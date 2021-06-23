/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

namespace alglin {

  template <typename t_Value>
  BorderedCR<t_Value>::BorderedCR( Utils::ThreadPool * TP )
  : m_baseValue("BorderedCR_values")
  , m_baseInteger("BorderedCR_integers")
  , m_basePointer("BorderedCR_pointers")
  , m_basePointerInteger("BorderedCR_pointer_integers")
  , m_superluValue("BorderedCR_superluValue")
  , m_superluInteger("BorderedCR_superluInteger")
  , m_number_of_blocks(0)
  , m_block_size(0)
  , m_qr(0)
  , m_qx(0)
  , m_nr(0)
  , m_nx(0)
  , m_Nr(0)
  , m_Nc(0)
  , n_x_2(0)
  , n_x_n(0)
  , n_x_nx(0)
  , nr_x_n(0)
  , nr_x_nx(0)
  , nr_x_qx(0)
  , m_last_selected(BORDERED_LAST_LU)
  , m_selected(BORDERED_LU)
  , m_last_must_use_PINV(false)
  , m_H0Nqp(nullptr)
  , m_Bmat(nullptr)
  , m_Cmat(nullptr)
  , m_Cqmat(nullptr)
  , m_Dmat(nullptr)
  , m_Emat(nullptr)
  , m_Tmat(nullptr)
  , m_Ttau(nullptr)
  , m_Work(nullptr)
  , m_Perm(nullptr)
  , m_Lwork(0)
  , m_Hmat(nullptr)
  , m_TP(TP)
  {
    m_available_thread = m_TP == nullptr ? 1 : m_TP->size();
    #ifdef LAPACK_WRAPPER_USE_OPENBLAS
    openblas_set_num_threads(1);
    goto_set_num_threads(1);
    #endif
  }

  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::allocate(
    integer nblock,
    integer n,
    integer qr,
    integer qx,
    integer nr,
    integer nx
  ) {

    if ( m_number_of_blocks == nblock &&
         m_block_size       == n      &&
         m_qr               == qr     &&
         m_qx               == qx     &&
         m_nr               == nr     &&
         m_nx               == nx ) return;

    m_number_of_blocks = nblock;
    m_block_size       = n;

    m_qr     = qr;
    m_qx     = qx;
    m_nr     = nr;
    m_nx     = nx;
    n_x_2    = n*2;
    n_x_n    = n*n;
    nr_x_n   = n*m_nr;
    n_x_nx   = n*m_nx;
    nr_x_nx  = m_nr*m_nx;
    nr_x_qx  = m_nr*m_qx;
    m_Nr     = n_x_2+m_nr+m_qr;
    m_Nc     = n_x_2+m_nx+m_qx;
    m_Tsize  = 2*n_x_n+n;

    integer N = std::max(m_Nr,m_Nc);
    m_Lwork = std::max(N,2*std::max(n_x_n,std::max(nr_x_n,n_x_nx)));

    real_type tmp; // get optimal allocation
    integer info = alglin::geqrf( m_Nr, m_Nc, nullptr, m_Nr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    if ( m_Lwork < integer(tmp) ) m_Lwork = integer(tmp);

    info = geqp3( m_Nr, m_Nc, nullptr, m_Nr, nullptr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( m_Lwork < integer(tmp) ) m_Lwork = integer(tmp);

    m_LworkT = 2*n_x_n;
    info = alglin::geqrf( n_x_2, n, nullptr, n_x_2, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    m_LworkQR = integer(tmp);

    info = alglin::geqp3( n_x_2, n, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( m_LworkQR < integer(tmp) ) m_LworkQR = integer(tmp);

    integer nnz = nblock*(n*(nx+2*n)+m_Tsize+n) +
                  (nblock+1)*nr*(n+qx) + (n+qr+m_Nr)*m_Nc +
                  m_Lwork + (m_LworkT+m_LworkQR+nr*(1+nx))*m_available_thread;
    integer innz = nblock*n + (3+n)*m_available_thread;

    m_baseValue.reallocate( size_t(nnz) );
    m_baseInteger.reallocate( size_t(innz) );
    m_basePointer.reallocate( size_t(4*m_available_thread) );
    m_basePointerInteger.reallocate( size_t(m_available_thread) );

    m_Bmat  = m_baseValue( size_t(nblock*n_x_nx) );
    m_Cmat  = m_baseValue( size_t((nblock+1)*nr_x_n) );
    m_Cqmat = m_baseValue( size_t(nr_x_qx) );
    m_Dmat  = m_baseValue( size_t(nblock*n_x_n) );
    m_Emat  = m_baseValue( size_t(nblock*n_x_n) );

    m_Tmat  = m_baseValue( size_t(nblock*m_Tsize) );
    m_Ttau  = m_baseValue( size_t(nblock*n) );
    m_Hmat  = m_baseValue( size_t(m_Nr*m_Nc) );
    m_H0Nqp = m_baseValue( size_t((n+qr)*m_Nc) );

    m_Work   = m_baseValue( size_t(m_Lwork) );
    m_Perm   = m_baseInteger( size_t(nblock*n) );
    m_iBlock = m_baseInteger( size_t(2*m_available_thread) );
    m_kBlock = m_baseInteger( size_t(m_available_thread) );

    m_perm_thread = m_basePointerInteger( size_t(m_available_thread) );
    m_xb_thread   = m_basePointer( size_t(m_available_thread) );
    m_WorkT       = m_basePointer( size_t(m_available_thread) );
    m_WorkQR      = m_basePointer( size_t(m_available_thread) );
    m_Fmat        = m_basePointer( size_t(m_available_thread) );

    // precompute partition for parallel computation
    for ( size_t nt = 0; nt < size_t(m_available_thread); ++nt ) {
      m_perm_thread[nt] = m_baseInteger( size_t(n) );
      m_xb_thread[nt]   = m_baseValue( size_t(nr) );
      m_WorkT[nt]       = m_baseValue( size_t(m_LworkT) );
      m_WorkQR[nt]      = m_baseValue( size_t(m_LworkQR) );
      m_Fmat[nt]        = m_baseValue( size_t(nr_x_nx) );
    }

    // calcolo partizionamento blocchi in modo che ogni thread
    // abbia almeno 2 blocchi di righe da processare.
    m_used_thread  = m_available_thread;
    m_reduced_nblk = 2*m_available_thread;
    while ( m_reduced_nblk > m_number_of_blocks )
      { m_reduced_nblk -= 2; --m_used_thread; }

    if ( m_used_thread <= 1 ) {
      m_used_thread = 1;
      m_iBlock[0]   = 0;
      m_iBlock[1]   = m_number_of_blocks;
    } else {
      m_iBlock[0] = 0;
      m_iBlock[1] = static_cast<integer>(m_number_of_blocks/m_used_thread);
      for ( integer nt = 1; nt < m_used_thread; ++nt ) {
        m_iBlock[2*nt+0] = m_iBlock[2*nt-1]+1;
        m_iBlock[2*nt+1] = static_cast<integer>(((nt+1)*m_number_of_blocks)/m_used_thread);
      }
      --m_reduced_nblk;
    }
  }

  /*\
   |       _
   |    __| |_   _ _ __
   |   / _` | | | | '_ \
   |  | (_| | |_| | |_) |
   |   \__,_|\__,_| .__/
   |              |_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::dup( BorderedCR const & M ) {
    m_available_thread = M.m_available_thread;
    m_TP               = M.m_TP;
    allocate( M.m_number_of_blocks, M.m_block_size, M.m_qr, M.m_qx, M.m_nr, M.m_nx );

    integer const & nblk = m_number_of_blocks;
    integer const & n    = m_block_size;

    alglin::copy( nblk*n_x_nx,     M.m_Bmat,    1, m_Bmat,    1 );
    alglin::copy( (nblk+1)*nr_x_n, M.m_Cmat,    1, m_Cmat,    1 );
    alglin::copy( nblk*n_x_n,      M.m_Dmat,    1, m_Dmat,    1 );
    alglin::copy( nblk*n_x_n,      M.m_Emat,    1, m_Emat,    1 );
    alglin::copy( nr_x_nx,         M.m_Fmat[0], 1, m_Fmat[0], 1 );
    alglin::copy( (n+m_qr)*m_Nc,   M.m_H0Nqp,   1, m_H0Nqp,   1 );
    alglin::copy( nr_x_qx,         M.m_Cqmat,   1, m_Cqmat,   1 );
  }

  /*\
   |   _        __
   |  (_)_ __  / _| ___
   |  | | '_ \| |_ / _ \
   |  | | | | |  _| (_) |
   |  |_|_| |_|_|  \___/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::info( ostream_type & stream ) const {
    fmt::print( stream,
      "rows   = {}\n"
      "cols   = {}\n"
      "nblock = {}\n"
      "n      = {}\n"
      "qr     = {}\n"
      "qx     = {}\n"
      "nr     = {}\n"
      "nx     = {}\n",
      numRows(), numCols(), m_number_of_blocks,
      m_block_size, m_qr, m_qx, m_nr, m_nx
    );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroD() {
    integer const & nblk = m_number_of_blocks;
    alglin::zero( nblk*n_x_n, m_Dmat, 1 );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroE() {
    integer const & nblk = m_number_of_blocks;
    alglin::zero( nblk*n_x_n, m_Emat, 1 );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroB() {
    integer const & nblk = m_number_of_blocks;
    alglin::zero( nblk*n_x_nx, m_Bmat, 1 );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroF()
  { alglin::zero( nr_x_nx, m_Fmat[0], 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroH() {
    integer const & n = m_block_size;
    alglin::zero( (n+m_qr)*m_Nc, m_H0Nqp, 1 );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroC() {
    integer const & nblk = m_number_of_blocks;
    alglin::zero( (nblk+1)*nr_x_n, m_Cmat, 1 );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroCq()
  { alglin::zero( nr_x_qx, m_Cqmat, 1 ); }

  /*\
   |   _                 _ ____        _   _
   |  | | ___   __ _  __| | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` |  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_|____/ \___/ \__|\__\___/|_| |_| |_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom(
    real_type const H0[], integer ld0,
    real_type const HN[], integer ldN,
    real_type const Hq[], integer ldQ,
    real_type const Hp[], integer ldP
  ) {
    integer const & n = m_block_size;
    // (n+qr) x ( n + n + qx + nx )
    integer     m = n + m_qr;
    real_type * H = m_H0Nqp;
    alglin::gecopy( m, n,    H0, ld0, H, m ); H += m * n;
    alglin::gecopy( m, n,    HN, ldN, H, m ); H += m * n;
    alglin::gecopy( m, m_qx, Hq, ldQ, H, m ); H += m * m_qx;
    alglin::gecopy( m, m_nx, Hp, ldP, H, m );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom(
    MatrixWrapper<real_type> const & H0,
    MatrixWrapper<real_type> const & HN,
    MatrixWrapper<real_type> const & Hq,
    MatrixWrapper<real_type> const & Hp
  ) {
    integer const & n = m_block_size;
    integer m = n + m_qr;

    UTILS_ASSERT(
      H0.numRows() == m && H0.numCols() == n,
      "loadBottom, bad dimension size(H0) = {} x {} expected {} x {}\n",
      H0.numRows(), H0.numCols(), m, n
    );

    UTILS_ASSERT(
      HN.numRows() == m && HN.numCols() == n,
      "loadBottom, bad dimension size(HN) = {} x {} expected {} x {}\n",
      HN.numRows(), HN.numCols(), m, n
    );

    UTILS_ASSERT(
      Hq.numRows() == m && Hq.numCols() == m_qx,
      "loadBottom, bad dimension size(Hq) = {} x {} expected {} x {}\n",
      Hq.numRows(), Hq.numCols(), m, m_qx
    );

    UTILS_ASSERT(
      Hp.numRows() == m && Hp.numCols() == m_nx,
      "loadBottom, bad dimension size(Hp) = {} x {} expected {} x {}\n",
      Hp.numRows(), Hp.numCols(), m, m_nx
    );

    loadBottom(
      H0.data(), H0.lDim(),
      HN.data(), HN.lDim(),
      Hq.data(), Hq.lDim(),
      Hp.data(), Hp.lDim()
    );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom( real_type const _H0Nqp[], integer ldH ) {
    integer const & n = m_block_size;
    integer nq = n + m_qr;
    alglin::gecopy( nq, m_Nc, _H0Nqp, ldH, m_H0Nqp, nq );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom( MatrixWrapper<real_type> const & H ) {
    integer const & n = m_block_size;
    integer m = n + m_qr;
    UTILS_ASSERT(
      H.numRows() == m && H.numCols() == m_Nc,
      "loadBottom, bad dimension size(H) = {} x {} expected {} x {}\n",
      H.numRows(), H.numCols(), m, m_Nc
    );
    loadBottom( H.data(), H.lDim() );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    real_type const C0[], integer ld0,
    real_type const CN[], integer ldN,
    real_type const Cq[], integer ldCq,
    real_type const F[],  integer ldF
  ) {
    integer const & nblk = m_number_of_blocks;
    integer const & n    = m_block_size;
    // (n+qr) x ( n + n + qx + nx )
    gecopy( m_nr, n,    C0, ld0,  m_Cmat,            m_nr );
    gecopy( m_nr, n,    CN, ldN,  m_Cmat+nblk*n_x_n, m_nr );
    gecopy( m_nr, m_qx, Cq, ldCq, m_Cqmat,           m_nr );
    gecopy( m_nr, m_nx, F,  ldF,  m_Fmat[0],         m_nr );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    MatrixWrapper<real_type> const & C0,
    MatrixWrapper<real_type> const & CN,
    MatrixWrapper<real_type> const & Cq,
    MatrixWrapper<real_type> const & F
  ) {
    integer const & n = m_block_size;

    UTILS_ASSERT(
      C0.numRows() == m_nr && C0.numCols() == n,
      "loadBottom2, bad dimension size(C0) = {} x {} expected {} x {}\n",
      C0.numRows(), C0.numCols(), m_nr, n
    );

    UTILS_ASSERT(
      CN.numRows() == m_nr && CN.numCols() == n,
      "loadBottom2, bad dimension size(CN) = {} x {} expected {} x {}\n",
      CN.numRows(), CN.numCols(), m_nr, n
    );

    UTILS_ASSERT(
      Cq.numRows() == m_nr && Cq.numCols() == m_qx,
      "loadBottom2, bad dimension size(Cq) = {} x {} expected {} x {}\n",
      Cq.numRows(), Cq.numCols(), m_nr, m_qx
    );

    UTILS_ASSERT(
      F.numRows() == m_nr && F.numCols() == m_nx,
      "loadBottom2, bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), m_nr, m_nx
    );

    loadBottom2(
      C0.data(), C0.lDim(),
      CN.data(), CN.lDim(),
      Cq.data(), Cq.lDim(),
      F.data(),  F.lDim()
    );

  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2( MatrixWrapper<real_type> const & H ) {
    integer const & nblock = m_number_of_blocks;
    real_type const * ptr = H.data();
    integer           ld  = H.lDim();
    this->loadC(  0,      ptr, ld ); ptr += nr_x_n;
    this->loadC(  nblock, ptr, ld ); ptr += nr_x_n;
    this->loadCq( ptr, ld ); ptr += m_nr * m_qx;
    this->loadF(  ptr, ld );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB( integer nbl, real_type const B[], integer ldB ) {
    integer const & n = m_block_size;
    alglin::gecopy( n, m_nx, B, ldB, m_Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB(
    integer                          nbl,
    MatrixWrapper<real_type> const & B
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      B.numRows() == n && B.numCols() == m_nx,
      "loadB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), n, m_nx
    );
    alglin::gecopy( n, m_nx, B.data(), B.lDim(), m_Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB( integer nbl, real_type const B[], integer ldB ) {
    integer const & n = m_block_size;
    real_type * BB = m_Bmat + nbl*n_x_nx;
    alglin::geadd( n, m_nx, 1.0, B, ldB, 1.0, BB, n, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB(
    integer                          nbl,
    MatrixWrapper<real_type> const & B
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      B.numRows() == n && B.numCols() == m_nx,
      "addtoB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), n, m_nx
    );
    real_type * BB = m_Bmat + nbl*n_x_nx;
    alglin::geadd( n, m_nx, 1.0, B.data(), B.lDim(), 1.0, BB, n, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC( integer nbl, real_type const C[], integer ldC ) {
    integer const & n = m_block_size;
    alglin::gecopy( m_nr, n, C, ldC, m_Cmat + nbl*nr_x_n, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC(
    integer                          nbl,
    MatrixWrapper<real_type> const & C
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == n,
      "loadC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, n
    );
    real_type * CC = m_Cmat + nbl*nr_x_n;
    alglin::gecopy( m_nr, n, C.data(), C.lDim(), CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC( integer nbl, real_type const C[], integer ldC ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      ldC >= m_nr, "addtoC( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    real_type * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n, 1.0, C, ldC, 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC(
    integer                          nbl,
    MatrixWrapper<real_type> const & C
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == n,
      "addtoC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, n
    );
    real_type * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n, 1.0, C.data(), C.lDim(), 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2( integer nbl, real_type const C[], integer ldC ) {
    UTILS_ASSERT(
      ldC >= m_nr, "addtoC2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    real_type * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n_x_2, 1.0, C, ldC, 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2(
    integer                          nbl,
    MatrixWrapper<real_type> const & C
  ) {
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == n_x_2,
      "addtoC2( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, n_x_2
    );
    real_type * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n_x_2, 1.0, C.data(), C.lDim(), 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD( integer nbl, real_type const D[], integer ldD ) {
    integer const & n = m_block_size;
    alglin::gecopy( n, n, D, ldD, m_Dmat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD(
    integer                          nbl,
    MatrixWrapper<real_type> const & D
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      D.numRows() == n && D.numCols() == n,
      "loadD( {}, D) bad dimension size(D) = {} x {} expected {} x {}\n",
      nbl, D.numRows(), D.numCols(), n, n
    );
    real_type * DD = m_Dmat + nbl*n_x_n;
    alglin::gecopy( n, n, D.data(), D.lDim(), DD, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE( integer nbl, real_type const E[], integer ldE ) {
    integer const & n = m_block_size;
    real_type * EE = m_Emat + nbl*n_x_n;
    alglin::gecopy( n, n, E, ldE, EE, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE(
    integer                          nbl,
    MatrixWrapper<real_type> const & E
  ) {
    integer const & n = m_block_size;
    UTILS_ASSERT(
      E.numRows() == n && E.numCols() == n,
      "loadE( {}, E) bad dimension size(E) = {} x {} expected {} x {}\n",
      nbl, E.numRows(), E.numCols(), n, n
    );
    real_type * EE = m_Emat + nbl*n_x_n;
    alglin::gecopy( n, n, E.data(), E.lDim(), EE, n );
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDE( integer nbl, real_type const DE[], integer ldDE ) {
    integer const & n = m_block_size;
    real_type * DD = m_Dmat + nbl*n_x_n;
    real_type * EE = m_Emat + nbl*n_x_n;
    alglin::gecopy( n, n, DE, ldDE, DD, n ); DE += n*ldDE;
    alglin::gecopy( n, n, DE, ldDE, EE, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDEB( integer nbl, real_type const DEB[], integer ldDEB ) {
    integer const & n = m_block_size;
    real_type * DD = m_Dmat + nbl*n_x_n;
    real_type * EE = m_Emat + nbl*n_x_n;
    real_type * BB = m_Bmat + nbl*n_x_nx;
    alglin::gecopy( n, n,    DEB, ldDEB, DD, n ); DEB += n*ldDEB;
    alglin::gecopy( n, n,    DEB, ldDEB, EE, n ); DEB += n*ldDEB;
    alglin::gecopy( n, m_nx, DEB, ldDEB, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( real_type const F[], integer ldF )
  { alglin::gecopy( m_nr, m_nx, F, ldF, m_Fmat[0], m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( MatrixWrapper<real_type> const & F ) {
    UTILS_ASSERT(
      F.numRows() == m_nr && F.numCols() == m_nx,
      "loadF(F) bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), m_nr, m_nx
    );
    alglin::gecopy( m_nr, m_nx, F.data(), F.lDim(), m_Fmat[0], m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( real_type const F[], integer ldF )
  { alglin::gecopy( m_nr, m_nx, F, ldF, m_Fmat[0], m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( MatrixWrapper<real_type> const & F ) {
    UTILS_ASSERT(
      F.numRows() == m_nr && F.numCols() == m_nx,
      "addtoF(F) bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), m_nr, m_nx
    );
    alglin::geadd(
      m_nr, m_nx,
      1.0, F.data(), F.lDim(),
      1.0, m_Fmat[0], m_nr,
      m_Fmat[0], m_nr
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( real_type const Cq[], integer ldC ) {
    alglin::gecopy( m_nr, m_qx, Cq, ldC, m_Cqmat, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( MatrixWrapper<real_type> const & Cq ) {
    UTILS_ASSERT(
      Cq.numRows() == m_nr && Cq.numCols() == m_qx,
      "loadCq(Cq) bad dimension size(Cq) = {} x {} expected {} x {}\n",
      Cq.numRows(), Cq.numCols(), m_nr, m_qx
    );
    alglin::gecopy( m_nr, m_qx, Cq.data(), Cq.lDim(), m_Cqmat, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCqF( real_type const CqF[], integer ldCF ) {
    integer const & nr = m_nr;
    integer const & nx = m_nx;
    integer const & qx = m_qx;
    alglin::gecopy( nr, qx, CqF, ldCF, m_Cqmat,   nr ); CqF += qx*ldCF;
    alglin::gecopy( nr, nx, CqF, ldCF, m_Fmat[0], nr );
  }

  /*\
   |   _____          _             _
   |  |  ___|_ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_CR() {
    integer const & nblock = m_number_of_blocks;

    if ( m_used_thread > 1 ) {
      // launch thread
      for ( integer nt = 1; nt < m_used_thread; ++nt ) {
        // fill zero F(...) only that of the extra threads
        if ( nr_x_nx > 0 ) alglin::zero( nr_x_nx, m_Fmat[size_t(nt)], 1 );
        m_TP->run( nt, &BorderedCR<t_Value>::factorize_block, this, nt );
      }
      factorize_block(0);
      // wait thread
      m_TP->wait_all();
      // accumulate F(...)
      if ( nr_x_nx > 0 )
        for ( integer nt = 1; nt < m_used_thread; ++nt )
          alglin::axpy( nr_x_nx, 1.0, m_Fmat[size_t(nt)], 1, m_Fmat[0], 1 );
      factorize_reduced();
    } else {
      m_iBlock[0] = 0;
      m_iBlock[1] = nblock;
      factorize_block(0);
    }
    UTILS_ASSERT0(
      load_and_factorize_last(),
      "BorderedCR<t_Value>::factorize_CR, load_and_factorize_last failed\n"
    );
  }

  /*\
    Compute   / L \ (U) = / LU \ = / TOP    \
              \ M /       \ MU /   \ BOTTOM /
  \*/
  template <typename t_Value>
  void
  BorderedCR<t_Value>::buildT(
    integer         nth,
    real_type const TOP[],
    real_type const BOTTOM[],
    real_type       T[],
    integer         iperm[]
  ) const {
    integer const & n = m_block_size;
    alglin::gecopy( n, n, TOP,    n, T,   n_x_2 );
    alglin::gecopy( n, n, BOTTOM, n, T+n, n_x_2 );
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      info = alglin::getrf( n_x_2, n, T, n_x_2, iperm );
      UTILS_ASSERT(
        info == 0, "BorderedCR::buildT getrf::INFO = {}\n", info
      );
      break;
    case BORDERED_QR:
      info = alglin::geqrf(
        n_x_2, n, T, n_x_2, T+2*n_x_n, m_WorkQR[size_t(nth)], m_LworkQR
      );
      UTILS_ASSERT(
        info == 0, "BorderedCR::buildT geqrf::INFO = {}\n", info
      );
      break;
    case BORDERED_QRP:
      {
        integer * P = m_perm_thread[size_t(nth)];
        std::fill_n( P, n, 0 );
        info = alglin::geqp3(
          n_x_2, n, T, n_x_2, P, T+2*n_x_n, m_WorkQR[size_t(nth)], m_LworkQR
        );
        if ( info == 0 ) permutation_to_exchange( n, P, iperm );
        UTILS_ASSERT(
          info == 0, "BorderedCR::buildT geqp3::INFO = {}\n", info
        );
      }
      break;
    case BORDERED_SUPERLU:
      UTILS_ERROR0( "BorderedCR::buildT, cannot be used with SUPERLU\n" );
      break;
    }
  }

  /*\
   APPLY
      / -M  I \ / L^(-1)   \ / TOP    \
      \  I  0 / \        I / \ BOTTOM /
      = / -M  I \ / L^(-1) TOP \
        \  I  0 / \  BOTTOM    /
      = / BOTTOM - M L^(-1) TOP \
        \ L^(-1) TOP            /

   ON
      / -M  I \ / L^(-1)   \ / L \ (U)
      \  I  0 / \        I / \ M /
        = / -M  I \ / I \ (U)
          \  I  0 / \ M /
        = / 0 \ (U)   / 0 \
          \ I /     = \ U /
  \*/
  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT(
    integer         nth,
    real_type const T[],
    integer   const iperm[],
    real_type       TOP[],
    integer         ldTOP,
    real_type       BOTTOM[],
    integer         ldBOTTOM,
    integer         ncol
  ) const {
    integer const & n = m_block_size;
    real_type     * W = m_WorkT[size_t(nth)];
    alglin::gecopy( n, ncol, TOP,    ldTOP,    W,   n_x_2 );
    alglin::gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, n_x_2 );
    // Apply row interchanges to the right hand sides.
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      info = alglin::swaps( ncol, W, n_x_2, 0, n-1, iperm, 1 );
      UTILS_ASSERT( info == 0, "BorderedCR::applyT INFO = {}\n", info );
      alglin::trsm(
        LEFT, LOWER, NO_TRANSPOSE, UNIT,
        n, ncol, 1.0, T, n_x_2, W, n_x_2
      );
      alglin::gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        n, ncol, n,
        -1.0, T+n, n_x_2,  // M
              W,   n_x_2,  // L^(-1) TOP
         1.0, W+n, n_x_2   // TOP = BOTTOM - M L^(-1) TOP
      );
      break;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = alglin::ormqr(
        LEFT, TRANSPOSE,
        n_x_2, ncol, // righe x colonne
        n,           // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        m_WorkQR[size_t(nth)], m_LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      UTILS_ERROR0( "BorderedCR::applyT, cannot be used with SUPERLU\n" );
      break;
    }
    alglin::gecopy( n, ncol, W+n, n_x_2, TOP,    ldTOP    );
    alglin::gecopy( n, ncol, W,   n_x_2, BOTTOM, ldBOTTOM );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT(
    integer         nth,
    real_type const T[],
    integer   const iperm[],
    real_type       TOP[],
    real_type       BOTTOM[]
  ) const {
    integer const & n = m_block_size;
    real_type * W = m_WorkT[size_t(nth)];
    size_t nn = size_t(n)*sizeof(real_type);
    memcpy( W,   TOP,    nn );
    memcpy( W+n, BOTTOM, nn );
    //copy( n, TOP,    1, W,   1 );
    //copy( n, BOTTOM, 1, W+n, 1 );
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      // Apply row interchanges to the right hand sides.
      info = swaps( 1, W, n_x_2, 0, n-1, iperm, 1 );
      if ( info == 0 ) {
        alglin::trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, n_x_2, W, 1 );
        alglin::gemv(
          NO_TRANSPOSE,
          n, n,
          -1.0, T+n, n_x_2,
                W,       1,
           1.0, W+n, 1
        );
      }
      break;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = alglin::ormqr(
        LEFT, TRANSPOSE,
        n_x_2, 1, // righe x colonne
        n,        // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        m_WorkQR[size_t(nth)], m_LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      UTILS_ERROR0( "BorderedCR::applyT, cannot be used with SUPERLU\n" );
      break;
    }
    UTILS_ASSERT( info == 0, "BorderedCR::applyT INFO = {}\n", info );
    memcpy( TOP,    W+n, nn );
    memcpy( BOTTOM, W,   nn );
    //copy( n, W+n, 1, TOP,    1 );
    //copy( n, W,   1, BOTTOM, 1 );
  }

  /*\
   |    __            _             _             _     _            _
   |   / _| __ _  ___| |_ ___  _ __(_)_______    | |__ | | ___   ___| | __
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \   | '_ \| |/ _ \ / __| |/ /
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/   | |_) | | (_) | (__|   <
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|___|_.__/|_|\___/ \___|_|\_\
   |                                        |_____|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_block( integer nth ) {
    integer const & n = m_block_size;

    integer iblock = m_iBlock[2*nth+0];
    integer eblock = m_iBlock[2*nth+1];
    integer nblk   = eblock - iblock;

    real_type * Bmat0 = m_Bmat + iblock*n_x_nx;
    real_type * Cmat0 = m_Cmat + iblock*nr_x_n;
    real_type * Dmat0 = m_Dmat + iblock*n_x_n;
    real_type * Emat0 = m_Emat + iblock*n_x_n;
    real_type * T0    = m_Tmat + iblock*m_Tsize;
    integer   * P0    = m_Perm + iblock*n;

    real_type * Fmat_th = m_Fmat[size_t(nth)];

    integer k = 1;
    while ( k < nblk ) {

      real_type * Bjp = Bmat0;
      real_type * Bj  = Bmat0 + k*n_x_nx;
      real_type * Cjp = Cmat0;
      real_type * Cj  = Cmat0 + k*nr_x_n;
      real_type * Djp = Dmat0;
      real_type * Dj  = Dmat0 + k*n_x_n;
      real_type * Ejp = Emat0;
      real_type * Ej  = Emat0 + k*n_x_n;
      real_type * T   = T0    + k*m_Tsize;
      integer   * P   = P0    + k*n;

      // -----------------------------------------------------------------------
      integer k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {

        buildT( nth, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( nth, T, P, Djp, n, Dj, n, n );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( nth, T, P, Ejp, n, Ej, n, n );

        if ( m_nx > 0 ) applyT( nth, T, P, Bjp, n, Bj, n, m_nx );

        if ( m_nr > 0 ) {

          if ( m_selected == BORDERED_QRP ) {
            integer i = n;
            do {
              --i;
              if ( P[i] > i ) swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            m_nr, n, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Dj,  n,
             1.0, Cjp, m_nr
          );

          integer     jpp = std::min(j+k,eblock);
          real_type * Cpp = m_Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Ej,  n,
             1.0, Cpp, m_nr
          );
        }

        if ( nr_x_nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_nx, n,
            -1.0, Cj,      m_nr,
                  Bj,      n,
             1.0, Fmat_th, m_nr // solo accumulo!
          );

        // NEXT STEP
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Djp += k_x_2*n_x_n;
        Dj  += k_x_2*n_x_n;
        Ejp += k_x_2*n_x_n;
        Ej  += k_x_2*n_x_n;
        Bj  += k_x_2*n_x_nx;
        Bjp += k_x_2*n_x_nx;
        Cj  += k_x_2*nr_x_n;
        Cjp += k_x_2*nr_x_n;
      }
      k *= 2;
    }
    m_kBlock[nth] = k;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_reduced() {
    integer const & n = m_block_size;
    integer k = 1;
    while ( k < m_reduced_nblk ) {
      // -----------------------------------------------------------------------
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer j  = m_iBlock[jj];
        integer jp = m_iBlock[jj-k];
        real_type * T   = m_Tmat + j*m_Tsize;
        integer   * P   = m_Perm + j*n;
        real_type * Djp = m_Dmat + jp*n_x_n;
        real_type * Dj  = m_Dmat + j*n_x_n;
        real_type * Ejp = m_Emat + jp*n_x_n;
        real_type * Ej  = m_Emat + j*n_x_n;

        buildT( 0, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( 0, T, P, Djp, n, Dj, n, n );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( 0, T, P, Ejp, n, Ej, n, n );

        if ( m_nx > 0 ) {
          real_type * Bj  = m_Bmat + j*n_x_nx;
          real_type * Bjp = m_Bmat + jp*n_x_nx;
          applyT( 0, T, P, Bjp, n, Bj, n, m_nx );
        }

        if ( m_nr > 0 ) {
          real_type * Cj  = m_Cmat + j*nr_x_n;
          real_type * Cjp = m_Cmat + jp*nr_x_n;
          if ( m_selected == BORDERED_QRP ) {
            integer i = n;
            do {
              --i;
              if ( P[i] > i ) swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            m_nr, n, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Dj,  n,
             1.0, Cjp, m_nr
          );

          integer     jpp = m_iBlock[std::min(jj+k,m_reduced_nblk)];
          real_type * Cpp = m_Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, n, n,
            -1.0, Cj,  m_nr,
                  Ej,  n,
             1.0, Cpp, m_nr
          );

        }

        if ( nr_x_nx > 0 ) {
          real_type * Cj = m_Cmat + j*nr_x_n;
          real_type * Bj = m_Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_nx, n,
            -1.0, Cj,        m_nr,
                  Bj,        n,
             1.0, m_Fmat[0], m_nr
          );
        }
      }
      k *= 2;
    }
  }

  /*\
   |   _                 _                   _
   |  | | ___   __ _  __| |   __ _ _ __   __| |
   |  | |/ _ \ / _` |/ _` |  / _` | '_ \ / _` |
   |  | | (_) | (_| | (_| | | (_| | | | | (_| |
   |  |_|\___/ \__,_|\__,_|  \__,_|_| |_|\__,_|
   |    __            _             _           _           _
   |   / _| __ _  ___| |_ ___  _ __(_)_______  | | __ _ ___| |_
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \ | |/ _` / __| __|
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/ | | (_| \__ \ |_
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___| |_|\__,_|___/\__|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::load_and_factorize_last() {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
   /*
    //    n   n  qx  nx
    //  +---+---+---+---+
    //  | D | E | 0 | B | n
    //  +===+===+===+===+
    //  |   |   |   |   |
    //  |H0 |HN |Hq |Hp | n+qr
    //  |   |   |   |   |
    //  +---+---+---+---+
    //  | C | C |Cq | F | nr
    //  +---+---+---+---+
    */
    real_type * Cnb = m_Cmat + nblock*nr_x_n;
    real_type * W0  = m_Hmat;
    real_type * WN  = W0+n*m_Nr;
    real_type * Wq  = WN+n*m_Nr;
    real_type * Wp  = Wq+m_qx*m_Nr;

    alglin::gecopy( n, n,    m_Dmat, n, W0, m_Nr );
    alglin::gecopy( n, n,    m_Emat, n, WN, m_Nr );
    alglin::gezero( n, m_qx,            Wq, m_Nr );
    alglin::gecopy( n, m_nx, m_Bmat, n, Wp, m_Nr );

    alglin::gecopy( n+m_qr, m_Nc, m_H0Nqp, n+m_qr, m_Hmat+n, m_Nr );

    integer offs = n_x_2+m_qr;

    alglin::gecopy( m_nr, n,    m_Cmat,    m_nr, W0+offs, m_Nr );
    alglin::gecopy( m_nr, n,    Cnb,       m_nr, WN+offs, m_Nr );
    alglin::gecopy( m_nr, m_qx, m_Cqmat,   m_nr, Wq+offs, m_Nr );
    alglin::gecopy( m_nr, m_nx, m_Fmat[0], m_nr, Wp+offs, m_Nr );

    bool ok = false;
    switch ( m_last_selected ) {
    case BORDERED_LAST_LU:
      ok = m_last_lu.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_LUPQ:
      ok = m_last_lupq.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_QR:
      ok = m_last_qr.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_QRP:
      ok = m_last_qrp.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_SVD:
      ok = m_last_svd.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_LSS:
      ok = m_last_lss.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_LSY:
      ok = m_last_lsy.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      break;
    case BORDERED_LAST_PINV:
      break;
    }
    if ( ok ) {
      m_last_must_use_PINV = false;
    } else {
      m_last_must_use_PINV = true;
      // try to facttorize using PINV
      if ( m_Nc > m_Nr )
        ok = m_last_pinv.t_factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
      else
        ok = m_last_pinv.factorize( m_Nr, m_Nc, m_Hmat, m_Nr );
    }
    return ok;
  }

  /*\
   |             _               _           _
   |   ___  ___ | |_   _____    | | __ _ ___| |_
   |  / __|/ _ \| \ \ / / _ \   | |/ _` / __| __|
   |  \__ \ (_) | |\ V /  __/   | | (_| \__ \ |_
   |  |___/\___/|_| \_/ \___|___|_|\__,_|___/\__|
   |                       |_____|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_last( real_type x[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    real_type * X = x + (nblock-1)*n;
    // sposto primo blocco rhs in fondo
    swap( n, X, 1, x, 1 ); // uso x stesso come temporaneo
    bool ok = false;
    if ( m_last_must_use_PINV ) {
      if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( X, 1, X, 1 );
      else               ok = m_last_pinv.mult_inv( X, 1, X, 1 );
    } else {
      switch ( m_last_selected ) {
        case BORDERED_LAST_LU:   ok = m_last_lu.solve( X );   break;
        case BORDERED_LAST_LUPQ: ok = m_last_lupq.solve( X ); break;
        case BORDERED_LAST_QR:   ok = m_last_qr.solve( X );   break;
        case BORDERED_LAST_QRP:  ok = m_last_qrp.solve( X );  break;
        case BORDERED_LAST_SVD:  ok = m_last_svd.solve( X );  break;
        case BORDERED_LAST_LSS:  ok = m_last_lss.solve( X );  break;
        case BORDERED_LAST_LSY:  ok = m_last_lsy.solve( X );  break;
        case BORDERED_LAST_PINV: break;
      }
    }
    if ( ok ) swap( n, X, 1, x, 1 );
    return ok;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_last(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    real_type * X = x + (nblock-1)*n;
    // sposto primo blocco rhs in fondo
    for ( integer i = 0; i < nrhs; ++i )
      swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    bool ok = false;
    if ( m_last_must_use_PINV ) {
      if ( m_Nc > m_Nr ) ok = m_last_pinv.t_mult_inv( nrhs, X, ldX, X, ldX );
      else               ok = m_last_pinv.mult_inv( nrhs, X, ldX, X, ldX );
    } else {
      switch ( m_last_selected ) {
        case BORDERED_LAST_LU:   ok = m_last_lu.solve( nrhs, X, ldX );   break;
        case BORDERED_LAST_LUPQ: ok = m_last_lupq.solve( nrhs, X, ldX ); break;
        case BORDERED_LAST_QR:   ok = m_last_qr.solve( nrhs, X, ldX );   break;
        case BORDERED_LAST_QRP:  ok = m_last_qrp.solve( nrhs, X, ldX );  break;
        case BORDERED_LAST_SVD:  ok = m_last_svd.solve( nrhs, X, ldX );  break;
        case BORDERED_LAST_LSS:  ok = m_last_lss.solve( nrhs, X, ldX );  break;
        case BORDERED_LAST_LSY:  ok = m_last_lsy.solve( nrhs, X, ldX );  break;
        case BORDERED_LAST_PINV: break;
      }
    }
    if ( ok )
      for ( integer i = 0; i < nrhs; ++i )
        swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    return ok;
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_CR( real_type x[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    real_type * xb = x + (nblock+1)*n + m_qr; // deve essere b!
    if ( m_used_thread > 1 ) {
      if ( m_nr > 0 ) {
        for ( integer nt = 1; nt < m_used_thread; ++nt ) {
          alglin::zero( m_nr, m_xb_thread[size_t(nt)], 1 );
          m_TP->run(
            nt, &BorderedCR<t_Value>::forward, this,
            nt, x, m_xb_thread[size_t(nt)]
          );
        }
        alglin::zero( m_nr, m_xb_thread[0], 1 );
        forward(0,x,xb);
        m_TP->wait_all();
        for ( integer nt = 1; nt < m_used_thread; ++nt )
          axpy( m_nr, 1.0, m_xb_thread[size_t(nt)], 1, xb, 1 );
      } else {
        for ( integer nt = 1; nt < m_used_thread; ++nt )
          m_TP->run(
            nt, &BorderedCR<t_Value>::forward, this,
            nt, x, m_xb_thread[size_t(nt)]
          );
        forward(0,x,xb);
        m_TP->wait_all();
      }
      forward_reduced(x,xb);
    } else {
      forward(0,x,xb);
    }

    if ( !solve_last( x ) ) return false;

    if ( m_used_thread > 1 ) {
      backward_reduced(x);
      for ( integer nt = 1; nt < m_used_thread; ++nt )
        m_TP->run( nt, &BorderedCR<t_Value>::backward, this, nt, x );
      backward(0,x);
      m_TP->wait_all();
    } else {
      backward(0,x);
    }
    return true;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_CR(
    integer   nrhs,
    real_type rhs[],
    integer   ldRhs
  ) const {
    if ( m_used_thread > 1 ) {
      for ( integer nt = 1; nt < m_used_thread; ++nt )
        m_TP->run(
          nt, &BorderedCR<t_Value>::forward_n, this, nt, nrhs, rhs, ldRhs
        );
      forward_n( 0, nrhs, rhs, ldRhs );
      m_TP->wait_all();
      forward_n_reduced( nrhs, rhs, ldRhs );
    } else {
      forward_n( 0, nrhs, rhs, ldRhs );
    }

    if ( !solve_last( nrhs, rhs, ldRhs ) ) return false;

    if ( m_used_thread > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < m_used_thread; ++nt )
        m_TP->run(
          nt, &BorderedCR<t_Value>::backward_n, this, nt, nrhs, rhs, ldRhs
        );
      backward_n( 0, nrhs, rhs, ldRhs );
      m_TP->wait_all();
    } else {
      backward_n( 0, nrhs, rhs, ldRhs );
    }
    return true;
  }

  /*\
   |    __                                  _
   |   / _| ___  _ ____      ____ _ _ __ __| |
   |  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
   |  |  _| (_) | |   \ V  V / (_| | | | (_| |
   |  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward(
    integer   nth,
    real_type x[],
    real_type xb[]
  ) const {
    integer const & n = m_block_size;

    integer iblock = m_iBlock[2*nth+0];
    integer eblock = m_iBlock[2*nth+1];
    integer nblk   = eblock - iblock;
    real_type * x0 = x      + iblock*n;
    real_type * T0 = m_Tmat + iblock*m_Tsize;
    integer   * P0 = m_Perm + iblock*n;
    real_type * C0 = m_Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      real_type * xj  = x0 + k*n;
      real_type * xjp = x0;
      real_type * T   = T0 + k*m_Tsize;
      integer   * P   = P0 + k*n;
      real_type * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;
      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, xj );
        if ( m_nr > 0 )
          gemv(
            NO_TRANSPOSE,
            m_nr, n, -1.0,
            Cj, m_nr, xj, 1,
            1.0, xb, 1
          ); // solo accumulato
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Cj  += k_x_2*nr_x_n;
      }
      k *= 2;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_n(
    integer   nth,
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    real_type * xb = x + (nblock+1)*n + m_qr;
    integer iblock = m_iBlock[2*nth+0];
    integer eblock = m_iBlock[2*nth+1];
    integer nblk   = eblock - iblock;
    real_type * x0 = x      + iblock*n;
    real_type * T0 = m_Tmat + iblock*m_Tsize;
    integer   * P0 = m_Perm + iblock*n;
    real_type * C0 = m_Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      real_type * xj  = x0 + k*n;
      real_type * xjp = x0;
      real_type * T   = T0 + k*m_Tsize;
      integer   * P   = P0 + k*n;
      real_type * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;

      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( m_nr > 0 ) {
          m_spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, nrhs, n,
            -1.0, Cj, m_nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          m_spin.unlock();
        }
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*m_Tsize;
        P   += k_x_2*n;
        Cj  += k_x_2*nr_x_n;
      }
      k *= 2;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_reduced(
    real_type x[],
    real_type xb[]
  ) const {
    integer const & n = m_block_size;

    integer k = 1;
    while ( k < m_reduced_nblk ) {
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer j  = m_iBlock[jj];
        integer jp = m_iBlock[jj-k];
        real_type const * T   = m_Tmat + j*m_Tsize;
        integer   const * P   = m_Perm + j*n;
        real_type       * xj  = x + j*n;
        real_type       * xjp = x + jp*n;
        applyT( 0, T, P, xjp, xj );
        if ( m_nr > 0 ) {
          real_type * Cj = m_Cmat + j*nr_x_n;
          alglin::gemv(
            NO_TRANSPOSE, m_nr, n,
            -1.0, Cj, m_nr, xj, 1, 1.0, xb, 1
          );
        }
      }
      k *= 2;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_n_reduced(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    real_type * xb = x + (nblock+1)*n + m_qr;
    integer k = 1;
    while ( k < m_reduced_nblk ) {
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer j  = m_iBlock[jj];
        integer jp = m_iBlock[jj-k];
        real_type const * T   = m_Tmat + j*m_Tsize;
        integer   const * P   = m_Perm + j*n;
        real_type       * xj  = x + j*n;
        real_type       * xjp = x + jp*n;
        applyT( 0, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( m_nr > 0 ) {
          real_type * Cj = m_Cmat + j*nr_x_n;
          m_spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, nrhs, n,
            -1.0, Cj, m_nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          m_spin.unlock();
        }
      }
      k *= 2;
    }
  }

  /*\
   |   _                _                           _
   |  | |__   __ _  ___| | ____      ____ _ _ __ __| |
   |  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
   |  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
   |  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward( integer nth, real_type x[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    real_type * xn = x + (nblock+1)*n + m_qx;
    integer iblock = m_iBlock[2*nth+0];
    integer eblock = m_iBlock[2*nth+1];
    real_type * x0 = x      + iblock*n;
    real_type * B0 = m_Bmat + iblock*n_x_nx;
    real_type * D0 = m_Dmat + iblock*n_x_n;
    real_type * E0 = m_Emat + iblock*n_x_n;
    real_type * T0 = m_Tmat + iblock*m_Tsize;
    integer k = m_kBlock[nth];
    while ( (k/=2) > 0 ) {
      real_type * xj = x0 + k*n;
      real_type * xp = x0;
      real_type * Bj = B0 + k*n_x_nx;
      real_type * Dj = D0 + k*n_x_n;
      real_type * Ej = E0 + k*n_x_n;
      real_type * T  = T0 + k*m_Tsize;
      integer   k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        real_type * xpp = x + jpp*n;
        alglin::gemv(
          NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 )
          gemv(
            NO_TRANSPOSE, n, m_nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1
          );
        alglin::trsv(
          UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*n;
        xp += k_x_2*n;
        T  += k_x_2*m_Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n(
    integer   nth,
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    real_type * xn = x + (nblock+1)*n + m_qx;
    integer iblock = m_iBlock[2*nth+0];
    integer eblock = m_iBlock[2*nth+1];
    real_type * x0 = x      + iblock*n;
    real_type * B0 = m_Bmat + iblock*n_x_nx;
    real_type * D0 = m_Dmat + iblock*n_x_n;
    real_type * E0 = m_Emat + iblock*n_x_n;
    real_type * T0 = m_Tmat + iblock*m_Tsize;
    integer k = m_kBlock[nth];
    while ( (k/=2) > 0 ) {
      real_type * xj = x0 + k*n;
      real_type * xp = x0;
      real_type * Bj = B0 + k*n_x_nx;
      real_type * Dj = D0 + k*n_x_n;
      real_type * Ej = E0 + k*n_x_n;
      real_type * T  = T0 + k*m_Tsize;
      integer  k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        real_type * xpp = x + jpp*n;
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Dj,  n,
                xp,  ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, m_nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        alglin::trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              swap( nrhs, xj+i, ldX, xj+P[i], ldX );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*n;
        xp += k_x_2*n;
        T  += k_x_2*m_Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_reduced( real_type x[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    real_type * xn = x + (nblock+1)*n + m_qx;
    integer k = 1;
    while ( k < m_reduced_nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer     j   = m_iBlock[jj];
        integer     jp  = m_iBlock[jj-k];
        integer     jpp = m_iBlock[std::min(jj+k,m_reduced_nblk)];
        real_type * Dj  = m_Dmat + j*n_x_n;
        real_type * Ej  = m_Emat + j*n_x_n;
        real_type * xj  = x + j*n;
        real_type * xp  = x + jp*n;
        real_type * xpp = x + jpp*n;
        alglin::gemv(
          NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 ) {
          real_type * Bj = m_Bmat + j*n_x_nx;
          gemv(
            NO_TRANSPOSE, n, m_nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1
          );
        }
        real_type const * T = m_Tmat + j*m_Tsize;
        alglin::trsv(
          UPPER, NO_TRANSPOSE, NON_UNIT,
          n, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] );
        }
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n_reduced(
    integer   nrhs,
    real_type x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    real_type * xn = x + (nblock+1)*n + m_qx;
    integer k = 1;
    while ( k < m_reduced_nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < m_reduced_nblk; jj += 2*k ) {
        integer     j   = m_iBlock[jj];
        integer     jp  = m_iBlock[jj-k];
        integer     jpp = m_iBlock[std::min(jj+k,m_reduced_nblk)];
        real_type * Dj  = m_Dmat + j*n_x_n;
        real_type * Ej  = m_Emat + j*n_x_n;
        real_type * xj  = x + j*n;
        real_type * xjp = x + jp*n;
        real_type * xpp = x + jpp*n;
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Dj,  n,
                xjp, ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 ) {
          real_type * Bj = m_Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, m_nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        real_type const * T = m_Tmat + j*m_Tsize;
        trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*n;
          for ( integer i = 0; i < n; ++i )
            if ( P[i] > i )
              alglin::swap( nrhs, xj+i, ldX, xj+P[i], ldX );
        }
      }
    }
  }

  /*\
   |             _     _ __  __
   |    __ _  __| | __| |  \/  |_   __
   |   / _` |/ _` |/ _` | |\/| \ \ / /
   |  | (_| | (_| | (_| | |  | |\ V /
   |   \__,_|\__,_|\__,_|_|  |_| \_/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::Mv( real_type const x[], real_type res[] ) const {
    alglin::zero( numRows(), res, 1 );
    addMv( 1.0, x, res );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addMv(
    real_type       alpha,
    real_type const x[],
    real_type       res[]
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    // internal blocks block
    t_Value const * D  = m_Dmat;
    t_Value const * E  = m_Emat;
    t_Value const * B  = m_Bmat;
    t_Value const * xx = x;
    t_Value const * xe = x  + nblock*n;
    t_Value const * xq = xe + n;
    t_Value const * xb = xq + m_qx;
    t_Value *       yy = res;
    for ( integer i = 0; i < nblock; ++i ) {
      alglin::gemv(
        NO_TRANSPOSE, n, n, alpha, D, n, xx, 1, 1.0, yy, 1
      );
      xx += n;
      alglin::gemv(
        NO_TRANSPOSE, n, n, alpha, E, n, xx, 1, 1.0, yy, 1
      );
      alglin::gemv(
        NO_TRANSPOSE, n, m_nx, alpha, B, n, xb, 1, 1.0, yy, 1
      );
      yy += n;
      D  += n_x_n;
      E  += n_x_n;
      B  += n_x_nx;
    }

    integer     m = n+m_qr;
    real_type * H = m_H0Nqp;
    alglin::gemv( NO_TRANSPOSE, m, n,    alpha, H, m, x,  1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( NO_TRANSPOSE, m, n,    alpha, H, m, xe, 1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( NO_TRANSPOSE, m, m_qx, alpha, H, m, xq, 1, 1.0, yy, 1 ); H += m * m_qx;
    alglin::gemv( NO_TRANSPOSE, m, m_nx, alpha, H, m, xb, 1, 1.0, yy, 1 );

    if ( m_nr > 0 ) {
      yy += m;
      alglin::gemv(
        NO_TRANSPOSE, m_nr, m_nx,
        alpha, m_Fmat[0], m_nr, xb, 1, 1.0, yy, 1
      );
      t_Value const * C = m_Cmat;
      xx = x;
      for ( integer i = 0; i <= nblock; ++i ) {
        alglin::gemv(
          NO_TRANSPOSE, m_nr, n,
          alpha, C, m_nr, xx, 1, 1.0, yy, 1
        );
        xx += n; C += nr_x_n;
      }
      alglin::gemv(
        NO_TRANSPOSE, m_nr, m_qx,
        alpha, m_Cqmat, m_nr, xx, 1, 1.0, yy, 1
      );
    }
  }

  /*\
   |  ____
   | / ___| _ __   __ _ _ __ ___  ___
   | \___ \| '_ \ / _` | '__/ __|/ _ \
   |  ___) | |_) | (_| | |  \__ \  __/
   | |____/| .__/ \__,_|_|  |___/\___|
   |       |_|
   |
  \*/

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternB(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    integer i0 = nbl*n + offs;
    integer j0 = (nblock+1)*n + m_qx + offs;
    for ( integer ij = 0; ij < n_x_nx; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesB( integer nbl, real_type V[] ) const {
    alglin::copy( n_x_nx, m_Bmat + nbl*n_x_nx, 1, V, 1 );
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternC(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    integer i0 = (nblock+1)*n + m_qr + offs;
    integer j0 = nbl*n + offs;
    for ( integer ij = 0; ij < nr_x_n; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesC( integer nbl, real_type V[] ) const {
    alglin::copy( nr_x_n, m_Cmat + nbl*nr_x_n, 1, V, 1 );
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternD(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & n = m_block_size;
    integer i0 = nbl*n + offs;
    integer j0 = nbl*n + offs;
    for ( integer ij = 0; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesD( integer nbl, real_type V[] ) const {
    alglin::copy( n_x_n, m_Dmat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternE(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & n = m_block_size;
    integer i0 = nbl*n + offs;
    integer j0 = (nbl+1)*n + offs;
    for ( integer ij = 0; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesE( integer nbl, real_type V[] ) const {
    alglin::copy( n_x_n, m_Emat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternF(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer i0 = (nblock+1)*n + m_qr + offs;
    integer j0 = (nblock+1)*n + m_qx + offs;
    for ( integer ij = 0; ij < nr_x_nx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesF( real_type V[] ) const {
    alglin::copy( nr_x_nx, m_Fmat[0], 1, V, 1 );
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternCq(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    integer i0 = (nblock+1)*n + m_qr + offs;
    integer j0 = (nblock+1)*n + offs;
    for ( integer ij = 0; ij < nr_x_qx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesCq( real_type V[] ) const {
    alglin::copy( nr_x_qx, m_Cqmat, 1, V, 1 );
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternH( integer I[], integer J[], integer offs ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer nqr = n + m_qr;
    integer nnz = nqr * m_Nc;
    integer i0 = nblock*n;
    integer j0 = i0 - n;
    for ( integer ij = 0; ij < nnz; ++ij ) {
      I[ij] = i0 + (ij % nqr) + offs;
      integer j = integer(ij/nqr); if ( j >= n ) j += j0;
      J[ij] = j + offs;
    }
    return nnz;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesH( real_type V[] ) const {
    integer const & n = m_block_size;
    integer nnz = (n + m_qr) * m_Nc;
    alglin::copy( nnz, m_H0Nqp, 1, V, 1 );
    return nnz;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparsePattern(
    integer I[],
    integer J[],
    integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;

    integer kkk = 0;
    for ( integer nbl = 0; nbl < nblock; ++nbl ) {
      kkk += this->patternD( nbl, I+kkk, J+kkk, offs );
      kkk += this->patternE( nbl, I+kkk, J+kkk, offs );
      kkk += this->patternB( nbl, I+kkk, J+kkk, offs );
    }

    // H
    kkk += this->patternH( I+kkk, J+kkk, offs );

    // C
    for ( integer nbl = 0; nbl <= nblock; ++nbl )
      kkk += this->patternC( nbl, I+kkk, J+kkk, offs );

    // F
    kkk += this->patternF( I+kkk, J+kkk, offs );

    // Cq
    kkk += this->patternCq( I+kkk, J+kkk, offs );

    UTILS_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparsePattern( I, J, offs ), inserted {} values, expected {}\n",
      kkk, sparseNnz()
    );

  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparseValues( real_type V[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer kkk = 0;
    for ( integer nbl = 0; nbl < nblock; ++nbl ) {
      kkk += this->valuesD( nbl, V+kkk );
      kkk += this->valuesE( nbl, V+kkk );
      kkk += this->valuesB( nbl, V+kkk );
    }

    // H
    kkk += this->valuesH( V+kkk );

    // C
    for ( integer nbl = 0; nbl <= nblock; ++nbl )
      kkk += this->valuesC( nbl, V+kkk );

    // F
    kkk += this->valuesF( V+kkk );

    // Cq
    kkk += this->valuesCq( V+kkk );

    UTILS_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparseValues( V ), inserted {} values, expected {}\n",
      kkk, sparseNnz()
    );

  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparseLoad(
    real_type const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz
  ) {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    integer const rH   = n*nblock;
    integer const rC   = rH + n + m_qr;
    integer const nrow = rC + m_nr;

    integer const cCq  = rH  + n;
    integer const cF   = cCq + m_qx;
    integer const ncol = cF  + m_nx;

    fillZero();

    for ( integer kkk = 0; kkk < M_nnz; ++kkk ) {
      integer   i = M_row[kkk] - r_offs;
      integer   j = M_col[kkk] - c_offs;
      real_type v = M_values[kkk];
      // cerca blocco
      bool ok = true;
      if ( i < rH ) {
        if ( j < cCq ) { // DE
          // cerca blocchi
          integer ib = i/n;
          integer jb = j/n;
          if ( ib == jb ) {
            D(ib,i%n,j%n) = v;
          } else if ( ib+1 == jb ) {
            E(ib,i%n,j%n) = v;
          } else {
            ok = false;
          }
        } else if ( j < cF ) { // Hq
          ok = false;
        } else if ( j < ncol ) { // B
          integer ib = i/n;
          B(ib,i%n,j-cF) = v;
        } else {
          ok = false;
        }
      } else if ( i < rC ) {
        if ( j < n ) { // H0
          H(i-rH,j) = v;
        } else if ( j < rH ) {
          ok = false;
        } else if ( j < ncol ) { // HN,Hq,Hp
          H(i-rH,j-rH+n) = v;
        } else {
          ok = false;
        }
      } else if ( i < nrow ) {
        if ( j < cCq ) {
          integer jb = j/n;
          C(jb,i-rC,j%n) = v;
        } else if ( j < cF ) {
          Cq(i-rC,j-cCq) =v;
        } else if ( j < ncol ) {
          F(i-rC,j-cF) = v;
        } else {
          ok = false;
        }
      } else {
        ok = false;
      }
      UTILS_ASSERT(
        ok,
        "in BorderedCR<t_Value>::sparseLoad, "
        "indices (i,j) = ( {}, {}) out of pattern!\n",
        M_row[kkk], M_col[kkk]
      );
    }
  }

  /*\
   |   ____                        _    _   _
   |  / ___| _   _ _ __   ___ _ __| |  | | | |
   |  \___ \| | | | '_ \ / _ \ '__| |  | | | |
   |   ___) | |_| | |_) |  __/ |  | |__| |_| |
   |  |____/ \__,_| .__/ \___|_|  |_____\___/
   |              |_|
  \*/

  static
  inline
  void
  Create_CompCol_Matrix(
    SuperMatrix * A,
    int           m,
    int           n,
    int           nnz,
    double      * nzval,
    int         * rowind,
    int         * colptr,
    Stype_t       stype,
    Dtype_t       dtype,
    Mtype_t       mtype
  ) {
    dCreate_CompCol_Matrix(
      A, m, n, nnz, nzval, rowind, colptr, stype, dtype, mtype
    );
  }

  static
  inline
  void
  Create_CompCol_Matrix(
    SuperMatrix * A,
    int           m,
    int           n,
    int           nnz,
    float       * nzval,
    int         * rowind,
    int         * colptr,
    Stype_t       stype,
    Dtype_t       dtype,
    Mtype_t       mtype
  ) {
    sCreate_CompCol_Matrix(
      A, m, n, nnz, nzval, rowind, colptr, stype, dtype, mtype
    );
  }

  static
  inline
  void
  Create_Dense_Matrix(
    SuperMatrix * X,
    int           m,
    int           n,
    double      * x,
    int           ldx,
    Stype_t       stype,
    Dtype_t       dtype,
    Mtype_t       mtype
  ) {
    dCreate_Dense_Matrix( X, m, n, x, ldx, stype, dtype, mtype );
  }

  static
  inline
  void
  Create_Dense_Matrix(
    SuperMatrix * X,
    int           m,
    int           n,
    float       * x,
    int           ldx,
    Stype_t       stype,
    Dtype_t       dtype,
    Mtype_t       mtype
  ) {
    sCreate_Dense_Matrix( X, m, n, x, ldx, stype, dtype, mtype );
  }

  template <typename T>
  class SuperLU {
  public:

    static
    void
    gstrf(
      superlu_options_t * options,
      SuperMatrix       * A,
      int                 relax,
      int                 panel_size,
      int               * etree,
      void              * work,
      int                 lwork,
      int               * perm_c,
      int               * perm_r,
      SuperMatrix       * L,
      SuperMatrix       * U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      GlobalLU_t        * Glu,
      #endif
      SuperLUStat_t     * stat,
      int               * info
    );

    static
    void
    gstrs(
      trans_t         trans,
      SuperMatrix   * L,
      SuperMatrix   * U,
      int           * perm_c,
      int           * perm_r,
      SuperMatrix   * B,
      SuperLUStat_t * stat,
      int           * info
    );

  };

  template <>
  void
  SuperLU<float>::gstrf(
    superlu_options_t * options,
    SuperMatrix       * A,
    int                 relax,
    int                 panel_size,
    int               * etree,
    void              * work,
    int                 lwork,
    int               * perm_c,
    int               * perm_r,
    SuperMatrix       * L,
    SuperMatrix       * U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t        * Glu,
    #endif
    SuperLUStat_t     * stat,
    int               * info
  ) {
    sgstrf(
      options, A, relax, panel_size, etree, work, lwork,
      perm_c, perm_r, L, U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      Glu,
      #endif
      stat, info
    );
  }

  template <>
  void
  SuperLU<double>::gstrf(
    superlu_options_t * options,
    SuperMatrix       * A,
    int                 relax,
    int                 panel_size,
    int               * etree,
    void              * work,
    int                 lwork,
    int               * perm_c,
    int               * perm_r,
    SuperMatrix       * L,
    SuperMatrix       * U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t        * Glu,
    #endif
    SuperLUStat_t     * stat,
    int               * info
  ) {
    dgstrf(
      options, A, relax, panel_size, etree, work, lwork,
      perm_c, perm_r, L, U,
      #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      Glu,
      #endif
      stat, info
    );
  }

  template <>
  void
  SuperLU<float>::gstrs(
    trans_t         trans,
    SuperMatrix   * L,
    SuperMatrix   * U,
    int           * perm_c,
    int           * perm_r,
    SuperMatrix   * B,
    SuperLUStat_t * stat,
    int           * info
  ) {
    sgstrs( trans, L, U, perm_c, perm_r, B, stat,  info );
  }

  template <>
  void
  SuperLU<double>::gstrs(
    trans_t         trans,
    SuperMatrix   * L,
    SuperMatrix   * U,
    int           * perm_c,
    int           * perm_r,
    SuperMatrix   * B,
    SuperLUStat_t * stat,
    int           * info
  ) {
    dgstrs( trans, L, U, perm_c, perm_r, B, stat,  info );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_SuperLU() {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    int neq = int(this->numRows());
    int nnz = nblock * ( 2 * n_x_n + m_nx * n ) +
              ( n * (nblock+1) + (m_qx+m_nx) ) * m_nr +
              ( n+m_qr ) * ( n_x_2+m_qx+m_nx );

    m_superluInteger.reallocate( size_t(nnz+5*neq+1) );
    m_superluValue.reallocate( size_t(nnz) );

    set_default_options(&m_slu_options);

    // Initialize the statistics variables.
    StatInit( &m_slu_stats );

    m_slu_perm_r = m_superluInteger( size_t(neq) ); /* row permutations from partial pivoting */
    m_slu_perm_c = m_superluInteger( size_t(neq) ); /* column permutation vector */
    m_slu_etree  = m_superluInteger( size_t(neq) );

    real_type * values = m_superluValue( size_t(nnz) );
    int       * rowind = m_superluInteger( size_t(nnz) );
    int       * colptr = m_superluInteger( size_t(neq+1) );

    // fill matrix
    int row1 = n*nblock;
    int row2 = row1 + n + m_qr;
    int kk   = 0;
    int jj   = 0;
    colptr[jj] = 0;
    for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < n; ++i ) {
        values[kk] = D( 0, i, j );
        rowind[kk] = i;
        ++kk;
      }
      for ( int i = 0; i < n + m_qr; ++i ) {
        values[kk] = H( i, j );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < m_nr; ++i ) {
        values[kk] = C( 0, i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    for ( int nbl = 1; nbl < nblock; ++nbl ) {
      int rown = nbl*n;
      for ( int j = 0; j < n; ++j ) {
        for ( int i = 0; i < n; ++i ) {
          values[kk] = E( nbl-1, i, j );
          rowind[kk] = i+rown-n;
          ++kk;
        }
        for ( int i = 0; i < n; ++i ) {
          values[kk] = D( nbl, i, j );
          rowind[kk] = i+rown;
          ++kk;
        }
        for ( int i = 0; i < m_nr; ++i ) {
          values[kk] = C( nbl, i, j );
          rowind[kk] = i+row2;
          ++kk;
        }
        colptr[++jj] = kk;
      }
    }

    for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < n; ++i ) {
        values[kk] = E( nblock-1, i, j );
        rowind[kk] = i+row1-n;
        ++kk;
      }
      for ( int i = 0; i < n+m_qr; ++i ) {
        values[kk] = H( i, j+n );;
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < m_nr; ++i ) {
        values[kk] = C( nblock, i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    for ( int j = 0; j < m_qx; ++j ) {
      for ( int i = 0; i < n+m_qr; ++i ) {
        values[kk] = H( i, j+2*n );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < m_nr; ++i ) {
        values[kk] = Cq( i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    for ( int j = 0; j < m_nx; ++j ) {
      for ( int nbl = 0; nbl < nblock; ++nbl ) {
        for ( int i = 0; i < n; ++i ) {
          values[kk] = B( nbl, i, j );
          rowind[kk] = i+nbl*n;
          ++kk;
        }
      }
      for ( int i = 0; i < n+m_qr; ++i ) {
        values[kk] = H( i, j+2*n+m_qx );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < m_nr; ++i ) {
        values[kk] = F( i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    UTILS_ASSERT(
      kk == nnz,
      "BABD_SuperLU::factorize -- dgstrf() error nnz = {} != {}\n", nnz, kk
    );

    // Create matrix A in the format expected by SuperLU.
    Create_CompCol_Matrix(
      &m_slu_A, neq, neq, nnz,
      values, rowind, colptr,
      SLU_NC, SLU_D, SLU_GE
    );

    /*\
     * Get column permutation vector perm_c[], according to permc_spec:
     *   ColPerm = 0: natural ordering
     *   ColPerm = 1: minimum degree on structure of A'*A
     *   ColPerm = 2: minimum degree on structure of A'+A
     *   ColPerm = 3: approximate minimum degree for unsymmetric matrices
    \*/
    //cout << "get_perm_c.\n";
    get_perm_c( m_slu_options.ColPerm, &m_slu_A, m_slu_perm_c );
    //cout << "sp_preorder.\n";
    sp_preorder( &m_slu_options, &m_slu_A, m_slu_perm_c, m_slu_etree, &m_slu_AC );

    int panel_size = sp_ienv(1);
    int relax      = sp_ienv(2);
    int info       = 0;
    //cout << "dgstrf.\n";
    SuperLU<t_Value>::gstrf(
      &m_slu_options, &m_slu_AC, relax, panel_size,
      m_slu_etree, nullptr, 0,
      m_slu_perm_c, m_slu_perm_r, &m_slu_L, &m_slu_U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &m_slu_glu,
    #endif
      &m_slu_stats, &info
    );

    // Free un-wanted storage
    Destroy_SuperMatrix_Store(&m_slu_A);
    Destroy_CompCol_Permuted(&m_slu_AC);
    StatFree(&m_slu_stats);

    UTILS_ASSERT(
      info == 0,
      "BABD_SuperLU::factorize -- [sd]gstrf() error returns INFO = {}\n", info
    );
    //cout << "done\n";
  }

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_SuperLU( real_type x[] ) const {
    int const   nrhs = 1;
    int         info;
    SuperMatrix B;

    trans_t trans = NOTRANS; // TRANS

    // Initialize the statistics variables.
    StatInit(&m_slu_stats) ;

    int nrow = m_slu_L.nrow;

    Create_Dense_Matrix(
      &B, nrow, nrhs,
      x, nrow,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    SuperLU<t_Value>::gstrs(
      trans, &m_slu_L, &m_slu_U, m_slu_perm_c, m_slu_perm_r, &B, &m_slu_stats, &info
    );

    Destroy_SuperMatrix_Store( &B ) ;
    StatFree(&m_slu_stats);

    return info == 0;
  }

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_SuperLU(
    integer   nrhs,
    real_type rhs[],
    integer   ldRhs
  ) const {
    int         info;
    SuperMatrix B;

    trans_t trans = NOTRANS; // TRANS

    // Initialize the statistics variables.
    StatInit(&m_slu_stats) ;

    int nrow = m_slu_L.nrow;

    Create_Dense_Matrix(
      &B, nrow, nrhs,
      rhs, ldRhs,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    SuperLU<t_Value>::gstrs(
      trans, &m_slu_L, &m_slu_U, m_slu_perm_c, m_slu_perm_r, &B, &m_slu_stats, &info
    );

    Destroy_SuperMatrix_Store( &B ) ;
    StatFree(&m_slu_stats);

    return info == 0;
  }

  // ---------------------------------------------------------------------------

  /*\
   |   __  __   _ _____ _      _   ___
   |  |  \/  | /_\_   _| |    /_\ | _ )
   |  | |\/| |/ _ \| | | |__ / _ \| _ \
   |  |_|  |_/_/ \_\_| |____/_/ \_\___/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::printMatlab( ostream_type & stream ) const {
    integer nnz = this->sparseNnz();
    Malloc<t_Value> mem(" BorderedCR::printMatlab real");
    Malloc<integer> memi(" BorderedCR::printMatlab integer");
    memi.allocate( size_t(2*nnz) );
    t_Value * V = mem.malloc( size_t(nnz) );
    integer * I = memi( size_t(nnz) );
    integer * J = memi( size_t(nnz) );

    this->sparsePattern( I, J, 1 );
    this->sparseValues( V );

    fmt::print( stream,
      "I = [ {}",
      I[0]
    );
    for ( integer i = 1; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", I[i] );
      else                 fmt::print( stream, ",  {}", I[i] );
    }
    fmt::print( stream, " ];\n\n" );

    fmt::print( stream,
      "J = [ {}",
      J[0]
    );
    for ( integer i = 1; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", J[i] );
      else                 fmt::print( stream, ",  {}", J[i] );
    }
    fmt::print( stream, " ];\n\n" );

    fmt::print( stream,
      "V = [ {}",
      V[0]
    );
    for ( integer i = 1; i < nnz; ++i ) {
      if ( (i % 20) == 0 ) fmt::print( stream, ", ...\n  {}", V[i] );
      else                 fmt::print( stream, ",  {}", V[i] );
    }
    fmt::print( stream,
      " ];\n\n"
      "nr  = {};\n"
      "nc  = {};\n\n"
      "MAT = sparse( I, J, V, nr, nc );\n",
      this->numRows(), this->numCols()
    );
  }

  // ---------------------------------------------------------------------------

  template class BorderedCR<float>;
  template class BorderedCR<double>;

}
