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
  BorderedCR<t_Value>::BorderedCR( Utils::ThreadPool * _TP )
  : m_baseValue("BorderedCR_values")
  , m_baseInteger("BorderedCR_integers")
  , m_superluValue("BorderedCR_superluValue")
  , m_superluInteger("BorderedCR_superluInteger")
  , m_number_of_blocks(0)
  , m_dim(0)
  , m_qr(0)
  , m_qx(0)
  , m_nr(0)
  , m_nx(0)
  , n_x_2(0)
  , n_x_n(0)
  , n_x_nx(0)
  , nr_x_n(0)
  , nr_x_nx(0)
  , nr_x_qx(0)
  , Nr(0)
  , Nc(0)
  , m_last_selected(BORDERED_LAST_LU)
  , m_selected(BORDERED_LU)
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
  , m_TP(_TP)
  {
    m_usedThread = m_TP == nullptr ? 1 : m_TP->size();
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
    integer _n,
    integer _qr,
    integer _qx,
    integer _nr,
    integer _nx
  ) {

    if ( m_number_of_blocks == nblock && m_dim == _n  &&
         m_qr     == _qr     && m_qx  == _qx &&
         m_nr     == _nr     && m_nx  == _nx ) return;

    m_number_of_blocks = nblock;
    m_dim    = _n;
    m_qr     = _qr;
    m_qx     = _qx;
    m_nr     = _nr;
    m_nx     = _nx;
    n_x_2    = m_dim*2;
    n_x_n    = m_dim*m_dim;
    nr_x_n   = m_dim*m_nr;
    n_x_nx   = m_dim*m_nx;
    nr_x_nx  = m_nr*m_nx;
    nr_x_qx  = m_nr*m_qx;
    Nr       = n_x_2+m_nr+m_qr;
    Nc       = n_x_2+m_nx+m_qx;
    Tsize    = 2*n_x_n+m_dim;

    integer N = std::max(Nr,Nc);
    m_Lwork = std::max(N,2*std::max(n_x_n,std::max(nr_x_n,n_x_nx)));

    valueType tmp; // get optimal allocation
    integer info = alglin::geqrf( Nr, Nc, nullptr, Nr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    if ( m_Lwork < integer(tmp) ) m_Lwork = integer(tmp);

    info = geqp3( Nr, Nc, nullptr, Nr, nullptr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( m_Lwork < integer(tmp) ) m_Lwork = integer(tmp);

    m_LworkT = 2*n_x_n;
    info = alglin::geqrf( n_x_2, m_dim, nullptr, n_x_2, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    m_LworkQR = integer(tmp);

    info = alglin::geqp3( n_x_2, m_dim, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 );
    UTILS_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( m_LworkQR < integer(tmp) ) m_LworkQR = integer(tmp);

    integer nnz = nblock*(n_x_nx+2*n_x_n+Tsize+m_dim) +
                  (nblock+1)*nr_x_n +
                  nr_x_qx +
                  Nr*Nc + (m_dim+m_qr)*Nc +
                  m_Lwork + (m_LworkT+m_LworkQR+nr_x_nx+m_nr)*m_usedThread;
    integer innz = nblock*m_dim + (3+m_dim)*m_usedThread;

    m_baseValue.allocate(size_t(nnz));
    m_baseInteger.allocate(size_t(innz));

    m_Bmat  = m_baseValue( size_t(nblock*n_x_nx) );
    m_Cmat  = m_baseValue( size_t((nblock+1)*nr_x_n) );
    m_Cqmat = m_baseValue( size_t(nr_x_qx) );
    m_Dmat  = m_baseValue( size_t(nblock*n_x_n) );
    m_Emat  = m_baseValue( size_t(nblock*n_x_n) );

    m_Tmat  = m_baseValue( size_t(nblock*Tsize) );
    m_Ttau  = m_baseValue( size_t(nblock*m_dim) );
    m_Hmat  = m_baseValue( size_t(Nr*Nc) );
    m_H0Nqp = m_baseValue( size_t((m_dim+m_qr)*Nc) );

    m_Work = m_baseValue( size_t(m_Lwork) );
    m_Perm = m_baseInteger( size_t(nblock*m_dim) );
    iBlock = m_baseInteger( size_t(2*m_usedThread) );
    kBlock = m_baseInteger( size_t(m_usedThread) );

    m_xb_thread.resize( size_t(m_usedThread) );
    m_perm_thread.resize( size_t(m_usedThread) );
    m_WorkT.resize( size_t(m_usedThread) );
    m_WorkQR.resize( size_t(m_usedThread) );
    m_Fmat.resize( size_t(m_usedThread) );

    // precompute partition for parallel computation
    for ( size_t nt = 0; nt < size_t(m_usedThread); ++nt ) {
      m_perm_thread[nt] = m_baseInteger( size_t(m_dim) );
      m_xb_thread[nt]   = m_baseValue( size_t(m_nr) );
      m_WorkT[nt]       = m_baseValue( size_t(m_LworkT) );
      m_WorkQR[nt]      = m_baseValue( size_t(m_LworkQR) );
      m_Fmat[nt]        = m_baseValue( size_t(nr_x_nx) );
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
    m_usedThread = M.m_usedThread;
    m_TP         = M.m_TP;
    allocate( M.m_number_of_blocks, M.m_dim, M.m_qr, M.m_qx, M.m_nr, M.m_nx );
    alglin::copy( m_number_of_blocks*n_x_nx,     M.m_Bmat,    1, m_Bmat,    1 );
    alglin::copy( (m_number_of_blocks+1)*nr_x_n, M.m_Cmat,    1, m_Cmat,    1 );
    alglin::copy( m_number_of_blocks*n_x_n,      M.m_Dmat,    1, m_Dmat,    1 );
    alglin::copy( m_number_of_blocks*n_x_n,      M.m_Emat,    1, m_Emat,    1 );
    alglin::copy( nr_x_nx,             M.m_Fmat[0], 1, m_Fmat[0], 1 );
    alglin::copy( (m_dim+m_qr)*Nc,     M.m_H0Nqp,   1, m_H0Nqp,   1 );
    alglin::copy( nr_x_qx,             M.m_Cqmat,   1, m_Cqmat,   1 );
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
      m_dim, m_qr, m_qx, m_nr, m_nx
    );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroD()
  { alglin::zero( m_number_of_blocks*n_x_n, m_Dmat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroE()
  { alglin::zero( m_number_of_blocks*n_x_n, m_Emat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroB()
  { alglin::zero( m_number_of_blocks*n_x_nx, m_Bmat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroF()
  { alglin::zero( nr_x_nx, m_Fmat[0], 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroH()
  { alglin::zero( (m_dim+m_qr)*Nc, m_H0Nqp, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroC()
  { alglin::zero( (m_number_of_blocks+1)*nr_x_n, m_Cmat, 1 ); }

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
    valueType const H0[], integer ld0,
    valueType const HN[], integer ldN,
    valueType const Hq[], integer ldQ,
    valueType const Hp[], integer ldP
  ) {
    // (n+qr) x ( n + n + qx + nx )
    integer     m = m_dim + m_qr;
    valueType * H = m_H0Nqp;
    alglin::gecopy( m, m_dim, H0, ld0, H, m ); H += m * m_dim;
    alglin::gecopy( m, m_dim, HN, ldN, H, m ); H += m * m_dim;
    alglin::gecopy( m, m_qx,  Hq, ldQ, H, m ); H += m * m_qx;
    alglin::gecopy( m, m_nx,  Hp, ldP, H, m );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom(
    MatrixWrapper<valueType> const & H0,
    MatrixWrapper<valueType> const & HN,
    MatrixWrapper<valueType> const & Hq,
    MatrixWrapper<valueType> const & Hp
  ) {

    integer m = m_dim + m_qr;

    UTILS_ASSERT(
      H0.numRows() == m && H0.numCols() == m_dim,
      "loadBottom, bad dimension size(H0) = {} x {} expected {} x {}\n",
      H0.numRows(), H0.numCols(), m, m_dim
    );

    UTILS_ASSERT(
      HN.numRows() == m && HN.numCols() == m_dim,
      "loadBottom, bad dimension size(HN) = {} x {} expected {} x {}\n",
      HN.numRows(), HN.numCols(), m, m_dim
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
  BorderedCR<t_Value>::loadBottom( valueType const _H0Nqp[], integer ldH ) {
    integer nq = m_dim + m_qr;
    alglin::gecopy( nq, Nc, _H0Nqp, ldH, m_H0Nqp, nq );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom( MatrixWrapper<valueType> const & H ) {
    integer m = m_dim + m_qr;
    UTILS_ASSERT(
      H.numRows() == m && H.numCols() == Nc,
      "loadBottom, bad dimension size(H) = {} x {} expected {} x {}\n",
      H.numRows(), H.numCols(), m, Nc
    );
    loadBottom( H.data(), H.lDim() );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    valueType const C0[], integer ld0,
    valueType const CN[], integer ldN,
    valueType const Cq[], integer ldCq,
    valueType const F[],  integer ldF
  ) {
    // (n+qr) x ( n + n + qx + nx )
    gecopy( m_nr, m_dim,  C0, ld0,  m_Cmat,                          m_nr );
    gecopy( m_nr, m_dim,  CN, ldN,  m_Cmat+m_number_of_blocks*n_x_n, m_nr );
    gecopy( m_nr, m_qx,   Cq, ldCq, m_Cqmat,                         m_nr );
    gecopy( m_nr, m_nx,   F,  ldF,  m_Fmat[0],                       m_nr );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    MatrixWrapper<valueType> const & C0,
    MatrixWrapper<valueType> const & CN,
    MatrixWrapper<valueType> const & Cq,
    MatrixWrapper<valueType> const & F
  ) {

    UTILS_ASSERT(
      C0.numRows() == m_nr && C0.numCols() == m_dim,
      "loadBottom2, bad dimension size(C0) = {} x {} expected {} x {}\n",
      C0.numRows(), C0.numCols(), m_nr, m_dim
    );

    UTILS_ASSERT(
      CN.numRows() == m_nr && CN.numCols() == m_dim,
      "loadBottom2, bad dimension size(CN) = {} x {} expected {} x {}\n",
      CN.numRows(), CN.numCols(), m_nr, m_dim
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
  BorderedCR<t_Value>::loadBottom2( MatrixWrapper<valueType> const & H ) {
    integer const & nblock = m_number_of_blocks;
    valueType const * ptr = H.data();
    integer           ld  = H.lDim();
    this->loadC(  0,      ptr, ld ); ptr += nr_x_n;
    this->loadC(  nblock, ptr, ld ); ptr += nr_x_n;
    this->loadCq( ptr, ld ); ptr += m_nr * m_qx;
    this->loadF(  ptr, ld );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB( integer nbl, valueType const B[], integer ldB )
  { alglin::gecopy( m_dim, m_nx, B, ldB, m_Bmat + nbl*n_x_nx, m_dim ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    UTILS_ASSERT(
      B.numRows() == m_dim && B.numCols() == m_nx,
      "loadB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), m_dim, m_nx
    );
    alglin::gecopy( m_dim, m_nx, B.data(), B.lDim(), m_Bmat + nbl*n_x_nx, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB( integer nbl, valueType const B[], integer ldB ) {
    valueType * BB = m_Bmat + nbl*n_x_nx;
    alglin::geadd( m_dim, m_nx, 1.0, B, ldB, 1.0, BB, m_dim, BB, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    UTILS_ASSERT(
      B.numRows() == m_dim && B.numCols() == m_nx,
      "addtoB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), m_dim, m_nx
    );
    valueType * BB = m_Bmat + nbl*n_x_nx;
    alglin::geadd( m_dim, m_nx, 1.0, B.data(), B.lDim(), 1.0, BB, m_dim, BB, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC( integer nbl, valueType const C[], integer ldC )
  { alglin::gecopy( m_nr, m_dim, C, ldC, m_Cmat + nbl*nr_x_n, m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == m_dim,
      "loadC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, m_dim
    );
    alglin::gecopy( m_nr, m_dim, C.data(), C.lDim(), m_Cmat + nbl*nr_x_n, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC( integer nbl, valueType const C[], integer ldC ) {
    UTILS_ASSERT(
      ldC >= m_nr, "addtoC( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    valueType * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, m_dim, 1.0, C, ldC, 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == m_dim,
      "addtoC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, m_dim
    );
    valueType * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, m_dim, 1.0, C.data(), C.lDim(), 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2( integer nbl, valueType const C[], integer ldC ) {
    UTILS_ASSERT(
      ldC >= m_nr, "addtoC2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    valueType * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n_x_2, 1.0, C, ldC, 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    UTILS_ASSERT(
      C.numRows() == m_nr && C.numCols() == n_x_2,
      "addtoC2( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), m_nr, n_x_2
    );
    valueType * CC = m_Cmat + nbl*nr_x_n;
    alglin::geadd( m_nr, n_x_2, 1.0, C.data(), C.lDim(), 1.0, CC, m_nr, CC, m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD( integer nbl, valueType const D[], integer ldD )
  { alglin::gecopy( m_dim, m_dim, D, ldD, m_Dmat + nbl*n_x_n, m_dim ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD(
    integer                          nbl,
    MatrixWrapper<valueType> const & D
  ) {
    UTILS_ASSERT(
      D.numRows() == m_dim && D.numCols() == m_dim,
      "loadD( {}, D) bad dimension size(D) = {} x {} expected {} x {}\n",
      nbl, D.numRows(), D.numCols(), m_dim, m_dim
    );
    alglin::gecopy( m_dim, m_dim, D.data(), D.lDim(), m_Dmat + nbl*n_x_n, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE( integer nbl, valueType const E[], integer ldE )
  { alglin::gecopy( m_dim, m_dim, E, ldE, m_Emat + nbl*n_x_n, m_dim ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE(
    integer                          nbl,
    MatrixWrapper<valueType> const & E
  ) {
    UTILS_ASSERT(
      E.numRows() == m_dim && E.numCols() == m_dim,
      "loadE( {}, E) bad dimension size(E) = {} x {} expected {} x {}\n",
      nbl, E.numRows(), E.numCols(), m_dim, m_dim
    );
    alglin::gecopy( m_dim, m_dim, E.data(), E.lDim(), m_Emat + nbl*n_x_n, m_dim );
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDE( integer nbl, valueType const DE[], integer ldDE ) {
    alglin::gecopy( m_dim, m_dim, DE, ldDE, m_Dmat + nbl*n_x_n, m_dim ); DE += m_dim*ldDE;
    alglin::gecopy( m_dim, m_dim, DE, ldDE, m_Emat + nbl*n_x_n, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDEB( integer nbl, valueType const DEB[], integer ldDEB ) {
    alglin::gecopy( m_dim, m_dim, DEB, ldDEB, m_Dmat + nbl*n_x_n,  m_dim ); DEB += m_dim*ldDEB;
    alglin::gecopy( m_dim, m_dim, DEB, ldDEB, m_Emat + nbl*n_x_n,  m_dim ); DEB += m_dim*ldDEB;
    alglin::gecopy( m_dim, m_nx,  DEB, ldDEB, m_Bmat + nbl*n_x_nx, m_dim );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( valueType const F[], integer ldF )
  { alglin::gecopy( m_nr, m_nx, F, ldF, m_Fmat[0], m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( MatrixWrapper<valueType> const & F ) {
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
  BorderedCR<t_Value>::addtoF( valueType const F[], integer ldF )
  { alglin::gecopy( m_nr, m_nx, F, ldF, m_Fmat[0], m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( MatrixWrapper<valueType> const & F ) {
    UTILS_ASSERT(
      F.numRows() == m_nr && F.numCols() == m_nx,
      "addtoF(F) bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), m_nr, m_nx
    );
    alglin::geadd( m_nr, m_nx, 1.0, F.data(), F.lDim(), 1.0, m_Fmat[0], m_nr, m_Fmat[0], m_nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( valueType const Cq[], integer ldC )
  { alglin::gecopy( m_nr, m_qx, Cq, ldC, m_Cqmat, m_nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( MatrixWrapper<valueType> const & Cq ) {
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
  BorderedCR<t_Value>::loadCqF( valueType const CqF[], integer ldCF ) {
    alglin::gecopy( m_nr, m_qx, CqF, ldCF, m_Cqmat, m_nr ); CqF += m_qx*ldCF;
    alglin::gecopy( m_nr, m_nx, CqF, ldCF, m_Fmat[0], m_nr );
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

    if ( m_usedThread > 1 ) {
      iBlock[0] = 0;
      iBlock[1] = static_cast<integer>(nblock/m_usedThread);
      for ( integer nt = 1; nt < m_usedThread; ++nt ) {
        iBlock[2*nt+0]  = iBlock[2*nt-1]+1;
        iBlock[2*nt+1]  = static_cast<integer>(((nt+1)*nblock)/m_usedThread);
      }
      // launch thread
      for ( integer nt = 1; nt < m_usedThread; ++nt ) {
        // fill zero F(...) only that of the extra threads
        if ( nr_x_nx > 0 ) alglin::zero( nr_x_nx, m_Fmat[size_t(nt)], 1 );
        m_TP->run( nt, &BorderedCR<t_Value>::factorize_block, this, nt );
      }
      factorize_block(0);
      // wait thread
      if ( m_usedThread > 1 ) { m_TP->wait_all();
        // accumulate F(...)
        if ( nr_x_nx > 0 )
          for ( integer nt = 1; nt < m_usedThread; ++nt )
            alglin::axpy( nr_x_nx, 1.0, m_Fmat[size_t(nt)], 1, m_Fmat[0], 1 );
      }
      factorize_reduced();
    } else {
      iBlock[0] = 0;
      iBlock[1] = nblock;
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
    valueType const TOP[],
    valueType const BOTTOM[],
    valueType       T[],
    integer         iperm[]
  ) const {
    alglin::gecopy( m_dim, m_dim, TOP,    m_dim, T,       n_x_2 );
    alglin::gecopy( m_dim, m_dim, BOTTOM, m_dim, T+m_dim, n_x_2 );
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      info = alglin::getrf( n_x_2, m_dim, T, n_x_2, iperm );
      break;
    case BORDERED_QR:
      info = alglin::geqrf( n_x_2, m_dim, T, n_x_2, T+2*n_x_n, m_WorkQR[size_t(nth)], m_LworkQR );
      break;
    case BORDERED_QRP:
      { integer * P = m_perm_thread[size_t(nth)];
        std::fill_n( P, m_dim, 0 );
        info = alglin::geqp3(
          n_x_2, m_dim, T, n_x_2, P, T+2*n_x_n, m_WorkQR[size_t(nth)], m_LworkQR
        );
        if ( info == 0 ) permutation_to_exchange( m_dim, P, iperm );
      }
      break;
    case BORDERED_SUPERLU:
      UTILS_ERROR0( "BorderedCR::buildT, cannot be used with SUPERLU\n" );
      break;
    }
    UTILS_ASSERT( info == 0, "BorderedCR::factorize INFO = {}\n", info );
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
    valueType const T[],
    integer   const iperm[],
    valueType       TOP[],
    integer         ldTOP,
    valueType       BOTTOM[],
    integer         ldBOTTOM,
    integer         ncol
  ) const {
    valueType * W = m_WorkT[size_t(nth)];
    alglin::gecopy( m_dim, ncol, TOP,    ldTOP,    W,       n_x_2 );
    alglin::gecopy( m_dim, ncol, BOTTOM, ldBOTTOM, W+m_dim, n_x_2 );
    // Apply row interchanges to the right hand sides.
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      info = alglin::swaps( ncol, W, n_x_2, 0, m_dim-1, iperm, 1 );
      UTILS_ASSERT( info == 0, "BorderedCR::applyT INFO = {}\n", info );
      alglin::trsm(
        LEFT, LOWER, NO_TRANSPOSE, UNIT,
        m_dim, ncol, 1.0, T, n_x_2, W, n_x_2
      );
      alglin::gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        m_dim, ncol, m_dim,
        -1.0, T+m_dim, n_x_2,  // M
              W,       n_x_2,  // L^(-1) TOP
         1.0, W+m_dim, n_x_2   // TOP = BOTTOM - M L^(-1) TOP
      );
      break;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = alglin::ormqr(
        LEFT, TRANSPOSE,
        n_x_2, ncol, // righe x colonne
        m_dim,       // numero riflettori usati nel prodotto Q
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
    alglin::gecopy( m_dim, ncol, W+m_dim, n_x_2, TOP,    ldTOP    );
    alglin::gecopy( m_dim, ncol, W,       n_x_2, BOTTOM, ldBOTTOM );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT(
    integer         nth,
    valueType const T[],
    integer   const iperm[],
    valueType       TOP[],
    valueType       BOTTOM[]
  ) const {
    valueType * W = m_WorkT[size_t(nth)];
    size_t nn = size_t(m_dim)*sizeof(valueType);
    memcpy( W,       TOP,    nn );
    memcpy( W+m_dim, BOTTOM, nn );
    //copy( n, TOP,    1, W,   1 );
    //copy( n, BOTTOM, 1, W+n, 1 );
    integer info = 0;
    switch ( m_selected ) {
    case BORDERED_LU:
      // Apply row interchanges to the right hand sides.
      info = swaps( 1, W, n_x_2, 0, m_dim-1, iperm, 1 );
      if ( info == 0 ) {
        alglin::trsv( LOWER, NO_TRANSPOSE, UNIT, m_dim, T, n_x_2, W, 1 );
        alglin::gemv(
          NO_TRANSPOSE,
          m_dim, m_dim,
          -1.0, T+m_dim, n_x_2,
                W,       1,
           1.0, W+m_dim, 1
        );
      }
      break;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = alglin::ormqr(
        LEFT, TRANSPOSE,
        n_x_2, 1, // righe x colonne
        m_dim,    // numero riflettori usati nel prodotto Q
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
    memcpy( TOP,    W+m_dim, nn );
    memcpy( BOTTOM, W,       nn );
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

    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;

    valueType * Bmat0 = m_Bmat + iblock*n_x_nx;
    valueType * Cmat0 = m_Cmat + iblock*nr_x_n;
    valueType * Dmat0 = m_Dmat + iblock*n_x_n;
    valueType * Emat0 = m_Emat + iblock*n_x_n;
    valueType * T0    = m_Tmat + iblock*Tsize;
    integer   * P0    = m_Perm + iblock*m_dim;

    valueType * Fmat_th = m_Fmat[size_t(nth)];

    integer k = 1;
    while ( k < nblk ) {

      valueType * Bjp = Bmat0;
      valueType * Bj  = Bmat0 + k*n_x_nx;
      valueType * Cjp = Cmat0;
      valueType * Cj  = Cmat0 + k*nr_x_n;
      valueType * Djp = Dmat0;
      valueType * Dj  = Dmat0 + k*n_x_n;
      valueType * Ejp = Emat0;
      valueType * Ej  = Emat0 + k*n_x_n;
      valueType * T   = T0    + k*Tsize;
      integer   * P   = P0    + k*m_dim;

      // -----------------------------------------------------------------------
      integer k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {

        buildT( nth, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( nth, T, P, Djp, m_dim, Dj, m_dim, m_dim );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( nth, T, P, Ejp, m_dim, Ej, m_dim, m_dim );

        if ( m_nx > 0 ) applyT( nth, T, P, Bjp, m_dim, Bj, m_dim, m_nx );

        if ( m_nr > 0 ) {

          if ( m_selected == BORDERED_QRP ) {
            integer i = m_dim;
            do {
              --i;
              if ( P[i] > i ) swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            m_nr, m_dim, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_dim, m_dim,
            -1.0, Cj,  m_nr,
                  Dj,  m_dim,
             1.0, Cjp, m_nr
          );

          integer     jpp = std::min(j+k,eblock);
          valueType * Cpp = m_Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_dim, m_dim,
            -1.0, Cj,  m_nr,
                  Ej,  m_dim,
             1.0, Cpp, m_nr
          );
        }

        if ( nr_x_nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_nx, m_dim,
            -1.0, Cj,      m_nr,
                  Bj,      m_dim,
             1.0, Fmat_th, m_nr // solo accumulo!
          );

        // NEXT STEP
        T   += k_x_2*Tsize;
        P   += k_x_2*m_dim;
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
    kBlock[nth] = k;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_reduced() {
    integer nblk = 2*m_usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      // -----------------------------------------------------------------------
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType * T   = m_Tmat + j*Tsize;
        integer   * P   = m_Perm + j*m_dim;
        valueType * Djp = m_Dmat + jp*n_x_n;
        valueType * Dj  = m_Dmat + j*n_x_n;
        valueType * Ejp = m_Emat + jp*n_x_n;
        valueType * Ej  = m_Emat + j*n_x_n;

        buildT( 0, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( 0, T, P, Djp, m_dim, Dj, m_dim, m_dim );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( 0, T, P, Ejp, m_dim, Ej, m_dim, m_dim );

        if ( m_nx > 0 ) {
          valueType * Bj  = m_Bmat + j*n_x_nx;
          valueType * Bjp = m_Bmat + jp*n_x_nx;
          applyT( 0, T, P, Bjp, m_dim, Bj, m_dim, m_nx );
        }

        if ( m_nr > 0 ) {
          valueType * Cj  = m_Cmat + j*nr_x_n;
          valueType * Cjp = m_Cmat + jp*nr_x_n;
          if ( m_selected == BORDERED_QRP ) {
            integer i = m_dim;
            do {
              --i;
              if ( P[i] > i ) swap( m_nr, Cj+i*m_nr, 1, Cj+P[i]*m_nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            m_nr, m_dim, 1.0, T, n_x_2, Cj, m_nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_dim, m_dim,
            -1.0, Cj,  m_nr,
                  Dj,  m_dim,
             1.0, Cjp, m_nr
          );

          integer     jpp = iBlock[std::min(jj+k,nblk)];
          valueType * Cpp = m_Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_dim, m_dim,
            -1.0, Cj,  m_nr,
                  Ej,  m_dim,
             1.0, Cpp, m_nr
          );

        }

        if ( nr_x_nx > 0 ) {
          valueType * Cj = m_Cmat + j*nr_x_n;
          valueType * Bj = m_Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, m_nx, m_dim,
            -1.0, Cj,        m_nr,
                  Bj,        m_dim,
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
    valueType * Cnb = m_Cmat + nblock*nr_x_n;
    valueType * W0  = m_Hmat;
    valueType * WN  = W0+m_dim*Nr;
    valueType * Wq  = WN+m_dim*Nr;
    valueType * Wp  = Wq+m_qx*Nr;

    alglin::gecopy( m_dim,  m_dim,  m_Dmat, m_dim, W0, Nr );
    alglin::gecopy( m_dim,  m_dim,  m_Emat, m_dim, WN, Nr );
    alglin::gezero( m_dim,  m_qx,                  Wq, Nr );
    alglin::gecopy( m_dim,  m_nx,   m_Bmat, m_dim, Wp, Nr );

    alglin::gecopy( m_dim+m_qr, Nc, m_H0Nqp, m_dim+m_qr, m_Hmat+m_dim, Nr );

    integer offs = n_x_2+m_qr;

    alglin::gecopy( m_nr, m_dim, m_Cmat,    m_nr, W0+offs, Nr );
    alglin::gecopy( m_nr, m_dim, Cnb,       m_nr, WN+offs, Nr );
    alglin::gecopy( m_nr, m_qx,  m_Cqmat,   m_nr, Wq+offs, Nr );
    alglin::gecopy( m_nr, m_nx,  m_Fmat[0], m_nr, Wp+offs, Nr );

    bool ok = false;
    switch ( m_last_selected ) {
    case BORDERED_LAST_LU:
      ok = last_lu.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_LUPQ:
      ok = last_lupq.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_QR:
      ok = last_qr.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_QRP:
      ok = last_qrp.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_SVD:
      ok = last_svd.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_LSS:
      ok = last_lss.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_LSY:
      ok = last_lsy.factorize( Nr, Nc, m_Hmat, Nr );
      break;
    case BORDERED_LAST_PINV:
      ok = last_pinv.factorize( Nr, Nc, m_Hmat, Nr );
      break;
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
  BorderedCR<t_Value>::solve_last( valueType x[] ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * X = x + (nblock-1)*m_dim;
    swap( m_dim, X, 1, x, 1 ); // uso x stesso come temporaneo
    bool ok = false;
    switch ( m_last_selected ) {
    case BORDERED_LAST_LU:   ok = last_lu.solve( X );   break;
    case BORDERED_LAST_LUPQ: ok = last_lupq.solve( X ); break;
    case BORDERED_LAST_QR:   ok = last_qr.solve( X );   break;
    case BORDERED_LAST_QRP:  ok = last_qrp.solve( X );  break;
    case BORDERED_LAST_SVD:  ok = last_svd.solve( X );  break;
    case BORDERED_LAST_LSS:  ok = last_lss.solve( X );  break;
    case BORDERED_LAST_LSY:  ok = last_lsy.solve( X );  break;
    case BORDERED_LAST_PINV: ok = last_pinv.solve( X ); break;
    }
    if ( ok ) swap( m_dim, X, 1, x, 1 );
    return ok;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_last(
    integer   nrhs,
    valueType x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * X = x + (nblock-1)*m_dim;
    for ( integer i = 0; i < nrhs; ++i ) swap( m_dim, X+i*ldX, 1, x+i*ldX, 1 );
    bool ok = false;
    switch ( m_last_selected ) {
    case BORDERED_LAST_LU:   ok = last_lu.solve( nrhs, X, ldX );   break;
    case BORDERED_LAST_LUPQ: ok = last_lupq.solve( nrhs, X, ldX ); break;
    case BORDERED_LAST_QR:   ok = last_qr.solve( nrhs, X, ldX );   break;
    case BORDERED_LAST_QRP:  ok = last_qrp.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_SVD:  ok = last_svd.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_LSS:  ok = last_lss.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_LSY:  ok = last_lsy.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_PINV: ok = last_pinv.solve( nrhs, X, ldX ); break;
    }
    if ( ok )
      for ( integer i = 0; i < nrhs; ++i )
        swap( m_dim, X+i*ldX, 1, x+i*ldX, 1 );
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
  BorderedCR<t_Value>::solve_CR( valueType x[] ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xb = x + (nblock+1)*m_dim + m_qr; // deve essere b!
    if ( m_usedThread > 1 ) {
      if ( m_nr > 0 ) {
        for ( integer nt = 1; nt < m_usedThread; ++nt ) {
          alglin::zero( m_nr, m_xb_thread[size_t(nt)], 1 );
          m_TP->run(
            nt, &BorderedCR<t_Value>::forward, this,
            nt, x, m_xb_thread[size_t(nt)]
          );
        }
        alglin::zero( m_nr, m_xb_thread[0], 1 );
        forward(0,x,xb);
        m_TP->wait_all();
        for ( integer nt = 1; nt < m_usedThread; ++nt )
          axpy( m_nr, 1.0, m_xb_thread[size_t(nt)], 1, xb, 1 );
      } else {
        for ( integer nt = 1; nt < m_usedThread; ++nt )
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

    if ( m_usedThread > 1 ) {
      backward_reduced(x);
      for ( integer nt = 1; nt < m_usedThread; ++nt )
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
    valueType rhs[],
    integer   ldRhs
  ) const {
    if ( m_usedThread > 1 ) {
      for ( integer nt = 1; nt < m_usedThread; ++nt )
        m_TP->run( nt, &BorderedCR<t_Value>::forward_n, this, nt, nrhs, rhs, ldRhs );
      forward_n( 0, nrhs, rhs, ldRhs );
      m_TP->wait_all();
      forward_n_reduced( nrhs, rhs, ldRhs );
    } else {
      forward_n( 0, nrhs, rhs, ldRhs );
    }

    if ( !solve_last( nrhs, rhs, ldRhs ) ) return false;

    if ( m_usedThread > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < m_usedThread; ++nt )
        m_TP->run( nt, &BorderedCR<t_Value>::backward_n, this, nt, nrhs, rhs, ldRhs );
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
    valueType x[],
    valueType xb[]
  ) const {
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;
    valueType * x0 = x      + iblock*m_dim;
    valueType * T0 = m_Tmat + iblock*Tsize;
    integer   * P0 = m_Perm + iblock*m_dim;
    valueType * C0 = m_Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      valueType * xj  = x0 + k*m_dim;
      valueType * xjp = x0;
      valueType * T   = T0 + k*Tsize;
      integer   * P   = P0 + k*m_dim;
      valueType * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;
      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, xj );
        if ( m_nr > 0 )
          gemv(
            NO_TRANSPOSE,
            m_nr, m_dim, -1.0,
            Cj, m_nr, xj, 1,
            1.0, xb, 1
          ); // solo accumulato
        xj  += k_x_2*m_dim;
        xjp += k_x_2*m_dim;
        T   += k_x_2*Tsize;
        P   += k_x_2*m_dim;
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
    valueType x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xb = x + (nblock+1)*m_dim + m_qr;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;
    valueType * x0 = x      + iblock*m_dim;
    valueType * T0 = m_Tmat + iblock*Tsize;
    integer   * P0 = m_Perm + iblock*m_dim;
    valueType * C0 = m_Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      valueType * xj  = x0 + k*m_dim;
      valueType * xjp = x0;
      valueType * T   = T0 + k*Tsize;
      integer   * P   = P0 + k*m_dim;
      valueType * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;

      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( m_nr > 0 ) {
          m_spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, nrhs, m_dim,
            -1.0, Cj, m_nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          m_spin.unlock();
        }
        xj  += k_x_2*m_dim;
        xjp += k_x_2*m_dim;
        T   += k_x_2*Tsize;
        P   += k_x_2*m_dim;
        Cj  += k_x_2*nr_x_n;
      }
      k *= 2;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_reduced(
    valueType x[],
    valueType xb[]
  ) const {
    integer nblk = 2*m_usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType const * T   = m_Tmat + j*Tsize;
        integer   const * P   = m_Perm + j*m_dim;
        valueType       * xj  = x + j*m_dim;
        valueType       * xjp = x + jp*m_dim;
        applyT( 0, T, P, xjp, xj );
        if ( m_nr > 0 ) {
          valueType * Cj = m_Cmat + j*nr_x_n;
          alglin::gemv(
            NO_TRANSPOSE, m_nr, m_dim,
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
    valueType x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xb = x + (nblock+1)*m_dim + m_qr;
    integer nblk = 2*m_usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType const * T   = m_Tmat + j*Tsize;
        integer   const * P   = m_Perm + j*m_dim;
        valueType       * xj  = x + j*m_dim;
        valueType       * xjp = x + jp*m_dim;
        applyT( 0, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( m_nr > 0 ) {
          valueType * Cj = m_Cmat + j*nr_x_n;
          m_spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_nr, nrhs, m_dim,
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
  BorderedCR<t_Value>::backward( integer nth, valueType x[] ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xn = x + (nblock+1)*m_dim + m_qx;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    valueType * x0 = x      + iblock*m_dim;
    valueType * B0 = m_Bmat + iblock*n_x_nx;
    valueType * D0 = m_Dmat + iblock*n_x_n;
    valueType * E0 = m_Emat + iblock*n_x_n;
    valueType * T0 = m_Tmat + iblock*Tsize;
    integer k = kBlock[nth];
    while ( (k/=2) > 0 ) {
      valueType * xj = x0 + k*m_dim;
      valueType * xp = x0;
      valueType * Bj = B0 + k*n_x_nx;
      valueType * Dj = D0 + k*n_x_n;
      valueType * Ej = E0 + k*n_x_n;
      valueType * T  = T0 + k*Tsize;
      integer   k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        valueType * xpp = x + jpp*m_dim;
        alglin::gemv(
          NO_TRANSPOSE, m_dim, m_dim,
          -1.0, Dj, m_dim, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          NO_TRANSPOSE, m_dim, m_dim,
          -1.0, Ej, m_dim, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 )
          gemv(
            NO_TRANSPOSE, m_dim, m_nx,
            -1.0, Bj, m_dim, xn, 1, 1.0, xj, 1
          );
        alglin::trsv(
          UPPER, NO_TRANSPOSE, NON_UNIT,
          m_dim, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*m_dim;
          for ( integer i = 0; i < m_dim; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*m_dim;
        xp += k_x_2*m_dim;
        T  += k_x_2*Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n(
    integer   nth,
    integer   nrhs,
    valueType x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xn = x + (nblock+1)*m_dim + m_qx;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    valueType * x0 = x      + iblock*m_dim;
    valueType * B0 = m_Bmat + iblock*n_x_nx;
    valueType * D0 = m_Dmat + iblock*n_x_n;
    valueType * E0 = m_Emat + iblock*n_x_n;
    valueType * T0 = m_Tmat + iblock*Tsize;
    integer k = kBlock[nth];
    while ( (k/=2) > 0 ) {
      valueType * xj = x0 + k*m_dim;
      valueType * xp = x0;
      valueType * Bj = B0 + k*n_x_nx;
      valueType * Dj = D0 + k*n_x_n;
      valueType * Ej = E0 + k*n_x_n;
      valueType * T  = T0 + k*Tsize;
      integer  k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        valueType * xpp = x + jpp*m_dim;
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          m_dim, nrhs, m_dim,
          -1.0, Dj,  m_dim,
                xp,  ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          m_dim, nrhs, m_dim,
          -1.0, Ej,  m_dim,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_dim, nrhs, m_nx,
            -1.0, Bj, m_dim,
                  xn, ldX,
             1.0, xj, ldX
          );
        alglin::trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          m_dim, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*m_dim;
          for ( integer i = 0; i < m_dim; ++i )
            if ( P[i] > i )
              swap( nrhs, xj+i, ldX, xj+P[i], ldX );
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx;
        Dj += k_x_2*n_x_n;
        Ej += k_x_2*n_x_n;
        xj += k_x_2*m_dim;
        xp += k_x_2*m_dim;
        T  += k_x_2*Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_reduced( valueType x[] ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xn = x + (nblock+1)*m_dim + m_qx;
    integer nblk = 2*m_usedThread-1;
    integer k = 1;
    while ( k < nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer     j   = iBlock[jj];
        integer     jp  = iBlock[jj-k];
        integer     jpp = iBlock[std::min(jj+k,nblk)];
        valueType * Dj  = m_Dmat + j*n_x_n;
        valueType * Ej  = m_Emat + j*n_x_n;
        valueType * xj  = x + j*m_dim;
        valueType * xp  = x + jp*m_dim;
        valueType * xpp = x + jpp*m_dim;
        alglin::gemv(
          NO_TRANSPOSE, m_dim, m_dim,
          -1.0, Dj, m_dim, xp,  1, 1.0, xj, 1
        );
        alglin::gemv(
          NO_TRANSPOSE, m_dim, m_dim,
          -1.0, Ej, m_dim, xpp, 1, 1.0, xj, 1
        );
        if ( m_nx > 0 ) {
          valueType * Bj = m_Bmat + j*n_x_nx;
          gemv(
            NO_TRANSPOSE, m_dim, m_nx,
            -1.0, Bj, m_dim, xn, 1, 1.0, xj, 1
          );
        }
        valueType const * T = m_Tmat + j*Tsize;
        alglin::trsv(
          UPPER, NO_TRANSPOSE, NON_UNIT,
          m_dim, T, n_x_2, xj, 1
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*m_dim;
          for ( integer i = 0; i < m_dim; ++i )
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
    valueType x[],
    integer   ldX
  ) const {
    integer const & nblock = m_number_of_blocks;
    valueType * xn = x + (nblock+1)*m_dim + m_qx;
    integer nblk = 2*m_usedThread-1;
    integer k = 1;
    while ( k < nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer     j   = iBlock[jj];
        integer     jp  = iBlock[jj-k];
        integer     jpp = iBlock[std::min(jj+k,nblk)];
        valueType * Dj  = m_Dmat + j*n_x_n;
        valueType * Ej  = m_Emat + j*n_x_n;
        valueType * xj  = x + j*m_dim;
        valueType * xjp = x + jp*m_dim;
        valueType * xpp = x + jpp*m_dim;
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          m_dim, nrhs, m_dim,
          -1.0, Dj,  m_dim,
                xjp, ldX,
           1.0, xj,  ldX
        );
        alglin::gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          m_dim, nrhs, m_dim,
          -1.0, Ej,  m_dim,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( m_nx > 0 ) {
          valueType * Bj = m_Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            m_dim, nrhs, m_nx,
            -1.0, Bj, m_dim,
                  xn, ldX,
             1.0, xj, ldX
          );
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        valueType const * T = m_Tmat + j*Tsize;
        trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          m_dim, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( m_selected == BORDERED_QRP ) {
          integer const * P = m_Perm + j*m_dim;
          for ( integer i = 0; i < m_dim; ++i )
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
  BorderedCR<t_Value>::Mv( valueType const x[], valueType res[] ) const {
    alglin::zero( numRows(), res, 1 );
    addMv( x, res );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addMv( valueType const x[], valueType res[] ) const {
    integer const & nblock = m_number_of_blocks;
    // internal blocks block
    t_Value const * D  = m_Dmat;
    t_Value const * E  = m_Emat;
    t_Value const * B  = m_Bmat;
    t_Value const * xx = x;
    t_Value const * xe = x  + nblock*m_dim;
    t_Value const * xq = xe + m_dim;
    t_Value const * xb = xq + m_qx;
    t_Value *       yy = res;
    for ( integer i = 0; i < nblock; ++i ) {
      alglin::gemv(
        NO_TRANSPOSE, m_dim, m_dim,
        1.0, D, m_dim, xx, 1, 1.0, yy, 1
      );
      xx += m_dim;
      alglin::gemv(
        NO_TRANSPOSE, m_dim, m_dim,
        1.0, E, m_dim, xx, 1, 1.0, yy, 1
      );
      alglin::gemv(
        NO_TRANSPOSE, m_dim, m_nx,
        1.0, B, m_dim, xb, 1, 1.0, yy, 1
      );
      yy += m_dim;
      D  += n_x_n;
      E  += n_x_n;
      B  += n_x_nx;
    }

    integer     m = m_dim+m_qr;
    valueType * H = m_H0Nqp;
    alglin::gemv( NO_TRANSPOSE, m, m_dim, 1.0, H, m, x,  1, 1.0, yy, 1 ); H += m * m_dim;
    alglin::gemv( NO_TRANSPOSE, m, m_dim, 1.0, H, m, xe, 1, 1.0, yy, 1 ); H += m * m_dim;
    alglin::gemv( NO_TRANSPOSE, m, m_qx,  1.0, H, m, xq, 1, 1.0, yy, 1 ); H += m * m_qx;
    alglin::gemv( NO_TRANSPOSE, m, m_nx,  1.0, H, m, xb, 1, 1.0, yy, 1 );

    if ( m_nr > 0 ) {
      yy += m;
      alglin::gemv(
        NO_TRANSPOSE, m_nr, m_nx,
        1.0, m_Fmat[0], m_nr, xb, 1, 1.0, yy, 1
      );
      t_Value const * C = m_Cmat;
      xx = x;
      for ( integer i = 0; i <= nblock; ++i ) {
        alglin::gemv(
          NO_TRANSPOSE, m_nr, m_dim,
          1.0, C, m_nr, xx, 1, 1.0, yy, 1
        );
        xx += m_dim; C += nr_x_n;
      }
      alglin::gemv(
        NO_TRANSPOSE, m_nr, m_qx,
        1.0, m_Cqmat, m_nr, xx, 1, 1.0, yy, 1
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
    integer i0 = nbl*m_dim + offs;
    integer j0 = (nblock+1)*m_dim + m_qx + offs;
    for ( integer ij = 0; ij < n_x_nx; ++ij ) {
      I[ij] = i0 + (ij % m_dim);
      J[ij] = j0 + integer(ij/m_dim);
    }
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesB( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_nx, m_Bmat + nbl*n_x_nx, 1, V, 1 );
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternC(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer i0 = (nblock+1)*m_dim + m_qr + offs;
    integer j0 = nbl*m_dim + offs;
    for ( integer ij = 0; ij < nr_x_n; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesC( integer nbl, valueType V[] ) const {
    alglin::copy( nr_x_n, m_Cmat + nbl*nr_x_n, 1, V, 1 );
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternD(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer i0 = nbl*m_dim + offs;
    integer j0 = nbl*m_dim + offs;
    for ( integer ij = 0; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % m_dim);
      J[ij] = j0 + integer(ij/m_dim);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesD( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_n, m_Dmat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternE(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer i0 = nbl*m_dim + offs;
    integer j0 = (nbl+1)*m_dim + offs;
    for ( integer ij = 0; ij < n_x_n; ++ij ) {
      I[ij] = i0 + (ij % m_dim);
      J[ij] = j0 + integer(ij/m_dim);
    }
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesE( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_n, m_Emat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternF(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer i0 = (nblock+1)*m_dim + m_qr + offs;
    integer j0 = (nblock+1)*m_dim + m_qx + offs;
    for ( integer ij = 0; ij < nr_x_nx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesF( valueType V[] ) const {
    alglin::copy( nr_x_nx, m_Fmat[0], 1, V, 1 );
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternCq(
    integer I[], integer J[], integer offs
  ) const {
    integer const & nblock = m_number_of_blocks;
    integer i0 = (nblock+1)*m_dim + m_qr + offs;
    integer j0 = (nblock+1)*m_dim + offs;
    for ( integer ij = 0; ij < nr_x_qx; ++ij ) {
      I[ij] = i0 + (ij % m_nr);
      J[ij] = j0 + integer(ij/m_nr);
    }
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesCq( valueType V[] ) const {
    alglin::copy( nr_x_qx, m_Cqmat, 1, V, 1 );
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternH( integer I[], integer J[], integer offs ) const {
    integer const & nblock = m_number_of_blocks;
    integer nqr = m_dim + m_qr;
    integer nnz = nqr * Nc;
    integer i0 = nblock*m_dim + offs;
    integer j0 = i0 - m_dim;
    for ( integer ij = 0; ij < nnz; ++ij ) {
      I[ij] = i0 + (ij % nqr);
      integer j = integer(ij/nqr); if ( j >= m_dim ) j += j0;
      J[ij] = j;
    }
    return nnz;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesH( valueType V[] ) const {
    integer nnz = (m_dim + m_qr) * Nc;
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
  BorderedCR<t_Value>::sparseValues( valueType V[] ) const {
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
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz
  ) {
    integer const & nblock = m_number_of_blocks;
    integer const rH   = m_dim*nblock;
    integer const rC   = rH + m_dim+m_qr;
    integer const nrow = rC + m_nr;

    integer const cCq  = rH  + m_dim;
    integer const cF   = cCq + m_qx;
    integer const ncol = cF  + m_nx;

    fillZero();

    for ( integer kkk = 0; kkk < M_nnz; ++kkk ) {
      integer   i = M_row[kkk] - r_offs;
      integer   j = M_col[kkk] - c_offs;
      valueType v = M_values[kkk];
      // cerca blocco
      bool ok = true;
      if ( i < rH ) {
        if ( j < cCq ) { // DE
          // cerca blocchi
          integer ib = i/m_dim;
          integer jb = j/m_dim;
          if ( ib == jb ) {
            D(ib,i%m_dim,j%m_dim) = v;
          } else if ( ib+1 == jb ) {
            E(ib,i%m_dim,j%m_dim) = v;
          } else {
            ok = false;
          }
        } else if ( j < cF ) { // Hq
          ok = false;
        } else if ( j < ncol ) { // B
          integer ib = i/m_dim;
          B(ib,i%m_dim,j-cF) = v;
        } else {
          ok = false;
        }
      } else if ( i < rC ) {
        if ( j < m_dim ) { // H0
          H(i-rH,j) = v;
        } else if ( j < rH ) {
          ok = false;
        } else if ( j < ncol ) { // HN,Hq,Hp
          H(i-rH,j-rH+m_dim) = v;
        } else {
          ok = false;
        }
      } else if ( i < nrow ) {
        if ( j < cCq ) {
          integer jb = j/m_dim;
          C(jb,i-rC,j%m_dim) = v;
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

    int neq = int(this->numRows());
    int nnz = nblock * ( 2 * n_x_n + m_nx * m_dim ) +
              ( m_dim * (nblock+1) + (m_qx+m_nx) ) * m_nr +
              ( m_dim+m_qr ) * ( n_x_2+m_qx+m_nx );

    m_superluInteger.allocate( size_t(nnz+5*neq+1) );
    m_superluValue.allocate( size_t(nnz) );

    set_default_options(&slu_options);

    // Initialize the statistics variables.
    StatInit( &slu_stats );

    slu_perm_r = m_superluInteger( size_t(neq) ); /* row permutations from partial pivoting */
    slu_perm_c = m_superluInteger( size_t(neq) ); /* column permutation vector */
    slu_etree  = m_superluInteger( size_t(neq) );

    valueType * values = m_superluValue( size_t(nnz) );
    int       * rowind = m_superluInteger( size_t(nnz) );
    int       * colptr = m_superluInteger( size_t(neq+1) );

    // fill matrix
    int row1 = m_dim*nblock;
    int row2 = row1+m_dim+m_qr;
    int kk   = 0;
    int jj   = 0;
    colptr[jj] = 0;
    for ( int j = 0; j < m_dim; ++j ) {
      for ( int i = 0; i < m_dim; ++i ) {
        values[kk] = D( 0, i, j );
        rowind[kk] = i;
        ++kk;
      }
      for ( int i = 0; i < m_dim+m_qr; ++i ) {
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
      int rown = nbl*m_dim;
      for ( int j = 0; j < m_dim; ++j ) {
        for ( int i = 0; i < m_dim; ++i ) {
          values[kk] = E( nbl-1, i, j );
          rowind[kk] = i+rown-m_dim;
          ++kk;
        }
        for ( int i = 0; i < m_dim; ++i ) {
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

    for ( int j = 0; j < m_dim; ++j ) {
      for ( int i = 0; i < m_dim; ++i ) {
        values[kk] = E( nblock-1, i, j );
        rowind[kk] = i+row1-m_dim;
        ++kk;
      }
      for ( int i = 0; i < m_dim+m_qr; ++i ) {
        values[kk] = H( i, j+m_dim );;
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
      for ( int i = 0; i < m_dim+m_qr; ++i ) {
        values[kk] = H( i, j+2*m_dim );
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
        for ( int i = 0; i < m_dim; ++i ) {
          values[kk] = B( nbl, i, j );
          rowind[kk] = i+nbl*m_dim;
          ++kk;
        }
      }
      for ( int i = 0; i < m_dim+m_qr; ++i ) {
        values[kk] = H( i, j+2*m_dim+m_qx );
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
      &slu_A, neq, neq, nnz,
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
    get_perm_c( slu_options.ColPerm, &slu_A, slu_perm_c );
    //cout << "sp_preorder.\n";
    sp_preorder( &slu_options, &slu_A, slu_perm_c, slu_etree, &slu_AC );

    int panel_size = sp_ienv(1);
    int relax      = sp_ienv(2);
    int info = 0;
    //cout << "dgstrf.\n";
    SuperLU<t_Value>::gstrf(
      &slu_options, &slu_AC, relax, panel_size,
      slu_etree, nullptr, 0,
      slu_perm_c, slu_perm_r, &slu_L, &slu_U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &slu_glu,
    #endif
      &slu_stats, &info
    );

    // Free un-wanted storage
    Destroy_SuperMatrix_Store(&slu_A);
    Destroy_CompCol_Permuted(&slu_AC);
    StatFree(&slu_stats);

    UTILS_ASSERT(
      info == 0,
      "BABD_SuperLU::factorize -- [sd]gstrf() error returns INFO = {}\n", info
    );
    //cout << "done\n";

  }

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_SuperLU( valueType x[] ) const {
    int const   nrhs = 1;
    int         info;
    SuperMatrix B;

    trans_t trans = NOTRANS; // TRANS

    // Initialize the statistics variables.
    StatInit(&slu_stats) ;

    int nrow = slu_L.nrow;

    Create_Dense_Matrix(
      &B, nrow, nrhs,
      x, nrow,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    SuperLU<t_Value>::gstrs(
      trans, &slu_L, &slu_U, slu_perm_c, slu_perm_r, &B, &slu_stats, &info
    );

    Destroy_SuperMatrix_Store( &B ) ;
    StatFree(&slu_stats);

    return info == 0;
  }

  template <typename t_Value>
  bool
  BorderedCR<t_Value>::solve_SuperLU(
    integer   nrhs,
    valueType rhs[],
    integer   ldRhs
  ) const {
    int         info;
    SuperMatrix B;

    trans_t trans = NOTRANS; // TRANS

    // Initialize the statistics variables.
    StatInit(&slu_stats) ;

    int nrow = slu_L.nrow;

    Create_Dense_Matrix(
      &B, nrow, nrhs,
      rhs, ldRhs,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    SuperLU<t_Value>::gstrs(
      trans, &slu_L, &slu_U, slu_perm_c, slu_perm_r, &B, &slu_stats, &info
    );

    Destroy_SuperMatrix_Store( &B ) ;
    StatFree(&slu_stats);

    return info == 0;

  }

  // ---------------------------------------------------------------------------

  template class BorderedCR<float>;
  template class BorderedCR<double>;

}
