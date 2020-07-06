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

#include "BABD_BorderedCR.hh"
#include "Alglin_Eigen.hh"

namespace alglin {

  template <typename t_Value>
  BorderedCR<t_Value>::BorderedCR( ThreadPool * _TP )
  : baseValue("BorderedCR_values")
  , baseInteger("BorderedCR_integers")
  , superluValue("BorderedCR_superluValue")
  , superluInteger("BorderedCR_superluInteger")
  , nblock(0)
  , n(0)
  , qr(0)
  , qx(0)
  , nr(0)
  , nx(0)
  , n_x_2(0)
  , n_x_n(0)
  , n_x_nx(0)
  , nr_x_n(0)
  , nr_x_nx(0)
  , nr_x_qx(0)
  , Nr(0)
  , Nc(0)
  , last_selected(BORDERED_LAST_LU)
  , selected(BORDERED_LU)
  , H0Nqp(nullptr)
  , Bmat(nullptr)
  , Cmat(nullptr)
  , Cqmat(nullptr)
  , Dmat(nullptr)
  , Emat(nullptr)
  , Tmat(nullptr)
  , Ttau(nullptr)
  , Work(nullptr)
  , Perm(nullptr)
  , Lwork(0)
  , Hmat(nullptr)
  , pTP(_TP)
  {
    usedThread = pTP == nullptr ? 1 : pTP->size();
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
    integer _nblock,
    integer _n,
    integer _qr,
    integer _qx,
    integer _nr,
    integer _nx
  ) {

    if ( nblock == _nblock && n  == _n  &&
         qr     == _qr     && qx == _qx &&
         nr     == _nr     && nx == _nx ) return;

    nblock  = _nblock;
    n       = _n;
    qr      = _qr;
    qx      = _qx;
    nr      = _nr;
    nx      = _nx;
    n_x_2   = n*2;
    n_x_n   = n*n;
    nr_x_n  = n*nr;
    n_x_nx  = n*nx;
    nr_x_nx = nr*nx;
    nr_x_qx = nr*qx;
    Nr      = n_x_2+nr+qr;
    Nc      = n_x_2+nx+qx;
    Tsize   = 2*n_x_n+n;

    integer N = std::max(Nr,Nc);
    Lwork = std::max(N,2*std::max(n_x_n,std::max(nr_x_n,n_x_nx)));

    valueType tmp; // get optimal allocation
    integer info = alglin::geqrf( Nr, Nc, nullptr, Nr, nullptr, &tmp, -1 );
    LW_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp);

    info = geqp3( Nr, Nc, nullptr, Nr, nullptr, nullptr, &tmp, -1 );
    LW_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp);

    LworkT = 2*n_x_n;
    info = alglin::geqrf( n_x_2, n, nullptr, n_x_2, nullptr, &tmp, -1 );
    LW_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = {}\n", info
    );
    LworkQR = integer(tmp);

    info = alglin::geqp3( n_x_2, n, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 );
    LW_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = {}\n", info
    );
    if ( LworkQR < integer(tmp) ) LworkQR = integer(tmp);

    integer nnz = nblock*(n_x_nx+2*n_x_n+Tsize+n) +
                  (nblock+1)*nr_x_n +
                  nr_x_qx +
                  Nr*Nc + (n+qr)*Nc +
                  Lwork + (LworkT+LworkQR+nr_x_nx+nr)*usedThread;
    integer innz = nblock*n + (3+n)*usedThread;

    baseValue.allocate(size_t(nnz));
    baseInteger.allocate(size_t(innz));

    Bmat  = baseValue( size_t(nblock*n_x_nx) );
    Cmat  = baseValue( size_t((nblock+1)*nr_x_n) );
    Cqmat = baseValue( size_t(nr_x_qx) );
    Dmat  = baseValue( size_t(nblock*n_x_n) );
    Emat  = baseValue( size_t(nblock*n_x_n) );

    Tmat  = baseValue( size_t(nblock*Tsize) );
    Ttau  = baseValue( size_t(nblock*n) );
    Hmat  = baseValue( size_t(Nr*Nc) );
    H0Nqp = baseValue( size_t((n+qr)*Nc) );

    Work   = baseValue( size_t(Lwork) );
    Perm   = baseInteger( size_t(nblock*n) );
    iBlock = baseInteger( size_t(2*usedThread) );
    kBlock = baseInteger( size_t(usedThread) );

    xb_thread   . resize( size_t(usedThread) );
    perm_thread . resize( size_t(usedThread) );
    WorkT       . resize( size_t(usedThread) );
    WorkQR      . resize( size_t(usedThread) );
    Fmat        . resize( size_t(usedThread) );

    // precompute partition for parallel computation
    for ( size_t nt = 0; nt < size_t(usedThread); ++nt ) {
      perm_thread[nt] = baseInteger( size_t(n) );
      xb_thread[nt]   = baseValue( size_t(nr) );
      WorkT[nt]       = baseValue( size_t(LworkT) );
      WorkQR[nt]      = baseValue( size_t(LworkQR) );
      Fmat[nt]        = baseValue( size_t(nr_x_nx) );
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
    usedThread = M.usedThread;
    pTP        = M.pTP;
    allocate( M.nblock, M.n, M.qr, M.qx, M.nr, M.nx );
    alglin::copy( nblock*n_x_nx,     M.Bmat,    1, Bmat,    1 );
    alglin::copy( (nblock+1)*nr_x_n, M.Cmat,    1, Cmat,    1 );
    alglin::copy( nblock*n_x_n,      M.Dmat,    1, Dmat,    1 );
    alglin::copy( nblock*n_x_n,      M.Emat,    1, Emat,    1 );
    alglin::copy( nr_x_nx,           M.Fmat[0], 1, Fmat[0], 1 );
    alglin::copy( (n+qr)*Nc,         M.H0Nqp,   1, H0Nqp,   1 );
    alglin::copy( nr_x_qx,           M.Cqmat,   1, Cqmat,   1 );
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
    stream
      << "rows   = " << numRows() << '\n'
      << "cols   = " << numCols() << '\n'
      << "nblock = " << nblock << '\n'
      << "n      = " << n << '\n'
      << "qr     = " << qr << '\n'
      << "qx     = " << qx << '\n'
      << "nr     = " << nr << '\n'
      << "nx     = " << nx << '\n';
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroD()
  { alglin::zero( nblock*n_x_n, Dmat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroE()
  { alglin::zero( nblock*n_x_n, Emat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroB()
  { alglin::zero( nblock*n_x_nx, Bmat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroF()
  { alglin::zero( nr_x_nx, Fmat[0], 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroH()
  { alglin::zero( (n+qr)*Nc, H0Nqp, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroC()
  { alglin::zero( (nblock+1)*nr_x_n, Cmat, 1 ); }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::zeroCq()
  { alglin::zero( nr_x_qx, Cqmat, 1 ); }

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
    integer     m = n + qr;
    valueType * H = H0Nqp;
    alglin::gecopy( m, n,  H0, ld0, H, m ); H += m * n;
    alglin::gecopy( m, n,  HN, ldN, H, m ); H += m * n;
    alglin::gecopy( m, qx, Hq, ldQ, H, m ); H += m * qx;
    alglin::gecopy( m, nx, Hp, ldP, H, m );
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

    integer m = n + qr;

    LW_ASSERT(
      H0.numRows() == m && H0.numCols() == n,
      "loadBottom, bad dimension size(H0) = {} x {} expected {} x {}\n",
      H0.numRows(), H0.numCols(), m, n
    );

    LW_ASSERT(
      HN.numRows() == m && HN.numCols() == n,
      "loadBottom, bad dimension size(HN) = {} x {} expected {} x {}\n",
      HN.numRows(), HN.numCols(), m, n
    );

    LW_ASSERT(
      Hq.numRows() == m && Hq.numCols() == qx,
      "loadBottom, bad dimension size(Hq) = {} x {} expected {} x {}\n",
      Hq.numRows(), Hq.numCols(), m, qx
    );

    LW_ASSERT(
      Hp.numRows() == m && Hp.numCols() == nx,
      "loadBottom, bad dimension size(Hp) = {} x {} expected {} x {}\n",
      Hp.numRows(), Hp.numCols(), m, nx
    );

    loadBottom(
      H0.get_data(), H0.lDim(),
      HN.get_data(), HN.lDim(),
      Hq.get_data(), Hq.lDim(),
      Hp.get_data(), Hp.lDim()
    );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom( valueType const _H0Nqp[], integer ldH ) {
    integer nq = n+qr;
    alglin::gecopy( nq, Nc, _H0Nqp, ldH, H0Nqp, nq );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom( MatrixWrapper<valueType> const & H ) {
    integer m = n + qr;
    LW_ASSERT(
      H.numRows() == m && H.numCols() == Nc,
      "loadBottom, bad dimension size(H) = {} x {} expected {} x {}\n",
      H.numRows(), H.numCols(), m, Nc
    );
    loadBottom( H.get_data(), H.lDim() );
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
    gecopy( nr, n,  C0, ld0,  this->Cmat,              nr );
    gecopy( nr, n,  CN, ldN,  this->Cmat+nblock*n_x_n, nr );
    gecopy( nr, qx, Cq, ldCq, this->Cqmat,             nr );
    gecopy( nr, nx, F,  ldF,  this->Fmat[0],           nr );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    MatrixWrapper<valueType> const & C0,
    MatrixWrapper<valueType> const & CN,
    MatrixWrapper<valueType> const & Cq,
    MatrixWrapper<valueType> const & F
  ) {

    LW_ASSERT(
      C0.numRows() == nr && C0.numCols() == n,
      "loadBottom2, bad dimension size(C0) = {} x {} expected {} x {}\n",
      C0.numRows(), C0.numCols(), nr, n
    );

    LW_ASSERT(
      CN.numRows() == nr && CN.numCols() == n,
      "loadBottom2, bad dimension size(CN) = {} x {} expected {} x {}\n",
      CN.numRows(), CN.numCols(), nr, n
    );

    LW_ASSERT(
      Cq.numRows() == nr && Cq.numCols() == qx,
      "loadBottom2, bad dimension size(Cq) = {} x {} expected {} x {}\n",
      Cq.numRows(), Cq.numCols(), nr, qx
    );

    LW_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "loadBottom2, bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), nr, nx
    );

    loadBottom2(
      C0.get_data(), C0.lDim(),
      CN.get_data(), CN.lDim(),
      Cq.get_data(), Cq.lDim(),
      F.get_data(),  F.lDim()
    );

  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2( MatrixWrapper<valueType> const & H ) {
    valueType const * ptr = H.get_data();
    integer           ld  = H.lDim();
    this->loadC(  0,      ptr, ld ); ptr += nr_x_n;
    this->loadC(  nblock, ptr, ld ); ptr += nr_x_n;
    this->loadCq( ptr, ld ); ptr += nr * qx;
    this->loadF(  ptr, ld );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB( integer nbl, valueType const B[], integer ldB )
  { alglin::gecopy( n, nx, B, ldB, Bmat + nbl*n_x_nx, n ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    LW_ASSERT(
      B.numRows() == n && B.numCols() == nx,
      "loadB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), n, nx
    );
    alglin::gecopy( n, nx, B.get_data(), B.lDim(), Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB( integer nbl, valueType const B[], integer ldB ) {
    valueType * BB = Bmat + nbl*n_x_nx;
    alglin::geadd( n, nx, 1.0, B, ldB, 1.0, BB, n, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    LW_ASSERT(
      B.numRows() == n && B.numCols() == nx,
      "addtoB( {}, B) bad dimension size(B) = {} x {} expected {} x {}\n",
      nbl, B.numRows(), B.numCols(), n, nx
    );
    valueType * BB = Bmat + nbl*n_x_nx;
    alglin::geadd( n, nx, 1.0, B.get_data(), B.lDim(), 1.0, BB, n, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC( integer nbl, valueType const C[], integer ldC )
  { alglin::gecopy( nr, n, C, ldC, Cmat + nbl*nr_x_n, nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    LW_ASSERT(
      C.numRows() == nr && C.numCols() == n,
      "loadC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), nr, n
    );
    alglin::gecopy( nr, n, C.get_data(), C.lDim(), Cmat + nbl*nr_x_n, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC( integer nbl, valueType const C[], integer ldC ) {
    LW_ASSERT(
      ldC >= nr, "addtoC( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    alglin::geadd( nr, n, 1.0, C, ldC, 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    LW_ASSERT(
      C.numRows() == nr && C.numCols() == n,
      "addtoC( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), nr, n
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    alglin::geadd( nr, n, 1.0, C.get_data(), C.lDim(), 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2( integer nbl, valueType const C[], integer ldC ) {
    LW_ASSERT(
      ldC >= nr, "addtoC2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    alglin::geadd( nr, n_x_2, 1.0, C, ldC, 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    LW_ASSERT(
      C.numRows() == nr && C.numCols() == n_x_2,
      "addtoC2( {}, C) bad dimension size(C) = {} x {} expected {} x {}\n",
      nbl, C.numRows(), C.numCols(), nr, n_x_2
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    alglin::geadd( nr, n_x_2, 1.0, C.get_data(), C.lDim(), 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD( integer nbl, valueType const D[], integer ldD )
  { alglin::gecopy( n, n, D, ldD, Dmat + nbl*n_x_n, n ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD(
    integer                          nbl,
    MatrixWrapper<valueType> const & D
  ) {
    LW_ASSERT(
      D.numRows() == n && D.numCols() == n,
      "loadD( {}, D) bad dimension size(D) = {} x {} expected {} x {}\n",
      nbl, D.numRows(), D.numCols(), n, n
    );
    alglin::gecopy( n, n, D.get_data(), D.lDim(), Dmat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE( integer nbl, valueType const E[], integer ldE )
  { alglin::gecopy( n, n, E, ldE, Emat + nbl*n_x_n, n ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE(
    integer                          nbl,
    MatrixWrapper<valueType> const & E
  ) {
    LW_ASSERT(
      E.numRows() == n && E.numCols() == n,
      "loadE( {}, E) bad dimension size(E) = {} x {} expected {} x {}\n",
      nbl, E.numRows(), E.numCols(), n, n
    );
    alglin::gecopy( n, n, E.get_data(), E.lDim(), Emat + nbl*n_x_n, n );
  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDE( integer nbl, valueType const DE[], integer ldDE ) {
    alglin::gecopy( n, n, DE, ldDE, Dmat + nbl*n_x_n, n ); DE += n*ldDE;
    alglin::gecopy( n, n, DE, ldDE, Emat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadDEB( integer nbl, valueType const DEB[], integer ldDEB ) {
    alglin::gecopy( n, n,  DEB, ldDEB, Dmat + nbl*n_x_n,  n ); DEB += n*ldDEB;
    alglin::gecopy( n, n,  DEB, ldDEB, Emat + nbl*n_x_n,  n ); DEB += n*ldDEB;
    alglin::gecopy( n, nx, DEB, ldDEB, Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( valueType const F[], integer ldF )
  { alglin::gecopy( nr, nx, F, ldF, Fmat[0], nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( MatrixWrapper<valueType> const & F ) {
    LW_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "loadF(F) bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), nr, nx
    );
    alglin::gecopy( nr, nx, F.get_data(), F.lDim(), Fmat[0], nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( valueType const F[], integer ldF )
  { alglin::gecopy( nr, nx, F, ldF, Fmat[0], nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( MatrixWrapper<valueType> const & F ) {
    LW_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "addtoF(F) bad dimension size(F) = {} x {} expected {} x {}\n",
      F.numRows(), F.numCols(), nr, nx
    );
    alglin::geadd( nr, nx, 1.0, F.get_data(), F.lDim(), 1.0, Fmat[0], nr, Fmat[0], nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( valueType const Cq[], integer ldC )
  { alglin::gecopy( nr, qx, Cq, ldC, Cqmat, nr ); }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( MatrixWrapper<valueType> const & Cq ) {
    LW_ASSERT(
      Cq.numRows() == nr && Cq.numCols() == qx,
      "loadCq(Cq) bad dimension size(Cq) = {} x {} expected {} x {}\n",
      Cq.numRows(), Cq.numCols(), nr, qx
    );
    alglin::gecopy( nr, qx, Cq.get_data(), Cq.lDim(), Cqmat, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCqF( valueType const CqF[], integer ldCF ) {
    alglin::gecopy( nr, qx, CqF, ldCF, Cqmat, nr ); CqF += qx*ldCF;
    alglin::gecopy( nr, nx, CqF, ldCF, Fmat[0], nr );
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
    if ( usedThread > 1 ) {
      iBlock[0] = 0;
      iBlock[1] = static_cast<integer>(nblock/usedThread);
      for ( integer nt = 1; nt < usedThread; ++nt ) {
        iBlock[2*nt+0]  = iBlock[2*nt-1]+1;
        iBlock[2*nt+1]  = static_cast<integer>(((nt+1)*nblock)/usedThread);
      }
      // launch thread
      for ( integer nt = 1; nt < usedThread; ++nt ) {
        // fill zero F(...) only that of the extra threads
        if ( nr_x_nx > 0 ) alglin::zero( nr_x_nx, Fmat[size_t(nt)], 1 );
        pTP->run( nt, &BorderedCR<t_Value>::factorize_block, this, nt );
      }
      factorize_block(0);
      // wait thread
      pTP->wait_all();
      // accumulate F(...)
      if ( nr_x_nx > 0 )
        for ( integer nt = 1; nt < usedThread; ++nt )
          alglin::axpy( nr_x_nx, 1.0, Fmat[size_t(nt)], 1, Fmat[0], 1 );
      factorize_reduced();
    } else {
      iBlock[0] = 0;
      iBlock[1] = nblock;
      factorize_block(0);
    }
    LW_ASSERT0(
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
    alglin::gecopy( n, n, TOP,    n, T,   n_x_2 );
    alglin::gecopy( n, n, BOTTOM, n, T+n, n_x_2 );
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      info = alglin::getrf( n_x_2, n, T, n_x_2, iperm );
      break;
    case BORDERED_QR:
      info = alglin::geqrf( n_x_2, n, T, n_x_2, T+2*n_x_n, WorkQR[size_t(nth)], LworkQR );
      break;
    case BORDERED_QRP:
      { integer * P = perm_thread[size_t(nth)];
        std::fill( P, P+n, 0 );
        info = alglin::geqp3( n_x_2, n, T, n_x_2, P, T+2*n_x_n, WorkQR[size_t(nth)], LworkQR );
        if ( info == 0 ) permutation_to_exchange( n, P, iperm );
      }
      break;
    case BORDERED_SUPERLU:
      LW_ERROR0( "BorderedCR::buildT, cannot be used with SUPERLU\n" );
      break;
    }
    LW_ASSERT( info == 0, "BorderedCR::factorize INFO = {}\n", info );
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
    valueType * W = WorkT[size_t(nth)];
    alglin::gecopy( n, ncol, TOP,    ldTOP,    W,   n_x_2 );
    alglin::gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, n_x_2 );
    // Apply row interchanges to the right hand sides.
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      info = alglin::swaps( ncol, W, n_x_2, 0, n-1, iperm, 1 );
      LW_ASSERT( info == 0, "BorderedCR::applyT INFO = {}\n", info );
      alglin::trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, n_x_2, W, n_x_2 );
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
        WorkQR[size_t(nth)], LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      LW_ERROR0( "BorderedCR::applyT, cannot be used with SUPERLU\n" );
      break;
    }
    alglin::gecopy( n, ncol, W+n, n_x_2, TOP,    ldTOP );
    alglin::gecopy( n, ncol, W,   n_x_2, BOTTOM, ldBOTTOM );
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
    valueType * W = WorkT[size_t(nth)];
    size_t nn = size_t(n)*sizeof(valueType);
    memcpy( W,   TOP,    nn );
    memcpy( W+n, BOTTOM, nn );
    //copy( n, TOP,    1, W,   1 );
    //copy( n, BOTTOM, 1, W+n, 1 );
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      // Apply row interchanges to the right hand sides.
      info = swaps( 1, W, n_x_2, 0, n-1, iperm, 1 );
      if ( info == 0 ) {
        alglin::trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, n_x_2, W, 1 );
        alglin::gemv(
          NO_TRANSPOSE,
          n, n,
          -1.0, T+n, n_x_2,
                W,   1,
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
        WorkQR[size_t(nth)], LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      LW_ERROR0( "BorderedCR::applyT, cannot be used with SUPERLU\n" );
      break;
    }
    LW_ASSERT( info == 0, "BorderedCR::applyT INFO = {}\n", info );
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

    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;

    valueType * Bmat0 = Bmat + iblock*n_x_nx;
    valueType * Cmat0 = Cmat + iblock*nr_x_n;
    valueType * Dmat0 = Dmat + iblock*n_x_n;
    valueType * Emat0 = Emat + iblock*n_x_n;
    valueType * T0    = Tmat + iblock*Tsize;
    integer   * P0    = Perm + iblock*n;

    valueType * Fmat_th = Fmat[size_t(nth)];

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
      integer   * P   = P0    + k*n;

      // -----------------------------------------------------------------------
      integer k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {

        buildT( nth, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( nth, T, P, Djp, n, Dj, n, n );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( nth, T, P, Ejp, n, Ej, n, n );

        if ( nx > 0 ) applyT( nth, T, P, Bjp, n, Bj, n, nx );

        if ( nr > 0 ) {

          if ( selected == BORDERED_QRP ) {
            integer i = n;
            do {
              --i;
              if ( P[i] > i ) swap( nr, Cj+i*nr, 1, Cj+P[i]*nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            nr, n, 1.0, T, n_x_2, Cj, nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Dj,  n,
             1.0, Cjp, nr
          );

          integer     jpp = std::min(j+k,eblock);
          valueType * Cpp = Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Ej,  n,
             1.0, Cpp, nr
          );
        }

        if ( nr_x_nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nx, n,
            -1.0, Cj,      nr,
                  Bj,      n,
             1.0, Fmat_th, nr // solo accumulo!
          );

        // NEXT STEP
        T   += k_x_2*Tsize;
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
    kBlock[nth] = k;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_reduced() {
    integer nblk = 2*usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      // -----------------------------------------------------------------------
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType * T   = Tmat + j*Tsize;
        integer   * P   = Perm + j*n;
        valueType * Djp = Dmat + jp*n_x_n;
        valueType * Dj  = Dmat + j*n_x_n;
        valueType * Ejp = Emat + jp*n_x_n;
        valueType * Ej  = Emat + j*n_x_n;

        buildT( 0, Ejp, Dj, T, P );

        alglin::zero( n_x_n, Dj, 1 );
        applyT( 0, T, P, Djp, n, Dj, n, n );

        alglin::zero( n_x_n, Ejp, 1 );
        applyT( 0, T, P, Ejp, n, Ej, n, n );

        if ( nx > 0 ) {
          valueType * Bj  = Bmat + j*n_x_nx;
          valueType * Bjp = Bmat + jp*n_x_nx;
          applyT( 0, T, P, Bjp, n, Bj, n, nx );
        }

        if ( nr > 0 ) {
          valueType * Cj  = Cmat + j*nr_x_n;
          valueType * Cjp = Cmat + jp*nr_x_n;
          if ( selected == BORDERED_QRP ) {
            integer i = n;
            do {
              --i;
              if ( P[i] > i ) swap( nr, Cj+i*nr, 1, Cj+P[i]*nr, 1 );
            } while ( i > 0 );
          }
          alglin::trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            nr, n, 1.0, T, n_x_2, Cj, nr
          );

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Dj,  n,
             1.0, Cjp, nr
          );

          integer     jpp = iBlock[std::min(jj+k,nblk)];
          valueType * Cpp = Cmat + jpp*nr_x_n;

          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Ej,  n,
             1.0, Cpp, nr
          );

        }

        if ( nr_x_nx > 0 ) {
          valueType * Cj = Cmat + j*nr_x_n;
          valueType * Bj = Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nx, n,
            -1.0, Cj,      nr,
                  Bj,      n,
             1.0, Fmat[0], nr
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
    valueType * Cnb = Cmat + nblock*nr_x_n;
    valueType * W0  = Hmat;
    valueType * WN  = W0+n*Nr;
    valueType * Wq  = WN+n*Nr;
    valueType * Wp  = Wq+qx*Nr;

    alglin::gecopy( n,  n,  Dmat, n, W0, Nr );
    alglin::gecopy( n,  n,  Emat, n, WN, Nr );
    alglin::gezero( n,  qx,          Wq, Nr );
    alglin::gecopy( n,  nx, Bmat, n, Wp, Nr );

    alglin::gecopy( n+qr, Nc, H0Nqp, n+qr, Hmat+n, Nr );

    integer offs = n_x_2+qr;

    alglin::gecopy( nr, n,  Cmat,    nr, W0+offs, Nr );
    alglin::gecopy( nr, n,  Cnb,     nr, WN+offs, Nr );
    alglin::gecopy( nr, qx, Cqmat,   nr, Wq+offs, Nr );
    alglin::gecopy( nr, nx, Fmat[0], nr, Wp+offs, Nr );

    bool ok = false;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      ok = last_lu.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_LUPQ:
      ok = last_lupq.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_QR:
      ok = last_qr.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_QRP:
      ok = last_qrp.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_SVD:
      ok = last_svd.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_LSS:
      ok = last_lss.factorize( Nr, Nc, Hmat, Nr );
      break;
    case BORDERED_LAST_LSY:
      ok = last_lsy.factorize( Nr, Nc, Hmat, Nr );
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
    valueType * X = x + (nblock-1)*n;
    swap( n, X, 1, x, 1 ); // uso x stesso come temporaneo
    bool ok = false;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:   ok = last_lu.solve( X );   break;
    case BORDERED_LAST_LUPQ: ok = last_lupq.solve( X ); break;
    case BORDERED_LAST_QR:   ok = last_qr.solve( X );   break;
    case BORDERED_LAST_QRP:  ok = last_qrp.solve( X );  break;
    case BORDERED_LAST_SVD:  ok = last_svd.solve( X );  break;
    case BORDERED_LAST_LSS:  ok = last_lss.solve( X );  break;
    case BORDERED_LAST_LSY:  ok = last_lsy.solve( X );  break;
    }
    if ( ok ) swap( n, X, 1, x, 1 );
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
    valueType * X = x + (nblock-1)*n;
    for ( integer i = 0; i < nrhs; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    bool ok = false;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:   ok = last_lu.solve( nrhs, X, ldX );   break;
    case BORDERED_LAST_LUPQ: ok = last_lupq.solve( nrhs, X, ldX ); break;
    case BORDERED_LAST_QR:   ok = last_qr.solve( nrhs, X, ldX );   break;
    case BORDERED_LAST_QRP:  ok = last_qrp.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_SVD:  ok = last_svd.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_LSS:  ok = last_lss.solve( nrhs, X, ldX );  break;
    case BORDERED_LAST_LSY:  ok = last_lsy.solve( nrhs, X, ldX );  break;
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
  BorderedCR<t_Value>::solve_CR( valueType x[] ) const {
    valueType * xb = x + (nblock+1)*n + qr; // deve essere b!
    if ( usedThread > 1 ) {
      if ( nr > 0 ) {
        for ( integer nt = 1; nt < usedThread; ++nt ) {
          alglin::zero( nr, xb_thread[size_t(nt)], 1 );
          pTP->run(
            nt, &BorderedCR<t_Value>::forward, this,
            nt, x, xb_thread[size_t(nt)]
          );
        }
        alglin::zero( nr, xb_thread[0], 1 );
        forward(0,x,xb);
        pTP->wait_all();
        for ( integer nt = 1; nt < usedThread; ++nt )
          axpy( nr, 1.0, xb_thread[size_t(nt)], 1, xb, 1 );
      } else {
        for ( integer nt = 1; nt < usedThread; ++nt )
          pTP->run(
            nt, &BorderedCR<t_Value>::forward, this,
            nt, x, xb_thread[size_t(nt)]
          );
        forward(0,x,xb);
        pTP->wait_all();
      }
      forward_reduced(x,xb);
    } else {
      forward(0,x,xb);
    }

    if ( !solve_last( x ) ) return false;

    if ( usedThread > 1 ) {
      backward_reduced(x);
      for ( integer nt = 1; nt < usedThread; ++nt )
        pTP->run( nt, &BorderedCR<t_Value>::backward, this, nt, x );
      backward(0,x);
      pTP->wait_all();
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
    if ( usedThread > 1 ) {
      for ( integer nt = 1; nt < usedThread; ++nt )
        pTP->run( nt, &BorderedCR<t_Value>::forward_n, this, nt, nrhs, rhs, ldRhs );
      forward_n( 0, nrhs, rhs, ldRhs );
      pTP->wait_all();
      forward_n_reduced( nrhs, rhs, ldRhs );
    } else {
      forward_n( 0, nrhs, rhs, ldRhs );
    }

    if ( !solve_last( nrhs, rhs, ldRhs ) ) return false;

    if ( usedThread > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < usedThread; ++nt )
        pTP->run( nt, &BorderedCR<t_Value>::backward_n, this, nt, nrhs, rhs, ldRhs );
      backward_n( 0, nrhs, rhs, ldRhs );
      pTP->wait_all();
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
    valueType * x0 = x    + iblock*n;
    valueType * T0 = Tmat + iblock*Tsize;
    integer   * P0 = Perm + iblock*n;
    valueType * C0 = Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      valueType * xj  = x0 + k*n;
      valueType * xjp = x0;
      valueType * T   = T0 + k*Tsize;
      integer   * P   = P0 + k*n;
      valueType * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;
      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, xj );
        if ( nr > 0 ) gemv( NO_TRANSPOSE, nr, n, -1.0, Cj, nr, xj, 1, 1.0, xb, 1 ); // solo accumulato
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*Tsize;
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
    valueType x[],
    integer   ldX
  ) const {
    valueType * xb = x + (nblock+1)*n + qr;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;
    valueType * x0 = x    + iblock*n;
    valueType * T0 = Tmat + iblock*Tsize;
    integer   * P0 = Perm + iblock*n;
    valueType * C0 = Cmat + iblock*nr_x_n;

    integer k = 1;
    while ( k < nblk ) {
      valueType * xj  = x0 + k*n;
      valueType * xjp = x0;
      valueType * T   = T0 + k*Tsize;
      integer   * P   = P0 + k*n;
      valueType * Cj  = C0 + k*nr_x_n;
      integer   k_x_2 = 2*k;

      for ( integer jj = k; jj < nblk; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( nr > 0 ) {
          spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nrhs, n,
            -1.0, Cj, nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          spin.unlock();
        }
        xj  += k_x_2*n;
        xjp += k_x_2*n;
        T   += k_x_2*Tsize;
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
    valueType x[],
    valueType xb[]
  ) const {
    integer nblk = 2*usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType const * T   = Tmat + j*Tsize;
        integer   const * P   = Perm + j*n;
        valueType       * xj  = x + j*n;
        valueType       * xjp = x + jp*n;
        applyT( 0, T, P, xjp, xj );
        if ( nr > 0 ) {
          valueType * Cj = Cmat + j*nr_x_n;
          alglin::gemv( NO_TRANSPOSE, nr, n, -1.0, Cj, nr, xj, 1, 1.0, xb, 1 );
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
    valueType * xb = x + (nblock+1)*n + qr;
    integer nblk = 2*usedThread-1;
    integer k = 1;
    while ( k < nblk ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer j  = iBlock[jj];
        integer jp = iBlock[jj-k];
        valueType const * T   = Tmat + j*Tsize;
        integer   const * P   = Perm + j*n;
        valueType       * xj  = x + j*n;
        valueType       * xjp = x + jp*n;
        applyT( 0, T, P, xjp, ldX, xj, ldX, nrhs );
        if ( nr > 0 ) {
          valueType * Cj = Cmat + j*nr_x_n;
          spin.lock();
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nrhs, n,
            -1.0, Cj, nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          spin.unlock();
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
    valueType * xn = x + (nblock+1)*n + qx;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    valueType * x0 = x    + iblock*n;
    valueType * B0 = Bmat + iblock*n_x_nx;
    valueType * D0 = Dmat + iblock*n_x_n;
    valueType * E0 = Emat + iblock*n_x_n;
    valueType * T0 = Tmat + iblock*Tsize;
    integer k = kBlock[nth];
    while ( (k/=2) > 0 ) {
      valueType * xj = x0 + k*n;
      valueType * xp = x0;
      valueType * Bj = B0 + k*n_x_nx;
      valueType * Dj = D0 + k*n_x_n;
      valueType * Ej = E0 + k*n_x_n;
      valueType * T  = T0 + k*Tsize;
      integer   k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        valueType * xpp = x + jpp*n;
        alglin::gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 );
        alglin::gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 );
        if ( nx > 0 ) gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 );
        alglin::trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 );
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n;
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
    valueType * xn = x + (nblock+1)*n + qx;
    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    valueType * x0 = x    + iblock*n;
    valueType * B0 = Bmat + iblock*n_x_nx;
    valueType * D0 = Dmat + iblock*n_x_n;
    valueType * E0 = Emat + iblock*n_x_n;
    valueType * T0 = Tmat + iblock*Tsize;
    integer k = kBlock[nth];
    while ( (k/=2) > 0 ) {
      valueType * xj = x0 + k*n;
      valueType * xp = x0;
      valueType * Bj = B0 + k*n_x_nx;
      valueType * Dj = D0 + k*n_x_n;
      valueType * Ej = E0 + k*n_x_n;
      valueType * T  = T0 + k*Tsize;
      integer  k_x_2 = 2*k;
      for ( integer j = iblock+k; j < eblock; j += k_x_2 ) {
        integer     jpp = std::min(j+k,eblock);
        valueType * xpp = x + jpp*n;
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
        if ( nx > 0 )
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        alglin::trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n;
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
        T  += k_x_2*Tsize;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_reduced( valueType x[] ) const {
    valueType * xn = x + (nblock+1)*n + qx;
    integer nblk = 2*usedThread-1;
    integer k = 1;
    while ( k < nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer     j   = iBlock[jj];
        integer     jp  = iBlock[jj-k];
        integer     jpp = iBlock[std::min(jj+k,nblk)];
        valueType * Dj  = Dmat + j*n_x_n;
        valueType * Ej  = Emat + j*n_x_n;
        valueType * xj  = x + j*n;
        valueType * xp  = x + jp*n;
        valueType * xpp = x + jpp*n;
        alglin::gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 );
        alglin::gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 );
        if ( nx > 0 ) {
          valueType * Bj = Bmat + j*n_x_nx;
          gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 );
        }
        valueType const * T = Tmat + j*Tsize;
        alglin::trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 );
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n;
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
    valueType x[],
    integer   ldX
  ) const {
    valueType * xn = x + (nblock+1)*n + qx;
    integer nblk = 2*usedThread-1;
    integer k = 1;
    while ( k < nblk ) k *= 2;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k; jj < nblk; jj += 2*k ) {
        integer     j   = iBlock[jj];
        integer     jp  = iBlock[jj-k];
        integer     jpp = iBlock[std::min(jj+k,nblk)];
        valueType * Dj  = Dmat + j*n_x_n;
        valueType * Ej  = Emat + j*n_x_n;
        valueType * xj  = x + j*n;
        valueType * xjp = x + jp*n;
        valueType * xpp = x + jpp*n;
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
        if ( nx > 0 ) {
          valueType * Bj = Bmat + j*n_x_nx;
          alglin::gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        valueType const * T = Tmat + j*Tsize;
        trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          n, nrhs, 1.0, T, n_x_2, xj, ldX
        );
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n;
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
  BorderedCR<t_Value>::Mv( valueType const x[], valueType res[] ) const {
    alglin::zero( numRows(), res, 1 );
    addMv( x, res );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addMv( valueType const x[], valueType res[] ) const {
    // internal blocks block
    t_Value const * D  = Dmat;
    t_Value const * E  = Emat;
    t_Value const * B  = Bmat;
    t_Value const * xx = x;
    t_Value const * xe = x  + nblock*n;
    t_Value const * xq = xe + n;
    t_Value const * xb = xq + qx;
    t_Value *       yy = res;
    for ( integer i = 0; i < nblock; ++i ) {
      alglin::gemv( NO_TRANSPOSE, n, n,  1.0, D, n, xx, 1, 1.0, yy, 1 );
      xx += n;
      alglin::gemv( NO_TRANSPOSE, n, n,  1.0, E, n, xx, 1, 1.0, yy, 1 );
      alglin::gemv( NO_TRANSPOSE, n, nx, 1.0, B, n, xb, 1, 1.0, yy, 1 );
      yy += n; D += n_x_n; E += n_x_n; B += n_x_nx;
    }

    integer     m = n+qr;
    valueType * H = H0Nqp;
    alglin::gemv( NO_TRANSPOSE, m, n,  1.0, H, m, x,  1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( NO_TRANSPOSE, m, n,  1.0, H, m, xe, 1, 1.0, yy, 1 ); H += m * n;
    alglin::gemv( NO_TRANSPOSE, m, qx, 1.0, H, m, xq, 1, 1.0, yy, 1 ); H += m * qx;
    alglin::gemv( NO_TRANSPOSE, m, nx, 1.0, H, m, xb, 1, 1.0, yy, 1 );

    if ( nr > 0 ) {
      yy += m;
      alglin::gemv( NO_TRANSPOSE, nr, nx, 1.0, Fmat[0], nr, xb, 1, 1.0, yy, 1 );
      t_Value const * C = Cmat;
      xx = x;
      for ( integer i = 0; i <= nblock; ++i ) {
        alglin::gemv( NO_TRANSPOSE, nr, n, 1.0, C, nr, xx, 1, 1.0, yy, 1 );
        xx += n; C += nr_x_n;
      }
      alglin::gemv( NO_TRANSPOSE, nr, qx, 1.0, Cqmat, nr, xx, 1, 1.0, yy, 1 );
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
    integer i0 = nbl*n + offs;
    integer j0 = (nblock+1)*n + qx + offs;
    for ( integer ij = 0; ij < n_x_nx; ++ij ) {
      I[ij] = i0 + (ij % n);
      J[ij] = j0 + integer(ij/n);
    }
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesB( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_nx, Bmat + nbl*n_x_nx, 1, V, 1 );
    return n_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternC(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
    integer i0 = (nblock+1)*n + qr + offs;
    integer j0 = nbl*n + offs;
    for ( integer ij = 0; ij < nr_x_n; ++ij ) {
      I[ij] = i0 + (ij % nr);
      J[ij] = j0 + integer(ij/nr);
    }
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesC( integer nbl, valueType V[] ) const {
    alglin::copy( nr_x_n, Cmat + nbl*nr_x_n, 1, V, 1 );
    return nr_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternD(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
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
  BorderedCR<t_Value>::valuesD( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_n, Dmat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternE(
    integer nbl, integer I[], integer J[], integer offs
  ) const {
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
  BorderedCR<t_Value>::valuesE( integer nbl, valueType V[] ) const {
    alglin::copy( n_x_n, Emat + nbl*n_x_n, 1, V, 1 );
    return n_x_n;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternF(
    integer I[], integer J[], integer offs
  ) const {
    integer i0 = (nblock+1)*n + qr + offs;
    integer j0 = (nblock+1)*n + qx + offs;
    for ( integer ij = 0; ij < nr_x_nx; ++ij ) {
      I[ij] = i0 + (ij % nr);
      J[ij] = j0 + integer(ij/nr);
    }
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesF( valueType V[] ) const {
    alglin::copy( nr_x_nx, Fmat[0], 1, V, 1 );
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternCq(
    integer I[], integer J[], integer offs
  ) const {
    integer i0 = (nblock+1)*n + qr + offs;
    integer j0 = (nblock+1)*n + offs;
    for ( integer ij = 0; ij < nr_x_qx; ++ij ) {
      I[ij] = i0 + (ij % nr);
      J[ij] = j0 + integer(ij/nr);
    }
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesCq( valueType V[] ) const {
    alglin::copy( nr_x_qx, Cqmat, 1, V, 1 );
    return nr_x_qx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternH( integer I[], integer J[], integer offs ) const {
    integer nqr = n + qr;
    integer nnz = nqr * Nc;
    integer i0 = nblock*n + offs;
    integer j0 = i0 - n;
    for ( integer ij = 0; ij < nnz; ++ij ) {
      I[ij] = i0 + (ij % nqr);
      integer j = integer(ij/nqr); if ( j >= n ) j += j0;
      J[ij] = j;
    }
    return nnz;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::valuesH( valueType V[] ) const {
    integer nnz = (n + qr) * Nc;
    alglin::copy( nnz, H0Nqp, 1, V, 1 );
    return nnz;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparsePattern(
    integer I[],
    integer J[],
    integer offs
  ) const {
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

    LW_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparsePattern( I, J, offs ), inserted {} values, expected {}\n",
      kkk, sparseNnz()
    );

  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparseValues( valueType V[] ) const {
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

    LW_ASSERT(
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
    integer const rH   = n*nblock;
    integer const rC   = rH + n+qr;
    integer const nrow = rC + nr;

    integer const cCq  = rH  + n;
    integer const cF   = cCq + qx;
    integer const ncol = cF  + nx;

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
      LW_ASSERT(
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

    int neq = int(this->numRows());
    int nnz = nblock * ( 2 * n_x_n + nx * n ) +
              ( n * (nblock+1) + (qx+nx) ) * nr +
              (n+qr) * (n_x_2+qx+nx);

    superluInteger.allocate( size_t(nnz+5*neq+1) );
    superluValue.allocate( size_t(nnz) );

    set_default_options(&slu_options);

    // Initialize the statistics variables.
    StatInit(&slu_stats);

    slu_perm_r = superluInteger( size_t(neq) ); /* row permutations from partial pivoting */
    slu_perm_c = superluInteger( size_t(neq) ); /* column permutation vector */
    slu_etree  = superluInteger( size_t(neq) );

    valueType * values = superluValue( size_t(nnz) );
    int       * rowind = superluInteger( size_t(nnz) );
    int       * colptr = superluInteger( size_t(neq+1) );

    // fill matrix
    int row1 = n*nblock;
    int row2 = row1+n+qr;
    int kk   = 0;
    int jj   = 0;
    colptr[jj] = 0;
    for ( int j = 0; j < n; ++j ) {
      for ( int i = 0; i < n; ++i ) {
        values[kk] = D( 0, i, j );
        rowind[kk] = i;
        ++kk;
      }
      for ( int i = 0; i < n+qr; ++i ) {
        values[kk] = H( i, j );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < nr; ++i ) {
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
        for ( int i = 0; i < nr; ++i ) {
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
      for ( int i = 0; i < n+qr; ++i ) {
        values[kk] = H( i, j+n );;
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < nr; ++i ) {
        values[kk] = C( nblock, i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    for ( int j = 0; j < qx; ++j ) {
      for ( int i = 0; i < n+qr; ++i ) {
        values[kk] = H( i, j+2*n );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < nr; ++i ) {
        values[kk] = Cq( i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    for ( int j = 0; j < nx; ++j ) {
      for ( int nbl = 0; nbl < nblock; ++nbl ) {
        for ( int i = 0; i < n; ++i ) {
          values[kk] = B( nbl, i, j );
          rowind[kk] = i+nbl*n;
          ++kk;
        }
      }
      for ( int i = 0; i < n+qr; ++i ) {
        values[kk] = H( i, j+2*n+qx );
        rowind[kk] = i+row1;
        ++kk;
      }
      for ( int i = 0; i < nr; ++i ) {
        values[kk] = F( i, j );
        rowind[kk] = i+row2;
        ++kk;
      }
      colptr[++jj] = kk;
    }

    LW_ASSERT(
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

    LW_ASSERT(
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
