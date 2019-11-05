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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

#ifdef ALGLIN_USE_SYSTEM_EIGEN
  #include <Eigen/Dense>
#else
  #include "Eigen/Dense"
#endif

namespace alglin {

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
    Nr      = n_x_2+nr+qr;
    Nc      = n_x_2+nx+qx;
    Tsize   = 2*n_x_n+n;

    integer N = std::max(Nr,Nc);
    Lwork = std::max(N,2*std::max(n_x_n,std::max(nr_x_n,n_x_nx)));

    valueType tmp; // get optimal allocation
    integer info = geqrf( Nr, Nc, nullptr, Nr, nullptr, &tmp, -1 );
    ALGLIN_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = " << info
    );
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp);

    info = geqp3( Nr, Nc, nullptr, Nr, nullptr, nullptr, &tmp, -1 );
    ALGLIN_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = " << info
    );
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp);

    LworkT = 2*n_x_n;
    info = geqrf( n_x_2, n, nullptr, n_x_2, nullptr, &tmp, -1 );
    ALGLIN_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqrf return info = " << info
    );
    LworkQR = integer(tmp);

    info = geqp3( n_x_2, n, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 );
    ALGLIN_ASSERT(
      info == 0,
      "BorderedCR::allocate call alglin::geqp3 return info = " << info
    );
    if ( LworkQR < integer(tmp) ) LworkQR = integer(tmp);

    usedThread = nblock >= 128*maxThread ? maxThread : nblock/128;
    if ( usedThread < 1 ) usedThread = 1;

    integer nnz = nblock*(n_x_nx+2*n_x_n+Tsize+n) +
                  (nblock+1)*nr_x_n +
                  nr*qx +
                  Nr*Nc + (n+qr)*Nc +
                  Lwork + (LworkT+LworkQR+nr_x_nx+nr)*usedThread;
    integer innz = nblock*n + (3+n)*usedThread;

    baseValue.allocate(size_t(nnz));
    baseInteger.allocate(size_t(innz));

    Bmat  = baseValue( size_t(nblock*n_x_nx) );
    Cmat  = baseValue( size_t((nblock+1)*nr_x_n) );
    Cqmat = baseValue( size_t(nr*qx) );
    Dmat  = baseValue( size_t(nblock*n_x_n) );
    Emat  = baseValue( size_t(nblock*n_x_n) );
    Fmat  = baseValue( size_t(nr_x_nx*usedThread) );

    Tmat  = baseValue( size_t(nblock*Tsize) );
    Ttau  = baseValue( size_t(nblock*n) );
    Hmat  = baseValue( size_t(Nr*Nc) );
    H0Nqp = baseValue( size_t((n+qr)*Nc) );

    Work   = baseValue( size_t(Lwork) );
    WorkT  = baseValue( size_t(LworkT*usedThread) );
    WorkQR = baseValue( size_t(LworkQR*usedThread) );

    xb_thread = baseValue( size_t(nr*usedThread) );

    Perm        = baseInteger( size_t(nblock*n) );
    perm_thread = baseInteger( size_t(n*usedThread) );
    iBlock      = baseInteger( size_t(2*usedThread) );
    kBlock      = baseInteger( size_t(usedThread) );

    // precompute partition for parallel computation
    if ( usedThread > 1 ) {
      iBlock[0] = 0;
      iBlock[1] = static_cast<integer>(nblock/usedThread);
      for ( integer nt = 1; nt < usedThread; ++nt ) {
        iBlock[2*nt+0] = iBlock[2*nt-1]+1;
        iBlock[2*nt+1] = static_cast<integer>(((nt+1)*nblock)/usedThread);
      }
    } else {
      iBlock[0] = 0;
      iBlock[1] = nblock;
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
    allocate( M.nblock, M.n, M.qr, M.qx, M.nr, M.nx );
    copy( nblock*n_x_nx,     M.Bmat,  1, Bmat,  1 );
    copy( (nblock+1)*nr_x_n, M.Cmat,  1, Cmat,  1 );
    copy( nblock*n_x_n,      M.Dmat,  1, Dmat,  1 );
    copy( nblock*n_x_n,      M.Emat,  1, Emat,  1 );
    copy( nr_x_nx,           M.Fmat,  1, Fmat,  1 );
    copy( (n+qr)*Nc,         M.H0Nqp, 1, H0Nqp, 1 );
    copy( nr*qx,             M.Cqmat, 1, Cqmat, 1 );
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
    gecopy( m, n,  H0, ld0, H, m ); H += m * n;
    gecopy( m, n,  HN, ldN, H, m ); H += m * n;
    gecopy( m, qx, Hq, ldQ, H, m ); H += m * qx;
    gecopy( m, nx, Hp, ldP, H, m );
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

    ALGLIN_ASSERT(
      H0.numRows() == m && H0.numCols() == n,
      "loadBottom, bad dimension size(H0) = " << H0.numRows() <<
      " x " << H0.numCols() << " expected " << m << " x " << n
    );

    ALGLIN_ASSERT(
      HN.numRows() == m && HN.numCols() == n,
      "loadBottom, bad dimension size(HN) = " << HN.numRows() <<
      " x " << HN.numCols() << " expected " << m << " x " << n
    );

    ALGLIN_ASSERT(
      Hq.numRows() == m && Hq.numCols() == qx,
      "loadBottom, bad dimension size(Hq) = " << Hq.numRows() <<
      " x " << Hq.numCols() << " expected " << m << " x " << qx
    );

    ALGLIN_ASSERT(
      Hp.numRows() == m && Hp.numCols() == nx,
      "loadBottom, bad dimension size(Hp) = " << Hp.numRows() <<
      " x " << Hp.numCols() << " expected " << m << " x " << nx
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
  BorderedCR<t_Value>::loadBottom( MatrixWrapper<valueType> const & H ) {
    integer m = n + qr;
    ALGLIN_ASSERT(
      H.numRows() == m && H.numCols() == Nc,
      "loadBottom, bad dimension size(H) = " << H.numRows() <<
      " x " << H.numCols() << " expected " << m << " x " << Nc
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
    gecopy( nr, nx, F,  ldF,  this->Fmat,              nr );
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadBottom2(
    MatrixWrapper<valueType> const & C0,
    MatrixWrapper<valueType> const & CN,
    MatrixWrapper<valueType> const & Cq,
    MatrixWrapper<valueType> const & F
  ) {

    ALGLIN_ASSERT(
      C0.numRows() == nr && C0.numCols() == n,
      "loadBottom2, bad dimension size(C0) = " << C0.numRows() <<
      " x " << C0.numCols() << " expected " << nr << " x " << n
    );

    ALGLIN_ASSERT(
      CN.numRows() == nr && CN.numCols() == n,
      "loadBottom2, bad dimension size(CN) = " << CN.numRows() <<
      " x " << CN.numCols() << " expected " << nr << " x " << n
    );

    ALGLIN_ASSERT(
      Cq.numRows() == nr && Cq.numCols() == qx,
      "loadBottom2, bad dimension size(Cq) = " << Cq.numRows() <<
      " x " << Cq.numCols() << " expected " << nr << " x " << qx
    );

    ALGLIN_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "loadBottom2, bad dimension size(F) = " << F.numRows() <<
      " x " << F.numCols() << " expected " << nr << " x " << nx
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
  BorderedCR<t_Value>::loadB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    ALGLIN_ASSERT(
      B.numRows() == n && B.numCols() == nx,
      "loadB( " << nbl << ", B) bad dimension size(B) = " << B.numRows() <<
      " x " << B.numCols() << " expected " << n << " x " << nx
    );
    gecopy( n, nx, B.get_data(), B.lDim(), Bmat + nbl*n_x_nx, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoB(
    integer                          nbl,
    MatrixWrapper<valueType> const & B
  ) {
    ALGLIN_ASSERT(
      B.numRows() == n && B.numCols() == nx,
      "addtoB( " << nbl << ", B) bad dimension size(B) = " << B.numRows() <<
      " x " << B.numCols() << " expected " << n << " x " << nx
    );
    valueType * BB = Bmat + nbl*n_x_nx;
    geadd( n, nx, 1.0, B.get_data(), B.lDim(), 1.0, BB, n, BB, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    ALGLIN_ASSERT(
      C.numRows() == nr && C.numCols() == n,
      "loadC( " << nbl << ", C) bad dimension size(C) = " << C.numRows() <<
      " x " << C.numCols() << " expected " << nr << " x " << n
    );
    gecopy( nr, n, C.get_data(), C.lDim(), Cmat + nbl*nr_x_n, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    ALGLIN_ASSERT(
      C.numRows() == nr && C.numCols() == n,
      "addtoC( " << nbl << ", C) bad dimension size(C) = " << C.numRows() <<
      " x " << C.numCols() << " expected " << nr << " x " << n
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    geadd( nr, n, 1.0, C.get_data(), C.lDim(), 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoC2(
    integer                          nbl,
    MatrixWrapper<valueType> const & C
  ) {
    ALGLIN_ASSERT(
      C.numRows() == nr && C.numCols() == n_x_2,
      "addtoC( " << nbl << ", C) bad dimension size(C) = " << C.numRows() <<
      " x " << C.numCols() << " expected " << nr << " x " << n_x_2
    );
    valueType * CC = Cmat + nbl*nr_x_n;
    geadd( nr, n_x_2, 1.0, C.get_data(), C.lDim(), 1.0, CC, nr, CC, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadD(
    integer                          nbl,
    MatrixWrapper<valueType> const & D
  ) {
    ALGLIN_ASSERT(
      D.numRows() == n && D.numCols() == n,
      "loadD( " << nbl << ", D) bad dimension size(D) = " << D.numRows() <<
      " x " << D.numCols() << " expected " << n << " x " << n
    );
    gecopy( n, n, D.get_data(), D.lDim(), Dmat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadE(
    integer                          nbl,
    MatrixWrapper<valueType> const & E
  ) {
    ALGLIN_ASSERT(
      E.numRows() == n && E.numCols() == n,
      "loadE( " << nbl << ", E) bad dimension size(E) = " << E.numRows() <<
      " x " << E.numCols() << " expected " << n << " x " << n
    );
    gecopy( n, n, E.get_data(), E.lDim(), Emat + nbl*n_x_n, n );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadF( MatrixWrapper<valueType> const & F ) {
    ALGLIN_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "loadF(F) bad dimension size(F) = " << F.numRows() <<
      " x " << F.numCols() << " expected " << nr << " x " << nx
    );
    gecopy( nr, nx, F.get_data(), F.lDim(), Fmat, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::addtoF( MatrixWrapper<valueType> const & F ) {
    ALGLIN_ASSERT(
      F.numRows() == nr && F.numCols() == nx,
      "addtoF(F) bad dimension size(F) = " << F.numRows() <<
      " x " << F.numCols() << " expected " << nr << " x " << nx
    );
    geadd( nr, nx, 1.0, F.get_data(), F.lDim(), 1.0, Fmat, nr, Fmat, nr );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  BorderedCR<t_Value>::loadCq( MatrixWrapper<valueType> const & Cq ) {
    ALGLIN_ASSERT(
      Cq.numRows() == nr && Cq.numCols() == qx,
      "loadCq(Cq) bad dimension size(Cq) = " << Cq.numRows() <<
      " x " << Cq.numCols() << " expected " << nr << " x " << qx
    );
    gecopy( nr, qx, Cq.get_data(), Cq.lDim(), Cqmat, nr );
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
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      // fill zero F(...)
      if ( nr_x_nx > 0 ) alglin::zero( nr_x_nx*(usedThread-1), Fmat + nr_x_nx, 1 );
      // launch thread
      for ( integer nt = 1; nt < usedThread; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::factorize_block, this, nt );
      factorize_block(0);
      // wait thread
      for ( integer nt = 1; nt < usedThread; ++nt ) threads[nt].join();
      // accumulate F(...)
      if ( nr_x_nx > 0 )
        for ( integer nt = 1; nt < usedThread; ++nt )
          alglin::axpy( nr_x_nx, 1.0, Fmat + nt*nr_x_nx, 1, Fmat, 1 );
      factorize_reduced();
    } else {
      factorize_block(0);
    }
    #else
    factorize_block(0);
    #endif
    load_and_factorize_last();
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
    gecopy( n, n, TOP,    n, T,   n_x_2 );
    gecopy( n, n, BOTTOM, n, T+n, n_x_2 );
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      info = getrf( n_x_2, n, T, n_x_2, iperm );
      break;
    case BORDERED_QR:
      info = geqrf( n_x_2, n, T, n_x_2, T+2*n_x_n, WorkQR+nth*LworkQR, LworkQR );
      break;
    case BORDERED_QRP:
      { integer * P = perm_thread+nth*n;
        std::fill( P, P+n, 0 );
        info = geqp3( n_x_2, n, T, n_x_2, P, T+2*n_x_n, WorkQR+nth*LworkQR, LworkQR );
        if ( info == 0 ) permutation_to_exchange( n, P, iperm );
      }
      break;
    case BORDERED_SUPERLU:
      ALGLIN_ERROR( "BorderedCR::buildT, cannot be used with SUPERLU" );
      break;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize INFO = " << info );
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
    valueType * W = WorkT + LworkT*nth;
    gecopy( n, ncol, TOP,    ldTOP,    W,   n_x_2 );
    gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, n_x_2 );
    // Apply row interchanges to the right hand sides.
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      info = swaps( ncol, W, n_x_2, 0, n-1, iperm, 1 );
      ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info );
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, n_x_2, W, n_x_2 );
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        n, ncol, n,
        -1.0, T+n, n_x_2,  // M
              W,   n_x_2,  // L^(-1) TOP
         1.0, W+n, n_x_2   // TOP = BOTTOM - M L^(-1) TOP
      );
      break;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = ormqr(
        LEFT, TRANSPOSE,
        n_x_2, ncol, // righe x colonne
        n,           // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        WorkQR+nth*LworkQR, LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      ALGLIN_ERROR( "BorderedCR::applyT, cannot be used with SUPERLU" );
      break;
    }
    gecopy( n, ncol, W+n, n_x_2, TOP,    ldTOP );
    gecopy( n, ncol, W,   n_x_2, BOTTOM, ldBOTTOM );
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
    valueType * W = WorkT + LworkT*nth;
    memcpy( W,   TOP,    n*sizeof(valueType) );
    memcpy( W+n, BOTTOM, n*sizeof(valueType) );
    //copy( n, TOP,    1, W,   1 );
    //copy( n, BOTTOM, 1, W+n, 1 );
    integer info = 0;
    switch ( selected ) {
    case BORDERED_LU:
      // Apply row interchanges to the right hand sides.
      info = swaps( 1, W, n_x_2, 0, n-1, iperm, 1 );
      if ( info == 0 ) {
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, n_x_2, W, 1 );
        gemv(
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
      info = ormqr(
        LEFT, TRANSPOSE,
        n_x_2, 1, // righe x colonne
        n,        // numero riflettori usati nel prodotto Q
        T, n_x_2,
        T+2*n_x_n,
        W, n_x_2,
        WorkQR+nth*LworkQR, LworkQR
      );
      break;
    case BORDERED_SUPERLU:
      ALGLIN_ERROR( "BorderedCR::applyT, cannot be used with SUPERLU" );
      break;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info );
    memcpy( TOP,    W+n, n*sizeof(valueType) );
    memcpy( BOTTOM, W,   n*sizeof(valueType) );
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

  //#define TRY_EIGEN

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_block( integer nth ) {

    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,Eigen::Dynamic> dmat_t;
    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,1>              dvec_t;

    integer iblock = iBlock[2*nth+0];
    integer eblock = iBlock[2*nth+1];
    integer nblk   = eblock - iblock;

    valueType * Bmat0 = Bmat + iblock*n_x_nx;
    valueType * Cmat0 = Cmat + iblock*nr_x_n;
    valueType * Dmat0 = Dmat + iblock*n_x_n;
    valueType * Emat0 = Emat + iblock*n_x_n;
    valueType * T0    = Tmat + iblock*Tsize;
    integer   * P0    = Perm + iblock*n;

    valueType * Fmat_th = Fmat + nth * nr_x_nx;

    #ifdef TRY_EIGEN
    Eigen::Map<dmat_t> F_mat(Fmat_th,nr,nx);
    #endif

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

        #ifdef TRY_EIGEN
        Eigen::Map<dmat_t> Bj_mat(Bj,n,nx);
        Eigen::Map<dmat_t> Cj_mat(Cj,nr,n);
        Eigen::Map<dmat_t> Cjp_mat(Cjp,nr,n);
        Eigen::Map<dmat_t> Dj_mat(Dj,n,n);
        Eigen::Map<dmat_t> Djp_mat(Djp,n,n);
        Eigen::Map<dmat_t> Ej_mat(Ej,n,n);
        Eigen::Map<dmat_t> Ejp_mat(Ejp,n,n);
        #endif

        buildT( nth, Ejp, Dj, T, P );

        #ifdef TRY_EIGEN
        Dj_mat.setZero();
        #else
        alglin::zero( n_x_n, Dj, 1 );
        #endif
        applyT( nth, T, P, Djp, n, Dj, n, n );

        #ifdef TRY_EIGEN
        Ejp_mat.setZero();
        #else
        alglin::zero( n_x_n, Ejp, 1 );
        #endif
        applyT( nth, T, P, Ejp, n, Ej, n, n );

        if ( nx > 0 ) applyT( nth, T, P, Bjp, n, Bj, n, nx );

        if ( nr > 0 ) {

          if ( selected == BORDERED_QRP ) {
            integer i = n;
            do {
              --i;
              if ( P[i] > i ) {
                #ifdef TRY_EIGEN
                Cj_mat.row(i).swap(Cj_mat.row(P[i]));
                #else
                swap( nr, Cj+i*nr, 1, Cj+P[i]*nr, 1 );
                #endif
              }
            } while ( i > 0 );
          }
          #ifdef TRY_EIGEN
          Eigen::Map<dmat_t,Eigen::Unaligned,Eigen::OuterStride<> >
            T_mat(T,n,n, Eigen::OuterStride<>(n_x_2));
          T_mat . template triangularView<Eigen::Upper>()
                . template solveInPlace<Eigen::OnTheRight>(Cj_mat);
          #else
          trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            nr, n, 1.0, T, n_x_2, Cj, nr
          );
          #endif

          #ifdef TRY_EIGEN
          Cjp_mat.noalias() -= Cj_mat * Dj_mat;
          #else
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Dj,  n,
             1.0, Cjp, nr
          );
          #endif

          integer     jpp = std::min(j+k,eblock);
          valueType * Cpp = Cmat + jpp*nr_x_n;

          #ifdef TRY_EIGEN
          Eigen::Map<dmat_t> Cpp_mat(Cpp,nr,n);
          Cpp_mat.noalias() -= Cj_mat * Ej_mat;
          #else
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Ej,  n,
             1.0, Cpp, nr
          );
          #endif
        }

        if ( nr_x_nx > 0 ) {
          #ifdef TRY_EIGEN
          F_mat.noalias() -= Cj_mat * Bj_mat;
          #else
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nx, n,
            -1.0, Cj,      nr,
                  Bj,      n,
             1.0, Fmat_th, nr // solo accumulo!
          );
          #endif
        }

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
          trsm(
            RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
            nr, n, 1.0, T, n_x_2, Cj, nr
          );

          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, n, n,
            -1.0, Cj,  nr,
                  Dj,  n,
             1.0, Cjp, nr
          );

          integer     jpp = iBlock[std::min(jj+k,nblk)];
          valueType * Cpp = Cmat + jpp*nr_x_n;

          gemm(
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
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nx, n,
            -1.0, Cj,   nr,
                  Bj,   n,
             1.0, Fmat, nr
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
  void
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

    gecopy( n,  n,  Dmat, n, W0, Nr );
    gecopy( n,  n,  Emat, n, WN, Nr );
    gezero( n,  qx,          Wq, Nr );
    gecopy( n,  nx, Bmat, n, Wp, Nr );

    gecopy( n+qr, Nc, H0Nqp, n+qr, Hmat+n, Nr );

    integer offs = n_x_2+qr;

    gecopy( nr, n,  Cmat,  nr, W0+offs, Nr );
    gecopy( nr, n,  Cnb,   nr, WN+offs, Nr );
    gecopy( nr, qx, Cqmat, nr, Wq+offs, Nr );
    gecopy( nr, nx, Fmat,  nr, Wp+offs, Nr );

    integer info = 0;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      last_lu.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_LUPQ:
      last_lupq.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_QR:
      last_qr.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_QRP:
      last_qrp.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_SVD:
      last_svd.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_LSS:
      last_lss.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    case BORDERED_LAST_LSY:
      last_lsy.factorize(
        "BorderedCR::load_and_factorize_last",
        Nr, Nc, Hmat, Nr
      );
      break;
    }
    ALGLIN_ASSERT(
      info == 0,
      "BorderedCR::factorize_last INFO = " << info <<
      " Nr = " << Nr << " Nc = " << Nc
    );
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
  void
  BorderedCR<t_Value>::solve_last( valueType x[] ) const {
    valueType * X = x + (nblock-1)*n;
    swap( n, X, 1, x, 1 ); // uso x stesso come temporaneo
    integer info = 0;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      last_lu.solve( X );
      break;
    case BORDERED_LAST_LUPQ:
      last_lupq.solve( X );
      break;
    case BORDERED_LAST_QR:
      last_qr.solve( X );
      break;
    case BORDERED_LAST_QRP:
      last_qrp.solve( X );
      break;
    case BORDERED_LAST_SVD:
      last_svd.solve( X );
      break;
    case BORDERED_LAST_LSS:
      last_lss.solve( X );
      break;
    case BORDERED_LAST_LSY:
      last_lsy.solve( X );
      break;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info );
    swap( n, X, 1, x, 1 );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last(
    integer   nrhs,
    valueType x[],
    integer   ldX
  ) const {
    valueType * X = x + (nblock-1)*n;
    for ( integer i = 0; i < nrhs; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 );
    integer info = 0;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      last_lu.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_LUPQ:
      last_lupq.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_QR:
      last_qr.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_QRP:
      last_qrp.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_SVD:
      last_svd.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_LSS:
      last_lss.solve( nrhs, X, ldX );
      break;
    case BORDERED_LAST_LSY:
      last_lsy.solve( nrhs, X, ldX );
      break;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info );
    for ( integer i = 0; i < nrhs; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 );
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_CR( valueType x[] ) const {
    valueType * xb = x + (nblock+1)*n + qr; // deve essere b!
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      if ( nr > 0 ) {
        alglin::zero( nr*usedThread, xb_thread, 1 );
        for ( integer nt = 1; nt < usedThread; ++nt )
          threads[nt] = std::thread( &BorderedCR<t_Value>::forward, this,
                                     nt, x, xb_thread+nr*nt );
        forward(0,x,xb);
        for ( integer nt = 1; nt < usedThread; ++nt ) {
          threads[nt].join();
          axpy( nr, 1.0, xb_thread+nr*nt, 1, xb, 1 );
        }
      } else {
        for ( integer nt = 1; nt < usedThread; ++nt )
          threads[nt] = std::thread( &BorderedCR<t_Value>::forward, this,
                                     nt, x, xb_thread+nr*nt );
        forward(0,x,xb);
        for ( integer nt = 1; nt < usedThread; ++nt ) threads[nt].join();
      }
      forward_reduced(x,xb);
    } else {
      forward(0,x,xb);
    }
    #else
      forward(0,x,xb);
    #endif

    solve_last( x );

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      backward_reduced(x);
      for ( integer nt = 1; nt < usedThread; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::backward, this, nt, x );
      backward(0,x);
      for ( integer nt = 1; nt < usedThread; ++nt ) threads[nt].join();
    } else {
      backward(0,x);
    }
    #else
    backward(0,x);
    #endif
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_CR(
    integer   nrhs,
    valueType rhs[],
    integer   ldRhs
  ) const {
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      for ( integer nt = 1; nt < usedThread; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::forward_n, this,
                                   nt, nrhs, rhs, ldRhs );
      forward_n( 0, nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < usedThread; ++nt ) threads[nt].join();
      forward_n_reduced( nrhs, rhs, ldRhs );
    } else {
      forward_n( 0, nrhs, rhs, ldRhs );
    }
    #else
      forward_n( 0, nrhs, rhs, ldRhs );
    #endif

    solve_last( nrhs, rhs, ldRhs );

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < usedThread; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::backward_n, this,
                                   nt, nrhs, rhs, ldRhs );
      backward_n( 0, nrhs, rhs, ldRhs );
      for ( integer nt = 1; nt < usedThread; ++nt ) threads[nt].join();
    } else {
      backward_n( 0, nrhs, rhs, ldRhs );
    }
    #else
      backward_n( 0, nrhs, rhs, ldRhs );
    #endif
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
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.lock();
          #endif
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nrhs, n,
            -1.0, Cj, nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.unlock();
          #endif
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
          gemv( NO_TRANSPOSE, nr, n, -1.0, Cj, nr, xj, 1, 1.0, xb, 1 );
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
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.lock();
          #endif
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            nr, nrhs, n,
            -1.0, Cj, nr,
                  xj, ldX,
             1.0, xb, ldX
          );
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.unlock();
          #endif
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
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 );
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 );
        if ( nx > 0 ) gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 );
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 );
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
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Dj,  n,
                xp,  ldX,
           1.0, xj,  ldX
        );
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( nx > 0 )
          gemm(
            NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, nx,
            -1.0, Bj, n,
                  xn, ldX,
             1.0, xj, ldX
          );
        trsm(
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
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 );
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 );
        if ( nx > 0 ) {
          valueType * Bj = Bmat + j*n_x_nx;
          gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 );
        }
        valueType const * T = Tmat + j*Tsize;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 );
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
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Dj,  n,
                xjp, ldX,
           1.0, xj,  ldX
        );
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, Ej,  n,
                xpp, ldX,
           1.0, xj,  ldX
        );
        if ( nx > 0 ) {
          valueType * Bj = Bmat + j*n_x_nx;
          gemm(
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
              swap( nrhs, xj+i, ldX, xj+P[i], ldX );
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
      gemv( NO_TRANSPOSE, n, n,  1.0, D, n, xx, 1, 1.0, yy, 1 );
      xx += n;
      gemv( NO_TRANSPOSE, n, n,  1.0, E, n, xx, 1, 1.0, yy, 1 );
      gemv( NO_TRANSPOSE, n, nx, 1.0, B, n, xb, 1, 1.0, yy, 1 );
      yy += n; D += n_x_n; E += n_x_n; B += n_x_nx;
    }

    integer     m = n+qr;
    valueType * H = H0Nqp;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, x,  1, 1.0, yy, 1 ); H += m * n;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, xe, 1, 1.0, yy, 1 ); H += m * n;
    gemv( NO_TRANSPOSE, m, qx, 1.0, H, m, xq, 1, 1.0, yy, 1 ); H += m * qx;
    gemv( NO_TRANSPOSE, m, nx, 1.0, H, m, xb, 1, 1.0, yy, 1 );

    if ( nr > 0 ) {
      yy += m;
      gemv( NO_TRANSPOSE, nr, nx, 1.0, Fmat, nr, xb, 1, 1.0, yy, 1 );
      t_Value const * C = Cmat;
      xx = x;
      for ( integer i = 0; i <= nblock; ++i ) {
        gemv( NO_TRANSPOSE, nr, n, 1.0, C, nr, xx, 1, 1.0, yy, 1 );
        xx += n; C += nr_x_n;
      }
      gemv( NO_TRANSPOSE, nr, qx, 1.0, Cqmat, nr, xx, 1, 1.0, yy, 1 );
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
    alglin::copy( nr_x_nx, Fmat, 1, V, 1 );
    return nr_x_nx;
  }

  template <typename t_Value>
  integer
  BorderedCR<t_Value>::patternCq(
    integer I[], integer J[], integer offs
  ) const {
    integer nr_x_qx = nr*qx;
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
    integer nr_x_qx = nr*qx;
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

#if 1
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

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparsePattern( I, J, offs ), inserted " << kkk <<
      " values, expected " << sparseNnz()
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

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparseValues( V ), inserted " << kkk <<
      " values, expected " << sparseNnz()
    );

  }

#else

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparsePattern(
    integer I[],
    integer J[],
    integer offs
  ) const {
    integer kkk = 0;

    // D, E, B
    integer je = nblock*n;
    integer jq = je+n;
    integer jb = jq+qx;
    for ( integer k = 0; k < nblock; ++k ) {
      integer ii = k*n;
      for ( integer i = 0; i < n; ++i ) {
        integer iii = ii+i;
        for ( integer j = 0; j < n; ++j ) {
          I[kkk] = iii; J[kkk] = ii+j;   ++kkk; // D[i+j*n]
          I[kkk] = iii; J[kkk] = ii+n+j; ++kkk; // E[i+j*n]
        }
        for ( integer j = 0; j < nx; ++j ) {
          I[kkk] = iii; J[kkk] = jb+j; ++kkk; // B[i+j*n]
        }
      }
    }

    // H0
    integer ie = nblock*n;
    for ( integer i = 0; i < n+qr; ++i ) {
      for ( integer j = 0; j < n; ++j ) {
        I[kkk] = ie+i; J[kkk] = j; ++kkk; // H[i+j*(n+qr)]
      }
    }

    // HNp
    for ( integer i = 0; i < n+qr; ++i ) {
      for ( integer j = 0; j < n+qx+nx; ++j ) {
        I[kkk] = ie+i; J[kkk] = je+j; ++kkk; // H[i+j*(n+qr)]
      }
    }

    // C
    ie += n+qr;
    for ( integer k = 0; k <= nblock; ++k ) {
      integer ii = k*n;
      for ( integer i = 0; i < nr; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          I[kkk] = ie+i; J[kkk] = ii+j; ++kkk; // C[i+j*nr]
        }
      }
    }

    // F
    for ( integer i = 0; i < nr; ++i ) {
      for ( integer j = 0; j < nx; ++j ) {
        I[kkk] = ie+i; J[kkk] = jb+j; ++kkk; // Fmat[i+j*nr]
      }
    }

    // Cq
    for ( integer i = 0; i < nr; ++i ) {
      for ( integer j = 0; j < qx; ++j ) {
        I[kkk] = ie+i; J[kkk] = jq+j; ++kkk; // Cqmat[i+j*nr]
      }
    }

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparsePattern( V ), inserted " << kkk <<
      " values, expected " << sparseNnz()
    );

    if ( offs != 0 ) {
      for ( integer kkk = 0; kkk < sparseNnz(); ++kkk ) {
        I[kkk] += offs;
        J[kkk] += offs;
      }
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::sparseValues( valueType V[] ) const {
    // D, E, B
    integer kkk = 0;
    for ( integer k = 0; k < nblock; ++k ) {
      //integer ii = offs+k*n;
      valueType * D = Dmat + k * n_x_n;
      valueType * E = Emat + k * n_x_n;
      valueType * B = Bmat + k * n_x_nx;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          V[kkk++] = D[i+j*n];
          V[kkk++] = E[i+j*n];
        }
        for ( integer j = 0; j < nx; ++j ) V[kkk++] = B[i+j*n];
      }
    }

    // H0
    valueType * H = this->H0Nqp;
    for ( integer i = 0; i < n+qr; ++i )
      for ( integer j = 0; j < n; ++j )
        V[kkk++] = H[i+j*(n+qr)];
    H += n*(n+qr);

    // HNp
    for ( integer i = 0; i < n+qr; ++i )
      for ( integer j = 0; j < n+qx+nx; ++j )
        V[kkk++] = H[i+j*(n+qr)];

    // C
    for ( integer k = 0; k <= nblock; ++k ) {
      valueType * C = Cmat + k * nr_x_n;
      for ( integer i = 0; i < nr; ++i )
        for ( integer j = 0; j < n; ++j )
          V[kkk++] = C[i+j*nr];
    }

    // F
    for ( integer i = 0; i < nr; ++i )
      for ( integer j = 0; j < nx; ++j )
        V[kkk++] = Fmat[i+j*nr];

    // Cq
    for ( integer i = 0; i < nr; ++i )
      for ( integer j = 0; j < qx; ++j )
        V[kkk++] = Cqmat[i+j*nr];

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BorderedCR::sparseValues( V ), inserted " << kkk <<
      " values, expected " << sparseNnz()
    );
  }
#endif

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
      ALGLIN_ASSERT(
        ok,
        "in BorderedCR<t_Value>::sparseLoad, "
        "indices (i,j) = ( " << M_row[kkk] <<
        ", " << M_col[kkk] << ") out on pattern!"
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

    ALGLIN_ASSERT(
      kk == nnz,
      "BABD_SuperLU::factorize -- dgstrf() error nnz = " << nnz <<
      " != " << kk
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

    ALGLIN_ASSERT( info == 0, "BABD_SuperLU::factorize -- [sd]gstrf() error returns INFO = " << info );
    //cout << "done\n";

  }

  template <typename t_Value>
  void
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

    ALGLIN_ASSERT( info == 0, "BABD_SuperLU::solve_SuperLU -- gstrs() error returns INFO= " << info );
  }

  template <typename t_Value>
  void
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

    ALGLIN_ASSERT( info == 0, "BABD_SuperLU::solve_SuperLU -- gstrs() error returns INFO= " << info );
  }

  // ---------------------------------------------------------------------------

  template class BorderedCR<float>;
  template class BorderedCR<double>;

}
