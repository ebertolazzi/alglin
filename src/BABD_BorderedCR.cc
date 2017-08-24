/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wweak-template-vtables"
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
  BorderedCR<t_Value>::allocate( integer _nblock,
                                 integer _n,
                                 integer _qr,
                                 integer _qx,
                                 integer _nr,
                                 integer _nx ) {

    if ( nblock == _nblock && n  == _n  &&
         qr     == _qr     && qx == _qx &&
         nr     == _nr     && nx == _nx ) return ;

    nblock  = _nblock ;
    n       = _n ;
    qr      = _qr ;
    qx      = _qx ;
    nr      = _nr ;
    nx      = _nx ;
    n_x_2   = n*2 ;
    n_x_n   = n*n ;
    nr_x_n  = n*nr ;
    n_x_nx  = n*nx ;
    nr_x_nx = nr*nx ;
    Nr      = n_x_2+nr+qr ;
    Nc      = n_x_2+nx+qx ;
    Tsize   = 2*n_x_n+n ;

    integer N = std::max(Nr,Nc) ;
    Lwork = std::max(N,2*std::max(n_x_n,std::max(nr_x_n,n_x_nx))) ;

    valueType tmp ; // get optimal allocation
    integer info = geqrf( Nr, Nc, nullptr, Nr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    info = geqp3( Nr, Nc, nullptr, Nr, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    LworkT = 2*n_x_n ;
    info = geqrf( n_x_2, n, nullptr, n_x_2, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    LworkQR = integer(tmp) ;

    info = geqp3( n_x_2, n, nullptr, n_x_2, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( LworkQR < integer(tmp) ) LworkQR = integer(tmp) ;

    usedThread = nblock >= 128*maxThread ? maxThread : nblock/128 ;
    if ( usedThread < 1 ) usedThread = 1 ;

    integer nnz = nblock*(n_x_nx+2*n_x_n+Tsize+n) +
                  (nblock+1)*nr_x_n +
                  nr*qx +
                  Nr*Nc + N + (n+qr)*Nc +
                  Lwork + (LworkT+LworkQR+nr_x_nx+nr)*usedThread ;
    integer innz = 2*N + nblock*n + (3+n)*usedThread ;

    baseValue.allocate(size_t(nnz)) ;
    baseInteger.allocate(size_t(innz)) ;

    Bmat  = baseValue(size_t(nblock*n_x_nx)) ;
    Cmat  = baseValue(size_t((nblock+1)*nr_x_n)) ;
    Cqmat = baseValue(size_t(nr*qx)) ;
    Dmat  = baseValue(size_t(nblock*n_x_n)) ;
    Emat  = baseValue(size_t(nblock*n_x_n)) ;
    Fmat  = baseValue(size_t(nr_x_nx*usedThread)) ;

    Tmat  = baseValue(size_t(nblock*Tsize)) ;
    Ttau  = baseValue(size_t(nblock*n)) ;
    Hmat  = baseValue(size_t(Nr*Nc)) ;
    Htau  = baseValue(size_t(N)) ;
    H0Nqp = baseValue(size_t((n+qr)*Nc)) ;

    Work   = baseValue(size_t(Lwork)) ;
    WorkT  = baseValue(size_t(LworkT*usedThread)) ;
    WorkQR = baseValue(size_t(LworkQR*usedThread)) ;

    xb_thread = baseValue(size_t(nr*usedThread)) ;

    Hperm  = baseInteger(size_t(N)) ;
    Hswaps = baseInteger(size_t(N)) ;
    Perm   = baseInteger(size_t(nblock*n)) ;

    perm_thread = baseInteger(size_t(n*usedThread)) ;
    iBlock      = baseInteger(size_t(2*usedThread)) ;
    kBlock      = baseInteger(size_t(usedThread)) ;

    // precompute partition for parallel computation
    if ( usedThread > 1 ) {
      iBlock[0] = 0 ;
      iBlock[1] = static_cast<integer>(nblock/usedThread) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt ) {
        iBlock[2*nt+0] = iBlock[2*nt-1]+1 ;
        iBlock[2*nt+1] = static_cast<integer>(((nt+1)*nblock)/usedThread) ;
      }
    } else {
      iBlock[0] = 0 ;
      iBlock[1] = nblock ;
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
    allocate( M.nblock, M.n, M.qr, M.qx, M.nr, M.nx ) ;
    copy( nblock*n_x_nx,     M.Bmat,  1, Bmat,  1 ) ;
    copy( (nblock+1)*nr_x_n, M.Cmat,  1, Cmat,  1 ) ;
    copy( nblock*n_x_n,      M.Dmat,  1, Dmat,  1 ) ;
    copy( nblock*n_x_n,      M.Emat,  1, Emat,  1 ) ;
    copy( nr_x_nx,           M.Fmat,  1, Fmat,  1 ) ;
    copy( (n+qr)*Nc,         M.H0Nqp, 1, H0Nqp, 1 ) ;
    copy( nr*qx,             M.Cqmat, 1, Cqmat, 1 ) ;
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
  BorderedCR<t_Value>::info( std::ostream & stream ) const {
    stream
      << "rows   = " << numRows() << '\n'
      << "cols   = " << numCols() << '\n'
      << "nblock = " << nblock << '\n'
      << "n      = " << n << '\n'
      << "qr     = " << qr << '\n'
      << "qx     = " << qx << '\n'
      << "nr     = " << nr << '\n'
      << "nx     = " << nx << '\n' ;
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
    valueConstPointer H0, integer ld0,
    valueConstPointer HN, integer ldN,
    valueConstPointer Hq, integer ldQ,
    valueConstPointer Hp, integer ldP
  ) {
    integer      m = n + qr ;
    valuePointer H = H0Nqp ;
    gecopy( m, n,  H0, ld0, H, m ) ; H += m * n ;
    gecopy( m, n,  HN, ldN, H, m ) ; H += m * n ;
    gecopy( m, qx, Hq, ldQ, H, m ) ; H += m * qx ;
    gecopy( m, nx, Hp, ldP, H, m ) ;
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
  BorderedCR<t_Value>::addMv( valueConstPointer x, valuePointer res ) const {
    // internal blocks block
    t_Value const * D  = Dmat ;
    t_Value const * E  = Emat ;
    t_Value const * B  = Bmat ;
    t_Value const * xx = x ;
    t_Value const * xe = x  + nblock*n ;
    t_Value const * xq = xe + n ;
    t_Value const * xb = xq + qx ;
    t_Value *       yy = res ;
    for ( integer i = 0 ; i < nblock ; ++i ) {
      gemv( NO_TRANSPOSE, n, n,  1.0, D, n, xx, 1, 1.0, yy, 1 ) ;
      xx += n ;
      gemv( NO_TRANSPOSE, n, n,  1.0, E, n, xx, 1, 1.0, yy, 1 ) ;
      gemv( NO_TRANSPOSE, n, nx, 1.0, B, n, xb, 1, 1.0, yy, 1 ) ;
      yy += n ; D += n_x_n ; E += n_x_n ; B += n_x_nx ;
    }

    integer      m = n+qr ;
    valuePointer H = H0Nqp ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, x,  1, 1.0, yy, 1 ) ; H += m * n ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, xe, 1, 1.0, yy, 1 ) ; H += m * n ;
    gemv( NO_TRANSPOSE, m, qx, 1.0, H, m, xq, 1, 1.0, yy, 1 ) ; H += m * qx ;
    gemv( NO_TRANSPOSE, m, nx, 1.0, H, m, xb, 1, 1.0, yy, 1 ) ;

    if ( nr > 0 ) {
      yy += m ;
      gemv( NO_TRANSPOSE, nr, nx, 1.0, Fmat, nr, xb, 1, 1.0, yy, 1 ) ;
      t_Value const * C = Cmat ;
      xx = x ;
      for ( integer i = 0 ; i <= nblock ; ++i ) {
        gemv( NO_TRANSPOSE, nr, n, 1.0, C, nr, xx, 1, 1.0, yy, 1 ) ;
        xx += n ; C += nr_x_n ;
      }
      gemv( NO_TRANSPOSE, nr, qx, 1.0, Cqmat, nr, xx, 1, 1.0, yy, 1 ) ;
    }
  }

  /*\
   |       _
   |    __| |_   _ _ __ ___  _ __     ___ ___ ___   ___  _ __
   |   / _` | | | | '_ ` _ \| '_ \   / __/ __/ _ \ / _ \| '__|
   |  | (_| | |_| | | | | | | |_) | | (_| (_| (_) | (_) | |
   |   \__,_|\__,_|_| |_| |_| .__/___\___\___\___/ \___/|_|
   |                        |_| |_____|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::dump_ccoord( std::ostream & stream ) const {
    integer ii ;
    integer nnz = nblock*(2*n_x_n+n_x_nx+nr_x_n) + nr_x_n + nr*(qx+nx) + (n+qr)*(2*n+qx+nx) ;
    stream << nnz << '\n' ;

    // bidiagonal
    integer je = nblock*n ;
    integer jq = je+n ;
    integer jb = jq+qx ;
    for ( integer k = 0 ; k < nblock ; ++k ) {
      ii = k*n ;
      valuePointer D = Dmat + k * n_x_n ;
      valuePointer E = Emat + k * n_x_n ;
      valuePointer B = Bmat + k * n_x_nx ;
      for ( integer i = 0 ; i < n ; ++i )
        for ( integer j = 0 ; j < n ; ++j )
          stream
            << ii+i << '\t' << ii+j   << '\t' << D[i+j*n] << '\n'
            << ii+i << '\t' << ii+n+j << '\t' << E[i+j*n] << '\n' ;
      for ( integer i = 0 ; i < n ; ++i )
        for ( integer j = 0 ; j < nx ; ++j )
          stream
            << ii+i << '\t' << jb+j << '\t' << B[i+j*n] << '\n' ;
    }
    integer ie = nblock*n ;
    valuePointer H = H0Nqp ;
    for ( integer i = 0 ; i < n+qr ; ++i )
      for ( integer j = 0 ; j < n ; ++j )
        stream << ie+i << '\t' << j << '\t' << H[i+j*(n+qr)] << '\n' ;

    H += n*(n+qr) ;
    for ( integer i = 0 ; i < n+qr ; ++i )
      for ( integer j = 0 ; j < n+qx+nx ; ++j )
        stream << ie+i << '\t' << je+j << '\t' << H[i+j*(n+qr)] << '\n' ;

    ie += n+qr ;
    for ( integer k = 0 ; k <= nblock ; ++k ) {
      ii = k*n ;
      valuePointer C = Cmat + k * nr_x_n ;
      for ( integer i = 0 ; i < nr ; ++i )
        for ( integer j = 0 ; j < n ; ++j )
          stream
            << ie+i << '\t' << ii+j << '\t' << C[i+j*nr] << '\n' ;
    }
    for ( integer i = 0 ; i < nr ; ++i )
      for ( integer j = 0 ; j < nx ; ++j )
        stream
          << ie+i << '\t' << jb+j << '\t' << Fmat[i+j*nr] << '\n' ;
    jb -= qx ;
    for ( integer i = 0 ; i < nr ; ++i )
      for ( integer j = 0 ; j < qx ; ++j )
        stream
          << ie+i << '\t' << jb+j << '\t' << Cqmat[i+j*nr] << '\n' ;
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
  BorderedCR<t_Value>::factorize() {
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      if ( nr_x_nx > 0 ) zero( nr_x_nx*(usedThread-1), Fmat + nr_x_nx, 1 ) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::factorize_block, this, nt ) ;
      factorize_block(0) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt ) threads[nt].join() ;
      if ( nr_x_nx > 0 )
        for ( integer nt = 1 ; nt < usedThread ; ++nt )
          axpy( nr_x_nx, 1.0, Fmat + nt*nr_x_nx, 1, Fmat, 1 ) ;
      factorize_reduced() ;
    } else {
      factorize_block(0) ;
    }
    #else
    factorize_block(0) ;
    #endif
    load_and_factorize_last() ;
  }

  /*\
    Compute   / L \ (U) = / LU \ = / TOP    \
              \ M /       \ MU /   \ BOTTOM /
  \*/
  template <typename t_Value>
  void
  BorderedCR<t_Value>::buildT( integer           nth,
                               valueConstPointer TOP,
                               valueConstPointer BOTTOM,
                               valuePointer      T,
                               integer *         iperm ) const {
    gecopy( n, n, TOP,    n, T,   n_x_2 ) ;
    gecopy( n, n, BOTTOM, n, T+n, n_x_2 ) ;
    integer info = 0 ;
    switch ( selected ) {
    case BORDERED_LU:
      info = getrf( n_x_2, n, T, n_x_2, iperm ) ;
      break ;
    case BORDERED_QR:
      info = geqrf( n_x_2, n, T, n_x_2, T+2*n_x_n, WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    case BORDERED_QRP:
      { integer * P = perm_thread+nth*n ;
        std::fill( P, P+n, 0 ) ;
        info = geqp3( n_x_2, n, T, n_x_2, P, T+2*n_x_n, WorkQR+nth*LworkQR, LworkQR ) ;
        if ( info == 0 ) permutation_to_exchange( n, P, iperm ) ;
      }
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize INFO = " << info ) ;
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
  BorderedCR<t_Value>::applyT( integer           nth,
                               valueConstPointer T,
                               integer const *   iperm,
                               valuePointer      TOP,
                               integer           ldTOP,
                               valuePointer      BOTTOM,
                               integer           ldBOTTOM,
                               integer           ncol ) const {
    valuePointer W = WorkT + LworkT*nth ;
    gecopy( n, ncol, TOP,    ldTOP,    W,   n_x_2 ) ;
    gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, n_x_2 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = 0 ;
    switch ( selected ) {
    case BORDERED_LU:
      info = swaps( ncol, W, n_x_2, 0, n-1, iperm, 1 ) ;
      ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info ) ;
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, n_x_2, W, n_x_2 ) ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE,
            n, ncol, n,
            -1.0, T+n, n_x_2,    // M
                  W,   n_x_2,    // L^(-1) TOP
             1.0, W+n, n_x_2 ) ; // TOP = BOTTOM - M L^(-1) TOP
      break ;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    n_x_2, ncol, // righe x colonne
                    n,         // numero riflettori usati nel prodotto Q
                    T, n_x_2,
                    T+2*n_x_n,
                    W, n_x_2,
                    WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    }
    gecopy( n, ncol, W+n, n_x_2, TOP,    ldTOP ) ;
    gecopy( n, ncol, W,   n_x_2, BOTTOM, ldBOTTOM ) ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT( integer           nth,
                               valueConstPointer T,
                               integer const *   iperm,
                               valuePointer      TOP,
                               valuePointer      BOTTOM ) const {
    valuePointer W = WorkT + LworkT*nth ;
    copy( n, TOP,    1, W,   1 ) ;
    copy( n, BOTTOM, 1, W+n, 1 ) ;
    integer info = 0 ;
    switch ( selected ) {
    case BORDERED_LU:
      // Apply row interchanges to the right hand sides.
      info = swaps( 1, W, n_x_2, 0, n-1, iperm, 1 ) ;
      if ( info == 0 ) {
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, n_x_2, W, 1 ) ;
        gemv( NO_TRANSPOSE,
              n, n,
              -1.0, T+n, n_x_2,
                    W,   1,
               1.0, W+n, 1 ) ;
      }
      break ;
    case BORDERED_QR:
    case BORDERED_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    n_x_2, 1, // righe x colonne
                    n,      // numero riflettori usati nel prodotto Q
                    T, n_x_2,
                    T+2*n_x_n,
                    W, n_x_2,
                    WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info ) ;
    copy( n, W+n, 1, TOP,    1 ) ;
    copy( n, W,   1, BOTTOM, 1 ) ;
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
    integer iblock = iBlock[2*nth+0] ;
    integer eblock = iBlock[2*nth+1] ;
    integer nblk   = eblock - iblock ;

    valuePointer Bmat0 = Bmat + iblock*n_x_nx ;
    valuePointer Cmat0 = Cmat + iblock*nr_x_n ;
    valuePointer Dmat0 = Dmat + iblock*n_x_n ;
    valuePointer Emat0 = Emat + iblock*n_x_n ;
    valuePointer T0    = Tmat + iblock*Tsize ;
    integer *    P0    = Perm + iblock*n ;

    valuePointer Fmat_th = Fmat + nth * nr_x_nx ;

    integer k = 1 ;
    while ( k < nblk ) {

      valuePointer Bjp = Bmat0 ;
      valuePointer Bj  = Bmat0 + k*n_x_nx ;
      valuePointer Cjp = Cmat0 ;
      valuePointer Cj  = Cmat0 + k*nr_x_n ;
      valuePointer Djp = Dmat0 ;
      valuePointer Dj  = Dmat0 + k*n_x_n ;
      valuePointer Ejp = Emat0 ;
      valuePointer Ej  = Emat0 + k*n_x_n ;
      valuePointer T   = T0    + k*Tsize ;
      integer *    P   = P0    + k*n ;

      // -----------------------------------------------------------------------
      integer k_x_2 = 2*k ;
      for ( integer j = iblock+k ; j < eblock ; j += k_x_2 ) {

        buildT( nth, Ejp, Dj, T, P ) ;

        zero( n_x_n, Dj, 1 ) ;
        applyT( nth, T, P, Djp, n, Dj, n, n ) ;

        zero( n_x_n, Ejp, 1 ) ;
        applyT( nth, T, P, Ejp, n, Ej, n, n ) ;

        if ( nx > 0 ) applyT( nth, T, P, Bjp, n, Bj, n, nx ) ;

        if ( nr > 0 ) {

          if ( selected == BORDERED_QRP ) {
            integer i = n ;
            do {
              --i ;
              if ( P[i] > i ) swap( nr, Cj+i*nr, 1, Cj+P[i]*nr, 1 ) ;
            } while ( i > 0 ) ;
          }
          trsm( RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
                nr, n, 1.0, T, n_x_2, Cj, nr ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, n, n,
                -1.0, Cj,  nr,
                      Dj,  n,
                 1.0, Cjp, nr ) ;

          integer      jpp = std::min(j+k,eblock) ;
          valuePointer Cpp = Cmat + jpp*nr_x_n ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, n, n,
                -1.0, Cj,  nr,
                      Ej,  n,
                 1.0, Cpp, nr ) ;
        }

        if ( nr_x_nx > 0 )
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, nx, n,
                -1.0, Cj,      nr,
                      Bj,      n,
                 1.0, Fmat_th, nr ) ; // solo accumulo!

        // NEXT STEP
        T   += k_x_2*Tsize ;
        P   += k_x_2*n ;
        Djp += k_x_2*n_x_n ;
        Dj  += k_x_2*n_x_n ;
        Ejp += k_x_2*n_x_n ;
        Ej  += k_x_2*n_x_n ;
        Bj  += k_x_2*n_x_nx ;
        Bjp += k_x_2*n_x_nx ;
        Cj  += k_x_2*nr_x_n ;
        Cjp += k_x_2*nr_x_n ;
      }
      k *= 2 ;
    }
    kBlock[nth] = k ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_reduced() {
    integer nblk = 2*usedThread-1 ;
    integer k = 1 ;
    while ( k < nblk ) {
      // -----------------------------------------------------------------------
      for ( integer jj = k ; jj < nblk ; jj += 2*k ) {
        integer j  = iBlock[jj] ;
        integer jp = iBlock[jj-k] ;
        valuePointer T   = Tmat + j*Tsize ;
        integer *    P   = Perm + j*n ;
        valuePointer Djp = Dmat + jp*n_x_n ;
        valuePointer Dj  = Dmat + j*n_x_n ;
        valuePointer Ejp = Emat + jp*n_x_n ;
        valuePointer Ej  = Emat + j*n_x_n ;

        buildT( 0, Ejp, Dj, T, P ) ;

        zero( n_x_n, Dj, 1 ) ;
        applyT( 0, T, P, Djp, n, Dj, n, n ) ;

        zero( n_x_n, Ejp, 1 ) ;
        applyT( 0, T, P, Ejp, n, Ej, n, n ) ;

        if ( nx > 0 ) {
          valuePointer Bj  = Bmat + j*n_x_nx ;
          valuePointer Bjp = Bmat + jp*n_x_nx ;
          applyT( 0, T, P, Bjp, n, Bj, n, nx ) ;
        }

        if ( nr > 0 ) {
          valuePointer Cj  = Cmat + j*nr_x_n ;
          valuePointer Cjp = Cmat + jp*nr_x_n ;
          if ( selected == BORDERED_QRP ) {
            integer i = n ;
            do {
              --i ;
              if ( P[i] > i ) swap( nr, Cj+i*nr, 1, Cj+P[i]*nr, 1 ) ;
            } while ( i > 0 ) ;
          }
          trsm( RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
                nr, n, 1.0, T, n_x_2, Cj, nr ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, n, n,
                -1.0, Cj,  nr,
                      Dj,  n,
                 1.0, Cjp, nr ) ;

          integer      jpp = iBlock[std::min(jj+k,nblk)] ;
          valuePointer Cpp = Cmat + jpp*nr_x_n ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, n, n,
                -1.0, Cj,  nr,
                      Ej,  n,
                 1.0, Cpp, nr ) ;

        }

        if ( nr_x_nx > 0 ) {
          valuePointer Cj = Cmat + j*nr_x_n ;
          valuePointer Bj = Bmat + j*n_x_nx ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, nx, n,
                -1.0, Cj,   nr,
                      Bj,   n,
                 1.0, Fmat, nr ) ;
        }
      }
      k *= 2 ;
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
    valuePointer Cnb = Cmat + nblock*nr_x_n ;
    valuePointer W0  = Hmat ;
    valuePointer WN  = W0+n*Nr ;
    valuePointer Wq  = WN+n*Nr ;
    valuePointer Wp  = Wq+qx*Nr ;

    gecopy( n,  n,  Dmat, n, W0, Nr ) ;
    gecopy( n,  n,  Emat, n, WN, Nr ) ;
    gezero( n,  qx,          Wq, Nr ) ;
    gecopy( n,  nx, Bmat, n, Wp, Nr ) ;

    gecopy( n+qr, Nc, H0Nqp, n+qr, Hmat+n, Nr ) ;

    integer offs = n_x_2+qr ;

    gecopy( nr, n,  Cmat,  nr, W0+offs, Nr ) ;
    gecopy( nr, n,  Cnb,   nr, WN+offs, Nr ) ;
    gecopy( nr, qx, Cqmat, nr, Wq+offs, Nr ) ;
    gecopy( nr, nx, Fmat,  nr, Wp+offs, Nr ) ;

    integer info = 0 ;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      info = getrf( Nr, Nc, Hmat, Nr, Hperm ) ;
      break ;
    case BORDERED_LAST_LUP:
      info = getc2( Nr, Hmat, Nr, Hperm, Hswaps ) ;
      break ;
    case BORDERED_LAST_QR:
      info = geqrf( Nr, Nc, Hmat, Nr, Htau, Work, Lwork ) ;
      break ;
    case BORDERED_LAST_QRP:
      std::fill( Hperm, Hperm+Nc, 0 ) ;
      info = geqp3( Nr, Nc, Hmat, Nr, Hperm, Htau, Work, Lwork ) ;
      // convert permutation to exchanges
      if ( info == 0 ) permutation_to_exchange( Nr, Hperm, Hswaps ) ;
      break ;
    }
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize_last INFO = " << info <<
                   " Nr = " << Nr << " Nc = " << Nc ) ;
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
  BorderedCR<t_Value>::solve_last( valuePointer x ) const {
    valuePointer X = x + (nblock-1)*n ;
    swap( n, X, 1, x, 1 ) ; // uso x stesso come temporaneo
    integer info = 0 ;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      info = getrs( NO_TRANSPOSE, Nr, 1, Hmat, Nr, Hperm, X, Nr ) ;
      break ;
    case BORDERED_LAST_LUP:
      {
        valueType scale = gesc2( Nr, Hmat, Nr, X, Hperm, Hswaps ) ;
        ALGLIN_ASSERT( isZero(scale-1), "BorderedCR::solve_last scale = " << scale ) ;
      }
      break ;
    case BORDERED_LAST_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    Nr, 1, // righe x colonne
                    Nr-1,  // numero riflettori usati nel prodotto Q
                    Hmat, Nr /*ldA*/,
                    Htau,
                    X, Nr,
                    Work, Lwork ) ;
      if ( info == 0 ) {
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, Nc, Hmat, Nr, X, 1 ) ;
        if ( last_selected == BORDERED_LAST_QRP ) {
          for ( integer i = 0 ; i < Nr ; ++i )
            if ( Hswaps[i] > i )
              std::swap( X[i], X[Hswaps[i]] ) ;
          //integer i = n ;
          //do {
          //  --i ;
          //  if ( Hswaps[i] > i ) std::swap( X[i], X[Hswaps[i]] ) ;
          //} while ( i > 0 ) ;
        }
      }
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info ) ;
    swap( n, X, 1, x, 1 ) ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last( integer      nrhs,
                                   valuePointer x,
                                   integer      ldX ) const {
    valuePointer X = x + (nblock-1)*n ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
    integer info = 0 ;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      info = getrs( NO_TRANSPOSE, Nr, nrhs, Hmat, Nr, Hperm, X, ldX ) ;
      break ;
    case BORDERED_LAST_LUP:
      for ( integer i = 0 ; i < nrhs && info == 0 ; ++i ) {
        valueType scale = gesc2( Nr, Hmat, Nr, X+i*ldX, Hperm, Hswaps ) ;
        ALGLIN_ASSERT( isZero(scale-1), "BorderedCR::solve_last scale = " << scale ) ;
      }
      break ;
    case BORDERED_LAST_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    Nr, nrhs, // righe x colonne
                    Nr-1,     // numero riflettori usati nel prodotto Q
                    Hmat, Nr,
                    Htau,
                    X, ldX,
                    Work, Lwork ) ;
      if ( info == 0 ) { // A P P^T x = b --> Q R P^T x = b --> x = P R^(-1) Q^T b
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
              Nr, nrhs, 1.0, Hmat, Nr, X, ldX ) ;
        if ( last_selected == BORDERED_LAST_QRP ) {
          for ( integer i = 0 ; i < Nr ; ++i )
            if ( Hswaps[i] > i )
              swap( nrhs, X+i, ldX, X+Hswaps[i], ldX ) ;
          //integer i = n ;
          //do {
          //  --i ;
          //  if ( Hswaps[i] > i ) swap( nrhs, X+i, ldX, X+Hswaps[i], ldX ) ;
          //} while ( i > 0 ) ;
        }
      }
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info ) ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
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
  BorderedCR<t_Value>::solve( valuePointer x ) const {
    valuePointer xb = x + (nblock+1)*n + qr ; // deve essere b!
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      if ( nr > 0 ) {
        zero( nr*usedThread, xb_thread, 1 ) ;
        for ( integer nt = 1 ; nt < usedThread ; ++nt )
          threads[nt] = std::thread( &BorderedCR<t_Value>::forward, this,
                                     nt, x, xb_thread+nr*nt ) ;
        forward(0,x,xb) ;
        for ( integer nt = 1 ; nt < usedThread ; ++nt ) {
          threads[nt].join() ;
          axpy( nr, 1.0, xb_thread+nr*nt, 1, xb, 1 ) ;
        }
      } else {
        for ( integer nt = 1 ; nt < usedThread ; ++nt )
          threads[nt] = std::thread( &BorderedCR<t_Value>::forward, this,
                                     nt, x, xb_thread+nr*nt ) ;
        forward(0,x,xb) ;
        for ( integer nt = 1 ; nt < usedThread ; ++nt ) threads[nt].join() ;
      }
      forward_reduced(x,xb) ;
    } else {
      forward(0,x,xb) ;
    }
    #else
      forward(0,x,xb) ;
    #endif

    solve_last( x ) ;

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      backward_reduced(x) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::backward, this, nt, x ) ;
      backward(0,x) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt ) threads[nt].join() ;
    } else {
      backward(0,x) ;
    }
    #else
    backward(0,x) ;
    #endif
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve( integer      nrhs,
                              valuePointer rhs,
                              integer      ldRhs ) const {
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      for ( integer nt = 1 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::forward_n, this,
                                   nt, nrhs, rhs, ldRhs ) ;
      forward_n( 0, nrhs, rhs, ldRhs ) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt ) threads[nt].join() ;
      forward_n_reduced( nrhs, rhs, ldRhs ) ;
    } else {
      forward_n( 0, nrhs, rhs, ldRhs ) ;
    }
    #else
      forward_n( 0, nrhs, rhs, ldRhs ) ;
    #endif

    solve_last( nrhs, rhs, ldRhs ) ;

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 1 ) {
      backward_n_reduced( nrhs, rhs, ldRhs ) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BorderedCR<t_Value>::backward_n, this,
                                   nt, nrhs, rhs, ldRhs ) ;
      backward_n( 0, nrhs, rhs, ldRhs ) ;
      for ( integer nt = 1 ; nt < usedThread ; ++nt ) threads[nt].join() ;
    } else {
      backward_n( 0, nrhs, rhs, ldRhs ) ;
    }
    #else
      backward_n( 0, nrhs, rhs, ldRhs ) ;
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
  BorderedCR<t_Value>::forward( integer nth, valuePointer x, valuePointer xb ) const {
    integer iblock = iBlock[2*nth+0] ;
    integer eblock = iBlock[2*nth+1] ;
    integer nblk   = eblock - iblock ;
    valuePointer x0 = x    + iblock*n ;
    valuePointer T0 = Tmat + iblock*Tsize ;
    integer *    P0 = Perm + iblock*n ;
    valuePointer C0 = Cmat + iblock*nr_x_n ;

    integer k = 1 ;
    while ( k < nblk ) {
      valuePointer xj  = x0 + k*n ;
      valuePointer xjp = x0 ;
      valuePointer T   = T0 + k*Tsize ;
      integer *    P   = P0 + k*n ;
      valuePointer Cj  = C0 + k*nr_x_n ;
      integer    k_x_2 = 2*k ;
      for ( integer jj = k ; jj < nblk ; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, xj ) ;
        if ( nr > 0 ) gemv( NO_TRANSPOSE, nr, n, -1.0, Cj, nr, xj, 1, 1.0, xb, 1 ) ; // solo accumulato
        xj  += k_x_2*n ;
        xjp += k_x_2*n ;
        T   += k_x_2*Tsize ;
        P   += k_x_2*n ;
        Cj  += k_x_2*nr_x_n ;
      }
      k *= 2 ;
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_n( integer      nth,
                                  integer      nrhs,
                                  valuePointer x,
                                  integer      ldX ) const {
    valuePointer xb = x + (nblock+1)*n + qr ;
    integer iblock = iBlock[2*nth+0] ;
    integer eblock = iBlock[2*nth+1] ;
    integer nblk   = eblock - iblock ;
    valuePointer x0 = x    + iblock*n ;
    valuePointer T0 = Tmat + iblock*Tsize ;
    integer *    P0 = Perm + iblock*n ;
    valuePointer C0 = Cmat + iblock*nr_x_n ;

    integer k = 1 ;
    while ( k < nblk ) {
      valuePointer xj  = x0 + k*n ;
      valuePointer xjp = x0 ;
      valuePointer T   = T0 + k*Tsize ;
      integer *    P   = P0 + k*n ;
      valuePointer Cj  = C0 + k*nr_x_n ;
      integer    k_x_2 = 2*k ;

      for ( integer jj = k ; jj < nblk ; jj += k_x_2 ) {
        applyT( nth, T, P, xjp, ldX, xj, ldX, nrhs ) ;
        if ( nr > 0 ) {
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.lock() ;
          #endif
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, nrhs, n,
                -1.0, Cj, nr,
                      xj, ldX,
                 1.0, xb, ldX ) ;
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.unlock() ;
          #endif
        }
        xj  += k_x_2*n ;
        xjp += k_x_2*n ;
        T   += k_x_2*Tsize ;
        P   += k_x_2*n ;
        Cj  += k_x_2*nr_x_n ;
      }
      k *= 2 ;
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_reduced( valuePointer x, valuePointer xb ) const {
    integer nblk = 2*usedThread-1 ;
    integer k = 1 ;
    while ( k < nblk ) {
      for ( integer jj = k ; jj < nblk ; jj += 2*k ) {
        integer j  = iBlock[jj] ;
        integer jp = iBlock[jj-k] ;
        valuePointer const T = Tmat + j*Tsize ;
        integer const *    P = Perm + j*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xjp = x + jp*n ;
        applyT( 0, T, P, xjp, xj ) ;
        if ( nr > 0 ) {
          valuePointer Cj = Cmat + j*nr_x_n ;
          gemv( NO_TRANSPOSE, nr, n, -1.0, Cj, nr, xj, 1, 1.0, xb, 1 ) ;
        }
      }
      k *= 2 ;
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::forward_n_reduced( integer      nrhs,
                                          valuePointer x,
                                          integer      ldX ) const {
    valuePointer xb = x + (nblock+1)*n + qr ;
    integer nblk = 2*usedThread-1 ;
    integer k = 1 ;
    while ( k < nblk ) {
      for ( integer jj = k ; jj < nblk ; jj += 2*k ) {
        integer j  = iBlock[jj] ;
        integer jp = iBlock[jj-k] ;
        valuePointer const T = Tmat + j*Tsize ;
        integer const *    P = Perm + j*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xjp = x + jp*n ;
        applyT( 0, T, P, xjp, ldX, xj, ldX, nrhs ) ;
        if ( nr > 0 ) {
          valuePointer Cj = Cmat + j*nr_x_n ;
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.lock() ;
          #endif
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nr, nrhs, n,
                -1.0, Cj, nr,
                      xj, ldX,
                 1.0, xb, ldX ) ;
          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          spin.unlock() ;
          #endif
        }
      }
      k *= 2 ;
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
  BorderedCR<t_Value>::backward( integer nth, valuePointer x ) const {
    valuePointer xn = x + (nblock+1)*n + qx ;
    integer iblock = iBlock[2*nth+0] ;
    integer eblock = iBlock[2*nth+1] ;
    valuePointer x0 = x    + iblock*n ;
    valuePointer B0 = Bmat + iblock*n_x_nx ;
    valuePointer D0 = Dmat + iblock*n_x_n ;
    valuePointer E0 = Emat + iblock*n_x_n ;
    valuePointer T0 = Tmat + iblock*Tsize ;
    integer k = kBlock[nth] ;
    while ( (k/=2) > 0 ) {
      valuePointer xj = x0 + k*n ;
      valuePointer xp = x0 ;
      valuePointer Bj = B0 + k*n_x_nx ;
      valuePointer Dj = D0 + k*n_x_n ;
      valuePointer Ej = E0 + k*n_x_n ;
      valuePointer T  = T0 + k*Tsize ;
      integer   k_x_2 = 2*k ;
      for ( integer j = iblock+k ; j < eblock ; j += k_x_2 ) {
        integer      jpp = std::min(j+k,eblock) ;
        valuePointer xpp = x + jpp*n ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 ) ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 ) ;
        if ( nx > 0 ) gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 ) ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 ) ;
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] ) ;
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx ;
        Dj += k_x_2*n_x_n ;
        Ej += k_x_2*n_x_n ;
        xj += k_x_2*n ;
        xp += k_x_2*n ;
        T  += k_x_2*Tsize ;
      }
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n( integer      nth,
                                   integer      nrhs,
                                   valuePointer x,
                                   integer      ldX ) const {
    valuePointer xn = x + (nblock+1)*n + qx ;
    integer iblock = iBlock[2*nth+0] ;
    integer eblock = iBlock[2*nth+1] ;
    valuePointer x0 = x    + iblock*n ;
    valuePointer B0 = Bmat + iblock*n_x_nx ;
    valuePointer D0 = Dmat + iblock*n_x_n ;
    valuePointer E0 = Emat + iblock*n_x_n ;
    valuePointer T0 = Tmat + iblock*Tsize ;
    integer k = kBlock[nth] ;
    while ( (k/=2) > 0 ) {
      valuePointer xj = x0 + k*n ;
      valuePointer xp = x0 ;
      valuePointer Bj = B0 + k*n_x_nx ;
      valuePointer Dj = D0 + k*n_x_n ;
      valuePointer Ej = E0 + k*n_x_n ;
      valuePointer T  = T0 + k*Tsize ;
      integer   k_x_2 = 2*k ;
      for ( integer j = iblock+k ; j < eblock ; j += k_x_2 ) {
        integer      jpp = std::min(j+k,eblock) ;
        valuePointer xpp = x + jpp*n ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Dj,  n,
                    xp,  ldX,
               1.0, xj,  ldX ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Ej,  n,
                    xpp, ldX,
               1.0, xj,  ldX ) ;
        if ( nx > 0 )
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                n, nrhs, nx,
                -1.0, Bj, n,
                      xn, ldX,
                 1.0, xj, ldX ) ;
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
              n, nrhs, 1.0, T, n_x_2, xj, ldX ) ;
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              swap( nrhs, xj+i, ldX, xj+P[i], ldX ) ;
        }
        // NEXT STEP
        Bj += k_x_2*n_x_nx ;
        Dj += k_x_2*n_x_n ;
        Ej += k_x_2*n_x_n ;
        xj += k_x_2*n ;
        xp += k_x_2*n ;
        T  += k_x_2*Tsize ;
      }
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_reduced( valuePointer x ) const {
    valuePointer xn = x + (nblock+1)*n + qx ;
    integer nblk = 2*usedThread-1 ;
    integer k = 1 ;
    while ( k < nblk ) k *= 2 ;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k ; jj < nblk ; jj += 2*k ) {
        integer      j   = iBlock[jj] ;
        integer      jp  = iBlock[jj-k] ;
        integer      jpp = iBlock[std::min(jj+k,nblk)] ;
        valuePointer Dj  = Dmat + j*n_x_n ;
        valuePointer Ej  = Emat + j*n_x_n ;
        valuePointer xj  = x + j*n ;
        valuePointer xp  = x + jp*n ;
        valuePointer xpp = x + jpp*n ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 ) ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 ) ;
        if ( nx > 0 ) {
          valuePointer Bj = Bmat + j*n_x_nx ;
          gemv( NO_TRANSPOSE, n, nx, -1.0, Bj, n, xn, 1, 1.0, xj, 1 ) ;
        }
        valuePointer const T = Tmat + j*Tsize ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, n_x_2, xj, 1 ) ;
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] ) ;
        }
      }
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::backward_n_reduced( integer      nrhs,
                                           valuePointer x,
                                           integer      ldX ) const {
    valuePointer xn = x + (nblock+1)*n + qx ;
    integer nblk = 2*usedThread-1 ;
    integer k = 1 ;
    while ( k < nblk ) k *= 2 ;
    while ( (k/=2) > 0 ) {
      for ( integer jj = k ; jj < nblk ; jj += 2*k ) {
        integer      j   = iBlock[jj] ;
        integer      jp  = iBlock[jj-k] ;
        integer      jpp = iBlock[std::min(jj+k,nblk)] ;
        valuePointer Dj  = Dmat + j*n_x_n ;
        valuePointer Ej  = Emat + j*n_x_n ;
        valuePointer xj  = x + j*n ;
        valuePointer xjp = x + jp*n ;
        valuePointer xpp = x + jpp*n ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Dj,  n,
                    xjp, ldX,
               1.0, xj,  ldX ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Ej,  n,
                    xpp, ldX,
               1.0, xj,  ldX ) ;
        if ( nx > 0 ) {
          valuePointer Bj = Bmat + j*n_x_nx ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                n, nrhs, nx,
                -1.0, Bj, n,
                      xn, ldX,
                 1.0, xj, ldX ) ;
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        valuePointer const T = Tmat + j*Tsize ;
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
              n, nrhs, 1.0, T, n_x_2, xj, ldX ) ;
        if ( selected == BORDERED_QRP ) {
          integer const * P = Perm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              swap( nrhs, xj+i, ldX, xj+P[i], ldX ) ;
        }
      }
    }
  }

  template class BorderedCR<float> ;
  template class BorderedCR<double> ;

}
