/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#include "BorderedCR.hh"

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
                                 integer _q,
                                 integer _nb ) {

    if ( nblock == _nblock && n == _n && q == _q && nb == _nb ) return ;

    nblock = _nblock ;
    n      = _n ;
    q      = _q ;
    nb     = _nb ;
    nx2    = n*2 ;
    nxn    = n*n ;
    nxnb   = n*nb ;
    N      = nx2+nb+q ;
    Tsize  = 2*nxn+n ;

    Lwork = max(N,2*max(nxn,nxnb)) ;

    valueType tmp ; // get optimal allocation
    integer info = geqrf( N, N, nullptr, N, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    info = geqp3( N, N, nullptr, N, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    LworkT = 2*nxn ;
    info = geqrf( nx2, n, nullptr, nx2, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    LworkQR = integer(tmp) ;

    info = geqp3( nx2, n, nullptr, nx2, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( LworkQR < integer(tmp) ) LworkQR = integer(tmp) ;

    integer nnz = Lwork+(LworkT+LworkQR)*numThread+N*(N+1)+nb*(nb+q)+nxnb+nblock*(2*nxnb+2*nxn+Tsize+n)+(n+q)*(nx2+nb+q) ;
    baseValue.allocate(size_t(nnz)) ;
    baseInteger.allocate(size_t(2*N+nblock*(3*n+1))) ;

    Bmat  = baseValue(size_t(nblock*nxnb)) ;
    Cmat  = baseValue(size_t((nblock+1)*nxnb)) ;
    Cqmat = baseValue(size_t(nb*q)) ;
    Dmat  = baseValue(size_t(nblock*nxn)) ;
    Emat  = baseValue(size_t(nblock*nxn)) ;
    Fmat  = baseValue(size_t(nb*nb)) ;
    
    Tmat  = baseValue(size_t(nblock*Tsize)) ;
    Ttau  = baseValue(size_t(nblock*n)) ;
    Hmat  = baseValue(size_t(N*N)) ;
    Htau  = baseValue(size_t(N)) ;
    H0Nqp = baseValue(size_t((n+q)*N)) ;

    Work   = baseValue(size_t(Lwork)) ;
    WorkT  = baseValue(size_t(LworkT*numThread)) ;
    WorkQR = baseValue(size_t(LworkQR*numThread)) ;

    Hperm  = baseInteger(size_t(N)) ;
    Hswaps = baseInteger(size_t(N)) ;
    Rperm  = baseInteger(size_t(nblock*n)) ;
    Cperm  = baseInteger(size_t(nblock*n)) ;
    
    task_done = baseInteger(size_t(nblock)) ;
    Tperm     = baseInteger(size_t(nblock*n)) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::fillZero() {
    zero( nblock*nxnb, Bmat, 1 ) ;
    zero( (nblock+1)*nxnb, Cmat, 1 ) ;
    zero( nb*q, Cqmat, 1 ) ;
    zero( nblock*nxn, Dmat, 1 ) ;
    zero( nblock*nxn, Emat, 1 ) ;
    zero( nb*nb, Fmat, 1 ) ;
    zero( N*N, Hmat, 1 ) ;
    zero( (n+q)*N, H0Nqp, 1 ) ;
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
    integer m = n + q ;
    valuePointer H = H0Nqp ;
    gecopy( m, n,  H0, ld0, H, m ) ; H += m * n ;
    gecopy( m, n,  HN, ldN, H, m ) ; H += m * n ;
    gecopy( m, q,  Hq, ldQ, H, m ) ; H += m * q ;
    gecopy( m, nb, Hp, ldP, H, m ) ;
  }

  /*\
   |   __  __
   |  |  \/  |_   __
   |  | |\/| \ \ / /
   |  | |  | |\ V /
   |  |_|  |_| \_/
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::Mv( valueConstPointer x, valuePointer res ) const {
    zero( n*(nblock+1)+nb+q, res, 1 ) ;
    // internal blocks block
    t_Value const * D  = Dmat ;
    t_Value const * E  = Emat ;
    t_Value const * B  = Bmat ;
    t_Value const * xx = x ;
    t_Value const * xe = x  + nblock*n ;
    t_Value const * xq = xe + n ;
    t_Value const * xb = xq + q ;
    t_Value *       yy = res ;
    for ( integer i = 0 ; i < nblock ; ++i ) {
      gemv( NO_TRANSPOSE, n, n,  1.0, D, n, xx, 1, 0.0, yy, 1 ) ;
      xx += n ;
      gemv( NO_TRANSPOSE, n, n,  1.0, E, n, xx, 1, 1.0, yy, 1 ) ;
      if ( nb > 0 ) gemv( NO_TRANSPOSE, n, nb, 1.0, B, n, xb, 1, 1.0, yy, 1 ) ;
      yy += n ; D += nxn ; E += nxn ; B += nxnb ;
    }

    integer      m = n+q ;
    valuePointer H = H0Nqp ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, x,  1, 0.0, yy, 1 ) ; H += m * n ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H, m, xe, 1, 1.0, yy, 1 ) ; H += m * n ;
    gemv( NO_TRANSPOSE, m, q,  1.0, H, m, xq, 1, 1.0, yy, 1 ) ; H += m * q ;
    gemv( NO_TRANSPOSE, m, nb, 1.0, H, m, xb, 1, 1.0, yy, 1 ) ;

    if ( nb > 0 ) {
      yy += m ;
      gemv( NO_TRANSPOSE, nb, nb, 1.0, Fmat, nb, xb, 1, 0.0, yy, 1 ) ;
      t_Value const * C = Cmat ;
      xx = x ;
      for ( integer i = 0 ; i <= nblock ; ++i ) {
        gemv( NO_TRANSPOSE, nb, n, 1.0, C, nb, xx, 1, 1.0, yy, 1 ) ;
        xx += n ; C += nxnb ;
      }
      gemv( NO_TRANSPOSE, nb, q, 1.0, Cqmat, nb, xx, 1, 1.0, yy, 1 ) ;
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
    integer nnz = 2*nblock*(nxn+nxnb) + nxnb + nb*nb + (n+q)*N ;
    stream << nnz << '\n' ;

    // bidiagonal
    integer je = nblock*n ;
    integer jq = je+n ;
    integer jb = jq+q ;
    for ( integer k = 0 ; k < nblock ; ++k ) {
      ii = k*n ;
      valuePointer D = Dmat + k * nxn ;
      valuePointer E = Emat + k * nxn ;
      valuePointer B = Bmat + k * nxnb ;
      for ( integer i = 0 ; i < n ; ++i )
        for ( integer j = 0 ; j < n ; ++j )
          stream
            << ii+i << '\t' << ii+j   << '\t' << D[i+j*n] << '\n'
            << ii+i << '\t' << ii+n+j << '\t' << E[i+j*n] << '\n' ;
      for ( integer i = 0 ; i < n ; ++i )
        for ( integer j = 0 ; j < nb ; ++j )
          stream
            << ii+i << '\t' << jb+j << '\t' << B[i+j*n] << '\n' ;
    }
    integer ie = nblock*n ;
    valuePointer H = H0Nqp ;
    for ( integer i = 0 ; i < n+q ; ++i )
      for ( integer j = 0 ; j < n ; ++j )
        stream << ie+i << '\t' << j << '\t' << H[i+j*(n+q)] << '\n' ;

    H += n*(n+q) ;
    for ( integer i = 0 ; i < n+q ; ++i )
      for ( integer j = 0 ; j < n+q+nb ; ++j )
        stream << ie+i << '\t' << je+j << '\t' << H[i+j*(n+q)] << '\n' ;

    ie += n+q ;
    for ( integer k = 0 ; k <= nblock ; ++k ) {
      ii = k*n ;
      valuePointer C = Cmat + k * nxnb ;
      for ( integer i = 0 ; i < nb ; ++i )
        for ( integer j = 0 ; j < n ; ++j )
          stream
            << ie+i << '\t' << ii+j << '\t' << C[i+j*nb] << '\n' ;
    }
    for ( integer i = 0 ; i < nb ; ++i )
      for ( integer j = 0 ; j < nb ; ++j )
        stream
          << ie+i << '\t' << jb+j << '\t' << Fmat[i+j*nb] << '\n' ;
    jb -= q ;
    for ( integer i = 0 ; i < nb ; ++i )
      for ( integer j = 0 ; j < q ; ++j )
        stream
          << ie+i << '\t' << jb+j << '\t' << Cqmat[i+j*nb] << '\n' ;
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
    std::fill( task_done, task_done + nblock, 0 ) ;
    barrier.setup(numThread) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt] = std::thread( &BorderedCR<t_Value>::factorize_mt, this, nt ) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt].join() ;
    #else
    factorize_mt(0) ;
    #endif
    load_last_block() ;
    factorize_last() ;
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
    gecopy( n, n, TOP,    n, T,   nx2 ) ;
    gecopy( n, n, BOTTOM, n, T+n, nx2 ) ;
    integer info = 0 ;
    switch ( selected ) {
    case BORDERED_LU:
      info = getrf( nx2, n, T, nx2, iperm ) ;
      break ;
    case BORDERED_QR:
      info = geqrf( nx2, n, T, nx2, T+2*nxn, WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    case BORDERED_LAST_QRP:
      { integer * P = Tperm+nth*n ;
        info = geqp3( nx2, n, T, nx2, P, T+2*nxn, WorkQR+nth*LworkQR, LworkQR ) ;
        // convert permutation to exchanges
        if ( info == 0 )
          for ( integer i = 0 ; i < n ; ++i ) {
            integer j = i ;
            while ( j < n ) { if ( P[j] == i+1 ) break ; ++j ; }
            std::swap( P[j], P[i] ) ;
            iperm[i] = j ;
          }
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
    gecopy( n, ncol, TOP,    ldTOP,    W,   nx2 ) ;
    gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, nx2 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = 0 ;
    switch ( selected ) {
    case BORDERED_LU:
      info = swaps( ncol, W, nx2, 0, n-1, iperm, 1 ) ;
      ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info ) ;
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, nx2, W, nx2 ) ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE,
            n, ncol, n,
            -1.0, T+n, nx2,    // M
                  W,   nx2,    // L^(-1) TOP
             1.0, W+n, nx2 ) ; // TOP = BOTTOM - M L^(-1) TOP
      break ;
    case BORDERED_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    nx2, ncol, // righe x colonne
                    n,         // numero riflettori usati nel prodotto Q
                    T, nx2,
                    T+2*nxn,
                    W, nx2,
                    WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    }
    gecopy( n, ncol, W+n, nx2, TOP,    ldTOP ) ;
    gecopy( n, ncol, W,   nx2, BOTTOM, ldBOTTOM ) ;
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
      info = swaps( 1, W, nx2, 0, n-1, iperm, 1 ) ;
      if ( info == 0 ) {
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, nx2, W, 1 ) ;
        gemv( NO_TRANSPOSE,
              n, n,
              -1.0, T+n, nx2,
                    W,   1,
               1.0, W+n, 1 ) ;
      }
      break ;
    case BORDERED_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    nx2, 1, // righe x colonne
                    n,      // numero riflettori usati nel prodotto Q
                    T, nx2,
                    T+2*nxn,
                    W, nx2,
                    WorkQR+nth*LworkQR, LworkQR ) ;
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::applyT INFO = " << info ) ;
    copy( n, W+n, 1, TOP,    1 ) ;
    copy( n, W,   1, BOTTOM, 1 ) ;
  }
  
  /*\
   |    __            _             _                       _
   |   / _| __ _  ___| |_ ___  _ __(_)_______     _ __ ___ | |_
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \   | '_ ` _ \| __|
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/   | | | | | | |_
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|___|_| |_| |_|\__|
   |                                        |_____|
  \*/
  
  #define CHECK_TASK_DONE(FLG)       \
  bool skip = false ;                \
  spin.lock() ;                      \
  skip = (task_done[j]&FLG) == FLG ; \
  if ( !skip ) task_done[j] |= FLG ; \
  spin.unlock() ;                    \
  if ( skip ) continue

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_mt( integer nth ) {
    for ( integer k = 1 ; k < nblock ; k = 2*k ) {
      // -----------------------------------------------------------------------
      for ( integer j = k ; j < nblock ; j += 2*k ) {

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        CHECK_TASK_DONE(0x01);
        #endif

        integer jp = j-k ;
        valuePointer T   = Tmat  + j*Tsize ;
        integer *    P   = Rperm + j*n ;
        valuePointer Djp = Dmat  + jp*nxn ;
        valuePointer Dj  = Dmat  + j*nxn ;
        valuePointer Ejp = Emat  + jp*nxn ;
        valuePointer Ej  = Emat  + j*nxn ;

        buildT( nth, Ejp, Dj, T, P ) ;

        zero( nxn, Dj, 1 ) ;
        applyT( nth, T, P, Djp, n, Dj, n, n ) ;

        zero( nxn, Ejp, 1 ) ;
        applyT( nth, T, P, Ejp, n, Ej, n, n ) ;

        if ( nb > 0 ) {

          valuePointer Bj  = Bmat + j*nxnb ;
          valuePointer Bjp = Bmat + jp*nxnb ;
          valuePointer Cj  = Cmat + j*nxnb ;
          valuePointer Cjp = Cmat + jp*nxnb ;

          applyT( nth, T, P, Bjp, n, Bj, n, nb ) ;
          
          // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          if ( selected == BORDERED_QRP ) {
            integer i = n ;
            do {
              --i ;
              if ( P[i] > i )
                swap( nb, Cj+i*nb, 1, Cj+P[i]*nb, 1 ) ;
            } while ( i > 0 ) ;
          }
          trsm( RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
                nb, n, 1.0, T, nx2, Cj, nb ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, n, n,
                -1.0, Cj,  nb,
                      Dj,  n,
                 1.0, Cjp, nb ) ;

          spin1.lock() ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, nb, n,
                -1.0, Cj,   nb,
                      Bj,   n,
                 1.0, Fmat, nb ) ;
          spin1.unlock() ;

        }
      }

      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      barrier.count_down_and_wait() ;
      #endif

      // -----------------------------------------------------------------------
      if ( nb > 0 ) {
        for ( integer j = k ; j < nblock ; j += 2*k ) {

          #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
          CHECK_TASK_DONE(0x02);
          #endif

          integer jpp = std::min(j+k,nblock) ;

          valueConstPointer Ej  = Emat + j*nxn ;
          valueConstPointer Cj  = Cmat + j*nxnb ;
          valuePointer      Cpp = Cmat + jpp*nxnb ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, n, n,
                -1.0, Cj,  nb,
                      Ej,  n,
                 1.0, Cpp, nb ) ;
        }

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        barrier.count_down_and_wait() ;
        #endif
      }
    }
  }
  
  /*\
   |   _                 _     _           _       _     _            _
   |  | | ___   __ _  __| |   | | __ _ ___| |_    | |__ | | ___   ___| | __
   |  | |/ _ \ / _` |/ _` |   | |/ _` / __| __|   | '_ \| |/ _ \ / __| |/ /
   |  | | (_) | (_| | (_| |   | | (_| \__ \ |_    | |_) | | (_) | (__|   <
   |  |_|\___/ \__,_|\__,_|___|_|\__,_|___/\__|___|_.__/|_|\___/ \___|_|\_\
   |                     |_____|             |_____|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::load_last_block() {
    /*
    //    n   n   q  nb
    //  +---+---+---+---+
    //  | D | E | 0 | B | n
    //  +===+===+===+===+
    //  |   |   |   |   |
    //  |H0 |HN |Hq |Hp | n+q
    //  |   |   |   |   |
    //  +---+---+---+---+
    //  | C | C |Cq | F | nb
    //  +---+---+---+---+
    */
    valuePointer Cnb = Cmat + nblock*nxnb ;
    valuePointer W0 = Hmat ;
    valuePointer WN = W0+n*N ;
    valuePointer Wq = WN+n*N ;
    valuePointer Wp = Wq+q*N ;

    gecopy( n,  n,  Dmat, n, W0, N ) ;
    gecopy( n,  n,  Emat, n, WN, N ) ;
    gezero( n,  q,           Wq, N ) ;
    gecopy( n,  nb, Bmat, n, Wp, N ) ;

    gecopy( n+q, N, H0Nqp, n+q, Hmat+n, N ) ;

    W0 += nx2+q ; WN += nx2+q ; Wq += nx2+q ; Wp += nx2+q ;

    gecopy( nb, n,  Cmat,  nb, W0, N ) ;
    gecopy( nb, n,  Cnb,   nb, WN, N ) ;
    gecopy( nb, q,  Cqmat, nb, Wq, N ) ;
    gecopy( nb, nb, Fmat,  nb, Wp, N ) ;

  }

  /*\
   |    __            _             _             _           _
   |   / _| __ _  ___| |_ ___  _ __(_)_______    | | __ _ ___| |_
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \   | |/ _` / __| __|
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/   | | (_| \__ \ |_
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|___|_|\__,_|___/\__|
   |                                        |_____|
  \*/

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_last() {
    integer info = 0 ;
    switch ( last_selected ) {
    case BORDERED_LAST_LU:
      info = getrf( N, N, Hmat, N, Hperm ) ;
      break ;
    case BORDERED_LAST_QR:
      info = geqrf( N, N, Hmat, N, Htau, Work, Lwork ) ;
      break ;
    case BORDERED_LAST_QRP:
      info = geqp3( N, N, Hmat, N, Hperm, Htau, Work, Lwork ) ;
      // convert permutation to exchanges
      if ( info == 0 )
        for ( integer i = 0 ; i < N ; ++i ) {
          integer j = i ;
          while ( j < N ) { if ( Hperm[j] == i+1 ) break ; ++j ; }
          std::swap( Hperm[j], Hperm[i] ) ;
          Hswaps[i] = j ;
        }
      break ;
    }
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize_last INFO = " << info << " N = " << N ) ;
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
      info = getrs( NO_TRANSPOSE, N, 1, Hmat, N, Hperm, X, N ) ;
      break ;
    case BORDERED_LAST_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    N, 1, // righe x colonne
                    N-1,  // numero riflettori usati nel prodotto Q
                    Hmat, N /*ldA*/,
                    Htau,
                    X, N,
                    Work, Lwork ) ;
      if ( info == 0 ) trsv( UPPER, NO_TRANSPOSE, NON_UNIT, N, Hmat, N, X, 1 ) ;
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info ) ;
    if ( last_selected == BORDERED_LAST_QRP )
      for ( integer i = 0 ; i < N ; ++i )
        if ( Hswaps[i] > i ) std::swap( X[i], X[Hswaps[i]] ) ;
    swap( n, X, 1, x, 1 ) ;
  }

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
      info = getrs( NO_TRANSPOSE, N, nrhs, Hmat, N, Hperm, X, ldX ) ;
      break ;
    case BORDERED_LAST_QR:
    case BORDERED_LAST_QRP:
      info = ormqr( LEFT, TRANSPOSE,
                    N, nrhs, // righe x colonne
                    N-1,  // numero riflettori usati nel prodotto Q
                    Hmat, N /*ldA*/,
                    Htau,
                    X, ldX,
                    Work, Lwork ) ;
      if ( info == 0 )
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
              N, nrhs, 1.0, Hmat, N, X, ldX ) ;
      break ;
    }
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last INFO = " << info ) ;
    if ( last_selected == BORDERED_LAST_QRP )
      for ( integer i = 0 ; i < N ; ++i )
        if ( Hswaps[i] > i )
          swap( nrhs, X+i, ldX, X+Hswaps[i], ldX ) ;
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
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    std::fill( task_done, task_done + nblock, 0 ) ;
    barrier.setup(numThread) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt] = std::thread( &BorderedCR<t_Value>::solve_mt, this, nt, x ) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt].join() ;
    #else
    solve_mt( 0, x ) ;
    #endif
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve( integer      nrhs,
                              valuePointer rhs,
                              integer      ldRhs ) const {
    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    std::fill( task_done, task_done + nblock, 0 ) ;
    barrier.setup(numThread) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt] = std::thread( &BorderedCR<t_Value>::solve_n_mt, this, nt, nrhs, rhs, ldRhs ) ;
    for ( integer nt = 0 ; nt < numThread ; ++nt )
      threads[nt].join() ;
    #else
    solve_mt( 0, x ) ;
    #endif
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_mt( integer nth, valuePointer x ) const {
    valuePointer xn = x + (nblock+1)*n + q ;
    integer k = 1 ;
    while ( k < nblock ) {
      for ( integer j = k ; j < nblock ; j += 2*k ) {

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        CHECK_TASK_DONE(0x01);
        #endif

        integer jp = j-k ;
        valuePointer const T = Tmat  + j*Tsize ;
        integer *    const P = Rperm + j*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xjp = x + jp*n ;
        applyT( nth, T, P, xjp, xj ) ;
        if ( nb > 0 ) {
          valuePointer Cj = Cmat + j*nxnb ;
          spin1.lock() ;
          gemv( NO_TRANSPOSE, nb, n, -1.0, Cj, nb, xj, 1, 1.0, xn, 1 ) ;
          spin1.unlock() ;
        }
      }

      k *= 2 ;

      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      barrier.count_down_and_wait() ;
      #endif
    }

    if ( nth == 0 ) solve_last( x ) ;

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    barrier.count_down_and_wait() ;
    #endif

    while ( (k/=2) > 0 ) {
      for ( integer j = k ; j < nblock ; j += 2*k ) {

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        CHECK_TASK_DONE(0x02);
        #endif

        integer      jp  = j-k ;
        integer      jpp = std::min(j+k,nblock) ;
        valuePointer Dj  = Dmat + j*nxn ;
        valuePointer Ej  = Emat + j*nxn ;
        valuePointer xj  = x + j*n ;
        valuePointer xp  = x + jp*n ;
        valuePointer xpp = x + jpp*n ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj, n, xp,  1, 1.0, xj, 1 ) ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej, n, xpp, 1, 1.0, xj, 1 ) ;
        if ( nb > 0 ) {
          valuePointer Bj = Bmat + j*nxnb ;
          gemv( NO_TRANSPOSE, n, nb, -1.0, Bj, n, xn, 1, 1.0, xj, 1 ) ;
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        valuePointer const T = Tmat + j*Tsize ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, nx2, xj, 1 ) ;
        if ( selected == BORDERED_QRP ) {
          integer * const P = Rperm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              std::swap( xj[i], xj[P[i]] ) ;
        }
      }

      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      barrier.count_down_and_wait() ;
      #endif
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_n_mt( integer      nth,
                                   integer      nrhs,
                                   valuePointer x,
                                   integer      ldX ) const {
    valuePointer xn = x + (nblock+1)*n + q ;
    integer k = 1 ;
    while ( k < nblock ) {
      for ( integer j = k ; j < nblock ; j += 2*k ) {

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        CHECK_TASK_DONE(0x01);
        #endif

        integer jp = j-k ;
        valuePointer const T = Tmat  + j*Tsize ;
        integer *    const P = Rperm + j*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xjp = x + jp*n ;
        applyT( nth, T, P, xjp, ldX, xj, ldX, nrhs ) ;
        if ( nb > 0 ) {
          valuePointer Cj = Cmat + j*nxnb ;
          spin1.lock() ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, nrhs, n,
                -1.0, Cj, nb,
                      xj, ldX,
                 1.0, xn, ldX ) ;
          spin1.unlock() ;
        }
      }

      k *= 2 ;

      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      barrier.count_down_and_wait() ;
      #endif
    }

    if ( nth == 0 ) solve_last( nrhs, x, ldX ) ;

    #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
    barrier.count_down_and_wait() ;
    #endif

    while ( (k/=2) > 0 ) {
      for ( integer j = k ; j < nblock ; j += 2*k ) {

        #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
        CHECK_TASK_DONE(0x02);
        #endif

        integer      jp  = j-k ;
        integer      jpp = std::min(j+k,nblock) ;
        valuePointer Dj  = Dmat + j*nxn ;
        valuePointer Ej  = Emat + j*nxn ;
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
        if ( nb > 0 ) {
          valuePointer Bj = Bmat + j*nxnb ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                n, nrhs, nb,
                -1.0, Bj, n,
                      xn, ldX,
                 1.0, xj, ldX ) ;
        }

        // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        valuePointer const T = Tmat + j*Tsize ;
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
              n, nrhs, 1.0, T, nx2, xj, ldX ) ;
        if ( selected == BORDERED_QRP ) {
          integer * const P = Rperm + j*n ;
          for ( integer i = 0 ; i < n ; ++i )
            if ( P[i] > i )
              swap( nrhs, xj+i, ldX, xj+P[i], ldX ) ;
        }
      }

      #ifdef BORDERED_CYCLIC_REDUCTION_USE_THREAD
      barrier.count_down_and_wait() ;
      #endif
    }
  }

  template class BorderedCR<float> ;
  template class BorderedCR<double> ;

}
