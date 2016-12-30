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

    Lwork = max(N,2*n*max(n,nb)) ;

    valueType tmp ; // get optimal allocation
    integer info = geqrf( N, N, nullptr, N, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    info = geqp3( N, N, nullptr, N, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::allocate call alglin::geqrf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;

    integer nnz = Lwork+N*(N+1)+nb*(nb+q)+nxnb+nblock*(2*nxnb+4*nxn+n)+(n+q)*(nx2+nb+q) ;
    baseValue.allocate(size_t(nnz)) ;
    baseInteger.allocate(size_t(2*(nblock*n+N))) ;

    Bmat  = baseValue(size_t(nblock*nxnb)) ;
    Cmat  = baseValue(size_t((nblock+1)*nxnb)) ;
    Cqmat = baseValue(size_t(nb*q)) ;
    Dmat  = baseValue(size_t(nblock*nxn)) ;
    Emat  = baseValue(size_t(nblock*nxn)) ;
    Fmat  = baseValue(size_t(nb*nb)) ;
    
    Tmat  = baseValue(size_t(2*nblock*nxn)) ;
    Ttau  = baseValue(size_t(nblock*n)) ;
    Hmat  = baseValue(size_t(N*N)) ;
    Htau  = baseValue(size_t(N)) ;
    H0Nqp = baseValue(size_t((n+q)*N)) ;

    Work  = baseValue(size_t(Lwork)) ;

    Hperm  = baseInteger(size_t(N)) ;
    Hswaps = baseInteger(size_t(N)) ;
    Rperm  = baseInteger(size_t(nblock*n)) ;
    Cperm  = baseInteger(size_t(nblock*n)) ;
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
   |
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
  /*\
    Compute   / L \ (U) = / LU \ = / TOP    \
              \ M /       \ MU /   \ BOTTOM /
  \*/
  template <typename t_Value>
  void
  BorderedCR<t_Value>::buildT( valueConstPointer TOP,
                               valueConstPointer BOTTOM,
                               valuePointer      T,
                               integer *         iperm ) const {
    gecopy( n, n, TOP,    n, T,   nx2 ) ;
    gecopy( n, n, BOTTOM, n, T+n, nx2 ) ;
    integer info = getrf( nx2, n, T, nx2, iperm ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize getrf INFO = " << info ) ;
  }

  /*\
   APPLY
      / I         \ / -M  I \ / L^(-1)   \ / TOP    \
      \    U^(-1) / \  I  0 / \        I / \ BOTTOM /
      = / I         \ / -M  I \ / L^(-1) TOP \
        \    U^(-1) / \  I  0 / \  BOTTOM    /
      = / BOTTOM - M L^(-1) TOP \
        \ U^(-1) L^(-1) TOP /

   ON
      / I         \ / -M  I \ / L^(-1)   \ / L \ (U)
      \    U^(-1) / \  I  0 / \        I / \ M /
        = / I         \ / -M  I \ / I \ (U)
          \    U^(-1) / \  I  0 / \ M /
        = / I         \ / 0 \ (U)   / 0 \
          \    U^(-1) / \ I /     = \ I /
  \*/
  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT( valueConstPointer T,
                               integer const *   iperm,
                               valuePointer      TOP,
                               integer           ldTOP,
                               valuePointer      BOTTOM,
                               integer           ldBOTTOM,
                               integer           ncol ) const {
    valuePointer W = Tmat ;
    gecopy( n, ncol, TOP,    ldTOP,    W,   nx2 ) ;
    gecopy( n, ncol, BOTTOM, ldBOTTOM, W+n, nx2 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = swaps( ncol, W, nx2, 0, n-1, iperm, 1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize swaps INFO = " << info ) ;
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, nx2, W, nx2 ) ;
    gecopy( n, ncol, W+n, nx2, TOP, ldTOP ) ;
    gemm( NO_TRANSPOSE, NO_TRANSPOSE,
          n, ncol, n,
          -1.0, T+n, nx2,      // M
                W,   nx2,      // L^(-1) TOP
           1.0, TOP, ldTOP ) ; // TOP = BOTTOM - M L^(-1) TOP
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, n, ncol, 1.0, T, nx2, W, nx2 ) ;
    gecopy( n, ncol, W, nx2, BOTTOM, ldBOTTOM ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT( valueConstPointer T,
                               integer const *   iperm,
                               valuePointer      TOP,
                               valuePointer      BOTTOM ) const {
    valuePointer W = Tmat ;
    copy( n, TOP,    1, W,   1 ) ;
    copy( n, BOTTOM, 1, W+n, 1 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = swaps( 1, W, nx2, 0, n-1, iperm, 1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize swaps INFO = " << info ) ;
    trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, nx2, W, 1 ) ;
    copy( n, W+n, 1, TOP, 1 ) ;
    gemv( NO_TRANSPOSE,
          n, n,
          -1.0, T+n, nx2,
                W,   1,
           1.0, TOP, 1 ) ;
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, nx2, W, 1 ) ;
    copy( n, W, 1, BOTTOM, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_LU() {
    for ( integer k = 1 ; k < nblock ; k = 2*k ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j = j1-k ;
        valuePointer T   = Tmat  + j1*2*nxn ;
        integer *    P   = Rperm + j1*n ;
        valuePointer Dj  = Dmat  + j*nxn ;
        valuePointer Dj1 = Dmat  + j1*nxn ;
        valuePointer Ej  = Emat  + j*nxn ;
        valuePointer Ej1 = Emat  + j1*nxn ;

        buildT( Ej, Dj1, T, P ) ;

        zero( nxn, Dj1, 1 ) ;
        applyT( T, P, Dj, n, Dj1, n, n ) ;

        zero( nxn, Ej, 1 ) ;
        applyT( T, P, Ej, n, Ej1, n, n ) ;

        if ( nb > 0 ) {
          integer j2 = j1+k ;
          if ( j2 >= nblock ) j2 = nblock ;

          valuePointer Bj  = Bmat + j*nxnb ;
          valuePointer Bj1 = Bmat + j1*nxnb ;
          valuePointer Cj  = Cmat + j*nxnb ;
          valuePointer Cj1 = Cmat + j1*nxnb ;
          valuePointer Cj2 = Cmat + j2*nxnb ;

          applyT( T, P, Bj, n, Bj1, n, nb ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, n, n,
                -1.0, Cj1, nb,
                      Dj1, n,
                 1.0, Cj,  nb ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, n, n,
                -1.0, Cj1, nb,
                      Ej1, n,
                 1.0, Cj2, nb ) ;

          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, nb, n,
                -1.0, Cj1, nb,
                      Bj1, n,
                 1.0, Fmat, nb ) ;
        }
      }
    }
  }


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

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_last_LU() {
    integer info = getrf( N, N, Hmat, N, Hperm ) ;
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize getrf INFO = " << info <<
                   " N = " << N ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_LU( valuePointer x ) const {
    valuePointer X = x + (nblock-1)*n ;
    swap( n, X, 1, x, 1 ) ; // uso x stesso come temporaneo
    integer info = getrs( NO_TRANSPOSE, N, 1, Hmat, N, Hperm, X, N ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_LU getrs INFO = " << info ) ;
    swap( n, X, 1, x, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_LU( integer      nrhs,
                                      valuePointer x,
                                      integer      ldX ) const {
    valuePointer X = x + (nblock-1)*n ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
    integer info = getrs( NO_TRANSPOSE, N, nrhs, Hmat, N, Hperm, X, ldX ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_LU getrs INFO = " << info ) ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_last_QR() {
    integer info = geqrf( N, N, Hmat, N, Htau, Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize geqrf INFO = " << info <<
                   " N = " << N ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_QR( valuePointer x ) const {
    valuePointer X = x + (nblock-1)*n ;
    swap( n, X, 1, x, 1 ) ; // uso x stesso come temporaneo
    integer info = ormqr( LEFT, TRANSPOSE,
                          N, 1, // righe x colonne
                          N-1,  // numero riflettori usati nel prodotto Q
                          Hmat, N /*ldA*/,
                          Htau,
                          X, N,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_QR ormqr INFO = " << info ) ;
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, N, Hmat, N, X, 1 ) ;
    swap( n, X, 1, x, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_QR( integer      nrhs,
                                      valuePointer x,
                                      integer      ldX ) const {
    valuePointer X = x + (nblock-1)*n ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
    integer info = ormqr( LEFT, TRANSPOSE,
                          N, nrhs, // righe x colonne
                          N-1,  // numero riflettori usati nel prodotto Q
                          Hmat, N /*ldA*/,
                          Htau,
                          X, ldX,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_QR ormqr INFO = " << info ) ;
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, N, nrhs, 1.0, Hmat, N, X, ldX ) ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize_last_QRP() {
    integer info = geqp3( N, N, Hmat, N, Hperm, Htau, Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize geqrf INFO = " << info <<
                   " N = " << N ) ;
    // convert permutation to exchanges
    for ( integer i = 0 ; i < N ; ++i ) {
      integer j = i ;
      while ( j < N ) { if ( Hperm[j] == i+1 ) break ; ++j ; }
      std::swap( Hperm[j], Hperm[i] ) ;
      Hswaps[i] = j ;
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_QRP( valuePointer x ) const {
    valuePointer X = x + (nblock-1)*n ;
    swap( n, X, 1, x, 1 ) ; // uso x stesso come temporaneo
    integer info = ormqr( LEFT, TRANSPOSE,
                          N, 1, // righe x colonne
                          N-1,  // numero riflettori usati nel prodotto Q
                          Hmat, N /*ldA*/,
                          Htau,
                          X, N,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_QR ormqr INFO = " << info ) ;
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, N, Hmat, N, X, 1 ) ;
    for ( integer i = 0 ; i < N ; ++i ) {
      if ( Hswaps[i] > i ) std::swap( X[i], X[Hswaps[i]] ) ;
    }
    swap( n, X, 1, x, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_last_QRP( integer      nrhs,
                                       valuePointer x,
                                       integer      ldX ) const {
    valuePointer X = x + (nblock-1)*n ;
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
    integer info = ormqr( LEFT, TRANSPOSE,
                          N, nrhs, // righe x colonne
                          N-1,  // numero riflettori usati nel prodotto Q
                          Hmat, N /*ldA*/,
                          Htau,
                          X, ldX,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve_last_QR ormqr INFO = " << info ) ;
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, N, nrhs, 1.0, Hmat, N, X, ldX ) ;
    for ( integer i = 0 ; i < N ; ++i ) {
      if ( Hswaps[i] > i )
        swap( nrhs, X+i, ldX, X+Hswaps[i], ldX ) ;
    }
    for ( integer i = 0 ; i < nrhs ; ++i ) swap( n, X+i*ldX, 1, x+i*ldX, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_LU( valuePointer x ) const {
    valuePointer xn = x + (nblock+1)*n + q ;
    integer k = 1 ;
    while ( k < nblock ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j = j1-k ;
        valuePointer const T = Tmat + j1*2*nxn ;
        integer *    const P = Rperm + j1*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xj1 = x + j1*n ;
        applyT( T, P, xj, xj1 ) ;
        if ( nb > 0 ) {
          valuePointer Cj1 = Cmat + j1*nxnb ;
          gemv( NO_TRANSPOSE, nb, n, -1.0, Cj1, nb, xj1, 1, 1.0, xn, 1 ) ;
        }
      }
      k *= 2 ;
    }

    switch ( last_selected ) {
      case BORDERED_LAST_LU:  solve_last_LU( x )  ; break ;
      case BORDERED_LAST_QR:  solve_last_QR( x )  ; break ;
      case BORDERED_LAST_QRP: solve_last_QRP( x ) ; break ;
    }

    while ( (k/=2) > 0 ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j  = j1-k ;
        integer j2 = j1+k ;
        if ( j2 >= nblock ) j2 = nblock ;
        valuePointer Dj1 = Dmat + j1*nxn ;
        valuePointer Ej1 = Emat + j1*nxn ;
        valuePointer xj  = x + j*n ;
        valuePointer xj1 = x + j1*n ;
        valuePointer xj2 = x + j2*n ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Dj1, n, xj,  1, 1.0, xj1, 1 ) ;
        gemv( NO_TRANSPOSE, n, n, -1.0, Ej1, n, xj2, 1, 1.0, xj1, 1 ) ;
        if ( nb > 0 ) {
          valuePointer Bj1 = Bmat + j1*nxnb ;
          gemv( NO_TRANSPOSE, n, nb, -1.0, Bj1, n, xn, 1, 1.0, xj1, 1 ) ;
        }
      }
    }
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve_LU( integer      nrhs,
                                 valuePointer x,
                                 integer      ldX ) const {
    valuePointer xn = x + nblock*n + q ;
    integer k = 1 ;
    while ( k < nblock ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j = j1-k ;
        valuePointer const T = Tmat + j1*2*nxn ;
        integer *    const P = Rperm + j1*n ;
        valuePointer xj  = x + j*n ;
        valuePointer xj1 = x + j1*n ;
        applyT( T, P, xj, ldX, xj1, ldX, nrhs ) ;
        if ( nb > 0 ) {
          valuePointer Cj1 = Cmat + j1*nxnb ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                nb, nrhs, n, 
                -1.0, Cj1, nb,
                      xj1, ldX,
                 1.0, xn,  ldX ) ;
        }
      }
      k *= 2 ;
    }

    switch ( last_selected ) {
      case BORDERED_LAST_LU:  solve_last_LU ( nrhs, x, ldX ) ; break ;
      case BORDERED_LAST_QR:  solve_last_QR ( nrhs, x, ldX ) ; break ;
      case BORDERED_LAST_QRP: solve_last_QRP( nrhs, x, ldX ) ; break ;
    }

    while ( (k/=2) > 0 ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j  = j1-k ;
        integer j2 = j1+k ;
        if ( j2 >= nblock ) j2 = nblock ;
        valuePointer Dj1 = Dmat + j1*nxn ;
        valuePointer Ej1 = Emat + j1*nxn ;
        valuePointer xj  = x + j*n ;
        valuePointer xj1 = x + j1*n ;
        valuePointer xj2 = x + j2*n ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Dj1, n,
                    xj,  ldX,
               1.0, xj1, ldX ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Ej1, n,
                    xj2, ldX,
               1.0, xj1, ldX ) ;
        if ( nb > 0 ) {
          valuePointer Bj1 = Bmat + j1*nxnb ;
          gemm( NO_TRANSPOSE, NO_TRANSPOSE,
                n, nrhs, nb,
                -1.0, Bj1, n,
                      xn,  ldX,
                 1.0, xj1, ldX ) ;
        }
      }
    }
  }

  template class BorderedCR<float> ;
  template class BorderedCR<double> ;

}
