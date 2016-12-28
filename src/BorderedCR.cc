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
/*
  string
  LastBlock_to_string( LASTBLOCK_Choice c ) {
    switch ( c ) {
      case LASTBLOCK_LU:  return "last block LU"  ;
      case LASTBLOCK_QR:  return "last block QR"  ;
      case LASTBLOCK_QRP: return "last block QRP" ;
      case LASTBLOCK_SVD: return "last block SVD" ;
    }
    return "last block not selected" ;
  }
*/
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
    valueConstPointer Hp, integer ldP,
    valueConstPointer Hq, integer ldQ
  ) {
    integer m = n + q ;
    gecopy( m, n,  H0, ld0, H0Npq,            m ) ;
    gecopy( m, n,  HN, ldN, H0Npq+m*n,        m ) ;
    gecopy( m, nb, Hp, ldP, H0Npq+m*nx2,      m ) ;
    gecopy( m, q,  Hq, ldQ, H0Npq+m*(nx2+nb), m ) ;
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
    t_Value const * xe = x+(nblock+1)*n;
    t_Value *       yy = res ;
    for ( integer i = 0 ; i < nblock ; ++i ) {
      gemv( NO_TRANSPOSE, n, n,  1.0, D, n, xx, 1, 0.0, yy, 1 ) ;
      xx += n ;
      gemv( NO_TRANSPOSE, n, n,  1.0, E, n, xx, 1, 1.0, yy, 1 ) ;
      if ( nb > 0 ) gemv( NO_TRANSPOSE, n, nb, 1.0, B, n, xe, 1, 1.0, yy, 1 ) ;
      yy += n ; D += nxn ; E += nxn ; B += nxnb ;
    }
    if ( nb > 0 ) {
      gemv( NO_TRANSPOSE, nb, nb, 1.0, Fmat, nb, xe, 1, 0.0, yy, 1 ) ;
      t_Value const * C = Cmat ;
      xx = x ;
      for ( integer i = 0 ; i <= nblock ; ++i ) {
        gemv( NO_TRANSPOSE, nb, n, 1.0, C, nb, xx, 1, 1.0, yy, 1 ) ;
        xx += n ; C += nxnb ;
      }
    }
    yy += nb ;

    integer m = n+q ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H0Npq,            m, x,     1, 0.0, yy, 1 ) ;
    gemv( NO_TRANSPOSE, m, n,  1.0, H0Npq+m*n,        m, xe-n,  1, 1.0, yy, 1 ) ;
    gemv( NO_TRANSPOSE, m, nb, 1.0, H0Npq+m*nx2,      m, xe,    1, 1.0, yy, 1 ) ;
    gemv( NO_TRANSPOSE, m, q,  1.0, H0Npq+m*(nx2+nb), m, xe+nb, 1, 1.0, yy, 1 ) ;
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
    integer nnz = 2*nblock*(nxn+nxnb) + nxnb + nb*nb + (n+q)*(2*n+nb+q) ;
    stream << nnz << '\n' ;

    // bidiagonal
    integer je = (nblock+1)*n ;
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
            << ii+i << '\t' << je+j << '\t' << B[i+j*n] << '\n' ;
    }
    integer ie = nblock*n ;
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
          << ie+i << '\t' << je+j << '\t' << Fmat[i+j*nb] << '\n' ;

    ie += nb ;
    je -= n ;
    valuePointer H0 = H0Npq ;
    for ( integer i = 0 ; i < n+q ; ++i )
      for ( integer j = 0 ; j < n ; ++j )
        stream << ie+i << '\t' << j << '\t' << H0[i+j*(n+q)] << '\n' ;

    valuePointer HNpq = H0Npq+n*(n+q) ;
    for ( integer i = 0 ; i < n+q ; ++i )
      for ( integer j = 0 ; j < n+q+nb ; ++j )
        stream << ie+i << '\t' << je+j << '\t' << HNpq[i+j*(n+q)] << '\n' ;

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
                               valuePointer      BOTTOM,
                               integer           ncol ) const {
    gecopy( n, ncol, TOP,    n, Work,   nx2 ) ;
    gecopy( n, ncol, BOTTOM, n, Work+n, nx2 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = swaps( ncol, Work, nx2, 0, n-1, iperm, 1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize swaps INFO = " << info ) ;
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, ncol, 1.0, T, nx2, Work, nx2 ) ;
    gecopy( n, ncol, Work+n, nx2, TOP, n ) ;
    gemm( NO_TRANSPOSE, NO_TRANSPOSE,
          n, ncol, n,
          -1.0, T+n,  nx2,
                Work, nx2,
           1.0, TOP,  n ) ;
    gecopy( n, ncol, Work, nx2, BOTTOM, n ) ;
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, n, ncol, 1.0, T, nx2, BOTTOM, n ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::applyT( valueConstPointer T,
                               integer const *   iperm,
                               valuePointer      TOP,
                               valuePointer      BOTTOM ) const {
    copy( n, TOP,    1, Work,   1 ) ;
    copy( n, BOTTOM, 1, Work+n, 1 ) ;
    // Apply row interchanges to the right hand sides.
    integer info = swaps( 1, Work, nx2, 0, n-1, iperm, 1 ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::factorize swaps INFO = " << info ) ;
    trsv( LOWER, NO_TRANSPOSE, UNIT, n, T, nx2, Work, 1 ) ;
    copy( n, Work+n, 1, TOP, 1 ) ;
    gemv( NO_TRANSPOSE,
          n, n,
          -1.0, T+n,  nx2,
                Work, 1,
           1.0, TOP,  1 ) ;
    copy( n, Work, 1, BOTTOM, 1 ) ;
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, T, nx2, BOTTOM, 1 ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::factorize() {
    for ( integer k = 1 ; k < nblock ; k = 2*k ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j = j1-k ;
        valuePointer T   = Tmat  + j1*2*nxn ;
        integer *    P   = Tperm + j1*n ;
        valuePointer Dj  = Dmat  + j*nxn ;
        valuePointer Ej  = Emat  + j*nxn ;
        valuePointer Dj1 = Dmat  + j1*nxn ;
        valuePointer Ej1 = Emat  + j1*nxn ;

        buildT( Ej, Dj1, T, P ) ;
        zero( nxn, Dj1, 1 ) ;
        applyT( T, P, Dj, Dj1, n ) ;
        zero( nxn, Ej, 1 ) ;
        applyT( T, P, Ej, Ej1, n  ) ;

        if ( nb > 0 ) {
          integer j2 = j1+k ;
          if ( j2 >= nblock ) j2 = nblock ;

          valuePointer Bj  = Bmat + j*nxnb ;
          valuePointer Cj  = Cmat + j*nxnb ;
          valuePointer Bj1 = Bmat + j1*nxnb ;
          valuePointer Cj1 = Cmat + j1*nxnb ;
          valuePointer Cj2 = Cmat + j2*nxnb ;

          applyT( T, P, Bj, Bj1, nb ) ;

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

    /*
    //  +-----+-----+---+--+
    //  |  D  |  E  | B |  | n
    //  +-----+-----+---+--+
    //  |  C  |  C  | F |  | nb
    //  +=====+=====+===+==+
    //  |     |     |   |  |
    //  |  H0 | HN  |Hp |Hq| n+q
    //  |     |     |   |  |
    //  +-----+-----+---+--+
    */
    integer N = nx2+nb+q ;
    valuePointer Cnb = Cmat + nblock*nxnb ;
    valuePointer W0 = Hmat ;
    valuePointer W1 = W0+n*N ;
    valuePointer W2 = W1+n*N ;

    gezero( n+nb, q, W2+nb*N, N ) ;
    gecopy( n,  n,  Dmat, n, W0, N ) ;
    gecopy( n,  n,  Emat, n, W1, N ) ;
    gecopy( n,  nb, Bmat, n, W2, N ) ;
    gecopy( nb, n,  Cmat, nb, W0+n, N ) ;
    gecopy( nb, n,  Cnb,  nb, W1+n, N ) ;
    gecopy( nb, nb, Fmat, nb, W2+n, N ) ;
    gecopy( n+q, N, H0Npq, n+q, W0+n+nb, N ) ;
    integer info = getrf( N, N, Hmat, N, Hperm ) ;
    ALGLIN_ASSERT( info == 0,
                   "BorderedCR::factorize getrf INFO = " << info <<
                   " N = " << N ) ;
  }

  template <typename t_Value>
  void
  BorderedCR<t_Value>::solve( valuePointer x ) const {
    valuePointer xn = x + (nblock+1)*n ;
    integer k = 1 ;
    while ( k < nblock ) {
      for ( integer j1 = k ; j1 < nblock ; j1 += 2*k ) {
        integer j = j1-k ;
        valuePointer T   = Tmat + j1*2*nxn ;
        integer *    P   = Tperm + j1*n ;
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

    integer N = nx2+nb+q ;
    valuePointer X = Work ;
    copy( n,      x,    1, X,   1 ) ;
    copy( n+nb+q, xn-n, 1, X+n, 1 ) ;
    integer info = getrs( NO_TRANSPOSE, N, 1, Hmat, N, Hperm, X, N ) ;
    ALGLIN_ASSERT( info == 0, "BorderedCR::solve getrs INFO = " << info ) ;
    copy( n,      X,   1, x,    1 ) ;
    copy( n+nb+q, X+n, 1, xn-n, 1 ) ;
    
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

  template class BorderedCR<float> ;
  template class BorderedCR<double> ;

}
