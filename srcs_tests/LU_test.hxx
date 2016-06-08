/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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

  #define NSIZE 10

  alglin::integer n    = NSIZE ;
  alglin::integer nblk = 500000 ;
  alglin::integer q    = 40 ;
  alglin::integer nq   = n+q ;
/*
  alglin::integer n    = 2 ;
  alglin::integer nblk = 11 ;
  alglin::integer q    = 0 ;
  alglin::integer nq   = n+q ;
*/
  alglin::integer N   = nblk*n + nq ;
  alglin::integer nnz = nblk*(2*n*n) + nq*(2*n+q) + 5*N ;

  alglin::Malloc<valueType>       baseValue("real") ;
  alglin::Malloc<alglin::integer> baseIndex("integer") ;

  baseValue.allocate(nnz) ;
  baseIndex.allocate(N) ;

  valueType * AdAu  = baseValue(2*n*n*nblk) ;
  valueType * H0    = baseValue(nq*n) ;
  valueType * HN    = baseValue(nq*n) ;
  valueType * Hq    = baseValue(nq*q) ;
  valueType * x     = baseValue(N) ;
  valueType * xref  = baseValue(N) ;
  valueType * xref1 = baseValue(N) ;
  valueType * rhs   = baseValue(N) ;
  valueType * resid = baseValue(N) ;
  
  for ( int i = 0 ; i < nq ; ++i ) {
    for ( int j = 0 ; j < n ; ++j ) {
      H0[i+j*nq] = rand(-1,0) ;
      HN[i+j*nq] = rand(-1,0) ;
    }
  }

  for ( int i = 0 ; i < nq ; ++i )
    for ( int j = 0 ; j < q ; ++j )
      Hq[i+j*nq] = rand(-1,0) ;

  // forzo diagonale dominanza
  valueType diag = 2*n ;
  for ( int k = 0 ; k < nblk ; ++k ) {
    valueType * AdAu_k = AdAu + 2*k*n*n ;
    for ( int i = 0 ; i < n ; ++i )
      for ( int j = 0 ; j < 2*n ; ++j )
        AdAu_k[i+j*n] = rand(-1,0) ;
    for ( int i = 0 ; i < n ; ++i )
      AdAu_k[i*(n+1)] += diag ;
  }

  for ( int j = 0 ; j < n ; ++j ) HN[j*(nq+1)]   += diag ;
  for ( int j = 0 ; j < q ; ++j ) Hq[n+j*(nq+1)] += diag ;

  cout << "N = " << N << '\n' ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  std::copy( x, x+N, xref ) ;
  alglin::babd_mv<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                              1.0, x, 1, 0, rhs, 1 ) ;

  alglin::babd_residue<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                                   rhs, 1, x, 1, resid, 1 ) ;

  cout << "Check residue |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;
