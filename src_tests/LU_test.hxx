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

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

  //#define NSIZE 10
  #define NSIZE 5

  alglin::integer n    = NSIZE;
  alglin::integer nblk = 1000000;
  alglin::integer q    = 40;
  //alglin::integer q    = 3;
  alglin::integer nq   = n+q;
/*
  alglin::integer n    = 2;
  alglin::integer nblk = 11;
  alglin::integer q    = 0;
  alglin::integer nq   = n+q;
*/
  alglin::integer N   = nblk*n + nq;
  alglin::integer nnz = nblk*(2*n*n) + nq*(2*n+q) + 5*N+9*N;

  alglin::Malloc<real_type>       base_value("real");
  alglin::Malloc<alglin::integer> base_index("integer");

  base_value.allocate(nnz);
  baseIndex.allocate(N);

  real_type * AdAu  = base_value(2*n*n*nblk);
  real_type * H0    = base_value(nq*n);
  real_type * HN    = base_value(nq*n);
  real_type * Hq    = base_value(nq*q);
  real_type * x     = base_value(N*10); // extra space per multiple rhs
  real_type * xref  = base_value(N);
  real_type * xref1 = base_value(N);
  real_type * rhs   = base_value(N);
  real_type * resid = base_value(N);

  for ( int i = 0; i < nq; ++i ) {
    for ( int j = 0; j < n; ++j ) {
      H0[i+j*nq] = rand(-1,0);
      HN[i+j*nq] = rand(-1,0);
    }
  }

  for ( int i = 0; i < nq; ++i )
    for ( int j = 0; j < q; ++j )
      Hq[i+j*nq] = rand(-1,0);

  // forzo diagonale dominanza
  real_type diag = 2*n;
  for ( int k = 0; k < nblk; ++k ) {
    real_type * AdAu_k = AdAu + 2*k*n*n;
    for ( int i = 0; i < n; ++i )
      for ( int j = 0; j < 2*n; ++j )
        AdAu_k[i+j*n] = rand(-1,0);
    for ( int i = 0; i < n; ++i )
      AdAu_k[i*(n+1)] += diag;
  }

  for ( int j = 0; j < n; ++j ) HN[j*(nq+1)]   += diag;
  for ( int j = 0; j < q; ++j ) Hq[n+j*(nq+1)] += diag;

  cout << "N = " << N << '\n';

  for ( alglin::integer i = 0; i < N; ++i ) x[i] = i;
  Copy_n( x, N, xref );
  alglin::babd_mv<real_type>( nblk, n, q, AdAu, H0, HN, Hq,
                              1.0, x, 1, 0, rhs, 1 );

  alglin::babd_residue<real_type>( nblk, n, q, AdAu, H0, HN, Hq,
                                   rhs, 1, x, 1, resid, 1 );

  cout << "Check residue |r|_inf = "
       << alglin::absmax( N, resid, 1 ) << '\n';
