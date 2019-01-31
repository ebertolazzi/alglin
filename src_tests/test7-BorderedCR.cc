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

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"
#include "BABD_BorderedCR.hh"

#if defined(__GCC__) || defined(__GNUC__) 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wc99-extensions"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif


using namespace std;
typedef double valueType;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

int
main() {

  alglin::BorderedCR<double> BCR, BCR_SAVED;

  //#define NSIZE 10
  #define NSIZE 6

  alglin::integer n      = NSIZE;
  alglin::integer nblock = 200000;
  alglin::integer qx     = 4;// 4+1;
  alglin::integer qr     = 4;// 4;
  alglin::integer nx     = 1;// 2-1;
  alglin::integer nr     = 1;//2;
  alglin::integer N      = (nblock+1)*n+nx+qx;
 
  BCR.allocate( nblock, n, qr, qx, nr, nx );

  alglin::Malloc<valueType>       baseValue("real");
  alglin::Malloc<alglin::integer> baseIndex("integer");
  
  baseValue.allocate( size_t(7*N) );
  
  valueType diag = 1.01*n;

  valueType * x     = baseValue(size_t(2*N)); // extra space per multiple rhs
  valueType * xref  = baseValue(size_t(N));
  valueType * xref1 = baseValue(size_t(N));
  valueType * rhs   = baseValue(size_t(2*N));
  valueType * resid = baseValue(size_t(N));

  BCR.select_LU();
  //BCR.select_QR();
  //BCR.select_QRP();

  //BCR.select_last_LU();
  //BCR.select_last_LUPQ();
  BCR.select_last_SVD();
  //BCR.select_last_QR();
  //BCR.select_last_QRP();
  //BCR.select_last_LSS();
  //BCR.select_last_LSY();

  for ( int i = 0; i < (n+qr); ++i ) {
    for ( int j = 0; j < (2*n+qx+nx); ++j ) {
      BCR.H(i,j) = rand(-1,0);
    }
    BCR.H(i,i+n) += diag; // force diagonal dominance
  }

  for ( int k = 0; k < nblock; ++k ) {
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < n; ++j ) {
        BCR.D(k,i,j) = rand(-1,0);
        BCR.E(k,i,j) = rand(-1,0);
      }
      BCR.D(k,i,i) += diag; // force diagonal dominance
    }
    for ( int i = 0; i < n; ++i ) {
      for ( int j = 0; j < nx; ++j ) BCR.B(k,i,j) = rand(-0.1,0.1);
      for ( int j = 0; j < nr; ++j ) BCR.C(k,j,i) = rand(-0.1,0.1);
    }
  }
  for ( int i = 0; i < nr; ++i ) {
    for ( int j = 0; j < nx; ++j ) {
      BCR.F(i,j) = rand(-0.1,0.1);
    }
    for ( int j = 0; j < qx; ++j ) {
      BCR.Cq(i,j) = rand(-0.1,0.1);
    }
    BCR.F(i,i) += diag; // force diagonal dominance
  }
  for ( int i = 0; i < n; ++i ) {
    for ( int j = 0; j < nr; ++j ) {
      BCR.C(nblock,j,i) = 1;
    }
  }

  for ( alglin::integer i = 0; i < N; ++i ) x[i] = 1+ (i % 100);
  std::copy( x, x+N, xref );
  BCR.Mv( x, rhs );
  BCR_SAVED.dup( BCR );

  cout << "N      = " << N      << '\n'
       << "nblock = " << nblock << '\n'
       << "n      = " << n      << '\n'
       << "nr     = " << nr     << '\n'
       << "nx     = " << nx     << '\n'
       << "qr     = " << qr     << '\n'
       << "qx     = " << qx     << '\n';
  /*
  ofstream file("mat.txt");
  file.precision(15);
  BCR.dump_ccoord( file );
  file.close();
  file.open("rhs.txt");
  file.precision(15);
  for ( int i = 0; i < N; ++i )
    file << rhs[i] << '\n';
  file.close();
  */

  TicToc tm;
  tm.reset();
  tm.tic();
  BCR.factorize();
  tm.toc();
  cout << "\nFactorize = " << tm.elapsedMilliseconds() << " [ms]\n\n";

  std::copy( rhs, rhs+N, x );
  std::copy( rhs, rhs+N, x+N );
  tm.tic();
  int ns = 1;
  BCR.solve( x );
  #if 1
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  std::copy( rhs, rhs+N, x );
  BCR.solve( x ); ++ns;
  #endif
  tm.toc();
  cout << "\nSolve = " << tm.elapsedMilliseconds()/ns << " [ms]\n\n";

  alglin::copy( N, xref, 1, xref1, 1 );
  alglin::axpy( N, -1.0, x, 1, xref1, 1 );
  valueType err = alglin::absmax( N, xref1, 1 );
  cout << "Check |err|_inf = " << err << '\n';
  ALGLIN_ASSERT( err < 1e-8, "test failed!" );
  err = alglin::asum( N, xref1, 1 )/N;
  cout << "Check |err|_1/N = " << err << '\n';
  ALGLIN_ASSERT( err < 1e-8, "test failed!" );

  std::copy( rhs, rhs+2*N, x );
  tm.tic();
  ns = 1;
  BCR.solve( 2, x, N );
  #if 1
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy( rhs, rhs+2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  #endif
  tm.toc();
  cout << "\nSolve2 = " << tm.elapsedMilliseconds()/ns << " [ms]\n\n";

  /*
  file.open("sol.txt");
  file.precision(15);
  for ( int i = 0; i < N; ++i )
    file << x[i] << '\n';
  file.close();
  */
  //for ( int i = 0; i < N; ++i )
  //  cout << i << " " << x[i] << '\n';

  alglin::copy( N, xref, 1, xref1, 1 );
  alglin::axpy( N, -1.0, x, 1, xref1, 1 );
  err = alglin::absmax( N, xref1, 1 );
  cout << "Check |err|_inf = " << err << '\n';
  ALGLIN_ASSERT( err < 1e-8, "test failed!" );
  err =  alglin::asum( N, xref1, 1 )/N;
  cout << "Check |err|_1/N = " << err << '\n';
  ALGLIN_ASSERT( err < 1e-8, "test failed!" );

  cout << "\n\ncheck residual\n\n";

  std::copy( rhs, rhs+N, resid );
  alglin::scal( N, -1.0, resid, 1 );
  std::copy( rhs, rhs+N, x );
  BCR.solve( x );
  BCR_SAVED.addMv( x, resid );

  valueType res = alglin::nrm2( BCR_SAVED.numRows(), resid, 1 );
  cout << "||res||_2   = " << res << '\n';
  ALGLIN_ASSERT( res < 1e-6, "test failed!" );
  res = alglin::asum( BCR_SAVED.numRows(), resid, 1 );
  cout << "||res||_1   = " << res << '\n';
  ALGLIN_ASSERT( res < 1e-6, "test failed!" );
  res = alglin::absmax( BCR_SAVED.numRows(), resid, 1 );
  cout << "||res||_inf = " << res << '\n';
  ALGLIN_ASSERT( res < 1e-6, "test failed!" );

  cout << "All done!\n";

  return 0;
}
