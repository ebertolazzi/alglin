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
#include <random>

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

static
void
fill_matrix(
  alglin::BorderedCR<double> & BCR,
  alglin::integer nblock,
  alglin::integer n,
  alglin::integer nr,
  alglin::integer nx,
  alglin::integer qr,
  alglin::integer qx
) {

  BCR.allocate( nblock, n, qr, qx, nr, nx );

  valueType diag = 1.01*n;

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
}

int
main() {

  alglin::integer nth = std::thread::hardware_concurrency();

  #ifdef LAPACK_WRAPPER_USE_OPENBLAS
  openblas_set_num_threads(1);
  goto_set_num_threads(1);
  #endif

  Utils::ThreadPool TP(nth);

  alglin::BorderedCR<double> BCR(&TP), BCR_SAVED(&TP);

  #define NSIZE 8

  alglin::integer n      = NSIZE;
  alglin::integer nblock = 10000;
  alglin::integer qx     = 4;// 4+1;
  alglin::integer qr     = 4;// 4;
  alglin::integer nx     = 1;// 2-1;
  alglin::integer nr     = 1;//2;
  alglin::integer N      = (nblock+1)*n+nx+qx;

  fill_matrix( BCR, nblock, n, nr, nx, qr, qx );

  alglin::Malloc<valueType>       baseValue("real");
  alglin::Malloc<alglin::integer> baseIndex("integer");

  baseValue.allocate( size_t(7*N) );
  valueType * x     = baseValue(size_t(2*N)); // extra space per multiple rhs
  valueType * xref  = baseValue(size_t(N));
  valueType * xref1 = baseValue(size_t(N));
  valueType * rhs   = baseValue(size_t(2*N));
  valueType * resid = baseValue(size_t(N));

  BCR.select_LU();
  //BCR.select_QR();
  //BCR.select_QRP();
  //BCR.select_SUPERLU();

  BCR.select_last_LU();
  //BCR.select_last_LUPQ();
  //BCR.select_last_SVD();
  //BCR.select_last_QR();
  //BCR.select_last_QRP();
  //BCR.select_last_LSS();
  //BCR.select_last_LSY();
  //BCR.select_last_PINV();

  for ( alglin::integer i = 0; i < N; ++i ) x[i] = 1+ (i % 100);
  std::copy_n( x, N, xref );
  BCR.Mv( x, rhs );
  BCR_SAVED.dup( BCR );

  fmt::print(
    "nthread (avilable) = {}\n"
    "nthread (used)     = {}\n",
    std::thread::hardware_concurrency(), nth
  );

  fmt::print(
    "N      = {}\n"
    "nblock = {}\n"
    "n      = {}\n"
    "nr     = {}\n"
    "nx     = {}\n"
    "qr     = {}\n"
    "qx     = {}\n",
    N, nblock, n, nr, nx, qr, qx
  );

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

  Utils::TicToc tm;
  tm.tic();
  BCR.factorize();
  tm.toc();
  fmt::print("\nFactorize = {:.5} [ms]\n\n", tm.elapsed_ms());

  std::copy_n( rhs, N, x );
  std::copy_n( rhs, N, x+N );
  tm.tic();
  int ns = 1;
  BCR.solve( x );
  #if 1
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  std::copy_n( rhs, N, x );
  BCR.solve( x ); ++ns;
  #endif
  tm.toc();
  fmt::print("\nSolve = {:.5} [ms]\n\n", tm.elapsed_ms()/ns);

  alglin::copy( N, xref, 1, xref1, 1 );
  alglin::axpy( N, -1.0, x, 1, xref1, 1 );
  valueType err = alglin::absmax( N, xref1, 1 );
  fmt::print("Check |err|_inf = {:.5}\n",err);
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );
  err = alglin::asum( N, xref1, 1 )/N;
  fmt::print("Check |err|_1/N = {:.5}\n",err);
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

  std::copy_n( rhs, 2*N, x );
  tm.tic();
  ns = 1;
  BCR.solve( 2, x, N );
  #if 1
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  std::copy_n( rhs, 2*N, x );
  BCR.solve( 2, x, N ); ++ns;
  #endif
  tm.toc();
  fmt::print("\nSolve2 = {:.5} [ms]\n",tm.elapsed_ms()/ns);

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
  fmt::print("Check |err|_inf = {:.5}\n",err);
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );
  err = alglin::asum( N, xref1, 1 )/N;
  fmt::print("Check |err|_1/N = {:.5}\n",err);
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

  fmt::print("\n\ncheck residual\n\n");

  std::copy_n( rhs, N, resid );
  alglin::scal( N, -1.0, resid, 1 );
  std::copy_n( rhs, N, x );
  BCR.solve( x );
  BCR_SAVED.addMv( x, resid );

  valueType res = alglin::nrm2( BCR_SAVED.numRows(), resid, 1 );
  fmt::print("||res||_2   = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
  res = alglin::asum( BCR_SAVED.numRows(), resid, 1 );
  fmt::print("||res||_1   = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
  res = alglin::absmax( BCR_SAVED.numRows(), resid, 1 );
  fmt::print("||res||_inf = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

  fmt::print("All done!\n");

  return 0;
}
