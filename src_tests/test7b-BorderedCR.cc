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

int
main() {

  try {

  alglin::integer nth = std::thread::hardware_concurrency();

  #ifdef LAPACK_WRAPPER_USE_OPENBLAS
  openblas_set_num_threads(1);
  goto_set_num_threads(1);
  #endif

  Utils::ThreadPool TP(nth);

  alglin::BorderedCR<double> BCR(&TP), BCR_SAVED(&TP);
  //alglin::BorderedCR<double> BCR(nullptr), BCR_SAVED(nullptr);
  //alglin::BorderedCR_eigen3<double> BCR(&TP), BCR_SAVED(&TP);

  #define NSIZE 8
  //#define NSIZE 40

  alglin::integer n      = NSIZE;
  alglin::integer nblock = 10000;
  //salglin::integer nblock = 200;
  alglin::integer qx     = 4;// 4+1;
  alglin::integer qr     = 4;// 4;
  alglin::integer nx     = 1;// 2-1;
  alglin::integer nr     = 1+1;//2;
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
  //BCR.select_SUPERLU();

  //BCR.select_last_LU();
  //BCR.select_last_LUPQ();
  //BCR.select_last_SVD();
  //BCR.select_last_QR();
  //BCR.select_last_QRP();
  //BCR.select_last_LSS();
  //BCR.select_last_LSY();
  BCR.select_last_PINV();

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
  std::copy_n( x, N, xref );
  BCR.Mv( x, rhs );
  BCR_SAVED.dup( BCR );

  cout << "nthread (avilable)= " << std::thread::hardware_concurrency() << '\n';
  cout << "nthread (used)    = " << nth << '\n';

  cout << "N      = " << N      << '\n'
       << "nblock = " << nblock << '\n'
       << "n      = " << n      << '\n'
       << "nr     = " << nr     << '\n'
       << "nx     = " << nx     << '\n'
       << "qr     = " << qr     << '\n'
       << "qx     = " << qx     << '\n';

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

  fmt::print("All done!\n");

  }
  catch ( std::exception const & err ) {
    std::cerr << "Error: " << err.what();
  }
  catch (...) {
    std::cerr << "Unknwn error\n";
  }

  return 0;
}
