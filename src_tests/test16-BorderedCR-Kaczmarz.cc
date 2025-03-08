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
#include "Alglin_Eigen.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

#ifdef LAPACK_WRAPPER_USE_OPENBLAS
#include <omp.h>
#endif

#include <random>

using namespace std;
typedef double real_type;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
real_type
rand( real_type xmin, real_type xmax ) {
  real_type const random{ static_cast<real_type>(generator())/generator.max() };
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
  BCR.fill_zero();

  real_type const diag{ 1.01*n };

  for ( int i{0}; i < (n+qr); ++i ) {
    for ( int j{0}; j < (2*n+qx+nx); ++j ) {
      BCR.H(i,j) = rand(-1,0);
    }
    BCR.H(i,i+n) += diag; // force diagonal dominance
  }

  for ( int k{0}; k < nblock; ++k ) {
    for ( int i{0}; i < n; ++i ) {
      for ( int j{0}; j < n; ++j ) {
        BCR.D(k,i,j) = rand(-1,0);
        BCR.E(k,i,j) = rand(-1,0);
      }
      BCR.D(k,i,i) += diag; // force diagonal dominance
    }
    for ( int i{0}; i < n; ++i ) {
      for ( int j{0}; j < nx; ++j ) BCR.B(k,i,j) = rand(-0.1,0.1);
      for ( int j{0}; j < nr; ++j ) BCR.C(k,j,i) = rand(-0.1,0.1);
    }
  }
  for ( int i{0}; i < nr; ++i ) {
    for ( int j{0}; j < nx; ++j ) {
      BCR.F(i,j) = rand(-0.1,0.1);
    }
    for ( int j{0}; j < qx; ++j ) {
      BCR.Cq(i,j) = rand(-0.1,0.1);
    }
    BCR.F(i,i) += diag; // force diagonal dominance
  }
  for ( int i{0}; i < n; ++i ) {
    for ( int j{0}; j < nr; ++j ) {
      BCR.C(nblock,j,i) = 1;
    }
  }
}

int
main( int argc, char *argv[] ) {

  try {

  alglin::integer nth{ static_cast<alglin::integer>(std::thread::hardware_concurrency()) };

  for ( int i{0}; i < argc; ++i )
    fmt::print( "arg[{}] = {}\n", i, argv[i] );

  if ( argc == 2 ) nth = atoi( argv[1] );

  Utils::ThreadPool1         TP(nth);
  alglin::BorderedCR<double> BCR(&TP);
  alglin::BorderedCR<double> BCR_SAVED(&TP);

  #define NSIZE 8

  alglin::integer n      { NSIZE };
  alglin::integer nblock { 10000 };
  alglin::integer qx     { 4 };// 4+1;
  alglin::integer qr     { 4 };// 4;
  alglin::integer nx     { 1 };// 2-1;
  alglin::integer nr     { 1 };//2;
  alglin::integer N      { (nblock+1) * n + nx + qx };

  fill_matrix( BCR, nblock, n, nr, nx, qr, qx );

  alglin::Malloc<real_type>       base_value("real");
  alglin::Malloc<alglin::integer> base_index("integer");

  base_value.allocate( 7*N );
  real_type * x     { base_value( 2*N) }; // extra space per multiple rhs
  real_type * xref  { base_value( N)   };
  real_type * xref1 { base_value( N)   };
  real_type * rhs   { base_value( 2*N) };
  real_type * resid { base_value( N)   };

  BCR.select_LU();
  //BCR.select_QR();
  //BCR.select_QRP();

  BCR.select_last_LU();
  //BCR.select_last_LUPQ();
  //BCR.select_last_SVD();
  //BCR.select_last_QR();
  //BCR.select_last_QRP();
  //BCR.select_last_LSS();
  //BCR.select_last_LSY();
  //BCR.select_last_PINV();

  for ( alglin::integer i{0}; i < N; ++i ) x[i] = 1+ (i % 100);
  alglin::Copy_n( x, N, xref );
  BCR.Mv( x, rhs );
  BCR_SAVED.dup( BCR );

  fmt::print(
    "nthread (available) = {}\n"
    "nthread (used)      = {}\n",
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

  BCR.set_factorize_use_thread( false );
  BCR.set_solve_use_thread( false );

  fmt::print( "{}\n", BCR.info() );

  Utils::TicToc tm;

  tm.tic();
  bool ok = BCR.factorize();
  tm.toc();
  if ( !ok ) fmt::print( "BCR.factorize failed, last error: {}\n",BCR.last_error() );
  fmt::print( "\nFactorize = {:.5} [ms]\n\n", tm.elapsed_ms() );

  alglin::Copy_n( rhs, N, x );
  alglin::Copy_n( rhs, N, x+N );
  tm.tic();
  int nrhs = 2;
  ok = BCR.solve( x );
  if ( !ok ) {
    fmt::print( "BCR.solve( x ) failed, last error: {}\n",BCR.last_error() );
  }
  tm.toc();
  fmt::print("\nSolve = {:.5} [ms]\n\n", tm.elapsed_ms());

  alglin::Copy_n( xref, N, xref1 );
  alglin::axpy( N, -1.0, x, 1, xref1, 1 );
  real_type err = alglin::absmax( N, xref1, 1 );
  msg.semaphore(
    err > 1e-8 ? 0 : 1,
    fmt::format("Check ‖err‖∞ = {:.5}\n",err)
  );
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );
  err = alglin::asum( N, xref1, 1 )/N;
  msg.semaphore(
    err > 1e-8 ? 0 : 1,
    fmt::format("Check ‖err‖_1/N = {:.5}\n",err)
  );
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

  alglin::Copy_n( rhs, 2*N, x );
  tm.tic();
  ok = BCR.solve( 2, x, N );
  tm.toc();
  if ( !ok ) fmt::print( "BCR.solve( nrhs, x, N ) failed, last error: {}\n",BCR.last_error() );
  fmt::print("\nSolve2 = {:.5} [ms]\n",tm.elapsed_ms());

  alglin::Copy_n( xref, N, xref1 );
  alglin::axpy( N, -1.0, x, 1, xref1, 1 );
  err = alglin::absmax( N, xref1, 1 );
  msg.semaphore(
    err > 1e-8 ? 0 : 1,
    fmt::format("Check ‖err‖∞ = {:.5}\n",err)
  );
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );
  err = alglin::asum( N, xref1, 1 )/N;
  msg.semaphore(
    err > 1e-8 ? 0 : 1,
    fmt::format("Check ‖err‖_1/N = {:.5}\n",err)
  );
  UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

  fmt::print("\n\ncheck residual\n\n");

  alglin::Copy_negate_n( rhs, N, resid );
  alglin::Copy_n( rhs, N, x );
  BCR.solve( x );
  BCR_SAVED.add_Mv( x, resid );

  real_type res = alglin::nrm2( BCR_SAVED.nrows(), resid, 1 );
  fmt::print("‖res‖₂ = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
  res = alglin::asum( BCR_SAVED.nrows(), resid, 1 );
  fmt::print("‖res‖₁ = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
  res = alglin::absmax( BCR_SAVED.nrows(), resid, 1 );
  fmt::print("‖res‖∞ = {:.5}\n",res);
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

  msg.green("All done!\n");

  }
  catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
  }
  catch (...) {
    msg.green( "Unknown error\n" );
  }

  return 0;
}
