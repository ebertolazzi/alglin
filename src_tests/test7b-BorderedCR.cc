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
  real_type random = real_type(generator())/generator.max();
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

  real_type diag = 1.01*n;

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
main( int argc, char *argv[] ) {

  try {

    alglin::integer nth = std::thread::hardware_concurrency();

    for ( int i = 0; i < argc; ++i )
      fmt::print( "arg[{}] = {}\n", i, argv[i] );

    if ( argc >= 2 ) nth = atoi( argv[1] );

    fmt::print(
      "nthread (available) = {}\n"
      "nthread (used)      = {}\n",
      std::thread::hardware_concurrency(), nth
    );

    for ( alglin::integer iii=0; iii<1; ++iii ) {
      for ( alglin::integer jjj=7; jjj<8; ++jjj ) {

        Utils::ThreadPool1         TP(nth);
        alglin::BorderedCR<double> BCR(&TP);
        alglin::BorderedCR<double> BCR_SAVED(&TP);

        #define NSIZE 2

        alglin::integer n      = NSIZE;
        alglin::integer nblock = 4;
        alglin::integer qx     = 1;// 4+1;
        alglin::integer qr     = 10;// 4;
        alglin::integer nx     = 0;// 2-1;
        alglin::integer nr     = 0;//2;
        alglin::integer Nr     = (nblock+1)*n+nr+qr;
        alglin::integer Nc     = (nblock+1)*n+nx+qx;
        alglin::integer N      = std::max( Nr, Nc );

        fill_matrix( BCR, nblock, n, nr, nx, qr, qx );

        alglin::Malloc<real_type>       base_value("real");
        alglin::Malloc<alglin::integer> base_index("integer");

        base_value.allocate( size_t(5*N) );

        real_type * x     = base_value(size_t(N)); // extra space per multiple rhs
        real_type * xref  = base_value(size_t(N));
        real_type * xref1 = base_value(size_t(N));
        real_type * rhs   = base_value(size_t(N));
        real_type * resid = base_value(size_t(N));

        std::cout << "\n\n\n\n\n\n";

        switch ( iii ) {
        case 0:
          BCR.select_LU();
          fmt::print("use LU\n");
          break;
        case 1:
          BCR.select_QR();
          fmt::print("use QR\n");
          break;
        case 2:
          BCR.select_QRP();
          fmt::print("use QRP\n");
          break;
        }

        switch ( jjj ) {
        case 0:
          BCR.select_last_LU();
          fmt::print("use last LU\n");
          break;
        case 1:
          BCR.select_last_LUPQ();
          fmt::print("use last LUPQ\n");
          break;
        case 2:
          BCR.select_last_SVD();
          fmt::print("use last SVD\n");
          break;
        case 3:
          BCR.select_last_QR();
          fmt::print("use last QR\n");
          break;
        case 4:
          BCR.select_last_QRP();
          fmt::print("use last QRP\n");
          break;
        case 5:
          BCR.select_last_LSS();
          fmt::print("use last LSS\n");
          break;
        case 6:
          BCR.select_last_LSY();
          fmt::print("use last LSY\n");
          break;
        case 7:
          BCR.select_last_PINV();
          fmt::print("use last PINV\n");
          break;
        }

        for ( alglin::integer i = 0; i < Nc; ++i ) x[i] = 1+(i % 100);
        alglin::Copy_n( x, Nc, xref );
        BCR.Mv( x, rhs );
        BCR_SAVED.dup( BCR );

        fmt::print(
          "x = {}"
          "b = {}",
          alglin::print_matrix( 1, Nc, x, 1 ),
          alglin::print_matrix( 1, Nr, rhs, 1 )
        );

        fmt::print(
          "nthread (avilable) = {}\n"
          "nthread (used)     = {}\n",
          std::thread::hardware_concurrency(), nth
        );

        fmt::print(
          "Nr     = {}\n"
          "Nc     = {}\n"
          "nblock = {}\n"
          "n      = {}\n"
          "nr     = {}\n"
          "nx     = {}\n"
          "qr     = {}\n"
          "qx     = {}\n",
          Nr, Nc, nblock, n, nr, nx, qr, qx
        );

        BCR.print_matlab_script( std::cout );

        Utils::TicToc tm;
        tm.tic();
        BCR.factorize();
        tm.toc();
        fmt::print("\nFactorize = {:.5} [ms]\n\n", tm.elapsed_ms());

        alglin::Copy_n( rhs, Nr, x );
        tm.tic();
        BCR.solve( x );
        tm.toc();

        fmt::print(
          "x = {}",
          alglin::print_matrix( 1, Nc, x, 1 )
        );

        fmt::print("\nSolve = {:.5} [ms]\n\n", tm.elapsed_ms() );

        BCR_SAVED.Mv( x, resid );
        alglin::axpy( Nr, -1.0, rhs, 1, resid, 1 );
        real_type err = alglin::asum( Nr, resid, 1 )/Nr;
        msg.semaphore(
          err > 1e-8 ? 0 : 1,
          fmt::format("Check |resid|_1/N = {:.5}\n",err)
        );

        alglin::Copy_n( xref, Nc, xref1 );
        alglin::axpy( Nc, -1.0, x, 1, xref1, 1 );
        err = alglin::absmax( Nc, xref1, 1 );
        msg.semaphore(
          err > 1e-8 ? 0 : 1,
          fmt::format("Check |err|_inf = {:.5}\n",err)
        );
        UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

        err = alglin::asum( Nc, xref1, 1 )/Nc;
        msg.semaphore(
          err > 1e-8 ? 0 : 1,
          fmt::format("Check |err|_1/N = {:.5}\n",err)
        );
        UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

        msg.green( "All done!\n" );
      }
    }
  }
  catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
  }
  catch (...) {
    msg.red( "Unknown error\n" );
  }

  return 0;
}
