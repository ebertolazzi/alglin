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
#endif

Utils::Console msg( &std::cout,4);

#ifdef LAPACK_WRAPPER_USE_OPENBLAS
#include <omp.h>
#endif

#include <random>

using namespace std;
typedef double real_type;

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

    Utils::ThreadPool1         TP(nth);
    alglin::BorderedCR<double> BCR(&TP);
    alglin::BorderedCR<double> BCR_SAVED(&TP);

    alglin::integer n      = 2;
    alglin::integer nblock = 2;
    alglin::integer qx     = 0;
    alglin::integer qr     = 1;
    alglin::integer nx     = 0;
    alglin::integer nr     = 0;
    alglin::integer Nc     = (nblock+1)*n+nx+qx;
    alglin::integer Nr     = (nblock+1)*n+nr+qr;

    BCR.allocate( nblock, n, qr, qx, nr, nx );

    alglin::Malloc<real_type>       base_value("real");
    alglin::Malloc<alglin::integer> base_index("integer");

    alglin::integer N = std::max(Nr,Nc);
    base_value.allocate( size_t(5*N) );

    real_type * x     = base_value(size_t(N)); // extra space per multiple rhs
    real_type * xref  = base_value(size_t(N));
    real_type * xref1 = base_value(size_t(N));
    real_type * rhs   = base_value(size_t(N));
    real_type * resid = base_value(size_t(N));

    BCR.D(0,0,0) = 1;
    BCR.D(0,0,1) = 2;
    BCR.D(0,1,0) = 3;
    BCR.D(0,1,1) = 3;

    BCR.E(0,0,0) = 4;
    BCR.E(0,0,1) = 2;
    BCR.E(0,1,0) = 3;
    BCR.E(0,1,1) = 3;

    BCR.D(1,0,0) = 4;
    BCR.D(1,0,1) = 0;
    BCR.D(1,1,0) = 3;
    BCR.D(1,1,1) = 3;

    BCR.E(1,0,0) = 0;
    BCR.E(1,0,1) = 1;
    BCR.E(1,1,0) = -1;
    BCR.E(1,1,1) = 3;

    BCR.H(0,0) = 1;
    BCR.H(0,1) = 1;
    BCR.H(1,0) = 1;
    BCR.H(1,1) = 0;
    BCR.H(2,0) = 0;
    BCR.H(2,1) = -1;

    BCR.H(0,2) = 2;
    BCR.H(0,3) = -1;
    BCR.H(1,2) = 1;
    BCR.H(1,3) = 1;
    BCR.H(2,2) = -3;
    BCR.H(2,3) = -1;

    std::cout << "\n\n\n\n\n\n";
    //BCR.select_LU();
    //BCR.select_QR();
    BCR.select_QRP();

    //BCR.select_last_LU();
    //BCR.select_last_LUPQ();
    //BCR.select_last_SVD();
    //BCR.select_last_QR();
    //BCR.select_last_QRP();
    //BCR.select_last_LSS();
    //BCR.select_last_LSY();
    BCR.select_last_PINV();

    for ( alglin::integer i = 0; i < Nc; ++i ) x[i] = i+1;
    alglin::Copy_n( x, Nc, xref );
    BCR.Mv( x, rhs );
    BCR_SAVED.dup( BCR );
    fmt::print(
      "Nr x Nc = {} x {}\n"
      "nblock  = {}\n"
      "n       = {}\n"
      "nr      = {}\n"
      "nx      = {}\n"
      "qr      = {}\n"
      "qx      = {}\n",
      Nr, Nc, nblock, n, nr, nx, qr, qx
    );

    fmt::print(
      "xe  = {}"
      "rhs = {}",
      alglin::print_matrix( 1, Nc, xref, 1 ),
      alglin::print_matrix( 1, Nr, rhs,  1 )
    );

    BCR.factorize();

    alglin::Copy_n( rhs, Nr, x );
    BCR.solve( x );

    fmt::print(
      "x   = {}",
      alglin::print_matrix( 1, Nc, x, 1 )
    );

    msg.green("All done!\n");
  }
  catch ( std::exception const & err ) {
    std::cerr << "Error: " << err.what();
  }
  catch (...) {
    std::cerr << "Unknown error\n";
  }

  return 0;
}
