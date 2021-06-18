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
typedef double real_type;

int
main() {

  try {

    alglin::integer nth = std::thread::hardware_concurrency();

    #ifdef LAPACK_WRAPPER_USE_OPENBLAS
    openblas_set_num_threads(1);
    goto_set_num_threads(1);
    #endif

    Utils::ThreadPool          TP(8);
    alglin::BorderedCR<double> BCR(&TP);
    alglin::BorderedCR<double> BCR_SAVED(&TP);
    //alglin::BorderedCR<double> BCR(nullptr);
    //alglin::BorderedCR<double> BCR_SAVED(nullptr);

    alglin::integer n      = 2;
    alglin::integer nblock = 2;
    alglin::integer qx     = 0;
    alglin::integer qr     = 1;
    alglin::integer nx     = 0;
    alglin::integer nr     = 0;
    alglin::integer Nc     = (nblock+1)*n+nx+qx;
    alglin::integer Nr     = (nblock+1)*n+nr+qr;

    BCR.allocate( nblock, n, qr, qx, nr, nx );

    alglin::Malloc<real_type>       baseValue("real");
    alglin::Malloc<alglin::integer> baseIndex("integer");

    alglin::integer N = std::max(Nr,Nc);
    baseValue.allocate( size_t(5*N) );

    real_type * x     = baseValue(size_t(N)); // extra space per multiple rhs
    real_type * xref  = baseValue(size_t(N));
    real_type * xref1 = baseValue(size_t(N));
    real_type * rhs   = baseValue(size_t(N));
    real_type * resid = baseValue(size_t(N));

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
    //BCR.select_SUPERLU();

    //BCR.select_last_LU();
    //BCR.select_last_LUPQ();
    //BCR.select_last_SVD();
    //BCR.select_last_QR();
    //BCR.select_last_QRP();
    //BCR.select_last_LSS();
    //BCR.select_last_LSY();
    BCR.select_last_PINV();

    for ( alglin::integer i = 0; i < Nc; ++i ) x[i] = i+1;
    std::copy_n( x, Nc, xref );
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
      lapack_wrapper::print_matrix( 1, Nc, xref, 1 ),
      lapack_wrapper::print_matrix( 1, Nr, rhs,  1 )
    );

    BCR.factorize();

    std::copy_n( rhs, Nr, x );
    BCR.solve( x );

    fmt::print(
      "x   = {}",
      lapack_wrapper::print_matrix( 1, Nc, x, 1 )
    );

    fmt::print("All done!\n");
  }
  catch ( std::exception const & err ) {
    std::cerr << "Error: " << err.what();
  }
  catch (...) {
    std::cerr << "Unknown error\n";
  }

  return 0;
}
