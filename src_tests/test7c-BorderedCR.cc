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

int
main( int argc, char *argv[] ) {

  try {

    alglin::integer nth{ static_cast<alglin::integer>(std::thread::hardware_concurrency()) };

    for ( int i{0}; i < argc; ++i )
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

    alglin::integer const n      { 2 };
    alglin::integer const nblock { 3 };
    alglin::integer const qx     { 1 };
    alglin::integer const qr     { 1 };
    alglin::integer const nx     { 2 };
    alglin::integer const nr     { 3 };
    alglin::integer const Nc     { (nblock+1)*n+nx+qx };
    alglin::integer const Nr     { (nblock+1)*n+nr+qr };

    BCR.allocate( nblock, n, qr, qx, nr, nx );

    alglin::Malloc<real_type>       base_value("real");
    alglin::Malloc<alglin::integer> base_index("integer");

    alglin::integer const N{ std::max(Nr,Nc) };
    base_value.allocate( 3 * N );

    real_type * x    { base_value( N ) }; // extra space per multiple rhs
    real_type * xref { base_value( N ) };
    real_type * rhs  { base_value( N ) };

    alglin::integer kkk{0};

    for ( alglin::integer k{0}; k < nblock; ++k ) {
      for ( alglin::integer j{0}; j < n; ++j ) {
        for ( alglin::integer i{0}; i < n; ++i ) {
          BCR.D(k,i,j) = ++kkk;
        }
      }
      for ( alglin::integer j{0}; j < n; ++j ) {
        for ( alglin::integer i{0}; i < n; ++i ) {
          BCR.E(k,i,j) = -(++kkk);
        }
      }
      for ( alglin::integer j{0}; j < nx; ++j ) {
        for ( alglin::integer i{0}; i < n; ++i ) {
          BCR.B(k,i,j) = ++kkk;
        }
      }
    }

    for ( alglin::integer j{0}; j < 2*n+qx+nx; ++j ) {
      for ( alglin::integer i{0}; i < n+qr; ++i ) {
        BCR.H(i,j) = ++kkk;
      }
    }

    for ( alglin::integer k{0}; k <= nblock; ++k ) {
      for ( alglin::integer j{0}; j < n; ++j ) {
        for ( alglin::integer i{0}; i < nr; ++i ) {
          BCR.C(k,i,j) = -(++kkk);
        }
      }
    }  

    for ( alglin::integer j{0}; j < nx; ++j ) {
      for ( alglin::integer i{0}; i < nr; ++i ) {
        BCR.F(i,j) = -(++kkk);
      }
    }

    for ( alglin::integer j{0}; j < qx; ++j ) {
      for ( alglin::integer i{0}; i < nr; ++i ) {
        BCR.Cq(i,j) = ++kkk;
      }
    }
    
    std::ofstream file("test7.m");
    BCR.print_matlab_script( file );
    file.close();

    alglin::real_type M_values[Nr*Nc];
    alglin::integer   M_row[Nr*Nc];
    alglin::integer   M_col[Nr*Nc];
    alglin::integer   M_nnz{Nr*Nc};
    
    kkk = 0;
    for ( alglin::integer j{0}; j < Nc; ++j ) {
      for ( alglin::integer i{0}; i < Nr; ++i ) {
        M_values[kkk] = kkk+1;
        M_row[kkk]    = i;
        M_col[kkk]    = j;
        ++kkk;
      }
    }
    
    BCR.sparse_load( M_values, M_row, 0, M_col, 0, M_nnz, false );

    file.open("test7b.m");
    BCR.print_matlab_script( file );
    file.close();

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

    for ( alglin::integer i{0}; i < Nc; ++i ) x[i] = i+1;
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
