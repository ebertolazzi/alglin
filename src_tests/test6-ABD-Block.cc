/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017c                                                     |
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
 |      Università degli Studi di Trento                                    |
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

#include <random>

typedef double real_type;

#include <random>
static unsigned seed1 = 2;
// std::chrono::system_clock::now().time_since_epoch().count();

static std::mt19937 generator(seed1);
static
real_type
rand( real_type const xmin, real_type const xmax ) {
  real_type const random{ static_cast<real_type>(generator())/generator.max() };
  return xmin + (xmax-xmin)*random;
}
#define TIC tm.tic()
#define TOC tm.toc()

using namespace std;

int
main() {

  alglin::integer NB{0};
  for (; NB < 100; NB = NB*2 + 1 ) {
    cout << "\n\n\nNB = " << NB << "\n\n\n\n";

    alglin::integer       dim      {100};
    alglin::integer       row0     {40};
    alglin::integer       col00    {2};
    alglin::integer       colNN    {10};
    alglin::integer       rowN     {(dim-row0)+(col00+colNN)};
    alglin::integer const numBlock {1000};

    alglin::integer       N   { row0 + rowN + numBlock*dim };
    alglin::integer const nnz { row0*(dim+col00) +
                                rowN*(dim+colNN) +
                                2*dim*dim + 13*(N+NB) + 2*N*NB + NB*NB };

    alglin::Malloc<real_type> base_value("real");
    base_value.allocate(nnz);

    real_type * block0 { base_value(row0*(dim+col00)) };
    real_type * blockN { base_value(rowN*(dim+colNN)) };
    real_type * AdAu   { base_value(2*dim*dim)        };

    real_type * B { base_value( N*NB ) };
    real_type * C { base_value( N*NB ) };
    real_type * D { base_value( NB*NB ) };

    real_type * x     { base_value(10*(N+NB) ) };
    real_type * xref  { base_value(N+NB ) };
    real_type * xref1 { base_value(N+NB ) };
    real_type * rhs   { base_value(N+NB ) };

    constexpr alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice ch[]{
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::LU,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::LUPQ,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::QR,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::QRP,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::SVD,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::LSS,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::LSY,
      alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice::PINV
    };
    char const * kind[]{ "LU", "LUPQ", "QR", "QRP", "SVD", "LSS", "LSY", "PINV" };

    alglin::BlockLU<real_type> LU;
    LU.allocate_top_bottom( numBlock, dim, row0, dim+col00, rowN, dim+colNN, NB );

    // carico matrice
    Utils::TicToc tm;

    for ( int test{0}; test < 8; ++test ) {
      fmt::print("\n\n\ntest N.{} NB = {} kind = {}\n", test, NB, kind[test]);
      real_type const diag{ static_cast<double>(2*dim) };

      alglin::integer const nn{ row0-col00 };

      for ( int k{0}; k < numBlock; ++k ) {
        for ( int i{0}; i < dim; ++i ) {
          for ( int j{0}; j < 2*dim; ++j )
            AdAu[i+j*dim] = rand(-1,1);
          AdAu[i*(dim+1)+nn*dim] += diag;
        }
        LU.load_block( k, AdAu, dim );
      }
      for ( int i{0}; i < row0; ++i ) {
        for ( int j{0}; j < dim+col00; ++j )
          block0[i+j*row0] = rand(-1,1);
        block0[i*(row0+1)] += diag;
      }
      for ( int i{0}; i < rowN; ++i ) {
        for ( int j{0}; j < dim+colNN; ++j )
          blockN[i+j*rowN] = rand(-1,1);
        blockN[i*(rowN+1)+nn*rowN] += diag;
      }
      for ( int j{0}; j < NB; ++j ) {
        for ( int i{0}; i < N; ++i ) B[i+j*N] = 1;
      }
      for ( int i{0}; i < NB; ++i ) {
        for ( int j{0}; j < N; ++j ) C[i+j*NB] = rand(-1,1);
        C[i*(NB+1)] += 10;
      }
      for ( int i{0}; i < NB; ++i )
        for ( int j{0}; j < NB; ++j )
          D[i+j*NB] = rand(-1,1)+(i==j?10:0);

      LU.load_right_blocks( B, N );
      LU.load_bottom_blocks( C, NB );
      LU.load_RB_block( D, NB );
      LU.load_top_bottom( block0, row0, blockN, rowN );
      LU.select_last_block_solver( ch[test] );

      //ofstream file("dump_mat.txt");
      //LU.dump_ccoord( file );
      //file.close();

      fmt::print(
        "N = {} n = {} col00 = {} colNN = {} row0 = {} rowN {}\n",
        N, dim, col00, colNN, row0, rowN
      );

      for ( alglin::integer i{0}; i < N+NB; ++i ) x[i] = 1+(i%4);
      alglin::Copy_n( x, N+NB, xref  );
      alglin::Copy_n( x, N+NB, xref1 );
      LU.Mv( x, rhs );

      TIC;
      LU.factorize_bordered();
      TOC;
      fmt::print("(Block {}) Factorize = {:.5} [ms]\n", kind[test], tm.elapsed_ms());

      alglin::Copy_n( rhs, N+NB, x );
      TIC;
      LU.solve_bordered( x );
      TOC;
      fmt::print("(Block {}) Solve = {:.5} [ms]\n", kind[test], tm.elapsed_ms());

      alglin::axpy( N+NB, -1.0, x, 1, xref, 1 );
      real_type err = alglin::absmax( N+NB, xref, 1 );
      msg.semaphore(
        err > 1e-8 ? 0 : 1,
        fmt::format("\nCheck ‖err‖∞ = {:.5}\n\n",err)
      );
      UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

      for ( alglin::integer i{0}; i < 10; ++i ) alglin::Copy_n( rhs, N+NB, x+i*(N+NB) );
      TIC;
      LU.solve_bordered( 1, x, N+NB );
      TOC;
      fmt::print("(Block {}) Solve = {:.5} [ms]\n", kind[test], tm.elapsed_ms());

      alglin::axpy( N+NB, -1.0, x, 1, xref1, 1 );
      err = alglin::absmax( N+NB, xref1, 1 );
      msg.semaphore(
        err > 1e-8 ? 0 : 1,
        fmt::format("\nCheck ‖err‖∞ = {:.5}\n\n",err)
      );
      UTILS_ASSERT0( err < 1e-8, "test failed!\n" );

    }
  }

  fmt::print("\n\nAll done!\n");

  return 0;
}
