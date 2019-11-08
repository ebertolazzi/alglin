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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <vector>

#include "Alglin_Config.hh"
#include "Alglin_aux.hh"
#include "ABD_Block.hh"

#ifdef __GNUC__ 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

typedef double valueType;

#include <random>
static unsigned seed1 = 2;
// std::chrono::system_clock::now().time_since_epoch().count();

static std::mt19937 generator(seed1);
static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}
#define TIC tm.tic()
#define TOC tm.toc()

using namespace std;

int
main() {

  alglin::integer NB = 0;
  for (; NB < 100; NB = NB*2 + 1 ) {
    cout << "\n\n\nNB = " << NB << "\n\n\n\n";

    alglin::integer dim      = 100;
    alglin::integer row0     = 40;
    alglin::integer col00    = 2;
    alglin::integer colNN    = 10;
    alglin::integer rowN     = (dim-row0)+(col00+colNN);
    alglin::integer numBlock = 1000;

    alglin::integer N   = row0 + rowN + numBlock*dim;
    alglin::integer nnz = row0*(dim+col00) +
                          rowN*(dim+colNN) +
                          2*dim*dim + 13*(N+NB) + 2*N*NB + NB*NB;

    alglin::Malloc<valueType> baseValue("real");
    baseValue.allocate(size_t(nnz));

    valueType * block0 = baseValue(size_t(row0*(dim+col00)));
    valueType * blockN = baseValue(size_t(rowN*(dim+colNN)));
    valueType * AdAu   = baseValue(size_t(2*dim*dim));

    valueType * B = baseValue(size_t(N*NB));
    valueType * C = baseValue(size_t(N*NB));
    valueType * D = baseValue(size_t(NB*NB));

    valueType * x      = baseValue(size_t(10*(N+NB)));
    valueType * xref   = baseValue(size_t(N+NB));
    valueType * xref1  = baseValue(size_t(N+NB));
    valueType * rhs    = baseValue(size_t(N+NB));

    alglin::LASTBLOCK_Choice ch[4] = { alglin::LASTBLOCK_LU,
                                       alglin::LASTBLOCK_QR,
                                       alglin::LASTBLOCK_QRP,
                                       alglin::LASTBLOCK_SVD };
    char const * kind[] = { "LU", "QR", "QRP", "SVD" };

    alglin::BlockLU<valueType> LU;
    LU.allocateTopBottom( numBlock, dim, row0, dim+col00, rowN, dim+colNN, NB );

    // carico matrice
    TicToc tm;

    for ( int test = 0; test < 3; ++test ) {
      fmt::print("\n\n\ntest N.{} NB = {} kind = {}\n", test, NB, kind[test]);
      valueType diag = 2*dim;

      alglin::integer nn = row0-col00;

      for ( int k = 0; k < numBlock; ++k ) {
        for ( int i = 0; i < dim; ++i ) {
          for ( int j = 0; j < 2*dim; ++j )
            AdAu[i+j*dim] = rand(-1,1);
          AdAu[i*(dim+1)+nn*dim] += diag;
        }
        LU.loadBlock( k, AdAu, dim );
      }
      for ( int i = 0; i < row0; ++i ) {
        for ( int j = 0; j < dim+col00; ++j )
          block0[i+j*row0] = rand(-1,1);
        block0[i*(row0+1)] += diag;
      }
      for ( int i = 0; i < rowN; ++i ) {
        for ( int j = 0; j < dim+colNN; ++j )
          blockN[i+j*rowN] = rand(-1,1);
        blockN[i*(rowN+1)+nn*rowN] += diag;
      }
      for ( int j = 0; j < NB; ++j ) {
        for ( int i = 0; i < N; ++i ) B[i+j*N] = 1;
      }
      for ( int i = 0; i < NB; ++i ) {
        for ( int j = 0; j < N; ++j ) C[i+j*NB] = rand(-1,1);
        C[i*(NB+1)] += 10;
      }
      for ( int i = 0; i < NB; ++i )
        for ( int j = 0; j < NB; ++j )
          D[i+j*NB] = rand(-1,1)+(i==j?10:0);

      LU.loadRightBlocks( B, N );
      LU.loadBottomBlocks( C, NB );
      LU.loadRBblock( D, NB );
      LU.loadTopBottom( block0, row0, blockN, rowN );
      LU.selectLastBlockSolver( ch[test] );

      //ofstream file("dump_mat.txt");
      //LU.dump_ccoord( file );
      //file.close();

      fmt::print(
        "N = {} n = {} col00 = {} colNN = {} row0 = {} rowN {}\n",
        N, dim, col00, colNN, row0, rowN
      );

      for ( alglin::integer i = 0; i < N+NB; ++i ) x[i] = 1+(i%4);
      std::copy( x, x+N+NB, xref  );
      std::copy( x, x+N+NB, xref1 );
      LU.Mv( x, rhs );

      TIC;
      LU.factorize_bordered();
      TOC;
      fmt::print("(Block {}) Factorize = {:10.5} [ms]\n", kind[test], tm.elapsed_ms());

      std::copy( rhs, rhs+N+NB, x );
      TIC;
      LU.solve_bordered( x );
      TOC;
      fmt::print("(Block {}) Solve     = {:10.5} [ms]\n", kind[test], tm.elapsed_ms());

      alglin::axpy( N+NB, -1.0, x, 1, xref, 1 );
      valueType err = alglin::absmax( N+NB, xref, 1 );
      fmt::print("\nCheck |err|_inf = {:10.5}\n\n",err);
      LW_ASSERT0( err < 1e-8, "test failed!" );

      for ( alglin::integer i = 0; i < 10; ++i ) std::copy( rhs, rhs+N+NB, x+i*(N+NB) );
      TIC;
      LU.solve_bordered( 1, x, N+NB );
      TOC;
      fmt::print("(Block {}) Solve     = {:10.5} [ms]\n", kind[test], tm.elapsed_ms());

      alglin::axpy( N+NB, -1.0, x, 1, xref1, 1 );
      err = alglin::absmax( N+NB, xref1, 1 );
      fmt::print("\nCheck |err|_inf = {:10.5}\n\n",err);
      LW_ASSERT0( err < 1e-8, "test failed!" );

    }
  }

  fmt::print("\n\nAll done!\n");

  return 0;
}
