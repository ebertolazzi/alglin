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
#include <fstream>
#include <vector>

#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"
#include "KKT_like.hh"

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
using alglin::integer;

static
void
test0() {

  alglin::KKT<valueType> kkt;
  integer const N = 2;
  integer const M = 1;

  valueType A[] = {
    1,      0,
    0,      1,
  };

  valueType B[] = { 0, 0 };
  valueType C[] = { 0, 0 };
  valueType D[] = { 2 };
  
  // 1 0 0 -> 1
  // 0 1 0 -> 2
  // 0 0 2 -> 6

  kkt.load( N, M,
            A, N, false,
            B, M, false,
            C, N, false,
            D, M, false );
  kkt.factorize();
  valueType x[N+M] = { 1, 2, 3 };
  valueType rhs[N+M];

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 );
  for ( integer i = 0; i < N+M; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n';
  kkt.solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';

}

static
void
test1() {

  alglin::KKT<valueType> kkt;
  integer const N = 2;
  integer const M = 1;

  valueType A[] = {
    1,      3,
   -2,      1,
  };

  valueType B[] = { 1, -1 };
  valueType C[] = { -1, 1 };
  valueType D[] = { 2 };
  
  //  1 -2  1 -> 0
  //  3  1 -1 -> 2
  // -1  1  2 -> 7

  kkt.load( N, M,
            A, N, false,
            B, M, false,
            C, N, false,
            D, M, false );
  kkt.factorize();
  valueType x[N+M] = { 1, 2, 3 };
  valueType rhs[N+M];

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 );
  for ( integer i = 0; i < N+M; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n';
  kkt.solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';

}

static
void
test2() {

  alglin::KKT<valueType> kkt;
  integer const N = 3;
  integer const M = 2;

  valueType A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  };

  valueType B[] = {
    0.001,      3,
    0.001,      -0.001,
    1,          2,
  };

  valueType C[] = {
    2,       3,
    -0.001,  -0.001,
    1,        2,
  };

  valueType D[] = {
    1,      4,
    -1,      1
  };
  
  kkt.load( N, M,
            A, N, false,
            B, N, false,
            C, M, false,
            D, M, false );
  kkt.factorize();
  valueType x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 };
  valueType rhs[2*(N+M)];

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 );
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 );
  std::copy( rhs, rhs + N+M, rhs+N+M );
  for ( integer i = 0; i < N+M; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n';
  //kkt.solve( rhs );
  kkt.solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';

}

static
void
test3() {

  alglin::KKT<valueType> kkt;
  integer const N = 3;
  integer const M = 2;

  valueType A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  };

  valueType B[] = {
    0.001,      3,
    0.001,     -0.001,
    1,          2,
  };

  valueType C[] = {
     2,       3,
    -0.001,  -0.001,
     1,       2,
  };

  valueType D[] = {
    1,      4,
    -1,      1
  };
  
  kkt.load( N, M,
            A, N, false,
            B, N, false,
            C, M, false,
            D, M, false );
  kkt.factorize();
  valueType x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 };
  valueType rhs[2*(N+M)];

  alglin::gemv(alglin::TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 );
  alglin::gemv(alglin::TRANSPOSE,
               M, N, 1.0, C, M,
               x+N, 1,
               1,
               rhs, 1 );

  alglin::gemv(alglin::TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               0,
               rhs+N, 1 );

  alglin::gemv(alglin::TRANSPOSE,
               N, M, 1.0, B, N,
               x, 1,
               1,
               rhs+N, 1 );

  std::copy( rhs, rhs+N+M, rhs+N+M );

  for ( integer i = 0; i < N+M; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n';
  kkt.t_solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';
}

static
void
test4() {

  alglin::KKT<valueType> kkt;
  integer const N = 10;
  integer const M = 2;
  alglin::BandedLU<valueType> bLU;

  bLU.setup( N, N, 3, 2 ); // number of upper diagonal
  bLU.zero();

  integer ldAA = 2;
  valueType AA[] = {
    3, 1,
    1, 3,
  };
  for ( integer i = 0; i < N; ++i ) bLU(i,i) = i+1;
  bLU.load_block( 2, 2, AA, ldAA, 0, 0 );
  bLU.load_block( 2, 2, AA, ldAA, 3, 3 );
  bLU.load_block( 2, 2, AA, ldAA, 6, 6 );
  bLU.load_block( 1, 1, AA, ldAA, 9, 9 );

  //bLU.dump( cout );

  valueType B[] = {
    0.001,      3,
    0.001,     -0.001,
    1,          2,
    1,          2,
    1,          2,
    0.001,      3,
    0.001,     -0.001,
    2,          3,
    0.001,      3,
    0.001,     -0.001,
  };

  valueType C[] = {
    1,          2,
    1,          2,
    0.001,      3,
    2,          3,
    -0.001,  -0.001,
    1,          2,
    1,          2,
    0.001,      3,
    1,          2,
    0.001,      3,
  };

  valueType D[] = {
    1,      4,
    -1,      1
  };

  valueType x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 };
  valueType rhs[2*(N+M)];

  for ( integer i = 0; i < N+M; ++i ) rhs[i] = 0;
  bLU.aAxpy( 1.0, x, rhs );

  alglin::gemv(alglin::TRANSPOSE,
               M, N, 1.0, C, M,
               x+N, 1,
               1,
               rhs, 1 );

  alglin::gemv(alglin::TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               0,
               rhs+N, 1 );

  alglin::gemv(alglin::TRANSPOSE,
               N, M, 1.0, B, N,
               x, 1,
               1,
               rhs+N, 1 );

  std::copy( rhs, rhs+N+M, rhs+N+M );

  for ( integer i = 0; i < N+M; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n';

  // must be factorized before to call kkt.factorize
  bLU.factorize();
  kkt.load( N, M,
            &bLU,
            B, N, false,
            C, M, false,
            D, M, false );
  kkt.factorize();
  kkt.t_solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';
}

static
void
test5() {

  alglin::KKT<valueType> kkt;
  integer const N   = 7;
  integer const M   = 2;
  integer const nnz = 22;

  integer ii[] = {
    1, 2, 3, 8, 2, 4, 9, 3, 5, 6, 7,
    8, 4, 9, 5, 8, 6, 9, 7, 9, 8, 9
  };

  integer jj[] = {
    1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3,
    3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9
  };

  valueType vals[] = {
    2, 1, -1, 1, 1, -1, 2, 2, 1, 2, 2,
    3, 2, 4, 3, 5, 3, 6, 1, 1, 1, 1
  };

  // must be factorized before to call kkt.factorize
  kkt.load_banded( N, M, 4, 4,
                   vals,
                   ii, -1,
                   jj, -1,
                   nnz, true );
  kkt.factorize();

  valueType x[]   = { 1, 2, 3, 4, 5, -1, -2, -3, -4, -5 };
  valueType rhs[] = { -2, -9, -5, -10, 3, -21, 0, 32, 8 };

  //kkt.t_solve( 1, rhs, N+M );
  kkt.t_solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
   cout << "x[" << i << "] = " << rhs[i] << '\n';
}

static
void
test6() {

  alglin::KKT<valueType> kkt;
  integer const N   = 7;
  integer const M   = 2;
  integer const nnz = 22;

  integer ii[] = {
    1, 2, 3, 8, 2, 4, 9, 3, 5, 6, 7,
    8, 4, 9, 5, 8, 6, 9, 7, 9, 8, 9
  };

  integer jj[] = {
    1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3,
    3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 9
  };

  valueType vals[] = {
    2, 1, -1, 1, 1, -1, 2, 2, 1, 2, 2,
    3, 2, 4, 3, 5, 3, 6, 1, 1, 1, 1
  };

  // must be factorized before to call kkt.factorize
  integer kblocks[] = { 0, 2, 4, 7 };
  kkt.load_triblock( N, M,
                     3, kblocks,
                     vals,
                     ii, -1,
                     jj, -1,
                     nnz, true );
  kkt.factorize();

  valueType x[]   = { 1, 2, 3, 4, 5, -1, -2, -3, -4, -5 };
  valueType rhs[] = { -2, -9, -5, -10, 3, -21, 0, 32, 8 };

  //kkt.t_solve( 1, rhs, N+M );
  kkt.t_solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
   cout << "x[" << i << "] = " << rhs[i] << '\n';
}

int
main() {

  try {
    cout << "\n\n\ntest0\n";
    test0();
    cout << "\n\n\ntest1\n";
    test1();
    cout << "\n\n\ntest2\n";
    test2();
    cout << "\n\n\ntest3\n";
    test3();
    cout << "\n\n\ntest4\n";
    test4();
    cout << "\n\n\ntest5\n";
    test5();
    cout << "\n\n\ntest6\n";
    test6();
  } catch ( exception const & exc ) {
    cerr << exc.what() << '\n';
  } catch ( ... ) {
    cerr << "Errore Sconosciuto!\n";
  }

  cout << "\nAll done!\n";

  return 0;
}
