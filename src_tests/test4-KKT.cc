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

using namespace std;
typedef double real_type;
using alglin::integer;
using alglin::Transposition;

static
void
test0() {

  alglin::KKT<real_type> kkt;
  integer const N = 2;
  integer const M = 1;

  real_type A[] = {
    1,      0,
    0,      1,
  };

  //
  //  A B
  //  C D
  //
  real_type B[] = { 0, 0 };
  real_type C[] = { 0, 0 };
  real_type D[] = { 2 };

  // 1 0 0 -> 1
  // 0 1 0 -> 2
  // 0 0 2 -> 6
  kkt.load(
    N, M,
    A, N, false, // N x N
    B, N, false, // N x M
    C, N, false, // M x N
    D, M, false  // M x M
  );
  kkt.factorize();
  real_type x[N+M] = { 1, 2, 3 };
  real_type rhs[N+M];

  alglin::gemv(
    Transposition::NO,
    N, N, 1.0, A, N,
    x, 1,
    0,
    rhs, 1
  );
  alglin::gemv(
    Transposition::NO,
    N, M, 1.0, B, N,
    x+N, 1,
    1,
    rhs, 1
  );
  alglin::gemv(
    Transposition::NO,
    M, N, 1.0, C, M,
    x, 1,
    0,
    rhs+N, 1
  );
  alglin::gemv(
    Transposition::NO,
    M, M, 1.0, D, M,
    x+N, 1,
    1,
    rhs+N, 1
  );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("rhs[{}] = {}\n",i,rhs[i]);
  kkt.solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);
}

static
void
test1() {

  alglin::KKT<real_type> kkt;
  integer const N = 2;
  integer const M = 1;

  real_type A[] = {
    1,      3,
   -2,      1,
  };

  real_type B[] = { 1, -1 };
  real_type C[] = { -1, 1 };
  real_type D[] = { 2 };

  //  1 -2  1 -> 0
  //  3  1 -1 -> 2
  // -1  1  2 -> 7

  kkt.load(
    N, M,
    A, N, false, // N x N
    B, N, false, // N x M
    C, N, false, // M x N
    D, M, false  // M x M
  );
  kkt.factorize();
  real_type x[N+M] = { 1, 2, 3 };
  real_type rhs[N+M];

  alglin::gemv(
    Transposition::NO,
    N, N, 1.0, A, N,
    x, 1,
    0,
    rhs, 1
  );
  alglin::gemv(
    Transposition::NO,
    N, M, 1.0, B, N,
    x+N, 1,
    1,
    rhs, 1
  );
  alglin::gemv(
    Transposition::NO,
    M, N, 1.0, C, M,
    x, 1,
    0,
    rhs+N, 1
  );
  alglin::gemv(
    Transposition::NO,
    M, M, 1.0, D, M,
    x+N, 1,
    1,
    rhs+N, 1
  );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("rhs[{}] = {}\n",i,rhs[i]);
  kkt.solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);
}

static
void
test2() {

  using namespace alglin;

  alglin::KKT<real_type> kkt;
  integer const N = 3;
  integer const M = 2;

  real_type A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  };

  real_type B[] = {
    0.001,      3,
    0.001,      -0.001,
    1,          2,
  };

  real_type C[] = {
    2,       3,
    -0.001,  -0.001,
    1,        2,
  };

  real_type D[] = {
    1,      4,
    -1,      1
  };

  kkt.load(
    N, M,
    A, N, false, // N x N
    B, N, false, // N x M
    C, N, false, // M x N
    D, M, false  // M x M
  );
  kkt.factorize();
  real_type x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 };
  real_type rhs[2*(N+M)];

  gemv(
    Transposition::NO,
    N, N, 1.0, A, N,
    x, 1,
    0,
    rhs, 1
  );
  gemv(
    Transposition::NO,
    N, M, 1.0, B, N,
    x+N, 1,
    1,
    rhs, 1
  );
  gemv(
    Transposition::NO,
    M, N, 1.0, C, M,
    x, 1,
    0,
    rhs+N, 1
  );
  gemv(
    Transposition::NO,
    M, M, 1.0, D, M,
    x+N, 1,
    1,
    rhs+N, 1
  );
  alglin::Copy_n( rhs, N+M, rhs+N+M );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("rhs[{}] = {}\n",i,rhs[i]);
  //kkt.solve( rhs );
  kkt.solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);

}

static
void
test3() {

  using namespace alglin;

  KKT<real_type> kkt;
  integer const N = 3;
  integer const M = 2;

  real_type A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  };

  real_type B[] = {
    0.001,      3,
    0.001,     -0.001,
    1,          2,
  };

  real_type C[] = {
     2,       3,
    -0.001,  -0.001,
     1,       2,
  };

  real_type D[] = {
    1,      4,
    -1,      1
  };

  kkt.load(
    N, M,
    A, N, false, // N x N
    B, N, false, // N x M
    C, N, false, // M x N
    D, M, false  // M x M
  );
  kkt.factorize();
  real_type x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 };
  real_type rhs[2*(N+M)];

  gemv(
    Transposition::YES,
    N, N, 1.0, A, N,
    x, 1,
    0,
    rhs, 1
  );
  gemv(
    Transposition::YES,
    M, N, 1.0, C, M,
    x+N, 1,
    1,
    rhs, 1
  );
  gemv(
    Transposition::YES,
    M, M, 1.0, D, M,
    x+N, 1,
    0,
    rhs+N, 1
  );
  gemv(
    Transposition::YES,
    N, M, 1.0, B, N,
    x, 1,
    1,
    rhs+N, 1
  );

  Copy_n( rhs, N+M, rhs+N+M );

  for ( integer i = 0; i < N+M; ++i )
    fmt::print("rhs[{}] = {}\n",i,rhs[i]);
  kkt.t_solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);
}

static
void
test4() {

  using namespace alglin;

  KKT<real_type> kkt;
  integer const N = 10;
  integer const M = 2;
  BandedLU<real_type> bLU;

  bLU.setup( N, N, 3, 2 ); // number of upper diagonal
  bLU.zero();

  integer ldAA = 2;
  real_type AA[] = {
    3, 1,
    1, 3,
  };
  for ( integer i = 0; i < N; ++i ) bLU(i,i) = i+1;
  bLU.load_block( 2, 2, AA, ldAA, 0, 0 );
  bLU.load_block( 2, 2, AA, ldAA, 3, 3 );
  bLU.load_block( 2, 2, AA, ldAA, 6, 6 );
  bLU.load_block( 1, 1, AA, ldAA, 9, 9 );

  //bLU.dump( cout );

  real_type B[] = {
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

  real_type C[] = {
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

  real_type D[] = {
    1,      4,
    -1,      1
  };

  real_type x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2 };
  real_type rhs[2*(N+M)];

  for ( integer i = 0; i < N+M; ++i ) rhs[i] = 0;
  bLU.aAxpy( 1.0, x, rhs );

  gemv(
    Transposition::YES,
    M, N, 1.0, C, M,
    x+N, 1,
    1,
    rhs, 1
  );

  gemv(
    Transposition::YES,
    M, M, 1.0, D, M,
    x+N, 1,
    0,
    rhs+N, 1
  );

  gemv(
    Transposition::YES,
    N, M, 1.0, B, N,
    x, 1,
    1,
    rhs+N, 1
  );

  Copy_n( rhs, N+M, rhs+N+M );

  for ( integer i = 0; i < N+M; ++i )
    fmt::print("rhs[{}] = {}\n",i,rhs[i]);

  // must be factorized before to call kkt.factorize
  bLU.factorize( "bLU" );
  kkt.load(
    N, M,
    &bLU,
    B, N, false,
    C, M, false,
    D, M, false
  );
  kkt.factorize();
  kkt.t_solve( 2, rhs, N+M );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);
}

static
void
test5() {

  using namespace alglin;

  KKT<real_type> kkt;
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

  real_type vals[] = {
    2, 1, -1, 1, 1, -1, 2, 2, 1, 2, 2,
    3, 2, 4, 3, 5, 3, 6, 1, 1, 1, 1
  };

  // must be factorized before to call kkt.factorize
  kkt.load_banded(
    N, M, 4, 4,
    vals,
    ii, -1,
    jj, -1,
    nnz, true
  );
  kkt.factorize();

  //real_type x[]   = { 1, 2, 3, 4, 5, -1, -2, -3, -4, -5 };
  real_type rhs[] = { -2, -9, -5, -10, 3, -21, 0, 32, 8 };

  //kkt.t_solve( 1, rhs, N+M );
  kkt.t_solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n';
}

static
void
test6() {

  using namespace alglin;

  KKT<real_type> kkt;
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

  real_type vals[] = {
    2, 1, -1, 1, 1, -1, 2, 2, 1, 2, 2,
    3, 2, 4, 3, 5, 3, 6, 1, 1, 1, 1
  };

  // must be factorized before to call kkt.factorize
  integer kblocks[] = { 0, 2, 4, 7 };
  kkt.load_triblock(
    N, M,
    3, kblocks,
    vals,
    ii, -1,
    jj, -1,
    nnz, true
  );
  kkt.factorize();

  //real_type x[]   = { 1, 2, 3, 4, 5, -1, -2, -3, -4, -5 };
  real_type rhs[] = { -2, -9, -5, -10, 3, -21, 0, 32, 8 };

  //kkt.t_solve( 1, rhs, N+M );
  kkt.t_solve( rhs );
  for ( integer i = 0; i < N+M; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);
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
    cerr << exc.what();
  } catch ( ... ) {
    cerr << "Errore Sconosciuto!\n";
  }

  cout << "\nAll done!\n";

  return 0;
}
