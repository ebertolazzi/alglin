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

using namespace std;
typedef double real_type;
using alglin::integer;
using alglin::Transposition;

static
void
test1() {
  alglin::QR<real_type> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  real_type A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  msg.green(
    fmt::format(
      "\n\n\nTest1:\n\nInitial A\n{}",
      alglin::print_matrix( M, N, A, M )
    )
  );

  fmt::print( "\nDo QR factorization of A^T\n" );
  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print( "\nR=\n{}\n", alglin::print_matrix( M, M, R, M ) );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(
    "\nLS solution of A x = b\n"
    "b^T = {}",
    alglin::print_matrix( 1, M, b, 1 )
  );

  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, M, x );
  alglin::Zero_n( x+3, N-M );
  qr.Q_mul( x );

  fmt::print( "x^T = {}", alglin::print_matrix( 1, M, x, 1 ) );

  alglin::gemv( Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1 );

  fmt::print(
    "residual = {}done test1\n",
    alglin::print_matrix( 1, M, b, 1 )
  );
}



static
void
test2() {
  alglin::QRP<real_type> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  real_type A[] = {
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  fmt::print(
    "\n\n\nTest2:\n\nInitial A\n"
    "{}"
    "\nDo QR factorization of A^T\n",
    alglin::print_matrix( M, N, A, M )
  );
  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print(
    "\nR=\n{}",
    alglin::print_matrix( M, M, R, M )
  );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  alglin::gemv(
    Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  alglin::gemv(
    Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  fmt::print(
    "\nLS solution of A x = b\n\n"
    "b^T = {}",
    alglin::print_matrix( 1, M, b, 1 )
  );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, M, x );
  alglin::Zero_n( x+3, N-M );
  qr.Q_mul( x );

  fmt::print( "x^T = {}", alglin::print_matrix( 1, M, x, 1 ) );

  alglin::gemv(
    Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1
  );
  fmt::print(
    "residual = {}done test2\n",
    alglin::print_matrix( 1, M, b, 1 )
  );
}

static
void
test3() {
  alglin::QRP<real_type> qr;
  integer const M   = 5;
  integer const N   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,     2, 3,
    0.001,  0.001,     0, 0.001, 1e-10,
    0,      0.001,     0, 0.001, 1e-12,
    0.001,     -1, 1e-12,    -1, -1e-12,
    0.000001,   5,     3,     5, 3
  };

  fmt::print(
    "\n\n\nTest3:\n\nInitial A\n"
    "{}"
    "\nDo QR factorization of A^T\n",
    alglin::print_matrix( M, N, A, M )
  );

  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  fmt::print(
    "\nR=\n{}",
    alglin::print_matrix( M, M, R, M )
  );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(
    "\nLS solution of A x = b\n\n"
    "b^T = {}",
    alglin::print_matrix( 1, M, b, 1 )
  );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, 3, x );
  alglin::Zero_n( x+3, 2 );
  qr.Q_mul( x );

  fmt::print(
    "x^T = {}",
    alglin::print_matrix( 1, M, x, 1 )
  );

  alglin::gemv( Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  fmt::print(
    "residual = {}done test3\n",
    alglin::print_matrix( 1, M, b, 1 )
  );
}

#define TEST4(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  alglin::Copy_n( rhs, M, x ); \
  alglin::Copy_n( rhs, M, b ); \
  /* L.solve( x ); */ \
  F.solve( 1, x, M); \
  cout << "x^T      =" << alglin::print_matrix( 1, M, x, 1 ); \
  \
  alglin::gemv( Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual =" << alglin::print_matrix( 1, M, b, 1 ); \
  res = alglin::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

static
void
test4() {
  alglin::LU<real_type>   lu;
  alglin::LUPQ<real_type> lupq;
  alglin::QR<real_type>   qr;
  alglin::QRP<real_type>  qrp;
  alglin::SVD<real_type>  svd;
  alglin::LSS<real_type>  lss;
  alglin::LSY<real_type>  lsy;
  alglin::PINV<real_type> pinv;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M], res;
  real_type x[M] = {1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(
    "\n\n\nTest4:\n\nInitial A\n{}",
    alglin::print_matrix( M, M, A, M )
  );

  TEST4("LU",lu);
  TEST4("LUPQ",lupq);
  TEST4("QR",qr);
  TEST4("QRP",qrp);
  TEST4("SVD",svd);
  TEST4("LSS",lss);
  TEST4("LSY",lsy);
  TEST4("PINV",pinv);

  fmt::print( "done test4\n" );
}

#define TEST5(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  alglin::Copy_n( rhs, M, x ); \
  alglin::Copy_n( rhs, M, b ); \
  /* F.t_solve( x ); */ \
  F.t_solve( 1, x, M); \
  cout << "x^T      = " << alglin::print_matrix( 1, M, x, 1 ); \
  \
  alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual = " << alglin::print_matrix( 1, M, b, 1 ); \
  res = alglin::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );


static
void
test5() {
  alglin::LU<real_type>   lu;
  alglin::LUPQ<real_type> lupq;
  alglin::QR<real_type>   qr;
  alglin::QRP<real_type>  qrp;
  alglin::SVD<real_type>  svd;
  alglin::LSS<real_type>  lss;
  alglin::LSY<real_type>  lsy;
  alglin::PINV<real_type> pinv;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,     3+1
  };

  real_type rhs[M], b[M], res;
  real_type x[M] = {1,2,3,4,5};
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(
    "\n\n\nTest5:\n\nInitial A\n{}",
    alglin::print_matrix( M, M, A, M )
  );

  TEST5("LU",lu);
  TEST5("LUPQ",lupq);
  TEST5("QR",qr);
  TEST5("QRP",qrp);
  TEST5("SVD",svd);
  TEST5("LSS",lss);
  TEST5("LSY",lsy);
  TEST5("PINV",pinv);

  fmt::print( "\ndone test5\n" );
}


static
void
test6() {
  alglin::TridiagonalLU<real_type> lu;
  alglin::TridiagonalQR<real_type> qr;

  integer const N = 5;
  real_type const D[] = { 1, 1, 2, -0.1, 0.1 };
  real_type const L[] = { -0.1, -1, -2, -0.1 };
  real_type const U[] = { -1, -10, -2, 0.1 };

  real_type rhs[N], b[N];
  real_type x[N] = {1,2,3,4,5};
  qr.axpy( N, 1.0, L, D, U, x, 0.0, rhs );
  qr.axpy( N, 1.0, L, D, U, x, 0.0, b   );

  msg.green( "\n\n\nTest6:\n\nInitial A\n" );

  msg.yellow( "\n\nDo (trid) LU  factorization of A\n" );
  lu.factorize( "lu", N, L, D, U );

  msg.yellow( "(trid) LU solution of A x = b\n" );
  alglin::Copy_n( rhs, N, x );
  alglin::Copy_n( rhs, N, b );
  lu.solve( x );
  //qr.t_solve( 1, x, M);
  fmt::print( "x^T = {}", alglin::print_matrix( 1, N, x, 1 ) );

  lu.axpy( N, -1.0, L, D, U, x, 1.0, b );
  fmt::format(
    "residual = {}"
    "\n\nDo (trid) QR  factorization of A\n",
    alglin::print_matrix( 1, N, b, 1 )
  );
  qr.factorize( "qr", N, L, D, U );

  msg.yellow( "(trid) QR solution of A x = b\n" );
  alglin::Copy_n( rhs, N, x );
  alglin::Copy_n( rhs, N, b );
  qr.solve( x );
  //qr.t_solve( 1, x, M);
  fmt::format(
    "x^T = {}",
    alglin::print_matrix( 1, N, x, 1 )
  );

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );
  fmt::format(
    "residual = {}\ndone test6\n",
    alglin::print_matrix( 1, N, b, 1 )
  );
}

static
void
test7() {
  alglin::QRP<real_type> qrp;

  integer const M   = 5;
  integer const LDA = 5;
  real_type A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M];
  real_type x[M] = {1,2,3,4,5};
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  msg.green(
    fmt::format(
      "\n\n\nTest7:\n\nInitial A\n"
      "{}"
      "\n\nDo QRP factorization of A\n",
      alglin::print_matrix( M, M, A, M )
    )
  );
  qrp.factorize( "qrp", M, M, A, LDA );

  fmt::print( "QRP solution of A x = b\n" );
  alglin::Copy_n( rhs, M, x );
  alglin::Copy_n( rhs, M, b );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  fmt::print(
    "x^T = {}",
    alglin::print_matrix( 1, M, x, 1 )
  );

  alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  fmt::print(
    "residual = {}"
    "\n\nDo QRP factorization of A\n",
    alglin::print_matrix( 1, M, b, 1 )
  );

  alglin::Matrix<real_type> mat;
  mat.setup( M, M );
  mat.load_block( 2, 5, A,   LDA, 0, 0 );
  mat.load_block( 3, 5, A+2, LDA, 2, 0 );
  qrp.factorize( "qrp", mat );

  fmt::print( "QRP solution of A x = b\n" );
  alglin::Copy_n( rhs, M, x );
  alglin::Copy_n( rhs, M, b );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  fmt::print(
    "x^T = {}",
    alglin::print_matrix( 1, M, x, 1 )
  );

  alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  fmt::print(
    "residual = {}\ndone test7\n",
    alglin::print_matrix( 1, M, b, 1 )
  );
}


int
main() {

  try {
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
  } catch ( exception const & exc ) {
    msg.red( exc.what() );
  } catch ( ... ) {
    msg.red( "Errore Sconosciuto!\n" );
  }

  msg.green( "\n\nAll done!\n" );

  return 0;
}
