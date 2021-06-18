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

using namespace std;
typedef double real_type;
using lapack_wrapper::integer;

static
void
test1() {
  lapack_wrapper::QR<real_type> qr;
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

  cout
    << rang::fg::green << "\n\n\nTest1:\n\nInitial A\n"
    << lapack_wrapper::print_matrix( M, N, A, M )
    << rang::fg::reset;

  cout << "\nDo QR factorization of A^T\n";
  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  cout
    << "\nR=\n"
    << lapack_wrapper::print_matrix( M, M, R, M );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\nLS solution of A x = b\n"
    << "b^T      = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 );

  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 )
    << "done test1\n";
}



static
void
test2() {
  lapack_wrapper::QRP<real_type> qr;
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

  cout
    << "\n\n\nTest2:\n\nInitial A\n"
    << lapack_wrapper::print_matrix( M, N, A, M )
    << "\nDo QR factorization of A^T\n";
  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  cout
    << "\nR=\n"
    << lapack_wrapper::print_matrix( M, M, R, M );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\nLS solution of A x = b\n\n"
    << "b^T =      "
    << lapack_wrapper::print_matrix( 1, M, b, 1 );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( M,   rhs, 1, x, 1 );
  lapack_wrapper::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  cout
    << "x^T =      "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 )
    << "done test2\n";
}

static
void
test3() {
  lapack_wrapper::QRP<real_type> qr;
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

  cout
    << "\n\n\nTest3:\n\nInitial A\n"
    << lapack_wrapper::print_matrix( M, N, A, M );

  cout << "\nDo QR factorization of A^T\n";
  qr.t_factorize( "qr", M, N, A, LDA );

  real_type R[M*M];
  qr.getR( R, M );
  cout
    << "\nR=\n"
    << lapack_wrapper::print_matrix( M, M, R, M );

  real_type rhs[M], b[M];
  real_type x[N] = {1,2,3,4,5};
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\nLS solution of A x = b\n\n"
    << "b^T =      "
    << lapack_wrapper::print_matrix( 1, M, b, 1 );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  lapack_wrapper::copy( 3, rhs, 1, x, 1 );
  lapack_wrapper::zero( 2, x+3, 1 );
  qr.Q_mul( x );

  cout
    << "x^T =      "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1
  );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 )
    << "done test3\n";
}

#define TEST4(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  lapack_wrapper::copy( M, rhs, 1, x, 1 ); \
  lapack_wrapper::copy( M, rhs, 1, b, 1 ); \
  /* L.solve( x ); */ \
  F.solve( 1, x, M); \
  cout << "x^T      =" << lapack_wrapper::print_matrix( 1, M, x, 1 ); \
  \
  lapack_wrapper::gemv( lapack_wrapper::NO_TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual =" << lapack_wrapper::print_matrix( 1, M, b, 1 ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );

static
void
test4() {
  lapack_wrapper::LU<real_type>   lu;
  lapack_wrapper::LUPQ<real_type> lupq;
  lapack_wrapper::QR<real_type>   qr;
  lapack_wrapper::QRP<real_type>  qrp;
  lapack_wrapper::SVD<real_type>  svd;
  lapack_wrapper::LSS<real_type>  lss;
  lapack_wrapper::LSY<real_type>  lsy;
  lapack_wrapper::PINV<real_type> pinv;

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
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\n\n\nTest4:\n\nInitial A\n"
   << lapack_wrapper::print_matrix( M, M, A, M );

  TEST4("LU",lu);
  TEST4("LUPQ",lupq);
  TEST4("QR",qr);
  TEST4("QRP",qrp);
  TEST4("SVD",svd);
  TEST4("LSS",lss);
  TEST4("LSY",lsy);
  TEST4("PINV",pinv);

  cout << "done test4\n";
}

#define TEST5(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( NAME, M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b\n"; \
  lapack_wrapper::copy( M, rhs, 1, x, 1 ); \
  lapack_wrapper::copy( M, rhs, 1, b, 1 ); \
  /* F.t_solve( x ); */ \
  F.t_solve( 1, x, M); \
  cout << "x^T      = " << lapack_wrapper::print_matrix( 1, M, x, 1 ); \
  \
  lapack_wrapper::gemv( lapack_wrapper::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual = " << lapack_wrapper::print_matrix( 1, M, b, 1 ); \
  res = lapack_wrapper::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  UTILS_ASSERT0( res < 1e-6, "test failed!\n" );


static
void
test5() {
  lapack_wrapper::LU<real_type>   lu;
  lapack_wrapper::LUPQ<real_type> lupq;
  lapack_wrapper::QR<real_type>   qr;
  lapack_wrapper::QRP<real_type>  qrp;
  lapack_wrapper::SVD<real_type>  svd;
  lapack_wrapper::LSS<real_type>  lss;
  lapack_wrapper::LSY<real_type>  lsy;
  lapack_wrapper::PINV<real_type> pinv;

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
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\n\n\nTest5:\n\nInitial A\n"
    << lapack_wrapper::print_matrix( M, M, A, M );

  TEST5("LU",lu);
  TEST5("LUPQ",lupq);
  TEST5("QR",qr);
  TEST5("QRP",qrp);
  TEST5("SVD",svd);
  TEST5("LSS",lss);
  TEST5("LSY",lsy);
  TEST5("PINV",pinv);

  cout << "\ndone test5\n";
}


static
void
test6() {
  lapack_wrapper::TridiagonalLU<real_type> lu;
  lapack_wrapper::TridiagonalQR<real_type> qr;

  integer const N = 5;
  real_type const D[] = { 1, 1, 2, -0.1, 0.1 };
  real_type const L[] = { -0.1, -1, -2, -0.1 };
  real_type const U[] = { -1, -10, -2, 0.1 };

  real_type rhs[N], b[N];
  real_type x[N] = {1,2,3,4,5};
  qr.axpy( N, 1.0, L, D, U, x, 0.0, rhs );
  qr.axpy( N, 1.0, L, D, U, x, 0.0, b   );

  cout << "\n\n\nTest6:\n\nInitial A\n";

  cout << "\n\nDo (trid) LU  factorization of A\n";
  lu.factorize( "lu", N, L, D, U );

  cout << "(trid) LU solution of A x = b\n";
  lapack_wrapper::copy( N, rhs, 1, x, 1 );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  lu.solve( x );
  //qr.t_solve( 1, x, M);
  cout
    << "x^T =      "
    << lapack_wrapper::print_matrix( 1, N, x, 1 );

  lu.axpy( N, -1.0, L, D, U, x, 1.0, b );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, N, b, 1 )
    << "\n\nDo (trid) QR  factorization of A\n";
  qr.factorize( "qr", N, L, D, U );

  cout << "(trid) QR solution of A x = b\n";
  lapack_wrapper::copy( N, rhs, 1, x, 1 );
  lapack_wrapper::copy( N, rhs, 1, b, 1 );
  qr.solve( x );
  //qr.t_solve( 1, x, M);
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, N, x, 1 );

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, N, b, 1 )
    << "\ndone test6\n";
}

static
void
test7() {
  lapack_wrapper::QRP<real_type> qrp;

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
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1
  );
  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1
  );

  cout
    << "\n\n\nTest7:\n\nInitial A\n"
    << lapack_wrapper::print_matrix( M, M, A, M )
    << "\n\nDo QRP factorization of A\n";
  qrp.factorize( "qrp", M, M, A, LDA );

  cout << "QRP solution of A x = b\n";
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv(
    lapack_wrapper::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1
  );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 )
    << "\n\nDo QRP factorization of A\n";

  lapack_wrapper::Matrix<real_type> mat;
  mat.setup( M, M );
  mat.load_block( 2, 5, A,   LDA, 0, 0 );
  mat.load_block( 3, 5, A+2, LDA, 2, 0 );
  qrp.factorize( "qrp", mat );

  cout << "QRP solution of A x = b\n";
  lapack_wrapper::copy( M, rhs, 1, x, 1 );
  lapack_wrapper::copy( M, rhs, 1, b, 1 );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  cout
    << "x^T      = "
    << lapack_wrapper::print_matrix( 1, M, x, 1 );

  lapack_wrapper::gemv( lapack_wrapper::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  cout
    << "residual = "
    << lapack_wrapper::print_matrix( 1, M, b, 1 )
    << "\ndone test7\n";

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
    cerr << exc.what();
  } catch ( ... ) {
    cerr << "Errore Sconosciuto!\n";
  }

  cout << "\n\nAll done!\n";

  return 0;
}
