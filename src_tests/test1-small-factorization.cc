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
#include <vector>
#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

using namespace std;
typedef double valueType;
using alglin::integer;

static
void
test1() {
  alglin::QR<valueType> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  valueType A[] = {
    0.001,      2,     3,
    0.001,  0.001, 0,
    0,      0.001, 0,
    0.001,     -1,     0,
    0.000001,      5,     3
  };

  cout << "\n\n\nTest1:\n\nInitial A\n";
  alglin::print_matrix( cout, M, N, A, M );

  cout << "Do QR factorization of A^T\n";
  qr.t_factorize( M, N, A, LDA );
  
  valueType R[M*M];
  qr.getR( R, M );
  cout << "R=\n";
  alglin::print_matrix( cout, M, M, R, M );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "LS solution of A x = b, \nb=\n";
  alglin::print_matrix( cout, M, 1, b, M );

  qr.invRt_mul( rhs, 1 );
  alglin::copy( M, rhs, 1, x, 1 );
  alglin::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  cout << "x=\n";
  alglin::print_matrix( cout, 5, 1, x, 5 );

  alglin::gemv( alglin::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  cout << "residual=\n";
  alglin::print_matrix( cout, M, 1, b, M );
  cout << "done test1\n";
}

static
void
test2() {
  alglin::QRP<valueType> qr;
  integer const M   = 3;
  integer const N   = 5;
  integer const LDA = 3;
  valueType A[] = {
    0.001,      2,     3,
    0.001,  0.001, 0,
    0,      0.001, 0,
    0.001,     -1,     0,
    0.000001,      5,     3
  };

  cout << "\n\n\nTest2:\n\nInitial A\n";
  alglin::print_matrix( cout, M, N, A, M );

  cout << "Do QR factorization of A^T\n";
  qr.t_factorize( M, N, A, LDA );
  
  valueType R[M*M];
  qr.getR( R, M );
  cout << "R=\n";
  alglin::print_matrix( cout, M, M, R, M );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "LS solution of A x = b, \nb=\n";
  alglin::print_matrix( cout, M, 1, b, M );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  alglin::copy( M,   rhs, 1, x, 1 );
  alglin::zero( N-M, x+3, 1 );
  qr.Q_mul( x );

  cout << "x=\n";
  alglin::print_matrix( cout, 5, 1, x, 5 );

  alglin::gemv( alglin::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  cout << "residual=\n";
  alglin::print_matrix( cout, M, 1, b, M );
  cout << "done test2\n";
}

static
void
test3() {
  alglin::QRP<valueType> qr;
  integer const M   = 5;
  integer const N   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,     2, 3,
    0.001,  0.001,     0, 0.001, 1e-10,
    0,      0.001,     0, 0.001, 1e-12,
    0.001,     -1, 1e-12,    -1, -1e-12,
    0.000001,   5,     3,     5, 3
  };

  cout << "\n\n\nTest3:\n\nInitial A\n";
  alglin::print_matrix( cout, M, N, A, M );

  cout << "Do QR factorization of A^T\n";
  qr.t_factorize( M, N, A, LDA );
  
  valueType R[M*M];
  qr.getR( R, M );
  cout << "R=\n";
  alglin::print_matrix( cout, M, M, R, M );

  valueType rhs[M], b[M];
  valueType x[N] = {1,2,3,4,5};
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::NO_TRANSPOSE, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "LS solution of A x = b, \nb=\n";
  alglin::print_matrix( cout, M, 1, b, M );

  qr.inv_permute( rhs ); // da aggiungere!
  qr.invRt_mul( rhs, 1 );
  alglin::copy( 3, rhs, 1, x, 1 );
  alglin::zero( 2, x+3, 1 );
  qr.Q_mul( x );

  cout << "x=\n";
  alglin::print_matrix( cout, 5, 1, x, 5 );

  alglin::gemv( alglin::NO_TRANSPOSE, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  cout << "residual=\n";
  alglin::print_matrix( cout, M, 1, b, M );
  cout << "done test3\n";
}

#define TEST4(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b"; \
  alglin::copy( M, rhs, 1, x, 1 ); \
  alglin::copy( M, rhs, 1, b, 1 ); \
  /* L.solve( x ); */ \
  F.solve( 1, x, M); \
  cout << "x=\n"; \
  alglin::print_matrix( cout, M, 1, x, M ); \
  \
  alglin::gemv( alglin::NO_TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual=\n"; \
  alglin::print_matrix( cout, M, 1, b, M ); \
  res = alglin::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  ALGLIN_ASSERT( res < 1e-6, "test failed!" );

static
void
test4() {
  alglin::LU<valueType>   lu;
  alglin::LUPQ<valueType> lupq;
  alglin::QR<valueType>   qr;
  alglin::QRP<valueType>  qrp;
  alglin::SVD<valueType>  svd;
  alglin::LSS<valueType>  lss;
  alglin::LSY<valueType>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,      3+1
  };

  valueType rhs[M], b[M], res;
  valueType x[M] = {1,2,3,4,5};
  alglin::gemv( alglin::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::NO_TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "\n\n\nTest4:\n\nInitial A\n";
  alglin::print_matrix( cout, M, M, A, M );

  TEST4("LU",lu);
  TEST4("LUPQ",lupq);
  TEST4("QR",qr);
  TEST4("QRP",qrp);
  TEST4("SVD",svd);
  TEST4("LSS",lss);
  TEST4("LSY",lsy);

  cout << "done test4\n";
}

#define TEST5(NAME,F) \
  cout << "\n\nDo " << NAME << " factorization of A\n"; \
  F.factorize( M, M, A, LDA ); \
  \
  cout << NAME << " solution of A x = b"; \
  alglin::copy( M, rhs, 1, x, 1 ); \
  alglin::copy( M, rhs, 1, b, 1 ); \
  /* F.t_solve( x ); */ \
  F.t_solve( 1, x, M); \
  cout << "x=\n"; \
  alglin::print_matrix( cout, M, 1, x, M ); \
  \
  alglin::gemv( alglin::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 ); \
  cout << "residual=\n"; \
  alglin::print_matrix( cout, M, 1, b, M ); \
  res = alglin::nrm2( M, b, 1 ); \
  cout << "||res||_2 = " << res << '\n'; \
  ALGLIN_ASSERT( res < 1e-6, "test failed!" );


static
void
test5() {
  alglin::LU<valueType>   lu;
  alglin::LUPQ<valueType> lupq;
  alglin::QR<valueType>   qr;
  alglin::QRP<valueType>  qrp;
  alglin::SVD<valueType>  svd;
  alglin::LSS<valueType>  lss;
  alglin::LSY<valueType>  lsy;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,      3+1
  };

  valueType rhs[M], b[M], res;
  valueType x[M] = {1,2,3,4,5};
  alglin::gemv( alglin::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "\n\n\nTest5:\n\nInitial A\n";
  alglin::print_matrix( cout, M, M, A, M );

  TEST5("LU",lu);
  TEST5("LUPQ",lupq);
  TEST5("QR",qr);
  TEST5("QRP",qrp);
  TEST5("SVD",svd);
  TEST5("LSS",lss);
  TEST5("LSY",lsy);

  cout << "done test5\n";

}


static
void
test6() {
  alglin::TridiagonalLU<valueType> lu;
  alglin::TridiagonalQR<valueType> qr;

  integer const N = 5;
  valueType const D[] = { 1, 1, 2, -0.1, 0.1 };
  valueType const L[] = { -0.1, -1, -2, -0.1 };
  valueType const U[] = { -1, -10, -2, 0.1 };

  valueType rhs[N], b[N];
  valueType x[N] = {1,2,3,4,5};
  qr.axpy( N, 1.0, L, D, U, x, 0.0, rhs );
  qr.axpy( N, 1.0, L, D, U, x, 0.0, b   );

  cout << "\n\n\nTest6:\n\nInitial A\n";

  cout << "\n\nDo (trid) LU  factorization of A\n";
  lu.factorize( N, L, D, U );

  cout << "(trid) LU solution of A x = b\n";
  alglin::copy( N, rhs, 1, x, 1 );
  alglin::copy( N, rhs, 1, b, 1 );
  lu.solve( x );
  //qr.t_solve( 1, x, M);
  cout << "x=\n";
  alglin::print_matrix( cout, N, 1, x, N );

  lu.axpy( N, -1.0, L, D, U, x, 1.0, b );
  cout << "residual=\n";
  alglin::print_matrix( cout, N, 1, b, N );



  cout << "\n\nDo (trid) QR  factorization of A\n";
  qr.factorize( N, L, D, U );

  cout << "(trid) QR solution of A x = b\n";
  alglin::copy( N, rhs, 1, x, 1 );
  alglin::copy( N, rhs, 1, b, 1 );
  qr.solve( x );
  //qr.t_solve( 1, x, M);
  cout << "x=\n";
  alglin::print_matrix( cout, N, 1, x, N );

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );
  cout << "residual=\n";
  alglin::print_matrix( cout, N, 1, b, N );
  cout << "done test6\n";

}


static
void
test7() {
  alglin::QRP<valueType> qrp;

  integer const M   = 5;
  integer const LDA = 5;
  valueType A[] = {
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,      3+1
  };

  valueType rhs[M], b[M];
  valueType x[M] = {1,2,3,4,5};
  alglin::gemv( alglin::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( alglin::TRANSPOSE, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  cout << "\n\n\nTest7:\n\nInitial A\n";
  alglin::print_matrix( cout, M, M, A, M );

  cout << "\n\nDo QRP factorization of A\n";
  qrp.factorize( M, M, A, LDA );

  cout << "QRP solution of A x = b";
  alglin::copy( M, rhs, 1, x, 1 );
  alglin::copy( M, rhs, 1, b, 1 );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  cout << "x=\n";
  alglin::print_matrix( cout, 5, 1, x, 5 );

  alglin::gemv( alglin::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  cout << "residual=\n";
  alglin::print_matrix( cout, M, 1, b, M );



  cout << "\n\nDo QRP factorization of A\n";
  qrp.allocate( M, M );
  qrp.load_block( 2, 5, A,   LDA, 0, 0 );
  qrp.load_block( 3, 5, A+2, LDA, 2, 0 );
  qrp.factorize();
  
  cout << "QRP solution of A x = b";
  alglin::copy( M, rhs, 1, x, 1 );
  alglin::copy( M, rhs, 1, b, 1 );
  qrp.t_solve( x );
  //qrp.t_solve( 1, x, M );
  cout << "x=\n";
  alglin::print_matrix( cout, 5, 1, x, 5 );

  alglin::gemv( alglin::TRANSPOSE, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  cout << "residual=\n";
  alglin::print_matrix( cout, M, 1, b, M );
  cout << "done test7\n";

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
    cerr << exc.what() << '\n';
  } catch ( ... ) {
    cerr << "Errore Sconosciuto!\n";
  }

  cout << "\n\nAll done!\n";

  return 0;
}
