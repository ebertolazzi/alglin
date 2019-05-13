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

#include "Alglin++.hh"
#include "KKT_like.hh"
#include <iostream>

using namespace std;

typedef alglin::doublereal real_type;

static
void
test1() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A[] = { 1, 2, 3, 4,
                    4, 4, 4, 4,
                    1, 2, 1, 2,
                    0, 0, 1, 1 };
  alglin::Eigenvalues<real_type> E;
  E.setup( 4, A, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";
}

static
void
test2() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A[] = { 1, 2, 3, 4,
                    4, 4, 4, 4,
                    1, 2, 1, 2,
                    0, 0, 1, 1 };
  real_type B[] = { 3, 2, 1, 1,
                    1, 1, 1, 1,
                    4, 3, 2, 1,
                    1,-1, 0, 1 };
  alglin::GeneralizedEigenvalues<real_type> E;
  E.setup( 4, A, 4, B, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";
}

static
void
test3() {
  /*
  A = [ 1, 2, 3, 4; 4, 4, 4, 4; 1, 2, 1, 2; 0, 0, 1, 1].';
  B = [ 3, 2, 1, 1; 1, 1, 1, 1; 4, 3, 2, 1; 1,-1, 0, 1].';
  [U,V,X,C,S] = gsvd(A,B)
  %[U,V,X,C,S] = gsvd(A,B,0)
  */
  real_type A_data[] = { 1, 2, 3, 4,
                         4, 4, 4, 4,
                         1, 2, 1, 2,
                         0, 0, 1, 1 };
  real_type B_data[] = { 3, 2, 1, 1,
                         1, 1, 1, 1,
                         4, 3, 2, 1,
                         1,-1, 0, 1 };
  real_type TMP_data[16], TMP1_data[16];
  alglin::MatrixWrapper<real_type> A(A_data,4,4,4);
  alglin::MatrixWrapper<real_type> B(B_data,4,4,4);
  alglin::MatrixWrapper<real_type> TMP(TMP_data,4,4,4);
  alglin::MatrixWrapper<real_type> TMP1(TMP1_data,4,4,4);

  alglin::GeneralizedSVD<real_type> E;
  E.setup( A, B );
  E.info( cout, 1e-10 ) ;

  alglin::gemm( 1.0,
                alglin::TRANSPOSE,    E.getU(),
                alglin::NO_TRANSPOSE, A,
                0.0, TMP );
  alglin::gemm( 1.0,
                alglin::NO_TRANSPOSE, TMP,
                alglin::NO_TRANSPOSE, E.getQ(),
                0.0, TMP1 );
  cout << "TMP1\n";
  TMP1.print0(cout,1e-12);

  alglin::gemm( 1.0,
                alglin::TRANSPOSE,    E.getV(),
                alglin::NO_TRANSPOSE, B,
                0.0, TMP );
  alglin::gemm( 1.0,
                alglin::NO_TRANSPOSE, TMP,
                alglin::NO_TRANSPOSE, E.getQ(),
                0.0, TMP1 );
  cout << "TMP1\n";
  TMP1.print0(cout,1e-12);
}

static
void
test4() {
  /*
  */
  real_type A_data[] = { 1, 2, 3, 4, 5,
                         6, 7, 8, 9, 10,
                         11,12,13,14,15 };
  real_type B_data[] = { 8, 3, 4,
                         1, 5, 9,
                         6, 7, 2 };
  real_type TMP_data[100], TMP1_data[100], TMP2_data[100], TMP3_data[100];
  alglin::MatrixWrapper<real_type> A(A_data,5,3,5);
  alglin::MatrixWrapper<real_type> B(B_data,3,3,3);
  alglin::MatrixWrapper<real_type> TMP(TMP_data,5,3,5);
  alglin::MatrixWrapper<real_type> TMP1(TMP1_data,5,3,5);
  alglin::MatrixWrapper<real_type> TMP2(TMP2_data,3,3,3);
  alglin::MatrixWrapper<real_type> TMP3(TMP3_data,3,3,3);

  alglin::GeneralizedSVD<real_type> E;
  E.setup( A, B );
  E.info( cout, 1e-12) ;

  alglin::gemm( 1.0,
                alglin::TRANSPOSE,    E.getU(),
                alglin::NO_TRANSPOSE, A,
                0.0, TMP );
  alglin::gemm( 1.0,
                alglin::NO_TRANSPOSE, TMP,
                alglin::NO_TRANSPOSE, E.getQ(),
                0.0, TMP1 );
  cout << "TMP1\n";
  TMP1.print0(cout,1e-12);

  alglin::gemm( 1.0,
                alglin::TRANSPOSE,    E.getV(),
                alglin::NO_TRANSPOSE, B,
                0.0, TMP2 );
  alglin::gemm( 1.0,
                alglin::NO_TRANSPOSE, TMP2,
                alglin::NO_TRANSPOSE, E.getQ(),
                0.0, TMP3 );
  cout << "TMP3\n";
  TMP3.print0(cout,1e-12);
}

static
void
test5() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[] = {
    2, 0, 2, 0,
    0, 2, 4, 0,
    2, 4, 4, 1,
    0, 0, 0, 4
  };
  alglin::MatrixWrapper<real_type> A(A_data,4,4,4);
  alglin::Eigenvectors<real_type> E;
  E.setup( 4, A_data, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";

  vector<vector< alglin::Eigenvectors<real_type>::complexType > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    cout << "vL[" << n << "] = ";
    for ( size_t i = 0; i < 4; ++i )
      cout << ' ' << vecs[n][i];
    cout << '\n';
  }
  E.getRightEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    cout << "vR[" << n << "] = ";
    for ( size_t i = 0; i < 4; ++i )
      cout << ' ' << vecs[n][i];
    cout << '\n';
  }
}

static
void
test6() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[] = { 1, 2, 3, 4,
                         4, 4, 4, 4,
                         1, 2, 1, 2,
                         0, 0, 1, 1 };
  real_type B_data[] = { 3, 2, 1, 1,
                         1, 1, 1, 1,
                         4, 3, 2, 1,
                         1,-1, 0, 1 };
  alglin::MatrixWrapper<real_type> A(A_data,4,4,4);
  alglin::MatrixWrapper<real_type> B(B_data,4,4,4);
  alglin::GeneralizedEigenvectors<real_type> E;
  E.setup( 4, A_data, 4, B_data, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  vector<complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";

  vector<vector< alglin::GeneralizedEigenvectors<real_type>::complexType > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    cout << "vL[" << n << "] = ";
    for ( size_t i = 0; i < 4; ++i )
      cout << ' ' << vecs[n][i];
    cout << '\n';
  }
  E.getRightEigenvector( vecs );
  for ( size_t n = 0; n < 4; ++n ) {
    cout << "vR[" << n << "] = ";
    for ( size_t i = 0; i < 4; ++i )
      cout << ' ' << vecs[n][i];
    cout << '\n';
  }
}

int
main() {
  cout << "test1\n";
  test1();
  cout << "\n\ntest2\n";
  test2();
  cout << "\n\ntest3\n";
  test3();
  cout << "\n\ntest4\n";
  test4();
  cout << "\n\ntest5\n";
  test5();
  cout << "\n\ntest6\n";
  test6();
  return 0;
}