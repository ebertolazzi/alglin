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

namespace fmt {
  template <typename Real> struct formatter<std::complex<Real>> : ostream_formatter {};
}

using namespace std;

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

using alglin::real_type;
using alglin::Transposition;

static
void
test1() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  constexpr real_type A[]{ 1, 2, 3, 4,
                           4, 4, 4, 4,
                           1, 2, 1, 2,
                           0, 0, 1, 1 };
  alglin::Eigenvalues<real_type> E;
  E.setup( 4, A, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  for ( auto const & it : e ) fmt::print( "{}\n", it );
}

static
void
test2() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  constexpr real_type A[]{ 1, 2, 3, 4,
                           4, 4, 4, 4,
                           1, 2, 1, 2,
                           0, 0, 1, 1 };
  constexpr real_type B[]{ 3, 2, 1, 1,
                           1, 1, 1, 1,
                           4, 3, 2, 1,
                           1,-1, 0, 1 };
  alglin::GeneralizedEigenvalues<real_type> E;
  E.setup( 4, A, 4, B, 4 );
  vector<complex<real_type> > e;
  E.getEigenvalues( e );

  for ( auto const & it : e ) fmt::print( "{}\n", it );
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
  real_type A_data[]{ 1, 2, 3, 4,
                      4, 4, 4, 4,
                      1, 2, 1, 2,
                      0, 0, 1, 1 };
  real_type B_data[]{ 3, 2, 1, 1,
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

  fmt::print( "E={}\n", E.info(1e-10) );

  alglin::gemm( 1.0,
                Transposition::YES, E.getU(),
                Transposition::NO,  A,
                0.0, TMP );
  alglin::gemm( 1.0,
                Transposition::NO, TMP,
                Transposition::NO, E.getQ(),
                0.0, TMP1 );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );

  alglin::gemm( 1.0,
                Transposition::YES, E.getV(),
                Transposition::NO,  B,
                0.0, TMP );
  alglin::gemm( 1.0,
                Transposition::NO, TMP,
                Transposition::NO, E.getQ(),
                0.0, TMP1 );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );
}

static
void
test4() {

  using namespace alglin;

  /*
  */
  real_type A_data[]{ 1, 2, 3, 4, 5,
                      6, 7, 8, 9, 10,
                      11,12,13,14,15 };
  real_type B_data[]{ 8, 3, 4,
                      1, 5, 9,
                      6, 7, 2 };
  real_type TMP_data[100], TMP1_data[100], TMP2_data[100], TMP3_data[100];
  MatrixWrapper<real_type> A(A_data,5,3,5);
  MatrixWrapper<real_type> B(B_data,3,3,3);
  MatrixWrapper<real_type> TMP(TMP_data,5,3,5);
  MatrixWrapper<real_type> TMP1(TMP1_data,5,3,5);
  MatrixWrapper<real_type> TMP2(TMP2_data,3,3,3);
  MatrixWrapper<real_type> TMP3(TMP3_data,3,3,3);

  alglin::GeneralizedSVD<real_type> E;
  E.setup( A, B );
  fmt::print( "E={}\n", E.info(1e-10) );

  gemm( 1.0,
        Transposition::YES, E.getU(),
        Transposition::NO,  A,
        0.0, TMP );
  gemm( 1.0,
        Transposition::NO, TMP,
        Transposition::NO, E.getQ(),
        0.0, TMP1 );
  fmt::print( "TMP1\n{}", TMP1.to_string(1e-12) );

  gemm( 1.0,
        Transposition::YES, E.getV(),
        Transposition::NO,  B,
        0.0, TMP2 );
  gemm( 1.0,
        Transposition::NO, TMP2,
        Transposition::NO, E.getQ(),
        0.0, TMP3 );
  fmt::print( "TMP3\n{}", TMP3.to_string(1e-12) );
}

static
void
test5() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[]{
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

  vector<vector< alglin::Eigenvectors<real_type>::complex_type > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n{0}; n < 4; ++n ) {
    fmt::print("vL[{}] = ",n);
    for ( size_t i{0}; i < 4; ++i )
      fmt::print(" {}",vecs[n][i]);
    cout << '\n';
  }
  E.getRightEigenvector( vecs );
  for ( size_t n{0}; n < 4; ++n ) {
    fmt::print("vR[{}] = ",n);
    for ( size_t i{0}; i < 4; ++i )
      fmt::print(" {}",vecs[n][i]);
    cout << '\n';
  }
}

static
void
test6() {
  // solution
  // 0, 0, 7/2+(1/2)*sqrt(73), 7/2-(1/2)*sqrt(73)])
  real_type A_data[]{ 1, 2, 3, 4,
                      4, 4, 4, 4,
                      1, 2, 1, 2,
                       0, 0, 1, 1 };
  real_type B_data[]{ 3, 2, 1, 1,
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

  vector<vector< alglin::GeneralizedEigenvectors<real_type>::complex_type > > vecs;
  E.getLeftEigenvector( vecs );
  for ( size_t n{0}; n < 4; ++n ) {
    fmt::print("vL[{}] = ",n);
    for ( size_t i{0}; i < 4; ++i )
      fmt::print(" {}",vecs[n][i]);
    cout << '\n';
  }
  E.getRightEigenvector( vecs );
  for ( size_t n{0}; n < 4; ++n ) {
    fmt::print("vR[{}] = ",n);
    for ( size_t i{0}; i < 4; ++i )
      fmt::print(" {}",vecs[n][i]);
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
  msg.green( "All done!\n" );
  return 0;
}
