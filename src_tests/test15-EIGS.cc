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
  alglin::Eigenvalues<double> E;
  E.setup( 4, A, 4 );
  std::vector<std::complex<real_type> > e;
  E.getEigenvalues( e );

  std::vector<std::complex<real_type> >::const_iterator it;
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
  alglin::GeneralizedEigenvalues<double> E;
  E.setup( 4, A, 4, B, 4 );
  std::vector<std::complex<real_type> > e;
  E.getEigenvalues( e );

  std::vector<std::complex<real_type> >::const_iterator it;
  for ( it = e.begin(); it != e.end(); ++it )
    cout << *it << "\n";
}

int
main() {
  cout << "test1\n";
  test1();
  cout << "\n\ntest2\n";
  test2();
  //cout << "\n\ntest3\n";
  //test3();
  //cout << "\n\ntest4\n";
  //test4();
  //cout << "All done!\n";
  return 0;
}
