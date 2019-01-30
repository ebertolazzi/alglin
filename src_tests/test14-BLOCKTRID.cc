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
#include <iostream>

using namespace std;


static
void
test1() {

  alglin::BlockTridiagonalSymmetic<alglin::doublereal> BT;
  alglin::integer rBlocks[] = { 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  alglin::doublereal D0[] = { 2, 1, 1, 1};
  alglin::doublereal D1[] = { 2, 0, 1,
                              0, 0, 0,
                              1, 0, 2};
  alglin::doublereal D2[] = { 2, 0,
                              0, 1};

  alglin::doublereal L0[] = { 1, 1, 0,
                              0, 1, 1 };
  alglin::doublereal L1[] = { 1, 1, 1,
                              1, 1, 1 };

  alglin::doublereal rhs[] = {
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.setD( 0, D0, 2 );
  BT.setD( 1, D1, 3 );
  BT.setD( 2, D2, 2 );
  BT.setL( 0, L0, 3 );
  BT.setL( 1, L1, 2 );

  BT.factorize();

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( alglin::integer k = 0; k < rBlocks[3]; ++k )
    cout << "x[ " << k << "] = " << rhs[7+k] << "\n";
}

static
void
test2() {

  alglin::BlockTridiagonalSymmetic<alglin::doublereal> BT;
  alglin::integer rBlocks[] = { 0, 2, 5, 7};
  BT.setup( 3, rBlocks );

  alglin::integer ii[] = {
    1, 2, 1, 2,
    3, 5, 3, 5,
    6, 7,
    3, 4, 4, 5,
    6, 6, 6, 7, 7, 7
  };

  alglin::integer jj[] = {
    1, 1, 2, 2,
    3, 3, 5, 5,
    6, 7,
    1, 1, 2, 2,
    3, 4, 5, 3, 4, 5
  };

  alglin::doublereal vals[] = {
    2, 1, 1, 1,
    2, 1, 1, 2,
    2, 1,
    1, 1, 1, 1,
    1, 1, 1, 1, 1, 1
  };

  alglin::doublereal rhs[] = {
    5, 4, 6, 4, 6, 5, 4,
    5, 4, 6, 4, 6, 5, 4
  };

  BT.zero();
  for ( int k = 0; k < 20; ++k )
    BT(ii[k]-1,jj[k]-1) += vals[k];

  BT.factorize();

  //BT.solve( rhs );
  BT.solve( 2, rhs, 7 );

  for ( alglin::integer k = 0; k < rBlocks[3]; ++k )
    cout << "x[ " << k << "] = " << rhs[7+k] << "\n";
}


int
main() {
  cout << "test1\n";
  test1();
  cout << "\n\ntest2\n";
  test2();

  cout << "All done!\n";
  return 0;
}
