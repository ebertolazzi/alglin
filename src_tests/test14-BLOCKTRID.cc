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


int
main() {

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

  cout << "All done!\n";
  return 0;
}
