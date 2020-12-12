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


int
main() {

  alglin::doublereal s1[] = {1,2,3};
  alglin::doublereal y1[] = {11./3.,-6./3.,13./3.};

  alglin::doublereal s2[] = {1,0,1};
  alglin::doublereal y2[] = {3,-2,3};

  alglin::doublereal s3[] = {0,-1,1};
  alglin::doublereal y3[] = {7/3.,-6/3.,8/3.};

  alglin::BFGS<alglin::doublereal> bfgs;

  alglin::doublereal epsi = 1e-8;

  bfgs.allocate(3);
  bfgs.init();
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y1, s1, epsi );
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y2, s2, epsi );
  cout << "\n\n";
  bfgs.print( cout );

  bfgs.update( y3, s3, epsi );
  cout << "\n\n";
  bfgs.print( cout );

  for ( int i = 0; i < 90 ; ++i ) {
    bfgs.update( y1, s1, epsi );
    bfgs.update( y2, s2, epsi );
    bfgs.update( y3, s3, epsi );
  }
  cout << "\n\n";
  bfgs.print( cout );

  cout << "All done!\n";
  return 0;
}
