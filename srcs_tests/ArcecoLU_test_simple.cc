/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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
#include "alglin.hh"
#include "ColrowLU.hh"

/*
  1 -2
  1 1 1 2
  2 3 1 1
      1 1 1 2
      2 3 1 1
          1 1 1
          1 2 3
*/

using namespace std ;
typedef double valueType ;

int
main() {
#if 1
  valueType block0[2] = { 1, -2 } ;
  valueType blocks[16] = {
    1, 2,
    1, 3,
    1, 1,
    2, 4,
    1, 2,
    1, 3,
    1, 1,
    2, -4,
  } ;
  valueType blockN[6] = {
    1, 2,
    1, 2,
    1, 3
  };
#else
  valueType block0[2] = { 1, 2 } ;
  valueType blocks[8] = {
    0, 0,
    1, 2,
    1, 1,
    2, 4,
  } ;
  valueType blockN[6] = {
    0, 0,
    1, 0,
    1, 3
  };
#endif
  /*
  
  N = 5
  Block 0
       1        1
  Block 1
       1        1        1        2
       2        3        1        4
  Block N
       0        0        0
       1        1        1
       1        2        3
  */

  alglin::integer numBlock = 2 ;
  alglin::integer dim      = 2 ;
  alglin::integer row0     = 1 ;
  alglin::integer rowN     = 2 ;
  
  alglin::ColrowLU<valueType> LU ;
  
  alglin::integer N = row0 + rowN + numBlock*dim ;
  cout << "N = " << N << '\n' ;

  valueType x[7] = {1,2,3,4,5,6,7}, y[7] ;

  alglin::print_colrow( cout, numBlock, dim, row0, rowN, block0, blocks, blockN ) ;
  alglin::print_colrow_to_maple( cout, numBlock, dim, row0, rowN, block0, blocks, blockN ) ;

  LU.factorize( numBlock, dim, row0, rowN, block0, blocks, blockN ) ;
  LU.print(cout) ;
  LU.mv( numBlock, dim, row0, rowN, block0, blocks, blockN, 1.0, x, 1, 0, y, 1 ) ;
  for ( int i = 0 ; i < N ; ++i )
    cout << "rhs[" << i << "] = " << y[i] << "\n" ;
  LU.solve( y ) ;
  for ( int i = 0 ; i < N ; ++i )
    cout << "y[" << i << "] = " << y[i] << "\n" ;
  
  cout << "All done!\n" ;

  return 0 ;
}
