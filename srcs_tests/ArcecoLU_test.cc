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
*/

using namespace std ;
typedef double valueType ;

int
main() {
  valueType block0[2*4] = {
     1,  3,  5,  7,
     2,  4,  6,  8
  } ;
  valueType blocks[3*4*8] = {
     9,  13, 17, 21, 10, 14, 18, 2,
     10, 16, 18, 22, 10, 14, 121, 22,
     11, 15, 19, 23, 10, 14, 18, 22,
     12, 16, 20, 24, 10, 14, 15, 22,
     9,  13, 17, 21, 10, 14, 18, 22,
     10, 14, -3345, 22, 11, 54, 13, 22,
     11, 15, 19, 23, 10, 14, 28, 22,
     12, 16, 20, 24, 12, 10, 18, 22,
     9,  13, 17, 21, -2345, 14, 18, 22,
     10, 14, 18, 22, 13, 11, 18, 12,
     11, 15, 19, 23, 10, 13, 12, 22,
     12, 16, 20, 24, 10, 14, 11, 22,
  } ;
  valueType blockN[4*6] = {
    1, 99,  3,  7, 11, 1,
    0,  0,  4,  8, 12, 2,
    0,  1,  5,  9, 234, 3,
    -1,  2,  6, 10, 14, 4
  };

  alglin::integer numBlock = 5-2 ;
  alglin::integer dim      = 4 ;
  alglin::integer row0     = 2 ;
  alglin::integer rowN     = 4 ;
  
  alglin::ColrowLU<valueType> LU ;
  
  alglin::integer N = row0 + rowN + numBlock*4 ;
  cout << "N = " << N << '\n' ;
  
  valueType x[20], y[20] ;

  alglin::ColrowLU<valueType>::print_colrow( cout, numBlock, dim, row0, rowN, block0, blocks, blockN ) ;
  alglin::ColrowLU<valueType>::print_colrow_to_maple( cout, numBlock, dim, row0, rowN, block0, blocks, blockN ) ;

  LU.factorize( numBlock, dim, row0, rowN, block0, blocks, blockN ) ;
  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  LU.mv( numBlock, dim, row0, rowN, block0, blocks, blockN, 2.0, x, 1, 0, y, 1 ) ;
  LU.solve( y ) ;
  for ( int i = 0 ; i < 20 ; ++i )
    cout << "y[" << i << "] = " << y[i] << "\n" ;
  
  cout << "All done!\n" ;

  return 0 ;
}
