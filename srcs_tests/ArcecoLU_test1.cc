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
  valueType block0[3*5] = {
     1,  1, 3, 1, 1,
     2,  2,  6,  8, 0,
     -2,  4,  -1,  8, 2
  } ;
  valueType block011[8] = {
    1, 2, 6, 8, -2, 4, 8, 2
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

  alglin::integer numBlock = 0 ; // 5-2 ;
  alglin::integer dim      = 4 ;
  alglin::integer row0     = 3 ;
  alglin::integer col0     = 5 ;
  alglin::integer rowN     = 4 ;
  alglin::integer colN     = 6 ;
  
  alglin::ColrowLU<valueType> LU ;
  
  alglin::integer N = row0 + rowN + numBlock*dim ;
  cout << "N = " << N << '\n' ;
  
  valueType x[20], y[20], resid[20] ;

  alglin::ColrowLU<valueType>::print( cout,
                                      row0, col0, block0,
                                      numBlock, dim, blocks,
                                      rowN, colN, blockN ) ;
  LU.factorize( row0, col0, block0,
                numBlock, dim, blocks,
                rowN, colN, blockN ) ;
  cout << "\n\n\n" ;
  LU.print(cout) ;
  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  //alglin::copy( N, x, 1, xref, 1 ) ;

  LU.mv( row0, col0, block0,
         numBlock, dim, blocks,
         rowN, colN, blockN,
         1.0, x, 1, 0, y, 1 ) ;

  for ( int i = 0 ; i < N ; ++i )
    cout << "rhs[" << i << "] = " << y[i] << "\n" ;

  LU.residue( row0, col0, block0,
              numBlock, dim, blocks,
              rowN, colN, blockN,
              y, 1, x, 1, resid, 1 ) ;

  cout << "Check |resid|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;

  LU.solve( y ) ;
  for ( int i = 0 ; i < 20 ; ++i )
    cout << "y[" << i << "] = " << y[i] << "\n" ;
  
  cout << "All done!\n" ;

  return 0 ;
}
