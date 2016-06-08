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
#include <vector>
#include <random>
#include "Alglin.hh"
#include "Alglin_aux.hh"
#include "TimeMeter.hh"
#include "LU_BABD_QR_N.hh"

using namespace std ;
typedef double valueType ;

unsigned seed1 = 2 ;
std::mt19937 generator(seed1);

valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

int
main() {

  #include "LU_test.hxx"

  alglin::BabdQR_N<valueType,NSIZE> LU ;

  //alglin::babd_print<valueType>( cout, nblk, n, q, AdAu, H0, HN, Hq ) ;

  TimeMeter tm ;
  tm.reset() ;

  tm.start() ;
  LU.factorize( nblk, q, AdAu, H0, HN, Hq ) ;
  tm.stop() ;
  cout << "Factorize (QRN) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  tm.start() ;
  //for ( int k = 0 ; k < 10 ; ++k ) {
    std::copy( rhs, rhs+N, x ) ;
    LU.solve( x ) ;
  //}
  tm.stop() ;
  cout << "Solve (QRN) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  //for ( alglin::integer i = 0 ; i < N ; ++i )
  //  cout << "x[" << i << "] = " << x[i] << '\n' ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  alglin::babd_residue<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                                   rhs, 1, x, 1, resid, 1 ) ;

  cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;

  LU.solve( resid ) ;
  alglin::axpy( N, +1.0, resid, 1, x, 1 ) ;

  cout << "All done!\n" ;

  return 0 ;
}
