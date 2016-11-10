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
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"
#include "BABD_Amodio.hh"

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

int
main() {

  #include "LU_test.hxx"

  cout << "nblk = " << nblk << "\n"
       << "n    = " << n << "\n"
       << "q    = " << q << "\n" ;

  alglin::AmodioLU<valueType> LU ;

  //alglin::babd_print<valueType>( cout, nblk, n, q, AdAu, H0, HN, Hq ) ;

  alglin::LASTBLOCK_Choice ch[4] = { alglin::LASTBLOCK_LU,
                                     alglin::LASTBLOCK_QR,
                                     alglin::LASTBLOCK_QRP,
                                     alglin::LASTBLOCK_SVD } ;
  char const * kind[] = { "LU", "QR", "QRP", "SVD" } ;

  for ( int test = 0 ; test < 3 ; ++test ) {
    cout << "\n\n\ntest N." << test << "\n" ;
    TicToc tm ;
    tm.reset() ;

    tm.tic() ;
    LU.allocate( nblk, n, q, 0 );
    LU.loadBlocks( AdAu, n ) ;
    LU.loadBottom( H0, n+q, HN, n+q, Hq, n+q ) ;
    LU.selectLastBlockSolver( ch[test] ) ;
    LU.factorize() ;

    tm.toc() ;
    cout << "Factorize (Amodio" << kind[test] << ") = " << tm.elapsedMilliseconds() << " [ms]\n" ;

    tm.tic() ;
    //for ( int k = 0 ; k < 10 ; ++k ) {
      std::copy( rhs, rhs+N, x ) ;
      LU.solve( 10, x, N ) ;
    //}
    tm.toc() ;
    cout << "Solve (Amodio" << kind[test] << ") = " << tm.elapsedMilliseconds() << " [ms]\n" ;

    //for ( alglin::integer i = 0 ; i < N ; ++i )
    //  cout << "x[" << i << "] = " << x[i] << '\n' ;

    alglin::copy( N, xref, 1, xref1, 1 ) ;
    alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
    cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

    alglin::babd_residue<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                                     rhs, 1, x, 1, resid, 1 ) ; 

    cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;

    //LU.solve( resid ) ;
    //alglin::axpy( N, +1.0, resid, 1, x, 1 ) ;
  }

  cout << "All done!\n" ;

  return 0 ;
}
