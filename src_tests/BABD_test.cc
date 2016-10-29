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
#include "BABD.hh"

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
static std::mt19937 generator(seed1);

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

  alglin::BABD<valueType> LU ;

  alglin::LASTBLOCK_Choice ch[4] = { alglin::LASTBLOCK_LU,
                                     alglin::LASTBLOCK_QR,
                                     alglin::LASTBLOCK_QRP,
                                     alglin::LASTBLOCK_SVD } ;

  alglin::BABD_Choice ch_sol[3] = { alglin::BABD_AMODIO,
                                    alglin::BABD_CYCLIC_REDUCTION_QR,
                                    alglin::BABD_CYCLIC_REDUCTION_QRP } ;

  for ( int solver = 0 ; solver < 3 ; ++solver ) {
    cout << "\n\n\nSOLVER: " << BABD_Choice_to_string(ch_sol[solver]) << '\n' ;
    LU.selectSolver(ch_sol[solver]) ;

    for ( int test = 0 ; test < 4 ; ++test ) {
      cout << "\n\nLAST_BLOCK: " << LastBlock_to_string(ch[test]) << '\n' ;

      LU.allocate( nblk, n );
      LU.loadBlocks( AdAu, n ) ;
      LU.loadBottom( q, H0, n+q, HN, n+q, Hq, n+q ) ;
      LU.selectLastBlockSolver( ch[test] ) ;

      TicToc tm ;
      tm.reset() ;
      tm.tic() ;
      LU.factorize() ;
      tm.toc() ;
      cout << "Factorize = " << tm.elapsedMilliseconds() << " [ms]\n" ;

      tm.tic() ;
      std::copy( rhs, rhs+N, x ) ;
      LU.solve( x ) ;
  
      tm.toc() ;
      cout << "Solve = " << tm.elapsedMilliseconds() << " [ms]\n" ;

      alglin::copy( N, xref, 1, xref1, 1 ) ;
      alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
      cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

      alglin::babd_residue<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                                       rhs, 1, x, 1, resid, 1 ) ;

      cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;
    }
  }

  cout << "All done!\n" ;

  return 0 ;
}
