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
#include "alglin.hh"
#include "ColrowLU.hh"
#include "TimeMeter.hh"
#include "LU_Arceco.hh"

using namespace std ;
typedef double valueType ;

unsigned seed1 = 2 ;
// std::chrono::system_clock::now().time_since_epoch().count();

std::mt19937 generator(seed1);

valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

int
main() {

  alglin::integer dim      = 100 ;
  alglin::integer row0     = 40 ;
  alglin::integer col0     = dim+2 ;
  alglin::integer colN     = dim+10 ;
  alglin::integer rowN     = (dim-row0)+(col0+colN-2*dim) ;
  alglin::integer numBlock = 1000 ;

  alglin::integer N   = row0 + rowN + numBlock*dim ;
  alglin::integer nnz = row0*col0 + rowN*colN + 2*dim*dim*numBlock + 5*N ;

  alglin::Malloc<valueType>       baseValue("real") ;
  alglin::Malloc<alglin::integer> baseIndex("integer") ;
  
  baseValue.allocate(nnz) ;
  baseIndex.allocate(N) ;

  valueType * block0 = baseValue(row0*col0) ;
  valueType * blocks = baseValue(2*dim*dim*numBlock) ;
  valueType * blockN = baseValue(rowN*colN) ;
  valueType * x      = baseValue(N) ;
  valueType * xref   = baseValue(N) ;
  valueType * xref1  = baseValue(N) ;
  valueType * rhs    = baseValue(N) ;
  valueType * resid  = baseValue(N) ;
  
  for ( int i = 0 ; i < row0 ; ++i )
    for ( int j = 0 ; j < col0 ; ++j )
      block0[i+j*row0] = rand(-1,1) ;
  
  for ( int i = 0 ; i < rowN ; ++i )
    for ( int j = 0 ; j < colN ; ++j )
      blockN[i+j*rowN] = rand(-1,1) ;

  for ( int k = 0 ; k < numBlock ; ++k )
    for ( int i = 0 ; i < dim ; ++i )
      for ( int j = 0 ; j < 2*dim ; ++j )
        blocks[2*k*dim*dim+i+j*dim] = rand(-1,1) ;
  
  // forzo diagonale dominanza
  valueType diag = 8.5 ;
  alglin::integer col00 = col0-dim ;
  alglin::integer row00 = row0-col00 ;
  for ( int i = 0 ; i < row0 ; ++i )
    block0[i*(row0+1)] += diag ;

  for ( int k = 0 ; k < numBlock ; ++k )
    for ( int i = 0 ; i < dim ; ++i )
      blocks[2*k*dim*dim+i+i*dim+row00*dim] += diag ;

  for ( int i = 0 ; i < rowN ; ++i )
    blockN[i*(rowN+1)+colN-rowN] += diag ;

  alglin::ColrowLU<valueType> LU(false) ;
  alglin::ColrowLU<valueType> LU_arceco(true) ;

  cout << "N = " << N << '\n' ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  std::copy( x, x+N, xref ) ;
  LU.mv( row0, col0, block0,
         numBlock, dim, blocks,
         rowN, colN, blockN,
         1.0, x, 1, 0, rhs, 1 ) ;

  TimeMeter tm ;
  tm.reset() ;

  tm.start() ;
  LU.factorize( row0, col0, block0,
                numBlock, dim, blocks,
                rowN, colN, blockN ) ;
  tm.stop() ;
  cout << "Factorize = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  tm.start() ;
  LU_arceco.factorize( row0, col0, block0,
                       numBlock, dim, blocks,
                       rowN, colN, blockN ) ;
  tm.stop() ;
  cout << "Factorize (arceco) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  std::copy( rhs, rhs+N, x ) ;
  tm.start() ;
  LU.solve( x ) ;
  tm.stop() ;
  cout << "Solve = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  LU.residue( row0, col0, block0,
              numBlock, dim, blocks,
              rowN, colN, blockN,
              rhs, 1, x, 1, resid, 1 ) ;

  cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;
  
  LU.solve( resid ) ;
  alglin::axpy( N, +1.0, resid, 1, x, 1 ) ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  LU.residue( row0, col0, block0,
              numBlock, dim, blocks,
              rowN, colN, blockN,
              rhs, 1, x, 1, resid, 1 ) ;
  cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;



  std::copy( rhs, rhs+N, x ) ;
  tm.start() ;
  LU_arceco.solve( x ) ;
  tm.stop() ;
  cout << "\n\nSolve (arceco) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check (arceco) |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  LU.residue( row0, col0, block0,
              numBlock, dim, blocks,
              rowN, colN, blockN,
              rhs, 1, x, 1, resid, 1 ) ;

  cout << "Check (arceco) |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;


  cout << "All done!\n" ;

  return 0 ;
}
