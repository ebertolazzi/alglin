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
  alglin::integer row0     = 5 ;
  alglin::integer col0     = dim+2 ;
  alglin::integer colN     = dim+10 ;
  alglin::integer rowN     = (dim-row0)+(col0+colN-2*dim) ;
  alglin::integer numBlock = 20 ;

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
      if ( i == j ) block0[i+j*row0] = 20+rand(-1,1) ;
      else          block0[i+j*row0] = rand(-1,1) ;
  
  for ( int i = 0 ; i < rowN ; ++i )
    for ( int j = 0 ; j < colN ; ++j )
      if ( i == j ) blockN[i+j*rowN] = 20+rand(-1,1) ;
      else          blockN[i+j*rowN] = rand(-1,1) ;

  for ( int k = 0 ; k < numBlock ; ++k )
    for ( int i = 0 ; i < dim ; ++i )
      for ( int j = 0 ; j < 2*dim ; ++j )
        if ( i == j ) blocks[2*k*dim*dim+i+j*dim] = 20+rand(-1,1) ;
        else          blocks[2*k*dim*dim+i+j*dim] = rand(-1,1) ;

  alglin::ColrowLU<valueType> LU ;

  cout << "N = " << N << '\n' ;

  //alglin::print_colrow( cout, numBlock, dim, row0, rowN, &block0.front(), &blocks.front(), &blockN.front() ) ;
  //alglin::print_colrow_to_maple( cout, numBlock, dim, row0, rowN, &block0.front(), &blocks.front(), &blockN.front() ) ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  std::copy( x, x+N, xref ) ;
  LU.mv( row0, col0, block0,
         numBlock, dim, blocks,
         rowN, colN, blockN,
         1.0, x, 1, 0, rhs, 1 ) ;

  std::copy( rhs, rhs+N, x ) ;

  TimeMeter tm ;
  tm.reset() ;
  tm.start() ;
  LU.factorize( row0, col0, block0,
                numBlock, dim, blocks,
                rowN, colN, blockN ) ;
  tm.stop() ;
  cout << "Factorize = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;
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

  int numInitialBc  = row0-(col0-dim) ;
  int numFinalBc    = dim-row0 ;
  int numInitialETA = 0 ;
  int numFinalETA   = 0 ;

  alglin::ArcecoLU LU1 ;
  LU1.load( numInitialBc,
            numFinalBc,
            numInitialETA,
            numFinalETA,
            numBlock,
            blocks, // AdAu,
            block0, // H0,
            blockN, // HN,
            blockN ) ;

  tm.start() ;
  LU1.factorize() ;
  tm.stop() ;
  cout << "Factorize (ARCECO) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  tm.start() ;
  LU1.solve(rhs) ;
  tm.stop() ;
  cout << "SOLVE (ARCECO) = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;

  //for ( int i = 0 ; i < N ; ++i )
  //  cout << "x[" << i << "] = " << x[i] << "\n" ;
  
  cout << "All done!\n" ;

  return 0 ;
}
