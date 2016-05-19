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
  alglin::integer row0     = 2  ;
  alglin::integer rowN     = dim-row0+2 ;
  alglin::integer numBlock = 50 ;

/*
  alglin::integer iscale   = 50 ;
  alglin::integer row0     = 2*iscale  ;
  alglin::integer rowN     = 10*iscale ;
  alglin::integer dim      = 4*iscale  ;
  alglin::integer numBlock = 1  ;
*/

  alglin::integer N   = row0 + rowN + numBlock*dim ;
  alglin::integer nnz = row0*dim +
                        2*dim*dim*numBlock +
                        (row0+rowN)*rowN +
                        3*N ;

  valueType * p_memory = new valueType[ nnz ] ;

  valueType * ptr    = p_memory ;
  valueType * block0 = ptr ; ptr += row0*dim ;
  valueType * blocks = ptr ; ptr += 2*dim*dim*numBlock ;
  valueType * blockN = ptr ; ptr += (row0+rowN)*rowN ;
  valueType * x      = ptr ; ptr += N ;
  valueType * xref   = ptr ; ptr += N ;
  valueType * rhs    = ptr ; ptr += N ;
  
  for ( int i = 0 ; i < nnz ; ++i ) p_memory[i] = rand(-1,1) ;
  for ( int j = 0 ; j < row0 ; ++j ) {
    block0[(row0+1)*j] += 20 ;
  }
  for ( int j = 0 ; j < rowN ; ++j ) {
    blockN[2*rowN+(rowN+1)*j] += 20 ;
  }
  //std::copy( blockN, blockN+rowN+row0, blockN+rowN+row0 ) ;
  //std::copy( blockN, blockN+rowN+row0, blockN+2*(rowN+row0) ) ;
  //std::fill( blockN, blockN+(rowN+row0)*(rowN-1), 0 ) ;
  
  for ( int i = 0 ; i < numBlock ; ++i ) {
    valueType * blks = blocks + 2*dim*dim*i ;
    for ( int j = 0 ; j < dim ; ++j ) {
      blks[(dim+1)*j]         += 10 ;
      blks[dim*dim+(dim+1)*j] -= 10 ;
    }
  }

  alglin::ColrowLU<valueType> LU ;

  cout << "N = " << N << '\n' ;

  //alglin::print_colrow( cout, numBlock, dim, row0, rowN, &block0.front(), &blocks.front(), &blockN.front() ) ;
  //alglin::print_colrow_to_maple( cout, numBlock, dim, row0, rowN, &block0.front(), &blocks.front(), &blockN.front() ) ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  std::copy( x, x+N, xref ) ;
  LU.mv( numBlock, dim, row0, rowN,
         block0, blocks, blockN,
         1.0, x, 1, 0, rhs, 1 ) ;
  std::copy( rhs, rhs+N, x ) ;

  TimeMeter tm ;
  tm.reset() ;
  tm.start() ;
  LU.factorize( numBlock, dim, row0, rowN, block0, blocks, blockN ) ;
  tm.stop() ;
  cout << "Factorize = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;
  tm.start() ;
  LU.solve( x ) ;
  tm.stop() ;
  cout << "Solve = " << tm.partialElapsedMilliseconds() << " [ms]\n" ;
  alglin::axpy( N, -1.0, x, 1, xref, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref, 1 ) << '\n' ;
  
  //for ( int i = 0 ; i < N ; ++i )
  //  cout << "x[" << i << "] = " << x[i] << "\n" ;
  
  cout << "All done!\n" ;
  
  delete [] p_memory ;

  return 0 ;
}
