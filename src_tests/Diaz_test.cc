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
#include "TicToc.hh"
#include "ABD_Arceco.hh"
#include "ABD_Diaz.hh"

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
// std::chrono::system_clock::now().time_since_epoch().count();

static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

int
main() {

  alglin::integer dim      = 100/10 ;
  alglin::integer row0     = 40/10 ;
  alglin::integer col0     = dim+2 ;
  alglin::integer colN     = dim+10 ;
  alglin::integer rowN     = (dim-row0)+(col0+colN-2*dim) ;
  alglin::integer numBlock = 100 ;

  alglin::integer N   = row0 + rowN + numBlock*dim ;
  alglin::integer nnz = row0*col0 + rowN*colN + 2*dim*dim*numBlock + 14*N ;

  alglin::Malloc<valueType>       baseValue("real") ;
  alglin::Malloc<alglin::integer> baseIndex("integer") ;
  
  baseValue.allocate(size_t(nnz)) ;
  baseIndex.allocate(size_t(N)) ;

  valueType * block0 = baseValue(size_t(row0*col0)) ;
  valueType * blocks = baseValue(size_t(2*dim*dim*numBlock)) ;
  valueType * blockN = baseValue(size_t(rowN*colN)) ;
  valueType * x      = baseValue(size_t(10*N)) ;
  valueType * xref   = baseValue(size_t(N)) ;
  valueType * xref1  = baseValue(size_t(N)) ;
  valueType * rhs    = baseValue(size_t(N)) ;
  valueType * resid  = baseValue(size_t(N)) ;
  
  for ( int i = 0 ; i < row0 ; ++i )
    for ( int j = 0 ; j < col0 ; ++j )
      block0[i+j*row0] = rand(-1,0) ;
  
  for ( int i = 0 ; i < rowN ; ++i )
    for ( int j = 0 ; j < colN ; ++j )
      blockN[i+j*rowN] = rand(-1,0) ;

  for ( int k = 0 ; k < numBlock ; ++k )
    for ( int i = 0 ; i < dim ; ++i )
      for ( int j = 0 ; j < 2*dim ; ++j )
        blocks[2*k*dim*dim+i+j*dim] = rand(-1,0) ;
  
  // forzo diagonale dominanza
  valueType diag = 2*dim ;
  alglin::integer col00 = col0-dim ;
  alglin::integer row00 = row0-col00 ;
  for ( int i = 0 ; i < row0 ; ++i )
    block0[i*(row0+1)] += diag ;

  for ( int k = 0 ; k < numBlock ; ++k )
    for ( int i = 0 ; i < dim ; ++i )
      blocks[2*k*dim*dim+i+i*dim+row00*dim] += diag ;

  for ( int i = 0 ; i < rowN ; ++i )
    blockN[i*(rowN+1)+colN-rowN] += diag ;

  alglin::DiazLU<valueType>   LU ;
  alglin::ArcecoLU<valueType> LU_arceco ;

  cout << "N = " << N << '\n' ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i ;
  std::copy( x, x+N, xref ) ;
  alglin::abd_mv<valueType>( row0, col0, block0,
                             numBlock, dim, blocks,
                             rowN, colN, blockN,
                             1.0, x, 1, 0, rhs, 1 ) ;

  // ruoto per mettere le BC in fondo
  alglin::integer neq = numBlock*dim+row0+rowN ;
  std::rotate( rhs, rhs + row0, rhs + neq ) ;

  // ruoto la soluzione se sborda
  std::rotate( xref, xref + col0-dim, xref + neq ) ;

  alglin::LASTBLOCK_Choice ch[4] = { alglin::LASTBLOCK_LU,
                                     alglin::LASTBLOCK_QR,
                                     alglin::LASTBLOCK_QRP,
                                     alglin::LASTBLOCK_SVD } ;
  char const * kind[] = { "LU", "QR", "QRP", "SVD" } ;

  TicToc tm ;
  tm.reset() ;

  for ( int test = 0 ; test < 3 ; ++test ) {
    cout << "\n\n\ntest N." << test << "\n" ;

    tm.tic() ;
    LU.allocate( numBlock, dim, row0+rowN-dim, 0 );
    LU.loadBlocks( blocks, dim ) ;
    LU.loadTopBottom( row0, col0, block0, row0,
                      rowN, colN, blockN, rowN ) ;
    LU.selectLastBlockSolver( ch[test] ) ;
    LU.factorize() ;
    tm.toc() ;
    cout << "(Diaz " << kind[test] << ") Factorize = " << tm.elapsedMilliseconds() << " [ms]\n" ;

    std::copy( rhs, rhs+N, x ) ;
    tm.tic() ;
    LU.solve( x ) ;
    tm.toc() ;
    cout << "(Diaz " << kind[test] << ") Solve = " << tm.elapsedMilliseconds() << " [ms]\n" ;

    std::copy( rhs, rhs+N, x ) ;
    tm.tic() ;
    LU.solve( 10, x, N ) ;
    tm.toc() ;
    cout << "(Diaz " << kind[test] << ") Solve = " << tm.elapsedMilliseconds() << " [ms]\n" ;

    alglin::copy( N, xref, 1, xref1, 1 ) ;
    alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
    cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

    LU.solve( resid ) ;
    alglin::axpy( N, +1.0, resid, 1, x, 1 ) ;

    alglin::copy( N, xref, 1, xref1, 1 ) ;
    alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
    cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;
  }

  // rimetto rhs e x a posto
  std::rotate( rhs, rhs + neq - row0, rhs + neq ) ;
  std::rotate( xref, xref + neq - col0 + dim, xref + neq ) ;

  tm.tic() ;
  LU_arceco.factorize( row0, col0, block0,
                       numBlock, dim, blocks,
                       rowN, colN, blockN ) ;
  tm.toc() ;
  cout << "(arceco) Factorize = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  std::copy( rhs, rhs+N, x ) ;
  tm.tic() ;
  LU_arceco.solve( x ) ;
  tm.toc() ;
  cout << "\n\n(arceco) Solve = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "(arceco) Check|err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  alglin::abd_residue<valueType>( row0, col0, block0,
                                  numBlock, dim, blocks,
                                  rowN, colN, blockN,
                                  rhs, 1, x, 1, resid, 1 ) ;

  cout << "(arceco) Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;


  cout << "All done!\n" ;

  return 0 ;
}
