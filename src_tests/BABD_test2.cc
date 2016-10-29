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

  alglin::integer dim      = 100/10 ;
  alglin::integer row0     = 40/10 ;
  alglin::integer col0     = dim+2 ;
  alglin::integer colN     = dim+10 ;
  alglin::integer rowN     = (dim-row0)+(col0+colN-2*dim) ;
  alglin::integer numBlock = 100 ;

  alglin::integer N   = row0 + rowN + numBlock*dim ;
  alglin::integer nnz = row0*col0 + rowN*colN + 2*dim*dim*numBlock + 5*N ;

  alglin::Malloc<valueType>       baseValue("real") ;
  alglin::Malloc<alglin::integer> baseIndex("integer") ;
  
  baseValue.allocate(size_t(nnz)) ;
  baseIndex.allocate(size_t(N)) ;

  valueType * block0 = baseValue(size_t(row0*col0)) ;
  valueType * blocks = baseValue(size_t(2*dim*dim*numBlock)) ;
  valueType * blockN = baseValue(size_t(rowN*colN)) ;
  valueType * x      = baseValue(size_t(N)) ;
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

  alglin::BABD<valueType> LU ;

  alglin::LASTBLOCK_Choice ch[4] = { alglin::LASTBLOCK_LU,
                                     alglin::LASTBLOCK_QR,
                                     alglin::LASTBLOCK_QRP,
                                     alglin::LASTBLOCK_SVD } ;

  alglin::BABD_Choice ch_sol[4] = { alglin::BABD_DIAZ,
                                    alglin::BABD_AMODIO,
                                    alglin::BABD_CYCLIC_REDUCTION_QR,
                                    alglin::BABD_CYCLIC_REDUCTION_QRP } ;

  for ( int solver = 0 ; solver < 4 ; ++solver ) {
    cout << "\n\n\nSOLVER: " << BABD_Choice_to_string(ch_sol[solver]) << '\n' ;
    LU.selectSolver(ch_sol[solver]) ;

    for ( int test = 0 ; test < 4 ; ++test ) {
      cout << "\n\nLAST_BLOCK: " << LastBlock_to_string(ch[test]) << '\n' ;

      LU.allocate( numBlock, dim );
      LU.loadBlocks( blocks, dim ) ;
      LU.loadTopBottom( row0, col0, block0, row0,
                        rowN, colN, blockN, rowN ) ;
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

    }
  }

  cout << "All done!\n" ;

  return 0 ;
}
