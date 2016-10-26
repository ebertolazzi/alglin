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
#include "CyclicReductionQR.hh"

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

  valueType LB[(2*n+q)*(2*n+q)], tmp[2*n+q] ;

  alglin::CyclicReductionQR<alglin::QRP<valueType> > CR ;

  TicToc tm ;
  tm.reset() ;
  tm.tic() ;
  //CR( nblk, n ) ;
  CR.reduce( nblk, n, AdAu ) ;

  // q, AdAu, H0, HN, Hq ) ;
  tm.toc() ;
  cout << "Reduction = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  tm.tic() ;
  alglin::integer nnq = 2*n+q ;
  alglin::gezero(nnq,nnq,LB,nnq) ;
  CR.getBlock_LR(LB,nnq) ;
  //alglin::print_matrix( cout, n, nnq, LB, nnq ) ;
  alglin::gecopy( nq, n, H0, nq, LB+n, nnq ) ;
  alglin::gecopy( nq, n, HN, nq, LB+nnq*n+n, nnq ) ;
  alglin::gecopy( nq, q, Hq, nq, LB+nnq*2*n+n, nnq ) ;
  //alglin::print_matrix( cout, nnq, nnq, LB, nnq ) ;
  //alglin::print_matrix( cout, nq, n, H0, nq ) ;
  //alglin::print_matrix( cout, nq, n, HN, nq ) ;
  //alglin::print_matrix( cout, nq, q, Hq, nq ) ;

  alglin::LU<valueType> lu ;
  lu.factorize( nnq, nnq, LB, nnq ) ;

  std::copy( rhs, rhs+N, x ) ;
  CR.forward( x ) ;
  alglin::copy( n,   x,        1, tmp,   1 ) ;
  alglin::copy( n+q, x+nblk*n, 1, tmp+n, 1 ) ;
  lu.solve(tmp) ;
  alglin::copy( n,   tmp,   1, x,        1 ) ;
  alglin::copy( n+q, tmp+n, 1, x+nblk*n, 1 ) ;
  CR.backward( x ) ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  alglin::babd_residue<valueType>( nblk, n, q, AdAu, H0, HN, Hq,
                                   rhs, 1, x, 1, resid, 1 ) ;

  cout << "Check |r|_inf = " << alglin::absmax( N, resid, 1 ) << '\n' ;
  cout << "All done!\n" ;

  return 0 ;
}