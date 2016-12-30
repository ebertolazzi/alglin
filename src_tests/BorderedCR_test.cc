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
#include <fstream>
#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"
#include "BorderedCR.hh"

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

  alglin::BorderedCR<double> BCR ;

  #define NSIZE 10
  //#define NSIZE 5

  alglin::integer n      = NSIZE ;
  alglin::integer nblock = 200000 ;
  alglin::integer q      = 4 ;
  alglin::integer nb     = 2 ;
  alglin::integer nq     = n+q ;
  alglin::integer N      = (nblock+1)*n+nb+q ;
 
  BCR.allocate( nblock, n, q, nb ) ;

  alglin::Malloc<valueType>       baseValue("real") ;
  alglin::Malloc<alglin::integer> baseIndex("integer") ;
  
  baseValue.allocate( size_t(6*N) ) ;
  
  valueType diag = 1.01*n ;

  valueType * x     = baseValue(size_t(2*N)) ; // extra space per multiple rhs
  valueType * xref  = baseValue(size_t(N)) ;
  valueType * xref1 = baseValue(size_t(N)) ;
  valueType * rhs   = baseValue(size_t(N)) ;
  valueType * resid = baseValue(size_t(N)) ;
  
  for ( int i = 0 ; i < nq ; ++i ) {
    for ( int j = 0 ; j < (nq+n+nb) ; ++j ) {
      BCR.H(i,j) = rand(-1,0) ;
    }
    BCR.H(i,i+n) += diag ; // force diagonal dominance
  }

  for ( int k = 0 ; k < nblock ; ++k ) {
    for ( int i = 0 ; i < n ; ++i ) {
      for ( int j = 0 ; j < n ; ++j ) {
        BCR.D(k,i,j) = rand(-1,0) ;
        BCR.E(k,i,j) = rand(-1,0) ;
      }
      BCR.D(k,i,i) += diag ; // force diagonal dominance
    }
    for ( int i = 0 ; i < n ; ++i ) {
      for ( int j = 0 ; j < nb ; ++j ) {
        BCR.B(k,i,j) = rand(-0.1,0.1) ;
        BCR.C(k,j,i) = rand(-0.1,0.1) ;
      }
    }
  }
  for ( int i = 0 ; i < nb ; ++i ) {
    for ( int j = 0 ; j < nb ; ++j ) {
      BCR.F(i,j) = rand(-0.1,0.1) ;
    }
    for ( int j = 0 ; j < q ; ++j ) {
      BCR.Cq(i,j) = rand(-0.1,0.1) ;
    }
    BCR.F(i,i) += diag ; // force diagonal dominance
  }
  for ( int i = 0 ; i < n ; ++i ) {
    for ( int j = 0 ; j < nb ; ++j ) {
      BCR.C(nblock,j,i) = 1 ;
    }
  }

  cout << "N = " << N << '\n' ;

  for ( alglin::integer i = 0 ; i < N ; ++i ) x[i] = i % 100 ;
  std::copy( x, x+N, xref ) ;
  BCR.Mv( x, rhs ) ;

  cout << "nblock = " << nblock << "\n"
       << "n      = " << n      << "\n"
       << "q      = " << q      << "\n"
       << "nb     = " << nb     << "\n" ;
  /*
  ofstream file("mat.txt") ;
  file.precision(15) ;
  BCR.dump_ccoord( file ) ;
  file.close() ;
  file.open("rhs.txt") ;
  file.precision(15) ;
  for ( int i = 0 ; i < N ; ++i )
    file << rhs[i] << '\n' ;
  file.close() ;
  */

  TicToc tm ;
  tm.reset() ;
  tm.tic() ;
  BCR.factorize() ;
  tm.toc() ;
  cout << "\nReduction = " << tm.elapsedMilliseconds() << " [ms]\n\n" ;

  std::copy( rhs, rhs+N, x ) ;
  std::copy( rhs, rhs+N, x+N ) ;
  tm.tic() ;
  //BCR.solve( 2, x, N ) ;
  BCR.solve( x ) ;
  tm.toc() ;
  cout << "\nSolve = " << tm.elapsedMilliseconds() << " [ms]\n\n" ;

  /*
  file.open("sol.txt") ;
  file.precision(15) ;
  for ( int i = 0 ; i < N ; ++i )
    file << x[i] << '\n' ;
  file.close() ;
  */
  //for ( int i = 0 ; i < N ; ++i )
  //  cout << i << " " << x[i] << '\n' ;

  alglin::copy( N, xref, 1, xref1, 1 ) ;
  alglin::axpy( N, -1.0, x, 1, xref1, 1 ) ;
  cout << "Check |err|_inf = " << alglin::absmax( N, xref1, 1 ) << '\n' ;

  cout << "All done!\n" ;

  return 0 ;
}
