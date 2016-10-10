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
#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_aux.hh"
#include "TicToc.hh"
#include "KKT_like.hh"

using namespace std ;
typedef double valueType ;
using alglin::integer ;

static
void
test0() {

  alglin::KKT<valueType> kkt ;
  integer const N   = 2 ;
  integer const M   = 1 ;

  valueType A[] = {
    1,      0,
    0,      1,
  } ;

  valueType B[] = { 0, 0 } ;
  valueType C[] = { 0, 0 } ;
  valueType D[] = { 2 } ;
  
  // 1 0 0 -> 1
  // 0 1 0 -> 2
  // 0 0 2 -> 6

  kkt.factorize( N, M,
                 A, N, false,
                 B, M, false,
                 C, N, false,
                 D, M, false ) ;
  valueType x[N+M] = { 1, 2, 3 } ;
  valueType rhs[N+M] ;

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n' ;
  kkt.solve( rhs ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n' ;

}

static
void
test1() {

  alglin::KKT<valueType> kkt ;
  integer const N   = 2 ;
  integer const M   = 1 ;

  valueType A[] = {
    1,      3,
   -2,      1,
  } ;

  valueType B[] = { 1, -1 } ;
  valueType C[] = { -1, 1 } ;
  valueType D[] = { 2 } ;
  
  //  1 -2  1 -> 0
  //  3  1 -1 -> 2
  // -1  1  2 -> 7

  kkt.factorize( N, M,
                 A, N, false,
                 B, M, false,
                 C, N, false,
                 D, M, false ) ;
  valueType x[N+M] = { 1, 2, 3 } ;
  valueType rhs[N+M] ;

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n' ;
  kkt.solve( rhs ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n' ;

}

static
void
test2() {

  alglin::KKT<valueType> kkt ;
  integer const N = 3 ;
  integer const M = 2 ;

  valueType A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  } ;

  valueType B[] = {
    0.001,      3,
    0.001,      -0.001,
    1,          2,
  } ;

  valueType C[] = {
    2,       3,
    -0.001,  -0.001,
    1,        2,
  } ;

  valueType D[] = {
    1,      4,
    -1,      1
  } ;
  
  kkt.factorize( N, M,
                 A, N, false,
                 B, N, false,
                 C, M, false,
                 D, M, false ) ;
  valueType x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 } ;
  valueType rhs[2*(N+M)] ;

  alglin::gemv(alglin::NO_TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               N, M, 1.0, B, N,
               x+N, 1,
               1,
               rhs, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, N, 1.0, C, M,
               x, 1,
               0,
               rhs+N, 1 ) ;
  alglin::gemv(alglin::NO_TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               1,
               rhs+N, 1 ) ;
  std::copy( rhs, rhs + N+M, rhs+N+M ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n' ;
  //kkt.solve( rhs ) ;
  kkt.solve( 2, rhs, N+M ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n' ;

}

static
void
test3() {

  alglin::KKT<valueType> kkt ;
  integer const N = 3 ;
  integer const M = 2 ;

  valueType A[] = {
    0.001,      2,     3,
    0.001,      0.001, 0,
    0,          0.001, 2,
  } ;

  valueType B[] = {
    0.001,      3,
    0.001,     -0.001,
    1,          2,
  } ;

  valueType C[] = {
     2,       3,
    -0.001,  -0.001,
     1,       2,
  } ;

  valueType D[] = {
    1,      4,
    -1,      1
  } ;
  
  kkt.factorize( N, M,
                 A, N, false,
                 B, N, false,
                 C, M, false,
                 D, M, false ) ;
  valueType x[] = { 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 } ;
  valueType rhs[2*(N+M)] ;

  alglin::gemv(alglin::TRANSPOSE,
               N, N, 1.0, A, N,
               x, 1,
               0,
               rhs, 1 ) ;
  alglin::gemv(alglin::TRANSPOSE,
               M, N, 1.0, C, M,
               x+N, 1,
               1,
               rhs, 1 ) ;

  alglin::gemv(alglin::TRANSPOSE,
               M, M, 1.0, D, M,
               x+N, 1,
               0,
               rhs+N, 1 ) ;

  alglin::gemv(alglin::TRANSPOSE,
               N, M, 1.0, B, N,
               x, 1,
               1,
               rhs+N, 1 ) ;

  std::copy( rhs, rhs+N+M, rhs+N+M ) ;

  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "rhs[" << i << "] = " << rhs[i] << '\n' ;
  kkt.t_solve( 2, rhs, N+M ) ;
  for ( integer i = 0 ; i < N+M ; ++i )
    cout << "x[" << i << "] = " << rhs[i] << '\n' ;
}


int
main() {

  try {
    test0() ;
    test1() ;
    test2() ;
    test3() ;
  } catch ( exception const & exc ) {
    cerr << exc.what() << '\n' ;
  } catch ( ... ) {
    cerr << "Errore Sconosciuto!\n" ;
  }
  return 0 ;
}
