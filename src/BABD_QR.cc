/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

#include "BABD_QR.hh"
#include <iostream>

/*\
 |  reduces a 2x3 block matrix / Ad  Au     \ to a 1x2 block matrix ( Ad' Au' )
 |                             \     Bd  Bu /
 |  
 |  by using the following algorithm:
 |
 |    P * / Au \ = / Au' \  = / M \ (U) = / G \ (L*U)
 |        \ Bd /   \ Bd' /    \ L /       \ I /
 |
 |  where G = M * L^(-1).
 |
 |    / I -G \ * P * / Au \ = / I -G \ / G \ (L*U) = /  0  \
 |    \    I /       \ Bd /   \    I / \ I /         \ L*U /
 |
 |
 |    P * / Ad Au    \ = / Ad'    0   Bu'  \
 |        \    Bd Bu /   \ Ad''  L*U  Bu'' /
 |
 |    / I -G \ * P * / Ad Au    \ = / E*    0   F*   \
 |    \    I /       \    Bd Bu /   \ Ad'' L*U  Bu'' /
 |
 |  where  E* = Ad' - G*Ad'',  F* = Bu' - G*Bu''
 |
 | ------------------------------------------------------------
 |
 |  Esempio
 |  +-----+-----+
 |  |  E  |  F  |
 |  +-----+-----+-----+
 |        |  E  |  F  |
 |        +-----+-----+-----+
 |              |  E  |  F  |
 |              +-----+-----+-----+
 |                    |  E  |  F  |
 |                    +-----+-----+-----+
 |                          |  E  |  F  |
 |                          +-----+-----+-----+
 |                                |  E  |  F  |
 |  +----+                        +-----+-----+--+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
 |  Applicazione riduzione
 |  ----------------------
 |
 |  Q^T * / E  F  0  \ = / E'  R  F' \
 |        \ 0  E1 F1 /   \ E*  0  F* /
 |
 |  Sistema ridotto e memorizzazione
 |  +-----+-----+
 |  | E*  | F* |
 |  +-----+-----+-----+
 |        | E*  | F*  |
 |        +-----+-----+-----+
 |              | E*  | F*  |
 |  +----+      +-----+-----+--+
 |  | H0 |            | HN  |  |
 |  |    |            |     |  |
 |  +----+            +-----+--+
\*/

namespace alglin {

  using namespace std ;

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename QR_type>
  void
  BabdQR<QR_type>::factorize() {
    this->reduce() ;
    this->last_block_factorize() ;
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  template <typename QR_type>
  void
  BabdQR<QR_type>::solve( valuePointer y ) const {

    integer const & nblock = this->nblock ;
    integer const & n      = this->n ;

    this->forward(y) ;
    
    valuePointer ye   = y + nblock * n ;
    valuePointer ytmp = ye - n ;

    swap( n, y, 1, ytmp, 1 ) ;
    this->la_factorization->solve( ytmp ) ;
    swap( n, y, 1, ytmp, 1 ) ;

    this->backward( y ) ;
  }

  template <typename QR_type>
  void
  BabdQR<QR_type>::solve( integer      nrhs,
                          valuePointer y,
                          integer      ldY ) const {

    integer const & nblock = this->nblock ;
    integer const & n      = this->n ;

    this->forward( nrhs, y, ldY ) ;
    
    valuePointer ye   = y + nblock * n ;
    valuePointer ytmp = ye - n ;

    for ( integer j = 0 ; j < nrhs ; ++j )
      swap( n, y + j*ldY, 1, ytmp + j*ldY, 1 ) ;
    this->la_factorization->solve( nrhs, ytmp, ldY ) ;
    for ( integer j = 0 ; j < nrhs ; ++j )
      swap( n, y + j*ldY, 1, ytmp + j*ldY, 1 ) ;

    this->backward( nrhs, y, ldY ) ;
  }

  template class BabdQR<QR<double> > ;
  template class BabdQR<QR<float> > ;
  template class BabdQR<QRP<double> > ;
  template class BabdQR<QRP<float> > ;

}
