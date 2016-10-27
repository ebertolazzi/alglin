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

#include "BABD_Amodio.hh"
#include <iostream>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

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
 |  |  Ad | Au  |
 |  +-----+-----+-----+
 |        | Ad  | Au  |
 |        +-----+-----+-----+
 |              | Ad  | Au  |
 |              +-----+-----+-----+
 |                    | Ad  | Au  |
 |                    +-----+-----+-----+
 |                          | Ad  | Au  |
 |                          +-----+-----+-----+
 |                                | Ad  | Au  |
 |  +----+                        +-----+-----+--+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
 |  Applicazione riduzione
 |  ----------------------
 |
 |  Au'' = selezione righe di Au --> / Au'' e Ad'' possono essere compattate
 |  Ad'' = selezione righe di Ad --> \ nella stessa matrice
 |
 |  P * / Ad  0 \ = / Ad'   Bu'  \ ===> P * / Ad \ = / Adu'  \
 |      \ 0  Bu /   \ Ad''  Bu'' /          \ Bu /   \ Adu'' /
 |
 |  where  E* = Ad' - G*Ad'',  F* = Bu' - G*Bu''
 |
 |  La matrice G è spachettata nel blocco
 |
 |  +-----+-----+ - - +
 |  | E*  |  0  | F*   <-- messo al posto dell "0"
 |  +-----+-----+-----+
 |    Ad' | L\U | Au' | <-- memorizzazione compatta
 |  + - - +-----+-----+-----+ - - +
 |              | E*  |  0  | F*
 |              +-----+-----+-----+
 |                Ad' | L\U | Au' |
 |              + - - +-----+-----+-----+ - - +
 |                          | E*  |  0  | F*
 |                          +-----+-----+-----+
 |                            Ad' | L\U | Au' |
 |  +----+                  + - - +-----+-----+--+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
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
 |
 |
 |  +-----+-----+
 |  | E*  | F*  |
 |  +-----+-----+-----+                            +===+
 |        | L\U | Adu'|                            | V | <-- memorizzazione compatta e fill in
 |        +-----+-----+-----+                      +===+
 |              | E*  | F*  |
 |              +-----+-----+-----+                +===+
 |                    | L\U | Adu'|                | V | <-- memorizzazione compatta e fill in
 |                    +-----+-----+-----+          +===+
 |                          | E*  | F*  |
 |                          +-----+-----+-----+    +===+
 |                                | L\U | Adu'|    | V | <-- memorizzazione compatta e fill in
 |  +----+                        +-----+-----+--+ +===+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
\*/

namespace alglin {

  using namespace std ;

  /*\
   |   _                 _ ____        _   _
   |  | | ___   __ _  __| | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` |  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_|____/ \___/ \__|\__\___/|_| |_| |_|
  \*/

  template <typename t_Value>
  void
  AmodioLU<t_Value>::loadBottom( integer           q,
                                 valueConstPointer H0, integer ld0,
                                 valueConstPointer HN, integer ldN,
                                 valueConstPointer Hq, integer ldQ ) {

    integer & n = this->n ;

    m = n+q ;

    baseValue.allocate( size_t((m+n)*(m+1)) ) ;
    tmpV = baseValue( size_t(n+m) ) ;
    H0Nq = baseValue( size_t(m*(n+m)) ) ;

    gecopy( m, n, H0, ld0, H0Nq,       m ) ;
    gecopy( m, n, HN, ldN, H0Nq+m*n,   m ) ;
    gecopy( m, q, Hq, ldQ, H0Nq+2*m*n, m ) ;
  }

  /*\
   |   _                 _ _____           ____        _   _
   |  | | ___   __ _  __| |_   _|__  _ __ | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` | | |/ _ \| '_ \|  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | | | (_) | |_) | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_| |_|\___/| .__/|____/ \___/ \__|\__\___/|_| |_| |_|
   |                                |_|
  \*/

  template <typename t_Value>
  void
  AmodioLU<t_Value>::loadTopBottom( // ----------------------------
                                    integer           row0,
                                    integer           col0,
                                    valueConstPointer block0,
                                    integer           ld0,
                                    // ----------------------------
                                    integer           rowN,
                                    integer           colN,
                                    valueConstPointer blockN,
                                    integer           ldN ) {

    integer & n = this->n ;

    m = col0+colN-n ;

    baseValue.allocate( size_t((m+n)*(m+1)) ) ;
    tmpV = baseValue( size_t(n+m) ) ;
    H0Nq = baseValue( size_t(m*(n+m)) ) ;
    
    /*\
     |  +----+-----+---+
     |  | H0 | HN  |Hq |
     |  |    |     |   |
     |  +----+-----+---+
     |  +----+------+---------+
     |  | 0  | blkN |blkN: 0  |
     |  |blk0|  0   |    :blk0|
     |  +----+------+---------+
    \*/

    zero( m*(n+m), H0Nq, 1 ) ;
    gecopy( rowN, colN, blockN, ldN, H0Nq+m*n, m ) ;

    valuePointer H0 = H0Nq+rowN ;
    valuePointer Hq = H0+m*(n+colN) ;
    gecopy( row0, col0-n, block0,              ld0, Hq, m ) ;
    gecopy( row0, n,      block0+(col0-n)*ld0, ld0, H0, m ) ;

  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  AmodioLU<t_Value>::factorize() {

    integer & n = this->n ;
    this->reduce() ;

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    this->la_factorization->allocate(n+m,n+m) ;
    this->la_factorization->load_block( n, 2*n, this->getPointer_LR(), n ) ;
    this->la_factorization->load_block( m, m+n, H0Nq, m, n, 0 ) ;
    if ( m > n ) this->la_factorization->zero_block( n, m-n, 0, 2*n ) ;

    // fattorizzazione ultimo blocco
    this->la_factorization->factorize() ;
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  template <typename t_Value>
  void
  AmodioLU<t_Value>::solve( valuePointer y ) const {

    integer const & n = this->n ;

    this->forward(y) ;
    
    valuePointer V0 = tmpV ;
    valuePointer V1 = tmpV+n ;
    
    valuePointer ye = y + this->getNblock() * n ;
    copy( n, y,  1, V0, 1 ) ;
    copy( m, ye, 1, V1, 1 ) ;

    this->la_factorization->solve( V0 ) ;

    copy( n, V0, 1, y,  1 ) ;
    copy( m, V1, 1, ye, 1 ) ;

    this->backward( y ) ;
  }

  template class AmodioLU<double> ;
  template class AmodioLU<float> ;

}
