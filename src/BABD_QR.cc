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
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
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

  /*  
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */

  template <typename QR_type>
  void
  BabdQR<QR_type>::factorize( BABD_QR_LASTBLOCK_Choice choice,
                              // ----------------------------
                              integer           nblock,
                              integer           _n,
                              integer           q,
                              valueConstPointer AdAu,
                              valueConstPointer H0,
                              valueConstPointer HN,
                              valueConstPointer Hq ) {

    n = _n ;
    m = n+q ;
    last_block = choice ;
    CR.reduce( nblock, n, AdAu, n ) ;

    // fattorizzazione ultimo blocco
    switch ( last_block ) {
      case BABD_QR_LASTBLOCK_LU:  factorization = & la_lu  ; break ;
      case BABD_QR_LASTBLOCK_QR:  factorization = & la_qr  ; break ;
      case BABD_QR_LASTBLOCK_QRP: factorization = & la_qrp ; break ;
      case BABD_QR_LASTBLOCK_SVD: factorization = & la_svd ; break ;
      ALGLIN_ERROR("AmodioLU<t_Value>::factorize -- no last block solver selected") ;
    }

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    tmpV.resize( n+m ) ;
    factorization->allocate(n+m,n+m) ;
    factorization->load_block( n, 2*n, CR.getPointer_LR(), n ) ;
    factorization->load_block( m, n, H0, m, n, 0 ) ;
    factorization->load_block( m, n, HN, m, n, n ) ;
    if ( q > 0 ) {
      factorization->zero_block( n, q, 0, 2*n ) ;
      factorization->load_block( m, q, Hq, m, n, 2*n ) ;
    }
    // fattorizzazione ultimo blocco
    factorization->factorize() ;
  }

  /*             _
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */

  template <typename QR_type>
  void
  BabdQR<QR_type>::solve( valuePointer y ) const {

    CR.forward(y) ;
    
    valuePointer V0 = &tmpV[0] ;
    valuePointer V1 = &tmpV[n] ;
    
    valuePointer ye = y + CR.getNblock() * n ;
    copy( n, y,  1, V0, 1 ) ;
    copy( m, ye, 1, V1, 1 ) ;

    factorization->solve( V0 ) ;

    copy( n, V0, 1, y,  1 ) ;
    copy( m, V1, 1, ye, 1 ) ;

    CR.backward( y ) ;
  }


  template class BabdQR<QR<double> > ;
  template class BabdQR<QR<float> > ;
  template class BabdQR<QRP<double> > ;
  template class BabdQR<QRP<float> > ;

}
