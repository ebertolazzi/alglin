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

#include "LU_BABD_Amodio.hh"

#include <iostream>
using namespace std ;

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
  
  /*
  //
  // +--------+
  // |        |
  // |   A    |
  // |        |
  // +--------+
  // |        |
  // |   B    |
  // |        |
  // +--------+
  //
  */
  template <typename t_Value>
  integer
  AmodioLU<t_Value>::LU_2_block( integer      n,
                                 valuePointer A,
                                 valuePointer B,
                                 integer      ipiv[] ) const {
    // LU DECOMPOSITION, ROW INTERCHANGES
    valuePointer Ajj = A ;
    valuePointer Bj  = B ;
    for ( integer j = 0 ; j < n ; Ajj += n+1, Bj += n  ) {
      integer MX1 = iamax( n-j, Ajj, 1 ) ;
      integer MX2 = iamax( n,   Bj,  1 ) ;
      if ( std::abs(Ajj[MX1]) < std::abs(Bj[MX2]) ) {
        ipiv[j] = MX2 + n ; // C-based
        swap( n, A + j, n, B + MX2, n ) ;
      } else {
        ipiv[j] = MX1 + j ; // C-based
        if ( MX1 > 0 ) swap( n, A + j, n, A + ipiv[j], n ) ;
      }
      if ( std::abs(Ajj[0]) == 0 ) return j+1 ;
      valueType ROWM = 1/Ajj[0] ;
      ++j ;
      scal(n-j, ROWM, Ajj+1, 1) ;
      scal(n,   ROWM, Bj,    1) ;
      ger(n-j, n-j, -1.0, Ajj+1, 1, Ajj+n, n, Ajj+n+1, n ) ;
      ger(n,   n-j, -1.0, Bj,    1, Ajj+n, n, Bj+n,    n ) ;
    }
    /*
    // compute G = B*L^(-1)
    //
    //  / L \ (U) = / I          \ (L*U) =  / I \ (L*U)
    //  \ M /       \ M * L^(-1) /          \ G /
    */
    trsm( SideMultiply::RIGHT,
          ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          n, n, 1.0, A, n, B, n ) ;
    return 0 ;
  }

  /*  
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */
  template <typename t_Value>
  void
  AmodioLU<t_Value>::factorize( integer           nblock,
                                integer           n,
                                integer           q,
                                valueConstPointer AdAu,
                                valueConstPointer H0,
                                valueConstPointer HN,
                                valueConstPointer Hq ) {

    this -> nblock = nblock ;
    this -> n      = n ;
    this -> m      = n+q ;

    integer nx2    = 2*n ;
    integer nxn    = n*n ;
    integer nxnx2  = nxn*2 ;
    integer nm     = n+m ;

    integer nnzG    = (nblock-1)*nxn ;
    integer nnzADAU = nblock*nxnx2 ;
    integer nnzLU   = nm*nm ;

    integer nv = nnzG + nnzADAU + nnzLU + nm ;
    integer ni = nblock*n + nm ;

    baseValue   . allocate(nv) ;
    baseInteger . allocate(ni) ;

    AdAu_blk = baseValue( nnzADAU ) ;
    G_blk    = baseValue( nnzG ) - nxn ; // 1 based
    LU_blk   = baseValue( nnzLU ) ;
    tmpV     = baseValue( nm ) ;

    ipiv_blk    = baseInteger( nblock*n ) ;
    LU_ipiv_blk = baseInteger( nm ) ;

    LU_rows_blk.resize( nblock ) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      LU_rows_blk[i].resize(nx2) ;

    alglin::copy( nnzADAU, AdAu, 1, AdAu_blk, 1 ) ;

    /*\
     |  Based on the algorithm 
     |  Pierluigi Amodio and Giuseppe Romanazzi
     |  Algorithm 859: BABDCR: a Fortran 90 package for the Solution of Bordered ABD Linear Systems
     |  ACM Transactions on Mathematical Software, 32, 4, 597—608 (2006)
    \*/

    // 0  1  2  3  4  5  6
    // 0  -  2  -  4  -  6* unchanged
    // 0  -  -  -  4  -  6* unchanged
    // 0  -  -  -  -  -  6* special
    // 0  -  -  -  -  -  -

    // 0  1  2  3  4  5  6  7  8  9 10 11 12
    // 0  -  2  -  4  -  6  -  8  - 10  - 12* unchanged
    // 0  -  -  -  4  -  -  -  8  -  -  - 12* unchanged
    // 0  -  -  -  -  -  -  -  8  -  -  - 12* unchanged
    // 0  -  -  -  -  -  -  -  -  -  -  - 12* special
    // 0  -  -  -  -  -  -  -  -  -  -  -

    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13
    // 0  -  2  -  4  -  6  -  8  - 10  - 12  -
    // 0  -  -  -  4  -  -  -  8  -  -  - 12  -
    // 0  -  -  -  -  -  -  -  8  -  -  - 12  -
    // 0  -  -  -  -  -  -  -  -  -  -  - 12  -
    // 0  -  -  -  -  -  -  -  -  -  -  -  -  -

    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
    // 0  -  2  -  4  -  6  -  8  - 10  - 12  - 14  -
    // 0  -  -  -  4  -  -  -  8  -  -  - 12  -  -  -
    // 0  -  -  -  -  -  -  -  8  -  -  -  -  -  -  -
    // 0  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    integer jump = 1 ;
    for ( ; jump < nblock ; jump *= 2 ) {
      integer k  = 0 ;
      integer k1 = jump ;
      while ( k1 < nblock ) {

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;
        std::vector<bool> & LU_rows = LU_rows_blk[k1] ;

        valuePointer E   = BLK0 ;       // <-- Ad'' - G*Ad'
        valuePointer F   = BLK0 + nxn ; // <-- Bu'' - G*Bu'
        valuePointer LU  = BLK1 ;
        valuePointer Adu = BLK1 + nxn ;

        /*\
         | factorize RS by means of the LU factorization
         | P * / Au \ = / G \ (L*U)
         |     \ Bd /   \ I /
        \*/
        copy( nxn, F, 1, G, 1 ) ;
        integer info = LU_2_block( n, LU, G, ipiv ) ;
        ALGLIN_ASSERT( info == 0,
                       "AmodioLU::factorize, at block N." << k <<
                       " singular matrix, info = " << info );

        /*\
         |   +-----+-----+ - - +
         |   | E*  |  0 <=  F*   <-- messo al posto dell "0"
         |   +-----+-----+-----+
         |     Ad''| L\U | Au''| <-- memorizzazione compatta
         |   + - - +-----+-----+
         |
         |   applico permutazione
         |   +-----+             +-----+-----+
         |   | Ad  |             |  E  |  F  |
         |   +-----+-----+  -->  +-----+-----+
         |         | Bu  |       | Ad''+Bu'' | compattato
         |         +-----+       +-----------+
        \*/

        // determine the permutation vector and select row vectors
        for ( integer i = 0 ; i < n ; ++i ) {
          LU_rows[i]   = true ;
          LU_rows[n+i] = false ;
        }
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            std::swap( LU_rows[i], LU_rows[ip] ) ;
            if ( ip < n ) swap( n, Adu + i, n, Adu + ip, n ) ; // scambia righe
            else          swap( n, Adu + i, n, E + ip-n, n ) ; // scambia righe
          }
        }
        valuePointer EE = LU_blk ;
        valuePointer FF = LU_blk+nxn ;
        zero( nxn, F, 1 ) ;
        zero( nxnx2, EE, 1 ) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            copy( n, E+i, n, F+i, n ) ;
            zero( n, E+i, n ) ;
          }
          if ( LU_rows[i] ) copy( n, Adu+i, n, FF+i, n ) ;
          else              copy( n, Adu+i, n, EE+i, n ) ;
        }
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              FF, n,
              1.0, F, n ) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
        k  += 2*jump ;
        k1 += 2*jump ;
      }
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
    gecopy( n, nx2, AdAu_blk, n, LU_blk, nm ) ;
    gecopy( m, n, H0, m, LU_blk+n,      nm ) ;
    gecopy( m, n, HN, m, LU_blk+n+n*nm, nm ) ;
    if ( m > n ) {
      gezero( n, q,        LU_blk+nx2*nm,   nm ) ;
      gecopy( m, q, Hq, m, LU_blk+nx2*nm+n, nm ) ;
    }

    integer INFO = getrf( nm, nm, LU_blk, nm, LU_ipiv_blk ) ;
    ALGLIN_ASSERT( INFO==0,
                   "AmodioLU::factorize(), singular matrix, getrf INFO = " << INFO ) ;

  }
  
  /*             _           
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */
  template <typename t_Value>
  void
  AmodioLU<t_Value>::solve( valuePointer y ) {
    /*\
     |
     |  Purpose
     |  =======
     |
     |  BABDCR_SOLV  solves a babd linear system whose coefficient matrix
     |  has been factorized by BABDCR_FACT. The algorithm consists of three
     |  phases: reduction (by using the subroutine REDUCE_RHS), solution of
     |  the 2 by 2 block system (by using the subroutine DGETRS) and 
     |  back-substitution (by using the subroutine SOLVE_BLOCK).
     |
     |  In input, BABDCR_SOLV requires the coefficient matrix factorized by
     |  the subroutine BABDCR_FACT (see details of this subroutine) and the
     |  right hand side which is stored in the block vector VECT_B in the
     |  following form
     |
     |  VECT_B  = [f(0), f(1), ...., f(NBLOKS)].
     |
     |  The first block element f(0) corresponds to the first row of the
     |  coefficient matrix containing Ba and Bb.
     |
     |  On exit, BABDCR_SOLV gives the solution of the babd linear system.
     |  The solution is stored in the block vector VECT_B.
     |
    \*/
    // some constanst
    integer nxn   = n*n ;
    integer nxnx2 = nxn*2 ;
    integer nm    = n+m ;
    
    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    integer jump = 1 ;
    for ( ; jump < nblock ; jump *= 2 ) {
      integer k  = 0 ;
      integer k1 = jump ;

      while ( k1 < nblock ) {

        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;

        /*\
         |   applico permutazione e moltiplico per / I -G \
         |                                         \    I /
        \*/
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            if ( ip < n ) std::swap( yk1[i], yk1[ip] ) ;
            else          std::swap( yk1[i], yk[ip-n] ) ;
          }
        }
        // yk -= G * yk1
        gemv( Transposition::NO_TRANSPOSE,
              n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
        // -----------------
        k  += 2*jump ;
        k1 += 2*jump ;
      }
    }

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! 2 by 2 block linear system solution !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    valuePointer ye = y + nblock * n ;
    copy( n, y,  1, tmpV,   1 ) ;
    copy( m, ye, 1, tmpV+n, 1 ) ;
    integer INFO = getrs( Transposition::NO_TRANSPOSE,
                          nm, 1, LU_blk, nm, LU_ipiv_blk, tmpV, nm ) ;
    ALGLIN_ASSERT( INFO==0,
                   "AmodioLU::solve(), singular matrix, getrs INFO = " << INFO ) ;
    copy( n, tmpV,   1, y,  1 ) ;
    copy( m, tmpV+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    for ( jump /= 2 ; jump > 0 ; jump /= 2 ) {

      integer k  = 0 ;
      integer k1 = jump ;

      while ( k1 < nblock ) {

        integer      k2   = min(k1+jump,nblock) ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        std::vector<bool> & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        valuePointer LU  = BLK1 ;
        valuePointer Adu = BLK1 + nxn ;

        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = LU_rows[i] ? yk2 : yk ; // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          yk1[i] -= dot( n, Adu+i, n, LR, 1 ) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        trsv( ULselect::LOWER,
              Transposition::NO_TRANSPOSE,
              DiagonalType::UNIT,
              n, LU, n, yk1, 1 ) ;
        trsv( ULselect::UPPER,
              Transposition::NO_TRANSPOSE,
              DiagonalType::NON_UNIT,
              n, LU, n, yk1, 1 ) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
        // -----------------
        k  += 2*jump ;
        k1 += 2*jump ;
      }
    }

  }

  template class AmodioLU<double> ;
  template class AmodioLU<float> ;

}
