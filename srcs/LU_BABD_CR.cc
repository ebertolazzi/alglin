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

#define EIGEN_NO_DEBUG
#define NDEBUG

#include "LU_BABD_CR.hh"
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

  #ifdef LU_BABD_CR_USE_THREAD
  template <typename t_Value, integer n>
  CyclicReductionLU<t_Value,n>::CyclicReductionLU( integer nth )
  : baseValue("CyclicReductionLU_value")
  , baseInteger("CyclicReductionLU_index")
  , NB(25)
  , numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= LU_BABD_CR_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << LU_BABD_CR_MAX_THREAD ) ;
  }
  #else
  template <typename t_Value, integer n>
  CyclicReductionLU<t_Value,n>::CyclicReductionLU()
  : NB(25)
  { }
  #endif

  template <typename t_Value, integer n>
  CyclicReductionLU<t_Value,n>::~CyclicReductionLU() {
  }

  /*
  //
  // +--------+
  // | \______|
  // | | A    |
  // | |      |
  // +--------+
  // | |      |
  // | | B    |
  // | |      |
  // +--------+
  //
  */
  template <typename t_Value, integer n>
  integer
  CyclicReductionLU<t_Value,n>::LU_2_block( mat_t  & A,
                                            mat_t  & B,
                                            ivec_t & ipiv ) const {
    // LU DECOMPOSITION, ROW INTERCHANGES
    //integer nb = std::min(NB,n) ;
    for ( integer j = 0 ; j < n ; ) {
      integer mx1 = j ;
      t_Value ax1 = std::abs( A(j,j) ) ;
      for ( integer i = j+1 ; i < n ; ++i ) {
        t_Value v = std::abs( A(i,j) ) ;
        if ( v > ax1 ) { ax1 = v ; mx1 = i ; }
      }
      integer mx2 = 0 ;
      t_Value ax2 = std::abs( B(0,j) ) ;
      for ( integer i = 1 ; i < n ; ++i ) {
        t_Value v = std::abs( B(i,j) ) ;
        if ( v > ax2 ) { ax2 = v ; mx2 = i ; }
      }
      if ( ax1 < ax2 ) {
        ipiv(j) = mx2 + n ;
        A.row(j).swap(B.row(mx2));
      } else {
        ipiv(j) = mx1 ;
        if ( mx1 > j ) A.row(j).swap(A.row(mx1)) ;
      }
      if ( A(j,j) == 0 ) return j+1 ;
      valueType ROWM = 1/A(j,j) ;
      A.block(j+1,j,n-j-1,1) *= ROWM ;
      B.block(0,j,n,1)       *= ROWM ;
      ++j ;
      A.block(j,j,n-j,n-j).noalias() -=
        A.block(j,j-1,n-j,1) * A.block(j-1,j,1,n-j) ;
      B.block(0,j,n,n-j).noalias() -=
        B.block(0,j-1,n,1) * A.block(j-1,j,1,n-j) ;
    }
    /*
    // compute G = B*L^(-1)
    //
    //  / L \ (U) = / I          \ (L*U) =  / I \ (L*U)
    //  \ M /       \ M * L^(-1) /          \ G /
    */
    A.template triangularView<Eigen::UnitLower>()
     .template solveInPlace<Eigen::OnTheRight>(B) ;
    return 0 ;
  }

  /*  
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */

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
  
  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::reduction() {
    do {

      if ( jump_block >= nblock ) break ;
      k_block = 0 ;

      do {

        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        ivec_t & ipiv = ipiv_blk[k1] ;
        std::vector<bool> & LU_rows = LU_rows_blk[k1] ;

        mat_t & G   = AdAuG_blk[3*k1+2] ;
        mat_t & E   = AdAuG_blk[3*k]    ; // <-- Ad'' - G*Ad'
        mat_t & F   = AdAuG_blk[3*k+1]  ; // <-- Bu'' - G*Bu'
        mat_t & LU  = AdAuG_blk[3*k1]   ;
        mat_t & Adu = AdAuG_blk[3*k1+1] ;

        /*\
         | factorize RS by means of the LU factorization
         | P * / Au \ = / G \ (L*U)
         |     \ Bd /   \ I /
        \*/
        G = F ;
        integer info = LU_2_block( LU, G, ipiv ) ;
        ALGLIN_ASSERT( info == 0,
                      "CyclicReductionLU::factorize, at block N." << k <<
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
            if ( ip < n ) Adu.row(i).swap(Adu.row(ip)) ;
            else          Adu.row(i).swap(E.row(ip-n)) ;
          }
        }
        F.setZero() ;
        EE.setZero() ;
        FF.setZero() ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            F.row(i).noalias() = E.row(i) ;
            E.row(i).setZero() ;
          }
          if ( LU_rows[i] ) FF.row(i).noalias() = Adu.row(i) ;
          else              EE.row(i).noalias() = Adu.row(i) ;
        }
        E -= G*EE ;
        F -= G*FF ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;
      jump_block *= 2 ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_CR_USE_THREAD
  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::reduction_mt( integer num_thread, integer nth ) {
    valuePointer EE = tmpM+nth*nxnx2 ;
    valuePointer FF = EE+nxn ;

    do {

      if ( jump_block >= nblock ) break ;

      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      do {

        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

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
                       "CyclicReductionLU::factorize, at block N." << k <<
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
      } while ( true ) ;

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
        } else {
          cond0.wait(lck);
        }
      }
    } while ( true ) ;
  }
  #endif

  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::factorize( integer           _nblock,
                                           integer           _q,
                                           valueConstPointer AdAu,
                                           valueConstPointer H0,
                                           valueConstPointer HN,
                                           valueConstPointer Hq ) {

    nblock = _nblock ;
    q      = _q ;

    #ifdef LU_BABD_CR_USE_THREAD
    integer nnzLU = numThread*nxnx2 ;
    if ( nnzLU < nm*nm ) nnzLU = nm*nm ;
    #endif

    AdAuG_blk.resize(3*nblock) ;
    ipiv_blk.resize( nblock ) ;
    LU_rows_blk.resize( nblock ) ;
    EE.resize(n,n) ;
    FF.resize(n,n) ;
    tmpV.resize(2*n+q) ;
    tmpV1.resize(2*n+q) ;
    integer nxn   = n*n ;
    integer nxnx2 = 2*n*n ;
    for ( integer i = 0 ; i < nblock ; ++i ) {
      LU_rows_blk[i].resize(2*n) ;
      AdAuG_blk[3*i].resize(n,n) ;
      AdAuG_blk[3*i+1].resize(n,n) ;
      AdAuG_blk[3*i+2].resize(n,n) ;
      ipiv_blk[i].resize(n) ;
      alglin::copy( nxn, AdAu + i*nxnx2,     1, AdAuG_blk[3*i].data(),   1 ) ;
      alglin::copy( nxn, AdAu + i*nxnx2+nxn, 1, AdAuG_blk[3*i+1].data(), 1 ) ;
    }

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef LU_BABD_CR_USE_THREAD
    integer usedThread = numThread ;
    #endif

    jump_block = 1 ;
    #ifdef LU_BABD_CR_USE_THREAD
    to_be_done = 0 ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &CyclicReductionLU<t_Value>::reduction_mt, this, usedThread, nt ) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #else
    reduction() ;
    #endif

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    integer nm = 2*n+q ;
    integer m  = n+q ;
    dmat_t tmp(nm,nm) ;
    tmp.block(0,0,n,n) = AdAuG_blk[0] ;
    tmp.block(0,n,n,n) = AdAuG_blk[1] ;
    tmp.block(n,0,m,n) = Eigen::Map<dmat_t const>(H0,m,n) ;
    tmp.block(n,n,m,n) = Eigen::Map<dmat_t const>(HN,m,n) ;
    if ( _q > 0 ) {
      tmp.block(0,2*n,n,_q).setZero() ;
      tmp.block(n,2*n,m,_q) = Eigen::Map<dmat_t const>(Hq,m,_q) ;
    }

    LU_last.compute(tmp) ;
    ALGLIN_ASSERT( LU_last.isInvertible(),
                   "CyclicReductionLU::factorize(), singular matrix" ) ;

  }

  /*
  //  __                                  _                  _
  // / _| ___  _ ____      ____ _ _ __ __| |    _ __ ___  __| |_   _  ___ ___
  // | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |   | '__/ _ \/ _` | | | |/ __/ _ \
  // |  _| (_) | |   \ V  V / (_| | | | (_| |   | | |  __/ (_| | |_| | (_|  __/
  // |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|___|_|  \___|\__,_|\__,_|\___\___|
  //                                       |_____|
  */

  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::forward_reduce( valuePointer y ) const {
    do {
      k_block = 0 ;
      if ( jump_block >= nblock ) break ;

      do {
        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;
        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        mat_t  const & G    = AdAuG_blk[3*k1+2] ;
        ivec_t const & ipiv = ipiv_blk[k1] ;

        Eigen::Map<vec_t> yk(y+k*n) ;
        Eigen::Map<vec_t> yk1(y+k1*n) ;

        /*\
         |   applico permutazione e moltiplico per / I -G \
         |                                         \    I /
        \*/
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            if ( ip < n ) std::swap( yk1(i), yk1(ip) ) ;
            else          std::swap( yk1(i), yk(ip-n) ) ;
          }
        }
        yk -= G * yk1 ;
      } while ( true ) ;

      // aspetta le altre thread
      jump_block *= 2 ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_CR_USE_THREAD
  template <typename t_Value>
  void
  CyclicReductionLU<t_Value>::forward_reduce_mt( integer num_thread, valuePointer y ) const {
    do {
      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      if ( jump_block >= nblock ) break ;

      do {
        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }
        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

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
      } while ( true ) ;

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
        } else {
          cond0.wait(lck);
        }
      }
    } while ( true ) ;
  }
  #endif

  /*
  //   _                _                  _         _   _ _         _
  //  | |__   __ _  ___| | __    ___ _   _| |__  ___| |_(_) |_ _   _| |_ ___
  //  | '_ \ / _` |/ __| |/ /   / __| | | | '_ \/ __| __| | __| | | | __/ _ \
  //  | |_) | (_| | (__|   <    \__ \ |_| | |_) \__ \ |_| | |_| |_| | ||  __/
  //  |_.__/ \__,_|\___|_|\_\___|___/\__,_|_.__/|___/\__|_|\__|\__,_|\__\___|
  //                       |_____|
  */

  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::back_substitute( valuePointer y ) const {
    do {
      jump_block /= 2 ;
      k_block    = 0 ;
      if ( jump_block == 0 ) break ;

      do {
        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;
      
        integer k2 = min(k1+jump_block,nblock) ;

        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;
        mat_t  const & LU  = AdAuG_blk[3*k1] ;
        mat_t  const & Adu = AdAuG_blk[3*k1+1] ;

        Eigen::Map<vec_t> yk(y+k*n,n) ;
        Eigen::Map<vec_t> yk1(y+k1*n,n) ;
        Eigen::Map<vec_t> yk2(y+k2*n,n) ;

        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i] ) yk1(i) -= yk2.dot(Adu.row(i)) ;
          else              yk1(i) -= yk.dot(Adu.row(i)) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        LU.template triangularView<Eigen::UnitLower>()
          .template solveInPlace<Eigen::OnTheLeft>(yk1) ;
        LU.template triangularView<Eigen::Upper>()
          .template solveInPlace<Eigen::OnTheLeft>(yk1) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_CR_USE_THREAD

  template <typename t_Value>
  void
  CyclicReductionLU<t_Value>::back_substitute_mt( integer num_thread, valuePointer y ) const {
    do {
      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          jump_block /= 2 ;
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      if ( jump_block == 0 ) break ;

      do {
        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;
      
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = LU_rows[i] ? yk2 : yk ;
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
      } while ( true ) ;
      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) cond0.notify_all() ; // wake up all tread
        else                     cond0.wait(lck);
      }
    } while ( true ) ;
  }
  #endif

  /*             _
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */

  template <typename t_Value, integer n>
  void
  CyclicReductionLU<t_Value,n>::solve( valuePointer y ) const {

    #ifdef LU_BABD_CR_USE_THREAD
    integer usedThread = numThread ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef LU_BABD_CR_USE_THREAD
    if ( usedThread == 1 ) {
      forward_reduce(y) ;
    } else {
      to_be_done = 0 ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &CyclicReductionLU<t_Value>::forward_reduce_mt, this, usedThread, y ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #else
    forward_reduce(y) ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! 2 by 2 block linear system solution !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    integer m = n+q ;
    Eigen::Map<dvec_t> yb(y, n) ;
    Eigen::Map<dvec_t> ye(y + nblock * n, m) ;
    tmpV.head(n) = yb ;
    tmpV.tail(m) = ye ;
    tmpV1 = LU_last.solve(tmpV) ;
    yb = tmpV1.head(n) ;
    ye = tmpV1.tail(m) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    #ifdef LU_BABD_CR_USE_THREAD
    if ( usedThread == 1 ) {
      back_substitute( y ) ;
    } else {
      to_be_done = 0 ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &CyclicReductionLU<t_Value>::back_substitute_mt, this, usedThread, y ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #else
    back_substitute( y ) ;
    #endif
  }

  template class CyclicReductionLU<double,2> ;
  template class CyclicReductionLU<double,4> ;
  template class CyclicReductionLU<double,6> ;
  template class CyclicReductionLU<float,2> ;

}
