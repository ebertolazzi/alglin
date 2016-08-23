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

#include "LU_BABD_QR.hh"
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

  #ifdef LU_BABD_QR_USE_THREAD
  template <typename t_Value>
  BabdQR<t_Value>::BabdQR( integer nth )
  : baseValue("BabdQR_value")
  , baseInteger("BabdQR_index")
  , numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= LU_BABD_QR_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << LU_BABD_QR_MAX_THREAD ) ;
  }
  #else
  template <typename t_Value>
  BabdQR<t_Value>::BabdQR()
  : baseValue("BabdQR_value")
  { }
  #endif

  template <typename t_Value>
  BabdQR<t_Value>::~BabdQR() {
    baseValue.free() ;
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
  
  template <typename t_Value>
  void
  BabdQR<t_Value>::reduction() {

    while ( jump_block < nblock ) {

      integer kstep = 2*jump_block ;
      integer kend  = nblock-jump_block ;

      for ( integer k = 0 ; k < kend ; k += kstep ) {
        integer k1 = k+jump_block ;
        /*\
         |  +-----+-----+
         |  |  E  |  F  |
         |  +-----+-----+-----+
         |        |  E1 |  F1 |
         |        +-----+-----+
        \*/

        valuePointer E  = AdAu_blk + k  * nxnx2 ;
        valuePointer E1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer F  = E  + nxn ;
        valuePointer F1 = E1 + nxn ;

        QR_type & QR = QR_blk[k1] ;

        /*\
         |  / F \ = Q / R \
         |  \ E1/     \ 0 /
        \*/
        gecopy( n, n, F,  n, M1_2n_n.data(),   nx2 ) ;
        gecopy( n, n, E1, n, M1_2n_n.data()+n, nx2 ) ;
        QR.compute(M1_2n_n) ;

        /*\
         |   / 0 I \ Q^T / F \ = / 0 \
         |   \ I 0 /     \ E1/   \ R /
        \*/

        gecopy( n, n, E, n, M1_2n_n.data(),   nx2 ) ;
        gezero( n, n,       M1_2n_n.data()+n, nx2 ) ;
        M2_2n_n.noalias() = QR.householderQ().transpose() * M1_2n_n ;
        gecopy( n, n, M2_2n_n.data(),   nx2, E1, n ) ; // memorizzo scambiate
        gecopy( n, n, M2_2n_n.data()+n, nx2, E,  n ) ;

        gezero( n, n,        M1_2n_n.data(),   nx2 ) ;
        gecopy( n, n, F1, n, M1_2n_n.data()+n, nx2 ) ;
        M2_2n_n.noalias() = QR.householderQ().transpose() * M1_2n_n ;
        gecopy( n, n, M2_2n_n.data(),   nx2, F1, n ) ; // memorizzo scambiate
        gecopy( n, n, M2_2n_n.data()+n, nx2, F,  n ) ;
      }
      jump_block *= 2 ;
    }
  }

  #ifdef LU_BABD_QR_USE_THREAD
  template <typename t_Value>
  void
  BabdQR<t_Value>::reduction_mt( integer nth ) {
    valuePointer EE = tmpM+nth*nxnx2 ;
    valuePointer FF = EE+nxn ;

    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;

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
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
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
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
        }
      }
    }
  }
  #endif

  template <typename t_Value>
  void
  BabdQR<t_Value>::factorize( integer           _nblock,
                              integer           _n,
                              integer           _q,
                              valueConstPointer AdAu,
                              valueConstPointer H0,
                              valueConstPointer HN,
                              valueConstPointer Hq ) {

    nblock = _nblock ;
    n      = _n ;
    m      = _n+_q ;

    nx2    = 2*n ;
    nxn    = n*n ;
    nxnx2  = nxn*2 ;
    nm     = n+m ;

    integer nnzADAU = nblock*nxnx2 ;
    baseValue.allocate(nnzADAU) ;
    AdAu_blk = baseValue(nnzADAU) ;

    v1_n    . resize(n) ;
    v2_n    . resize(n) ;
    v1_nx2  . resize(nx2) ;
    v2_nx2  . resize(nx2) ;
    v1_nm   . resize(nm) ;
    v2_nm   . resize(nm) ;
    M1_2n_n . resize(nx2,n) ;
    M2_2n_n . resize(nx2,n) ;
    QR_blk  . resize(nblock) ;

    alglin::copy( nnzADAU, AdAu, 1, AdAu_blk, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef LU_BABD_QR_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    #endif

    jump_block = 1 ;
    #ifdef LU_BABD_QR_USE_THREAD
    to_be_done = usedThread ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &BabdQR<t_Value>::reduction_mt, this, nt ) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #endif
    reduction() ;

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    matType M(nm,nm) ;
    gecopy( n, nx2, AdAu_blk, n, M.data(), nm ) ;
    gecopy( m, n, H0, m, M.data()+n,      nm ) ;
    gecopy( m, n, HN, m, M.data()+n+n*nm, nm ) ;
    if ( _q > 0 ) {
      gezero( n, _q,        M.data()+nx2*nm,   nm ) ;
      gecopy( m, _q, Hq, m, M.data()+nx2*nm+n, nm ) ;
    }
    QR_last_blk.compute(M) ;
  }

  /*
  //  __                                  _                  _
  // / _| ___  _ ____      ____ _ _ __ __| |    _ __ ___  __| |_   _  ___ ___
  // | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |   | '__/ _ \/ _` | | | |/ __/ _ \
  // |  _| (_) | |   \ V  V / (_| | | | (_| |   | | |  __/ (_| | |_| | (_|  __/
  // |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|___|_|  \___|\__,_|\__,_|\___\___|
  //                                       |_____|
  */

  template <typename t_Value>
  void
  BabdQR<t_Value>::forward_reduce( valuePointer y ) const {
    while ( jump_block < nblock ) {

      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        QR_type const & QR = QR_blk[k1] ;
        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;

        copy( n, yk,  1, v1_nx2.data(),   1 ) ;
        copy( n, yk1, 1, v1_nx2.data()+n, 1 ) ;

        /*\
         |   applico Q^T e scambio
        \*/
        v2_nx2.noalias() = QR.householderQ().transpose() * v1_nx2 ;
        copy( n, v2_nx2.data(),   1, yk1, 1 ) ;
        copy( n, v2_nx2.data()+n, 1, yk,  1 ) ;
      }

      // aspetta le altre thread
      jump_block *= 2 ;
    }
  }

  #ifdef LU_BABD_QR_USE_THREAD
  template <typename t_Value>
  void
  BabdQR<t_Value>::forward_reduce_mt( integer nth ) const {
    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer     k1 = k+jump_block ;
        valuePointer G = G_blk    + k1 * nxn ;
        integer * ipiv = ipiv_blk + k1 * n ;

        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;

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
        gemv( NO_TRANSPOSE, n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
        }
      }

    }
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

  template <typename t_Value>
  void
  BabdQR<t_Value>::back_substitute( valuePointer y, integer jump_block_min ) const {
    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer      k1 = k+jump_block ;
        integer      k2 = min(k1+jump_block,nblock) ;
        valuePointer E1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer F1 = E1 + nxn ;

        QR_type const & QR = QR_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;
        
        // io -= M*io1
        gemv( NO_TRANSPOSE,
              n, n,
              -1.0, E1, n,
              yk, 1,
              1, yk1, 1 ) ;
        gemv( NO_TRANSPOSE,
              n, n,
              -1.0, F1, n,
              yk2, 1,
              1, yk1, 1 ) ;
        copy( n, yk1, 1, v1_n.data(), 1 ) ;
        QR.matrixQR()
          .topLeftCorner(n,n)
          .template triangularView<Eigen::Upper>()
          .template solveInPlace<Eigen::OnTheLeft>(v1_n) ;

        v2_n.noalias() = QR.colsPermutation()*v1_n;
        copy( n, v2_n.data(), 1, yk1, 1 ) ;

      }
      jump_block /= 2 ;
    }
  }

  #ifdef LU_BABD_QR_USE_THREAD

  template <typename t_Value>
  void
  BabdQR<t_Value>::back_substitute_mt( integer nth ) const {

    while ( jump_block > 0 ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk1 = y_thread + k1 * n ;
        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = y_thread + ( LU_rows[i] ? k2 : k ) * n ;
          yk1[i] -= dot( n, Adu+i, n, LR, 1 ) ;
        }

        // (LU)^(-1) = U^(-1) L^(-1)
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, LU, n, yk1, 1 ) ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, LU, n, yk1, 1 ) ;

        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      }
      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block /= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
        }
      }
    }
  }
  #endif

  /*             _
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */

  template <typename t_Value>
  void
  BabdQR<t_Value>::solve( valuePointer y ) const {

    #ifdef LU_BABD_QR_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    y_thread          = y ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef LU_BABD_QR_USE_THREAD
    if ( usedThread > 0 ) {
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR<t_Value>::forward_reduce_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #endif
    forward_reduce(y) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! 2 by 2 block linear system solution !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    valuePointer ye = y + nblock * n ;
    copy( n, y,  1, v1_nm.data(),   1 ) ;
    copy( m, ye, 1, v1_nm.data()+n, 1 ) ;
    v2_nm = QR_last_blk.solve(v1_nm) ;
    copy( n, v2_nm.data(),   1, y,  1 ) ;
    copy( m, v2_nm.data()+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block /= 2 ;
    #ifdef LU_BABD_QR_USE_THREAD
    if ( usedThread > 0 ) {
      back_substitute( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR<t_Value>::back_substitute_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
      back_substitute( y, 0 ) ;
    }
    #else
    back_substitute( y, 0 ) ;
    #endif

  }

  template class BabdQR<double> ;
  template class BabdQR<float> ;

}
