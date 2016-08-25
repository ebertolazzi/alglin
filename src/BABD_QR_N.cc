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

#include "Alglin_tmpl.hh"
#include "BABD_QR_N.hh"
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

  #ifdef BABD_QR_N_USE_THREAD
  template <typename t_Value, integer n>
  BabdQR_N<t_Value,n>::BabdQR_N( integer nth )
  : numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= BABD_QR_N_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << BABD_QR_N_MAX_THREAD ) ;
  }
  #else
  template <typename t_Value, integer n>
  BabdQR_N<t_Value,n>::BabdQR_N()
  { }
  #endif

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
  BabdQR_N<t_Value,n>::reduction() {
    matType2NN  M_2n_n ;
    matType2N2N M1_2n_2n, M2_2n_2n ;
    M1_2n_2n.setZero() ;

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

        matTypeNN & E  = AdAu_blk[2*k] ;
        matTypeNN & F  = AdAu_blk[2*k+1] ;
        matTypeNN & E1 = AdAu_blk[2*k1] ;
        matTypeNN & F1 = AdAu_blk[2*k1+1]  ;

        QR_type_2NN & QR = QR_blk[k1] ;

        /*\
         |  / F \ = Q / R \
         |  \ E1/     \ 0 /
        \*/
        M_2n_n.template block<n,n>(0,0) = F ;
        M_2n_n.template block<n,n>(n,0) = E1 ;
        QR.compute(M_2n_n) ;

        /*\
         |   / 0 I \ Q^T / F \ = / 0 \
         |   \ I 0 /     \ E1/   \ R /
        \*/
        M1_2n_2n.template block<n,n>(0,0) = E  ;
        //M1_2n_2n.block<n,n>(n,0).setZero() ;
        M1_2n_2n.template block<n,n>(n,n) = F1 ;
        //M1_2n_2n.block<n,n>(0,n).setZero() ;
        M2_2n_2n.noalias() = QR.householderQ().transpose() * M1_2n_2n ;
        E1 = M2_2n_2n.template block<n,n>(0,0) ;
        E  = M2_2n_2n.template block<n,n>(n,0) ;
        F1 = M2_2n_2n.template block<n,n>(0,n) ;
        F  = M2_2n_2n.template block<n,n>(n,n) ;
      }
      jump_block *= 2 ;
    }
  }

  #ifdef BABD_QR_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::reduction_mt( integer nth ) {
    matType2NN  M_2n_n ;
    matType2N2N M1_2n_2n, M2_2n_2n ;
    M1_2n_2n.setZero() ;

    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;

        /*\
         |  +-----+-----+
         |  |  E  |  F  |
         |  +-----+-----+-----+
         |        |  E1 |  F1 |
         |        +-----+-----+
        \*/

        matTypeNN & E  = AdAu_blk[2*k] ;
        matTypeNN & F  = AdAu_blk[2*k+1] ;
        matTypeNN & E1 = AdAu_blk[2*k1] ;
        matTypeNN & F1 = AdAu_blk[2*k1+1]  ;

        QR_type_2NN & QR = QR_blk[k1] ;

        /*\
         |  / F \ = Q / R \
         |  \ E1/     \ 0 /
        \*/
        M_2n_n.template block<n,n>(0,0) = F ;
        M_2n_n.template block<n,n>(n,0) = E1 ;
        QR.compute(M_2n_n) ;

        /*\
         |   / 0 I \ Q^T / F \ = / 0 \
         |   \ I 0 /     \ E1/   \ R /
        \*/
        M1_2n_2n.template block<n,n>(0,0) = E  ;
        //M1_2n_2n.block<n,n>(n,0).setZero() ;
        M1_2n_2n.template block<n,n>(n,n) = F1 ;
        //M1_2n_2n.block<n,n>(0,n).setZero() ;
        M2_2n_2n.noalias() = QR.householderQ().transpose() * M1_2n_2n ;
        E1 = M2_2n_2n.template block<n,n>(0,0) ;
        E  = M2_2n_2n.template block<n,n>(n,0) ;
        F1 = M2_2n_2n.template block<n,n>(0,n) ;
        F  = M2_2n_2n.template block<n,n>(n,n) ;
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

  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::factorize( integer           _nblock,
                                  integer           _q,
                                  valueConstPointer AdAu,
                                  valueConstPointer H0,
                                  valueConstPointer HN,
                                  valueConstPointer Hq ) {
    integer const nx2 = n*2 ;
    integer const nxn = n*n ;
    integer const m   = n+_q ;
    integer const nm  = m+n ;

    nblock = _nblock ;
    q      = _q ;

    AdAu_blk . resize(size_t(2*nblock)) ;
    QR_blk   . resize(size_t(nblock)) ;
    v1_nm    . resize(size_t(nm)) ;
    v2_nm    . resize(size_t(nm)) ;
    M_nm_nm  . resize(nm,nm) ;
    
    for ( integer k = 0 ; k < 2*nblock ; ++k, AdAu += nxn )
      Copy<t_Value,nxn,1,1>::eval( AdAu, AdAu_blk[k].data() ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef BABD_QR_N_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    #endif

    jump_block = 1 ;
    #ifdef BABD_QR_N_USE_THREAD
    to_be_done = usedThread ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &BabdQR_N<t_Value,n>::reduction_mt, this, nt ) ;
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
    M_nm_nm.template block<n,n>(0,0) = AdAu_blk[0] ;
    M_nm_nm.template block<n,n>(0,n) = AdAu_blk[1] ;
    gecopy( m, n, H0, m, M_nm_nm.data()+n,      nm ) ;
    gecopy( m, n, HN, m, M_nm_nm.data()+n+n*nm, nm ) ;
    if ( _q > 0 ) {
      gezero( n, _q,        M_nm_nm.data()+nx2*nm,   nm ) ;
      gecopy( m, _q, Hq, m, M_nm_nm.data()+nx2*nm+n, nm ) ;
    }
    QR_last_blk.compute(M_nm_nm) ;
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
  BabdQR_N<t_Value,n>::forward_reduce( valuePointer y ) const {
    vecType2N v1_nx2, v2_nx2 ;
    while ( jump_block < nblock ) {

      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        QR_type_2NN const & QR = QR_blk[k1] ;
        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;

        /*\
         |   applico Q^T e scambio
        \*/
        Copy<t_Value,n,1,1>::eval( yk,  v1_nx2.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk1, v1_nx2.data()+n ) ;
        v2_nx2.noalias() = QR.householderQ().transpose() * v1_nx2 ;
        Copy<t_Value,n,1,1>::eval( v2_nx2.data(),   yk1 ) ;
        Copy<t_Value,n,1,1>::eval( v2_nx2.data()+n, yk  ) ;
      }

      // aspetta le altre thread
      jump_block *= 2 ;
    }
  }

  #ifdef BABD_QR_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::forward_reduce_mt( integer nth ) const {
    vecType2N v1_nx2, v2_nx2 ;
    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer     k1 = k+jump_block ;

        QR_type_2NN const & QR = QR_blk[k1] ;
        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;

        /*\
         |   applico Q^T e scambio
        \*/
        Copy<t_Value,n,1,1>::eval( yk,  v1_nx2.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk1, v1_nx2.data()+n ) ;
        v2_nx2.noalias() = QR.householderQ().transpose() * v1_nx2 ;
        Copy<t_Value,n,1,1>::eval( v2_nx2.data(),   yk1 ) ;
        Copy<t_Value,n,1,1>::eval( v2_nx2.data()+n, yk  ) ;
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

  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::back_substitute( valuePointer y, integer jump_block_min ) const {
    vecTypeN v0_n, v1_n, v2_n ;
    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        integer k2 = min(k1+jump_block,nblock) ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        Copy<t_Value,n,1,1>::eval( yk,  v0_n.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk1, v1_n.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk2, v2_n.data() ) ;

        matTypeNN   const & E1 = AdAu_blk[2*k1] ;
        matTypeNN   const & F1 = AdAu_blk[2*k1+1] ;
        QR_type_2NN const & QR = QR_blk[k1] ;

        v1_n -= E1*v0_n + F1*v2_n ;
        QR.matrixQR()
          .topLeftCorner(n,n)
          .template triangularView<Eigen::Upper>()
          .template solveInPlace<Eigen::OnTheLeft>(v1_n) ;

        #ifdef BABD_QR_N_USE_PIVOTING
        v2_n.noalias() = QR.colsPermutation()*v1_n;
        Copy<t_Value,n,1,1>::eval( v2_n.data(), yk1 ) ;
        #else
        Copy<t_Value,n,1,1>::eval( v1_n.data(), yk1 ) ;
        #endif
      }
      jump_block /= 2 ;
    }
  }

  #ifdef BABD_QR_N_USE_THREAD

  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::back_substitute_mt( integer nth ) const {
    vecTypeN v0_n, v1_n, v2_n ;
    while ( jump_block > 0 ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        integer k2 = min(k1+jump_block,nblock) ;

        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;
        valuePointer yk2 = y_thread + k2 * n ;

        Copy<t_Value,n,1,1>::eval( yk,  v0_n.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk1, v1_n.data() ) ;
        Copy<t_Value,n,1,1>::eval( yk2, v2_n.data() ) ;

        matTypeNN   const & E1 = AdAu_blk[2*k1] ;
        matTypeNN   const & F1 = AdAu_blk[2*k1+1] ;
        QR_type_2NN const & QR = QR_blk[k1] ;

        v1_n -= E1*v0_n + F1*v2_n ;
        QR.matrixQR()
          .topLeftCorner(n,n)
          .template triangularView<Eigen::Upper>()
          .template solveInPlace<Eigen::OnTheLeft>(v1_n) ;
        #ifdef BABD_QR_N_USE_PIVOTING
        v2_n.noalias() = QR.colsPermutation()*v1_n;
        Copy<t_Value,n,1,1>::eval( v2_n.data(), yk1 ) ;
        #else
        Copy<t_Value,n,1,1>::eval( v1_n.data(), yk1 ) ;
        #endif
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

  template <typename t_Value, integer n>
  void
  BabdQR_N<t_Value,n>::solve( valuePointer y ) const {
  
    integer const m = n+q ;

    #ifdef BABD_QR_N_USE_THREAD
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
    #ifdef BABD_QR_N_USE_THREAD
    if ( usedThread > 0 ) {
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR_N<t_Value,n>::forward_reduce_mt, this, nt ) ;
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
    Copy<t_Value,n,1,1>::eval( y, v1_nm.data() ) ;
    copy( m, ye, 1, v1_nm.data()+n, 1 ) ;
    v2_nm = QR_last_blk.solve(v1_nm) ;
    Copy<t_Value,n,1,1>::eval(v2_nm.data(), y ) ;
    copy( m, v2_nm.data()+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block /= 2 ;
    #ifdef BABD_QR_N_USE_THREAD
    if ( usedThread > 0 ) {
      back_substitute( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR_N<t_Value,n>::back_substitute_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
      back_substitute( y, 0 ) ;
    }
    #else
    back_substitute( y, 0 ) ;
    #endif

  }

  // speedup for windows
  #ifndef ALGLIN_OS_WINDOWS
  template class BabdQR_N<double,2> ;
  template class BabdQR_N<double,3> ;
  template class BabdQR_N<double,4> ;
  template class BabdQR_N<double,5> ;
  template class BabdQR_N<double,6> ;
  template class BabdQR_N<double,7> ;
  template class BabdQR_N<double,8> ;
  template class BabdQR_N<double,9> ;
  template class BabdQR_N<double,10> ;
  /*
  template class BabdQR_N<double,11> ;
  template class BabdQR_N<double,12> ;
  template class BabdQR_N<double,13> ;
  template class BabdQR_N<double,14> ;
  template class BabdQR_N<double,15> ;
  template class BabdQR_N<double,16> ;
  */
  //template class BabdQR_N<float> ;
  #endif
}
