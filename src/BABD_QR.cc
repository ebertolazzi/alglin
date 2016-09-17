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

  #ifdef BABD_QR_USE_THREAD
  template <typename QR_type>
  BabdQR<QR_type>::BabdQR( integer nth )
  : baseValue("BabdQR_value")
  , numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= BABD_QR_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << BABD_QR_MAX_THREAD ) ;
  }
  #else
  template <typename QR_type>
  BabdQR<QR_type>::BabdQR()
  : baseValue("BabdQR_value")
  { }
  #endif

  template <typename QR_type>
  BabdQR<QR_type>::~BabdQR() {
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
  
  template <typename QR_type>
  void
  BabdQR<QR_type>::reduction() {

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

        QR_type & QR = *QR_blk[k1] ;

        /*\
         |  / F \ = Q / R \
         |  \ E1/     \ 0 /
        \*/
        QR.block_load( n, n, F,  n, 0, 0 ) ;
        QR.block_load( n, n, E1, n, n, 0 ) ;
        QR.factorize() ;

        /*\
         |   / 0 I \ Q^T / F \ = / 0 \
         |   \ I 0 /     \ E1/   \ R /
        \*/
        gecopy( n, n, E, n,  M_2n_2n,         nx2 ) ;
        gezero( n, n,        M_2n_2n+n,       nx2 ) ;
        gezero( n, n,        M_2n_2n+n*nx2,   nx2 ) ;
        gecopy( n, n, F1, n, M_2n_2n+n*nx2+n, nx2 ) ;

        QR.Qt_mul( nx2, nx2, M_2n_2n, nx2 ) ;

        gecopy( n, n, M_2n_2n,         nx2, E1, n ) ; // memorizzo scambiate
        gecopy( n, n, M_2n_2n+n,       nx2, E,  n ) ;
        gecopy( n, n, M_2n_2n+n*nx2,   nx2, F1, n ) ; // memorizzo scambiate
        gecopy( n, n, M_2n_2n+n*nx2+n, nx2, F,  n ) ;
      }
      jump_block *= 2 ;
    }
  }
  
#ifdef BABD_QR_USE_THREAD
  template <typename QR_type>
  void
  BabdQR<QR_type>::reduction_mt( integer nth ) {
    valuePointer M = M_2n_2n_mt[nth] ;
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

        valuePointer E  = AdAu_blk + k  * nxnx2 ;
        valuePointer E1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer F  = E  + nxn ;
        valuePointer F1 = E1 + nxn ;

        QR_type & QR = *QR_blk[k1] ;

        /*\
         |  / F \ = Q / R \
         |  \ E1/     \ 0 /
        \*/
        QR.block_load( n, n, F,  n, 0, 0 ) ;
        QR.block_load( n, n, E1, n, n, 0 ) ;
        QR.factorize() ;

        /*\
         |   / 0 I \ Q^T / F \ = / 0 \
         |   \ I 0 /     \ E1/   \ R /
        \*/
        gecopy( n, n, E, n,  M,         nx2 ) ;
        gezero( n, n,        M+n,       nx2 ) ;
        gezero( n, n,        M+n*nx2,   nx2 ) ;
        gecopy( n, n, F1, n, M+n*nx2+n, nx2 ) ;

        QR.Qt_mul( nx2, nx2, M, nx2 ) ;

        gecopy( n, n, M,         nx2, E1, n ) ; // memorizzo scambiate
        gecopy( n, n, M+n,       nx2, E,  n ) ;
        gecopy( n, n, M+n*nx2,   nx2, F1, n ) ; // memorizzo scambiate
        gecopy( n, n, M+n*nx2+n, nx2, F,  n ) ;

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

  template <typename QR_type>
  void
  BabdQR<QR_type>::allocate( integer _nblock, integer _n, integer _q ) {

    if ( _nblock == nblock && n == _n && _q == m-n ) return ;

    nblock = _nblock ;
    n      = _n ;
    m      = _n+_q ;

    nx2    = 2*n ;
    nxn    = n*n ;
    nxnx2  = nxn*2 ;
    nm     = n+m ;

    integer mem = nblock*nxnx2+4*nxn+nx2+nm ;
    #ifdef BABD_QR_USE_THREAD
    mem += (4*nxn+nx2+nm)* BABD_QR_MAX_THREAD ;
    #endif
    baseValue.allocate(size_t(mem)) ;
    AdAu_blk = baseValue(size_t(nblock*nxnx2)) ;
    M_2n_2n  = baseValue(size_t(4*nxn)) ;
    v_nx2    = baseValue(size_t(nx2)) ;
    v_nm     = baseValue(size_t(nm)) ;

    #ifdef BABD_QR_USE_THREAD
    for ( integer nth = 0 ; nth < BABD_QR_MAX_THREAD ; ++nth ) {
      M_2n_2n_mt[nth]  = baseValue(size_t(4*nxn)) ;
      v_nx2_mt[nth]    = baseValue(size_t(nx2)) ;
      v_nm_mt[nth]     = baseValue(size_t(nm)) ;
    }
    #endif

    QR_blk.resize(nblock) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      QR_blk[i] = new QR_type( nx2, n ) ;

    QR_last_blk.allocate(nm,nm) ;

  }

  template <typename QR_type>
  void
  BabdQR<QR_type>::factorize( integer           _nblock,
                              integer           _n,
                              integer           _q,
                              valueConstPointer AdAu,
                              valueConstPointer H0,
                              valueConstPointer HN,
                              valueConstPointer Hq ) {

    allocate( _nblock, _n, _q ) ;
    alglin::copy( nblock*nxnx2, AdAu, 1, AdAu_blk, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef BABD_QR_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    #endif

    jump_block = 1 ;
    #ifdef BABD_QR_USE_THREAD
    to_be_done = usedThread ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &BabdQR<QR_type>::reduction_mt, this, nt ) ;
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
    QR_last_blk.block_load( n, nx2, AdAu_blk, n, 0, 0 ) ;
    QR_last_blk.block_load( m, n,         H0, m, n, 0 ) ;
    QR_last_blk.block_load( m, n,         HN, m, n, n ) ;
    if ( _q > 0 ) {
      QR_last_blk.block_zero( n, _q,        0, nx2 ) ;
      QR_last_blk.block_load( m, _q, Hq, m, n, nx2 ) ;
    }
    QR_last_blk.factorize() ;
  }

  /*
  //  __                                  _                  _
  // / _| ___  _ ____      ____ _ _ __ __| |    _ __ ___  __| |_   _  ___ ___
  // | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |   | '__/ _ \/ _` | | | |/ __/ _ \
  // |  _| (_) | |   \ V  V / (_| | | | (_| |   | | |  __/ (_| | |_| | (_|  __/
  // |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|___|_|  \___|\__,_|\__,_|\___\___|
  //                                       |_____|
  */

  template <typename QR_type>
  void
  BabdQR<QR_type>::forward_reduce( valuePointer y ) const {
    while ( jump_block < nblock ) {

      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        QR_type const & QR = *QR_blk[k1] ;
        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        /*\
         |   applico Q^T e scambio
        \*/
        copy( n, yk,  1, v_nx2,   1 ) ;
        copy( n, yk1, 1, v_nx2+n, 1 ) ;
        QR.Qt_mul( v_nx2 ) ;
        copy( n, v_nx2,   1, yk1, 1 ) ;
        copy( n, v_nx2+n, 1, yk,  1 ) ;
      }

      // aspetta le altre thread
      jump_block *= 2 ;
    }
  }

  #ifdef BABD_QR_USE_THREAD
  template <typename QR_type>
  void
  BabdQR<QR_type>::forward_reduce_mt( integer nth, valuePointer y ) const {
    valuePointer v_tmp = v_nx2_mt[nth] ;
    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        QR_type const & QR = *QR_blk[k1] ;
        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        /*\
         |   applico Q^T e scambio
        \*/
        copy( n, yk,  1, v_tmp,   1 ) ;
        copy( n, yk1, 1, v_tmp+n, 1 ) ;
        QR.Qt_mul( v_tmp ) ;
        copy( n, v_tmp,   1, yk1, 1 ) ;
        copy( n, v_tmp+n, 1, yk,  1 ) ;
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

  template <typename QR_type>
  void
  BabdQR<QR_type>::back_substitute( valuePointer y, integer jump_block_min ) const {
    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer      k1 = k+jump_block ;
        integer      k2 = min(k1+jump_block,nblock) ;
        valuePointer E1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer F1 = E1 + nxn ;

        QR_type const & QR = *QR_blk[k1] ;

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
        QR.invR_mul( yk1 ) ;
        QR.permute( yk1 ) ;
      }
      jump_block /= 2 ;
    }
  }

  #ifdef BABD_QR_USE_THREAD
  template <typename QR_type>
  void
  BabdQR<QR_type>::back_substitute_mt( integer nth, valuePointer y, integer jump_block_min ) const {

    while ( jump_block > 0 ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer      k1 = k+jump_block ;
        integer      k2 = min(k1+jump_block,nblock) ;
        valuePointer E1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer F1 = E1 + nxn ;

        QR_type const & QR = *QR_blk[k1] ;

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
        QR.invR_mul( yk1 ) ;
        QR.permute( yk1 ) ;
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

  template <typename QR_type>
  void
  BabdQR<QR_type>::solve( valuePointer y ) const {

    #ifdef BABD_QR_USE_THREAD
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
    #ifdef BABD_QR_USE_THREAD
    if ( usedThread > 0 ) {
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR<QR_type>::forward_reduce_mt, this, nt, y ) ;
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
    copy( n, y,  1, v_nm,   1 ) ;
    copy( m, ye, 1, v_nm+n, 1 ) ;
    QR_last_blk.solve(v_nm) ;
    copy( n, v_nm,   1, y,  1 ) ;
    copy( m, v_nm+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block /= 2 ;
    #ifdef BABD_QR_USE_THREAD
    if ( usedThread > 0 ) {
      back_substitute( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &BabdQR<QR_type>::back_substitute_mt, this, nt, y, jump_block_max_mt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
      back_substitute( y, 0 ) ;
    }
    #else
    back_substitute( y, 0 ) ;
    #endif

  }

  template class BabdQR<QR<double> > ;
  template class BabdQR<QR<float> > ;
  template class BabdQR<QRP<double> > ;
  template class BabdQR<QRP<float> > ;

}
