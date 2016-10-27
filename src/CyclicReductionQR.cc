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

#include "CyclicReductionQR.hh"
#include <iostream>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wpadded"
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
 |                                +-----+-----+
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
 |              +-----+-----+
 |
 |
\*/

namespace alglin {

  using namespace std ;

  /*\
   |   ____           _ _        ____          _            _   _
   |  / ___|   _  ___| (_) ___  |  _ \ ___  __| |_   _  ___| |_(_) ___  _ __
   | | |  | | | |/ __| | |/ __| | |_) / _ \/ _` | | | |/ __| __| |/ _ \| '_ \
   | | |__| |_| | (__| | | (__  |  _ <  __/ (_| | |_| | (__| |_| | (_) | | | |
   |  \____\__, |\___|_|_|\___| |_| \_\___|\__,_|\__,_|\___|\__|_|\___/|_| |_|
   |       |___/
  \*/

  template <typename QR_type>
  #ifdef CYCLIC_REDUCTION_USE_THREAD
  CyclicReductionQR<QR_type>::CyclicReductionQR( integer nth )
  #else
  CyclicReductionQR<QR_type>::CyclicReductionQR()
  #endif
  : baseValue("CyclicReductionQR_value")
  #ifdef CYCLIC_REDUCTION_USE_THREAD
  , numThread(nth)
  #endif
  , NB(25)
  {
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    ALGLIN_ASSERT( numThread > 0 && numThread <= CYCLIC_REDUCTION_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << CYCLIC_REDUCTION_MAX_THREAD ) ;
    #endif
  }

  template <typename QR_type>
  CyclicReductionQR<QR_type>::~CyclicReductionQR() {
    baseValue.free() ;
  }
  
  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/

  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::allocate( integer _nblock, integer _n ) {

    if ( _nblock == nblock && n == _n ) return ;

    nblock = _nblock ;
    n      = _n ;
    nx2    = 2*n ;
    nxn    = n*n ;
    nxnx2  = nxn*2 ;

    integer mem = nblock*nxnx2+4*nxn+nx2 ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    mem += (4*nxn+nx2)*CYCLIC_REDUCTION_MAX_THREAD ;
    #endif
    baseValue.allocate(size_t(mem)) ;
    AdAu_blk = baseValue(size_t(nblock*nxnx2)) ;
    M_2n_2n  = baseValue(size_t(4*nxn)) ;
    v_nx2    = baseValue(size_t(nx2)) ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    for ( integer nth = 0 ; nth < CYCLIC_REDUCTION_MAX_THREAD ; ++nth ) {
      M_2n_2n_mt[nth] = baseValue(size_t(4*nxn)) ;
      v_nx2_mt[nth]   = baseValue(size_t(nx2)) ;
    }
    #endif

    QR_blk.resize(nblock) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      QR_blk[i] = new QR_type( nx2, n ) ;

  }

  /*\
   |                _
   |   _ __ ___  __| |_   _  ___ ___
   |  | '__/ _ \/ _` | | | |/ __/ _ \
   |  | | |  __/ (_| | |_| | (_|  __/
   |  |_|  \___|\__,_|\__,_|\___\___|
  \*/
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
  CyclicReductionQR<QR_type>::reduce() {
    jump_block = 1 ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    to_be_done = usedThread = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &CyclicReductionQR<QR_type>::reduce_mt, this, nt ) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #endif

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
        QR.load_block( n, n, F,  n, 0, 0 ) ;
        QR.load_block( n, n, E1, n, n, 0 ) ;
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

  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::reduce_mt( integer nth ) {
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
        QR.load_block( n, n, F,  n, 0, 0 ) ;
        QR.load_block( n, n, E1, n, n, 0 ) ;
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

  /*\
   |    __                                  _
   |   / _| ___  _ ____      ____ _ _ __ __| |
   |  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
   |  |  _| (_) | |   \ V  V / (_| | | | (_| |
   |  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
  \*/
  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::forward( valuePointer y ) const {
    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 0 ) {
      y_thread = y ;
      to_be_done = usedThread = numThread ;
      jump_block_max_mt = nblock>>(usedThread-1) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &CyclicReductionQR<QR_type>::forward_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #endif

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

  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::forward_mt( integer nth ) const {
    valuePointer v_tmp = v_nx2_mt[nth] ;
    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;
        QR_type const & QR = *QR_blk[k1] ;
        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;
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

  /*\
   |   _                _                           _
   |  | |__   __ _  ___| | ____      ____ _ _ __ __| |
   |  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
   |  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
   |  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
  \*/
  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::backward( valuePointer y, integer jump_block_min ) const {
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

  #ifdef CYCLIC_REDUCTION_USE_THREAD

  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::backward_mt( integer nth ) const {

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

        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;
        valuePointer yk2 = y_thread + k2 * n ;
        
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

  template <typename QR_type>
  void
  CyclicReductionQR<QR_type>::backward( valuePointer y ) const {
    for ( jump_block = 1 ; jump_block < nblock ; jump_block *= 2 ) {}
    jump_block /= 2 ;
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 0 ) {
      backward( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &CyclicReductionQR<QR_type>::backward_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
      backward( y, 0 ) ;
    }
    #else
    backward( y, 0 ) ;
    #endif
  }

  template class CyclicReductionQR<QR<double> > ;
  template class CyclicReductionQR<QR<float> > ;
  template class CyclicReductionQR<QRP<double> > ;
  template class CyclicReductionQR<QRP<float> > ;

}
