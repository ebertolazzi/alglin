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

#ifndef CYCLIC_REDUCTION_QR_HH
#define CYCLIC_REDUCTION_QR_HH

#include "Alglin.hh"
#include "Alglin++.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#endif

#include <vector>

#ifdef CYCLIC_REDUCTION_USE_THREAD
  #include <thread>
  #include <mutex>
  #include <condition_variable>
  #include <atomic>
  #ifndef CYCLIC_REDUCTION_MAX_THREAD
    #define CYCLIC_REDUCTION_MAX_THREAD 256
  #endif
#endif

namespace alglin {

  /*\
   |   ____           _ _        ____          _            _   _
   |  / ___|   _  ___| (_) ___  |  _ \ ___  __| |_   _  ___| |_(_) ___  _ __
   | | |  | | | |/ __| | |/ __| | |_) / _ \/ _` | | | |/ __| __| |/ _ \| '_ \
   | | |__| |_| | (__| | | (__  |  _ <  __/ (_| | |_| | (__| |_| | (_) | | | |
   |  \____\__, |\___|_|_|\___| |_| \_\___|\__,_|\__,_|\___|\__|_|\___/|_| |_|
   |       |___/
  \*/

  //! Cyclic reduction of a block bidiagonal matrix
  /*!
   * 
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
   */
  template <typename QR_type>
  class CyclicReductionQR {
  private:

    typedef typename QR_type::valueType         valueType ;
    typedef typename QR_type::valueType*        valuePointer ;
    typedef typename QR_type::valueType const * valueConstPointer ;

    Malloc<valueType> baseValue ;

    CyclicReductionQR(CyclicReductionQR const &) ;
    CyclicReductionQR const & operator = (CyclicReductionQR const &) ;

    integer nblock ; //!< total number of blocks
    integer n      ; //!< size of square blocks

    /*
    //
    //  Matrix structure
    //
    //                 (n+1) * nblock
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n   
    //  +-----+-----+-----+----.........-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |   
    //  |                                                |
    //  :                                                 > n * nblock
    //  :                                                | 
    //  :                                                |
    //  :                                                |
    //  :                              +-----+-----+     |
    //  :                              | Au  |  0  |     |
    //  :                        +-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  | n   |
    //  +-----+-----+---......---+-----+-----+=====+    /
    //
    */

    ///////////////////////////////////////////////////////
    std::vector<QR_type*> QR_blk ;
    valuePointer          AdAu_blk ;

    mutable valuePointer  M_2n_2n, v_nx2 ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    mutable mutex              mtx0 ;
    mutable condition_variable cond0 ;
    mutable std::thread        threads[CYCLIC_REDUCTION_MAX_THREAD] ;
    mutable integer            to_be_done ;
            integer const      numThread ;
    mutable integer            usedThread ;
    mutable integer            jump_block_max_mt ;
    mutable valuePointer       y_thread ;
    
    mutable valuePointer  M_2n_2n_mt[CYCLIC_REDUCTION_MAX_THREAD],
                          v_nx2_mt[CYCLIC_REDUCTION_MAX_THREAD] ;

    void forward_mt( integer nth ) const ;
    void backward_mt( integer nth ) const ;
    void reduce_mt( integer nth ) ;
    #endif

    void backward( valuePointer y, integer jump_block_min ) const ;

    mutable integer jump_block ;

    // some derived constants
    integer nx2 ;
    integer nxn ;
    integer nxnx2 ;

    integer NB ; // blocking factor

  public:

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    explicit CyclicReductionQR( integer nth = integer(std::thread::hardware_concurrency()) ) ;
    #else
    explicit CyclicReductionQR() ;
    #endif

    ~CyclicReductionQR() ;

    //! load matrix in the class
    void allocate( integer nblk, integer n ) ;

    integer getNblock() const { return nblock ; }
    integer getN() const { return n ; }

    void
    loadBlock( integer           nbl,
               valueConstPointer AdAu,
               integer           ldA ) {
      gecopy( n, nx2, AdAu, ldA, AdAu_blk + nbl*nxnx2, n ) ;
    }

    void
    loadBlockLeft( integer           nbl,
                   valueConstPointer Ad,
                   integer           ldA ) {
      gecopy( n, n, Ad, ldA, AdAu_blk + nbl*nxnx2, n ) ;
    }
    
    void
    loadBlockRight( integer           nbl,
                    valueConstPointer Au,
                    integer           ldA ) {
      gecopy( n, n, Au, ldA, AdAu_blk + nbl*nxnx2 + nxn, n ) ;
    }

    valueConstPointer getPointer_LR() const { return AdAu_blk ; }

    void
    getBlock_LR( valuePointer LR, integer ldA ) const
    { gecopy( n, nx2, AdAu_blk, n, LR, ldA ) ; }

    void
    getBlock_L( valuePointer L, integer ldA ) const
    { gecopy( n, n, AdAu_blk, n, L, ldA ) ; }
    
    void
    getBlock_R( valuePointer R, integer ldA ) const
    { gecopy( n, n, AdAu_blk+nxn, n, R, ldA ) ; }

    //! load matrix in the class
    void
    reduce( integer           _nblk,
            integer           _n,
            valueConstPointer AdAu ) {
      allocate( _nblk, _n ) ;
      copy( nblock*nxnx2, AdAu, 1, AdAu_blk, 1 ) ;
      reduce() ;
    }

    /*
    //
    //  Apply reduction to the coeff matric of the linear system
    //
    //                 (n+1) * nblock
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n   
    //  +-----+-----+-----+----.........-----+-----+    +----+     +----+ \
    //  |  Ad | Au  |  0  |                  |  0  | n  | y1 |     | b1 | |
    //  +-----+-----+-----+             -----+-----+    +----+     +----+ |
    //  |  0  | Ad  | Au  |  0               |  0  | n  | y2 |     | b2 | |
    //  +-----+-----+-----+-----+       -----+-----+    +----+     +----+ |
    //  |  0  |  0  | Ad  | Au  |            |  0  | n  | y3 |     | b3 | |
    //  +-----+-----+-----+-----+       -----+-----+    +----+     +----+ |
    //  |                                                                 |
    //  :                                                      =          > n * nblock
    //  :                                                                 |
    //  :                                                                 |
    //  :                                                                 |
    //  :                              +-----+-----+    +----+     +----+ |
    //  :                              | Au  |  0  |    | .. |     | .. | |
    //  :                        +-----+-----+-----+    +----+     +----+ |
    //  :                        |  0  | Ad  | Au  | n  | yn |     | bn | |
    //  +-----+-----+---......---+-----+-----+=====+    +----+     +----+ /
    //
    //  After reduction the linear system reduce to
    //
    //     n     n
    //  +-----+-----+  +----+
    //  |  L  | R   |  | y1 |
    //  +-----+-----+  +----+
    //                 | yn |
    //                 +----+
    */
    void reduce() ;

    /*
    //  Apply reduction to the RHS of the linear system.
    //  After reduction the linear system becomes:
    //
    //     n     n
    //  +-----+-----+  +----+   +----+
    //  |  L  | R   |  | y1 |   | c1 |
    //  +-----+-----+  +----+ = +----+
    //                 | yn |   | cn |
    //                 +----+   +----+
    */
    void forward( valuePointer rhs ) const ;

    /*
    //  Given y1 and yn of the reduced linear system compute y2, y3, ... y(n-1)
    */
    void backward( valuePointer y ) const ;

  } ;
}

#endif
