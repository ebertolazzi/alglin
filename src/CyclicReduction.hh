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

#ifndef CYCLIC_REDUCTION_HH
#define CYCLIC_REDUCTION_HH

#include "BlockBidiagonal.hh"
#include <vector>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#endif

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

  template <typename t_Value>
  integer
  LU_2_block( integer n,
              t_Value A[],
              t_Value B[],
              integer ipiv[],
              integer NB ) ;

  template <typename t_Value, integer n>
  integer
  LU_2_block( t_Value A[],
              t_Value B[],
              integer ipiv[] ) ;

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
  template <typename t_Value>
  class CyclicReduction : public BlockBidiagonal<t_Value> {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    std::vector<std::vector<bool> > LU_rows_blk ;

    CyclicReduction(CyclicReduction const &) ;
    CyclicReduction const & operator = (CyclicReduction const &) ;

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
    valuePointer G_blk ;
    valuePointer tmpM ;

    integer * ipiv_blk ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    mutable mutex              mtx0 ;
    mutable condition_variable cond0 ;
    mutable std::thread        threads[CYCLIC_REDUCTION_MAX_THREAD] ;
    mutable integer            to_be_done ;
            integer const      numThread ;
    mutable integer            jump_block_max_mt ;
    mutable integer            usedThread ;
    mutable valuePointer       y_thread ;
    mutable integer            nrhs_thread ;
    mutable integer            ldY_thread ;
    void forward_mt( integer nth ) const ;
    void backward_mt( integer nth ) const ;
    void reduce_mt( integer nth ) ;

    void forward_nrhs_mt( integer nth ) const ;
    void backward_nrhs_mt( integer nth ) const ;
    #endif

    integer NB ; // blocking factor

    #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
    template <integer N>
    class FixedSize {
      CyclicReduction * const pCR ;
    public:
      FixedSize( CyclicReduction * _pCR ) ;
      void reduce() ;
      void forward( valuePointer y ) const ;
      void backward( valuePointer y, integer jump_block_min ) const ;
      #ifdef CYCLIC_REDUCTION_USE_THREAD
      void reduce_mt( integer nth ) ;
      void forward_mt( integer nth ) const ;
      void backward_mt( integer nth ) const ;
      #endif
    } ;
    
    FixedSize<2>  fixed2 ;
    FixedSize<3>  fixed3 ;
    FixedSize<4>  fixed4 ;
    FixedSize<5>  fixed5 ;
    FixedSize<6>  fixed6 ;
    FixedSize<7>  fixed7 ;
    FixedSize<8>  fixed8 ;
    FixedSize<9>  fixed9 ;
    FixedSize<10> fixed10 ;
    #endif

    void backward( valuePointer y, integer jump_block_min ) const ;
    void backward( integer      nrhs,
                   valuePointer y,
                   integer      ldY,
                   integer      jump_block_min ) const ;

    mutable integer jump_block ;

  public:

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    explicit
    CyclicReduction( integer nth = integer(std::thread::hardware_concurrency()) ) ;
    #else
    explicit
    CyclicReduction() ;
    #endif

    ~CyclicReduction() ;

    //! load matrix in the class
    virtual
    void
    allocate( integer nblk, integer n, integer q ) ;

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
    //  :                                                       =         > n * nblock
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
    void forward( integer nrhs, valuePointer rhs, integer ldRhs ) const ;

    /*
    //  Given y1 and yn of the reduced linear system compute y2, y3, ... y(n-1)
    */
    void backward( valuePointer y ) const ;
    void backward( integer nrhs, valuePointer y, integer ldY ) const ;

  } ;
}

#endif
