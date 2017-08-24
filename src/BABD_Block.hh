/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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

#ifndef LU_BABD_BLOCK_HH
#define LU_BABD_BLOCK_HH

#include "BlockBidiagonal.hh"
#include <vector>

//! Various LU decomposition classes
namespace alglin {

  /*
  //   ____  ____  _            _    _    _   _
  //  | __ )| __ )| | ___   ___| | _| |  | | | |
  //  |  _ \|  _ \| |/ _ \ / __| |/ / |  | | | |
  //  | |_) | |_) | | (_) | (__|   <| |__| |_| |
  //  |____/|____/|_|\___/ \___|_|\_\_____\___/
  */
  //! LU decomposition of a BABD matrix
  /*!
   * 
   * \date     May 30, 2006
   * \version  1.0
   * \note     first release May 30, 2006
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento<br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   */
  template <typename t_Value>
  class BBlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

  private:

    BBlockLU(BBlockLU const &) ;
    BBlockLU const & operator = (BBlockLU const &) ;

    integer m   ; //!< number final rows (m>=n)
    integer N   ; //!< n * (nblock+1) + q
    integer nnz ; //!< total number of non zeros

    /*!
    //
    //  Matrix structure
    //
    //                  n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n                        n     q
    //  +-----+-----+-----+----.........-----+-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |   
    //  |                                                :     |
    //  :                                                :      > n * nblock
    //  :                                                :     | 
    //  :                                                :     |
    //  :                                                :     |
    //  :                              +-----+-----+-----+     |
    //  :                              | Au  |  0  |  0  |     |
    //  :                        +-----+-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  |  0  | n   |
    //  +-----+-----+---......---+-----+-----+=====+=====+    /
    //  |     |     |            |     |     !     |     !
    //  | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
    //  |     |     |            |     |     !     |     !
    //  +-----+-----+---......---+-----+-----+=====+=====+
    //
    */
    
    /*!
    //
    //  Working block AdH_blk [ size = nblock * ( n * (n+m) ) ]
    //
    //                 n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n                      
    //  +-----+-----+-----+----........+-----+
    //  |  Ad | Ad  | Ad  |               Ad | n
    //  +-----+-----+-----+            +-----+
    //  |     |     |            |     |     !
    //  | H0  |  0  |            |     |  0  ! m
    //  |     |     |            |     |     !
    //  +-----+-----+---......---+-----+-----+
    //
    */
    valuePointer AdH_blk ;

    /*!
    //
    //  Working block Au_blk [ size = nblock * (n*n) ]
    //
    //                  n * nblock + q
    //    ________________^________________________
    //   /                                         \
    //     n     n     n                  n      q
    //  +-----+-----+-----+----.........-----+------+
    //  |  Au | Au  | Au  |               Au | fill | n
    //  +-----+-----+-----+----.........-----+------+
    //
    */
    valuePointer Au_blk ;

    /*!
    //
    //  Working block DD_blk [ size = m * m ]
    //
    //     n     q     (n+q=m)
    //  +=====+=====+
    //  !     |     !
    //  ! HN  | Hq  ! m
    //  !     |     !
    //  +=====+=====+
    //
    */    
    valuePointer DD_blk ;

    /*!
    //
    //  Working block FF_blk [ size = nblock * (n*m) ]
    //
    //     n     q        (n+q=m)
    //  +-----+-----+   \
    //  |fill |fill | n |
    //  +-----+-----+   |
    //  |fill |fill | n |
    //  +-----+-----+   |
    //  |fill |fill | n |
    //  +-----+-----+   :
    //  |           |   :  n * (nblock-1)
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  |fill |fill | n |
    //  +-----+-----+   /
    //
    */
    valuePointer FF_blk ;

    // pivot vector
    integer * ipiv_blk ;

  public:
  
    using BlockBidiagonal<t_Value>::factorize ;
    using BlockBidiagonal<t_Value>::dump_ccoord ;

    explicit BBlockLU() { }

    ~BBlockLU() { }

    virtual
    void
    factorize() ALGLIN_OVERRIDE ;

    //! solve linear system previously factorized
    virtual
    void
    solve( valuePointer in_out ) const ALGLIN_OVERRIDE ;

  } ;

  // explicit instantiation declaration to suppress warnings

  #ifdef __GCC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma GCC diagnostic ignored "-Wweak-template-vtables"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class BBlockLU<float> ;
  extern template class BBlockLU<double> ;

  #ifdef __GCC__
  #pragma GCC diagnostic pop
  #endif
  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

}

#endif
