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

#ifndef LU_BABD_BLOCK_HH
#define LU_BABD_BLOCK_HH

#include "BlockBidiagonal.hh"
#include <vector>

//! Various LU decomposition classes
namespace alglin {

  /*
  //   ____  _            _    _    _   _ 
  //  | __ )| | ___   ___| | _| |  | | | |
  //  |  _ \| |/ _ \ / __| |/ / |  | | | |
  //  | |_) | | (_) | (__|   <| |__| |_| |
  //  |____/|_|\___/ \___|_|\_\_____\___/ 
  //
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
  class BlockLU : public BlockBidiagonal<t_Value> {

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    BlockLU(BlockLU const &) ;
    BlockLU const & operator = (BlockLU const &) ;

    integer m   ; //!< number final rows (m>=n)
    integer N   ; //!< n * (nblock+1) + q
    integer nnz ; //!< total number of non zeros
    integer __padding ;

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

    explicit BlockLU() { }

    ~BlockLU() { }

    virtual
    void
    allocate( integer nblock, integer n, integer q ) ;

    virtual
    void
    factorize() ;

    //! solve linear system previously factorized
    virtual
    void
    solve( valuePointer in_out ) const ;

    //! load BABD linear system to class
    void
    factorize( integer           nblk,
               integer           n,
               integer           q,
               valueConstPointer AdAu,
               valueConstPointer H0,
               valueConstPointer HN,
               valueConstPointer Hq ) {
      this->allocate( nblk, n, q );
      this->loadBlocks( AdAu, n ) ;
      this->loadBottom( H0, n+q, HN, n+q, Hq, n+q ) ;
      factorize() ;
    }

  } ;

}

#endif
