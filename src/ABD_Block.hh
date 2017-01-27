/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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

///
/// file: ABD_Block.hh
///

#ifndef ABD_BLOCK_HH
#define ABD_BLOCK_HH

#include "BlockBidiagonal.hh"
#include <iostream>

namespace alglin {

  //! LU decomposition of a ABD matrix
  /*!
   * 
   * \date     May, 2016
   * \version  1.0
   * \note     first release May, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           `enrico.bertolazzi\@unitn.it`
   * 
   */
  /*
      Matrix NNZ structure
        col0-col00
     |    +----+----+                        |
     |    |    |    |                        | <- sizeBlock
     |    +----+----+----+                   |
     |         |    |    |                   |
     |         +----+----+----+              |
     |              |    |    |              |
     |              +----+----+----+         |
     |                   |    |    |         |
     |                   +----+----+---+     |
     |                        |        |     | <-- rowN
     |                        | BOTTOM |     |
     |    +----+              +--------+--+  |
     |    | TOP|                       |  |  | <-- row0
     |    +----+                       +--+  |
                                colN  col00

          col0
       |       |
     / +-------+....+                    \
     | |  TOP  |  F :                    | <-- row0
     | +--+----+----+....+               |
     |    |    |    |  F :               | <- sizeBlock
     |    +----+----+----+....+          |
     |         |    |    |  F :          |
     |         +----+----+----+....+     |
     |              |    |    |  F :     |
     |              +----+----+----+...+ |
     |                   |    |    | E | |
     |                   +----+----+---+ |
     |                        |        | | <-- rowN
     |                        | BOTTOM | |
     \                        +--------+ /
                              |        |
                                 colN
  */
  template <typename t_Value>
  class BlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

  private:

    BlockLU( BlockLU const & ) ;
    BlockLU const & operator = ( BlockLU const & ) ;

    mutable integer nblk ;

    integer * swapR_blks ;
    t_Value * F_mat ;

  public:

    using BlockBidiagonal<valueType>::factorize ;
    using BlockBidiagonal<valueType>::dump_ccoord ;

    explicit BlockLU() {}
    ~BlockLU() {}

    virtual
    void
    allocate( integer /* nblock */,
              integer /* n */,
              integer /* nb */,
              // ----------------------
              integer /* numInitialBC */,
              integer /* numFinalBC */,
              integer /* numCyclicBC */,
              // ----------------------
              integer /* numInitialOMEGA */,
              integer /* numFinalOMEGA */,
              integer /* numCyclicOMEGA */ ) {
    /*
      integer Fnnz = (_nblock-1)*_n*_n ;
      integer innz = _nblock*_n ;
      BlockBidiagonal<t_Value>::allocate(_nblock, _n, _q, _nb, Fnnz, innz ) ;
      swapR_blks = this->baseInteger(size_t(innz)) ;
      F_mat      = this->baseReal(size_t(Fnnz)) ;
    */
    }

    virtual
    void
    factorize() ;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( valuePointer in_out ) const ;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( integer nrhs, valuePointer in_out, integer ldRhs ) const ;

  } ;
}

#endif

///
/// eof: ABD_Block.hh
///

