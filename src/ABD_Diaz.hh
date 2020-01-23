/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
/// file: ABD_Diaz.hh
///

#pragma once

#ifndef ABD_DIAZ_HH
#define ABD_DIAZ_HH

#include "BlockBidiagonal.hh"
#include <iostream>

namespace alglin {

  //! LU decomposition of a ABD matrix
  /*!
   * 
   * \date     May, 2016
   * \version  1.0
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

      Matrix NNZ internal structure
          col0
       |       |
     / +-------+                         \
     | |  TOP  |                         | <-- row0
     | +--+----+----+                    |
     |    |    |    |                    | <- sizeBlock
     |    +----+----+----+               |
     |         |    |    |               |
     |         +----+----+----+          |
     |              |    |    |          |
     |              +----+----+----+     |
     |                   |    |    |     |
     |                   +----+----+---+ |
     |                        |        | | <-- rowN
     |                        | BOTTOM | |
     \                        +--------+ /
                              |        |
                                 colN
  */

  template <typename t_Value>
  class DiazLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value valueType;

  private:

    DiazLU( DiazLU const & );
    DiazLU const & operator = ( DiazLU const & );

    integer NB;

    mutable integer nblk;

    integer * swapRC_blks;

    void
    LU_left_right(
      integer nrA,
      integer ncA,
      integer ncL,
      integer ncR,
      t_Value * A, integer ldA,
      integer swapR[]
    );

    void
    LU_top_bottom(
      integer nrT,
      integer nrA,
      integer ncA,
      t_Value * A, integer ldA,
      integer nrB,
      t_Value * B, integer ldB,
      integer swapC[]
    );

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal( bool do_permute, valueType in_out[] ) const;

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal(
      bool      do_permute,
      integer   nrhs,
      valueType in_out[],
      integer   ldRhs
    ) const;

  public:

    using BlockBidiagonal<valueType>::factorize;
    using BlockBidiagonal<valueType>::dump_ccoord;

    explicit ALGLIN_CONSTEXPR DiazLU() : NB(25) {}

    virtual
    ~DiazLU() ALGLIN_OVERRIDE
    {}

    virtual
    void
    allocate(
      integer /* nblock */,
      integer /* n      */,
      integer /* nb     */,
      // ----------------------
      integer /* numInitialBC */,
      integer /* numFinalBC   */,
      integer /* numCyclicBC  */,
      // ----------------------
      integer /* numInitialOMEGA */,
      integer /* numFinalOMEGA   */,
      integer /* numCyclicOMEGA  */
    ) ALGLIN_OVERRIDE
    { LW_ERROR0("DiazLU::allocate() not defined!\n"); }

    virtual
    void
    allocateTopBottom(
      integer _nblock,
      integer _n,
      integer _row0,
      integer _col0,
      integer _rowN,
      integer _colN,
      integer _nb
    ) ALGLIN_OVERRIDE {
      integer inv = _nblock*_n+(_col0+_colN-2*_n);
      BlockBidiagonal<t_Value>::allocateTopBottom(
        _nblock, _n,
        _row0, _col0,
        _rowN, _colN,
        _nb, 0, inv
      );
      swapRC_blks = this->baseInteger(size_t(inv));
    }

    virtual
    void
    factorize() ALGLIN_OVERRIDE;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( valueType in_out[] ) const ALGLIN_OVERRIDE
    { solve_internal( true, in_out ); }

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( integer nrhs, valueType in_out[], integer ldRhs ) const ALGLIN_OVERRIDE
    { solve_internal( true, nrhs, in_out, ldRhs ); }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD( valueType in_out[] ) const
    { solve_internal( false, in_out ); }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD( integer nrhs, valueType in_out[], integer ldRhs ) const
    { solve_internal( false, nrhs, in_out, ldRhs ); }

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class DiazLU<float>;
  extern template class DiazLU<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

#endif

///
/// eof: ABD_Diaz.hh
///

