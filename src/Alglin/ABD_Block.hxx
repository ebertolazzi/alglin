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
/// file: ABD_Block.hxx
///

namespace alglin {

  //! LU decomposition of a ABD matrix
  /*!
   *
   * \date     January, 2017
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
  \*/
  /*\
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
  \*/
  template <typename t_Value>
  class BlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value valueType;

  private:

    BlockLU( BlockLU const & );
    BlockLU const & operator = ( BlockLU const & );

    mutable integer nblk;

    integer * swap0;
    integer * swapR_blks;
    integer   Work_lda, F_size, F_lda;
    t_Value * Work_mat;
    t_Value * Work_mat1;
    t_Value * F_mat;

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

    explicit BlockLU() {}
    // ~BlockLU() UTILS_OVERRIDE {}

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
    ) UTILS_OVERRIDE
    { UTILS_ERROR0("BlockLU::allocate() not defined!\n"); }

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
    ) UTILS_OVERRIDE;

    virtual
    void
    factorize() UTILS_OVERRIDE;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( valueType in_out[] ) const UTILS_OVERRIDE
    { solve_internal( true, in_out ); }

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve(
      integer   nrhs,
      valueType in_out[],
      integer   ldRhs
    ) const UTILS_OVERRIDE {
      this->solve_internal( true, nrhs, in_out, ldRhs );
    }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD( valueType in_out[] ) const
    { this->solve_internal( false, in_out ); }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD(
      integer   nrhs,
      valueType in_out[],
      integer   ldRhs
    ) const {
      this->solve_internal( false, nrhs, in_out, ldRhs );
    }

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class BlockLU<float>;
  extern template class BlockLU<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Block.hxx
///
