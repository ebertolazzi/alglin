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

    mutable integer m_nblk;

    integer * m_swap0;
    integer * m_swapR_blks;
    integer   m_Work_lda;
    integer   m_F_size;
    integer   m_F_lda;
    t_Value * m_Work_mat;
    t_Value * m_Work_mat1;
    t_Value * m_F_mat;

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

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    using BlockBidiagonal<t_Value>::m_n;
    using BlockBidiagonal<t_Value>::m_nxnx2;
    using BlockBidiagonal<t_Value>::m_nxn;
    using BlockBidiagonal<t_Value>::m_nblock;

    using BlockBidiagonal<t_Value>::m_baseInteger;
    using BlockBidiagonal<t_Value>::m_baseValue;

    using BlockBidiagonal<t_Value>::m_numInitialBC;
    using BlockBidiagonal<t_Value>::m_numFinalBC;
    using BlockBidiagonal<t_Value>::m_numCyclicBC;
    using BlockBidiagonal<t_Value>::m_numInitialOMEGA;
    using BlockBidiagonal<t_Value>::m_numFinalOMEGA;
    using BlockBidiagonal<t_Value>::m_numCyclicOMEGA;

    using BlockBidiagonal<t_Value>::m_DE_blk;
    using BlockBidiagonal<t_Value>::m_H0Nq;
    using BlockBidiagonal<t_Value>::m_block0;
    using BlockBidiagonal<t_Value>::m_blockN;
    using BlockBidiagonal<t_Value>::m_Bmat;
    using BlockBidiagonal<t_Value>::m_Cmat;
    using BlockBidiagonal<t_Value>::m_Dmat;

    using BlockBidiagonal<t_Value>::m_la_matrix;
    using BlockBidiagonal<t_Value>::m_la_factorization;

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

