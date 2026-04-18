/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: ABD_Block.hxx
///

/*!
 * \file ABD_Block.hxx
 * \brief Block LU solver for almost block diagonal (ABD) systems.
 *
 * This header defines `alglin::BlockLU`, a specialized implementation for ABD
 * matrices described by a top block, a sequence of internal bidiagonal blocks,
 * and a final bottom block. The class inherits the loading and sparse export
 * machinery from `BlockBidiagonal` and adds LU factorization together with
 * solvers for both the full bordered system and the ABD part alone.
 */

namespace alglin {

  /*!
   * \brief LU factorization for ABD systems with top/bottom blocks.
   *
   * The class assumes the usual top/bottom geometry of almost block diagonal
   * systems. After the blocks are loaded through the base-class API,
   * `factorize()` builds the internal decomposition and `solve()` handles the
   * full bordered system.
   *
   * The `solve_ABD()` overloads instead apply the decomposition to the ABD part
   * only, without the permutation associated with the bordered formulation.
   *
   * \date     January, 2017
   * \version  1.0
   * \author   Enrico Bertolazzi
   */
  /*\
      Matrix NNZ structure
        col0-col00
     |    +----+----+                        |
     |    |    |    |                        | <- size_block
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
     |    |    |    |  F :               | <- size_block
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

    typedef t_Value real_type;

    BlockLU( BlockLU const & ) = delete;
    BlockLU const & operator = ( BlockLU const & ) = delete;

  private:

    mutable integer m_nblk{0};

    integer * m_swap0{nullptr};
    integer * m_swapR_blks{nullptr};
    integer   m_Work_ldim{0};
    integer   m_F_size{0};
    integer   m_F_ldim{0};
    t_Value * m_Work_mat{nullptr};
    t_Value * m_Work_mat1{nullptr};
    t_Value * m_F_mat{nullptr};

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal( bool do_permute, real_type in_out[] ) const;

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal(
      bool      do_permute,
      integer   nrhs,
      real_type in_out[],
      integer   ldRhs
    ) const;

  public:

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;
    using BlockBidiagonal<t_Value>::n_x_n;

    using BlockBidiagonal<t_Value>::m_mem;
    using BlockBidiagonal<t_Value>::m_mem_int;

    using BlockBidiagonal<t_Value>::m_num_initial_BC;
    using BlockBidiagonal<t_Value>::m_num_final_BC;
    using BlockBidiagonal<t_Value>::m_num_cyclic_BC;
    using BlockBidiagonal<t_Value>::m_num_initial_OMEGA;
    using BlockBidiagonal<t_Value>::m_num_final_OMEGA;
    using BlockBidiagonal<t_Value>::m_num_cyclic_OMEGA;

    using BlockBidiagonal<t_Value>::m_DE_blk;
    using BlockBidiagonal<t_Value>::m_H0Nq;
    using BlockBidiagonal<t_Value>::m_block0;
    using BlockBidiagonal<t_Value>::m_blockN;
    using BlockBidiagonal<t_Value>::m_Bmat;
    using BlockBidiagonal<t_Value>::m_Cmat;
    using BlockBidiagonal<t_Value>::m_Dmat;

    using BlockBidiagonal<t_Value>::m_la_matrix;
    using BlockBidiagonal<t_Value>::m_la_factorization;

    //! \brief Builds a solver that has not been allocated yet.
    explicit BlockLU() {}
    // ~BlockLU() override {}

    //! \brief Not supported by this specialization: use \ref allocate_top_bottom.
    virtual
    void
    allocate(
      integer /* nblock */,
      integer /* n      */,
      integer /* nb     */,
      // ----------------------
      integer /* num_initial_BC */,
      integer /* num_final_BC   */,
      integer /* num_cyclic_BC  */,
      // ----------------------
      integer /* num_initial_OMEGA */,
      integer /* num_final_OMEGA   */,
      integer /* num_cyclic_OMEGA  */
    ) override
    { UTILS_ERROR0("BlockLU::allocate() not defined!\n"); }

    /*!
     * \brief Allocates the solver for an ABD structure with top and bottom blocks.
     * \param nblock number of internal blocks.
     * \param n size of the square internal blocks.
     * \param row0 number of rows in the top block.
     * \param col0 number of columns in the top block.
     * \param rowN number of rows in the bottom block.
     * \param colN number of columns in the bottom block.
     * \param nb size of the right border, if present.
     */
    virtual
    void
    allocate_top_bottom(
      integer nblock,
      integer n,
      integer row0,
      integer col0,
      integer rowN,
      integer colN,
      integer nb
    ) override;

    //! \brief Factorizes the already loaded ABD/bordered system.
    virtual
    void
    factorize() override;

    //! \brief Solves the bordered system for a single right-hand side.
    virtual
    void
    solve( real_type in_out[] ) const override
    { solve_internal( true, in_out ); }

    //! \brief Solves the bordered system for multiple right-hand sides.
    virtual
    void
    solve(
      integer   nrhs,
      real_type in_out[],
      integer   ldRhs
    ) const override {
      this->solve_internal( true, nrhs, in_out, ldRhs );
    }

    /*!
     * \brief Solves only the ABD part of the factorization.
     * \param in_out right-hand side on input, solution on output.
     *
     * This method avoids the permutation induced by the bordered formulation and
     * is useful when the decomposition must be reused on the ABD structure alone.
     */
    void
    solve_ABD( real_type in_out[] ) const
    { this->solve_internal( false, in_out ); }

    //! \brief Multi-RHS version of \ref solve_ABD(real_type[]) const.
    void
    solve_ABD(
      integer   nrhs,
      real_type in_out[],
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
  extern template class BlockLU<float>;
  extern template class BlockLU<double>;
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Block.hxx
///
