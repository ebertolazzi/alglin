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
/// file: ABD_Diaz.hh
///

/*!
 * \file ABD_Diaz.hxx
 * \brief Diaz implementation of the solver for almost block diagonal systems.
 *
 * This header declares `alglin::DiazLU`, an ABD solver variant based on the
 * Diaz algorithm for systems with top and bottom blocks. The class shares the
 * same data layout as `BlockBidiagonal`, but uses a dedicated elimination
 * procedure together with internally stored permutations during the solve phase.
 */

namespace alglin {

  /*!
   * \brief ABD solver based on the Diaz factorization.
   *
   * The class is designed for almost block diagonal matrices described by an
   * initial block, a sequence of internal blocks, and a final block. The
   * algorithm uses row and column permutations to keep the decomposition stable
   * and reuses the base-class infrastructure for system loading and sparse
   * export.
   *
   * \date     May, 2016
   * \version  1.0
   * \author   Enrico Bertolazzi
   */
  /*


    Matrix NNZ structure
          col0-col00
     ┌    ┌────┼────┐                       ┐
     │    │    │    │                       │ <- size_block
     │    └────┼────┼────┐                  │
     │         │    │    │                  │
     │         └────┼────┼────┐             │
     │              │    │    │             │
     │              └────┼────┼────┐        │
     │                   │    │    │        │
     │                   └───-┼────┴───┐    │
     │                        │        │    │ <-- rowN
     │                        │ BOTTOM │    │
     │    ┌────┐              └────────┼─┐  │
     │    │TOP │                       │ │  │ <- row0
     │    └────┘                       └─┘  │
     └                        
                              colN  col00

    Matrix NNZ structure
          col0
       │       │
     ┌ ┌───────┐                         ┐
     │ │  TOP  │                         │ <-- row0
     │ └──┬────┼────┐                    │
     │    │    │    │                    │ <- size_block
     │    └────┼────┼────┐               │
     │         │    │    │               │
     │         └────┼────┼────┐          │
     │              │    │    │          │
     │              └────┼────┼────┐     │
     │                   │    │    │     │
     │                   └───-┼────┴───┐ │
     │                        │        │ │ <-- rowN
     │                        │ BOTTOM │ │
     └                        └────────┘
                              │        │
                                 colN

  */

  template <typename t_Value>
  class DiazLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value real_type;

    DiazLU( DiazLU const & ) = delete;
    DiazLU const & operator = ( DiazLU const & ) = delete;

  private:

    integer m_NB{25};

    mutable integer m_nblk{0};

    integer * m_swapRC_blks{nullptr};

    void
    LU_left_right(
      integer nrA,
      integer ncA,
      integer ncL,
      integer ncR,
      t_Value A[], integer ldA,
      integer swapR[]
    );

    void
    LU_top_bottom(
      integer nrT,
      integer nrA,
      integer ncA,
      t_Value A[], integer ldA,
      integer nrB,
      t_Value B[], integer ldB,
      integer swapC[]
    );

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

    using BlockBidiagonal<t_Value>::m_mem_int;

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;
    using BlockBidiagonal<t_Value>::n_x_n;

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

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    using BlockBidiagonal<t_Value>::m_la_factorization;

    //! \brief Builds a Diaz solver that has not been allocated yet.
    explicit constexpr DiazLU() = default;

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
    { UTILS_ERROR0("DiazLU::allocate() not defined!\n"); }

    /*!
     * \brief Allocates the solver for a top/bottom block structure.
     * \param nblock number of internal blocks.
     * \param n size of the square internal block.
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
    ) override {
      integer inv{ nblock*n+(col0+colN-2*n) };
      BlockBidiagonal<t_Value>::allocate_top_bottom(
        nblock, n,
        row0, col0,
        rowN, colN,
        nb, 0, inv
      );
      m_swapRC_blks = m_mem_int( inv );
    }

    //! \brief Factorizes the already loaded system.
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
    solve( integer nrhs, real_type in_out[], integer ldRhs ) const override
    { solve_internal( true, nrhs, in_out, ldRhs ); }

    /*!
     * \brief Solves only the ABD part of the factorization.
     * \param in_out right-hand side on input, solution on output.
     */
    void
    solve_ABD( real_type in_out[] ) const
    { solve_internal( false, in_out ); }

    //! \brief Multi-RHS version of \ref solve_ABD(real_type[]) const.
    void
    solve_ABD( integer nrhs, real_type in_out[], integer ldRhs ) const
    { solve_internal( false, nrhs, in_out, ldRhs ); }

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  extern template class DiazLU<float>;
  extern template class DiazLU<double>;
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Diaz.hh
///
