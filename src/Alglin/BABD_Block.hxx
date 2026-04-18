/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                       |
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
/// file: BABD_Block.hxx
///

/*!
 * \file BABD_Block.hxx
 * \brief Block LU solver for cyclic BABD systems.
 *
 * This header defines `alglin::BBlockLU`, a specialization of
 * `BlockBidiagonal` for block-almost-bidiagonal systems without a top block,
 * but with a lower border and cyclic coupling. The class uses dedicated block
 * workspaces to build and factorize the final reduced system.
 */

//! Various LU decomposition classes
namespace alglin {

  /*
  //   ____  ____  _            _    _    _   _
  //  | __ )| __ )| | ___   ___| | _| |  | | | |
  //  |  _ \|  _ \| |/ _ \ / __| |/ / |  | | | |
  //  | |_) | |_) | | (_) | (__|   <| |__| |_| |
  //  |____/|____/|_|\___/ \___|_|\_\_____\___/
  */
  /*!
   * \brief LU factorization for cyclic BABD systems.
   *
   * This class is the "block LU" counterpart for BABD systems whose structure
   * is already expressed through the base-class conventions. It provides a full
   * allocation/factorization/solve workflow for systems with a lower border and
   * final coupling.
   *
   * \date     May 30, 2006
   * \version  1.0
   * \note     first release May 30, 2006
   * \author   Enrico Bertolazzi
   */
  template <typename t_Value>
  class BBlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value real_type;

  private:

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;
    using BlockBidiagonal<t_Value>::m_extra_bc;
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

    BBlockLU( BBlockLU const & ) = delete;
    BBlockLU const & operator = ( BBlockLU const & ) = delete;

    //!
    //!  Matrix structure
    //!
    //!  \verbatim
    //!
    //!                  n * nblock
    //!   ┌────────────────┴──────────────────┐
    //!     n     n     n                        n     q
    //!  ╔═════╦═════╦═════╦═════........═════╦═════╦═════╗    ─┐
    //!  ║  Ad ║ Au  ║  0  ║                  ║  0  ║  0  ║ n   │
    //!  ╠═════╬═════╬═════╬═            ═════╬═════╬═════╣     │
    //!  ║  0  ║ Ad  ║ Au  ║  0               ║  0  ║  0  ║ n   │
    //!  ╠═════╬═════╬═════╬═════╬═      ═════╬═════╬═════╣     │
    //!  |                                                :     │
    //!  :                                                :     ├─ n * nblock
    //!  :                                                :     │
    //!  :                                                :     │
    //!  :                                                :     │
    //!  :                               ═════╬═════╬═════╣     │
    //!  :                              ║ Au  ║  0  ║  0  ║     │
    //!  :                        ╬═════╬═════╬═════╬═════╣     │
    //!  :                        ║  0  ║ Ad  ║ Au  ║  0  ║ n   │
    //!  ╠═════╬═════╬════.....═══╬═════╬═════╬═════╬═════╣    ─┘
    //!  ║     ║     ║            ║     ║     ║     ║     ║
    //!  ║ H0  ║  0  ║            ║     ║  0  ║ HN  ║ Hq  ║ m
    //!  ║     ║     ║            ║     ║     ║     ║     ║
    //!  ╚═════╩═════╩════.....═══╩═════╩═════╩═════╩═════╝
    //!
    //!  \endverbatim
    //!
    //!  Working block AdH_blk [ size = nblock * ( n * (n+m) ) ]
    //!
    //!  \verbatim
    //!
    //!                  n * nblock
    //!   ┌────────────────┴──────────────────┐
    //!     n     n     n
    //!  ╔═════╦═════╦═════╦════........╦═════╗
    //!  ║  Ad ║ Ad  ║ Ad  ║            ║  Ad ║ n
    //!  ╠═════╬═════╬═════╬═           ╬═════╣
    //!  ║     ║     ║            ║     ║     ║
    //!  ║ H0  ║  0  ║            ║     ║  0  ║ m
    //!  ║     ║     ║            ║     ║     ║
    //!  ╚═════╩═════╩════.....═══╩═════╩═════╝
    //!
    //!  \endverbatim
    //!
    real_type * m_AdH_blk{nullptr};

    //!
    //!
    //!  Working block Au_blk [ size = nblock * (n*n) ]
    //!
    //!  \verbatim
    //!                  n * nblock+ q
    //!  ┌─────────────────┴─────────────────────────┐
    //!     n     n     n                  n      q
    //!  ╔═════╦═════╦═════╦════........╦═════╦══════╗
    //!  ║  Au ║ Au  ║ Au  ║               Au ║ fill ║ n
    //!  ╚═════╩═════╩═════╩════.....═══╩═════╩══════╝
    //!
    //!  \endverbatim
    //!
    real_type * m_Au_blk{nullptr};

    //!
    //!  Working block DD_blk [ size = m * m ]
    //!
    //!  \verbatim
    //!
    //!     n     q     (n+q=m)
    //!  +=====+=====+
    //!  !     |     !
    //!  ! HN  | Hq  ! m
    //!  !     |     !
    //!  +=====+=====+
    //!
    //!  \endverbatim
    //!
    real_type * m_DD_blk{nullptr};

    //!
    //!  Working block FF_blk [ size = nblock * (n*m) ]
    //!
    //!  \verbatim
    //!
    //!     n     q        (n+q=m)
    //!  ┌─────┬─────┐  ─┐
    //!  │fill │fill │ n │
    //!  ├─────┼─────┤   │
    //!  │fill │fill │ n │
    //!  ├─────┼─────┤   │
    //!  │fill │fill │ n │
    //!  ├─────┼─────┤   :
    //!  |           |   :  n * (nblock-1)
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  │fill │fill │ n │
    //!  └─────┴─────┘  ─┘
    //!
    //!  \endverbatim
    //!
    real_type * m_FF_blk{nullptr};

    // pivot vector
    integer * m_ipiv_blk{nullptr};

  public:

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    //! \brief Builds a solver that has not been allocated yet.
    explicit BBlockLU() = default;

    /*!
     * \brief Allocates the solver for the general BABD geometry.
     * \param nblock number of internal blocks.
     * \param n size of the square internal block.
     * \param nb size of the right border.
     * \param num_initial_BC number of initial border rows.
     * \param num_final_BC number of final border rows.
     * \param num_cyclic_BC number of cyclic border rows.
     * \param num_initial_OMEGA number of initial extra columns.
     * \param num_final_OMEGA number of final extra columns.
     * \param num_cyclic_OMEGA number of cyclic extra columns.
     */
    virtual
    void
    allocate(
      integer nblock,
      integer n,
      integer nb,
      integer num_initial_BC,
      integer num_final_BC,
      integer num_cyclic_BC,
      integer num_initial_OMEGA,
      integer num_final_OMEGA,
      integer num_cyclic_OMEGA
    ) override;

    /*!
     * \brief Allocation variant for top/bottom-style interfaces.
     *
     * Specialization `BBlockLU` only supports cases compatible with the cyclic
     * BABD structure; this method makes explicit the same convention used by
     * the other solvers in the module.
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

    //! \brief Factorizes the loaded BABD system.
    virtual
    void
    factorize() override;

    //! \brief Solves the system for a single right-hand side.
    virtual
    void
    solve( real_type in_out[] ) const override;

    //! \brief Solves the system for multiple right-hand sides.
    virtual
    void
    solve(
      integer   nrhs,
      real_type in_out[],
      integer   ldRhs
    ) const override;

  };

  // explicit instantiation declaration to suppress warnings
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  extern template class BBlockLU<float>;
  extern template class BBlockLU<double>;
  #pragma clang diagnostic pop
  #endif

}

///
/// eof: BABD_Block.hxx
///
