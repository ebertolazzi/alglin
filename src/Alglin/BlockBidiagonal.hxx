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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BlockBidiagonal.hxx
///

/*!
 * \file BlockBidiagonal.hxx
 * \brief Abstract base for block bidiagonal and bordered systems.
 *
 * This header collects the common base class used by the ABD/BABD solvers in
 * the library. The class manages:
 * - allocation of the block data structure;
 * - loading of bidiagonal blocks and border terms;
 * - selection of the solver applied to the last block or reduced system;
 * - sparse export and matrix-vector products.
 *
 * Derived classes implement the actual numerical factorization and solve
 * strategy.
 */

namespace alglin {

  /*\
   |  ____  _            _      ____  _     _ _                               _
   | | __ )| | ___   ___| | __ | __ )(_) __| (_) __ _  __ _  ___  _ __   __ _| |
   | |  _ \| |/ _ \ / __| |/ / |  _ \| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
   | | |_) | | (_) | (__|   <  | |_) | | (_| | | (_| | (_| | (_) | | | | (_| | |
   | |____/|_|\___/ \___|_|\_\ |____/|_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
   |                                                  |___/
  \*/

  /*!
   * \brief Base class for block bidiagonal matrices with an optional border.
   *
   * The class encapsulates problem geometry and the buffers representing:
   * - the internal diagonal/bidiagonal structure;
   * - the right and lower border blocks;
   * - the final reduced block produced by elimination.
   *
   * Derived implementations must provide methods \ref factorize and \ref solve,
   * but can reuse the whole loading, sparse export, and final-solver management
   * infrastructure.
   *
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   * \author   Enrico Bertolazzi
   */
  template <typename t_Value>
  class BlockBidiagonal {
  public:

    typedef t_Value real_type;

    //! \brief Solvers available for the final reduced block.
    using BB_LASTBLOCK_Choice = enum class BB_LASTBLOCK_Choice : integer {
      LU   = 0, //!< Standard LU factorization.
      LUPQ = 1, //!< LU factorization with complete pivoting.
      QR   = 2, //!< QR factorization.
      QRP  = 3, //!< QR factorization with pivoting.
      SVD  = 4, //!< Singular value decomposition.
      LSS  = 5, //!< Least-squares standard.
      LSY  = 6, //!< Least-squares with regularization/symmetry support.
      PINV = 7  //!< Pseudoinverse.
    };

    static
    string
    //! \brief Converts a final-solver choice into a readable string.
    //! \param c choice to convert.
    //! \return Solver description.
    LastBlock_to_string( BB_LASTBLOCK_Choice c ) {
      switch ( c ) {
        case BB_LASTBLOCK_Choice::LU:   return "last block LU";
        case BB_LASTBLOCK_Choice::LUPQ: return "last block LUPQ";
        case BB_LASTBLOCK_Choice::QR:   return "last block QR";
        case BB_LASTBLOCK_Choice::QRP:  return "last block QRP";
        case BB_LASTBLOCK_Choice::SVD:  return "last block SVD";
        case BB_LASTBLOCK_Choice::LSS:  return "last block LSS";
        case BB_LASTBLOCK_Choice::LSY:  return "last block LSY";
        case BB_LASTBLOCK_Choice::PINV: return "last block PINV";
      }
      return "last block not selected";
    }

  public:

    BlockBidiagonal( BlockBidiagonal const & ) = delete;
    BlockBidiagonal const & operator = ( BlockBidiagonal const & ) = delete;

  protected:

    Malloc<real_type> m_mem{"BlockBidiagonal::m_mem"};
    Malloc<integer>   m_mem_int{"BlockBidiagonal::m_mem_int"};

    integer m_number_of_blocks{0}; //!< total number of blocks
    integer m_block_size{0};       //!< size of square blocks
    integer m_extra_bc{0};         //!< extra BC
    integer m_border_size{0};      //!< border size
    integer m_num_equations{0};

    // some derived constants
    integer n_x_n{0};
    integer n_x_nb{0};

    integer m_num_initial_BC{0};
    integer m_num_final_BC{0};
    integer m_num_cyclic_BC{0};
    integer m_num_initial_OMEGA{0};
    integer m_num_final_OMEGA{0};
    integer m_num_cyclic_OMEGA{0};

    Matrix<t_Value>               m_la_matrix;
    Matrix<t_Value>               m_bb_matrix;
    LinearSystemSolver<t_Value> * m_la_factorization{nullptr};
    LinearSystemSolver<t_Value> * m_bb_factorization{nullptr};

    /*
    //
    //  Matrix structure
    //
    //                 (n+1) * nblock
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n
    //  ╒═════╦═════╦═════╦----.........═════╦═════╗    \
    //  ║  D  ║  E  ║  0  ║                  ║  0  ║ n   |
    //  ╠═════╬═════╬═════╬             ═════╬═════╣     |
    //  ║  0  ║  D  ║  E  ║  0               ║  0  ║ n   |
    //  ╠═════╬═════╬═════╬═════╬       ═════╬═════╣     |
    //  ║  0  ║  0  ║  D  ║  E  ║            ║  0  ║ n   |
    //  ╠═════╬═════╬═════╬═════╝       ═════╬═════╣     |
    //  ║                                                |
    //  :                                                 > n * nblock
    //  :                                                |
    //  :                                                |
    //  :                                                |
    //  :                              ╔═════╬═════╣     |
    //  :                              ║  E  ║  0  ║     |
    //  :                        ╔═════╬═════╬═════╣     |
    //  :                        ║  0  ║  D  ║  E  ║ n   |
    //  ╠═════╬═════+═══......═══╬═════╬═════╬═════╬══╣  /
    //  ║     ║                              ║     ║  ║  \
    //  ║ H0  ║                              ║ HN  ║Hq║  |
    //  ║     ║                              ║     ║  ║  | n+q
    //  ╚═════╩═════╩═══......═══╩═════╩═════╩═════╩══╝  /
    //                                               q
    //
    //  Bordered matrix
    //  / A  B \
    //  \ C  D /
    //
    */

    real_type * m_DE_blk { nullptr };
    real_type * m_H0Nq   { nullptr };
    real_type * m_block0 { nullptr };
    real_type * m_blockN { nullptr };
    real_type * m_Bmat   { nullptr };
    real_type * m_Cmat   { nullptr };
    real_type * m_Dmat   { nullptr };

  private:

    LU<real_type>   m_la_lu;
    LUPQ<real_type> m_la_lupq;
    QR<real_type>   m_la_qr;
    QRP<real_type>  m_la_qrp;
    SVD<real_type>  m_la_svd;
    LSS<real_type>  m_la_lss;
    LSY<real_type>  m_la_lsy;
    PINV<real_type> m_la_pinv;

    LU<real_type>   m_bb_lu;
    LUPQ<real_type> m_bb_lupq;
    QR<real_type>   m_bb_qr;
    QRP<real_type>  m_bb_qrp;
    SVD<real_type>  m_bb_svd;
    LSS<real_type>  m_bb_lss;
    LSY<real_type>  m_bb_lsy;
    PINV<real_type> m_bb_pinv;

  public:

    explicit
    BlockBidiagonal()
    : m_la_factorization(&m_la_lu)
    , m_bb_factorization(&m_bb_lu)
    {}

    virtual ~BlockBidiagonal() = default;

    /*!
     * \brief Allocates or resizes the problem using the general structure.
     *
     * This method initializes system geometry, border-block sizes, and any
     * extra storage required by the derived algorithm.
     *
     * \param nblock number of internal bidiagonal blocks.
     * \param n size of each square internal block.
     * \param nb size of the right border.
     * \param num_initial_BC number of rows in the initial lower block.
     * \param num_final_BC number of rows in the final lower block.
     * \param num_cyclic_BC number of rows in the cyclic lower block.
     * \param num_initial_OMEGA number of extra columns associated with the first block.
     * \param num_final_OMEGA number of extra columns associated with the last block.
     * \param num_cyclic_OMEGA number of cyclic extra columns.
     * \param num_extra_r extra real storage requested by the derived class.
     * \param num_extra_i extra integer storage requested by the derived class.
     */
    void
    allocate(
      integer nblock,
      integer n,
      integer nb,
      // ----------------------
      integer num_initial_BC,
      integer num_final_BC,
      integer num_cyclic_BC,
      // ----------------------
      integer num_initial_OMEGA,
      integer num_final_OMEGA,
      integer num_cyclic_OMEGA,
      // ----------------------
      integer num_extra_r,
      integer num_extra_i
    );

    /*!
     * \brief Convenience variant for systems with only top and bottom blocks.
     *
     * It converts the dimensions of `top` and `bottom` blocks into the general
     * class convention and forwards to \ref allocate.
     *
     * \param nblock number of internal blocks.
     * \param n size of the square internal blocks.
     * \param row0 number of rows in the top block.
     * \param col0 number of columns in the top block.
     * \param rowN number of rows in the bottom block.
     * \param colN number of columns in the bottom block.
     * \param nb size of the right border.
     * \param num_extra_r extra real storage requested by the derived class.
     * \param num_extra_i extra integer storage requested by the derived class.
     */
    void
    allocate_top_bottom(
      integer nblock,
      integer n,
      integer row0,
      integer col0,
      integer rowN,
      integer colN,
      integer nb,
      integer num_extra_r,
      integer num_extra_i
    ) {
      allocate(
        nblock, n, nb,
        row0, rowN, 0,
        col0-n, colN-n, 0,
        num_extra_r, num_extra_i
      );
    }

    //! \name Loading the bidiagonal part
    //! @{

    //! \brief Loads all internal `[Ad | Au]` blocks.
    void load_blocks( real_type const AdAu[], integer ldA );
    //! \brief Loads full internal block `nbl` as `[Ad | Au]`.
    void load_block( integer nbl, real_type const AdAu[], integer ldA );

    //! \brief Loads only the left diagonal part `Ad` of block `nbl`.
    void load_block_left( integer nbl, real_type const Ad[], integer ldA );
    //! \brief Loads only the right superdiagonal part `Au` of block `nbl`.
    void load_block_right( integer nbl, real_type const Au[], integer ldA );

    //! @}
    //! \name Lower-border blocks
    //! @{

    //! \brief Zeros all lower-border blocks.
    void set_zero_bottom_blocks();
    //! \brief Loads all lower-border blocks.
    void load_bottom_blocks( real_type const C[], integer ldC );
    //! \brief Loads the lower block associated with internal block `nbl`.
    void load_bottom_block( integer nbl, real_type const C[], integer ldC );
    //! \brief Adds a contribution to lower block `nbl`.
    void add_to_bottom_block( integer nbl, real_type const C[], integer ldC );

    //! \brief Adds a contribution distributed across blocks `nbl` and `nbl+1`.
    void add_to_bottom_block2( integer nbl, real_type const C[], integer ldC );
    //! \brief Loads the last lower-border block.
    void load_bottom_last_block( real_type const C[], integer ldC );

    //! @}
    //! \name Right-border blocks
    //! @{

    //! \brief Zeros all right-border blocks.
    void set_zero_right_blocks();
    //! \brief Loads all right-border blocks.
    void load_right_blocks( real_type const B[], integer ldB );
    //! \brief Loads right-border block `nbl`.
    void loadRightBlock( integer nbl, real_type const B[], integer ldB );
    //! \brief Loads the last right-border block.
    void loadRightLastBlock( real_type const B[], integer ldB );

    //! @}
    //! \name Bottom-right corner block
    //! @{

    //! \brief Zeros the bottom-right corner block.
    void setZeroRBblock();
    //! \brief Loads the bottom-right corner block.
    void load_RB_block( real_type const D[], integer ldD );

    //! @}
    //! \name Extraction of the final reduced block
    //! @{

    //! \brief Returns the pointer to the concatenated final reduced block.
    real_type const * getPointer_LR() const { return m_DE_blk; }

    //! \brief Copies the concatenated final block `[L | R]`.
    void getBlock_LR( real_type LR[], integer ldA ) const;
    //! \brief Copies the left final block `L`.
    void getBlock_L ( real_type L[],  integer ldA ) const;
    //! \brief Copies the right final block `R`.
    void getBlock_R ( real_type R[],  integer ldA ) const;
    //! \brief Copies reduced-lower-border block `H0`.
    void getBlock_H0( real_type H0[], integer ld0 ) const;
    //! \brief Copies reduced-lower-border block `HN`.
    void getBlock_HN( real_type HN[], integer ldN ) const;
    //! \brief Copies reduced-lower-border block `Hq`.
    void getBlock_Hq( real_type Hq[], integer ldQ ) const;

    //! @}

    /*!
     * \brief Specialized allocation implemented by the derived class.
     *
     * This overload receives only the problem geometry; base storage is already
     * managed by the current class.
     */
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
    ) = 0;

    /*!
     * \brief Specialized allocation for systems described by top/bottom blocks.
     *
     * Derived classes may use this signature to implement more natural
     * allocation paths for their specific algorithms.
     */
    virtual
    void
    allocate_top_bottom(
      integer /* nblock */,
      integer /* n      */,
      integer /* row0   */,
      integer /* col0   */,
      integer /* rowN   */,
      integer /* colN   */,
      integer /* nb     */
    ) = 0;

    //! \brief Factorizes the system loaded in the internal structure.
    virtual
    void
    factorize() = 0;

    //! \brief Solves the system for a single right-hand side.
    virtual
    void
    solve( real_type [] ) const = 0;

    //! \brief Solves the system for multiple right-hand sides.
    virtual
    void
    solve(
      integer      /* nrhs  */,
      real_type [] /* rhs   */,
      integer      /* ldRhs */
    ) const = 0;

    //! \brief Number of internal blocks.
    integer number_of_blocks() const { return m_number_of_blocks; }
    //! \brief Size of the square internal block.
    integer block_size()       const { return m_block_size; }
    //! \brief Total number of equations in the internal bidiagonal block structure.
    integer num_equations()    const { return m_num_equations; }
    //! \brief Size of the right / bottom-right border.
    integer border_size()      const { return m_border_size; }

    /*!
     * \brief Loads the reduced lower border through blocks `H0`, `HN`, `Hq`.
     * \param H0 contribution connected to the first block.
     * \param ld0 leading dimension di `H0`.
     * \param HN contribution connected to the last block.
     * \param ldN leading dimension di `HN`.
     * \param Hq contribution associated with the extra variables.
     * \param ldQ leading dimension di `Hq`.
     */
    void
    load_bottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ
    );

    /*!
     * \brief Loads `top` and `bottom` blocks in the top/bottom parameterization.
     * \param block0 top block.
     * \param ld0 leading dimension di `block0`.
     * \param blockN bottom block.
     * \param ldN leading dimension di `blockN`.
     *
     * Block `top` has size `row0 x col0`, while `bottom` has size `rowN x colN`
     * according to the geometry fixed during allocation.
     */
    void
    load_top_bottom(
      real_type const block0[], integer ld0,
      real_type const blockN[], integer ldN
    );

    /*!
     * \brief Selects the solver to use on the final reduced block.
     * \param choice desired numerical strategy.
     */
    void
    select_last_block_solver( BB_LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case BB_LASTBLOCK_Choice::LU:   m_la_factorization = &m_la_lu;   break;
        case BB_LASTBLOCK_Choice::LUPQ: m_la_factorization = &m_la_lupq; break;
        case BB_LASTBLOCK_Choice::QR:   m_la_factorization = &m_la_qr;   break;
        case BB_LASTBLOCK_Choice::QRP:  m_la_factorization = &m_la_qrp;  break;
        case BB_LASTBLOCK_Choice::SVD:  m_la_factorization = &m_la_svd;  break;
        case BB_LASTBLOCK_Choice::LSS:  m_la_factorization = &m_la_lss;  break;
        case BB_LASTBLOCK_Choice::LSY:  m_la_factorization = &m_la_lsy;  break;
        case BB_LASTBLOCK_Choice::PINV: m_la_factorization = &m_la_pinv; break;
      }
    }

    //! \brief Selects LU for the final reduced block.
    void select_last_block_solver_LU()   { m_la_factorization = &m_la_lu;   }
    //! \brief Selects LUPQ for the final reduced block.
    void select_last_block_solver_LUPQ() { m_la_factorization = &m_la_lupq; }
    //! \brief Selects QR for the final reduced block.
    void select_last_block_solver_QR()   { m_la_factorization = &m_la_qr;   }
    //! \brief Selects QRP for the final reduced block.
    void select_last_block_solver_QRP()  { m_la_factorization = &m_la_qrp;  }
    //! \brief Selects SVD for the final reduced block.
    void select_last_block_solver_SVD()  { m_la_factorization = &m_la_svd;  }
    //! \brief Selects LSS for the final reduced block.
    void select_last_block_solver_LSS()  { m_la_factorization = &m_la_lss;  }
    //! \brief Selects LSY for the final reduced block.
    void select_last_block_solver_LSY()  { m_la_factorization = &m_la_lsy;  }
    //! \brief Selects PINV for the final reduced block.
    void select_last_block_solver_PINV() { m_la_factorization = &m_la_pinv; }

    /*!
     * \brief Selects the solver to use on the final bordered system.
     * \param choice desired numerical strategy.
     */
    void
    select_last_border_block_solver( BB_LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case BB_LASTBLOCK_Choice::LU:   m_bb_factorization = &m_bb_lu;   break;
        case BB_LASTBLOCK_Choice::LUPQ: m_bb_factorization = &m_bb_lupq; break;
        case BB_LASTBLOCK_Choice::QR:   m_bb_factorization = &m_bb_qr;   break;
        case BB_LASTBLOCK_Choice::QRP:  m_bb_factorization = &m_bb_qrp;  break;
        case BB_LASTBLOCK_Choice::SVD:  m_bb_factorization = &m_bb_svd;  break;
        case BB_LASTBLOCK_Choice::LSS:  m_bb_factorization = &m_bb_lss;  break;
        case BB_LASTBLOCK_Choice::LSY:  m_bb_factorization = &m_bb_lsy;  break;
        case BB_LASTBLOCK_Choice::PINV: m_bb_factorization = &m_bb_pinv; break;
      }
    }

    //! \brief Selects LU for the final bordered system.
    void select_last_border_block_solver_LU()   { m_bb_factorization = &m_bb_lu;   }
    //! \brief Selects LUPQ for the final bordered system.
    void select_last_border_block_solver_LUPQ() { m_bb_factorization = &m_bb_lupq; }
    //! \brief Selects QR for the final bordered system.
    void select_last_border_block_solver_QR()   { m_bb_factorization = &m_bb_qr;   }
    //! \brief Selects QRP for the final bordered system.
    void select_last_border_block_solver_QRP()  { m_bb_factorization = &m_bb_qrp;  }
    //! \brief Selects SVD for the final bordered system.
    void select_last_border_block_solver_SVD()  { m_bb_factorization = &m_bb_svd;  }
    //! \brief Selects LSS for the final bordered system.
    void select_last_border_block_solver_LSS()  { m_bb_factorization = &m_bb_lss;  }
    //! \brief Selects LSY for the final bordered system.
    void select_last_border_block_solver_LSY()  { m_bb_factorization = &m_bb_lsy;  }
    //! \brief Selects PINV for the final bordered system.
    void select_last_border_block_solver_PINV() { m_bb_factorization = &m_bb_pinv; }

    //! \brief Factorizes the final reduced block without border terms.
    void
    last_block_factorize();

    //! \brief Factorizes the reduced system with right/lower border terms.
    void
    factorize_bordered();

    //! \brief Solves the bordered system for a single right-hand side.
    void
    solve_bordered( real_type [] ) const;

    //! \brief Solves the bordered system for multiple right-hand sides.
    void
    solve_bordered(
      integer      /* nrhs  */,
      real_type [] /* rhs   */,
      integer      /* ldRhs */
    ) const;

    /*!
     * \brief Loads and factorizes a bordered system in a single step.
     *
     * Blocks are provided in dense column-major form, consistent with the
     * structure of the already allocated problem.
     */
    void
    factorize(
      real_type const AdAu[],
      real_type const B[],
      real_type const C[],
      real_type const D[],
      real_type const H0[],
      real_type const HN[],
      real_type const Hq[]
    ) {
      integer const & n   = m_block_size;
      integer const & nb  = m_border_size;
      integer const & q   = m_extra_bc;
      integer const & neq = m_num_equations;

      this->load_blocks( AdAu, n );
      integer nq = n + q;
      this->load_bottom( H0, nq, HN, nq, Hq, nq );
      if ( nb > 0 ) {
        this->load_right_blocks( B, neq );
        this->load_bottom_blocks( C, nb );
        this->load_RB_block( D, nb );
        this->factorize_bordered();
      } else {
        this->factorize();
      }
    }

    /*!
     * \brief Loads and factorizes a system without a right border in a single step.
     */
    void
    factorize(
      real_type const AdAu[],
      real_type const H0[],
      real_type const HN[],
      real_type const Hq[]
    ) {
      integer const & n  = m_block_size;
      integer const & q  = m_extra_bc;
      integer const & nb = m_border_size;

      UTILS_ASSERT0( nb == 0, "factorize nb > 0 and no border assigned\n" );
      integer nq = n + q;
      this->load_blocks( AdAu, n );
      this->load_bottom( H0, nq, HN, nq, Hq, nq );
      this->factorize();
    }

    /*!
     * \brief Computes the matrix-vector product `res = A*x`.
     * \param x input vector.
     * \param res output buffer.
     */
    void
    Mv( real_type const x[], real_type res[] ) const;

    /*\
     |   ____
     |  |  _ \ _   _ _ __ ___  _ __
     |  | | | | | | | '_ ` _ \| '_ \
     |  | |_| | |_| | | | | | | |_) |
     |  |____/ \__,_|_| |_| |_| .__/
     |                        |_|
    \*/

    //! \brief Exports the matrix in readable column-coordinate form.
    void
    dump_ccoord( ostream_type & stream ) const;

    //! \brief Prints the matrix in a Maple-compatible format.
    void
    dump_to_Maple( ostream_type & stream ) const;

    /*\
     |   ___ _ __   __ _ _ __ ___  ___
     |  / __| '_ \ / _` | '__/ __|/ _ \
     |  \__ \ |_) | (_| | |  \__ \  __/
     |  |___/ .__/ \__,_|_|  |___/\___|
     |      |_|
    \*/

    //! \brief Number of nonzero entries in the full matrix.
    integer
    sparse_nnz() const;

    //! \brief Exports the sparse pattern of the full matrix.
    void
    sparse_pattern( integer I[], integer J[] ) const;

    //! \brief Exports the values of the full matrix in sparse order.
    void
    sparse_values( real_type vals[] ) const;

  };
}

///
/// eof: BlockBidiagonal.hxx
///
