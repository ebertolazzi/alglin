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
 |      Universita' degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BABD.hxx
///

/*!
 * \file BABD.hxx
 * \brief High-level wrapper for the BABD solvers available in the library.
 *
 * This header exposes class `alglin::BABD`, a lightweight dispatcher that
 * explicitly selects one of the real backends for block-almost-bidiagonal
 * systems:
 * - `DiazLU`, suitable for ABD systems with top/bottom blocks;
 * - `BorderedCR`, suitable for systems solved by cyclic reduction.
 *
 * The wrapper keeps only the operations that are genuinely common to the
 * backends (`factorize_system`, `solve_system`, sparse export, and dimension
 * queries) and leaves system loading to the active backend through `diaz()` or
 * `cr()`.
 */

#include <variant>

namespace alglin {

  /*!
   * \brief Dedicated wrapper around the real backends for BABD systems.
   *
   * The class does not try to hide the interface differences between `DiazLU`
   * and `BorderedCR`: it only exposes common factorization/solve operations
   * and, for backend-specific system setup, provides explicit access to the
   * active backend.
   *
   * The backend is selected through \ref BABD_Choice. When the active backend
   * is `BorderedCR`, the choice of the final elimination algorithm is
   * automatically propagated to the underlying solver.
   */
  template <typename t_Value>
  class BABD {
  public:

    using real_type = t_Value;
    using DiazBackend = DiazLU<t_Value>;
    using CRBackend   = BorderedCR<t_Value>;
    using block_choice = typename BlockBidiagonal<t_Value>::BB_LASTBLOCK_Choice;

    //! \brief High-level solution strategy selectable in the wrapper.
    using BABD_Choice = enum class BABD_Choice : integer {
      DIAZ                 = 1, //!< Backend `DiazLU`.
      CYCLIC_REDUCTION_LU  = 2, //!< `BorderedCR` backend with cyclic reduction and LU last block.
      CYCLIC_REDUCTION_QR  = 3, //!< `BorderedCR` backend with cyclic reduction and QR last block.
      CYCLIC_REDUCTION_QRP = 4  //!< `BorderedCR` backend with cyclic reduction and QRP last block.
    };

    static
    string
    //! \brief Converts a wrapper choice into a readable string.
    //! \param c choice to convert.
    //! \return Symbolic name also used in diagnostic strings.
    Choice_to_string( BABD_Choice c ) {
      switch ( c ) {
        case BABD_Choice::DIAZ:                 return "Diaz";
        case BABD_Choice::CYCLIC_REDUCTION_LU:  return "CyclicReduction+LU";
        case BABD_Choice::CYCLIC_REDUCTION_QR:  return "CyclicReduction+QR";
        case BABD_Choice::CYCLIC_REDUCTION_QRP: return "CyclicReduction+QRP";
      }
      return "none";
    }

  private:

    using backend_type = std::variant<DiazBackend, CRBackend>;

    BABD( BABD<t_Value> const & ) = delete;
    BABD<t_Value> const & operator = ( BABD<t_Value> const & ) = delete;

    BABD_Choice             m_choice{ BABD_Choice::CYCLIC_REDUCTION_LU };
    Utils::ThreadPoolBase * m_thread_pool{ nullptr };
    integer                 m_num_parallel_block{ 0 };
    backend_type            m_backend;

    static
    constexpr
    bool
    is_cr_choice( BABD_Choice choice ) {
      return choice != BABD_Choice::DIAZ;
    }

    void
    configure_cr_choice() {
      switch ( m_choice ) {
        case BABD_Choice::CYCLIC_REDUCTION_LU:
          this->cr().select_LU();
          break;
        case BABD_Choice::CYCLIC_REDUCTION_QR:
          this->cr().select_QR();
          break;
        case BABD_Choice::CYCLIC_REDUCTION_QRP:
          this->cr().select_QRP();
          break;
        case BABD_Choice::DIAZ:
          break;
      }
    }

  public:

    /*!
     * \brief Builds the wrapper by choosing the initial backend.
     * \param choice initial backend/algorithm choice.
     * \param TP optional thread pool to forward to `BorderedCR`.
     * \param num_parallel_block maximum number of blocks processed in parallel
     *   by the `BorderedCR` backend.
     */
    explicit
    BABD(
      BABD_Choice             choice = BABD_Choice::CYCLIC_REDUCTION_LU,
      Utils::ThreadPoolBase * TP     = nullptr,
      integer                 num_parallel_block = 0
    )
    : m_thread_pool(TP)
    , m_num_parallel_block(num_parallel_block)
    , m_backend(std::in_place_type<CRBackend>, TP, num_parallel_block) {
      this->selectSolver( choice );
    }

    ~BABD() = default;

    //! \brief Returns the currently selected strategy.
    BABD_Choice selected_solver() const { return m_choice; }

    //! \brief Returns `true` when the active backend is `DiazLU`.
    bool using_diaz()              const { return std::holds_alternative<DiazBackend>(m_backend); }
    //! \brief Returns `true` when the active backend is `BorderedCR`.
    bool using_cyclic_reduction()  const { return std::holds_alternative<CRBackend>(m_backend); }

    /*!
     * \brief Access to backend `DiazLU`.
     * \return Reference to the active backend.
     * \note Triggers an assertion if the active backend is not `DiazLU`.
     */
    DiazBackend &
    diaz() {
      UTILS_ASSERT0( this->using_diaz(), "BABD::diaz(), active backend is not DiazLU\n" );
      return std::get<DiazBackend>( m_backend );
    }

    /*!
     * \brief Const access to backend `DiazLU`.
     * \return Const reference to the active backend.
     * \note Triggers an assertion if the active backend is not `DiazLU`.
     */
    DiazBackend const &
    diaz() const {
      UTILS_ASSERT0( this->using_diaz(), "BABD::diaz() const, active backend is not DiazLU\n" );
      return std::get<DiazBackend>( m_backend );
    }

    /*!
     * \brief Access to backend `BorderedCR`.
     * \return Reference to the active backend.
     * \note Triggers an assertion if the active backend is not `BorderedCR`.
     */
    CRBackend &
    cr() {
      UTILS_ASSERT0(
        this->using_cyclic_reduction(),
        "BABD::cr(), active backend is not BorderedCR\n"
      );
      return std::get<CRBackend>( m_backend );
    }

    /*!
     * \brief Const access to backend `BorderedCR`.
     * \return Const reference to the active backend.
     * \note Triggers an assertion if the active backend is not `BorderedCR`.
     */
    CRBackend const &
    cr() const {
      UTILS_ASSERT0(
        this->using_cyclic_reduction(),
        "BABD::cr() const, active backend is not BorderedCR\n"
      );
      return std::get<CRBackend>( m_backend );
    }

    /*!
     * \brief Changes the backend or the final factorization algorithm.
     *
     * If the new choice requires backend `BorderedCR`, the wrapper creates the
     * backend if needed and realigns the LU/QR/QRP selection. If strategy
     * `DIAZ` is requested instead, the active backend is replaced with a
     * `DiazLU` object.
     *
     * \param choice new solver choice.
     */
    void
    selectSolver( BABD_Choice choice ) {
      bool const need_cr{ is_cr_choice(choice) };
      bool const has_cr { this->using_cyclic_reduction() };

      m_choice = choice;
      if ( need_cr ) {
        if ( !has_cr )
          m_backend.template emplace<CRBackend>( m_thread_pool, m_num_parallel_block );
        this->configure_cr_choice();
      } else if ( has_cr ) {
        m_backend.template emplace<DiazBackend>();
      }
    }

    /*!
     * \brief Compact description of the active algorithm.
     * \return Backend name or descriptive string for the final solver.
     */
    string
    info() const {
      if ( this->using_diaz() ) return Choice_to_string( m_choice );
      return this->cr().info_algo();
    }

    /*!
     * \brief Number of rows of the global linear system.
     * \return Size of the residual vector associated with the active backend.
     */
    integer
    nrows() const {
      if ( this->using_diaz() ) {
        DiazBackend const & s{ this->diaz() };
        return s.num_equations() + s.border_size();
      }
      return this->cr().nrows();
    }

    /*!
     * \brief Number of columns of the global linear system.
     * \return Size of the global unknown associated with the active backend.
     */
    integer
    ncols() const {
      if ( this->using_diaz() ) {
        DiazBackend const & s{ this->diaz() };
        return s.num_equations() + s.border_size();
      }
      return this->cr().ncols();
    }

    /*!
     * \brief Factorizes the system already loaded in the active backend.
     * \return `true` if the active backend exposes a boolean status, or `true`
     *   after execution of backend `DiazLU`.
     */
    bool
    factorize_system() {
      if ( this->using_diaz() ) {
        this->diaz().factorize_bordered();
        return true;
      }
      return this->cr().factorize();
    }

    /*!
     * \brief Solves the system for a single right-hand side.
     * \param in_out Contains the right-hand side on input and the solution on output.
     * \return `true` if the solve completes successfully.
     */
    bool
    solve_system( real_type in_out[] ) const {
      if ( this->using_diaz() ) {
        this->diaz().solve_bordered( in_out );
        return true;
      }
      return this->cr().solve( in_out );
    }

    /*!
     * \brief Solves the system for multiple right-hand sides.
     * \param nrhs number of columns in the right-hand-side block.
     * \param in_out column-major matrix containing right-hand sides on input and
     *   the corresponding solutions on output.
     * \param ldIO leading dimension of `in_out`.
     * \return `true` if the solve completes successfully.
     */
    bool
    solve_system(
      integer   nrhs,
      real_type in_out[],
      integer   ldIO
    ) const {
      if ( this->using_diaz() ) {
        this->diaz().solve_bordered( nrhs, in_out, ldIO );
        return true;
      }
      return this->cr().solve( nrhs, in_out, ldIO );
    }

    //! \brief Short alias for \ref factorize_system.
    bool
    factorize() {
      return this->factorize_system();
    }

    //! \brief Short alias for \ref solve_system(real_type[]) const.
    bool
    solve( real_type in_out[] ) const {
      return this->solve_system( in_out );
    }

    //! \brief Short alias for \ref solve_system(integer,real_type[],integer) const.
    bool
    solve(
      integer   nrhs,
      real_type in_out[],
      integer   ldIO
    ) const {
      return this->solve_system( nrhs, in_out, ldIO );
    }

    /*!
     * \brief Number of nonzero entries in the sparse representation of the system.
     * \return Structural nonzero count of the active backend.
     */
    integer
    sparse_nnz() const {
      if ( this->using_diaz() ) return this->diaz().sparse_nnz();
      return this->cr().sparse_nnz();
    }

    /*!
     * \brief Exports the sparse pattern of the global system.
     * \param I row-index array for the nonzeros.
     * \param J column-index array for the nonzeros.
     */
    void
    sparse_pattern( integer I[], integer J[] ) const {
      if ( this->using_diaz() ) this->diaz().sparse_pattern( I, J );
      else                      this->cr().sparse_pattern( I, J, 0 );
    }

    /*!
     * \brief Exports the values of the global system in the same order as \ref sparse_pattern.
     * \param vals output buffer for the nonzero values.
     */
    void
    sparse_values( real_type vals[] ) const {
      if ( this->using_diaz() ) this->diaz().sparse_values( vals );
      else                      this->cr().sparse_values( vals );
    }
  };
}

///
/// eof: BABD.hxx
///
