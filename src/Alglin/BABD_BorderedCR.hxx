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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BABD_BorderedCR.hxx
///

/*!
 * \file BABD_BorderedCR.hxx
 * \brief Cyclic-reduction-based solver for bordered BABD systems.
 *
 * This header defines `alglin::BorderedCR`, the most feature-complete class in
 * the module for block almost bidiagonal systems with border terms on both the
 * right and lower sides. The solver supports several numerical strategies:
 * - cyclic reduction with LU, QR, or QRP block eliminations;
 * - multiple solvers for the final reduced block;
 * - an optional fallback solver for ill-conditioned or singular cases.
 *
 * The class provides fine-grained APIs to load individual blocks, export the
 * system in sparse form, perform matrix-vector products, and control internal
 * parallelism.
 */

namespace alglin {

  using std::atomic;
  using std::swap;
  using std::min;
  using std::max;
  using std::string_view;

  /*\
   |   ___             _                _    ___ ___
   |  | _ ) ___ _ _ __| |___ _ _ ___ __| |  / __| _ \
   |  | _ \/ _ \ '_/ _` / -_) '_/ -_) _` | | (__|   /
   |  |___/\___/_| \__,_\___|_| \___\__,_|  \___|_|_\
  \*/

  /*!
   * \brief Bordered BABD solver based on cyclic reduction.
   *
   * The handled system has an internal block-bidiagonal structure augmented by
   * right and lower border blocks. The class separates loading of the various
   * matrix components, selection of the final reduced-block solvers, and the
   * factorization/solve phase.
   *
   * The solver supports both single-thread use and execution on an external
   * thread pool, with explicit control over the maximum number of reduced
   * blocks processed in parallel.
   *
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   * \author   Enrico Bertolazzi
   *
   * \par Matrix structure
   * \verbatim
   *                 n * (nblock+1)
   *  ┌───────────────────┴─────────────────────┐
   *   n   n   n                              n  qx  nx
   * ╔═══╦═══╦═══╦═════................═════╦═══╦═══╦═══╗   ─┐
   * ║ D ║ E ║   ║                          ║   ║   ║ B ║ n  │
   * ╠═══╬═══╬═══╬═══                  ═════╬═══╬═══╬═══╣    │
   * ║   ║ D ║ E ║                          ║   ║   ║ B ║ n  │
   * ╠═══╬═══╬═══╬═══╬                 ═════╬═══╬═══╬═══╣    │
   * ║   ║   ║ D ║ E ║                      ║   ║   ║ B ║ n  │
   * ╠═══╬═══╬═══╬═══╬                 ═════╬═══╬═══╬═══╣    │
   * :                                                  :    │
   * :                                                  :    │
   * :                                                  :    ├─ n * nblock
   * :                                                  :    │
   * :                                                  :    │
   * :                              ╬═══╬═══╬═══╬═══╬═══╣    │
   * :                              ║ D ║ E ║   ║   ║ B ║ n  │
   * :                              ╬═══╬═══╬═══╬═══╬═══╣    │
   * :                                  ║ D ║ E ║   ║ B ║ n  │
   * ╠═══╬═══╬═══................═══╬═══╬═══╬═══╬═══╬═══╣   ─┤
   * ║   ║   ║                          ║   ║   ║   ║   ║    │
   * ║H0 ║ 0 ║                          ║ 0 ║HN ║ Hq║ Hp║    │ n+qr
   * ║   ║   ║                          ║   ║   ║   ║   ║    │
   * ╠═══╬═══╬═══................═══╬═══╬═══╬═══╬═══╬═══╣   ─┤
   * ║ C ║ C ║                      ║ C ║ C ║ C ║ Cq║ F ║    │ nr
   * ╚═══╩═══╩═══................═══╩═══╩═══╩═══╩═══╩═══╝   ─┘
   *                                            nr*qx
   *
   * \endverbatim
   */
  template <typename t_Value>
  class BorderedCR : public LinearSystemSolver<t_Value> {
  public:
    using real_type = t_Value;
    using MatW      = MatrixWrapper<t_Value>;

    //! \brief Solvers available for blocks eliminated during cyclic reduction.
    using BORDERED_Choice = enum class BORDERED_Choice : integer {
      LU  = 0, //!< Reduction using LU factorization.
      QR  = 1, //!< Reduction using QR factorization.
      QRP = 2  //!< Reduction using QR factorization with pivoting.
    };

    //! \brief Primary solver applied to the final reduced system.
    using BORDERED_LAST_Choice = enum class BORDERED_LAST_Choice : integer {
      LU   = 0, //!< Standard LU.
      LUPQ = 1, //!< LU with complete pivoting.
      QR   = 2, //!< Standard QR.
      QRP  = 3, //!< QR with pivoting.
      SVD  = 4, //!< Singular value decomposition.
      LSS  = 5, //!< Least-squares standard.
      LSY  = 6, //!< Least-squares with alternate solver.
      PINV = 7  //!< Pseudoinverse.
    };

    //! \brief Fallback solver for the final block, used only when needed.
    using BORDERED_LAST_Choice2 = enum class BORDERED_LAST_Choice2 : integer {
      NONE = 0, //!< No fallback.
      SVD  = 1, //!< Fallback SVD.
      LSS  = 2, //!< Fallback least-squares standard.
      LSY  = 3, //!< Fallback alternate least-squares.
      PINV = 4  //!< Fallback using the pseudoinverse.
    };

    // block copy constructor
    BorderedCR() = delete;
    BorderedCR( BorderedCR const & ) = delete;
    BorderedCR const & operator = ( BorderedCR const & ) = delete;

    //! \brief Enables internal debug prints and checks.
    bool m_debug{false};

  protected:

    Malloc<real_type> m_mem{"BorderedCR_values"};
    Malloc<integer>   m_mem_int{"BorderedCR_integers"};

    mutable Malloc<real_type> m_work_mem{"BorderedCR_work_mem"};
    mutable string            m_last_error{"no error"};

    integer m_number_of_blocks{0}; //!< total number of blocks
    integer m_block_size{0};       //!< size of square blocks
    integer m_qr{0};               //!< extra BC
    integer m_qx{0};               //!< extra BC
    integer m_nr{0};               //!< border size
    integer m_nx{0};               //!< border size

    integer m_Nr{0};
    integer m_Nc{0};
    integer m_Tsize{0};

    // some derived constants
    integer n_x_2{0};
    integer n_x_n{0};
    integer n_x_nx{0};
    integer nr_x_n{0};
    integer nr_x_nx{0};
    integer nr_x_qx{0};

    bool m_factorize_use_thread{true};
    bool m_solve_use_thread{true};
    bool m_matrix_is_factorized{false};

    BORDERED_Choice       m_selected{BORDERED_Choice::LU};
    BORDERED_LAST_Choice  m_last_block_selected{BORDERED_LAST_Choice::LU};
    BORDERED_LAST_Choice2 m_last_block_selected2{BORDERED_LAST_Choice2::NONE};
    bool                  m_last_block_use_solver2{false};

    real_type * m_H0Nqp{nullptr};
    real_type * m_Bmat{nullptr};
    real_type * m_Cmat{nullptr};
    real_type * m_Cqmat{nullptr};
    real_type * m_Dmat{nullptr};
    real_type * m_Emat{nullptr};
    real_type * m_Fmat{nullptr};

    vector<real_type*> m_xb_thread;

    integer            m_Work_Lapack_size{0};
    vector<real_type*> m_Work_Lapack_thread;

    mutable vector<vector<real_type>> m_Work_T_thread;

    // working block
    real_type * m_Tmat{nullptr};
    real_type * m_Ttau{nullptr};
    integer *   m_Perm{nullptr};


    // last block
    real_type * m_Hmat{nullptr};

    LU<real_type>   m_last_lu;
    LUPQ<real_type> m_last_lupq;
    QR<real_type>   m_last_qr;
    QRP<real_type>  m_last_qrp;
    SVD<real_type>  m_last_svd;
    LSS<real_type>  m_last_lss;
    LSY<real_type>  m_last_lsy;
    PINV<real_type> m_last_pinv;

    integer * m_iBlock{nullptr};
    integer * m_kBlock{nullptr};

    // used also with a unique thread
    integer m_max_parallel_block{0};
    integer m_used_parallel_block{0};
    integer m_reduced_nblk{0};
    mutable vector<integer*> m_perm_thread;

    mutable UTILS_SPINLOCK m_spin;
    mutable atomic<bool>   m_ok_thread;

    Utils::ThreadPool0      m_TP_fake{1};
    Utils::ThreadPoolBase * m_TP{&m_TP_fake};

    bool
    buildT(
      integer   nth,
      real_type T[],
      integer   iperm[]
    ) const;

    bool
    applyT(
      integer         nth,
      real_type const T[],
      integer   const iperm[],
      real_type       TOP[],
      integer         ldTOP,
      real_type       BOTTOM[],
      integer         ldBOTTOM,
      integer         ncol
    ) const;

    bool
    applyT(
      integer         nth,
      real_type const T[],
      integer   const iperm[],
      real_type       TOP[],
      real_type       BOTTOM[]
    ) const;

    real_type *
    Work_T( integer n_thread, integer size ) const {
      vector<real_type> & W{ m_Work_T_thread[n_thread] };
      if ( integer(W.size()) < size ) W.resize(size);
      return W.data();
    }

    // convert permutation to exchanges
    void
    permutation_to_exchange( integer nn, integer P[], integer S[] ) const {
      for ( integer i{0}; i < nn; ++i ) {
        integer j = i;
        while ( j < nn ) { if ( P[j] == i+1 ) break; ++j; }
        //UTILS_ASSERT0( j < nn, "permutation_to_exchange error!\n" );
        std::swap( P[j], P[i] );
        S[i] = j;
      }
    }

    /*
    //    __         _           _
    //   / _|__ _ __| |_ ___ _ _(_)______
    //  |  _/ _` / _|  _/ _ \ '_| |_ / -_)
    //  |_| \__,_\__|\__\___/_| |_/__\___|
    */

    void factorize_block( integer nth );
    bool factorize_reduced();

    /*
    //    __                            _
    //   / _|___ _ ___ __ ____ _ _ _ __| |
    //  |  _/ _ \ '_\ V  V / _` | '_/ _` |
    //  |_| \___/_|  \_/\_/\__,_|_| \__,_|
    */

    bool forward( integer nth, real_type x[] ) const;

    void
    forward_n(
      integer   nth,
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs,
      real_type work[]
    ) const;

    bool
    forward_reduced( real_type x[], real_type xb[] ) const;

    bool
    forward_n_reduced(
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs
    ) const;

    /*
    //   _             _                       _
    //  | |__  __ _ __| |____ __ ____ _ _ _ __| |
    //  | '_ \/ _` / _| / /\ V  V / _` | '_/ _` |
    //  |_.__/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    */

    void backward( integer nth, real_type x[] ) const;
    void backward_reduced( real_type x[] ) const;

    void
    backward_n(
      integer   nth,
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs
    ) const;

    void
    backward_n_reduced(
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs
    ) const;

    bool
    load_and_factorize_last();

    /*
    //   _         _
    //  | |__ _ __| |_
    //  | / _` (_-<  _|
    //  |_\__,_/__/\__|
    */

    bool
    solve_last( real_type [] ) const;

    bool
    solve_last(
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs
    ) const;

  public:

    using LinearSystemSolver<t_Value>::factorize;
    using LinearSystemSolver<t_Value>::solve;
    using LinearSystemSolver<t_Value>::t_solve;

    /*!
     * \brief Builds the solver.
     * \param TP optional external thread pool.
     * \param num_parallel_block maximum number of reduced blocks processed in parallel.
     *
     * If `TP` is null, a local single-thread pool is created.
     */
    explicit
    BorderedCR( Utils::ThreadPoolBase * TP, integer num_parallel_block = 0 ) {
      if ( TP != nullptr ) {
        m_TP = TP;
        m_max_parallel_block = num_parallel_block == 0 ?
                               TP->thread_count() :
                               min( integer(TP->thread_count()), num_parallel_block );
      } else {
        m_TP                 = &m_TP_fake;
        m_max_parallel_block = 1;
      }
    }

    virtual
    ~BorderedCR() override
    {}

    //! \brief Enables or disables parallelism in factorization.
    void set_factorize_use_thread( bool yes_no ) { m_factorize_use_thread = yes_no; }
    //! \brief Enables or disables parallelism in the solve phase.
    void set_solve_use_thread( bool yes_no )     { m_solve_use_thread     = yes_no; }

    //! \brief Returns `true` if factorization uses threads.
    bool factorize_use_thread() const { return m_factorize_use_thread; }
    //! \brief Returns `true` if the solve phase uses threads.
    bool solve_use_thread()     const { return m_solve_use_thread; }

    /*!
     * \brief Sets the maximum number of reduced blocks processed in parallel.
     * \param num_parallel_block new upper bound.
     *
     * If the problem has already been allocated, the internal geometry is
     * reallocated so that the new choice becomes effective.
     */
    void
    set_num_parallel_block( integer num_parallel_block ) {
      // Clamp the maximum number of parallel blocks to the number of threads.
      m_max_parallel_block = max(
        integer(1),
        min( integer(m_TP->thread_count()), num_parallel_block )
      );
      if ( m_number_of_blocks > 0 )
        this->allocate(
          m_number_of_blocks,
          m_block_size,
          m_qr,
          m_qx,
          m_nr,
          m_nx
        );
    }

    /*!
     * \brief Allocates the system and initializes internal buffers.
     * \param nblock number of internal blocks.
     * \param n size of the square internal blocks.
     * \param qr number of extra rows in the upper-lower border block.
     * \param qx number of extra columns in the lower border.
     * \param nr number of rows in the lower `C/F` border.
     * \param nx number of columns in the right border.
     */
    void
    allocate(
      integer nblock,
      integer n,
      integer qr,
      integer qx,
      integer nr,
      integer nx
    );

    //! \brief Last internally recorded error message.
    string_view last_error() const { return m_last_error; }

    /*!
     * \brief Fully duplicates the state of another solver.
     * \param other source object to copy from.
     */
    void dup( BorderedCR const & );

    //!
    //! \name Select Linar Algebra solver
    //! @{
    //!

    //! \brief Selects the reduction-phase solver from a string.
    void
    select( string_view s ) {
      if      ( s == "LU"  ) select_LU();
      else if ( s == "QR"  ) select_QR();
      else if ( s == "QRP" ) select_QRP();
      else {
        UTILS_ERROR(
          "BorderedCR::select( choice='{}' )\n"
          "unknown/unsupported factorization type ('LU','QR','QRP')", s
        );
      }
    }

    //! \brief Selects LU in the reduction phase.
    void select_LU()  { m_selected = BORDERED_Choice::LU; }
    //! \brief Selects QR in the reduction phase.
    void select_QR()  { m_selected = BORDERED_Choice::QR; }
    //! \brief Selects QRP in the reduction phase.
    void select_QRP() { m_selected = BORDERED_Choice::QRP; }

    //! \brief Selects the primary last-block solver from a string.
    void
    select_last( string_view s ) {
      if      ( s == "LU"   ) select_last_LU();
      else if ( s == "LUPQ" ) select_last_LUPQ();
      else if ( s == "QR"   ) select_last_QR();
      else if ( s == "QRP"  ) select_last_QRP();
      else if ( s == "SVD"  ) select_last_SVD();
      else if ( s == "LSS"  ) select_last_LSS();
      else if ( s == "LSY"  ) select_last_LSY();
      else if ( s == "PINV" ) select_last_PINV();
      else {
        UTILS_ERROR(
          "BorderedCR::select_last( choice='{}' )\n"
          "unknown/unsupported factorization type ('LU','LUPQ','QR','QRP','SVD','LSS','LSY','PINV')", s
        );
      }
    }

    //! \brief Selects LU as the primary final solver.
    void select_last_LU()   { m_last_block_selected = BORDERED_LAST_Choice::LU;   }
    //! \brief Selects LUPQ as the primary final solver.
    void select_last_LUPQ() { m_last_block_selected = BORDERED_LAST_Choice::LUPQ; }
    //! \brief Selects QR as the primary final solver.
    void select_last_QR()   { m_last_block_selected = BORDERED_LAST_Choice::QR;   }
    //! \brief Selects QRP as the primary final solver.
    void select_last_QRP()  { m_last_block_selected = BORDERED_LAST_Choice::QRP;  }
    //! \brief Selects SVD as the primary final solver.
    void select_last_SVD()  { m_last_block_selected = BORDERED_LAST_Choice::SVD;  }
    //! \brief Selects LSS as the primary final solver.
    void select_last_LSS()  { m_last_block_selected = BORDERED_LAST_Choice::LSS;  }
    //! \brief Selects LSY as the primary final solver.
    void select_last_LSY()  { m_last_block_selected = BORDERED_LAST_Choice::LSY;  }
    //! \brief Selects PINV as the primary final solver.
    void select_last_PINV() { m_last_block_selected = BORDERED_LAST_Choice::PINV; }

    //! \brief Selects the fallback last-block solver from a string.
    void
    select_last2( string_view s ) {
      if      ( s == "NONE" ) select_last2_NONE();
      else if ( s == "SVD"  ) select_last2_SVD();
      else if ( s == "LSS"  ) select_last2_LSS();
      else if ( s == "LSY"  ) select_last2_LSY();
      else if ( s == "PINV" ) select_last2_PINV();
      else {
        UTILS_ERROR(
          "BorderedCR::select_last2( choice='{}' )\n"
          "unknown/unsupported factorization type ('NONE','SVD','LSS','LSY','PINV')", s
        );
      }
    }

    //! \brief Disables fallback on the final block.
    void select_last2_NONE() { m_last_block_selected2 = BORDERED_LAST_Choice2::NONE;  }
    //! \brief Selects SVD as final-block fallback.
    void select_last2_SVD()  { m_last_block_selected2 = BORDERED_LAST_Choice2::SVD;  }
    //! \brief Selects LSS as final-block fallback.
    void select_last2_LSS()  { m_last_block_selected2 = BORDERED_LAST_Choice2::LSS;  }
    //! \brief Selects LSY as final-block fallback.
    void select_last2_LSY()  { m_last_block_selected2 = BORDERED_LAST_Choice2::LSY;  }
    //! \brief Selects PINV as final-block fallback.
    void select_last2_PINV() { m_last_block_selected2 = BORDERED_LAST_Choice2::PINV; }

    //! \brief Returns the textual description of the reduction solver.
    static
    string
    choice_to_string( BORDERED_Choice c ) {
      string res{"none"};
      switch ( c ) {
      case BORDERED_Choice::LU:  res = "CyclicReduction+LU";         break;
      case BORDERED_Choice::QR:  res = "CyclicReduction+QR";         break;
      case BORDERED_Choice::QRP: res = "CyclicReduction+QRP";        break;
      }
      return res;
    }

    //! \brief Returns the textual description of the primary final solver.
    static
    string
    choice_to_string( BORDERED_LAST_Choice c ) {
      string res{"LastBlock not selected"};
      switch ( c ) {
      case BORDERED_LAST_Choice::LU:   res = "LastBlock LU";   break;
      case BORDERED_LAST_Choice::LUPQ: res = "LastBlock LUPQ"; break;
      case BORDERED_LAST_Choice::QR:   res = "LastBlock QR";   break;
      case BORDERED_LAST_Choice::QRP:  res = "LastBlock QRP";  break;
      case BORDERED_LAST_Choice::SVD:  res = "LastBlock SVD";  break;
      case BORDERED_LAST_Choice::LSS:  res = "LastBlock LSS";  break;
      case BORDERED_LAST_Choice::LSY:  res = "LastBlock LSY";  break;
      case BORDERED_LAST_Choice::PINV: res = "LastBlock PINV"; break;
      }
      return res;
    }

    //! \brief Returns the textual description of the fallback final solver.
    static
    string
    choice_to_string( BORDERED_LAST_Choice2 c ) {
      string res{"LastBlock not selected"};
      switch ( c ) {
      case BORDERED_LAST_Choice2::NONE: res = "LastBlock NONE"; break;
      case BORDERED_LAST_Choice2::SVD:  res = "LastBlock SVD";  break;
      case BORDERED_LAST_Choice2::LSS:  res = "LastBlock LSS";  break;
      case BORDERED_LAST_Choice2::LSY:  res = "LastBlock LSY";  break;
      case BORDERED_LAST_Choice2::PINV: res = "LastBlock PINV"; break;
      }
      return res;
    }

    //! \brief Compact description of the selected full algorithm.
    string
    info_algo() const {
      return fmt::format(
        "{} and {} and {}",
        choice_to_string(m_selected),
        choice_to_string(m_last_block_selected),
        choice_to_string(m_last_block_selected2)
      );
    }

    //! \brief Description of the reduction phase only.
    string info_algo_block()       const { return choice_to_string(m_selected); }
    //! \brief Description of the primary final solver.
    string info_algo_last_block()  const { return choice_to_string(m_last_block_selected); }
    //! \brief Description of the fallback final solver.
    string info_algo_last_block2() const { return choice_to_string(m_last_block_selected2); }

    //! \brief Returns an extended diagnostic string describing solver state.
    string info( string_view indent = "" ) const;

    //! \brief Prints the extended diagnostics to a stream.
    void info( ostream_type & stream ) const { stream << info(); }

    //!
    //! @}
    //!
    //! \brief Number of rows of the linear system
    //!
    integer
    nrows() const {
      integer const & nblk = m_number_of_blocks;
      integer const & n    = m_block_size;
      return n * (nblk+1) + m_qr + m_nr;
    }

    //!
    //! \brief Number of columns of the linear system
    //!
    integer
    ncols() const {
      integer const & nblk{ m_number_of_blocks};
      integer const & n{ m_block_size };
      return n * (nblk+1) + m_qx + m_nx;
    }

    //! \brief Number of internal blocks.
    integer number_of_blocks() const { return m_number_of_blocks; }
    //! \brief Size of the square internal block.
    integer block_size()       const { return m_block_size; }

    //! \brief Number of extra rows in the upper-lower border block.
    integer dim_qr() const { return m_qr; }
    //! \brief Number of extra columns in the lower border block.
    integer dim_qx() const { return m_qx; }
    //! \brief Number of rows in the lower `C/F` border.
    integer dim_nr() const { return m_nr; }
    //! \brief Number of columns in the right `B/F` border.
    integer dim_nx() const { return m_nx; }
    //! \brief Number of rows in the final reduced system.
    integer Nr()     const { return m_Nr; }
    //! \brief Number of columns in the final reduced system.
    integer Nc()     const { return m_Nc; }

    //!
    //! \name Filling all or part of the linear system with zero
    //! @{
    //!

    //! \brief Zeros all `D` blocks.
    void zero_D();
    //! \brief Zeros all `E` blocks.
    void zero_E();
    //! \brief Zeros all `B` blocks.
    void zero_B();
    //! \brief Zeros block `F`.
    void zero_F();
    //! \brief Zeros the upper-lower block `[H0, HN, Hq, Hp]`.
    void zero_H();
    //! \brief Zeros all `C` blocks.
    void zero_C();
    //! \brief Zeros block `Cq`.
    void zero_Cq();

    //! \brief Zeros all system blocks.
    void
    fill_zero() {
      zero_B();
      zero_C();
      zero_Cq();
      zero_D();
      zero_E();
      zero_F();
      zero_H();
    }

    //!
    //!  @}
    //!
    //!  \name Access to single block
    //!
    //!  Matrix structure
    //!
    //!  \verbatim
    //!
    //!                 n * (nblock+1)
    //!   ┌───────────────────┴──────────────────────┐
    //!    n   n   n                              n  qx  nx
    //!  ╔═══╦═══╦═══╦════.................═════╦═══╦═══╦═══╗   ─┐
    //!  ║ D ║ E ║   ║                          ║   ║   ║ B ║ n  │
    //!  ╠═══╬═══╬═══╣                     ═════╬═══╬═══╬═══╣    │
    //!  ║   ║ D ║ E ║                          ║   ║   ║ B ║ n  │
    //!  ╠═══╬═══╬═══╬═══╗                 ═════╬═══╬═══╬═══╣    │
    //!  ║   ║   ║ D ║ E ║                      ║   ║   ║ B ║ n  │
    //!  ╠═══╬═══╬═══╬═══╬                 ═════╬═══╬═══╬═══╣    │
    //!  :                                                  :    │
    //!  :                                                  :    │
    //!  :                                                  :    ├─ n * nblock
    //!  :                                                  :    │
    //!  :                                                  :    │
    //!  :                              ╔═══╦═══╦═══╦═══╦═══╣    │
    //!  :                              ║ D ║ E ║   ║   ║ B ║ n  │
    //!  :                              ╚═══╬═══╬═══╬═══╬═══╣    │
    //!  :                                  ║ D ║ E ║   ║ B ║ n  │
    //!  ╠═══╬═══╬═══................═══╬═══╬═══╬═══╬═══╬═══╣   ─┤
    //!  ║   ║   ║                          ║   ║   ║   ║   ║    │
    //!  ║H0 ║ 0 ║                          ║ 0 ║HN ║ Hq║ Hp║    │ n+qr
    //!  ║   ║   ║                          ║   ║   ║   ║   ║    │
    //!  ╠═══╬═══╬═══................═══╬═══╬═══╬═══╬═══╬═══╣   ─┤
    //!  ║ C ║ C ║                      ║ C ║ C ║ C ║ Cq║ F ║    │ nr
    //!  ╚═══╩═══╩═══................═══╩═══╩═══╩═══╩═══╩═══╝   ─┘
    //!                                             nr*qx
    //!  \endverbatim
    //!
    //!  @{
    //!

    //! \name Right-border blocks `B`
    //! @{

    //! \brief Loads block `B` associated with internal block `nbl`.
    /*\
     |  ____
     | | __ )
     | |  _ \
     | | |_) |
     | |____/
    \*/
    void load_B( integer nbl, real_type const B[], integer ldB );
    //! \brief Overload of \ref load_B(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_B( integer nbl, MatW const & B );
    //! \brief Adds a contribution to block `B` associated with `nbl`.
    void add_to_B( integer nbl, real_type const B[], integer ldB );
    //! \brief Overload of \ref add_to_B(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void add_to_B( integer nbl, MatW const & B );

    //! \brief Exports the sparse pattern of block `B` number `nbl`.
    integer pattern_B( integer nbl, integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `B` number `nbl`.
    integer values_B( integer nbl, real_type V[] ) const;

    //! @}
    //! \name Lower-border blocks `C`
    //! @{

    //! \brief Loads block `C` associated with internal block `nbl`.
    /*\
     |   ____
     |  / ___|
     | | |
     | | |___
     |  \____|
    \*/
    void load_C( integer nbl, real_type const C[], integer ldC );
    //! \brief Overload of \ref load_C(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_C( integer nbl, MatW const & C );
    //! \brief Adds a contribution to block `C` associated with `nbl`.
    void add_to_C( integer nbl, real_type const C[], integer ldC );
    //! \brief Overload of \ref add_to_C(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void add_to_C( integer nbl, MatW const & C );

    //! \brief Exports the sparse pattern of block `C` number `nbl`.
    integer pattern_C( integer nbl, integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `C` number `nbl`.
    integer values_C( integer nbl, real_type V[] ) const;

    //! \brief Adds a contribution distributed across blocks `C_nbl` and `C_{nbl+1}`.
    void add_to_C2( integer nbl, real_type const C[], integer ldC );
    //! \brief Overload of \ref add_to_C2(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void add_to_C2( integer nbl, MatW const & C );

    //! \brief Adds a joint contribution to blocks `C_{nbl}`, `C_{nbl+1}`, and `F`.
    void add_to_C2F( integer nbl, real_type const C2F[], integer ldC );
    //! \brief Overload of \ref add_to_C2F(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void add_to_C2F( integer nbl, MatW const & C2F );

    //! @}
    //! \name Diagonal blocks `D`
    //! @{

    /*\
     |  ____
     | |  _ \
     | | | | |
     | | |_| |
     | |____/
    \*/
    //! \brief Loads diagonal block `D` number `nbl`.
    void load_D( integer nbl, real_type const D[], integer ldD );
    //! \brief Overload of \ref load_D(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_D( integer nbl, MatW const & D );

    //! \brief Exports the sparse pattern of block `D` number `nbl`.
    integer pattern_D( integer nbl, integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `D` number `nbl`.
    integer values_D( integer nbl, real_type V[] ) const;

    //! @}
    //! \name Superdiagonal blocks `E`
    //! @{

    /*\
     |  _____
     | | ____|
     | |  _|
     | | |___
     | |_____|
    \*/
    //! \brief Loads superdiagonal block `E` number `nbl`.
    void load_E( integer nbl, real_type const E[], integer ldE );
    //! \brief Overload of \ref load_E(integer,const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_E( integer nbl, MatW const & E );

    //! \brief Exports the sparse pattern of block `E` number `nbl`.
    integer pattern_E( integer nbl, integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `E` number `nbl`.
    integer values_E( integer nbl, real_type V[] ) const;

    //! \brief Loads the concatenated block `[D | E]` in a single call.
    void load_DE( integer nbl, real_type const DE[], integer ldDE );
    //! \brief Loads the concatenated block `[D | E | B]` in a single call.
    void load_DEB( integer nbl, real_type const DEB[], integer ldDEB );

    //! @}
    //! \name Final lower-right block `F`
    //! @{

    /*\
     |  _____
     | |  ___|
     | | |_
     | |  _|
     | |_|
    \*/
    //! \brief Loads block `F`.
    void load_F( real_type const F[], integer ldF );
    //! \brief Overload of \ref load_F(const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_F( MatW const & F );
    //! \brief Adds a contribution to block `F`.
    void add_to_F( real_type const F[], integer ldF );
    //! \brief Overload of \ref add_to_F(const real_type*,integer) accepting a `MatrixWrapper` view.
    void add_to_F( MatW const & F );

    //! \brief Exports the sparse pattern of block `F`.
    integer pattern_F( integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `F`.
    integer values_F( real_type V[] ) const;

    //! @}
    //! \name Blocks `Cq`
    //! @{

    /*\
     |   ____
     |  / ___|__ _
     | | |   / _` |
     | | |__| (_| |
     |  \____\__, |
     |          |_|
    \*/

    //! \brief Loads block `Cq`.
    void load_Cq( real_type const Cq[], integer ldC );
    //! \brief Overload of \ref load_Cq(const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_Cq( MatW const & Cq );
    //! \brief Loads the concatenated block `[Cq | F]` in a single call.
    void load_CqF( real_type const CqF[], integer ldCF );

    //! \brief Exports the sparse pattern of block `Cq`.
    integer pattern_Cq( integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `Cq`.
    integer values_Cq( real_type V[] ) const;

    //! @}
    //! \name Upper-lower block `H`
    //! @{

    /*\
     |  _   _
     | | | | |
     | | |_| |
     | |  _  |
     | |_| |_|
    \*/
    //! \brief Exports the sparse pattern of block `H = [H0 | HN | Hq | Hp]`.
    integer pattern_H( integer I[], integer J[], integer offs ) const;
    //! \brief Exports the values of block `H`.
    integer values_H( real_type V[] ) const;

    /*!
     * \brief Loads the sub-blocks of block `H` separately.
     * \param H0 block connected to the first internal block.
     * \param ld0 leading dimension di `H0`.
     * \param HN block connected to the last internal block.
     * \param ldN leading dimension di `HN`.
     * \param Hq block associated with extra variables `qx`.
     * \param ldQ leading dimension di `Hq`.
     * \param Hp block associated with extra variables `nx`.
     * \param ldP leading dimension di `Hp`.
     */
    void
    load_bottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ,
      real_type const Hp[], integer ldP
    );

    //! \brief Overload of \ref load_bottom(const real_type*,integer,const real_type*,integer,const real_type*,integer,const real_type*,integer) accepting `MatrixWrapper` views.
    void
    load_bottom(
      MatW const & H0,
      MatW const & HN,
      MatW const & Hq,
      MatW const & Hp
    );

    //! \brief Loads block `H` already concatenated into a single dense matrix.
    void load_bottom( real_type const _H0Nqp[], integer ldH );
    //! \brief Overload of \ref load_bottom(const real_type*,integer) accepting a `MatrixWrapper` view.
    void load_bottom( MatW const & H );

    /*\
     |  ┌───┬───┬───┬───┐
     |  │ C0│ CN│ Cq│ F │
     |  └───┴───┴───┴───┘
    \*/

    /*!
     * \brief Loads the sub-blocks of the final lower border `[C0 | CN | Cq | F]` separately.
     */
    void
    load_bottom2(
      real_type const C0[], integer ld0,
      real_type const CN[], integer ldN,
      real_type const Cq[], integer ldCq,
      real_type const F[],  integer ldF
    );

    //! \brief Overload of \ref load_bottom2(const real_type*,integer,const real_type*,integer,const real_type*,integer,const real_type*,integer) accepting `MatrixWrapper` views.
    void
    load_bottom2(
      MatW const & C0,
      MatW const & CN,
      MatW const & Cq,
      MatW const & F
    );

    //! \brief Loads the final lower border already concatenated into a single matrix.
    void
    load_bottom2( MatW const & H );

    //!
    //! @}
    //!

    //!
    //!
    //! \name Access to blocks by element
    //!
    //! @{
    //!

    //! \brief Writable access to element `(i,j)` of block `B_nbl`.
    t_Value       & B ( integer nbl, integer i, integer j )       { return m_Bmat [ nbl*n_x_nx + i + j*m_block_size        ]; }
    //! \brief Read-only access to element `(i,j)` of block `B_nbl`.
    t_Value const & B ( integer nbl, integer i, integer j ) const { return m_Bmat [ nbl*n_x_nx + i + j*m_block_size        ]; }
    //! \brief Writable access to element `(i,j)` of block `C_nbl`.
    t_Value       & C ( integer nbl, integer i, integer j )       { return m_Cmat [ nbl*nr_x_n + i + j*m_nr                ]; }
    //! \brief Read-only access to element `(i,j)` of block `C_nbl`.
    t_Value const & C ( integer nbl, integer i, integer j ) const { return m_Cmat [ nbl*nr_x_n + i + j*m_nr                ]; }
    //! \brief Writable access to element `(i,j)` of block `D_nbl`.
    t_Value       & D ( integer nbl, integer i, integer j )       { return m_Dmat [ nbl*n_x_n  + i + j*m_block_size        ]; }
    //! \brief Read-only access to element `(i,j)` of block `D_nbl`.
    t_Value const & D ( integer nbl, integer i, integer j ) const { return m_Dmat [ nbl*n_x_n  + i + j*m_block_size        ]; }
    //! \brief Writable access to element `(i,j)` of block `E_nbl`.
    t_Value       & E ( integer nbl, integer i, integer j )       { return m_Emat [ nbl*n_x_n  + i + j*m_block_size        ]; }
    //! \brief Read-only access to element `(i,j)` of block `E_nbl`.
    t_Value const & E ( integer nbl, integer i, integer j ) const { return m_Emat [ nbl*n_x_n  + i + j*m_block_size        ]; }
    //! \brief Writable access to element `(i,j)` of block `F`.
    t_Value       & F (              integer i, integer j )       { return m_Fmat [              i + j*m_nr                ]; }
    //! \brief Read-only access to element `(i,j)` of block `F`.
    t_Value const & F (              integer i, integer j ) const { return m_Fmat [              i + j*m_nr                ]; }
    //! \brief Writable access to element `(i,j)` of block `Cq`.
    t_Value       & Cq(              integer i, integer j )       { return m_Cqmat[              i + j*m_nr                ]; }
    //! \brief Read-only access to element `(i,j)` of block `Cq`.
    t_Value const & Cq(              integer i, integer j ) const { return m_Cqmat[              i + j*m_nr                ]; }
    //! \brief Writable access to element `(i,j)` of block `H`.
    t_Value       & H (              integer i, integer j )       { return m_H0Nqp[              i + j*(m_block_size+m_qr) ]; }
    //! \brief Read-only access to element `(i,j)` of block `H`.
    t_Value const & H (              integer i, integer j ) const { return m_H0Nqp[              i + j*(m_block_size+m_qr) ]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //! \brief Returns a `MatrixWrapper` view of block `B_nbl`.
    void B( integer nbl, MatW & B_wrap  ) const { B_wrap.setup( m_Bmat + nbl*n_x_nx, m_block_size, m_nx, m_block_size ); }
    //! \brief Returns a `MatrixWrapper` view of block `C_nbl`.
    void C( integer nbl, MatW & C_wrap  ) const { C_wrap.setup( m_Cmat + nbl*nr_x_n, m_nr, m_block_size, m_nr ); }
    //! \brief Returns a `MatrixWrapper` view of block `D_nbl`.
    void D( integer nbl, MatW & D_wrap  ) const { D_wrap.setup( m_Dmat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }
    //! \brief Returns a `MatrixWrapper` view of block `E_nbl`.
    void E( integer nbl, MatW & E_wrap  ) const { E_wrap.setup( m_Emat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }
    //! \brief Returns a `MatrixWrapper` view of block `F`.
    void F(              MatW & F_wrap  ) const { F_wrap.setup( m_Fmat, m_nr, m_nx, m_nr ); }
    //! \brief Returns a `MatrixWrapper` view of block `Cq`.
    void Cq(             MatW & Cq_wrap ) const { Cq_wrap.setup( m_Cqmat, m_nr, m_qx, m_nr ); }
    //! \brief Returns a `MatrixWrapper` view of block `H`.
    void H(              MatW & H_wrap  ) const { H_wrap.setup( m_H0Nqp, m_block_size + m_qr, m_Nc, m_block_size + m_qr ); }

    //!
    //! @}
    //!

    /*\
     |         _      _               _
     |  __   _(_)_ __| |_ _   _  __ _| |___
     |  \ \ / / | '__| __| | | |/ _` | / __|
     |   \ V /| | |  | |_| |_| | (_| | \__ \
     |    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    //! \brief Factorizes the loaded bordered system.
    bool factorize();
    //! \brief Solves the system for a single right-hand side.
    bool solve( real_type x[] ) const override;
    //! \brief Solves the system for multiple right-hand sides.
    bool solve( integer nrhs, real_type rhs[], integer ldRhs ) const override;

    //! \brief Transposed variant not implemented for a single RHS.
    virtual bool t_solve( real_type [] ) const override;
    //! \brief Transposed variant not implemented for multiple RHS.
    virtual bool t_solve( integer, real_type [], integer ) const override;

    /*\
     |     _
     |    / \  _   ___  __
     |   / _ \| | | \ \/ /
     |  / ___ \ |_| |>  <
     | /_/   \_\__,_/_/\_\
     |
    \*/
    /*!
     * \brief Computes the matrix-vector product `res = A*x`.
     * \param x input vector.
     * \param res output buffer.
     */
    void Mv( real_type const x[], real_type res[] ) const;
    /*!
     * \brief Accumulates `alpha*A*x` into `res`.
     * \param alpha scalar coefficient.
     * \param x input vector.
     * \param res output accumulator.
     */
    void add_Mv( real_type alpha, real_type const x[], real_type res[] ) const;

    //! \brief Shorthand for `add_Mv(1,x,res)`.
    void
    add_Mv( real_type const x[], real_type res[] ) const {
      add_Mv( 1.0, x, res );
    }

    /*\
     |  ____
     | / ___| _ __   __ _ _ __ ___  ___
     | \___ \| '_ \ / _` | '__/ __|/ _ \
     |  ___) | |_) | (_| | |  \__ \  __/
     | |____/| .__/ \__,_|_|  |___/\___|
     |       |_|
     |
    \*/

    //! \brief Number of nonzero entries in the full system.
    integer
    sparse_nnz() const {
      integer const & nblock { m_number_of_blocks };
      integer const & n      { m_block_size };
      integer const & nr     { m_nr };
      integer const & nx     { m_nx };
      integer const & qr     { m_qr };
      integer const & qx     { m_qx };
      return n*nblock*(2*n+nx+nr) + nr*(n+qx+nx) + (n+qr)*(2*n+qx+nx);
    }

    //! \brief Exports the sparse pattern of the full matrix.
    void
    sparse_pattern( integer I[], integer J[], integer offs ) const;

    //! \brief Exports system values in the same order as \ref sparse_pattern.
    void
    sparse_values( real_type V[] ) const;

    /*!
     * \brief Loads the matrix from an external sparse representation.
     * \param M_values nonzero values.
     * \param M_row row indices.
     * \param r_offs offset applied to row indices.
     * \param M_col column indices.
     * \param c_offs offset applied to column indices.
     * \param M_nnz number of nonzero entries.
     * \param do_error if `true`, incompatible elements are reported as errors.
     */
    void
    sparse_load(
      real_type const M_values[],
      integer   const M_row[], integer r_offs,
      integer   const M_col[], integer c_offs,
      integer         M_nnz,
      bool            do_error = true
    );

    /*\
     |   __  __   _ _____ _      _   ___
     |  |  \/  | /_\_   _| |    /_\ | _ )
     |  | |\/| |/ _ \| | | |__ / _ \| _ \
     |  |_|  |_/_/ \_\_| |____/_/ \_\___/
    \*/

    //! \brief Prints a MATLAB/Octave script that reconstructs the matrix.
    void
    print_matlab_script( ostream_type & stream ) const;

    /*\
     |   ___        _                       _
     |  | _ \___ __| |_ __ _ _ _  __ _ _  _| |__ _ _ _
     |  |   / -_) _|  _/ _` | ' \/ _` | || | / _` | '_|
     |  |_|_\___\__|\__\__,_|_||_\__, |\_,_|_\__,_|_|
     |                           |___/
    \*/
    #if 0
    bool
    mult_inv(
      real_type const b[],
      integer         incb,
      real_type       x[],
      integer         incx
    ) const;

    bool
    t_mult_inv(
      real_type const b[],
      integer         incb,
      real_type       x[],
      integer         incx
    ) const;

    bool
    mult_inv(
      integer         nrhs,
      real_type const B[],
      integer         ldB,
      real_type       X[],
      integer         ldX
    ) const;

    bool
    t_mult_inv(
      integer         nrhs,
      real_type const B[],
      integer         ldB,
      real_type       X[],
      integer         ldX
    ) const;
    #endif

    //! \brief Enables or disables internal debugging.
    void set_debug( bool debug ) { m_debug = debug; }
    //! \brief Checks consistency and basic properties of the loaded matrix.
    void check_matrix();

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  extern template class BorderedCR<float>;
  extern template class BorderedCR<double>;
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: BABD_BorderedCR.hxx
///
