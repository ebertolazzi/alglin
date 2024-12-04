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
/// file: BABD_BorderedCR.hxx
///

namespace alglin {

  using std::atomic;
  using std::swap;
  using std::min;
  using std::max;

  /*\
   |   ___             _                _    ___ ___
   |  | _ ) ___ _ _ __| |___ _ _ ___ __| |  / __| _ \
   |  | _ \/ _ \ '_/ _` / -_) '_/ -_) _` | | (__|   /
   |  |___/\___/_| \__,_\___|_| \___\__,_|  \___|_|_\
  \*/

  //!
  //! Cyclic reduction of a block bidiagonal matrix
  //!
  //!
  //! \date     October 25, 2016
  //! \version  1.0
  //! \note     October 25, 2016
  //!
  //! \author   Enrico Bertolazzi
  //!
  //! \par      Affiliation:
  //!           Department of Industrial Engineering<br>
  //!           University of Trento <br>
  //!           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
  //!           enrico.bertolazzi\@unitn.it
  //!
  //!
  //!
  //!  Matrix structure
  //!
  //!  \verbatim
  //!                 n * (nblock+1)
  //!    ___________________^____________________
  //!   /                                        \
  //!    n   n   n                              n  qx  nx
  //!  +---+---+---+----.................-----+---+---+---+   -+
  //!  | D | E |   |                          |   |   | B | n  |
  //!  +---+---+---+                     -----+---+---+---+    |
  //!  |   | D | E |                          |   |   | B | n  |
  //!  +---+---+---+---+                 -----+---+---+---+    |
  //!  |   |   | D | E |                      |   |   | B | n  |
  //!  +---+---+---+---+                 -----+---+---+---+    |
  //!  :                                                  :    |
  //!  :                                                  :    |
  //!  :                                                  :     > n * nblock
  //!  :                                                  :    |
  //!  :                                                  :    |
  //!  :                              +---+---+---+---+---+    |
  //!  :                              | D | E |   |   | B | n  |
  //!  :                              +---+---+---+---+---+    |
  //!  :                                  | D | E |   | B | n  |
  //!  +---+---+---................---+---+---+---+---+---+   -+
  //!  |   |   |                          |   |   |   |   |    |
  //!  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
  //!  |   |   |                          |   |   |   |   |    |
  //!  +---+---+---................---+---+---+---+---+---+   -+
  //!  | C | C |                      | C | C | C | Cq| F |    | nr
  //!  +---+---+---................---+---+---+---+---+---+   -+
  //!                                             nr*qx
  //!
  //!  \endverbatim
  //!
  template <typename t_Value>
  class BorderedCR : public LinearSystemSolver<t_Value> {
  public:
    using real_type       = t_Value;
    using MatW            = MatrixWrapper<t_Value>;
    using BORDERED_Choice = enum class BORDERED_Choice : integer {
      LU  = 0,
      QR  = 1,
      QRP = 3
    };

    using BORDERED_LAST_Choice = enum class BORDERED_LAST_Choice : integer {
      LU   = 0,
      LUPQ = 1,
      QR   = 2,
      QRP  = 3,
      SVD  = 4,
      LSS  = 5,
      LSY  = 6,
      PINV = 7
    };

    // block copy constructor
    BorderedCR() = delete;
    BorderedCR( BorderedCR const & ) = delete;
    BorderedCR const & operator = ( BorderedCR const & ) = delete;

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

    BORDERED_LAST_Choice m_last_selected{BORDERED_LAST_Choice::LU};
    BORDERED_Choice      m_selected{BORDERED_Choice::LU};
    bool                 m_last_can_use_PINV{false};
    bool                 m_last_use_PINV{false};

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

    Utils::ThreadPool0      m_TP_fake{0};
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

    void set_factorize_use_thread( bool yes_no ) { m_factorize_use_thread = yes_no; }
    void set_solve_use_thread( bool yes_no )     { m_solve_use_thread     = yes_no; }
    void set_can_use_pinv( bool yes_no )         { m_last_can_use_PINV    = yes_no; }
    bool can_use_pinv() const                    { return m_last_can_use_PINV; }

    bool factorize_use_thread() const { return m_factorize_use_thread; }
    bool solve_use_thread()     const { return m_solve_use_thread; }

    void
    set_num_parallel_block( integer num_parallel_block ) {
      // taglia il massimo numero di blocchi paralleli al numero delle thread
      m_max_parallel_block = m_TP->thread_count();
      if ( m_max_parallel_block > num_parallel_block )
        m_max_parallel_block = num_parallel_block;
    }

    //! load matrix in the class
    void
    allocate(
      integer nblock,
      integer n,
      integer qr,
      integer qx,
      integer nr,
      integer nx
    );

    string const & last_error() const { return m_last_error; }

    void dup( BorderedCR const & );

    //!
    //! \name Select Linar Algebra solver
    //! @{
    //!

    void select_LU()  { m_selected = BORDERED_Choice::LU; }
    void select_QR()  { m_selected = BORDERED_Choice::QR; }
    void select_QRP() { m_selected = BORDERED_Choice::QRP; }

    void select_last_LU()   { m_last_selected = BORDERED_LAST_Choice::LU;   }
    void select_last_LUPQ() { m_last_selected = BORDERED_LAST_Choice::LUPQ; }
    void select_last_QR()   { m_last_selected = BORDERED_LAST_Choice::QR;   }
    void select_last_QRP()  { m_last_selected = BORDERED_LAST_Choice::QRP;  }
    void select_last_SVD()  { m_last_selected = BORDERED_LAST_Choice::SVD;  }
    void select_last_LSS()  { m_last_selected = BORDERED_LAST_Choice::LSS;  }
    void select_last_LSY()  { m_last_selected = BORDERED_LAST_Choice::LSY;  }
    void select_last_PINV() { m_last_selected = BORDERED_LAST_Choice::PINV; }

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

    string
    info_algo() const {
      return choice_to_string(m_selected)+
             " and "+
             choice_to_string(m_last_selected);
    }

    string info_algo_block()      const { return choice_to_string(m_selected); }
    string info_algo_last_block() const { return choice_to_string(m_last_selected); }

    string info( char const * indent = "" ) const;

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

    integer number_of_blocks() const { return m_number_of_blocks; }
    integer block_size()       const { return m_block_size; }

    integer dim_qr() const { return m_qr; }
    integer dim_qx() const { return m_qx; }
    integer dim_nr() const { return m_nr; }
    integer dim_nx() const { return m_nx; }
    integer Nr()     const { return m_Nr; }
    integer Nc()     const { return m_Nc; }

    //!
    //! \name Filling all or part of the linear system with zero
    //! @{
    //!

    void zero_D();
    void zero_E();
    void zero_B();
    void zero_F();
    void zero_H();
    void zero_C();
    void zero_Cq();

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
    //!    ___________________^____________________
    //!   /                                        \
    //!    n   n   n                              n  qx  nx
    //!  +---+---+---+----.................-----+---+---+---+   -+
    //!  | D | E |   |                          |   |   | B | n  |
    //!  +---+---+---+                     -----+---+---+---+    |
    //!  |   | D | E |                          |   |   | B | n  |
    //!  +---+---+---+---+                 -----+---+---+---+    |
    //!  |   |   | D | E |                      |   |   | B | n  |
    //!  +---+---+---+---+                 -----+---+---+---+    |
    //!  :                                                  :    |
    //!  :                                                  :    |
    //!  :                                                  :     > n * nblock
    //!  :                                                  :    |
    //!  :                                                  :    |
    //!  :                              +---+---+---+---+---+    |
    //!  :                              | D | E |   |   | B | n  |
    //!  :                              +---+---+---+---+---+    |
    //!  :                                  | D | E |   | B | n  |
    //!  +---+---+---................---+---+---+---+---+---+   -+
    //!  |   |   |                          |   |   |   |   |    |
    //!  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
    //!  |   |   |                          |   |   |   |   |    |
    //!  +---+---+---................---+---+---+---+---+---+   -+
    //!  | C | C |                      | C | C | C | Cq| F |    | nr
    //!  +---+---+---................---+---+---+---+---+---+   -+
    //!                                             nr*qx
    //!  \endverbatim
    //!
    //!  @{
    //!

    // Border Right blocks
    /*\
     |  ____
     | | __ )
     | |  _ \
     | | |_) |
     | |____/
    \*/
    void load_B( integer nbl, real_type const B[], integer ldB );
    void load_B( integer nbl, MatW const & B );
    void add_to_B( integer nbl, real_type const B[], integer ldB );
    void add_to_B( integer nbl, MatW const & B );

    integer pattern_B( integer nbl, integer I[], integer J[], integer offs ) const;
    integer values_B( integer nbl, real_type V[] ) const;

    // Border Bottom blocks
    /*\
     |   ____
     |  / ___|
     | | |
     | | |___
     |  \____|
    \*/
    void load_C( integer nbl, real_type const C[], integer ldC );
    void load_C( integer nbl, MatW const & C );
    void add_to_C( integer nbl, real_type const C[], integer ldC );
    void add_to_C( integer nbl, MatW const & C );

    integer pattern_C( integer nbl, integer I[], integer J[], integer offs ) const;
    integer values_C( integer nbl, real_type V[] ) const;

    // add to block nbl and nbl+1
    void add_to_C2( integer nbl, real_type const C[], integer ldC );
    void add_to_C2( integer nbl, MatW const & C );

    void add_to_C2F( integer nbl, real_type const C2F[], integer ldC );
    void add_to_C2F( integer nbl, MatW const & C2F );

    // -------------------------------------------------------------------------
    /*\
     |  ____
     | |  _ \
     | | | | |
     | | |_| |
     | |____/
    \*/
    void load_D( integer nbl, real_type const D[], integer ldD );
    void load_D( integer nbl, MatW const & D );

    integer pattern_D( integer nbl, integer I[], integer J[], integer offs ) const;
    integer values_D( integer nbl, real_type V[] ) const;

    /*\
     |  _____
     | | ____|
     | |  _|
     | | |___
     | |_____|
    \*/
    void load_E( integer nbl, real_type const E[], integer ldE );
    void load_E( integer nbl, MatW const & E );

    integer pattern_E( integer nbl, integer I[], integer J[], integer offs ) const;
    integer values_E( integer nbl, real_type V[] ) const;

    void load_DE( integer nbl, real_type const DE[], integer ldDE );
    void load_DEB( integer nbl, real_type const DEB[], integer ldDEB );

    // -------------------------------------------------------------------------
    /*\
     |  _____
     | |  ___|
     | | |_
     | |  _|
     | |_|
    \*/
    void load_F( real_type const F[], integer ldF );
    void load_F( MatW const & F );
    void add_to_F( real_type const F[], integer ldF );
    void add_to_F( MatW const & F );

    integer pattern_F( integer I[], integer J[], integer offs ) const;
    integer values_F( real_type V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |   ____
     |  / ___|__ _
     | | |   / _` |
     | | |__| (_| |
     |  \____\__, |
     |          |_|
    \*/

    void load_Cq( real_type const Cq[], integer ldC );
    void load_Cq( MatW const & Cq );
    void load_CqF( real_type const CqF[], integer ldCF );

    integer pattern_Cq( integer I[], integer J[], integer offs ) const;
    integer values_Cq( real_type V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |  _   _
     | | | | |
     | | |_| |
     | |  _  |
     | |_| |_|
    \*/
    integer pattern_H( integer I[], integer J[], integer offs ) const;
    integer values_H( real_type V[] ) const;

    void
    load_bottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ,
      real_type const Hp[], integer ldP
    );

    void
    load_bottom(
      MatW const & H0,
      MatW const & HN,
      MatW const & Hq,
      MatW const & Hp
    );

    void load_bottom( real_type const _H0Nqp[], integer ldH );
    void load_bottom( MatW const & H );

    /*\
     |  +---+---+---+---+
     |  | C0| CN| Cq| F |
     |  +---+---+---+---+
    \*/

    void
    load_bottom2(
      real_type const C0[], integer ld0,
      real_type const CN[], integer ldN,
      real_type const Cq[], integer ldCq,
      real_type const F[],  integer ldF
    );

    void
    load_bottom2(
      MatW const & C0,
      MatW const & CN,
      MatW const & Cq,
      MatW const & F
    );

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

    t_Value &
    B( integer nbl, integer i, integer j )
    { return m_Bmat[ nbl*n_x_nx + i + j*m_block_size ]; }

    t_Value const &
    B( integer nbl, integer i, integer j ) const
    { return m_Bmat[ nbl*n_x_nx + i + j*m_block_size ]; }

    t_Value &
    C( integer nbl, integer i, integer j )
    { return m_Cmat[ nbl*nr_x_n + i + j*m_nr ]; }

    t_Value const &
    C( integer nbl, integer i, integer j ) const
    { return m_Cmat[ nbl*nr_x_n + i + j*m_nr ]; }

    t_Value &
    D( integer nbl, integer i, integer j )
    { return m_Dmat[ nbl*n_x_n + i + j*m_block_size]; }

    t_Value const &
    D( integer nbl, integer i, integer j ) const
    { return m_Dmat[ nbl*n_x_n + i + j*m_block_size]; }

    t_Value &
    E( integer nbl, integer i, integer j )
    { return m_Emat[ nbl*n_x_n + i + j*m_block_size]; }

    t_Value const &
    E( integer nbl, integer i, integer j ) const
    { return m_Emat[ nbl*n_x_n + i + j*m_block_size]; }

    t_Value &
    F( integer i, integer j )
    { return m_Fmat[ i + j*m_nr ]; }

    t_Value const &
    F( integer i, integer j ) const
    { return m_Fmat[ i + j*m_nr ]; }

    t_Value &
    Cq( integer i, integer j )
    { return m_Cqmat[ i + j*m_nr ]; }

    t_Value const &
    Cq( integer i, integer j ) const
    { return m_Cqmat[ i + j*m_nr ]; }

    t_Value &
    H( integer i, integer j )
    { return m_H0Nqp[ i + j*(m_block_size+m_qr) ]; }

    t_Value const &
    H( integer i, integer j ) const
    { return m_H0Nqp[ i + j*(m_block_size+m_qr) ]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    B( integer nbl, MatW & B_wrap )
    { B_wrap.setup( m_Bmat + nbl*n_x_nx, m_block_size, m_nx, m_block_size ); }

    void
    C( integer nbl, MatW & C_wrap )
    { C_wrap.setup( m_Cmat + nbl*nr_x_n, m_nr, m_block_size, m_nr ); }

    void
    D( integer nbl, MatW & D_wrap )
    { D_wrap.setup( m_Dmat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    E( integer nbl, MatW & E_wrap )
    { E_wrap.setup( m_Emat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    F( MatW & F_wrap )
    { F_wrap.setup( m_Fmat, m_nr, m_nx, m_nr ); }

    void
    Cq( MatW & Cq_wrap )
    { Cq_wrap.setup( m_Cqmat, m_nr, m_qx, m_nr ); }

    void
    H( MatW & H_wrap ) const
    { H_wrap.setup( m_H0Nqp, m_block_size + m_qr, m_Nc, m_block_size + m_qr ); }

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

    bool factorize();
    bool solve( real_type x[] ) const override;
    bool solve( integer nrhs, real_type rhs[], integer ldRhs ) const override;

    virtual bool t_solve( real_type [] ) const override;
    virtual bool t_solve( integer, real_type [], integer ) const override;

    /*\
     |     _
     |    / \  _   ___  __
     |   / _ \| | | \ \/ /
     |  / ___ \ |_| |>  <
     | /_/   \_\__,_/_/\_\
     |
    \*/
    void Mv( real_type const x[], real_type res[] ) const;
    void add_Mv( real_type alpha, real_type const x[], real_type res[] ) const;

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

    void
    sparse_pattern( integer I[], integer J[], integer offs ) const;

    void
    sparse_values( real_type V[] ) const;

    void
    sparse_load(
      real_type const M_values[],
      integer   const M_row[], integer r_offs,
      integer   const M_col[], integer c_offs,
      integer         M_nnz
    );

    /*\
     |   __  __   _ _____ _      _   ___
     |  |  \/  | /_\_   _| |    /_\ | _ )
     |  | |\/| |/ _ \| | | |__ / _ \| _ \
     |  |_|  |_/_/ \_\_| |____/_/ \_\___/
    \*/

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

    void set_debug( bool debug ) { m_debug = debug; }
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
