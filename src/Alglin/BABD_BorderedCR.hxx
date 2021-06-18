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

  /*\
   |   ___             _                _    ___ ___
   |  | _ ) ___ _ _ __| |___ _ _ ___ __| |  / __| _ \
   |  | _ \/ _ \ '_/ _` / -_) '_/ -_) _` | | (__|   /
   |  |___/\___/_| \__,_\___|_| \___\__,_|  \___|_|_\
  \*/

  //! Cyclic reduction of a block bidiagonal matrix
  /*!
   *
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
  \*/

  /*\
   |
   |  Matrix structure
   |
   |                 n * (nblock+1)
   |    ___________________^____________________
   |   /                                        \
   |    n   n   n                              n  qx  nx
   |  +---+---+---+----.................-----+---+---+---+   -+
   |  | D | E |   |                          |   |   | B | n  |
   |  +---+---+---+                     -----+---+---+---+    |
   |  |   | D | E |                          |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  |   |   | D | E |                      |   |   | B | n  |
   |  +---+---+---+---+                 -----+---+---+---+    |
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                                                  :     > n * nblock
   |  :                                                  :    |
   |  :                                                  :    |
   |  :                              +---+---+---+---+---+    |
   |  :                              | D | E |   |   | B | n  |
   |  :                              +---+---+---+---+---+    |
   |  :                                  | D | E |   | B | n  |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  |   |   |                          |   |   |   |   |    |
   |  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
   |  |   |   |                          |   |   |   |   |    |
   |  +---+---+---................---+---+---+---+---+---+   -+
   |  | C | C |                      | C | C | C | Cq| F |    | nr
   |  +---+---+---................---+---+---+---+---+---+   -+
   |                                             nr*qx
  \*/

  template <typename t_Value>
  class BorderedCR : public LinearSystemSolver<t_Value> {
  public:
    typedef t_Value real_type;

    typedef enum {
      BORDERED_LU      = 0,
      BORDERED_QR      = 1,
      BORDERED_QRP     = 3,
      BORDERED_SUPERLU = 4
    } BORDERED_Choice;

    typedef enum {
      BORDERED_LAST_LU   = 0,
      BORDERED_LAST_LUPQ = 1,
      BORDERED_LAST_QR   = 2,
      BORDERED_LAST_QRP  = 3,
      BORDERED_LAST_SVD  = 4,
      BORDERED_LAST_LSS  = 5,
      BORDERED_LAST_LSY  = 6,
      BORDERED_LAST_PINV = 7
    } BORDERED_LAST_Choice;

  private:

    // block copy constructor
    BorderedCR() = delete;
    BorderedCR( BorderedCR const & ) = delete;
    BorderedCR const & operator = ( BorderedCR const & ) = delete;

  protected:

    Malloc<real_type>  m_baseValue;
    Malloc<integer>    m_baseInteger;
    Malloc<real_type*> m_basePointer;
    Malloc<integer*>   m_basePointerInteger;

    Malloc<real_type>  m_superluValue;
    Malloc<int>        m_superluInteger;

    integer m_number_of_blocks; //!< total number of blocks
    integer m_block_size;       //!< size of square blocks
    integer m_qr, m_qx;         //!< extra BC
    integer m_nr, m_nx;         //!< border size

    integer m_Nr;
    integer m_Nc;
    integer m_Tsize;

    // some derived constants
    integer n_x_2;
    integer n_x_n;
    integer n_x_nx;
    integer nr_x_n;
    integer nr_x_nx;
    integer nr_x_qx;

    // for SuperLU =====================
    int * m_slu_perm_r; // row permutations from partial pivoting
    int * m_slu_perm_c; // column permutation vector
    int * m_slu_etree;

    superlu_options_t     m_slu_options;
    mutable SuperLUStat_t m_slu_stats;
    mutable SuperMatrix   m_slu_A;  // messo mutable per zittire warning
    mutable SuperMatrix   m_slu_AC; // messo mutable per zittire warning
    mutable SuperMatrix   m_slu_L;  // messo mutable per zittire warning
    mutable SuperMatrix   m_slu_U;  // messo mutable per zittire warning

    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    mutable GlobalLU_t    m_slu_glu;
    #endif

    // for SuperLU ===================== END

    BORDERED_LAST_Choice m_last_selected;
    BORDERED_Choice      m_selected;
    bool                 m_last_must_use_PINV;

    real_type * m_H0Nqp;
    real_type * m_Bmat;
    real_type * m_Cmat;
    real_type * m_Cqmat;
    real_type * m_Dmat;
    real_type * m_Emat;

    real_type** m_Fmat;
    real_type** m_WorkT;
    real_type** m_WorkQR;

    // working block
    real_type * m_Tmat;
    real_type * m_Ttau;
    real_type * m_Work;
    integer   * m_Perm;
    integer     m_Lwork;
    integer     m_LworkT;
    integer     m_LworkQR;

    // last block
    real_type * m_Hmat;

    LU<real_type>   m_last_lu;
    LUPQ<real_type> m_last_lupq;
    QR<real_type>   m_last_qr;
    QRP<real_type>  m_last_qrp;
    SVD<real_type>  m_last_svd;
    LSS<real_type>  m_last_lss;
    LSY<real_type>  m_last_lsy;
    PINV<real_type> m_last_pinv;

    integer *m_iBlock;
    integer *m_kBlock;

    // used also with a unique thread
    integer                     m_available_thread;
    integer                     m_used_thread;
    integer                     m_reduced_nblk;
    mutable integer**           m_perm_thread;
    mutable real_type**         m_xb_thread;
    mutable Utils::ThreadPool * m_TP;
    mutable Utils::SpinLock     m_spin;

    void
    buildT(
      integer         nth,
      real_type const TOP[],
      real_type const BOTTOM[],
      real_type       T[],
      integer         iperm[]
    ) const;

    void
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

    void
    applyT(
      integer         nth,
      real_type const T[],
      integer   const iperm[],
      real_type       TOP[],
      real_type       BOTTOM[]
    ) const;

    // convert permutation to exchanges
    void
    permutation_to_exchange( integer nn, integer P[], integer S[] ) const {
      for ( integer i = 0; i < nn; ++i ) {
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
    void factorize_reduced();

    /*
    //    __                            _
    //   / _|___ _ ___ __ ____ _ _ _ __| |
    //  |  _/ _ \ '_\ V  V / _` | '_/ _` |
    //  |_| \___/_|  \_/\_/\__,_|_| \__,_|
    */

    void
    forward( integer nth, real_type x[], real_type xb[] ) const;

    void
    forward_n(
      integer   nth,
      integer   nrhs,
      real_type rhs[],
      integer   ldRhs
    ) const;

    void
    forward_reduced( real_type x[], real_type xb[] ) const;

    void
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
    BorderedCR( Utils::ThreadPool * TP );

    virtual
    ~BorderedCR() override
    {}

    void
    setThreadPool( Utils::ThreadPool * TP ) {
      m_TP               = TP;
      m_available_thread = m_TP == nullptr ? 1 : m_TP->size();
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

    void
    dup( BorderedCR const & );

    /*!
     | \name Select Linar Algebra solver
     | @{
    \*/

    void select_LU()      { m_selected = BORDERED_LU; }
    void select_QR()      { m_selected = BORDERED_QR; }
    void select_QRP()     { m_selected = BORDERED_QRP; }
    void select_SUPERLU() { m_selected = BORDERED_SUPERLU; }

    void select_last_LU()   { m_last_selected = BORDERED_LAST_LU;   }
    void select_last_LUPQ() { m_last_selected = BORDERED_LAST_LUPQ; }
    void select_last_QR()   { m_last_selected = BORDERED_LAST_QR;   }
    void select_last_QRP()  { m_last_selected = BORDERED_LAST_QRP;  }
    void select_last_SVD()  { m_last_selected = BORDERED_LAST_SVD;  }
    void select_last_LSS()  { m_last_selected = BORDERED_LAST_LSS;  }
    void select_last_LSY()  { m_last_selected = BORDERED_LAST_LSY;  }
    void select_last_PINV() { m_last_selected = BORDERED_LAST_PINV; }

    static
    std::string
    choice_to_string( BORDERED_Choice c ) {
      std::string res = "none";
      switch ( c ) {
      case BORDERED_LU:      res = "CyclicReduction+LU";         break;
      case BORDERED_QR:      res = "CyclicReduction+QR";         break;
      case BORDERED_QRP:     res = "CyclicReduction+QRP";        break;
      case BORDERED_SUPERLU: res = "SuperLU(LastBlock ignored)"; break;
      }
      return res;
    }

    static
    std::string
    choice_to_string( BORDERED_LAST_Choice c ) {
      std::string res = "LastBlock not selected";
      switch ( c ) {
      case BORDERED_LAST_LU:   res = "LastBlock LU";   break;
      case BORDERED_LAST_LUPQ: res = "LastBlock LUPQ"; break;
      case BORDERED_LAST_QR:   res = "LastBlock QR";   break;
      case BORDERED_LAST_QRP:  res = "LastBlock QRP";  break;
      case BORDERED_LAST_SVD:  res = "LastBlock SVD";  break;
      case BORDERED_LAST_LSS:  res = "LastBlock LSS";  break;
      case BORDERED_LAST_LSY:  res = "LastBlock LSY";  break;
      case BORDERED_LAST_PINV: res = "LastBlock PINV"; break;
      }
      return res;
    }

    std::string
    info_algo() const {
      std::string a = choice_to_string(m_selected);
      std::string b = choice_to_string(m_last_selected);
      return a+" and "+b;
    }

    void info( ostream_type & stream ) const;

    /*!
     | @}
    \*/

    //! \brief Number of rows of the linear system
    integer
    numRows() const {
      integer const & nblk = m_number_of_blocks;
      integer const & n    = m_block_size;
      return n * (nblk+1) + m_qr + m_nr;
    }

    //! \brief Number of columns of the linear system
    integer
    numCols() const {
      integer const & nblk = m_number_of_blocks;
      integer const & n    = m_block_size;
      return n * (nblk+1) + m_qx + m_nx;
    }

    /*!
     | \name Filling all or part of the linear system with zero
     | @{
    \*/

    void zeroD();
    void zeroE();
    void zeroB();
    void zeroF();
    void zeroH();
    void zeroC();
    void zeroCq();

    void
    fillZero()
    { zeroB(); zeroC(); zeroCq(); zeroD(); zeroE(); zeroF(); zeroH(); }

    /*!
     |  @}
     |
     |  \name Access to single block
     |
     |  Matrix structure
     |
     |  \verbatim
     |
     |                 n * (nblock+1)
     |    ___________________^____________________
     |   /                                        \
     |    n   n   n                              n  qx  nx
     |  +---+---+---+----.................-----+---+---+---+   -+
     |  | D | E |   |                          |   |   | B | n  |
     |  +---+---+---+                     -----+---+---+---+    |
     |  |   | D | E |                          |   |   | B | n  |
     |  +---+---+---+---+                 -----+---+---+---+    |
     |  |   |   | D | E |                      |   |   | B | n  |
     |  +---+---+---+---+                 -----+---+---+---+    |
     |  :                                                  :    |
     |  :                                                  :    |
     |  :                                                  :     > n * nblock
     |  :                                                  :    |
     |  :                                                  :    |
     |  :                              +---+---+---+---+---+    |
     |  :                              | D | E |   |   | B | n  |
     |  :                              +---+---+---+---+---+    |
     |  :                                  | D | E |   | B | n  |
     |  +---+---+---................---+---+---+---+---+---+   -+
     |  |   |   |                          |   |   |   |   |    |
     |  |H0 | 0 |                          | 0 |HN | Hq| Hp|    | n+qr
     |  |   |   |                          |   |   |   |   |    |
     |  +---+---+---................---+---+---+---+---+---+   -+
     |  | C | C |                      | C | C | C | Cq| F |    | nr
     |  +---+---+---................---+---+---+---+---+---+   -+
     |                                             nr*qx
     |  \endverbatim
     |
     |  @{
    \*/

    // Border Right blocks
    /*\
     |  ____
     | | __ )
     | |  _ \
     | | |_) |
     | |____/
    \*/
    void loadB( integer nbl, real_type const B[], integer ldB );
    void loadB( integer nbl, MatrixWrapper<real_type> const & B );
    void addtoB( integer nbl, real_type const B[], integer ldB );
    void addtoB( integer nbl, MatrixWrapper<real_type> const & B );

    integer patternB( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesB( integer nbl, real_type V[] ) const;

    // Border Bottom blocks
    /*\
     |   ____
     |  / ___|
     | | |
     | | |___
     |  \____|
    \*/
    void loadC( integer nbl, real_type const C[], integer ldC );
    void loadC( integer nbl, MatrixWrapper<real_type> const & C );
    void addtoC( integer nbl, real_type const C[], integer ldC );
    void addtoC( integer nbl, MatrixWrapper<real_type> const & C );

    integer patternC( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesC( integer nbl, real_type V[] ) const;

    // add to block nbl and nbl+1
    void addtoC2( integer nbl, real_type const C[], integer ldC );
    void addtoC2( integer nbl, MatrixWrapper<real_type> const & C );

    // -------------------------------------------------------------------------
    /*\
     |  ____
     | |  _ \
     | | | | |
     | | |_| |
     | |____/
    \*/
    void loadD( integer nbl, real_type const D[], integer ldD );
    void loadD( integer nbl, MatrixWrapper<real_type> const & D );

    integer patternD( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesD( integer nbl, real_type V[] ) const;

    /*\
     |  _____
     | | ____|
     | |  _|
     | | |___
     | |_____|
    \*/
    void loadE( integer nbl, real_type const E[], integer ldE );
    void loadE( integer nbl, MatrixWrapper<real_type> const & E );

    integer patternE( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesE( integer nbl, real_type V[] ) const;

    void loadDE( integer nbl, real_type const DE[], integer ldDE );
    void loadDEB( integer nbl, real_type const DEB[], integer ldDEB );

    // -------------------------------------------------------------------------
    /*\
     |  _____
     | |  ___|
     | | |_
     | |  _|
     | |_|
    \*/
    void loadF( real_type const F[], integer ldF );
    void loadF( MatrixWrapper<real_type> const & F );
    void addtoF( real_type const F[], integer ldF );
    void addtoF( MatrixWrapper<real_type> const & F );

    integer patternF( integer I[], integer J[], integer offs ) const;
    integer valuesF( real_type V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |   ____
     |  / ___|__ _
     | | |   / _` |
     | | |__| (_| |
     |  \____\__, |
     |          |_|
    \*/

    void loadCq( real_type const Cq[], integer ldC );
    void loadCq( MatrixWrapper<real_type> const & Cq );
    void loadCqF( real_type const CqF[], integer ldCF );

    integer patternCq( integer I[], integer J[], integer offs ) const;
    integer valuesCq( real_type V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |  _   _
     | | | | |
     | | |_| |
     | |  _  |
     | |_| |_|
    \*/
    integer patternH( integer I[], integer J[], integer offs ) const;
    integer valuesH( real_type V[] ) const;

    void
    loadBottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ,
      real_type const Hp[], integer ldP
    );

    void
    loadBottom(
      MatrixWrapper<real_type> const & H0,
      MatrixWrapper<real_type> const & HN,
      MatrixWrapper<real_type> const & Hq,
      MatrixWrapper<real_type> const & Hp
    );

    void loadBottom( real_type const _H0Nqp[], integer ldH );
    void loadBottom( MatrixWrapper<real_type> const & H );

    /*\
     |  +---+---+---+---+
     |  | C | C | Cq| F |
     |  +---+---+---+---+
    \*/

    void
    loadBottom2(
      real_type const C0[], integer ld0,
      real_type const CN[], integer ldN,
      real_type const Cq[], integer ldCq,
      real_type const F[],  integer ldF
    );

    void
    loadBottom2(
      MatrixWrapper<real_type> const & C0,
      MatrixWrapper<real_type> const & CN,
      MatrixWrapper<real_type> const & Cq,
      MatrixWrapper<real_type> const & F
    );

    void
    loadBottom2( MatrixWrapper<real_type> const & H );

    /*!
     | @}
    \*/

    /*!
     |
     | \name Access to blocks by element
     |
     | @{
    \*/

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
    { return m_Fmat[0][ i + j*m_nr ]; }

    t_Value const &
    F( integer i, integer j ) const
    { return m_Fmat[0][ i + j*m_nr ]; }

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
    B( integer nbl, MatrixWrapper<real_type> & B_wrap )
    { B_wrap.setup( m_Bmat + nbl*n_x_nx, m_block_size, m_nx, m_block_size ); }

    void
    C( integer nbl, MatrixWrapper<real_type> & C_wrap )
    { C_wrap.setup( m_Cmat + nbl*nr_x_n, m_nr, m_block_size, m_nr ); }

    void
    D( integer nbl, MatrixWrapper<real_type> & D_wrap )
    { D_wrap.setup( m_Dmat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    E( integer nbl, MatrixWrapper<real_type> & E_wrap )
    { E_wrap.setup( m_Emat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    F( MatrixWrapper<real_type> & F_wrap )
    { F_wrap.setup( m_Fmat[0], m_nr, m_nx, m_nr ); }

    void
    Cq( MatrixWrapper<real_type> & Cq_wrap )
    { Cq_wrap.setup( m_Cqmat, m_nr, m_qx, m_nr ); }

    void
    H( MatrixWrapper<real_type> & H_wrap ) const
    { H_wrap.setup( m_H0Nqp, m_block_size + m_qr, m_Nc, m_block_size + m_qr ); }

    /*!
     | @}
    \*/

    void factorize_SuperLU();
    void factorize_CR();

    void
    factorize() {
      if ( m_selected == BORDERED_SUPERLU ) {
        this->factorize_SuperLU();
      } else {
        this->factorize_CR();
      }
    }

    bool solve_SuperLU( real_type x[] ) const;
    bool solve_CR( real_type x[] ) const;
    bool solve_SuperLU( integer nrhs, real_type rhs[], integer ldRhs ) const;
    bool solve_CR( integer nrhs, real_type rhs[], integer ldRhs ) const;

    /*\
     |         _      _               _
     |  __   _(_)_ __| |_ _   _  __ _| |___
     |  \ \ / / | '__| __| | | |/ _` | / __|
     |   \ V /| | |  | |_| |_| | (_| | \__ \
     |    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( real_type x[] ) const override {
      if ( m_selected == BORDERED_SUPERLU ) {
        return solve_SuperLU( x );
      } else {
        return solve_CR( x );
      }
    }

    virtual
    bool
    solve( integer nrhs, real_type rhs[], integer ldRhs ) const override {
      if ( m_selected == BORDERED_SUPERLU ) {
        return solve_SuperLU( nrhs, rhs, ldRhs );
      } else {
        return solve_CR( nrhs, rhs, ldRhs );
      }
    }

    virtual
    bool
    t_solve( real_type [] ) const override {
      UTILS_ERROR0( "BorderedCR::t_solve() not defined\n" );
    }

    virtual
    bool
    t_solve( integer, real_type [], integer ) const override {
      UTILS_ERROR0( "BorderedCR::t_solve() not defined\n" );
    }

    /*\
     |     _
     |    / \  _   ___  __
     |   / _ \| | | \ \/ /
     |  / ___ \ |_| |>  <
     | /_/   \_\__,_/_/\_\
     |
    \*/
    void Mv( real_type const x[], real_type res[] ) const;
    void addMv( real_type alpha, real_type const x[], real_type res[] ) const;

    void
    addMv( real_type const x[], real_type res[] ) const {
      addMv( 1.0, x, res );
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
    sparseNnz() const {
      integer const & nblock = m_number_of_blocks;
      integer const & n      = m_block_size;
      integer const & nr     = m_nr;
      integer const & nx     = m_nx;
      integer const & qr     = m_qr;
      integer const & qx     = m_qx;
      return n*nblock*(2*n+nx+nr) + nr*(n+qx+nx) + (n+qr)*(2*n+qx+nx);
    }

    void
    sparsePattern( integer I[], integer J[], integer offs ) const;

    void
    sparseValues( real_type V[] ) const;

    void
    sparseLoad(
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
    printMatlab( ostream_type & stream ) const;

    /*\
     |   ___        _                       _
     |  | _ \___ __| |_ __ _ _ _  __ _ _  _| |__ _ _ _
     |  |   / -_) _|  _/ _` | ' \/ _` | || | / _` | '_|
     |  |_|_\___\__|\__\__,_|_||_\__, |\_,_|_\__,_|_|
     |                           |___/
    \*/
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

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class BorderedCR<float>;
  extern template class BorderedCR<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: BABD_BorderedCR.hxx
///
