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
    typedef t_Value valueType;

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

    Malloc<valueType> m_baseValue;
    Malloc<integer>   m_baseInteger;

    Malloc<valueType> m_superluValue;
    Malloc<int>       m_superluInteger;

    integer m_number_of_blocks; //!< total number of blocks
    integer m_block_size;       //!< size of square blocks
    integer m_qr, m_qx;         //!< extra BC
    integer m_nr, m_nx;         //!< border size

    // some derived constants
    integer n_x_2;
    integer n_x_n;
    integer n_x_nx;
    integer nr_x_n;
    integer nr_x_nx;
    integer nr_x_qx;
    integer Nr;
    integer Nc;
    integer Tsize;

    // for SuperLU =====================
    int * slu_perm_r; // row permutations from partial pivoting
    int * slu_perm_c; // column permutation vector
    int * slu_etree;

    superlu_options_t     slu_options;
    mutable SuperLUStat_t slu_stats;
    mutable SuperMatrix   slu_A;  // messo mutable per zittire warning
    mutable SuperMatrix   slu_AC; // messo mutable per zittire warning
    mutable SuperMatrix   slu_L;  // messo mutable per zittire warning
    mutable SuperMatrix   slu_U;  // messo mutable per zittire warning

    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    mutable GlobalLU_t    slu_glu;
    #endif

    // for SuperLU ===================== END

    BORDERED_LAST_Choice m_last_selected;
    BORDERED_Choice      m_selected;

    void
    buildT(
      integer         nth,
      valueType const TOP[],
      valueType const BOTTOM[],
      valueType       T[],
      integer         iperm[]
    ) const;

    void
    applyT(
      integer         nth,
      valueType const T[],
      integer   const iperm[],
      valueType       TOP[],
      integer         ldTOP,
      valueType       BOTTOM[],
      integer         ldBOTTOM,
      integer         ncol
    ) const;

    void
    applyT(
      integer         nth,
      valueType const T[],
      integer   const iperm[],
      valueType       TOP[],
      valueType       BOTTOM[]
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
    forward( integer nth, valueType x[], valueType xb[] ) const;

    void
    forward_n(
      integer   nth,
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    void
    forward_reduced( valueType x[], valueType xb[] ) const;

    void
    forward_n_reduced(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    /*
    //   _             _                       _
    //  | |__  __ _ __| |____ __ ____ _ _ _ __| |
    //  | '_ \/ _` / _| / /\ V  V / _` | '_/ _` |
    //  |_.__/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    */

    void backward( integer nth, valueType x[] ) const;
    void backward_reduced( valueType x[] ) const;

    void
    backward_n(
      integer   nth,
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    void
    backward_n_reduced(
      integer   nrhs,
      valueType rhs[],
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
    solve_last( valueType [] ) const;

    bool
    solve_last(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const;

    valueType * m_H0Nqp;
    valueType * m_Bmat;
    valueType * m_Cmat;
    valueType * m_Cqmat;
    valueType * m_Dmat;
    valueType * m_Emat;

    std::vector<valueType*> m_Fmat;
    std::vector<valueType*> m_WorkT;
    std::vector<valueType*> m_WorkQR;

    // working block
    valueType * m_Tmat;
    valueType * m_Ttau;
    valueType * m_Work;
    integer   * m_Perm;
    integer     m_Lwork;
    integer     m_LworkT;
    integer     m_LworkQR;

    // last block
    valueType * m_Hmat;

    LU<valueType>   last_lu;
    LUPQ<valueType> last_lupq;
    QR<valueType>   last_qr;
    QRP<valueType>  last_qrp;
    SVD<valueType>  last_svd;
    LSS<valueType>  last_lss;
    LSY<valueType>  last_lsy;
    PINV<valueType> last_pinv;

    integer *iBlock, *kBlock;

    // used also with a unique thread
    integer                         m_usedThread;
    mutable std::vector<integer*>   m_perm_thread;
    mutable std::vector<valueType*> m_xb_thread;
    mutable Utils::ThreadPool *     m_TP;
    mutable Utils::SpinLock         m_spin;

  public:

    using LinearSystemSolver<t_Value>::factorize;

    explicit
    BorderedCR( Utils::ThreadPool * _TP = nullptr );

    virtual
    ~BorderedCR() UTILS_OVERRIDE
    {}

    void
    setThreadPool( Utils::ThreadPool * _TP = nullptr ) {
      m_TP         = _TP;
      m_usedThread = m_TP == nullptr ? 1 : m_TP->size();
    }

    //! load matrix in the class
    void
    allocate(
      integer _nblock,
      integer _n,
      integer _qr,
      integer _qx,
      integer _nr,
      integer _nx
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
      case BORDERED_LAST_LU:   res = "LastBlock LU";    break;
      case BORDERED_LAST_LUPQ: res = "LastBlock LUPQ";  break;
      case BORDERED_LAST_QR:   res = "LastBlock QR";    break;
      case BORDERED_LAST_QRP:  res = "LastBlock QRP";   break;
      case BORDERED_LAST_SVD:  res = "LastBlock SVD";   break;
      case BORDERED_LAST_LSS:  res = "LastBlock LSS";   break;
      case BORDERED_LAST_LSY:  res = "LastBlock LSY";   break;
      case BORDERED_LAST_PINV: res = "LastBlock PINV";  break;
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
      return n * (nblk+1) + m_qx + m_nx;
    }

    //! \brief Number of columns of the linear system
    integer
    numCols() const {
      integer const & nblk = m_number_of_blocks;
      integer const & n    = m_block_size;
      return n * (nblk+1) + m_qr + m_nr;
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
    void loadB( integer nbl, valueType const B[], integer ldB );
    void loadB( integer nbl, MatrixWrapper<valueType> const & B );
    void addtoB( integer nbl, valueType const B[], integer ldB );
    void addtoB( integer nbl, MatrixWrapper<valueType> const & B );

    integer patternB( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesB( integer nbl, valueType V[] ) const;

    // Border Bottom blocks
    /*\
     |   ____
     |  / ___|
     | | |
     | | |___
     |  \____|
    \*/
    void loadC( integer nbl, valueType const C[], integer ldC );
    void loadC( integer nbl, MatrixWrapper<valueType> const & C );
    void addtoC( integer nbl, valueType const C[], integer ldC );
    void addtoC( integer nbl, MatrixWrapper<valueType> const & C );

    integer patternC( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesC( integer nbl, valueType V[] ) const;

    // add to block nbl and nbl+1
    void addtoC2( integer nbl, valueType const C[], integer ldC );
    void addtoC2( integer nbl, MatrixWrapper<valueType> const & C );

    // -------------------------------------------------------------------------
    /*\
     |  ____
     | |  _ \
     | | | | |
     | | |_| |
     | |____/
    \*/
    void loadD( integer nbl, valueType const D[], integer ldD );
    void loadD( integer nbl, MatrixWrapper<valueType> const & D );

    integer patternD( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesD( integer nbl, valueType V[] ) const;

    /*\
     |  _____
     | | ____|
     | |  _|
     | | |___
     | |_____|
    \*/
    void loadE( integer nbl, valueType const E[], integer ldE );
    void loadE( integer nbl, MatrixWrapper<valueType> const & E );

    integer patternE( integer nbl, integer I[], integer J[], integer offs ) const;
    integer valuesE( integer nbl, valueType V[] ) const;

    void loadDE( integer nbl, valueType const DE[], integer ldDE );
    void loadDEB( integer nbl, valueType const DEB[], integer ldDEB );

    // -------------------------------------------------------------------------
    /*\
     |  _____
     | |  ___|
     | | |_
     | |  _|
     | |_|
    \*/
    void loadF( valueType const F[], integer ldF );
    void loadF( MatrixWrapper<valueType> const & F );
    void addtoF( valueType const F[], integer ldF );
    void addtoF( MatrixWrapper<valueType> const & F );

    integer patternF( integer I[], integer J[], integer offs ) const;
    integer valuesF( valueType V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |   ____
     |  / ___|__ _
     | | |   / _` |
     | | |__| (_| |
     |  \____\__, |
     |          |_|
    \*/

    void loadCq( valueType const Cq[], integer ldC );
    void loadCq( MatrixWrapper<valueType> const & Cq );
    void loadCqF( valueType const CqF[], integer ldCF );

    integer patternCq( integer I[], integer J[], integer offs ) const;
    integer valuesCq( valueType V[] ) const;

    // -------------------------------------------------------------------------
    /*\
     |  _   _
     | | | | |
     | | |_| |
     | |  _  |
     | |_| |_|
    \*/
    integer patternH( integer I[], integer J[], integer offs ) const;
    integer valuesH( valueType V[] ) const;

    void
    loadBottom(
      valueType const H0[], integer ld0,
      valueType const HN[], integer ldN,
      valueType const Hq[], integer ldQ,
      valueType const Hp[], integer ldP
    );

    void
    loadBottom(
      MatrixWrapper<valueType> const & H0,
      MatrixWrapper<valueType> const & HN,
      MatrixWrapper<valueType> const & Hq,
      MatrixWrapper<valueType> const & Hp
    );

    void loadBottom( valueType const _H0Nqp[], integer ldH );
    void loadBottom( MatrixWrapper<valueType> const & H );

    /*\
     |  +---+---+---+---+
     |  | C | C | Cq| F |
     |  +---+---+---+---+
    \*/

    void
    loadBottom2(
      valueType const C0[], integer ld0,
      valueType const CN[], integer ldN,
      valueType const Cq[], integer ldCq,
      valueType const F[],  integer ldF
    );

    void
    loadBottom2(
      MatrixWrapper<valueType> const & C0,
      MatrixWrapper<valueType> const & CN,
      MatrixWrapper<valueType> const & Cq,
      MatrixWrapper<valueType> const & F
    );

    void
    loadBottom2( MatrixWrapper<valueType> const & H );

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
    B( integer nbl, MatrixWrapper<valueType> & B_wrap )
    { B_wrap.setup( m_Bmat + nbl*n_x_nx, m_block_size, m_nx, m_block_size ); }

    void
    C( integer nbl, MatrixWrapper<valueType> & C_wrap )
    { C_wrap.setup( m_Cmat + nbl*nr_x_n, m_nr, m_block_size, m_nr ); }

    void
    D( integer nbl, MatrixWrapper<valueType> & D_wrap )
    { D_wrap.setup( m_Dmat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    E( integer nbl, MatrixWrapper<valueType> & E_wrap )
    { E_wrap.setup( m_Emat + nbl*n_x_n, m_block_size, m_block_size, m_block_size ); }

    void
    F( MatrixWrapper<valueType> & F_wrap )
    { F_wrap.setup( m_Fmat[0], m_nr, m_nx, m_nr ); }

    void
    Cq( MatrixWrapper<valueType> & Cq_wrap )
    { Cq_wrap.setup( m_Cqmat, m_nr, m_qx, m_nr ); }

    void
    H( MatrixWrapper<valueType> & H_wrap ) const
    { H_wrap.setup( m_H0Nqp, m_block_size + m_qr, Nc, m_block_size + m_qr ); }

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

    bool solve_SuperLU( valueType x[] ) const;
    bool solve_CR( valueType x[] ) const;
    bool solve_SuperLU( integer nrhs, valueType rhs[], integer ldRhs ) const;
    bool solve_CR( integer nrhs, valueType rhs[], integer ldRhs ) const;

    /*\
     |         _      _               _
     |  __   _(_)_ __| |_ _   _  __ _| |___
     |  \ \ / / | '__| __| | | |/ _` | / __|
     |   \ V /| | |  | |_| |_| | (_| | \__ \
     |    \_/ |_|_|   \__|\__,_|\__,_|_|___/
    \*/

    virtual
    bool
    solve( valueType x[] ) const UTILS_OVERRIDE {
      if ( m_selected == BORDERED_SUPERLU ) {
        return solve_SuperLU( x );
      } else {
        return solve_CR( x );
      }
    }

    virtual
    bool
    solve( integer nrhs, valueType rhs[], integer ldRhs ) const UTILS_OVERRIDE {
      if ( m_selected == BORDERED_SUPERLU ) {
        return solve_SuperLU( nrhs, rhs, ldRhs );
      } else {
        return solve_CR( nrhs, rhs, ldRhs );
      }
    }

    virtual
    bool
    t_solve( valueType [] ) const UTILS_OVERRIDE {
      UTILS_ERROR0( "BorderedCR::t_solve() not defined\n" );
    }

    virtual
    bool
    t_solve( integer, valueType [], integer ) const UTILS_OVERRIDE {
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
    void Mv( valueType const x[], valueType res[] ) const;
    void addMv( valueType const x[], valueType res[] ) const;

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
      return nblock*(2*n_x_n+n_x_nx+nr_x_n) + nr_x_n +
             m_nr*(m_qx+m_nx) + (n+m_qr)*(2*n+m_qx+m_nx);
    }

    void
    sparsePattern( integer I[], integer J[], integer offs ) const;

    void
    sparseValues( valueType V[] ) const;

    void
    sparseLoad(
      valueType const M_values[],
      integer   const M_row[], integer r_offs,
      integer   const M_col[], integer c_offs,
      integer         M_nnz
    );

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
