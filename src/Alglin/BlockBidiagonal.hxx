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
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: BlockBidiagonal.hxx
///

namespace alglin {

  /*\
   |  ____  _            _      ____  _     _ _                               _
   | | __ )| | ___   ___| | __ | __ )(_) __| (_) __ _  __ _  ___  _ __   __ _| |
   | |  _ \| |/ _ \ / __| |/ / |  _ \| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
   | | |_) | | (_) | (__|   <  | |_) | | (_| | | (_| | (_| | (_) | | | | (_| | |
   | |____/|_|\___/ \___|_|\_\ |____/|_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
   |                                                  |___/
  \*/

  //!
  //!  Cyclic reduction of a block bidiagonal matrix
  //!
  //!  \date     October 25, 2016
  //!  \version  1.0
  //!  \note     October 25, 2016
  //!
  //!  \author   Enrico Bertolazzi
  //!
  //!  \par      Affiliation:
  //!            Department of Industrial Engineering<br>
  //!            University of Trento <br>
  //!            Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
  //!            enrico.bertolazzi\@unitn.it
  //!
  template <typename t_Value>
  class BlockBidiagonal {
  public:

    typedef t_Value real_type;

    //! available LU factorization code
    using BB_LASTBLOCK_Choice = enum class BB_LASTBLOCK_Choice : integer {
      LU   = 0,
      LUPQ = 1,
      QR   = 2,
      QRP  = 3,
      SVD  = 4,
      LSS  = 5,
      LSY  = 6,
      PINV = 7
    };

    static
    string
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
    //  +-----+-----+-----+----.........-----+-----+    \
    //  |  D  |  E  |  0  |                  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+     |
    //  |  0  |  D  |  E  |  0               |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |
    //  |  0  |  0  |  D  |  E  |            |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |
    //  |                                                |
    //  :                                                 > n * nblock
    //  :                                                |
    //  :                                                |
    //  :                                                |
    //  :                              +-----+-----+     |
    //  :                              |  E  |  0  |     |
    //  :                        +-----+-----+-----+     |
    //  :                        |  0  |  D  |  E  | n   |
    //  +-----+-----+---......---+-----+-----+=====+--+  /
    //  |     |                              |     |  |  \
    //  | H0  |                              | HN  |Hq|  |
    //  |     |                              |     |  |  | n+q
    //  +-----+-----+---......---+-----+-----+=====+--+  /
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

    //! allocatew and resize the problem
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

    // filling bidiagonal part of the matrix
    void load_blocks( real_type const AdAu[], integer ldA );
    void load_block( integer nbl, real_type const AdAu[], integer ldA );

    void load_block_left( integer nbl, real_type const Ad[], integer ldA );
    void load_block_right( integer nbl, real_type const Au[], integer ldA );

    // Border Bottom blocks
    void set_zero_bottom_blocks();
    void load_bottom_blocks( real_type const C[], integer ldC );
    void load_bottom_block( integer nbl, real_type const C[], integer ldC );
    void add_to_bottom_block( integer nbl, real_type const C[], integer ldC );

    // add to bottom block nbl and nbl+1
    void add_to_bottom_block2( integer nbl, real_type const C[], integer ldC );
    void load_bottom_last_block( real_type const C[], integer ldC );

    // Border Right blocks
    void set_zero_right_blocks();
    void load_right_blocks( real_type const B[], integer ldB );
    void loadRightBlock( integer nbl, real_type const B[], integer ldB );
    void loadRightLastBlock( real_type const B[], integer ldB );

    // Border RBblock
    void setZeroRBblock();
    void load_RB_block( real_type const D[], integer ldD );

    // final blocks after cyclic reduction
    real_type const * getPointer_LR() const { return m_DE_blk; }

    void getBlock_LR( real_type LR[], integer ldA ) const;
    void getBlock_L ( real_type L[],  integer ldA ) const;
    void getBlock_R ( real_type R[],  integer ldA ) const;
    void getBlock_H0( real_type H0[], integer ld0 ) const;
    void getBlock_HN( real_type HN[], integer ldN ) const;
    void getBlock_Hq( real_type Hq[], integer ldQ ) const;

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

    virtual
    void
    factorize() = 0;

    virtual
    void
    solve( real_type [] ) const = 0;

    virtual
    void
    solve(
      integer      /* nrhs  */,
      real_type [] /* rhs   */,
      integer      /* ldRhs */
    ) const = 0;

    void
    load_bottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ
    );

    // block0 = row0 * col0
    // blockN = rowN * colN
    void
    load_top_bottom(
      real_type const block0[], integer ld0,
      real_type const blockN[], integer ldN
    );

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

    void select_last_block_solver_LU()   { m_la_factorization = &m_la_lu;   }
    void select_last_block_solver_LU_Q() { m_la_factorization = &m_la_lupq; }
    void select_last_block_solver_QR()   { m_la_factorization = &m_la_qr;   }
    void select_last_block_solver_QRP()  { m_la_factorization = &m_la_qrp;  }
    void select_last_block_solver_SVD()  { m_la_factorization = &m_la_svd;  }
    void select_last_block_solver_LSS()  { m_la_factorization = &m_la_lss;  }
    void select_last_block_solver_LSY()  { m_la_factorization = &m_la_lsy;  }
    void select_last_block_solver_PINV() { m_la_factorization = &m_la_pinv; }

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

    void select_last_border_block_solver_LU()   { m_bb_factorization = &m_bb_lu;   }
    void select_last_border_block_solver_LU_Q() { m_bb_factorization = &m_bb_lupq; }
    void select_last_border_block_solver_QR()   { m_bb_factorization = &m_bb_qr;   }
    void select_last_border_block_solver_QRP()  { m_bb_factorization = &m_bb_qrp;  }
    void select_last_border_block_solver_SVD()  { m_bb_factorization = &m_bb_svd;  }
    void select_last_border_block_solver_LSS()  { m_bb_factorization = &m_bb_lss;  }
    void select_last_border_block_solver_LSY()  { m_bb_factorization = &m_bb_lsy;  }
    void select_last_border_block_solver_PINV() { m_bb_factorization = &m_bb_pinv; }

    void
    last_block_factorize();

    void
    factorize_bordered();

    void
    solve_bordered( real_type [] ) const;

    void
    solve_bordered(
      integer      /* nrhs  */,
      real_type [] /* rhs   */,
      integer      /* ldRhs */
    ) const;

    // All in one
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

    // All in one
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

    // aux function
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

    void
    dump_ccoord( ostream_type & stream ) const;

    void
    dump_to_Maple( ostream_type & stream ) const;

    /*\
     |   ___ _ __   __ _ _ __ ___  ___
     |  / __| '_ \ / _` | '__/ __|/ _ \
     |  \__ \ |_) | (_| | |  \__ \  __/
     |  |___/ .__/ \__,_|_|  |___/\___|
     |      |_|
    \*/

    integer
    sparse_nnz() const;

    void
    sparse_pattern( integer I[], integer J[] ) const;

    void
    sparse_values( real_type vals[] ) const;

  };
}

///
/// eof: BlockBidiagonal.hxx
///
