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
   */

  template <typename t_Value>
  class BlockBidiagonal {
  public:

    typedef t_Value real_type;

    //! available LU factorization code
    typedef enum {
      BB_LASTBLOCK_LU   = 0,
      BB_LASTBLOCK_LUPQ = 1,
      BB_LASTBLOCK_QR   = 2,
      BB_LASTBLOCK_QRP  = 3,
      BB_LASTBLOCK_SVD  = 4,
      BB_LASTBLOCK_LSS  = 5,
      BB_LASTBLOCK_LSY  = 6,
      BB_LASTBLOCK_PINV = 7
    } BB_LASTBLOCK_Choice;

    static
    std::string
    LastBlock_to_string( BB_LASTBLOCK_Choice c ) {
      switch ( c ) {
        case BB_LASTBLOCK_LU:   return "last block LU";
        case BB_LASTBLOCK_LUPQ: return "last block LUPQ";
        case BB_LASTBLOCK_QR:   return "last block QR";
        case BB_LASTBLOCK_QRP:  return "last block QRP";
        case BB_LASTBLOCK_SVD:  return "last block SVD";
        case BB_LASTBLOCK_LSS:  return "last block LSS";
        case BB_LASTBLOCK_LSY:  return "last block LSY";
        case BB_LASTBLOCK_PINV: return "last block PINV";
      }
      return "last block not selected";
    }

  private:

    BlockBidiagonal(BlockBidiagonal const &) = delete;
    BlockBidiagonal const & operator = (BlockBidiagonal const &) = delete;

  protected:

    Malloc<real_type> m_baseValue;
    Malloc<integer>   m_baseInteger;

    integer m_number_of_blocks; //!< total number of blocks
    integer m_block_size;       //!< size of square blocks
    integer m_extra_bc;         //!< extra BC
    integer m_border_size;      //!< border size
    integer m_num_equations;

    // some derived constants
    integer n_x_n;
    integer n_x_nb;

    integer m_numInitialBC;
    integer m_numFinalBC;
    integer m_numCyclicBC;
    integer m_numInitialOMEGA;
    integer m_numFinalOMEGA;
    integer m_numCyclicOMEGA;

    Matrix<t_Value>               m_la_matrix;
    Matrix<t_Value>               m_bb_matrix;
    LinearSystemSolver<t_Value> * m_la_factorization;
    LinearSystemSolver<t_Value> * m_bb_factorization;

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

    real_type * m_DE_blk;
    real_type * m_H0Nq;
    real_type * m_block0;
    real_type * m_blockN;
    real_type * m_Bmat;
    real_type * m_Cmat;
    real_type * m_Dmat;

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
    : m_baseValue("BlockBidiagonal_values")
    , m_baseInteger("BlockBidiagonal_integers")
    , m_number_of_blocks(0)
    , m_block_size(0)
    , m_extra_bc(0)
    , m_border_size(0)
    , m_num_equations(0)
    , n_x_n(0)
    , n_x_nb(0)
    , m_numInitialBC(0)
    , m_numFinalBC(0)
    , m_numCyclicBC(0)
    , m_numInitialOMEGA(0)
    , m_numFinalOMEGA(0)
    , m_numCyclicOMEGA(0)
    , m_la_factorization(&m_la_lu)
    , m_bb_factorization(&m_bb_lu)
    , m_DE_blk(nullptr)
    , m_H0Nq(nullptr)
    , m_block0(nullptr)
    , m_blockN(nullptr)
    , m_Bmat(nullptr)
    , m_Cmat(nullptr)
    , m_Dmat(nullptr)
    {}

    virtual ~BlockBidiagonal()
    {}

    //! allocatew and resize the problem
    void
    allocate(
      integer nblock,
      integer n,
      integer nb,
      // ----------------------
      integer numInitialBC,
      integer numFinalBC,
      integer numCyclicBC,
      // ----------------------
      integer numInitialOMEGA,
      integer numFinalOMEGA,
      integer numCyclicOMEGA,
      // ----------------------
      integer num_extra_r,
      integer num_extra_i
    );

    void
    allocateTopBottom(
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
    void
    loadBlocks( real_type const AdAu[], integer ldA ) {
      integer const & nblock = m_number_of_blocks;
      integer const & n      = m_block_size;
      gecopy( n, 2*n*nblock, AdAu, ldA, m_DE_blk, n );
    }

    void
    loadBlock( integer nbl, real_type const AdAu[], integer ldA ) {
      integer const & n = m_block_size;
      gecopy( n, 2*n, AdAu, ldA, m_DE_blk + 2*nbl*n_x_n, n );
    }

    void
    loadBlockLeft( integer nbl, real_type const Ad[], integer ldA ) {
      integer const & n = m_block_size;
      gecopy( n, n, Ad, ldA, m_DE_blk + 2*nbl*n_x_n, n );
    }

    void
    loadBlockRight( integer nbl, real_type const Au[], integer ldA ) {
      integer const & n = m_block_size;
      gecopy( n, n, Au, ldA, m_DE_blk + (2*nbl+1)*n_x_n, n );
    }

    // Border Bottom blocks
    void
    setZeroBottomBlocks() {
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      alglin::zero( nb*neq, m_Cmat, 1 );
    }

    void
    loadBottomBlocks( real_type const C[], integer ldC ) {
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      gecopy( nb, neq, C, ldC, m_Cmat, nb );
    }

    void
    loadBottomBlock( integer nbl, real_type const C[], integer ldC ) {
      integer const & n  = m_block_size;
      integer const & nb = m_border_size;
      UTILS_ASSERT(
        ldC >= nb, "loadBottomBlock( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
      );
      real_type * CC = m_Cmat + nbl*n_x_nb;
      gecopy( nb, n, C, ldC, CC, nb );
    }

    void
    addtoBottomBlock( integer nbl, real_type const C[], integer ldC ) {
      integer const & n  = m_block_size;
      integer const & nb = m_border_size;
      UTILS_ASSERT(
        ldC >= nb, "addtoBottomBlock( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
      );
      real_type * CC = m_Cmat + nbl*n_x_nb;
      geadd( nb, n, 1.0, C, ldC, 1.0, CC, nb, CC, nb );
    }

    // add to bottom block nbl and nbl+1
    void
    addtoBottomBlock2( integer nbl, real_type const C[], integer ldC ) {
      integer const & n  = m_block_size;
      integer const & nb = m_border_size;
      UTILS_ASSERT(
        ldC >= nb, "addtoBottomBlock2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC
      );
      real_type * CC = m_Cmat + nbl*n_x_nb;
      geadd( nb, 2*n, 1.0, C, ldC, 1.0, CC, nb, CC, nb );
    }

    void
    loadBottomLastBlock( real_type const C[], integer ldC ) {
      integer const & nblock = m_number_of_blocks;
      integer const & nb     = m_border_size;
      integer const & q      = m_extra_bc;
      gecopy( nb, q, C, ldC, m_Cmat + (nblock+1)*n_x_nb, nb );
    }

    // Border Right blocks
    void
    setZeroRightBlocks() {
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      alglin::zero( neq*nb, m_Bmat, 1 );
    }

    void
    loadRightBlocks( real_type const B[], integer ldB ) {
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      gecopy( neq, nb, B, ldB, m_Bmat, neq );
    }

    void
    loadRightBlock( integer nbl, real_type const B[], integer ldB ) {
      integer const & n   = m_block_size;
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      alglin::gecopy( n, nb, B, ldB, m_Bmat + nbl*n, neq );
    }

    void
    loadRightLastBlock( real_type const B[], integer ldB ) {
      integer const & n   = m_block_size;
      integer const & q   = m_extra_bc;
      integer const & nb  = m_border_size;
      integer const & neq = m_num_equations;
      alglin::gecopy( n+q, nb, B, ldB, m_Bmat + neq-n-q, neq );
    }

    // Border RBblock
    void
    setZeroRBblock() {
      integer const & nb = m_border_size;
      alglin::zero( nb*nb, m_Dmat, 1 );
    }

    void
    loadRBblock( real_type const D[], integer ldD ) {
      integer const & nb = m_border_size;
      alglin::gecopy( nb, nb, D, ldD, m_Dmat, nb );
    }

    // final blocks after cyclic reduction
    real_type const *
    getPointer_LR() const
    { return m_DE_blk; }

    void
    getBlock_LR( real_type LR[], integer ldA ) const {
      integer const & n = m_block_size;
      gecopy( n, 2*n, m_DE_blk, n, LR, ldA );
    }

    void
    getBlock_L( real_type L[], integer ldA ) const {
      integer const & n = m_block_size;
      gecopy( n, n, m_DE_blk, n, L, ldA );
    }

    void
    getBlock_R( real_type R[], integer ldA ) const {
      integer const & n = m_block_size;
      gecopy( n, n, m_DE_blk+n_x_n, n, R, ldA );
    }

    // BC blocks
    void
    getBlock_H0( real_type H0[], integer ld0 ) const {
      integer const & n = m_block_size;
      integer const & q = m_extra_bc;
      gecopy( n+q, n, m_H0Nq, n+q, H0, ld0 );
    }

    void
    getBlock_HN( real_type HN[], integer ldN ) const {
      integer const & n = m_block_size;
      integer const & q = m_extra_bc;
      integer nq = n + q;
      gecopy( nq, n, m_H0Nq+n*nq, nq, HN, ldN );
    }

    void
    getBlock_Hq( real_type Hq[], integer ldQ ) const {
      integer const & n = m_block_size;
      integer const & q = m_extra_bc;
      integer nq = n + q;
      gecopy( nq, q, m_H0Nq+2*n*nq, nq, Hq, ldQ );
    }

    virtual
    void
    allocate(
      integer /* nblock */,
      integer /* n      */,
      integer /* nb     */,
      // ----------------------
      integer /* numInitialBC */,
      integer /* numFinalBC   */,
      integer /* numCyclicBC  */,
      // ----------------------
      integer /* numInitialOMEGA */,
      integer /* numFinalOMEGA   */,
      integer /* numCyclicOMEGA  */
    ) = 0;

    virtual
    void
    allocateTopBottom(
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
    loadBottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ
    );

    // block0 = row0 * col0
    // blockN = rowN * colN
    void
    loadTopBottom(
      real_type const block0[], integer ld0,
      real_type const blockN[], integer ldN
    );

    void
    selectLastBlockSolver( BB_LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case BB_LASTBLOCK_LU:   m_la_factorization = &m_la_lu;   break;
        case BB_LASTBLOCK_LUPQ: m_la_factorization = &m_la_lupq; break;
        case BB_LASTBLOCK_QR:   m_la_factorization = &m_la_qr;   break;
        case BB_LASTBLOCK_QRP:  m_la_factorization = &m_la_qrp;  break;
        case BB_LASTBLOCK_SVD:  m_la_factorization = &m_la_svd;  break;
        case BB_LASTBLOCK_LSS:  m_la_factorization = &m_la_lss;  break;
        case BB_LASTBLOCK_LSY:  m_la_factorization = &m_la_lsy;  break;
        case BB_LASTBLOCK_PINV: m_la_factorization = &m_la_pinv; break;
      }
    }

    void selectLastBlockSolver_LU()   { m_la_factorization = &m_la_lu;   }
    void selectLastBlockSolver_LU_Q() { m_la_factorization = &m_la_lupq; }
    void selectLastBlockSolver_QR()   { m_la_factorization = &m_la_qr;   }
    void selectLastBlockSolver_QRP()  { m_la_factorization = &m_la_qrp;  }
    void selectLastBlockSolver_SVD()  { m_la_factorization = &m_la_svd;  }
    void selectLastBlockSolver_LSS()  { m_la_factorization = &m_la_lss;  }
    void selectLastBlockSolver_LSY()  { m_la_factorization = &m_la_lsy;  }
    void selectLastBlockSolver_PINV() { m_la_factorization = &m_la_pinv; }

    void
    selectLastBorderBlockSolver( BB_LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case BB_LASTBLOCK_LU:   m_bb_factorization = &m_bb_lu;   break;
        case BB_LASTBLOCK_LUPQ: m_bb_factorization = &m_bb_lupq; break;
        case BB_LASTBLOCK_QR:   m_bb_factorization = &m_bb_qr;   break;
        case BB_LASTBLOCK_QRP:  m_bb_factorization = &m_bb_qrp;  break;
        case BB_LASTBLOCK_SVD:  m_bb_factorization = &m_bb_svd;  break;
        case BB_LASTBLOCK_LSS:  m_bb_factorization = &m_bb_lss;  break;
        case BB_LASTBLOCK_LSY:  m_bb_factorization = &m_bb_lsy;  break;
        case BB_LASTBLOCK_PINV: m_bb_factorization = &m_bb_pinv; break;
      }
    }

    void selectLastBorderBlockSolver_LU()   { m_bb_factorization = &m_bb_lu;   }
    void selectLastBorderBlockSolver_LU_Q() { m_bb_factorization = &m_bb_lupq; }
    void selectLastBorderBlockSolver_QR()   { m_bb_factorization = &m_bb_qr;   }
    void selectLastBorderBlockSolver_QRP()  { m_bb_factorization = &m_bb_qrp;  }
    void selectLastBorderBlockSolver_SVD()  { m_bb_factorization = &m_bb_svd;  }
    void selectLastBorderBlockSolver_LSS()  { m_bb_factorization = &m_bb_lss;  }
    void selectLastBorderBlockSolver_LSY()  { m_bb_factorization = &m_bb_lsy;  }
    void selectLastBorderBlockSolver_PINV() { m_bb_factorization = &m_bb_pinv; }

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

      this->loadBlocks( AdAu, n );
      integer nq = n + q;
      this->loadBottom( H0, nq, HN, nq, Hq, nq );
      if ( nb > 0 ) {
        this->loadRightBlocks( B, neq );
        this->loadBottomBlocks( C, nb );
        this->loadRBblock( D, nb );
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
      this->loadBlocks( AdAu, n );
      this->loadBottom( H0, nq, HN, nq, Hq, nq );
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
    sparseNnz() const;

    void
    sparsePattern( integer I[], integer J[] ) const;

    void
    sparseValues( real_type vals[] ) const;

  };
}

///
/// eof: BlockBidiagonal.hxx
///
