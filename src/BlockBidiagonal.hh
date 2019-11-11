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

#ifndef BLOCK_BIDIAGONAL_HH
#define BLOCK_BIDIAGONAL_HH

#include "Alglin_Config.hh"

#ifdef __GNUC__ 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

namespace alglin {

  //! available LU factorization code
  typedef enum {
    LASTBLOCK_LU  = 0,
    LASTBLOCK_QR  = 1,
    LASTBLOCK_QRP = 2,
    LASTBLOCK_SVD = 3
  } LASTBLOCK_Choice;

  extern std::string LastBlock_to_string( LASTBLOCK_Choice c );

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

    typedef t_Value valueType;

  private:

    BlockBidiagonal(BlockBidiagonal const &);
    BlockBidiagonal const & operator = (BlockBidiagonal const &);

  protected:

    Malloc<valueType> baseValue;
    Malloc<integer>   baseInteger;

    integer nblock; //!< total number of blocks
    integer n;      //!< size of square blocks
    integer q;      //!< extra BC
    integer nb;     //!< border size

    // some derived constants
    integer neq;
    integer nx2;
    integer nxn;
    integer nxnx2;
    integer nxnb;

    integer numInitialBC;
    integer numFinalBC;
    integer numCyclicBC;
    integer numInitialOMEGA;
    integer numFinalOMEGA;
    integer numCyclicOMEGA;

    Matrix<t_Value> la_matrix;
    Matrix<t_Value> bb_matrix;
    LinearSystemSolver<t_Value> * la_factorization;
    LinearSystemSolver<t_Value> * bb_factorization;

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

    valueType * DE_blk;
    valueType * H0Nq;
    valueType * block0;
    valueType * blockN;
    valueType * Bmat;
    valueType * Cmat;
    valueType * Dmat;

  private:

    LU<t_Value>  la_lu,  bb_lu;
    QR<t_Value>  la_qr,  bb_qr;
    QRP<t_Value> la_qrp, bb_qrp;
    SVD<t_Value> la_svd, bb_svd;

  public:

    explicit
    BlockBidiagonal()
    : baseValue("BlockBidiagonal_values")
    , baseInteger("BlockBidiagonal_integers")
    , nblock(0)
    , n(0)
    , q(0)
    , nb(0)
    , neq(0)
    , nx2(0)
    , nxn(0)
    , nxnx2(0)
    , nxnb(0)
    , numInitialBC(0)
    , numFinalBC(0)
    , numCyclicBC(0)
    , numInitialOMEGA(0)
    , numFinalOMEGA(0)
    , numCyclicOMEGA(0)
    , la_factorization(&la_lu)
    , bb_factorization(&bb_lu)
    , DE_blk(nullptr)
    , H0Nq(nullptr)
    , block0(nullptr)
    , blockN(nullptr)
    , Bmat(nullptr)
    , Cmat(nullptr)
    , Dmat(nullptr)
    {}

    virtual ~BlockBidiagonal()
    {}

    //! allocatew and resize the problem
    void
    allocate(
      integer _nblock,
      integer _n,
      integer _nb,
      // ----------------------
      integer _numInitialBC,
      integer _numFinalBC,
      integer _numCyclicBC,
      // ----------------------
      integer _numInitialOMEGA,
      integer _numFinalOMEGA,
      integer _numCyclicOMEGA,
      // ----------------------
      integer num_extra_r,
      integer num_extra_i
    );

    void
    allocateTopBottom(
      integer _nblock,
      integer _n,
      integer _row0,
      integer _col0,
      integer _rowN,
      integer _colN,
      integer _nb,
      integer num_extra_r,
      integer num_extra_i
    ) {
      allocate(
        _nblock, _n, _nb,
        _row0, _rowN, 0,
        _col0-_n, _colN-_n, 0,
        num_extra_r, num_extra_i
      );
    }

    // filling bidiagonal part of the matrix
    void
    loadBlocks( valueType const AdAu[], integer ldA )
    { gecopy( n, nblock * nx2, AdAu, ldA, DE_blk, n ); }

    void
    loadBlock( integer nbl, valueType const AdAu[], integer ldA )
    { gecopy( n, nx2, AdAu, ldA, DE_blk + nbl*nxnx2, n ); }

    void
    loadBlockLeft( integer nbl, valueType const Ad[], integer ldA )
    { gecopy( n, n, Ad, ldA, DE_blk + nbl*nxnx2, n ); }

    void
    loadBlockRight( integer nbl, valueType const Au[], integer ldA )
    { gecopy( n, n, Au, ldA, DE_blk + nbl*nxnx2 + nxn, n ); }

    // Border Bottom blocks
    void
    setZeroBottomBlocks()
    { alglin::zero( nb*neq, Cmat, 1 ); }

    void
    loadBottomBlocks( valueType const C[], integer ldC )
    { gecopy( nb, neq, C, ldC, Cmat, nb ); }

    void
    loadBottomBlock( integer nbl, valueType const C[], integer ldC ) {
      LW_ASSERT(
        ldC >= nb, "loadBottomBlock( {}, C, ldC = {} ) bad ldC", nbl, ldC
      );
      valueType * CC = Cmat + nbl*nxnb;
      gecopy( nb, n, C, ldC, CC, nb );
    }

    void
    addtoBottomBlock( integer nbl, valueType const C[], integer ldC ) {
      LW_ASSERT(
        ldC >= nb, "addtoBottomBlock( {}, C, ldC = {} )bad ldC", nbl, ldC
      );
      valueType * CC = Cmat + nbl*nxnb;
      geadd( nb, n, 1.0, C, ldC, 1.0, CC, nb, CC, nb );
    }

    // add to bottom block nbl and nbl+1
    void
    addtoBottomBlock2( integer nbl, valueType const C[], integer ldC ) {
      LW_ASSERT(
        ldC >= nb, "addtoBottomBlock2( {}, C, ldC = {} ) bad ldC", nbl, ldC
      );
      valueType * CC = Cmat + nbl*nxnb;
      geadd( nb, nx2, 1.0, C, ldC, 1.0, CC, nb, CC, nb );
    }

    void
    loadBottomLastBlock( valueType const C[], integer ldC )
    { gecopy( nb, q, C, ldC, Cmat + (nblock+1)*nxnb, nb ); }

    // Border Right blocks
    void
    setZeroRightBlocks()
    { alglin::zero( neq*nb, Bmat, 1 ); }

    void
    loadRightBlocks( valueType const B[], integer ldB )
    { gecopy( neq, nb, B, ldB, Bmat, neq ); }

    void
    loadRightBlock( integer nbl, valueType const B[], integer ldB )
    { alglin::gecopy( n, nb, B, ldB, Bmat + nbl*n, neq ); }

    void
    loadRightLastBlock( valueType const B[], integer ldB )
    { alglin::gecopy( n+q, nb, B, ldB, Bmat + neq-n-q, neq ); }

    // Border RBblock
    void
    setZeroRBblock()
    { alglin::zero( nb*nb, Dmat, 1 ); }

    void
    loadRBblock( valueType const D[], integer ldD )
    { alglin::gecopy( nb, nb, D, ldD, Dmat, nb ); }

    // final blocks after cyclic reduction
    valueType const *
    getPointer_LR() const
    { return DE_blk; }

    void
    getBlock_LR( valueType LR[], integer ldA ) const
    { gecopy( n, nx2, DE_blk, n, LR, ldA ); }

    void
    getBlock_L( valueType L[], integer ldA ) const
    { gecopy( n, n, DE_blk, n, L, ldA ); }

    void
    getBlock_R( valueType R[], integer ldA ) const
    { gecopy( n, n, DE_blk+nxn, n, R, ldA ); }

    // BC blocks
    void
    getBlock_H0( valueType H0[], integer ld0 ) const
    { gecopy( n+q, n, H0Nq, n+q, H0, ld0 ); }

    void
    getBlock_HN( valueType HN[], integer ldN ) const
    { gecopy( n+q, n, H0Nq+n*(n+q), n+q, HN, ldN ); }

    void
    getBlock_Hq( valueType Hq[], integer ldQ ) const
    { gecopy( n+q, q, H0Nq+nx2*(n+q), n+q, Hq, ldQ ); }

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
    ) ALGLIN_PURE_VIRTUAL;

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
    ) ALGLIN_PURE_VIRTUAL;

    virtual
    void
    factorize() ALGLIN_PURE_VIRTUAL;

    virtual
    void
    solve( valueType [] ) const ALGLIN_PURE_VIRTUAL;

    virtual
    void
    solve(
      integer      /* nrhs  */,
      valueType [] /* rhs   */,
      integer      /* ldRhs */
    ) const ALGLIN_PURE_VIRTUAL;

    void
    loadBottom(
      valueType const H0[], integer ld0,
      valueType const HN[], integer ldN,
      valueType const Hq[], integer ldQ
    );

    // block0 = row0 * col0
    // blockN = rowN * colN
    void
    loadTopBottom(
      valueType const block0[], integer ld0,
      valueType const blockN[], integer ldN
    );

    void
    selectLastBlockSolver( LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case LASTBLOCK_LU:  la_factorization = &la_lu;  break;
        case LASTBLOCK_QR:  la_factorization = &la_qr;  break;
        case LASTBLOCK_QRP: la_factorization = &la_qrp; break;
        case LASTBLOCK_SVD: la_factorization = &la_svd; break;
      }
    }

    void selectLastBlockSolver_LU()  { la_factorization = &la_lu; }
    void selectLastBlockSolver_QR()  { la_factorization = &la_qr; }
    void selectLastBlockSolver_QRP() { la_factorization = &la_qrp; }
    void selectLastBlockSolver_SVD() { la_factorization = &la_svd; }

    void
    selectLastBorderBlockSolver( LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case LASTBLOCK_LU:  bb_factorization = &bb_lu;  break;
        case LASTBLOCK_QR:  bb_factorization = &bb_qr;  break;
        case LASTBLOCK_QRP: bb_factorization = &bb_qrp; break;
        case LASTBLOCK_SVD: bb_factorization = &bb_svd; break;
      }
    }

    void selectLastBorderBlockSolver_LU()  { la_factorization = &la_lu;  }
    void selectLastBorderBlockSolver_QR()  { la_factorization = &la_qr;  }
    void selectLastBorderBlockSolver_QRP() { la_factorization = &la_qrp; }
    void selectLastBorderBlockSolver_SVD() { la_factorization = &la_svd; }

    void
    last_block_factorize();

    void
    factorize_bordered();

    void
    solve_bordered( valueType [] ) const;

    void
    solve_bordered(
      integer      /* nrhs  */,
      valueType [] /* rhs   */,
      integer      /* ldRhs */
    ) const;

    // All in one
    void
    factorize(
      valueType const AdAu[],
      valueType const B[],
      valueType const C[],
      valueType const D[],
      valueType const H0[],
      valueType const HN[],
      valueType const Hq[]
    ) {
      this->loadBlocks( AdAu, n );
      this->loadBottom( H0, n+q, HN, n+q, Hq, n+q );
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
      valueType const AdAu[],
      valueType const H0[],
      valueType const HN[],
      valueType const Hq[]
    ) {
      LW_ASSERT0( nb == 0, "factorize nb > 0 and no border assigned" );
      this->loadBlocks( AdAu, n );
      this->loadBottom( H0, n+q, HN, n+q, Hq, n+q );
      this->factorize();
    }

    // aux function
    void
    Mv( valueType const x[], valueType res[] ) const;

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
    sparseValues( valueType vals[] ) const;

  };
}

#endif
