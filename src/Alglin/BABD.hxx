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
/// file: BABD.hxx
///

namespace alglin {

  using namespace ::std;

  //! available LU factorization code
  typedef enum {
    BABD_DIAZ                 = 1, // no CR
    BABD_CYCLIC_REDUCTION_LU  = 2, // LU+QR
    BABD_CYCLIC_REDUCTION_QR  = 3, // CR+QR
    BABD_CYCLIC_REDUCTION_QRP = 4  // CR+QR
  } BABD_Choice;

  extern string BABD_Choice_to_string( BABD_Choice c );

  //! LU decomposition of a ABD matrix
  /*!
   *
   * \date     June 30, 2007
   * \version  1.0
   * \note     first release June 30, 2007
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi@\unitn.it
   *
  \*/
  template <typename t_Value>
  class BABD {
  private:

    typedef t_Value valueType;

    BABD( BABD<t_Value> const & );
    BABD<t_Value> const & operator = ( BABD<t_Value> const & );

    DiazLU<t_Value>            m_diaz_LU;
    BorderedCR<t_Value>        m_bordered;

    BlockBidiagonal<t_Value> * m_babd_solver;

  public:

    explicit UTILS_CONSTEXPR BABD()
    : m_babd_solver(&m_bordered)
    {}

    ~BABD() {}

    // filling bidiagonal part of the matrix
    void
    loadBlocks( valueType const AdAu[], integer ldA )
    { m_babd_solver->loadBlocks( AdAu, ldA ); }

    void
    loadBlock( integer nbl, valueType const AdAu[], integer ldA )
    { m_babd_solver->loadBlock( nbl, AdAu, ldA ); }

    void
    loadBlockLeft( integer nbl, valueType const Ad[], integer ldA )
    { m_babd_solver->loadBlockLeft( nbl, Ad, ldA ); }

    void
    loadBlockRight( integer nbl, valueType const Au[], integer ldA )
    { m_babd_solver->loadBlockRight( nbl, Au, ldA ); }

    // Border Bottom blocks
    void
    setZeroBottomBlocks()
    { m_babd_solver->setZeroBottomBlocks(); }

    void
    loadBottomBlocks( valueType const C[], integer ldC )
    { m_babd_solver->loadBottomBlocks( C, ldC ); }

    void
    loadBottomBlock( integer nbl, valueType const C[], integer ldC )
    { m_babd_solver->loadBottomBlock( nbl, C, ldC ); }

    void
    addtoBottomBlock( integer nbl, valueType const C[], integer ldC )
    { m_babd_solver->addtoBottomBlock( nbl, C, ldC ); }

    void
    addtoBottomBlock2( integer nbl, valueType const C[], integer ldC )
    { m_babd_solver->addtoBottomBlock2( nbl, C, ldC ); }

    void
    loadBottomLastBlock( valueType const C[], integer ldC )
    { m_babd_solver->loadBottomLastBlock( C, ldC ); }

    // Border Right blocks
    void
    setZeroRightBlocks()
    { m_babd_solver->setZeroRightBlocks(); }

    void
    loadRightBlocks( valueType const B[], integer ldB )
    { m_babd_solver->loadRightBlocks( B, ldB ); }

    void
    loadRightBlock( integer nbl, valueType const B[], integer ldB )
    { m_babd_solver->loadRightBlock( nbl, B, ldB ); }

    void
    loadRightLastBlock( valueType const B[], integer ldB )
    { m_babd_solver->loadRightLastBlock( B, ldB ); }

    // Border RBblock
    void
    setZeroRBblock()
    { m_babd_solver->setZeroRBblock(); }

    void
    loadRBblock( valueType const D[], integer ldD )
    { m_babd_solver->loadRBblock( D, ldD ); }

    // Bottom BC
    void
    loadBottom(
      valueType const H0[], integer ld0,
      valueType const HN[], integer ldN,
      valueType const Hq[], integer ldQ
    ) {
      m_babd_solver->loadBottom( H0, ld0, HN, ldN, Hq, ldQ );
    }

    void
    loadTopBottom(
      integer         row0,
      integer         col0,
      valueType const block0[],
      integer         ld0,
      // ----------------------------
      integer         rowN,
      integer         colN,
      valueType const blockN[],
      integer         ldN
    ) {
      m_babd_solver->loadTopBottom(
        row0, col0, block0, ld0,
        rowN, colN, blockN, ldN
      );
    }

    void
    selectLastBlockSolver( LASTBLOCK_Choice choice )
    { m_babd_solver->selectLastBlockSolver( choice ); }

    void
    selectLastBorderBlockSolver( LASTBLOCK_Choice choice )
    { m_babd_solver->selectLastBorderBlockSolver( choice ); }

    void
    allocate( integer nblk, integer n, integer q, integer nb )
    { m_babd_solver->allocate( nblk, n, q, nb ); }

    void
    selectSolver( BABD_Choice choice ) {
      switch ( choice ) {
        case BABD_DIAZ:
          m_babd_solver = &m_diaz_LU;
          break;
        case BABD_CYCLIC_REDUCTION_LU:
        case BABD_CYCLIC_REDUCTION_QR:
        case BABD_CYCLIC_REDUCTION_QRP:
          m_babd_solver = &m_bordered;
          break;
      };
    }

    /*\
     |   _                 _ ____   ____
     |  | | ___   __ _  __| | __ ) / ___|
     |  | |/ _ \ / _` |/ _` |  _ \| |
     |  | | (_) | (_| | (_| | |_) | |___
     |  |_|\___/ \__,_|\__,_|____/ \____|
    \*/
    void
    loadBC( // ----------------------
      integer numInitialBC,
      integer numFinalBC,
      integer numCyclicBC,
      // ----------------------
      integer numInitialOMEGA,
      integer numFinalOMEGA,
      integer numCyclicOMEGA,
      // ----------------------
      valueType H0[], integer ld0,
      valueType HN[], integer ldN,
      valueType Hq[], integer ldq
    ) {
      m_babd_solver->loadBC(
        numInitialBC,    numFinalBC,    numCyclicBC,
        numInitialOMEGA, numFinalOMEGA, numCyclicOMEGA,
        H0, ld0, HN, ldN, Hq, ldq
      );
    }

    /*\
     |    __            _             _
     |   / _| __ _  ___| |_ ___  _ __(_)_______
     |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
     |  |  _| (_| | (__| || (_) | |  | |/ /  __/
     |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
    \*/
    void
    factorize()
    { m_babd_solver->factorize(); }

    void
    factorize_bordered()
    { m_babd_solver->factorize_bordered(); }

    /*\
     |             _
     |   ___  ___ | |_   _____
     |  / __|/ _ \| \ \ / / _ \
     |  \__ \ (_) | |\ V /  __/
     |  |___/\___/|_| \_/ \___|
    \*/
    //! solve linear sistem using internal factorized matrix
    void
    solve( valueType in_out[] ) const
    { m_babd_solver->solve( in_out ); }

    void
    solve( integer nrhs, valueType in_out[], integer ldIO ) const
    { m_babd_solver->solve( nrhs, in_out, ldIO ); }

    void
    solve_bordered( valueType in_out[] ) const
    { m_babd_solver->solve_bordered( in_out ); }

    void
    solve_bordered(
      integer   nrhs,
      valueType rhs[],
      integer   ldRhs
    ) const {
      m_babd_solver->solve_bordered( nrhs, rhs, ldRhs );
    }

    /*\
     |   ____
     |  |  _ \ _   _ _ __ ___  _ __
     |  | | | | | | | '_ ` _ \| '_ \
     |  | |_| | |_| | | | | | | |_) |
     |  |____/ \__,_|_| |_| |_| .__/
     |                        |_|
    \*/

    void
    dump_to_Maple( ostream_type & stream ) const
    { m_babd_solver->dump_to_Maple( stream ); }

    void
    dump_ccoord( ostream_type & stream ) const
    { m_babd_solver->dump_ccoord( stream ); }

    /*\
     |   ___ _ __   __ _ _ __ ___  ___
     |  / __| '_ \ / _` | '__/ __|/ _ \
     |  \__ \ |_) | (_| | |  \__ \  __/
     |  |___/ .__/ \__,_|_|  |___/\___|
     |      |_|
    \*/

    integer
    sparseNnz() const
    { return m_babd_solver->sparseNnz(); }

    void
    sparsePattern( integer I[], integer J[] ) const
    { m_babd_solver->sparsePattern(I,J); }

    void
    sparseValues( valueType vals[] ) const
    { m_babd_solver->sparseValues(vals); }

  };
}

///
/// eof: BABD.hxx
///
