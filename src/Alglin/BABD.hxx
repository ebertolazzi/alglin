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
/// file: BABD.hxx
///

namespace alglin {

  using namespace ::std;

  //!
  //!  LU decomposition of a ABD matrix
  //!
  //!
  //!  \date     June 30, 2007
  //!  \version  1.0
  //!  \note     first release June 30, 2007
  //!
  //!  \author   Enrico Bertolazzi
  //!
  //!  \par      Affiliation:
  //!            Department of Industrial Engineering<br>
  //!            University of Trento <br>
  //!            Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
  //!            enrico.bertolazzi@\unitn.it
  //!
  //!
  template <typename t_Value>
  class BABD {
  public:
    //!
    //! available LU factorization code
    //!
    using BABD_Choice = enum class BABD_Choice : integer {
      DIAZ                 = 1, // no CR
      CYCLIC_REDUCTION_LU  = 2, // LU+QR
      CYCLIC_REDUCTION_QR  = 3, // CR+QR
      CYCLIC_REDUCTION_QRP = 4  // CR+QR
    };

    static
    string
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

    typedef t_Value real_type;

    BABD( BABD<t_Value> const & ) = delete;
    BABD<t_Value> const & operator = ( BABD<t_Value> const & ) = delete;

    DiazLU<t_Value>            m_diaz_LU;
    BorderedCR<t_Value>        m_bordered;

    BlockBidiagonal<t_Value> * m_babd_solver{nullptr};

  public:

    explicit constexpr
    BABD() : m_babd_solver(&m_bordered) {}

    ~BABD() {}

    // filling bidiagonal part of the matrix
    void
    load_blocks( real_type const AdAu[], integer ldA )
    { m_babd_solver->load_blocks( AdAu, ldA ); }

    void
    load_block( integer nbl, real_type const AdAu[], integer ldA )
    { m_babd_solver->load_block( nbl, AdAu, ldA ); }

    void
    load_block_left( integer nbl, real_type const Ad[], integer ldA )
    { m_babd_solver->load_block_left( nbl, Ad, ldA ); }

    void
    load_block_right( integer nbl, real_type const Au[], integer ldA )
    { m_babd_solver->load_block_right( nbl, Au, ldA ); }

    // Border Bottom blocks
    void
    set_zero_bottom_blocks()
    { m_babd_solver->set_zero_bottom_blocks(); }

    void
    load_bottom_blocks( real_type const C[], integer ldC )
    { m_babd_solver->load_bottom_blocks( C, ldC ); }

    void
    load_bottom_block( integer nbl, real_type const C[], integer ldC )
    { m_babd_solver->load_bottom_block( nbl, C, ldC ); }

    void
    add_to_bottom_block( integer nbl, real_type const C[], integer ldC )
    { m_babd_solver->add_to_bottom_block( nbl, C, ldC ); }

    void
    add_to_bottom_block2( integer nbl, real_type const C[], integer ldC )
    { m_babd_solver->add_to_bottom_block2( nbl, C, ldC ); }

    void
    load_bottom_last_block( real_type const C[], integer ldC )
    { m_babd_solver->load_bottom_last_block( C, ldC ); }

    // Border Right blocks
    void
    set_zero_right_blocks()
    { m_babd_solver->set_zero_right_blocks(); }

    void
    load_right_blocks( real_type const B[], integer ldB )
    { m_babd_solver->load_right_blocks( B, ldB ); }

    void
    loadRightBlock( integer nbl, real_type const B[], integer ldB )
    { m_babd_solver->loadRightBlock( nbl, B, ldB ); }

    void
    loadRightLastBlock( real_type const B[], integer ldB )
    { m_babd_solver->loadRightLastBlock( B, ldB ); }

    // Border RBblock
    void
    setZeroRBblock()
    { m_babd_solver->setZeroRBblock(); }

    void
    load_RB_block( real_type const D[], integer ldD )
    { m_babd_solver->load_RB_block( D, ldD ); }

    // Bottom BC
    void
    load_bottom(
      real_type const H0[], integer ld0,
      real_type const HN[], integer ldN,
      real_type const Hq[], integer ldQ
    ) {
      m_babd_solver->load_bottom( H0, ld0, HN, ldN, Hq, ldQ );
    }

    void
    load_top_bottom(
      integer         row0,
      integer         col0,
      real_type const block0[],
      integer         ld0,
      // ----------------------------
      integer         rowN,
      integer         colN,
      real_type const blockN[],
      integer         ldN
    ) {
      m_babd_solver->load_top_bottom(
        row0, col0, block0, ld0,
        rowN, colN, blockN, ldN
      );
    }

    void
    select_last_block_solver(
      typename BlockBidiagonal<t_Value>::BB_LASTBLOCK_Choice choice
    ) {
      m_babd_solver->select_last_block_solver( choice );
    }

    void
    select_last_border_block_solver(
      typename BlockBidiagonal<t_Value>::BB_LASTBLOCK_Choice choice
    ) {
      m_babd_solver->select_last_border_block_solver( choice );
    }

    void
    allocate( integer nblk, integer n, integer q, integer nb )
    { m_babd_solver->allocate( nblk, n, q, nb ); }

    void
    selectSolver( BABD_Choice choice ) {
      switch ( choice ) {
        case BABD_Choice::DIAZ:
          m_babd_solver = &m_diaz_LU;
          break;
        case BABD_Choice::CYCLIC_REDUCTION_LU:
        case BABD_Choice::CYCLIC_REDUCTION_QR:
        case BABD_Choice::CYCLIC_REDUCTION_QRP:
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
    load_BC( // ----------------------
      integer num_initial_BC,
      integer num_final_BC,
      integer num_cyclic_BC,
      // ----------------------
      integer num_initial_OMEGA,
      integer num_final_OMEGA,
      integer num_cyclic_OMEGA,
      // ----------------------
      real_type H0[], integer ld0,
      real_type HN[], integer ldN,
      real_type Hq[], integer ldq
    ) {
      m_babd_solver->load_BC(
        num_initial_BC,    num_final_BC,    num_cyclic_BC,
        num_initial_OMEGA, num_final_OMEGA, num_cyclic_OMEGA,
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
    solve( real_type in_out[] ) const
    { m_babd_solver->solve( in_out ); }

    void
    solve( integer nrhs, real_type in_out[], integer ldIO ) const
    { m_babd_solver->solve( nrhs, in_out, ldIO ); }

    void
    solve_bordered( real_type in_out[] ) const
    { m_babd_solver->solve_bordered( in_out ); }

    void
    solve_bordered(
      integer   nrhs,
      real_type rhs[],
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
    sparse_nnz() const
    { return m_babd_solver->sparse_nnz(); }

    void
    sparse_pattern( integer I[], integer J[] ) const
    { m_babd_solver->sparse_pattern(I,J); }

    void
    sparse_values( real_type vals[] ) const
    { m_babd_solver->sparse_values(vals); }

  };
}

///
/// eof: BABD.hxx
///
