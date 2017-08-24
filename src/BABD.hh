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

#ifndef BABD_HH
#define BABD_HH

#include "Alglin.hh"
#include "Alglin++.hh"

#include "ABD_Diaz.hh"
#include "ABD_Arceco.hh"
#include "BABD_BorderedCR.hh"
#include "BABD_Block.hh"

namespace alglin {

  using namespace ::std ;

  //! available LU factorization code
  typedef enum {
    BABD_DIAZ                 = 1, // no CR
    BABD_CYCLIC_REDUCTION_LU  = 2, // LU+QR
    BABD_CYCLIC_REDUCTION_QR  = 3, // CR+QR
    BABD_CYCLIC_REDUCTION_QRP = 4  // CR+QR
  } BABD_Choice;

  extern string BABD_Choice_to_string( BABD_Choice c ) ;

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
   */
  template <typename t_Value>
  class BABD {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    BABD( BABD<t_Value> const & ) ;
    BABD<t_Value> const & operator = ( BABD<t_Value> const & ) ;

    DiazLU<t_Value>       diaz_LU ;
    BorderedCR<t_Value>   bordered ;

    BlockBidiagonal<t_Value> * babd_solver ;

  public:

    explicit ALGLIN_CONSTEXPR BABD()
    : babd_solver(&bordered)
    {}

    ~BABD() {}

    // filling bidiagonal part of the matrix
    void
    loadBlocks( valueConstPointer AdAu, integer ldA )
    { babd_solver->loadBlocks( AdAu, ldA ) ; }

    void
    loadBlock( integer nbl, valueConstPointer AdAu, integer ldA )
    { babd_solver->loadBlock( nbl, AdAu, ldA ) ; }

    void
    loadBlockLeft( integer nbl, valueConstPointer Ad, integer ldA )
    { babd_solver->loadBlockLeft( nbl, Ad, ldA ) ; }

    void
    loadBlockRight( integer nbl, valueConstPointer Au, integer ldA )
    { babd_solver->loadBlockRight( nbl, Au, ldA ) ; }

    // Border Bottom blocks
    void
    setZeroBottomBlocks()
    {  babd_solver->setZeroBottomBlocks() ; }

    void
    loadBottomBlocks( valueConstPointer C, integer ldC )
    {  babd_solver->loadBottomBlocks( C, ldC ) ; }

    void
    loadBottomBlock( integer nbl, valueConstPointer C, integer ldC )
    { babd_solver->loadBottomBlock( nbl, C, ldC ) ; }

    void
    addtoBottomBlock( integer nbl, valueConstPointer C, integer ldC )
    { babd_solver->addtoBottomBlock( nbl, C, ldC ) ; }

    void
    addtoBottomBlock2( integer nbl, valueConstPointer C, integer ldC )
    { babd_solver->addtoBottomBlock2( nbl, C, ldC ) ; }

    void
    loadBottomLastBlock( valueConstPointer C, integer ldC )
    { babd_solver->loadBottomLastBlock( C, ldC ) ; }

    // Border Right blocks
    void
    setZeroRightBlocks()
    { babd_solver->setZeroRightBlocks() ; }

    void
    loadRightBlocks( valueConstPointer B, integer ldB )
    { babd_solver->loadRightBlocks( B, ldB ) ; }

    void
    loadRightBlock( integer nbl, valueConstPointer B, integer ldB )
    { babd_solver->loadRightBlock( nbl, B, ldB ) ; }

    void
    loadRightLastBlock( valueConstPointer B, integer ldB )
    { babd_solver->loadRightLastBlock( B, ldB ) ; }

    // Border RBblock
    void
    setZeroRBblock()
    { babd_solver->setZeroRBblock() ; }

    void
    loadRBblock( valueConstPointer D, integer ldD )
    { babd_solver->loadRBblock( D, ldD ) ; }

    // Bottom BC
    void
    loadBottom( valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hq, integer ldQ )
    { babd_solver->loadBottom( H0, ld0, HN, ldN, Hq, ldQ ) ; }

    void
    loadTopBottom( integer           _row0,
                   integer           _col0,
                   valueConstPointer _block0,
                   integer           _ld0,
                   // ----------------------------
                   integer           _rowN,
                   integer           _colN,
                   valueConstPointer _blockN,
                   integer           _ldN )
    { babd_solver->loadTopBottom( _row0, _col0, _block0, _ld0,
                                  _rowN, _colN, _blockN, _ldN ) ; }

    void
    selectLastBlockSolver( LASTBLOCK_Choice choice )
    { babd_solver->selectLastBlockSolver( choice ) ; }

    void
    selectLastBorderBlockSolver( LASTBLOCK_Choice choice )
    { babd_solver->selectLastBorderBlockSolver( choice ) ; }

    void
    allocate( integer nblk, integer n, integer q, integer nb )
    { babd_solver->allocate( nblk, n, q, nb ) ; }

    void
    selectSolver( BABD_Choice choice ) {
      switch ( choice ) {
        case BABD_DIAZ:
          babd_solver = &diaz_LU ;
          break ;
        case BABD_CYCLIC_REDUCTION_LU:
        case BABD_CYCLIC_REDUCTION_QR:
        case BABD_CYCLIC_REDUCTION_QRP:
          babd_solver = &bordered ;
          break ;
      } ;
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
            integer      numInitialBC,
            integer      numFinalBC,
            integer      numCyclicBC,
            // ----------------------
            integer      numInitialOMEGA,
            integer      numFinalOMEGA,
            integer      numCyclicOMEGA,
            // ----------------------
            valuePointer H0, integer ld0,
            valuePointer HN, integer ldN,
            valuePointer Hq, integer ldq ) {
      babd_solver->loadBC( numInitialBC,  numFinalBC,  numCyclicBC,
                           numInitialOMEGA, numFinalOMEGA, numCyclicOMEGA,
                           H0, ld0, HN, ldN, Hq, ldq ) ;
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
    { babd_solver->factorize() ; }

    void
    factorize_bordered()
    { babd_solver->factorize_bordered() ; }

    /*\
     |             _
     |   ___  ___ | |_   _____
     |  / __|/ _ \| \ \ / / _ \
     |  \__ \ (_) | |\ V /  __/
     |  |___/\___/|_| \_/ \___|
    \*/
    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) const
    { babd_solver->solve( in_out ) ; }

    void
    solve( integer nrhs, valuePointer in_out, integer ldIO ) const
    { babd_solver->solve( nrhs, in_out, ldIO ) ; }

    void
    solve_bordered( valuePointer in_out ) const
    { babd_solver->solve_bordered( in_out ) ; }

    void
    solve_bordered( integer      nrhs,
                    valuePointer rhs,
                    integer      ldRhs ) const
    { babd_solver->solve_bordered( nrhs, rhs, ldRhs ) ; }

    /*\
     |   ____
     |  |  _ \ _   _ _ __ ___  _ __
     |  | | | | | | | '_ ` _ \| '_ \
     |  | |_| | |_| | | | | | | |_) |
     |  |____/ \__,_|_| |_| |_| .__/
     |                        |_|
    \*/
    void
    dump_to_Maple( basic_ostream<char> & stream ) const
    { babd_solver->dump_to_Maple( stream ) ; }

    void
    dump_ccoord( basic_ostream<char> & stream ) const
    { babd_solver->dump_ccoord( stream ) ; }

  } ;
}

#endif
