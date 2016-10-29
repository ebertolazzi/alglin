/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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
#include "BABD_Amodio.hh"
#include "BABD_Block.hh"
#include "BABD_QR.hh"

namespace alglin {

  using namespace ::std ;

  //! available LU factorization code
  typedef enum {
    BABD_DIAZ                 = 1, // no CR
    BABD_AMODIO               = 2, // CR_LU
    BABD_CYCLIC_REDUCTION_QR  = 3, // CR+QR
    BABD_CYCLIC_REDUCTION_QRP = 4  // CR+QR
  } BABD_Choice;

  extern string BABD_Choice_to_string( BABD_Choice c ) ;

  /*
   * NB: prima le condizioni finali, poi quelle iniziali.
   *
   */

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
    AmodioLU<t_Value>     amodio_LU ;
    BabdQR<QR<t_Value> >  babd_QR ;
    BabdQR<QRP<t_Value> > babd_QRP ;
    
    BlockBidiagonal<t_Value> * babd_solver ;

  public:

    explicit BABD()
    : babd_solver(&amodio_LU)
    {}

    ~BABD() {}

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

    void
    loadBottom( integer           q,
                valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hq, integer ldQ )
    { babd_solver->loadBottom( q, H0, ld0, HN, ldN, Hq, ldQ ) ; }

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
    allocate( integer nblk, integer n )
    { babd_solver->allocate( nblk, n ) ; }

    void
    selectSolver( BABD_Choice choice ) {
      switch ( choice ) {
        case BABD_DIAZ:                 babd_solver = &diaz_LU   ; break ;
        case BABD_AMODIO:               babd_solver = &amodio_LU ; break ;
        case BABD_CYCLIC_REDUCTION_QR:  babd_solver = &babd_QR   ; break ;
        case BABD_CYCLIC_REDUCTION_QRP: babd_solver = &babd_QRP  ; break ;
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
            integer numInitialBc,
            integer numFinalBc,
            integer numCyclicBC,
            // ----------------------
            integer numInitialETA,
            integer numFinalETA,
            integer numCyclicOMEGA,
            // ----------------------
            valuePointer H0, integer ld0,
            valuePointer HN, integer ldN,
            valuePointer Hq, integer ldq ) {
      babd_solver->loadBC( numInitialBc,  numFinalBc,  numCyclicBC,
                           numInitialETA, numFinalETA, numCyclicOMEGA,
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

  } ;
}

#endif
