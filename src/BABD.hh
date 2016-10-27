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
    BABD_AUTOMATIC  = 0,
    // ----------------
    BABD_COLROW_LU  = 1,
    BABD_COLROW_QR  = 2,
    BABD_COLROW_QRP = 3,
    BABD_COLROW_SVD = 4,
    // ----------------
    BABD_ARCECO     = 5,
    // ----------------
    BABD_AMODIO_LU  = 6,
    BABD_AMODIO_QR  = 7,
    BABD_AMODIO_QRP = 8,
    BABD_AMODIO_SVD = 9,
    // ----------------
    BABD_BLOCK_LU             = 10,
    BABD_CYCLIC_REDUCTION_QR  = 11,
    BABD_CYCLIC_REDUCTION_QRP = 12
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

    // allocate temporary
    std::vector<valueType> block0 ;
    std::vector<valueType> blockN ;

    DiazLU<t_Value>   diaz_LU ;
    ArcecoLU<t_Value> arceco_LU ;

    // -------------------------
    AmodioLU<t_Value>     amodio_LU ;
    BlockLU<t_Value>      block_LU ;
    BabdQR<QR<t_Value> >  babd_QR ;
    BabdQR<QRP<t_Value> > babd_QRP ;

    integer row0, col0 ;
    integer rowN, colN ;
    integer numEquations ;

    integer numInitialBc ;
    integer numFinalBc ;
    integer numCyclicBC ;
    integer numInitialETA ;
    integer numFinalETA ;
    integer numCyclicOMEGA ;

    BABD_Choice solver_used ;

    void shift( valuePointer in_out ) const ;
    void unshift( valuePointer in_out ) const ;

  public:

    explicit BABD()
    : solver_used(BABD_AMODIO_LU)
    {}
    
    ~BABD() {}

    /*
    //    __            _             _
    //   / _| __ _  ___| |_ ___  _ __(_)_______
    //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
    //  |  _| (_| | (__| || (_) | |  | |/ /  __/
    //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
    */
    void
    factorize( BABD_Choice  solver,
               // ----------------------
               integer      numInitialBc,
               integer      numFinalBc,
               integer      numCyclicBC,
               // ----------------------
               integer      numInitialETA,
               integer      numFinalETA,
               integer      numCyclicOMEGA,
               // ----------------------
               integer      numBlock,
               valuePointer AdAu,
               valuePointer H0,
               valuePointer HN,
               valuePointer Hq ) ;
    /*             _
    //   ___  ___ | |_   _____
    //  / __|/ _ \| \ \ / / _ \
    //  \__ \ (_) | |\ V /  __/
    //  |___/\___/|_| \_/ \___|
    */
    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) const ;

  } ;
}

#endif
