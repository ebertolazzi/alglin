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

#ifndef LU_ARCECO_HH
#define LU_ARCECO_HH

#include "Alglin.hh"
#include "Alglin++.hh"
#include "LU_ArcecoSolver.hh"

namespace alglin {

  using namespace ::std ;
  
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
   *           `enrico.bertolazzi@unitn.it`
   * 
   */
  template <typename t_Value>
  class ArcecoLU {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;

    ArcecoLU( ArcecoLU<t_Value> const & ) ;
    ArcecoLU<t_Value> const & operator = ( ArcecoLU<t_Value> const & ) ;

    integer      nRow0, nCol0 ;
    integer      nRowN, nColN ;
    integer      numEquations ;
    integer      numInitialETA ;

    valuePointer AR, X  ;
    integer *    PIVOT  ;
    integer *    MTR    ;
    integer      NBLOCK ;
    
    Arceco<t_Value> arcecoSolver ;

  public:

    explicit ArcecoLU() ;
    ~ArcecoLU() ;

    //! load matrix in the class
    /*!
      \param numInitialBc   number of initial boundary condition
      \param numFinalBc     number of final boundary condition
      \param numInitialETA  initial bc blocks
      \param numFinalETA    final bc blocks
      \param numBlock       number of diagonal blocks
      \param AdAu           diagonal blocks
      \param H0             pointer to the block \f$ H_0 \f$
      \param HN             pointer to the block \f$ H_N \f$
      \param Hq             pointer to the block \f$ H_q \f$

      \code
      
      nq = numInitialBc + numFinalBc
      q  = numInitialETA + numFinalETA
      n  = nq - q
      
      compatibility conditions
      numFinalETA + numInitialETA == numInitialBc + numFinalBc - n

      Matrix structure
    
       +                     +
       |  Ad Au              |
       |     Ad Au           |
       |          ...        |
       |            Ad Au    |
       |  H0           HN Hq |
       +                     +
    
              nq x n         nq x n         nq x q
            +-------+       +-------+       +---+
            |  0    |       |  NZ   |       |y|0|
       H0 = +-------+  HN = +-------+  Hq = +---+
            |  NZ   |       |  0    |       |0|x|
            |  NZ   |       |       |       |0|x|
            +-------+       +-------+       +---+
    
       +-+
       |y| numFinalBc x numFinalETA
       +-+

       +-+
       |x| numInitialBc x numInitialETA
       |x|
       +-+
      \endcode
    */

    void
    load( integer      numInitialBc,
          integer      numFinalBc,
          integer      numInitialETA,
          integer      numFinalETA,
          integer      numBlock,
          valuePointer AdAu,
          valuePointer H0,
          valuePointer HN,
          valuePointer Hq ) ;

    /*
    //    __            _             _
    //   / _| __ _  ___| |_ ___  _ __(_)_______
    //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
    //  |  _| (_| | (__| || (_) | |  | |/ /  __/
    //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
    */
    //! factorize the matrix
    void
    factorize()
    { arcecoSolver.factorize() ; }

    /*             _
    //   ___  ___ | |_   _____
    //  / __|/ _ \| \ \ / / _ \
    //  \__ \ (_) | |\ V /  __/
    //  |___/\___/|_| \_/ \___|
    */
    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) {
      alglin::copy( numEquations - nRow0, in_out, 1, X + nRow0, 1 ) ;
      alglin::copy( nRow0, in_out + numEquations - nRow0, 1, X, 1 ) ;
      arcecoSolver.solve ( X ) ;
      alglin::copy( numEquations - numInitialETA, X + numInitialETA, 1, in_out, 1 ) ;
      alglin::copy( numInitialETA, X, 1, in_out + numEquations - numInitialETA, 1 ) ;
    }
  } ;
}

#endif
