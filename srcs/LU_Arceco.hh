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

#ifndef MECHATRONIX_CORE_ARCECO_LU_HH
#define MECHATRONIX_CORE_ARCECO_LU_HH

#include "Alglin.hh"
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
  typedef double  valueType ;
  typedef int     indexType ;
  typedef double* valuePointer ;
  typedef int*    indexPointer ;
  
  class ArcecoLU {
  private:

    Malloc<valueType> baseValue ;
    Malloc<indexType> baseIndex ;

    ArcecoLU( ArcecoLU const & ) ;
    ArcecoLU const & operator = ( ArcecoLU const & ) ;

    indexType    nRow0, nCol0 ;
    indexType    nRowN, nColN ;
    indexType    numEquations ;
    indexType    numInitialETA ;

    valuePointer AR, X  ;
    indexPointer PIVOT  ;
    indexPointer MTR    ;
    indexType    NBLOCK ;
    
    Arceco arcecoSolver ;

  public:

    explicit ArcecoLU() ;
    ~ArcecoLU() ;

    //! load matrix in the class
    /*!
      \param numInitialBc   number of initial boundary condition
      \param numFinalBc     number of final boundary condition
      \param numInitialETA
      \param numFinalETA
      \param numBlock
      \param AdAu
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
    load( indexType    numInitialBc,
          indexType    numFinalBc,
          indexType    numInitialETA,
          indexType    numFinalETA,
          indexType    numBlock,
          valuePointer AdAu,
          valuePointer H0,
          valuePointer HN,
          valuePointer Hq ) ;

    //! factorize the matrix
    void factorize() ;

    //! solve linear sistem using internal factorized matrix
    void solve( valuePointer in_out ) ;

  } ;
}

#endif
