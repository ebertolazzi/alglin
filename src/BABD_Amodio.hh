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

#ifndef BABD_AMODIO_HH
#define BABD_AMODIO_HH

#include "Alglin.hh"
#include "Alglin++.hh"
#include "CyclicReduction.hh"
#include <vector>

namespace alglin {

  //! available LU factorization code
  typedef enum {
    AMODIO_LASTBLOCK_LU  = 0,
    AMODIO_LASTBLOCK_QR  = 1,
    AMODIO_LASTBLOCK_QRP = 2,
    AMODIO_LASTBLOCK_SVD = 3
  } AMODIO_LASTBLOCK_Choice;

  /*
  //      _                        _ _       _    _   _ 
  //     / \   _ __ ___   ___   __| (_) ___ | |  | | | |
  //    / _ \ | '_ ` _ \ / _ \ / _` | |/ _ \| |  | | | |
  //   / ___ \| | | | | | (_) | (_| | | (_) | |__| |_| |
  //  /_/   \_\_| |_| |_|\___/ \__,_|_|\___/|_____\___/ 
  */

  //! LU decomposition of a BABD matrix
  /*!
   * 
   * \date     May 30, 2006
   * \version  1.0
   * \note     first release May 30, 2006
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
  class AmodioLU {
  private:
  
    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    AmodioLU(AmodioLU const &) ;
    AmodioLU const & operator = (AmodioLU const &) ;

    mutable vector<t_Value>  tmpV ;
    CyclicReduction<t_Value> CR ;
    LU<t_Value>  la_lu ;
    QR<t_Value>  la_qr ;
    QRP<t_Value> la_qrp ;
    SVD<t_Value> la_svd ;

    integer n ;
    integer m ;
    
    Factorization<t_Value> * factorization ;

    AMODIO_LASTBLOCK_Choice last_block ;

  public:

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    explicit
    AmodioLU( integer nth = integer(std::thread::hardware_concurrency()) )
    : CR(nth)
    {}
    #else
    explicit AmodioLU() : CR() {}
    #endif

    ~AmodioLU() {}

    //! load matrix in the class
    /*!
      \param nblk number of (square) blocks
      \param n    size of the blocks
      \param q    extra bc
      \param AdAu pointer to the blocks diagonal ad upper diagonal
      \param H0   pointer to the block \f$ H_0 \f$
      \param HN   pointer to the block \f$ H_N \f$
      \param Hq   pointer to the block \f$ H_q \f$
      \code
      Matrix structure
      
                      n * nblock
        ________________^_________________
       /                                  \
         n     n     n                        n     q
      +-----+-----+-----+----.........-----+-----+-----+    \
      |  Ad | Au  |  0  |                  |  0  |  0  | n   |
      +-----+-----+-----+             -----+-----+-----+     |
      |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |
      |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |   
      |                                                :     |
      :                                                :      > n * nblock
      :                                                :     | 
      :                                                :     |
      :                                                :     |
      :                              +-----+-----+-----+     |
      :                              | Au  |  0  |  0  |     |
      :                        +-----+-----+-----+-----+     |
      :                        |  0  | Ad  | Au  |  0  | n   |
      +-----+-----+---......---+-----+-----+=====+=====+    /
      |     |     |            |     |     !     |     !
      | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
      |     |     |            |     |     !     |     !
      +-----+-----+---......---+-----+-----+=====+=====+
      \endcode
    */
    void
    factorize( AMODIO_LASTBLOCK_Choice choice,
               // ----------------------------
               integer           nblk,
               integer           n,
               integer           q,
               valueConstPointer AdAu,
               valueConstPointer H0,
               valueConstPointer HN,
               valueConstPointer Hq ) ;

    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) const ;

  } ;
}

#endif
