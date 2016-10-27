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

#ifndef LU_BABD_QR_HH
#define LU_BABD_QR_HH

#include "CyclicReductionQR.hh"
#include <vector>

namespace alglin {

  using namespace std ;

  //! available LU factorization code
  typedef enum {
    BABD_QR_LASTBLOCK_LU  = 0,
    BABD_QR_LASTBLOCK_QR  = 1,
    BABD_QR_LASTBLOCK_QRP = 2,
    BABD_QR_LASTBLOCK_SVD = 3
  } BABD_QR_LASTBLOCK_Choice;

  /*
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
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
   *           University of Trento<br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
   */
  template <typename QR_type>
  class BabdQR {
  private:

    typedef typename QR_type::valueType         valueType ;
    typedef typename QR_type::valueType*        valuePointer ;
    typedef typename QR_type::valueType const * valueConstPointer ;

    BabdQR(BabdQR const &) ;
    BabdQR const & operator = (BabdQR const &) ;

    /*
    //
    //  Matrix structure
    //
    //                  n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n                        n     q
    //  +-----+-----+-----+----.........-----+-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |   
    //  |                                                :     |
    //  :                                                :      > n * nblock
    //  :                                                :     | 
    //  :                                                :     |
    //  :                                                :     |
    //  :                              +-----+-----+-----+     |
    //  :                              | Au  |  0  |  0  |     |
    //  :                        +-----+-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  |  0  | n   |
    //  +-----+-----+---......---+-----+-----+=====+=====+    /
    //  |     |     |            |     |     !     |     !
    //  | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
    //  |     |     |            |     |     !     |     !
    //  +-----+-----+---......---+-----+-----+=====+=====+
    //
    */

    ///////////////////////////////////////////////////////

    mutable vector<valueType>  tmpV ;
    CyclicReductionQR<QR_type> CR ;

    LU<valueType>  la_lu ;
    QR<valueType>  la_qr ;
    QRP<valueType> la_qrp ;
    SVD<valueType> la_svd ;

    integer n ;
    integer m ;
    
    Factorization<valueType> * factorization ;

    BABD_QR_LASTBLOCK_Choice last_block ;

  public:

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    explicit
    BabdQR( integer nth = integer(std::thread::hardware_concurrency()) )
    : CR(nth)
    {}
    #else
    explicit BabdQR() : CR() {}
    #endif

    ~BabdQR() {}

    //! load matrix in the class
    /*!
      \param nblock number of (square) blocks
      \param n      size of the blocks
      \param q      extra bc
      \param AdAu   pointer to the blocks diagonal ad upper diagonal
      \param H0     pointer to the block \f$ H_0 \f$
      \param HN     pointer to the block \f$ H_N \f$
      \param Hq     pointer to the block \f$ H_q \f$
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
    factorize( BABD_QR_LASTBLOCK_Choice choice,
               // ----------------------------
               integer           nblock,
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
