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
  class BabdQR : public CyclicReductionQR<QR_type> {
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

    Malloc<valueType> baseValue ;

    mutable valuePointer tmpV ;
    valuePointer H0Nq ;

    integer m ;
    
  public:

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    explicit
    BabdQR( integer nth = integer(std::thread::hardware_concurrency()) )
    : CyclicReductionQR<QR_type>(nth)
    , baseValue("BabdQR_values")
    {}
    #else
    explicit BabdQR() : CyclicReductionQR<QR_type>() {}
    #endif

    ~BabdQR() {}

    //! load matrix in the class
    /*!
      \param q  extra bc
      \param H0 pointer to the block \f$ H_0 \f$
      \param HN pointer to the block \f$ H_N \f$
      \param Hq pointer to the block \f$ H_q \f$
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
    virtual
    void
    loadBottom( integer           q,
                valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hq, integer ldQ ) ;

    virtual
    void
    loadBC( integer /*numInitialBc*/,
            integer /*numFinalBc*/,
            integer /*numCyclicBC*/,
            // ----------------------
            integer numInitialETA,
            integer numFinalETA,
            integer numCyclicOMEGA,
            // ----------------------
            valueConstPointer H0, integer ld0,
            valueConstPointer HN, integer ldN,
            valueConstPointer Hq, integer ldQ ) {
      this->loadBottom( numInitialETA + numFinalETA + numCyclicOMEGA,
                        H0, ld0, HN, ldN, Hq, ldQ ) ;
    }

    virtual
    void
    loadTopBottom( // ----------------------------
                   integer           row0,
                   integer           col0,
                   valueConstPointer block0,
                   integer           ld0,
                   // ----------------------------
                   integer           rowN,
                   integer           colN,
                   valueConstPointer blockN,
                   integer           ldN ) ;

    virtual
    void
    factorize() ;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( valuePointer in_out ) const ;
    
    void
    factorize( LASTBLOCK_Choice  choice,
               integer           nblk,
               integer           n,
               integer           q,
               valueConstPointer AdAu,
               valueConstPointer H0,
               valueConstPointer HN,
               valueConstPointer Hq ) {
      this->allocate( nblk, n );
      this->loadBlocks( AdAu, n ) ;
      loadBottom( q, H0, n+q, HN, n+q, Hq, n+q ) ;
      this->selectLastBlockSolver( choice ) ;
      factorize() ;
    }

  } ;
}

#endif
