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

#ifndef BLOCK_BIDIAGONAL_HH
#define BLOCK_BIDIAGONAL_HH

#include "Alglin.hh"
#include "Alglin++.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace alglin {

  //! available LU factorization code
  typedef enum {
    LASTBLOCK_LU  = 0,
    LASTBLOCK_QR  = 1,
    LASTBLOCK_QRP = 2,
    LASTBLOCK_SVD = 3
  } LASTBLOCK_Choice;

  //! available LU factorization code
  typedef enum {
    CYCLIC_REDUCTION_LU  = 0,
    CYCLIC_REDUCTION_QR  = 1,
    CYCLIC_REDUCTION_QRP = 2
  } CYCLIC_REDUCTION_Choice;

  /*\
   |  ____  _            _      ____  _     _ _                               _
   | | __ )| | ___   ___| | __ | __ )(_) __| (_) __ _  __ _  ___  _ __   __ _| |
   | |  _ \| |/ _ \ / __| |/ / |  _ \| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
   | | |_) | | (_) | (__|   <  | |_) | | (_| | | (_| | (_| | (_) | | | | (_| | |
   | |____/|_|\___/ \___|_|\_\ |____/|_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
   |                                                  |___/
  \*/

  //! Cyclic reduction of a block bidiagonal matrix
  /*!
   * 
   * \date     October 25, 2016
   * \version  1.0
   * \note     October 25, 2016
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
  class BlockBidiagonal {
  private:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    Malloc<valueType> baseValue ;

    BlockBidiagonal(BlockBidiagonal const &) ;
    BlockBidiagonal const & operator = (BlockBidiagonal const &) ;

  protected:
    integer nblock ; //!< total number of blocks
    integer n      ; //!< size of square blocks

    /*
    //
    //  Matrix structure
    //
    //                 (n+1) * nblock
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n   
    //  +-----+-----+-----+----.........-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+     |   
    //  |                                                |
    //  :                                                 > n * nblock
    //  :                                                | 
    //  :                                                |
    //  :                                                |
    //  :                              +-----+-----+     |
    //  :                              | Au  |  0  |     |
    //  :                        +-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  | n   |
    //  +-----+-----+---......---+-----+-----+=====+    /
    //
    */

    ///////////////////////////////////////////////////////
    valuePointer AdAu_blk ;

    // some derived constants
    integer nx2 ;
    integer nxn ;
    integer nxnx2 ;

  public:

    explicit
    BlockBidiagonal()
    : baseValue("BlockBidiagonal_values")
    , nblock(0)
    , n(0)
    , nx2(0)
    , nxn(0)
    , nxnx2(0)
    {}

    virtual ~BlockBidiagonal()
    {}

    //! load matrix in the class
    virtual
    void
    allocate( integer _nblock, integer _n ) {
      if ( _nblock != nblock || n != _n ) {
        nblock = _nblock ;
        n      = _n ;
        nx2    = 2*n ;
        nxn    = n*n ;
        nxnx2  = nxn*2 ;
        baseValue.allocate(size_t(nblock*nxnx2)) ;
        AdAu_blk = baseValue(size_t(nblock*nxnx2)) ;
      }
    }

    integer getNblock() const { return nblock ; }
    integer getN() const { return n ; }

    void
    loadBlocks( valueConstPointer AdAu, integer ldA ) {
      gecopy( n, nblock * nx2, AdAu, ldA, AdAu_blk, n ) ;
    }

    void
    loadBlock( integer           nbl,
               valueConstPointer AdAu,
               integer           ldA ) {
      gecopy( n, nx2, AdAu, ldA, AdAu_blk + nbl*nxnx2, n ) ;
    }

    void
    loadBlockLeft( integer           nbl,
                   valueConstPointer Ad,
                   integer           ldA ) {
      gecopy( n, n, Ad, ldA, AdAu_blk + nbl*nxnx2, n ) ;
    }
    
    void
    loadBlockRight( integer           nbl,
                    valueConstPointer Au,
                    integer           ldA ) {
      gecopy( n, n, Au, ldA, AdAu_blk + nbl*nxnx2 + nxn, n ) ;
    }

    valueConstPointer
    getPointer_LR() const
    { return AdAu_blk ; }

    void
    getBlock_LR( valuePointer LR, integer ldA ) const
    { gecopy( n, nx2, AdAu_blk, n, LR, ldA ) ; }

    void
    getBlock_L( valuePointer L, integer ldA ) const
    { gecopy( n, n, AdAu_blk, n, L, ldA ) ; }
    
    void
    getBlock_R( valuePointer R, integer ldA ) const
    { gecopy( n, n, AdAu_blk+nxn, n, R, ldA ) ; }

    //! load matrix in the class
    void
    reduce( integer           _nblk,
            integer           _n,
            valueConstPointer AdAu,
            integer           ldAdDu ) {
      allocate( _nblk, _n ) ;
      gecopy( n, nblock*nx2, AdAu, ldAdDu, AdAu_blk, n ) ;
      reduce() ;
    }

    /*
    //
    //  Apply reduction to the coeff matric of the linear system
    //
    //                 (n+1) * nblock
    //    ___________________^____________________
    //   /                                        \
    //     n     n     n                        n   
    //  +-----+-----+-----+----.........-----+-----+    +----+     +----+ \
    //  |  Ad | Au  |  0  |                  |  0  | n  | y1 |     | b1 | |
    //  +-----+-----+-----+             -----+-----+    +----+     +----+ |
    //  |  0  | Ad  | Au  |  0               |  0  | n  | y2 |     | b2 | |
    //  +-----+-----+-----+-----+       -----+-----+    +----+     +----+ |
    //  |  0  |  0  | Ad  | Au  |            |  0  | n  | y3 |     | b3 | |
    //  +-----+-----+-----+-----+       -----+-----+    +----+     +----+ |
    //  |                                                                 |
    //  :                                                      =          > n * nblock
    //  :                                                                 |
    //  :                                                                 |
    //  :                                                                 |
    //  :                              +-----+-----+    +----+     +----+ |
    //  :                              | Au  |  0  |    | .. |     | .. | |
    //  :                        +-----+-----+-----+    +----+     +----+ |
    //  :                        |  0  | Ad  | Au  | n  | yn |     | bn | |
    //  +-----+-----+---......---+-----+-----+=====+    +----+     +----+ /
    //
    //  After reduction the linear system reduce to
    //
    //     n     n
    //  +-----+-----+  +----+
    //  |  L  | R   |  | y1 |
    //  +-----+-----+  +----+
    //                 | yn |
    //                 +----+
    */
    virtual void reduce()
    { ALGLIN_ERROR("BlockBidiagonal::reduce() not defined!") ; }

    /*
    //  Apply reduction to the RHS of the linear system.
    //  After reduction the linear system becomes:
    //
    //     n     n
    //  +-----+-----+  +----+   +----+
    //  |  L  | R   |  | y1 |   | c1 |
    //  +-----+-----+  +----+ = +----+
    //                 | yn |   | cn |
    //                 +----+   +----+
    */
    virtual void forward( valuePointer ) const
    { ALGLIN_ERROR("BlockBidiagonal::forward( valuePointer ) not defined!") ; }

    /*
    //  Given y1 and yn of the reduced linear system compute y2, y3, ... y(n-1)
    */
    virtual void backward( valuePointer ) const
    { ALGLIN_ERROR("BlockBidiagonal::backward( valuePointer ) not defined!") ; }

    virtual void factorize()
    { ALGLIN_ERROR("BlockBidiagonal::factorize( valuePointer ) not defined!") ; }

    virtual void solve( valuePointer ) const
    { ALGLIN_ERROR("BlockBidiagonal::solve( valuePointer ) not defined!") ; }

  } ;
}

#endif
