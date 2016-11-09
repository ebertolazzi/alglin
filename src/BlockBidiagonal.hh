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

  extern string LastBlock_to_string( LASTBLOCK_Choice c ) ;

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

    BlockBidiagonal(BlockBidiagonal const &) ;
    BlockBidiagonal const & operator = (BlockBidiagonal const &) ;

  protected:

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;

    integer nblock ; //!< total number of blocks
    integer n      ; //!< size of square blocks
    integer q      ; //!< extra BC
    integer nb     ; //!< border size

    // some derived constants
    integer neq ;
    integer nx2 ;
    integer nxn ;
    integer nxnx2 ;
    integer nxnb ;

    integer numInitialBc ;
    integer numFinalBc ;
    integer numCyclicBC ;
    integer numInitialOMEGA ;
    integer numFinalOMEGA ;
    integer numCyclicOMEGA ;

    Factorization<t_Value> * la_factorization ;

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
    //  +-----+-----+---......---+-----+-----+=====+--+  /
    //  |     |                              |     |  |  \
    //  | H0  |                              | HN  |Hq|  |
    //  |     |                              |     |  |  | n+q
    //  +-----+-----+---......---+-----+-----+=====+--+  /
    //                                               q
    //
    //  Bordered matrix
    //  / A  B \
    //  \ C  D /
    //
    */

    valuePointer AdAu_blk ;
    valuePointer H0Nq ;
    valuePointer block0 ;
    valuePointer blockN ;
    valuePointer Bmat, Cmat, Dmat ;

  private:

    LU<t_Value>  la_lu ;
    QR<t_Value>  la_qr ;
    QRP<t_Value> la_qrp ;
    SVD<t_Value> la_svd ;

  public:

    explicit
    BlockBidiagonal()
    : baseValue("BlockBidiagonal_values")
    , baseInteger("BlockBidiagonal_integers")
    , nblock(0)
    , n(0)
    , q(0)
    , nb(0)
    , neq(0)
    , nx2(0)
    , nxn(0)
    , nxnx2(0)
    , nxnb(0)
    , numInitialBc(0)
    , numFinalBc(0)
    , numCyclicBC(0)
    , numInitialOMEGA(0)
    , numFinalOMEGA(0)
    , numCyclicOMEGA(0)
    , la_factorization(&la_lu)
    {}

    virtual ~BlockBidiagonal()
    {}

    //! load matrix in the class
    void
    allocate( integer _nblock,
              integer _n,
              integer _q,
              integer _nb,
              integer num_extra_r,
              integer num_extra_i ) {
      nblock = _nblock ;
      n      = _n ;
      q      = _q ;
      nb     = _nb ;
      neq    = (nblock+1)*n+q ;
      nx2    = n*2 ;
      nxn    = n*n ;
      nxnx2  = nxn*2 ;
      nxnb   = n*nb ;
      integer AdAu_size = nblock*nxnx2;
      integer H0Nq_size = (n+q)*(nx2+q);
      integer BC_size   = nb*neq ;
      baseValue.allocate(size_t(AdAu_size+H0Nq_size+2*BC_size+nb*nb+num_extra_r)) ;
      baseInteger.allocate(size_t(num_extra_i)) ;
      AdAu_blk = baseValue(size_t(AdAu_size)) ;
      H0Nq     = baseValue(size_t(H0Nq_size)) ;
      Bmat     = baseValue(size_t(BC_size)) ;
      Cmat     = baseValue(size_t(BC_size)) ;
      Dmat     = baseValue(size_t(nb*nb)) ;
      block0   = nullptr ;
      blockN   = nullptr ;
    }

    // filling bidiagonal part of the matrix
    void
    loadBlocks( valueConstPointer AdAu, integer ldA )
    { gecopy( n, nblock * nx2, AdAu, ldA, AdAu_blk, n ) ; }

    void
    loadBlock( integer nbl, valueConstPointer AdAu, integer ldA )
    { gecopy( n, nx2, AdAu, ldA, AdAu_blk + nbl*nxnx2, n ) ; }

    void
    loadBlockLeft( integer nbl, valueConstPointer Ad, integer ldA )
    { gecopy( n, n, Ad, ldA, AdAu_blk + nbl*nxnx2, n ) ; }

    void
    loadBlockRight( integer nbl, valueConstPointer Au, integer ldA )
    { gecopy( n, n, Au, ldA, AdAu_blk + nbl*nxnx2 + nxn, n ) ; }

    // Border Bottom blocks
    void
    loadBottomBlocks( valueConstPointer C, integer ldC )
    { gecopy( nb, neq, C, ldC, Cmat, nb ) ; }

    void
    loadBottomBlock( integer nbl, valueConstPointer C, integer ldC )
    { gecopy( nb, n, C, ldC, Cmat + nbl*nxnb, nb ) ; }

    void
    loadBottomLastBlock( valueConstPointer C, integer ldC )
    { gecopy( nb, q, C, ldC, Cmat + (nblock+1)*nxnb, nb ) ; }

    // Border Right blocks
    void
    loadRightBlocks( valueConstPointer B, integer ldB )
    { gecopy( neq, nb, B, ldB, Bmat, neq ) ; }

    void
    loadRightBlock( integer nbl, valueConstPointer B, integer ldB )
    { gecopy( n, nb, B, ldB, Bmat + nbl*n, n ) ; }

    void
    loadRightLastBlock( valueConstPointer B, integer ldB )
    { gecopy( n+q, nb, B, ldB, Bmat + neq-n-q, n+q ) ; }

    // final blocks after cyclic reduction
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

    // BC blocks
    void
    getBlock_H0( valuePointer H0, integer ld0 ) const
    { gecopy( n+q, n, H0Nq, n+q, H0, ld0 ) ; }

    void
    getBlock_HN( valuePointer HN, integer ldN ) const
    { gecopy( n+q, n, H0Nq+n*(n+q), n+q, HN, ldN ) ; }

    void
    getBlock_Hq( valuePointer Hq, integer ldQ ) const
    { gecopy( n+q, q, H0Nq+nx2*(n+q), n+q, Hq, ldQ ) ; }

    virtual
    void
    allocate( integer /*nblock*/, integer /*n*/, integer /*q*/ )
    { ALGLIN_ERROR("BlockBidiagonal::allocate() not defined!") ; }

    virtual
    void
    factorize()
    { ALGLIN_ERROR("BlockBidiagonal::factorize() not defined!") ; }

    virtual
    void
    solve( valuePointer ) const
    { ALGLIN_ERROR("BlockBidiagonal::solve() not defined!") ; }

    virtual
    void
    solve( integer      /* nrhs*/,
           valuePointer /*rhs*/,
           integer      /*ldRhs*/ ) const
    { ALGLIN_ERROR("BlockBidiagonal::solve() not defined!") ; }

    void
    loadBottom( valueConstPointer H0, integer ld0,
                valueConstPointer HN, integer ldN,
                valueConstPointer Hq, integer ldQ ) ;

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

    void
    loadBC( integer numInitialBc,
            integer numFinalBc,
            integer numCyclicBC,
            // ----------------------
            integer numInitialOMEGA,
            integer numFinalOMEGA,
            integer numCyclicOMEGA,
            // ----------------------
            valueConstPointer H0, integer ld0,
            valueConstPointer HN, integer ldN,
            valueConstPointer Hq, integer ldQ ) ;

    void
    selectLastBlockSolver( LASTBLOCK_Choice choice ) {
      switch ( choice ) {
        case LASTBLOCK_LU:  la_factorization = &la_lu  ; break ;
        case LASTBLOCK_QR:  la_factorization = &la_qr  ; break ;
        case LASTBLOCK_QRP: la_factorization = &la_qrp ; break ;
        case LASTBLOCK_SVD: la_factorization = &la_svd ; break ;
      }
    }

    void
    last_block_factorize() ;

    void
    dumpMatrix ( basic_ostream<char> & stream ) const ;

  } ;
}

#endif
