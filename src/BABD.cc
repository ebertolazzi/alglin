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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "BABD.hh"

namespace alglin {

  string
  BABD_Choice_to_string( BABD_Choice c ) {
    string res ;
    switch ( c ) {
    case alglin::BABD_AUTOMATIC:  res = "Automatic"  ; break ;
    case alglin::BABD_COLROW_LU:  res = "Colrow+LU"  ; break ;
    case alglin::BABD_COLROW_QR:  res = "Colrow+QR"  ; break ;
    case alglin::BABD_COLROW_QRP: res = "Colrow+QRP" ; break ;
    case alglin::BABD_COLROW_SVD: res = "Colrow+SVD" ; break ;
    case alglin::BABD_ARCECO:     res = "Arceco"     ; break ;
    case alglin::BABD_AMODIO_LU:  res = "Amodio+LU"  ; break ;
    case alglin::BABD_AMODIO_QR:  res = "Amodio+QR"  ; break ;
    case alglin::BABD_AMODIO_QRP: res = "Amodio+QRP" ; break ;
    case alglin::BABD_AMODIO_SVD: res = "Amodio+SVD" ; break ;
    case alglin::BABD_BLOCK_LU:   res = "BlockLU"    ; break ;
    case alglin::BABD_CYCLIC_REDUCTION_QR:
      res = "babdQR" ;
      break ;
    case alglin::BABD_CYCLIC_REDUCTION_QRP:
      res = "babdQRP" ;
      break ;
    }
    return res ;
  }

  template <typename t_Value>
  void
  BABD<t_Value>::shift( valuePointer in_out ) const {
    std::rotate( in_out, in_out + numEquations - row0,in_out + numEquations ) ;
  }

  template <typename t_Value>
  void
  BABD<t_Value>::unshift( valuePointer in_out ) const {
    std::rotate( in_out, in_out + numInitialETA, in_out + numEquations ) ;
  }

  /*
  //    __            _             _
  //   / _| __ _  ___| |_ ___  _ __(_)_______
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */
  /*
  //  +                     +            nq x n         nq x n        nq x q
  //  |  Ad Au              |          +-------+       +-------+       +---+
  //  |     Ad Au           |          |  0    |       |  NZ   |       |y|0|
  //  |          ...        |     H0 = +-------+  HN = +-------+  Hq = +---+
  //  |            Ad Au    |          |  NZ   |       |  0    |       |0|x|
  //  |  H0           HN Hq |          |  NZ   |       |       |       |0|x|
  //  +                     +          +-------+       +-------+       +---+
  //
  //  Primo blocco per ARCECO               Ultimo blocco per ARCECO 
  //  (numInitialBc) x (n+numInitialETA)    (numFinalBc) x (n+numFinalETA)
  //  +-+-------+                           +-------+-+
  //  |x|  NZ   |                           |  NZ   |y|
  //  |x|  NZ   |                           +-------+-+
  //  +-+-------+
  //
  //  +-+                                   +-+
  //  |x| numInitialBc x numInitialETA      |y| numFinalBc x numFinalETA
  //  |x|                                   +-+
  //  +-+
  */
  template <typename t_Value>
  void
  BABD<t_Value>::factorize( BABD_Choice  _solver,
                            // ----------------------
                            integer      _numInitialBc,
                            integer      _numFinalBc,
                            integer      _numCyclicBC,
                            // ----------------------
                            integer      _numInitialETA,
                            integer      _numFinalETA,
                            integer      _numCyclicOMEGA,
                            // ----------------------
                            integer      numBlock,
                            valuePointer AdAu,
                            valuePointer H0,
                            valuePointer HN,
                            valuePointer Hq ) {

    solver_used    = _solver ;
    numInitialBc   = _numInitialBc ;
    numFinalBc     = _numFinalBc ;
    numCyclicBC    = _numCyclicBC ;
    numInitialETA  = _numInitialETA ;
    numFinalETA    = _numFinalETA ;
    numCyclicOMEGA = _numCyclicOMEGA ;

    integer nq = numInitialBc  + numFinalBc  + numCyclicBC ;
    integer q  = numInitialETA + numFinalETA + numCyclicOMEGA ;
    integer n  = nq - q ;
    
    if ( solver_used == BABD_AUTOMATIC ) solver_used = BABD_ARCECO ;

    switch ( solver_used ) {
    case BABD_COLROW_LU:
    case BABD_COLROW_QR:
    case BABD_COLROW_QRP:
    case BABD_COLROW_SVD:
    case BABD_ARCECO:
      if ( numCyclicOMEGA == 0 && numCyclicBC == 0 ) {
        row0         = numInitialBc ;
        col0         = n + numInitialETA ;
        rowN         = numFinalBc ;
        colN         = n + numFinalETA ;
        numEquations = numBlock * n + row0 + rowN ;

        // allocate temporary
        block0.resize(size_t(row0*col0)) ;
        blockN.resize(size_t(rowN*colN)) ;
    
        #define IDX0(I,J) ((I)+(J)*row0)
        #define IDXN(I,J) ((I)+(J)*rowN)
        std::fill( block0.begin(), block0.end(), 0 ) ;
        std::fill( blockN.begin(), blockN.end(), 0 ) ;
        //  +-+-------+
        //  |x|  NZ   |
        //  |x|  NZ   |
        //  +-+-------+
        gecopy( row0, numInitialETA, Hq+rowN+numFinalETA*nq, nq, &block0[IDX0(0,0)],             row0 ) ;
        gecopy( row0, n,             H0+rowN,                nq, &block0[IDX0(0,numInitialETA)], row0 ) ;
        //  +-------+-+
        //  |  NZ   |y|
        //  +-------+-+
        alglin::gecopy( rowN, n,           HN, nq, &blockN[IDXN(0,0)], rowN ) ;
        alglin::gecopy( rowN, numFinalETA, Hq, nq, &blockN[IDXN(0,n)], rowN ) ;
        #undef IDX0
        #undef IDXN
      } else {
        solver_used = BABD_AMODIO_LU ; // non si puo usare ABD
      }
      break;
    case BABD_AMODIO_LU:
    case BABD_AMODIO_QR:
    case BABD_AMODIO_QRP:
    case BABD_AMODIO_SVD:
    case BABD_BLOCK_LU:
    case BABD_CYCLIC_REDUCTION_QR:
    case BABD_CYCLIC_REDUCTION_QRP:
    case BABD_AUTOMATIC:
      break ;
    }

    switch ( solver_used ) {
    /*
    case BABD_COLROW_LU:
      colrow_LU.factorize( LASTBLOCK_LU,
                           row0, col0, &block0.front(),
                           numBlock, n, AdAu,
                           rowN, colN, &blockN.front() ) ;
      break ;
    case BABD_COLROW_QR:
      colrow_LU.factorize( LASTBLOCK_QR,
                           row0, col0, &block0.front(),
                           numBlock, n, AdAu,
                           rowN, colN, &blockN.front() ) ;
      break ;
    case BABD_COLROW_QRP:
      colrow_LU.factorize( LASTBLOCK_QRP,
                           row0, col0, &block0.front(),
                           numBlock, n, AdAu,
                           rowN, colN, &blockN.front() ) ;
      break ;
    case BABD_COLROW_SVD:
      colrow_LU.factorize( LASTBLOCK_SVD,
                           row0, col0, &block0.front(),
                           numBlock, n, AdAu,
                           rowN, colN, &blockN.front() ) ;
      break ;*/
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    case BABD_ARCECO:
      arceco_LU.factorize( row0, col0, &block0.front(),
                           numBlock, n, AdAu,
                           rowN, colN, &blockN.front() ) ;
      break ;
    case BABD_AMODIO_LU:
      amodio_LU.factorize( AMODIO_LASTBLOCK_LU, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_AMODIO_QR:
      amodio_LU.factorize( AMODIO_LASTBLOCK_QR, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_AMODIO_QRP:
      amodio_LU.factorize( AMODIO_LASTBLOCK_QRP, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_AMODIO_SVD:
      amodio_LU.factorize( AMODIO_LASTBLOCK_SVD, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_BLOCK_LU:
      block_LU.factorize( numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_CYCLIC_REDUCTION_QR:
      babd_QR.factorize( BABD_QR_LASTBLOCK_QR, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_CYCLIC_REDUCTION_QRP:
      babd_QRP.factorize( BABD_QR_LASTBLOCK_QRP, numBlock, n, q, AdAu, H0, HN, Hq ) ;
      break ;
    case BABD_AUTOMATIC:
      ALGLIN_ERROR("BABD<t_Value>::factorize -- no solver selected") ;
    }
  }

  /*             _
  //   ___  ___ | |_   _____
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */
  //! solve linear sistem using internal factorized matrix
  template <typename t_Value>
  void
  BABD<t_Value>::solve( valuePointer in_out ) const {
    switch ( solver_used ) {
    case BABD_COLROW_LU:
    case BABD_COLROW_QR:
    case BABD_COLROW_QRP:
    case BABD_COLROW_SVD:
      shift(in_out) ;
      diaz_LU.solve( in_out ) ;
      unshift(in_out) ;
      break ;
    case BABD_ARCECO:
      shift(in_out) ;
      arceco_LU.solve( in_out ) ;
      unshift(in_out) ;
      break ;
    case BABD_AMODIO_LU:
    case BABD_AMODIO_QR:
    case BABD_AMODIO_QRP:
    case BABD_AMODIO_SVD:
      amodio_LU.solve( in_out ) ;
      break ;
    case BABD_BLOCK_LU:
      block_LU.solve( in_out ) ;
      break ;
    case BABD_CYCLIC_REDUCTION_QR:
      babd_QR.solve( in_out ) ;
      break ;
    case BABD_CYCLIC_REDUCTION_QRP:
      babd_QRP.solve( in_out ) ;
      break ;
    case BABD_AUTOMATIC:
      ALGLIN_ERROR( "BABD<t_Value>::solve -- no solver selected" ) ;
    }
  }
  
  template class BABD<float> ;
  template class BABD<double> ;

}
