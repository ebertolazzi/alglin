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

#include "LU_Arceco.hh"
#include "Alglin.hh"

namespace alglin {

  using namespace std ;

  ArcecoLU::ArcecoLU()
  : baseValue("ArcecoLU_value")
  , baseIndex("ArcecoLU_index")
  { }

  ArcecoLU::~ArcecoLU() {
    baseValue . free() ;
    baseIndex . free() ;
  }

  void
  ArcecoLU::load( indexType    numInitialBc,
                  indexType    numFinalBc,
                  indexType    numInitialETA,
                  indexType    numFinalETA,
                  indexType    numBlock,
                  valuePointer AdAu,
                  valuePointer H0,
                  valuePointer HN,
                  valuePointer Hq ) {

    this -> numInitialETA = numInitialETA ;

    indexType nq = numInitialBc + numFinalBc ;
    indexType q  = numInitialETA + numFinalETA ;
    indexType n  = nq - q ;
    
    nRow0        = numInitialBc ;
    nCol0        = n + numInitialETA ;
    nRowN        = numFinalBc ;
    nColN        = n + numFinalETA ;
    numEquations = n * numBlock + nRow0 + nRowN ;
    NBLOCK       = numBlock + 2 ;

    baseValue . allocate( 2*n*n*numBlock + nRow0 * nCol0 + nRowN * nColN + numEquations ) ;
    baseIndex . allocate( numEquations + 3*NBLOCK) ;

    AR    = baseValue( 2*n*n*numBlock + nRow0 * nCol0 + nRowN * nColN ) ;
    X     = baseValue( numEquations ) ;
    PIVOT = baseIndex( numEquations ) ;
    MTR   = baseIndex( 3*NBLOCK ) ;

    /*
    //  +                     +
    //  |  Ad Au              |
    //  |     Ad Au           |
    //  |          ...        |
    //  |            Ad Au    |
    //  |  H0           HN Hq |
    //  +                     +
    //
    //  H0 = nq x n
    //  HN = nq x n
    //  Hq = nq x q
    //
    //       +-------+       +-------+       +---+
    //       |  0    |       |  NZ   |       |y|0|
    //  H0 = +-------+  HN = +-------+  Hq = +---+
    //       |  NZ   |       |  0    |       |0|x|
    //       |  NZ   |       |       |       |0|x|
    //       +-------+       +-------+       +---+
    //
    //  Primo blocco per ARCECO (numInitialBc) x (n+numInitialETA) 
    //
    //  +-+-------+
    //  |x|  NZ   |
    //  |x|  NZ   |
    //  +-+-------+
    //
    //  Ultimo blocco per ARCECO (numFinalBc) x (n+numFinalETA) 
    //
    //  +-------+-+
    //  |  NZ   |y|
    //  +-------+-+
    //
    */
    
    #define IDX0(I,J) ((I)+(J)*nRow0)
    #define IDXN(I,J) ((I)+(J)*nRowN)

    valuePointer blk0 = AR ;
    valuePointer blkN = AR + 2*n*n*numBlock + nRow0 * nCol0 ;
    
    // Fill structures
    alglin::zero( nRow0 * nCol0, blk0, 1 ) ;
    alglin::zero( nRowN * nColN, blkN, 1 ) ;
    alglin::copy( 2*n*n*numBlock, AdAu, 1, AR + nRow0 * nCol0, 1 ) ;

    indexPointer mtr = MTR ;
    *mtr++ = nRow0 ;
    *mtr++ = nCol0 ;
    *mtr++ = n ;
    for ( indexType i = 0 ; i < numBlock ; ++i ) {
      *mtr++ = n ;
      *mtr++ = 2*n ;
      *mtr++ = n ;
    }
    *mtr++ = nRowN ;
    *mtr++ = nColN ;
    *mtr++ = 0 ;
    
    std::fill( blk0, blk0 + nRow0 * nCol0, 0 ) ;
    std::fill( blkN, blkN + nRowN * nColN, 0 ) ;
    
    //  +-+-------+
    //  |x|  NZ   |
    //  |x|  NZ   |
    //  +-+-------+
    alglin::gecopy( nRow0, numInitialETA, Hq+nRowN+numFinalETA*nq, nq, blk0 + IDX0(0,0),             nRow0 ) ;
    alglin::gecopy( nRow0, n,             H0+nRowN,                nq, blk0 + IDX0(0,numInitialETA), nRow0 ) ;

    //  +-------+-+
    //  |  NZ   |y|
    //  +-------+-+
    alglin::gecopy( nRowN, n,           HN, nq, blkN + IDXN(0,0), nRowN ) ;
    alglin::gecopy( nRowN, numFinalETA, Hq, nq, blkN + IDXN(0,n), nRowN ) ;
    
    arcecoSolver . loadByRef ( NBLOCK, MTR, AR, PIVOT ) ;

  }

  /*
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */
  void
  ArcecoLU::factorize() {
  
    #ifdef USE_F90_ARCECO
    indexType IFLAG ;
    F77NAME(arcedc)(numEquations, AR, MTR, NBLOCK, PIVOT, IFLAG) ;
    ASSERT(IFLAG == 0, "ArcecoLU::factorize()" ) ;
    #else
    arcecoSolver . factorize() ;
    #endif
  }

  /*             _           
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */                       
  void
  ArcecoLU::solve( valuePointer in_out ) {

    alglin::copy( numEquations - nRow0, in_out, 1, X + nRow0, 1 ) ;
    alglin::copy( nRow0, in_out + numEquations - nRow0, 1, X, 1 ) ;

    #ifdef USE_F90_ARCECO
    F77NAME(arcesl)( AR, MTR, NBLOCK, PIVOT, X ) ;
    #else
    arcecoSolver . solve ( X ) ;
    #endif
    
    alglin::copy( numEquations - numInitialETA, X + numInitialETA, 1, in_out, 1 ) ;
    alglin::copy( numInitialETA, X, 1, in_out + numEquations - numInitialETA, 1 ) ;
  }
}
