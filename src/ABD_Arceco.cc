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

#include "ABD_Arceco.hh"
#include "Alglin.hh"

namespace alglin {

  /*\
   |   _                 _ ____        ____       __
   |  | | ___   __ _  __| | __ ) _   _|  _ \ ___ / _|
   |  | |/ _ \ / _` |/ _` |  _ \| | | | |_) / _ \ |_
   |  | | (_) | (_| | (_| | |_) | |_| |  _ <  __/  _|
   |  |_|\___/ \__,_|\__,_|____/ \__, |_| \_\___|_|
   |                             |___/
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::loadByRef ( integer      _numberOfBlocks,
                                 integer *    _matrixStructure,
                                 valuePointer _array,
                                 integer *    _pivot ) {

    numberOfBlocks  = _numberOfBlocks  ;
    matrixStructure = _matrixStructure ;
    array           = _array           ;
    pivot_array     = _pivot           ;
  }

  /*\
   |        _               _     ____  _                   _
   |    ___| |__   ___  ___| | __/ ___|| |_ _ __ _   _  ___| |_ _   _ _ __ ___
   |   / __| '_ \ / _ \/ __| |/ /\___ \| __| '__| | | |/ __| __| | | | '__/ _ \
   |  | (__| | | |  __/ (__|   <  ___) | |_| |  | |_| | (__| |_| |_| | | |  __/
   |   \___|_| |_|\___|\___|_|\_\|____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::checkStructure( integer neq ) {
  
    ALGLIN_ASSERT( numOverlap(numberOfBlocks-1) == 0,
                   "Arceco::checkStructure: numOverlap(" << numberOfBlocks-1 << ") = " <<
                    numOverlap(numberOfBlocks-1) << " expected zero!" ) ;

    // check index
    for ( integer k = 0 ; k < numberOfBlocks ; ++k ) {
      ALGLIN_ASSERT( numCols(k)    >= 1, "Arceco::checkStructure: numCols(" << k << ") = " << numCols(k) << " < 1 " ) ;
      ALGLIN_ASSERT( numRows(k)    >= 1, "Arceco::checkStructure: numRows(" << k << ") = " << numRows(k) << " < 1 " ) ;
      ALGLIN_ASSERT( numOverlap(k) >= 0, "Arceco::checkStructure: numOverlap(" << k << ") = " << numOverlap(k) << " < 0 " ) ;
      ALGLIN_ASSERT( numCols(k)    >= numOverlap(k),
              "Arceco::checkStructure: numCols(" << k << ") = " << numCols(k) << 
              " < numOverlap(" << k << ") = " << numOverlap(k) << "!" ) ; 
    }

    // check ovelapping
    for ( integer k = 1 ; k < numberOfBlocks ; ++k )
      ALGLIN_ASSERT( numOverlap(k-1) + numOverlap(k) <= numCols(k),
                     "Arceco::checkStructure: at block " << k << " three consecutive block overlap" ) ;

    // controlla che i blocchi atraversino la diagonale
    //             c                  c+numCol
    // r           +------------------+
    //             |                  |
    // r+numRow    +------------------+
    integer r = numRows(0), c = numCols(0) - numOverlap(0) ;
    for ( integer k = 1 ; k < numberOfBlocks ; ++k ) {
      ALGLIN_ASSERT( c <= r && c+numCols(k) >= r+numRows(k),
              "Arceco::checkStructure: block n. " << k << " do not cross the diagonal" ) ;
      r += numRows(k) ;
      c += numCols(k)-numOverlap(k) ;
    }

    integer isum1 = 0 ;
    integer isum2 = 0 ;
    for ( integer k = 0 ; k < numberOfBlocks ; ++k ) {
      isum1 += numRows(k) ;
      isum2 += numCols(k) - numOverlap(k) ;
    }
    ALGLIN_ASSERT( isum1 == isum2, "Arceco::checkStructure: matrix not squared!\nrow sum = " << isum1 << " column sum = " << isum2 ) ;
    ALGLIN_ASSERT( isum1 == neq,   "Arceco::checkStructure: block dimension = " << isum1 << " different from expected dimension = " << neq ) ;
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  //! factorize the matrix
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::factorize( integer           row0,
                                integer           col0,
                                valueConstPointer block0,
                                integer           nblk,
                                integer           n,
                                valueConstPointer blocks,
                                integer           rowN,
                                integer           colN,
                                valueConstPointer blockN ) {

    numberOfBlocks = nblk + 2 ;

    integer size0        = row0 * col0 ;
    integer sizeN        = rowN * colN ;
    integer BLK_size     = 2*n*n*nblk ;
    integer numEquations = nblk * n + row0 + rowN ;

    // allocazione dinamica
    baseValue   . allocate( size_t( BLK_size + size0 + sizeN ) ) ;
    baseInteger . allocate( size_t( 3*numberOfBlocks + numEquations ) ) ;

    array           = baseValue( size_t( BLK_size + size0 + sizeN ) ) ;
    pivot_array     = baseInteger( size_t( numEquations ) ) ;
    matrixStructure = baseInteger( size_t( 3*numberOfBlocks ) ) ;

    // Fill structures
    alglin::copy( size0,    block0, 1, array, 1 ) ;
    alglin::copy( BLK_size, blocks, 1, array + size0, 1 ) ;
    alglin::copy( sizeN,    blockN, 1, array + size0 + BLK_size, 1 ) ;

    integer * mtr = matrixStructure ;
    *mtr++ = row0 ;
    *mtr++ = col0 ;
    *mtr++ = n ;
    for ( integer i = 0 ; i < nblk ; ++i ) {
      *mtr++ = n ;
      *mtr++ = 2*n ;
      *mtr++ = n ;
    }
    *mtr++ = rowN ;
    *mtr++ = colN ;
    *mtr++ = 0 ;

    factorize() ;
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::factorize () {

    integer index1 = 0 ;
    integer indpiv = 0 ;
    integer numRowsBlock   = numRows(0) ;
    integer numColsBlock   = numCols(0) ;
    integer numOverlapCols = numOverlap(0) ;
    integer numRowsPivot   = numColsBlock - numOverlapCols ;

    /*\
     |  +------------------+
     |  |                  |
     |  |                  |
     |  +---------+--------+---------+
     |            |                  |
     |            +------------------+
    \*/

    rowElimination ( array + index1,
                     numRowsBlock, numColsBlock, numRowsPivot,
                     pivot_array + indpiv ) ;

    for ( integer k = 1 ; k < numberOfBlocks ; ++k ) {
      indpiv += numRowsPivot ;
      integer index2        = index1 + numRowsBlock * numRowsPivot ;
      integer index3        = index2 + numRowsBlock * numOverlapCols ;
      integer numColsPivot  = numRowsBlock - numRowsPivot ;
      integer numRowsBlock2 = numRows(k) ;

      columnElimination ( array + index2, numRowsBlock,  numOverlapCols, 
                          array + index3, numRowsBlock2, numColsPivot,
                          pivot_array + indpiv ) ;

      numRowsBlock   = numRowsBlock2 ;
      index1         = index3 + numRowsBlock * numColsPivot ;
      numColsBlock   = numCols(k) - numColsPivot ;
      numOverlapCols = numOverlap(k) ;
      numRowsPivot   = numColsBlock - numOverlapCols ;
      indpiv        += numColsPivot ;

      rowElimination ( array + index1,
                       numRowsBlock, numColsBlock, numRowsPivot,
                       pivot_array + indpiv );

    }
  }
  
  /*\
   |   _____ _ _           _             _   _
   |  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __
   |  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \
   |  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
   |  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::rowElimination ( valuePointer block,
                                      integer      numRowsBlock,
                                      integer      numColsBlock,
                                      integer      numRowsPivot,
                                      integer *    pivot ) {

    #define BLOCK(I,J) block[(I) + (J) * numRowsBlock]

    for ( integer j = 0 ; j < numRowsPivot ; ++j ) {
      integer   jplus1 = j + 1 ;
      integer   jmax   = j ;
      valueType rowmax = std::abs(BLOCK(j,j)) ;
      for ( integer i1 = jplus1 ; i1 < numRowsBlock ; ++i1 ) {
        valueType tempiv = std::abs(BLOCK(i1,j)) ;
        if ( tempiv > rowmax ) { rowmax = tempiv ; jmax = i1 ; }
      }

      ALGLIN_ASSERT( rowmax > 0, "Arceco::rowElimination, singular matrix" ) ;

      pivot[j] = jmax ;
      if ( j != jmax )
        for ( integer j1 = j ; j1 < numColsBlock ; ++j1 )
          std::swap( BLOCK(jmax,j1), BLOCK(j,j1) ) ;

      valueType rowpiv = BLOCK(j,j) ;
      for ( integer i1 = jplus1 ; i1 < numRowsBlock ; ++i1 ) {
        valueType rowmlt = ( BLOCK(i1,j) /= rowpiv ) ;
        for ( integer j1 = jplus1 ; j1 < numColsBlock ; ++j1 )
          BLOCK(i1,j1) -= rowmlt * BLOCK(j,j1) ;
      }
    }

    #undef BLOCK

  }

  template <typename t_Value>
  void
  ArcecoLU<t_Value>::columnElimination ( valuePointer topblk,
                                         integer      numRowsTopBlock,
                                         integer      numOverlapCols,
                                         valuePointer botblk,
                                         integer      numRowsBottomBlock,
                                         integer      numColsPivot,
                                         integer *    pivot ) {

    #define TOPBLK(I,J) topblk[(I) + (J) * numRowsTopBlock]
    #define BOTBLK(I,J) botblk[(I) + (J) * numRowsBottomBlock]

    for ( integer j = 0 ; j < numColsPivot ; ++j ) {
      integer jplus1 = j + 1 ;
      integer i      = numRowsTopBlock - numColsPivot + j ;
      integer jmax   = j ;
      valueType colmax = std::abs(TOPBLK(i,j)) ;
      for ( integer j1 = jplus1 ; j1 < numOverlapCols; ++j1 ) {
        valueType tempiv = std::abs(TOPBLK(i,j1)) ;
        if ( tempiv > colmax) { colmax = tempiv ; jmax = j1 ; }
      }
      
      ALGLIN_ASSERT( colmax > 0, "Arceco::columnElimination, singular matrix" ) ;

      pivot[j] = jmax ;
      if ( j != jmax ) {
        for ( integer k = i ; k < numRowsTopBlock    ; ++k ) std::swap(TOPBLK(k,j),TOPBLK(k,jmax)) ;
        for ( integer k = 0 ; k < numRowsBottomBlock ; ++k ) std::swap(BOTBLK(k,j),BOTBLK(k,jmax)) ;
      }
      valueType colpiv = TOPBLK(i,j) ;
      for ( integer j1 = jplus1 ; j1 < numOverlapCols ; ++j1 ) {
        valueType colmlt = (TOPBLK(i,j1) /= colpiv) ;
        for ( integer k = i+1 ; k < numRowsTopBlock    ; ++k ) TOPBLK(k,j1) -= colmlt * TOPBLK(k,j) ;
        for ( integer k = 0   ; k < numRowsBottomBlock ; ++k ) BOTBLK(k,j1) -= colmlt * BOTBLK(k,j) ;
      }
    }

    #undef TOPBLK
    #undef BOTBLK

  }
  
  /*\
   |   ____        _
   |  / ___|  ___ | |_   _____
   |  \___ \ / _ \| \ \ / / _ \
   |   ___) | (_) | |\ V /  __/
   |  |____/ \___/|_| \_/ \___|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::solve( valuePointer b ) const {
    integer indpiv         = 0 ;
    integer indexa         = 0 ;
    integer numRowsBlock   = numRows(0) ;
    integer numColsBlock   = numCols(0) ;
    integer numOverlapCols = numOverlap(0) ;
    integer numRowsPivot   = numColsBlock - numOverlapCols ;

    forwardElimination( array + indexa,
                        numRowsBlock, numRowsPivot,
                        pivot_array + indpiv, b + indpiv ) ;

    integer numColsPivot = 0 ;
    for ( integer k = 1 ; k < numberOfBlocks ; ++k ) {
      indexa      += numRowsBlock * numRowsPivot ;
      numColsPivot = numRowsBlock - numRowsPivot ;
      indpiv      += numRowsPivot ;

      //if ( numColsPivot > 0 ) 
      forwardSolution ( array + indexa, numRowsBlock, numColsPivot, numOverlapCols, b + indpiv ) ;

      indexa      += numOverlapCols * numRowsBlock ;
      numRowsBlock = numRows(k) ;
      
      //if ( numColsPivot > 0 ) 
      forwardModification ( array + indexa, numRowsBlock, numColsPivot, b + indpiv ) ;
      
      indexa        += numRowsBlock * numColsPivot ;
      numColsBlock   = numCols(k) - numColsPivot ;
      numOverlapCols = numOverlap(k) ;
      numRowsPivot   = numColsBlock - numOverlapCols ;
      indpiv        += numColsPivot;
      
      //if ( numRowsPivot > 0 )  
      forwardElimination ( array + indexa,
                           numRowsBlock, numRowsPivot,
                           pivot_array + indpiv, b + indpiv ) ;
    }
    // BACKWARD LOOP
    for ( integer k = numberOfBlocks - 2 ; k >= 0 ; --k ) {

      if ( numRowsPivot != 0 ) {
        if ( numRowsPivot != numColsBlock ) backwardModification ( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv ) ;
        backwardSolution ( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv ) ;
      }

      indexa        -= numRowsBlock * numColsPivot ;
      numRowsBlock   = numRows(k) ;
      numOverlapCols = numOverlap(k) ;
      indexa        -= numRowsBlock * numOverlapCols ;
      indpiv        -= numColsPivot ;

      //if ( numColsPivot > 0 ) 
      backwardElimination ( array + indexa,
                            numRowsBlock, numColsPivot, numOverlapCols,
                            pivot_array + indpiv, b + indpiv ) ;

      numRowsPivot = numRowsBlock - numColsPivot ;
      numColsBlock = numOverlapCols + numRowsPivot ;
      indexa      -= numRowsBlock * numRowsPivot ;
      indpiv      -= numRowsPivot ;
      numColsPivot = numCols(k) - numColsBlock ;
    }

    if ( numRowsPivot > 0 ) {
      if ( numRowsPivot != numColsBlock ) backwardModification( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv ) ;
      backwardSolution( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv ) ;
    }
  }

  /*\
   |    __                                  _
   |   / _| ___  _ ____      ____ _ _ __ __| |
   |  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
   |  |  _| (_) | |   \ V  V / (_| | | | (_| |
   |  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
   |
   |
   |   _____ _ _           _             _   _
   |  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __
   |  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \
   |  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
   |  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::forwardElimination ( valuePointer block,
                                          integer      numRowsBlock,
                                          integer      numRowsPivot,
                                          integer *    pivot,
                                          valuePointer b ) const {
    valueConstPointer blockI = block ;
    for ( integer i = 0 ; i < numRowsPivot ; ++i, blockI += numRowsBlock ) {
      integer pivoti = pivot[i];
      if ( pivoti != i ) std::swap(b[pivoti],b[i]) ;
      valueType bi = b[i] ;
      for ( integer l = i+1 ; l < numRowsBlock ; ++l ) b[l] -= blockI[l] * bi;
    }
  }

  /*\
   |   ____        _       _   _
   |  / ___|  ___ | |_   _| |_(_) ___  _ __
   |  \___ \ / _ \| | | | | __| |/ _ \| '_ \
   |   ___) | (_) | | |_| | |_| | (_) | | | |
   |  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
   |
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::forwardSolution ( valuePointer block,
                                       integer      numRowsBlock,
                                       integer      numColsPivot,
                                       integer      /* numOverlapCols */,
                                       valuePointer b ) const {
    integer      kk      = numRowsBlock - numColsPivot ;
    valuePointer blockJS = block + kk ;
    for ( integer j = 0 ; j < numColsPivot ; ++j, blockJS += numRowsBlock ) {
      valueType xj = (b[j] /= blockJS[j]) ;
      for ( integer l = j+1 ; l < numColsPivot ; ++l )
        b[l] -= blockJS[l] * xj ;
    }
  }

  /*\
   |   __  __           _ _  __ _           _   _
   |  |  \/  | ___   __| (_)/ _(_) ___ __ _| |_(_) ___  _ __
   |  | |\/| |/ _ \ / _` | | |_| |/ __/ _` | __| |/ _ \| '_ \
   |  | |  | | (_) | (_| | |  _| | (_| (_| | |_| | (_) | | | |
   |  |_|  |_|\___/ \__,_|_|_| |_|\___\__,_|\__|_|\___/|_| |_|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::forwardModification ( valuePointer block,
                                           integer      numRowsBlock,
                                           integer      numColsPivot,
                                           valuePointer b ) const {
    valuePointer blockJ = block ;
    for ( integer j = 0 ; j < numColsPivot ; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j];
      for ( integer l = 0 ; l < numRowsBlock ; ++l )
        b[numColsPivot + l] -= blockJ[l] * xj ;
    }
  }
  /*\
   |   _                _                           _
   |  | |__   __ _  ___| | ____      ____ _ _ __ __| |
   |  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
   |  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
   |  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
   |   __  __           _ _  __ _           _   _
   |  |  \/  | ___   __| (_)/ _(_) ___ __ _| |_(_) ___  _ __
   |  | |\/| |/ _ \ / _` | | |_| |/ __/ _` | __| |/ _ \| '_ \
   |  | |  | | (_) | (_| | |  _| | (_| (_| | |_| | (_) | | | |
   |  |_|  |_|\___/ \__,_|_|_| |_|\___\__,_|\__|_|\___/|_| |_|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::backwardModification ( valuePointer block,
                                            integer      numRowsBlock,
                                            integer      numColsBlock,
                                            integer      numRowsPivot,
                                            valuePointer b ) const {
    valuePointer blockJ = block + numRowsPivot * numRowsBlock ;
    for ( integer j = numRowsPivot ; j < numColsBlock ; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j] ;
      for ( integer l = 0 ; l < numRowsPivot ; ++l )
        b[l] -= blockJ[l] * xj ;
    }
  }
  /*\
   |   ____        _       _   _
   |  / ___|  ___ | |_   _| |_(_) ___  _ __
   |  \___ \ / _ \| | | | | __| |/ _ \| '_ \
   |   ___) | (_) | | |_| | |_| | (_) | | | |
   |  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
   |
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::backwardSolution ( valuePointer block,
                                        integer      numRowsBlock,
                                        integer      /* numColsBlock */,
                                        integer      numRowsPivot,
                                        valuePointer b ) const {
    for ( integer j = numRowsPivot - 1 ; j >= 0 ; --j ) {
      valuePointer blockJ = block + j * numRowsBlock ;
      valueType xj = (b[j] /= blockJ[j]) ;
      for ( integer l = 0 ; l < j ; ++l )
        b[l] -= blockJ[l] * xj ;
    }
  }
  /*\
   |   _____ _ _           _             _   _
   |  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __
   |  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \
   |  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
   |  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  \*/
  template <typename t_Value>
  void
  ArcecoLU<t_Value>::backwardElimination( valuePointer block,
                                          integer      numRowsBlock,
                                          integer      numColsPivot,
                                          integer      numOverlapCols,
                                          integer *    pivot,
                                          valuePointer x )  const {
    integer kk = numRowsBlock - numColsPivot ;
    integer j1 = numColsPivot ;
    while ( j1 > 0 ) {
      valueConstPointer blockS = block + j1-1 + kk + j1 * numRowsBlock ;
      integer n = numOverlapCols - j1 ;
      // Kahan summation algorithm
      valueType dotprd = dot(n,x+j1,1,blockS,numRowsBlock) ;
      x[--j1] -= dotprd ;
      integer pivotj = pivot[j1] ;
      if ( pivotj != j1 ) std::swap( x[pivotj], x[j1] ) ;
    }
  }

  #ifdef __GCC__
  #pragma GCC diagnostic ignored "-Wpadded"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic ignored "-Wpadded"
  #endif

  template class ArcecoLU<float> ;
  template class ArcecoLU<double> ;

}
