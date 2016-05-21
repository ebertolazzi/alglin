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

#include "LU_ArcecoSolver.hh"
#include "Alglin.hh"

namespace alglin {

  /*
  //   _                 _ ____        ____       __ 
  //  | | ___   __ _  __| | __ ) _   _|  _ \ ___ / _|
  //  | |/ _ \ / _` |/ _` |  _ \| | | | |_) / _ \ |_ 
  //  | | (_) | (_| | (_| | |_) | |_| |  _ <  __/  _|
  //  |_|\___/ \__,_|\__,_|____/ \__, |_| \_\___|_|  
  //                             |___/               
  */
  void
  Arceco::loadByRef ( indexType    numberOfBlocks,
                      indexPointer matrixStructure,
                      valuePointer array,
                      indexPointer pivot ) {

    this -> numberOfBlocks  = numberOfBlocks  ;
    this -> matrixStructure = matrixStructure ;
    this -> array           = array           ;
    this -> pivot           = pivot           ;
  }

  /*  
  //        _               _     ____  _                   _                  
  //    ___| |__   ___  ___| | __/ ___|| |_ _ __ _   _  ___| |_ _   _ _ __ ___ 
  //   / __| '_ \ / _ \/ __| |/ /\___ \| __| '__| | | |/ __| __| | | | '__/ _ \
  //  | (__| | | |  __/ (__|   <  ___) | |_| |  | |_| | (__| |_| |_| | | |  __/
  //   \___|_| |_|\___|\___|_|\_\|____/ \__|_|   \__,_|\___|\__|\__,_|_|  \___|
  */
  void
  Arceco::checkStructure( indexType neq ) {
  
    ALGLIN_ASSERT( numOverlap(numberOfBlocks-1) == 0,
            "Arceco::checkStructure: numOverlap(" << numberOfBlocks-1 << ") = " <<
            numOverlap(numberOfBlocks-1) << " expected zero!" ) ;

    // check index
    for ( indexType k = 0 ; k < numberOfBlocks ; ++k ) {
      ALGLIN_ASSERT( numCols(k)    >= 1, "Arceco::checkStructure: numCols(" << k << ") = " << numCols(k) << " < 1 " ) ;
      ALGLIN_ASSERT( numRows(k)    >= 1, "Arceco::checkStructure: numRows(" << k << ") = " << numRows(k) << " < 1 " ) ;
      ALGLIN_ASSERT( numOverlap(k) >= 0, "Arceco::checkStructure: numOverlap(" << k << ") = " << numOverlap(k) << " < 0 " ) ;
      ALGLIN_ASSERT( numCols(k)    >= numOverlap(k),
              "Arceco::checkStructure: numCols(" << k << ") = " << numCols(k) << 
              " < numOverlap(" << k << ") = " << numOverlap(k) << "!" ) ; 
    }

    // check ovelapping
    for ( indexType k = 1 ; k < numberOfBlocks ; ++k )
      ALGLIN_ASSERT( numOverlap(k-1) + numOverlap(k) <= numCols(k),
              "Arceco::checkStructure: at block " << k << " three consecutive block overlap" ) ;

    // controlla che i blocchi atraversino la diagonale
    //             c                  c+numCol
    // r           +------------------+
    //             |                  |
    // r+numRow    +------------------+
    indexType r = numRows(0), c = numCols(0) - numOverlap(0) ;
    for ( indexType k = 1 ; k < numberOfBlocks ; ++k ) {
      ALGLIN_ASSERT( c <= r && c+numCols(k) >= r+numRows(k),
              "Arceco::checkStructure: block n. " << k << " do not cross the diagonal" ) ;
      r += numRows(k) ;
      c += numCols(k)-numOverlap(k) ;
    }

    indexType isum1 = 0 ;
    indexType isum2 = 0 ;
    for ( indexType k = 0 ; k < numberOfBlocks ; ++k ) {
      isum1 += numRows(k) ;
      isum2 += numCols(k) - numOverlap(k) ;
    }
    ALGLIN_ASSERT( isum1 == isum2, "Arceco::checkStructure: matrix not squared!\nrow sum = " << isum1 << " column sum = " << isum2 ) ;
    ALGLIN_ASSERT( isum1 == neq,   "Arceco::checkStructure: block dimension = " << isum1 << " different from expected dimension = " << neq ) ;
  }

  /*
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */                                        
  void
  Arceco::factorize () {

    indexType index1 = 0 ;
    indexType indpiv = 0 ;
    indexType numRowsBlock   = numRows(0) ;
    indexType numColsBlock   = numCols(0) ;
    indexType numOverlapCols = numOverlap(0) ;
    indexType numRowsPivot   = numColsBlock - numOverlapCols ;

    /*
    //  +------------------+
    //  |                  |
    //  |                  |
    //  +---------+--------+---------+
    //            |                  |
    //            +------------------+
    */

    rowElimination ( array + index1, numRowsBlock, numColsBlock, numRowsPivot, pivot + indpiv ) ;

    //CHECK_NAN( array + index1,
    //           "Arceco::rowElimination, first block",
    //           numRowsBlock*numColsBlock ) ;

    for ( indexType k = 1 ; k < numberOfBlocks ; ++k ) {
      indpiv += numRowsPivot ;
      indexType index2        = index1 + numRowsBlock * numRowsPivot ;
      indexType index3        = index2 + numRowsBlock * numOverlapCols ;
      indexType numColsPivot  = numRowsBlock - numRowsPivot ;
      indexType numRowsBlock2 = numRows(k) ;

      columnElimination ( array + index2, numRowsBlock,  numOverlapCols, 
                          array + index3, numRowsBlock2, numColsPivot,
                          pivot + indpiv ) ;
      //CHECK_NAN( array + index2,
      //           "Arceco::columnElimination, top block",
      //           numRowsBlock*numOverlapCols ) ;
      //CHECK_NAN( array + index3,
      //          "Arceco::columnElimination, bottom block",
      //           numRowsBlock2*numColsPivot ) ;

      numRowsBlock   = numRowsBlock2 ;
      index1         = index3 + numRowsBlock * numColsPivot ;
      numColsBlock   = numCols(k) - numColsPivot ;
      numOverlapCols = numOverlap(k) ;
      numRowsPivot   = numColsBlock - numOverlapCols ;
      indpiv        += numColsPivot ;

      rowElimination ( array + index1, numRowsBlock, numColsBlock, numRowsPivot, pivot + indpiv );
      //CHECK_NAN( array + index1,
      //           "Arceco::rowElimination",
      //           numRowsBlock*numColsBlock ) ;
    }
  }
  
  /*
  //   _____ _ _           _             _   _             
  //  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __  
  //  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \ 
  //  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
  //  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  */
  
  void
  Arceco::rowElimination ( valuePointer block,
                           indexType    numRowsBlock,
                           indexType    numColsBlock,
                           indexType    numRowsPivot,
                           indexPointer pivot ) {

    #define BLOCK(I,J) block[(I) + (J) * numRowsBlock]

    for ( indexType j = 0 ; j < numRowsPivot ; ++j ) {
      indexType jplus1 = j + 1 ;
      indexType jmax   = j ;
      valueType rowmax = std::abs(BLOCK(j,j)) ;
      for ( indexType i1 = jplus1 ; i1 < numRowsBlock ; ++i1 ) {
        valueType tempiv = std::abs(BLOCK(i1,j)) ;
        if ( tempiv > rowmax ) { rowmax = tempiv ; jmax = i1 ; }
      }

      ALGLIN_ASSERT( rowmax > 0, "Arceco::rowElimination, singular matrix" ) ;

      pivot[j] = jmax ;
      if ( j != jmax )
        for ( indexType j1 = j ; j1 < numColsBlock ; ++j1 )
          std::swap( BLOCK(jmax,j1), BLOCK(j,j1) ) ;

      valueType rowpiv = BLOCK(j,j) ;
      for ( indexType i1 = jplus1 ; i1 < numRowsBlock ; ++i1 ) {
        valueType rowmlt = ( BLOCK(i1,j) /= rowpiv ) ;
        for ( indexType j1 = jplus1 ; j1 < numColsBlock ; ++j1 )
          BLOCK(i1,j1) -= rowmlt * BLOCK(j,j1) ;
      }
    }

    #undef BLOCK

  }

  void
  Arceco::columnElimination ( valuePointer   topblk,
                              indexType      numRowsTopBlock,
                              indexType      numOverlapCols,
                              valuePointer   botblk,
                              indexType      numRowsBottomBlock,
                              indexType      numColsPivot,
                              indexPointer   pivot ) {

    #define TOPBLK(I,J) topblk[(I) + (J) * numRowsTopBlock]
    #define BOTBLK(I,J) botblk[(I) + (J) * numRowsBottomBlock]

    for ( indexType j = 0 ; j < numColsPivot ; ++j ) {
      indexType jplus1 = j + 1 ;
      indexType i      = numRowsTopBlock - numColsPivot + j ;
      indexType jmax   = j ;
      valueType colmax = std::abs(TOPBLK(i,j)) ;
      for ( indexType j1 = jplus1 ; j1 < numOverlapCols; ++j1 ) {
        valueType tempiv = std::abs(TOPBLK(i,j1)) ;
        if ( tempiv > colmax) { colmax = tempiv ; jmax = j1 ; }
      }
      
      ALGLIN_ASSERT( colmax > 0, "Arceco::columnElimination, singular matrix" ) ;

      pivot[j] = jmax ;
      if ( j != jmax ) {
        for ( indexType k = i ; k < numRowsTopBlock    ; ++k ) std::swap(TOPBLK(k,j),TOPBLK(k,jmax)) ;
        for ( indexType k = 0 ; k < numRowsBottomBlock ; ++k ) std::swap(BOTBLK(k,j),BOTBLK(k,jmax)) ;
      }
      valueType colpiv = TOPBLK(i,j) ;
      for ( indexType j1 = jplus1 ; j1 < numOverlapCols ; ++j1 ) {
        valueType colmlt = (TOPBLK(i,j1) /= colpiv) ;
        for ( indexType k = i+1 ; k < numRowsTopBlock    ; ++k ) TOPBLK(k,j1) -= colmlt * TOPBLK(k,j) ;
        for ( indexType k = 0   ; k < numRowsBottomBlock ; ++k ) BOTBLK(k,j1) -= colmlt * BOTBLK(k,j) ;
      }
    }

    #undef TOPBLK
    #undef BOTBLK

  }
  
  /*
  //   ____        _           
  //  / ___|  ___ | |_   _____ 
  //  \___ \ / _ \| \ \ / / _ \
  //   ___) | (_) | |\ V /  __/
  //  |____/ \___/|_| \_/ \___|
  */                         
  void
  Arceco::solve( valuePointer b ) const {
    indexType indpiv         = 0 ;
    indexType indexa         = 0 ;
    indexType numRowsBlock   = numRows(0) ;
    indexType numColsBlock   = numCols(0) ;
    indexType numOverlapCols = numOverlap(0) ;
    indexType numRowsPivot   = numColsBlock - numOverlapCols ;

    forwardElimination( array + indexa, numRowsBlock, numRowsPivot, pivot + indpiv, b + indpiv ) ;

    indexType numColsPivot = 0 ;
    for ( indexType k = 1 ; k < numberOfBlocks ; ++k ) {
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
      forwardElimination ( array + indexa, numRowsBlock, numRowsPivot, pivot + indpiv, b + indpiv ) ;
    }
    // BACKWARD LOOP
    for ( indexType k = numberOfBlocks - 2 ; k >= 0 ; --k ) {

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
      backwardElimination ( array + indexa, numRowsBlock, numColsPivot, numOverlapCols, pivot + indpiv, b + indpiv ) ;

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

  /*  
  //    __                                  _ 
  //   / _| ___  _ ____      ____ _ _ __ __| |
  //  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
  //  |  _| (_) | |   \ V  V / (_| | | | (_| |
  //  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
  */
  /*  
  //   _____ _ _           _             _   _             
  //  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __  
  //  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \ 
  //  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
  //  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  */
  void
  Arceco::forwardElimination ( valuePointer block,
                               indexType    numRowsBlock,
                               indexType    numRowsPivot,
                               indexPointer pivot,
                               valuePointer b ) const {
    valueConstPointer blockI = block ;
    for ( indexType i = 0 ; i < numRowsPivot ; ++i, blockI += numRowsBlock ) {
      indexType pivoti = pivot[i];
      if ( pivoti != i ) std::swap(b[pivoti],b[i]) ;
      valueType bi = b[i] ;
      for ( indexType l = i+1 ; l < numRowsBlock ; ++l ) b[l] -= blockI[l] * bi;
    }
  }
  /*
  //   ____        _       _   _             
  //  / ___|  ___ | |_   _| |_(_) ___  _ __  
  //  \___ \ / _ \| | | | | __| |/ _ \| '_ \ 
  //   ___) | (_) | | |_| | |_| | (_) | | | |
  //  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
  //                                       
  */                                                         
  void
  Arceco::forwardSolution ( valuePointer block,
                            indexType    numRowsBlock,
                            indexType    numColsPivot,
                            indexType    /* numOverlapCols */,
                            valuePointer b ) const {
    indexType    kk      = numRowsBlock - numColsPivot ;
    valuePointer blockJS = block + kk ;
    for ( indexType j = 0 ; j < numColsPivot ; ++j, blockJS += numRowsBlock ) {
      valueType xj = (b[j] /= blockJS[j]) ;
      for ( indexType l = j+1 ; l < numColsPivot ; ++l )
        b[l] -= blockJS[l] * xj ;
    }
  }
  /*
  //   __  __           _ _  __ _           _   _             
  //  |  \/  | ___   __| (_)/ _(_) ___ __ _| |_(_) ___  _ __  
  //  | |\/| |/ _ \ / _` | | |_| |/ __/ _` | __| |/ _ \| '_ \ 
  //  | |  | | (_) | (_| | |  _| | (_| (_| | |_| | (_) | | | |
  //  |_|  |_|\___/ \__,_|_|_| |_|\___\__,_|\__|_|\___/|_| |_|
  */
  void
  Arceco::forwardModification ( valuePointer block,
                                indexType    numRowsBlock,
                                indexType    numColsPivot,
                                valuePointer b ) const {
    valuePointer blockJ = block ;
    for ( indexType j = 0 ; j < numColsPivot ; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j];
      for ( indexType l = 0 ; l < numRowsBlock ; ++l )
        b[numColsPivot + l] -= blockJ[l] * xj ;
    }
  }
  /*  
  //   _                _                           _ 
  //  | |__   __ _  ___| | ____      ____ _ _ __ __| |
  //  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
  //  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
  //  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
  */  
  /*
  //   __  __           _ _  __ _           _   _             
  //  |  \/  | ___   __| (_)/ _(_) ___ __ _| |_(_) ___  _ __  
  //  | |\/| |/ _ \ / _` | | |_| |/ __/ _` | __| |/ _ \| '_ \ 
  //  | |  | | (_) | (_| | |  _| | (_| (_| | |_| | (_) | | | |
  //  |_|  |_|\___/ \__,_|_|_| |_|\___\__,_|\__|_|\___/|_| |_|
  */
  void
  Arceco::backwardModification ( valuePointer block,
                                 indexType    numRowsBlock,
                                 indexType    numColsBlock,
                                 indexType    numRowsPivot,
                                 valuePointer b ) const {
    valuePointer blockJ = block + numRowsPivot * numRowsBlock ;
    for ( indexType j = numRowsPivot ; j < numColsBlock ; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j] ;
      for ( indexType l = 0 ; l < numRowsPivot ; ++l )
        b[l] -= blockJ[l] * xj ;
    }
  }
  /*
  //   ____        _       _   _             
  //  / ___|  ___ | |_   _| |_(_) ___  _ __  
  //  \___ \ / _ \| | | | | __| |/ _ \| '_ \ 
  //   ___) | (_) | | |_| | |_| | (_) | | | |
  //  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
  //                                       
  */                                                         
  void
  Arceco::backwardSolution ( valuePointer block,
                             indexType    numRowsBlock,
                             indexType    /* numColsBlock */,
                             indexType    numRowsPivot,
                             valuePointer b ) const {
    for ( indexType j = numRowsPivot - 1 ; j >= 0 ; --j ) {
      valuePointer blockJ = block + j * numRowsBlock ;
      valueType xj = (b[j] /= blockJ[j]) ;
      for ( indexType l = 0 ; l < j ; ++l )
        b[l] -= blockJ[l] * xj ;
    }
  }
  /*  
  //   _____ _ _           _             _   _             
  //  | ____| (_)_ __ ___ (_)_ __   __ _| |_(_) ___  _ __  
  //  |  _| | | | '_ ` _ \| | '_ \ / _` | __| |/ _ \| '_ \ 
  //  | |___| | | | | | | | | | | | (_| | |_| | (_) | | | |
  //  |_____|_|_|_| |_| |_|_|_| |_|\__,_|\__|_|\___/|_| |_|
  */
  void
  Arceco::backwardElimination( valuePointer block,
                               indexType    numRowsBlock,
                               indexType    numColsPivot,
                               indexType    numOverlapCols,
                               indexPointer pivot,
                               valuePointer x )  const {
    indexType kk = numRowsBlock - numColsPivot ;
    indexType j1 = numColsPivot ;
    while ( j1 > 0 ) {
      valueConstPointer blockS = block + j1-1 + kk + j1 * numRowsBlock ;
      valueConstPointer xS     = x + j1 ;
      indexType n = numOverlapCols - j1 ;
      // Kahan summation algorithm
      valueType dotprd = 0 ;
      valueType comp   = 0 ;
      for ( indexType k = 0 ; k < n ; ++k, ++xS, blockS += numRowsBlock ) {
        valueType y = xS[0] * blockS[0] - comp ;
        valueType t = dotprd + y ;
        comp = (t-dotprd) - y ;
        dotprd = t ;
      }
      --j1 ;
      x[j1] -= dotprd ;
      indexType pivotj = pivot[j1] ;
      if ( pivotj != j1 ) std::swap( x[pivotj], x[j1] ) ;
    }
  }

}
