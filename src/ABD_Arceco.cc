/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
  ArcecoLU<t_Value>::loadByRef(
    integer   _numberOfBlocks,
    integer   _matrixStructure[],
    valueType _array[],
    integer   _pivot[]
  ) {
    numberOfBlocks  = _numberOfBlocks;
    matrixStructure = _matrixStructure;
    array           = _array;
    pivot_array     = _pivot;
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

    UTILS_ASSERT(
      numOverlap(numberOfBlocks-1) == 0,
      "Arceco::checkStructure: numOverlap({}) = {} expected zero!\n",
      numberOfBlocks-1, numOverlap(numberOfBlocks-1)
    );

    // check index
    for ( integer k = 0; k < numberOfBlocks; ++k ) {
      UTILS_ASSERT(
        numCols(k) >= 1,
        "ArcecoLU::checkStructure: numCols({}) = {} < 1\n", k, numCols(k)
      );
      UTILS_ASSERT(
        numRows(k) >= 1,
        "ArcecoLU::checkStructure: numRows({}) = {} < 1\n", k, numRows(k)
      );
      UTILS_ASSERT(
        numOverlap(k) >= 0,
        "ArcecoLU::checkStructure: numOverlap({}) = {} < 0\n", k, numOverlap(k)
      );
      UTILS_ASSERT(
        numCols(k) >= numOverlap(k),
        "ArcecoLU::checkStructure: numCols({}) = {} < numOverlap({}) = {}\n",
        k, numCols(k), k, numOverlap(k)
      );
    }

    // check ovelapping
    for ( integer k = 1; k < numberOfBlocks; ++k )
      UTILS_ASSERT(
        numOverlap(k-1) + numOverlap(k) <= numCols(k),
        "Arceco::checkStructure: at block {} three consecutive block overlap\n", k
      );

    // controlla che i blocchi atraversino la diagonale
    //             c                  c+numCol
    // r           +------------------+
    //             |                  |
    // r+numRow    +------------------+
    integer r = numRows(0), c = numCols(0) - numOverlap(0);
    for ( integer k = 1; k < numberOfBlocks; ++k ) {
      UTILS_ASSERT(
        c <= r && c+numCols(k) >= r+numRows(k),
        "ArcecoLU::checkStructure: block n. {} do not cross the diagonal\n", k
      );
      r += numRows(k);
      c += numCols(k)-numOverlap(k);
    }

    integer isum1 = 0;
    integer isum2 = 0;
    for ( integer k = 0; k < numberOfBlocks; ++k ) {
      isum1 += numRows(k);
      isum2 += numCols(k) - numOverlap(k);
    }
    UTILS_ASSERT(
      isum1 == isum2,
      "ArcecoLU::checkStructure: matrix not squared!\n"
      "row sum = {} column sum = {}\n",
      isum1, isum2
    );
    UTILS_ASSERT(
      isum1 == neq,
      "ArcecoLU::checkStructure: block dimension = {}"
      " different from expected dimension = {}\n",
      isum1, neq
    );
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
  ArcecoLU<t_Value>::factorize(
    integer         row0,
    integer         col0,
    valueType const block0[],
    integer         nblk,
    integer         n,
    valueType const blocks[],
    integer         rowN,
    integer         colN,
    valueType const blockN[]
  ) {

    numberOfBlocks = nblk + 2;

    integer size0        = row0 * col0;
    integer sizeN        = rowN * colN;
    integer BLK_size     = 2*n*n*nblk;
    integer numEquations = nblk * n + row0 + rowN;

    // allocazione dinamica
    baseValue   . allocate( size_t( BLK_size + size0 + sizeN ) );
    baseInteger . allocate( size_t( 3*numberOfBlocks + numEquations ) );

    array           = baseValue( size_t( BLK_size + size0 + sizeN ) );
    pivot_array     = baseInteger( size_t( numEquations ) );
    matrixStructure = baseInteger( size_t( 3*numberOfBlocks ) );

    // Fill structures
    alglin::copy( size0,    block0, 1, array, 1 );
    alglin::copy( BLK_size, blocks, 1, array + size0, 1 );
    alglin::copy( sizeN,    blockN, 1, array + size0 + BLK_size, 1 );

    integer * mtr = matrixStructure;
    *mtr++ = row0;
    *mtr++ = col0;
    *mtr++ = n;
    for ( integer i = 0; i < nblk; ++i ) {
      *mtr++ = n;
      *mtr++ = 2*n;
      *mtr++ = n;
    }
    *mtr++ = rowN;
    *mtr++ = colN;
    *mtr++ = 0;

    factorize();
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
  ArcecoLU<t_Value>::factorize() {

    integer index1 = 0;
    integer indpiv = 0;
    integer numRowsBlock   = numRows(0);
    integer numColsBlock   = numCols(0);
    integer numOverlapCols = numOverlap(0);
    integer numRowsPivot   = numColsBlock - numOverlapCols;

    /*\
     |  +------------------+
     |  |                  |
     |  |                  |
     |  +---------+--------+---------+
     |            |                  |
     |            +------------------+
    \*/

    rowElimination( array + index1,
                    numRowsBlock, numColsBlock, numRowsPivot,
                    pivot_array + indpiv );

    for ( integer k = 1; k < numberOfBlocks; ++k ) {
      indpiv += numRowsPivot;
      integer index2        = index1 + numRowsBlock * numRowsPivot;
      integer index3        = index2 + numRowsBlock * numOverlapCols;
      integer numColsPivot  = numRowsBlock - numRowsPivot;
      integer numRowsBlock2 = numRows(k);

      columnElimination(
        array + index2, numRowsBlock,  numOverlapCols,
        array + index3, numRowsBlock2, numColsPivot,
        pivot_array + indpiv
      );

      numRowsBlock   = numRowsBlock2;
      index1         = index3 + numRowsBlock * numColsPivot;
      numColsBlock   = numCols(k) - numColsPivot;
      numOverlapCols = numOverlap(k);
      numRowsPivot   = numColsBlock - numOverlapCols;
      indpiv        += numColsPivot;

      rowElimination(
        array + index1,
        numRowsBlock, numColsBlock, numRowsPivot,
        pivot_array + indpiv
      );

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
  ArcecoLU<t_Value>::rowElimination(
    valueType block[],
    integer   numRowsBlock,
    integer   numColsBlock,
    integer   numRowsPivot,
    integer   pivot[]
  ) {

    #define BLOCK(I,J) block[(I) + (J) * numRowsBlock]

    for ( integer j = 0; j < numRowsPivot; ++j ) {
      integer   jplus1 = j + 1;
      integer   jmax   = j;
      valueType rowmax = std::abs(BLOCK(j,j));
      for ( integer i1 = jplus1; i1 < numRowsBlock; ++i1 ) {
        valueType tempiv = std::abs(BLOCK(i1,j));
        if ( tempiv > rowmax ) { rowmax = tempiv; jmax = i1; }
      }

      UTILS_ASSERT0( rowmax > 0, "Arceco::rowElimination, singular matrix\n" );

      pivot[j] = jmax;
      if ( j != jmax )
        for ( integer j1 = j; j1 < numColsBlock; ++j1 )
          std::swap( BLOCK(jmax,j1), BLOCK(j,j1) );

      valueType rowpiv = BLOCK(j,j);
      for ( integer i1 = jplus1; i1 < numRowsBlock; ++i1 ) {
        valueType rowmlt = ( BLOCK(i1,j) /= rowpiv );
        for ( integer j1 = jplus1; j1 < numColsBlock; ++j1 )
          BLOCK(i1,j1) -= rowmlt * BLOCK(j,j1);
      }
    }

    #undef BLOCK

  }

  template <typename t_Value>
  void
  ArcecoLU<t_Value>::columnElimination(
    valueType topblk[],
    integer   numRowsTopBlock,
    integer   numOverlapCols,
    valueType botblk[],
    integer   numRowsBottomBlock,
    integer   numColsPivot,
    integer   pivot[]
  ) {

    #define TOPBLK(I,J) topblk[(I) + (J) * numRowsTopBlock]
    #define BOTBLK(I,J) botblk[(I) + (J) * numRowsBottomBlock]

    for ( integer j = 0; j < numColsPivot; ++j ) {
      integer jplus1 = j + 1;
      integer i      = numRowsTopBlock - numColsPivot + j;
      integer jmax   = j;
      valueType colmax = std::abs(TOPBLK(i,j));
      for ( integer j1 = jplus1; j1 < numOverlapCols; ++j1 ) {
        valueType tempiv = std::abs(TOPBLK(i,j1));
        if ( tempiv > colmax) { colmax = tempiv; jmax = j1; }
      }

      UTILS_ASSERT0( colmax > 0, "Arceco::columnElimination, singular matrix\n" );

      pivot[j] = jmax;
      if ( j != jmax ) {
        for ( integer k = i; k < numRowsTopBlock;    ++k ) std::swap(TOPBLK(k,j),TOPBLK(k,jmax));
        for ( integer k = 0; k < numRowsBottomBlock; ++k ) std::swap(BOTBLK(k,j),BOTBLK(k,jmax));
      }
      valueType colpiv = TOPBLK(i,j);
      for ( integer j1 = jplus1; j1 < numOverlapCols; ++j1 ) {
        valueType colmlt = (TOPBLK(i,j1) /= colpiv);
        for ( integer k = i+1; k < numRowsTopBlock;    ++k ) TOPBLK(k,j1) -= colmlt * TOPBLK(k,j);
        for ( integer k = 0;   k < numRowsBottomBlock; ++k ) BOTBLK(k,j1) -= colmlt * BOTBLK(k,j);
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
  ArcecoLU<t_Value>::solve( valueType b[] ) const {
    integer indpiv         = 0;
    integer indexa         = 0;
    integer numRowsBlock   = numRows(0);
    integer numColsBlock   = numCols(0);
    integer numOverlapCols = numOverlap(0);
    integer numRowsPivot   = numColsBlock - numOverlapCols;

    forwardElimination(
      array + indexa,
      numRowsBlock, numRowsPivot,
      pivot_array + indpiv, b + indpiv
    );

    integer numColsPivot = 0;
    for ( integer k = 1; k < numberOfBlocks; ++k ) {
      indexa      += numRowsBlock * numRowsPivot;
      numColsPivot = numRowsBlock - numRowsPivot;
      indpiv      += numRowsPivot;

      //if ( numColsPivot > 0 ) 
      forwardSolution(
        array + indexa, numRowsBlock, numColsPivot, numOverlapCols, b + indpiv
      );

      indexa      += numOverlapCols * numRowsBlock;
      numRowsBlock = numRows(k);
      
      //if ( numColsPivot > 0 ) 
      forwardModification(
        array + indexa, numRowsBlock, numColsPivot, b + indpiv
      );
      
      indexa        += numRowsBlock * numColsPivot;
      numColsBlock   = numCols(k) - numColsPivot;
      numOverlapCols = numOverlap(k);
      numRowsPivot   = numColsBlock - numOverlapCols;
      indpiv        += numColsPivot;
      
      //if ( numRowsPivot > 0 )  
      forwardElimination(
        array + indexa,
        numRowsBlock, numRowsPivot,
        pivot_array + indpiv, b + indpiv
      );
    }
    // BACKWARD LOOP
    for ( integer k = numberOfBlocks - 2; k >= 0; --k ) {

      if ( numRowsPivot != 0 ) {
        if ( numRowsPivot != numColsBlock ) backwardModification ( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv );
        backwardSolution ( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv );
      }

      indexa        -= numRowsBlock * numColsPivot;
      numRowsBlock   = numRows(k);
      numOverlapCols = numOverlap(k);
      indexa        -= numRowsBlock * numOverlapCols;
      indpiv        -= numColsPivot;

      //if ( numColsPivot > 0 ) 
      backwardElimination(
        array + indexa,
        numRowsBlock, numColsPivot, numOverlapCols,
        pivot_array + indpiv, b + indpiv
      );

      numRowsPivot = numRowsBlock - numColsPivot;
      numColsBlock = numOverlapCols + numRowsPivot;
      indexa      -= numRowsBlock * numRowsPivot;
      indpiv      -= numRowsPivot;
      numColsPivot = numCols(k) - numColsBlock;
    }

    if ( numRowsPivot > 0 ) {
      if ( numRowsPivot != numColsBlock ) backwardModification( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv );
      backwardSolution( array + indexa, numRowsBlock, numColsBlock, numRowsPivot, b + indpiv );
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
  ArcecoLU<t_Value>::forwardElimination(
    valueType block[],
    integer   numRowsBlock,
    integer   numRowsPivot,
    integer   pivot[],
    valueType b[]
  ) const {
    valueType const * blockI = block;
    for ( integer i = 0; i < numRowsPivot; ++i, blockI += numRowsBlock ) {
      integer pivoti = pivot[i];
      if ( pivoti != i ) std::swap(b[pivoti],b[i]);
      valueType bi = b[i];
      for ( integer l = i+1; l < numRowsBlock; ++l ) b[l] -= blockI[l] * bi;
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
  ArcecoLU<t_Value>::forwardSolution(
    valueType block[],
    integer   numRowsBlock,
    integer   numColsPivot,
    integer   /* numOverlapCols */,
    valueType b[]
  ) const {
    integer     kk      = numRowsBlock - numColsPivot;
    valueType * blockJS = block + kk;
    for ( integer j = 0; j < numColsPivot; ++j, blockJS += numRowsBlock ) {
      valueType xj = (b[j] /= blockJS[j]);
      for ( integer l = j+1; l < numColsPivot; ++l )
        b[l] -= blockJS[l] * xj;
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
  ArcecoLU<t_Value>::forwardModification(
    valueType block[],
    integer   numRowsBlock,
    integer   numColsPivot,
    valueType b[]
  ) const {
    valueType * blockJ = block;
    for ( integer j = 0; j < numColsPivot; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j];
      for ( integer l = 0; l < numRowsBlock; ++l )
        b[numColsPivot + l] -= blockJ[l] * xj;
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
  ArcecoLU<t_Value>::backwardModification(
    valueType block[],
    integer   numRowsBlock,
    integer   numColsBlock,
    integer   numRowsPivot,
    valueType b[]
  ) const {
    valueType * blockJ = block + numRowsPivot * numRowsBlock;
    for ( integer j = numRowsPivot; j < numColsBlock; ++j, blockJ += numRowsBlock ) {
      valueType xj = b[j];
      for ( integer l = 0; l < numRowsPivot; ++l )
        b[l] -= blockJ[l] * xj;
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
  ArcecoLU<t_Value>::backwardSolution(
    valueType block[],
    integer   numRowsBlock,
    integer   /* numColsBlock */,
    integer   numRowsPivot,
    valueType b[]
  ) const {
    for ( integer j = numRowsPivot - 1; j >= 0; --j ) {
      valueType * blockJ = block + j * numRowsBlock;
      valueType xj = (b[j] /= blockJ[j]);
      for ( integer l = 0; l < j; ++l )
        b[l] -= blockJ[l] * xj;
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
  ArcecoLU<t_Value>::backwardElimination(
    valueType block[],
    integer   numRowsBlock,
    integer   numColsPivot,
    integer   numOverlapCols,
    integer   pivot[],
    valueType x[]
  ) const {
    integer kk = numRowsBlock - numColsPivot;
    integer j1 = numColsPivot;
    while ( j1 > 0 ) {
      valueType const * blockS = block + j1-1 + kk + j1 * numRowsBlock;
      integer n = numOverlapCols - j1;
      // Kahan summation algorithm
      valueType dotprd = dot(n,x+j1,1,blockS,numRowsBlock);
      x[--j1] -= dotprd;
      integer pivotj = pivot[j1];
      if ( pivotj != j1 ) std::swap( x[pivotj], x[j1] );
    }
  }

  template class ArcecoLU<float>;
  template class ArcecoLU<double>;

}
