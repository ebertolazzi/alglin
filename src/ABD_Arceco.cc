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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"
#include "Alglin_Eigen.hh"

namespace alglin {

  using std::abs;
  using std::swap;

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
  ArcecoLU<t_Value>::load_by_ref(
    integer const number_of_blocks,
    integer       matrix_structure[],
    real_type     array[],
    integer       pivot[]
  ) {
    m_number_of_blocks = number_of_blocks;
    m_matrix_structure = matrix_structure;
    m_array            = array;
    m_pivot_array      = pivot;
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
  ArcecoLU<t_Value>::checkStructure( integer const neq ) const {
    integer const & nblocks{ m_number_of_blocks };

    UTILS_ASSERT(
      noverlap(nblocks-1) == 0,
      "Arceco::checkStructure: noverlap({}) = {} expected zero!\n",
      nblocks-1, noverlap(nblocks-1)
    );

    // check index
    for ( integer k{0}; k < nblocks; ++k ) {
      UTILS_ASSERT(
        ncols(k) >= 1,
        "ArcecoLU::checkStructure: ncols({}) = {} < 1\n", k, ncols(k)
      );
      UTILS_ASSERT(
        nrows(k) >= 1,
        "ArcecoLU::checkStructure: nrows({}) = {} < 1\n", k, nrows(k)
      );
      UTILS_ASSERT(
        noverlap(k) >= 0,
        "ArcecoLU::checkStructure: noverlap({}) = {} < 0\n", k, noverlap(k)
      );
      UTILS_ASSERT(
        ncols(k) >= noverlap(k),
        "ArcecoLU::checkStructure: ncols({}) = {} < noverlap({}) = {}\n",
        k, ncols(k), k, noverlap(k)
      );
    }

    // check ovelapping
    for ( integer k{1}; k < nblocks; ++k )
      UTILS_ASSERT(
        noverlap(k-1) + noverlap(k) <= ncols(k),
        "Arceco::checkStructure: at block {} three consecutive block overlap\n", k
      );

    // controlla che i blocchi atraversino la diagonale
    //             c                  c+numCol
    // r           ┌──────────────────┐
    //             │                  │
    // r+numRow    └──────────────────┘
    integer r{ nrows(0) };
    integer c{ ncols(0) - noverlap(0) };
    for ( integer k{1}; k < nblocks; ++k ) {
      UTILS_ASSERT(
        c <= r && c+ncols(k) >= r+nrows(k),
        "ArcecoLU::checkStructure: block n. {} do not cross the diagonal\n", k
      );
      r += nrows(k);
      c += ncols(k)-noverlap(k);
    }

    integer isum1{0};
    integer isum2{0};
    for ( integer k{0}; k < nblocks; ++k ) {
      isum1 += nrows(k);
      isum2 += ncols(k) - noverlap(k);
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
    integer   const row0,
    integer   const col0,
    real_type const block0[],
    integer         nblk,
    integer         n,
    real_type const blocks[],
    integer   const rowN,
    integer   const colN,
    real_type const blockN[]
  ) {

    m_number_of_blocks = nblk + 2;

    integer const size0        { row0 * col0 };
    integer const sizeN        { rowN * colN };
    integer const BLK_size     { 2 * n * n * nblk };
    integer const numEquations { nblk * n + row0 + rowN };

    // allocazione dinamica
    m_mem.reallocate( BLK_size + size0 + sizeN );
    m_mem_int.reallocate( 3*m_number_of_blocks + numEquations );

    m_array            = m_mem( BLK_size + size0 + sizeN );
    m_pivot_array      = m_mem_int( numEquations );
    m_matrix_structure = m_mem_int( 3*m_number_of_blocks );

    // Fill structures
    Copy_n( block0, size0,    m_array );
    Copy_n( blocks, BLK_size, m_array + size0 );
    Copy_n( blockN, sizeN,    m_array + size0 + BLK_size );

    integer * mtr = m_matrix_structure;
    *mtr++ = row0;
    *mtr++ = col0;
    *mtr++ = n;
    for ( integer i{0}; i < nblk; ++i ) {
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
    integer nrows_block   = nrows(0);
    integer ncols_block   = ncols(0);
    integer noverlap_cols = noverlap(0);
    integer nrows_pivot   = ncols_block - noverlap_cols;

    /*\
     |  ┌──────────────────┐
     |  │                  │
     |  │                  │
     |  └─────────┬────────┴─────────┐
     |            │                  │
     |            └──────────────────┘
    \*/

    row_elimination(
      m_array + index1,
      nrows_block, ncols_block, nrows_pivot,
      m_pivot_array + indpiv
    );

    for ( integer k{1}; k < m_number_of_blocks; ++k ) {
      indpiv += nrows_pivot;
      integer index2       = index1 + nrows_block * nrows_pivot;
      integer index3       = index2 + nrows_block * noverlap_cols;
      integer ncols_pivot  = nrows_block - nrows_pivot;
      integer nrows_block2 = nrows(k);

      column_elimination(
        m_array + index2, nrows_block,  noverlap_cols,
        m_array + index3, nrows_block2, ncols_pivot,
        m_pivot_array + indpiv
      );

      nrows_block   = nrows_block2;
      index1        = index3 + nrows_block * ncols_pivot;
      ncols_block   = ncols(k) - ncols_pivot;
      noverlap_cols = noverlap(k);
      nrows_pivot   = ncols_block - noverlap_cols;
      indpiv        += ncols_pivot;

      row_elimination(
        m_array + index1,
        nrows_block, ncols_block, nrows_pivot,
        m_pivot_array + indpiv
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
  ArcecoLU<t_Value>::row_elimination(
    real_type     block[],
    integer const nrows_block,
    integer const ncols_block,
    integer const nrows_pivot,
    integer       pivot[]
  ) {

    #define BLOCK(I,J) block[(I) + (J) * nrows_block]

    for ( integer j{0}; j < nrows_pivot; ++j ) {
      integer const jplus1 { j + 1 };
      integer       jmax   { j };
      real_type     rowmax { abs(BLOCK(j,j)) };
      for ( integer i1{jplus1}; i1 < nrows_block; ++i1 ) {
        real_type tempiv{ abs(BLOCK(i1,j)) };
        if ( tempiv > rowmax ) { rowmax = tempiv; jmax = i1; }
      }

      UTILS_ASSERT0( rowmax > 0, "Arceco::row_elimination, singular matrix\n" );

      pivot[j] = jmax;
      if ( j != jmax )
        for ( integer j1{j}; j1 < ncols_block; ++j1 )
          std::swap( BLOCK(jmax,j1), BLOCK(j,j1) );

      real_type rowpiv{ BLOCK(j,j) };
      for ( integer i1{jplus1}; i1 < nrows_block; ++i1 ) {
        real_type rowmlt{ BLOCK(i1,j) /= rowpiv };
        for ( integer j1{jplus1}; j1 < ncols_block; ++j1 )
          BLOCK(i1,j1) -= rowmlt * BLOCK(j,j1);
      }
    }

    #undef BLOCK

  }

  template <typename t_Value>
  void
  ArcecoLU<t_Value>::column_elimination(
    real_type     topblk[],
    integer const nrows_top_block,
    integer const noverlap_cols,
    real_type     botblk[],
    integer const nrows_bottom_block,
    integer const ncols_pivot,
    integer       pivot[]
  ) {

    #define TOPBLK(I,J) topblk[(I) + (J) * nrows_top_block]
    #define BOTBLK(I,J) botblk[(I) + (J) * nrows_bottom_block]

    for ( integer j{0}; j < ncols_pivot; ++j ) {
      integer const jplus1 { j + 1 };
      integer const i      { nrows_top_block - ncols_pivot + j };
      integer       jmax   { j };
      real_type colmax{ abs(TOPBLK(i,j)) };
      for ( integer j1{jplus1}; j1 < noverlap_cols; ++j1 ) {
        real_type tempiv{ abs(TOPBLK(i,j1)) };
        if ( tempiv > colmax) { colmax = tempiv; jmax = j1; }
      }

      UTILS_ASSERT0( colmax > 0, "Arceco::column_elimination, singular matrix\n" );

      pivot[j] = jmax;
      if ( j != jmax ) {
        for ( integer k{i}; k < nrows_top_block;    ++k ) std::swap(TOPBLK(k,j),TOPBLK(k,jmax));
        for ( integer k{0}; k < nrows_bottom_block; ++k ) std::swap(BOTBLK(k,j),BOTBLK(k,jmax));
      }
      real_type colpiv{ TOPBLK(i,j) };
      for ( integer j1{jplus1}; j1 < noverlap_cols; ++j1 ) {
        real_type colmlt{ TOPBLK(i,j1) /= colpiv };
        for ( integer k{i+1}; k < nrows_top_block;    ++k ) TOPBLK(k,j1) -= colmlt * TOPBLK(k,j);
        for ( integer k{0};   k < nrows_bottom_block; ++k ) BOTBLK(k,j1) -= colmlt * BOTBLK(k,j);
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
  ArcecoLU<t_Value>::solve( real_type b[] ) const {
    integer indpiv        { 0 };
    integer indexa        { 0 };
    integer nrows_block   { nrows(0) };
    integer ncols_block   { ncols(0) };
    integer noverlap_cols { noverlap(0) };
    integer nrows_pivot   { ncols_block - noverlap_cols };

    forward_elimination(
      m_array + indexa,
      nrows_block, nrows_pivot,
      m_pivot_array + indpiv, b + indpiv
    );

    integer ncols_pivot{0};
    for ( integer k{1}; k < m_number_of_blocks; ++k ) {
      indexa     += nrows_block * nrows_pivot;
      ncols_pivot = nrows_block - nrows_pivot;
      indpiv     += nrows_pivot;

      //if ( ncols_pivot > 0 )
      forward_solution(
        m_array + indexa,
        nrows_block, ncols_pivot, noverlap_cols,
        b + indpiv
      );

      indexa     += noverlap_cols * nrows_block;
      nrows_block = nrows(k);

      //if ( ncols_pivot > 0 )
      forward_modification(
        m_array + indexa,
        nrows_block, ncols_pivot,
        b + indpiv
      );

      indexa       += nrows_block * ncols_pivot;
      ncols_block   = ncols(k) - ncols_pivot;
      noverlap_cols = noverlap(k);
      nrows_pivot   = ncols_block - noverlap_cols;
      indpiv       += ncols_pivot;

      //if ( nrows_pivot > 0 )
      forward_elimination(
        m_array + indexa,
        nrows_block, nrows_pivot,
        m_pivot_array + indpiv, b + indpiv
      );
    }
    // BACKWARD LOOP
    for ( integer k{m_number_of_blocks - 2}; k >= 0; --k ) {

      if ( nrows_pivot != 0 ) {
        if ( nrows_pivot != ncols_block )
          backward_modification(
            m_array + indexa,
            nrows_block,
            ncols_block,
            nrows_pivot,
            b + indpiv
          );
        backward_solution(
          m_array + indexa,
          nrows_block,
          ncols_block,
          nrows_pivot,
          b + indpiv
        );
      }

      indexa       -= nrows_block * ncols_pivot;
      nrows_block   = nrows(k);
      noverlap_cols = noverlap(k);
      indexa       -= nrows_block * noverlap_cols;
      indpiv       -= ncols_pivot;

      //if ( ncols_pivot > 0 )
      backward_elimination(
        m_array + indexa,
        nrows_block, ncols_pivot, noverlap_cols,
        m_pivot_array + indpiv, b + indpiv
      );

      nrows_pivot = nrows_block - ncols_pivot;
      ncols_block = noverlap_cols + nrows_pivot;
      indexa     -= nrows_block * nrows_pivot;
      indpiv     -= nrows_pivot;
      ncols_pivot = ncols(k) - ncols_block;
    }

    if ( nrows_pivot > 0 ) {
      if ( nrows_pivot != ncols_block )
        backward_modification(
          m_array + indexa,
          nrows_block, ncols_block, nrows_pivot,
          b + indpiv
        );
      backward_solution(
        m_array + indexa,
        nrows_block, ncols_block, nrows_pivot,
        b + indpiv
      );
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
  ArcecoLU<t_Value>::forward_elimination(
    real_type const block[],
    integer   const nrows_block,
    integer   const nrows_pivot,
    integer   const pivot[],
    real_type       b[]
  ) const {
    real_type const * blockI{block};
    for ( integer i{0}; i < nrows_pivot; ++i, blockI += nrows_block ) {
      integer pivoti{pivot[i]};
      if ( pivoti != i ) std::swap( b[pivoti], b[i] );
      real_type bi{b[i]};
      for ( integer l{i+1}; l < nrows_block; ++l ) b[l] -= blockI[l] * bi;
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
  ArcecoLU<t_Value>::forward_solution(
    real_type     block[],
    integer const nrows_block,
    integer const ncols_pivot,
    integer   /* noverlap_cols */,
    real_type     b[]
  ) const {
    integer     kk      { nrows_block - ncols_pivot };
    real_type * blockJS { block + kk };
    for ( integer j{0}; j < ncols_pivot; ++j, blockJS += nrows_block ) {
      real_type xj{ b[j] /= blockJS[j] };
      for ( integer l{j+1}; l < ncols_pivot; ++l )
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
  ArcecoLU<t_Value>::forward_modification(
    real_type     block[],
    integer const nrows_block,
    integer const ncols_pivot,
    real_type     b[]
  ) const {
    real_type * blockJ{ block };
    for ( integer j{0}; j < ncols_pivot; ++j, blockJ += nrows_block ) {
      real_type xj{ b[j] };
      for ( integer l{0}; l < nrows_block; ++l )
        b[ncols_pivot + l] -= blockJ[l] * xj;
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
  ArcecoLU<t_Value>::backward_modification(
    real_type     block[],
    integer const nrows_block,
    integer const ncols_block,
    integer const nrows_pivot,
    real_type     b[]
  ) const {
    real_type * blockJ{ block + nrows_pivot * nrows_block };
    for ( integer j{nrows_pivot}; j < ncols_block; ++j, blockJ += nrows_block ) {
      real_type xj{b[j]};
      for ( integer l{0}; l < nrows_pivot; ++l )
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
  ArcecoLU<t_Value>::backward_solution(
    real_type     block[],
    integer const nrows_block,
    integer       /* ncols_block */,
    integer const nrows_pivot,
    real_type     b[]
  ) const {
    for ( integer j{nrows_pivot - 1}; j >= 0; --j ) {
      real_type * blockJ{ block + j * nrows_block };
      real_type xj{ b[j] /= blockJ[j] };
      for ( integer l{0}; l < j; ++l )
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
  ArcecoLU<t_Value>::backward_elimination(
    real_type     block[],
    integer const nrows_block,
    integer const ncols_pivot,
    integer const noverlap_cols,
    integer const pivot[],
    real_type     x[]
  ) const {
    integer kk{ nrows_block - ncols_pivot };
    integer j1{ ncols_pivot };
    while ( j1 > 0 ) {
      real_type const * blockS { block + j1-1 + kk + j1 * nrows_block };
      integer   const   n      { noverlap_cols - j1 };
      // Kahan summation algorithm
      real_type const dotprd { dot(n,x+j1,1,blockS,nrows_block) };
      x[--j1] -= dotprd;
      integer pivotj{ pivot[j1] };
      if ( pivotj != j1 ) swap( x[pivotj], x[j1] );
    }
  }

  template class ArcecoLU<float>;
  template class ArcecoLU<double>;

}
