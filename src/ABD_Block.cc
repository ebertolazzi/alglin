/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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

///
/// file: ABD_Block.cc
///

#include "ABD_Block.hh"
#include <iomanip>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace alglin {

  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/
  template <typename t_Value>
  void
  BlockLU<t_Value>::allocateTopBottom( integer _nblock,
                                       integer _n,
                                       integer _row0,
                                       integer _col0,
                                       integer _rowN,
                                       integer _colN,
                                       integer _nb ) {
    F_lda    = _n + _row0 - _col0 ;
    F_size   = _n*F_lda ;
    Work_lda = F_lda + _n ;
    integer Fnnz = _nblock*F_size ;
    integer innz = _nblock*_n ;
    BlockBidiagonal<t_Value>::allocateTopBottom( _nblock, _n,
                                                 _row0, _col0,
                                                 _rowN, _colN,
                                                 _nb,
                                                 Fnnz+2*Work_lda*_n,
                                                 innz+_row0) ;
    swap0      = this->baseInteger(size_t(_row0)) ;
    swapR_blks = this->baseInteger(size_t(innz)) ;
    F_mat      = this->baseValue(size_t(Fnnz)) ;
    Work_mat   = this->baseValue(size_t(Work_lda*_n)) ;
    Work_mat1  = this->baseValue(size_t(Work_lda*_n)) ;
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
  BlockLU<t_Value>::factorize() {

    ALGLIN_ASSERT( this->numCyclicOMEGA == 0 && this->numCyclicBC == 0,
                   "BlockLU cannot manage cyclic BC" ) ;

    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;

    integer const & col00  = this->numInitialOMEGA ;
    integer const & colNN  = this->numFinalOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valuePointer & block0 = this->block0 ;
    valuePointer & blockN = this->blockN ;
    
    if ( nblock > 0 ) {

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U | M       |
      //  |     |         |
      //  +-----+---------+
      */
      t_Value * block00 = block0 + col00*row0 ;
      if ( col00 > 0 ) {
        integer info = getrf( row0, col00, block0, row0, swap0 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, getrf INFO = " << info ) ;
        // applico permutazioni
        info = swaps( n, block00, row0, 0, col00-1, swap0, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, swaps INFO = " << info ) ;
        trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, col00, n, 1.0, block0, row0, block00, row0 ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              row00, n, col00,
              -1.0, block0+col00,  row0,
                    block00,       row0,
               1.0, block00+col00, row0 ) ;
      }

      /*
      // +------------+
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |      A     |     B      |
      // |            |            |
      // |            |            |
      // +------------+------------+------------+
      //              |            |            |
      //              |            |            |
      //
      */

      // copy from
      integer info = gecopy( row00, n, block00 + col00, row0, Work_mat, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;
      gezero( row00, n, Work_mat1, Work_lda ) ;

      info = gecopy( n, n, this->AdAu_blk, n, Work_mat+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      info = gecopy( n, n, this->AdAu_blk + nxn, n, Work_mat1+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      // factorize
      info = getrf( n+row00, n, Work_mat, Work_lda, swapR_blks ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, getrf INFO = " << info ) ;

      // apply permutation
      info = swaps( n, Work_mat1, Work_lda, 0, n-1, swapR_blks, 1 ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, swaps INFO = " << info ) ;
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, n, 1.0, Work_mat, Work_lda, Work_mat1, Work_lda ) ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE,
            row00, n, n,
            -1.0, Work_mat+n,  Work_lda,
                  Work_mat1,   Work_lda,
             1.0, Work_mat1+n, Work_lda ) ;

      // copy to
      info = gecopy( row00, n, Work_mat, Work_lda, block00 + col00, row0 ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      info = gecopy( row00, n, Work_mat1, Work_lda, F_mat, F_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      info = gecopy( n, n, Work_mat+row00, Work_lda, this->AdAu_blk, n ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      info = gecopy( n, n, Work_mat1+row00, Work_lda, this->AdAu_blk + nxn, n ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

      /*
      // |    A       |
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |            |            |
      // |     B      |    C       |
      // |            |            |
      // +------------+------------+------------+
      //              |            |            |
      //              |            |            |
      //
      */

      for ( integer k = 1 ; k < nblock ; ++k ) {

        valuePointer B = this->AdAu_blk + nxnx2 * k ;
        valuePointer A = B - nxn ;
        valuePointer C = B + nxn ;
        valuePointer F = F_mat + k*F_size ;
        integer * swps = swapR_blks + k * n ;

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, Work_mat, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;
        gezero( row00, n, Work_mat1, Work_lda ) ;

        info = gecopy( n, n, B, n, Work_mat+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n, n, C, n, Work_mat1+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        // factorize
        info = getrf( n+row00, n, Work_mat, Work_lda, swps ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, getrf INFO = " << info ) ;

        // apply permutation
        info = swaps( n, Work_mat1, Work_lda, 0, n-1, swps, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, swaps INFO = " << info ) ;
        trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, n, 1.0, Work_mat, Work_lda, Work_mat1, Work_lda ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              row00, n, n,
              -1.0, Work_mat+n,  Work_lda,
                    Work_mat1,   Work_lda,
               1.0, Work_mat1+n, Work_lda ) ;

        // copy to
        info = gecopy( row00, n, Work_mat, Work_lda, A + n_m_row00, n ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( row00, n, Work_mat1, Work_lda, F, F_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n, n, Work_mat+row00, Work_lda, B, n ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n, n, Work_mat1+row00, Work_lda, C, n ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;
      }
      /*
      // |    A       |      B     |
      // |            | ********** |
      // +------------+------------+-+
      //              |              |
      //              |   BlockN     |
      //              |              |
      //              +--------------+
      //
      */

      // fattorizzazione ultimo blocco
      valuePointer B = this->AdAu_blk + nxnx2 * nblock - nxn ;
      integer N = row00+rowN ;
      this->la_factorization->allocate(N,N) ;
      this->la_factorization->zero_block(row00,N-n,0,n) ;
      this->la_factorization->load_block(row00,n,B+n_m_row00,n,0,0) ;
      this->la_factorization->load_block(rowN,N,blockN,rowN,row00,0) ;
      this->la_factorization->factorize() ;
    } else { // case nblock == 0
      integer N = row0+rowN ;
      this->la_factorization->allocate(N,N) ;
      this->la_factorization->zero_block(N,N,0,0) ;
      this->la_factorization->load_block(row0,n+col00,block0,row0,0,0) ;
      this->la_factorization->load_block(rowN,n+colNN,blockN,rowN,row0,col00) ;
      this->la_factorization->factorize() ;
    }
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
  BlockLU<t_Value>::solve( valuePointer in_out ) const {

    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;
    integer const & col00  = this->numInitialOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valueConstPointer const & block0 = this->block0 ;

    // permuto le x
    integer neq = nblock*n+row0+rowN ;
    std::rotate( in_out, in_out + neq - row0, in_out + neq ) ;

    if ( nblock > 0 ) {

      valuePointer io = in_out ;
      /*
      //             _             _
      //   ___  ___ | |_   _____  | |
      //  / __|/ _ \| \ \ / / _ \ | |
      //  \__ \ (_) | |\ V /  __/ | |___
      //  |___/\___/|_| \_/ \___| |_____|
      */
      /*
      // primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        // apply permutation
        integer info = swaps( 1, io, row0, 0, col00-1, swap0, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;
        trsv( LOWER, NO_TRANSPOSE, UNIT, col00, block0, row0, in_out, 1 ) ;
        io += col00 ;
        gemv( NO_TRANSPOSE,
              row00, col00,
              -1.0, block0+col00, row0,
                    in_out, 1,
               1.0, io, 1 ) ;
      }

      /*
      // +-+------------+
      // | |            | * * * * * * <-- fill in
      // +-+------------+------------+
      //   |            |            |
      //   |            |            |
      //   |     A      |     B      |
      //   |            |            |
      //   +------------+------------+------------+
      //                |            |            |
      //                |            |            |
      //
      */
      valueConstPointer block00 = block0 + col00*row0 ;

      // apply permutation
      integer info = swaps( 1, io, n, 0, n-1, swapR_blks, 1 ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, Work_mat, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, gecopy INFO = " << info ) ;

      info = gecopy( n, n, this->AdAu_blk, n, Work_mat+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, gecopy INFO = " << info ) ;

      trsv( LOWER, NO_TRANSPOSE, UNIT, n, Work_mat, Work_lda, io, 1 ) ;
      gemv( NO_TRANSPOSE,
            row00, n,
            -1.0, Work_mat+n, Work_lda,
                  io, 1,
             1.0, io+n, 1 ) ;
      io += n ;

      /*
      // |            |
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |            |            |
      // |            |    A       |
      // |            |            |
      // +------------+------------+------------+
      //              |            |            |
      //              |    B       |            |
      //
      */
      integer k = 0 ;
      while ( ++k < nblock ) {

        valuePointer B = this->AdAu_blk + nxnx2 * k ;
        valuePointer A = B - nxn ;
        integer * swps = swapR_blks + k * n ;

        // apply permutation
        info = swaps( 1, io, n, 0, n-1, swps, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, Work_mat, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n, n, B, n, Work_mat+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        trsv( LOWER, NO_TRANSPOSE, UNIT, n, Work_mat, Work_lda, io, 1 ) ;
        gemv( NO_TRANSPOSE,
              row00, n,
              -1.0, Work_mat+n, Work_lda,
                    io, 1,
              1.0, io+n, 1 ) ;
        io += n ;
      }

      /*
      // |    A       |      B     |
      // |            |            |
      // +------------+------------+-+
      //              |              |
      //              |   BlockN     |
      //              |              |
      //              +--------------+
      //
      */

      // risolvo ultimo blocco
      this->la_factorization->solve( io ) ;
    
      /*
      //             _             _   _
      //   ___  ___ | |_   _____  | | | |
      //  / __|/ _ \| \ \ / / _ \ | | | |
      //  \__ \ (_) | |\ V /  __/ | |_| |
      //  |___/\___/|_| \_/ \___|  \___/
      */
      /*
      // |     A      |
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |            |            |
      // |     B      |     C      |
      // |            |            | * * * * * * <- fill in
      // +------------+------------+------------+
      //              |            |            |
      //              |            |            |
      //
      */
      while ( --k > 0 ) {

        valuePointer B = this->AdAu_blk + nxnx2 * k ;
        valuePointer C = B + nxn ;
        valuePointer A = B - nxn ;
        valuePointer F = F_mat + k*F_size ;
      
        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, Work_mat, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

        info = gecopy( row00, n, F, F_lda, Work_mat1, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

        info = gecopy( n_m_row00, n, B, n, Work_mat+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n_m_row00, n, C, n, Work_mat1+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        io -= n ;
        gemv( NO_TRANSPOSE,
              n, n,
              -1.0, Work_mat1, Work_lda,
                    io+n,      1,
               1.0, io,        1 ) ;

        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, Work_mat, Work_lda, io, 1 ) ;

      }

      valuePointer B = this->AdAu_blk ;
      valuePointer C = B + nxn ;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, Work_mat, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

      info = gecopy( row00, n, F_mat, F_lda, Work_mat1, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

      info = gecopy( n_m_row00, n, B, n, Work_mat+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = " << info ) ;

      info = gecopy( n_m_row00, n, C, n, Work_mat1+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = " << info ) ;

      io -= n ;
      gemv( NO_TRANSPOSE,
            n, n,
            -1.0, Work_mat1, Work_lda,
                  io+n,      1,
             1.0, io,        1 ) ;

      trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, Work_mat, Work_lda, io, 1 ) ;

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        io -= col00 ;
        gemv( NO_TRANSPOSE,
              col00, n,
              -1.0, block00,  row0,
                    io+col00, 1,
               1.0, io,       1 ) ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, col00, block0, row0, io, 1 ) ;
      }
    } else {  // case nblock = 0
      this->la_factorization->solve( in_out ) ;
    }

    // permuto le x
    std::rotate( in_out, in_out + col00, in_out + neq ) ;
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
  BlockLU<t_Value>::solve( integer      nrhs,
                           valuePointer in_out,
                           integer      ldRhs ) const {

    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;
    integer const & col00  = this->numInitialOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valueConstPointer const & block0 = this->block0 ;

    integer neq = nblock*n+row0+rowN ;

    // permuto le x
    valuePointer io = in_out ;
    for ( integer k = 0 ; k < nrhs ; ++k ) {
      std::rotate( io, io + neq - row0, io + neq ) ;
      io += ldRhs ;
    }


    if ( nblock > 0 ) {

      io = in_out ;
      /*
      //             _             _
      //   ___  ___ | |_   _____  | |
      //  / __|/ _ \| \ \ / / _ \ | |
      //  \__ \ (_) | |\ V /  __/ | |___
      //  |___/\___/|_| \_/ \___| |_____|
      */
      /*
      // primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        // apply permutation
        integer info = swaps( 1, io, row0, 0, col00-1, swap0, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;
        trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, col00, nrhs, 1.0, block0, row0, in_out, ldRhs ) ;
        io += col00 ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              row00, nrhs, col00,
              -1.0, block0+col00, row0,
                    in_out, ldRhs,
               1.0, io, ldRhs ) ;
      }

      /*
      // +-+------------+
      // | |            | * * * * * * <-- fill in
      // +-+------------+------------+
      //   |            |            |
      //   |            |            |
      //   |     A      |     B      |
      //   |            |            |
      //   +------------+------------+------------+
      //                |            |            |
      //                |            |            |
      //
      */
      valueConstPointer block00 = block0 + col00*row0 ;

      // apply permutation
      integer info = swaps( ldRhs, io, ldRhs, 0, n-1, swapR_blks, 1 ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, Work_mat, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, gecopy INFO = " << info ) ;

      info = gecopy( n, n, this->AdAu_blk, n, Work_mat+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, gecopy INFO = " << info ) ;

      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, nrhs, 1.0, Work_mat, Work_lda, io, ldRhs ) ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE,
            row00, nrhs, n,
            -1.0, Work_mat+n, Work_lda,
                  io, ldRhs,
             1.0, io+n, ldRhs ) ;
      io += n ;

      /*
      // |            |
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |            |            |
      // |            |    A       |
      // |            |            |
      // +------------+------------+------------+
      //              |            |            |
      //              |    B       |            |
      //
      */
      integer k = 0 ;
      while ( ++k < nblock ) {

        valuePointer B = this->AdAu_blk + nxnx2 * k ;
        valuePointer A = B - nxn ;
        integer * swps = swapR_blks + k * n ;

        // apply permutation
        info = swaps( nrhs, io, ldRhs, 0, n-1, swps, 1 ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = " << info ) ;

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, Work_mat, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n, n, B, n, Work_mat+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, nrhs, 1.0, Work_mat, Work_lda, io, ldRhs ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              row00, nrhs, n,
              -1.0, Work_mat+n, Work_lda,
                    io, ldRhs,
               1.0, io+n, ldRhs ) ;
        io += n ;
      }

      /*
      // |    A       |      B     |
      // |            |            |
      // +------------+------------+-+
      //              |              |
      //              |   BlockN     |
      //              |              |
      //              +--------------+
      //
      */

      // risolvo ultimo blocco
      this->la_factorization->solve( nrhs, io, ldRhs ) ;
    
      /*
      //             _             _   _
      //   ___  ___ | |_   _____  | | | |
      //  / __|/ _ \| \ \ / / _ \ | | | |
      //  \__ \ (_) | |\ V /  __/ | |_| |
      //  |___/\___/|_| \_/ \___|  \___/
      */
      /*
      // |     A      |
      // |            | * * * * * * <-- fill in
      // +------------+------------+
      // |            |            |
      // |            |            |
      // |     B      |     C      |
      // |            |            | * * * * * * <- fill in
      // +------------+------------+------------+
      //              |            |            |
      //              |            |            |
      //
      */
      while ( --k > 0 ) {

        valuePointer B = this->AdAu_blk + nxnx2 * k ;
        valuePointer C = B + nxn ;
        valuePointer A = B - nxn ;
        valuePointer F = F_mat + k*F_size ;
      
        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, Work_mat, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

        info = gecopy( row00, n, F, F_lda, Work_mat1, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

        info = gecopy( n_m_row00, n, B, n, Work_mat+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        info = gecopy( n_m_row00, n, C, n, Work_mat1+row00, Work_lda ) ;
        ALGLIN_ASSERT( info == 0, "BlockLU::factorize, first block, gecopy INFO = " << info ) ;

        io -= n ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, nrhs, n,
              -1.0, Work_mat1, Work_lda,
                    io+n,      ldRhs,
               1.0, io,        ldRhs ) ;

        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, n, nrhs, 1.0, Work_mat, Work_lda, io, ldRhs ) ;

      }

      valuePointer B = this->AdAu_blk ;
      valuePointer C = B + nxn ;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, Work_mat, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

      info = gecopy( row00, n, F_mat, F_lda, Work_mat1, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = " << info ) ;

      info = gecopy( n_m_row00, n, B, n, Work_mat+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = " << info ) ;

      info = gecopy( n_m_row00, n, C, n, Work_mat1+row00, Work_lda ) ;
      ALGLIN_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = " << info ) ;

      io -= n ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE,
            n, nrhs, n,
            -1.0, Work_mat1, Work_lda,
                  io+n,      ldRhs,
             1.0, io,        ldRhs ) ;

      trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, n, nrhs, 1.0, Work_mat, Work_lda, io, ldRhs ) ;

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        io -= col00 ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              col00, nrhs, n,
              -1.0, block00,  row0,
                    io+col00, ldRhs,
               1.0, io,       ldRhs ) ;
        trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, col00, nrhs, 1.0, block0, row0, io, ldRhs ) ;
      }
    } else {  // case nblock = 0
      this->la_factorization->solve( nrhs, in_out, ldRhs ) ;
    }

    // permuto le x
    io = in_out ;
    for ( integer k = 0 ; k < nrhs ; ++k ) {
      std::rotate( io, io + col00, io + neq ) ;
      io += ldRhs ;
    }
  }

  template class BlockLU<double> ;
  template class BlockLU<float> ;

}

///
/// eof: ABD_Block.cc
///

