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

///
/// file: ABD_Block.cc
///

#include "Alglin.hh"

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
  BlockLU<t_Value>::allocateTopBottom(
    integer _nblock,
    integer _n,
    integer _row0,
    integer _col0,
    integer _rowN,
    integer _colN,
    integer _nb
  ) {
    m_F_lda    = _n + _row0 - _col0;
    m_F_size   = _n*m_F_lda;
    m_Work_lda = m_F_lda + _n;
    integer Fnnz = _nblock*m_F_size;
    integer innz = _nblock*_n;
    BlockBidiagonal<t_Value>::allocateTopBottom(
      _nblock, _n,
      _row0, _col0,
      _rowN, _colN,
      _nb,
      Fnnz+2*m_Work_lda*_n,
      innz+_row0
    );
    m_swap0      = m_baseInteger( size_t(_row0) );
    m_swapR_blks = m_baseInteger( size_t(innz) );
    m_F_mat      = m_baseValue( size_t(Fnnz) );
    m_Work_mat   = m_baseValue( size_t(m_Work_lda*_n) );
    m_Work_mat1  = m_baseValue( size_t(m_Work_lda*_n) );
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

    UTILS_ASSERT0(
      this->numCyclicOMEGA == 0 && this->numCyclicBC == 0,
      "BlockLU cannot manage cyclic BC\n"
    );

    integer const & n      = this->n;
    integer const & nxnx2  = this->nxnx2;
    integer const & nxn    = this->nxn;
    integer const & nblock = this->nblock;

    integer const & col00  = this->numInitialOMEGA;
    integer const & colNN  = this->numFinalOMEGA;
    integer const & row0   = this->numInitialBC;
    integer const & rowN   = this->numFinalBC;

    integer const row00     = row0 - col00;
    integer const n_m_row00 = n - row00;
    
    if ( nblock > 0 ) {

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U | M       |
      //  |     |         |
      //  +-----+---------+
      */
      t_Value * block00 = m_block0 + col00*row0;
      if ( col00 > 0 ) {
        integer info = getrf( row0, col00, m_block0, row0, m_swap0 );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, getrf INFO = {}\n", info
        );
        // applico permutazioni
        info = swaps( n, block00, row0, 0, col00-1, m_swap0, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, swaps INFO = {}\n", info
        );
        trsm(
          LEFT, LOWER, NO_TRANSPOSE, UNIT,
          col00, n, 1.0, m_block0, row0, block00, row0
        );
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          row00, n, col00,
          -1.0, m_block0+col00, row0,
                block00,        row0,
           1.0, block00+col00,  row0
        );
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
      integer info = gecopy( row00, n, block00 + col00, row0, m_Work_mat, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );
      gezero( row00, n, m_Work_mat1, m_Work_lda );

      info = gecopy( n, n, m_DE_blk, n, m_Work_mat+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( n, n, m_DE_blk + nxn, n, m_Work_mat1+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0,
        "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

      // factorize
      info = getrf( n+row00, n, m_Work_mat, m_Work_lda, m_swapR_blks );
      UTILS_ASSERT(
        info == 0,
        "BlockLU::factorize, first block, getrf INFO = {}\n", info
      );

      // apply permutation
      info = swaps( n, m_Work_mat1, m_Work_lda, 0, n-1, m_swapR_blks, 1 );
      UTILS_ASSERT(
        info == 0,
        "BlockLU::factorize, first block, swaps INFO = {}\n", info
      );
      trsm(
        LEFT, LOWER, NO_TRANSPOSE, UNIT,
        n, n, 1.0, m_Work_mat, m_Work_lda, m_Work_mat1, m_Work_lda
      );
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        row00, n, n,
        -1.0, m_Work_mat+n,  m_Work_lda,
              m_Work_mat1,   m_Work_lda,
         1.0, m_Work_mat1+n, m_Work_lda
      );

      // copy to
      info = gecopy( row00, n, m_Work_mat, m_Work_lda, block00 + col00, row0 );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( row00, n, m_Work_mat1, m_Work_lda, m_F_mat, m_F_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( n, n, m_Work_mat+row00, m_Work_lda, m_DE_blk, n );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( n, n, m_Work_mat1+row00, m_Work_lda, m_DE_blk + nxn, n );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
      );

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

      for ( integer k = 1; k < nblock; ++k ) {

        valueType * B    = m_DE_blk + nxnx2 * k;
        valueType * A    = B - nxn;
        valueType * C    = B + nxn;
        valueType * F    = m_F_mat + k*m_F_size;
        integer   * swps = m_swapR_blks + k * n;

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, m_Work_mat, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );
        gezero( row00, n, m_Work_mat1, m_Work_lda );

        info = gecopy( n, n, B, n, m_Work_mat+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n, n, C, n, m_Work_mat1+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        // factorize
        info = getrf( n+row00, n, m_Work_mat, m_Work_lda, swps );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, getrf INFO = {}\n", info
        );

        // apply permutation
        info = swaps( n, m_Work_mat1, m_Work_lda, 0, n-1, swps, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, swaps INFO = {}\n", info
        );
        trsm(
          LEFT, LOWER, NO_TRANSPOSE, UNIT,
          n, n, 1.0, m_Work_mat, m_Work_lda, m_Work_mat1, m_Work_lda
        );
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          row00, n, n,
          -1.0, m_Work_mat+n,  m_Work_lda,
                m_Work_mat1,   m_Work_lda,
           1.0, m_Work_mat1+n, m_Work_lda
        );

        // copy to
        info = gecopy( row00, n, m_Work_mat, m_Work_lda, A + n_m_row00, n );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( row00, n, m_Work_mat1, m_Work_lda, F, m_F_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n, n, m_Work_mat+row00, m_Work_lda, B, n );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n, n, m_Work_mat1+row00, m_Work_lda, C, n );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );
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
      valueType * B = m_DE_blk + nxnx2 * nblock - nxn;
      integer N = row00+rowN;
      this->la_matrix.setup( N, N );
      this->la_matrix.zero_block(row00,N-n,0,n);
      this->la_matrix.load_block(row00,n,B+n_m_row00,n,0,0);
      this->la_matrix.load_block(rowN,N,m_blockN,rowN,row00,0);
      this->la_factorization->factorize( "BlockLU::factorize", this->la_matrix );
    } else { // case nblock == 0
      integer N = row0+rowN;
      this->la_matrix.setup( N, N );
      this->la_matrix.zero_block(N,N,0,0);
      this->la_matrix.load_block(row0,n+col00,m_block0,row0,0,0);
      this->la_matrix.load_block(rowN,n+colNN,m_blockN,rowN,row0,col00);
      this->la_factorization->factorize( "BlockLU::factorize[nblock=0]", this->la_matrix );
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
  BlockLU<t_Value>::solve_internal(
    bool      do_permute,
    valueType in_out[]
  ) const {

    integer const & n      = this->n;
    integer const & nxnx2  = this->nxnx2;
    integer const & nxn    = this->nxn;
    integer const & nblock = this->nblock;
    integer const & col00  = this->numInitialOMEGA;
    integer const & row0   = this->numInitialBC;
    integer const & rowN   = this->numFinalBC;

    integer const row00     = row0 - col00;
    integer const n_m_row00 = n - row00;

    // permuto le x
    integer neq = nblock*n+row0+rowN;
    if ( do_permute ) std::rotate( in_out, in_out + neq - row0, in_out + neq );

    if ( nblock > 0 ) {

      valueType * io = in_out;
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
        integer info = swaps( 1, io, row0, 0, col00-1, m_swap0, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
        );
        trsv( LOWER, NO_TRANSPOSE, UNIT, col00, m_block0, row0, in_out, 1 );
        io += col00;
        gemv(
          NO_TRANSPOSE,
          row00, col00,
          -1.0, m_block0+col00, row0,
                in_out, 1,
           1.0, io, 1
        );
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
      valueType const * block00 = m_block0 + col00*row0;

      // apply permutation
      integer info = swaps( 1, io, n, 0, n-1, m_swapR_blks, 1 );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
      );

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, m_Work_mat, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( n, n, m_DE_blk, n, m_Work_mat+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, gecopy INFO = {}\n", info
      );

      trsv( LOWER, NO_TRANSPOSE, UNIT, n, m_Work_mat, m_Work_lda, io, 1 );
      gemv(
        NO_TRANSPOSE,
        row00, n,
        -1.0, m_Work_mat+n, m_Work_lda,
              io, 1,
         1.0, io+n, 1
      );
      io += n;

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
      integer k = 0;
      while ( ++k < nblock ) {

        valueType * B    = m_DE_blk + nxnx2 * k;
        valueType * A    = B - nxn;
        integer   * swps = m_swapR_blks + k * n;

        // apply permutation
        info = swaps( 1, io, n, 0, n-1, swps, 1 );
        UTILS_ASSERT(
          info == 0,
          "BlockLU::solve, first block, swaps INFO = {}\n", info
        );

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, m_Work_mat, m_Work_lda );
        UTILS_ASSERT(
          info == 0,
          "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n, n, B, n, m_Work_mat+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0,
          "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        trsv( LOWER, NO_TRANSPOSE, UNIT, n, m_Work_mat, m_Work_lda, io, 1 );
        gemv(
          NO_TRANSPOSE,
          row00, n,
          -1.0, m_Work_mat+n, m_Work_lda,
                io, 1,
           1.0, io+n, 1
        );
        io += n;
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
      this->la_factorization->solve( io );
    
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

        valueType * B = m_DE_blk + nxnx2 * k;
        valueType * C = B + nxn;
        valueType * A = B - nxn;
        valueType * F = m_F_mat + k*m_F_size;
      
        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, m_Work_mat, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
        );

        info = gecopy( row00, n, F, m_F_lda, m_Work_mat1, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
        );

        info = gecopy( n_m_row00, n, B, n, m_Work_mat+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n_m_row00, n, C, n, m_Work_mat1+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        io -= n;
        gemv(
          NO_TRANSPOSE,
          n, n,
          -1.0, m_Work_mat1, m_Work_lda,
                io+n,        1,
           1.0, io,          1
        );

        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, m_Work_mat, m_Work_lda, io, 1 );

      }

      valueType * B = m_DE_blk;
      valueType * C = B + nxn;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, m_Work_mat, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
      );

      info = gecopy( row00, n, m_F_mat, m_F_lda, m_Work_mat1, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
      );

      info = gecopy( n_m_row00, n, B, n, m_Work_mat+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, U block, gecopy INFO = {}\n", info
      );

      info = gecopy( n_m_row00, n, C, n, m_Work_mat1+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::factorize, U block, gecopy INFO = {}\n", info
      );

      io -= n;
      gemv(
        NO_TRANSPOSE,
        n, n,
        -1.0, m_Work_mat1, m_Work_lda,
              io+n,        1,
         1.0, io,          1
      );

      trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, m_Work_mat, m_Work_lda, io, 1 );

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        io -= col00;
        gemv(
          NO_TRANSPOSE,
          col00, n,
          -1.0, block00,  row0,
                io+col00, 1,
           1.0, io,       1
        );
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, col00, m_block0, row0, io, 1 );
      }
    } else {  // case nblock = 0
      this->la_factorization->solve( in_out );
    }

    // permuto le x
    if ( do_permute ) std::rotate( in_out, in_out + col00, in_out + neq );
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
  BlockLU<t_Value>::solve_internal(
    bool      do_permute,
    integer   nrhs,
    valueType in_out[],
    integer   ldRhs
  ) const {

    integer const & n      = this->n;
    integer const & nxnx2  = this->nxnx2;
    integer const & nxn    = this->nxn;
    integer const & nblock = this->nblock;
    integer const & col00  = this->numInitialOMEGA;
    integer const & row0   = this->numInitialBC;
    integer const & rowN   = this->numFinalBC;

    integer const row00     = row0 - col00;
    integer const n_m_row00 = n - row00;

    integer neq = nblock*n+row0+rowN;

    // permuto le x
    valueType * io = in_out;
    if ( do_permute ) {
      for ( integer k = 0; k < nrhs; ++k ) {
        std::rotate( io, io + neq - row0, io + neq );
        io += ldRhs;
      }
      io = in_out;
    }

    if ( nblock > 0 ) {

      io = in_out;
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
        integer info = swaps( 1, io, row0, 0, col00-1, m_swap0, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
        );
        trsm(
          LEFT, LOWER, NO_TRANSPOSE, UNIT,
          col00, nrhs, 1.0, m_block0, row0, in_out, ldRhs
        );
        io += col00;
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          row00, nrhs, col00,
          -1.0, m_block0+col00, row0,
                in_out, ldRhs,
           1.0, io, ldRhs
        );
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
      valueType const * block00 = m_block0 + col00*row0;

      // apply permutation
      integer info = swaps( nrhs, io, ldRhs, 0, n-1, m_swapR_blks, 1 );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
      );

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, m_Work_mat, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, gecopy INFO = {}\n", info
      );

      info = gecopy( n, n, m_DE_blk, n, m_Work_mat+row00, m_Work_lda );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, gecopy INFO = {}\n", info
      );

      trsm(
        LEFT, LOWER, NO_TRANSPOSE, UNIT,
        n, nrhs, 1.0, m_Work_mat, m_Work_lda, io, ldRhs
      );
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        row00, nrhs, n,
        -1.0, m_Work_mat+n, m_Work_lda,
              io, ldRhs,
         1.0, io+n, ldRhs
      );
      io += n;

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
      integer k = 0;
      while ( ++k < nblock ) {

        valueType * B    = m_DE_blk + nxnx2 * k;
        valueType * A    = B - nxn;
        integer   * swps = m_swapR_blks + k * n;

        // apply permutation
        info = swaps( nrhs, io, ldRhs, 0, n-1, swps, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
        );

        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, m_Work_mat, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n, n, B, n, m_Work_mat+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        trsm(
          LEFT, LOWER, NO_TRANSPOSE, UNIT,
          n, nrhs, 1.0, m_Work_mat, m_Work_lda, io, ldRhs
        );
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          row00, nrhs, n,
          -1.0, m_Work_mat+n, m_Work_lda,
                io, ldRhs,
           1.0, io+n, ldRhs
        );
        io += n;
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
      this->la_factorization->solve( nrhs, io, ldRhs );
    
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

        valueType * B = m_DE_blk + nxnx2 * k;
        valueType * C = B + nxn;
        valueType * A = B - nxn;
        valueType * F = m_F_mat + k*m_F_size;
      
        // copy from
        info = gecopy( row00, n, A + n_m_row00, n, m_Work_mat, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
        );

        info = gecopy( row00, n, F, m_F_lda, m_Work_mat1, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info
        );

        info = gecopy( n_m_row00, n, B, n, m_Work_mat+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        info = gecopy( n_m_row00, n, C, n, m_Work_mat1+row00, m_Work_lda );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, gecopy INFO = {}\n", info
        );

        io -= n;
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          n, nrhs, n,
          -1.0, m_Work_mat1, m_Work_lda,
                io+n,        ldRhs,
           1.0, io,          ldRhs
        );

        trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          n, nrhs, 1.0, m_Work_mat, m_Work_lda, io, ldRhs
        );

      }

      valueType * B = m_DE_blk;
      valueType * C = B + nxn;

      // copy from
      info = gecopy( row00, n, block00 + col00, row0, m_Work_mat, m_Work_lda );
      UTILS_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info );

      info = gecopy( row00, n, m_F_mat, m_F_lda, m_Work_mat1, m_Work_lda );
      UTILS_ASSERT( info == 0, "BlockLU::solve, U block, gecopy INFO = {}\n", info );

      info = gecopy( n_m_row00, n, B, n, m_Work_mat+row00, m_Work_lda );
      UTILS_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = {}\n", info );

      info = gecopy( n_m_row00, n, C, n, m_Work_mat1+row00, m_Work_lda );
      UTILS_ASSERT( info == 0, "BlockLU::factorize, U block, gecopy INFO = {}\n", info );

      io -= n;
      gemm(
        NO_TRANSPOSE, NO_TRANSPOSE,
        n, nrhs, n,
        -1.0, m_Work_mat1, m_Work_lda,
              io+n,        ldRhs,
         1.0, io,          ldRhs
      );

      trsm(
        LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
        n, nrhs, 1.0, m_Work_mat, m_Work_lda, io, ldRhs
      );

      /*
      // fattorizzo primo blocco
      //
      //  +-----+---------+
      //  | L\U |         |
      //  |  M  |         |
      //  +-----+---------+
      */
      if ( col00 > 0 ) {
        io -= col00;
        gemm(
          NO_TRANSPOSE, NO_TRANSPOSE,
          col00, nrhs, n,
          -1.0, block00,  row0,
                io+col00, ldRhs,
           1.0, io,       ldRhs
        );
        trsm(
          LEFT, UPPER, NO_TRANSPOSE, NON_UNIT,
          col00, nrhs, 1.0, m_block0, row0, io, ldRhs
        );
      }
    } else {  // case nblock = 0
      this->la_factorization->solve( nrhs, in_out, ldRhs );
    }

    // permuto le x
    if ( do_permute ) {
      io = in_out;
      for ( integer k = 0; k < nrhs; ++k ) {
        std::rotate( io, io + col00, io + neq );
        io += ldRhs;
      }
    }
  }

  template class BlockLU<double>;
  template class BlockLU<float>;

}

///
/// eof: ABD_Block.cc
///

