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

///
/// file: ABD_Block.cc
///

#include "Alglin.hh"
#include "Alglin_Eigen.hh"

namespace alglin {

  using std::rotate;

  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/
  template <typename t_Value>
  void
  BlockLU<t_Value>::allocate_top_bottom(
    integer nblock,
    integer n,
    integer row0,
    integer col0,
    integer rowN,
    integer colN,
    integer nb
  ) {
    m_F_ldim    = n + row0 - col0;
    m_F_size    = n*m_F_ldim;
    m_Work_ldim = m_F_ldim + n;
    integer Fnnz = nblock*m_F_size;
    integer innz = nblock*n;
    BlockBidiagonal<t_Value>::allocate_top_bottom(
      nblock , n,
      row0, col0,
      rowN, colN,
      nb,
      Fnnz+2*m_Work_ldim*n,
      innz+row0
    );
    m_swap0      = m_mem_int( row0 );
    m_swapR_blks = m_mem_int( innz );
    m_F_mat      = m_mem( Fnnz );
    m_Work_mat   = m_mem( m_Work_ldim*n );
    m_Work_mat1  = m_mem( m_Work_ldim*n );
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
      m_num_cyclic_OMEGA == 0 && m_num_cyclic_BC == 0,
      "BlockLU cannot manage cyclic BC\n"
    );

    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;

    integer const & col00  = m_num_initial_OMEGA;
    integer const & colNN  = m_num_final_OMEGA;
    integer const & row0   = m_num_initial_BC;
    integer const & rowN   = m_num_final_BC;

    integer const row00     = row0 - col00;
    integer const n_m_row00 = n - row00;

    if ( nblock > 0 ) {

      /*
      // fattorizzo primo blocco
      //
      //  ┌─────┬─────────┐
      //  │ L\U │ M       │
      //  │     │         │
      //  └─────┴─────────┘
      */
      t_Value * block00 = m_block0 + col00*row0;
      if ( col00 > 0 ) {
        integer info = getrf( row0, col00, m_block0, row0, m_swap0 );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, getrf INFO = {}\n", info
        );
        // applico permutazioni
        info = alglin::swaps( n, block00, row0, 0, col00-1, m_swap0, 1 );
        UTILS_ASSERT(
          info == 0, "BlockLU::factorize, first block, swaps INFO = {}\n", info
        );
        trsm(
          SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          col00, n, 1.0, m_block0, row0, block00, row0
        );
        gemm(
          Transposition::NO,
          Transposition::NO,
          row00, n, col00,
          -1.0, m_block0+col00, row0,
                block00,        row0,
           1.0, block00+col00,  row0
        );
      }

      /*
      // ┌────────────┐
      // │            │ * * * * * * <-- fill in
      // ├────────────┼────────────┐
      // │            │            │
      // │      A     │     B      │
      // │            │            │
      // │            │            │
      // └────────────┼────────────┼────────────┐
      //              │            │            │
│
      //
      */

      // copy from
      GEcopy( row00, n, block00 + col00,  row0, m_Work_mat,        m_Work_ldim );
      GEzero( row00, n,                         m_Work_mat1,       m_Work_ldim );
      GEcopy( n,     n, m_DE_blk,         n,    m_Work_mat+row00,  m_Work_ldim );
      GEcopy( n,     n, m_DE_blk + n_x_n, n,    m_Work_mat1+row00, m_Work_ldim );

      // factorize
      integer info = getrf( n+row00, n, m_Work_mat, m_Work_ldim, m_swapR_blks );
      UTILS_ASSERT(
        info == 0,
        "BlockLU::factorize, first block, getrf INFO = {}\n", info
      );

      // apply permutation
      info = alglin::swaps( n, m_Work_mat1, m_Work_ldim, 0, n-1, m_swapR_blks, 1 );
      UTILS_ASSERT(
        info == 0,
        "BlockLU::factorize, first block, swaps INFO = {}\n", info
      );
      trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, n, 1.0, m_Work_mat,  m_Work_ldim,
                   m_Work_mat1, m_Work_ldim
      );
      gemm(
        Transposition::NO,
        Transposition::NO,
        row00, n, n,
        -1.0, m_Work_mat+n,  m_Work_ldim,
              m_Work_mat1,   m_Work_ldim,
         1.0, m_Work_mat1+n, m_Work_ldim
      );

      // copy to
      GEcopy( row00, n, m_Work_mat,        m_Work_ldim, block00 + col00,  row0 );
      GEcopy( row00, n, m_Work_mat1,       m_Work_ldim, m_F_mat,          m_F_ldim );
      GEcopy( n,     n, m_Work_mat+row00,  m_Work_ldim, m_DE_blk,         n );
      GEcopy( n,     n, m_Work_mat1+row00, m_Work_ldim, m_DE_blk + n_x_n, n );

      /*
      // │    A       │
      // │            │ * * * * * * <-- fill in
      // ├────────────┼────────────┬────────────┤
      // │            │            │
      // │            │            │
      // │     B      │    C       │
      // │            │            │
      // ├────────────┼────────────┼────────────┤
      //              │            │            │
      //              │            │            │
      //
      */

      for ( integer k{1}; k < nblock; ++k ) {

        real_type * B    { m_DE_blk + (2*k)*n_x_n };
        real_type * A    { B - n_x_n              };
        real_type * C    { B + n_x_n              };
        real_type * F    { m_F_mat + k*m_F_size   };
        integer   * swps { m_swapR_blks + k * n   };

        // copy from
        GEcopy( row00, n, A + n_m_row00, n, m_Work_mat,        m_Work_ldim );
        GEzero( row00, n,                   m_Work_mat1,       m_Work_ldim );
        GEcopy( n,     n, B,             n, m_Work_mat+row00,  m_Work_ldim );
        GEcopy( n,     n, C,             n, m_Work_mat1+row00, m_Work_ldim );

        // factorize
        info = getrf( n+row00, n, m_Work_mat, m_Work_ldim, swps );
        UTILS_ASSERT( info == 0, "BlockLU::factorize, first block, getrf INFO = {}\n", info );

        // apply permutation
        info = alglin::swaps( n, m_Work_mat1, m_Work_ldim, 0, n-1, swps, 1 );
        UTILS_ASSERT( info == 0, "BlockLU::factorize, first block, swaps INFO = {}\n", info );
        trsm(
          SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          n, n, 1.0, m_Work_mat,  m_Work_ldim,
                     m_Work_mat1, m_Work_ldim
        );
        gemm(
          Transposition::NO,
          Transposition::NO,
          row00, n, n,
          -1.0, m_Work_mat+n,  m_Work_ldim,
                m_Work_mat1,   m_Work_ldim,
           1.0, m_Work_mat1+n, m_Work_ldim
        );

        // copy to
        GEcopy( row00, n, m_Work_mat,        m_Work_ldim, A + n_m_row00, n );
        GEcopy( row00, n, m_Work_mat1,       m_Work_ldim, F,             m_F_ldim );
        GEcopy( n,     n, m_Work_mat+row00,  m_Work_ldim, B,             n );
        GEcopy( n,     n, m_Work_mat1+row00, m_Work_ldim, C,             n );
      }
      /*
      // │    A       │      B     │
      // │            │ ********** │
      // ├────────────┼────────────┴─┐
      //              │              │
      //              │   BlockN     │
      //              │              │
      //              └──────────────┘
      //
      */

      // fattorizzazione ultimo blocco
      real_type * B{ m_DE_blk + (2*nblock-1)*n_x_n };
      integer     N{ row00 + rowN };
      m_la_matrix.setup( N, N );
      m_la_matrix.zero_block(row00,N-n,0,n);
      m_la_matrix.load_block(row00,n,B+n_m_row00,n,0,0);
      m_la_matrix.load_block(rowN,N,m_blockN,rowN,row00,0);
      m_la_factorization->factorize( "BlockLU::factorize", m_la_matrix );
    } else { // case nblock == 0
      integer N = row0+rowN;
      m_la_matrix.setup( N, N );
      m_la_matrix.zero_block(N,N,0,0);
      m_la_matrix.load_block(row0,n+col00,m_block0,row0,0,0);
      m_la_matrix.load_block(rowN,n+colNN,m_blockN,rowN,row0,col00);
      m_la_factorization->factorize( "BlockLU::factorize[nblock=0]", m_la_matrix );
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
    bool const do_permute,
    real_type  in_out[]
  ) const {

    integer const & nblock  { m_number_of_blocks };
    integer const & n       { m_block_size };

    integer const & col00   { m_num_initial_OMEGA };
    integer const & row0    { m_num_initial_BC };
    integer const & rowN    { m_num_final_BC };

    integer const row00     { row0 - col00 };
    integer const n_m_row00 { n - row00 };

    // permuto le x
    integer neq{ nblock * n + row0 + rowN };
    if ( do_permute ) rotate( in_out, in_out + neq - row0, in_out + neq );

    if ( nblock > 0 ) {

      real_type * io = in_out;
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
      //  ┌─────┬─────────┐
      //  │ L\U │         │
      //  │  M  │         │
      //  └─────┴─────────┘
      */
      if ( col00 > 0 ) {
        // apply permutation
        integer info{ alglin::swaps( 1, io, row0, 0, col00-1, m_swap0, 1 ) };
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
        );
        trsv(
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          col00, m_block0, row0, in_out, 1
        );
        io += col00;
        gemv(
          Transposition::NO,
          row00, col00,
          -1.0, m_block0+col00, row0,
                in_out, 1,
           1.0, io, 1
        );
      }

      /*
      // ┌─┬────────────┐
      // │ │            │ * * * * * * <-- fill in
      // └─┼────────────┼────────────┐
      //   │            │            │
      //   │            │            │
      //   │     A      │     B      │
      //   │            │            │
      //   └────────────┼────────────┼────────────┐
      //                │            │            │
      //                │            │            │
      //
      */
      real_type const * block00 { m_block0 + col00 * row0 };

      // apply permutation
      integer info = alglin::swaps( 1, io, n, 0, n-1, m_swapR_blks, 1 );
      UTILS_ASSERT(
        info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
      );

      // copy from
      GEcopy( row00, n, block00 + col00, row0, m_Work_mat,       m_Work_ldim );
      GEcopy( n,     n, m_DE_blk,        n,    m_Work_mat+row00, m_Work_ldim );

      trsv(
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, m_Work_mat, m_Work_ldim, io, 1
      );
      gemv(
        Transposition::NO,
        row00, n,
        -1.0, m_Work_mat+n, m_Work_ldim,
              io, 1,
         1.0, io+n, 1
      );
      io += n;

      /*
      // │            │
      // │            │ * * * * * * <-- fill in
      // ├────────────┼────────────┐
      // │            │            │
      // │            │            │
      // │            │    A       │
      // │            │            │
      // └────────────┼────────────┼────────────┐
      //              │            │            │
      //              │    B       │            │
      //
      */
      integer k{0};
      while ( ++k < nblock ) {

        real_type * B    { m_DE_blk + (2*k)*n_x_n };
        real_type * A    { B - n_x_n };
        integer   * swps { m_swapR_blks + k * n };

        // apply permutation
        info = alglin::swaps( 1, io, n, 0, n-1, swps, 1 );
        UTILS_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info );

        // copy from
        GEcopy( row00, n, A + n_m_row00, n, m_Work_mat,       m_Work_ldim );
        GEcopy( n,     n, B,             n, m_Work_mat+row00, m_Work_ldim );

        trsv(
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          n, m_Work_mat, m_Work_ldim, io, 1
        );
        gemv(
          Transposition::NO,
          row00, n,
          -1.0, m_Work_mat+n, m_Work_ldim,
                io, 1,
           1.0, io+n, 1
        );
        io += n;
      }

      /*
      // │    A       │      B     │
      // │            │            │
      // └────────────┼────────────┴─┐
      //              │              │
      //              │   BlockN     │
      //              │              │
      //              └──────────────┘
      //
      */

      // risolvo ultimo blocco
      m_la_factorization->solve( io );

      /*
      //             _             _   _
      //   ___  ___ | |_   _____  | | | |
      //  / __|/ _ \| \ \ / / _ \ | | | |
      //  \__ \ (_) | |\ V /  __/ | |_| |
      //  |___/\___/|_| \_/ \___|  \___/
      */
      /*
      // │     A      │
      // │            │ * * * * * * <-- fill in
      // ├────────────┼────────────┐
      // │            │            │
      // │            │            │
      // │     B      │     C      │
      // │            │            │ * * * * * * <- fill in
      // └────────────┼────────────┼────────────┐
      //              │            │            │
      //              │            │            │
      //
      */
      while ( --k > 0 ) {

        real_type * B { m_DE_blk + (2*k)*n_x_n };
        real_type * C { B + n_x_n };
        real_type * A { B - n_x_n };
        real_type * F { m_F_mat + k*m_F_size };

        // copy from
        GEcopy( row00,     n, A + n_m_row00, n,        m_Work_mat,        m_Work_ldim );
        GEcopy( row00,     n, F,             m_F_ldim, m_Work_mat1,       m_Work_ldim );
        GEcopy( n_m_row00, n, B,             n,        m_Work_mat+row00,  m_Work_ldim );
        GEcopy( n_m_row00, n, C,             n,        m_Work_mat1+row00, m_Work_ldim );

        io -= n;
        gemv(
          Transposition::NO,
          n, n,
          -1.0, m_Work_mat1, m_Work_ldim,
                io+n,        1,
           1.0, io,          1
        );

        trsv(
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, m_Work_mat, m_Work_ldim, io, 1
        );

      }

      real_type * B { m_DE_blk };
      real_type * C { B + n_x_n };

      // copy from
      GEcopy( row00,     n, block00 + col00, row0,     m_Work_mat,        m_Work_ldim );
      GEcopy( row00,     n, m_F_mat,         m_F_ldim, m_Work_mat1,       m_Work_ldim );
      GEcopy( n_m_row00, n, B,               n,        m_Work_mat+row00,  m_Work_ldim );
      GEcopy( n_m_row00, n, C,               n,        m_Work_mat1+row00, m_Work_ldim );

      io -= n;
      gemv(
        Transposition::NO,
        n, n,
        -1.0, m_Work_mat1, m_Work_ldim,
              io+n,        1,
         1.0, io,          1
      );

      trsv(
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        n, m_Work_mat, m_Work_ldim, io, 1
      );

      /*
      // fattorizzo primo blocco
      //
      //  ┌─────┬─────────┐
      //  │ L\U │         │
      //  │  M  │         │
      //  └─────┴─────────┘
      */
      if ( col00 > 0 ) {
        io -= col00;
        gemv(
          Transposition::NO,
          col00, n,
          -1.0, block00,  row0,
                io+col00, 1,
           1.0, io,       1
        );
        trsv(
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          col00, m_block0, row0, io, 1
        );
      }
    } else {  // case nblock = 0
      m_la_factorization->solve( in_out );
    }

    // permuto le x
    if ( do_permute ) rotate( in_out, in_out + col00, in_out + neq );
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
    bool const do_permute,
    integer    nrhs,
    real_type  in_out[],
    integer    ldRhs
  ) const {

    integer const & nblock  { m_number_of_blocks };
    integer const & n       { m_block_size };

    integer const & col00   { m_num_initial_OMEGA };
    integer const & row0    { m_num_initial_BC };
    integer const & rowN    { m_num_final_BC };

    integer const row00     { row0 - col00 };
    integer const n_m_row00 { n - row00 };

    integer neq{ nblock * n + row0 + rowN };

    // permuto le x
    real_type * io{ in_out };
    if ( do_permute ) {
      for ( integer k{0}; k < nrhs; ++k ) {
        rotate( io, io + neq - row0, io + neq );
        io += ldRhs;
      }
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
      //  ┌─────┬─────────┐
      //  │ L\U │         │
      //  │  M  │         │
      //  └─────┴─────────┘
      */
      if ( col00 > 0 ) {
        // apply permutation
        integer info{ alglin::swaps( 1, io, row0, 0, col00-1, m_swap0, 1 ) };
        UTILS_ASSERT(
          info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info
        );
        trsm(
          SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          col00, nrhs, 1.0, m_block0, row0, in_out, ldRhs
        );
        io += col00;
        gemm(
          Transposition::NO,
          Transposition::NO,
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
      real_type const * block00 { m_block0 + col00*row0 };

      // apply permutation
      integer info { alglin::swaps( nrhs, io, ldRhs, 0, n-1, m_swapR_blks, 1 ) };
      UTILS_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info );

      // copy from
      GEcopy( row00, n, block00 + col00, row0, m_Work_mat,       m_Work_ldim );
      GEcopy( n,     n, m_DE_blk,        n,    m_Work_mat+row00, m_Work_ldim );

      trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, nrhs, 1.0, m_Work_mat, m_Work_ldim,
                      io,         ldRhs
      );
      gemm(
        Transposition::NO,
        Transposition::NO,
        row00, nrhs, n,
        -1.0, m_Work_mat+n, m_Work_ldim,
              io,           ldRhs,
         1.0, io+n,         ldRhs
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
      integer k{0};
      while ( ++k < nblock ) {

        real_type * B    { m_DE_blk + (2*k)*n_x_n };
        real_type * A    { B - n_x_n };
        integer   * swps { m_swapR_blks + k * n };

        // apply permutation
        info = alglin::swaps( nrhs, io, ldRhs, 0, n-1, swps, 1 );
        UTILS_ASSERT( info == 0, "BlockLU::solve, first block, swaps INFO = {}\n", info );

        // copy from
        GEcopy( row00, n, A + n_m_row00, n, m_Work_mat,       m_Work_ldim );
        GEcopy( n,     n, B,             n, m_Work_mat+row00, m_Work_ldim );

        trsm(
          SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO,
          DiagonalType::UNIT,
          n, nrhs, 1.0, m_Work_mat, m_Work_ldim,
                        io,         ldRhs
        );
        gemm(
          Transposition::NO,
          Transposition::NO,
          row00, nrhs, n,
          -1.0, m_Work_mat+n, m_Work_ldim,
                io,           ldRhs,
           1.0, io+n,         ldRhs
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
      m_la_factorization->solve( nrhs, io, ldRhs );

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

        real_type * B { m_DE_blk + (2*k)*n_x_n };
        real_type * C { B + n_x_n };
        real_type * A { B - n_x_n };
        real_type * F { m_F_mat + k*m_F_size };

        // copy from
        GEcopy( row00,     n, A + n_m_row00, n,        m_Work_mat,        m_Work_ldim );
        GEcopy( row00,     n, F,             m_F_ldim, m_Work_mat1,       m_Work_ldim );
        GEcopy( n_m_row00, n, B,             n,        m_Work_mat+row00,  m_Work_ldim );
        GEcopy( n_m_row00, n, C,             n,        m_Work_mat1+row00, m_Work_ldim );

        io -= n;
        gemm(
          Transposition::NO,
          Transposition::NO,
          n, nrhs, n,
          -1.0, m_Work_mat1, m_Work_ldim,
                io+n,        ldRhs,
           1.0, io,          ldRhs
        );

        trsm(
          SideMultiply::LEFT,
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          n, nrhs, 1.0, m_Work_mat, m_Work_ldim,
                        io,         ldRhs
        );

      }

      real_type * B { m_DE_blk  };
      real_type * C { B + n_x_n };

      // copy from
      GEcopy( row00,     n, block00 + col00, row0,     m_Work_mat,        m_Work_ldim );
      GEcopy( row00,     n, m_F_mat,         m_F_ldim, m_Work_mat1,       m_Work_ldim );
      GEcopy( n_m_row00, n, B,               n,        m_Work_mat+row00,  m_Work_ldim );
      GEcopy( n_m_row00, n, C,               n,        m_Work_mat1+row00, m_Work_ldim );

      io -= n;
      gemm(
        Transposition::NO,
        Transposition::NO,
        n, nrhs, n,
        -1.0, m_Work_mat1, m_Work_ldim,
              io+n,        ldRhs,
         1.0, io,          ldRhs
      );

      trsm(
        SideMultiply::LEFT,
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        n, nrhs, 1.0, m_Work_mat, m_Work_ldim,
                      io,         ldRhs
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
          Transposition::NO,
          Transposition::NO,
          col00, nrhs, n,
          -1.0, block00,  row0,
                io+col00, ldRhs,
           1.0, io,       ldRhs
        );
        trsm(
          SideMultiply::LEFT,
          ULselect::UPPER,
          Transposition::NO,
          DiagonalType::NON_UNIT,
          col00, nrhs, 1.0, m_block0, row0, io, ldRhs
        );
      }
    } else {  // case nblock = 0
      m_la_factorization->solve( nrhs, in_out, ldRhs );
    }

    // permuto le x
    if ( do_permute ) {
      io = in_out;
      for ( integer k{0}; k < nrhs; ++k ) {
        rotate( io, io + col00, io + neq );
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

