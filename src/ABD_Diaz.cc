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
/// file: ABD_Diaz.cc
///

#include "Alglin.hh"

/*
//  Blocco algoritmo di Diaz
//
//   ┌─────┬─────┐
//   │  A  │  B  │
//   ├─────┼─────┼───────────┐
//   │     │     │           │
//   │  C  │  D  │     E     │
//   │     │     │           │
//   │     │     │           │
//   └─────┴─────┼───────────┴───────────┐
//               │                       │
//               │                       │
//               │                       │
//               │                       │
//               └───────────────────────┘
//
*/

namespace alglin {

  using std::rotate;
  using std::swap;

  /*
  //
  // ┌───┬─────┬────────┐
  // │   │     │        │
  // │ L │  A  │    R   │
  // │   │     │        │
  // └───┴─────┴────────┘
  //
  //
  // ┌───┬─────┬────────┐
  // │   │ \ U │        │
  // │ L │   \ │    R   │
  // │   │  L  │        │
  // └───┴─────┴────────┘
  //
  */

  template <typename t_Value>
  void
  DiazLU<t_Value>::LU_left_right(
    integer   nrA,
    integer   ncA,
    integer   ncL,
    integer   ncR,
    t_Value * A, integer ldA,
    integer   swapR[]
  ) {

    integer ierr;
    if ( 2*m_NB < nrA ) ierr = getry<t_Value>( nrA, ncA, A, ldA, swapR, m_NB );
    else                ierr = gty( nrA, ncA, A, ldA, swapR );

    UTILS_ASSERT( ierr == 0, "DiazLU::LU_left_right, found ierr: {}\n", ierr );
    // applico permutazione al blocco
    t_Value * R { A + ncA * ldA };
    t_Value * L { A - ncL * ldA };
    for ( integer i{0}; i < ncA; ++i ) {
      integer ip = swapR[i];
      if ( ip > i ) {
        alglin::swap( ncR, R + i, ldA, R + ip, ldA );
        alglin::swap( ncL, L + i, ldA, L + ip, ldA );
      }
    }
    alglin::trsm(
      SideMultiply::LEFT,
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      ncA, ncR, 1.0,
      A, ldA,
      R, ldA
    );
    alglin::gemm(
      Transposition::NO,
      Transposition::NO,
      nrA-ncA, ncR, ncA,
      -1.0,  A+ncA, ldA,
      R, ldA,
      1.0,  R+ncA, ldA
    );
  }

  /*
  //                     columns permuted
  //   +-----------+      +-----------+
  //   |     T     |      |     T     |
  //   +-----------+      +-----------+
  //   |     A     |      | L \   U   |
  //   +===========+  ==> +====+======+
  //   |           |      |    |      |
  //   |     B     |      | L  |  U   |
  //   |           |      |    |      |
  //   +-----------+      +----+------+
  */

  template <typename t_Value>
  void
  DiazLU<t_Value>::LU_top_bottom(
    integer   nrT,
    integer   nrA,
    integer   ncA,
    t_Value * A, integer ldA,
    integer   nrB,
    t_Value * B, integer ldB,
    integer   swapC[]
  ) {

    integer ierr;
    if ( 2*m_NB < ncA ) ierr = getrx( nrA, ncA, A, ldA, swapC, m_NB );
    else                ierr = gtx( nrA, ncA, A, ldA, swapC );

    UTILS_ASSERT( ierr == 0, "DiazLU::LU_top_bottom, found ierr: {}\n", ierr );
    // applico permutazione al blocco
    t_Value * T{A - nrT};
    for ( integer i{0}; i < nrA; ++i ) {
      if ( integer const ip{swapC[i]}; ip > i ) {
        alglin::swap( nrT, T + i*ldA, 1, T + ip*ldA, 1 );
        alglin::swap( nrB, B + i*ldB, 1, B + ip*ldB, 1 );
      }
    }
    alglin::trsm(
      SideMultiply::RIGHT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      nrB, nrA, 1.0,
      A, ldA,
      B, ldB
    );
    alglin::gemm(
      Transposition::NO,
      Transposition::NO,
      nrB, ncA-nrA, nrA,
      -1.0, B, ldB,
      A+nrA*ldA, ldA,
      1.0,  B+nrA*ldB, ldB
    );
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
  DiazLU<t_Value>::factorize() {

    UTILS_ASSERT0(
      m_num_cyclic_OMEGA == 0 && m_num_cyclic_BC == 0,
      "DiazLU cannot manage cyclic BC\n"
    );

    integer const & nblock  { m_number_of_blocks };
    integer const & n       { m_block_size };

    integer const & col00   { m_num_initial_OMEGA };
    integer const & row0    { m_num_initial_BC };
    integer const & rowN    { m_num_final_BC };

    integer const col0      { n + col00 };
    integer const row00     { row0 - col00 };
    integer const n_m_row00 { n - row00 };

    // primo blocco
    integer * swapRC{ m_swapRC_blks };
    if ( col00 > 0 ) {
      LU_left_right( row0, col00, 0, col0-col00, m_block0, row0, swapRC );
      swapRC += col00;
    }

    /*
    //    dimA   n
    //  +--///--+----------+
    //  |  ---  | B0       | ldA-row00
    //  +       +----------+
    //  |  A1   | B1       | row00
    //  +--///--+----+-----+----------+
    //          |    |     |          |
    //          | C  | D   |    E     |
    //          |    |     |          |
    //          +----+-----+----------+
    //          row00
    */

    // blocchi intermedi (n-1)
    if ( nblock > 0 ) {
      real_type * B1 = m_block0 + (row0+1)*col00;
      LU_top_bottom(
        col00, row00,
        n, B1, row0,
        n, m_DE_blk,
        n, swapRC
      );
      swapRC += row00;
      real_type * D{ m_DE_blk + row00 * n };
      //                 NR          NC            L       R
      LU_left_right( n, n_m_row00, row00, n, D, n, swapRC );
      swapRC += n_m_row00;
    }

    real_type * C{ m_DE_blk };
    for ( m_nblk = 1; m_nblk < nblock; ++m_nblk ) {
      C += 2*n_x_n;
      real_type * B1 { C - n_x_n + n_m_row00 };
      LU_top_bottom(
        n_m_row00, row00,
        n, B1, n,
        n,  C, n, swapRC
      );
      swapRC += row00;
      real_type * D { C + row00 * n };
      //                 NR          NC            L       R
      LU_left_right( n, n_m_row00, row00, n, D, n, swapRC );
      swapRC += n_m_row00;
    }

    // ultimo blocco
    /*
    //    dimA   n
    //  +--///--+----------+
    //  |  ---  | B0       | row11 = ldA-row00 o col00
    //  +       +----------+
    //  |  A1   | B1       |
    //  +--///--+----+-----+-----+
    //          | C0 | D0  |  E0 |
    //          |    |     |     | rowN
    //          | C1 | D1  |  E1 |
    //          +----+-----+-----+
    //          row00       colNN
    */

    swapRC = m_swapRC_blks + ( col00 + nblock*n );

    if ( nblock == 0 ) {
      real_type * B1 { m_block0 + (row0+1) * col00 };
      LU_top_bottom(
        col00, row00, n,
        B1, row0,
        rowN, m_blockN, rowN,
        swapRC
      );
    } else {
      real_type * B1 { m_DE_blk + (2*nblock-1) * n_x_n + n_m_row00 };
      LU_top_bottom(
        n_m_row00, row00, n,
        B1, n,
        rowN, m_blockN, rowN,
        swapRC
      );
    }

    // fattorizzazione ultimo blocco
    real_type * D0 { m_blockN + row00 * rowN };
    m_la_factorization->factorize(
      "DiazLU::factorize", rowN, rowN, D0, rowN
    );
  }

  /*\
   |             _               _       _                        _
   |   ___  ___ | |_   _____    (_)_ __ | |_ ___ _ __ _ __   __ _| |
   |  / __|/ _ \| \ \ / / _ \   | | '_ \| __/ _ \ '__| '_ \ / _` | |
   |  \__ \ (_) | |\ V /  __/   | | | | | ||  __/ |  | | | | (_| | |
   |  |___/\___/|_| \_/ \___|___|_|_| |_|\__\___|_|  |_| |_|\__,_|_|
   |                       |_____|
  \*/
  // ---------------------------------------------------------------------------
  template <typename t_Value>
  void
  DiazLU<t_Value>::solve_internal(
    bool const do_permute,
    real_type  in_out[]
  ) const {

    integer const & nblock  { m_number_of_blocks };
    integer const & n       { m_block_size };

    integer const & col00   { m_num_initial_OMEGA };
    integer const & colNN   { m_num_final_OMEGA };
    integer const & row0    { m_num_initial_BC };
    integer const & rowN    { m_num_final_BC };

    integer const col0      { n + col00 };
    integer const colN      { n + colNN };
    integer const row00     { row0 - col00 };
    integer const n_m_row00 { n - row00 };

    integer neq { nblock * n + row0 + rowN };
    if ( do_permute ) rotate( in_out, in_out + neq - row0, in_out + neq );

    // applico permutazione alla RHS
    integer const * swapR { m_swapRC_blks };
    real_type * io { in_out };
    for ( integer k{0}; k < col00; ++k ) {
      integer k1{ swapR[k] }; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] );
    }
    io    += row0;
    swapR += row0;
    for ( m_nblk = 0; m_nblk < nblock; ++m_nblk )  {
      for ( integer k{0}; k < n_m_row00; ++k ) {
        integer k1 = swapR[k]; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] );
      }
      io    += n;
      swapR += n;
    }

    // primo blocco
    io = in_out;
    trsv(
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      row0,
      m_block0, row0, io, 1
    );

    io += row0;
    // blocchi intermedi
    m_nblk = 0;
    while ( m_nblk < nblock ) {

      /*
      //
      //  +----+---+---+---+
      //  |  M | L |   :   |
      //  |  M | L | L :   |
      //  +----+---+---+---+
      //  row00
      */

      real_type * io1 { io - row00 };
      real_type * M   { m_DE_blk + (2*m_nblk) * n_x_n };
      real_type * L   { M + row00 * n };

      // io -= M*io1
      gemv(
        Transposition::NO,
        n, row00,
        -1, M, n,
        io1, 1,
        1, io, 1
      );

      trsv(
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, L, n, io, 1
      );

      io += n;
      ++m_nblk;
    }

    // soluzione ultimo blocco
    integer ncol { colN - rowN };
    gemv(
      Transposition::NO,
      rowN, ncol,
      -1, m_blockN, rowN,
      io-ncol, 1,
      1, io, 1
    );

    m_la_factorization->solve(io);

    while ( m_nblk > 0 ) {
      --m_nblk;
      io -= n;

      /*
      //  +---------+---+-------+
      //  |       U |   |       |
      //  |         | U |   M   |
      //  +---------+---+-------+
      //  row00
      */

      real_type * io1 { io + n };
      real_type * U   { m_DE_blk + (2*m_nblk)*n_x_n + row00 * n };
      real_type * M   { U + n_x_n };

      gemv(
        Transposition::NO,
        n, n_m_row00,
        -1, M, n,
        io1, 1,
        1, io, 1
      );

      trsv(
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        n, U, n, io, 1
      );
    }

    // primo blocco
    io -= row0;

    // soluzione primo blocco
    gemv(
      Transposition::NO,
      row0, col0-row0,
      -1, m_block0+row0*row0, row0,
      io+row0, 1,
      1, io, 1
    );

    trsv(
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      row0,
      m_block0, row0, io, 1
    );

    // applico permutazione alla Soluzione
    integer const * swapC { m_swapRC_blks + ( nblock * n + col00 ) };
    io = in_out + nblock*n + col00;
    integer k{ row00 };
    while ( k > 0 ) {
      integer k1 { swapC[--k] }; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] );
    }
    for ( m_nblk = 0; m_nblk < nblock; ++m_nblk )  {
      io    -= n;
      swapC -= n;
      k      = row00;
      while ( k > 0 ) {
        integer k1{ swapC[--k] }; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] );
      }
    }

    // permuto le x
    if ( do_permute ) rotate( in_out, in_out + col00, in_out + neq );
  }

  template <typename t_Value>
  void
  DiazLU<t_Value>::solve_internal(
    bool const do_permute,
    integer    nrhs,
    real_type  in_out[],
    integer    ldRhs
  ) const {

    integer const & nblock  { m_number_of_blocks };
    integer const & n       { m_block_size };

    integer const & col00   { m_num_initial_OMEGA };
    integer const & colNN   { m_num_final_OMEGA };
    integer const & row0    { m_num_initial_BC };
    integer const & rowN    { m_num_final_BC };

    integer const col0      { n + col00 };
    integer const colN      { n + colNN };
    integer const row00     { row0 - col00 };
    integer const n_m_row00 { n - row00 };

    integer neq { nblock * n + row0 + rowN };

    // permuto le x
    real_type * io { in_out };
    if ( do_permute ) {
      for ( integer k{0}; k < nrhs; ++k ) {
        rotate( io, io + neq - row0, io + neq );
        io += ldRhs;
      }
      io = in_out;
    }

    // applico permutazione alla RHS
    integer const * swapR { m_swapRC_blks };
    for ( integer k{0}; k < col00; ++k ) {
      integer k1 { swapR[k] }; // 0 based
      if ( k1 > k )
        alglin::swap( nrhs, io+k, ldRhs, io+k1, ldRhs );
    }
    io    += row0;
    swapR += row0;
    for ( m_nblk = 0; m_nblk < nblock; ++m_nblk )  {
      for ( integer k{0}; k < n_m_row00; ++k ) {
        integer k1 { swapR[k] }; // 0 based
        if ( k1 > k )
          alglin::swap( nrhs, io+k, ldRhs, io+k1, ldRhs );
      }
      io    += n;
      swapR += n;
    }

    // primo blocco
    io = in_out;
    alglin::trsm(
      SideMultiply::LEFT,
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      row0, nrhs,
      1.0, m_block0, row0,
      io, ldRhs
    );

    io += row0;
    // blocchi intermedi
    m_nblk = 0;
    while ( m_nblk < nblock ) {

      /*
      //
      //  +----+---+---+---+
      //  |  M | L |   :   |
      //  |  M | L | L :   |
      //  +----+---+---+---+
      //  row00
      */

      real_type * io1 { io - row00 };
      real_type * M   { m_DE_blk + (2*m_nblk) * n_x_n };
      real_type * L   { M + row00 * n };

      // io -= M*io1
      alglin::gemm(
        Transposition::NO,
        Transposition::NO,
        n, nrhs, row00,
        -1, M, n,
        io1, ldRhs,
        1, io, ldRhs
      );

      alglin::trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, nrhs,
        1.0, L, n,
        io, ldRhs
      );

      io += n;
      ++m_nblk;
    }

    // soluzione ultimo blocco
    integer ncol { colN - rowN };
    alglin::gemm(
      Transposition::NO,
      Transposition::NO,
      rowN, nrhs, ncol,
      -1, m_blockN, rowN,
      io-ncol, ldRhs,
      1, io, ldRhs
    );

    m_la_factorization->solve(nrhs,io,ldRhs);

    while ( m_nblk > 0 ) {
      --m_nblk;
      io -= n;

      /*
      //  +---------+---+-------+
      //  |       U |   |       |
      //  |         | U |   M   |
      //  +---------+---+-------+
      //  row00
      */

      real_type * io1 { io + n };
      real_type * U   { m_DE_blk + (2*m_nblk) * n_x_n + row00 * n };
      real_type * M   { U + n_x_n };

      alglin::gemm(
        Transposition::NO,
        Transposition::NO,
        n, nrhs, n_m_row00,
        -1, M, n,
        io1, ldRhs,
        1, io, ldRhs
      );

      alglin::trsm(
        SideMultiply::LEFT,
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        n, nrhs,
        1.0, U, n,
        io, ldRhs
      );
    }

    // primo blocco
    io -= row0;

    // soluzione primo blocco
    alglin::gemm(
      Transposition::NO,
      Transposition::NO,
      row0, nrhs, col0-row0,
      -1, m_block0+row0*row0, row0,
      io+row0, ldRhs,
      1, io, ldRhs
    );

    alglin::trsm(
      SideMultiply::LEFT,
      ULselect::UPPER,
      Transposition::NO,
      DiagonalType::NON_UNIT,
      row0, nrhs,
      1.0, m_block0, row0,
      io, ldRhs
    );

    // applico permutazione alla Soluzione
    integer const * swapC { m_swapRC_blks + ( nblock * n + col00 ) };
    io = in_out + nblock*n + col00;
    integer k { row00 };
    while ( k > 0 ) {
      integer k1 { swapC[--k] }; // 0 based
      if ( k1 > k )
        alglin::swap( nrhs, io+k, ldRhs, io+k1, ldRhs );
    }
    for ( m_nblk = 0; m_nblk < nblock; ++m_nblk )  {
      io    -= n;
      swapC -= n;
      k      = row00;
      while ( k > 0 ) {
        integer k1 { swapC[--k] }; // 0 based
        if ( k1 > k )
          alglin::swap( nrhs, io+k, ldRhs, io+k1, ldRhs );
      }
    }

    // permuto le x
    if ( do_permute ) {
      io = in_out;
      for ( k = 0; k < nrhs; ++k ) {
        rotate( io, io + col00, io + neq );
        io += ldRhs;
      }
    }
  }

  template class DiazLU<double>;
  template class DiazLU<float>;

}

///
/// eof: ABD_Diaz.cc
///

