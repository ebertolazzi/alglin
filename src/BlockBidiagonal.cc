/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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
#include "Alglin_Eigen.hh"

namespace alglin {

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_blocks( real_type const AdAu[], integer ldA ) {
    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    GEcopy( n, 2*n*nblock, AdAu, ldA, m_DE_blk, n );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_block( integer const nbl, real_type const AdAu[], integer const ldA ) {
    integer const & n{m_block_size};
    GEcopy( n, 2*n, AdAu, ldA, m_DE_blk + 2*nbl*n_x_n, n );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_block_left( integer const nbl, real_type const Ad[], integer const ldA ) {
    integer const & n{m_block_size};
    GEcopy( n, n, Ad, ldA, m_DE_blk + 2*nbl*n_x_n, n );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_block_right( integer const nbl, real_type const Au[], integer const ldA ) {
    integer const & n{m_block_size};
    GEcopy( n, n, Au, ldA, m_DE_blk + (2*nbl+1)*n_x_n, n );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::set_zero_bottom_blocks() {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    alglin::Zero_n( m_Cmat, nb*neq );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_bottom_blocks( real_type const C[], integer ldC ) {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    GEcopy( nb, neq, C, ldC, m_Cmat, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_bottom_block( integer nbl, real_type const C[], integer ldC ) {
    integer const & n{m_block_size};
    integer const & nb{m_border_size};
    UTILS_ASSERT( ldC >= nb, "load_bottom_block( {}, C, ldC = {} ) bad ldC\n", nbl, ldC );
    real_type * CC{m_Cmat + nbl*n_x_nb};
    GEcopy( nb, n, C, ldC, CC, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::add_to_bottom_block( integer nbl, real_type const C[], integer ldC ) {
    integer const & n{m_block_size};
    integer const & nb{m_border_size};
    UTILS_ASSERT( ldC >= nb, "add_to_bottom_block( {}, C, ldC = {} ) bad ldC\n", nbl, ldC );
    real_type * CC{m_Cmat + nbl*n_x_nb};
    GEadd( nb, n, C, ldC, CC, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::add_to_bottom_block2( integer nbl, real_type const C[], integer ldC ) {
    integer const & n{m_block_size};
    integer const & nb{m_border_size};
    UTILS_ASSERT( ldC >= nb, "add_to_bottom_block2( {}, C, ldC = {} ) bad ldC\n", nbl, ldC );
    real_type * CC{m_Cmat + nbl*n_x_nb};
    GEadd( nb, 2*n, C, ldC, CC, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_bottom_last_block( real_type const C[], integer ldC ) {
    integer const & nblock{m_number_of_blocks};
    integer const & nb{m_border_size};
    integer const & q{m_extra_bc};
    GEcopy( nb, q, C, ldC, m_Cmat + (nblock+1)*n_x_nb, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::set_zero_right_blocks() {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    alglin::Zero_n( m_Bmat, neq*nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_right_blocks( real_type const B[], integer ldB ) {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    GEcopy( neq, nb, B, ldB, m_Bmat, neq );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::loadRightBlock( integer const nbl, real_type const B[], integer const ldB ) {
    integer const & n{m_block_size};
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    GEcopy( n, nb, B, ldB, m_Bmat + nbl*n, neq );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::loadRightLastBlock( real_type const B[], integer ldB ) {
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    GEcopy( n+q, nb, B, ldB, m_Bmat + neq-n-q, neq );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::setZeroRBblock() {
    integer const & nb{m_border_size};
    alglin::Zero_n( m_Dmat, nb*nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_RB_block( real_type const D[], integer ldD ) {
    integer const & nb{m_border_size};
    GEcopy( nb, nb, D, ldD, m_Dmat, nb );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_LR( real_type LR[], integer ldA ) const {
    integer const & n{m_block_size};
    GEcopy( n, 2*n, m_DE_blk, n, LR, ldA );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_L( real_type L[], integer ldA ) const {
    integer const & n{m_block_size};
    GEcopy( n, n, m_DE_blk, n, L, ldA );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_R( real_type R[], integer ldA ) const {
    integer const & n{m_block_size};
    GEcopy( n, n, m_DE_blk+n_x_n, n, R, ldA );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_H0( real_type H0[], integer ld0 ) const {
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};
    GEcopy( n+q, n, m_H0Nq, n+q, H0, ld0 );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_HN( real_type HN[], integer ldN ) const {
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};
    integer nq{n + q};
    GEcopy( nq, n, m_H0Nq+n*nq, nq, HN, ldN );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::getBlock_Hq( real_type Hq[], integer ldQ ) const {
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};
    integer nq{n + q};
    GEcopy( nq, q, m_H0Nq+2*n*nq, nq, Hq, ldQ );
  }

  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::allocate(
    integer const nblock,
    integer const n,
    integer const nb,
    // ----------------------
    integer const num_initial_BC,
    integer const num_final_BC,
    integer const num_cyclic_BC,
    // ----------------------
    integer const num_initial_OMEGA,
    integer const num_final_OMEGA,
    integer const num_cyclic_OMEGA,
    // ----------------------
    integer const num_extra_r,
    integer const num_extra_i
  ) {

    UTILS_ASSERT(
      num_initial_BC  >= 0 && num_final_BC      >= 0 &&
      num_cyclic_BC   >= 0 && num_initial_OMEGA >= 0 &&
      num_final_OMEGA >= 0 && num_cyclic_OMEGA  >= 0,
      "Bad BC specification:"
      "\nnum_initial_BC    = {}"
      "\nnum_final_BC      = {}"
      "\nnum_cyclic_BC     = {}"
      "\nnum_initial_OMEGA = {}"
      "\nnum_final_OMEGA   = {}"
      "\nnum_cyclic_OMEGA  = {}\n",
      num_initial_BC,    num_final_BC,    num_cyclic_BC,
      num_initial_OMEGA, num_final_OMEGA, num_cyclic_OMEGA
    );

    m_extra_bc = num_initial_OMEGA + num_final_OMEGA + num_cyclic_OMEGA;
    integer const & q = m_extra_bc;

    UTILS_ASSERT(
      num_initial_BC + num_final_BC + num_cyclic_BC == n+q,
      "Bad BC specification:"
      "\nnum_initial_BC    = {}"
      "\nnum_final_BC      = {}"
      "\nnum_cyclic_BC     = {}"
      "\nnum_initial_OMEGA + num_final_OMEGA + num_cyclic_OMEGA must be = {}\n",
      num_initial_BC, num_final_BC, num_cyclic_BC, n+q
    );

    m_num_initial_BC    = num_initial_BC;
    m_num_final_BC      = num_final_BC;
    m_num_cyclic_BC     = num_cyclic_BC;
    m_num_initial_OMEGA = num_initial_OMEGA;
    m_num_final_OMEGA   = num_final_OMEGA;
    m_num_cyclic_OMEGA  = num_cyclic_OMEGA;

    m_number_of_blocks = nblock;
    m_block_size       = n;
    m_border_size      = nb;
    m_num_equations    = (nblock+1)*n+q;
    n_x_n              = n*n;
    n_x_nb             = n*nb;

    integer const DE_size   { 2 * nblock * n_x_n  };
    integer const H0Nq_size { (n+q) * (2*n+q)      };
    integer const BC_size   { nb * m_num_equations };

    m_mem.reallocate( DE_size + H0Nq_size + 2*BC_size + nb*nb + num_extra_r );
    m_mem_int.reallocate( num_extra_i );
    m_DE_blk = m_mem( DE_size );
    m_H0Nq   = m_mem( H0Nq_size );
    m_Bmat   = m_mem( BC_size );
    m_Cmat   = m_mem( BC_size );
    m_Dmat   = m_mem( nb*nb );
    m_block0 = nullptr;
    m_blockN = nullptr;
  }

  /*\
   |   _                 _ ____        _   _
   |  | | ___   __ _  __| | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` |  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_|____/ \___/ \__|\__\___/|_| |_| |_|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_bottom(
    real_type const H0[], integer ld0,
    real_type const HN[], integer ldN,
    real_type const Hq[], integer ldQ
  ) {
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};

    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {
      /*\
       |  +----+-----+---+
       |  | H0 | HN  |Hq |
       |  |    |     |   |
       |  +----+-----+---+
       |  +----+------+---------+
       |  | 0  | blkN |blkN: 0  |
       |  |blk0|  0   |    :blk0|
       |  +----+------+---------+
      \*/
      integer row0  = m_num_initial_BC;
      integer rowN  = m_num_final_BC;
      integer col00 = m_num_initial_OMEGA;
      integer colNN = m_num_final_OMEGA;

      m_block0 = m_H0Nq;
      m_blockN = m_H0Nq+row0*(n+col00);

      GEcopy( rowN, n,     HN, ldN, m_blockN,        rowN );
      GEcopy( rowN, colNN, Hq, ldQ, m_blockN+n*rowN, rowN );

      GEcopy( row0, col00, Hq+rowN+colNN*ldQ, ldQ, m_block0,            row0 );
      GEcopy( row0, n,     H0+rowN,           ld0, m_block0+col00*row0, row0 );

    } else {
      integer m = n + q;
      GEcopy( m, n, H0, ld0, m_H0Nq,       m );
      GEcopy( m, n, HN, ldN, m_H0Nq+m*n,   m );
      GEcopy( m, q, Hq, ldQ, m_H0Nq+2*m*n, m );
    }

  }

  /*
  // ---------------------------------------------------------------------------
  //             col0
  //       +---+----------+
  // col00 |              |
  //       +   + - - - -  |  row0
  // row00 |   :          |
  //       |   :          |
  //       +---+----------+----------+
  //     col00 |                     |
  //           |                     | n
  //           |                     |
  //           +----------+----------+
  //                 2*n
  //
  //     +----------+----------+
  //     |                     |
  //     |                     |
  //     |                     | colNN
  //     +----------+----------+-----+
  //                |                |
  //                |                | rowN
  //                |                |
  //                |                |
  //                +----------------+
  //                      colN
  //
  */
  /*\
   |   _                 _ _____           ____        _   _
   |  | | ___   __ _  __| |_   _|__  _ __ | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` | | |/ _ \| '_ \|  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | | | (_) | |_) | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_| |_|\___/| .__/|____/ \___/ \__|\__\___/|_| |_| |_|
   |                                |_|
  \*/
  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::load_top_bottom(
    real_type const block0_in[], integer ld0,
    real_type const blockN_in[], integer ldN
  ) {
    integer const & n{m_block_size};

    UTILS_ASSERT(
      m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0,
      "in load_top_bottom num_cyclic_BC = {} and num_cyclic_OMEGA = {} "
      "must be both zero!\n",
      m_num_cyclic_BC, m_num_cyclic_OMEGA
    );

    integer row0 = m_num_initial_BC;
    integer rowN = m_num_final_BC;
    integer col0 = n + m_num_initial_OMEGA;
    integer colN = n + m_num_final_OMEGA;

    m_block0 = m_H0Nq;
    m_blockN = m_H0Nq+row0*col0;

    GEcopy( row0, col0, block0_in, ld0, m_block0, row0 );
    GEcopy( rowN, colN, blockN_in, ldN, m_blockN, rowN );

  }

  /*\
   |   _           _     _     _            _
   |  | | __ _ ___| |_  | |__ | | ___   ___| | __
   |  | |/ _` / __| __| | '_ \| |/ _ \ / __| |/ /
   |  | | (_| \__ \ |_  | |_) | | (_) | (__|   <
   |  |_|\__,_|___/\__| |_.__/|_|\___/ \___|_|\_\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::last_block_factorize() {
    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};
    integer const   nx2{n*2};

    integer m{n+q};
    integer mn{m+n};

    m_la_matrix.setup( mn, mn );
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {
      /*\
       |  +----+-----+---+
       |  | H0 | HN  |Hq |
       |  |    |     |   |
       |  +----+-----+---+
       |  +----+------+---------+
       |  | 0  | blkN |blkN: 0  |
       |  |blk0|  0   |    :blk0|
       |  +----+------+---------+
      \*/
      integer const row0{m_num_initial_BC};
      integer const rowN{m_num_final_BC};
      integer const col00{m_num_initial_OMEGA};
      integer const colNN{m_num_final_OMEGA};
      m_la_matrix.load_block( n, nx2, m_DE_blk, n );
      m_la_matrix.zero_block( m, mn, n, 0 );
      m_la_matrix.load_block( rowN, n+colNN, m_blockN,            rowN, n,      n );
      m_la_matrix.load_block( row0, n,       m_block0+col00*row0, row0, n+rowN, 0 );
      m_la_matrix.load_block( row0, col00,   m_block0,            row0, n+rowN, nx2+colNN );
    } else {
      m_la_matrix.load_block( n, nx2, m_DE_blk, n );
      m_la_matrix.load_block( m, mn,  m_H0Nq,   m, n, 0 );
      if ( m > n ) m_la_matrix.zero_block( n, q, 0, nx2 );
    }

    // fattorizzazione ultimo blocco
    m_la_factorization->factorize( "BlockBidiagonal::last_block_factorize", m_la_matrix );
  }

  /*\
   |   _                   _                   _
   |  | |__   ___  _ __ __| | ___ _ __ ___  __| |
   |  | '_ \ / _ \| '__/ _` |/ _ \ '__/ _ \/ _` |
   |  | |_) | (_) | | | (_| |  __/ | |  __/ (_| |
   |  |_.__/ \___/|_|  \__,_|\___|_|  \___|\__,_|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::factorize_bordered() {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};

    this->factorize(); // factorize top left block
    if ( nb > 0 ) {
      // Compute aux matrix
      // Z = A^(-1)*B
      // W = C*Z - D
      real_type * Zmat = m_Bmat;
      this->solve( nb, Zmat, neq );
      gemm(
        Transposition::NO,
        Transposition::NO,
        nb, nb, neq,
        1,
        m_Cmat, nb,
        Zmat, neq,
        -1,
        m_Dmat, nb
      );
      m_bb_factorization->factorize(
        "BlockBidiagonal::factorize_bordered",
        nb, nb, m_Dmat, nb
      );
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::solve_bordered( real_type xb[] ) const {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    // a' = A^(-1)*a
    this->solve( xb );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemv(
        Transposition::NO,
        nb, neq,
        1, m_Cmat, nb,
        xb, 1,
        -1, xb+neq, 1
      );
      // y = W^(-1) * b'
      m_bb_factorization->solve( xb+neq );
      // x = a' - Z*y
      real_type * Zmat = m_Bmat;
      gemv(
        Transposition::NO,
        neq, nb,
        -1, Zmat, neq,
        xb+neq, 1,
        1, xb, 1
      );
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::solve_bordered(
    integer   nrhs,
    real_type xb[],
    integer   ldRhs
  ) const {
    integer const & nb{m_border_size};
    integer const & neq{m_num_equations};
    // a' = A^(-1)*a
    this->solve( nrhs, xb, ldRhs );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemm(
        Transposition::NO,
        Transposition::NO,
        nb, nrhs, neq,
        1, m_Cmat, nb,
        xb, ldRhs,
        -1, xb+neq, ldRhs
      );
      // y = W^(-1) * b'
      m_bb_factorization->solve( nrhs, xb+neq, ldRhs );
      // x = a' - Z*y
      real_type * Zmat = m_Bmat;
      gemm(
        Transposition::NO,
        Transposition::NO,
        neq, nrhs, nb,
        -1, Zmat, neq,
        xb+neq, ldRhs,
        1, xb, ldRhs
      );
    }
  }

  /*\
   |       _                       __  __       _        _
   |    __| |_   _ _ __ ___  _ __ |  \/  | __ _| |_ _ __(_)_  __
   |   / _` | | | | '_ ` _ \| '_ \| |\/| |/ _` | __| '__| \ \/ /
   |  | (_| | |_| | | | | | | |_) | |  | | (_| | |_| |  | |>  <
   |   \__,_|\__,_|_| |_| |_| .__/|_|  |_|\__,_|\__|_|  |_/_/\_\
   |                        |_|
  \*/

  template <typename T>
  static
  void
  dumpOneMatrix(
    ostream_type & stream,
    string_view    name,
    T const        M[],
    integer        numRow,
    integer        numCol
  ) {
    fmt::print( stream,
      "# {} Size: {} x {}\n"
      "{} := <",
      name, numRow, numCol, name
    );
    for ( integer nc{0}; nc < numCol; ++nc ) {
      fmt::print( stream, "<" );
      for ( integer nr{0}; nr < numRow; ++nr ) {
        fmt::print( stream, "{}", M[ nr + nc * numRow ] );
        if ( nr+1 < numRow ) fmt::print( stream, "," );
        else                 fmt::print( stream, ">" );
      }
      if ( nc+1 < numCol ) fmt::print( stream, "|\n"  );
      else                 fmt::print( stream, ">;\n" );
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::dump_to_Maple( ostream_type & stream ) const {

    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};

    fmt::print( stream, "interface( rtablesize = 40 );\n" );
    for ( integer row{0}; row < nblock; ++row ) {
      real_type const * Ad{m_DE_blk + 2*row*n*n};
      real_type const * Au{Ad + n*n};
      dumpOneMatrix( stream, "Ad", Ad, n, n );
      dumpOneMatrix( stream, "Au", Au, n, n );
    }

    dumpOneMatrix( stream, "H0Nq", m_H0Nq, n+q, 2*n+q );

    fmt::print( stream,
      "with(LinearAlgebra):\n"
      "Determinant(Ad);\n"
      "Determinant(Au);\n"
      "Rank(H0Np);\n"
      "Rank(<H0N|Hp>);\n"
    );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::Mv(
    real_type const x[],
    real_type       res[]
  ) const {

    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & q      {m_extra_bc};
    integer const & nb     {m_border_size};
    integer const & neq    {m_num_equations};
    integer const   nx2    {n*2};

    alglin::Zero_n( res, neq+nb );
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {
      integer row0  {m_num_initial_BC};
      integer rowN  {m_num_final_BC};
      integer col00 {m_num_initial_OMEGA};
      integer colNN {m_num_final_OMEGA};

      real_type const * xe{x+neq-(n+colNN+col00)};
      real_type       * rese{res+neq-(row0+rowN)};

      gemv(
        Transposition::NO, rowN, n+colNN,
        1.0, m_blockN, rowN,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        Transposition::NO, row0, n,
        1.0, m_block0+col00*row0, row0,
        x, 1,
        1.0, rese+rowN, 1
      );

      gemv(
        Transposition::NO, row0, col00,
        1.0, m_block0, row0,
        xe+n+colNN, 1,
        1.0, rese+rowN, 1
      );
    } else {
      integer m{n+q};
      real_type const * xe   { x+neq-(n+m_num_initial_OMEGA+m_num_final_OMEGA+m_num_cyclic_OMEGA) };
      real_type       * rese { res+neq-(m_num_initial_BC+m_num_final_BC+m_num_cyclic_BC) };

      gemv(
        Transposition::NO, m, n,
        1.0, m_H0Nq, m,
        x, 1,
        1.0, rese, 1
      );

      gemv(
        Transposition::NO, m, n,
        1.0, m_H0Nq+m*n, m,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        Transposition::NO, m, q,
        1.0, m_H0Nq+m*nx2, m,
        xe+n, 1,
        1.0, rese, 1
      );
    }

    // internal blocks block
    t_Value const * xx{x};
    t_Value *       yy{res};
    t_Value const * DE{m_DE_blk};
    for ( integer i{0}; i < nblock; ++i ) {
      gemv(
        Transposition::NO, n, nx2,
        1.0, DE, n,
        xx, 1,
        1.0, yy, 1
      );
      xx += n;
      yy += n;
      DE += 2*n_x_n;
    }
    if ( nb > 0 ) {
      gemv(
        Transposition::NO, neq, nb,
        1.0, m_Bmat, neq,
        x+neq, 1,
        1.0, res, 1
      );
      gemv(
        Transposition::NO, nb, neq,
        1.0, m_Cmat, nb,
        x, 1,
        1.0, res+neq, 1
      );
      gemv(
        Transposition::NO, nb, nb,
        1.0, m_Dmat, nb,
        x+neq, 1,
        1.0, res+neq, 1
      );
    }
  }

  /*\
   |
   |       _
   |    __| |_   _ _ __ ___  _ __     ___ ___ ___   ___  _ __
   |   / _` | | | | '_ ` _ \| '_ \   / __/ __/ _ \ / _ \| '__|
   |  | (_| | |_| | | | | | | |_) | | (_| (_| (_) | (_) | |
   |   \__,_|\__,_|_| |_| |_| .__/___\___\___\___/ \___/|_|
   |                        |_| |_____|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::dump_ccoord( ostream_type & stream ) const {
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & q      {m_extra_bc};
    integer const & nb     {m_border_size};
    integer const   nx2    {n*2};

    integer nnz{ (2*nblock)*n_x_n+ 2*(nblock+1)*n*nb + nb*nb };

    // BC
    integer ii;
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {

      integer const row0  = m_num_initial_BC;
      integer const rowN  = m_num_final_BC;
      integer const col00 = m_num_initial_OMEGA;
      integer const colNN = m_num_final_OMEGA;

      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
      fmt::print( stream, "{}\n", nnz );

      ii = nblock*n;
      for ( integer i{0}; i < rowN; ++i )
        for ( integer j{0}; j < n+colNN; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n", ii+i, ii+j,
            m_blockN[i+j*rowN]
          );

      for ( integer i{0}; i < row0; ++i )
        for ( integer j{0}; j < n; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+rowN+i, j, m_block0[i+(j+col00)*row0]
          );

      for ( integer i{0}; i < row0; ++i )
        for ( integer j{0}; j < col00; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+rowN+i, ii+n+colNN+j, m_block0[i+j*row0]
          );

    } else {

      integer const nq{n+q};
      nnz += nq*(2*nq);
      fmt::print( stream, "{}\n", nnz );

      real_type * H0{m_H0Nq};
      ii = nblock*n;
      for ( integer i{0}; i < nq; ++i )
        for ( integer j{0}; j < n; ++j )
          fmt::print( stream,
            "{}\t{}\t{}\n",
            ii+i, j, H0[i+j*nq]
          );

      real_type * HNq{m_H0Nq+n*nq};
      for ( integer i{0}; i < nq; ++i )
        for ( integer j{0}; j < nq; ++j )
          fmt::print( stream,
            "{}\t{}\t{}\n",
            ii+i, ii+j, HNq[i+j*nq]
          );
    }

    // bidiagonal
    for ( integer k{0}; k < nblock; ++k ) {
      ii = k*n;
      real_type * DE{ m_DE_blk + (2*k) * n_x_n };
      for ( integer i{0}; i < n; ++i )
        for ( integer j{0}; j < nx2; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+i, ii+j, DE[i+j*n]
          );
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i{0}; i < nb; ++i )
      for ( integer j{0}; j < ii; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n{}\t{}\t{}\n",
            ii+i, j, m_Cmat[i+j*nb],
            j, ii+i, m_Bmat[j+i*ii]
          );

    for ( integer i{0}; i < nb; ++i )
      for ( integer j{0}; j < nb; ++j )
        fmt::print(
          stream, "{}\t{}\t{}\n",
          ii+i, ii+j, m_Dmat[i+j*nb]
        );

  }

  /*\
   |   ___ _ __   __ _ _ __ ___  ___
   |  / __| '_ \ / _` | '__/ __|/ _ \
   |  \__ \ |_) | (_| | |  \__ \  __/
   |  |___/ .__/ \__,_|_|  |___/\___|
   |      |_|
  \*/

  template <typename t_Value>
  integer
  BlockBidiagonal<t_Value>::sparse_nnz() const {
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & q      {m_extra_bc};
    integer const & nb     {m_border_size};

    integer nnz{ (2*nblock)*n_x_n + 2*(nblock+1)*n*nb + nb*nb };
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {
      integer const row0  {m_num_initial_BC};
      integer const rowN  {m_num_final_BC};
      integer const col00 {m_num_initial_OMEGA};
      integer const colNN {m_num_final_OMEGA};
      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
    } else {
      nnz += (n+q)*(2*n+q);
    }
    return nnz;
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparse_pattern( integer I[], integer J[] ) const {
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & q      {m_extra_bc};
    integer const & nb     {m_border_size};
    integer const   nx2    {n*2};
    integer kkk = 0;

    // BC
    integer ii;
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {

      integer const row0  {m_num_initial_BC};
      integer const rowN  {m_num_final_BC};
      integer const col00 {m_num_initial_OMEGA};
      integer const colNN {m_num_final_OMEGA};

      ii = nblock*n;
      for ( integer i{0}; i < rowN; ++i ) {
        for ( integer j{0}; j < n+colNN; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }

      for ( integer i{0}; i < row0; ++i ) {
        for ( integer j{0}; j < n; ++j ) {
          I[kkk] = ii+rowN+i;
          J[kkk] = j;
          ++kkk;
        }
      }

      for ( integer i{0}; i < row0; ++i ) {
        for ( integer j{0}; j < col00; ++j ) {
          I[kkk] = ii+rowN+i;
          J[kkk] = ii+n+colNN+j;
          ++kkk;
        }
      }

    } else {

      integer const nq{n + q};

      ii = nblock*n;
      for ( integer i{0}; i < nq; ++i ) {
        for ( integer j{0}; j < n; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = j;
          ++kkk;
        }
      }

      for ( integer i{0}; i < nq; ++i ) {
        for ( integer j{0}; j < nq; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }
    }

    // bidiagonal
    for ( integer k{0}; k < nblock; ++k ) {
      ii = k*n;
      for ( integer i{0}; i < n; ++i ) {
        for ( integer j{0}; j < nx2; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i{0}; i < nb; ++i ) {
      for ( integer j{0}; j < ii; ++j ) {
        I[kkk] = ii+i;
        J[kkk] = j;
        ++kkk;
        I[kkk] = j;
        J[kkk] = ii+i;
        ++kkk;
      }
    }

    for ( integer i{0}; i < nb; ++i ) {
      for ( integer j{0}; j < nb; ++j ) {
        I[kkk] = ii+i;
        J[kkk] = ii+j;
        ++kkk;
      }
    }

    UTILS_ASSERT(
      kkk == sparse_nnz(),
      "BlockBidiagonal::sparse_pattern( V ), inserted {} values, expected {}\n",
      kkk, sparse_nnz()
    );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparse_values( real_type V[] ) const {
    integer const & nblock {m_number_of_blocks};
    integer const & n      {m_block_size};
    integer const & q      {m_extra_bc};
    integer const & nb     {m_border_size};
    integer const   nx2    {n*2};
    integer kkk{0};

    // BC
    integer ii;
    if ( m_num_cyclic_BC == 0 && m_num_cyclic_OMEGA == 0 ) {

      integer const row0  { m_num_initial_BC };
      integer const rowN  { m_num_final_BC };
      integer const col00 { m_num_initial_OMEGA };
      integer const colNN { m_num_final_OMEGA };

      ii = nblock*n;
      for ( integer i{0}; i < rowN; ++i )
        for ( integer j{0}; j < n+colNN; ++j )
          V[kkk++] = m_blockN[i+j*rowN];

      for ( integer i{0}; i < row0; ++i )
        for ( integer j{0}; j < n; ++j )
          V[kkk++] = m_block0[i+(j+col00)*row0];

      for ( integer i{0}; i < row0; ++i )
        for ( integer j{0}; j < col00; ++j )
          V[kkk++] = m_block0[i+j*row0];

    } else {

      integer const nq{n + q};

      real_type * H0{m_H0Nq};
      ii = nblock*n;
      for ( integer i{0}; i < nq; ++i )
        for ( integer j{0}; j < n; ++j )
          V[kkk++] = H0[i+j*nq];

      real_type * HNq{m_H0Nq+n*nq};
      for ( integer i{0}; i < nq; ++i )
        for ( integer j{0}; j < nq; ++j )
          V[kkk++] = HNq[i+j*nq];
    }

    // bidiagonal
    for ( integer k{0}; k < nblock; ++k ) {
      ii = k*n;
      real_type * DE{m_DE_blk + (2*k) * n_x_n};
      for ( integer i{0}; i < n; ++i )
        for ( integer j{0}; j < nx2; ++j )
          V[kkk++] = DE[i+j*n];
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i{0}; i < nb; ++i ) {
      for ( integer j{0}; j < ii; ++j ) {
        V[kkk++] = m_Cmat[i+j*nb];
        V[kkk++] = m_Bmat[j+i*ii];
      }
    }

    for ( integer i{0}; i < nb; ++i )
      for ( integer j{0}; j < nb; ++j )
        V[kkk++] = m_Dmat[i+j*nb];

    UTILS_ASSERT(
      kkk == sparse_nnz(),
      "BlockBidiagonal::sparse_values( V ), inserted {} values, expected {}\n",
      kkk, sparse_nnz()
    );

  }

  template class BlockBidiagonal<float>;
  template class BlockBidiagonal<double>;

}
