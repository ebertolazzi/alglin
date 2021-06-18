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
  BlockBidiagonal<t_Value>::allocate(
    integer nblock,
    integer n,
    integer nb,
    // ----------------------
    integer numInitialBC,
    integer numFinalBC,
    integer numCyclicBC,
    // ----------------------
    integer numInitialOMEGA,
    integer numFinalOMEGA,
    integer numCyclicOMEGA,
    // ----------------------
    integer num_extra_r,
    integer num_extra_i
  ) {

    UTILS_ASSERT(
      numInitialBC  >= 0 && numFinalBC      >= 0 &&
      numCyclicBC   >= 0 && numInitialOMEGA >= 0 &&
      numFinalOMEGA >= 0 && numCyclicOMEGA  >= 0,
      "Bad BC specification:"
      "\nnumInitialBC    = {}"
      "\nnumFinalBC      = {}"
      "\nnumCyclicBC     = {}"
      "\nnumInitialOMEGA = {}"
      "\nnumFinalOMEGA   = {}"
      "\nnumCyclicOMEGA  = {}\n",
      numInitialBC,    numFinalBC,    numCyclicBC,
      numInitialOMEGA, numFinalOMEGA, numCyclicOMEGA
    );

    m_extra_bc = numInitialOMEGA + numFinalOMEGA + numCyclicOMEGA;
    integer const & q = m_extra_bc;

    UTILS_ASSERT(
      numInitialBC + numFinalBC + numCyclicBC == n+q,
      "Bad BC specification:"
      "\nnumInitialBC    = {}"
      "\nnumFinalBC      = {}"
      "\nnumCyclicBC     = {}"
      "\nnumInitialOMEGA + numFinalOMEGA + numCyclicOMEGA must be = {}\n",
      numInitialBC, numFinalBC, numCyclicBC, n+q
    );

    m_numInitialBC    = numInitialBC;
    m_numFinalBC      = numFinalBC;
    m_numCyclicBC     = numCyclicBC;
    m_numInitialOMEGA = numInitialOMEGA;
    m_numFinalOMEGA   = numFinalOMEGA;
    m_numCyclicOMEGA  = numCyclicOMEGA;

    m_number_of_blocks = nblock;
    m_block_size       = n;
    m_border_size      = nb;
    m_num_equations    = (nblock+1)*n+q;
    n_x_n              = n*n;
    n_x_nb             = n*nb;
    integer DE_size    = 2*nblock*n_x_n;
    integer H0Nq_size  = (n+q)*(2*n+q);
    integer BC_size    = nb*m_num_equations;
    m_baseValue.reallocate( size_t(DE_size+H0Nq_size+2*BC_size+nb*nb+num_extra_r) );
    m_baseInteger.reallocate( size_t(num_extra_i) );
    m_DE_blk = m_baseValue( size_t(DE_size) );
    m_H0Nq   = m_baseValue( size_t(H0Nq_size) );
    m_Bmat   = m_baseValue( size_t(BC_size) );
    m_Cmat   = m_baseValue( size_t(BC_size) );
    m_Dmat   = m_baseValue( size_t(nb*nb) );
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
  BlockBidiagonal<t_Value>::loadBottom(
    real_type const H0[], integer ld0,
    real_type const HN[], integer ldN,
    real_type const Hq[], integer ldQ
  ) {
    integer const & n = m_block_size;
    integer const & q = m_extra_bc;

    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {
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
      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;

      m_block0 = m_H0Nq;
      m_blockN = m_H0Nq+row0*(n+col00);

      gecopy( rowN, n,     HN, ldN, m_blockN,        rowN );
      gecopy( rowN, colNN, Hq, ldQ, m_blockN+n*rowN, rowN );

      gecopy( row0, col00, Hq+rowN+colNN*ldQ, ldQ, m_block0,            row0 );
      gecopy( row0, n,     H0+rowN,           ld0, m_block0+col00*row0, row0 );

    } else {
      integer m = n + q;
      gecopy( m, n, H0, ld0, m_H0Nq,       m );
      gecopy( m, n, HN, ldN, m_H0Nq+m*n,   m );
      gecopy( m, q, Hq, ldQ, m_H0Nq+2*m*n, m );
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
  BlockBidiagonal<t_Value>::loadTopBottom(
    real_type const block0_in[], integer ld0,
    real_type const blockN_in[], integer ldN
  ) {
    integer const & n = m_block_size;

    UTILS_ASSERT(
      m_numCyclicBC == 0 && m_numCyclicOMEGA == 0,
      "in loadTopBottom numCyclicBC = {} and numCyclicOMEGA = {} "
      "must be both zero!\n",
      m_numCyclicBC, m_numCyclicOMEGA
    );

    integer row0 = m_numInitialBC;
    integer rowN = m_numFinalBC;
    integer col0 = n + m_numInitialOMEGA;
    integer colN = n + m_numFinalOMEGA;

    m_block0 = m_H0Nq;
    m_blockN = m_H0Nq+row0*col0;

    gecopy( row0, col0, block0_in, ld0, m_block0, row0 );
    gecopy( rowN, colN, blockN_in, ldN, m_blockN, rowN );

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
    integer const & n   = m_block_size;
    integer const & q   = m_extra_bc;
    integer const   nx2 = n*2;

    integer m  = n+q;
    integer mn = m+n;

    m_la_matrix.setup( mn, mn );
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {
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
      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;
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
    integer const & nb  = m_border_size;
    integer const & neq = m_num_equations;

    this->factorize(); // factorize top left block
    if ( nb > 0 ) {
      // Compute aux matrix
      // Z = A^(-1)*B
      // W = C*Z - D
      real_type * Zmat = m_Bmat;
      this->solve( nb, Zmat, neq );
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
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
    integer const & nb  = m_border_size;
    integer const & neq = m_num_equations;
    // a' = A^(-1)*a
    this->solve( xb );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemv(
        NO_TRANSPOSE,
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
        NO_TRANSPOSE,
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
    integer const & nb  = m_border_size;
    integer const & neq = m_num_equations;
    // a' = A^(-1)*a
    this->solve( nrhs, xb, ldRhs );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
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
        NO_TRANSPOSE,
        NO_TRANSPOSE,
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
    char const *   name,
    T const        M[],
    integer        numRow,
    integer        numCol
  ) {
    stream
      << "# " << name
      << " Size: " << numRow << " x " << numCol << '\n'
      << name << " := <";
    for ( integer nc = 0; nc < numCol; ++nc ) {
      stream << '<';
      for ( integer nr = 0; nr < numRow; ++nr ) {
        stream << M[ nr + nc * numRow ];
        if ( nr+1 < numRow ) stream << ',';
        else                 stream << '>';
      }
      if ( nc+1 < numCol ) stream << "|\n";
      else                 stream << ">;\n";
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::dump_to_Maple( ostream_type & stream ) const {

    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;

    stream << "interface( rtablesize = 40 );\n";
    for ( integer row = 0; row < nblock; ++row ) {
      real_type const * Ad = m_DE_blk + 2*row*n*n;
      real_type const * Au = Ad + n*n;
      dumpOneMatrix( stream, "Ad", Ad, n, n );
      dumpOneMatrix( stream, "Au", Au, n, n );
    }

    dumpOneMatrix( stream, "H0Nq", m_H0Nq, n+q, 2*n+q );

    stream << "with(LinearAlgebra):\n";
    stream << "Determinant(Ad);\n";
    stream << "Determinant(Au);\n";
    stream << "Rank(H0Np);\n";
    stream << "Rank(<H0N|Hp>);\n";

  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::Mv(
    real_type const x[],
    real_type       res[]
  ) const {

    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;
    integer const & nb     = m_border_size;
    integer const & neq    = m_num_equations;
    integer const   nx2    = n*2;

    alglin::zero( neq+nb, res, 1 );
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {
      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;

      real_type const * xe   = x+neq-(n+colNN+col00);
      real_type       * rese = res+neq-(row0+rowN);

      gemv(
        NO_TRANSPOSE, rowN, n+colNN,
        1.0, m_blockN, rowN,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, row0, n,
        1.0, m_block0+col00*row0, row0,
        x, 1,
        1.0, rese+rowN, 1
      );

      gemv(
        NO_TRANSPOSE, row0, col00,
        1.0, m_block0, row0,
        xe+n+colNN, 1,
        1.0, rese+rowN, 1
      );
    } else {
      integer m = n+q;
      real_type const * xe   = x+neq-(n+m_numInitialOMEGA+m_numFinalOMEGA+m_numCyclicOMEGA);
      real_type       * rese = res+neq-(m_numInitialBC+m_numFinalBC+m_numCyclicBC);

      gemv(
        NO_TRANSPOSE, m, n,
        1.0, m_H0Nq, m,
        x, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, m, n,
        1.0, m_H0Nq+m*n, m,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, m, q,
        1.0, m_H0Nq+m*nx2, m,
        xe+n, 1,
        1.0, rese, 1
      );
    }

    // internal blocks block
    t_Value const * xx = x;
    t_Value *       yy = res;
    t_Value const * DE = m_DE_blk;
    for ( integer i = 0; i < nblock; ++i ) {
      gemv(
        NO_TRANSPOSE, n, nx2,
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
        NO_TRANSPOSE, neq, nb,
        1.0, m_Bmat, neq,
        x+neq, 1,
        1.0, res, 1
      );
      gemv(
        NO_TRANSPOSE, nb, neq,
        1.0, m_Cmat, nb,
        x, 1,
        1.0, res+neq, 1
      );
      gemv(
        NO_TRANSPOSE, nb, nb,
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
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;
    integer const & nb     = m_border_size;
    integer const   nx2    = n*2;

    integer nnz = (2*nblock)*n_x_n+ 2*(nblock+1)*n*nb + nb*nb;

    // BC
    integer ii;
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {

      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;

      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
      fmt::print( stream, "{}\n", nnz );

      ii = nblock*n;
      for ( integer i = 0; i < rowN; ++i )
        for ( integer j = 0; j < n+colNN; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n", ii+i, ii+j,
            m_blockN[i+j*rowN]
          );

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < n; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+rowN+i, j, m_block0[i+(j+col00)*row0]
          );

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < col00; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+rowN+i, ii+n+colNN+j, m_block0[i+j*row0]
          );

    } else {

      integer nq = n+q;
      nnz += nq*(2*nq);
      stream << nnz << '\n';

      real_type * H0 = m_H0Nq;
      ii = nblock*n;
      for ( integer i = 0; i < nq; ++i )
        for ( integer j = 0; j < n; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+i, j, H0[i+j*nq]
          );

      real_type * HNq = m_H0Nq+n*nq;
      for ( integer i = 0; i < nq; ++i )
        for ( integer j = 0; j < nq; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+i, ii+j, HNq[i+j*nq]
          );
    }

    // bidiagonal
    for ( integer k = 0; k < nblock; ++k ) {
      ii = k*n;
      real_type * DE = m_DE_blk + (2*k) * n_x_n;
      for ( integer i = 0; i < n; ++i )
        for ( integer j = 0; j < nx2; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n",
            ii+i, ii+j, DE[i+j*n]
          );
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < ii; ++j )
          fmt::print(
            stream, "{}\t{}\t{}\n{}\t{}\t{}\n",
            ii+i, j, m_Cmat[i+j*nb],
            j, ii+i, m_Bmat[j+i*ii]
          );

    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < nb; ++j )
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
  BlockBidiagonal<t_Value>::sparseNnz() const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;
    integer const & nb     = m_border_size;

    integer nnz = (2*nblock)*n_x_n + 2*(nblock+1)*n*nb + nb*nb;
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {
      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;
      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
    } else {
      nnz += (n+q)*(2*n+q);
    }
    return nnz;
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparsePattern( integer I[], integer J[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;
    integer const & nb     = m_border_size;
    integer const   nx2    = n*2;
    integer kkk = 0;

    // BC
    integer ii;
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {

      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;

      ii = nblock*n;
      for ( integer i = 0; i < rowN; ++i ) {
        for ( integer j = 0; j < n+colNN; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }

      for ( integer i = 0; i < row0; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          I[kkk] = ii+rowN+i;
          J[kkk] = j;
          ++kkk;
        }
      }

      for ( integer i = 0; i < row0; ++i ) {
        for ( integer j = 0; j < col00; ++j ) {
          I[kkk] = ii+rowN+i;
          J[kkk] = ii+n+colNN+j;
          ++kkk;
        }
      }

    } else {

      integer nq = n + q;

      ii = nblock*n;
      for ( integer i = 0; i < nq; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = j;
          ++kkk;
        }
      }

      for ( integer i = 0; i < nq; ++i ) {
        for ( integer j = 0; j < nq; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }
    }

    // bidiagonal
    for ( integer k = 0; k < nblock; ++k ) {
      ii = k*n;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < nx2; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = ii+j;
          ++kkk;
        }
      }
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i = 0; i < nb; ++i ) {
      for ( integer j = 0; j < ii; ++j ) {
        I[kkk] = ii+i;
        J[kkk] = j;
        ++kkk;
        I[kkk] = j;
        J[kkk] = ii+i;
        ++kkk;
      }
    }

    for ( integer i = 0; i < nb; ++i ) {
      for ( integer j = 0; j < nb; ++j ) {
        I[kkk] = ii+i;
        J[kkk] = ii+j;
        ++kkk;
      }
    }

    UTILS_ASSERT(
      kkk == sparseNnz(),
      "BlockBidiagonal::sparsePattern( V ), inserted {} values, expected {}\n",
      kkk, sparseNnz()
    );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparseValues( real_type V[] ) const {
    integer const & nblock = m_number_of_blocks;
    integer const & n      = m_block_size;
    integer const & q      = m_extra_bc;
    integer const & nb     = m_border_size;
    integer const   nx2    = n*2;
    integer kkk = 0;

    // BC
    integer ii;
    if ( m_numCyclicBC == 0 && m_numCyclicOMEGA == 0 ) {

      integer row0  = m_numInitialBC;
      integer rowN  = m_numFinalBC;
      integer col00 = m_numInitialOMEGA;
      integer colNN = m_numFinalOMEGA;

      ii = nblock*n;
      for ( integer i = 0; i < rowN; ++i )
        for ( integer j = 0; j < n+colNN; ++j )
          V[kkk++] = m_blockN[i+j*rowN];

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < n; ++j )
          V[kkk++] = m_block0[i+(j+col00)*row0];

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < col00; ++j )
          V[kkk++] = m_block0[i+j*row0];

    } else {

      integer nq = n + q;

      real_type * H0 = m_H0Nq;
      ii = nblock*n;
      for ( integer i = 0; i < nq; ++i )
        for ( integer j = 0; j < n; ++j )
          V[kkk++] = H0[i+j*nq];

      real_type * HNq = m_H0Nq+n*nq;
      for ( integer i = 0; i < nq; ++i )
        for ( integer j = 0; j < nq; ++j )
          V[kkk++] = HNq[i+j*nq];
    }

    // bidiagonal
    for ( integer k = 0; k < nblock; ++k ) {
      ii = k*n;
      real_type * DE = m_DE_blk + (2*k) * n_x_n;
      for ( integer i = 0; i < n; ++i )
        for ( integer j = 0; j < nx2; ++j )
          V[kkk++] = DE[i+j*n];
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i = 0; i < nb; ++i ) {
      for ( integer j = 0; j < ii; ++j ) {
        V[kkk++] = m_Cmat[i+j*nb];
        V[kkk++] = m_Bmat[j+i*ii];
      }
    }

    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < nb; ++j )
        V[kkk++] = m_Dmat[i+j*nb];

    UTILS_ASSERT(
      kkk == sparseNnz(),
      "BlockBidiagonal::sparseValues( V ), inserted {} values, expected {}\n",
      kkk, sparseNnz()
    );

  }

  template class BlockBidiagonal<float>;
  template class BlockBidiagonal<double>;

}
