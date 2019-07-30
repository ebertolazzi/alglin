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

#include "BlockBidiagonal.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace alglin {

  std::string
  LastBlock_to_string( LASTBLOCK_Choice c ) {
    switch ( c ) {
      case LASTBLOCK_LU:  return "last block LU";
      case LASTBLOCK_QR:  return "last block QR";
      case LASTBLOCK_QRP: return "last block QRP";
      case LASTBLOCK_SVD: return "last block SVD";
    }
    return "last block not selected";
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
    integer _nblock,
    integer _n,
    integer _nb,
    // ----------------------
    integer _numInitialBC,
    integer _numFinalBC,
    integer _numCyclicBC,
    // ----------------------
    integer _numInitialOMEGA,
    integer _numFinalOMEGA,
    integer _numCyclicOMEGA,
    // ----------------------
    integer num_extra_r,
    integer num_extra_i
  ) {

    ALGLIN_ASSERT(
      _numInitialBC  >= 0 && _numFinalBC      >= 0 &&
      _numCyclicBC   >= 0 && _numInitialOMEGA >= 0 &&
      _numFinalOMEGA >= 0 && _numCyclicOMEGA  >= 0,
      "Bad BC specification:" <<
      "\nnumInitialBC    = " << _numInitialBC <<
      "\nnumFinalBC      = " << _numFinalBC <<
      "\nnumCyclicBC     = " << _numCyclicBC <<
      "\nnumInitialOMEGA = " << _numInitialOMEGA <<
      "\nnumFinalOMEGA   = " << _numFinalOMEGA <<
      "\nnumCyclicOMEGA  = " << _numCyclicOMEGA
    );

    q = _numInitialOMEGA + _numFinalOMEGA + _numCyclicOMEGA;

    ALGLIN_ASSERT(
      _numInitialBC + _numFinalBC + _numCyclicBC == _n+q,
      "Bad BC specification:" <<
      "\nnumInitialBC    = " << _numInitialBC <<
      "\nnumFinalBC      = " << _numFinalBC   <<
      "\nnumCyclicBC     = " << _numCyclicBC  <<
      "\nnumInitialOMEGA + numFinalOMEGA + numCyclicOMEGA must be = " << _n+q
    );

    numInitialBC    = _numInitialBC;
    numFinalBC      = _numFinalBC;
    numCyclicBC     = _numCyclicBC;
    numInitialOMEGA = _numInitialOMEGA;
    numFinalOMEGA   = _numFinalOMEGA;
    numCyclicOMEGA  = _numCyclicOMEGA;

    nblock = _nblock;
    n      = _n;
    nb     = _nb;
    neq    = (nblock+1)*n+q;
    nx2    = n*2;
    nxn    = n*n;
    nxnx2  = nxn*2;
    nxnb   = n*nb;
    integer DE_size   = nblock*nxnx2;
    integer H0Nq_size = (n+q)*(nx2+q);
    integer BC_size   = nb*neq;
    baseValue.allocate(size_t(DE_size+H0Nq_size+2*BC_size+nb*nb+num_extra_r));
    baseInteger.allocate(size_t(num_extra_i));
    DE_blk = baseValue(size_t(DE_size));
    H0Nq   = baseValue(size_t(H0Nq_size));
    Bmat   = baseValue(size_t(BC_size));
    Cmat   = baseValue(size_t(BC_size));
    Dmat   = baseValue(size_t(nb*nb));
    block0 = nullptr;
    blockN = nullptr;
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
    valueType const H0[], integer ld0,
    valueType const HN[], integer ldN,
    valueType const Hq[], integer ldQ
  ) {

    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {
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
      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;

      block0 = H0Nq;
      blockN = H0Nq+row0*(n+col00);

      gecopy( rowN, n,     HN, ldN, blockN,        rowN );
      gecopy( rowN, colNN, Hq, ldQ, blockN+n*rowN, rowN );

      gecopy( row0, col00, Hq+rowN+colNN*ldQ, ldQ, block0,            row0 );
      gecopy( row0, n,     H0+rowN,           ld0, block0+col00*row0, row0 );

    } else {
      integer m = n + q;
      gecopy( m, n, H0, ld0, H0Nq,       m );
      gecopy( m, n, HN, ldN, H0Nq+m*n,   m );
      gecopy( m, q, Hq, ldQ, H0Nq+2*m*n, m );
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
    valueType const block0_in[], integer ld0,
    valueType const blockN_in[], integer ldN
  ) {

    ALGLIN_ASSERT(
      numCyclicBC == 0 && numCyclicOMEGA == 0,
      "in loadTopBottom numCyclicBC = " << numCyclicBC <<
      " and numCyclicOMEGA = " << numCyclicOMEGA << " must be both zero!"
    );

    integer row0 = numInitialBC;
    integer rowN = numFinalBC;
    integer col0 = n + numInitialOMEGA;
    integer colN = n + numFinalOMEGA;

    block0 = H0Nq;
    blockN = H0Nq+row0*col0;

    gecopy( row0, col0, block0_in, ld0, block0, row0 );
    gecopy( rowN, colN, blockN_in, ldN, blockN, rowN );

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
    integer m  = n+q;
    integer mn = m+n;

    la_factorization->allocate( mn, mn );
    if ( this->numCyclicBC == 0 && this->numCyclicOMEGA == 0 ) {
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
      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;
      la_factorization->load_block( n, nx2, DE_blk, n );
      la_factorization->zero_block( m, mn, n, 0 );
      la_factorization->load_block( rowN, n+colNN, blockN, rowN, n, n );
      la_factorization->load_block( row0, n,     block0+col00*row0, row0, n+rowN, 0 );
      la_factorization->load_block( row0, col00, block0,            row0, n+rowN, nx2+colNN );
    } else {
      la_factorization->load_block( n, nx2, DE_blk, n );
      la_factorization->load_block( m, mn, H0Nq, m, n, 0 );
      if ( m > n ) la_factorization->zero_block( n, q, 0, nx2 );
    }

    // fattorizzazione ultimo blocco
    this->la_factorization->factorize( "BlockBidiagonal::last_block_factorize" );
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

    this->factorize(); // factorize top left block
    if ( nb > 0 ) {
      // Compute aux matrix
      // Z = A^(-1)*B
      // W = C*Z - D
      valueType * Zmat = Bmat;
      this->solve( nb, Zmat, neq );
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
        nb, nb, neq,
        1,
        Cmat, nb,
        Zmat, neq,
        -1,
        Dmat, nb
      );
      bb_factorization->factorize(
        "BlockBidiagonal::factorize_bordered",
        nb, nb, Dmat, nb
      );
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::solve_bordered( valueType xb[] ) const {
    // a' = A^(-1)*a
    this->solve( xb );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemv(
        NO_TRANSPOSE,
        nb, neq,
        1, Cmat, nb,
        xb, 1,
        -1, xb+neq, 1
      );
      // y = W^(-1) * b'
      bb_factorization->solve( xb+neq );
      // x = a' - Z*y
      valueType * Zmat = this->Bmat;
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
    valueType xb[],
    integer   ldRhs
  ) const {
    // a' = A^(-1)*a
    this->solve( nrhs, xb, ldRhs );
    if ( nb > 0 ) {
      // b' = C*a' - b
      gemm(
        NO_TRANSPOSE,
        NO_TRANSPOSE,
        nb, nrhs, neq,
        1, Cmat, nb,
        xb, ldRhs,
        -1, xb+neq, ldRhs
      );
      // y = W^(-1) * b'
      bb_factorization->solve( nrhs, xb+neq, ldRhs );
      // x = a' - Z*y
      valueType * Zmat = this->Bmat;
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

    stream << "interface( rtablesize = 40 );\n";
    for ( integer row = 0; row < nblock; ++row ) {
      valueType const * Ad = DE_blk + 2*row*n*n;
      valueType const * Au = Ad + n*n;
      dumpOneMatrix( stream, "Ad", Ad, n, n );
      dumpOneMatrix( stream, "Au", Au, n, n );
    }

    dumpOneMatrix( stream, "H0Nq", H0Nq, n+q, 2*n+q );

    stream << "with(LinearAlgebra):\n";
    stream << "Determinant(Ad);\n";
    stream << "Determinant(Au);\n";
    stream << "Rank(H0Np);\n";
    stream << "Rank(<H0N|Hp>);\n";

  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::Mv(
    valueType const x[],
    valueType       res[]
  ) const {
    alglin::zero( neq+nb, res, 1 );
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {
      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;

      valueType const * xe   = x+neq-(n+colNN+col00);
      valueType       * rese = res+neq-(row0+rowN);

      gemv(
        NO_TRANSPOSE, rowN, n+colNN,
        1.0, blockN, rowN,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, row0, n,
        1.0, block0+col00*row0, row0,
        x, 1,
        1.0, rese+rowN, 1
      );

      gemv(
        NO_TRANSPOSE, row0, col00,
        1.0, block0, row0,
        xe+n+colNN, 1,
        1.0, rese+rowN, 1
      );
    } else {
      integer m = n+q;
      valueType const * xe   = x+neq-(n+numInitialOMEGA+numFinalOMEGA+numCyclicOMEGA);
      valueType       * rese = res+neq-(numInitialBC+numFinalBC+numCyclicBC);

      gemv(
        NO_TRANSPOSE, m, n,
        1.0, H0Nq, m,
        x, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, m, n,
        1.0, H0Nq+m*n, m,
        xe, 1,
        1.0, rese, 1
      );

      gemv(
        NO_TRANSPOSE, m, q,
        1.0, H0Nq+m*nx2, m,
        xe+n, 1,
        1.0, rese, 1
      );
    }

    // internal blocks block
    t_Value const * xx = x;
    t_Value *       yy = res;
    t_Value const * DE = DE_blk;
    for ( integer i = 0; i < nblock; ++i ) {
      gemv(
        NO_TRANSPOSE, n, nx2,
        1.0, DE, n,
        xx, 1,
        1.0, yy, 1
      );
      xx += n;
      yy += n;
      DE += nxnx2;
    }
    if ( nb > 0 ) {
      gemv(
        NO_TRANSPOSE, neq, nb,
        1.0, Bmat, neq,
        x+neq, 1,
        1.0, res, 1
      );
      gemv(
        NO_TRANSPOSE, nb, neq,
        1.0, Cmat, nb,
        x, 1,
        1.0, res+neq, 1
      );
      gemv(
        NO_TRANSPOSE, nb, nb,
        1.0, Dmat, nb,
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
    integer nnz = nblock*nxnx2 + 2*(nblock+1)*n*nb + nb*nb;

    // BC
    integer ii;
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {

      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;

      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
      stream << nnz << '\n';

      ii = nblock*n;
      for ( integer i = 0; i < rowN; ++i )
        for ( integer j = 0; j < n+colNN; ++j )
          stream
            << ii+i << '\t'
            << ii+j << '\t'
            << blockN[i+j*rowN] << '\n';

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < n; ++j )
          stream
            << ii+rowN+i << '\t'
            << j         << '\t'
            << block0[i+(j+col00)*row0] << '\n';

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < col00; ++j )
          stream
            << ii+rowN+i    << '\t'
            << ii+n+colNN+j << '\t'
            << block0[i+j*row0] << '\n';

    } else {

      nnz += (n+q)*(2*n+q);
      stream << nnz << '\n';

      valueType * H0 = H0Nq;
      ii = nblock*n;
      for ( integer i = 0; i < n+q; ++i )
        for ( integer j = 0; j < n; ++j )
          stream << ii+i << '\t' << j << '\t' << H0[i+j*(n+q)] << '\n';

      valueType * HNq = H0Nq+n*(n+q);
      for ( integer i = 0; i < n+q; ++i )
        for ( integer j = 0; j < n+q; ++j )
          stream << ii+i << '\t' << ii+j << '\t' << HNq[i+j*(n+q)] << '\n';
    }

    // bidiagonal
    for ( integer k = 0; k < nblock; ++k ) {
      ii = k*n;
      valueType * DE = DE_blk + k * nxnx2;
      for ( integer i = 0; i < n; ++i )
        for ( integer j = 0; j < nx2; ++j )
          stream << ii+i << '\t' << ii+j << '\t' << DE[i+j*n] << '\n';
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < ii; ++j )
        stream
          << ii+i << '\t' << j << '\t' << Cmat[i+j*nb] << '\n'
          << j << '\t' << ii+i << '\t' << Bmat[j+i*ii] << '\n';

    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < nb; ++j )
        stream << ii+i << '\t' << ii+j << '\t' << Dmat[i+j*nb] << '\n';

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
    integer nnz = nblock*nxnx2 + 2*(nblock+1)*n*nb + nb*nb;
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {
      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;
      nnz += row0 * ( n + col00 ) + rowN * ( n + colNN );
    } else {
      nnz += (n+q)*(2*n+q);
    }
    return nnz;
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparsePattern( integer I[], integer J[] ) const {
    integer kkk = 0;

    // BC
    integer ii;
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {

      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;

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

      ii = nblock*n;
      for ( integer i = 0; i < n+q; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          I[kkk] = ii+i;
          J[kkk] = j;
          ++kkk;
        }
      }

      for ( integer i = 0; i < n+q; ++i ) {
        for ( integer j = 0; j < n+q; ++j ) {
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

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BlockBidiagonal::sparsePattern( V ), inserted " << kkk <<
      " values, expected " << sparseNnz()
    );
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::sparseValues( valueType V[] ) const {
    integer kkk = 0;

    // BC
    integer ii;
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {

      integer row0  = numInitialBC;
      integer rowN  = numFinalBC;
      integer col00 = numInitialOMEGA;
      integer colNN = numFinalOMEGA;

      ii = nblock*n;
      for ( integer i = 0; i < rowN; ++i )
        for ( integer j = 0; j < n+colNN; ++j )
          V[kkk++] = blockN[i+j*rowN];

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < n; ++j )
          V[kkk++] = block0[i+(j+col00)*row0];

      for ( integer i = 0; i < row0; ++i )
        for ( integer j = 0; j < col00; ++j )
          V[kkk++] = block0[i+j*row0];

    } else {

      valueType * H0 = H0Nq;
      ii = nblock*n;
      for ( integer i = 0; i < n+q; ++i )
        for ( integer j = 0; j < n; ++j )
          V[kkk++] = H0[i+j*(n+q)];

      valueType * HNq = H0Nq+n*(n+q);
      for ( integer i = 0; i < n+q; ++i )
        for ( integer j = 0; j < n+q; ++j )
          V[kkk++] = HNq[i+j*(n+q)];
    }

    // bidiagonal
    for ( integer k = 0; k < nblock; ++k ) {
      ii = k*n;
      valueType * DE = DE_blk + k * nxnx2;
      for ( integer i = 0; i < n; ++i )
        for ( integer j = 0; j < nx2; ++j )
          V[kkk++] = DE[i+j*n];
    }

    // border
    ii = n*(nblock+1)+q;
    for ( integer i = 0; i < nb; ++i ) {
      for ( integer j = 0; j < ii; ++j ) {
        V[kkk++] = Cmat[i+j*nb];
        V[kkk++] = Bmat[j+i*ii];
      }
    }

    for ( integer i = 0; i < nb; ++i )
      for ( integer j = 0; j < nb; ++j )
        V[kkk++] = Dmat[i+j*nb];

    ALGLIN_ASSERT(
      kkk == sparseNnz(),
      "BlockBidiagonal::sparseValues( V ), inserted " << kkk <<
      " values, expected " << sparseNnz()
    );

  }

  template class BlockBidiagonal<float>;
  template class BlockBidiagonal<double>;

}
