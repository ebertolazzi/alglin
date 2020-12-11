/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
/// file: Alglin_aux.hxx
///

namespace alglin {

  /*
      Matrix NNZ structure
          col0
       |       |
     / +-------+                         \
     | |  TOP  |                         | <-- row0
     | +-------+----+                    |
     |    |    |    |                    | <- sizeBlock
     |    +----+----+----+               |
     |         |    |    |               |
     |         +----+----+----+          |
     |              |    |    |          |
     |              +----+----+----+     |
     |                   |    |    |     |
     |                   +----+----+---+ |
     |                        |        | | <-- rowN
     |                        | BOTTOM | |
     \                        +--------+ /
                              |        |
                                 colN
  */

  //! compute y = alpha*A*x+beta*y
  template <typename t_Value>
  inline
  void
  abd_mv(
    integer         row0,
    integer         col0,
    t_Value const * block0,
    integer         numBlock,
    integer         dimBlock,
    t_Value const * blocks,
    integer         rowN,
    integer         colN,
    t_Value const * blockN,
    t_Value         alpha,
    t_Value const * x,
    integer         incx,
    t_Value         beta,
    t_Value *       y,
    integer         incy
  ) {

    // first block y = alpha * _block0 * x + beta * y
    gemv(
      NO_TRANSPOSE, row0, col0,
      alpha, block0, row0,
      x, incx,
      beta, y, incy
    );

    // internal blocks block
    t_Value const * xx   = x+(col0-dimBlock)*incx;
    t_Value *       yy   = y+row0*incy;
    t_Value const * blks = blocks;
    for ( integer i = 0; i < numBlock; ++i ) {
      gemv(
        NO_TRANSPOSE, dimBlock, 2*dimBlock,
        alpha, blks, dimBlock,
        xx, incx,
        beta, yy, incy
      );
      xx   += dimBlock*incx;
      yy   += dimBlock*incy;
      blks += 2*dimBlock*dimBlock;
    }

    // last block
    gemv(
      NO_TRANSPOSE, rowN, colN,
      alpha, blockN, rowN,
      xx, incx,
      beta, yy, incy
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //! compute r = b-A*x
  template <typename t_Value>
  inline
  void
  abd_residue(
    integer         row0,
    integer         col0,
    t_Value const * block0,
    integer         numBlock,
    integer         dimBlock,
    t_Value const * blocks,
    integer         rowN,
    integer         colN,
    t_Value const * blockN,
    t_Value const * b,
    integer         incb,
    t_Value const * x,
    integer         incx,
    t_Value *       res,
    integer         incr
  ) {
    copy( numBlock*dimBlock+row0+rowN, b, incb, res, incr );
    abd_mv(
      row0, col0, block0,
      numBlock, dimBlock, blocks,
      rowN, colN, blockN,
      t_Value(-1.0), x, incx, t_Value(1.0), res, incr
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  abd_print(
    ostream_type &  stream,
    integer         row0,
    integer         col0,
    t_Value const * block0,
    integer         numBlock,
    integer         dimBlock,
    t_Value const * blocks,
    integer         rowN,
    integer         colN,
    t_Value const * blockN
  ) {
    integer sizeBlock = 2*dimBlock*dimBlock;
    stream << "Block 0\n";
    for ( integer i = 0; i < row0; ++i ) {
      stream << std::setw(8) << block0[i];
      for ( integer j = 1; j < col0; ++j )
        stream << ' ' << std::setw(8) << block0[i+j*row0];
      stream << '\n';
    }
    for ( integer k = 0; k < numBlock; ++k ) {
      stream << "Block " << k+1 << '\n';
      t_Value const * blk = blocks+k*sizeBlock;
      for ( integer i = 0; i < dimBlock; ++i ) {
        stream << std::setw(8) << blk[i];
        for ( integer j = 1; j < 2*dimBlock; ++j )
          stream << ' ' << std::setw(8) << blk[i+j*dimBlock];
        stream << '\n';
      }
    }
    stream << "Block N\n";
    for ( integer i = 0; i < rowN; ++i ) {
      stream << std::setw(8) << blockN[i];
      for ( integer j = 1; j < colN; ++j )
        stream << ' ' << std::setw(8) << blockN[i+j*rowN];
      stream << '\n';
    }
  }

  /*!
    Matrix structure

                      n * nblock
        ________________^_________________
       /                                  \
         n     n     n                        n     q
      +-----+-----+-----+----.........-----+-----+-----+    \
      |  Ad | Au  |  0  |                  |  0  |  0  | n   |
      +-----+-----+-----+             -----+-----+-----+     |
      |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |
      |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |
      |                                                :     |
      :                                                :      > n * nblock
      :                                                :     |
      :                                                :     |
      :                                                :     |
      :                              +-----+-----+-----+     |
      :                              | Au  |  0  |  0  |     |
      :                        +-----+-----+-----+-----+     |
      :                        |  0  | Ad  | Au  |  0  | n   |
      +-----+-----+---......---+-----+-----+=====+=====+    /
      |     |     |            |     |     !     |     !
      | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
      |     |     |            |     |     !     |     !
      +-----+-----+---......---+-----+-----+=====+=====+
  */
  //! compute y = alpha*A*x+beta*y
  template <typename t_Value>
  inline
  void
  babd_mv(
    integer         nblk,
    integer         n,
    integer         q,
    t_Value const * AdAu,
    t_Value const * H0,
    t_Value const * HN,
    t_Value const * Hq,
    t_Value         alpha,
    t_Value const * x,
    integer         incx,
    t_Value         beta,
    t_Value *       y,
    integer         incy
  ) {

    // internal blocks block
    t_Value const * xx   = x;
    t_Value *       yy   = y;
    t_Value const * blks = AdAu;
    for ( integer i = 0; i < nblk; ++i ) {
      gemv(
        NO_TRANSPOSE, n, 2*n,
        alpha, blks, n,
        xx, incx,
        beta, yy, incy
      );
      xx   += n*incx;
      yy   += n*incy;
      blks += 2*n*n;
    }

    // last blocks
    integer nq = n+q;
    gemv( NO_TRANSPOSE, nq, n, alpha, H0, nq, x, incx, beta, yy, incy );
    gemv( NO_TRANSPOSE, nq, n, alpha, HN, nq, xx, incx, t_Value(1), yy, incy );

    xx += n*incx;
    gemv( NO_TRANSPOSE, nq, q, alpha, Hq, nq, xx, incx, t_Value(1), yy, incy );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  babd_residue(
    integer         nblk,
    integer         n,
    integer         q,
    t_Value const * AdAu,
    t_Value const * H0,
    t_Value const * HN,
    t_Value const * Hq,
    t_Value const * b,
    integer         incb,
    t_Value const * x,
    integer         incx,
    t_Value *       res,
    integer         incr
  ) {
    copy( nblk*n+n+q, b, incb, res, incr );
    babd_mv(
      nblk, n, q, AdAu, H0, HN, Hq,
      t_Value(-1.0), x, incx, t_Value(1.0), res, incr
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  babd_print(
    ostream_type  & stream,
    integer         nblk,
    integer         n,
    integer         q,
    t_Value const * AdAu,
    t_Value const * H0,
    t_Value const * HN,
    t_Value const * Hq
  ) {
    integer sizeBlock = 2*n*n;
    for ( integer k = 0; k < nblk; ++k ) {
      stream << "Block " << k+1 << '\n';
      t_Value const * blk = AdAu+k*sizeBlock;
      for ( integer i = 0; i < n; ++i ) {
        stream << std::setw(8) << blk[i];
        for ( integer j = 1; j < 2*n; ++j )
          stream << ' ' << std::setw(8) << blk[i+j*n];
        stream << '\n';
      }
    }
    integer nq = n+q;
    stream << "Block H0\n";
    for ( integer i = 0; i < nq; ++i ) {
      stream << std::setw(8) << H0[i];
      for ( integer j = 1; j < n; ++j )
        stream << ' ' << std::setw(8) << H0[i+j*nq];
      stream << '\n';
    }
    stream << "Block HN\n";
    for ( integer i = 0; i < nq; ++i ) {
      stream << std::setw(8) << HN[i];
      for ( integer j = 1; j < n; ++j )
        stream << ' ' << std::setw(8) << HN[i+j*nq];
      stream << '\n';
    }
    if ( q > 0 ) {
      stream << "Block Hq\n";
      for ( integer i = 0; i < nq; ++i ) {
        stream << std::setw(8) << Hq[i];
        for ( integer j = 1; j < q; ++j )
          stream << ' ' << std::setw(8) << Hq[i+j*nq];
        stream << '\n';
      }
    }
  }

  /*
  //              _                 __
  //    __ _  ___| |_ _ ____  __   / /   _
  //   / _` |/ _ \ __| '__\ \/ /  / / | | |
  //  | (_| |  __/ |_| |   >  <  / /| |_| |
  //   \__, |\___|\__|_|  /_/\_\/_/  \__, |
  //   |___/                         |___/
  */
  template <typename REAL>
  integer
  getrx(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[],
    integer NB
  );

  template <typename REAL>
  integer
  getry(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[],
    integer NB
  );

  template <typename REAL>
  integer
  gtx(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[]
  );

  template <typename REAL>
  integer
  gty(
    integer M,
    integer N,
    REAL    A[],
    integer LDA,
    integer IPIV[]
  );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  integer
  equilibrate(
    integer M,
    integer N,
    T const A[],
    integer LDA,
    T       R[],
    T       C[],
    integer maxIter,
    T       epsi
  );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  void
  triTikhonov(
    integer N,
    T const Tmat[],
    integer LDT,
    integer nrhs,
    T       RHS[],
    integer ldRHS,
    T       lambda
  );

  inline
  bool
  outMATRIXcheck( MatrixType const & MT, integer i, integer j ) {
    bool ok = MT == FULL_MATRIX ||
              ( MT == LOWER_TRIANGULAR_MATRIX && i >= j ) ||
              ( MT == UPPER_TRIANGULAR_MATRIX && i <= j );
    return ok;
  }

  template <typename T>
  inline
  void
  outMATRIX(
    MatrixType const & MT,
    integer            NR,
    integer            NC,
    T const            A[],
    integer            LDA,
    ostream_type     & s,
    integer            prec = 4,
    integer            rperm[] = nullptr,
    integer            cperm[] = nullptr
  ) {
    integer j0 = cperm == nullptr ? 0 : cperm[0]-1;
    for ( integer i = 0; i < NR; ++i ) {
      integer ii = rperm == nullptr ? i : rperm[i]-1;
      if ( outMATRIXcheck(MT,i,0) )
        s << std::setprecision(prec) << std::setw(prec+6) << A[ii+j0*LDA];
      else
        s << std::setw(prec+6) << " ";
      for ( integer j = 1; j < NC; ++j ) {
        integer jj = cperm == nullptr ? j : cperm[j]-1;
        if ( outMATRIXcheck(MT,i,j) )
          s << " " << std::setprecision(prec) << std::setw(prec+6)
            << A[ii+jj*LDA];
        else
          s << " " << std::setw(prec+6) << " ";
      }
      s << '\n';
    }
  }

  template <typename T>
  inline
  void
  outMAPLE(
    integer        NR,
    integer        NC,
    T const        A[],
    integer        LDA,
    ostream_type & s
  ) {
    s << "<";
    for ( integer j = 0; j < NC; ++j ) {
      s << "<" << std::setprecision(20) << A[j*LDA];
      for ( integer i = 1; i < NR; ++i )
        s << "," << std::setprecision(20) << A[i+j*LDA];
      if ( j < NC-1 ) s << ">|\n";
      else            s << ">>;\n";
    }
  }

} // end namespace alglin

///
/// eof: Alglin_aux.hxx
///
