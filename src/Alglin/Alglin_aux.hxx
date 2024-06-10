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
     |    |    |    |                    | <- size_block
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

  //!
  //! compute y = alpha*A*x+beta*y
  //!
  template <typename t_Value>
  inline
  void
  abd_mv(
    integer       row0,
    integer       col0,
    t_Value const block0[],
    integer       num_block,
    integer       dim_block,
    t_Value const blocks[],
    integer       rowN,
    integer       colN,
    t_Value const blockN[],
    t_Value       alpha,
    t_Value const x[],
    integer       incx,
    t_Value       beta,
    t_Value       y[],
    integer       incy
  ) {

    // first block y = alpha * _block0 * x + beta * y
    gemv(
      Transposition::NO,
      row0, col0,
      alpha, block0, row0,
      x, incx,
      beta, y, incy
    );

    // internal blocks block
    t_Value const * xx{ x+(col0-dim_block)*incx};
    t_Value *       yy{ y+row0*incy };
    t_Value const * blks{ blocks };
    for ( integer i{0}; i < num_block; ++i ) {
      gemv(
        Transposition::NO,
        dim_block, 2*dim_block,
        alpha, blks, dim_block,
        xx, incx,
        beta, yy, incy
      );
      xx   += dim_block*incx;
      yy   += dim_block*incy;
      blks += 2*dim_block*dim_block;
    }

    // last block
    gemv(
      Transposition::NO, rowN, colN,
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
    integer       row0,
    integer       col0,
    t_Value const block0[],
    integer       num_block,
    integer       dim_block,
    t_Value const blocks[],
    integer       rowN,
    integer       colN,
    t_Value const blockN[],
    t_Value const b[],
    integer       incb,
    t_Value const x[],
    integer       incx,
    t_Value       res[],
    integer       incr
  ) {
    alglin::copy( num_block*dim_block+row0+rowN, b, incb, res, incr );
    abd_mv(
      row0, col0, block0,
      num_block, dim_block, blocks,
      rowN, colN, blockN,
      t_Value(-1.0), x, incx, t_Value(1.0), res, incr
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  abd_print(
    ostream_type & stream,
    integer        row0,
    integer        col0,
    t_Value const  block0[],
    integer        num_block,
    integer        dim_block,
    t_Value const  blocks[],
    integer        rowN,
    integer        colN,
    t_Value const  blockN[]
  ) {
    integer size_block{ 2*dim_block*dim_block };
    fmt::print( stream, "Block 0\n");
    for ( integer i{0}; i < row0; ++i ) {
      fmt::print( stream, "{:8}", block0[i] );
      for ( integer j{1}; j < col0; ++j )
        fmt::print( stream, " {:8}", block0[i+j*row0] );
      fmt::print( stream, "\n" );
    }
    for ( integer k{0}; k < num_block; ++k ) {
      fmt::print( stream, "Block {}\n", k+1 );
      t_Value const * blk{blocks+k*size_block};
      for ( integer i{0}; i < dim_block; ++i ) {
        fmt::print( stream, "{:8}", blk[i] );
        for ( integer j{1}; j < 2*dim_block; ++j )
          fmt::print( stream, " {:8}", blk[i+j*dim_block] );
        fmt::print( stream, "\n" );
      }
    }
    fmt::print( stream, "Block N\n" );
    for ( integer i{0}; i < rowN; ++i ) {
      fmt::print( stream, "{:8}", blockN[i] );
      for ( integer j{1}; j < colN; ++j )
        fmt::print( stream, " {:8}", blockN[i+j*rowN] );
      fmt::print( stream, "\n" );
    }
  }

  //!
  //! Matrix structure
  //!
  //! \verbatim
  //!                   n * nblock
  //!     ________________^_________________
  //!    /                                  \
  //!      n     n     n                        n     q
  //!   +-----+-----+-----+----.........-----+-----+-----+    \
  //!   |  Ad | Au  |  0  |                  |  0  |  0  | n   |
  //!   +-----+-----+-----+             -----+-----+-----+     |
  //!   |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
  //!   +-----+-----+-----+-----+       -----+-----+-----+     |
  //!   |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
  //!   +-----+-----+-----+-----+       -----+-----+-----+     |
  //!   |                                                :     |
  //!   :                                                :      > n * nblock
  //!   :                                                :     |
  //!   :                                                :     |
  //!   :                                                :     |
  //!   :                              +-----+-----+-----+     |
  //!   :                              | Au  |  0  |  0  |     |
  //!   :                        +-----+-----+-----+-----+     |
  //!   :                        |  0  | Ad  | Au  |  0  | n   |
  //!   +-----+-----+---......---+-----+-----+=====+=====+    /
  //!   |     |     |            |     |     !     |     !
  //!   | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
  //!   |     |     |            |     |     !     |     !
  //!   +-----+-----+---......---+-----+-----+=====+=====+
  //!
  //! \endverbatim
  //! compute y = alpha*A*x+beta*y
  //!
  template <typename t_Value>
  inline
  void
  babd_mv(
    integer       nblk,
    integer       n,
    integer       q,
    t_Value const AdAu[],
    t_Value const H0[],
    t_Value const HN[],
    t_Value const Hq[],
    t_Value       alpha,
    t_Value const x[],
    integer       incx,
    t_Value       beta,
    t_Value       y[],
    integer       incy
  ) {

    // internal blocks block
    t_Value const * xx{x};
    t_Value *       yy{y};
    t_Value const * blks{AdAu};
    for ( integer i{0}; i < nblk; ++i ) {
      gemv(
        Transposition::NO, n, 2*n,
        alpha, blks, n,
        xx, incx,
        beta, yy, incy
      );
      xx   += n*incx;
      yy   += n*incy;
      blks += 2*n*n;
    }

    // last blocks
    integer nq{n+q};
    gemv(Transposition:: NO, nq, n, alpha, H0, nq, x, incx, beta, yy, incy );
    gemv( Transposition::NO, nq, n, alpha, HN, nq, xx, incx, t_Value(1), yy, incy );

    xx += n*incx;
    gemv( Transposition::NO, nq, q, alpha, Hq, nq, xx, incx, t_Value(1), yy, incy );

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  inline
  void
  babd_residue(
    integer       nblk,
    integer       n,
    integer       q,
    t_Value const AdAu[],
    t_Value const H0[],
    t_Value const HN[],
    t_Value const Hq[],
    t_Value const b[],
    integer       incb,
    t_Value const x[],
    integer       incx,
    t_Value       res[],
    integer       incr
  ) {
    alglin::copy( nblk*n+n+q, b, incb, res, incr );
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
    ostream_type & stream,
    integer        nblk,
    integer        n,
    integer        q,
    t_Value const  AdAu[],
    t_Value const  H0[],
    t_Value const  HN[],
    t_Value const  Hq[]
  ) {
    integer size_block{ 2*n*n };
    for ( integer k{0}; k < nblk; ++k ) {
      fmt::print( stream, "Block {}\n", k+1 );
      t_Value const * blk{ AdAu+k*size_block };
      for ( integer i{0}; i < n; ++i ) {
        fmt::print( stream, "{:8}", blk[i] );
        for ( integer j{1}; j < 2*n; ++j )
         fmt::print( stream, " {:8}", blk[i+j*n] );
        fmt::print( stream, "\n" );
      }
    }
    integer nq{n+q};
    fmt::print( stream, "Block H0\n" );
    for ( integer i{0}; i < nq; ++i ) {
      fmt::print( stream, "{:8}", H0[i] );
      for ( integer j{1}; j < n; ++j )
        fmt::print( stream, " {:8}", H0[i+j*nq] );
      fmt::print( stream, "\n" );
    }
    fmt::print( stream, "Block HN\n" );
    for ( integer i{0}; i < nq; ++i ) {
      fmt::print( stream, "{:8}", HN[i] );
      for ( integer j{1}; j < n; ++j )
        fmt::print( stream, " {:8}", HN[i+j*nq] );
      fmt::print( stream, "\n" );
    }
    if ( q > 0 ) {
      fmt::print( stream, "Block Hq\n" );
      for ( integer i{0}; i < nq; ++i ) {
        fmt::print( stream, "{:8}", Hq[i] );
        for ( integer j{1}; j < q; ++j )
          fmt::print( stream, " {:8}", Hq[i+j*nq] );
        fmt::print( stream, "\n" );
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

  extern template integer gtx( integer M, integer N, float  A[], integer LDA, integer IPIV[] );
  extern template integer gtx( integer M, integer N, double A[], integer LDA, integer IPIV[] );
  extern template integer gty( integer M, integer N, float  A[], integer LDA, integer IPIV[] );
  extern template integer gty( integer M, integer N, double A[], integer LDA, integer IPIV[] );

  extern template integer getrx( integer M, integer N, float  A[], integer LDA, integer IPIV[], integer MB );
  extern template integer getrx( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );
  extern template integer getry( integer M, integer N, float  A[], integer LDA, integer IPIV[], integer MB  );
  extern template integer getry( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  );

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
    bool ok = MT == MatrixType::FULL ||
              ( MT == MatrixType::LOWER_TRIANGULAR && i >= j ) ||
              ( MT == MatrixType::UPPER_TRIANGULAR && i <= j );
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
    using std::setprecision;
    using std::setw;
    integer j0{ cperm == nullptr ? 0 : cperm[0]-1 };
    for ( integer i{0}; i < NR; ++i ) {
      integer ii{ rperm == nullptr ? i : rperm[i]-1 };
      if ( outMATRIXcheck(MT,i,0) )
        s << setprecision(prec) << setw(prec+6) << A[ii+j0*LDA];
      else
        s << setw(prec+6) << " ";
      for ( integer j{1}; j < NC; ++j ) {
        integer jj{cperm == nullptr ? j : cperm[j]-1};
        if ( outMATRIXcheck(MT,i,j) )
          s << " " << setprecision(prec) << setw(prec+6)
            << A[ii+jj*LDA];
        else
          s << " " << setw(prec+6) << " ";
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
    using std::setprecision;
    using std::setw;
    fmt::print( s, "<" );
    for ( integer j{0}; j < NC; ++j ) {
      fmt::print( s, "<{:20}", A[j*LDA] );
      for ( integer i{1}; i < NR; ++i )
        fmt::print( s, ",{:20}", A[i+j*LDA] );
      if ( j < NC-1 ) fmt::print( s, ">|\n" );
      else            fmt::print( s, ">>;\n" );
    }
  }

} // end namespace alglin

///
/// eof: Alglin_aux.hxx
///
