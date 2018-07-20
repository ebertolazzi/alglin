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
/// file: Alglin_aux.hh
///

#ifndef ALGLIN_AUX_HH
#define ALGLIN_AUX_HH

/*
*/

#include "Alglin.hh"

namespace alglin {

  /*\
   *  Purpose
   *  =======
   *
   *  DLARNV returns a vector of n random real numbers from a uniform or
   *  normal distribution.
   *
   *  Arguments
   *  =========
   *
   *  IDIST   (input) INTEGER
   *          Specifies the distribution of the random numbers:
   *          = 1:  uniform (0,1)
   *          = 2:  uniform (-1,1)
   *          = 3:  normal (0,1)
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry, the seed of the random number generator; the array
   *          elements must be between 0 and 4095, and ISEED(4) must be
   *          odd.
   *          On exit, the seed is updated.
   *
   *  N       (input) INTEGER
   *          The number of random numbers to be generated.
   *
   *  X       (output) DOUBLE PRECISION array, dimension (N)
   *          The generated random numbers.
  \*/

  // use standard Lapack routine
  #if defined(ALGLIN_USE_LAPACK) || defined(ALGLIN_USE_OPENBLAS) || defined(ALGLIN_USE_ATLAS)
  extern "C" {

    integer
    LAPACK_F77NAME(slarnv)( integer * IDIST,
                            integer * ISEED,
                            integer * N,
                            real      X[] );

    integer
    LAPACK_F77NAME(dlarnv)( integer  * IDIST,
                            integer  * ISEED,
                            integer  * N,
                            doublereal X[] );
  }
  #endif

  inline
  integer
  larnv( integer IDIST,
         integer ISEED[4],
         integer N,
         real    X[] ) {
    #if defined(ALGLIN_USE_LAPACK) || defined(ALGLIN_USE_OPENBLAS) || defined(ALGLIN_USE_ATLAS)
    LAPACK_F77NAME(slarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(ALGLIN_USE_MKL)
    slarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(ALGLIN_USE_ACCELERATE)
    return CLAPACKNAME(slarnv)( &IDIST, ISEED, &N, X );
    #else
    return LAPACKNAME(slarnv)( &IDIST, ISEED, &N, X );
    #endif
  }

  inline
  integer
  lanrv( integer    IDIST,
         integer    ISEED[4],
         integer    N,
         doublereal X[] ) {
    #if defined(ALGLIN_USE_LAPACK) || defined(ALGLIN_USE_OPENBLAS) || defined(ALGLIN_USE_ATLAS)
    LAPACK_F77NAME(dlarnv)( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(ALGLIN_USE_MKL)
    dlarnv( &IDIST, ISEED, &N, X ); // return void in openblas
    return 0;
    #elif defined(ALGLIN_USE_ACCELERATE)
    return CLAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X );
    #else
    return LAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X );
    #endif
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DLARGE pre- and post-multiplies a real general n by n matrix A
   *  with a random orthogonal matrix: A = U*D*U'.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the original n by n matrix A.
   *          On exit, A is overwritten by U*A*U' for some random
   *          orthogonal matrix U.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= N.
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry, the seed of the random number generator; the array
   *          elements must be between 0 and 4095, and ISEED(4) must be
   *          odd.
   *          On exit, the seed is updated.
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/
  
  template <typename T>
  inline
  integer
  large( integer N,
         T       A[],
         integer LDA,
         integer ISEED[4],
         T       WORK[],
         integer LWORK ) {
    ALGLIN_ASSERT( LWORK >= 2*N,
                   "large, LWORK = " << LWORK << " must be >= " << 2*N );
    // Test the input arguments
    if ( N < 0 ) return -1;
    if ( LDA < max_index( 1, N ) ) return -3;
    // pre- and post-multiply A by random orthogonal matrix
    for ( integer i = N-1; i >= 0; --i ) {
      // generate random reflection
      integer Ni = N-i;
      lanrv( 3, ISEED, N, WORK );
      T WN = nrm2( Ni, WORK, 1 );
      T WA = WN;
      if ( WORK[0] < 0 ) WA = -WN;
      T TAU = 0;
      if ( WN > 0 ) {
        T WB = WORK[0] + WA;
        scal( Ni-1, 1.0/WB, WORK+1, 1 );
        WORK[0] = 1;
        TAU = WB / WA;
      }
      // multiply A(i:n,1:n) by random reflection from the left
      gemv( TRANSPOSE, Ni, N, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+N, 1 );
      ger( Ni, N, -TAU, WORK, 1, WORK+N, 1, A+i, LDA );
      // multiply A(1:n,i:n) by random reflection from the right
      gemv( NO_TRANSPOSE, N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 );
      ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA );
    }
    return 0;
  }
  
  /*\
   *  Purpose
   *  =======
   *
   *  SLAROR pre- or post-multiplies an M by N matrix A by a random
   *  orthogonal matrix U, overwriting A.  A may optionally be initialized
   *  to the identity matrix before multiplying by U.  U is generated using
   *  the method of G.W. Stewart (SIAM J. Numer. Anal. 17, 1980, 403-409).
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          Specifies whether A is multiplied on the left or right by U.
   *          = 'L':         Multiply A on the left (premultiply) by U
   *          = 'R':         Multiply A on the right (postmultiply) by U'
   *          = 'C' or 'T':  Multiply A on the left by U and the right
   *                          by U' (Here, U' means U-transpose.)
   *
   *  INIT    (input) CHARACTER*1
   *          Specifies whether or not A should be initialized to the
   *          identity matrix.
   *          = 'I':  Initialize A to (a section of) the identity matrix
   *                   before applying U.
   *          = 'N':  No initialization.  Apply U to the input matrix A.
   *
   *          INIT = 'I' may be used to generate square or rectangular
   *          orthogonal matrices:
   *
   *          For M = N and SIDE = 'L' or 'R', the rows will be orthogonal
   *          to each other, as will the columns.
   *
   *          If M < N, SIDE = 'R' produces a dense matrix whose rows are
   *          orthogonal and whose columns are not, while SIDE = 'L'
   *          produces a matrix whose rows are orthogonal, and whose first
   *          M columns are orthogonal, and whose remaining columns are
   *          zero.
   *
   *          If M > N, SIDE = 'L' produces a dense matrix whose columns
   *          are orthogonal and whose rows are not, while SIDE = 'R'
   *          produces a matrix whose columns are orthogonal, and whose
   *          first M rows are orthogonal, and whose remaining rows are
   *          zero.
   *
   *  M       (input) INTEGER
   *          The number of rows of A.
   *
   *  N       (input) INTEGER
   *          The number of columns of A.
   *
   *  A       (input/output) REAL array, dimension (LDA, N)
   *          On entry, the array A.
   *          On exit, overwritten by U A ( if SIDE = 'L' ),
   *           or by A U ( if SIDE = 'R' ),
   *           or by U A U' ( if SIDE = 'C' or 'T').
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  ISEED   (input/output) INTEGER array, dimension (4)
   *          On entry ISEED specifies the seed of the random number
   *          generator. The array elements should be between 0 and 4095;
   *          if not they will be reduced mod 4096.  Also, ISEED(4) must
   *          be odd.  The random number generator uses a linear
   *          congruential sequence limited to small integers, and so
   *          should produce machine independent random numbers. The
   *          values of ISEED are changed on exit, and can be used in the
   *          next call to SLAROR to continue the same random number
   *          sequence.
   *
   *  X       (workspace) REAL array, dimension (3*MAX( M, N ))
   *          Workspace of length
   *              2*M + N if SIDE = 'L',
   *              2*N + M if SIDE = 'R',
   *              3*N     if SIDE = 'C' or 'T'.
   *
   *  INFO    (output) INTEGER
   *          An error flag.  It is set to:
   *          = 0:  normal return
   *          < 0:  if INFO = -k, the k-th argument had an illegal value
   *          = 1:  if the random numbers generated by SLARND are bad.
   *
   *  =====================================================================
  \*/
  
  template <typename T>
  inline
  integer
  laror( SideMultiply const & LR,
         bool                 init,
         integer              M,
         integer              N,
         T                    A[],
         integer              LDA,
         integer              ISEED[4],
         T                    WORK[],
         integer              LWORK ) {

    if ( LR == LEFT ) {
      ALGLIN_ASSERT( LWORK >= 2*M + N,
                     "laror, LWORK = " << LWORK << " must be >= " << 2*M + N );
    } else if ( LR == RIGHT ) {
      ALGLIN_ASSERT( LWORK >= 2*N + M,
                     "laror, LWORK = " << LWORK << " must be >= " << 2*N + M );
    } else {
      return -1;
    }

    if ( N == 0 || M == 0 ) return 0;

    // Check for argument errors.
    if ( M   < 0 ) return -3;
    if ( N   < 0 ) return -4;
    if ( LDA < M ) return -6;

    // Initialize A to the identity matrix if desired
    if ( init ) geid( M, N, A, LDA );

    ALGLIN_ASSERT( LWORK >= 2*N,
                   "large, LWORK = " << LWORK << " must be >= " << 2*N );
    // Test the input arguments
    if ( N < 0 ) return -1;
    if ( LDA < max_index( 1, N ) ) return -3;
    
    if ( LR == LEFT ) {
      // pre-multiply A by random orthogonal matrix
      for ( integer i = M-1; i >= 0; --i ) {
        // generate random reflection
        integer Mi = M-i;
        lanrv( 1, ISEED, M, WORK );
        T WN = nrm2( Mi, WORK, 1 );
        T WA = WN;
        if ( WORK[0] < 0 ) WA = -WN;
        T TAU = 0;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA;
          scal( Mi-1, 1.0/WB, WORK+1, 1 );
          WORK[0] = 1;
          TAU = WB / WA;
        }
        // multiply A(i:n,1:n) by random reflection from the left
        gemv( TRANSPOSE, Mi, M, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+M, 1 );
        ger( Mi, N, -TAU, WORK, 1, WORK+M, 1, A+i, LDA );
      }
    } else {
      // post-multiply A by random orthogonal matrix
      for ( integer i = N-1; i >= 0; --i ) {
        // generate random reflection
        integer Ni = N-i;
        lanrv( 1, ISEED, N, WORK );
        T WN = nrm2( Ni, WORK, 1 );
        T WA = WN;
        if ( WORK[0] < 0 ) WA = -WN;
        T TAU = 0;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA;
          scal( Ni-1, 1.0/WB, WORK+1, 1 );
          WORK[0] = 1;
          TAU = WB / WA;
        }
        // multiply A(1:n,i:n) by random reflection from the right
        gemv( NO_TRANSPOSE, N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 );
        ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA );
      }
    }

    return 0;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  inline
  void
  print_matrix( std::basic_ostream<char> & stream,
                integer                    nr,
                integer                    nc,
                t_Value const              A[],
                integer                    ldA ) {
    for ( integer i = 0; i < nr; ++i ) {
      for ( integer j = 0; j < nc; ++j )
        stream << std::setw(14) << A[i+j*ldA] << " ";
      stream << '\n';
    }
  }

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
  abd_mv( integer         row0,
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
          integer         incy ) {

    // first block y = alpha * _block0 * x + beta * y
    gemv( NO_TRANSPOSE, row0, col0,
          alpha, block0, row0,
          x, incx,
          beta, y, incy );

    // internal blocks block
    t_Value const * xx   = x+(col0-dimBlock)*incx;
    t_Value *       yy   = y+row0*incy;
    t_Value const * blks = blocks;
    for ( integer i = 0; i < numBlock; ++i ) {
      gemv( NO_TRANSPOSE, dimBlock, 2*dimBlock,
            alpha, blks, dimBlock,
            xx, incx,
            beta, yy, incy );
      xx   += dimBlock*incx;
      yy   += dimBlock*incy;
      blks += 2*dimBlock*dimBlock;
    }

    // last block
    gemv( NO_TRANSPOSE, rowN, colN,
          alpha, blockN, rowN,
          xx, incx,
          beta, yy, incy );
  }

  //! compute r = b-A*x
  template <typename t_Value>
  inline
  void
  abd_residue( integer         row0,
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
               integer         incr ) {
    copy( numBlock*dimBlock+row0+rowN, b, incb, res, incr );
    abd_mv( row0, col0, block0,
            numBlock, dimBlock, blocks,
            rowN, colN, blockN,
            t_Value(-1.0), x, incx, t_Value(1.0), res, incr );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  inline
  void
  abd_print( std::basic_ostream<char> & stream,
             integer                    row0,
             integer                    col0,
             t_Value const *            block0,
             integer                    numBlock,
             integer                    dimBlock,
             t_Value const *            blocks,
             integer                    rowN,
             integer                    colN,
             t_Value const *            blockN ) {
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
  babd_mv( integer         nblk,
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
           integer         incy ) {

    // internal blocks block
    t_Value const * xx   = x;
    t_Value *       yy   = y;
    t_Value const * blks = AdAu;
    for ( integer i = 0; i < nblk; ++i ) {
      gemv( NO_TRANSPOSE, n, 2*n,
            alpha, blks, n,
            xx, incx,
            beta, yy, incy );
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

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  inline
  void
  babd_residue( integer         nblk,
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
                integer         incr ) {
    copy( nblk*n+n+q, b, incb, res, incr );
    babd_mv( nblk, n, q, AdAu, H0, HN, Hq,
             t_Value(-1.0), x, incx, t_Value(1.0), res, incr );
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  inline
  void
  babd_print( std::basic_ostream<char> & stream,
              integer                    nblk,
              integer                    n,
              integer                    q,
              t_Value const *            AdAu,
              t_Value const *            H0,
              t_Value const *            HN,
              t_Value const *            Hq ) {
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
  getrx( integer M,
         integer N,
         REAL    A[],
         integer LDA,
         integer IPIV[],
         integer NB );

  template <typename REAL>
  integer
  getry( integer M,
         integer N,
         REAL    A[],
         integer LDA,
         integer IPIV[],
         integer NB );

  template <typename REAL>
  integer
  gtx( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] );

  template <typename REAL>
  integer
  gty( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  integer
  equilibrate( integer M,
               integer N,
               T const A[],
               integer LDA,
               T       R[],
               T       C[],
               integer maxIter,
               T       epsi );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  void
  triTikhonov( integer N,
               T const Tmat[],
               integer LDT,
               integer nrhs,
               T       RHS[],
               integer ldRHS,
               T       lambda );

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
  outMATRIX( MatrixType const &         MT,
             integer                    NR,
             integer                    NC,
             T const                    A[],
             integer                    LDA,
             std::basic_ostream<char> & s,
             integer                    prec = 4,
             integer                    rperm[] = nullptr,
             integer                    cperm[] = nullptr ) {
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
  outMAPLE( integer                    NR,
            integer                    NC,
            T const                    A[],
            integer                    LDA,
            std::basic_ostream<char> & s ) {
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

#endif

///
/// eof: Alglin_aux.hh
///

