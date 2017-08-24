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
/// file: alglin_blas_ge.hxx
///

namespace alglin {

  /*
  //    __ _  ___ _ __
  //   / _` |/ _ \ '__|
  //  | (_| |  __/ |
  //   \__, |\___|_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGER   performs the rank 1 operation
   *
   *     A := alpha*x*y' + A,
   *
   *  where alpha is a scalar, x is an m element vector, y is an n element
   *  vector and A is an m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of the matrix A.
   *           M must be at least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( m - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the m
   *           element vector x.
   *           Unchanged on exit.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  Y      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCY ) ).
   *           Before entry, the incremented array Y must contain the n
   *           element vector y.
   *           Unchanged on exit.
   *
   *  INCY   - INTEGER.
   *           On entry, INCY specifies the increment for the elements of
   *           Y. INCY must not be zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry, the leading m by n part of the array A must
   *           contain the matrix of coefficients. On exit, A is
   *           overwritten by the updated matrix.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Level 2 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    using namespace blas_type ;
    void
    BLASNAME(sger)( integer          const * M,
                    integer          const * N,
                    single_precision const * ALPHA,
                    single_precision const   X[],
                    integer          const * INCX,
                    single_precision const   Y[],
                    integer          const * INCY,
                    single_precision         A[],
                    integer          const * LDA ) ;
    void
    BLASNAME(dger)( integer          const * M,
                    integer          const * N,
                    double_precision const * ALPHA,
                    double_precision const   X[],
                    integer          const * INCX,
                    double_precision const   Y[],
                    integer          const * INCY,
                    double_precision         A[],
                    integer          const * LDA ) ;
  }
  #endif

  inline
  void
  ger( integer    M,
       integer    N,
       real       ALPHA,
       real const X[],
       integer    INCX,
       real const Y[],
       integer    INCY,
       real       A[],
       integer    LDA )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(sger)( CblasColMajor, M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) ; }
  #else
  { BLASNAME(sger)( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ) ; }
  #endif

  inline
  void
  ger( integer          M,
       integer          N,
       doublereal       ALPHA,
       doublereal const X[],
       integer          INCX,
       doublereal const Y[],
       integer          INCY,
       doublereal       A[],
       integer          LDA )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dger)( CblasColMajor, M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) ; }
  #else
  { BLASNAME(dger)( &M, &N, &ALPHA, X, &INCX, Y, &INCY, A, &LDA ) ; }
  #endif

  /*
  //    __ _  ___ _ __ _____   __
  //   / _` |/ _ \ '_ ` _ \ \ / /
  //  | (_| |  __/ | | | | \ V /
  //   \__, |\___|_| |_| |_|\_/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGEMV  performs one of the matrix-vector operations
   *
   *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
   *
   *  where alpha and beta are scalars, x and y are vectors and A is an
   *  m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the operation to be performed as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
   *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
   *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of the matrix A.
   *           M must be at least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry, the leading m by n part of the array A must
   *           contain the matrix of coefficients.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of DIMENSION at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
   *           and at least
   *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
   *           Before entry, the incremented array X must contain the
   *           vector x.
   *           Unchanged on exit.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  BETA   - DOUBLE PRECISION.
   *           On entry, BETA specifies the scalar beta. When BETA is
   *           supplied as zero then Y need not be set on input.
   *           Unchanged on exit.
   *
   *  Y      - DOUBLE PRECISION array of DIMENSION at least
   *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
   *           and at least
   *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
   *           Before entry with BETA non-zero, the incremented array Y
   *           must contain the vector y. On exit, Y is overwritten by the
   *           updated vector y.
   *
   *  INCY   - INTEGER.
   *           On entry, INCY specifies the increment for the elements of
   *           Y. INCY must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    using namespace blas_type ;
    void
    BLASNAME(sgemv)( character        const   TRANS[],
                     integer          const * M,
                     integer          const * N,
                     single_precision const * ALPHA,
                     single_precision const   A[],
                     integer          const * LDA,
                     single_precision const   X[],
                     integer          const * INCX,
                     single_precision const * BETA,
                     single_precision         Y[],
                     integer          const * INCY ) ;
    void
    BLASNAME(dgemv)( character        const   TRANS[],
                     integer          const * M,
                     integer          const * N,
                     double_precision const * ALPHA,
                     double_precision const   A[],
                     integer          const * LDA,
                     double_precision const   X[],
                     integer          const * INCX,
                     double_precision const * BETA,
                     double_precision         Y[],
                     integer          const * INCY ) ;
    }
  #endif

  inline
  void
  gemv( Transposition const & TRANS,
        integer               M,
        integer               N,
        real                  ALPHA,
        real const            A[],
        integer               LDA,
        real const            X[],
        integer               INCX,
        real                  BETA,
        real                  Y[],
        integer               INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(sgemv)( CblasColMajor, trans_cblas[TRANS],
                      M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) ; }
  #else
  { BLASNAME(sgemv)( trans_blas[TRANS],
                     &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ) ; }
  #endif

  inline
  void
  gemv( Transposition const & TRANS,
        integer               M,
        integer               N,
        doublereal            ALPHA,
        doublereal const      A[],
        integer               LDA,
        doublereal const      X[],
        integer               INCX,
        doublereal            BETA,
        doublereal            Y[],
        integer               INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dgemv)( CblasColMajor, trans_cblas[TRANS],
                      M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) ; }
  #else
  { BLASNAME(dgemv)( trans_blas[TRANS],
                     &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY ) ; }
  #endif

  /*
  //    __ _  ___ _ __ ___  _ __ ___
  //   / _` |/ _ \ '_ ` _ \| '_ ` _ \
  //  | (_| |  __/ | | | | | | | | | |
  //   \__, |\___|_| |_| |_|_| |_| |_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGEMM  performs one of the matrix-matrix operations
   *
   *     C := alpha*op( A )*op( B ) + beta*C,
   *
   *  where  op( X ) is one of
   *
   *     op( X ) = X   or   op( X ) = X',
   *
   *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
   *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   *
   *  Parameters
   *  ==========
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n',  op( A ) = A.
   *              TRANSA = 'T' or 't',  op( A ) = A'.
   *              TRANSA = 'C' or 'c',  op( A ) = A'.
   *
   *           Unchanged on exit.
   *
   *  TRANSB - CHARACTER*1.
   *           On entry, TRANSB specifies the form of op( B ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSB = 'N' or 'n',  op( B ) = B.
   *              TRANSB = 'T' or 't',  op( B ) = B'.
   *              TRANSB = 'C' or 'c',  op( B ) = B'.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry,  M  specifies  the number  of rows  of the  matrix
   *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry,  N  specifies the number  of columns of the matrix
   *           op( B ) and the number of columns of the matrix C. N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  K      - INTEGER.
   *           On entry,  K  specifies  the number of columns of the matrix
   *           op( A ) and the number of rows of the matrix op( B ). K must
   *           be at least  zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry, ALPHA specifies the scalar alpha.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
   *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
   *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
   *           part of the array  A  must contain the matrix  A,  otherwise
   *           the leading  k by m  part of the array  A  must contain  the
   *           matrix A.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
   *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
   *           least  max( 1, k ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
   *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
   *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
   *           part of the array  B  must contain the matrix  B,  otherwise
   *           the leading  n by k  part of the array  B  must contain  the
   *           matrix B.
   *           Unchanged on exit.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
   *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
   *           least  max( 1, n ).
   *           Unchanged on exit.
   *
   *  BETA   - DOUBLE PRECISION.
   *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
   *           supplied as zero then C need not be set on input.
   *           Unchanged on exit.
   *
   *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
   *           Before entry, the leading  m by n  part of the array  C must
   *           contain the matrix  C,  except when  beta  is zero, in which
   *           case C need not be set on entry.
   *           On exit, the array  C  is overwritten by the  m by n  matrix
   *           ( alpha*op( A )*op( B ) + beta*C ).
   *
   *  LDC    - INTEGER.
   *           On entry, LDC specifies the first dimension of C as declared
   *           in  the  calling  (sub)  program.   LDC  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *
   *  Level 3 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    using namespace blas_type ;
    void
    BLASNAME(sgemm)( character        const   TRANSA[],
                     character        const   TRANSB[],
                     integer          const * M,
                     integer          const * N,
                     integer          const * K,
                     single_precision const * ALPHA,
                     single_precision const   A[],
                     integer          const * LDA,
                     single_precision const   B[],
                     integer          const * LDB,
                     single_precision const * BETA,
                     single_precision         C[],
                     integer          const * LDC ) ;
    void
    BLASNAME(dgemm)( character        const   TRANSA[],
                     character        const   TRANSB[],
                     integer          const * M,
                     integer          const * N,
                     integer          const * K,
                     double_precision const * ALPHA,
                     double_precision const   A[],
                     integer          const * LDA,
                     double_precision const   B[],
                     integer          const * LDB,
                      double_precision const * BETA,
                     double_precision         C[],
                     integer          const * LDC ) ;
    }
  #endif

  inline
  void
  gemm( Transposition const & TRANSA,
        Transposition const & TRANSB,
        integer               M,
        integer               N,
        integer               K,
        real                  ALPHA,
        real            const A[],
        integer               LDA,
        real            const B[],
        integer               LDB,
        real                  BETA,
        real                  C[],
        integer               LDC )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(sgemm)( CblasColMajor,
                      trans_cblas[TRANSA],
                      trans_cblas[TRANSB],
                      M, N, K,
                      ALPHA, A, LDA,
                      B, LDB,
                      BETA, C, LDC ) ; }
  #else
  { BLASNAME(sgemm)( trans_blas[TRANSA],
                     trans_blas[TRANSB],
                     &M, &N, &K,
                     &ALPHA, A, &LDA,
                     B, &LDB,
                     &BETA, C, &LDC ) ; }
  #endif

  inline
  void
  gemm( Transposition const & TRANSA,
        Transposition const & TRANSB,
        integer               M,
        integer               N,
        integer               K,
        doublereal            ALPHA,
        doublereal const      A[],
        integer               LDA,
        doublereal const      B[],
        integer               LDB,
        doublereal            BETA,
        doublereal            C[],
        integer               LDC )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dgemm)( CblasColMajor,
                      trans_cblas[TRANSA],
                      trans_cblas[TRANSB],
                      M, N, K,
                      ALPHA, A, LDA,
                      B, LDB,
                      BETA, C, LDC ) ; }
  #else
  { BLASNAME(dgemm)( const_cast<character*>(trans_blas[TRANSA]),
                     const_cast<character*>(trans_blas[TRANSB]),
                     &M, &N, &K,
                     &ALPHA, A, &LDA,
                     B, &LDB,
                     &BETA, C, &LDC ) ; }
  #endif

}
