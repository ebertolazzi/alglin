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
/// file: alglin_blas.hxx
///

namespace alglin {

  /*
  //                              __   __ _ _ _
  //    ___ ___  _ __  _   _     / /  / _(_) | |
  //   / __/ _ \| '_ \| | | |   / /  | |_| | | |
  //  | (_| (_) | |_) | |_| |  / /   |  _| | | |
  //   \___\___/| .__/ \__, | /_/    |_| |_|_|_|
  //            |_|    |___/
  */
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(scopy)( integer const * N,
                           real    const   X[],
                           integer const * INCX,
                           real            Y[],
                           integer const * INCY ) ;
    void
    LAPACK_F77NAME(dcopy)( integer    const * N,
                           doublereal const   X[],
                           integer    const * INCX,
                           doublereal         Y[],
                           integer    const * INCY ) ;
  }
  #endif

  inline
  void
  copy( integer    N,
        real const X[],
        integer    INCX,
        real       Y[],
        integer    INCY )
  #if defined(ALGLIN_USE_MKL)
  { scopy( &N, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(scopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(scopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  copy( integer          N,
        doublereal const X[],
        integer          INCX,
        doublereal       Y[],
        integer          INCY )
  #if defined(ALGLIN_USE_MKL)
  { dcopy( &N, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(dcopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(dcopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  fill( integer N,
        real    Y[],
        integer INCY,
        real    V )
  { copy( N, &V, 0, Y, INCY ) ; }

  inline
  void
  fill( integer    N,
        doublereal Y[],
        integer    INCY,
        doublereal V )
  { copy( N, &V, 0, Y, INCY ) ; }

  /*
  //   _____      ____ _ _ __
  //  / __\ \ /\ / / _` | '_ \
  //  \__ \\ V  V / (_| | |_) |
  //  |___/ \_/\_/ \__,_| .__/
  //                    |_|
  */
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(sswap)( integer const * N,
                           real            X[],
                           integer const * INCX,
                           real            Y[],
                           integer const * INCY ) ;
    void
    LAPACK_F77NAME(dswap)( integer const * N,
                           doublereal      X[],
                           integer const * INCX,
                           doublereal      Y[],
                           integer const * INCY ) ;
    void
    LAPACK_F77NAME(slaswp)( integer const * NCOL,
                            real            A[],
                            integer const * LDA,
                            integer const * K1,
                            integer const * K2,
                            integer const   IPIV[],
                            integer const * INC ) ;
    void
    LAPACK_F77NAME(dlaswp)( integer const * NCOL,
                            doublereal      A[],
                            integer const * LDA,
                            integer const * K1,
                            integer const * K2,
                            integer const   IPIV[],
                            integer const * INC ) ;
  }
  #endif

  inline
  void
  swap( integer N,
        real    X[],
        integer INCX,
        real    Y[],
        integer INCY )
  #if defined(ALGLIN_USE_MKL)
  { sswap( &N, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(sswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(sswap)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  swap( integer    N,
        doublereal X[],
        integer    INCX,
        doublereal Y[],
        integer    INCY )
  #if defined(ALGLIN_USE_MKL)
  { dswap( &N, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(dswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(dswap)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLASWP performs a series of row interchanges on the matrix A.
   *  One row interchange is initiated for each of rows K1 through K2 of A.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the matrix of column dimension N to which the row
   *          interchanges will be applied.
   *          On exit, the permuted matrix.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *
   *  K1      (input) INTEGER
   *          The first element of IPIV for which a row interchange will
   *          be done.
   *
   *  K2      (input) INTEGER
   *          The last element of IPIV for which a row interchange will
   *          be done.
   *
   *  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
   *          The vector of pivot indices.  Only the elements in positions
   *          K1 through K2 of IPIV are accessed.
   *          IPIV(K) = L implies rows K and L are to be interchanged.
   *
   *  INCX    (input) INTEGER
   *          The increment between successive values of IPIV.  If IPIV
   *          is negative, the pivots are applied in reverse order.
  \*/

  inline
  integer
  swaps( integer       NCOL,
         real          A[],
         integer       LDA,
         integer       I1,
         integer       I2,
         integer const IPIV[],
         integer       INC ) {
    integer info(0);
    #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_MKL)
      for ( integer i = I1 ; i <= I2 ; i += INC )
        swap( NCOL, A+i, LDA, A+IPIV[i]-1, LDA ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      BLASFUNC(slaswp)( &NCOL, A, &LDA, &K1, &K2,
                        const_cast<integer*>(IPIV), &INC ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      info = CLAPACKNAME(slaswp)( &NCOL, A, &LDA, &K1, &K2,
                                  const_cast<integer*>(IPIV), &INC ) ;
    #else
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      LAPACKNAME(slaswp)( &NCOL, A, &LDA, &K1, &K2,
                          const_cast<integer*>(IPIV), &INC ) ;
    #endif
    return info;
  }

  inline
  integer
  swaps( integer       NCOL,
         doublereal    A[],
         integer       LDA,
         integer       I1,
         integer       I2,
         integer const IPIV[],
         integer       INC ) {
    integer info(0);
    #if  defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_MKL)
      for ( integer i = I1 ; i <= I2 ; i += INC )
        swap( NCOL, A+i, LDA, A+IPIV[i]-1, LDA ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      BLASFUNC(dlaswp)( &NCOL, A, &LDA, &K1,&K2,
                        const_cast<integer*>(IPIV), &INC ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      info = CLAPACKNAME(dlaswp)( &NCOL, A, &LDA, &K1, &K2,
                                  const_cast<integer*>(IPIV), &INC ) ;
    #else
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      LAPACKNAME(dlaswp)( &NCOL, A, &LDA, &K1, &K2,
                          const_cast<integer*>(IPIV), &INC ) ;
    #endif
    return info ;
  }

  /*
  //                 _
  //   ___  ___ __ _| |
  //  / __|/ __/ _` | |
  //  \__ \ (_| (_| | |
  //  |___/\___\__,_|_|
  */
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(sscal)( integer const * N,
                           real    const * S,
                           real            X[],
                           integer const * INCX ) ;
    void
    LAPACK_F77NAME(dscal)( integer    const * N,
                           doublereal const * S,
                           doublereal         X[],
                           integer    const * INCX ) ;
    void
    LAPACK_F77NAME(srscl)( integer const * N,
                           real    const * SA,
                           real          * SX,
                           integer const * INCX ) ;
    void
    LAPACK_F77NAME(drscl)( integer    const * N,
                           doublereal const * SA,
                           doublereal       * SX,
                           integer    const * INCX ) ;
  }
  #endif

  inline
  void
  scal( integer N,
        real    S,
        real    X[],
        integer INCX )
  #if defined(ALGLIN_USE_MKL)
  { sscal( &N, &S, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(sscal)( N, S, X, INCX ) ; }
  #else
  { LAPACK_F77NAME(sscal)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  scal( integer    N,
        doublereal S,
        doublereal X[],
        integer    INCX )
  #if defined(ALGLIN_USE_MKL)
  { dscal( &N, &S, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(dscal)( N, S, X, INCX ) ; }
  #else
  { LAPACK_F77NAME(dscal)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  rscal( integer N,
         real    S,
         real    X[],
         integer INCX )
  #if defined(ALGLIN_USE_MKL)
  { srscl( &N, &S, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { real rS = 1/S ; CBLASNAME(sscal)( N, rS, X, INCX ) ; }
  #else
  { LAPACK_F77NAME(srscl)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  rscal( integer    N,
         doublereal S,
         doublereal X[],
         integer    INCX )
  #if defined(ALGLIN_USE_MKL)
  { drscl( &N, &S, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { doublereal rS = 1/S ; CBLASNAME(dscal)( N, rS, X, INCX ) ; }
  #else
  { LAPACK_F77NAME(drscl)( &N, &S, X, &INCX ) ; }
  #endif

  /*
  //    __ ___  ___ __  _   _
  //   / _` \ \/ / '_ \| | | |
  //  | (_| |>  <| |_) | |_| |
  //   \__,_/_/\_\ .__/ \__, |
  //             |_|    |___/
  */
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(saxpy)( integer const * N,
                           real    const * A,
                           real    const   X[],
                           integer const * INCX,
                           real            Y[],
                           integer const * INCY ) ;
    void
    LAPACK_F77NAME(daxpy)( integer    const * N,
                           doublereal const * A,
                           doublereal const   X[],
                           integer    const * INCX,
                           doublereal         Y[],
                           integer    const * INCY ) ;
  }
  #endif
  inline
  void
  axpy( integer    N,
        real       A,
        real const X[],
        integer    INCX,
        real       Y[],
        integer    INCY )
  #if defined(ALGLIN_USE_MKL)
  { saxpy( &N, &A, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(saxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(saxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  axpy( integer          N,
        doublereal       A,
        doublereal const X[],
        integer          INCX,
        doublereal       Y[],
        integer          INCY )
  #if defined(ALGLIN_USE_MKL)
  { daxpy( &N, &A, X, &INCX, Y, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(daxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { LAPACK_F77NAME(daxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
  #endif

  /*
  //   _______ _ __ ___
  //  |_  / _ \ '__/ _ \
  //   / /  __/ | | (_) |
  //  /___\___|_|  \___/
  */
  inline
  void
  zero( integer N,
        real    X[],
        integer INCX )
  #if defined(ALGLIN_USE_MKL)
  { real z = 0 ; integer iz = 0 ; scopy( &N, &z, &iz, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { real z = 0 ; CBLASNAME(scopy)( N, &z, 0, X, INCX ) ; }
  #else
  { real z = 0 ; integer iz = 0 ; LAPACK_F77NAME(scopy)( &N, &z, &iz, X, &INCX ) ; }
  #endif

  inline
  void
  zero( integer    N,
        doublereal X[],
        integer    INCX )
  #if defined(ALGLIN_USE_MKL)
  { doublereal z = 0 ; integer iz = 0 ; dcopy( &N, &z, &iz, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { doublereal z = 0 ; CBLASNAME(dcopy)( N, &z, 0, X, INCX ) ; }
  #else
  { doublereal z = 0 ; integer iz = 0 ; LAPACK_F77NAME(dcopy)( &N, &z, &iz, X, &INCX ) ; }
  #endif

  /*
  //             _
  //   _ __ ___ | |_
  //  | '__/ _ \| __|
  //  | | | (_) | |_
  //  |_|  \___/ \__|
  */
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(srot)( integer const * N,
                          real            DX[],
                          integer const * INCX,
                          real            DY[],
                          integer const * INCY,
                          real    const * C,
                          real    const * S ) ;
    void
    LAPACK_F77NAME(drot)( integer    const * N,
                          doublereal         DX[],
                          integer    const * INCX,
                          doublereal         DY[],
                          integer    const * INCY,
                          doublereal const * C,
                          doublereal const * S ) ;
  }
  #endif

  inline
  void
  rot( integer N,
       real    DX[],
       integer INCX,
       real    DY[],
       integer INCY,
       real    C,
       real    S )
  #if defined(ALGLIN_USE_MKL)
  { srot( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(srot)( N, DX, INCX, DY, INCY, C, S ) ; }
  #else
  { LAPACK_F77NAME(srot)( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
  #endif

  inline
  void
  rot( integer    N,
       doublereal DX[],
       integer    INCX,
       doublereal DY[],
       integer    INCY,
       doublereal C,
       doublereal S )
  #if defined(ALGLIN_USE_MKL)
  { drot( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(drot)( N, DX, INCX, DY, INCY, C, S ) ; }
  #else
  { LAPACK_F77NAME(drot)( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
  #endif

  /*
  //             _
  //   _ __ ___ | |_ __ _
  //  | '__/ _ \| __/ _` |
  //  | | | (_) | || (_| |
  //  |_|  \___/ \__\__, |
  //                |___/
  */
  /*\
   *   Construct the Givens transformation
   *
   *         ( DC  DS )
   *     G = (        ) ,    DC**2 + DS**2 = 1 ,
   *         (-DS  DC )
   *
   *     which zeros the second entry of the 2-vector  (DA,DB)**T .
   *
   *     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
   *     storage.  The value of DB is overwritten by a value Z which
   *     allows DC and DS to be recovered by the following algorithm.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    LAPACK_F77NAME(srotg)( real * DX,
                           real * DY,
                           real * C,
                           real * S ) ;
    void
    LAPACK_F77NAME(drotg)( doublereal * DX,
                           doublereal * DY,
                           doublereal * C,
                           doublereal * S ) ;
  }
  #endif

  inline
  void
  rotg( real & DX,
        real & DY,
        real & C,
        real & S )
  #if defined(ALGLIN_USE_MKL)
  { srotg( &DX, &DY, &C, &S ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(srotg)( &DX, &DY, &C, &S ) ; }
  #else
  { LAPACK_F77NAME(srotg)( &DX, &DY, &C, &S ) ; }
  #endif

  inline
  void
  rotg( doublereal & DX,
        doublereal & DY,
        doublereal & C,
        doublereal & S )
  #if defined(ALGLIN_USE_MKL)
  { drotg( &DX, &DY, &C, &S ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { CBLASNAME(drotg)( &DX, &DY, &C, &S ) ; }
  #else
  { LAPACK_F77NAME(drotg)( &DX, &DY, &C, &S ) ; }
  #endif

  /*
  //                        ____
  //   _ __  _ __ _ __ ___ |___ \
  //  | '_ \| '__| '_ ` _ \  __) |
  //  | | | | |  | | | | | |/ __/
  //  |_| |_|_|  |_| |_| |_|_____|
  */
  /*\
   *  DNRM2 returns the euclidean norm of a vector via the function
   *  name, so that
   *
   *     DNRM2 := sqrt( x'*x )
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    real
    LAPACK_F77NAME(snrm2)( integer const * N,
                           real    const   X[],
                           integer const * INCX ) ;

    doublereal
    LAPACK_F77NAME(dnrm2)( integer    const * N,
                           doublereal const   X[],
                           integer    const * INCX ) ;
  }
  #endif

  inline
  real
  nrm2( integer    N,
        real const X[],
        integer    INCX )
  #if defined(ALGLIN_USE_MKL)
  { return snrm2( &N, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(snrm2)( N, X, INCX ) ; }
  #else
  { return LAPACK_F77NAME(snrm2)( &N, X, &INCX ) ; }
  #endif

  inline
  doublereal
  nrm2( integer          N,
        doublereal const X[],
        integer          INCX )
  #if defined(ALGLIN_USE_MKL)
  { return dnrm2( &N, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(dnrm2)( N, X, INCX ) ; }
  #else
  { return LAPACK_F77NAME(dnrm2)( &N, X, &INCX ) ; }
  #endif

  /*
  //    __ _ ___ _   _ _ __ ___
  //   / _` / __| | | | '_ ` _ \
  //  | (_| \__ \ |_| | | | | | |
  //   \__,_|___/\__,_|_| |_| |_|
  */
  /*\
   * Purpose
   * =======
   *
   *  DASUM takes the sum of the absolute values.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    real
    LAPACK_F77NAME(sasum)( integer const * N,
                           real    const   X[],
                           integer const * INCX ) ;

    doublereal
    LAPACK_F77NAME(dasum)( integer    const * N,
                           doublereal const   X[],
                           integer    const * INCX ) ;
  }
  #endif

  inline
  real
  asum( integer    N,
        real const X[],
        integer    INCX)
  #if defined(ALGLIN_USE_MKL)
  { return sasum( &N, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(sasum)( N, X, INCX ) ; }
  #else
  { return LAPACK_F77NAME(sasum)( &N, X, &INCX ) ; }
  #endif

  inline
  doublereal
  asum( integer          N,
        doublereal const X[],
        integer          INCX)
  #if defined(ALGLIN_USE_MKL)
  { return dasum( &N, X, &INCX ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(dasum)( N, X, INCX ) ; }
  #else
  { return LAPACK_F77NAME(dasum)( &N, X, &INCX ) ; }
  #endif

  /*
  //    __ _ _ __ ___   __ ___  __
  //   / _` | '_ ` _ \ / _` \ \/ /
  //  | (_| | | | | | | (_| |>  <
  //   \__,_|_| |_| |_|\__,_/_/\_\
  */
  /*\
   *  Purpose
   *  =======
   *
   *     IDAMAX finds the index of element having max. absolute value.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    integer
    LAPACK_F77NAME(isamax)( integer const * N,
                            real    const   X[],
                            integer const * INCX ) ;

    integer
    LAPACK_F77NAME(idamax)( integer    const * N,
                            doublereal const   X[],
                            integer    const * INCX ) ;
  }
  #endif

  inline
  integer
  iamax( integer    N,
         real const X[],
         integer    INCX )
  #if defined(ALGLIN_USE_MKL)
  { return isamax( &N, X, &INCX )-1 ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return integer(CBLASNAME(isamax)( N, X, INCX )) ; }
  #else
  { return LAPACK_F77NAME(isamax)( &N, X, &INCX )-1 ; }
  #endif

  inline
  integer
  iamax( integer          N,
         doublereal const X[],
         integer          INCX )
  #if defined(ALGLIN_USE_MKL)
  { return idamax( &N, X, &INCX )-1 ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return integer(CBLASNAME(idamax)( N, X, INCX )) ; }
  #else
  { return LAPACK_F77NAME(idamax)( &N, X, &INCX )-1 ; }
  #endif

  inline
  real
  absmax( integer    N,
          real const X[],
          integer    INCX )
  { real tmp = X[iamax(N,X,INCX)] ; return tmp > 0 ? tmp : -tmp ; }

  inline
  doublereal
  absmax( integer          N,
          doublereal const X[],
          integer          INCX )
  { doublereal tmp = X[iamax(N,X,INCX)] ; return tmp > 0 ? tmp : -tmp ; }

  /*
  //       _       _
  //    __| | ___ | |_
  //   / _` |/ _ \| __|
  //  | (_| | (_) | |_
  //   \__,_|\___/ \__|
  */
  /*\
   *  Purpose
   *  =======
   *
   *     DDOT forms the dot product of two vectors.
   *     uses unrolled loops for increments equal to one.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    real
    LAPACK_F77NAME(sdot)( integer const * N,
                          real    const   SX[],
                          integer const * INCX,
                          real    const   SY[],
                          integer const * INCY ) ;

    doublereal
    LAPACK_F77NAME(ddot)( integer    const * N,
                          doublereal const   SX[],
                          integer    const * INCX,
                          doublereal const   SY[],
                          integer    const * INCY ) ;
  }
  #endif

  inline
  real
  dot( integer    N,
       real const SX[],
       integer    INCX,
       real const SY[],
       integer    INCY )
  #if defined(ALGLIN_USE_MKL)
  { return sdot( &N, SX, &INCX, SY, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(sdot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return LAPACK_F77NAME(sdot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif

  inline
  doublereal
  dot( integer          N,
       doublereal const SX[],
       integer          INCX,
       doublereal const SY[],
       integer          INCY )
  #if defined(ALGLIN_USE_MKL)
  { return ddot( &N, SX, &INCX, SY, &INCY ) ; }
  #elif defined(ALGLIN_USE_CBLAS)
  { return CBLASNAME(ddot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return LAPACK_F77NAME(ddot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif

}
