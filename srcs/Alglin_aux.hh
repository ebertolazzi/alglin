/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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
  #ifndef ALGLIN_USE_ACCELERATE
  extern "C" {
  
    integer
    LAPACKNAME(slarnv)( integer * IDIST,
                        integer * ISEED,
                        integer * N,
                        real      X[] ) ;

    integer
    LAPACKNAME(dlarnv)( integer  * IDIST,
                        integer  * ISEED,
                        integer  * N,
                        doublereal X[]) ;
  }
  #endif

  inline
  integer
  larnv( integer IDIST,
         integer ISEED[4],
         integer N,
         real    X[] ) {
    #if defined(ALGLIN_USE_OPENBLAS)
    return LAPACK_NAME(slarnv)( &IDIST, ISEED, &N, X ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    return CLAPACKNAME(slarnv)( &IDIST, ISEED, &N, X ) ;
    #else
    return LAPACKNAME(slarnv)( &IDIST, ISEED, &N, X ) ;
    #endif
  }

  inline
  integer
  lanrv( integer    IDIST,
         integer    ISEED[4],
         integer    N,
         doublereal X[] ) {
    #if defined(ALGLIN_USE_OPENBLAS)
    return LAPACK_NAME(dlarnv)( &IDIST, ISEED, &N, X ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    return CLAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X ) ;
    #else
    return LAPACKNAME(dlarnv)( &IDIST, ISEED, &N, X ) ;
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
    ASSERT( LWORK >= 2*N, "large, LWORK = " << LWORK << " must be >= " << 2*N ) ;
    // Test the input arguments
    if ( N < 0 ) return -1 ;
    if ( LDA < std::max( 1, N ) ) return -3 ;
    // pre- and post-multiply A by random orthogonal matrix
    for ( integer i = N-1 ; i >= 0 ; --i ) {
      // generate random reflection
      integer Ni = N-i ;
      lanrv( 3, ISEED, N, WORK ) ;
      T WN = nrm2( Ni, WORK, 1 ) ;
      T WA = WN ;
      if ( WORK[0] < 0 ) WA = -WN  ;
      T TAU = 0 ;
      if ( WN > 0 ) {
        T WB = WORK[0] + WA ;
        scal( Ni-1, 1.0/WB, WORK+1, 1 ) ;
        WORK[0] = 1 ;
        TAU = WB / WA ;
      }
      // multiply A(i:n,1:n) by random reflection from the left
      gemv( Transposition::TRANSPOSE,
            Ni, N, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+N, 1 ) ;
      ger( Ni, N, -TAU, WORK, 1, WORK+N, 1, A+i, LDA ) ;
      // multiply A(1:n,i:n) by random reflection from the right
      gemv( Transposition::NO_TRANSPOSE,
            N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 ) ;
      ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA ) ;
    }
    return 0 ;
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

    if ( LR == SideMultiply::LEFT ) {
      ASSERT( LWORK >= 2*M + N , "laror, LWORK = " << LWORK << " must be >= " << 2*M + N ) ;
    } else if ( LR == SideMultiply::RIGHT ) {
      ASSERT( LWORK >= 2*N + M, "laror, LWORK = " << LWORK << " must be >= " << 2*N + M ) ;
    } else {
      return -1 ;
    }

    if ( N == 0 || M == 0 ) return 0 ;

    // Check for argument errors.
    if ( M   < 0 ) return -3 ;
    if ( N   < 0 ) return -4 ;
    if ( LDA < M ) return -6 ;

    // Initialize A to the identity matrix if desired
    if ( init ) geid( M, N, A, LDA ) ;

    ASSERT( LWORK >= 2*N, "large, LWORK = " << LWORK << " must be >= " << 2*N ) ;
    // Test the input arguments
    if ( N < 0 ) return -1 ;
    if ( LDA < std::max( 1, N ) ) return -3 ;
    
    if ( LR == SideMultiply::LEFT ) {
      // pre-multiply A by random orthogonal matrix
      for ( integer i = M-1 ; i >= 0 ; --i ) {
        // generate random reflection
        integer Mi = M-i ;
        lanrv( 1, ISEED, M, WORK ) ;
        T WN = nrm2( Mi, WORK, 1 ) ;
        T WA = WN ;
        if ( WORK[0] < 0 ) WA = -WN  ;
        T TAU = 0 ;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA ;
          scal( Mi-1, 1.0/WB, WORK+1, 1 ) ;
          WORK[0] = 1 ;
          TAU = WB / WA ;
        }
        // multiply A(i:n,1:n) by random reflection from the left
        gemv( Transposition::TRANSPOSE,
              Mi, M, 1.0, A+i, LDA, WORK, 1, 0.0, WORK+M, 1 ) ;
        ger( Mi, N, -TAU, WORK, 1, WORK+M, 1, A+i, LDA ) ;
      }
    } else {
      // post-multiply A by random orthogonal matrix
      for ( integer i = N-1 ; i >= 0 ; --i ) {
        // generate random reflection
        integer Ni = N-i ;
        lanrv( 1, ISEED, N, WORK ) ;
        T WN = nrm2( Ni, WORK, 1 ) ;
        T WA = WN ;
        if ( WORK[0] < 0 ) WA = -WN  ;
        T TAU = 0 ;
        if ( WN > 0 ) {
          T WB = WORK[0] + WA ;
          scal( Ni-1, 1.0/WB, WORK+1, 1 ) ;
          WORK[0] = 1 ;
          TAU = WB / WA ;
        }
        // multiply A(1:n,i:n) by random reflection from the right
        gemv( Transposition::NO_TRANSPOSE,
              N, Ni, 1.0, A+i*LDA, LDA, WORK, 1, 0.0, WORK+N, 1 ) ;
        ger( N, Ni, -TAU, WORK+N, 1, WORK, 1, A+i*LDA, LDA ) ;
      }
    }

    return 0 ;
  }

} // end namespace alglin

#endif

///
/// eof: Alglin_aux.hh
///

