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

#include "Alglin.hh"
#include <vector>
#include <limits>

// workaround for windows
#ifdef _MSC_VER
  #ifdef max
    #undef max
  #endif
  #ifdef min
    #undef min
  #endif
#endif

#ifdef ALGLIN_OS_WINDOWS
#include <cstdlib>
static
char * basename(char *path) {
   char *drive = nullptr;
   char *dir   = nullptr;
   char *fname = nullptr;
   char *ext   = nullptr;
   _splitpath(path, drive, dir, fname, ext);
  return dir;
}
#else
#include <libgen.h>
#endif

#include <algorithm>

namespace alglin {

  #ifdef ALGLIN_USE_CBLAS
    CBLAS_TRANSPOSE trans_cblas[3] = { CblasNoTrans, CblasTrans, CblasConjTrans } ;
    CBLAS_UPLO      uplo_cblas[2]  = { CblasUpper, CblasLower } ;
    CBLAS_DIAG      diag_cblas[2]  = { CblasUnit, CblasNonUnit } ;
    CBLAS_SIDE      side_cblas[2]  = { CblasLeft, CblasRight } ;
  #endif

  character const *trans_blas[3]   = { "N", "T", "C" } ;
  character const *uplo_blas[2]    = { "Upper", "Lower" } ;
  character const *diag_blas[2]    = { "Unit",  "NonUnit" } ;
  character const *side_blas[2]    = { "Left",  "Right" } ;
  character const *balance_blas[4] = { "N",
                                       "Permute Only",
                                       "Scale Only",
                                       "Both permute and scale" } ;

  character const *job_blas[4]   = { "A", "S", "O", "N" } ;
  character const *sense_blas[4] = { "N", "E", "V", "B" } ;

  character const *direct_blas[2] = { "Forward", "Backward" } ;
  character const *store_blas[2]  = { "Row", "Column" } ;

  character const *mtype_blas[5]  = {
    "G", // A is a full matrix.
    "L", // A is a lower triangular matrix.
    "U", // A is an upper triangular matrix.
    "H", // A is an upper Hessenberg matrix.
    "B"  // A is a symmetric band matrix with lower bandwidth KL
  } ;

  character const *equilibrate_blas[4]  = { "N", "R", "C", "B" } ;

  /*
  //    ____                _              _       
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___ 
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
  */

  //! \cond NODOC
  #ifdef ALGLIN_OS_WINDOWS
    unsigned _int64 lNaN = ((unsigned _int64) 1 << 63) - 1;
    doublereal const NaN = (*(doublereal*)&lNaN) ;
  #else
    doublereal const NaN = __builtin_nan("") ;
  #endif
  //! \endcond

  doublereal const machineEps     = std::numeric_limits<doublereal>::epsilon() ;
  doublereal const sqrtMachineEps = sqrt(std::numeric_limits<doublereal>::epsilon()) ;
  doublereal const maximumValue   = std::numeric_limits<doublereal>::max() ;
  doublereal const minimumValue   = std::numeric_limits<doublereal>::min() ;

  //============================================================================

  /*    __                       _ _   _       _   _ 
  //   / _| ___  _   _ _ __   __| | \ | | __ _| \ | |
  //  | |_ / _ \| | | | '_ \ / _` |  \| |/ _` |  \| |
  //  |  _| (_) | |_| | | | | (_| | |\  | (_| | |\  |
  //  |_|  \___/ \__,_|_| |_|\__,_|_| \_|\__,_|_| \_|
  */
  //! check if the vector `pv` os size `DIM` contains only regular floats
  bool
  foundNaN( doublereal const pv[], integer DIM ) {
    for ( integer i = 0 ; i < DIM ; ++i )
      if ( !isRegular(pv[i]) )
        return true ;
    return false ;
  }

  bool
  foundNaN( real const pv[], integer DIM ) {
    for ( integer i = 0 ; i < DIM ; ++i )
      if ( !isRegular(pv[i]) )
        return true ;
    return false ;
  }

  /*       _               _    _   _       _   _ 
  //   ___| |__   ___  ___| | _| \ | | __ _| \ | |
  //  / __| '_ \ / _ \/ __| |/ /  \| |/ _` |  \| |
  // | (__| | | |  __/ (__|   <| |\  | (_| | |\  |
  //  \___|_| |_|\___|\___|_|\_\_| \_|\__,_|_| \_|
  */
  
  #define LINE_LINE_LINE_LINE "--------------------------------------------------------------------------------"

  //! check if the vector `pv` os size `DIM` contains only regular floats. If not an error is issued
  void
  checkNaN( doublereal const pv[],
            char       const v_name[],
            integer          DIM,
            integer          line,
            char       const file[] ) {
    for ( integer i = 0 ; i < DIM ; ++i ) {
      if ( isInfinite(pv[i]) ) {
        ALGLIN_ERROR( LINE_LINE_LINE_LINE <<
                      "\n(" << basename(const_cast<char*>(file)) <<
                      ':' << line <<
                      ") found Infinity at " << v_name << "[" << i << "]\n" <<
                      LINE_LINE_LINE_LINE )
      } else if ( isNaN(pv[i]) ) {
        ALGLIN_ERROR( LINE_LINE_LINE_LINE <<
                      "\n(" << basename(const_cast<char*>(file)) <<
                      ':' << line <<
                      ") found NaN at " << v_name << "[" << i << "]\n" <<
                      LINE_LINE_LINE_LINE ) ;
      }
    }
  }

  void
  checkNaN( real const pv[],
            char const v_name[],
            integer    DIM,
            integer    line,
            char const file[] ) {
    for ( integer i = 0 ; i < DIM ; ++i ) {
      if ( isInfinite(pv[i]) ) {
        ALGLIN_ERROR( LINE_LINE_LINE_LINE <<
                      "\n(" << basename(const_cast<char*>(file)) <<
                      ':' << line <<
                      ") found Infinity at " << v_name << "[" << i << "]\n" <<
                      LINE_LINE_LINE_LINE )
      } else if ( isNaN(pv[i]) ) {
        ALGLIN_ERROR( LINE_LINE_LINE_LINE <<
                      "\n(" << basename(const_cast<char*>(file)) <<
                      ':' << line <<
                      ") found NaN at " << v_name << "[" << i << "]\n" <<
                      LINE_LINE_LINE_LINE ) ;
      }
    }
  }

  bool
  foundNaN( doublereal pv[], integer DIM ) ;

  bool
  foundNaN( real pv[], integer DIM ) ;

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

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
  gtx( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] ) {
    // LU DECOMPOSITION, COLUMN INTERCHANGES
    REAL * Ajj = A ;
    for ( int j = 0 ; j < M ; Ajj += LDA+1 ) {
      integer MX = iamax( N-j, Ajj, LDA ) ;
      IPIV[j] = MX + j ; // C-based
      if ( j < IPIV[j] ) swap( M, A + j*LDA, 1, A + IPIV[j]*LDA, 1 ) ;
      if ( isZero(Ajj[0]) ) return j ;
      REAL COLM = 1/Ajj[0] ;
      ++j ;
      scal(M-j, COLM, Ajj+1, 1) ;
      ger(M-j, N-j, -1.0, Ajj+1, 1, Ajj+LDA, LDA, Ajj+LDA+1, LDA) ;
    }
    return 0 ;
  }

  template <typename REAL>
  integer
  getrx( integer M,      // NUMBER OF ROWS IN A, UNCHANGED
         integer N,      // NUMBER OF COLUMNS IN A, UNCHANGED
         REAL    A[],    // A - ARRAY OF DIMENSION (LDA, N), OVERWRITTEN BY FACTORIZATION
                         // L - UNIT LOVER TRIANGULAR
                         // U - UPPER TRIANGULAR
         integer LDA,    // FIRST DIMENSION OF A AS DECLARED IN THE CALLING (SUB)PROGRAM, UNCHANGED
         integer IPIV[], // ARRAY OF DIMENSION (M) ON EXIT CONTAINS PIVOT INDICES
         integer MB
  ) {
    // COMPUTES LU FACTORIZATION OF M-BY-N MATRIX;
    // PARTIAL PIVOTING COLUMN INTERCHANGES
    // INFO - ERROR INDICATOR, 0 - NORMAL RETURN, POSITIVE VALUE (K) INDICATE
    // THAT U(K, K) = 0.E0 EXACTLY, ALWAYS CHECK AFTER CALL
    if ( M == 0 || N == 0 ) return 0 ;
    REAL * Ajj = A ;
    for ( integer j = 0 ; j < M ; j += MB, Ajj += MB*(LDA+1) ) {
      integer JB = min_index(M-j, MB) ;
      // FACTORIZE DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR SINGULARITY
      integer INFO = gtx( JB, N-j, Ajj, LDA, IPIV+j ) ;
      if ( INFO != 0 ) return INFO + j ;
      // APPLY INTERCHANGES TO PREVIOUS BLOCKS
      integer jjB = j + JB ;
      REAL * Aj = A+jjB ;
      for ( integer i = j ; i < jjB ; ++i ) {
        integer IP = (IPIV[i] += j) ; ;
        if ( i < IP ) {
          swap( j,     A  + i*LDA, 1, A  + IP*LDA, 1 ) ;
          swap( M-jjB, Aj + i*LDA, 1, Aj + IP*LDA, 1 ) ;
        }
      }
      // COMPUTE SUPERDIAGONAL BLOCK OF U
      trsm( RIGHT,
            UPPER,
            NO_TRANSPOSE,
            NON_UNIT,
            M-jjB, JB, 1, Ajj, LDA, Ajj+JB, LDA ) ;
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm( NO_TRANSPOSE,
            NO_TRANSPOSE,
            M-jjB, N-jjB, JB,
            -1.0, Ajj+JB,         LDA,
                  Ajj+JB*LDA,     LDA,
            1.0,  Ajj+JB*(LDA+1), LDA ) ;
    }
    return 0 ;
  }

  template <typename REAL>
  integer
  gty( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] ) {
    // LU DECOMPOSITION, ROW INTERCHANGES
    REAL * Ajj = A ;
    for ( int j = 0 ; j < N ; Ajj += LDA+1 ) {
      integer MX = iamax( M-j, Ajj, 1 ) ;
      IPIV[j] = MX + j ; // C-based
      if ( j < IPIV[j] ) swap( N, A + j, LDA, A + IPIV[j], LDA ) ;
      if ( isZero(Ajj[0]) ) return j ;
      REAL ROWM = 1/Ajj[0] ;
      ++j ;
      scal(M-j, ROWM, Ajj+1, 1) ;
      ger(M-j, N-j, -1.0, Ajj+1, 1, Ajj+LDA, LDA, Ajj+LDA+1, LDA ) ;
    }
    return 0 ;
  }

  template <typename REAL>
  integer
  getry( integer M,      // NUMBER OF ROWS IN A, UNCHANGED
         integer N,      // NUMBER OF COLUMNS IN A, UNCHANGED
         REAL    A[],    // A - ARRAY OF DIMENSION (LDA, N), OVERWRITTEN BY FACTORIZATION
                         // L - UNIT LOVER TRIANGULAR
                         // U - UPPER TRIANGULAR
         integer LDA,    // FIRST DIMENSION OF A AS DECLARED IN THE CALLING (SUB)PROGRAM, UNCHANGED
         integer IPIV[], // ARRAY OF DIMENSION (M) ON EXIT CONTAINS PIVOT INDICES
         integer NB
  ) {
    // COMPUTES LU FACTORIZATION OF M-BY-N MATRIX;
    // PARTIAL PIVOTING ROW INTERCHANGES
    // INFO - ERROR INDICATOR, 0 - NORMAL RETURN, POSITIVE VALUE (K) INDICATE
    // THAT U(K, K) = 0.E0 EXACTLY, ALWAYS CHECK AFTER CALL
    if ( M == 0 || N == 0 ) return 0 ;
    REAL * Ajj = A ;
    for ( integer j = 0 ; j < N ; j += NB, Ajj += NB*(LDA+1) ) {
      integer JB = min_index(N-j, NB) ;
      // FACTORIZE DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR SINGULARITY
      integer INFO = gty( M-j, JB, Ajj, LDA, IPIV+j ) ;
      // APPLY INTERCHANGES TO PREVIOUS BLOCKS
      integer jjB = j+JB ;
      REAL * Aj = A+jjB*LDA ;
      for ( integer i = j ; i < jjB ; ++i ) {
        integer IP = (IPIV[i] += j) ;
        if ( i < IP ) {
          swap( j,     A + i,  LDA, A + IP,  LDA ) ;
          swap( N-jjB, Aj + i, LDA, Aj + IP, LDA ) ;
        }
      }
      if ( INFO != 0 ) return INFO + j ;
      // COMPUTE SUPERDIAGONAL BLOCK OF U
      trsm( LEFT,
            LOWER,
            NO_TRANSPOSE,
            UNIT,
            JB, N-jjB, 1.0, Ajj, LDA, Ajj+JB*LDA, LDA ) ;
      // UPDATE DIAGONAL AND SUBDIAGONAL BLOCKS
      gemm( NO_TRANSPOSE,
            NO_TRANSPOSE,
            M-jjB, N-jjB, JB,
            -1.0, Ajj+JB,         LDA,
                  Ajj+JB*LDA,     LDA,
            1.0,  Ajj+JB*(LDA+1), LDA ) ;
    }
    return 0 ;
  }
  
  template integer gtx( integer M, integer N, float A[],  integer LDA, integer IPIV[] ) ;
  template integer gtx( integer M, integer N, double A[], integer LDA, integer IPIV[] ) ;
  template integer gty( integer M, integer N, float A[],  integer LDA, integer IPIV[] ) ;
  template integer gty( integer M, integer N, double A[], integer LDA, integer IPIV[] ) ;

  template integer getrx( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB ) ;
  template integer getrx( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  ) ;
  template integer getry( integer M, integer N, float A[],  integer LDA, integer IPIV[], integer MB  ) ;
  template integer getry( integer M, integer N, double A[], integer LDA, integer IPIV[], integer MB  ) ;

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
               T       epsi ) {
    #if 0
    alglin::fill( M, R, 1, 1 ) ;
    for ( integer k = 0 ; k < maxIter ; ++k ) {
      for ( integer j = 0 ; j < N ; ++j ) {
        T bf(0) ;
        for ( integer i = 0 ; i < M ; ++i )
          bf += R[i]*std::abs(A[i+j*LDA]) ;
        if ( bf < M*machineEps ) bf = M*machineEps ;
        C[j] = sqrt(M/bf) ;
      }
      for ( integer i = 0 ; i < M ; ++i ) {
        T bf(0) ;
        for ( integer j = 0 ; j < N ; ++j ) {
          T tmp = std::abs(A[i+j*LDA]*C[j]) ;
          if ( bf < tmp ) bf = tmp ;
        }
        if ( bf < machineEps ) bf = machineEps ;
        R[i] = 1/bf ;
      }
    }
    T r = sqrt(alglin::absmax(M,R,1)/alglin::absmax(N,C,1)) ;
    alglin::scal(  M, r, C, 1 ) ;
    alglin::rscal( N, r, R, 1 ) ;

    return 0 ;
    #else

    for ( integer j = 0 ; j < M ; ++j ) {
      T bf(0) ;
      for ( integer i = 0 ; i < N ; ++i )
        bf = std::max( bf, std::abs(A[i+j*LDA]) ) ;
      C[j] = 1/bf ;
    }
    for ( integer i = 0 ; i < N ; ++i ) {
      T bf(0) ;
      for ( integer j = 0 ; j < M ; ++j )
        bf = std::max( bf, C[j]*std::abs(A[i+j*LDA]) ) ;
      R[i] = 1/bf ;
    }
    for ( integer k = 0 ; k < maxIter ; ++k ) {
      for ( integer j = 0 ; j < M ; ++j ) {
        T bf(0) ;
        for ( integer i = 0 ; i < N ; ++i )
          bf += R[i]*std::abs(A[i+j*LDA]) ;
        C[j] = 1/bf ;
      }
      for ( integer i = 0 ; i < N ; ++i ) {
        T bf(0) ;
        for ( integer j = 0 ; j < M ; ++j )
          bf += C[j]*std::abs(A[i+j*LDA]) ;
        R[i] = 1/bf ;
      }
    }
    T normCinf(0) ;
    for ( integer j = 0 ; j < M ; ++j ) {
      T bf(0) ;
      for ( integer i = 0 ; i < N ; ++i )
        bf = std::max( bf, R[i]*std::abs(A[i+j*LDA]) ) ;
      C[j] = 1/bf ;
      normCinf = std::max( normCinf, std::abs(C[j]) ) ;
    }
    T normRinf(0) ;
    for ( integer i = 0 ; i < N ; ++i ) {
      T bf(0) ;
      for ( integer j = 0 ; j < M ; ++j )
        bf = std::max( bf, C[j]*std::abs(A[i+j*LDA]) ) ;
      R[i] = 1/bf ;
      normRinf = std::max( normRinf, std::abs(R[i]) ) ;
    }
    T alpha = sqrt(normRinf/normCinf) ;
    alglin::scal(  M, alpha, C, 1 ) ;
    alglin::rscal( N, alpha, R, 1 ) ;
    return 0 ;
    #endif
  }

  /*\
   *
   *   RELEASE 3.0, WGS COPYRIGHT 1997.
   *
   *   PURPOSE
   *
   *   To compute (optionally) a rank-revealing QR factorization of a 
   *   real general M-by-N matrix  A,  which may be rank-deficient,
   *   and estimate its effective rank using incremental condition 
   *   estimation.
   *
   *   The routine uses a QR factorization with column pivoting:
   *      A * P = Q * R,  where  R = [ R11 R12 ],
   *                                 [  0  R22 ]
   *   with R11 defined as the largest leading submatrix whose estimated
   *   condition number is less than 1/RCOND.  The order of R11, RANK,
   *   is the effective rank of A.
   *
   *   MB03OD  does not perform any scaling of the matrix A.
   *
   *   ARGUMENTS 
   *
   *   Mode Parameters
   *
   *   JOBQR   CHARACTER*1
   *           = 'Q':  Perform a QR factorization with column pivoting;
   *           = 'N':  Do not perform the QR factorization (but assume
   *                   that it has been done outside).  
   *
   *   Input/Output Parameters
   *
   *   M       (input) INTEGER
   *           The number of rows of the matrix A.  M >= 0.
   *          
   *   N       (input) INTEGER
   *           The number of columns of the matrix A.  N >= 0.
   *          
   *   A       (input/output) DOUBLE PRECISION array, dimension 
   *           ( LDA, N )
   *           On entry with JOBQR = 'Q', the leading M by N part of this
   *           array must contain the given matrix A.
   *           On exit with JOBQR = 'Q', the leading min(M,N) by N upper
   *           triangular part of A contains the triangular factor R, 
   *           and the elements below the diagonal, with the array TAU, 
   *           represent the orthogonal matrix Q as a product of 
   *           min(M,N) elementary reflectors.
   *           On entry and on exit with JOBQR = 'N', the leading
   *           min(M,N) by N upper triangular part of A contains the
   *           triangular factor R, as determined by the QR factorization
   *           with pivoting.  The elements below the diagonal of A are 
   *           not referenced.
   *          
   *   LDA     INTEGER
   *           The leading dimension of the array A.  LDA >= max(1,M).
   *
   *   RCOND   (input) DOUBLE PRECISION
   *           RCOND is used to determine the effective rank of A, which
   *           is defined as the order of the largest leading triangular
   *           submatrix R11 in the QR factorization with pivoting of A,
   *           whose estimated condition number is less than 1/RCOND.
   *           RCOND >= 0.
   *           NOTE that when SVLMAX > 0, the estimated rank could be
   *           less than that defined above (see SVLMAX).
   *
   *   TAU     (output) DOUBLE PRECISION array, dimension ( MIN( M, N ) )
   *           On exit with JOBQR = 'Q', the leading min(M,N) elements of
   *           TAU contain the scalar factors of the elementary 
   *           reflectors.
   *           Array TAU is not referenced when JOBQR = 'N'.
   *          
   *   RANK    (output) INTEGER
   *           The effective (estimated) rank of A, i.e. the order of 
   *           the submatrix R11.
   *          
   *   SVAL    (output) DOUBLE PRECISION array, dimension ( 3 )
   *           The estimates of some of the singular values of the 
   *           triangular factor R:
   *           SVAL(1): largest singular value of R(1:RANK,1:RANK);
   *           SVAL(2): smallest singular value of R(1:RANK,1:RANK);
   *           SVAL(3): smallest singular value of R(1:RANK+1,1:RANK+1),
   *                    if RANK < MIN( M, N ), or of R(1:RANK,1:RANK),
   *                    otherwise.
   *           If the triangular factorization is a rank-revealing one
   *           (which will be the case if the leading columns were well-
   *           conditioned), then SVAL(1) will also be an estimate for
   *           the largest singular value of A, and SVAL(2) and SVAL(3)
   *           will be estimates for the RANK-th and (RANK+1)-st singular
   *           values of A, respectively.
   *           By examining these values, one can confirm that the rank
   *           is well defined with respect to the chosen value of RCOND.
   *           The ratio SVAL(1)/SVAL(2) is an estimate of the condition
   *           number of R(1:RANK,1:RANK).
   *          
   *   Workspace
   *
   *   DWORK   DOUBLE PRECISION array, dimension ( LDWORK )
   *           where LDWORK = max( 1, 2*min( M, N ) )
   *
   *   Error Indicator
   *
   *   INFO    INTEGER
   *           = 0:  successful exit
   *           < 0:  if INFO = -i, the i-th argument had an illegal 
   *                 value.
   *
   *   METHOD
   *
   *   The routine computes or uses a QR factorization with column 
   *   pivoting of A,  A * P = Q * R,  with  R  defined above, and then
   *   finds the largest leading submatrix whose estimated condition
   *   number is less than 1/RCOND, taking the possible positive value of
   *   SVLMAX into account.  This is performed using the LAPACK
   *   incremental condition estimation scheme and a slightly modified
   *   rank decision test.
   *
   *   CONTRIBUTOR
   *
   *   V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996.
   *
   *  ******************************************************************
  \*/

  template <typename T>
  integer
  rankEstimate( integer   M,
                integer   N,
                T         A[],
                integer   LDA,
                T         RCOND,
                integer & RANK,
                T         SVAL[3] ) {

    integer MN = std::min(M, N) ;
    std::vector<T> Wmin(MN), Wmax(MN) ;

    // Test the input scalar arguments.
    if ( M < 0 )          return -1 ;
    if ( N < 0 )          return -2 ;
    if ( LDA < max(1,M) ) return -4 ;
    if ( RCOND < 0 )      return -5 ;

    // Quick return if possible
    RANK    = 0 ;
    SVAL[0] = 0 ;
    SVAL[1] = 0 ;
    SVAL[2] = 0 ;
    if ( MN == 0 ) return 0 ;

    // Determine RANK using incremental condition estimation
    T SMAX = std::abs( A[0] ) ;
    if ( SMAX > 0 ) {
      T SMIN   = SMAX ;
      T SMINPR = SMIN ;
      Wmin[0] = Wmax[0] = 1 ;
      while ( ++RANK < MN ) {
        T SMAXPR, S1, C1, S2, C2 ;
        T * pA0r = A + RANK * LDA ;
        T & Arr  = pA0r[RANK] ;
        laic1( 2, RANK, &Wmin.front(), SMIN, pA0r, Arr, SMINPR, S1, C1 ) ;
        laic1( 1, RANK, &Wmax.front(), SMAX, pA0r, Arr, SMAXPR, S2, C2 ) ;

        if ( SMAXPR*RCOND > SMINPR ) break ;

        for ( integer i=0 ; i < RANK ; ++i )
          { Wmin[i] *= S1 ; Wmax[i] *= S2 ; }

        Wmin[RANK] = C1 ;
        Wmax[RANK] = C2 ;
        SMIN = SMINPR ;
        SMAX = SMAXPR ;
      } 
      SVAL[0] = SMAX ;
      SVAL[1] = SMIN ;
      SVAL[2] = SMINPR ;
    }
    return 0 ;
  }

  /*\
   *
   *  Solve
   *
   *  min || T x - b ||^2 + lambda ||x||^2
   *
   *
   *  / * * * * * * \
   *  |   * * * * * |
   *  |     * * * * |
   *  |       * * * |
   *  |         * * |
   *  |           * |
   *  | x - - - - - |
   *  |   x         |
   *  |     x       |
   *  |       x     |
   *  |         x   |
   *  \           x /
   *
  \*/

  template <typename T>
  void
  triTikhonov( integer N,
               T const Amat[],
               integer ldA,
               integer nrhs,
               T       RHS[],
               integer ldRHS,
               T       lambda ) {

    std::vector<T> Tmat(N*N), line(N), r(nrhs) ;
    T C, S ;
    
    gecopy( N, N, Amat, ldA, &Tmat.front(), N ) ;

    for ( integer i = 0 ; i < N ; ++i ) {
      std::fill( line.begin()+i, line.end(), T(0) ) ;
      std::fill( r.begin(), r.end(), T(0) ) ;
      line[i] = lambda ;
      for ( integer j = i ; j < N ; ++j ) {
        T * pTjj = &Tmat[j*(N+1)] ;
        rotg( *pTjj, line[j], C, S ) ;
        if ( N-j-1 > 0 ) rot( N-j-1, pTjj+N, N, &line[j+1], 1, C, S ) ;
        rot( nrhs, RHS+j, ldRHS, &r.front(), 1, C, S ) ;
      }
    }
    // risolvo R
    trsm( LEFT,
          UPPER,
          NO_TRANSPOSE,
          NON_UNIT,
          N, nrhs, 1.0, &Tmat.front(), N, RHS, ldRHS ) ;
  }
  
  #ifdef ALGLIN_USE_OPENBLAS

  template <typename T>
  integer
  getc2_tmpl( integer N,
              T       A[],
              integer LDA,
              integer IPIV[],
              integer JPIV[] ) {
    ALGLIN_ERROR("NOT YET IMPLEMENTED" ) ;
  }

  template <typename T>
  T
  gesc2_tmpl( integer       N,
              T       const A[],
              integer       LDA,
              T             RHS[],
              integer const IPIV[],
              integer const JPIV[] ) {
    ALGLIN_ERROR("NOT YET IMPLEMENTED" ) ;
  }

  template <typename T>
  void
  laqge_tmpl( integer             M,
              integer             N,
              T                   A[],
              integer             LDA,
              T const             R[],
              T const             C[],
              T                   ROWCND,
              T                   COLCND,
              T                   AMAX,
              EquilibrationType & equ ) {
    ALGLIN_ERROR("NOT YET IMPLEMENTED" ) ;
  }

  #endif

  template integer equilibrate( integer    M,
                                integer    N,
                                real const A[],
                                integer    LDA,
                                real       R[],
                                real       C[],
                                integer    maxIter,
                                real       epsi ) ;

  template integer equilibrate( integer          M,
                                integer          N,
                                doublereal const A[],
                                integer          LDA,
                                doublereal       R[],
                                doublereal       C[],
                                integer          maxIter,
                                doublereal       epsi ) ;

  template integer rankEstimate( integer   M,
                                 integer   N,
                                 real      A[],
                                 integer   LDA,
                                 real      RCOND,
                                 integer & RANK,
                                 real      SVAL[3] ) ;
  
  template integer rankEstimate( integer    M,
                                 integer    N,
                                 doublereal A[],
                                 integer    LDA,
                                 doublereal RCOND,
                                 integer  & RANK,
                                 doublereal SVAL[3] ) ;

  template void triTikhonov( integer    N,
                             real const Tmat[],
                             integer    LDT,
                             integer    nrhs,
                             real       RHS[],
                             integer    ldRHS,
                             real       lambda ) ;

  template void triTikhonov( integer          N,
                             doublereal const Tmat[],
                             integer          LDT,
                             integer          nrhs,
                             doublereal       RHS[],
                             integer          ldRHS,
                             doublereal       lambda ) ;
  #ifdef ALGLIN_USE_OPENBLAS
  
  template integer getc2_tmpl( integer N,
                               real    A[],
                               integer LDA,
                               integer IPIV[],
                               integer JPIV[] ) ;

  template integer getc2_tmpl( integer    N,
                               doublereal A[],
                               integer    LDA,
                               integer    IPIV[],
                               integer    JPIV[] ) ;

  template real gesc2_tmpl( integer       N,
                            real    const A[],
                            integer       LDA,
                            real          RHS[],
                            integer const IPIV[],
                            integer const JPIV[] ) ;

  template doublereal gesc2_tmpl( integer          N,
                                  doublereal const A[],
                                  integer          LDA,
                                  doublereal       RHS[],
                                  integer    const IPIV[],
                                  integer    const JPIV[] ) ;

  template void laqge_tmpl( integer             M,
                            integer             N,
                            real                A[],
                            integer             LDA,
                            real const          R[],
                            real const          C[],
                            real                ROWCND,
                            real                COLCND,
                            real                AMAX,
                            EquilibrationType & equ ) ;

  #endif
} // end namespace alglin

///
/// eof: alglin.cc
///

