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
/// file: alglin.hh
///

#ifndef ALGLIN_HH
#define ALGLIN_HH

/*
// [sd]lamch
// [sd]copy
// [sd]swap [sd]swaps
// [sd]scal
// [sd]axpy
// [sd]zero
// [sd]rot, [sd]rotg, [sd]lartg
// [sd]lapy
// [sd]ger
// [sd]gemv
// [sd]gemm
// [sd]nrm2
// [sd]asum
// [sd]amax
// [sd]dot
// [sd]trmv, [sd]trsv, [sd]trmm, [sd]trsm
// [sd]tbmv, [sd]tbsv
// [sd]gttrf, [sd]gttrs
// [sd]geadd
// [sd]getrf, [sd]getrs, [sd]gesv, [sd]getc2, [sd]gesc2
// [sd]gecon, [sd]geequ, [sd]laqge
// [sd]gecopy, [sd]gezero, [sd]geid
// [sd]normInf [sd]norm1, [sd]normF, [sd]maxabs
// [sd]lascl
// [sd]gev, [sd]ggev, [sd]ggevx (autovalori)
// [sd]gesvd, [sd]gesdd
// [sd]gelsd, [sd]gelss, [sd]gelsy
// [sd]ormqr, [sd]orm2r
// [sd]larft, [sd]larfg, [sd]larfb
// [sd]geqrf, [sd]geqr2, [sd]geqp3
// [sd]tzrzf, [sd]ormrz
// [sd]laic1
// rankEstimate
*/

#include <cmath>
#include <ctgmath>
#include <cstring>

// osx architecture
#if defined(__APPLE__) && defined(__MACH__)
  #define ALGLIN_OS_OSX 1
  #if defined(__i386__)
    #define ALGLIN_ARCH32 1
  #elif defined(__x86_64__)
    #define ALGLIN_ARCH64 1
  #endif
  #if !defined(ALGLIN_USE_ACCELERATE) && \
      !defined(ALGLIN_USE_ATLAS)      && \
      !defined(ALGLIN_USE_OPENBLAS)   && \
      !defined(ALGLIN_USE_LAPACK)
    #define ALGLIN_USE_ACCELERATE 1
  #endif
#elif defined(__unix__)
  #define ALGLIN_OS_LINUX 1
  #if defined(__i386__)
    #define ALGLIN_ARCH32 1
  #elif defined(__x86_64__)
    #define ALGLIN_ARCH64 1
  #endif
  #if !defined(ALGLIN_USE_ACCELERATE) && \
      !defined(ALGLIN_USE_ATLAS)      && \
      !defined(ALGLIN_USE_OPENBLAS)   && \
      !defined(ALGLIN_USE_LAPACK)
    #define ALGLIN_USE_ATLAS 1
  #endif
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define ALGLIN_OS_WINDOWS 1
  #if defined(_M_X64) || defined(_M_AMD64)
    #define ALGLIN_ARCH64 1
  #else
    #define ALGLIN_ARCH32 1
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #if !defined(ALGLIN_USE_ACCELERATE) && \
      !defined(ALGLIN_USE_ATLAS)      && \
      !defined(ALGLIN_USE_OPENBLAS)   && \
      !defined(ALGLIN_USE_LAPACK)
    #define ALGLIN_USE_LAPACK 1
  #endif
#else
  #error "unsupported OS!"
#endif

// if C++11
#ifdef ALGLIN_OS_WINDOWS
  #if _MSC_VER >= 1900
    #define ALGLIN_USE_CXX11
  #endif
#else
  #if __cplusplus > 199711L
    #ifndef DO_NOT_USE_CXX11
      #define ALGLIN_USE_CXX11
    #endif
  #endif
#endif

// find Headers for Lapack/Blas

#define BLASNAME(A)   A##_
#define LAPACKNAME(A) A##_

#if defined(ALGLIN_USE_ACCELERATE)

  #include <Accelerate/Accelerate.h>
  #define CBLASNAME(A)   cblas_##A
  #define CLAPACKNAME(A) A##_

#elif defined(ALGLIN_USE_ATLAS)

  // atlas 3.6.0
  extern "C" {
    #include <cblas.h>
    #include <clapack.h>
  }
  #define CBLASNAME(A)   cblas_##A
  #define CLAPACKNAME(A) clapack_##A

#elif defined(ALGLIN_USE_OPENBLAS)

  #define LAPACK_COMPLEX_CUSTOM
  #include <complex>
  #define lapack_complex_float  std::complex<float>
  #define lapack_complex_double std::complex<double>

  #include <cblas.h>
  #include <lapacke.h>

  #define CBLASNAME(A)    cblas_##A
  #define LAPACK_NAME(A)  LAPACK_##A
  #define LAPACKE_NAME(A) LAPACKE_##A

#elif defined(ALGLIN_USE_LAPACK)

  #include <cstdint>
  namespace blas_type {
    typedef char     character ;
    #ifdef LAPACK_32_BIT
    typedef int32_t  integer ; // 32 bit integer lapack
    #else
    typedef int64_t  integer ; // 64 bit integer lapack (default)
    #endif
    typedef float    single_precision ;
    typedef double   double_precision ;
  }

#else
  #error "You must select the linear algebra packages used!"
#endif

#ifndef ALGLIN_OS_WINDOWS
  #include <cstdint>
#endif

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#ifndef ALGLIN_ERROR
  #define ALGLIN_ERROR(MSG) { \
    std::ostringstream ost ; ost << MSG << '\n' ; \
    throw std::runtime_error(ost.str()) ; \
  }
#endif

#ifndef ALGLIN_ASSERT
  #define ALGLIN_ASSERT(COND,MSG) if ( !(COND) ) ALGLIN_ERROR( "in alglin::" << MSG )
#endif

namespace alglin {

  // Standard types
  #ifdef ALGLIN_OS_WINDOWS
    typedef          __int8  int8_t   ;
    typedef          __int16 int16_t  ;
    typedef          __int32 int32_t  ;
    typedef          __int64 int64_t  ;
    typedef unsigned __int8  uint8_t  ;
    typedef unsigned __int16 uint16_t ;
    typedef unsigned __int32 uint32_t ;
    typedef unsigned __int64 uint64_t ;
  #endif

  using namespace std ;

  typedef enum { NO_TRANSPOSE        = 0,
                 TRANSPOSE           = 1,
                 CONJUGATE_TRANSPOSE = 2 } Transposition ;

  typedef enum { UPPER = 0,
                 LOWER = 1 } ULselect ;

  typedef enum { UNIT     = 0,
                 NON_UNIT = 1 } DiagonalType ;

  typedef enum { LEFT  = 0,
                 RIGHT = 1 } SideMultiply ;

  typedef enum { NO_BALANCE        = 0, // 'N'
                 PERMUTE_ONLY      = 1, // 'P'
                 SCALE_ONLY        = 2, // 'S'
                 PERMUTE_AND_SCALE = 3  // 'B'
               } BalanceType ;

  typedef enum { ALL     = 0, // 'A'
                 REDUCED = 1, // 'S'
                 INPLACE = 2, // 'O'
                 NO_JOB  = 3  // 'N'
               } JobType ;

  typedef enum { NONE                         = 0, // 'N'
                 EIGENVALUES_ONLY             = 1, // 'E'
                 EIGENVECTORS_ONLY            = 2, // 'V'
                 EIGENVALUES_AND_EIGENVECTORS = 3  // 'B'
               } SenseType ;

  typedef enum { FORWARD  = 0,
                 BACKWARD = 1 } DirectionType ;

  typedef enum { COLUMNWISE = 0,
                 ROWWISE    = 1 } StorageType ;

  typedef enum { FULL_MATRIX             = 0,
                 LOWER_TRIANGULAR_MATRIX = 1,
                 UPPER_TRIANGULAR_MATRIX = 2,
                 HESSENBERG_MATRIX       = 3,
                 BANDED_MATRIX           = 4 } MatrixType ;

  typedef enum { NO_EQUILIBRATE      = 0,
                 EQUILIBRATE_ROWS    = 1,
                 EQUILIBRATE_COLUMNS = 2,
                 EQUILIBRATE_BOTH    = 3 } EquilibrationType ;

  #ifdef ALGLIN_USE_ACCELERATE
    typedef __CLPK_integer    integer ;
    typedef __CLPK_real       real ;
    typedef __CLPK_doublereal doublereal ;
    typedef char              character ;
    typedef doublereal        return_precision ;
    #define ALGLIN_USE_CBLAS
  #elif defined(ALGLIN_USE_ATLAS)
    typedef int        integer ;
    typedef float      real ;
    typedef double     doublereal ;
    typedef char       character ;
    typedef doublereal return_precision ;
    #define ALGLIN_USE_CBLAS
  #elif defined(ALGLIN_USE_OPENBLAS)
    typedef blasint  integer ;
    typedef float    real ;
    typedef double   doublereal ;
    typedef char     character ;
    typedef double   return_precision ;
    #define ALGLIN_USE_CBLAS
  #elif defined(ALGLIN_USE_LAPACK)
    typedef blas_type::integer          integer ;
    typedef blas_type::single_precision real ;
    typedef blas_type::double_precision doublereal ;
    typedef blas_type::character        character ;
    typedef doublereal return_precision ;
  #else
    #error "You must select the linear algebra packages used!"
  #endif

  #ifdef ALGLIN_USE_CBLAS
    extern CBLAS_TRANSPOSE trans_cblas[3] ;
    extern CBLAS_UPLO      uplo_cblas[2] ;
    extern CBLAS_DIAG      diag_cblas[2] ;
    extern CBLAS_SIDE      side_cblas[2] ;
  #endif

  extern character const *trans_blas[3] ;
  extern character const *uplo_blas[2] ;
  extern character const *diag_blas[2] ;
  extern character const *side_blas[2] ;
  extern character const *balance_blas[4] ;
  extern character const *job_blas[4] ;
  extern character const *sense_blas[4] ;
  extern character const *direct_blas[2] ;
  extern character const *store_blas[2] ;
  extern character const *mtype_blas[5] ;
  extern character const *equilibrate_blas[4] ;

  /*
  //    ____                _              _       
  //   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___ 
  //  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
  //  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
  //   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
  */

  extern doublereal const NaN            ; //!< Not a number constant
  extern doublereal const machineEps     ; //!< machine espilon
  extern doublereal const sqrtMachineEps ; //!< square root of machine espilon
  extern doublereal const maximumValue   ; //!< maximum floating point value
  extern doublereal const minimumValue   ; //!< minimum floating point value

  static
  inline
  bool isZero( doublereal x )
  { return FP_ZERO == std::fpclassify(x) ; }

  static
  inline
  bool isZero( real x )
  { return FP_ZERO == std::fpclassify(x) ; }

  static
  inline
  bool isInfinite( doublereal x )
  { return FP_INFINITE == std::fpclassify(x) ; }

  static
  inline
  bool isInfinite( real x )
  { return FP_INFINITE == std::fpclassify(x) ; }

  static
  inline
  bool isNaN( doublereal x )
  { return FP_NAN == std::fpclassify(x) ; }

  static
  inline
  bool isNaN( real x )
  { return FP_NAN == std::fpclassify(x) ; }

  static
  inline
  bool isRegular( doublereal x )
  { return !( FP_INFINITE == std::fpclassify(x) || FP_NAN == std::fpclassify(x) ) ; }

  static
  inline
  bool isRegular( real x )
  { return !( FP_INFINITE == std::fpclassify(x) || FP_NAN == std::fpclassify(x) ) ; }

  static
  inline
  bool isInteger( doublereal x )
  { return isZero(x-floor(x)) ; }

  static
  inline
  bool isInteger( real x )
  { return isZero(x-floor(x)) ; }

  static
  inline
  bool isUnsigned( doublereal x )
  { return isInteger(x) && x >= 0 ; }

  static
  inline
  bool isUnsigned( real x )
  { return isInteger(x) && x >= 0 ; }

  //============================================================================

  bool
  foundNaN( doublereal const pv[], integer DIM ) ;

  bool
  foundNaN( real const pv[], integer DIM ) ;

  void
  checkNaN( doublereal const pv[],
            char       const v_name[],
            integer          DIM,
            integer          line,
            char       const file[] ) ;

  void
  checkNaN( real const pv[],
            char const v_name[],
            integer    DIM,
            integer    line,
            char const file[] ) ;

  //============================================================================

  //! `m_e` the value of \f$ e \f$.
  static doublereal const m_e = 2.718281828459045235360287471352662497757 ;

  //! `m_pi` the value of \f$ \pi \f$.
  static doublereal const m_pi = 3.141592653589793238462643383279502884197 ;

  //! `m_2pi` the value of \f$ 2\pi \f$.
  static doublereal const m_2pi = 6.283185307179586476925286766559005768394 ;

  //! `m_pi_2` the value of \f$ \pi/2 \f$.
  static doublereal const m_pi_2 = 1.570796326794896619231321691639751442098 ;

  //! `m_pi_4` the value of \f$ \pi/4 \f$.
  static doublereal const m_pi_4 = 0.7853981633974483096156608458198757210492 ;

  //! `m_1_pi` the value of \f$ 1/\pi \f$.
  static doublereal const m_1_pi = 0.3183098861837906715377675267450287240689 ;

  //! `m_2_pi` the value of \f$ 2/\pi \f$.
  static doublereal const m_2_pi = 0.6366197723675813430755350534900574481378 ;
  
  //! `m_sqrtpi` the value of \f$ \sqrt{\pi} \f$.
  static doublereal const m_sqrtpi = 1.772453850905516027298167483341145182798 ;

  //! `m_2_sqrtpi` the value of \f$ 2/\sqrt{\pi} \f$.
  static doublereal const m_2_sqrtpi = 1.128379167095512573896158903121545171688 ;

  //! `m_sqrt2` the value of \f$ \sqrt{2} \f$.
  static doublereal const m_sqrt2 = 1.414213562373095048801688724209698078570 ;

  //! `m_1_sqrt2` the value of \f$ 1/\sqrt{2} \f$.
  static doublereal const m_1_sqrt2 = 0.7071067811865475244008443621048490392850 ;

  //============================================================================

  static
  inline  
  integer
  min_index( integer a, integer b)
  { return a < b ? a : b ; }

  static
  inline  
  integer
  max_index( integer a, integer b)
  { return a > b ? a : b ; }

  /*
  //   _                      _
  //  | | __ _ _ __ ___   ___| |__
  //  | |/ _` | '_ ` _ \ / __| '_ \
  //  | | (_| | | | | | | (__| | | |
  //  |_|\__,_|_| |_| |_|\___|_| |_|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {
    real
    LAPACKNAME(slamch)( character * what ) ;

    doublereal
    LAPACKNAME(dlamch)( character * what ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLAMCH determines double precision machine parameters.
   *
   *  Arguments
   *  =========
   *
   *  CMACH   (input) CHARACTER*1
   *          Specifies the value to be returned by DLAMCH:
   *          = 'E' or 'e',   DLAMCH := eps
   *          = 'S' or 's ,   DLAMCH := sfmin
   *          = 'B' or 'b',   DLAMCH := base
   *          = 'P' or 'p',   DLAMCH := eps*base
   *          = 'N' or 'n',   DLAMCH := t
   *          = 'R' or 'r',   DLAMCH := rnd
   *          = 'M' or 'm',   DLAMCH := emin
   *          = 'U' or 'u',   DLAMCH := rmin
   *          = 'L' or 'l',   DLAMCH := emax
   *          = 'O' or 'o',   DLAMCH := rmax
   *
   *          where
   *
   *          eps   = relative machine precision
   *          sfmin = safe minimum, such that 1/sfmin does not overflow
   *          base  = base of the machine
   *          prec  = eps*base
   *          t     = number of (base) digits in the mantissa
   *          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
   *          emin  = minimum exponent before (gradual) underflow
   *          rmin  = underflow threshold - base**(emin-1)
   *          emax  = largest exponent before overflow
   *          rmax  = overflow threshold  - (base**emax)*(1-eps)
   *
   * =====================================================================
  \*/
  template <typename T> T lamch( character const WHAT[] ) ;

  template <>
  inline
  real
  lamch( character const WHAT[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(slamch)( const_cast<character*>(WHAT) ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return real(CLAPACKNAME(slamch)( const_cast<character*>(WHAT) )) ; }
  #else
  { return LAPACKNAME(slamch)( const_cast<character*>(WHAT) ) ; }
  #endif

  template <>
  inline
  doublereal
  lamch( character const WHAT[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(dlamch)( const_cast<character*>(WHAT) ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlamch)( const_cast<character*>(WHAT) ) ; }
  #else
  { return LAPACKNAME(dlamch)( const_cast<character*>(WHAT) ) ; }
  #endif

  /*
  //   _ _
  //  (_) | __ _  ___ _ ____   __
  //  | | |/ _` |/ _ \ '_ \ \ / /
  //  | | | (_| |  __/ | | \ V /
  //  |_|_|\__,_|\___|_| |_|\_/
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {
    integer
    LAPACKNAME(ilaenv)( integer   * ISPEC,
                        character   NAME[],
                        character   OPTS[],
                        integer   * N1,
                        integer   * N2,
                        integer   * N3,
                        integer   * N4,
                        size_t    * len_NAME,
                        size_t    * len_OPTS ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  ILAENV is called from the LAPACK routines to choose problem-dependent
   *  parameters for the local environment.  See ISPEC for a description of
   *  the parameters.
   *
   *  ILAENV returns an INTEGER
   *  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
   *  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
   *
   *  This version provides a set of parameters which should give good,
   *  but not optimal, performance on many of the currently available
   *  computers.  Users are encouraged to modify this subroutine to set
   *  the tuning parameters for their particular machine using the option
   *  and problem size information in the arguments.
   *
   *  This routine will not function correctly if it is converted to all
   *  lower case.  Converting it to all upper case is allowed.
   *
   *  Arguments
   *  =========
   *
   *  ISPEC   (input) INTEGER
   *          Specifies the parameter to be returned as the value of
   *          ILAENV.
   *          = 1: the optimal blocksize; if this value is 1, an unblocked
   *               algorithm will give the best performance.
   *          = 2: the minimum block size for which the block routine
   *               should be used; if the usable block size is less than
   *               this value, an unblocked routine should be used.
   *          = 3: the crossover point (in a block routine, for N less
   *               than this value, an unblocked routine should be used)
   *          = 4: the number of shifts, used in the nonsymmetric
   *               eigenvalue routines (DEPRECATED)
   *          = 5: the minimum column dimension for blocking to be used;
   *               rectangular blocks must have dimension at least k by m,
   *               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
   *          = 6: the crossover point for the SVD (when reducing an m by n
   *               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
   *               this value, a QR factorization is used first to reduce
   *               the matrix to a triangular form.)
   *          = 7: the number of processors
   *          = 8: the crossover point for the multishift QR method
   *               for nonsymmetric eigenvalue problems (DEPRECATED)
   *          = 9: maximum size of the subproblems at the bottom of the
   *               computation tree in the divide-and-conquer algorithm
   *               (used by xGELSD and xGESDD)
   *          =10: ieee NaN arithmetic can be trusted not to trap
   *          =11: infinity arithmetic can be trusted not to trap
   *          12 <= ISPEC <= 16:
   *               xHSEQR or one of its subroutines,
   *               see IPARMQ for detailed explanation
   *
   *  NAME    (input) CHARACTER*(*)
   *          The name of the calling subroutine, in either upper case or
   *          lower case.
   *
   *  OPTS    (input) CHARACTER*(*)
   *          The character options to the subroutine NAME, concatenated
   *          into a single character string.  For example, UPLO = 'U',
   *          TRANS = 'T', and DIAG = 'N' for a triangular routine would
   *          be specified as OPTS = 'UTN'.
   *
   *  N1      (input) INTEGER
   *  N2      (input) INTEGER
   *  N3      (input) INTEGER
   *  N4      (input) INTEGER
   *          Problem dimensions for the subroutine NAME; these may not all
   *          be required.
   *
   *  Further Details
   *  ===============
   *
   *  The following conventions have been used when calling ILAENV from the
   *  LAPACK routines:
   *  1)  OPTS is a concatenation of all of the character options to
   *      subroutine NAME, in the same order that they appear in the
   *      argument list for NAME, even if they are not used in determining
   *      the value of the parameter specified by ISPEC.
   *  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
   *      that they appear in the argument list for NAME.  N1 is used
   *      first, N2 second, and so on, and unused problem dimensions are
   *      passed a value of -1.
   *  3)  The parameter value returned by ILAENV is checked for validity in
   *      the calling subroutine.  For example, ILAENV is used to retrieve
   *      the optimal blocksize for STRTRI as follows:
   *
   *      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
   *      IF( NB.LE.1 ) NB = MAX( 1, N )
   *
   *  =====================================================================
  \*/

  inline
  integer
  ilaenv( integer         ISPEC,
          character const NAME[],
          character const OPTS[],
          integer         N1,
          integer         N2,
          integer         N3,
          integer         N4 ) {
    #ifdef ALGLIN_USE_OPENBLAS
    cerr << "Call of ilaenv in OPEN BLAS\n" ;
    exit(1) ;
    return 0 ;
    #elif !defined(ALGLIN_USE_ACCELERATE)
    size_t len_NAME = strlen( NAME ) ;
    size_t len_OPTS = strlen( OPTS ) ;
    return LAPACKNAME(ilaenv)( &ISPEC,
                               const_cast<character*>(NAME),
                               const_cast<character*>(OPTS),
                               &N1, &N2, &N3, &N4, &len_NAME, &len_OPTS ) ;
    #else
    return CLAPACKNAME(ilaenv)( &ISPEC,
                                const_cast<character*>(NAME),
                                const_cast<character*>(OPTS),
                                &N1, &N2, &N3, &N4 ) ;
    #endif
  }

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
    using namespace blas_type ;
    void
    BLASNAME(scopy)( integer          const * N,
                     single_precision const   X[],
                     integer          const * INCX,
                     single_precision         Y[],
                     integer          const * INCY ) ;
    void
    BLASNAME(dcopy)( integer          const * N,
                     double_precision const   X[],
                     integer          const * INCX,
                     double_precision         Y[],
                     integer          const * INCY ) ;
    }
  #endif

  inline
  void
  copy( integer    N,
        real const X[],
        integer    INCX,
        real       Y[],
        integer    INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(scopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(scopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  copy( integer          N,
        doublereal const X[],
        integer          INCX,
        doublereal       Y[],
        integer          INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dcopy)( N, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(dcopy)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  fill( integer    N,
        real       Y[],
        integer    INCY,
        real       V )
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
  using namespace blas_type ;
  void
  BLASNAME(sswap)( integer          const * N, 
                   single_precision         X[], 
                   integer          const * INCX, 
                   single_precision         Y[], 
                   integer          const * INCY ) ;
  void
  BLASNAME(dswap)( integer          const * N, 
                   double_precision         X[], 
                   integer          const * INCX, 
                   double_precision         Y[], 
                   integer          const * INCY ) ;
  void
  BLASNAME(slaswp)( integer          const * NCOL, 
                    single_precision         A[], 
                    integer          const * LDA, 
                    integer          const * K1, 
                    integer          const * K2,
                    integer          const   IPIV[],
                    integer          const * INC ) ;
  void
  BLASNAME(dlaswp)( integer          const * NCOL, 
                    double_precision         A[], 
                    integer          const * LDA, 
                    integer          const * K1, 
                    integer          const * K2,
                    integer          const   IPIV[],
                    integer          const * INC ) ;
  }
  #endif

  inline
  void
  swap( integer N,
        real    X[],
        integer INCX,
        real    Y[],
        integer INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(sswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(sswap)( &N, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  swap( integer    N,
        doublereal X[],
        integer    INCX,
        doublereal Y[],
        integer    INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dswap)( N, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(dswap)( &N, X, &INCX, Y, &INCY ) ; }
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
    #ifdef ALGLIN_USE_ATLAS
      for ( integer i = I1 ; i <= I2 ; i += INC ) swap( NCOL, A+i, LDA, A+IPIV[i]-1, LDA ) ;
      return 0 ;
    #elif defined(ALGLIN_USE_OPENBLAS)
      integer K1 = I1+1;
      integer K2 = I2+1;
      LAPACK_NAME(slaswp)( &NCOL, A, &LDA, &K1, &K2, IPIV, &INC ) ;
      return 0;
    #elif defined(ALGLIN_USE_ACCELERATE)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      return CLAPACKNAME(slaswp)( &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC ) ;
    #else
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      LAPACKNAME(slaswp)( &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC ) ;
      return 0;
    #endif
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
    #ifdef ALGLIN_USE_ATLAS
      for ( integer i = I1 ; i <= I2 ; i += INC ) swap( NCOL, A+i, LDA, A+IPIV[i]-1, LDA ) ;
      return 0 ;
    #elif defined(ALGLIN_USE_OPENBLAS)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      LAPACK_NAME(dlaswp)( &NCOL, A, &LDA, &K1,&K2, IPIV, &INC ) ;
      return 0;
    #elif defined(ALGLIN_USE_ACCELERATE)
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      return CLAPACKNAME(dlaswp)( &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC ) ;
    #else
      integer K1 = I1+1 ;
      integer K2 = I2+1 ;
      LAPACKNAME(dlaswp)( &NCOL, A, &LDA, &K1, &K2, const_cast<integer*>(IPIV), &INC ) ;
      return 0;
    #endif
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
  using namespace blas_type ;
  void
  BLASNAME(sscal)( integer          const * N,
                   single_precision const * S,
                   single_precision         X[],
                   integer          const * INCX ) ;
  void
  BLASNAME(dscal)( integer          const * N,
                   double_precision const * S,
                   double_precision         X[],
                   integer          const * INCX ) ;
  }
  extern "C" {
    void
    BLASNAME(srscl)( integer * N, real * SA, real * SX, integer * INCX ) ;
    void
    BLASNAME(drscl)( integer * N, doublereal * SA, doublereal * SX, integer * INCX ) ;
  }
  #endif

  inline
  void
  scal( integer N,
        real    S,
        real    X[],
        integer INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(sscal)( N, S, X, INCX ) ; }
  #else
  { BLASNAME(sscal)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  scal( integer    N,
        doublereal S,
        doublereal X[],
        integer    INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dscal)( N, S, X, INCX ) ; }
  #else
  { BLASNAME(dscal)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  rscal( integer N,
         real    S,
         real    X[],
         integer INCX )
  #ifdef ALGLIN_USE_CBLAS
  { real rS = 1/S ; CBLASNAME(sscal)( N, rS, X, INCX ) ; }
  #else
  { BLASNAME(srscl)( &N, &S, X, &INCX ) ; }
  #endif

  inline
  void
  rscal( integer    N,
         doublereal S,
         doublereal X[],
         integer    INCX )
  #ifdef ALGLIN_USE_CBLAS
  { doublereal rS = 1/S ; CBLASNAME(dscal)( N, rS, X, INCX ) ; }
  #else
  { BLASNAME(drscl)( &N, &S, X, &INCX ) ; }
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
  using namespace blas_type ;
  void
  BLASNAME(saxpy)( integer          const * N, 
                   single_precision const * A, 
                   single_precision const   X[], 
                   integer          const * INCX, 
                   single_precision         Y[], 
                   integer          const * INCY ) ;
  void
  BLASNAME(daxpy)( integer          const * N, 
                   double_precision const * A, 
                   double_precision const   X[], 
                   integer          const * INCX, 
                   double_precision         Y[], 
                   integer          const * INCY ) ;
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
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(saxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(saxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
  #endif

  inline
  void
  axpy( integer          N,
        doublereal       A,
        doublereal const X[],
        integer          INCX,
        doublereal       Y[],
        integer          INCY )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(daxpy)( N, A, X, INCX, Y, INCY ) ; }
  #else
  { BLASNAME(daxpy)( &N, &A, X, &INCX, Y, &INCY ) ; }
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
  #ifdef ALGLIN_USE_CBLAS
  { real z = 0 ; CBLASNAME(scopy)( N, &z, 0, X, INCX ) ; }
  #else
  { real z = 0 ; integer iz = 0 ; BLASNAME(scopy)( &N, &z, &iz, X, &INCX ) ; }
  #endif

  inline
  void
  zero( integer    N,
        doublereal X[],
        integer    INCX )
  #ifdef ALGLIN_USE_CBLAS
  { doublereal z = 0 ; CBLASNAME(dcopy)( N, &z, 0, X, INCX ) ; }
  #else
  { doublereal z = 0 ; integer iz = 0 ; BLASNAME(dcopy)( &N, &z, &iz, X, &INCX ) ; }
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
    BLASNAME(srot)( integer * N,
                    real      DX[],
                    integer * INCX,
                    real      DY[],
                    integer * INCY,
                    real    * C,
                    real    * S ) ;
    void
    BLASNAME(drot)( integer    * N,
                    doublereal   DX[],
                    integer    * INCX,
                    doublereal   DY[],
                    integer    * INCY,
                    doublereal * C,
                    doublereal * S  ) ;
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
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(srot)( N, DX, INCX, DY, INCY, C, S ) ; }
  #else
  { BLASNAME(srot)( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
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
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(drot)( N, DX, INCX, DY, INCY, C, S ) ; }
  #else
  { BLASNAME(drot)( &N, DX, &INCX, DY, &INCY, &C, &S ) ; }
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
    BLASNAME(srotg)( real * DX,
                     real * DY,
                     real * C,
                     real * S ) ;
    void
    BLASNAME(drotg)( doublereal * DX,
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
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(srotg)( &DX, &DY, &C, &S ) ; }
  #else
  { BLASNAME(srotg)( &DX, &DY, &C, &S ) ; }
  #endif

  inline
  void
  rotg( doublereal & DX,
        doublereal & DY,
        doublereal & C,
        doublereal & S )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(drotg)( &DX, &DY, &C, &S ) ; }
  #else
  { BLASNAME(drotg)( &DX, &DY, &C, &S ) ; }
  #endif

  /*
  //   _            _
  //  | | __ _ _ __| |_ __ _
  //  | |/ _` | '__| __/ _` |
  //  | | (_| | |  | || (_| |
  //  |_|\__,_|_|   \__\__, |
  //                   |___/
  */
  #ifndef ALGLIN_USE_ACCELERATE
  extern "C" {

    real
    LAPACKNAME(slartg)( real * F,
                        real * G,
                        real * C,
                        real * S,
                        real * R ) ;

    doublereal
    LAPACKNAME(dlartg)( doublereal * F,
                        doublereal * G,
                        doublereal * C,
                        doublereal * S,
                        doublereal * R ) ;

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLARTG generate a plane rotation so that
   *
   *     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
   *     [ -SN  CS  ]     [ G ]     [ 0 ]
   *
   *  This is a slower, more accurate version of the BLAS1 routine DROTG,
   *  with the following other differences:
   *     F and G are unchanged on return.
   *     If G=0, then CS=1 and SN=0.
   *     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
   *        floating point operations (saves work in DBDSQR when
   *        there are zeros on the diagonal).
   *
   *  If F exceeds G in magnitude, CS will be positive.
   *
   *  Arguments
   *  =========
   *
   *  F       (input) DOUBLE PRECISION
   *          The first component of vector to be rotated.
   *
   *  G       (input) DOUBLE PRECISION
   *          The second component of vector to be rotated.
   *
   *  CS      (output) DOUBLE PRECISION
   *          The cosine of the rotation.
   *
   *  SN      (output) DOUBLE PRECISION
   *          The sine of the rotation.
   *
   *  R       (output) DOUBLE PRECISION
   *          The nonzero component of the rotated vector.
   *
  \*/

  inline
  real
  lartg( real   F,
         real   G,
         real & C,
         real & S,
         real & R )
  #ifdef ALGLIN_USE_ACCELERATE
  { return CLAPACKNAME(slartg)( &F, &G, &C, &S, &R ) ; }
  #else
  { return LAPACKNAME(slartg)( &F, &G, &C, &S, &R ) ; }
  #endif

  inline
  doublereal
  lartg( doublereal   F,
         doublereal   G,
         doublereal & C,
         doublereal & S,
         doublereal & R )
  #ifdef ALGLIN_USE_ACCELERATE
  { return CLAPACKNAME(dlartg)( &F, &G, &C, &S, &R ) ; }
  #else
  { return LAPACKNAME(dlartg)( &F, &G, &C, &S, &R ) ; }
  #endif

  /*
  //  _
  // | | __ _ _ __  _   _
  // | |/ _` | '_ \| | | |
  // | | (_| | |_) | |_| |
  // |_|\__,_| .__/ \__, |
  //         |_|    |___/
  */

  #ifdef ALGLIN_USE_ATLAS
  extern "C" {

    real
    LAPACKNAME(slapy2)( real * X, real * Y ) ;

    doublereal
    LAPACKNAME(dlapy2)( doublereal * X, doublereal * Y ) ;

    real
    LAPACKNAME(slapy3)( real * X, real * Y, real * Z ) ;

    doublereal
    LAPACKNAME(dlapy3)( doublereal * X, doublereal * Y, doublereal * Z ) ;

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
   *  overflow.
   *
   *  Arguments
   *  =========
   *
   *  X       (input) DOUBLE PRECISION
   *  Y       (input) DOUBLE PRECISION
   *          X and Y specify the values x and y.
   *  =====================================================================
   *
   *  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
   *  unnecessary overflow.
   *
   *  Arguments
   *  =========
   *
   *  X       (input) DOUBLE PRECISION
   *  Y       (input) DOUBLE PRECISION
   *  Z       (input) DOUBLE PRECISION
   *          X, Y and Z specify the values x, y and z.
   *  =====================================================================
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {

    real
    LAPACKNAME(slapy2)( real * x, real * y ) ;

    real
    LAPACKNAME(slapy3)( real * x, real * y, real * z  ) ;

    doublereal
    LAPACKNAME(dlapy2)( doublereal * x, doublereal * y ) ;

    doublereal
    LAPACKNAME(dlapy3)( doublereal * x, doublereal * y, doublereal * z ) ;
  }
  #endif


  inline
  real
  lapy( real x, real y )
  #if defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(slapy2)( &x, &y ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return real(CLAPACKNAME(slapy2)( &x, &y )) ; }
  #else
  { return LAPACKNAME(slapy2)( &x, &y ) ; }
  #endif

  inline
  doublereal
  lapy( doublereal x, doublereal y, doublereal z )
  #if defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(dlapy3)( &x, &y, &z ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlapy3)( &x, &y, &z ) ; }
  #else
  { return LAPACKNAME(dlapy3)( &x, &y, &z ) ; }
  #endif

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
   *
   *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
   *
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
   *
   *              TRANSA = 'T' or 't',  op( A ) = A'.
   *
   *              TRANSA = 'C' or 'c',  op( A ) = A'.
   *
   *           Unchanged on exit.
   *
   *  TRANSB - CHARACTER*1.
   *           On entry, TRANSB specifies the form of op( B ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSB = 'N' or 'n',  op( B ) = B.
   *
   *              TRANSB = 'T' or 't',  op( B ) = B'.
   *
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
  using namespace blas_type ;
  single_precision
  BLASNAME(snrm2)( integer const * N, single_precision const X[], integer const * INCX ) ;

  double_precision
  BLASNAME(dnrm2)( integer const * N, double_precision const X[], integer const * INCX ) ;
  }
  #endif

  inline
  real
  nrm2( integer    N,
        real const X[],
        integer    INCX )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(snrm2)( N, X, INCX ) ; }
  #else
  { return BLASNAME(snrm2)( &N, X, &INCX ) ; }
  #endif

  inline
  doublereal
  nrm2( integer          N,
        doublereal const X[],
        integer          INCX )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(dnrm2)( N, X, INCX ) ; }
  #else
  { return BLASNAME(dnrm2)( &N, X, &INCX ) ; }
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
  using namespace blas_type ;
  single_precision
  BLASNAME(sasum)( integer const * N, single_precision const X[], integer const * INCX ) ;

  double_precision
  BLASNAME(dasum)( integer const * N, double_precision const X[], integer const * INCX ) ;
  }
  #endif

  inline
  real
  asum( integer    N,
        real const X[],
        integer    INCX)
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(sasum)( N, X, INCX ) ; }
  #else
  { return BLASNAME(sasum)( &N, X, &INCX ) ; }
  #endif

  inline
  doublereal
  asum( integer          N,
        doublereal const X[],
        integer          INCX)
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(dasum)( N, X, INCX ) ; }
  #else
  { return BLASNAME(dasum)( &N, X, &INCX ) ; }
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
  using namespace blas_type ;
  integer
  BLASNAME(isamax)( integer const * N, single_precision const X[], integer const * INCX ) ;

  integer
  BLASNAME(idamax)( integer const * N, double_precision const X[], integer const * INCX ) ;
  }
  #endif
  
  inline
  integer
  iamax( integer    N,
         real const X[],
         integer    INCX )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(isamax)( N, X, INCX ) ; }
  #else
  { return BLASNAME(isamax)( &N, X, &INCX )-1 ; }
  #endif

  inline
  integer
  iamax( integer          N,
         doublereal const X[],
         integer          INCX )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(idamax)( N, X, INCX ) ; }
  #else
  { return BLASNAME(idamax)( &N, X, &INCX )-1 ; }
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
  using namespace blas_type ;
  single_precision
  BLASNAME(sdot)( integer          const * N, 
                  single_precision const   SX[], 
                  integer          const * INCX, 
                  single_precision const   SY[], 
                  integer          const * INCY ) ;

  double_precision
  BLASNAME(ddot)( integer          const * N, 
                  double_precision const   SX[], 
                  integer          const * INCX, 
                  double_precision const   SY[], 
                  integer          const * INCY ) ;
  }
  #endif
  
  inline
  real
  dot( integer    N,
       real const SX[],
       integer    INCX,
       real const SY[],
       integer    INCY )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(sdot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return BLASNAME(sdot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif

  inline
  doublereal
  dot( integer          N,
       doublereal const SX[],
       integer          INCX,
       doublereal const SY[],
       integer          INCY )
  #ifdef ALGLIN_USE_CBLAS
  { return CBLASNAME(ddot)( N, SX, INCX, SY, INCY ) ; }
  #else
  { return BLASNAME(ddot)( &N, SX, &INCX, SY, &INCY ) ; }
  #endif

  /*
  //   _
  //  | |_ _ __ _ __ _____   __
  //  | __| '__| '_ ` _ \ \ / /
  //  | |_| |  | | | | | \ V /
  //   \__|_|  |_| |_| |_|\_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRMV  performs one of the matrix-vector operations
   *
   *     x := A*x,   or   x := A**T*x,
   *
   *  where x is an n element vector and  A is an n by n unit, or non-unit,
   *  upper or lower triangular matrix.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the operation to be performed as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   x := A*x.
   *
   *              TRANS = 'T' or 't'   x := A**T*x.
   *
   *              TRANS = 'C' or 'c'   x := A**T*x.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular matrix and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular matrix and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
   *           A are not referenced either, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, n ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element vector x. On exit, X is overwritten with the
   *           tranformed vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 2 Blas routine.
   *  The vector and matrix arguments are not referenced when N = 0, or M = 0
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  BLASNAME(strmv)( character        const   UPLO[],
                   character        const   TRANS[],
                   character        const   DIAG[],
                   integer          const * N,
                   single_precision const   A[],
                   integer          const * LDA,
                   single_precision         X[],
                   integer          const * INCX ) ;
  void
  BLASNAME(dtrmv)( character        const   UPLO[],
                   character        const   TRANS[],
                   character        const   DIAG[],
                   integer          const * N,
                   double_precision const   A[],
                   integer          const * LDA,
                   double_precision         X[],
                   integer          const * INCX ) ;
  }
  #endif

  inline
  void
  trmv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               N,
        real          const   A[],
        integer               LDA,
        real                  X[],
        integer               INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(strmv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, A, LDA, X, INCX ) ; }
  #else
  { BLASNAME(strmv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, A, &LDA, X, &INCX ) ; }
  #endif

  inline
  void
  trmv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               N,
        doublereal    const   A[],
        integer               LDA,
        doublereal            X[],
        integer               INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dtrmv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, A, LDA, X, INCX ) ; }
  #else
  { BLASNAME(dtrmv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, A, &LDA, X, &INCX ) ; }
  #endif

  /*
  //   _
  //  | |_ _ __ _____   __
  //  | __| '__/ __\ \ / /
  //  | |_| |  \__ \\ V /
  //   \__|_|  |___/ \_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRSV  solves one of the systems of equations
   *
   *     A*x = b,   or   A**T*x = b,
   *
   *  where b and x are n element vectors and A is an n by n unit, or
   *  non-unit, upper or lower triangular matrix.
   *
   *  No test for singularity or near-singularity is included in this
   *  routine. Such tests must be performed before calling this routine.
   *
   *  Arguments
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the equations to be solved as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   A*x = b.
   *
   *              TRANS = 'T' or 't'   A**T*x = b.
   *
   *              TRANS = 'C' or 'c'   A**T*x = b.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with  UPLO = 'U' or 'u', the leading n by n
   *           upper triangular part of the array A must contain the upper
   *           triangular matrix and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry with UPLO = 'L' or 'l', the leading n by n
   *           lower triangular part of the array A must contain the lower
   *           triangular matrix and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u', the diagonal elements of
   *           A are not referenced either, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           max( 1, n ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element right-hand side vector b. On exit, X is overwritten
   *           with the solution vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  BLASNAME(strsv)( character        const   UPLO[], 
                   character        const   TRANS[], 
                   character        const   DIAG[], 
                   integer          const * N, 
                   single_precision const   A[], 
                   integer          const * LDA, 
                   single_precision         X[], 
                   integer          const * INCX ) ;
  void
  BLASNAME(dtrsv)( character        const   UPLO[], 
                   character        const   TRANS[], 
                   character        const   DIAG[], 
                   integer          const * N, 
                   double_precision const   A[], 
                   integer          const * LDA, 
                   double_precision         X[], 
                   integer          const * INCX ) ;
  }
  #endif

  inline
  void
  trsv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               N,
        real          const   A[],
        integer               LDA,
        real                  X[],
        integer               INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(strsv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, A, LDA, X, INCX ) ; }
  #else
  { BLASNAME(strsv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, A, &LDA, X, &INCX ) ; }
  #endif

  inline
  void
  trsv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               N,
        doublereal    const   A[],
        integer               LDA,
        doublereal            X[],
        integer               INCX )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dtrsv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, A, LDA, X, INCX  ) ; }
  #else
  { BLASNAME(dtrsv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, A, &LDA, X, &INCX ) ; }
  #endif

  /*
  //   _
  //  | |_ _ __ _ __ ___  _ __ ___
  //  | __| '__| '_ ` _ \| '_ ` _ \
  //  | |_| |  | | | | | | | | | | |
  //   \__|_|  |_| |_| |_|_| |_| |_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRMM  performs one of the matrix-matrix operations
   *
   *     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
   *
   *  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
   *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   *
   *     op( A ) = A   or   op( A ) = A**T.
   *
   *  Arguments
   *  ==========
   *
   *  SIDE   - CHARACTER*1.
   *           On entry,  SIDE specifies whether  op( A ) multiplies B from
   *           the left or right as follows:
   *
   *              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
   *
   *              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
   *
   *           Unchanged on exit.
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix A is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n'   op( A ) = A.
   *
   *              TRANSA = 'T' or 't'   op( A ) = A**T.
   *
   *              TRANSA = 'C' or 'c'   op( A ) = A**T.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit triangular
   *           as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of B. M must be at
   *           least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of B.  N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
   *           zero then  A is not referenced and  B need not be set before
   *           entry.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
   *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
   *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
   *           upper triangular part of the array  A must contain the upper
   *           triangular matrix  and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
   *           lower triangular part of the array  A must contain the lower
   *           triangular matrix  and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
   *           A  are not referenced either,  but are assumed to be  unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
   *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
   *           then LDA must be at least max( 1, n ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
   *           Before entry,  the leading  m by n part of the array  B must
   *           contain the matrix  B,  and  on exit  is overwritten  by the
   *           transformed matrix.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in  the  calling  (sub)  program.   LDB  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 3 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  BLASNAME(strmm) ( character        const   SIDE[],     // "L" or "R" 
                    character        const   UPLO[],     // "U" or "L"
                    character        const   TRANSA[],   // "N", "T", "C"
                    character        const   DIAG[],     // "N", "U"
                    integer          const * M,
                    integer          const * N,
                    single_precision const * ALPHA,
                    single_precision const   A[],
                    integer          const * LDA,
                    single_precision         B[],
                    integer          const * LDB ) ;
  void
  BLASNAME(dtrmm) ( character        const   SIDE[],     // "L" or "R" 
                    character        const   UPLO[],     // "U" or "L"
                    character        const   TRANSA[],   // "N", "T", "C"
                    character        const   DIAG[],     // "N", "U"
                    integer          const * M,
                    integer          const * N,
                    double_precision const * ALPHA,
                    double_precision const   A[],
                    integer          const * LDA,
                    double_precision         B[],
                    integer          const * LDB ) ;
  }
  #endif

  inline
  void
  trmm( SideMultiply  const & SIDE,
        ULselect      const & UPLO,
        Transposition const & TRANSA,
        DiagonalType  const & DIAG,
        integer               M,
        integer               N,
        real                  alpha,
        real          const   A[],
        integer               LDA,
        real                  B[],
        integer               LDB )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(strmm) ( CblasColMajor,
                       side_cblas[SIDE],
                       uplo_cblas[UPLO],
                       trans_cblas[TRANSA],
                       diag_cblas[DIAG],
                       M, N, alpha, A, LDA, B, LDB ) ; }
  #else
  { BLASNAME(strmm)( side_blas[SIDE],
                     uplo_blas[UPLO],
                     trans_blas[TRANSA],
                     diag_blas[DIAG],
                     &M, &N, &alpha, A, &LDA, B, &LDB ) ; }
  #endif

  inline
  void
  trmm( SideMultiply  const & SIDE,
        ULselect      const & UPLO,
        Transposition const & TRANSA,
        DiagonalType  const & DIAG,
        integer               M,
        integer               N,
        doublereal            alpha,
        doublereal    const   A[],
        integer               LDA,
        doublereal            B[],
        integer               LDB )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dtrmm) ( CblasColMajor,
                       side_cblas[SIDE],
                       uplo_cblas[UPLO],
                       trans_cblas[TRANSA],
                       diag_cblas[DIAG],
                       M, N, alpha, A, LDA, B, LDB ) ;  }
  #else
  { BLASNAME(dtrmm)( side_blas[SIDE],
                     uplo_blas[UPLO],
                     trans_blas[TRANSA],
                     diag_blas[DIAG],
                     &M, &N, &alpha, A, &LDA, B, &LDB ) ; }
  #endif

  /*
  //   _
  //  | |_ _ __ ___ _ __ ___
  //  | __| '__/ __| '_ ` _ \
  //  | |_| |  \__ \ | | | | |
  //   \__|_|  |___/_| |_| |_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DTRSM  solves one of the matrix equations
   *
   *     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   *
   *  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
   *  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   *
   *     op( A ) = A   or   op( A ) = A**T.
   *
   *  The matrix X is overwritten on B.
   *
   *  Arguments
   *  ==========
   *
   *  SIDE   - CHARACTER*1.
   *           On entry, SIDE specifies whether op( A ) appears on the left
   *           or right of X as follows:
   *
   *              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
   *
   *              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
   *
   *           Unchanged on exit.
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix A is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANSA - CHARACTER*1.
   *           On entry, TRANSA specifies the form of op( A ) to be used in
   *           the matrix multiplication as follows:
   *
   *              TRANSA = 'N' or 'n'   op( A ) = A.
   *
   *              TRANSA = 'T' or 't'   op( A ) = A**T.
   *
   *              TRANSA = 'C' or 'c'   op( A ) = A**T.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit triangular
   *           as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  M      - INTEGER.
   *           On entry, M specifies the number of rows of B. M must be at
   *           least zero.
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the number of columns of B.  N must be
   *           at least zero.
   *           Unchanged on exit.
   *
   *  ALPHA  - DOUBLE PRECISION.
   *           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
   *           zero then  A is not referenced and  B need not be set before
   *           entry.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
   *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
   *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
   *           upper triangular part of the array  A must contain the upper
   *           triangular matrix  and the strictly lower triangular part of
   *           A is not referenced.
   *           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
   *           lower triangular part of the array  A must contain the lower
   *           triangular matrix  and the strictly upper triangular part of
   *           A is not referenced.
   *           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
   *           A  are not referenced either,  but are assumed to be  unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
   *           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
   *           then LDA must be at least max( 1, n ).
   *           Unchanged on exit.
   *
   *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
   *           Before entry,  the leading  m by n part of the array  B must
   *           contain  the  right-hand  side  matrix  B,  and  on exit  is
   *           overwritten by the solution matrix  X.
   *
   *  LDB    - INTEGER.
   *           On entry, LDB specifies the first dimension of B as declared
   *           in  the  calling  (sub)  program.   LDB  must  be  at  least
   *           max( 1, m ).
   *           Unchanged on exit.
   *
   *  Further Details
   *  ===============
   *
   *  Level 3 Blas routine.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  BLASNAME(strsm)( character        const   SIDE[], 
                   character        const   UPLO[], 
                   character        const   TRANSA[], 
                   character        const   DIAG[], 
                   integer          const * M,
                   integer          const * N,
                   single_precision const * alpha, 
                   single_precision const   A[], 
                   integer          const * LDA, 
                   single_precision         B[], 
                   integer          const * LDB ) ;
  void
  BLASNAME(dtrsm)( character        const   SIDE[], 
                   character        const   UPLO[], 
                   character        const   TRANSA[], 
                   character        const   DIAG[], 
                   integer          const * M,
                   integer          const * N, 
                   double_precision const * alpha, 
                   double_precision const   A[], 
                   integer          const * LDA, 
                   double_precision         B[], 
                   integer          const * LDB ) ;
  }
  #endif

  inline
  void
  trsm( SideMultiply  const & SIDE,
        ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               M,
        integer               N,
        real                  alpha,
        real          const   A[],
        integer               LDA,
        real                  B[],
        integer               LDB )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(strsm)( CblasColMajor,
                      side_cblas[SIDE],
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      M, N, alpha, A, LDA, B, LDB ) ; }
  #else
  { BLASNAME(strsm)( side_blas[SIDE],
                     uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &M, &N, &alpha, A, &LDA, B, &LDB ) ; }
  #endif

  inline
  void
  trsm( SideMultiply  const & SIDE,
        ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer               M,
        integer               N,
        doublereal            alpha,
        doublereal    const   A[],
        integer               LDA,
        doublereal            B[],
        integer               LDB )
  #ifdef ALGLIN_USE_CBLAS
  { CBLASNAME(dtrsm)( CblasColMajor,
                      side_cblas[SIDE],
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      M, N, alpha, A, LDA, B, LDB ) ; }
  #else
  { BLASNAME(dtrsm)( side_blas[SIDE],
                     uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &M, &N, &alpha, A, &LDA, B, &LDB ) ; }
  #endif

  /*
  //   _____     _ ____                  _          _
  //  |_   _| __(_) __ )  __ _ _ __   __| | ___  __| |
  //    | || '__| |  _ \ / _` | '_ \ / _` |/ _ \/ _` |
  //    | || |  | | |_) | (_| | | | | (_| |  __/ (_| |
  //    |_||_|  |_|____/ \__,_|_| |_|\__,_|\___|\__,_|
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DTBMV  performs one of the matrix-vector operations
   *
   *     x := A*x,   or   x := A'*x,
   *
   *  where x is an n element vector and  A is an n by n unit, or non-unit,
   *  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
   *
   *  Parameters
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the operation to be performed as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   x := A*x.
   *
   *              TRANS = 'T' or 't'   x := A'*x.
   *
   *              TRANS = 'C' or 'c'   x := A'*x.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  K      - INTEGER.
   *           On entry with UPLO = 'U' or 'u', K specifies the number of
   *           super-diagonals of the matrix A.
   *           On entry with UPLO = 'L' or 'l', K specifies the number of
   *           sub-diagonals of the matrix A.
   *           K must satisfy  0 .le. K.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
   *           by n part of the array A must contain the upper triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row
   *           ( k + 1 ) of the array, the first super-diagonal starting at
   *           position 2 in row k, and so on. The top left k by k triangle
   *           of the array A is not referenced.
   *           The following program segment will transfer an upper
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = K + 1 - J
   *                    DO 10, I = MAX( 1, J - K ), J
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
   *           by n part of the array A must contain the lower triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row 1 of
   *           the array, the first sub-diagonal starting at position 1 in
   *           row 2, and so on. The bottom right k by k triangle of the
   *           array A is not referenced.
   *           The following program segment will transfer a lower
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = 1 - J
   *                    DO 10, I = J, MIN( N, J + K )
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Note that when DIAG = 'U' or 'u' the elements of the array A
   *           corresponding to the diagonal elements of the matrix are not
   *           referenced, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           ( k + 1 ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element vector x. On exit, X is overwritten with the
   *           tranformed vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
  \*/

  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    BLASNAME(stbmv)( char const * UPLO,
                     char const * TRANS,
                     char const * DIAG,
                     integer    * N,
                     integer    * K,
                     real         A[],
                     integer    * LDA,
                     real         X[],
                     integer    * INCXU );
    void
    BLASNAME(dtbmv)( char const * UPLO,
                     char const * TRANS,
                     char const * DIAG,
                     integer    * N,
                     integer    * K,
                     doublereal   A[],
                     integer    * LDA,
                     doublereal   X[],
                     integer    * INCXU );
  }
  #endif

  inline
  void
  tbmv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer             N,
        integer             K,
        real          const A[],
        integer             ldA,
        real                xb[],
        integer             incx ) {
    #ifdef ALGLIN_USE_CBLAS
    CBLASNAME(stbmv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, K, A, ldA, xb, incx ) ;
    #else
    BLASNAME(stbmv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, &K,
                     const_cast<real*>(A), &ldA, xb, &incx ) ;
    #endif
  }

  inline
  void
  tbmv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer             N,
        integer             K,
        doublereal    const A[],
        integer             ldA,
        doublereal          xb[],
        integer             incx ) {
    #ifdef ALGLIN_USE_CBLAS
    CBLASNAME(dtbmv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, K, A, ldA, xb, incx ) ;
    #else
    BLASNAME(dtbmv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, &K,
                     const_cast<doublereal*>(A), &ldA, xb, &incx ) ;
    #endif
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DTBSV  solves one of the systems of equations
   *
   *     A*x = b,   or   A'*x = b,
   *
   *  where b and x are n element vectors and A is an n by n unit, or
   *  non-unit, upper or lower triangular band matrix, with ( k + 1 )
   *  diagonals.
   *
   *  No test for singularity or near-singularity is included in this
   *  routine. Such tests must be performed before calling this routine.
   *
   *  Parameters
   *  ==========
   *
   *  UPLO   - CHARACTER*1.
   *           On entry, UPLO specifies whether the matrix is an upper or
   *           lower triangular matrix as follows:
   *
   *              UPLO = 'U' or 'u'   A is an upper triangular matrix.
   *
   *              UPLO = 'L' or 'l'   A is a lower triangular matrix.
   *
   *           Unchanged on exit.
   *
   *  TRANS  - CHARACTER*1.
   *           On entry, TRANS specifies the equations to be solved as
   *           follows:
   *
   *              TRANS = 'N' or 'n'   A*x = b.
   *
   *              TRANS = 'T' or 't'   A'*x = b.
   *
   *              TRANS = 'C' or 'c'   A'*x = b.
   *
   *           Unchanged on exit.
   *
   *  DIAG   - CHARACTER*1.
   *           On entry, DIAG specifies whether or not A is unit
   *           triangular as follows:
   *
   *              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
   *
   *              DIAG = 'N' or 'n'   A is not assumed to be unit
   *                                  triangular.
   *
   *           Unchanged on exit.
   *
   *  N      - INTEGER.
   *           On entry, N specifies the order of the matrix A.
   *           N must be at least zero.
   *           Unchanged on exit.
   *
   *  K      - INTEGER.
   *           On entry with UPLO = 'U' or 'u', K specifies the number of
   *           super-diagonals of the matrix A.
   *           On entry with UPLO = 'L' or 'l', K specifies the number of
   *           sub-diagonals of the matrix A.
   *           K must satisfy  0 .le. K.
   *           Unchanged on exit.
   *
   *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
   *           by n part of the array A must contain the upper triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row
   *           ( k + 1 ) of the array, the first super-diagonal starting at
   *           position 2 in row k, and so on. The top left k by k triangle
   *           of the array A is not referenced.
   *           The following program segment will transfer an upper
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = K + 1 - J
   *                    DO 10, I = MAX( 1, J - K ), J
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
   *           by n part of the array A must contain the lower triangular
   *           band part of the matrix of coefficients, supplied column by
   *           column, with the leading diagonal of the matrix in row 1 of
   *           the array, the first sub-diagonal starting at position 1 in
   *           row 2, and so on. The bottom right k by k triangle of the
   *           array A is not referenced.
   *           The following program segment will transfer a lower
   *           triangular band matrix from conventional full matrix storage
   *           to band storage:
   *
   *                 DO 20, J = 1, N
   *                    M = 1 - J
   *                    DO 10, I = J, MIN( N, J + K )
   *                       A( M + I, J ) = matrix( I, J )
   *              10    CONTINUE
   *              20 CONTINUE
   *
   *           Note that when DIAG = 'U' or 'u' the elements of the array A
   *           corresponding to the diagonal elements of the matrix are not
   *           referenced, but are assumed to be unity.
   *           Unchanged on exit.
   *
   *  LDA    - INTEGER.
   *           On entry, LDA specifies the first dimension of A as declared
   *           in the calling (sub) program. LDA must be at least
   *           ( k + 1 ).
   *           Unchanged on exit.
   *
   *  X      - DOUBLE PRECISION array of dimension at least
   *           ( 1 + ( n - 1 )*abs( INCX ) ).
   *           Before entry, the incremented array X must contain the n
   *           element right-hand side vector b. On exit, X is overwritten
   *           with the solution vector x.
   *
   *  INCX   - INTEGER.
   *           On entry, INCX specifies the increment for the elements of
   *           X. INCX must not be zero.
   *           Unchanged on exit.
   *
   *
   *  Level 2 Blas routine.
   *
   *  -- Written on 22-October-1986.
   *     Jack Dongarra, Argonne National Lab.
   *     Jeremy Du Croz, Nag Central Office.
   *     Sven Hammarling, Nag Central Office.
   *     Richard Hanson, Sandia National Labs.
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
    void
    BLASNAME(stbsv)( char const * UPLO,
                     char const * TRANS,
                     char const * DIAG,
                     integer    * N,
                     integer    * K,
                     real         A[],
                     integer    * LDA,
                     real         X[],
                     integer    * INCXU );
    void
    BLASNAME(dtbsv)( char const * UPLO,
                     char const * TRANS,
                     char const * DIAG,
                     integer    * N,
                     integer    * K,
                     doublereal   A[],
                     integer    * LDA,
                     doublereal   X[],
                     integer    * INCXU );
  }
  #endif

  inline
  void
  tbsv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer             N,
        integer             K,
        real          const A[],
        integer             ldA,
        real                xb[],
        integer             incx ) {
    #ifdef ALGLIN_USE_CBLAS
    CBLASNAME(stbsv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, K, A, ldA, xb, incx ) ;
    #else
    BLASNAME(stbsv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, &K,
                     const_cast<real*>(A), &ldA, xb, &incx ) ;
    #endif
  }

  inline
  void
  tbsv( ULselect      const & UPLO,
        Transposition const & TRANS,
        DiagonalType  const & DIAG,
        integer             N,
        integer             K,
        doublereal    const A[],
        integer             ldA,
        doublereal          xb[],
        integer             incx ) {
    #ifdef ALGLIN_USE_CBLAS
    CBLASNAME(dtbsv)( CblasColMajor,
                      uplo_cblas[UPLO],
                      trans_cblas[TRANS],
                      diag_cblas[DIAG],
                      N, K, A, ldA, xb, incx ) ;
    #else
    BLASNAME(dtbsv)( uplo_blas[UPLO],
                     trans_blas[TRANS],
                     diag_blas[DIAG],
                     &N, &K,
                     const_cast<doublereal*>(A), &ldA, xb, &incx ) ;
    #endif
  }

  /*
  //   _____     _     _ _                               _
  //  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |
  //    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
  //    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |
  //    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
  //                             |___/
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DGTTRF computes an LU factorization of a real tridiagonal matrix A
   *  using elimination with partial pivoting and row interchanges.
   *
   *  The factorization has the form
   *     A = L * U
   *  where L is a product of permutation and unit lower bidiagonal
   *  matrices and U is upper triangular with nonzeros in only the main
   *  diagonal and first two superdiagonals.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.
   *
   *  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DL must contain the (n-1) sub-diagonal elements of
   *          A.
   *
   *          On exit, DL is overwritten by the (n-1) multipliers that
   *          define the matrix L from the LU factorization of A.
   *
   *  D       (input/output) DOUBLE PRECISION array, dimension (N)
   *          On entry, D must contain the diagonal elements of A.
   *
   *          On exit, D is overwritten by the n diagonal elements of the
   *          upper triangular matrix U from the LU factorization of A.
   *
   *  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DU must contain the (n-1) super-diagonal elements
   *          of A.
   *
   *          On exit, DU is overwritten by the (n-1) elements of the first
   *          super-diagonal of U.
   *
   *  DU2     (output) DOUBLE PRECISION array, dimension (N-2)
   *          On exit, DU2 is overwritten by the (n-2) elements of the
   *          second super-diagonal of U.
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -k, the k-th argument had an illegal value
   *          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
   *                has been completed, but the factor U is exactly
   *                singular, and division by zero will occur if it is used
   *                to solve a system of equations.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgttrf)( integer * N,
                        real    * DL,
                        real    * D,
                        real    * DU,
                        real    * DU2,
                        integer * IPIV,
                        integer * INFO ) ;

    void
    LAPACKNAME(dgttrf)( integer    * N,
                        doublereal * DL,
                        doublereal * D,
                        doublereal * DU,
                        doublereal * DU2,
                        integer    * IPIV,
                        integer    * INFO ) ;
  }
  #endif

  inline
  integer
  gttrf( integer N,
         real    DL[],
         real    D[],
         real    DU[],
         real    DU2[],
         integer IPIV[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer INFO ; LAPACK_NAME(sgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO ; CLAPACKNAME(sgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #else
  { integer INFO ; LAPACKNAME(sgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #endif

  inline
  integer
  gttrf( integer    N,
         doublereal DL[],
         doublereal D[],
         doublereal DU[],
         doublereal DU2[],
         integer    IPIV[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer INFO ; LAPACK_NAME(dgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO ; CLAPACKNAME(dgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #else
  { integer INFO ; LAPACKNAME(dgttrf)( &N, DL, D, DU, DU2, IPIV, &INFO ) ; return INFO ; }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGTTRS solves one of the systems of equations
   *     A*X = B  or  A**T*X = B,
   *  with a tridiagonal matrix A using the LU factorization computed
   *  by DGTTRF.
   *
   *  Arguments
   *  =========
   *
   *  TRANS   (input) CHARACTER*1
   *          Specifies the form of the system of equations.
   *          = 'N':  A * X = B  (No transpose)
   *          = 'T':  A**T* X = B  (Transpose)
   *          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  DL      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) multipliers that define the matrix L from the
   *          LU factorization of A.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the upper triangular matrix U from
   *          the LU factorization of A.
   *
   *  DU      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) elements of the first super-diagonal of U.
   *
   *  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
   *          The (n-2) elements of the second super-diagonal of U.
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the matrix of right hand side vectors B.
   *          On exit, B is overwritten by the solution vectors X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgttrs)( character * TRANS,
                        integer   * N,
                        integer   * NRHS,
                        real      * DL,
                        real      * D,
                        real      * DU,
                        real      * DU2,
                        integer   * IPIV,
                        real      * B,
                        integer   * LDB,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgttrs)( character  * TRANS,
                        integer    * N,
                        integer    * NRHS,
                        doublereal * DL,
                        doublereal * D,
                        doublereal * DU,
                        doublereal * DU2,
                        integer    * IPIV,
                        doublereal * B,
                        integer    * LDB,
                        integer    * INFO ) ;
  }
  #endif

  inline
  integer
  gttrs( Transposition const & TRANS,
         integer               N,
         integer               NRHS,
         real          const   DL[],
         real          const   D[],
         real          const   DU[],
         real          const   DU2[],
         integer       const   IPIV[],
         real                  B[],
         integer               LDB )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgttrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgttrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS,
                         const_cast<real*>(DL),
                         const_cast<real*>(D),
                         const_cast<real*>(DU),
                         const_cast<real*>(DU2),
                         const_cast<integer*>(IPIV), B, &LDB, &INFO ) ;
    #else
    LAPACKNAME(sgttrs)( const_cast<character*>(trans_blas[TRANS]),
                        &N, &NRHS,
                        const_cast<real*>(DL),
                        const_cast<real*>(D),
                        const_cast<real*>(DU),
                        const_cast<real*>(DU2),
                        const_cast<integer*>(IPIV),
                        B, &LDB, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gttrs( Transposition const & TRANS,
         integer               N,
         integer               NRHS,
         doublereal    const   DL[],
         doublereal    const   D[],
         doublereal    const   DU[],
         doublereal    const   DU2[],
         integer       const   IPIV[],
         doublereal            B[],
         integer               LDB )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgttrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgttrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS,
                         const_cast<doublereal*>(DL),
                         const_cast<doublereal*>(D),
                         const_cast<doublereal*>(DU),
                         const_cast<doublereal*>(DU2),
                         const_cast<integer*>(IPIV),
                         B, &LDB, &INFO ) ;
    #else
    LAPACKNAME(dgttrs)( const_cast<character*>(trans_blas[TRANS]),
                        &N, &NRHS,
                        const_cast<doublereal*>(DL),
                        const_cast<doublereal*>(D),
                        const_cast<doublereal*>(DU),
                        const_cast<doublereal*>(DU2),
                        const_cast<integer*>(IPIV),
                        B, &LDB, &INFO ) ;
    #endif
    return INFO ;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGTSV  solves the equation
   *
   *     A*X = B,
   *
   *  where A is an n by n tridiagonal matrix, by Gaussian elimination with
   *  partial pivoting.
   *
   *  Note that the equation  A**T*X = B  may be solved by interchanging the
   *  order of the arguments DU and DL.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DL must contain the (n-1) sub-diagonal elements of
   *          A.
   *
   *          On exit, DL is overwritten by the (n-2) elements of the
   *          second super-diagonal of the upper triangular matrix U from
   *          the LU factorization of A, in DL(1), ..., DL(n-2).
   *
   *  D       (input/output) DOUBLE PRECISION array, dimension (N)
   *          On entry, D must contain the diagonal elements of A.
   *
   *          On exit, D is overwritten by the n diagonal elements of U.
   *
   *  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
   *          On entry, DU must contain the (n-1) super-diagonal elements
   *          of A.
   *
   *          On exit, DU is overwritten by the (n-1) elements of the first
   *          super-diagonal of U.
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N by NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N by NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
   *               has not been computed.  The factorization has not been
   *               completed unless i = N.
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgtsv)( integer * N,
                       integer * NRHS,
                       real    * DL,
                       real    * D,
                       real    * DU,
                       real    * B,
                       integer * LDB,
                       integer * INFO ) ;

    void
    LAPACKNAME(dgtsv)( integer    * N,
                       integer    * NRHS,
                       doublereal * DL,
                       doublereal * D,
                       doublereal * DU,
                       doublereal * B,
                       integer    * LDB,
                       integer    * INFO ) ;
  }
  #endif

  inline
  integer
  gtsv( integer N,
        integer NRHS,
        real    DL[],
        real    D[],
        real    DU[],
        real    B[],
        integer LDB )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #else
    LAPACKNAME(sgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gtsv( integer    N,
        integer    NRHS,
        doublereal DL[],
        doublereal D[],
        doublereal DU[],
        doublereal B[],
        integer    LDB )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #else
    LAPACKNAME(dgtsv)( &N, &NRHS, DL, D, DU, B, &LDB, &INFO ) ;
    #endif
    return INFO ;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGTCON estimates the reciprocal of the condition number of a real
   *  tridiagonal matrix A using the LU factorization as computed by
   *  DGTTRF.
   *
   *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
   *  condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
   *
   *  Arguments
   *  =========
   *
   *  NORM    (input) CHARACTER*1
   *          Specifies whether the 1-norm condition number or the
   *          infinity-norm condition number is required:
   *          = '1' or 'O':  1-norm;
   *          = 'I':         Infinity-norm.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  DL      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) multipliers that define the matrix L from the
   *          LU factorization of A as computed by DGTTRF.
   *
   *  D       (input) DOUBLE PRECISION array, dimension (N)
   *          The n diagonal elements of the upper triangular matrix U from
   *          the LU factorization of A.
   *
   *  DU      (input) DOUBLE PRECISION array, dimension (N-1)
   *          The (n-1) elements of the first superdiagonal of U.
   *
   *  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
   *          The (n-2) elements of the second superdiagonal of U.
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices; for 1 <= i <= n, row i of the matrix was
   *          interchanged with row IPIV(i).  IPIV(i) will always be either
   *          i or i+1; IPIV(i) = i indicates a row interchange was not
   *          required.
   *
   *  ANORM   (input) DOUBLE PRECISION
   *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
   *          If NORM = 'I', the infinity-norm of the original matrix A.
   *
   *  RCOND   (output) DOUBLE PRECISION
   *          The reciprocal of the condition number of the matrix A,
   *          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
   *          estimate of the 1-norm of inv(A) computed in this routine.
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
   *
   *  IWORK   (workspace) INTEGER array, dimension (N)
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgtcon)( character * NORM,
                        integer   * N,
                        real      * DL,
                        real      * D,
                        real      * DU,
                        real      * DU2,
                        integer   * IPIV,
                        real      * ANORM,
                        real      * RCOND,
                        real      * WORK,
                        integer   * IWORK,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgtcon)( character  * NORM,
                        integer    * N,
                        doublereal * DL,
                        doublereal * D,
                        doublereal * DU,
                        doublereal * DU2,
                        integer    * IPIV,
                        doublereal * ANORM,
                        doublereal * RCOND,
                        doublereal * WORK,
                        integer    * IWORK,
                        integer    * INFO ) ;
  }
  #endif

  inline
  integer
  gtcon1( integer       N,
          real    const DL[],
          real    const D[],
          real    const DU[],
          real    const DU2[],
          integer const IPIV[],
          real          ANORM,
          real        & RCOND,
          real          WORK[],
          integer       IWORK[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgtcon)( const_cast<character*>("1"), &N,
                         const_cast<real*>(DL),
                         const_cast<real*>(D),
                         const_cast<real*>(DU),
                         const_cast<real*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgtcon)( const_cast<character*>("1"), &N,
                         const_cast<real*>(DL),
                         const_cast<real*>(D),
                         const_cast<real*>(DU),
                         const_cast<real*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #else
    LAPACKNAME(sgtcon)( const_cast<character*>("1"), &N,
                        const_cast<real*>(DL),
                        const_cast<real*>(D),
                        const_cast<real*>(DU),
                        const_cast<real*>(DU2),
                        const_cast<integer*>(IPIV),
                        &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gtcon1( integer          N,
          doublereal const DL[],
          doublereal const D[],
          doublereal const DU[],
          doublereal const DU2[],
          integer    const IPIV[],
          doublereal       ANORM,
          doublereal     & RCOND,
          doublereal       WORK[],
          integer          IWORK[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgtcon)( const_cast<character*>("1"), &N,
                         const_cast<doublereal*>(DL),
                         const_cast<doublereal*>(D),
                         const_cast<doublereal*>(DU),
                         const_cast<doublereal*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgtcon)( const_cast<character*>("1"), &N,
                         const_cast<doublereal*>(DL),
                         const_cast<doublereal*>(D),
                         const_cast<doublereal*>(DU),
                         const_cast<doublereal*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #else
    LAPACKNAME(dgtcon)( const_cast<character*>("1"), &N,
                        const_cast<doublereal*>(DL),
                        const_cast<doublereal*>(D),
                        const_cast<doublereal*>(DU),
                        const_cast<doublereal*>(DU2),
                        const_cast<integer*>(IPIV),
                        &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gtconInf( integer       N,
            real    const DL[],
            real    const D[],
            real    const DU[],
            real    const DU2[],
            integer const IPIV[],
            real          ANORM,
            real        & RCOND,
            real          WORK[],
            integer       IWORK[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgtcon)( const_cast<character*>("I"), &N,
                         const_cast<real*>(DL),
                         const_cast<real*>(D),
                         const_cast<real*>(DU),
                         const_cast<real*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgtcon)( const_cast<character*>("I"), &N,
                         const_cast<real*>(DL),
                         const_cast<real*>(D),
                         const_cast<real*>(DU),
                         const_cast<real*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #else
    LAPACKNAME(sgtcon)( const_cast<character*>("I"), &N,
                        const_cast<real*>(DL),
                        const_cast<real*>(D),
                        const_cast<real*>(DU),
                        const_cast<real*>(DU2),
                        const_cast<integer*>(IPIV),
                        &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gtconInf( integer          N,
            doublereal const DL[],
            doublereal const D[],
            doublereal const DU[],
            doublereal const DU2[],
            integer    const IPIV[],
            doublereal       ANORM,
            doublereal     & RCOND,
            doublereal       WORK[],
            integer          IWORK[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgtcon)( const_cast<character*>("I"), &N,
                         const_cast<doublereal*>(DL),
                         const_cast<doublereal*>(D),
                         const_cast<doublereal*>(DU),
                         const_cast<doublereal*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgtcon)( const_cast<character*>("I"), &N,
                         const_cast<doublereal*>(DL),
                         const_cast<doublereal*>(D),
                         const_cast<doublereal*>(DU),
                         const_cast<doublereal*>(DU2),
                         const_cast<integer*>(IPIV),
                         &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #else
    LAPACKNAME(dgtcon)( const_cast<character*>("I"), &N,
                        const_cast<doublereal*>(DL),
                        const_cast<doublereal*>(D),
                        const_cast<doublereal*>(DU),
                        const_cast<doublereal*>(DU2),
                        const_cast<integer*>(IPIV),
                        &ANORM, &RCOND, WORK, IWORK, &INFO ) ;
    #endif
    return INFO ;
  }

  /*
  //              _         __
  //    __ _  ___| |_ _ __ / _|
  //   / _` |/ _ \ __| '__| |_
  //  | (_| |  __/ |_| |  |  _|
  //   \__, |\___|\__|_|  |_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETRF computes an LU factorization of a general M-by-N matrix A
   *  using partial pivoting with row interchanges.
   *
   *  The factorization has the form
   *     A = P * L * U
   *  where P is a permutation matrix, L is lower triangular with unit
   *  diagonal elements (lower trapezoidal if m > n), and U is upper
   *  triangular (upper trapezoidal if m < n).
   *
   *  This is the right-looking Level 3 BLAS version of the algorithm.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix to be factored.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  IPIV    (output) INTEGER array, dimension (min(M,N))
   *          The pivot indices; for 1 <= i <= min(M,N), row i of the
   *          matrix was interchanged with row IPIV(i).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
   *                has been completed, but the factor U is exactly
   *                singular, and division by zero will occur if it is used
   *                to solve a system of equations.
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {  
  using namespace blas_type ;
  void
  LAPACKNAME(sgetrf)( integer          const * N, 
                      integer          const * M, 
                      single_precision         A[], 
                      integer          const * LDA, 
                      integer                  IPIV[], 
                      integer                * INFO ) ;
  void
  LAPACKNAME(dgetrf)( integer          const * N, 
                      integer          const * M, 
                      double_precision         A[], 
                      integer          const * LDA, 
                      integer                  IPIV[], 
                      integer                * INFO ) ;
  }
  #endif

  inline
  integer
  getrf( integer N,
         integer M,
         real    A[],
         integer LDA,
         integer IPIV[] )
  #if defined(ALGLIN_USE_ATLAS)
  { return CLAPACKNAME(sgetrf)( CblasColMajor, N, M, A, LDA, IPIV ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0; CLAPACKNAME(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ; return INFO ; }
  #else
  { integer INFO = 0; LAPACKNAME(sgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ; return INFO ; }
  #endif

  inline
  integer
  getrf( integer    N,
         integer    M,
         doublereal A[],
         integer    LDA,
         integer    IPIV[] )
  #if defined(ALGLIN_USE_ATLAS)
  { return CLAPACKNAME(dgetrf)( CblasColMajor, N, M, A, LDA, IPIV ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO  = 0; CLAPACKNAME(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ; return INFO ; }
  #else
  { integer INFO = 0; LAPACKNAME(dgetrf)( &N, &M, A, &LDA, IPIV, &INFO ) ; return INFO ; }
  #endif

  /*
  //              _
  //    __ _  ___| |_ _ __ ___
  //   / _` |/ _ \ __| '__/ __|
  //  | (_| |  __/ |_| |  \__ \
  //   \__, |\___|\__|_|  |___/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETRS solves a system of linear equations
   *     A * X = B  or  A**T * X = B
   *  with a general N-by-N matrix A using the LU factorization computed
   *  by DGETRF.
   *
   *  Arguments
   *  =========
   *
   *  TRANS   (input) CHARACTER*1
   *          Specifies the form of the system of equations:
   *          = 'N':  A * X = B  (No transpose)
   *          = 'T':  A**T* X = B  (Transpose)
   *          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The factors L and U from the factorization A = P*L*U
   *          as computed by DGETRF.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (input) INTEGER array, dimension (N)
   *          The pivot indices from DGETRF; for 1<=i<=N, row i of the
   *          matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the right hand side matrix B.
   *          On exit, the solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  LAPACKNAME(sgetrs)( character        const   TRANS[], 
                      integer          const * N, 
                      integer          const * NRHS, 
                      single_precision const   A[], 
                      integer          const * LDA, 
                      integer          const   IPIV[], 
                      single_precision         B[], 
                      integer          const * LDB, 
                      integer                * INFO ) ;
  void
  LAPACKNAME(dgetrs)( character        const   TRANS[],
                      integer          const * N, 
                      integer          const * NRHS, 
                      double_precision const   A[], 
                      integer          const * LDA, 
                      integer          const   IPIV[], 
                      double_precision         B[], 
                      integer          const * LDB, 
                      integer                * INFO ) ;
  }
  #endif

  inline
  integer
  getrs( Transposition const & TRANS,
         integer               N,
         integer               NRHS,
         real          const   A[],
         integer               LDA,
         integer       const   IPIV[],
         real                  B[],
         integer               LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { CLAPACKNAME(sgetrs)( CblasColMajor, trans_cblas[TRANS], N, NRHS, A, LDA, IPIV, B, LDB ) ; return 0 ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(sgetrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(sgetrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS,
                         const_cast<real*>(A), &LDA,
                         const_cast<integer*>(IPIV),
                         B, &LDB, &INFO ) ;
    return INFO ; }
  #else
  { integer INFO = 0;
    LAPACKNAME(sgetrs)( trans_blas[TRANS], &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #endif

  inline
  integer
  getrs( Transposition const & TRANS,
         integer               N,
         integer               NRHS,
         doublereal    const   A[],
         integer               LDA,
         integer       const   IPIV[],
         doublereal            B[],
         integer               LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { CLAPACKNAME(dgetrs)( CblasColMajor, trans_cblas[TRANS], N, NRHS, A, LDA, IPIV, B, LDB ) ; return 0 ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(dgetrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(dgetrs)( const_cast<character*>(trans_blas[TRANS]),
                         &N, &NRHS,
                         const_cast<doublereal*>(A), &LDA,
                         const_cast<integer*>(IPIV),
                         B, &LDB, &INFO ) ;
    return INFO ; }
  #else
  { integer INFO = 0;
    LAPACKNAME(dgetrs)( trans_blas[TRANS], &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #endif

  /*
  //    __ _  ___  _____   __
  //   / _` |/ _ \/ __\ \ / /
  //  | (_| |  __/\__ \\ V /
  //   \__, |\___||___/ \_/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESV computes the solution to a real system of linear equations
   *     A * X = B,
   *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   *
   *  The LU decomposition with partial pivoting and row interchanges is
   *  used to factor A as
   *     A = P * L * U,
   *  where P is a permutation matrix, L is unit lower triangular, and U is
   *  upper triangular.  The factored form of A is then used to solve the
   *  system of equations A * X = B.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void 
  LAPACKNAME(sgesv)( integer          const * N, 
                     integer          const * NRHS, 
                     single_precision         A[], 
                     integer          const * LDA, 
                     integer                  IPIV[], 
                     single_precision         B[], 
                     integer          const * LDB, 
                     integer                * INFO ) ;
  void
  LAPACKNAME(dgesv)( integer          const * N, 
                     integer          const * NRHS, 
                     double_precision         A[], 
                     integer          const * LDA, 
                     integer                  IPIV[], 
                     double_precision         B[], 
                     integer          const * LDB, 
                     integer                * INFO ) ;
  }
  #endif

  inline
  integer
  gesv( integer N,
        integer NRHS,
        real    A[],
        integer LDA,
        integer IPIV[],
        real    B[],
        integer LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { integer INFO = CLAPACKNAME(sgetrf)( CblasColMajor, N, N, A, LDA, IPIV ) ;
    if ( INFO == 0 ) CLAPACKNAME(sgetrs)( CblasColMajor, CblasNoTrans, N, NRHS, A, LDA, IPIV, B, LDB ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0; CLAPACKNAME(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ; return INFO ; }
  #else
  { integer INFO = 0; LAPACKNAME(sgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ; return INFO ; }
  #endif

  inline
  integer
  gesv( integer    N,
        integer    NRHS,
        doublereal A[],
        integer    LDA,
        integer    IPIV[],
        doublereal B[],
        integer    LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { integer INFO = CLAPACKNAME(dgetrf)( CblasColMajor, N, N, A, LDA, IPIV );
    if ( INFO == 0 ) CLAPACKNAME(dgetrs)( CblasColMajor, CblasNoTrans, N, NRHS, A, LDA, IPIV, B, LDB ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { integer INFO = 0;
    LAPACK_NAME(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ;
    return INFO ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0; CLAPACKNAME(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ; return INFO ; }
  #else
  { integer INFO = 0; LAPACKNAME(dgesv)( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO ) ; return INFO ; }
  #endif

  /*
  //    ____ _____ _____ ____ ____
  //   / ___| ____|_   _/ ___|___ \
  //  | |  _|  _|   | || |     __) |
  //  | |_| | |___  | || |___ / __/
  //   \____|_____| |_| \____|_____|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGETC2 computes an LU factorization with complete pivoting of the
   *  n-by-n matrix A. The factorization has the form A = P * L * U * Q,
   *  where P and Q are permutation matrices, L is lower triangular with
   *  unit diagonal elements and U is upper triangular.
   *
   *  This is the Level 2 BLAS algorithm.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A. N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the n-by-n matrix A to be factored.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U*Q; the unit diagonal elements of L are not stored.
   *          If U(k, k) appears to be less than SMIN, U(k, k) is given the
   *          value of SMIN, i.e., giving a nonsingular perturbed system.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension(N).
   *          The pivot indices; for 1 <= i <= N, row i of the
   *          matrix has been interchanged with row IPIV(i).
   *
   *  JPIV    (output) INTEGER array, dimension(N).
   *          The pivot indices; for 1 <= j <= N, column j of the
   *          matrix has been interchanged with column JPIV(j).
   *
   *  INFO    (output) INTEGER
   *           = 0: successful exit
   *           > 0: if INFO = k, U(k, k) is likely to produce owerflow if
   *                we try to solve for x in Ax = b. So U is perturbed to
   *                avoid the overflow.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
   *     Umea University, S-901 87 Umea, Sweden.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgetc2)( integer   * N,
                        real      * A,
                        integer   * LDA,
                        integer   * IPIV,
                        integer   * JPIV,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgetc2)( integer    * N,
                        doublereal * A,
                        integer    * LDA,
                        integer    * IPIV,
                        integer    * JPIV,
                        integer    * INFO ) ;
  }
  #endif

  #ifdef ALGLIN_USE_OPENBLAS
  template <typename T>
  integer
  getc2_tmpl( integer N,
              T       A[],
              integer LDA,
              integer IPIV[],
              integer JPIV[] ) ;
  #endif

  inline
  integer
  getc2( integer   N,
         real      A[],
         integer   LDA,
         integer   IPIV[],
         integer   JPIV[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    INFO = getc2_tmpl<real>(N,A,LDA,IPIV,JPIV) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO ) ;
    #else
    LAPACKNAME(sgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  getc2( integer    N,
         doublereal A[],
         integer    LDA,
         integer    IPIV[],
         integer    JPIV[] )
  { integer INFO = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    INFO = getc2_tmpl<doublereal>(N,A,LDA,IPIV,JPIV) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO ) ;
    #else
    LAPACKNAME(dgetc2)( &N, A, &LDA, IPIV, JPIV, &INFO ) ;
    #endif
    return INFO ;
  }

  /*
  //    ____ _____ ____   ____ ____
  //   / ___| ____/ ___| / ___|___ \
  //  | |  _|  _| \___ \| |     __) |
  //  | |_| | |___ ___) | |___ / __/
  //   \____|_____|____/ \____|_____|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESC2 solves a system of linear equations
   *
   *            A * X = scale* RHS
   *
   *  with a general N-by-N matrix A using the LU factorization with
   *  complete pivoting computed by DGETC2.
   *
   *  Arguments
   *  =========
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the  LU part of the factorization of the n-by-n
   *          matrix A computed by DGETC2:  A = P * L * U * Q
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1, N).
   *
   *  RHS     (input/output) DOUBLE PRECISION array, dimension (N).
   *          On entry, the right hand side vector b.
   *          On exit, the solution vector X.
   *
   *  IPIV    (input) INTEGER array, dimension (N).
   *          The pivot indices; for 1 <= i <= N, row i of the
   *          matrix has been interchanged with row IPIV(i).
   *
   *  JPIV    (input) INTEGER array, dimension (N).
   *          The pivot indices; for 1 <= j <= N, column j of the
   *          matrix has been interchanged with column JPIV(j).
   *
   *  SCALE    (output) DOUBLE PRECISION
   *           On exit, SCALE contains the scale factor. SCALE is chosen
   *           0 <= SCALE <= 1 to prevent owerflow in the solution.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
   *     Umea University, S-901 87 Umea, Sweden.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgesc2)( integer   * N,
                        real      * A,
                        integer   * LDA,
                        real      * RHS,
                        integer   * IPIV,
                        integer   * JPIV,
                        real      * SCALE ) ;

    void
    LAPACKNAME(dgesc2)( integer    * N,
                        doublereal * A,
                        integer    * LDA,
                        doublereal * RHS,
                        integer    * IPIV,
                        integer    * JPIV,
                        doublereal * SCALE ) ;
  }
  #endif

  #ifdef ALGLIN_USE_OPENBLAS
  template <typename T>
  T
  gesc2_tmpl( integer       N,
              T       const A[],
              integer       LDA,
              T             RHS[],
              integer const IPIV[],
              integer const JPIV[] ) ;
  #endif

  inline
  real
  gesc2( integer       N,
         real    const A[],
         integer       LDA,
         real          RHS[],
         integer const IPIV[],
         integer const JPIV[] )
  { real SCALE = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    SCALE = gesc2_tmpl<real>( N, A, LDA, RHS, IPIV, JPIV ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgesc2)( &N,
                         const_cast<real*>(A), &LDA, RHS,
                         const_cast<integer*>(IPIV),
                         const_cast<integer*>(JPIV), &SCALE ) ;
    #else
    LAPACKNAME(sgesc2)( &N,
                        const_cast<real*>(A), &LDA, RHS,
                        const_cast<integer*>(IPIV),
                        const_cast<integer*>(JPIV), &SCALE ) ;
    #endif
    return SCALE ;
  }

  inline
  doublereal
  gesc2( integer          N,
         doublereal const A[],
         integer          LDA,
         doublereal       RHS[],
         integer    const IPIV[],
         integer    const JPIV[] )
  { doublereal SCALE = 0;
    #if defined(ALGLIN_USE_OPENBLAS)
    SCALE = gesc2_tmpl<doublereal>( N, A, LDA, RHS, IPIV, JPIV ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgesc2)( &N,
                         const_cast<doublereal*>(A), &LDA, RHS,
                         const_cast<integer*>(IPIV),
                         const_cast<integer*>(JPIV),
                         &SCALE ) ;
    #else
    LAPACKNAME(dgesc2)( &N,
                        const_cast<doublereal*>(A), &LDA, RHS,
                        const_cast<integer*>(IPIV),
                        const_cast<integer*>(JPIV),
                        &SCALE ) ;
    #endif
    return SCALE ;
  }

  /*
  //    ____ _____ ____ ___  _   _
  //   / ___| ____/ ___/ _ \| \ | |
  //  | |  _|  _|| |  | | | |  \| |
  //  | |_| | |__| |__| |_| | |\  |
  //   \____|_____\____\___/|_| \_|
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGECON estimates the reciprocal of the condition number of a general
   *  real matrix A, in either the 1-norm or the infinity-norm, using
   *  the LU factorization computed by DGETRF.
   *
   *  An estimate is obtained for norm(inv(A)), and the reciprocal of the
   *  condition number is computed as
   *     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
   *
   *  Arguments
   *  =========
   *
   *  NORM    (input) CHARACTER*1
   *          Specifies whether the 1-norm condition number or the
   *          infinity-norm condition number is required:
   *          = '1' or 'O':  1-norm;
   *          = 'I':         Infinity-norm.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The factors L and U from the factorization A = P*L*U
   *          as computed by DGETRF.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  ANORM   (input) DOUBLE PRECISION
   *          If NORM = '1' or 'O', the 1-norm of the original matrix A.
   *          If NORM = 'I', the infinity-norm of the original matrix A.
   *
   *  RCOND   (output) DOUBLE PRECISION
   *          The reciprocal of the condition number of the matrix A,
   *          computed as RCOND = 1/(norm(A) * norm(inv(A))).
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
   *
   *  IWORK   (workspace) INTEGER array, dimension (N)
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgecon)( character * NORM,
                        integer   * N,
                        real      * A,
                        integer   * ldA,
                        real      * ANORM,
                        real      * RCOND,
                        real      * WORK,
                        integer   * IWORK,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgecon)( character  * NORM,
                        integer    * N,
                        doublereal * A,
                        integer    * ldA,
                        doublereal * ANORM,
                        doublereal * RCOND,
                        doublereal * WORK,
                        integer    * IWORK,
                        integer    * INFO ) ;
  }
  #endif

  inline
  integer
  gecon1( integer    N,
          real const A[],
          integer    LDA,
          real       anorm,
          real     & rcond,
          real       work[],
          integer    iwork[] ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(sgecon)( const_cast<character*>("1"), &N,
                        const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgecon)( const_cast<character*>("1"), &N,
                         const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgecon)( const_cast<character*>("1"), &N,
                         const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #else
    LAPACKNAME(sgecon)( const_cast<character*>("1"), &N,
                        const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO) ;
    #endif
    return INFO ;
  }

  inline
  integer
  gecon1( integer          N,
          doublereal const A[],
          integer          LDA,
          doublereal       anorm,
          doublereal     & rcond,
          doublereal       work[],
          integer          iwork[] ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(dgecon)( const_cast<character*>("1"), &N,
                        const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgecon)( const_cast<character*>("1"), &N,
                         const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgecon)( const_cast<character*>("1"), &N,
                         const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #else
    LAPACKNAME(dgecon)( const_cast<character*>("1"), &N,
                        const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO) ;
    #endif
    return INFO ;
  }

  inline
  integer
  geconInf( integer    N,
            real const A[],
            integer    LDA,
            real       anorm,
            real     & rcond,
            real       work[],
            integer    iwork[] ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(sgecon)( const_cast<character*>("I"), &N,
                        const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgecon)( const_cast<character*>("I"), &N,
                         const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgecon)( const_cast<character*>("I"), &N,
                         const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #else
    LAPACKNAME(sgecon)( const_cast<character*>("I"), &N,
                        const_cast<real*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO) ;
    #endif
    return INFO ;
  }

  inline
  integer
  geconInf( integer          N,
            doublereal const A[],
            integer          LDA,
            doublereal       anorm,
            doublereal     & rcond,
            doublereal       work[],
            integer          iwork[] ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(dgecon)( const_cast<character*>("I"), &N,
                        const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgecon)( const_cast<character*>("I"), &N,
                         const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgecon)( const_cast<character*>("I"), &N,
                         const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO ) ;
    #else
    LAPACKNAME(dgecon)( const_cast<character*>("I"), &N,
                        const_cast<doublereal*>(A), &LDA, &anorm, &rcond, work, iwork, &INFO) ;
    #endif
    return INFO ;
  }

  /*\
   *  Purpose
   *  =======
   *
   *  DGEEQU computes row and column scalings intended to equilibrate an
   *  M-by-N matrix A and reduce its condition number.  R returns the row
   *  scale factors and C the column scale factors, chosen to try to make
   *  the largest element in each row and column of the matrix B with
   *  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
   *
   *  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
   *  number and BIGNUM = largest safe number.  Use of these scaling
   *  factors is not guaranteed to reduce the condition number of A but
   *  works well in practice.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The M-by-N matrix whose equilibration factors are
   *          to be computed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  R       (output) DOUBLE PRECISION array, dimension (M)
   *          If INFO = 0 or INFO > M, R contains the row scale factors
   *          for A.
   *
   *  C       (output) DOUBLE PRECISION array, dimension (N)
   *          If INFO = 0,  C contains the column scale factors for A.
   *
   *  ROWCND  (output) DOUBLE PRECISION
   *          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
   *          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
   *          AMAX is neither too large nor too small, it is not worth
   *          scaling by R.
   *
   *  COLCND  (output) DOUBLE PRECISION
   *          If INFO = 0, COLCND contains the ratio of the smallest
   *          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
   *          worth scaling by C.
   *
   *  AMAX    (output) DOUBLE PRECISION
   *          Absolute value of largest matrix element.  If AMAX is very
   *          close to overflow or very close to underflow, the matrix
   *          should be scaled.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i,  and i is
   *                <= M:  the i-th row of A is exactly zero
   *                >  M:  the (i-M)-th column of A is exactly zero
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgeequ)( integer   * M,
                        integer   * N,
                        real      * A,
                        integer   * ldA,
                        real      * R,
                        real      * C,
                        real      * ROWCND,
                        real      * COLCND,
                        real      * AMAX,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgeequ)( integer    * M,
                        integer    * N,
                        doublereal * A,
                        integer    * ldA,
                        doublereal * R,
                        doublereal * C,
                        doublereal * ROWCND,
                        doublereal * COLCND,
                        doublereal * AMAX,
                        integer    * INFO ) ;
  }
  #endif

  inline
  integer
  geequ( integer    M,
         integer    N,
         real const A[],
         integer    LDA,
         real       R[],
         real       C[],
         real     & ROWCND,
         real     & COLCND,
         real     & AMAX ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(sgeequ)( &M, &N, const_cast<real*>(A), &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgeequ)( &M, &N, const_cast<real*>(A), &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgeequ)( &M, &N, const_cast<real*>(A), &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #else
    LAPACKNAME(sgeequ)( &M, &N, const_cast<real*>(A), &LDA, R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #endif
    return INFO ;
  }

  inline
  integer
  geequ( integer          M,
         integer          N,
         doublereal const A[],
         integer          LDA,
         doublereal       R[],
         doublereal       C[],
         doublereal     & ROWCND,
         doublereal     & COLCND,
         doublereal     & AMAX ) {
    integer INFO = 0 ;
    #if defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(dgeequ)( &M, &N,
                        const_cast<doublereal*>(A), &LDA,
                        R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgeequ)( &M, &N,
                         const_cast<doublereal*>(A), &LDA,
                         R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgeequ)( &M, &N,
                        const_cast<doublereal*>(A), &LDA,
                        R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #else
    LAPACKNAME(dgeequ)( &M, &N,
                        const_cast<doublereal*>(A), &LDA,
                        R, C, &ROWCND, &COLCND, &AMAX, &INFO ) ;
    #endif
    return INFO ;
  }

  /*\
   * SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED )
   *
   *  Purpose
   *  =======
   *
   *  DLAQGE equilibrates a general M by N matrix A using the row and
   *  column scaling factors in the vectors R and C.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M by N matrix A.
   *          On exit, the equilibrated matrix.  See EQUED for the form of
   *          the equilibrated matrix.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(M,1).
   *
   *  R       (input) DOUBLE PRECISION array, dimension (M)
   *          The row scale factors for A.
   *
   *  C       (input) DOUBLE PRECISION array, dimension (N)
   *          The column scale factors for A.
   *
   *  ROWCND  (input) DOUBLE PRECISION
   *          Ratio of the smallest R(i) to the largest R(i).
   *
   *  COLCND  (input) DOUBLE PRECISION
   *          Ratio of the smallest C(i) to the largest C(i).
   *
   *  AMAX    (input) DOUBLE PRECISION
   *          Absolute value of largest matrix entry.
   *
   *  EQUED   (output) CHARACTER*1
   *          Specifies the form of equilibration that was done.
   *          = 'N':  No equilibration
   *          = 'R':  Row equilibration, i.e., A has been premultiplied by
   *                  diag(R).
   *          = 'C':  Column equilibration, i.e., A has been postmultiplied
   *                  by diag(C).
   *          = 'B':  Both row and column equilibration, i.e., A has been
   *                  replaced by diag(R) * A * diag(C).
   *
   *  Internal Parameters
   *  ===================
   *
   *  THRESH is a threshold value used to decide if row or column scaling
   *  should be done based on the ratio of the row or column scaling
   *  factors.  If ROWCND < THRESH, row scaling is done, and if
   *  COLCND < THRESH, column scaling is done.
   *
   *  LARGE and SMALL are threshold values used to decide if row scaling
   *  should be done based on the absolute size of the largest matrix
   *  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(slaqge)( integer   * M,
                        integer   * N,
                        real      * A,
                        integer   * ldA,
                        real      * R,
                        real      * C,
                        real      * ROWCND,
                        real      * COLCND,
                        real      * AMAX,
                        character * EQUED ) ;

    void
    LAPACKNAME(dlaqge)( integer    * M,
                        integer    * N,
                        doublereal * A,
                        integer    * ldA,
                        doublereal * R,
                        doublereal * C,
                        doublereal * ROWCND,
                        doublereal * COLCND,
                        doublereal * AMAX,
                        character  * EQUED ) ;
  }
  #endif


  #ifdef ALGLIN_USE_OPENBLAS
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
              EquilibrationType & equ ) ;
  #endif

  inline
  void
  laqge( integer             M,
         integer             N,
         real                A[],
         integer             LDA,
         real const          R[],
         real const          C[],
         real                ROWCND,
         real                COLCND,
         real                AMAX,
         EquilibrationType & equ ) {
    #if defined(ALGLIN_USE_OPENBLAS)
    laqge_tmpl<real>( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, equ ) ;
    #elif defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(slaqge)( &M, &N, A, &LDA,
                        const_cast<real*>(R),
                        const_cast<real*>(C),
                        &ROWCND, &COLCND, &AMAX,
                        const_cast<character*>(equilibrate_blas[equ]) ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(slaqge)( &M, &N, A, &LDA,
                         const_cast<real*>(R),
                         const_cast<real*>(C),
                         &ROWCND, &COLCND, &AMAX,
                         const_cast<character*>(equilibrate_blas[equ]) ) ;
    #else
    LAPACKNAME(slaqge)( &M, &N, A, &LDA,
                        const_cast<real*>(R),
                        const_cast<real*>(C),
                        &ROWCND, &COLCND, &AMAX,
                        const_cast<character*>(equilibrate_blas[equ]) ) ;
    #endif
  }

  inline
  void
  laqge( integer             M,
         integer             N,
         doublereal          A[],
         integer             LDA,
         doublereal const    R[],
         doublereal const    C[],
         doublereal          ROWCND,
         doublereal          COLCND,
         doublereal          AMAX,
         EquilibrationType & equ ) {
    #if defined(ALGLIN_USE_OPENBLAS)
    laqge_tmpl<doublereal>( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, equ ) ;
    #elif defined(ALGLIN_USE_ATLAS)
    LAPACKNAME(dlaqge)( &M, &N, A, &LDA,
                        const_cast<doublereal*>(R),
                        const_cast<doublereal*>(C),
                        &ROWCND, &COLCND, &AMAX,
                        const_cast<character*>(equilibrate_blas[equ]) ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dlaqge)( &M, &N, A, &LDA,
                         const_cast<doublereal*>(R),
                         const_cast<doublereal*>(C),
                         &ROWCND, &COLCND, &AMAX,
                         const_cast<character*>(equilibrate_blas[equ]) ) ;
    #else
    LAPACKNAME(dlaqge)( &M, &N, A, &LDA,
                        const_cast<doublereal*>(R),
                        const_cast<doublereal*>(C),
                        &ROWCND, &COLCND, &AMAX,
                        const_cast<character*>(equilibrate_blas[equ]) ) ;
    #endif
  }

  /*
  //  #          #    ######     #     #####  #    #
  //  #         # #   #     #   # #   #     # #   #
  //  #        #   #  #     #  #   #  #       #  #
  //  #       #     # ######  #     # #       ###
  //  #       ####### #       ####### #       #  #
  //  #       #     # #       #     # #     # #   #
  //  ####### #     # #       #     #  #####  #    #
  */

  /*
  //    __ _  ___  ___ ___  _ __  _   _
  //   / _` |/ _ \/ __/ _ \| '_ \| | | |
  //  | (_| |  __/ (_| (_) | |_) | |_| |
  //   \__, |\___|\___\___/| .__/ \__, |
  //   |___/               |_|    |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DLACPY copies all or part of a two-dimensional matrix A to another
   *  matrix B.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          Specifies the part of the matrix A to be copied to B.
   *          = 'U':      Upper triangular part
   *          = 'L':      Lower triangular part
   *          Otherwise:  All of the matrix A
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The m by n matrix A.  If UPLO = 'U', only the upper triangle
   *          or trapezoid is accessed; if UPLO = 'L', only the lower
   *          triangle or trapezoid is accessed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
   *          On exit, B = A in the locations specified by UPLO.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,M).
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  LAPACKNAME(slacpy)( character        const   UPLO[], 
                      integer          const * M, 
                      integer          const * N, 
                      single_precision const   A[], 
                      integer          const * LDA, 
                      single_precision         B[], 
                      integer          const * LDB ) ;
  void
  LAPACKNAME(dlacpy)( character        const   UPLO[], 
                      integer          const * M, 
                      integer          const * N, 
                      double_precision const   A[], 
                      integer          const * LDA, 
                      double_precision         B[], 
                      integer          const * LDB ) ;
  }
  #endif

  inline
  integer
  gecopy( integer    M,
          integer    N,
          real const A[],
          integer    LDA,
          real       B[],
          integer    LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { for ( integer i = 0 ; i < N ; ++i ) alglin::copy( M, A+i*LDA, 1, B+i*LDB, 1 ) ; return 0 ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(slacpy)( const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB ) ; return 0 ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(slacpy)( const_cast<character*>("A"), &M, &N,
                                const_cast<real*>(A), &LDA, B, &LDB ) ; }
  #else
  { LAPACKNAME(slacpy)( const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB ) ; return 0 ; }
  #endif

  inline
  integer
  gecopy( integer          M,
          integer          N,
          doublereal const A[],
          integer          LDA,
          doublereal       B[],
          integer          LDB )
  #if defined(ALGLIN_USE_ATLAS)
  { for ( integer i = 0 ; i < N ; ++i ) alglin::copy( M, A+i*LDA, 1, B+i*LDB, 1 ) ; return 0 ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(dlacpy)( const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB ) ; return 0 ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlacpy)( const_cast<character*>("A"), &M, &N,
                                const_cast<doublereal*>(A), &LDA, B, &LDB ) ; }
  #else
  { LAPACKNAME(dlacpy)( const_cast<character*>("A"), &M, &N, A, &LDA, B, &LDB ) ; return 0 ; }
  #endif

  /*
  //    __ _  ___ _______ _ __ ___
  //   / _` |/ _ \_  / _ \ '__/ _ \
  //  | (_| |  __// /  __/ | | (_) |
  //   \__, |\___/___\___|_|  \___/
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
   *  ALPHA on the offdiagonals.
   *
   *  Arguments
   *  =========
   *
   *  UPLO    (input) CHARACTER*1
   *          Specifies the part of the matrix A to be set.
   *          = 'U':      Upper triangular part is set; the strictly lower
   *                      triangular part of A is not changed.
   *          = 'L':      Lower triangular part is set; the strictly upper
   *                      triangular part of A is not changed.
   *          Otherwise:  All of the matrix A is set.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  ALPHA   (input) DOUBLE PRECISION
   *          The constant to which the offdiagonal elements are to be set.
   *
   *  BETA    (input) DOUBLE PRECISION
   *          The constant to which the diagonal elements are to be set.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On exit, the leading m-by-n submatrix of A is set as follows:
   *
   *          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
   *          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
   *          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
   *
   *          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
  \*/
  #ifdef ALGLIN_USE_LAPACK
  extern "C" {
  using namespace blas_type ;
  void
  LAPACKNAME(slaset)( character        const   UPLO[], 
                      integer          const * M, 
                      integer          const * N, 
                      single_precision const * ALPHA, 
                      single_precision const * BETA, 
                      single_precision         A[], 
                      integer          const * LDA ) ;
  void
  LAPACKNAME(dlaset)( character        const   UPLO[], 
                      integer          const * M, 
                      integer          const * N, 
                      double_precision const * ALPHA, 
                      double_precision const * BETA, 
                      double_precision         A[], 
                      integer          const * LDA ) ;
  }
  #endif

  inline
  void
  gezero( integer M,
          integer N,
          real    A[],
          integer LDA )
  #if defined(ALGLIN_USE_ATLAS)
  { real zero = 0 ; for ( integer i = 0 ; i < N ; ++i ) alglin::copy( M, &zero, 0, A+i*LDA, 1 ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { real zero = 0 ; LAPACK_NAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { real zero = 0 ; CLAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #else
  { real zero = 0 ; LAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #endif

  inline
  void
  gezero( integer    M,
          integer    N,
          doublereal A[],
          integer    LDA )
  #if defined(ALGLIN_USE_ATLAS)
  { doublereal zero = 0 ; for ( integer i = 0 ; i < N ; ++i ) alglin::copy( M, &zero, 0, A+i*LDA, 1 ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { doublereal zero = 0 ; LAPACK_NAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { doublereal zero = 0 ; CLAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #else
  { doublereal zero = 0 ; LAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &zero, A, &LDA ) ; }
  #endif

  /*
  //             _     _
  //   __ _  ___(_) __| |
  //  / _` |/ _ \ |/ _` |
  // | (_| |  __/ | (_| |
  //  \__, |\___|_|\__,_|
  //  |___/
  */
  inline
  void
  geid( integer M,
        integer N,
        real    A[],
        integer LDA,
        real    diag = 1 )
  #if defined(ALGLIN_USE_ATLAS)
  { alglin::gezero( M, N, A, LDA ) ;
    integer MN = min_index(M,N) ;
    alglin::copy( MN, &diag, 0, A, LDA+1 ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { real zero = 0 ; LAPACK_NAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { real zero = 0 ; CLAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #else
  { real zero = 0 ; LAPACKNAME(slaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #endif

  inline
  void
  geid( integer    M,
        integer    N,
        doublereal A[],
        integer    LDA,
        doublereal diag = 1 )
  #if defined(ALGLIN_USE_ATLAS)
  { alglin::gezero( M, N, A, LDA ) ;
    integer MN = min_index(M,N) ;
    alglin::copy( MN, &diag, 0, A, LDA+1 ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { doublereal zero = 0 ; LAPACK_NAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { doublereal zero = 0 ; CLAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #else
  { doublereal zero = 0 ; LAPACKNAME(dlaset)( const_cast<character*>("A"), &M, &N, &zero, &diag, A, &LDA ) ; }
  #endif

  /*
  //                        _     _
  //    __ _  ___  __ _  __| | __| |
  //   / _` |/ _ \/ _` |/ _` |/ _` |
  //  | (_| |  __/ (_| | (_| | (_| |
  //   \__, |\___|\__,_|\__,_|\__,_|
  //   |___/
  */
  /*
   * C <- alpha * A + beta * B
   */

  #ifdef __GCC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wpartial-availability"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wpartial-availability"
  #endif

  inline
  void
  geadd( integer    M,
         integer    N,
         real       alpha,
         real const A[], integer LDA,
         real       beta,
         real const B[], integer LDB,
         real       C[], integer LDC )
  #if defined(ALGLIN_USE_ACCELERATE) && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_10
  { appleblas_sgeadd( CblasColMajor, CblasNoTrans, CblasNoTrans,
                      M, N, alpha, A, LDA, beta, B, LDB, C, LDC ) ; }
  #else
  {
    integer ierr = gecopy( M, N, B, LDB, C, LDC ) ;
    ALGLIN_ASSERT( ierr == 0, "geadd, ierr = " << ierr ) ;
    real const * Aj = A ;
    real       * Cj = C ;
    for ( integer j = 0 ; j < N ; ++j, Aj += LDA, Cj += LDC ) {
      alglin::scal( M, beta,  Cj, 1 ) ;
      alglin::axpy( M, alpha, Aj, 1, Cj, 1 ) ;
    }
  }
  #endif

  inline
  void
  geadd( integer          M,
         integer          N,
         doublereal       alpha,
         doublereal const A[], integer LDA,
         doublereal       beta,
         doublereal const B[], integer LDB,
         doublereal       C[], integer LDC )
  #if defined(ALGLIN_USE_ACCELERATE) && MAC_OS_X_VERSION_MIN_REQUIRED >= MAC_OS_X_VERSION_10_10
  { appleblas_dgeadd( CblasColMajor, CblasNoTrans, CblasNoTrans,
                      M, N, alpha, A, LDA, beta, B, LDB, C, LDC ) ; }
  #else
  {
    integer ierr = gecopy( M, N, B, LDB, C, LDC ) ;
    ALGLIN_ASSERT( ierr == 0, "geadd, ierr = " << ierr ) ;
    doublereal const * Aj = A ;
    doublereal       * Cj = C ;
    for ( integer j = 0 ; j < N ; ++j, Aj += LDA, Cj += LDC ) {
      alglin::scal( M, beta,  Cj, 1 ) ;
      alglin::axpy( M, alpha, Aj, 1, Cj, 1 ) ;
    }
  }
  #endif

  #ifdef __GCC__
  #pragma GCC diagnostic pop
  #endif
  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif

  /*
  //   __  __       _        _        _   _
  //  |  \/  | __ _| |_ _ __(_)_  __ | \ | | ___  _ __ _ __ ___
  //  | |\/| |/ _` | __| '__| \ \/ / |  \| |/ _ \| '__| '_ ` _ \
  //  | |  | | (_| | |_| |  | |>  <  | |\  | (_) | |  | | | | | |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\ |_| \_|\___/|_|  |_| |_| |_|
  */
  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    doublereal
    LAPACKNAME(dlange)( character  const   NORM[],
                        integer    const * M,
                        integer    const * N,
                        doublereal const * A,
                        integer    const * LDA,
                        doublereal       * WORK ) ;

    real
    LAPACKNAME(slange)( character const   NORM[],
                        integer   const * M,
                        integer   const * N,
                        real      const * A,
                        integer   const * LDA,
                        real            * WORK ) ;

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
   *  the  infinity norm,  or the  element of  largest absolute value  of a
   *  real matrix A.
   *
   *  Description
   *  ===========
   *
   *  DLANGE returns the value
   *
   *     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
   *              (
   *              ( norm1(A),         NORM = '1', 'O' or 'o'
   *              (
   *              ( normI(A),         NORM = 'I' or 'i'
   *              (
   *              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
   *
   *  where  norm1  denotes the  one norm of a matrix (maximum column sum),
   *  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
   *  normF  denotes the  Frobenius norm of a matrix (square root of sum of
   *  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
   *
   *  Arguments
   *  =========
   *
   *  NORM    (input) CHARACTER*1
   *          Specifies the value to be returned in DLANGE as described
   *          above.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.  When M = 0,
   *          DLANGE is set to zero.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.  When N = 0,
   *          DLANGE is set to zero.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          The m by n matrix A.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(M,1).
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
   *          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
   *          referenced.
  \*/

  inline
  real
  normInf( integer    N,
           integer    M,
           real const A[],
           integer    LDA )
  { real res = 0 ;
    for ( integer i=0 ; i < N ; ++i ) {
      real normcol = alglin::asum(M,A+i,LDA) ;
      if ( normcol > res ) res = normcol ;
    }
    return res ;
  }

  inline
  doublereal
  normInf( integer          N,
           integer          M,
           doublereal const A[],
           integer          LDA )
  { doublereal res = 0 ;
    for ( integer i=0 ; i < N ; ++i ) {
      doublereal normcol = alglin::asum(M,A+i,LDA) ;
      if ( normcol > res ) res = normcol ;
    }
    return res ;
  }

  ///////////////////

  inline
  real
  norm1( integer    N,
         integer    M,
         real const A[],
         integer    LDA )
  #if defined(ALGLIN_USE_ATLAS)
  { real res = 0 ;
    for ( integer i=0 ; i < M ; ++i ) {
      real normcol = alglin::asum(N,A+i*LDA,1) ;
      if ( normcol > res ) res = normcol ;
    }
    return res ;
  }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(slange)( const_cast<character*>("1"), &N, &M, A, &LDA, nullptr ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return real(CLAPACKNAME(slange)( const_cast<character*>("1"), &N, &M,
                                     const_cast<real*>(A), &LDA, nullptr ) ) ; }
  #else
  { return LAPACKNAME(slange)( const_cast<character*>("1"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  inline
  doublereal
  norm1( integer          N,
         integer          M,
         doublereal const A[],
         integer          LDA )
  #if defined(ALGLIN_USE_ATLAS)
  { return_precision res = 0 ;
    for ( integer i=0 ; i < M ; ++i ) {
      return_precision normcol = alglin::asum(N,A+i*LDA,1) ;
      if ( normcol > res ) res = normcol ;
    }
    return res ;
  }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(dlange)( const_cast<character*>("1"), &N, &M, A, &LDA, nullptr ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlange)( const_cast<character*>("1"), &N, &M,
                                const_cast<doublereal*>(A), &LDA, nullptr ) ; }
  #else
  { return LAPACKNAME(dlange)( const_cast<character*>("1"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  ///////////////////

  inline
  return_precision
  normF( integer    N,
         integer    M,
         real const A[],
         integer    LDA )
  #if defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(slange)( const_cast<character*>("F"), &N, &M,
                                const_cast<real*>(A), &LDA, nullptr ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(slange)( const_cast<character*>("F"), &N, &M, A, &LDA, nullptr ) ; }
  #else
  { return LAPACKNAME(slange)( const_cast<character*>("F"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  inline
  doublereal
  normF( integer          N,
         integer          M,
         doublereal const A[],
         integer          LDA )
  #if defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlange)( const_cast<character*>("F"), &N, &M,
                                const_cast<doublereal*>(A), &LDA, nullptr ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(dlange)( const_cast<character*>("F"), &N, &M, A, &LDA, nullptr ) ; }
  #else
  { return LAPACKNAME(dlange)( const_cast<character*>("F"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  ///////////////////

  inline
  real
  maxabs( integer    N,
          integer    M,
          real const A[],
          integer    LDA )
  #if defined(ALGLIN_USE_ACCELERATE)
  { return real(CLAPACKNAME(slange)( const_cast<character*>("M"), &N, &M,
                                     const_cast<real*>(A), &LDA, nullptr )) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(slange)( const_cast<character*>("M"), &N, &M, A, &LDA, nullptr ) ; }
  #else
  { return LAPACKNAME(slange)( const_cast<character*>("M"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  inline
  doublereal
  maxabs( integer          N,
          integer          M,
          doublereal const A[],
          integer          LDA )
  #if defined(ALGLIN_USE_ACCELERATE)
  { return CLAPACKNAME(dlange)( const_cast<character*>("M"), &N, &M,
                                const_cast<doublereal*>(A), &LDA, nullptr ) ; }
  #elif defined(ALGLIN_USE_OPENBLAS)
  { return LAPACK_NAME(dlange)( const_cast<character*>("M"), &N, &M, A, &LDA, nullptr ) ; }
  #else
  { return LAPACKNAME(dlange)( const_cast<character*>("M"), &N, &M, A, &LDA, nullptr ) ; }
  #endif

  /*
  //   _                _
  //  | | __ _ ___  ___| |
  //  | |/ _` / __|/ __| |
  //  | | (_| \__ \ (__| |
  //  |_|\__,_|___/\___|_|
  */
  /*\
   *
   *  Purpose
   *  =======
   *
   *  DLASCL multiplies the M by N real matrix A by the real scalar
   *  CTO/CFROM.  This is done without over/underflow as long as the final
   *  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
   *  A may be full, upper triangular, lower triangular, upper Hessenberg,
   *  or banded.
   *
   *  Arguments
   *  =========
   *
   *  TYPE    (input) CHARACTER*1
   *          TYPE indices the storage type of the input matrix.
   *          = 'G':  A is a full matrix.
   *          = 'L':  A is a lower triangular matrix.
   *          = 'U':  A is an upper triangular matrix.
   *          = 'H':  A is an upper Hessenberg matrix.
   *          = 'B':  A is a symmetric band matrix with lower bandwidth KL
   *                  and upper bandwidth KU and with the only the lower
   *                  half stored.
   *          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
   *                  and upper bandwidth KU and with the only the upper
   *                  half stored.
   *          = 'Z':  A is a band matrix with lower bandwidth KL and upper
   *                  bandwidth KU.
   *
   *  KL      (input) INTEGER
   *          The lower bandwidth of A.  Referenced only if TYPE = 'B',
   *          'Q' or 'Z'.
   *
   *  KU      (input) INTEGER
   *          The upper bandwidth of A.  Referenced only if TYPE = 'B',
   *          'Q' or 'Z'.
   *
   *  CFROM   (input) DOUBLE PRECISION
   *  CTO     (input) DOUBLE PRECISION
   *          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
   *          without over/underflow if the final result CTO*A(I,J)/CFROM
   *          can be represented without over/underflow.  CFROM must be
   *          nonzero.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
   *          storage type.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  INFO    (output) INTEGER
   *          0  - successful exit
   *          <0 - if INFO = -i, the i-th argument had an illegal value.
   *
   *  =====================================================================
  \*/
  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {
    void
    LAPACKNAME(dlascl)( character    TYPE[],
                        integer    * KL,
                        integer    * KU,
                        doublereal * FROM,
                        doublereal * TO,
                        integer    * M,
                        integer    * N,
                        doublereal * A,
                        integer    * LDA,
                        integer    * info ) ;

    void
    LAPACKNAME(slascl)( character    TYPE[],
                        integer    * KL,
                        integer    * KU,
                        real       * FROM,
                        real       * TO,
                        integer    * M,
                        integer    * N,
                        real       * A,
                        integer    * LDA,
                        integer    * info ) ;
  }
  #endif

  inline
  integer
  lascl( MatrixType const & TYPE,
         integer    KL,
         integer    KU,
         real       FROM,
         real       TO,
         integer    M,
         integer    N,
         real       A[],
         integer    LDA ) {
    integer info = 0 ;
    #if defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(slascl)( const_cast<character*>(mtype_blas[TYPE]),
                         &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    #else
    LAPACKNAME(slascl)( const_cast<character*>(mtype_blas[TYPE]),
                        &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  lascl( MatrixType const & TYPE,
         integer    KL,
         integer    KU,
         doublereal FROM,
         doublereal TO,
         integer    M,
         integer    N,
         doublereal A[],
         integer    LDA ) {
    integer info = 0 ;
    #if defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dlascl)( const_cast<character*>(mtype_blas[TYPE]),
                         &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info ) ;
    #elif defined(ALGLIN_USE_OPENBLAS)
    #else
    LAPACKNAME(dlascl)( const_cast<character*>(mtype_blas[TYPE]),
                        &KL, &KU, &FROM, &TO, &M, &N, A, &LDA, &info ) ;
    #endif
    return info ;
  }

  /*
  //   _____ _                            _
  //  | ____(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
  //  |  _| | |/ _` |/ _ \ '_ \ \ / / _` | | | | |/ _ \/ __|
  //  | |___| | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
  //  |_____|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
  //           |___/
  */
  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    // DGEEV - compute for an N-by-N real nonsymmetric matrix A,
    //         the eigenvalues and, optionally, the left and/or right
    //          eigenvectors

    void
    LAPACKNAME(dgeev)( character    jobvl[],
                       character    jobvr[],
                       integer    * N,
                       doublereal * A,
                       integer    * LDA,
                       doublereal * wr,
                       doublereal * wi,
                       doublereal * vl,
                       integer    * Lvl,
                       doublereal * vr,
                       integer    * Lvr,
                       doublereal * WORK,
                       integer    * Lwork,
                       integer    * info ) ;

    void
    LAPACKNAME(sgeev)( character   jobvl[],
                       character   jobvr[],
                       integer   * N,
                       real      * A,
                       integer   * LDA,
                       real      * wr,
                       real      * wi,
                       real      * vl,
                       integer   * Lvl,
                       real      * vr,
                       integer   * Lvr,
                       real      * WORK,
                       integer   * Lwork,
                       integer   * info ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
   *  eigenvalues and, optionally, the left and/or right eigenvectors.
   *
   *  The right eigenvector v(j) of A satisfies
   *                   A * v(j) = lambda(j) * v(j)
   *  where lambda(j) is its eigenvalue.
   *  The left eigenvector u(j) of A satisfies
   *                u(j)**T * A = lambda(j) * u(j)**T
   *  where u(j)**T denotes the transpose of u(j).
   *
   *  The computed eigenvectors are normalized to have Euclidean norm
   *  equal to 1 and largest component real.
   *
   *  Arguments
   *  =========
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N': left eigenvectors of A are not computed;
   *          = 'V': left eigenvectors of A are computed.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N': right eigenvectors of A are not computed;
   *          = 'V': right eigenvectors of A are computed.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A. N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N matrix A.
   *          On exit, A has been overwritten.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  WR      (output) DOUBLE PRECISION array, dimension (N)
   *  WI      (output) DOUBLE PRECISION array, dimension (N)
   *          WR and WI contain the real and imaginary parts,
   *          respectively, of the computed eigenvalues.  Complex
   *          conjugate pairs of eigenvalues appear consecutively
   *          with the eigenvalue having the positive imaginary part
   *          first.
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order
   *          as their eigenvalues.
   *          If JOBVL = 'N', VL is not referenced.
   *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
   *          the j-th column of VL.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
   *          u(j+1) = VL(:,j) - i*VL(:,j+1).
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the array VL.  LDVL >= 1; if
   *          JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order
   *          as their eigenvalues.
   *          If JOBVR = 'N', VR is not referenced.
   *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
   *          the j-th column of VR.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
   *          v(j+1) = VR(:,j) - i*VR(:,j+1).
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the array VR.  LDVR >= 1; if
   *          JOBVR = 'V', LDVR >= N.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
   *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
   *          performance, LWORK must generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if INFO = i, the QR algorithm failed to compute all the
   *                eigenvalues, and no eigenvectors have been computed;
   *                elements i+1:N of WR and WI contain eigenvalues which
   *                have converged.
  \*/

  inline
  integer
  geev( bool    jobvl, // false = do not compute the left generalized eigenvectors
        bool    jobvr, // false = do not compute the right generalized eigenvectors
        integer n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
        real    A[],
        integer lda,
        real    wr[],
        real    wi[],
        real    vl[],
        integer ldvl,
        real    vr[],
        integer ldvr,
        real    work[],
        integer lwork )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(sgeev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n, A, &lda, wr, wi,
                        vl, &ldvl,
                        vr, &ldvr, work, &lwork, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0 ;
    CLAPACKNAME(sgeev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n,
                        const_cast<real*>(A), &lda, wr, wi,
                        vl, &ldvl,
                        vr, &ldvr,
                        work, &lwork, &INFO ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(sgeev)( const_cast<character*>(jobvl?"V":"N"),
                       const_cast<character*>(jobvr?"V":"N"),
                       &n, A, &lda, wr, wi,
                       vl, &ldvl,
                       vr, &ldvr,
                       work, &lwork, &INFO ) ;
    return INFO ;
  }
  #endif

  inline
  integer
  geev( bool       jobvl, // false = do not compute the left generalized eigenvectors
        bool       jobvr, // false = do not compute the right generalized eigenvectors
        integer    n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
        doublereal A[],
        integer    lda,
        doublereal wr[],
        doublereal wi[],
        doublereal vl[],
        integer    ldvl,
        doublereal vr[],
        integer    ldvr,
        doublereal work[],
        integer    lwork )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(dgeev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n, A, &lda, wr, wi,
                        vl, &ldvl,
                        vr, &ldvr, work, &lwork, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(dgeev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n,
                        const_cast<doublereal*>(A), &lda, wr, wi,
                        vl, &ldvl,
                        vr, &ldvr,
                        work, &lwork, &INFO ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(dgeev)( const_cast<character*>(jobvl?"V":"N"),
                       const_cast<character*>(jobvr?"V":"N"),
                       &n, A, &lda, wr, wi,
                       vl, &ldvl,
                       vr, &ldvr,
                       work, &lwork, &INFO ) ;
    return INFO ;
  }
  #endif

  //////////////////////////////////////////////////////////////////////////////

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    //  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
    //  the generalized eigenvalues, and optionally, the left and/or right
    //  generalized eigenvectors.
    //
    //  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
    //  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
    //  singular. It is usually represented as the pair (alpha,beta), as
    //  there is a reasonable interpretation for beta=0, and even for both
    //  being zero.


    void
    LAPACKNAME(dggev)( character    jobvl[],
                       character    jobvr[],
                       integer    * N,
                       doublereal * A,
                       integer    * LDA,
                       doublereal * B,
                       integer    * LDB,
                       doublereal * alphar,
                       doublereal * alphai,
                       doublereal * beta,
                       doublereal * vl,
                       integer    * Lvl,
                       doublereal * vr,
                       integer    * Lvr,
                       doublereal * WORK,
                       integer    * Lwork,
                       integer    * info ) ;

    void
    LAPACKNAME(sggev)( character   jobvl[],
                       character   jobvr[],
                       integer   * N,
                       real      * A,
                       integer   * LDA,
                       real      * B,
                       integer   * LDB,
                       real      * alphar,
                       real      * alphai,
                       real      * beta,
                       real      * vl,
                       integer   * Lvl,
                       real      * vr,
                       integer   * Lvr,
                       real      * WORK,
                       integer   * Lwork,
                       integer   * info ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j).
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B .
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          alpha/beta.  However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector is scaled so the largest component has
   *          abs(real part)+abs(imag. part)=1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector is scaled so the largest component has
   *          abs(real part)+abs(imag. part)=1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,8*N).
   *          For good performance, LWORK must generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
  \*/

  inline
  integer
  ggev( bool    jobvl, // false = do not compute the left generalized eigenvectors
        bool    jobvr, // false = do not compute the right generalized eigenvectors
        integer n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
        real    A[],
        integer lda,
        real    B[],
        integer ldb,
        real    alphar[],
        real    alphai[],
        real    beta[],
        real    vl[],
        integer ldvl,
        real    vr[],
        integer ldvr,
        real    work[],
        integer lwork )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(sggev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n, A, &lda, B, &ldb, alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr, work, &lwork, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(sggev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n,
                        const_cast<real*>(A), &lda,
                        const_cast<real*>(B), &ldb,
                        alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr, work, &lwork, &INFO ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(sggev)( const_cast<character*>(jobvl?"V":"N"),
                       const_cast<character*>(jobvr?"V":"N"),
                       &n, A, &lda, B, &ldb, alphar, alphai, beta,
                       vl, &ldvl, vr, &ldvr, work, &lwork, &INFO ) ;
    return INFO ;
  }
  #endif


  inline
  integer
  ggev( bool       jobvl, // false = do not compute the left generalized eigenvectors
        bool       jobvr, // false = do not compute the right generalized eigenvectors
        integer    n,     // The order of the matrices A, B, VL, and VR.  N >= 0.
        doublereal A[],
        integer	   lda,
        doublereal B[],
        integer    ldb,
        doublereal alphar[],
        doublereal alphai[],
        doublereal beta[],
        doublereal vl[],
        integer    ldvl,
        doublereal vr[],
        integer    ldvr,
        doublereal work[],
        integer    lwork )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(dggev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n, A, &lda, B, &ldb, alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr, work, &lwork, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(dggev)( const_cast<character*>(jobvl?"V":"N"),
                        const_cast<character*>(jobvr?"V":"N"),
                        &n,
                        const_cast<doublereal*>(A), &lda,
                        const_cast<doublereal*>(B), &ldb,
                        alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr, work, &lwork, &INFO ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(dggev)( const_cast<character*>(jobvl?"V":"N"),
                       const_cast<character*>(jobvr?"V":"N"),
                       &n, A, &lda, B, &ldb, alphar, alphai, beta,
                       vl, &ldvl, vr, &ldvr, work, &lwork, &INFO ) ;
    return INFO ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  Optionally also, it computes a balancing transformation to improve
   *  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
   *  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
   *  the eigenvalues (RCONDE), and reciprocal condition numbers for the
   *  right eigenvectors (RCONDV).
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j) .
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B.
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  BALANC  (input) CHARACTER*1
   *          Specifies the balance option to be performed.
   *          = 'N':  do not diagonally scale or permute;
   *          = 'P':  permute only;
   *          = 'S':  scale only;
   *          = 'B':  both permute and scale.
   *          Computed reciprocal condition numbers will be for the
   *          matrices after permuting and/or balancing. Permuting does
   *          not change condition numbers (in exact arithmetic), but
   *          balancing does.
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  SENSE   (input) CHARACTER*1
   *          Determines which reciprocal condition numbers are computed.
   *          = 'N': none are computed;
   *          = 'E': computed for eigenvalues only;
   *          = 'V': computed for eigenvectors only;
   *          = 'B': computed for eigenvalues and eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then A contains the first part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then B contains the second part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  ILO     (output) INTEGER
   *  IHI     (output) INTEGER
   *          ILO and IHI are integer values such that on exit
   *          A(i,j) = 0 and B(i,j) = 0 if i > j and
   *          j = 1,...,ILO-1 or i = IHI+1,...,N.
   *          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
   *
   *  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the left side of A and B.  If PL(j) is the index of the
   *          row interchanged with row j, and DL(j) is the scaling
   *          factor applied to row j, then
   *            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
   *                      = DL(j)  for j = ILO,...,IHI
   *                      = PL(j)  for j = IHI+1,...,N.
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the right side of A and B.  If PR(j) is the index of the
   *          column interchanged with column j, and DR(j) is the scaling
   *          factor applied to column j, then
   *            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
   *                      = DR(j)  for j = ILO,...,IHI
   *                      = PR(j)  for j = IHI+1,...,N
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  ABNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix A.
   *
   *  BBNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix B.
   *
   *  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'E' or 'B', the reciprocal condition numbers of
   *          the eigenvalues, stored in consecutive elements of the array.
   *          For a complex conjugate pair of eigenvalues two consecutive
   *          elements of RCONDE are set to the same value. Thus RCONDE(j),
   *          RCONDV(j), and the j-th columns of VL and VR all correspond
   *          to the j-th eigenpair.
   *          If SENSE = 'N or 'V', RCONDE is not referenced.
   *
   *  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'V' or 'B', the estimated reciprocal condition
   *          numbers of the eigenvectors, stored in consecutive elements
   *          of the array. For a complex eigenvector two consecutive
   *          elements of RCONDV are set to the same value. If the
   *          eigenvalues cannot be reordered to compute RCONDV(j),
   *          RCONDV(j) is set to 0; this can only occur when the true
   *          value would be very small anyway.
   *          If SENSE = 'N' or 'E', RCONDV is not referenced.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= max(1,2*N).
   *          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',
   *          LWORK >= max(1,6*N).
   *          If SENSE = 'E' or 'B', LWORK >= max(1,10*N).
   *          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  IWORK   (workspace) INTEGER array, dimension (N+6)
   *          If SENSE = 'E', IWORK is not referenced.
   *
   *  BWORK   (workspace) LOGICAL array, dimension (N)
   *          If SENSE = 'N', BWORK is not referenced.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
   *  Further Details
   *  ===============
   *
   *  Balancing a matrix pair (A,B) includes, first, permuting rows and
   *  columns to isolate eigenvalues, second, applying diagonal similarity
   *  transformation to the rows and columns to make the rows and columns
   *  as close in norm as possible. The computed reciprocal condition
   *  numbers correspond to the balanced matrix. Permuting rows and columns
   *  will not change the condition numbers (in exact arithmetic) but
   *  diagonal scaling will.  For further explanation of balancing, see
   *  section 4.11.1.2 of LAPACK Users' Guide.
   *
   *  An approximate error bound on the chordal distance between the i-th
   *  computed generalized eigenvalue w and the corresponding exact
   *  eigenvalue lambda is
   *
   *       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
   *
   *  An approximate error bound for the angle between the i-th computed
   *  eigenvector VL(i) or VR(i) is given by
   *
   *       EPS * norm(ABNRM, BBNRM) / DIF(i).
   *
   *  For further explanation of the reciprocal condition numbers RCONDE
   *  and RCONDV, see section 4.11 of LAPACK User's Guide.
   *
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    // 'M' or 'm'            maxabs
    // '1', 'O' or 'o'       norm1
    // 'I' or 'i'            normI
    // 'F', 'f', 'E' or 'e'  normF

    void
    LAPACKNAME(dggevx)(character  const   balanc[],
                       character  const   jobvl[],
                       character  const   jobvr[],
                       character  const   sense[],
                       integer    const * N,
                       doublereal       * A,
                       integer    const * LDA,
                       doublereal       * B,
                       integer    const * LDB,
                       doublereal       * alphar,
                       doublereal       * alphai,
                       doublereal       * beta,
                       doublereal       * vl,
                       integer    const * Lvl,
                       doublereal       * vr,
                       integer    const * Lvr,
                       integer          * ILO,
                       integer          * IHI,
                       doublereal       * LSCALE,
                       doublereal       * RSCALE,
                       doublereal       * ABNRM,
                       doublereal       * BBNRM,
                       doublereal       * RCONDE,
                       doublereal       * RCONDV,
                       doublereal       * WORK,
                       integer    const * Lwork,
                       integer          * Iwork,
                       integer          * Bwork,
                       integer          * info ) ;

    void
    LAPACKNAME(sggevx)(character  const   balanc[],
                       character  const   jobvl[],
                       character  const   jobvr[],
                       character  const   sense[],
                       integer    const  * N,
                       real              * A,
                       integer    const  * LDA,
                       real              * B,
                       integer    const  * LDB,
                       real              * alphar,
                       real              * alphai,
                       real              * beta,
                       real              * vl,
                       integer    const  * Lvl,
                       real              * vr,
                       integer    const  * Lvr,
                       integer           * ILO,
                       integer           * IHI,
                       real              * LSCALE,
                       real              * RSCALE,
                       real              * ABNRM,
                       real              * BBNRM,
                       real              * RCONDE,
                       real              * RCONDV,
                       real              * WORK,
                       integer     const * Lwork,
                       integer           * Iwork,
                       integer           * Bwork,
                       integer           * info  ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGGEVX computes for a pair of N-by-N real nonsymmetric matrices (A,B)
   *  the generalized eigenvalues, and optionally, the left and/or right
   *  generalized eigenvectors.
   *
   *  Optionally also, it computes a balancing transformation to improve
   *  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
   *  LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
   *  the eigenvalues (RCONDE), and reciprocal condition numbers for the
   *  right eigenvectors (RCONDV).
   *
   *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
   *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
   *  singular. It is usually represented as the pair (alpha,beta), as
   *  there is a reasonable interpretation for beta=0, and even for both
   *  being zero.
   *
   *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   A * v(j) = lambda(j) * B * v(j) .
   *
   *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
   *  of (A,B) satisfies
   *
   *                   u(j)**H * A  = lambda(j) * u(j)**H * B.
   *
   *  where u(j)**H is the conjugate-transpose of u(j).
   *
   *
   *  Arguments
   *  =========
   *
   *  BALANC  (input) CHARACTER*1
   *          Specifies the balance option to be performed.
   *          = 'N':  do not diagonally scale or permute;
   *          = 'P':  permute only;
   *          = 'S':  scale only;
   *          = 'B':  both permute and scale.
   *          Computed reciprocal condition numbers will be for the
   *          matrices after permuting and/or balancing. Permuting does
   *          not change condition numbers (in exact arithmetic), but
   *          balancing does.
   *
   *  JOBVL   (input) CHARACTER*1
   *          = 'N':  do not compute the left generalized eigenvectors;
   *          = 'V':  compute the left generalized eigenvectors.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N':  do not compute the right generalized eigenvectors;
   *          = 'V':  compute the right generalized eigenvectors.
   *
   *  SENSE   (input) CHARACTER*1
   *          Determines which reciprocal condition numbers are computed.
   *          = 'N': none are computed;
   *          = 'E': computed for eigenvalues only;
   *          = 'V': computed for eigenvectors only;
   *          = 'B': computed for eigenvalues and eigenvectors.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A, B, VL, and VR.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the matrix A in the pair (A,B).
   *          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then A contains the first part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the matrix B in the pair (A,B).
   *          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
   *          or both, then B contains the second part of the real Schur
   *          form of the "balanced" versions of the input A and B.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of B.  LDB >= max(1,N).
   *
   *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
   *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
   *  BETA    (output) DOUBLE PRECISION array, dimension (N)
   *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
   *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
   *          the j-th eigenvalue is real; if positive, then the j-th and
   *          (j+1)-st eigenvalues are a complex conjugate pair, with
   *          ALPHAI(j+1) negative.
   *
   *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
   *          may easily over- or underflow, and BETA(j) may even be zero.
   *          Thus, the user should avoid naively computing the ratio
   *          ALPHA/BETA. However, ALPHAR and ALPHAI will be always less
   *          than and usually comparable with norm(A) in magnitude, and
   *          BETA always less than and usually comparable with norm(B).
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVL = 'N'.
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the matrix VL. LDVL >= 1, and
   *          if JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order as
   *          their eigenvalues. If the j-th eigenvalue is real, then
   *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
   *          (j+1)-th eigenvalues form a complex conjugate pair, then
   *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
   *          Each eigenvector will be scaled so the largest component have
   *          abs(real part) + abs(imag. part) = 1.
   *          Not referenced if JOBVR = 'N'.
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the matrix VR. LDVR >= 1, and
   *          if JOBVR = 'V', LDVR >= N.
   *
   *  ILO     (output) INTEGER
   *  IHI     (output) INTEGER
   *          ILO and IHI are integer values such that on exit
   *          A(i,j) = 0 and B(i,j) = 0 if i > j and
   *          j = 1,...,ILO-1 or i = IHI+1,...,N.
   *          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
   *
   *  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the left side of A and B.  If PL(j) is the index of the
   *          row interchanged with row j, and DL(j) is the scaling
   *          factor applied to row j, then
   *            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
   *                      = DL(j)  for j = ILO,...,IHI
   *                      = PL(j)  for j = IHI+1,...,N.
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
   *          Details of the permutations and scaling factors applied
   *          to the right side of A and B.  If PR(j) is the index of the
   *          column interchanged with column j, and DR(j) is the scaling
   *          factor applied to column j, then
   *            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
   *                      = DR(j)  for j = ILO,...,IHI
   *                      = PR(j)  for j = IHI+1,...,N
   *          The order in which the interchanges are made is N to IHI+1,
   *          then 1 to ILO-1.
   *
   *  ABNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix A.
   *
   *  BBNRM   (output) DOUBLE PRECISION
   *          The one-norm of the balanced matrix B.
   *
   *  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'E' or 'B', the reciprocal condition numbers of
   *          the eigenvalues, stored in consecutive elements of the array.
   *          For a complex conjugate pair of eigenvalues two consecutive
   *          elements of RCONDE are set to the same value. Thus RCONDE(j),
   *          RCONDV(j), and the j-th columns of VL and VR all correspond
   *          to the j-th eigenpair.
   *          If SENSE = 'N or 'V', RCONDE is not referenced.
   *
   *  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
   *          If SENSE = 'V' or 'B', the estimated reciprocal condition
   *          numbers of the eigenvectors, stored in consecutive elements
   *          of the array. For a complex eigenvector two consecutive
   *          elements of RCONDV are set to the same value. If the
   *          eigenvalues cannot be reordered to compute RCONDV(j),
   *          RCONDV(j) is set to 0; this can only occur when the true
   *          value would be very small anyway.
   *          If SENSE = 'N' or 'E', RCONDV is not referenced.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= max(1,2*N).
   *          If BALANC = 'S' or 'B', or JOBVL = 'V', or JOBVR = 'V',
   *          LWORK >= max(1,6*N).
   *          If SENSE = 'E' or 'B', LWORK >= max(1,10*N).
   *          If SENSE = 'V' or 'B', LWORK >= 2*N*N+8*N+16.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  IWORK   (workspace) INTEGER array, dimension (N+6)
   *          If SENSE = 'E', IWORK is not referenced.
   *
   *  BWORK   (workspace) LOGICAL array, dimension (N)
   *          If SENSE = 'N', BWORK is not referenced.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          = 1,...,N:
   *                The QZ iteration failed.  No eigenvectors have been
   *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
   *                should be correct for j=INFO+1,...,N.
   *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
   *                =N+2: error return from DTGEVC.
   *
   *  Further Details
   *  ===============
   *
   *  Balancing a matrix pair (A,B) includes, first, permuting rows and
   *  columns to isolate eigenvalues, second, applying diagonal similarity
   *  transformation to the rows and columns to make the rows and columns
   *  as close in norm as possible. The computed reciprocal condition
   *  numbers correspond to the balanced matrix. Permuting rows and columns
   *  will not change the condition numbers (in exact arithmetic) but
   *  diagonal scaling will.  For further explanation of balancing, see
   *  section 4.11.1.2 of LAPACK Users' Guide.
   *
   *  An approximate error bound on the chordal distance between the i-th
   *  computed generalized eigenvalue w and the corresponding exact
   *  eigenvalue lambda is
   *
   *       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
   *
   *  An approximate error bound for the angle between the i-th computed
   *  eigenvector VL(i) or VR(i) is given by
   *
   *       EPS * norm(ABNRM, BBNRM) / DIF(i).
   *
   *  For further explanation of the reciprocal condition numbers RCONDE
   *  and RCONDV, see section 4.11 of LAPACK User's Guide.
   *
  \*/

  inline
  integer
  ggevx( BalanceType const & balanc,
         bool                jobvl,
         bool                jobvr,
         SenseType   const & sense,
         integer             n,
         real                A[],
         integer             lda,
         real                B[],
         integer             ldb,
         real                alphar[],
         real                alphai[],
         real                beta[],
         real                vl[],
         integer             ldvl,
         real                vr[],
         integer             ldvr,
         integer           & ilo,
         integer           & ihi,
         real                lscale[],
         real                rscale[],
         real              & abnrm,
         real              & bbnrm,
         real                rconde[],
         real                rcondv[],
         real                work[],
         integer             lwork,
         integer             iwork[],
         integer             bwork[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(sggevx)( const_cast<character*>(balance_blas[balanc]),
                         const_cast<character*>(jobvl?"V":"N"),
                         const_cast<character*>(jobvr?"V":"N"),
                         const_cast<character*>(sense_blas[sense]),
                         &n,
                         A, &lda,
                         B, &ldb,
                         alphar, alphai, beta,
                         vl, &ldvl, vr, &ldvr,
                         &ilo, &ihi,
                         lscale, rscale,
                         &abnrm, &bbnrm,
                         rconde, rcondv,
                         work, &lwork, iwork, bwork, &info) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(sggevx)( const_cast<character*>(balance_blas[balanc]),
                         const_cast<character*>(jobvl?"V":"N"),
                         const_cast<character*>(jobvr?"V":"N"),
                         const_cast<character*>(sense_blas[sense]),
                         &n,
                         const_cast<real*>(A), &lda,
                         const_cast<real*>(B), &ldb,
                         alphar, alphai, beta,
                         vl, &ldvl, vr, &ldvr,
                         &ilo, &ihi,
                         lscale, rscale,
                         &abnrm, &bbnrm,
                         rconde, rcondv,
                         work, &lwork, iwork, bwork, &INFO ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(sggevx)( balance_blas[balanc],
                        (jobvl?"V":"N"),
                        (jobvr?"V":"N"),
                        sense_blas[sense],
                        &n, A, &lda, B, &ldb, alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr,
                        &ilo, &ihi,
                        lscale, rscale,
                        &abnrm, &bbnrm,
                        rconde, rcondv,
                        work, &lwork, iwork, bwork, &INFO ) ;
    return INFO ;
  }
  #endif

  inline
  integer
  ggevx(BalanceType const & balanc,
        bool                jobvl,
        bool                jobvr,
        SenseType   const & sense,
        integer             n,
        doublereal          A[],
        integer	            lda,
        doublereal          B[],
        integer             ldb,
        doublereal          alphar[],
        doublereal          alphai[],
        doublereal          beta[],
        doublereal          vl[],
        integer             ldvl,
        doublereal          vr[],
        integer             ldvr,
        integer           & ilo,
        integer           & ihi,
        doublereal          lscale[],
        doublereal          rscale[],
        doublereal        & abnrm,
        doublereal        & bbnrm,
        doublereal          rconde[],
        doublereal          rcondv[],
        doublereal          work[],
        integer             lwork,
        integer             iwork[],
        integer             bwork[])
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(dggevx)( const_cast<character*>(balance_blas[balanc]),
                         const_cast<character*>(jobvl?"V":"N"),
                         const_cast<character*>(jobvr?"V":"N"),
                         const_cast<character*>(sense_blas[sense]),
                         &n,
                         A, &lda,
                         B, &ldb,
                         alphar, alphai, beta,
                         vl, &ldvl, vr, &ldvr,
                         &ilo, &ihi,
                         lscale, rscale,
                         &abnrm, &bbnrm,
                         rconde, rcondv,
                         work, &lwork, iwork, bwork, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer INFO = 0;
    CLAPACKNAME(dggevx)( const_cast<character*>(balance_blas[balanc]),
                         const_cast<character*>(jobvl?"V":"N"),
                         const_cast<character*>(jobvr?"V":"N"),
                         const_cast<character*>(sense_blas[sense]),
                         &n,
                         const_cast<doublereal*>(A), &lda,
                         const_cast<doublereal*>(B), &ldb,
                         alphar, alphai, beta,
                         vl, &ldvl, vr, &ldvr,
                         &ilo, &ihi,
                         lscale, rscale,
                         &abnrm, &bbnrm,
                         rconde, rcondv,
                         work, &lwork, iwork, bwork, &INFO  ) ;
    return INFO ;
  }
  #else
  { integer INFO = 0;
    LAPACKNAME(dggevx)( balance_blas[balanc],
                        (jobvl?"V":"N"),
                        (jobvr?"V":"N"),
                        sense_blas[sense],
                        &n, A, &lda, B, &ldb, alphar, alphai, beta,
                        vl, &ldvl, vr, &ldvr,
                        &ilo, &ihi,
                        lscale, rscale,
                        &abnrm, &bbnrm,
                        rconde, rcondv,
                        work, &lwork, iwork, bwork, &INFO ) ;
    return INFO ;
  }
  #endif

  /*
  //   ____  ______     __
  //  / ___||  _ \ \   / /
  //  \___ \| | | \ \ / /
  //   ___) | |_| |\ V /
  //  |____/|____/  \_/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESVD computes the singular value decomposition (SVD) of a real
   *  M-by-N matrix A, optionally computing the left and/or right singular
   *  vectors. The SVD is written
   *
   *       A = U * SIGMA * transpose(V)
   *
   *  where SIGMA is an M-by-N matrix which is zero except for its
   *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   *  are the singular values of A; they are real and non-negative, and
   *  are returned in descending order.  The first min(m,n) columns of
   *  U and V are the left and right singular vectors of A.
   *
   *  Note that the routine returns V**T, not V.
   *
   *  Arguments
   *  =========
   *
   *  JOBU    (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix U:
   *          = 'A':  all M columns of U are returned in array U:
   *          = 'S':  the first min(m,n) columns of U (the left singular
   *                  vectors) are returned in the array U;
   *          = 'O':  the first min(m,n) columns of U (the left singular
   *                  vectors) are overwritten on the array A;
   *          = 'N':  no columns of U (no left singular vectors) are
   *                  computed.
   *
   *  JOBVT   (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix
   *          V**T:
   *          = 'A':  all N rows of V**T are returned in the array VT;
   *          = 'S':  the first min(m,n) rows of V**T (the right singular
   *                  vectors) are returned in the array VT;
   *          = 'O':  the first min(m,n) rows of V**T (the right singular
   *                  vectors) are overwritten on the array A;
   *          = 'N':  no rows of V**T (no right singular vectors) are
   *                  computed.
   *
   *          JOBVT and JOBU cannot both be 'O'.
   *
   *  M       (input) INTEGER
   *          The number of rows of the input matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the input matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit,
   *          if JOBU = 'O',  A is overwritten with the first min(m,n)
   *                          columns of U (the left singular vectors,
   *                          stored columnwise);
   *          if JOBVT = 'O', A is overwritten with the first min(m,n)
   *                          rows of V**T (the right singular vectors,
   *                          stored rowwise);
   *          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
   *                          are destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A, sorted so that S(i) >= S(i+1).
   *
   *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   *          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
   *          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
   *          if JOBU = 'S', U contains the first min(m,n) columns of U
   *          (the left singular vectors, stored columnwise);
   *          if JOBU = 'N' or 'O', U is not referenced.
   *
   *  LDU     (input) INTEGER
   *          The leading dimension of the array U.  LDU >= 1; if
   *          JOBU = 'S' or 'A', LDU >= M.
   *
   *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   *          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
   *          V**T;
   *          if JOBVT = 'S', VT contains the first min(m,n) rows of
   *          V**T (the right singular vectors, stored rowwise);
   *          if JOBVT = 'N' or 'O', VT is not referenced.
   *
   *  LDVT    (input) INTEGER
   *          The leading dimension of the array VT.  LDVT >= 1; if
   *          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   *          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
   *          superdiagonal elements of an upper bidiagonal matrix B
   *          whose diagonal is in S (not necessarily sorted). B
   *          satisfies A = U * B * VT, so it has the same singular values
   *          as A, and singular vectors related by U and VT.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if DBDSQR did not converge, INFO specifies how many
   *                superdiagonals of an intermediate bidiagonal form B
   *                did not converge to zero. See the description of WORK
   *                above for details.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(dgesvd)( character const   JOBU[],
                        character const   JOBVT[],
                        integer   const * M,
                        integer   const * N,
                        doublereal      * A,
                        integer   const * LDA,
                        doublereal      * S,
                        doublereal      * U,
                        integer   const * LDU,
                        doublereal      * VT,
                        integer   const * LDVT,
                        doublereal      * WORK,
                        integer   const * LWORK,
                        integer         * info ) ;

    void
    LAPACKNAME(sgesvd)( character const   JOBU[],
                        character const   JOBVT[],
                        integer   const * M,
                        integer   const * N,
                        real            * A,
                        integer   const * LDA,
                        real            * S,
                        real            * U,
                        integer   const * LDU,
                        real            * VT,
                        integer   const * LDVT,
                        real            * WORK,
                        integer   const * LWORK,
                        integer         * info  ) ;
  }
  #endif

  inline
  integer
  gesvd( JobType const & JOBU,
         JobType const & JOBVT,
         integer         M,
         integer         N,
         real            A[],
         integer         LDA,
         real            S[],
         real            U[],
         integer         LDU,
         real            VT[],
         integer         LDVT,
         real            WORK[],
         integer         LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgesvd)( const_cast<character*>(job_blas[JOBU]),
                         const_cast<character*>(job_blas[JOBVT]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgesvd)( const_cast<character*>(job_blas[JOBU]),
                         const_cast<character*>(job_blas[JOBVT]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(sgesvd)( const_cast<character*>(job_blas[JOBU]),
                        const_cast<character*>(job_blas[JOBVT]),
                        &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  gesvd( JobType const & JOBU,
         JobType const & JOBVT,
         integer         M,
         integer         N,
         doublereal      A[],
         integer         LDA,
         doublereal      S[],
         doublereal      U[],
         integer         LDU,
         doublereal      VT[],
         integer         LDVT,
         doublereal      WORK[],
         integer         LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgesvd)( const_cast<character*>(job_blas[JOBU]),
                         const_cast<character*>(job_blas[JOBVT]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgesvd)( const_cast<character*>(job_blas[JOBU]),
                         const_cast<character*>(job_blas[JOBVT]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(dgesvd)( const_cast<character*>(job_blas[JOBU]),
                        const_cast<character*>(job_blas[JOBVT]),
                        &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  /*
  //                       _     _
  //    __ _  ___  ___  __| | __| |
  //   / _` |/ _ \/ __|/ _` |/ _` |
  //  | (_| |  __/\__ \ (_| | (_| |
  //   \__, |\___||___/\__,_|\__,_|
  //   |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGESDD computes the singular value decomposition (SVD) of a real
   *  M-by-N matrix A, optionally computing the left and right singular
   *  vectors.  If singular vectors are desired, it uses a
   *  divide-and-conquer algorithm.
   *
   *  The SVD is written
   *
   *       A = U * SIGMA * transpose(V)
   *
   *  where SIGMA is an M-by-N matrix which is zero except for its
   *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   *  are the singular values of A; they are real and non-negative, and
   *  are returned in descending order.  The first min(m,n) columns of
   *  U and V are the left and right singular vectors of A.
   *
   *  Note that the routine returns VT = V**T, not V.
   *
   *  The divide and conquer algorithm makes very mild assumptions about
   *  floating point arithmetic. It will work on machines with a guard
   *  digit in add/subtract, or on those binary machines without guard
   *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
   *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
   *  without guard digits, but we know of none.
   *
   *  Arguments
   *  =========
   *
   *  JOBZ    (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix U:
   *          = 'A':  all M columns of U and all N rows of V**T are
   *                  returned in the arrays U and VT;
   *          = 'S':  the first min(M,N) columns of U and the first
   *                  min(M,N) rows of V**T are returned in the arrays U
   *                  and VT;
   *          = 'O':  If M >= N, the first N columns of U are overwritten
   *                  on the array A and all rows of V**T are returned in
   *                  the array VT;
   *                  otherwise, all columns of U are returned in the
   *                  array U and the first M rows of V**T are overwritten
   *                  in the array A;
   *          = 'N':  no columns of U or rows of V**T are computed.
   *
   *  M       (input) INTEGER
   *          The number of rows of the input matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the input matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit,
   *          if JOBZ = 'O',  A is overwritten with the first N columns
   *                          of U (the left singular vectors, stored
   *                          columnwise) if M >= N;
   *                          A is overwritten with the first M rows
   *                          of V**T (the right singular vectors, stored
   *                          rowwise) otherwise.
   *          if JOBZ .ne. 'O', the contents of A are destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A, sorted so that S(i) >= S(i+1).
   *
   *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
   *          UCOL = min(M,N) if JOBZ = 'S'.
   *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
   *          orthogonal matrix U;
   *          if JOBZ = 'S', U contains the first min(M,N) columns of U
   *          (the left singular vectors, stored columnwise);
   *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
   *
   *  LDU     (input) INTEGER
   *          The leading dimension of the array U.  LDU >= 1; if
   *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
   *
   *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
   *          N-by-N orthogonal matrix V**T;
   *          if JOBZ = 'S', VT contains the first min(M,N) rows of
   *          V**T (the right singular vectors, stored rowwise);
   *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
   *
   *  LDVT    (input) INTEGER
   *          The leading dimension of the array VT.  LDVT >= 1; if
   *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
   *          if JOBZ = 'S', LDVT >= min(M,N).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 1.
   *          If JOBZ = 'N',
   *            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
   *          If JOBZ = 'O',
   *            LWORK >= 3*min(M,N)*min(M,N) +
   *                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
   *          If JOBZ = 'S' or 'A'
   *            LWORK >= 3*min(M,N)*min(M,N) +
   *                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
   *          For good performance, LWORK should generally be larger.
   *          If LWORK = -1 but other input arguments are legal, WORK(1)
   *          returns the optimal LWORK.
   *
   *  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  DBDSDC did not converge, updating process failed.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Ming Gu and Huan Ren, Computer Science Division, University of
   *     California at Berkeley, USA
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(dgesdd)( character const   JOBZ[],
                        integer   const * M,
                        integer   const * N,
                        doublereal      * A,
                        integer   const * LDA,
                        doublereal      * S,
                        doublereal      * U,
                        integer   const * LDU,
                        doublereal      * VT,
                        integer   const * LDVT,
                        doublereal      * WORK,
                        integer   const * LWORK,
                        integer   const * IWORK, // (8*min(M,N))
                        integer         * info ) ;

    void
    LAPACKNAME(sgesdd)( character const   JOBZ[],
                        integer   const * M,
                        integer   const * N,
                        real            * A,
                        integer   const * LDA,
                        real            * S,
                        real            * U,
                        integer   const * LDU,
                        real            * VT,
                        integer   const * LDVT,
                        real            * WORK,
                        integer   const * LWORK,
                        integer   const * IWORK, // (8*min(M,N))
                        integer         * info  ) ;
  }
  #endif

  inline
  integer
  gesdd( JobType const & JOBZ,
         integer         M,
         integer         N,
         real            A[],
         integer         LDA,
         real            S[],
         real            U[],
         integer         LDU,
         real            VT[],
         integer         LDVT,
         real            WORK[],
         integer         LWORK,
         integer         IWORK[] ) // (8*min(M,N))
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgesdd)( const_cast<character*>(job_blas[JOBZ]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgesdd)( const_cast<character*>(job_blas[JOBZ]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #else
    LAPACKNAME(sgesdd)( const_cast<character*>(job_blas[JOBZ]),
                        &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  gesdd( JobType const & JOBZ,
         integer         M,
         integer         N,
         doublereal      A[],
         integer         LDA,
         doublereal      S[],
         doublereal      U[],
         integer         LDU,
         doublereal      VT[],
         integer         LDVT,
         doublereal      WORK[],
         integer         LWORK,
         integer         IWORK[] ) // (8*min(M,N))
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgesdd)( const_cast<character*>(job_blas[JOBZ]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgesdd)( const_cast<character*>(job_blas[JOBZ]),
                         &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #else
    LAPACKNAME(dgesdd)( const_cast<character*>(job_blas[JOBZ]),
                        &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, IWORK, &info ) ;
    #endif
    return info ;
  }


  /*
  //    ____ _____ _     ____  ____
  //   / ___| ____| |   / ___||  _ \
  //  | |  _|  _| | |   \___ \| | | |
  //  | |_| | |___| |___ ___) | |_| |
  //   \____|_____|_____|____/|____/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGELSD computes the minimum-norm solution to a real linear least
   *  squares problem:
   *      minimize 2-norm(| b - A*x |)
   *  using the singular value decomposition (SVD) of A. A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
   *  matrix X.
   *
   *  The problem is solved in three steps:
   *  (1) Reduce the coefficient matrix A to bidiagonal form with
   *      Householder transformations, reducing the original problem
   *      into a "bidiagonal least squares problem" (BLS)
   *  (2) Solve the BLS using a divide and conquer approach.
   *  (3) Apply back all the Householder tranformations to solve
   *      the original least squares problem.
   *
   *  The effective rank of A is determined by treating as zero those
   *  singular values which are less than RCOND times the largest singular
   *  value.
   *
   *  The divide and conquer algorithm makes very mild assumptions about
   *  floating point arithmetic. It will work on machines with a guard
   *  digit in add/subtract, or on those binary machines without guard
   *  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
   *  Cray-2. It could conceivably fail on hexadecimal or decimal machines
   *  without guard digits, but we know of none.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of A. N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrices B and X. NRHS >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, A has been destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, B is overwritten by the N-by-NRHS solution
   *          matrix X.  If m >= n and RANK = n, the residual
   *          sum-of-squares for the solution in the i-th column is given
   *          by the sum of squares of elements n+1:m in that column.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A in decreasing order.
   *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
   *
   *  RCOND   (input) DOUBLE PRECISION
   *          RCOND is used to determine the effective rank of A.
   *          Singular values S(i) <= RCOND*S(1) are treated as zero.
   *          If RCOND < 0, machine precision is used instead.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the number of singular values
   *          which are greater than RCOND*S(1).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK must be at least 1.
   *          The exact minimum amount of workspace needed depends on M,
   *          N and NRHS. As long as LWORK is at least
   *              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
   *          if M is greater than or equal to N or
   *              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
   *          if M is less than N, the code will execute correctly.
   *          SMLSIZ is returned by ILAENV and is equal to the maximum
   *          size of the subproblems at the bottom of the computation
   *          tree (usually about 25), and
   *             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *     Ming Gu and Ren-Cang Li, Computer Science Division, University of
   *       California at Berkeley, USA
   *     Osni Marques, LBNL/NERSC, USA
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

   void
   LAPACKNAME(sgelsd)( integer * m,
                       integer * n,
                       integer * nrhs,
                       real      a[],
	                     integer * lda,
                       real      b[],
                       integer * ldb,
                       real      s[],
                       real    * rcond,
                       integer * rank,
                       real      work[],
                       integer * lwork,
                       integer   iwork[],
                       integer * info ) ;

   void
   LAPACKNAME(dgelsd)( integer    * m,
                       integer    * n,
                       integer    * nrhs,
                       doublereal   a[],
	                     integer    * lda,
                       doublereal   b[],
                       integer    * ldb,
                       doublereal   s[],
                       doublereal * rcond,
                       integer    * rank,
                       doublereal   work[],
                       integer    * lwork,
                       integer      iwork[],
                       integer    * info ) ;
  }
  #endif

  inline
  integer
  gelsd( integer   m,
         integer   n,
         integer   nrhs,
         real      a[],
	       integer   lda,
         real      b[],
         integer   ldb,
         real      s[],
         real      rcond,
         integer & rank,
         real      work[],
         integer   lwork,
         integer   iwork[] )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info ) ;
    #else
    LAPACKNAME(sgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info) ;
    #endif
    return info ;
  }

  inline
  integer
  gelsd( integer    m,
         integer    n,
         integer    nrhs,
         doublereal a[],
	       integer    lda,
         doublereal b[],
         integer    ldb,
         doublereal s[],
         doublereal rcond,
         integer  & rank,
         doublereal work[],
         integer    lwork,
         integer    iwork[] )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info ) ;
    #else
    LAPACKNAME(dgelsd)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info) ;
    #endif
    return info ;
  }

  /*
  //    ____ _____ _     ____ ____
  //   / ___| ____| |   / ___/ ___|
  //  | |  _|  _| | |   \___ \___ \
  //  | |_| | |___| |___ ___) |__) |
  //   \____|_____|_____|____/____/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  DGELSS computes the minimum norm solution to a real linear least
   *  squares problem:
   *
   *  Minimize 2-norm(| b - A*x |).
   *
   *  using the singular value decomposition (SVD) of A. A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
   *  X.
   *
   *  The effective rank of A is determined by treating as zero those
   *  singular values which are less than RCOND times the largest singular
   *  value.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A. N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrices B and X. NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the first min(m,n) rows of A are overwritten with
   *          its right singular vectors, stored rowwise.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, B is overwritten by the N-by-NRHS solution
   *          matrix X.  If m >= n and RANK = n, the residual
   *          sum-of-squares for the solution in the i-th column is given
   *          by the sum of squares of elements n+1:m in that column.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,max(M,N)).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A in decreasing order.
   *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
   *
   *  RCOND   (input) DOUBLE PRECISION
   *          RCOND is used to determine the effective rank of A.
   *          Singular values S(i) <= RCOND*S(1) are treated as zero.
   *          If RCOND < 0, machine precision is used instead.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the number of singular values
   *          which are greater than RCOND*S(1).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 1, and also:
   *          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
   *          For good performance, LWORK should generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  the algorithm for computing the SVD failed to converge;
   *                if INFO = i, i off-diagonal elements of an intermediate
   *                bidiagonal form did not converge to zero.
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

   void
   LAPACKNAME(sgelss)( integer * m,
                       integer * n,
                       integer * nrhs,
                       real      a[],
	                     integer * lda,
                       real      b[],
                       integer * ldb,
                       real      s[],
                       real    * rcond,
                       integer * rank,
                       real    * work,
                       integer * lwork,
                       integer * info ) ;

   void
   LAPACKNAME(dgelss)( integer    * m,
                       integer    * n,
                       integer    * nrhs,
                       doublereal   a[],
	                     integer    * lda,
                       doublereal   b[],
                       integer    * ldb,
                       doublereal   s[],
                       doublereal * rcond,
                       integer    * rank,
                       doublereal * work,
                       integer    * lwork,
                       integer    * info ) ;
  }
  #endif

  inline
  integer
  gelss( integer   m,
         integer   n,
         integer   nrhs,
         real      a[],
	       integer   lda,
         real      b[],
         integer   ldb,
         real      s[],
         real      rcond,
         integer & rank,
         real      work[],
         integer   lwork )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgelss)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info ) ;
    #else
    LAPACKNAME(sgelss)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info) ;
    #endif
    return info ;
  }

  inline
  integer
  gelss( integer    m,
         integer    n,
         integer    nrhs,
         doublereal a[],
	       integer    lda,
         doublereal b[],
         integer    ldb,
         doublereal s[],
         doublereal rcond,
         integer  & rank,
         doublereal work[],
         integer    lwork )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgelss)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgelss)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info ) ;
    #else
    LAPACKNAME(dgelss)( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info) ;
    #endif
    return info ;
  }

  /*
  //              _
  //    __ _  ___| |___ _   _
  //   / _` |/ _ \ / __| | | |
  //  | (_| |  __/ \__ \ |_| |
  //   \__, |\___|_|___/\__, |
  //   |___/            |___/
  */
  /*\
   *  Purpose
   *  =======
   *
   *  SGELSY computes the minimum-norm solution to a real linear least
   *  squares problem:
   *      minimize || A * X - B ||
   *  using a complete orthogonal factorization of A.  A is an M-by-N
   *  matrix which may be rank-deficient.
   *
   *  Several right hand side vectors b and solution vectors x can be
   *  handled in a single call; they are stored as the columns of the
   *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
   *  matrix X.
   *
   *  The routine first computes a QR factorization with column pivoting:
   *      A * P = Q * [ R11 R12 ]
   *                  [  0  R22 ]
   *  with R11 defined as the largest leading submatrix whose estimated
   *  condition number is less than 1/RCOND.  The order of R11, RANK,
   *  is the effective rank of A.
   *
   *  Then, R22 is considered to be negligible, and R12 is annihilated
   *  by orthogonal transformations from the right, arriving at the
   *  complete orthogonal factorization:
   *     A * P = Q * [ T11 0 ] * Z
   *                 [  0  0 ]
   *  The minimum-norm solution is then
   *     X = P * Z' [ inv(T11)*Q1'*B ]
   *                [        0       ]
   *  where Q1 consists of the first RANK columns of Q.
   *
   *  This routine is basically identical to the original xGELSX except
   *  three differences:
   *    o The call to the subroutine xGEQPF has been substituted by the
   *      the call to the subroutine xGEQP3. This subroutine is a Blas-3
   *      version of the QR factorization with column pivoting.
   *    o Matrix B (the right hand side) is updated with Blas-3.
   *    o The permutation of matrix B (the right hand side) is faster and
   *      more simple.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of
   *          columns of matrices B and X. NRHS >= 0.
   *
   *  A       (input/output) REAL array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, A has been overwritten by details of its
   *          complete orthogonal factorization.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  B       (input/output) REAL array, dimension (LDB,NRHS)
   *          On entry, the M-by-NRHS right hand side matrix B.
   *          On exit, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B. LDB >= max(1,M,N).
   *
   *  JPVT    (input/output) INTEGER array, dimension (N)
   *          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
   *          to the front of AP, otherwise column i is a free column.
   *          On exit, if JPVT(i) = k, then the i-th column of AP
   *          was the k-th column of A.
   *
   *  RCOND   (input) REAL
   *          RCOND is used to determine the effective rank of A, which
   *          is defined as the order of the largest leading triangular
   *          submatrix R11 in the QR factorization with pivoting of A,
   *          whose estimated condition number < 1/RCOND.
   *
   *  RANK    (output) INTEGER
   *          The effective rank of A, i.e., the order of the submatrix
   *          R11.  This is the same as the order of the submatrix T11
   *          in the complete orthogonal factorization of A.
   *
   *  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          The unblocked strategy requires that:
   *             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
   *          where MN = min( M, N ).
   *          The block algorithm requires that:
   *             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
   *          where NB is an upper bound on the blocksize returned
   *          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR,
   *          and SORMRZ.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: If INFO = -i, the i-th argument had an illegal value.
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

   void
   LAPACKNAME(sgelsy)( integer * m,
                       integer * n,
                       integer * nrhs,
                       real      a[],
	                     integer * lda,
                       real      b[],
                       integer * ldb,
                       integer   jpvt[],
                       real    * rcond,
                       integer * rank,
                       real    * work,
                       integer * lwork,
                       integer * info ) ;

   void
   LAPACKNAME(dgelsy)( integer    * m,
                       integer    * n,
                       integer    * nrhs,
                       doublereal   a[],
	                     integer    * lda,
                       doublereal   b[],
                       integer    * ldb,
                       integer      jpvt[],
                       doublereal * rcond,
                       integer    * rank,
                       doublereal * work,
                       integer    * lwork,
                       integer    * info ) ;
  }
  #endif

  inline
  integer
  gelsy( integer   m,
         integer   n,
         integer   nrhs,
         real      a[],
	       integer   lda,
         real      b[],
         integer   ldb,
         integer   jpvt[],
         real      rcond,
         integer & rank,
         real      work[],
         integer   lwork )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgelsy)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgelsy)( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info ) ;
    #else
    LAPACKNAME(sgelsy)( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info) ;
    #endif
    return info ;
  }

  inline
  integer
  gelsy( integer    m,
         integer    n,
         integer    nrhs,
         doublereal a[],
	       integer    lda,
         doublereal b[],
         integer    ldb,
         integer    jpvt[],
         doublereal rcond,
         integer  & rank,
         doublereal work[],
         integer    lwork )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgelsy)( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgelsy)( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info ) ;
    #else
    LAPACKNAME(dgelsy)( &m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &rcond, &rank, work, &lwork, &info) ;
    #endif
    return info ;
  }

  /*
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

   void
   LAPACKNAME(sormqr)( character * SIDE,
                       character * TRANS,
                       integer   * M,
                       integer   * N,
                       integer   * K,
                       real        A[],
                       integer   * LDA,
                       real        TAU[],
                       real        C[],
                       integer   * LDC,
                       real        WORK[],
                       integer   * LWORK,
                       integer   * INFO ) ;

   void
   LAPACKNAME(dormqr)( character * SIDE,
                       character * TRANS,
                       integer   * M,
                       integer   * N,
                       integer   * K,
                       doublereal  A[],
                       integer   * LDA,
                       doublereal  TAU[],
                       doublereal  C[],
                       integer   * LDC,
                       doublereal  WORK[],
                       integer   * LWORK,
                       integer   * INFO ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DORMQR overwrites the general real M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'T':      Q**T * C       C * Q**T
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q**T from the Left;
   *          = 'R': apply Q or Q**T from the Right.
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N':  No transpose, apply Q;
   *          = 'T':  Transpose, apply Q**T.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
   *          The i-th column must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DGEQRF in the first k columns of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *          If SIDE = 'L', LDA >= max(1,M);
   *          if SIDE = 'R', LDA >= max(1,N).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DGEQRF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the M-by-N matrix C.
   *          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If SIDE = 'L', LWORK >= max(1,N);
   *          if SIDE = 'R', LWORK >= max(1,M).
   *          For optimum performance LWORK >= N*NB if SIDE = 'L', and
   *          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
   *          blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
  \*/

  inline
  integer
  ormqr( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         real                  A[],
         integer               LDA,
         real                  TAU[],
         real                  C[],
         integer               LDC,
         real                  WORK[],
         integer               LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sormqr)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sormqr)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(sormqr)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  ormqr( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         doublereal            A[],
         integer               LDA,
         doublereal            TAU[],
         doublereal            C[],
         integer               LDC,
         doublereal            WORK[],
         integer               LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dormqr)( const_cast<character*>(side_cblas[SIDE]),
                         const_cast<character*>(trans_cblas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dormqr)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(dormqr)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  #ifndef ALGLIN_USE_ACCELERATE
  // use standard Lapack routine
  extern "C" {

   void
   LAPACKNAME(sorm2r)( character * SIDE,
                       character * TRANS,
                       integer   * M,
                       integer   * N,
                       integer   * K,
                       real        A[],
                       integer   * LDA,
                       real        TAU[],
                       real        C[],
                       integer   * LDC,
                       real        WORK[],
                       integer   * INFO ) ;

   void
   LAPACKNAME(dorm2r)( character * SIDE,
                       character * TRANS,
                       integer   * M,
                       integer   * N,
                       integer   * K,
                       doublereal  A[],
                       integer   * LDA,
                       doublereal  TAU[],
                       doublereal  C[],
                       integer   * LDC,
                       doublereal  WORK[],
                       integer   * INFO ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DORM2R overwrites the general real m by n matrix C with
   *
   *        Q * C  if SIDE = 'L' and TRANS = 'N', or
   *
   *        Q'* C  if SIDE = 'L' and TRANS = 'T', or
   *
   *        C * Q  if SIDE = 'R' and TRANS = 'N', or
   *
   *        C * Q' if SIDE = 'R' and TRANS = 'T',
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q' from the Left
   *          = 'R': apply Q or Q' from the Right
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N': apply Q  (No transpose)
   *          = 'T': apply Q' (Transpose)
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
   *          The i-th column must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DGEQRF in the first k columns of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.
   *          If SIDE = 'L', LDA >= max(1,M);
   *          if SIDE = 'R', LDA >= max(1,N).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DGEQRF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the m by n matrix C.
   *          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace) DOUBLE PRECISION array, dimension
   *                                   (N) if SIDE = 'L',
   *                                   (M) if SIDE = 'R'
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit
   *          < 0: if INFO = -i, the i-th argument had an illegal value
   *
  \*/

  inline
  integer
  orm2r( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         real                  A[],
         integer               LDA,
         real                  TAU[],
         real                  C[],
         integer               LDC,
         real                  WORK[] )
  { integer info = 0;
    #ifdef ALGLIN_USE_ACCELERATE
    CLAPACKNAME(sorm2r)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &info ) ;
    #else
    LAPACKNAME(sorm2r)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  orm2r( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         doublereal            A[],
         integer               LDA,
         doublereal            TAU[],
         doublereal            C[],
         integer               LDC,
         doublereal            WORK[] )
  { integer info = 0;
    #ifdef ALGLIN_USE_ACCELERATE
    CLAPACKNAME(dorm2r)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &info ) ;
    #else
    LAPACKNAME(dorm2r)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &info ) ;
    #endif
    return info ;
  }

  /*
  //   _             __ _
  //  | | __ _ _ __ / _| |_
  //  | |/ _` | '__| |_| __|
  //  | | (_| | |  |  _| |_
  //  |_|\__,_|_|  |_|  \__|
  */
  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(slarft)( character   DIRECT[],
                        character   STOREV[],
                        integer   * N,
                        integer   * K,
                        real        V[],
                        integer   * LDV,
                        real        TAU[],
                        real        T[],
                        integer   * LDT ) ;

    void
    LAPACKNAME(dlarft)( character    DIRECT[],
                        character    STOREV[],
                        integer    * N,
                        integer    * K,
                        doublereal   V[],
                        integer    * LDV,
                        doublereal   TAU[],
                        doublereal   T[],
                        integer    * LDT ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  forms the triangular factor T of a real block reflector H
   *  of order n, which is defined as a product of k elementary reflectors.
   *
   *  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
   *
   *  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
   *
   *  If STOREV = 'C', the vector which defines the elementary reflector
   *  H(i) is stored in the i-th column of the array V, and
   *
   *     H  =  I - V * T * V'
   *
   *  If STOREV = 'R', the vector which defines the elementary reflector
   *  H(i) is stored in the i-th row of the array V, and
   *
   *     H  =  I - V' * T * V
   *
   *  Arguments
   *  =========
   *
   *  DIRECT = Specifies the order in which the elementary reflectors are
   *           multiplied to form the block reflector:
   *         = 'F': H = H(1) H(2) . . . H(k) (Forward)
   *         = 'B': H = H(k) . . . H(2) H(1) (Backward)
   *
   *  STOREV = Specifies how the vectors which define the elementary
   *           reflectors are stored (see also Further Details):
   *         = 'C': columnwise
   *         = 'R': rowwise
   *
   *  N = The order of the block reflector H. N >= 0.
   *  K = The order of the triangular factor T (= the number of elementary reflectors). K >= 1.
   *  V = (LDV,K) if STOREV = 'C'
   *      (LDV,N) if STOREV = 'R'
   *      The matrix V. See further details.
   *
   *  LDV = The leading dimension of the array V.
   *        If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
   *
   *  TAU = TAU(i) must contain the scalar factor of the elementary reflector H(i).
   *
   *  T = The k by k triangular factor T of the block reflector.
   *      If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
   *      lower triangular. The rest of the array is not used.
   *
   *  LDT = The leading dimension of the array T. LDT >= K.
   *
   *  Further Details
   *  ===============
   *
   *  The shape of the matrix V and the storage of the vectors which define
   *  the H(i) is best illustrated by the following example with n = 5 and
   *  k = 3. The elements equal to 1 are not stored; the corresponding
   *  array elements are modified but restored on exit. The rest of the
   *  array is not used.
   *
   *  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
   *
   *               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
   *                   ( v1  1    )                     (     1 v2 v2 v2 )
   *                   ( v1 v2  1 )                     (        1 v3 v3 )
   *                   ( v1 v2 v3 )
   *                   ( v1 v2 v3 )
   *
   *  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
   *
   *               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
   *                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
   *                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
   *                   (     1 v3 )
   *                   (        1 )
   *
  \*/

  inline
  void
  larft( DirectionType & DIRECT, // direct_blas[DIRECT]
         StorageType   & STOREV, // store_blas[STOREV]
         integer         N,
         integer         K,
         real            V[],
         integer         LDV,
         real            TAU[],
         real            T[],
         integer         LDT )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(slarft)( const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(slarft)( const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #else
  { LAPACKNAME(slarft)( const_cast<character*>(direct_blas[DIRECT]),
                        const_cast<character*>(store_blas[STOREV]),
                        &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #endif

  inline
  void
  larft( DirectionType & DIRECT, // direct_blas[DIRECT]
         StorageType   & STOREV, // store_blas[STOREV]
         integer         N,
         integer         K,
         doublereal      V[],
         integer         LDV,
         doublereal      TAU[],
         doublereal      T[],
         integer         LDT )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(dlarft)( const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(dlarft)( const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #else
  { LAPACKNAME(dlarft)( const_cast<character*>(direct_blas[DIRECT]),
                        const_cast<character*>(store_blas[STOREV]),
                        &N, &K, V, &LDV, TAU, T, &LDT ) ; }
  #endif

  /*
  //   _             __
  //  | | __ _ _ __ / _| __ _
  //  | |/ _` | '__| |_ / _` |
  //  | | (_| | |  |  _| (_| |
  //  |_|\__,_|_|  |_|  \__, |
  //                    |___/
  */
  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(slarfg)( integer * N,
                        real    * ALPHA,
                        real      X[],
                        integer * INCX,
                        real      TAU[] ) ;

    void
    LAPACKNAME(dlarfg)( integer    * N,
                        doublereal * ALPHA,
                        doublereal   X[],
                        integer    * INCX,
                        doublereal   TAU[] ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  generates a real elementary reflector H of order n, such that
   *
   *        H * ( alpha ) = ( beta ),   H' * H = I.
   *            (   x   )   (   0  )
   *
   *  where alpha and beta are scalars, and x is an (n-1)-element real
   *  vector. H is represented in the form
   *
   *        H = I - tau * ( 1 ) * ( 1 v' ) ,
   *                      ( v )
   *
   *  where tau is a real scalar and v is a real (n-1)-element vector.
   *
   *  If the elements of x are all zero, then tau = 0 and H is taken to be
   *  the unit matrix.
   *
   *  Otherwise  1 <= tau <= 2.
   *
   *  Arguments
   *  =========
   *
   *  N     = The order of the elementary reflector.
   *  ALPHA = On entry, the value alpha.
   *          On exit, it is overwritten with the value beta.
   *  X     = On entry, the vector x.
   *          On exit, it is overwritten with the vector v.
   *  INCX = The increment between elements of X. INCX <> 0.
   *  TAU  = The value tau.
  \*/

  inline
  void
  larfg( integer N,
         real  & ALPHA,
         real    X[],
         integer INCX,
         real    TAU[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(slarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(slarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #else
  { LAPACKNAME(slarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #endif

  inline
  void
  larfg( integer      N,
         doublereal & ALPHA,
         doublereal   X[],
         integer      INCX,
         doublereal   TAU[] )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(dlarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(dlarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #else
  { LAPACKNAME(dlarfg)( &N, &ALPHA, X, &INCX, TAU ) ; }
  #endif

  /*
  //   _             __ _
  //  | | __ _ _ __ / _| |__
  //  | |/ _` | '__| |_| '_ \
  //  | | (_| | |  |  _| |_) |
  //  |_|\__,_|_|  |_| |_.__/
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(slarfb)( character   SIDE[],
                        character   TRANS[],
                        character   DIRECT[],
                        character   STOREV[],
                        integer   * M,
                        integer   * N,
                        integer   * K,
                        real        V[],
                        integer   * LDV,
                        real        T[],
                        integer   * LDT,
                        real        C[],
                        integer   * LDC,
                        real        WORK[],
                        integer   * LDWORK ) ;

    void
    LAPACKNAME(dlarfb)( character    SIDE[],
                        character    TRANS[],
                        character    DIRECT[],
                        character    STOREV[],
                        integer    * M,
                        integer    * N,
                        integer    * K,
                        doublereal   V[],
                        integer    * LDV,
                        doublereal   T[],
                        integer    * LDT,
                        doublereal   C[],
                        integer    * LDC,
                        doublereal   WORK[],
                        integer    * LDWORK ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  applies a real block reflector H or its transpose H' to a
   *  real m by n matrix C, from either the left or the right.
   *
   *  Arguments
   *  =========
   *
   *  SIDE = 'L': apply H or H' from the Left
   *       = 'R': apply H or H' from the Right
   *
   *  TRANS = 'N': apply H (No transpose)
   *        = 'T': apply H' (Transpose)
   *
   *  DIRECT Indicates how H is formed from a product of elementary reflectors
   *       = 'F': H = H(1) H(2) . . . H(k) (Forward)
   *       = 'B': H = H(k) . . . H(2) H(1) (Backward)
   *
   *  STOREV Indicates how the vectors which define the elementary reflectors are stored:
   *       = 'C': Columnwise
   *       = 'R': Rowwise
   *
   *  M = The number of rows of the matrix C.
   *  N = The number of columns of the matrix C.
   *  K = The order of the matrix T (= the number of elementary
   *      reflectors whose product defines the block reflector).
   *  V = (LDV,K) if STOREV = 'C'
   *      (LDV,M) if STOREV = 'R' and SIDE = 'L'
   *      (LDV,N) if STOREV = 'R' and SIDE = 'R'
   *      The matrix V. See further details.
   *
   *  LDV = The leading dimension of the array V.
   *          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
   *          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
   *          if STOREV = 'R', LDV >= K.
   *
   *  T = (LDT,K) The triangular k by k matrix T in the representation of the block reflector.
   *  LDT = The leading dimension of the array T. LDT >= K.
   *  C  = (LDC,N)
   *     On entry, the m by n matrix C.
   *     On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
   *  LDC = The leading dimension of the array C. LDA >= max(1,M).
   *  WORK = (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
   *  LDWORK = The leading dimension of the array WORK.
   *     If SIDE = 'L', LDWORK >= max(1,N);
   *     if SIDE = 'R', LDWORK >= max(1,M).
   *
   *  =====================================================================
  \*/

  inline
  void
  larfb( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         DirectionType const & DIRECT,
         StorageType   const & STOREV,
         integer               M,
         integer               N,
         integer               K,
         real          const   V[],
         integer               LDV,
         real          const   T[],
         integer               LDT,
         real                  C[],
         integer               LDC,
         real                  WORK[],
         integer               LDWORK )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(slarfb)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &M, &N, &K, V, &LDV, T, &LDT, C, &LDC, WORK, &LDWORK ) ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(slarfb)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &M, &N, &K,
                         const_cast<real*>(V), &LDV,
                         const_cast<real*>(T), &LDT,
                         C, &LDC, WORK, &LDWORK ) ;
  }
  #else
  { LAPACKNAME(slarfb)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        const_cast<character*>(direct_blas[DIRECT]),
                        const_cast<character*>(store_blas[STOREV]),
                        &M, &N, &K,
                        const_cast<real*>(V), &LDV,
                        const_cast<real*>(T), &LDT,
                        C, &LDC, WORK, &LDWORK ) ;
  }
  #endif

  inline
  void
  larfb( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         DirectionType const & DIRECT,
         StorageType   const & STOREV,
         integer               M,
         integer               N,
         integer               K,
         doublereal    const   V[],
         integer               LDV,
         doublereal    const   T[],
         integer               LDT,
         doublereal            C[],
         integer               LDC,
         doublereal            WORK[],
         integer               LDWORK )
  #if defined(ALGLIN_USE_OPENBLAS)
  { LAPACK_NAME(dlarfb)(const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        const_cast<character*>(direct_blas[DIRECT]),
                        const_cast<character*>(store_blas[STOREV]),
                        &M, &N, &K, V, &LDV, T, &LDT, C, &LDC, WORK, &LDWORK ) ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { CLAPACKNAME(dlarfb)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         const_cast<character*>(direct_blas[DIRECT]),
                         const_cast<character*>(store_blas[STOREV]),
                         &M, &N, &K,
                         const_cast<doublereal*>(V), &LDV,
                         const_cast<doublereal*>(T), &LDT,
                         C, &LDC, WORK, &LDWORK ) ;
  }
  #else
  { LAPACKNAME(dlarfb)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        const_cast<character*>(direct_blas[DIRECT]),
                        const_cast<character*>(store_blas[STOREV]),
                        &M, &N, &K,
                        const_cast<doublereal*>(V), &LDV,
                        const_cast<doublereal*>(T), &LDT,
                        C, &LDC, WORK, &LDWORK ) ;
  }
  #endif

  /*
  //                          __
  //    __ _  ___  __ _ _ __ / _|
  //   / _` |/ _ \/ _` | '__| |_
  //  | (_| |  __/ (_| | |  |  _|
  //   \__, |\___|\__, |_|  |_|
  //   |___/         |_|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgeqrf)( integer  * M,
                        integer  * N,
                        real       A[],
                        integer  * LDA,
                        real       TAU[],
                        real       WORK[],
                        integer  * LWORK,
                        integer  * INFO ) ;

    void
    LAPACKNAME(dgeqrf)( integer    * M,
                        integer    * N,
                        doublereal   A[],
                        integer    * LDA,
                        doublereal   TAU[],
                        doublereal   WORK[],
                        integer    * LWORK,
                        integer    * INFO ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEQRF computes a QR factorization of a real M-by-N matrix A:
   *  A = Q * R.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the elements on and above the diagonal of the array
   *          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
   *          upper triangular if m >= n); the elements below the diagonal,
   *          with the array TAU, represent the orthogonal matrix Q as a
   *          product of min(m,n) elementary reflectors (see Further
   *          Details).
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The scalar factors of the elementary reflectors (see Further
   *          Details).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,N).
   *          For optimum performance LWORK >= N*NB, where NB is
   *          the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real scalar, and v is a real vector with
   *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
   *  and tau in TAU(i).
   *
   *  =====================================================================
  \*/

  inline
  integer
  geqrf( integer M,
         integer N,
         real    A[],
         integer LDA,
         real    TAU[],
         real    WORK[],
         integer LWORK )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(sgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer info = 0;
    CLAPACKNAME(sgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #else
  { integer info = 0;
    LAPACKNAME(sgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #endif

  inline
  integer
  geqrf( integer    M,
         integer    N,
         doublereal A[],
         integer    LDA,
         doublereal TAU[],
         doublereal WORK[],
         integer    LWORK )
  #if defined(ALGLIN_USE_OPENBLAS)
  { integer info = 0 ;
    LAPACK_NAME(dgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #elif defined(ALGLIN_USE_ACCELERATE)
  { integer info = 0;
    CLAPACKNAME(dgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #else
  { integer info = 0;
    LAPACKNAME(dgeqrf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    return info ;
  }
  #endif

  /*
  //                        ____
  //    __ _  ___  __ _ _ _|___ \
  //   / _` |/ _ \/ _` | '__|__) |
  //  | (_| |  __/ (_| | |  / __/
  //   \__, |\___|\__, |_| |_____|
  //   |___/         |_|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgeqr2)( integer * M,
                        integer * N,
                        real      A[],
                        integer * LDA,
                        real      TAU[],
                        real      WORK[],
                        integer * INFO ) ;

    void
    LAPACKNAME(dgeqr2)( integer    * M,
                        integer    * N,
                        doublereal   A[],
                        integer    * LDA,
                        doublereal   TAU[],
                        doublereal   WORK[],
                        integer    * INFO ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  computes a QR factorization of a real m by n matrix A:
   *  A = Q * R.
   *
   *  Arguments
   *  =========
   *
   *  M = The number of rows of the matrix A.  M >= 0.
   *  N = The number of columns of the matrix A.  N >= 0.
   *  A = (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the m by n matrix A.
   *          On exit, the elements on and above the diagonal of the array
   *          contain the min(m,n) by n upper trapezoidal matrix R (R is
   *          upper triangular if m >= n); the elements below the diagonal,
   *          with the array TAU, represent the orthogonal matrix Q as a
   *          product of elementary reflectors (see Further Details).
   *
   *  LDA = The leading dimension of the array A.  LDA >= max(1,M).
   *  TAU = The scalar factors of the elementary reflectors (see Further Details).
   *  WORK = (workspace) DOUBLE PRECISION array, dimension (N)
   *  INFO = 0: successful exit
   *       < 0: if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real scalar, and v is a real vector with
   *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
   *  and tau in TAU(i).
  \*/

  inline
  integer
  geqr2( integer M,
         integer N,
         real    A[],
         integer LDA,
         real    TAU[],
         real    WORK[] )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #else
    LAPACKNAME(sgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  geqr2( integer    M,
         integer    N,
         doublereal A[],
         integer    LDA,
         doublereal TAU[],
         doublereal WORK[] )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #else
    LAPACKNAME(dgeqr2)( &M, &N, A, &LDA, TAU, WORK, &info ) ;
    #endif
    return info ;
  }

  /*
  //                         _
  //    __ _  ___  __ _ _ __| |_
  //   / _` |/ _ \/ _` | '__| __|
  //  | (_| |  __/ (_| | |  | |_
  //   \__, |\___|\__, |_|   \__|
  //   |___/         |_|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  extern "C" {

    void
    LAPACKNAME(sgeqrt)( integer * M,
                        integer * N,
                        integer * NB,
                        real      A[],
                        integer * LDA,
                        real      T[],
                        integer * LDT,
                        real      WORK[],
                        integer * INFO ) ;

    void
    LAPACKNAME(dgeqrt)( integer    * M,
                        integer    * N,
                        integer    * NB,
                        doublereal   A[],
                        integer    * LDA,
                        doublereal   T[],
                        integer    * LDT,
                        doublereal   WORK[],
                        integer    * INFO ) ;
  }
  #endif

  /*
  //                         _____
  //    __ _  ___  __ _ _ __|___ /
  //   / _` |/ _ \/ _` | '_ \ |_ \
  //  | (_| |  __/ (_| | |_) |__) |
  //   \__, |\___|\__, | .__/____/
  //   |___/         |_|_|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(sgeqp3)( integer   * M,
                        integer   * N,
                        real        A[],
                        integer   * LDA,
                        integer     JPVT[],
                        real        TAU[],
                        real        WORK[],
                        integer   * LWORK,
                        integer   * INFO ) ;

    void
    LAPACKNAME(dgeqp3)( integer   * M,
                        integer   * N,
                        doublereal  A[],
                        integer   * LDA,
                        integer     JPVT[],
                        doublereal  TAU[],
                        doublereal  WORK[],
                        integer   * LWORK,
                        integer   * INFO ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DGEQP3 computes a QR factorization with column pivoting of a
   *  matrix A:  A*P = Q*R  using Level 3 BLAS.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit, the upper triangle of the array contains the
   *          min(M,N)-by-N upper trapezoidal matrix R; the elements below
   *          the diagonal, together with the array TAU, represent the
   *          orthogonal matrix Q as a product of min(M,N) elementary
   *          reflectors.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A. LDA >= max(1,M).
   *
   *  JPVT    (input/output) INTEGER array, dimension (N)
   *          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
   *          to the front of A*P (a leading column); if JPVT(J)=0,
   *          the J-th column of A is a free column.
   *          On exit, if JPVT(J)=K, then the J-th column of A*P was the
   *          the K-th column of A.
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The scalar factors of the elementary reflectors.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 3*N+1.
   *          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
   *          is the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0: successful exit.
   *          < 0: if INFO = -i, the i-th argument had an illegal value.
   *
   *  Further Details
   *  ===============
   *
   *  The matrix Q is represented as a product of elementary reflectors
   *
   *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
   *
   *  Each H(i) has the form
   *
   *     H(i) = I - tau * v * v'
   *
   *  where tau is a real/complex scalar, and v is a real/complex vector
   *  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
   *  A(i+1:m,i), and tau in TAU(i).
   *
   *  Based on contributions by
   *    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
   *    X. Sun, Computer Science Dept., Duke University, USA
   *
  \*/

  inline
  integer
  geqp3( integer   M,
         integer   N,
         real      A[],
         integer   LDA,
         integer   JPVT[],
         real      TAU[],
         real      WORK[],
         integer   LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(sgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  geqp3( integer    M,
         integer    N,
         doublereal A[],
         integer    LDA,
         integer    JPVT[],
         doublereal TAU[],
         doublereal WORK[],
         integer    LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(dgeqp3)( &M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  /*
  //   _                 __
  //  | |_ _____ __ ____/ _|
  //  | __|_  / '__|_  / |_
  //  | |_ / /| |   / /|  _|
  //   \__/___|_|  /___|_|
  */

  /*\
   *  Purpose
   *  =======
   *
   *  DTZRZF reduces the M-by-N ( M<=N ) real upper trapezoidal matrix A
   *  to upper triangular form by means of orthogonal transformations.
   *
   *  The upper trapezoidal matrix A is factored as
   *
   *     A = ( R  0 ) * Z,
   *
   *  where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
   *  triangular matrix.
   *
   *  Arguments
   *  =========
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix A.  N >= M.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the leading M-by-N upper trapezoidal part of the
   *          array A must contain the matrix to be factorized.
   *          On exit, the leading M-by-M upper triangular part of A
   *          contains the upper triangular matrix R, and elements M+1 to
   *          N of the first M rows of A, with the array TAU, represent the
   *          orthogonal matrix Z as a product of M elementary reflectors.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  TAU     (output) DOUBLE PRECISION array, dimension (M)
   *          The scalar factors of the elementary reflectors.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,M).
   *          For optimum performance LWORK >= M*NB, where NB is
   *          the optimal blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
   *
   *  The factorization is obtained by Householder's method.  The kth
   *  transformation matrix, Z( k ), which is used to introduce zeros into
   *  the ( m - k + 1 )th row of A, is given in the form
   *
   *     Z( k ) = ( I     0   ),
   *              ( 0  T( k ) )
   *
   *  where
   *
   *     T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
   *                                                 (   0    )
   *                                                 ( z( k ) )
   *
   *  tau is a scalar and z( k ) is an ( n - m ) element vector.
   *  tau and z( k ) are chosen to annihilate the elements of the kth row
   *  of X.
   *
   *  The scalar tau is returned in the kth element of TAU and the vector
   *  u( k ) in the kth row of A, such that the elements of z( k ) are
   *  in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned in
   *  the upper triangular part of A.
   *
   *  Z is given by
   *
   *     Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
   *
   *  =====================================================================
  \*/

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(stzrzf)( integer * M,
                        integer * N,
                        real    * A,
                        integer * LDA,
                        real    * TAU,
                        real    * WORK,
                        integer * LWORK,
                        integer * INFO ) ;

    void
    LAPACKNAME(dtzrzf)( integer    * M,
                        integer    * N,
                        doublereal * A,
                        integer    * LDA,
                        doublereal * TAU,
                        doublereal * WORK,
                        integer    * LWORK,
                        integer    * INFO ) ;

  }
  #endif

  inline
  integer
  tzrzf( integer M,
         integer N,
         real    A[],
         integer LDA,
         real    TAU[],
         real    WORK[],
         integer LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(stzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(stzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(stzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  tzrzf( integer    M,
         integer    N,
         doublereal A[],
         integer    LDA,
         doublereal TAU[],
         doublereal WORK[],
         integer    LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dtzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dtzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(dtzrzf)( &M, &N, A, &LDA, TAU, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  /*
  //    ___  _ __ _ __ ___  _ __ ____
  //   / _ \| '__| '_ ` _ \| '__|_  /
  //  | (_) | |  | | | | | | |   / /
  //   \___/|_|  |_| |_| |_|_|  /___|
  */

  #if defined(ALGLIN_USE_ATLAS) || defined(ALGLIN_USE_LAPACK)
  // use standard Lapack routine
  extern "C" {

    void
    LAPACKNAME(sormrz)( character * SIDE,
                        character * TRANS,
                        integer   * M,
                        integer   * N,
                        integer   * K,
                        integer   * L,
                        real        A[],
                        integer   * LDA,
                        real        TAU[],
                        real        C[],
                        integer   * LDC,
                        real        WORK[],
                        integer   * LWORK,
                        integer   * INFO  ) ;

    void
    LAPACKNAME(dormrz)( character  * SIDE,
                        character  * TRANS,
                        integer    * M,
                        integer    * N,
                        integer    * K,
                        integer    * L,
                        doublereal A[],
                        integer    * LDA,
                        doublereal   TAU[],
                        doublereal   C[],
                        integer    * LDC,
                        doublereal   WORK[],
                        integer    * LWORK,
                        integer    * INFO ) ;

  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DORMRZ overwrites the general real M-by-N matrix C with
   *
   *                  SIDE = 'L'     SIDE = 'R'
   *  TRANS = 'N':      Q * C          C * Q
   *  TRANS = 'T':      Q**T * C       C * Q**T
   *
   *  where Q is a real orthogonal matrix defined as the product of k
   *  elementary reflectors
   *
   *        Q = H(1) H(2) . . . H(k)
   *
   *  as returned by DTZRZF. Q is of order M if SIDE = 'L' and of order N
   *  if SIDE = 'R'.
   *
   *  Arguments
   *  =========
   *
   *  SIDE    (input) CHARACTER*1
   *          = 'L': apply Q or Q**T from the Left;
   *          = 'R': apply Q or Q**T from the Right.
   *
   *  TRANS   (input) CHARACTER*1
   *          = 'N':  No transpose, apply Q;
   *          = 'T':  Transpose, apply Q**T.
   *
   *  M       (input) INTEGER
   *          The number of rows of the matrix C. M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the matrix C. N >= 0.
   *
   *  K       (input) INTEGER
   *          The number of elementary reflectors whose product defines
   *          the matrix Q.
   *          If SIDE = 'L', M >= K >= 0;
   *          if SIDE = 'R', N >= K >= 0.
   *
   *  L       (input) INTEGER
   *          The number of columns of the matrix A containing
   *          the meaningful part of the Householder reflectors.
   *          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
   *
   *  A       (input) DOUBLE PRECISION array, dimension
   *                               (LDA,M) if SIDE = 'L',
   *                               (LDA,N) if SIDE = 'R'
   *          The i-th row must contain the vector which defines the
   *          elementary reflector H(i), for i = 1,2,...,k, as returned by
   *          DTZRZF in the last k rows of its array argument A.
   *          A is modified by the routine but restored on exit.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A. LDA >= max(1,K).
   *
   *  TAU     (input) DOUBLE PRECISION array, dimension (K)
   *          TAU(i) must contain the scalar factor of the elementary
   *          reflector H(i), as returned by DTZRZF.
   *
   *  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
   *          On entry, the M-by-N matrix C.
   *          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
   *
   *  LDC     (input) INTEGER
   *          The leading dimension of the array C. LDC >= max(1,M).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If SIDE = 'L', LWORK >= max(1,N);
   *          if SIDE = 'R', LWORK >= max(1,M).
   *          For optimum performance LWORK >= N*NB if SIDE = 'L', and
   *          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
   *          blocksize.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *
   *  Further Details
   *  ===============
   *
   *  Based on contributions by
   *    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
   *
   *  =====================================================================
  \*/

  inline
  integer
  ormrz( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         integer               L,
         real                  A[],
         integer               LDA,
         real                  TAU[],
         real                  C[],
         integer               LDC,
         real                  WORK[],
         integer               LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(sormrz)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(sormrz)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(sormrz)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  inline
  integer
  ormrz( SideMultiply  const & SIDE,
         Transposition const & TRANS,
         integer               M,
         integer               N,
         integer               K,
         integer               L,
         doublereal            A[],
         integer               LDA,
         doublereal            TAU[],
         doublereal            C[],
         integer               LDC,
         doublereal            WORK[],
         integer               LWORK )
  { integer info = 0 ;
    #if defined(ALGLIN_USE_OPENBLAS)
    LAPACK_NAME(dormrz)( const_cast<character*>(side_cblas[SIDE]),
                         const_cast<character*>(trans_cblas[TRANS]),
                         &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info);
    #elif defined(ALGLIN_USE_ACCELERATE)
    CLAPACKNAME(dormrz)( const_cast<character*>(side_blas[SIDE]),
                         const_cast<character*>(trans_blas[TRANS]),
                         &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #else
    LAPACKNAME(dormrz)( const_cast<character*>(side_blas[SIDE]),
                        const_cast<character*>(trans_blas[TRANS]),
                        &M, &N, &K, &L, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &info ) ;
    #endif
    return info ;
  }

  /*\
   *  _        _    ___ ____ _
   * | |      / \  |_ _/ ___/ |
   * | |     / _ \  | | |   | |
   * | |___ / ___ \ | | |___| |
   * |_____/_/   \_\___\____|_|
   *
  \*/

  #ifndef ALGLIN_USE_ACCELERATE
  extern "C" {

    void
    LAPACKNAME(slaic1)( integer * JOB,
                        integer * J,
                        real      X[],
                        real    * SEST,
                        real      W[],
                        real    * GAMMA,
                        real    * SESTPR,
                        real    * S,
                        real    * C ) ;

    void
    LAPACKNAME(dlaic1)( integer    * JOB,
                        integer    * J,
                        doublereal   X[],
                        doublereal * SEST,
                        doublereal   W[],
                        doublereal * GAMMA,
                        doublereal * SESTPR,
                        doublereal * S,
                        doublereal * C ) ;
  }
  #endif

  /*\
   *  Purpose
   *  =======
   *
   *  DLAIC1 applies one step of incremental condition estimation in
   *  its simplest version:
   *
   *  Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
   *  lower triangular matrix L, such that
   *           twonorm(L*x) = sest
   *  Then DLAIC1 computes sestpr, s, c such that
   *  the vector
   *                  [ s*x ]
   *           xhat = [  c  ]
   *  is an approximate singular vector of
   *                  [ L     0  ]
   *           Lhat = [ w' gamma ]
   *  in the sense that
   *           twonorm(Lhat*xhat) = sestpr.
   *
   *  Depending on JOB, an estimate for the largest or smallest singular
   *  value is computed.
   *
   *  Note that [s c]' and sestpr**2 is an eigenpair of the system
   *
   *      diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
   *                                            [ gamma ]
   *
   *  where  alpha =  x'*w.
   *
   *  Arguments
   *  =========
   *
   *  JOB     (input) INTEGER
   *          = 1: an estimate for the largest singular value is computed.
   *          = 2: an estimate for the smallest singular value is computed.
   *
   *  J       (input) INTEGER
   *          Length of X and W
   *
   *  X       (input) DOUBLE PRECISION array, dimension (J)
   *          The j-vector x.
   *
   *  SEST    (input) DOUBLE PRECISION
   *          Estimated singular value of j by j matrix L
   *
   *  W       (input) DOUBLE PRECISION array, dimension (J)
   *          The j-vector w.
   *
   *  GAMMA   (input) DOUBLE PRECISION
   *          The diagonal element gamma.
   *
   *  SESTPR  (output) DOUBLE PRECISION
   *          Estimated singular value of (j+1) by (j+1) matrix Lhat.
   *
   *  S       (output) DOUBLE PRECISION
   *          Sine needed in forming xhat.
   *
   *  C       (output) DOUBLE PRECISION
   *          Cosine needed in forming xhat.
   *
  \*/

  inline
  void
  laic1( integer    JOB,
         integer    J,
         real const X[],
         real       SEST,
         real const W[],
         real       GAMMA,
         real     & SESTPR,
         real     & S,
         real     & C )
  #ifdef ALGLIN_USE_ACCELERATE
  { CLAPACKNAME(slaic1)( &JOB, &J,
                         const_cast<real*>(X), &SEST,
                         const_cast<real*>(W), &GAMMA,
                         &SESTPR, &S, &C  ) ; }
  #else
  { LAPACKNAME(slaic1)( &JOB, &J,
                        const_cast<real*>(X), &SEST,
                        const_cast<real*>(W), &GAMMA,
                        &SESTPR, &S, &C ) ; }
  #endif

  inline
  void
  laic1( integer          JOB,
         integer          J,
         doublereal const X[],
         doublereal       SEST,
         doublereal const W[],
         doublereal       GAMMA,
         doublereal     & SESTPR,
         doublereal     & S,
         doublereal     & C )
  #ifdef ALGLIN_USE_ACCELERATE
  { CLAPACKNAME(dlaic1)( &JOB, &J,
                         const_cast<doublereal*>(X), &SEST,
                         const_cast<doublereal*>(W), &GAMMA,
                         &SESTPR, &S, &C  ) ; }
  #else
  { LAPACKNAME(dlaic1)( &JOB, &J,
                        const_cast<doublereal*>(X), &SEST,
                        const_cast<doublereal*>(W), &GAMMA,
                        &SESTPR, &S, &C ) ; }
  #endif

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
         integer NB ) ;

  template <typename REAL>
  integer
  getry( integer M,
         integer N,
         REAL    A[],
         integer LDA,
         integer IPIV[],
         integer NB ) ;

  template <typename REAL>
  integer
  gtx( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] ) ;

  template <typename REAL>
  integer
  gty( integer M,
       integer N,
       REAL    A[],
       integer LDA,
       integer IPIV[] ) ;

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
               T       epsi ) ;

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  integer
  rankEstimate( integer   M,
                integer   N,
                T         A[],
                integer   LDA,
                T         RCOND,
                integer & RANK,
                T         SVAL[3] ) ;

  template <typename T>
  void
  triTikhonov( integer N,
               T const Tmat[],
               integer LDT,
               integer nrhs,
               T       RHS[],
               integer ldRHS,
               T       lambda ) ;

  inline
  bool
  outMATRIXcheck( MatrixType const & MT, integer i, integer j ) {
    bool ok = MT == FULL_MATRIX ||
              ( MT == LOWER_TRIANGULAR_MATRIX && i >= j ) ||
              ( MT == UPPER_TRIANGULAR_MATRIX && i <= j ) ;
    return ok ;
  }

  template <typename T>
  inline
  void
  outMATRIX( MatrixType const & MT,
             integer   NR,
             integer   NC,
             T const   A[],
             integer   LDA,
             ostream & s,
             integer   prec = 4,
             integer   rperm[] = nullptr,
             integer   cperm[] = nullptr ) {
    integer j0 = cperm == nullptr ? 0 : cperm[0]-1 ;
    for ( integer i = 0 ; i < NR ; ++i ) {
      integer ii = rperm == nullptr ? i : rperm[i]-1 ;
      if ( outMATRIXcheck(MT,i,0) )
        s << std::setprecision(prec) << std::setw(prec+6) << A[ii+j0*LDA] ;
      else
        s << std::setw(prec+6) << " " ;
      for ( integer j = 1 ; j < NC ; ++j ) {
        integer jj = cperm == nullptr ? j : cperm[j]-1 ;
        if ( outMATRIXcheck(MT,i,j) )
          s << " " << std::setprecision(prec) << std::setw(prec+6) << A[ii+jj*LDA] ;
        else
          s << " " << std::setw(prec+6) << " " ;
      }
      s << '\n' ;
    }
  }

  template <typename T>
  inline
  void
  outMAPLE( integer   NR,
            integer   NC,
            T const   A[],
            integer   LDA,
            ostream & s ) {
    s << "<" ;
    for ( integer j = 0 ; j < NC ; ++j ) {
      s << "<" << std::setprecision(20) << A[j*LDA] ;
      for ( integer i = 1 ; i < NR ; ++i )
        s << "," << std::setprecision(20) << A[i+j*LDA] ;
      if ( j < NC-1 ) s << ">|\n" ;
      else            s << ">>;\n" ;
    }
  }

} // end namespace alglin

#endif

///
/// eof: alglin.hh
///
