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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wpadded"
#endif

#include "Alglin++.hh"
#include <iomanip>
#include <vector>

namespace alglin {

  /*
  //   ____                             ____ ____ ___   ___  ____
  //  / ___| _ __   __ _ _ __ ___  ___ / ___/ ___/ _ \ / _ \|  _ \
  //  \___ \| '_ \ / _` | '__/ __|/ _ \ |  | |  | | | | | | | |_) |
  //   ___) | |_) | (_| | |  \__ \  __/ |__| |__| |_| | |_| |  _ <
  //  |____/| .__/ \__,_|_|  |___/\___|\____\____\___/ \___/|_| \_\
  //        |_|
  */

  //! Sparse Matrix Structure
  template <typename T>
  bool
  SparseCCOOR<T>::foundNaN() const
  { return alglin::foundNaN( vals, nnz ); }

  /*
  //   __  __       _        _     __        __
  //  |  \/  | __ _| |_ _ __(_)_  _\ \      / / __ __ _ _ __  _ __   ___ _ __
  //  | |\/| |/ _` | __| '__| \ \/ /\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
  //  | |  | | (_| | |_| |  | |>  <  \ V  V /| | | (_| | |_) | |_) |  __/ |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\  \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
  //                                                   |_|   |_|
  */

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  MatrixWrapper<T>::MatrixWrapper( valueType * _data,
                                   integer     nr,
                                   integer     nc,
                                   integer     ld )
  : numRows(nr)
  , numCols(nc)
  , ldData(ld)
  , data(_data)
  {
  #ifndef ALGLIN_NO_DEBUG
    ALGLIN_ASSERT( nr >= 0 && nc >= 0 && ldData >= nr,
                   "MatrixWrapper( data, nr=" << nr <<
                   ", nc=" << nc << ", ld=" << ld <<
                   ") bad dimensions" );
  #endif
  }

  template <typename T>
  void
  MatrixWrapper<T>::setup(
    valueType * _data,
    integer     nr,
    integer     nc,
    integer     ld
  ) {
    data    = _data;
    numRows = nr;
    numCols = nc;
    ldData  = ld;
    #ifndef ALGLIN_NO_DEBUG
    ALGLIN_ASSERT( nr >= 0 && nc >= 0 && ldData >= nr,
                   "MatrixWrapper( data, nr=" << nr <<
                   ", nc=" << nc << ", ld=" << ld <<
                   ") bad dimensions" );
    #endif
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( MatW const & A ) const {
    ALGLIN_ASSERT( A.numRows == numRows && A.numCols == numCols,
                   "MatrixWrapper::check(A) size(A) = " <<
                   A.numRows << " x " << A.numRows << " expected " <<
                   numRows << " x " << numCols );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::check( Sparse const & sp ) const {
    ALGLIN_ASSERT( sp.numRows <= numRows && sp.numCols <= numCols,
                   "MatrixWrapper::check(sp) size(sp) = " <<
                   sp.numRows << " x " << sp.numRows <<
                   " mus be contained in " << numRows << " x " << numCols );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( valueType const data_in[], integer ldData_in ) {
    integer info = gecopy( numRows, numCols, data_in, ldData_in, data, ldData );
    ALGLIN_ASSERT( info == 0,
                   "MatrixWrapper::load call alglin::gecopy return info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( MatW const & A ) {
    #ifndef ALGLIN_NO_DEBUG
    check(A);
    #endif
    integer info = gecopy( A.numRows, A.numCols, A.data, A.ldData, data, ldData );
    ALGLIN_ASSERT( info == 0,
                   "MatrixWrapper::load call alglin::gecopy return info = " << info );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load0( Sparse const & sp ) {
    #ifndef ALGLIN_NO_DEBUG
    check(sp);
    #endif
    zero();
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx],sp.cols[idx])] = sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( Sparse const & sp ) {
    #ifndef ALGLIN_NO_DEBUG
    check(sp);
    #endif
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx],sp.cols[idx])] = sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( valueType const data_in[], integer ldData_in ) {
    geadd( numRows, numCols,
           1.0, data_in, ldData_in,
           1.0, data, ldData,
           data, ldData );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( Sparse const & sp ) {
    #ifndef ALGLIN_NO_DEBUG
    check(sp);
    #endif
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx],sp.cols[idx])] += sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::add( valueType alpha, Sparse const & sp ) {
    #ifndef ALGLIN_NO_DEBUG
    check(sp);
    #endif
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx],sp.cols[idx])] += alpha * sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  void
  MatrixWrapper<T>::load( Sparse const & sp,
                          integer        i_offs,
                          integer        j_offs ) {
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] = sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  void
  MatrixWrapper<T>::add( Sparse const & sp,
                         integer        i_offs,
                         integer        j_offs ) {
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += sp.vals[idx];
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  template <typename T>
  void
  MatrixWrapper<T>::add( valueType      alpha,
                         Sparse const & sp,
                         integer        i_offs,
                         integer        j_offs ) {
    for ( integer idx = 0; idx < sp.nnz; ++idx )
      data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += alpha * sp.vals[idx];
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
                T         SVAL[3] ) {

    integer MN = std::min(M, N);
    std::vector<T> Wmin( MN ), Wmax( MN );

    // Test the input scalar arguments.
    ALGLIN_ASSERT( M >= 0 && N >= 0,
                   "rankEstimate, bad size matrix " << M << " x " << N );
    ALGLIN_ASSERT( LDA >= max_index(1,M),
                   "rankEstimate, bad leading dimension ldA = " << LDA );
    ALGLIN_ASSERT( RCOND >= 0,
                   "rankEstimate, bad condision number rcond = " << RCOND );

    // Quick return if possible
    SVAL[0] = 0;
    SVAL[1] = 0;
    SVAL[2] = 0;
    if ( MN == 0 ) return 0;

    // Determine RANK using incremental condition estimation
    integer RANK = 0;
    T SMAX = std::abs( A[0] );
    if ( SMAX > 0 ) {
      T SMIN   = SMAX;
      T SMINPR = SMIN;
      Wmin[0] = Wmax[0] = 1;
      while ( ++RANK < MN ) {
        T SMAXPR, S1, C1, S2, C2;
        T * pA0r = A + RANK * LDA;
        T & Arr  = pA0r[RANK];
        laic1( 2, RANK, &Wmin.front(), SMIN, pA0r, Arr, SMINPR, S1, C1 );
        laic1( 1, RANK, &Wmax.front(), SMAX, pA0r, Arr, SMAXPR, S2, C2 );

        if ( SMAXPR*RCOND > SMINPR ) break;

        for ( integer i=0; i < RANK; ++i )
          { Wmin[i] *= S1; Wmax[i] *= S2; }

        Wmin[RANK] = C1;
        Wmax[RANK] = C2;
        SMIN = SMINPR;
        SMAX = SMAXPR;
      }
      SVAL[0] = SMAX;
      SVAL[1] = SMIN;
      SVAL[2] = SMINPR;
    }
    return RANK;
  }

  /*\
   |   _    _   _
   |  | |  | | | |
   |  | |  | | | |
   |  | |__| |_| |
   |  |_____\___/
   |
  \*/

  template <typename T>
  LU<T>::LU()
  : Factorization<T>()
  , allocReals("allocReals")
  , allocIntegers("allocIntegers")
  {}

  template <typename T>
  LU<T>::~LU() {
    allocReals.free();
    allocIntegers.free();
  }

  template <typename T>
  void
  LU<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol+2*(nRow+nCol)) );
      allocIntegers.allocate( size_t(2*nRow) );
      this -> Amat    = allocReals( size_t(nRow*nCol) );
      this -> Work    = allocReals( size_t(2*(nRow+nCol)) );
      this -> i_pivot = allocIntegers( size_t(nRow) );
      this -> Iwork   = allocIntegers( size_t(nRow) );
    }
  }

  template <typename T>
  void
  LU<T>::factorize() {
    integer info = getrf( nRow, nCol, Amat, nRow, i_pivot );
    ALGLIN_ASSERT( info == 0, "LU::factorize getrf INFO = " << info );
  }

  template <typename T>
  void
  LU<T>::factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    ALGLIN_ASSERT( info == 0, "LU::factorize gecopy INFO = " << info );
    factorize();
  }

  template <typename T>
  void
  LU<T>::check_ls( char const who[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "LU<T>::" << who << ", rectangular matrix " <<
                   nRow << " x " << nCol );
  }

  template <typename T>
  void
  LU<T>::solve( valueType xb[] ) const {
    check_ls("solve");
    integer info = getrs( NO_TRANSPOSE,
                          nRow, 1, Amat, nRow, i_pivot,
                          xb, nRow );
    ALGLIN_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  template <typename T>
  void
  LU<T>::t_solve( valueType xb[] ) const {
    check_ls("t_solve");
    integer info = getrs( TRANSPOSE,
                          nRow, 1, Amat, nRow, i_pivot,
                          xb, nRow );
    ALGLIN_ASSERT( info == 0, "LU::t_solve getrs INFO = " << info );
  }

  template <typename T>
  void
  LU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("solve");
    integer info = getrs( NO_TRANSPOSE,
                          nRow, nrhs, Amat, nRow, i_pivot,
                          B, ldB );
    ALGLIN_ASSERT( info == 0, "LU::solve getrs INFO = " << info );
  }

  template <typename T>
  void
  LU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("t_solve");
    integer info = getrs( TRANSPOSE,
                          nRow, nrhs, Amat, nRow, i_pivot,
                          B, ldB );
    ALGLIN_ASSERT( info >= 0, "LU::t_solve getrs INFO = " << info );
  }

  template <typename T>
  typename LU<T>::valueType
  LU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gecon1( nRow, Amat, nRow,
                           norm1, rcond, Work, Iwork );
    ALGLIN_ASSERT( info == 0, "LU::cond1, gecon1 return info = " << info );
    return rcond;
  }

  template <typename T>
  typename LU<T>::valueType
  LU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = geconInf( nRow, Amat, nRow,
                             normInf, rcond, Work, Iwork );
    ALGLIN_ASSERT( info == 0, "LU::condInf, geconInf return info = " << info );
    return rcond;
  }

  /*\
   |   _    _   _ ____   ___
   |  | |  | | | |  _ \ / _ \
   |  | |  | | | | |_) | | | |
   |  | |__| |_| |  __/| |_| |
   |  |_____\___/|_|    \__\_\
   |
  \*/

  template <typename T>
  LUPQ<T>::LUPQ()
  : Factorization<T>()
  , allocReals("LUPQ-allocReals")
  , allocIntegers("LUPQ-allocIntegers")
  {}

  template <typename T>
  LUPQ<T>::~LUPQ() {
    allocReals.free();
    allocIntegers.free();
  }

  template <typename T>
  void
  LUPQ<T>::allocate( integer NR, integer NC ) {
    ALGLIN_ASSERT( NR == NC,
                   "LUPQ<T>::allocate, cannot allocate rectangular matrix " <<
                   NR << " x " << NC );
    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      allocReals.allocate( size_t(nRow*nCol) );
      Amat = allocReals( size_t(nRow*nCol) );
      allocIntegers.allocate( size_t(2*nRow) );
      ipiv = allocIntegers( size_t(nRow) );
      jpiv = allocIntegers( size_t(nRow) );
    }
  }

  template <typename T>
  void
  LUPQ<T>::factorize() {
    integer info = getc2( nRow, Amat, nRow, ipiv, jpiv );
    ALGLIN_ASSERT( info == 0, "LUPQ::factorize getrf INFO = " << info );
  }

  template <typename T>
  void
  LUPQ<T>::factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
    ALGLIN_ASSERT( NR == NC,
                   "LUPQ<T>::factorize, cannot factorize rectangular matrix " <<
                   NR << " x " << NC );
    allocate( NR, NC );
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow );
    ALGLIN_ASSERT( info == 0, "LUPQ::factorize gecopy INFO = " << info );
    factorize();
  }

  template <typename T>
  void
  LUPQ<T>::solve( valueType xb[] ) const {
    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsv( LOWER, NO_TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for U part
    trsv( UPPER, NO_TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, -1 );
  }

  template <typename T>
  void
  LUPQ<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, 1 );

    // Solve for L part
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for U part
    trsm( LEFT, UPPER, NO_TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, -1 );
  }

  template <typename T>
  void
  LUPQ<T>::t_solve( valueType xb[] ) const {
    // Apply permutations JPIV to the solution (RHS)
    swaps( 1, xb, nRow, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsv( UPPER, TRANSPOSE, NON_UNIT, nRow, Amat, nRow, xb, 1 );

    // Solve for L part
    trsv( LOWER, TRANSPOSE, UNIT, nRow, Amat, nRow, xb, 1 );

    // Apply permutations IPIV to RHS
    swaps( 1, xb, nRow, 0, nRow-2, ipiv, -1 );

  }

  template <typename T>
  void
  LUPQ<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {

    // Apply permutations JPIV to the solution (RHS)
    swaps( nrhs, B, ldB, 0, nRow-2, jpiv, 1 );

    // Solve for U part
    trsm( LEFT, UPPER, TRANSPOSE, NON_UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Solve for L part
    trsm( LEFT, LOWER, TRANSPOSE, UNIT, nRow, nrhs, 1.0, Amat, nRow, B, ldB );

    // Apply permutations IPIV to RHS
    swaps( nrhs, B, ldB, 0, nRow-2, ipiv, -1 );
  }

  /*\
   |    ___  ____
   |   / _ \|  _ \
   |  | | | | |_) |
   |  | |_| |  _ <
   |   \__\_\_| \_\
   |
  \*/
  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC, integer Lwrk ) {
    nRow       = NR;
    nCol       = NC;
    nReflector = std::min(nRow,nCol);
    Lwork      = Lwrk;
    allocReals.allocate( size_t(nRow*nCol+Lwork+nReflector) );
    Amat = allocReals( size_t(nRow*nCol) );
    Work = allocReals( size_t(Lwork) );
    Tau  = allocReals( size_t(nReflector) );
  }

  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC ) {
      valueType tmp; // get optimal allocation
      integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 );
      ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::geqrf return info = " << info );
      integer L = integer(tmp);
      if ( L < NR ) L = NR;
      if ( L < NC ) L = NC;
      allocate( NR, NC, L );
    }
  }

  template <typename T>
  void
  QR<T>::applyQ( SideMultiply  SIDE,
                 Transposition TRANS,
                 integer       nRefl,
                 integer       NR,
                 integer       NC,
                 valueType     C[],
                 integer       ldC ) const {
    ALGLIN_ASSERT( (SIDE == alglin::LEFT  && NR == nRow) ||
                   (SIDE == alglin::RIGHT && NC == nRow),
                   "QR::applyQ NR = " << NR << " NC = " << NC <<
                   " nRow = " << nRow );
    integer info = ormqr( SIDE, TRANS,
                          NR, NC,
                          nRefl,  // numero riflettori usati nel prodotto Q
                          Amat, nRow /*ldA*/,
                          Tau,
                          C, ldC,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "QR::applyQ call alglin::ormqr return info = " << info <<
                   " Lwork = " << Lwork );
  }

  template <typename T>
  void
  QR<T>::getR( valueType R[], integer ldR ) const {
    integer minRC = std::min( nRow, nCol );
    gezero( minRC, minRC, R, ldR );
    for ( integer i = 0; i < minRC; ++i )
      for ( integer j = i; j < minRC; ++j )
        R[i+j*ldR] = Amat[ i+j*nRow];
  }

  template <typename T>
  void
  QR<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve, factored matrix must be square" );
    Qt_mul(xb);
    invR_mul(xb);
  }

  template <typename T>
  void
  QR<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve_t, factored matrix must be square" );
    invRt_mul(xb);
    Q_mul(xb);
  }

  template <typename T>
  void
  QR<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve, factored matrix must be square" );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
  }

  template <typename T>
  void
  QR<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve_t, factored matrix must be square" );
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
  }

  /*\
   |    ___  ____  ____
   |   / _ \|  _ \|  _ \
   |  | | | | |_) | |_) |
   |  | |_| |  _ <|  __/
   |   \__\_\_| \_\_|
   |
  \*/
  template <typename T>
  void
  QRP<T>::permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[JPVT[i]-1] = x[i];
    copy( nCol, Work, 1, x, 1 );
  }

  template <typename T>
  void
  QRP<T>::inv_permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0; i < nCol; ++i ) Work[i] = x[JPVT[i]-1];
    copy( nCol, Work, 1, x, 1 );
  }

  template <typename T>
  void
  QRP<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QRP::solve, factored matrix must be square" );
    Qt_mul(xb);
    invR_mul(xb);
    permute(xb); // da aggiungere!
  }

  template <typename T>
  void
  QRP<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QRP::solve_t, factored matrix must be square" );
    inv_permute(xb); // da aggiungere!
    invRt_mul(xb);
    Q_mul(xb);
  }

  template <typename T>
  void
  QRP<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QRP::solve, factored matrix must be square" );
    Qt_mul( nRow, nrhs, XB, ldXB );
    invR_mul( nRow, nrhs, XB, ldXB );
    permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
  }

  template <typename T>
  void
  QRP<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QRP::solve_t, factored matrix must be square" );
    inv_permute_rows( nRow, nrhs, XB, ldXB ); // da aggiungere!
    invRt_mul( nRow, nrhs, XB, ldXB );
    Q_mul( nRow, nrhs, XB, ldXB );
  }

  /*\
   |   ______     ______
   |  / ___\ \   / /  _ \
   |  \___ \\ \ / /| | | |
   |   ___) |\ V / | |_| |
   |  |____/  \_/  |____/
   |
  \*/
  template <typename T>
  void
  SVD<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      minRC      = std::min(NR,NC);
      valueType tmp;
      integer info = gesvd( REDUCED, REDUCED,
                            NR, NC,
                            nullptr, NR,
                            nullptr,
                            nullptr, NR,
                            nullptr, minRC,
                            &tmp, -1 );
      ALGLIN_ASSERT( info == 0,
                     "alglin::SVD::allocate, in gesvd info = " << info );
      Lwork = integer(tmp);
      info = gesdd( REDUCED,
                    NR, NC,
                    nullptr, NR,
                    nullptr,
                    nullptr, NR,
                    nullptr, minRC,
                    &tmp, -1, nullptr );
       if ( integer(tmp) > Lwork ) Lwork = integer(tmp);
       allocReals.allocate( size_t(nRow*nCol+minRC*(nRow+nCol+1)+Lwork) );
       Amat = allocReals( size_t(nRow*nCol));
       Svec  = allocReals( size_t(minRC) );
       Umat  = allocReals( size_t(minRC*nRow) );
       VTmat = allocReals( size_t(minRC*nCol) );
       Work  = allocReals( size_t(Lwork) );
       allocIntegers.allocate( size_t(8*minRC) );
       IWork = allocIntegers( size_t(8*minRC) );
    }
  }

  template <typename T>
  void
  SVD<T>::factorize() {
    integer info;
    switch ( svd_used ) {
    case USE_GESVD:
      info = gesvd( REDUCED,
                    REDUCED,
                    nRow, nCol, Amat, nRow,
                    Svec,
                    Umat, nRow,
                    VTmat, minRC,
                    Work, Lwork );
      ALGLIN_ASSERT( info == 0,
                     "SVD::factorize call alglin::gesvd return info = " << info );
      break;
    case USE_GESDD:
      info = gesdd( REDUCED,
                    nRow, nCol, Amat, nRow,
                    Svec,
                    Umat, nRow,
                    VTmat, minRC,
                    Work, Lwork, IWork );
      ALGLIN_ASSERT( info == 0,
                     "SVD::factorize call alglin::gesdd return info = " << info );
      break;
    }
  }

  template <typename T>
  void
  SVD<T>::solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    Ut_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= Svec[i];
    V_mul( 1.0, Work, 1, 0.0, xb, 1 );
  }

  template <typename T>
  void
  SVD<T>::t_solve( valueType xb[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    Vt_mul( 1.0, xb, 1, 0.0, Work, 1 );
    for ( integer i = 0; i < minRC; ++i ) Work[i] /= Svec[i];
    U_mul( 1.0, Work, 1, 0.0, xb, 1 );
  }

  /*\
   |   _     ____ ____
   |  | |   / ___/ ___|
   |  | |   \___ \___ \
   |  | |___ ___) |__) |
   |  |_____|____/____/
   |
  \*/
  template <typename T>
  void
  LSS<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelss( NR, NC, 1, nullptr, NR, nullptr, NR, nullptr,
                            rcond, rank, &tmp, -1 );
      ALGLIN_ASSERT( info == 0,
                     "alglin::LSS::allocate, in gelss info = " << info );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelss( NC, NR, 1, nullptr, NC, nullptr, NC, nullptr,
                      rcond, rank, &tmp, -1 );
        ALGLIN_ASSERT( info == 0,
                       "alglin::LSS::allocate, in gelss info = " << info );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      integer minRC = std::min(NR,NC);

      allocReals.allocate( size_t(2*NR*NC+Lwork+minRC) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      sigma    = allocReals( size_t(minRC) );
      AmatWork = Amat+nRow*nCol;
    }
  }

  template <typename T>
  void
  LSS<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelss( nRow, nCol, 1,
                          AmatWork, nRow,
                          xb, nRow,
                          sigma, rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::solve (rhs=1), in gelss info = " << info );
  }

  template <typename T>
  void
  LSS<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss( nCol, nRow, 1,
                          AmatWork, nCol,
                          xb, nCol,
                          sigma, rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::t_solve (rhs=1), in gelss info = " << info );
  }

  template <typename T>
  void
  LSS<T>::solve( integer   nrhs,
                 valueType B[],
                 integer   ldB ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelss( nRow, nCol, nrhs,
                          AmatWork, nRow,
                          B, ldB,
                          sigma, rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::solve (rhs=" << nrhs <<
                   "), in gelss info = " << info );
  }

  template <typename T>
  void
  LSS<T>::t_solve( integer   nrhs,
                   valueType B[],
                   integer   ldB ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelss( nCol, nRow, nrhs,
                          AmatWork, nCol,
                          B, ldB,
                          sigma, rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::t_solve (rhs=" << nrhs <<
                   "), in gelss info = " << info );
  }

  /*\
   |  _     ______   __
   | | |   / ___\ \ / /
   | | |   \___ \\ V /
   | | |___ ___) || |
   | |_____|____/ |_|
   |
  \*/
  template <typename T>
  void
  LSY<T>::allocate( integer NR, integer NC ) {

    if ( nRow != NR || nCol != NC ) {
      nRow = NR;
      nCol = NC;
      valueType tmp;
      integer info = gelsy( NR, NC, 1, nullptr, NR, nullptr, NR, nullptr,
                            rcond, rank, &tmp, -1 );
      ALGLIN_ASSERT( info == 0,
                     "alglin::LSY::allocate, in gelss info = " << info );

      Lwork = integer(tmp);
      if ( NR != NC ) {
        info = gelsy( NC, NR, 1, nullptr, NC, nullptr, NC, nullptr,
                      rcond, rank, &tmp, -1 );
        ALGLIN_ASSERT( info == 0,
                       "alglin::LSY::allocate, in gelss info = " << info );
        if ( Lwork < integer(tmp) ) Lwork = integer(tmp);
      }

      allocReals.allocate( size_t(2*NR*NC+Lwork) );
      Amat     = allocReals( size_t(2*nRow*nCol) );
      Work     = allocReals( size_t(Lwork) );
      AmatWork = Amat+nRow*nCol;
      allocInts.allocate( size_t(NC) );
      jpvt = allocInts( size_t(NC) );
    }
  }

  template <typename T>
  void
  LSY<T>::solve( valueType xb[] ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1);
    integer info = gelsy( nRow, nCol, 1,
                          AmatWork, nRow,
                          xb, nRow, jpvt,
                          rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::solve (rhs=1), in gelss info = " << info );
  }

  template <typename T>
  void
  LSY<T>::t_solve( valueType xb[] ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy( nCol, nRow, 1,
                          AmatWork, nCol,
                          xb, nCol, jpvt,
                          rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSS::t_solve (rhs=1), in gelss info = " << info );
  }

  template <typename T>
  void
  LSY<T>::solve( integer   nrhs,
                 valueType B[],
                 integer   ldB ) const {
    // save matrix
    copy( nRow*nCol, Amat, 1, AmatWork, 1 );
    integer info = gelsy( nRow, nCol, nrhs,
                          AmatWork, nRow,
                          B, ldB, jpvt,
                          rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSD::solve (rhs=" << nrhs <<
                   "), in gelsd info = " << info );
  }

  template <typename T>
  void
  LSY<T>::t_solve( integer   nrhs,
                   valueType B[],
                   integer   ldB ) const {
    // save matrix
    for ( integer i = 0; i < nCol; ++i )
      copy( nRow, Amat+i*nRow, 1, AmatWork+i, nCol );
    integer info = gelsy( nCol, nRow, nrhs,
                          AmatWork, nCol,
                          B, ldB, jpvt,
                          rcond, rank,
                          Work, Lwork );
    ALGLIN_ASSERT( info == 0,
                   "alglin::LSD::t_solve (rhs=" << nrhs <<
                   "), in gelsd info = " << info );
  }

  //============================================================================

  template <typename valueType>
  inline
  void
  tridiag_axpy( integer         N,
                valueType       alpha,
                valueType const L[],
                valueType const D[],
                valueType const U[],
                valueType const x[],
                valueType       beta,
                valueType       y[] ) {
    if ( isZero(beta) ) {
      y[0] = alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    } else {
      y[0] = beta*y[0] + alpha*(D[0]*x[0] + U[0] * x[1]);
      for ( integer i = 1; i < N-1; ++i )
        y[i] = beta*y[i] + alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]);
      y[N-1] = beta*y[N-1] + alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]);
    }
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ ____  ____  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| / ___||  _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | \___ \| |_) | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |___) |  __/| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|____/|_|   |____/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalSPD<T>::factorize( integer         N,
                               valueType const _L[],
                               valueType const _D[] ) {
    if ( this -> nRC != N ) {
      this -> nRC = N;
      allocReals.allocate(3*N);
      L    = allocReals(N);
      D    = allocReals(N);
      WORK = allocReals(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    integer info = pttrf( N, L, D );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalSPD::factorize, return info = " << info );
  }

  template <typename T>
  T
  TridiagonalSPD<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = ptcon1( nRC, D, L, norm1, rcond, WORK );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalSPD::cond1, return info = " << info );
    return rcond;
  }

  template <typename T>
  void
  TridiagonalSPD<T>::solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalSPD::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( valueType xb[] ) const {
    integer info = pttrs( nRC, 1, D, L, xb, nRC );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalSPD::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalSPD<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalSPD::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalSPD<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = pttrs( nRC, nrhs, D, L, xb, ldXB );
    ALGLIN_ASSERT( info == 0, "TridiagonalSPD::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalSPD<T>::axpy( integer         N,
                           valueType       alpha,
                           valueType const _L[],
                           valueType const _D[],
                           valueType const x[],
                           valueType       beta,
                           valueType       y[] ) const {
    tridiag_axpy( N, alpha, _L, _D, _L, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _ _    _   _
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
   |                             |___/
  \*/
  template <typename T>
  void
  TridiagonalLU<T>::factorize( integer         N,
                               valueType const _L[],
                               valueType const _D[],
                               valueType const _U[] ) {
    if ( this -> nRC != N ) {
      this -> nRC = N;
      allocReals.allocate(6*N);
      allocIntegers.allocate(2*N);
      L     = allocReals(N);
      D     = allocReals(N);
      U     = allocReals(N);
      U2    = allocReals(N);
      WORK  = allocReals(2*N);
      IPIV  = allocIntegers(N);
      IWORK = allocIntegers(N);
    }
    copy( N, _L, 1, L, 1 );
    copy( N, _D, 1, D, 1 );
    copy( N, _U, 1, U, 1 );
    integer info = gttrf( N, L, D, U, U2, IPIV );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::factorize, return info = " << info );
  }

  template <typename T>
  T
  TridiagonalLU<T>::cond1( valueType norm1 ) const {
    valueType rcond;
    integer info = gtcon1( nRC, L, D, U, U2, IPIV, norm1, rcond, WORK, IWORK );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::cond1, return info = " << info );
    return rcond;
  }

  template <typename T>
  T
  TridiagonalLU<T>::condInf( valueType normInf ) const {
    valueType rcond;
    integer info = gtconInf( nRC, L, D, U, U2, IPIV, normInf, rcond, WORK, IWORK );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::cond1, return info = " << info );
    return rcond;
  }

  template <typename T>
  void
  TridiagonalLU<T>::solve( valueType xb[] ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( valueType xb[] ) const {
    integer info = gttrs( TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalLU<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB );
    ALGLIN_ASSERT( info == 0,
                   "TridiagonalLU::solve, return info = " << info );
  }

  template <typename T>
  void
  TridiagonalLU<T>::axpy( integer         N,
                          valueType       alpha,
                          valueType const _L[],
                          valueType const _D[],
                          valueType const _U[],
                          valueType const x[],
                          valueType       beta,
                          valueType       y[] ) const {
    tridiag_axpy( N, alpha, _L, _D, _U, x, beta, y );
  }

  //============================================================================
  /*\
   |   _____     _     _ _                               _  ___  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
   |                             |___/
  \*/

  template <typename T>
  void
  TridiagonalQR<T>::factorize( integer         N,
                               valueType const L[],
                               valueType const D[],
                               valueType const U[] ) {
    allocReals.allocate(size_t(5*(N-1)));
    this -> nRC = N;
    this -> C   = allocReals(size_t(N-1));
    this -> S   = allocReals(size_t(N-1));
    this -> BD  = allocReals(size_t(N));
    this -> BU  = allocReals(size_t(N-1));
    this -> BU2 = allocReals(size_t(N-2));

    /*\
      | d u       | d u @     | d u @     | d u @     | d u @     |
      | l d u     | 0 d u     | 0 d u @   | 0 d u @   | 0 d u @   |
      |   l d u   |   l d u   |   0 d u   |   0 d u @ |   0 d u @ |
      |     l d u |     l d u |     l d u |     0 d u |     0 d u |
      |       l d |       l d |       l d |       l d |       0 d |
    \*/

    alglin::copy( N,   D, 1, BD, 1 );
    alglin::copy( N-1, U, 1, BU, 1 );
    alglin::zero( N-2, BU2, 1 );

    normInfA = 0;
    integer i = 0;
    for (; i < N-2; ++i ) {
      valueType Li = L[i];
      rotg( BD[i], Li, C[i], S[i] );
      rot( 1, &BU[i],  1, &BD[i+1], 1, C[i], S[i] );
      rot( 1, &BU2[i], 1, &BU[i+1], 1, C[i], S[i] );
      valueType sum = std::abs(BD[i]) + std::abs(BU[i]) + std::abs(BU2[i]);
      if ( sum > normInfA ) normInfA = sum;
    }
    valueType Li = L[i];
    rotg( BD[i], Li, C[i], S[i] );
    rot( 1, &BU[i], 1, &BD[i+1], 1, C[i], S[i] );

    valueType sum = std::abs(BD[i]) + std::abs(BU[i]);
    if ( sum > normInfA ) normInfA = sum;
    sum = std::abs(BD[i+1]);
    if ( sum > normInfA ) normInfA = sum;

    // Q A = R
  }

  template <typename T>
  void
  TridiagonalQR<T>::Rsolve( valueType xb[] ) const {
    xb[nRC-1] /= BD[nRC-1];
    xb[nRC-2] = (xb[nRC-2]-BU[nRC-2]*xb[nRC-1])/BD[nRC-2];
    for ( integer i = nRC-3; i >= 0; --i )
      xb[i] = (xb[i]-BU[i]*xb[i+1]-BU2[i]*xb[i+2])/BD[i];
  }

  template <typename T>
  void
  TridiagonalQR<T>::RsolveTransposed( valueType xb[] ) const {
    xb[0] /= BD[0];
    xb[1] = (xb[1]-BU[0]*xb[0])/BD[1];
    for ( integer i = 2; i < nRC; ++i )
      xb[i] = (xb[i]-BU[i]*xb[i-1]-BU2[i]*xb[i-2])/BD[i];
  }

  template <typename T>
  void
  TridiagonalQR<T>::solve( valueType xb[] ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], S[i] );
    Rsolve( xb );
  }

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( valueType xb[] ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    RsolveTransposed( xb );
    // applico Q^T b
    for ( integer i = nRC-2; i >= 0; --i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], -S[i] );
  }

  template <typename T>
  void
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], S[i] );
    for ( integer i = 0; i < nrhs; ++i )
      Rsolve( xb+i*ldXB );
  }

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    for ( integer i = 0; i < nrhs; ++i )
      RsolveTransposed(xb+i*ldXB);
    for ( integer i = nRC-2; i >= 0; --i )
      rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], -S[i] );
  }

  template <typename T>
  void
  TridiagonalQR<T>::axpy( integer         N,
                          valueType       alpha,
                          valueType const L[],
                          valueType const D[],
                          valueType const U[],
                          valueType const x[],
                          valueType       beta,
                          valueType       y[] ) const {
    tridiag_axpy( N, alpha, L, D, U, x, beta, y );
  }

  /*\
   *
   *  Solve
   *
   *  min || T x - b ||^2 + lambda ||x||^2
   *
   *
   *  / * * *       \
   *  |   * * *     |
   *  |     * * *   |
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
  TridiagonalQR<T>::lsq( integer nrhs,
                         T       RHS[],
                         integer ldRHS,
                         T       lambda_in) const {

    valueType lambda = normInfA * lambda_in;
    std::vector<T> D(nRC),
                   U(nRC-1),
                   U2(nRC-2),
                   tmp(nrhs);
    T CC, SS;

    for ( integer i = 0; i < nRC-1; ++i )
      rot( nrhs, RHS+i, ldRHS, RHS+i+1, ldRHS, C[i], S[i] );

    copy( nRC,   BD,  1, &D.front(),  1 );
    copy( nRC-1, BU,  1, &U.front(),  1 );
    copy( nRC-2, BU2, 1, &U2.front(), 1 );
    T line[3];
    integer i = 0;
    while ( i < nRC-1 ) {
      line[0] = line[2] = 0; line[1] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      integer j = i;
      while ( j < nRC-2 ) {
        line[0] = line[1];
        line[1] = line[2];
        line[2] = 0;
        rotg( D[j], line[0], CC, SS );
        rot( 1, &U[j],  1, &line[1], 1, CC, SS );
        rot( 1, &U2[j], 1, &line[2], 1, CC, SS );
        rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
        ++j;
      }
      // penultima
      rotg( D[j], line[1], CC, SS );
      rot( 1, &U[j],  1, &line[2], 1, CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
      ++j;

      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      // ultima
      line[2] = lambda;
      rotg( D[j], line[2], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );

      ++i; // next lambda
    }
    if ( nRC > 0 ) {
      integer j = nRC-1;
      line[0] = lambda;
      std::fill( tmp.begin(), tmp.end(), T(0) );
      rotg( D[j], line[0], CC, SS );
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS );
    }

    for ( integer j = 0; j < nrhs; ++j ) {
      T * xb = RHS + j*ldRHS;
      xb[nRC-1] /= D[nRC-1];
      xb[nRC-2] = (xb[nRC-2]-U[nRC-2]*xb[nRC-1])/D[nRC-2];
      for ( integer k = nRC-3; k >= 0; --k )
        xb[k] = (xb[k]-U[k]*xb[k+1]-U2[k]*xb[k+2])/D[k];
    }
  }

  /*\
   |   ____                  _          _ __  __       _        _
   |  | __ )  __ _ _ __   __| | ___  __| |  \/  | __ _| |_ _ __(_)_  __
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |\/| |/ _` | __| '__| \ \/ /
   |  | |_) | (_| | | | | (_| |  __/ (_| | |  | | (_| | |_| |  | |>  <
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|_|  |_|\__,_|\__|_|  |_/_/\_\
  \*/

  template <typename T>
  BandedLU<T>::BandedLU()
  : allocReals("_BandedLU_reals")
  , allocIntegers("_BandedLU_integers")
  , m(0)
  , n(0)
  , nL(0)
  , nU(0)
  , ldAB(0)
  , is_factorized(false)
  {}

  template <typename T>
  BandedLU<T>::~BandedLU()
  {}

  //! base class for linear system solver
  template <typename T>
  void
  BandedLU<T>::setup( integer _m,
                      integer _n,
                      integer _nL,
                      integer _nU ) {
    m    = _m;
    n    = _n;
    nL   = _nL;
    nU   = _nU;
    ldAB = 2*nL+nU+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    allocIntegers.allocate(m);
    AB   = allocReals( nnz );
    ipiv = allocIntegers( m );
    is_factorized = false;
  }

  template <typename T>
  void
  BandedLU<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    ALGLIN_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    ALGLIN_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  template <typename T>
  void
  BandedLU<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    ALGLIN_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, 1, AB, ldAB, ipiv, xb, m );
    ALGLIN_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  template <typename T>
  void
  BandedLU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    ALGLIN_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    ALGLIN_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( NO_TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    ALGLIN_ASSERT( info == 0, "BandedLU::solve, info = " << info );
  }

  template <typename T>
  void
  BandedLU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    ALGLIN_ASSERT( is_factorized, "BandedLU::solve, matrix not yet factorized" );
    ALGLIN_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrs( TRANSPOSE, m, nL, nU, nrhs, AB, ldAB, ipiv, B, ldB );
    ALGLIN_ASSERT( info == 0, "BandedLU::t_solve, info = " << info );
  }

  template <typename T>
  void
  BandedLU<T>::factorize() {
    ALGLIN_ASSERT( !is_factorized, "BandedLU::solve, matrix yet factorized" );
    ALGLIN_ASSERT( m == n, "BandedLU::solve, matrix must be square" );
    integer info = gbtrf( m, n, nL, nU, AB, ldAB, ipiv );
    ALGLIN_ASSERT( info == 0, "BandedLU::factorize, info = " << info );
    is_factorized = true;
  }

  template <typename T>
  void
  BandedLU<T>::zero() {
    integer nnz = m*(2*nL+nU+1);
    alglin::zero( nnz, AB, 1 );
    is_factorized = false;
  }

  template <typename T>
  void
  BandedLU<T>::iaddr_check( integer i, integer j ) const {
    ALGLIN_ASSERT( i >= 0 && i < m && j >= 0 && j < n,
                   "BandedLU:iaddr_check( " << i << " , " << j <<
                   " ) out of range" );
    ALGLIN_ASSERT( j >= i-nL && j <= i+nU,
                   "BandedLU:iaddr_check( " << i << " , " << j <<
                   " ) out of band" );
  }

  template <typename T>
  void
  BandedLU<T>::load_block( integer         nr,
                           integer         nc,
                           valueType const B[],
                           integer         ldB,
                           integer         irow,
                           integer         icol ) {

    ALGLIN_ASSERT( !is_factorized,
                   "BandedLU::load_block, matrix is factorized" );

    #if 1
    for ( integer r = 0; r < nr; ++r )
      for ( integer c = 0; c < nc; ++c )
        AB[iaddr( irow+r, icol+c )] = B[ r + c * ldB ];
    #else
    // must be checked
    iaddr_check( irow,      icol      );
    iaddr_check( irow+nr-1, icol+nc-1 );
    iaddr_check( irow,      icol+nc-1 );
    iaddr_check( irow+nr-1, icol      );

    // copy by diagonal
    for ( integer r = 0; r < nr; ++r ) {
      integer ia = iaddr( irow+r, icol );
      copy( std::min(nr-r,nc), B+r, ldB+1, AB+ia, ldAB );
    }
    for ( integer c = 1; c < nc; ++c ) {
      integer ia = iaddr( irow, icol+c );
      copy( std::min(nc-c,nr), B+c*ldB, ldB+1, AB+ia, ldAB );
    }
    #endif
  }

  // y <- beta*y + alpha*A*x
  /*
    +---------+
    | \       |
    |  \      |
    +---+-----+
  */
  template <typename T>
  void
  BandedLU<T>::aAxpy( valueType       alpha,
                      valueType const x[],
                      valueType       y[] ) const {

    valueType const * col = AB + nL;
    for ( integer j = 0; j < n; ++j, col += ldAB ) {
      integer imin  = j-nU;
      integer imax  = std::min(j+nL,m-1);
      integer imin0 = imin > 0 ? imin : 0;
      alglin::axpy( imax-imin0+1,
                    alpha*x[j],
                    col+imin0-imin, 1,
                    y+imin0,        1 );
    }
  }

  template <typename T>
  void
  BandedLU<T>::dump( std::ostream & stream ) const {
    for ( integer i = 0; i <= nL+nU; ++i ) {
      valueType const * col = AB + nL + i;
      for ( integer j = 0; j < n; ++j, col += ldAB )
        stream << std::setw(10) << col[0] << ' ';
      stream << '\n';
    }
  }

  // --------------------------------------------------------------

  template <typename T>
  BandedSPD<T>::BandedSPD()
  : allocReals("_BandedSPD_reals")
  , n(0)
  , nD(0)
  , ldAB(0)
  , is_factorized(false)
  {}

  template <typename T>
  BandedSPD<T>::~BandedSPD()
  {}

  //! base class for linear system solver
  template <typename T>
  void
  BandedSPD<T>::setup( ULselect _UPLO,
                       integer  _N,
                       integer  _nD ) {
    UPLO = _UPLO;
    n    = _N;
    nD   = _nD;
    ldAB = nD+1;
    integer nnz = n*ldAB;
    allocReals.allocate( nnz );
    AB   = allocReals( nnz );
    is_factorized = false;
  }

  template <typename T>
  void
  BandedSPD<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( is_factorized,
                   "BandedSPD::solve, matrix not yet factorized" );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    ALGLIN_ASSERT( info == 0,
                   "BandedSPD::solve, info = " << info );
  }

  template <typename T>
  void
  BandedSPD<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( is_factorized,
                   "BandedSPD::solve, matrix not yet factorized" );
    integer info = pbtrs( UPLO, n, nD, 1, AB, ldAB, xb, n );
    ALGLIN_ASSERT( info == 0,
                   "BandedSPD::t_solve, info = " << info );
  }

  template <typename T>
  void
  BandedSPD<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    ALGLIN_ASSERT( is_factorized,
                   "BandedSPD::solve, matrix not yet factorized" );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    ALGLIN_ASSERT( info == 0,
                   "BandedSPD::solve, info = " << info );
  }

  template <typename T>
  void
  BandedSPD<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    ALGLIN_ASSERT( is_factorized,
                   "BandedSPD::solve, matrix not yet factorized" );
    integer info = pbtrs( UPLO, n, nD, nrhs, AB, ldAB, B, ldB );
    ALGLIN_ASSERT( info == 0,
                   "BandedSPD::t_solve, info = " << info );
  }

  template <typename T>
  void
  BandedSPD<T>::factorize() {
    ALGLIN_ASSERT( !is_factorized,
                   "BandedSPD::solve, matrix yet factorized" );
    integer info = pbtrf( UPLO, n, nD, AB, ldAB );
    ALGLIN_ASSERT( info == 0,
                   "BandedSPD::factorize, info = " << info );
    is_factorized = true;
  }

  template <typename T>
  void
  BandedSPD<T>::zero() {
    alglin::zero( n*ldAB, AB, 1 );
    is_factorized = false;
  }

  template integer rankEstimate( integer   M,
                                 integer   N,
                                 real      A[],
                                 integer   LDA,
                                 real      RCOND,
                                 real      SVAL[3] );

  template integer rankEstimate( integer    M,
                                 integer    N,
                                 doublereal A[],
                                 integer    LDA,
                                 doublereal RCOND,
                                 doublereal SVAL[3] );


  template class SparseCCOOR<real>;
  template class SparseCCOOR<doublereal>;

  template class MatrixWrapper<real>;
  template class MatrixWrapper<doublereal>;

  template class LU<real>;
  template class LU<doublereal>;

  template class LUPQ<real>;
  template class LUPQ<doublereal>;

  template class QR<real>;
  template class QR<doublereal>;

  template class QRP<real>;
  template class QRP<doublereal>;

  template class SVD<real>;
  template class SVD<doublereal>;

  template class LSS<real>;
  template class LSS<doublereal>;

  template class LSY<real>;
  template class LSY<doublereal>;

  template class TridiagonalSPD<real>;
  template class TridiagonalSPD<doublereal>;

  template class TridiagonalLU<real>;
  template class TridiagonalLU<doublereal>;

  template class TridiagonalQR<real>;
  template class TridiagonalQR<doublereal>;

  template class BandedLU<real>;
  template class BandedLU<doublereal>;

  template class BandedSPD<real>;
  template class BandedSPD<doublereal>;

} // end namespace alglin

///
/// eof: Alglin++.cc
///
