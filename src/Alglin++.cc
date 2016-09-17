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

#include "Alglin++.hh"
#include <iomanip>
#include <vector>

namespace alglin {
  /*
  //   _    _   _
  //  | |  | | | |
  //  | |  | | | |
  //  | |__| |_| |
  //  |_____\___/
  */
  template <typename T>
  LU<T>::LU()
  : nRow(0)
  , nCol(0)
  , Lwork(0)
  , allocReals("allocReals")
  , allocIntegers("allocIntegers")
  {}

  template <typename T>
  LU<T>::~LU() {
    allocReals.free() ;
    allocIntegers.free() ;
  }

  template <typename T>
  void
  LU<T>::factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
    if ( nRow != NR || nCol != NC ) {
      nRow = NR ;
      nCol = NC ;
      allocReals.allocate(size_t(nRow*nCol+2*(nRow+nCol))) ;
      allocIntegers.allocate(size_t(2*nRow)) ;
      this -> Amat    = allocReals(size_t(nRow*nCol)) ;
      this -> Work    = allocReals(size_t(2*(nRow+nCol))) ;
      this -> i_pivot = allocIntegers(size_t(nRow)) ;
      this -> Iwork   = allocIntegers(size_t(nRow)) ;
    }
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "LU::factorize gecopy INFO = " << info ) ;
    info = getrf( nRow, nCol, Amat, nRow, i_pivot ) ;
    ALGLIN_ASSERT( info == 0, "LU::factorize getrf INFO = " << info ) ;
  }

  template <typename T>
  void
  LU<T>::check_ls( char const who[] ) const {
    ALGLIN_ASSERT( nRow == nCol, "LU<T>::" << who << ", rectangular matrix " << nRow << " x " << nCol ) ;
  }

  template <typename T>
  void
  LU<T>::solve( valueType xb[] ) const {
    check_ls("solve") ;
    integer info = getrs( NO_TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow ) ;
    ALGLIN_ASSERT( info == 0, "LU::solve getrs INFO = " << info ) ;
  }

  template <typename T>
  void
  LU<T>::t_solve( valueType xb[] ) const {
    check_ls("t_solve") ;
    integer info = getrs( TRANSPOSE, nRow, 1, Amat, nRow, i_pivot, xb, nRow ) ;
    ALGLIN_ASSERT( info == 0, "LU::t_solve getrs INFO = " << info ) ;
  }

  template <typename T>
  void
  LU<T>::solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("solve") ;
    integer info = getrs( NO_TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB ) ;
    ALGLIN_ASSERT( info == 0, "LU::solve getrs INFO = " << info ) ;
  }

  template <typename T>
  void
  LU<T>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    check_ls("t_solve") ;
    integer info = getrs( TRANSPOSE, nRow, nrhs, Amat, nRow, i_pivot, B, ldB ) ;
    ALGLIN_ASSERT( info >= 0, "LU::t_solve getrs INFO = " << info ) ;
  }

  template <typename T>
  T
  LU<T>::cond1( valueType norm1 ) const {
    valueType rcond ;
    integer info = gecon1( nRow, Amat, nRow, norm1, rcond, Work, Iwork ) ;
    ALGLIN_ASSERT( info == 0, "LU::cond1, gecon1 return info = " << info ) ;
    return rcond ;
  }

  template <typename T>
  T
  LU<T>::condInf( valueType normInf ) const {
    valueType rcond ;
    integer info = geconInf( nRow, Amat, nRow, normInf, rcond, Work, Iwork ) ;
    ALGLIN_ASSERT( info == 0, "LU::condInf, geconInf return info = " << info ) ;
    return rcond ;
  }

  /*
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
  */
  template <typename T>
  void
  QR<T>::allocate( integer NR, integer NC ) {
    if ( nRow != NR || nCol != NC ) {
      nRow       = NR ;
      nCol       = NC ;
      nReflector = min(nRow,nCol) ;
      valueType tmp ; // get optimal allocation
      integer info = geqrf( NR, NC, nullptr, NR, nullptr, &tmp, -1 ) ;
      ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::geqrf return info = " << info ) ;
      Lwork = integer(tmp) ;
      allocReals.allocate(size_t(nRow*nCol+Lwork+nReflector)) ;
      Amat = allocReals(size_t(nRow*nCol)) ;
      Work = allocReals(size_t(Lwork)) ;
      Tau  = allocReals(size_t(nReflector)) ;
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
                   "QR::applyQ NR = " << NR << " NC = " << NC << " nRow = " << nRow ) ;
    integer info = ormqr( SIDE, TRANS,
                          NR, NC,
                          nRefl,  // numero riflettori usati nel prodotto Q
                          Amat, nRow /*ldA*/,
                          Tau,
                          C, ldC,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "QR::applyQ call alglin::ormqr return info = " << info ) ;
  }
      
  template <typename T>
  void
  QR<T>::getR( valueType R[], integer ldR ) const {
    integer minRC = std::min( nRow, nCol ) ;
    gezero( minRC, minRC, R, ldR ) ;
    for ( integer i = 0 ; i < minRC ; ++i )
      for ( integer j = i ; j < minRC ; ++j )
        R[i+j*ldR] = Amat[ i+j*nRow] ;
  }

  template <typename T>
  void
  QR<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve, factored matrix must be square" ) ;
    Qt_mul(xb) ;
    invR_mul(xb) ;
  }

  template <typename T>
  void
  QR<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve_t, factored matrix must be square" ) ;
    invRt_mul(xb) ;
    Q_mul(xb) ;
  }

  template <typename T>
  void
  QR<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve, factored matrix must be square" ) ;
    Qt_mul( nRow, nrhs, XB, ldXB ) ;
    invR_mul( nRow, nrhs, XB, ldXB ) ;
  }

  template <typename T>
  void
  QR<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( nRow == nCol,
                   "in QR::solve_t, factored matrix must be square" ) ;
    invRt_mul( nRow, nrhs, XB, ldXB ) ;
    Q_mul( nRow, nrhs, XB, ldXB ) ;
  }

  /*
  //    ___  ____  ____
  //   / _ \|  _ \|  _ \
  //  | | | | |_) | |_) |
  //  | |_| |  _ <|  __/
  //   \__\_\_| \_\_|
  */

  template <typename T>
  void
  QRP<T>::permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0 ; i < this->nCol ; ++i ) this->Work[JPVT[i]-1] = x[i] ;
    copy( this->nCol, this->Work, 1, x, 1 ) ;
  }

  template <typename T>
  void
  QRP<T>::inv_permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0 ; i < this->nCol ; ++i ) this->Work[i] = x[JPVT[i]-1] ;
    copy( this->nCol, this->Work, 1, x, 1 ) ;
  }

  template <typename T>
  void
  QRP<T>::solve( valueType xb[] ) const {
    ALGLIN_ASSERT( this->nRow == this->nCol,
                   "in QRP::solve, factored matrix must be square" ) ;
    this->Qt_mul(xb) ;
    this->invR_mul(xb) ;
    this->permute(xb) ; // da aggiungere!
  }

  template <typename T>
  void
  QRP<T>::t_solve( valueType xb[] ) const {
    ALGLIN_ASSERT( this->nRow == this->nCol,
                   "in QRP::solve_t, factored matrix must be square" ) ;
    this->inv_permute(xb) ; // da aggiungere!
    this->invRt_mul(xb) ;
    this->Q_mul(xb) ;
  }

  template <typename T>
  void
  QRP<T>::solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( this->nRow == this->nCol,
                   "in QRP::solve, factored matrix must be square" ) ;
    this->Qt_mul( this->nRow, nrhs, XB, ldXB ) ;
    this->invR_mul( this->nRow, nrhs, XB, ldXB ) ;
    this->permute_rows( this->nRow, nrhs, XB, ldXB ) ; // da aggiungere!
  }

  template <typename T>
  void
  QRP<T>::t_solve( integer nrhs, valueType XB[], integer ldXB ) const {
    ALGLIN_ASSERT( this->nRow == this->nCol,
                   "in QRP::solve_t, factored matrix must be square" ) ;
    this->inv_permute_rows( this->nRow, nrhs, XB, ldXB ) ; // da aggiungere!
    this->invRt_mul( this->nRow, nrhs, XB, ldXB ) ;
    this->Q_mul( this->nRow, nrhs, XB, ldXB ) ;
  }

  /*
  //   ______     ______
  //  / ___\ \   / /  _ \
  //  \___ \\ \ / /| | | |
  //   ___) |\ V / | |_| |
  //  |____/  \_/  |____/
  */
  template <typename T>
  void
  SVD<T>::allocate( integer NR, integer NC, integer LDA ) {
  
    if ( nRow != NR || nCol != NC ) {
      nRow  = NR ;
      nCol  = NC ;
      minRC = min(NR,NC) ;
      valueType tmp ;
      integer info = gesvd( REDUCED, REDUCED,
                            NR, NC,
                            nullptr, LDA,
                            nullptr,
                            nullptr, NR,
                            nullptr, minRC,
                            &tmp, -1 ) ;
      ALGLIN_ASSERT( info == 0,
                     "alglin::SVD::allocate, in gesvd info = " << info ) ;
      Lwork = integer(tmp) ;
      info = gesdd( REDUCED,
                    NR, NC,
                    nullptr, LDA,
                    nullptr,
                    nullptr, NR,
                    nullptr, minRC,
                    &tmp, -1, nullptr ) ;
       if ( integer(tmp) > Lwork ) Lwork = integer(tmp) ;
       allocReals.allocate(nRow*nCol+minRC*(nRow+nCol+1)+Lwork) ;
       Amat    = allocReals(nRow*nCol) ;
       Svec    = allocReals(minRC) ;
       Umat    = allocReals(minRC*nRow) ;
       VTmat   = allocReals(minRC*nCol) ;
       Work    = allocReals(Lwork) ;
       allocIntegers.allocate(8*minRC) ;
       IWork   = allocIntegers(8*minRC) ;
    }
  }

  template <typename T>
  void
  SVD<T>::factorize( integer         NR,
                     integer         NC,
                     valueType const A[],
                     integer         LDA ) {

    allocate( NR, NC, LDA ) ;
    integer info = gecopy( nRow, nCol, A, LDA, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gecopy return info = " << info ) ;
    if ( use_gesvd ) {
      info = gesvd( REDUCED,
                    REDUCED,
                    nRow, nCol, Amat, nRow,
                    Svec,
                    Umat, nRow,
                    VTmat, minRC,
                    Work, Lwork ) ;
      ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gesvd return info = " << info ) ;
    } else {
      info = gesdd( REDUCED,
                    nRow, nCol, Amat, nRow,
                    Svec,
                    Umat, nRow,
                    VTmat, minRC,
                    Work, Lwork, IWork ) ;
      ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gesdd return info = " << info ) ;
    }
  }

  template <typename T>
  void
  SVD<T>::solve( valueType const in[], valueType out[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    Ut_mul( 1.0, in, 1, 0.0, Work, 1 ) ;
    for ( integer i = 0 ; i < minRC ; ++i ) Work[i] /= Svec[i] ;
    V_mul( 1.0, Work, 1, 0.0, out, 1 ) ;
  }

  template <typename T>
  void
  SVD<T>::solve( integer         nrhs,
                 valueType const rhs[], integer ldRHS,
                 valueType       x[],   integer ldX ) const {
    for ( integer i = 0 ; i < nrhs ; ++i ) solve( rhs + i*ldRHS, x + i*ldX ) ; 
  }

  template <typename T>
  void
  SVD<T>::t_solve( valueType const in[], valueType out[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    Vt_mul( 1.0, in, 1, 0.0, Work, 1 ) ;
    for ( integer i = 0 ; i < minRC ; ++i ) Work[i] /= Svec[i] ;
    U_mul( 1.0, Work, 1, 0.0, out, 1 ) ;
  }

  template <typename T>
  void
  SVD<T>::t_solve( integer         nrhs,
                   valueType const rhs[], integer ldRHS,
                   valueType       x[],   integer ldX ) const {
    for ( integer i = 0 ; i < nrhs ; ++i ) t_solve( rhs + i*ldRHS, x + i*ldX ) ;
  }

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
      y[0] = alpha*(D[0]*x[0] + U[0] * x[1]) ;
      for ( integer i = 1 ; i < N-1 ; ++i )
        y[i] = alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]) ;
      y[N-1] = alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]) ;
    } else {
      y[0] = beta*y[0] + alpha*(D[0]*x[0] + U[0] * x[1]) ;
      for ( integer i = 1 ; i < N-1 ; ++i )
        y[i] = beta*y[i] + alpha*(D[i]*x[i] + U[i] * x[i+1] + L[i-1] * x[i-1]) ;
      y[N-1] = beta*y[N-1] + alpha*(D[N-1]*x[N-1] + L[N-2] * x[N-2]) ;
    }
  }

  //============================================================================
  /*
  //   _____     _     _ _                               _ _    _   _
  //  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| | |  | | | |
  //    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | |  | | | |
  //    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |__| |_| |
  //    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|_____\___/
  //                             |___/
  */
  template <typename T>
  void
  TridiagonalLU<T>::factorize( integer         N,
                               valueType const _L[],
                               valueType const _D[],
                               valueType const _U[] ) {
    if ( this -> nRC != N ) {
      this -> nRC = N ;
      allocReals.allocate(6*N) ;
      allocIntegers.allocate(2*N) ;
      L     = allocReals(N) ;
      D     = allocReals(N) ;
      U     = allocReals(N) ;
      U2    = allocReals(N) ;
      WORK  = allocReals(2*N) ;
      IPIV  = allocIntegers(N) ;
      IWORK = allocIntegers(N) ;
    }
    copy( N, _L, 1, L, 1 ) ;
    copy( N, _D, 1, D, 1 ) ;
    copy( N, _U, 1, U, 1 ) ;
    integer info = gttrf( N, L, D, U, U2, IPIV ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::factorize, return info = " << info ) ;
  }
  
  template <typename T>
  T
  TridiagonalLU<T>::cond1( valueType norm1 ) const {
    valueType rcond ;
    integer info = gtcon1( nRC, L, D, U, U2, IPIV, norm1, rcond, WORK, IWORK ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::cond1, return info = " << info ) ;
    return rcond ;
  }

  template <typename T>
  T
  TridiagonalLU<T>::condInf( valueType normInf ) const {
    valueType rcond ;
    integer info = gtconInf( nRC, L, D, U, U2, IPIV, normInf, rcond, WORK, IWORK ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::cond1, return info = " << info ) ;
    return rcond ;
  }

  template <typename T>
  void
  TridiagonalLU<T>::solve( valueType xb[] ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
  }

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( valueType xb[] ) const {
    integer info = gttrs( TRANSPOSE, nRC, 1, L, D, U, U2, IPIV, xb, nRC ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
  }

  template <typename T>
  void
  TridiagonalLU<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( NO_TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
  }

  template <typename T>
  void
  TridiagonalLU<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    integer info = gttrs( TRANSPOSE, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB ) ;
    ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
  }

  template <typename T>
  void
  TridiagonalLU<T>::axpy( integer         N,
                          valueType       alpha,
                          valueType const L[],
                          valueType const D[],
                          valueType const U[],
                          valueType const x[],
                          valueType       beta,
                          valueType       y[] ) const {
    tridiag_axpy( N, alpha, L, D, U, x, beta, y ) ;
  }

  //============================================================================
  /*
  //   _____     _     _ _                               _  ___  ____
  //  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
  //    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
  //    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
  //    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
  //                             |___/
  */

  template <typename T>
  void
  TridiagonalQR<T>::factorize( integer         N,
                               valueType const L[],
                               valueType const D[],
                               valueType const U[] ) {
    allocReals.allocate(size_t(5*(N-1))) ;
    this -> nRC = N ;
    this -> C   = allocReals(size_t(N-1)) ;
    this -> S   = allocReals(size_t(N-1)) ;
    this -> BD  = allocReals(size_t(N)) ;
    this -> BU  = allocReals(size_t(N-1)) ;
    this -> BU2 = allocReals(size_t(N-2)) ;

    /*\
      | d u       | d u @     | d u @     | d u @     | d u @     |
      | l d u     | 0 d u     | 0 d u @   | 0 d u @   | 0 d u @   |
      |   l d u   |   l d u   |   0 d u   |   0 d u @ |   0 d u @ |
      |     l d u |     l d u |     l d u |     0 d u |     0 d u |
      |       l d |       l d |       l d |       l d |       0 d |
    \*/

    alglin::copy( N,   D, 1, BD, 1 ) ;
    alglin::copy( N-1, U, 1, BU, 1 ) ;
    alglin::zero( N-2, BU2, 1 ) ;

    normInfA = 0 ;
    integer i = 0 ;
    for ( ; i < N-2 ; ++i ) {
      valueType Li = L[i] ;
      rotg( BD[i], Li, C[i], S[i] ) ;
      rot( 1, &BU[i],  1, &BD[i+1], 1, C[i], S[i] ) ;
      rot( 1, &BU2[i], 1, &BU[i+1], 1, C[i], S[i] ) ;
      valueType sum = std::abs(BD[i]) + std::abs(BU[i]) + std::abs(BU2[i]) ;
      if ( sum > normInfA ) normInfA = sum ;
    }
    valueType Li = L[i] ;
    rotg( BD[i], Li, C[i], S[i] ) ;
    rot( 1, &BU[i], 1, &BD[i+1], 1, C[i], S[i] ) ;

    valueType sum = std::abs(BD[i]) + std::abs(BU[i]) ;
    if ( sum > normInfA ) normInfA = sum ;
    sum = std::abs(BD[i+1]) ;
    if ( sum > normInfA ) normInfA = sum ;

    // Q A = R
  }

  template <typename T>
  void
  TridiagonalQR<T>::Rsolve( valueType xb[] ) const {
    xb[nRC-1] /= BD[nRC-1] ;
    xb[nRC-2] = (xb[nRC-2]-BU[nRC-2]*xb[nRC-1])/BD[nRC-2] ;
    for ( integer i = nRC-3 ; i >= 0 ; --i )
      xb[i] = (xb[i]-BU[i]*xb[i+1]-BU2[i]*xb[i+2])/BD[i] ;
  }

  template <typename T>
  void
  TridiagonalQR<T>::RsolveTransposed( valueType xb[] ) const {
    xb[0] /= BD[0] ;
    xb[1] = (xb[1]-BU[0]*xb[0])/BD[1] ;
    for ( integer i = 2 ; i < nRC ; ++i )
      xb[i] = (xb[i]-BU[i]*xb[i-1]-BU2[i]*xb[i-2])/BD[i] ;
  }

  template <typename T>
  void
  TridiagonalQR<T>::solve( valueType xb[] ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0 ; i < nRC-1 ; ++i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], S[i] ) ;
    Rsolve( xb ) ;
  }

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( valueType xb[] ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    RsolveTransposed( xb ) ;
    // applico Q^T b
    for ( integer i = nRC-2 ; i >= 0 ; --i )
      rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], -S[i] ) ;
  }

  template <typename T>
  void
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A x = b --> Q A x = Q b --> R x = Q b
    // applico Q b
    for ( integer i = 0 ; i < nRC-1 ; ++i ) rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], S[i] ) ;
    for ( integer i = 0 ; i < nrhs  ; ++i ) Rsolve( xb+i*ldXB ) ;
  }

  template <typename T>
  void
  TridiagonalQR<T>::t_solve( integer nrhs, valueType xb[], integer ldXB ) const {
    // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
    for ( integer i = 0     ; i < nrhs ; ++i ) RsolveTransposed(xb+i*ldXB) ;
    for ( integer i = nRC-2 ; i >= 0   ; --i ) rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], -S[i] ) ;
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
    tridiag_axpy( N, alpha, L, D, U, x, beta, y ) ;
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

    valueType lambda = normInfA * lambda_in ;
    std::vector<T> D(nRC),
                   U(nRC-1),
                   U2(nRC-2),
                   tmp(nrhs) ;
    T CC, SS ;

    for ( integer i = 0 ; i < nRC-1 ; ++i )
      rot( nrhs, RHS+i, ldRHS, RHS+i+1, ldRHS, C[i], S[i] ) ;

    copy( nRC,   BD,  1, &D.front(),  1 ) ;
    copy( nRC-1, BU,  1, &U.front(),  1 ) ;
    copy( nRC-2, BU2, 1, &U2.front(), 1 ) ;
    T line[3] ;
    integer i = 0 ;
    while ( i < nRC-1 ) {
      line[0] = line[2] = 0 ; line[1] = lambda ;
      std::fill( tmp.begin(), tmp.end(), T(0) ) ;
      integer j = i ;
      while ( j < nRC-2 ) {
        line[0] = line[1] ;
        line[1] = line[2] ;
        line[2] = 0 ;
        rotg( D[j], line[0], CC, SS ) ;
        rot( 1, &U[j],  1, &line[1], 1, CC, SS ) ;
        rot( 1, &U2[j], 1, &line[2], 1, CC, SS ) ;
        rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS ) ;
        ++j ;
      }
      // penultima
      rotg( D[j], line[1], CC, SS ) ;
      rot( 1, &U[j],  1, &line[2], 1, CC, SS ) ;
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS ) ;
      ++j ;

      rotg( D[j], line[2], CC, SS ) ;
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS ) ;

      // ultima
      line[2] = lambda ;
      rotg( D[j], line[2], CC, SS ) ;
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS ) ;

      ++i ; // next lambda
    }
    if ( nRC > 0 ) {
      integer j = nRC-1 ;
      line[0] = lambda ;
      std::fill( tmp.begin(), tmp.end(), T(0) ) ;
      rotg( D[j], line[0], CC, SS ) ;
      rot( nrhs, RHS+j, ldRHS, &tmp.front(), 1, CC, SS ) ;
    }

    for ( integer j = 0 ; j < nrhs ; ++j ) {
      T * xb = RHS + j*ldRHS ;
      xb[nRC-1] /= D[nRC-1] ;
      xb[nRC-2] = (xb[nRC-2]-U[nRC-2]*xb[nRC-1])/D[nRC-2] ;
      for ( integer k = nRC-3 ; k >= 0 ; --k )
        xb[k] = (xb[k]-U[k]*xb[k+1]-U2[k]*xb[k+2])/D[k] ;
    }
  }

  /*
  //   _                   _   ____
  //  | |    ___  __ _ ___| |_/ ___|  __ _ _   _  __ _ _ __ ___  ___
  //  | |   / _ \/ _` / __| __\___ \ / _` | | | |/ _` | '__/ _ \/ __|
  //  | |__|  __/ (_| \__ \ |_ ___) | (_| | |_| | (_| | | |  __/\__ \
  //  |_____\___|\__,_|___/\__|____/ \__, |\__,_|\__,_|_|  \___||___/
  //                                    |_|
  */
  template <typename T>
  integer
  LeastSquares<T>::factorize( integer NR,
                              integer NC,
                              valueType const A[], integer LDA,
                              valueType lambda ) {
    /*
    // calcolo fattorizzazione QR della matrice
    //
    //  / A        \ ( P )   / Q1 Q2 \ / R \
    //  |          |       = |       | |   |
    //  \ lambda I /         \ Q3 Q4 / \ 0 /
    //
    //  A = NR x NC
    //  I = NC x NC
    //  R = NC x NC
    */
    //ASSERT( lambda > 0, "LeastSquares::factorize( NR = " << NR <<
    //                    ", NC = " << NC << ", A, LDA = " << LDA <<
    //                    ", lambda = " << lambda << " ) lambda must be positive" ) ;
    nRow = NR+NC ;
    nCol = NC ;
    integer N = max(nRow,nCol) ;
    Lwork = 2*N+(N+1)*N ;
    allocReals.allocate(size_t(nRow*nCol+Lwork+N+nRow)) ;
    allocIntegers.allocate(size_t(N)) ;
    Amat = allocReals(size_t(nRow*nCol)) ;
    Work = allocReals(size_t(Lwork)) ;
    Tau  = allocReals(size_t(N)) ;
    tmp  = allocReals(size_t(nRow)) ;
    JPVT = allocIntegers(size_t(N)) ;
    integer info = gecopy( NR, NC, A, LDA, Amat, nRow ) ;
    if ( info == 0 ) {
      geid( NC, NC, Amat + NR, nRow, lambda ) ;
      info = geqp3( nRow, nCol, Amat, nRow, JPVT, Tau, Work, Lwork ) ;
    }
    return info ;
  }

  template <typename T>
  void
  LeastSquares<T>::solve( valueType const in[], valueType out[] ) const {
    /*
    // calcolo soluzione ai minimi quadrati di
    //
    //  || Ax - b ||^2 + lambda^2 ||x||^2
    //
    //  / A        \       / Q1 Q2 \ / R \ (P^T) ( x )     / b \
    //  |          | (x) = |       | |   |              ~= |   |
    //  \ lambda I /       \ Q3 Q4 / \ 0 /                 \ 0 /
    //
    //
    //  ( x ) = ( P R^(-1) 0 ) / Q1 Q2 \^T  / b \
    //                         |       |    |   |
    //                         \ Q3 Q4 /    \ 0 /
    //
    */
    // copio `in` in vettore temporaneo (b 0)
    copy( nRow-nCol, in, 1, tmp, 1 ) ;
    zero( nCol, tmp + nCol, 1 ) ;
    // moltiplico per Q
    integer info = ormqr( LEFT,
                          TRANSPOSE,
                          nRow, 1, // dimensione tmp
                          nCol,  // numero riflettori Q
                          Amat, nRow,
                          Tau,
                          tmp, nRow,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "LeastSquares::solve, ormqr return info = " << info ) ;
    // risolvo R
    trsv( UPPER,
          NO_TRANSPOSE,
          NON_UNIT,
          nCol,
          Amat, nRow,
          tmp, 1 ) ;
    // applico permutazione
    info = swaps( 1, tmp, nRow, 0, nRow-1, JPVT, 1 ) ;
    ALGLIN_ASSERT( info == 0, "LeastSquares::solve, swaps return info = " << info ) ;
    // copio soluzione
    copy( nCol, tmp, 1, out, 1 ) ;
  }

  template class LU<real> ;
  template class LU<doublereal> ;

  template class QR<real> ;
  template class QR<doublereal> ;

  template class QRP<real> ;
  template class QRP<doublereal> ;

  template class SVD<real> ;
  template class SVD<doublereal> ;

  template class TridiagonalLU<real> ;
  template class TridiagonalLU<doublereal> ;

  template class TridiagonalQR<real> ;
  template class TridiagonalQR<doublereal> ;

  template class LeastSquares<real> ;
  template class LeastSquares<doublereal> ;

} // end namespace alglin

///
/// eof: Alglin++.cc
///

