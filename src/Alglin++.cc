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
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
  */

  template <typename T>
  void
  QR<T>::factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
    mxabs = maxabs( NR, NC, A, LDA ) ;
    ALGLIN_ASSERT( mxabs > 0, "QR::factorize, mxabs = " << mxabs ) ;

    // calcolo fattorizzazione QR della matrice A
    nRow       = NR ;
    nCol       = NC ;
    nReflector = min(nRow,nCol) ;
    integer mm = max(nRow,nCol) ;
    valueType tmp ; // get optimal allocation
    integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::geqp3 return info = " << info ) ;
    Lwork = integer(tmp) ;
    info = tzrzf( NR, NC, nullptr, NR, nullptr, &tmp, -1 ) ;
    ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::tzrzf return info = " << info ) ;
    if ( Lwork < integer(tmp) ) Lwork = integer(tmp) ;
    if ( Lwork < nRow*nCol    ) Lwork = nRow*nCol ;
    allocReals.allocate(nRow*nCol+Lwork+nCol+2*nReflector+mm) ;
    allocIntegers.allocate(nCol) ;
    Amat   = allocReals(nRow*nCol) ;
    Work   = allocReals(Lwork) ;
    Tau    = allocReals(nReflector) ;
    TauTop = allocReals(nReflector) ;
    JPVT   = allocIntegers(nCol) ;
    xbtmp  = allocReals(mm) ;
    info   = gecopy( NR, NC, A, LDA, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::gecopy return info = " << info ) ;
    info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nRow, nCol, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::lascl return info = " << info ) ;
    std::fill( JPVT, JPVT + nCol, 0 ) ;
    ///outMAPLE( nRow, nCol, Amat, nRow, cout ) ;
    info = geqp3( nRow, nCol, Amat, nRow, JPVT, Tau, Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::geqp3 return info = " << info ) ;
    top_factorized = false ;
    rank           = nReflector ;
  }

  template <typename T>
  void
  QR<T>::applyQ( SideMultiply  const & SIDE,
                 Transposition const & TRANS,
                 integer NR, integer NC,
                 valueType C[], integer ldC ) const {
    ALGLIN_ASSERT( (SIDE == alglin::LEFT  && NR == nRow) ||
                   (SIDE == alglin::RIGHT && NC == nRow),
                   "QR::applyQ NR = " << NR << " NC = " << NC << " nRow = " << nRow ) ;
    integer info = ormqr( SIDE, TRANS,
                          NR, NC,
                          nReflector,  // numero riflettori Q
                          Amat, nRow /*ldA*/,
                          Tau,
                          C, ldC,
                          Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "QR::applyQ call alglin::ormqr return info = " << info ) ;
  }
    
  template <typename T>
  void
  QR<T>::solveT( Transposition const & TRANS,
                 integer rk, valueType x[], integer incx ) const {
    // risolvo R
    trsv( UPPER, TRANS, NON_UNIT, rk, Amat, nRow, x, incx ) ;
  }
    
  template <typename T>
  void
  QR<T>::solveT( SideMultiply  const & SIDE,
                 Transposition const & TRANS,
                 integer N, integer M,
                 valueType alpha, valueType Bmat[], integer ldB ) const {
    // risolvo R
    trsm( SIDE, UPPER, TRANS, NON_UNIT, N, M, alpha, Amat, nRow, Bmat, ldB ) ;
  }

  template <typename T>
  void
  QR<T>::permute( valueType x[] ) const {
    // applico permutazione
    for ( integer i = 0 ; i < nCol ; ++i ) Work[JPVT[i]-1] = x[i] ;
    copy( nCol, Work, 1, x, 1 ) ;
  }

  template <typename T>
  void
  QR<T>::permute( integer ncolC, valueType C[], integer ldC ) const {
    // applico permutazione
    for ( integer j = 0 ; j < ncolC ; ++j ) permute( C + ldC*j ) ;
  }

  template <typename T>
  void
  QR<T>::topFactorize() const {
    integer info = tzrzf( rank, nCol, Amat, nRow, TauTop, Work, Lwork ) ;
    ALGLIN_ASSERT( info == 0, "QR::topFactorize call alglin::tzrzf return info = " << info ) ;
    top_factorized = true ;
  }

  template <typename T>
  void
  QR<T>::solve( valueType const rhs[], valueType x[] ) const {
    copy( nRow, rhs, 1, xbtmp, 1 ) ;
    applyQ( LEFT, TRANSPOSE, nRow, 1, xbtmp, nRow ) ;
    solveT( NO_TRANSPOSE, rank, xbtmp, 1 ) ;
    for ( integer i = 0 ; i < nCol ; ++i ) x[JPVT[i]-1] = xbtmp[i] ;
    integer info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nCol, 1, x, nCol ) ;
    ALGLIN_ASSERT( info == 0, "QR::minq call alglin::lascl return info = " << info ) ;
  }

  template <typename T>
  void
  QR<T>::minq( valueType const rhs[], valueType x[], valueType lambda ) const {

    // se NR > NC |x| < |rhs|
    copy( nRow, rhs, 1, xbtmp, 1 ) ;
    applyQ( LEFT, TRANSPOSE, nRow, 1, xbtmp, nRow ) ;

    if ( rank < nCol && !top_factorized ) topFactorize() ;

    if ( lambda > 0 ) triTikhonov( rank, Amat, nRow, 1, xbtmp, rank, lambda ) ;
    else              solveT( NO_TRANSPOSE, rank, xbtmp, 1 ) ;

    if ( rank < nCol ) {
      for ( integer i = rank ; i < nCol ; ++i ) xbtmp[i] = 0 ;
      integer info = ormrz( LEFT, TRANSPOSE,
                            nCol, 1, rank, nCol-rank,
                            Amat, nRow,
                            TauTop,
                            xbtmp, nCol,
                            Work, Lwork ) ;
      ALGLIN_ASSERT( info == 0, "QR::minq call alglin::ormrz return info = " << info ) ;
    }
    for ( integer i = 0 ; i < nCol ; ++i ) x[JPVT[i]-1] = xbtmp[i] ;

    integer info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nCol, 1, x, nCol ) ;
    ALGLIN_ASSERT( info == 0, "QR::minq call alglin::lascl return info = " << info ) ;

  }

  template <typename T>
  void
  QR<T>::minq( integer nrhs,
               valueType const rhs[], integer ldRHS,
               valueType       x[],   integer ldX,
               valueType lambda ) const {
    ALGLIN_ASSERT( ldRHS >= nRow, "QR::minq ldRHS = " << ldRHS << " < " << nRow << " = nRow" ) ;
    ALGLIN_ASSERT( ldX   >= nCol, "QR::minq ldX = " << ldX << " < " << nCol << " = nCol" ) ;

    // se NR > NC |x| < |rhs|
    for ( integer j = 0 ; j < nrhs ; ++j ) {
      copy( nRow, rhs+j*ldRHS, 1, xbtmp, 1 ) ;
      applyQ( LEFT, TRANSPOSE, nRow, 1, xbtmp, nRow ) ;
      copy( nRow, xbtmp, 1, x+j*ldX, 1 ) ;
    }

    if ( rank < nCol && !top_factorized ) topFactorize() ;

    if ( lambda > 0 ) triTikhonov( rank, Amat, nRow, nrhs, x, rank, lambda ) ;
    else              solveT( LEFT, NO_TRANSPOSE, rank, nrhs, 1.0, x, ldX ) ;

    if ( rank < nCol ) {
      for ( integer j = 0 ; j < nrhs ; ++j ) {
        valueType * xj = x+j*ldX ;
        copy( rank, xj, 1, xbtmp, 1 ) ;
        zero( nCol-rank, xbtmp+rank, 1 ) ;
        integer info = ormrz( LEFT, TRANSPOSE,
                              nCol, 1, rank, nCol-rank,
                              Amat, nRow,
                              TauTop,
                              xbtmp, nCol,
                              Work, Lwork ) ;
        ALGLIN_ASSERT( info == 0, "QR::minq call alglin::ormrz return info = " << info ) ;
        for ( integer i = 0 ; i < nCol ; ++i ) xj[JPVT[i]-1] = xbtmp[i] ;
      }
    }

    integer info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nCol, nrhs, x, ldX ) ;
    ALGLIN_ASSERT( info == 0, "QR::minq call alglin::lascl return info = " << info ) ;

  }

  template <typename T>
  void
  QR<T>::print( ostream & stream, integer prec ) const {
    stream << "R " << nRow << " x " << nCol << " =\n" ;
    outMATRIX( UPPER_TRIANGULAR_MATRIX, nRow, nCol, Amat, nRow, stream, prec ) ;
    stream << "A " << nRow << " x " << nCol << " =\n" ;
    outMATRIX( FULL_MATRIX, nRow, nCol, Amat, nRow, stream, prec ) ;
    stream << "TAU 1 x " << nReflector << " =\n" ;
    outMATRIX( FULL_MATRIX, 1, nReflector, Tau, 1, stream, prec ) ;
    stream << "JPVT 1 x " << nReflector << " =\n" ;
    outMATRIX( FULL_MATRIX, 1, nReflector, JPVT, 1, stream, prec ) ;
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
  SVD<T>::factorize( integer NR,
                     integer NC,
                     T const A[],
                     integer LDA,
                     T       threshold ) {
    mxabs = maxabs( NR, NC, A, LDA ) ;
    ALGLIN_ASSERT( mxabs > 0, "SVD::factorize, mxabs = " << mxabs ) ;
    nRow  = NR ;
    nCol  = NC ;
    minRC = min(NR,NC) ;
    valueType tmp ;
    integer info = use_gesvd ? gesvd( REDUCED,
                                      REDUCED,
                                      NR, NC,
                                      nullptr, LDA,
                                      nullptr,
                                      nullptr, NR,
                                      nullptr, minRC,
                                      &tmp, -1 ) :
                               gesdd( REDUCED,
                                      NR, NC,
                                      nullptr, LDA,
                                      nullptr,
                                      nullptr, NR,
                                      nullptr, minRC,
                                      &tmp, -1, nullptr ) ;
    Lwork = integer(tmp) ;
    allocReals.allocate(nRow*nCol+minRC*(nRow+nCol+2)+Lwork) ;
    Amat    = allocReals(nRow*nCol) ;
    Svec    = allocReals(minRC) ;
    invSvec = allocReals(minRC) ;
    Umat    = allocReals(minRC*nRow) ;
    VTmat   = allocReals(minRC*nCol) ;
    Work    = allocReals(Lwork) ;
    info    = gecopy( NR, NC, A, LDA, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gecopy return info = " << info ) ;
    info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nRow, nCol, Amat, nRow ) ;
    ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::lascl return info = " << info ) ;
    if ( use_gesvd ) {
      info = gesvd( REDUCED,
                    REDUCED,
                    NR, NC, Amat, NR,
                    Svec,
                    Umat, NR,
                    VTmat, minRC,
                    Work, Lwork ) ;
      ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gesvd return info = " << info ) ;
    } else {
      std::vector<integer> IWork(8*minRC) ;
      info = gesdd( REDUCED,
                    NR, NC, Amat, NR,
                    Svec,
                    Umat, NR,
                    VTmat, minRC,
                    Work, Lwork, &IWork.front() ) ;
      ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gesdd return info = " << info ) ;
    }
    tmp = threshold*Svec[0] ;
    for ( _rank = 0 ; _rank < minRC && Svec[_rank] > tmp; ++_rank ) invSvec[_rank] = 1/Svec[_rank] ;
    for ( integer i = _rank+1 ; i < minRC ; ++i ) invSvec[i] = 2/hypot(tmp,Svec[i]) ;
  }

  template <typename T>
  void
  SVD<T>::solve( valueType const in[], valueType out[] ) const {
    // A = U*S*VT
    // U*S*VT*x=b --> VT^T S^+ U^T b
    // U  nRow x minRC
    // VT minRC x nCol
    gemv( TRANSPOSE,
          nRow,
          minRC,
          1.0, Umat, nRow,
          in, 1,
          0.0, Work, 1 ) ;

    for ( integer i = 0 ; i < minRC ; ++i ) Work[i] *= invSvec[i] ;

    gemv( TRANSPOSE,
          minRC,
          nCol,
          1.0, VTmat, minRC,
          Work, 1,
          0.0, out, 1 ) ;

    integer info = lascl( FULL_MATRIX, 0, 0, mxabs, 1, nCol, 1, out, nCol ) ;
    ALGLIN_ASSERT( info == 0, "SVD::solve call alglin::lascl return info = " << info ) ;
  }

  template <typename T>
  void
  SVD<T>::solve( integer         nrhs,
                 valueType const rhs[], integer ldRHS,
                 valueType       x[],   integer ldX ) const {
    for ( integer i = 0 ; i < nrhs ; ++i ) solve( rhs + i*ldRHS, x + i*ldX ) ; 
  }

  /*
  //   _____     _     _ _                               _
  //  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |
  //    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | |
  //    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | |
  //    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|
  //                             |___/
  */

  template <typename T>
  void
  TridiagonalQR<T>::factorize( integer         N,
                               valueType const L[],
                               valueType const D[],
                               valueType const U[] ) {
    allocReals.allocate(5*(N-1)) ;
    this -> nRC = N ;
    this -> C   = allocReals(N-1) ;
    this -> S   = allocReals(N-1) ;
    this -> BD  = allocReals(N) ;
    this -> BU  = allocReals(N-1) ;
    this -> BU2 = allocReals(N-2) ;

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
  TridiagonalQR<T>::solve( valueType xb[], Transposition const & TRANS ) const {
    if ( TRANS == NO_TRANSPOSE ) {
      // A x = b --> Q A x = Q b --> R x = Q b
      // applico Q b
      for ( integer i = 0 ; i < nRC-1 ; ++i )
        rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], S[i] ) ;
      Rsolve( xb ) ;
    } else {
      // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
      RsolveTransposed( xb ) ;
      // applico Q^T b
      for ( integer i = nRC-2 ; i >= 0 ; --i )
        rot( 1, &xb[i], 1, &xb[i+1], 1, C[i], -S[i] ) ;
    }
  }

  template <typename T>
  void
  TridiagonalQR<T>::solve( integer nrhs, valueType xb[], integer ldXB,
                           Transposition const & TRANS ) const {
    if ( TRANS == NO_TRANSPOSE ) {
      // A x = b --> Q A x = Q b --> R x = Q b
      // applico Q b
      for ( integer i = 0 ; i < nRC-1 ; ++i ) rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], S[i] ) ;
      for ( integer i = 0 ; i < nrhs  ; ++i ) Rsolve( xb+i*ldXB ) ;
    } else {
      // A^T x = b --> A^T Q^T Q x = b --> R^T Q x = b --> R^T y = b  x = Q^T y
      for ( integer i = 0     ; i < nrhs ; ++i ) RsolveTransposed(xb+i*ldXB) ;
      for ( integer i = nRC-2 ; i >= 0   ; --i ) rot( nrhs, &xb[i], ldXB, &xb[i+1], ldXB, C[i], -S[i] ) ;
    }
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
  TridiagonalQR<T>::minq( integer nrhs,
                          T       RHS[],
                          integer ldRHS,
                          T       lambda_in) const {

    valueType lambda = normInfA * lambda_in ;
    std::vector<T> D(nRC), U(nRC-1), U2(nRC-2), tmp(nrhs) ;
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
      for ( integer i = nRC-3 ; i >= 0 ; --i )
        xb[i] = (xb[i]-U[i]*xb[i+1]-U2[i]*xb[i+2])/D[i] ;
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
    allocReals.allocate(nRow*nCol+Lwork+N+nRow) ;
    allocIntegers.allocate(N) ;
    Amat = allocReals(nRow*nCol) ;
    Work = allocReals(Lwork) ;
    Tau  = allocReals(N) ;
    tmp  = allocReals(nRow) ;
    JPVT = allocIntegers(N) ;
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

  template class QR<real> ;
  template class QR<doublereal> ;

  template class TridiagonalQR<real> ;
  template class TridiagonalQR<doublereal> ;

  template class SVD<real> ;
  template class SVD<doublereal> ;

  template class LeastSquares<real> ;
  template class LeastSquares<doublereal> ;

} // end namespace alglin

///
/// eof: Alglin++.cc
///

