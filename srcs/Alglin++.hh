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
/// file: Alglin++.hh
///

#ifndef ALGLINPP_HH
#define ALGLINPP_HH

#include "Alglin.hh"

/*
//    ____            _       _             __
//   / ___| _     _  (_)_ __ | |_ ___ _ __ / _| __ _  ___ ___
//  | |   _| |_ _| |_| | '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
//  | |__|_   _|_   _| | | | | ||  __/ |  |  _| (_| | (_|  __/
//   \____||_|   |_| |_|_| |_|\__\___|_|  |_|  \__,_|\___\___|
*/

namespace alglin {

  //============================================================================

  template <typename T>
  class LU {
  public:
    typedef T valueType ;
  
  private:

    integer     nRC, Lwork ;
    valueType * Amat ;
    valueType * Work ;
    integer   * Iwork ;
    integer   * i_pivot ;
  
    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

  public:

    LU()
    : nRC(0)
    , Lwork(0)
    , allocReals("allocReals")
    , allocIntegers("allocIntegers")
    {}

    ~LU() {
      allocReals.free() ;
      allocIntegers.free() ;
    }

    integer
    factorize( integer N, valueType const A[], integer LDA ) {
      nRC = N ;
      allocReals.allocate(nRC*nRC+4*nRC) ;
      allocIntegers.allocate(2*nRC) ;
      this -> Amat    = allocReals(nRC*nRC) ;
      this -> Work    = allocReals(4*nRC) ;
      this -> i_pivot = allocIntegers(nRC) ;
      this -> Iwork   = allocIntegers(nRC) ;
      integer info = gecopy( nRC, nRC, A, LDA, Amat, nRC ) ;
      ALGLIN_ASSERT( info == 0, "LU::factorize gecopy INFO = " << info ) ;
      info = getrf( nRC, nRC, Amat, nRC, i_pivot ) ;
      ALGLIN_ASSERT( info >= 0, "LU::factorize getrf INFO = " << info ) ;
      return info ;
    }

    void
    solve( valueType xb[], Transposition const & TRANS = NO_TRANSPOSE ) const {
      integer info = getrs( TRANS, nRC, 1, Amat, nRC, i_pivot, xb, nRC ) ;
      ALGLIN_ASSERT( info >= 0, "LU::solve getrs INFO = " << info ) ;
    }

    void
    solve( integer   nrhs,
           valueType B[],
           integer   ldB,
           Transposition const & TRANS = NO_TRANSPOSE ) const {
      integer info = getrs( TRANS, nRC, nrhs, Amat, nRC, i_pivot, B, ldB ) ;
      ALGLIN_ASSERT( info >= 0, "LU::solve getrs INFO = " << info ) ;
    }

    valueType
    cond1( valueType norm1 ) const {
      valueType rcond ;
      integer info = gecon1( nRC, Amat, nRC, norm1, rcond, Work, Iwork ) ;
      ALGLIN_ASSERT( info == 0, "LU::cond1, gecon1 return info = " << info ) ;
      return rcond ;
    }

    valueType
    condInf( valueType normInf ) const {
      valueType rcond ;
      integer info = geconInf( nRC, Amat, nRC, normInf, rcond, Work, Iwork ) ;
      ALGLIN_ASSERT( info == 0, "LU::condInf, geconInf return info = " << info ) ;
      return rcond ;
    }

    void print( ostream & stream, integer prec = 4 ) const ;

  } ;

  //============================================================================

  template <typename T>
  class QR {
  public:
    typedef T valueType ;
  
  private:

    integer     rank ;
    integer     nRow, nCol, nReflector, Lwork ;
    valueType * Amat ;
    valueType * Work ;
    valueType * Tau ;
    valueType * TauTop ;
    integer   * JPVT ;
    valueType * xbtmp ;
    valueType   SVAL[3] ;
    valueType   mxabs ;
  
    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    mutable bool top_factorized ;

    void topFactorize() const ;

  public:

    QR()
    : rank(-1)
    , nRow(0)
    , nCol(0)
    , nReflector(0)
    , Lwork(0)
    , allocReals("allocReals")
    , allocIntegers("allocIntegers")
    , top_factorized(false)
    {}

    ~QR() {
      allocReals.free() ;
      allocIntegers.free() ;
    }

    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    void
    evalNumericRank( valueType rcond ) {
      integer info = rankEstimate<valueType>( nRow, nCol, Amat, nRow, rcond, rank, SVAL ) ;
      ALGLIN_ASSERT( info == 0, "QR::evalNumericRank, return info = " << info ) ;
    }

    integer numericRank() const { return rank ; }

    void
    minq( valueType const rhs[], valueType x[], valueType lambda = 0 ) const ;

    void
    minq( integer nrhs,
          valueType const rhs[], integer ldRHS,
          valueType       x[],   integer ldX,
          valueType lambda = 0 ) const ;

    /*\
     *  overwrites the general real M-by-N matrix C with
     *
     *                  SIDE = 'L'     SIDE = 'R'
     *  TRANS = 'N':      Q * C          C * Q
     *  TRANS = 'T':      Q**T * C       C * Q**T
    \*/
    void
    applyQ( SideMultiply  const & SIDE,
            Transposition const & TRANS,
            integer rank,
            integer ncol,
            valueType C[], integer ldC ) const ;

    void
    solveT( Transposition const & TRANS,
            integer rank,
            valueType x[], integer incx ) const ;
    
    void
    solveT( SideMultiply  const & SIDE,
            Transposition const & TRANS,
            integer rank, integer ncol,
            valueType alpha,
            valueType B[], integer ldB ) const ;

    void permute( valueType x[] ) const ;
    void permute( integer ncolC, valueType C[], integer ldC ) const ;

    void print( ostream & stream, integer prec = 4 ) const ;

  } ;

  //============================================================================

  template <typename T>
  class SVD {
  public:
    typedef T valueType ;
  
  private:
    integer     nRow, nCol, minRC, _rank, Lwork ;
    valueType * Amat ;
    valueType * Work ;
    valueType * Umat ;
    valueType * VTmat ;
    valueType * Svec ;
    valueType * invSvec ;
    valueType   mxabs ;
    bool        use_gesvd ;

    Malloc<valueType> allocReals ;

  public:

    SVD()
    : nRow(0)
    , nCol(0)
    , _rank(0)
    , Lwork(0)
    , mxabs(0)
    , use_gesvd(true)
    , allocReals("allocReals")
    {}

    ~SVD() { allocReals.free() ; }

    void
    factorize( integer         NR,
               integer         NC,
               valueType const A[],
               integer         LDA,
               valueType       threshold ) ;

    integer numRow() const { return nRow ; }
    integer numCol() const { return nCol ; }
    integer rank()   const { return _rank ; }

    valueType U( integer i, integer j ) const { return Umat[i+j*nRow] ; }
    valueType V( integer i, integer j ) const { return VTmat[j+i*nCol] ; }
    valueType sigma( integer i ) const { return Svec[i] ; }

    void
    solve( valueType const in[], valueType out[] ) const ;

    void
    solve( integer         nrhs,
           valueType const rhs[], integer ldRHS,
           valueType       x[],   integer ldX ) const ;

  } ;

  //============================================================================

  template <typename T>
  class TridiagonalLU {
  public:
    typedef T valueType ;
  
  private:

    integer     nRC ;
    valueType * L ;
    valueType * D ;
    valueType * U ;
    valueType * U2 ;
    valueType * WORK ;
    integer   * IPIV ;
    integer   * IWORK ;
  
    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

  public:

    TridiagonalLU()
    : nRC(0)
    , allocReals("allocReals")
    , allocIntegers("allocIntegers")
    {}

    ~TridiagonalLU() {
      allocReals.free() ;
      allocIntegers.free() ;
    }

    void
    factorize( integer         N,
               valueType const L[],
               valueType const D[],
               valueType const U[] ) {
      this -> nRC = N ;
      allocReals.allocate(6*N) ;
      allocIntegers.allocate(2*N) ;
      this -> L     = allocReals(N) ;
      this -> D     = allocReals(N) ;
      this -> U     = allocReals(N) ;
      this -> U2    = allocReals(N) ;
      this -> WORK  = allocReals(2*N) ;
      this -> IPIV  = allocIntegers(N) ;
      this -> IWORK = allocIntegers(N) ;
      copy( N, L, 1, this->L, 1 ) ;
      copy( N, D, 1, this->D, 1 ) ;
      copy( N, U, 1, this->U, 1 ) ;
      integer info = gttrf( N, this->L, this->D, this->U, this->U2, IPIV ) ;
      ALGLIN_ASSERT( info == 0, "Tridiagonal::factorize, return info = " << info ) ;
    }

    valueType
    cond1( valueType norm1 ) const {
      valueType rcond ;
      integer info = gtcon1( nRC, L, D, U, U2, IPIV, norm1, rcond, WORK, IWORK ) ;
      ALGLIN_ASSERT( info == 0, "Tridiagonal::cond1, return info = " << info ) ;
      return rcond ;
    }

    valueType
    condInf( valueType normInf ) const {
      valueType rcond ;
      integer info = gtconInf( nRC, L, D, U, U2, IPIV, normInf, rcond, WORK, IWORK ) ;
      ALGLIN_ASSERT( info == 0, "Tridiagonal::cond1, return info = " << info ) ;
      return rcond ;
    }

    void
    solve( valueType xb[], Transposition const & TRANS = NO_TRANSPOSE ) const {
      integer info = gttrs( TRANS, nRC, 1, L, D, U, U2, IPIV, xb, nRC ) ;
      ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
    }

    void
    solve( integer nrhs, valueType xb[], integer ldXB, Transposition const & TRANS = NO_TRANSPOSE ) const {
      integer info = gttrs( TRANS, nRC, nrhs, L, D, U, U2, IPIV, xb, ldXB ) ;
      ALGLIN_ASSERT( info == 0, "Tridiagonal::solve, return info = " << info ) ;
    }
  } ;

  //============================================================================

  template <typename T>
  class TridiagonalQR {
  public:
    typedef T valueType ;
  
  private:

    integer     nRC ;
    valueType * C   ; // rotazioni givens
    valueType * S   ;
    valueType * BD  ; // band triangular matrix
    valueType * BU  ; // band triangular matrix
    valueType * BU2 ; // band triangular matrix
    
    valueType   normInfA ;
  
    Malloc<valueType> allocReals ;

    void Rsolve( valueType xb[] ) const ;
    void RsolveTransposed( valueType xb[] ) const ;

  public:

    TridiagonalQR()
    : nRC(0)
    , allocReals("allocReals")
    {}

    ~TridiagonalQR() {
      allocReals.free() ;
    }

    void
    factorize( integer         N,
               valueType const L[],
               valueType const D[],
               valueType const U[] ) ;

    void
    solve( valueType xb[],
           Transposition const & TRANS = NO_TRANSPOSE ) const ;

    void
    solve( integer nrhs, valueType xb[], integer ldXB,
           Transposition const & TRANS = NO_TRANSPOSE ) const ;

    void
    minq( integer nrhs,
          T       RHS[],
          integer ldRHS,
          T       lambda ) const ;

  } ;

  //============================================================================

  template <typename T>
  class LeastSquares {
  public:
    typedef T valueType ;

  private:
    integer     nRow, nCol, Lwork ;
    valueType * Amat ;
    valueType * Work ;
    valueType * Tau ;
    valueType * tmp ;
    integer   * JPVT ;
    valueType   mxabs ;

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

  public:

    LeastSquares()
    : nRow(0)
    , nCol(0)
    , Lwork(0)
    , mxabs(0)
    , allocReals("allocReals")
    , allocIntegers("allocIntegers")
    {}

    ~LeastSquares() {
      allocReals.free() ;
      allocIntegers.free() ;
    }

    integer
    factorize( integer NR, integer NC,
               valueType const A[], integer LDA,
               valueType lambda ) ;

    void
    solve( valueType const in[], valueType out[] ) const ;

  } ;

} // end namespace alglin

#endif

///
/// eof: Alglin++.hh
///

