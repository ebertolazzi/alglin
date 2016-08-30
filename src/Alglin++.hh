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

  /*
  //   __  __       _ _
  //  |  \/  | __ _| | | ___   ___
  //  | |\/| |/ _` | | |/ _ \ / __|
  //  | |  | | (_| | | | (_) | (__
  //  |_|  |_|\__,_|_|_|\___/ \___|
  */

  //! Allocate memory
  template <typename T>
  class Malloc {
  public:
    typedef T                 valueType           ;
    typedef valueType*        valuePointer        ;
    typedef const valueType*  valueConstPointer   ;
    typedef valueType&        valueReference      ;
    typedef const valueType&  valueConstReference ;

  private:

    std::string  _name ;
    size_t       numTotValues ;
    size_t       numTotReserved ;
    size_t       numAllocated ;
    valuePointer pMalloc ;

    Malloc(Malloc<T> const &) ; // blocco costruttore di copia
    Malloc<T> const & operator = (Malloc<T> &) const ; // blocco copia

  public:

    //! malloc object constructor
    explicit Malloc( std::string const & __name )
    : _name(__name)
    , numTotValues(0)
    , numTotReserved(0)
    , numAllocated(0)
    , pMalloc(nullptr)
    { }

    //! malloc object destructor
    ~Malloc() { free() ; }

    //! allocate memory for `n` objects
    void
    allocate( size_t n ) {
      try {
        if ( n > numTotReserved ) {
          delete [] pMalloc ;
          numTotValues   = n ;
          numTotReserved = n + (n>>3) ; // 12% more values
          pMalloc = new T[numTotReserved] ;
        }
      }
      catch ( exception const & exc ) {
        cerr << "Memory allocation failed: " << exc.what()
             << "\nTry to allocate " << n << " bytes for " << _name
             << '\n' ;
        exit(0) ;
      }
      catch (...) {
        cerr << "Malloc allocation failed for " << _name << ": memory exausted\n" ;
        exit(0) ;
      }
      numTotValues = n ;
      numAllocated = 0 ;
    }

    //! free memory
    void
    free(void) {
      if ( pMalloc != nullptr ) {
        delete [] pMalloc ; pMalloc = nullptr;
        numTotValues   = 0 ;
        numTotReserved = 0 ;
        numAllocated   = 0 ;
      }
    }

    //! number of objects allocated
    size_t size(void) const { return numTotValues ; }

    //! get pointer of allocated memory for `sz` objets
    T * operator () ( size_t sz ) {
      size_t offs = numAllocated ;
      numAllocated += sz ;
      if ( numAllocated > numTotValues ) {
        cerr << "\nMalloc<" << _name << ">::operator () (" << sz << ") -- Malloc EXAUSTED\n" ;
        exit(0) ;
      }
      return pMalloc + offs ;
    }
  } ;

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
      allocReals.allocate(size_t(nRC*nRC+4*nRC)) ;
      allocIntegers.allocate(size_t(2*nRC)) ;
      this -> Amat    = allocReals(size_t(nRC*nRC)) ;
      this -> Work    = allocReals(size_t(4*nRC)) ;
      this -> i_pivot = allocIntegers(size_t(nRC)) ;
      this -> Iwork   = allocIntegers(size_t(nRC)) ;
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

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    integer     nRow, nCol, nReflector, Lwork ;
    valueType * Amat ;
    valueType * Work ;
    valueType * Tau ;
    valueType * TauTop ;
    integer   * JPVT ;
    valueType * xbtmp ;
    valueType   SVAL[3] ;
    valueType   mxabs ;
    integer     rank ;

    mutable bool top_factorized ;
    bool __pad[3] ;

    void topFactorize() const ;

  public:

    QR()
    : allocReals("allocReals")
    , allocIntegers("allocIntegers")
    , nRow(0)
    , nCol(0)
    , nReflector(0)
    , Lwork(0)
    , rank(-1)
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
    solve( valueType const rhs[], valueType x[] ) const ;

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

    Malloc<valueType> allocReals ;

    valueType * Amat ;
    valueType * Work ;
    valueType * Umat ;
    valueType * VTmat ;
    valueType * Svec ;
    valueType * invSvec ;

    valueType   mxabs ;

    integer     nRow, nCol, minRC, _rank ;
    integer     Lwork ;
    bool        use_gesvd ;

  public:

    SVD()
    : allocReals("allocReals")
    , mxabs(0)
    , nRow(0)
    , nCol(0)
    , _rank(0)
    , Lwork(0)
    , use_gesvd(true)
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
               valueType const _L[],
               valueType const _D[],
               valueType const _U[] ) {
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
      copy( N, _L, 1, L, 1 ) ;
      copy( N, _D, 1, D, 1 ) ;
      copy( N, _U, 1, U, 1 ) ;
      integer info = gttrf( N, L, D, U, U2, IPIV ) ;
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

    Malloc<valueType> allocReals ;

    valueType * C   ; // rotazioni givens
    valueType * S   ;
    valueType * BD  ; // band triangular matrix
    valueType * BU  ; // band triangular matrix
    valueType * BU2 ; // band triangular matrix

    valueType   normInfA ;
    integer     nRC ;

    void Rsolve( valueType xb[] ) const ;
    void RsolveTransposed( valueType xb[] ) const ;

  public:

    TridiagonalQR()
    : allocReals("allocReals")
    , nRC(0)
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

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    valueType * Amat ;
    valueType * Work ;
    valueType * Tau ;
    valueType * tmp ;
    integer   * JPVT ;
    valueType   mxabs ;

    integer     nRow, nCol, Lwork ;

  public:

    LeastSquares()
    : allocReals("allocReals")
    , allocIntegers("allocIntegers")
    , mxabs(0)
    , nRow(0)
    , nCol(0)
    , Lwork(0)
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

