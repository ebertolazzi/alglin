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
  /*
  //   _    _   _
  //  | |  | | | |
  //  | |  | | | |
  //  | |__| |_| |
  //  |_____\___/
  */
  template <typename T>
  class LU {
  public:
    typedef T valueType ;
  
  private:

    integer     nRow, nCol, Lwork, __padding ;
    valueType * Amat ;
    valueType * Work ;
    integer   * Iwork ;
    integer   * i_pivot ;
  
    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;
    
    void check_ls( char const who[] ) const ;

  public:

    LU() ;
    ~LU() ;

    integer numRow() const { return nRow ; }
    integer numCol() const { return nCol ; }

    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    void
    t_factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    void
    solve( valueType xb[] ) const ;

    void
    t_solve( valueType xb[] ) const ;

    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const ;

    valueType cond1( valueType norm1 ) const ;
    valueType condInf( valueType normInf ) const ;
    void print( ostream & stream, integer prec = 4 ) const ;
  } ;

  //============================================================================
  /*
  //    ___  ____
  //   / _ \|  _ \
  //  | | | | |_) |
  //  | |_| |  _ <
  //   \__\_\_| \_\
  */
  template <typename T>
  class QR {
  public:
    typedef T valueType ;
  
  private:

    Malloc<valueType> allocReals ;

  protected:

    valueType * Amat ;
    valueType * Work ;
    valueType * Tau ;

    integer nRow, nCol, nReflector, Lwork ;

    void allocate( integer nr, integer nc ) ;

  public:

    QR()
    : allocReals("QR-allocReals")
    , nRow(0)
    , nCol(0)
    , nReflector(0)
    , Lwork(0)
    {}

    ~QR() {
      allocReals.free() ;
    }

    integer numRow() const { return nRow ; }
    integer numCol() const { return nCol ; }

    /*!
      Do QR factorization of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    /*!
      Do QR factorization of the transpose of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    t_factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ QR x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    void
    solve( valueType xb[] ) const ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ (QR)^T x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    void
    t_solve( valueType xb[] ) const ;

    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const ;

    /*\
     *  overwrites the general real M-by-N matrix C with
     *
     *                  SIDE = 'L'     SIDE = 'R'
     *  TRANS = 'N':      Q * C          C * Q
     *  TRANS = 'T':      Q**T * C       C * Q**T
    \*/
    void
    applyQ( SideMultiply  SIDE,
            Transposition TRANS,
            integer       nRefl,
            integer       nr,
            integer       nc,
            valueType     C[],
            integer       ldC ) const ;
    //! x <- Q*x
    void
    Q_mul( valueType x[] ) const
    { applyQ( LEFT, NO_TRANSPOSE, nReflector, nRow, 1, x, nRow ) ; }
    
    //! x <- Q'*x
    void
    Qt_mul( valueType x[] ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, nRow, 1, x, nRow ) ; }

    //! C <- Q*C
    void
    Q_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, NO_TRANSPOSE, nReflector, nr, nc, C, ldC ) ; }

    //! C <- Q'*C
    void
    Qt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, nr, nc, C, ldC ) ; }
    
    //! C <- C*Q
    void
    mul_Q( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, NO_TRANSPOSE, nReflector, nr, nc, C, ldC ) ; }
    
    //! C <- C*Q'
    void
    mul_Qt( integer nr, integer nc, valueType C[], integer ldC ) const
    { applyQ( RIGHT, TRANSPOSE, nReflector, nr, nc, C, ldC ) ; }

    // -------------------------------------------------------------------------

    void
    Rsolve( Transposition TRANS,
            integer       rk,
            valueType     x[],
            integer       incx ) const {
      trsv( UPPER, TRANS, NON_UNIT, rk, Amat, nRow, x, incx ) ;
    }

    void
    Rsolve( SideMultiply  SIDE,
            Transposition TRANS,
            integer       nr,
            integer       nc,
            valueType     alpha,
            valueType     Bmat[],
            integer       ldB ) const {
      trsm( SIDE, UPPER, TRANS, NON_UNIT, nr, nc, alpha, Amat, nRow, Bmat, ldB ) ;
    }
    
    //! x <- R^(-1) * x
    void
    invR_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, NO_TRANSPOSE, NON_UNIT, nReflector, Amat, nRow, x, incx ) ; }

    //! x <- R^(-T) * x
    void
    invRt_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, TRANSPOSE, NON_UNIT, nReflector, Amat, nRow, x, incx ) ; }

    //! C <- R^(-1) * C
    void
    invR_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ) ; }

    //! C <- R^(-T) * C
    void
    invRt_mul( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( LEFT, TRANSPOSE, nr, nc, 1.0, C, ldC ) ; }

    //! C <- C * R^(-1)
    void
    mul_invR( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, NO_TRANSPOSE, nr, nc, 1.0, C, ldC ) ; }

    void
    mul_invRt( integer nr, integer nc, valueType C[], integer ldC ) const
    { Rsolve( RIGHT, TRANSPOSE,  nr, nc, 1.0, C, ldC ) ; }
    
    void
    getR( valueType R[], integer ldR ) const ;

  } ;

  /*
  //    ___  ____  ____
  //   / _ \|  _ \|  _ \
  //  | | | | |_) | |_) |
  //  | |_| |  _ <|  __/
  //   \__\_\_| \_\_|
  */
  template <typename T>
  class QRP : public QR<T> {
  public:
    typedef typename QR<T>::valueType valueType ;

  private:
    Malloc<integer> allocIntegers ;
    integer * JPVT ;

  public:

    QRP()
    : allocIntegers("QRP-allocIntegers")
    {}

    ~QRP() {
      allocIntegers.free() ;
    }

    integer numRow() const { return this->nRow ; }
    integer numCol() const { return this->nCol ; }

    /*!
      Do QR factorization with column pivoting of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    /*!
      Do QR factorization with column pivoting of the transpose of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    t_factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ QR x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    void
    solve( valueType xb[] ) const ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ (QR)^T x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    void
    t_solve( valueType xb[] ) const ;

    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const ;

    // -------------------------------------------------------------------------
    void
    permute( valueType x[] ) const ;

    void
    inv_permute( valueType x[] ) const ;
    
    void
    permute_rows( integer nr, integer nc, valueType C[], integer ldC ) const {
      for ( integer j = 0 ; j < nc ; ++j ) permute( C + ldC*j ) ;
    }
    
    void
    inv_permute_rows( integer nr, integer nc, valueType C[], integer ldC ) const {
      for ( integer j = 0 ; j < nc ; ++j ) inv_permute( C + ldC*j ) ;
    }
  } ;

  //============================================================================
  /*
  //   ______     ______
  //  / ___\ \   / /  _ \
  //  \___ \\ \ / /| | | |
  //   ___) |\ V / | |_| |
  //  |____/  \_/  |____/
  */
  template <typename T>
  class SVD {
  public:
    typedef T valueType ;
  
  protected:

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    valueType * Amat ;
    valueType * Work ;
    valueType * Umat ;
    valueType * VTmat ;
    valueType * Svec ;
    integer   * IWork ;

    integer     nRow, nCol, minRC, Lwork ;
    bool        use_gesvd ;
    
    void allocate( integer NR, integer NC, integer LDA ) ;

  public:

    SVD()
    : allocReals("SVD-allocReals")
    , allocIntegers("SVD-allocIntegers")
    , nRow(0)
    , nCol(0)
    , Lwork(0)
    , use_gesvd(true)
    {}

    ~SVD() { allocReals.free() ; }

    integer numRow() const { return nRow ; }
    integer numCol() const { return nCol ; }

    /*!
      Do SVD factorization of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    valueType U( integer i, integer j ) const { return Umat[i+j*nRow] ; }
    valueType V( integer i, integer j ) const { return VTmat[j+i*nCol] ; }
    valueType sigma( integer i ) const { return Svec[i] ; }

    void
    solve( valueType const in[], valueType out[] ) const ;

    void
    solve( integer         nrhs,
           valueType const rhs[], integer ldRHS,
           valueType       x[],   integer ldX ) const ;

    void
    t_solve( valueType const in[], valueType out[] ) const ;

    void
    t_solve( integer         nrhs,
             valueType const rhs[], integer ldRHS,
             valueType       x[],   integer ldX ) const ;

    //! y <- alpha * U * x + beta * y
    void
    U_mul( valueType alpha, valueType const x[], integer incx,
           valueType beta,  valueType       y[], integer incy ) const {
      gemv( NO_TRANSPOSE,
            nRow, minRC,
            alpha, Umat, nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * U' * x + beta * y
    void
    Ut_mul( valueType alpha, valueType const x[], integer incx,
            valueType beta,  valueType       y[], integer incy ) const {
      gemv( TRANSPOSE,
            nRow, minRC,
            alpha, Umat, nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * V * x + beta * y
    void
    V_mul( valueType alpha, valueType const x[], integer incx,
           valueType beta,  valueType       y[], integer incy ) const {
      gemv( TRANSPOSE,
            minRC, nCol,
            alpha, VTmat, nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * V' * x + beta * y
    void
    Vt_mul( valueType alpha, valueType const x[], integer incx,
            valueType beta,  valueType       y[], integer incy ) const {
      gemv( NO_TRANSPOSE,
            minRC, nCol,
            alpha, VTmat, nRow,
            x, incx,
            beta, y, incy ) ;
    }

  } ;

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
               valueType const _U[] ) ;

    valueType cond1( valueType norm1 ) const ;
    valueType condInf( valueType normInf ) const ;

    void solve( valueType xb[] ) const ;
    void t_solve( valueType xb[] ) const ;
    void solve( integer nrhs, valueType xb[], integer ldXB ) const ;
    void t_solve( integer nrhs, valueType xb[], integer ldXB ) const ;

    void
    axpy( integer         N,
          valueType       alpha,
          valueType const L[],
          valueType const D[],
          valueType const U[],
          valueType const x[],
          valueType       beta,
          valueType       y[] ) const ;
  } ;

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

    void solve( valueType xb[] ) const ;
    void t_solve( valueType xb[] ) const ;
    void solve( integer nrhs, valueType xb[], integer ldXB ) const ;
    void t_solve( integer nrhs, valueType xb[], integer ldXB ) const ;

    void
    axpy( integer         N,
          valueType       alpha,
          valueType const L[],
          valueType const D[],
          valueType const U[],
          valueType const x[],
          valueType       beta,
          valueType       y[] ) const ;

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

