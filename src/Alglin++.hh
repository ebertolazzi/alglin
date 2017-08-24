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

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wc++98-compat"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

/*\
 |    ____            _       _             __
 |   / ___| _     _  (_)_ __ | |_ ___ _ __ / _| __ _  ___ ___
 |  | |   _| |_ _| |_| | '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
 |  | |__|_   _|_   _| | | | | ||  __/ |  |  _| (_| | (_|  __/
 |   \____||_|   |_| |_|_| |_|\__\___|_|  |_|  \__,_|\___\___|
\*/

namespace alglin {

  /*\
   |   __  __       _ _
   |  |  \/  | __ _| | | ___   ___
   |  | |\/| |/ _` | | |/ _ \ / __|
   |  | |  | | (_| | | | (_) | (__
   |  |_|  |_|\__,_|_|_|\___/ \___|
  \*/

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
    explicit
    Malloc( std::string const & __name )
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
      catch ( std::exception const & exc ) {
        std::cerr << "Memory allocation failed: " << exc.what()
             << "\nTry to allocate " << n << " bytes for " << _name
             << '\n' ;
        std::exit(0) ;
      }
      catch (...) {
        std::cerr << "Malloc allocation failed for " << _name << ": memory exausted\n" ;
        std::exit(0) ;
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
        std::cerr << "\nMalloc<" << _name << ">::operator () (" << sz << ") -- Malloc EXAUSTED\n" ;
        std::exit(0) ;
      }
      return pMalloc + offs ;
    }
  } ;

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  template <typename T>
  integer
  rankEstimate( integer   M,
                integer   N,
                T         A[],
                integer   LDA,
                T         RCOND,
                T         SVAL[3] ) ;

  //! base class for linear systema solver
  template <typename T>
  class LinearSystemSolver {
  public:
    typedef T valueType ;

  public:

    LinearSystemSolver()
    {}

    virtual ~LinearSystemSolver()
    {}

    virtual
    void
    solve( valueType xb[] ) const ALGLIN_PURE_VIRTUAL ;

    virtual
    void
    t_solve( valueType xb[] ) const ALGLIN_PURE_VIRTUAL ;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0 ; i < nrhs ; ++i ) solve( B + i*ldB ) ;
    }

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const {
      for ( integer i = 0 ; i < nrhs ; ++i ) t_solve( B + i*ldB ) ;
    }

  } ;
  
  /*\
   |   _____          _             _          _   _
   |  |  ___|_ _  ___| |_ ___  _ __(_)______ _| |_(_) ___  _ __
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _` | __| |/ _ \| '_ \
   |  |  _| (_| | (__| || (_) | |  | |/ / (_| | |_| | (_) | | | |
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\__,_|\__|_|\___/|_| |_|
  \*/

  template <typename T>
  class Factorization : public LinearSystemSolver<T> {
  public:
    typedef T valueType ;
  
  protected:

    valueType * Amat ;
    integer     nRow ;
    integer     nCol ;

  public:

    Factorization()
    : nRow(0)
    , nCol(0)
    {}

    virtual
    ~Factorization()
    {}

    integer           numRow()   const { return nRow ; }
    integer           numCol()   const { return nCol ; }
    valueType *       Apointer()       { return Amat ; }
    valueType const * Apointer() const { return Amat ; }

    virtual
    void
    allocate( integer NR, integer NC ) ALGLIN_PURE_VIRTUAL;

    virtual
    void
    factorize() ALGLIN_PURE_VIRTUAL;

    virtual
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ALGLIN_PURE_VIRTUAL ;

    /*!
      Zeroes a rectangular block of the stored matrix staring at `(irow,icol)` position

      \param[in] nr    number of rows of the block to be zeroed
      \param[in] nc    number of columns of the block to be zeroed
      \param[in] irow  starting row
      \param[in] icol  stating column
    \*/
    void
    zero_block( integer nr,
                integer nc,
                integer irow,
                integer icol ) {
      valueType * Ablk = Amat + irow + icol * nRow ;
      gezero( nr, nc, Ablk, nRow ) ;
    }

    /*!
      Copy a matrix to a rectangular block of the stored matrix staring at `(irow,icol)` position

      \param[in] nr    number of rows of the block to be zeroed
      \param[in] nc    number of columns of the block to be zeroed
      \param[in] B     pointer to memory storing the input matrix `B`
      \param[in] ldB   leading dimension of the matrix `B`
      \param[in] irow  starting row
      \param[in] icol  stating column
    \*/
    void
    load_block( integer         nr,
                integer         nc,
                valueType const B[],
                integer         ldB,
                integer         irow = 0,
                integer         icol = 0 ) {
      valueType * Ablk = Amat + irow + icol * nRow ;
      integer info = gecopy( nr, nc, B, ldB, Ablk, nRow ) ;
      ALGLIN_ASSERT( info == 0, "block_load call alglin::gecopy return info = " << info ) ;
    }

    /*!
      Copy vector `column` to the `icol`th column of the internal stored matrix
      \param[in] column the column vector
      \param[in] icol   the column to be changed
    \*/
    void
    load_column( valueType const column[], integer icol ) {
      copy( nRow, column, 1, Amat + icol * nRow, 1 ) ;
    }

    /*!
      Copy vector `row` to the `irow`th row of the internal stored matrix
      \param[in] row  the row vector
      \param[in] irow the row to be changed
    \*/
    void
    load_row( valueType const row[], integer irow ) {
      copy( nRow, row, 1, Amat + irow, nRow ) ;
    }

    /*!
      Copy vector element of a sparse vector to a column of the internal stored matrix
      \param[in] nnz        number of nonzeros of the columns
      \param[in] values     the values of the sparse vector
      \param[in] row        index position of the values of the sparse vector
      \param[in] icol       the column to be changed
      \param[in] index_base offset for the index, 0 for C based vector 1 for FORTRAN based vector
    \*/
    void
    load_sparse_column( integer         nnz,
                        valueType const values[],
                        integer   const row[],
                        integer         icol,
                        integer         index_base = 0 ) {
      valueType * Acol = Amat + icol * nRow ;
      zero( nRow, Acol, 1 ) ;
      for ( integer i = 0 ; i < nnz ; ++i ) Acol[row[i]-index_base] = values[i] ;
    }

    /*!
      Copy vector element of a sparse vector to a row of the internal stored matrix
      \param[in] nnz        number of nonzeros of the columns
      \param[in] values     the values of the sparse vector
      \param[in] col        index position of the values of the sparse vector
      \param[in] irow       the column to be changed
      \param[in] index_base offset for the index, 0 for C based vector 1 for FORTRAN based vector
    \*/
    void
    load_sparse_row( integer         nnz,
                     valueType const values[],
                     integer   const col[],
                     integer         irow,
                     integer         index_base = 0 ) {
      valueType * Arow = Amat + irow ;
      zero( nRow, Arow, nRow ) ;
      for ( integer i = 0 ; i < nnz ; ++i ) Arow[col[i]-index_base] = values[i] ;
    }

    /*!
      Copy a sparse matrix into the internal stored matrix
      \param[in] nnz        number of nonzeros of the columns
      \param[in] values     the values of the sparse vector
      \param[in] row        index row position of the values of the sparse vector
      \param[in] col        index column position of the values of the sparse vector
      \param[in] index_base offset for the index, 0 for C based vector 1 for FORTRAN based vector
    \*/
    void
    load_sparse( integer         nnz,
                 valueType const values[],
                 integer   const row[],
                 integer   const col[],
                 integer         index_base = 0 ) {
      zero( nRow*nCol, Amat, 1 ) ;
      for ( integer i = 0 ; i < nnz ; ++i )
        Amat[row[i]-index_base + (col[i]-index_base) * nRow] = values[i] ;
    }

  } ;

  //============================================================================
  /*\
   |   _    _   _
   |  | |  | | | |
   |  | |  | | | |
   |  | |__| |_| |
   |  |_____\___/
  \*/
  template <typename T>
  class LU : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType ;
  
  private:

    valueType * Work ;
    integer   * Iwork ;
    integer   * i_pivot ;
  
    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;
    
    void check_ls( char const who[] ) const ;

  public:

    LU() ;
    virtual ~LU();

    virtual
    void
    allocate( integer NR, integer NC ) ;

    virtual
    void
    factorize() ;

    virtual
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ;

    virtual
    void
    solve( valueType xb[] ) const ;

    virtual
    void
    t_solve( valueType xb[] ) const ;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const ;

    valueType cond1( valueType norm1 ) const ;
    valueType condInf( valueType normInf ) const ;

  } ;

  //============================================================================
  /*\
   |    ___  ____
   |   / _ \|  _ \
   |  | | | | |_) |
   |  | |_| |  _ <
   |   \__\_\_| \_\
  \*/
  template <typename T>
  class QR : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType ;
  
  private:

    Malloc<valueType> allocReals ;

  protected:

    valueType * Work ;
    valueType * Tau ;

    integer nReflector, Lwork ;

  public:

    QR()
    : Factorization<T>()
    , allocReals("QR-allocReals")
    , nReflector(0)
    , Lwork(0)
    {}

    QR( integer nr, integer nc )
    : Factorization<T>()
    , allocReals("QR-allocReals")
    , nReflector(0)
    , Lwork(0)
    { allocate(nr,nc) ; }

    virtual
    ~QR() {
      allocReals.free() ;
    }

    void
    allocate( integer nr, integer nc, integer Lwrk ) ;

    virtual
    void
    allocate( integer nr, integer nc ) ;

    virtual
    void
    factorize() {
      integer info = geqrf( this->nRow, this->nCol, this->Amat, this->nRow, Tau, Work, Lwork ) ;
      ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::geqrf return info = " << info ) ;
    }

    /*!
      Do QR factorization of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
      allocate( NR, NC ) ;
      integer info = gecopy( NR, NC, A, LDA, this->Amat, this->nRow ) ;
      ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::gecopy return info = " << info ) ;
      factorize() ;
    }

    /*!
      Do QR factorization of the transpose of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    t_factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
      allocate( NC, NR ) ;
      for ( integer i = 0 ; i < NR ; ++i )
        copy( NC, A+i, LDA, this->Amat + i*this->nRow, 1 ) ;
      factorize() ;
    }

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ QR x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    solve( valueType xb[] ) const ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ (QR)^T x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    t_solve( valueType xb[] ) const ;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    virtual
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
    { applyQ( LEFT, NO_TRANSPOSE, nReflector, this->nRow, 1, x, this->nRow ) ; }
    
    //! x <- Q'*x
    void
    Qt_mul( valueType x[] ) const
    { applyQ( LEFT, TRANSPOSE, nReflector, this->nRow, 1, x, this->nRow ) ; }

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
      trsv( UPPER, TRANS, NON_UNIT, rk, this->Amat, this->nRow, x, incx ) ;
    }

    void
    Rsolve( SideMultiply  SIDE,
            Transposition TRANS,
            integer       nr,
            integer       nc,
            valueType     alpha,
            valueType     Bmat[],
            integer       ldB ) const {
      trsm( SIDE, UPPER, TRANS, NON_UNIT, nr, nc, alpha, this->Amat, this->nRow, Bmat, ldB ) ;
    }
    
    //! x <- R^(-1) * x
    void
    invR_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, NO_TRANSPOSE, NON_UNIT, nReflector, this->Amat, this->nRow, x, incx ) ; }

    //! x <- R^(-T) * x
    void
    invRt_mul( valueType x[], integer incx = 1 ) const
    { trsv( UPPER, TRANSPOSE, NON_UNIT, nReflector, this->Amat, this->nRow, x, incx ) ; }

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
    
    // -------------------------------------------------------------------------
    // dummy routines
    void permute( valueType [] ) const {}
    void inv_permute( valueType [] ) const {}
    void permute_rows( integer, integer, valueType [], integer ) const { }
    void inv_permute_rows( integer, integer, valueType [], integer ) const { }
  } ;

  //============================================================================
  /*\
   |    ___  ____  ____
   |   / _ \|  _ \|  _ \
   |  | | | | |_) | |_) |
   |  | |_| |  _ <|  __/
   |   \__\_\_| \_\_|
  \*/
  template <typename T>
  class QRP : public QR<T> {
  public:
    typedef typename QR<T>::valueType valueType ;

  private:
    Malloc<integer> allocIntegers ;
    integer * JPVT ;

  public:

    QRP()
    : QR<T>(), allocIntegers("QRP-allocIntegers")
    {}

    QRP( integer nr, integer nc )
    : QR<T>(), allocIntegers("QRP-allocIntegers")
    { allocate(nr,nc) ; }

    virtual
    ~QRP() {
      allocIntegers.free() ;
    }
    
    virtual
    void
    allocate( integer NR, integer NC ) {
      if ( this->nRow != NR || this->nCol != NC ) {
        valueType tmp ; // get optimal allocation
        integer info = geqp3( NR, NC, nullptr, NR, nullptr, nullptr, &tmp, -1 ) ;
        ALGLIN_ASSERT( info == 0, "QRP::factorize call alglin::geqp3 return info = " << info ) ;
        QR<T>::allocate( NR, NC, integer(tmp) );
      }
      //QR<T>::allocate( NR, NC ) ;
      allocIntegers.allocate(size_t(NC)) ;
      JPVT = allocIntegers(size_t(NC)) ;
    }

    virtual
    void
    factorize() {
      std::fill( JPVT, JPVT + this->nCol, 0 ) ;
      integer info = geqp3( this->nRow, this->nCol,
                            this->Amat, this->nRow,
                            JPVT,
                            this->Tau,
                            this->Work, this->Lwork ) ;
      ALGLIN_ASSERT( info == 0, "QRP::factorize call alglin::geqrf return info = " << info ) ;
    }

    /*!
      Do QR factorization with column pivoting of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
      // calcolo fattorizzazione QR della matrice A
      allocate( NC, NR ) ;
      integer info = gecopy( NR, NC, A, LDA, this->Amat, this->nRow ) ;
      ALGLIN_ASSERT( info == 0, "QR::factorize call alglin::gecopy return info = " << info ) ;
      factorize() ;
    }

    /*!
      Do QR factorization with column pivoting of the transpose of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    t_factorize( integer NR, integer NC, valueType const A[], integer LDA ) {
      // calcolo fattorizzazione QR della matrice A
      allocate( NC, NR ) ;
      for ( integer i = 0 ; i < NR ; ++i )
        copy( NC, A+i, LDA, this->Amat + i*this->nRow, 1 ) ;
      factorize() ;
    }

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ QR x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    solve( valueType xb[] ) const ;

    /*!
      In case of QR factorization of a square matrix solve the 
      linear system \f$ (QR)^T x = b \f$
      \param xb on input the rhs of linear system on output the solution
    \*/
    virtual
    void
    t_solve( valueType xb[] ) const ;

    virtual
    void
    solve( integer nrhs, valueType B[], integer ldB ) const ;

    virtual
    void
    t_solve( integer nrhs, valueType B[], integer ldB ) const ;

    // -------------------------------------------------------------------------
    void
    permute( valueType x[] ) const ;

    void
    inv_permute( valueType x[] ) const ;
    
    void
    permute_rows( integer nr, integer nc, valueType C[], integer ldC ) const {
      ALGLIN_ASSERT( nr == this->nRow,
                     "QRP::permute_rows, bad number of row, expected " << this->nRow <<
                     " find " << nr ) ;
      for ( integer j = 0 ; j < nc ; ++j ) permute( C + ldC*j ) ;
    }
    
    void
    inv_permute_rows( integer nr, integer nc, valueType C[], integer ldC ) const {
      ALGLIN_ASSERT( nr == this->nRow,
                     "QRP::permute_rows, bad number of row, expected " << this->nRow <<
                     " find " << nr ) ;
      for ( integer j = 0 ; j < nc ; ++j ) inv_permute( C + ldC*j ) ;
    }
    
    integer
    rankEstimate( valueType rcond ) const {
      valueType SVAL[3] ;
      return alglin::rankEstimate( this->nRow, this->nCol,
                                   this->Amat, this->nRow,
                                   rcond, SVAL ) ;
    }
  } ;

  //============================================================================
  /*\
   |   ______     ______
   |  / ___\ \   / /  _ \
   |  \___ \\ \ / /| | | |
   |   ___) |\ V / | |_| |
   |  |____/  \_/  |____/
  \*/
  template <typename T>
  class SVD : public Factorization<T> {
  public:
    typedef typename Factorization<T>::valueType valueType ;
  
  protected:

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    valueType * Work ;
    valueType * Umat ;
    valueType * VTmat ;
    valueType * Svec ;
    integer   * IWork ;

    integer     minRC, Lwork ;
    
    typedef enum { USE_GESVD = 0, USE_GESDD = 1 } SVD_USED ;
    SVD_USED    svd_used ;

  public:
  
    using Factorization<T>::solve ;
    using Factorization<T>::t_solve ;

    SVD( SVD_USED _svd_used = USE_GESVD )
    : Factorization<T>()
    , allocReals("SVD-allocReals")
    , allocIntegers("SVD-allocIntegers")
    , Lwork(0)
    , svd_used(_svd_used)
    {}

    virtual
    ~SVD() ALGLIN_OVERRIDE
    { allocReals.free() ; }

    virtual
    void
    allocate( integer NR, integer NC ) ALGLIN_OVERRIDE ;

    virtual
    void
    factorize() ALGLIN_OVERRIDE ;

    /*!
      Do SVD factorization of a rectangular matrix
      \param NR  number of rows of the matrix
      \param NC  number of columns of the matrix
      \param A   pointer to the matrix
      \param LDA Leading dimension of the matrix
    \*/
    virtual
    void
    factorize( integer NR, integer NC, valueType const A[], integer LDA ) ALGLIN_OVERRIDE {
      allocate( NR, NC ) ;
      integer info = gecopy( this->nRow, this->nCol, A, LDA, this->Amat, this->nRow ) ;
      ALGLIN_ASSERT( info == 0, "SVD::factorize call alglin::gecopy return info = " << info ) ;
      factorize() ;
    }

    valueType U( integer i, integer j ) const { return Umat[i+j*this->nRow] ; }
    valueType V( integer i, integer j ) const { return VTmat[j+i*this->nCol] ; }
    valueType sigma( integer i ) const { return Svec[i] ; }

    virtual
    void
    solve( valueType xb[] ) const ALGLIN_OVERRIDE ;

    virtual
    void
    t_solve( valueType xb[] ) const ALGLIN_OVERRIDE ;

    //! y <- alpha * U * x + beta * y
    void
    U_mul( valueType alpha, valueType const x[], integer incx,
           valueType beta,  valueType       y[], integer incy ) const {
      gemv( NO_TRANSPOSE,
            this->nRow, minRC,
            alpha, Umat, this->nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * U' * x + beta * y
    void
    Ut_mul( valueType alpha, valueType const x[], integer incx,
            valueType beta,  valueType       y[], integer incy ) const {
      gemv( TRANSPOSE,
            this->nRow, minRC,
            alpha, Umat, this->nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * V * x + beta * y
    void
    V_mul( valueType alpha, valueType const x[], integer incx,
           valueType beta,  valueType       y[], integer incy ) const {
      gemv( TRANSPOSE,
            minRC, this->nCol,
            alpha, VTmat, this->nRow,
            x, incx,
            beta, y, incy ) ;
    }

    //! y <- alpha * V' * x + beta * y
    void
    Vt_mul( valueType alpha, valueType const x[], integer incx,
            valueType beta,  valueType       y[], integer incy ) const {
      gemv( NO_TRANSPOSE,
            minRC, this->nCol,
            alpha, VTmat, this->nRow,
            x, incx,
            beta, y, incy ) ;
    }

  } ;

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
  class TridiagonalSPD : public LinearSystemSolver<T> {
  public:
    typedef T valueType ;
  
  private:
  
    Malloc<valueType> allocReals ;

    valueType * L ;
    valueType * D ;
    valueType * WORK ;
    integer     nRC ;

  public:

    TridiagonalSPD()
    : allocReals("allocReals")
    , nRC(0)
    {}

    virtual
    ~TridiagonalSPD() {
      allocReals.free() ;
    }

    void
    factorize( integer         N,
               valueType const _L[],
               valueType const _D[] ) ;

    valueType cond1( valueType norm1 ) const ;

    virtual void solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void t_solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;
    virtual void t_solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;

    void
    axpy( integer         N,
          valueType       alpha,
          valueType const L[],
          valueType const D[],
          valueType const x[],
          valueType       beta,
          valueType       y[] ) const ;
  } ;

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
  class TridiagonalLU : public LinearSystemSolver<T> {
  public:
    typedef T valueType ;
  
  private:

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;

    valueType * L ;
    valueType * D ;
    valueType * U ;
    valueType * U2 ;
    valueType * WORK ;
    integer   * IPIV ;
    integer   * IWORK ;

    integer     nRC ;

  public:

    TridiagonalLU()
    : allocReals("allocReals")
    , allocIntegers("allocIntegers")
    , nRC(0)
    {}

    virtual
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

    virtual void solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void t_solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;
    virtual void t_solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;

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
  /*\
   |   _____     _     _ _                               _  ___  ____
   |  |_   _| __(_) __| (_) __ _  __ _  ___  _ __   __ _| |/ _ \|  _ \
   |    | || '__| |/ _` | |/ _` |/ _` |/ _ \| '_ \ / _` | | | | | |_) |
   |    | || |  | | (_| | | (_| | (_| | (_) | | | | (_| | | |_| |  _ <
   |    |_||_|  |_|\__,_|_|\__,_|\__, |\___/|_| |_|\__,_|_|\__\_\_| \_\
   |                             |___/
  \*/
  template <typename T>
  class TridiagonalQR : public LinearSystemSolver<T> {
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

    virtual
    ~TridiagonalQR() {
      allocReals.free() ;
    }

    void
    factorize( integer         N,
               valueType const L[],
               valueType const D[],
               valueType const U[] ) ;

    virtual void solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void t_solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;
    virtual void t_solve( integer nrhs, valueType xb[], integer ldXB ) const ALGLIN_OVERRIDE;

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
    lsq( integer nrhs,
         T       RHS[],
         integer ldRHS,
         T       lambda ) const ;

  } ;

  //============================================================================
  /*\
   |   ____                  _          _ _    _   _
   |  | __ )  __ _ _ __   __| | ___  __| | |  | | | |
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` | |  | | | |
   |  | |_) | (_| | | | | (_| |  __/ (_| | |__| |_| |
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|_____\___/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedLU : public LinearSystemSolver<T> {
  public:
    typedef T valueType ;

    Malloc<valueType> allocReals ;
    Malloc<integer>   allocIntegers ;
    
    integer     m, n, nL, nU, ldAB ;
    integer   * ipiv ;
    valueType * AB ;

    bool is_factorized ;

  public:

    BandedLU() ;
    virtual ~BandedLU();

    void
    setup( integer M,    // number of rows
           integer N,    // number of columns
           integer nL,   // number of lower diagonal
           integer nU ); // number of upper diagonal

    virtual void solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void t_solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void solve( integer nrhs, valueType B[], integer ldB ) const ALGLIN_OVERRIDE;
    virtual void t_solve( integer nrhs, valueType B[], integer ldB ) const ALGLIN_OVERRIDE;

    integer
    iaddr( integer i, integer j ) const {
      integer d = (i-j+nL+nU) ;
      return d+j*ldAB ;
    }

    void
    iaddr_check( integer i, integer j ) const ;

    valueType const &
    operator () ( integer i, integer j ) const { return AB[iaddr(i,j)] ; }

    valueType &
    operator () ( integer i, integer j ) { return AB[iaddr(i,j)] ; }

    void zero() ;

    void
    load_block( integer         nr,
                integer         nc,
                valueType const B[],
                integer         ldB,
                integer         irow,
                integer         icol ) ;

    // do internal factorization, to be executed (only once) before to call solve or t_solve
    void factorize() ;

    // y <- beta*y + alpha*A*x
    void
    aAxpy( valueType       alpha,
           valueType const x[],
           valueType       y[] ) const ;

    void
    dump( std::ostream & stream ) const ;

  } ;

  //============================================================================
  /*\
   |   ____                  _          _ ____  ____  ____
   |  | __ )  __ _ _ __   __| | ___  __| / ___||  _ \|  _ \
   |  |  _ \ / _` | '_ \ / _` |/ _ \/ _` \___ \| |_) | | | |
   |  | |_) | (_| | | | | (_| |  __/ (_| |___) |  __/| |_| |
   |  |____/ \__,_|_| |_|\__,_|\___|\__,_|____/|_|   |____/
  \*/
  //! base class for linear systema solver
  template <typename T>
  class BandedSPD : public LinearSystemSolver<T> {
  public:
    typedef T valueType ;

    Malloc<valueType> allocReals ;

    integer     n, nD, ldAB ;
    valueType * AB ;
    ULselect    UPLO ;
    bool is_factorized ;

  public:

    BandedSPD() ;
    virtual ~BandedSPD();

    void
    setup( ULselect UPLO,
           integer  N,    // numbe of rows and columns
           integer  nD ); // number of upper diagonal

    virtual void solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void t_solve( valueType xb[] ) const ALGLIN_OVERRIDE;
    virtual void solve( integer nrhs, valueType B[], integer ldB ) const ALGLIN_OVERRIDE;
    virtual void t_solve( integer nrhs, valueType B[], integer ldB ) const ALGLIN_OVERRIDE;

    valueType const &
    operator () ( integer i, integer j ) const { return AB[i+j*ldAB] ; }

    valueType &
    operator () ( integer i, integer j ) { return AB[i+j*ldAB] ; }

    void zero() ;

    // do internal fatcorization, to be executed (only once) before to call solve or t_solve
    void factorize() ;

    /* not yet available
    // y <- beta*y + alpha*A*x
    void
    aAxpy( valueType       alpha,
           valueType const x[],
           valueType       y[] ) const ;

    void
    dump( ostream & stream ) const ;
    */

  } ;

  // explicit instantiation declaration to suppress warnings

  #ifdef __GCC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class LU<real> ;
  extern template class LU<doublereal> ;

  extern template class QR<real> ;
  extern template class QR<doublereal> ;

  extern template class QRP<real> ;
  extern template class QRP<doublereal> ;

  extern template class SVD<real> ;
  extern template class SVD<doublereal> ;

  extern template class TridiagonalSPD<real> ;
  extern template class TridiagonalSPD<doublereal> ;

  extern template class TridiagonalLU<real> ;
  extern template class TridiagonalLU<doublereal> ;

  extern template class TridiagonalQR<real> ;
  extern template class TridiagonalQR<doublereal> ;

  extern template class BandedLU<real> ;
  extern template class BandedLU<doublereal> ;

  extern template class BandedSPD<real> ;
  extern template class BandedSPD<doublereal> ;

} // end namespace alglin

#ifdef __GCC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

///
/// eof: Alglin++.hh
///
