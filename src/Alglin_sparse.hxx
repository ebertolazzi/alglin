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

///
/// file: alglin_sparse.hxx
///

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
  class SparseCCOOR {
  public:
    typedef T valueType;
    integer     numRows; //!< Number of rows
    integer     numCols; //!< Number of columns
    integer     nnz;     //!< Total number of nonzeros
    integer     __dummy;
    valueType * vals;    //!< pointer to the values of the sparse matrix
    integer   * rows;    //!< pointer to the rows index
    integer   * cols;    //!< pointer to the columns index

    bool
    foundNaN() const
    { return alglin::foundNaN( vals, nnz ); }

  };

  /*
  //   __  __       _        _     __        __
  //  |  \/  | __ _| |_ _ __(_)_  _\ \      / / __ __ _ _ __  _ __   ___ _ __
  //  | |\/| |/ _` | __| '__| \ \/ /\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
  //  | |  | | (_| | |_| |  | |>  <  \ V  V /| | | (_| | |_) | |_) |  __/ |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\  \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
  //                                                   |_|   |_|
  */

  //! Sparse Matrix Structure
  template <typename T>
  class MatrixWrapper {
    typedef T                valueType;
    typedef MatrixWrapper<T> MatW ;
    typedef SparseCCOOR<T>   Sparse ;

    integer
    iaddr( integer i,  integer j ) const
    { return i + j*ldData; }

  public:

    // public access to data, use with caution
    integer   numRows; //!< Number of rows
    integer   numCols; //!< Number of columns
    integer   ldData;  //!< Leadind dimension
    integer   __dummy; // padding
    valueType *data;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    explicit
    MatrixWrapper( )
    : numRows(0)
    , numCols(0)
    , ldData(0)
    , data(nullptr)
    {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    explicit
    MatrixWrapper( valueType * _data,
                   integer     nr,
                   integer     nc,
                   integer     ld )
    : numRows(nr)
    , numCols(nc)
    , ldData(ld)
    , data(_data)
    {
    #if defined(_DEBUG) || defined(DEBUG)
      ALGLIN_ASSERT( nr >= 0 && nc >= 0 && ldData >= nr,
                     "MatrixWrapper( data, nr=" << nr <<
                     ", nc=" << nc << ", ld=" << ld <<
                     ") bad dimensions" );
    #endif
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    setup( valueType * _data,
           integer     nr,
           integer     nc,
           integer     ld ) {
      data    = _data;
      numRows = nr;
      numCols = nc;
      ldData  = ld;
      #if defined(_DEBUG) || defined(DEBUG)
      ALGLIN_ASSERT( nr >= 0 && nc >= 0 && ldData >= nr,
                     "MatrixWrapper( data, nr=" << nr <<
                     ", nc=" << nc << ", ld=" << ld <<
                     ") bad dimensions" );
      #endif
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    valueType const &
    operator () ( integer i,  integer j ) const
    { return data[iaddr(i,j)]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    valueType &
    operator () ( integer i,  integer j )
    { return data[iaddr(i,j)]; }

    //! fill matrix with 0
    void zero() { gezero( numRows, numCols, data, ldData ); }
    void id( valueType dg ) { geid( numRows, numCols, data, ldData, dg ); }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    load( MatW const & A ) {
      integer info = gecopy( A.numRows, A.numCols, A.data, A.ldData,
                             data, ldData );
      ALGLIN_ASSERT( info == 0,
                     "MatrixWrapper::load call alglin::gecopy return info = " << info );

    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    load0( Sparse const & sp ) {
      zero();
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] = sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    load( Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] = sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] += sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( valueType alpha, Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] += alpha * sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    load( Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] = sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += sp.vals[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( valueType alpha, Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += alpha * sp.vals[idx] ;
    }

    // alpha*A + beta*B -> C
    friend
    void
    add( valueType    alpha,
         MatW const & A,
         valueType    beta,
         MatW const & B,
         MatW       & C ) {
      geadd( C.numRows, C.numCols,
             alpha, A.data, A.ldData,
             beta,  B.data, B.ldData,
             C.data, C.ldData );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    block( integer i_offs,
           integer j_offs,
           integer nrow,
           integer ncol,
           MatW & to ) {
      to.setup( data+iaddr(i_offs,j_offs), nrow, ncol, ldData ) ;
    }

  };

  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm( char const               where[],
        T                        alpha,
        MatrixWrapper<T> const & A,
        MatrixWrapper<T> const & B,
        T                        beta,
        MatrixWrapper<T>       & C ) {
    ALGLIN_ASSERT( A.numCols == B.numRows &&
                   A.numRows == C.numRows &&
                   B.numCols == C.numCols,
                   "gemm, inconsistent dimensions: " << where <<
                   "\nA = " << A.numRows << " x " << A.numCols <<
                   "\nB = " << B.numRows << " x " << B.numCols <<
                   "\nC = " << C.numRows << " x " << C.numCols );
    alglin::gemm( NO_TRANSPOSE,
                  NO_TRANSPOSE,
                  A.numRows, B.numCols, A.numCols,
                  alpha,
                  A.data, A.ldData,
                  B.data, B.ldData,
                  beta,
                  C.data, C.ldData );
  }

  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm( char const               where[],
        T                        alpha,
        Transposition    const & TRANSA,
        MatrixWrapper<T> const & A,
        Transposition    const & TRANSB,
        MatrixWrapper<T> const & B,
        T                        beta,
        MatrixWrapper<T>       & C ) {

    integer Ar = TRANSA == NO_TRANSPOSE ? A.numRows : A.numCols;
    integer Ac = TRANSA == NO_TRANSPOSE ? A.numCols : A.numRows;
    integer Br = TRANSB == NO_TRANSPOSE ? B.numRows : B.numCols;
    integer Bc = TRANSB == NO_TRANSPOSE ? B.numCols : B.numRows;

    ALGLIN_ASSERT( C.numRows == Ar && C.numCols == Bc && Ac == Br,
                   "gemm, inconsistent dimensions: " << where <<
                   "\nA = " << A.numRows << " x " << A.numCols <<
                   "\nB = " << B.numRows << " x " << B.numCols <<
                   "\nC = " << C.numRows << " x " << C.numCols <<
                   "\nA " << (NO_TRANSPOSE?"NO":"") << " transposed" <<
                   "\nB " << (NO_TRANSPOSE?"NO":"") << " transposed" );
    alglin::gemm( TRANSA,
                  TRANSB,
                  C.numRows, C.numCols, Ac,
                  alpha,
                  A.data, A.ldData,
                  B.data, B.ldData,
                  beta,
                  C.data, C.ldData );
  }

  // explicit instantiation declaration to suppress warnings

  #ifdef ALGLIN_USE_CXX11

  #ifdef __GCC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
  #endif
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class SparseCCOOR<real>;
  extern template class SparseCCOOR<doublereal>;

  extern template class MatrixWrapper<real>;
  extern template class MatrixWrapper<doublereal>;

  #endif

}


