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
    valueType * vals;    //!< pointer to the values of the sparse matrix
    integer   * rows;    //!< pointer to the rows index
    integer   * cols;    //!< pointer to the columns index
    integer     numRows; //!< Number of rows
    integer     numCols; //!< Number of columns
    integer     nnz;     //!< Total number of nonzeros
  };

  /*
  //   __  __       _        _     __        __
  //  |  \/  | __ _| |_ _ __(_)_  _\ \      / / __ __ _ _ __  _ __   ___ _ __
  //  | |\/| |/ _` | __| '__| \ \/ /\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
  //  | |  | | (_| | |_| |  | |>  <  \ V  V /| | | (_| | |_) | |_) |  __/ |
  //  |_|  |_|\__,_|\__|_|  |_/_/\_\  \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
  //                                                 |_|   |_|
  */

  //! Sparse Matrix Structure
  template <typename T>
  class MatrixWrapper {
    typedef T                valueType;
    typedef MatrixWrapper<T> MatW ;
    typedef SparseCCOOR<T>   Sparse ;

    valueType *data;
    integer   numRows; //!< Number of rows
    integer   numCols; //!< Number of columns
    integer   ldData;  //!< Leadind dimension

    integer
    iaddr( integer i,  integer j ) const
    { return i + j*ldData; }

  public:

    explicit
    MatrixWrapper( valueType * _data,
                   integer     nr,
                   integer     nc,
                   integer     ld )
    : data(_data)
    , numRows(nr)
    , numCols(nc)
    , ldData(ld)
    {
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
    load( Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] = sp.data[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] += sp.data[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( valueType alpha, Sparse const & sp ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx],sp.cols[idx])] += alpha * sp.data[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    load( Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] = sp.data[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += sp.data[idx] ;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void
    add( valueType alpha, Sparse const & sp, integer i_offs, integer j_offs ) {
      for ( integer idx = 0 ; idx < sp.nnz ; ++idx )
        data[iaddr(sp.rows[idx]+i_offs,sp.cols[idx]+j_offs)] += alpha * sp.data[idx] ;
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

  };

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


