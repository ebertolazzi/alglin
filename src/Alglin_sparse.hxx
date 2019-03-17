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
/// file: Alglin_sparse.hxx
///

#include <vector>

namespace alglin {

  template <typename T>
  class SparseMatrixBase {
  public:
    typedef T valueType;
  protected:

    /*!
     * \brief SparseMatrixBase:
     *        Protected Constructor of the class SparseMatrixBase.
    \*/
    SparseMatrixBase() {}

    void
    y_manage(
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const {
      if ( isZero(beta) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] = 0;
      } else if ( !isZero(beta-1) ) {
        for ( integer i = 0; i < DimY; ++i ) y[i*incY] *= beta;
      }
    }

  public:

    /*!
     * \brief ~SparseMatrixBase:
     *        Virtual Destructor of the class SparseMatrixBase.
    \*/
    virtual
    ~SparseMatrixBase() {}

    virtual
    bool
    FORTRAN_indexing() const
    { return false; }

    virtual
    integer
    get_number_of_rows() const ALGLIN_PURE_VIRTUAL;

    virtual
    integer
    get_number_of_cols() const ALGLIN_PURE_VIRTUAL;

    virtual
    integer
    get_nnz() const ALGLIN_PURE_VIRTUAL;

    /*!
     * \brief get_info:
     *        Returns the number of nonzeroes and dimension of the sparse matrix.
     * \param[out] numRows Row dimension of the matrix
     * \param[out] numCols Column dimension of the matrix
     * \param[out] nnz     the number of nonzeroes
     *
    \*/
    virtual
    void
    get_info(
      integer & numRows,
      integer & numCols,
      integer & nnz
    ) const {
      numRows = this->get_number_of_rows();
      numCols = this->get_number_of_cols();
      nnz     = this->get_nnz();
    }

    /*!
     * \brief get_data:
     *        Returns the sparse matrix.
     *
     * Example of usage
     *
     * \code
     * int_type nnz;
     * int_type const * rows;
     * int_type const * cols;
     * real     const * values;
     * ptr->get_data( nnz, &rows, &cols, &values );
     * \endcode
     *
     * \param[out] pRows   vector of pointers of row indices where to data will be copied.
     * \param[out] pCols   vector of pointers where to store column indices.
     * \param[out] pValues vector of pointers where to store values.
     *
    \*/

    virtual
    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      valueType const * & pValues
    ) const ALGLIN_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha A x + y\f$).
     *
     * \param[in]     alpha         Scalar \f$\alpha \f$.
     * \param[in]     DimX          Dimension of the vector x.
     * \param[in]     x             Vector x.
     * \param[in]     incX          Stride of vector x.
     * \param[in]     DimY          Dimension of the vector y.
     * \param[in]     beta          Scalar \f$\beta \f$.
     * \param[in,out] y             In-/Output vector y.
     * \param[in]     incY          Stride of vector y.
     *
    \*/
    virtual
    void
    gemv(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha A^T x + y\f$).
     *
     * \param[in]     alpha         Scalar \f$\alpha \f$.
     * \param[in]     DimX          Dimension of the vector x.
     * \param[in]     x             Vector x.
     * \param[in]     incX          Stride of vector x.
     * \param[in]     DimY          Dimension of the vector y.
     * \param[in]     beta          Scalar \f$\beta \f$.
     * \param[in,out] y             In-/Output vector y.
     * \param[in]     incY          Stride of vector y.
     *
    \*/

    virtual
    void
    gemv_Transposed(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_PURE_VIRTUAL;

    /*!
     * \brief gemv:
     *        Calls the blas-Routine (\f$y = \beta y + \alpha (A+A^T-\textrm{diag}(A)) x + y\f$).
     *
     * \param[in]     alpha         Scalar \f$\alpha \f$.
     * \param[in]     DimX          Dimension of the vector x.
     * \param[in]     x             Vector x.
     * \param[in]     incX          Stride of vector x.
     * \param[in]     DimY          Dimension of the vector y.
     * \param[in]     beta          Scalar \f$\beta \f$.
     * \param[in,out] y             In-/Output vector y.
     * \param[in]     incY          Stride of vector y.
     *
    \*/
    virtual
    void
    gemv_Symmetric(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_PURE_VIRTUAL;

  };

  /*
  //   ____                             ____ ____ ___   ___  ____
  //  / ___| _ __   __ _ _ __ ___  ___ / ___/ ___/ _ \ / _ \|  _ \
  //  \___ \| '_ \ / _` | '__/ __|/ _ \ |  | |  | | | | | | | |_) |
  //   ___) | |_) | (_| | |  \__ \  __/ |__| |__| |_| | |_| |  _ <
  //  |____/| .__/ \__,_|_|  |___/\___|\____\____\___/ \___/|_| \_\
  //        |_|
  */

  template <typename T> class MatrixWrapper;

  //! Sparse Matrix Structure
  template <typename T>
  class SparseCCOOR : public SparseMatrixBase<T> {
  public:
    typedef SparseMatrixBase<T>        Sparse;
    typedef MatrixWrapper<T>           MatW;
    typedef typename Sparse::valueType valueType;

  protected:
    integer                nRows; //!< Number of rows
    integer                nCols; //!< Number of columns
    integer                nnz;   //!< Total number of nonzeros
    std::vector<valueType> vals;  //!< the values of the sparse matrix
    std::vector<integer>   rows;  //!< the rows index
    std::vector<integer>   cols;  //!< the columns index

    bool fortran_indexing;

  public:

    SparseCCOOR()
    : nRows(0)
    , nCols(0)
    , nnz(0)
    , fortran_indexing(false)
    {}

    SparseCCOOR(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi
    );

    virtual
    ~SparseCCOOR() ALGLIN_OVERRIDE
    {}

    virtual
    bool
    FORTRAN_indexing() const ALGLIN_OVERRIDE
    { return this->fortran_indexing; }

    virtual
    integer
    get_number_of_rows() const ALGLIN_OVERRIDE
    { return this->nRows; }

    virtual
    integer
    get_number_of_cols() const ALGLIN_OVERRIDE
    { return this->nCols; }

    virtual
    integer
    get_nnz() const ALGLIN_OVERRIDE
    { return this->nnz; }

    virtual
    void
    get_data(
      integer   const * & pRows,
      integer   const * & pCols,
      valueType const * & pValues
    ) const ALGLIN_OVERRIDE {
      pRows   = &this->rows.front();
      pCols   = &this->cols.front();
      pValues = &this->vals.front();
    }

    template <typename i_type, typename r_type>
    void
    export_data(
      integer NNZ,
      i_type  i[],
      i_type  j[],
      r_type  val[],
      bool    fi
    ) {
      ALGLIN_ASSERT(
        NNZ == this->rows.size() &&
        NNZ == this->cols.size() &&
        NNZ == this->vals.size(),
        "export_data, bad dimension"
      );
      integer offs = 0;
      if ( fi ) ++offs;
      if ( this->fortran_indexing ) --offs;
      for ( integer index = 0; index < NNZ; ++index ) {
        i[index]   = i_type(this->rows[index]+offs);
        j[index]   = i_type(this->cols[index]+offs);
        val[index] = r_type(this->vals[index]);
      }
    }

    template <typename i_type, typename r_type>
    void
    export_data(
      std::vector<i_type> & i,
      std::vector<i_type> & j,
      std::vector<r_type> & v,
      bool                  fi
    ) {
      integer offs = 0;
      if ( fi ) ++offs;
      if ( this->fortran_indexing ) --offs;
      i.clear(); i.reserve( this->nnz );
      j.clear(); j.reserve( this->nnz );
      v.clear(); v.reserve( this->nnz );
      for ( integer index = 0; index < this->nnz; ++index ) {
        i.push_back( i_type(this->rows[index]+offs) );
        j.push_back( i_type(this->cols[index]+offs) );
        v.push_back( r_type(this->vals[index]) );
      }
    }

    virtual
    void
    gemv(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_OVERRIDE;

    virtual
    void
    gemv_Transposed(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_OVERRIDE;

    virtual
    void
    gemv_Symmetric(
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      valueType       beta,
      integer         DimY,
      valueType       y[],
      integer         incY
    ) const ALGLIN_OVERRIDE;

    /*\
    :|:           _    _ _ _   _               _             _   _            _
    :|:   __ _ __| |__| (_) |_(_)___ _ _  __ _| |  _ __  ___| |_| |_  ___  __| |___
    :|:  / _` / _` / _` | |  _| / _ \ ' \/ _` | | | '  \/ -_)  _| ' \/ _ \/ _` (_-<
    :|:  \__,_\__,_\__,_|_|\__|_\___/_||_\__,_|_| |_|_|_\___|\__|_||_\___/\__,_/__/
    \*/

    void clear();

    void
    init(
      integer N,
      integer M,
      integer reserve_nnz,
      bool    fi = false
    );

    void
    transpose() {
      this->rows.swap(this->cols);
      std::swap( this->nRows, this->nCols);
    }

    void to_FORTRAN_indexing();
    void to_C_indexing();

    void
    setup_as_full_row_major( integer N, integer M, bool fi = false );

    void
    setup_as_full_column_major( integer N, integer M, bool fi = false );

    void
    fill( valueType const V[], integer M );

    void
    fill( std::vector<valueType> const & V );

    void
    reserve( integer reserve_nnz );

    void
    push_value_C(
      integer   row,
      integer   col,
      valueType val
    );

    void
    push_value_F(
      integer   row,
      integer   col,
      valueType val
    );

    void
    push_matrix(
      integer        row_offs,
      integer        col_offs,
      Sparse const & Matrix,
      bool           transpose = false
    );

    void
    get_matrix( MatW & M ) const;

    void
    get_matrix_symmetric( MatW & M ) const;

    void
    get_matrix_transposed( MatW & M ) const;

    bool foundNaN() const;
  };

  // explicit instantiation declaration to suppress warnings

  #ifdef ALGLIN_USE_CXX11
  extern template class SparseCCOOR<real>;
  extern template class SparseCCOOR<doublereal>;
  #endif

}
