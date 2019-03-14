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
     *        Calls the blas-Routine (\f$y = \beta y + \alpha \mathrm{op}(A) x + y\f$).
     *
     * \param[in]     TransposedA   True if the matrix A is transposed.
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
      bool            TransposedA,
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      integer         DimY,
      valueType       beta,
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

  //! Sparse Matrix Structure
  template <typename T>
  class SparseCCOOR : public SparseMatrixBase<T> {
  public:
    typedef typename SparseMatrixBase<T>::valueType valueType;
  private:
    integer     nRows; //!< Number of rows
    integer     nCols; //!< Number of columns
    integer     nnz;   //!< Total number of nonzeros
    valueType * vals;  //!< pointer to the values of the sparse matrix
    integer   * rows;  //!< pointer to the rows index
    integer   * cols;  //!< pointer to the columns index

  public:

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
      pRows   = this->rows;
      pCols   = this->cols;
      pValues = this->vals;
    }

    virtual
    void
    gemv(
      bool            TransposedA,
      valueType       alpha,
      integer         DimX,
      valueType const x[],
      integer         incX,
      integer         DimY,
      valueType       beta,
      valueType       y[],
      integer         incY
    ) const ALGLIN_OVERRIDE;

    bool foundNaN() const;
  };

  // explicit instantiation declaration to suppress warnings

  #ifdef ALGLIN_USE_CXX11
  extern template class SparseCCOOR<real>;
  extern template class SparseCCOOR<doublereal>;
  #endif

}


