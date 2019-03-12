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
/// file: Alglin_wrapper.hxx
///

namespace alglin {

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
    typedef MatrixWrapper<T> MatW;
    typedef SparseCCOOR<T>   Sparse;

    #if defined(DEBUG) || defined(_DEBUG)
    integer
    iaddr( integer i,  integer j ) const {
      ALGLIN_ASSERT(
        i >= 0 && i < numRows && j >= 0 && j < numCols,
        "iaddr(" << i << ", " << j << ") out of range [0," <<
        numRows << ") x [0," << numCols << ")"
      );
      return i + j*ldData;
    }
    #else
    integer
    iaddr( integer i,  integer j ) const
    { return i + j*ldData; }
    #endif

    void
    check( Sparse const & sp ) const;

    void
    check( MatW const & M ) const;

  public:

    // public access to data, use with caution
    integer     numRows; //!< Number of rows
    integer     numCols; //!< Number of columns
    integer     ldData;  //!< Leadind dimension
    valueType * data;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: build an empy wrapper
    \*/
    explicit
    MatrixWrapper( )
    : numRows(0)
    , numCols(0)
    , ldData(0)
    , data(nullptr)
    {
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data  pointer of the memory to be mapped as a matrix
    :|: \param nr     number of rows of the mapped matrix
    :|: \param nc     number of columns of the mapped matrix
    :|: \param ld     leading dimension of the matrix (Fortran addressing 0 based)
    \*/
    explicit
    MatrixWrapper(
      valueType * _data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Map a piece of memory into a matrix object
    :|:
    :|: \param _data  pointer of the memory to be mapped as a matrix
    :|: \param nr     number of rows of the mapped matrix
    :|: \param nc     number of columns of the mapped matrix
    :|: \param ld     leading dimension of the matrix (Fortran addressing 0 based)
    \*/
    void
    setup(
      valueType * _data,
      integer     nr,
      integer     nc,
      integer     ld
    );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Access element (i,j) of the matrix
    :|:
    :|: \param i row of the element
    :|: \param j column of the element
    \*/
    valueType const &
    operator () ( integer i,  integer j ) const
    { return data[iaddr(i,j)]; }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Access element (i,j) of the matrix
    :|:
    :|: \param i row of the element
    :|: \param j column of the element
    \*/
    valueType &
    operator () ( integer i,  integer j )
    { return data[iaddr(i,j)]; }

    /*!
    :|: Fill the matrix with zeros
    \*/
    void zero() { gezero( numRows, numCols, data, ldData ); }

    /*!
    :|: Initialize the matrix as an identity matrix
    \*/
    void id( valueType dg ) { geid( numRows, numCols, data, ldData, dg ); }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Initialize matrix
    :|:
    :|: \param data   pointer of memory with data to be copied
    :|: \param ldData leading dimension of the memory to be copied
    :|:
    \*/
    void load( valueType const data[], integer ldData );

    /*!
    :|: Initialize matrix
    :|:
    :|: \param A initialize matrix with the matrix of the object A
    :|:
    \*/
    void load( MatW const & A );

    /*!
    :|: Initialize matrix to 0 and next copy sparse matrix
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void load0( Sparse const & sp );

    /*!
    :|: Copy sparse matrix into the object
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void load( Sparse const & sp );

    /*!
    :|: Copy sparse matrix into the object
    :|:
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void load( Sparse const & sp, integer i_offs, integer j_offs );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Add a block of memory to the matrix
    :|:
    :|: \param data   pointer of memory with data to be copied
    :|: \param ldData leading dimension of the memory to be copied
    :|:
    \*/
    void add( valueType const data[], integer ldData );

    /*!
    :|: Add sparse matrix to the object
    :|:
    :|: \param sp sparse matrix to be copied
    :|:
    \*/
    void add( Sparse const & sp );

    /*!
    :|: Add sparse matrix to the object
    :|:
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void add( Sparse const & sp, integer i_offs, integer j_offs );

    /*!
    :|: Add sparse multiplied by `alpha` matrix to the object
    :|:
    :|: \param alpha  scalar used to multiply the sparse matrix
    :|: \param sp     sparse matrix to be copied
    :|:
    \*/
    void add( valueType alpha, Sparse const & sp );

    /*!
    :|: Add sparse multiplied by `alpha` matrix to the object
    :|:
    :|: \param alpha  scalar used to multiply the sparse matrix
    :|: \param sp     sparse matrix to be copied
    :|: \param i_offs offset added to to the row indices
    :|: \param j_offs offset added to to the column indices
    :|:
    \*/
    void add( valueType alpha, Sparse const & sp, integer i_offs, integer j_offs );

    // alpha*A + beta*B -> C
    /*!
    :|: Add a linear combination of two matrix and assign to a third one
    :|:
    :|: \param alpha scalar used to multiply the matrix `A`
    :|: \param A     first mapped matrix
    :|: \param beta  scalar used to multiply the matrix `B`
    :|: \param B     second mapped matrix
    :|: \param C     mapped matrix which store the result
    :|:
    \*/
    friend
    void
    add(
      valueType    alpha,
      MatW const & A,
      valueType    beta,
      MatW const & B,
      MatW       & C
    ) {
      geadd( C.numRows, C.numCols,
             alpha, A.data, A.ldData,
             beta,  B.data, B.ldData,
             C.data, C.ldData );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    /*!
    :|: Extract a block
    :|:
    :|: \param i_offs initial row of the block to be extracted
    :|: \param j_offs initial column of the block to be extracted
    :|: \param nrow   number of rows of the block to be extracted
    :|: \param ncol   number of columns of the block to be extracted
    :|: \param to     mapped matrix which store the result
    :|:
    \*/
    void
    block(
      integer i_offs,
      integer j_offs,
      integer nrow,
      integer ncol,
      MatW &  to
    ) {
      #if defined(DEBUG) || defined(_DEBUG)
      ALGLIN_ASSERT(
        i_offs >= 0 && i_offs+nrow <= numRows &&
        j_offs >= 0 && j_offs+ncol <= numCols,
        "block(i_offs=" << i_offs << ", j_offs=" << j_offs <<
        ", nrow=" << nrow << ", ncol=" << ncol <<
        ") will be out of range [0," <<
        numRows << ") x [0," << numCols << ")"
      );
      #endif
      to.setup( data+iaddr(i_offs,j_offs), nrow, ncol, ldData );
    }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  /*!
  :|: Perform matrix vector multiplication `c = beta*c + alpha*A*v`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose if multiply with the transpose of `A`
  :|: \param A      matrix used in the multiplication
  :|: \param v      vector to be multiplied
  :|: \param incv   stride of the vector `v`
  :|: \param beta   scalar used to multiply `c`
  :|: \param c      result vector
  :|: \param incc   stride of the vector `c`
  :|:
  \*/
  template <typename T>
  inline
  void
  gemv(
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    T const                  v[],
    integer                  incv,
    T                        beta,
    T                        c[],
    integer                  incc
  ) {
    alglin::gemv( TRANSA,
                  A.numRows, A.numCols,
                  alpha,
                  A.data, A.ldData,
                  v, incv,
                  beta,
                  c, incc );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // c = beta*c + alpha*A*v
  /*!
  :|: Perform matrix vector multiplication `c = beta*c + alpha*A*v`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param v      vector to be multiplied
  :|: \param incv   stride of the vector `v`
  :|: \param beta   scalar used to multiply `c`
  :|: \param c      result vector
  :|: \param incc   stride of the vector `c`
  :|:
  \*/
  template <typename T>
  inline
  void
  gemv(
    T                        alpha,
    MatrixWrapper<T> const & A,
    T const                  v[],
    integer                  incv,
    T                        beta,
    T                        c[],
    integer                  incc
  ) {
    alglin::gemv( NO_TRANSPOSE,
                  A.numRows, A.numCols,
                  alpha,
                  A.data, A.ldData,
                  v, incv,
                  beta,
                  c, incc );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param where  message added in case of error
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm(
    char const               where[],
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    ALGLIN_ASSERT(
      A.numCols == B.numRows &&
      A.numRows == C.numRows &&
      B.numCols == C.numCols,
      "gemm, at `" << where << "' inconsistent dimensions: " <<
      "\nA = " << A.numRows << " x " << A.numCols <<
      "\nB = " << B.numRows << " x " << B.numCols <<
      "\nC = " << C.numRows << " x " << C.numCols
    );
    alglin::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows, B.numCols, A.numCols,
      alpha,
      A.data, A.ldData,
      B.data, B.ldData,
      beta,
      C.data, C.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param A      matrix used in the multiplication
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
 // C = beta*C + alpha*A*B (NO DEBUG)
  template <typename T>
  inline
  void
  gemm(
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    alglin::gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      A.numRows, B.numCols, A.numCols,
      alpha,
      A.data, A.ldData,
      B.data, B.ldData,
      beta,
      C.data, C.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param where  message added in case of error
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose `A` transposed or not
  :|: \param A      matrix used in the multiplication
  :|: \param TRANSB choose `B` transposed or not
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  // C = beta*C + alpha*A*B
  template <typename T>
  inline
  void
  gemm(
    char const               where[],
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    Transposition    const & TRANSB,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {

    integer Ar = TRANSA == NO_TRANSPOSE ? A.numRows : A.numCols;
    integer Ac = TRANSA == NO_TRANSPOSE ? A.numCols : A.numRows;
    integer Br = TRANSB == NO_TRANSPOSE ? B.numRows : B.numCols;
    integer Bc = TRANSB == NO_TRANSPOSE ? B.numCols : B.numRows;

    ALGLIN_ASSERT(
      C.numRows == Ar && C.numCols == Bc && Ac == Br,
      "gemm, at `" << where << "' inconsistent dimensions: " <<
      "\nA = " << A.numRows << " x " << A.numCols <<
      "\nB = " << B.numRows << " x " << B.numCols <<
      "\nC = " << C.numRows << " x " << C.numCols <<
      "\nA " << (NO_TRANSPOSE?"NO":"") << " transposed" <<
      "\nB " << (NO_TRANSPOSE?"NO":"") << " transposed"
    );
    alglin::gemm(
      TRANSA,
      TRANSB,
      C.numRows, C.numCols, Ac,
      alpha,
      A.data, A.ldData,
      B.data, B.ldData,
      beta,
      C.data, C.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // C = beta*C + alpha*A*B (NO DEBUG)
  /*!
  :|: Perform matrix matrix multiplication `C = beta*C + alpha*A*B`
  :|:
  :|: \param alpha  matrix `A` is multiplied by `alpha`
  :|: \param TRANSA choose `A` transposed or not
  :|: \param A      matrix used in the multiplication
  :|: \param TRANSB choose `B` transposed or not
  :|: \param B      matrix used in the multiplication
  :|: \param beta   scalar used to multiply `C`
  :|: \param C      result matrix
  :|:
  \*/
  template <typename T>
  inline
  void
  gemm(
    T                        alpha,
    Transposition    const & TRANSA,
    MatrixWrapper<T> const & A,
    Transposition    const & TRANSB,
    MatrixWrapper<T> const & B,
    T                        beta,
    MatrixWrapper<T>       & C
  ) {
    integer Ac = TRANSA == NO_TRANSPOSE ? A.numCols : A.numRows;
    alglin::gemm(
      TRANSA,
      TRANSB,
      C.numRows, C.numCols, Ac,
      alpha,
      A.data, A.ldData,
      B.data, B.ldData,
      beta,
      C.data, C.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trmv(
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trmv, matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trmv(
      UPLO, TRANS, DIAG, A.numRows, A.data, A.ldData, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A * x or x = A^T * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param where  message added in case of error
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trmv(
    char const               where[],
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trmv at `" << where << "` matrix is " <<
      A.numRows << " x " << A.numRows << " expected square"
    );
    alglin::trmv(
      UPLO, TRANS, DIAG, A.numRows, A.data, A.ldData, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
  :|: Perform matrix vector multiplication `x = A * x` or `x = A^T * x`
  :|:
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trsv(
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trsv, matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trsv(
      UPLO, TRANS, DIAG, A.numRows, A.data, A.ldData, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // x = A^(-1) * x or x = A^(-T) * x
  /*!
  :|: Perform matrix vector multiplication `x = A^(-1) * x` or `x = A^(-T) * x`
  :|:
  :|: \param where  message added in case of error
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param A      matrix to multiplied
  :|: \param x      vector used in the multiplcation and result
  :|: \param incx   stride of vector `x`
  :|:
  \*/
  template <typename T>
  inline
  void
  trsv(
    char const               where[],
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    MatrixWrapper<T> const & A,
    T                        x[],
    integer                  incx
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trsv at `" << where << "` matrix is " <<
      A.numRows << " x " << A.numRows << " expected square"
    );
    alglin::trsv(
      UPLO, TRANS, DIAG, A.numRows, A.data, A.ldData, x, incx
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trmm(
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trmm, matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows, B.numCols, alpha,
      A.data, A.ldData,
      B.data, B.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param where  message added in case of error
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trmm(
    char const               where[],
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trmm, at `" << where <<
      "` matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trmm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows, B.numCols, alpha,
      A.data, A.ldData,
      B.data, B.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trsm(
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trsm, matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows, B.numCols, alpha,
      A.data, A.ldData,
      B.data, B.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /*!
  :|: Perform matrix matrix multiplication
  :|:
  :|:    B := alpha*op( A^(-1) )*B,   or   B := alpha*B*op( A^(-1) ),
  :|:    op( A ) = A   or   op( A ) = A**T.
  :|:
  :|: \param where  message added in case of error
  :|: \param SIDE   multiply on left or right
  :|: \param UPLO   choose upper or lower part
  :|: \param TRANS  choose `A` transposed or not
  :|: \param DIAG   use or not diagonal elements to 1
  :|: \param alpha  scalar used in the multiplication
  :|: \param A      matrix to be multiplied
  :|: \param B      matrix to be multiplied
  :|:
  \*/
  template <typename T>
  inline
  void
  trsm(
    char const               where[],
    SideMultiply     const & SIDE,
    ULselect         const & UPLO,
    Transposition    const & TRANS,
    DiagonalType     const & DIAG,
    T                        alpha,
    MatrixWrapper<T> const & A,
    MatrixWrapper<T>       & B
  ) {
    ALGLIN_ASSERT(
      A.numRows == A.numCols,
      "trmm, at `" << where <<
      "` matrix is " << A.numRows << " x " <<
      A.numRows << " expected square"
    );
    alglin::trsm(
      SIDE, UPLO, TRANS, DIAG,
      B.numRows, B.numCols, alpha,
      A.data, A.ldData,
      B.data, B.ldData
    );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normInf( MatrixWrapper<T> const & A ) {
    return alglin::normInf( A.numRows, A.numCols, A.data, A.ldData );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  norm1( MatrixWrapper<T> const & A ) {
    return alglin::norm1( A.numRows, A.numCols, A.data, A.ldData );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  normF( MatrixWrapper<T> const & A ) {
    return alglin::normF( A.numRows, A.numCols, A.data, A.ldData );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename T>
  inline
  T
  maxabs( MatrixWrapper<T> const & A ) {
    return alglin::maxabs( A.numRows, A.numCols, A.data, A.ldData );
  }

  // explicit instantiation declaration to suppress warnings

  #ifdef ALGLIN_USE_CXX11
  extern template class MatrixWrapper<real>;
  extern template class MatrixWrapper<doublereal>;
  #endif

}
