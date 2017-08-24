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

#ifndef KKT_LIKE_N_HH
#define KKT_LIKE_N_HH

#include "Alglin.hh"
#include "Alglin++.hh"
#include <vector>

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#endif

namespace alglin {

  using namespace std ;

  /*
  //   _  ___  _______
  //  | |/ / |/ /_   _|
  //  | ' /| ' /  | |
  //  | . \| . \  | |
  //  |_|\_\_|\_\ |_|
  */
  //! LU decomposition of a KKT like matrix
  /*!
   *
   *  Solve the linear system of the form
   *  -----------------------------------
   *
   *  / A B \ / x \   / a \
   *  |     | |   | = |   |
   *  \ C D / \ y /   \ b /
   *
   *  A*x + B*y = a
   *  C*x + D*y = b
   *
   *  x = A^(-1) * ( a - B*y )
   *  C*A^(-1)*a -b = C*A^(-1)*B*y -D*y
   *
   *  Solution procedure
   *
   *  Compute aux matrix
   *  Z = A^(-1)*B
   *  W = C*Z - D
   *
   *  solve the system
   *  a' = A^(-1)*a
   *  b' = C*a' - b
   *
   *  y = W^(-1) * b'
   *  x = a' - Z*y
   *
   *  Solve the linear system of the form
   *  -----------------------------------
   *
   *  / A^T C^T \ / x \   / a \
   *  |         | |   | = |   |
   *  \ B^T D^T / \ y /   \ b /
   *
   *  A^T*x + C^T*y = a
   *  B^T*x + D^T*y = b
   *
   *  x = A^(-T) * ( a - C^T*y )
   *  B^T*A^(-T)*a -b = B^T*A^(-T)*C^T*y - D^T*y
   *                  = (A^(-1)B)^T*C^T*y - D^T*y
   *                  = (C*(A^(-1)B))^T*y - D^T*y
   *                  = W^T*y
   *  (A^(-1)*B)^T*a -b = W^T*y
   *           Z^T*a -b = W^T*y
   *
   *  Solution procedure
   *
   *  b' = Z^T*a -b
   *  y  = W^(-T)*b'
   *  a' = a - C^T*y
   *  x = A^(-T) a'
  \*/
  template <typename t_Value>
  class KKT : public LinearSystemSolver<t_Value> {
  public:
    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    typedef LinearSystemSolver<t_Value> LSS ;

    Malloc<valueType> allocReals ;

    LSS const * pAsolver ;
    LU<t_Value> A_lu ;
    LU<t_Value> W_lu ;

    valuePointer Zmat, Cmat ;

    BandedLU<valueType> banded_LU;

  private:

    // A is n x n
    // B is n x m
    // C is m x n
    // D is m x m
    // m << n
    integer n ;
    integer m ;

    /*
    //
    //  Matrix structure
    //
    //     n     m
    //  +-----+-----+
    //  |  A  |  B  | n
    //  +-----+-----+
    //  |  C  |  D  | m
    //  +-----+-----+
    */

  public:

    explicit KKT() ;
    ~KKT();

    void
    allocate( integer n, integer m ) ;

    /*!
      Load A block as a factorized Matrix class
      \param[in] Asystem the pointer to the object containing the factorized matrix
    \*/
    void
    load_A( LSS const * Asystem ) ;

    /*!
      Load A block in sparse form

      \param[in] A_values     vector with the values
      \param[in] A_row        vector with the row index of the values
      \param[in] r_offs       offset to be applied, `A_row[j]+r_offs` is the index
      \param[in] A_col        vector with the column index of the values
      \param[in] c_offs       offset to be applied, `A_col[j]+c_offs` is the index
      \param[in] A_nnz        number of element of the sparse matrix
      \param[in] is_symmetric true if matrix is symmetric and only lower or upper parte is passed
    \*/
    void
    load_A( valueConstPointer A_values,
            integer const *   A_row, integer r_offs,
            integer const *   A_col, integer c_offs,
            integer           A_nnz,
            bool  is_symmetric = false ) ;
    
    /*!
      Load A block as a FORTRAN matrix
      \param[in] A          the pointer to the matrix
      \param[in] ldA        leading dimension of `A`
      \param[in] transposed true if the matrix must be loaded transposed
    \*/
    void
    load_A( valueConstPointer A, integer ldA, bool transposed = false ) ;

    /*!
      Load B block in sparse form

      \param[in] B_values vector with the values
      \param[in] B_row    vector with the row index of the values
      \param[in] r_offs   offset to be applied, `B_row[j]+r_offs` is the index
      \param[in] B_col    vector with the column index of the values
      \param[in] c_offs   offset to be applied, `B_col[j]+c_offs` is the index
      \param[in] B_nnz    number of element of the sparse matrix
    \*/
    void
    load_B( valueConstPointer B_values,
            integer const *   B_row, integer r_offs,
            integer const *   B_col, integer c_offs,
            integer           B_nnz ) ;
    
    /*!
      Load B block as a FORTRAN matrix
      \param[in] B          the pointer to the matrix
      \param[in] ldB        leading dimension of `B`
      \param[in] transposed true if the matrix must be loaded transposed
    \*/
    void
    load_B( valueConstPointer B, integer ldB, bool transposed = false ) ;

    /*!
      Load C block in sparse form

      \param[in] C_values vector with the values
      \param[in] C_row    vector with the row index of the values
      \param[in] r_offs   offset to be applied, `C_row[j]+r_offs` is the index
      \param[in] C_col    vector with the column index of the values
      \param[in] c_offs   offset to be applied, `C_col[j]+c_offs` is the index
      \param[in] C_nnz    number of element of the sparse matrix
    \*/
    void
    load_C( valueConstPointer C_values,
            integer const *   C_row, integer r_offs,
            integer const *   C_col, integer c_offs,
            integer           C_nnz ) ;
    
    /*!
      Load C block as a FORTRAN matrix
      \param[in] C          the pointer to the matrix
      \param[in] ldC        leading dimension of `C`
      \param[in] transposed true if the matrix must be loaded transposed
    \*/
    void
    load_C( valueConstPointer C, integer ldC, bool transposed = false ) ;

    /*!
      Load D block in sparse form

      \param[in] D_values     vector with the values
      \param[in] D_row        vector with the row index of the values
      \param[in] r_offs       offset to be applied, `D_row[j]+r_offs` is the index
      \param[in] D_col        vector with the column index of the values
      \param[in] c_offs       offset to be applied, `D_col[j]+c_offs` is the index
      \param[in] D_nnz        number of element of the sparse matrix
      \param[in] is_symmetric true if matrix `D` is symmetric and only lower or upper part is passed
    \*/
    void
    load_D( valueConstPointer D_values,
            integer const *   D_row, integer r_offs,
            integer const *   D_col, integer c_offs,
            integer           D_nnz,
            bool is_symmetric = false ) ;
    
    /*!
      Load D block as a FORTRAN matrix
      \param[in] D          the pointer to the matrix
      \param[in] ldD        leading dimension of `D`
      \param[in] transposed true if the matrix must be loaded transposed
    \*/
    void
    load_D( valueConstPointer D, integer ldD, bool transposed = false ) ;

    void factorize() ;

    //! load matrix in the class
    /*!
      \param _n  number of row and column of the first block
      \param _m  number of rows of block `C` and columns of block `B`

      \param[in] A_values elements `Aij` of the matrix `A`
      \param[in] A_row    row index of the corresponding element in `A_values`
      \param[in] A_col    column index of the corresponding element in `A_values`
      \param[in] A_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Ar_offs  offset to be applied, `A_row[j]+Ar_offs` is the index
      \param[in] Ac_offs  offset to be applied, `A_col[j]+Ac_offs` is the index
      \param[in] A_is_symmetric true if matrix `A` is symmetric and only lower or upper part is passed

      \param[in] B_values elements `Bij` of the matrix `B`
      \param[in] B_row    row index of the corresponding element in `B_values`
      \param[in] B_col    column index of the corresponding element in `B_values`
      \param[in] B_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Br_offs  offset to be applied, `B_row[j]+Br_offs` is the index
      \param[in] Bc_offs  offset to be applied, `B_col[j]+Bc_offs` is the index

      \param[in] C_values elements `Cij` of the matrix `C`
      \param[in] C_row    row index of the corresponding element in `C_values`
      \param[in] C_col    column index of the corresponding element in `C_values`
      \param[in] C_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Cr_offs  offset to be applied, `C_row[j]+Cr_offs` is the index
      \param[in] Cc_offs  offset to be applied, `C_col[j]+Cc_offs` is the index

      \param[in] D_values elements `Dij` of the matrix `D`
      \param[in] D_row    row index of the corresponding element in `D_values`
      \param[in] D_col    column index of the corresponding element in `D_values`
      \param[in] D_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Dr_offs  offset to be applied, `D_row[j]+Dr_offs` is the index
      \param[in] Dc_offs  offset to be applied, `D_col[j]+Dc_offs` is the index
      \param[in] D_is_symmetric true if matrix `D` is symmetric and only lower or upper part is passed
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               valueConstPointer A_values,
               integer const *   A_row, integer Ar_offs,
               integer const *   A_col, integer Ac_offs,
               integer           A_nnz,
               bool              A_is_symmetric,
               // -----------------------
               valueConstPointer B_values,
               integer const *   B_row, integer Br_offs,
               integer const *   B_col, integer Bc_offs,
               integer           B_nnz,
               // -----------------------
               valueConstPointer C_values,
               integer const *   C_row, integer Cr_offs,
               integer const *   C_col, integer Cc_offs,
               integer           C_nnz,
               // -----------------------
               valueConstPointer D_values,
               integer const *   D_row, integer Dr_offs,
               integer const *   D_col, integer Dc_offs,
               integer           D_nnz,
               bool              D_is_symmetric ) ;

    //! load matrix in the class
    /*!
      \param _n           # of row and column of the first block
      \param _m           # of rows of block C and columns o B

      \param A_values     elements `Aij` of the matrix `A`
      \param ldA          leading dimension of matrix `A`
      \param A_transposed true if matrix `A` must be loaded transposed

      \param B_values     elements `Bij` of the matrix `B`
      \param ldB          leading dimension of matrix `B`
      \param B_transposed true if matrix `B` must be loaded transposed

      \param C_values     elements `Cij` of the matrix `C`
      \param ldC          leading dimension of matrix `C`
      \param C_transposed true if matrix `C` must be loaded transposed

      \param D_values     elements `Dij` of the matrix `D`
      \param ldD          leading dimension of matrix `D`
      \param D_transposed true if matrix `D` must be loaded transposed
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               valueConstPointer A_values,
               integer           ldA,
               bool              A_transposed,
               // -----------------------
               valueConstPointer B_values,
               integer           ldB,
               bool              B_transposed,
               // -----------------------
               valueConstPointer C_values,
               integer           ldC,
               bool              C_transposed,
               // -----------------------
               valueConstPointer D_values,
               integer           ldD,
               bool              D_transposed ) ;

    //! load matrix in the class
    /*!
      \param[in] _n number row and column of the first block
      \param[in] _m number of rows of block `C` and columns of block `B`

      \param[in] Asystem  pointer to a class with `solve` and `t_solve` method

      \param[in] B_values elements `Bij` of the matrix `B`
      \param[in] B_row    row index of the corresponding element in `B_values`
      \param[in] B_col    column index of the corresponding element in `B_values`
      \param[in] B_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Br_offs  offset to be applied, `B_row[j]+Br_offs` is the index
      \param[in] Bc_offs  offset to be applied, `B_col[j]+Bc_offs` is the index

      \param[in] C_values elements `Cij` of the matrix `C`
      \param[in] C_row    row index of the corresponding element in `C_values`
      \param[in] C_col    column index of the corresponding element in `C_values`
      \param[in] C_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Cr_offs  offset to be applied, `C_row[j]+Cr_offs` is the index
      \param[in] Cc_offs  offset to be applied, `C_col[j]+Cc_offs` is the index

      \param[in] D_values elements `Dij` of the matrix `D`
      \param[in] D_row    row index of the corresponding element in `D_values`
      \param[in] D_col    column index of the corresponding element in `D_values`
      \param[in] D_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] Dr_offs  offset to be applied, `D_row[j]+Dr_offs` is the index
      \param[in] Dc_offs  offset to be applied, `D_col[j]+Dc_offs` is the index
      \param[in] D_is_symmetric true if matrix `D` is symmetric and only lower or upper part is passed
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               LSS     const *   Asystem,
               // -----------------------
               valueConstPointer B_values,
               integer const *   B_row, integer Br_offs,
               integer const *   B_col, integer Bc_offs,
               integer           B_nnz,
               // -----------------------
               valueConstPointer C_values,
               integer const *   C_row, integer Cr_offs,
               integer const *   C_col, integer Cc_offs,
               integer           C_nnz,
               // -----------------------
               valueConstPointer D_values,
               integer const *   D_row, integer Dr_offs,
               integer const *   D_col, integer Dc_offs,
               integer           D_nnz,
               bool              D_is_symmetric ) ;

    //! load matrix in the class
    /*!
      \param _n number of row and column of the first block
      \param _m number of rows of block `C` and columns of block `B`

      \param Asystem  pointer to a class with `solve` and `t_solve` method

      \param B_values     elements `Bij` of the matrix `B`
      \param ldB          leading dimension of matrix `B`
      \param B_transposed true if matrix `B` must be loaded transposed

      \param C_values     elements `Cij` of the matrix `C`
      \param ldC          leading dimension of matrix `C`
      \param C_transposed true if matrix `C` must be loaded transposed

      \param D_values     elements `Dij` of the matrix `D`
      \param ldD          leading dimension of matrix `D`
      \param D_transposed true if matrix `D` must be loaded transposed
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               LSS const *       Asystem,
               // -----------------------
               valueConstPointer B_values,
               integer           ldB,
               bool              B_transposed,
               // -----------------------
               valueConstPointer C_values,
               integer           ldC,
               bool              C_transposed,
               // -----------------------
               valueConstPointer D_values,
               integer           ldD,
               bool              D_transposed ) ;

    //! load matrix in the class
    /*!
      \param _n  number of row and column of the first block
      \param _m  number of rows of block `C` and columns of block `B`
      \param _nL number of extra lower diagonals
      \param _nU number of extra upper diagonals

      \param[in] M_values elements `Aij` of the matrix `A`
      \param[in] M_row    row index of the corresponding element in `A_values`
      \param[in] M_col    column index of the corresponding element in `A_values`
      \param[in] M_nnz    number of entries in the vectors `values`, `row` and `col`
      \param[in] r_offs  offset to be applied, `A_row[j]+Ar_offs` is the index
      \param[in] c_offs  offset to be applied, `A_col[j]+Ac_offs` is the index
      \param[in] M_is_symmetric true if matrix `A` is symmetric and only lower or upper part is passed
    */
    void
    factorize( integer           _n,
               integer           _m,
               integer           _nL,
               integer           _nU,
               // -----------------------
               valueConstPointer M_values,
               integer const *   M_row, integer r_offs,
               integer const *   M_col, integer c_offs,
               integer           M_nnz,
               bool              M_is_symmetric ) ;

    // -------------------------------------------------------------------------
    // virtuals redefined

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

  } ;

  // explicit instantiation declaration to suppress warnings
  extern template class KKT<float> ;
  extern template class KKT<double> ;

}

#ifdef __GCC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
