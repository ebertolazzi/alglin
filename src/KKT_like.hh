/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#pragma clang diagnostic ignored "-Wc++98-compat"
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

    void
    load_A( LSS const * Asystem ) ;

    void
    load_A( valueConstPointer A_values,
            integer const *   A_row,
            integer const *   A_col,
            integer           A_nnz ) ;
    
    void
    load_A( valueConstPointer A, integer ldA, bool transposed = false ) ;

    void
    load_B( valueConstPointer B_values,
            integer const *   B_row,
            integer const *   B_col,
            integer           B_nnz ) ;
    
    void
    load_B( valueConstPointer B, integer ldB, bool transposed = false ) ;

    void
    load_C( valueConstPointer C_values,
            integer const *   C_row,
            integer const *   C_col,
            integer           C_nnz ) ;
    
    void
    load_C( valueConstPointer C, integer ldC, bool transposed = false ) ;

    void
    load_D( valueConstPointer D_values,
            integer const *   D_row,
            integer const *   D_col,
            integer           D_nnz ) ;
    
    void
    load_D( valueConstPointer D, integer ldD, bool transposed = false ) ;

    void factorize() ;

    //! load matrix in the class
    /*!
      \param _n       # of row and column of the first block
      \param _m       # of rows of block C and columns o B

      \param A_values elements `Aij` of the matrix `A`
      \param A_row    row index of the corresponding element in `A_values`
      \param A_col    column index of the corresponding element in `A_values`
      \param A_nnz    number of entries in the vectors `values`, `row` and `col`

      \param B_values elements `Bij` of the matrix `B`
      \param B_row    row index of the corresponding element in `B_values`
      \param B_col    column index of the corresponding element in `B_values`
      \param B_nnz    number of entries in the vectors `values`, `row` and `col`

      \param C_values elements `Cij` of the matrix `C`
      \param C_row    row index of the corresponding element in `C_values`
      \param C_col    column index of the corresponding element in `C_values`
      \param C_nnz    number of entries in the vectors `values`, `row` and `col`

      \param D_values elements `Dij` of the matrix `D`
      \param D_row    row index of the corresponding element in `D_values`
      \param D_col    column index of the corresponding element in `D_values`
      \param D_nnz    number of entries in the vectors `values`, `row` and `col`
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               valueConstPointer A_values,
               integer const *   A_row,
               integer const *   A_col,
               integer           A_nnz,
               // -----------------------
               valueConstPointer B_values,
               integer const *   B_row,
               integer const *   B_col,
               integer           B_nnz,
               // -----------------------
               valueConstPointer C_values,
               integer const *   C_row,
               integer const *   C_col,
               integer           C_nnz,
               // -----------------------
               valueConstPointer D_values,
               integer const *   D_row,
               integer const *   D_col,
               integer           D_nnz ) ;

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
      \param _n        # of row and column of the first block
      \param _m        # of rows of block C and columns o B

      \param Asystem  pointer to a class with `solve` and `t_solve` method

      \param B_values elements `Bij` of the matrix `B`
      \param B_row    row index of the corresponding element in `B_values`
      \param B_col    column index of the corresponding element in `B_values`
      \param B_nnz    number of entries in the vectors `values`, `row` and `col`

      \param C_values elements `Cij` of the matrix `C`
      \param C_row    row index of the corresponding element in `C_values`
      \param C_col    column index of the corresponding element in `C_values`
      \param C_nnz    number of entries in the vectors `values`, `row` and `col`

      \param D_values elements `Dij` of the matrix `D`
      \param D_row    row index of the corresponding element in `D_values`
      \param D_col    column index of the corresponding element in `D_values`
      \param D_nnz    number of entries in the vectors `values`, `row` and `col`
    */
    void
    factorize( integer           _n,
               integer           _m,
               // -----------------------
               LSS     const *   Asystem,
               // -----------------------
               valueConstPointer B_values,
               integer const *   B_row,
               integer const *   B_col,
               integer           B_nnz,
               // -----------------------
               valueConstPointer C_values,
               integer const *   C_row,
               integer const *   C_col,
               integer           C_nnz,
               // -----------------------
               valueConstPointer D_values,
               integer const *   D_row,
               integer const *   D_col,
               integer           D_nnz ) ;

    //! load matrix in the class
    /*!
      \param _n        # of row and column of the first block
      \param _m        # of rows of block C and columns o B

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
   *  / L  0 \ / U  U' \
   *  |      | |       |
   *  \ L' I / \ 0  0  /
   *
   *  L*U*x + B*y = a
   *  C*x + D*y = b
   * 
  \*/

}

#ifdef __GCC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
