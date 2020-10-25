/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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

#include "Alglin.hh"

namespace alglin {

  using namespace std;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::allocate( integer _n, integer _m ) {
    UTILS_ASSERT(
      _n > 0 && _m > 0, "KKT::allocate( {}, {} ) bad dimension\n", _n, _m
    );
    if ( n != _n || m != _m ) {
      n = _n;
      m = _m;
      allocReals.allocate( size_t(2*n*m) );
      Zmat = allocReals( size_t(n*m) );
      Cmat = allocReals( size_t(n*m) );
      W_LU.allocate(m,m);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_A( LSS const * Asystem ) {
    pAsolver = Asystem;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_A(
    valueType const A_values[],
    integer   const A_row[], integer r_offs,
    integer   const A_col[], integer c_offs,
    integer         A_nnz,
    bool            is_symmetric
  ) {
    pAsolver = &A_LU;
    A_LU_working.setup(n,n);
    if ( is_symmetric )
      A_LU_working.load_symmetric( r_offs, c_offs, A_row, A_col, A_values, A_nnz );
    else
      A_LU_working.load( r_offs, c_offs, A_row, A_col, A_values, A_nnz );
    A_LU.factorize( "KKT::load_A", A_LU_working );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_A(
    valueType const A[],
    integer         ldA,
    bool            transposed
  ) {
    UTILS_ASSERT(
      ldA >= n, "KKT::load_A bad ldA = {} must be >= {}\n", ldA, n
    );
    pAsolver = &A_LU;
    A_LU_working.setup(n,n);
    if ( transposed ) {
      for ( integer i = 0; i < n; ++i )
        A_LU_working.load_row(A+i*ldA,i);
    } else {
      A_LU_working.load_block(n,n,A,ldA);
    }
    A_LU.factorize( "KKT::load_A\n", A_LU_working );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // n x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_B(
    valueType const B_values[],
    integer   const B_row[], integer r_offs,
    integer   const B_col[], integer c_offs,
    integer         B_nnz
  ) {
     gezero( n, m, Zmat, n );
     for ( integer k = 0; k < B_nnz; ++k ) {
       integer i = B_row[k]+r_offs;
       integer j = B_col[k]+c_offs;
       UTILS_ASSERT(
         i >= 0 && i < n && j >= 0 && j < m,
         "KKT::load_B bad index (i,j) = ({},{}) at position {}\n", i, j, k
       );
       Zmat[ i + n * j ] += B_values[k];
     }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // n x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_B(
    valueType const B[],
    integer         ldB,
    bool            transposed
  ) {
    if ( transposed ) {
      UTILS_ASSERT(
        ldB >= m, "KKT::load_B bad ldB = {} must be >= {}\n", ldB,  m
      );
      for ( integer i = 0; i < n; ++i )
        copy( m, B+i, n, Zmat+i*n, 1 );
    } else {
      integer info = gecopy( n, m, B, ldB, Zmat, n );
      UTILS_ASSERT( info == 0, "KKT::load_B bad call gecopy, info = {}\n", info );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // m x n
  template <typename t_Value>
  void
  KKT<t_Value>::load_C(
    valueType const C_values[],
    integer   const C_row[], integer r_offs,
    integer   const C_col[], integer c_offs,
    integer         C_nnz
  ) {
    gezero( m, n, Cmat, m );
    for ( integer k = 0; k < C_nnz; ++k ) {
      integer i = C_row[k]+r_offs;
      integer j = C_col[k]+c_offs;
      UTILS_ASSERT(
        i >= 0 && i < m && j >= 0 && j < n,
        "KKT::load_C bad index (i,j) = ({},{}) at position {}\n", i, j, k
      );
      Cmat[ i + m * j ] += C_values[k];
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // m x n
  template <typename t_Value>
  void
  KKT<t_Value>::load_C(
    valueType const C[],
    integer         ldC,
    bool            transposed
  ) {
    if ( transposed ) {
      UTILS_ASSERT(
        ldC >= n, "KKT::load_C bad ldC = {} must be >= {}\n", ldC, n
      );
      for ( integer i = 0; i < m; ++i )
        copy( n, C+i, m, Cmat+i*m, 1 );
    } else {
      integer info = gecopy( m, n, C, ldC, Cmat, m );
      UTILS_ASSERT(
        info == 0, "KKT::load_C bad call gecopy, info = {}\n", info
      );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // m x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_D(
    valueType const D_values[],
    integer   const D_row[], integer r_offs,
    integer   const D_col[], integer c_offs,
    integer         D_nnz,
    bool            is_symmetric_D
  ) {
    W_LU_working.setup(m,m);
    W_LU_working.zero_fill();
    for ( integer k = 0; k < D_nnz; ++k ) {
      integer i = D_row[k]+r_offs;
      integer j = D_col[k]+c_offs;
      UTILS_ASSERT(
        i >= 0 && i < m && j >= 0 && j < m,
        "KKT::load_D bad index (i,j) = ({},{}) at position {}\n", i, j, k
      );
      W_LU_working(i,j) += D_values[k];
      if ( is_symmetric_D && i != j ) W_LU_working(j,i) += D_values[k];
     }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // m x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_D(
    valueType const D[],
    integer         ldD,
    bool            transposed
  ) {
    UTILS_ASSERT( ldD >= m, "KKT::load_D bad ldD = {} must be >= {}\n", ldD, m );
    W_LU_working.setup(m,m);
    if ( transposed ) W_LU_working.load_transposed( D, ldD );
    else              W_LU_working.load( D, ldD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::factorize() {
    // Compute aux matrix
    // Z = A^(-1)*B
    // W = C*Z - D
    pAsolver->solve( m, Zmat, n );
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      m, m, n,
      1,
      Cmat, m,
      Zmat, n,
      -1,
      W_LU_working.data(), m
    );
    W_LU.factorize( "KKT<::factorize<W>", W_LU_working );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer         _n,
    integer         _m,
    // -----------------------
    valueType const A_values[],
    integer   const A_row[], integer Ar_offs,
    integer   const A_col[], integer Ac_offs,
    integer         A_nnz,
    bool            A_is_symmetric,
    // -----------------------
    valueType const B_values[],
    integer   const B_row[], integer Br_offs,
    integer   const B_col[], integer Bc_offs,
    integer         B_nnz,
    // -----------------------
    valueType const C_values[],
    integer   const C_row[], integer Cr_offs,
    integer   const C_col[], integer Cc_offs,
    integer         C_nnz,
    // -----------------------
    valueType const D_values[],
    integer   const D_row[], integer Dr_offs,
    integer   const D_col[], integer Dc_offs,
    integer         D_nnz,
    bool            D_is_symmetric
  ) {
    allocate( _n, _m );
    load_A( A_values, A_row, Ar_offs, A_col, Ac_offs, A_nnz, A_is_symmetric );
    load_B( B_values, B_row, Br_offs, B_col, Bc_offs, B_nnz );
    load_C( C_values, C_row, Cr_offs, C_col, Cc_offs, C_nnz );
    load_D( D_values, D_row, Dr_offs, D_col, Dc_offs, D_nnz, D_is_symmetric );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer         _n,
    integer         _m,
    // -----------------------
    valueType const A_values[],
    integer         ldA,
    bool            A_transposed,
    // -----------------------
    valueType const B_values[],
    integer         ldB,
    bool            B_transposed,
    // -----------------------
    valueType const C_values[],
    integer         ldC,
    bool            C_transposed,
    // -----------------------
    valueType const D_values[],
    integer         ldD,
    bool            D_transposed
  ) {
    allocate( _n, _m );
    load_A( A_values, ldA, A_transposed );
    load_B( B_values, ldB, B_transposed );
    load_C( C_values, ldC, C_transposed );
    load_D( D_values, ldD, D_transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer _n,
    integer _m,
    // -----------------------
    LSS const * Asystem,
    // -----------------------
    valueType const B_values[],
    integer   const B_row[], integer Br_offs,
    integer   const B_col[], integer Bc_offs,
    integer         B_nnz,
    // -----------------------
    valueType const C_values[],
    integer   const C_row[], integer Cr_offs,
    integer   const C_col[], integer Cc_offs,
    integer         C_nnz,
    // -----------------------
    valueType const D_values[],
    integer   const D_row[], integer Dr_offs,
    integer   const D_col[], integer Dc_offs,
    integer         D_nnz,
    bool            D_is_symmetric
  ) {
    allocate( _n, _m );
    load_A( Asystem );
    load_B( B_values, B_row, Br_offs, B_col, Bc_offs, B_nnz );
    load_C( C_values, C_row, Cr_offs, C_col, Cc_offs, C_nnz );
    load_D( D_values, D_row, Dr_offs, D_col, Dc_offs, D_nnz, D_is_symmetric );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer _n,
    integer _m,
    // -----------------------
    LSS const * Asystem,
    // -----------------------
    valueType const B_values[],
    integer         ldB,
    bool            B_transposed,
    // -----------------------
    valueType const C_values[],
    integer         ldC,
    bool            C_transposed,
    // -----------------------
    valueType const D_values[],
    integer         ldD,
    bool            D_transposed
  ) {
    allocate( _n, _m );
    load_A( Asystem );
    load_B( B_values, ldB, B_transposed );
    load_C( C_values, ldC, C_transposed );
    load_D( D_values, ldD, D_transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_banded(
    integer         _n,
    integer         _m,
    integer         _nL,
    integer         _nU,
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( _n, _m );
    A_banded_LU.setup( _n, _n, _nL, _nU );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      A_banded_LU
    );
    pAsolver = &A_banded_LU;
    A_banded_LU.factorize( "KKT::load_banded" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_triblock(
    integer         _n,
    integer         _m,
    // ---- BLOCK TRIDIAGONAL STRUCTURE ----
    integer         _nblocks,
    integer const   rBlocks[],
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( _n, _m );
    A_strid_LDL.setup( _nblocks, rBlocks );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      A_strid_LDL
    );
    pAsolver = &A_strid_LDL;
    A_strid_LDL.factorize( "KKT::load_triblock" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_triblock(
    integer _n,
    integer _m,
    // ---- BLOCK TRIDIAGONAL STRUCTURE ----
    integer _nblocks,
    integer _block_size,
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( _n, _m );
    A_strid_LDL.setup( _nblocks, _block_size );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      A_strid_LDL
    );
    pAsolver = &A_strid_LDL;
    A_strid_LDL.factorize( "KKT::load_triblock" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::solve( valueType xb[] ) const {
    // a' = A^(-1)*a
    pAsolver->solve( xb );
    // b' = C*a' - b
    gemv(
      NO_TRANSPOSE,
      m, n,
      1, Cmat, m,
      xb, 1,
      -1, xb+n, 1
    );
    // y = W^(-1) * b'
    W_LU.solve( xb+n );
    // x = a' - Z*y
    gemv(
      NO_TRANSPOSE,
      n, m,
      -1, Zmat, n,
      xb+n, 1,
      1, xb, 1
    );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::t_solve( valueType xb[] ) const {
    // b' = Z^T*a -b
    gemv(
      TRANSPOSE,
      n, m,
      1, Zmat, n,
      xb, 1,
      -1, xb+n, 1
    );
    // y  = W^(-T)*b'
    W_LU.t_solve( xb+n );
    // a' = a - C^T*y
    gemv(
      TRANSPOSE,
      m, n,
      -1, Cmat, m,
      xb+n, 1,
      1, xb, 1
    );
    // x = A^(-T) a'
    pAsolver->t_solve( xb );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // a' = A^(-1)*a
    pAsolver->solve( nrhs, B, ldB );
    // b' = C*a' - b
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      m, nrhs, n,
      1, Cmat, m,
      B, ldB,
      -1, B+n, ldB
    );
    // y = W^(-1) * b'
    W_LU.solve( nrhs, B+n, ldB );
    // x = a' - Z*y
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      n, nrhs, m,
      -1, Zmat, n,
      B+n, ldB,
      1, B, ldB
    );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    // b' = Z^T*a -b
    gemm(
      TRANSPOSE,
      NO_TRANSPOSE,
      m, nrhs, n,
      1, Zmat, n,
      B, ldB,
      -1, B+n, ldB
    );
    // y  = W^(-T)*b'
    W_LU.t_solve( nrhs, B+n, ldB );
    // a' = a - C^T*y
    gemm(
      TRANSPOSE,
      NO_TRANSPOSE,
      n, nrhs, m,
      -1, Cmat, m,
      B+n, ldB,
      1, B, ldB
    );
    // x = A^(-T) a'
    pAsolver->t_solve( nrhs, B, ldB );
    return true;
  }

  template class KKT<float>;
  template class KKT<double>;

}
