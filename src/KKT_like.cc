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
  KKT<t_Value>::allocate( integer n, integer m ) {
    UTILS_ASSERT(
      n > 0 && m > 0, "KKT::allocate( {}, {} ) bad dimension\n", n, m
    );
    if ( n != m_dim1 || m != m_dim2 ) {
      m_dim1 = n;
      m_dim2 = m;
      m_allocReals.reallocate( size_t(2*n*m) );
      m_Zmat = m_allocReals( size_t(n*m) );
      m_Cmat = m_allocReals( size_t(n*m) );
      m_W_LU.allocate(m,m);
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_A( LSS const * Asystem ) {
    m_Asolver = Asystem;
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
    integer const & n = m_dim1;
    m_Asolver = &m_A_LU;
    m_A_LU_working.setup(n,n);
    if ( is_symmetric )
      m_A_LU_working.load_symmetric( r_offs, c_offs, A_row, A_col, A_values, A_nnz );
    else
      m_A_LU_working.load( r_offs, c_offs, A_row, A_col, A_values, A_nnz );
    m_A_LU.factorize( "KKT::load_A", m_A_LU_working );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_A(
    valueType const A[],
    integer         ldA,
    bool            transposed
  ) {
    integer const & n = m_dim1;
    UTILS_ASSERT(
      ldA >= n, "KKT::load_A bad ldA = {} must be >= {}\n", ldA, n
    );
    m_Asolver = &m_A_LU;
    m_A_LU_working.setup(n,n);
    if ( transposed ) {
      for ( integer i = 0; i < n; ++i )
        m_A_LU_working.load_row(A+i*ldA,i);
    } else {
      m_A_LU_working.load_block(n,n,A,ldA);
    }
    m_A_LU.factorize( "KKT::load_A\n", m_A_LU_working );
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
    integer const & n = m_dim1;
    integer const & m = m_dim2;
     gezero( n, m, m_Zmat, n );
     for ( integer k = 0; k < B_nnz; ++k ) {
       integer i = B_row[k]+r_offs;
       integer j = B_col[k]+c_offs;
       UTILS_ASSERT(
         i >= 0 && i < n && j >= 0 && j < m,
         "KKT::load_B bad index (i,j) = ({},{}) at position {}\n", i, j, k
       );
       m_Zmat[ i + n * j ] += B_values[k];
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
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    if ( transposed ) {
      UTILS_ASSERT(
        ldB >= m, "KKT::load_B bad ldB = {} must be >= {}\n", ldB,  m
      );
      for ( integer i = 0; i < n; ++i )
        copy( m, B+i, n, m_Zmat+i*n, 1 );
    } else {
      integer info = gecopy( n, m, B, ldB, m_Zmat, n );
      UTILS_ASSERT(
        info == 0,
        "KKT::load_B bad call gecopy, info = {}\n", info
      );
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
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    gezero( m, n, m_Cmat, m );
    for ( integer k = 0; k < C_nnz; ++k ) {
      integer i = C_row[k]+r_offs;
      integer j = C_col[k]+c_offs;
      UTILS_ASSERT(
        i >= 0 && i < m && j >= 0 && j < n,
        "KKT::load_C bad index (i,j) = ({},{}) at position {}\n", i, j, k
      );
      m_Cmat[ i + m * j ] += C_values[k];
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
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    if ( transposed ) {
      UTILS_ASSERT(
        ldC >= n, "KKT::load_C bad ldC = {} must be >= {}\n", ldC, n
      );
      for ( integer i = 0; i < m; ++i )
        copy( n, C+i, m, m_Cmat+i*m, 1 );
    } else {
      integer info = gecopy( m, n, C, ldC, m_Cmat, m );
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
    integer const & m = m_dim2;
    m_W_LU_working.setup(m,m);
    m_W_LU_working.zero_fill();
    for ( integer k = 0; k < D_nnz; ++k ) {
      integer i = D_row[k]+r_offs;
      integer j = D_col[k]+c_offs;
      UTILS_ASSERT(
        i >= 0 && i < m && j >= 0 && j < m,
        "KKT::load_D bad index (i,j) = ({},{}) at position {}\n", i, j, k
      );
      m_W_LU_working(i,j) += D_values[k];
      if ( is_symmetric_D && i != j ) m_W_LU_working(j,i) += D_values[k];
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
    integer const & m = m_dim2;
    UTILS_ASSERT( ldD >= m, "KKT::load_D bad ldD = {} must be >= {}\n", ldD, m );
    m_W_LU_working.setup(m,m);
    if ( transposed ) m_W_LU_working.load_transposed( D, ldD );
    else              m_W_LU_working.load( D, ldD );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::factorize() {
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    // Compute aux matrix
    // Z = A^(-1)*B
    // W = C*Z - D
    m_Asolver->solve( m, m_Zmat, n );
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      m, m, n,
      1,
      m_Cmat, m,
      m_Zmat, n,
      -1,
      m_W_LU_working.data(), m
    );
    m_W_LU.factorize( "KKT<::factorize<W>", m_W_LU_working );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer         n,
    integer         m,
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
    allocate( n, m );
    load_A( A_values, A_row, Ar_offs, A_col, Ac_offs, A_nnz, A_is_symmetric );
    load_B( B_values, B_row, Br_offs, B_col, Bc_offs, B_nnz );
    load_C( C_values, C_row, Cr_offs, C_col, Cc_offs, C_nnz );
    load_D( D_values, D_row, Dr_offs, D_col, Dc_offs, D_nnz, D_is_symmetric );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer         n,
    integer         m,
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
    allocate( n, m );
    load_A( A_values, ldA, A_transposed );
    load_B( B_values, ldB, B_transposed );
    load_C( C_values, ldC, C_transposed );
    load_D( D_values, ldD, D_transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer n,
    integer m,
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
    allocate( n, m );
    load_A( Asystem );
    load_B( B_values, B_row, Br_offs, B_col, Bc_offs, B_nnz );
    load_C( C_values, C_row, Cr_offs, C_col, Cc_offs, C_nnz );
    load_D( D_values, D_row, Dr_offs, D_col, Dc_offs, D_nnz, D_is_symmetric );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load(
    integer n,
    integer m,
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
    allocate( n, m );
    load_A( Asystem );
    load_B( B_values, ldB, B_transposed );
    load_C( C_values, ldC, C_transposed );
    load_D( D_values, ldD, D_transposed );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_banded(
    integer         n,
    integer         m,
    integer         nL,
    integer         nU,
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( n, m );
    m_A_banded_LU.setup( n, n, nL, nU );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      m_A_banded_LU
    );
    m_Asolver = &m_A_banded_LU;
    m_A_banded_LU.factorize( "KKT::load_banded" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_triblock(
    integer         n,
    integer         m,
    // ---- BLOCK TRIDIAGONAL STRUCTURE ----
    integer         nblocks,
    integer const   rBlocks[],
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( n, m );
    m_A_strid_LDL.setup( nblocks, rBlocks );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      m_A_strid_LDL
    );
    m_Asolver = &m_A_strid_LDL;
    m_A_strid_LDL.factorize( "KKT::load_triblock" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  void
  KKT<t_Value>::load_triblock(
    integer n,
    integer m,
    // ---- BLOCK TRIDIAGONAL STRUCTURE ----
    integer nblocks,
    integer block_size,
    // -----------------------
    valueType const M_values[],
    integer   const M_row[], integer r_offs,
    integer   const M_col[], integer c_offs,
    integer         M_nnz,
    bool            M_is_symmetric
  ) {
    allocate( n, m );
    m_A_strid_LDL.setup( nblocks, block_size );
    this->load(
      M_values, M_row, r_offs, M_col, c_offs, M_nnz, M_is_symmetric,
      m_A_strid_LDL
    );
    m_Asolver = &m_A_strid_LDL;
    m_A_strid_LDL.factorize( "KKT::load_triblock" );
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::solve( valueType xb[] ) const {
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    // a' = A^(-1)*a
    m_Asolver->solve( xb );
    // b' = C*a' - b
    gemv(
      NO_TRANSPOSE,
      m, n,
      1, m_Cmat, m,
      xb, 1,
      -1, xb+n, 1
    );
    // y = W^(-1) * b'
    m_W_LU.solve( xb+n );
    // x = a' - Z*y
    gemv(
      NO_TRANSPOSE,
      n, m,
      -1, m_Zmat, n,
      xb+n, 1,
      1, xb, 1
    );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::t_solve( valueType xb[] ) const {
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    // b' = Z^T*a -b
    gemv(
      TRANSPOSE,
      n, m,
      1, m_Zmat, n,
      xb, 1,
      -1, xb+n, 1
    );
    // y  = W^(-T)*b'
    m_W_LU.t_solve( xb+n );
    // a' = a - C^T*y
    gemv(
      TRANSPOSE,
      m, n,
      -1, m_Cmat, m,
      xb+n, 1,
      1, xb, 1
    );
    // x = A^(-T) a'
    m_Asolver->t_solve( xb );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::solve( integer nrhs, valueType B[], integer ldB ) const {
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    // a' = A^(-1)*a
    m_Asolver->solve( nrhs, B, ldB );
    // b' = C*a' - b
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      m, nrhs, n,
      1, m_Cmat, m,
      B, ldB,
      -1, B+n, ldB
    );
    // y = W^(-1) * b'
    m_W_LU.solve( nrhs, B+n, ldB );
    // x = a' - Z*y
    gemm(
      NO_TRANSPOSE,
      NO_TRANSPOSE,
      n, nrhs, m,
      -1, m_Zmat, n,
      B+n, ldB,
      1, B, ldB
    );
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename t_Value>
  bool
  KKT<t_Value>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    integer const & n = m_dim1;
    integer const & m = m_dim2;
    // b' = Z^T*a -b
    gemm(
      TRANSPOSE,
      NO_TRANSPOSE,
      m, nrhs, n,
      1, m_Zmat, n,
      B, ldB,
      -1, B+n, ldB
    );
    // y  = W^(-T)*b'
    m_W_LU.t_solve( nrhs, B+n, ldB );
    // a' = a - C^T*y
    gemm(
      TRANSPOSE,
      NO_TRANSPOSE,
      n, nrhs, m,
      -1, m_Cmat, m,
      B+n, ldB,
      1, B, ldB
    );
    // x = A^(-T) a'
    m_Asolver->t_solve( nrhs, B, ldB );
    return true;
  }

  template class KKT<float>;
  template class KKT<double>;
}
