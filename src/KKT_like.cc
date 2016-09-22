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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "KKT_like.hh"

namespace alglin {

  using namespace std ;

  template <typename t_Value>
  void
  KKT<t_Value>::allocate( integer _n, integer _m ) {
    ALGLIN_ASSERT( _n > 0 && _m > 0,
                   "KKT::allocate( " << _n << "," << _m << ")  bad dimension" ) ;
    if ( n != _n || m != _m ) {
      n = _n ;
      m = _m ;
      allocReals.allocate( 2*n*m ) ;
      Zmat = allocReals( n*m ) ;
      Cmat = allocReals( n*m ) ;
      W_lu.allocate(m,m) ;
    }
  }

  template <typename t_Value>
  void
  KKT<t_Value>::load_A( LSS const * Asystem ) {
    pAsolver = Asystem ;
  }

  template <typename t_Value>
  void
  KKT<t_Value>::load_A( valueConstPointer A_values,
                        integer const *   A_row,
                        integer const *   A_col,
                        integer           A_nnz ) {
    pAsolver = &A_lu ;
    A_lu.allocate(n,n) ;
    A_lu.load_sparse( A_nnz,
                      A_values,
                      A_row,
                      A_col ) ;
    A_lu.factorize() ;
  }
    
  template <typename t_Value>
  void
  KKT<t_Value>::load_A( valueConstPointer A,
                        integer           ldA,
                        bool transposed ) {
    ALGLIN_ASSERT( ldA >= n,
                   "KKT::load_A bad ldA = " << ldA << " must be >= " << n ) ;
    pAsolver = &A_lu ;
    A_lu.allocate(n,n) ;
    if ( transposed ) {
      for ( integer i = 0 ; i < n ; ++i )
        A_lu.load_row(A+i*ldA,i) ;
    } else {
      A_lu.load_block(n,n,A,ldA) ;
    }
    A_lu.factorize() ;
  }

  // n x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_B( valueConstPointer B_values,
                        integer const *   B_row,
                        integer const *   B_col,
                        integer           B_nnz ) {
     gezero( n, m, Zmat, n ) ;
     for ( integer k = 0 ; k < B_nnz ; ++k ) {
       integer i = B_row[k] ;
       integer j = B_col[k] ;
       ALGLIN_ASSERT( i >= 0 && i < n && j >= 0 && j < m,
                      "KKT::load_B bad index (i,j) = (" << i << "," << j <<
                      ") at position " << k ) ;
       Zmat[ i + n * j ] += B_values[k] ;
     }
  }

  // n x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_B( valueConstPointer B,
                        integer           ldB,
                        bool              transposed ) {
    if ( transposed ) {
      ALGLIN_ASSERT( ldB >= m,
                     "KKT::load_B bad ldB = " << ldB << " must be >= " << m ) ;
      for ( integer i = 0 ; i < n ; ++i )
        copy( m, B+i, n, Zmat+i*n, 1 ) ;
    } else {
     integer info = gecopy( n, m, B, ldB, Zmat, n ) ;
     ALGLIN_ASSERT( info == 0,
                    "KKT::load_B bad call gecopy, info = " << info ) ;
    }
  }

  // m x n
  template <typename t_Value>
  void
  KKT<t_Value>::load_C( valueConstPointer C_values,
                        integer const *   C_row,
                        integer const *   C_col,
                        integer           C_nnz ) {
     gezero( m, n, Cmat, m ) ;
     for ( integer k = 0 ; k < C_nnz ; ++k ) {
       integer i = C_row[k] ;
       integer j = C_col[k] ;
       ALGLIN_ASSERT( i >= 0 && i < m && j >= 0 && j < n,
                      "KKT::load_C bad index (i,j) = (" << i << "," << j <<
                      ") at position " << k ) ;
       Cmat[ i + m * j ] += C_values[k] ;
     }
  }

  // m x n
  template <typename t_Value>
  void
  KKT<t_Value>::load_C( valueConstPointer C,
                        integer           ldC,
                        bool              transposed ) {
    if ( transposed ) {
      ALGLIN_ASSERT( ldC >= n,
                     "KKT::load_C bad ldC = " << ldC << " must be >= " << n ) ;
      for ( integer i = 0 ; i < m ; ++i )
        copy( n, C+i, m, Cmat+i*m, 1 ) ;
    } else {
     integer info = gecopy( m, n, C, ldC, Cmat, m ) ;
     ALGLIN_ASSERT( info == 0,
                    "KKT::load_C bad call gecopy, info = " << info ) ;
    }
  }

  // m x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_D( valueConstPointer D_values,
                        integer const *   D_row,
                        integer const *   D_col,
                        integer           D_nnz ) {
     valuePointer Wmat = W_lu.Apointer() ;
     gezero( m, m, Wmat, m ) ;
     for ( integer k = 0 ; k < D_nnz ; ++k ) {
       integer i = D_row[k] ;
       integer j = D_col[k] ;
       ALGLIN_ASSERT( i >= 0 && i < m && j >= 0 && j < m,
                      "KKT::load_D bad index (i,j) = (" << i << "," << j <<
                      ") at position " << k ) ;
       Wmat[ i + m * j ] += D_values[k] ;
     }
  }

  // m x m
  template <typename t_Value>
  void
  KKT<t_Value>::load_D( valueConstPointer D,
                        integer           ldD,
                        bool              transposed ) {
    ALGLIN_ASSERT( ldD >= m,
                   "KKT::load_D bad ldD = " << ldD << " must be >= " << m ) ;
    valuePointer Wmat = W_lu.Apointer() ;
    if ( transposed ) {
      for ( integer i = 0 ; i < m ; ++i )
        copy( m, D+i, m, Wmat+i*m, 1 ) ;
    } else {
     integer info = gecopy( m, m, D, ldD, Wmat, m ) ;
     ALGLIN_ASSERT( info == 0,
                    "KKT::load_C bad call gecopy, info = " << info ) ;
    }
  }

  template <typename t_Value>
  void
  KKT<t_Value>::factorize() {
    // Compute aux matrix
    // Z = A^(-1)*B
    // W = C*Z - D
    pAsolver->solve( m, Zmat, n ) ;
    valuePointer Wmat = W_lu.Apointer() ;
    gemm( NO_TRANSPOSE,
          NO_TRANSPOSE,
          m, m, n,
          1,
          Cmat, m,
          Zmat, n,
          -1,
          Wmat, m ) ;
    W_lu.factorize() ;
  }
  
  template <typename t_Value>
  void
  KKT<t_Value>::solve( valueType xb[] ) const {
    // a' = A^(-1)*a
    pAsolver->solve( xb ) ;
    // b' = C*a' - b
    gemv( NO_TRANSPOSE,
          m, n,
          1, Cmat, m,
          xb, 1,
          -1, xb+n, 1 ) ;
    // y = W^(-1) * b'
    W_lu.solve( xb+n ) ;
    // x = a' - Z*y
    gemv( NO_TRANSPOSE,
          n, m,
          -1, Zmat, n,
          xb+n, 1,
          1, xb, 1 ) ;
  }

  template <typename t_Value>
  void
  KKT<t_Value>::t_solve( valueType xb[] ) const {
    // b' = Z^T*a -b
    gemv( TRANSPOSE,
          n, m,
          1, Zmat, n,
          xb, 1,
          -1, xb+n, 1 ) ;
    // y  = W^(-T)*b'
    W_lu.t_solve( xb+n ) ;
    // a' = a - C^T*y
    gemv( TRANSPOSE,
          m, n,
          -1, Cmat, m,
          xb+n, 1,
          1, xb, 1 ) ;
    // x = A^(-T) a'
    pAsolver->t_solve( xb ) ;
  }

  template <typename t_Value>
  void
  KKT<t_Value>::solve( integer nrhs, valueType B[], integer ldB ) const {
    // a' = A^(-1)*a
    pAsolver->solve( nrhs, B, ldB ) ;
    // b' = C*a' - b
    gemm( NO_TRANSPOSE,
          NO_TRANSPOSE,
          m, nrhs, n,
          1, Cmat, m,
          B, ldB,
          -1, B+n, ldB ) ;
    // y = W^(-1) * b'
    W_lu.solve( nrhs, B+n, ldB ) ;
    // x = a' - Z*y
    gemm( NO_TRANSPOSE,
          NO_TRANSPOSE,
          n, nrhs, m,
          -1, Zmat, n,
          B+n, ldB,
          1, B, ldB ) ;
  }

  template <typename t_Value>
  void
  KKT<t_Value>::t_solve( integer nrhs, valueType B[], integer ldB ) const {
    // b' = Z^T*a -b
    gemm( TRANSPOSE,
          NO_TRANSPOSE,
          m, nrhs, n,
          1, Zmat, n,
          B, ldB,
          -1, B+n, ldB ) ;
    // y  = W^(-T)*b'
    W_lu.t_solve( nrhs, B+n, ldB ) ;
    // a' = a - C^T*y
    gemm( TRANSPOSE,
          NO_TRANSPOSE,
          n, nrhs, m,
          -1, Cmat, m,
          B+n, ldB,
          1, B, ldB ) ;
    // x = A^(-T) a'
    pAsolver->t_solve( nrhs, B, ldB ) ;
  }

  template class KKT<float> ;
  template class KKT<double> ;

}
