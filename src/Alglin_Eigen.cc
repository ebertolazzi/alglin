/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
/// file: Alglin_Eigen.cc
///

#include "Alglin.hh"
#include "Alglin_Eigen.hh"

namespace alglin {

  template <typename t_Value>
  void
  GEcopy(
    integer         N,
    integer         M,
    t_Value const * Amat, integer ldA,
    t_Value       * Bmat, integer ldB
  ) {
    MapMatrix<t_Value> A( const_cast<t_Value*>(Amat), ldA, M );
    MapMatrix<t_Value> B( Bmat, ldB, M );
    integer kk = (ldA == N?0:1)|(ldB == N?0:2);
    switch ( kk ) {
    case 0: B.noalias()            = A;            break;
    case 1: B.noalias()            = A.topRows(N); break;
    case 2: B.topRows(N).noalias() = A;            break;
    case 3: B.topRows(N).noalias() = A.topRows(N); break;
    }
  }

  template
  void
  GEcopy(
    integer        N,
    integer        M,
    double const * Amat, integer ldA,
    double       * Bmat, integer ldB
  );

  template
  void
  GEcopy(
    integer       N,
    integer       M,
    float const * Amat, integer ldA,
    float       * Bmat, integer ldB
  );

  // --------------------------------------

  template <typename t_Value>
  void
  GEzero(
    integer   N,
    integer   M,
    t_Value * Amat, integer ldA
  ) {
    if ( N == ldA ) {
      Zero_n( Amat, N*M );
    } else {
      MapMatrix<t_Value> A( Amat, ldA, M );
      A.topRows(N).setZero();
    }
  }

  template void GEzero( integer N, integer M, double * Amat, integer ldA );
  template void GEzero( integer N, integer M, float * Amat, integer ldA );

  // --------------------------------------

  template <typename t_Value>
  void
  GEadd(
    integer         N,
    integer         M,
    t_Value const * Amat, integer ldA,
    t_Value       * Bmat, integer ldB
  ) {
    MapMatrix<t_Value> A( const_cast<t_Value*>(Amat), ldA, M );
    MapMatrix<t_Value> B( Bmat, ldB, M );
    integer kk = (ldA == N?0:1)|(ldB == N?0:2);
    switch ( kk ) {
    case 0: B.noalias()            += A;            break;
    case 1: B.noalias()            += A.topRows(N); break;
    case 2: B.topRows(N).noalias() += A;            break;
    case 3: B.topRows(N).noalias() += A.topRows(N); break;
    }
  }

  template
  void
  GEadd(
    integer        N,
    integer        M,
    double const * Amat, integer ldA,
    double       * Bmat, integer ldB
  );

  template
  void
  GEadd(
    integer       N,
    integer       M,
    float const * Amat, integer ldA,
    float       * Bmat, integer ldB
  );

} // end namespace alglin

///
/// eof: Alglin_Eigen.cc
///
