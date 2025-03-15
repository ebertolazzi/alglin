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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Alglin_Eigen.hh
///

#pragma once

#ifndef ALGLIN_EIGEN_dot_HH
#define ALGLIN_EIGEN_dot_HH

#ifdef NO_SYSTEM_UTILS
  #include "Utils_eigen.hh"
#else
  #include <Utils_eigen.hh>
#endif

namespace alglin {

  template <typename t_Value> using EigenVector = Eigen::Matrix<t_Value,Eigen::Dynamic,1>;
  template <typename t_Value> using EigenMatrix = Eigen::Matrix<t_Value,Eigen::Dynamic,Eigen::Dynamic>;
  template <typename t_Value> using MapVector   = Eigen::Map<EigenVector<t_Value>>;
  template <typename t_Value> using MapMatrix   = Eigen::Map<EigenMatrix<t_Value>>;

  template <typename t_Value>
  inline
  void
  Copy_n( t_Value const * from, integer n, t_Value * to ) {
    MapVector<t_Value> TO( to, n );
    MapVector<t_Value> FROM( const_cast<t_Value*>(from), n );
    TO.noalias() = FROM;
  }

  template <typename t_Value>
  inline
  void
  Copy_negate_n( t_Value const * from, integer n, t_Value * to ) {
    MapVector<t_Value> TO( to, n );
    MapVector<t_Value> FROM( const_cast<t_Value*>(from), n );
    TO.noalias() = -FROM;
  }

  template <typename t_Value>
  inline
  void
  Axpy_n( integer n, t_Value a, t_Value const * from, t_Value * to ) {
    MapVector<t_Value> TO( to, n );
    MapVector<t_Value> FROM( const_cast<t_Value*>(from), n );
    TO.noalias() += a*FROM;
  }

  template <typename t_Value>
  inline
  void
  Fill_n( t_Value * to, integer n, t_Value v )
  { MapVector<t_Value>( to, n ).fill(v); }

  template <typename t_Value>
  inline
  void
  Zero_n( t_Value * to, integer n )
  { MapVector<t_Value> ( to, n ).setZero(); }

  template <typename t_Value>
  void
  GEcopy(
    integer         N,
    integer         M,
    t_Value const * Amat, integer ldA,
    t_Value       * Bmat, integer ldB
  );

  extern template void GEcopy(
    integer       N,
    integer       M,
    float const * Amat, integer ldA,
    float       * Bmat, integer ldB
  );
  extern template void GEcopy(
    integer        N,
    integer        M,
    double const * Amat, integer ldA,
    double       * Bmat, integer ldB
  );

  template <typename t_Value>
  void
  GEzero(
    integer   N,
    integer   M,
    t_Value * Amat, integer ldA
  );

  extern template void GEzero(
    integer N, integer M, float * Amat, integer ldA
  );
  extern template void GEzero(
    integer N, integer M, double * Amat, integer ldA
  );

  template <typename t_Value>
  void
  GEadd(
    integer         N,
    integer         M,
    t_Value const * Amat, integer ldA,
    t_Value       * Bmat, integer ldB
  );

  extern template void GEadd(
    integer       N,
    integer       M,
    float const * Amat, integer ldA,
    float       * Bmat, integer ldB
  );
  extern template void GEadd(
    integer        N,
    integer        M,
    double const * Amat, integer ldA,
    double       * Bmat, integer ldB
  );

} // end namespace alglin

#endif

///
/// eof: Alglin_Eigen.hh
///
