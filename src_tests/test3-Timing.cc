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


#include <iostream>
#include <vector>
#include <random>

#include "Alglin_Config.hh"
#include "Alglin_tmpl.hh"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wdeprecated"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wunused-template"
#endif

#ifdef ALGLIN_USE_SYSTEM_EIGEN
  #include <Eigen/Dense>
#else
  #include "Eigen/Dense"
#endif

using namespace std;
typedef double valueType;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

using namespace alglin;
typedef Eigen::Matrix<valueType,Eigen::Dynamic,Eigen::Dynamic> dmat_t;
typedef Eigen::Matrix<valueType,Eigen::Dynamic,1>              dvec_t;

template <int N>
void
testMM() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  typedef Eigen::Matrix<valueType,N,N> matN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<valueType>       baseValue("real");
  Malloc<alglin::integer> baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  valueType * M1 = baseValue(N*N);
  valueType * M2 = baseValue(N*N);
  valueType * M3 = baseValue(N*N);

  matN_t m1, m2, m3;
  dmat_t dm1, dm2, dm3;

  dm1.resize(N,N);
  dm2.resize(N,N);
  dm3.resize(N,N);

  for ( int i = 0; i < N; ++i ) {
    for ( int j = 0; j < N; ++j ) {
      m1(i,j) = dm1(i,j) = M1[i+j*N] = rand(-1,1);
      m2(i,j) = dm2(i,j) = M2[i+j*N] = rand(-1,1);
      m3(i,j) = dm3(i,j) = M3[i+j*N] = rand(-1,1);
    }
  }

  TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    gemm(
      NO_TRANSPOSE, NO_TRANSPOSE,
      N, N, N,
      -1.0, M1, N,
      M2, N,
      1.0, M3, N
    );
    copy( N*N, M3, 1, M2, 1);
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (lapack)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    dm3.noalias() -= dm1*dm2;
    dm2 = dm3;
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (eigen dynamic)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Eigen::Map<dmat_t> mm1(M1,N,N);
    Eigen::Map<dmat_t> mm2(M2,N,N);
    Eigen::Map<dmat_t> mm3(M3,N,N);
    mm3.noalias() -= mm1*mm2;
    mm2 = mm3;
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (eigen map dynamic)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    m3.noalias() -= m1*m2;
    m2 = m3;
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (eigen fixed)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Eigen::Map<matN_t> mm1(M1);
    Eigen::Map<matN_t> mm2(M2);
    Eigen::Map<matN_t> mm3(M3);
    mm3.noalias() -= mm1*mm2;
    mm2 = mm3;
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (eigen fixed map)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    MM<valueType,N,N,N,N,N,N>::subTo(M1,M2,M3);
    memcpy( M2, M3, N*N*sizeof(valueType) );
    //Vec2<valueType,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MM) MULT = {:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  fmt::print("All done!\n");
}


template <int N>
void
testMv() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  typedef Eigen::Matrix<valueType,N,N> matN_t;
  typedef Eigen::Matrix<valueType,N,1> vecN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<valueType>       baseValue("real");
  Malloc<alglin::integer> baseIndex("integer");

  baseValue.allocate(N*N*10);
  baseIndex.allocate(N*10);

  valueType * M = baseValue(N*N);
  valueType * V = baseValue(N);
  valueType * R = baseValue(N);

  matN_t m;
  dmat_t dm;

  vecN_t v,  r;
  dvec_t dv, dr;

  dm.resize(N,N);
  dv.resize(N);
  dr.resize(N);

  for ( int i = 0; i < N; ++i ) {
    dv(i) = v(i) = V[i] = rand(-1,1);
    dr(i) = r(i) = R[i] = rand(-1,1);
    for ( int j = 0; j < N; ++j ) {
      m(i,j) = dm(i,j) = M[i+j*N] = rand(-1,1);
    }
  }

  TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    gemv(
      NO_TRANSPOSE,
      N, N,
      -1.0, M, N,
      V, 1,
      1.0, R, 1
    );
    copy( N, R, 1, V, 1);
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (lapack)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    dr.noalias() -= dm*dv;
    dv = dr;
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (eigen dynamic)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Eigen::Map<dmat_t> mm(M,N,N);
    Eigen::Map<dvec_t> vv(V,N);
    Eigen::Map<dvec_t> rr(R,N);
    rr.noalias() -= mm*vv;
    vv = rr;
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (eigen map dynamic)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    r.noalias() -= m*v;
    v = r;
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (eigen fixed)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Eigen::Map<matN_t> mm(M);
    Eigen::Map<vecN_t> vv(V);
    Eigen::Map<vecN_t> rr(R);
    rr.noalias() -= mm*vv;
    vv = rr;
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (eigen fixed map)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    Mv<valueType,N,N,N,N,N>::subTo(M,V,R);
    memcpy( V, R, N*sizeof(valueType) );
    //Vec2<valueType,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MV) MULT = {:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  fmt::print("All done!\n");
}



int
main() {

  testMv<2>();
  testMv<3>();
  testMv<4>();
  testMv<5>();
  testMv<6>();
  testMv<7>();
  testMv<8>();
  testMv<8>();
  testMv<10>();
  testMv<11>();
  testMv<12>();
  testMv<13>();
  testMv<14>();
  testMv<15>();
  testMv<16>();
  testMv<17>();
  testMv<18>();
  testMv<19>();
  testMv<20>();

  testMM<2>();
  testMM<3>();
  testMM<4>();
  testMM<5>();
  testMM<6>();
  testMM<7>();
  testMM<8>();
  testMM<8>();
  testMM<10>();
  testMM<11>();
  testMM<12>();
  testMM<13>();
  testMM<14>();
  testMM<15>();
  testMM<16>();
  testMM<17>();
  testMM<18>();
  testMM<19>();
  testMM<20>();

  fmt::print("\n\nAll done!\n");

  return 0;
}
