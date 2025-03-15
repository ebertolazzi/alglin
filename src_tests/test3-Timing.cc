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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"
#include "Alglin_Eigen.hh"
#include <random>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

using namespace std;
typedef double real_type;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
real_type
rand( real_type xmin, real_type xmax ) {
  real_type random = static_cast<real_type>(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

using namespace alglin;
typedef Eigen::Matrix<real_type,Eigen::Dynamic,Eigen::Dynamic> dmat_t;
typedef Eigen::Matrix<real_type,Eigen::Dynamic,1>              dvec_t;

real_type lapack_time,
          eigen_dynamic_time,
          eigen_map_dynamic_time,
          eigen_fixed_time,
          eigen_fixed_map_time;

static
void
message() {
  real_type to_perc = lapack_time;
  //to_perc = min( to_perc, eigen_dynamic_time );
  //to_perc = min( to_perc, eigen_map_dynamic_time );
  //to_perc = min( to_perc, eigen_fixed_time );
  //to_perc = min( to_perc, eigen_fixed_map_time );
  to_perc = 100/to_perc;
  fmt::print(
    "{:8.5} [ps] {:8.5}% (lapack)\n"
    "{:8.5} [ps] {:8.5}% (eigen dynamic)\n"
    "{:8.5} [ps] {:8.5}% (eigen map dynamic)\n"
    "{:8.5} [ps] {:8.5}% (eigen fixed)\n"
    "{:8.5} [ps] {:8.5}% (eigen fixed map)\n",
    lapack_time,            to_perc*lapack_time,
    eigen_dynamic_time,     to_perc*eigen_dynamic_time,
    eigen_map_dynamic_time, to_perc*eigen_map_dynamic_time,
    eigen_fixed_time,       to_perc*eigen_fixed_time,
    eigen_fixed_map_time,   to_perc*eigen_fixed_map_time
  );
}

/*\
===============================================================================
===============================================================================
===============================================================================
===============================================================================
\*/

template <int N>
void
testCopy() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  typedef Eigen::Matrix<real_type,N,N> matN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<real_type>       base_value("real");
  Malloc<alglin::integer> base_index("integer");

  base_value.allocate(N*N*10);
  base_index.allocate(N*10);

  real_type * M1 = base_value(N*N);
  real_type * M2 = base_value(N*N);
  real_type * M3 = base_value(N*N);

  matN_t m1, m2, m3;
  dmat_t dm1, dm2, dm3;

  dm1.resize(N,N);
  dm2.resize(N,N);
  dm3.resize(N,N);

  for ( int i{0}; i < N; ++i ) {
    for ( int j{0}; j < N; ++j ) {
      m1(i,j) = dm1(i,j) = M1[i+j*N] = rand(-1,1);
      m2(i,j) = dm2(i,j) = M2[i+j*N] = rand(-1,1);
      m3(i,j) = dm3(i,j) = M3[i+j*N] = rand(-1,1);
    }
  }

  Utils::TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    alglin::gecopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
    alglin::gecopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
    alglin::gecopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
    alglin::gecopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
    alglin::gecopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
  }
  tm.toc();
  lapack_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    dm2.noalias() = dm1; dm1(0,0) += dm2(0,0);
    dm2.noalias() = dm1; dm1(0,0) += dm2(0,0);
    dm2.noalias() = dm1; dm1(0,0) += dm2(0,0);
    dm2.noalias() = dm1; dm1(0,0) += dm2(0,0);
    dm2.noalias() = dm1; dm1(0,0) += dm2(0,0);
  }
  tm.toc();
  eigen_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================
  {
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      alglin::GEcopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
      alglin::GEcopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
      alglin::GEcopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
      alglin::GEcopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
      alglin::GEcopy( N, N, M1, N, M2, N ); M1[0] += M2[0];
    }
    tm.toc();
  }
  eigen_map_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    m2.noalias() = m1; m1(0,0) += m2(0,0);
    m2.noalias() = m1; m1(0,0) += m2(0,0);
    m2.noalias() = m1; m1(0,0) += m2(0,0);
    m2.noalias() = m1; m1(0,0) += m2(0,0);
    m2.noalias() = m1; m1(0,0) += m2(0,0);
  }
  tm.toc();
  eigen_fixed_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    Eigen::Map<matN_t> mm1(M1);
    Eigen::Map<matN_t> mm2(M2);
    mm2.noalias() = mm1; mm1(0,0) += mm2(0,0);
    mm2.noalias() = mm1; mm1(0,0) += mm2(0,0);
    mm2.noalias() = mm1; mm1(0,0) += mm2(0,0);
    mm2.noalias() = mm1; mm1(0,0) += mm2(0,0);
    mm2.noalias() = mm1; mm1(0,0) += mm2(0,0);
  }
  tm.toc();
  eigen_fixed_map_time = to_ps*tm.elapsed_ms();

  fmt::print("(MM) COPY\n" );
  message();

  msg.green( "All done!\n" );
}

/*\
===============================================================================
===============================================================================
===============================================================================
===============================================================================
\*/

template <int N>
void
testVV() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  typedef Eigen::Matrix<real_type,N,1> vecN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<real_type>       base_value("real");
  Malloc<alglin::integer> base_index("integer");

  base_value.allocate(N*10);
  base_index.allocate(N*10);

  real_type * V1 = base_value(N);
  real_type * V2 = base_value(N);
  real_type * V3 = base_value(N);

  vecN_t v1, v2, v3;
  dvec_t dv1, dv2, dv3;

  dv1.resize(N);
  dv2.resize(N);
  dv3.resize(N);

  for ( int i{0}; i < N; ++i ) {
    v1(i) = dv1(i) = V1[i] = rand(-1,1);
    v2(i) = dv2(i) = V2[i] = rand(-1,1);
    v3(i) = dv3(i) = V3[i] = rand(-1,1);
  }

  Utils::TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    real_type alpha = 1.0/i;
    alglin::copy( N, V2, 1, V3, 1 );
    alglin::axpy( N, alpha, V1, 1, V3, 1 );
    alglin::Copy_n( V3, N, V1 );
  }
  tm.toc();
  lapack_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    real_type alpha = 1.0/i;
    dv3.noalias() = dv2 + alpha*dv1;
    dv1.noalias() = dv3;
  }
  tm.toc();
  eigen_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      real_type alpha = 1.0/i;
      alglin::Copy_n( V2, N, V3 );
      alglin::Axpy_n( N, alpha, V1, V3 );
      alglin::Copy_n( V3, N, V1 );
    }
    tm.toc();
  }
  eigen_map_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    real_type alpha = 1.0/i;
    v3.noalias() = v2 + alpha*v1;
    v2.noalias() = v3;
  }
  tm.toc();
  eigen_fixed_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    Eigen::Map<vecN_t> vv1(nullptr), vv2(nullptr), vv3(nullptr);
    new (&vv1) Eigen::Map<vecN_t>(V1);
    new (&vv2) Eigen::Map<vecN_t>(V2);
    new (&vv3) Eigen::Map<vecN_t>(V3);
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      real_type alpha = 1.0/i;
      vv3.noalias() = vv2 + alpha*vv1;
      vv1.noalias() = vv3;
    }
    tm.toc();
  }
  eigen_fixed_map_time = to_ps*tm.elapsed_ms();

  fmt::print( "(VV) AXPY\n" );

  message();

  msg.green( "All done!\n" );
}

/*\
===============================================================================
===============================================================================
===============================================================================
===============================================================================
\*/

template <int N>
void
testMv() {

  int     N_TIMES = (1000000/N);
  double  to_ps   = 1000000.0/N_TIMES;

  typedef Eigen::Matrix<real_type,N,N> matN_t;
  typedef Eigen::Matrix<real_type,N,1> vecN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<real_type>       base_value("real");
  Malloc<alglin::integer> base_index("integer");

  base_value.allocate(N*N*10);
  base_index.allocate(N*10);

  real_type * M = base_value(N*N);
  real_type * V = base_value(N);
  real_type * R = base_value(N);

  matN_t m;
  dmat_t dm;

  vecN_t v,  r;
  dvec_t dv, dr;

  dm.resize(N,N);
  dv.resize(N);
  dr.resize(N);

  for ( int i{0}; i < N; ++i ) {
    dv(i) = v(i) = V[i] = rand(-1,1);
    dr(i) = r(i) = R[i] = rand(-1,1);
    for ( int j{0}; j < N; ++j ) {
      m(i,j) = dm(i,j) = M[i+j*N] = rand(-1,1);
    }
  }

  Utils::TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    gemv(
      Transposition::NO,
      N, N,
      -1.0, M, N,
      V, 1,
      1.0, R, 1
    );
    alglin::Copy_n( R, N, V );
  }
  tm.toc();
  lapack_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    dr.noalias() -= dm*dv;
    dv.noalias() = dr;
  }
  tm.toc();
  eigen_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    //Eigen::Index Z(0);
    //Eigen::Map<dmat_t> mm(NULL,Z,Z);
    //Eigen::Map<dvec_t> vv(NULL,Z), rr(R,Z);
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      //new (&mm) Eigen::Map<dmat_t>(M,Eigen::Index(N),Eigen::Index(N));
      //new (&vv) Eigen::Map<dvec_t>(V,Eigen::Index(N));
      //new (&rr) Eigen::Map<dvec_t>(R,Eigen::Index(N));
      Eigen::Map<dmat_t> mm(M,static_cast<Eigen::Index>(N),static_cast<Eigen::Index>(N));
      Eigen::Map<dvec_t> vv(V,static_cast<Eigen::Index>(N));
      Eigen::Map<dvec_t> rr(R,static_cast<Eigen::Index>(N));
      rr.noalias() -= mm*vv;
      vv.noalias() = rr;
    }
    tm.toc();
  }
  eigen_map_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    r.noalias() -= m*v;
    v.noalias() = r;
  }
  tm.toc();
  eigen_fixed_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    Eigen::Map<matN_t> mm(nullptr);
    Eigen::Map<vecN_t> vv(nullptr), rr(nullptr);
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      new (&mm) Eigen::Map<matN_t>(M);
      new (&vv) Eigen::Map<vecN_t>(V);
      new (&rr) Eigen::Map<vecN_t>(R);
      rr.noalias() -= mm*vv;
      vv = rr;
    }
    tm.toc();
  }
  eigen_fixed_map_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    alglin::Mv<real_type,N,N,N,N,N>::subTo(M,V,R);
    memcpy( V, R, N*sizeof(real_type) );
    //Vec2<real_type,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MV) MULT\n");
  message();
  fmt::print("{:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  msg.green( "All done!\n" );
}

/*\
===============================================================================
===============================================================================
===============================================================================
===============================================================================
\*/

template <int N>
void
testMM() {

  int     N_TIMES = (10000/N);
  double  to_ps   = 10000.0/N_TIMES;

  typedef Eigen::Matrix<real_type,N,N> matN_t;

  fmt::print("\nSize N = {}\n",N);

  Malloc<real_type>       base_value("real");
  Malloc<alglin::integer> base_index("integer");

  base_value.allocate(N*N*10);
  base_index.allocate(N*10);

  real_type * M1 { base_value( N*N ) };
  real_type * M2 { base_value( N*N ) };
  real_type * M3 { base_value( N*N ) };

  matN_t m1, m2, m3;
  dmat_t dm1, dm2, dm3;

  dm1.resize( static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N) );
  dm2.resize( static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N) );
  dm3.resize( static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N) );

  for ( int i{0}; i < N; ++i ) {
    for ( int j{0}; j < N; ++j ) {
      m1(i,j) = dm1(i,j) = M1[i+j*N] = rand(-1,1);
      m2(i,j) = dm2(i,j) = M2[i+j*N] = rand(-1,1);
      m3(i,j) = dm3(i,j) = M3[i+j*N] = rand(-1,1);
    }
  }

  Utils::TicToc tm;

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    gemm(
      Transposition::NO,
      Transposition::NO,
      N, N, N,
      -1.0, M1, N,
      M2, N,
      1.0, M3, N
    );
    alglin::Copy_n( M3, N*N, M2 );
  }
  tm.toc();
  lapack_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    dm3.noalias() -= dm1*dm2;
    dm2.noalias()  = dm3;
  }
  tm.toc();
  eigen_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    //Eigen::Index Z(0);
    //Eigen::Map<dmat_t> mm1(NULL,Z,Z), mm2(NULL,Z,Z), mm3(NULL,Z,Z);
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      //new (&mm1) Eigen::Map<dmat_t>(M1,Eigen::Index(N),Eigen::Index(N));
      //new (&mm2) Eigen::Map<dmat_t>(M2,Eigen::Index(N),Eigen::Index(N));
      //new (&mm3) Eigen::Map<dmat_t>(M3,Eigen::Index(N),Eigen::Index(N));
      Eigen::Map<dmat_t> mm1(M1,static_cast<Eigen::Index>(N),static_cast<Eigen::Index>(N));
      Eigen::Map<dmat_t> mm2(M2,static_cast<Eigen::Index>(N),static_cast<Eigen::Index>(N));
      Eigen::Map<dmat_t> mm3(M3,static_cast<Eigen::Index>(N),static_cast<Eigen::Index>(N));
      mm3.noalias() -= mm1*mm2;
      mm2.noalias()  = mm3;
    }
    tm.toc();
  }
  eigen_map_dynamic_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    m3.noalias() -= m1*m2;
    m2.noalias() = m3;
  }
  tm.toc();
  eigen_fixed_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  {
    Eigen::Map<matN_t> mm1(nullptr), mm2(nullptr), mm3(nullptr);
    tm.tic();
    for ( int i{0}; i < N_TIMES; ++i ) {
      new (&mm1) Eigen::Map<matN_t>(M1);
      new (&mm2) Eigen::Map<matN_t>(M2);
      new (&mm3) Eigen::Map<matN_t>(M3);
      mm3.noalias() -= mm1*mm2;
      mm2.noalias()  = mm3;
    }
    tm.toc();
  }
  eigen_fixed_map_time = to_ps*tm.elapsed_ms();

  // ===========================================================================

  tm.tic();
  for ( int i{0}; i < N_TIMES; ++i ) {
    MM<real_type,N,N,N,N,N,N>::subTo(M1,M2,M3);
    memcpy( M2, M3, N*N*sizeof(real_type) );
    //Vec2<real_type,N*N,1,1>::copy(M3,M2);
  }
  tm.toc();
  fmt::print("(MM) MULT\n" );

  message();

  fmt::print("(MM) MULT = {:8.4} [ps] (hand unrolled)\n", to_ps*tm.elapsed_ms() );

  // ===========================================================================

  msg.green( "All done!\n" );
}


static
void
testVVall() {
  testVV<2>();
  //testCopy<3>();
  testVV<4>();
  //testCopy<5>();
  testVV<6>();
  //testCopy<7>();
  testVV<8>();
  //testCopy<9>();
  //testCopy<10>();
  //testCopy<11>();
  testVV<12>();
  //testCopy<13>();
  //testCopy<14>();
  //testCopy<15>();
  testVV<16>();
  //testCopy<17>();
  //testCopy<18>();
  //testCopy<19>();
  //testCopy<20>();
  testVV<32>();
  testVV<64>();
  testVV<100>();
  testVV<1000>();
  testVV<10000>();
}

static
void
testMvAll() {
  testMv<2>();
  //testMv<3>();
  testMv<4>();
  //testMv<5>();
  testMv<6>();
  //testMv<7>();
  testMv<8>();
  //testMv<9>();
  //testMv<10>();
  //testMv<11>();
  testMv<12>();
  //testMv<13>();
  //testMv<14>();
  //testMv<15>();
  testMv<16>();
  //testMv<17>();
  //testMv<18>();
  //testMv<19>();
  //testMv<20>();
  testMv<32>();
  testMv<64>();
  testMv<100>();
}

static
void
testMMall() {
  testMM<2>();
  //testMM<3>();
  testMM<4>();
  //testMM<5>();
  testMM<6>();
  //testMM<7>();
  testMM<8>();
  //testMM<9>();
  //testMM<10>();
  //testMM<11>();
  testMM<12>();
  //testMM<13>();
  //testMM<14>();
  //testMM<15>();
  testMM<16>();
  //testMM<17>();
  //testMM<18>();
  //testMM<19>();
  //testMM<20>();
  testMM<100>();
}

static
void
testCopyAll() {
  testCopy<2>();
  //testCopy<3>();
  testCopy<4>();
  //testCopy<5>();
  testCopy<6>();
  //testCopy<7>();
  testCopy<8>();
  //testCopy<9>();
  //testCopy<10>();
  //testCopy<11>();
  testCopy<12>();
  //testCopy<13>();
  //testCopy<14>();
  //testCopy<15>();
  testCopy<16>();
  //testCopy<17>();
  //testCopy<18>();
  //testCopy<19>();
  //testCopy<20>();
  testCopy<100>();
  //testCopy<1000>();
  //testCopy<10000>();
}

int
main() {
  Eigen::initParallel();
  Eigen::setNbThreads(4);

  //testCopyAll();
  //testVVall();
  //testMvAll();
  testMMall();

  fmt::print("\n\nNUM THREAD {}\n",Eigen::nbThreads());
  fmt::print("\n\nAll done!\n");
  return 0;
}
