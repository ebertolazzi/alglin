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
#include "Alglin.hh"
#include "Alglin++.hh"
#include "Alglin_tmpl.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wc++98-compat-pedantic"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wdocumentation-unknown-command"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wdocumentation"
#pragma GCC diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma GCC diagnostic ignored "-Wused-but-marked-unused"
#pragma GCC diagnostic ignored "-Wdeprecated"
#pragma GCC diagnostic ignored "-Wundefined-func-template"
#pragma GCC diagnostic ignored "-Wextra-semi"
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
#endif

#ifndef ALGLIN_USE_CXX11

using namespace std ;

int
main() {
  cout << "To test timing you must compile with a C++11 capable compiler\nAll done!\n" ;
  return 0 ;
}

#else
#include "TicToc.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wc99-extensions"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif

#ifdef ALGLIN_USE_SYSTEM_EIGEN
  #include <Eigen/Dense>
#else
  #include "Eigen/Dense"
#endif

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
static std::mt19937 generator(seed1);

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

using namespace alglin ;
typedef Eigen::Matrix<valueType,Eigen::Dynamic,Eigen::Dynamic> dmat_t ;

#define N_TIMES 1000000

template <int N>
void
testN() {

  typedef Eigen::Matrix<valueType,N,N> matN_t ;

  cout << "\nSize N = " << N << "\n" ;

  Malloc<valueType>       baseValue("real") ;
  Malloc<alglin::integer> baseIndex("integer") ;

  baseValue.allocate(N*N*10) ;
  baseIndex.allocate(N*10) ;

  valueType * M1 = baseValue(N*N) ;
  valueType * M2 = baseValue(N*N) ;
  valueType * M3 = baseValue(N*N) ;
  
  matN_t m1, m2, m3 ;
  dmat_t dm1, dm2, dm3 ;
  
  dm1.resize(N,N) ;
  dm2.resize(N,N) ;
  dm3.resize(N,N) ;
  
  for ( int i = 0 ; i < N ; ++i ) {
    for ( int j = 0 ; j < N ; ++j ) {
      m1(i,j) = dm1(i,j) = M1[i+j*N] = rand(-1,1) ;
      m2(i,j) = dm2(i,j) = M2[i+j*N] = rand(-1,1) ;
      m3(i,j) = dm3(i,j) = M3[i+j*N] = rand(-1,1) ;
    }
  }

  TicToc tm ;
  tm.reset() ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    gemm( NO_TRANSPOSE, NO_TRANSPOSE,
          N, N, N,
          -1.0, M1, N,
          M2, N,
          1.0, M3, N ) ;
    copy( N*N, M3, 1, M2, 1) ;
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (lapack)\n" ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    dm3.noalias() -= dm1*dm2 ;
    dm2 = dm3 ;
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (eigen dynamic)\n" ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    Eigen::Map<dmat_t> mm1(M1,N,N) ;
    Eigen::Map<dmat_t> mm2(M2,N,N) ;
    Eigen::Map<dmat_t> mm3(M3,N,N) ;
    mm3.noalias() -= mm1*mm2 ;
    mm2 = mm3 ;
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (eigen map dynamic)\n" ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    m3.noalias() -= m1*m2 ;
    m2 = m3 ;
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (eigen fixed)\n" ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    Eigen::Map<matN_t> mm1(M1) ;
    Eigen::Map<matN_t> mm2(M2) ;
    Eigen::Map<matN_t> mm3(M3) ;
    mm3.noalias() -= mm1*mm2 ;
    mm2 = mm3 ;
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (eigen fixed map)\n" ;

  // ===========================================================================

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    MM<valueType,N,N,N,N,N,N>::subTo(M1,M2,M3) ;
    memcpy( M2, M3, N*N*sizeof(valueType) ) ;
    //Vec2<valueType,N*N,1,1>::copy(M3,M2);
  }
  tm.toc() ;
  cout << "MULT = " << tm.elapsedMilliseconds() << " [ms] (hand unrolled)\n" ;

  // ===========================================================================

  cout << "All done!\n" ;
}



int
main() {

  testN<2>() ;
  testN<3>() ;
  testN<4>() ;
  testN<5>() ;
  testN<6>() ;
  testN<7>() ;
  testN<8>() ;

  cout << "\n\nAll done!\n" ;

  return 0 ;
}

#endif

