/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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
#include "TicToc.hh"

#ifdef USE_MECHATRONIX_EIGEN
  #include <MechatronixCore/Eigen/Dense>
#else
  #include <Eigen/Dense>
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

int
main() {

  #define N 2
  #define N_TIMES 100000

  Malloc<valueType>       baseValue("real") ;
  Malloc<alglin::integer> baseIndex("integer") ;

  baseValue.allocate(N*N*10) ;
  baseIndex.allocate(N*10) ;

  valueType * M1 = baseValue(N*N) ;
  valueType * M2 = baseValue(N*N) ;
  valueType * M3 = baseValue(N*N) ;
  
  //typedef Eigen::Matrix<valueType,N,1> vec_t ;
  typedef Eigen::Matrix<valueType,N,N> mat_t ;
  //typedef Eigen::Matrix<valueType,Eigen::Dynamic,1>              dvec_t ;
  typedef Eigen::Matrix<valueType,Eigen::Dynamic,Eigen::Dynamic> dmat_t ;
  
  mat_t m1, m2, m3 ;
  dmat_t dm1, dm2, dm3 ;
  
  dm1.resize(N,N) ;
  dm2.resize(N,N) ;
  dm3.resize(N,N) ;
  
  for ( int i = 0 ; i < N ; ++i ) {
    for ( int j = 0 ; j < N ; ++j ) {
      dm1(i,j) = m1(i,j) = M1[i+j*N] = rand(-1,1) ;
      dm2(i,j) = m2(i,j) = M2[i+j*N] = rand(-1,1) ;
      dm3(i,j) = m3(i,j) = M3[i+j*N] = rand(-1,1) ;
    }
  }

  TicToc tm ;
  tm.reset() ;

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
  cout << "MULT (lapack) = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    m3.noalias() -= m1*m2 ;
    m2 = m3 ;
  }
  tm.toc() ;
  cout << "MULT (eigen) = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    dm3.noalias() -= dm1*dm2 ;
    dm2 = dm3 ;
  }
  tm.toc() ;
  cout << "MULT (deigen) = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    Eigen::Map<mat_t> mm1(M1) ;
    Eigen::Map<mat_t> mm2(M2) ;
    Eigen::Map<mat_t> mm3(M3) ;
    mm3.noalias() -= mm1*mm2 ;
    mm2 = mm3 ;
  }
  tm.toc() ;
  cout << "MULT (eigen map) = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    Eigen::Map<dmat_t> mm1(M1,N,N) ;
    Eigen::Map<dmat_t> mm2(M2,N,N) ;
    Eigen::Map<dmat_t> mm3(M3,N,N) ;
    mm3.noalias() -= mm1*mm2 ;
    mm2 = mm3 ;
  }
  tm.toc() ;
  cout << "MULT (deigen map) = " << tm.elapsedMilliseconds() << " [ms]\n" ;


  tm.tic() ;
  for ( int i = 0 ; i < N_TIMES ; ++i ) {
    M3[0+0*N] -= M1[0+0*N]*M2[0+0*N]+M1[0+1*N]*M2[1+0*N] ;
    M3[0+1*N] -= M1[0+0*N]*M2[0+1*N]+M1[0+1*N]*M2[1+1*N] ;
    M3[1+0*N] -= M1[1+0*N]*M2[0+0*N]+M1[1+1*N]*M2[1+0*N] ;
    M3[1+1*N] -= M1[1+0*N]*M2[0+1*N]+M1[1+1*N]*M2[1+1*N] ;
    M2[0+0*N] = M2[0+0*N];
    M2[0+1*N] = M2[0+1*N];
    M2[1+0*N] = M2[1+0*N];
    M2[1+1*N] = M2[1+1*N];
  }
  tm.toc() ;
  cout << "MULT (hand) = " << tm.elapsedMilliseconds() << " [ms]\n" ;

  cout << "All done!\n" ;

  return 0 ;
}
