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

///
/// file: Alglin_Config.hh
///

#ifndef ALGLIN_CONFIG_HH
#define ALGLIN_CONFIG_HH

#ifdef ALGLIN_USE_SYSTEM_LAPACK_WRAPPER
  #include <lapack_wrapper/lapack_wrapper.hh>
  #include <lapack_wrapper/lapack_wrapper++.hh>
#else
  #include "lapack_wrapper/lapack_wrapper.hh"
  #include "lapack_wrapper/lapack_wrapper++.hh"
#endif

#ifdef LAPACK_WRAPPER_USE_THREAD
  #ifndef LAPACK_WRAPPER_USE_CXX11
    #error "Alglin libray compiled with c++11 support, cannot use c++ < c++11"
  #endif
  #define CYCLIC_REDUCTION_USE_THREAD
  #define CYCLIC_REDUCTION_USE_FIXED_SIZE
  #define BORDERED_CYCLIC_REDUCTION_USE_THREAD
  #define BABD_AMODIO_N_USE_THREAD
  #define BABD_QR_USE_THREAD
  #define BABD_QR_N_USE_THREAD
  #define BABD_QR_N_USE_PIVOTING
#endif

#define ALGLIN_PURE_VIRTUAL = 0
#if defined(LAPACK_WRAPPER_USE_CXX11) && !defined(LAPACK_WRAPPER_OS_WINDOWS)
  #define ALGLIN_OVERRIDE  override
  #define ALGLIN_CONSTEXPR constexpr
  #ifdef __clang__
    #pragma clang diagnostic ignored "-Wc++98-compat"
  #endif
#else
  #define ALGLIN_OVERRIDE
  #define ALGLIN_CONSTEXPR
#endif

#ifndef ALGLIN_ERROR
  #define ALGLIN_ERROR(MSG) {                    \
    std::ostringstream ost;                      \
    ost << "in file: " << __FILE__ << "\nline: " \
        << __LINE__ << '\n' << MSG << '\n';      \
    throw std::runtime_error(ost.str());         \
  }
#endif

#ifndef ALGLIN_ASSERT
  #define ALGLIN_ASSERT(COND,MSG) \
    if ( !(COND) ) ALGLIN_ERROR( "in alglin::" << MSG )
#endif

#ifdef LAPACK_WRAPPER_USE_CXX11
  #define ALGLIN_USE_CXX11 1
#endif

#ifdef LAPACK_WRAPPER_USE_ACCELERATE
  #define ALGLIN_USE_ACCELERATE 1
#endif
#ifdef LAPACK_WRAPPER_USE_ATLAS
  #define ALGLIN_USE_ATLAS 1
#endif
#ifdef LAPACK_WRAPPER_USE_OPENBLAS
  #define ALGLIN_USE_OPENBLAS 1
#endif
#ifdef LAPACK_WRAPPER_USE_LAPACK
  #define ALGLIN_USE_LAPACK 1
#endif
#ifdef LAPACK_WRAPPER_USE_MKL
  #define ALGLIN_USE_MKL 1
#endif

namespace alglin {
  using namespace lapack_wrapper;
}



#endif

///
/// eof: Alglin_Config.hh
///