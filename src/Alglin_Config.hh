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

#define ALGLIN_PURE_VIRTUAL = 0
#ifndef LAPACK_WRAPPER_OS_WINDOWS
  #define ALGLIN_OVERRIDE  override
  #define ALGLIN_CONSTEXPR constexpr
  #ifdef __clang__
    #pragma clang diagnostic ignored "-Wc++98-compat"
  #endif
#else
  #define ALGLIN_OVERRIDE
  #define ALGLIN_CONSTEXPR
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


#ifdef LAPACK_WRAPPER_OS_WINDOWS
  #define ALGLIN_OS_WINDOWS 1
#endif
#ifdef LAPACK_WRAPPER_OS_OSX
  #define ALGLIN_OS_OSX 1
#endif
#ifdef LAPACK_WRAPPER_OS_LINUX
  #define ALGLIN_OS_LINUX 1
#endif
#ifdef LAPACK_WRAPPER_ARCH32
  #define ALGLIN_ARCH32 1
#endif
#ifdef LAPACK_WRAPPER_ARCH64
  #define ALGLIN_ARCH64 1
#endif

namespace alglin {
  using namespace lapack_wrapper;
}

#endif

///
/// eof: Alglin_Config.hh
///
