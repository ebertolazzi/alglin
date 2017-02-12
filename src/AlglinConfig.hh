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

///
/// file: AlglinConfig.hh
///

#ifndef ALGLIN_CONFIG_HH
#define ALGLIN_CONFIG_HH

//#define ALGLIN_USE_ACCELERATE 1
//#define ALGLIN_USE_ATLAS 1
//#define ALGLIN_USE_OPENBLAS 1
//#define ALGLIN_USE_LAPACK 1

//#define ALGLIN_USE_THREAD
//#define ALGLIN_DO_NOT_USE_CXX11

// select computer architecture
#if defined(__APPLE__) && defined(__MACH__)
  // osx architecture
  #define ALGLIN_OS_OSX 1
  #if defined(__i386__)
    #define ALGLIN_ARCH32 1
  #elif defined(__x86_64__)
    #define ALGLIN_ARCH64 1
  #endif
#elif defined(__unix__)
  // linux architecture
  #define ALGLIN_OS_LINUX 1
  #if defined(__i386__)
    #define ALGLIN_ARCH32 1
  #elif defined(__x86_64__)
    #define ALGLIN_ARCH64 1
  #endif
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define ALGLIN_OS_WINDOWS 1
  #if defined(_M_X64) || defined(_M_AMD64)
    #define ALGLIN_ARCH64 1
  #else
    #define ALGLIN_ARCH32 1
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else
  #error "unsupported OS!"
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
  #ifndef ALGLIN_DO_NOT_USE_CXX11
    #define ALGLIN_USE_CXX11
  #endif
#else
  // not C++11 compiler
  #ifndef nullptr
    #define nullptr NULL
  #endif
#endif


#ifdef ALGLIN_USE_THREAD
  #ifndef ALGLIN_USE_CXX11
    #error "Alglin libray compiled without c++11 support, cannot use thread"
  #endif
  #define CYCLIC_REDUCTION_USE_THREAD
  #define CYCLIC_REDUCTION_USE_FIXED_SIZE
  #define BORDERED_CYCLIC_REDUCTION_USE_THREAD
  #define BABD_AMODIO_N_USE_THREAD
  #define BABD_QR_USE_THREAD
  #define BABD_QR_N_USE_THREAD
  #define BABD_QR_N_USE_PIVOTING
#endif

#endif

///
/// eof: AlglinConfig.hh
///
