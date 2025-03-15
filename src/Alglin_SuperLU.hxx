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

///
/// file: Alglin_SuperLU.hh
///

// Eigen3
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wpadded"
#endif
#ifdef _MSC_VER
#pragma warning( disable : 4200 )
#endif

// workaround for SUPERLU INCLUSION
#define dgemm_  dgemm_BUGGED
#define dtrsv_  dtrsv_BUGGED
#define dtrsm_  dtrsm_BUGGED
#define dgemv_  dgemv_BUGGED
#define sgemm_  sgemm_BUGGED
#define strsv_  strsv_BUGGED
#define strsm_  strsm_BUGGED
#define sgemv_  sgemv_BUGGED
#define xerbla_ xerbla_BUGGED

#ifdef ALGLIN_USE_SYSTEM_SUPERLU
  #include <superlu/slu_sdefs.h>
  #include <superlu/slu_ddefs.h>
#else
  #include "superlu/slu_sdefs.h"
  #include "superlu/slu_ddefs.h"
#endif

#undef dgemm_
#undef dtrsv_
#undef dtrsm_
#undef dgemv_
#undef sgemm_
#undef strsv_
#undef strsm_
#undef sgemv_
#undef xerbla_

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

///
/// eof: Alglin_SuperLU.hxx
///
