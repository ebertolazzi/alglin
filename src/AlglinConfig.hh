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

#define ALGLIN_USE_THREAD
//#define ALGLIN_DO_NOT_USE_CXX11

#ifdef ALGLIN_USE_THREAD
  #define BABD_AMODIO_USE_THREAD
  #define BABD_AMODIO_N_USE_THREAD
  //#define BABD_QR_USE_THREAD
  #define BABD_QR_N_USE_THREAD
  #define BABD_QR_N_USE_PIVOTING
#endif


#endif

///
/// eof: AlglinConfig.hh
///
