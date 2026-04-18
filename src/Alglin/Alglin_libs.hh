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
/// file: Alglin_libs.hh
///

/*!
 * \file Alglin_libs.hh
 * \brief Auto-linking helper for Microsoft Visual C++ toolchains.
 *
 * This header does not expose numerical APIs; it centralizes the
 * `#pragma comment(lib, ...)` directives used on Windows to automatically link
 * the correct static-library variant according to architecture and build type.
 */

#pragma once

#ifndef ALGLIN_LIBS_HH
#define ALGLIN_LIBS_HH

#include "Alglin_Config.hh"

// automatic library inclusion
#if defined(UTILS_OS_WINDOWS) && defined(_MSC_VER)
  #if defined(_DEBUG) || defined(DEBUG)
    #ifdef ALGLIN_ARCH64
      #pragma comment(lib, "Alglin_win_x64_static_debug.lib")
    #else
      #pragma comment(lib, "Alglin_win_x86_static_debug.lib")
    #endif
  #else
    #ifdef ALGLIN_ARCH64
      #pragma comment(lib, "Alglin_win_x64_static.lib")
    #else
      #pragma comment(lib, "Alglin_win_x86_static.lib")
    #endif
  #endif
#endif

#endif

///
/// eof: Alglin_libs.hh
///
