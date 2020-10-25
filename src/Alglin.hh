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
/// file: Alglin.hh
///

#pragma once

#ifndef ALGLIN_dot_HH
#define ALGLIN_dot_HH

#include "Alglin_Config.hh"

#include <algorithm>
#include <vector>
#include <map>

#include <thread>
#include <mutex>

#include <string>
#include <cmath>

#include "Alglin_Eigen.hxx"
#include "Alglin/Alglin_aux.hxx"
#include "Alglin/Alglin_SuperLU.hxx"
#include "Alglin/Alglin_FD.hxx"

#include "Alglin/BlockBidiagonal.hxx"
#include "Alglin/ABD_Block.hxx"
#include "Alglin/ABD_Arceco.hxx"
#include "Alglin/ABD_Diaz.hxx"
#include "Alglin/BABD_Block.hxx"
#include "Alglin/BABD_BorderedCR.hxx"
#include "Alglin/BABD_SuperLU.hxx"
#include "Alglin/KKT_like.hxx"
#include "Alglin/Simplex.hxx"

#include "Alglin/BABD.hxx"

#endif

///
/// eof: Alglin.hh
///
