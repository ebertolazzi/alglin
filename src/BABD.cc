/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

#include "BABD.hh"

namespace alglin {

  string
  BABD_Choice_to_string( BABD_Choice c ) {
    string res = "none" ;
    switch ( c ) {
    case BABD_DIAZ:                 res = "Diaz" ;                break ;
    case BABD_AMODIO:               res = "Amodio" ;              break ;
    case BABD_CYCLIC_REDUCTION_QR:  res = "CyclicReduction+QR" ;  break ;
    case BABD_CYCLIC_REDUCTION_QRP: res = "CyclicReduction+QRP" ; break ;
    }
    return res ;
  }

}
