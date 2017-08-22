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

#include "BABD_C_interface.h"
#include <stdio.h>

int
main() {
  int i ;
  ABD_intType n      = 2 ;
  ABD_intType nblock = 3 ;
  ABD_intType row0   = 2 ;
  ABD_intType col0   = 3 ;
  ABD_intType rowN   = 3 ;
  ABD_intType colN   = 4 ;

  ABD_intType ldTOP = 2;
  ABD_realType const TOP[] = {
    1, 2,
    2, 3,
    3, -1
  };

  ABD_intType ldBOTTOM = 3;
  ABD_realType const BOTTOM[] = {
    1,  2,  1,
    2,  3,  0,
    3, -1, -1,
    1, -1,  2
  };

  ABD_intType ldDE = 2;
  ABD_realType const DE[] = {
     1, 2,
     3, 4,
     3, -1,
     2, -1,

     1, 2,
     3, 4,
     3, -1,
     2, -1,

     1, 2,
     3, 4,
     3, -1,
     2, -1
  };

  ABD_realType rhs[] = { 6, 4, 9, 4, 9, 4, 9, 4, 7, 3, 2 } ;

  ABD_intType mat_id = 1;

  int ok = ABD_factorize( mat_id,
                          row0, col0, TOP, ldTOP,
                          nblock, n, DE, ldDE,
                          rowN, colN, BOTTOM, ldBOTTOM ) ;
  if ( ok != 0 ) printf("ERR = %s\n", ABD_get_last_error() ) ;

  ok = ABD_solve( mat_id, rhs ) ;
  if ( ok != 0 ) printf("ERR = %s\n", ABD_get_last_error() ) ;
  for ( i=0 ; i < 11 ; ++i )
    printf("x[%d] = %lf\n", i, rhs[i] ) ;

  printf("All done!\n") ;

  return 0 ;
}
