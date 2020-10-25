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

#include "Alglin/BABD_C_interface.h"
#include <stdio.h>

static
void
test1() {

  int i;

  BABD_intType mat_id          = 1;
	BABD_intType mat_fact        = 0;
	BABD_intType last_block_fact = 0;
	BABD_intType nblock          = 2;
	BABD_intType n               = 2;
	BABD_intType qr              = 1;
	BABD_intType qx              = 1;

  BABD_intType ldDE = 2;
  BABD_realType const DE[] = {
    10,2,
    2,40,
    2,1,
    3,-1,
        10,2,
        2,40,
        2,1,
        3,-1
  };

  BABD_intType ldH0 = 3;
  BABD_realType const H0[] = {
    0,1,1,
    1,0,2
  };

  BABD_intType ldHN = 3;
  BABD_realType const HN[] = {
    20,0,2,
    2,-10,2
  };

  BABD_intType ldHq = 3;
  BABD_realType const Hq[] = {
    1,1,30
  };

  /*
     Matrix
     10   2   2   3   .   .   .
      2  40   1  -1   .   .   .
      .   .  10   2   2   3   .
      .   .   2  40   1  -1   .
      0   1   .   .  20   2   1
      1   0   .   .   0 -10   1
      1   2   .   .   2   2  30
  */

  /* BABD_realType sol[] = { 1, 1, 1, 1, 1, 1, 1 }; */
  /* BABD_realType rhs[] = { 17, 42, 17, 42, 24, -8, 37 }; */

  /* BABD_realType sol[] = { 1, 2, 3, 4, 5, 6, 7 }; */
  BABD_realType rhs[] = { 32, 81, 66, 165, 121, -52, 237 };

  int ok = BABD_factorize( mat_id, mat_fact, last_block_fact, nblock, n, qr, qx,
                           DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq );
  if ( ok != 0 ) printf("ERR = %s\n", BABD_get_last_error() );

  for ( i=0; i < 7; ++i )
    printf("rhs[%d] = %lf\n", i, rhs[i] );

  ok = BABD_solve( mat_id, rhs );
  if ( ok != 0 ) printf("ERR = %s\n", BABD_get_last_error() );
  for ( i=0; i < 7; ++i )
    printf("x[%d] = %lf\n", i, rhs[i] );
}

static
void
test2() {
  int i;
  /*
     Matrix
     10   2   2   3   .   .   .   1
      2  40   1  -1   .   .   .   1
      .   .  10   2   2   3   .   1
      .   .   2  40   1  -1   .   1
      0   1   .   .  20   2   1   1
      1   0   .   .   0 -10   1   1
      1   2   .   .   2   2  30   1
      1  -1   1  -1   1  -1   1  -1
  */

  BABD_intType mat_id          = 1;
	BABD_intType mat_fact        = 0;
	BABD_intType last_block_fact = 0;
	BABD_intType nblock          = 2;
	BABD_intType n               = 2;
	BABD_intType qr              = 1;
	BABD_intType qx              = 1;

  BABD_intType ldDE = 2;
  BABD_realType const DE[] = {
    10,2,
    2,40,
    2,1,
    3,-1,
        10,2,
        2,40,
        2,1,
        3,-1
  };

  BABD_intType ldH0 = 3;
  BABD_realType const H0[] = {
    0,1,1,
    1,0,2
  };

  BABD_intType ldHN = 3;
  BABD_realType const HN[] = {
    20,0,2,
    2,-10,2
  };

  BABD_intType ldHq = 3;
  BABD_realType const Hq[] = {
    1,1,30
  };


	BABD_intType nr = 1;
	BABD_intType nx = 1;

  BABD_intType ldB = 7;
  BABD_realType const B[] = { 1,1,1,1,1,1,1 };

  BABD_intType ldC = 1;
  BABD_realType const C[] = { 1,-1,1,-1,1,-1,1 };

  BABD_intType ldD = 1;
  BABD_realType const D[] = { -1 };

  BABD_realType rhs1[] = { 32, 81, 66, 165, 121, -52, 237, 4 };

  int ok = BABD_factorize_bordered( mat_id, mat_fact, last_block_fact,
                                    nblock, n, qr, nr, qx, nx,
                                    DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq,
                                    B, ldB, C, ldC, D, ldD );
  if ( ok != 0 ) printf("ERR = %s\n", BABD_get_last_error() );

  for ( i=0; i < 8; ++i )
    printf("rhs[%d] = %lf\n", i, rhs1[i] );

  ok = BABD_solve( mat_id, rhs1 );
  if ( ok != 0 ) printf("ERR = %s\n", BABD_get_last_error() );
  for ( i=0; i < 8; ++i )
    printf("x[%d] = %lf\n", i, rhs1[i] );
}

int
main() {

  test1();
  test2();
  printf("All done!\n");
  return 0;
}
