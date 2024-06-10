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

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

using namespace std;

int
main() {

  alglin::BandedLU<alglin::real_type> BLU;

  alglin::integer const N  = 10;
  alglin::real_type D[]   = {1,2,3,4,5,6,7,8,9,10};
  alglin::real_type L0[]  = { 1,   1, 1,  1, 1,  1, 1, 1, 1};
  alglin::real_type L1[]  = { 1,  -1, 1, -1, 1, -1, 1, 1};
  alglin::real_type L2[]  = { 1,-1.4, 1, -1, 1, -1, 1 };
  alglin::real_type U0[]  = {-1,  -1,-1, -1,-1, -1,-1, -1,-1};
  alglin::real_type U1[]  = { 1,  -1, 1, -1, 1, -1, 1, -1};
  alglin::real_type rhs[N];

  BLU.setup( N, N, 3, 2 );
  BLU.zero();

  for ( int i = 0; i < N; ++i )
    rhs[i] = BLU(i,i) = D[i];

  for ( int i = 0; i < N-1; ++i ) {
    rhs[i]   += U0[i]; BLU(i,i+1) = U0[i];
    rhs[i+1] += L0[i]; BLU(i+1,i) = L0[i];
  }

  for ( int i = 0; i < N-2; ++i ) {
    rhs[i]   += U1[i]; BLU(i,i+2) = U1[i];
    rhs[i+2] += L1[i]; BLU(i+2,i) = L1[i];
  }

  for ( int i = 0; i < N-3; ++i ) {
    rhs[i+3] += L2[i]; BLU(i+3,i) = L2[i];
  }

  BLU.factorize("BLU");
  BLU.solve(rhs);

  for ( int i = 0; i < N; ++i )
    fmt::print("x[{}] = {}\n",i,rhs[i]);

  msg.green("All done!\n");
  return 0;
}
