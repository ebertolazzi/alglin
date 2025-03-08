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
  constexpr alglin::real_type s1[]{1,2,3};
  constexpr alglin::real_type y1[]{11./3.,-6./3.,13./3.};

  constexpr alglin::real_type s2[]{1,0,1};
  constexpr alglin::real_type y2[]{3,-2,3};

  constexpr alglin::real_type s3[]{0,-1,1};
  constexpr alglin::real_type y3[]{7/3.,-6/3.,8/3.};

  alglin::BFGS<alglin::real_type> bfgs;

  constexpr alglin::real_type epsi{1e-8};

  bfgs.allocate(3);
  bfgs.init();
  fmt::print( "\n\n{}", bfgs.to_string() );

  bfgs.update( y1, s1, epsi );
  fmt::print( "\n\n{}", bfgs.to_string() );

  bfgs.update( y2, s2, epsi );
  fmt::print( "\n\n{}", bfgs.to_string() );

  bfgs.update( y3, s3, epsi );
  fmt::print( "\n\n{}", bfgs.to_string() );

  for ( int i{0}; i < 90 ; ++i ) {
    bfgs.update( y1, s1, epsi );
    bfgs.update( y2, s2, epsi );
    bfgs.update( y3, s3, epsi );
  }
  fmt::print( "\n\n{}", bfgs.to_string() );

  msg.green("All done!\n");
  return 0;
}
