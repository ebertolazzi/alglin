/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright 2016                                                          |
 |                                                                          |
 |  Enrico Bertolazzi^(*)  and  Matthias Gerdts^(**) (Ingenieurmathematik)  |
 |                                                                          |
 |  (*) Department of Industrial Engineering                                |
 |      University of Trento                                                |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 | (**) Institut fuer Mathematik und Rechneranwendung                       |
 |      Fakultaet fuer Luftund Raumfahrttechnik                             |
 |      Universitaet der Bundeswehr Muenchen                                |
 |      email: matthias.gerdts@unibw.de                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

Utils::Console msg( &std::cout,4);

int
main() {

  using Simplex::infinity;

  // Cycling pivot test
  Simplex::integer   m    = 3;
  Simplex::integer   n    = 3;
  Simplex::real_type A[]  = {  1,  2,  0,
                               1,  1, -1,
                               1, -1,  1 };
  Simplex::real_type c[]  = { -2, -3, -1 };
  Simplex::real_type L[]  = { 0, 0, 0, -infinity, 10, 10 };
  Simplex::real_type U[]  = { infinity, infinity, infinity, 40, infinity, infinity };
  //Simplex::real_type x[]  = { 0, 0, 0, 0, 0, 1 };
  //Simplex::integer   IB[] = { 0, 1, 2 };

  Simplex::Problem        simplex_problem;
  Simplex::AuxProblem     simplex_problem_aux;
  Simplex::StandardSolver simplex("simplex");

  try {
    //simplex.solve( &simplex_problem, x, IB );
    Simplex::real_type xd[100], xdd[100];
    Simplex::integer   IBd[100];

    simplex_problem.setup( m, n, A, m, c, L, U );

    Simplex::StandardProblemAdaptor simplex_problem_adaptor(simplex_problem);

    fmt::print( "simplex_problem\n" );
    simplex_problem.info(std::cout);
    fmt::print( "\n\nsimplex_problem_adaptor\n" );
    simplex_problem_adaptor.info(std::cout);

    simplex_problem_aux.setup( &simplex_problem_adaptor );

    simplex_problem_aux.feasible_point( xd, IBd );
    simplex.solve( &simplex_problem_aux, xd, IBd );

    fmt::print( "\n\n\n\n\n\n\n\n\n\n" );

    simplex_problem_aux.to_primal( xd, xdd, IBd );
    simplex.solve( &simplex_problem_adaptor, xdd, IBd );
  }
  catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
  }
  catch (...) {
    msg.red( "Unknwn error\n" );
  }

  msg.green( "\nAll Done Folks!\n" );
  return 0;
}
