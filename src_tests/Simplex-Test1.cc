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

int
main() {

  using Simplex::infinity;
  Simplex::integer m = 3;
  Simplex::integer n = 5;
  Simplex::real_type A[] = { 1, 40, 6,
                             1, 120, 12,
                             1, 0, 0,
                             0, 1, 0,
                             0, 0, 1 };
  Simplex::real_type b[]  = { 40, 2400, 312 };
  Simplex::real_type c[]  = { -100, -250, 0, 0, 0 };
  Simplex::real_type L[]  = { 0, 0, 0, 0, 0 };
  Simplex::real_type U[]  = { infinity, infinity, infinity, infinity, infinity };
  Simplex::real_type x[]  = { 0, 0, 40, 2400, 312 };
  Simplex::integer   IB[] = { 2, 3, 4 };

  Simplex::StandardProblem simplex_problem;
  Simplex::AuxProblem      simplex_problem_aux;
  Simplex::StandardSolver  simplex("simplex");

  try {
    simplex_problem.setup( m, n, A, m, b, c, L, U );
    simplex_problem_aux.setup( &simplex_problem );
    simplex.solve( &simplex_problem, x, IB );

    Simplex::real_type xd[100], xdd[100];
    Simplex::integer   IBd[100];

    fmt::print( "\n\n\n\n\n\n\n\n\n\n" );

    simplex_problem_aux.feasible_point( xd, IBd );
    simplex.solve( &simplex_problem_aux, xd, IBd );
    simplex_problem_aux.to_primal( xd, xdd, IBd );

    fmt::print(
      "xdd = {} {} {} {} {}\n",
      xdd[0], xdd[1], xdd[2], xdd[3], xdd[4]
    );
    fmt::print( "IBd = {} {} {}\n", IBd[0], IBd[1], IBd[2]);
    simplex.solve( &simplex_problem, xdd, IBd );
  }
  catch ( std::exception const & err ) {
    fmt::print( "Error: {}\n", err.what() );
  }
  catch (...) {
     fmt::print( "Unknwn error\n" );
  }

  fmt::print( "\nAll Done Folks!\n" );
  return 0;
}
