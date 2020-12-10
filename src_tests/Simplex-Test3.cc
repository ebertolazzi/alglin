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

  // Cycling pivot test
  Simplex::integer   m    = 3;
  Simplex::integer   n    = 6;
  Simplex::valueType A[]  = {  1,  2,  0,
                               1,  1, -1,
                               1, -1,  1,
                               1,  0,  0,
                               0, -1,  0,
                               0,  0, -1};
  Simplex::valueType b[]  = { 40, 10, 10 };
  Simplex::valueType c[]  = { -2, -3, -1, 0, 0, 0 };
  Simplex::valueType L[]  = { 0, 0, 0, 0, 0, 0 };
  Simplex::valueType U[]  = { infinity, infinity, infinity, infinity, infinity, infinity };
  //Simplex::valueType x[]  = { 0, 0, 0, 0, 0, 1 };
  //Simplex::integer   IB[] = { 0, 1, 2 };
  
  Simplex::StandardProblem simplex_problem;
  Simplex::AuxProblem      simplex_problem_aux;
  Simplex::StandardSolver  simplex("simplex");
  
  try {
    //simplex.solve( &simplex_problem, x, IB );
    Simplex::valueType xd[100], xdd[100];
    Simplex::integer   IBd[100];
    
    simplex_problem.setup( m, n, A, m, b, c, L, U );
    simplex_problem_aux.setup( &simplex_problem );
    simplex_problem_aux.feasible_point( xd, IBd );
    simplex.solve( &simplex_problem_aux, xd, IBd );

    std::cout << "\n\n\n\n\n\n\n\n\n\n";
    
    simplex_problem_aux.to_primal( xd, xdd, IBd );
    simplex.solve( &simplex_problem, xdd, IBd );
  }
  catch ( std::exception const & err ) {
    std::cerr << "Error: " << err.what();
  }
  catch (...) {
    std::cerr << "Unknwn error\n";
  }

  std::cout << "\nAll Done Folks!\n";
  return 0;
}
