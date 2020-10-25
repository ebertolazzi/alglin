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
 |  Licensed under the EUPL, Version 1.1 or – as soon they will be          |
 |  approved by the European Commission - subsequent versions of the EUPL   |
 |  (the "Licence"); You may not use this work except in compliance with    |
 |  the Licence.                                                            |
 |  You may obtain a copy of the Licence at:                                |
 |                                                                          |
 |  http://ec.europa.eu/idabc/eupl5                                         |
 |                                                                          |
 |  Unless required by applicable law or agreed to in writing, software     |
 |  distributed under the Licence is distributed on an "AS IS" basis,       |
 |  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or         |
 |  implied.                                                                |
 |  See the Licence for the specific language governing permissions and     |
 |  limitations under the Licence.                                          |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

int
main() {

  using Simplex::infinity;

  // Cycling pivot test
  Simplex::integer   m    = 3;
  Simplex::integer   n    = 7;
  Simplex::valueType A[]  = {  0.5,  0.5,  1,
                              -5.5, -1.5,  0,
                              -2.5, -0.5,  0,
                                 9,    1,  0,
                                 1,    0,  0,
                                 0,    1,  0,
                                 0,    0,  1};
  Simplex::valueType b[]  = { 0, 0, 1 };
  Simplex::valueType c[]  = { 10, 57, 9, 24, 0, 0, 0 };
  Simplex::valueType L[]  = { 0, 0, 0, 0, 0, 0, 0 };
  Simplex::valueType U[]  = { infinity, infinity, infinity, infinity, infinity, infinity, infinity };
  Simplex::valueType x[]  = { 0, 0, 0, 0, 0, 0, 1 };
  Simplex::integer   IB[] = { 4, 5, 6 };
  
  Simplex::StandardProblem simplex_problem;
  Simplex::AuxProblem      simplex_problem_aux;
  Simplex::StandardSolver  simplex("simplex");
  try {
    simplex_problem.setup( m, n, A, m, b, c, L, U );
    simplex_problem_aux.setup( &simplex_problem );
    simplex.solve( &simplex_problem, x, IB );

    Simplex::valueType xd[100], xdd[100];
    Simplex::integer   IBd[100];
    
    std::cout << "\n\n\n\n\n\n\n\n\n\n";
    
    simplex_problem_aux.feasible_point( xd, IBd );
    simplex.solve( &simplex_problem_aux, xd, IBd );
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
