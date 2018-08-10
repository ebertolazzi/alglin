/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright 2016                                                          |                                                                          |
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

#include "Simplex.hh"
#include <vector>
#include <iostream>

using namespace std;

int
main() {
  using Simplex::infinity;
  Simplex::integer m = 3;
  Simplex::integer n = 5;
  Simplex::valueType A[] = { 1, 40, 6,
                             1, 120, 12,
                             1, 0, 0,
                             0, 1, 0,
                             0, 0, 1 };
  Simplex::valueType b[]  = { 40, 2400, 312 };
  Simplex::valueType c[]  = { -100, -250, 0, 0, 0 };
  Simplex::valueType L[]  = { 0, 0, 0, 0, 0 };
  Simplex::valueType U[]  = { infinity, infinity, infinity, infinity, infinity };
  Simplex::valueType x[]  = { 0, 0, 40, 2400, 312 };
  Simplex::integer   IB[] = { 2, 3, 4 };

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

    cout << "xdd = "
         << xdd[0] << " "
         << xdd[1] << " "
         << xdd[2] << " "
         << xdd[3] << " "
         << xdd[4] << "\n";
    cout << "IBd = "
         << IBd[0] << " "
         << IBd[1] << " "
         << IBd[2] << "\n";
    simplex.solve( &simplex_problem, xdd, IBd );
  }
  catch (  exception const & err ) {
    cerr << "Error: " << err.what() << "\n";
  }
  catch (...) {
    cerr << "Unknwn error\n";
  }

  cout << "\nAll Done Folks!\n";
  return 0;
}
