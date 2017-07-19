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

#include "Alglin.hh"
#include "Alglin++.hh"

#include "ABD_Diaz.hh"
//#include "ABD_Arceco.hh"
//#include "BorderedCR.hh"
//#include "BABD_Block.hh"

#include "BABD_C_interface.h"
#include <map>
#include <string>

namespace alglin {

  static std::map<ABD_intType,DiazLU<ABD_realType> > abd_database ;
  static string abd_last_error = "no error" ;

  extern "C"
  int
  ABD_factorize( ABD_intType        mat_id,
                 ABD_intType        row0,
                 ABD_intType        col0,
                 ABD_realType const TOP[], ABD_intType ldTOP,
                 ABD_intType        nblock,
                 ABD_intType        n,
                 ABD_realType const DE[], ABD_intType ldDE,
                 ABD_intType        rowN,
                 ABD_intType        colN,
                 ABD_realType const BOTTOM[], ABD_intType ldBOTTOM ) {
    try {
      DiazLU<ABD_realType> & lu = abd_database[mat_id] ; // find or create
      lu.allocateTopBottom( nblock, n, row0, col0, rowN, colN, 0 );
      lu.loadTopBottom( TOP, ldTOP, BOTTOM, ldBOTTOM ) ;
      lu.loadBlocks( DE, ldDE ) ;
      lu.factorize();
    }
    catch ( exception const & err ) {
      abd_last_error = err.what() ;
      return -1;
    }
    catch ( ... ) {
      abd_last_error = "ABD_factorize unknown error" ;
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_solve( ABD_intType mat_id, ABD_realType rhs_sol[] ) {
    try {
      std::map<ABD_intType,DiazLU<ABD_realType> >::const_iterator it = abd_database.find(mat_id);
      if ( it == abd_database.end() ) {
        abd_last_error = "ABD_solve mat_id do not correspond to any factorization" ;
        return -3;
      }
      it->second.solve_ABD( rhs_sol ) ;
    }
    catch ( exception const & err ) {
      abd_last_error = err.what() ;
      return -1;
    }
    catch ( ... ) {
      abd_last_error = "ABD_solve unknown error" ;
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_solve_nrhs( ABD_intType  mat_id,
                  ABD_intType  nrhs,
                  ABD_realType rhs_sol[],
                  ABD_intType  ldRhs )  {
    try {
      std::map<ABD_intType,DiazLU<ABD_realType> >::const_iterator it = abd_database.find(mat_id);
      if ( it == abd_database.end() ) {
        abd_last_error = "ABD_solve_nrhs mat_id do not correspond to any factorization" ;
        return -3;
      }
      it->second.solve_ABD( nrhs, rhs_sol, ldRhs ) ;
    }
    catch ( exception const & err ) {
      abd_last_error = err.what() ;
      return -1;
    }
    catch ( ... ) {
      abd_last_error = "ABD_factorize unknown error" ;
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_free( ABD_intType mat_id ) {
    abd_database.erase(mat_id) ;
    return 0 ;
  }

  extern "C"
  char const *
  ABD_get_last_error( )
  { return abd_last_error.c_str() ; }

}
