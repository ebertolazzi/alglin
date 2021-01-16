/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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
#include "Alglin/BABD_C_interface.h"

namespace alglin {

  /*\
   |     _    ____  ____
   |    / \  | __ )|  _ \
   |   / _ \ |  _ \| | | |
   |  / ___ \| |_) | |_| |
   | /_/   \_\____/|____/
   |
  \*/

  typedef std::map<ABD_intType,std::shared_ptr<DiazLU<ABD_realType> > > MAP_DIAZ;
  using Utils::ThreadPool;

  static ThreadPool  * TP             = nullptr;
  static MAP_DIAZ    * abd_database   = new MAP_DIAZ(); // workaround da indagare
  static std::string * abd_last_error = new std::string("no error");

  extern "C"
  int
  ABD_factorize(
    ABD_intType        mat_id,
    ABD_intType        row0,
    ABD_intType        col0,
    ABD_realType const TOP[], ABD_intType ldTOP,
    ABD_intType        nblock,
    ABD_intType        n,
    ABD_realType const DE[], ABD_intType ldDE,
    ABD_intType        rowN,
    ABD_intType        colN,
    ABD_realType const BOTTOM[], ABD_intType ldBOTTOM
  ) {
    try {
      MAP_DIAZ::iterator lu = abd_database->find(mat_id); // find or create
      if ( lu == abd_database->end() ) {
        MAP_DIAZ::mapped_type ptr( new DiazLU<ABD_realType>() );
        MAP_DIAZ::value_type  P( mat_id, ptr );
        pair<MAP_DIAZ::iterator,bool> res = abd_database->insert( P );
        lu = res.first;
      }

      lu->second->allocateTopBottom( nblock, n, row0, col0, rowN, colN, 0 );
      lu->second->loadTopBottom( TOP, ldTOP, BOTTOM, ldBOTTOM );
      lu->second->loadBlocks( DE, ldDE );
      lu->second->factorize();
    }
    catch ( std::exception const & err ) {
      (*abd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*abd_last_error) = "ABD_factorize unknown error";
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_solve( ABD_intType mat_id, ABD_realType rhs_sol[] ) {
    try {
      MAP_DIAZ::const_iterator it = abd_database->find(mat_id);
      if ( it == abd_database->end() ) {
        (*abd_last_error) = "ABD_solve mat_id do not correspond to any factorization";
        return -3;
      }
      it->second->solve_ABD( rhs_sol );
    }
    catch ( std::exception const & err ) {
      (*abd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*abd_last_error) = "ABD_solve unknown error";
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_solve_nrhs(
    ABD_intType  mat_id,
    ABD_intType  nrhs,
    ABD_realType rhs_sol[],
    ABD_intType  ldRhs
  )  {
    try {
      MAP_DIAZ::const_iterator it = abd_database->find(mat_id);
      if ( it == abd_database->end() ) {
        (*abd_last_error) = "ABD_solve_nrhs mat_id do not correspond to any factorization";
        return -3;
      }
      it->second->solve_ABD( nrhs, rhs_sol, ldRhs );
    }
    catch ( std::exception const & err ) {
      (*abd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*abd_last_error) = "ABD_factorize unknown error";
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  ABD_free( ABD_intType mat_id ) {
    abd_database->erase(mat_id);
    return 0;
  }

  extern "C"
  char const *
  ABD_get_last_error( )
  { return abd_last_error->c_str(); }

  extern "C"
  void
  ABD_get_last_error_f90( char res[], long len ) {
    #ifdef UTILS_OS_WINDOWS
    errno_t e = strncpy_s( res, len, abd_last_error->c_str(), abd_last_error->length() );
    #else
    strncpy( res, abd_last_error->c_str(), size_t(len) );
    #endif
  }

  /*\
   |  ____    _    ____  ____
   | | __ )  / \  | __ )|  _ \
   | |  _ \ / _ \ |  _ \| | | |
   | | |_) / ___ \| |_) | |_| |
   | |____/_/   \_\____/|____/
   |
  \*/

  typedef std::map<BABD_intType,std::shared_ptr<BorderedCR<BABD_realType> > > MAP_BABD;

  static MAP_BABD    * babd_database   = new MAP_BABD();
  static std::string * babd_last_error = new std::string("no error");

  extern "C"
  int
  BABD_factorize(
    BABD_intType        mat_id,
    BABD_intType        mat_fact,
    BABD_intType        last_block_fact,
    BABD_intType        nblock,
    BABD_intType        n,
    BABD_intType        qr,
    BABD_intType        qx,
    BABD_realType const DE[], BABD_intType ldDE,
    BABD_realType const H0[], BABD_intType ldH0,
    BABD_realType const HN[], BABD_intType ldHN,
    BABD_realType const Hq[], BABD_intType ldHq
  ) {
    try {

      MAP_BABD::iterator lu = babd_database->find(mat_id); // find or create
      if ( lu == babd_database->end() ) {
        if ( TP == nullptr ) TP = new ThreadPool( std::thread::hardware_concurrency() );
        MAP_BABD::mapped_type ptr( new BorderedCR<BABD_realType>( TP ) );
        MAP_BABD::value_type  P( mat_id, ptr );
        pair<MAP_BABD::iterator,bool> res = babd_database->insert( P );
        lu = res.first;
      }

      lu->second->allocate( nblock, n, qr, qx, 0, 0 );
      switch ( mat_fact ) {
        case 0: lu->second->select_LU();  break;
        case 1: lu->second->select_QR();  break;
        case 2: lu->second->select_QRP(); break;
      }
      switch ( last_block_fact ) {
        case 0: lu->second->select_last_LU();   break;
        case 1: lu->second->select_last_LUPQ(); break;
        case 2: lu->second->select_last_QR();   break;
        case 3: lu->second->select_last_QRP();  break;
        case 4: lu->second->select_last_SVD();  break;
        case 5: lu->second->select_last_LSS();  break;
        case 6: lu->second->select_last_LSY();  break;
      }
      lu->second->loadBottom( H0, ldH0, HN, ldHN, Hq, ldHq, nullptr, 0 );
      for ( BABD_intType nbl = 0; nbl < nblock; ++nbl )
        lu->second->loadDE( nbl, DE + 2*nbl*n*ldDE, ldDE );
      lu->second->factorize();
    }
    catch ( std::exception const & err ) {
      (*babd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*babd_last_error) = "BABD_factorize unknown error";
      return -2;
    }
    return 0;

  }

  extern "C"
  int
  BABD_factorize_bordered(
    BABD_intType        mat_id,
    BABD_intType        mat_fact,
    BABD_intType        last_block_fact,
    BABD_intType        nblock,
    BABD_intType        n,
    BABD_intType        qr,
    BABD_intType        nr,
    BABD_intType        qx,
    BABD_intType        nx,
    BABD_realType const DE[], BABD_intType ldDE,  // n x (2*n*nblock)
    BABD_realType const H0[], BABD_intType ldH0,  // (n+qr) x n
    BABD_realType const HN[], BABD_intType ldHN,  // (n+qr) x n
    BABD_realType const Hq[], BABD_intType ldHq,  // (n+qr) x qx
    BABD_realType const B[],  BABD_intType ldB,   // (n*(nblock+1)+qr) x nx
    BABD_realType const C[],  BABD_intType ldC,   // nr x (n*(nblock+1)+qx)
    BABD_realType const D[],  BABD_intType ldD    // nr x nx
  ) {
    try {

      MAP_BABD::iterator lu = babd_database->find(mat_id); // find or create
      if ( lu == babd_database->end() ) {
        if ( TP == nullptr ) TP = new ThreadPool( std::thread::hardware_concurrency() );
        MAP_BABD::mapped_type ptr( new BorderedCR<BABD_realType>( TP ) );
        MAP_BABD::value_type  P( mat_id, ptr );
        pair<MAP_BABD::iterator,bool> res = babd_database->insert( P );
        lu = res.first;
      }
      lu->second->allocate( nblock, n, qr, qx, nr, nx );
      lu->second->loadBottom( H0, ldH0, HN, ldHN, Hq, ldHq, B+(nblock*n), ldB );
      switch ( mat_fact ) {
        case 0: lu->second->select_LU(); break;
        case 1: lu->second->select_QR(); break;
        case 2: lu->second->select_QRP(); break;
      }
      switch ( last_block_fact ) {
        case 0: lu->second->select_last_LU();   break;
        case 1: lu->second->select_last_LUPQ(); break;
        case 2: lu->second->select_last_QR();   break;
        case 3: lu->second->select_last_QRP();  break;
        case 4: lu->second->select_last_SVD();  break;
        case 5: lu->second->select_last_LSS();  break;
        case 6: lu->second->select_last_LSY();  break;
      }
      for ( BABD_intType nbl = 0; nbl < nblock; ++nbl ) {
        lu->second->loadB( nbl, B + nbl*n, ldB );
        lu->second->loadC( nbl, C + nbl*n*ldC, ldC );
        lu->second->loadDE( nbl, DE + 2*nbl*n*ldDE, ldDE );
      }
      lu->second->loadC( nblock, C + nblock*n*ldC, ldC );
      lu->second->loadCq( C + (nblock+1)*n*ldC, ldC );
      lu->second->loadF( D, ldD );
      lu->second->factorize();
    }
    catch ( std::exception const & err ) {
      (*abd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*abd_last_error) = "BABD_factorize unknown error";
      return -2;
    }
    return 0;

  }

  extern "C"
  int
  BABD_solve( BABD_intType mat_id, BABD_realType rhs_sol[] ) {
    try {
      MAP_BABD::const_iterator it = babd_database->find(mat_id);
      if ( it == babd_database->end() ) {
        (*babd_last_error) = "BABD_solve mat_id do not correspond to any factorization";
        return -3;
      }
      it->second->solve( rhs_sol );
    }
    catch ( std::exception const & err ) {
      (*babd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*babd_last_error) = "BABD_solve unknown error";
      return -2;
    }
    return 0;
  }

  extern "C"
  int
  BABD_solve_nrhs(
    BABD_intType  mat_id,
    BABD_intType  nrhs,
    BABD_realType rhs_sol[],
    BABD_intType  ldRhs
  ) {
    try {
      MAP_BABD::const_iterator it = babd_database->find(mat_id);
      if ( it == babd_database->end() ) {
        (*abd_last_error) = "BABD_solve_nrhs mat_id do not correspond to any factorization";
        return -3;
      }
      it->second->solve( nrhs, rhs_sol, ldRhs );
    }
    catch ( std::exception const & err ) {
      (*babd_last_error) = err.what();
      return -1;
    }
    catch ( ... ) {
      (*babd_last_error) = "BABD_factorize unknown error";
      return -2;
    }
    return 0;

  }

  extern "C"
  int
  BABD_free( BABD_intType mat_id ) {
    babd_database->erase(mat_id);
    return 0;
  }

  extern "C"
  char const *
  BABD_get_last_error( )
  { return babd_last_error->c_str(); }

  extern "C"
  void
  BABD_get_last_error_f90( char res[], long len ) {
    #ifdef UTILS_OS_WINDOWS
    errno_t e = strncpy_s( res, len, babd_last_error->c_str(), babd_last_error->length() );
    #else
    strncpy( res, babd_last_error->c_str(), size_t(len) );
    #endif
  }

}
