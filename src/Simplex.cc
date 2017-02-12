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
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

//! namespace for nonlinear systems and nonlinearsolver
namespace Simplex {

  valueType const epsilon        = std::numeric_limits<valueType>::epsilon() ;
  valueType const relaxedEpsilon = pow(epsilon,0.8) ;
  valueType const infinity       = std::numeric_limits<valueType>::max() ;

 /*\
  |     _             ___         _    _
  |    /_\ _  ___ __ | _ \_ _ ___| |__| |___ _ __
  |   / _ \ || \ \ / |  _/ '_/ _ \ '_ \ / -_) '  \
  |  /_/ \_\_,_/_\_\ |_| |_| \___/_.__/_\___|_|_|_|
 \*/
  void
  AuxProblem::setup( StandardProblemBase * _pBase ) {

    pBase = _pBase ;

    n = pBase->dim_x() ;
    m = pBase->dim_g() ;

    // compute size of variable z, w, and p
    nz = nw = np = 0 ;
    for ( integer i = 0 ; i < n ; ++i ) {
      integer icase = 0 ;
      if ( pBase->Upper_is_free(i) ) icase  = 1 ;
      if ( pBase->Lower_is_free(i) ) icase += 2 ;
      switch ( icase ) {
      case 0: ++nz ; ++nw ; break ; // l <= x <= u
      case 1: ++nz ;        break ; // l <= x
      case 2: ++nw ;        break ; // x <= u
      case 3: ++np ;        break ;
      }
    }
    // allocate memory
    baseInteger.allocate( nz+nw+np+n+m ) ;
    map_z    = baseInteger( nz ) ;
    map_w    = baseInteger( nw ) ;
    map_p    = baseInteger( np ) ;
    map_case = baseInteger( n  ) ;
    i_row    = baseInteger( m  ) ;

    baseReals.allocate( 2*m ) ;
    d      = baseReals( m ) ;
    values = baseReals( m ) ;

    // fill mapping and vector q
    nz = nw = np = 0 ;
    pBase->load_b( d ) ;
    for ( integer i = 0 ; i < n ;  ++i ) {
      valueType q = 0 ;
      integer icase = 0 ;
      if ( pBase->Upper_is_free(i) ) icase  = 1 ;
      if ( pBase->Lower_is_free(i) ) icase += 2 ;
      map_case[i] = icase ;
      switch ( icase ) {
      case 0: // l <= x <= u
        map_z[nz++] = i ;
        map_w[nw++] = i ;
        q = (pBase->Upper(i) + pBase->Lower(i))/2 ;
        break ;
      case 1: // l <= x
        map_z[nz++] = i ;
        q = pBase->Lower(i) ;
        break ;
      case 2: // x <= u
        map_w[nw++] = i ;
        q = pBase->Upper(i) ;
        break ;
      case 3:
        map_p[np++] = i ;
        break ;
      }
      integer nnz = pBase->load_A_column( i, values, i_row ) ;
      for ( integer k = 0 ; k < nnz ; ++k )
        d[i_row[k]] -= q*values[k] ;
    }

    b_max_abs = alglin::absmax( m, d, 1) ;
    c_max_abs = 1 ;
    A_max_abs = pBase->get_A_max_abs() ;

  }

  bool
  AuxProblem::Lower_is_free( integer i ) const
  { return i >= nz+nw && i < nz+nw+np ; }

  bool
  AuxProblem::Upper_is_free( integer   ) const
  { return true ; }

  valueType
  AuxProblem::Lower( integer i ) const {
    if ( Lower_is_free(i) ) return infinity ;
    else                    return 0 ;
  }
  
  valueType
  AuxProblem::Upper( integer ) const
  { return infinity ; }

  void
  AuxProblem::feasible_point( valueType x[], integer IB[] ) const {
    integer nn = nz+nw+np ;
    alglin::zero( nn, x, 1 ) ;
    for ( integer i = 0 ; i < m ; ++i ) {
      x[nn+i] = std::abs(d[i]) ;
      IB[i] = nn+i ;
    }
  }

  /*\
   !           +
   !           | (z_i-w_i+(l_i+u_i))/2 if -infinity < l_i and u_i < infinity
   !           |
   !           | p_i                   if |u_i| = |l_i| = infinity
   !     x_i = |
   !           | u_i-w_i               if l_i = -infinity
   !           |
   !           | z_i+l_i               if u_i = infinity
   !           +
  \*/
  void
  AuxProblem::to_primal( valueType const x[],
                         valueType       xo[],
                         integer         IBo[] ) const {
    integer ib = 0 ;
    integer nzz = 0 ;
    integer nww = 0 ;
    integer npp = 0 ;
    valueConstPointer z = x ;
    valueConstPointer w = z+nz ;
    valueConstPointer p = w+nw ;
    for ( integer i = 0 ; i < n ;  ++i ) {
      integer icase = 0 ;
      if ( pBase->Upper_is_free(i) ) icase  = 1 ;
      if ( pBase->Lower_is_free(i) ) icase += 2 ;
      switch ( icase ) {
      case 0: // l <= x <= u
        xo[i] = (z[map_z[nzz]]-w[map_w[nww]]+(pBase->Upper(i)+pBase->Lower(i)))/2 ;
        IBo[ib++] = i ;
        ++nzz ; ++nww ;
        break ;
      case 1: // l <= x
        xo[i] = pBase->Lower(i) + z[map_z[nzz++]] ;
        break ;
      case 2: // x <= u
        xo[i] = pBase->Upper(i) - w[map_w[nww++]] ;
        break ;
      case 3:
        xo[i] = p[map_p[npp++]] ;
        break ;
      }
    }
  }

  integer
  AuxProblem::load_A_column( integer jcol, valueType vals[], integer irow[] ) const {
    integer nnz ;
    if ( jcol < nz ) { // Z matrix
      integer j = map_z[jcol] ;
      nnz = pBase->load_A_column( j, vals, irow ) ;
      for ( integer i = 0 ; i < nnz ; ++i )
        if ( d[irow[i]] < 0 ) vals[i] = -vals[i] ;
    } else {
      jcol -= nz ;
      if ( jcol < nw ) { // W matrix
        integer j = map_w[jcol] ;
        nnz = pBase->load_A_column( j, vals, irow ) ;
        for ( integer i = 0 ; i < nnz ; ++i )
          if ( d[irow[i]] >= 0 ) vals[i] = -vals[i] ; // reverse sign
      } else {
        jcol -= nw ;
        if ( jcol < np ) { // P matrix
          integer j = map_z[jcol] ;
          nnz = pBase->load_A_column( j, vals, irow ) ;
          for ( integer i = 0 ; i < nnz ; ++i )
            if ( d[irow[i]] < 0 ) vals[i] = -vals[i] ;
        } else { // I matrix
          jcol   -= np ;
          nnz     = 1 ;
          vals[0] = 1 ;
          irow[0] = jcol ;
        }
      }
    }
    return nnz ;
  }

  void
  AuxProblem::subtract_Ax( valueType const x[], valueType res[] ) const {
    for ( integer i = 0 ; i < dim_x() ; ++i ) {
      integer nnz = load_A_column( i, values, i_row ) ;
      for ( integer k = 0 ; k < nnz ; ++k )
        res[i_row[k]] -= x[i]*values[k] ;
    }
  }

  /*\
   |   ___ _                _             _   ___         _    _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  void
  StandardProblem::setup( integer           _m,
                          integer           _n,
                          valueConstPointer _A,
                          integer           _ldA,
                          valueConstPointer _b,
                          valueConstPointer _c,
                          valueConstPointer _L,
                          valueConstPointer _U ) {
    m   = _m ;
    n   = _n ;
    c   = _c ;
    A   = _A ;
    b   = _b ;
    L   = _L ;
    U   = _U ;
    ldA = _ldA ;

    // some check
    SIMPLEX_ASSERT( m >= 0 && n > 0,
                    "Bad problem dimensions, m = " << m << " n = " << n ) ;
  
    SIMPLEX_ASSERT( ldA >= m,
                    "Bad leading dimension of matrix A, ldA = " << ldA << " A is " << m << " x " << n ) ;

    SIMPLEX_ASSERT( n >= m,
                    "Dimension of x (" << n <<
                    ") must be greater than number of equality constraints (" << m << ")" ) ;
    
    L_free.resize(n) ;
    U_free.resize(n) ;
    for ( integer i = 0 ; i < n ; ++i ) {
      L_free[i] = L[i] <= -infinity ;
      U_free[i] = U[i] >=  infinity ;
      if ( !( L_free[i] || U_free[i]) ) {
        SIMPLEX_ASSERT( L[i] <= U[i],
                       "L[i] = " << L[i] << " must be less or equal to U[i] = " << U[i] ) ;
      }
    }

    b_max_abs = alglin::absmax( m, b, 1) ;
    c_max_abs = alglin::absmax( n, c, 1) ;
    A_max_abs = alglin::maxabs( m, n, A, ldA ) ;
  }

  /*\
   |   ___         _    _
   |  | _ \_ _ ___| |__| |___ _ __
   |  |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |_| |_| \___/_.__/_\___|_|_|_|
  \*/

  void
  Problem::setup( integer           _m,
                  integer           _n,
                  valueConstPointer _A,
                  integer           _ldA,
                  valueConstPointer _c,
                  valueConstPointer _L,
                  valueConstPointer _U ) {
    m   = _m ;
    n   = _n ;
    c   = _c ;
    A   = _A ;
    L   = _L ;
    U   = _U ;
    ldA = _ldA ;

    // some check
    SIMPLEX_ASSERT( m >= 0 && n > 0,
                    "Bad problem dimensions, m = " << m << " n = " << n ) ;
  
    SIMPLEX_ASSERT( ldA >= m,
                    "Bad leading dimension of matrix A, ldA = " << ldA << " A is " << m << " x " << n ) ;
    
    L_free.resize(n) ;
    U_free.resize(n) ;
    for ( integer i = 0 ; i < n+m ; ++i ) {
      L_free[i] = L[i] <= -infinity ;
      U_free[i] = U[i] >=  infinity ;
      if ( !( L_free[i] || U_free[i]) ) {
        SIMPLEX_ASSERT( L[i] <= U[i],
                       "L[i] = " << L[i] << " must be less or equal to U[i] = " << U[i] ) ;
      }
    }

    c_max_abs = alglin::absmax( n, c, 1) ;
    A_max_abs = alglin::maxabs( m, n, A, ldA ) ;
  }

  /*\
   |   ___ _                _             _   ___      _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | / __| ___| |_ _____ _ _
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | \__ \/ _ \ \ V / -_) '_|
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |___/\___/_|\_/\___|_|
  \*/
  StandardSolver::StandardSolver( std::string const & name )
  : _name(name)
  , baseReals(name+"-SimplexFull-baseReals")
  , baseIntegers(name+"-SimplexFull-baseInteger")
  , pStream(&std::cout)
  , maxIter(1000)
  {
  }
  
  static
  inline
  void
  writeIter( std::ostream * pStream,
             integer        iter,
             integer const  IB[],
             integer        m ) {
    (*pStream) << "\n====================================================\n"
               << "                    TABLEAU  " << iter
               << "\n====================================================\n"
               << "Basis index set :" ;
    for ( integer i = 0 ; i < m ; ++i ) (*pStream) << " " << IB[i] ;
    (*pStream) << "\n" ;
  }

  std::string
  StandardSolver::Lstring( integer i ) const {
    if ( L_bounded(i) ) {
      std::ostringstream msg ;
      msg << L(i) ;
      return msg.str() ;
    } else {
      return "-Infinity" ;
    }
  }
  
  std::string
  StandardSolver::Ustring( integer i ) const {
    if ( U_bounded(i) ) {
      std::ostringstream msg ;
      msg << U(i) ;
      return msg.str() ;
    } else {
      return "+Infinity" ;
    }
  }

  void
  StandardSolver::solve( StandardProblemBase * _problem,
                         valueType x[],
                         integer   IB[],
                         valueType eps ) {

    problem = _problem ;

    n = problem->dim_x() ;
    m = problem->dim_g() ;

    // memory allocation
    baseReals.allocate(4*m+2*n) ;
    valueType * gamma  = baseReals(m) ;
    valueType * beta   = baseReals(m) ;
    valueType * ceta   = baseReals(n-m) ;
    valueType * y      = baseReals(m) ;
    valueType * c      = baseReals(n) ;
    valueType * values = baseReals(m) ;

    baseIntegers.allocate(n) ;
    integer * IN    = baseIntegers(n-m) ;
    integer * i_row = baseIntegers(m) ;

    problem->load_c( c ) ;

    // no equality constraints, solution is on the border of the box
    if ( m == 0 ) { // trivial solution
      for ( integer i = 0 ; i < n ; ++i ) {
        if      ( c[i] > 0 ) x[i] = L(i) ;
        else if ( c[i] < 0 ) x[i] = U(i) ;
        else {
          if      ( L(i) > 0 ) x[i] = L(i) ;
          else if ( U(i) < 0 ) x[i] = U(i) ;
          else                 x[i] = 0 ;
        }
      }
      return ;
    }
    
    // build non bases index set IN
    std::sort( IB, IB+m ) ;
    SIMPLEX_ASSERT( IB[0] >= 0 && IB[m-1] < n,
                    "Error in Simplex: Wrong Basis index set" ) ;
    for ( integer i = 1 ; i < m ; ++i ) {
      SIMPLEX_ASSERT( IB[i-1] != IB[i],
                      "Error in Simplex: Duplicated index in basis index set" ) ;
    }

    integer ib = 0 ;
    integer in = 0 ;
    for ( integer i = 0 ; i < n ; ++i ) {
      while ( ib < m && IB[ib] < i ) ++ib ;
      if ( ib == m || IB[ib] != i ) IN[in++] = i ;
      SIMPLEX_ASSERT( in <= n-m, "Bad non basis index construction" ) ;
    }

    if ( n == m ) { // trivial solution
      alglin::LU<valueType> lu ; // qr decomposition manage class
      // select ALL the column of matrix A
      for ( integer j = 0 ; j < m ; ++j ) {
        integer nnz = problem->load_A_column( j, values, i_row ) ;
        lu.load_sparse_column( nnz, values, i_row, j ) ;
      }
      lu.factorize() ;
      problem->load_b( x ) ;
      lu.t_solve( x ) ;
      // check solution
      for ( integer i = 0 ; i < n ; ++i ) {
        bool ok = true ;
        if ( U_bounded(i) ) ok = x[i] <= U(i)+eps ;
        if ( L_bounded(i) ) ok = ok && x[i] >= L(i)-eps ;
        SIMPLEX_ASSERT( ok,
                        "Error in Simplex: Infeasible solution x[" << i <<
                        "] = " << x[i] << " not in [" <<
                        Lstring(i) << "," << Ustring(i) << "]" ) ;
      }
    }

    alglin::QRP<valueType> qr ; // qr decomposition manage class
    qr.allocate( m, m ) ;
    
    // check basis
    for ( integer i = 0 ; i < m ; ++i ) {
      integer ibi = IB[i] ;
      bool ok = true ;
      if ( U_bounded(ibi) ) ok = x[ibi] <= U(ibi)+eps ;
      if ( L_bounded(ibi) ) ok = ok && x[ibi] >= L(ibi)-eps ;
      SIMPLEX_ASSERT( ok,
                      "Error in Simplex: Infeasible starting point, violating basis variable x[" << ibi <<
                      "] = " << x[ibi] << " not in [" <<
                      Lstring(ibi) << ", " << Lstring(ibi) << "]" ) ;
    }

    // check non basis
    for ( integer i = 0 ; i < n-m ; ++i ) {
      integer ini = IN[i] ;
      bool ok = false ;
      if ( U_bounded(ini) && x[ini] >= U(ini)-eps ) {
        x[ini] = U(ini) ;
        ok     = true ;
      } else if ( L_bounded(ini) && x[ini] <= L(ini)+eps ) {
        x[ini] = L(ini) ;
        ok     = true ;
      }
      SIMPLEX_ASSERT( ok, "Error in Simplex: Infeasible starting point, violating non basis variable x[" << ini <<
                          "] = " << x[ini] << " not equal to " <<
                          Lstring(ini) << " or " << Ustring(ini) ) ;
    }

    // check residual
    valueType scale_residual = 1 ;
    if ( problem->get_A_max_abs() > scale_residual ) scale_residual = problem->get_A_max_abs() ;
    if ( problem->get_b_max_abs() > scale_residual ) scale_residual = problem->get_b_max_abs() ;

    problem->load_b( y ) ;
    problem->subtract_Ax( x, y ) ;
    for ( integer i = 0 ; i < m ; ++i ) {
      valueType error = std::abs(y[i]) ;
      SIMPLEX_ASSERT( error < eps*scale_residual,
                      "Error in Simplex: Infeasible starting point, violating " << i <<
                      "th equality constraint, error = " << y[i] ) ;
    }
    /*
    //  _
    // | |___  ___ _ __
    // | / _ \/ _ \ '_ \
    // |_\___/\___/ .__/
    //            |_|
    */
    for ( integer iter = 0 ; iter < maxIter ; ++iter ) {
      // select the column of matrix A
      for ( integer i = 0 ; i < m ; ++i ) {
        integer ibi = IB[i] ;
        integer nnz = problem->load_A_column( ibi, values, i_row ) ;
        qr.load_sparse_column( nnz, values, i_row, i ) ;
        y[i] = c[ibi] ;
      }
      if ( pStream != nullptr ) {
        writeIter( pStream, iter, IB, m ) ;
        (*pStream) << "Objective c'x = " << alglin::dot( n, x, 1, c, 1 ) << "\n" ;
      }

      /*
      // determine pivot column:
      */
      // 1.  QR decomposition of AB (without pivoting)
      // 2.  compute dual variable Y from  AB' * Y = CB  <=>  R' * Z = CB, Z = Q' * Y
      // 2B. compute Y = Q * Z
      qr.factorize() ;
      integer rank = qr.rankEstimate(relaxedEpsilon) ;
      SIMPLEX_ASSERT( rank == m,
                      "Error in Simplex: Basis Matrix Singular, rank = " <<
                      rank << " should be " << m ) ;
      qr.t_solve( y ) ;
      // 3. compute pivot column gamma
      integer   inmin      = n ;
      integer   Q          = -1 ;
      valueType ceta_max_L = 0 ;
      valueType ceta_min_U = 0 ;
      // select the non bases (pivot) columns (n-m)
      for ( integer i = 0 ; i < n-m ; ++i ) {
        integer ini = IN[i] ;
        // 3A.  compute negative reduced costs Y'*A_N - C_N
        //problem->load_A_column( column, ini ) ;
        //ceta[i] = alglin::dot(m,y,1,column,1) - c[ini] ;
        integer nnz = problem->load_A_column( ini, values, i_row ) ;
        ceta[i] = 0 ;
        for ( integer j = 0 ; j < nnz ; ++j ) ceta[i] += y[i_row[j]] * values[j] ;
        ceta[i] -= c[ini] ;
        if ( alglin::isZero(x[ini]-L(ini)) ) {
          if ( ceta_max_L < ceta[i] ) ceta_max_L = ceta[i] ;
          if ( ini < inmin && ceta[i] >  eps ) { inmin = ini ; Q = i ; }
        } else if ( alglin::isZero(x[ini]-U(ini)) ) {
          if ( ceta_min_U > ceta[i] ) ceta_min_U = ceta[i] ;
          if ( ini < inmin && ceta[i] < -eps ) { inmin = ini ; Q = i ; }
        } else {
          if ( pStream != nullptr )
            (*pStream) << "Non basis variable " << ini << " not at the boundary\n" ;
        }
      }
      // 4.  compute basic solution
      problem->load_b( beta ) ;
      for ( integer j = 0 ; j < n-m ; ++j ) {
        integer inj = IN[j] ;
        integer nnz = problem->load_A_column( inj, values, i_row ) ;
        for ( integer k = 0 ; k < nnz ; ++k )
          beta[i_row[k]] -= x[inj] * values[k] ;
      }
      qr.solve( beta ) ;
      // check for optimality (no new pivot selected)
      if ( Q < 0 ) {
        /*
        // return basic solution
        */
        for ( integer i = 0 ; i < m ; ++i ) x[IB[i]] = beta[i] ;
        /*
        // output
        */
        if ( pStream != nullptr ) {
          bool unique = std::abs(ceta_max_L) <= eps ||
                        std::abs(ceta_min_U) <= eps ;
          std::vector<bool>      B(n) ;
          std::vector<valueType> D(n) ;
          std::fill( B.begin(), B.end(), false ) ;
          for ( integer i = 0 ; i < m ; ++i ) {
            B[IB[i]] = true ;
            D[IB[i]] = y[i] ;
          }
          for ( integer i = 0 ; i < n-m ; ++i ) {
            integer ini = IN[i] ;
            if ( alglin::isZero(x[ini]-L(ini)) ) {
              D[ini] = ceta[i] ;
            } else if ( alglin::isZero(x[ini]-U(ini)) ) {
              D[ini] = -ceta[i] ;
            } else {
              SIMPLEX_ERROR( "Something wrong in Simplex, " <<
                         " x[" << ini << "] = " << x[ini] << " != L[" << ini << "] = " << L(ini) << " and " <<
                         " x[" << ini << "] = " << x[ini] << " != U[" << ini << "] = " << U(ini) <<
                         " for non basis variable") ;
            }
          }
          //writeIter( pStream, iter, IB, m ) ;
          (*pStream)
            << "\n\n"
            << std::setw(4)  << "TYP"
            << std::setw(4)  << "I"
            << std::setw(14) << "Lower"
            << std::setw(14) << "X"
            << std::setw(14) << "Upper"
            << std::setw(14) << "Dual"
            << '\n' ;
          for ( integer i = 0 ; i < n ; ++i )
            (*pStream) << std::setw(4)  << (B[i]?"B":"N")
                       << std::setw(4)  << i
                       << std::setw(14) << Lstring(i)
                       << std::setw(14) << x[i]
                       << std::setw(14) << Ustring(i)
                       << std::setw(14) << D[i]
                       << (B[i]?" Ax=b\n":"\n") ;
          (*pStream) << "\nObjective c'x = " << alglin::dot( n, x, 1, c, 1 ) << "\n" ;
          if ( !unique ) (*pStream) << "\nThe solution is NOT unique\n\n" ;
        }
        return ;
      }
      /*
      // pivot row
      */
      integer nnz = problem->load_A_column( IN[Q], values, i_row ) ;
      alglin::zero( m, gamma, 1 ) ;
      for ( integer j = 0 ; j < nnz ; ++j ) gamma[i_row[j]] = values[j] ;
      //problem->load_A_column( gamma, IN[Q] ) ;
      /*
      // 5.  compute pivot column
      */
      qr.solve( gamma ) ;
      /*
      // determine pivot row
      */
      integer   IBmin    = n ; // maximum number of column of A (size of x)
      integer   P        = -1 ;
      valueType beta_min = infinity ;
      if ( ceta[Q] > 0 ) {
        /*
        // case 1: ceta[Q] > 0 (resp. X[IN[Q]]=L[IN[Q]])
        */
        for ( integer i = 0 ; i < m ; ++i ) {
          integer   IBi  = IB[i] ;
          valueType QUOT = infinity ;
          if ( gamma[i] < -eps ) {
            if ( U_bounded(IBi) ) QUOT = (beta[i]-U(IBi))/gamma[i] ;
          } else if ( gamma[i] > eps ) {
            if ( L_bounded(IBi) ) QUOT = (beta[i]-L(IBi))/gamma[i] ;
          }
          if ( QUOT < beta_min ) {
            beta_min = QUOT ; IBmin = IBi ; P = i ;
          }
          if ( std::abs(QUOT-beta_min) < eps && IBi < IBmin ) {
            IBmin = IBi ; P = i ;
          }
        }
      } else if ( ceta[Q] < 0 ) {
        /*
        // case 2: ceta[Q] < 0 (resp. x[IN[Q]]=U[IN[Q]])
        */
        for ( integer i = 0 ; i < m ; ++i ) {
          integer   IBi  = IB[i] ;
          valueType QUOT = infinity ;
          if ( gamma[i] > eps ) {
            if ( U_bounded(IBi) ) QUOT = (U(IBi)-beta[i])/gamma[i] ;
          } else if ( gamma[i] < -infinity ) {
            if ( L_bounded(IBi) ) QUOT = (L(IBi)-beta[i])/gamma[i] ;
          }
          if ( QUOT < beta_min ) {
            beta_min = QUOT ; IBmin = IBi ; P = i ;
          }
          if ( std::abs(QUOT-beta_min) < eps && IBi < IBmin ) {
            IBmin = IBi ; P = i ;
          }
        }
      }
      /*
      // stopping criteria
      */
      valueType TQ  = infinity ;
      integer   INq = IN[Q] ;
      integer   IBp = IB[P] ;
      if ( U_bounded(INq) && L_bounded(INq) ) TQ = U(INq)-L(INq) ;

      SIMPLEX_ASSERT( TQ < infinity || beta_min < infinity, "Problem is NOT solvable!" ) ;
      if ( TQ < infinity && TQ <= beta_min ) {
        if      ( ceta[Q] > eps  ) x[INq] = U(INq) ;
        else if ( ceta[Q] < -eps ) x[INq] = L(INq) ;
        SIMPLEX_ASSERT( ceta[Q] <= eps || ceta[Q] >= -eps,
                        "Warning in Simplex: Could not find a new pivot row in iteration " << iter ) ;
      } else if ( beta_min < infinity && beta_min < TQ ) {
        if ( ceta[Q] > 0 ) {
          if      ( gamma[P] < 0 ) x[IBp] = U(IBp);
          else if ( gamma[P] > 0 ) x[IBp] = L(IBp);
        } else if ( ceta[Q] < 0 ) {
          if      ( gamma[P] < 0 ) x[IBp] = L(IBp);
          else if ( gamma[P] > 0 ) x[IBp] = U(IBp);
        }
        if ( pStream != nullptr )
          (*pStream)
            << "Pivot element gamma(" << IB[P] << "," << IN[Q] << ") = " << gamma[P] << "\n" ;
        /*
        // perform basis change
        */
        std::swap( IB[P], IN[Q] ) ;
      } else {
        SIMPLEX_ERROR( "Error in Simplex: No continuation defined!" ) ;
      }
    }
  }

}