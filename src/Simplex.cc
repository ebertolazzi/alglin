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
#include "Alglin_Eigen.hh"

//! namespace for nonlinear systems and nonlinearsolver
namespace Simplex {

  using std::string;
  using std::abs;
  using std::numeric_limits;
  using std::sort;
  using std::swap;
  using std::cout;

  real_type const epsilon        = numeric_limits<real_type>::epsilon();
  real_type const relaxedEpsilon = pow(epsilon,0.8);
  real_type const infinity       = numeric_limits<real_type>::max();

 /*\
  |     _             ___         _    _
  |    /_\ _  ___ __ | _ \_ _ ___| |__| |___ _ __
  |   / _ \ || \ \ / |  _/ '_/ _ \ '_ \ / -_) '  \
  |  /_/ \_\_,_/_\_\ |_| |_| \___/_.__/_\___|_|_|_|
 \*/
  void
  AuxProblem::setup( StandardProblemBase * _pBase ) {

    m_problem_base = _pBase;

    m_n = m_problem_base->dim_x();
    m_m = m_problem_base->dim_g();

    // compute size of variable z, w, and p
    m_nz = m_nw = m_np = 0;
    for ( integer i{0}; i < m_n; ++i ) {
      integer icase{0};
      if ( m_problem_base->Upper_is_free(i) ) icase  = 1;
      if ( m_problem_base->Lower_is_free(i) ) icase += 2;
      switch ( icase ) {
      case 0: ++m_nz; ++m_nw; break; // l <= x <= u
      case 1: ++m_nz;         break; // l <= x
      case 2: ++m_nw;         break; // x <= u
      case 3: ++m_np;         break;
      }
    }
    // allocate memory
    m_mem_int.reallocate( size_t(m_nz+m_nw+m_np+m_n+m_m) );
    m_map_z    = m_mem_int( size_t(m_nz) );
    m_map_w    = m_mem_int( size_t(m_nw) );
    m_map_p    = m_mem_int( size_t(m_np) );
    m_map_case = m_mem_int( size_t(m_n)  );
    m_i_row    = m_mem_int( size_t(m_m)  );

    m_mem.reallocate( size_t(2*m_m) );
    m_d      = m_mem( size_t(m_m) );
    m_values = m_mem( size_t(m_m) );

    // fill mapping and vector q
    m_nz = m_nw = m_np = 0;
    m_problem_base->load_b( m_d );
    for ( integer i{0}; i < m_n;  ++i ) {
      real_type q{0};
      integer icase{0};
      if ( m_problem_base->Upper_is_free(i) ) icase  = 1;
      if ( m_problem_base->Lower_is_free(i) ) icase += 2;
      m_map_case[i] = icase;
      switch ( icase ) {
      case 0: // l <= x <= u
        m_map_z[m_nz++] = i;
        m_map_w[m_nw++] = i;
        q = (m_problem_base->Upper(i) + m_problem_base->Lower(i))/2;
        break;
      case 1: // l <= x
        m_map_z[m_nz++] = i;
        q = m_problem_base->Lower(i);
        break;
      case 2: // x <= u
        m_map_w[m_nw++] = i;
        q = m_problem_base->Upper(i);
        break;
      case 3:
        m_map_p[m_np++] = i;
        break;
      }
      integer nnz = m_problem_base->load_A_column( i, m_values, m_i_row );
      for ( integer k = 0; k < nnz; ++k )
        m_d[m_i_row[k]] -= q*m_values[k];
    }

    m_b_max_abs = alglin::absmax( m_m, m_d, 1);
    m_c_max_abs = 1;
    m_A_max_abs = m_problem_base->get_A_max_abs();

  }

  bool
  AuxProblem::Lower_is_free( integer i ) const
  { return i >= m_nz+m_nw && i < m_nz+m_nw+m_np; }

  bool
  AuxProblem::Upper_is_free( integer   ) const
  { return true; }

  real_type
  AuxProblem::Lower( integer i ) const {
    if ( Lower_is_free(i) ) return infinity;
    else                    return 0;
  }

  real_type
  AuxProblem::Upper( integer ) const
  { return infinity; }

  void
  AuxProblem::feasible_point( real_type x[], integer IB[] ) const {
    integer nn = m_nz+m_nw+m_np;
    alglin::Zero_n( x, nn );
    for ( integer i{0}; i < m_m; ++i ) {
      x[nn+i] = abs(m_d[i]);
      IB[i] = nn+i;
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
  AuxProblem::to_primal(
    real_type const x[],
    real_type       xo[],
    integer         IBo[]
  ) const {
    integer ib{0};
    integer nzz{0};
    integer nww{0};
    integer npp{0};
    real_type const * z{x};
    real_type const * w{z+m_nz};
    real_type const * p{w+m_nw};
    for ( integer i{0}; i < m_n;  ++i ) {
      integer icase{0};
      if ( m_problem_base->Upper_is_free(i) ) icase  = 1;
      if ( m_problem_base->Lower_is_free(i) ) icase += 2;
      switch ( icase ) {
      case 0: // l <= x <= u
        xo[i] = ( z[m_map_z[nzz]] - w[m_map_w[nww]] +
                  (m_problem_base->Upper(i)+m_problem_base->Lower(i)) )/2;
        IBo[ib++] = i;
        ++nzz; ++nww;
        break;
      case 1: // l <= x
        xo[i] = m_problem_base->Lower(i) + z[m_map_z[nzz++]];
        break;
      case 2: // x <= u
        xo[i] = m_problem_base->Upper(i) - w[m_map_w[nww++]];
        break;
      case 3:
        xo[i] = p[m_map_p[npp++]];
        break;
      }
    }
  }

  integer
  AuxProblem::load_A_column( integer jcol, real_type vals[], integer irow[] ) const {
    integer nnz;
    if ( jcol < m_nz ) { // Z matrix
      integer j{m_map_z[jcol]};
      nnz = m_problem_base->load_A_column( j, vals, irow );
      for ( integer i{0}; i < nnz; ++i )
        if ( m_d[irow[i]] < 0 ) vals[i] = -vals[i];
    } else {
      jcol -= m_nz;
      if ( jcol < m_nw ) { // W matrix
        integer j{m_map_w[jcol]};
        nnz = m_problem_base->load_A_column( j, vals, irow );
        for ( integer i{0}; i < nnz; ++i )
          if ( m_d[irow[i]] >= 0 ) vals[i] = -vals[i]; // reverse sign
      } else {
        jcol -= m_nw;
        if ( jcol < m_np ) { // P matrix
          integer j{m_map_z[jcol]};
          nnz = m_problem_base->load_A_column( j, vals, irow );
          for ( integer i{0}; i < nnz; ++i )
            if ( m_d[irow[i]] < 0 ) vals[i] = -vals[i];
        } else { // I matrix
          jcol   -= m_np;
          nnz     = 1;
          vals[0] = 1;
          irow[0] = jcol;
        }
      }
    }
    return nnz;
  }

  void
  AuxProblem::subtract_Ax( real_type const x[], real_type res[] ) const {
    for ( integer i{0}; i < dim_x(); ++i ) {
      integer nnz{ load_A_column( i, m_values, m_i_row )} ;
      for ( integer k{0}; k < nnz; ++k )
        res[m_i_row[k]] -= x[i]*m_values[k];
    }
  }

  /*\
   |   ___ _                _             _   ___         _    _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  void
  StandardProblem::setup(
    integer         m,
    integer         n,
    real_type const A[],
    integer         ldA,
    real_type const b[],
    real_type const c[],
    real_type const L[],
    real_type const U[]
  ) {
    m_m      = m;
    m_n      = n;
    m_c      = c;
    m_A      = A;
    m_b      = b;
    m_L      = L;
    m_U      = U;
    m_A_ldim = ldA;

    // some check
    UTILS_ASSERT(
      m_m >= 0 && m_n > 0,
      "StandardProblem::Bad problem dimensions, m={} n={}\n",
      m_m, m_n
    );

    UTILS_ASSERT(
      m_A_ldim >= m_m,
      "StandardProblem::Bad leading dimension of matrix A, ldA={} A is {} x {}\n",
      m_A_ldim, m_m, m_n
    );

    UTILS_ASSERT(
      m_n >= m_m,
      "StandardProblem::Dimension of x ({}) must be greater than\n"
      "number of equality constraints ({})\n",
      m_m, m_n
    );

    m_L_free.resize( size_t(m_n) );
    m_U_free.resize( size_t(m_n) );
    for ( integer i{0}; i < m_n; ++i ) {
      m_L_free[i] = m_L[i] <= -infinity;
      m_U_free[i] = m_U[i] >=  infinity;
      if ( !( m_L_free[i] || m_U_free[i]) ) {
        UTILS_ASSERT(
          m_L[i] <= m_U[i],
          "StandardProblem, L[i]={} must be less or equal to U[i]={}\n",
          m_L[i], m_U[i]
        );
      }
    }

    m_b_max_abs = alglin::absmax( m_m, m_b, 1);
    m_c_max_abs = alglin::absmax( m_n, m_c, 1);
    m_A_max_abs = alglin::maxabs( m_m, m_n, m_A, m_A_ldim );
  }

  /*\
   |   ___         _    _
   |  | _ \_ _ ___| |__| |___ _ __
   |  |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |_| |_| \___/_.__/_\___|_|_|_|
  \*/

  void
  Problem::setup(
    integer         m,
    integer         n,
    real_type const A[],
    integer         ldA,
    real_type const c[],
    real_type const L[],
    real_type const U[]
  ) {
    m_m      = m;
    m_n      = n;
    m_c      = c;
    m_A      = A;
    m_L      = L;
    m_U      = U;
    m_A_ldim = ldA;

    // some check
    UTILS_ASSERT(
      m_m >= 0 && m_n > 0,
      "Problem::Bad problem dimensions, m={} n={}\n",
      m_m, m_n
    );

    UTILS_ASSERT(
      m_A_ldim >= m_m,
      "Problem::Bad leading dimension of matrix A, ldA={} A is {} x {}\n",
      m_A_ldim, m_m, m_n
    );

    m_L_free.resize(m_n);
    m_U_free.resize(m_n);
    for ( integer i{0}; i < m_n+m_m; ++i ) {
      m_L_free[i] = m_L[i] <= -infinity;
      m_U_free[i] = m_U[i] >=  infinity;
      if ( !( m_L_free[i] || m_U_free[i]) ) {
        UTILS_ASSERT(
          m_L[i] <= m_U[i],
          "Problem, L[i]={} must be less or equal to U[i]={}\n",
          m_L[i], m_U[i]
        );
      }
    }

    m_c_max_abs = alglin::absmax( m_n, m_c, 1);
    m_A_max_abs = alglin::maxabs( m_m, m_n, m_A, m_A_ldim );
  }

  /*\
   |   ___ _                _             _   ___      _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | / __| ___| |_ _____ _ _
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | \__ \/ _ \ \ V / -_) '_|
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |___/\___/_|\_/\___|_|
  \*/
  StandardSolver::StandardSolver( string const & name )
  : m_name(name)
  , m_mem(fmt::format("StandardSolver[{}]::m_mem",name))
  , m_mem_int(fmt::format("StandardSolver[{}]::m_mem_int",name))
  , pStream(&cout)
  {
  }

  static
  inline
  void
  writeIter(
    ostream_type * pStream,
    integer        iter,
    integer const  IB[],
    integer        m
  ) {
    fmt::print( *pStream,
      "\n"
      "====================================================\n"
      "                    TABLEAU {}\n"
      "====================================================\n"
      "Basis index set :",
      iter
    );
    for ( integer i{0}; i < m; ++i ) fmt::print( *pStream, "{} ", IB[i] );
    fmt::print( *pStream, "\n" );
  }

  string
  StandardSolver::Lstring( integer i ) const {
    if ( L_bounded(i) ) {
      return fmt::format("{}",L(i));
    } else {
      return "-Infinity";
    }
  }

  string
  StandardSolver::Ustring( integer i ) const {
    if ( U_bounded(i) ) {
      return fmt::format("{}",U(i));
    } else {
      return "+Infinity";
    }
  }

  // windows workaround!
  #ifdef IN
    #undef IN
  #endif
  #ifdef IB
    #undef IB
  #endif

  void
  StandardSolver::solve(
    StandardProblemBase * _problem,
    real_type x[],
    integer   IB[],
    real_type eps
  ) {

    m_problem = _problem;

    m_n = m_problem->dim_x();
    m_m = m_problem->dim_g();
    integer nm{m_n-m_m};

    // memory allocation
    m_mem.reallocate(4*m_m+2*m_n);
    real_type * gamma  { m_mem(m_m) };
    real_type * beta   { m_mem(m_m) };
    real_type * ceta   { m_mem(nm)  };
    real_type * y      { m_mem(m_m) };
    real_type * c      { m_mem(m_n) };
    real_type * values { m_mem(m_m) };

    m_mem_int.reallocate(m_n);
    integer * IN    { m_mem_int(nm)  };
    integer * i_row { m_mem_int(m_m) };

    m_problem->load_c( c );

    // no equality constraints, solution is on the border of the box
    if ( m_m == 0 ) { // trivial solution
      for ( integer i{0}; i < m_n; ++i ) {
        if      ( c[i] > 0 ) x[i] = L(i);
        else if ( c[i] < 0 ) x[i] = U(i);
        else {
          if      ( L(i) > 0 ) x[i] = L(i);
          else if ( U(i) < 0 ) x[i] = U(i);
          else                 x[i] = 0;
        }
      }
      return;
    }

    // build non bases index set IN
    sort( IB, IB+m_m );
    UTILS_ASSERT0(
      IB[0] >= 0 && IB[m_m-1] < m_n,
      "StandardSolver. Wrong Basis index set\n"
    );
    for ( integer i{1}; i < m_m; ++i ) {
      UTILS_ASSERT0(
        IB[i-1] != IB[i],
        "StandardSolver. Duplicated index in basis index set\n"
      );
    }

    integer ib{0};
    integer in{0};
    for ( integer i{0}; i < m_n; ++i ) {
      while ( ib < m_m && IB[ib] < i ) ++ib;
      if ( ib == m_m || IB[ib] != i ) IN[in++] = i;
      UTILS_ASSERT0( in <= nm, "Bad non basis index construction\n" );
    }

    if ( m_n == m_m ) { // trivial solution
      alglin::Matrix<real_type> mat;
      alglin::LU<real_type>     lu; // qr decomposition manage class
      mat.setup( m_m, m_m );
      mat.zero_fill();
      // select ALL the column of matrix A
      for ( integer j{0}; j < m_m; ++j ) {
        integer nnz{ m_problem->load_A_column( j, values, i_row ) };
        mat.load_sparse_column( nnz, values, i_row, j );
      }
      lu.factorize( "StandardSolver::solve<lu>", mat );

      m_problem->load_b( x );
      lu.t_solve( x );
      // check solution
      for ( integer i{0}; i < m_n; ++i ) {
        bool ok{true};
        if ( U_bounded(i) ) ok = x[i] <= U(i)+eps;
        if ( L_bounded(i) ) ok = ok && x[i] >= L(i)-eps;
        UTILS_ASSERT(
          ok,
          "StandardSolver. Infeasible solution x[{}]={} not in [{},{}]\n",
          i, x[i], Lstring(i), Ustring(i)
        );
      }
    }

    alglin::Matrix<real_type> qr_mat;
    alglin::QRP<real_type>    qr; // qr decomposition manage class
    qr_mat.setup( m_m, m_m );
    qr_mat.zero_fill();

    // check basis
    for ( integer i{0}; i < m_m; ++i ) {
      integer ibi = IB[i];
      bool ok = true;
      if ( U_bounded(ibi) ) ok = x[ibi] <= U(ibi)+eps;
      if ( L_bounded(ibi) ) ok = ok && x[ibi] >= L(ibi)-eps;
      UTILS_ASSERT(
        ok,
        "StandardSolver. Infeasible starting point,\n"
        "violating basis variable x[{}] = {} not in [{},{}]\n",
        ibi, x[ibi], Lstring(ibi), Lstring(ibi)
      );
    }

    // check non basis
    for ( integer i{0}; i < nm; ++i ) {
      integer ini = IN[i];
      bool ok{false};
      if ( U_bounded(ini) && x[ini] >= U(ini)-eps ) {
        x[ini] = U(ini);
        ok     = true;
      } else if ( L_bounded(ini) && x[ini] <= L(ini)+eps ) {
        x[ini] = L(ini);
        ok     = true;
      }
      UTILS_ASSERT(
        ok,
        "StandardSolver. Infeasible starting point,\n"
        "violating non basis variable x[{}]={} not equal to {} or {}\n",
        ini, x[ini], Lstring(ini), Ustring(ini)
      );
    }

    // check residual
    real_type scale_residual{1};
    if ( m_problem->get_A_max_abs() > scale_residual ) scale_residual = m_problem->get_A_max_abs();
    if ( m_problem->get_b_max_abs() > scale_residual ) scale_residual = m_problem->get_b_max_abs();

    m_problem->load_b( y );
    m_problem->subtract_Ax( x, y );
    for ( integer i{0}; i < m_m; ++i ) {
      real_type error{abs(y[i])};
      UTILS_ASSERT(
        error < eps*scale_residual,
        "StandardSolver. Infeasible starting point,\n"
        "violating {}th equality constraint, error={}\n",
        i, y[i]
      );
    }
    /*
    //  _
    // | |___  ___ _ __
    // | / _ \/ _ \ '_ \
    // |_\___/\___/ .__/
    //            |_|
    */
    for ( integer iter{0}; iter < m_max_iter; ++iter ) {
      // select the column of matrix A
      for ( integer i{0}; i < m_m; ++i ) {
        integer ibi{IB[i]};
        integer nnz{m_problem->load_A_column( ibi, values, i_row )};
        qr_mat.load_sparse_column( nnz, values, i_row, i );
        y[i] = c[ibi];
      }
      if ( pStream != nullptr ) {
        writeIter( pStream, iter, IB, m_m );
        fmt::print( *pStream,
          "Objective c'x = {}\n", alglin::dot( m_n, x, 1, c, 1 )
        );
      }

      /*
      // determine pivot column:
      */
      // 1.  QR decomposition of AB (without pivoting)
      // 2.  compute dual variable Y from  AB' * Y = CB  <=>  R' * Z = CB, Z = Q' * Y
      // 2B. compute Y = Q * Z
      qr.factorize( "StandardSolver::solve<qr>", qr_mat );
      integer rank{qr.rankEstimate(relaxedEpsilon)};
      UTILS_ASSERT(
        rank == m_m,
        "StandardSolver. Basis Matrix Singular, rank={} should be {}\n",
        rank, m_m
      );
      qr.t_solve( y );
      // 3. compute pivot column gamma
      integer   inmin{m_n};
      integer   Q{-1};
      real_type ceta_max_L{0};
      real_type ceta_min_U{0};
      // select the non bases (pivot) columns (n-m)
      for ( integer i{0}; i < nm; ++i ) {
        integer ini{IN[i]};
        // 3A.  compute negative reduced costs Y'*A_N - C_N
        //problem->load_A_column( column, ini );
        //ceta[i] = alglin::dot(m,y,1,column,1) - c[ini];
        integer nnz{m_problem->load_A_column( ini, values, i_row )};
        ceta[i] = 0;
        for ( integer j{0}; j < nnz; ++j ) ceta[i] += y[i_row[j]] * values[j];
        ceta[i] -= c[ini];
        if ( Utils::is_zero(x[ini]-L(ini)) ) {
          if ( ceta_max_L < ceta[i] ) ceta_max_L = ceta[i];
          if ( ini < inmin && ceta[i] >  eps ) { inmin = ini; Q = i; }
        } else if ( Utils::is_zero(x[ini]-U(ini)) ) {
          if ( ceta_min_U > ceta[i] ) ceta_min_U = ceta[i];
          if ( ini < inmin && ceta[i] < -eps ) { inmin = ini; Q = i; }
        } else {
          if ( pStream != nullptr )
            fmt::print( *pStream,  "Non basis variable {} not at the boundary\n", ini );
        }
      }
      // 4.  compute basic solution
      m_problem->load_b( beta );
      for ( integer j{0}; j < nm; ++j ) {
        integer inj{IN[j]};
        integer nnz{m_problem->load_A_column( inj, values, i_row )};
        for ( integer k{0}; k < nnz; ++k )
          beta[i_row[k]] -= x[inj] * values[k];
      }
      qr.solve( beta );
      // check for optimality (no new pivot selected)
      if ( Q < 0 ) {
        /*
        // return basic solution
        */
        for ( integer i{0}; i < m_m; ++i ) x[IB[i]] = beta[i];
        /*
        // output
        */
        if ( pStream != nullptr ) {
          bool unique{ abs(ceta_max_L) <= eps || abs(ceta_min_U) <= eps };
          vector<bool>      B(m_n);
          vector<real_type> D(m_n);
          fill( B.begin(), B.end(), false );
          for ( integer i{0}; i < m_m; ++i ) {
            B[IB[i]] = true;
            D[IB[i]] = y[i];
          }
          for ( integer i{0}; i < nm; ++i ) {
            integer ini{IN[i]};
            if ( Utils::is_zero(x[ini]-L(ini)) ) {
              D[ini] = ceta[i];
            } else if ( Utils::is_zero(x[ini]-U(ini)) ) {
              D[ini] = -ceta[i];
            } else {
              UTILS_ERROR(
                "StandardSolver. Something wrong\n"
                " x[{}]={} != L[{}]={} and x[{}]={} != U[{}]={}\n"
                "for non basis variable\n",
                ini, x[ini], ini, L(ini), ini, x[ini], ini, U(ini)
              );
            }
          }
          //writeIter( pStream, iter, IB, m );
          fmt::print( *pStream,
            "\n\n{:4}{:4}{:14}{:14}{:14}\n",
            "TYP", "I", "Lower", "X", "Upper", "Dual"
          );
          for ( integer i{0}; i < m_n; ++i )
            fmt::print( *pStream,
              "\n\n{:4}{:4}{:14}{:14}{:14}{}",
               (B[i]?"B":"N"),
               i, Lstring(i), x[i], Ustring(i), D[i],
               (B[i]?" Ax=b\n":"\n")
            );
          fmt::print( *pStream,
            "\nObjective c'x = {}\n",
            alglin::dot( m_n, x, 1, c, 1 )
          );
          if ( !unique ) fmt::print( *pStream, "\nThe solution is NOT unique\n\n" );
        }
        return;
      }
      /*
      // pivot row
      */
      integer nnz{m_problem->load_A_column( IN[Q], values, i_row )};
      alglin::Zero_n( gamma, m_m );
      for ( integer j{0}; j < nnz; ++j ) gamma[i_row[j]] = values[j];
      //problem->load_A_column( gamma, IN[Q] );
      /*
      // 5.  compute pivot column
      */
      qr.solve( gamma );
      /*
      // determine pivot row
      */
      integer   IBmin{m_n}; // maximum number of column of A (size of x)
      integer   P{-1};
      real_type beta_min{infinity};
      if ( ceta[Q] > 0 ) {
        /*
        // case 1: ceta[Q] > 0 (resp. X[IN[Q]]=L[IN[Q]])
        */
        for ( integer i{0}; i < m_m; ++i ) {
          integer   IBi{IB[i]};
          real_type QUOT{infinity};
          if ( gamma[i] < -eps ) {
            if ( U_bounded(IBi) ) QUOT = (beta[i]-U(IBi))/gamma[i];
          } else if ( gamma[i] > eps ) {
            if ( L_bounded(IBi) ) QUOT = (beta[i]-L(IBi))/gamma[i];
          }
          if ( QUOT < beta_min ) {
            beta_min = QUOT; IBmin = IBi; P = i;
          }
          if ( abs(QUOT-beta_min) < eps && IBi < IBmin ) {
            IBmin = IBi; P = i;
          }
        }
      } else if ( ceta[Q] < 0 ) {
        /*
        // case 2: ceta[Q] < 0 (resp. x[IN[Q]]=U[IN[Q]])
        */
        for ( integer i{0}; i < m_m; ++i ) {
          integer   IBi{IB[i]};
          real_type QUOT{infinity};
          if ( gamma[i] > eps ) {
            if ( U_bounded(IBi) ) QUOT = (U(IBi)-beta[i])/gamma[i];
          } else if ( gamma[i] < -infinity ) {
            if ( L_bounded(IBi) ) QUOT = (L(IBi)-beta[i])/gamma[i];
          }
          if ( QUOT < beta_min ) {
            beta_min = QUOT; IBmin = IBi; P = i;
          }
          if ( abs(QUOT-beta_min) < eps && IBi < IBmin ) {
            IBmin = IBi; P = i;
          }
        }
      }
      /*
      // stopping criteria
      */
      real_type TQ{infinity};
      integer   INq{IN[Q]};
      integer   IBp{IB[P]};
      if ( U_bounded(INq) && L_bounded(INq) ) TQ = U(INq)-L(INq);

      UTILS_ASSERT0(
        TQ < infinity || beta_min < infinity,
        "StandardSolver. Problem is NOT solvable!\n"
      );
      if ( TQ < infinity && TQ <= beta_min ) {
        if      ( ceta[Q] > eps  ) x[INq] = U(INq);
        else if ( ceta[Q] < -eps ) x[INq] = L(INq);
        UTILS_ASSERT(
          ceta[Q] <= eps || ceta[Q] >= -eps,
          "StandardSolver. Could not find a new pivot row in iteration {}\n", iter
        );
      } else if ( beta_min < infinity && beta_min < TQ ) {
        if ( ceta[Q] > 0 ) {
          if      ( gamma[P] < 0 ) x[IBp] = U(IBp);
          else if ( gamma[P] > 0 ) x[IBp] = L(IBp);
        } else if ( ceta[Q] < 0 ) {
          if      ( gamma[P] < 0 ) x[IBp] = L(IBp);
          else if ( gamma[P] > 0 ) x[IBp] = U(IBp);
        }
        if ( pStream != nullptr )
          fmt::print( *pStream,
            "Pivot element gamma({},{}) = {}\n",
            IB[P], IN[Q], gamma[P]
          );
        /*
        // perform basis change
        */
        std::swap( IB[P], IN[Q] );
      } else {
        UTILS_ERROR0( "StandardSolver. No continuation defined!" );
      }
    }
  }

}
