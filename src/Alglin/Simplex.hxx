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

///
/// file: Simplex.hxx
///

#ifndef SIMPLEX_API_DLL
  #ifdef UTILS_OS_WINDOWS
    #ifdef SIMPLEX_EXPORT
      #define SIMPLEX_API_DLL __declspec(dllexport)
    #elif defined(SIMPLEX_IMPORT)
      #define SIMPLEX_API_DLL __declspec(dllimport)
    #else
      #define SIMPLEX_API_DLL
    #endif
  #else
    #define SIMPLEX_API_DLL
  #endif
#endif

//!
//! Namespace for nonlinear systems and nonlinearsolver
//!
namespace Simplex {

  using alglin::ostream_type;
  using alglin::integer;
  using alglin::real_type; //!< double value

  using std::abs;
  using std::fill;
  using std::fill_n;
  using std::copy_n;
  using std::vector;
  using std::string;
  using std::string_view;

  extern real_type const epsilon; // machine epsilon
  extern real_type const relaxedEpsilon;
  extern real_type const infinity;

  /*\
   |   ___ _                _             _   ___         _    _             ___
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __   | _ ) __ _ ___ ___
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \  | _ \/ _` (_-</ -_)
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_| |___/\__,_/__/\___|
  \*/
  /*--------------------------------------------------------------------------*\
   !                                                                          !
   !     Base class for the definition of Primal simplex method               !
   !     for solving linear programs of type                                  !
   !                                                                          !
   !     Minimise     c'x                                                     !
   !                                                                          !
   !     subject to   l <= x <= u                                             !
   !                  Ax = b                                                  !
   !                                                                          !
   !     where A is a m by n matrix, rank(A)=m, m<=n, l<=u.                   !
   !                                                                          !
  \*--------------------------------------------------------------------------*/
  class StandardProblemBase {
  private:

    // block copy constructor
    StandardProblemBase(StandardProblemBase const &);
    StandardProblemBase const & operator = (StandardProblemBase const &);

  public:

    SIMPLEX_API_DLL
    explicit
    StandardProblemBase()
    {}

    virtual ~StandardProblemBase() = default;

    virtual real_type get_b_max_abs() const = 0;
    virtual real_type get_c_max_abs() const = 0;
    virtual real_type get_A_max_abs() const = 0;

    virtual integer dim_x() const = 0; //!< dimension of x
    virtual integer dim_g() const = 0; //!< number of equality constraints

    virtual void load_c( real_type c[] ) const = 0;

    //!
    //! Fill the vector b with the rhs of the constraints \f$ Ax = b \$f
    //!
    virtual void load_b( real_type b[] ) const = 0;

    //!
    //! Fill the sparse vector with the column `j_col` of matrix `A` in a sparse form as values
    //! a index of the row.
    //!
    //! \param  j_col  colum to be extracted
    //! \param  values the nonzeros elements of the column
    //! \param  i_row  the index of the row of the corresponding nonzeros element
    //! \return number of nonzeros elements of the column
    //!
    virtual integer load_A_column( integer j_col, real_type values[], integer i_row[] ) const = 0;

    //!
    //! Subtract to `res` the product `Ax`.
    //!
    virtual void subtract_Ax( real_type const x[], real_type res[] ) const = 0;

    //!
    //! Lower bound of `x_i`.
    //!
    virtual real_type Lower( integer i ) const = 0;

    //!
    //! Upper bound of `x_i`.
    //!
    virtual real_type Upper( integer i ) const = 0;

    //!
    //! Return true if lower bound of `x_i` is unlimited.
    //!
    virtual bool Lower_is_free( integer i ) const = 0;

    //!
    //! Return true if upper bound of `x_i` is unlimited.
    //!
    virtual bool Upper_is_free( integer i ) const = 0;

    SIMPLEX_API_DLL
    string
    info() const {
      string res{ fmt::format( "{:14} {:14} {:14} {:14}\n", "Flag", "Lower", "Upper", "Flag" ) };
      for ( integer i{0}; i < dim_x(); ++i ) {
        res += fmt::format(
          "{:14} {:14.5} {:14.5} {:14}\n",
          (Lower_is_free(i)?"Free":"Bounded"),
          Lower(i), Upper(i),
          (Upper_is_free(i)?"Free":"Bounded")
        );
      }
      return res;
    }

    SIMPLEX_API_DLL
    void info( ostream_type & stream ) const { stream << info(); }

  };

  /*\
   |   ___         _    _             ___
   |  | _ \_ _ ___| |__| |___ _ __   | _ ) __ _ ___ ___
   |  |  _/ '_/ _ \ '_ \ / -_) '  \  | _ \/ _` (_-</ -_)
   |  |_| |_| \___/_.__/_\___|_|_|_| |___/\__,_/__/\___|
  \*/
  /*--------------------------------------------------------------------------*\
   !                                                                          !
   !     Base class for the definition of Primal simplex method               !
   !     for solving linear programs of type                                  !
   !                                                                          !
   !     Minimise     c'x                                                     !
   !                                                                          !
   !     subject to   l <= / x  \ <= u                                        !
   !                       \ Ax /                                             !
   !                                                                          !
   !     where A is a m by n matrix, rank(A)=m, m<=n, l<=u.                   !
   !                                                                          !
  \*--------------------------------------------------------------------------*/
  class ProblemBase {
  private:

    // block copy constructor
    ProblemBase(ProblemBase const &);
    ProblemBase const & operator = (ProblemBase const &);

  public:

    SIMPLEX_API_DLL
    explicit
    ProblemBase()
    {}

    ~ProblemBase()
    {}

    virtual real_type get_b_max_abs() const = 0;
    virtual real_type get_c_max_abs() const = 0;
    virtual real_type get_A_max_abs() const = 0;

    virtual integer dim_x() const = 0; //!< dimension of x
    virtual integer dim_g() const = 0; //!< number of equality constraints

    virtual void load_c( real_type c[] ) const = 0;

    //!
    //! Fill the sparse vector with the column `j_col` of matrix `A` in a sparse form as values
    //! a index of the row.
    //!
    //! \param  j_col  colum to be extracted
    //! \param  values the nonzeros elements of the column
    //! \param  i_row  the index of the row of the corresponding nonzeros element
    //! \return number of nonzeros elements of the column
    //!
    virtual integer load_A_column( integer j_col, real_type values[], integer i_row[] ) const = 0;

    //!
    //! Subtract to `res` the product `Ax`.
    //!
    virtual void subtract_Ax( real_type const x[], real_type res[] ) const = 0;

    //!
    //! Lower bound of `x_i`.
    //!
    virtual real_type Lower( integer i ) const = 0;

    //!
    //! Upper bound of `x_i`.
    //!
    virtual real_type Upper( integer i ) const = 0;

    //!
    //! Return true if lower bound of `x_i` is unlimited.
    //!
    virtual bool Lower_is_free( integer i ) const = 0;

    //!
    //! Return true if upper bound of `x_i` is unlimited.
    //!
    virtual bool Upper_is_free( integer i ) const = 0;

    SIMPLEX_API_DLL
    string
    info() const {
      string res{ fmt::format( "{:14} {:14} {:14} {:14}\n", "Flag", "Lower", "Upper", "Flag" ) };
      for ( integer i{0}; i < dim_x()+dim_x(); ++i ) {
        res += fmt::format(
          "{:14} {:14.5} {:14.5} {:14}\n",
          (Lower_is_free(i)?"Free":"Bounded"),
          Lower(i), Upper(i),
          (Upper_is_free(i)?"Free":"Bounded")
        );
      }
      return res;
    }

    SIMPLEX_API_DLL
    void info( ostream_type & stream ) const { stream << info(); }

  };

  /*\
   |   ___ _                _             _   ___         _    _               _      _           _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __     /_\  __| |__ _ _ __| |_ ___ _ _
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \   / _ \/ _` / _` | '_ \  _/ _ \ '_|
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_| /_/ \_\__,_\__,_| .__/\__\___/_|
   |                                                                                       |_|
  \*/
  /*--------------------------------------------------------------------------*\
   !                                                                          !
   !     Base class for the definition of Primal simplex method               !
   !     for solving linear programs of type                                  !
   !                                                                          !
   !     Minimise     c'x                                                     !
   !                                                                          !
   !     subject to   l <= / x  \ <= u                                        !
   !                       \ Ax /                                             !
   !                                                                          !
   !     where A is a m by n matrix, rank(A)=m, m<=n, l<=u.                   !
   !     The initial x is supposed to be a feasible basic solution.           !
   !                                                                          !
   !     Tranformed to Standard Problem                                       !
   !                                                                          !
   !     Minimise     c'x                                                     !
   !                                                                          !
   !     subject to   l <= / x  \ <= u                                        !
   !                       \ z  /                                             !
   !                                                                          !
   !                       [A -I] / x \ = 0                                   !
   !                              \ z /                                       !
  \*--------------------------------------------------------------------------*/
  class StandardProblemAdaptor : public StandardProblemBase {
  private:

    ProblemBase * m_problem{nullptr};

  public:

    SIMPLEX_API_DLL
    explicit
    StandardProblemAdaptor( ProblemBase & problem_reference )
    : StandardProblemBase()
    , m_problem(&problem_reference)
    { }

    real_type get_b_max_abs() const override { return 0; }
    real_type get_c_max_abs() const override { return m_problem->get_c_max_abs(); }
    real_type get_A_max_abs() const override { return m_problem->get_A_max_abs(); }

    integer dim_x() const override { return m_problem->dim_x()+m_problem->dim_g(); }
    integer dim_g() const override { return m_problem->dim_g(); }

    void
    load_c( real_type c[] ) const override {
      m_problem->load_c( c );
      fill_n( c + m_problem->dim_x(), m_problem->dim_g(), 0 );
    }

    void
    load_b( real_type b[] ) const override {
      fill_n( b, m_problem->dim_g(), 0 );
    }

    integer
    load_A_column(
      integer   j_col,
      real_type values[],
      integer   i_row[]
    ) const override {
      if ( j_col < m_problem->dim_x() )
        return m_problem->load_A_column( j_col, values, i_row );
      values[0] = -1;
      i_row[0]  = j_col - m_problem->dim_x();
      return 1;
    }

    void
    subtract_Ax( real_type const x[], real_type res[] ) const override {
      // [ A - I ] x
      m_problem->subtract_Ax( x, res );
      alglin::axpy( m_problem->dim_g(), 1.0, x + m_problem->dim_x(), 1, res, 1 );
    }

    real_type Lower( integer i )         const override { return m_problem->Lower(i); }
    real_type Upper( integer i )         const override { return m_problem->Upper(i); }
    bool      Lower_is_free( integer i ) const override { return m_problem->Lower_is_free(i); }
    bool      Upper_is_free( integer i ) const override { return m_problem->Upper_is_free(i); }
  };

  /*\
   |     _             ___         _    _
   |    /_\ _  ___ __ | _ \_ _ ___| |__| |___ _ __
   |   / _ \ || \ \ / |  _/ '_/ _ \ '_ \ / -_) '  \
   |  /_/ \_\_,_/_\_\ |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  /*--------------------------------------------------------------------------*\
   !                                                                          !
   !     Base class for the definition of Primal simplex method               !
   !     for solving linear programs of type                                  !
   !                                                                          !
   !     Minimise     e's                                                     !
   !                                                                          !
   !     subject to  ( Z W P I ) / z \                  / z \                 !
   !                             | w | = Dd >= 0   0 <= | w | <= infinity     !
   !                             | p |                  \ s /                 !
   !                             \ s /       -infinity <= p <= infinity       !
   !                                                                          !
   !     Feasible initial point z = w = q = 0, s = Dd                         !
   !                                                                          !
   !     Connection with the initial problem                                  !
   !                                                                          !
   !     z_i = x_i - l_i  if l_i > -infinity                                  !
   !                                                                          !
   !     w_i = u_i - x_i  if u_i < infinity                                   !
   !                                                                          !
   !     p_i = x_i        if l_i = -infinity   and   u_i = infinity           !
   !                                                                          !
   !     recontruction of x from the solution of aux problem                  !
   !                                                                          !
   !           +                                                              !
   !           | (z_i-w_i+(l_i+u_i))/2 if -infinity < l_i and u_i < infinity  !
   !           |                                                              !
   !           | p_i                   if |u_i| = |l_i| = infinity            !
   !     x_i = |                                                              !
   !           | u_i-w_i               if l_i = -infinity                     !
   !           |                                                              !
   !           | z_i+l_i               if u_i = infinity                      !
   !           +                                                              !
   !                                                                          !
   !     d = b-A*q                                                            !
   !                                                                          !
   !           +                                                              !
   !           | (l_i+u_i)/2   if -infinity < l_i and u_i < infinity          !
   !           |                                                              !
   !           | 0             if |u_i| = |l_i| = infinity                    !
   !     q_i = |                                                              !
   !           | u_i           if l_i = -infinity                             !
   !           |                                                              !
   !           | l_i           if u_i = infinity                              !
   !           +                                                              !
   !                                                                          !
  \*--------------------------------------------------------------------------*/
  class AuxProblem : public StandardProblemBase {
  private:

    StandardProblemBase * m_problem_base{nullptr};

    alglin::Malloc<real_type> m_mem{"Simplex::AuxProblem::m_mem"};
    alglin::Malloc<integer>   m_mem_int{"Simplex::AuxProblem::m_mem_int"};

    real_type * m_d{nullptr};
    real_type * m_values{nullptr};

    integer m_n{0};
    integer m_m{0};
    integer m_nz{0};
    integer m_nw{0};
    integer m_np{0};
    integer * m_map_z{nullptr};
    integer * m_map_w{nullptr};
    integer * m_map_p{nullptr};
    integer * m_map_case{nullptr};
    integer * m_i_row{nullptr};

    real_type m_b_max_abs{0};
    real_type m_c_max_abs{0};
    real_type m_A_max_abs{0};

  public:

    SIMPLEX_API_DLL
    explicit
    AuxProblem()
    : StandardProblemBase()
    {}

    ~AuxProblem() override {}

    SIMPLEX_API_DLL
    void
    setup( StandardProblemBase * _pBase );

    integer dim_x() const override { return m_nz+m_nw+m_np+m_m; }
    integer dim_g() const override { return m_m; }

    real_type get_b_max_abs() const override { return m_b_max_abs; }
    real_type get_c_max_abs() const override { return m_c_max_abs; }
    real_type get_A_max_abs() const override { return m_A_max_abs; }

    void
    load_c( real_type c[] ) const override {
      integer nn{ m_nz+m_nw+m_np };
      fill_n( c, nn, 0 );
      fill_n( c + nn, m_m, 1 );
    }

    void
    load_b( real_type b[] ) const override {
      for ( integer i{0}; i < m_m; ++i )
        b[i] = abs(m_d[i]);
    }

    integer load_A_column( integer j_col, real_type values[], integer i_row[] ) const override;

    real_type Lower( integer ) const override;
    real_type Upper( integer ) const override;

    bool Lower_is_free( integer ) const override;
    bool Upper_is_free( integer ) const override;

    void subtract_Ax( real_type const x[], real_type res[] ) const override;

    //!
    //! Get initial feasible point for the solution of Simplex problem.
    //!
    SIMPLEX_API_DLL
    void
    feasible_point( real_type x[], integer IB[] ) const;

    //!
    //! Get the solution of the Aux problem and transform
    //! to initial point of primal problem,
    //!
    SIMPLEX_API_DLL
    void
    to_primal( real_type const x[], real_type xo[], integer IBo[] ) const;
  };

  /*\
   |   ___ _                _             _   ___         _    _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  class StandardProblem : public StandardProblemBase {
  private:
    integer           m_n{0};
    integer           m_m{0};
    real_type const * m_c{nullptr};
    real_type const * m_A{nullptr};
    real_type const * m_b{nullptr};
    real_type const * m_L{nullptr};
    real_type const * m_U{nullptr};
    integer           m_A_ldim{0};

    real_type m_b_max_abs{0};
    real_type m_c_max_abs{0};
    real_type m_A_max_abs{0};

    vector<bool> m_L_free, m_U_free;

  public:

    SIMPLEX_API_DLL
    explicit
    StandardProblem()
    : StandardProblemBase()
    { }

    ~StandardProblem() override
    {}

    //!
    //! Setup linear programs of type
    //!
    //! Minimise \f$ c'x \f$
    //!
    //! subject to
    //! \f[ l \leq x \leq u \f]
    //! \f[ Ax = b \f]
    //!
    //! where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
    //! \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
    //! The initial \f$ x \f$ is supposed to be a feasible basic solution.
    //!
    //! \param m   Number of rows of A (number of linear constraints)
    //! \param n   Number of optimization variables (dimension of x)
    //! \param A   The matrix A stored columnwise (Fortran storage)
    //! \param ldA Leading dimension of A (size m x n)
    //! \param b   r.h.s of equality constraints (size m)
    //! \param c   vector of objective function (size n)
    //! \param L   lower bound of x
    //! \param U   upper bound of x
    //!
    SIMPLEX_API_DLL
    void
    setup(
      integer         m,
      integer         n,
      real_type const A[],
      integer         ldA,
      real_type const b[],
      real_type const c[],
      real_type const L[],
      real_type const U[]
    );

    real_type get_b_max_abs() const override { return m_b_max_abs; }
    real_type get_c_max_abs() const override { return m_c_max_abs; }
    real_type get_A_max_abs() const override { return m_A_max_abs; }

    integer dim_x() const override { return m_n; }
    integer dim_g() const override { return m_m; }

    void load_c( real_type c[] ) const override { copy_n( m_c, m_n, c ); }
    void load_b( real_type b[] ) const override { copy_n( m_b, m_m, b ); }

    integer
    load_A_column(
      integer   j_col,
      real_type values[],
      integer   i_row[]
    ) const override {
      copy_n( m_A+j_col*m_A_ldim, m_m, values );
      for ( integer i{0}; i < m_m; ++i ) i_row[i] = i;
      return m_m;
    }

    //! subtract to `res` the product `Ax`
    void
    subtract_Ax( real_type const x[], real_type res[] ) const override {
      using namespace alglin;
      gemv( Transposition::NO, m_m, m_n, -1.0, m_A, m_A_ldim, x, 1, 1.0, res, 1);
    }

    real_type Lower( integer i ) const override { return m_L[i]; }
    real_type Upper( integer i ) const override { return m_U[i]; }

    bool Lower_is_free( integer i ) const override { return m_L_free[i]; }
    bool Upper_is_free( integer i ) const override { return m_U_free[i]; }

  };

  /*\
   |   ___         _    _
   |  | _ \_ _ ___| |__| |___ _ __
   |  |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  class Problem : public ProblemBase {
  private:
    integer           m_n{0};
    integer           m_m{0};
    real_type const * m_c{nullptr};
    real_type const * m_A{nullptr};
    real_type const * m_L{nullptr};
    real_type const * m_U{nullptr};
    integer           m_A_ldim{0};

    real_type m_c_max_abs{0};
    real_type m_A_max_abs{0};

    vector<bool> m_L_free, m_U_free;

  public:

    SIMPLEX_API_DLL
    explicit
    Problem()
    : ProblemBase()
    { }

    ~Problem()
    {}

    //!
    //! Setup linear programs of type
    //!
    //! Minimise \f$ c'x \f$
    //!
    //! subject to
    //! \f[ l \leq x \leq u \f]
    //! \f[ Ax = b \f]
    //!
    //! where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
    //! \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
    //! The initial \f$ x \f$ is supposed to be a feasible basic solution.
    //!
    //! \param m   Number of rows of A (number of linear constraints)
    //! \param n   Number of optimization variables (dimension of x)
    //! \param A   The matrix A stored columnwise (Fortran storage)
    //! \param ldA Leading dimension of A (size m x n)
    //! \param c   vector of objective function (size n)
    //! \param L   lower bound of x
    //! \param U   upper bound of x
    //!
    SIMPLEX_API_DLL
    void
    setup(
      integer         m,
      integer         n,
      real_type const A[],
      integer         ldA,
      real_type const c[],
      real_type const L[],
      real_type const U[]
    );

    real_type get_b_max_abs() const override { return 0; }
    real_type get_c_max_abs() const override { return m_c_max_abs; }
    real_type get_A_max_abs() const override { return m_A_max_abs; }

    integer dim_x() const override { return m_n; }
    integer dim_g() const override { return m_m; }

    void load_c( real_type c[] ) const override { copy_n( m_c, m_n, c ); }

    integer
    load_A_column(
      integer   j_col,
      real_type values[],
      integer   i_row[]
    ) const override {
      copy_n( m_A+j_col*m_A_ldim, m_m, values );
      for ( integer i{0}; i < m_m; ++i ) i_row[i] = i;
      return m_m;
    }

    //!
    //! Subtract to `res` the product `Ax`.
    //!
    void
    subtract_Ax( real_type const x[], real_type res[] ) const override {
      using namespace alglin;
      gemv( Transposition::NO, m_m, m_n, -1.0, m_A, m_A_ldim, x, 1, 1.0, res, 1);
    }

    real_type Lower( integer i ) const override { return m_L[i]; }
    real_type Upper( integer i ) const override { return m_U[i]; }

    bool Lower_is_free( integer i ) const override { return m_L_free[i]; }
    bool Upper_is_free( integer i ) const override { return m_U_free[i]; }

  };


  /*\
   |   ___ _                _             _   ___      _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | / __| ___| |_ _____ _ _
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | \__ \/ _ \ \ V / -_) '_|
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |___/\___/_|\_/\___|_|
  \*/
  /*--------------------------------------------------------------------------*\
   !                                                                          !
   !     Primal simplex method for solving linear programs of type            !
   !                                                                          !
   !     Minimise     c'x                                                     !
   !                                                                          !
   !     subject to   l <= x <= u                                             !
   !                  Ax = b                                                  !
   !                                                                          !
   !     where A is a m by n matrix, rank(A)=m, m<=n, l<=u.                   !
   !     The initial x is supposed to be a feasible basic solution.           !
   !                                                                          !
  \*--------------------------------------------------------------------------*/
  class StandardSolver {
  public:
    // block copy constructor
    StandardSolver() = delete;
    StandardSolver( StandardSolver const & ) = delete;
    StandardSolver const & operator = ( StandardSolver const & ) = delete;

  private:
    string const m_name; //!< name of the NLP problem defined

    alglin::Malloc<real_type> m_mem;
    alglin::Malloc<integer>   m_mem_int;

    ostream_type        * pStream{nullptr};
    StandardProblemBase * m_problem{nullptr};

    real_type L( integer i ) const { return m_problem->Lower(i); }
    real_type U( integer i ) const { return m_problem->Upper(i); }
    bool L_free( integer i ) const { return m_problem->Lower_is_free(i); }
    bool U_free( integer i ) const { return m_problem->Upper_is_free(i); }
    bool L_bounded( integer i ) const { return !m_problem->Lower_is_free(i); }
    bool U_bounded( integer i ) const { return !m_problem->Upper_is_free(i); }

    string Lstring( integer i ) const;
    string Ustring( integer i ) const;

  protected:

    integer m_m{0}; //!< number of equality constraints
    integer m_n{0}; //!< dimension of x (dimension of the solution)
    integer m_max_iter{1000};

  public:

    SIMPLEX_API_DLL
    explicit
    StandardSolver( string_view n );

    virtual
    ~StandardSolver() = default;

    //!
    //! The name of the class.
    //!
    SIMPLEX_API_DLL string_view name(void) const { return m_name; }

    //!
    //! Solve linear programs of type
    //!
    //! Minimise \f$ c'x \f$
    //!
    //! subject to
    //! \f[ l \leq x \leq u \f]
    //! \f[ Ax = b \f]
    //!
    //! where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
    //! \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
    //! The initial \f$ x \f$ is supposed to be a feasible basic solution.
    //!
    //! \param _problem Pointer to a class instance describing the problem
    //! \param x   in input feasible starting point for simplex algorithm.
    //!            On outpout the computed solution
    //! \param IB  index set defining the selection of basic variables
    //! \param eps
    //!
    void
    solve(
      StandardProblemBase * problem,
      real_type             x[],
      integer               IB[],
      real_type             eps = relaxedEpsilon
    );

  };

}

///
/// eof: Simplex.hxx
///
