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

#ifndef SIMPLEX_HH
#define SIMPLEX_HH

#include "Alglin_Config.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#ifndef SIMPLEX_API_DLL
  #ifdef ALGLIN_OS_WINDOWS
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

#ifndef SIMPLEX_VIRTUAL
  #define SIMPLEX_VIRTUAL SIMPLEX_API_DLL virtual
#endif

#ifndef SIMPLEX_ERROR
  #define SIMPLEX_ERROR(MSG) { \
    std::ostringstream ost; ost << MSG << '\n'; \
    throw std::runtime_error(ost.str()); \
  }
#endif

#ifndef SIMPLEX_ASSERT
  #define SIMPLEX_ASSERT(COND,MSG) if ( !(COND) ) SIMPLEX_ERROR(MSG);
#endif

#ifdef __GNUC__ 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-vtables"
#pragma clang diagnostic ignored "-Wpadded"
#endif

//! namespace for nonlinear systems and nonlinearsolver
namespace Simplex {

  using alglin::ostream_type;
  using alglin::integer;
  using std::vector;

  typedef alglin::doublereal valueType; //!< double value

  extern valueType const epsilon; // machine epsilon
  extern valueType const relaxedEpsilon;
  extern valueType const infinity;

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

    SIMPLEX_VIRTUAL
    ~StandardProblemBase()
    {}

    SIMPLEX_VIRTUAL valueType get_b_max_abs() const ALGLIN_PURE_VIRTUAL;
    SIMPLEX_VIRTUAL valueType get_c_max_abs() const ALGLIN_PURE_VIRTUAL;
    SIMPLEX_VIRTUAL valueType get_A_max_abs() const ALGLIN_PURE_VIRTUAL;

    SIMPLEX_VIRTUAL integer dim_x() const ALGLIN_PURE_VIRTUAL; //!< dimension of x
    SIMPLEX_VIRTUAL integer dim_g() const ALGLIN_PURE_VIRTUAL; //!< number of equality constraints

    SIMPLEX_VIRTUAL void load_c( valueType c[] ) const ALGLIN_PURE_VIRTUAL;

    //! fill the vector b with the rhs of the constraints \f$ Ax = b \$f
    SIMPLEX_VIRTUAL void load_b( valueType b[] ) const ALGLIN_PURE_VIRTUAL;

    /*!
      fill the sparse vector with the column `j_col` of matrix `A` in a sparse form as values
      a index of the row.
      \param  j_col  colum to be extracted
      \param  values the nonzeros elements of the column
      \param  i_row  the index of the row of the corresponding nonzeros element
      \return number of nonzeros elements of the column
    \*/
    SIMPLEX_VIRTUAL integer load_A_column( integer j_col, valueType values[], integer i_row[] ) const ALGLIN_PURE_VIRTUAL;

    //! subtract to `res` the product `Ax`
    SIMPLEX_VIRTUAL void subtract_Ax( valueType const x[], valueType res[] ) const ALGLIN_PURE_VIRTUAL;

    //! lower bound of `x_i`
    SIMPLEX_VIRTUAL valueType Lower( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! upper bound of `x_i`
    SIMPLEX_VIRTUAL valueType Upper( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! return true if lower bound of `x_i` is unlimited
    SIMPLEX_VIRTUAL bool Lower_is_free( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! return true if upper bound of `x_i` is unlimited
    SIMPLEX_VIRTUAL bool Upper_is_free( integer i ) const ALGLIN_PURE_VIRTUAL;

    SIMPLEX_API_DLL
    void
    info( ostream_type & stream ) {
      fmt::print(
        stream, "{:14} {:14} {:14} {:14}\n",
        "Flag", "Lower", "Upper", "Flag"
      );
      for ( integer i = 0; i < dim_x(); ++i )
        fmt::print(
          stream, "{:14} {:14.5} {:14.5} {:14}\n",
          (Lower_is_free(i)?"Free":"Bounded"),
          Lower(i), Upper(i),
          (Upper_is_free(i)?"Free":"Bounded")
        );
    }
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

    SIMPLEX_VIRTUAL
    ~ProblemBase()
    {}

    SIMPLEX_VIRTUAL valueType get_c_max_abs() const ALGLIN_PURE_VIRTUAL;
    SIMPLEX_VIRTUAL valueType get_A_max_abs() const ALGLIN_PURE_VIRTUAL;

    SIMPLEX_VIRTUAL integer dim_x() const ALGLIN_PURE_VIRTUAL; //!< dimension of x
    SIMPLEX_VIRTUAL integer dim_g() const ALGLIN_PURE_VIRTUAL; //!< number of equality constraints

    SIMPLEX_VIRTUAL void load_c( valueType c[] ) const ALGLIN_PURE_VIRTUAL;

    /*!
      fill the sparse vector with the column `j_col` of matrix `A` in a sparse form as values
      a index of the row.
      \param  j_col  colum to be extracted
      \param  values the nonzeros elements of the column
      \param  i_row  the index of the row of the corresponding nonzeros element
      \return number of nonzeros elements of the column
    \*/
    SIMPLEX_VIRTUAL integer load_A_column( integer j_col, valueType values[], integer i_row[] ) const ALGLIN_PURE_VIRTUAL;

    //! subtract to `res` the product `Ax`
    SIMPLEX_VIRTUAL void subtract_Ax( valueType const x[], valueType res[] ) const ALGLIN_PURE_VIRTUAL;

    //! lower bound of `x_i`
    SIMPLEX_VIRTUAL valueType Lower( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! upper bound of `x_i`
    SIMPLEX_VIRTUAL valueType Upper( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! return true if lower bound of `x_i` is unlimited
    SIMPLEX_VIRTUAL bool Lower_is_free( integer i ) const ALGLIN_PURE_VIRTUAL;

    //! return true if upper bound of `x_i` is unlimited
    SIMPLEX_VIRTUAL bool Upper_is_free( integer i ) const ALGLIN_PURE_VIRTUAL;

    SIMPLEX_API_DLL
    void
    info( ostream_type & stream ) {
      fmt::print(
        stream, "{:14} {:14} {:14} {:14}\n",
        "Flag", "Lower", "Upper", "Flag"
      );
      for ( integer i = 0; i < dim_x()+dim_x(); ++i )
        fmt::print(
          stream, "{:14} {:14.5} {:14.5} {:14}\n",
          (Lower_is_free(i)?"Free":"Bounded"),
          Lower(i), Upper(i),
          (Upper_is_free(i)?"Free":"Bounded")
        );
    }

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

    ProblemBase * problem;

  public:

    SIMPLEX_API_DLL
    explicit
    StandardProblemAdaptor( ProblemBase & problem_reference )
    : StandardProblemBase()
    , problem(&problem_reference)
    { }

    SIMPLEX_VIRTUAL
    ~StandardProblemAdaptor()
    {}

    SIMPLEX_API_DLL valueType get_b_max_abs() const { return 0; }
    SIMPLEX_API_DLL valueType get_c_max_abs() const { return problem->get_c_max_abs(); }
    SIMPLEX_API_DLL valueType get_A_max_abs() const { return problem->get_A_max_abs(); }

    SIMPLEX_VIRTUAL integer dim_x()  const { return problem->dim_x()+problem->dim_g(); }
    SIMPLEX_VIRTUAL integer dim_g()  const { return problem->dim_g(); }

    SIMPLEX_VIRTUAL
    void
    load_c( valueType c[] ) const {
      problem->load_c( c );
      alglin::zero( problem->dim_g(), c + problem->dim_x(), 1 );
    }

    SIMPLEX_VIRTUAL
    void
    load_b( valueType b[] ) const {
      alglin::zero( problem->dim_g(), b, 1 );
    }

    SIMPLEX_VIRTUAL
    integer
    load_A_column( integer j_col, valueType values[], integer i_row[] ) const {
      if ( j_col < problem->dim_x() )
        return problem->load_A_column( j_col, values, i_row );
      values[0] = -1;
      i_row[0]  = j_col - problem->dim_x();
      return 1;
    }

    SIMPLEX_VIRTUAL
    void
    subtract_Ax( valueType const x[], valueType res[] ) const {
      // [ A - I ] x
      problem->subtract_Ax( x, res );
      alglin::axpy( problem->dim_g(), 1.0, x + problem->dim_x(), 1, res, 1 );
    }

    SIMPLEX_VIRTUAL valueType Lower( integer i )         const { return problem->Lower(i); }
    SIMPLEX_VIRTUAL valueType Upper( integer i )         const { return problem->Upper(i); }
    SIMPLEX_VIRTUAL bool      Lower_is_free( integer i ) const { return problem->Lower_is_free(i); }
    SIMPLEX_VIRTUAL bool      Upper_is_free( integer i ) const { return problem->Upper_is_free(i); }
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

    StandardProblemBase *pBase;

    alglin::Malloc<valueType> baseReals;
    alglin::Malloc<integer>   baseInteger;

    valueType * d;
    valueType * values;

    integer n;
    integer m;
    integer nz;
    integer nw;
    integer np;
    integer *map_z;
    integer *map_w;
    integer *map_p;
    integer *map_case;
    integer *i_row;

    valueType b_max_abs;
    valueType c_max_abs;
    valueType A_max_abs;

  public:

    SIMPLEX_API_DLL
    explicit
    AuxProblem()
    : StandardProblemBase()
    , baseReals("Simplex::AuxProblem_reals")
    , baseInteger("Simplex::AuxProblem_integers")
    {}

    SIMPLEX_VIRTUAL
    ~AuxProblem()
    {}

    SIMPLEX_API_DLL
    void
    setup( StandardProblemBase * _pBase );

    SIMPLEX_VIRTUAL integer dim_x() const { return nz+nw+np+m; }
    SIMPLEX_VIRTUAL integer dim_g() const { return m; }

    SIMPLEX_API_DLL valueType get_b_max_abs() const { return b_max_abs; }
    SIMPLEX_API_DLL valueType get_c_max_abs() const { return c_max_abs; }
    SIMPLEX_API_DLL valueType get_A_max_abs() const { return A_max_abs; }

    SIMPLEX_VIRTUAL
    void
    load_c( valueType c[] ) const {
      integer nn = nz+nw+np;
      alglin::zero( nn, c, 1 );
      alglin::fill( m, c + nn, 1, 1.0 );
    }

    SIMPLEX_VIRTUAL
    void
    load_b( valueType b[] ) const {
      for ( integer i = 0; i < m; ++i )
        b[i] = std::abs(d[i]);
    }

    SIMPLEX_VIRTUAL integer load_A_column( integer j_col, valueType values[], integer i_row[] ) const;

    SIMPLEX_VIRTUAL valueType Lower( integer ) const;
    SIMPLEX_VIRTUAL valueType Upper( integer ) const;

    SIMPLEX_VIRTUAL bool Lower_is_free( integer ) const;
    SIMPLEX_VIRTUAL bool Upper_is_free( integer ) const;

    SIMPLEX_VIRTUAL void subtract_Ax( valueType const x[], valueType res[] ) const;

    //! get initial feasible point for the solution of Simplex problem
    SIMPLEX_API_DLL
    void
    feasible_point( valueType x[], integer IB[] ) const;

    //! get the solution of the Aux problem and transform to initial point of primal problem
    SIMPLEX_API_DLL
    void
    to_primal( valueType const x[], valueType xo[], integer IBo[] ) const;
  };

  /*\
   |   ___ _                _             _   ___         _    _
   |  / __| |_ __ _ _ _  __| |__ _ _ _ __| | | _ \_ _ ___| |__| |___ _ __
   |  \__ \  _/ _` | ' \/ _` / _` | '_/ _` | |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |___/\__\__,_|_||_\__,_\__,_|_| \__,_| |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  class StandardProblem : public StandardProblemBase {
  private:
    integer           n;
    integer           m;
    valueType const * c;
    valueType const * A;
    valueType const * b;
    valueType const * L;
    valueType const * U;
    integer           ldA;

    valueType b_max_abs;
    valueType c_max_abs;
    valueType A_max_abs;

    std::vector<bool> L_free, U_free;

  public:

    SIMPLEX_API_DLL
    explicit
    StandardProblem()
    : StandardProblemBase()
    , n(0)
    , m(0)
    , c(nullptr)
    , A(nullptr)
    , b(nullptr)
    , L(nullptr)
    , U(nullptr)
    , ldA(0)
    { }

    SIMPLEX_VIRTUAL
    ~StandardProblem()
    {}

    /*!
      Setup linear programs of type

      Minimise \f$ c'x \f$

      subject to
      \f[ l \leq x \leq u \f]
      \f[ Ax = b \f]

      where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
      \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
      The initial \f$ x \f$ is supposed to be a feasible basic solution.

      \param m   Number of rows of A (number of linear constraints)
      \param n   Number of optimization variables (dimension of x)
      \param A   The matrix A stored columnwise (Fortran storage)
      \param ldA Leading dimension of A (size m x n)
      \param b   r.h.s of equality constraints (size m)
      \param c   vector of objective function (size n)
      \param L   lower bound of x
      \param U   upper bound of x
    \*/

    SIMPLEX_API_DLL
    void
    setup(
      integer         m,
      integer         n,
      valueType const A[],
      integer         ldA,
      valueType const b[],
      valueType const c[],
      valueType const L[],
      valueType const U[]
    );

    SIMPLEX_API_DLL valueType get_b_max_abs() const { return b_max_abs; }
    SIMPLEX_API_DLL valueType get_c_max_abs() const { return c_max_abs; }
    SIMPLEX_API_DLL valueType get_A_max_abs() const { return A_max_abs; }

    SIMPLEX_VIRTUAL integer dim_x() const { return n; }
    SIMPLEX_VIRTUAL integer dim_g() const { return m; }

    SIMPLEX_VIRTUAL
    void
    load_c( valueType _c[] ) const {
      alglin::copy( n, c, 1, _c, 1 );
    }

    SIMPLEX_VIRTUAL
    void
    load_b( valueType _b[] ) const {
      alglin::copy( m, b, 1, _b, 1 );
    }

    SIMPLEX_VIRTUAL
    integer
    load_A_column( integer j_col, valueType values[], integer i_row[] ) const {
      alglin::copy( m, A+j_col*ldA, 1, values, 1 );
      for ( integer i = 0; i < m; ++i ) i_row[i] = i;
      return m;
    }

    //! subtract to `res` the product `Ax`
    SIMPLEX_VIRTUAL
    void
    subtract_Ax( valueType const x[], valueType res[] ) const {
      alglin::gemv( alglin::NO_TRANSPOSE, m, n, -1.0, A, ldA, x, 1, 1.0, res, 1);
    }

    SIMPLEX_VIRTUAL valueType Lower( integer i ) const { return L[i]; }
    SIMPLEX_VIRTUAL valueType Upper( integer i ) const { return U[i]; }

    SIMPLEX_VIRTUAL bool Lower_is_free( integer i ) const { return L_free[i]; }
    SIMPLEX_VIRTUAL bool Upper_is_free( integer i ) const { return U_free[i]; }

  };

  /*\
   |   ___         _    _
   |  | _ \_ _ ___| |__| |___ _ __
   |  |  _/ '_/ _ \ '_ \ / -_) '  \
   |  |_| |_| \___/_.__/_\___|_|_|_|
  \*/
  class Problem : public ProblemBase {
  private:
    integer           n;
    integer           m;
    valueType const * c;
    valueType const * A;
    valueType const * L;
    valueType const * U;
    integer           ldA;

    valueType c_max_abs;
    valueType A_max_abs;

    std::vector<bool> L_free, U_free;

  public:

    SIMPLEX_API_DLL
    explicit
    Problem()
    : ProblemBase()
    , n(0)
    , m(0)
    , c(nullptr)
    , A(nullptr)
    , L(nullptr)
    , U(nullptr)
    , ldA(0)
    { }

    SIMPLEX_VIRTUAL
    ~Problem()
    {}

    /*!
      Setup linear programs of type

      Minimise \f$ c'x \f$

      subject to
      \f[ l \leq x \leq u \f]
      \f[ Ax = b \f]

      where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
      \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
      The initial \f$ x \f$ is supposed to be a feasible basic solution.

      \param m   Number of rows of A (number of linear constraints)
      \param n   Number of optimization variables (dimension of x)
      \param A   The matrix A stored columnwise (Fortran storage)
      \param ldA Leading dimension of A (size m x n)
      \param c   vector of objective function (size n)
      \param L   lower bound of x
      \param U   upper bound of x
    \*/

    SIMPLEX_API_DLL
    void
    setup(
      integer         m,
      integer         n,
      valueType const A[],
      integer         ldA,
      valueType const c[],
      valueType const L[],
      valueType const U[]
    );

    SIMPLEX_API_DLL valueType get_b_max_abs() const { return 0; }
    SIMPLEX_API_DLL valueType get_c_max_abs() const { return c_max_abs; }
    SIMPLEX_API_DLL valueType get_A_max_abs() const { return A_max_abs; }

    SIMPLEX_VIRTUAL integer dim_x() const { return n; }
    SIMPLEX_VIRTUAL integer dim_g() const { return m; }

    SIMPLEX_VIRTUAL
    void
    load_c( valueType _c[] ) const {
      alglin::copy( n, c, 1, _c, 1 );
    }

    SIMPLEX_VIRTUAL
    integer
    load_A_column( integer j_col, valueType values[], integer i_row[] ) const {
      alglin::copy( m, A+j_col*ldA, 1, values, 1 );
      for ( integer i = 0; i < m; ++i ) i_row[i] = i;
      return m;
    }

    //! subtract to `res` the product `Ax`
    SIMPLEX_VIRTUAL
    void
    subtract_Ax( valueType const x[], valueType res[] ) const {
      alglin::gemv( alglin::NO_TRANSPOSE, m, n, -1.0, A, ldA, x, 1, 1.0, res, 1);
    }

    SIMPLEX_VIRTUAL valueType Lower( integer i ) const { return L[i]; }
    SIMPLEX_VIRTUAL valueType Upper( integer i ) const { return U[i]; }

    SIMPLEX_VIRTUAL bool Lower_is_free( integer i ) const { return L_free[i]; }
    SIMPLEX_VIRTUAL bool Upper_is_free( integer i ) const { return U_free[i]; }

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
  private:
    std::string const _name; //!< name of the NLP problem defined

    // block copy constructor
    StandardSolver();
    StandardSolver(StandardSolver const &);
    StandardSolver const & operator = (StandardSolver const &);

    alglin::Malloc<valueType> baseReals;
    alglin::Malloc<integer>   baseIntegers;

    ostream_type        * pStream;
    StandardProblemBase * problem;

    valueType L( integer i ) const { return problem->Lower(i); }
    valueType U( integer i ) const { return problem->Upper(i); }
    bool L_free( integer i ) const { return problem->Lower_is_free(i); }
    bool U_free( integer i ) const { return problem->Upper_is_free(i); }
    bool L_bounded( integer i ) const { return !problem->Lower_is_free(i); }
    bool U_bounded( integer i ) const { return !problem->Upper_is_free(i); }

    std::string Lstring( integer i ) const;
    std::string Ustring( integer i ) const;

  protected:

    integer m; //!< number of equality constraints
    integer n; //!< dimension of x (dimension of the solution)
    integer maxIter;

  public:

    SIMPLEX_API_DLL
    explicit
    StandardSolver( std::string const & n );

    SIMPLEX_VIRTUAL
    ~StandardSolver()
    {}

    //! The name of the class
    SIMPLEX_API_DLL std::string const & name(void) const { return _name; }

    /*!
      Solve linear programs of type

      Minimise \f$ c'x \f$

      subject to
      \f[ l \leq x \leq u \f]
      \f[ Ax = b \f]

      where \f$ A \f$ is a \f$ m \f$ by \f$ n\f$  matrix,
      \f$ \textrm{rank}(A)=m \f$ , \f$ m \leq n\f$ , \f$ l \leq u\f$ .
      The initial \f$ x \f$ is supposed to be a feasible basic solution.

      \param _problem Pointer to a class instance describing the problem
      \param x   in input feasible starting point for simplex algorithm.
                 On outpout the computed solution
      \param IB  index set defining the selection of basic variables
      \param eps
    \*/
    void
    solve(
      StandardProblemBase * _problem,
      valueType             x[],
      integer               IB[],
      valueType             eps = relaxedEpsilon
    );

  };

}

#ifdef __GNUC__ 
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
