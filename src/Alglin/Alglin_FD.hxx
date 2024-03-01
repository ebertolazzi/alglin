/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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

///
/// file: Alglin_FD.hxx
///

#include <tuple>
#include <string>
#include <cmath>

namespace alglin {

  using std::abs;
  using std::max;
  using std::min;
  using std::pow;
  //using std::cbrt;
  using std::string;
  using std::numeric_limits;
  using std::vector;
  using std::tuple;
  using std::tie;

  template <typename Number>
  class finite_difference_epsilon {
    Number const m_epsilon{numeric_limits<Number>::epsilon()};
    Number const m_epsilon1{sqrt(numeric_limits<Number>::epsilon())};
    Number const m_epsilon2{pow(numeric_limits<Number>::epsilon(),Number(0.75))};
    Number const m_epsilon3{pow(numeric_limits<Number>::epsilon(),Number(0.25))};
  public:
    Number epsilon1( Number v ) const { return (abs(v)+1)*m_epsilon1; }
    Number epsilon2( Number v ) const { return (abs(v)+1)*m_epsilon2; }
    Number epsilon3( Number v ) const { return (abs(v)+1)*m_epsilon3; }
  };

  template <typename Number>
  static
  Number
  finite_difference_side(
    Number f0,
    Number f1,
    Number f2,
    Number h1,
    Number h2
  ) {
    Number df1{f1-f0};
    Number df2{f2-f0};
    Number dH{h1-h2};
    return ((h1/h2)*df2-(h2/h1)*df1)/dH;
  }

  template <typename Number>
  static
  bool
  finite_difference_side(
    Number const f0[],
    Number const f1[],
    Number const f2[],
    Number       h1,
    Number       h2,
    Number       Df[],
    integer      dim
  ) {
    bool ok{true};
    Number dH{h1-h2};
    Number t1{-(h2/h1)/dH};
    Number t2{(h1/h2)/dH};
    for ( integer i{0}; ok && i < dim; ++i ) {
      Number df1{f1[i]-f0[i]};
      Number df2{f2[i]-f0[i]};
      Df[i] = df1*t1+df2*t2;
      ok = Utils::is_finite(Df[i]);
    }
    return ok;
  }

  //
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //

  #define ALGLIN_FINITE_DIFFERENCE_WEIGHT                 \
    Number w1{abs(D1)+numeric_limits<Number>::epsilon()}; \
    Number w2{abs(D2)+numeric_limits<Number>::epsilon()}; \
    Number ww{max(w1,w2)};                                \
    w1 /= ww; w1 = sqrt(w1);                              \
    w2 /= ww; w2 = sqrt(w2)

  #define ALGLIN_FINITE_DIFFERENCE_AVERAGE (D1*w2+D2*w1)/(w1+w2)

  template <typename Number>
  static
  Number
  finite_difference_centered(
    Number fL,
    Number fC,
    Number fR,
    Number h
  ) {
    Number D1{(fR-fC)/h};
    Number D2{(fC-fL)/h};
    ALGLIN_FINITE_DIFFERENCE_WEIGHT;
    return ALGLIN_FINITE_DIFFERENCE_AVERAGE;
  }

  template <typename Number>
  static
  integer
  finite_difference_centered(
    Number   fL, bool ok_L,
    Number   fC, bool ok_C,
    Number   fR, bool ok_R,
    Number   h,
    Number & res
  ) {
    if ( ok_R && ok_L && ok_C ) {
      Number D1{(fR-fC)/h};
      Number D2{(fC-fL)/h};
      ALGLIN_FINITE_DIFFERENCE_WEIGHT;
      res = ALGLIN_FINITE_DIFFERENCE_AVERAGE;
      return 0;
    } else if ( ok_C && ok_R ) {
      res = (fR-fC)/h;
      return 1;
    } else if ( ok_C && ok_L ) {
      res = (fC-fL)/h;
      return -1;
    }
    res = 0;
    return -2;
  }

  template <typename Number>
  static
  integer
  finite_difference_centered(
    Number const fL[], bool ok_L,
    Number const fC[], bool ok_C,
    Number const fR[], bool ok_R,
    Number       h,
    Number       res[],
    integer      f_dim
  ) {
    if ( ok_R && ok_L && ok_C  ) {
      for ( integer i{0}; i < f_dim; ++i ) {
        Number D1{(fR[i]-fC[i])/h};
        Number D2{(fC[i]-fL[i])/h};
        ALGLIN_FINITE_DIFFERENCE_WEIGHT;
        res[i] = ALGLIN_FINITE_DIFFERENCE_AVERAGE;
      }
      return 0;
    } else if ( ok_C && ok_R ) {
      for ( integer i{0}; i < f_dim; ++i ) res[i] = (fR[i]-fC[i])/h;
      return 1;
    } else if ( ok_C && ok_L ) {
      for ( integer i{0}; i < f_dim; ++i ) res[i] = (fC[i]-fL[i])/h;
      return -1;
    }
    for ( integer i{0}; i < f_dim; ++i ) res[i] = 0;
    return -2;
  }

  /*
  //    ____               _ _            _
  //   / ___|_ __ __ _  __| (_) ___ _ __ | |_
  //  | |  _| '__/ _` |/ _` | |/ _ \ '_ \| __|
  //  | |_| | | | (_| | (_| | |  __/ | | | |_
  //   \____|_|  \__,_|\__,_|_|\___|_| |_|\__|
  */
  template <typename FUNCTION, typename Number>
  static
  bool
  finite_difference_gradient(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number       grad[]
  ) {

    static finite_difference_epsilon<Number> EPS;

    auto FUN = [&fun]( Number const x[], Number & f ) -> bool {
      bool ok{fun( x, f )};
      if ( ok ) ok = Utils::is_finite(f);
      return ok;
    };

    Number vC{0}, vR{0}, vL{0}; // only to stop warning
    bool ok_C = FUN( x, vC );
    if ( !ok_C ) return false;

    Number * X = const_cast<Number*>(x);

    for ( integer i{0}; i < dim_x; ++i ) {
      Number temp{x[i]};
      Number h1{EPS.epsilon1(temp)};
      Number h2{EPS.epsilon2(temp)};
      X[i] = temp+h1; bool ok_R{FUN( X, vR )};
      X[i] = temp-h1; bool ok_L{FUN( X, vL )};
      integer ic = finite_difference_centered( vL, ok_L, vC, ok_C, vR, ok_R, h1, grad[i] );
      switch ( ic ) {
      case  0:
        break;
      case  1:
        {
          Number & vRR = vL;
          X[i] = temp+h2; // modify the vector only at i position
          bool ok_RR{FUN( X, vRR )};
          if ( ok_RR ) grad[i] = finite_difference_side( vC, vR, vRR, h1, h2 );
          if ( ! (ok_RR&&Utils::is_finite(grad[i])) ) grad[i] = (vR-vC)/h1; // low precision FD
        }
        break;
      case -1:
        {
          Number & vLL = vR;
          X[i] = temp-h2; // modify the vector only at i position
          bool ok_LL{FUN( X, vLL )};
          if ( ok_LL ) grad[i] = finite_difference_side( vC, vL, vLL, -h1, -h2 );
          if ( ! (ok_LL&&Utils::is_finite(grad[i])) ) grad[i] = (vC-vL)/h1; // low precision FD
        }
        break;
      case -2:
        return false;
        break;
      }
      X[i] = temp; // restore i position
      if ( !Utils::is_finite(grad[i]) ) return false;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  static
  void
  finite_difference_check_gradient(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number const grad[],
    Number       epsi,
    //            i       A      FD     err   scale
    vector<tuple<integer,Number,Number,Number,Number>> & err_list
  ) {

    static finite_difference_epsilon<Number> EPS;

    auto FUN = [&fun]( Number const x[], Number & f ) -> bool {
      bool ok{fun( x, f )};
      if ( ok ) ok = Utils::is_finite(f);
      return ok;
    };

    Number vC{0}, vR{0}, vL{0}, gradi{0}; // only to stop warning
    bool ok_C{FUN( x, vC )};
    if ( !ok_C ) {
      err_list.emplace_back( -1, 0, 0, 0, 0 );
      return;
    }

    Number * X{const_cast<Number*>(x)};

    for ( integer i{0}; i < dim_x; ++i ) {
      Number temp{x[i]};
      Number h1{EPS.epsilon1(temp)};
      Number h2{EPS.epsilon2(temp)};
      X[i] = temp+h1; bool ok_R{FUN( X, vR )};
      X[i] = temp-h1; bool ok_L{FUN( X, vL )};
      integer ic = finite_difference_centered( vL, ok_L, vC, ok_C, vR, ok_R, h1, gradi );
      switch ( ic ) {
      case  0:
        break;
      case  1:
        {
          Number & vRR = vL;
          X[i] = temp+h2; // modify the vector only at i position
          bool ok_RR{FUN( X, vRR )};
          if ( ok_RR ) gradi = finite_difference_side( vC, vR, vRR, h1, h2 );
          if ( ! (ok_RR&&Utils::is_finite(gradi)) ) gradi = (vR-vC)/h1; // low precision FD
        }
        break;
      case -1:
        {
          Number & vLL = vR;
          X[i] = temp-h2; // modify the vector only at i position
          bool ok_LL{FUN( X, vLL )};
          if ( ok_LL ) gradi = finite_difference_side( vC, vL, vLL, -h1, -h2 );
          if ( ! (ok_LL&&Utils::is_finite(gradi)) ) gradi = (vC-vL)/h1; // low precision FD
        }
        break;
      case -2:
        break;
      }
      X[i] = temp; // restore i position
      Number scale = max( Number(1), max(abs(gradi), abs(grad[i])) );
      Number err   = abs(gradi - grad[i]);
      if ( !Utils::is_finite(gradi)   ||
           !Utils::is_finite(grad[i]) ||
           err > epsi * scale )
        err_list.emplace_back( i, grad[i], gradi, err, scale );
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  static
  string
  finite_difference_check_gradient(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number const grad[],
    Number       epsi
  ) {
    vector<tuple<integer,Number,Number,Number,Number>> err_list;
    finite_difference_check_gradient( x, dim_x, fun, grad, epsi, err_list );
    string res{""};
    for ( auto & e : err_list ) {
      integer i;
      Number A, FD, err, scale;
      tie ( i, A, FD, err, scale ) = e;
      res += fmt::format(
        "grad[{:3}] = {:>12} : {:<12} [A:FD] --- {:>12} : {:<12} [err:err(%)]\n",
        i,
        fmt::format("{:.5}", A),
        fmt::format("{:.5}", FD),
        fmt::format("{:.5}", err),
        fmt::format("{:.5}", 100 * err / scale)
      );
    }
    return res;
  }

  /*
  //       _                 _     _
  //      | | __ _  ___ ___ | |__ (_) __ _ _ __
  //   _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \
  //  | |_| | (_| | (_| (_) | |_) | | (_| | | | |
  //   \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|
  */
  template <typename FUNCTION, typename Number>
  static
  bool
  finite_difference_jacobian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    integer      dim_f,
    Number       Jac[],
    integer      ldJ,
    Number *     work,
    integer      lwork
  ) {

    static finite_difference_epsilon<Number> EPS;

    UTILS_ASSERT(
      lwork >= 3*dim_f,
      "finite_difference_jacobian(...,dim_f={},...,lwork={}), lwork must be >= 3*dim_f\n",
      dim_f, lwork
    );

    Number * vC{work};
    Number * vR{work+dim_f};
    Number * vL{work+2*dim_f};

    auto FUN = [&fun,dim_f]( Number const x[], Number f[] ) -> bool {
      bool ok{fun( x, f )};
      for ( integer i{0}; ok && i < dim_f; ++i ) ok = Utils::is_finite(f[i]);
      return ok;
    };

    bool ok_C{FUN( x, vC )};
    if ( !ok_C ) return false;

    Number * X{const_cast<Number*>(x)};
    Number * pjac{Jac};

    for ( integer j{0}; j < dim_x; ++j ) {
      Number temp{x[j]};
      Number h1{EPS.epsilon1(temp)};
      Number h2{EPS.epsilon2(temp)};

      X[j] = temp+h1; bool ok_R{FUN( X, vR )};
      X[j] = temp-h1; bool ok_L{FUN( X, vL )};

      integer ic = finite_difference_centered( vL, ok_L, vC, ok_C, vR, ok_R, h1, pjac, dim_f );

      switch ( ic ) {
      case  0:
        break;
      case  1:
        {
          Number * vRR{vL};
          X[j] = temp+h2; // modify the vector only at i position
          bool ok_RR{FUN( X, vRR )};
          if ( ok_RR ) ok_RR = finite_difference_side( vC, vR, vRR, h1, h2, pjac, dim_f );
          if ( !ok_RR ) {
            for ( integer i{0}; i < dim_f; ++i ) {
              pjac[i] = (vR[i]-vC[i])/h1; // low precision FD
              if ( !Utils::is_finite( pjac[i] ) ) return false;
            }
          }
        }
        break;
      case -1:
        {
          Number * vLL{vR};
          X[j] = temp-h2; // modify the vector only at i position
          bool ok_LL{FUN( X, vLL )};
          if ( ok_LL ) ok_LL = finite_difference_side( vC, vL, vLL, -h1, -h2, pjac, dim_f );
          if ( !ok_LL ) {
            for ( integer i{0}; i < dim_f; ++i ) {
              pjac[i] = (vC[i]-vL[i])/h1; // low precision FD
              if ( !Utils::is_finite( pjac[i] ) ) return false;
            }
          }
        }
        break;
      case -2:
        return false;
        break;
      }
      X[j] = temp; // restore i position
      pjac += ldJ;
    }
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  static
  void
  finite_difference_check_jacobian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    integer      dim_f,
    Number const Jac[],
    integer      ldJ,
    Number       epsi,
    Number *     work,
    integer      lwork,
    vector<tuple<integer,integer,Number,Number,Number,Number>> & err_list
  ) {

    static finite_difference_epsilon<Number> EPS;

    UTILS_ASSERT(
      lwork >= 4*dim_f,
      "finite_difference_check_jacobian(...,dim_f={},...,lwork={}), lwork must be >= 3*dim_f\n",
      dim_f, lwork
    );

    Number * vC{work};
    Number * vR{work+dim_f};
    Number * vL{work+2*dim_f};
    Number * Jcol{work+3*dim_f};

    auto FUN = [&fun,dim_f]( Number const x[], Number f[] ) -> bool {
      bool ok{fun( x, f )};
      for ( integer i{0}; ok && i < dim_f; ++i ) ok = Utils::is_finite(f[i]);
      return ok;
    };

    bool ok_C{FUN( x, vC )};
    if ( !ok_C ) { err_list.emplace_back(-1,-1,0,0,0,0); return; }

    Number       * X{const_cast<Number*>(x)};
    Number const * Ajac{Jac};

    for ( integer j{0}; j < dim_x; ++j ) {
      Number temp{x[j]};
      Number h1{EPS.epsilon1(temp)};
      Number h2{EPS.epsilon2(temp)};

      X[j] = temp+h1; bool ok_R{FUN( X, vR )};
      X[j] = temp-h1; bool ok_L{FUN( X, vL )};

      integer ic = finite_difference_centered( vL, ok_L, vC, ok_C, vR, ok_R, h1, Jcol, dim_f );

      switch ( ic ) {
      case  0:
        break;
      case  1:
        {
          Number * vRR{vL};
          X[j] = temp+h2; // modify the vector only at i position
          bool ok_RR{FUN( X, vRR )};
          if ( ok_RR ) ok_RR = finite_difference_side( vC, vR, vRR, h1, h2, Jcol, dim_f );
          if ( !ok_RR )
            for ( integer i = 0; i < dim_f; ++i )
              Jcol[i] = (vR[i]-vC[i])/h1; // low precision FD
        }
        break;
      case -1:
        {
          Number * vLL{vR};
          X[j] = temp-h2; // modify the vector only at i position
          bool ok_LL{FUN( X, vLL )};
          if ( ok_LL ) ok_LL = finite_difference_side( vC, vL, vLL, -h1, -h2, Jcol, dim_f );
          if ( !ok_LL )
            for ( integer i = 0; i < dim_f; ++i )
              Jcol[i] = (vC[i]-vL[i])/h1; // low precision FD
        }
        break;
      case -2:
        break;
      }

      X[j] = temp; // restore i position
      for ( integer i{0}; i < dim_f; ++i ) {
        Number A{Ajac[i]};
        Number FD{Jcol[i]};
        Number scale{max( Number(1), max(abs(A), abs(FD)))};
        Number err{abs(A - FD)};
        if ( !Utils::is_finite(A)  ||
             !Utils::is_finite(FD) ||
             err > epsi * scale )
          err_list.emplace_back( i, j, A, FD, err, scale );
      }
      Ajac += ldJ;
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  static
  string
  finite_difference_check_jacobian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    integer      dim_f,
    Number const Jac[],
    integer      ldJ,
    Number       epsi,
    Number *     work,
    integer      lwork
  ) {
    vector<tuple<integer,integer,Number,Number,Number,Number>> err_list;
    finite_difference_check_jacobian( x, dim_x, fun, dim_f, Jac, ldJ, epsi, work, lwork, err_list );
    string res{""};
    for ( auto & e : err_list ) {
      integer i, j;
      Number A, FD, err, scale;
      tie ( i, j, A, FD, err, scale ) = e;
      res += fmt::format(
        "jac[{:3},{:3}] = {:>12} : {:<12} [A:FD] --- {:>12} : {:<12} [err:err(%)]\n",
        i, j,
        fmt::format("{:.5}", A),
        fmt::format("{:.5}", FD),
        fmt::format("{:.5}", err),
        fmt::format("{:.5}", 100 * err / scale)
      );
    }
    return res;
  }

  /*
  //   _   _               _
  //  | | | | ___  ___ ___(_) __ _ _ __
  //  | |_| |/ _ \/ __/ __| |/ _` | '_ \
  //  |  _  |  __/\__ \__ \ | (_| | | | |
  //  |_| |_|\___||___/___/_|\__,_|_| |_|
  */

  template <typename FUNCTION, typename Number>
  static
  bool
  finite_difference_hessian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number       Hess[],
    integer      ldH
  ) {
    static finite_difference_epsilon<Number> EPS;

    auto FUN = [&fun]( Number const x[], Number & f ) -> bool {
      bool ok{fun( x, f )};
      if ( ok ) ok = Utils::is_finite(f);
      return ok;
    };

    Number fp{0}, fm{0}, fc{0}, hij{0}, tempi{0}, tempj{0}, hi{0}, hj{0},
           fpp{0}, fpm{0}, fmp{0}, fmm{0};
    bool ok{true};

    Number * X{const_cast<Number*>(x)};
    for ( integer j{0}; j < dim_x && ok; ++j ) {
      tempj = x[j];
      hj    = EPS.epsilon3(tempj);
      ok    = FUN( X, fc );
      if ( !ok ) goto skip;
      X[j] = tempj+hj;
      ok   = FUN( X, fp );
      if ( !ok ) goto skip;
      X[j] = tempj-hj;
      ok   = FUN( X, fm );
      if ( !ok ) goto skip;
      Hess[j*(ldH+1)] = ((fp+fm)-2*fc)/(hj*hj);
      for ( integer i{j+1}; i < dim_x && ok; ++i ) {
        tempi = X[i];
        hi    = EPS.epsilon3(tempi);
        X[i]  = tempi+hi;
        X[j]  = tempj+hj;
        ok    = FUN( X, fpp );
        if ( !ok ) goto skip2;
        X[i] = tempi-hi;
        ok   = FUN( X, fmp );
        if ( !ok ) goto skip2;
        X[j] = tempj-hj;
        ok   = FUN( X, fmm );
        if ( !ok ) goto skip2;
        X[i] = tempi+hi;
        ok   = FUN( X, fpm );
        if ( !ok ) goto skip2;
        hij = 4*hi*hj;
        Hess[j+i*ldH] = Hess[i+j*ldH] = ( (fpp+fmm) - (fpm+fmp) )/hij;
      skip2:
        X[i] = tempi;
      }
    skip:
      X[j] = tempj;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  static
  void
  finite_difference_check_hessian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number const Hess[],
    integer      ldH,
    Number       epsi,
    vector<tuple<integer,integer,Number,Number,Number,Number>> err_list
  ) {
    static finite_difference_epsilon<Number> EPS;

    auto FUN = [&fun]( Number const x[], Number & f ) -> bool {
      bool ok{fun( x, f )};
      if ( ok ) ok = Utils::is_finite(f);
      return ok;
    };

    Number fc{0}, fp{0}, fm{0}, fpp{0}, fpm{0}, fmp{0}, fmm{0},
           ddji{0}, ddij{0}, dde{0}, dd{0};

    Number * X{const_cast<Number*>(x)};
    for ( integer j{0}; j < dim_x; ++j ) {
      Number tempj{x[j]};
      Number hj{EPS.epsilon3(tempj)};

      bool ok{FUN( X, fc )};
      if ( ok ) { X[j] = tempj+hj; ok = FUN( X, fp ); }
      if ( ok ) { X[j] = tempj-hj; ok = FUN( X, fm ); }
      if ( ok ) {
        dde = Hess[j*(ldH+1)];
        dd  = ((fp+fm)-2*fc)/(hj*hj);
        ok  = Utils::is_finite(dd);
      }
      if ( ok ) {
        Number scale{1+max(abs(dd),abs(dde))};
        Number err{abs(dd-dde)};
        if ( err > epsi*scale ) err_list.emplace_back( j, j, dd, dde, err, scale );
        for ( integer i{j+1}; i < dim_x && ok; ++i ) {
          Number tempi{X[i]};
          Number hi{EPS.epsilon3(tempi)};
          X[i] = tempi+hi;
          X[j] = tempj+hj;
          ok   = FUN( X, fpp );
          if ( ok ) { X[i] = tempi-hi; ok = FUN( X, fmp ); }
          if ( ok ) { X[j] = tempj-hj; ok = FUN( X, fmm ); }
          if ( ok ) { X[i] = tempi+hi; ok = FUN( X, fpm ); }
          if ( ok ) {
            Number hij{4*hi*hj};
            ddji = Hess[j+i*ldH];
            ddij = Hess[i+j*ldH];
            dd   = ( (fpp+fmm) - (fpm+fmp) )/hij;
            ok   = Utils::is_finite(dd);
          }
          X[i] = tempi;
          if ( ok ) {
            scale = 1+max(abs(dd),abs(ddij));
            err   = abs(dd-ddij);
            if ( err > epsi*scale ) err_list.emplace_back( i, j, dd, ddij, err, scale );
            scale = 1+max(abs(dd),abs(ddji));
            err   = abs(dd-ddji);
            if ( err > epsi*scale ) err_list.emplace_back( j, i, dd, ddji, err, scale );
          } else {
            err_list.emplace_back( i, j, dd, ddij, 0, 1 );
          }
        }
      } else {
        err_list.emplace_back( -1, -1, 0, 0, 0, 0 );
      }
      X[j] = tempj;
    }
  }

  template <typename FUNCTION, typename Number>
  static
  string
  finite_difference_check_hessian(
    Number const x[],
    integer      dim_x,
    FUNCTION     fun,
    Number const Hess[],
    integer      ldH,
    Number       epsi
  ) {
    vector<tuple<integer,integer,Number,Number,Number,Number>> err_list;
    finite_difference_check_hessian( x, dim_x, fun, Hess, ldH, epsi, err_list );
    string res{""};
    for ( auto & e : err_list ) {
      integer i,j;
      Number A, FD, err, scale;
      tie ( i, j, A, FD, err, scale ) = e;
      res += fmt::format(
        "hessian[{:3},{:3}] = {:>12} : {:<12} [A:FD] --- {:>12} : {:<12} [err:err(%)]\n",
        i, j,
        fmt::format("{:.5}", A),
        fmt::format("{:.5}", FD),
        fmt::format("{:.5}", err),
        fmt::format("{:.5}", 100 * err / scale)
      );
    }
    return res;
  }

} // end namespace alglin

///
/// eof: Alglin_FD.hxx
///
