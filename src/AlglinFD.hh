/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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
/// file: alglinFD.hh
///

#include "Alglin.hh"
#include <iostream>
#include <iomanip>

#ifndef ALGLIN_FD_HH
#define ALGLIN_FD_HH

namespace alglin {

  /*
  //    ____               _ _            _
  //   / ___|_ __ __ _  __| (_) ___ _ __ | |_
  //  | |  _| '__/ _` |/ _` | |/ _ \ '_ \| __|
  //  | |_| | | | (_| | (_| | |  __/ | | | |_
  //   \____|_|  \__,_|\__,_|_|\___|_| |_|\__|
  */
  template <typename FUNCTION, typename Number>
  inline
  bool
  finite_difference_gradient( integer      fd_gradient,
                              Number const x[],
                              integer      dim_x,
                              FUNCTION   * fun,
                              Number       grad[] ) {

    Number const eps = fd_gradient == 0 ?
                       cbrt(std::numeric_limits<Number>::epsilon()):
                       sqrt(std::numeric_limits<Number>::epsilon());

    Number value0=0, value1=0; // only to stop warning
    bool ok = true;
    if ( fd_gradient != 0 ) {
      ok = (*fun)( x, value0 ) && isRegular(value0);
      if ( !ok ) return false;
    }

    Number * X = const_cast<Number*>(x);

    for ( integer i = 0; i < dim_x && ok; ++i ) {
      Number temp = x[i];
      Number h    = std::max( eps*std::abs(temp), eps );
      switch ( fd_gradient ) {
      case 0:
        X[i] = temp+h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && isRegular(value1);
        if ( !ok ) break;

        X[i] = temp-h; // modify the vector only at i position
        ok = (*fun)( X, value0 ) && isRegular(value0);
        if ( ok ) grad[i] = (value1-value0)/(2*h);
        break;

      case 1:
        X[i] = temp+h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && isRegular(value1);
        if ( ok ) grad[i] = (value1-value0)/h;
        break;

      case -1:
        X[i] = temp-h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && isRegular(value1);
        if ( ok ) grad[i] = (value0-value1)/h;
        break;
      default:
        ALGLIN_ERROR( "finite_difference_gradient::fd_gradient = " <<
                      fd_gradient << " must be in {-1, 0, 1}" );
      }
      X[i] = temp; // restore i position
      ok = isRegular(grad[i]);
    }
    return ok;
  }

  template <typename FUNCTION, typename Number>
  inline
  void
  finite_difference_check_gradient( Number const   x[],
                                    integer        dim_x,
                                    FUNCTION     * fun,
                                    Number const   grad[],
                                    Number         epsi,
                                    ostream_type & stream ) {

    Number const eps = cbrt(std::numeric_limits<Number>::epsilon());
    Number * X = const_cast<Number*>(x);

    bool ok = true;
    Number value0, value1, gradi=0;
    for ( integer i = 0; i < dim_x && ok; ++i ) {
      Number temp = x[i];
      Number h    = std::max( eps*std::abs(temp), eps );
      X[i] = temp+h; // modify the vector only at i position
      ok = (*fun)( X, value1 ) && isRegular(value1);
      if ( !ok ) break;

      X[i] = temp-h; // modify the vector only at i position
      ok = (*fun)( X, value0 ) && isRegular(value0);
      if ( ok ) {
        gradi = (value1-value0)/(2*h);
        ok    = isRegular(gradi);
      }

      X[i] = temp; // restore i position
      if ( ok ) {
        Number scale = std::max(eps,std::max(std::abs(gradi),std::abs(grad[i])));
        Number err   = std::abs(gradi-grad[i]);
        if ( err > epsi*std::max(Number(1),scale) ) {
          stream << "grad[" << std::setw(3) << i << "] = "
                 << std::setw(14) << grad[i] << " [A] --- "
                 << std::setw(14) << gradi   << " [FD]  err = "
                 << std::setw(14) << err << "  err (%) = "
                 << 100*err/scale << '\n';
        }
      }
    }
    if ( !ok )
      stream << "failed finite difference in `finite_difference_check_gradient`\n";
  }

  // ---------------------------------------------------------------------------

  template <typename Number>
  inline
  bool
  isRegular( Number const v[], integer n ) {
    bool ok = true;
    for ( integer i = 0; i < n && ok; ++i ) ok = alglin::isRegular(v[i]);
    return ok;
  }

  /*
  //       _                 _     _
  //      | | __ _  ___ ___ | |__ (_) __ _ _ __
  //   _  | |/ _` |/ __/ _ \| '_ \| |/ _` | '_ \
  //  | |_| | (_| | (_| (_) | |_) | | (_| | | | |
  //   \___/ \__,_|\___\___/|_.__/|_|\__,_|_| |_|
  */
  template <typename FUNCTION, typename Number>
  bool
  finite_difference_jacobian( integer      fd_jacobian,
                              Number const x[],
                              integer      dim_x,
                              FUNCTION   * fun,
                              integer      dim_f,
                              Number       Jac[],
                              integer      ldJ ) {

    Number const eps = fd_jacobian == 0 ?
                       cbrt(std::numeric_limits<Number>::epsilon()):
                       sqrt(std::numeric_limits<Number>::epsilon());

    std::vector<Number> tmp(size_t(2*dim_f));
    Number * g0 = &tmp[size_t(0)];
    Number * g1 = &tmp[size_t(dim_f)];
    bool ok = true;

    if ( fd_jacobian != 0 ) {
      ok = (*fun)( x, g0 ) && isRegular(g0,dim_f);
      if ( !ok ) return false;
    }

    Number * X = const_cast<Number*>(x);
    Number * pjac = Jac;

    for ( integer j = 0; j < dim_x && ok; ++j ) {
      Number temp = x[j];
      Number h    = std::max( eps*std::abs(temp), eps );
      switch ( fd_jacobian ) {
      case 0:
        X[j] = temp+h; // modify the vector only at j position
        ok = (*fun)( X, g1 ) && isRegular(g1,dim_f);
        if ( !ok ) break;

        X[j] = temp-h; // modify the vector only at j position
        ok = (*fun)( X, g0 ) && isRegular(g0,dim_f);
        if ( ok )
          for ( integer i = 0; i < dim_f; ++i )
            pjac[i] = (g1[i]-g0[i])/(2*h);
        break;

      case 1:
        X[j] = temp+h; // modify the vector only at j position
        ok = (*fun)( X, g1 ) && isRegular(g1,dim_f);
        if ( ok )
          for ( integer i = 0; i < dim_f; ++i )
            pjac[i] = (g1[i]-g0[i])/h;
        break;

      case -1:
        X[j] = temp-h; // modify the vector only at i position
        ok = (*fun)( X, g1 ) && isRegular(g1,dim_f);
        if ( ok )
          for ( integer i = 0; i < dim_f; ++i )
            pjac[j] = (g0[i]-g1[i])/h;
        break;
      default:
        ALGLIN_ERROR( "finite_difference_jacobian::fd_gradient = " <<
                      fd_jacobian << " must be in {-1, 0, 1}" );
      }
      X[j] = temp; // modify the vector only at i position
      ok = isRegular(pjac,dim_f);
      pjac += ldJ;
    }
    return ok;
  }

  template <typename FUNCTION, typename Number>
  void
  finite_difference_check_jacobian( Number const   x[],
                                    integer        dim_x,
                                    FUNCTION    *  fun,
                                    integer        dim_f,
                                    Number const   Jac[],
                                    integer        ldJ,
                                    Number         epsi,
                                    ostream_type & stream  ) {

    Number const eps = cbrt(std::numeric_limits<Number>::epsilon());

    std::vector<Number> tmp(size_t(2*dim_f));
    Number * g0 = &tmp[size_t(0)];
    Number * g1 = &tmp[size_t(dim_f)];
    bool ok = true;

    Number       * X    = const_cast<Number*>(x);
    Number const * pjac = Jac;

    for ( integer j = 0; j < dim_x && ok; ++j ) {
      Number temp = x[j];
      Number h    = std::max( eps*std::abs(temp), eps );

      X[j] = temp+h; // modify the vector only at j position
      ok = (*fun)( X, g1 ) && isRegular(g1,dim_f);
      if ( !ok ) break;

      X[j] = temp-h; // modify the vector only at j position
      ok = (*fun)( X, g0 ) && isRegular(g0,dim_f);

      for ( integer i = 0; i < dim_f && ok; ++i ) {
        ok = isRegular(pjac,dim_f);
        Number d     = (g1[i]-g0[i])/(2*h);
        Number scale = std::max(eps,std::max(std::abs(d),std::abs(pjac[i])));
        Number err   = std::abs(d-pjac[i]);
        if ( err > epsi*std::max(Number(1),scale) ) {
          stream << "jac[" << std::setw(3) << i << ", "
                 << std::setw(3) << j << "] = "
                 << std::setw(14) << pjac[i] << " [A] --- "
                 << std::setw(14) << d << " [FD]  err = "
                 << std::setw(14) << err << "  err (%) = "
                 << 100*err/scale << '\n';
        }
      }

      X[j] = temp; // modify the vector only at i position
      pjac += ldJ;
    }
    if ( !ok )
      stream << "failed finite difference in `finite_difference_check_jacobian`\n";
  }

  /*
  //   _   _               _
  //  | | | | ___  ___ ___(_) __ _ _ __
  //  | |_| |/ _ \/ __/ __| |/ _` | '_ \
  //  |  _  |  __/\__ \__ \ | (_| | | | |
  //  |_| |_|\___||___/___/_|\__,_|_| |_|
  */
  template <typename FUNCTION, typename Number>
  bool
  finite_difference_hessian( Number const x[],
                             integer      dim_x,
                             FUNCTION   * fun,
                             Number       Hess[],
                             integer      ldH ) {

    Number const eps = pow(std::numeric_limits<Number>::epsilon(),0.25);
    bool ok = true;

    Number * X = const_cast<Number*>(x);
    Number fpp, fpm, fmp, fmm;
    for ( integer j = 0; j < dim_x && ok; ++j ) {
      Number tempj = x[j];
      Number hj    = std::max( eps*std::abs(tempj), eps );

      Number fp, fm, fc;
      ok = (*fun)( X, fc ) && alglin::isRegular(fc);
      if ( !ok ) break;
      X[j] = tempj+hj;
      ok = (*fun)( X, fp ) && alglin::isRegular(fp);
      if ( !ok ) break;
      X[j] = tempj-hj;
      ok = (*fun)( X, fm ) && alglin::isRegular(fm);
      if ( !ok ) break;
      Hess[j*(ldH+1)] = ((fp+fm)-2*fc)/(hj*hj);
      for ( integer i = j+1; i < dim_x && ok; ++i ) {
        Number tempi = X[i];
        Number hi    = std::max( eps*std::abs(tempi), eps );
        X[i] = tempi+hi;
        X[j] = tempj+hj;
        ok = (*fun)( X, fpp ) && alglin::isRegular(fpp);
        if ( !ok ) break;
        X[i] = tempi-hi;
        ok = (*fun)( X, fmp ) && alglin::isRegular(fmp);
        if ( !ok ) break;
        X[j] = tempj-hj;
        ok = (*fun)( X, fmm ) && alglin::isRegular(fmm);
        if ( !ok ) break;
        X[i] = tempi+hi;
        ok = (*fun)( X, fpm ) && alglin::isRegular(fpm);
        if ( !ok ) break;
        Number hij = 4*hi*hj;
        Hess[j+i*ldH] = Hess[i+j*ldH] = ( (fpp+fmm) - (fpm+fmp) )/hij;
        X[i] = tempi;
      }
      X[j] = tempj;
    }
    return ok;
  }


  template <typename FUNCTION, typename Number>
  bool
  finite_difference_check_hessian( Number const   x[],
                                   integer        dim_x,
                                   FUNCTION     * fun,
                                   Number const   Hess[],
                                   integer        ldH,
                                   Number         epsi,
                                   ostream_type & stream ) {

    Number const eps = pow(std::numeric_limits<Number>::epsilon(),0.25);
    bool ok = true;

    Number * X = const_cast<Number*>(x);
    Number fpp, fpm, fmp, fmm;
    for ( integer j = 0; j < dim_x && ok; ++j ) {
      Number tempj = x[j];
      Number hj    = std::max( eps*std::abs(tempj), eps );

      Number fp, fm, fc;
      ok = (*fun)( X, fc ) && alglin::isRegular(fc);
      if ( !ok ) break;
      X[j] = tempj+hj;
      ok = (*fun)( X, fp ) && alglin::isRegular(fp);
      if ( !ok ) break;
      X[j] = tempj-hj;
      ok = (*fun)( X, fm ) && alglin::isRegular(fm);
      if ( !ok ) break;

      Number dde = Hess[j*(ldH+1)];
      Number dd  = ((fp+fm)-2*fc)/(hj*hj);
      ok = alglin::isRegular(dd);
      if ( !ok ) break;
      Number scale = std::max(eps,std::max(std::abs(dd),std::abs(dde)));
      Number err   = std::abs(dd-dde);
      if ( err > epsi*std::max(Number(1),scale) ) {
        stream << "Hess[" << std::setw(3) << j << ", "
               << std::setw(3) << j << "] = "
               << std::setw(14) << dde << " [A] --- "
               << std::setw(14) << dd << " [FD]  err = "
               << std::setw(14) << err << "  err (%) = "
               << 100*err/scale << '\n';
      }

      for ( integer i = j+1; i < dim_x && ok; ++i ) {
        Number tempi = X[i];
        Number hi    = std::max( eps*std::abs(tempi), eps );
        X[i] = tempi+hi;
        X[j] = tempj+hj;
        ok = (*fun)( X, fpp ) && alglin::isRegular(fpp);
        if ( !ok ) break;
        X[i] = tempi-hi;
        ok = (*fun)( X, fmp ) && alglin::isRegular(fmp);
        if ( !ok ) break;
        X[j] = tempj-hj;
        ok = (*fun)( X, fmm ) && alglin::isRegular(fmm);
        if ( !ok ) break;
        X[i] = tempi+hi;
        ok = (*fun)( X, fpm ) && alglin::isRegular(fpm);
        if ( !ok ) break;
        Number hij  = 4*hi*hj;
        Number ddji = Hess[j+i*ldH];
        Number ddij = Hess[i+j*ldH];
        dd   = ( (fpp+fmm) - (fpm+fmp) )/hij;
        ok = alglin::isRegular(dd);
        if ( !ok ) break;
        scale = std::max(eps,std::max(std::abs(dd),std::abs(ddij)));
        err   = std::abs(dd-ddij);
        if ( err > epsi*std::max(Number(1),scale) ) {
          stream << "Hess[" << std::setw(3) << i << ", " << std::setw(3) << j
                 << "] = " << std::setw(14) << ddij << " [A] --- "
                 << std::setw(14) << dd << " [FD]  err = "
                 << std::setw(14) << err << "  err (%) = "
                 << 100*err/scale << '\n';
        }
        scale = std::max(eps,std::max(std::abs(dd),std::abs(ddji)));
        err   = std::abs(dd-ddji);
        if ( err > epsi*std::max(Number(1),scale) ) {
          stream << "Hess[" << std::setw(3) << j << ", " << std::setw(3) << i
                 << "] = " << std::setw(14) << ddij << " [A] --- "
                 << std::setw(14) << dd << " [FD]  err = "
                 << std::setw(14) << err << "  err (%) = "
                 << 100*err/scale << '\n';
        }
        X[i] = tempi;
      }
      X[j] = tempj;
    }
    return ok;
  }

} // end namespace alglin

#endif

///
/// eof: AlglinFD.hh
///
