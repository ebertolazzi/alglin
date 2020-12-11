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
  finite_difference_gradient(
    integer      fd_gradient,
    Number const x[],
    integer      dim_x,
    FUNCTION   * fun,
    Number       grad[]
  ) {

    Number const eps = fd_gradient == 0 ?
                       cbrt(std::numeric_limits<Number>::epsilon()):
                       sqrt(std::numeric_limits<Number>::epsilon());

    Number value0=0, value1=0; // only to stop warning
    bool ok = true;
    if ( fd_gradient != 0 ) {
      ok = (*fun)( x, value0 ) && Utils::isRegular(value0);
      if ( !ok ) return false;
    }

    Number * X = const_cast<Number*>(x);

    for ( integer i = 0; i < dim_x && ok; ++i ) {
      Number temp = x[i];
      Number h    = std::max( eps*std::abs(temp), eps );
      switch ( fd_gradient ) {
      case 0:
        X[i] = temp+h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && Utils::isRegular(value1);
        if ( !ok ) break;

        X[i] = temp-h; // modify the vector only at i position
        ok = (*fun)( X, value0 ) && Utils::isRegular(value0);
        if ( ok ) grad[i] = (value1-value0)/(2*h);
        break;

      case 1:
        X[i] = temp+h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && Utils::isRegular(value1);
        if ( ok ) grad[i] = (value1-value0)/h;
        break;

      case -1:
        X[i] = temp-h; // modify the vector only at i position
        ok = (*fun)( X, value1 ) && Utils::isRegular(value1);
        if ( ok ) grad[i] = (value0-value1)/h;
        break;

      default:
        UTILS_ERROR(
          "finite_difference_gradient::fd_gradient = {} must be in {{-1, 0, 1}}\n",
          fd_gradient
        );
      }
      X[i] = temp; // restore i position
      ok   = Utils::isRegular(grad[i]);
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  inline
  void
  finite_difference_check_gradient(
    Number const   x[],
    integer        dim_x,
    FUNCTION     * fun,
    Number const   grad[],
    Number         epsi,
    ostream_type & stream
  ) {

    Number const eps = cbrt(std::numeric_limits<Number>::epsilon());
    Number * X = const_cast<Number*>(x);

    bool ok = true;
    Number value0, value1, gradi=0;
    for ( integer i = 0; i < dim_x && ok; ++i ) {
      Number temp = x[i];
      Number h    = std::max( eps*std::abs(temp), eps );
      X[i] = temp+h; // modify the vector only at i position
      ok = (*fun)( X, value1 ) && Utils::isRegular(value1);
      if ( !ok ) break;

      X[i] = temp-h; // modify the vector only at i position
      ok = (*fun)( X, value0 ) && Utils::isRegular(value0);
      if ( ok ) {
        gradi = (value1-value0)/(2*h);
        ok    = Utils::isRegular(gradi);
      }

      X[i] = temp; // restore i position
      if ( ok ) {
        Number scale = std::max(eps,std::max(std::abs(gradi),std::abs(grad[i])));
        Number err   = std::abs(gradi-grad[i]);
        if ( err > epsi*std::max(Number(1),scale) ) {
          fmt::print( stream,
            "grad[{:3}] = {:14.5} [A] --- {:14.5} [FD]  err = {:14.5} err (%) {:8.4}\n",
            i, grad[i], gradi, err, 100*err/scale
          );
        }
      }
    }
    if ( !ok )
      stream << "failed finite difference in `finite_difference_check_gradient`\n";
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename Number>
  inline
  bool
  isRegular( Number const v[], integer n ) {
    bool ok = true;
    for ( integer i = 0; i < n && ok; ++i ) ok = Utils::isRegular(v[i]);
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
  finite_difference_jacobian(
    integer      fd_jacobian,
    Number const x[],
    integer      dim_x,
    FUNCTION   * fun,
    integer      dim_f,
    Number       Jac[],
    integer      ldJ
  ) {

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
        UTILS_ERROR(
          "finite_difference_jacobian::fd_gradient = {} must be in {{-1, 0, 1}}\n",
          fd_jacobian
        );
      }
      X[j] = temp; // modify the vector only at i position
      if ( ok ) ok = isRegular(pjac,dim_f);
      pjac += ldJ;
    }
    return ok;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  template <typename FUNCTION, typename Number>
  void
  finite_difference_check_jacobian(
    Number const   x[],
    integer        dim_x,
    FUNCTION    *  fun,
    integer        dim_f,
    Number const   Jac[],
    integer        ldJ,
    Number         epsi,
    ostream_type & stream
  ) {

    Number const eps = cbrt(std::numeric_limits<Number>::epsilon());

    std::vector<Number> tmp(size_t(2*dim_f));
    Number * g0 = &tmp[size_t(0)];
    Number * g1 = &tmp[size_t(dim_f)];
    bool ok = true;

    Number       * X    = const_cast<Number*>(x);
    Number const * pjac = Jac;

    for ( integer j = 0; j < dim_x && ok; ++j ) {
      ok = isRegular(pjac,dim_f);
      if ( !ok ) break;

      Number temp = x[j];
      Number h    = std::max( eps*std::abs(temp), eps );

      X[j] = temp+h; // modify the vector only at j position
      ok = (*fun)( X, g1 ) && isRegular(g1,dim_f);
      if ( !ok ) { X[j] = temp; break; }

      X[j] = temp-h; // modify the vector only at j position
      ok = (*fun)( X, g0 ) && isRegular(g0,dim_f);
      X[j] = temp; // modify the vector only at i position
      if ( !ok ) break;

      for ( integer i = 0; i < dim_f; ++i ) {
        Number d     = (g1[i]-g0[i])/(2*h);
        Number scale = std::max(eps,std::max(std::abs(d),std::abs(pjac[i])));
        Number err   = std::abs(d-pjac[i]);
        if ( err > epsi*std::max(Number(1),scale) ) {
          fmt::print( stream,
            "jac[{:3},{:3}] = {:14.5} [A] --- {:14.5} [FD]  err = {:14.5} err (%) {:8.4}\n",
            i, j, pjac[i], d, err, 100*err/scale
          );
        }
      }

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
  finite_difference_hessian(
    Number const x[],
    integer      dim_x,
    FUNCTION   * fun,
    Number       Hess[],
    integer      ldH
  ) {
    Number fp, fm, fc, hij, tempi, tempj, hi, hj;
    Number const eps = pow(std::numeric_limits<Number>::epsilon(),0.25);
    bool ok = true;

    Number * X = const_cast<Number*>(x);
    Number fpp, fpm, fmp, fmm;
    for ( integer j = 0; j < dim_x && ok; ++j ) {
      tempj = x[j];
      hj    = std::max( eps*std::abs(tempj), eps );
      ok = (*fun)( X, fc ) && Utils::isRegular(fc);
      if ( !ok ) goto skip;
      X[j] = tempj+hj;
      ok = (*fun)( X, fp ) && Utils::isRegular(fp);
      if ( !ok ) goto skip;
      X[j] = tempj-hj;
      ok = (*fun)( X, fm ) && Utils::isRegular(fm);
      if ( !ok ) goto skip;
      Hess[j*(ldH+1)] = ((fp+fm)-2*fc)/(hj*hj);
      for ( integer i = j+1; i < dim_x && ok; ++i ) {
        tempi = X[i];
        hi    = std::max( eps*std::abs(tempi), eps );
        X[i] = tempi+hi;
        X[j] = tempj+hj;
        ok = (*fun)( X, fpp ) && Utils::isRegular(fpp);
        if ( !ok ) goto skip2;
        X[i] = tempi-hi;
        ok = (*fun)( X, fmp ) && Utils::isRegular(fmp);
        if ( !ok ) goto skip2;
        X[j] = tempj-hj;
        ok = (*fun)( X, fmm ) && Utils::isRegular(fmm);
        if ( !ok ) goto skip2;
        X[i] = tempi+hi;
        ok = (*fun)( X, fpm ) && Utils::isRegular(fpm);
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
  bool
  finite_difference_check_hessian(
    Number const   x[],
    integer        dim_x,
    FUNCTION     * fun,
    Number const   Hess[],
    integer        ldH,
    Number         epsi,
    ostream_type & stream
  ) {
    Number fc, fp, fm, fpp, fpm, fmp, fmm, tempi, tempj, hi, hj, hij;
    Number ddji, ddij, dde, dd, scale, err;
    Number const eps = pow(std::numeric_limits<Number>::epsilon(),0.25);
    bool ok = true;

    Number * X = const_cast<Number*>(x);
    for ( integer j = 0; j < dim_x && ok; ++j ) {
      tempj = x[j];
      hj    = std::max( eps*std::abs(tempj), eps );

      ok = (*fun)( X, fc ) && Utils::isRegular(fc);
      if ( !ok ) goto skip;
      X[j] = tempj+hj;
      ok = (*fun)( X, fp ) && Utils::isRegular(fp);
      if ( !ok ) goto skip;
      X[j] = tempj-hj;
      ok = (*fun)( X, fm ) && Utils::isRegular(fm);
      if ( !ok ) goto skip;

      dde = Hess[j*(ldH+1)];
      dd  = ((fp+fm)-2*fc)/(hj*hj);
      ok = Utils::isRegular(dd);
      if ( !ok ) goto skip;
      scale = std::max(eps,std::max(std::abs(dd),std::abs(dde)));
      err   = std::abs(dd-dde);
      if ( err > epsi*std::max(Number(1),scale) ) {
        fmt::print( stream,
          "Hess[{:3},{:3}] = {:14.5} [A] --- {:14.5} [FD]  err = {:14.5}  err (%) = {:8.4}\n",
          j, j, dde, dd, err, 100*err/scale
        );
      }

      for ( integer i = j+1; i < dim_x && ok; ++i ) {
        tempi = X[i];
        hi    = std::max( eps*std::abs(tempi), eps );
        X[i] = tempi+hi;
        X[j] = tempj+hj;
        ok = (*fun)( X, fpp ) && Utils::isRegular(fpp);
        if ( !ok ) goto skip2;
        X[i] = tempi-hi;
        ok = (*fun)( X, fmp ) && Utils::isRegular(fmp);
        if ( !ok ) goto skip2;
        X[j] = tempj-hj;
        ok = (*fun)( X, fmm ) && Utils::isRegular(fmm);
        if ( !ok ) goto skip2;
        X[i] = tempi+hi;
        ok = (*fun)( X, fpm ) && Utils::isRegular(fpm);
        if ( !ok ) goto skip2;
        hij  = 4*hi*hj;
        ddji = Hess[j+i*ldH];
        ddij = Hess[i+j*ldH];
        dd   = ( (fpp+fmm) - (fpm+fmp) )/hij;
        ok = Utils::isRegular(dd);
        if ( !ok ) goto skip2;
        scale = std::max(eps,std::max(std::abs(dd),std::abs(ddij)));
        err   = std::abs(dd-ddij);
        if ( err > epsi*std::max(Number(1),scale) ) {
          fmt::print( stream,
            "Hess[{:3},{:3}] = {:14.5} [A] --- {:14.5} [FD]  err = {:14.5}  err (%) = {:8.4}\n",
            i, j, ddij, dd, err, 100*err/scale
          );
        }
        scale = std::max(eps,std::max(std::abs(dd),std::abs(ddji)));
        err   = std::abs(dd-ddji);
        if ( err > epsi*std::max(Number(1),scale) ) {
          fmt::print( stream,
            "Hess[{:3},{:3}] = {:14.5} [A] --- {:14.5} [FD]  err = {:14.5}  err (%) = {:8.4}\n",
            j, i, ddij, dd, err, 100*err/scale
          );
        }
      skip2:
        X[i] = tempi;
      }
    skip:
      X[j] = tempj;
    }
    return ok;
  }

} // end namespace alglin

///
/// eof: Alglin_FD.hxx
///
