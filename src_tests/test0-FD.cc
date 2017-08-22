/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2008-2015                                                 |
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

#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include "AlglinFD.hh"

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wc99-extensions"
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

using namespace alglin ;

static
bool
fun( doublereal const x[4], doublereal & res ) {
  res = x[0] + sin(x[1]) + exp(x[2]*x[3]) + x[0]*x[1]*x[2] ;
  return true ;
}

static
bool
fun_grad( doublereal const x[4], doublereal grad[4] ) {
  grad[0] = 1+x[1]*x[2] ;
  grad[1] = cos(x[1])+ x[0]*x[2]  ;
  grad[2] = exp(x[2]*x[3])*x[3]+x[0]*x[1] ;
  grad[3] = exp(x[2]*x[3])*x[2] ;
  return true ;
}

static
bool
fun_jac( doublereal const x[4], doublereal jac[4*4] ) {

  jac[0+0*4] = 0 ;
  jac[0+1*4] = x[2] ;
  jac[0+2*4] = x[1] ;
  jac[0+3*4] = 0 ;
  
  jac[1+0*4] = x[2] ;
  jac[1+1*4] = -sin(x[1]) ;
  jac[1+2*4] = x[0] ;
  jac[1+3*4] = 0 ;

  jac[2+0*4] = x[1] ;
  jac[2+1*4] = x[0] ;
  jac[2+2*4] = exp(x[2]*x[3])*x[3]*x[3] ;
  jac[2+3*4] = exp(x[2]*x[3]) + exp(x[2]*x[3])*x[2]*x[3]  ;

  jac[3+0*4] = 0 ;
  jac[3+1*4] = 0 ;
  jac[3+2*4] = exp(x[2]*x[3])+exp(x[2]*x[3])*x[2]*x[3] ;
  jac[3+3*4] = exp(x[2]*x[3])*x[2]*x[2] ;

  return true ;
}

class fun_class {
public:
  bool operator () ( doublereal const x[4], doublereal & res ) {
    return fun(x,res) ;
  }
} ;

class fun_class1 {
public:
  bool operator () ( doublereal const x[4], doublereal res[4] ) {
    return fun_grad( x, res ) ;
  }
} ;

int
main() {

  using namespace std ;

  integer    dim_x = 4 ;
  doublereal x[4] = { 1, 2, 1, 4 } ;
  doublereal gradFD[4], grad[4], jac[4*4], jacFD[4*4] ;

  fun_class  ff ;
  fun_class1 gg ;

  bool ok = finite_difference_gradient( 0, x, dim_x, ff, gradFD ) ;
  ok = fun_grad( x, grad ) ;
  doublereal epsi = 1e-6 ;
  cout << "\n\nCheck Gradient\n" ;
  finite_difference_check_gradient( x, dim_x, ff, grad, epsi, cout ) ;
  cout << "Done\n" ;

  cout << "diff grad(FD)[0] = " << gradFD[0] - grad[0] << '\n' ;
  cout << "diff grad(FD)[1] = " << gradFD[1] - grad[1] << '\n' ;
  cout << "diff grad(FD)[2] = " << gradFD[2] - grad[2] << '\n' ;
  cout << "diff grad(FD)[3] = " << gradFD[3] - grad[3] << '\n' ;
  cout << "diff grad(FD)[0] = " << (gradFD[0] - grad[0])/max(1.0,abs(grad[0])) << '\n' ;
  cout << "diff grad(FD)[1] = " << (gradFD[1] - grad[1])/max(1.0,abs(grad[1])) << '\n' ;
  cout << "diff grad(FD)[2] = " << (gradFD[2] - grad[2])/max(1.0,abs(grad[2])) << '\n' ;
  cout << "diff grad(FD)[3] = " << (gradFD[3] - grad[3])/max(1.0,abs(grad[3])) << '\n' ;

  ok = fun_jac( x, jac ) ;
  ok = finite_difference_jacobian( 0, x, dim_x, gg, dim_x, jacFD, dim_x ) ;
  cout << "\n\nCheck Jacobian\n" ;
  finite_difference_check_jacobian( x, dim_x, gg, dim_x, jac, dim_x, epsi, cout ) ;
  cout << "Done\n" ;

  for ( int i = 0 ; i < dim_x ; ++i ) {
    for ( int j = 0 ; j < dim_x ; ++j ) {
      cout << "jac[" << i << "," << j << "] = " << setw(14) << jac[i+j*dim_x]
           << " err = " << abs(jac[i+j*dim_x]-jacFD[i+j*dim_x])
           << "\n" ;
    }
  }

  ok = finite_difference_hessian( x, dim_x, ff, jacFD, dim_x ) ;
  cout << "ok = " << (ok?"TRUE\n":"FALSE\n") ;
  for ( int i = 0 ; i < dim_x ; ++i ) {
    for ( int j = 0 ; j < dim_x ; ++j ) {
      cout << "Hess[" << i << "," << j << "] = " << setw(14) << jac[i+j*dim_x]
           << " HessFD[" << i << "," << j << "] = " << setw(14) << jacFD[i+j*dim_x]
           << " err = " << abs(jac[i+j*dim_x]-jacFD[i+j*dim_x])
           << "\n" ;
    }
  }

  cout << "\n\nCheck Hessian\n" ;
  finite_difference_check_hessian( x, dim_x, ff, jac, dim_x, epsi, cout ) ;
  cout << "Done\n" ;

  cout << "All done!\n" ;

  return 0 ;
}
