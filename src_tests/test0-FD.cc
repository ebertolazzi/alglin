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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

using namespace alglin;

static Utils::Console msg( &std::cout,4);

static
bool
fun( real_type const x[4], real_type & res ) {
  res = x[0] + sin(x[1]) + exp(x[2]*x[3]) + x[0]*x[1]*x[2];
  return true;
}

static
bool
fun_grad( real_type const x[4], real_type grad[4] ) {
  grad[0] = 1+x[1]*x[2];
  grad[1] = cos(x[1])+ x[0]*x[2];
  grad[2] = exp(x[2]*x[3])*x[3]+x[0]*x[1];
  grad[3] = exp(x[2]*x[3])*x[2];
  return true;
}

static
bool
fun_jac( real_type const x[4], real_type jac[4*4] ) {

  jac[0+0*4] = 0;
  jac[0+1*4] = x[2];
  jac[0+2*4] = x[1];
  jac[0+3*4] = 0;

  jac[1+0*4] = x[2];
  jac[1+1*4] = -sin(x[1]);
  jac[1+2*4] = x[0];
  jac[1+3*4] = 0;

  jac[2+0*4] = x[1];
  jac[2+1*4] = x[0];
  jac[2+2*4] = exp(x[2]*x[3])*x[3]*x[3];
  jac[2+3*4] = exp(x[2]*x[3]) + exp(x[2]*x[3])*x[2]*x[3];

  jac[3+0*4] = 0;
  jac[3+1*4] = 0;
  jac[3+2*4] = exp(x[2]*x[3])+exp(x[2]*x[3])*x[2]*x[3];
  jac[3+3*4] = exp(x[2]*x[3])*x[2]*x[2];

  return true;
}

class fun_class {
public:
  bool operator () ( real_type const x[4], real_type & res ) {
    return fun(x,res);
  }
};

class fun_class1 {
public:
  bool operator () ( real_type const x[4], real_type res[4] ) {
    return fun_grad( x, res );
  }
};

int
main() {

  using namespace std;

  integer   const dim_x{4};
  real_type const x[4]{ 1, 2, 1, 4 };
  real_type       gradFD[4], grad[4], jac[4*4], jacFD[4*4];

  fun_class  ff;
  fun_class1 gg;

  bool ok{ finite_difference_gradient( x, dim_x, ff, gradFD ) };
  if (ok) ok = fun_grad( x, grad );
  if (!ok) { fmt::print("test failed!\n"); return -1; }
  real_type epsi{1e-6};

  fmt::print(
    "\n\nCheck Gradient\n{}Done\n",
     finite_difference_check_gradient( x, dim_x, ff, grad, epsi )
  );

  fmt::print("diff grad(FD)[0] = {:>12.6}\n", gradFD[0] - grad[0]);
  fmt::print("diff grad(FD)[1] = {:>12.6}\n", gradFD[1] - grad[1]);
  fmt::print("diff grad(FD)[2] = {:>12.6}\n", gradFD[2] - grad[2]);
  fmt::print("diff grad(FD)[3] = {:>12.6}\n", gradFD[3] - grad[3]);
  fmt::print("diff grad(FD)[0] = {:>12.6}\n", (gradFD[0] - grad[0])/max(1.0,abs(grad[0])));
  fmt::print("diff grad(FD)[1] = {:>12.6}\n", (gradFD[1] - grad[1])/max(1.0,abs(grad[1])));
  fmt::print("diff grad(FD)[2] = {:>12.6}\n", (gradFD[2] - grad[2])/max(1.0,abs(grad[2])));
  fmt::print("diff grad(FD)[3] = {:>12.6}\n", (gradFD[3] - grad[3])/max(1.0,abs(grad[3])));

  std::vector<real_type> work(4*dim_x);

  ok = fun_jac( x, jac );
  if ( ok ) ok = finite_difference_jacobian(
    x, dim_x, gg, dim_x, jacFD, dim_x,
    work.data(), static_cast<integer>(work.size())
  );
  if (!ok) { fmt::print("test failed!\n"); return -1; }
  fmt::print(
    "\n\nCheck Jacobian\n{}Done\n",
    finite_difference_check_jacobian(
      x, dim_x, gg, dim_x, jac, dim_x,
      epsi, work.data(), static_cast<integer>(work.size())
    )
  );

  for ( int i{0}; i < dim_x; ++i ) {
    for ( int j{0}; j < dim_x; ++j ) {
      fmt::print(
        "jac[{},{}] = {:<12.6}   err = {:<12.6}\n",
        i, j, jac[i+j*dim_x], abs(jac[i+j*dim_x]-jacFD[i+j*dim_x])
      );
    }
  }

  ok = finite_difference_hessian( x, dim_x, ff, jacFD, dim_x );
  fmt::print("ok = {}\n",ok);
  for ( int i{0}; i < dim_x; ++i ) {
    for ( int j{0}; j < dim_x; ++j ) {
      fmt::print(
        "Hess[{0},{1}] = {2:<12.6} HessFD[{0},{1}] = {3:<12.6} err = {4:<12.6}\n",
        i, j, jac[i+j*dim_x], jacFD[i+j*dim_x], abs(jac[i+j*dim_x]-jacFD[i+j*dim_x])
      );
    }
  }

  fmt::print(
    "\n\nCheck Hessian\n{}Done\n",
    finite_difference_check_hessian( x, dim_x, ff, jac, dim_x, epsi )
  );

  msg.green( "All done!\n" );

  return 0;
}
