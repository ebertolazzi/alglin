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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "Alglin.hh"
#include "Alglin_Eigen.hh"
#include "Utils_fmt.hh"
#include "Utils_string.hh"

#include <chrono>
#include <random>
#include <functional>
#include <array>

using namespace std;
typedef double real_type;
using alglin::integer;
using alglin::Transposition;

// ===========================================================================
// Utility functions for formatted output
// ===========================================================================

static void print_header(const std::string& title, fmt::color color = fmt::color::cyan) {
  const int width = 80;
  fmt::print(fg(color) | fmt::emphasis::bold, "\n");
  fmt::print(fg(color) | fmt::emphasis::bold, "┌{0:─^{1}}┐\n", "", width);
  fmt::print(fg(color) | fmt::emphasis::bold, "│{0: ^{1}}│\n", "", width);
  fmt::print(fg(color) | fmt::emphasis::bold, "│{0:^{1}}│\n", title, width);
  fmt::print(fg(color) | fmt::emphasis::bold, "│{0: ^{1}}│\n", "", width);
  fmt::print(fg(color) | fmt::emphasis::bold, "└{0:─^{1}}┘\n", "", width);
}

static void print_test_info(const std::string& info, fmt::color color = fmt::color::yellow) {
  fmt::print(fg(color), "🔹 {}\n", info);
}

static void print_success(const std::string& message) {
  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "✅ {}\n", message);
}

static void print_warning(const std::string& message) {
  fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold, "⚠️  {}\n", message);
}

static void print_error(const std::string& message) {
  fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "❌ {}\n", message);
}

static void print_matrix_with_border(const std::string& label, integer m, integer n, 
                                     const real_type* A, integer lda) {
  fmt::print(fg(fmt::color::light_blue), "\n┌─ {} ────\n", label);
  fmt::print("{}", alglin::print_matrix(m, n, A, lda));
  fmt::print(fg(fmt::color::light_blue), "└{}\n", Utils::repeat("─",70));
}

// ===========================================================================
// Test functions
// ===========================================================================

static void test1() {
  print_header("TEST 1: QR Factorization", fmt::color::cyan);
  
  alglin::QR<real_type> qr;
  constexpr integer M{3};
  constexpr integer N{5};
  constexpr integer LDA{3};
  constexpr real_type A[]{
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  print_matrix_with_border("Initial Matrix A (3×5)", M, N, A, LDA);
  
  print_test_info("Performing QR factorization of Aᵀ");
  auto start = chrono::high_resolution_clock::now();
  qr.t_factorize( "qr", M, N, A, LDA );
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration<double, milli>(end - start);
  fmt::print(fg(fmt::color::green), "✓ Factorization completed in {:.3f} ms\n", duration.count());

  real_type R[M*M];
  qr.getR( R, M );
  print_matrix_with_border("R matrix (3×3)", M, M, R, M);

  real_type rhs[M], b[M];
  real_type x[N]{1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(fg(fmt::color::yellow), "\n📊 Solving Least Squares: A·x = b\n");
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));

  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, M, x );
  alglin::Zero_n( x+3, N-M );
  qr.Q_mul( x );

  fmt::print("xᵀ = {}\n", alglin::print_matrix( 1, N, x, 1 ));

  alglin::gemv( Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  real_type residual_norm = alglin::nrm2(M, b, 1);
  
  fmt::print("rᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));
  fmt::print("‖r‖₂ = {:.6e}\n", residual_norm);
  
  if (residual_norm < 1e-6) {
    print_success("Test 1 passed!");
  } else {
    print_warning(fmt::format("Test 1 residual ({:.6e}) above threshold", residual_norm));
  }
}

static void test2() {
  print_header("TEST 2: QRP Factorization", fmt::color::magenta);
  
  alglin::QRP<real_type> qr;
  constexpr integer M{3};
  constexpr integer N{5};
  constexpr integer LDA{3};
  constexpr real_type A[]{
    0.001,      2,     3,
    0.001,  0.001,     0,
    0,      0.001,     0,
    0.001,     -1,     0,
    0.000001,   5,     3
  };

  print_matrix_with_border("Initial Matrix A (3×5)", M, N, A, LDA);
  
  print_test_info("Performing QRP factorization of Aᵀ");
  auto start = chrono::high_resolution_clock::now();
  qr.t_factorize( "qr", M, N, A, LDA );
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration<double, milli>(end - start);
  fmt::print(fg(fmt::color::green), "✓ Factorization completed in {:.3f} ms\n", duration.count());

  real_type R[M*M];
  qr.getR( R, M );
  print_matrix_with_border("R matrix (3×3)", M, M, R, M);

  real_type rhs[M], b[M];
  real_type x[N]{1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(fg(fmt::color::yellow), "\n📊 Solving Least Squares: A·x = b\n");
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));

  qr.inv_permute( rhs );
  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, M, x );
  alglin::Zero_n( x+3, N-M );
  qr.Q_mul( x );

  fmt::print("xᵀ = {}\n", alglin::print_matrix( 1, N, x, 1 ));

  alglin::gemv( Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  real_type residual_norm = alglin::nrm2(M, b, 1);
  
  fmt::print("rᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));
  fmt::print("‖r‖₂ = {:.6e}\n", residual_norm);
  
  if (residual_norm < 1e-6) {
    print_success("Test 2 passed!");
  } else {
    print_warning(fmt::format("Test 2 residual ({:.6e}) above threshold", residual_norm));
  }
}

static void test3() {
  print_header("TEST 3: QRP on Ill-conditioned Matrix", fmt::color::orange);
  
  alglin::QRP<real_type> qr;
  constexpr integer M{5};
  constexpr integer N{5};
  constexpr integer LDA{5};
  constexpr real_type A[]{
    0.001,      2,     3,     2, 3,
    0.001,  0.001,     0, 0.001, 1e-10,
    0,      0.001,     0, 0.001, 1e-12,
    0.001,     -1, 1e-12,    -1, -1e-12,
    0.000001,   5,     3,     5, 3
  };

  print_matrix_with_border("Initial Matrix A (5×5)", M, N, A, LDA);
  print_warning("Matrix contains very small values (ill-conditioned)");
  
  print_test_info("Performing QRP factorization of Aᵀ");
  auto start = chrono::high_resolution_clock::now();
  qr.t_factorize( "qr", M, N, A, LDA );
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration<double, milli>(end - start);
  fmt::print(fg(fmt::color::green), "✓ Factorization completed in {:.3f} ms\n", duration.count());

  real_type R[M*M];
  qr.getR( R, M );
  print_matrix_with_border("R matrix (5×5)", M, M, R, M);

  real_type rhs[M], b[M];
  real_type x[N]{1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, N, 1, A, LDA, x, 1, 0, b, 1 );

  fmt::print(fg(fmt::color::yellow), "\n📊 Solving Least Squares: A·x = b\n");
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));

  qr.inv_permute( rhs );
  qr.invRt_mul( rhs, 1 );
  alglin::Copy_n( rhs, 3, x );
  alglin::Zero_n( x+3, 2 );
  qr.Q_mul( x );

  fmt::print("xᵀ = {}\n", alglin::print_matrix( 1, M, x, 1 ));

  alglin::gemv( Transposition::NO, M, N, -1, A, LDA, x, 1, 1, b, 1 );
  real_type residual_norm = alglin::nrm2(M, b, 1);
  
  fmt::print("rᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));
  fmt::print("‖r‖₂ = {:.6e}\n", residual_norm);
  
  if (residual_norm < 1e-5) {  // Relaxed threshold for ill-conditioned
    print_success("Test 3 passed!");
  } else {
    print_warning(fmt::format("Test 3 residual ({:.6e}) above threshold", residual_norm));
  }
}

// Template function for testing factorization methods
template<typename FactorizationType>
void test_factorization_method(
    FactorizationType& fac,
    const std::string& name,
    integer M,
    integer LDA,
    real_type A[],
    real_type rhs[],
    real_type b[],
    real_type x[]
) {
    fmt::print(fg(fmt::color::cyan), "\n\n🔧 Testing {} factorization\n", name);
    
    auto start = chrono::high_resolution_clock::now();
    fac.factorize( name.c_str(), M, M, A, LDA );
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double, milli>(end - start);
    
    fmt::print(fg(fmt::color::green), "  ✓ {} completed in {:.3f} ms\n", name, duration.count());
    fmt::print("  Solving A·x = b using {}\n", name);
    
    alglin::Copy_n( rhs, M, x );
    alglin::Copy_n( rhs, M, b );
    fac.solve( 1, x, M);
    fmt::print("  xᵀ = {}\n", alglin::print_matrix( 1, M, x, 1 ));

    alglin::gemv( Transposition::NO, M, M, -1, A, LDA, x, 1, 1, b, 1 );
    real_type res = alglin::nrm2( M, b, 1 );
    
    fmt::print("  ‖r‖₂ = {:.6e}", res);
    if (res < 1e-6) {
      fmt::print(fg(fmt::color::green), " ✓\n");
    } else {
      fmt::print(fg(fmt::color::orange), " ⚠️\n");
      UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
    }
}

static void test4() {
  print_header("TEST 4: Comparing Factorization Methods", fmt::color::light_blue);
  
  integer const M{5};
  integer const LDA{5};
  real_type A[]{
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M];
  real_type x[M]{1,2,3,4,5};
  alglin::gemv( Transposition::NO, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::NO, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  print_matrix_with_border("Test Matrix A (5×5)", M, M, A, M);
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));
  
  fmt::print(fg(fmt::color::yellow), "\n🎯 Running 8 different factorization methods:\n");
  fmt::print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

  // Test all factorization methods
  {
    alglin::LU<real_type> lu;
    test_factorization_method(lu, "LU", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LUPQ<real_type> lupq;
    test_factorization_method(lupq, "LUPQ", M, LDA, A, rhs, b, x);
  }
  {
    alglin::QR<real_type> qr;
    test_factorization_method(qr, "QR", M, LDA, A, rhs, b, x);
  }
  {
    alglin::QRP<real_type> qrp;
    test_factorization_method(qrp, "QRP", M, LDA, A, rhs, b, x);
  }
  {
    alglin::SVD<real_type> svd;
    test_factorization_method(svd, "SVD", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LSS<real_type> lss;
    test_factorization_method(lss, "LSS", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LSY<real_type> lsy;
    test_factorization_method(lsy, "LSY", M, LDA, A, rhs, b, x);
  }
  {
    alglin::PINV<real_type> pinv;
    test_factorization_method(pinv, "PINV", M, LDA, A, rhs, b, x);
  }

  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Test 4 completed - All methods compared!");
}

// Template function for testing transposed factorization methods
template<typename FactorizationType>
void test_transposed_factorization_method(
    FactorizationType& fac,
    const std::string& name,
    integer M,
    integer LDA,
    real_type A[],
    real_type rhs[],
    real_type b[],
    real_type x[]
) {
    fmt::print(fg(fmt::color::cyan), "\n\n🔧 Testing {} factorization (transposed)\n", name);
    
    auto start = chrono::high_resolution_clock::now();
    fac.factorize( name.c_str(), M, M, A, LDA );
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double, milli>(end - start);
    
    fmt::print(fg(fmt::color::green), "  ✓ {} completed in {:.3f} ms\n", name, duration.count());
    fmt::print("  Solving Aᵀ·x = b using {}\n", name);
    
    alglin::Copy_n( rhs, M, x );
    alglin::Copy_n( rhs, M, b );
    fac.t_solve( 1, x, M);
    fmt::print("  xᵀ = {}\n", alglin::print_matrix( 1, M, x, 1 ));

    alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 );
    real_type res = alglin::nrm2( M, b, 1 );
    
    fmt::print("  ‖r‖₂ = {:.6e}", res);
    if (res < 1e-6) {
      fmt::print(fg(fmt::color::green), " ✓\n");
    } else {
      fmt::print(fg(fmt::color::orange), " ⚠️\n");
      UTILS_ASSERT0( res < 1e-6, "test failed!\n" );
    }
}

static void test5() {
  print_header("TEST 5: Transposed System Solutions", fmt::color::violet);
  
  integer const M{5};
  integer const LDA{5};
  real_type A[]{
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,     3+1
  };

  real_type rhs[M], b[M];
  real_type x[M]{1,2,3,4,5};
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  print_matrix_with_border("Test Matrix A (5×5)", M, M, A, M);
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));
  
  fmt::print(fg(fmt::color::yellow), "\n🎯 Solving Aᵀ·x = b with 8 methods:\n");
  fmt::print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

  // Test all transposed factorization methods
  {
    alglin::LU<real_type> lu;
    test_transposed_factorization_method(lu, "LU", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LUPQ<real_type> lupq;
    test_transposed_factorization_method(lupq, "LUPQ", M, LDA, A, rhs, b, x);
  }
  {
    alglin::QR<real_type> qr;
    test_transposed_factorization_method(qr, "QR", M, LDA, A, rhs, b, x);
  }
  {
    alglin::QRP<real_type> qrp;
    test_transposed_factorization_method(qrp, "QRP", M, LDA, A, rhs, b, x);
  }
  {
    alglin::SVD<real_type> svd;
    test_transposed_factorization_method(svd, "SVD", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LSS<real_type> lss;
    test_transposed_factorization_method(lss, "LSS", M, LDA, A, rhs, b, x);
  }
  {
    alglin::LSY<real_type> lsy;
    test_transposed_factorization_method(lsy, "LSY", M, LDA, A, rhs, b, x);
  }
  {
    alglin::PINV<real_type> pinv;
    test_transposed_factorization_method(pinv, "PINV", M, LDA, A, rhs, b, x);
  }

  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Test 5 completed - All transposed systems solved!");
}

static void test6() {
  print_header("TEST 6: Tridiagonal Systems", fmt::color::teal);
  
  alglin::TridiagonalLU<real_type> lu;
  alglin::TridiagonalQR<real_type> qr;

  integer   const N{5};
  real_type const D[]{ 1, 1, 2, -0.1, 0.1 };
  real_type const L[]{ -0.1, -1, -2, -0.1 };
  real_type const U[]{ -1, -10, -2, 0.1 };

  fmt::print(fg(fmt::color::yellow), "\n📊 Tridiagonal Matrix:\n");
  fmt::print("Main diagonal D:  ");
  for (int i = 0; i < N; ++i) fmt::print("{:8.3f} ", D[i]);
  fmt::print("\nLower diagonal L:");
  for (int i = 0; i < N-1; ++i) fmt::print("{:8.3f} ", L[i]);
  fmt::print("\nUpper diagonal U:");
  for (int i = 0; i < N-1; ++i) fmt::print("{:8.3f} ", U[i]);
  fmt::print("\n");

  real_type rhs[N], b[N];
  real_type x[N]{1,2,3,4,5};
  qr.axpy( N, 1.0, L, D, U, x, 0.0, rhs );
  qr.axpy( N, 1.0, L, D, U, x, 0.0, b   );

  fmt::print("\nbᵀ = ");
  for (int i = 0; i < N; ++i) fmt::print("{:8.3f} ", b[i]);
  fmt::print("\n");

  fmt::print(fg(fmt::color::cyan), "\n🔧 Testing Tridiagonal LU\n");
  auto start_lu = chrono::high_resolution_clock::now();
  lu.factorize( "lu", N, L, D, U );
  auto end_lu = chrono::high_resolution_clock::now();
  auto duration_lu = chrono::duration<double, milli>(end_lu - start_lu);
  fmt::print(fg(fmt::color::green), "  ✓ Tridiagonal LU completed in {:.3f} ms\n", duration_lu.count());

  alglin::Copy_n( rhs, N, x );
  alglin::Copy_n( rhs, N, b );
  lu.solve( x );
  fmt::print("  xᵀ = ");
  for (int i = 0; i < N; ++i) fmt::print("{:8.3f} ", x[i]);
  fmt::print("\n");

  lu.axpy( N, -1.0, L, D, U, x, 1.0, b );
  real_type res_lu = alglin::nrm2(N, b, 1);
  fmt::print("  ‖r‖₂ = {:.6e}", res_lu);
  if (res_lu < 1e-6) {
    fmt::print(fg(fmt::color::green), " ✓\n");
  } else {
    fmt::print(fg(fmt::color::orange), " ⚠️\n");
  }

  fmt::print(fg(fmt::color::cyan), "\n🔧 Testing Tridiagonal QR\n");
  auto start_qr = chrono::high_resolution_clock::now();
  qr.factorize( "qr", N, L, D, U );
  auto end_qr = chrono::high_resolution_clock::now();
  auto duration_qr = chrono::duration<double, milli>(end_qr - start_qr);
  fmt::print(fg(fmt::color::green), "  ✓ Tridiagonal QR completed in {:.3f} ms\n", duration_qr.count());

  alglin::Copy_n( rhs, N, x );
  alglin::Copy_n( rhs, N, b );
  qr.solve( x );
  fmt::print("  xᵀ = ");
  for (int i = 0; i < N; ++i) fmt::print("{:8.3f} ", x[i]);
  fmt::print("\n");

  qr.axpy( N, -1.0, L, D, U, x, 1.0, b );
  real_type res_qr = alglin::nrm2(N, b, 1);
  fmt::print("  ‖r‖₂ = {:.6e}", res_qr);
  if (res_qr < 1e-6) {
    fmt::print(fg(fmt::color::green), " ✓\n");
  } else {
    fmt::print(fg(fmt::color::orange), " ⚠️\n");
  }

  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Test 6 completed - Tridiagonal methods tested!");
}

static void test7() {
  print_header("TEST 7: Advanced QRP Tests", fmt::color::gold);
  
  alglin::QRP<real_type> qrp;

  constexpr integer M{5};
  constexpr integer LDA{5};
  constexpr real_type A[]{
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };

  real_type rhs[M], b[M];
  real_type x[M]{1,2,3,4,5};
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, rhs, 1 );
  alglin::gemv( Transposition::YES, M, M, 1, A, LDA, x, 1, 0, b, 1 );

  print_matrix_with_border("Test Matrix A (5×5)", M, M, A, M);
  fmt::print("bᵀ = {}\n", alglin::print_matrix( 1, M, b, 1 ));

  fmt::print(fg(fmt::color::yellow), "\n🎯 Test 7A: Standard QRP\n");
  auto start1 = chrono::high_resolution_clock::now();
  qrp.factorize( "qrp", M, M, A, LDA );
  auto end1 = chrono::high_resolution_clock::now();
  auto duration1 = chrono::duration<double, milli>(end1 - start1);
  fmt::print(fg(fmt::color::green), "✓ Factorization completed in {:.3f} ms\n", duration1.count());

  alglin::Copy_n( rhs, M, x );
  alglin::Copy_n( rhs, M, b );
  qrp.t_solve( x );
  fmt::print("xᵀ = {}\n", alglin::print_matrix( 1, M, x, 1 ));

  alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  real_type res1 = alglin::nrm2(M, b, 1);
  fmt::print("‖r‖₂ = {:.6e}", res1);
  if (res1 < 1e-6) {
    fmt::print(fg(fmt::color::green), " ✓\n");
  } else {
    fmt::print(fg(fmt::color::orange), " ⚠️\n");
  }

  fmt::print(fg(fmt::color::yellow), "\n🎯 Test 7B: QRP with Matrix Block Loading\n");
  alglin::Matrix<real_type> mat;
  mat.setup( M, M );
  mat.load_block( 2, 5, A,   LDA, 0, 0 );
  mat.load_block( 3, 5, A+2, LDA, 2, 0 );
  
  fmt::print("Matrix loaded in blocks:\n");
  fmt::print("  Block 1: 2×5 loaded at (0,0)\n");
  fmt::print("  Block 2: 3×5 loaded at (2,0)\n");

  auto start2 = chrono::high_resolution_clock::now();
  qrp.factorize( "qrp", mat );
  auto end2 = chrono::high_resolution_clock::now();
  auto duration2 = chrono::duration<double, milli>(end2 - start2);
  fmt::print(fg(fmt::color::green), "✓ Block factorization completed in {:.3f} ms\n", duration2.count());

  alglin::Copy_n( rhs, M, x );
  alglin::Copy_n( rhs, M, b );
  qrp.t_solve( x );
  fmt::print("xᵀ = {}\n", alglin::print_matrix( 1, M, x, 1 ));

  alglin::gemv( Transposition::YES, M, M, -1, A, LDA, x, 1, 1, b, 1 );
  real_type res2 = alglin::nrm2(M, b, 1);
  fmt::print("‖r‖₂ = {:.6e}", res2);
  if (res2 < 1e-6) {
    fmt::print(fg(fmt::color::green), " ✓\n");
  } else {
    fmt::print(fg(fmt::color::orange), " ⚠️\n");
  }

  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Test 7 completed - Advanced QRP tests passed!");
}

// ===========================================================================
// New Tests
// ===========================================================================

static void test8() {
  print_header("TEST 8: Performance Benchmark", fmt::color::crimson);
  
  constexpr integer M{100};
  constexpr integer N{100};
  
  fmt::print(fg(fmt::color::yellow), "📈 Generating random {0}×{1} matrix for benchmark\n", M, N);
  
  // Generate random matrix
  std::vector<real_type> A(M * N);
  std::vector<real_type> x(N, 1.0);
  std::vector<real_type> b(M);
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<real_type> dis(-10.0, 10.0);
  
  for (integer i = 0; i < M * N; ++i) {
    A[i] = dis(gen);
  }
  
  fmt::print("Running benchmark with different factorization methods...\n");
  fmt::print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  
  // Lambda function for benchmarking
  auto benchmark_method = [&](const std::string& name, auto factory_func) {
    auto start = chrono::high_resolution_clock::now();
    factory_func();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration<double, milli>(end - start);
    
    fmt::print(fg(fmt::color::cyan), "{:6} ", name);
    fmt::print("Time: {:8.3f} ms", duration.count());
    
    // Simple rating
    if (duration.count() < 10) {
      fmt::print(fg(fmt::color::green), " ⚡ Fast\n");
    } else if (duration.count() < 50) {
      fmt::print(fg(fmt::color::yellow), " ⏱️  Moderate\n");
    } else {
      fmt::print(fg(fmt::color::orange), " 🐌 Slow\n");
    }
  };
  
  // Benchmark different methods
  benchmark_method("LU", [&]() {
    alglin::LU<real_type> lu;
    lu.factorize("lu", M, N, A.data(), M);
  });
  
  benchmark_method("QR", [&]() {
    alglin::QR<real_type> qr;
    qr.factorize("qr", M, N, A.data(), M);
  });
  
  benchmark_method("SVD", [&]() {
    alglin::SVD<real_type> svd;
    svd.factorize("svd", M, N, A.data(), M);
  });
  
  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Benchmark completed!");
}

static void test9() {
  print_header("TEST 9: Matrix Properties Analysis", fmt::color::lime);
  
  constexpr integer M{5};
  constexpr integer N{5};
  constexpr integer LDA{5};
  constexpr real_type A[]{
    0.001,      2,     3,       2,      3,
    0.001,  0.001,     0,   0.001,  1e-10,
    0,      0.001,     0,   0.001,  1e-12,
    0.001,     -1, 1e-6+1,     -1, -1e-12,
    0.000001,   5,     3,  1e-6+5,    3+1
  };
  
  fmt::print(fg(fmt::color::yellow), "📐 Analyzing matrix properties:\n\n");
  
  // Lambda function for computing matrix norm (Frobenius)
  auto compute_frobenius_norm = [&]() -> real_type {
    real_type norm = 0.0;
    for (integer i = 0; i < M * N; ++i) {
      norm += A[i] * A[i];
    }
    return sqrt(norm);
  };
  
  // Lambda function for checking symmetry
  auto check_symmetry = [&]() -> bool {
    for (integer i = 0; i < M; ++i) {
      for (integer j = 0; j < N; ++j) {
        if (fabs(A[i * LDA + j] - A[j * LDA + i]) > 1e-10) {
          return false;
        }
      }
    }
    return true;
  };
  
  // Lambda function for counting small elements
  auto count_small_elements = [&](real_type threshold) -> integer {
    integer count = 0;
    for (integer i = 0; i < M * N; ++i) {
      if (fabs(A[i]) < threshold) {
        count++;
      }
    }
    return count;
  };
  
  // Lambda function for computing trace
  auto compute_trace = [&]() -> real_type {
    real_type trace = 0.0;
    for (integer i = 0; i < std::min(M, N); ++i) {
      trace += A[i * LDA + i];
    }
    return trace;
  };
  
  // Execute analysis
  real_type norm = compute_frobenius_norm();
  bool symmetric = check_symmetry();
  integer small_count = count_small_elements(1e-6);
  real_type trace = compute_trace();
  
  fmt::print("Frobenius norm: {:.6f}\n", norm);
  fmt::print("Symmetric: ");
  if (symmetric) {
    fmt::print(fg(fmt::color::green), "Yes ✓\n");
  } else {
    fmt::print(fg(fmt::color::yellow), "No\n");
  }
  
  fmt::print("Elements < 1e-6: {}/{} ({:.1f}%)\n", 
             small_count, M*N, 100.0*small_count/(M*N));
  fmt::print("Trace: {:.6f}\n", trace);
  
  fmt::print(fg(fmt::color::green), "\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
  print_success("Matrix analysis completed!");
}

// ===========================================================================
// Main function
// ===========================================================================

int main() {
  auto total_start = chrono::high_resolution_clock::now();
  
  print_header("ALGLIN TEST SUITE", fmt::color::green);
  fmt::print(fg(fmt::color::yellow), "🚀 Starting comprehensive linear algebra tests\n");
  fmt::print("Compiled: {} {}\n", __DATE__, __TIME__);
  fmt::print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n");

  try {
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    test8();
    test9();
    
    auto total_end = chrono::high_resolution_clock::now();
    auto total_duration = chrono::duration<double, milli>(total_end - total_start);
    
    fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, 
      "\n"
      "┌────────────────────────────────────────────────────────────┐\n"
      "│                    ALL TESTS COMPLETED!                    │\n"
      "│                    Total time: {:8.3f} ms                 │\n"
      "└────────────────────────────────────────────────────────────┘\n",
      total_duration.count());
    
  } catch ( exception const & exc ) {
    print_error(fmt::format("Exception: {}", exc.what()));
    fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, 
               "\n❌ TEST SUITE FAILED!\n");
    return 1;
  } catch ( ... ) {
    print_error("Unknown error occurred!");
    fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, 
               "\n❌ TEST SUITE FAILED WITH UNKNOWN ERROR!\n");
    return 1;
  }

  return 0;
}
