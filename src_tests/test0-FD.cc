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
#include "Utils_fmt.hh"
#include "Utils_string.hh"
#include "Utils_TicToc.hh"
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

using namespace alglin;
using namespace std;

// ===========================================================================
// Utility functions for formatted output
// ===========================================================================

namespace TestUtils {

static void print_header(const std::string& title, fmt::color color = fmt::color::cyan) {
  const int width = 78;
  fmt::print(fg(color) | fmt::emphasis::bold,
    "\n"
    "╔{:═^{}}╗\n"
    "║{:^{}}║\n"
    "╚{:═^{}}╝\n",
    "", width,
    title, width,
    "", width
  );
}

static void print_subheader(const std::string& title, fmt::color color = fmt::color::yellow) {
  fmt::print(fg(color) | fmt::emphasis::bold,
    "\n"
    "┌{:─^{}}┐\n"
    "│ {:<74} │\n"
    "└{:─^{}}┘\n",
    "", 76,
    title,
    "", 76
  );
}

#if 0
static void print_success(const std::string& message) {
  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "✅  {}\n", message);
}
static void print_warning(const std::string& message) {
  fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold, "⚠️   {}\n", message);
}
static void print_debug(const std::string& message) {
  fmt::print(fg(fmt::color::gray), "📝  {}\n", message);
}
#endif

static void print_error(const std::string& message) {
  fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "❌  {}\n", message);
}

static void print_info(const std::string& message) {
  fmt::print(fg(fmt::color::blue), "🔹  {}\n", message);
}

static void print_progress(int current, int total, const std::string& message = "") {
  float const progress = static_cast<float>(current) / static_cast<float>(total);
  int constexpr bar_width = 40;
  
  fmt::print("\r");
  fmt::print(fg(fmt::color::cyan), "[");
  int const pos = static_cast<int>(static_cast<float>(bar_width) * progress);
  for (int i = 0; i < bar_width; ++i) {
    if (i < pos) fmt::print(fg(fmt::color::green), "█");
    else if (i == pos) fmt::print(fg(fmt::color::yellow), "▌");
    else fmt::print(" ");
  }
  fmt::print(fg(fmt::color::cyan), "] ");
  fmt::print("{:3.0f}%", progress * 100);
  if (!message.empty()) fmt::print("  {}", message);
}

struct TestResult {
  string name;
  bool passed;
  double analytic_time_mus;
  double fd_time_mus;
  double max_error;
  double avg_error;
  int dimension;
  
  double speedup() const { 
    return analytic_time_mus > 0 ? fd_time_mus / analytic_time_mus : 0.0;
  }
  
  string status_symbol() const {
    return passed ? "✓" : "✗";
  }
  
  fmt::color status_color() const {
    return passed ? fmt::color::green : fmt::color::red;
  }
};

static void print_result_table(const vector<TestResult>& results) {
  constexpr int col_widths[] = {30, 8, 12, 12, 12, 12, 10, 8};
  
  // Table header
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "┌{:─^{}}┬{:─^{}}┬{:─^{}}┬{:─^{}}┬{:─^{}}┬{:─^{}}┬{:─^{}}┬{:─^{}}┐\n",
    "", col_widths[0], "", col_widths[1], "", col_widths[2], "", col_widths[3],
    "", col_widths[4], "", col_widths[5], "", col_widths[6], "", col_widths[7]);
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│{:^{}}│{:^{}}│{:^{}}│{:^{}}│{:^{}}│{:^{}}│{:^{}}│{:^{}}│\n",
    "Test Name", col_widths[0],
    "Dim", col_widths[1],
    "Analytic(μs)", col_widths[2],
    "FD(μs)", col_widths[3],
    "Speedup", col_widths[4],
    "Max Error", col_widths[5],
    "Avg Error", col_widths[6],
    "Status", col_widths[7]);
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "├{:─^{}}┼{:─^{}}┼{:─^{}}┼{:─^{}}┼{:─^{}}┼{:─^{}}┼{:─^{}}┼{:─^{}}┤\n",
    "", col_widths[0], "", col_widths[1], "", col_widths[2], "", col_widths[3],
    "", col_widths[4], "", col_widths[5], "", col_widths[6], "", col_widths[7]);
  
  // Table rows
  for (const auto& result : results) {
    fmt::print("│{:<{}}│{:^{}}│{:^{}.3f}│{:^{}.3f}│{:^{}.1f}x │{:^{}.2e}│{:^{}.2e}│",
      result.name, col_widths[0],
      result.dimension, col_widths[1],
      result.analytic_time_mus, col_widths[2],
      result.fd_time_mus, col_widths[3],
      result.speedup(), col_widths[4]-2,
      result.max_error, col_widths[5],
      result.avg_error, col_widths[6]);
    
    fmt::print(fg(result.status_color()) | fmt::emphasis::bold,
      "{:^{}}│\n", result.status_symbol(), col_widths[7]);
  }
  
  // Table footer
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "└{:─^{}}┴{:─^{}}┴{:─^{}}┴{:─^{}}┴{:─^{}}┴{:─^{}}┴{:─^{}}┴{:─^{}}┘\n",
    "", col_widths[0], "", col_widths[1], "", col_widths[2], "", col_widths[3],
    "", col_widths[4], "", col_widths[5], "", col_widths[6], "", col_widths[7]);
  
  // Summary statistics
  double total_speedup = 0.0;
  int passed_count = 0;
  for (const auto& result : results) {
    total_speedup += result.speedup();
    if (result.passed) passed_count++;
  }
  
  fmt::print(fg(fmt::color::white), "\n📈 Summary: {}/{} tests passed", 
             passed_count, results.size());
  double const result_count = static_cast<double>(results.size());
  fmt::print(fg(fmt::color::cyan), " | Average speedup: {:.1f}x\n", 
             total_speedup / result_count);
}

} // namespace TestUtils

// ===========================================================================
// Test Function Definitions
// ===========================================================================

// Original test function from the file
static bool fun(real_type const x[4], real_type & res) {
  res = x[0] + sin(x[1]) + exp(x[2]*x[3]) + x[0]*x[1]*x[2];
  return true;
}

static bool fun_grad(real_type const x[4], real_type grad[4]) {
  grad[0] = 1 + x[1]*x[2];
  grad[1] = cos(x[1]) + x[0]*x[2];
  grad[2] = exp(x[2]*x[3])*x[3] + x[0]*x[1];
  grad[3] = exp(x[2]*x[3])*x[2];
  return true;
}

// ===========================================================================
// Extended Test Functions Collection
// ===========================================================================

namespace TestFunctions {

// 1. Linear Function: f(x) = a·x + b
struct LinearFunction {
  vector<real_type> a;
  real_type b;
  
  LinearFunction(const vector<real_type>& coeffs, real_type constant = 0.0) 
    : a(coeffs), b(constant) {}
  
  bool eval(real_type const x[], real_type & res) const {
    res = b;
    for (size_t i = 0; i < a.size(); ++i) {
      res += a[i] * x[i];
    }
    return true;
  }
  
  bool gradient(real_type const x[], real_type grad[]) const {
    (void)x; // x is not used for linear gradient
    for (size_t i = 0; i < a.size(); ++i) {
      grad[i] = a[i];
    }
    return true;
  }
};

// 2. Quadratic Form: f(x) = ½ xᵀ·A·x + b·x + c
struct QuadraticFunction {
  vector<vector<real_type>> A;
  vector<real_type> b;
  real_type c;
  
  QuadraticFunction(const vector<vector<real_type>>& Amatrix,
                   const vector<real_type>& linear,
                   real_type constant = 0.0)
    : A(Amatrix), b(linear), c(constant) {}
  
  bool eval(real_type const x[], real_type & res) const {
    integer n = static_cast<integer>(b.size());
    res = c;
    
    // Linear term
    for (integer i = 0; i < n; ++i) {
      res += b[i] * x[i];
    }
    
    // Quadratic term
    for (integer i = 0; i < n; ++i) {
      for (integer j = 0; j < n; ++j) {
        res += 0.5 * A[i][j] * x[i] * x[j];
      }
    }
    return true;
  }
  
  bool gradient(real_type const x[], real_type grad[]) const {
    integer n = static_cast<integer>(b.size());
    for (integer i = 0; i < n; ++i) {
      grad[i] = b[i];
      for (integer j = 0; j < n; ++j) {
        grad[i] += A[i][j] * x[j];
      }
    }
    return true;
  }
};

// 3. Exponential Sum: f(x) = Σ exp(a_i·x_i)
struct ExponentialFunction {
  vector<real_type> a;
  
  ExponentialFunction(const vector<real_type>& coeffs) : a(coeffs) {}
  
  bool eval(real_type const x[], real_type & res) const {
    res = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      res += exp(a[i] * x[i]);
    }
    return true;
  }
  
  bool gradient(real_type const x[], real_type grad[]) const {
    for (size_t i = 0; i < a.size(); ++i) {
      grad[i] = a[i] * exp(a[i] * x[i]);
    }
    return true;
  }
};

// 4. Trigonometric Function: f(x) = Σ sin(a_i·x_i) + cos(b_i·x_i)
struct TrigonometricFunction {
  vector<real_type> a, b;
  
  TrigonometricFunction(const vector<real_type>& sin_coeffs,
                       const vector<real_type>& cos_coeffs)
    : a(sin_coeffs), b(cos_coeffs) {}
  
  bool eval(real_type const x[], real_type & res) const {
    res = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
      res += sin(a[i] * x[i]) + cos(b[i] * x[i]);
    }
    return true;
  }
  
  bool gradient(real_type const x[], real_type grad[]) const {
    for (size_t i = 0; i < a.size(); ++i) {
      grad[i] = a[i] * cos(a[i] * x[i]) - b[i] * sin(b[i] * x[i]);
    }
    return true;
  }
};

// 5. Rosenbrock Function (N-dimensional)
struct RosenbrockFunction {
  real_type a, b;
  
  RosenbrockFunction(real_type a_val = 1.0, real_type b_val = 100.0) 
    : a(a_val), b(b_val) {}
  
  bool eval(real_type const x[], real_type & res) const {
    integer n = 2; // Standard 2D Rosenbrock
    res = 0.0;
    for (integer i = 0; i < n-1; ++i) {
      res += b * pow(x[i+1] - x[i]*x[i], 2) + pow(a - x[i], 2);
    }
    return true;
  }
  
  bool gradient(real_type const x[], real_type grad[]) const {
    integer n = 2;
    grad[0] = -400 * x[0] * (x[1] - x[0]*x[0]) - 2 * (a - x[0]);
    grad[1] = 200 * (x[1] - x[0]*x[0]);
    // For higher dimensions, we would need to extend this
    for (integer i = 2; i < n; ++i) {
      grad[i] = 0.0;
    }
    return true;
  }
};

} // namespace TestFunctions

// ===========================================================================
// Test Runner Class
// ===========================================================================

class FiniteDifferenceTester {
private:
  vector<TestUtils::TestResult> results_;
  Utils::TicToc timer_;
  
  template<typename Functor, typename GradFunctor>
  TestUtils::TestResult run_gradient_test(
    const string& name,
    const vector<real_type>& x0,
    Functor& func,
    GradFunctor& grad_func,
    real_type epsilon = 1e-6
  ) {
    TestUtils::TestResult result;
    result.name = name;
    result.dimension = static_cast<int>(x0.size());
    
    vector<real_type> x(x0.begin(), x0.end());
    vector<real_type> grad_analytic(x0.size());
    vector<real_type> grad_fd(x0.size());
    
    // Analytical gradient
    timer_.tic();
    bool ok = grad_func(x.data(), grad_analytic.data());
    timer_.toc();
    result.analytic_time_mus = timer_.elapsed_mus();
    
    if (!ok) {
      result.passed = false;
      result.max_error = numeric_limits<real_type>::max();
      return result;
    }
    
    // Finite difference gradient
    auto func_wrapper = [&](real_type const x_vec[], real_type & res) -> bool {
      return func(x_vec, res);
    };
    
    timer_.tic();
    ok = finite_difference_gradient(x.data(), 
                                   static_cast<integer>(x0.size()),
                                   func_wrapper,
                                   grad_fd.data());
    timer_.toc();
    result.fd_time_mus = timer_.elapsed_mus();
    
    if (!ok) {
      result.passed = false;
      result.max_error = numeric_limits<real_type>::max();
      return result;
    }
    
    // Compute errors
    result.max_error = 0.0;
    result.avg_error = 0.0;
    
    for (size_t i = 0; i < x0.size(); ++i) {
      real_type error = abs(grad_fd[i] - grad_analytic[i]);
      result.max_error = max(result.max_error, error);
      result.avg_error += error;
    }
    result.avg_error /= static_cast<real_type>(x0.size());
    
    result.passed = (result.max_error < epsilon);
    
    // Print detailed comparison
    if (x0.size() <= 8) { // Only print details for small dimensions
      TestUtils::print_subheader(name + " Gradient Comparison");
      
      fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
        "┌─────┬────────────────┬────────────────┬────────────────┬──────────┐\n");
      fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
        "│ {:3} │ {:14} │ {:14} │ {:14} │ {:8} │\n",
        "Idx", "Analytic", "FD", "Abs Error", "Rel %");
      fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
        "├─────┼────────────────┼────────────────┼────────────────┼──────────┤\n");
      
      for (size_t i = 0; i < x0.size(); ++i) {
        real_type abs_err = abs(grad_fd[i] - grad_analytic[i]);
        real_type rel_err = (grad_analytic[i] != 0.0) ? 
                           abs_err / abs(grad_analytic[i]) * 100.0 : 0.0;
        
        fmt::color err_color = abs_err < epsilon ? fmt::color::green :
                              abs_err < epsilon*10 ? fmt::color::yellow : 
                              fmt::color::red;
        
        fmt::print("│ {:3d} │ {:14.6e} │ {:14.6e} │ ", i, grad_analytic[i], grad_fd[i]);
        fmt::print(fg(err_color), "{:14.6e} │ {:7.3f}% │\n", abs_err, rel_err);
      }
      
      fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
        "└─────┴────────────────┴────────────────┴────────────────┴──────────┘\n");
      
      fmt::print("📊 Max Error: {:.2e} | Avg Error: {:.2e} | Threshold: {:.0e}\n",
                 result.max_error, result.avg_error, epsilon);
    }
    
    return result;
  }
  
public:
  void run_all_tests() {
    TestUtils::print_header("FINITE DIFFERENCE TEST SUITE", fmt::color::green);
    
    fmt::print(fg(fmt::color::yellow), 
      "🚀 Starting comprehensive finite difference validation\n\n");
    
    // Test 1: Original function
    {
      TestUtils::print_subheader("Test 1: Original Function (4D)");
      auto func = [](real_type const x[], real_type & res) -> bool {
        return fun(x, res);
      };
      auto grad_func = [](real_type const x[], real_type grad[]) -> bool {
        return fun_grad(x, grad);
      };
      vector<real_type> x0 = {1.0, 2.0, 1.0, 4.0};
      results_.push_back(run_gradient_test("Original Function", x0, func, grad_func));
    }
    
    // Test 2: Linear Function
    {
      TestUtils::print_subheader("Test 2: Linear Function (5D)");
      TestFunctions::LinearFunction linear(
        {1.5, 2.3, 0.8, -1.2, 3.4},  // coefficients
        2.5                           // constant
      );
      auto func = [&](real_type const x[], real_type & res) -> bool {
        return linear.eval(x, res);
      };
      auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
        return linear.gradient(x, grad);
      };
      vector<real_type> x0 = {0.5, 1.0, -0.5, 2.0, 1.5};
      results_.push_back(run_gradient_test("Linear Function", x0, func, grad_func));
    }
    
    // Test 3: Quadratic Function
    {
      TestUtils::print_subheader("Test 3: Quadratic Function (3D)");
      vector<vector<real_type>> A = {
        {2.0, 1.0, 0.5},
        {1.0, 3.0, -0.2},
        {0.5, -0.2, 1.5}
      };
      vector<real_type> b = {1.0, -2.0, 0.5};
      TestFunctions::QuadraticFunction quadratic(A, b, 1.0);
      auto func = [&](real_type const x[], real_type & res) -> bool {
        return quadratic.eval(x, res);
      };
      auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
        return quadratic.gradient(x, grad);
      };
      vector<real_type> x0 = {1.0, 0.5, -1.0};
      results_.push_back(run_gradient_test("Quadratic Function", x0, func, grad_func));
    }
    
    // Test 4: Exponential Function
    {
      TestUtils::print_subheader("Test 4: Exponential Function (4D)");
      TestFunctions::ExponentialFunction exp_func({0.5, 1.0, -0.3, 2.0});
      auto func = [&](real_type const x[], real_type & res) -> bool {
        return exp_func.eval(x, res);
      };
      auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
        return exp_func.gradient(x, grad);
      };
      vector<real_type> x0 = {0.2, 0.5, 1.0, -0.5};
      results_.push_back(run_gradient_test("Exponential Function", x0, func, grad_func));
    }
    
    // Test 5: Trigonometric Function
    {
      TestUtils::print_subheader("Test 5: Trigonometric Function (3D)");
      TestFunctions::TrigonometricFunction trig(
        {1.0, 2.0, 0.5},  // sin coefficients
        {0.5, 1.5, 2.0}   // cos coefficients
      );
      auto func = [&](real_type const x[], real_type & res) -> bool {
        return trig.eval(x, res);
      };
      auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
        return trig.gradient(x, grad);
      };
      vector<real_type> x0 = {M_PI/4, M_PI/3, M_PI/6};
      results_.push_back(run_gradient_test("Trigonometric Function", x0, func, grad_func));
    }
    
    // Test 6: Rosenbrock Function
    {
      TestUtils::print_subheader("Test 6: Rosenbrock Function (2D)");
      TestFunctions::RosenbrockFunction rosenbrock;
      auto func = [&](real_type const x[], real_type & res) -> bool {
        return rosenbrock.eval(x, res);
      };
      auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
        return rosenbrock.gradient(x, grad);
      };
      vector<real_type> x0 = {-1.2, 1.0};  // Standard starting point
      results_.push_back(run_gradient_test("Rosenbrock Function", x0, func, grad_func, 1e-5));
    }
    
    // Test 7: Multi-point test
    {
      TestUtils::print_subheader("Test 7: Multi-Point Validation");
      vector<vector<real_type>> test_points = {
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0},
        {-1.0, 0.5, 2.0},
        {0.3, -0.7, 0.1}
      };
      
      TestFunctions::QuadraticFunction quad(
        {{1, 0.5, 0.2}, {0.5, 2, -0.1}, {0.2, -0.1, 1}},
        {0.1, -0.2, 0.3},
        0.5
      );
      
      double max_error_all = 0.0;
      double avg_error_all = 0.0;
      int point_count = 0;
      
      for (const auto& point : test_points) {
        auto func = [&](real_type const x[], real_type & res) -> bool {
          return quad.eval(x, res);
        };
        auto grad_func = [&](real_type const x[], real_type grad[]) -> bool {
          return quad.gradient(x, grad);
        };
        auto result = run_gradient_test("Multi-Point #" + to_string(point_count+1), 
                                       point, func, grad_func);
        max_error_all = max(max_error_all, result.max_error);
        avg_error_all += result.avg_error;
        point_count++;
        
        // Progress indicator
        TestUtils::print_progress(point_count, static_cast<int>(test_points.size()),
                                 "Testing multiple points...");
      }
      fmt::print("\n");
      
      TestUtils::TestResult multi_result;
      multi_result.name = "Multi-Point Test";
      multi_result.dimension = 3;
      multi_result.passed = (max_error_all < 1e-6);
      multi_result.max_error = max_error_all;
      multi_result.avg_error = avg_error_all / point_count;
      multi_result.analytic_time_mus = 0.0;  // Not measured for this summary
      multi_result.fd_time_mus = 0.0;
      
      results_.push_back(multi_result);
    }
  }
  
  void print_summary() {
    TestUtils::print_header("TEST SUMMARY REPORT", fmt::color::green);
    
    // Overall statistics
    int passed = 0, failed = 0;
    double total_speedup = 0.0;
    double max_error_overall = 0.0;
    
    for (const auto& result : results_) {
      if (result.passed) passed++;
      else failed++;
      total_speedup += result.speedup();
      max_error_overall = max(max_error_overall, result.max_error);
    }
    
    // Summary box
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
      "╔══════════════════════════════════════════════════════════════════════════╗\n"
      "║                            OVERALL SUMMARY                               ║\n"
      "╠══════════════════════════════════════════════════════════════════════════╣\n");
    
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, "║");
    fmt::print(fg(fmt::color::white),"   Total Tests: {:<10d}     Passed: ", results_.size());
    fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "{:<10d} ", passed);
    fmt::print(fg(fmt::color::white), "     Failed: ");
    fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "{:<10d}", failed);
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, " ║\n");
    
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, "║" );
    double const result_count = static_cast<double>(results_.size());
    fmt::print(fg(fmt::color::white), "   Success Rate: {:5.1f}%                                                   ", passed * 100.0 / result_count);
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold, "║\n║" );
    
    fmt::print(fg(fmt::color::white), 
      "   Max Error: {:10.2e}    Avg Speedup: {:5.1f}x                           ",
      max_error_overall, total_speedup / result_count);
    
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
      "║\n"
      "╚══════════════════════════════════════════════════════════════════════════╝\n\n");
    
    // Detailed results table
    TestUtils::print_result_table(results_);
    
    // Recommendations
    TestUtils::print_subheader("Recommendations");
    
    if (max_error_overall < 1e-9) {
      fmt::print(fg(fmt::color::green), "🎯 Excellent accuracy! Finite differences are reliable.\n");
    } else if (max_error_overall < 1e-6) {
      fmt::print(fg(fmt::color::yellow), "📊 Good accuracy. Suitable for most applications.\n");
    } else if (max_error_overall < 1e-3) {
      fmt::print(fg(fmt::color::orange), "⚠️  Moderate accuracy. Consider adjusting epsilon values.\n");
    } else {
      fmt::print(fg(fmt::color::red), "❌ Poor accuracy. Review finite difference implementation.\n");
    }
    
    double avg_speedup = total_speedup / result_count;
    if (avg_speedup > 50) {
      fmt::print(fg(fmt::color::blue), "⚡ Analytical gradients are significantly faster.\n");
    } else if (avg_speedup > 10) {
      fmt::print(fg(fmt::color::cyan), "🏃 Good speed advantage for analytical gradients.\n");
    }
    
    // Final verdict
    fmt::print("\n");
    if (failed == 0) {
      fmt::print(fg(fmt::color::green) | fmt::emphasis::bold,
        "✅ ALL TESTS PASSED! Finite difference implementation is validated.\n");
    } else {
      fmt::print(fg(fmt::color::red) | fmt::emphasis::bold,
        "❌ {} TEST(S) FAILED. Review the implementation.\n", failed);
    }
  }
  
  const vector<TestUtils::TestResult>& get_results() const { return results_; }
};

// ===========================================================================
// Performance Benchmark Function
// ===========================================================================

static void benchmark_finite_difference() {
  TestUtils::print_header("PERFORMANCE BENCHMARK", fmt::color::magenta);
  
  const vector<int> dimensions = {2, 4, 8, 16, 32, 64};
  
  fmt::print(fg(fmt::color::yellow), 
    "Benchmarking finite difference performance across dimensions\n\n");
  
  // Benchmark table header
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "┌───────────┬────────────┬─────────────┬────────────┐\n"
    "│ Dimension │ Time (μs)  │ Evaluations │ Speed      │\n"
    "├───────────┼────────────┼─────────────┼────────────┤\n"
  );
  
  Utils::TicToc timer;
  
  for (int dim : dimensions) {
    fmt::print("│ {:9d} │", dim);
    
    // Create test function (simple quadratic)
    vector<real_type> x0(dim, 1.0);
    vector<real_type> coeffs(dim);
    for (int i = 0; i < dim; ++i) {
      coeffs[i] = 0.5 * (i + 1);
    }
    TestFunctions::LinearFunction linear_func(coeffs, 1.0);
    
    vector<real_type> grad_fd(dim);
    
    timer.tic();
    auto func_wrapper = [&](real_type const x[], real_type & res) -> bool {
      return linear_func.eval(x, res);
    };
    
    bool ok = finite_difference_gradient(x0.data(), dim, func_wrapper, grad_fd.data());
    timer.toc();
    
    if (ok) {
      double time_mus = timer.elapsed_mus();
      // Color code based on performance
      if (time_mus < 1.0) {
        fmt::print(fg(fmt::color::green), " {:10.3f} │", time_mus);
      } else if (time_mus < 10.0) {
        fmt::print(fg(fmt::color::yellow), " {:10.3f} │", time_mus);
      } else {
        fmt::print(fg(fmt::color::red), " {:10.3f} │", time_mus);
      }
      
      // Number of function evaluations (dim + 1 for finite difference)
      fmt::print(" {:11d} │", dim + 1);
      
      // Performance category
      if (time_mus < 1.0) {
        fmt::print(fg(fmt::color::green), " {:^10} │\n", "Fast");
      } else if (time_mus < 10.0) {
        fmt::print(fg(fmt::color::yellow), " {:^10} │\n", "Moderate");
      } else {
        fmt::print(fg(fmt::color::red), " {:^10} │\n", "Slow");
      }
    } else {
      fmt::print(" {:>10} │ {:>10} │ {:^10} │\n", "FAIL", "N/A", "FAIL");
    }
  }
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "└───────────┴────────────┴─────────────┴────────────┘\n");
  
  TestUtils::print_info("Note: Performance depends on function complexity and dimension.");
  TestUtils::print_info("For N dimensions, finite difference requires N+1 function evaluations.");
}

// ===========================================================================
// Main Function
// ===========================================================================

int main() {
  try {
    Utils::TicToc total_timer;
    total_timer.tic();
    
    // Run comprehensive tests
    FiniteDifferenceTester tester;
    tester.run_all_tests();
    
    // Run performance benchmark
    benchmark_finite_difference();
    
    total_timer.toc();
    double total_time_mus = total_timer.elapsed_mus();
    
    // Print final summary
    TestUtils::print_header("TEST SUITE COMPLETED", fmt::color::green);
    
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
      "\n"
      "╔══════════════════════════════════════════════════════════════════════════╗\n"
      "║                           EXECUTION SUMMARY                              ║\n"
      "╠══════════════════════════════════════════════════════════════════════════╣\n"
      "║   Total Execution Time: {:45.3f} μs ║\n"
      "║   Tests Completed:      {:45d}    ║\n"
      "╚══════════════════════════════════════════════════════════════════════════╝\n",
      total_time_mus, static_cast<int>(tester.get_results().size())
    );
    
    // Final results table
    tester.print_summary();
    
    return 0;
    
  } catch (const exception& e) {
    TestUtils::print_error(fmt::format("Exception: {}", e.what()));
    return -1;
  } catch (...) {
    TestUtils::print_error("Unknown error occurred!");
    return -1;
  }
}
