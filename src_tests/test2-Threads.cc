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
#include "Utils_TicToc.hh"
#include <random>
#include <vector>
#include <numeric>
#include <algorithm>
#include <atomic>

using namespace std;
using namespace alglin;

// ===========================================================================
// Utility functions for formatted output
// ===========================================================================

namespace ThreadTestUtils {

static void print_header(const std::string& title, fmt::color color = fmt::color::cyan) {
  const int width = 78;
  fmt::print("\n");
  fmt::print(fg(color) | fmt::emphasis::bold, "╔{:═^{}}╗\n", "", width);
  fmt::print(fg(color) | fmt::emphasis::bold, "║{:^{}}║\n", title, width);
  fmt::print(fg(color) | fmt::emphasis::bold, "╚{:═^{}}╝\n", "", width);
}

static void print_subheader(const std::string& title, fmt::color color = fmt::color::yellow) {
  fmt::print(fg(color) | fmt::emphasis::bold, "\n┌{0:─^{2}}┐\n", "", title, 76);
  fmt::print(fg(color) | fmt::emphasis::bold, "│ {:<74} │\n", title);
  fmt::print(fg(color) | fmt::emphasis::bold, "└{:─^{}}┘\n", "", 76);
}

static void print_success(const std::string& message) {
  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "✅  {}\n", message);
}

static void print_warning(const std::string& message) {
  fmt::print(fg(fmt::color::orange) | fmt::emphasis::bold, "⚠️   {}\n", message);
}

static void print_error(const std::string& message) {
  fmt::print(fg(fmt::color::red) | fmt::emphasis::bold, "❌  {}\n", message);
}

static void print_info(const std::string& message) {
  fmt::print(fg(fmt::color::blue), "🔹  {}\n", message);
}

} // namespace ThreadTestUtils

// ===========================================================================
// Thread Statistics Structure
// ===========================================================================

struct ThreadStatistics {
  int thread_id;
  double creation_time_ms;
  double start_time_ms;
  double phase1_time_ms;
  double barrier1_wait_ms;
  double phase2_time_ms;
  double barrier2_wait_ms;
  double completion_time_ms;
  double total_time_ms;
  
  // For concurrent access
  std::atomic<bool> phase1_completed{false};
  std::atomic<bool> phase2_completed{false};
};

// ===========================================================================
// Test Functions with different workloads
// ===========================================================================

class ThreadWorkload {
private:
  static std::mutex s_output_mutex;
  static std::atomic<int> s_active_threads;
  
public:
  static void simple_workload(int thread_id, ThreadStatistics& stats, 
                              Utils::Barrier& barrier, std::mt19937& gen) {
    Utils::TicToc timer;
    
    // Record creation time
    stats.creation_time_ms = timer.elapsed_ms();
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::cyan), "🧵 Thread {:3d} started\n", thread_id);
    }
    
    // Phase 1: Random computation
    timer.tic();
    int work_ms = static_cast<int>(gen() % 500 + 100); // 100-600 ms
    Utils::sleep_for_milliseconds(work_ms);
    timer.toc();
    stats.phase1_time_ms = timer.elapsed_ms();
    stats.phase1_completed = true;
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::green), "📊 Thread {:3d} completed Phase 1 ({:6.1f} ms)\n", 
                 thread_id, stats.phase1_time_ms);
    }
    
    // Wait at barrier 1
    timer.tic();
    barrier.count_down_and_wait();
    timer.toc();
    stats.barrier1_wait_ms = timer.elapsed_ms();
    
    // Phase 2: Different random computation
    timer.tic();
    work_ms = static_cast<int>(gen() % 800 + 200); // 200-1000 ms
    Utils::sleep_for_milliseconds(work_ms);
    timer.toc();
    stats.phase2_time_ms = timer.elapsed_ms();
    stats.phase2_completed = true;
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::yellow), "📈 Thread {:3d} completed Phase 2 ({:6.1f} ms)\n", 
                 thread_id, stats.phase2_time_ms);
    }
    
    // Wait at barrier 2
    timer.tic();
    barrier.count_down_and_wait();
    timer.toc();
    stats.barrier2_wait_ms = timer.elapsed_ms();
    
    stats.completion_time_ms = stats.creation_time_ms + 
                              stats.phase1_time_ms + 
                              stats.barrier1_wait_ms +
                              stats.phase2_time_ms + 
                              stats.barrier2_wait_ms;
    stats.total_time_ms = stats.completion_time_ms;
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::magenta), "🏁 Thread {:3d} completed all work ({:6.1f} ms total)\n", 
                 thread_id, stats.total_time_ms);
    }
    
    s_active_threads--;
  }
  
  static void matrix_computation(int thread_id, ThreadStatistics& stats,
                                 Utils::Barrier& barrier, std::mt19937& gen) {
    Utils::TicToc timer;
    
    // Record start time
    stats.start_time_ms = timer.elapsed_ms();
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::blue), "🔢 Thread {:3d} starting matrix computation\n", thread_id);
    }
    
    // Phase 1: Matrix multiplication simulation
    timer.tic();
    int size = 100 + (gen() % 100); // 100-200
    vector<vector<real_type>> A(size, vector<real_type>(size));
    vector<vector<real_type>> B(size, vector<real_type>(size));
    vector<vector<real_type>> C(size, vector<real_type>(size));
    
    // Fill with random values
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        A[i][j] = static_cast<real_type>(gen()) / gen.max();
        B[i][j] = static_cast<real_type>(gen()) / gen.max();
      }
    }
    
    // Simple matrix multiplication
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        real_type sum = 0;
        for (int k = 0; k < size; ++k) {
          sum += A[i][k] * B[k][j];
        }
        C[i][j] = sum;
      }
    }
    timer.toc();
    stats.phase1_time_ms = timer.elapsed_ms();
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::green), "🧮 Thread {:3d} completed matrix {}x{} ({:6.1f} ms)\n", 
                 thread_id, size, size, stats.phase1_time_ms);
    }
    
    // Barrier
    timer.tic();
    barrier.count_down_and_wait();
    timer.toc();
    stats.barrier1_wait_ms = timer.elapsed_ms();
    
    // Phase 2: Vector operations
    timer.tic();
    int vec_size = 100000 + (gen() % 900000); // 100k-1M
    vector<real_type> v1(vec_size), v2(vec_size);
    
    for (int i = 0; i < vec_size; ++i) {
      v1[i] = static_cast<real_type>(gen()) / gen.max();
      v2[i] = static_cast<real_type>(gen()) / gen.max();
    }
    
    real_type dot_product = 0;
    for (int i = 0; i < vec_size; ++i) {
      dot_product += v1[i] * v2[i];
    }
    
    // Some additional computation
    real_type norm1 = 0, norm2 = 0;
    for (int i = 0; i < vec_size; ++i) {
      norm1 += v1[i] * v1[i];
      norm2 += v2[i] * v2[i];
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);
    timer.toc();
    stats.phase2_time_ms = timer.elapsed_ms();
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::yellow), "📐 Thread {:3d} completed vector ops (size: {}, dot: {:.6f}) ({:6.1f} ms)\n", 
                 thread_id, vec_size, dot_product, stats.phase2_time_ms);
    }
    
    // Final barrier
    timer.tic();
    barrier.count_down_and_wait();
    timer.toc();
    stats.barrier2_wait_ms = timer.elapsed_ms();
    
    stats.total_time_ms = stats.phase1_time_ms + stats.phase2_time_ms +
                         stats.barrier1_wait_ms + stats.barrier2_wait_ms;
    
    {
      std::lock_guard<std::mutex> lock(s_output_mutex);
      fmt::print(fg(fmt::color::magenta), "🎯 Thread {:3d} finished computation ({:6.1f} ms total)\n", 
                 thread_id, stats.total_time_ms);
    }
    
    s_active_threads--;
  }
};

// Static member initialization
std::mutex ThreadWorkload::s_output_mutex;
std::atomic<int> ThreadWorkload::s_active_threads{0};

// ===========================================================================
// Test Scenarios
// ===========================================================================

void test_simple_threads(int num_threads) {
  ThreadTestUtils::print_header("TEST 1: Simple Threads with Barriers", fmt::color::green);
  
  vector<ThreadStatistics> stats(num_threads);
  vector<thread> threads;
  threads.reserve(num_threads);
  
  Utils::Barrier barrier;
  barrier.setup(num_threads);
  std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
  
  ThreadTestUtils::print_info(fmt::format("Starting {} threads with barrier synchronization", num_threads));
  
  Utils::TicToc total_timer;
  total_timer.tic();
  
  // Launch threads
  for (int i = 0; i < num_threads; ++i) {
    stats[i].thread_id = i;
    threads.emplace_back([i, &stats, &barrier, &gen]() {
      ThreadWorkload::simple_workload(i, stats[i], barrier, gen);
    });
  }
  
  // Wait for all threads to complete
  for (auto& t : threads) {
    t.join();
  }
  
  total_timer.toc();
  double total_time_ms = total_timer.elapsed_ms();
  
  // Print statistics table
  ThreadTestUtils::print_subheader("Thread Performance Statistics");
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "┌─────┬──────────┬──────────┬──────────┬────────────┬──────────┬────────────┬──────────┐\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│ ID  │ Phase 1  │ Barrier1 │ Phase 2  │ Barrier 2  │ Total    │ Efficiency │ Status   │\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│     │ (ms)     │ Wait(ms) │ (ms)     │ Wait(ms)   │ (ms)     │ (%)        │          │\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "├─────┼──────────┼──────────┼──────────┼────────────┼──────────┼────────────┼──────────┤\n");
  
  double avg_phase1 = 0, avg_phase2 = 0, avg_total = 0;
  double max_phase1 = 0, max_phase2 = 0, max_total = 0;
  
  for (int i = 0; i < num_threads; ++i) {
    const auto& s = stats[i];
    double compute_time = s.phase1_time_ms + s.phase2_time_ms;
    double efficiency = (compute_time / s.total_time_ms) * 100.0;
    
    avg_phase1 += s.phase1_time_ms;
    avg_phase2 += s.phase2_time_ms;
    avg_total += s.total_time_ms;
    
    max_phase1 = max(max_phase1, s.phase1_time_ms);
    max_phase2 = max(max_phase2, s.phase2_time_ms);
    max_total = max(max_total, s.total_time_ms);
    
    fmt::color status_color = efficiency > 80 ? fmt::color::green :
                             efficiency > 60 ? fmt::color::yellow : fmt::color::red;
    
    fmt::print("│ {:3d} │ {:8.1f} │ {:8.1f} │ {:8.1f} │ {:10.1f} │ {:8.1f} │ {:10.1f} │ ",
               i, s.phase1_time_ms, s.barrier1_wait_ms, s.phase2_time_ms,
               s.barrier2_wait_ms, s.total_time_ms, efficiency);
    fmt::print(fg(status_color) | fmt::emphasis::bold, "{:^8} │\n",
               efficiency > 80 ? "Good" : efficiency > 60 ? "Fair" : "Poor");
  }
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "├─────┼──────────┼──────────┼──────────┼────────────┼──────────┼────────────┼──────────┤\n");
  
  avg_phase1 /= num_threads;
  avg_phase2 /= num_threads;
  avg_total /= num_threads;
  
  fmt::print("│ Avg │ {:8.1f} │          │ {:8.1f} │            │ {:8.1f} │            │          │\n",
             avg_phase1, avg_phase2, avg_total);
  fmt::print("│ Max │ {:8.1f} │          │ {:8.1f} │            │ {:8.1f} │            │          │\n",
             max_phase1, max_phase2, max_total);
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "└─────┴──────────┴──────────┴──────────┴────────────┴──────────┴────────────┴──────────┘\n");
  
  // Summary
  ThreadTestUtils::print_subheader("Test Summary");
  fmt::print("Total test time: {:8.1f} ms\n", total_time_ms);
  fmt::print("Number of threads: {}\n", num_threads);
  fmt::print("Average thread time: {:8.1f} ms\n", avg_total);
  
  double speedup_ideal = max_phase1 + max_phase2;
  double speedup_actual = total_time_ms;
  double parallel_efficiency = (speedup_ideal / speedup_actual) * 100.0;
  
  fmt::print("Theoretical best time: {:8.1f} ms\n", speedup_ideal);
  fmt::print("Parallel efficiency: {:6.1f}%\n", parallel_efficiency);
  
  if (parallel_efficiency > 70) {
    ThreadTestUtils::print_success("Excellent parallel efficiency!");
  } else if (parallel_efficiency > 50) {
    ThreadTestUtils::print_warning("Moderate parallel efficiency");
  } else {
    ThreadTestUtils::print_warning("Low parallel efficiency - consider thread pool optimization");
  }
}

void test_computation_threads(int num_threads) {
  ThreadTestUtils::print_header("TEST 2: Computational Workload", fmt::color::blue);
  
  vector<ThreadStatistics> stats(num_threads);
  vector<thread> threads;
  threads.reserve(num_threads);
  
  Utils::Barrier barrier;
  barrier.setup(num_threads);
  std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
  
  ThreadTestUtils::print_info(fmt::format("Starting {} threads with computational workload", num_threads));
  
  Utils::TicToc total_timer;
  total_timer.tic();
  
  // Launch threads with matrix computation
  for (int i = 0; i < num_threads; ++i) {
    stats[i].thread_id = i;
    threads.emplace_back([i, &stats, &barrier, &gen]() {
      ThreadWorkload::matrix_computation(i, stats[i], barrier, gen);
    });
  }
  
  // Wait for all threads to complete
  for (auto& t : threads) {
    t.join();
  }
  
  total_timer.toc();
  double total_time_ms = total_timer.elapsed_ms();
  
  // Print summary table
  ThreadTestUtils::print_subheader("Computation Threads Summary");
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "┌─────┬────────────────┬────────────────┬────────────────┬────────────────┐\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│ ID  │ Matrix Time    │ Barrier Wait   │ Vector Time    │ Total Time     │\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│     │ (ms)           │ (ms)           │ (ms)           │ (ms)           │\n");
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "├─────┼────────────────┼────────────────┼────────────────┼────────────────┤\n");
  
  double total_compute_time = 0;
  double total_wait_time = 0;
  
  for (int i = 0; i < num_threads; ++i) {
    const auto& s = stats[i];
    double compute_time = s.phase1_time_ms + s.phase2_time_ms;
    double wait_time = s.barrier1_wait_ms + s.barrier2_wait_ms;
    
    total_compute_time += compute_time;
    total_wait_time += wait_time;
    
    fmt::print("│ {:3d} │ {:14.1f} │ {:14.1f} │ {:14.1f} │ {:14.1f} │\n",
               i, s.phase1_time_ms, s.barrier1_wait_ms, s.phase2_time_ms, s.total_time_ms);
  }
  
  fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "└─────┴────────────────┴────────────────┴────────────────┴────────────────┘\n");
  
  // Statistics
  ThreadTestUtils::print_subheader("Performance Analysis");
  
  double avg_compute_time = total_compute_time / num_threads;
  double avg_wait_time = total_wait_time / num_threads;
  double cpu_utilization = (total_compute_time / (num_threads * total_time_ms)) * 100.0;
  
  fmt::print("Total test duration: {:10.1f} ms\n", total_time_ms);
  fmt::print("Average compute time per thread: {:8.1f} ms\n", avg_compute_time);
  fmt::print("Average wait time per thread: {:8.1f} ms\n", avg_wait_time);
  fmt::print("Total CPU time (sum of all threads): {:8.1f} ms\n", total_compute_time);
  fmt::print("CPU utilization: {:6.1f}%\n", cpu_utilization);
  
  if (cpu_utilization > 80) {
    ThreadTestUtils::print_success("High CPU utilization - good thread distribution");
  } else if (cpu_utilization > 50) {
    ThreadTestUtils::print_warning("Moderate CPU utilization");
  } else {
    ThreadTestUtils::print_warning("Low CPU utilization - threads may be waiting too much");
  }
}

void test_thread_scalability() {
  ThreadTestUtils::print_header("TEST 3: Thread Scalability Analysis", fmt::color::magenta);
  
  vector<int> thread_counts = {1, 2, 4, 8, 16};
  vector<double> execution_times;
  vector<double> efficiencies;
  
  std::string table;
  
  table = fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "┌──────────┬────────────────┬────────────────┬────────────────┬────────────────┐\n");
  table += fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│ Threads  │ Serial Time    │ Parallel Time  │ Speedup        │ Efficiency     │\n");
  table += fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "│          │ (ms)           │ (ms)           │ (ideal = N)    │ (%)            │\n");
  table += fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "├──────────┼────────────────┼────────────────┼────────────────┼────────────────┤\n");
  
  double serial_time = 0;
  
  for (int num_threads : thread_counts) {
    Utils::TicToc timer;
    timer.tic();
    
    vector<thread> threads;
    Utils::Barrier barrier;
    barrier.setup(num_threads);
    std::mt19937 gen(static_cast<unsigned>(time(nullptr) + num_threads));
    
    // Launch threads
    for (int i = 0; i < num_threads; ++i) {
      threads.emplace_back([i, &barrier, &gen]() {
        ThreadStatistics stats;
        ThreadWorkload::simple_workload(i, stats, barrier, gen);
      });
    }
    
    // Wait for completion
    for (auto& t : threads) {
      t.join();
    }
    
    timer.toc();
    double time_ms = timer.elapsed_ms();
    execution_times.push_back(time_ms);
    
    if (num_threads == 1) {
      serial_time = time_ms;
    }
    
    double speedup = serial_time / time_ms;
    double efficiency = (speedup / num_threads) * 100.0;
    efficiencies.push_back(efficiency);
    
    fmt::color speedup_color = speedup > num_threads * 0.7 ? fmt::color::green :
                              speedup > num_threads * 0.5 ? fmt::color::yellow : fmt::color::red;
    
    table += fmt::format("│ {:8d} │ {:14.1f} │ {:14.1f} │ ", num_threads, serial_time, time_ms);
    table += fmt::format(fg(speedup_color), "{:14.2f} │ {:14.1f} │\n", speedup, efficiency);
  }
  
  table += fmt::format(fg(fmt::color::cyan) | fmt::emphasis::bold,
    "└──────────┴────────────────┴────────────────┴────────────────┴────────────────┘\n");
    
  std::cout << table;
  
  // Scalability analysis
  ThreadTestUtils::print_subheader("Scalability Insights");
  
  // Find the optimal thread count
  auto max_efficiency = max_element(efficiencies.begin(), efficiencies.end());
  int optimal_threads = thread_counts[distance(efficiencies.begin(), max_efficiency)];
  
  fmt::print("Optimal thread count: {} (efficiency: {:.1f}%)\n", 
             optimal_threads, *max_efficiency);
  
  if (*max_efficiency > 80) {
    ThreadTestUtils::print_success("Excellent scalability up to optimal thread count");
  } else if (*max_efficiency > 60) {
    ThreadTestUtils::print_warning("Moderate scalability - diminishing returns after optimal count");
  } else {
    ThreadTestUtils::print_warning("Poor scalability - consider optimizing thread synchronization");
  }
  
  // Recommendation
  fmt::print("\n");
  ThreadTestUtils::print_info("Recommendation:");
  if (optimal_threads <= 4) {
    fmt::print("🔹 Consider using {} threads for best performance/efficiency tradeoff\n", optimal_threads);
  } else {
    fmt::print("🔹 System shows good scalability up to {} threads\n", optimal_threads);
  }
}

// ===========================================================================
// Main Test Function
// ===========================================================================

int main() {
  try {
    ThreadTestUtils::print_header("THREAD CONCURRENCY TEST SUITE", fmt::color::green);
    
    fmt::print(fg(fmt::color::yellow), "🚀 Starting comprehensive thread concurrency tests\n\n");
    
    Utils::TicToc total_timer;
    total_timer.tic();
    
    // Test 1: Simple threads with barriers
    test_simple_threads(8);
    
    // Test 2: Computational workload
    test_computation_threads(4);
    
    // Test 3: Scalability analysis
    test_thread_scalability();
    
    total_timer.toc();
    double total_time_ms = total_timer.elapsed_ms();
    
    // Final summary
    ThreadTestUtils::print_header("TEST SUITE COMPLETED", fmt::color::green);
    
    fmt::print(fg(fmt::color::cyan), "\n");
    fmt::print(fg(fmt::color::cyan) | fmt::emphasis::bold,
      "╔══════════════════════════════════════════════════════════════════════════╗\n"
      "║                           EXECUTION SUMMARY                              ║\n"
      "╠══════════════════════════════════════════════════════════════════════════╣\n"
      "║   Total Execution Time: {:45.3f} ms ║\n"
      "║   Tests Completed:      3                                                ║\n"
      "╚══════════════════════════════════════════════════════════════════════════╝\n",
      total_time_ms
    );
    
    fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, 
               "\n✅ All thread concurrency tests completed successfully!\n");
    
    return 0;
    
  } catch (const exception& e) {
    ThreadTestUtils::print_error(fmt::format("Exception: {}", e.what()));
    return -1;
  } catch (...) {
    ThreadTestUtils::print_error("Unknown error occurred!");
    return -1;
  }
}
