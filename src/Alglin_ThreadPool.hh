#pragma once

#ifndef THREADPOOL_HH
#define THREADPOOL_HH

#include "Alglin_Config.hh"

#include <algorithm>
#include <utility>
#include <vector>
#include <type_traits>

#include <thread>
#include <condition_variable>
#include <mutex>
#include <functional>

namespace alglin {

  class SpinLock {
    #ifdef ALGLIN_OS_WINDOWS
    std::atomic_flag locked;
    #else
    std::atomic_flag locked = ATOMIC_FLAG_INIT;
    #endif
  public:

    #ifdef ALGLIN_OS_WINDOWS
    SpinLock() { locked.clear(); }
    #else
    SpinLock() {}
    #endif

    void
    lock() {
      while (locked.test_and_set(std::memory_order_acquire)) { ; }
    }

    void
    unlock() {
      locked.clear(std::memory_order_release);
    }

  };

  class SpinLock_barrier {
  private:
    std::atomic<unsigned> m_count;
    std::atomic<unsigned> m_generation;
    unsigned int m_count_reset_value;
  public:
    SpinLock_barrier(const SpinLock_barrier&) = delete;
    SpinLock_barrier& operator=(const SpinLock_barrier&) = delete;

    explicit
    SpinLock_barrier()
    : m_generation(0)
    {}

    void
    setup( unsigned count ) {
      m_count_reset_value = m_count = count ;
    }

    void
    count_down() {
      unsigned gen = m_generation.load();
      if ( --m_count == 0 ) {
        if ( m_generation.compare_exchange_weak(gen, gen + 1) )
          m_count = m_count_reset_value;
        return;
      }
    }

    void
    wait() {
      unsigned gen = m_generation.load();
      while ((gen == m_generation) && (m_count != 0))
        std::this_thread::yield();
    }

    void
    count_down_and_wait() {
      unsigned gen = m_generation.load();
      if ( --m_count == 0 ) {
        if ( m_generation.compare_exchange_weak(gen, gen + 1) )
          m_count = m_count_reset_value;
        return;
      }
      while ((gen == m_generation) && (m_count != 0))
        std::this_thread::yield();
    }
  };

  class Barrier {
    int to_be_done, usedThread;
    std::mutex              mtx;
    std::condition_variable cond;
  public:
    Barrier() : to_be_done(0) {}

    void
    setup( int nthreads )
    { usedThread = to_be_done = nthreads ; }

    void
    count_down() {
      std::unique_lock<std::mutex> lck(mtx);
      if ( --to_be_done <= 0 ) cond.notify_all() ; // wake up all tread
    }

    void
    wait() {
      std::unique_lock<std::mutex> lck(mtx);
      cond.wait(lck);
    }

    void
    count_down_and_wait() {
      std::unique_lock<std::mutex> lck(mtx);
      if ( --to_be_done <= 0 ) {
        cond.notify_all() ; // wake up all tread
        to_be_done = usedThread ;
      } else {
        cond.wait(lck);
      }
    }
  };

  class SimpleSemaphore {
  private:
    bool                    m_go;
    std::mutex              m_mutex;
    std::condition_variable m_cv;
  public:
    SimpleSemaphore() noexcept : m_go(true) {}

    void
    green() noexcept {
      { std::unique_lock<std::mutex> lock(m_mutex); m_go = true; }
      m_cv.notify_one();
    }

    void
    red() noexcept {
      { std::unique_lock<std::mutex> lock(m_mutex); m_go = false; }
      m_cv.notify_one();
    }

    void
    wait() noexcept {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_cv.wait(lock, [&]()->bool { return m_go; });
    }

  };

  /*\
   |  __        __         _
   |  \ \      / /__  _ __| | _____ _ __
   |   \ \ /\ / / _ \| '__| |/ / _ \ '__|
   |    \ V  V / (_) | |  |   <  __/ |
   |     \_/\_/ \___/|_|  |_|\_\___|_|
  \*/

  class Worker {
    bool                  active;
    SimpleSemaphore       is_running, job_done;
    std::thread           running_thread;
    std::function<void()> job;

    //disable copy
    Worker( Worker const & ) = delete;
    Worker& operator = ( Worker const & ) = delete;

    void
    loop() {
      while ( active ) {
        is_running.wait();
        if ( active ) job();
        is_running.red();
        job_done.green();
      }
    }

  public:

    Worker() : active(false) { start(); }
    ~Worker() { stop(); }

    Worker( Worker && rhs ) {
      active         = rhs.active;
      job            = rhs.job;
      running_thread = std::move(rhs.running_thread);
    }

    void
    start() {
      if ( !active ) {
        active = true;
        is_running.red();
        job_done.green();
        running_thread = std::thread( [&] () -> void { loop(); } );
      }
    }

    void
    stop() {
      if ( active ) {
        active = false;        // deactivate computation
        is_running.green();    // for exiting from the loop
        running_thread.join(); // wait thread for exiting
      }
    }

    void wait() { job_done.wait(); }

    template < class Func, class... Args >
    void
    run( Func && func, Args && ... args ) {
      //launch( std::bind(std::forward<Func>(func), std::forward<Args>(args)...) );
      job_done.wait(); // se gia occupato in task aspetta
      job = std::bind(std::forward<Func>(func), std::forward<Args>(args)...);
      job_done.red();
      is_running.green(); // activate computation
    }

  };

  /*\
   |   _____ _                        _ ____             _
   |  |_   _| |__  _ __ ___  __ _  __| |  _ \ ___   ___ | |
   |    | | | '_ \| '__/ _ \/ _` |/ _` | |_) / _ \ / _ \| |
   |    | | | | | | | |  __/ (_| | (_| |  __/ (_) | (_) | |
   |    |_| |_| |_|_|  \___|\__,_|\__,_|_|   \___/ \___/|_|
  \*/

  class ThreadPool {
    // need to keep track of threads so we can join them
    std::vector<Worker> workers;

    //disable copy
    ThreadPool() = delete;
    ThreadPool( ThreadPool const & ) = delete;
    ThreadPool& operator = ( ThreadPool const & ) = delete;

  public:

    ThreadPool(
      integer nthread = std::max(
        integer(1),
        integer(std::thread::hardware_concurrency()-1)
      )
    ) {
      workers.resize( size_t( nthread ) );
    }

    //! Submit a job to be run by the thread pool.
    template <typename Func, typename... Args>
    void
    run( integer nt, Func && func, Args && ... args ) {
      workers[ size_t(nt) ].run( func, args...);
    }

    void wait_all()  { for ( auto && w : workers ) w.wait(); }
    void start_all() { for ( auto && w : workers ) w.start(); }
    void stop_all()  { for ( auto && w : workers ) w.stop(); }

    integer size() const { return integer(workers.size()); }

    void
    resize( integer numThreads ) {
      wait_all();
      stop_all();
      workers.resize( size_t(numThreads) );
      start_all();
    }

  };

}

#endif
