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
/// file: Alglin_threads.hh
///

#ifndef ALGLIN_THREADS_HH
#define ALGLIN_THREADS_HH

/*
*/

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#if defined(__GCC__) || defined(__GNUC__) 
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++98-compat"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++98-compat"
#endif

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
      while (locked.test_and_set(std::memory_order_acquire)) {; }
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
      m_count_reset_value = m_count = count;
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class Barrier {
    int to_be_done, usedThread;
    std::mutex              mtx;
    std::condition_variable cond;
  public:
    Barrier() : to_be_done(0) {}

    void
    setup( int nthreads )
    { usedThread = to_be_done = nthreads; }

    void
    count_down() {
      std::unique_lock<std::mutex> lck(mtx);
      if ( --to_be_done <= 0 ) cond.notify_all(); // wake up all tread
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
        cond.notify_all(); // wake up all tread
        to_be_done = usedThread;
      } else {
        cond.wait(lck);
      }
    }
  };
} // end namespace alglin

#if defined(__GCC__) || defined(__GNUC__) 
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif

///
/// eof: Alglin_threads.hh
///
