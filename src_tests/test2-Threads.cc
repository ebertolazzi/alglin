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
#include <random>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

static Utils::Console msg( &std::cout,4);

using namespace std;
typedef double real_type;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static Utils::TicToc tictoc;
static std::mutex mtx;

static Utils::Barrier bar;

static
void
test( int nth ) {

  int ms{ static_cast<int>( generator() % 2000 ) };

  mtx.lock();
  fmt::print( "Thread N.{} ms = {}\n", nth, ms );
  mtx.unlock();

  Utils::sleep_for_milliseconds( ms );

  mtx.lock();
  fmt::print( "Thread N.{} done_and_wait\n", nth );
  mtx.unlock();

  bar.count_down_and_wait();

  ms = generator() % 2000;

  mtx.lock();
  fmt::print( "Thread N.{} second part ms = {}\n", nth, ms );
  mtx.unlock();

   Utils::sleep_for_milliseconds( ms );

  mtx.lock();
  fmt::print( "Thread N.{} done second part\n", nth );
  mtx.unlock();
  bar.count_down_and_wait();
}

int
main() {

  std::thread threads[100];
  int usedThread{10};

  bar.setup(usedThread);
  for ( int nt{0}; nt < usedThread; ++nt )
    threads[nt] = std::thread( &test, nt );

  for ( int nt{0}; nt < usedThread; ++nt ) {
    threads[nt].join();
    mtx.lock();
    fmt::print( "Thread N.{} joined\n", nt );
    mtx.unlock();
  }

  msg.green( "All done!\n" );
  return 0;
}
