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
#include <vector>
#include <random>
#include <fstream>
#include "Alglin_threads.hh"
#include "TicToc.hh"

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

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
static std::mt19937 generator(seed1);

static TicToc tictoc ;
static std::mutex mtx ;

static alglin::Barrier bar ;

static
void
test( int nth ) {

  int ms = generator() % 2000 ;

  mtx.lock() ;
  cout << "Thread N." << nth << " ms = " << ms << "\n" ;
  mtx.unlock() ;

  tictoc.sleep_for_milliseconds( ms ) ;

  mtx.lock() ;
  cout << "Thread N." << nth << " done_and_wait\n" ;
  mtx.unlock() ;

  bar.count_down_and_wait() ;
  
  ms = generator() % 2000 ;

  mtx.lock() ;
  cout << "Thread N." << nth << " second part ms = " << ms << "\n" ;
  mtx.unlock() ;

  tictoc.sleep_for_milliseconds( ms ) ;

  mtx.lock() ;
  cout << "Thread N." << nth << " done second part\n" ;
  mtx.unlock() ;
  bar.count_down_and_wait() ;
}

int
main() {
  std::thread threads[100] ;
  int usedThread = 10 ;

  bar.setup(usedThread) ;
  for ( int nt = 0 ; nt < usedThread ; ++nt )
    threads[nt] = std::thread( &test, nt ) ;

  for ( int nt = 0 ; nt < usedThread ; ++nt ) {
    threads[nt].join() ;
    mtx.lock() ;
    cout << "Thread N." << nt << " joined\n" ;
    mtx.unlock() ;
  }

  cout << "All done!\n" ;

  return 0 ;
}
