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

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wdocumentation"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wc99-extensions"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wdeprecated"
#pragma GCC diagnostic ignored "-Wused-but-marked-unused"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wswitch-enum"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include "Alglin_threads.hh"
#include "TicToc.hh"

using namespace std ;
typedef double valueType ;

static unsigned seed1 = 2 ;
static std::mt19937 generator(seed1);

TicToc tictoc ;
std::mutex mtx ;

static
valueType
rand( valueType xmin, valueType xmax ) {
  valueType random = valueType(generator())/generator.max();
  return xmin + (xmax-xmin)*random ;
}

alglin::Barrier bar ;

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
