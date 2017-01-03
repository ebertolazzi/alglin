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

  bar.done_and_wait() ;
  
  ms = generator() % 2000 ;

  mtx.lock() ;
  cout << "Thread N." << nth << " second part ms = " << ms << "\n" ;
  mtx.unlock() ;

  tictoc.sleep_for_milliseconds( ms ) ;

  mtx.lock() ;
  cout << "Thread N." << nth << " done second part\n" ;
  mtx.unlock() ;
  bar.done_and_wait() ;
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
