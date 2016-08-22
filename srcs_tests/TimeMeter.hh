/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010                                                      |
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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 |      version: 0.2 07-03-2011                                             |
 |                                                                          |
 \*--------------------------------------------------------------------------*/

#ifndef TIME_METER_HH
#define TIME_METER_HH

#include <chrono>

class TimeMeter {

  std::chrono::system_clock::time_point start_time ;
  std::chrono::system_clock::time_point stop_time ;

  double elapsedTotal, elapsedPartial ;

  TimeMeter( TimeMeter const & ) ;
  TimeMeter const & operator = ( TimeMeter const & ) const ;

public:

  TimeMeter() : elapsedTotal(0) { start() ; }
  ~TimeMeter() {}

  void reset() { elapsedTotal = 0 ; }

  void
  start() {
    start_time = std::chrono::system_clock::now();
  }

  void
  stop() {
    stop_time = std::chrono::system_clock::now();
    elapsedPartial = double(1e-6*std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time).count());
    elapsedTotal  += elapsedPartial ;
  }

  double totalElapsedSeconds()      const { return elapsedTotal ; }
  double totalElapsedMilliseconds() const { return elapsedTotal*1000 ; }

  double partialElapsedSeconds()      const { return elapsedPartial ; }
  double partialElapsedMilliseconds() const { return elapsedPartial*1000 ; }

} ;

#endif
