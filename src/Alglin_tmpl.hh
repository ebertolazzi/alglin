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

///
/// file: Alglin_tmpl.hh
///

#ifndef ALGLIN_TMPL_HH
#define ALGLIN_TMPL_HH

/*
*/

#include "Alglin.hh"

namespace alglin {

  /*
  //   _____
  //  |__  /___ _ __ ___
  //    / // _ \ '__/ _ \
  //   / /|  __/ | | (_) |
  //  /____\___|_|  \___/
  */
  template <typename t_Value, int N, int STEP>
  class Zero {
  public:
    static inline void eval( t_Value * A ) {
      A[0] = 0 ; Zero<t_Value,N-1,STEP>::eval( A+STEP ) ;
    }
  } ;

  template <typename t_Value, int STEP>
  class Zero<t_Value,0,STEP> {
  public:
    static inline void eval( t_Value * A ) {}
  } ;
  
  /*
  //    ____
  //   / ___|___  _ __  _   _
  //  | |   / _ \| '_ \| | | |
  //  | |__| (_) | |_) | |_| |
  //   \____\___/| .__/ \__, |
  //             |_|    |___/
  */
  template <typename t_Value, int N, int STEPA, int STEPB>
  class Copy {
  public:
    static inline void eval( t_Value const * A, t_Value * B ) {
      B[0] = A[0] ;
      Copy<t_Value,N-1,STEPA,STEPB>::eval( A+STEPA, B+STEPB ) ;
    }
  } ;

  template <typename t_Value, int STEPA, int STEPB>
  class Copy<t_Value,0,STEPA,STEPB> {
  public:
    static inline void eval( t_Value const * A, t_Value * B ) {}
  } ;
  
  /*
  //   ____
  //  / ___|_      ____ _ _ __
  //  \___ \ \ /\ / / _` | '_ \
  //   ___) \ V  V / (_| | |_) |
  //  |____/ \_/\_/ \__,_| .__/
  //                     |_|
  */
  template <typename t_Value, int N, int STEPA, int STEPB>
  class Swap {
  public:
    static inline void eval( t_Value * A, t_Value * B ) {
      std::swap(B[0],A[0]) ; Swap<t_Value,N-1,STEPA,STEPB>::eval( A+STEPA, B+STEPB ) ;
    }
  } ;

  template <typename t_Value, int STEPA, int STEPB>
  class Swap<t_Value,0,STEPA,STEPB> {
  public:
    static inline void eval( t_Value * A, t_Value * B ) {}
  } ;
  
  /*
  //   ____        _
  //  |  _ \  ___ | |_
  //  | | | |/ _ \| __|
  //  | |_| | (_) | |_
  //  |____/ \___/ \__|
  */
  template <typename t_Value, int N, int STEPA, int STEPB>
  class Dot {
  public:
    static inline t_Value eval( t_Value const * A, t_Value const * B ) {
      return B[0]*A[0]+Dot<t_Value,N-1,STEPA,STEPB>::eval( A+STEPA, B+STEPB ) ;
    }
  } ;

  template <typename t_Value, int STEPA, int STEPB>
  class Dot<t_Value,0,STEPA,STEPB> {
  public:
    static inline t_Value eval( t_Value const * A, t_Value const * B )
    { return t_Value(0) ; }
  } ;
  
  /*
  //   ____            _
  //  / ___|  ___ __ _| |
  //  \___ \ / __/ _` | |
  //   ___) | (_| (_| | |
  //  |____/ \___\__,_|_|
  */
  template <typename t_Value, int N, int STEP>
  class Scal {
  public:
    static inline void eval( t_Value s, t_Value * x ) {
      x[0] *= s ;
      Scal<t_Value,N-1,STEP>::eval(s,x+STEP) ;
    }
  } ;

  template <typename t_Value, int STEP>
  class Scal<t_Value,0,STEP> {
  public:
    static inline void eval( t_Value s, t_Value * x ) {}
  } ;
  
  /*
  //      _
  //     / \   _ __ ___   __ ___  __
  //    / _ \ | '_ ` _ \ / _` \ \/ /
  //   / ___ \| | | | | | (_| |>  <
  //  /_/   \_\_| |_| |_|\__,_/_/\_\
  */
  template <typename t_Value, int N, int STEP>
  class Amax {
  public:
    static inline void eval( t_Value const * x, t_Value & am, integer & ipos ) {
      integer ipos1 ;
      t_Value am1 ;
      am = std::abs(x[0]) ;
      Amax<t_Value,N-1,STEP>::eval(x+STEP,am1,ipos1) ;
      if ( am < am1 ) { am = am1 ; ipos = ipos1+1 ; }
      else ipos = 0 ;
    }
  } ;

  template <typename t_Value, int STEP>
  class Amax<t_Value,1,STEP> {
  public:
    static inline void eval( t_Value const * x, t_Value & am, integer & ipos )
    { am = std::abs(x[0]) ; ipos = 0 ; }
  } ;

  /*
  //  __  __        __   __
  //  \ \/ /___  _ _\ \ / /
  //   \  // _ \| '_ \ V /
  //   /  \ (_) | |_) | |
  //  /_/\_\___/| .__/|_|
  //            |_|
  //
  //  y  = a*x ;
  //  y += a*x ;
  //  y  = b*y+a*x ;
  //  y -= x ;
  //  y += x ;
  //
  */
  template <typename t_Value, int N, int STEPX, int STEPY>
  class XopY {
  public:
    static inline void xpy( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] += x[0] ;
      XopY<t_Value,N-1,STEPX,STEPY>::xpy(x+STEPX,y+STEPY) ;
    }
    static inline void xmy( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] -= x[0] ;
      XopY<t_Value,N-1,STEPX,STEPY>::xmy(x+STEPX,y+STEPY) ;
    }
    static inline void axy( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] = a*x[0] ;
      XopY<t_Value,N-1,STEPX,STEPY>::axy(a,x+STEPX,y+STEPY) ;
    }
    static inline void axpy( t_Value a, t_Value const * x, t_Value * y ) {
      y[0] += a*x[0] ;
      XopY<t_Value,N-1,STEPX,STEPY>::axpy(a,x+STEPX,y+STEPY) ;
    }
    static inline void axpby( t_Value a, t_Value const * x, t_Value b, t_Value * y ) {
      y[0] = b*y[0] + a*x[0] ;
      XopY<t_Value,N-1,STEPX,STEPY>::axpby(a,x+STEPX,b,y+STEPY) ;
    }
  } ;

  template <typename t_Value, int STEPX, int STEPY>
  class XopY<t_Value,0,STEPX,STEPY> {
  public:
    static inline void xpy( t_Value a, t_Value const * x, t_Value * y ) {}
    static inline void xmy( t_Value a, t_Value const * x, t_Value * y ) {}
    static inline void axy( t_Value a, t_Value const * x, t_Value * y ) {}
    static inline void axpy( t_Value a, t_Value const * x, t_Value * y ) {}
    static inline void axpby( t_Value a, t_Value const * x, t_Value b, t_Value * y ) {}
  } ;
  
  /*
  //   _               _
  //  | |    ___  ___ | |_   _____
  //  | |   / __|/ _ \| \ \ / / _ \
  //  | |___\__ \ (_) | |\ V /  __/
  //  |_____|___/\___/|_| \_/ \___|
  */
  template <typename t_Value, int N, int LDL>
  class LsolveUnit {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) {
      XopY<t_Value,N-1,1,1>::axpy(-x[0],L+1,x+1) ;
      LsolveUnit<t_Value,N-1,LDL>::eval(L+LDL+1,x+1) ;
    }
  } ;

  template <typename t_Value, int LDL>
  class LsolveUnit<t_Value,1,LDL> {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) {}
  } ;

  template <typename t_Value, int N, int LDL>
  class Lsolve {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) {
      x[0] /= L[0] ;
      XopY<t_Value,N-1,1,1>::axpy(-x[0],L+1,x+1) ;
      Lsolve<t_Value,N-1,LDL>::eval(L+LDL+1,x+1) ;
    }
  } ;

  template <typename t_Value, int LDL>
  class Lsolve<t_Value,1,LDL> {
  public:
    static inline void eval( t_Value const * L, t_Value * x ) { x[0] /= L[0] ; }
  } ;

  /*
  //   _   _           _
  //  | | | |___  ___ | |_   _____
  //  | | | / __|/ _ \| \ \ / / _ \
  //  | |_| \__ \ (_) | |\ V /  __/
  //   \___/|___/\___/|_| \_/ \___|
  */
  template <typename t_Value, int N, int LDU>
  class UsolveUnit {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      XopY<t_Value,N-1,1,1>::axpy(-x[N-1],U+(N-1)*LDU,x) ;
      UsolveUnit<t_Value,N-1,LDU>::eval(U,x) ;
    }
  } ;

  template <typename t_Value, int LDU>
  class UsolveUnit<t_Value,1,LDU> {
  public:
    static void eval( t_Value const * U, t_Value * x ) {}
  } ;

  template <typename t_Value, int N, int LDU>
  class Usolve {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      x[N-1] /= U[(N-1)*(LDU+1)] ;
      XopY<t_Value,N-1,1,1>::axpy(-x[N-1],U+(N-1)*LDU,x) ;
      Usolve<t_Value,N-1,LDU>::eval(U,x) ;
    }
  } ;

  template <typename t_Value, int LDU>
  class Usolve<t_Value,1,LDU> {
  public:
    static inline void eval( t_Value const * U, t_Value * x ) {
      x[0] /= U[0] ;
    }
  } ;
  
  /*
  //   __  __
  //  |  \/  |_   __
  //  | |\/| \ \ / /
  //  | |  | |\ V /
  //  |_|  |_| \_/
  */
  template <typename t_Value, int M, int N, int LDM, int INCX, int INCR>
  class Mv {
  public:
    static inline void ass( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      res[0] = Dot<t_Value,N,LDM,INCX>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM,INCX,INCR>::ass( Mat+1, x, res+INCR ) ;
    }
    static inline void addTo( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      res[0] += Dot<t_Value,N,LDM,INCX>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM,INCX,INCR>::addTo( Mat+1, x, res+INCR ) ;
    }
    static inline void subTo( t_Value const * Mat, t_Value const * x, t_Value * res ) {
      res[0] -= Dot<t_Value,N,LDM,INCX>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM,INCX,INCR>::subTo( Mat+1, x, res+INCR ) ;
    }
    static inline void aMxpby( t_Value a, t_Value const * Mat, t_Value const * x,
                               t_Value b, t_Value * res ) {
      res[0] = b*res[0] + a*Dot<t_Value,N,LDM,INCX>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM,INCX,INCR>::subTo( Mat+1, x, res+INCR ) ;
    }
  } ;

  template <typename t_Value, int N, int LDM, int INCX, int INCR>
  class Mv<t_Value,0,N,LDM,INCX,INCR> {
  public:
    static inline void ass( t_Value const * Mat, t_Value const * x, t_Value * res ) { }
    static inline void addTo( t_Value const * Mat, t_Value const * x, t_Value * res ) { }
    static inline void subTo( t_Value const * Mat, t_Value const * x, t_Value * res ) { }
    static inline void aMxpby( t_Value a, t_Value const * Mat, t_Value const * x,
                               t_Value b, t_Value * res ) { }
  } ;
  
  /*
  //   __  __ __  __
  //  |  \/  |  \/  |
  //  | |\/| | |\/| |
  //  | |  | | |  | |
  //  |_|  |_|_|  |_|
  //
  //  C  = A*B
  //  C += A*B
  //  C -= A*B
  //  C  = a*C + b*A*B
  */
  template <typename t_Value, int M, int N, int K, int LDA, int LDB, int LDC>
  class MM {
  public:
    static inline void ass( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0 ; i < N ; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::ass( A, B, C ) ;
    }
    static inline void addTo( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0 ; i < N ; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::addTo( A, B, C ) ;
    }
    static inline void subTo( t_Value const * A, t_Value const * B, t_Value * C ) {
      for ( integer i = 0 ; i < N ; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::subTo( A, B, C ) ;
    }
    static inline void aMxpby( t_Value a, t_Value const * A, t_Value const * B,
                               t_Value b, t_Value * C ) {
      for ( integer i = 0 ; i < N ; ++i, B += LDB, C += LDC )
        Mv<t_Value,M,K,LDA,1,1>::aMxpby( a, A, B, b, C ) ;
    }
  } ;

  /*
  //   ____             _    _
  //  |  _ \ __ _ _ __ | | _/ |
  //  | |_) / _` | '_ \| |/ / |
  //  |  _ < (_| | | | |   <| |
  //  |_| \_\__,_|_| |_|_|\_\_|
  //
  //  M = u*v^T
  //  M += u*v^T
  //  M -= u*v^T
  //  M  = b*M + a*u*v^T
  //
  */
  template <typename t_Value, int M, int N, int LDM, int INCX, int INCY>
  class Rank1 {
  public:
    static inline void ass( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      XopY<t_Value,M,INCX,1>::axy( y[0], x, Mat ) ;
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::ass(x,y+INCY,Mat+LDM) ;
    }
    static inline void addTo( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      XopY<t_Value,M,INCX,1>::axpy( y[0], x, Mat ) ;
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::addTo(x,y+INCY,Mat+LDM) ;
    }
    static inline void subTo( t_Value const * x, t_Value const * y, t_Value * Mat ) {
      XopY<t_Value,M,INCX,1>::axpy( -y[0], x, Mat ) ;
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::subTo(x,y+INCY,Mat+LDM) ;
    }
    static inline void auvpbM( t_Value a, t_Value const * x, t_Value const * y,
                               t_Value b, t_Value * Mat ) {
      XopY<t_Value,M,INCX,1>::axpby( a*y[0], x, b, Mat ) ;
      Rank1<t_Value,M,N-1,LDM,INCX,INCY>::auvpbM(x,y+INCY,Mat+LDM) ;
    }
  } ;

  template <typename t_Value, int M, int LDM, int INCX, int INCY>
  class Rank1<t_Value,M,0,LDM,INCX,INCY> {
  public:
    static inline void ass( t_Value const * x, t_Value const * y, t_Value * Mat ) { }
    static inline void addTo( t_Value const * x, t_Value const * y, t_Value * Mat ) { }
    static inline void subTo( t_Value const * x, t_Value const * y, t_Value * Mat ) { }
    static inline void auvpbM( t_Value a, t_Value const * x, t_Value const * y,
                               t_Value b, t_Value * Mat ) { }
  } ;


} // end namespace alglin

#endif

///
/// eof: Alglin_tmpl.hh
///

