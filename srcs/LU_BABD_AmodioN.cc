/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2003                                                      |
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

#include "LU_BABD_AmodioN.hh"
#include <iostream>
#include <algorithm>

/*\
 |  reduces a 2x3 block matrix / Ad  Au     \ to a 1x2 block matrix ( Ad' Au' )
 |                             \     Bd  Bu /
 |  
 |  by using the following algorithm:
 |
 |    P * / Au \ = / Au' \  = / M \ (U) = / G \ (L*U)
 |        \ Bd /   \ Bd' /    \ L /       \ I /
 |
 |  where G = M * L^(-1).
 |
 |    / I -G \ * P * / Au \ = / I -G \ / G \ (L*U) = /  0  \
 |    \    I /       \ Bd /   \    I / \ I /         \ L*U /
 |
 |
 |    P * / Ad Au    \ = / Ad'    0   Bu'  \
 |        \    Bd Bu /   \ Ad''  L*U  Bu'' /
 |
 |    / I -G \ * P * / Ad Au    \ = / E*    0   F*   \
 |    \    I /       \    Bd Bu /   \ Ad'' L*U  Bu'' /
 |
 |  where  E* = Ad' - G*Ad'',  F* = Bu' - G*Bu''
 |
 | ------------------------------------------------------------
 |
 |  Esempio
 |  +-----+-----+
 |  |  Ad | Au  |
 |  +-----+-----+-----+
 |        | Ad  | Au  |
 |        +-----+-----+-----+
 |              | Ad  | Au  |
 |              +-----+-----+-----+
 |                    | Ad  | Au  |
 |                    +-----+-----+-----+
 |                          | Ad  | Au  |
 |                          +-----+-----+-----+
 |                                | Ad  | Au  |
 |  +----+                        +-----+-----+--+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
 |  Applicazione riduzione
 |  ----------------------
 |
 |  Au'' = selezione righe di Au --> / Au'' e Ad'' possono essere compattate
 |  Ad'' = selezione righe di Ad --> \ nella stessa matrice
 |
 |  P * / Ad  0 \ = / Ad'   Bu'  \ ===> P * / Ad \ = / Adu'  \
 |      \ 0  Bu /   \ Ad''  Bu'' /          \ Bu /   \ Adu'' /
 |
 |  where  E* = Ad' - G*Ad'',  F* = Bu' - G*Bu''
 |
 |  La matrice G è spachettata nel blocco
 |
 |  +-----+-----+ - - +
 |  | E*  |  0  | F*   <-- messo al posto dell "0"
 |  +-----+-----+-----+
 |    Ad' | L\U | Au' | <-- memorizzazione compatta
 |  + - - +-----+-----+-----+ - - +
 |              | E*  |  0  | F*
 |              +-----+-----+-----+
 |                Ad' | L\U | Au' |
 |              + - - +-----+-----+-----+ - - +
 |                          | E*  |  0  | F*
 |                          +-----+-----+-----+
 |                            Ad' | L\U | Au' |
 |  +----+                  + - - +-----+-----+--+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
 |  Sistema ridotto e memorizzazione
 |  +-----+-----+
 |  | E*  | F* |
 |  +-----+-----+-----+
 |        | E*  | F*  |
 |        +-----+-----+-----+
 |              | E*  | F*  |
 |  +----+      +-----+-----+--+
 |  | H0 |            | HN  |  |
 |  |    |            |     |  |
 |  +----+            +-----+--+
 |
 |
 |  +-----+-----+
 |  | E*  | F*  |
 |  +-----+-----+-----+                            +===+
 |        | L\U | Adu'|                            | V | <-- memorizzazione compatta e fill in
 |        +-----+-----+-----+                      +===+
 |              | E*  | F*  |
 |              +-----+-----+-----+                +===+
 |                    | L\U | Adu'|                | V | <-- memorizzazione compatta e fill in
 |                    +-----+-----+-----+          +===+
 |                          | E*  | F*  |
 |                          +-----+-----+-----+    +===+
 |                                | L\U | Adu'|    | V | <-- memorizzazione compatta e fill in
 |  +----+                        +-----+-----+--+ +===+
 |  | H0 |                              | HN  |  |
 |  |    |                              |     |  |
 |  +----+                              +-----+--+
 |
\*/

namespace alglin {

  using namespace std ;

  #ifdef LU_BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::AmodioN( integer nth )
  : baseValue("AmodioLU_value")
  , baseInteger("AmodioLU_index")
  , numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= LU_BABD_AMODIO_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << LU_BABD_AMODIO_MAX_THREAD ) ;
  }
  #else
  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::AmodioN()
  : baseValue("AmodioLU_value")
  , baseInteger("AmodioLU_index")
  { }
  #endif

  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::~AmodioN() {
    baseValue   . free() ;
    baseInteger . free() ;
  }

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
  template <typename t_Value, int M, int N, int LDM>
  class Mv {
  public:
    static inline void ass( t_Value const * Mat, t_Value * x, t_Value * res ) {
      res[0] = Dot<t_Value,N,LDM,1>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM>::ass( Mat+1, x, res+1 ) ;
    }
    static inline void addTo( t_Value const * Mat, t_Value * x, t_Value * res ) {
      res[0] += Dot<t_Value,N,LDM,1>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM>::addTo( Mat+1, x, res+1 ) ;
    }
    static inline void subTo( t_Value const * Mat, t_Value * x, t_Value * res ) {
      res[0] -= Dot<t_Value,N,LDM,1>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM>::subTo( Mat+1, x, res+1 ) ;
    }
    static inline void aMxpby( t_Value a, t_Value const * Mat, t_Value * x,
                               t_Value b, t_Value * res ) {
      res[0] = b*res[0] + a*Dot<t_Value,N,LDM,1>::eval( Mat, x ) ;
      Mv<t_Value,M-1,N,LDM>::subTo( Mat+1, x, res+1 ) ;
    }
  } ;

  template <typename t_Value, int N, int LDM>
  class Mv<t_Value,0,N,LDM> {
  public:
    static inline void ass( t_Value const * Mat, t_Value * x, t_Value * res ) { }
    static inline void addTo( t_Value const * Mat, t_Value * x, t_Value * res ) { }
    static inline void subTo( t_Value const * Mat, t_Value * x, t_Value * res ) { }
    static inline void aMxpby( t_Value a, t_Value const * Mat, t_Value * x,
                               t_Value b, t_Value * res ) { }
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

  /*
  //
  // +--------+
  // | \______|
  // | | A    |
  // | |      |
  // +--------+
  // | |      |
  // | | B    |
  // | |      |
  // +--------+
  //
  */
  template <typename t_Value, integer N, integer LEVEL>
  class LU_2_block_inline {
    enum { J = N-LEVEL, NJ = LEVEL, NJ1 = LEVEL-1 } ;
  public:
    static inline integer
    eval( t_Value * A,
          t_Value * B,
          integer   ipiv[] ) {
      integer MX1, MX2 ;
      t_Value absA, absB ;
      t_Value * Ajj = A + J*(N+1) ;
      t_Value * Bj  = B + J*N ;
      Amax<t_Value,NJ,1>::eval( Ajj, absA, MX1 ) ;
      Amax<t_Value,N,1>::eval( Bj, absB, MX2 ) ;
      if ( absA < absB ) {
        ipiv[J] = MX2 + N ; // C-based
        Swap<t_Value,N,N,N>::eval( A+J, B+MX2 ) ;
      } else {
        ipiv[J] = MX1 + J ; // C-based
        if ( MX1 > 0 ) Swap<t_Value,N,N,N>::eval( A+J, A+ipiv[J] ) ;
      }
      if ( Ajj[0] == 0 ) return J+1 ;
      t_Value ROWM = 1/Ajj[0] ;
      Scal<t_Value,NJ1,1>::eval(ROWM,Ajj+1) ;
      Scal<t_Value,N,1>::eval(ROWM,Bj) ;
      //   N, M, LdM, INCX, INCY
      Rank1<t_Value,NJ1,NJ1,N,1,N>::subTo(Ajj+1,Ajj+N,Ajj+N+1) ;
      Rank1<t_Value,N,NJ1,N,1,N>::subTo(Bj,Ajj+N,Bj+N) ;
      return LU_2_block_inline<t_Value,N,LEVEL-1>::eval(A,B,ipiv) ;
    }

  } ;

  template <typename t_Value, integer N>
  class LU_2_block_inline<t_Value,N,0> {
  public:
    static inline integer eval( t_Value * A,
                                t_Value * B,
                                integer   ipiv[] ) {
      return 0 ;
    }
  } ;

  template <typename t_Value, integer n>
  integer
  AmodioN<t_Value,n>::LU_2_block( valuePointer A,
                                  valuePointer B,
                                  integer      ipiv[] ) const {
    
    integer ierr = LU_2_block_inline<t_Value,n,n>::eval(A,B,ipiv) ;
    #if 0
    // LU DECOMPOSITION, ROW INTERCHANGES
    valuePointer Ajj = A ;
    valuePointer Bj  = B ;
    integer MX1, MX2 ;
    t_Value absB ;
    for ( integer j = 0 ; j < n ; Ajj += n+1, Bj += n  ) {
      MX1 = iamax( n-j, Ajj, 1 ) ;
      Amax<t_Value,n,1>::eval( Bj, absB, MX2 ) ;
      if ( std::abs(Ajj[MX1]) < absB ) {
        ipiv[j] = MX2 + n ; // C-based
        Swap<t_Value,n,n,n>::eval( A + j, B + MX2 ) ;
      } else {
        ipiv[j] = MX1 + j ; // C-based
        if ( MX1 > 0 ) Swap<t_Value,n,n,n>::eval( A + j, A + ipiv[j] ) ;
      }
      if ( std::abs(Ajj[0]) == 0 ) return j+1 ;
      valueType ROWM = 1/Ajj[0] ;
      ++j ;
      scal(n-j, ROWM, Ajj+1, 1) ;
      Scal<t_Value,n,1>::eval(ROWM,Bj) ;
      ger(n-j, n-j, -1.0, Ajj+1, 1, Ajj+n, n, Ajj+n+1, n ) ;
      ger(n,   n-j, -1.0, Bj,    1, Ajj+n, n, Bj+n,    n ) ;
    }
    #endif
    /*
    // compute G = B*L^(-1)
    //
    //  / L \ (U) = / I          \ (L*U) =  / I \ (L*U)
    //  \ M /       \ M * L^(-1) /          \ G /
    */
    trsm( SideMultiply::RIGHT,
          ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          n, n, 1.0, A, n, B, n ) ;
    return ierr ;
  }

  /*  
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */

  /*\
   |  Based on the algorithm 
   |  Pierluigi Amodio and Giuseppe Romanazzi
   |  Algorithm 859: BABDCR: a Fortran 90 package for the Solution of Bordered ABD Linear Systems
   |  ACM Transactions on Mathematical Software, 32, 4, 597—608 (2006)
  \*/

  // 0  1  2  3  4  5  6
  // 0  -  2  -  4  -  6* unchanged
  // 0  -  -  -  4  -  6* unchanged
  // 0  -  -  -  -  -  6* special
  // 0  -  -  -  -  -  -

  // 0  1  2  3  4  5  6  7  8  9 10 11 12
  // 0  -  2  -  4  -  6  -  8  - 10  - 12* unchanged
  // 0  -  -  -  4  -  -  -  8  -  -  - 12* unchanged
  // 0  -  -  -  -  -  -  -  8  -  -  - 12* unchanged
  // 0  -  -  -  -  -  -  -  -  -  -  - 12* special
  // 0  -  -  -  -  -  -  -  -  -  -  -

  // 0  1  2  3  4  5  6  7  8  9 10 11 12 13
  // 0  -  2  -  4  -  6  -  8  - 10  - 12  -
  // 0  -  -  -  4  -  -  -  8  -  -  - 12  -
  // 0  -  -  -  -  -  -  -  8  -  -  - 12  -
  // 0  -  -  -  -  -  -  -  -  -  -  - 12  -
  // 0  -  -  -  -  -  -  -  -  -  -  -  -  -

  // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  // 0  -  2  -  4  -  6  -  8  - 10  - 12  - 14  -
  // 0  -  -  -  4  -  -  -  8  -  -  - 12  -  -  -
  // 0  -  -  -  -  -  -  -  8  -  -  -  -  -  -  -
  // 0  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
  
  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::reduction() {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
  
    valuePointer EE = tmpM ;
    valuePointer FF = tmpM+nxn ;

    do {

      if ( jump_block >= nblock ) break ;
      k_block = 0 ;

      do {

        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;
        std::vector<bool> & LU_rows = LU_rows_blk[k1] ;

        valuePointer E   = BLK0 ;       // <-- Ad'' - G*Ad'
        valuePointer F   = BLK0 + nxn ; // <-- Bu'' - G*Bu'
        valuePointer LU  = BLK1 ;
        valuePointer Adu = BLK1 + nxn ;

        /*\
         | factorize RS by means of the LU factorization
         | P * / Au \ = / G \ (L*U)
         |     \ Bd /   \ I /
        \*/
        Copy<t_Value,nxn,1,1>::eval(F,G) ;
        integer info = LU_2_block( LU, G, ipiv ) ;
        ALGLIN_ASSERT( info == 0,
                       "AmodioLU::factorize, at block N." << k <<
                       " singular matrix, info = " << info );


        /*\
         |   +-----+-----+ - - +
         |   | E*  |  0 <=  F*   <-- messo al posto dell "0"
         |   +-----+-----+-----+
         |     Ad''| L\U | Au''| <-- memorizzazione compatta
         |   + - - +-----+-----+
         |
         |   applico permutazione
         |   +-----+             +-----+-----+
         |   | Ad  |             |  E  |  F  |
         |   +-----+-----+  -->  +-----+-----+
         |         | Bu  |       | Ad''+Bu'' | compattato
         |         +-----+       +-----------+
        \*/

        // determine the permutation vector and select row vectors
        for ( integer i = 0 ; i < n ; ++i ) {
          LU_rows[i]   = true ;
          LU_rows[n+i] = false ;
        }
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            std::swap( LU_rows[i], LU_rows[ip] ) ;
            if ( ip < n ) Swap<t_Value,n,n,n>::eval( Adu + i, Adu + ip ) ; // scambia righe
            else          Swap<t_Value,n,n,n>::eval( Adu + i, E + ip-n ) ; // scambia righe
          }
        }
        Zero<t_Value,nxn,1>::eval( F ) ;
        Zero<t_Value,nxnx2,1>::eval( EE ) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            Copy<t_Value,n,n,n>::eval(E+i,F+i) ;
            Zero<t_Value,n,n>::eval(E+i) ;
          }
          if ( LU_rows[i] ) Copy<t_Value,n,n,n>::eval( Adu+i, FF+i ) ;
          else              Copy<t_Value,n,n,n>::eval( Adu+i, EE+i ) ;
        }
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              FF, n,
              1.0, F, n ) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;
      jump_block *= 2 ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::reduction_mt( integer num_thread, integer nth ) {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;

    valuePointer EE = tmpM+nth*nxnx2 ;
    valuePointer FF = EE+nxn ;

    do {

      if ( jump_block >= nblock ) break ;

      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      do {

        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;
        std::vector<bool> & LU_rows = LU_rows_blk[k1] ;

        valuePointer E   = BLK0 ;       // <-- Ad'' - G*Ad'
        valuePointer F   = BLK0 + nxn ; // <-- Bu'' - G*Bu'
        valuePointer LU  = BLK1 ;
        valuePointer Adu = BLK1 + nxn ;

        /*\
         | factorize RS by means of the LU factorization
         | P * / Au \ = / G \ (L*U)
         |     \ Bd /   \ I /
        \*/
        Copy<t_Value,nxn,1,1>::eval( F, G ) ;
        integer info = LU_2_block( LU, G, ipiv ) ;
        ALGLIN_ASSERT( info == 0,
                       "AmodioLU::factorize, at block N." << k <<
                       " singular matrix, info = " << info );


        /*\
         |   +-----+-----+ - - +
         |   | E*  |  0 <=  F*   <-- messo al posto dell "0"
         |   +-----+-----+-----+
         |     Ad''| L\U | Au''| <-- memorizzazione compatta
         |   + - - +-----+-----+
         |
         |   applico permutazione
         |   +-----+             +-----+-----+
         |   | Ad  |             |  E  |  F  |
         |   +-----+-----+  -->  +-----+-----+
         |         | Bu  |       | Ad''+Bu'' | compattato
         |         +-----+       +-----------+
        \*/

        // determine the permutation vector and select row vectors
        for ( integer i = 0 ; i < n ; ++i ) {
          LU_rows[i]   = true ;
          LU_rows[n+i] = false ;
        }
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            std::swap( LU_rows[i], LU_rows[ip] ) ;
            if ( ip < n ) Swap<t_Value,n,n,n>::eval( Adu + i, Adu + ip ) ; // scambia righe
            else          Swap<t_Value,n,n,n>::eval( Adu + i, E + ip-n ) ; // scambia righe
          }
        }
        Zero<t_Value,nxn,1>::eval(F) ;
        Zero<t_Value,nxnx2,1>::eval(EE) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            Copy<t_Value,n,n,n>::eval( E+i, F+i ) ;
            Zero<t_Value,n,n>::eval(E+i) ;
          }
          if ( LU_rows[i] ) Copy<t_Value,n,n,n>::eval( Adu+i, FF+i ) ;
          else              Copy<t_Value,n,n,n>::eval( Adu+i, EE+i ) ;
        }
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( Transposition::NO_TRANSPOSE,
              Transposition::NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              FF, n,
              1.0, F, n ) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
        } else {
          cond0.wait(lck);
        }
      }
    } while ( true ) ;
  }
  #endif

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::factorize( integer           _nblock,
                                 integer           _q,
                                 valueConstPointer AdAu,
                                 valueConstPointer H0,
                                 valueConstPointer HN,
                                 valueConstPointer Hq ) {

    integer const nx2   = n*2 ;
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;

    nblock = _nblock ;
    q      = _q ;
    m      = n+q ;
    nm     = n+m ;

    integer nnzG    = (nblock-1)*nxn ;
    integer nnzADAU = nblock*nxnx2 ;

    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    integer nnzLU = numThread*nxnx2 ;
    if ( nnzLU < nm*nm ) nnzLU = nm*nm ;
    #else
    integer nnzLU = nm*nm ;
    #endif

    integer nv = nnzG + nnzADAU + nnzLU + nm ;
    integer ni = nblock*n + nm ;

    baseValue   . allocate(nv) ;
    baseInteger . allocate(ni) ;

    AdAu_blk = baseValue( nnzADAU ) ;
    G_blk    = baseValue( nnzG ) - nxn ; // 1 based
    LU_blk   = baseValue( nnzLU ) ;
    tmpM     = LU_blk ;
    tmpV     = baseValue( nm ) ;

    ipiv_blk    = baseInteger( nblock*n ) ;
    LU_ipiv_blk = baseInteger( nm ) ;

    LU_rows_blk.resize( nblock ) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      LU_rows_blk[i].resize(nx2) ;

    alglin::copy( nnzADAU, AdAu, 1, AdAu_blk, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    integer usedThread = numThread ;
    #endif

    jump_block = 1 ;
    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    to_be_done = 0 ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &AmodioN<t_Value,n>::reduction_mt, this, usedThread, nt ) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #else
    reduction() ;
    #endif

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    gecopy( n, nx2, AdAu_blk, n, LU_blk, nm ) ;
    gecopy( m, n, H0, m, LU_blk+n,      nm ) ;
    gecopy( m, n, HN, m, LU_blk+n+n*nm, nm ) ;
    if ( _q > 0 ) {
      gezero( n, _q,        LU_blk+nx2*nm,   nm ) ;
      gecopy( m, _q, Hq, m, LU_blk+nx2*nm+n, nm ) ;
    }

    integer INFO = getrf( nm, nm, LU_blk, nm, LU_ipiv_blk ) ;
    ALGLIN_ASSERT( INFO==0,
                   "AmodioLU::factorize(), singular matrix, getrf INFO = " << INFO ) ;

  }

  /*
  //  __                                  _                  _
  // / _| ___  _ ____      ____ _ _ __ __| |    _ __ ___  __| |_   _  ___ ___
  // | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |   | '__/ _ \/ _` | | | |/ __/ _ \
  // |  _| (_) | |   \ V  V / (_| | | | (_| |   | | |  __/ (_| | |_| | (_|  __/
  // |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|___|_|  \___|\__,_|\__,_|\___\___|
  //                                       |_____|
  */

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::forward_reduce( valuePointer y ) const {
    integer const nxn = n*n ;
    do {
      k_block = 0 ;
      if ( jump_block >= nblock ) break ;

      do {
        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;
        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;

        /*\
         |   applico permutazione e moltiplico per / I -G \
         |                                         \    I /
        \*/
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            if ( ip < n ) std::swap( yk1[i], yk1[ip] ) ;
            else          std::swap( yk1[i], yk[ip-n] ) ;
          }
        }
        // yk -= G * yk1
        Mv<t_Value,n,n,n>::subTo(G,yk1,yk) ;
        //gemv( Transposition::NO_TRANSPOSE,
        //      n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
      } while ( true ) ;

      // aspetta le altre thread
      jump_block *= 2 ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::forward_reduce_mt( integer num_thread, valuePointer y ) const {
    integer const nxn = n*n ;
    do {
      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      if ( jump_block >= nblock ) break ;

      do {
        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }
        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;

        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;

        /*\
         |   applico permutazione e moltiplico per / I -G \
         |                                         \    I /
        \*/
        for ( integer i = 0 ; i < n ; ++i ) {
          integer ip = ipiv[i] ;
          if ( ip > i ) {
            if ( ip < n ) std::swap( yk1[i], yk1[ip] ) ;
            else          std::swap( yk1[i], yk[ip-n] ) ;
          }
        }
        // yk -= G * yk1
        Mv<t_Value,n,n,n>::subTo(G,yk1,yk) ;
        //gemv( Transposition::NO_TRANSPOSE,
        //      n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
      } while ( true ) ;

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
        } else {
          cond0.wait(lck);
        }
      }
    } while ( true ) ;
  }
  #endif

  /*
  //   _                _                  _         _   _ _         _
  //  | |__   __ _  ___| | __    ___ _   _| |__  ___| |_(_) |_ _   _| |_ ___
  //  | '_ \ / _` |/ __| |/ /   / __| | | | '_ \/ __| __| | __| | | | __/ _ \
  //  | |_) | (_| | (__|   <    \__ \ |_| | |_) \__ \ |_| | |_| |_| | ||  __/
  //  |_.__/ \__,_|\___|_|\_\___|___/\__,_|_.__/|___/\__|_|\__|\__,_|\__\___|
  //                       |_____|
  */

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::back_substitute( valuePointer y ) const {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    do {
      jump_block /= 2 ;
      k_block    = 0 ;
      if ( jump_block == 0 ) break ;

      do {
        integer k, k1 ;
        // update
        k        = k_block ;
        k_block += 2*jump_block ;

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;
      
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = LU_rows[i] ? yk2 : yk ;
          yk1[i] -= Dot<t_Value,n,n,1>::eval( Adu+i, LR ) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        LsolveUnit<t_Value,n,n>::eval(LU,yk1) ;
        Usolve<t_Value,n,n>::eval(LU,yk1) ;

        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;
    } while ( true ) ;
  }

  #ifdef LU_BABD_AMODIO_N_USE_THREAD

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::back_substitute_mt( integer num_thread, valuePointer y ) const {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    do {
      { unique_lock<mutex> lck(mtx0);
        if ( to_be_done == 0 ) { // prima thread che passa
          jump_block /= 2 ;
          k_block    = 0 ;
          to_be_done = num_thread ;
        }
      }

      if ( jump_block == 0 ) break ;

      do {
        integer k, k1 ;
        // update
        {
          unique_lock<mutex> lck(mtx1);
          k        = k_block ;
          k_block += 2*jump_block ;
        }

        // update....
        k1 = k+jump_block ;
        if ( k1 >= nblock ) break ;
      
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = LU_rows[i] ? yk2 : yk ;
          yk1[i] -= Dot<t_Value,n,n,1>::eval( Adu+i, LR ) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        LsolveUnit<t_Value,n,n>::eval(LU,yk1) ;
        Usolve<t_Value,n,n>::eval(LU,yk1) ;

        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      } while ( true ) ;
      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) cond0.notify_all() ; // wake up all tread
        else                     cond0.wait(lck);
      }
    } while ( true ) ;
  }
  #endif

  /*             _
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::solve( valuePointer y ) const {

    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    integer usedThread = numThread ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    if ( usedThread == 1 ) {
      forward_reduce(y) ;
    } else {
      to_be_done = 0 ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &AmodioN<t_Value,n>::forward_reduce_mt, this, usedThread, y ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #else
    forward_reduce(y) ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! 2 by 2 block linear system solution !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    valuePointer ye = y + nblock * n ;
    Copy<t_Value,n,1,1>::eval( y, tmpV ) ;
    copy( m, ye, 1, tmpV+n, 1 ) ;
    integer INFO = getrs( Transposition::NO_TRANSPOSE,
                          nm, 1, LU_blk, nm, LU_ipiv_blk, tmpV, nm ) ;
    ALGLIN_ASSERT( INFO==0,
                   "AmodioLU::solve(), singular matrix, getrs INFO = " << INFO ) ;
    Copy<t_Value,n,1,1>::eval( tmpV, y ) ;
    copy( m, tmpV+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    #ifdef LU_BABD_AMODIO_N_USE_THREAD
    if ( usedThread == 1 ) {
      back_substitute( y ) ;
    } else {
      to_be_done = 0 ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &AmodioN<t_Value,n>::back_substitute_mt, this, usedThread, y ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #else
    back_substitute( y ) ;
    #endif
  }

  template class AmodioN<double,2> ;
  template class AmodioN<double,3> ;
  template class AmodioN<double,4> ;
  template class AmodioN<double,5> ;
  template class AmodioN<double,6> ;
  template class AmodioN<double,7> ;
  template class AmodioN<double,8> ;
  template class AmodioN<double,9> ;
  template class AmodioN<double,10> ;
  template class AmodioN<double,11> ;
  template class AmodioN<double,12> ;
  template class AmodioN<double,13> ;
  template class AmodioN<double,14> ;
  template class AmodioN<double,15> ;
  template class AmodioN<double,16> ;

}
