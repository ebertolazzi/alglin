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

#include "CyclicReduction.hh"

#ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE

#include "Alglin_tmpl.hh"

#include <iostream>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace alglin {

  using namespace std ;
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
      if ( isZero(Ajj[0]) ) return J+1 ;
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
    static inline integer
    eval( t_Value *,
          t_Value *,
          integer [] ) {
      return 0 ;
    }
  } ;

  template <typename t_Value, integer n>
  integer
  LU_2_block( t_Value A[],
              t_Value B[],
              integer ipiv[] ) {
    integer ierr = LU_2_block_inline<t_Value,n,n>::eval(A,B,ipiv) ;
    /*
    // compute G = B*L^(-1)
    //
    //  / L \ (U) = / I          \ (L*U) =  / I \ (L*U)
    //  \ M /       \ M * L^(-1) /          \ G /
    */
    if ( ierr == 0 ) trsm( RIGHT, LOWER, NO_TRANSPOSE, UNIT,
                           n, n, 1.0, A, n, B, n ) ;
    return ierr ;
  }

  template integer LU_2_block<float,1>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,2>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,3>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,4>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,5>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,6>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,7>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,8>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,9>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,10>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,11>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,12>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,13>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,14>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,15>( float A[], float B[], integer ipiv[] ) ;
  template integer LU_2_block<float,16>( float A[], float B[], integer ipiv[] ) ;

  template integer LU_2_block<double,1>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,2>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,3>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,4>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,5>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,6>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,7>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,8>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,9>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,10>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,11>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,12>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,13>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,14>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,15>( double A[], double B[], integer ipiv[] ) ;
  template integer LU_2_block<double,16>( double A[], double B[], integer ipiv[] ) ;

  template <typename t_Value>
  template <integer n>
  CyclicReduction<t_Value>::FixedSize<n>::FixedSize( CyclicReduction * const _pCR )
  : pCR(_pCR) {}

  /*\
   |                _
   |   _ __ ___  __| |_   _  ___ ___
   |  | '__/ _ \/ _` | | | |/ __/ _ \
   |  | | |  __/ (_| | |_| | (_|  __/
   |  |_|  \___|\__,_|\__,_|\___\___|
  \*/
  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::reduce() {

    integer const n     = N ;
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    
    integer      & jump_block = pCR->jump_block ;
    integer      & nblock     = pCR->nblock ;
    valuePointer & AdAu_blk   = pCR->AdAu_blk ;
    valuePointer & G_blk      = pCR->G_blk ;
    integer *    & ipiv_blk   = pCR->ipiv_blk ;
  
    valuePointer EE = pCR->tmpM ;
    valuePointer FF = pCR->tmpM+nxn ;

    while ( jump_block < nblock ) {
      integer kstep = 2*jump_block ;
      integer kend  = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += kstep ) {
        integer k1 = k+jump_block ;

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;
        std::vector<bool> & LU_rows = pCR->LU_rows_blk[k1] ;

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
        integer info = LU_2_block<t_Value,N>( LU, G, ipiv ) ;
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

        MM<t_Value,n,n,n,n,n,n>::subTo( G, EE, E ) ;
        MM<t_Value,n,n,n,n,n,n>::subTo( G, FF, F ) ;

        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      }
      jump_block *= 2 ;
    }

  }

  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::reduce_mt( integer nth ) {
    integer const n     = N ;
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;

    integer      & usedThread        = pCR->usedThread ;
    integer      & jump_block        = pCR->jump_block ;
    integer      & jump_block_max_mt = pCR->jump_block_max_mt ;
    integer      & nblock            = pCR->nblock ;
    valuePointer & AdAu_blk          = pCR->AdAu_blk ;
    valuePointer & G_blk             = pCR->G_blk ;
    integer *    & ipiv_blk          = pCR->ipiv_blk ;

    valuePointer EE = pCR->tmpM+nth*nxnx2 ;
    valuePointer FF = EE+nxn ;

    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;

        valuePointer BLK0 = AdAu_blk + k  * nxnx2 ;
        valuePointer BLK1 = AdAu_blk + k1 * nxnx2 ;
        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;
        std::vector<bool> & LU_rows = pCR->LU_rows_blk[k1] ;

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
        integer info = LU_2_block<t_Value,N>( LU, G, ipiv ) ;
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

        MM<t_Value,n,n,n,n,n,n>::subTo( G, EE, E ) ;
        MM<t_Value,n,n,n,n,n,n>::subTo( G, FF, F ) ;

        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(pCR->mtx0);
        if ( --pCR->to_be_done == 0 ) {
          pCR->cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          pCR->to_be_done = usedThread ;
        } else {
          pCR->cond0.wait(lck);
        }
      }
    }
  }
  #endif

  /*\
   |    __                                  _
   |   / _| ___  _ ____      ____ _ _ __ __| |
   |  | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |
   |  |  _| (_) | |   \ V  V / (_| | | | (_| |
   |  |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|
  \*/
  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::forward( valuePointer y ) const {

    integer const n   = N ;
    integer const nxn = n*n ;
    
    integer      & jump_block = pCR->jump_block ;
    integer      & nblock     = pCR->nblock ;
    valuePointer & G_blk      = pCR->G_blk ;
    integer *    & ipiv_blk   = pCR->ipiv_blk ;

    while ( jump_block < nblock ) {

      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;

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
        Mv<t_Value,n,n,n,1,1>::subTo(G,yk1,yk) ;
      }

      // aspetta le altre thread
      jump_block *= 2 ;
    }
  }
  
  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::forward_mt( integer nth ) const {

    integer const n   = N ;
    integer const nxn = n*n ;
    
    integer      & jump_block        = pCR->jump_block ;
    integer      & jump_block_max_mt = pCR->jump_block_max_mt ;
    integer      & nblock            = pCR->nblock ;
    valuePointer & G_blk             = pCR->G_blk ;
    integer *    & ipiv_blk          = pCR->ipiv_blk ;
    valuePointer & y_thread          = pCR->y_thread ;
    integer      & usedThread        = pCR->usedThread ;
    integer      & to_be_done        = pCR->to_be_done ;

    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {
        integer k1 = k+jump_block ;

        valuePointer G    = G_blk    + k1 * nxn ;
        integer *    ipiv = ipiv_blk + k1 * n ;

        valuePointer yk  = y_thread + k  * n ;
        valuePointer yk1 = y_thread + k1 * n ;

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
        Mv<t_Value,n,n,n,1,1>::subTo(G,yk1,yk) ;
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(pCR->mtx0);
        if ( --to_be_done == 0 ) {
          pCR->cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          to_be_done = usedThread ;
        } else {
          pCR->cond0.wait(lck);
        }
      }
    }
  }
  #endif

  /*\
   |   _                _                           _
   |  | |__   __ _  ___| | ____      ____ _ _ __ __| |
   |  | '_ \ / _` |/ __| |/ /\ \ /\ / / _` | '__/ _` |
   |  | |_) | (_| | (__|   <  \ V  V / (_| | | | (_| |
   |  |_.__/ \__,_|\___|_|\_\  \_/\_/ \__,_|_|  \__,_|
  \*/
  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::backward( valuePointer y, integer jump_block_min ) const {

    integer const n     = N ;
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    
    integer      & jump_block = pCR->jump_block ;
    integer      & nblock     = pCR->nblock ;
    valuePointer & AdAu_blk   = pCR->AdAu_blk ;

    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer      k1  = k+jump_block ;
        integer      k2  = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = pCR->LU_rows_blk[k1] ;

        valuePointer yk1 = y + k1 * n ;
        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = y + ( LU_rows[i] ? k2 : k ) * n ;
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
      }
      jump_block /= 2 ;
    }
  }

  #ifdef CYCLIC_REDUCTION_USE_THREAD

  template <typename t_Value>
  template <integer N>
  void
  CyclicReduction<t_Value>::FixedSize<N>::backward_mt( integer nth ) const {

    integer const n     = N ;
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    
    integer      & jump_block = pCR->jump_block ;
    integer      & nblock     = pCR->nblock ;
    integer      & usedThread = pCR->usedThread ;
    integer      & to_be_done = pCR->to_be_done ;
    valuePointer & AdAu_blk   = pCR->AdAu_blk ;
    valuePointer & y_thread   = pCR->y_thread ;

    while ( jump_block > 0 ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = pCR->LU_rows_blk[k1] ;

        valuePointer yk1 = y_thread + k1 * n ;
        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = y_thread + ( LU_rows[i] ? k2 : k ) * n ;
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
      }
      // aspetta le altre thread
      { unique_lock<mutex> lck(pCR->mtx0);
        if ( --to_be_done == 0 ) {
          pCR->cond0.notify_all() ; // wake up all tread
          jump_block /= 2 ;
          pCR->to_be_done = usedThread ;
        } else {
          pCR->cond0.wait(lck);
        }
      }
    }  
  }
  #endif

  template class CyclicReduction<float>::FixedSize<1> ;
  template class CyclicReduction<double>::FixedSize<1> ;

  template class CyclicReduction<float>::FixedSize<2> ;
  template class CyclicReduction<double>::FixedSize<2> ;

  template class CyclicReduction<float>::FixedSize<3> ;
  template class CyclicReduction<double>::FixedSize<3> ;

  template class CyclicReduction<float>::FixedSize<4> ;
  template class CyclicReduction<double>::FixedSize<4> ;

  template class CyclicReduction<float>::FixedSize<5> ;
  template class CyclicReduction<double>::FixedSize<5> ;

  template class CyclicReduction<float>::FixedSize<6> ;
  template class CyclicReduction<double>::FixedSize<6> ;

  template class CyclicReduction<float>::FixedSize<7> ;
  template class CyclicReduction<double>::FixedSize<7> ;

  template class CyclicReduction<float>::FixedSize<8> ;
  template class CyclicReduction<double>::FixedSize<8> ;

  template class CyclicReduction<float>::FixedSize<9> ;
  template class CyclicReduction<double>::FixedSize<9> ;

  template class CyclicReduction<float>::FixedSize<10> ;
  template class CyclicReduction<double>::FixedSize<10> ;
}

#endif
