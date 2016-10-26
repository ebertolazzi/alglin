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
 |                                +-----+-----+
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
 |                          + - - +-----+-----+
 |
 |  Sistema ridotto e memorizzazione
 |  +-----+-----+
 |  | E*  | F* |
 |  +-----+-----+-----+
 |        | E*  | F*  |
 |        +-----+-----+-----+
 |              | E*  | F*  |
 |              +-----+-----+
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
 |                                +-----+-----+    +===+
 |
\*/

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

namespace alglin {

  using namespace std ;

  /*\
   |   +--------+
   |   | \______|
   |   | | A    |
   |   | |      |
   |   +--------+
   |   | |      |
   |   | | B    |
   |   | |      |
   |   +--------+
  \*/
  template <typename t_Value>
  integer
  LU_2_block( integer N,
              t_Value A[],
              t_Value B[],
              integer ipiv[],
              integer NB ) {
    // LU DECOMPOSITION, ROW INTERCHANGES
    t_Value * Ajj = A ;
    t_Value * Bj  = B ;
    integer nb = min_index(NB,N) ;
    for ( integer j = 0 ; j < N ; Ajj += N+1, Bj += N  ) {
      integer MX1 = iamax( N-j, Ajj, 1 ) ;
      integer MX2 = iamax( N,   Bj,  1 ) ;
      if ( std::abs(Ajj[MX1]) < std::abs(Bj[MX2]) ) {
        ipiv[j] = MX2 + N ; // C-based
        swap( N, A + j, N, B + MX2, N ) ;
      } else {
        ipiv[j] = MX1 + j ; // C-based
        if ( MX1 > 0 ) swap( N, A + j, N, A + ipiv[j], N ) ;
      }
      if ( isZero(Ajj[0]) ) return j+1 ;
      t_Value ROWM = 1/Ajj[0] ;
      ++j ;
      scal(N-j, ROWM, Ajj+1, 1) ;
      scal(N,   ROWM, Bj,    1) ;
      if ( nb == j ) { // applico al prox blocco
        /*
        // / L 0 \-1   /  L^(-1)      \
        // \ M I /   = \ -M*L^(-1)  I /
        */
        nb += NB ;
        if ( nb > N-NB/2 ) nb = N ;
        trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT,
              j, nb-j, 1.0, A, N, A+j*N, N ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              N-j, nb-j, j,
              -1.0, A+j,       N,
                    A+j*N,     N,
               1.0, A+j*(N+1), N ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              N, nb-j, j,
              -1.0, B,     N,
                    A+j*N, N,
               1.0, B+j*N, N ) ;
      } else {
        ger(N-j, nb-j, -1.0, Ajj+1, 1, Ajj+N, N, Ajj+N+1, N ) ;
        ger(N,   nb-j, -1.0, Bj,    1, Ajj+N, N, Bj+N,    N ) ;
      }
    }
    /*
    // compute G = B*L^(-1)
    //
    //  / L \ (U) = / I          \ (L*U) =  / I \ (L*U)
    //  \ M /       \ M * L^(-1) /          \ G /
    */
    trsm( RIGHT, LOWER, NO_TRANSPOSE, UNIT,
          N, N, 1.0, A, N, B, N ) ;
    return 0 ;
  }

  /*\
   |   ____           _ _        ____          _            _   _
   |  / ___|   _  ___| (_) ___  |  _ \ ___  __| |_   _  ___| |_(_) ___  _ __
   | | |  | | | |/ __| | |/ __| | |_) / _ \/ _` | | | |/ __| __| |/ _ \| '_ \
   | | |__| |_| | (__| | | (__  |  _ <  __/ (_| | |_| | (__| |_| | (_) | | | |
   |  \____\__, |\___|_|_|\___| |_| \_\___|\__,_|\__,_|\___|\__|_|\___/|_| |_|
   |       |___/
  \*/

  template <typename t_Value>
  #ifdef CYCLIC_REDUCTION_USE_THREAD
  CyclicReduction<t_Value>::CyclicReduction( integer nth )
  #else
  CyclicReduction<t_Value>::CyclicReduction()
  #endif
  : baseValue("CyclicReduction_value")
  , baseInteger("CyclicReduction_index")
  #ifdef CYCLIC_REDUCTION_USE_THREAD
  , numThread(nth)
  #endif
  #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
  , fixed2(this)
  , fixed3(this)
  , fixed4(this)
  , fixed5(this)
  , fixed6(this)
  , fixed7(this)
  , fixed8(this)
  , fixed9(this)
  , fixed10(this)
  #endif
  , NB(25)
  {
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    ALGLIN_ASSERT( numThread > 0 && numThread <= CYCLIC_REDUCTION_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << CYCLIC_REDUCTION_MAX_THREAD ) ;
    #endif
  }

  template <typename t_Value>
  CyclicReduction<t_Value>::~CyclicReduction() {
    baseValue   . free() ;
    baseInteger . free() ;
  }
  
  /*\
   |         _ _                 _
   |    __ _| | | ___   ___ __ _| |_ ___
   |   / _` | | |/ _ \ / __/ _` | __/ _ \
   |  | (_| | | | (_) | (_| (_| | ||  __/
   |   \__,_|_|_|\___/ \___\__,_|\__\___|
  \*/

  template <typename t_Value>
  void
  CyclicReduction<t_Value>::allocate( integer _nblock, integer _n ) {

    nblock = _nblock ;
    n      = _n ;
    nx2    = 2*n ;
    nxn    = n*n ;
    nxnx2  = nxn*2 ;

    integer nnzG    = (nblock-1)*nxn ;
    integer nnzADAU = nblock*nxnx2 ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    integer nnzLU = numThread*nxnx2 ;
    #else
    integer nnzLU = n*n ;
    #endif

    integer nv = nnzG + nnzADAU + nnzLU ;
    integer ni = nblock*n ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    nv += numThread*nxnx2 ;
    #else
    nv += nxnx2 ;
    #endif

    baseValue   . allocate(size_t(nv)) ;
    baseInteger . allocate(size_t(ni)) ;

    AdAu_blk = baseValue(size_t( nnzADAU )) ;
    G_blk    = baseValue(size_t( nnzG )) - nxn ; // 1 based

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    tmpM = baseValue(size_t( numThread*nxnx2 )) ;
    #else
    tmpM = baseValue(size_t( nxnx2 )) ;
    #endif

    ipiv_blk = baseInteger(size_t( nblock*n )) ;

    LU_rows_blk.resize( nblock ) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      LU_rows_blk[i].resize(nx2) ;
  }

  /*\
   |                _
   |   _ __ ___  __| |_   _  ___ ___
   |  | '__/ _ \/ _` | | | |/ __/ _ \
   |  | | |  __/ (_| | |_| | (_|  __/
   |  |_|  \___|\__,_|\__,_|\___\___|
  \*/
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

  template <typename t_Value>
  void
  CyclicReduction<t_Value>::reduce() {
    jump_block = 1 ;

    #ifdef CYCLIC_REDUCTION_USE_THREAD
    to_be_done = usedThread = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt ) {
      #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
      switch ( n ) {
      case 2:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<2>::reduce_mt, &fixed2, nt ) ; break ;
      case 3:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<3>::reduce_mt, &fixed3, nt ) ; break ;
      case 4:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<4>::reduce_mt, &fixed4, nt ) ; break ;
      case 5:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<5>::reduce_mt, &fixed5, nt ) ; break ;
      case 6:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<6>::reduce_mt, &fixed6, nt ) ; break ;
      case 7:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<7>::reduce_mt, &fixed7, nt ) ; break ;
      case 8:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<8>::reduce_mt, &fixed8, nt ) ; break ;
      case 9:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<9>::reduce_mt, &fixed9, nt ) ; break ;
      case 10: threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<10>::reduce_mt, &fixed10, nt ) ; break ;
      default: threads[nt] = std::thread( &CyclicReduction<t_Value>::reduce_mt, this, nt ) ;
      }
      #else
      threads[nt] = std::thread( &CyclicReduction<t_Value>::reduce_mt, this, nt ) ;
      #endif
    }
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #endif

    #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
    switch ( n ) {
    case 2:  fixed2.reduce() ; return ;
    case 3:  fixed3.reduce() ; return ;
    case 4:  fixed4.reduce() ; return ;
    case 5:  fixed5.reduce() ; return ;
    case 6:  fixed6.reduce() ; return ;
    case 7:  fixed7.reduce() ; return ;
    case 8:  fixed8.reduce() ; return ;
    case 9:  fixed9.reduce() ; return ;
    case 10: fixed10.reduce() ; return ;
    }
    #endif

    valuePointer EE = tmpM ;
    valuePointer FF = tmpM+nxn ;

    while ( jump_block < nblock ) {

      integer kstep = 2*jump_block ;
      integer kend  = nblock-jump_block ;

      for ( integer k = 0 ; k < kend ; k += kstep ) {

        integer k1 = k+jump_block ;

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
        copy( nxn, F, 1, G, 1 ) ;
        integer info = LU_2_block( n, LU, G, ipiv, NB ) ;
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
            if ( ip < n ) swap( n, Adu + i, n, Adu + ip, n ) ; // scambia righe
            else          swap( n, Adu + i, n, E + ip-n, n ) ; // scambia righe
          }
        }
        zero( nxn, F, 1 ) ;
        zero( nxnx2, EE, 1 ) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            copy( n, E+i, n, F+i, n ) ;
            zero( n, E+i, n ) ;
          }
          if ( LU_rows[i] ) copy( n, Adu+i, n, FF+i, n ) ;
          else              copy( n, Adu+i, n, EE+i, n ) ;
        }
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
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
      }
      jump_block *= 2 ;
    }
  }

  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename t_Value>
  void
  CyclicReduction<t_Value>::reduce_mt( integer nth ) {
    valuePointer EE = tmpM+nth*nxnx2 ;
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
        copy( nxn, F, 1, G, 1 ) ;
        integer info = LU_2_block( n, LU, G, ipiv, NB ) ;
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
            if ( ip < n ) swap( n, Adu + i, n, Adu + ip, n ) ; // scambia righe
            else          swap( n, Adu + i, n, E + ip-n, n ) ; // scambia righe
          }
        }
        zero( nxn, F, 1 ) ;
        zero( nxnx2, EE, 1 ) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          if ( LU_rows[i+n] ) {
            copy( n, E+i, n, F+i, n ) ;
            zero( n, E+i, n ) ;
          }
          if ( LU_rows[i] ) copy( n, Adu+i, n, FF+i, n ) ;
          else              copy( n, Adu+i, n, EE+i, n ) ;
        }
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
              n, n, n,
              -1.0, G, n,
              EE, n,
              1.0, E, n ) ;
        gemm( NO_TRANSPOSE, NO_TRANSPOSE,
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
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
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
  void
  CyclicReduction<t_Value>::forward( valuePointer y ) const {
    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 0 ) {
      y_thread = y ;
      to_be_done = usedThread = numThread ;
      jump_block_max_mt = nblock>>(usedThread-1) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt ) {
        #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
        switch ( n ) {
        case 2:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<2>::forward_mt, &fixed2, nt ) ; break ;
        case 3:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<3>::forward_mt, &fixed3, nt ) ; break ;
        case 4:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<4>::forward_mt, &fixed4, nt ) ; break ;
        case 5:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<5>::forward_mt, &fixed5, nt ) ; break ;
        case 6:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<6>::forward_mt, &fixed6, nt ) ; break ;
        case 7:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<7>::forward_mt, &fixed7, nt ) ; break ;
        case 8:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<8>::forward_mt, &fixed8, nt ) ; break ;
        case 9:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<9>::forward_mt, &fixed9, nt ) ; break ;
        case 10: threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<10>::forward_mt, &fixed10, nt ) ; break ;
        default: threads[nt] = std::thread( &CyclicReduction<t_Value>::forward_mt, this, nt ) ;
        }
        #else
        threads[nt] = std::thread( &CyclicReduction<t_Value>::forward_mt, this, nt ) ;
        #endif
      }
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #endif
 
    #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
    switch ( n ) {
    case 2:  fixed2.forward(y) ; return ;
    case 3:  fixed3.forward(y) ; return ;
    case 4:  fixed4.forward(y) ; return ;
    case 5:  fixed5.forward(y) ; return ;
    case 6:  fixed6.forward(y) ; return ;
    case 7:  fixed7.forward(y) ; return ;
    case 8:  fixed8.forward(y) ; return ;
    case 9:  fixed9.forward(y) ; return ;
    case 10: fixed10.forward(y) ; return ;
    }
    #endif

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
        gemv( NO_TRANSPOSE, n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
      }

      // aspetta le altre thread
      jump_block *= 2 ;
    }
  }

  #ifdef CYCLIC_REDUCTION_USE_THREAD
  template <typename t_Value>
  void
  CyclicReduction<t_Value>::forward_mt( integer nth ) const {
    while ( jump_block < jump_block_max_mt ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer     k1 = k+jump_block ;
        valuePointer G = G_blk    + k1 * nxn ;
        integer * ipiv = ipiv_blk + k1 * n ;

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
        gemv( NO_TRANSPOSE, n, n, -1.0, G, n, yk1, 1, 1.0, yk, 1 ) ;
      }

      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block *= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
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
  void
  CyclicReduction<t_Value>::backward( valuePointer y, integer jump_block_min ) const {
    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer      k1  = k+jump_block ;
        integer      k2  = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk  = y + k  * n ;
        valuePointer yk1 = y + k1 * n ;
        valuePointer yk2 = y + k2 * n ;

        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = LU_rows[i] ? yk2 : yk ;
          yk1[i] -= dot( n, Adu+i, n, LR, 1 ) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, LU, n, yk1, 1 ) ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, LU, n, yk1, 1 ) ;
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
  void
  CyclicReduction<t_Value>::backward_mt( integer nth ) const {

    while ( jump_block > 0 ) {

      integer k_step = 2*usedThread*jump_block ;
      integer k0     = 2*nth*jump_block ;
      integer kend   = nblock-jump_block ;

      for ( integer k = k0 ; k < kend ; k += k_step ) {

        integer k1 = k+jump_block ;
        integer k2 = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

        valuePointer yk1 = y_thread + k1 * n ;
        for ( integer i = 0 ; i < n ; ++i ) {
          valuePointer LR = y_thread + ( LU_rows[i] ? k2 : k ) * n ;
          yk1[i] -= dot( n, Adu+i, n, LR, 1 ) ;
        }
        // (LU)^(-1) = U^(-1) L^(-1)
        trsv( LOWER, NO_TRANSPOSE, UNIT, n, LU, n, yk1, 1 ) ;
        trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, LU, n, yk1, 1 ) ;
        /*\
         |  +-----+-----+ - - +
         |  | E*  |  0  | F*    <-- messo al posto dell "0"
         |  +-----+-----+-----+
         |    Ad' | L\U | Bu' | <-- memorizzazione compatta
         |  + - - +-----+-----+-----+ - - +
        \*/
      }
      // aspetta le altre thread
      { unique_lock<mutex> lck(mtx0);
        if ( --to_be_done == 0 ) {
          cond0.notify_all() ; // wake up all tread
          jump_block /= 2 ;
          to_be_done = usedThread ;
        } else {
          cond0.wait(lck);
        }
      }
    }
  }
  #endif

  template <typename t_Value>
  void
  CyclicReduction<t_Value>::backward( valuePointer y ) const {
    for ( jump_block = 1 ; jump_block < nblock ; jump_block *= 2 ) {}
    jump_block /= 2 ;
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    if ( usedThread > 0 ) {
      backward( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt ) {
        #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
        switch ( n ) {
        case 2:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<2>::backward_mt, &fixed2, nt ) ; break ;
        case 3:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<3>::backward_mt, &fixed3, nt ) ; break ;
        case 4:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<4>::backward_mt, &fixed4, nt ) ; break ;
        case 5:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<5>::backward_mt, &fixed5, nt ) ; break ;
        case 6:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<6>::backward_mt, &fixed6, nt ) ; break ;
        case 7:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<7>::backward_mt, &fixed7, nt ) ; break ;
        case 8:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<8>::backward_mt, &fixed8, nt ) ; break ;
        case 9:  threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<9>::backward_mt, &fixed9, nt ) ; break ;
        case 10: threads[nt] = std::thread( &CyclicReduction<t_Value>::FixedSize<10>::backward_mt, &fixed10, nt ) ; break ;
        default: threads[nt] = std::thread( &CyclicReduction<t_Value>::backward_mt, this, nt ) ;
        }
        #else
        threads[nt] = std::thread( &CyclicReduction<t_Value>::backward_mt, this, nt ) ;
        #endif
      }
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
    #endif
      #ifdef CYCLIC_REDUCTION_USE_FIXED_SIZE
      switch ( n ) {
      case 2:  fixed2.backward( y, 0 ) ; break ;
      case 3:  fixed3.backward( y, 0 ) ; break ;
      case 4:  fixed4.backward( y, 0 ) ; break ;
      case 5:  fixed5.backward( y, 0 ) ; break ;
      case 6:  fixed6.backward( y, 0 ) ; break ;
      case 7:  fixed7.backward( y, 0 ) ; break ;
      case 8:  fixed8.backward( y, 0 ) ; break ;
      case 9:  fixed9.backward( y, 0 ) ; break ;
      case 10: fixed10.backward( y, 0 ) ; break ;
      default: backward( y, 0 ) ;
      }
      #else
      backward( y, 0 ) ;
      #endif
    #ifdef CYCLIC_REDUCTION_USE_THREAD
    }
    #endif
  }

  template class CyclicReduction<double> ;
  template class CyclicReduction<float> ;
  
  template integer LU_2_block( integer N, double A[], double B[], integer ipiv[], integer NB ) ;
  template integer LU_2_block( integer N, float A[], float B[], integer ipiv[], integer NB ) ;

}
