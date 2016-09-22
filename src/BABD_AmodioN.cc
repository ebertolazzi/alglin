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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#include "BABD_AmodioN.hh"
#include "Alglin_tmpl.hh"
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

  #ifdef BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::AmodioN( integer nth )
  : baseValue("AmodioN_values")
  , baseInteger("AmodioN_integers")
  , numThread(nth)
  {
    ALGLIN_ASSERT( numThread > 0 && numThread <= BABD_AMODIO_MAX_THREAD,
                   "Bad number of thread specification [" << numThread << "]\n"
                   "must be a number > 0 and <= " << BABD_AMODIO_MAX_THREAD ) ;
  }
  #else
  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::AmodioN()
  : baseValue("AmodioN_values")
  , baseInteger("AmodioN_integers")
  { }
  #endif

  template <typename t_Value, integer n>
  AmodioN<t_Value,n>::~AmodioN() {
    baseValue   . free() ;
    baseInteger . free() ;
  }
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
  AmodioN<t_Value,n>::LU_2_block( valuePointer A,
                                  valuePointer B,
                                  integer      ipiv[] ) const {
    
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

  #ifdef BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::reduction_mt( integer nth ) {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;

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

    #ifdef BABD_AMODIO_N_USE_THREAD
    integer nnzLU = numThread*nxnx2 ;
    if ( nnzLU < nm*nm ) nnzLU = nm*nm ;
    #else
    integer nnzLU = nm*nm ;
    #endif

    integer nv = nnzG + nnzADAU + nnzLU + nm ;
    integer ni = nblock*n + nm ;

    baseValue   . allocate(size_t(nv)) ;
    baseInteger . allocate(size_t(ni)) ;

    AdAu_blk = baseValue(size_t( nnzADAU )) ;
    G_blk    = baseValue(size_t( nnzG )) - nxn ; // 1 based
    LU_blk   = baseValue(size_t( nnzLU )) ;
    tmpM     = LU_blk ;
    tmpV     = baseValue(size_t( nm )) ;

    ipiv_blk    = baseInteger(size_t( nblock*n )) ;
    LU_ipiv_blk = baseInteger(size_t( nm )) ;

    LU_rows_blk.resize( nblock ) ;
    for ( integer i = 0 ; i < nblock ; ++i )
      LU_rows_blk[i].resize(nx2) ;

    alglin::copy( nnzADAU, AdAu, 1, AdAu_blk, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/

    #ifdef BABD_AMODIO_N_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    #endif

    jump_block = 1 ;
    #ifdef BABD_AMODIO_N_USE_THREAD
    // esegue un numero di livelli in parallelo
    to_be_done = usedThread ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt] = std::thread( &AmodioN<t_Value,n>::reduction_mt, this, nt ) ;
    for ( integer nt = 0 ; nt < usedThread ; ++nt )
      threads[nt].join() ;
    #endif
    reduction() ;

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
  //  __                                   _                  _
  // / _|  ___  _ ____      ____ _ _ __ __| |    _ __ ___  __| |_   _  ___ ___
  // | |_ / _ \| '__\ \ /\ / / _` | '__/ _` |   | '__/ _ \/ _` | | | |/ __/ _ \
  // |  _| (_) | |   \ V  V / (_| | | | (_| |   | | |  __/ (_| | |_| | (_|  __/
  // |_|  \___/|_|    \_/\_/ \__,_|_|  \__,_|___|_|  \___|\__,_|\__,_|\___\___|
  //                                       |_____|
  */

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::forward_reduce( valuePointer y ) const {
    integer const nxn = n*n ;

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

  #ifdef BABD_AMODIO_N_USE_THREAD
  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::forward_reduce_mt( integer nth ) const {
    integer const nxn = n*n ;
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
  AmodioN<t_Value,n>::back_substitute( valuePointer y, integer jump_block_min ) const {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;
    while ( jump_block > jump_block_min ) {
      integer k_step = 2*jump_block ;
      integer kend   = nblock-jump_block ;
      for ( integer k = 0 ; k < kend ; k += k_step ) {
        integer      k1  = k+jump_block ;
        integer      k2  = min(k1+jump_block,nblock) ;
        valuePointer LU  = AdAu_blk + k1 * nxnx2 ;
        valuePointer Adu = LU + nxn ;
        std::vector<bool> const & LU_rows = LU_rows_blk[k1] ;

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

  #ifdef BABD_AMODIO_N_USE_THREAD

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::back_substitute_mt( integer nth ) const {
    integer const nxn   = n*n ;
    integer const nxnx2 = n*n*2 ;

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

  /*             _
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */

  template <typename t_Value, integer n>
  void
  AmodioN<t_Value,n>::solve( valuePointer y ) const {

    #ifdef BABD_AMODIO_N_USE_THREAD
    usedThread        = numThread ;
    jump_block_max_mt = nblock>>(usedThread-1) ;
    y_thread          = y ;
    #endif

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! reduction phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block = 1 ;
    #ifdef BABD_AMODIO_N_USE_THREAD
    if ( usedThread > 0 ) {
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &AmodioN<t_Value,n>::forward_reduce_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    }
    #endif
    forward_reduce(y) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! 2 by 2 block linear system solution !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    valuePointer ye = y + nblock * n ;
    Copy<t_Value,n,1,1>::eval( y, tmpV ) ;
    copy( m, ye, 1, tmpV+n, 1 ) ;
    integer INFO = getrs( NO_TRANSPOSE, nm, 1, LU_blk, nm, LU_ipiv_blk, tmpV, nm ) ;
    ALGLIN_ASSERT( INFO==0,
                   "AmodioLU::solve(), singular matrix, getrs INFO = " << INFO ) ;
    Copy<t_Value,n,1,1>::eval( tmpV, y ) ;
    copy( m, tmpV+n, 1, ye, 1 ) ;

    /*\
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     | !!!! back-substitution phase !!!!
     | !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    \*/
    jump_block /= 2 ;
    #ifdef BABD_AMODIO_N_USE_THREAD
    if ( usedThread > 0 ) {
      back_substitute( y, jump_block_max_mt ) ;
      to_be_done = usedThread ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt] = std::thread( &AmodioN<t_Value,n>::back_substitute_mt, this, nt ) ;
      for ( integer nt = 0 ; nt < usedThread ; ++nt )
        threads[nt].join() ;
    } else {
      back_substitute( y, 0 ) ;
    }
    #else
      back_substitute( y, 0 ) ;
    #endif
  }

  // speedup for windows
  #ifndef ALGLIN_OS_WINDOWS
  template class AmodioN<double,1> ;
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
  //template class AmodioN<double,17> ;
  //template class AmodioN<double,18> ;
  //template class AmodioN<double,19> ;
  //template class AmodioN<double,20> ;
  //template class AmodioN<double,21> ;
  //template class AmodioN<double,22> ;
  //template class AmodioN<double,23> ;
  //template class AmodioN<double,24> ;
  //template class AmodioN<double,25> ;
  #endif

}
