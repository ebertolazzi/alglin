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

#include "LU_BABD_Amodio.hh"

#include <iostream>
using namespace std ;


namespace alglin {

  using namespace std ;

  /*  
  //    __            _             _         
  //   / _| __ _  ___| |_ ___  _ __(_)_______ 
  //  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
  //  |  _| (_| | (__| || (_) | |  | |/ /  __/
  //  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  */
  template <typename t_Value>
  void
  AmodioLU<t_Value>::factorize( integer           nblock,
                                integer           n,
                                integer           q,
                                valueConstPointer AdAu,
                                valueConstPointer H0,
                                valueConstPointer HN,
                                valueConstPointer Hq ) {

    this -> nblock = nblock ;
    this -> n      = n ;
    this -> m      = n+q ;

    integer nnzF   = nblock*n*n ;
    integer nnzSR  = 2*nnzF ;
    integer nnzD   = (n+m)*(n+m) ;
    integer nnzTMP = 2*n*n ;

    integer nv = nnzF + nnzSR + nnzD + nnzTMP ;
    integer ni = 2*n*nblock + 2*n+m ;

    baseValue   . allocate(nv) ;
    baseInteger . allocate(ni) ;

    SR_blk    = baseValue( nnzSR  ) ;
    D_blk     = baseValue( nnzD   ) ;
    F_blk     = baseValue( nnzF   ) ;
    TMP_blk   = baseValue( nnzTMP ) ;

    perm_blk  = baseInteger( 2*n*nblock ) ;
    ipiv_blk  = baseInteger( n+m ) ;
    ipiv_work = baseInteger( n ) ;

    alglin::copy( nnzSR, AdAu, 1, SR_blk, 1 ) ;

    /*\
     |  Based on the algorithm 
     |  Pierluigi Amodio and Giuseppe Romanazzi
     |  Algorithm 859: BABDCR: a Fortran 90 package for the Solution of Bordered ABD Linear Systems
     |  ACM Transactions on Mathematical Software, 32, 4, 597—608 (2006)
    \*/
    
    // some constanst
    integer const nxn = n*n ;
    integer const nx2 = n*2 ;

    // initialize indices
    valuePointer F    = F_blk ;
    integer *    perm = perm_blk ;
    
    // 0  1  2  3  4  5  6  7  8  9 10 11 12
    // 0  *  2  *  4  *  6  *  8  * 10  * 12
    // 0  -  *  -  4  -  *  -  8  -  *  - 12
    // 0  -  -  -  *  -  -  -  8  -  -  - 12
    // 0  -  -  -  -  -  -  -  *  -  -  - 12    

    // !!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! reduction phase !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!
    for ( integer jump = 1 ; jump < nblock ; jump *= 2 ) {

      for ( integer k = 1 ; k*jump < nblock ; k += 2 ) {
      
        integer      jC = k*jump ;
        integer      jL = jC-jump ;
        integer      jR = min(jC + jump,nblock) ;
        valuePointer S  = SR_blk +  2*jL*nxn ;
        valuePointer RS = SR_blk + (2*jC-1)*nxn ;
        valuePointer R  = SR_blk + (2*jR-1)*nxn ;

        // reshape ( Rb  Sb ) as ( Rb ) and save in ( , ) columns by columns
        //                       ( Sb )
        copy( 2*nxn, RS, 1, TMP_blk, 1 ) ;
        gecopy( n, n, TMP_blk,     n, RS,   nx2 ) ;
        gecopy( n, n, TMP_blk+nxn, n, RS+n, nx2 ) ;

        /*\
         |  reduces a 2x3 block matrix / Sa  Rb     \  to a 1x2 block matrix ( Sa' Rc' )
         |                             \     Sb  Rc /
         |  
         |  by using the following algorithm:
         |
         |    P * / Rb \ = / Rb' \  = / L          \ (U) = / I              \ (L*U) = / I \ (L*U)
         |        \ Sb /   \ Sb' /    \ Sb' U^(-1) /       \ Sb' (L*U)^(-1) /         \ G /
         |
         |  where  and  G = Sb'*U^(-1)*L^(-1).
         |
         |    / -G  I \ * P * / Rb \ = / -G  I \ / I \ (L*U) = /  0  \
         |    \  I  0 /       \ Sb /   \  I  0 / \ G /         \ L*U /
         |
         | ------------------------------------------------------------
         |
         |    / -G  I \ * P * / Sa   0  \ = / -G  I \ / Sa'  Rc'  \ = / Sa''' Rc''' \ 
         |    \  I  0 /       \ 0    Rc /   \  I  0 / \ Sa'' Rc'' /   \ Sa'   Rc'   /
         |
         |  where  Sa''' = Sa'' - G*Sa', Rc''' = Rc'' - G*Rc' and thus
         |
         |    / -G  I \ * P * / Sa  Rb     \ = / Sa'''  0    Rc''' \
         |    \  I  0 /       \     Sb  Rc /   \ Sa'    L*U  Rc'   /
         |  
        \*/
        /*
        // factorize RS by means of the LU factorization
        // P * / Rb \ = / L          \ (U)
        //     \ Sb /   \ Sb' U^(-1) /
        */
        integer INFO = getrf( nx2, n, RS, nx2, ipiv_work ) ;
        ALGLIN_ASSERT( INFO==0, "AmodioLU::factorize, singular matrix" );

        /*
        // compute G = Sb' U^(-1) L^(-1)
        //  / L          \ (U) = / I                 \ (L*U) =  / I \ (L*U)
        //  \ Sb' U^(-1) /       \ Sb' U^(-1) L^(-1) /          \ G /
        */
        trsm( SideMultiply::RIGHT,
              ULselect::LOWER,
              Transposition::NO_TRANSPOSE,
              DiagonalType::UNIT,
              n, n, 1.0, RS, nx2, RS+n, nx2 ) ;

        // determine the permutation vector
        for ( integer i = 0 ; i < nx2 ; ++i ) perm[i] = i ;
        for ( integer i = 0 ; i < n   ; ++i ) {
          integer ip = ipiv_work[i]-1 ;
          if ( ip != i ) std::swap( perm[i], perm[ip] ) ;
        }

        /*
        // P * / Sa   0  \ = / Sa'  Rc'  \ 
        //     \ 0    Rc /   \ Sa'' Rc'' /
        // S = Sa'' and R = Rc''
        */
        valuePointer Sa1 = TMP_blk ;
        valuePointer Rc1 = TMP_blk + nxn ;
        copy( nxn, S, 1, Sa1, 1 ) ; zero( nxn, S, 1 ) ;
        copy( nxn, R, 1, Rc1, 1 ) ; zero( nxn, R, 1 ) ;
        for ( integer i = 0 ; i < n ; ++i ) {
          integer pi = perm[i+n] ;
          if ( pi < n ) copy( n, Sa1+pi,   n, S+i, n ) ;
          else          copy( n, Rc1+pi-n, n, R+i, n ) ;
        }

        /*
        // / -G  I \ / Sa'  Rc'  \ = / Sa''' Rc''' \ 
        // \  I  0 / \ Sa'' Rc'' /   \ Sa'   Rc'   /
        //
        // where  Sa''' = Sa'' - G*Sa', Rc''' = Rc'' - G*Rc' 
        // where the nonnull rows of Sa' and Rc' are saved in F
        */
        for ( integer i = 0 ; i < n ; ++i ) {
          integer      pi = perm[i] ;
          valuePointer Fi = F  + n*i ;
          valuePointer Gi = RS + (2*i+1)*n ; // i-th column of G
          valuePointer R_or_S, Ra_or_Sc ;
          if ( pi < n ) { // save nonnull row of Sa' in the column of F
            Ra_or_Sc = Sa1 + pi ;
            R_or_S   = S ;
          } else { // save nonnull row of Rc' in the column of F
            Ra_or_Sc = Rc1 + pi - n ;
            R_or_S   = R ;
          }
          copy( n, Ra_or_Sc, n, Fi, 1 ) ;
          ger( n, n, -1.0, Gi, 1, Fi, 1, R_or_S, n ) ;
        }

        perm += nx2 ;
        F    += nxn ;
      }
    }

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last 2 by 2 matrix !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    valuePointer S = SR_blk ;
    valuePointer R = SR_blk + (2*nblock-1)*nxn ;

    integer nm = n+m ;
    gecopy( n, n, S,  n, D_blk,        nm ) ;
    gecopy( n, n, R,  n, D_blk+n*nm,   nm ) ;

    gecopy( m, n, H0, m, D_blk+n,      nm ) ;
    gecopy( m, n, HN, m, D_blk+n+n*nm, nm ) ;

    if ( m > n ) {
      gezero( n, m-n,        D_blk+nx2*nm,   nm ) ;
      gecopy( m, m-n, Hq, m, D_blk+nx2*nm+n, nm ) ;
    }

    integer INFO = getrf( nm, nm, D_blk, nm, ipiv_blk ) ;
    ALGLIN_ASSERT( INFO==0, "AmodioLU::factorize(), singular matrix" ) ;

  }
  
  /*             _           
  //   ___  ___ | |_   _____ 
  //  / __|/ _ \| \ \ / / _ \
  //  \__ \ (_) | |\ V /  __/
  //  |___/\___/|_| \_/ \___|
  */
  template <typename t_Value>
  void
  AmodioLU<t_Value>::solve( valuePointer y ) {
    /*\
     |
     |  Purpose
     |  =======
     |
     |  BABDCR_SOLV  solves a babd linear system whose coefficient matrix
     |  has been factorized by BABDCR_FACT. The algorithm consists of three
     |  phases: reduction (by using the subroutine REDUCE_RHS), solution of
     |  the 2 by 2 block system (by using the subroutine DGETRS) and 
     |  back-substitution (by using the subroutine SOLVE_BLOCK).
     |
     |  In input, BABDCR_SOLV requires the coefficient matrix factorized by
     |  the subroutine BABDCR_FACT (see details of this subroutine) and the
     |  right hand side which is stored in the block vector VECT_B in the
     |  following form
     |
     |  VECT_B  = [f(0), f(1), ...., f(NBLOKS)].
     |
     |  The first block element f(0) corresponds to the first row of the
     |  coefficient matrix containing Ba and Bb.
     |
     |  On exit, BABDCR_SOLV gives the solution of the babd linear system.
     |  The solution is stored in the block vector VECT_B.
     |
    \*/
    // some constanst
    integer const nxn = n*n ;
    integer const nx2 = 2*n ;

    // initialize indices
    valuePointer F    = F_blk ;
    integer *    perm = perm_blk ;
    
    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! reduction phase !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!
    */
    integer jump = 1 ;
    for ( ; jump < nblock ; jump *=2 ) {

      for ( integer k = 1 ; k*jump < nblock ; k += 2 ) {

        integer      jC = k * jump ;
        integer      jL = jC - jump ;
        valuePointer RS = SR_blk + (2*jC-1)*nxn ;
        valuePointer Va = y + jL*n ;
        valuePointer Vb = y + jC*n ;

        /*\
         |   P * / Va \ = / Vp \
         |       \ Vb /   \ Vq /
        \*/
        copy( n, Va, 1, TMP_blk,   1 ) ; 
        copy( n, Vb, 1, TMP_blk+n, 1 ) ;

        for ( integer i = 0 ; i < n ; ++i ) {
          Va[i] = TMP_blk[perm[i+n]] ;
          Vb[i] = TMP_blk[perm[i]] ;
        }
        /*\
         |  / -G  I \ / Vp \ = / Vq-G*Vp \
         |  \  I  0 / \ Vq /   \ Vp      /
        \*/
        gemv( Transposition::NO_TRANSPOSE,
              n, n, -1.0, RS + n, nx2, Vb, 1, 1.0, Va, 1 ) ;

        perm += nx2 ;
        F    += nxn ;
      }
    }

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! 2 by 2 block linear system solution !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    integer nm = n+m ;

    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    valuePointer ye = y + (nblock-1) * n ;
    
    // Apply row interchanges to the right hand sides.
    swap( n, y, 1, ye, 1 ) ;
    swaps( 1, ye, nm, 0, nm-1, ipiv_blk, 1 ) ;
    trsv( alglin::ULselect::LOWER, alglin::Transposition::NO_TRANSPOSE, alglin::DiagonalType::UNIT,     nm, D_blk, nm, ye, 1 ) ;
    trsv( alglin::ULselect::UPPER, alglin::Transposition::NO_TRANSPOSE, alglin::DiagonalType::NON_UNIT, nm, D_blk, nm, ye, 1 ) ;
    swap( n, y, 1, ye, 1 ) ;

    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! back-substitution phase !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    for ( jump /= 2 ; jump > 0 ; jump /= 2 ) {

      integer k = (nblock-1)/jump ;
      if ( (k & 0x01) == 0 ) --k ; 
    
      for ( ; k > 0 ; k -= 2 ) {
        integer      jC = k*jump ;
        integer      jL = jC-jump ;
        integer      jR = min(jC + jump,nblock) ;
        valuePointer Xm = y + jL*n ;
        valuePointer Xs = y + jR*n ;
        valuePointer V  = y + jC*n ;
        valuePointer RS = SR_blk + (2*jC-1)*nxn ;

        perm -= nx2 ;
        F    -= nxn ;

        /*\
         |  computes the solution Xn (of length n) of the linear system
         |
         |          L*U Xn = Vp - Sp*Xm - Rp*Xs
         |
         |  obtained from the subroutines reduceBlock and reduceRHS applied
         |  to the system (see the subroutine reduceBlock for further details)
         |
         |   ( Sa  Rb     )  ( Xm )   ( Va ).
         |   (     Sb  Rc )  ( Xn ) = ( Vb )
         |                   ( Xs )
         |
         |   (  I    ) * P * ( Sa  Rb     ) ( Xm )   ( Sp  L*U  Rp  ) ( Xm )   ( Vp )
         |   ( -G  I )       (     Sb  Rc ) ( Xn ) = ( Sa'  0   Rc' ) ( Xn ) = ( V' )
         |                                  ( Xs )                    ( Xs )
        \*/
     
        // compute V = Vp - Sp*Xm - Rp*Xs
        for ( integer i = 0 ; i < n ; ++i )
          V[i] -= dot( n, F + i*n, 1, perm[i] < n ? Xm : Xs, 1 ) ;

        // solve the system L*U Xn = V
        trsv( ULselect::LOWER, Transposition::NO_TRANSPOSE, DiagonalType::UNIT,     n, RS, nx2, V, 1 ) ;
        trsv( ULselect::UPPER, Transposition::NO_TRANSPOSE, DiagonalType::NON_UNIT, n, RS, nx2, V, 1 ) ;
      }
    }
  }

  template class AmodioLU<double> ;
  template class AmodioLU<float> ;

}
