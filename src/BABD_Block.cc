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

#include "BABD_Block.hh"
#include "Alglin.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace alglin {

  using namespace std ;

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  BBlockLU<t_Value>::factorize() {

    integer const & nblock = this->nblock ;
    integer const & n      = this->n ;
    integer const & q      = this->q ;
    integer const & nx2    = this->nx2 ;
    integer const & nxn    = this->nxn ;
    integer const & nxnx2  = this->nxnx2 ;

    // fill matrix
    integer nm = n+m ;
    for ( integer k = 0 ; k < nblock ; ++k ) {
      valueConstPointer Ad = this->AdAu_blk + k*nxnx2 ;
      valueConstPointer Au = Ad + nxn ;
      gecopy( n, n, Ad, n, AdH_blk + k*nm*n, nm ) ;
      gecopy( n, n, Au, n, Au_blk  + k*nxn,  n  ) ;
    }

    gecopy( m, n, this->H0Nq,       m, AdH_blk + n,  nm ) ;
    gecopy( m, n, this->H0Nq+n*m,   m, DD_blk,       m  ) ;
    gecopy( m, q, this->H0Nq+nx2*m, m, DD_blk + m*n, m  ) ;

    integer rowFF = (nblock-1)*n ;
    integer INFO ;

    integer *  ipivk = ipiv_blk ;
    valuePointer AdH = AdH_blk  ;
    valuePointer Au  = Au_blk   ;
    valuePointer FF  = FF_blk   ;

    for ( integer k = 0 ; k < nblock-1 ; ++k, ipivk += n, AdH += nm*n, Au += n*n, FF += n ) {

      INFO = getrf( nm, n, AdH, nm, ipivk ) ; // LU factorization
      ALGLIN_ASSERT( INFO==0, "BlockLU::factorize(), matrix singolar" ) ;

      valuePointer H  = AdH + n ;
      valuePointer CC = AdH + nm*n + n ;

      for ( integer i = 0 ; i < n ; ++i ) {
        integer ip = ipivk[i]-1 ;
        if ( ip != i ) { // detected row exchange
          if ( ip < n ) { // exchange row
            swap( n, Au + i, n,     Au + ip, n     ) ;
            swap( m, FF + i, rowFF, FF + ip, rowFF ) ; // last column block
          } else {
            swap( n, Au + i, n,     CC     + (ip-n), nm ) ;
            swap( m, FF + i, rowFF, DD_blk + (ip-n), m  ) ; // last column block
          }
        }
      }

      //
      //  +---------+---------+ ........ +--------+
      //  | \    U  |         |          |        |
      //  |    \    |L^(-1)Au |          |L^(-1)FF|
      //  |   L   \ |         |          |        |
      //  |=========|=========|          +========+
      //  |         |         |          |        |
      //  |    H    |   CC*   |          |   DD*  |
      //  |         |         |          |        |
      //  +---------+---------+ ........ +--------+
      //
      //  CC* = CC - H (LU)^(-1) Au
      //  DD* = DD - H (LU)^(-1) FF
      //
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, n, 1, AdH, nm, Au, n ) ;
      trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, m, 1, AdH, nm, FF, rowFF ) ;

      gemm( NO_TRANSPOSE, NO_TRANSPOSE, m, n, n, -1, H, nm, Au, n, 1, CC, nm ) ;
      gemm( NO_TRANSPOSE, NO_TRANSPOSE, m, m, n, -1, H, nm, FF, rowFF, 1, DD_blk, m ) ;
    }

    // factorize last two block
    INFO = getrf( nm, n, AdH, nm, ipivk ) ; // LU factorization
    ALGLIN_ASSERT( INFO==0, "BlockLU::factorize(), matrix singolar" ) ;

    for ( integer i = 0 ; i < n ; ++i ) {
      integer ip = ipivk[i]-1 ;
      if ( ip != i ) { // detected row exchange
        if ( ip < n ) { // exchange row
          swap( m, Au + i, n, Au + ip, n ) ; // last column block
        } else {
          swap( m, Au + i, n, DD_blk + (ip-n), m ) ; // last column block
        }
      }
    }

    //
    // +---------+---------+
    // | \    U  |         |
    // |    \    |L^(-1)Au |
    // |   L   \ |         |
    // |=========|=========|
    // |         |         |
    // |    H    |   CC*   |
    // |         |         |
    // +---------+---------+
    //
    //  CC* = CC - H (LU)^(-1) Au
    //
    valuePointer H = AdH + n ;
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT, n, m, 1, AdH, nm, Au, n ) ;
    gemm( NO_TRANSPOSE, NO_TRANSPOSE, m, m, n, -1, H, nm, Au, n, 1, DD_blk, m ) ;
    
    // factorize last block
    ipivk += n ;
    INFO = getrf( m, m, DD_blk, m, ipivk ) ; // LU factorization
    ALGLIN_ASSERT( INFO==0, "BlockLU::factorize(), singular matrix" ) ;
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/
  template <typename t_Value>
  void
  BBlockLU<t_Value>::solve( valuePointer y ) const {

    integer const & n      = this->n ;
    integer const & nblock = this->nblock ;

    // solve L
    integer nm    = n+m ;
    integer rowFF = (nblock-1) * n ;
    valuePointer ye = y + nblock * n ;

    for ( integer k = 0 ; k < nblock ; ++k ) {
      integer const * ipivk = ipiv_blk + k * n ;
      valueConstPointer AdH = AdH_blk  + k * (nm*n) ;
      valuePointer      yk  = y        + k * n ;

      // apply permutation
      for ( integer i = 0 ; i < n ; ++i ) {
        integer ip = ipivk[i]-1 ;
        if ( ip != i ) { // detected row exchange
          if ( ip < n ) std::swap( yk[i], yk[ip] ) ;
          else          std::swap( yk[i], ye[ip-n] ) ;
        }
      }
      trsv( LOWER, NO_TRANSPOSE, UNIT, n, AdH, nm, yk, 1 ) ;
      gemv( NO_TRANSPOSE, m, n, -1, AdH + n, nm, yk, 1, 1, ye, 1 ) ;
    }

    integer const * ipive = ipiv_blk + nblock * n ;
    integer            ok = getrs( NO_TRANSPOSE, m, 1, DD_blk, m, ipive, ye, m ) ;

    ALGLIN_ASSERT( ok == 0, "BlockLU::solve(...) failed" ) ;

    if ( rowFF > 0 ) gemv( NO_TRANSPOSE, rowFF, m, -1, FF_blk, rowFF, ye, 1, 1, y, 1 ) ;
    if (     m > n ) gemv( NO_TRANSPOSE, n, m-n, -1, Au_blk + nblock*n*n, n, ye+n, 1, 1, ye-n, 1 ) ;

    integer k = nblock ;
    do {
      --k ;
      valueConstPointer AdH = AdH_blk + k*nm*n ;
      valueConstPointer Au  = Au_blk  + k*n*n  ;
      valuePointer      yk  = y       + k*n    ;

      gemv( NO_TRANSPOSE, n, n, -1, Au, n, yk + n, 1, 1, yk, 1 ) ;
      trsv( UPPER, NO_TRANSPOSE, NON_UNIT, n, AdH, nm, yk, 1 ) ;

    } while ( k > 0 ) ;
  }

  template class BBlockLU<double> ;
  template class BBlockLU<float> ;

}
