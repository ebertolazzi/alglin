/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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

#include "Alglin.hh"
#include "Alglin_Eigen.hh"

namespace alglin {

  using std::swap;
  using std::max;
  using std::min;

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

    integer const & nblock{m_number_of_blocks};
    integer const & n{m_block_size};
    integer const & q{m_extra_bc};

    // fill matrix
    integer m{n+q};
    integer nm{n+m};
    for ( integer k{0}; k < nblock; ++k ) {
      real_type const * Ad{m_DE_blk + (2*k)*n_x_n};
      real_type const * Au{Ad + n_x_n};
      GEcopy( n, n, Ad, n, m_AdH_blk + k*nm*n,  nm );
      GEcopy( n, n, Au, n, m_Au_blk  + k*n_x_n, n  );
    }

    GEcopy( m, n, m_H0Nq,       m, m_AdH_blk + n,  nm );
    GEcopy( m, n, m_H0Nq+n*m,   m, m_DD_blk,       m  );
    GEcopy( m, q, m_H0Nq+2*n*m, m, m_DD_blk + m*n, m  );

    integer rowFF{(nblock-1)*n};
    integer INFO;

    integer   * ipivk {m_ipiv_blk};
    real_type * AdH   {m_AdH_blk};
    real_type * Au    {m_Au_blk};
    real_type * FF    {m_FF_blk};

    for ( integer k{0};
          k < nblock-1;
          ++k, ipivk += n, AdH += nm*n, Au += n*n, FF += n ) {

      INFO = getrf( nm, n, AdH, nm, ipivk ); // LU factorization
      UTILS_ASSERT0( INFO==0, "BlockLU::factorize(), matrix singular\n" );

      real_type * H  { AdH + n };
      real_type * CC { AdH + nm*n + n };

      for ( integer i{0}; i < n; ++i ) {
        integer ip{ipivk[i]-1};
        if ( ip != i ) { // detected row exchange
          if ( ip < n ) { // exchange row
            alglin::swap( n, Au + i, n,     Au + ip, n     );
            alglin::swap( m, FF + i, rowFF, FF + ip, rowFF ); // last column block
          } else {
            alglin::swap( n, Au + i, n,     CC       + (ip-n), nm );
            alglin::swap( m, FF + i, rowFF, m_DD_blk + (ip-n), m  ); // last column block
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
      trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, n, 1, AdH, nm, Au, n
      );
      trsm(
        SideMultiply::LEFT,
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, m, 1, AdH, nm, FF, rowFF
      );

      gemm(
        Transposition::NO,
        Transposition::NO,
        m, n, n, -1, H, nm, Au, n, 1, CC, nm
      );
      gemm(
        Transposition::NO,
        Transposition::NO,
        m, m, n, -1, H, nm, FF, rowFF, 1, m_DD_blk, m
      );
    }

    // factorize last two block
    INFO = getrf( nm, n, AdH, nm, ipivk ); // LU factorization
    UTILS_ASSERT0( INFO==0, "BlockLU::factorize(), matrix singular\n" );

    for ( integer i{0}; i < n; ++i ) {
      integer ip{ipivk[i]-1};
      if ( ip != i ) { // detected row exchange
        if ( ip < n ) { // exchange row
          alglin::swap( m, Au + i, n, Au + ip, n ); // last column block
        } else {
          alglin::swap( m, Au + i, n, m_DD_blk + (ip-n), m ); // last column block
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
    real_type * H{AdH + n};
    trsm(
      SideMultiply::LEFT,
      ULselect::LOWER,
      Transposition::NO,
      DiagonalType::UNIT,
      n, m, 1, AdH, nm, Au, n
    );
    gemm(
      Transposition::NO,
      Transposition::NO,
      m, m, n, -1, H, nm, Au, n, 1, m_DD_blk, m
    );

    // factorize last block
    ipivk += n;
    INFO = getrf( m, m, m_DD_blk, m, ipivk ); // LU factorization
    UTILS_ASSERT0( INFO==0, "BlockLU::factorize(), singular matrix\n" );
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
  BBlockLU<t_Value>::solve( real_type y[] ) const {

    integer const & nblock { m_number_of_blocks };
    integer const & n      { m_block_size };
    integer const & q      { m_extra_bc };

    // solve L
    integer     m     {n+q};
    integer     nm    {n+m};
    integer     rowFF {(nblock-1) * n};
    real_type * ye    {y + nblock * n};

    for ( integer k{0}; k < nblock; ++k ) {
      integer   const * ipivk {m_ipiv_blk + k * n};
      real_type const * AdH   {m_AdH_blk  + k * (nm*n)};
      real_type       * yk    {y          + k * n};

      // apply permutation
      for ( integer i{0}; i < n; ++i ) {
        integer ip{ipivk[i]-1};
        if ( ip != i ) { // detected row exchange
          if ( ip < n ) std::swap( yk[i], yk[ip] );
          else          std::swap( yk[i], ye[ip-n] );
        }
      }
      trsv(
        ULselect::LOWER,
        Transposition::NO,
        DiagonalType::UNIT,
        n, AdH, nm, yk, 1
      );
      gemv(
        Transposition::NO,
        m, n, -1, AdH + n, nm, yk, 1, 1, ye, 1
      );
    }

    integer const * ipive { m_ipiv_blk + nblock * n };
    integer const   ok    { getrs( Transposition::NO, m, 1, m_DD_blk, m, ipive, ye, m ) };

    UTILS_ASSERT0( ok == 0, "BlockLU::solve(...) failed\n" );

    if ( rowFF > 0 ) gemv( Transposition::NO, rowFF, m, -1, m_FF_blk, rowFF, ye, 1, 1, y, 1 );
    if (     m > n ) gemv( Transposition::NO, n, m-n, -1, m_Au_blk + nblock*n*n, n, ye+n, 1, 1, ye-n, 1 );

    integer k = nblock;
    do {
      --k;
      real_type const * AdH { m_AdH_blk + k*nm*n };
      real_type const * Au  { m_Au_blk  + k*n*n  };
      real_type       * yk  { y         + k*n    };

      gemv(
        Transposition::NO,
        n, n, -1, Au, n, yk + n, 1, 1, yk, 1
      );
      trsv(
        ULselect::UPPER,
        Transposition::NO,
        DiagonalType::NON_UNIT,
        n, AdH, nm, yk, 1
      );

    } while ( k > 0 );
  }

  template class BBlockLU<double>;
  template class BBlockLU<float>;

}
