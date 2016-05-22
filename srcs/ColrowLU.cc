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
/// file: ColrowLU.cc
///

#include "ColrowLU.hh"
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>

/*
//  Blocco algoritmo di Diaz
//
//   +-----+-----+
//   |  A  |  B  |
//   +-----+-----+-----------+
//   |     |     |           |
//   |  C  |  D  |     E     |
//   |     |     |           |
//   |     |     |           |
//   +-----+-----+-----------+-----------+
//               |                       |
//               |                       |
//               |                       |
//               |                       |
//               +-----------------------+
//
*/

namespace alglin {

  template <> float  const ColrowLU<float>::epsi  = 100*std::numeric_limits<float>::epsilon() ;
  template <> double const ColrowLU<double>::epsi = 1000*std::numeric_limits<double>::epsilon() ;

  template <typename t_Value>
  ColrowLU<t_Value>::ColrowLU( bool _use_arceco )
  : baseValue("ColrowLU reals")
  , baseIndex("ColrowLU ints")
  , last_block(ColrowLU_QR0)
  //, last_block(ColrowLU_QR1)
  //, last_block(ColrowLU_QR2)
  //, last_block(ColrowLU_LU0)
  //, last_block(ColrowLU_LU1)
  , use_arceco(_use_arceco)
  {
  }
  
  template <typename t_Value>
  ColrowLU<t_Value>::~ColrowLU() {
    baseValue.free() ;
    baseIndex.free() ;
  }
  
  // ---------------------------------------------------------------------------

  //! compute y += alpha*A*x
  template <typename t_Value>
  void
  ColrowLU<t_Value>::mv( integer           _row0,
                         integer           _col0,
                         valueConstPointer _block0,
                        
                         integer           _numBlock,
                         integer           _dimBlock,
                         valueConstPointer _blocks,
                        
                         integer           _rowN,
                         integer           _colN,
                         valueConstPointer _blockN,
                        
                         valueType         alpha,
                         valueConstPointer x,
                         integer           incx,
                        
                         valueType         beta,
                         valuePointer      y,
                         integer           incy) {

    // first block y = alpha * _block0 * x + beta * y
    gemv( Transposition::NO_TRANSPOSE, _row0, _col0,
          alpha, _block0, _row0,
          x, incx,
          beta, y, incy ) ;

    // internal blocks block
    valueConstPointer xx   = x+(_col0-_dimBlock)*incx ;
    valuePointer      yy   = y+_row0*incy ;
    valueConstPointer blks = _blocks ;
    for ( integer i = 0 ; i < _numBlock ; ++i ) {
      gemv( Transposition::NO_TRANSPOSE, _dimBlock, 2*_dimBlock,
            alpha, blks, _dimBlock,
            xx, incx,
            beta, yy, incy ) ;
      xx   += _dimBlock*incx ;
      yy   += _dimBlock*incy ;
      blks += 2*_dimBlock*_dimBlock ;
    }

    // last block
    gemv( Transposition::NO_TRANSPOSE, _rowN, _colN,
          alpha, _blockN, _rowN,
          xx, incx,
          beta, yy, incy ) ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::residue(  integer           _row0,
                               integer           _col0,
                               valueConstPointer _block0,
                               integer           _numBlock,
                               integer           _dimBlock,
                               valueConstPointer _blocks,
                               integer           _rowN,
                               integer           _colN,
                               valueConstPointer _blockN,
                               valueConstPointer b,
                               integer           incb,
                               valueConstPointer x,
                               integer           incx,
                               valuePointer      res,
                               integer           incr ) {

    copy( _numBlock*_dimBlock+_row0+_rowN, b, incb, res, incr ) ;
    mv( _row0, _col0, _block0,
        _numBlock, _dimBlock, _blocks,
        _rowN, _colN, _blockN,
        -1.0, x, incx, 1.0, res, incr ) ;
  }

  /*
  // ---------------------------------------------------------------------------
  //             col0
  //       +---+----------+
  // col00 |              |
  //       +   + - - - -  |  row0
  // row00 |   :          |
  //       |   :          |
  //       +---+----------+----------+
  //     col00 |                     |
  //           |                     | dimBlock
  //           |                     |
  //           +----------+----------+
  //                 2*dimBlock
  //
  //     +----------+----------+
  //     |                     |
  //     |                     |
  //     |                     | colNN
  //     +----------+----------+-----+
  //                |                |
  //                |                | rowN
  //                |                |
  //                |                |
  //                +----------------+
  //                      colN
  //
  */
  template <typename t_Value>
  void
  ColrowLU<t_Value>::setup() {

    hSizeBlock = dimBlock*dimBlock ;
    sizeBlock  = 2*hSizeBlock ;

    col00 = col0 - dimBlock ;
    colNN = colN - dimBlock ;
    row00 = row0 - col00 ;
    rowNN = rowN - colNN ;
    neq   = numBlock * dimBlock + row0 + rowN ;

    integer neq1 = (numBlock+1) * dimBlock + col00 + colNN ;

    dimBlock_m_row00 = dimBlock - row00 ;

    ALGLIN_ASSERT( row0 > 0 && rowN > 0 && row00 >= 0 && col00 >= 0 &&
                   col0 >= dimBlock && colN >= dimBlock && dimBlock > 0,
                   "Bad parameter(s):" <<
                   "\nrow0     = " << row0 <<
                   "\ncol0     = " << col0 <<
                   "\ncol00    = " << col00 <<
                   "\nrowN     = " << rowN <<
                   "\ncolN     = " << colN <<
                   "\ncolNN    = " << colNN <<
                   "\ndimBlock = " << dimBlock ) ;

    ALGLIN_ASSERT( neq == neq1,
                   "Bad parameter(s):" <<
                   "\nrow0     = " << row0 <<
                   "\ncol0     = " << col0 <<
                   "\ncol00    = " << col00 <<
                   "\nrowN     = " << rowN <<
                   "\ncolN     = " << colN <<
                   "\ncolNN    = " << colNN <<
                   "\ndimBlock = " << dimBlock <<
                   "\nneq      = " << neq <<
                   "\nneq1     = " << neq1 <<
                   "\nneq and neq1 = (numBlock+1) * dimBlock + col00 + colNN must be equal" ) ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize( integer           _row0,
                                integer           _col0,
                                valueConstPointer _block0,
                                integer           _numBlock,
                                integer           _dimBlock,
                                valueConstPointer _blocks,
                                integer           _rowN,
                                integer           _colN,
                                valueConstPointer _blockN ) {
    numBlock = _numBlock ;
    dimBlock = _dimBlock ;
    row0     = _row0 ;
    col0     = _col0 ;
    rowN     = _rowN ;
    colN     = _colN ;
    
    setup() ;
    
    if ( use_arceco ) {

      baseIndex.allocate( 3*(numBlock+2) + neq ) ;
      baseValue.allocate( sizeBlock*numBlock + row0*col0 + rowN*colN ) ;

      matrixStructure = baseIndex( 3*(numBlock+2) ) ;
      array           = baseValue( sizeBlock*numBlock + row0*col0 + rowN*colN ) ;
      pivot           = baseIndex( neq ) ;

      integer * mS = matrixStructure ;
      valuePointer      a = array ;
      valueConstPointer b = _blocks ;
      *mS++ = row0 ;
      *mS++ = col0 ;
      *mS++ = dimBlock ;
      std::copy( _block0, _block0 + row0*col0, a ) ; a += row0*col0 ;
      for ( integer nb = 0 ; nb < numBlock ; ++nb ) {
        *mS++ = dimBlock ;
        *mS++ = 2*dimBlock ;
        *mS++ = dimBlock ;
        std::copy( b, b+sizeBlock, a ) ;
        a += sizeBlock ;
        b += sizeBlock ;
      }
      *mS++ = rowN ;
      *mS++ = colN ;
      *mS++ = 0 ;
      std::copy( _blockN, _blockN+rowN*colN, a ) ;

      LU_arceco.loadByRef( numBlock+2, matrixStructure, array, pivot ) ;
      LU_arceco.factorize() ;

    } else {

      // allocate
      baseIndex.allocate( 2*(numBlock*dimBlock+row0) ) ;
      baseValue.allocate( sizeBlock*numBlock + row0*col0 + rowN*colN ) ;
      block0      = baseValue( row0*col0 ) ;
      blocks      = baseValue( sizeBlock * numBlock ) ;
      blockN      = baseValue( rowN*colN ) ;
      swapRC_blks = baseIndex( 2*(numBlock*dimBlock+row0) ) ;

      // copy block
      copy( row0*col0,          _block0, 1, block0, 1 ) ;
      copy( sizeBlock*numBlock, _blocks, 1, blocks, 1 ) ;
      copy( rowN*colN,          _blockN, 1, blockN, 1 ) ;

      factorize() ;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_inplace(  integer      _row0,
                                         integer      _col0,
                                         valuePointer _block0,
                                         integer      _numBlock,
                                         integer      _dimBlock,
                                         valuePointer _blocks,
                                         integer      _rowN,
                                         integer      _colN,
                                         valuePointer _blockN ) {
    numBlock = _numBlock ;
    dimBlock = _dimBlock ;
    row0     = _row0 ;
    col0     = _col0 ;
    rowN     = _rowN ;
    colN     = _colN ;

    setup() ;

    baseIndex.allocate( 2*(numBlock*dimBlock+row0) ) ;
    baseValue.allocate( sizeBlock*numBlock + row0*col0 + rowN*colN ) ;
    swapRC_blks = baseIndex( 2*(numBlock*dimBlock+row0) ) ;

    factorize() ;
  }

  /*
  //
  // +---+-----+--------+
  // |   |     |        |
  // | L |  A  |    R   |
  // |   |     |        |
  // +---+-----+--------+
  //
  //
  // +---+-----+--------+
  // |   | \ U |        |
  // | L |   \ |    R   |
  // |   |  L  |        |
  // +---+-----+--------+
  //
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::LU_left_right( integer nrA,
                                    integer ncA,
                                    integer ncL,
                                    integer ncR,
                                    t_Value * A, integer ldA,
                                    integer swapR[],
                                    integer swapC[] ) {

    t_Value * Ak  = A ;
    t_Value * Akk = A ;
    t_Value * L   = A - ncL*ldA ;

    integer ncc = ncA+ncL+ncR ;

    // LU full pivoting
    integer rr = nrA ;
    integer cc = ncA+ncR ;
    for ( integer k = 0 ; k < ncA ; ++k ) {
      // cerco pivot
      integer ipiv = k, jpiv = k ;
      t_Value amax = std::abs(Akk[0]) ;
      for ( integer i = k ; i < nrA ; ++i ) {
        for ( integer j = k ; j < ncA ; ++j ) {
          t_Value amax1 = std::abs(A[i+j*ldA]) ;
          if ( amax1 > amax ) { ipiv = i ; jpiv = j ; amax = amax1 ; }
        }
      }
      // controllo pivot non zero
      ALGLIN_ASSERT( amax > epsi,
                     "ColrowLU::LU_left_right, found pivot: " << Akk[0] <<
                     " at block nblk = " << nblk ) ;
      // memorizzo scambio
      swapR[k] = ipiv ;
      swapC[k] = jpiv ;
      // applico permutazione
      if ( ipiv > k ) swap( ncc, L + k, ldA, L + ipiv, ldA ) ; // scambio righe
      if ( jpiv > k ) swap( nrA, A + k*ldA, 1, A + jpiv*ldA, 1 ) ; // scambio colonne
      // memorizzo L^(-1)
      t_Value * X = Akk + 1 ;
      t_Value * Y = Akk + ldA ;
      // aggiorno k
      --rr ;
      --cc ;
      rscal( rr, Akk[0], X, 1 ) ;
      // applico update rank 1 ==> A := A - x*y'
      Ak  += ldA ;
      Akk += ldA+1 ;
      ger( rr, cc, -1, X, 1, Y, ldA, Akk, ldA ) ;
    }
  }

  /*
  //                     columns permuted
  //   +-----------+      +-----------+
  //   |     T     |      |     T     |
  //   +-----------+      +-----------+
  //   |     A     |      | L \   U   |
  //   +===========+  ==> +====+======+
  //   |           |      |    |      |
  //   |     B     |      | L  |  U   |
  //   |           |      |    |      |
  //   +-----------+      +----+------+
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::LU_top_bottom( integer nrT,
                                    integer nrA,
                                    integer ncA,
                                    t_Value * A, integer ldA,
                                    integer nrB,
                                    t_Value * B, integer ldB,
                                    integer swapR[],
                                    integer swapC[] ) {
    t_Value * Ak  = A ;
    t_Value * Akk = A ;
    t_Value * Bk  = B ;
    t_Value * T   = A-nrT ;

    // LU full pivoting
    integer rr   = nrA ;
    integer cc   = ncA ;
    integer nrAT = nrA + nrT ;
    for ( integer k = 0 ; k < nrA ; ++k ) {
      // cerco pivot
      integer ipiv = k, jpiv = k ;
      t_Value amax = std::abs(Akk[0]) ;
      for ( integer i = k ; i < nrA ; ++i ) {
        for ( integer j = k ; j < ncA ; ++j ) {
          t_Value amax1 = std::abs(A[i+j*ldA]) ;
          if ( amax1 > amax ) { ipiv = i ; jpiv = j ; amax = amax1 ; }
        }
      }
      // controllo pivot non zero
      ALGLIN_ASSERT( amax > epsi,
                     "ColrowLU::LU_top_bottom, found pivot: " << Akk[0] <<
                     " at block nblk = " << nblk ) ;
      // memorizzo scambio
      swapR[k] = ipiv ;
      swapC[k] = jpiv ;
      // applico permutazione
      if ( ipiv > k ) swap( ncA, A + k, ldA, A + ipiv, ldA ) ; // scambio righe
      if ( jpiv > k ) {
        swap( nrAT, T + k*ldA, 1, T + jpiv*ldA, 1 ) ; // scambio colonne
        swap( nrB,  B + k*ldB, 1, B + jpiv*ldB, 1 ) ; // scambio colonne
      }
      // memorizzo L^(-1)
      t_Value * X = Akk + 1 ;
      t_Value * Y = Akk + ldA ;
      t_Value * Z = Bk ;
      // aggiorno k
      --rr ;
      --cc ;
      t_Value rval = 1/Akk[0] ;
      scal( rr,  rval, X,  1 ) ;
      scal( nrB, rval, Bk, 1 ) ;
      // applico update rank 1 ==> A := A - x*y'
      Ak  += ldA ;
      Akk += ldA+1 ;
      Bk  += ldB ;
      ger( rr,  cc, -1, X, 1, Y, ldA, Akk, ldA ) ;
      ger( nrB, cc, -1, Z, 1, Y, ldA, Bk,  ldB ) ;
    }
  }

  /*
  //  col00
  //  +----+----------+
  //  | A0 |    B0    |  \
  //  +----+----------+  ldA = row0
  //  | A1 |    B1    |  /
  //  +----+----------+
  //
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_first_block() {
    if ( col00 > 0 ) {
      integer * swapC = swapRC_blks ;
      integer * swapR = swapRC_blks+col00 ;
      LU_left_right( row0, col00, 0, col0-col00,
                     block0, row0, swapR, swapC ) ;
    }
  }

  /*
  //    dimA   dimBlock
  //  +--///--+----------+
  //  |  ---  | B0       | ldA-row00
  //  +       +----------+
  //  |  A1   | B1       | row00
  //  +--///--+----+-----+----------+
  //          |    |     |          |
  //          | C  | D   |    E     |
  //          |    |     |          |
  //          +----+-----+----------+
  //          row00
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_block() {

    valuePointer C = blocks + nblk * sizeBlock ;
    valuePointer D = C + row00 * dimBlock ;

    integer      ldA, dimA, rowT ;
    valuePointer B0, B1, A1 ;
    if ( nblk == 0 ) {
      ldA  = row0 ;
      dimA = col00 ;
      rowT = col00 ;
      A1   = block0 + col00 ;
      B0   = block0 + row0*col00 ;
      B1   = B0     + col00 ;
    } else {
      ldA  = dimBlock ;
      dimA = dimBlock ;
      rowT = dimBlock_m_row00 ;
      B0   = C - hSizeBlock ;
      B1   = B0 + dimBlock_m_row00 ;
      A1   = B1 - hSizeBlock ;
    }

    integer * swapC = swapRC_blks + 2*(col00+nblk*dimBlock) ;
    integer * swapR = swapC + row00 ;

    LU_top_bottom( rowT, row00, dimBlock, B1, ldA, dimBlock, C, dimBlock, swapR, swapC ) ;

    // propago permutazioni al blocco a sinistra
    for ( integer i = 0 ; i < row00 ; ++i ) {
      integer i1 = swapR[i] ;
      if ( i1 > i ) swap( dimA, A1 + i, ldA, A1 + i1, ldA ) ;
    }

    swapC = swapR + row00 ;
    swapR = swapC + dimBlock_m_row00 ;
    //                 NR          NC            L       R
    LU_left_right( dimBlock, dimBlock_m_row00, row00, dimBlock,
                   D, dimBlock, swapR, swapC ) ;

    // propago permutazioni ai vari blocchi
    for ( integer i = 0 ; i < dimBlock_m_row00 ; ++i ) {
      integer i1 = swapC[i] ;
      if ( i1 > i ) swap( ldA, B0 + (i+row00) * ldA, 1, B0 + (i1+row00) * ldA, 1 ) ;
    }

  }

  // ---------------------------------------------------------------------------

  /*
  //    dimA   dimBlock
  //  +--///--+----------+
  //  |  ---  | B0       | row11 = ldA-row00 o col00
  //  +       +----------+
  //  |  A1   | B1       |
  //  +--///--+----+-----+-----+
  //          | C0 | D0  |  E0 |
  //          |    |     |     | rowN
  //          | C1 | D1  |  E1 |
  //          +----+-----+-----+
  //          row00       colNN
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_last_block() {

    valuePointer C0 = blockN ;
    valuePointer D0 = C0 + row00 * rowN ;
    valuePointer B0, B1, A1 ;
    integer ldA, dimA, rowT ;

    if ( numBlock == 0 ) {
      A1   = block0 + col00 ;
      B0   = block0 + row0 * col00 ;
      B1   = B0 + col00 ;
      ldA  = row0 ;
      dimA = col00 ;
      rowT = col00 ;
    } else {
      B0   = blocks + numBlock * sizeBlock - hSizeBlock ;
      B1   = B0 + dimBlock_m_row00 ;
      A1   = B1 - hSizeBlock ;
      ldA  = dimBlock ;
      dimA = dimBlock ;
      rowT = dimBlock_m_row00 ;
    }

    integer * swapC = swapRC_blks + 2*( col00 + numBlock*dimBlock ) ;
    integer * swapR = swapC + row00 ;

    LU_top_bottom( rowT, row00, dimBlock,
                   B1, ldA,
                   rowN, blockN, rowN,
                   swapR, swapC ) ;

    // propago permutazioni ai vari blocchi
    for ( integer i = 0 ; i < row00 ; ++i ) {
      integer i1 = swapR[i] ;
      if ( i1 > i ) swap( dimA, A1 + i, ldA, A1 + i1, ldA ) ;
    }

    // fattorizzazione ultimo blocco
    switch ( last_block ) {
      case ColrowLU_QR0: la_solve0.compute(Eigen::Map<mat>(D0,rowN,rowN)) ; break ;
      case ColrowLU_QR1: la_solve1.compute(Eigen::Map<mat>(D0,rowN,rowN)) ; break ;
      case ColrowLU_QR2: la_solve2.compute(Eigen::Map<mat>(D0,rowN,rowN)) ; break ;
      case ColrowLU_LU0: la_solve3.compute(Eigen::Map<mat>(D0,rowN,rowN)) ; break ;
      case ColrowLU_LU1: la_solve4.compute(Eigen::Map<mat>(D0,rowN,rowN)) ; break ;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize() {
    // primo blocco
    nblk = -1 ;
    factorize_first_block() ;

    // blocchi intermedi (n-1)
    for ( nblk = 0 ; nblk < numBlock ; ++nblk ) factorize_block() ;

    // copio in ultimo blocco
    factorize_last_block() ;
  }

  /*
  //
  //  +----+---+---+---+
  //  |  M | L |   :   |
  //  |  M | L | L :   |
  //  +----+---+---+---+
  //  row00
  */
  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve_block_L( valuePointer io ) const {

    valuePointer io1 = io - row00 ;
    valuePointer M   = blocks + nblk * sizeBlock ;
    valuePointer L   = M + row00 * dimBlock ;

    // io -= M*io1
    gemv( Transposition::NO_TRANSPOSE,
          dimBlock, row00,
          -1, M, dimBlock,
          io1, 1,
          1, io, 1 ) ;

    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          dimBlock,
          L, dimBlock, io, 1 ) ;

  }

  /*
  //  +---------+---+-------+
  //  |       U |   |       |
  //  |         | U |   M   |
  //  +---------+---+-------+
  //  row00
  */
  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve_block_U( valuePointer io ) const {

    valuePointer io1 = io + dimBlock ;
    valuePointer U   = blocks + nblk * sizeBlock + row00 * dimBlock ;
    valuePointer M   = U + hSizeBlock ;

    gemv( Transposition::NO_TRANSPOSE,
          dimBlock, dimBlock_m_row00,
          -1, M, dimBlock,
          io1, 1,
          1, io, 1 ) ;

    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          dimBlock,
          U, dimBlock, io, 1 ) ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve_last_block( valuePointer io ) const {
    // soluzione ultimo blocco
    integer ncol = colN-rowN ;
    gemv( Transposition::NO_TRANSPOSE,
          rowN, ncol,
          -1, blockN, rowN,
          io-ncol, 1,
          1, io, 1 ) ;
    switch ( last_block ) {
      case ColrowLU_QR0: Eigen::Map<vec>(io,rowN) = la_solve0.solve(Eigen::Map<vec>(io,rowN)) ; break ;
      case ColrowLU_QR1: Eigen::Map<vec>(io,rowN) = la_solve1.solve(Eigen::Map<vec>(io,rowN)) ; break ;
      case ColrowLU_QR2: Eigen::Map<vec>(io,rowN) = la_solve2.solve(Eigen::Map<vec>(io,rowN)) ; break ;
      case ColrowLU_LU0: Eigen::Map<vec>(io,rowN) = la_solve3.solve(Eigen::Map<vec>(io,rowN)) ; break ;
      case ColrowLU_LU1: Eigen::Map<vec>(io,rowN) = la_solve4.solve(Eigen::Map<vec>(io,rowN)) ; break ;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve( valuePointer in_out ) const {
    if ( use_arceco ) {
      LU_arceco.solve(in_out) ;
      return ;
    }

    // applico permutazione alla RHS
    // col00, col00, < row00, row00 > dimBlock_m_row00, dimBlock_m_row00
    // C      R        C      R       C                 R
    integer const * swapRC = swapRC_blks ;
    valuePointer io = in_out ;

    integer const * swapR = swapRC+col00 ;
    for ( integer k = 0 ; k < col00 ; ++k ) {
      integer k1 = swapR[k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    io     += col00 ;
    swapRC += 2*col00 ;
    for ( nblk = 0 ; nblk < numBlock ; ++nblk )  {
      swapR = swapRC + row00 ;
      for ( integer k = 0 ; k < row00 ; ++k ) {
        integer k1 = swapR[k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
      io     += row00 ;
      swapRC += 2*row00 ;
      swapR   = swapRC + dimBlock_m_row00 ;
      for ( integer k = 0 ; k < dimBlock_m_row00 ; ++k ) {
        integer k1 = swapR[k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
      io     += dimBlock_m_row00 ;
      swapRC += 2*dimBlock_m_row00 ;
    }
    // ultimo blocco
    swapR = swapRC + row00 ;
    for ( integer k = 0 ; k < row00 ; ++k ) {
      integer k1 = swapR[k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;

    // primo blocco
    io = in_out ;
    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          row0,
          block0, row0, io, 1 ) ;

    io += row0 ;
    // blocchi intermedi
    nblk = 0 ;
    while ( nblk < numBlock ) {
      solve_block_L( io ) ;
      io += dimBlock ;
      ++nblk ;
    }

    // soluzione ultimo blocco
    solve_last_block( io ) ;

    while ( nblk > 0 ) {
      --nblk ;
      io -= dimBlock ;
      solve_block_U( io ) ;
    }
    
    // primo blocco
    io -= row0 ;
  
    // soluzione primo blocco
    gemv( Transposition::NO_TRANSPOSE,
          row0, col0-row0,
          -1, block0+row0*row0, row0,
          io+row0, 1,
          1, io, 1 ) ;
    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          row0,
          block0, row0, io, 1 ) ;

    // applico permutazione alla Soluzione
    // col00, col00, < row00, row00 > dimBlock_m_row00, dimBlock_m_row00
    // C      R        C      R       C                 R

    integer const * swapC = swapRC_blks+2*(numBlock*dimBlock+col00) ;
    io = in_out + numBlock*dimBlock + col00 ;
    integer k = row00 ;
    while ( k > 0 ) {
      integer k1 = swapC[--k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    for ( nblk = 0 ; nblk < numBlock ; ++nblk )  {
      io    -= dimBlock_m_row00 ;
      swapC -= 2*dimBlock_m_row00 ;
      k      = dimBlock_m_row00 ;
      while ( k > 0 ) {
        integer k1 = swapC[--k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
      io    -= row00 ;
      swapC -= 2*row00 ;
      k      = row00 ;
      while ( k > 0 ) {
        integer k1 = swapC[--k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
    }
    io    -= col00 ;
    swapC -= 2*col00 ;
    k      = col00 ;
    while ( k > 0 ) {
      integer k1 = swapC[--k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::print( std::ostream & stream ) const {
    cout << "\nneq         = " << neq
         << "\nnumBlock    = " << numBlock
         << "\ndimBlock    = " << dimBlock
         << "\n(row0,col0) = (" << row0 << "," << col0 << ")"
         << "\n(rowN,colN) = (" << rowN << "," << colN << ")"
         << "\ncol00       = " << col00
         << "\ncolNN       = " << colNN ;
    stream << "\nBlock 0\n" ;
    for ( integer i = 0 ; i < row0 ; ++i ) {
      stream << setw(8) << block0[i] ;
      for ( integer j = 1 ; j < col0 ; ++j )
        stream << ' ' << setw(8) << block0[i+j*row0] ;
      stream << '\n' ;
    }
    for ( integer k = 0 ; k < numBlock ; ++k ) {
      stream << "Block " << k+1 << '\n' ;
      valueConstPointer blk = blocks+k*sizeBlock ;
      for ( integer i = 0 ; i < dimBlock ; ++i ) {
        stream << setw(8) << blk[i] ;
        for ( integer j = 1 ; j < 2*dimBlock ; ++j )
          stream << ' ' << setw(8) << blk[i+j*dimBlock] ;
        stream << '\n' ;
      }
    }
    stream << "Block N\n" ;
    for ( integer i = 0 ; i < rowN ; ++i ) {
      stream << setw(8) << blockN[i] ;
      for ( integer j = 1 ; j < colN ; ++j )
        stream << ' ' << setw(8) << blockN[i+j*rowN] ;
      stream << '\n' ;
    }
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::print( std::ostream & stream,
                            integer         row0,
                            integer         col0,
                            t_Value const * block0,
                            integer         numBlock,
                            integer         dimBlock,
                            t_Value const * blocks,
                            integer         rowN,
                            integer         colN,
                            t_Value const * blockN ) {
    integer sizeBlock = 2*dimBlock*dimBlock ;
    stream << "Block 0\n" ;
    for ( integer i = 0 ; i < row0 ; ++i ) {
      stream << setw(8) << block0[i] ;
      for ( integer j = 1 ; j < col0 ; ++j )
        stream << ' ' << setw(8) << block0[i+j*row0] ;
      stream << '\n' ;
    }
    for ( integer k = 0 ; k < numBlock ; ++k ) {
      stream << "Block " << k+1 << '\n' ;
      t_Value const * blk = blocks+k*sizeBlock ;
      for ( integer i = 0 ; i < dimBlock ; ++i ) {
        stream << setw(8) << blk[i] ;
        for ( integer j = 1 ; j < 2*dimBlock ; ++j )
          stream << ' ' << setw(8) << blk[i+j*dimBlock] ;
        stream << '\n' ;
      }
    }
    stream << "Block N\n" ;
    for ( integer i = 0 ; i < rowN ; ++i ) {
      stream << setw(8) << blockN[i] ;
      for ( integer j = 1 ; j < colN ; ++j )
        stream << ' ' << setw(8) << blockN[i+j*rowN] ;
      stream << '\n' ;
    }
  }

  template class ColrowLU<double> ;
  //template class ColrowLU<float> ;

}

///
/// eof: ColrowLU.cc
///

