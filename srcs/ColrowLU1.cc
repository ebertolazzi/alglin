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
  ColrowLU<t_Value>::ColrowLU()
  : baseValue("ColrowLU reals")
  , baseIndex("ColrowLU ints")
  , last_block(ColrowLU_QR0)
  //, last_block(ColrowLU_QR1)
  //, last_block(ColrowLU_QR2)
  //, last_block(ColrowLU_LU0)
  //, last_block(ColrowLU_LU1)
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
  ColrowLU<t_Value>::mv( integer           _numBlock,
                         integer           _dimBlock,
                         integer           _row0,
                         integer           _rowN,
                         valueConstPointer _block0,
                         valueConstPointer _blocks,
                         valueConstPointer _blockN,
                         valueType         alpha,
                         valueConstPointer x,
                         integer           incx,
                         valueType         beta,
                         valuePointer      y,
                         integer           incy ) {

    // first block y = alpha * _block0 * x + beta * y
    gemv( Transposition::NO_TRANSPOSE, _row0, _dimBlock,
          alpha, _block0, _row0,
          x, incx,
          beta, y, incy ) ;

    // internal blocks block
    valueConstPointer xx   = x ;
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
    gemv( Transposition::NO_TRANSPOSE, _rowN, _row0+_rowN,
          alpha, _blockN, _rowN,
          xx, incx,
          beta, yy, incy ) ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::residue( integer           _numBlock,
                              integer           _dimBlock,
                              integer           _row0,
                              integer           _rowN,
                              valueConstPointer _block0,
                              valueConstPointer _blocks,
                              valueConstPointer _blockN,
                              valueConstPointer b,
                              integer           incb,
                              valueConstPointer x,
                              integer           incx,
                              valuePointer      res,
                              integer           incr ) {

    copy( _numBlock*_dimBlock+_row0+_rowN, b, incb, res, incr ) ;
    // first block res -= _block0 * x
    gemv( Transposition::NO_TRANSPOSE, _row0, _dimBlock,
          -1.0, _block0, _row0,
          x, incx,
          1.0, res, incr ) ;

    // internal blocks block
    valueConstPointer xx   = x ;
    valuePointer      rr   = res+_row0*incr ;
    valueConstPointer blks = _blocks ;
    for ( integer i = 0 ; i < _numBlock ; ++i ) {
      gemv( Transposition::NO_TRANSPOSE, _dimBlock, 2*_dimBlock,
            -1.0, blks, _dimBlock,
            xx, incx,
            1.0, rr, incr ) ;
      xx   += _dimBlock*incx ;
      rr   += _dimBlock*incr ;
      blks += 2*_dimBlock*_dimBlock ;
    }

    // last block
    gemv( Transposition::NO_TRANSPOSE, _rowN, _row0+_rowN,
          -1.0, _blockN, _rowN,
          xx, incx,
          1.0, rr, incr ) ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize( integer           _numBlock,
                                integer           _dimBlock,
                                integer           _row0,
                                integer           _rowN,
                                valueConstPointer _block0,
                                valueConstPointer _blocks,
                                valueConstPointer _blockN  ) {
    numBlock        = _numBlock ;
    dimBlock        = _dimBlock ;
    dimBlock_m_row0 = _dimBlock-_row0 ;
    row0            = _row0 ;
    rowN            = _rowN ;
    sizeBlock       = 2*dimBlock*dimBlock ;
    offs            = dimBlock*dimBlock - dimBlock_m_row0 ;
    Nlast           = row0+rowN ;

    // allocate
    baseIndex.allocate( numBlock*(dimBlock+row0)+2*Nlast ) ; // @@@@@@@
    baseValue.allocate( sizeBlock*numBlock + dimBlock*row0 + Nlast*Nlast ) ;
    block0      = baseValue( dimBlock * row0 ) ;
    blocks      = baseValue( sizeBlock * numBlock ) ;
    blockNN     = baseValue( Nlast*Nlast ) ;
    swapRC_blks = baseIndex( numBlock*(dimBlock+row0) +2*Nlast ) ; // @@@@@@@@
    
    // copy block
    copy( dimBlock*row0,      _block0, 1, block0, 1 ) ;
    copy( sizeBlock*numBlock, _blocks, 1, blocks, 1 ) ;

    gezero( row0, Nlast, blockNN, Nlast ) ;
    integer ierr = gecopy( rowN, Nlast, _blockN, rowN, blockNN+row0, Nlast ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize, in gecopy ierr = " << ierr ) ;

    factorize() ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_inplace( integer      _numBlock,
                                        integer      _dimBlock,
                                        integer      _row0,
                                        integer      _rowN,
                                        valuePointer _block0,
                                        valuePointer _blocks,
                                        valuePointer _blockN ) {
    numBlock        = _numBlock ;
    dimBlock        = _dimBlock ;
    dimBlock_m_row0 = _dimBlock-_row0 ;
    row0            = _row0 ;
    rowN            = _rowN ;
    sizeBlock       = 2*dimBlock*dimBlock ;
    offs            = dimBlock*dimBlock - dimBlock_m_row0 ;
    Nlast           = row0+rowN ;

    baseIndex.allocate( numBlock*(dimBlock+row0)+2*Nlast ) ; // @@@@@@@
    baseValue.allocate( Nlast*Nlast ) ;

    blockNN     = baseValue( Nlast*Nlast ) ;
    swapRC_blks = baseIndex( numBlock*(dimBlock+row0) +2*Nlast ) ; // @@@@@@@

    integer ierr = gecopy( rowN, Nlast, _blockN, rowN, blockNN+row0, Nlast ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_inplace, in gecopy ierr = " << ierr ) ;
    factorize() ;
  }
  
  /*
  //  +------+             +-------+             +-------+
  //  |  A   |             | L\ U  |             | L\ U  |
  //  +------+------+ ===> +---+---+------+ ===> +---+---+------+
  //  |             |      |   |          |      |   | \ |      |
  //  |      B      |      | 0 |          |      | 0 | 0 |      |
  //  +-------------+      +---+----------+      +---+----------+
  //
  //  applico permutazione colonne (da applicare alla rovescia sulla soluzione)
  //
  //  +----+          +---+-----+              +---+-----+
  //  |  L |          | U |  N  |              |LU | LN  |
  //  +----+------+   +---+-----+--------+ P = +---+-----+--------+
  //  | M1 | I    |   | 0 |     /M1\     |     |M1U|              |
  //  | M2 |    I |   | 0 |  C- \M2/ N   |     |M2U|     C        |
  //  +----+------+   +---+--------------+     +---+--------------+
  //
  //  applico permutazione righe (da applicare alla rovescia sulla rhs)
  //
  //    +----+          +---+          +---+-----+
  //    |  L |          | I |          | U |  N  |
  //  P +----+------+   +---+---+---+  +---+-----+--------+ P
  //    |  M | I    |   |   | L |   |  |   |  U  |        |
  //    |  M |    I |   |   | W | I |  |   |     |    D   |
  //    +----+------+   +---+---+---+  +---+-----+--------+
  //
  //  A = row0 x dimBlock
  //  B = ldB x ( 2*dimBlock )
  //
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_block( bool first ) {

    // applico row0 passi di Gauss con full pivoting
    valuePointer Bk  = Bmat ;
    valuePointer Ak  = Amat ;
    valuePointer Akk = Amat ;

    integer * swapC = swapRC ;
    integer * swapR = swapRC+row0 ;

    // LU full pivoting
    integer k = 0 ;
    bool do_loop ;
    do {
      // cerco pivot
      integer ipiv = k, jpiv = k ;
      valueType amax = std::abs(Akk[0]) ;
      for ( integer i = k ; i < row0 ; ++i ) {
        //integer i = k ;
        for ( integer j = k ; j < dimBlock ; ++j ) {
          valueType amax1 = std::abs(Amat[i+j*ldA]) ;
          if ( amax1 > amax ) { ipiv = i ; jpiv = j ; amax = amax1 ; }
        }
      }
      swapR[k] = ipiv ;
      swapC[k] = jpiv ; // memorizzo scambio
      if ( ipiv > k ) swap( dimBlock, Amat + k,   ldA, Amat + ipiv,   ldA ) ;
      if ( jpiv > k ) swap( row0,     Amat + k*ldA, 1, Amat + jpiv*ldA, 1 ) ;

      // controllo pivot non zero
      ALGLIN_ASSERT( amax > epsi,
                     "ColrowLU::factorize_block, found pivot: " << Akk[0] << " at block N." <<
                     nblk << " row N." << k << " of the block" ) ;
      // memorizzo L^(-1)
      valuePointer X = Akk + 1 ;
      valuePointer Y = Akk + ldA ;
      // aggiorno k
      ++k ;
      do_loop = k < row0 ;
      if ( do_loop ) {
        rscal( row0-k, Akk[0], X,  1 ) ;
        // applico update rank 1 ==> A := A - x*y'
        Ak  += ldA ;
        Bk  += dimBlock ;
        Akk += ldA+1 ;
        ger( row0-k, dimBlock-k, -1, X, 1, Y, ldA, Akk, ldA ) ;
      }
    } while ( do_loop ) ;

    // propago permutazioni ai vari blocchi
    if ( first ) {
      for ( integer j = 0 ; j < row0 ; ++j ) {
        integer j1 = swapC[j] ;
        if ( j1 > j ) swap( dimBlock, Bmat + j*dimBlock,  1,
                                      Bmat + j1*dimBlock, 1 ) ;
      }
    } else {
      valuePointer Tmat = Amat - dimBlock_m_row0 ;
      valuePointer Lmat = Amat - dimBlock*dimBlock ;
      for ( integer j = 0 ; j < row0 ; ++j ) {
        integer j1 = swapC[j] ;
        if ( j1 > j ) {
          integer offs1 = j*dimBlock ;
          integer offs2 = j1*dimBlock ;
          swap( dimBlock,        Bmat + offs1, 1, Bmat + offs2, 1 ) ;
          swap( dimBlock_m_row0, Tmat + offs1, 1, Tmat + offs2, 1 ) ;
        }
        j1 = swapR[j] ;
        if ( j1 > j ) swap( dimBlock, Lmat + j, dimBlock, Lmat + j1, dimBlock ) ;
      }
    }

    /*
    //  +---------+
    //  | A1   A2 |
    //  +---------+--------+
    //  | B1   C1 |  D1    |
    //  | B2   C2 |  D2    |
    //  +---------+--------+
    */

    valuePointer A1 = Amat ;
    valuePointer A2 = Amat+ldA*row0 ;

    valuePointer B  = Bmat ;

    valuePointer C  = Bmat+row0*dimBlock ;
    valuePointer C1 = C ;
    valuePointer C2 = C1+dimBlock_m_row0 ;

    valuePointer D  = Bmat+dimBlock*dimBlock ;
    valuePointer D1 = D ;
    valuePointer D2 = D+dimBlock_m_row0 ;

    // B * Upper(A1)^(-1)
    trsm( SideMultiply::RIGHT,
          ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          dimBlock, row0, 1.0,
          A1, ldA,
          B,  dimBlock ) ;
    // C -= B * A2
    gemm( Transposition::NO_TRANSPOSE,
          Transposition::NO_TRANSPOSE,
          dimBlock, dimBlock_m_row0, row0,
          -1.0,
          B,  dimBlock,
          A2, ldA,
          1.0,
          C, dimBlock ) ;
    // ----------
    // C = L*U
    swapR = swapRC + 2*row0 ;
    integer ierr = getrf( dimBlock, dimBlock_m_row0, C, dimBlock, swapR ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_block, found ierr = " << ierr <<
                   " at LU factorizatioon on block N." << nblk+1 ) ;
    ierr = swaps( dimBlock, D, dimBlock, 0, dimBlock_m_row0-1, swapR, 1 ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_block, found ierr = " << ierr <<
                   " at LU swaps rows on block N." << nblk+1 ) ;
    ierr = swaps( row0, B, dimBlock, 0, dimBlock_m_row0-1, swapR, 1 ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_block, found ierr = " << ierr <<
                   " at LU swaps rows on block N." << nblk+1 ) ;
    // Lower(C1)^(-1) * D1
    trsm( SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          dimBlock_m_row0, dimBlock, 1.0,
          C1, dimBlock,
          D1, dimBlock ) ;
    // D2 -= C2*D1
    gemm( Transposition::NO_TRANSPOSE,
          Transposition::NO_TRANSPOSE,
          row0, dimBlock, dimBlock_m_row0,
          -1.0,
          C2, dimBlock,
          D1, dimBlock,
          1.0,
          D2, dimBlock ) ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize() {
    // primo blocco
    swapRC = swapRC_blks ;
    Amat   = block0 ;
    ldA    = row0 ;
    Bmat   = blocks ;
    nblk   = 0 ;
    factorize_block(true) ;

    // blocchi intermedi (n-1)
    ldA = dimBlock ;
    while ( ++nblk < numBlock ) {
      swapRC += dimBlock+row0 ;
      Bmat   += sizeBlock ;
      Amat    = Bmat-offs ;
      factorize_block(false) ;
    }

    swapRC += dimBlock+row0 ;
    Amat    = Bmat + sizeBlock - offs ;

    // copio in ultimo blocco
    integer N = row0+rowN ;
    integer ierr = gecopy( row0, dimBlock, Amat, dimBlock, blockNN, N ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize, in gecopy ierr = " << ierr ) ;

    // fattorizzazione ultimo blocco
    switch ( last_block ) {
      case ColrowLU_QR0: la_solve0.compute(Eigen::Map<mat>(blockNN,N,N)) ; break ;
      case ColrowLU_QR1: la_solve1.compute(Eigen::Map<mat>(blockNN,N,N)) ; break ;
      case ColrowLU_QR2: la_solve2.compute(Eigen::Map<mat>(blockNN,N,N)) ; break ;
      case ColrowLU_LU0: la_solve3.compute(Eigen::Map<mat>(blockNN,N,N)) ; break ;
      case ColrowLU_LU1: la_solve4.compute(Eigen::Map<mat>(blockNN,N,N)) ; break ;
    }

    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize, failed to factorize last block ierr = " << ierr ) ;
  }

  /*
  //  applico permutazione righe (da applicare alla rovescia sulla rhs)
  //
  //  +----+          +---+          +---+-----+
  //  |  L |          | I |          | U |  N  |
  //  +----+------+   +---+---+---+  +---+-----+--------+
  //  |  M | I    |   |   | L |   |  |   |  U  |        |
  //  |  M |    I |   |   | W | I |  |   |     |    D   |
  //  +----+------+   +---+---+---+  +---+-----+--------+
  //
  */
  template <typename t_Value>
  void
  ColrowLU<t_Value>::solver_block_L( valuePointer io ) const {

    valuePointer io1 = io+row0 ;

    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          row0,
          Amat, ldA, io, 1 ) ;
    gemv( Transposition::NO_TRANSPOSE,
          dimBlock, row0,
          -1, Bmat, dimBlock,
          io, 1,
          1, io1, 1 ) ;
    valuePointer io2 = io+dimBlock ;
    valuePointer L   = Bmat+dimBlock*row0 ;
    valuePointer W   = L+dimBlock_m_row0 ;

    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          dimBlock_m_row0,
          L, dimBlock, io1, 1 ) ;
    gemv( Transposition::NO_TRANSPOSE,
          row0, dimBlock_m_row0,
          -1, W, dimBlock,
          io1, 1,
          1, io2, 1 ) ;
  }

  /*
  //  applico permutazione righe (da applicare alla rovescia sulla rhs)
  //
  //                                  io   io1   io2
  //  +----+          +---+          +---+-----+
  //  |  L |          | I |          | A |  N  |
  //  +----+------+   +---+---+---+  +---+-----+--------+
  //  |  M | I    |   |   | L |   |  |   |  U  |   D    |
  //  |  M |    I |   |   | W | I |  |   |     |        |
  //  +----+------+   +---+---+---+  +---+-----+--------+
  //
  */

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solver_block_U( valuePointer in_out ) const {

    // y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    valuePointer      io1 = in_out+row0 ;
    valuePointer      io2 = in_out+dimBlock ;
    valueConstPointer D   = Bmat+dimBlock*dimBlock ;
    valueConstPointer U   = Bmat+row0*dimBlock ;
    valueConstPointer N   = Amat+row0*ldA ;

    // io1 -= D*io2
    gemv( Transposition::NO_TRANSPOSE,
          dimBlock_m_row0, dimBlock,
          -1, D, dimBlock,
          io2, 1,
          1, io1, 1 ) ;

    // io1 = U^(-1)*io1
    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          dimBlock_m_row0,
          U, dimBlock, io1, 1 ) ;

    // io1 -= N*io1
    gemv( Transposition::NO_TRANSPOSE,
          row0, dimBlock_m_row0,
          -1, N, ldA,
          io1, 1,
          1, in_out, 1 ) ;

    // io = U(A)^(-1)*io
    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          row0,
          Amat, ldA, in_out, 1 ) ;
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve( valuePointer in_out ) const {

    // applico permutazione alla RHS
    for ( integer nb = 0 ; nb < numBlock ; ++nb )  {
      integer * swapRC = swapRC_blks+nb*(dimBlock+row0) ;
      valuePointer io = in_out+nb*dimBlock ;

      integer const * swapR = swapRC+row0 ;
      for ( integer k = 0 ; k < row0 ; ++k ) {
        integer k1 = swapR[k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;

      // applico permutazione
      swaps( 1, io+row0, dimBlock, 0, dimBlock_m_row0-1, swapRC+2*row0, 1 ) ;
    }

    valuePointer io = in_out ;

    // primo blocco
    swapRC = swapRC_blks ;
    Amat   = block0 ;
    ldA    = row0 ;
    Bmat   = blocks ;
    solver_block_L( io ) ;

    // blocchi intermedi
    ldA = dimBlock ;
    for ( integer i = 1 ; i < numBlock ; ++i ) {
      swapRC += dimBlock+row0 ;
      Bmat   += sizeBlock ;
      Amat    = Bmat-offs ;
      io     += dimBlock ;
      solver_block_L( io ) ;
    }

    // risolve rispetto ultimo blocco
    integer N = row0+rowN ;
    swapRC += dimBlock+row0 ;
    Bmat   += sizeBlock ;
    io     += dimBlock ;

    // fattorizzazione ultimo blocco
    switch ( last_block ) {
      case ColrowLU_QR0: Eigen::Map<vec>(io,N) = la_solve0.solve(Eigen::Map<vec>(io,N)) ; break ;
      case ColrowLU_QR1: Eigen::Map<vec>(io,N) = la_solve1.solve(Eigen::Map<vec>(io,N)) ; break ;
      case ColrowLU_QR2: Eigen::Map<vec>(io,N) = la_solve2.solve(Eigen::Map<vec>(io,N)) ; break ;
      case ColrowLU_LU0: Eigen::Map<vec>(io,N) = la_solve3.solve(Eigen::Map<vec>(io,N)) ; break ;
      case ColrowLU_LU1: Eigen::Map<vec>(io,N) = la_solve4.solve(Eigen::Map<vec>(io,N)) ; break ;
    }

    for ( integer i = 1 ; i < numBlock ; ++i ) {
      swapRC -= dimBlock+row0 ;
      Bmat   -= sizeBlock ;
      Amat    = Bmat-offs ;
      io     -= dimBlock ;
      solver_block_U( io ) ;
    }

    swapRC = swapRC_blks ;
    Amat   = block0 ;
    ldA    = row0 ;
    Bmat   = blocks ;
    solver_block_U( in_out ) ;

    for ( integer nb = 0 ; nb < numBlock ; ++nb )  {
      integer const * swapC = swapRC_blks+nb*(dimBlock+row0) ;
      valuePointer io = in_out+nb*dimBlock ;
      integer k = row0 ;
      do {
        integer k1 = swapC[--k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } while ( k > 0 ) ;
    }

  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::print( std::ostream & stream ) const {
    stream << "Block 0\n" ;
    for ( integer i = 0 ; i < row0 ; ++i ) {
      stream << setw(8) << block0[i] ;
      for ( integer j = 1 ; j < dimBlock ; ++j )
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
    integer N = rowN+row0 ;
    for ( integer i = 0 ; i < N ; ++i ) {
      stream << setw(8) << blockNN[i] ;
      for ( integer j = 1 ; j < N ; ++j )
        stream << ' ' << setw(8) << blockNN[i+j*N] ;
      stream << '\n' ;
    }
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::print_colrow( std::ostream & stream,
                                   integer         numBlock,
                                   integer         dimBlock,
                                   integer         row0,
                                   integer         rowN,
                                   t_Value const * block0,
                                   t_Value const * blocks,
                                   t_Value const * blockN ) {
    integer sizeBlock = 2*dimBlock*dimBlock ;
    stream << "Block 0\n" ;
    for ( integer i = 0 ; i < row0 ; ++i ) {
      stream << setw(8) << block0[i] ;
      for ( integer j = 1 ; j < dimBlock ; ++j )
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
      for ( integer j = 1 ; j < row0+rowN ; ++j )
        stream << ' ' << setw(8) << blockN[i+j*rowN] ;
      stream << '\n' ;
    }
  }
  
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::print_colrow_to_maple( std::ostream & stream,
                                            integer         numBlock,
                                            integer         dimBlock,
                                            integer         row0,
                                            integer         rowN,
                                            t_Value const * block0,
                                            t_Value const * blocks,
                                            t_Value const * blockN ) {

    integer N = row0+rowN+numBlock*dimBlock ;
    std::vector<t_Value> mat(N*N) ;
    std::fill(mat.begin(),mat.end(),0) ;

    integer sizeBlock = 2*dimBlock*dimBlock ;
    for ( integer i = 0 ; i < row0 ; ++i ) {
      for ( integer j = 0 ; j < dimBlock ; ++j )
        mat[i+j*N] = block0[i+j*row0] ;
    }
    for ( integer k = 0 ; k < numBlock ; ++k ) {
      integer offs = k*dimBlock ;
      t_Value const * blk = blocks+k*sizeBlock ;
      for ( integer i = 0 ; i < dimBlock ; ++i ) {
        for ( integer j = 0 ; j < 2*dimBlock ; ++j )
          mat[i+offs+row0+(j+offs)*N] = blk[i+j*dimBlock] ;
      }
    }
    integer offs = numBlock*dimBlock ;
    for ( integer i = 0 ; i < rowN ; ++i ) {
      for ( integer j = 0 ; j < row0+rowN ; ++j )
        mat[i+offs+row0+(j+offs)*N] = blockN[i+j*rowN] ;
    }
    stream << "M :=\n<" ;
    for ( integer j = 0 ; j < N ; ++j ) {
      stream << "<" << mat[j*N] ;
      for ( integer i = 1 ; i < N ; ++i ) stream << "," << mat[i+j*N] ;
      if ( j < N-1 ) stream << ">|\n" ;
      else           stream << ">>;\n" ;
    }
  }

  template class ColrowLU<double> ;
  template class ColrowLU<float> ;

}

///
/// eof: ColrowLU.cc
///

