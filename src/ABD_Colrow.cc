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
/// file: ABD_Colrow.cc
///

#include "ABD_Colrow.hh"
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
  : baseValue("ColrowLU_values")
  , baseIndex("ColrowLU_integers")
  , last_block(COLROW_LASTBLOCK_LU)
  //, last_block(COLROW_LASTBLOCK_QR)
  //, last_block(COLROW_LASTBLOCK_SVD)
  , NB(25)
  {
  }
  
  template <typename t_Value>
  ColrowLU<t_Value>::~ColrowLU() {
    baseValue.free() ;
    baseIndex.free() ;
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
  ColrowLU<t_Value>::factorize( COLROW_LASTBLOCK_Choice choice,
                                // ----------------------------
                                integer           _row0,
                                integer           _col0,
                                valueConstPointer _block0,
                                // ----------------------------
                                integer           _numBlock,
                                integer           _dimBlock,
                                valueConstPointer _blocks,
                                // ----------------------------
                                integer           _rowN,
                                integer           _colN,
                                valueConstPointer _blockN ) {
    last_block = choice ;

    numBlock = _numBlock ;
    dimBlock = _dimBlock ;
    row0     = _row0 ;
    col0     = _col0 ;
    rowN     = _rowN ;
    colN     = _colN ;

    setup() ;

    // allocate
    baseIndex.allocate(size_t( numBlock*dimBlock+row0 )) ;
    baseValue.allocate(size_t( sizeBlock*numBlock + row0*col0 + rowN*colN )) ;
    block0      = baseValue(size_t( row0*col0 )) ;
    blocks      = baseValue(size_t( sizeBlock * numBlock )) ;
    blockN      = baseValue(size_t( rowN*colN )) ;
    swapRC_blks = baseIndex(size_t( numBlock*dimBlock+row0 )) ;

    // copy block
    copy( row0*col0,          _block0, 1, block0, 1 ) ;
    copy( sizeBlock*numBlock, _blocks, 1, blocks, 1 ) ;
    copy( rowN*colN,          _blockN, 1, blockN, 1 ) ;

    factorize() ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize_inplace( COLROW_LASTBLOCK_Choice choice,
                                        // ----------------------------
                                        integer      _row0,
                                        integer      _col0,
                                        valuePointer _block0,
                                        // ----------------------------
                                        integer      _numBlock,
                                        integer      _dimBlock,
                                        valuePointer _blocks,
                                        // ----------------------------
                                        integer      _rowN,
                                        integer      _colN,
                                        valuePointer _blockN ) {
    last_block = choice ;

    numBlock = _numBlock ;
    dimBlock = _dimBlock ;
    row0     = _row0 ;
    col0     = _col0 ;
    rowN     = _rowN ;
    colN     = _colN ;
    
    block0   = _block0 ;
    blocks   = _blocks ;
    blockN   = _blockN ;

    setup() ;

    baseIndex.allocate(size_t( numBlock*dimBlock+row0 )) ;
    baseValue.allocate(size_t( sizeBlock*numBlock + row0*col0 + rowN*colN )) ;
    swapRC_blks = baseIndex(size_t( numBlock*dimBlock+row0 )) ;

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
                                    integer swapR[] ) {
#if 1
    integer ierr ;
    if ( 2*NB < nrA ) ierr = getry( nrA, ncA, A, ldA, swapR, NB ) ;
    else              ierr = gty( nrA, ncA, A, ldA, swapR ) ;

    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::LU_left_right, found ierr: " << ierr ) ;
    // applico permutazione al blocco
    t_Value * R = A + ncA * ldA ;
    t_Value * L = A - ncL * ldA ;
    for ( integer i = 0 ; i < ncA ; ++i ) {
      integer ip = swapR[i] ;
      if ( ip > i ) {
        swap( ncR, R + i, ldA, R + ip, ldA ) ;
        swap( ncL, L + i, ldA, L + ip, ldA ) ;
      }
    }
    trsm( LEFT, LOWER, NO_TRANSPOSE, UNIT,
          ncA, ncR, 1.0,
          A, ldA,
          R, ldA ) ;
    gemm( NO_TRANSPOSE, NO_TRANSPOSE,
          nrA-ncA, ncR, ncA,
          -1.0,  A+ncA, ldA,
          R, ldA,
          1.0,  R+ncA, ldA ) ;

#else
    t_Value * Ak  = A ;
    t_Value * Akk = A ;
    t_Value * L   = A - ncL*ldA ;

    integer ncc = ncA+ncL+ncR ;

    // LU row pivoting
    integer rr = nrA ;
    integer cc = ncA+ncR ;
    for ( integer k = 0 ; k < ncA ; ++k ) {
      // cerco pivot
      integer ipiv = k ;
      t_Value amax = std::abs(Akk[0]) ;
      // solo scambi di riga
      for ( integer i = k ; i < nrA ; ++i ) {
        t_Value amax1 = std::abs(Ak[i]) ;
        if ( amax1 > amax ) { ipiv = i ; amax = amax1 ; }
      }
      // controllo pivot non zero
      ALGLIN_ASSERT( amax > epsi,
                     "ColrowLU::LU_left_right, found pivot: " << Akk[0] <<
                     " at block nblk = " << nblk ) ;
      // memorizzo scambio
      swapR[k] = ipiv ;
      // applico permutazione
      if ( ipiv > k ) swap( ncc, L + k, ldA, L + ipiv, ldA ) ; // scambio righe
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
  #endif
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
                                    integer swapC[] ) {
#if 1
    integer ierr ;
    if ( 2*NB < ncA ) ierr = getrx( nrA, ncA, A, ldA, swapC, NB ) ;
    else              ierr = gtx( nrA, ncA, A, ldA, swapC ) ;

    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::LU_top_bottom, found ierr: " << ierr ) ;
    // applico permutazione al blocco
    t_Value * T = A - nrT ;
    for ( integer i = 0 ; i < nrA ; ++i ) {
      integer ip = swapC[i] ;
      if ( ip > i ) {
        swap( nrT, T + i*ldA, 1, T + ip*ldA, 1 ) ;
        swap( nrB, B + i*ldB, 1, B + ip*ldB, 1 ) ;
      }
    }
    trsm( RIGHT, UPPER, NO_TRANSPOSE, NON_UNIT,
          nrB, nrA, 1.0,
          A, ldA,
          B, ldB ) ;
    gemm( NO_TRANSPOSE, NO_TRANSPOSE,
          nrB, ncA-nrA, nrA,
          -1.0, B, ldB,
          A+nrA*ldA, ldA,
          1.0,  B+nrA*ldB, ldB ) ;

#else
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
      integer jpiv = k ;
      t_Value amax = std::abs(Akk[0]) ;
      for ( integer j = k ; j < ncA ; ++j ) {
        t_Value amax1 = std::abs(A[k+j*ldA]) ;
        if ( amax1 > amax ) { jpiv = j ; amax = amax1 ; }
      }
      // controllo pivot non zero
      ALGLIN_ASSERT( amax > epsi,
                     "ColrowLU::LU_top_bottom, found pivot: " << Akk[0] <<
                     " at block nblk = " << nblk ) ;
      // memorizzo scambio
      swapC[k] = jpiv ;
      // applico permutazione
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
#endif
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize() {
    // primo blocco
    integer * swapRC = swapRC_blks ;
    if ( col00 > 0 ) {
      LU_left_right( row0, col00, 0, col0-col00, block0, row0, swapRC ) ;
      swapRC += col00 ;
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

    // blocchi intermedi (n-1)
    if ( numBlock > 0 ) {
      valuePointer B1 = block0 + (row0+1)*col00 ;
      LU_top_bottom( col00, row00,
                     dimBlock, B1, row0,
                     dimBlock, blocks, dimBlock, swapRC ) ;
      swapRC += row00 ;
      valuePointer D = blocks + row00 * dimBlock ;
      //                 NR          NC            L       R
      LU_left_right( dimBlock, dimBlock_m_row00, row00, dimBlock,
                     D, dimBlock, swapRC ) ;
      swapRC += dimBlock_m_row00 ;
    }
    
    valuePointer C = blocks ;
    for ( nblk = 1 ; nblk < numBlock ; ++nblk ) {
      C += sizeBlock;
      valuePointer B1 = C - hSizeBlock + dimBlock_m_row00 ;
      LU_top_bottom( dimBlock_m_row00, row00,
                     dimBlock, B1, dimBlock,
                     dimBlock,  C, dimBlock, swapRC ) ;
      swapRC += row00 ;
      valuePointer D = C + row00 * dimBlock ;
      //                 NR          NC            L       R
      LU_left_right( dimBlock, dimBlock_m_row00, row00, dimBlock,
                     D, dimBlock, swapRC ) ;
      swapRC += dimBlock_m_row00 ;
    }

    // ultimo blocco
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

    swapRC = swapRC_blks + ( col00 + numBlock*dimBlock ) ;

    if ( numBlock == 0 ) {
      valuePointer B1 = block0 + (row0+1) * col00 ;
      LU_top_bottom( col00, row00, dimBlock,
                     B1, row0,
                     rowN, blockN, rowN,
                     swapRC ) ;
    } else {
      valuePointer B1 = blocks + numBlock * sizeBlock - hSizeBlock + dimBlock_m_row00 ;
      LU_top_bottom( dimBlock_m_row00, row00, dimBlock,
                     B1, dimBlock,
                     rowN, blockN, rowN,
                     swapRC ) ;
    }

    // fattorizzazione ultimo blocco
    valuePointer D0 = blockN + row00 * rowN;
    switch ( last_block ) {
      case COLROW_LASTBLOCK_LU:  la_lu.factorize(rowN,rowN,D0,rowN)  ; break ;
      case COLROW_LASTBLOCK_QR:  la_qr.factorize(rowN,rowN,D0,rowN)  ; break ;
      case COLROW_LASTBLOCK_SVD: la_svd.factorize(rowN,rowN,D0,rowN) ; break ;
    }
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve( valuePointer in_out ) const {

    // applico permutazione alla RHS
    integer const * swapR = swapRC_blks ;
    valuePointer io = in_out ;
    for ( integer k = 0 ; k < col00 ; ++k ) {
      integer k1 = swapR[k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    io    += row0 ;
    swapR += row0 ;
    for ( nblk = 0 ; nblk < numBlock ; ++nblk )  {
      for ( integer k = 0 ; k < dimBlock_m_row00 ; ++k ) {
        integer k1 = swapR[k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
      io    += dimBlock ;
      swapR += dimBlock ;
    }

    // primo blocco
    io = in_out ;
    trsv( LOWER, NO_TRANSPOSE, UNIT,
          row0,
          block0, row0, io, 1 ) ;

    io += row0 ;
    // blocchi intermedi
    nblk = 0 ;
    while ( nblk < numBlock ) {

      /*
      //
      //  +----+---+---+---+
      //  |  M | L |   :   |
      //  |  M | L | L :   |
      //  +----+---+---+---+
      //  row00
      */

      valuePointer io1 = io - row00 ;
      valuePointer M   = blocks + nblk * sizeBlock ;
      valuePointer L   = M + row00 * dimBlock ;

      // io -= M*io1
      gemv( NO_TRANSPOSE,
            dimBlock, row00,
            -1, M, dimBlock,
            io1, 1,
            1, io, 1 ) ;

      trsv( LOWER, NO_TRANSPOSE, UNIT,
            dimBlock,
            L, dimBlock, io, 1 ) ;

      io += dimBlock ;
      ++nblk ;
    }

    // soluzione ultimo blocco
    integer ncol = colN-rowN ;
    gemv( NO_TRANSPOSE,
          rowN, ncol,
          -1, blockN, rowN,
          io-ncol, 1,
          1, io, 1 ) ;

    switch ( last_block ) {
      case COLROW_LASTBLOCK_LU:  la_lu.solve(io)     ; break ;
      case COLROW_LASTBLOCK_QR:  la_qr.solve(io)     ; break ;
      case COLROW_LASTBLOCK_SVD: la_svd.solve(io,io) ; break ;
    }

    while ( nblk > 0 ) {
      --nblk ;
      io -= dimBlock ;

      /*
      //  +---------+---+-------+
      //  |       U |   |       |
      //  |         | U |   M   |
      //  +---------+---+-------+
      //  row00
      */

      valuePointer io1 = io + dimBlock ;
      valuePointer U   = blocks + nblk * sizeBlock + row00 * dimBlock ;
      valuePointer M   = U + hSizeBlock ;

      gemv( NO_TRANSPOSE,
            dimBlock, dimBlock_m_row00,
            -1, M, dimBlock,
            io1, 1,
            1, io, 1 ) ;

      trsv( UPPER, NO_TRANSPOSE, NON_UNIT,
            dimBlock,
            U, dimBlock, io, 1 ) ;
    }
    
    // primo blocco
    io -= row0 ;
  
    // soluzione primo blocco
    gemv( NO_TRANSPOSE,
          row0, col0-row0,
          -1, block0+row0*row0, row0,
          io+row0, 1,
          1, io, 1 ) ;

    trsv( UPPER, NO_TRANSPOSE, NON_UNIT,
          row0,
          block0, row0, io, 1 ) ;

    // applico permutazione alla Soluzione
    integer const * swapC = swapRC_blks+(numBlock*dimBlock+col00) ;
    io = in_out + numBlock*dimBlock + col00 ;
    integer k = row00 ;
    while ( k > 0 ) {
      integer k1 = swapC[--k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    for ( nblk = 0 ; nblk < numBlock ; ++nblk )  {
      io    -= dimBlock ;
      swapC -= dimBlock ;
      k      = row00 ;
      while ( k > 0 ) {
        integer k1 = swapC[--k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
    }
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
  template class ColrowLU<float> ;

}

///
/// eof: ColrowLU.cc
///

