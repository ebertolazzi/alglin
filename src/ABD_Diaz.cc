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

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

///
/// file: ABD_Diaz.cc
///

#include "ABD_Diaz.hh"
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

  template <> float  const DiazLU<float>::epsi  = 100*std::numeric_limits<float>::epsilon() ;
  template <> double const DiazLU<double>::epsi = 1000*std::numeric_limits<double>::epsilon() ;

  template <typename t_Value>
  DiazLU<t_Value>::DiazLU()
  : baseValue("DiazLU_values")
  , baseIndex("DiazLU_integers")
  , NB(25)
  {
  }
  
  template <typename t_Value>
  DiazLU<t_Value>::~DiazLU() {
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
  //           |                     | n
  //           |                     |
  //           +----------+----------+
  //                 2*n
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
  /*\
   |   _                 _ _____           ____        _   _
   |  | | ___   __ _  __| |_   _|__  _ __ | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` | | |/ _ \| '_ \|  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | | | (_) | |_) | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_| |_|\___/| .__/|____/ \___/ \__|\__\___/|_| |_| |_|
   |                                |_|
  \*/
  template <typename t_Value>
  void
  DiazLU<t_Value>::loadTopBottom( integer           _row0,
                                  integer           _col0,
                                  valueConstPointer _block0,
                                  integer           ld0,
                                  // ----------------------------
                                  integer           _rowN,
                                  integer           _colN,
                                  valueConstPointer _blockN,
                                  integer           ldN ) {

    integer & n      = this->n ;
    integer & nblock = this->nblock ;

    row0 = _row0 ;
    col0 = _col0 ;
    rowN = _rowN ;
    colN = _colN ;

    col00 = col0 - n ;
    colNN = colN - n ;
    row00 = row0 - col00 ;
    rowNN = rowN - colNN ;

    n_m_row00 = n - row00 ;

    ALGLIN_ASSERT( row0 > 0 && rowN > 0 && row00 >= 0 && col00 >= 0 &&
                   col0 >= n && colN >= n && n > 0 &&
                   row0 + rowN + n == col0 + colN,
                   "Bad parameter(s):" <<
                   "\nrow0     = " << row0 <<
                   "\ncol0     = " << col0 <<
                   "\ncol00    = " << col00 <<
                   "\nrowN     = " << rowN <<
                   "\ncolN     = " << colN <<
                   "\ncolNN    = " << colNN <<
                   "\nn        = " << n ) ;

    // allocate
    baseIndex.allocate(size_t( nblock*n+row0 )) ;
    baseValue.allocate(size_t( row0*col0 + rowN*colN )) ;
    block0      = baseValue(size_t( row0*col0 )) ;
    blockN      = baseValue(size_t( rowN*colN )) ;
    swapRC_blks = baseIndex(size_t( nblock*n+row0 )) ;

    // copy block
    gecopy( row0, col0, _block0, ld0, block0, row0 ) ;
    gecopy( rowN, colN, _blockN, ldN, blockN, rowN ) ;
  }
  
  /*\
   |  _                 _ ____   ____
   | | | ___   __ _  __| | __ ) / ___|
   | | |/ _ \ / _` |/ _` |  _ \| |
   | | | (_) | (_| | (_| | |_) | |___
   | |_|\___/ \__,_|\__,_|____/ \____|
  \*/
  template <typename t_Value>
  void
  DiazLU<t_Value>::loadBC( integer numInitialBc,
                           integer numFinalBc,
                           integer numCyclicBC,
                           // ----------------------
                           integer numInitialETA,
                           integer numFinalETA,
                           integer numCyclicOMEGA,
                           // ----------------------
                           valueConstPointer H0, integer ld0,
                           valueConstPointer HN, integer ldN,
                           valueConstPointer Hq, integer ldQ ) {
    
    integer & n      = this->n ;
    integer & nblock = this->nblock ;

    integer nq = numInitialBc  + numFinalBc  + numCyclicBC ;
    integer q  = numInitialETA + numFinalETA + numCyclicOMEGA ;

    ALGLIN_ASSERT( n+q == nq,
                   "Bad parameter(s):" <<
                   "\nn              = " << n <<
                   "\nnumInitialBc   = " << numInitialBc   <<
                   "\nnumFinalBc     = " << numFinalBc     <<
                   "\nnumCyclicBC    = " << numCyclicBC    <<
                   "\nnumInitialETA  = " << numInitialETA  <<
                   "\nnumFinalETA    = " << numFinalETA    <<
                   "\nnumCyclicOMEGA = " << numCyclicOMEGA ) ;
    
    ALGLIN_ASSERT( numCyclicOMEGA == 0 && numCyclicBC == 0,
                   "DiazLU cannot manage cyclic BC" ) ;

    col00 = numInitialETA ;
    colNN = numFinalETA ;

    row0 = numInitialBc ;
    col0 = n + col00 ;
    rowN = numFinalBc ;
    colN = n + colNN ;

    row00 = row0 - col00 ;
    rowNN = rowN - colNN ;

    // allocate
    baseIndex.allocate(size_t( nblock*n+row0 )) ;
    baseValue.allocate(size_t( row0*col0 + rowN*colN )) ;
    block0      = baseValue(size_t( row0*col0 )) ;
    blockN      = baseValue(size_t( rowN*colN )) ;
    swapRC_blks = baseIndex(size_t( nblock*n+row0 )) ;

    // zeros block
    zero( row0*col0 + rowN*colN, block0, 1 ) ;

    #define IDX0(I,J) ((I)+(J)*row0)
    #define IDXN(I,J) ((I)+(J)*rowN)
    //  +-+-------+
    //  |x|  NZ   |
    //  |x|  NZ   |
    //  +-+-------+
    gecopy( row0, col00, Hq+rowN+colNN*nq, ldQ, block0+IDX0(0,0),     row0 ) ;
    gecopy( row0, n,     H0+rowN,          ld0, block0+IDX0(0,col00), row0 ) ;
    //  +-------+-+
    //  |  NZ   |y|
    //  +-------+-+
    gecopy( rowN, n,     HN, ldN, blockN+IDXN(0,0), rowN ) ;
    gecopy( rowN, colNN, Hq, ldQ, blockN+IDXN(0,n), rowN ) ;
    #undef IDX0
    #undef IDXN
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
  DiazLU<t_Value>::LU_left_right( integer nrA,
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
                   "DiazLU::LU_left_right, found ierr: " << ierr ) ;
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
                     "DiazLU::LU_left_right, found pivot: " << Akk[0] <<
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
  DiazLU<t_Value>::LU_top_bottom( integer nrT,
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
                   "DiazLU::LU_top_bottom, found ierr: " << ierr ) ;
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
                     "DiazLU::LU_top_bottom, found pivot: " << Akk[0] <<
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

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  DiazLU<t_Value>::factorize() {

    integer & n      = this->n ;
    integer & nxnx2  = this->nxnx2 ;
    integer & nxn    = this->nxn      ;
    integer & nblock = this->nblock ;

    // primo blocco
    integer * swapRC = swapRC_blks ;
    if ( col00 > 0 ) {
      LU_left_right( row0, col00, 0, col0-col00, block0, row0, swapRC ) ;
      swapRC += col00 ;
    }

    /*
    //    dimA   n
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
    if ( nblock > 0 ) {
      valuePointer B1 = block0 + (row0+1)*col00 ;
      LU_top_bottom( col00, row00,
                     n, B1, row0,
                     n, this->AdAu_blk, n, swapRC ) ;
      swapRC += row00 ;
      valuePointer D = this->AdAu_blk + row00 * n ;
      //                 NR          NC            L       R
      LU_left_right( n, n_m_row00, row00, n, D, n, swapRC ) ;
      swapRC += n_m_row00 ;
    }
    
    valuePointer C = this->AdAu_blk ;
    for ( nblk = 1 ; nblk < nblock ; ++nblk ) {
      C += nxnx2;
      valuePointer B1 = C - nxn + n_m_row00 ;
      LU_top_bottom( n_m_row00, row00,
                     n, B1, n,
                     n,  C, n, swapRC ) ;
      swapRC += row00 ;
      valuePointer D = C + row00 * n ;
      //                 NR          NC            L       R
      LU_left_right( n, n_m_row00, row00, n, D, n, swapRC ) ;
      swapRC += n_m_row00 ;
    }

    // ultimo blocco
    /*
    //    dimA   n
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

    swapRC = swapRC_blks + ( col00 + nblock*n ) ;

    if ( nblock == 0 ) {
      valuePointer B1 = block0 + (row0+1) * col00 ;
      LU_top_bottom( col00, row00, n,
                     B1, row0,
                     rowN, blockN, rowN,
                     swapRC ) ;
    } else {
      valuePointer B1 = this->AdAu_blk + nblock * nxnx2 - nxn + n_m_row00 ;
      LU_top_bottom( n_m_row00, row00, n,
                     B1, n,
                     rowN, blockN, rowN,
                     swapRC ) ;
    }

    // fattorizzazione ultimo blocco
    valuePointer D0 = blockN + row00 * rowN;
    this->la_factorization->factorize(rowN,rowN,D0,rowN) ;
  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  DiazLU<t_Value>::solve( valuePointer in_out ) const {
  
    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;

    integer neq = nblock*n+row0+rowN ;
    std::rotate( in_out, in_out + neq - row0, in_out + neq ) ;

    // applico permutazione alla RHS
    integer const * swapR = swapRC_blks ;
    valuePointer io = in_out ;
    for ( integer k = 0 ; k < col00 ; ++k ) {
      integer k1 = swapR[k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    io    += row0 ;
    swapR += row0 ;
    for ( nblk = 0 ; nblk < nblock ; ++nblk )  {
      for ( integer k = 0 ; k < n_m_row00 ; ++k ) {
        integer k1 = swapR[k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
      io    += n ;
      swapR += n ;
    }

    // primo blocco
    io = in_out ;
    trsv( LOWER, NO_TRANSPOSE, UNIT,
          row0,
          block0, row0, io, 1 ) ;

    io += row0 ;
    // blocchi intermedi
    nblk = 0 ;
    while ( nblk < nblock ) {

      /*
      //
      //  +----+---+---+---+
      //  |  M | L |   :   |
      //  |  M | L | L :   |
      //  +----+---+---+---+
      //  row00
      */

      valuePointer io1 = io - row00 ;
      valuePointer M   = this->AdAu_blk + nblk * nxnx2 ;
      valuePointer L   = M + row00 * n ;

      // io -= M*io1
      gemv( NO_TRANSPOSE,
            n, row00,
            -1, M, n,
            io1, 1,
            1, io, 1 ) ;

      trsv( LOWER, NO_TRANSPOSE, UNIT,
            n, L, n, io, 1 ) ;

      io += n ;
      ++nblk ;
    }

    // soluzione ultimo blocco
    integer ncol = colN-rowN ;
    gemv( NO_TRANSPOSE,
          rowN, ncol,
          -1, blockN, rowN,
          io-ncol, 1,
          1, io, 1 ) ;

    this->la_factorization->solve(io) ;

    while ( nblk > 0 ) {
      --nblk ;
      io -= n ;

      /*
      //  +---------+---+-------+
      //  |       U |   |       |
      //  |         | U |   M   |
      //  +---------+---+-------+
      //  row00
      */

      valuePointer io1 = io + n ;
      valuePointer U   = this->AdAu_blk + nblk * nxnx2 + row00 * n ;
      valuePointer M   = U + nxn ;

      gemv( NO_TRANSPOSE,
            n, n_m_row00,
            -1, M, n,
            io1, 1,
            1, io, 1 ) ;

      trsv( UPPER, NO_TRANSPOSE, NON_UNIT,
            n, U, n, io, 1 ) ;
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
    integer const * swapC = swapRC_blks+(nblock*n+col00) ;
    io = in_out + nblock*n + col00 ;
    integer k = row00 ;
    while ( k > 0 ) {
      integer k1 = swapC[--k] ; // 0 based
      if ( k1 > k ) std::swap( io[k], io[k1] ) ;
    } ;
    for ( nblk = 0 ; nblk < nblock ; ++nblk )  {
      io    -= n ;
      swapC -= n ;
      k      = row00 ;
      while ( k > 0 ) {
        integer k1 = swapC[--k] ; // 0 based
        if ( k1 > k ) std::swap( io[k], io[k1] ) ;
      } ;
    }

    // permuto le x
    std::rotate( in_out, in_out + col00, in_out + neq ) ;
  }
 
  template class DiazLU<double> ;
  template class DiazLU<float> ;

}

///
/// eof: DiazLU.cc
///

