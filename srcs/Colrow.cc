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
/// file: Colrow.cc
///

#include "Colrow.hh"
#include <iomanip>
#include <vector>
#include <limits>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <stdexcept>

#ifndef ALGLIN_ASSERT
  #define ALGLIN_ASSERT(COND,MSG) if ( !(COND) ) { \
    std::ostringstream ost ; ost << "in alglin::" << MSG << '\n' ; \
    throw std::runtime_error(ost.str()) ; \
  }
#endif


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

  template <> float  const Colrow<float>::epsi  = 100*std::numeric_limits<float>::epsilon() ;
  template <> double const Colrow<double>::epsi = 1000*std::numeric_limits<double>::epsilon() ;

  template <typename t_Value>
  Colrow<t_Value>::Colrow()
  //, last_block(ColrowLU_LU)
  //, last_block(ColrowLU_fullLU)
  : last_block(Colrow_QR0)
  {
  }
  
  template <typename t_Value>
  Colrow<t_Value>::~Colrow() {
  }

  //! compute y += alpha*A*x
  template <typename t_Value>
  void
  Colrow<t_Value>::mv( integer     _numBlock,
                       integer     _dimBlock,
                       integer     _row0,
                       integer     _rowN,
                       mat const & _block0,
                       mat const & _blocks,
                       mat const & _blockN,
                       valueType   alpha,
                       vec const & x,
                       valueType   beta,
                       vec       & y ) {

    // first block y = alpha * _block0 * x + beta * y
    y *= beta ;
    y.head(_row0) += alpha * _block0 * x.head(_dimBlock) ;
    
    // internal blocks block
    for ( integer i = 0 ; i < _numBlock ; ++i )
      y.segment(_row0+i*_dimBlock,_dimBlock) +=
        alpha * _blocks.block(i*_dimBlock,0,_dimBlock,2*_dimBlock) * x.segment(i*_dimBlock,2*_dimBlock) ;

    // last block
    y.tail(_rowN) += alpha * _blockN * x.tail(_rowN+_row0) ;
  }

  template <typename t_Value>
  void
  Colrow<t_Value>::factorize( integer     _numBlock,
                              integer     _dimBlock,
                              integer     _row0,
                              integer     _rowN,
                              mat const & _block0,
                              mat const & _blocks,
                              mat const & _blockN ) {
    numBlock        = _numBlock ;
    dimBlock        = _dimBlock ;
    dimBlock_m_row0 = _dimBlock-_row0 ;
    row0            = _row0 ;
    rowN            = _rowN ;
    sizeBlock       = 2*dimBlock*dimBlock ;
    offs            = dimBlock*dimBlock - dimBlock_m_row0 ;
    Nlast           = row0+rowN ;

    // allocate
    block0      = _block0 ;
    blocks      = _blocks ;
    blockNN.resize(Nlast,Nlast) ;
    swapRC_blks.resize( numBlock*dimBlock+2*Nlast ) ;
    
    blockNN.block( 0, 0, row0, Nlast ).setZero() ;
    blockNN.block( row0, 0, Nlast-row0, Nlast ) = _blockN ;

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
  Colrow<t_Value>::factorize_block() {
    // applico row0 passi di Gauss con pivoting di colonna
    valuePointer Bk  = Bmat ;
    valuePointer Ak  = Amat ;
    valuePointer Akk = Amat ;
    integer k = 0 ;
    bool do_loop ;
    do {
      // cerco pivot
      integer ipiv = iamax( dimBlock-k, Akk, ldA ) ;
      // scambio colonne
      if ( ipiv > 0 ) {
        swap( row0,     Ak, 1, Ak + ipiv*ldA,      1 ) ;
        swap( dimBlock, Bk, 1, Bk + ipiv*dimBlock, 1 ) ;
      }
      swapRC[k] = k+ipiv ; // memorizzo scambio
      // controllo pivot non zero
      ALGLIN_ASSERT( std::abs(Akk[0]) > epsi,
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
    // Apply transform to block B
    trsm( SideMultiply::RIGHT,
          ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          dimBlock, row0, 1.0,
          Amat, ldA,
          Bmat, dimBlock ) ;
    Bk = Bmat+dimBlock*row0 ;
    gemm( Transposition::NO_TRANSPOSE,
          Transposition::NO_TRANSPOSE,
          dimBlock, dimBlock_m_row0, row0,
          -1.0,
          Bmat, dimBlock,
          Amat+ldA*row0, ldA,
          1.0,
          Bk, dimBlock ) ;
    // ----------
    integer ierr = getrf( dimBlock, dimBlock_m_row0, Bk, dimBlock, swapRC+row0 ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_block, found ierr = " << ierr <<
                   " at LU factorizatioon on block N." << nblk+1 ) ;
    valuePointer Cmat = Bmat+dimBlock*dimBlock ;
    ierr = swaps( dimBlock, Cmat, dimBlock,
                  0, dimBlock_m_row0-1,
                  swapRC+row0, 1 ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize_block, found ierr = " << ierr <<
                   " at LU swaps rows on block N." << nblk+1 ) ;
    trsm( SideMultiply::LEFT,
          ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          dimBlock_m_row0, dimBlock, 1.0,
          Bk, dimBlock,
          Cmat, dimBlock ) ;
    gemm( Transposition::NO_TRANSPOSE,
          Transposition::NO_TRANSPOSE,
          row0, dimBlock, dimBlock_m_row0,
          -1.0,
          Bk+dimBlock_m_row0, dimBlock,
          Cmat, dimBlock,
          1.0,
          Cmat+dimBlock_m_row0, dimBlock ) ;
  }

  template <typename t_Value>
  void
  ColrowLU<t_Value>::factorize() {
    // primo blocco
    swapRC = swapRC_blks ;
    Amat   = block0 ;
    ldA    = row0 ;
    Bmat   = blocks ;
    nblk   = 0 ;
    factorize_block() ;

    // blocchi intermedi (n-1)
    ldA = dimBlock ;
    while ( ++nblk < numBlock ) {
      swapRC += dimBlock ;
      Bmat   += sizeBlock ;
      Amat    = Bmat-offs ;
      factorize_block() ;
    }

    swapRC += dimBlock ;
    Amat    = Bmat + sizeBlock - offs ;

    // copio in ultimo blocco
    integer N = row0+rowN ;
    integer ierr = gecopy( row0, dimBlock, Amat, dimBlock, blockNN, N ) ;
    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::factorize, in gecopy ierr = " << ierr ) ;

    // fattorizzazione ultimo blocco
    switch ( last_block ) {
      case ColrowLU_LU:
        ierr = getrf( N, N, blockNN, N, swapRC ) ;
        break ;
      case ColrowLU_fullLU:
        ierr = getc2( N, blockNN, N, swapRC, swapRC+N ) ;
        break ;
      case ColrowLU_QR:
        //#define USE_EIGEN
        #ifndef USE_EIGEN
        //ierr = geqp3( N, N, blockNN, N, swapRC, Tau, Work, Lwork ) ;
        ierr = geqr2( N, N, blockNN, N, Tau, Work ) ;
        //ierr = geqrf( N, N, blockNN, N, Tau, Work ) ;
        #else
        qr.compute(Eigen::Map<mat>(blockNN,N,N)) ;
        #endif
        break ;
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
  ColrowLU<t_Value>::solver_block_L( valuePointer in_out ) const {
    valuePointer io1 = in_out+row0 ;

    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          row0,
          Amat, ldA, in_out, 1 ) ;
    gemv( Transposition::NO_TRANSPOSE,
          dimBlock, row0,
          -1, Bmat, dimBlock,
          in_out, 1,
          1, io1, 1 ) ;
    valuePointer io2 = in_out+dimBlock ;
    valuePointer Bk  = Bmat+dimBlock*row0 ;

    // applico permutazione
    for ( integer k = 0 ; k < dimBlock_m_row0 ; ++k ) {
      integer k1 = swapRC[row0+k]-1 ; // base 1 da lapack
      if ( k1 > k ) std::swap( io1[k], io1[k1] ) ;
    }

    trsv( ULselect::LOWER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::UNIT,
          dimBlock_m_row0,
          Bk, dimBlock, io1, 1 ) ;
    gemv( Transposition::NO_TRANSPOSE,
          row0, dimBlock_m_row0,
          -1, Bk+dimBlock_m_row0, dimBlock,
          io1, 1,
          1, io2, 1 ) ;
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
  ColrowLU<t_Value>::solver_block_U( valuePointer in_out ) const {
    // y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
    valuePointer io1  = in_out+row0 ;
    valuePointer io2  = in_out+dimBlock ;
    valuePointer Cmat = Bmat+dimBlock*dimBlock ;
    gemv( Transposition::NO_TRANSPOSE,
          dimBlock_m_row0, dimBlock,
          -1, Cmat, dimBlock,
          io2, 1,
          1, io1, 1 ) ;
    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          dimBlock_m_row0,
          Bmat+row0*dimBlock, dimBlock, io1, 1 ) ;
    gemv( Transposition::NO_TRANSPOSE,
          row0, dimBlock_m_row0,
          -1, Amat+row0*ldA, ldA,
          io1, 1,
          1, in_out, 1 ) ;
    trsv( ULselect::UPPER,
          Transposition::NO_TRANSPOSE,
          DiagonalType::NON_UNIT,
          row0,
          Amat, ldA, in_out, 1 ) ;

    // applico permutazione
    for ( integer k = row0-1 ; k >= 0 ; --k ) {
      integer k1 = swapRC[k] ; // 0 based
      if ( k1 > k ) std::swap( in_out[k], in_out[k1] ) ;
    }
  }

  template <typename t_Value>
  void
  ColrowLU<t_Value>::solve( valuePointer in_out ) const {

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
      swapRC += dimBlock ;
      Bmat   += sizeBlock ;
      Amat    = Bmat-offs ;
      io     += dimBlock ;
      solver_block_L( io ) ;
    }

    // risolve rispetto ultimo blocco
    integer N = row0+rowN ;
    swapRC += dimBlock ;
    Bmat   += sizeBlock ;
    io     += dimBlock ;

    integer ierr = 0 ;
    switch ( last_block ) {
      case ColrowLU_LU:
        ierr = getrs( Transposition::NO_TRANSPOSE, N, 1, blockNN, N, swapRC, io, N ) ;
        break ;
      case ColrowLU_fullLU:
        {
          valueType scale = gesc2( N, blockNN, N, io, swapRC, swapRC+N ) ;
          scal( N, scale, io, 1 ) ;
        }
        break ;
      case ColrowLU_QR:
        #ifndef USE_EIGEN
        ierr = ormqr( SideMultiply::LEFT,
                      Transposition::TRANSPOSE,
                      N, 1, N, blockNN, N,
                      Tau, io, N, Work, Lwork ) ;
        trsv( ULselect::UPPER,
              Transposition::NO_TRANSPOSE,
              DiagonalType::NON_UNIT,
              N, blockNN, N, io, 1 ) ;
        // applico permutazione
        //for ( integer i = 0 ; i < N ; ++i ) Work[swapRC[i]-1] = io[i] ;
        //copy( N, Work, 1, io, 1 ) ;
        #else
        Eigen::Map<vec>(io,N) = qr.solve(Eigen::Map<vec>(io,N)) ;
        #endif
        break ;
    }

    ALGLIN_ASSERT( ierr == 0,
                   "ColrowLU::solve, failed to solve on last block ierr = " << ierr ) ;

    for ( integer i = 1 ; i < numBlock ; ++i ) {
      swapRC -= dimBlock ;
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
  }

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


  template <typename t_Value>
  void
  print_colrow( std::ostream & stream,
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

  template <typename t_Value>
  void
  print_colrow_to_maple( std::ostream & stream,
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

  template void print_colrow( std::ostream & stream,
                              integer         numBlock,
                              integer         dimBlock,
                              integer         row0,
                              integer         rowN,
                              double const *  block0,
                              double const *  blocks,
                              double const *  blockN ) ;

  template void print_colrow( std::ostream & stream,
                              integer         numBlock,
                              integer         dimBlock,
                              integer         row0,
                              integer         rowN,
                              float const *   block0,
                              float const *   blocks,
                              float const *   blockN ) ;

  template void print_colrow_to_maple( std::ostream & stream,
                                       integer         numBlock,
                                       integer         dimBlock,
                                       integer         row0,
                                       integer         rowN,
                                       double const *  block0,
                                       double const *  blocks,
                                       double const *  blockN ) ;

  template void print_colrow_to_maple( std::ostream & stream,
                                       integer         numBlock,
                                       integer         dimBlock,
                                       integer         row0,
                                       integer         rowN,
                                       float const *   block0,
                                       float const *   blocks,
                                       float const *   blockN ) ;

  template class ColrowLU<double> ;
  template class ColrowLU<float> ;

}

///
/// eof: ColrowLU.cc
///

