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

#include "BlockBidiagonal.hh"

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace alglin {

  string
  LastBlock_to_string( LASTBLOCK_Choice c ) {
    switch ( c ) {
      case LASTBLOCK_LU:  return "last block LU"  ;
      case LASTBLOCK_QR:  return "last block QR"  ;
      case LASTBLOCK_QRP: return "last block QRP" ;
      case LASTBLOCK_SVD: return "last block SVD" ;
    }
    return "last block not selected" ;
  }
  
  /*\
   |   _                 _ ____        _   _
   |  | | ___   __ _  __| | __ )  ___ | |_| |_ ___  _ __ ___
   |  | |/ _ \ / _` |/ _` |  _ \ / _ \| __| __/ _ \| '_ ` _ \
   |  | | (_) | (_| | (_| | |_) | (_) | |_| || (_) | | | | | |
   |  |_|\___/ \__,_|\__,_|____/ \___/ \__|\__\___/|_| |_| |_|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::loadBottom(
    valueConstPointer H0, integer ld0,
    valueConstPointer HN, integer ldN,
    valueConstPointer Hq, integer ldQ
  ) {

    numInitialBc    = 0 ;
    numFinalBc      = 0 ;
    numCyclicBC     = n+q ;
    numInitialOMEGA = 0 ;
    numFinalOMEGA   = 0 ;
    numCyclicOMEGA  = q ;

    integer m = n + q ;
    gecopy( m, n, H0, ld0, H0Nq,       m ) ;
    gecopy( m, n, HN, ldN, H0Nq+m*n,   m ) ;
    gecopy( m, q, Hq, ldQ, H0Nq+2*m*n, m ) ;

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
  BlockBidiagonal<t_Value>::loadTopBottom(
    integer           row0,
    integer           col0,
    valueConstPointer block0_in,
    integer           ld0,
    // ----------------------------
    integer           rowN,
    integer           colN,
    valueConstPointer blockN_in,
    integer           ldN
  ) {

    numInitialBc    = row0 ;
    numFinalBc      = rowN ;
    numCyclicBC     = 0 ;
    numInitialOMEGA = col0 - n ;
    numFinalOMEGA   = colN - n ;
    numCyclicOMEGA  = 0 ;

    block0 = H0Nq ;
    blockN = H0Nq+row0*col0 ;

    gecopy( row0, col0, block0_in, ld0, block0, row0 ) ;
    gecopy( rowN, colN, blockN_in, ldN, blockN, rowN ) ;

  }

  /*\
   |   _                 _ ____   ____
   |  | | ___   __ _  __| | __ ) / ___|
   |  | |/ _ \ / _` |/ _` |  _ \| |
   |  | | (_) | (_| | (_| | |_) | |___
   |  |_|\___/ \__,_|\__,_|____/ \____|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::loadBC(
    integer _numInitialBc,
    integer _numFinalBc,
    integer _numCyclicBC,
    // ----------------------
    integer _numInitialOMEGA,
    integer _numFinalOMEGA,
    integer _numCyclicOMEGA,
    // ----------------------
    valueConstPointer H0, integer ld0,
    valueConstPointer HN, integer ldN,
    valueConstPointer Hq, integer ldQ
  ) {

    integer nBC = _numInitialBc + _numFinalBc + _numCyclicBC ;
    integer nOM = _numInitialOMEGA + _numFinalOMEGA + _numCyclicOMEGA ;

    ALGLIN_ASSERT( n+q == nBC && nOM == q,
                   "BlockBidiagonal::loadBC bad parameter(s):" <<
                   "\nn               = " << n <<
                   "\nq               = " << q <<
                   "\nnumInitialBc    = " << _numInitialBc    <<
                   "\nnumFinalBc      = " << _numFinalBc      <<
                   "\nnumCyclicBC     = " << _numCyclicBC     <<
                   "\nnumInitialOMEGA = " << _numInitialOMEGA <<
                   "\nnumFinalOMEGA   = " << _numFinalOMEGA   <<
                   "\nnumCyclicOMEGA  = " << _numCyclicOMEGA ) ;

    numInitialBc    = _numInitialBc ;
    numFinalBc      = _numFinalBc ;
    numCyclicBC     = _numCyclicBC ;
    numInitialOMEGA = _numInitialOMEGA ;
    numFinalOMEGA   = _numFinalOMEGA ;
    numCyclicOMEGA  = _numCyclicOMEGA ;
    
    if ( numCyclicBC == 0 && numCyclicOMEGA == 0 ) {
      /*\
       |  +----+-----+---+
       |  | H0 | HN  |Hq |
       |  |    |     |   |
       |  +----+-----+---+
       |  +----+------+---------+
       |  | 0  | blkN |blkN: 0  |
       |  |blk0|  0   |    :blk0|
       |  +----+------+---------+
      \*/
      integer row0  = numInitialBc ;
      integer rowN  = numFinalBc ;
      integer col00 = numInitialOMEGA ;
      integer colNN = numFinalOMEGA ;

      block0 = H0Nq ;
      blockN = H0Nq+row0*(n+col00) ;

      gecopy( rowN, n,     HN, ldN, blockN,        rowN ) ;
      gecopy( rowN, colNN, Hq, ldQ, blockN+n*rowN, rowN ) ;

      gecopy( row0, col00, Hq+rowN+colNN*ldQ, ldQ, block0,            row0 ) ;
      gecopy( row0, n,     H0+rowN,           ld0, block0+col00*row0, row0 ) ;

    } else {
      integer m = n + q ;
      gecopy( m, n, H0, ld0, H0Nq,       m ) ;
      gecopy( m, n, HN, ldN, H0Nq+m*n,   m ) ;
      gecopy( m, q, Hq, ldQ, H0Nq+2*m*n, m ) ;
    }
  }
  
  /*\
   |   _           _     _     _            _
   |  | | __ _ ___| |_  | |__ | | ___   ___| | __
   |  | |/ _` / __| __| | '_ \| |/ _ \ / __| |/ /
   |  | | (_| \__ \ |_  | |_) | | (_) | (__|   <
   |  |_|\__,_|___/\__| |_.__/|_|\___/ \___|_|\_\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::last_block_factorize() {
    /*
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!! factorization of the last block !!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    */
    /*
    // / S  R  0  \ /x(0)\  = b(0)
    // \ H0 HN Hq / \x(N)/  = b(N)
    */
    integer m  = n+q ;
    integer mn = m+n ;

    la_factorization->allocate( mn, mn ) ;
    if ( this->numCyclicBC == 0 && this->numCyclicOMEGA == 0 ) {
      /*\
       |  +----+-----+---+
       |  | H0 | HN  |Hq |
       |  |    |     |   |
       |  +----+-----+---+
       |  +----+------+---------+
       |  | 0  | blkN |blkN: 0  |
       |  |blk0|  0   |    :blk0|
       |  +----+------+---------+
      \*/
      integer row0  = numInitialBc ;
      integer rowN  = numFinalBc ;
      integer col00 = numInitialOMEGA ;
      integer colNN = numFinalOMEGA ;
      la_factorization->load_block( n, nx2, AdAu_blk, n ) ;
      la_factorization->zero_block( m, mn, n, 0 ) ;
      la_factorization->load_block( rowN, n+colNN, blockN, rowN, n, n ) ;
      la_factorization->load_block( row0, n,     block0+col00*row0, row0, n+rowN, 0 ) ;
      la_factorization->load_block( row0, col00, block0,            row0, n+rowN, nx2+colNN ) ;
    } else {
      la_factorization->load_block( n, nx2, AdAu_blk, n ) ;
      la_factorization->load_block( m, mn, H0Nq, m, n, 0 ) ;
      if ( m > n ) la_factorization->zero_block( n, q, 0, nx2 ) ;
    }

    // fattorizzazione ultimo blocco
    this->la_factorization->factorize() ;
  }

  template <typename T>
  static
  void
  dumpOneMatrix ( basic_ostream<char> & stream,
                  char const *          name,
                  T const               M[],
                  integer               numRow,
                  integer               numCol ) {
    stream << "# " << name << " Size: " << numRow << " x " << numCol << '\n' ;
    stream << name << " := <" ;
    for ( integer nc = 0 ; nc < numCol ; ++nc ) {
      stream << '<' ;
      for ( integer nr = 0 ; nr < numRow ; ++nr ) {
        stream << M[ nr + nc * numRow ] ;
        if ( nr+1 < numRow ) stream << ',' ;
        else                 stream << '>' ;
      }
      if ( nc+1 < numCol ) stream << "|\n" ;
      else                 stream << "> ;\n" ;
    }
  }

  template <typename t_Value>
  void
  BlockBidiagonal<t_Value>::dumpMatrix( basic_ostream<char> & stream ) const {

    stream << "interface( rtablesize = 40 ) ;\n" ;
    for ( integer row = 0 ; row < nblock ; ++row ) {
      valueConstPointer Ad = AdAu_blk + 2*row*n*n ;
      valueConstPointer Au = Ad + n*n ;
      dumpOneMatrix( stream, "Ad", Ad, n, n ) ;
      dumpOneMatrix( stream, "Au", Au, n, n ) ;
    }

    dumpOneMatrix( stream, "H0Nq", H0Nq, n+q, 2*n+q ) ;

    stream << "with(LinearAlgebra):\n" ;
    stream << "Determinant(Ad);\n" ;
    stream << "Determinant(Au);\n" ;
    stream << "Rank(H0Np);\n" ;
    stream << "Rank(<H0N|Hp>);\n" ;

  }

  template class BlockBidiagonal<float> ;
  template class BlockBidiagonal<double> ;

}
