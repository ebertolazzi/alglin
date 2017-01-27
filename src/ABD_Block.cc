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
/// file: ABD_Block.cc
///

#include "ABD_Block.hh"
#include <iomanip>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

#if 0

namespace alglin {

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/

  template <typename t_Value>
  void
  BlockLU<t_Value>::factorize() {

    ALGLIN_ASSERT( this->numCyclicOMEGA == 0 && this->numCyclicBC == 0,
                   "DiazLU cannot manage cyclic BC" ) ;

    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;

    integer const & col00  = this->numInitialOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const col0      = n + col00 ;
    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valuePointer & block0 = this-> block0 ;
    valuePointer & blockN = this-> blockN ;

  }

  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BlockLU<t_Value>::solve( valuePointer in_out ) const {
  
    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;
    integer const & col00  = this->numInitialOMEGA ;
    integer const & colNN  = this->numFinalOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const col0      = n + col00 ;
    integer const colN      = n + colNN ;
    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valueConstPointer const & block0 = this->block0 ;
    valueConstPointer const & blockN = this->blockN ;

    integer neq = nblock*n+row0+rowN ;
    std::rotate( in_out, in_out + neq - row0, in_out + neq ) ;


    // permuto le x
    std::rotate( in_out, in_out + col00, in_out + neq ) ;
  }
 
 
  // ---------------------------------------------------------------------------

  template <typename t_Value>
  void
  BlockLU<t_Value>::solve( integer      nrhs,
                          valuePointer in_out,
                          integer      ldRhs ) const {
  
    integer const & n      = this->n ;
    integer const & nxnx2  = this->nxnx2 ;
    integer const & nxn    = this->nxn      ;
    integer const & nblock = this->nblock ;
    integer const & col00  = this->numInitialOMEGA ;
    integer const & colNN  = this->numFinalOMEGA ;
    integer const & row0   = this->numInitialBC ;
    integer const & rowN   = this->numFinalBC ;

    integer const col0      = n + col00 ;
    integer const colN      = n + colNN ;
    integer const row00     = row0 - col00 ;
    integer const n_m_row00 = n - row00 ;

    valueConstPointer const & block0 = this->block0 ;
    valueConstPointer const & blockN = this->blockN ;

    integer neq = nblock*n+row0+rowN ;

    // permuto le x
    valuePointer io = in_out ;
    for ( integer k = 0 ; k < nrhs ; ++k ) {
      std::rotate( io, io + neq - row0, io + neq ) ;
      io += ldRhs ;
    }


    // permuto le x
    io = in_out ;
    for ( integer k = 0 ; k < nrhs ; ++k ) {
      std::rotate( io, io + col00, io + neq ) ;
      io += ldRhs ;
    }
  }

  template class BlockLU<double> ;
  template class BlockLU<float> ;

}

#endif

///
/// eof: ABD_Block.cc
///

