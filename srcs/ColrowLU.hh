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
/// file: ColrowLU.hh
///

#ifndef COLROW_LU_HH
#define COLROW_LU_HH

#include "alglin.hh"
#include <iostream>

namespace alglin {

  template <typename t_Value>
  void
  print_colrow( std::ostream & stream,
                integer         numBlock,
                integer         dimBlock,
                integer         row0,
                integer         rowN,
                t_Value const * block0,
                t_Value const * blocks,
                t_Value const * blockN ) ;

  template <typename t_Value>
  void
  print_colrow_to_maple( std::ostream & stream,
                         integer         numBlock,
                         integer         dimBlock,
                         integer         row0,
                         integer         rowN,
                         t_Value const * block0,
                         t_Value const * blocks,
                         t_Value const * blockN ) ;

  //! LU decomposition of a ABD matrix
  /*!
   * 
   * \date     May, 2016
   * \version  1.0
   * \note     first release May, 2016
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Dipartimento di Ingegneria Industriale <br>
   *           University of Trento <br>
   *           via Mesiano 77, I -- 38050 Trento, Italy <br>
   *           enrico.bertolazzi@ing.unitn.it
   * 
   */
  template <typename t_Value>
  class ColrowLU {
  public:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

  private:

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseIndex ;

    ColrowLU( ColrowLU const & ) ;
    ColrowLU const & operator = ( ColrowLU const & ) ;
    
    static valueType const epsi ;
    
    integer      numBlock ;
    integer      dimBlock ;
    integer      dimBlock_m_row0 ;
    integer      Nlast ;
    integer      row0 ;
    integer      rowN ;
    valuePointer block0 ;
    valuePointer blocks ;
    valuePointer blockNN ;
    integer *    swapRC_blks ;
    bool full_LU ;

    // internal parameters
    mutable valuePointer Amat ;
    mutable integer      ldA ;
    mutable valuePointer Bmat ;
    mutable integer *    swapRC ;
    mutable integer      sizeBlock ;
    mutable integer      offs ;
    mutable integer      nblk ;

    void factorize() ;
    void factorize_block() ;
    void solver_block_L( valuePointer in_out ) const ;
    void solver_block_U( valuePointer in_out ) const ;

  public:

    explicit ColrowLU() ;
    ~ColrowLU() ;

    //! compute y = alpha*A*x+beta*y
    static
    void
    mv( integer           _numBlock,
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
        integer           incy ) ;

    //! factorize the matrix
    void
    factorize( integer           _numBlock,
               integer           _dimBlock,
               integer           _row0,
               integer           _rowN,
               valueConstPointer _block0,
               valueConstPointer _blocks,
               valueConstPointer _blockN ) ;

    //! factorize the matrix
    void
    factorize_inplace( integer      _numBlock,
                       integer      _dimBlock,
                       integer      _row0,
                       integer      _rowN,
                       valuePointer _block0,
                       valuePointer _blocks,
                       valuePointer _blockN ) ;

    //! solve linear sistem using internal factorized matrix
    void solve( valuePointer in_out ) const ;

    void print( std::ostream & stream ) const ;

  } ;
}

#endif

///
/// eof: ColrowLU.hh
///

