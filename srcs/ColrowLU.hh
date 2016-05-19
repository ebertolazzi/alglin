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

// Eigen3
#include <Eigen/Dense>

namespace alglin {

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

    typedef enum { ColrowLU_QR0 = 0,
                   ColrowLU_QR1 = 1,
                   ColrowLU_QR2 = 2,
                   ColrowLU_LU0 = 3,
                   ColrowLU_LU1 = 4 } LAST_BLOCK ;

    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,1>              vec ;
    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,Eigen::Dynamic> mat ;

    Eigen::HouseholderQR<mat>        la_solve0 ;
    Eigen::ColPivHouseholderQR<mat>  la_solve1 ;
    Eigen::FullPivHouseholderQR<mat> la_solve2 ;
    Eigen::PartialPivLU<mat>         la_solve3 ;
    Eigen::FullPivLU<mat>            la_solve4 ;

    static
    void
    print_colrow( std::ostream & stream,
                  integer         numBlock,
                  integer         dimBlock,
                  integer         row0,
                  integer         rowN,
                  t_Value const * block0,
                  t_Value const * blocks,
                  t_Value const * blockN ) ;

    static
    void
    print_colrow_to_maple( std::ostream & stream,
                           integer         numBlock,
                           integer         dimBlock,
                           integer         row0,
                           integer         rowN,
                           t_Value const * block0,
                           t_Value const * blocks,
                           t_Value const * blockN ) ;

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
    LAST_BLOCK   last_block ;

    // internal parameters
    mutable valuePointer Amat ;
    mutable integer      ldA ;
    mutable valuePointer Bmat ;
    mutable integer *    swapRC ;
    mutable integer      sizeBlock ;
    mutable integer      offs ;
    mutable integer      nblk ;

    void factorize() ;
    void factorize_block( bool first ) ;
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

    //! compute r = b-A*x
    static
    void
    residue( integer           _numBlock,
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
             integer           incr ) ;

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

