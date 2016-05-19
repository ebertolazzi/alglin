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
/// file: Colrow.hh
///

#ifndef COLROW_HH
#define COLROW_HH

#include <iostream>

// Eigen3
#include <Eigen/Dense>

namespace alglin {

  typedef int integer ;

  template <typename t_Value>
  void
  print_colrow( std::ostream &  stream,
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
  class Colrow {
  public:

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    typedef enum { Colrow_QR0 = 0,
                   Colrow_QR1 = 1,
                   Colrow_QR2 = 2,
                   Colrow_LU0 = 3,
                   Colrow_LU1 = 4 } LAST_BLOCK ;

    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,1>              vec ;
    typedef Eigen::Matrix<t_Value,Eigen::Dynamic,Eigen::Dynamic> mat ;

    typedef Eigen::Matrix<integer,Eigen::Dynamic,1>              ivec ;
    typedef Eigen::Matrix<integer,Eigen::Dynamic,Eigen::Dynamic> imat ;

    Eigen::HouseholderQR<mat>        la_solve0 ;
    Eigen::ColPivHouseholderQR<mat>  la_solve1 ;
    Eigen::FullPivHouseholderQR<mat> la_solve2 ;
    Eigen::PartialPivLU<mat>         la_solve3 ;
    Eigen::FullPivLU<mat>            la_solve4 ;

  private:

    Colrow( Colrow const & ) ;
    Colrow const & operator = ( Colrow const & ) ;

    static valueType const epsi ;
    
    integer    numBlock ;
    integer    dimBlock ;
    integer    dimBlock_m_row0 ;
    integer    Nlast ;
    integer    row0 ;
    integer    rowN ;
    mat        block0 ;
    mat        blocks ;
    mat        blockNN ;
    ivec       swapRC_blks ;
    LAST_BLOCK last_block ;

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
    void solver_block_L( vec & in_out ) const ;
    void solver_block_U( vec & in_out ) const ;

  public:

    explicit Colrow() ;
    ~Colrow() ;

    //! compute y = alpha*A*x+beta*y
    static
    void
    mv( integer     _numBlock,
        integer     _dimBlock,
        integer     _row0,
        integer     _rowN,
        mat const & _block0,
        mat const & _blocks,
        mat const & _blockN,
        valueType   alpha,
        vec const & x,
        valueType   beta,
        vec       & y ) ;

    //! factorize the matrix
    void
    factorize( integer     _numBlock,
               integer     _dimBlock,
               integer     _row0,
               integer     _rowN,
               mat const & _block0,
               mat const & _blocks,
               mat const & _blockN ) ;

    //! solve linear sistem using internal factorized matrix
    void solve( vec & in_out ) const ;

    void print( std::ostream & stream ) const ;

  } ;
}

#endif

///
/// eof: Colrow.hh
///
