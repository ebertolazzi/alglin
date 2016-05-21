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

#include "Alglin.hh"
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
  /*
      Matrix NNZ structure
          col0
       |       |
     / +-------+                         \
     | |  TOP  |                         | <-- row0
     | +-------+----+                    |
     |    |    |    |                    | <- sizeBlock
     |    +----+----+----+               |
     |         |    |    |               |
     |         +----+----+----+          |
     |              |    |    |          |
     |              +----+----+----+     |
     |                   |    |    |     |
     |                   +----+----+---+ |
     |                        |        | | <-- rowN
     |                        | BOTTOM | |
     \                        +--------+ /
                              |        |
                                 colN
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
    print( std::ostream & stream,
           integer         row0,
           integer         col0,
           t_Value const * block0,
           integer         numBlock,
           integer         dimBlock,
           t_Value const * blocks,
           integer         rowN,
           integer         colN,
           t_Value const * blockN ) ;

  private:

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseIndex ;

    ColrowLU( ColrowLU const & ) ;
    ColrowLU const & operator = ( ColrowLU const & ) ;
    
    static valueType const epsi ;
    
    integer      neq ;
    integer      numBlock ;
    integer      dimBlock ;
    integer      Nlast ;
    integer      row0, col0 ;
    integer      rowN, colN ;
    valuePointer block0 ;
    valuePointer blocks ;
    valuePointer blockN ;
    integer *    swapRC_blks ;
    LAST_BLOCK   last_block ;

    integer      sizeBlock ;
    integer      hSizeBlock ;
    integer      row00, col00 ;
    integer      rowNN, colNN ;
    integer      dimBlock_m_row00 ;
    
    mutable integer nblk ;

    void setup() ;

    void factorize() ;
    void factorize_first_block() ;
    void factorize_block() ;
    void factorize_last_block() ;

    void solve_block_L( valuePointer in_out ) const ;
    void solve_block_U( valuePointer in_out ) const ;
    void solve_last_block( valuePointer in_out ) const ;
    void LU_full( integer nr, integer nc,
                  t_Value * A, integer ldA,
                  integer swapR[], integer swapC[] ) ;

  public:

    explicit ColrowLU() ;
    ~ColrowLU() ;

    //! compute y = alpha*A*x+beta*y
    static
    void
    mv( integer           _row0,
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
        integer           incy ) ;

    //! compute r = b-A*x
    static
    void
    residue( integer           _row0,
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
             integer           incr ) ;

    //! factorize the matrix
    void
    factorize( integer           _row0,
               integer           _col0,
               valueConstPointer _block0,
               integer           _numBlock,
               integer           _dimBlock,
               valueConstPointer _blocks,
               integer           _rowN,
               integer           _colN,
               valueConstPointer _blockN ) ;

    //! factorize the matrix
    void
    factorize_inplace( integer      _row0,
                       integer      _col0,
                       valuePointer _block0,
                       integer      _numBlock,
                       integer      _dimBlock,
                       valuePointer _blocks,
                       integer      _rowN,
                       integer      _colN,
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

