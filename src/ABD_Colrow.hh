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
/// file: ABD_Colrow.hh
///

#ifndef ABD_COLROW_HH
#define ABD_COLROW_HH

#include "Alglin.hh"
#include "Alglin++.hh"
#include <iostream>

namespace alglin {

  //! available LU factorization code
  typedef enum {
    COLROW_LASTBLOCK_LU  = 0,
    COLROW_LASTBLOCK_QR  = 1,
    COLROW_LASTBLOCK_SVD = 2
  } COLROW_LASTBLOCK_Choice;

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
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           `enrico.bertolazzi@unitn.it`
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

    LU<t_Value>  la_lu ;
    QR<t_Value>  la_qr ;
    SVD<t_Value> la_svd ;

    COLROW_LASTBLOCK_Choice last_block ;

    ColrowLU( ColrowLU const & ) ;
    ColrowLU const & operator = ( ColrowLU const & ) ;

    static valueType const epsi ;

    integer      neq ;
    integer      numBlock ;
    integer      dimBlock ;
    integer      Nlast ;
    integer      row0, col0 ;
    integer      rowN, colN ;

    integer      sizeBlock ;
    integer      hSizeBlock ;
    integer      row00, col00 ;
    integer      rowNN, colNN ;
    integer      dimBlock_m_row00 ;
    integer      NB ;

    mutable integer nblk ;

    valuePointer block0 ;
    valuePointer blocks ;
    valuePointer blockN ;
    integer *    swapRC_blks ;

    void setup() ;
    void factorize() ;

    void
    LU_left_right( integer nrA,
                   integer ncA,
                   integer ncL,
                   integer ncR,
                   t_Value * A, integer ldA,
                   integer swapR[] ) ;

    void
    LU_top_bottom( integer nrT,
                   integer nrA,
                   integer ncA,
                   t_Value * A, integer ldA,
                   integer nrB,
                   t_Value * B, integer ldB,
                   integer swapC[] ) ;

    integer *    matrixStructure ; //!< structure of the matrix
    integer *    pivot           ; //!< permutation array
    valuePointer array           ; //!< the matrix data

  public:

    explicit ColrowLU() ;
    ~ColrowLU() ;

    //! factorize the matrix
    void
    factorize( COLROW_LASTBLOCK_Choice choice,
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
               valueConstPointer _blockN ) ;

    //! factorize the matrix
    void
    factorize_inplace( COLROW_LASTBLOCK_Choice choice,
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
                       valuePointer _blockN ) ;

    //! solve linear sistem using internal factorized matrix
    void solve( valuePointer in_out ) const ;

    void print( std::ostream & stream ) const ;

  } ;
}

#endif

///
/// eof: ABD_Colrow.hh
///

