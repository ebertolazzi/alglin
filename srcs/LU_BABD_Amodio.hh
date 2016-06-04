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

#ifndef LU_BABD_AMODIO_HH
#define LU_BABD_AMODIO_HH

#include "Alglin.hh"
#include <vector>

namespace alglin {

  /*
  //      _                        _ _       _    _   _ 
  //     / \   _ __ ___   ___   __| (_) ___ | |  | | | |
  //    / _ \ | '_ ` _ \ / _ \ / _` | |/ _ \| |  | | | |
  //   / ___ \| | | | | | (_) | (_| | | (_) | |__| |_| |
  //  /_/   \_\_| |_| |_|\___/ \__,_|_|\___/|_____\___/ 
  */

  //! LU decomposition of a BABD matrix
  /*!
   * 
   * \date     May 30, 2006
   * \version  1.0
   * \note     first release May 30, 2006
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Mechanics and Structures Engineering <br>
   *             University of Trento <br>
   *           via Mesiano 77, I -- 38050 Trento, Italy <br>
   *           enrico.bertolazzi@ing.unitn.it
   *
   */
  template <typename t_Value>
  class AmodioLU {
  private:
  
    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;
    std::vector<std::vector<bool> > LU_rows_blk ;

    AmodioLU(AmodioLU const &) ;
    AmodioLU const & operator = (AmodioLU const &) ;

    integer nblock ; //!< total number of blocks
    integer n      ; //!< size of square blocks
    integer m      ; //!< number final rows (m>=n)
    integer nnz    ; //!< total number of non zeros

    // some derived constanst
    integer nx2 ;
    integer nxn ;
    integer nxnx2 ;
    integer nm ;

    /*
    //
    //  Matrix structure
    //
    //                  n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n                        n     q
    //  +-----+-----+-----+----.........-----+-----+-----+    \
    //  |  Ad | Au  |  0  |                  |  0  |  0  | n   |
    //  +-----+-----+-----+             -----+-----+-----+     |
    //  |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |
    //  |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
    //  +-----+-----+-----+-----+       -----+-----+-----+     |   
    //  |                                                :     |
    //  :                                                :      > n * nblock
    //  :                                                :     | 
    //  :                                                :     |
    //  :                                                :     |
    //  :                              +-----+-----+-----+     |
    //  :                              | Au  |  0  |  0  |     |
    //  :                        +-----+-----+-----+-----+     |
    //  :                        |  0  | Ad  | Au  |  0  | n   |
    //  +-----+-----+---......---+-----+-----+=====+=====+    /
    //  |     |     |            |     |     !     |     !
    //  | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
    //  |     |     |            |     |     !     |     !
    //  +-----+-----+---......---+-----+-----+=====+=====+
    //
    */

    ///////////////////////////////////////////////////////

    valuePointer G_blk ;
    valuePointer AdAu_blk ;
    valuePointer LU_blk ; // last LU and working space
    valuePointer tmpV ;

    integer * ipiv_blk ;
    integer * LU_ipiv_blk ;
    integer   NB ; // blocking factor

    integer
    LU_2_block( integer      n,
                valuePointer A,
                valuePointer B,
                integer      ipiv[] ) const ;

    void
    reduce_block( integer      k,
                  integer      k1,
                  valuePointer y ) const ;

    void
    back_substitute( integer      k,
                     integer      k1,
                     integer      k2,
                     valuePointer y ) const ;

  public:

    explicit AmodioLU()
    : baseValue("AmodioLU_value")
    , baseInteger("AmodioLU_index")
    , NB(25)
    { }

    ~AmodioLU() {
      baseValue   . free() ;
      baseInteger . free() ;
    }

    //! load matrix in the class
    /*!
      \param nblk number of (square) blocks
      \param n    size of the blocks
      \param q    extra bc
      \param AdAu pointer to the blocks diagonal ad upper diagonal
      \param H0   pointer to the block \f$ H_0 \f$
      \param HN   pointer to the block \f$ H_N \f$
      \param Hq   pointer to the block \f$ H_q \f$
      \code
      Matrix structure
      
                      n * nblock
        ________________^_________________
       /                                  \
         n     n     n                        n     q
      +-----+-----+-----+----.........-----+-----+-----+    \
      |  Ad | Au  |  0  |                  |  0  |  0  | n   |
      +-----+-----+-----+             -----+-----+-----+     |
      |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |
      |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
      +-----+-----+-----+-----+       -----+-----+-----+     |   
      |                                                :     |
      :                                                :      > n * nblock
      :                                                :     | 
      :                                                :     |
      :                                                :     |
      :                              +-----+-----+-----+     |
      :                              | Au  |  0  |  0  |     |
      :                        +-----+-----+-----+-----+     |
      :                        |  0  | Ad  | Au  |  0  | n   |
      +-----+-----+---......---+-----+-----+=====+=====+    /
      |     |     |            |     |     !     |     !
      | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
      |     |     |            |     |     !     |     !
      +-----+-----+---......---+-----+-----+=====+=====+
      \endcode
    */
    void
    factorize( integer           nblk,
               integer           n,
               integer           q,
               valueConstPointer AdAu,
               valueConstPointer H0,
               valueConstPointer HN,
               valueConstPointer Hq ) ;

    //! solve linear sistem using internal factorized matrix
    void
    solve( valuePointer in_out ) const ;

  } ;
}

#endif
