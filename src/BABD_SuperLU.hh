/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
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

#ifndef BABD_SUPERLU_HH
#define BABD_SUPERLU_HH

#include "AlglinConfig.hh"
#include "Alglin.hh"
#include "Alglin++.hh"
#include "AlglinSuperLU.hh"

#ifdef __GCC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpadded"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpadded"
#endif

namespace alglin {

  /*
  //   ____                        _    _   _
  //  / ___| _   _ _ __   ___ _ __| |  | | | |
  //  \___ \| | | | '_ \ / _ \ '__| |  | | | |
  //   ___) | |_| | |_) |  __/ |  | |__| |_| |
  //  |____/ \__,_| .__/ \___|_|  |_____\___/
  //              |_|
  */

  //! LU decomposition of a BABD matrix
  /*!
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Dipartimento di Ingegneria Industriale<br>
   *           Universita di Trento<br>
   *           via Sommarive 9, I -- 38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   */
  class BABD_SuperLU {
  private:

    typedef double            valueType;
    typedef valueType *       valuePointer;
    typedef valueType const * valueConstPointer;

    Malloc<valueType> baseValue;
    Malloc<int>       baseInteger;

    BABD_SuperLU(BABD_SuperLU const &);
    BABD_SuperLU const & operator = (BABD_SuperLU const &);

    integer nblock; //!< total number of blocks
    integer n;      //!< size of square blocks
    integer m;      //!< number final rows (m>=n)
    integer nnz;    //!< total number of non zeros
    integer neq;

    //int               info;
    //mem_usage_t       mem_usage;
    int         *perm_r; // row permutations from partial pivoting
    int         *perm_c; // column permutation vector
    int         *etree;
    valueType   one_norm_A, inf_norm_A;

    superlu_options_t     slu_options;
    mutable SuperLUStat_t slu_stats;
    mutable SuperMatrix   A, L, U; // messo mutable per zittire warning

    //! factorize the matrix
    void factorize();

  public:

    explicit BABD_SuperLU();

    ~BABD_SuperLU();

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
               valueConstPointer Hq );

    //! solve linear sistem using internal factorized matrix
    void solve( valuePointer in_out ) const;

    //! get condition number
    void cond( valueType & rcond_1, valueType & rcond_inf ) const;
  };
}

#ifdef __GCC__
#pragma GCC diagnostic pop
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#endif
