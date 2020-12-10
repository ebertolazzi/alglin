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

///
/// file: BABD_SuperLU.hxx
///

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

    typedef double valueType;

    Malloc<valueType> m_baseValue;
    Malloc<int>       m_baseInteger;

    BABD_SuperLU(BABD_SuperLU const &);
    BABD_SuperLU const & operator = (BABD_SuperLU const &);

    integer m_number_of_blocks; //!< total number of blocks
    integer m_block_size;       //!< size of square blocks
    integer m;                  //!< number final rows (m>=m_block_size)
    integer m_nnz;              //!< total number of non zeros
    integer m_neq;

    int     * m_perm_r; // row permutations from partial pivoting
    int     * m_perm_c; // column permutation vector
    int     * m_etree;
    valueType m_one_norm_A;
    valueType m_inf_norm_A;

    superlu_options_t     m_slu_options;
    mutable SuperLUStat_t m_slu_stats;
    mutable SuperMatrix   m_A; // messo mutable per zittire warning
    mutable SuperMatrix   m_L; // messo mutable per zittire warning
    mutable SuperMatrix   m_U; // messo mutable per zittire warning

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
    factorize(
      integer         nblk,
      integer         n,
      integer         q,
      valueType const AdAu[],
      valueType const H0[],
      valueType const HN[],
      valueType const Hq[]
    );

    //! solve linear sistem using internal factorized matrix
    void solve( valueType in_out[] ) const;

    //! get condition number
    void cond( valueType & rcond_1, valueType & rcond_inf ) const;
  };
}

///
/// eof: BABD_SuperLU.hxx
///
