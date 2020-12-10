/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                       |
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
/// file: BABD_Block.hxx
///

//! Various LU decomposition classes
namespace alglin {

  /*
  //   ____  ____  _            _    _    _   _
  //  | __ )| __ )| | ___   ___| | _| |  | | | |
  //  |  _ \|  _ \| |/ _ \ / __| |/ / |  | | | |
  //  | |_) | |_) | | (_) | (__|   <| |__| |_| |
  //  |____/|____/|_|\___/ \___|_|\_\_____\___/
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
   *           Department of Industrial Engineering<br>
   *           University of Trento<br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   */
  template <typename t_Value>
  class BBlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value valueType;

  private:

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;
    using BlockBidiagonal<t_Value>::m_extra_bc;

    using BlockBidiagonal<t_Value>::m_nx2;
    using BlockBidiagonal<t_Value>::m_nxnx2;
    using BlockBidiagonal<t_Value>::m_nxn;

    using BlockBidiagonal<t_Value>::m_numInitialBC;
    using BlockBidiagonal<t_Value>::m_numFinalBC;
    using BlockBidiagonal<t_Value>::m_numCyclicBC;
    using BlockBidiagonal<t_Value>::m_numInitialOMEGA;
    using BlockBidiagonal<t_Value>::m_numFinalOMEGA;
    using BlockBidiagonal<t_Value>::m_numCyclicOMEGA;

    using BlockBidiagonal<t_Value>::m_DE_blk;
    using BlockBidiagonal<t_Value>::m_H0Nq;
    using BlockBidiagonal<t_Value>::m_block0;
    using BlockBidiagonal<t_Value>::m_blockN;
    using BlockBidiagonal<t_Value>::m_Bmat;
    using BlockBidiagonal<t_Value>::m_Cmat;
    using BlockBidiagonal<t_Value>::m_Dmat;

    BBlockLU(BBlockLU const &);
    BBlockLU const & operator = (BBlockLU const &);

    integer m;   //!< number final rows (m>=n)
    integer N;   //!< n * (nblock+1) + q
    integer nnz; //!< total number of non zeros

    /*!
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

    /*!
    //
    //  Working block AdH_blk [ size = nblock * ( n * (n+m) ) ]
    //
    //                 n * nblock
    //    ________________^_________________
    //   /                                  \
    //     n     n     n
    //  +-----+-----+-----+----........+-----+
    //  |  Ad | Ad  | Ad  |               Ad | n
    //  +-----+-----+-----+            +-----+
    //  |     |     |            |     |     !
    //  | H0  |  0  |            |     |  0  ! m
    //  |     |     |            |     |     !
    //  +-----+-----+---......---+-----+-----+
    //
    */
    valueType * m_AdH_blk;

    /*!
    //
    //  Working block Au_blk [ size = nblock * (n*n) ]
    //
    //                  n * nblock + q
    //    ________________^________________________
    //   /                                         \
    //     n     n     n                  n      q
    //  +-----+-----+-----+----.........-----+------+
    //  |  Au | Au  | Au  |               Au | fill | n
    //  +-----+-----+-----+----.........-----+------+
    //
    */
    valueType * m_Au_blk;

    /*!
    //
    //  Working block DD_blk [ size = m * m ]
    //
    //     n     q     (n+q=m)
    //  +=====+=====+
    //  !     |     !
    //  ! HN  | Hq  ! m
    //  !     |     !
    //  +=====+=====+
    //
    */
    valueType * m_DD_blk;

    /*!
    //
    //  Working block FF_blk [ size = nblock * (n*m) ]
    //
    //     n     q        (n+q=m)
    //  +-----+-----+   \
    //  |fill |fill | n |
    //  +-----+-----+   |
    //  |fill |fill | n |
    //  +-----+-----+   |
    //  |fill |fill | n |
    //  +-----+-----+   :
    //  |           |   :  n * (nblock-1)
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  :           :   :
    //  |fill |fill | n |
    //  +-----+-----+   /
    //
    */
    valueType * m_FF_blk;

    // pivot vector
    integer * m_ipiv_blk;

  public:

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    explicit BBlockLU() { }

    virtual
    ~BBlockLU() UTILS_OVERRIDE
    { }

    virtual
    void
    factorize() UTILS_OVERRIDE;

    //! solve linear system previously factorized
    virtual
    void
    solve( valueType in_out[] ) const UTILS_OVERRIDE;

  };

  // explicit instantiation declaration to suppress warnings

  extern template class BBlockLU<float>;
  extern template class BBlockLU<double>;

}

///
/// eof: BABD_Block.hxx
///
