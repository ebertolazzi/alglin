/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                       |
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
  //!
  //!  LU decomposition of a BABD matrix
  //!
  //!  \date     May 30, 2006
  //!  \version  1.0
  //!  \note     first release May 30, 2006
  //!
  //!  \author   Enrico Bertolazzi
  //!
  //!  \par      Affiliation:
  //!            Department of Industrial Engineering<br>
  //!            University of Trento<br>
  //!            Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
  //!            enrico.bertolazzi\@unitn.it
  //!
  template <typename t_Value>
  class BBlockLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value real_type;

  private:

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;
    using BlockBidiagonal<t_Value>::m_extra_bc;
    using BlockBidiagonal<t_Value>::n_x_n;

    using BlockBidiagonal<t_Value>::m_num_initial_BC;
    using BlockBidiagonal<t_Value>::m_num_final_BC;
    using BlockBidiagonal<t_Value>::m_num_cyclic_BC;
    using BlockBidiagonal<t_Value>::m_num_initial_OMEGA;
    using BlockBidiagonal<t_Value>::m_num_final_OMEGA;
    using BlockBidiagonal<t_Value>::m_num_cyclic_OMEGA;

    using BlockBidiagonal<t_Value>::m_DE_blk;
    using BlockBidiagonal<t_Value>::m_H0Nq;
    using BlockBidiagonal<t_Value>::m_block0;
    using BlockBidiagonal<t_Value>::m_blockN;
    using BlockBidiagonal<t_Value>::m_Bmat;
    using BlockBidiagonal<t_Value>::m_Cmat;
    using BlockBidiagonal<t_Value>::m_Dmat;

    BBlockLU( BBlockLU const & ) = delete;
    BBlockLU const & operator = ( BBlockLU const & ) = delete;

    //!
    //!  Matrix structure
    //!
    //!  \verbatim
    //!
    //!                  n * nblock
    //!    ________________^_________________
    //!   /                                  \
    //!     n     n     n                        n     q
    //!  +-----+-----+-----+----.........-----+-----+-----+    \
    //!  |  Ad | Au  |  0  |                  |  0  |  0  | n   |
    //!  +-----+-----+-----+             -----+-----+-----+     |
    //!  |  0  | Ad  | Au  |  0               |  0  |  0  | n   |
    //!  +-----+-----+-----+-----+       -----+-----+-----+     |
    //!  |  0  |  0  | Ad  | Au  |            |  0  |  0  | n   |
    //!  +-----+-----+-----+-----+       -----+-----+-----+     |
    //!  |                                                :     |
    //!  :                                                :      > n * nblock
    //!  :                                                :     |
    //!  :                                                :     |
    //!  :                                                :     |
    //!  :                              +-----+-----+-----+     |
    //!  :                              | Au  |  0  |  0  |     |
    //!  :                        +-----+-----+-----+-----+     |
    //!  :                        |  0  | Ad  | Au  |  0  | n   |
    //!  +-----+-----+---......---+-----+-----+=====+=====+    /
    //!  |     |     |            |     |     !     |     !
    //!  | H0  |  0  |            |     |  0  ! HN  | Hq  ! m
    //!  |     |     |            |     |     !     |     !
    //!  +-----+-----+---......---+-----+-----+=====+=====+
    //!
    //!  \endverbatim
    //!
    //!  Working block AdH_blk [ size = nblock * ( n * (n+m) ) ]
    //!
    //!  \verbatim
    //!
    //!                 n * nblock
    //!    ________________^_________________
    //!   /                                  \
    //!     n     n     n
    //!  +-----+-----+-----+----........+-----+
    //!  |  Ad | Ad  | Ad  |               Ad | n
    //!  +-----+-----+-----+            +-----+
    //!  |     |     |            |     |     !
    //!  | H0  |  0  |            |     |  0  ! m
    //!  |     |     |            |     |     !
    //!  +-----+-----+---......---+-----+-----+
    //!
    //!  \endverbatim
    //!
    real_type * m_AdH_blk{nullptr};

    //!
    //!
    //!  Working block Au_blk [ size = nblock * (n*n) ]
    //!
    //!  \verbatim
    //!
    //!                  n * nblock + q
    //!    ________________^________________________
    //!   /                                         \
    //!     n     n     n                  n      q
    //!  +-----+-----+-----+----.........-----+------+
    //!  |  Au | Au  | Au  |               Au | fill | n
    //!  +-----+-----+-----+----.........-----+------+
    //!
    //!  \endverbatim
    //!
    real_type * m_Au_blk{nullptr};

    //!
    //!  Working block DD_blk [ size = m * m ]
    //!
    //!  \verbatim
    //!
    //!     n     q     (n+q=m)
    //!  +=====+=====+
    //!  !     |     !
    //!  ! HN  | Hq  ! m
    //!  !     |     !
    //!  +=====+=====+
    //!
    //!  \endverbatim
    //!
    real_type * m_DD_blk{nullptr};

    //!
    //!  Working block FF_blk [ size = nblock * (n*m) ]
    //!
    //!  \verbatim
    //!
    //!     n     q        (n+q=m)
    //!  +-----+-----+   \
    //!  |fill |fill | n |
    //!  +-----+-----+   |
    //!  |fill |fill | n |
    //!  +-----+-----+   |
    //!  |fill |fill | n |
    //!  +-----+-----+   :
    //!  |           |   :  n * (nblock-1)
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  :           :   :
    //!  |fill |fill | n |
    //!  +-----+-----+   /
    //!
    //!  \endverbatim
    //!
    real_type * m_FF_blk{nullptr};

    // pivot vector
    integer * m_ipiv_blk{nullptr};

  public:

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    explicit BBlockLU() = default;

    virtual
    void
    factorize() override;

    //!
    //! solve linear system previously factorized
    //!
    virtual
    void
    solve( real_type in_out[] ) const override;

  };

  // explicit instantiation declaration to suppress warnings
  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  extern template class BBlockLU<float>;
  extern template class BBlockLU<double>;
  #pragma clang diagnostic pop
  #endif

}

///
/// eof: BABD_Block.hxx
///
