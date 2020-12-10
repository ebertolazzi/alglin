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
/// file: ABD_Diaz.hh
///

namespace alglin {

  //! LU decomposition of a ABD matrix
  /*!
   *
   * \date     May, 2016
   * \version  1.0
   *
   * \author   Enrico Bertolazzi
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           `enrico.bertolazzi\@unitn.it`
   *
   */
  /*
      Matrix NNZ structure
        col0-col00
     |    +----+----+                        |
     |    |    |    |                        | <- sizeBlock
     |    +----+----+----+                   |
     |         |    |    |                   |
     |         +----+----+----+              |
     |              |    |    |              |
     |              +----+----+----+         |
     |                   |    |    |         |
     |                   +----+----+---+     |
     |                        |        |     | <-- rowN
     |                        | BOTTOM |     |
     |    +----+              +--------+--+  |
     |    | TOP|                       |  |  | <-- row0
     |    +----+                       +--+  |
                                colN  col00

      Matrix NNZ internal structure
          col0
       |       |
     / +-------+                         \
     | |  TOP  |                         | <-- row0
     | +--+----+----+                    |
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
  class DiazLU : public BlockBidiagonal<t_Value> {
  public:

    typedef t_Value valueType;

  private:

    DiazLU( DiazLU const & ) = delete;
    DiazLU const & operator = ( DiazLU const & ) = delete;

    integer m_NB;

    mutable integer m_nblk;

    integer * m_swapRC_blks;

    void
    LU_left_right(
      integer nrA,
      integer ncA,
      integer ncL,
      integer ncR,
      t_Value * A, integer ldA,
      integer swapR[]
    );

    void
    LU_top_bottom(
      integer nrT,
      integer nrA,
      integer ncA,
      t_Value * A, integer ldA,
      integer nrB,
      t_Value * B, integer ldB,
      integer swapC[]
    );

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal( bool do_permute, valueType in_out[] ) const;

    //! solve linear sistem using internal factorized matrix
    void
    solve_internal(
      bool      do_permute,
      integer   nrhs,
      valueType in_out[],
      integer   ldRhs
    ) const;

  public:

    using BlockBidiagonal<t_Value>::m_baseInteger;

    using BlockBidiagonal<t_Value>::m_number_of_blocks;
    using BlockBidiagonal<t_Value>::m_block_size;

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

    using BlockBidiagonal<t_Value>::factorize;
    using BlockBidiagonal<t_Value>::dump_ccoord;

    using BlockBidiagonal<t_Value>::m_la_factorization;

    explicit UTILS_CONSTEXPR DiazLU() : m_NB(25) {}

    virtual
    ~DiazLU() UTILS_OVERRIDE
    {}

    virtual
    void
    allocate(
      integer /* nblock */,
      integer /* n      */,
      integer /* nb     */,
      // ----------------------
      integer /* numInitialBC */,
      integer /* numFinalBC   */,
      integer /* numCyclicBC  */,
      // ----------------------
      integer /* numInitialOMEGA */,
      integer /* numFinalOMEGA   */,
      integer /* numCyclicOMEGA  */
    ) UTILS_OVERRIDE
    { UTILS_ERROR0("DiazLU::allocate() not defined!\n"); }

    virtual
    void
    allocateTopBottom(
      integer _nblock,
      integer _n,
      integer _row0,
      integer _col0,
      integer _rowN,
      integer _colN,
      integer _nb
    ) UTILS_OVERRIDE {
      integer inv = _nblock*_n+(_col0+_colN-2*_n);
      BlockBidiagonal<t_Value>::allocateTopBottom(
        _nblock, _n,
        _row0, _col0,
        _rowN, _colN,
        _nb, 0, inv
      );
      m_swapRC_blks = m_baseInteger(size_t(inv));
    }

    virtual
    void
    factorize() UTILS_OVERRIDE;

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( valueType in_out[] ) const UTILS_OVERRIDE
    { solve_internal( true, in_out ); }

    //! solve linear sistem using internal factorized matrix
    virtual
    void
    solve( integer nrhs, valueType in_out[], integer ldRhs ) const UTILS_OVERRIDE
    { solve_internal( true, nrhs, in_out, ldRhs ); }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD( valueType in_out[] ) const
    { solve_internal( false, in_out ); }

    //! solve linear sistem using internal factorized matrix
    void
    solve_ABD( integer nrhs, valueType in_out[], integer ldRhs ) const
    { solve_internal( false, nrhs, in_out, ldRhs ); }

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
  #endif

  extern template class DiazLU<float>;
  extern template class DiazLU<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Diaz.hh
///

