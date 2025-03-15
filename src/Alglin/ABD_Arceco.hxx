/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2020                                                      |
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
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: ABD_Arceco.hxx
///

namespace alglin {

  /*\
   *  A R C E C O
  \*/

  //!
  //!  \date    June 11, 2010
  //!  \version 1.0
  //!
  //!  \author  Enrico Bertolazzi and Daniele Bosetti
  //!
  //!  \par     Affiliation:
  //!           Department of Industrial Engineering<br>
  //!           University of Trento <br>
  //!           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
  //!           enrico.bertolazzi\@unitn.it
  //!
  //!  \par Abstract
  //!  This program solves the linear system A*X = B where A is
  //!  an almost block diagonal matrix. The method implemented is
  //!  based on Gauss elimination with alternate row and column
  //!  eliminaion with partial pivoting, which produces a stable
  //!  decomposition of the matrix A without introducing fill-in-
  //!
  //!  This class is an implementation of the Alternate Row and
  //!  Column Elimination in the C++ language.
  //!  Being almost block diagonal, the matrix is given by the
  //!  3-tuple (number_of_blocks,m_matrix_structure,array), where
  //!  number_of_blocks is the number of the blocks forming the matrix,
  //!  m_matrix_structure is an array which describes the structure
  //!  of the matrix, and array contains the data of the matrix.
  //!
  template <typename t_Value>
  class ArcecoLU {

    typedef t_Value real_type;

    Malloc<real_type> m_mem{"ArcecoLU_values"};
    Malloc<integer>   m_mem_int{"ArcecoLU_integers"};

    integer   * m_matrix_structure{nullptr}; //!< structure of the matrix
    integer   * m_pivot_array{nullptr};     //!< permutation array
    real_type * m_array{nullptr};           //!< the matrix data

    integer m_number_of_blocks{0};  //!< total number of blocks of the matrix A

    //!
    //!  row_elimination performs nrows_pivot row elimination on the matrix block.
    //!
    //!  \param block       pointer to the first element of the block
    //!  \param nrows_block number of rows of the block
    //!  \param ncols_block number of columns of the block
    //!  \param nrows_pivot number of rows to eliminate
    //!  \param pivot       pointer to a pivot array
    //!
    void
    row_elimination(
      real_type block[],
      integer   nrows_block,
      integer   ncols_block,
      integer   nrows_pivot,
      integer   pivot[]
    );

    //!
    //!  ColumnElimination performs ncols_pivot column elimination on the matrix
    //!  top-block and bottom-block.
    //!
    //!  \param topblk             pointer to the first element of the top block
    //!  \param nrows_top_block    number of rows of the top block
    //!  \param noverlap_cols      number of overlapping columns
    //!  \param botblk             pointer to the first element of the bottom block
    //!  \param nrows_bottom_block number of rows of the bottom block
    //!  \param ncols_pivot        number of columns to eliminate
    //!  \param pivot              pointer to a pivot array
    //!
    void
    column_elimination(
      real_type topblk[],
      integer   nrows_top_block,
      integer   noverlap_cols,
      real_type botblk[],
      integer   nrows_bottom_block,
      integer   ncols_pivot,
      integer   pivot[]
    );

    //!
    //! Performs the forward elimination step in the solution phase of solveByRef
    //!
    void
    forward_elimination(
      real_type const block[],
      integer         nrows_block,
      integer         nrows_pivot,
      integer   const pivot[],
      real_type       b[]
    ) const;

    //!
    //! Performs the forward solution step in the solution phase of solveByRef
    //!
    void
    forward_solution(
      real_type block[],
      integer   nrows_block,
      integer   ncols_pivot,
      integer   noverlap_cols,
      real_type b[]
    ) const;

    //!
    //! Performs the forward modification step in the solution phase of solve
    //!
    void
    forward_modification(
      real_type block[],
      integer   nrows_block,
      integer   ncols_pivot,
      real_type b[]
    ) const;

    //!
    //! Performs the backward modification step in the solution phase of solve
    //!
    void
    backward_modification(
      real_type block[],
      integer   nrows_block,
      integer   ncols_block,
      integer   nrows_pivot,
      real_type b[]
    ) const;

    //!
    //! Performs the backward substitution step in the solution phase of solve
    //!
    void
    backward_solution(
      real_type block[],
      integer   nrows_block,
      integer   ncols_block,
      integer   nrows_pivot,
      real_type b[]
    ) const;

    //!
    //! Performs the backward elimination step in the solution phase of solve
    //!
    void
    backward_elimination(
      real_type     block[],
      integer       nrows_block,
      integer       ncols_pivot,
      integer       noverlap_cols,
      integer const pivot[],
      real_type     b[]
    ) const;

    [[nodiscard]]
    integer
    nrows( integer num_block ) const
    { return m_matrix_structure[num_block*3+0]; }

    [[nodiscard]]
    integer
    ncols( integer num_block ) const
    { return m_matrix_structure[num_block*3+1]; }

    [[nodiscard]]
    integer
    noverlap( integer num_block ) const
    { return m_matrix_structure[num_block*3+2]; }

  public:

    explicit constexpr ArcecoLU() = default;
    ~ArcecoLU() = default;

    ArcecoLU(ArcecoLU<t_Value> const &) = delete;
    ArcecoLU<t_Value> const &operator = (ArcecoLU<t_Value> const &) = delete;

    //!
    //!  load_by_ref function gives to the class the sizes of the
    //!  problem and the pointers to the memory locations the class works on.
    //!
    //!  \param number_of_blocks Total number of blocks in A
    //!  \param pivot pointer to an array with n elements
    //!  \param matrix_structure pointer to an array with
    //!         3*number_of_blocks elements. Describes the block structure of A:
    //!         matrix_structure[3*i] = number of rows in the i-th block,
    //!         matrix_structure[3*i+1] = number of columns in the i-th block,
    //!         matrix_structure[3*i+2] = number of columns overlapped by
    //!         block i and block (i+1).
    //!  \param array pointer to an array with
    //!         SUM(matrix_structure[3*i]*matrix_structure[3*i+1]),
    //!         i=0...number_of_blocks-1, elements.
    //!         Contains the entries of the almost block diagonal system A whose
    //!         structure is given by the integer array matrix_structure.
    //!         The elements of A are stored by columns, in blocks corresponding
    //!         to the given structure.
    //!         The class will use this space to store the matrix decomposition.
    //!
    void
    load_by_ref(
      integer   number_of_blocks,
      integer   matrix_structure[],
      real_type array[],
      integer   pivot[]
    );

    //!
    //! \param neq the order of the linear system, and
    //!            n = SUM(matrix_structure[3*k],K=0,number_of_blocks-1)
    //!
    void
    checkStructure( integer neq ) const;

    //!
    //!  Decompose supervises the modified alternate row and column
    //!  decomposition with partial pivoting of the almost block
    //!  diagonal matrix A stored in the arrays array and matrix_structure.
    //!
    void
    factorize();

    //!
    //! factorize the matrix
    //!
    void
    factorize(
      integer         row0,
      integer         col0,
      real_type const block0[],
      // ----------------
      integer         num_block,
      integer         dim_block,
      real_type const blocks[],
      // ----------------
      integer         rowN,
      integer         colN,
      real_type const blockN[]
    );

    //!
    //!  Solve supervises the solution of the linear system
    //!
    //!                          A*X=B
    //!
    //!  using the decomposition of the matrix  A  already generated
    //!  in Decompose. It involves two loops, the forward loop,
    //!  consisting of forward solution, forward modification, and
    //!  forward elimination, and the backward loop, consisting of
    //!  backward solution, backward modification, and backward
    //!  elimination.
    //!
    void
    solve( real_type b[] ) const;

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  extern template class ArcecoLU<float>;
  extern template class ArcecoLU<double>;
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Arceco.hxx
///
