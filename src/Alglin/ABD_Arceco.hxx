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
 |      Universita` degli Studi di Trento                                   |
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

  /*!
   *
   *  \date    June 11, 2010
   *  \version 1.0
   *
   *  \author  Enrico Bertolazzi and Daniele Bosetti
   *
   *  \par     Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi\@unitn.it
   *
   *  \par Abstract
   *  This program solves the linear system A*X = B where A is
   *  an almost block diagonal matrix. The method implemented is
   *  based on Gauss elimination with alternate row and column
   *  eliminaion with partial pivoting, which produces a stable
   *  decomposition of the matrix A without introducing fill-in-
   *
   *  This class is an implementation of the Alternate Row and
   *  Column Elimination in the C++ language.
   *  Being almost block diagonal, the matrix is given by the
   *  3-tuple (numberOfBlocks,m_matrix_structure,array), where
   *  numberOfBlocks is the number of the blocks forming the matrix,
   *  m_matrix_structure is an array which describes the structure
   *  of the matrix, and array contains the data of the matrix.
  \*/
  template <typename t_Value>
  class ArcecoLU {

    typedef t_Value real_type;

    Malloc<real_type> m_baseValue;
    Malloc<integer>   m_baseInteger;

    ArcecoLU(ArcecoLU<t_Value> const &);
    ArcecoLU<t_Value> const &operator = (ArcecoLU<t_Value> const &);

    integer   * m_matrix_structure; //!< structure of the matrix
    integer   * m_pivot_array;     //!< permutation array
    real_type * m_array;           //!< the matrix data

    integer m_number_of_blocks;  //!< total number of blocks of the matrix A

    /*!
     *  RowElimination performs numRowsPivot row elimination on the matrix block.
     *  \param block        pointer to the first element of the block
     *  \param numRowsBlock number of rows of the block
     *  \param numColsBlock number of columns of the block
     *  \param numRowsPivot number of rows to eliminate
     *  \param pivot        pointer to a pivot array
     */
    void
    rowElimination(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsBlock,
      integer   numRowsPivot,
      integer   pivot[]
    );

    /*!
     *  ColumnElimination performs numColsPivot column elimination on the matrix top-block and bottom-block.
     *  \param topblk             pointer to the first element of the top block
     *  \param numRowsTopBlock    number of rows of the top block
     *  \param numOverlapCols     number of overlapping columns
     *  \param botblk             pointer to the first element of the bottom block
     *  \param numRowsBottomBlock number of rows of the bottom block
     *  \param numColsPivot       number of columns to eliminate
     *  \param pivot              pointer to a pivot array
     */
    void
    columnElimination(
      real_type topblk[],
      integer   numRowsTopBlock,
      integer   numOverlapCols,
      real_type botblk[],
      integer   numRowsBottomBlock,
      integer   numColsPivot,
      integer   pivot[]
    );

    //! Performs the forward elimination step in the solution phase of solveByRef
    void
    forwardElimination(
      real_type block[],
      integer   numRowsBlock,
      integer   numRowsPivot,
      integer   pivot[],
      real_type b[]
    ) const;

    //! Performs the forward solution step in the solution phase of solveByRef
    void
    forwardSolution(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsPivot,
      integer   numOverlapCols,
      real_type b[]
    ) const;

    //! Performs the forward modification step in the solution phase of solve
    void
    forwardModification(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsPivot,
      real_type b[]
    ) const;

    //! Performs the backward modification step in the solution phase of solve
    void
    backwardModification(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsBlock,
      integer   numRowsPivot,
      real_type b[]
    ) const;

    //! Performs the backward substitution step in the solution phase of solve
    void
    backwardSolution(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsBlock,
      integer   numRowsPivot,
      real_type b[]
    ) const;

    //! Performs the backward elimination step in the solution phase of solve
    void
    backwardElimination(
      real_type block[],
      integer   numRowsBlock,
      integer   numColsPivot,
      integer   numOverlapCols,
      integer   pivot[],
      real_type b[]
    ) const;

    integer
    numRows( integer numBlock ) const
    { return m_matrix_structure[numBlock*3+0]; }

    integer
    numCols( integer numBlock ) const
    { return m_matrix_structure[numBlock*3+1]; }

    integer
    numOverlap( integer numBlock ) const
    { return m_matrix_structure[numBlock*3+2]; }

  public:

    explicit constexpr ArcecoLU()
    : m_baseValue("ArcecoLU_values")
    , m_baseInteger("ArcecoLU_integers")
    {}

    ~ArcecoLU()
    {}

    /*!
     *  loadByRef function gives to the class the sizes of the
     *  problem and the pointers to the memory locations the class works on.
     *
     *  \param numberOfBlocks Total number of blocks in A
     *  \param pivot pointer to an array with n elements
     *  \param matrixStructure pointer to an array with
     *         3*numberOfBlocks elements. Describes the block structure of A:
     *         matrixStructure[3*i] = number of rows in the i-th block,
     *         matrixStructure[3*i+1] = number of columns in the i-th block,
     *         matrixStructure[3*i+2] = number of columns overlapped by
     *         block i and block (i+1).
     *  \param array pointer to an array with
     *         SUM(matrixStructure[3*i]*matrixStructure[3*i+1]),
     *         i=0...numberOfBlocks-1, elements.
     *         Contains the entries of the almost block diagonal system A whose
     *         structure is given by the integer array matrixStructure.
     *         The elements of A are stored by columns, in blocks corresponding
     *         to the given structure.
     *         The class will use this space to store the matrix decomposition.
     */
    void
    loadByRef(
      integer   numberOfBlocks,
      integer   matrixStructure[],
      real_type array[],
      integer   pivot[]
    );

    //! \@param neq the order of the linear system, and n = SUM(matrixStructure[3*k],K=0,numberOfBlocks-1)
    void
    checkStructure( integer neq );

    /*!
     *  Decompose supervises the modified alternate row and column
     *  decomposition with partial pivoting of the almost block
     *  diagonal matrix A stored in the arrays array  and
     *  matrixStructure.
     */
    void
    factorize();

    //! factorize the matrix
    void
    factorize(
      integer         row0,
      integer         col0,
      real_type const block0[],
      // ----------------
      integer         numBlock,
      integer         dimBlock,
      real_type const blocks[],
      // ----------------
      integer         rowN,
      integer         colN,
      real_type const blockN[]
    );

    /*!
     *  Solve supervises the solution of the linear system
     *                          A*X=B
     *  using the decomposition of the matrix  A  already generated
     *  in Decompose. It involves two loops, the forward loop,
     *  consisting of forward solution, forward modification, and
     *  forward elimination, and the backward loop, consisting of
     *  backward solution, backward modification, and backward
     *  elimination.
     */
    void
    solve( real_type b[] ) const;

  };

  // explicit instantiation declaration to suppress warnings

  #ifdef __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
  #endif

  extern template class ArcecoLU<float>;
  extern template class ArcecoLU<double>;

  #ifdef __clang__
  #pragma clang diagnostic pop
  #endif
}

///
/// eof: ABD_Arceco.hxx
///
