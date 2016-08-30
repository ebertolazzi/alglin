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

#ifndef ABD_ARCECO_HH
#define ABD_ARCECO_HH

#include "Alglin.hh"
#include "Alglin++.hh"

namespace alglin {

  /*
   *  A R C E C O
   */

  /*!
   *
   * \date     June 11, 2010
   * \version  1.0
   *
   * \author   Enrico Bertolazzi and Daniele Bosetti
   *
   * \par      Affiliation:
   *           Department of Industrial Engineering<br>
   *           University of Trento <br>
   *           Via Sommarive 9, I-38123 Povo, Trento, Italy<br>
   *           enrico.bertolazzi@unitn.it
   *
   * \par Abstract
   *  This program solves the linear system A*X = B where A is
   *  an almost block diagonal matrix. The method implemented is
   *  based on Gauss elimination with alternate row and column
   *  eliminaion with partial pivoting, which produces a stable
   *  decomposition of the matrix A without introducing fill-in-
   *
   *  This class is an implementation of the Alternate Row and Column Elimination in the C++ language.
   *  Being almost block diagonal, the matrix is given by the 3-tuple (numberOfBlocks,matrixStructure,array), where
   *  numberOfBlocks is the number of the blocks forming the matrix, matrixStructure is an array which describes the structure
   *  of the matrix, and array contains the data of the matrix.
   */
  template <typename t_Value>
  class ArcecoLU {

    typedef t_Value         valueType ;
    typedef t_Value*        valuePointer ;
    typedef t_Value const * valueConstPointer ;

    Malloc<valueType> baseValue ;
    Malloc<integer>   baseInteger ;

    ArcecoLU(ArcecoLU<t_Value> const &);
    ArcecoLU<t_Value> const &operator = (ArcecoLU<t_Value> const &) ;

    integer *    matrixStructure ; //!< structure of the matrix
    integer *    pivot           ; //!< permutation array
    valuePointer array           ; //!< the matrix data

    integer      numberOfBlocks  ; //!< total number of blocks of the matrix A
    integer      __padding ;

    /*!
     *  RowElimination performs numRowsPivot row elimination on the matrix block.
     *  \param block        pointer to the first element of the block
     *  \param numRowsBlock number of rows of the block
     *  \param numColsBlock number of columns of the block
     *  \param numRowsPivot number of rows to eliminate
     *  \param pivot        pointer to a pivot array
     */
    void
    rowElimination ( valuePointer block,
                     integer      numRowsBlock,
                     integer      numColsBlock,
                     integer      numRowsPivot,
                     integer    * pivot ) ;

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
    columnElimination ( valuePointer topblk,
                        integer      numRowsTopBlock,
                        integer      numOverlapCols,
                        valuePointer botblk,
                        integer      numRowsBottomBlock,
                        integer      numColsPivot,
                        integer *    pivot ) ;

    //! Performs the forward elimination step in the solution phase of solveByRef
    void
    forwardElimination ( valuePointer block,
                         integer      numRowsBlock,
                         integer      numRowsPivot,
                         integer *    pivot,
                         valuePointer b ) const ;

    //! Performs the forward solution step in the solution phase of solveByRef
    void
    forwardSolution ( valuePointer block,
                      integer      numRowsBlock,
                      integer      numColsPivot,
                      integer      numOverlapCols,
                      valuePointer b ) const ;

    //! Performs the forward modification step in the solution phase of solve
    void
    forwardModification ( valuePointer block,
                          integer      numRowsBlock,
                          integer      numColsPivot,
                          valuePointer b ) const ;

    //! Performs the backward modification step in the solution phase of solve
    void
    backwardModification ( valuePointer block,
                           integer      numRowsBlock,
                           integer      numColsBlock,
                           integer      numRowsPivot,
                           valuePointer b ) const ;

    //! Performs the backward substitution step in the solution phase of solve
    void
    backwardSolution ( valuePointer block,
                       integer      numRowsBlock,
                       integer      numColsBlock,
                       integer      numRowsPivot,
                       valuePointer b ) const ;

    //! Performs the backward elimination step in the solution phase of solve
    void
    backwardElimination ( valuePointer block,
                          integer      numRowsBlock,
                          integer      numColsPivot,
                          integer      numOverlapCols,
                          integer *    pivot,
                          valuePointer b ) const ;

    integer numRows    ( integer numBlock ) const { return matrixStructure[numBlock*3+0] ; }
    integer numCols    ( integer numBlock ) const { return matrixStructure[numBlock*3+1] ; }
    integer numOverlap ( integer numBlock ) const { return matrixStructure[numBlock*3+2] ; }

  public:
  
    explicit ArcecoLU()
    : baseValue("ArcecoLU_values")
    , baseInteger("ArcecoLU_integers")
    {}

    ~ArcecoLU()
    {}

    /*!
     *  loadByRef function gives to the class the sizes of the problem and the pointers to the memory
     *  locations the class works on.
     *
     *  \param numberOfBlocks Total number of blocks in A
     *  \param pivot pointer to an array with n elements
     *  \param matrixStructure pointer to an array with 3*numberOfBlocks elements. Describes the block structure of A:
     *         matrixStructure[3*i] = number of rows in the i-th block,
     *         matrixStructure[3*i+1] = number of columns in the i-th block,
     *         matrixStructure[3*i+2] = number of columns overlapped by block i and block (i+1).
     *  \param array pointer to an array with SUM(matrixStructure[3*i]*matrixStructure[3*i+1]),i=0...numberOfBlocks-1, elements.
     *  Contains the entries of the almost block diagonal system A whose structure is given
     *  by the integer array matrixStructure. The elements of A are stored by columns, in blocks corresponding
     *  to the given structure. The class will use this space to store the matrix decomposition.
     */
    void
    loadByRef ( integer      numberOfBlocks,
                integer *    matrixStructure,
                valuePointer array,
                integer *    pivot ) ;

    //! \@param neq the order of the linear system, and n = SUM(matrixStructure[3*k],K=0,numberOfBlocks-1)
    void
    checkStructure( integer neq ) ;

    /*!
     *  Decompose supervises the modified alternate row and column
     *  decomposition with partial pivoting of the almost block
     *  diagonal matrix A stored in the arrays array  and
     *  matrixStructure.
     */
    void
    factorize() ;

    //! factorize the matrix
    void
    factorize( integer           _row0,
               integer           _col0,
               valueConstPointer _block0,
               // ----------------
               integer           _numBlock,
               integer           _dimBlock,
               valueConstPointer _blocks,
               // ----------------
               integer           _rowN,
               integer           _colN,
               valueConstPointer _blockN ) ;

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
    solve( valuePointer b ) const ;

  };
}

#endif
