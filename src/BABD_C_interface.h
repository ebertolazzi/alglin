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

#ifndef BABD_C_INTERFACE_HH
#define BABD_C_INTERFACE_HH

/*
  Use -lstdc++ when linking with C code.
*/

#ifdef __cplusplus
extern "C" {
#endif

typedef int    ABD_intType ;
typedef double ABD_realType ;

/*!
 *  perform Diaz factorization of ABD matrix defined by blocks TOP, BOTTOM, D and E.
 *  Matrix structure is the following:
 *
 *           col0
 *        |       |
 *      / +-------+                         \
 *      | |  TOP  |                         | <-- row0
 *      | +--+----+----+                    |
 *      |    | D1 | E1 |                    | <- sizeBlock
 *      |    +----+----+----+               |
 *      |         | D2 | E2 |               |
 *      |         +----+----+----+          |
 *      |               ...........         |
 *      |              +----+----+----+     |
 *      |                   | DN | EN |     |
 *      |                   +----+----+---+ |
 *      |                        |        | | <-- rowN
 *      |                        | BOTTOM | |
 *      \                        +--------+ /
 *                               |        |
 *                                  colN
 *
 *  \param mat_id   identifier for the factorization, used in the subsequent `ABD_solve`
 *  \param row0     number of rows of the `TOP` block
 *  \param col0     number of cols of the `TOP` block
 *  \param TOP      pointer of the `TOP` block stored by column (FORTRAN STORAGE)
 *  \param ldTOP    leading dimension of matrix `TOP`
 *  \param nblock   number of blocks `D` and `E`
 *  \param n        dimension of the blocks `D` and `E` (size `n` x `n`)
 *  \param DE       pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
 *                  The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
 *  \param ldDE     leading dimension of matrices `D` and `E`
 *  \param rowN     number of rows of the `BOTTOM` block
 *  \param colN     number of cols of the `BOTTOM` block
 *  \param BOTTOM   pointer of the `BOTTOM` block stored by column (FORTRAN STORAGE)
 *  \param ldBOTTOM leading dimension of matrix `BOTTOM`
 *
 *  \return  0 no error found
 */

int
ABD_factorize( ABD_intType        mat_id,
               ABD_intType        row0,
               ABD_intType        col0,
               ABD_realType const TOP[], ABD_intType ldTOP,
               ABD_intType        nblock,
               ABD_intType        n,
               ABD_realType const DE[], ABD_intType ldDE,
               ABD_intType        rowN,
               ABD_intType        colN,
               ABD_realType const BOTTOM[], ABD_intType ldBOTTOM ) ;

/*!
 *  solve linear ABD system using factorization of `ABD_factorize` call
 *
 *  \param mat_id   identifier for the factorization
 *  \param rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
 *
 *  \return  0 no error found
 */

int
ABD_solve( ABD_intType mat_id, ABD_realType rhs_sol[] ) ;


/*!
 *  solve linear ABD system using factorization of `ABD_factorize` call
 *
 *  \param mat_id   identifier for the factorization
 *  \param nrhs     number of rhs
 *  \param rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
 *  \param ldRhs    leadind dimension of `rhs_sol`
 *
 *  \return  0 no error found
 */

int
ABD_solve_nrhs( ABD_intType  mat_id,
                ABD_intType  nrhs,
                ABD_realType rhs_sol[],
                ABD_intType  ldRhs ) ;

/*!
 *  destroy (and free mempry) of a factorization of ABD matrix.
 *  \param mat_id identifier for the factorization
 *
 *  \return  0 no error found
 */

int
ABD_free( ABD_intType mat_id ) ;

/*!
 *  \return pointer to a string with the last error found
 */

char const *
ABD_get_last_error( ) ;

#ifdef __cplusplus
}
#endif

#endif
