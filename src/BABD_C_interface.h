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

#ifndef BABD_C_INTERFACE_HH
#define BABD_C_INTERFACE_HH

/*
  Use -lstdc++ when linking with C code.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*\
 |     _    ____  ____
 |    / \  | __ )|  _ \
 |   / _ \ |  _ \| | | |
 |  / ___ \| |_) | |_| |
 | /_/   \_\____/|____/
 |
\*/

typedef int    ABD_intType;
typedef double ABD_realType;

/*!
 *  perform Diaz factorization of ABD matrix defined by blocks TOP, BOTTOM, D and E.
 *  Matrix structure is the following:
 *
 *  \verbatim
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
 *  \endverbatim
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
ABD_factorize(
  ABD_intType        mat_id,
  ABD_intType        row0,
  ABD_intType        col0,
  ABD_realType const TOP[], ABD_intType ldTOP,
  ABD_intType        nblock,
  ABD_intType        n,
  ABD_realType const DE[], ABD_intType ldDE,
  ABD_intType        rowN,
  ABD_intType        colN,
  ABD_realType const BOTTOM[], ABD_intType ldBOTTOM
);

/*!
 *  solve linear ABD system using factorization of `ABD_factorize` call
 *
 *  \param mat_id   identifier for the factorization
 *  \param rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
 *
 *  \return  0 no error found
 */

int
ABD_solve( ABD_intType mat_id, ABD_realType rhs_sol[] );


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
ABD_solve_nrhs(
  ABD_intType  mat_id,
  ABD_intType  nrhs,
  ABD_realType rhs_sol[],
  ABD_intType  ldRhs
);

/*!
 *  destroy (and free mempry) of a factorization of ABD matrix.
 *  \param mat_id identifier for the factorization
 *
 *  \return  0 no error found
 */

int
ABD_free( ABD_intType mat_id );

/*!
 *  \return pointer to a string with the last error found
 */

char const *
ABD_get_last_error(void);

void
ABD_get_last_error_f90( char [], long len );

/*\
 |  ____    _    ____  ____
 | | __ )  / \  | __ )|  _ \
 | |  _ \ / _ \ |  _ \| | | |
 | | |_) / ___ \| |_) | |_| |
 | |____/_/   \_\____/|____/
 |
\*/

typedef int    BABD_intType;
typedef double BABD_realType;

/*!
 *  Perform factorization of BABD matrix defined by blocks H, D and E.
 *  Matrix can be bordered.
 *
 *  \verbatim
 *  Matrix structure
 *
 *                     (n+1) * nblock
 *        ___________________^____________________
 *       /                                        \
 *         n     n     n                        n
 *      +-----+-----+-----+----.........-----+-----+    \
 *      |  D  |  E  |  0  |                  |  0  | n   |
 *      +-----+-----+-----+             -----+-----+     |
 *      |  0  |  D  |  E  |  0               |  0  | n   |
 *      +-----+-----+-----+-----+       -----+-----+     |
 *      |  0  |  0  |  D  |  E  |            |  0  | n   |
 *      +-----+-----+-----+-----+       -----+-----+     |
 *      |                                                |
 *  A = :                                                 > n * nblock
 *      :                                                |
 *      :                                                |
 *      :                                                |
 *      :                              +-----+-----+     |
 *      :                              |  E  |  0  |     |
 *      :                        +-----+-----+-----+     |
 *      :                        |  0  |  D  |  E  | n   |
 *      +-----+-----+---......---+-----+-----+=====+--+  /
 *      |     |                              |     |  |  \
 *      | H0  |                              | HN  |Hq|  |
 *      |     |                              |     |  |  | n+qr
 *      +-----+-----+---......---+-----+-----+=====+--+  /
 *                                                  qx
 *  \endverbatim
 *
 *  \param mat_id           identifier for the factorization, used in the subsequent `ABD_solve`
 *  \param mat_fact         factorization of cyclic reduction intermediate blocks 0 - LU, 1 - QR, 2 - QRP
 *  \param last_block_fact  last block factorization type 0 - LU, 1 - QR, 2 - QRP, 3 - SVD
 *  \param nblock           number of blocks `D` and `E`
 *  \param n                dimension of the blocks `D` and `E` (size `n` x `n`)
 *  \param qr               number of extra equation
 *  \param qx               number of extra unknown
 *  \param DE               pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
 *                          The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
 *  \param ldDE             leading dimension of matrices `D` and `E`
 *  \param H0               pointer of the `H0` block stored by column (FORTRAN STORAGE)
 *  \param ldH0             leading dimension of matrix `H0`
 *  \param HN               pointer of the `HN` block stored by column (FORTRAN STORAGE)
 *  \param ldHN             leading dimension of matrix `HN`
 *  \param Hq               pointer of the `Hq` block stored by column (FORTRAN STORAGE)
 *  \param ldHq             leading dimension of matrix `Hq`
 *
 *  \return  0 no error found
 *
\*/

int
BABD_factorize(
  BABD_intType        mat_id,
  BABD_intType        mat_fact,
  BABD_intType        last_block_fact,
  BABD_intType        nblock,
  BABD_intType        n,
  BABD_intType        qr,
  BABD_intType        qx,
  BABD_realType const DE[], BABD_intType ldDE,
  BABD_realType const H0[], BABD_intType ldH0,
  BABD_realType const HN[], BABD_intType ldHN,
  BABD_realType const Hq[], BABD_intType ldHq
);

/*!
 *  Perform factorization of BABD matrix defined by blocks H, D and E.
 *  Matrix can be bordered.
 *
 *  \verbatim
 *  Matrix structure
 *
 *                 n * (nblock+1)
 *    ___________________^____________________
 *   /                                        \
 *    n   n   n                              n  qx  nx
 *  +---+---+---+----.................-----+---+---+---+   -+
 *  | D | E |   |                          |   |   | B | n  |
 *  +---+---+---+                     -----+---+---+---+    |
 *  |   | D | E |                          |   |   | B | n  |
 *  +---+---+---+---+                 -----+---+---+---+    |
 *  |   |   | D | E |                      |   |   | B | n  |
 *  +---+---+---+---+                 -----+---+---+---+    |
 *  :                                                  :    |
 *  :                                                  :    |
 *  :                                                  :     > n * nblock
 *  :                                                  :    |
 *  :                                                  :    |
 *  :                              +---+---+---+---+---+    |
 *  :                              | D | E |   |   | B | n  |
 *  :                              +---+---+---+---+---+    |
 *  :                                  | D | E |   | B | n  |
 *  +---+---+---................---+---+---+---+---+---+   -+
 *  |   |   |                          |   |   |   |   |    |
 *  |H0 | 0 |                          | 0 |HN | Hq| Bp|    | n+qr
 *  |   |   |                          |   |   |   |   |    |
 *  +---+---+---................---+---+---+---+---+---+   -+
 *  | C | C |                      | C | C | C | Cq| F |    | nr
 *  +---+---+---................---+---+---+---+---+---+   -+
 *                                             nr*qx
 *
 *  \endverbatim
 *
 *  \param mat_id           identifier for the factorization, used in the subsequent `ABD_solve`
 *  \param mat_fact         factorization of cyclic reduction intermediate blocks 0 - LU, 1 - QR, 2 - QRP
 *  \param last_block_fact  last block factorization type 0 - LU, 1 - QR, 2 - QRP, 3 - SVD
 *  \param nblock           number of blocks `D` and `E`
 *  \param n                dimension of the blocks `D` and `E` (size `n` x `n`)
 *  \param qr               integer
 *  \param nr               integer
 *  \param qx               integer
 *  \param nx               integer
 *  \param DE               pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
 *                          The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
 *  \param ldDE             leading dimension of matrices `D` and `E`
 *  \param H0               pointer of the `H0` block stored by column (FORTRAN STORAGE)
 *  \param ldH0             leading dimension of matrix `H0`
 *  \param HN               pointer of the `HN` block stored by column (FORTRAN STORAGE)
 *  \param ldHN             leading dimension of matrix `HN`
 *  \param Hq               pointer of the `Hq` block stored by column (FORTRAN STORAGE)
 *  \param ldHq             leading dimension of matrix `Hq`
 *  \param B                pointer of the `B` block stored by column (FORTRAN STORAGE)
 *  \param ldB              leading dimension of matrix `B`
 *  \param C                pointer of the `C` block stored by column (FORTRAN STORAGE)
 *  \param ldC              leading dimension of matrix `C`
 *  \param D                pointer of the `D` block stored by column (FORTRAN STORAGE)
 *  \param ldD              leading dimension of matrix `D`
 *
 *  \return 0 no error found
 *
\*/

int
BABD_factorize_bordered(
  BABD_intType        mat_id,
  BABD_intType        mat_fact,
  BABD_intType        last_block_fact,
  BABD_intType        nblock,
  BABD_intType        n,
  BABD_intType        qr,
  BABD_intType        nr,
  BABD_intType        qx,
  BABD_intType        nx,
  BABD_realType const DE[], BABD_intType ldDE,  /* n x (2*n*nblock) */
  BABD_realType const H0[], BABD_intType ldH0,  /* (n+qr) x n */
  BABD_realType const HN[], BABD_intType ldHN,  /* (n+qr) x n */
  BABD_realType const Hq[], BABD_intType ldHq,  /* (n+qr) x qx */
  BABD_realType const B[],  BABD_intType ldB,   /* (n*(nblock+1)+qr) x nx */
  BABD_realType const C[],  BABD_intType ldC,   /* nr x (n*(nblock+1)+qx) */
  BABD_realType const D[],  BABD_intType ldD    /* nr x nx */
);

/*!
 *  solve linear ABD system using factorization of `ABD_factorize` call
 *
 *  \param mat_id   identifier for the factorization
 *  \param rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
 *
 *  \return  0 no error found
 */

int
BABD_solve( BABD_intType mat_id, BABD_realType rhs_sol[] );

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
BABD_solve_nrhs(
  BABD_intType  mat_id,
  BABD_intType  nrhs,
  BABD_realType rhs_sol[],
  BABD_intType  ldRhs
);

/*!
 *  destroy (and free mempry) of a factorization of ABD matrix.
 *  \param mat_id identifier for the factorization
 *
 *  \return  0 no error found
 */

int
BABD_free( BABD_intType mat_id );

/*!
 *  \return pointer to a string with the last error found
 */

char const *
BABD_get_last_error(void);

void
BABD_get_last_error_f90( char [], long len );

#ifdef __cplusplus
}
#endif

#endif
