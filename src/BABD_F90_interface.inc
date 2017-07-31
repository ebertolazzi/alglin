!use iso_c_binding, only: C_CHAR, C_NULL_CHAR

interface
	!  
	!      _    ____  ____
	!     / \  | __ )|  _ \
	!    / _ \ |  _ \| | | |
	!   / ___ \| |_) | |_| |
	!  /_/   \_\____/|____/
	!  
	!  
	! 
	! 
	!  perform Diaz factorization of ABD matrix defined by blocks TOP, BOTTOM, D and E.
	!  Matrix structure is the following:
	! 
	!           col0
	!        |       |
	!      / +-------+                         \
	!      | |  TOP  |                         | <-- row0
	!      | +--+----+----+                    |
	!      |    | D1 | E1 |                    | <- sizeBlock
	!      |    +----+----+----+               |
	!      |         | D2 | E2 |               |
	!      |         +----+----+----+          |
	!      |               ...........         |
	!      |              +----+----+----+     |
	!      |                   | DN | EN |     |
	!      |                   +----+----+---+ |
	!      |                        |        | | <-- rowN
	!      |                        | BOTTOM | |
	!      \                        +--------+ /
	!                               |        |
	!                                  colN
	! 
	!  - id       identifier for the factorization, used in the subsequent `ABD_solve`
	!  - row0     number of rows of the `TOP` block
	!  - col0     number of cols of the `TOP` block
	!  - TOP      pointer of the `TOP` block stored by column (FORTRAN STORAGE)
	!  - ldTOP    leading dimension of matrix `TOP`
	!  - nblock   number of blocks `D` and `E`
	!  - n        dimension of the blocks `D` and `E` (size `n` x `n`)
	!  - DE       pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
	!             The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
	!  - ldDE     leading dimension of matrices `D` and `E`
	!  - rowN     number of rows of the `BOTTOM` block
	!  - colN     number of cols of the `BOTTOM` block
	!  - BOTTOM   pointer of the `BOTTOM` block stored by column (FORTRAN STORAGE)
	!  - ldBOTTOM leading dimension of matrix `BOTTOM`
	! 
	!  return  0 no error found
	! 
  integer &
  function &
  ABD_factorize( id, row0, col0, TOP, ldTOP, nblock, n, &
  	             DE, ldDE, rowN, colN, BOTTOM, ldBOTTOM ) &
  							 bind(C, name="ABD_factorize")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id, row0, col0, ldTOP, nblock, n, ldDE, rowN, colN, ldBOTTOM
  	real(kind=C_DOUBLE) :: TOP(*), DE(*), BOTTOM(*)
  end function ABD_factorize
  !
  !  solve linear ABD system using factorization of `ABD_factorize` call
  ! 
  !  - id       identifier for the factorization
  !  - rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
  ! 
  !  return  0 no error found
  ! 
  integer &
  function &
  ABD_solve( id, rhs_sol ) bind(C, name="ABD_solve")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id
  	real(kind=C_DOUBLE) :: rhs_sol(*)
  end function ABD_solve
  ! 
  !  solve linear ABD system using factorization of `ABD_factorize` call
  ! 
  !  - id       identifier for the factorization
  !  - nrhs     number of rhs
  !  - rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
  !  - ldRhs    leadind dimension of `rhs_sol`
  ! 
  !  return 0 no error found
  ! 
  integer &
  function &
  ABD_solve_nrhs( id, nrhs, rhs_sol, ldRhs ) bind(C, name="ABD_solve_nrhs")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id, nrhs, ldRhs
  	real(kind=C_DOUBLE) :: rhs_sol(*)
  end function ABD_solve_nrhs
  !
  ! destroy (and free mempry) of a factorization of ABD matrix.
  ! \param mat_id identifier for the factorization
  !
  ! \return  0 no error found
  !
  integer &
  function &
  ABD_free( id ) bind(C, name="ABD_free")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id
  end function ABD_free
 	!
 	!
 	subroutine &
 	ABD_get_last_error_f90( err_message, len ) bind(C,name="ABD_get_last_error_f90")
  	use, intrinsic :: iso_c_binding
   	character :: err_message(*)
    integer(kind=C_INT), intent(in), value :: len
  end subroutine ABD_get_last_error_f90
	!  
	!   ____    _    ____  ____
	!  | __ )  / \  | __ )|  _ \
	!  |  _ \ / _ \ |  _ \| | | |
	!  | |_) / ___ \| |_) | |_| |
	!  |____/_/   \_\____/|____/
	!
  ! 
	!  Perform factorization of BABD matrix defined by blocks.
	!  Matrix structure
	!
	!                 n * (nblock+1)
  !    ___________________^____________________
  !   /                                        \
  !    n   n   n                              n  qx  nx
  !  +---+---+---+----.................-----+---+---+   -+
  !  | D | E |   |                          |   |   | n  |
  !  +---+---+---+                     -----+---+---+    |
  !  |   | D | E |                          |   |   | n  |
  !  +---+---+---+---+                 -----+---+---+    |
  !  |   |   | D | E |                      |   |   | n  |
  !  +---+---+---+---+                 -----+---+---+    |
  !  :                                                   |
  !  :                                                   |
  !  :                                                    > n * nblock
  !  :                                                   |
  !  :                                                   |
  !  :                              +---+---+---+---+    |
  !  :                              | D | E |   |   | n  |
  !  :                              +---+---+---+---+    |
  !  :                                  | D | E |   | n  |
  !  +---+---+---................---+---+---+---+---+   -+
  !  |   |   |                          |   |   |   |    |
  !  |H0 | 0 |                          | 0 |HN | Hq|    | n+qr
  !  |   |   |                          |   |   |   |    |
  !  +---+---+---................---+---+---+---+---+   -+
  !
  ! 
  !  - mat_id           identifier for the factorization, used in the subsequent `ABD_solve`
  !  - mat_fact         factorization of cyclic reduction intermediate blocks 0 - LU, 1 - QR, 2 - QRP
  !  - last_block_fact  last block factorization type 0 - LU, 1 - QR, 2 - QRP, 3 - SVD
  !  - nblock           number of blocks `D` and `E`
  !  - n                dimension of the blocks `D` and `E` (size `n` x `n`)
  !  - qr               extra variable
  !  - qx               extra equations
  !  - DE               pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
  !                     The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
  !  - ldDE             leading dimension of matrices `D` and `E`
  !  - H0               pointer of the `H0` block stored by column (FORTRAN STORAGE)
  !  - ldH0             leading dimension of matrix `H0`
  !  - HN               pointer of the `HN` block stored by column (FORTRAN STORAGE)
  !  - ldHN             leading dimension of matrix `HN`
  !  - Hq               pointer of the `Hq` block stored by column (FORTRAN STORAGE)
  !  - ldHq             leading dimension of matrix `Hq`
  !  - B                pointer of the `B` block stored by column (FORTRAN STORAGE)
  !  - ldB              leading dimension of matrix `B`
  !  - C                pointer of the `C` block stored by column (FORTRAN STORAGE)
  !  - ldC              leading dimension of matrix `C`
  !  - D                pointer of the `D` block stored by column (FORTRAN STORAGE)
  !  - ldD              leading dimension of matrix `D`
  !
  integer &
  function &
  BABD_factorize( id, mat_fact, last_block_fact, nblock, n, qr, qx, &
                  DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq ) &
  		            bind(C, name="BABD_factorize")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id, mat_fact, last_block_fact, &
		                                          nblock, n, qr, qx, ldDE, ldH0, ldHN, ldHq
  	real(kind=C_DOUBLE) :: DE(*), H0(*), HN(*), Hq(*)
  end function
	!
  ! 
	!  Perform factorization of BABD matrix defined by blocks.
	!  Matrix structure
	!
	!                 n * (nblock+1)
  !    ___________________^____________________
  !   /                                        \
  !    n   n   n                              n  qx  nx
  !  +---+---+---+----.................-----+---+---+---+   -+
  !  | D | E |   |                          |   |   | B | n  |
  !  +---+---+---+                     -----+---+---+---+    |
  !  |   | D | E |                          |   |   | B | n  |
  !  +---+---+---+---+                 -----+---+---+---+    |
  !  |   |   | D | E |                      |   |   | B | n  |
  !  +---+---+---+---+                 -----+---+---+---+    |
  !  :                                                  :    |
  !  :                                                  :    |
  !  :                                                  :     > n * nblock
  !  :                                                  :    |
  !  :                                                  :    |
  !  :                              +---+---+---+---+---+    |
  !  :                              | D | E |   |   | B | n  |
  !  :                              +---+---+---+---+---+    |
  !  :                                  | D | E |   | B | n  |
  !  +---+---+---................---+---+---+---+---+---+   -+
  !  |   |   |                          |   |   |   |   |    |
  !  |H0 | 0 |                          | 0 |HN | Hq| Bp|    | n+qr
  !  |   |   |                          |   |   |   |   |    |
  !  +---+---+---................---+---+---+---+---+---+   -+
  !  | C | C |                      | C | C | C | Cq| F |    | nr
  !  +---+---+---................---+---+---+---+---+---+   -+
  !                                             nr*qx
  ! 
  !  - mat_id           identifier for the factorization, used in the subsequent `ABD_solve`
  !  - mat_fact         factorization of cyclic reduction intermediate blocks 0 - LU, 1 - QR, 2 - QRP
  !  - last_block_fact  last block factorization type 0 - LU, 1 - QR, 2 - QRP, 3 - SVD
  !  - nblock           number of blocks `D` and `E`
  !  - n                dimension of the blocks `D` and `E` (size `n` x `n`)
  !  - qr               integer
  !  - nr               integer
  !  - qx               integer
  !  - nx               integer
  !  - DE               pointer to the blocks  `D` and `E` stored by column (FORTRAN STORAGE).
  !                     The blocks are ordered as [D1,E1,D2,E2,...,DN,EN]
  !  - ldDE             leading dimension of matrices `D` and `E`
  !  - H0               pointer of the `H0` block stored by column (FORTRAN STORAGE)
  !  - ldH0             leading dimension of matrix `H0`
  !  - HN               pointer of the `HN` block stored by column (FORTRAN STORAGE)
  !  - ldHN             leading dimension of matrix `HN`
  !  - Hq               pointer of the `Hq` block stored by column (FORTRAN STORAGE)
  !  - ldHq             leading dimension of matrix `Hq`
  !  - B                pointer of the `B` block stored by column (FORTRAN STORAGE)
  !  - ldB              leading dimension of matrix `B`
  !  - C                pointer of the `C` block stored by column (FORTRAN STORAGE)
  !  - ldC              leading dimension of matrix `C`
  !  - D                pointer of the `D` block stored by column (FORTRAN STORAGE)
  !  - ldD              leading dimension of matrix `D`
  !
  integer &
  function &
  BABD_factorize_bordered( id, mat_fact, last_block_fact, nblock, n, qr, nr, qx, nx, &
                           DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq, B, ldB, C, ldC, D, ldD ) &
  		                     bind(C, name="BABD_factorize_bordered")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id, mat_fact, last_block_fact, &
		                                          nblock, n, qr, qx, nr, nx, ldDE, ldH0, ldHN, &
																							ldHq, ldB, ldC, ldD
  	real(kind=C_DOUBLE) :: DE(*), H0(*), HN(*), Hq(*), B(*), C(*), D(*)
  end function
  !
  !  solve linear BABD system using factorization of `BABD_factorize` call
  ! 
  !  - id       identifier for the factorization
  !  - rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
  ! 
  !  return  0 no error found
  ! 
  integer &
  function &
  BABD_solve( id, rhs_sol ) bind(C, name="BABD_solve")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id
  	real(kind=C_DOUBLE) :: rhs_sol(*)
  end function BABD_solve
  ! 
  !  solve linear BABD system using factorization of `BABD_factorize` call
  ! 
  !  - id       identifier for the factorization
  !  - nrhs     number of rhs
  !  - rhs_sol  rhs (INPUT) and solution (OUTOUT) of linear system
  !  - ldRhs    leadind dimension of `rhs_sol`
  ! 
  !  return 0 no error found
  ! 
  integer &
  function &
  BABD_solve_nrhs( id, nrhs, rhs_sol, ldRhs ) bind(C, name="BABD_solve_nrhs")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id, nrhs, ldRhs
  	real(kind=C_DOUBLE) :: rhs_sol(*)
  end function BABD_solve_nrhs
  !
  ! destroy (and free mempry) of a factorization of ABD matrix.
  ! \param mat_id identifier for the factorization
  !
  ! \return  0 no error found
  !
  integer &
  function &
  BABD_free( id ) bind(C, name="BABD_free")
  	use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in), value :: id
  end function BABD_free
 	!
 	!
 	subroutine &
 	BABD_get_last_error_f90( err_message, len ) bind(C,name="BABD_get_last_error_f90")
  	use, intrinsic :: iso_c_binding
   	character :: err_message(*)
    integer(kind=C_INT), intent(in), value :: len
  end subroutine BABD_get_last_error_f90


end interface
