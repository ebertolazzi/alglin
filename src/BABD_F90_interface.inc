!use iso_c_binding, only: C_CHAR, C_NULL_CHAR

interface
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
 
end interface
