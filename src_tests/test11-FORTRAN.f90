program main
  include 'BABD_F90_interface.inc'
  DOUBLE PRECISION DE, H0, HN, Hq, X, X1, B, C, D
  DOUBLE PRECISION ERROR, ERR
  DIMENSION DE(2,8), H0(3,2), HN(3,2), Hq(3,1), X(7), X1(8), B(7,1), C(1,7), D(1,1)

  !
	!  1 2 2  3
	!  2 4 1 -1
	!      1 2 2  3
	!      2 4 1 -1
	!  0 1     2  2 1
	!  1 0     0 -1 1
	!  1 2     2 -2 3
	!
  integer :: id, mat_fact, last_block_fact, nblock, n, qr, nr, qx, nx, ok
  character(LEN=1000) :: err_message
  !
  id              = 1
	mat_fact        = 0
	last_block_fact = 0
	nblock          = 2
	n               = 2
	qr              = 1
	qx              = 1
  ldDE            = 2
	ldH0            = 3
	ldHN            = 3
	ldHq            = 3
	DE = reshape((/10,2,2,40,2,1,3,-1,10,2,2,40,2,1,3,-1/), shape(DE))
	H0 = reshape((/0,1,1,1,0,2/), shape(H0))
	HN = reshape((/20,0,2,2,-10,2/), shape(HN))
	Hq = reshape((/1,1,30/), shape(Hq))
	X  = reshape((/32, 81, 66, 165, 121, -52, 237/), shape(X))

  print *, 'DE = ', DE
	
  ok = BABD_factorize( id, mat_fact, last_block_fact, nblock, n, qr, qx, &
                       DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq ) !, B, ldB, C, ldC, D, ldD ) ;
  print *, 'ok = ', ok
  ok = BABD_solve( 1, X ) ;
  !ok = BABD_solve_nrhs( id, 1, X, 12 ) ;
  print *, 'ok = ', ok
  print *, 'X = ', X
  call BABD_get_last_error_f90(err_message,len(err_message))
  print *,'err=',err_message

  ldB = 7
  ldC = 1
  ldD = 1
	nr  = 1
	nx  = 1
	
	B = reshape((/1,1,1,1,1,1,1/), shape(B))
	C = reshape((/1,-1,1,-1,1,-1,1/), shape(C))
	D = reshape((/-1/), shape(D))

  ok = BABD_factorize_bordered( id, mat_fact, last_block_fact, nblock, n, qr, nr, qx, nx, &
                                DE, ldDE, H0, ldH0, HN, ldHN, Hq, ldHq, &
																B, ldB, C, ldC, D, ldD ) ;
  print *, 'ok = ', ok
	X1  = reshape((/32, 81, 66, 165, 121, -52, 237, 4/), shape(X1))
  ok = BABD_solve( 1, X1 ) ;
  !ok = BABD_solve_nrhs( id, 1, X, 12 ) ;
  print *, 'ok = ', ok
  print *, 'X1 = ', X1
  call BABD_get_last_error_f90(err_message,len(err_message))
  print *,'err=',err_message
  print *, 'All done folks!'


end
