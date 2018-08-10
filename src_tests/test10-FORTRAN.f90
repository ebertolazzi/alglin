program main
  include 'BABD_F90_interface.inc'
  !     ALGORITHM 603, COLLECTED ALGORITHMS FROM ACM.
  !     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 3,
  !     SEP., 1983, P. 376-380.
  !***********************************************************************
  !
  !   THIS DRIVER RUNS A TEST PROBLEM USING C O L R O W.
  !   THE ARRAYS TOP, AR, BOT, AND B ARE SET UP IN DATA
  !   STATEMENTS, AS ARE THE PARAMETERS NRWTOP, NRWBLK, NCLBLK,
  !   NBLOKS, NRWBOT, WHICH DEFINE THE STRUCTURE OF THE COEFF-
  !   ICIENT MATRIX.
  !
  !***********************************************************************
  !
  DOUBLE PRECISION TOP, AR, BOT, B, X
  DOUBLE PRECISION ERROR, ERR
  DIMENSION TOP(2,4), AR(4,8,2), BOT(2,4), B(12), X(12)
  DATA TOP(1,1), TOP(1,2), TOP(1,3), TOP(1,4), TOP(2,1), TOP(2,2), TOP(2,3), TOP(2,4) &
       /0.0D0,-0.98D0,-0.79D0,-0.15D0,-1.00D0,0.25D0,-0.87D0,0.35D0/
  DATA AR(1,1,1), AR(1,2,1), AR(1,3,1), AR(1,4,1), AR(1,5,1), AR(1,6,1), AR(1,7,1), AR(1,8,1) &
       /0.78D0,0.31D0,-0.85D0,0.89D0,-0.69D0,-0.98D0,-0.76D0,-0.82D0/
  DATA AR(2,1,1), AR(2,2,1), AR(2,3,1), AR(2,4,1), AR(2,5,1), AR(2,6,1), AR(2,7,1), AR(2,8,1) &
       /0.12D0,-0.01D0,0.75D0,0.32D0,-1.00D0,-0.53D0,-0.83D0,-0.98D0/
  DATA AR(3,1,1), AR(3,2,1), AR(3,3,1), AR(3,4,1), AR(3,5,1), AR(3,6,1), AR(3,7,1), AR(3,8,1) &
       /-0.58D0,0.04D0,0.87D0,0.38D0,-1.00D0,-0.21D0,-0.93D0,-0.84D0/
  DATA AR(4,1,1), AR(4,2,1), AR(4,3,1), AR(4,4,1), AR(4,5,1), AR(4,6,1), AR(4,7,1), AR(4,8,1) &
       /-0.21D0,-0.91D0,-0.09D0,-0.62D0,-1.99D0,-1.12D0,-1.21D0,0.07D0/
  DATA AR(1,1,2), AR(1,2,2), AR(1,3,2), AR(1,4,2), AR(1,5,2), AR(1,6,2), AR(1,7,2), AR(1,8,2) &
       /0.78D0,-0.93D0,-0.76D0,0.48D0,-0.87D0,-0.14D0,-1.00D0,-0.59D0/
  DATA AR(2,1,2), AR(2,2,2), AR(2,3,2), AR(2,4,2), AR(2,5,2), AR(2,6,2), AR(2,7,2), AR(2,8,2) &
       /-0.99D0,0.21D0,-0.73D0,-0.48D0,-0.93D0,-0.91D0,0.10D0,-0.89D0/
  DATA AR(3,1,2), AR(3,2,2), AR(3,3,2), AR(3,4,2), AR(3,5,2), AR(3,6,2), AR(3,7,2), AR(3,8,2) &
       /-0.68D0,-0.09D0,-0.58D0,-0.21D0,0.85D0,-0.39D0,0.79D0,-0.71D0/
  DATA AR(4,1,2), AR(4,2,2), AR(4,3,2), AR(4,4,2), AR(4,5,2), AR(4,6,2), AR(4,7,2), AR(4,8,2) &
       /0.39D0,-0.99D0,-0.12D0,-0.75D0,-0.68D0,-0.99D0,0.50D0,-0.88D0/
  DATA BOT(1,1), BOT(1,2), BOT(1,3), BOT(1,4), BOT(2,1), BOT(2,2), BOT(2,3), BOT(2,4) &
       /0.71D0,-0.64D0,0.0D0,0.48D0,0.08D0,100.0D0,50.00D0,15.00D0/
  DATA B(1), B(2), B(3), B(4), B(5), B(6), B(7), B(8), B(9), B(10), B(11), B(12) &
       /-1.92D0,-1.27D0,-2.12D0,-2.16D0,-2.27D0,-6.08D0, -3.03D0,-4.62D0,-1.02D0,-3.52D0,.55D0,165.08D0/
  DATA X /1,1,1,1,1,1,1,1,1,1,1,1/

  integer :: row0, col0, ldTOP, nblock, n, ldDE, rowN, colN, ldBOTTOM, ok
  character(LEN=1000) :: err_message

  !
  !***********************************************************************
  !
  !   THE INPUT MATRIX IS GIVEN BY:
  !
  !  0.0  -0.98 -0.79 -0.15
  ! -1.00  0.25 -0.87  0.35
  !  0.78  0.31 -0.85  0.89 -0.69 -0.98 -0.76 -0.82
  !  0.12 -0.01  0.75  0.32 -1.00 -0.53 -0.83 -0.98
  ! -0.58  0.04  0.87  0.38 -1.00 -0.21 -0.93 -0.84
  ! -0.21 -0.91 -0.09 -0.62 -1.99 -1.12 -1.21  0.07
  !                          0.78 -0.93 -0.76  0.48 -0.87 -0.14 -1.00 -0.5
  !                         -0.99  0.21 -0.73 -0.48 -0.93 -0.91  0.10 -0.8
  !                         -0.68 -0.09 -0.58 -0.21  0.85 -0.39  0.79 -0.7
  !                          0.39 -0.99 -0.12 -0.75 -0.68 -0.99  0.50 -0.8
  !                                                  0.71 -0.64  0.0   0.4
  !                                                  0.08 100.0 50.00 15.0
  !
  !       THE RIGHT HAND SIDE IS GIVEN BY:
  !
  !         B = (-1.92,-1.27,-2.12,-2.16,-2.27,-6.08,-3.03,-4.62,
  !              -1.02,-3.52,0.55,165.08)
  !
  !       THE SOLUTION OF THIS SYSTEM IS GIVEN BY;
  !
  !          X = (1,1,1,1,1,1,1,1,1,1,1,1)
  !
  !***********************************************************************
  !
  row0     = 2
  col0     = 4
  ldTOP    = 2
  nblock   = 2
  n        = 4
  ldDE     = 4
  rowN     = 2
  colN     = 4
  ldBOTTOM = 2

  ok = ABD_factorize( 1, row0, col0, TOP, ldTOP, nblock, n, AR, ldDE, rowN, colN, BOT, ldBOTTOM );
  print *, 'ok = ', ok
  !ok = ABD_solve( 1, B );
  ok = ABD_solve_nrhs( 1, 1, B, 12 );
  print *, 'ok = ', ok
  print *, 'B = ', B
  print *, 'X = ', X
  print *, 'All done folks!'
  call ABD_get_last_error_f90(err_message,len(err_message))
  print *,'err=',err_message
end
