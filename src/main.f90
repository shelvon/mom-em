PROGRAM main
  USE interface
  USE io
  USE common
  USE data
#ifdef _OPENMP
  USE omp_lib
#endif
  IMPLICIT NONE

  TYPE(model_type)    :: model
  INTEGER             :: max_threads
  INTEGER             :: narg, iarg
  CHARACTER (LEN=256) :: varg

  WRITE(*,*) 'Method of Moments solver'
  WRITE(*,*) 'Laboratory of Photonics, Tampere University of Technology'
  WRITE(*,*) '---'

  max_threads = OMP_GET_MAX_THREADS()
  CALL OMP_SET_NUM_THREADS(max_threads)
  WRITE(*,'(A,I0,:)') ' Number of threads for OpenMP: ', max_threads

  ! read parameters
  narg = IARGC()
  IF (narg .EQ. 0) THEN
    ! method 1: mom < input
    WRITE(*,*) 'Read modelling parameters through input redirection in command line.'
    CALL msgloop()
  ELSE IF (narg .GT. 0) THEN
    ! method 2: mom input.json
    ! DO iarg = 1, narg
      iarg = 1
      CALL GETARG(iarg, varg)

      WRITE(*,*) 'Read modelling parameters through JSON file "', TRIM(varg), '".'
      CALL json2model(TRIM(varg), model)

      CALL solve(model)
    ! END DO
  END IF
END PROGRAM main
