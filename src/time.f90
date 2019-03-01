! MODULE: time
! AUTHOR: Jouni Makitalo
! Modified: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines for measuring time and time intervals.
MODULE time
  USE constants
  IMPLICIT NONE

  REAL (KIND=r8), PRIVATE, SAVE :: time_start, time_start_cpu

CONTAINS
  SUBROUTINE show_progress(p, i)
  ! SUBROUTINE show_progress(p, n, i)
  ! prinit command line progress bar
    REAL(KIND=dp), INTENT(IN) :: p ! finished percentage, belongs to [0, 1]
    ! INTEGER, INTENT(IN)       :: n ! show n times
    INTEGER, PARAMETER        :: n = 50 ! show 50 times by default
    INTEGER, INTENT(INOUT)    :: i ! a counter

    INTEGER                   :: ierr
    INTEGER                   :: nc ! total number of the characters
    CHARACTER(LEN=8)          :: nstr, ncstr, pstr ! convert number to string
    CHARACTER(LEN=64)         :: format_line
    CHARACTER(LEN=256)        :: progress_string
    CHARACTER(LEN=2)          :: end_string

    nc = n + 2 + 5 + 1 ! n: =/# symbols, 2: [], 5: ff.ff, 1: %
    WRITE(nstr,'(I0)') n
    WRITE(ncstr,'(I0)') nc
    WRITE(pstr,'(F6.2)') p*100
    format_line = '(A,A1,A' // TRIM(ADJUSTL(nstr)) // ',A1,F6.2,A1)'

    IF (FLOOR(p*n) >= i) THEN
      IF ( p == 1.0_dp) THEN
        end_string = '\n'
      ELSE
        end_string = '\r'
      END IF
      !
      ! method 1
      ! progress_string = REPEAT('=', i) // REPEAT(' ', n-i)
      ! WRITE(*, TRIM(ADJUSTL(format_line)),  ADVANCE='NO') 'In progress:','[',progress_string,']',p*100.0_dp, '%'
      ! WRITE(*, TRIM(ADJUSTL(format_line))) 'In progress:','[',progress_string,']',p*100.0_dp, '%'
      !
      ! method 2
      progress_string = 'printf "%s' // end_string // '" " ['// REPEAT('=', i) // REPEAT(' ', n-i) // ']' &
                        // TRIM(ADJUSTL(pstr)) // '%' // '"'
      ! WRITE(*,*) progress_string
      CALL SYSTEM(progress_string)
      i = i +1
    END IF
  END SUBROUTINE show_progress

  FUNCTION get_time() RESULT(sec)
    IMPLICIT NONE
    INTEGER, DIMENSION(8) :: t
    REAL (KIND=r8) :: sec

    CALL DATE_AND_TIME(values=t)
    sec = t(3)*86400.0_r8+t(5)*3600.0_r8+t(6)*60.0_r8+t(7)+t(8)/1000.0_r8
  END FUNCTION get_time

  SUBROUTINE timer_start()
    time_start = get_time()
  END SUBROUTINE timer_start

  FUNCTION timer_end() RESULT(sec)
    IMPLICIT NONE
    REAL (KIND=r8) :: sec
    sec = get_time() - time_start
  END FUNCTION timer_end

  FUNCTION sec_to_str(sec) RESULT(str)
    REAL (KIND=r8), INTENT(IN) :: sec
    CHARACTER (LEN=128) :: str
    INTEGER :: hours, minutes, seconds

    hours = FLOOR(sec/3600.0_r8)
    minutes = FLOOR((sec-hours*3600.0_r8)/60.0_r8)
    seconds = FLOOR((sec-hours*3600.0_r8-minutes*60.0_r8))
    write(str, '(I0,A,I0,A,I0,A)') hours, ' hours ', minutes, ' minutes ',&
         seconds, ' seconds '
  END FUNCTION sec_to_str

  SUBROUTINE cpu_timer_start()
    CALL cpu_time(time_start_cpu)
  END SUBROUTINE cpu_timer_start

  FUNCTION cpu_timer_end() RESULT(sec)
    IMPLICIT NONE
    REAL (KIND=r8) :: sec

    CALL cpu_time(sec)
    sec = sec - time_start_cpu
  END FUNCTION cpu_timer_end

END MODULE time
