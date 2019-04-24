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

  ! prinit progress bar in command line
  SUBROUTINE show_progress(ido, ndo, nt_set)
    INTEGER, INTENT(IN)           :: ido, ndo
    INTEGER, OPTIONAL, INTENT(IN) :: nt_set! user specified number of frames

    INTEGER                   :: it, nt=20 ! update nt times the percentage in the progress bar
    CHARACTER(LEN=8)          :: pctstr ! convert number to string
    CHARACTER(LEN=64)         :: format_line
    CHARACTER(LEN=256)        :: progress_string
    CHARACTER(LEN=2)          :: end_string

    IF (PRESENT(nt_set)) THEN
      nt =  nt_set
    END IF
    it = ido*nt/ndo
    IF ( it > (ido-1)*nt/ndo ) THEN
      WRITE(pctstr,'(F6.2)') REAL(ido*100)/REAL(ndo)
      format_line = '(A,A1,A' // TRIM(ADJUSTL(pctstr)) // ',A1,F6.2,A1)'

      IF ( ndo >= ido) THEN
        IF ( ndo == ido ) THEN
          end_string = '\n'
        ELSE
          end_string = '\r'
        END IF

        progress_string = 'printf "%s' // end_string // '" "     ['// REPEAT('===', it) // REPEAT('   ', nt-it) // ']' &
                          // TRIM(ADJUSTL(pctstr)) // '%' // '"'
        ! WRITE(*,*) progress_string
        CALL SYSTEM(progress_string)
      END IF
    END IF
    RETURN
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
