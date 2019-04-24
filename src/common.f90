! MODULE: common
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Implements various data aggregation types for the description of scattering
! problems, including media and domains. These form the basis for the high-level
! code that is used to define certain types of problems via the user interface.
MODULE common
  USE data
  USE source
  USE focal
  USE greenprd
  USE nlsurf
  USE nlbulk

  IMPLICIT NONE

CONTAINS
  SUBROUTINE batch_defaults(b)
    TYPE(batch), INTENT(INOUT) :: b

    b%name = 'unnamed'
    b%mesh_file = 'default.msh'
    b%scale = 1d-9
    b%nwl = 0
    b%src%type = 0

    b%qd_tri = tri_quad_data('tri_gl13')
    b%qd_tetra = tetra_quad_data('tetra_gl4')
  END SUBROUTINE batch_defaults


  ! Checks whether the batch contains nonlinear centrosymmetric materials.
  FUNCTION is_nl_centrosym(b) RESULT(res)
    TYPE(batch), INTENT(IN) :: b
    LOGICAL :: res
    INTEGER :: n

    res = .FALSE.

    DO n=1,SIZE(b%media)
       IF(b%media(n)%type==mtype_nls .OR. b%media(n)%type==mtype_nlb_nonlocal) THEN
          res = .TRUE.
          RETURN
       END IF
    END DO
  END FUNCTION is_nl_centrosym

  ! Checks whether the batch contains nonlinear non-centrosymmetric materials.
  FUNCTION is_nl_noncentrosym(b) RESULT(res)
    TYPE(batch), INTENT(IN) :: b
    LOGICAL :: res
    INTEGER :: n

    res = .FALSE.

    DO n=1,SIZE(b%media)
       IF(b%media(n)%type==mtype_nlb_dipole) THEN
          res = .TRUE.
          RETURN
       END IF
    END DO
  END FUNCTION is_nl_noncentrosym

  ! Checks whether the batch contains nonlinear materials for third harmonic generation (thg).
  FUNCTION is_nl_thg(b) RESULT(res)
    TYPE(batch), INTENT(IN) :: b
    LOGICAL :: res
    INTEGER :: n

    res = .FALSE.

    DO n=1,SIZE(b%media)
       IF(b%media(n)%type==mtype_thg_bulk) THEN
          res = .TRUE.
          RETURN
       END IF
    END DO
  END FUNCTION is_nl_thg

  SUBROUTINE delete_batch(b)
    TYPE(batch), INTENT(INOUT) :: b
    INTEGER :: i

    CALL delete_quad_data(b%qd_tri)
    CALL delete_quad_data(b%qd_tetra)

    CALL delete_mesh(b%mesh)

    IF(ALLOCATED(b%sols)) THEN
       DO i=1,SIZE(b%sols)
          CALL delete_solution(b%sols(i))
       END DO

       DEALLOCATE(b%sols)
    END IF

    IF(ALLOCATED(b%ga)) THEN
       DO i=1,SIZE(b%ga)
          DEALLOCATE(b%ga(i)%ef)
       END DO

       DEALLOCATE(b%ga)
    END IF

    IF(ALLOCATED(b%domains)) THEN
       DO i=1,SIZE(b%domains)
          CALL delete_domain(b%domains(i))
       END DO

       DEALLOCATE(b%domains)
    END IF

    IF(ALLOCATED(b%media)) THEN
       DEALLOCATE(b%media)
    END IF

    IF(ASSOCIATED(b%prd)) THEN
       DO i=1,SIZE(b%prd)
          CALL clear_prd(b%prd(i))
       END DO

       DEALLOCATE(b%prd)
    END IF

    IF(ALLOCATED(b%media)) THEN
       DO i=1,SIZE(b%media)
          DEALLOCATE(b%media(i)%prop)
       END DO

       DEALLOCATE(b%media)
    END IF

    IF(ALLOCATED(b%src)) THEN
       DEALLOCATE(b%src)
    END IF
  END SUBROUTINE delete_batch

  SUBROUTINE delete_solution(s)
    TYPE(solution), INTENT(INOUT) :: s

    IF(ALLOCATED(s%x)) THEN
       DEALLOCATE(s%x)
    END IF

    IF(ALLOCATED(s%nlx)) THEN
       DEALLOCATE(s%nlx)
    END IF

    IF(ALLOCATED(s%src_coef)) THEN
       DEALLOCATE(s%src_coef)
    END IF

    IF(ALLOCATED(s%eigvec)) THEN
       DEALLOCATE(s%eigvec)
    END IF

    IF(ALLOCATED(s%eigval)) THEN
       DEALLOCATE(s%eigval)
    END IF
  END SUBROUTINE delete_solution

  SUBROUTINE delete_domain(d)
    TYPE(domain), INTENT(INOUT) :: d

    CALL delete_mesh(d%mesh)
  END SUBROUTINE delete_domain
END MODULE common
