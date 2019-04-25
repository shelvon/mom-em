! MODULE: modes
! AUTHOR: Xiaorun ZANG
! DESCRIPTION:
! Routines to solve both linear and nonlinear wave scattering problems.
MODULE scattering
  USE sysmat
  USE bc
  USE aux
  USE time
  USE material
  USE nfpost

  IMPLICIT NONE

CONTAINS
  SUBROUTINE scattering_linear( domain, ga, media, source, zwl, base )
    IMPLICIT NONE
    TYPE(domain_type), DIMENSION(-1:), INTENT(IN) :: domain
    TYPE(group_action), DIMENSION(:), INTENT(IN)  :: ga
    TYPE(media_type), DIMENSION(0:), INTENT(IN)   :: media
    TYPE(source_type), DIMENSION(:), INTENT(IN)   :: source
    COMPLEX(KIND=dp), INTENT(IN)                  :: zwl
    TYPE(base_type), INTENT(INOUT)                :: base

    COMPLEX(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: A
    TYPE(group_action), DIMENSION(1:1)              :: ga_tmp
    COMPLEX(KIND=dp)    :: ri, eta, phdx, phdy
    REAL(KIND=dp)       :: wl
    INTEGER             :: nbasis, nga, nsrc, isrc, ndom

    nbasis = domain(-1)%elements%nedges ! The total number of edges
    !nga = SIZE(ga)
    ! Temporarily, no group action.
    ! The 1st one is just identity group action.
    nga = 1
    ga_tmp = ga(1)
    nsrc = SIZE(source)
    ndom = UBOUND(domain,1)
!    phdx = CMPLX(1, 0, KIND=dp)
!    phdy = CMPLX(1, 0, KIND=dp)

    ALLOCATE( A(nbasis*2,nbasis*2,nga) )

    CALL sysmat_pmchwt_linear(domain, ndom, ga_tmp, nga, media, zwl, A, nbasis)

    WRITE(*,*) 'Computing source vectors'
    CALL timer_start()

    wl = REAL(zwl, KIND=dp)
    ri = get_ri(media(domain(0)%media)%ri, wl)
    ALLOCATE( base%x(nbasis*2, nga, nsrc) )
    DO isrc = 1, nsrc
      CALL srcvec(domain(0)%elements, domain(0)%elements%nedges, wl, ri, &
        ga_tmp, source(isrc), domain(0)%qd_tri, base%x(:,:,isrc))
!         WRITE(*,*) base%x(:,:,isrc)
!         STOP
    END DO
    WRITE(*,*) sec_to_str(timer_end())

    ! Solve the linear system of equations for each representation,
    ! enforcing possible boundary conditions.
    ! This is O(N^3) but often in practice very fast.
    WRITE(*,*) 'Solving system'
    CALL timer_start()
    CALL cpu_timer_start()
    CALL solve_systems(domain(-1)%elements, ga_tmp, phdx, phdy, A, base%x)
    WRITE(*,*) 'Wall-clock time:'
    WRITE(*,*) sec_to_str(timer_end())
    WRITE(*,*) 'CPU time:'
    WRITE(*,*) sec_to_str(cpu_timer_end())

  END SUBROUTINE scattering_linear

END MODULE scattering
