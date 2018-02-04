! MODULE: modes
! AUTHOR: Xiaorun ZANG
! DESCRIPTION:
! Routines to find intrinsic resonance modes of single nanoparticle.
MODULE modes 
  USE sysmat
  USE bc
  USE aux
  USE time

  IMPLICIT NONE

CONTAINS
  ! Calculate modes of nanoparticles for a range of fixed real angular frequencies
  ! as specified in input structure b.
  SUBROUTINE modes_mueller(b)
    TYPE(batch), INTENT(INOUT) :: b
    ! INTEGER, INTENT(IN) :: type
    INTEGER :: n, m, l, k, nbasis, nf, nind, neig
    REAL (KIND=dp) :: wl, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, Aadj, F, invF
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec, eigdata
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    INTEGER, DIMENSION(:), ALLOCATABLE :: eigind

    COMPLEX (KIND=dp) :: ri, eta
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: epsp
    INTEGER, DIMENSION(b%mesh%nedges) :: ind
    TYPE(medium_prop) :: mprop

    ! Check that necessary input data exists.
    IF(ALLOCATED(b%sols)==.FALSE.) THEN
       WRITE(*,*) 'Setup wavelengths prior to solving!'
       RETURN
    END IF

    IF(ALLOCATED(b%mesh%nodes)==.FALSE.) THEN
       WRITE(*,*) 'Load mesh prior to solving!'
       RETURN
    END IF

!    IF(ALLOCATED(b%domains)==.FALSE.) THEN
!       WRITE(*,*) 'Set up domains prior to solving!'
!       RETURN
!    END IF

    IF(ALLOCATED(b%media)==.FALSE.) THEN
       WRITE(*,*) 'Set up media prior to solving!'
       RETURN
    END IF

!    IF(SIZE(b%domains)/=2) THEN
!       WRITE(*,*) 'Only single scatter is supported!'
!       RETURN
!    END IF

    nbasis = b%mesh%nedges
    neig = nbasis*2

    ALLOCATE(A(nbasis*2,nbasis*2), Aadj(nbasis*2,nbasis*2),&
            F(nbasis,nbasis), invF(nbasis,nbasis))
    ALLOCATE(eigval(nbasis*2), eigind(neig), eigvec(nbasis*2,nbasis*2))
    ALLOCATE(eigdata(neig,b%nwl))

    WRITE(*,*) 'Name: ', TRIM(b%name)
    WRITE(*,*) 'Preparing to solve the eigenmodes (Mueller)'
    WRITE(*,*) '---'
    WRITE(*,*) '--- Begin wavelength batch ---'
    ! Go through all given wavelengths.

    CALL rwg_moments(b%mesh, b%qd_tri, F)
    CALL matrix_inverse(F, invF)

    DO n=1, b%nwl
       wl = b%sols(n)%wl

       ! Print some information.
       WRITE(*,'(A,F6.1,A,I0,A,I0,A)') ' Wavelength: ', wl*1d9, ' nm (', n, ' of ', b%nwl, ')'

       ! 1,2 two media for a single scatter
       DO m=1,2
          ri = b%media(m)%prop(n)%ri
          WRITE(*,'(A,I0,A,"(",F6.3,",",F6.3,")")') ' Refractive index of medium ', m, ': ', ri
       END DO
       ! eta=eta0/b%media(b%domains(1)%medium_index)%prop(n)%ri

       WRITE(*,*) 'Building matrices ...'
       WRITE(*,*) 'Triangle quadrature: ', TRIM(b%qd_tri%description)

       ! Compute the system matrix for all group representations.
       ! This is O(N^2), but in practice the most time consuming part.
       CALL timer_start()
       CALL cpu_timer_start()
       CALL sysmat_mueller(b, n, A, b%mesh%nedges)

       ! Compute the adjoint.
!       Aadj = A
!       Aadj = TRANSPOSE(CONJG(Aadj))
!       CALL matmul_blockrt(Aadj, F, nbasis)
!       CALL matmul_blocklt(MATMUL(invF, invF), Aadj, nbasis)

       ! Divide out the self-moment matrix on the RHS.
       CALL matmul_blocklt(invF, A, nbasis)

       ! Eigenvalue problem for Fredholm operator I + A.
       ! Add the identity part.
       ! IF(type==muller_fredholm) THEN
          DO m=1,(nbasis*2)
             A(m,m) = A(m,m) + 1.0_dp

             ! Aadj(m,m) = Aadj(m,m) + 1.0_dp
          END DO
          ! SVD for compact operator A.
       ! ELSE IF(type==muller_svd) THEN
       !   A = MATMUL(A, Aadj)
       ! END IF

       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())


       ! Solve the eigenvalue problem for each representation,
       ! enforcing possible boundary conditions.
       WRITE(*,*) 'Solving eigenvalues ...'
       CALL timer_start()
       CALL cpu_timer_start()

       ! Allocate memory for eigenvectors, group representation is not used
       ALLOCATE(b%sols(n)%eigval(neig))
       ALLOCATE(b%sols(n)%eigvec(nbasis*2,neig))
       ALLOCATE(b%sols(n)%x(nbasis*2,1,neig))

       ! Solve the eigenvalue problem of matrix A.
       ! No group representation.
       A(1:nbasis,1:nbasis) = A(1:nbasis,1:nbasis) + F
       A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) = A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) + F
       CALL matrix_eigenvalues(A, eigval, eigvec)

       ! eigind = find_smallest(CMPLX((0.0_dp, AIMAG(eigval)), KIND=dp), nbasis*2, neig)
       ! eigind = find_smallest(CMPLX(REAL(eigval, KIND=dp), KIND=dp), nbasis*2, neig)
       ! eigind = find_smallest(eigval, nbasis*2, neig)
       eigind = find_eigind(eigval, nbasis*2, neig)
       b%sols(n)%eigval(:)=eigval(eigind(:))
       b%sols(n)%eigvec(:,:)=eigvec(:,eigind(:))
       ! b%sols(n)%x(:,1,:)=eigvec(:,eigind(:))
       ! b%sols(n)%x(1:nbasis,1,:)=eigvec(1:nbasis,eigind(:))/eta
       b%sols(n)%x(1:nbasis,1,:)=eigvec(1:nbasis,eigind(:))
       b%sols(n)%x(nbasis+1:2*nbasis,1,:)=eigvec(nbasis+1:2*nbasis,eigind(:))
       eigdata(:,n) = eigval(eigind(:))

       ! Normalize modes to <fn,fn> = 1.
       DO k=1,nbasis*2
          CALL normalize_eigvec(b%sols(n)%eigvec(:,k), F, nbasis)
       END DO

       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

    END DO
    WRITE(*,*) '--- End wavelength batch ---'

    CALL write_data(TRIM(b%name) // '-real_eig.txt', REAL(eigdata))
    CALL write_data(TRIM(b%name) // '-imag_eig.txt', AIMAG(eigdata))

    DEALLOCATE(A, Aadj, F, invF, eigval, eigvec)
  END SUBROUTINE modes_mueller

  FUNCTION find_eigind(eigval, eigdim, neig) RESULT(ind)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: eigval
    INTEGER, INTENT(IN) :: eigdim, neig

    INTEGER, DIMENSION(neig) :: ind
    COMPLEX (KIND=dp), DIMENSION(eigdim) :: tmpeigval
    INTEGER :: n
    INTEGER, DIMENSION(1) :: loc
    REAL (KIND=dp) :: maxeig

    tmpeigval(:) = eigval(:)

    maxeig = MAXVAL(ABS(eigval))

    DO n=1,neig
       loc = MINLOC(ABS(tmpeigval))
       tmpeigval(loc(1)) = maxeig

       ind(n) = loc(1)
    END DO
  END FUNCTION find_eigind

  SUBROUTINE normalize_eigvec(v, F, eigdim)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(INOUT) :: v
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: F
    INTEGER, INTENT(IN) :: eigdim
    REAL (KIND=dp) :: fnorm
    INTEGER :: m, l

    fnorm = 0.0_dp

    DO m=1,eigdim
       DO l=(m+1),eigdim
          fnorm = fnorm + 2*F(m,l)*REAL(CONJG(v(m))*v(l) +&
               CONJG(v(m + eigdim))*v(l + eigdim))
       END DO

       fnorm = fnorm + F(m,m)*(ABS(v(m))**2 + ABS(v(m + eigdim))**2)
    END DO

    v = v/SQRT(fnorm)

  END SUBROUTINE normalize_eigvec

END MODULE modes
