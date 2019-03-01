! MODULE: modes
! AUTHOR: Xiaorun ZANG
! DESCRIPTION:
! Routines to find intrinsic resonance modes of single nanoparticle.
MODULE modes 
  USE sysmat
  USE bc
  USE aux
  USE time
  USE material
  USE nfpost

  IMPLICIT NONE

CONTAINS

  ! near field of the mode on the given mesh surface
  SUBROUTINE modes_nfms( name, geom, ga, media, mode )
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: name
    TYPE(geom_type), INTENT(IN)         :: geom
    TYPE(group_action), DIMENSION(:), INTENT(IN)  :: ga
    TYPE(media_type), DIMENSION(0:), INTENT(IN)   :: media
    TYPE(mode_type), INTENT(IN)         :: mode

    INTEGER                             :: nbasis, imode
    COMPLEX(KIND=dp)                    :: zomega, ri, eps
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: x
    CHARACTER(LEN=256)                  :: filename

    nbasis = geom%domain(-1)%elements%nedges ! total number of edges
    ALLOCATE( x(2*nbasis, 1) )
    DO imode = 1, SIZE(mode%modeidx)
      filename = TRIM(name)//'-mode'//TRIM(num2str(imode))//'.msh'
      x(:,1) = mode%eigvec(:,mode%modeidx(imode))
      zomega = 2.0_dp*pi*c0/mode%zwl
      eps = media_eps( media(geom%domain(0)%media), mode%zwl )*eps0
      ri = SQRT(eps/eps0)
      CALL field_mesh( TRIM(filename), geom%domain(0)%elements, &
                       geom%mesh%scale, nbasis, x, ga, zomega, ri)
    END DO

  END SUBROUTINE modes_nfms

  SUBROUTINE pole_cauchy( domain, ga, media, zwl, mode )
    IMPLICIT NONE
    TYPE(domain_type), DIMENSION(-1:), INTENT(IN) :: domain
    TYPE(group_action), DIMENSION(:), INTENT(IN)  :: ga
    TYPE(media_type), DIMENSION(0:), INTENT(IN)   :: media
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)    :: zwl
    TYPE(mode_type), INTENT(INOUT)                :: mode

    COMPLEX(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: B, invB
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, C1, C2, C3, F, invF
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: U, VH
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: Sigma, invSigma
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigval, omega
    INTEGER, DIMENSION(:), ALLOCATABLE            :: eigidx
    INTEGER           :: m, k, nbasis, neig, nwl, iwl, r
    COMPLEX(KIND=dp)  :: ri, eta
    REAL(KIND=dp)     :: time_begin, time_end, time_loop

    nbasis = domain(-1)%elements%nedges ! The total number of edges
    neig = nbasis*2
    nwl = SIZE(zwl)

    ! Allocate memory for matrix/vectors, group representation is not used.
    ALLOCATE( omega(nwl) )
    ALLOCATE( B(1:nbasis*2,1:nbasis*2, 1:nwl), C3(1:nbasis*2,1:nbasis*2), &
              invB(1:nbasis*2,1:nbasis*2, 1:nwl), &
              F(1:nbasis,1:nbasis), invF(1:nbasis,1:nbasis),  &
              C1(1:nbasis*2,1:nbasis*2), C2(1:nbasis*2,1:nbasis*2), &
              U(1:nbasis*2,1:nbasis*2), VH(1:nbasis*2,1:nbasis*2), &
              Sigma(1:nbasis*2,1:nbasis*2), invSigma(1:nbasis*2,1:nbasis*2) )

    ! initialization
    B(:,:,:) = CMPLX(0.0_dp, 0.0_dp)
    invB(:,:,:) = CMPLX(0.0_dp, 0.0_dp)
    C1(:,:) = CMPLX(0.0_dp, 0.0_dp)
    C2(:,:) = CMPLX(0.0_dp, 0.0_dp)
    C3(:,:) = CMPLX(0.0_dp, 0.0_dp)
    U(:,:) = CMPLX(0.0_dp, 0.0_dp)
    VH(:,:) = CMPLX(0.0_dp, 0.0_dp)
    Sigma(:,:) = 0.0_dp
    invSigma(:,:) = 0.0_dp

    WRITE(*,*) 'Building matrices ...'
    WRITE(*,*) 'Triangle quadrature: ', TRIM(domain(-1)%quad(2)%qd%description)
    ! Compute the system matrix for all group representations.
    ! This is O(N^2), but in practice the most time consuming part.
    CALL rwg_moments(domain(-1)%elements, domain(-1)%quad(2)%qd, F)
    CALL matrix_inverse(F, invF)

    time_loop = 0.0_dp
    DO iwl = 1, nwl
      WRITE(*,*) '  -- wavelength: ', zwl(iwl)
      WRITE(*,*) '    ( ',iwl,' of ',nwl,' )'
      ! CALL timer_start()
      CALL CPU_TIME( time_begin )

      CALL sysmat_mueller( domain, ga, media, zwl(iwl), B(:,:,iwl), nbasis )
      ! Divide out the self-moment matrix on the RHS.
      CALL matmul_blocklt(invF, B(:,:,iwl), nbasis)
      DO m=1,(nbasis*2)
        B(m,m,iwl)= B(m,m,iwl) + 1.0_dp
      END DO
      B(1:nbasis,1:nbasis,iwl) = B(1:nbasis,1:nbasis,iwl) + F
      B(nbasis+1:2*nbasis,nbasis+1:2*nbasis,iwl) =  &
        B(nbasis+1:2*nbasis,nbasis+1:2*nbasis,iwl) + F

      ! inverse to get the scattering matrix form with a few poles
      ! CALL write_matrix( 'B-wl'//TRIM(num2str(iwl))//'.txt', B(:,:,iwl) )
      CALL matrix_inverse(B(:,:,iwl), invB(:,:,iwl))

      IF ( iwl >= 2 ) THEN
        omega(iwl) = 2.0_dp*pi*c0*( 1.0_dp/zwl(iwl) + 1.0_dp/zwl(iwl-1) )/2.0_dp
        C1 = C1 + (zwl(iwl)-zwl(iwl-1))/2.0_dp &
            * (invB(:,:,iwl) + invB(:,:,iwl-1))
        C2 = C2 + omega(iwl)*(zwl(iwl)-zwl(iwl-1))/2.0_dp &
            * (invB(:,:,iwl) + invB(:,:,iwl-1))
      END IF

      IF ( iwl == nwl ) THEN
        omega(iwl) = 2.0_dp*pi*c0*( 1.0_dp/zwl(1) + 1.0_dp/zwl(iwl) )/2.0_dp
        C1 = C1 + (zwl(1)-zwl(iwl))/2.0_dp &
            * (invB(:,:,1) + invB(:,:,iwl))
        C2 = C2 + omega(iwl)*(zwl(1)-zwl(iwl))/2.0_dp &
            * (invB(:,:,1) + invB(:,:,iwl))
      END IF

      CALL CPU_TIME( time_end )

      time_loop = time_loop + time_end-time_begin
      WRITE(*,*) 'Wall-clock time:'
      ! WRITE(*,*) sec_to_str(timer_end())
      WRITE(*,*) sec_to_str( time_end-time_begin )

      WRITE(*,*) 'Wall-clock time (loop):'
      WRITE(*,*) sec_to_str( time_loop )

    END DO
    C1 = C1/(2.0_dp*pi*(0,1.0_dp))
    C2 = C2/(2.0_dp*pi*(0,1.0_dp))

    ! SVD to C1 matrix
    ! CALL GESVD( C1, Sigma U, VH, [,ww] [,job] [,info])
    CALL zsvd( C1, U, Sigma, VH )
    DO r = 1, neig
      IF ( Sigma(r,r) > tol*1E5 ) THEN ! Sigma(r,r) must be real, positive
        C3(r,:) = MATMUL( CONJG(U(:,r)), C2(:,:) )

        C3(r,r) = dotc( C3(r,:), CONJG(VH(r,:)) )
        C3(r,r) = C3(r,r)/Sigma(r,r)
      ELSE
        EXIT
      END IF
    END DO
    r = r-1
    WRITE(*,*) 'r=',r
    ! CALL matrix_inverse(Sigma, invSigma)
    ! C3 = MATMUL( CONJG(TRANSPOSE(U)), C2 )
    ! C3 = MATMUL( C3, CONJG(TRANSPOSE(VH)) )
    ! C3 = MATMUL( C3, invSigma )

    ALLOCATE( A(r,r), eigval(r), eigvec(r,r), eigidx(r) )
    A(1:r,1:r) = C3(1:r,1:r)
    WRITE(*,*) 'Solving eigen values ...'
    CALL matrix_eigenvalues( A, eigval, eigvec )

    WRITE(*,*) 'Wall-clock time:'
    WRITE(*,*) sec_to_str(timer_end())
    WRITE(*,*) 'CPU time:'
    WRITE(*,*) sec_to_str(cpu_timer_end())

    ALLOCATE( mode%eigval(r) )
    ALLOCATE( mode%eigvec(r,r) )
    mode%eigval = eigval
    mode%eigvec = eigvec

    DEALLOCATE(A, B, invB, F, invF, omega, eigval, eigvec, eigidx)
    DEALLOCATE( C1, C2, C3, U, Sigma, invSigma, VH)
  END SUBROUTINE pole_cauchy

  SUBROUTINE modes_mueller( domain, ga, media, zwl, mode )
    IMPLICIT NONE
    TYPE(domain_type), DIMENSION(-1:), INTENT(IN) :: domain
    TYPE(group_action), DIMENSION(:), INTENT(IN)  :: ga
    TYPE(media_type), DIMENSION(0:), INTENT(IN)   :: media
    COMPLEX(KIND=dp), INTENT(IN)                  :: zwl
    TYPE(mode_type), INTENT(INOUT)                :: mode

    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, Aadj, F, invF
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigval
    INTEGER, DIMENSION(:), ALLOCATABLE            :: eigidx
    INTEGER           :: m, k, nbasis, neig
    COMPLEX(KIND=dp)  :: ri, eta, omega

    nbasis = domain(-1)%elements%nedges ! The total number of edges
    neig = nbasis*2

    ! Allocate memory for matrix/vectors, group representation is not used.
    ALLOCATE( mode%eigval(1:neig) )
    ALLOCATE( mode%eigvec(1:nbasis*2,1:neig) )
    ALLOCATE( A(1:nbasis*2,1:nbasis*2), Aadj(1:nbasis*2,1:nbasis*2),  &
              F(1:nbasis,1:nbasis), invF(1:nbasis,1:nbasis) )
    ALLOCATE( eigval(1:nbasis*2), eigidx(1:neig), &
              eigvec(1:nbasis*2,1:nbasis*2) )

    ! initialization
    A(:,:) = CMPLX(0.0_dp, 0.0_dp)
    Aadj(:,:) = CMPLX(0.0_dp, 0.0_dp)
    eigval(:) = CMPLX(0.0_dp, 0.0_dp)
    eigvec(:,:) = CMPLX(0.0_dp, 0.0_dp)

    WRITE(*,*) 'Building matrices ...'
    WRITE(*,*) 'Triangle quadrature: ', TRIM(domain(-1)%quad(2)%qd%description)

    CALL rwg_moments(domain(-1)%elements, domain(-1)%quad(2)%qd, F)
    CALL matrix_inverse(F, invF)

    ! Compute the system matrix for all group representations.
    ! This is O(N^2), but in practice the most time consuming part.
    CALL sysmat_mueller( domain, ga, media, zwl, A, nbasis )

    ! Divide out the self-moment matrix on the RHS.
    CALL matmul_blocklt( invF, A, nbasis )

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

    ! Solve the eigenvalue problem for each representation,
    ! enforcing possible boundary conditions.
    WRITE(*,*) 'Solving eigenvalues ...'
    ! Solve the eigenvalue problem of matrix A.
    ! No group representation.
    A(1:nbasis,1:nbasis) = A(1:nbasis,1:nbasis) + F
    A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) =  &
      A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) + F

    CALL matrix_eigenvalues(A, eigval, eigvec)

    ! eigidx = find_smallest(CMPLX((0.0_dp, AIMAG(eigval)), KIND=dp), nbasis*2, neig)
    ! eigidx = find_smallest(CMPLX(REAL(eigval, KIND=dp), KIND=dp), nbasis*2, neig)
    ! eigidx = find_smallest(eigval, nbasis*2, neig)
    eigidx(1:neig) = find_eigidx(eigval, nbasis*2, neig)
    mode%eigval(1:neig)=eigval(eigidx(1:neig))
    mode%eigvec(1:nbasis*2,1:neig)=eigvec(1:nbasis*2,eigidx(1:neig))

    ! Normalize modes to <fn,fn> = 1.
    DO k = 1 , nbasis*2
      CALL normalize_eigvec( mode%eigvec(1:nbasis*2,k), F, nbasis)
    END DO

    DEALLOCATE(A, Aadj, F, invF, eigval, eigvec, eigidx)
  END SUBROUTINE modes_mueller

  ! Calculate modes of a scatter in a range of fixed real angular frequencies
  ! as specified in input structure b.
  SUBROUTINE modes_mueller_b(b)
    TYPE(batch), INTENT(INOUT) :: b

    ! INTEGER, INTENT(IN) :: type
    INTEGER :: n, m, l, k, nbasis, nf, nind, neig
    REAL (KIND=dp) :: wl, omega
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: A, Aadj, F, invF
    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec, eigdata
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    INTEGER, DIMENSION(:), ALLOCATABLE :: eigidx

    COMPLEX (KIND=dp) :: ri, eta
    ! COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: epsp
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

    ALLOCATE(A(1:nbasis*2,1:nbasis*2), Aadj(1:nbasis*2,1:nbasis*2),&
            F(1:nbasis,1:nbasis), invF(1:nbasis,1:nbasis))
    ALLOCATE(eigval(1:nbasis*2), eigidx(1:neig), eigvec(1:nbasis*2,1:nbasis*2))
    ALLOCATE(eigdata(1:b%nwl,1:neig))

    WRITE(*,*) 'Name: ', TRIM(b%name)
    WRITE(*,*) 'Preparing to solve the eigenmodes (Mueller)'
    WRITE(*,*) '---'
    WRITE(*,*) '--- Begin wavelength batch ---'
    ! Go through all given wavelengths.

    CALL rwg_moments(b%mesh, b%qd_tri, F)
    CALL matrix_inverse(F, invF)

    ! DO n=b%nwl,1,-1
    DO n=1,b%nwl
       wl = b%sols(n)%wl

       ! initialize in each loop
       A(:,:) = CMPLX(0.0_dp, 0.0_dp)
       Aadj(:,:) = CMPLX(0.0_dp, 0.0_dp)
       eigval(:) = CMPLX(0.0_dp, 0.0_dp)
       eigvec(:,:) = CMPLX(0.0_dp, 0.0_dp)

       ! Allocate memory for eigenvectors, group representation is not used
       ALLOCATE(b%sols(n)%eigval(1:neig))
       ALLOCATE(b%sols(n)%eigvec(1:nbasis*2,1:neig))
       ALLOCATE(b%sols(n)%x(1:nbasis*2,1,1:neig))
       b%sols(n)%eigval(1:neig) = CMPLX(0.0_dp, 0.0_dp)
       b%sols(n)%eigvec(1:nbasis*2,1:neig) = CMPLX(0.0_dp, 0.0_dp)
       b%sols(n)%x(1:nbasis*2,1,1:neig) = CMPLX(0.0_dp, 0.0_dp)

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
       CALL sysmat_mueller_b(b, n, A, b%mesh%nedges)

!       ! Compute the adjoint.
!       Aadj = A
!       Aadj = TRANSPOSE(CONJG(Aadj))
!       CALL matmul_blockrt(Aadj, F, nbasis)
!       CALL matmul_blocklt(MATMUL(invF, invF), Aadj, nbasis)

       ! Divide out the self-moment matrix on the RHS.
       CALL matmul_blocklt(invF, A, nbasis)

       ! Eigenvalue problem for Fredholm operator I + A.
       ! Add the identity part.
!       IF(type==muller_fredholm) THEN
         DO m=1,(nbasis*2)
            A(m,m) = A(m,m) + 1.0_dp
             ! Aadj(m,m) = Aadj(m,m) + 1.0_dp
          END DO
       ! SVD for compact operator A.
!       ELSE IF(type==muller_svd) THEN
!         A = MATMUL(A, Aadj)
!       END IF

       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())


       ! Solve the eigenvalue problem for each representation,
       ! enforcing possible boundary conditions.
       WRITE(*,*) 'Solving eigenvalues ...'
       CALL timer_start()
       CALL cpu_timer_start()

       ! Solve the eigenvalue problem of matrix A.
       ! No group representation.
       A(1:nbasis,1:nbasis) = A(1:nbasis,1:nbasis) + F
       A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) = A(nbasis+1:2*nbasis,nbasis+1:2*nbasis) + F

       CALL matrix_eigenvalues(A, eigval, eigvec)

       ! eigidx = find_smallest(CMPLX((0.0_dp, AIMAG(eigval)), KIND=dp), nbasis*2, neig)
       ! eigidx = find_smallest(CMPLX(REAL(eigval, KIND=dp), KIND=dp), nbasis*2, neig)
       ! eigidx = find_smallest(eigval, nbasis*2, neig)
       eigidx(1:neig) = find_eigidx(eigval, nbasis*2, neig)
       b%sols(n)%eigval(1:neig)=eigval(eigidx(1:neig))
       b%sols(n)%eigvec(1:nbasis*2,1:neig)=eigvec(1:nbasis*2,eigidx(1:neig))
       ! b%sols(n)%x(:,1,:)=eigvec(:,eigidx(:))
       ! b%sols(n)%x(1:nbasis,1,:)=eigvec(1:nbasis,eigidx(:))/eta
       b%sols(n)%x(1:nbasis,1,1:neig)=eigvec(1:nbasis,eigidx(1:neig))
       b%sols(n)%x(nbasis+1:2*nbasis,1,1:neig)=eigvec(nbasis+1:2*nbasis,eigidx(1:neig))
       eigdata(n,1:neig) = eigval(eigidx(1:neig))

       ! Normalize modes to <fn,fn> = 1.
       DO k=1,nbasis*2
          CALL normalize_eigvec(b%sols(n)%eigvec(1:nbasis*2,k), F, nbasis)
       END DO

       WRITE(*,*) 'Wall-clock time:'
       WRITE(*,*) sec_to_str(timer_end())
       WRITE(*,*) 'CPU time:'
       WRITE(*,*) sec_to_str(cpu_timer_end())

       ! incrementally write result in each loop
       CALL write_data(TRIM(b%name) // '-real_eig.txt', REAL(eigdata))
       CALL write_data(TRIM(b%name) // '-imag_eig.txt', AIMAG(eigdata))
    END DO
    WRITE(*,*) '--- End wavelength batch ---'

    DEALLOCATE(A, Aadj, F, invF, eigval, eigvec)
  END SUBROUTINE modes_mueller_b

  FUNCTION find_eigidx(eigval, eigdim, neig) RESULT(ind)
    IMPLICIT NONE
    COMPLEX (KIND=dp), DIMENSION(1:neig), INTENT(IN) :: eigval
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
  END FUNCTION find_eigidx

  SUBROUTINE normalize_eigvec(v, F, eigdim)
    IMPLICIT NONE
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
