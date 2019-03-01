! MODULE: focal
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines to focal field definition and calculation.
MODULE focal
  USE data
  USE pupil
  USE quad
  USE aux
  USE linalg
  USE time
  USE bessel
  ! USE hankel_api

  IMPLICIT NONE

CONTAINS

  ! determine the equivalent pupil function for the focused beam
  FUNCTION focus2pupil(focus) RESULT(pupil)
    TYPE(focus_type), INTENT(IN)      :: focus
    TYPE(pupil_type)                  :: pupil

    ! initialize pupil_type::pupil
    pupil%type = 'analytical'
    pupil%aperture%type = 'circ'
    pupil%aperture%circ%r = 1
    pupil%phase%type = focus%type
    pupil%phase%input = 'hg00'
    pupil%phase%thetaRot=0.0_dp ! counter clockwise rotation
    pupil%phase%vJones(1)=CMPLX(1.0_dp, 0.0_dp)
    pupil%phase%vJones(2)=CMPLX(0.0_dp, 0.0_dp)
    pupil%phase%nMax = 1

    ! modification for each focus beam
    IF ( focus%type == 'hg10' .OR.  &
         focus%type == 'hg01' ) THEN

    ELSE IF ( focus%type == 'hg00' ) THEN
      pupil%phase%nMax = 0

    ELSE IF ( focus%type == 'rad' ) THEN
      ! rad = hg10*nx + hg10*ny, local SOP in radial direction
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

    ELSE IF ( focus%type == 'azi' ) THEN
      ! azi = hg01*nx + hg01*ny, local SOP in azimuthal direction, clockwise
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

    ELSE IF ( focus%type == 'm45' ) THEN ! m45 stands for minus 45
      ! general CVB, local SOP in the direction of -45 [deg], minus 45
      ! with respect to the positive radial axis
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

    ELSE IF ( focus%type == 'p45' ) THEN! p45 stands for plus 45
      ! general CVB, local SOP in the direction of +45 [deg]
      ! with respect to the positive radial axis
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

    ELSE IF ( focus%type == 'rcp' ) THEN! rcp stands for right-circular pol.
      ! general CVB, local SOP is right-circular polarization
      ! with respect to the positive radial axis
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

    ELSE IF ( focus%type == 'lcp' ) THEN! lcp stands for left-circular pol.
      ! general CVB, local SOP is left-circular polarization
      ! with respect to the positive radial axis
      pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
    ELSE IF ( focus%type == 'bessel' ) THEN
      pupil%phase%bessel = focus%bessel

      IF ( TRIM(pupil%phase%bessel%input) == 'rad' ) THEN
        ! rad = hg10*nx + hg10*ny
        pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
        pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)

      ELSE IF ( TRIM(pupil%phase%bessel%input) == 'azi' ) THEN
        ! azi = hg01*nx + hg01*ny, clockwise
        pupil%phase%vJones(1) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
        pupil%phase%vJones(2) = CMPLX(1.0_dp/SQRT(2.0_dp), 0.0_dp)
      END IF

    ELSE
      WRITE(*,*) 'The focused beam is not supported!'
    END IF

    ! get expansion coefficients
    ALLOCATE( pupil%cn(-pupil%phase%nMax:pupil%phase%nMax) )
    CALL pupil_coeff(pupil)
  END FUNCTION focus2pupil

  SUBROUTINE focal_field(focus, wl, ri1, ri2, grid, field)
    TYPE(focus_type), INTENT(IN)            :: focus
    REAL(KIND=dp), INTENT(IN)               :: wl
    COMPLEX(KIND=dp), INTENT(IN)            :: ri1, ri2
    TYPE(grid3d_type), INTENT(INOUT)        :: grid
    TYPE(field3d_type), INTENT(INOUT)       :: field

    TYPE(pupil_type)        :: pupil0 ! equivalent pupil function
    REAL(KIND=dp)           :: E0, omega, k, f, w0, na, theta_max, f0
    COMPLEX(KIND=dp)        :: Cprefix
    REAL(KIND=dp)           :: x, y, z, rho, varphi, varphiRot
    INTEGER                 :: ix, iy, iz, ibar
    INTEGER                 :: j, s, c, n, nMax, nbasis
    CHARACTER(LEN=32)       :: phase_type
    CHARACTER(LEN=12)       :: polyName ! name of basis function
    REAL(KIND=dp), DIMENSION(1:2)     :: thetaRot
    COMPLEX(KIND=dp), DIMENSION(1:2)  :: vJones

    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)   ::  &
    intJjS0C0, intJjS0C1, intJjS1C0, intJjS00C01
    LOGICAL, ALLOCATABLE, DIMENSION(:)            ::  &
    intJjS0C0_done, intJjS0C1_done, intJjS1C0_done, intJjS00C01_done
    LOGICAL, ALLOCATABLE, DIMENSION(:)            ::  &
    intJjS0C0_need, intJjS0C1_need, intJjS1C0_need, intJjS00C01_need
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: En
    COMPLEX(KIND=dp), DIMENSION(1:3)              :: Etilde

    ! determine pupil function based on focused beam type
    pupil0 = focus2pupil(focus)
    nMax = pupil0%phase%nMax
    nbasis = focus%bessel%nbasis
    polyName = focus%bessel%polyName

    E0 = 1.0_dp ! amplitude of the fundamental mode, hg00
    omega = 2.0_dp*pi*c0/wl
    k = REAL(ri2, KIND=dp)*omega/c0
    f = focus%focal
    w0 = focus%waist
    na = focus%na
    theta_max = ASIN(na/REAL(ri2,KIND=dp))
    Cprefix = (0,1)*k*f*EXP(-(0,1)*k*f)*(1.0_dp/(2.0_dp*pi))    &
             *(1.0_dp/2.0_dp)*(SQRT(ri1/ri2))

    ALLOCATE( grid%cj(nbasis*2,0:nMax+2,grid%ny,grid%nx,grid%nz), &
              grid%cjs0c0(nbasis*2,0:nMax+2,grid%ny,grid%nx,grid%nz), &
              grid%cjs0c1(nbasis*2,0:nMax+2,grid%ny,grid%nx,grid%nz), &
              grid%cjs1c0(nbasis*2,0:nMax+2,grid%ny,grid%nx,grid%nz) )

    ALLOCATE( intJjS0C0(-(nMax+2):nMax+2), &
              intJjS0C1(-(nMax+2):nMax+2), &
              intJjS1C0(-(nMax+2):nMax+1) )
    ALLOCATE( intJjS0C0_done(-(nMax+2):nMax+2), &
              intJjS0C1_done(-(nMax+2):nMax+2), &
              intJjS1C0_done(-(nMax+2):nMax+1) )
    ALLOCATE( intJjS00C01(-(nMax+2):nMax+2), &
              intJjS00C01_done(-(nMax+2):nMax+2) )
    ALLOCATE( En(-nMax:nMax,1:3) )

    ! get filling factor, only works for circ aperture type.
    IF ( pupil0%aperture%type=='circ' ) THEN
      f0 = REAL(pupil0%aperture%circ%r, KIND=dp)
    ELSE
      ! filling factor is set to 1 for all other apertures.
      !!!! to be extended.
      f0 = 1.0_dp
    END IF

    phase_type = TRIM(pupil0%phase%type)
    thetaRot = pupil0%phase%thetaRot
    vJones = pupil0%phase%vJones
    ! loop each point in the 3D grid
    ibar = 0 ! counter for the progress bar
    WRITE(*,*) 'Focal field calculation is running ...'
    DO iz = 1, grid%nz
      DO ix = 1, grid%nx
        DO iy = 1, grid%ny

          x = grid%x(ix)
          y = grid%y(iy)
          z = grid%z(iz)

          rho = SQRT(x**2+y**2)
          varphi = ATAN2(y,x)

          intJjS0C0 = CMPLX(0.0_dp,0.0_dp)
          intJjS0C1 = CMPLX(0.0_dp,0.0_dp)
          intJjS1C0 = CMPLX(0.0_dp,0.0_dp)
          intJjS00C01 = CMPLX(0.0_dp,0.0_dp)
          En = CMPLX(0.0_dp, 0.0_dp)
          intJjS0C0_done = .FALSE.
          intJjS0C1_done = .FALSE.
          intJjS1C0_done = .FALSE.
          intJjS00C01_done = .FALSE.

          DO n = nMax, -nMax, -1
            IF ( pupil0%cn(n) == CMPLX(0.0_dp, 0.0_dp) ) THEN
              CYCLE
            END IF

            DO j = n+2, n-2, -2
              IF ( intJjS0C0_done(j) .EQV. .FALSE. ) THEN
                IF ( intJjS0C0_done(-j) .EQV. .TRUE. ) THEN
                  intJjS0C0(j)=(-1)**(-j)*intJjS0C0(-j)
                ELSE
                  s=0; c=0;
                  ! intJjS0C0(j)=IntegralJjSsCc(j, s, c)
                  CALL IntegralJjSsCc(j, s, c, intJjS0C0(j), grid%cjs0c0(:,j,iy,ix,iz))
                END IF
              END IF

              IF ( intJjS0C1_done(j) .EQV. .FALSE. ) THEN
                IF ( intJjS0C1_done(-j) .EQV. .TRUE. ) THEN
                  intJjS0C1(j)=(-1)**(-j)*intJjS0C1(-j)
                ELSE
                  s=0; c=1;
                  ! intJjS0C1(j)=IntegralJjSsCc(j, s, c)
                  CALL IntegralJjSsCc(j, s, c, intJjS0C1(j), grid%cjs0c1(:,j,iy,ix,iz))
                END IF
              END IF

              !IF ( intJjS00C01_done(j) == .FALSE. ) THEN
              !  IF ( intJjS00C01_done(-j) == .TRUE. ) THEN
              !    intJjS00C01(j)=(-1)**(-j)*intJjS00C01(-j)
              !  ELSE
              !    s1=0; s2=0; c1=0; c2=1
              !    intJjS00C01(j)=IntegralJjSssCcc(l, s1, s2, c1, c2)
              !  END IF
              !END IF

              intJjS0C0_done(j) = .TRUE.
              intJjS0C1_done(j) = .TRUE.
            END DO

            DO j = n+1, n-1, -2
              IF ( intJjS1C0_done(j) .EQV. .FALSE. ) THEN
                IF ( intJjS1C0_done(-j) .EQV. .TRUE. ) THEN
                  intJjS1C0(j)=(-1)**(-j)*intJjS1C0(-j)
                ELSE
                  s=1; c=0;
                  !intJjS1C0(j)=IntegralJjSsCc(j, s, c)
                  CALL IntegralJjSsCc(j, s, c, intJjS1C0(j), grid%cjs1c0(:,j,iy,ix,iz))
                END IF
              END IF

              intJjS1C0_done(j) = .TRUE.
            END DO

            ! calc. field from x-polarized hg00 input beam
            varphiRot = varphi - thetaRot(1)*deg2rad

            Etilde(1) = 2*EXP((0,1)*n*varphiRot)*(intJjS0C0(n)+intJjS0C1(n))  &
              +EXP((0,1)*(n+2)*varphiRot)*(intJjS0C0(n+2)-intJjS0C1(n+2))     &
              +EXP((0,1)*(n-2)*varphiRot)*(intJjS0C0(n-2)-intJjS0C1(n-2))
            Etilde(2) = (0,1)*( &
              EXP((0,1)*(n-2)*varphiRot)*(intJjS0C0(n-2)-intJjS0C1(n-2))      &
              -EXP((0,1)*(n+2)*varphiRot)*(intJjS0C0(n+2)-intJjS0C1(n+2)))
            Etilde(3) = 2*(0,1)*( EXP((0,1)*(n-1)*varphiRot)*intJjS1C0(n-1)   &
                                 -EXP((0,1)*(n+1)*varphiRot)*intJjS1C0(n+1) )
            ! rotate Etilde(1:2) in thetaRot(1) [deg]
            Etilde = MATMUL( matrix_rz(thetaRot(1)*deg2rad), Etilde )
            En(n,:) = En(n,:) + pupil0%cn(n)*(0,1)**n*Etilde(:)*vJones(1)

            ! calc. field from y-polarized hg00 input beam
            varphiRot = varphi - pi/2 - thetaRot(2)*deg2rad

            Etilde(1) = 2*EXP((0,1)*n*varphiRot)*(intJjS0C0(n)+intJjS0C1(n))  &
              +EXP((0,1)*(n+2)*varphiRot)*(intJjS0C0(n+2)-intJjS0C1(n+2))     &
              +EXP((0,1)*(n-2)*varphiRot)*(intJjS0C0(n-2)-intJjS0C1(n-2))
            Etilde(2) = (0,1)*( &
              EXP((0,1)*(n-2)*varphiRot)*(intJjS0C0(n-2)-intJjS0C1(n-2))      &
              -EXP((0,1)*(n+2)*varphiRot)*(intJjS0C0(n+2)-intJjS0C1(n+2)))
            Etilde(3) = 2*(0,1)*( EXP((0,1)*(n-1)*varphiRot)*intJjS1C0(n-1)   &
                                 -EXP((0,1)*(n+1)*varphiRot)*intJjS1C0(n+1) )

            ! rotate pi/2 in order to convert hg00 from x-pol to y-pol
            Etilde = MATMUL( matrix_rz(pi/2), Etilde )
            ! rotate Etilde(1:2) in thetaRot(2) [deg]
            Etilde = MATMUL( matrix_rz(thetaRot(2)*deg2rad), Etilde )

            En(n,:) = En(n,:) + pupil0%cn(n)*(0,1)**n*Etilde(:)*vJones(2)
          END DO

          field%ex(iy,ix,iz) = Cprefix*E0*pi*SUM(En(:,1), DIM=1)
          field%ey(iy,ix,iz) = Cprefix*E0*pi*SUM(En(:,2), DIM=1)
          field%ez(iy,ix,iz) = Cprefix*E0*pi*SUM(En(:,3), DIM=1)

          ! show progress
          CALL show_progress( REAL(iy*ix*iz, KIND=dp)&
                              /REAL(grid%ny*grid%nx*grid%nz, KIND=dp), ibar)
        END DO! iy=
      END DO! ix=
    END DO! iz=
    WRITE(*,*) 'Focal field calculation is finished.'

    DEALLOCATE( intJjS0C0, intJjS0C1, intJjS1C0, intJjS00C01, &
      intJjS0C0_done, intJjS0C1_done, intJjS1C0_done, intJjS00C01_done )
    DEALLOCATE( En )

  CONTAINS

    !> The integral is integrated over diverging angle \f$\theta\f$.
    !> \f[
    !! I_{JjSmCn} = \int_0^{\theta_{max}} g(\theta)
    !!          J_j(k\rho sin(\theta)) sin^s(\theta) cos^c(\theta) d\theta,
    !! \f]
    !> with \f$g\left(\theta\right) \equiv
    !! \sqrt{\cos\theta}\sin\theta e^{ikz\cos\theta}f_{w}\left(\theta\right)\f$,
    !! where \f$f_{w}\left(\theta\right) \equiv
    !! e^{-\frac{1}{f_{0}^{2}}\frac{\sin^{2}\theta}{\sin^{2}\theta_{\max}}}\f$
    !! is a common phase factor that is shared by the Gaussian laser modes.
    !> The meaning of the subscript,
    !> J: Bessel function of the first kind
    !> j: Integer order of the first kind Bessel function
    !> S: Sine function
    !> s: Power of sine function excluding the common part
    !> C: Cosine function
    !> c: Power of cosine function excluding the common part
    SUBROUTINE IntegralJjSsCc(j, s, c, integ, cj)
    ! FUNCTION IntegralJjSsCc(j,s,c) RESULT(integ)
      INTEGER, INTENT(IN)         :: j,s,c
      COMPLEX(KIND=dp), INTENT(OUT) :: integ
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT) :: cj

      REAL(KIND=dp)               :: theta_bessel, delta_bessel, &
                                     theta0, theta1, theta2
      REAL(KIND=dp)               :: maxerr, krho_tol=1D1, kz_tol=1D0
      ! REAL(KIND=dp)               :: maxerr, krho_tol=1D3, kz_tol=1D3
      INTEGER                     :: maxDepth, ierr, n

      ! accuracy control for the Adaptive Simpson's method
      maxerr = 1.0D-9
      maxDepth = 30
      ! n = nbasis ! n basis in bessel_transform
      integ = CMPLX(0.0_dp, 0.0_dp)
      ! WRITE(*,*) j, s, c
      theta0 = 0.0_dp

      IF ( TRIM(phase_type) == 'bessel' ) THEN
        theta_bessel = pupil0%phase%bessel%theta*deg2rad
        delta_bessel = pupil0%phase%bessel%delta*deg2rad
        IF ( delta_bessel == 0.0_dp ) THEN
          ! infinitesimal thin ring, no integration
          integ = integrandTheta(theta_bessel)
        ELSE
          theta1 = MAX(0.0_dp, theta_bessel-delta_bessel/2.0_dp)
          theta2 = theta_bessel+delta_bessel/2.0_dp

          integ = asqz(integrandTheta, theta1, theta2, maxerr, maxDepth)
        END IF
      ELSE
        theta1 = 0.0_dp
        theta2 = theta_max

        integ = asqz(integrandTheta, theta1, theta2, maxerr, maxDepth)
      END IF

      ! WRITE(*,*) x,y,z
!      WRITE(*,*) '-------------------------------------------------------------'
!      WRITE(*,*) cyl_bessel_zeros(MAX(0,j-1),1)
!      WRITE(*,*) '-------------------------------------------------------------'
!      WRITE(*,*) j, k*rho, k*z, theta0*rad2deg, theta1*rad2deg, theta2*rad2deg
      ! The 1st method for numerical integration
      ! CALL zhankel(gx_ht, SIN(theta1), SIN(theta2), j, k*rho, integ, ierr)
      ! The subroutine zhankel is not working well, as it requires writing
      ! the integrand function in Fortran 77 in order to use COMMON blocks.
      ! Could check this method later.

      ! The 2nd method for numerical integration.
      ! Split into 1 or 2 subintervals for Bessel (or Hankel) transfrom.
      !IF ( ( k*rho*SIN(theta1) < cyl_bessel_zeros(MAX(0,j-1),1) ) &
      !    .AND. ( cyl_bessel_zeros(MAX(0,j-1),1) < k*rho*SIN(theta2) ) ) THEN
      !  theta0 = theta1
      !  theta1 = ASIN(cyl_bessel_zeros(MAX(0,j-1),1)/(k*rho))
      !  integ = besselj_transfrom(gx,SIN(theta0),SIN(theta1),j+1,k*z,k*rho,n)
      !END IF
      !integ = integ + &
      !        besselj_transfrom(gx,SIN(theta1),SIN(theta2),j+1,k*z,k*rho,n)

      ! When k*rho and k*z is small, use Simpson's quadrature rule.
!      IF ( ( k*rho < krho_tol) .AND. ( ABS(k*z) < kz_tol ) ) THEN
!        integ = asqz(integrandTheta, theta1, theta2, maxerr, maxDepth)
!      ELSE
!        ! compare lower bound with the 1st maxima of bessel function
!        ! IF ( ( k*rho*SIN(theta1) < cyl_bessel_zeros(MAX(0,j-1),1) ) &
!        !     .AND. ( cyl_bessel_zeros(MAX(0,j-1),1) < k*rho*SIN(theta2) ) ) THEN
!        !   theta0 = theta1
!        !   theta1 = ASIN(cyl_bessel_zeros(MAX(0,j-1),1)/(k*rho))
!        !   integ = besselj_transfrom(gx,SIN(theta0),SIN(theta1),j+1,k*z,k*rho,n)
!        ! END IF
!        ! integ = integ + &
!        !         besselj_transfrom(gx,SIN(theta1),SIN(theta2),j+1,k*z,k*rho,n)
!        ! integ = besselj_transfrom(gx,SIN(theta1),SIN(theta2),j+1,k*z,k*rho,n)
!        ! CALL besselj_transfrom_detail( gx, SIN(theta1), SIN(theta2), j+1, &
!        !                                k*z, k*rho, n, polyName, integ, cj )
!        ! CALL besselj_transfrom_detail( gx, SIN(theta1), SIN(theta2), j+1, &! debug
!        !                                0.0_dp, 4.267823982235190_dp, n, polyName, integ, cj )! debug
!      END IF
    !END FUNCTION IntegralJjSsCc
    END SUBROUTINE IntegralJjSsCc

    FUNCTION integrandTheta(theta) RESULT(res)
      REAL(KIND=dp), INTENT(IN)   :: theta
      COMPLEX(KIND=dp)            :: res

      res = besselj( j, k*rho*SIN(theta) ) * &
            g(theta)*(SIN(theta)**s)*(COS(theta)**c)
    END FUNCTION integrandTheta

    FUNCTION gx(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)   :: x ! x = sin(theta)
      COMPLEX(KIND=dp)            :: res

      ! Host Association: k, z, phase_type, l, w0
      ! basis: HG00 mode
      ! res = fx_w(x) * (1-x**2)**((2.0_dp*c-1.0_dp)/4.0_dp) * x**(s+1) &
      !       * EXP( (0,1)*k*z*(1-x**2)**0.5_dp )
      res = fx_w(x) * (1-x**2)**((2.0_dp*REAL(c, KIND=dp)-1.0_dp)/4.0_dp) * x**(s+1)
      ! EXP( (0,1)*k*z*(1-x**2)**0.5_dp ) is considered in bessel_transform
      ! as the factor in term exp(i*r1*(1-x^2)^0.5)*J_nu(r_2 * x)
      ! Check bessel_transform for more details.

      ! multiply extra factor for each type of phase
      IF ( TRIM(phase_type) == 'petal' .OR. &
           TRIM(phase_type) == 'petal_rect' ) THEN
        res = res * ((SQRT(2.0_dp)*f0*x/w0)**ABS(j))
      ELSE IF ( TRIM(phase_type) == 'bessel' ) THEN
        res = res * x
      ELSE IF ( TRIM(phase_type) == 'azi' ) THEN
        res = res * x
      ELSE IF ( TRIM(phase_type) == 'rad' ) THEN
        res = res * x
      ELSE IF ( TRIM(phase_type) == 'hg01' ) THEN
        res = res * x
      END IF
    END FUNCTION gx

    ! Definition for the integrand FUN(X) in the following hankel transform
    !
    !                 ( infty
    ! HT(FUN,NU,XI) = |      (XI * X)**(1/2) JNUX(XI * X) FUN(X) dX
    !                 ) 0
    FUNCTION gx_ht(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)   :: x! x = sin(theta)
      COMPLEX(KIND=dp)            :: res

      ! Host Association: k, z, phase_type, l, w0
      ! basis: HG00 mode
      WRITE(*,*) rho
      res = fx_w(x) * (1-x**2)**((2*c+1)/4) &
            * x**(s+1) * EXP( (0,1)*k*z*SQRT(1-x**2) ) &
            / (k*rho*x)**(0.5)
      ! multiply extra factor for each type of phase
      IF ( TRIM(phase_type) == 'petal' .OR. &
           TRIM(phase_type) == 'petal_rect' ) THEN
        res = res * ((SQRT(2.0_dp)*f0*x/w0)**ABS(j))
      ELSE IF ( TRIM(phase_type) == 'bessel' ) THEN
        res = res * x
      ELSE IF ( TRIM(phase_type) == 'hg01' ) THEN
        res = res * x
      END IF
    END FUNCTION gx_ht

    ! inlcude the filling factor
    FUNCTION fx_w(x) RESULT(res)
      REAL (KIND=dp), INTENT(IN)  :: x
      REAL (KIND=dp)              :: res

      ! Host Association: f0, theta_max
      res = EXP( -( x**2 / (f0**2 * SIN(theta_max)**2) ) )
    END FUNCTION fx_w

    FUNCTION g(theta) RESULT(res)
      REAL (KIND=dp), INTENT(IN)  :: theta
      COMPLEX (KIND=dp)           :: res

      ! Host Association: k, z, phase_type, l, w0
      ! basis: HG00 mode
      res = f_w(theta) * SQRT( COS(theta) ) &
            * SIN(theta) * EXP( (0,1)*k*z*COS(theta) )
      ! multiply extra factor for each type of phase
      IF ( TRIM(phase_type) == 'petal' .OR. &
           TRIM(phase_type) == 'petal_rect' ) THEN
        res = res * ((SQRT(2.0_dp)*f0*SIN(theta)/w0)**ABS(j))
      ELSE IF ( TRIM(phase_type) == 'bessel' ) THEN
        res = res * SIN(theta)
      ELSE IF ( TRIM(phase_type) == 'azi' ) THEN
        res = res * SIN(theta)
      ELSE IF ( TRIM(phase_type) == 'rad' ) THEN
        res = res *  SIN(theta)
      ELSE IF ( TRIM(phase_type) == 'hg01' ) THEN
        res = res * SIN(theta)
      END IF
    END FUNCTION g

    ! inlcude the filling factor
    FUNCTION f_w(theta) RESULT(res)
      REAL (KIND=dp), INTENT(IN)  :: theta
      REAL (KIND=dp)              :: res

      ! Host Association: f0, theta_max
      res = EXP( -( SIN(theta)**2 / (f0**2 * SIN(theta_max)**2) ) )
    END FUNCTION f_w

  END SUBROUTINE focal_field

END MODULE focal
