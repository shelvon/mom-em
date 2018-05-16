! MODULE: focal
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines to focal field definition and calculation.
MODULE focal
  USE constants
  USE pupil
  USE source
  USE quad
  USE aux

  IMPLICIT NONE

  TYPE points
    REAL (KIND=dp), DIMENSION(3)        :: xyz
  END TYPE points

  TYPE fields
    COMPLEX (KIND=dp), DIMENSION(3)     :: data ! data(1:3) for Ex, Ey, Ez components
  END TYPE fields

  TYPE data_focal
    INTEGER                     :: nx, ny, nz
    INTEGER                     :: jMax
    REAL (KIND=dp)              :: xa, xb, ya, yb, za, zb
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)         :: x, y, z
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridx
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridy
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridz
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:,:)   :: grid
    ! TYPE(points), ALLOCATABLE, DIMENSION(:,:,:)     :: grid
    TYPE(fields), ALLOCATABLE, DIMENSION(:,:,:)     :: e
    TYPE(fields), ALLOCATABLE, DIMENSION(:,:,:)     :: h
  END TYPE data_focal

CONTAINS
  SUBROUTINE meshgrid(focal)
    TYPE(data_focal), INTENT(INOUT) :: focal

    INTEGER :: ix, iy, iz, ixyz
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:)         :: x, y, z
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridx
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridy
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)     :: gridz

    ALLOCATE(x(1:focal%nx), y(1:focal%ny), z(1:focal%nz))
    ALLOCATE(gridx(1:focal%ny, 1:focal%nx, 1:focal%nz))
    ALLOCATE(gridy(1:focal%ny, 1:focal%nx, 1:focal%nz))
    ALLOCATE(gridz(1:focal%ny, 1:focal%nx, 1:focal%nz))

    ixyz = 0
    ! dimension 3
    IF (focal%nz>1) THEN
      z=(/(focal%za+(iz-1)*(focal%zb-focal%za)/(focal%nz-1), iz=1, focal%nz, 1)/)
    ELSE IF(focal%nz==1) THEN
      z=focal%za
    END IF
    ! dimension 2
    IF (focal%nx>1) THEN
      x=(/(focal%xa+(ix-1)*(focal%xb-focal%xa)/(focal%nx-1), ix=1, focal%nx, 1)/)
    ELSE IF(focal%nx==1) THEN
      x=focal%xa
    END IF
    ! dimension 1
    IF (focal%ny>1) THEN
      y=(/(focal%ya+(iy-1)*(focal%yb-focal%ya)/(focal%ny-1), iy=1, focal%ny, 1)/)
    ELSE IF(focal%ny==1) THEN
      y=focal%ya
    END IF

    gridx=(SPREAD(SPREAD(x,1,focal%ny),3,focal%nz))
    gridy=(SPREAD(SPREAD(y,2,focal%nx),3,focal%nz))
    gridz=(SPREAD(SPREAD(z,1,focal%nx),1,focal%ny))
    focal%grid(1,:,:,:)=gridx
    focal%grid(2,:,:,:)=gridy
    focal%grid(3,:,:,:)=gridz
    focal%gridx=gridx
    focal%gridy=gridy
    focal%gridz=gridz
    focal%x=x
    focal%y=y
    focal%z=z

    DEALLOCATE(x,y,z,gridx,gridy,gridz)

  END SUBROUTINE meshgrid

  SUBROUTINE field_focal(wl,ri1,ri2,src,pupil,focal,name)
    REAL (KIND=dp), INTENT(IN)          :: wl
    COMPLEX (KIND=dp), INTENT(IN)       :: ri1, ri2
    TYPE(srcdata), INTENT(IN)           :: src
    TYPE(data_pupil), INTENT(INOUT)     :: pupil
    TYPE(data_focal), INTENT(INOUT)     :: focal
    CHARACTER (LEN=*), INTENT(IN)       :: name

    REAL (KIND=dp)      :: E0
    REAL (KIND=dp)      :: omega, k, f, w0, na, theta_max, f0
    COMPLEX (KIND=dp)   :: Cprefix, prefixHG00, prefixLG0l
    REAL (KIND=dp)      :: x,y,z, rho, varphi
    INTEGER             :: ix, iy, iz, jExtra
    COMPLEX (KIND=dp), DIMENSION(0:focal%jMax)      :: intJjS0C0, intJjS0C1, intJjS1C0, intJjS1C1, intJjS2C0
    REAL (KIND=dp), DIMENSION(1:2,0:focal%jMax+2)   :: dj
    COMPLEX (KIND=dp), DIMENSION(-(focal%jMax+2):focal%jMax+2)    :: cj
    COMPLEX (KIND=dp), DIMENSION(-focal%jMax:focal%jMax,1:3)      :: Ej
    CHARACTER (LEN=256) :: filename
    CHARACTER (LEN=32)  :: phase_type

    INTEGER :: charge,l,j,s,c, focustype

    E0 = 1.0_dp

    ! get focustype
    IF(src%type==src_focus_x) THEN
      focustype = focustype_x
    END IF

    ! assign extra orders in Fourier series expansion for each beam
    IF(focustype==focustype_x) THEN
      jExtra = 2
    END IF

    ALLOCATE(focal%e(1:focal%ny, 1:focal%nx, 1:focal%nz))
    ! ALLOCATE(focal%h(1:focal%ny, 1:focal%nx, 1:focal%nz))
    ALLOCATE(pupil%phase%rect%dj(1:2,0:focal%jMax+jExtra))

    omega = 2.0_dp*pi*c0/wl
    k = REAL(ri2,KIND=dp)*omega/c0

    f = src%focal
    w0 = src%waist
    na = src%napr
    theta_max = ASIN(na/REAL(ri2,KIND=dp))

    Cprefix = (0,1)*k*f*EXP(-(0,1)*k*f)*(1.0_dp/(2.0_dp*pi))*(1.0_dp/2.0_dp)*(SQRT(ri1/ri2))

    phase_type = pupil%phase%type

    ! get filling factor, only works for circ aperture type.
    IF(pupil%aperture%type=='circ') THEN
      f0 = pupil%aperture%circ%r
    ELSE
      f0 = 1! filling factor is set to 1 for all other apertures.
    END IF

    ! get the coefficients of Fourier series expansion for the exponent of the phase factor
    !DO j=0,focal%jMax+jExtra
    !  dj(:,j) = coeffRectPhase(pupil,j)
    !  ! WRITE(*,*) dj(1,j), dj(2,j)
    !END DO
    !pupil%phase%rect%dj=dj

    ! get the complex coefficients of complex Fourier series expansion of the
    ! entire pupil function, but only for the one doesn't depend on radial variable.
    IF(pupil%aperture%type=='circ' .AND. &
    (phase_type=='rect' &
    .OR. phase_type=='vortex' &
    .OR. phase_type=='petal'  &
    .OR. phase_type=='bessel' &
    .OR. phase_type=='petal_rect')) THEN
      IF(phase_type=='vortex') THEN
        cj = CMPLX(0.0_dp,0.0_dp)
        charge = pupil%phase%vortex%charge
        cj(charge) = coeffPupil(pupil,charge)
      ELSE IF(phase_type=='petal') THEN
        cj = CMPLX(0.0_dp,0.0_dp)
        l = pupil%phase%petal%l
        cj(l) = coeffPupil(pupil,l)
        cj(-l) = coeffPupil(pupil,-l)
      ELSE IF(phase_type=='bessel') THEN
        cj = CMPLX(1.0_dp,0.0_dp)
      ELSE
        DO j=-(focal%jMax+jExtra),focal%jMax+jExtra
          cj(j) = coeffPupil(pupil,j)
        END DO
      END IF
    ELSE
      WRITE(*,*) 'series expansion approach is not supported for this phase mask!'
      RETURN
    END IF

    ! loop each point in the 3D grid
    DO iz=1,focal%nz
      DO ix=1,focal%nx
        DO iy=1,focal%ny
          x=focal%grid(1,iy,ix,iz)
          y=focal%grid(2,iy,ix,iz)
          z=focal%grid(3,iy,ix,iz)

          rho=SQRT(x**2+y**2)
          varphi=ATAN2(y,x)

          intJjS0C0 = CMPLX(0.0_dp,0.0_dp)
          intJjS0C1 = CMPLX(0.0_dp,0.0_dp)
          intJjS1C0 = CMPLX(0.0_dp,0.0_dp)
          intJjS1C1 = CMPLX(0.0_dp,0.0_dp)
          intJjS2C0 = CMPLX(0.0_dp,0.0_dp)
          Ej = CMPLX(0.0_dp,0.0_dp)

          ! parallel computation of series of integrals
          IF(focustype==focustype_x) THEN

            IF(phase_type=='rect' .OR. phase_type=='petal_rect') THEN
              DO j=-focal%jMax,focal%jMax
                IF(cj(j)==0 .AND. cj(j+2)==0 .AND. cj(j-2)==0 .AND. cj(j+1)==0 .AND. cj(j-1)==0) THEN
                  CONTINUE
                END IF

                s=0;c=0;    intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
                s=0;c=1;    intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
                s=1;c=0;    intJjS1C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)

                ! field of j-order as j-index with respect to I_{JjSsCc}
                Ej(j,1) = (2*cj(j)-cj(j+2)-cj(j-2))*intJjS0C0(j)+&
                          (2*cj(j)+cj(j+2)+cj(j-2))*intJjS0C1(j)
                Ej(j,2) = -(0,1)*(cj(j+2)-cj(j-2))*(intJjS0C0(j)-intJjS0C1(j))
                Ej(j,3) = -2*(cj(j+1)+cj(j-1))*intJjS1C0(j)

                prefixHG00 = E0*pi*((0,1)**j)*EXP((0,1)*j*varphi)
                Ej(j,:) = prefixHG00*Ej(j,:)
              END DO! j=
            ! vortex beam
            ELSE IF(phase_type=='vortex') THEN
              charge = pupil%phase%vortex%charge

              j=charge;s=0;c=0;     intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge+2;           intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge-2;           intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge;s=0;c=1;     intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge+2;           intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge-2;           intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge+1;s=1;c=0;   intJjS1C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=charge-1;           intJjS1C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)

              ! calc. field
              Ej(charge,1) = 2*(intJjS0C0(charge)+intJjS0C1(charge))+&
                EXP((0,1)*(2)*varphi)*(intJjS0C0(charge+2)-intJjS0C1(charge+2))+&
                EXP((0,1)*(-2)*varphi)*(intJjS0C0(charge-2)-intJjS0C1(charge-2))

              Ej(charge,2)=(0,1)*&
                (EXP((0,1)*(-2)*varphi)*(intJjS0C0(charge-2)-intJjS0C1(charge-2))-&
                EXP((0,1)*2*varphi)*(intJjS0C0(charge+2)-intJjS0C1(charge+2)))

              Ej(charge,3) = 2*(0,1)*&
                (EXP((0,1)*(-1)*varphi)*intJjS1C0(charge-1)-&
                EXP((0,1)*1*varphi)*intJjS1C0(charge+1))

              prefixHG00 = cj(charge)*E0*pi*((0,1)**charge)*EXP((0,1)*charge*varphi)
              Ej(charge,:) = prefixHG00*Ej(charge,:)

            ! bessel beam
            ELSE IF(phase_type=='bessel') THEN
              ! temporarily consider only the o-order of Bessel beam, i.e. J_0
              charge=0
              theta_max = pupil%phase%bessel%theta

              j=charge;s=2;c=0;     intJjS2C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)

              ! calc. field
              Ej(charge,3) = intJjS2C0(charge)

              prefixHG00 = cj(charge)*E0*pi*((0,1)**charge)*EXP((0,1)*charge*varphi)
              Ej(charge,:) = prefixHG00*Ej(charge,:)

            ! petal beam
            ELSE IF(phase_type=='petal') THEN
              l = pupil%phase%petal%l

              j=l;s=0;c=0;     intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l+2;           intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l-2;           intJjS0C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l;s=0;c=1;     intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l+2;           intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l-2;           intJjS0C1(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l+1;s=1;c=0;   intJjS1C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)
              j=l-1;           intJjS1C0(j)=IntegralJjSsCc(j,s,c,theta_max,k,rho,z)

              ! calc. field, l
              Ej(l,1) = 2*(intJjS0C0(l)+intJjS0C1(l))+&
                EXP((0,1)*(2)*varphi)*(intJjS0C0(l+2)-intJjS0C1(l+2))+&
                EXP((0,1)*(-2)*varphi)*(intJjS0C0(l-2)-intJjS0C1(l-2))

              Ej(l,2)=(0,1)*&
                (EXP((0,1)*(-2)*varphi)*(intJjS0C0(l-2)-intJjS0C1(l-2))-&
                EXP((0,1)*(2)*varphi)*(intJjS0C0(l+2)-intJjS0C1(l+2)))

              Ej(l,3) = 2*(0,1)*(EXP((0,1)*(-1)*varphi)*intJjS1C0(l-1)-&
                EXP((0,1)*1*varphi)*intJjS1C0(l+1))

              prefixLG0l = cj(l)*E0*SQRT(2/pi/factorial_n(ABS(l)))*(1/w0)*pi*((0,1)**l)*EXP((0,1)*l*varphi)
              Ej(l,:) = prefixLG0l*Ej(l,:)

              ! calc. field, -l
              Ej(-l,1) = 2*((-1)**l)*(intJjS0C0(l)+intJjS0C1(l))+&
                (EXP((0,1)*(2)*varphi)+((-1)**(l+2))*EXP((0,1)*(-2)*varphi))*(intJjS0C0(l+2)-intJjS0C1(l+2))+&
                (EXP((0,1)*(-2)*varphi)+((-1)**(l-2))*EXP((0,1)*(2)*varphi))*(intJjS0C0(l-2)-intJjS0C1(l-2))

              Ej(-l,2)=(0,1)*&
                (EXP((0,1)*(-2)*varphi)*((-1)**(l+2))*(intJjS0C0(l+2)-intJjS0C1(l+2))-&
                EXP((0,1)*(2)*varphi)*((-1)**(l-2))*(intJjS0C0(l-2)-intJjS0C1(l-2)))

              Ej(-l,3) = 2*(0,1)*(EXP((0,1)*(-1)*varphi)*((-1)**(l+1))*intJjS1C0(l+1)-&
                EXP((0,1)*1*varphi)*((-1)**(l-1))*intJjS1C0(l-1))

              prefixLG0l = cj(-l)*E0*SQRT(2/pi/factorial_n(ABS(l)))*(1/w0)*pi*((0,1)**(-l))*EXP((0,1)*(-l)*varphi)
              Ej(-l,:) = prefixLG0l*Ej(-l,:)
            END IF
            focal%e(iy,ix,iz)%data(1) = Cprefix*SUM(Ej(:,1))
            focal%e(iy,ix,iz)%data(2) = Cprefix*SUM(Ej(:,2))
            focal%e(iy,ix,iz)%data(3) = Cprefix*SUM(Ej(:,3))

          END IF!(focustype==focustype_x) THEN

        END DO! iy=
        WRITE(*,'(A,I0,A)') 'Computed ', NINT(100*REAL(ix)/REAL(focal%nx)), ' percent of sources'
      END DO! ix=

      ! write to file
      WRITE(filename, '(A,I0)') '-iz', iz
      filename = name // filename
      CALL write_data(TRIM(filename) // '-ex-re.dat', REAL(focal%e(:,:,iz)%data(1)))
      CALL write_data(TRIM(filename) // '-ey-re.dat', REAL(focal%e(:,:,iz)%data(2)))
      CALL write_data(TRIM(filename) // '-ez-re.dat', REAL(focal%e(:,:,iz)%data(3)))

      CALL write_data(TRIM(filename) // '-ex-im.dat', AIMAG(focal%e(:,:,iz)%data(1)))
      CALL write_data(TRIM(filename) // '-ey-im.dat', AIMAG(focal%e(:,:,iz)%data(2)))
      CALL write_data(TRIM(filename) // '-ez-im.dat', AIMAG(focal%e(:,:,iz)%data(3)))

!      CALL write_data(TRIM(filename) // '-ex-re.dat', REAL(RESHAPE(focal%e(:,:,iz)%data(1),(/focal%ny,focal%nx/)), KIND=dp))
!      CALL write_data(TRIM(filename) // '-ey-re.dat', REAL(RESHAPE(focal%e(:,:,iz)%data(2),(/focal%ny,focal%nx/)), KIND=dp))
!      CALL write_data(TRIM(filename) // '-ez-re.dat', REAL(RESHAPE(focal%e(:,:,iz)%data(3),(/focal%ny,focal%nx/)), KIND=dp))
!
!      CALL write_data(TRIM(filename) // '-ex-im.dat', AIMAG(RESHAPE(focal%e(:,:,iz)%data(1),(/focal%ny,focal%nx/)), KIND=dp))
!      CALL write_data(TRIM(filename) // '-ey-im.dat', AIMAG(RESHAPE(focal%e(:,:,iz)%data(2),(/focal%ny,focal%nx/)), KIND=dp))
!      CALL write_data(TRIM(filename) // '-ez-im.dat', AIMAG(RESHAPE(focal%e(:,:,iz)%data(3),(/focal%ny,focal%nx/)), KIND=dp))
    END DO! iz=

!      ! write to file
!      iy = 1
!      filename = name // filename
!      CALL write_data(TRIM(filename) // '-ex-re.dat', REAL(focal%e(iy,:,:)%data(1)))
!      CALL write_data(TRIM(filename) // '-ey-re.dat', REAL(focal%e(iy,:,:)%data(2)))
!      CALL write_data(TRIM(filename) // '-ez-re.dat', REAL(focal%e(iy,:,:)%data(3)))
!
!      CALL write_data(TRIM(filename) // '-ex-im.dat', AIMAG(focal%e(iy,:,:)%data(1)))
!      CALL write_data(TRIM(filename) // '-ey-im.dat', AIMAG(focal%e(iy,:,:)%data(2)))
!      CALL write_data(TRIM(filename) // '-ez-im.dat', AIMAG(focal%e(iy,:,:)%data(3)))

    WRITE(*,*) 'focal field calculation is finished.'

  CONTAINS! contains functions inside "SUBROUTINE field_focal"
    ! The general integral to be integed over beaming diverging angle
    ! theta, from 0 to theta_max;
    ! I_JjSmCn = \int_0^{\theta_max} * Icommon *
    !          * J_j(k\rho sin(\theta)) sin^m(\theta) cos^n(\theta) d\theta,
    ! with Icommon = f_w(\theta) exp(ikzcos(\theta)) \sqrt(cos(\theta)) sin(\theta)
    ! The meaning of the subscript,
    ! J: Bessel function of the first kind
    ! j: Integer order of the first kind Bessel function
    ! S: Sine function
    ! m: Power of sine function excluding the common part
    ! C: Cosine function
    ! n: Power of cosine function excluding the common part
    FUNCTION IntegralJjSsCc(j,s,c,theta_max,k,rho,z) RESULT(integ)
    ! FUNCTION IntegralJjSsCc(j,s,c) RESULT(integ)
    INTEGER, INTENT(IN)         :: j,s,c
    REAL (KIND=dp), INTENT(IN)  :: theta_max,k,rho,z
    REAL (KIND=dp)              :: integ

    REAL (KIND=dp)  :: maxerr
    INTEGER         :: maxDepth

    ! accuracy for the Adaptive Simpson's method
    maxerr = 1D-4
    maxDepth = 20

    IF(phase_type=='bessel') THEN
      integ = integrandTheta(theta_max) ! Here borrow the variable name theta_max->theta
    ELSE
      integ = asqz(integrandTheta, 0.0_dp, theta_max, maxerr, maxDepth)
    END IF
    END FUNCTION IntegralJjSsCc

    FUNCTION integrandTheta(theta) RESULT(res)
    REAL (KIND=dp), INTENT(IN)  :: theta
    COMPLEX (KIND=dp)           :: res
    REAL (KIND=dp)              :: x, cyr, cyi, zr, zi, fnu
    INTEGER                     :: kode, n, nz, ierr
    ! INTEGER :: besm

    ! Host Association: k, rho, s, c
    x = k*rho*SIN(theta)
    zr = x
    zi = 0.0_dp
    fnu = ABS(REAL(j,dp))! take the non-negative order
    kode = 1
    n = 1
    call ZBESJ(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
    ! cyr = BESSEL_JN(j,x) ! requires Fortran 2008 features
    ! CALL besselj(besm,j,maxerr,k*rho*SIN(theta),cyr)! subroutine in bessel.f90
    ! IF(ierr==0) THEN
      ! WRITE(*,*) 'ZBESJ, COMPUTATION COMPLETED'
    ! ELSE
    IF(ierr/=0) THEN
      WRITE(*,*) 'ZBESJ, COMPUTATION FAILED, failed code: ierr = ', ierr, '(focal.f90:integrandTheta)'
      STOP
    END IF

    IF(cyi>tol) THEN
      WRITE(*,*) 'ZBESJ, INTEGER ORDER, IMAGINARY PART, NON-ZERO (focal.f90:integrandTheta)'
      STOP
    END IF

    !for negative integer order, use the following relation
    !  J(-FNU,Z) = ((-1)**FNU)*J(FNU,Z)
    IF(j<0) THEN
      cyr = ((-1)**fnu)*cyr
    END IF

    ! cyr=cyr*EXP(ABS(zi)) ! if kode = 2
    res = CMPLX(cyr, cyi)*integrandCommon(theta)*(SIN(theta)**s)*(COS(theta)**c)
    END FUNCTION integrandTheta

    FUNCTION integrandCommon(theta) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: theta
      COMPLEX (KIND=dp) :: res

      ! Host Association: k, z, phase_type, l, w0
      IF(phase_type=='petal' .OR. phase_type=='petal_rect') THEN
        res = ((SQRT(2.0_dp)*f0*SIN(theta)/w0)**ABS(l))*f_w(theta)*SQRT(COS(theta))*EXP((0,1)*k*z*COS(theta))
      ELSE
        res = f_w(theta)*SQRT(COS(theta))*EXP((0,1)*k*z*COS(theta))
      END IF
    END FUNCTION integrandCommon

    ! inlcude the filling factor
    FUNCTION f_w(theta) RESULT(res)
      REAL (KIND=dp), INTENT(IN) :: theta
      REAL (KIND=dp) :: res

      ! Host Association: f0, theta_max
      res = EXP(-1/(f0**2)*(SIN(theta)**2)/(SIN(theta_max)**2))
    END FUNCTION f_w

  END SUBROUTINE field_focal
END MODULE focal
