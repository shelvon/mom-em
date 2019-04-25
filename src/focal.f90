! MODULE: focal
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines for focal field calculation.
MODULE focal
  USE data
  USE quad
  USE aux
  USE linalg
  USE time
  USE bessel
  USE pupil

  IMPLICIT NONE

CONTAINS

  SUBROUTINE focal_field(focus, wl, ri, pt, Einc, Hinc)
    TYPE(focus_type), INTENT(IN)            :: focus
    REAL(KIND=dp), INTENT(IN)               :: wl
    COMPLEX(KIND=dp), INTENT(IN)            :: ri
    REAL(KIND=dp), DIMENSION(3), INTENT(IN) :: pt
    ! The focal fields Einc and Hinc are calculated in cylindrical coordinate.
    ! [E|H]inc(1): radial field component;
    ! [E|H]inc(2): azimuthal field component;
    ! [E|H]inc(3): z-field component.
    COMPLEX(KIND=dp), DIMENSION(3), INTENT(OUT) :: Einc, Hinc

    !--some beam parameters
    TYPE(pupil_type)  :: pupilE, pupilH
    REAL(KIND=dp)     :: E0, omega, k, f, w0, na, theta_max, f0, ldn, zeta
    COMPLEX(KIND=dp)  :: Cf

    ! spatial and related variables
    ! (x,y,z) the location where the focal field is to be calculated.
    ! (rho,varphi,z) the corresponding cylindrical coordinate representation
    REAL(KIND=dp)     :: x, y, z, rho, varphi

    !--initialization
    f0 = 1.0_dp ! filling factor is set to 1 for all other apertures.
    E0 = focus%E0
    w0 = focus%waist
    na = focus%na
    omega = 2.0_dp*pi*c0/wl
    k = REAL(ri, KIND=dp)*omega/c0
    ldn = wl/REAL(ri, KIND=dp)
    theta_max = ASIN(na/REAL(ri,KIND=dp))
    ! f = focus%focal ! !!!focal length is not an independent parameter!!!
    f = w0*f0/(SIN(theta_max))
    Cf = -(0,1)*k*f/(2.0_dp*pi)*SQRT(1.0_dp/ri)
    zeta = REAL(ri, KIND=dp)/eta0

    ! Generate a general and equivalent pupil function for a specific input beam
    pupilE = generate_pupil(focus, theta_max)
    pupilH = pupil_e2h(pupilE)

    x = pt(1)
    y = pt(2)
    z = pt(3)
    rho = SQRT(x**2+y**2)
    varphi = ATAN2(y,x)
    
    IF (focus%paraxial) THEN
      Einc = E0*vecP(rho, varphi, z, pupilE, ldn, f, w0)
      Hinc = zeta*E0*vecP(rho, varphi, z, pupilH, ldn, f, w0)

    ELSE
      Einc = E0*Cf*vecI(rho, varphi, z, pupilE, k, f, w0, f0, theta_max)
      Hinc = zeta*E0*Cf*vecI(rho, varphi, z, pupilH, k, f, w0, f0, theta_max)

    END IF
  END SUBROUTINE focal_field

  ! calculate the vector form of the focal field under paraxial approximation
  FUNCTION vecP(rho, varphi, z, pupil0, ldn, f, w0)
    TYPE(pupil_type), INTENT(IN)      :: pupil0
    REAL(KIND=dp), INTENT(IN)         :: rho, varphi, z, ldn, f, w0
    COMPLEX(KIND=dp), DIMENSION(1:3)  :: vecP

    REAL(KIND=dp)     :: zR, wz, Rz, gaussian, phig, x, theta
    COMPLEX(KIND=dp)  :: Gouy
    INTEGER           :: n, ir, l, nMin, nMax
    COMPLEX(KIND=dp), DIMENSION(1:3)                :: vecR
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: vecPhi_cn

    nMin = LBOUND(pupil0%Phi(1)%cn,1)
    nMax = UBOUND(pupil0%Phi(1)%cn,1)
    ALLOCATE( vecPhi_cn(1:3,nMin:nMax) )
    vecPhi_cn = CMPLX(0,0, KIND=dp)
    vecP = CMPLX(0,0, KIND=dp)

    zR = pi*w0**2/ldn
    wz = w0*SQRT(1.0_dp+(z/zR)**2)
    gaussian = EXP(-(rho/wz)**2)

    theta = ASIN(rho/f)

    !--vecR(theta) and vecPhi_cn, loop for all three field components
    DO ir = 1, 3
      vecPhi_cn(ir,nMin:nMax) = pupil0%Phi(ir)%cn

      SELECT CASE ( TRIM(pupil0%R(ir)%fun) )
        CASE ('const')
          vecR(ir) = pupil0%R(ir)%const

        CASE ('sin')
          vecR(ir) = SIN(pupil0%R(ir)%sin%x*theta)**pupil0%R(ir)%sin%y

        CASE ('cos')
          vecR(ir) = COS(pupil0%R(ir)%cos%x*theta)**pupil0%R(ir)%cos%y

        CASE ('H1')
          ! f*sin(theta) = rho under paraxial approx.
          vecR(ir) = (2.0_dp/w0)*rho

        CASE ('LG')
          n = pupil0%R(ir)%lg%n
          l = pupil0%R(ir)%lg%l

          x = SQRT(2.0_dp)*rho/wz
          phig = ATAN(z/zR)
          Gouy = EXP(-(0,1)*(2.0_dp*n+ABS(l)+1)*phig)
          Rz = z*(1.0_dp+(zR/z)**2)

          vecR(ir) = Gouy*SQRT(2.0_dp*factorial_n(n)/(pi*factorial_n(n+ABS(l)))) &
                    * (1.0_dp/wz) * (x**ABS(l)) * gLaguerrePoly(n, ABS(l), x**2)

        CASE DEFAULT
          WRITE(*,*) 'Not supported Radial function type!'
          STOP
      END SELECT

!      WRITE(*,*) vecPhi_cn(ir, nMin:nMax)*(/(EXP((0,1)*n*varphi),n=nMin, nMax)/)
!      STOP
      vecP(ir) = vecR(ir)*SUM(vecPhi_cn(ir, nMin:nMax)*(/(EXP((0,1)*n*varphi),n=nMin, nMax)/))
    END DO! ir

    vecP = vecP*gaussian

  END FUNCTION vecP

  ! calculate the vector form of the tightly focal field through the vector
  ! diffraction theory developed by B. Richards and E. Wolf in 1959
  FUNCTION vecI(rho, varphi, z, pupil0, k, f, w0, f0, theta_max)
    TYPE(pupil_type), INTENT(IN)      :: pupil0
    REAL(KIND=dp), INTENT(IN)         :: rho, varphi, z, k, f, w0, f0, theta_max
    COMPLEX(KIND=dp), DIMENSION(1:3)  :: vecI

    REAL(KIND=dp)   :: theta1, theta2, maxErr
    INTEGER         :: maxDepth, nMin, nMax

    ! accuracy control for the Adaptive Simpson's method
    maxErr = 1.0D-6
    maxDepth = 50

    nMin = LBOUND(pupil0%Phi(1)%cn,1)
    nMax = UBOUND(pupil0%Phi(1)%cn,1)

    ! get three field components
    theta1 = pupil0%R(1)%thetaBounds(1)
    theta2 = pupil0%R(1)%thetaBounds(2)
    IF ( theta1 == theta2 ) THEN
    ! infinitesimal thin ring, no integration
      vecI(1) = K_r(theta1)
    ELSE
      vecI(1) = asqz(K_r, theta1, theta2, maxErr, maxDepth)! adaptive simpson
    END IF

    theta1 = pupil0%R(2)%thetaBounds(1)
    theta2 = pupil0%R(2)%thetaBounds(2)
    IF ( theta1 == theta2 ) THEN
      vecI(2) = K_phi(theta1)
    ELSE
      vecI(2) = asqz(K_phi, theta1, theta2, maxErr, maxDepth)
    END IF

    theta1 = pupil0%R(3)%thetaBounds(1)
    theta2 = pupil0%R(3)%thetaBounds(2)
    IF ( theta1 == theta2 ) THEN
      vecI(3) = K_z(theta1)
    ELSE
      vecI(3) = asqz(K_z, theta1, theta2, maxErr, maxDepth)
    END IF

  CONTAINS
    !----K_q----, the integrand over theta
    FUNCTION K_r(theta)
      REAL(KIND=dp), INTENT(IN)     :: theta
      COMPLEX(KIND=dp)              :: K_r

      K_r = g(theta) * SUM(matB(theta, (/1/), (/1, 2, 3/)))
    END FUNCTION K_r

    FUNCTION K_phi(theta)
      REAL(KIND=dp), INTENT(IN)     :: theta
      COMPLEX(KIND=dp)              :: K_phi

      K_phi = g(theta) * SUM(matB(theta, (/2/), (/1, 2, 3/)))
    END FUNCTION K_phi

    FUNCTION K_z(theta)
      REAL(KIND=dp), INTENT(IN)     :: theta
      COMPLEX(KIND=dp)              :: K_z

      K_z = g(theta) * SUM(matB(theta, (/3/), (/1, 2, 3/)))
    END FUNCTION K_z
    !----K_q----

    !----g(theta)----, the common part of the function over theta
    COMPLEX(KIND=dp) FUNCTION g(theta)
      REAL(KIND=dp), INTENT(IN)   :: theta

      g = SQRT(COS(theta))*SIN(theta)*fw(theta)*EXP((0,1)*k*z*COS(theta))
    END FUNCTION g

    ! the filling factor
    REAL(KIND=dp) FUNCTION fw(theta)
      REAL(KIND=dp), INTENT(IN)   :: theta

      fw = EXP(-(SIN(theta)**2/(f0**2*SIN(theta_max)**2)))
    END FUNCTION fw
    !----g(theta)----

    !----matB(theta, q1, q2)----
    ! e.g., matB(theta, 2, 1) represents the phi (2nd) focal field component that is
    ! due to the radial (1st) input field component
    FUNCTION matB(theta, q1, q2)
      REAL(KIND=dp), INTENT(IN)             :: theta
      INTEGER, DIMENSION(:), INTENT(IN)     :: q1, q2
      COMPLEX(KIND=dp), DIMENSION(1:3,1:3)  :: matB

      COMPLEX(KIND=dp), DIMENSION(1:3,1:3)            :: matJ, matP
      REAL(KIND=dp), DIMENSION(1:3,1:3)               :: matA
      REAL(KIND=dp)                                   :: x
      INTEGER                                         :: n, m, ir, l
      COMPLEX(KIND=dp), DIMENSION(1:3)                :: vecR
      COMPLEX(KIND=dp), DIMENSION(1:3,nMin:nMax)      :: vecPhi_cn
      LOGICAL, DIMENSION(nMin-1:nMax+1)               :: Jn_done
      COMPLEX(KIND=dp), DIMENSION(nMin-1:nMax+1)      :: Jn

      !--vecR(theta) and vecPhi_cn, loop for all three field components
      DO ir = 1, 3
        vecPhi_cn(ir,nMin:nMax) = pupil0%Phi(ir)%cn

        SELECT CASE ( TRIM(pupil0%R(ir)%fun) )
          CASE ('const')
            vecR(ir) = pupil0%R(ir)%const

          CASE ('sin')
            vecR(ir) = SIN(pupil0%R(ir)%sin%x*theta)**pupil0%R(ir)%sin%y

          CASE ('cos')
            vecR(ir) = COS(pupil0%R(ir)%cos%x*theta)**pupil0%R(ir)%cos%y

          CASE ('H1')
            vecR(ir) = (2.0_dp*f/w0)*SIN(theta)

          CASE ('LG')
            n = pupil0%R(ir)%lg%n
            l = pupil0%R(ir)%lg%l
            x = SQRT(2.0_dp)*f*SIN(theta)/w0
            vecR(ir) = SQRT(2.0_dp*factorial_n(n)/(pi*factorial_n(n+ABS(l)))) &
                      * (1.0_dp/w0) * (x**ABS(l)) * gLaguerrePoly(n, ABS(l), x**2)

          CASE DEFAULT
            WRITE(*,*) 'Not supported Radial function type!'
            STOP
        END SELECT
      END DO! ir

      !--A(theta): matrix representation for the aplanatic system
      matA = RESHAPE( (/COS(theta),   0.0_dp,   SIN(theta), &
                        0.0_dp,       1.0_dp,   0.0_dp,     &
                        -SIN(theta),  0.0_dp,   COS(theta)/),(/3,3/))

      matB = CMPLX(0,0,KIND=dp)
      matP = CMPLX(0,0,KIND=dp)

      !--Jn(k*rho*SIN(theta))
      Jn_done = .FALSE.
      Jn = CMPLX(0,0,KIND=dp)
      DO n = nMax, nMin, -1
        ! skip insignificant contributions
!        IF (  ALL((/vecPhi_cn(1,n)==CMPLX(0,0,KIND=dp),  &
!                    vecPhi_cn(2,n)==CMPLX(0,0,KIND=dp),  &
!                    vecPhi_cn(3,n)==CMPLX(0,0,KIND=dp)/)) ) THEN
        IF (  ALL((/ABS(vecPhi_cn(1,n)) .LT. epsilon_dp,  &
                    ABS(vecPhi_cn(2,n)) .LT. epsilon_dp,  &
                    ABS(vecPhi_cn(3,n)) .LT. epsilon_dp/)) ) THEN
          CYCLE
        END IF

        ! The Bessel function J_(n-1) and J_(n+1) are needed as well.
        DO m = n+1, n-1, -1
          ! Whether or not the Bessel function Jn(m) or Jn(-m) has been calculated.
          IF ( Jn_done(m) .EQV. .FALSE. ) THEN
            ! While inquiring Jn(-m), don't forget about legal the array bounds.
            IF ( (-m .GT. nMin) .AND. (-m .LT. nMax) ) THEN
              IF (Jn_done(-m) .EQV. .TRUE.) THEN
                Jn(m)=(-1)**(-m)*Jn(-m)
              ELSE ! Jn(-m) is legal, but not calculated
                ! new calculation for the Bessel function
                Jn(m) = besselj( m, k*rho*SIN(theta) )
              END IF
            ELSE ! Jn(-m) is out of bound
              ! new calculation for the Bessel function
              Jn(m) = besselj( m, k*rho*SIN(theta) )
            END IF

          END IF
          ! mark the Bessel function of mth-order Jn(m) as calculated
          Jn_done(m) = .TRUE.
        END DO! m

        !--when \vec{k} \vec{r} =kzcos(theta)-k\rho sin(theta)cos(\phi-\varphi)
        matJ = pi*(-(0,1)*EXP((0,1)*varphi))**n * RESHAPE( &
          (/(0,1)*(Jn(n-1)-Jn(n+1)), -(Jn(n-1)+Jn(n+1)),  CMPLX(0,0,KIND=dp), &
            (Jn(n-1)+Jn(n+1)),  (0,1)*(Jn(n-1)-Jn(n+1)),  CMPLX(0,0,KIND=dp), &
            CMPLX(0,0,KIND=dp),      CMPLX(0,0,KIND=dp),  2*Jn(n)/), (/3,3/) )

        !--when \vec{k} \vec{r} =kzcos(theta)+k\rho sin(theta)cos(\phi-\varphi)
!        matJ = pi*((0,1)*EXP((0,1)*varphi))**n * RESHAPE( &
!          (/-(0,1)*(Jn(n-1)-Jn(n+1)),  (Jn(n-1)+Jn(n+1)), CMPLX(0,0,KIND=dp), &
!            -(Jn(n-1)+Jn(n+1)), -(0,1)*(Jn(n-1)-Jn(n+1)), CMPLX(0,0,KIND=dp), &
!            CMPLX(0,0,KIND=dp),      CMPLX(0,0,KIND=dp),  2*Jn(n)/), (/3,3/) )

        matP(1,q2) = vecR(q2)*vecPhi_cn(q2,n)
        matP(2,q2) = matP(1,q2)
        matP(3,q2) = matP(1,q2)

        matB(q1,q2) = matB(q1,q2) + MATMUL( matJ(q1,:), matA(:,q2)*matP(:,q2) )
      END DO! n

    END FUNCTION matB
    !----matB(theta)----

  END FUNCTION vecI

END MODULE focal
