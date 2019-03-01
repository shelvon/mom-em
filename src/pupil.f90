! MODULE: pupil
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines to define a pupil function inserted in between the incident
! beam and focusing objective, i.e. right before focusing
! Here, fun_pupil = fun_aperture x fun_phase.

! Define the pupil function in cylindrical coordinate as P(theta, phi)

! For an analytical pupil function, the pupil function writes as
! P(theta, phi)=A(theta, phi)*Phi(theta, phi) with
! A(theta, phi) the aperture/amplitude function;
! Phi(theta, phi) the phase modulation function.

!!!! to be done:
! For numerical pupil function:  the pupil function is a
! two-dimensional matrix with row->theta; column->phi.
MODULE pupil
  USE data
  USE constants

  IMPLICIT NONE

CONTAINS

  FUNCTION circ(rho) RESULT(res)
    REAL (KIND=dp) :: rho
    REAL :: res

    IF ( (rho <= 1.0) .AND. (rho >= 0.0) ) THEN
      res = 1.0
    ELSE
      res = 0.0
    END IF
  END FUNCTION circ

  FUNCTION rect(t) RESULT(res)
    REAL (KIND=dp) :: t
    REAL :: res

    IF ( (t > 0.5) .OR. (t < -0.5) ) THEN
      res = 0
    ELSE IF ( (t == 0.5) .OR. (t == -0.5) ) THEN
      res = 0.5
    ELSE IF ( (t < 0.5) .AND. (t > -0.5) ) THEN
      res = 1
    END IF

  END FUNCTION rect

  FUNCTION rectGeneral(t, a, b) RESULT(res)
    REAL (KIND=dp) :: t, a, b
    REAL :: res

    res = rect((t-a)/(b-a)-0.5)
  END FUNCTION rectGeneral

  FUNCTION rectPhase(t,intervals) RESULT(res)
    REAL (KIND=dp) :: t
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: intervals
    INTEGER :: n, intervals_n
    REAL :: res

    res = 0.0

    intervals_n=SIZE(intervals,2)
    DO n=1,intervals_n
      res = res + rectGeneral(t, intervals(n,1), intervals(n,2))
    END DO
  END FUNCTION rectPhase

  ! Get fourier series coefficients from equations (3.79) and (3.80) in
  ! Oppenheim, Alan V., Alan S. Willsky and Syed Hamid Nawab.
  ! 1997. Signals & Systems. 2nd ed. Prentice-Hall Signal Processing Series.
  ! Upper Saddle River, N.J: Prentice Hall.
  FUNCTION coeffSquareWave(T, T1, t0, k) RESULT(ak)
    REAL (KIND=dp)      :: T, T1, t0
    INTEGER             :: k
    COMPLEX (KIND=dp)   :: ak

    IF(k==0) THEN
      ak = 2*T1/T
    ELSE
      ak = SIN(k*2*pi*T1/T)/(k*pi)
    END IF
    ! time shift
    ak = ak*EXP(-(0,1)*k*2*pi*t0/T)
  END FUNCTION coeffSquareWave

  FUNCTION coeffRectPhase(pupil,j) RESULT(dj)
    TYPE(pupil_type), INTENT(IN)    :: pupil
    INTEGER, INTENT(IN)             :: j
    INTEGER                         :: n
    REAL (KIND=dp)                  :: aj, bj, delta, Ak
    COMPLEX (KIND=dp)               :: cj
    REAL (KIND=dp), DIMENSION(1:2)  :: dj
    REAL (KIND=dp), DIMENSION(1:2,pupil%phase%rect%intervals_n) :: intervals

    intervals = pi*pupil%phase%rect%intervals ! convert to unit in [\pi]
    IF(pupil%phase%type=='rect') THEN
    ! coefficients for 'rect' type phase mask
      !ck = (0.0_dp, 0.0_dp)
      !DO n=1,pupil%phase%rect%intervals_n
      !  aj = intervals(1,n)
      !  bj = intervals(2,n)
      !  cj = cj + coeffSquareWave(2*pi,bj-aj,(bj+aj)/2,j)
      !END DO
      !IF(j/=0) THEN
      !  delta = ATAN2(AIMAG(cj),REAL(cj))
      !  Ak = ABS(cj)
      !  aj = 2*Ak*COS(delta)
      !  bj = -2*Ak*SIN(delta)
      !END IF
      IF(j==0) THEN
        aj = SUM(intervals(2,:) - intervals(1,:))
        bj = 0.0_dp
      ELSE
        aj = SUM(SIN(j*intervals(2,:)) - SIN(j*intervals(1,:)))
        aj = (1.0_dp/REAL(j,dp))*aj
        bj = SUM(COS(j*intervals(2,:)) - COS(j*intervals(1,:)))
        bj = -(1.0_dp/REAL(j,dp))*bj
      END IF
    END IF

    ! WRITE(*,*) tol
    IF(ABS(aj) .LT. tol) THEN
      aj = 0.0_dp
    END IF
    IF(ABS(bj) .LT. tol) THEN
      bj = 0.0_dp
    END IF

    dj=(1.0_dp/pi)*(/aj,bj/)
  END FUNCTION coeffRectPhase

  SUBROUTINE pupil_coeff(pupil)
    TYPE(pupil_type), INTENT(INOUT)     :: pupil

    COMPLEX(KIND=dp)                    :: prefix
    INTEGER                             :: n, charge, l, id, intervals_n
    REAL(KIND=dp), DIMENSION(1:2)       :: delta ! amount of phase shift [\pi]
    REAL(KIND=dp), DIMENSION(1:2, pupil%phase%rect%intervals_n)  &
                                        :: intervals, intervalsLower
    REAL(KIND=dp), DIMENSION(1:2, pupil%phase%petal_rect%intervals_n)  &
                                        :: intervalsPetal

    pupil%cn = CMPLX(0.0_dp,0.0_dp) ! initialized to zeros
    DO n = -pupil%phase%nMax, pupil%phase%nMax
    ! Cylindrical Vector Beams  : CVBs
    ! State of Polarization     : SOP

      ! Hermite Gaussian (HG00) beam, i.e. the basis
      IF ( TRIM(pupil%phase%type) == 'hg00') THEN
        IF ( n == 0 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, 0.0_dp)
        END IF

      ! Hermite Gaussian (HG01) beam
      ELSE IF ( TRIM(pupil%phase%type) == 'hg01') THEN
        IF ( n == 1 ) THEN
          pupil%cn(n) = CMPLX(0.0_dp, -1.0_dp/2.0_dp)
        ELSE IF ( n == -1 ) THEN
          pupil%cn(n) = CMPLX(0.0_dp, 1.0_dp/2.0_dp)
        END IF

      ! Hermite Gaussian (HG10) beam
      ELSE IF ( TRIM(pupil%phase%type) == 'hg10' ) THEN
        IF ( ABS(n) == 1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp/2.0_dp, 0.0_dp)
        END IF

      ! CVBs, local SOP: azimuthally polarized
      ELSE IF ( TRIM(pupil%phase%type) == 'azi') THEN
        IF ( n == 1 ) THEN
          pupil%cn(n) = CMPLX(0.0_dp, -1.0_dp/2.0_dp)
        ELSE IF ( n == -1 ) THEN
          pupil%cn(n) = CMPLX(0.0_dp, 1.0_dp/2.0_dp)
        END IF

      ! CVBs, local SOP: radially polarized
      ELSE IF ( TRIM(pupil%phase%type) == 'rad') THEN
        IF ( ABS(n) == 1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp/2.0_dp, 0.0_dp)
        END IF

      ! CVBs, local SOP: m45, i.e. minus 45 [deg] polarized
      ELSE IF ( TRIM(pupil%phase%type) == 'm45') THEN
        IF ( n == 1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, -1.0_dp)/2.0_dp
        ELSE IF ( n == -1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, 1.0_dp)/2.0_dp
        END IF

      ! CVBs, local SOP: p45, i.e. plus 45 [deg] polarized
      ELSE IF ( TRIM(pupil%phase%type) == 'p45') THEN
        IF ( n == 1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, 1.0_dp)/2.0_dp
        ELSE IF ( n == -1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, -1.0_dp)/2.0_dp
        END IF

      ! CVBs, local SOP: right-circular polarization
      ELSE IF ( TRIM(pupil%phase%type) == 'rcp') THEN
        IF ( n == 1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, 0.0_dp)
        END IF

      ! CVBs, local SOP: left-circular polarization
      ELSE IF ( TRIM(pupil%phase%type) == 'lcp') THEN
        IF ( n == -1 ) THEN
          pupil%cn(n) = CMPLX(1.0_dp, 0.0_dp)
        END IF

      ! vortex beam
      ELSE IF ( TRIM(pupil%phase%type) == 'vortex' ) THEN
        IF ( n == pupil%phase%vortex%charge ) THEN
          pupil%cn(n) = CMPLX(1.0_dp,0.0_dp)
        END IF

      ELSE IF ( TRIM(pupil%phase%type) == 'petal' ) THEN
        IF ( ABS(n) == pupil%phase%petal%l ) THEN
          pupil%cn(n) = CMPLX(1.0_dp/2.0_dp,0.0_dp)
        END IF

      ! bessel beam
      ELSE IF ( TRIM(pupil%phase%type) == 'bessel' ) THEN
        IF ( TRIM(pupil%phase%bessel%input) == 'azi' ) THEN
          IF ( n == 1 ) THEN
            pupil%cn(n) = CMPLX(0.0_dp, -1.0_dp/2.0_dp)
          ELSE IF ( n == -1 ) THEN
            pupil%cn(n) = CMPLX(0.0_dp, 1.0_dp/2.0_dp)
          END IF
        ELSE IF ( TRIM(pupil%phase%bessel%input) == 'rad' ) THEN
          IF ( ABS(n) == 1 ) THEN
            pupil%cn(n) = CMPLX(1.0_dp/2.0_dp, 0.0_dp)
          END IF
        END IF

      ! test for petal with binary phase mask
      ELSE IF ( TRIM(pupil%phase%type) == 'petal_rect' ) THEN
        l = pupil%phase%petal_rect%l
        intervals_n = pupil%phase%petal_rect%intervals_n
        intervalsPetal = pupil%phase%petal_rect%intervals
        IF ( n == l ) THEN
          pupil%cn(n) = (1.0_dp/(2*pi*2))*  &
            ( ( intervalsPetal(2,1) - intervalsPetal(1,1) ) - &
            ( EXP(-(0,1)*(l+n)*intervalsPetal(2,1)) -         &
              EXP(-(0,1)*(l+n)*intervalsPetal(1,1)) )/((0,1)*(l+n)) )
        ELSE IF ( n == -l ) THEN
          pupil%cn(n) = (1.0_dp/(2*pi*2))*  &
            ( ( intervalsPetal(2,1) - intervalsPetal(1,1) ) + &
            ( EXP((0,1)*(l-n)*intervalsPetal(2,1)) -       &
              EXP((0,1)*(l-n)*intervalsPetal(1,1)) )/((0,1)*(l-n)) )
        ELSE
          pupil%cn(n) = (1.0_dp/(2*pi*2))*  &
          ( -( EXP(-(0,1)*(l+n)*intervalsPetal(2,1)) -  &
               EXP(-(0,1)*(l+n)*intervalsPetal(1,1)) )/((0,1)*(l+n)) &
            +( EXP((0,1)*(l-n)*intervalsPetal(2,1)) -   &
               EXP((0,1)*(l-n)*intervalsPetal(1,1)) )/((0,1)*(l-n)) )
        END IF

      ELSE IF ( TRIM(pupil%phase%type) == 'rect' ) THEN
        intervals_n = pupil%phase%rect%intervals_n
        intervals = pi*pupil%phase%rect%intervals ! convert to unit in [\pi]
        delta = pi*pupil%phase%rect%delta ! convert to unit in [\pi]

        ! phase modulation higher level
        prefix = (1.0_dp/2/pi)*EXP((0,1)*delta(2))
        IF( n == 0 ) THEN
          pupil%cn(n) = prefix*SUM(intervals(2,:) - intervals(1,:))
        ELSE
          pupil%cn(n) = prefix*(0,1)/n * &
            SUM(EXP(-(0,1)*n*intervals(2,:)) - EXP(-(0,1)*n*intervals(1,:)))
        END IF

        ! phase modulation lower level
        prefix = (1.0_dp/2/pi)*EXP((0,1)*delta(1))
        intervalsLower(1,:) = (/intervals(2,1:intervals_n)/)
        intervalsLower(2,:) = (/intervals(1,2:intervals_n), &
                                intervals(1,1)+2*pi/)
        IF ( n == 0 ) THEN
          pupil%cn(n) = pupil%cn(n) + &
            prefix*SUM(intervalsLower(2,:) - intervalsLower(1,:))
        ELSE
          pupil%cn(n) = pupil%cn(n) + prefix*(0,1)/n * &
            SUM(  EXP(-(0,1)*n*intervalsLower(2,:)) -     &
                  EXP(-(0,1)*n*intervalsLower(1,:)) )
        END IF

      END IF

      ! remove insignificant orders
      IF ( ABS(REAL(pupil%cn(n))) .LT. tol ) THEN
        pupil%cn(n) = CMPLX(0.0_dp, AIMAG(pupil%cn(n)))
      END IF
      IF ( ABS(AIMAG(pupil%cn(n))) .LT. tol ) THEN
        pupil%cn(n) = CMPLX(REAL(pupil%cn(n)), 0.0_dp)
      END IF
    END DO
  END SUBROUTINE pupil_coeff

END MODULE pupil
