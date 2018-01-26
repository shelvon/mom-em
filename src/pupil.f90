! MODULE: pupil
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines to define a pupile function inserted in between
! the incident beam and focusing objective.
MODULE pupil
  USE constants

  IMPLICIT NONE

  ! Define the pupil function in cylindrical coordinate as P(theta, phi)

  ! For an analytical pupile function, the pupil function writes as
  ! P(theta, phi)=A(theta, phi)*exp(i*Phi(theta, phi)) with
  ! A(theta, phi) the aperture/amplitude function;
  ! Phi(theta, phi) the phase modulation function.

  ! For numerical pupile function:  the pupil function is a
  ! two-dimensional matrix with row->theta; column->phi.

  TYPE data_circ
    REAL (KIND=dp)      :: r=1  ! radius of circ function, normalized to
                                ! the waist radius of the incident beam
  END TYPE data_circ

  ! intervals where function is equal to 1, which is also
  ! where the phase modulation applies.
  TYPE data_rect
    REAL (KIND=dp), DIMENSION(1:2)              :: delta! amount of phase shift [\pi]
    INTEGER                                     :: intervals_n
    ! REAL (KIND=dp), DIMENSION(2) :: intervals
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: intervals
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: dj
    ! with (1,:) stores aj and (2,:) bj
    ! aj for cos(j\phi) and bj for sin(j\phi)
  END TYPE data_rect

  ! topological charge of vortex beam
  TYPE data_petal
    INTEGER :: l ! l is the integer subindex in Petal_{p,l} = LG_{p,l} + LG_{p,-l}
  END TYPE data_petal

  ! topological charge of vortex beam
  TYPE data_vortex
    INTEGER :: charge! 'charge x 2\pi' phase change over one rotation in azimuthal angle
  END TYPE data_vortex

  TYPE data_aperture
    ! type of pupil function: either analytical or numerical.
    CHARACTER (LEN=32)  :: type
    TYPE(data_circ)     :: circ
  END TYPE data_aperture

  TYPE data_phase
    ! type of pupil function: either analytical or numerical.
    CHARACTER (LEN=32)      :: type
    TYPE(data_rect)         :: rect
    TYPE(data_vortex)       :: vortex !spiral phase plate
    TYPE(data_petal)        :: petal  !petal beam phase term
  END TYPE data_phase

  TYPE data_pupil
    ! type of pupil function: either analytical or numerical.
    CHARACTER (LEN=32)  :: type

    ! aperture function
    TYPE(data_aperture) :: aperture

    ! phase function
    TYPE(data_phase)    :: phase

    ! complex-valued Fourier expansion coefficient
    ! by using complex Fourier series
    ! Works only for pupil function that is independent of radial variable, \rho
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: cj
    ! cj = 1/(2*pi)*\Int_{-pi}^{pi} e^{i\Phi(\phi)}e^{-im\phi} d\phi
  END TYPE data_pupil

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

!  FUNCTION circAperture(rho,waist) RESULT(res)
!    REAL (KIND=dp) :: rho
!    REAL (KIND=dp) :: waist
!    REAL :: res
!
!    res = circ(rho/waist)
!  END FUNCTION circAperture

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

  ! get fourier series coefficients from equations (3.79) and (3.80) in
  ! Oppenheim, Alan V., Alan S. Willsky and Syed Hamid Nawab. 1997. Signals & Systems.
  ! 2nd ed. Prentice-Hall Signal Processing Series. Upper Saddle River, N.J: Prentice Hall.
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
    TYPE(data_pupil), INTENT(IN)    :: pupil
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

  FUNCTION coeffPupil(pupil,j) RESULT(cj)
    TYPE(data_pupil), INTENT(IN)    :: pupil
    INTEGER, INTENT(IN)             :: j
    COMPLEX (KIND=dp)               :: cj

    COMPLEX (KIND=dp)               :: prefix
    INTEGER                         :: charge, l, id, intervals_n
    REAL (KIND=dp), DIMENSION(1:2)  :: delta! amount of phase shift [\pi]
    REAL (KIND=dp), DIMENSION(1:2,pupil%phase%rect%intervals_n) :: intervals,intervalsLower

    intervals_n = pupil%phase%rect%intervals_n
    intervals = pi*pupil%phase%rect%intervals ! convert to unit in [\pi]

    IF(pupil%phase%type=='vortex') THEN
      charge = pupil%phase%vortex%charge
      IF(j==charge) THEN
        cj = CMPLX(1.0_dp,0.0_dp)
      ELSE
        cj = CMPLX(0.0_dp,0.0_dp)
      END IF
    ELSE IF(pupil%phase%type=='petal') THEN
      l = pupil%phase%petal%l
      IF(ABS(j)==l) THEN
        cj = CMPLX(1.0_dp,0.0_dp)
      ELSE
        cj = CMPLX(0.0_dp,0.0_dp)
      END IF
    ELSE IF(pupil%phase%type=='rect') THEN
      delta = pi*pupil%phase%rect%delta ! convert to unit in [\pi]

      ! phase modulation higher level
      prefix = (1.0_dp/2/pi)*EXP((0,1)*delta(2))
      IF(j==0) THEN
        cj = prefix*SUM(intervals(2,:) - intervals(1,:))
      ELSE
        cj = prefix/(-(0,1)*j)*SUM(EXP(-(0,1)*j*intervals(2,:)) - EXP(-(0,1)*j*intervals(1,:)))
      END IF

      ! phase modulation lower level
      prefix = (1.0_dp/2/pi)*EXP((0,1)*delta(1))
      intervalsLower(1,:) = (/intervals(2,1:intervals_n)/)
      intervalsLower(2,:) = (/intervals(1,2:intervals_n), intervals(1,1)+2*pi/)
      IF(j==0) THEN
        cj = cj+prefix*SUM(intervalsLower(2,:) - intervalsLower(1,:))
      ELSE
        cj = cj+prefix/(-(0,1)*j)*SUM(EXP(-(0,1)*j*intervalsLower(2,:)) - EXP(-(0,1)*j*intervalsLower(1,:)))
      END IF
    END IF

    IF(REAL(cj) .LT. tol) THEN
      cj = CMPLX(0.0_dp,AIMAG(cj))
    END IF
    IF(AIMAG(cj) .LT. tol) THEN
      cj = CMPLX(REAL(cj), 0.0_dp)
    END IF
  END FUNCTION coeffPupil

END MODULE pupil
