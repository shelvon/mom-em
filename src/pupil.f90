! MODULE: pupil
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Routines to define the pupil function for creating general input beam
! from the basic input beam profile.
! The radial and azimuthal function with index 1->r, 2->phi, 3->z
! are applied to E_r,phi,z separately.
MODULE pupil
  USE data
  USE constants

  IMPLICIT NONE

CONTAINS

  FUNCTION pupil_e2h(pE) RESULT(pH)
    TYPE(pupil_type), INTENT(IN)      :: pE

    TYPE(pupil_type)                  :: pH

    pH = pE

    pH%Phi(2) = pE%Phi(1)
    pH%R(2) = pE%R(1)
    pH%Phi(1) = pE%Phi(2)
    pH%R(1) = pE%R(2)
    pH%Phi(1)%cn = -pH%Phi(1)%cn

!    WRITE(*,*) LBOUND(pH%Phi(1)%cn), pH%Phi(1)%cn
!    WRITE(*,*) LBOUND(pE%Phi(2)%cn), pE%Phi(2)%cn
!    STOP

  END FUNCTION pupil_e2h

  ! generate a general and equivalent pupil function
  FUNCTION generate_pupil(focus, theta_max) RESULT(p0)
    TYPE(focus_type), INTENT(IN)      :: focus
    REAL(KIND=dp), INTENT(IN)         :: theta_max
    TYPE(pupil_type)                  :: p0

    INTEGER                           :: nMin, nMax

    ! initialize pupil_type::pupil
    SELECT CASE ( TRIM(focus%input) )
      CASE ('pupil')
        ! get pupil function from user input directly
        p0 = focus%pupil
        !--generate cn values if the functional form is not series
        CALL cn_series(p0%Phi(1))
        CALL cn_series(p0%Phi(2))
        CALL cn_series(p0%Phi(3))

!      CASE ('invbasis')
!        p0 = invbasis2pupil(focus%invbasis%pol, focus%invbasis%n)

      CASE ('basis')
        p0 = basis2pupil(focus%basis%pol, focus%basis%n)

      CASE ('cos')
        p0 = cos2pupil(focus%cos%pol, focus%cos%l)

      CASE ('sin')
        p0 = sin2pupil(focus%sin%pol, focus%sin%l)

      CASE ('hg')
        p0 = hg2pupil(focus%hg%pol, focus%hg%m, focus%hg%n)

      CASE ('petal')
        p0 = petal2pupil(focus%petal%pol, focus%petal%n, focus%petal%l)

      CASE ('lg')
        p0 = lg2pupil(focus%lg%pol, focus%lg%n, focus%lg%l)

      CASE ('rect')
        p0 = rect2pupil(focus%rect)

      CASE ('bessel')
        p0 = basis2pupil(focus%bessel%pol, focus%bessel%j)

      CASE DEFAULT
        WRITE(*,*) "Not supported input beam for focal field calculation,"
        WRITE(*,*) "please try to use 'pupil' for the input beam!"
        STOP
    END SELECT

    !--get the lower and upper bounds of the harmonics in all three Phi functions
    nMin = MINVAL((/LBOUND(p0%Phi(1)%cn),LBOUND(p0%Phi(2)%cn),LBOUND(p0%Phi(3)%cn)/))
    nMax = MAXVAL((/UBOUND(p0%Phi(1)%cn),UBOUND(p0%Phi(2)%cn),UBOUND(p0%Phi(3)%cn)/))
    ! synchronize the cn's indices in all three Phi functions
    CALL cn_pad(p0%Phi(1)%cn, nMin, nMax)
    CALL cn_pad(p0%Phi(2)%cn, nMin, nMax)
    CALL cn_pad(p0%Phi(3)%cn, nMin, nMax)

    !--bounds of the integration over theta
    IF ( TRIM(focus%input)=='bessel' ) THEN
      p0%R(1)%thetaBounds = focus%bessel%theta*deg2rad  &
                          + focus%bessel%delta*(/-0.5_dp, 0.5_dp/)*deg2rad
      p0%R(2)%thetaBounds = p0%R(1)%thetaBounds
      p0%R(3)%thetaBounds = p0%R(1)%thetaBounds
    ELSE
      p0%R(1)%thetaBounds = (/0.0_dp, theta_max/)
      p0%R(2)%thetaBounds = (/0.0_dp, theta_max/)
      p0%R(3)%thetaBounds = (/0.0_dp, theta_max/)
    END IF

  END FUNCTION generate_pupil

  !--general pupil function for input beam of basic types
  RECURSIVE FUNCTION basis2pupil(pol, n) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: n
    TYPE(pupil_type)              :: res

    SELECT CASE (TRIM(pol))
      CASE ('r', 'phi')
        res%R(:)%fun = 'const'
        res%R(:)%const = CMPLX(1, 0, KIND=dp)
        res%Phi(:)%fun = 'series'
        ALLOCATE(res%Phi(1)%cn(n:n))
        ALLOCATE(res%Phi(2)%cn(n:n))
        ALLOCATE(res%Phi(3)%cn(n:n))
        res%Phi(1)%cn = CMPLX(0,0,KIND=dp)
        res%Phi(2)%cn = CMPLX(0,0,KIND=dp)
        res%Phi(3)%cn = CMPLX(0,0,KIND=dp)
        IF (TRIM(pol)=='r') THEN
          res%Phi(1)%cn = CMPLX(1,0,KIND=dp)
        ELSE IF (TRIM(pol)=='phi') THEN
          res%Phi(2)%cn = CMPLX(1,0,KIND=dp)
        END IF

      CASE ('+')
        ! |+> = 1/sqrt(2)*(|r> + i|phi>)
        res = basis2pupil('r',n)
        res = pupil_add( res, basis2pupil('phi',n), &
                        (1.0,0.0)/SQRT(2.0), (0.0,1.0)/SQRT(2.0) )

      CASE ('-')
        ! |-> = 1/sqrt(2)*(|r> - i|phi>)
        res = basis2pupil('r',n)
        res = pupil_add( res, basis2pupil('phi',n), &
                        (1.0,0.0)/SQRT(2.0), (0.0,-1.0)/SQRT(2.0) )

      CASE ('x')
        res = basis2pupil('r',n-1)
        res = pupil_add( res, basis2pupil('r',n+1), (0.5,0.0), (0.5,0.0) )
        res = pupil_add( res, basis2pupil('phi',n-1), (1.0,0.0), (0.0,-0.5) )
        res = pupil_add( res, basis2pupil('phi',n+1), (1.0,0.0), (0.0,0.5) )

      CASE ('y')
        res = basis2pupil('r',n-1)
        res = pupil_add( res, basis2pupil('r',n+1), (0.0,0.5), (0.0,-0.5) )
        res = pupil_add( res, basis2pupil('phi',n-1), (1.0,0.0), (0.5,0.0) )
        res = pupil_add( res, basis2pupil('phi',n+1), (1.0,0.0), (0.5,0.0) )

      CASE DEFAULT
        WRITE(*,*) "Not supported pol for the basis input beam!"
        WRITE(*,*) "please try to use 'pupil' for the input beam!"
        STOP
    END SELECT

  END FUNCTION basis2pupil

!  FUNCTION invbasis2pupil(pol, l) RESULT(res)
!    CHARACTER(LEN=*), INTENT(IN)  :: pol
!    INTEGER, INTENT(IN)           :: l
!    TYPE(pupil_type)              :: res
!
!    res = basis2pupil(pol,l)
!    res%R(1)
!  END FUNCTION invbasis2pupil

  FUNCTION cos2pupil(pol, l) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: l
    TYPE(pupil_type)              :: res

    res = pupil_add( basis2pupil(pol,l),basis2pupil(pol,-l),(0.5,0.0),(0.5,0.0) )
  END FUNCTION cos2pupil

  FUNCTION sin2pupil(pol, l) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: l
    TYPE(pupil_type)              :: res

    res = pupil_add( basis2pupil(pol,l),basis2pupil(pol,-l),(0.0,-0.5),(0.0,-0.5) )
  END FUNCTION sin2pupil

  FUNCTION lg2pupil(pol, n, l) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: n,l
    TYPE(pupil_type)              :: res

    res = basis2pupil(pol,-l)! orbital angular momentum EXP(-il\varphi)
    res%R(:)%fun = 'LG'
    res%R(:)%lg%n = n
    res%R(:)%lg%l = l
  END FUNCTION lg2pupil

  FUNCTION petal2pupil(pol, n, l) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: n,l
    TYPE(pupil_type)              :: res

    res = pupil_add( basis2pupil(pol,l),basis2pupil(pol,-l),(0.5,0.0),(0.5,0.0) )
    res%R(:)%fun = 'LG'
    res%R(:)%lg%n = n
    res%R(:)%lg%l = l
  END FUNCTION petal2pupil

  FUNCTION hg2pupil(pol, m, n) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: pol
    INTEGER, INTENT(IN)           :: m,n
    TYPE(pupil_type)              :: res

    IF ( m==0 .AND. n==0 ) THEN
      res = basis2pupil(pol,0)
    ELSE IF ( m==1 .AND. n==0 ) THEN
      res = cos2pupil(pol, 1)
      res%R(:)%fun = 'H1'
    ELSE IF ( m==0 .AND. n==1 ) THEN
      res = sin2pupil(pol, 1)
      res%R(:)%fun = 'H1'
    ELSE
      WRITE(*,*) "Higher order HG modes are not supported!"
      WRITE(*,*) "please try to use 'pupil' for the input beam!"
      STOP
    END IF
  END FUNCTION hg2pupil

  RECURSIVE FUNCTION rect2pupil(rect) RESULT(res)
    TYPE(rect_type), INTENT(IN)   :: rect
    TYPE(pupil_type)              :: res

    TYPE(rect_type)                   :: rect_r, rect_phi
    TYPE(pupil_type), DIMENSION(-1:1) :: pupil_r, pupil_phi

    SELECT CASE ( TRIM(rect%pol) )
      CASE ('r')
        res = basis2pupil('r',0)
        res%Phi(1)%fun = 'rect'
        res%Phi(1)%rect = rect
        CALL cn_series(res%Phi(1))

      CASE ('phi')
        res = basis2pupil('phi',0)
        res%Phi(2)%fun = 'rect'
        res%Phi(2)%rect = rect
        CALL cn_series(res%Phi(2))

      CASE ('x','y')
        rect_r = rect; rect_r%pol = 'r';
        pupil_r(0) = rect2pupil(rect_r)
        pupil_r(-1) = pupil_r(0)
        pupil_r(1) = pupil_r(0)
        CALL cn_shift( pupil_r(1)%Phi(1)%cn, +1 )
        CALL cn_shift( pupil_r(-1)%Phi(1)%cn, -1 )

        rect_phi = rect; rect_phi%pol = 'phi';
        pupil_phi(0) = rect2pupil(rect_phi)
        pupil_phi(-1) = pupil_phi(0)
        pupil_phi(1) = pupil_phi(0)
        CALL cn_shift( pupil_phi(1)%Phi(2)%cn, +1 )
        CALL cn_shift( pupil_phi(-1)%Phi(2)%cn, -1 )

        IF ( TRIM(rect%pol)=='x' ) THEN
          res = pupil_add( pupil_add(pupil_r(1),pupil_r(-1),(0.5,0.0),(0.5,0.0)), &
                  pupil_add(pupil_phi(1),pupil_phi(-1),(0.0,0.5),(0.0,-0.5)) )
        ELSE IF ( TRIM(rect%pol)=='y' ) THEN
          res = pupil_add( pupil_add(pupil_r(1),pupil_r(-1),(0.0,-0.5),(0.0,0.5)), &
                  pupil_add(pupil_phi(1),pupil_phi(-1),(0.5,0.0),(0.5,0.0)) )
        END IF

      CASE DEFAULT
        WRITE(*,*) "Not supported pol for the rect input beam!"
        WRITE(*,*) "please try to use 'pupil' for the input beam!"
        STOP
    END SELECT
  END FUNCTION rect2pupil

  RECURSIVE FUNCTION pupil_add(b1, b2, w1_set, w2_set) RESULT(b0)
    TYPE(pupil_type), INTENT(IN)  :: b1, b2
    COMPLEX, INTENT(IN), OPTIONAL :: w1_set, w2_set! corresponding weights
    TYPE(pupil_type)              :: b0

    INTEGER, DIMENSION(1:2)       :: n1, n2, n0
    INTEGER                       :: ir
    COMPLEX                       :: w1, w2

    IF ( PRESENT(w1_set) ) THEN
      w1 = w1_set
    ELSE
      w1 = (1.0,0.0)
    END IF

    IF ( PRESENT(w2_set) ) THEN
      w2 = w2_set
    ELSE
      w2 = (1.0,0.0)
    END IF

    DO ir = 1, 3
      ! the R%fun and Phi%fun must be the same
      IF ( ( TRIM(b1%R(ir)%fun)==TRIM(b2%R(ir)%fun) ) .AND. &
           ( TRIM(b1%Phi(ir)%fun)==TRIM(b2%Phi(ir)%fun) ) ) THEN

        b0%R(ir) = b1%R(ir)
        b0%Phi(ir)%fun = b1%Phi(ir)%fun

        n1 = (/LBOUND(b1%Phi(ir)%cn,1),UBOUND(b1%Phi(ir)%cn,1)/)
        n2 = (/LBOUND(b2%Phi(ir)%cn,1),UBOUND(b2%Phi(ir)%cn,1)/)
        n0 = (/MIN(n1(1),n2(1)), MAX(n1(2),n2(2))/)

        ALLOCATE( b0%Phi(ir)%cn(n0(1):n0(2)) )
        b0%Phi(ir)%cn = CMPLX(0,0,KIND=dp)
        b0%Phi(ir)%cn(n1(1):n1(2)) = b0%Phi(ir)%cn(n1(1):n1(2)) &
                                   + w1*b1%Phi(ir)%cn(n1(1):n1(2))
        b0%Phi(ir)%cn(n2(1):n2(2)) = b0%Phi(ir)%cn(n2(1):n2(2)) &
                                   + w2*b2%Phi(ir)%cn(n2(1):n2(2))

      ELSE
        WRITE(*,*) 'illegal addition of the pupil functions!'
        STOP
      END IF

    END DO! ir

  END FUNCTION pupil_add

  ! calculate cn-values for Phi_q when its functional form is not 'series'.
  SUBROUTINE cn_series(Phi)
    TYPE(funPhi_type), INTENT(INOUT)  :: Phi

    COMPLEX(KIND=dp)                    :: prefix
    INTEGER                             :: i, n
    REAL(KIND=dp)                       :: delta, a, b

!    IF ( TRIM(Phi%fun)=='series' ) THEN
!      WRITE(*,*) '    Phi%fun==series, cn-values is taken from user input.'
!
!    ELSE IF ( TRIM(Phi%fun)=='rect' ) THEN
    IF ( TRIM(Phi%fun)=='rect' ) THEN
      IF ( ALLOCATED(Phi%cn) ) THEN
        DEALLOCATE(Phi%cn)
      END IF
      ALLOCATE( Phi%cn(-Phi%rect%series_n:Phi%rect%series_n) )
      Phi%cn = CMPLX(0,0,KIND=dp)

      DO i = 1, SIZE(Phi%rect%delta)
        ! multiply pi to convert unit: [\pi] --> [rad]
        delta = pi * Phi%rect%delta(i)
        a = pi * Phi%rect%intervals(1,i)
        b = pi * Phi%rect%intervals(2,i)
        prefix = ( EXP((0,1)*delta)-1.0_dp )/(2.0_dp*pi)

        DO n = -Phi%rect%series_n, Phi%rect%series_n
          IF( n == 0 ) THEN
            Phi%cn(n) = Phi%cn(n) + prefix*(b-a) + 1.0_dp
          ELSE
            Phi%cn(n) = Phi%cn(n) + &
              prefix*(0,1)/n*(EXP(-(0,1)*n*b)-EXP(-(0,1)*n*a))
          END IF
        END DO!n

      END DO! i

    END IF

    ! remove insignificant orders
    DO n = LBOUND(Phi%cn,1), UBOUND(Phi%cn,1)
      IF ( ABS(REAL(Phi%cn(n))) .LT. epsilon_dp ) THEN
        Phi%cn(n) = CMPLX(0.0_dp, AIMAG(Phi%cn(n)), KIND=dp)
      END IF
      IF ( ABS(AIMAG(Phi%cn(n))) .LT. epsilon_dp ) THEN
        Phi%cn(n) = CMPLX(REAL(Phi%cn(n)), 0.0_dp, KIND=dp)
      END IF
    END DO! n

    ! update fun type to series
    Phi%fun='series'

  END SUBROUTINE cn_series

  !--extend cn size by padding zero
  SUBROUTINE cn_pad(cn, nMin, nMax)
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: cn
    INTEGER, INTENT(IN)                           :: nMin, nMax

    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: cn_tmp

    IF ( LBOUND(cn,1) .LT. nMin ) THEN
      WRITE(*,*) 'The lower bound of cn cannot be smaller than nMin.'
      STOP
    END IF
    IF ( UBOUND(cn,1) .GT. nMax ) THEN
      WRITE(*,*) 'The upper bound of cn cannot be larger than nMax.'
      STOP
    END IF

    ALLOCATE(cn_tmp(LBOUND(cn,1):UBOUND(cn,1)))
    cn_tmp = cn
    DEALLOCATE(cn)
    ALLOCATE(cn(nMin:nMax))
    cn = CMPLX(0,0,KIND=dp)
    cn(LBOUND(cn_tmp,1):UBOUND(cn_tmp,1)) = cn_tmp
    DEALLOCATE(cn_tmp)

  END SUBROUTINE cn_pad

  !--shift cn values from index n to index nShift for each element
  SUBROUTINE cn_shift(cn, nShift)
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: cn
    INTEGER, INTENT(IN)                           :: nShift

    INTEGER                                       :: nMin, nMax
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: cn_tmp

    nMin=LBOUND(cn,1)
    nMax=UBOUND(cn,1)
    ALLOCATE( cn_tmp(nMin:nMax) )
    cn_tmp = cn
    DEALLOCATE(cn)
    ALLOCATE( cn(nMin+nShift:nMax+nShift) )
    cn(nMin+nShift:nMax+nShift) = cn_tmp(nMin:nMax)
    DEALLOCATE(cn_tmp)

  END SUBROUTINE cn_shift

END MODULE pupil
