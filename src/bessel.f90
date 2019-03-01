! MODULE: bessel
! AUTHOR: Jouni Makitalo
! Modified: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Implements the Bessel function of the first kind for real arguments.
! Modification: Use as many as possible the intrinsic functions
! or the libraries in literature.
MODULE bessel
  USE constants
  USE linalg
  USE io

  IMPLICIT NONE

    INTEGER, PARAMETER                  :: norders=50, nzeros=1
    REAL(KIND=dp), DIMENSION(0:norders,nzeros), PARAMETER :: &
      cyl_bessel_zeros = (/ &
      2.404825557695772768621631879326454643124244909146_dp, &
      3.8317059702075123156144358863081607665645452742878_dp, &
      5.1356223018406825563014016901377654569737723475006_dp, &
      6.3801618959239835062366146419427033053263036919031_dp, &
      7.588342434503804385069630007985617417369977901313_dp, &
      8.7714838159599540191228671334095605629810770148974_dp, &
      9.9361095242176848946930891269651919315561767744089_dp, &
      11.086370019245083845762764435929999140272717259475_dp, &
      12.225092264004655175612804769107398951208625498161_dp, &
      13.354300477435331066419924883491922176258557055643_dp, &
      14.47550068655454123845163765541315197630481262997_dp, &
      15.589847884455484680876296733454892607710570604222_dp, &
      16.698249933848246473203530230023595665738447769195_dp, &
      17.801435153282442168759357635981491425605478162412_dp, &
      18.899997953174023642932004510201750561812019425666_dp, &
      19.994430629816384585761082823694914662734458804394_dp, &
      21.085146113064718937846237094033074997697154586175_dp, &
      22.172494618826326205595850739307328423172829133849_dp, &
      23.256776085110037252890681563869726827092514303773_dp, &
      24.338249623407174552892989123587853235529038844182_dp, &
      25.417140814072523580431614414564773476565763475563_dp, &
      26.493647416019034039110240980824074886219498170784_dp, &
      27.567943891262234251642067737827524362018502986702_dp, &
      28.640185030763960626682088201909932183314951656221_dp, &
      29.710508889811232903428314290105243945547827848484_dp, &
      30.779039186567266165349287634595565057805936661822_dp, &
      31.845887278687318383888340006136908642719815161925_dp, &
      32.911153804984183797142849878194721929278654413294_dp, &
      33.974930058748687799845843757667220598604537471066_dp, &
      35.0372991442601950733500900167237520659009345974_dp, &
      36.098336956747724799926814797787611255237471588045_dp, &
      37.158113017536591270821315128294208223209700519944_dp, &
      38.216691189603830972763238495178457206266078377386_dp, &
      39.274130293745929565978978222762487068705776168605_dp, &
      40.330484641659080746260316354811778994977414744729_dp, &
      41.385804499172417080970530316796222134642942739139_dp, &
      42.440136490457815319593257421986207919418638265557_dp, &
      43.493523952117245974984770428241585984341041629033_dp, &
      44.546007244508866383861287730660076174581244075032_dp, &
      45.597624026432090522996531982029164361723758769649_dp, &
      46.648409498285736446144028740195805763953326276045_dp, &
      47.698396617992997613940821722256042964937880441021_dp, &
      48.747616293314518642134488939014092540648166430104_dp, &
      49.796097553616815295277680453905751293090359701853_dp, &
      50.843867703704602618032720082459070118444264224084_dp, &
      51.890952461944124593639769071063442417468368190031_dp, &
      52.937376084585930176808515550312783490565076032591_dp, &
      53.983161477928257904969383382307960375947829258943_dp, &
      55.028330299737101743494256596552294445779769739072_dp, &
      56.072903051148752617426323922727085339024827707947_dp, &
      57.116899160119174119362278697052681975274983051978_dp/)


  INTERFACE besselj
    MODULE PROCEDURE xbesselj_int, xbesselj_real, zbesselj_real
  END INTERFACE

  INTERFACE besseljd
    MODULE PROCEDURE xbesseljd_int
  END INTERFACE

CONTAINS

  REAL(KIND=dp) FUNCTION xbesselj_int(n, x)
    INTEGER, INTENT(IN)           :: n
    REAL(KIND=dp), INTENT(IN)     :: x

    xbesselj_int = BESSEL_JN(ABS(n), x) ! requires f2008
    IF ( n<0 ) THEN
      xbesselj_int = (-1)**(-n)*xbesselj_int
    END IF
  END FUNCTION xbesselj_int

  REAL(KIND=dp) FUNCTION xbesselj_real(nu, x)
    REAL(KIND=dp), INTENT(IN)     :: nu
    REAL(KIND=dp), INTENT(IN)     :: x
    COMPLEX(KIND=dp)              :: z

    z = CMPLX(x, 0.0_dp, KIND=dp)
    xbesselj_real = REAL( zbesselj_real(nu, z) )
  END FUNCTION xbesselj_real

  COMPLEX(KIND=dp) FUNCTION zbesselj_real(nu, z)
    REAL(KIND=dp), INTENT(IN)     :: nu
    COMPLEX(KIND=dp), INTENT(IN)  :: z

    REAL(KIND=dp)                 :: cyr, cyi, zr, zi, fnu
    INTEGER                       :: kode, n, nz, ierr

    cyr = 0.0_dp
    cyi = 0.0_dp
    zr = REAL(z)
    zi = AIMAG(z)
    fnu = ABS(nu) ! take the non-negative order
    kode = 1
    n = 1

    CALL ZBESJ(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
    IF(ierr/=0) THEN
      WRITE(*,*) 'ZBESJ, COMPUTATION FAILED, failed code: ierr = ', ierr
      STOP ERR1
    END IF
    ! assemble the real and imaginary parts
    zbesselj_real=CMPLX(cyr, cyi, KIND=dp)

    IF ( nu<0 ) THEN
      ! FOR NEGATIVE ORDERS, USE THE FORMULA
      ! J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
      CALL ZBESY(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
      IF(ierr/=0) THEN
        WRITE(*,*) 'ZBESY, COMPUTATION FAILED, failed code: ierr = ', ierr
        STOP ERR1
      END IF

      zbesselj_real = zbesselj_real*COS(pi*fnu) &
                    - CMPLX(cyr, cyi, KIND=dp)*SIN(pi*fnu)
    END IF

    IF ( kode == 2) THEN
      ! ON KODE=2, ZBESJ RETURNS THE SCALED FUNCTIONS
      ! CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
      !
      ! ON KODE=2, ZBESY RETURNS THE SCALED FUNCTIONS
      ! CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
      !
      ! All together, we have the scaled value as
      zbesselj_real=EXP(-ABS(zi))*zbesselj_real
    END IF
  END FUNCTION zbesselj_real

  !> Numerical evaluation (Levin's method [1]#Example 1) of Bessel transfrom,
  !! \f$\int_a^b g(x) J_n(rx) dx\f$
  !
  !> 1. D. Levin, "Fast integration of rapidly oscillatory functions,"
  !> J. Comput. Appl. Math. 67, 95–101 (1996).
  !
  !> 2. S. Xiang, "Numerical analysis of a fast integration method for 
  !> highly oscillatory functions," BIT Numer. Math. 47, 469–482 (2007).
  FUNCTION besselj_transfrom(extfun, a, b, nu, r1, r2, n) RESULT(res)
    COMPLEX(KIND=dp), EXTERNAL      :: extfun
    REAL(KIND=dp), INTENT(IN)       :: a, b, r1, r2
    !> The interval is limited to 0<=a<b<1 for focal field calculation
    INTEGER, INTENT(IN)             :: n, nu

    COMPLEX(KIND=dp)                :: res
    REAL(KIND=dp)                   :: d, h, xj, xtol
    INTEGER                         :: j
    COMPLEX(KIND=dp), DIMENSION(2)  :: pb, pa, wb, wa
    COMPLEX(KIND=dp), DIMENSION(2*n):: f, c, u, up
    COMPLEX(KIND=dp), DIMENSION(2*n, 2*n) :: lp
    COMPLEX(KIND=dp), DIMENSION(2, 2)  :: AT

    d = (a + b)/2.0_dp
    h = b - a
    xtol = 0.01_dp*h/REAL(n-1);

    ! f(1:n) = CMPLX(0.0_dp, 0.0_dp)
    f(1+n:2*n) = CMPLX(0.0_dp, 0.0_dp)
    DO j = 1, n
      xj = a +REAL(j-1, KIND=dp)*h/REAL(n-1)
      IF ( xj < xtol ) THEN
      ! avoid xj ~ 0
        xj = xj + xtol
      ELSE IF ( ABS(xj-b) < xtol ) THEN
      ! avoid xj ~ b
        xj = xj - xtol
      ELSE IF ( (ABS(xj-d) < xtol) .AND. (xj >= d) ) THEN
      ! avoid xj ~ d, see FUNCTION up_vec, u_vec
        xj = xj + xtol
      ELSE IF ( (ABS(xj-d) < xtol) .AND. (xj < d) ) THEN
        xj = xj - xtol
      END IF
      AT = TRANSPOSE(A_matrix(xj))
      ! f(n+j) = extfun(xj)
      f(j) = extfun(xj)
      lp(j,1:n) = ub1d(xj, n, a, b, 'monomial') &
                  + AT(1,1)*ub(xj, n, a, b, 'monomial')
      lp(j,n+1:2*n) = AT(1,2)*ub(xj, n, a, b, 'monomial')
      lp(j+n,1:n) = ub1d(xj, n, a, b, 'monomial') &
                    + AT(2,1)*ub(xj, n, a, b, 'monomial')
      lp(j+n,n+1:2*n) = AT(2,2)*ub(xj, n, a, b, 'monomial')
    END DO

    CALL solve_linsys(lp, f)
    c = f

    ! Least squares problem, to be done
    ! c = f
    ! CALL solve_lsqr(lp, c)
    ! CALL write_matrix('c.txt', c)
    ! STOP

    pb = pn_vec(b)
    pa = pn_vec(a)
    wb = w_vec(b)
    wa = w_vec(a)
    res = dotc(wb, pb) - dotc(wa, pa)

  CONTAINS
    FUNCTION pn_vec(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: x
      COMPLEX(KIND=dp), DIMENSION(2)  :: res

      INTEGER                         :: k

      res(1) = dotrc( ub(x, n, a, b, 'monomial'), c(1:n) )
      res(2) = dotrc( ub(x, n, a, b, 'monomial'), c(1+n:2*n) )
    END FUNCTION pn_vec

    FUNCTION w_vec(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: x
      COMPLEX(KIND=dp), DIMENSION(2)  :: res

      ! nu, r1, r2 from host association
      ! res = (/besselj(nu-1, r2*x), besselj(nu, r2*x)/)
      res = CMPLX( (/besselj(nu-1, r2*x), besselj(nu, r2*x)/), 0.0_dp) *  &
            EXP( (0,1)*r1*(1.0_dp-x**2)**(0.5_dp) )
    END FUNCTION w_vec

    FUNCTION A_matrix(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)           :: x
      COMPLEX(KIND=dp), DIMENSION(2, 2)   :: res

      COMPLEX(KIND=dp)                    :: b11, b12
      COMPLEX(KIND=dp), DIMENSION(2, 2)   :: b0, a0

      ! nu, r1, r2 from host association
      ! case1, J_nu(r_2 * x)
      ! res = RESHAPE((/(nu-1)/x, r2, -r2, -nu/x/), (/2, 2/))
      ! case2, exp(i*r1*(1-x^2)^0.5)*J_nu(r_2 * x)
      b11 = -(0,1)*r1*x*(1.0_dp-x**2)**(-0.5_dp)
      b12 = CMPLX(0.0_dp,0.0_dp)
      b0 = RESHAPE((/b11, b12, b12, b11/), (/2, 2/))
      a0 = RESHAPE((/REAL(nu-1)/x, r2, -r2, -REAL(nu)/x/), (/2, 2/))

      res = a0 + b0
    END FUNCTION A_matrix

  END FUNCTION besselj_transfrom

  SUBROUTINE besselj_transfrom_detail(extfun, a, b, nu, r1, r2, n, polyName, bt, cj)
    COMPLEX(KIND=dp), EXTERNAL      :: extfun
    REAL(KIND=dp), INTENT(IN)       :: a, b, r1, r2
    !> The interval is limited to 0<=a<b<1 for focal field calculation
    INTEGER, INTENT(IN)             :: n, nu
    CHARACTER(LEN=*), INTENT(IN)    :: polyName
    COMPLEX(KIND=dp), INTENT(OUT)   :: bt
    COMPLEX(KIND=dp), DIMENSION(2*n), INTENT(INOUT) :: cj

    REAL(KIND=dp)                   :: d, h, xj, xtol
    INTEGER                         :: j
    COMPLEX(KIND=dp), DIMENSION(2)  :: pb, pa, wb, wa
    COMPLEX(KIND=dp), DIMENSION(2*n):: f, u, up
    COMPLEX(KIND=dp), DIMENSION(2*n, 2*n) :: lp
    COMPLEX(KIND=dp), DIMENSION(2, 2)  :: AT

    ! polyName = 'monomial'
    ! polyName = 'Legendre'
    ! polyName = 'cancelA1/x'
    d = (a + b)/2.0_dp
    h = b - a
    xtol = 0.0001_dp*h/REAL(n-1);

    f(1+n:2*n) = CMPLX(0.0_dp, 0.0_dp)
    DO j = 1, n
      ! 1. equally spaced nodes
      ! xj = a +REAL(j-1, KIND=dp)*h/REAL(n-1)
      ! 2. Chebyshev nodes
      xj = d + h/2.0_dp * COS( pi * (2.0_dp*REAL(j)-1.0_dp)/(2.0_dp*REAL(n)) )
      IF ( xj < xtol ) THEN
      ! avoid xj ~ 0
        xj = xj + xtol
      ELSE IF ( ABS(xj-b) < xtol ) THEN
      ! avoid xj ~ b
        xj = xj - xtol
      ELSE IF ( (ABS(xj-d) < xtol) .AND. (xj >= d) ) THEN
      ! avoid xj ~ d, see FUNCTION up_vec, u_vec
        xj = xj + xtol
      ELSE IF ( (ABS(xj-d) < xtol) .AND. (xj < d) ) THEN
        xj = xj - xtol
      END IF
      AT = TRANSPOSE(A_matrix(xj))
      ! f(n+j) = extfun(xj)
      f(j) = extfun(xj)
      lp(j,1:n) = ub1d(xj, n, a, b,polyName) &
                  + AT(1,1)*ub(xj, n, a, b, polyName)
      lp(j,n+1:2*n) = AT(1,2)*ub(xj, n, a, b, polyName)
      lp(j+n,1:n) = AT(2,1)*ub(xj, n, a, b, polyName)
      lp(j+n,n+1:2*n) = ub1d(xj, n, a, b, polyName) &
                    + AT(2,2)*ub(xj, n, a, b, polyName)
    END DO

    ! WRITE(*,*) nu, a, b, r1, r2 !debug
    CALL solve_linsys(lp, f)
    cj = f
    ! CALL write_matrix('c.txt', cj)! debug
    ! STOP! debug

    pb = pn_vec(b)
    pa = pn_vec(a)
    wb = w_vec(b)
    wa = w_vec(a)
    bt = dotc(wb, pb) - dotc(wa, pa)
  CONTAINS
    FUNCTION pn_vec(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: x
      COMPLEX(KIND=dp), DIMENSION(2)  :: res

      ! WRITE(*,*) x, a, b
      res(1) = dotrc( ub(x, n, a, b, polyName), cj(1:n) )
      res(2) = dotrc( ub(x, n, a, b, polyName), cj(1+n:2*n) )
    END FUNCTION pn_vec

    FUNCTION w_vec(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: x
      COMPLEX(KIND=dp), DIMENSION(2)  :: res

      res = CMPLX( (/besselj(nu-1, r2*x), besselj(nu, r2*x)/), 0.0_dp) *  &
            EXP( (0,1)*r1*(1.0_dp-x**2)**(0.5_dp) )
    END FUNCTION w_vec

    FUNCTION A_matrix(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)           :: x
      COMPLEX(KIND=dp), DIMENSION(2, 2)   :: res

      COMPLEX(KIND=dp)                    :: b11, b12
      COMPLEX(KIND=dp), DIMENSION(2, 2)   :: b0, a0

      b11 = -(0,1)*r1*x*(1.0_dp-x**2)**(-0.5_dp)
      b12 = CMPLX(0.0_dp,0.0_dp)
      b0 = RESHAPE((/b11, b12, b12, b11/), (/2, 2/))
      a0 = RESHAPE((/REAL(nu-1)/x, r2, -r2, -REAL(nu)/x/), (/2, 2/))
      res = a0 + b0
    END FUNCTION A_matrix
  END SUBROUTINE besselj_transfrom_detail

  ! directives of bessel function of the first kind
  ! with integer order, real argument
  RECURSIVE FUNCTION xbesseljd_int(n, x, d) RESULT(ybesselj)
    INTEGER, INTENT(IN)           :: n, d
    REAL(KIND=dp), INTENT(IN)     :: x
    REAL(KIND=dp)                 :: ybesselj

    IF ( d==1 ) THEN
      ybesselj = besselj(n-1,x) - n/x*besselj(n, x)
    ELSE IF ( d==2 ) THEN
      ybesselj = xbesseljd_int(n-1, x, 1) + n/(x**2)*besselj(n, x) &
               - n/x*xbesseljd_int(n, x, 1)
    END IF
  END FUNCTION xbesseljd_int

END MODULE bessel
