! MODULE: constants
! AUTHOR: Jouni Makitalo
! Modified: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Various constants for general use.
MODULE constants
  IMPLICIT NONE
  ! INTRINSIC SQRT

  ! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
  ! single, double and quadruple precision,
  ! by "Metcalf, Michael, et.al., 2011. Modern Fortran Explained. 4th Edition
  INTEGER, PARAMETER ::                             &
    sp = KIND(1.0),                                 &
    dp = SELECTED_REAL_KIND(2*PRECISION(1.0_sp)),   &
    qp = SELECTED_REAL_KIND(2*PRECISION(1.0_dp))

  INTEGER, PARAMETER :: r8=SELECTED_REAL_KIND(12,100) ! for time.f90

  REAL(KIND=dp), PARAMETER ::       &
    tol = 1.0E-12_dp,               &
    pi = 3.141592653589E0_dp,       &
    eps0 = 8.8541878176E-12_dp,     &
    mu0 = 4*pi*1E-7_dp,             &
    c0 = 2.997924580003E8_dp,       &
    eta0 = 3.767303134622E2_dp,     &
    hplank = 4.13566766225E-15_dp,  &
    epsilon_dp = EPSILON(pi)

  REAL(KIND=dp), PARAMETER ::       &
    rad2deg = 180.0_dp/pi,          &
    deg2rad = pi/180.0_dp

  INTEGER, PARAMETER ::             &
    prd_none = 1,                   &
    prd_2d = 2

  CHARACTER(LEN=32), PARAMETER  ::  &
    ERR0 = 'Required json node missing!',   &
    ERR1 = 'Procedure fails to run!',       &
    fmt_cmplx = '(SP, ES19.12, SP, ES19.12, "i")'
END MODULE constants
