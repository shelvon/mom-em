! MODULE: aux
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Auxiliary routines for various purposes, such as string manipulation, file access and simple
! mathematical tools. This module is a good candidate for refactorization in the future.
MODULE aux
  USE constants
  USE data

  IMPLICIT NONE

CONTAINS

  FUNCTION num2str(num) RESULT(str)
    INTEGER, INTENT(IN)   :: num
    CHARACTER(LEN=15)     :: str

    WRITE(str, *) num
    str = TRIM(ADJUSTL(str))
  END FUNCTION num2str

  FUNCTION str2num(str) RESULT(num)
    CHARACTER(LEN=*), INTENT(IN)  :: str
    INTEGER                       :: num
    CHARACTER(LEN=15)             :: str_trim

    str_trim = TRIM(ADJUSTL(str))
    READ(str_trim, *) num
  END FUNCTION str2num
  ! Finds nval smallest-in-magnitude values in val, whose size is dim.
  ! Returns the linear indices to these values.
  FUNCTION find_smallest(val, dim, nval) RESULT(ind)
    COMPLEX (KIND=dp), DIMENSION(:), INTENT(IN) :: val
    INTEGER, INTENT(IN) :: dim, nval

    INTEGER, DIMENSION(nval) :: ind
    COMPLEX (KIND=dp), DIMENSION(dim) :: tmpval
    INTEGER :: n
    INTEGER, DIMENSION(1) :: loc
    REAL (KIND=dp) :: mval

    tmpval(:) = val(:)

    mval = MAXVAL(ABS(val))

    DO n=1,nval
       loc = MINLOC(ABS(tmpval))
       tmpval(loc(1)) = mval

       ind(n) = loc(1)
    END DO
  END FUNCTION find_smallest

  ! Returns the three-letter filename extension of the filename.
  FUNCTION getext(filename) RESULT(ext)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    CHARACTER (LEN=3) :: ext
    INTEGER :: n

    n = LEN_TRIM(filename)

    ext = filename((n-2):n)
  END FUNCTION getext

  ! Replaces the three-letter filename extension in filename with ext.
  SUBROUTINE replace_ext(filename, ext)
    CHARACTER (LEN=*), INTENT(INOUT) :: filename
    CHARACTER (LEN=3), INTENT(IN) :: ext
    INTEGER :: n

    n = LEN_TRIM(filename)
    filename((n-2):n) = ext
  END SUBROUTINE replace_ext

  FUNCTION factorial_n(n) RESULT(res)
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp) :: res
    INTEGER :: i

    res = 1.0_dp

    DO i=2,n
       res = res*REAL(i,KIND=dp)
    END DO
  END FUNCTION factorial_n

  ! Linear interpolation between start and end with t in [0,1].
  FUNCTION linterp(start, end, t) RESULT(res)
    REAL (KIND=dp), INTENT(IN) :: start, end, t
    REAL (KIND=dp) :: res

    res = start + t*(end - start)
  END FUNCTION linterp

  FUNCTION linspace(start, end, n) RESULT(res)
    REAL(KIND=dp), INTENT(IN) :: start, end
    INTEGER, INTENT(IN)       :: n
    REAL(KIND=dp)             :: res(1:n)

    INTEGER                   :: i

    IF ( n .LT. 1 ) THEN
      WRITE(*,*) 'n must be larger than 1.'
    ELSE IF ( n .EQ. 1 ) THEN
      res = start
    ELSE
      DO i = 1, n
        res(i) = linterp( start, end, REAL(i-1, KIND=dp)/REAL(n-1, KIND=dp) )
      END DO
    END IF
  END FUNCTION linspace

  FUNCTION zlinspace(za, zb, nx, ny) RESULT(res)
    COMPLEX(KIND=dp), INTENT(IN)  :: za, zb
    INTEGER, INTENT(IN)           :: nx, ny
    COMPLEX(KIND=dp)              :: res(1:nx*ny)

    REAL(KIND=dp)                 :: xa, ya, xb, yb, x, y
    INTEGER                       :: i, ix, iy

    xa = REAL(za)
    ya = AIMAG(za)
    xb = REAL(zb)
    yb = AIMAG(zb)
    IF ( (nx .LT. 1) .OR. (ny .LT. 1) ) THEN
      WRITE(*,*) 'n must be larger than or equal to 1.'
    ELSE
      i = 0
      DO ix = 1, nx
        IF ( nx .EQ. 1 ) THEN
          x = xa
        ELSE
          x = xa + REAL(ix-1, KIND=dp)*(xb-xa)/REAL(nx-1, KIND=dp)
        END IF
        DO iy = 1, ny
          IF ( ny .EQ. 1 ) THEN
            y = ya
          ELSE
            y = ya + REAL(iy-1, KIND=dp)*(yb-ya)/REAL(ny-1, KIND=dp)
          END IF
          i = i+1
          res(i) = CMPLX(x, y, KIND=dp)
        END DO
      END DO
    END IF
  END FUNCTION zlinspace

  ! Generate 3d matrix (i.e. gridx, gridy, gridz) from 1d vector x, y, z.
  ! Dimensions are in the order of y, x, z.
  SUBROUTINE meshgrid(x, y, z, gridx, gridy, gridz)
    REAL (KIND=dp), DIMENSION(:), INTENT(IN)                    :: x, y, z
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT)  :: gridx
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT)  :: gridy
    REAL (KIND=dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(OUT)  :: gridz

    INTEGER                                                     :: nx, ny, nz

    nx = SIZE(x)
    ny = SIZE(y)
    nz = SIZE(z)
    ALLOCATE(gridx(1:ny, 1:nx, 1:nz))
    ALLOCATE(gridy(1:ny, 1:nx, 1:nz))
    ALLOCATE(gridz(1:ny, 1:nx, 1:nz))

    gridx=(SPREAD(SPREAD(x,1,ny),3,nz))
    gridy=(SPREAD(SPREAD(y,2,nx),3,nz))
    gridz=(SPREAD(SPREAD(z,1,nx),1,ny))
  END SUBROUTINE meshgrid

  ! Reads a matrix from file.
  SUBROUTINE get_matrix(filename, mat, nrows, ncols)
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i
    INTEGER, PARAMETER :: nattempts = 10
    REAL (KIND=dp), DIMENSION(nrows,ncols), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nrows, ncols

    ! Attempt to open the file for nattempts times before resorting
    ! to failure. This is so that multiple instances of the program
    ! can access the same refractive index files simultaneously.

    DO i=1,nattempts
       OPEN(fid, FILE=TRIM(filename), ACTION='READ', IOSTAT=iovar)
       IF(iovar>0) THEN
          CALL SLEEP(1)
       END IF
    END DO

    IF(i==nattempts .AND. iovar>0) THEN
       WRITE(*,*) 'Could not open matrix data file!'
       STOP
    END IF

    iovar = 0
    DO i=1,nrows
       READ(fid,*) mat(i,1:ncols)
    END DO

    CLOSE(fid)

  END SUBROUTINE get_matrix

  ! Finds from list the value that is closest to val and returns index to it.
  FUNCTION find_closest(val, list) RESULT(ind)
    REAL (KIND=dp), INTENT(IN) :: val
    REAL (KIND=dp), DIMENSION(:), INTENT(IN) :: list
    INTEGER :: ind, n
    REAL (KIND=dp) :: dist

    ind = 1
    dist = ABS(val-list(1))

    IF(SIZE(list)==1) THEN
       RETURN
    END IF

    DO n=2,SIZE(list)
       IF(ABS(val-list(n))<dist) THEN
          ind = n
          dist = ABS(val-list(n))
       END IF
    END DO
  END FUNCTION find_closest

  ! Writes matrix data to file.
  SUBROUTINE write_data(filename, data)
    REAL (KIND=dp), DIMENSION(:,:), INTENT(IN) :: data
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i, j

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1,SIZE(data,1)
       DO j=1,SIZE(data,2)
          WRITE(fid, '(EN15.3)', ADVANCE='NO') data(i,j)
       END DO
       WRITE(fid, '(/)')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_data

  ! Writes complex matrix data to file. Complex numbers are written in the form x+iy.
  SUBROUTINE write_cdata(filename, data)
    COMPLEX (KIND=dp), DIMENSION(:,:), INTENT(IN) :: data
    CHARACTER (LEN=*), INTENT(IN) :: filename
    INTEGER :: fid = 10, iovar, i, j

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1,SIZE(data,1)
       DO j=1,SIZE(data,2)
          WRITE(fid, '(EN15.3,"+",EN15.3,"i")', ADVANCE='NO') REAL(data(i,j)), IMAG(data(i,j))
       END DO
       WRITE(fid, '(/)')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_cdata

  ! Returns the direction vector for spherical polar angles theta and phi.
  FUNCTION get_dir(theta, phi) RESULT(dir)
    REAL (KIND=dp), INTENT(IN) :: theta, phi
    REAL (KIND=dp), DIMENSION(3) :: dir
    
    dir = (/SIN(theta)*COS(phi), SIN(theta)*SIN(phi), COS(theta)/)
  END FUNCTION get_dir

  ! Returns a 'polarization' vector of plane-wave propagating in direction given
  ! by get_dir. The angle psi is the rotation of the polarization along the direction of propagation.
  FUNCTION get_pol(theta, phi, psi) RESULT(pol)
    REAL (KIND=dp), INTENT(IN) :: theta, phi, psi
    REAL (KIND=dp), DIMENSION(3) :: pol

    pol = (/COS(psi)*SIN(phi) - SIN(psi)*COS(theta)*COS(phi),&
         -COS(psi)*COS(phi) - SIN(psi)*COS(theta)*SIN(phi),&
         SIN(psi)*SIN(theta)/)
  END FUNCTION get_pol

  ! 1st derivative of basis function defined on interval (a,b)
  FUNCTION ub1d(x, n, a, b, polyName) RESULT(res)
    REAL(KIND=dp), INTENT(IN)       :: x
    INTEGER, INTENT(IN)             :: n
    REAL(KIND=dp), INTENT(IN)       :: a
    REAL(KIND=dp), INTENT(IN)       :: b
    CHARACTER(LEN=*), INTENT(IN)    :: polyName
    REAL(KIND=dp), DIMENSION(n)     :: res

    INTEGER                         :: k
    REAL(KIND=dp)                   :: d, h, t

    d = (a+b)/2.0_dp
    h = b-a
    t = 2.0_dp*(x-d)/h
    DO k = 1, n
      IF ( TRIM(polyName) == 'monomial' ) THEN
        ! 1. Monomial polynomials
        res(k) = REAL(k-1)*(x-d)**(k-2)
      ELSE IF ( TRIM(polyName) == 'Legendre' ) THEN
        ! 2. Legendre polynomials
        IF ( k == 1 ) THEN
          res(k) = 0.0_dp
        ELSE
          res(k) = ( t*polyLegendre(k-1, t)  -polyLegendre(k-2, t) )  &
                   * REAL(k-1)/(t**2-1) * h/2.0_dp
        END IF
      ELSE IF ( TRIM(polyName) == 'cancelA1/x' ) THEN
        ! 3. Polynomials to cancel 1/x in A
        res(k) = REAL(k)*(x)**(k-1)
      ELSE IF ( TRIM(polyName) == 'test1' ) THEN
        ! 4. test polynomial
        res(k) = (x-d)**(k-1)+x*(k-1)*(x-d)**(k-2)
      ELSE IF ( TRIM(polyName) == 'test2' ) THEN
        ! 5. test polynomial
        res(k) = REAL(k)*(x)**(k-1)*SQRT(1-x**2) - x**(k+1)/SQRT(1-x**2)
      END IF
    END DO
  END FUNCTION ub1d

  ! basis function defined on interval (a,b)
  FUNCTION ub(x, n, a, b, polyName) RESULT(res)
    REAL(KIND=dp), INTENT(IN)       :: x
    INTEGER, INTENT(IN)             :: n
    REAL(KIND=dp), INTENT(IN)       :: a
    REAL(KIND=dp), INTENT(IN)       :: b
    CHARACTER(LEN=*), INTENT(IN)    :: polyName
    REAL(KIND=dp), DIMENSION(n)     :: res

    INTEGER                         :: k
    REAL(KIND=dp)                   :: d, h, t

    d = (a+b)/2.0_dp
    h = b-a
    t = 2.0_dp*(x-d)/h
    DO k = 1, n
      IF ( TRIM(polyName) == 'monomial' ) THEN
        ! 1. Monomial polynomials
        res(k) = (x-d)**(k-1)
      ELSE IF ( TRIM(polyName) == 'Legendre' ) THEN
        ! 2. Legendre polynomials
        res(k) = polyLegendre(k-1, t)
      ELSE IF ( TRIM(polyName) == 'cancelA1/x' ) THEN
        ! 3. Polynomials to cancel 1/x in A
        res(k) = (x)**(k)
      ELSE IF ( TRIM(polyName) == 'test1' ) THEN
        ! 4. test polynomial
        res(k) = x*(x-d)**(k-1)
      ELSE IF ( TRIM(polyName) == 'test2' ) THEN
        ! 5. test polynomial
        res(k) = x*SQRT(1-x**2)*(x)**(k-1)
      END IF
    END DO
  END FUNCTION ub

  FUNCTION polyLegendre(n, x) RESULT(res)
    INTEGER, INTENT(IN)       :: n
    REAL(KIND=dp), INTENT(IN) :: x
    REAL(KIND=dp)             :: res

    INTEGER :: k

    res = 0.0_dp
    IF ( n>=0 ) THEN
      DO k = 0, n
        res = res + factorial_n(n)/(factorial_n(k)*factorial_n(n-k))  &
                    * (x-1)**(n-k) * (x+1)**k
      END DO
      res = res/(2.0_dp**n)
    ELSE
      WRITE(*,*) 'n must not be negative in polyLegendre(n,x)!'
      STOP
    END IF

  END FUNCTION polyLegendre

!  FUNCTION PolyLagrange(n, x) RESULT(res)
!    INTEGER, INTENT(IN)       :: n
!    REAL(KIND=dp), INTENT(IN) :: x
!    REAL(KIND=dp)             :: res
!
!    INTEGER :: k
!
!    res = 0.0_dp
!    DO k = 0, n
!      res = res + factorial_n(n)/(factorial_n(k)*factorial_n(n-k))  &
!                  * (x-1)**(n-k) * (x+1)**k
!    END DO
!    res = res/(2.0_dp**n)
!
!  END FUNCTION PolyLagrange
END MODULE aux
