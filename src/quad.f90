! MODULE: quad
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Quadrature routines for integartion over triangles and tetrahedra (Gauss-Legendre).
! Also routines for Simpson quadrature over interval and rectangle.
MODULE quad
  USE data
  USE mesh

  IMPLICIT NONE

  REAL (KIND=dp), DIMENSION(4), PARAMETER :: qw4 = (/-0.56250000_dp,&
       0.52083333_dp,0.52083333_dp,0.52083333_dp/)

  REAL (KIND=dp), DIMENSION(4,3), PARAMETER :: quadFactors4 = RESHAPE((/&
       1.0_dp/3.0_dp,&
       0.6_dp,&
       0.2_dp,&
       0.2_dp,&

       1.0_dp/3.0_dp,&
       0.2_dp,&
       0.6_dp,&
       0.2_dp,&

       1.0_dp/3.0_dp,&
       0.2_dp,&
       0.2_dp,&
       0.6_dp/), (/4,3/))

  REAL (KIND=dp), DIMENSION(7), PARAMETER :: qw7 = (/0.22500000_dp,&
       0.13239415_dp, 0.13239415_dp, 0.13239415_dp, 0.12593918_dp,&
       0.12593918_dp, 0.12593918_dp/)

  REAL (KIND=dp), DIMENSION(7,3), PARAMETER :: quadFactors7 = RESHAPE((/&
       0.33333333_dp,&
       0.05971587_dp,&
       0.47014206_dp,&
       0.47014206_dp,&
       0.79742698_dp,&
       0.10128650_dp,&
       0.10128650_dp,&

       0.33333333_dp,&
       0.47014206_dp,&
       0.05971587_dp,&
       0.47014206_dp,&
       0.10128650_dp,&
       0.79742698_dp,&
       0.10128650_dp,&

       0.33333333_dp,&
       0.47014206_dp,&
       0.47014206_dp,&
       0.05971587_dp,&
       0.10128650_dp,&
       0.10128650_dp,&
       0.79742698_dp&
       /), (/7,3/))

  REAL (KIND=dp), DIMENSION(6), PARAMETER :: qw6 = (/&
       0.109951743655322_dp,&
       0.109951743655322_dp,&
       0.109951743655322_dp,&
       0.223381589678011_dp,&
       0.223381589678011_dp,&
       0.223381589678011_dp/)

  REAL (KIND=dp), DIMENSION(6,3), PARAMETER :: quadFactors6 = RESHAPE((/&
       0.816847572980459_dp,&
       0.091576213509771_dp,&
       0.091576213509771_dp,&
       0.108103018168070_dp,&
       0.445948490915965_dp,&
       0.445948490915965_dp,&

       0.091576213509771_dp,&
       0.816847572980459_dp,&
       0.091576213509771_dp,&
       0.445948490915965_dp,&
       0.108103018168070_dp,&
       0.445948490915965_dp,&

       0.091576213509771_dp,&
       0.091576213509771_dp,&
       0.816847572980459_dp,&
       0.445948490915965_dp,&
       0.445948490915965_dp,&
       0.108103018168070_dp/),(/6,3/))

  REAL (KIND=dp), DIMENSION(1), PARAMETER :: qw1 = 1.0_dp

  REAL (KIND=dp), DIMENSION(1,3), PARAMETER :: quadFactors1 = RESHAPE((/&
       0.333333333333333_dp,&
       0.333333333333333_dp,&
       0.333333333333333_dp/),(/1,3/))

  REAL (KIND=dp), DIMENSION(3), PARAMETER :: qw3 = (/&
       0.333333333333333_dp,&
       0.333333333333333_dp,&
       0.333333333333333_dp/)

  REAL (KIND=dp), DIMENSION(3,3), PARAMETER :: quadFactors3 = RESHAPE((/&
       0.666666666666667,&
       0.166666666666667,&
       0.166666666666667,&

       0.166666666666667,&
       0.666666666666667,&
       0.166666666666667,&

       0.166666666666667,&
       0.166666666666667,&
       0.666666666666667/),(/3,3/))

  REAL (KIND=dp), DIMENSION(13), PARAMETER :: qw13 = (/&
       -0.149570044467670,&
       0.175615257433204,&
       0.175615257433204,&
       0.175615257433204,&
       0.053347235608839,&
       0.053347235608839,&
       0.053347235608839,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257,&
       0.077113760890257/)

  REAL (KIND=dp), DIMENSION(13,3), PARAMETER :: quadFactors13 = RESHAPE((/&
       0.333333333333333,&
       0.479308067841923,&
       0.260345966079038,&
       0.260345966079038,&
       0.869739794195568,&
       0.065130102902216,&
       0.065130102902216,&
       0.638444188569809,&
       0.638444188569809,&
       0.312865496004875,&
       0.312865496004875,&
       0.048690315425316,&
       0.048690315425316,&

       0.333333333333333,&
       0.260345966079038,&
       0.260345966079038,&
       0.479308067841923,&
       0.065130102902216,&
       0.065130102902216,&
       0.869739794195568,&
       0.048690315425316,&
       0.312865496004875,&
       0.048690315425316,&
       0.638444188569809,&
       0.312865496004875,&
       0.638444188569809,&

       0.333333333333333,&
       0.260345966079038,&
       0.479308067841923,&
       0.260345966079038,&
       0.065130102902216,&
       0.869739794195568,&
       0.065130102902216,&
       0.312865496004875,&
       0.048690315425316,&
       0.638444188569809,&
       0.048690315425316,&
       0.638444188569809,&
       0.312865496004875/),(/13,3/))

  REAL (KIND=dp), DIMENSION(1), PARAMETER :: volQw1 = (/1.0_dp/)

  REAL (KIND=dp), DIMENSION(1,4), PARAMETER :: volQuadFactors1 = RESHAPE((/&
       0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp/), (/1,4/))

  REAL (KIND=dp), DIMENSION(4), PARAMETER :: volQw4 = (/&
       0.25_dp,&
       0.25_dp,&
       0.25_dp,&
       0.25_dp/)

  REAL (KIND=dp), DIMENSION(4,4), PARAMETER :: volQuadFactors4 = RESHAPE((/&
       0.585410196624969,&
       0.138196601125011,&
       0.138196601125011,&
       0.138196601125011,&

       0.138196601125011,&
       0.585410196624969,&
       0.138196601125011,&
       0.138196601125011,&

       0.138196601125011,&
       0.138196601125011,&
       0.585410196624969,&
       0.138196601125011,&

       0.138196601125011,&
       0.138196601125011,&
       0.138196601125011,&
       0.585410196624969/),(/4,4/))

  REAL (KIND=dp), DIMENSION(11), PARAMETER :: volQw11 = (/&
       -0.013155555555556,&

       0.007622222222222,&
       0.007622222222222,&
       0.007622222222222,&
       0.007622222222222,&

       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889,&
       0.024888888888889/)

  REAL (KIND=dp), DIMENSION(11,4), PARAMETER :: volQuadFactors11 = RESHAPE((/&
       0.250000000000000,&
       0.785714285714286,&
       0.071428571428571,&
       0.071428571428571,&
       0.071428571428571,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&

       0.250000000000000,&
       0.071428571428571,&
       0.785714285714286,&
       0.071428571428571,&
       0.071428571428571,&
       0.399403576166799,&
       0.100596423833201,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&

       0.250000000000000,&
       0.071428571428571,&
       0.071428571428571,&
       0.785714285714286,&
       0.071428571428571,&
       0.100596423833201,&
       0.399403576166799,&
       0.399403576166799,&
       0.100596423833201,&
       0.100596423833201,&
       0.399403576166799,&

       0.250000000000000,&
       0.071428571428571,&
       0.071428571428571,&
       0.071428571428571,&
       0.785714285714286,&
       0.100596423833201,&
       0.399403576166799,&
       0.100596423833201,&
       0.399403576166799,&
       0.399403576166799,&
       0.100596423833201/),(/11,4/))

CONTAINS

  SUBROUTINE prepare_quad( domain )
    TYPE(domain_type), DIMENSION(-1:), INTENT(INOUT) :: domain

    INTEGER   :: idom, iquad

    domain(-1)%quad(2)%qd = tri_quad_data( TRIM(domain(-1)%quad(2)%rule) )
    domain(-1)%quad(3)%qd = tetra_quad_data( TRIM(domain(-1)%quad(3)%rule) )
    DO idom = 0, (SIZE(domain(0:))-1)
      DO iquad = 1, SIZE(domain(idom)%quad)
        IF ( domain(idom)%quad(iquad)%dim == 2 ) THEN
          domain(idom)%quad(iquad)%qd = &
            tri_quad_data( TRIM(domain(idom)%quad(iquad)%rule) )
        ELSE IF ( domain(idom)%quad(iquad)%dim == 3 ) THEN
          domain(idom)%quad(iquad)%qd = &
            tetra_quad_data( TRIM(domain(idom)%quad(iquad)%rule) )
        END IF
      END DO
    END DO
  END SUBROUTINE prepare_quad

  FUNCTION tri_quad_data(name) RESULT(qd)
    CHARACTER (LEN=*) :: name
    TYPE(quad_data) :: qd

    IF(name=='tri_gl1') THEN
       qd%description = '1-point Gauss-Legendre'
       qd%num_nodes = 1
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,3))
       qd%weights = qw1
       qd%nodes = quadFactors1
    ELSE IF(name=='tri_gl3') THEN
       qd%description = '3-point Gauss-Legendre'
       qd%num_nodes = 3
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,3))
       qd%weights = qw3
       qd%nodes = quadFactors3
    ELSE IF(name=='tri_gl4') THEN
       qd%description = '4-point Gauss-Legendre'
       qd%num_nodes = 4
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,3))
       qd%weights = qw4
       qd%nodes = quadFactors4
    ELSE IF(name=='tri_gl7') THEN
       qd%description = '7-point Gauss-Legendre'
       qd%num_nodes = 7
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,3))
       qd%weights = qw7
       qd%nodes = quadFactors7
    ELSE IF(name=='tri_gl13') THEN
       qd%description = '13-point Gauss-Legendre'
       qd%num_nodes = 13
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,3))
       qd%weights = qw13
       qd%nodes = quadFactors13
    ELSE
       WRITE(*,*) 'Unrecognized quadrature type ', name, '!'
       STOP
    END IF
  END FUNCTION tri_quad_data

  FUNCTION tetra_quad_data(name) RESULT(qd)
    CHARACTER (LEN=*) :: name
    TYPE(quad_data) :: qd

    IF(name=='tetra_gl1') THEN
       qd%description = '1-point Gauss-Legendre'
       qd%num_nodes = 1
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,4))
       qd%weights = volQw1
       qd%nodes = volQuadFactors1
    ELSE IF(name=='tetra_gl4') THEN
       qd%description = '4-point Gauss-Legendre'
       qd%num_nodes = 4
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,4))
       qd%weights = volQw4
       qd%nodes = volQuadFactors4
    ELSE IF(name=='tetra_gl11') THEN
       qd%description = '11-point Gauss-Legendre'
       qd%num_nodes = 11
       ALLOCATE(qd%weights(qd%num_nodes), qd%nodes(qd%num_nodes,4))
       qd%weights = volQw11
       qd%nodes = volQuadFactors11
    ELSE
       WRITE(*,*) 'Unrecognized quadrature type ', name, '!'
       STOP
    END IF
  END FUNCTION tetra_quad_data

  SUBROUTINE delete_quad_data(qd)
    TYPE(quad_data), INTENT(INOUT) :: qd

    DEALLOCATE(qd%weights)
    DEALLOCATE(qd%nodes)
  END SUBROUTINE delete_quad_data

  FUNCTION quad_tri_points(qd, faceind, mesh) RESULT(res)
    TYPE(quad_data), INTENT(IN) :: qd
    INTEGER, INTENT(IN) :: faceind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: res
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3
    INTEGER :: n

    p1 = mesh%nodes(mesh%faces(faceind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%faces(faceind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%faces(faceind)%node_indices(3))%p

    DO n=1,qd%num_nodes
       res(:,n) = qd%nodes(n,3)*p1 + qd%nodes(n,1)*p2 + qd%nodes(n,2)*p3
    END DO
  END FUNCTION quad_tri_points

  FUNCTION quad_tetra_points(qd, solidind, mesh) RESULT(res)
    TYPE(quad_data), INTENT(IN) :: qd
    INTEGER, INTENT(IN) :: solidind
    TYPE(mesh_container), INTENT(IN) :: mesh
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: res
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, p3, p4
    INTEGER :: n

    p1 = mesh%nodes(mesh%solids(solidind)%node_indices(1))%p
    p2 = mesh%nodes(mesh%solids(solidind)%node_indices(2))%p
    p3 = mesh%nodes(mesh%solids(solidind)%node_indices(3))%p
    p4 = mesh%nodes(mesh%solids(solidind)%node_indices(4))%p

    DO n=1,qd%num_nodes
       res(:,n) = qd%nodes(n,1)*p1 + qd%nodes(n,2)*p2 + qd%nodes(n,3)*p3 + qd%nodes(n,4)*p4
    END DO
  END FUNCTION quad_tetra_points

  ! n is the number of subintervals -> number of nodes is n+1.
  ! n must be even.
  SUBROUTINE get_simpsons_weights(a, b, n, weights)
    REAL (KIND=dp), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp), DIMENSION(n+1), INTENT(OUT) :: weights
    REAL (KIND=dp) :: h

    IF(MOD(n,2)/=0) THEN
       WRITE(*,*) "Number of subintervals for Simpson's rule must be even!"
       STOP
    END IF

    h = (b-a)/n

    weights((/1,n+1/)) = h/3.0_dp
    weights(3:(n-1):2) = 2.0_dp*h/3.0_dp
    weights(2:n:2) = 4.0_dp*h/3.0_dp
  END SUBROUTINE get_simpsons_weights

  ! n is the number of subintervals -> number of nodes is n+1.
  ! n must be even.
  SUBROUTINE get_simpsons_points(a, b, n, points)
    REAL (KIND=dp), INTENT(IN) :: a, b
    INTEGER, INTENT(IN) :: n
    REAL (KIND=dp), DIMENSION(n+1), INTENT(OUT) :: points
    REAL (KIND=dp) :: h
    INTEGER :: i

    IF(MOD(n,2)/=0) THEN
       WRITE(*,*) "Number of subintervals for Simpson's rule must be even!"
       STOP
    END IF

    h = (b-a)/n

    points(1:(n+1)) = a + h*(/(i,i=0,n)/)
  END SUBROUTINE get_simpsons_points

  RECURSIVE FUNCTION asqz_aux(f, a, b, eps, s, fa, fb, fm, level, maxDepth) RESULT(res)
    COMPLEX (KIND=dp), EXTERNAL :: f
    REAL (KIND=dp), INTENT(IN) :: a, b, eps
    COMPLEX (KIND=dp), INTENT(IN) :: s, fa, fb, fm
    INTEGER, INTENT(IN) :: level, maxDepth

    COMPLEX (KIND=dp) :: res, flm, frm, sleft, sright, s2
    REAL(KIND=dp) :: m, h, lm, rm

    m = (a + b)/2
    h = b - a

    lm = (a + m)/2
    rm = (m + b)/2
    
    IF ( (eps/2.0_dp==eps) .OR. (a==lm) ) THEN
      WRITE(*,*) 'Serious numerical trouble: it won''','t converge in asqz.'
      STOP
    ELSE
      ! continue calculation
      flm = f(lm)
      frm = f(rm)

      sleft = (h/12)*(fa + 4*flm + fm)
      sright = (h/12)*(fm + 4*frm + fb)

      s2 = sleft + sright

      IF(level>=maxDepth .OR. (ABS(s2 - s)<=15*eps .AND. level>1) ) THEN
        res = s2 + (s2 - s)/15
        IF ( level>=maxDepth ) THEN
          WRITE(*,*) 'Numerical integration fails as the setup accuracy is not'// &
                      'achieved for maximum (',maxDepth,') steps of iteration.'
          STOP
        END IF
      ELSE
         res = asqz_aux(f, a, m, eps/2, sleft, fa, fm, flm, level+1, maxDepth) +&
              asqz_aux(f, m, b, eps/2, sright, fm, fb, frm, level+1, maxDepth)
      END IF
    END IF
  END FUNCTION asqz_aux

  ! Adaptive Simpson's method based on listing at
  ! http://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method#C
  FUNCTION asqz(f, a, b, eps, maxDepth) RESULT(res)
    COMPLEX(KIND=dp), EXTERNAL    :: f
    REAL(KIND=dp), INTENT(IN)     :: a, b, eps
    INTEGER, INTENT(IN)           :: maxDepth

    COMPLEX (KIND=dp) :: fa, fb, fm, s, res
    REAL (KIND=dp) :: m, h

    m = (a + b)/2
    h = b - a

    IF (h <= eps) THEN
      WRITE(*,*) 'The interval is too narrow for numerical integration!'
      STOP
    END IF

    fa = f(a)
    fb = f(b)
    fm = f(m)

    s = (h/6)*(fa + 4*fm + fb)

    res = asqz_aux(f, a, b, eps, s, fa, fb, fm, 0, maxDepth)
  END FUNCTION asqz

  ! Integrates f(x,y) over [x1,x2]x[y1,y2].
  SUBROUTINE asqz2(fInteg, x1, x2, y1, y2, eps, maxDepth, res)
    COMPLEX (KIND=dp), EXTERNAL :: fInteg
    REAL (KIND=dp), INTENT(IN) :: x1, x2, y1, y2, eps
    INTEGER, INTENT(IN) :: maxDepth
    COMPLEX (KIND=dp), INTENT(INOUT) :: res

    ! Auxiliary variable to make f(x,y) to appear a
    ! function of single variable for two nested
    ! integration routines.
    REAL (KIND=dp) :: gy

    res = asqz(fnested, y1, y2, eps, maxDepth)

  CONTAINS
    ! Evaluates f with y fixed to global value gy.
    FUNCTION fproxy(x) RESULT(z)
      REAL (KIND=dp), INTENT(IN) :: x
      COMPLEX (KIND=dp) :: z

      z = fInteg(x,gy)
      ! Here, 'f(x,gy)' is referring to 'integ(theta, phi)'
    END FUNCTION fproxy

    FUNCTION fnested(y) RESULT(z)
      REAL (KIND=dp), INTENT(IN) :: y
      COMPLEX (KIND=dp) :: z

      gy = y

      z = asqz(fproxy, x1, x2, eps, maxDepth)
    END FUNCTION fnested

  END SUBROUTINE asqz2
END MODULE quad
