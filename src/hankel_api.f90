! MODULE: hankel_api
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! This module is a caller to the driver of toms794/hankel_transform.f
! THis is a wrapper for the Hankel transform in toms794/src.f
! which is originally written by THOMAS WIEDER
! Refer to toms794/*.f for more details.
MODULE hankel_api
  USE constants
  IMPLICIT NONE

  INTERFACE
    ! subroutine in toms794/hankel_transform.f
    SUBROUTINE HANKEL_TRANSFORM(HANKEL_KERNEL,gNU,gXI,gXLIMIT,HT,IERR)
      DOUBLE PRECISION              :: HANKEL_KERNEL
      EXTERNAL                      :: HANKEL_KERNEL
      DOUBLE PRECISION, INTENT(IN)  :: gNU, gXI, gXLIMIT
      DOUBLE PRECISION, INTENT(OUT) :: HT
      INTEGER, INTENT(OUT)          :: IERR
!C     MOST IMPORTANT VARIABLES:
!C
!C     XI = POSITION AT WHICH THE HANKEL TRANSFORM *** HT ***
!C          OF FUNCTION *** FUN *** HAS TO BE EVALUATED.
!C     NU = ORDER OF HANKEL TRANSFORM
!C        = ORDER OF BESSEL FUNCTION *** DHJNUX ***.
!C
!C                     ( infty
!C     HT(FUN,NU,XI) =  |        (XI * X)**(1/2) JNUX(XI * X) FUN(X) dX
!C                   0 )
!C
!C     HT = HANKEL TRANSFORM OF FUNCTION *** FUN ***.
!C
!C     THE FUNCTION *** FUN *** IS OF THE FORM: Y = FUN(X)
!C     IERR: ERROR FLAG
    END SUBROUTINE HANKEL_TRANSFORM
  END INTERFACE

CONTAINS
  ! The zhankel is a Hankel transfrom for complex-valued f(x).
  ! The Hankel transform, in toms794, only supports real-valued f(x).
  ! But, x must be real-valued in all cases.
  SUBROUTINE zhankel(zhankel_kernel, xa, xb, j, gxi, ht, ierr)
    COMPLEX(KIND=dp), EXTERNAL      :: zhankel_kernel
    INTEGER, INTENT(IN)             :: j
    REAL(KIND=dp), INTENT(IN)       :: gxi, xa, xb
    COMPLEX(KIND=dp), INTENT(OUT)   :: ht
    INTEGER, INTENT(OUT)            :: ierr

    REAL(KIND=dp)                   :: gnu, ht_real, ht_imag

    gnu = REAL(j, KIND=dp)
    IF ( gxi == 0.0_dp ) THEN
      ht = CMPLX(0.0_dp, 0.0_dp)
    ELSE
      ! get the hankel transform of the real part of f(x)
      CALL HANKEL_TRANSFORM(HANKEL_KERNEL_REAL, gNU, gXI, XB, HT_REAL, IERR)

      ! get the hankel transform of the imaginary part of f(x)
      CALL HANKEL_TRANSFORM(HANKEL_KERNEL_IMAG, gNU, gXI, XB, HT_IMAG, IERR)

      ht = CMPLX(ht_real, ht_imag)
    END IF

  CONTAINS
    FUNCTION hankel_kernel_real(xi) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: xi

      REAL(KIND=dp)                   :: res

      res = REAL(zhankel_kernel(xi))
    END FUNCTION hankel_kernel_real

    FUNCTION hankel_kernel_imag(xi) RESULT(res)
      REAL(KIND=dp), INTENT(IN)       :: xi

      REAL(KIND=dp)                   :: res

      res = AIMAG(zhankel_kernel(xi))
    END FUNCTION hankel_kernel_imag

  END SUBROUTINE zhankel

  FUNCTION hankel_kernel(xi) RESULT(res)
    REAL(KIND=dp), INTENT(IN)       :: xi

    REAL(KIND=dp)                   :: res

    ! res = REAL(zhankel(xi))
    res = xi**2.0_dp
  END FUNCTION hankel_kernel


END MODULE hankel_api
