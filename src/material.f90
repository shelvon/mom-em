! MODULE: material
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! The optical properties of the material
MODULE material
  USE constants
  USE data
  USE aux

  IMPLICIT NONE

CONTAINS
  FUNCTION media_eps( media, zwl ) RESULT(value)
    TYPE(media_type), INTENT(IN)  :: media
    COMPLEX(KIND=dp), INTENT(IN)  :: zwl

    COMPLEX(KIND=dp)              :: value
    IF ( media%ri%active ) THEN
      value = get_ri( media%ri, REAL(zwl) )
      value = CMPLX( REAL(value)**2-AIMAG(value)**2,  &
                     2.0_dp*REAL(value)*AIMAG(value), KIND=dp )
    ELSE IF ( media%eps%active ) THEN
      value = get_eps( media%eps, TRIM(media%name), zwl )
    ELSE
      WRITE(*,*) '  --Neither media%ri nor media%eps is set up!'
      STOP
    END IF

  END FUNCTION media_eps

!  FUNCTION media_ri( media, zwl) RESULT(value)
!    TYPE(media_type), INTENT(IN)  :: media
!    COMPLEX(KIND=dp), INTENT(IN)  :: zwl
!
!    COMPLEX(KIND=dp)              :: value
!
!    IF ( media%ri%active ) THEN
!      value = get_ri( media%ri, REAL(zwl) )
!    ELSE IF ( media%eps%active ) THEN
!      value = get_eps( media%eps, TRIM(media%name), zwl )
!    END IF
!  END FUNCTION media_ri

  !> #epsilon
  FUNCTION get_eps( eps, name, zwl ) RESULT(value)
    CHARACTER(*), INTENT(IN)            :: name
    TYPE(material_property), INTENT(IN) :: eps
    COMPLEX(KIND=dp),INTENT(IN)         :: zwl
    COMPLEX(KIND=dp)                    :: value

    IF ( TRIM(eps%method) == 'model') THEN
      value = eps_model( TRIM(eps%model), name, 2.0_dp*pi/zwl*c0 )
    ELSE IF ( TRIM(eps%method) == 'value') THEN
      value = eps%z(1)
    ELSE IF ( TRIM(eps%method) == 'table' ) THEN
      WRITE(*,*) '  --Obtain epsilon by table method is not supported!'
!      value = get_ri_table( TRIM(eps%ref_file)//'.ref', REAL(zwl) )
!      value = CMPLX( REAL(value)**2-AIMAG(value)**2,  &
!                     2.0_dp*REAL(value)*AIMAG(value), KIND=dp )
    END IF
  END FUNCTION get_eps

  FUNCTION eps_model( model, name, omega ) RESULT(value)
    CHARACTER(*), INTENT(IN)            :: model
    CHARACTER(*), INTENT(IN)            :: name
    COMPLEX(KIND=dp),INTENT(IN)         :: omega
    COMPLEX(KIND=dp)                    :: value

    COMPLEX(KIND=dp)                          :: omega_eV
    REAL(KIND=dp)                             :: eps_inf, sigma2eps0
    INTEGER                                   :: p,ip
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: AL, BL, CL
    REAL(KIND=dp)                             :: omegaD, gammaD
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: ACP, phiCP, OmegaCP, GammaCP

    value = (1.0_dp, 0.0_dp)
    IF ( model == 'L4') THEN
      omega_eV = hplank*omega/(2.0_dp*pi)
      p = 4
      ! 4 Lorentzian-pole pairs (L4) model, with parameters taken from
      ! F. Hao and P. Nordlander, Chemical Physics Letters 446, 115 (2007).
      ! Note: use j = -i = (0, -1) in Eq. (1)
      ALLOCATE(AL(1:p), BL(1:p), CL(1:p))
      IF ( (name == 'gold') .OR. (name == 'Au') ) THEN
        eps_inf = 1.0_dp
        sigma2eps0 = 1355.01_dp
        AL = (/-8.577E4_dp, -2.875_dp, -997.6_dp, -1.630_dp/)
        BL = (/-1.156E4_dp, 0.0_dp, -3090.0_dp, -4.409_dp/)
        CL = (/5.557E7_dp, 2.079E3_dp, 6.921E5_dp, 26.15_dp/)
      ELSE IF ( (name == 'silver') .OR. (name == 'Ag') ) THEN
        eps_inf = 1.0_dp
        sigma2eps0 = 3157.56_dp
        AL = (/-1.160E5_dp, -4.252_dp, -0.496_dp, -2.118_dp/)
        BL = (/-3050.0_dp, -0.8385_dp, -13.85_dp, -10.23_dp/)
        CL = (/3.634E8_dp, 112.2_dp, 1.815_dp, 14.31_dp/)
      ELSE
        WRITE(*,*) '  --',model,' model for ',name,' is not implemented yet.'
        value = value/0.0_dp
      END IF

      value = eps_inf + sigma2eps0/((0,-1)*omega_eV)
      DO ip = 1, p
        value = value + CL(ip)/(omega_eV**2 + AL(ip)*(0,-1)*omega_eV + BL(ip))
      END DO! ip
      DEALLOCATE( AL, BL, CL )
    ELSE IF ( model == 'D+2CP') THEN
      p = 2
      ! Drude+2 Critical points (2CP) model, with parameters taken from
      ! A. Vial and T. Laroche, Applied Physics B 93, 139 (2008).
      ALLOCATE(ACP(1:p), phiCP(1:p), OmegaCP(1:p), GammaCP(1:p))
      IF ( (name == 'gold') .OR. (name == 'Au') ) THEN
        eps_inf = 1.1431_dp
        omegaD = 1.3202E16_dp
        gammaD = 1.0805E14_dp
        ACP = (/0.26698_dp, 3.0834_dp/)
        phiCP = (/-1.2371_dp, -1.0968_dp/)
        OmegaCP = (/3.8711E15_dp, 4.1684E15_dp/)
        GammaCP = (/4.4642E14_dp, 2.3555E15_dp/)
      ELSE IF ( (name == 'silver') .OR. (name == 'Ag') ) THEN
        eps_inf = 15.833_dp
        omegaD = 1.3861E16_dp
        gammaD = 4.5841E13_dp
        ACP = (/1.0171_dp, 15.797_dp/)
        phiCP = (/-0.93935_dp, 1.8087_dp/)
        OmegaCP = (/6.6327E15_dp, 9.2726E17_dp/)
        GammaCP = (/1.6666E15_dp, 2.3716E17_dp/)
      ELSE
        WRITE(*,*) '  --',model,' model for ',name,' is not implemented yet.'
        value = value/0.0_dp
      END IF

      value = eps_inf - Drude(omega, omegaD, gammaD) +  &
              CP(omega, p, ACP, phiCP, OmegaCP, GammaCP)
      DEALLOCATE( ACP, phiCP, OmegaCP, GammaCP )
    ELSE IF ( model == 'Drude') THEN
      ! Drude model, with parameters taken from
      ! A. Vial and T. Laroche, Applied Physics B 93, 139 (2008).
      IF ( (name == 'gold') .OR. (name == 'Au') ) THEN
        eps_inf = 1.1431_dp
        omegaD = 1.3202E16_dp
        gammaD = 1.0805E14_dp
      ELSE IF ( (name == 'silver') .OR. (name == 'Ag') ) THEN
        eps_inf = 15.833_dp
        omegaD = 1.3861E16_dp
        gammaD = 4.5841E13_dp
      ELSE
        WRITE(*,*) '  --',model,' model for ',name,' is not implemented yet.'
        value = value/0.0_dp
      END IF

      value = eps_inf - Drude(omega, omegaD, gammaD)
    ELSE
      WRITE(*,*) '  --',model,' model is not implemented yet for this media.'
      value = value/0.0_dp
    END IF
  END FUNCTION eps_model

  ! Evaluates the Drude model of noble metal dispersion at frequency omega.
  PURE FUNCTION L4(omega, omegaD, gammaD) RESULT(value)
    COMPLEX(KIND=dp), INTENT(IN)  :: omega
    REAL(KIND=dp), INTENT(IN)     :: omegaD, gammaD
    COMPLEX(KIND=dp)              :: value

    value = omegaD**2/(omega*(omega + (0,1)*gammaD))
  END FUNCTION L4

  ! Evaluates the Drude model of noble metal dispersion at frequency omega.
  PURE FUNCTION Drude(omega, omegaD, gammaD) RESULT(value)
    COMPLEX(KIND=dp), INTENT(IN)  :: omega
    REAL(KIND=dp), INTENT(IN)     :: omegaD, gammaD
    COMPLEX(KIND=dp)              :: value

    value = omegaD**2/(omega*(omega + (0,1)*gammaD))
  END FUNCTION Drude

  ! Evaluates the Critical points of noble metal dispersion at frequency omega.
  PURE FUNCTION CP(omega, p, ACP, phiCP, OmegaCP, GammaCP) RESULT(value)
    COMPLEX(KIND=dp), INTENT(IN)  :: omega
    INTEGER, INTENT(IN)           :: p
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: ACP, phiCP, OmegaCP, GammaCP
    COMPLEX(KIND=dp)              :: value

    INTEGER                       :: ip

    value = (0.0_dp, 0.0_dp)
    DO ip = 1, p
      value = value + ACP(ip)*OmegaCP(ip)*( &
        EXP((0,1)*phiCP(ip))/(OmegaCP(ip)-omega-(0,1)*GammaCP(ip)) + &
        EXP((0,-1)*phiCP(ip))/(OmegaCP(ip)+omega+(0,1)*GammaCP(ip)) )
    END DO! ip
  END FUNCTION CP

  !> #ri
  FUNCTION get_ri(ri, wl) RESULT(res)
    TYPE(material_property), INTENT(IN) :: ri
    REAL(KIND=dp), INTENT(IN)           :: wl
    COMPLEX(KIND=dp)                    :: res

    CHARACTER(LEN=32)                   :: ref_file

    IF ( ri%method == 'value' ) THEN
      res = ri%z(1)
    ELSE IF ( ri%method == 'table' ) THEN
      res = get_ri_table(TRIM(ri%ref_file) // '.ref', wl)
    ELSE IF ( ri%method == 'model' ) THEN
      WRITE(*,*) '  --Obtain refractive index by model method is not supported!'
      STOP
    ELSE
      WRITE(*,*) 'ri%method not supported, stop!'
      STOP
    END IF
  END FUNCTION get_ri

  ! Returns a complex refractive index value from the file at given wavelength.
  ! Linear interpolation is used.
  FUNCTION get_ri_table(filename, wl) RESULT(res)
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    REAL(KIND=dp), INTENT(IN)     :: wl
    COMPLEX(KIND=dp)              :: res
    REAL(KIND=dp)                 :: wl1, wl2, nr1, ni1, nr2, ni2, t
    INTEGER                       :: fid = 10, iovar, i
    INTEGER, PARAMETER            :: nattempts = 10

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
       WRITE(*,*) 'Could not open refractive index lookup file!'
       STOP
    END IF

    iovar = 0
    READ(fid,*) wl1, nr1, ni1

    IF(wl<wl1) THEN
      WRITE(*,*) 'Wavelength below refractive index sample range!'
      WRITE(*,*) wl, ' < ', wl1
      WRITE(*,*) 'Using unity refractive index.'
      res = 1.0_dp
      CLOSE(fid)
      RETURN
    END IF

    DO WHILE(iovar==0)
      READ(fid,*) wl2, nr2, ni2
      ! IF( (wl1<=wl) .AND. (wl<=wl2) .AND. (wl1/=wl2)) THEN
      IF ( (wl1<=wl) .AND. (wl<=wl2) .AND. (ABS(wl1-wl2)>=epsilon_dp) ) THEN
        t = (wl-wl1)/(wl2 - wl1)
        res = CMPLX(linterp(nr1, nr2, t), linterp(ni1, ni2, t), KIND=dp)
        CLOSE(fid)
        RETURN
      ELSE IF(wl==wl1) THEN
        res = CMPLX(nr1, ni1, KIND=dp)
        CLOSE(fid)
        RETURN
      ELSE IF(wl==wl2) THEN
        res = CMPLX(nr2, ni2, KIND=dp)
        CLOSE(fid)
        RETURN
      ELSE IF(wl1>wl2) THEN
        WRITE(*,*) 'Refractive index data disordered!'
        CLOSE(fid)
        STOP
      END IF

      wl1 = wl2
      nr1 = nr2
      ni1 = ni2
    END DO

    CLOSE(fid)

    WRITE(*,*) 'Wavelength above refractive index sample range!'
    WRITE(*,*) 'Using unity refractive index.'
    res = 1.0_dp
  END FUNCTION get_ri_table

  FUNCTION media_mu( media, zwl ) RESULT(value)
    TYPE(media_type), INTENT(IN)  :: media
    COMPLEX(KIND=dp), INTENT(IN)  :: zwl

    COMPLEX(KIND=dp)              :: value

    ! nonmagnetic material
    value = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
  END FUNCTION

END MODULE material
