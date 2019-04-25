! MODULE: io
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! Do the input/output operations.
MODULE io
  USE fson
  ! Functions for accessing data as an array
  USE fson_value_m, only: fson_value_count, fson_value_get
  ! Typical usage should only require an explicit use of the fson module.
  ! The other modules will be used privatley by fson as required.
  USE data
  USE aux
!#ifdef _USE_HDF5
  USE h5_wrapper
!#endif

  IMPLICIT NONE

  INTERFACE write_matrix
    MODULE PROCEDURE write_matrix1d, write_matrix2d
  END INTERFACE

CONTAINS
  ! fson_value to model
  SUBROUTINE json2model(json_file, model)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)        :: json_file
    TYPE(model_type), INTENT(INOUT)     :: model

    TYPE(fson_value), POINTER           :: json_root, json_item

    ! parse the json file
    json_root => fson_parse(json_file)
    ! print the parsed data to the console
    !CALL fson_print(json_root)

    ! fson_value to model%name
    CALL fson_get(json_root, "name", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "name", model%name) ! required
      WRITE(*,*) 'Node model%name is set to "', TRIM(model%name), '".'
    ELSE
      WRITE(*,*) 'Node model%name is not found!'
      STOP ERR0
    END IF

    ! fson_value to model%geom
    CALL fson_get(json_root,"geom", json_item) ! optional, not needed in source calculation
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2geom(json_item, model%geom)
      WRITE(*,*) 'Node model%geom is loaded.'
    END IF

    ! fson_value to model%physics
    CALL fson_get(json_root,"physics", json_item) ! partially optional, only media(0:1) is required for source calculation
    IF ( ASSOCIATED(json_item) ) THEN
      WRITE(*,*) 'Node model%physics is loading ...'
      CALL json2physics(json_item, model%physics)
      WRITE(*,*) 'Node model%physics is loaded.'
    ELSE
      WRITE(*,*) 'Node model%physics is not found!'
      STOP ERR0
    END IF

    ! fson_value to model%solver
    CALL fson_get(json_root,"simulation", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2simulation(json_item, model%simulation)
      WRITE(*,*) 'Node model%simulation is loaded.'
    ELSE
      WRITE(*,*) 'Node model%simulation is not found!'
      STOP ERR0
    END IF

    ! fson_value to model%solution
    CALL fson_get(json_root,"solution", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2solution(json_item, model%solution, SIZE(model%simulation%zwl) )
      WRITE(*,*) 'Node model%solution is loaded.'
    ELSE
      WRITE(*,*) 'Node model%solution is not found!'
      STOP ERR0
    END IF

    ! clean up
    CALL fson_destroy(json_root)
  END SUBROUTINE json2model

  ! fson_value to geom
  SUBROUTINE json2geom(json_root, geom)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(geom_type), INTENT(INOUT)            :: geom

    TYPE(fson_value), POINTER                 :: json_item, json_array
    INTEGER                                   :: ndom, idom
    CHARACTER(LEN=16)                         :: qd_rule

    ! fson_value to geom%mesh
    CALL fson_get(json_root,"mesh", json_item)
    CALL json2mesh(json_item, geom%mesh)

    ! fson_value to geom%domain
    CALL fson_get(json_root, "domain", json_array)
    ndom = fson_value_count(json_array)
    ! index starts from -1, i.e.
    ! -1: all domains, i.e. the entire system;
    ! 0: air/vacuum (or the surrounding media);
    ! 1: the first scatter;
    ! etc.
    ALLOCATE(geom%domain(-1:ndom-1))

    ! set global [domain(-1)] quadrature rule
    CALL fson_get(json_root, "qd_tri", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "qd_tri", geom%domain(-1)%qdrule_tri)
    ELSE
      geom%domain(-1)%qdrule_tri = 'tri_gl13'
    END IF

    CALL fson_get(json_root, "qd_tetra", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "qd_tetra", geom%domain(-1)%qdrule_tetra)
    ELSE
      geom%domain(-1)%qdrule_tetra = 'tetra_gl4'
    END IF

    DO idom = 1, ndom
      ! Get the array item (this is an associative array)
      json_item => fson_value_get(json_array, idom)
      ! array in json file is indexed from 0,
      ! but fson_value_get function is indexed from 1

      CALL json2domain(json_item, geom%domain(idom-1), geom%domain(-1))
    END DO
  END SUBROUTINE json2geom

  ! fson_value to mesh
  SUBROUTINE json2mesh(json_root, mesh)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(mesh_type), INTENT(INOUT)            :: mesh

    ! initial type and file string to null
    mesh%type = ''
    mesh%file = ''
    CALL fson_get(json_root, "type", mesh%type)
    CALL fson_get(json_root, "file", mesh%file)
    CALL fson_get(json_root, "scale", mesh%scale)

  END SUBROUTINE json2mesh

  ! fson_value to domain
  SUBROUTINE json2domain(json_root, domain, domainAll)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(domain_type), INTENT(INOUT)          :: domain
    TYPE(domain_type), INTENT(IN)             :: domainAll

    TYPE(fson_value), POINTER                 :: json_item, json_array
    INTEGER                                   :: nsurface, nvolume
    CHARACTER(LEN=16)                         :: qd_rule

    ! get media index
    CALL fson_get(json_root, "media", domain%media)

    ! get surface labels
    CALL fson_get(json_root, "surface", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nsurface=fson_value_count(json_array)
      ALLOCATE(domain%surface(1:nsurface))
      CALL fson_get(json_array, "", domain%surface)
      WRITE(*,*) '  Node geom%domain%surface is loaded.'
    ELSE
      WRITE(*,*) '  Node geom%domain%surface is not found!'
      STOP ERR0
    END IF

    ! get volume labels
    CALL fson_get(json_root, "volume", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nvolume=fson_value_count(json_array)
      ALLOCATE(domain%volume(1:nvolume))
      CALL fson_get(json_array, "", domain%volume)
    ELSE
      WRITE(*,*) '  Node geom%domain%volume is not found!'
    END IF

    ! set global [domain(-1)] quadrature rule
    CALL fson_get(json_root, "qd_tri", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "qd_tri", domain%qdrule_tri)
    ELSE
      domain%qdrule_tri = domainAll%qdrule_tri
    END IF

    CALL fson_get(json_root, "qd_tetra", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "qd_tetra", domain%qdrule_tetra)
    ELSE
      domain%qdrule_tetra = domainAll%qdrule_tetra
    END IF

  END SUBROUTINE json2domain

  ! json_value to physics
  SUBROUTINE json2physics(json_root, physics)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(physics_type), INTENT(INOUT)         :: physics

    TYPE(fson_value), POINTER                 :: json_item, json_array
    INTEGER     :: nmedia, imedia, nsource, isource, nga

    ! json_value to physics%symmetry
    CALL fson_get(json_root,"group", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nga=fson_value_count(json_array)
      ALLOCATE(physics%group%name(1:nga))
      CALL fson_get(json_array, "", physics%group%name)
      WRITE(*,*) '  Node physics%group is loaded.'
    ELSE
      ALLOCATE(physics%group%name(1))
      physics%group%name(1) = 'id'
    END IF

    ! json_value to physics%media
    CALL fson_get(json_root, "media",json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nmedia=fson_value_count(json_array)
      ALLOCATE(physics%media(0:nmedia-1))
      DO imedia = 1, nmedia
        ! Get the array item (this is an associative array)
        json_item => fson_value_get(json_array, imedia)
        ! array in json file is indexed from 0,
        ! but fson_value_get function is indexed from 1

        CALL json2media(json_item, physics%media(imedia-1))
      END DO! imedia
      WRITE(*,*) '  Node physics%media is loaded.'
    ELSE
      WRITE(*,*) '  Node physics%media is not found!'
      STOP ERR0
    END IF

    ! json_value to physics%source
    CALL fson_get(json_root, "source",json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nsource=fson_value_count(json_array)
      ALLOCATE(physics%source(1:nsource)) ! index from 1
      DO isource = 1, nsource
        json_item => fson_value_get(json_array, isource)
        CALL json2source(json_item, physics%source(isource))
      END DO! isource
      WRITE(*,*) '  Node physics%source is loaded.'
    ELSE
      WRITE(*,*) '  Node physics%source is not found! Only for mode solver!'
    END IF

  END SUBROUTINE json2physics

  ! json_value to symmetry
!  SUBROUTINE json2symmetry(json_root, symmetry)
!    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
!    TYPE(symmetry_type), INTENT(INOUT)        :: symmetry
!
!    CALL fson_get(json_root, "nsubgroups", symmetry%nsubgroups)
!    ALLOCATE(symmetry%group_name(1:symmetry%nsubgroups))
!    CALL fson_get(json_root, "group_name", symmetry%group_name)
!  END SUBROUTINE json2symmetry

  ! json_value to media
  SUBROUTINE json2media(json_root, media)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(media_type), INTENT(INOUT)           :: media

    ! get media name
    CALL fson_get(json_root, "name", media%name)

    ! get epsilon
    CALL json2property(json_root, "epsilon", media%eps)
    IF (media%eps%active) THEN
      media%ri%active = .FALSE.
    END IF

    ! get ri, i.e. refractive index
    CALL json2property(json_root, "ri", media%ri)
    IF (media%ri%active) THEN
      media%eps%active = .FALSE.
    END IF

    ! set default value for permittivity through refractive index
    IF ( (.NOT. media%ri%active) .AND. (.NOT. media%eps%active) ) THEN
      media%ri%active = .TRUE.
      media%ri%method = 'value'
      ! non-dispersive, i.e. constant value
      ALLOCATE(media%ri%re(1))
      ALLOCATE(media%ri%im(1))
      ALLOCATE(media%ri%z(1))
      media%ri%re(1) = 1.0_dp
      media%ri%im(1) = 0.0_dp
      media%ri%z(1) = CMPLX(0.0, 0.0, KIND=dp)
    END IF

    ! get mu
    CALL json2property(json_root, "mu", media%mu)

  END SUBROUTINE json2media

  ! json_value to source
  SUBROUTINE json2source(json_root, source)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    ! TYPE(srcdata), INTENT(INOUT)              :: source
    TYPE(source_type), INTENT(INOUT)          :: source

    TYPE(fson_value), POINTER                 :: json_item

    ! get media name
    CALL fson_get(json_root, "type", source%type)
    IF (source%type == 'pw') THEN
      CALL fson_get(json_root, "pw",json_item)
      IF ( ASSOCIATED(json_item) ) THEN
        CALL json2pw(json_item, source%pw)
        WRITE(*,*) '  Node physics%source%pw is loaded.'
      ELSE
        WRITE(*,*) '  Node physics%source%pw is not found!'
        STOP ERR0
      END IF

    ELSE IF (source%type == 'focus') THEN
      CALL fson_get(json_root, "focus",json_item)
      IF ( ASSOCIATED(json_item) ) THEN
        CALL json2focus(json_item, source%focus)
        WRITE(*,*) '  Node physics%source%focus is loaded.'
      ELSE
        WRITE(*,*) '  Node physics%source%focus is not found!'
        STOP ERR0
      END IF

    END IF
  END SUBROUTINE json2source

  SUBROUTINE json2pw(json_root, pw)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(pw_type), INTENT(INOUT)              :: pw
    TYPE(fson_value), POINTER                 :: json_item

    ! get incident angle, required
    CALL fson_get(json_root, "theta", pw%theta)
    CALL fson_get(json_root, "phi", pw%phi)
    ! get relative phase difference
    CALL fson_get(json_root, "phase", pw%phase)
    IF ( ABS(pw%phase) .EQ. 90) THEN
      ! circularly polarized light
      pw%A = 1.0_dp/SQRT(2.0_dp)
      pw%B = 1.0_dp/SQRT(2.0_dp)
      WRITE(*,*) 'circularly polarized plane wave!'
    ELSE
      ! get either psi or A, B values
      CALL fson_get(json_root, "psi", json_item)
      IF (associated(json_item)) THEN
        ! polarization direction reads from psi angle!
        CALL fson_get(json_root, "psi", pw%psi)
      ELSE
        CALL fson_get(json_root, "A", json_item)
        IF (associated(json_item)) THEN
          ! get electric strength of two orthogonal components
          CALL fson_get(json_root, "A", pw%A)
          CALL fson_get(json_root, "B", pw%B)
        ELSE
          WRITE(*,*) 'error input parameters, either psi or A/B should be defined!'
        END IF
      END IF
    END IF

  END SUBROUTINE json2pw

  SUBROUTINE json2focus(json_root, focus)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(focus_type), INTENT(INOUT)           :: focus
    TYPE(fson_value), POINTER                 :: json_item, json_item_sub

    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: pos

    CALL fson_get(json_root, "paraxial", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "paraxial", focus%paraxial)
    ELSE
      focus%paraxial = .FALSE.
    END IF

    ! required nodes
    CALL fson_get(json_root, "input", focus%input)! input beam
    CALL fson_get(json_root, "focal", focus%focal)
    CALL fson_get(json_root, "waist", focus%waist)
    CALL fson_get(json_root, "na", focus%na)
    CALL fson_get(json_root, "norm", focus%norm)

    ! optional nodes, inquire the existence first
    CALL fson_get(json_root, "E0", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "E0", focus%E0)
    ELSE
      focus%E0 = 1.0_dp
    END IF

    CALL fson_get(json_root, "pos", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "pos", pos)
      focus%pos = pos
    ELSE
      focus%pos = (/0.0_dp, 0.0_dp, 0.0_dp/)
    END IF

    ! read the input beam parameters
    SELECT CASE (focus%input)
      CASE ('pupil')
        ! json_value to focus%pupil
        CALL fson_get(json_root, "pupil", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL json2pupil(json_item, focus%pupil)
          WRITE(*,*) '  Node physics%source%focus%pupil is loaded.'
        END IF

        ! Theoretically, it seems this inversed beam doesn't work!
!      ! inversed basis
!      ! Assume the basis polarization states in the focal region,
!      ! then deduce back the incident beam parameters.
!      CASE ('invbasis')
!        CALL fson_get(json_root, "invbasis", json_item)
!        IF ( ASSOCIATED(json_item) ) THEN
!          CALL fson_get(json_item, "n", focus%invbasis%n)
!          CALL fson_get(json_item, "pol", focus%invbasis%pol)
!          WRITE(*,*) '  Node physics%source%focus%invbasis is loaded.'
!        END IF

      CASE ('basis')
        CALL fson_get(json_root, "basis", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "n", focus%basis%n)
          CALL fson_get(json_item, "pol", focus%basis%pol)
          WRITE(*,*) '  Node physics%source%focus%basis is loaded.'
        END IF

      CASE ('cos')
        CALL fson_get(json_root, "cos", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "l", focus%cos%l)
          CALL fson_get(json_item, "pol", focus%cos%pol)
          WRITE(*,*) '  Node physics%source%focus%cos is loaded.'
        END IF

      CASE ('sin')
        CALL fson_get(json_root, "sin", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "l", focus%sin%l)
          CALL fson_get(json_item, "pol", focus%sin%pol)
          WRITE(*,*) '  Node physics%source%focus%sin is loaded.'
        END IF

      CASE ('petal')
        CALL fson_get(json_root, "petal", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "pol", focus%petal%pol)
          CALL fson_get(json_item, "n", focus%petal%n)
          CALL fson_get(json_item, "l", focus%petal%l)
          WRITE(*,*) '  Node physics%source%focus%petal is loaded.'
        END IF

      CASE ('lg')
        CALL fson_get(json_root, "lg", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "pol", focus%lg%pol)
          CALL fson_get(json_item, "n", focus%lg%n)
          CALL fson_get(json_item, "l", focus%lg%l)
          WRITE(*,*) '  Node physics%source%focus%lg is loaded.'
        END IF

      CASE ('hg')
        CALL fson_get(json_root, "hg", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "m", focus%hg%m)
          CALL fson_get(json_item, "n", focus%hg%n)
          CALL fson_get(json_item, "pol", focus%hg%pol)
          WRITE(*,*) '  Node physics%source%focus%hg is loaded.'
        END IF

      CASE ('rect')
        CALL fson_get(json_root, "rect", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL json2rect(json_item, focus%rect)
          WRITE(*,*) '  Node physics%source%focus%rect is loaded.'
        END IF

      CASE ('bessel')
        CALL fson_get(json_root, "bessel", json_item)
        IF ( ASSOCIATED(json_item) ) THEN
          CALL fson_get(json_item, "pol", focus%bessel%pol)
          CALL fson_get(json_item, "j", json_item_sub)
          IF ( ASSOCIATED(json_item_sub) ) THEN
            CALL fson_get(json_item, "j", focus%bessel%j)
          ELSE
            focus%bessel%j = 0
          END IF
          CALL fson_get(json_item, "delta", focus%bessel%delta)
          CALL fson_get(json_item, "theta", focus%bessel%theta)
          WRITE(*,*) '  Node physics%source%focus%bessel is loaded.'
        END IF

    END SELECT
  END SUBROUTINE json2focus

  SUBROUTINE json2pupil(json_root, pupil)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(pupil_type), INTENT(INOUT)           :: pupil

    TYPE(fson_value), POINTER                 :: json_item

    !--read radial modulation, R_q
    CALL fson_get(json_root, "R_r", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funR(json_item, pupil%R(1))
      WRITE(*,*) '  Node physics%source%focus%pupil%R_r is loaded.'
    ELSE
      pupil%R(1)%fun = 'const'
      pupil%R(1)%const = CMPLX(1, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%R_r is not found, set to 1.'
    END IF

    CALL fson_get(json_root, "R_phi", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funR(json_item, pupil%R(2))
      WRITE(*,*) '  Node physics%source%focus%pupil%R_phi is loaded.'
    ELSE
      pupil%R(2)%fun = 'const'
      pupil%R(2)%const = CMPLX(1, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%R_phi is not found, set to 1.'
    END IF

    CALL fson_get(json_root, "R_z", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funR(json_item, pupil%R(3))
      WRITE(*,*) '  Node physics%source%focus%pupil%R_z is loaded.'
    ELSE
      ! By default, force Ez = 0, i.e., the input beam is pure transversal.
      pupil%R(3)%fun = 'const'
      pupil%R(3)%const = CMPLX(1, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%R_z is not found, set to 1.'
    END IF

    !--read azimuthal modulation, Phi_q
    CALL fson_get(json_root, "Phi_r", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funPhi(json_item, pupil%Phi(1))
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_r is loaded.'
    ELSE
      pupil%Phi(1)%fun = 'series'
      ALLOCATE( pupil%Phi(1)%cn(0:0) )
      pupil%Phi(1)%cn = CMPLX(0, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_r is not found'
    END IF

    CALL fson_get(json_root, "Phi_phi", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funPhi(json_item, pupil%Phi(2))
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_phi is loaded.'
    ELSE
      pupil%Phi(2)%fun = 'series'
      ALLOCATE( pupil%Phi(2)%cn(0:0) )
      pupil%Phi(2)%cn = CMPLX(0, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_phi is not found'
    END IF

    CALL fson_get(json_root, "Phi_z", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2funPhi(json_item, pupil%Phi(3))
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_z is loaded.'
    ELSE
      pupil%Phi(3)%fun = 'series'
      ALLOCATE( pupil%Phi(3)%cn(0:0) )
      ! By default, force Ez = 0, i.e., the input beam is pure transversal.
      pupil%Phi(3)%cn = CMPLX(0, 0, KIND=dp)
      WRITE(*,*) '  Node physics%source%focus%pupil%Phi_z is not found'
    END IF
  END SUBROUTINE json2pupil

  SUBROUTINE json2funR(json_root, funR)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(funR_type), INTENT(INOUT)            :: funR

    REAL(KIND=dp)                             :: Rconst_re, Rconst_im
    CALL fson_get(json_root, "fun", funR%fun)
    IF ( funR%fun == 'const' ) THEN
      CALL fson_get(json_root, "const.re", Rconst_re)
      CALL fson_get(json_root, "const.im", Rconst_im)
      funR%const = CMPLX(Rconst_re, Rconst_im, KIND=dp)
    ELSE IF ( funR%fun == 'sin' ) THEN
      CALL fson_get(json_root, "sin.x", funR%sin%x)
      CALL fson_get(json_root, "sin.y", funR%sin%y)
    ELSE IF ( funR%fun == 'cos' ) THEN
      CALL fson_get(json_root, "cos.x", funR%cos%x)
      CALL fson_get(json_root, "cos.y", funR%cos%y)
    END IF

  END SUBROUTINE json2funR

  SUBROUTINE json2funPhi(json_root, funPhi)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(funPhi_type), INTENT(INOUT)          :: funPhi

    TYPE(fson_value), POINTER                 :: json_item
    INTEGER                                   :: icn
    INTEGER, ALLOCATABLE, DIMENSION(:)        :: cn_n
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: cn_re, cn_im

    CALL fson_get(json_root, "fun", funPhi%fun)
    IF ( funPhi%fun == 'series' ) THEN
      ! Get expansion coefficients 'cn' from user input directly.
      CALL fson_get(json_root, "cn.n", cn_n)
      CALL fson_get(json_root, "cn.re", cn_re)
      CALL fson_get(json_root, "cn.im", cn_im)

      IF ( (SIZE(cn_n) .NE. SIZE(cn_re)) &
          .OR. (SIZE(cn_n) .NE. SIZE(cn_im)) ) THEN
        WRITE(*,*) ' Inconsistent number of series terms, STOP!'
        STOP ERR0
      END IF
      ALLOCATE( funPhi%cn(MINVAL(cn_n):MAXVAL(cn_n)) )
      DO icn = 1, SIZE(cn_n)
        funPhi%cn(cn_n(icn)) = CMPLX(cn_re(icn),cn_im(icn), KIND=dp)
      END DO!

    ! Otherwise, calculate expansion coefficients 'cn' later on.
    ELSE IF ( funPhi%fun == 'rect' ) THEN
      CALL fson_get(json_root, "rect", json_item)
      IF ( ASSOCIATED(json_item) ) THEN
        CALL json2rect(json_item, funPhi%rect)
        WRITE(*,*) '  Node funPhi%rect is loaded.'
      ELSE
        WRITE(*,*) '  Node funPhi%rect is not found, STOP!'
        STOP ERR0
      END IF
    END IF
  END SUBROUTINE json2funPhi

  SUBROUTINE json2rect(json_root, rect)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(rect_type), INTENT(INOUT)            :: rect

    TYPE(fson_value), POINTER                 :: json_array,json_item
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: intervals
    INTEGER                                   :: narray,iarray

    CALL fson_get(json_root, "pol", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "pol", rect%pol)
      WRITE(*,*) '  Node rect%pol is loaded.'
    END IF

    CALL fson_get(json_root, "delta", rect%delta)
    rect%intervals_n = SIZE(rect%delta)
    CALL fson_get(json_root, "series_n", rect%series_n)
    ALLOCATE( rect%intervals(2, rect%intervals_n) )
    CALL fson_get(json_root, "intervals", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      narray = fson_value_count(json_array)
      IF ( narray .NE. rect%intervals_n) THEN
        WRITE(*,*) ' Inconsistent number of intervals, STOP!'
        STOP ERR0
      END IF
      DO iarray = 1, narray
        ! Get the array item (this is an associative array)
        json_item => fson_value_get(json_array, iarray)
        ! array in json file is indexed from 0,
        ! but fson_value_get function is indexed from 1

        CALL fson_get(json_item, "", intervals)
        rect%intervals(:,iarray) = intervals
      END DO! iarray
    END IF

  END SUBROUTINE json2rect

  ! json_value to property by name
  SUBROUTINE json2property(json_root, name, property)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    CHARACTER(*), INTENT(IN)                  :: name
    TYPE(material_property), INTENT(INOUT)    :: property

    TYPE(fson_value), POINTER                 :: json_item, json_item_sub
    CHARACTER(32)                             :: json_path

    ! inquiry first into the existence of the corresponding property
    json_path = TRIM(name)
    CALL fson_get(json_root, json_path, json_item)
    IF (associated(json_item)) THEN
      property%active = .TRUE.

      CALL fson_get(json_item, "method", property%method)
      IF ( property%method == 'value' ) THEN
        ! non-dispersive, i.e. constant value
        ALLOCATE(property%re(1))
        ALLOCATE(property%im(1))
        ALLOCATE(property%z(1))
        property%re(1) = 0.0_dp
        property%im(1) = 0.0_dp
        property%z(1) = 0.0_dp
        CALL fson_get(json_item, "real", json_item_sub)
        IF (associated(json_item_sub)) THEN
          CALL fson_get(json_item, "real", property%re(1))
        ELSE
          property%re(1) = 1.0_dp
        END IF

        CALL fson_get(json_item, "imag", json_item_sub)
        IF (associated(json_item_sub)) THEN
          CALL fson_get(json_item, "imag", property%im(1))
        ELSE
          property%im(1) = 0.0_dp
        END IF
        property%z = CMPLX(property%re, property%im, KIND=dp)

      ELSE IF ( property%method == 'table' ) THEN
        CALL fson_get(json_item, "ref_file", property%ref_file)

      ELSE IF (property%method == 'model') THEN
        CALL fson_get(json_item, "model", property%model)
      END IF

    ELSE
      property%active = .FALSE.
    END IF

  END SUBROUTINE json2property

  ! json_value to simulation
  SUBROUTINE json2simulation(json_root, simulation)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(simulation_type), INTENT(INOUT)      :: simulation

    TYPE(fson_value), POINTER                 :: json_item

    CALL fson_get(json_root, "solver", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2solver( json_item, simulation%solver )
      WRITE(*,*) '  Node simulation%solver is loaded.'
    ELSE
      WRITE(*,*) '  Node simulation%solver is not found!'
      STOP ERR0
    END IF

    ! get frquency at which the simulation is running
    CALL fson_get(json_root, "frequency", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2frequency( json_item, simulation )
    ELSE
      WRITE(*,*) '  Node simulation%frequency is not found!'
    END IF

    ! get wavelength at which the simulation is running
    ! json_value to simulation%wl
    CALL fson_get(json_root, "wavelength", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2wavelength( json_item, simulation )
    ELSE
      WRITE(*,*) '  Node simulation%wl is not found!'
    END IF

  END SUBROUTINE json2simulation

  SUBROUTINE json2solver( json_root, solver )
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(solver_type), INTENT(INOUT)          :: solver

    TYPE(fson_value), POINTER                 :: json_item

    CALL fson_get(json_root, "name", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "name", solver%name)
      WRITE(*,*) '  Node solver%name is loaded.'
    ELSE
      WRITE(*,*) '  Node solver%name is not found!'
      STOP ERR0
    END IF

    CALL fson_get(json_root, "parallel", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "parallel", solver%parallel)
      WRITE(*,*) '  Node solver%parallel is loaded.'
    ELSE
      solver%parallel = .FALSE.
    END IF

    CALL fson_get(json_root, "memory", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "memory", solver%memory)
      WRITE(*,*) '  Node solver%memory is loaded.'
    ELSE
      solver%memory = .FALSE.
    END IF

    IF ( TRIM(solver%name) == 'mode' ) THEN
      CALL fson_get(json_root, "poles", json_item)
      IF ( ASSOCIATED(json_item) ) THEN
        CALL fson_get(json_root, "poles", solver%poles)
        WRITE(*,*) '  Node solver%poles is loaded.'
      ELSE
        WRITE(*,*) '  Node solver%poles is not found!'
        solver%poles = 'cauchy'! Defualt poles searching method: cauchy integral.
      END IF
    END IF

  END SUBROUTINE json2solver

  SUBROUTINE json2frequency( json_root, simulation )
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(simulation_type), INTENT(INOUT)      :: simulation

    TYPE(fson_value), POINTER                 :: json_item
    CHARACTER(32)                             :: freq_method
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: freq_value
    COMPLEX(KIND=dp)                          :: zfreqa, zfreqb, zfreqc, zfreqd
    INTEGER                                   :: nfreq, ia, ib, nfreqx, nfreqy
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zfreq_value

    CALL fson_get(json_root, "method", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "method", freq_method)
      WRITE(*,*) '  Loading method for node simulation%frequency is set to'// &
                 '"', TRIM(freq_method),'".'
    ELSE
      WRITE(*,*) '  Loading method for node simulation%frequency not found!'
      STOP ERR0
    END IF

    CALL fson_get(json_root, "value", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "value", freq_value)

      WRITE(*,*) '  Node simulation%frequency is loaded.'
    ELSE
      WRITE(*,*) '  Loading value for node simulation%frequency is not found!'
      STOP ERR0
    END IF

    ! generate frequency/wavelength values
    IF ( TRIM(freq_method) == 'list' ) THEN
      ALLOCATE(simulation%omega(1:SIZE(freq_value,1)))
      simulation%omega = freq_value
    ELSE IF ( TRIM(freq_method) == 'range' ) THEN
      nfreq = NINT(freq_value(1))
      ALLOCATE(simulation%omega(1:nfreq))
      simulation%omega = linspace(freq_value(2), freq_value(3), nfreq)
    ELSE IF ( TRIM(freq_method) == 'zrange' ) THEN
      nfreqx = NINT(freq_value(1))
      nfreqy = NINT(freq_value(2))
      ! corners of a rectangular region
      zfreqa = CMPLX(freq_value(3), freq_value(4), KIND=dp)
      zfreqb = CMPLX(freq_value(5), freq_value(4), KIND=dp)
      zfreqc = CMPLX(freq_value(5), freq_value(6), KIND=dp)
      zfreqd = CMPLX(freq_value(3), freq_value(6), KIND=dp)
      ! frequencies for cauchy's line integral
      IF ( (TRIM(simulation%solver%name) == 'mode') &
        .AND. (TRIM(simulation%solver%poles) == 'cauchy') ) THEN

        ! Allocate points along lines (counterclockwise) for cauchy's integral.
        ! But, there are 4 points at corners that are counted twice.
        ALLOCATE( simulation%zomega(1:(nfreqx*2+nfreqy*2-4)) )
        ia = 1; ib = ia+nfreqx-1;
        simulation%zomega(ia:ib) = zlinspace( zfreqa, zfreqb, nfreqx, 1)
        ia = ib; ib = ia+nfreqy-1;
        simulation%zomega(ia:ib) = zlinspace( zfreqb, zfreqc, 1, nfreqy)
        ia = ib; ib = ia+nfreqx-1;
        simulation%zomega(ia:ib) = zlinspace( zfreqc, zfreqd, nfreqx, 1)
        ia = ib; ib = ia+nfreqy-1;
        ALLOCATE( zfreq_value(ia:ib) )
        zfreq_value(ia:ib) = zlinspace( zfreqd, zfreqa, 1, nfreqy)
        simulation%zomega(ia:ib-1) = zfreq_value(ia:ib-1)
      ELSE
        ALLOCATE(simulation%zomega(1:nfreqx*nfreqy))
        simulation%zomega = zlinspace(  &
                            CMPLX(freq_value(3), freq_value(4), KIND=dp), &
                            CMPLX(freq_value(5), freq_value(6), KIND=dp), &
                            nfreqx, nfreqy)
      END IF
      ! copy to simulation%zwl
      ALLOCATE( simulation%zwl(1:SIZE(simulation%zomega)) )
      simulation%zwl = 2.0_dp*pi*c0/simulation%zomega
    ELSE
      WRITE(*,*) '  Loading method for node simulation%frequency is not supported!'
      STOP ERR0
    END IF

    CALL fson_get(json_root, "starter", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "starter", simulation%iwl)
    ELSE
      simulation%iwl = 1 ! start from 1 by default
    END IF
  END SUBROUTINE json2frequency

  SUBROUTINE json2wavelength( json_root, simulation )
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(simulation_type), INTENT(INOUT)      :: simulation

    TYPE(fson_value), POINTER                 :: json_item
    CHARACTER(32)                             :: wl_method
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: wl_value
    INTEGER                                   :: nwl, nwlx, nwly

    CALL fson_get(json_root, "method", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "method", wl_method)
      WRITE(*,*)  '  Loading method for node simulation%wl is set to "',  &
                  TRIM(wl_method),'".'
    ELSE
      WRITE(*,*) '  Loading method for node simulation%wl is not found!'
      STOP ERR0
    END IF

    CALL fson_get(json_root, "value", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "value", wl_value)
      ! CALL fson_get(json_root, "wavelength.value", wl_json)
      ! ALLOCATE( wl_value(SIZE(wl_json)) )
      ! wl_value = REAL(wl_json, KIND=dp)
      ! WRITE(*,*) wl_json, wl_value

      IF ( TRIM(wl_method) == 'list' ) THEN
        ALLOCATE(simulation%wl(1:SIZE(wl_value,1)))
        simulation%wl = wl_value

        ALLOCATE(simulation%zwl(1:SIZE(wl_value,1)))
        simulation%zwl = CMPLX(wl_value,0.0, KIND=dp)
        WRITE(*,*) simulation%zwl

      ELSE IF ( TRIM(wl_method) == 'range' ) THEN
        nwl = INT(wl_value(1))
        ALLOCATE(simulation%wl(1:nwl))
        simulation%wl = linspace(wl_value(2), wl_value(3), nwl)

        ALLOCATE(simulation%zwl(1:nwl))
        simulation%zwl = zlinspace(CMPLX(wl_value(2), 0.0_dp, KIND=dp),&
                                   CMPLX(wl_value(3), 0.0_dp, KIND=dp),&
                                   nwl, 1)
        WRITE(*,*) simulation%zwl

      ELSE IF ( TRIM(wl_method) == 'zrange' ) THEN
        nwlx = NINT(wl_value(1))
        nwly = NINT(wl_value(2))
        ALLOCATE(simulation%zwl(1:nwlx*nwly))
        simulation%zwl = zlinspace(CMPLX(wl_value(3), wl_value(4), KIND=dp),&
                                   CMPLX(wl_value(5), wl_value(6), KIND=dp),&
                                   nwlx, nwly)
      ELSE
        WRITE(*,*) '  Loading method for node simulation%wl is not supported!'
        STOP ERR0
      END IF

      WRITE(*,*) '  Node simulation%wl is loaded.'

    ELSE
      WRITE(*,*) '  Loading value for node simulation%wl is not found!'
      STOP ERR0
    END IF

    CALL fson_get(json_root, "starter", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "starter", simulation%iwl)
    ELSE
      simulation%iwl = 1 ! start from 1 by default
    END IF
  END SUBROUTINE json2wavelength

  SUBROUTINE json2solution(json_root, solution, nwl)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(solution_type), INTENT(INOUT)        :: solution
    INTEGER, INTENT(IN)                       :: nwl

    TYPE(fson_value), POINTER                 :: json_item, json_item_sub, json_array
    INTEGER                                   :: nfocal, ifocal, isave, nsave
    INTEGER, DIMENSION(:), ALLOCATABLE        :: saveidx

    ALLOCATE( solution%base(nwl) )
    solution%base(:)%save = .FALSE.

    CALL fson_get(json_root, "base", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "base.save", json_item_sub)
      IF ( ASSOCIATED(json_item_sub) ) THEN
        CALL fson_get(json_root, "base.save", saveidx)

        nsave = SIZE(saveidx)
        DO isave = 1, nsave
          IF ( (saveidx(isave) > 0) .AND. (saveidx(isave) <= nwl) ) THEN
            solution%base(saveidx(isave))%save = .TRUE.
          ELSE IF ( (saveidx(isave) < 0) .AND. (saveidx(isave) >= -nwl) ) THEN
            solution%base(nwl+1+saveidx(isave))%save = .TRUE.
          END IF

        END DO ! isave

        WRITE(*,*) '    Node solution%base%save is loaded.'
      END IF
    END IF

    CALL fson_get(json_root, "focal", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nfocal=fson_value_count(json_array)
      ALLOCATE(solution%focal(1:nfocal))
      DO ifocal = 1, nfocal
        ! Get the array item (this is an associative array)
        json_item => fson_value_get(json_array, ifocal)
        ! array in json file is indexed from 0,
        ! but fson_value_get function is indexed from 1

        CALL json2focal(json_item, solution%focal(ifocal))
      END DO! ifocal
      WRITE(*,*) '  Node solution%focal is loaded.'
    END IF

  END SUBROUTINE json2solution

  SUBROUTINE json2focal(json_root, focal)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(focal3d_type), INTENT(INOUT)         :: focal

    TYPE(fson_value), POINTER                 :: json_item
    REAL(KIND=dp)                 :: x1a, x1b, x2a, x2b, x3a, x3b
    REAL(KIND=dp), DIMENSION(2)   :: x1, x2, x3
    INTEGER                       :: n1, n2, n3

    CALL fson_get(json_root, "label", focal%label)
    CALL fson_get(json_root, "isrc", focal%isrc)
    CALL fson_get(json_root, "wl", focal%wl)

    CALL fson_get(json_root, "coordinate", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "coordinate", focal%grid%coord)
    ELSE
      focal%grid%coord = 'car'
    END IF

    CALL fson_get(json_root, "x1[1]", n1)
    CALL fson_get(json_root, "x2[1]", n2)
    CALL fson_get(json_root, "x3[1]", n3)

    CALL fson_get(json_root, "x1[2]", x1a)
    CALL fson_get(json_root, "x2[2]", x2a)
    CALL fson_get(json_root, "x3[2]", x3a)

    IF ( n1 .EQ. 1 ) THEN
      x1b = x1a
    ELSE
      CALL fson_get(json_root, "x1[3]", x1b)
    END IF

    IF ( n2 .EQ. 1 ) THEN
      x2b = x2a
    ELSE
      CALL fson_get(json_root, "x2[3]", x2b)
    END IF

    IF ( n3 .EQ. 1 ) THEN
      x3b = x3a
    ELSE
      CALL fson_get(json_root, "x3[3]", x3b)
    END IF

    x1 = (/x1a, x1b/)
    x2 = (/x2a, x2b/)
    x3 = (/x3a, x3b/)

    ALLOCATE( focal%nodes(1:n1*n2*n3) )

    CALL generate_nodes(n1, n2, n3, x1, x2, x3, focal%nodes, focal%grid%coord)

    focal%grid%n1 = n1
    focal%grid%n2 = n2
    focal%grid%n3 = n3
    focal%grid%x1 = x1
    focal%grid%x2 = x2
    focal%grid%x3 = x3

  END SUBROUTINE json2focal

  SUBROUTINE write_matrix1d(filename, vector)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)                :: filename
    COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector
    INTEGER                                     :: fid = 10, iovar, i

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1, SIZE(vector)
!      WRITE(fid, fmt_cmplx, ADVANCE='NO') &
!            REAL(vector(i)), AIMAG(vector(i)), 'i'
      WRITE(fid, fmt_cmplx, ADVANCE='NO') vector(i)
      WRITE(fid, '(/)', ADVANCE='NO')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_matrix1d

  SUBROUTINE write_matrix2d(filename, matrix)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)                :: filename
    COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN):: matrix
    INTEGER                                     :: fid = 10, iovar, i, j

    OPEN(fid, FILE=TRIM(filename), ACTION='WRITE', IOSTAT=iovar)
    IF(iovar>0) THEN
       WRITE(*,*) 'Could not open output file' // filename // '!'
       STOP
    END IF

    DO i=1,SIZE(matrix,1)
      DO j=1,SIZE(matrix,2)
!        WRITE(fid, fmt_cmplx, ADVANCE='NO') &
!              REAL(matrix(i,j)), AIMAG(matrix(i,j)), 'i'
        WRITE(fid, fmt_cmplx, ADVANCE='NO') matrix(i,j)
        WRITE(fid, '(4X)', ADVANCE='NO')
      END DO
      WRITE(fid, '(/)', ADVANCE='NO')
    END DO

    CLOSE(fid)
  END SUBROUTINE write_matrix2d

  SUBROUTINE write_field_h5(filename, grid, field, coord_set)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)      :: filename
    TYPE(grid_type), INTENT(IN)       :: grid
    TYPE(field_type), INTENT(IN), DIMENSION(:) :: field
    INTEGER(HID_T)                    :: file_id
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL  :: coord_set! ['car', 'cyl', 'sph']

    CHARACTER(LEN=3)                  :: coord
    REAL(KIND=dp), DIMENSION(grid%n1, grid%n2, grid%n3) :: e1_re, e1_im, e2_re, e2_im, e3_re, e3_im
    REAL(KIND=dp), DIMENSION(grid%n1) :: vx1
    REAL(KIND=dp), DIMENSION(grid%n2) :: vx2
    REAL(KIND=dp), DIMENSION(grid%n3) :: vx3

    IF (PRESENT(coord_set)) THEN
      coord = coord_set
    ELSE
      coord = 'car'
    END IF
    ! create and open h5 file, get the file_id
    CALL create_h5(filename, file_id)

    e1_re = RESHAPE(REAL(field(:)%e(1)), (/grid%n1, grid%n2, grid%n3/))
    e1_im = RESHAPE(AIMAG(field(:)%e(1)), (/grid%n1, grid%n2, grid%n3/))
    e2_re = RESHAPE(REAL(field(:)%e(2)), (/grid%n1, grid%n2, grid%n3/))
    e2_im = RESHAPE(AIMAG(field(:)%e(2)), (/grid%n1, grid%n2, grid%n3/))
    e3_re = RESHAPE(REAL(field(:)%e(3)), (/grid%n1, grid%n2, grid%n3/))
    e3_im = RESHAPE(AIMAG(field(:)%e(3)), (/grid%n1, grid%n2, grid%n3/))

    SELECT CASE (TRIM(coord))
      CASE ('car')
        CALL write_h5_d3(file_id, 'ex_re', e1_re)
        CALL write_h5_d3(file_id, 'ex_im', e1_im)
        CALL write_h5_d3(file_id, 'ey_re', e2_re)
        CALL write_h5_d3(file_id, 'ey_im', e2_im)
      CASE ('cyl')
        CALL write_h5_d3(file_id, 'erho_re', e1_re)
        CALL write_h5_d3(file_id, 'erho_im', e1_im)
        CALL write_h5_d3(file_id, 'ephi_re', e2_re)
        CALL write_h5_d3(file_id, 'ephi_im', e2_im)
    END SELECT
    CALL write_h5_d3(file_id, 'ez_re', e3_re)
    CALL write_h5_d3(file_id, 'ez_im', e3_im)

    vx1 = linspace(grid%x1(1), grid%x1(2), grid%n1)
    vx2 = linspace(grid%x2(1), grid%x2(2), grid%n2)
    vx3 = linspace(grid%x3(1), grid%x3(2), grid%n3)
    SELECT CASE ( grid%coord )
      CASE ('car')
        CALL write_h5_d1(file_id, 'x', vx1)
        CALL write_h5_d1(file_id, 'y', vx2)
        CALL write_h5_d1(file_id, 'z', vx3)
      CASE ('cyl')
        CALL write_h5_d1(file_id, 'rho', vx1)
        CALL write_h5_d1(file_id, 'theta', vx2)
        CALL write_h5_d1(file_id, 'z', vx3)
      CASE ('sph')
        CALL write_h5_d1(file_id, 'r', vx1)
        CALL write_h5_d1(file_id, 'theta', vx2)
        CALL write_h5_d1(file_id, 'phi', vx3)
    END SELECT

    ! close h5 file by file_id
    CALL close_h5(file_id)

    WRITE(*,*) '    Save field to "', TRIM(filename), '".'
  END SUBROUTINE write_field_h5

!  SUBROUTINE save(model)
!    TYPE(model_type), INTENT(INOUT)   :: model
!  END SUBROUTINE save

END MODULE io
