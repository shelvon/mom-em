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
  USE h5_wrapper

  IMPLICIT NONE

  INTERFACE write_matrix
    MODULE PROCEDURE write_matrix1d, write_matrix2d
  END INTERFACE

CONTAINS
  ! fson_value to model
  SUBROUTINE json2model(json_file, model)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN)       :: json_file
    TYPE(model_type), INTENT(INOUT)     :: model

    TYPE(fson_value), POINTER           :: json_root, json_array, json_item

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
      CALL json2solution(json_item, model%solution)
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
    INTEGER                                   :: ndom, idom, nquad, iquad, dim

    ! fson_value to geom%mesh
    CALL fson_get(json_root,"mesh", json_item)
    CALL json2mesh(json_item, geom%mesh)

    ! set global quadrature rule
    geom%quad(2)%dim = 2
    geom%quad(2)%rule = 'tri_gl13'  ! by default
    geom%quad(3)%dim = 3
    geom%quad(3)%rule = 'tetra_gl4' ! by default
    ! load global qd rules if they exist
    CALL fson_get(json_root, "quad", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nquad = fson_value_count(json_array)
      DO iquad = 1, nquad
        ! Get the array item (this is an associative array)
        json_item => fson_value_get(json_array, iquad)
        ! array in json file is indexed from 0,
        ! but fson_value_get function is indexed from 1

        CALL fson_get(json_item, "dim", dim)
        IF ( (dim == 2) .OR. (dim == 3) ) THEN
          CALL json2quad( json_item, geom%quad(dim) )
        END IF
      END DO
    END IF

    ! fson_value to geom%domain
    CALL fson_get(json_root, "domain", json_array)
    ndom = fson_value_count(json_array)
    ! index starts from -1, i.e.
    ! -1: all domains, i.e. the entire system;
    ! 0: air/vacuum (or the surrounding media);
    ! 1: the first scatter;
    ! etc.
    ALLOCATE(geom%domain(-1:ndom-1))

    ! copy global qd rule to domain(-1), i.e. the entire system
    ALLOCATE(geom%domain(-1)%quad(2:3))
    geom%domain(-1)%quad(2) = geom%quad(2)
    geom%domain(-1)%quad(3) = geom%quad(3)

    DO idom = 1, ndom
      ! Get the array item (this is an associative array)
      json_item => fson_value_get(json_array, idom)
      ! array in json file is indexed from 0,
      ! but fson_value_get function is indexed from 1

      CALL json2domain(json_item, geom%domain(idom-1))
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
  SUBROUTINE json2domain(json_root, domain)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(domain_type), INTENT(INOUT)          :: domain

    TYPE(fson_value), POINTER                 :: json_item, json_array
    INTEGER                                   :: nsurface, isurface,  &
                                                 nvolume, ivolume,    &
                                                 nquad, iquad

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

    ! get quadrature rules for surface/volume, specified by surface/volume ids
    CALL fson_get(json_root, "quad", json_array)
    IF ( ASSOCIATED(json_array) ) THEN
      nquad = fson_value_count(json_array)
      ALLOCATE(domain%quad(1:nquad))
      DO iquad = 1, nquad
        ! Get the array item (this is an associative array)
        json_item => fson_value_get(json_array, iquad)
        ! array in json file is indexed from 0,
        ! but fson_value_get function is indexed from 1

        CALL json2quad( json_item, domain%quad(iquad) )
      END DO! imedia
    END IF

  END SUBROUTINE json2domain

  SUBROUTINE json2quad(json_root, quad)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(quad_type), INTENT(INOUT)            :: quad

    CALL fson_get(json_root, "dim", quad%dim)
    CALL fson_get(json_root, "id", quad%id)
    CALL fson_get(json_root, "rule", quad%rule)

  END SUBROUTINE json2quad
  
  ! json_value to physics
  SUBROUTINE json2physics(json_root, physics)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(physics_type), INTENT(INOUT)         :: physics

    TYPE(fson_value), POINTER                 :: json_item, json_array
    INTEGER                                   :: nmedia, imedia,  &
                                                 nsource, isource,&
                                                 npupil, ipupil, nga

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

    ! json_value to physics%pupil
    CALL fson_get(json_root, "pupil",json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL json2pupil(json_item, physics%pupil)
      WRITE(*,*) '  Node physics%pupil is loaded.'
    END IF
!    CALL fson_get(json_root, "pupil",json_array)
!    IF ( associated(json_array) ) THEN
!      npupil=fson_value_count(json_array)
!      ALLOCATE(physics%pupil(1:npupil))
!      DO ipupil = 1, npupil
!        json_item => fson_value_get(json_array, ipupil)
!        CALL json2pupil(json_item, physics%pupil(ipupil))
!      END DO! ipupil
!    END IF
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
    ! get ri, i.e. refractive index
    CALL json2property(json_root, "ri", media%ri)
    ! get epsilon
    CALL json2property(json_root, "epsilon", media%eps)
    ! get mu
    CALL json2property(json_root, "mu", media%mu)
  END SUBROUTINE json2media

  ! json_value to source
  SUBROUTINE json2source(json_root, source)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    ! TYPE(srcdata), INTENT(INOUT)              :: source
    TYPE(source_type), INTENT(INOUT)          :: source

    TYPE(fson_value), POINTER                 :: json_item, json_array
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: pos
    CHARACTER(LEN=16)                         :: bessel_input

    ! get media name
    CALL fson_get(json_root, "type", source%type)
    CALL fson_get(json_root, "norm", source%norm)
    IF (source%type == 'pw') THEN
      ! get incident angle, required
      CALL fson_get(json_root, "theta", source%pw%theta)
      CALL fson_get(json_root, "phi", source%pw%phi)
      ! get relative phase difference
      CALL fson_get(json_root, "phase", source%pw%phase)
      IF ( ABS(source%pw%phase) .EQ. 90) THEN
        ! circularly polarized light
        source%pw%A = 1.0_dp/SQRT(2.0_dp)
        source%pw%B = 1.0_dp/SQRT(2.0_dp)
        WRITE(*,*) 'circularly polarized plane wave!'
      ELSE
        ! get either psi or A, B values
        CALL fson_get(json_root, "psi", json_item)
        IF (associated(json_item)) THEN
          ! polarization direction reads from psi angle!
          CALL fson_get(json_root, "psi", source%pw%psi)
        ELSE
          CALL fson_get(json_root, "A", json_item)
          IF (associated(json_item)) THEN
            ! get electric strength of two orthogonal components
            CALL fson_get(json_root, "A", source%pw%A)
            CALL fson_get(json_root, "B", source%pw%B)
          ELSE
            WRITE(*,*) 'error input parameters, either psi or A/B should be defined!'
          END IF
        END IF
      END IF

    ELSE IF (source%type == 'focus') THEN
      CALL fson_get(json_root, "focustype", source%focus%type)
      CALL fson_get(json_root, "focal", source%focus%focal)
      CALL fson_get(json_root, "waist", source%focus%waist)
      CALL fson_get(json_root, "na", source%focus%na)

      CALL fson_get(json_root, "pos", pos)
      source%focus%pos = pos

      IF (source%focus%type == 'bessel') THEN
        CALL fson_get(json_root, "bessel.theta", source%focus%bessel%theta)
        CALL fson_get(json_root, "bessel.delta", source%focus%bessel%delta)
        CALL fson_get(json_root, "bessel.j", source%focus%bessel%j)
        CALL fson_get(json_root, "bessel.input", bessel_input)
        CALL fson_get(json_root, "bessel.nbasis", source%focus%bessel%nbasis)
        CALL fson_get(json_root, "bessel.polyName", source%focus%bessel%polyName)
        IF ( LEN_TRIM(bessel_input) .NE. 3 ) THEN
          WRITE(*,*) 'Invalid input beam type for bessel beam calculation.'
          STOP ERR0
        ELSE
          source%focus%bessel%input = TRIM(bessel_input)
        END IF
      END IF
    END IF
  END SUBROUTINE json2source

  SUBROUTINE json2pupil(json_root, pupil)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(pupil_type), INTENT(INOUT)           :: pupil

    TYPE(fson_value), POINTER                 :: json_item
    CALL fson_get(json_root, "type", pupil%type)
    CALL fson_get(json_root, "aperture", json_item)
    CALL json2aperture(json_item, pupil%aperture)
    CALL fson_get(json_root, "phase", json_item)
    CALL json2phase(json_item, pupil%phase)
  END SUBROUTINE json2pupil

  SUBROUTINE json2aperture(json_root, aperture)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(aperture_type), INTENT(INOUT)        :: aperture

    CALL fson_get(json_root, "type", aperture%type)
    IF ( aperture%type == 'circ' ) THEN
      CALL fson_get(json_root, "circ.r", aperture%circ%r)
    ELSE
      WRITE(*,*) 'Aperture type of ', aperture%type, ' is not supported!'
      RETURN
    END IF
  END SUBROUTINE json2aperture

  SUBROUTINE json2phase(json_root, phase)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(phase_type), INTENT(INOUT)           :: phase

    CALL fson_get(json_root, "type", phase%type)
    IF ( phase%type == 'bessel' ) THEN
      CALL fson_get(json_root, "bessel.j", phase%bessel%j)
      CALL fson_get(json_root, "bessel.theta", phase%bessel%theta)
    END IF
  END SUBROUTINE json2phase

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
      CALL fson_get(json_item, 'method', property%method)
      IF ( property%method=='value' ) THEN
        ! non-dispersive, i.e. constant value
        ALLOCATE(property%re(1))
        ALLOCATE(property%im(1))
        ALLOCATE(property%z(1))
        property%re(1) = 0.0_dp
        property%im(1) = 0.0_dp
        property%z(1) = 0.0_dp
        CALL fson_get(json_item, 'real', json_item_sub)
        IF (associated(json_item_sub)) THEN
          CALL fson_get(json_item, 'real', property%re(1))
        END IF

        CALL fson_get(json_item, 'imag', json_item_sub)
        IF (associated(json_item_sub)) THEN
          CALL fson_get(json_item, 'imag', property%im(1))
        END IF
        property%z = CMPLX(property%re, property%im)

      ELSE IF ( property%method=='table' ) THEN
        CALL fson_get(json_item, 'ref_file', property%ref_file)
        ! CALL fson_get(json_item, 'ref_file', json_item_sub)
        ! IF (associated(json_item_sub)) THEN
        !   CALL fson_get(json_item, 'ref_file', property%ref_file)
        ! END IF
      ELSE IF (property%method=='model') THEN
        CALL fson_get(json_item, 'model', property%model)
      END IF
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
      WRITE(*,*) '  Loading method for node simulation%frequency is set to &
                 "', TRIM(freq_method),'".'
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
      nfreq = freq_value(1)
      ALLOCATE(simulation%omega(1:nfreq))
      simulation%omega = linspace(freq_value(2), freq_value(3), nfreq)
    ELSE IF ( TRIM(freq_method) == 'zrange' ) THEN
      nfreqx = freq_value(1)
      nfreqy = freq_value(2)
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
    INTEGER                                   :: nwl, iwl, nwlx, nwly

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
      ELSE IF ( TRIM(wl_method) == 'range' ) THEN
        nwl = wl_value(1)
        ALLOCATE(simulation%wl(1:nwl))
        simulation%wl = linspace(wl_value(2), wl_value(3), nwl)

        ALLOCATE(simulation%zwl(1:nwl))
        simulation%zwl = zlinspace(CMPLX(wl_value(2), 0.0_dp, KIND=dp),&
                                   CMPLX(wl_value(3), 0.0_dp, KIND=dp),&
                                   nwl, 1)
        WRITE(*,*) simulation%zwl
      ELSE IF ( TRIM(wl_method) == 'zrange' ) THEN
        nwlx = wl_value(1)
        nwly = wl_value(2)
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

  SUBROUTINE json2solution(json_root, solution)
    IMPLICIT NONE
    TYPE(fson_value), POINTER, INTENT(IN)     :: json_root
    TYPE(solution_type), INTENT(INOUT)        :: solution

    TYPE(fson_value), POINTER                 :: json_item, json_item_sub, json_array
    INTEGER                                   :: nfocal, ifocal

    CALL fson_get(json_root, "base", json_item)
    IF ( ASSOCIATED(json_item) ) THEN
      CALL fson_get(json_root, "base.save", json_item_sub)
      IF ( ASSOCIATED(json_item) ) THEN
        CALL fson_get(json_root, "base.save", solution%base%save)
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

    CALL fson_get(json_root, "label", focal%label)
    CALL fson_get(json_root, "isrc", focal%isrc)
    CALL fson_get(json_root, "wl", focal%wl)
    CALL fson_get(json_root, "x[1]", focal%grid%nx)
    CALL fson_get(json_root, "y[1]", focal%grid%ny)
    CALL fson_get(json_root, "z[1]", focal%grid%nz)

    CALL fson_get(json_root, "x[2]", focal%grid%xa)
    CALL fson_get(json_root, "y[2]", focal%grid%ya)
    CALL fson_get(json_root, "z[2]", focal%grid%za)

    IF ( focal%grid%nx .EQ. 1 ) THEN
      focal%grid%xb = focal%grid%xa
    ELSE
      CALL fson_get(json_root, "x[3]", focal%grid%xb)
    END IF

    IF ( focal%grid%ny .EQ. 1 ) THEN
      focal%grid%yb = focal%grid%ya
    ELSE
      CALL fson_get(json_root, "y[3]", focal%grid%yb)
    END IF

    IF ( focal%grid%nz .EQ. 1 ) THEN
      focal%grid%zb = focal%grid%za
    ELSE
      CALL fson_get(json_root, "z[3]", focal%grid%zb)
    END IF
    ALLOCATE( focal%grid%x(1:focal%grid%nx), focal%grid%y(1:focal%grid%ny), &
              focal%grid%z(1:focal%grid%nz))
    focal%grid%x = linspace(focal%grid%xa, focal%grid%xb, focal%grid%nx)
    focal%grid%y = linspace(focal%grid%ya, focal%grid%yb, focal%grid%ny)
    focal%grid%z = linspace(focal%grid%za, focal%grid%zb, focal%grid%nz)

    CALL meshgrid( focal%grid%x, focal%grid%y, focal%grid%z, &
                   focal%grid%gridx, focal%grid%gridy, focal%grid%gridz)

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

  SUBROUTINE write_field_h5(filename, grid, field)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)      :: filename
    TYPE(grid3d_type), INTENT(IN)     :: grid
    TYPE(field3d_type), INTENT(IN)    :: field

    CALL create_h5(filename)
    CALL write_h5_d3(filename, 'ex_re', REAL(field%ex))
    CALL write_h5_d3(filename, 'ex_im', AIMAG(field%ex))
    CALL write_h5_d3(filename, 'ey_re', REAL(field%ey))
    CALL write_h5_d3(filename, 'ey_im', AIMAG(field%ey))
    CALL write_h5_d3(filename, 'ez_re', REAL(field%ez))
    CALL write_h5_d3(filename, 'ez_im', AIMAG(field%ez))

    CALL write_h5_d1(filename, 'x', grid%x)
    CALL write_h5_d1(filename, 'y', grid%y)
    CALL write_h5_d1(filename, 'z', grid%z)

    ! For test purpose, save curve interpolation coefficients, cj
    ! CALL write_h5_z5(filename, 'cj', grid%cj) ! cj: double complex, rank 5
    CALL write_h5_z5(filename, 'cjs0c0', grid%cjs0c0) ! cj: double complex, rank 5
    CALL write_h5_z5(filename, 'cjs0c1', grid%cjs0c1) ! cj: double complex, rank 5
    CALL write_h5_z5(filename, 'cjs1c0', grid%cjs1c0) ! cj: double complex, rank 5

    WRITE(*,*) 'Save field to "', TRIM(filename), '".'
  END SUBROUTINE write_field_h5

!  SUBROUTINE save(model)
!    TYPE(model_type), INTENT(INOUT)   :: model
!  END SUBROUTINE save

END MODULE io
