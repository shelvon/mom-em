! MODULE: data
! AUTHOR: Xiaorun (Shelvon) ZANG
! DESCRIPTION:
! The module for the definition and related operations of
! the essenstial data structures of the model.
MODULE data
  USE constants

  IMPLICIT NONE

  ! electric and magnetic fields in 3d
  TYPE field3d_type
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)   :: ex, ey, ez
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)   :: hx, hy, hz
  END TYPE field3d_type

  ! data structure in the focal region (3d)
  TYPE grid3d_type
    INTEGER                                           :: nx, ny, nz
    REAL(KIND=dp)                                     :: xa, xb, ya, yb, za, zb
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: x, y, z
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)      :: gridx
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)      :: gridy
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:)      :: gridz
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: cj
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: cjs0c0
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: cjs0c1
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: cjs1c0
  END TYPE grid3d_type

  TYPE focal3d_type
    CHARACTER(LEN=256)                                :: label
    INTEGER                                           :: isrc
    ! REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)          :: wl
    REAL(KIND=dp)                                     :: wl
    TYPE(grid3d_type)                                 :: grid
    TYPE(field3d_type)                                :: field
    ! TYPE(field3d_type), ALLOCATABLE, DIMENSION(:)     :: field ! field(nwl)
  END TYPE focal3d_type

  !----data type for greenprd----
  TYPE prdcoef
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: samples
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: samplesz
     REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: rho ! Jouni's theis, Eq. (3.28)
     REAL (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: kt ! Jouni's theis, Eq. (3.29)
     REAL (KIND=dp) :: E, k0x, k0y, wl, range ! E, the splitting parameter
     INTEGER :: ic, jc, kc, n, m, np, npz
     COMPLEX (KIND=dp) :: ri, ri0, k
  END TYPE prdcoef

  ! Stores data for pre-computed periodic GF. This data is passed down to
  ! routines that calculate potential integrals.
  TYPE prdnfo
     CHARACTER (LEN=256) :: filename
     INTEGER :: type
     INTEGER :: nwl, cwl
     REAL (KIND=dp) :: dx, dy, dz, phi, phasex1, phasex2, phasey
     REAL (KIND=dp) :: cp, sp, pwtheta, pwphi
     LOGICAL :: oblique ! oblique incidence
     TYPE(prdcoef), DIMENSION(:), ALLOCATABLE :: coef
  END TYPE prdnfo
  !----data type for greenprd----

  !----data type for symmetry----
  TYPE group_action
     ! Action for electric field. Action for H-field is ef*detj.
     ! There can be multiple sub-problems with different field actions.
     COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: ef

     ! Jacobian the group action for points in R^3.
     REAL (KIND=dp), DIMENSION(3,3) :: j

     ! Determinant of the jacobian.
     REAL (KIND=dp) :: detj

     ! Identifier of a special action.
     INTEGER :: id

     ! Bit flags for group element generators.
     INTEGER :: genbits
  END type group_action
  !----data type for symmetry----

  !----data type for pupil----
  ! circ function for apodization in raidal direction
  TYPE circ_type
    REAL (KIND=dp)      :: r=1  ! radius of circ function, normalized to
                                ! the waist radius of the incident beam
  END TYPE circ_type

  ! Intervals where function is equal to 1, which is also
  ! where the phase modulation applies.
  TYPE rect_type
    REAL (KIND=dp), DIMENSION(1:2)              :: delta! phase shift [\pi]
    INTEGER                                     :: intervals_n
    ! REAL (KIND=dp), DIMENSION(2) :: intervals
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: intervals
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: dj
    ! with (1,:) stores aj and (2,:) bj
    ! aj for cos(j\phi) and bj for sin(j\phi)
  END TYPE rect_type

  ! topological charge of petal beam
  TYPE petal_type
    INTEGER :: l  ! l is the integer subindex in
                  ! Petal_{p,l} = LG_{p,l} + LG_{p,-l}
  END TYPE petal_type

  ! petal beam + rect phase mask
  TYPE petal_rect_type
    INTEGER :: l  ! l is the integer subindex in
                  ! Petal_{p,l} = LG_{p,l} + LG_{p,-l}
    INTEGER :: ln ! ln-petal is retained and the other petals are blocked
    INTEGER                                     :: intervals_n
    REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: intervals
  END TYPE petal_rect_type

  TYPE bessel_type
    CHARACTER(LEN=3)  :: input ! input beam: rad/azi
    INTEGER           :: j = 0 ! jth-order of Bessel function
    REAL(KIND=dp)     :: theta ! the diverging angle (of axicon)
                               ! going into J_j(k \rho sin (\theta))
    REAL(KIND=dp)     :: delta ! arising from a finite thickness of the ring
    INTEGER           :: nbasis ! number of basis functions in bessel_transform
    CHARACTER(LEN=12) :: polyName ! name of basis function
  END TYPE bessel_type

  TYPE vortex_type
    INTEGER :: charge ! 'topological charge x 2\pi' phase change
                      ! over one rotation in azimuthal angle
  END TYPE vortex_type

  TYPE aperture_type
    CHARACTER (LEN=32)  :: type ! name of the aperture function, chosen below
    TYPE(circ_type)     :: circ
  END TYPE aperture_type

  TYPE phase_type
    CHARACTER(LEN=32)       :: type   ! name of the phase function below
    CHARACTER(LEN=8)        :: input  ! input beam: hg00, hg01, hg10, rad/azi
    INTEGER                 :: nMax   ! series expansion up to nMax order
    TYPE(rect_type)         :: rect   ! binary type phase modulation
    TYPE(vortex_type)       :: vortex ! spiral phase plate
    TYPE(petal_type)        :: petal  ! petal beam phase term
    TYPE(bessel_type)       :: bessel ! bessel beam
    TYPE(petal_rect_type)   :: petal_rect ! petal beam phase + rect phase
    REAL(KIND=dp), DIMENSION(1:2) :: thetaRot ! rotation angles[deg], x-, y-pols
    COMPLEX(KIND=dp), DIMENSION(1:2) :: vJones ! Jones vector, x-, y-pols
  END TYPE phase_type

  TYPE pupil_type
    ! either analytical or numerical
    CHARACTER (LEN=10)  :: type = 'analytical'

    ! aperture function
    TYPE(aperture_type) :: aperture

    ! phase function
    TYPE(phase_type)    :: phase

    ! The complex-valued Fourier expansion coefficient, works only for
    ! pupil function that is independent of radial variable.
    ! cn = 1/(2*pi)*\Int_{-pi}^{pi} e^{i\Phi(\phi)}e^{-i n \phi} d\phi
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: cn
  END TYPE pupil_type
  !----data type for pupil----


  !----data type for source----
  !
  ! Plane wave in general form:
  ! E = (E1*e1 + E2*e2)*exp(i*k*r)*exp(-i*omega*t),
  ! where direction of wave propagation vector is k;
  ! e1, e2 are two orthogonal basis unit vector and
  ! (e1, e2, k) is a right-handed orthogonal triad of unit vector;
  ! e1 is along k x z direction and e2 is in the plane of k, z.
  ! Note that even in the special cases where k//z where phi angle is irrelevant,
  ! angle phi is still used to define e1, e2.
  ! E1 = A*exp(i*delta_phi1), E2 = B*exp(i*delta_phi2);
  ! delta_phi2 - delta_phi1 = phase, relative phase difference between E1 and E2.
  ! The eccentricity and orientation of the ellipse depend on the
  ! phase difference and the amplitude ratio B/A.
  !
  ! --three cases:
  ! a) phase = pi*n : linearly polarized
  !    then, psi or A/B value defines the polarization direction
  ! b) phase = pi/2*n & A = B: circularly porlarized
  ! c) arbitrary phase and A, B values: elliptically polarized
  !
  ! --k wave vector
  ! theta: angle between +z axis and wave vector k;
  ! phi: angle between +x axis and k projected on to xy plane
  ! when theta = 0 or 180, then phi value does not matter
  ! kz = k*cos(theta), kx = k*sin(theta)*cos(phi), ky = k*sin(theta)*sin(phi)
  TYPE pw_type
    REAL(KIND=dp)     :: theta=180, phi=0! define propagation direction
    REAL(KIND=dp)     :: phase=0
    REAL(KIND=dp)     :: psi=0! angle from e1 to e2, i.e. B/A = tan(psi)
    ! psi: polarization direction with respect to +x axis for linearly polarized case
    REAL(KIND=dp)     :: A=1, B=0! electric field strength, to be normalized
  END TYPE pw_type

  ! focus beams are catagorized as:
  ! -- hg00, hg01, hg10
  ! -- linear, rad, azi
  ! -- bessel, petal beams
  ! -- others, to be done
  TYPE focus_type
    CHARACTER(LEN=32)           :: type
    REAL(KIND=dp)               :: focal, waist, na
    REAL(KIND=dp), DIMENSION(3) :: pos
    TYPE(bessel_type)           :: bessel
    TYPE(petal_type)            :: petal
  END TYPE focus_type

  TYPE srcdata
     ! src_type: Identifier for the excitation source.
     INTEGER :: type

     ! Angles for a plane-wave excitation.
     REAL (KIND=dp) :: theta, phi, psi

     ! Phase factor between Ex and Ey, where (Ex,Ey,k) form a right-handed system.
     COMPLEX (KIND=dp) :: phase

     ! Focused beam excitation parameters.
     REAL (KIND=dp) :: focal, waist, napr

     ! Normalize beam fields to maximum at focal plane?
     LOGICAL :: nfocus

     ! Excitation source position for, e.g., beams and a dipole.
     REAL (KIND=dp), DIMENSION(3) :: pos

     ! In which plane of the source locate
     CHARACTER (LEN=256) :: sxyz

     ! Dipole moment of a dipole source.
     COMPLEX (KIND=dp), DIMENSION(3) :: dmom

  END TYPE srcdata
  !----data type for source----

  !----data type for nlsurf----
  ! Description of the surface nonlinearity of material.
  TYPE medium_nls
     COMPLEX (KIND=dp) :: chi2_nnn
     COMPLEX (KIND=dp) :: chi2_ntt
     COMPLEX (KIND=dp) :: chi2_ttn
     INTEGER :: nsurf_ids
     INTEGER, DIMENSION(:), ALLOCATABLE :: surf_ids
  END TYPE medium_nls
  !----data type for nlsurf----

  !----data type for nlbulk----
  ! Properties of material bulk nonlinearity.
  TYPE medium_nlb
     COMPLEX (KIND=dp) :: delta_prime
     COMPLEX (KIND=dp) :: gamma
     !COMPLEX (KIND=dp) :: chi2zzz ! deprecated

     ! Maps (Ex**2,Ey**2,Ez**2,EyEz,ExEz,ExEy) to (Px,Py,Pz) in crystal frame.
     ! chi2 is in D-matrix form with Kleinman's symmetry
     ! d_il=(d11,d12,...,d16;
     !       d21,d22,...,d26;
     !       d31,d32,...,d36)
     ! ijk  11  22  33  23,32   31,13   12,21
     !  l    1   2   3      4       5       6
     COMPLEX (KIND=dp), DIMENSION(3,6) :: chi2

     ! Maps (ExExEx,...) to (Px,Py,Pz) in crystal frame.
     ! chi3 is in D-matrix form with Kleinman's symmetry
     ! d_il=(d11,d12,...,d110;
     !       d21,d22,...,d210;
     !       d31,d32,...,d310)
     ! ijk  111 222 333 233,323,332 322,232,223 331,313,133 113,131,311 122,212,221 112,121,211 123,132,213,231,312,321
     !  l    1   2   3       4           5           6           7           8           9                 10
     COMPLEX (KIND=dp), DIMENSION(3,10) :: chi3

     ! T maps from crystal frame to lab frame. invT = inverse(T).
     REAL (KIND=dp), DIMENSION(3,3) :: T, invT
  END type medium_nlb
  !----data type for nlbulk----

  !----data type for quad----
  TYPE quad_data
     CHARACTER (LEN=256) :: description
     INTEGER :: num_nodes
     REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: weights
     REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: nodes
  END type quad_data
  !----data type for quad----

  !----data type for nfields----
  TYPE nfield_plane
     REAL (KIND=dp), DIMENSION(3) :: origin, v1, v2
     REAL (KIND=dp) :: d1, d2
     INTEGER :: n1, n2
  END TYPE nfield_plane
  !----data type for nfields----

  !----data type for mesh_container----
  TYPE node
    REAL (KIND=dp), DIMENSION(3)    :: p
    ! INTEGER, DIMENSION(:), POINTER :: face_indices, node_indices
    INTEGER                         :: parent_index
    ! INTEGER                         :: bnd, nbnd
  END TYPE node

  TYPE face
    REAL (KIND=dp), DIMENSION(3) :: n, cp
    REAL (KIND=dp), DIMENSION(3,3) :: s, m
    INTEGER, DIMENSION(3) :: node_indices
    INTEGER, DIMENSION(3) :: edge_indices
    REAL (KIND=dp) :: area, pd
    INTEGER :: id
    INTEGER :: parent_index
  END TYPE face

  TYPE line
    INTEGER, DIMENSION(2) :: node_indices
    INTEGER               :: id
  END TYPE line

  TYPE edge
    INTEGER, DIMENSION(2) :: node_indices, bnode_indices, face_indices
    REAL (KIND=dp)        :: length
    ! REAL (KIND=dp), DIMENSION(2) :: rwgDiv
    INTEGER               :: bnd
    INTEGER               :: parent_index
    INTEGER               :: couple_index
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: child_indices ! (:,1)=submesh, (:,2)=local edge
  END TYPE edge

  ! Interior face.
  TYPE solid_face
    INTEGER, DIMENSION(3) :: node_indices
    INTEGER, DIMENSION(2) :: solid_indices, bnode_indices
    INTEGER               :: face_index ! -1 if not a boundary
    REAL (KIND=dp)        :: area
  END TYPE solid_face

  TYPE solid
    INTEGER, DIMENSION(4) :: node_indices
    INTEGER, DIMENSION(4) :: solid_face_indices
    REAL (KIND=dp)        :: volume
    INTEGER               :: id
  END TYPE solid

  TYPE mesh_container
    TYPE(node), DIMENSION(:), ALLOCATABLE        :: nodes
    TYPE(face), DIMENSION(:), ALLOCATABLE        :: faces
    TYPE(line), DIMENSION(:), ALLOCATABLE        :: lines
    TYPE(edge), DIMENSION(:), ALLOCATABLE        :: edges
    TYPE(solid), DIMENSION(:), ALLOCATABLE       :: solids
    TYPE(solid_face), DIMENSION(:), ALLOCATABLE  :: solid_faces
    INTEGER        :: nnodes, nfaces, nlines, nedges, nsolids, nsolid_faces
    REAL (KIND=dp) :: avelen
  END TYPE mesh_container
  !----data type for mesh_container----

  !----data type for batch----
  TYPE medium_prop
     COMPLEX (KIND=dp) :: ri
     COMPLEX (KIND=dp) :: shri
     COMPLEX (KIND=dp) :: thri
     TYPE(medium_nls) :: nls
     TYPE(medium_nlb) :: nlb
  END type medium_prop

  TYPE medium
     INTEGER :: type

     ! Range of wavelengths.
     TYPE(medium_prop), DIMENSION(:), ALLOCATABLE :: prop
  END type medium

  ! Medium types.
  INTEGER, PARAMETER :: mtype_linear = 1,&
       mtype_nls = 2,&
       mtype_nlb_nonlocal = 3,&
       mtype_nlb_dipole = 4,&
       mtype_thg_bulk = 31

  ! Solution data.
  TYPE solution
     ! Solution vector. First dimension denotes basis coefficients of (J,M).
     ! Second dimension denotes sub-problems related to group representations.
     ! Third dimension denotes excitation source.
     COMPLEX (KIND=dp), DIMENSION(:,:,:), ALLOCATABLE :: x, nlx

     ! Basis coefficients for jumps in M and J due to surface sources.
     ! 1st dim: coefficients of M and J jump expansions.
     ! 2nd dim: domain, 3rd dim: group representation, 4th dim: excitation source.
     COMPLEX (KIND=dp), DIMENSION(:,:,:,:), ALLOCATABLE :: src_coef

     ! Eigenvectors of a spectral problem. First dimension denotes basis coefficients.
     ! Second dimension denotes eigenvalue index.
     COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
     COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: adjeigvec

     ! Eigenvalues of a spectral problem.
     COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval

     ! Wavelength.
     REAL (KIND=dp) :: wl
  END TYPE solution

  ! Data for computational domain or a sub-domain.
  TYPE domain
     ! Index of precomputed Green's function data. Index -1 denotes non-periodic GF.
     INTEGER :: gf_index

     ! Index for medium associated with the domain.
     INTEGER :: medium_index

     ! Mesh of the domain boundary. May contain also a volume mesh.
     TYPE(mesh_container) :: mesh
  END TYPE domain

  ! Data for a computation batch, which way involve a set of wavelengths.
  TYPE batch
     ! Name of the problem (used for output file naming).
     CHARACTER (LEN=256) :: name

     ! Name of the mesh file.
     CHARACTER (LEN=256) :: mesh_file

     ! Mesh data (may contain multiple domains).
     TYPE(mesh_container) :: mesh

     ! Scale factor for the mesh.
     REAL (KIND=dp) :: scale

     ! nwl: Number of wavelengths.
     INTEGER :: nwl

     ! Pupil function for focused beam
     TYPE (pupil_type)    :: pupil

     ! data of focused beam
     ! TYPE (grid3d_type), DIMENSION(:), ALLOCATABLE :: focal ! multiple wavelengths if needed
     TYPE (grid3d_type)  :: focal

     ! Source data.
     TYPE(srcdata), DIMENSION(:), ALLOCATABLE :: src

     ! Solution data for each wavelength.
     TYPE(solution), DIMENSION(:), ALLOCATABLE :: sols

     ! Periodic Green's function data.
     TYPE(prdnfo), DIMENSION(:), POINTER :: prd

     ! Array of group actions.
     TYPE(group_action), DIMENSION(:), ALLOCATABLE :: ga

     ! Domains of the problem.
     TYPE(domain), DIMENSION(:), ALLOCATABLE :: domains

     ! Media.
     TYPE(medium), DIMENSION(:), ALLOCATABLE :: media

     ! Quadrature data.
     TYPE(quad_data) :: qd_tri
     TYPE(quad_data) :: qd_tetra
  END TYPE batch
  !----data type for batch----

  ! model->geom->mesh
  TYPE mesh_type
    CHARACTER(LEN=32)     :: type ! mesh file type, gmsh/
    CHARACTER(LEN=256)    :: file ! mesh filename.
    REAL (KIND=dp)        :: scale ! Scale factor for the mesh.
    ! Mesh data for the entire system (may contain multiple domains).
    ! A submesh under domain node stores the mesh for that specific domain.
    TYPE(mesh_container)  :: elements
  END TYPE mesh_type

  TYPE quad_type
    INTEGER               :: dim, id
    CHARACTER(LEN=16)     :: rule
    TYPE(quad_data)       :: qd
  END TYPE quad_type

  ! model->geom->domain
  TYPE domain_type
    INTEGER                             :: media    ! media index
    INTEGER, DIMENSION(:), ALLOCATABLE  :: surface  ! surface indices
    INTEGER, DIMENSION(:), ALLOCATABLE  :: volume   ! volume indices
    ! The submesh of the domain. May contain also a volume mesh.
    TYPE(mesh_container)                :: elements
    ! Surface/volume specific quadrature data
    ! In principle, possible to set separate quadrature rule for each surface.
    ! In principle, possible to set separate quadrature rule for each volume.
    TYPE(quad_type), DIMENSION(:), ALLOCATABLE  :: quad
  END TYPE domain_type

  ! model->geom
  ! store information for discretization, such as mesh and quadrature rules
  TYPE geom_type
    TYPE(mesh_type)                               :: mesh   ! loaded mesh
    TYPE(domain_type), DIMENSION(:), ALLOCATABLE  :: domain ! contains submesh
    ! Global quadrature data.
    ! Triangular surface quadrature rule for surface mesh elements.
    ! Tetrahedral volume quadrature rule for volume mesh elements.
    ! Use global quadrature rule except for the surface/volume specific ones.
    TYPE(quad_type), DIMENSION(2:3)               :: quad ! qd_tri, qd_tetra
  END TYPE geom_type

  ! model->physics->group
  TYPE group_type
    CHARACTER(LEN=3) , DIMENSION(:), ALLOCATABLE  :: name
    TYPE(group_action), DIMENSION(:), ALLOCATABLE :: action
  END TYPE group_type

  ! model->physics->media->material_property
  TYPE material_property
    LOGICAL                                       :: active
    CHARACTER(LEN=32)                             :: method
    CHARACTER(LEN=32)                             :: ref_file! for "[ref_file].ref" file
    CHARACTER(LEN=32)                             :: model
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE      :: re
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE      :: im
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: z
  END TYPE material_property

  ! model->physics->media
  TYPE media_type
    CHARACTER(LEN=32)             :: name
    TYPE(material_property)       :: eps ! relative permittivity
    TYPE(material_property)       :: mu ! relative permeability
    TYPE(material_property)       :: ri ! refractive index, i.e. ri = n + i*k
  END TYPE media_type

  ! model->physics->source
  TYPE source_type
    CHARACTER(LEN=32)             :: type
    LOGICAL                       :: norm ! normalize beam to maximum in focal plane
    TYPE(pw_type)                 :: pw
    TYPE(focus_type)              :: focus
  END TYPE source_type

  ! model->physics
  TYPE physics_type
    TYPE(group_type)                              :: group
    TYPE(media_type), DIMENSION(:), ALLOCATABLE   :: media
    TYPE(source_type), DIMENSION(:), ALLOCATABLE  :: source
    ! TYPE(pupil_type), DIMENSION(:), ALLOCATABLE   :: pupil
    TYPE(pupil_type)                              :: pupil
  END TYPE physics_type

  ! model->simulation->solver
  TYPE solver_type
    CHARACTER(LEN=32)                             :: name
    LOGICAL                                       :: parallel
    LOGICAL                                       :: memory
    CHARACTER(LEN=12)                             :: poles

  END TYPE solver_type
  ! model->simulation
  TYPE simulation_type
    TYPE(solver_type)                             :: solver
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)      :: wl
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)      :: omega
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)   :: zwl ! z means complex
    INTEGER                                       :: iwl ! starter for zwl
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)   :: zomega
  END TYPE simulation_type

  ! model->solution->base
  TYPE base_type
    LOGICAL               :: save

  END TYPE base_type

  ! model->solution->derived->cs
  TYPE cs_type
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: sca, ext, abs
  END TYPE cs_type

  ! model->solution->derived
  TYPE derived_type
    TYPE(cs_type)         :: cs
  END TYPE derived_type

  ! model->solution->mode
  TYPE mode_type
    COMPLEX(KIND=dp)                              :: zwl ! at what wavelength
    INTEGER, DIMENSION(:), ALLOCATABLE            :: modeidx ! modes indices
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE   :: eigval
    COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
  END TYPE mode_type

  ! model->solution
  TYPE solution_type
    TYPE(base_type)       :: base
    TYPE(derived_type)    :: derived
    TYPE(focal3d_type), DIMENSION(:), ALLOCATABLE :: focal ! array of planes
    TYPE(mode_type), DIMENSION(:), ALLOCATABLE    :: mode ! array of zwl
  END TYPE solution_type

  ! model
  TYPE model_type
     CHARACTER(LEN=120)     :: name ! model name (also the prefix for output).

     TYPE(geom_type)        :: geom

     TYPE(physics_type)     :: physics

     TYPE(simulation_type)  :: simulation

     TYPE(solution_type)    :: solution

  END TYPE model_type

CONTAINS

END MODULE data
