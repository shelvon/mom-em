! MODULE: sysmat
! AUTHOR: Jouni Makitalo
! DESCRIPTION:
! Routines for calculating approximate matrix representations of
! three boundary integral operators in space spanned by RWG basis.
! Also contains routines for calculating PMCHWT system matrix from
! these operator representations.
MODULE sysmat
  USE srcint
  USE common

  IMPLICIT NONE

CONTAINS
  ! Computes the PMCHWT system matrix for a problem described by
  ! struct b and for wavelength of given index wlind. Matrices are
  ! computed for all group representations.
  SUBROUTINE sysmat_pmchwt(b, wlind, A)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlind

    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: A
    INTEGER :: N, nf, ns, nd, nind
    INTEGER, DIMENSION(1:b%mesh%nedges) :: ind
    REAL (KIND=dp) :: omega, detj
    COMPLEX (KIND=dp) :: Zsi, ri, gae
    TYPE(prdnfo), POINTER :: prd

    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: M

    N = b%mesh%nedges

    omega = 2.0_dp*pi*c0/b%sols(wlind)%wl
    
    A(:,:,:) = 0.0_dp

    ALLOCATE(M(N,N))

    ! Loop through domains.
    DO nd=1,SIZE(b%domains)
       ri = b%media(b%domains(nd)%medium_index)%prop(wlind)%ri
       Zsi = (ri**2)/(eta0**2)

       ! Determine global edges indices to index basis functions.
       nind = b%domains(nd)%mesh%nedges
       ind(1:nind) = b%domains(nd)%mesh%edges(:)%parent_index

       ! Choose Green's function.
       IF(b%domains(nd)%gf_index==-1) THEN
          prd => NULL()
       ELSE
          prd => b%prd(b%domains(nd)%gf_index)
       END IF

       ! Loop through group actions.
       DO ns=1,SIZE(b%ga)
          IF(admissible_ga(b%domains(nd)%mesh, b%ga(ns), nd==1)==.FALSE.) THEN
             CYCLE
          END IF

          detj = b%ga(ns)%detj
       
          ! Compute block matrix using action of index ns.
          ! Matrices are independent of field actions.

          CALL computeD(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          ! Add the block matrix to different sub-problems multiplied by proper
          ! field actions.
          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)
          END DO

          CALL computeH(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)
          END DO

          CALL computeK(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A((ind(1:nind)+N),ind(1:nind),nf) = A((ind(1:nind)+N),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A(ind(1:nind),(ind(1:nind)+N),nf) = A(ind(1:nind),(ind(1:nind)+N),nf)&
                  - gae*detj*M(1:nind,1:nind)
          END DO

       END DO
    END DO

    DEALLOCATE(M)
  END SUBROUTINE sysmat_pmchwt

  ! Computes the PMCHWT system matrix for second/third-harmonic problem,
  ! Adds products of src_coef and the submatrices to src_vec.
  ! Note that omega and ri must correspond to the SH/TH case.
  ! Also the electric field group actions (gae) are the squares of
  ! actions in the linear problem.
  SUBROUTINE sysmat_pmchwt_nls(b, wlind, A, src_coef, src_vec)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlind
    COMPLEX (KIND=dp), DIMENSION(:,:,:,:), INTENT(IN) :: src_coef
    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: src_vec

    COMPLEX (KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: A
    INTEGER :: N, nf, ns, nd, nind, nc, nsources
    INTEGER, DIMENSION(1:b%mesh%nedges) :: ind
    REAL (KIND=dp) :: omega, detj
    COMPLEX (KIND=dp) :: Zsi, ri, gae
    TYPE(prdnfo), POINTER :: prd

    COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: M

    N = b%mesh%nedges

    nsources = SIZE(src_coef,4)

    ! Second-harmonic or third-harmonic frequency.
    IF(is_nl_thg(b)) THEN
        omega = 3.0_dp*2.0_dp*pi*c0/b%sols(wlind)%wl
    ELSE
        omega = 2.0_dp*2.0_dp*pi*c0/b%sols(wlind)%wl
    END IF

    A(:,:,:) = 0.0_dp

    ALLOCATE(M(N,N))

    ! Loop through domains.
    DO nd=1,SIZE(b%domains)
       IF(is_nl_thg(b)) THEN
           ri = b%media(b%domains(nd)%medium_index)%prop(wlind)%thri
       ELSE
           ri = b%media(b%domains(nd)%medium_index)%prop(wlind)%shri
       END IF
       Zsi = (ri**2)/(eta0**2)

       ! Determine global edges indices to index basis functions.
       nind = b%domains(nd)%mesh%nedges
       ind(1:nind) = b%domains(nd)%mesh%edges(:)%parent_index

       ! Choose Green's function.
       IF(b%domains(nd)%gf_index==-1) THEN
          prd => NULL()
       ELSE
          prd => b%prd(b%domains(nd)%gf_index)
       END IF

       ! Loop through actions.
       DO ns=1,SIZE(b%ga)
          IF(admissible_ga(b%domains(nd)%mesh, b%ga(ns), nd==1)==.FALSE.) THEN
             CYCLE
          END IF

          detj = b%ga(ns)%detj
       
          ! Compute block matrix using action of index ns.
          ! Matrices are independent of field actions.

          CALL computeD(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          ! Add the block matrix to different sub-problems multiplied by proper
          ! field actions.
          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)

             ! Loop through different sources.
             DO nc=1,nsources
                src_vec(ind(1:nind),nf,nc) = src_vec(ind(1:nind),nf,nc) - &
                     gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf,nc))

                src_vec((ind(1:nind)+N),nf,nc) = src_vec((ind(1:nind)+N),nf,nc) - &
                     gae*detj*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf,nc))*Zsi
             END DO
          END DO

          CALL computeH(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A(ind(1:nind),ind(1:nind),nf) = A(ind(1:nind),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A((ind(1:nind)+N),(ind(1:nind)+N),nf) = A((ind(1:nind)+N),(ind(1:nind)+N),nf)&
                  + gae*detj*Zsi*M(1:nind,1:nind)

             DO nc=1,nsources
                src_vec(ind(1:nind),nf,nc) = src_vec(ind(1:nind),nf,nc) - &
                     gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf,nc))
                
                src_vec((ind(1:nind)+N),nf,nc) = src_vec((ind(1:nind)+N),nf,nc) - &
                     gae*detj*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf,nc))*Zsi
             END DO
          END DO

          CALL computeK(omega, ri, b%domains(nd)%mesh, b%ga(ns), prd, b%qd_tri, M(1:nind,1:nind))

          DO nf=1,SIZE(b%ga)
             gae = b%ga(ns)%ef(nf)

             A((ind(1:nind)+N),ind(1:nind),nf) = A((ind(1:nind)+N),ind(1:nind),nf)&
                  + gae*M(1:nind,1:nind)
             A(ind(1:nind),(ind(1:nind)+N),nf) = A(ind(1:nind),(ind(1:nind)+N),nf)&
                  - gae*detj*M(1:nind,1:nind)

             DO nc=1,nsources
                src_vec(ind(1:nind),nf,nc) = src_vec(ind(1:nind),nf,nc) + &
                     gae*detj*MATMUL(M(1:nind,1:nind), src_coef(1:nind,nd,nf,nc))

                src_vec((ind(1:nind)+N),nf,nc) = src_vec((ind(1:nind)+N),nf,nc) - &
                     gae*MATMUL(M(1:nind,1:nind), src_coef((nind+1):(2*nind),nd,nf,nc))
             END DO
          END DO

       END DO
    END DO

    DEALLOCATE(M)
  END SUBROUTINE sysmat_pmchwt_nls


  ! Computes the system matrix for eigen modes problem of single particle
  ! described by struct b and for wavelength of given index wlind.
  SUBROUTINE sysmat_mueller(b, wlind, A, N)
    TYPE(batch), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: wlind
    INTEGER, INTENT(IN) :: N
    COMPLEX (KIND=dp), DIMENSION(1:2*N,1:2*N), INTENT(INOUT) :: A

    COMPLEX (KIND=dp), DIMENSION(1:N,1:N) :: nxD1, nxD2, nxHline1, nxHline2, nxK1, nxK2
    INTEGER :: nd, nind
    INTEGER, DIMENSION(1:b%mesh%nedges) :: ind
    REAL (KIND=dp) :: omega
    COMPLEX (KIND=dp) :: ri1, ri2, k1, k2, eps1, eps2, mu1, mu2
    TYPE(mesh_container) :: mesh
    TYPE(prdnfo), POINTER :: prd

    ! no group representation
    IF(SIZE(b%ga)/=1) THEN
       WRITE(*,*) 'group action is not supported, exit!'
       RETURN
    END IF
    ! no periodi green function
    prd => NULL()

    omega = 2.0_dp*pi*c0/b%sols(wlind)%wl

    ! Loop through media 1, 2.
    ! 1: the infinite surrounding medium;
    ! 2: the single particle
    ri1 = b%media(1)%prop(wlind)%ri
    ri2 = b%media(2)%prop(wlind)%ri

    k1 = ri1*omega/c0
    k2 = ri2*omega/c0
    eps1 = (ri1**2)*eps0
    eps2 = (ri2**2)*eps0
    mu1 = mu0
    mu2 = mu0
    mesh = b%mesh
    ! mesh = b%domains(1)%mesh

!    !$OMP PARALLEL DEFAULT(SHARED)
!
!    !$OMP SECTIONS
!    !$OMP SECTION
    CALL computenxD(omega, ri1, mesh, b%ga(1), prd, b%qd_tri, nxD1)
!    !$OMP SECTION
    CALL computenxD(omega, ri2, mesh, b%ga(1), prd, b%qd_tri, nxD2)
!    !$OMP SECTION
    CALL computenxHline(omega, ri1, mesh, b%ga(1), prd, b%qd_tri, nxHline1)
!    !$OMP SECTION
    CALL computenxHline(omega, ri2, mesh, b%ga(1), prd, b%qd_tri, nxHline2)
!    !$OMP SECTION
    CALL computenxK(omega, ri1, mesh, b%ga(1), prd, b%qd_tri, nxK1)
!    !$OMP SECTION
    CALL computenxK(omega, ri2, mesh, b%ga(1), prd, b%qd_tri, nxK2)
!    !$OMP END SECTIONS
!
!    !$OMP END PARALLEL

!    !$OMP SECTIONS
!    !$OMP SECTION
    A(1:N,1:N) = A(1:N,1:N) + nxK1*2.0_dp*mu1/(mu1+mu2)
    A(1:N,1:N) = A(1:N,1:N) - nxK2*2.0_dp*mu2/(mu1+mu2)
!    !$OMP SECTION
    A(1:N,(N+1):(2*N)) = A(1:N,(N+1):(2*N)) - (nxD1+nxHline1)*2.0_dp*eps1/(mu1+mu2)
    A(1:N,(N+1):(2*N)) = A(1:N,(N+1):(2*N)) + (nxD2+nxHline2)*2.0_dp*eps2/(mu1+mu2)
!    !$OMP SECTION
    A((N+1):(2*N),1:N) = A((N+1):(2*N),1:N) + (nxD1+nxHline1)*2.0_dp*eps1/(eps1+eps2)
    A((N+1):(2*N),1:N) = A((N+1):(2*N),1:N) - (nxD2+nxHline2)*2.0_dp*eps2/(eps1+eps2)
!    !$OMP SECTION
    A((N+1):(2*N),(N+1):(2*N)) = A((N+1):(2*N),(N+1):(2*N)) + nxK1*2.0_dp*eps1/(eps1+eps2)
    A((N+1):(2*N),(N+1):(2*N)) = A((N+1):(2*N),(N+1):(2*N)) - nxK2*2.0_dp*eps2/(eps1+eps2)
!    !$OMP END SECTIONS
!
!    !$OMP END PARALLEL

  END SUBROUTINE sysmat_mueller

  ! Determine if faces of indices n and m are considered near,
  ! so as to use singularity subtraction, which is time-consuming.
  FUNCTION near_faces(mesh, prd, n, m, ga) RESULT(res)
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    INTEGER, INTENT(IN) :: n, m
    TYPE(group_action), INTENT(IN) :: ga

    LOGICAL :: res
    REAL (KIND=dp), DIMENSION(3) :: diff
    REAL (KIND=dp) :: distsq, threshold

    threshold = mesh%avelen*3

    diff = mesh%faces(m)%cp - MATMUL(ga%j, mesh%faces(n)%cp)
    distsq = SUM(diff*diff)

    IF(distsq<threshold**2) THEN
       res = .TRUE.
    ELSE IF(ASSOCIATED(prd)) THEN
       IF(prd%type==prd_2d .AND. distsq>(MINVAL((/prd%dx,prd%dy/))-threshold)**2) THEN
          res = .TRUE.
       ELSE
          res = .FALSE.
       END IF
    ELSE
       res = .FALSE.
    END IF

  END FUNCTION near_faces

  ! Computes the D-matrix elements
  ! D_mng = -i*omega_mu*int_Sm dS fm(r). int_Sn' dS' J_g*fn(r')*O'_gG(r,r').
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! D: The D-matrix.
  SUBROUTINE computeD(omega, ri, mesh, ga, prd, qd, D)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: D
    COMPLEX (KIND=dp) :: c1, k
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    COMPLEX (KIND=dp) :: int1
    COMPLEX (KIND=dp), DIMENSION(3,3,qd%num_nodes,mesh%nfaces) :: intaux
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes,3) :: fmv
    LOGICAL :: near

    WRITE(*,*) 'Building a D-matrix'

    nweights = qd%num_nodes

    D(:,:) = 0.0_dp

    k = ri*omega/c0

    ! Coefficients of partial integrals.
    c1 = -(0,1)*omega*mu0

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p)) ! f_n(r), vector basis function
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,prd,m,qd)&
       !$OMP PRIVATE(n,near,r)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO r=1,nweights
             intaux(:,:,r,n) = intK2(qpm(:,r), n, mesh, k, ga, prd, near, qd)
          END DO                   
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       DO n=1,mesh%nfaces

          DO p=1,3
             DO q=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qd%weights(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), intaux(:,q,r,n))
                END DO
                int1 = int1*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                D(index1,index2) = D(index1,index2) + c1*int1
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeD

  ! Computes the H-matrix elements
  ! H_mng = -1/(i*omega*epsilon)*int_Sm dS div_S(fm(r)) int_Sn' dS' div'_S(fn(r'))*O'_gG(r,r').
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! H: The H-matrix.
  SUBROUTINE computeH(omega, ri, mesh, ga, prd, qd, H)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: H
    COMPLEX (KIND=dp) :: c2, k
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    COMPLEX (KIND=dp) :: int1
    REAL (KIND=dp), DIMENSION(3) :: fmDiv
    COMPLEX (KIND=dp), DIMENSION(3,qd%num_nodes,mesh%nfaces) :: intaux
    LOGICAL :: near

    WRITE(*,*) 'Building an H-matrix'

    nweights = qd%num_nodes

    H(:,:) = 0.0_dp

    k = ri*omega/c0

    ! Coefficients of partial integrals.
    c2 = -1.0_dp/((0,1)*omega*eps0*(ri**2))

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          fmDiv(p) = rwgDiv(m, p, mesh)
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,prd,m,qd)&
       !$OMP PRIVATE(r,n,near)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO r=1,nweights
             intaux(:,r,n) = intK1(qpm(:,r), n, mesh, k, ga, prd, near, qd)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    
       DO n=1,mesh%nfaces

          DO p=1,3
             DO q=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qd%weights(r)*intaux(q,r,n)
                END DO
                int1 = int1*Am*fmDiv(p)

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                H(index1,index2) = H(index1,index2) + c2*int1
             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeH

  ! Computes the K-matrix elements
  ! K_mng = int_Sm dS fm(r). int_Sn' dS' [O'_g grad'G(r,r')]x(J_g*fn(r')).
  ! IN:
  ! omega: The angular frequency.
  ! ri: Refractive index of the domain.
  ! mesh: Surface mesh of the scatterer.
  ! ga: Group action identifiers.
  ! prd: The periodic Green's function.
  ! OUT:
  ! A: The K-matrix.
  SUBROUTINE computeK(omega, ri, mesh, ga, prd, qd, A)
    ! Input variables
    TYPE(mesh_container), INTENT(IN) :: mesh
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    REAL (KIND=dp), INTENT(IN) :: omega
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd

    ! Internal variables
    COMPLEX (KIND=dp), INTENT(INOUT), DIMENSION(mesh%nedges,mesh%nedges) :: A
    COMPLEX (KIND=dp) :: k, int1
    INTEGER :: nweights, n, m, p, q, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    COMPLEX (KIND=dp), DIMENSION(3,3,qd%num_nodes,mesh%nfaces) :: intaux
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes,3) :: fmv
    LOGICAL :: near

    WRITE(*,*) 'Building a K-matrix'

    nweights = qd%num_nodes

    A(:,:) = 0.0_dp

    k = ri*omega/c0

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area

       DO p=1,3
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intaux,qpm,mesh,k,ga,m,prd,qd)&
       !$OMP PRIVATE(r,n,near)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO r=1,nweights
             intaux(:,:,r,n) = intK4(qpm(:,r), n, mesh, k, ga, m, prd, near, qd)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL
    
       DO n=1,mesh%nfaces

          DO q=1,3
             DO p=1,3

                int1 = 0.0_dp
                DO r=1,nweights
                   int1 = int1 + qd%weights(r)*dotc(CMPLX(fmv(:,r,p),KIND=dp), intaux(:,q,r,n))
                END DO
                int1 = int1*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                A(index1,index2) = A(index1,index2) + int1

             END DO
          END DO
       END DO

    END DO

  END SUBROUTINE computeK

  ! Compute the nxD-matrix elements,
  ! D_mng = int_Sm dS fm(r).n(r) x int_Sn' dS' J_g*fn(r')*O'_gG(r,r'),
  ! where "x" means cross product and "." means dot product.
  SUBROUTINE computenxD(omega, ri, mesh, ga, prd, qd, nxD)
    ! Input variables
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd
    ! Input and output variables
    COMPLEX (KIND=dp), DIMENSION(1:mesh%nedges,1:mesh%nedges), INTENT(INOUT) :: nxD

    ! Internal variables
    COMPLEX (KIND=dp) :: k, eps, mu, coeff, intnxD
    INTEGER :: n, m, p, q, t, r, index1, index2, nweights
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3) :: nor, diff
    LOGICAL :: near
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes,3) :: fmv
    COMPLEX (KIND=dp), DIMENSION(3,3,qd%num_nodes,mesh%nfaces) :: intnxDaux

    WRITE(*,*) 'Building the nxD-matrix for eigen mode solver'

    nweights = qd%num_nodes

    nxD(:,:) = CMPLX(0.0_dp,0.0_dp)
    k = ri*omega/c0
    eps = (ri**2)*eps0
    mu = mu0
    coeff = (0,1)*omega*mu

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area
       nor = mesh%faces(m)%n

       DO p=1,3
          ! output: f_n(r), vector basis function
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intnxDaux,qpm,mesh,k,ga,prd,m,qd)&
       !$OMP PRIVATE(n,near,r)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO r=1,nweights
             ! intnxDaux = int_Sn' dS' J_g*fn(r')*O'_gG(r,r')
             intnxDaux(:,:,r,n) = intK2(qpm(:,r), n, mesh, k, ga, prd, near, qd)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       DO n=1,mesh%nfaces
          DO q=1,3
             DO p=1,3

                intnxD = 0.0_dp
                DO r=1,nweights
                   ! int_Sm dS fm(r).n(r) x intnxDaux.
                   intnxD = intnxD + qd%weights(r)*&
                   dotc(CMPLX(fmv(:,r,p),KIND=dp), crossc(CMPLX(nor,KIND=dp), intnxDaux(:,q,r,n)))
                END DO!r
                intnxD = intnxD*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                nxD(index1,index2) = nxD(index1,index2) + intnxD
             END DO!p
          END DO!q
       END DO!n

    END DO!m

    nxD(:,:) = coeff*nxD(:,:)
  END SUBROUTINE computenxD

  ! Compute the nxH-matrix elements,
  ! H_mng = int_Sm dS fm(r).n(r) x grad int_Sn' dS' div'_S(fn(r'))*O'_gG(r,r'),
  ! which is  done via close line integration as
  ! H_mng = -sum_(+-) int_(Tm+-) dC m(r).fm(r) int_Sn' dS' div'_S(fn(r'))*O'_gG(r,r'),
  ! where m(r) is a unit vector tangent to close line C which has induced orientation from that of S.
  SUBROUTINE computenxHline(omega, ri, mesh, ga, prd, qd, nxHline)
    ! Input variables
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd
    ! Input and output variables
    COMPLEX (KIND=dp), DIMENSION(1:mesh%nedges,1:mesh%nedges), INTENT(INOUT) :: nxHline

    ! Internal variables
    INTEGER, PARAMETER :: nweights = 6
    COMPLEX (KIND=dp) :: k, eps, mu, coeff
    INTEGER :: n, m, p, q, t, r, index1, index2
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3) :: nor, diff
    LOGICAL :: near
    REAL (KIND=dp), DIMENSION(3) :: p1, p2, ptm
    REAL (KIND=dp), DIMENSION(3,3) :: tv
    REAL (KIND=dp), DIMENSION(3,nweights*3) :: ptv
    REAL (KIND=dp), DIMENSION(nweights) :: glp, glw
    REAL (KIND=dp), DIMENSION(3,nweights*3,3) :: fmv
    COMPLEX (KIND=dp), DIMENSION(3,nweights*3,mesh%nfaces) :: intnxHaux
    COMPLEX (KIND=dp), DIMENSION(3) :: intnxH

    WRITE(*,*) 'Building the nxHline-matrix for eigen mode solver'

    ! Need further improvement.
    ! Points for n=6 Gauss-Legendre quadrature.
    ! For close line integration.
    glp = (/-0.238619186083197_dp, 0.238619186083197_dp,&
         -0.661209386466265_dp, 0.661209386466265_dp,&
         -0.932469514203152_dp, 0.932469514203152_dp/)
    glw = (/0.467913934572691_dp, 0.467913934572691_dp,&
         0.360761573048139_dp, 0.360761573048139_dp,&
         0.171324492379170_dp, 0.171324492379170_dp/)

    nxHline(:,:) = CMPLX(0.0_dp,0.0_dp)
    k = ri*omega/c0
    eps = (ri**2)*eps0
    mu = mu0
    coeff = -1.0_dp/((0,1)*omega*eps)

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area
       nor = mesh%faces(m)%n

       DO t=1,3
          p1 = mesh%nodes( mesh%faces(m)%node_indices(indexrot3(t)) )%p
          p2 = mesh%nodes( mesh%faces(m)%node_indices(indexrot3(t+1)) )%p
          tv(:,t) = p2 - p1

          DO r=1,nweights
             ptm = p1 + tv(:,t)*(glp(r)*0.5_dp + 0.5_dp)
             ptv(:,r+nweights*(t-1)) = ptm(:)
          END DO
       END DO

       DO p=1,3
          ! output: f_n(r), but for n(r) x H(fn(r)) calculation
          ! and a different quadrature is used.
          CALL vrwg(ptv(:,:),m,p,mesh,fmv(:,:,p))
       END DO

!       !$OMP PARALLEL DEFAULT(NONE)&
!       !$OMP SHARED(nweights,intnxHaux,ptv,mesh,k,ga,prd,m,qd)&
!       !$OMP PRIVATE(n,near,r,t)
!       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO t=1,3
             ! Integration points.
             DO r=1,nweights
                ! intnxHaux = int_Sn' dS' div'_S(fn(r'))*O'_gG(r,r'), where the singularity
                ! is cancelled out in the combined kernel of G(r,r') in Mueller formulation.
                intnxHaux(:,r+nweights*(t-1),n) = intK1Mueller(ptv(:,r+nweights*(t-1)), n, mesh, k, ga, prd, near, qd)
             END DO
          END DO
!       END DO
!       !$OMP END DO
!       !$OMP END PARALLEL
!
!       DO n=1,mesh%nfaces
          DO q=1,3
             ! Edges of the test triangle in contour integration.
             intnxH(:) = 0.0_dp
             DO t=1,3
                ! Integration points.
                DO r=1,nweights
                   ! Edges of the test triangle
                   DO p=1,3
                      intnxH(p) = intnxH(p) + glw(r)*dotr(tv(:,t),fmv(:,r+nweights*(t-1),p))*intnxHaux(q,r+nweights*(t-1),n)
                   END DO
                END DO
             END DO

             ! Jacobian of [0,1] -> [-1,1]
             intnxH(:) = intnxH(:)*0.5_dp

             index2 = mesh%faces(n)%edge_indices(q)

             DO p=1,3
                index1 = mesh%faces(m)%edge_indices(p)

                nxHline(index1,index2) = nxHline(index1,index2) + intnxH(p)
             END DO!p

          END DO!q
       END DO!n

    END DO!m

    nxHline(:,:) = coeff*nxHline(:,:)
  END SUBROUTINE computenxHline


  ! Compute the nxK-matrix elements,
  ! K_mng = int_Sm dS fm(r).n(r)x int_Sn' dS' [O'_g grad'G(r,r')]x(J_g*fn(r')),
  ! where "x" means cross product and "." means dot product.
  SUBROUTINE computenxK(omega, ri, mesh, ga, prd, qd, nxK)
    ! Input variables
    REAL (KIND=dp), INTENT(IN) :: omega
    COMPLEX (KIND=dp), INTENT(IN) :: ri
    TYPE(mesh_container), INTENT(IN) :: mesh
    TYPE(group_action), INTENT(IN) :: ga
    TYPE(prdnfo), POINTER, INTENT(IN) :: prd
    TYPE(quad_data), INTENT(IN) :: qd
    ! Input and output variables
    COMPLEX (KIND=dp), DIMENSION(1:mesh%nedges,1:mesh%nedges), INTENT(INOUT) :: nxK

    ! Internal variables
    COMPLEX (KIND=dp) :: k, eps, mu, coeff, intnxK
    INTEGER :: n, m, p, q, t, r, index1, index2, nweights
    REAL (KIND=dp) :: Am
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes) :: qpm
    REAL (KIND=dp), DIMENSION(3) :: nor, diff
    LOGICAL :: near
    REAL (KIND=dp), DIMENSION(3,qd%num_nodes,3) :: fmv
    COMPLEX (KIND=dp), DIMENSION(3,3,qd%num_nodes,mesh%nfaces) :: intnxKaux

    WRITE(*,*) 'Building the nxK-matrix for eigen mode solver'

    nweights = qd%num_nodes

    nxK(:,:) = CMPLX(0.0_dp,0.0_dp)
    k = ri*omega/c0

    DO m=1,mesh%nfaces

       qpm = quad_tri_points(qd, m, mesh)
       Am = mesh%faces(m)%area
       nor = mesh%faces(m)%n

       DO p=1,3
          ! output: f_n(r), vector basis function
          CALL vrwg(qpm,m,p,mesh,fmv(:,:,p))
       END DO

       !$OMP PARALLEL DEFAULT(NONE)&
       !$OMP SHARED(nweights,intnxKaux,qpm,mesh,k,ga,prd,m,qd)&
       !$OMP PRIVATE(n,near,r)
       !$OMP DO SCHEDULE(STATIC)
       DO n=1,mesh%nfaces
          near = near_faces(mesh, prd, n, m, ga)

          DO r=1,nweights
             ! intnxKaux = int_Sn' dS' [O'_g grad'G(r,r')]x(J_g*fn(r'))
             intnxKaux(:,:,r,n) = intK4(qpm(:,r), n, mesh, k, ga, m, prd, near, qd)
          END DO
       END DO
       !$OMP END DO
       !$OMP END PARALLEL

       DO n=1,mesh%nfaces
          DO q=1,3
             DO p=1,3

                intnxK = 0.0_dp
                DO r=1,nweights
                   ! int_Sm dS fm(r).n(r) x intnxKaux.
                   intnxK = intnxK + qd%weights(r)*&
                   dotc(CMPLX(fmv(:,r,p),KIND=dp), crossc(CMPLX(nor,KIND=dp), intnxKaux(:,q,r,n)))
                END DO!r
                intnxK = intnxK*Am

                index1 = mesh%faces(m)%edge_indices(p)
                index2 = mesh%faces(n)%edge_indices(q)

                nxK(index1,index2) = nxK(index1,index2) + intnxK
             END DO!p
          END DO!q
       END DO!n

    END DO!m

  END SUBROUTINE computenxK

END MODULE sysmat
