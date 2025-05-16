!***********************************************************************************************************************************************************
!  CODE: VortoNeX - The Full Non-linear Vortex Tube-Vorton Method (FTVM)                                                                                   !
!  PURPOSE: to solve detached flow/fluid past a shell-body (pre and post-stall condition)                                                                  !
!  VERSION: 1.1                                                                                                                                            !
!  AUTHOR AND DEVELOPER: Jesús Carlos Pimentel-García                                                                                                      !
!  CONTACT: pimentel_garcia@yahoo.com.mx                                                                                                                   !
!  CONTRIBUTOR: Enrique Ortega A. (UPC)                                                                                                                    !
!                                                                                                                                                          !
!  Note: Since the main objective of the current code is to explicitly show the procedure used to apply the FTVM to solve fluid flow past a thin body,     !
!        some aspects such as simplification, optimization, and robustness are left for future development. Of course, the current code (and its concepts) !
!        can be modified/improved/extended for non-commercial purposes; however, the licensee or licensor reserves the right to claim for future           !
!        commercial uses based on the protection under an international patent application:                                                                ! 
!        "Full-surface detached vorticity method and system to solve fluid dynamics", International patent application: WO/2024/136634,                    !
!        World Intellectual Property Organization (2024). https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2024136634&_cid=P21-LXXE19-45988-1     !    
!                                                                                                                                                          !              
!  Control version:                                                                                                                                        !
!       v1.0 (Sep/2023): Only implements the Kutta-Zhukovski (KJ) force calculation (validated for the pre-stall condition): https://link.springer.com/content/pdf/10.1186/s42774-023-00168-8.pdf                                           
!       v1.1 (May/2025): Implements both KJ (pre-stall) and an alternative method (post-stall) for force calculations: https://www.researchgate.net/publication/390250217_The_Full_Nonlinear_Vortex_Tube-Vorton_Method_flat_plate_at_post-stall_condition    
!***********************************************************************************************************************************************************
   
Program VortoNeX ! Vorton-Next (generation); JCPG/EOA
    Use ARRAYS, only : pr,mygridcon
    Use Mesh_ops, only : genmesh_st
    Use OMP_LIB
    Implicit None
    real(kind=pr) :: circ,point(3),tot_area,circ_vec(3),pos_ivort(3),pos_kvort(3),vind_total(3),circ_vec_2(3),pos_kvort_2(3)
    integer :: j,n,step,op,fstlen,stret,i,k,den_aux,option,option_2
    
    call READ_PARAM ! reads input parameters from input.txt (alpha, dt, vortex core radius, etc.)
    call GEOM ( ) ! reads the GiD's generated mesh from .dat file
    call GENMESH_ST ( mygridcon )  ! mesh helpers
    call KUTTA_EDGES ( )  ! kutta edges' data 
        !call COUPLING_BODY ( circ, j, point, tot_area ) ! assemblies the shell-body's influence coefficient matrix (ICM) with bounded vortex rings
    call V2_COUPLING_BODY ( circ, j, point, tot_area ) ! assemblies the shell-body's influence coefficient matrix (ICM) with multiple vorton rings
    call WAKE_BUFFER ( ) ! generates the first wake row for starting solution; it has been modified from the original version to avoid remaining wake rows generation
    call LE_DETEC ( ) ! detects the leading edge's (LE) bounded vortex rings (BVR)/"panels" (thereinafter BVR and panels are used as synonyms)
        !call STARTING_SOL ( circ, n, point, tot_area,j,op,fstlen,stret,step ) ! solves a starting solution (could be the steady-state one) to initialize the unsteady case (with bounded vortex rings)
    call V2_STARTING_SOL ( circ, n, point, tot_area, op, step ) ! solves a starting solution (could be the steady-state one) to initialize the unsteady case (with multiple vorton rings)
    call TIME_STEPS ( point,step,circ_vec,pos_ivort,pos_kvort,i,k,tot_area,den_aux,option,vind_total,circ,j,option_2,circ_vec_2,pos_kvort_2 ) ! main (unsteady) loop
End Program VortoNeX
    
Subroutine READ_PARAM ( ) ! reads input parameters (alpha, dt, vortex core radius, kinematic viscosity, etc.); JCPG
    Use ARRAYS,only : pr,input
    Implicit none
    integer :: stat
    
    open(unit=1, file='input.txt', action='read', status='old', access='sequential')
    read(1,*);read(1,*) input%mesh_num         ! mesh's number (mesh: quadrangular: 1, 2, 3, 4; rectangular: 5; swept-back: 6, 7)
    read(1,*);read(1,*) input%nsteps           ! number of total unsteady iterations (iter.; e.g., 15, 80, 150,...)
    read(1,*);read(1,*) input%alpha            ! angle of attack (alpha; positive values, e.g., 5, 17.5, 40,...)
    read(1,*);read(1,*) input%dt               ! time step (Delta t; e.g., 0.0625)
    read(1,*);read(1,*) input%core_rad_init    ! vortex core radius (sigma_zero; e.g., 0.0442)
    read(1,*);read(1,*) input%eps              ! distance from plate to put the nascent vortons (epsilon; e.g., 0.0442)
    read(1,*);read(1,*) input%nthreads         ! number of threads (threads; e.g., serial: 1, parallel: 2, 4, 8)
    read(1,*);read(1,*) input%detach_model     ! flow detachment model (pre-stall: 1, post-stall (testing): 2)
    read(1,*);read(1,*) input%dens             ! flow density (rho; e.g., unitary)
    read(1,*);read(1,*) input%q_inf            ! free-stream velocity (q_inf; e.g., unitary)
    read(1,*);read(1,*) input%char_len         ! characteristic length ("chord"; e.g., unitary)
    read(1,*);read(1,*) input%x_pos_mom_fact   ! pitching moment coefficient factor (x_cm; e.g., 0.25)
    read(1,*);read(1,*) input%wake_len         ! straight wake length (PHI; e.g., unitary for pre-stall or a fraction (e.g., 0.01) for post-stall)
    read(1,*);read(1,*) input%fst_wake_factor  ! first wake row length factor (phi; e.g., unitary)
    read(1,*);read(1,*) input%tol_rad          ! tolerance to avoid mathematical indetermination for induced velocities (tol_rad; e.g., 1e-8)
    read(1,*);read(1,*) input%tol_vort         ! tolerance to avoid mathematical indetermination for vortex stretching (tol_vort; e.g., 1e-8 or higher)
    read(1,*);read(1,*) input%kin_visc         ! fluid kinematic viscosity (nu; e.g., 0)
    read(1,*);read(1,*) input%regul_function   ! regularization function (g_sigma; 0: for high-order, 1: for 2nd-order Gaussian, 2: for Gaussian erf)
    read(1,*);read(1,*) input%vx_stretch       ! vortex stretching (vx_stretch; 0: constant volumes, 1: variable volumes)
    read(1,*);read(1,*) input%pres_inf         ! pressure at infinite (unitary)
    read(1,*);read(1,*) input%thick            ! plate's thickness (e.g., 0.04); for double layer pressure calculation; left 0.0 for single layer (default VLM-based solution)
    close(1)
     
    open(4, iostat=stat, file='aero_coef.dat') 
    if (stat==0) then 
        close(4, status='delete') ! it deletes file if exists
    else; end if
    
End Subroutine READ_PARAM
    
Subroutine GEOM ( ) ! creates geometry from a mesh file (GiD); JCPG
    Use ARRAYS, only : pr,grid,trailing_edges,kuttaedges,internal_LE,plate_LE,plate_lat,plate_TE,input,solid_angle
    Implicit none
    ! local variables
    integer :: i,j,k,nsize
    
    ! GEOMETRY FILES (select one for testing)
    select case(input%mesh_num)
    case(1) ! quadrangular AR=1
        open(2,file='meshes\mesh_2x2_AWD.dat')
            !open(2,file='meshes\mesh_2x2_AWD_inv.dat')
    case(2) ! quadrangular AR=1  
        open(2,file='meshes\mesh_4x4_AWD.dat')
    case(3) ! quadrangular AR=1
        open(2,file='meshes\mesh_10x10_AWD.dat')
        !open(2,file='mesh_2x2_AWD_inv.dat')
    case(4) ! quadrangular AR=1
        open(2,file='meshes\mesh_16x16_AWD.dat') 
    case(5) ! rectangular AR=0.5
        open(2,file='meshes\mesh_16x8_AWD.dat') 
    case(6) ! swept-back AR=1
        open(2,file='meshes\mesh_10x10_swept_AWD.dat') 
    case(7) ! swept-back AR=1
        open(2,file='meshes\mesh_16x16_swept_AWD.dat') 
    end select
    
    read(2,*);read(2,*);read(2,*)  
    read(2,*) grid%nnod,grid%nelem  ! #nodes y #elements
    read(2,*)
    nsize=max(grid%nnod,grid%nelem) ! max size
    allocate(grid%panel(4,grid%nelem),grid%coord(3,grid%nnod),grid%pnod(grid%nelem),grid%coord_start(3,grid%nnod))

    grid%panel(1:4,1:grid%nelem) = 0 ! mandatory
    
    do i=1,grid%nnod ! read nodes coordinates
        read(2,*) j, grid%coord(:,i)
    end do
    
    read(2,*)
    
    do i=1,grid%nelem ! read type of element (quad or tri) and panel nodes
        read(2,*) j,grid%pnod(i),(grid%panel(k,i),k=1,MIN(4,grid%pnod(i)+1))
    end do
     
    grid%pnod(:)=grid%pnod(:)+1 ! nodes number
    read(2,*);read(2,*)
    read(2,*);read(2,*);read(2,*);read(2,*) 
    read(2,*) kuttaedges%nte 
    allocate(trailing_edges(2,kuttaedges%nte))
    do i=1,kuttaedges%nte 
        read(2,*) trailing_edges(1:2,i)
    end do
    read(2,*);read(2,*)
    read(2,*) grid%nle ! number of internal LE
    allocate(internal_LE(1:2,grid%nle))
    do i=1,grid%nle
        read(2,*) internal_LE(1:2,i)
    end do
    
    read(2,*);read(2,*)
    read(2,*) grid%plate_nle ! number of plate's LE
    allocate(plate_LE(1:2,grid%plate_nle))
    do i=1,grid%plate_nle
        read(2,*) plate_LE(1:2,i)
    end do
    
    read(2,*);read(2,*)
    read(2,*) grid%plate_lat_nodes ! number of plate's lateral edges (it does not apply for the current code but only for the previous one: UFVLM)
    allocate(plate_lat(grid%plate_lat_nodes))
    do i=1,grid%plate_lat_nodes
        read(2,*) plate_lat(i)
    end do
    
    read(2,*);read(2,*)
    read(2,*) grid%plate_nte ! number of plate's trailing edges
    allocate(plate_TE(1:2,grid%plate_nte))
    do i=1,grid%plate_nte
        read(2,*) plate_TE(1:2,i)
    end do 
    
    allocate(solid_angle(grid%nelem,grid%nelem))
    
    close(2)
End Subroutine GEOM

Subroutine V2_COUPLING_BODY( circ, j, point, tot_area ) ! generates the body's influence coefficient matrix (ICM) for the multi-vorton scheme (V2); JCPG
    Use ARRAYS, only : pr,ctrl_pts,nor_vec,A_body,grid,rhs,pivot,A_solv,bound_circ,vdir,input,tan_vec1,tan_vec2,area,bound_circ_old,deriv_gamma,rhs_freeflow,solid_angle
    Use MATH, only : DOTPRODUCT,VECTORNORM,CROSSPRODUCT
    Implicit none
    integer :: i,op
    integer, intent(out) :: j
    real(kind=pr), intent(out) :: circ,point(3),tot_area
    real(kind=pr) :: vind_ring(3),pco(3),geometry(10),equiv_vcr
    real(kind=pr) :: r_pa(3),r_pb(3),r_pc(3),norm_r_pa,norm_r_pb,norm_r_pc,rpb_x_rpc(3)
    real(kind=pr) :: p(3),v(3,4),sa
    
    allocate(ctrl_pts(3,grid%nelem),nor_vec(3,grid%nelem),A_body(grid%nelem,grid%nelem),rhs(grid%nelem),tan_vec1(3,grid%nelem),tan_vec2(3,grid%nelem), area(grid%nelem),rhs_freeflow(grid%nelem))
    
    vdir(1:3) = (/ input%q_inf*cosd(input%alpha), 0._pr , input%q_inf*sind(input%alpha) /)  ! total velocity ( Vwind-Vbody, viewed from the body)
    grid%elemtype=grid%pnod(1)
    
    do i=1, grid%nelem
        call PANELGEO(geometry,pco,i)
        ctrl_pts(:,i) = pco ! panel's control points on the body
        
        tan_vec1(:,i) = geometry(1:3) ! first tangential vector
        tan_vec2(:,i) = geometry(4:6) ! second tangential vector
        nor_vec(:,i)  = geometry(7:9) ! normal vector
        area(i)       = geometry(10)  ! panel's area
    end do
    
    tot_area = 0._pr
    do i=1, grid%nelem
        tot_area = tot_area + area(i) ! total body's area
    end do
    
    ! new lines for precise solid angles calculation
    do i=1, grid%nelem ! on control points
        do j=1, grid%nelem ! on panels
            p(:) = ctrl_pts(:,i)
            v(:,1) = grid%coord(:,grid%panel(1,j))
            v(:,2) = grid%coord(:,grid%panel(2,j))
            v(:,3) = grid%coord(:,grid%panel(3,j))
            v(:,4) = grid%coord(:,grid%panel(4,j))
            call polygon_solid_angle_3d ( 4, v, p, sa )
            solid_angle(i,j) = sa ! solid angle            
        end do
    end do

    ! Shell-body's influence coefficients' matrix (A_body; only bounded VR)
    do i=1, grid%nelem ! over control points
        point(:) = ctrl_pts(:,i)
        do j=1, grid%nelem ! over bounded vorton rings
            circ = 1.0_pr ! unitary VR's circulation (for starting solution)
            op = 0 ! for plate case
            call VORTONRING ( j,circ,point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),A_body(i,j)) ! ICM's elements
        end do
    end do ! ends j
        
    allocate( A_solv(grid%nelem,grid%nelem),pivot(grid%nelem),bound_circ(grid%nelem),bound_circ_old(grid%nelem),deriv_gamma(grid%nelem) )

End Subroutine V2_COUPLING_BODY
    
Subroutine WAKE_BUFFER ( ) ! generates the first wake row (for bounded vortex rings-based starting solution); EOA
    Use ARRAYS, only : pr,vdir,kuttaedges,nor_vec,mygridcon,grid,input,wake
    Use MATH, only : DOTPRODUCT,NORMALIZE
    Implicit none
    integer :: j,nnew,pnew,iedg,jnod,node_id,pan_id,nwn
    integer :: int_aux_1(input%nsteps*grid%nnod),addnods(4)
    real(kind=pr) :: tnvec(3),tmp,wdir(3)
    
    int_aux_1(:) = -1 ;  ! auxiliar
    nnew = 0 ; pnew = 0
    grid%coord_start(:,:) = grid%coord(:,:)
    do iedg = 1,kuttaedges%nte  ! only first layer
        tnvec(1:3) = 0.0_pr  ! edge averaged directions
        do j = kuttaedges%esed2(iedg) + 1, kuttaedges%esed2(iedg+1)
            tnvec(1:3) = tnvec(1:3) + nor_vec(1:3,kuttaedges%esed1(j)) ! tangential vector
        end do  
        
        call NORMALIZE (tnvec(1:3), tmp)
        call DOTPRODUCT(vdir,tnvec(1:3), tmp)
        wdir(1:3) = vdir(1:3)
        call NORMALIZE( wdir(1:3), tmp)
        
        addnods(1:2) = mygridcon%inpoed(1:2,kuttaedges%edges(iedg))   ! addnods is an 4 components integer (local variable)
        
        do j = 1,2  ! edge's nodes
            jnod = addnods(j) 
            if ( int_aux_1(jnod).lt.0 )  then  ! add a node
                nnew = nnew + 1 ; node_id = grid%nnod + nnew 
                grid%coord(:,node_id) = grid%coord(:,jnod) + input%dt*wdir(1:3) ! new position
                grid%coord_start(:,node_id) = grid%coord_start(:,jnod) + input%wake_len*wdir(1:3) ! new position for starting (or steady-state) solution
                addnods(j+2) = node_id ; int_aux_1(jnod) = node_id 
                kuttaedges%flnods(nnew) = jnod
            else  ! read a node
                addnods(j+2) = int_aux_1(jnod) 
            end if
        end do  ! edge's nodes
        
        pnew = pnew + 1; 
        pan_id = grid%nelem + pnew  ! add panel (it has the same edge orientation)
        
        grid%panel(1,pan_id) = addnods(1)
        grid%panel(2,pan_id) = addnods(2)
        grid%panel(3,pan_id) = addnods(4)
        grid%panel(4,pan_id) = addnods(3)
        
        kuttaedges%edge2wpan(iedg) = pan_id 
    end do  ! kutta edges
    
    nwn = nnew ; wake%nwp = pnew   ! set number of wake nods&pans (ONLY! first row at start)
    ! when ends, nnew and pnew is the total new nodes and panels of the wake; nnew+nnod=node_id is the total nodes (body and wake), and nbp+pnew=pan_id is the total number of panels
End Subroutine WAKE_BUFFER
    
Subroutine LE_DETEC ( ) ! Leading edges detection (requires additional mesh data; see .dat file); JCPG
    USE Arrays, only: internal_LE,grid,wake,internal_LE_detected,plate_LE_detected,plate_LE,plate_TE,plate_TE_detected,spanwise_wake_elem,chordwise_wake_elem
    Implicit none
    integer :: i,j,cont
    logical :: log
    
    allocate(internal_LE_detected(grid%nle),plate_LE_detected(grid%plate_nle),plate_TE_detected(grid%plate_nte),spanwise_wake_elem(grid%nle + grid%plate_nle + grid%plate_nte),chordwise_wake_elem(wake%nwp - (grid%nle + grid%plate_nle + grid%plate_nte)) )
    
    do i=1, grid%nle ! over LE 1->12; shown ranges are used for 4x4 discretization example
        do j=grid%nelem+1, grid%nelem + wake%nwp ! 17->56
            if ( (internal_LE(1,i)==grid%panel(2,j)) .and. (internal_LE(2,i)==grid%panel(1,j)) ) then
                internal_LE_detected(i) = j ! j-wake
            else
            end if  
        end do ! ends j
    end do ! ends i
    
    do i=1, grid%plate_nle ! over LE 1->4
        do j=grid%nelem+1, grid%nelem + wake%nwp ! 17->56
            if ( (plate_LE(1,i)==grid%panel(2,j)) .and. (plate_LE(2,i)==grid%panel(1,j)) ) then
                plate_LE_detected(i) = j ! j-wake
            else
            end if  
        end do ! ends j
    end do ! ends i
    
    do i=1, grid%plate_nte ! over TE 1->4 ! mandatory to recognize lateral edges (for Gutnikov's force calculation)
        do j=grid%nelem+1, grid%nelem + wake%nwp ! 17->56
            if ( (plate_TE(1,i)==grid%panel(2,j)) .and. (plate_TE(2,i)==grid%panel(1,j)) ) then
                plate_TE_detected(i) = j ! j-wake
            else
            end if  
        end do ! ends j
    end do ! ends i    
    
    do i=1, grid%nle + grid%plate_nle + grid%plate_nte ! concatenation for Gutnikov's force calculation
        if (i<=grid%nle) then
            spanwise_wake_elem(i) = internal_LE_detected(i)
        else if  (i<=grid%nle + grid%plate_nle) then
            spanwise_wake_elem(i) = plate_LE_detected(i-grid%nle)
        else ! (i<=grid%nle+grid%plate_nle+grid%plate_nte) then
            spanwise_wake_elem(i) = plate_TE_detected(i-grid%nle-grid%plate_nle)
        end if
    end do
    
    cont=0
    do i=grid%nelem+1, grid%nelem+wake%nwp ! 17->56
        log=.false.
        do j=1, grid%nle + grid%plate_nle + grid%plate_nte ! 1->20
            if ( i==spanwise_wake_elem(j) ) then
                log=.true.
                cont=cont+1
            else
            end if
        end do
        if (log==.false.) then
            chordwise_wake_elem(i-grid%nelem-cont) = i ! lateral wake elements (for Gutnikov's force calculation)
        else 
        end if
    end do
    
End Subroutine LE_DETEC

Subroutine V2_STARTING_SOL( circ, n, point, tot_area, op, step ) ! solves the system of equations to obtain a starting solution for the multi-vorton scheme (V2); JCPG
    Use UTILS, only : RESIZE_REAL_1
    Use ARRAYS, only : grid,ctrl_pts,wake,kuttaedges,pr,A_body,nor_vec,A_first,vdir,rhs,A_solv,pivot,bound_circ,A_total,gamma_wake,A_multiwake,esed1_reduc,orien_reduc,internal_LE_detected,plate_LE_detected,fstwake_len,gamma_wake_start,esed1_modif_Gutnikov,del_pres,del_cp,del_force,gamma_wake_mod,input,rhs_freeflow,del_mom,del_pres_down,del_force_down
    Use MATH, only : DOTPRODUCT
    Implicit none
    integer :: i,j,k,trailpan,ok,m
    real(kind=pr) :: vind_ring(3),mid_point(3),equiv_vcr!,sum_circ_plate
    real(kind=pr), intent(inout) :: tot_area
    real(kind=pr), intent(out) :: circ,point(3)
    integer, intent(out) :: n,op,step
    real(kind=pr) :: dot
    logical :: log1,log2

    allocate( A_first(grid%nelem,grid%nelem), A_total(grid%nelem,grid%nelem), A_multiwake(grid%nelem,grid%nelem) )
    allocate(gamma_wake(wake%nwp),esed1_reduc(wake%nwp),orien_reduc(wake%nwp), fstwake_len(4,wake%nwp), gamma_wake_start(wake%nwp),esed1_modif_Gutnikov(4,grid%nelem),del_pres(grid%nelem),del_cp(grid%nelem),del_force(3,grid%nelem),gamma_wake_mod(wake%nwp),del_mom(3,grid%nelem),del_pres_down(grid%nelem),del_force_down(3,grid%nelem) )
    
    ! Multi-wake effect to the A (total) matrix; A_total = A_body + A_multiwake
    A_multiwake(:,:) = 0._pr
    do i = 1, grid%nelem  ! over all control points
        point(:) = ctrl_pts(:,i)
        do j=1, wake%nwp ! wake panels
            circ = 1._pr ! unitary circulation
            n = kuttaedges%edge2wpan(j) ! wake panel numeration (sequential)
            op = 1 ! for wake case
            call VORTONRING ( n,circ,point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),dot)
                
            do k = kuttaedges%esed2(j)+1, kuttaedges%esed2(j+1)  ! matrix entries
                trailpan = kuttaedges%esed1(k) ! body's panel sharing the separation edge
                          
                log1 = .false.; log2 = .false.
                do m=1, grid%nle
                    if ( n==internal_LE_detected(m) ) then ! determines if "n-wake" is an internal LE
                        log1 = .true.
                    else
                    end if
                end do
                do m=1, grid%plate_nle
                    if ( n==plate_LE_detected(m) ) then ! determines if "n-wake" is a plate's LE
                        log2 = .true.
                    else
                    end if
                end do   
                                   
                if ( kuttaedges%orien(k) == .true. ) then ! panel and wake have the same orientation
                    if (log1 == .true.) then ! for internal LEs
                        A_first(i,trailpan) = 0._pr
                    else if (log2 == .true.) then ! for plate's LEs
                        select case(input%detach_model)
                        case(1)
                            A_first(i,trailpan) = dot
                        case(2) ! for positive LE wake (Kutta condition)
                            A_first(i,trailpan) = -dot ! opposite sign for LE
                            !A_first(i,trailpan) = dot ! maintains the same sign during the assembly, but changes it during advection!
                        end select
                    else ! for trailing (and lateral) wakes
                        A_first(i,trailpan) = -dot 
                    end if
                else ! .false., for positive (different orientation) wakes respect its emiting panel
                    A_first(i,trailpan) = dot 
                end if
                    
            A_multiwake(i,trailpan) = A_multiwake(i,trailpan) + A_first(i,trailpan)
            end do  ! ends panels sharing the separation edge
        end do ! ends wake panels
        call DOTPRODUCT( -vdir(:), nor_vec(:,i), rhs(i) ) ! system of equations' RHS (only one wake row)
    end do  ! body panels
    
    rhs_freeflow(:) = rhs(:) ! for Kelvin's condition
    A_total(:,:) = A_body(:,:) + A_multiwake(:,:) ! A total (body + multiwake)
        
    A_solv(:,:) = A_total(:,:) ! to avoid undesired modification to original matrix after solving the system
    
    if (pr .eq. 4) then ! system solution
        call sgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! for single precision
    else
        call dgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! for double precision
    end if 
    
    bound_circ(:) = rhs(:) ! new panels' circulation

    !sum_circ_plate = sum(abs(bound_circ)) ! total plate's circulation
    
    step=0 ! for starting solution case (wake rings)
    call DETACHVORT ( step ) ! calculates the wake circulation strenght between the BVR (internal and external wakes)
    !call V2_START_FORCES_KJ ( circ, mid_point, tot_area ) ! calls to starting (or steady) force calculation
    
End Subroutine V2_STARTING_SOL

Subroutine V2_START_FORCES_KJ ( circ, mid_point, tot_area ) ! hydrodynamic coefficients calculation via Kutta-Zhukovski (KJ) approach for the multi-vorton scheme; JCPG/EOA
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,orien_reduc,internal_LE_detected,plate_LE_detected,gamma_wake_start,x_cm,ctrl_pts,half_chordpan
    Use MATH, only : CROSSPRODUCT
    Implicit none

    integer :: h,j,k,L,fst_nod,sec_nod,op,m
    real(kind=pr) :: sum_mom,den,cfx,cfy,cfz,cl,cd,cy,cn,equiv_vcr,dist,cm
    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta
    real(kind=pr) :: tot_force_steady(3),pan_force_steady(3),vind_body(3),vind_ring(3),ra(3),rb(3),vind_wakeonpoint(3),vind_total(3),del_gamma(3),segm_force(3),tot_force(3)
    real(kind=pr), intent(in) :: tot_area
    real(kind=pr), intent(out) :: mid_point(3),circ
    logical :: log1,log2

    tot_force_steady(:) = 0._pr
    sum_mom = 0._pr
    x_cm = input%x_pos_mom_fact*input%char_len ! moment axis position
    
    select case(input%mesh_num) ! for CM calculation (to displace vertical force component to panel's LE)
    case(1)     ! 2x2
        half_chordpan=0.25_pr
    case(2)     ! 4x4
        half_chordpan=0.125_pr
    case(3,6)   ! 10x10
        half_chordpan=0.05_pr
    case(4,5,7) ! 16x16, 16x8
        half_chordpan=0.03125_pr
    end select 
    
    do k=1, grid%nelem ! over BVR/panels
        pan_force_steady(:) = 0._pr
        
        fst_nod = grid%panel(4,k) ! first node (it depends on mesh numeration's starting point)
        sec_nod = grid%panel(3,k) ! second node (it depends on mesh numeration's starting point)
                
        log1 = .false.; log2 = .false.
        do m=1, grid%nle
            if ( fst_nod==grid%panel(2,internal_LE_detected(m)) .and. sec_nod==grid%panel(1,internal_LE_detected(m))) then
                log1 = .true.
            else
            end if
        end do
        
        do m=1, grid%plate_nle
            if ( fst_nod==grid%panel(2,plate_LE_detected(m)) .and. sec_nod==grid%panel(1,plate_LE_detected(m)) ) then
                select case(input%detach_model)
                case(1)
                    log2 = .true.
                case(2) ! for positive LE wake (Kutta condition)
                    log2 = .false.
                end select
            else
            end if
        end do 
        
        if (log1==.true. .or. log2==.true.) then
            vind_body(:) = 0._pr
            mid_point(:) = ((grid%coord(:,grid%panel(4,k)) + grid%coord(:,grid%panel(3,k)))/2._pr)
        
            do L=1, grid%nelem ! over panels
                circ = bound_circ(L)
                op = 0 ! for plate case
                call VORTONRING ( L,circ,mid_point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
                vind_body(:) = vind_body(:) + vind_ring(:)
            end do !ends L
        
            vind_wakeonpoint(:) = 0._pr
            h=1
            do j = grid%nelem + 1, grid%nelem + wake%nwp ! wake panel number, for 4x4 case: 17->56,96,136...
                circ = gamma_wake_start(j - grid%nelem) ! for 4x4 case: 41->80,120...
                op = 1 ! for wake case
                call VORTONRING ( j,circ,mid_point,vind_ring,op,equiv_vcr ) ! induced velocity by a bounded vorton ring
            
                if ( orien_reduc(h) == .true. ) then ! determines the correct orientation of the detached wakes
                    vind_ring(:) = -vind_ring(:)
                else
                end if

                vind_wakeonpoint(:) = vind_wakeonpoint(:) + vind_ring(:)
                h=h+1
                if (h==wake%nwp+1) then
                    h=1
                else
                end if
            end do ! ends j
            vind_total(:) = vdir(:) + vind_wakeonpoint(:) + vind_body(:) ! total local flow velocity
        
            ! Only LE for FMVLM-based scheme
            ra(:) = grid%coord(:,grid%panel(4,k)) ! bounded segment's first node position
            rb(:) = grid%coord(:,grid%panel(3,k)) ! bounded segment's second node position
         
            del_gamma(:) = (rb(:)-ra(:))*bound_circ(k) ! panel's remaining segments
            call CROSSPRODUCT(vind_total(:), del_gamma(:), segm_force(:)) ! forces on vortex segment
        else
            segm_force(:) = 0._pr ! for Kutta's type separation edge
        end if
        
        segm_force(:) = input%dens*segm_force(:)
        pan_force_steady(:) = pan_force_steady(:) + segm_force(:) ! perpendicular force on the k-panel
        tot_force_steady(:) = tot_force_steady(:) + pan_force_steady(:) ! body's steady force
        
        dist = ctrl_pts(1,k) - x_cm - half_chordpan !-0.03125_pr
        sum_mom = sum_mom + ( pan_force_steady(3) * dist ) !((ctrl_pts(1,k)-(char_len/real(trail_lines))) - x_cm)  )

    end do ! ends k
        
    tot_force(:) = tot_force_steady(:) ! body's total force (only steady)

    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
    cfx = 2._pr*tot_force(1)/den; cfy = 2._pr*tot_force(2)/den; cfz = 2._pr*tot_force(3)/den  ! body axes coefficients
    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !input%beta*(pi/180._pr)  ! angles in radians

    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)

    cl  = cfz*cosalpha - cfx*sinalpha                                ! lift coef.            (CL)
    cd = cfx*cosalpha*cosbeta + cfy*sinbeta + cfz*sinalpha*cosbeta   ! drag coef.            (CD)
    cy  = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef.   (CY)
    cn  = cl*cosalpha + cd*sinalpha                                  ! normal force coef.    (CN)
    cm = -sum_mom/(0.5_pr*input%dens*input%q_inf*input%q_inf*tot_area*input%char_len) ! pitching moment coefficient (is multiplied by -1 to represent its value in the aerodynamical convention; positive: clockwise)
    
400 format(f12.6)
    print 400, cl; print 400, cd; print 400, cm; print 400, cn
    print *,'---'

    500 format(f12.6,f12.6,f12.6,f12.6,f12.6)
    open(4, file='aero_coef.dat', access='append') ! writes output file
    write(4,500) cl, cd, cm, cn
    close(4)

End Subroutine V2_START_FORCES_KJ

Subroutine TIME_STEPS (point, step,circ_vec,pos_ivort,pos_kvort,i,k,tot_area,den_aux,option,vind_total,circ,j,option_2,circ_vec_2,pos_kvort_2) ! unsteady calculation; JCPG
    Use UTILS, only : RESIZE_REAL_1,RESIZE_REAL_2
    Use ARRAYS, only : grid,input,wake,pr,vdir,bound_circ,fouroverthree,pi,ctrl_pts,nor_vec,rhs,bound_circ_old,a_body,a_solv,pivot,vel_fst_point_old,vel_sec_point_old,vind_body,vind_cloud,oneoverthree,vorticity,vorticity_old,sum_circ_plate,seg2vortchild,elem2seg,sum_circ_plate_total,vxcore_tube_var,dL,dL_old,vort_mag,vort_old_mag,fourpi,orien_reduc_steps,orien_reduc,vel_vorton_avg,third_term,fourth_term,second_term,second_term_down,fourth_term_down
    Use MATH, only : DOTPRODUCT,VECTORNORM
    Use OMP_LIB
    Implicit none
    integer :: ok,s,cont,cont2
    integer, intent(out) :: step,i,k,den_aux,option,j,option_2
    real(kind=pr) :: vind_v2p(3),vind_total_fst_point(3),vind_total_sec_point(3),vind_v2p_2(3),dlen_dt_vec(3),del_vort(3),volume
    real(kind=pr) :: dlen,vxcore_rad3,gammav_size,dif_circ_total,volume_cyl,vort_visc,zero_vort
    real(kind=pr), intent(out) :: point(3),circ_vec(3),pos_ivort(3),pos_kvort(3),circ,circ_vec_2(3),pos_kvort_2(3)
    real(kind=pr), intent(inout) :: vind_total(3),tot_area

    allocate ( wake%pos_old(3,wake%nwp*input%nsteps),wake%pos_old_fst(3,wake%nwp*input%nsteps),wake%pos_old_sec(3,wake%nwp*input%nsteps) )
    allocate ( wake%pos(3,wake%nwp*(input%nsteps+1)),wake%gammav(3,wake%nwp*(input%nsteps+1)),wake%r_dir(3,wake%nwp*input%nsteps),wake%gamma_mag(wake%nwp*(input%nsteps+1)),wake%volume(wake%nwp*(input%nsteps+1)),wake%vxcore_rad(wake%nwp*(input%nsteps+1)),wake%vxcore_tube(wake%nwp*(input%nsteps+1)),wake%vxcore_rad_old(wake%nwp*input%nsteps),wake%vxcore_tube_old(wake%nwp*input%nsteps),wake%pos_multifix2plate(3,wake%nwp),wake%multicirc_vec(3,wake%nwp),wake%multivxcore_rad(wake%nwp),wake%partxfil(wake%nwp),wake%multicirc_mag(wake%nwp),wake%gammav_Gutnikov(3,4),seg2vortchild(wake%nwp),elem2seg(2,wake%nwp) )
    allocate ( wake%pos_fst_point(3,wake%nwp*(input%nsteps+1)),wake%pos_sec_point(3,wake%nwp*(input%nsteps+1)),wake%pos_fst_point_modif(3,wake%nwp*(input%nsteps+1)),wake%pos_sec_point_modif(3,wake%nwp*(input%nsteps+1)) )
    allocate ( vel_fst_point_old(3,wake%nwp*input%nsteps),vel_sec_point_old(3,wake%nwp*input%nsteps),wake%length(wake%nwp*(input%nsteps+1)),wake%length_old(wake%nwp*input%nsteps),vorticity_old(3,wake%nwp*(input%nsteps+1)),vorticity(3,wake%nwp*(input%nsteps+1)),vxcore_tube_var(wake%nwp*(input%nsteps+1)),dL(3,wake%nwp*(input%nsteps+1)),dL_old(3,wake%nwp*input%nsteps) )
    allocate (orien_reduc_steps(wake%nwp*(input%nsteps+1)))
    ! new allocations
    allocate ( vel_vorton_avg(3,wake%nwp*input%nsteps),third_term(grid%nelem),fourth_term(grid%nelem),second_term(grid%nelem),fourth_term_down(grid%nelem),second_term_down(grid%nelem) )
    
    dif_circ_total = 0._pr
    sum_circ_plate_total = 0._pr
    orien_reduc_steps(:) = .false.

    do step=1, input%nsteps ! main unsteady loop
        wake%vxcore_rad(1:wake%nwp) = input%core_rad_init ! bounded vortons' core radius (constant)
        wake%volume(1:wake%nwp) = fouroverthree*pi*input%core_rad_init*input%core_rad_init*input%core_rad_init  ! constant vortons' volume
        call NEW_VORTONS( s,cont,den_aux,step,cont2 ) ! generates the nascent vortons over the plate at a prescribed distance (epsilon) 

        circ=0._pr;point(1:3)=0._pr        

        wake%num_vort = wake%nwp*step ! number of vortons up the current time step; 40,80,120,... (4x4 example)
        
        if (step==1) then
            wake%num_vort = wake%nwp
        else
        end if
        
        ! ADVECTION STEP
        wake%pos_old(:,:) = wake%pos(:,:) ! copies the previous vortons' positions to convect them
        wake%pos_old_fst(:,:) = wake%pos_fst_point(:,:) ! copies the previous vortons first point's (hypothetical vortex filament) positions to convect them
        wake%pos_old_sec(:,:) = wake%pos_sec_point(:,:) ! copies the previous vortons second point's (hypothetical vortex filament) positions to convect them
        
        if (step==1) then ! only for the first iteration (step=1); Euler (one step) advection as (STRICKLAND, 2002)
            ! FOR FILAMENT'S FIRST POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,vdir,input) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,vind_body,point,option,vind_total_fst_point,vind_total,vel_fst_point_old) !& ! (declaration of private variables)
            !!!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons (inviscid particles at this point)
                pos_ivort(:) = wake%pos_fst_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body/plate
                point(:) = pos_ivort(:) ! vorton location
                option=0 ! for plate case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do  
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                
                vind_total_fst_point(:) = vind_body(:) + vdir(:)
                vind_total(:) = vind_total_fst_point(:)
                
                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons
                else
                end if
                
                wake%pos_fst_point(:,i) = wake%pos_old_fst(:,i) + vind_total(:)*input%fst_wake_factor*input%dt ! one step forward convection
                vel_fst_point_old(:,i) = vind_total(:) ! total velocity over i-vorton (for first AB calculation at step=2)
            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            
            ! FOR FILAMENT'S SECOND POINT (separated to avoid to use 3-dimensional arrays)
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,den_aux,vdir,input,vel_sec_point_old) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,point,option,k,circ_vec,pos_kvort,vind_v2p,vind_total_sec_point,vind_total) & ! (declaration of private variables)
            !!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons
                pos_ivort(:) = wake%pos_sec_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body/plate
                point(:) = pos_ivort(:) ! vorton's location
                option=0 ! for plate case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
               
                vind_total_sec_point(:) = vind_body(:) + vdir(:)
                vind_total(:) = vind_total_sec_point(:)
                
                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons
                else
                end if
                
                wake%pos_sec_point(:,i) = wake%pos_old_sec(:,i) + vind_total(:)*input%fst_wake_factor*input%dt ! one step forward convection
                vel_sec_point_old(:,i) = vind_total(:) ! total velocity over i-vorton (for first AB calculation at step=2)
                
                vel_vorton_avg(:,i) = (vel_fst_point_old(:,i) + vel_sec_point_old(:,i)) / 2._pr ! average velocity (on vorton); NEW LINE for pressure calculation
            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
                        
        else ! (step>1) ! 2nd order Adams-Bashforth convection as in (STRICKLAND, 2002).
            ! FOR FILAMENT'S FIRST POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,vdir,input,vind_body,vind_cloud,vel_fst_point_old,vind_total) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,point,option,option_2,vind_total_fst_point) !& ! (declaration of private variables)
            !!!$OMP REDUCTION(+: vind_body,vind_cloud) ! (declaration of reduction variables)
            !!$OMP DO ! the issue to parallelize this DO seems to be related to ADABASH2_TRAJ_FST subroutine, where wake%pos_fst_point(:,i) must be declared as REDUCTION type variable...(enable omp_set_nested...)
            do i=1, wake%num_vort ! over the number of vortons
                point(:) = wake%pos_fst_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body/plate
                option=0 ! for plate case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux ! for multi-vortons per vortex filament case
                    circ_vec(:) = wake%multicirc_vec(:,k) ! k-child vorton circulation
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-child vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_cloud(:) = 0._pr ! clears old values
                option_2=1 ! for wake case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(i,point,option_2,WAKE) PRIVATE(j,vind_v2p_2,circ_vec_2,pos_kvort_2) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
                !$OMP DO
                do j = wake%nwp+1, wake%num_vort ! over the free vortons (avoids nascent ones)
                    !if (j/=i) then ! avoids to induce velocity over itself (it does not allow to obtain a div-free grid)
                        circ_vec_2(:) = wake%gammav(:,j) ! vorton's circulation
                        pos_kvort_2(:) = wake%pos(:,j)
                        call VEL_VORT2POINT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,j,option_2) ! calculates the induced velocity by a vorton over another one
                        !call VEL_VORT2POINTVORT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,i,j,option_2) ! calculates the induced velocity by a vorton over another one
                        vind_cloud(:) = vind_cloud(:) + vind_v2p_2(:) ! induced velocity by the vorton cloud over a single one
                    !else
                    !end if
                end do ! ends j
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                    
                vind_total_fst_point(:) = vind_body(:) + vind_cloud(:) + vdir(:)
                vind_total(:) = vind_total_fst_point(:)

                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING ( i,vind_total ) ! forces to tangential velocity (by nulifying the negative normal velocity component) ONLY to nascent vortons
                else
                end if
                call ADABASH2_TRAJ_FST (i,vind_total) ! trajectory integration
                vel_fst_point_old(:,i) = vind_total(:) ! for next AB calculation 
            end do ! ends i
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            wake%pos_fst_point(:,:) = wake%pos_fst_point_modif(:,:) ! copies new modified positions for the next iterations
                        
            ! FOR FILAMENT'S SECOND POINT 
            
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(wake,den_aux,vdir,input,vel_sec_point_old) & ! (declaration of shared variables)
            !!$OMP PRIVATE(i,pos_ivort,point,option,option_2,k,j,circ_vec,circ_vec_2,pos_kvort,pos_kvort_2,vind_v2p,vind_v2p_2,vind_total_sec_point,vind_total) & ! (declaration of private variables)
            !!$OMP REDUCTION(+: vind_body,vind_cloud) ! (declaration of reduction variables)
            !!$OMP DO
            do i=1, wake%num_vort ! over the number of vortons
                point(:) = wake%pos_sec_point(:,i) ! position vector of the i-vorton
                vind_body(:) = 0._pr ! clears old values; induced velocity by the body/plate
                option=0 ! for plate case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(den_aux,i,point,option,WAKE) PRIVATE(k,vind_v2p,circ_vec,pos_kvort) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
                !$OMP DO
                do k=1, den_aux ! over the fixed vortons on plate
                    circ_vec(:) = wake%multicirc_vec(:,k)
                    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                    !call VEL_VORT2POINTVORT(circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one
                    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
                end do 
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_cloud(:) = 0._pr ! clears old values
                option_2=1 ! for wake case
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(i,point,option_2,WAKE) PRIVATE(j,vind_v2p_2,circ_vec_2,pos_kvort_2) & ! (declaration of variables)
                !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
                !$OMP DO
                do j = wake%nwp+1, wake%num_vort ! over the free vortons (avoids nascent ones)
                    !if (j/=i) then ! avoids to induce velocity over itself (it does not allow to obtain a div-free grid)
                        circ_vec_2(:) = wake%gammav(:,j) ! vorton's circulation
                        pos_kvort_2(:) = wake%pos(:,j)
                        call VEL_VORT2POINT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,j,option_2) ! calculates the induced velocity by a vorton over another one
                        !call VEL_VORT2POINTVORT(circ_vec_2,point,pos_kvort_2,vind_v2p_2,i,j,option_2) ! calculates the induced velocity by a vorton over another one
                        vind_cloud(:) = vind_cloud(:) + vind_v2p_2(:) ! induced velocity by the vorton cloud over a single one
                    !else
                    !end if
                end do ! ends j
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                vind_total_sec_point(:) = vind_body(:) + vind_cloud(:) + vdir(:)
                vind_total(:) = vind_total_sec_point(:)

                if (i<=wake%nwp) then
                    call TANG_VEL_FORCING ( i,vind_total ) ! forces to tangential velocity (by nulifying the normal velocity component) ONLY to nascent vortons. Is it physically justifiable?
                else
                end if
                    
                call ADABASH2_TRAJ_SEC (i,vind_total) ! trajectory integration
                vel_sec_point_old(:,i) = vind_total(:) ! for next AB calculation 
                
                vel_vorton_avg(:,i) = (vel_fst_point_old(:,i) + vel_sec_point_old(:,i)) / 2._pr ! average velocity (on vorton); NEW LINE for pressure calculation
                
            end do ! ends i 
            !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
            wake%pos_sec_point(:,:) = wake%pos_sec_point_modif(:,:) ! copies new modified positions for the next iteration
        end if
        
        orien_reduc_steps(1:wake%nwp) = orien_reduc(:)
        
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(wake,orien_reduc_steps,dL_old) PRIVATE(i,dL) ! (declaration of variables)
            !$OMP DO        
            do i=1, wake%num_vort ! over the number of vortons
                if (orien_reduc_steps(i) == .true.) then
                    dL_old(:,i) = wake%pos_old_fst(:,i) - wake%pos_old_sec(:,i) ! vectorial length
                    dL(:,i) = wake%pos_fst_point(:,i) - wake%pos_sec_point(:,i) ! vectorial length
                    wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_sec_point(:,i) ! new vorton position due vortex strain
                else !(orien_reduc_steps(i) == .false.)
                    dL_old(:,i) = wake%pos_old_sec(:,i) - wake%pos_old_fst(:,i) ! vectorial length
                    dL(:,i) = wake%pos_sec_point(:,i) - wake%pos_fst_point(:,i) ! vectorial length
                    wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_fst_point(:,i) ! new vorton position due vortex strain
                end if
            end do
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL            
            
            ! SIMPLIFIED WAKE PENETRATION AVOIDANCE (only valid for flat plate cases)
            call POINT_POLYGON (step) ! comment it to maintain a div-free grid (scheme 1: vorticity crosses the plate without restrictions)
        
           ! VORTEX STRETCHING/SQUEEZING CALCULATION
            wake%length_old(:) = wake%length(:) ! previous vortex filament-vorton length
            wake%vxcore_rad_old(:) = wake%vxcore_rad(:) ! previous vortex core radius
            vorticity_old(:,:) = vorticity(:,:)        
                !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
                !$OMP SHARED(wake,orien_reduc_steps,input,dL_old,step,vorticity_old,vorticity) PRIVATE(zero_vort,vort_visc,volume,volume_cyl,i,dL,dlen_dt_vec,dlen,vxcore_rad3,gammav_size,del_vort,vort_mag,vort_old_mag) & ! (declaration of variables)
                !$OMP REDUCTION(+: vxcore_tube_var) ! (declaration of reduction variables)
                !$OMP DO
                do i=1, wake%num_vort ! over the number of current wake vortons
                    if (orien_reduc_steps(i) == .true.) then
                        dL(:,i) = wake%pos_fst_point(:,i) - wake%pos_sec_point(:,i) ! vectorial length
                        wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_sec_point(:,i) ! new vorton position due vortex strain
                    else ! (orien_reduc_steps(i) == .false.)
                        dL(:,i) = wake%pos_sec_point(:,i) - wake%pos_fst_point(:,i) ! vectorial length
                        wake%pos(:,i) = ( dL(:,i)/2._pr ) + wake%pos_fst_point(:,i) ! new vorton position due vortex strain
                    end if
                
                    call VECTORNORM(dL(:,i),wake%length(i)) ! filament's length (for next iteration)
                    
                    ! NEW VOLUME DUE TO VORTEX FILAMENT/TUBE STRETCHING (inviscid and viscous scheme)
                    ! In this section, the equation numbers correspond to the referenced in the main publication (PIMENTEL, 2023)
                    dlen_dt_vec(:) = ( dL(:,i) - dL_old(:,i) ) / input%dt ! EQ. 9; position vector's time derivative
                    dlen = wake%length(i) - wake%length_old(i) ! filament length variation
                    vxcore_rad3 = wake%vxcore_rad(i)*wake%vxcore_rad(i)*wake%vxcore_rad(i) ! cubic vortex core radius (sphere)
                    wake%vxcore_tube_old(i) = sqrt(fouroverthree*(vxcore_rad3/wake%length_old(i))) ! EQ. 10; equivalent vortex filament radius; vol_vorton=vol_tube
                
                    call VECTORNORM(wake%gammav(:,i),gammav_size) ! vortex filament scalar circulation
                    volume_cyl = pi*wake%vxcore_tube_old(i)*wake%vxcore_tube_old(i)*wake%length_old(i) ! EQ. 11; new vortex tube's volume
                    del_vort(:) = ( gammav_size/volume_cyl ) * dlen_dt_vec(:) ! EQ. 12; vorticity vector's time derivative
                
                    if (dlen >= 0._pr) then
                        vorticity(:,i) = vorticity(:,i) + del_vort(:)*input%dt  ! EQ. 13a; vorticity variation due to vortex stretching (inviscid); (STRICKLAND, 2002; eq. 56)
                    else !(dlen<0._pr)
                        vorticity(:,i) = vorticity(:,i) - del_vort(:)*input%dt  ! EQ. 13b; vorticity variation due to vortex squeezing (inviscid); (STRICKLAND, 2002; eq. 56)
                    end if     
           
                    call VECTORNORM(vorticity(:,i),vort_mag) ! vorticity's magnitude
                    call VECTORNORM(vorticity_old(:,i),vort_old_mag) ! old vorticity's magnitude
                    
                    if (vort_mag >= input%tol_vort) then ! avoids dividing by zero (NaN); higher values for tol_vort can be selected to avoid instability
                        select case (input%vx_stretch)
                        case (0) ! constant volumes' scheme
                            volume = wake%volume(i) ! the vortex tube volume remains the same according to (STRICKLAND, 2002)
                            wake%vxcore_tube(i) = wake%vxcore_tube_old(i) * sqrt(wake%length_old(i)/wake%length(i)) ! corresponding vortex tube core radius
                            vxcore_tube_var(i) = (wake%vxcore_tube(i) - wake%vxcore_tube_old(i)) / input%dt ! vortex tube's core variation (dsigma/dt) due to pure advection
                        case (1) ! variable volumes' scheme
                            volume = (vort_old_mag/vort_mag)*wake%volume(i) ! EQ. 14; new vortex tube's volume after pure advection (inviscid)
                            wake%vxcore_tube(i) = sqrt(volume/(pi*wake%length(i))) ! EQ. 15; new vortex tube's core radius
                            vxcore_tube_var(i) = (wake%vxcore_tube(i) - wake%vxcore_tube_old(i)) / input%dt ! EQ. 16; vortex tube's core radius variation due to pure advection
                        end select
                                                
                        select case(input%regul_function)
                        case (1) ! for second-order Gaussian regularization function (STRICKLAND, 2002)
                            vxcore_tube_var(i) = vxcore_tube_var(i) + ((2._pr*input%kin_visc)/wake%vxcore_tube_old(i)) ! EQ. 17 with k=2; vxcore_tube_old instead vxcore_tube (BARBA, 2004)
                        case (2) ! for Gaussian error function (ALVAREZ, 2022)
                            vxcore_tube_var(i) = vxcore_tube_var(i) + (input%kin_visc/wake%vxcore_tube_old(i)) ! EQ. 17 with k=1
                        end select
                        
                        wake%vxcore_tube(i) = wake%vxcore_tube_old(i) + vxcore_tube_var(i)*input%dt ! EQ. 18; new vortex tube's core radius after pure advection and viscous diffusion
                        wake%volume(i) = pi*wake%vxcore_tube(i)*wake%vxcore_tube(i)*wake%length(i) ! EQ. 19; new volume due to pure advection and viscous diffusion
                        wake%vxcore_rad(i) = ((3._pr*wake%volume(i))/fourpi)**oneoverthree ! EQ. 20; new vorton's core radius
                        vort_visc = vort_mag*(volume/wake%volume(i)) ! EQ. 21; new vorticity magnitude
                        vorticity(:,i) = vorticity(:,i)*(vort_visc/vort_mag) ! EQ. 22; new vorticity vector
                        zero_vort = vort_old_mag*volume_cyl - vort_visc*wake%volume(i) ! for testing; it is non-zero for constant volumes' scheme and 'numerically zero' for variable volumes' one
                    else
                    end if
                end do
                !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
                !$OMP end PARALLEL
                
        ! COPYING (MOVING DATA INTO ARRAYS) FOR NEXT ITERATION
        if (step<=input%nsteps) then 
            !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !!$OMP SHARED(step,wake) PRIVATE(i) ! (declaration of variables)
            !!$OMP DO
                !!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(wake, vorticity, vorticity_old, orien_reduc_steps, step) ! COPILOT
            do i=(step+1)*wake%nwp, wake%nwp+1, -1 ! wake vorton's number ! 80 -> 41
                wake%pos(:,i) = wake%pos(:,i-wake%nwp) ! vortons positions
                wake%pos_fst_point(:,i) = wake%pos_fst_point(:,i-wake%nwp) ! first point hypothetical filament position
                wake%pos_sec_point(:,i) = wake%pos_sec_point(:,i-wake%nwp) ! second point hypothetical filament position
                wake%length(i) = wake%length(i-wake%nwp) ! vortex filament length
                wake%gammav(:,i) = wake%gammav(:,i-wake%nwp) ! copies first layers vortons' circulations to the following one
                wake%vxcore_rad(i) = wake%vxcore_rad(i-wake%nwp) ! vortex core radius
                wake%vxcore_tube(i) = wake%vxcore_tube(i-wake%nwp)
                wake%volume(i) = wake%volume(i-wake%nwp) ! vorton volumes
                wake%gamma_mag(i) = wake%gamma_mag(i-wake%nwp) ! circulation strengths
                vorticity(:,i) = vorticity(:,i-wake%nwp)
                vorticity_old(:,i) = vorticity_old(:,i-wake%nwp)
                orien_reduc_steps(i) = orien_reduc_steps(i-wake%nwp) ! for correct direction of dL (vortex stretching) 
            end do
                !!$OMP END PARALLEL DO ! COPILOT
            !!$OMP END DO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !!$OMP end PARALLEL
        else; end if 
        
        ! VORTEX MERGING (optional future improvement)
        ! VISIBILITY (optional future improvement to avoid vorticity crossing/penetration for more complex geometries; it have been applied for the bidimensional case)
            
        ! FREE VORTONS CONTRIBUTION TO RHS
        do i = 1, grid%nelem  ! over all control points
            vind_cloud(:) = 0._pr
            point(:) = ctrl_pts(:,i)
            option=1 ! for wake case
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(wake,step,point,option,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
            !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
            !$OMP DO
            do k = wake%nwp+1, wake%nwp*(step+1) ! from 41 to 80 (4x4 example) 
                circ_vec(:) = wake%gammav(:,k) ! k-vorton's vectorial circulation
                pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over the control point
                vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over the control point
            end do ! end wake panels
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL
            call DOTPRODUCT( -(vdir(:) + vind_cloud(:)), nor_vec(:,i), rhs(i) ) ! new (unsteady) RHS
                !rhs(i) = dot_product( -(vdir(:) + vind_cloud(:)), nor_vec(:,i) ) ! same as the previous line
        end do ! ends control points 
      
        bound_circ_old(:) = bound_circ(:) ! for force calculation
        
        ! SOLVES THE SYSTEM OF EQUATIONS 
        A_solv(:,:) = A_body(:,:) ! to avoid A_body matrix modification after solving the system

        if (pr .eq. 4) then ! system solution
            call sgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! single precision
        else
            call dgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! double precision
        end if 
        bound_circ(:) = rhs(:) ! new panels' circulation
        
        !sum_circ_plate = sum(abs(bound_circ)) ! total plate's circulation
    
        call DETACHVORT ( step ) ! Calculates the new wake circulation strength between the bounded vortex rings (internal and external wakes)

        print *, step
        print *,'---'
        
        ! FORCE CALCULATION
            select case (input%detach_model)
        case(1) ! pre-stall
            call FORCES_KJ ( circ_vec, point, tot_area, pos_kvort, den_aux, option )
        case(2) ! post-stall (it should work the same way for pre-stall by a more complex vortons scheme (detached from upside and downside faces of the plate))
            call FORCES_PRES ( circ_vec, point, tot_area, pos_kvort, den_aux, option )
        end select
            !call FORCES_GUTNIKOV ( point, tot_area, circ_vec, den_aux, option, pos_kvort) ! abandoned method

        ! SMOOTHING PRESSURES (optional for a future improvement to the pressure distribution calculation)
        
        ! WRITES OUTPUT FILES (Paraview); see the output_vtk folder for State and Windows Arrangement loads. 
        call OUTPUT_PLATE (step,den_aux)
        call OUTPUT_WAKE (step)
        call OUTPUT_FIL (step)
        !!call OUTPUT_TRAJ (step)

        ! WRITES OUTPUT FILES (GiD's postprocessor; it needs license for more than 10k elements)
        !call OUTPUT_MSH (step,den_aux) ! mesh file
        !call OUTPUT_RES (step,den_aux) ! results file
        
    end do ! end steps   
    pause

End Subroutine TIME_STEPS
        
!!! SECONDARY SUBROUTINES    
Subroutine NEW_VORTONS ( i,cont,den_aux,step,cont2 ) ! generates the nascent vortons at a prescribed distance (epsilon) from the plate; JCPG
    Use ARRAYS, only : pr,grid,wake,esed1_reduc,nor_vec,input,orien_reduc,gamma_wake,ds,pi,fouroverthree,vorticity,gamma_wake_mod,esed1_modif_Gutnikov,elem2seg,sum_circ_plate_step,sum_circ_plate_total
    USE MATH, only : VECTORNORM,NORMALIZE
    implicit none
    integer :: j,k
    integer,intent(in) :: step
    integer, intent(out) :: i
    integer, intent(inout) :: cont,den_aux,cont2
    
    if (step==1) then
        cont=1; den_aux = 0
    else
    end if
    cont2=1
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,grid,input,esed1_reduc,nor_vec,orien_reduc,gamma_wake,gamma_wake_start,step,den_aux,cont,cont2) &
    !!$OMP PRIVATE(i,ds,vorticity) ! (declaration of variables)
    !!$OMP DO
    do i=1, wake%nwp ! over all detached wakes from surface
        ds(:) = grid%coord(:,grid%panel(1,i+grid%nelem)) - grid%coord(:,grid%panel(2,i+grid%nelem)) ! segment's vector
        call VECTORNORM(ds,wake%length(i)) ! segment's distance
        
        wake%mid_segm(:) = ( ds(:)/2._pr ) + grid%coord(:,grid%panel(2,i+grid%nelem)) ! segment's midpoint
        wake%pos(:,i) = wake%mid_segm(:) + input%eps*nor_vec(:,esed1_reduc(i)) ! i-nascent vortons' positions
        
        wake%pos_fst_point(:,i) = grid%coord(:,grid%panel(1,i+grid%nelem)) + input%eps*nor_vec(:,esed1_reduc(i)) ! i-nascent vortons' first node position
        wake%pos_sec_point(:,i) = grid%coord(:,grid%panel(2,i+grid%nelem)) + input%eps*nor_vec(:,esed1_reduc(i)) ! i-nascent vortons' second node position
        
        if (orien_reduc(i)==.true.) then ! assigns the orientation to the vorton (depending on merged wake direction)
            wake%r_dir(:,i) = ds(:) ! maintains the direction (+)
        else ! orien_reduc(i)==.false.
            wake%r_dir(:,i) = -ds(:) ! inverts the direction (-)
        end if
            
        select case(input%detach_model)
            case(1) ! negative LE wake
                wake%gammav(:,i) = gamma_wake_mod(i)*wake%r_dir(:,i) ! circulation vector
            case(2) ! positive LE wake
                wake%gammav(:,i) = gamma_wake(i)*wake%r_dir(:,i) ! circulation vector
        end select
                    
        if (step==1) then
            wake%volume(i) = fouroverthree*pi*input%core_rad_init*input%core_rad_init*input%core_rad_init ! constant vorton volumes
            vorticity(:,i) = wake%gammav(:,i) / wake%volume(i) ! vorticity vector; omega vector
            call MULTI_VORTONS_FIL ( i, cont,den_aux ) ! calculates the number and sizes of multi-vortons per vortex filament (only once cause it remains the same along the entire simulation)
        else
        end if
        
        call MULTI_VORTONS_CIRC ( i, cont2 ) ! calculates the bounded vorton circulation strengths (changes their values at each time step according to the previous solution)
    end do ! end wakes
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    
    elem2seg(:,:) = 0
    if (step==1) then ! for Gutnikov's force calculation
        do i=1,grid%nelem ! over number of elements/panels
            do j=1,4 ! over 4 edges per element/panel
                do k=1,wake%nwp ! over detachment segments
                    if (esed1_modif_Gutnikov(j,i)==k .and. elem2seg(1,k)==0) then
                        elem2seg(1,k) = i ! puts on first space the i-panel corresponding to k-segment
                    else if (esed1_modif_Gutnikov(j,i)==k .and. elem2seg(1,k)/=0) then ! for internal/shared segments
                        elem2seg(2,k) = i ! puts on second space the i-panel corresponding to k-segment
                    else
                    end if
                end do ! ends 
            end do ! ends 
        end do ! ends 
    else
    end if
    
    sum_circ_plate_step = 0._pr
    do i=1, wake%nwp ! for conservation of circulation (Kelvin's theorem) verification
        call VECTORNORM( wake%gammav(:,i), wake%gamma_mag(i) ) ! vortex circulation strength; Gamma
        sum_circ_plate_step = sum_circ_plate_step + wake%gamma_mag(i) ! plate's circulation per time step
    end do
    
    sum_circ_plate_total = sum_circ_plate_total + sum_circ_plate_step ! total (all steps) plate's circulations
    
End Subroutine NEW_VORTONS
    
Subroutine MULTI_VORTONS_FIL ( i,cont,den_aux ) ! calculates the number of vortons per vortex filament (based on the nominal vortex core radius; sigma_zero) for the unsteady algorithm; JCPG
    Use ARRAYS, only : pr,input,wake,ds,grid,seg2vortchild
    Use MATH, only : VECTORNORM
    Use UTILS, only : RESIZE_REAL_2,RESIZE_REAL_1,RESIZE_INT_1
    Implicit none
    real(kind=pr) :: del_pos(3),den,pos_kvort(3),equiv_vcr
    integer, intent(in) :: i
    integer, intent(inout) :: den_aux,cont
    integer :: vps
        
    den = real(ceiling(wake%length(i)/input%core_rad_init)) + 1._pr ! denominator: number of vortons per segment
        !den=1._pr ! for testing with only one vorton per segment
    wake%partxfil(i) = den
    del_pos(:) = ds(:) / den  ! delta position 
    equiv_vcr = input%core_rad_init / ( den**(1._pr/3._pr) ) ! equivalent vortex core radius (to maintain a constant volume)
    
    den_aux = den_aux + int(den) ! auxiliar denominator to resize matrices (at last is the total number of fixed vortons; i.e. 160 for 4x4 discretization: 4 vortons x 40 vortex filaments/detached vortons)
    call RESIZE_REAL_2 (3 , den_aux, wake%pos_multifix2plate) ! resizes vorton positions' matrix
    call RESIZE_REAL_2 (3 , den_aux, wake%multicirc_vec) ! resizes circulation's matrix (here instead MULTI_VORTONS_CIRC to save computation)
    call RESIZE_REAL_1 (den_aux, wake%multivxcore_rad) ! resizes circulation's matrix
    call RESIZE_REAL_1 (den_aux, wake%multicirc_mag)
    call RESIZE_INT_1 (den_aux, seg2vortchild)
    
    do vps=1, int(den) ! vortons per segment
        pos_kvort(:) = grid%coord(:,grid%panel(2,i+grid%nelem)) + del_pos(:)/2._pr + (vps-1)*del_pos(:) ! extra k-vorton position
        wake%pos_multifix2plate(:,cont) = pos_kvort(:) ! child vorton position
        seg2vortchild(cont) = i
        wake%multivxcore_rad(cont) = equiv_vcr ! child vorton vortex core radius
        cont=cont+1 ! allows to save positions serially into matrices
    end do

End subroutine MULTI_VORTONS_FIL 
    
Subroutine MULTI_VORTONS_CIRC ( i,cont2 ) ! calculates the circulation strength for the child vortons; JCPG
    Use ARRAYS, only : pr,wake,gamma_wake
    Use MATH, only : VECTORNORM
    Use UTILS, only : RESIZE_REAL_2,RESIZE_REAL_1
    Implicit none
    real(kind=pr) :: den,fac(3)
    integer, intent(in) :: i
    integer, intent(inout) :: cont2
    integer :: vps
        
    den = wake%partxfil(i) ! denominator: number of vortons per segment
    fac(:) = wake%r_dir(:,i) / den  ! delta position 
    
    do vps=1, int(den) ! vortons per segment
        wake%multicirc_vec(:,cont2) = gamma_wake(i)*fac(:) ! child vorton circulation vector
        call VECTORNORM(wake%multicirc_vec(:,cont2), wake%multicirc_mag(cont2))
        cont2=cont2+1 ! allows to save positions serially into matrices
    end do

    End subroutine MULTI_VORTONS_CIRC
    
Subroutine TANG_VEL_FORCING ( i,vind_total ) ! eliminates the vertical component of the velocity --> forces it to be tagentially to the surface; JCPG
    Use ARRAYS, only : pr,nor_vec,esed1_reduc
    Use MATH, only : DOTPRODUCT
    Implicit none
    integer, intent(in) :: i
    real(kind=pr) :: vn
    real(kind=pr), intent(inout) :: vind_total(3)
        
    call DOTPRODUCT(vind_total(:),nor_vec(:,esed1_reduc(i)),vn) ! normal velocity component (on body axes)
    if (vn<0._pr) then ! normal velocity component is negative ("pushes" towards the plate) 
        vind_total(:) = vind_total(:) - nor_vec(:,esed1_reduc(i))*vn ! puts normal velocity's component to zero
    else
    end if

End Subroutine TANG_VEL_FORCING
    
Subroutine ADABASH2_TRAJ_FST (i,vind_total) ! second order Adams-Bashforth trajectory integration for the first filament's point; JCPG
    Use ARRAYS, only : pr,wake,input,vel_fst_point_old,vel_now
    implicit none
    integer, intent(in) :: i
    real(kind=pr), intent(inout) :: vind_total(3)

    vel_now(:) = vind_total(:) ! current total velocity

    if (i>wake%nwp) then ! for free vortons
        wake%pos_fst_point_modif(:,i) = wake%pos_old_fst(:,i) + (input%dt/2._pr)*(3._pr*vel_now(:) - vel_fst_point_old(:,i)) ! first point vortex filament's position
    else ! for first wake row
        wake%pos_fst_point_modif(:,i) = wake%pos_old_fst(:,i) + ( (input%fst_wake_factor*input%dt) / 2._pr)*(3._pr*vel_now(:) - vel_fst_point_old(:,i)) 
    end if
End Subroutine ADABASH2_TRAJ_FST

Subroutine ADABASH2_TRAJ_SEC (i,vind_total) ! second order Adams-Bashforth trajectory integration for the second filament's point; JCPG
    Use ARRAYS, only : pr,wake,input,vel_sec_point_old,vel_now
    implicit none
    integer, intent(in) :: i
    real(kind=pr), intent(inout) :: vind_total(3)
    
    vel_now(:) = vind_total(:) ! current total velocity

    if (i>wake%nwp) then
        wake%pos_sec_point_modif(:,i) = wake%pos_old_sec(:,i) + (input%dt/2._pr)*(3._pr*vel_now(:) - vel_sec_point_old(:,i)) ! first point vortex filament's position
    else ! for first wake row
        wake%pos_sec_point_modif(:,i) = wake%pos_old_sec(:,i) + ( (input%fst_wake_factor*input%dt) / 2._pr)*(3._pr*vel_now(:) - vel_sec_point_old(:,i)) 
    end if
End Subroutine ADABASH2_TRAJ_SEC
            
Subroutine VEL_VORT2POINT (circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over a point; JCPG
    Use ARRAYS, only : pr,input,wake,fiveovertwo,fourpi,fouroverpi,e_num,twopi,pi,twooverpi
    Use MATH, only : VECTORNORM,CROSSPRODUCT,NORMALIZE
    implicit none
    real(kind=pr) :: r_mag,r_mag2,r_mag3,fac1,fac2,fac3,cons,g_winck,rho,rho2,rho3!,sym_core_rad
    real(kind=pr) :: r(3),winck(3)
    real(kind=pr), intent(in) :: circ_vec(3),point(3),pos_kvort(3)
    real(kind=pr), intent(out) :: vind_v2p(3)
    integer, intent(in) :: k,option

    r(:) = point(:) - pos_kvort(:) ! radio vector between k-vorton and an evaluation point
    CALL VECTORNORM(r(:), r_mag) ! vector norm (distance between points)
    
    if (option==0) then ! for plate
        wake%core_rad2 = wake%multivxcore_rad(k)*wake%multivxcore_rad(k) ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*wake%multivxcore_rad(k) ! cubic vortex core radius
        rho = r_mag / wake%multivxcore_rad(k)
    else if (option==1) then ! for wake (acts on points)
        !sym_core_rad = sqrt( ( wake%vxcore_rad(i)*wake%vxcore_rad(i) + wake%vxcore_rad(k)*wake%vxcore_rad(k) ) / 2._pr ) ! symmetrized vortex core radius according to (HE, 2009); could it be applied to the filament/tube's endpoints? Some authors applied it on direct vorton-vorton interaction but not on vorton-tube one. Thus, should it be applied to the half volume per endpoint?
        wake%core_rad2 = wake%vxcore_rad(k)*wake%vxcore_rad(k) ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*wake%vxcore_rad(k) ! cubic vortex core radius
        rho = r_mag / wake%vxcore_rad(k) 
    else
        pause
    end if
        
    r_mag2 = r_mag*r_mag ! squared distance
    r_mag3 = r_mag2*r_mag ! cubic distance
    rho2 = rho*rho; rho3=rho*rho*rho
        
    if (r_mag<=input%tol_rad) then ! vortons are too close
        vind_v2p(:) = 0._pr ! induced velocity is zero
    else if (input%core_rad_init == 0._pr) then ! singularized vortex core; it does not apply
    else ! regularized vortex core
        select case (input%regul_function)
            case (0) ! High-order RF
                g_winck = (rho3 * (rho2 + fiveovertwo)) / ((rho2 + 1._pr)**fiveovertwo) ! High-order RF (WINCKELMANS, 1989)
            case (1) ! 2nd-order Gaussian RF (BERDOWSKY, 2015) 
                g_winck = 1._pr - e_num**(-rho3)
            case (2) ! Gaussian error function (used by several authors)
                cons = 0.5_pr*rho2 !r_mag2/(2._pr*wake%core_rad2)     
                !g_winck = erf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons)
                g_winck = derf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons) ! double precision erf

            end select          
                
        fac1 = -1._pr/fourpi
        fac2 = g_winck / r_mag3 ! (BIRD, 2021)
        fac3 = fac1 * fac2 ! first times second factors
        winck(:) = fac3 * r(:) ! it works for (WINCKELMANS, 2005) and (ÁLVAREZ, 2018); same result
        call CROSSPRODUCT(winck(:), circ_vec(:), vind_v2p(:)) ! induced velocity by a vorton over another one (or over an evaluation point)
    end if

End Subroutine VEL_VORT2POINT

Subroutine DETACHVORT ( step ) ! calculates the vortex strength and orientation of each detached vorton (internal and external wakes); JCPG
    USE ARRAYS, only : pr,kuttaedges,grid,bound_circ,wake,gamma_wake,esed1_reduc,orien_reduc,plate_LE_detected,gamma_wake_start,esed1_modif_Gutnikov,gamma_wake_mod,input
    Implicit none
    integer :: i,j,n,cont,detachBVR,detachBVR_fst,k,m
    integer, intent(in) :: step
    real(kind=pr) :: fst_wake_vort,sec_wake_vort!,sum_circ_plate

    if (step == 0) then
        esed1_modif_Gutnikov(:,:) = 0 ! to save computation
    end if
    
    k=0
    do i = 1, wake%nwp ! wake VR, for 4x4 case: 1->40
        n = kuttaedges%edge2wpan(i) ! wake VR num., for 4x4 case: 17->56
        cont=1
        do j = kuttaedges%esed2(i)+1 , kuttaedges%esed2(i+1) ! 1 (external bound VR) or 2 (internal bound VRs) wakes  
            detachBVR = kuttaedges%esed1(j) ! for 4x4 case: 1 to 16 (bounded VRs)
            if (cont==2) then ! for two wakes case (internal)
                
                sec_wake_vort = bound_circ(detachBVR) ! circulation for second (internal) wake
                if ( (abs(fst_wake_vort)) >= (abs(sec_wake_vort)) ) then
                    gamma_wake(i) = fst_wake_vort - sec_wake_vort ! difference between both bounded VRs
                    esed1_reduc(i) = detachBVR_fst ! writes previous BVR number
                    orien_reduc(i) = kuttaedges%orien(k) 
                    k=k+1
                    
                    if (step == 0) then
                        do m=1, 4
                            if ( esed1_modif_Gutnikov(m,detachBVR)==0 ) then
                                esed1_modif_Gutnikov(m,detachBVR) = i
                                exit
                            else
                            end if
                        end do
                    else
                    end if
                    
                else 
                    gamma_wake(i) = sec_wake_vort - fst_wake_vort ! difference between both bounded VRs
                    esed1_reduc(i) = detachBVR ! writes current BVR number
                    orien_reduc(i) = kuttaedges%orien(k+1)
                    k=k+1
                    
                    if (step == 0) then
                        do m=1, 4
                            if ( esed1_modif_Gutnikov(m,detachBVR)==0 ) then
                                esed1_modif_Gutnikov(m,detachBVR) = i
                                exit
                            else
                            end if
                        end do
                    else
                    end if
                    
                end if
            else ! for first (or only one; external) wake
                fst_wake_vort = bound_circ(detachBVR) ! circulation for first (or only one) wake
                gamma_wake(i) = fst_wake_vort
                detachBVR_fst = detachBVR ! used only for two-wake case
                esed1_reduc(i) = detachBVR ! for two-wake case is rewritted
                k=k+1
                orien_reduc(i) = kuttaedges%orien(k) ! orientation of internal-reduced (and external) wakes
                
                if (step == 0) then ! to avoid extra computation
                    do m=1, 4 ! over four edges' element/panel
                        if ( esed1_modif_Gutnikov(m,detachBVR)==0 ) then
                            esed1_modif_Gutnikov(m,detachBVR) = i ! assigns the wake number to a modified "edges to panel" matrix
                            exit
                        else
                        end if
                    end do
                else
                end if
                
            end if
            cont=cont+1
        end do ! ends j
    end do ! ends i        
    
        gamma_wake_start(:) = gamma_wake(:)
        gamma_wake_mod(:) = gamma_wake(:)
        do i = 1, wake%nwp ! wake VR, for 4x4 case: 1->40
            n = kuttaedges%edge2wpan(i) ! wake VR num., for 4x4 case: 17->56            
            do j=1, grid%plate_nle
                if (n==plate_LE_detected(j)) then
                    
                    select case(input%detach_model)
                    case(1) ! negative LE wake
                        gamma_wake_mod(i) = -gamma_wake(i)
                        gamma_wake_start(i) = -gamma_wake_start(i) ! inverted LE wake circulation condition (only for vortex ring wake, not for vorton method!)
                    case(2) ! this is not necessary (avoids rewriting)
                    end select
                else
                end if
            end do
        end do
        
    !sum_circ_plate = sum(abs(gamma_wake)) ! total plate's circulation    

End Subroutine DETACHVORT

Subroutine FORCES_KJ ( circ_vec, point, tot_area, pos_kvort, den_aux, option ) ! unsteady hydrodynamic coefficient calculation via Kutta-Zhukovski (KJ); JCPG/EOA
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,area,deriv_gamma,bound_circ_old,nor_vec,internal_LE_detected,plate_LE_detected,del_cp,del_pres,ctrl_pts,x_cm,half_chordpan,del_force
    Use MATH, only : CROSSPRODUCT
    Implicit none

    integer :: i,k,fst_nod,sec_nod,m
    integer, intent(in) :: den_aux
    integer, intent(out) :: option
    real(kind=pr) :: sum_mom,den,cfx,cfy,cfz,cl,cd,cy,cn,dist,cm
    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta,force_x,force_y,force_z
    real(kind=pr) :: tot_force_steady(3),pan_force_steady(3),vind_body(3),ra(3),rb(3),vind_cloud(3),vind_total(3),del_gamma(3),segm_force(3),tot_force_unst(3),pan_force_unst(3),tot_force(3),mid_point(3),vind_v2p(3)
    real(kind=pr), intent(in) :: tot_area
    real(kind=pr), intent(out) :: point(3),circ_vec(3),pos_kvort(3)
    logical :: log1,log2

    tot_force_steady(:) = 0._pr
    tot_force_unst(:) = 0._pr
    sum_mom = 0._pr
    x_cm = input%x_pos_mom_fact*input%char_len ! moment axis position

    do i=1, grid%nelem ! over BVR/panels
        deriv_gamma(i) = ( bound_circ(i) - bound_circ_old(i) ) / input%dt ! Gamma time derivative
        pan_force_steady(:) = 0._pr
        pan_force_unst(:) = 0._pr
    
        fst_nod = grid%panel(4,i) ! first node (it depends on mesh numeration's starting point)
        sec_nod = grid%panel(3,i) ! second node (it depends on mesh numeration's starting point)
                
        log1 = .false.; log2 = .false.
        do m=1, grid%nle
            if ( fst_nod==grid%panel(2,internal_LE_detected(m)) .and. sec_nod==grid%panel(1,internal_LE_detected(m))) then
                log1 = .true.
            else
            end if
        end do
        
        do m=1, grid%plate_nle
            if ( fst_nod==grid%panel(2,plate_LE_detected(m)) .and. sec_nod==grid%panel(1,plate_LE_detected(m)) ) then
                select case(input%detach_model)
                case(1) ! negative LE wake
                    log2 = .true.
                case(2) ! positive LE wake (Kutta condition)
                    log2 = .false.
                end select
            else
            end if
        end do 
        
        if (log1==.true. .or. log2==.true.) then
            vind_body(:) = 0._pr
            mid_point(:) = ((grid%coord(:,grid%panel(4,i)) + grid%coord(:,grid%panel(3,i))) / 2._pr) ! vortex segment's midpoint 
        
            point(:) = mid_point(:)  ! position vector of the i-vorton
            
            vind_body(:) = 0._pr ! induced velocity by the body/plate
            option=0 ! for plate case
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(den_aux,point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
            !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
            !$OMP DO
            do k=1, den_aux
                circ_vec(:) = wake%multicirc_vec(:,k)
                pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
            end do 
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL
                    
            vind_cloud(:) = 0._pr ! clears old values
            option=1 ! for wake case
            !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
            !$OMP SHARED(point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
            !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
            !$OMP DO
            do k = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
                circ_vec(:) = wake%gammav(:,k) ! vorton's circulation vector
                pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
                call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
                vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single vorton
            end do ! ends k     
            !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
            !$OMP end PARALLEL
        
            vind_total(:) = vdir(:) + vind_body(:) + vind_cloud(:)  ! total local flow velocity (freestream + plate + wake)
        
            ra(:) = grid%coord(:,grid%panel(4,i)) ! bounded segment's first node position
            rb(:) = grid%coord(:,grid%panel(3,i)) ! bounded segment's second node position
 
            del_gamma(:) = (rb(:)-ra(:)) * bound_circ(i) ! panel's remaining segments
            call CROSSPRODUCT(vind_total(:), del_gamma(:), segm_force(:)) ! force on segment
        else
            segm_force(:) = 0._pr ! for Kutta type separation edge
        end if
        
        segm_force(:) = input%dens*segm_force(:)
        pan_force_steady(:) = pan_force_steady(:) + segm_force(:) ! perpendicular force over the i-panel
        pan_force_unst(:) = input%dens*deriv_gamma(i)*area(i)*nor_vec(:,i) ! perpendicular unsteady force over the i-panel
        tot_force_steady(:) = tot_force_steady(:) + pan_force_steady(:) ! body's steady force contribution
        tot_force_unst(:) = tot_force_unst(:) + pan_force_unst(:) ! body's unsteady force contribution
        
        dist = ctrl_pts(1,i) - x_cm - half_chordpan
        sum_mom = sum_mom + (  (pan_force_steady(3) + pan_force_unst(3)) * dist )
        del_pres(i) = -( pan_force_steady(3) + pan_force_unst(3) ) / area(i) ! -1 to obtain the correct sign according to aerodynamical convention
        del_cp(i) = 2._pr*del_pres(i) / (input%dens * input%q_inf*input%q_inf) ! pressure coefficient jump across the panel            

    end do ! ends i
  
!!!! Pressure integration    
!    alpharad = input%alpha*(pi/180._pr)!; betarad  = input%beta*(pi/180._pr)  ! angles in radians
!    sinalpha = sin(alpharad); cosalpha = cos(alpharad)!; sinbeta = sin(betarad); cosbeta = cos(betarad)
!    !
!    !cl  = cfz*cosalpha - cfx*sinalpha                                ! lift coef.            (CL)
!    !cd = cfx*cosalpha*cosbeta + cfy*sinbeta + cfz*sinalpha*cosbeta   ! drag coef.            (CD)
!    !cy  = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef.   (CY)
!    !cn  = cl*cosalpha + cd*sinalpha                                  ! normal force coef.    (CN)
!    !!cm = sum_mom/(0.5_pr*rho*q_inf*q_inf*(env*char_len)*char_len)   ! pitching moment coef. (CM)
!
!    force_x=0._pr; force_y=0._pr; force_z=0._pr
!    !sum_mom=0._pr
!    do k=1, grid%nelem
!        del_force(:,k) = -(del_pres(k)*area(k))*nor_vec(:,k) ! force per panel
!        !force_x = force_x - del_force(1,k) ! longitudinal force !!! original (-)
!        !!force_y = force_y - del_force(2,k) ! lateral force
!        !force_z = force_z - del_force(3,k) ! vertical force
!        force_x = force_x + del_force(1,k) ! longitudinal force
!        !force_y = force_y + del_force(2,k) ! lateral force
!        force_z = force_z + del_force(3,k) ! vertical force     
!        
!        !dist_seg = ctrl_pts(1,k) !- ((0.5_pr*char_len)/real(trail_lines))
!        dist = ctrl_pts(1,k) - x_cm ! -0.125_pr
!        sum_mom = sum_mom + (  del_force(3,k) * dist ) !((ctrl_pts(1,k)-(char_len/real(trail_lines))) - x_cm)  )
!    end do
!
!    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
!    cfx = 2._pr*force_x / den
!    !cfy = 2._pr*force_y / den
!    cfz = 2._pr*force_z / den
!    
!    cfx = (1._pr/input%char_len)*cfx  ! longitudinal force coefficent
!    !cfy = (1._pr/input%char_len)*cfy ! lateral force coefficent
!    cfz = (1._pr/input%char_len)*cfz  ! vertical force coefficient
!    
!    cl = cfz*cosalpha - cfx*sinalpha ! lift coefficient, CL
!    cd = cfz*sinalpha + cfx*cosalpha ! drag coefficient, CD
!    cm = -sum_mom/(0.5_pr*input%dens*input%q_inf*input%q_inf*tot_area*input%char_len) ! pitching moment coefficient (is multiplied by -1 to represent its value in the aerodynamical sign convention; positive: clockwise) 
!!!!    
   
!!! Direct Kutta-Zhukovski (KJ)
    !TOT_FORCE_UNST(:) = 0._pr ! for testing
    tot_force(:) = tot_force_steady(:) + tot_force_unst(:) ! body's total force

    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
    cfx = 2._pr*tot_force(1)/den; cfy = 2._pr*tot_force(2)/den; cfz = 2._pr*tot_force(3)/den  ! body axes coefficients
    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !input%beta*(pi/180._pr)  ! angles in radians

    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)

    cl  = cfz*cosalpha - cfx*sinalpha                                ! lift coef. (CL)
    cd = cfx*cosalpha*cosbeta + cfy*sinbeta + cfz*sinalpha*cosbeta   ! drag coef. (CD)
    cm = -sum_mom/(0.5_pr*input%dens*input%q_inf*input%q_inf*tot_area*input%char_len) ! pitching moment coefficient (is multiplied by -1 to represent its value in the aerodynamical sign convention; positive: clockwise)
    cy  = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef. (CY)
    cn  = cl*cosalpha + cd*sinalpha                                  ! normal force coef. (CN)
!!!

    400 format(f12.6)
    print 400, cl; print 400, cd; print 400, cm; print 400, cn
    print *,'---'

    500 format(f12.6,f12.6,f12.6,f12.6,f12.6,f12.6)
    open(5, file='aero_coef.dat', access='append') ! writes output file
    write(5,500) cl, cd, cm, cn
    close(5)

End Subroutine FORCES_KJ

Subroutine FORCES_PRES ( circ_vec, point, tot_area, pos_kvort, den_aux, option ) ! Dynnikova's method; JCPG
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,area,deriv_gamma,bound_circ_old,nor_vec,del_cp,del_pres,ctrl_pts,del_force,third_term,inv_4pi,solid_angle,fourth_term,vel_vorton_avg,second_term,fourpi,x_cm,del_mom,del_pres_down,del_force_down,second_term_down,fourth_term_down
    Use MATH, only : CROSSPRODUCT,DOTPRODUCT,VECTORNORM
    Implicit none

    integer :: i,k,j
    integer, intent(in) :: den_aux
    integer, intent(out) :: option
    real(kind=pr), intent(in) :: tot_area
    real(kind=pr) :: den,cfx,cfy,cfz,cl,cd,cy,cn,vel_mag,vel_mag_cp,cm
    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta,force_x,force_y,force_z,vel_total,stat_pres
real(kind=pr) :: force_x_down,force_y_down,force_z_down,force_x_tot,force_y_tot,force_z_tot
    real(kind=pr) :: vind_body(3),vind_cloud(3),vind_total(3),vind_v2p(3),vind_total_cp(3),r_0(3),cross(3),sum_mom_vec(3)
    real(kind=pr), intent(out) :: circ_vec(3),pos_kvort(3),point(3)
   
    do i=1, grid%nelem
            !deriv_gamma(i) = ( bound_circ(i) - bound_circ_old(i) ) / input%dt ! Gamma time derivative (original)
        deriv_gamma(i) = ( bound_circ_old(i) - bound_circ(i) ) / input%dt ! inverted sign Gamma time derivative (Shcheglov, 2008) 
        !deriv_gamma(i) = deriv_gamma(i)/2._pr ! ; uncomment this line for double layer calculation
    end do
    
    ! THIRD TERM
    do i=1, grid%nelem ! on control points
        third_term(i) = 0._pr
        do j=1, grid%nelem ! on panels 
            third_term(i) = third_term(i) + solid_angle(i,j)*deriv_gamma(j)
        end do ! ends panels
    end do ! ends control points
    third_term(:) = inv_4pi*third_term(:) ! third term of eq. 3.9 (Dergachev, 2018)
    !third_term(:) = 2._pr*third_term(:) ! uncomment this line for double layer calculation
    
    ! SECOND AND FOURTH TERMS
    do i=1, grid%nelem ! on control points
        ! SINGLE SURFACE (VLM-BASED) PRESSURE CALCULATION
        second_term(i) = 0._pr
        point(:) = ctrl_pts(:,i) + (/0._pr,0._pr,input%thick/2._pr/) ! control point's position
        
        !!!!--- commented intentionally for testing purposes
        !
        !vind_body(:) = 0._pr ! induced velocity by the body
        !option = 0 ! for body case
        !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !!$OMP SHARED(den_aux,point,option,wake,i) PRIVATE(k,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
        !!$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
        !!$OMP DO
        !do k=1, den_aux
        !    circ_vec(:) = wake%multicirc_vec(:,k)
        !    pos_kvort(:) = wake%pos_multifix2plate(:,k) ! k-vorton's position
        !    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
        !    vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
        !end do 
        !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !!$OMP end PARALLEL
        !
        !!!!---
                    
        vind_cloud(:) = 0._pr ! clears old values
        fourth_term(i) = 0._pr
        option = 1 ! for wake case
        !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !$OMP SHARED(point,option,wake,i,vel_vorton_avg,vdir) PRIVATE(k,circ_vec,pos_kvort,vind_v2p,vel_total) & ! (declaration of variables)
        !$OMP REDUCTION(+: vind_cloud,fourth_term) ! (declaration of reduction variables)
        !$OMP DO
        do k = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
            circ_vec(:) = wake%gammav(:,k) ! vorton's circulation vector
            pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
            call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
            vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single vorton
            call DOTPRODUCT(vind_v2p(:),vel_vorton_avg(:,k - wake%nwp) + vdir(:),vel_total) ! +vdir(:) (see VM2D code, Marchevsky et al.: Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + V0; (file: MeasureVP2D.cpp, line 288))           
            fourth_term(i) = fourth_term(i) + vel_total
        end do ! ends k     
        !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !$OMP end PARALLEL
        
        vind_total(:) = vdir(:) + vind_cloud(:)  ! local flow velocity (freestream + wake)
            !vind_total(:) = vdir(:) + vind_cloud(:) + vind_body(:)  ! local flow velocity (freestream + wake + body) ! for testing
            !vind_total(:) = vind_cloud(:)  ! local flow velocity (freestream + wake) ! for testing
        
        call VECTORNORM(vind_total(:),vel_mag)
        second_term(i) = vel_mag

        del_pres(i) = input%pres_inf + input%dens * ( ((input%q_inf * input%q_inf)/2._pr) - ((second_term(i)*second_term(i))/2._pr) - third_term(i) + fourth_term(i) )

        !!!!--- DOWNSIDE ELEMENT'S PRESSURE CALCULATION; uncomment this section for double layer calculation
        !second_term_down(i) = 0._pr
        !point(:) = ctrl_pts(:,i) - (/0._pr,0._pr,input%thick/2._pr/)  ! control point's position
        !                
        !vind_cloud(:) = 0._pr ! clears old values
        !fourth_term_down(i) = 0._pr
        !option = 1 ! for wake case
        !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !!$OMP SHARED(point,option,wake,i,vel_vorton_avg,vdir) PRIVATE(k,circ_vec,pos_kvort,vind_v2p,vel_total) & ! (declaration of variables)
        !!$OMP REDUCTION(+: vind_cloud,fourth_term_down) ! (declaration of reduction variables)
        !!$OMP DO
        !do k = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
        !    circ_vec(:) = wake%gammav(:,k) ! vorton's circulation vector
        !    pos_kvort(:) = wake%pos(:,k) ! k-vorton's position
        !    call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,k,option) ! calculates the induced velocity by a vorton over another one
        !    vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single vorton
        !    call DOTPRODUCT(vind_v2p(:),vel_vorton_avg(:,k - wake%nwp) + vdir(:),vel_total) ! +vdir(:) (see VM2D code, Marchevsky et al.: Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + V0; (file: MeasureVP2D.cpp, line 288))
        !    fourth_term_down(i) = fourth_term_down(i) + vel_total
        !end do ! ends k     
        !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !!$OMP end PARALLEL
        !
        !vind_total(:) = vdir(:) + vind_cloud(:)  ! local flow velocity (freestream + wake)
        !call VECTORNORM(vind_total(:),vel_mag)
        !second_term_down(i) = vel_mag
        !
        !del_pres_down(i) = input%pres_inf + input%dens * ( ((input%q_inf * input%q_inf)/2._pr) - ((second_term_down(i)*second_term_down(i))/2._pr) - third_term(i) + fourth_term_down(i) )
        !!!!---

        
        !!!--- for testing (pending)
        !vind_total_cp(:) = vdir(:) + vind_cloud(:) + vind_body(:) ! for cp calculation
        !call VECTORNORM(vind_total_cp(:),vel_mag_cp)
        !del_cp(i) = 2._pr*del_pres(i) / (input%dens * input%q_inf*input%q_inf) ! pressure coefficient
        !     !del_cp(i) = 1._pr - ( (vel_mag_cp*vel_mag_cp) / (input%q_inf*input%q_inf) )
        !     !stat_pres = del_pres(i) - (0.5_pr*input%dens * input%q_inf*input%q_inf)
        !     !del_cp(i) = 2._pr*(stat_pres - input%pres_inf) / (input%dens * input%q_inf*input%q_inf)
        !     !del_cp(i) = 2._pr*stat_pres / (input%dens * input%q_inf*input%q_inf)
        !!!---
    end do ! ends i    

    del_pres(:) = -2._pr*del_pres(:) ! for single layer calculation    
    !del_pres(:) = -del_pres(:); del_pres_down(:) = -del_pres_down(:) ! uncomment this line for double layer calculation      
        
    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !betarad  = input%beta*(pi/180._pr) ! angles in radians
    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)

    force_x=0._pr; force_y=0._pr; force_z=0._pr ! for single layer calculation  
    !force_x_down=0._pr; force_y_down=0._pr; force_z_down=0._pr ! uncomment this line for double layer calculation

    x_cm = input%x_pos_mom_fact * input%char_len ! moment axis position
    r_0(:) = (/x_cm, 0.5_pr, 0._pr/) ! moment's position vector
    sum_mom_vec(:) = 0._pr

    do k=1, grid%nelem
        del_force(:,k) = -del_pres(k) * area(k) * nor_vec(:,k) ! force per panel
        !del_force_down(:,k) = -del_pres_down(k) * area(k) * -nor_vec(:,k) ! uncomment this line for double layer calculation

        force_x = force_x + del_force(1,k); force_y = force_y + del_force(2,k); force_z = force_z + del_force(3,k) ! force's components 
        !force_x_down = force_x_down + del_force_down(1,k); force_y_down = force_y_down + del_force_down(2,k); force_z_down = force_z_down + del_force_down(3,k) ! uncomment this line for double layer calculation
        
        call CROSSPRODUCT(ctrl_pts(:,k) - r_0(:), nor_vec(:,k), cross(:)) 
        del_mom(:,k) = del_pres(k) * area(k) * cross(:)
        sum_mom_vec(:) = sum_mom_vec(:) + del_mom(:,k)
    end do

    !force_x_tot = force_x - force_x_down; force_y_tot = force_y - force_y_down; force_z_tot = force_z - force_z_down ! uncomment this line for double layer calculation

    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
    
    cfx = 2._pr*force_x / den; cfy = 2._pr*force_y / den; cfz = 2._pr*force_z / den ! for single layer calculation  
    !cfx = 2._pr*force_x_tot / den; cfy = 2._pr*force_y_tot / den; cfz = 2._pr*force_z_tot / den ! uncomment this line for double layer calculation

    cfx = (1._pr/input%char_len)*cfx; cfy = (1._pr/input%char_len)*cfy; cfz = (1._pr/input%char_len)*cfz ! longitudinal, lateral and vertical force coefficent
    
    cl = cfz*cosalpha - cfx*sinalpha ! lift coefficient, CL
    cd = cfz*sinalpha + cfx*cosalpha ! drag coefficient, CD
    cm = -sum_mom_vec(2) / (0.5_pr*input%dens*input%q_inf*input%q_inf*tot_area*input%char_len) ! pitching moment coefficient (is multiplied by -1 to represent its value in the aerodynamical sign convention; positive: clockwise)
    cy = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef. (CY)
    cn = cl*cosalpha + cd*sinalpha ! normal force coefficient, CN
    
    400 format(f12.6)
    print 400, cl; print 400, cd; print 400, cm; print 400, cn
    print *,'---'

    500 format(f12.6,f12.6,f12.6,f12.6)
    open(5, file='aero_coef.dat', access='append') ! writes output file
    write(5,500) cl, cd, cm, cn
    close(5)
    
End Subroutine FORCES_PRES
    
Subroutine VORTONRING ( n,circ,point,vind_ring,op,equiv_vcr ) ! calculates the induced velocity by a multi-vorton ring; JCPG
    Use ARRAYS, only : pr,grid,input
    Use MATH, only : VECTORNORM
    Implicit none
    real(kind=pr), intent(in) :: circ
    real(kind=pr), intent(out) :: vind_ring(3),equiv_vcr
    real(kind=pr) :: ra(3),rb(3),ds(3),pos_kvort(3),circ_vec(3),vind_v2p(3),ds_len,del_pos(3),den,vind_segm(3)!,equiv_vcr
    real(kind=pr), intent(inout) :: point(3)
    integer, intent(in) :: n,op
    integer :: nodes(5),k,vps
    
    vind_ring(:) = 0._pr
    nodes(:) = [ grid%panel(1:grid%elemtype,n) , grid%panel(1,n) ] ! panel's number

    do k=1, grid%elemtype ! over the 3 or 4 vortex segments
        !if (k<grid%elemtype) then ! first 2 or 3 segments (clockwise direction; KATZ, 2001)
        !    ra = grid%coord(:,grid%panel(k+1,n))  
        !    rb = grid%coord(:,grid%panel(k,n))   
        !else ! last segment
        !    ra = grid%coord(:,grid%panel(1,n))   
        !    rb = grid%coord(:,grid%panel(grid%elemtype,n))   
        !end if
        if (op==1) then ! vorton wake's case
            ra(:) = grid%coord_start(:,nodes(k+1)) ! clockwise direction
            rb(:) = grid%coord_start(:,nodes(k)) 
        else ! bounded vorton rings' case
            ra(:) = grid%coord(:,nodes(k+1)) ! clockwise direction  
            rb(:) = grid%coord(:,nodes(k)) 
        end if
        
        ds(:) = rb(:) - ra(:) ! vortex segment's vector
        call VECTORNORM(ds(:),ds_len) ! segment's lenght
        
        den = real(ceiling(ds_len/input%core_rad_init)) + 1._pr ! number of vortons per segment
             !den=1._pr ! for testing
        equiv_vcr = input%core_rad_init / ( den**(1._pr/3._pr) )

        del_pos(:) = ds(:) / den  ! delta position
        circ_vec(:) = circ*del_pos(:)
        
        vind_segm(:) = 0._pr
        do vps=1, int(den) ! vortons per segment
            pos_kvort(:) = ra(:) + del_pos(:)/2._pr + (vps-1)*del_pos(:) ! extra k-vorton position
            call VEL_VORT2POINT_START(circ_vec,point,pos_kvort,vind_v2p,equiv_vcr) ! calculates the induced velocity by a vorton over a point 
            vind_segm(:) = vind_segm(:) + vind_v2p(:)
        end do

        vind_ring(:) = vind_ring(:) + vind_segm(:) ! induced velocity by a bounded vorton ring over a point 
    end do
  
End Subroutine VORTONRING
    
Subroutine VEL_VORT2POINT_START (circ_vec,point,pos_kvort,vind_v2p,equiv_vcr) ! calculates the induced velocity for a vorton over another one (or control point); JCPG
    Use ARRAYS, only : pr,input,wake,fiveovertwo,fourpi,fouroverpi,e_num,twopi,pi,twooverpi
    Use MATH, only : VECTORNORM,CROSSPRODUCT,NORMALIZE
    implicit none
    real(kind=pr) :: r_mag,r_mag2,r_mag3,fac1,fac2,fac3,cons,g_winck,rho,rho2,rho3
    real(kind=pr) :: r(3),winck(3)
    real(kind=pr), intent(in) :: circ_vec(3),point(3),pos_kvort(3),equiv_vcr
    real(kind=pr), intent(out) :: vind_v2p(3)

    wake%core_rad2 = equiv_vcr*equiv_vcr
    wake%core_rad3 = wake%core_rad2*equiv_vcr
            
    r(:) = point(:) - pos_kvort(:) ! radio vector between k-vorton and an evaluation point
    CALL VECTORNORM(r(:), r_mag) ! vector norm (distance between points)
    
    r_mag2 = r_mag*r_mag ! squared distance
    r_mag3 = r_mag2*r_mag ! cubic distance
    rho = r_mag / equiv_vcr; rho2=rho*rho; rho3=rho*rho*rho
        
    if (r_mag<=input%tol_rad) then ! particles/vortons are too close
        vind_v2p(:) = 0._pr ! induced velocity is zero
    else if (input%core_rad_init == 0._pr) then ! singularized vortex core radius; it does not apply
    else ! regularized vortex core
      
        select case (input%regul_function)
        case (0) ! High-order RF
            g_winck = (rho3 * (rho2 + fiveovertwo)) / ((rho2 + 1._pr)**fiveovertwo) ! High-order RF (WINCKELMANS, 1989)
        case (1) ! 2nd-order Gaussian RF (BERDOWSKY, 2015) 
            g_winck = 1._pr - e_num**(-rho3)
        case (2)
            cons = 0.5_pr*rho2 !r_mag2/(2._pr*wake%core_rad2)     
            !g_winck = erf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons)
            g_winck = derf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons) ! double precision erf

        end select
            
        fac1 = -1._pr/fourpi
        fac2 = g_winck / r_mag3 ! (BIRD, 2021)
        fac3 = fac1 * fac2 ! first times second factors
        winck(:) = fac3 * r(:)
        call CROSSPRODUCT(winck(:), circ_vec(:), vind_v2p(:)) ! induced velocity by a vorton over another one (or over an evaluation point)
    end if

End Subroutine VEL_VORT2POINT_START

Subroutine POINT_POLYGON (step) ! Allows to relocate the vortex elements that cross/penetrate the plate; JCPG
    Use ARRAYS, only : pr,input,wake,vel_fst_point_old,vel_sec_point_old,fst_plate_vertex_pos, sec_plate_vertex_pos,nod2fil_extend,kuttaedges,nod2fil_reduc,grid,spanwise_wake_elem,chordwise_wake_elem
    Use MATH, only: DOTPRODUCT,VECTORNORM
    Use UTILS, only : RESIZE_INT_2, RESIZE_LOG_1
    Implicit none
    integer :: i,k,cont,j,node_1,node_2,m,aa
    integer, intent(in) :: step
    real(kind=pr) :: or_1,or_2,or_3,or_4,p1_x,p1_y,p2_x,p2_y,q1_x,q1_y,q2_x,q2_y
    logical :: aux_or(4)
    logical(1), allocatable :: correc_pos(:)

    if (step==1) then
        allocate ( nod2fil_reduc(2,kuttaedges%nte), nod2fil_extend(2,kuttaedges%nte),correc_pos(wake%nwp) )
        nod2fil_reduc(:,:) = 0

        do m=1, grid%nle + grid%plate_nle + grid%plate_nte 
            aa=spanwise_wake_elem(m) - grid%nelem ! wake element (1,3,5,7,10,12 for 2x2 mesh)
            nod2fil_reduc(1,aa) = grid%panel(1,spanwise_wake_elem(m))
            nod2fil_reduc(2,aa) = grid%panel(2,spanwise_wake_elem(m))
        end do ! ends m
            
        do m=1, wake%nwp - (grid%nle + grid%plate_nle + grid%plate_nte)
            aa=chordwise_wake_elem(m) - grid%nelem ! wake element (2,4,6,8,9,11 for 2x2 mesh)
                nod2fil_reduc(1,aa) = grid%panel(1,chordwise_wake_elem(m))
                nod2fil_reduc(2,aa) = grid%panel(2,chordwise_wake_elem(m))
        end do ! ends m
    else
    end if

    if (step==1) nod2fil_extend(:,1:wake%nwp) = nod2fil_reduc(:,1:wake%nwp)

    call RESIZE_INT_2 ( 2, wake%nwp*step, nod2fil_extend) ! all nodes id (body + multiple wakes)
    call RESIZE_LOG_1 (wake%nwp*step, correc_pos)

    if(step>1) then
        do i=wake%nwp*step, wake%nwp + 1, -1 ! wake vorton's number ! 80 -> 41
            nod2fil_extend(:,i) = nod2fil_extend(:,i-wake%nwp ) + grid%nnod
        end do
    else
    end if

    aux_or(:) = .false.
    correc_pos(:) = .false.
    if (step==1) then
        if (input%mesh_num == 1 .or. input%mesh_num == 2 .or. input%mesh_num == 3 .or. input%mesh_num == 4) then ! for quadrangular (AR=1 case)
            fst_plate_vertex_pos(1,1) = 0.999_pr  ; fst_plate_vertex_pos(2,1) = 0.001_pr ! with a small increment (0.001) to avoid numerical issues (quadrangular shape is converted into an octagonal one)
            fst_plate_vertex_pos(1,2) = 0.5_pr ; fst_plate_vertex_pos(2,2) = -0.001_pr
            fst_plate_vertex_pos(1,3) = 0.001_pr  ; fst_plate_vertex_pos(2,3) = 0.001_pr
            fst_plate_vertex_pos(1,4) = -0.001_pr  ; fst_plate_vertex_pos(2,4) = 0.5_pr
            fst_plate_vertex_pos(1,5) = -0.001_pr  ; fst_plate_vertex_pos(2,5) = 0.5_pr
            fst_plate_vertex_pos(1,6) = 0.001_pr  ; fst_plate_vertex_pos(2,6) = 0.999_pr
            fst_plate_vertex_pos(1,7) = 0.5_pr ; fst_plate_vertex_pos(2,7) = 1.001_pr
            fst_plate_vertex_pos(1,8) = 0.999_pr  ; fst_plate_vertex_pos(2,8) = 0.999_pr
    
            sec_plate_vertex_pos(1,1) = 1.001_pr  ; sec_plate_vertex_pos(2,1) = 0.5_pr
            sec_plate_vertex_pos(1,2) = 0.999_pr  ; sec_plate_vertex_pos(2,2) = 0.001_pr
            sec_plate_vertex_pos(1,3) = 0.5_pr ; sec_plate_vertex_pos(2,3) = -0.001_pr
            sec_plate_vertex_pos(1,4) = 0.001_pr  ; sec_plate_vertex_pos(2,4) = 0.001_pr
            sec_plate_vertex_pos(1,5) = 0.001_pr  ; sec_plate_vertex_pos(2,5) = 0.999_pr
            sec_plate_vertex_pos(1,6) = 0.5_pr ; sec_plate_vertex_pos(2,6) = 1.001_pr
            sec_plate_vertex_pos(1,7) = 0.999_pr  ; sec_plate_vertex_pos(2,7) = 0.999_pr
            sec_plate_vertex_pos(1,8) = 1.001_pr  ; sec_plate_vertex_pos(2,8) = 0.5_pr
        else if (input%mesh_num == 5) then ! for AR=0.5 flat plate
            fst_plate_vertex_pos(1,1) = 0.999_pr  ; fst_plate_vertex_pos(2,1) = 0.001_pr
            fst_plate_vertex_pos(1,2) = 0.5_pr ; fst_plate_vertex_pos(2,2) = -0.001_pr
            fst_plate_vertex_pos(1,3) = 0.001_pr  ; fst_plate_vertex_pos(2,3) = 0.001_pr
            fst_plate_vertex_pos(1,4) = -0.001_pr  ; fst_plate_vertex_pos(2,4) = 0.25_pr
            fst_plate_vertex_pos(1,5) = -0.001_pr  ; fst_plate_vertex_pos(2,5) = 0.25_pr
            fst_plate_vertex_pos(1,6) = 0.001_pr  ; fst_plate_vertex_pos(2,6) = 0.499_pr
            fst_plate_vertex_pos(1,7) = 0.5_pr ; fst_plate_vertex_pos(2,7) = 0.501_pr
            fst_plate_vertex_pos(1,8) = 0.999_pr  ; fst_plate_vertex_pos(2,8) = 0.499_pr
    
            sec_plate_vertex_pos(1,1) = 1.001_pr  ; sec_plate_vertex_pos(2,1) = 0.25_pr
            sec_plate_vertex_pos(1,2) = 0.999_pr  ; sec_plate_vertex_pos(2,2) = 0.001_pr
            sec_plate_vertex_pos(1,3) = 0.5_pr ; sec_plate_vertex_pos(2,3) = -0.001_pr
            sec_plate_vertex_pos(1,4) = 0.001_pr  ; sec_plate_vertex_pos(2,4) = 0.001_pr
            sec_plate_vertex_pos(1,5) = 0.001_pr  ; sec_plate_vertex_pos(2,5) = 0.499_pr
            sec_plate_vertex_pos(1,6) = 0.5_pr ; sec_plate_vertex_pos(2,6) = 0.501_pr
            sec_plate_vertex_pos(1,7) = 0.999_pr  ; sec_plate_vertex_pos(2,7) = 0.499_pr
            sec_plate_vertex_pos(1,8) = 1.001_pr  ; sec_plate_vertex_pos(2,8) = 0.25_pr
        else if (input%mesh_num == 6 .or. input%mesh_num == 7) then ! for 45 deg. swept-back (AR=1) flat plate
            fst_plate_vertex_pos(1,1) = 1.499_pr ; fst_plate_vertex_pos(2,1) = 0.001_pr
            fst_plate_vertex_pos(1,2) = 1._pr  ; fst_plate_vertex_pos(2,2) = -0.001_pr
            fst_plate_vertex_pos(1,3) = 0.5_pr ; fst_plate_vertex_pos(2,3) = 0.001_pr
            fst_plate_vertex_pos(1,4) = 0._pr  ; fst_plate_vertex_pos(2,4) = 0.5_pr
            fst_plate_vertex_pos(1,5) = 0._pr  ; fst_plate_vertex_pos(2,5) = 0.5_pr
            fst_plate_vertex_pos(1,6) = 0.5_pr ; fst_plate_vertex_pos(2,6) = 0.999_pr
            fst_plate_vertex_pos(1,7) = 1._pr  ; fst_plate_vertex_pos(2,7) = 1.001_pr
            fst_plate_vertex_pos(1,8) = 1.499_pr ; fst_plate_vertex_pos(2,8) = 0.999_pr
    
            sec_plate_vertex_pos(1,1) = 1._pr  ; sec_plate_vertex_pos(2,1) = 0.5_pr
            sec_plate_vertex_pos(1,2) = 1.499_pr ; sec_plate_vertex_pos(2,2) = 0.001_pr
            sec_plate_vertex_pos(1,3) = 1._pr  ; sec_plate_vertex_pos(2,3) = -0.001_pr
            sec_plate_vertex_pos(1,4) = 0.5_pr ; sec_plate_vertex_pos(2,4) = 0.001_pr
            sec_plate_vertex_pos(1,5) = 0.5_pr ; sec_plate_vertex_pos(2,5) = 0.999_pr
            sec_plate_vertex_pos(1,6) = 1._pr  ; sec_plate_vertex_pos(2,6) = 1.001_pr
            sec_plate_vertex_pos(1,7) = 1.499_pr ; sec_plate_vertex_pos(2,7) = 0.999_pr
            sec_plate_vertex_pos(1,8) = 1._pr  ; sec_plate_vertex_pos(2,8) = 0.5_pr
        else
        end if
    else ! for next steps
    end if

    !!! RELOCATES ONLY THE VORTONS THAT CROSS/PENETRATE THE PLATE (scheme 2: non-div-free grid)
    !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !$OMP SHARED(input,wake,fst_plate_vertex_pos,sec_plate_vertex_pos,correc_pos) PRIVATE(i,k,p1_x,p1_y,q1_x,q1_y,cont,p2_x,p2_y,q2_x,q2_y,or_1,or_2,or_3,or_4,aux_or) !& ! (declaration of variables)
    !!$OMP SHARED(input,wake,fst_plate_vertex_pos,sec_plate_vertex_pos,correc_pos,vel_fst_point_old,vel_sec_point_old) PRIVATE(i,k,p1_x,p1_y,q1_x,q1_y,cont,p2_x,p2_y,q2_x,q2_y,or_1,or_2,or_3,or_4,aux_or) !& ! (declaration of variables)
    !$OMP DO
    do i=1, wake%num_vort ! over the number of current wake vortons
        !if ( wake%pos(3,i) <= wake%vxcore_rad(i)*0.99_pr ) then ! if vertical position is below its vortex core radius value (vorton penetration)
        if ( wake%pos(3,i) <= input%eps*0.99_pr ) then ! if vertical position is below 'almost' the nascent distance epsilon
            p1_x = wake%pos(1,i); p1_y = wake%pos(2,i) ! vorton's initial point
            q1_x = 100._pr;       q1_y = wake%pos(2,i) ! evaluation final point (far downstream)
            cont = 0
            do k=1, 8 ! for all edges
                if (k==1 .or. k==8) then ! for rightside edges
                    p2_x = fst_plate_vertex_pos(1,k) + wake%vxcore_rad(i); p2_y = fst_plate_vertex_pos(2,k) ! edge segment's starting point
                    q2_x = sec_plate_vertex_pos(1,k) + wake%vxcore_rad(i); q2_y = sec_plate_vertex_pos(2,k) ! edge segment's final point
                else if (k==2 .or. k==3) then ! for lowside edges
                    p2_x = fst_plate_vertex_pos(1,k);                      p2_y = fst_plate_vertex_pos(2,k) - wake%vxcore_rad(i) ! edge segment's starting point
                    q2_x = sec_plate_vertex_pos(1,k);                      q2_y = sec_plate_vertex_pos(2,k) - wake%vxcore_rad(i) ! edge segment's final point
                else if (k==4 .or. k==5) then ! for leftside edges
                    p2_x = fst_plate_vertex_pos(1,k) - wake%vxcore_rad(i); p2_y = fst_plate_vertex_pos(2,k) ! edge segment's starting point
                    q2_x = sec_plate_vertex_pos(1,k) - wake%vxcore_rad(i); q2_y = sec_plate_vertex_pos(2,k) ! edge segment's final point
                else if (k==6 .or. k==7) then ! for upside edges
                    p2_x = fst_plate_vertex_pos(1,k);                      p2_y = fst_plate_vertex_pos(2,k) + wake%vxcore_rad(i) ! edge segment's starting point
                    q2_x = sec_plate_vertex_pos(1,k);                      q2_y = sec_plate_vertex_pos(2,k) + wake%vxcore_rad(i) ! edge segment's final point
                end if

                or_1 = -(p2_y-q1_y)*(q1_x-p1_x) ! orientation with respect to the starting point of the panel's segment
                or_2 = -(q2_y-q1_y)*(q1_x-p1_x) ! orientation with respect to the final point of the panel's segment
                or_3 =  (q2_y-p2_y)*(p1_x-q2_x) - (p1_y-q2_y)*(q2_x-p2_x) ! orientation with respect to the starting point of the vortex segment
                or_4 =  (q2_y-p2_y)*(q1_x-q2_x) - (q1_y-q2_y)*(q2_x-p2_x) ! orientation with respect to the final point of the vortex segment
                
                if (or_1 > 0._pr) then  ! panel's start point (clockwise direction)
                    aux_or(1) = .true.
                else
                    aux_or(1) = .false. ! panel's start point (counter-clockwise direction)
                end if
                
                if (or_2 > 0._pr) then  ! panel's final point (clockwise direction)
                    aux_or(2) = .true.
                else
                    aux_or(2) = .false. ! panel's final point (counter-clockwise direction)
                end if
                
                if (or_3 > 0._pr) then  ! vortex segment's start point (clockwise direction)
                    aux_or(3) = .true.
                else
                    aux_or(3) = .false. ! vortex segment's start point (counter-clockwise direction)
                end if
                
                if (or_4 > 0._pr) then  ! vortex segment's final point (clockwise direction)
                    aux_or(4) = .true.
                else
                    aux_or(4) = .false. ! vortex segment's final point (counter-clockwise direction)
                end if
                
                if ((aux_or(1) /= aux_or(2)) .and. (aux_or(3) /= aux_or(4))) then ! the point crosses this edge
                    cont = cont + 1
                else ! the point does not cross this edge
                end if
            end do
        
            if (mod(cont,2) /= 0._pr) then ! odd, the point is over the plate
                wake%pos_fst_point(3,i) = input%eps ! first point's new vertical position
                wake%pos_sec_point(3,i) = input%eps ! second point's new vertical position
                    !wake%pos_fst_point(3,i) = wake%vxcore_rad(i) ! first point's new vertical position
                    !wake%pos_sec_point(3,i) = wake%vxcore_rad(i) ! second point's new vertical position
                !vel_fst_point_old(3,i) = 0._pr ! first point's new vertical velocity due to the wall boundary; can be physically justified to put the vertical velocity to zero after its relocation?
                !vel_sec_point_old(3,i) = 0._pr ! second point's new vertical velocity due to the wall boundary
                correc_pos(i) = .true.
            else ! even (=0)  the point is outside the plate
            end if
        
        else ! even (=0)  the point is far the plate
            
!!!---improvement v1.1 (only for quadrangular AR=1 case); modify for rectangular AR=0.5 and swept-back plate cases!   
            if ( wake%pos_fst_point(3,i) <= input%eps*0.99_pr ) then ! z-direction
                if ( wake%pos_fst_point(2,i) >= 0._pr .and. wake%pos_fst_point(2,i) <= 1._pr ) then ! y-direction
                    if ( wake%pos_fst_point(1,i) >= 0._pr .and. wake%pos_fst_point(1,i) <= 1._pr ) then ! x-direction
                        wake%pos_fst_point(3,i) = input%eps
                    end if
                end if
            else if ( wake%pos_sec_point(3,i) <= input%eps*0.99_pr ) then ! z-direction
                if ( wake%pos_sec_point(2,i) >= 0._pr .and. wake%pos_sec_point(2,i) <= 1._pr ) then ! y-direction
                    if ( wake%pos_sec_point(1,i) >= 0._pr .and. wake%pos_sec_point(1,i) <= 1._pr ) then ! x-direction
                        wake%pos_sec_point(3,i) = input%eps
                    end if
                end if
            end if
!!!---ends improvement
            
        end if
    end do
    !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !$OMP end PARALLEL
    !!!
            
    !!!! UNCOMMENT (AND MAINTAIN ACTIVE THE PREVIOUS SECTION) TO ENFORCE A PERFECT DIV-FREE GRID WHEN VORTONS IMPACT TO THE PLATE (scheme 3: it adjust the surrounding filaments/tubes' nodes)
    !!$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
    !!$OMP SHARED(wake,nod2fil_extend,step,correc_pos) PRIVATE(node_1,node_2,i,j,m) ! (declaration of variables)
    !!$OMP DO
    !do i=1, wake%num_vort ! over the number of current wake vortons
    !    if (correc_pos(i) == .true.) then
    !        node_1 = nod2fil_extend(1,i)
    !        node_2 = nod2fil_extend(2,i)
    !        do j=1, 2 ! over first and second position
    !            do m=1, wake%nwp*step 
    !                if (nod2fil_extend(j,m) == node_1 .and. i/=m) then ! search for first node
    !                    if (j==1) then
    !                        wake%pos_fst_point(:,m) = wake%pos_fst_point(:,i)
    !                        !vel_fst_point_old(3,m) = 0._pr ! first point's new vertical velocity due to the wall boundary
    !                    else ! (j==2)
    !                        wake%pos_sec_point(:,m) = wake%pos_fst_point(:,i)
    !                        !vel_sec_point_old(3,m) = 0._pr ! second point's new vertical velocity due to the wall boundary
    !                    end if
    !                else if (nod2fil_extend(j,m) == node_2 .and. i/=m) then ! search for second node
    !                    if (j==1) then
    !                        wake%pos_fst_point(:,m) = wake%pos_sec_point(:,i)
    !                        !vel_fst_point_old(3,m) = 0._pr ! first point's new vertical velocity due to the wall boundary
    !                    else ! (j==2)
    !                        wake%pos_sec_point(:,m) = wake%pos_sec_point(:,i)
    !                        !vel_sec_point_old(3,m) = 0._pr ! second point's new vertical velocity due to the wall boundary
    !                    end if
    !                else
    !                end if
    !            end do !ends m
    !        end do !ends j
    !    else
    !    end if
    !end do
    !!$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
    !!$OMP end PARALLEL
    !!!!

End Subroutine POINT_POLYGON

Subroutine PANELGEO (geometry,pco,i) ! panel geometry (tangential and normal vectors, control points and areas); EOA
    Use ARRAYS, only : pr,grid,oneoverthree
    Use MATH, only : DOTPRODUCT,VECTORNORM,CROSSPRODUCT,NORMALIZE
    Implicit none
    integer,intent(in) :: i
    real(kind=pr),intent(out) :: geometry(10),pco(3)!,dcarac
   
    ! local variables
    integer :: k
    real(kind=pr) :: p1(3),p2(3),p3(3),p4(3),vmod,t1(3),t2(3),t3(3),d1,d2  
   
    ! ·······················································································      
    !   calculates normal and tangential vectors, control points and panel's area (tria/quad)
    ! ·······················································································      
    
    select case (grid%elemtype) ! pnod
        case (4) ! quadrilateral
   
            do k = 1,3 ! panel coordinates
                p1(k) = grid%coord(k,grid%panel(1,i))
                p2(k) = grid%coord(k,grid%panel(2,i))
                p3(k) = grid%coord(k,grid%panel(3,i))
                p4(k) = grid%coord(k,grid%panel(4,i))
            end do
   
            pco = (p1+p2+p3+p4)*0.25_pr ! control point  
            t1 = (p2+p3)*0.50_pr - pco ! characteristic size (^2)
            t2 = (p3+p4)*0.50_pr - pco
    
            call VECTORNORM(t1,d1) ! vector norm of t1
            call VECTORNORM(t2,d2) ! vector norm of t2
        
            t1 = p3 - p1 ! normal vector (n)
            t2 = p4 - p2
    
            call CROSSPRODUCT(t1,t2,t3)
            call NORMALIZE(t3,vmod)
    
            geometry(10) = vmod*0.50 ! panel area
            t2 = (p3+p4)*0.50_pr - pco ! second tangential vector (m)
    
            call NORMALIZE(t2,vmod)
            call CROSSPRODUCT (t2,t3,t1) ! first tangential vector (l)
    
            geometry(1:3) = t1(:) ! l              
            geometry(4:6) = t2(:) ! m              
            geometry(7:9) = t3(:) ! n

        case (3) ! triangle (it does not apply)

            !p1 = grid%coord(:,grid%panel(1,i))
            !p2 = grid%coord(:,grid%panel(2,i))
            !p3 = grid%coord(:,grid%panel(3,i))
            !
            !pco = (p1+p2+p3)*oneoverthree    
            !
            !t1 = (p1+p2)*0.50_pr - pco 
            !t2 = (p3+p2)*0.50_pr - pco
            !t3 = (p3+p1)*0.50_pr - pco
            !
            !call DOTPRODUCT(t1,t1,d1)
            !call DOTPRODUCT(t2,t2,d2)
            !call DOTPRODUCT(t3,t3,d3)
            !
            !!dcarac = max(d1,d2,d3)*4.0_pr  
            !
            !t1 = p2 - p1 
            !t2 = p3 - p1
            !
            !call CROSSPRODUCT(t1,t2,t3)
            !call NORMALIZE(t3,vmod)
            !
            !geometry(10) = vmod*0.50_pr 
            !
            !t2 = (p3+p1)*0.50_pr - pco 
            !
            !call NORMALIZE(t2,vmod)
            !call CROSSPRODUCT(t2,t3,t1)  
            !
            !geometry(1:3) = t1(:) ! l              
            !geometry(4:6) = t2(:) ! m              
            !geometry(7:9) = t3(:) ! n
    
    end select
   
    !dcarac = dcarac*farfield*farfield ! panel's characteristic distance factor (^2)   
end Subroutine PANELGEO

Subroutine OUTPUT_PLATE (step,den_aux) ! writes the output results file (for the plate's vortons) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
    Implicit none
    integer :: i,j,totnod,totpan
    integer, intent(in) :: step,den_aux
    character(12) :: it,it2,it3,st
          
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(3,file='C:\Users\pimen\Documents\Visual Studio 2017\VORTONEX_v1\VortoNeX\output_vtk\platevort_'//trim(st)//'.vtk')
        write(3,'(a)')'# vtk DataFile Version 2.0'
        write(3,'(a)') 'This file was created by Fortran'
        write(3,'(a)') 'ASCII'
        write(3,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(it, '(i12)') grid%nnod + den_aux
        it=adjustl(it)        
        write(3,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (4)
            do i = 1,grid%nnod
                write(3,'(3(1x,e14.6))') grid%coord(1:3,i)
            end do
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(3(1x,e14.6))') wake%pos_multifix2plate(1,i), wake%pos_multifix2plate(2,i), wake%pos_multifix2plate(3,i)
            end do
            write(3,'(a)') ''
            write(it2, '(i12)') grid%nelem + den_aux
            write(it3, '(i12)') (grid%elemtype+1)*grid%nelem + 2*den_aux
            it2=adjustl(it2); it3=adjustl(it3)        
            write(3,'(a)') 'CELLS  '//trim(it2)//'  '//trim(it3)//''
            do i = 1,grid%nelem
                write(3,'(5(i7,i7))') grid%elemtype, (grid%panel(j,i)-1,j=1,4)  ! -1 to take into account 0 position
            end do
            do i = grid%nnod, den_aux + grid%nnod-1 ! -1 to take into account 0 position
                write(3,'(i7,i7)') 1, i
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'CELL_TYPES  '//trim(it2)//''
            do i=1, grid%nelem
                write(3,'(i7)') 9 ! 9 for quadrilaterals
            end do
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(i7)') 1 ! 1 for vertex
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'POINT_DATA  '//trim(it)//''
            write(3,'(a)') 'SCALARS  vx_strength  double'
            write(3,'(a)') 'LOOKUP_TABLE default'
            do i = 1,grid%nnod
                write(3,'(1(1x,e14.6))') 0._pr
            end do
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(1(1x,e14.6))') wake%multicirc_mag(i)
            end do
            write(3,'(a)') ''
            write(3,'(a)') 'SCALARS  vx_core_radius  double'
            write(3,'(a)') 'LOOKUP_TABLE default'
            do i = 1,grid%nnod
                write(3,'(1(1x,e14.6))') 0._pr
            end do
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(1(1x,e14.6))') wake%multivxcore_rad(i) 
            end do
            end select
        close(3)  
    
End Subroutine OUTPUT_PLATE

Subroutine OUTPUT_WAKE (step) ! writes the output results file (for the vortons' wake) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
    Implicit none
    integer :: i,totnod,totpan
    integer, intent(in) :: step
    character(12) :: it,it2,st
    
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(4,file='C:\Users\pimen\Documents\Visual Studio 2017\VORTONEX_v1\VortoNeX\output_vtk\wakevort_'//trim(st)//'.vtk')
        write(4,'(a)')'# vtk DataFile Version 2.0'
        write(4,'(a)') 'This file was created by Fortran'
        write(4,'(a)') 'ASCII'
        write(4,'(a)') 'DATASET UNSTRUCTURED_GRID'
        write(it, '(i12)') wake%nwp*step
        it=adjustl(it)        
        write(4,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (4)
            do i = 1, wake%nwp*step ! from 1 to 40,80,120... 
                write(4,'(3(1x,e14.6))') wake%pos(1,i+wake%nwp), wake%pos(2,i+wake%nwp), wake%pos(3,i+wake%nwp)
            end do
            write(4,'(a)') ''
            write(it2, '(i12)') 2*wake%nwp*step
            it2=adjustl(it2)        
            write(4,'(a)') 'CELLS  '//trim(it)//'  '//trim(it2)//''
            
            do i = 1, wake%nwp*step ! from 1 to 40,80,120...
                write(4,'(i7,i7)') 1, i-1
            end do
            
            write(4,'(a)') ''
            write(4,'(a)') 'CELL_TYPES  '//trim(it)//''
            
            do i = 1, wake%nwp*step ! from 1 to 40,80,120...
                write(4,'(i7)') 1 ! 1 for vertex
            end do
            
            write(4,'(a)') ''
            write(4,'(a)') 'POINT_DATA  '//trim(it)//''
            write(4,'(a)') 'SCALARS  vx_strength  double'
            write(4,'(a)') 'LOOKUP_TABLE default'
            
            do i=wake%nwp+1, wake%nwp*(step+1)
                write(4,'(1(1x,e14.6))') wake%gamma_mag(i)
            end do   
            
            write(4,'(a)') ''
            write(4,'(a)') 'SCALARS  vx_core_radius  double'
            write(4,'(a)') 'LOOKUP_TABLE default'
            
            do i=wake%nwp+1, wake%nwp*(step+1)               
                write(4,'(1(1x,e14.6))') wake%vxcore_rad(i) 
            end do       
            
            end select
        close(4)  
    
End Subroutine OUTPUT_WAKE

Subroutine OUTPUT_FIL (step) ! writes the output results file (for the wake grid) for Paraview; JCPG
    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
    Implicit none
    integer :: i,totnod,totpan,cont
    integer, intent(in) :: step
    character(12) :: it,it2,it3,st
    
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    write(st,'(i12)') step
    st = adjustl(st)
        open(5,file='C:\Users\pimen\Documents\Visual Studio 2017\VORTONEX_v1\VortoNeX\output_vtk\wakefil_'//trim(st)//'.vtk')
        write(5,'(a)')'# vtk DataFile Version 2.0'
        write(5,'(a)') 'This file was created by Fortran'
        write(5,'(a)') 'ASCII'
        write(5,'(a)') 'DATASET POLYDATA'
        write(it, '(i12)') 2*wake%nwp*step
        it=adjustl(it)        
        write(5,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'

        select case(grid%elemtype)
        case (4)
            do i = 1, wake%nwp*step ! from 1 to 40,80,120... 
                write(5,'(3(1x,e14.6))') wake%pos_fst_point(1,i+wake%nwp), wake%pos_fst_point(2,i+wake%nwp), wake%pos_fst_point(3,i+wake%nwp)
                write(5,'(3(1x,e14.6))') wake%pos_sec_point(1,i+wake%nwp), wake%pos_sec_point(2,i+wake%nwp), wake%pos_sec_point(3,i+wake%nwp)
            end do
            
            write(5,'(a)') ''
            write(it2, '(i12)') wake%nwp*step
            write(it3, '(i12)') 3*wake%nwp*step
            it2=adjustl(it2); it3=adjustl(it3)
            write(5,'(a)') 'LINES  '//trim(it2)//'  '//trim(it3)//''
            
            cont=1
            do i = 1, wake%nwp*step ! from 1 to 40,80,120...
                write(5,'(i7,i7,i7)') 2, cont-1, cont
                cont = cont + 2
            end do
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'CELL_TYPES  '//trim(it2)//''
            !
            !do i = 1, wake%nwp*step ! from 1 to 40,80,120...
            !    write(5,'(i7)') 2 ! 1 for vertex
            !end do
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'POINT_DATA  '//trim(it2)//''
            !write(5,'(a)') 'SCALARS  vx_strength  double'
            !write(5,'(a)') 'LOOKUP_TABLE default'
            !
            !do i=wake%nwp+1, wake%nwp*(step+1)
            !    write(5,'(1(1x,e14.6))') wake%gamma_mag(i)
            !end do   
            
            !write(5,'(a)') ''
            !write(5,'(a)') 'SCALARS  vx_core_radius  double'
            !write(5,'(a)') 'LOOKUP_TABLE default'
            !
            !do i=wake%nwp+1, wake%nwp*(step+1)               
            !    write(5,'(1(1x,e14.6))') wake%vxcore_rad(i) 
            !end do       
            !
            end select
        close(5)  
    
End Subroutine OUTPUT_FIL

! UNUSED SUBROUTINES (they have not been deleted for reference or ongoing development)
Subroutine VEL_VORT2POINTVORT (circ_vec,point,pos_kvort,vind_v2p,i,k,option) ! calculates the induced velocity by a vorton over another one (for symmetrized advection); it does not apply cause the advection is performed on points (filament's endpoints); JCPG
    Use ARRAYS, only : pr,input,wake,fiveovertwo,fourpi,fouroverpi,e_num,twopi,pi,oneoverthree,twooverpi
    Use MATH, only : VECTORNORM,CROSSPRODUCT,NORMALIZE
    implicit none
    real(kind=pr) :: r_mag,r_mag2,r_mag3,fac1,fac2,fac3,cons,g_winck,rho,rho2,rho3,minor_radius,sym_core_rad
    real(kind=pr) :: r(3),winck(3)
    real(kind=pr), intent(in) :: circ_vec(3),point(3),pos_kvort(3)
    real(kind=pr), intent(out) :: vind_v2p(3)
    integer, intent(in) :: i,k,option

    r(:) = point(:) - pos_kvort(:) ! radio vector between k-vorton and an evaluation point
    CALL VECTORNORM(r(:), r_mag) ! vector norm (distance between points)
    minor_radius = wake%vxcore_rad(i) / (2._pr**oneoverthree) ! equivalent radius for half volume
    
    if (option==0) then ! for plate
        sym_core_rad = sqrt( ( minor_radius*minor_radius + wake%multivxcore_rad(k)*wake%multivxcore_rad(k) ) / 2._pr ) ! symmetrized vortex core radius according to (HE, 2009).
        wake%core_rad2 = sym_core_rad*sym_core_rad ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*sym_core_rad ! cubic vortex core radius
        rho = r_mag / sym_core_rad
    else if (option==1) then ! for wake (acts on points)
        sym_core_rad = sqrt( ( minor_radius*minor_radius + wake%vxcore_rad(k)*wake%vxcore_rad(k) ) / 2._pr ) ! symmetrized vortex core radius according to (HE, 2009).
        wake%core_rad2 = sym_core_rad*sym_core_rad ! squared vortex core radius
        wake%core_rad3 = wake%core_rad2*sym_core_rad ! cubic vortex core radius
        rho = r_mag / sym_core_rad 
    else
        pause
    end if
        
    r_mag2 = r_mag*r_mag ! squared distance
    r_mag3 = r_mag2*r_mag ! cubic distance
    rho2 = rho*rho; rho3=rho*rho*rho
        
    if (r_mag<=input%tol_rad) then ! vortons are too close
        vind_v2p(:) = 0._pr ! induced velocity is zero
    else if (input%core_rad_init == 0._pr) then ! singularized vortex core radius; it does not apply
    else ! regularized vortex core
        select case (input%regul_function)
            case (0) ! High-order RF
                g_winck = (rho3 * (rho2 + fiveovertwo)) / ((rho2 + 1._pr)**fiveovertwo) ! High-order RF; Winckelmans
            case (1) ! 2nd-order Gaussian RF (BERDOWSKY, 2015) 
                g_winck = 1._pr - e_num**(-rho3)
            case (2) ! Gaussian error function (used by several authors)
                cons = 0.5_pr*rho2 !r_mag2/(2._pr*wake%core_rad2)     
                !g_winck = erf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons)
                g_winck = derf(rho/sqrt(2._pr)) - rho*sqrt(twooverpi)*e_num**(-cons) ! double precision erf
        end select
                
        fac1 = -1._pr/fourpi
        fac2 = g_winck / r_mag3 !  (BIRD, 2021)
        fac3 = fac1 * fac2 ! first times second factors
        winck(:) = fac3 * r(:) ! it works for (WINCKELMANS, 2005) and (ÁLVAREZ, 2018); same result
        call CROSSPRODUCT(winck(:), circ_vec(:), vind_v2p(:)) ! induced velocity by a vorton over another one (or over an evaluation point)
    end if

End Subroutine VEL_VORT2POINTVORT

Subroutine FORCES_GUTNIKOV ( point, tot_area, circ_vec, den_aux, option, pos_kvort) ! hydrodynamic coefficients calculation via Gutnikov's scheme (GUTNIKOV, 2006); TESTING YET!; JCPG
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,ctrl_pts,area,deriv_gamma,nor_vec,del_cp,del_pres,del_force,spanwise_wake_elem,chordwise_wake_elem,bound_circ_old,elem2seg
    Use MATH, only : CROSSPRODUCT,DOTPRODUCT
    Implicit none

    integer :: k,fst_nod,sec_nod,m,i,aa,s,w,cont
    integer, intent(in) :: den_aux
    integer, intent(out) :: option
    real(kind=pr) :: den,cfx,cfy,cfz,cl,cd,force_x,force_y,force_z!,cd
    real(kind=pr) :: veltimesdeltavel,alpharad,sinalpha,cosalpha
    real(kind=pr) :: vind_body(3),vind_total(3),gamma_vector(3),parenth(3),vind_v2p(3),vind_cloud(3),ds(3)
    real(kind=pr), intent(in) :: tot_area
    real(kind=pr), intent(out) :: point(3),circ_vec(3),pos_kvort(3)
    logical :: logic(wake%nwp)
    
    logic(:) = .false.
    do k=1, grid%nelem ! over BVR/panels
        deriv_gamma(k) = ( bound_circ(k) - bound_circ_old(k) ) / input%dt ! Gamma time derivative
        point(:) = ctrl_pts(:,k) ! over control point (panel's center)

        vind_body(:) = 0._pr
        option=0 ! for plate case
        !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !$OMP SHARED(den_aux,point,option,wake,k) PRIVATE(s,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
        !$OMP REDUCTION(+: vind_body) ! (declaration of reduction variables)
        !$OMP DO
        do s=1, den_aux
            circ_vec(:) = wake%multicirc_vec(:,s)
            pos_kvort(:) = wake%pos_multifix2plate(:,s) ! k-vorton's position
            call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,s,option) ! calculates the induced velocity by a vorton over another one
            vind_body(:) = vind_body(:) + vind_v2p(:) ! induced velocity by the fixed vortons over a single one                        
        end do 
        !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !$OMP end PARALLEL
            
        vind_cloud(:) = 0._pr ! clears old values
        option=1 ! for wake case
        !$OMP PARALLEL IF(input%nthreads>1) NUM_THREADS(input%nthreads) DEFAULT(none) & ! begins a parallel section
        !$OMP SHARED(point,option,wake,k) PRIVATE(s,circ_vec,pos_kvort,vind_v2p) & ! (declaration of variables)
        !$OMP REDUCTION(+: vind_cloud) ! (declaration of reduction variables)
        !$OMP DO
        do s = wake%nwp+1, wake%nwp + wake%num_vort ! over the free vortons
            circ_vec(:) = wake%gammav(:,s) ! vorton's circulation vector
            pos_kvort(:) = wake%pos(:,s) ! k-vorton's position
            call VEL_VORT2POINT(circ_vec,point,pos_kvort,vind_v2p,s,option) ! calculates the induced velocity by a vorton over another one
            vind_cloud(:) = vind_cloud(:) + vind_v2p(:) ! induced velocity by the vorton cloud over a single one
        end do ! ends k  
        !$OMP ENDDO !(!$OMP ENDDO NOWAIT allows a thread to continue without others finish)
        !$OMP end PARALLEL
            
        vind_total(:) = vdir(:) + vind_cloud(:) + vind_body(:) ! total local flow velocity

        do i=1, 4 ! for each panel's edge (four filaments)
            if (i==3 .or. i==4) then
                if (i==3) then
                    fst_nod = grid%panel(i+1,k) ! first node (it depends on mesh numeration's starting point)
                    sec_nod = grid%panel(i,k) ! second node (it depends on mesh numeration's starting point)
                else ! i==4
                    fst_nod = grid%panel(1,k) ! first node (it depends on mesh numeration's starting point)
                    sec_nod = grid%panel(4,k) ! second node (it depends on mesh numeration's starting point)
                end if
            else !(i==1 .or. i==2)
                fst_nod = grid%panel(i,k) ! second node (it depends on mesh numeration's starting point)
                sec_nod = grid%panel(i+1,k) ! first node (it depends on mesh numeration's starting point)
            end if
                
            if (i==1 .or. i==3) then ! for spanwise edges 
                cont=0 ! to save computation
                do m=1, wake%nwp - (grid%nle + grid%plate_nle + grid%plate_nte)
                    aa=spanwise_wake_elem(m) - grid%nelem ! wake element (1,3,5,7,10,12 for 2x2 mesh)
                    select case(input%mesh_num)
                    case(1,2,4,5,7) ! 4x4 or 16x16
                        if ( fst_nod == grid%panel(2,spanwise_wake_elem(m)) .and. sec_nod == grid%panel(1,spanwise_wake_elem(m)) .and. cont==0 ) then !.and. logic(aa)==.false.) then ! for 4x4 and 16x16 case
                            cont=cont+1
                            if (i==1) then
                                ds(:) = grid%coord(:,grid%panel(2,aa+grid%nelem)) - grid%coord(:,grid%panel(1,aa+grid%nelem)) ! segment's vector
                            else !(i==3)
                                ds(:) = grid%coord(:,grid%panel(1,aa+grid%nelem)) - grid%coord(:,grid%panel(2,aa+grid%nelem)) ! segment's vector
                            end if
                            do w=1,2
                                if (elem2seg(w,aa)/=0) then
                                    wake%gammav_Gutnikov(:,i) = bound_circ(elem2seg(w,aa))*ds(:) ! inverts circulation for lateral (internal and external) edges
                                    EXIT
                                else ! is zero (segment does not belongs to this panel)
                                end if
                            end do ! ends w
                        else
                        end if
                    case(3,6) ! 10x10
                        if ( fst_nod == grid%panel(2,spanwise_wake_elem(m)) .and. sec_nod == grid%panel(1,spanwise_wake_elem(m)) .and. cont==0 ) then !.and. logic(aa)==.false.) then ! for 10x10 mesh
                            cont=cont+1
                            if (i==1) then
                                ds(:) = grid%coord(:,grid%panel(2,aa+grid%nelem)) - grid%coord(:,grid%panel(1,aa+grid%nelem)) ! segment's vector
                            else !(i==3)
                                ds(:) = grid%coord(:,grid%panel(1,aa+grid%nelem)) - grid%coord(:,grid%panel(2,aa+grid%nelem)) ! segment's vector
                            end if
                            do w=1,2
                                if (elem2seg(w,aa)/=0) then
                                    wake%gammav_Gutnikov(:,i) = bound_circ(elem2seg(w,aa))*ds(:) ! inverts circulation for lateral (internal and external) edges
                                    EXIT
                                else ! is zero (segment does not belongs to this panel)
                                end if
                            end do ! ends w
                        else
                        end if
                    end select
                end do ! ends m

            else if (i==2 .or. i==4) then ! for chordwise edges
                cont=0
                do m=1, grid%nle + grid%plate_nle + grid%plate_nte
                    aa=chordwise_wake_elem(m) - grid%nelem ! wake element (2,4,6,8,9,11 for 2x2 mesh)
                    select case(input%mesh_num)
                    case(1,2,4,5,7) ! 4x4 or 16x16
                        if ( fst_nod == grid%panel(2,chordwise_wake_elem(m)) .and. sec_nod == grid%panel(1,chordwise_wake_elem(m)) .and. cont==0 ) then !.and. logic(aa)==.false.) then ! for 4x4 and 16x16 case
                            cont=cont+1
                            if (i==2) then
                                ds(:) = grid%coord(:,grid%panel(2,aa+grid%nelem)) - grid%coord(:,grid%panel(1,aa+grid%nelem)) ! segment's vector
                            else !(i==4)
                                ds(:) = grid%coord(:,grid%panel(1,aa+grid%nelem)) - grid%coord(:,grid%panel(2,aa+grid%nelem)) ! segment's vector
                            end if
                            do w=1,2
                                if (elem2seg(w,aa)/=0) then
                                    wake%gammav_Gutnikov(:,i) = bound_circ(elem2seg(w,aa))*ds(:) ! inverts circulation for lateral (internal and external) edges
                                    EXIT
                                else ! is zero (segment does not belongs to this panel)
                                end if
                            end do ! ends w
                        else
                        end if
                    case(3,6) ! 10x10
                        if ( fst_nod == grid%panel(1,chordwise_wake_elem(m)) .and. sec_nod == grid%panel(2,chordwise_wake_elem(m)) .and. cont==0 ) then !.and. logic(aa)==.false.) then ! for 10x10 mesh
                            cont=cont+1
                            if (i==2) then
                                ds(:) = grid%coord(:,grid%panel(1,aa+grid%nelem)) - grid%coord(:,grid%panel(2,aa+grid%nelem)) ! segment's vector
                            else !(i==4)
                                ds(:) = grid%coord(:,grid%panel(2,aa+grid%nelem)) - grid%coord(:,grid%panel(1,aa+grid%nelem)) ! segment's vector
                            end if
                            do w=1,2
                                if (elem2seg(w,aa)/=0) then
                                    wake%gammav_Gutnikov(:,i) = bound_circ(elem2seg(w,aa))*ds(:) ! inverts circulation for lateral (internal and external) edges
                                    EXIT
                                else ! is zero (segment does not belongs to this panel)
                                end if
                            end do ! ends w
                        else
                        end if
                    end select
                end do ! ends m
            end if
        end do ! ends i
        
        gamma_vector(:) = ( wake%gammav_Gutnikov(:,1) - wake%gammav_Gutnikov(:,2) + wake%gammav_Gutnikov(:,3) - wake%gammav_Gutnikov(:,4) ) / area(k) ! total circulation divided by the panel's area
        call CROSSPRODUCT(vind_total(:), gamma_vector(:), parenth(:)) ! parenthesis of eq. 3.4 in (GUTNIKOV, 2006)
        call DOTPRODUCT(parenth(:), nor_vec(:,k), veltimesdeltavel) ! velocity times delta velocity (W*deltaW)
        del_pres(k) = input%dens*veltimesdeltavel - input%dens*deriv_gamma(k) ! jump pressure across the panel
        del_cp(k) = 2._pr*del_pres(k) / (input%dens * input%q_inf*input%q_inf) ! pressure coefficient jump across the panel            
    end do ! ends k

    alpharad = input%alpha*(pi/180._pr)!; betarad  = input%beta*(pi/180._pr)  ! angles in radians
    sinalpha = sin(alpharad); cosalpha = cos(alpharad)!; sinbeta = sin(betarad); cosbeta = cos(betarad)

    force_x=0._pr; force_y=0._pr; force_z=0._pr
    do k=1, grid%nelem
        del_force(:,k) = -(del_pres(k)*area(k))*nor_vec(:,k) ! force per panel
        force_x = force_x + del_force(1,k) ! longitudinal force
        force_y = force_y + del_force(2,k) ! lateral force
        force_z = force_z + del_force(3,k) ! vertical force        
    end do

    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
    cfx = 2._pr*force_x / den; cfy = 2._pr*force_y / den; cfz = 2._pr*force_z / den
    
    cfx = (1._pr/input%char_len)*cfx  ! longitudinal force coefficient
    cfy = (1._pr/input%char_len)*cfy  ! lateral force coefficient
    cfz = (1._pr/input%char_len)*cfz  ! vertical force coefficient
    
    cl = cfz*cosalpha - cfx*sinalpha  ! lift coefficient (CL)
    cd = cfz*sinalpha + cfx*cosalpha  ! drag coefficient (CD)
    !cy  = -cfx*cosalpha + cfz*sinalpha ! lateral force coef. (CY)   
    
    400 format(f12.6)
    print 400, cl; print 400, cd!; print 400, cm; !print*, cn
    print *,'---'
    
    500 format(f12.6,f12.6,f12.6)
    open(4, file='aero_coef.dat', access='append') ! writes output file
    write(4,500) cl, cd !, cm
    close(4)

End Subroutine FORCES_GUTNIKOV
    
Subroutine OUTPUT_MSH (step,den_aux) ! writes output mesh file for GiD's postprocessor; JCPG/EOA
    Use ARRAYS, only : grid,kuttaedges,input,wake
    Implicit none
    integer :: i,j,totnod,totpan
    integer, intent(in) :: step,den_aux
    character(12) :: it
    
    totnod = grid%sizen
    totpan = grid%nelem + kuttaedges%nte*input%nsteps
    
    if (step==1) then ! for first time step
        open(3,file='test.post.msh')
        write(3,'(a,/)')'# encoding utf-8'
        write(3,*) ''
        write(3,*) 'Group "unsteady_1"'
        write(3,*) ''
    
        select case(grid%elemtype)
        case (4)
            write(3,*)'MESH body dimension 3 ElemType Quadrilateral Nnode 4' ! body panels
            write(3,'(a)')'Coordinates'
            do i = 1,grid%nnod
                write(3,'(i7,3(1x,e14.6))') i, grid%coord(1:3,i)
            end do
            write(3,'(a)')'End coordinates'
            write(3,'(a)')'Elements'
            do i = 1,grid%nelem
                write(3,'(5(1x,i7))') i , (grid%panel(j,i),j=1,4)  
            end do
            write(3,'(a)')'End elements'
            write(3,*) ''
            
            write(3,*)'MESH bodyblobs dimension 2 Elemtype Sphere nnode 1' ! multi-bounded blobs
            write(3,'(a)')'Coordinates'
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(i7,3(1x,e14.6))') i+grid%nnod, wake%pos_multifix2plate(1,i), wake%pos_multifix2plate(2,i), wake%pos_multifix2plate(3,i)
            end do
            write(3,'(a)')'End coordinates'
            write(3,'(a)')'Elements'
            do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                write(3,'(i7,i7,1(1x,e14.6))') grid%nnod+i, grid%nnod+i, wake%multivxcore_rad(i) 
            end do
            write(3,'(a)')'End elements'
            write(3,*) ''
    
            write(3,*)'MESH wake dimension 2 Elemtype Sphere nnode 1'
            write(3,'(a)')'Coordinates'
            do i = 1, wake%nwp*step ! from 1 to 40 
                write(3,'(i7,3(1x,e14.6))') i+grid%nnod+den_aux, wake%pos(1,i+wake%nwp), wake%pos(2,i+wake%nwp), wake%pos(3,i+wake%nwp)
            end do
            write(3,'(a)')'End coordinates'
            write(3,'(a)')'Elements'
            do i = 1, wake%nwp*step ! from 1 to 40 
                write(3,'(i7,i7,1(1x,e14.6))') grid%nnod+i+den_aux, grid%nnod+i+den_aux, wake%vxcore_rad(i) ! 17,18,...,56   
            end do
            write(3,'(a)')'End elements'
            write(3,*) ''
            write(3,*) 'End Group'
            write(3,*) ''
        end select
        
        else ! step>=2
            open(3,file='test.post.msh', access = 'append')
            write(it, '(i12)') step
            it=adjustl(it)
            write(3,*) 'group '//'"unsteady_'//trim(it)//'"'
            write(3,*) ''
    
            select case(grid%elemtype)
            case (4)
                write(3,*)'MESH body dimension 3 ElemType Quadrilateral Nnode 4' ! body panels
                write(3,'(a)')'Coordinates'
                do i = 1,grid%nnod
                    write(3,'(i7,3(1x,e14.6))') i, grid%coord(1:3,i)
                end do
                write(3,'(a)')'End coordinates'
                write(3,'(a)')'Elements'
                do i = 1,grid%nelem
                    write(3,'(5(1x,i7))') i , (grid%panel(j,i),j=1,4)
                end do
                write(3,'(a)')'End elements'
                write(3,*) ''
                
                write(3,*)'MESH bodyblobs dimension 2 Elemtype Sphere nnode 1' ! multi-bounded blobs
                write(3,'(a)')'Coordinates'
                !!$OMP DO
                do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                    write(3,'(i7,3(1x,e14.6))') i+grid%nnod, wake%pos_multifix2plate(1,i), wake%pos_multifix2plate(2,i), wake%pos_multifix2plate(3,i)
                end do
                !!$OMP END DO
                write(3,'(a)')'End coordinates'
                write(3,'(a)')'Elements'
                do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
                    write(3,'(i7,i7,1(1x,e14.6))') grid%nnod+i, grid%nnod+i, wake%multivxcore_rad(i) 
                end do
                write(3,'(a)')'End elements'
                write(3,*) ''
    
                write(3,*)'MESH wake dimension 2 Elemtype Sphere nnode 1'
                write(3,'(a)')'Coordinates'
                do i = 1, wake%nwp*step ! from 1 to 40,80,120... 
                    write(3,'(i7,3(1x,e14.6))') i+grid%nnod+den_aux, wake%pos(1,i+wake%nwp), wake%pos(2,i+wake%nwp), wake%pos(3,i+wake%nwp)
                end do
                write(3,'(a)')'End coordinates'
                write(3,'(a)')'Elements'
                do i=wake%nwp+1, wake%nwp*(step+1)
                    write(3,'(i7,i7,1(1x,e14.6))') grid%nnod+i-wake%nwp+den_aux,grid%nnod+i-wake%nwp+den_aux, wake%vxcore_rad(i) ! 17,18,...,56   
                end do
                write(3,'(a)')'End elements'
                write(3,*) ''
                write(3,*) 'End Group'
                write(3,*) ''        
            end select
        end if
                
        close(3)  
    
End Subroutine OUTPUT_MSH

Subroutine OUTPUT_RES (step,den_aux) ! writes output results file for GiD's postprocessor; JCPG/EOA
    Use ARRAYS, only: pr,grid,input,wake
    implicit none
    integer, intent(in) :: step,den_aux
    integer :: i
    character(24) :: group_name
    character(12) :: it    
    write(it,'(i12)') step ! num_step
    it = adjustl(it)   
    group_name = 'Group '//'"unsteady'//'_'//trim(it)//'"' 
   
    if (step==1) then
        open (4,file='test.post.res')
        write(4,*) 'GiD Post Results File 1.0'
        write(4,'(a,/)')'# encoding utf-8'
        write(4,*) ''
        write(4,'(a,/)')'On'//group_name
        write(4,*) ''
        write(4,*) 'result "vx strength" "load analysis"',input%dt*step,'scalar onnodes'
        write(4,*) 'values'
        do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
            write(4,'(i7,1(1x,e14.6))') grid%nnod+i, wake%multicirc_mag(i)
        end do        
        do i=wake%nwp+1, wake%nwp*(step+1) ! 41 -> 80
            write(4,'(i7,3(1x,e14.6))') grid%nnod+i-wake%nwp+den_aux, wake%gamma_mag(i) ! circulation strength
        end do
        write(4,*) 'end values'
        write(4,*) ''        
        write(4,*) 'result "vx core radius" "load analysis"',input%dt*step,'scalar onnodes'
        write(4,*) 'values'
        do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 votons per vortex segment)
            write(4,'(i7,1(1x,e14.6))') grid%nnod+i, wake%multivxcore_rad(i) ! value remains constant along the simulation
        end do
        do i=wake%nwp+1, wake%nwp*(step+1)
            write(4,'(i7,3(1x,e14.6))') grid%nnod+i-wake%nwp+den_aux, wake%vxcore_rad(i) ! vortex core radius
        end do
        write(4,*) 'end values'
        write(4,*) ''        
        write(4,*) 'end ongroup'
        
    else
        open (4,file='test.post.res', access = 'append')
        write(it,'(i12)') step
        it=adjustl(it)
        write(4,'(a,/)')'On'//group_name
        write(4,*) ''
        write(4,*) 'result "vx strength" "load analysis"',input%dt*step,'scalar onnodes'
        write(4,*) 'values'
        do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
            write(4,'(i7,1(1x,e14.6))') grid%nnod+i, wake%multicirc_mag(i)
        end do 
        do i=wake%nwp+1, wake%nwp*(step+1)
            write(4,'(i7,3(1x,e14.6))') grid%nnod+i-wake%nwp+den_aux, wake%gamma_mag(i) ! circulation strength
        end do
        write(4,*) 'end values'
        write(4,*) ''        
        write(4,*) 'result "vx core radius" "load analysis"',input%dt*step,'scalar onnodes'
        write(4,*) 'values'
        do i = 1, den_aux ! from 1 to 120 (4x4 discretization; 3 vortons per vortex segment)
            write(4,'(i7,1(1x,e14.6))') grid%nnod+i, wake%multivxcore_rad(i) 
        end do
        do i=wake%nwp+1, wake%nwp*(step+1)
            write(4,'(i7,3(1x,e14.6))') grid%nnod+i-wake%nwp+den_aux, wake%vxcore_rad(i) ! vortex core radius 
        end do
        write(4,*) 'end values'
        write(4,*) ''   
        write(4,*) 'end ongroup'
        
    end if
    
    close(4)
End Subroutine OUTPUT_RES

!Subroutine OUTPUT_TRAJ (step) ! writes the output results file (for the trajectories) for Paraview (not finished yet!); JCPG
!    Use ARRAYS, only : grid,kuttaedges,input,wake,pr
!    Implicit none
!    integer :: i,totnod,totpan,cont
!    integer, intent(in) :: step
!    character(12) :: it,it2,it3,st
!    
!    totnod = grid%sizen
!    totpan = grid%nelem + kuttaedges%nte*input%nsteps
!    write(st,'(i12)') step
!    st = adjustl(st)
!        open(5,file='C:\Users\pimen\Documents\Visual Studio 2017\VORTONEX_v1\VortoNeX\output_vtk\wakevort_'//trim(st)//'.vtk')
!        write(5,'(a)')'# vtk DataFile Version 2.0'
!        write(5,'(a)') 'This file was created by Fortran'
!        write(5,'(a)') 'ASCII'
!        write(5,'(a)') 'DATASET POLYDATA'
!        write(it, '(i12)') 2*wake%nwp*step*2
!        it=adjustl(it)        
!        write(5,'(a)') 'POINTS  '//trim(it)//'  DOUBLE'
!
!        select case(grid%elemtype)
!        case (4)
!            
!            do i = 1, wake%nwp*step ! from 1 to 40,80,120... 
!                write(5,'(3(1x,e14.6))') wake%pos_old_fst(1,i), wake%pos_old_fst(2,i), wake%pos_old_fst(3,i)
!                write(5,'(3(1x,e14.6))') wake%pos_old_sec(1,i), wake%pos_old_sec(2,i), wake%pos_old_sec(3,i)
!            end do
!            
!            do i = 1, wake%nwp*step ! from 1 to 40,80,120... 
!                write(5,'(3(1x,e14.6))') wake%pos_fst_point(1,i+wake%nwp), wake%pos_fst_point(2,i+wake%nwp), wake%pos_fst_point(3,i+wake%nwp)
!                write(5,'(3(1x,e14.6))') wake%pos_sec_point(1,i+wake%nwp), wake%pos_sec_point(2,i+wake%nwp), wake%pos_sec_point(3,i+wake%nwp)
!            end do
!            
!            write(5,'(a)') ''
!            write(it2, '(i12)') wake%nwp*step*2
!            write(it3, '(i12)') 3*wake%nwp*step*2
!            it2=adjustl(it2); it3=adjustl(it3)
!            write(5,'(a)') 'LINES  '//trim(it2)//'  '//trim(it3)//''
!            
!            !!cont=1
!            !do i = 1, wake%nwp*step ! from 1 to 40,80,120...
!            !    write(5,'(i7,i7,i7)') 2, cont-1, cont !i+wake%nwp*step
!            !    !cont = cont + 2
!            !end do
!            
!            cont = 1
!            do i = 1, wake%nwp*step*2
!                !write(5,'(i7,i7,i7)') 2, i-1, i-1+wake%nwp
!                write(5,'(i7,i7,i7)') 2, cont-1, cont + wake%nwp*step*2 - 1
!                cont = cont + 1
!            end do
!            
!            write(5,'(a)') ''
!            write(5,'(a)') 'CELL_TYPES  '//trim(it2)//''
!            
!            do i = 1, wake%nwp*step*2 ! from 1 to 40,80,120...
!                write(5,'(i7)') 2 ! 1 for vertex
!            end do
!            
!            !write(4,'(a)') ''
!            !write(4,'(a)') 'POINT_DATA  '//trim(it2)//''
!            !write(4,'(a)') 'SCALARS  vx_strength  double'
!            !write(4,'(a)') 'LOOKUP_TABLE default'
!            !
!            !do i=wake%nwp+1, wake%nwp*(step+1)
!            !    write(4,'(1(1x,e14.6))') wake%gamma_mag(i)
!            !end do   
!            
!            !write(4,'(a)') ''
!            !write(4,'(a)') 'SCALARS  vx_core_radius  double'
!            !write(4,'(a)') 'LOOKUP_TABLE default'
!            !
!            !do i=wake%nwp+1, wake%nwp*(step+1)               
!            !    write(4,'(1(1x,e14.6))') wake%vxcore_rad(i) 
!            !end do       
!            !
!            end select
!        close(5)  
!    
!End Subroutine OUTPUT_TRAJ

Subroutine COUPLING_BODY( circ, j, point, tot_area ) ! generates the body's influence coefficient matrix (ICM) for the BVRs scheme (V1); JCPG
    Use ARRAYS, only : pr,ctrl_pts,nor_vec,A_body,grid,rhs,pivot,A_solv,bound_circ,vdir,input,tan_vec1,tan_vec2,area,bound_circ_old,deriv_gamma,rhs_freeflow
    Use MATH, only : DOTPRODUCT
    Implicit none
    integer :: i,op,step,k,r,fstlen,stret
    integer, intent(out) :: j
    real(kind=pr), intent(out) :: circ,point(3),tot_area
    real(kind=pr) :: vind_ring(3),ra(3),rb(3),pco(3),geometry(10)!,vind_v2p(3)
    
    allocate(ctrl_pts(3,grid%nelem),nor_vec(3,grid%nelem),A_body(grid%nelem,grid%nelem),rhs(grid%nelem),tan_vec1(3,grid%nelem),tan_vec2(3,grid%nelem), area(grid%nelem),rhs_freeflow(grid%nelem))
    
    vdir(1:3) = (/ input%q_inf*cosd(input%alpha), 0._pr , input%q_inf*sind(input%alpha) /)  ! total velocity ( Vwind-Vbody, viewed from the body)
    grid%elemtype=grid%pnod(1)
    
    do i=1, grid%nelem
        call PANELGEO(geometry,pco,i)
        ctrl_pts(:,i) = pco ! panel's control points on the body
        
        tan_vec1(:,i) = geometry(1:3) ! first tangential vector
        tan_vec2(:,i) = geometry(4:6) ! second tangential vector
        nor_vec(:,i)  = geometry(7:9) ! normal vector
        area(i)       = geometry(10)  ! panel's area
    end do
    
    tot_area = 0._pr
    do i=1, grid%nelem
        tot_area = tot_area + area(i) ! total body's area
    end do
    
    ! Shell-body's influence coefficients' matrix (A_body; only bounded VR)
    do i=1, grid%nelem ! on control points
        point(:) = ctrl_pts(:,i)
        do j=1, grid%nelem ! on bounded VR
            circ = 1.0_pr ! unitary VR's circulation for starting solution
            k=0; r=0; step=0; op=0; fstlen=0; stret=0 ! such values comes from previous UFVLM code
            call VORING ( j,point,circ,vind_ring, ra, rb,k,r,step,op,fstlen,stret) ! vortex ring induced velocity calculation
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),A_body(i,j)) ! ICM's elements
        end do ! ends j
    end do ! ends i
        
    allocate( A_solv(grid%nelem,grid%nelem),pivot(grid%nelem),bound_circ(grid%nelem),bound_circ_old(grid%nelem),deriv_gamma(grid%nelem) )

End Subroutine COUPLING_BODY 

Subroutine STARTING_SOL( circ, n, point, tot_area,j,op,fstlen,stret,step ) ! solves the system of equations to obtain a starting solution for the BVRs' scheme (V1); JCPG
    Use UTILS, only : RESIZE_REAL_1
    Use ARRAYS, only : grid,ctrl_pts,wake,kuttaedges,pr,A_body,nor_vec,A_first,vdir,rhs,A_solv,pivot,bound_circ,A_total,gamma_wake,A_multiwake,esed1_reduc,orien_reduc,internal_LE_detected,plate_LE_detected,fstwake_len,gamma_wake_start,esed1_modif_Gutnikov,del_pres,del_cp,del_force,gamma_wake_mod,input,rhs_freeflow,del_pres_down
    Use MATH, only : DOTPRODUCT
    Implicit none
    integer :: i,k,trailpan,ok,p,m
    integer, intent(out) :: j
    integer :: staux,jaux
    real(kind=pr) :: vind_ring(3),ra(3),rb(3),mid_point(3)!,sum_circ_plate
    real(kind=pr), intent(inout) :: tot_area
    real(kind=pr), intent(out) :: circ,point(3)
    integer, intent(out) :: n,op,fstlen,stret,step
    real(kind=pr) :: dot
    logical :: log1,log2

    allocate( A_first(grid%nelem,grid%nelem), A_total(grid%nelem,grid%nelem), A_multiwake(grid%nelem,grid%nelem) )
    allocate( gamma_wake(wake%nwp),esed1_reduc(wake%nwp),orien_reduc(wake%nwp), fstwake_len(4,wake%nwp), gamma_wake_start(wake%nwp),esed1_modif_Gutnikov(4,grid%nelem),del_pres(grid%nelem),del_cp(grid%nelem),del_force(3,grid%nelem),gamma_wake_mod(wake%nwp),del_pres_down(grid%nelem) )

    ! Multi-wake effect to the A (total) matrix; A_total = A_body + A_multiwake
    A_multiwake(:,:) = 0._pr
    do i = 1, grid%nelem  ! over all control points
        point(:) = ctrl_pts(:,i)
        do j=1, wake%nwp ! wake panels
            circ = 1._pr ! unitary circulation
            n = kuttaedges%edge2wpan(j) ! wake panel numeration (sequential)
            p=0; op=1 ! op=1 only for wake rings case
            jaux=j; staux=0
            fstlen=0; stret=0
            call VORING ( n,point,circ,vind_ring, ra, rb,p,jaux,staux,op,fstlen,stret) ! induced velocity at control point due to the j-panel wake
            call DOTPRODUCT(vind_ring(:),nor_vec(:,i),dot)
                
            do k = kuttaedges%esed2(j)+1, kuttaedges%esed2(j+1)  ! matrix entries
                trailpan = kuttaedges%esed1(k) ! body's panel sharing the separation edge
                          
                log1 = .false.; log2 = .false.
                do m=1, grid%nle
                    if ( n==internal_LE_detected(m) ) then ! determines if "n-wake" is an internal LE
                        log1 = .true.
                    else
                    end if
                end do
                do m=1, grid%plate_nle
                    if ( n==plate_LE_detected(m) ) then ! determines if "n-wake" is a plate's LE
                        log2 = .true.
                    else
                    end if
                end do   
                                   
                if ( kuttaedges%orien(k) == .true. ) then ! panel and wake have the same orientation
                    if (log1 == .true.) then ! for internal LEs
                        A_first(i,trailpan) = 0._pr
                    else if (log2 == .true.) then ! for plate's LEs
                        select case(input%detach_model)
                        case(1)
                            A_first(i,trailpan) = dot
                        case(2) ! for positive LE wake (Kutta condition)
                            A_first(i,trailpan) = -dot
                            !A_first(i,trailpan) = dot ! maintains the same sign during the assembly, but changes it during advection!
                        end select
                    else ! for trailing (and lateral) wakes
                        A_first(i,trailpan) = -dot 
                    end if
                else ! .false., for positive (different orientation) wakes respect its emiting panel
                    A_first(i,trailpan) = dot 
                end if
                    
            A_multiwake(i,trailpan) = A_multiwake(i,trailpan) + A_first(i,trailpan)
            end do  ! ends panels sharing the separation edge
        end do ! ends wake panels
        call DOTPRODUCT( -vdir(:), nor_vec(:,i), rhs(i) ) ! system of equations' RHS (only one wake row)
    end do  ! body panels
    
    rhs_freeflow(:) = rhs(:) ! for Kelvin's condition
    A_total(:,:) = A_body(:,:) + A_multiwake(:,:) ! A total (body + multiwake)
        
    A_solv(:,:) = A_total(:,:) ! to avoid undesired modification to original matrix after solving the system
    
    if (pr .eq. 4) then ! system solution
        call sgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! for single precision
    else
        call dgesv(grid%nelem,1,A_solv,grid%nelem,pivot,rhs,grid%nelem,ok) ! for double precision
    end if 
    bound_circ(:) = rhs(:) ! new panels' circulation
    
    !sum_circ_plate = sum(abs(bound_circ)) ! total plate's circulation
    
    step=0 ! for starting solution case (wake rings)
    call DETACHVORT ( step ) ! calculates the wake circulation strenght between the BVRs (internal and external wakes)
    !call START_FORCES_KJ( circ, mid_point, tot_area ) ! calls to starting (or steady) force calculation

End Subroutine STARTING_SOL

Subroutine START_FORCES_KJ ( circ, mid_point, tot_area ) ! hydrodynamic coefficients calculation through Kutta-Zhukovski (KJ) approach; JCPG/EOA
    Use ARRAYS, only : pr,input,grid,bound_circ,wake,vdir,pi,orien_reduc,internal_LE_detected,plate_LE_detected,gamma_wake_start
    Use MATH, only : CROSSPRODUCT
    Implicit none

    integer :: h,j,k,L,fst_nod,sec_nod,kaux,jaux,r,staux,op,m,fstlen,stret
    real(kind=pr) :: sum_mom,x_cm,den,cfx,cfy,cfz,cl,cd,cy,cn
    real(kind=pr) :: alpharad,betarad,sinalpha,cosalpha,sinbeta,cosbeta
    real(kind=pr) :: tot_force_steady(3),pan_force_steady(3),vind_body(3),vind_ring(3),ra(3),rb(3),vind_wakeonpoint(3),vind_total(3),del_gamma(3),segm_force(3),tot_force(3)
    real(kind=pr), intent(in) :: tot_area
    real(kind=pr), intent(out) :: mid_point(3),circ
    logical :: log1,log2

    tot_force_steady(:) = 0._pr
    sum_mom = 0._pr
    x_cm = input%x_pos_mom_fact*input%char_len ! moment axis position

    do k=1, grid%nelem ! over BVR/panels
        pan_force_steady(:) = 0._pr
        
        fst_nod = grid%panel(4,k) ! first node (it depends on mesh numeration's starting point)
        sec_nod = grid%panel(3,k) ! second node (it depends on mesh numeration's starting point)
                
        log1 = .false.; log2 = .false.
        do m=1, grid%nle
            if ( fst_nod==grid%panel(2,internal_LE_detected(m)) .and. sec_nod==grid%panel(1,internal_LE_detected(m))) then
                log1 = .true.
            else
            end if
        end do
        
        do m=1, grid%plate_nle
            if ( fst_nod==grid%panel(2,plate_LE_detected(m)) .and. sec_nod==grid%panel(1,plate_LE_detected(m)) ) then
                select case(input%detach_model)
                case(1)
                    log2 = .true.
                case(2)
                    log2 = .false.
                end select
            else
            end if
        end do 
        
        if (log1==.true. .or. log2==.true.) then
            vind_body(:) = 0._pr
            mid_point(:) = ((grid%coord(:,grid%panel(4,k)) + grid%coord(:,grid%panel(3,k)))/2._pr)
        
            do L=1, grid%nelem ! over panels
                circ = bound_circ(L)
                kaux=0; r=0; 
                staux=0; op=0 ! for bounded VR
                fstlen=0; stret=0
                call VORING (L,mid_point,circ,vind_ring, ra, rb,kaux,r,staux,op,fstlen,stret) ! bounded vortex ring's induced velocity
                vind_body(:) = vind_body(:) + vind_ring(:)
            end do !ends L
        
            vind_wakeonpoint(:) = 0._pr
            h=1
            do j = grid%nelem + 1, grid%nelem + wake%nwp ! wake panel number, for 4x4 case: 17->56,96,136...
                circ = gamma_wake_start(j - grid%nelem) ! for 4x4 case: 41->80,120...
                kaux=0; jaux=h; 
                staux=0; op=1 ! for wake rings
                fstlen=0; stret=0
                call VORING ( j,mid_point,circ,vind_ring, ra, rb,kaux, jaux ,staux,op,fstlen,stret) ! wake vortex ring's induced velocity
            
                if ( orien_reduc(h) == .true. ) then ! determines the correct orientation of the detached wakes
                    vind_ring(:) = -vind_ring(:)
                else
                end if

                vind_wakeonpoint(:) = vind_wakeonpoint(:) + vind_ring(:)
                h=h+1
                if (h==wake%nwp+1) then
                    h=1
                else
                end if
            end do ! ends j
            vind_total(:) = vdir(:) + vind_wakeonpoint(:) + vind_body(:) ! total local flow velocity
        
            ! Only LE for FMVLM-based scheme
            ra(:) = grid%coord(:,grid%panel(4,k)) ! bounded segment's first node position
            rb(:) = grid%coord(:,grid%panel(3,k)) ! bounded segment's second node position
         
            del_gamma(:) = (rb(:)-ra(:))*bound_circ(k) ! panel's remaining segments
            call CROSSPRODUCT(vind_total(:), del_gamma(:), segm_force(:)) ! forces on vortex segment
        else
            segm_force(:) = 0._pr ! for Kutta's type separation edge
        end if
        
        segm_force(:) = input%dens*segm_force(:)
        pan_force_steady(:) = pan_force_steady(:) + segm_force(:) ! perpendicular force on the k-panel
        tot_force_steady(:) = tot_force_steady(:) + pan_force_steady(:) ! body's steady force
    end do ! ends k
        
    tot_force(:) = tot_force_steady(:) ! body's total force (only steady)

    den = input%dens * input%q_inf * input%q_inf * tot_area ! denominator
    cfx = 2._pr*tot_force(1)/den; cfy = 2._pr*tot_force(2)/den; cfz = 2._pr*tot_force(3)/den  ! body axes coefficients
    alpharad = input%alpha*(pi/180._pr); betarad  = 0._pr !input%beta*(pi/180._pr)  ! angles in radians

    sinalpha = sin(alpharad); cosalpha = cos(alpharad); sinbeta = sin(betarad); cosbeta = cos(betarad)

    cl  = cfz*cosalpha - cfx*sinalpha                                ! lift coef.            (CL)
    cd = cfx*cosalpha*cosbeta + cfy*sinbeta + cfz*sinalpha*cosbeta   ! drag coef.            (CD)
    cy  = -cfx*cosalpha*sinbeta + cfy*cosbeta + cfz*sinalpha*sinbeta ! lateral force coef.   (CY)
    cn  = cl*cosalpha + cd*sinalpha                                  ! normal force coef.    (CN)
    !cm = sum_mom/(0.5_pr*rho*q_inf*q_inf*(env*char_len)*char_len)   ! pitching moment coef. (CM)

    400 format(f12.4)
    print 400, cl; print 400, cd; !print 400, cm; !print*, cn
    print *,'---'

    500 format(f12.4,f12.4,f12.4)
    open(4, file='aero_coef.dat', access='append') ! writes output file
    write(4,500) cl, cd !, cm
    close(4)

End Subroutine START_FORCES_KJ

Subroutine VORING ( n,point,circ,vind_ring, ra, rb, k,j,step,op,fstlen,stret) ! calculates induced velocity by a vortex ring; JCPG
    Use ARRAYS, only : pr,grid
    Implicit none
    real(kind=pr), intent(in) :: circ
    real(kind=pr), intent(out) :: vind_ring(3),ra(3),rb(3)
    real(kind=pr), intent(inout) :: point(3)
    real(kind=pr) :: v_ind(3)
    integer, intent(in) :: n
    integer, intent(inout) :: k,step,j,op,fstlen,stret
    integer :: nodes(5)
    
    vind_ring(:) = 0._pr
    if (step>=1 .and. op==1) then ! for wake rings (for unsteady case)
        nodes(:) = [ grid%panel_wake(1:grid%elemtype,n) , grid%panel_wake(1,n) ]
    else ! for bounded and wake rings (for starting solution)
        nodes(:) = [ grid%panel(1:grid%elemtype,n) , grid%panel(1,n) ]
    end if

    do k=1, grid%elemtype ! over the 3 or 4 vortex segments
    
        if (step==0 .and. op==1) then ! for starting (or steady) solution (only for wakes)
            ra = grid%coord_start(:,nodes(k+1)) ! clockwise direction (KATZ, 2001)  
            rb = grid%coord_start(:,nodes(k)) 
        else if (step>=1 .and. op==1) then ! for unsteady case (only for wakes)
            ra = grid%coord_wake(:,nodes(k+1)) ! clockwise direction 
            rb = grid%coord_wake(:,nodes(k)) 
        else ! for bounded rings
            ra = grid%coord(:,nodes(k+1)) ! clockwise direction 
            rb = grid%coord(:,nodes(k)) 
        end if
    
        call VATISTAS_CM(ra,rb,point,circ,v_ind,k,j,step,op,fstlen,stret) ! single vortex filament induced velocity 
        vind_ring(:) = vind_ring(:) + v_ind(:) ! vortex ring's induced velocity
    end do
  
End Subroutine VORING

Subroutine VATISTAS_CM(ra,rb,point,circ,v_ind,k,j,step,op,fstlen,stret) ! calculates regularized induced velocity by a single vortex ring segment; JCPG
    Use ARRAYS, only : pr,fourpi,input,pi,fstwake_len,inv_4pi,fouroverthree
    Use MATH, only : DOTPRODUCT,VECTORNORM,CROSSPRODUCT
    Implicit none 
    integer, intent(in) :: k,j,step,op,fstlen,stret
    
    real(kind=pr), intent(in) :: ra(3),rb(3),point(3),circ
    real(kind=pr), intent(out) :: v_ind(3)
    real(kind=pr) :: r0(3),r1(3),r2(3),r1xr2(3),c_v_par,c_v,s_v,a_v,xxx,exp_int,q,fac,rel_stretch,eff_rad
    real(kind=pr) :: norm_r1,norm_r2,dot_r0_r1,dot_r1_r2,norm_r1xr2,dot_r0_r2,norm_r0,vcr3,vcr_cyl
    
    r0 = rb - ra ! inductor segment's local vector
    r1 = point(:) - ra ! local vector between the inductor vector's first node and the evaluation point
    r2 = point(:) - rb ! local vector between the inductor vector's second node and the evaluation point
    
    call DOTPRODUCT(r0,r1,dot_r0_r1)     ! scalar product between r0 and r1
    call DOTPRODUCT(r0,r2,dot_r0_r2)     ! scalar product between r0 and r2
    call DOTPRODUCT(r1,r2,dot_r1_r2)     ! scalar product between r1 and r2
    call VECTORNORM(r0,norm_r0)          ! vector norm of r0
    call VECTORNORM(r1,norm_r1)          ! vector norm of r1
    call VECTORNORM(r2,norm_r2)          ! vector norm of r2
    call CROSSPRODUCT(r1,r2,r1xr2)       ! cross product between r1 and r2
    call VECTORNORM(r1xr2,norm_r1xr2)    ! vector norm of r1xr2
    
    vcr3 = input%core_rad_init*input%core_rad_init*input%core_rad_init
    vcr_cyl = sqrt((fouroverthree*vcr3)/norm_r0) ! equivalent vortex core radius (for vortex filament/tube)
    
    if (op==1 .and. fstlen==1 .and. vcr_cyl/=0._pr) then ! vortex stretching (wake rings)
    
        fstwake_len(k,j)=norm_r0 ! determines the segment's lenght of the first wake rings row; k-segment, j-wake ring
    else; end if
    
    if (vcr_cyl==0._pr) then
        !*** Cut-off model
        if (norm_r1<=input%tol_rad .or. norm_r2<=input%tol_rad .or. (norm_r1xr2*norm_r1xr2)<=input%tol_rad) then 
            v_ind(:) = 0._pr ! to avoid overshooted velocity
        else ! for distances longer than tolerance  
            fac = norm_r1xr2*norm_r1xr2 ; fac = inv_4pi / fac
            q = (dot_r0_r1/norm_r1 - dot_r0_r2/norm_r2)*circ*fac 
            v_ind(:) = q*r1xr2(:) ! induced velocity vector
        end if
        !***
    else ! with a vortex core radius (vcr)
        !*** Vatistas' Core Model (LEUTHOLD, 2015)
        s_v = 0.5_pr*( norm_r1 + norm_r2 + norm_r0 )
        a_v = (2._pr*sqrt( s_v*(s_v-norm_r1)*(s_v-norm_r2)*(s_v-norm_r0) )) / (norm_r0) ! perpendicular distance between the segment and the evaluation point
             
        if (a_v>=input%tol_rad) then ! for distances longer than tolerance
            
            !!! Kauffman/Scully (better fit than Lamb-Osen acording to (HOMMES, 2015))
            !c_v_par = ( ((norm_r0*norm_r1)*(norm_r0*norm_r1)) - (dot_r0_r1*dot_r0_r1) ) / (norm_r0*norm_r0)  ! main parenthesis of c_v
            !
            !if (input%vx_stretch==1 .and. op==1 .and. stret==1 .and. step>=2) then ! for vortex stretching induced velocity correction
            !    rel_stretch = norm_r0/fstwake_len(k,j) ! -56, -112
            !!    eff_rad = input%core_rad_init*rel_stretch**-0.5_pr ! effective radius 
            !    eff_rad = vcr_cyl*rel_stretch**-0.5_pr ! effective radius 
            
            !    xxx = eff_rad**2._pr + c_v_par ! for Kaufmman/Scully VCM
            !else ! w/o vortex stretching, for bounded VR and all iterations
            !!    xxx =  input%core_rad_init**2._pr + c_v_par
            !    xxx =  vcr_cyl**2._pr + c_v_par                     
            !end if
            !   
            !exp_int = 1._pr/xxx ! only for Kauffman/Scully VCM
            !c_v = c_v_par*exp_int ! effective viscosity parameter
            !!!
                
            !! Lamb-Osen with vortex stretching (SEBASTIAN, 2012)
            c_v_par = ( ((norm_r0*norm_r1)*(norm_r0*norm_r1)) - (dot_r0_r1*dot_r0_r1) ) / (norm_r0*norm_r0)  ! main parenthesis of c_v
            
            if (input%vx_stretch==1 .and. op==1 .and. step>=2) then ! for vortex stretching induced velocity correction
                rel_stretch = norm_r0/fstwake_len(k,j) ! -56, -112
                !eff_rad = input%core_rad_init*rel_stretch**-0.5_pr ! effective radius 
                eff_rad = vcr_cyl*rel_stretch**-0.5_pr ! effective radius 
                xxx = eff_rad**4._pr + c_v_par**2._pr ! for Lamb-Osen VCM
            else ! w/o vortex stretching, for bounded VR and all iterations
                !xxx = input%core_rad_init**4._pr + c_v_par**2._pr      
                xxx = vcr_cyl**4._pr + c_v_par**2._pr                      
            end if
                   
            exp_int = xxx**-0.5_pr   
            c_v = c_v_par*exp_int ! effective viscosity parameter
            !!
            v_ind(:) = ( c_v*circ*(norm_r1 + norm_r2)*(r1xr2) ) / ( fourpi*norm_r1*norm_r2 * (norm_r1*norm_r2 + dot_r1_r2) ) ! induced velocity vector
        else
            v_ind(:) = 0._pr
        end if
        !***
    end if
         
    End Subroutine VATISTAS_CM