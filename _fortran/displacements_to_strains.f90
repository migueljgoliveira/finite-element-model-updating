! -----------------------------------------------------------------------------
!
!   Table of Contents
!
!       . displacements_to_strains()
!
!       . update_lagrangian_hughes_winget()
!       . update_lagrangian()
!       . total_lagrangian()
!       
!       . element_quad4_updated_lagrangian()
!       . element_quad4_total_lagrangian()
!
!       . gauss_jordan()
!
!       . polar_decomposition()
!       . spectral_decomposition()
!       . jacobi()
!
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       displacements_to_strains()
!
! -----------------------------------------------------------------------------
!
!   Description
!
!       . The displacement field is converted in strains through the 
!    deformation gradient. In addition, the element area for each element
!    is computed.
!       . To calculate the deformation gradient the shape functions of the 
!   chosen element are used with the local coordinates of the centre of the
!   the element [0,0].
!       . The displacement field is given always relative to the undeformed 
!   configuration.
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Computatation of the deformation gradient
!       . Polar Decomposition
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           nelems      : number of elements
!           nodpelem    : number of nodes per element
!           incmax      : maximum experimental increments
!           nnodes      : number of nodes
!           elems       : element connectivity
!           displ       : displacement field
!           nodes       : nodes coordinates
!
!       . Output:
!           dfgrd       : deformation gradient
!           indfgrd     : inverse of deformation gradient
!           rot         : rotation tensor
!           strain      : strain tensor
!           detdfgrd    : determinant of deformation gradient
!           area        : area of elements
!
!       . Auxiliar:
!           ldispl      : displacement of each node of element
!           dudx        : derivatives of displacement field relative to 
!                         undeformed configuration
!           Ident       : identity matrix
!           grot        : global rotations
!           incrot      : incremental rotations
!
! ----------------------------------------------------------------------

    subroutine displacements_to_strains(nelems,nodpelem,incmax,nnodes,&
                                        elems,nodes,displ,flag,&
                                        strain,dfgrd,indfgrd,detdfgrd,rot,area)
                                        
        implicit none

        ! Auxiliar
        integer(4),parameter                   :: ndof=2,ndi=2,nshr=1,&
                                                  ntens=3,ncomp=4
        real(8),dimension(ndof,ndof)           :: ident
        logical                                :: grot,incrot

        ! Input
        integer(4)                             :: nelems,nodpelem,incmax,nnodes
        integer(4),dimension(nelems,nodpelem)  :: elems
        real(8),dimension(nnodes,ndof)         :: nodes
        real(8),dimension(incmax,nnodes,ndof)  :: displ
        logical                                :: flag

        ! Output
        real(8),dimension(nelems)              :: area
        real(8),dimension(incmax,nelems)       :: detdfgrd
        real(8),dimension(incmax,nelems,ncomp) :: dfgrd,rot,indfgrd
        real(8),dimension(incmax,nelems,ntens) :: strain

        ! --------------------------------------------------------------

!f2py intent(in) nelems,nodpelem,incmax,nnodes,elems,displ,nodes,flag
!f2py intent(out) strain,dfgrd,indfgrd,detdfgrd,rot,area
!f2py depend(nelems,nodpelem) elems
!f2py depend(nnodes,ndof) nodes
!f2py depend(incmax,nnodes,ndof) displ
!f2py depend(nelems) area
!f2py depend(incmax,nelems,ncomp) dfgrd,rot,indfgrd
!f2py depend(incmax,nelems,ntens) strain
!f2py depend(incmax,nelems) detdfgrd

        ! --------------------------------------------------------------

        ident = 0.0d0
        ident(1,1) = 1.0d0
        ident(2,2) = 1.0d0

        grot = .False.
        incrot = .False.

        if ( incrot .eqv. .True. ) then
            call update_lagrangian_hughes_winget(incmax,nelems,nodpelem,&
                                                 displ,ident,elems,nnodes,&
                                                 nodes,flag,dfgrd,rot,indfgrd,&
                                                 strain,detdfgrd,area)
        else if ( grot .eqv. .True. ) then
            call update_lagrangian(incmax,nelems,nodpelem,displ,ident,elems,&
                                   nnodes,nodes,flag,dfgrd,rot,indfgrd,strain,&
                                   detdfgrd,area)
        else
            call total_lagrangian(incmax,nelems,nodpelem,displ,ident,elems,&
                                  nnodes,nodes,flag,dfgrd,rot,indfgrd,strain,&
                                  detdfgrd,area)
        end if
    end subroutine

! ----------------------------------------------------------------------
!
!   Subroutine
!
!       update_lagrangian_hughes_winget()
!
! ----------------------------------------------------------------------

    subroutine update_lagrangian_hughes_winget(incmax,nelems,nodpelem,displ,&
                                               ident,elems,nnodes,nodes,flag,&
                                               dfgrd,rot,indfgrd,strain,&
                                               detdfgrd,area)

        implicit none

        ! Auxiliar
        integer(4)                             :: ierr,i,j
        integer(4),parameter                   :: ndof=2,ndi=2,nshr=1,ntens=3,&
                                                  ncomp=4
        real(8)                                :: auxarea,fdeterm2
        real(8),dimension(ndi)                 :: auxB
        real(8),dimension(ndi,ndi)             :: auxdfgrd,mrot,auxstrain,&
                                                  auxstrainmid,auxdfgrdt,&
                                                  auxdfgrdf,mrott,mrotf
        real(8),dimension(ndof,ndof)           :: dudx
        real(8),dimension(nodpelem,ndof)       :: ldispl,ldisplu

        ! Input
        integer(4)                             :: incmax,nelems,nodpelem,nnodes
        integer(4),dimension(nelems,nodpelem)  :: elems
        real(8),dimension(ndof,ndof)           :: ident
        real(8),dimension(nnodes,ndof)         :: nodes
        real(8),dimension(incmax,nnodes,ndof)  :: displ
        logical                                :: flag

        ! Output
        real(8),dimension(nelems)              :: area
        real(8),dimension(incmax,nelems)       :: detdfgrd
        real(8),dimension(incmax,nelems,ncomp) :: dfgrd,rot,indfgrd
        real(8),dimension(incmax,nelems,ntens) :: strain

        ! --------------------------------------------------------------

        do i = 1,incmax
            do j = 1,nelems
                if ( i == 1 ) then
                    ldispl(:,:) = displ(i,elems(j,:),:)
                    ldisplu(:,:) = displ(i,elems(j,:),:)
                else
                    ldispl(:,:) = displ(i-1,elems(j,:),:)
                    ldisplu(:,:) = displ(i,elems(j,:),:) - &
                                   displ(i-1,elems(j,:),:)
                end if

                call element_quad4_updated_lagrangian(j,ndof,ndi,nnodes,&
                                                      nelems,nodpelem,nodes,&
                                                      elems,ldispl,ldisplu,&
                                                      dudx,auxarea)
                
                if ( i == 1 ) area(j) = auxarea

                auxdfgrd = dudx + ident
                call polar_decomposition(auxdfgrd,mrot,auxstrain)
                auxstrain = 0.5d0*(dudx+transpose(dudx))

                if ( i == 1 ) then
                    ! STRAIN TENSOR
                    strain(i,j,1) = auxstrain(1,1)
                    strain(i,j,2) = auxstrain(2,2)
                    strain(i,j,3) = auxstrain(1,2)*2

                    if ( flag ) then
                        ! DEFORMATION GRADIENT
                        dfgrd(i,j,1) = auxdfgrd(1,1)
                        dfgrd(i,j,2) = auxdfgrd(1,2)
                        dfgrd(i,j,3) = auxdfgrd(2,1)
                        dfgrd(i,j,4) = auxdfgrd(2,2)

                        ! ROTATION TENSOR
                        rot(i,j,1) = mrot(1,1)
                        rot(i,j,2) = mrot(1,2)
                        rot(i,j,3) = mrot(2,1)
                        rot(i,j,4) = mrot(2,2)

                        ! DETERMINANT OF DEFORMATION GRADIENT
                        detdfgrd(i,j) = fdeterm2(auxdfgrd)

                        ! INVERSE DEFORMATION GRADIENT
                        call gauss_jordan(auxdfgrd,ndi,auxB,ierr)

                        indfgrd(i,j,1) = auxdfgrd(1,1)
                        indfgrd(i,j,2) = auxdfgrd(1,2)
                        indfgrd(i,j,3) = auxdfgrd(2,1)
                        indfgrd(i,j,4) = auxdfgrd(2,2)
                    end if
                else
                    ! STRAIN TENSOR
                    strain(i,j,1) = strain(i-1,j,1) + auxstrain(1,1)
                    strain(i,j,2) = strain(i-1,j,2) + auxstrain(2,2)
                    strain(i,j,3) = strain(i-1,j,3) + auxstrain(1,2)*2

                    if ( flag ) then
                        ! DEFORMATION GRADIENT
                        auxdfgrdt(1,1) = dfgrd(i-1,j,1)
                        auxdfgrdt(1,2) = dfgrd(i-1,j,2)
                        auxdfgrdt(2,1) = dfgrd(i-1,j,3)
                        auxdfgrdt(2,2) = dfgrd(i-1,j,4)

                        auxdfgrdf = matmul(auxdfgrd,auxdfgrdt)

                        dfgrd(i,j,1) = auxdfgrdf(1,1)
                        dfgrd(i,j,2) = auxdfgrdf(1,2)
                        dfgrd(i,j,3) = auxdfgrdf(2,1)
                        dfgrd(i,j,4) = auxdfgrdf(2,2)

                        
                        ! ROTATION TENSOR
                        mrott(1,1) = rot(i-1,j,1)
                        mrott(1,2) = rot(i-1,j,2)
                        mrott(2,1) = rot(i-1,j,3)
                        mrott(2,2) = rot(i-1,j,4)

                        Mrotf = matmul(Mrot,Mrott)

                        rot(i,j,1) = mrotf(1,1)
                        rot(i,j,2) = mrotf(1,2)
                        rot(i,j,3) = mrotf(2,1)
                        rot(i,j,4) = mrotf(2,2)

                        ! DETERMINANT OF DEFORMATION GRADIENT
                        detdfgrd(i,j) = fdeterm2(auxdfgrdf)

                        ! INVERSE OF DEFORMATION GRADIENT
                        call gauss_jordan(auxdfgrdf,ndi,auxB,ierr)

                        indfgrd(i,j,1) = auxdfgrdf(1,1)
                        indfgrd(i,j,2) = auxdfgrdf(1,2)
                        indfgrd(i,j,3) = auxdfgrdf(2,1)
                        indfgrd(i,j,4) = auxdfgrdf(2,2)
                    end if
                end if
            end do
        end do

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       update_lagrangian()
!
! -----------------------------------------------------------------------------

    subroutine update_lagrangian(incmax,nelems,nodpelem,displ,ident,elems,&
                                 nnodes,nodes,flag,dfgrd,rot,indfgrd,strain,&
                                 detdfgrd,area)

        implicit none
        
        ! Auxiliar
        integer(4)                             :: ierr,i,j
        integer(4),parameter                   :: ndof=2,ndi=2,nshr=1,ntens=3,&
                                                  ncomp=4
        real(8)                                :: auxarea,fdeterm2
        real(8),dimension(ndi)                 :: auxB
        real(8),dimension(ndi,ndi)             :: auxdfgrd,mrot,auxstrain,&
                                                  auxstrainmid,auxdfgrdt,&
                                                  auxdfgrdf,mrott,mrotf
        real(8),dimension(ndof,ndof)           :: dudx
        real(8),dimension(nodpelem,ndof)       :: ldispl,ldisplu

        ! Input
        integer(4)                             :: incmax,nelems,nodpelem,nnodes
        integer(4),dimension(nelems,nodpelem)  :: elems
        real(8),dimension(ndof,ndof)           :: ident
        real(8),dimension(nnodes,ndof)         :: nodes
        real(8),dimension(incmax,nnodes,ndof)  :: displ
        logical                                :: flag

        ! Output
        real(8),dimension(nelems)              :: area
        real(8),dimension(incmax,nelems)       :: detdfgrd
        real(8),dimension(incmax,nelems,ncomp) :: dfgrd,rot,indfgrd
        real(8),dimension(incmax,nelems,ntens) :: strain

        ! ---------------------------------------------------------------------

        do i = 1,incmax
            do j = 1,nelems
                if ( i == 1 ) then
                    ldispl(:,:) = 0.0d0
                    ldisplu(:,:) = displ(i,elems(j,:),:)
                else
                    ldispl(:,:) = displ(i-1,elems(j,:),:)
                    ldisplu(:,:) = displ(i,elems(j,:),:) - &
                                   displ(i-1,elems(j,:),:)
                end if

                call element_quad4_updated_lagrangian(j,ndof,ndi,nnodes,&
                                                      nelems,nodpelem,nodes,&
                                                      elems,ldispl,ldisplu,dudx,&
                                                      auxarea)

                if ( i == 1 ) area(j) = auxarea

                auxdfgrd = dudx + ident
                call polar_decomposition(auxdfgrd,mrot,auxstrain)

                if ( i == 1 ) then
                    ! STRAIN TENSOR
                    auxstrainmid = matmul(transpose(mrot),auxstrain)
                    auxstrain = matmul(auxstrainmid,mrot)

                    strain(i,j,1) = auxstrain(1,1)
                    strain(i,j,2) = auxstrain(2,2)
                    strain(i,j,3) = auxstrain(1,2)*2

                    if ( flag ) then
                        ! DEFORMATION GRADIENT
                        dfgrd(i,j,1) = auxdfgrd(1,1)
                        dfgrd(i,j,2) = auxdfgrd(1,2)
                        dfgrd(i,j,3) = auxdfgrd(2,1)
                        dfgrd(i,j,4) = auxdfgrd(2,2)

                        ! ROTATION TENSOR
                        rot(i,j,1) = mrot(1,1)
                        rot(i,j,2) = mrot(1,2)
                        rot(i,j,3) = mrot(2,1)
                        rot(i,j,4) = mrot(2,2)

                        ! DETERMINANT OF DEFORMATION GRADIENT
                        detdfgrd(i,j) = fdeterm2(auxdfgrd)

                        ! INVERSE DEFORMATION GRADIENT
                        call gauss_jordan(auxdfgrd,ndi,auxB,ierr)

                        indfgrd(i,j,1) = auxdfgrd(1,1)
                        indfgrd(i,j,2) = auxdfgrd(1,2)
                        indfgrd(i,j,3) = auxdfgrd(2,1)
                        indfgrd(i,j,4) = auxdfgrd(2,2)
                    end if
                else
                    ! STRAIN TENSOR
                    mrott(1,1) = rot(i-1,j,1)
                    mrott(1,2) = rot(i-1,j,2)
                    mrott(2,1) = rot(i-1,j,3)
                    mrott(2,2) = rot(i-1,j,4)

                    mrotf = matmul(mrot,mrott)
                    auxstrainmid = matmul(transpose(mrotf),auxstrain)
                    auxstrain = matmul(auxstrainmid,mrotf)

                    strain(i,j,1) = strain(i-1,j,1) + auxstrain(1,1)
                    strain(i,j,2) = strain(i-1,j,2) + auxstrain(2,2)
                    strain(i,j,3) = strain(i-1,j,3) + auxstrain(1,2)*2

                    if ( flag ) then
                        ! DEFORMATION GRADIENT
                        auxdfgrdt(1,1) = dfgrd(i-1,j,1)
                        auxdfgrdt(1,2) = dfgrd(i-1,j,2)
                        auxdfgrdt(2,1) = dfgrd(i-1,j,3)
                        auxdfgrdt(2,2) = dfgrd(i-1,j,4)
                        
                        auxdfgrdf = matmul(auxdfgrd,auxdfgrdt)

                        dfgrd(i,j,1) = auxdfgrdf(1,1)
                        dfgrd(i,j,2) = auxdfgrdf(1,2)
                        dfgrd(i,j,3) = auxdfgrdf(2,1)
                        dfgrd(i,j,4) = auxdfgrdf(2,2)

                        ! ROTATION TENSOR
                        mrotf = matmul(mrot,mrott)
                        rot(i,j,1) = mrotf(1,1)
                        rot(i,j,2) = mrotf(1,2)
                        rot(i,j,3) = mrotf(2,1)
                        rot(i,j,4) = mrotf(2,2)

                        ! DETERMINANT OF DEFORMATION GRADIENT
                        detdfgrd(i,j) = fdeterm2(auxdfgrdf)

                        ! INVERSE DEFORMATION GRADIENT
                        call gauss_jordan(auxdfgrdf,ndi,auxB,ierr)

                        indfgrd(i,j,1) = auxdfgrdf(1,1)
                        indfgrd(i,j,2) = auxdfgrdf(1,2)
                        indfgrd(i,j,3) = auxdfgrdf(2,1)
                        indfgrd(i,j,4) = auxdfgrdf(2,2)
                    end if
                end if
            end do
        end do

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       total_lagrangian()
!
! -----------------------------------------------------------------------------

    subroutine total_lagrangian(incmax,nelems,nodpelem,displ,ident,elems,&
                                nnodes,nodes,flag,dfgrd,rot,indfgrd,strain,&
                                detdfgrd,area)

        implicit none

        ! Auxiliar
        integer(4)                             :: ierr,i,j
        integer(4),parameter                   :: ndof=2,ndi=2,nshr=1,ntens=3,&
                                                  ncomp=4
        real(8)                                :: auxarea,fdeterm2
        real(8),dimension(ndi)                 :: auxB
        real(8),dimension(ndi,ndi)             :: auxdfgrd,mrot,auxstrain,&
                                                  auxstrainmid
        real(8),dimension(ndof,ndof)           :: dudx
        real(8),dimension(nodpelem,ndof)       :: ldispl

        ! Input
        integer(4)                             :: incmax,nelems,nodpelem,nnodes
        integer(4),dimension(nelems,nodpelem)  :: elems
        real(8),dimension(ndof,ndof)           :: ident
        real(8),dimension(nnodes,ndof)         :: nodes
        real(8),dimension(incmax,nnodes,ndof)  :: displ
        logical                                :: flag

        ! Output
        real(8),dimension(nelems)              :: area
        real(8),dimension(incmax,nelems)       :: detdfgrd
        real(8),dimension(incmax,nelems,ncomp) :: dfgr                    d,rot,indfgrd
        real(8),dimension(incmax,nelems,ntens) :: strain

        ! ---------------------------------------------------------------------

        do i = 1,incmax
            do j = 1,nelems
                ldispl(1,:) = displ(i,elems(j,1),:)
                ldispl(2,:) = displ(i,elems(j,2),:)
                ldispl(3,:) = displ(i,elems(j,3),:)
                ldispl(4,:) = displ(i,elems(j,4),:)

                call element_quad4_total_lagrangian(j,nnodes,nelems,nodpelem,&
                                                    nodes,elems,ldispl,dudx,&
                                                    auxarea)

                if ( i == 1 ) area(j) = auxarea
                
                auxdfgrd = dudx + ident
                call polar_decomposition(auxdfgrd,mrot,auxstrain)

                ! STRAIN TENSOR
                auxstrainmid = matmul(transpose(mrot),auxstrain)
                auxstrain = matmul(auxstrainmid,mrot)

                strain(i,j,1) = auxstrain(1,1)
                strain(i,j,2) = auxstrain(2,2)
                strain(i,j,3) = auxstrain(1,2)*2

                if ( flag ) then
                    ! DEFORMATION GRADIENT
                    dfgrd(i,j,1) = auxdfgrd(1,1)
                    dfgrd(i,j,2) = auxdfgrd(1,2)
                    dfgrd(i,j,3) = auxdfgrd(2,1)
                    dfgrd(i,j,4) = auxdfgrd(2,2)

                    ! ROTATION TENSOR
                    rot(i,j,1) = mrot(1,1)
                    rot(i,j,2) = mrot(1,2)
                    rot(i,j,3) = mrot(2,1)
                    rot(i,j,4) = mrot(2,2)

                    ! DETERMINANT OF DEFORMATION GRADIENT
                    detdfgrd(i,j) = fdeterm2(auxdfgrd)

                    ! INVERSE DEFORMATION GRADIENT
                    call gauss_jordan(auxdfgrd,ndi,auxB,ierr)

                    indfgrd(i,j,1) = auxdfgrd(1,1)
                    indfgrd(i,j,2) = auxdfgrd(1,2)
                    indfgrd(i,j,3) = auxdfgrd(2,1)
                    indfgrd(i,j,4) = auxdfgrd(2,2)
                end if
            end do
        end do

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       element_quad4_total_lagrangian()
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Compute the gradient of the displacement field relative to the
!       undeformed configuration for a total Lagrangian formulation
!       . Compute the element's area
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           nel         : element number
!           nnodes      : number of nodes
!           nelems      : number of elements
!           nodpelem    : number of nodes per element
!           nodes       : coordinates of nodes
!           elems       : element connectivity
!           ldispl      : displacement of each node of element
!
!       . Output:
!           dudx        : gradient of displacement field
!           area        : element's area
!
!       . Auxiliar:
!           jac         : jacobian matrix
!           jacinv      : inverse of jacobian matrix
!           duden       : 
!           dNde        : 
!           dNdn        :
!
! -----------------------------------------------------------------------------

    subroutine element_quad4_total_lagrangian(nel,nnodes,nelems,nodpelem,&
                                              nodes,elems,ldispl,dudx,area)

        implicit none

        ! Auxiliar
        integer(4)                            :: i,ierr
        integer(4),parameter                  :: inodes=4,ndof=2,ndi=2
        real(8)                               :: xi,eta,fdeterm2
        real(8),dimension(inodes)             :: dNde,dNdn,a1,a2
        real(8),dimension(ndof)               :: temp2
        real(8),dimension(ndof,ndof)          :: jac,jacinv,duden
        
        ! Input
        integer(4)                            :: nel,nelems,nnodes,nodpelem
        integer(4),dimension(nelems,nodpelem) :: elems
        real(8),dimension(inodes,ndof)        :: ldispl
        real(8),dimension(nnodes,ndof)        :: nodes

        ! Output
        real(8)                               :: area
        real(8),dimension(ndi,ndi)            :: dudx

        !----------------------------------------------------------------------

        xi  = 0.0d0
        eta = 0.0d0
        a1 = (/ -1.0d0, -1.0d0,  1.0d0, 1.0d0 /)
        a2 = (/  1.0d0, -1.0d0, -1.0d0, 1.0d0 /)

        ! DERIVATIVES OF SHAPE FUNCTIONS (NATURAL COORDINATES)
        dNde = 0.0d0
        dNdn = 0.0d0
        do i = 1,inodes
            dNde(i) = a1(i)*(1+a2(i)*eta)/4.0d0
            dNdn(i) = a2(i)*(1+a1(i)*xi )/4.0d0
        end do

        ! JACOBIAN MATRIX AND ITS INVERSE
        jac = 0.0d0
        do i = 1,inodes
            jac(1,1) = jac(1,1) + dNde(i)*nodes(elems(nel,i),1)
            jac(2,1) = jac(2,1) + dNde(i)*nodes(elems(nel,i),2)
            jac(1,2) = jac(1,2) + dNdn(i)*nodes(elems(nel,i),1)
            jac(2,2) = jac(2,2) + dNdn(i)*nodes(elems(nel,i),2)
        end do
        area = 4.0d0*fdeterm2(jac)
        temp2 = 0.0d0
        jacinv = 0.0d0
        jacinv = jac
        call gauss_jordan(jacinv,ndof,temp2,ierr)

        ! TRANSPOSE OF DISPLACEMENT GRADIENT
        dudx = 0.0d0
        duden = 0.0d0
        do i = 1,inodes
            duden(1,1) = duden(1,1) + dNde(i)*ldispl(i,1)
            duden(2,1) = duden(2,1) + dNde(i)*ldispl(i,2)
            duden(1,2) = duden(1,2) + dNdn(i)*ldispl(i,1)
            duden(2,2) = duden(2,2) + dNdn(i)*ldispl(i,2)
        end do
        dudx = 0.0d0
        dudx = matmul(duden,jacinv)

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       element_quad4_updated_lagrangian()
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Compute the gradient of the displacement field relative to the
!       undeformed configuration for an updated Lagrangian formulation
!       . Compute the element's area
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           nel         : element number
!           nnodes      : number of nodes
!           nelems      : number of elements
!           nodpelem    : number of nodes per element
!           nodes       : coordinates of nodes
!           elems       : element connectivity
!           ldispl      : displacement of each node of element
!           ldisplu     : displacement increment of each node of element
!
!       . Output:
!           dudx        : gradient of displacement field
!           area        : element's area
!
!       . Auxiliar:
!           jac         : jacobian matrix
!           jacinv      : inverse of jacobian matrix
!           duden       : 
!           dNde        : 
!           dNdn        :
!
! -----------------------------------------------------------------------------

    subroutine element_quad4_updated_lagrangian(nel,nnodes,nelems,nodpelem,&
                                                nodes,elems,ldispl,ldisplu,&
                                                dudx,area)

        implicit none

        ! Auxiliar
        integer(4),parameter                  :: inodes=4,ndof=2,ndi=2
        integer(4)                            :: i,ierr
        real(8)                               :: xi,eta,fdeterm2
        real(8),dimension(inodes)             :: dNde,dNdn,a1,a2
        real(8),dimension(ndof)               :: temp2
        real(8),dimension(ndof,ndof)          :: jac,jacinv,duden
        
        ! Input
        integer(4)                            :: nel,nnodes,nodpelem,nelems
        integer(4),dimension(nelems,nodpelem) :: elems
        real(8),dimension(inodes,ndof)        :: ldispl,ldisplu
        real(8),dimension(nnodes,ndof)        :: nodes

        ! Output
        real(8),dimension(ndi,ndi)            :: dudx
        real(8)                               :: area

        !----------------------------------------------------------------------

        xi  = 0.0d0
        eta = 0.0d0
        a1 = (/ -1.0d0, -1.0d0,  1.0d0, 1.0d0 /)
        a2 = (/  1.0d0, -1.0d0, -1.0d0, 1.0d0 /)

        ! DERIVATIVES OF SHAPE FUNCTIONS (NATURAL COORDINATES)
        dNde = 0.0d0
        dNdn = 0.0d0
        do i = 1,inodes
            dNde(i) = a1(i)*(1+a2(i)*eta)/4.0d0
            dNdn(i) = a2(i)*(1+a1(i)*xi )/4.0d0
        end do

        ! JACOBIAN MATRIX AND ITS INVERSE
        jac = 0.0d0
        do i = 1,inodes
            jac(1,1) = jac(1,1) + dNde(i)*(nodes(elems(nel,i),1)+ldispl(i,1))
            jac(2,1) = jac(2,1) + dNde(i)*(nodes(elems(nel,i),2)+ldispl(i,2))
            jac(1,2) = jac(1,2) + dNdn(i)*(nodes(elems(nel,i),1)+ldispl(i,1))
            jac(2,2) = jac(2,2) + dNdn(i)*(nodes(elems(nel,i),2)+ldispl(i,2))
        end do
        area = 4.0d0*fdeterm2(jac)
        temp2 = 0.0d0
        jacinv = 0.0d0
        jacinv = jac
        call gauss_jordan(jacinv,ndof,temp2,ierr)

        ! TRANSPOSE OF DISPLACEMENT GRADIENT
        dudx = 0.0d0
        duden = 0.0d0
        do i = 1,inodes
            duden(1,1) = duden(1,1) + dNde(i)*ldisplu(i,1)
            duden(2,1) = duden(2,1) + dNde(i)*ldisplu(i,2)
            duden(1,2) = duden(1,2) + dNdn(i)*ldisplu(i,1)
            duden(2,2) = duden(2,2) + dNdn(i)*ldisplu(i,2)
        end do
        dudx = 0.0d0
        dudx = matmul(duden,jacinv)

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       gauss_jordan()
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Compute the inverse of a matrix
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           n           : number of dimensions
!           F           : deformation gradient
!
!       . Output:
!           Mrot        : rotation tensor
!           auxStrain   : logarithmic strain tensor
!
!       . Auxiliar:
!           U           : right Cauchy strech tensor
!           V           : left Cauchy strech tensor
!           eigc        : eigenvalues of V obtained by spectral decomposition
!           eigprj      : eigenvectors of V obtained by spectral decomposition
!
! -----------------------------------------------------------------------------

    subroutine gauss_jordan(A,n,B,ierr)

        implicit none

        ! Auxiliar
        integer(4)                 :: i,j,k,l,ll,irow,icol,ierr
        integer(4),parameter       :: nmax=5000
        integer(4),dimension(nmax) :: indxc,indxr,ipiv
        real(8)                    :: dum,pivinv,big
        
        ! Input
        integer(4)                 :: n
        real(8),dimension(n)       :: B
        real(8),dimension(n,n)     :: A

        ! Output
        ! real(8),dimension(n,n)   :: A


        ! ---------------------------------------------------------------------

        ierr = 0   
        if (n .gt. nmax ) then
            ierr = 2
            return
        end if
        ipiv = 0
        do i=1,n
            big = 0.0d0
            do j=1,n
                if ( ipiv(j) .ne. 1 ) then
                    do k = 1,n
                        if ( ipiv(k) .eq. 0 ) then
                            if ( abs(A(j,k)) .ge. big ) then
                                big = abs(A(j,k))
                                irow = j
                                icol = k
                            end if
                        else if ( ipiv(k) .gt. 1 ) then
                            ierr = 1
                        return
                        end if
                    end do
                end if
            end do
            ipiv(icol) = ipiv(icol)+1
            if ( irow .ne. icol) then
                do l = 1,n
                    dum = A(irow,l)
                    A(irow,l) = A(icol,l)
                    A(icol,l) = dum
                end do
                dum = B(irow)
                B(irow) = B(icol)
                B(icol) = dum
            end if
            indxr(i) = irow
            indxc(i) = icol
            if ( A(icol,icol) .eq. 0.0d0 ) then
                ierr = 1
            end if
            pivinv = 1.0d0/A(icol,icol)
            A(icol,icol) = 1.0d0
            do l = 1,n
                A(icol,L) = A(icol,L)*pivinv
            end do
            B(icol) = B(icol)*pivinv
            do ll = 1,n
                if ( ll .ne. icol ) then
                    dum = A(ll,icol)
                    A(ll,icol) = 0.0d0
                    do l = 1,n
                        A(ll,l) = A(ll,l)-A(icol,l)*dum
                    end do
                    B(ll) = B(ll)-B(icol)*dum
                end if
            end do
        end do
        do l = n,1,-1
            if ( indxr(l) .ne. indxc(l) ) then
                do k = 1,n
                    dum = A(k,indxr(l))
                    A(k,indxr(l)) = A(k,indxc(l))
                    A(k,indxc(l)) = dum
                end do
            end if
        end do

    end subroutine

! -----------------------------------------------------------------------------

    function fdeterm2(amatr) result(deter)

        real(8),dimension(2,2) :: amatr
        real(8)                :: deter

        deter = 0.0d0
        deter = amatr(1,1)*amatr(2,2) - amatr(1,2)*amatr(2,1)

    end function

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       polar_decomposition()
!   
! -----------------------------------------------------------------------------
!
!   Description
!
!       . The polar decomposition of the Deformation Gradient tensor is
!    performed. The Left and Right Cauchy strech tensor and the rotation
!    tensor are determined.
!
!       . The respective procedure used her is described in the book
!    "Computational Methods for Plasticity" Souza Neto, Peric and Owen, pp.735.
!
!       . The subroutine has been tested with the examples of
!    "www.continuummechanics.org".
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Polar decomposition
!       . Determination of the rotation tensor
!       . Determination of the logarithmic strain tensor
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           ndi         : number of dimensions = 2
!           F           : deformation gradient
!
!       . Output:
!           Mrot        : rotation tensor
!           auxStrain   : logarithmic strain tensor
!
!       . Auxiliar:
!           U           : right Cauchy strech tensor
!           V           : left Cauchy strech tensor
!           eigc        : eigenvalues of V obtained by spectral decomposition
!           eigprj      : eigenvectors of V obtained by spectral decomposition
!
! -----------------------------------------------------------------------------

    subroutine polar_decomposition(F,Mrot,auxStrain)

        implicit none

        ! Auxiliar
        integer(4)                   :: i,j,k,idim,icomp
        integer(4),parameter         :: ndim=2
        real(8)                      :: ueig,um1eig
        real(8),dimension(3)         :: um1vec(3),uvec(3)
        real(8),dimension(4)         :: cvec
        real(8),dimension(ndim)      :: eigc
        real(8),dimension(4,ndim)    :: eigprj
        real(8),dimension(ndim,ndim) :: um1,C,U,V,Vd
        logical                      :: dummy
        
        ! Input
        real(8),dimension(ndim,ndim) :: F

        ! Output
        real(8),dimension(ndim,ndim) :: Mrot,auxStrain

        ! ---------------------------------------------------------------------

        ! COMPUTE C := F^T F

        C=0.d0
        do i = 1,ndim
            do j = 1,ndim
                do k = 1,ndim
                    C(i,j)=C(i,j) + F(k,i)*F(k,j)
                end do
            end do
        end do

        ! PERFORM SPECTRAL DECOMPOSITION OF C

        cvec(1) = C(1,1)
        cvec(2) = C(2,2)
        cvec(3) = C(1,2)
        eigprj = 0.d0
        eigc = 0.d0
        call spectral_decomposition(eigprj,eigc,dummy,cvec)

        ! COMPUTE U:=(C)^1/2 and U^-1
        
        ! vector form
        uvec = 0.0d0
        um1vec = 0.0d0
        do idim = 1,ndim
            ueig = sqrt(eigc(idim))
            um1eig = 1.0d0/ueig
            do icomp = 1,3
                uvec(icomp) = uvec(icomp) + ueig*eigprj(icomp,idim)
                um1vec(icomp) = um1vec(icomp) + um1eig*eigprj(icomp,idim)
            end do
        end do
        
        ! matrix form
        U(1,1) = uvec(1)
        U(2,2) = uvec(2)
        U(1,2) = uvec(3)
        U(2,1) = uvec(3)
        um1(1,1) = UM1VEC(1)
        um1(2,2) = UM1VEC(2)
        um1(1,2) = UM1VEC(3)
        um1(2,1) = UM1VEC(3)
        
        ! COMPUTE ROTATION R := F U^-1

        Mrot = 0.d0
        do i = 1,ndim
            do j = 1,ndim
                do k = 1,ndim
                    Mrot(i,j)=Mrot(i,j) + F(i,k)*UM1(k,j)
                end do
            end do
        end do

        ! COMPUTE LEFT STRETCH TENSOR V := F R^T
        
        !Vd=matmul(U,transpose(Mrot))
        !V=matmul(Mrot,Vd)
        V = matmul(F,transpose(Mrot))
        
        ! PERFORM SPECTRAL DECOMPOSITION OF V

        cvec(1) = V(1,1)
        cvec(2) = V(2,2)
        cvec(3) = V(1,2)
        call spectral_decomposition(eigprj,eigc,dummy,cvec)

        uvec = 0.d0
        UM1VEC = 0.d0
        do idim = 1,ndim
            ueig = log(eigc(idim))
            do icomp = 1,3
                uvec(icomp) = uvec(icomp) + ueig*eigprj(icomp,idim)
            end do
        end do

        ! LOGARITHMIC STRAIN TENSOR (MATRIX FORM)

        auxStrain(1,1) = uvec(1)
        auxStrain(2,2) = uvec(2)
        auxStrain(1,2) = uvec(3)
        auxStrain(2,1) = uvec(3)

        return

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       spectral_decomposition()
!
! -----------------------------------------------------------------------------
!
!   Description
!
!       . The spectral decompostion of 2x2 symetric tensor is performed.
!
!       . The respective procedure used her is described in the book
!    "Computational Methods for Plasticity" Souza Neto, Peric and Owen, pp.735.
!
!       . To calculate the Deformation Gradient the shape functions of the 
!    chosen element are used with the local coordinates of the nodes [-1,1].
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Spectral decompostion
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           X           : tensor to compute spectral decomposition
!
!       . Output:
!           eigx        : eigenvalues of X
!           eigprj      : eigenvectors of X
!
!       . Auxiliar:
!           small       : 
!           repeat      :
!           small       :
!           eigvec      :
!
! -----------------------------------------------------------------------------

    subroutine spectral_decomposition(eigprj,eigx,repeat,X)

        implicit none

        ! Auxiliar
        integer(4)                    :: idir
        integer(4),parameter          :: mcomp=4,ndim=2
        real(8)                       :: trx,auxB,auxC,differ,amxeig
        real(8),parameter             :: small=1.0e-5
        real(8),dimension(ndim,ndim)  :: auxmtx,eigvec
        logical                       :: repeat

        ! Input
        real(8),dimension(mcomp)      :: X

        ! Output
        real(8),dimension(ndim)       :: eigx
        real(8),dimension(mcomp,ndim) :: eigprj

        ! ---------------------------------------------------------------------

        repeat = .false.

        ! COMPUTE EIGENVALUES OF X
        trx = X(1) + X(2)
        auxB = sqrt((X(1)-X(2))**2 + 4.0d0*X(3)*X(3))
        eigx(1) = 0.5d0*(trx+auxB)
        eigx(2) = 0.5d0*(trx-auxB)

        ! COMPUTE EIGENVECTORS TENSORS
        differ = abs(eigx(1)-eigx(2))
        amxeig = max( abs(eigx(1)),abs(eigx(2))) 
        if ( amxeig .ne. 0.0d0 ) differ = differ/amxeig

        if ( differ .lt. small ) then

            ! write(*,*) 'There is a problem with the polar decomposition'

            repeat = .true.

            auxmtx(1,1) = X(1)
            auxmtx(2,2) = X(2)
            auxmtx(1,2) = X(3)
            auxmtx(2,1) = auxmtx(1,2)
    
            call jacobi(auxmtx,eigx,eigvec,2)
            do idir = 1,2
                eigprj(1,idir) = eigvec(1,idir)*eigvec(1,idir)
                eigprj(2,idir) = eigvec(2,idir)*eigvec(2,idir)
                eigprj(3,idir) = eigvec(1,idir)*eigvec(2,idir)
                eigprj(4,idir) = 0.0d0
            end do
        else
            ! use closed formula to compute eigenvectors tensors
            do idir = 1,2
                auxB = eigx(idir)-trx
                auxC = 1.0d0/(eigx(idir)+auxB)
                eigprj(1,idir) = auxC*(X(1)+auxB)
                eigprj(2,idir) = auxC*(X(2)+auxB)
                eigprj(3,idir) = auxC*X(3)
                eigprj(4,idir) = 0.0d0
            end do
        end if

        return

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       jacobi()
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Spectral decomposition of a n-dimensional symmetric matrix
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           A           : matrix form of tensor
!           D           : eigenvalues of tensor 
!           n           : n-dimension of matrix
!
!       . Output:
!           V           : eigenvalues of X
!
!       . Auxiliar:
!
! -----------------------------------------------------------------------------

    subroutine jacobi(A,D,V,n)

        implicit double precision (A-H,O-Z)

        ! Auxiliar
        integer(4),parameter    :: mjiter=50,nmax=100
        real(8)                 :: SM
        real(8),parameter       :: toler=1.0e-12
        real(8),dimension(nmax) :: B,Z

        ! Input
        integer(4)              :: n
        real(8),dimension(n)    :: D
        real(8),dimension(n,n)  :: A

        ! Output
        real(8),dimension(n,n)  :: V

        ! ---------------------------------------------------------------------

        if ( n .gt. nmax ) then
            ! call ERRPRT('EI0025')
        end if

        V = 0.0d0
        do ip = 1,n
            V(ip,ip) = 1.0d0
        end do

        do ip = 1,n
            B(ip) = A(ip,ip)
            D(ip) = B(ip)
            Z(ip) = 0.0d0
        end do

        do i = 1,mjiter
            SM = 0.0d0
            do ip = 1,n-1
                do iq = ip+1,n
                    SM = SM + abs(A(ip,iq))
                end do
            end do
            if ( SM .lt. toler ) goto 999
            if ( I .lt. 4 ) then
                tresh = 0.2d0*SM/dble(n**2)
            else
                tresh = 0.0d0
            end if
            do ip = 1,n-1
                do iq = ip+1,n
                    G = 100.0d0*abs(A(ip,iq))
                    if ( (I .gt. 4) .and. (abs(D(ip))+G .eq. abs(D(ip)))&
                     .and. (abs(D(iq))+G .eq. abs(D(iq))) ) then
                        A(ip,iq) = 0.0d0
                    else if ( abs(A(ip,iq)) .gt. tresh ) then
                        H = D(iq) - D(ip)
                        if ( abs(H)+G .eq. abs(H)) then
                            T = A(ip,iq)/H
                        else
                            theta = 0.5d0*H/A(ip,iq)
                            T = 1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
                            if ( theta .lt. 0.0d0 ) T = -T
                        end if
                        C = 1.0d0/sqrt(1.0d0+T**2)
                        S = T*C
                        tau = S/(1.0d0+C)
                        H = T*A(ip,iq)
                        Z(ip) = Z(ip)-H
                        Z(iq) = Z(iq)+H
                        D(ip) = D(ip)-H
                        D(iq) = D(iq)+H
                        A(ip,iq) = 0.0d0
                        DO J = 1,ip-1
                            G = A(J,ip)
                            H = A(J,iq)
                            A(J,ip) = G-S*(H+G*tau)
                            A(J,iq) = H+S*(G-H*tau)
                        end do
                        do J = ip+1,iq-1
                            G = A(ip,J)
                            H = A(J,iq)
                            A(ip,J) = G-S*(H+G*tau)
                            A(J,iq) = H+S*(G-H*tau)
                        end do
                        do J = iq+1,N
                            G = A(ip,J)
                            H = A(iq,J)
                            A(ip,J) = G-S*(H+G*tau)
                            A(iq,J) = H+S*(G-H*tau)
                        end do
                        do J=1,N
                            G=V(J,ip)
                            H=V(J,iq)
                            V(J,ip)=G-S*(H+G*tau)
                            V(J,iq)=H+S*(G-H*tau)
                        end do
                    end if
                end do
            end do
            do ip = 1,N
                B(ip) = B(ip)+Z(ip)
                D(ip) = B(ip)
                Z(ip) = 0.0d0
            end do
        end do
        ! call ERRPRT('EE0005')
        999 continue

        return

    end subroutine