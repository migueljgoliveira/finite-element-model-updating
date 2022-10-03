! -----------------------------------------------------------------------------
!
!   Table of Contents
!
!       . vfm()
!
!       . rotate_tensor()
!       . out_of_plane_strain()
!       . stress_pullback()
!
!       . gauss_jordan()
!
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       vfm()
!   
! -----------------------------------------------------------------------------
!
!   Description
!
!       . The virtual fields method is executed.
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Performs the virtual fields method.
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!
!       . Output:
!
!       . Auxiliar:
!
! -----------------------------------------------------------------------------

    subroutine vfm(n,m,incmax,nelems,ntests,vars,name,elprops,thickness,&
                   length,norm,residual,expforce,expstrain,maxexpforce,&
                   maxexpstrain,area,detdfgrd,dfgrd,invdfgrd,rot,flag,&
                   fvec,objfunc)

        implicit none

        ! Auxiliar
        integer(4)                               :: i,j,k,iel,iprops
        integer(4),parameter                     :: ndi=2,nshr=1,ntens=3,&
                                                    ncomp=4
        real(8)                                  :: SY0,H,SSat,CY,KA,lambda,&
                                                    phi1,phi2,phi,eqstress,&
                                                    q11,q22,q12,q66,weight
        real(8),dimension(n)                     :: x
        real(8),dimension(m)                     :: quadfvec
        real(8),dimension(:),allocatable         :: props
        real(8),dimension(3,1)                   :: vstrain,astress
        real(8),dimension(ntens,ntens)           :: matq
        real(8),dimension(incmax,nelems)         :: eqplas,strain33,hard
        real(8),dimension(incmax,nelems,ntens)   :: dstrain,estrain,pstrain,&
                                                    stress,dstress
        real(8),dimension(incmax,nelems,ntens+1) :: fpkstress
        character(200)                           :: path
        logical,parameter                        :: nlgeom=.True.

        ! Input
        integer(4)                               :: n,m,incmax,nelems,&
                                                    residual,ntests
        real(8)                                  :: thickness,length,&
                                                    maxexpforce,maxexpstrain
        real(8),dimension(n)                     :: vars
        real(8),dimension(2)                     :: elprops
        real(8),dimension(incmax)                :: expforce
        real(8),dimension(nelems)                :: area
        real(8),dimension(incmax,nelems)         :: detdfgrd
        real(8),dimension(incmax,nelems,ncomp)   :: dfgrd,invdfgrd,rot
        real(8),dimension(incmax,nelems,ntens)   :: expstrain
        character(200)                           :: name
        logical                                  :: norm,flag
        

        ! Output
        real(8)                                  :: objfunc
        real(8),dimension(m)                     :: fvec


        ! ---------------------------------------------------------------------

!f2py intent(in) n,m,incmax,nelems,ntests,vars,name,elprops,thickness,length
!f2py intent(in) norm,residual,expforce,expstrain,maxexpforce,maxexpstrain,area
!f2py intent(in) detdfgrd,dfgrd,invdfgrd,rot,flag
!f2py intent(out) objfunc,fvec
!f2py depend(n) vars
!f2py depend(incmax) expforce
!f2py depend(nelems) area
!f2py depend(incmax,nelems) detdfgrd
!f2py depend(incmax,nelems,ncomp) dfgrd,invdfgrd,rot
!f2py depend(incmax,nelems,ntens) expstrain
!f2py depend(2) elprops
!f2py depend(m) fvec

        ! ---------------------------------------------------------------------

        x = vars

        iprops = 8
        allocate(props(iprops))

        ! INITIALISE VARIABLES
        stress  = 0.0d0
        dstress = 0.0d0
        hard    = 0.0d0
        eqplas  = 0.0d0
        lambda= 0.0d0
        estrain= 0.0d0
        pstrain= 0.0d0
        strain33 = 0.0d0
        fpkstress = 0.0d0

        SY0  = (x(2)/x(1))**(1.0d0/x(3))
        H    = 0.0d0
        SSat = x(1)
        CY   = x(3)
        KA   = 0.0d0

        props(1) = elprops(1)
        props(2) = elprops(2)
        props(3) = SY0
        props(4) = thickness
        props(5) = H
        props(6) = SSat
        props(7) = CY
        props(8) = KA

        weight = 1.0d0/dble(incmax*ntests)

        ! INTERNAL AND EXTERNAL WORK
        phi1 = 0.0d0
        do i = 1,incmax
            if ( nlgeom .eqv. .False. ) then
                continue
            else if ( nlgeom .eqv. .True. ) then
                if ( i == 1 ) then
                    dstrain(i,:,:) = expstrain(i,:,:)
                else
                    dstrain(i,:,:) = expstrain(i,:,:) - expstrain(i-1,:,:)
                end if
            end if

            do j = 1,nelems
                ! UNDERSTAND ???
                ! iel = j
                ! if ( nlgeom == .False. ) then
                !     call find_element(j,iel,i)
                ! end if
                ! if ( iel == -1 )then
                !     write(*,*)'Data point skipped'
                !     !goto 100
                ! end if

                ! set elastic constitutive matrix
                q11 = props(1)/(1.0d0-props(2)**2.0d0)
                q22 = q11
                q12 = props(2)*props(1)/(1.0d0-props(2)**2.0d0)
                q66 = (q11-q12)/2.0d0

                matq = 0.0d0
                matq(1,1) = q11
                matq(2,2) = q22
                matq(1,2) = q12
                matq(2,1) = q12
                matq(3,3) = q66
                
                if ( i > 1 ) then
                    stress(i,j,:) = stress(i,j,:) + stress(i-1,j,:)
                    hard(i,j) = hard(i,j) + hard(i-1,j)
                    eqplas(i,j) = eqplas(i,j) + eqplas(i-1,j)
                end if
                if ( nlgeom .eqv. .True. ) then
                    if ( i /= 1 ) then
                        call rotate_tensor(stress(i,j,:),rot(i-1,j,:),0,0)
                    end if
                end if
                call ludwik_isotropic_hardening(ndi,nshr,ntens,iprops,props,&
                                                stress(i,j,:),dstress(i,j,:),&
                                                expstrain(i,j,:),&
                                                dstrain(i,j,:),hard(i,j),&
                                                eqplas(i,j),matq,eqstress)
                if ( nlgeom .eqv. .True. ) then
                    if ( i > 1 ) then
                        lambda = eqplas(i,j) - eqplas(i-1,j)
                        estrain(i,j,:) = estrain(i-1,j,:)
                        pstrain(i,j,:) = pstrain(i-1,j,:)
                        strain33(i,j) = strain33(i-1,j)
                    else
                        lambda = eqplas(i,j)
                    end if
                    call out_of_plane_strain(props(2),lambda,rot(i,j,:),&
                                             eqstress,stress(i,j,:),&
                                             estrain(i,j,:),pstrain(i,j,:),&
                                             dstrain(i,j,:),strain33(i,j))
                    call rotate_tensor(stress(i,j,:),rot(i,j,:),1,0)
                    call stress_pullback(stress(i,j,:),detdfgrd(i,j),&
                                         invdfgrd(i,j,:),strain33(i,j),&
                                         fpkstress(i,j,:))
                end if
            end do

            ! Internal Forces
            phi1 = 0.0d0
            do j = 1,nelems
                vstrain = 0.0d0
                vstrain(1,1) = 1.0d0
                do k = 1,ntens
                    astress(k,1) = stress(i,j,k)
                end do

                if ( nlgeom .eqv. .True. ) then
                    phi1 = phi1 + area(j) * ( &
                           fpkstress(i,j,1)*vstrain(1,1) + &
                           fpkstress(i,j,2)*vstrain(3,1) + &
                           fpkstress(i,j,3)*vstrain(3,1) + &
                           fpkstress(i,j,4)*vstrain(2,1) )
                else
                    do k = 1,ntens
                        phi1 = phi1 + vstrain(k,1)*astress(k,1)*area(j)
                    end do
                end if
            end do

            ! External Forces
            phi2 = (expforce(i)/thickness)*(length/2.0d0)

            ! Compute Residuals
            if ( norm .eqv. .False. ) then
                phi = phi1-phi2
            else
                phi = (phi1-phi2)/phi2
            end if

            if ( residual == 0 ) then
                fvec(i) = phi
            else if ( residual == 1 ) then
                fvec(i) = sqrt(weight)*phi
            else if ( residual == 2 ) then
                fvec(i) = weight*phi**2
            end if
        end do

        ! STORE OBJECTIVE FUNCTION VALUE
        if ( residual == 0 ) then
            quadfvec = fvec**2
            objfunc = objfunc + weight*sum(quadfvec)
        else if ( residual == 1 ) then
            quadfvec = fvec**2
            objfunc = objfunc + sum(quadfvec)
        else if ( residual == 2 ) then
            objfunc = objfunc + sum(fvec)
        end if

        objfunc = 0.5d0*objfunc

        ! SAVE OUTPUT FIELDS
        if ( flag ) then
            path = 'Output/tmp/'//'VFM_NumStress_'//trim(name)//'.dat'
            open(unit=6,file=trim(path),status='NEW')

            do i = 1,incmax
                do j = 1,nelems
                    write(6,*) (stress(i,j,k), k = 1,3)
                end do
            end do

            close(6)
        end if

        return
      
    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       rotate_tensor()
!
! -----------------------------------------------------------------------------
!
!   Description
!
!       . 2nd order tensor is rotated for the global or material reference
!       system.
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Rotates 2nd order tensor.
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           . tensor: 2nd order tensor
!           . auxrot: rotation tensor
!           . dir: coordinates system (0: global / 1: material)
!           . type: stress or strain tensor (0: stress / 1: strain)
!
!       . Output:
!           . tensor: rotated 2nd order tensor
!
!       . Auxiliar:
!           . auxtensor: matrix form of 2nd order tensor
!           . auxrot: matrix form of rotation tensor
!           . auxmatrix: auxiliar matrix
!
! -----------------------------------------------------------------------------

    subroutine rotate_tensor(tensor,rot,dir,tipo)

        implicit none

        ! Auxiliar
        integer(4),parameter     :: ntens=3,ncomp=4
        real(8),dimension(2,2)   :: auxtensor,auxrot,auxmatrix

        ! Input
        integer(4)               :: dir,tipo
        real(8),dimension(ntens) :: tensor
        real(8),dimension(ncomp) :: rot
        
        ! Output
        !real(8),dimension(ntens) :: tensor

        ! ---------------------------------------------------------------------
        
        auxtensor(1,1) = tensor(1)
        auxtensor(2,2) = tensor(2)
        auxtensor(1,2) = tensor(3)
        auxtensor(2,1) = tensor(3)

        if ( tipo == 0 ) then
            auxtensor(1,2) = tensor(3)
            auxtensor(2,1) = tensor(3)
        else if ( tipo == 1 ) then
            auxtensor(1,2) = tensor(3)/2.d0
            auxtensor(2,1) = tensor(3)/2.d0
        end if

        if ( dir == 0 ) then
            auxrot(1,1) = rot(1)
            auxrot(2,1) = rot(2)
            auxrot(1,2) = rot(3)
            auxrot(2,2) = rot(4)
        else if ( dir == 1 ) then
            auxrot(1,1) = rot(1)
            auxrot(1,2) = rot(2)
            auxrot(2,1) = rot(3)
            auxrot(2,2) = rot(4)
        end if

        auxmatrix = matmul(auxrot,auxtensor)
        auxtensor = matmul(auxmatrix,transpose(auxrot))

        tensor(1) = auxtensor(1,1)
        tensor(2) = auxtensor(2,2)

        if (tipo == 0 ) then
            tensor(3) = auxtensor(1,2)
        else if (tipo == 1 ) then
            tensor(3) = auxtensor(1,2)*2.d0
        end if

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       out_of_plane_strain()
!
! -----------------------------------------------------------------------------
!
!   Description
!
!       . Computes the out-of-plane normal strain component.
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           . nu: Poisson's coefficient
!           . lambda: plastic strain-rate multiplier
!           . eqstress: equivalent stress
!           . stress: Cauchy stress tensor in the material frame
!           . estrain: elastic deformation tensor in the global frame
!           . pstrain: plastic deformation tensor in the global frame
!           . dstrain: increment deformation tensor
!
!       . Output:
!           . strain33: out-of-plane strain normal component
!
!       . Auxiliar:
!           . ntens: number of tensor3 components
!           . ncomp: number of tensor4 components
!
! -----------------------------------------------------------------------------

    subroutine out_of_plane_strain(nu,lambda,rot,eqstress,stress,estrain,&
                                   pstrain,dstrain,strain33)

        ! Auxiliar
        integer(4),parameter     :: ntens=3,ncomp=4

        ! Input
        real(8)                  :: nu,eqstress,lambda
        real(8),dimension(ntens) :: stress,estrain,pstrain,dstrain,&
                                    destrain,dpstrain
        real(8),dimension(ncomp) :: rot

        ! Output
        real(8) :: strain33

        ! ---------------------------------------------------------------------
        if ( lambda > 0.0d0 ) then
            dpstrain = 0.0d0
            destrain = 0.0d0

            ! PLASTIC STRAIN-RATE TENSOR
            dpstrain(1) = lambda/(2.0d0*eqstress)*(2.0d0*stress(1)-stress(2))
            dpstrain(2) = lambda/(2.0d0*eqstress)*(2.0d0*stress(2)-stress(1))
            dpstrain(3) = lambda/(2.0d0*eqstress)*(6.0d0*stress(3))
            
            ! ELASTIC STRAIN-RATE TENSOR
            destrain = dstrain(:) - dpstrain

            ! ROTATE TENSORS TO MATERIAL FRAME
            call rotate_tensor(destrain,rot,1,1)
            call rotate_tensor(dpstrain,rot,1,1)
            call rotate_tensor(estrain,rot,1,1)
            call rotate_tensor(pstrain,rot,1,1)
            estrain(:) = estrain(:) + destrain
            pstrain(:) = pstrain(:) + dpstrain

            ! OUT-OF-PLANE STRAIN
            strain33 = strain33 + (-nu/(1-nu))*(destrain(1)+destrain(2)) - &
                       (dpstrain(1)+dpstrain(2))

            ! ROTATE TENSORS TO GLOBAL FRAME
            call rotate_tensor(estrain,rot,0,1)
            call rotate_tensor(pstrain,rot,0,1)
        else
            ! ROTATE TENSORS TO MATERIAL FRAME
            call rotate_tensor(dstrain,rot,1,1)
            call rotate_tensor(estrain,rot,1,1)
            estrain(:) = estrain(:) + dstrain(:)

            ! OUT-OF-PLANE STRAIN
            strain33 = strain33 + (-nu/(1-nu))*(dstrain(1)+dstrain(2))
        
            ! ROTATE TENSORS TO GLOBAL FRAME
            call rotate_tensor(estrain,rot,0,1)
        end if
    
    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       stress_pullback()
!
! -----------------------------------------------------------------------------

    subroutine stress_pullback(stress,detdfgrd,invdfgrd,strain33,fpkstress)

        ! Auxiliar
        real(8)              :: volr

        ! Input
        real(8)              :: detdfgrd,strain33
        real(8),dimension(3) :: stress
        real(8),dimension(4) :: invdfgrd
        
        ! Output
        real(8),dimension(4) :: fpkstress

        ! ---------------------------------------------------------------------
        
        volr = detdfgrd*(1.0d0+strain33)

        fpkstress(1) = volr * (stress(1)*invdfgrd(1) + stress(3)*invdfgrd(2))
        fpkstress(2) = volr * (stress(1)*invdfgrd(3) + stress(3)*invdfgrd(4))
        fpkstress(3) = volr * (stress(3)*invdfgrd(1) + stress(2)*invdfgrd(2))
        fpkstress(4) = volr * (stress(3)*invdfgrd(3) + stress(2)*invdfgrd(4))
    
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
