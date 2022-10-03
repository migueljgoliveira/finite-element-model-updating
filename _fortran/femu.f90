! -----------------------------------------------------------------------------
!
!   Table of Contents
!
!       . femu()
!
!       . strain_difference()
!       . force_difference()
!
! -----------------------------------------------------------------------------

! -----------------------------------------------------------------------------
!   Subroutine
!
!       femu()
!
! -----------------------------------------------------------------------------
!
!   Description
!
!       . The finite element model updating method is executed, using Abaqus
!       software. 
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . Performs the finite element model updating method
!
! -----------------------------------------------------------------------------
!
!   Variables
!   
!       . Input:
!           n            : dimension of design variables
!           m            : dimension of residuals
!           incmax       : maximum experimental increments
!           nelems       : number of elements
!           ntests       : number of tests
!           cpus         : number of cpus for Abaqus job
!           vars         : design variables
!           name         : name of test
!           norm         : objective function normalised type
!           residual     : objective function residual type
!           femuresidual : femu residual type
!           exptime      : experimental time
!           expforce     : experimental force
!           expstrain    : experimental strain
!           flag         : flag to store output fields
!
!       . Output:
!           objfunc      : objective function value
!           fvec         : residuals
!
!       . Auxiliar:
!           incnum       : number of numerical increments
!           nstrain      : number of strain residuals
!           nforce       : number of force residuals
!           ntens        : number of tensor components
!           x            : copy of design variables
!           quadfvec     : quadratic residuals
!           wquadfvec    : weighted quadratic residuals
!           numforce     : interpolated numerical force
!           numstrain    : interpolated numerical strain
!           auxnumforce  : numerical force
!           auxnumstrain : numerical strain
!
! -----------------------------------------------------------------------------

    subroutine femu(m,incmax,incnum,nelems,ntests,norm,residual,femuresidual,&
                    exptime,expforce,expstrain,numtime,auxnumforce,&
                    auxnumstrain,maxexpforce,maxexpstrain,fvec,objfunc)

        implicit none

        ! Auxiliar
        integer(4)                             :: i,j,k,mc,nstrain,nforce
        integer(4),parameter                   :: ntens=3
        real(8)                                :: auxb,auxa,auxc,weight1,&
                                                  weight2,phi,phi1,phi2
        real(8),dimension(m)                   :: quadfvec,wquadfvec
        real(8),dimension(incmax)              :: numforce
        real(8),dimension(incmax,nelems,ntens) :: numstrain
        logical                                :: equal

        ! Input
        integer(4)                             :: m,incmax,incnum,nelems,&
                                                  ntests,residual,femuresidual
        real(8)                                :: maxexpforce,maxexpstrain
        real(8),dimension(incmax)              :: exptime,expforce
        real(8),dimension(incnum)              :: numtime,auxnumforce
        real(8),dimension(incmax,nelems,ntens) :: expstrain
        real(8),dimension(incnum,nelems,ntens) :: auxnumstrain
        logical                                :: norm

        ! Output
        real(8)                                :: objfunc
        real(8),dimension(m)                   :: fvec

        ! ---------------------------------------------------------------------

!f2py intent(in) m,incmax,incnum,nelems,ntests,norm,residual,femuresidual
!f2py intent(in) exptime,expforce,expstrain,numtime,auxnumforce,auxnumstrain
!f2py intent(in) maxexpforce,maxexpstrain
!f2py intent(out) fvec,objfunc
!f2py depend(m) fvec
!f2py depend(incmax) exptime,expforce
!f2py depend(incnum) numtime,auxnumforce
!f2py depend(incmax,nelems,ntens) expstrain
!f2py depend(incnum,nelems,ntens) auxnumstrain

        ! ---------------------------------------------------------------------

        ! LINEAR INTERPOLATION OF NUMERICAL DATA
        if ( incnum .ne. incmax ) then
            do i = 1,incmax
                equal = .False.
                do j = 1,incnum
                    if ( exptime(i) == numtime(j) ) then
                        equal = .True.
                        goto 10
                    end if
                end do
                10 continue
                if ( .not. equal) then
                    k = 1
                    do while ( exptime(i) .gt. numtime(k) )
                        k = k + 1
                    end do
                    auxb = numtime(k) - numtime(k-1)
                    auxa = exptime(i) - numtime(k-1)
                    auxc = auxa/auxb
                    numstrain(i,:,:) = auxnumstrain(k-1,:,:) + &
                                       (auxnumstrain(k,:,:) - &
                                       auxnumstrain(k-1,:,:)) * auxc
                    numforce(i) = auxnumforce(k-1) + &
                                  (auxnumforce(k) - auxnumforce(k-1)) * auxc
                else
                    numstrain(i,:,:) = auxnumstrain(j,:,:)
                    numforce(i) = auxnumforce(j)
                end if
            end do
        else
            numforce = auxnumforce
            numstrain = auxnumstrain
        end if

        ! COMPUTE STRAIN AND FORCE QUADRATIC DIFFERENCES
        fvec = 0.0d0

        nstrain = ntens*nelems*incmax*ntests
        weight1 = 1.0d0/dble(nstrain)
        nforce = incmax*ntests
        weight2 = 1.0d0/dble(nforce)

        mc = 0

        ! Residuals for Each Increment (Strain and Force)
        if ( femuresidual == 0 ) then
            do i = 1,incmax
                phi = 0.0d0
                do j = 1,nelems
                    do k = 1,ntens
                        call strain_difference(residual,norm,weight1,&
                                               numstrain(i,j,k),&
                                               expstrain(i,j,k),&
                                               maxexpstrain,phi1)
                        phi = phi + phi1
                    end do
                end do
                call force_difference(residual,norm,weight2,numforce(i),&
                                      expforce(i),maxexpforce,phi1)
                phi = phi + phi1
                mc = mc + 1
                fvec(mc) = phi
            end do
        else
        ! Residuals for Each Increment (Strain)
            if ( femuresidual == 1 ) then
                do i = 1,incmax
                    phi = 0.0d0
                    do j = 1,nelems
                        do k = 1,ntens
                            call strain_difference(residual,norm,weight1,&
                                                   numstrain(i,j,k),&
                                                   expstrain(i,j,k),&
                                                   maxexpstrain,phi1)
                            phi = phi + phi1
                        end do
                    end do
                    mc = mc + 1
                    fvec(mc) = phi
                end do
        ! Residuals for Each Increment, Point and Component (Strain)
            else if ( femuresidual == 2 ) then
                do i = 1,incmax
                    do j = 1,nelems
                        do k = 1,ntens
                            call strain_difference(residual,norm,weight1,&
                                                   numstrain(i,j,k),&
                                                   expstrain(i,j,k),&
                                                   maxexpstrain,phi1)
                            mc = mc + 1
                            fvec(mc) = phi1
                        end do
                    end do
                end do
            end if
        ! Residuals for Each Increment (Force)
            do i = 1,incmax
                call force_difference(residual,norm,weight2,numforce(i),&
                                      expforce(i),maxexpforce,phi2)
                mc = mc + 1
                fvec(mc) = phi2
            end do
        end if

        ! STORE OBJECTIVE FUNCTION VALUE
        if ( residual == 0 ) then
            quadfvec = fvec**2
            do i = 1,nstrain
                wquadfvec(i) = weight1*quadfvec(i)
            end do
            do i = 1,nforce
                wquadfvec(nstrain+i) = weight2*quadfvec(nstrain+i)
            end do
            objfunc = objfunc + sum(wquadfvec)
        else if ( residual == 1 ) then
            quadfvec = fvec**2
            objfunc = objfunc + sum(quadfvec)
        else if ( residual == 2 ) then
            objfunc = objfunc + sum(fvec)
        end if

        objfunc = 0.5d0*objfunc

        return
      
    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       strain_difference()
!
! -----------------------------------------------------------------------------

    subroutine strain_difference(residual,norm,weight1,auxnumstrain,&
                                 auxexpstrain,maxexpstrain,auxphi)

        ! Auxiliar
        real(8)    :: difstrain

        ! Input
        real(8)    :: weight1,auxnumstrain,auxexpstrain,maxexpstrain
        integer(4) :: residual
        logical    :: norm

        ! Output
        real(8)    :: auxphi

        ! ---------------------------------------------------------------------
        
        auxphi = 0.0d0
        difstrain = auxnumstrain - auxexpstrain

        ! DIFFERENCE
        if ( residual == 0 ) then
            if ( norm .eqv. .False. ) then
                auxphi = difstrain
            else
                auxphi = difstrain/maxexpstrain
            end if
        ! WEIGHTED DIFFERENCE
        else if ( residual == 1 ) then
            if ( norm .eqv. .False. ) then
                auxphi = sqrt(weight1) * difstrain
            else
                auxphi = sqrt(weight1) * (difstrain/maxexpstrain)
            end if
        ! WEIGHTED QUADRATIC DIFFERENCE
        else if ( residual == 2 ) then
            if ( norm .eqv. .False. ) then
                auxphi = weight1 * difstrain**2
            else
                auxphi = weight1 * (difstrain/maxexpstrain)**2
            end if
        end if

    end subroutine

! -----------------------------------------------------------------------------
!
!   Subroutine
!
!       force_difference()
!
! -----------------------------------------------------------------------------

    subroutine force_difference(residual,norm,weight2,auxnumforce,&
                                auxexpforce,maxexpforce,auxphi)

        ! Auxiliar
        real(8)    :: difforce

        ! Input
        real(8)    :: weight2,auxnumforce,auxexpforce,maxexpforce
        integer(4) :: residual
        logical    :: norm

        ! Output
        real(8)    :: auxphi

        ! ---------------------------------------------------------------------
        
        auxphi = 0.0d0
        difforce = auxnumforce - auxexpforce

        ! DIFFERENCE
        if ( residual == 0 ) then
            if ( norm .eqv. .False. ) then
                auxphi = difforce
            else
                auxphi = difforce/maxexpforce
            end if
        ! WEIGHTED  DIFFERENCE
        else if ( residual == 1 ) then
            if ( norm .eqv. .False. ) then
                auxphi = sqrt(weight2) * difforce
            else
                auxphi = sqrt(weight2) * (difforce/maxexpforce)
            end if
        ! WEIGHTED QUADRATIC DIFFERENCE
        else if ( residual == 2 ) then
            if ( norm .eqv. .False. ) then
                auxphi = weight2 * difforce**2
            else
                auxphi = weight2 * (difforce/maxexpforce)**2
            end if
        end if
    
    end subroutine

! -----------------------------------------------------------------------------