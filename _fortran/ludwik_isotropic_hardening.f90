!
!   Subroutine
!
!       ludwik_isotropic_hardening()
!   
! -----------------------------------------------------------------------------
!
!   Description
!
!       . 
!
! -----------------------------------------------------------------------------
!
!   Operations
!
!       . 
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

    subroutine ludwik_isotropic_hardening(ndi,nshr,ntens,iprops,props,stress,dstress,strain,dstrain,hard,eqplas,matD,SVM)
        
        implicit none
        
        real(8),parameter                                       ::YieldTol=1.0d-5
        integer(4),intent(IN)                                   ::iprops,ndi,nshr,ntens
        real(8),dimension(iprops),intent(IN)                    ::props
        real(8),dimension(1,ntens),intent(INOUT)                ::stress,dstress
        real(8),dimension(1,ntens),intent(IN)                   ::strain,dstrain
        real(8),dimension(1,1),intent(INOUT)                    ::hard,eqplas
        real(8),dimension(ntens,ntens),intent(INOUT)            ::matD
        real(8),intent(OUT)                                     ::SVM
        integer(4)                                              ::j,k1,k2,errorflag,counter
        real(8),dimension(ntens,1)                              ::stressold,stressin
        real(8),dimension(ntens)                                ::flow,temp1
        real(8)                                                 ::EMod,ve,esp,SY0,H
        real(8)                                                 ::eg,ek,YieldFunc,SYield,dp,yieldout
        real(8),dimension(ntens,1)                              ::DEP,rr,stressrate
        real(8),dimension(ntens,ntens)                          ::dadS,matQ,Qinv,matR
        real(8),dimension(ntens,ntens)                          ::matDInc,MatDCons
        real(8)                                                 ::param,dpr,YFS,deltadp
        real(8),dimension(1,1)                                  ::num,den
        real(8),dimension(1,3)                                  ::aux13,aux13b
        real(8),dimension(3,1)                                  ::aux31,aux31b
        real(8),dimension(3,3)                                  ::aux33,aux33b
        real(8),dimension(ntens,1)                              ::a
        real(8),dimension(ntens,ntens)                          ::matdelas
        real(8)                                                 ::kappa,expn
        !For consistent tangent stiffness
        real(8),dimension(ntens,ntens)                          ::Ident
        
        !Initialise variables -----
        Stressold=0.0d0
        MatDElas = MatD
        ident=0.0d0
        do k1=1,ntens
            ident(k1,k1)=1.0d0
        end do
        !Material Properties ------------------------
        EMod  = props(1)
        ve    = props(2)
        SY0   = abs(props(3))
        esp   = abs(props(4))
        H     = abs(props(5))
        kappa = abs(props(6))
        expn  = abs(props(7))
        !Shear Modulus ---------------------------------
        EG = Emod/(2.0d0*(1.0d0+ve))
        !Bulk Modulus for Plane Stress------------------
        ek = EMod/(2.0d0*(1.0d0-ve))
        !Save stress at the begining -------------------
        do j=1,ntens
            stressold(j,1)=stress(1,j)
        end do
        !Increment in strain --------------------------
        do j=1,ntens
            aux31(j,1) = dstrain(1,j)
        end do
        !Elastic increment in stress ------------------
        aux31b=matmul(MatD,aux31)
        do j=1,ntens
            dstress(1,j) = aux31b(j,1)
        end do
        !Trial Stress ---------------------------------
        do j=1,ntens
            stress(1,j) = stress(1,j) + dstress(1,j)
        end do
        !Store stress before returning to yield surface
        do j=1,ntens
            stressold(j,1) = stress(1,j)
        end do
        !Equivalent von Mises Stress ------------------
        do j=1,ntens
            stressin(j,1) = stress(1,j)
        end do
        call von_mises(stressin,ntens,SVM)
        !Determine if Actively Yielding ---------------
        !SYield = SY0 + kappa*eqplas(1,1)**(expn)
        SYield = kappa*(SY0+eqplas(1,1))**(expn)
        YieldFunc = SVM - SYield
        continue
        !Yielding Condition NOT satisfied -------------
        counter=0
        if (YieldFunc .gt. YieldTol)then
            deltadp=0.0d0
            dpr=0.0d0
            dp=0.0d0
            YFS=YieldFunc
            YieldOUT=YieldFunc/SVM*100.0d0
            !Flow direction (from Crisfield)
            flow(1)=2.0d0*stress(1,1)-stress(1,2)
            flow(2)=2.0d0*stress(1,2)-stress(1,1)
            flow(3)=6.0d0*stress(1,3)
            flow=flow/(2.0d0*SVM)
            do while(abs(YieldFunc) .Gt. YieldTol)
                counter=counter+1
                continue
                !Determine initial plastic multiplier and Yield stress
                if (counter==1)then
                    param=0.0d0
                    do k1=1,ntens
                        do k2=1,ntens
                            param=param + flow(k1)*matdelas(k1,k2)*flow(k2)
                        end do
                    end do
                    !Determination of plastic multiplier and new yield stress ------------
                    dp = 0.0d0
                    !H = kappa*expn*(eqplas(1,1)+dp)**(expn - 1.0d0)
                    H = kappa*expn*(SY0+eqplas(1,1))**(expn - 1.0d0)
                    !hard(1,1) = kappa*(eqplas(1,1)+dp)**(expn)
                    hard(1,1) = kappa*( SY0 + eqplas(1,1))**(expn)
                    dp= YieldFunc/(param+hard(1,1))
                    !Determine backward stress -------------------------------------------
                    DEP=0.0d0
                    do k1=1,ntens
                        do k2=1,ntens
                            DEP(k1,1)=DEP(k1,1) + matdelas(k1,k2)*flow(k2)
                        end do
                    end do
                    DEP=DEP*dp
                    !Update stress -------------------------------------------------------
                    do j=1,ntens
                        stress(1,j)=stress(1,j) - DEP(j,1)
                    end do
                    !Yield surface verification ------------------------------------------
                    do j=1,ntens
                        stressin(j,1)=stress(1,j)
                    end do
                    call von_mises(stressin,ntens,SVM)
                    SYield = kappa*(SY0+eqplas(1,1)+dp)**(expn)
                    YieldFunc = SVM - SYield
                    H = kappa*expn*(SY0 + eqplas(1,1) + dp)**(expn - 1.0d0)
                    hard(1,1) = kappa*(SY0 + eqplas(1,1) + dp)**(expn)
                end if
                flow(1)=2.0d0*stress(1,1)-stress(1,2)
                flow(2)=2.0d0*stress(1,2)-stress(1,1)
                flow(3)=6.0d0*stress(1,3)
                flow=flow/(2.0d0*SVM)
                DEP=0.0d0
                do k1=1,ntens
                    do k2=1,ntens
                        DEP(k1,1)=DEP(k1,1) + matdelas(k1,k2)*flow(k2)
                    end do
                end do
                DEP=DEP*dp
                !Difference between current stress and bakward-Euler stress -----
                do j=1,ntens
                    rr(j,1)=stress(1,j) - (stressold(j,1) - DEP(j,1))
                end do
                dadS=0.0d0
                dadS(1,1)=2.0d0
                dadS(1,2)=-1.0d0
                dadS(2,1)=-1.0d0
                dadS(2,2)=2.0d0
                dadS(3,3)=6.0d0
                dadS=dadS/(2.0d0*SVM)
                a=0.0d0
                do k1=1,ntens
                    do k2=1,ntens
                        dadS(k1,k2)=dadS(k1,k2)-flow(k1)*flow(k2)/svm
                    end do
                    a(k1,1)=flow(k1)
                end do
                !Change in stress -----
                matQ=(ident + dp*matmul(matdelas,dadS))
                Qinv= MatQ
                temp1 = 0.0d0
                call gauss_jordan(Qinv,ntens,temp1,errorflag)
                num=0.0d0
                den=0.0d0
                aux13=matmul(transpose(a),Qinv)
                num=matmul(aux13,rr)
                aux13b=matmul(aux13,matdelas)
                den=matmul(aux13b,a)
                !Plastic multiplier rate ------------------------------------
                dpr=0.0d0
                dpr=(YieldFunc-num(1,1))/(den(1,1) + H)
                aux31 = matmul(Qinv,rr)
                aux31b = matmul(matmul(Qinv,matdelas),a)
                !Stress rate ------------------------------------------------
                stressrate=-1.0d0*aux31 - dpr*aux31b
                do k1=1,ntens
                    stress(1,k1)=stress(1,k1) + stressrate(k1,1)
                end do
                !end if
                !---------------------------------------------------------------------
                !Yield surface verification ------------------------------------------
                do j=1,ntens
                    stressin(j,1)=stress(1,j)
                end do
                call von_mises(stressin,ntens,SVM)
                dp= dp + dpr
                !H = kappa*expn*(eqplas(1,1)+dp)**(expn - 1.0d0)
                H = kappa*expn*(SY0+eqplas(1,1)+dp)**(expn - 1.0d0)
                !Update yield surface ----
                !hard(1,1) = kappa*(eqplas(1,1)+dp)**(expn)
                hard(1,1) = kappa*(SY0+eqplas(1,1)+dp)**(expn)
                !SYield = SY0 + kappa*(eqplas(1,1)+dp)**(expn)
                SYield = kappa*(SY0+eqplas(1,1)+dp)**(expn)
                !Yield Function ----
                YieldFunc=SVM-SYield
                YieldOUT=abs(SVM-SYield)/SVM*100.0d0
                continue
                !---------------------------------------------------------------------
            end do
            eqplas(1,1) = eqplas(1,1) + dp
            !Inconsistent tangent matrix taken from Crisfield -------
            aux33=matmul(matmul(a,transpose(a)),matdelas)
            den=matmul(matmul(transpose(a),matdelas),a)
            aux33b=0.0d0
            aux33b=Ident - aux33/(den(1,1) + H)
            matDInc=matmul(matdelas,aux33b)
            !Consistent tangent matrix taken from Crisfield -------
            dadS=0.0d0
            dadS(1,1)=2.0d0
            dadS(1,2)=-1.0d0
            dadS(2,1)=-1.0d0
            dadS(2,2)=2.0d0
            dadS(3,3)=6.0d0
            dadS=dadS/(2.0d0*SVM)
            a=0.0d0
            do k1=1,ntens
                do k2=1,ntens
                    dadS(k1,k2)=dadS(k1,k2)-flow(k1)*flow(k2)/svm
                end do
                a(k1,1)=flow(k1)
            end do
            matQ=(ident + dp*matmul(matdelas,dadS))
            Qinv= MatQ
            temp1 = 0.0d0
            call gauss_jordan(Qinv,ntens,temp1,errorflag)
            matR=0.0d0
            matR=matmul(Qinv,matdelas)
            aux33 = matmul(a,transpose(a))
            aux33b = matmul(aux33,matR)
            den=0.0d0
            den=matmul(matmul(transpose(a),matR),a)
            aux33=0.0d0
            aux33=Ident - aux33b/(den(1,1) + H)
            matDCons=matmul(matR,aux33)
            matd=matDCons
            continue
        endif
        continue
    
    end subroutine

    subroutine von_mises(stress,sdim,SVM)
        implicit none
        
        integer(4),intent(IN)::sdim
        real(8),dimension(sdim,1)::stress
        real(8),intent(OUT)::SVM
        
        SVM=0.0d0
        
        if(sdim==3)then
            SVM=sqrt(stress(1,1)**2.0d0+stress(2,1)**2.0d0-stress(1,1)*stress(2,1)+3.0d0*stress(3,1)**2.0d0)
        elseif(sdim==6)then
            SVM=(stress(1,1)-stress(2,1))**2+(stress(2,1)-stress(3,1))**2+(stress(1,1)-stress(3,1))**2
            SVM=SVM+6.0d0*(stress(4,1)**2+stress(5,1)**2+stress(6,1)**2)
            SVM=SVM/2.0d0
            SVM=sqrt(SVM)
        end if

        continue

    end subroutine