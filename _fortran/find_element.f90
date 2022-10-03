!
!   Subroutine
!
!       find_element()
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

    subroutine fin_element(ipt,inc)
        
        implicit none
        integer(4),intent(IN)               ::ipt,inc
        integer(4),intent(OUT)              ::iel
        integer(4)                          ::i
        real(8)                             ::ElArea, PTArea, A1, A2, A3, A4
        real(8),dimension(3,2)              ::v1,v2,v3,v4
        
        
        do i = 1,nelems
            ! AREA OF VIRTUAL ELEMENT
            ElArea = 0.0d0
            v1 = 0.0d0
            v1(1,:) = nodes(elems(i,1),:)
            v1(2,:) = nodes(elems(i,2),:)
            v1(3,:) = nodes(elems(i,3),:)
            call AreaTri(v1,A1)
            v2 = 0.0d0
            v2(1,:)=nodes(elems(i,1),:)
            v2(2,:)=nodes(elems(i,4),:)
            v2(3,:)=nodes(elems(i,3),:)
            call AreaTri(v2,A2)
            ElArea = A1 + A2
            !--------------------------------------------------------------------------------
            ! Area of the triangles built from each edge of the element and the data point
            !--------------------------------------------------------------------------------
            PTArea = 0.0d0
            iel = i
            v1 = 0.0d0
            v1(1,:)=nodes(elems(i,1),:)
            v1(2,:)=nodes(elems(i,2),:)
            v1(3,:)=centroid(inc,ipt,:)
            call AreaTri(v1,A1)
            v2 = 0.0d0
            v2(1,:)=nodes(elems(i,2),:)
            v2(2,:)=nodes(elems(i,3),:)
            v2(3,:)=centroid(inc,ipt,:)
            call AreaTri(v2,A2)
            v3 = 0.0d0
            v3(1,:)=nodes(elems(i,3),:)
            v3(2,:)=nodes(elems(i,4),:)
            v3(3,:)=centroid(inc,ipt,:)
            call AreaTri(v3,A3)
            v4 = 0.0d0
            v4(1,:)=nodes(elems(i,4),:)
            v4(2,:)=nodes(elems(i,1),:)
            v4(3,:)=centroid(inc,ipt,:)
            call AreaTri(v4,A4)
            PTArea = A1 + A2 + A3 + A4
            if(PTArea <= ElArea*1.01d0)then
                iel = i
                goto 1
            elseif(PTArea .gt. ElArea .and. i==nelems)then
                write(*,*)'Warning: data point does not belong to any element of the mesh'
                iel = -1
            end if
            continue
        end do
    1   continue
    end subroutine

! -----------------------------------------------------------------------------
