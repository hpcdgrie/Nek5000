to set actual boundary conditions
in usrdat2()


       do iel=1,nelt
        do ifc=1,2*ndim 
          id_face = bc(5,ifc,iel,1)
          bc(5,ifc,iel,2) = bc(5,ifc,iel,1)
          if (id_face.eq.1) then           ! for non-interface sideset
             cbc(ifc,iel,1) = 'W  '
             cbc(ifc,iel,2) = 't  '
          elseif (id_face.eq.2) then       ! for interface sideset 
             cbc(ifc,iel,1) = 'W  '
             cbc(ifc,iel,2) = 'E  '        ! this is optional
          endif
        enddo
        enddo
