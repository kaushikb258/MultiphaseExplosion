                    subroutine primitive_gas(cons,prim,rhocut,nspec)

          implicit none
          integer :: nspec
          real(kind=8) :: cons(1:4+nspec), prim(1:5+nspec), rhocut
          real(kind=8) :: ke, sumyk         
          integer :: ns

           if(cons(1).le.rhocut) then
            print*, 'rho too small ', cons(:)
            call flush()
            stop
           endif

          prim(1) = cons(1)
          prim(2) = cons(2)/prim(1)
          prim(3) = cons(3)/prim(1)
  
          ke = 0.5d0*(prim(2)**2.0d0 + prim(3)**2.0d0)

          prim(4) = cons(4)/prim(1) - ke


          do ns = 1, nspec
           prim(5+ns) = cons(4+ns)/prim(1)
           prim(5+ns) = min(max(prim(5+ns),0.0d0),1.0d0)
          enddo


          sumyk = sum(prim(6:5+nspec))  
          

          do ns = 1, nspec
           prim(5+ns) = prim(5+ns)/sumyk
          enddo


!                if(prim(4).le.0.0d0) then
!                 print*, 'prim4 < 0 ', prim, cons
!                 call flush()
!                 stop
!                endif


           ! pressure
           call compute_pres(prim(1), prim(6:5+nspec), prim(4), prim(5))


                if(prim(5).le.0.0d0) then
                 print*, 'prim5 < 0 ', prim, cons
                 call flush()
                 stop
                endif


                    return
                    end subroutine

!-----------------------------------------------------------------------

                    subroutine primitive_solid(cons,prim,nvar2,rho2cut)

          implicit none
          integer :: nvar2, ns
          real(kind=8) :: cons(1:nvar2), prim(1:nvar2), rho2cut
          real(kind=8) :: ke         

           if(cons(1).le.rho2cut) then
            prim = 0.0d0 
            return
           endif

          prim(1) = cons(1)
          prim(2) = cons(2)/prim(1)
          prim(3) = cons(3)/prim(1)
  
          ke = 0.5d0*(prim(2)**2.0d0 + prim(3)**2.0d0)

          prim(4) = cons(4)/prim(1) - ke

          prim(5) = cons(5)/prim(1)
          prim(6) = cons(6)/prim(1)


                    return
                    end subroutine

!-----------------------------------------------------------------------
