                               program euler2d

           USE THERMO

              integer, parameter :: imax = 1001
              integer, parameter :: jmax = 1001
              integer, parameter :: nmax =  1 !50000 ! number of time steps
              integer, parameter :: ofile = 1 !250 ! how often to output vtk file             

              real(kind=8) :: dx, dy, x(1:imax), y(1:jmax), xlen, ylen
              real(kind=8) :: time, dt, cfl
              integer :: i, j, n, ns, irestart

              real(kind=8), parameter :: pi = acos(-1.0d0)
              integer, parameter :: nghost = 5 ! number of ghose cells
              integer, parameter :: nspec = 2 ! number of species
              integer, parameter :: ncons = 4+nspec ! number of conservative variables
              integer, parameter :: nprim = 5+nspec ! number of primitive variables
              integer, parameter :: nvar2 = 6 ! number of solid variables

              real(kind=8) :: cons(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: prim(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: xflux(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: yflux(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: xflux2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: yflux2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: cons2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: augmented_prim(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim+nvar2)

              real(kind=8) :: prim_intxl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: prim_intxr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intxl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: cons_intxr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: prim_intyl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: prim_intyr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nprim)
              real(kind=8) :: cons_intyl(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
              real(kind=8) :: cons_intyr(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)

              real(kind=8) :: prim_intxl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intxr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intxl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intxr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intyl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: prim_intyr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intyl2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)
              real(kind=8) :: cons_intyr2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2) 

              real(kind=8) :: d_ip12, d_im12, d_im32, eps1, kappa, phi1, phi2, d1, d2
              real(kind=8) :: spl, spr, wavespeed, terml, termr, termrl, etotl, etotr

              real(kind=8) :: gam, vel, cc, ke, xc, yc, rc, etot, theta, rc1
              real(kind=8) :: rhocut 
              real(kind=8) :: temp, molwt(nspec), molmix
              real(kind=8) :: Rgas
              real(kind=8) :: rhomin, rhomax, rhomid, pres
              integer :: iter
              real(kind=8) :: dist
              character ( len = 100 ) filename
              real(kind=8) :: gravity, source(1:ncons), omegadot(1:nspec)
              real(kind=8) :: maxvel, maxrho, maxp, minvel, minrho, minp
              real(kind=8) :: hfk(1:nspec), enth, ie, conc(1:nspec)   
              real(kind=8) :: rand1, rand2, rand3, rand4, rand5, rand6, rand7, rand8 

              integer :: rkstep 
              integer, parameter :: rk = 3
              real(kind=8) :: k_rk(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons,1:rk)
              real(kind=8) :: cons0(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,ncons)
            
              real(kind=8) :: k_rk2(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2,1:rk)
              real(kind=8) :: cons02(1-nghost:imax-1+nghost,1-nghost:jmax-1+nghost,nvar2)

              real(kind=8) :: dummy, T1, hfs
 
              integer, parameter :: space_scheme = 4 ! enter 1/2/3/4

              ! boundary conditions (1: periodic; 2: outflow; 3: slip wall)
              integer, parameter :: xbc_low = 3
              integer, parameter :: ybc_low = 3
              integer, parameter :: xbc_high = 2
              integer, parameter :: ybc_high = 2

              !--------------------------------------
              ! Variables for paraview 
              integer, parameter :: output_unit = 29
              character (len=100) title
              real(kind=8) :: xyz(1:imax,1:jmax,1:3)   
              !--------------------------------------


              real(kind=8) :: cvsolid, dia, rho2cut, rhosol
              real(kind=8) :: num, dia1, cd, nuss, sigmadot, T2, reyn, mu, mach 
              real(kind=8) :: fdrag(2), qheat, cpmix, cpk(1:nspec), lambda, prandtl  
              real(kind=8) :: Jkin, Jdiff, Zhyb, Ea_R, taup, muc, source_fign, tign
 
              real(kind=8) :: gaml, gamr, rhol, rhor, pl, pr, ul, ur, vl, vr, iel,   &
                   ier, ysl(1:nspec), ysr(1:nspec), flux(1:ncons)
              integer :: sweep
 
                 call thermo_data()


              ! conservative variables
              ! 1: rho
              ! 2: rho*u
              ! 3: rho*v 
              ! 4: rho*(e + u^2/2 + v^2/2)
              ! 5: rho*Y1
              ! 6: rho*Y2

              ! primitive variables
              ! 1: rho
              ! 2: u
              ! 3: v
              ! 4: e
              ! 5: p
              ! 6: Y1
              ! 7: Y2

              !----------------------

              ! Set values here
              xlen = 1.5d0    
              ylen = 1.5d0    
              cfl = 0.05d0 
              rhocut = 1.0d-4 ! cutoff rho
              temp = 293.0d0 ! temperature
              gravity = 0.0d0 !-1.0d0 !-9.8d0 ! gravity

              molwt = MW

              !----------------------
              ! MUSCL PARAMETERS
              
               eps1 = 1.0d0 
               kappa = 1.0d0/3.0d0 
              !----------------------

              ! SOLID PARAMETERS

              dia = 10.0d-6
              cvsolid = 766.0d0  
              rho2cut = 1.0d-4
              rhosol = 2500.0d0               
              prandtl = 0.71d0  

              Zhyb = 7.5d4
              Ea_R = 8500.0d0 

              !----------------------

 
              dx = xlen/dble(imax-1)
              dy = ylen/dble(jmax-1)
              print*, 'dx, dy: ', dx, dy 
              

              do i = 1, imax
               x(i) = dx*dble(i-1)
              enddo
              do j = 1, jmax
               y(j) = dy*dble(j-1)
              enddo

              ! initialize arrays
              cons = 0.0d0
              prim = 0.0d0


               open(4,file='restart.data',form='formatted')
                read(4,*) irestart
               close(4)


                   n = 0
                   time = 0.0d0

                 if(irestart.eq.0) then

             do i = 1, imax-1 !1-nghost, imax-1+nghost
              do j = 1, jmax-1 !1-nghost, jmax-1+nghost

                 xc = 0.5d0*(x(i+1)+x(i))
                 yc = 0.5d0*(y(j+1)+y(j)) 

                 theta = atan(yc/xc) 

                 rc = sqrt(xc**2.0d0 + yc**2.0d0)

               prim2(i,j,1) = 1.0d-3

               temp = 293.0d0 ! SURROUNDING TEMPERATURE

               prim(i,j,2) = 0.0d0 ! x-velocity
               prim(i,j,3) = 0.0d0 ! y-velocity

                !call random_number(dummy)
                !rc1 = 0.1d0*(1.0d0 + dummy*0.05d0*sin(128.0d0*theta/3.14159265d0))

                call random_number(rand1) 
                call random_number(rand2) 
                call random_number(rand3) 
                call random_number(rand4)
                call random_number(rand5) 
                call random_number(rand6) 
                call random_number(rand7) 
                call random_number(rand8)

                rc1 =  rand1*0.005d0*sin(128.0d0*theta/pi) 
                rc1 = rc1 + rand2*0.007d0*cos(128.0d0*theta/pi)
                rc1 = rc1 + rand3*0.01d0*cos(256.0d0*theta/pi) 
                rc1 = rc1 + rand4*0.008d0*cos(256.0d0*theta/pi)
                rc1 = rc1 + rand5*0.009d0*cos(512.0d0*theta/pi) 
                rc1 = rc1 + rand6*0.011d0*cos(512.0d0*theta/pi)
                rc1 = rc1 + rand7*0.0125d0*cos(1024.0d0*theta/pi)
                rc1 = rc1 + rand8*0.0083d0*cos(1024.0d0*theta/pi)
                rc1 = 0.1d0*(1.0d0 + rc1)
                



                if(rc.le.rc1) then 
               prim(i,j,5) = 100.0d5 ! pressure
               ! species mass fraction
               prim(i,j,6) = 1.0d0 ! AIR
               prim(i,j,7) = 0.0d0 ! AIR
                else
               prim(i,j,5) = 1.0d5 ! pressure
               ! species mass fraction
               prim(i,j,6) = 0.0d0 ! AIR
               prim(i,j,7) = 1.0d0 ! AIR
                endif


               ! mol wt of mixture
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 


               ! density
               prim(i,j,1) = prim(i,j,5)/Rgas/temp 

               ! internal energy
               call COMPUTE_HFK (TEMP, HFK) 
               enth = 0.0d0
               do ns = 1, nspec
                 enth = enth + prim(i,j,5+ns)*hfk(ns)
               enddo
               ie = enth - prim(i,j,5)/prim(i,j,1)
               prim(i,j,4) = ie

              
               prim2(i,j,2) = 0.0d0 ! x-velocity
               prim2(i,j,3) = 0.0d0 ! y-velocity
               if(prim2(i,j,1).gt.rho2cut) then
                prim2(i,j,4) = cvsolid*(293.0d0)
                num = prim2(i,j,1)/rhosol*6.0d0/3.14159265d0/(dia**3.0d0)
                prim2(i,j,6) = num ! number density
               else
                prim2(i,j,4) = 0.0d0 ! e2
                prim2(i,j,6) = 0.0d0 ! number density
               endif
               prim2(i,j,5) = 0.0d0 ! fign
 
             enddo
            enddo

                   else

              write(filename,'("output/output_",I5.5,".dat")'),irestart
              open(34,file=filename,form='formatted') 
               read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 
                 read(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4), prim2(i,j,5), prim2(i,j,6) 
                   prim2(i,j,4) = prim2(i,j,4)*cvsolid

               ! internal energy
               call COMPUTE_HFK (T1, HFK) 
               enth = 0.0d0
               do ns = 1, nspec
                 enth = enth + prim(i,j,5+ns)*hfk(ns)
               enddo
               ie = enth - prim(i,j,5)/prim(i,j,1)
               prim(i,j,4) = ie

                enddo 
               enddo 
              close(34)


                   endif


                  !---------------------------------------
                  !---------------------------------------
                  ! set primitive variables in ghost cells

                      if(xbc_low.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(imax-i,j,1:nprim)
                 prim2(1-i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)
                enddo
               enddo
                        else if(xbc_low.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(1,j,1:nprim)
                 prim2(1-i,j,1:nvar2) = prim2(1,j,1:nvar2)
                enddo
               enddo
                        else if(xbc_low.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(i,j,1:nprim)
                 prim2(1-i,j,1:nvar2) = prim2(i,j,1:nvar2)

                 ! u velocity should be opposite
                 prim(1-i,j,2) = -prim(i,j,2)
                 prim2(1-i,j,2) = -prim2(i,j,2)
                enddo
               enddo
                        endif


                      if(xbc_high.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(imax-1+i,j,1:nprim) = prim(i,j,1:nprim)
                 prim2(imax-1+i,j,1:nvar2) = prim2(i,j,1:nvar2)
                enddo
               enddo
                        else if(xbc_high.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(imax-1+i,j,1:nprim) = prim(imax-1,j,1:nprim)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-1,j,1:nvar2)
                enddo
               enddo
                        else if(xbc_high.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(imax-1+i,j,1:nprim) = prim(imax-i,j,1:nprim)
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)

                 ! u velocity should be opposite
                 prim(imax-1+i,j,2) = -prim(imax-i,j,2)
                 prim2(1-i,j,2) = -prim2(i,j,2)
                 prim2(imax-1+i,j,2) = -prim2(imax-i,j,2)
                enddo
               enddo
                        endif


                !-----------------------------------------------------

                        if(ybc_low.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,jmax-j,1:nprim)
                 prim2(i,1-j,1:nvar2) = prim2(i,jmax-j,1:nvar2)
                enddo
               enddo
                        else if(ybc_low.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,1,1:nprim)
                 prim2(i,1-j,1:nvar2) = prim2(i,1,1:nvar2)
                enddo
               enddo
                        else if(ybc_low.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,1-j,1:nprim) = prim(i,j,1:nprim)
                 prim2(i,1-j,1:nvar2) = prim2(i,j,1:nvar2)

                 ! v velocity should be opposite
                 prim(i,1-j,3) = -prim(i,j,3)
                 prim2(i,1-j,3) = -prim2(i,j,3)
                enddo
               enddo
                        endif
                      

                        if(ybc_high.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,jmax-1+j,1:nprim) = prim(i,j,1:nprim)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,j,1:nvar2)
                enddo
               enddo
                        else if(ybc_high.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-1,1:nprim)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-1,1:nvar2)
                enddo
               enddo
                        else if(ybc_high.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-j,1:nprim)
                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-j,1:nvar2)

                 ! v velocity should be opposite
                 prim(i,jmax-1+j,3) = -prim(i,jmax-j,3)
                 prim2(i,jmax-1+j,3) = -prim2(i,jmax-j,3)
                enddo
               enddo
                        endif

                  !---------------------------------------
                  !---------------------------------------


             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               cons(i,j,1) = prim(i,j,1)
               cons(i,j,2) = prim(i,j,1)*prim(i,j,2)
               cons(i,j,3) = prim(i,j,1)*prim(i,j,3)
               cons(i,j,4) = prim(i,j,1)*prim(i,j,4) + prim(i,j,1)*(prim(i,j,2)*prim(i,j,2) &
               + prim(i,j,3)*prim(i,j,3))/2.0d0
                 
                 do ns = 1, nspec 
               cons(i,j,4+ns) = prim(i,j,1)*prim(i,j,5+ns)
                 enddo


               cons2(i,j,1) = prim2(i,j,1)
               cons2(i,j,2) = prim2(i,j,1)*prim2(i,j,2)
               cons2(i,j,3) = prim2(i,j,1)*prim2(i,j,3)
               cons2(i,j,4) = prim2(i,j,1)*prim2(i,j,4) + prim2(i,j,1)*(prim2(i,j,2)*prim2(i,j,2) &
               + prim2(i,j,3)*prim2(i,j,3))/2.0d0

               cons2(i,j,5) = prim2(i,j,1)*prim2(i,j,5)
               cons2(i,j,6) = prim2(i,j,1)*prim2(i,j,6)

              enddo
             enddo
             

                       if(irestart.eq.0) then
              print*, 'writing vtk file for paraview '
              ! write paraview vtk file
              title = 'vtk_initial'       
              n = 0
              write(filename,'("vtk/output_",I5.5,".vtk")'),n   



             do i = 1, imax
              xyz(i,:,1) = x(i) + dx 
             enddo
             do j = 1, jmax
              xyz(:,j,2) = y(j) + dy
             enddo
             xyz(:,:,3) = 0.0d0 ! no z-direction

             augmented_prim(:,:,1:nprim) = prim(:,:,1:nprim) 
             augmented_prim(:,:,nprim+1:nprim+nvar2) = prim2(:,:,1:nvar2) 

              call vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim+nvar2    & 
                    ,augmented_prim(1:imax,1:jmax,1:nprim+nvar2))

              n = 0
              time = 0.0d0
              write(filename,'("output/output_",I5.5,".dat")'),n
              open(34,file=filename,form='formatted') 
               write(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 
                 call compute_temperature(prim(i,j,1), prim(i,j,5), prim(i,j,6:5+nspec), T1)
                 write(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4)/cvsolid, prim2(i,j,5), prim2(i,j,6) 
                enddo 
               enddo 
              close(34)

                      endif

!------------------------------------------------------------
              print*, 'beginning time loop '
 
              ! begin time loop
              do n = irestart, nmax
               
               vel = 0.0d0
               do i = 1, imax-1
                do j = 1, jmax-1
                       if(prim(i,j,5).le.0.0d0.or.prim(i,j,1).le.0.0d0) then
                         print*, 'something wrong: p, rho: ', prim(i,j,5), prim(i,j,1), i, j
                         stop
                       endif
                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc) 
                 vel = max(vel,sqrt(prim(i,j,2)**2.0d0 +prim(i,j,3)**2.0d0)+cc) 
                 vel = max(vel,sqrt(prim2(i,j,2)**2.0d0 + prim2(i,j,3)**2.0d0))
                enddo 
               enddo 
                dt = cfl*min(dx,dy)/vel
                time = time + dt
                print*, ' '
                print*, 'iteration, dt, time: ', n, dt, time
                call flush() 


                ! compute fluxes at cell interfaces
                ! here, flux(i) refers to flux(i-1/2) 

                k_rk = 0.0d0                  
                cons0 = cons ! cons0 is the conservative variables (gas) at time step n (beginning of RK)   
 
                k_rk2 = 0.0d0                  
                cons02 = cons2 ! cons02 is the conservative variables (solid) at time step n (beginning of RK)   

         do rkstep = 1, rk

                !-----------------------------------------------------
                ! apply boundary conditions 

                        if(xbc_low.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(1-i,j,1:nprim) = prim(imax-i,j,1:nprim)
                 prim2(1-i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)

                 cons(1-i,j,1:ncons) = cons(imax-i,j,1:ncons)          
                 cons2(1-i,j,1:nvar2) = cons2(imax-i,j,1:nvar2)        

                 k_rk(1-i,j,1:ncons,:) = k_rk(imax-i,j,1:ncons,:)       
                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(imax-i,j,1:nvar2,:)    
                enddo
               enddo
                        else if(xbc_low.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(1-i,j,1:nprim) = prim(1,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(1,j,1:nvar2)

                 cons(1-i,j,1:ncons) = cons(1,j,1:ncons)          

                 cons2(1-i,j,1:nvar2) = cons2(1,j,1:nvar2)          

                 k_rk(1-i,j,1:ncons,:) = k_rk(1,j,1:ncons,:)       

                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(1,j,1:nvar2,:)       
                enddo
               enddo
                        else if(xbc_low.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(1-i,j,1:nprim) = prim(i,j,1:nprim)

                 prim2(1-i,j,1:nvar2) = prim2(i,j,1:nvar2)

                 cons(1-i,j,1:ncons) = cons(i,j,1:ncons)          

                 cons2(1-i,j,1:nvar2) = cons2(i,j,1:nvar2)          

                 k_rk(1-i,j,1:ncons,:) = k_rk(i,j,1:ncons,:)       

                 k_rk2(1-i,j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)       

                 ! u velocity should be opposite
                 prim(1-i,j,2) = -prim(i,j,2)

                 prim2(1-i,j,2) = -prim2(i,j,2)

                 cons(1-i,j,2) = -cons(i,j,2)          

                 cons2(1-i,j,2) = -cons2(i,j,2)          

                 k_rk(1-i,j,2,:) = -k_rk(i,j,2,:)       

                 k_rk2(1-i,j,2,:) = -k_rk2(i,j,2,:)       

                enddo
               enddo
                        endif


                        if(xbc_high.eq.1) then
                 ! BC in x-direction (periodic)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(imax-1+i,j,1:nprim) = prim(i,j,1:nprim)        
                 prim2(imax-1+i,j,1:nvar2) = prim2(i,j,1:nvar2)        

                 cons(imax-1+i,j,1:ncons) = cons(i,j,1:ncons)        
                 cons2(imax-1+i,j,1:nvar2) = cons2(i,j,1:nvar2)        

                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(i,j,1:ncons,:)    
                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)    
                enddo
               enddo
                        else if(xbc_high.eq.2) then
                 ! BC in x-direction (outflow)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost 
                 prim(imax-1+i,j,1:nprim) = prim(imax-1,j,1:nprim)     
                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-1,j,1:nvar2)   

                 cons(imax-1+i,j,1:ncons) = cons(imax-1,j,1:ncons)     
                 cons2(imax-1+i,j,1:nvar2) = cons2(imax-1,j,1:nvar2)   

                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(imax-1,j,1:ncons,:) 
                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(imax-1,j,1:nvar2,:)
                enddo
               enddo
                        else if(xbc_high.eq.3) then
                 ! BC in x-direction (slip wall)
               do i = 1, nghost
                do j = 1-nghost,jmax-1+nghost
                 prim(imax-1+i,j,1:nprim) = prim(imax-i,j,1:nprim)

                 prim2(imax-1+i,j,1:nvar2) = prim2(imax-i,j,1:nvar2)

                 cons(imax-1+i,j,1:ncons) = cons(imax-i,j,1:ncons)     

                 cons2(imax-1+i,j,1:nvar2) = cons2(imax-i,j,1:nvar2)   

                 k_rk(imax-1+i,j,1:ncons,:) = k_rk(imax-i,j,1:ncons,:) 

                 k_rk2(imax-1+i,j,1:nvar2,:) = k_rk2(imax-i,j,1:nvar2,:)

                 ! u velocity should be opposite
                 prim(imax-1+i,j,2) = -prim(imax-i,j,2)

                 prim2(imax-1+i,j,2) = -prim2(imax-i,j,2)

                 cons(imax-1+i,j,2) = -cons(imax-i,j,2)     

                 cons2(imax-1+i,j,2) = -cons2(imax-i,j,2)     

                 k_rk(imax-1+i,j,2,:) = -k_rk(imax-i,j,2,:) 

                 k_rk2(imax-1+i,j,2,:) = -k_rk2(imax-i,j,2,:) 

                enddo
               enddo
                        endif

                !-----------------------------------------------------

                        if(ybc_low.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,jmax-j,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,jmax-j,1:nvar2)

                 cons(i,1-j,1:ncons) = cons(i,jmax-j,1:ncons)          

                 cons2(i,1-j,1:nvar2) = cons2(i,jmax-j,1:nvar2)        

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,jmax-j,1:ncons,:)       

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,jmax-j,1:nvar2,:)     
                enddo
               enddo
                        else if(ybc_low.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,1,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,1,1:nvar2)

                 cons(i,1-j,1:ncons) = cons(i,1,1:ncons)          

                 cons2(i,1-j,1:nvar2) = cons2(i,1,1:nvar2)          

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,1,1:ncons,:)       

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,1,1:nvar2,:)       
                enddo
               enddo
                        else if(ybc_low.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,1-j,1:nprim) = prim(i,j,1:nprim)

                 prim2(i,1-j,1:nvar2) = prim2(i,j,1:nvar2)

                 cons(i,1-j,1:ncons) = cons(i,j,1:ncons)          

                 cons2(i,1-j,1:nvar2) = cons2(i,j,1:nvar2)          

                 k_rk(i,1-j,1:ncons,:) = k_rk(i,j,1:ncons,:)       

                 k_rk2(i,1-j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)       

                 ! v velocity should be opposite
                 prim(i,1-j,3) = -prim(i,j,3)

                 prim2(i,1-j,3) = -prim2(i,j,3)

                 cons(i,1-j,3) = -cons(i,j,3)          

                 cons2(i,1-j,3) = -cons2(i,j,3)          

                 k_rk(i,1-j,3,:) = -k_rk(i,j,3,:)       

                 k_rk2(i,1-j,3,:) = -k_rk2(i,j,3,:)       

                enddo
               enddo
                        endif


                        if(ybc_high.eq.1) then
                 ! BC in y-direction (periodic)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,jmax-1+j,1:nprim) = prim(i,j,1:nprim)        

                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,j,1:nvar2)        

                 cons(i,jmax-1+j,1:ncons) = cons(i,j,1:ncons)        

                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,j,1:nvar2)        

                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,j,1:ncons,:)    

                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,j,1:nvar2,:)    
                enddo
               enddo
                        else if(ybc_high.eq.2) then
                 ! BC in y-direction (outflow)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-1,1:nprim)     

                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-1,1:nvar2)    

                 cons(i,jmax-1+j,1:ncons) = cons(i,jmax-1,1:ncons)     

                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,jmax-1,1:nvar2)   

                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,jmax-1,1:ncons,:)  

                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,jmax-1,1:nvar2,:)
                enddo
               enddo
                        else if(ybc_high.eq.3) then
                 ! BC in y-direction (slip wall)
               do j = 1, nghost
                do i = 1-nghost,imax-1+nghost 
                 prim(i,jmax-1+j,1:nprim) = prim(i,jmax-j,1:nprim)     

                 prim2(i,jmax-1+j,1:nvar2) = prim2(i,jmax-j,1:nvar2)     

                 cons(i,jmax-1+j,1:ncons) = cons(i,jmax-j,1:ncons)     

                 cons2(i,jmax-1+j,1:nvar2) = cons2(i,jmax-j,1:nvar2)     

                 k_rk(i,jmax-1+j,1:ncons,:) = k_rk(i,jmax-j,1:ncons,:)  

                 k_rk2(i,jmax-1+j,1:nvar2,:) = k_rk2(i,jmax-j,1:nvar2,:)  

                 ! v velocity should be opposite
                 prim(i,jmax-1+j,3) = -prim(i,jmax-j,3)     

                 prim2(i,jmax-1+j,3) = -prim2(i,jmax-j,3)     

                 cons(i,jmax-1+j,3) = -cons(i,jmax-j,3)     

                 cons2(i,jmax-1+j,3) = -cons2(i,jmax-j,3)     

                 k_rk(i,jmax-1+j,3,:) = -k_rk(i,jmax-j,3,:)  

                 k_rk2(i,jmax-1+j,3,:) = -k_rk2(i,jmax-j,3,:)  

                enddo
               enddo
                        endif


                !-----------------------------------------------------
                ! Interpolate in x-direction 
                do i = -1, imax+2 
                 do j = -1, jmax+2 

                 ! compute primitive and conservative variables at the cell interfaces
                 ! i refers to i-1/2 for interface variables


             if(space_scheme.eq.1) then

              do ns = 1, ncons
               d_im32 = cons(i-1,j,ns) - cons(i-2,j,ns)
               d_im12 = cons(i,j,ns) - cons(i-1,j,ns)
               d_ip12 = cons(i+1,j,ns) - cons(i,j,ns)

               cons_intxl(i,j,ns) = cons(i-1,j,ns)   &
               + eps1/4.0d0*((1.0d0-kappa)*d_im32+(1.0d0+kappa)*d_im12)
               cons_intxr(i,j,ns) = cons(i,j,ns)     &
               - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo 


              do ns = 1, nvar2
               d_im32 = cons2(i-1,j,ns) - cons2(i-2,j,ns)
               d_im12 = cons2(i,j,ns) - cons2(i-1,j,ns)
               d_ip12 = cons2(i+1,j,ns) - cons2(i,j,ns)

               cons_intxl2(i,j,ns) = cons2(i-1,j,ns)   &
               + eps1/4.0d0*((1.0d0-kappa)*d_im32+(1.0d0+kappa)*d_im12)
               cons_intxr2(i,j,ns) = cons2(i,j,ns)     &
               - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo 

             else if(space_scheme.eq.2) then
                  
               cons_intxl(i,j,:) = cons(i-1,j,:)
               cons_intxr(i,j,:) = cons(i,j,:)

               cons_intxl2(i,j,:) = cons2(i-1,j,:)
               cons_intxr2(i,j,:) = cons2(i,j,:)

             else if(space_scheme.eq.3) then

             else if(space_scheme.eq.4) then

              ! WENO5 USES PRIMITIVE VARIABLES
 
              do ns = 1, nprim
                call weno5(prim(i-2:i+2,j,ns),prim_intxl(i+1,j,ns),prim_intxr(i,j,ns))
              enddo

         call compute_pres(prim_intxl(i+1,j,1), prim_intxl(i+1,j,6:5+nspec), prim_intxl(i+1,j,4), prim_intxl(i+1,j,5)) 

                       if(prim_intxl(i+1,j,5).le.0.0d0.or.prim_intxl(i+1,j,1).le.0.0d0) then
                         print*, 'something wrong: p, rho (intxl): ', prim_intxl(i+1,j,5), prim_intxl(i+1,j,1), i, j
                         stop
                       endif

         call compute_pres(prim_intxr(i,j,1), prim_intxr(i,j,6:5+nspec), prim_intxr(i,j,4), prim_intxr(i,j,5)) 

                       if(prim_intxr(i,j,5).le.0.0d0.or.prim_intxr(i,j,1).le.0.0d0) then
                         print*, 'something wrong: p, rho (intxr): ', prim_intxr(i,j,5), prim_intxr(i,j,1), i, j
                         stop
                       endif


              prim_intxl(i+1,j,5+nspec) = min(max(1.0d0 - sum(prim_intxl(i+1,j,6:5+nspec-1)),0.0d0),1.0d0)
              prim_intxr(i,j,5+nspec) = min(max(1.0d0 - sum(prim_intxr(i,j,6:5+nspec-1)),0.0d0),1.0d0)


              do ns = 1, nvar2
                call weno5(prim2(i-2:i+2,j,ns),prim_intxl2(i+1,j,ns),prim_intxr2(i,j,ns))
              enddo


             else
               
               print*, 'wrong entry for space_scheme ', space_scheme

             endif

                enddo 
               enddo 





                       if(space_scheme.ne.4) then

                do i = 0, imax+1 
                 do j = 0, jmax+1 


                          ! left
                 call primitive_gas(cons_intxl(i,j,1:ncons),prim_intxl(i,j,1:nprim),rhocut,nspec)  
                 call primitive_solid(cons_intxl2(i,j,1:nvar2),prim_intxl2(i,j,1:nvar2),nvar2,rho2cut)


                          ! right
                 call primitive_gas(cons_intxr(i,j,1:ncons),prim_intxr(i,j,1:nprim),rhocut,nspec)  
                 call primitive_solid(cons_intxr2(i,j,1:nvar2),prim_intxr2(i,j,1:nvar2),nvar2,rho2cut)


                 enddo
                enddo  

                        endif

                !-----------------------------------------------------
                ! Interpolate in y-direction 
                do j = -1, jmax+2 
                 do i = -1, imax+2 

                 ! compute primitive and conservative variables at the cell interfaces
                 ! i refers to i-1/2 for interface variables


                  if(space_scheme.eq.1) then

              do ns = 1, ncons
               d_im32 = cons(i,j-1,ns) - cons(i,j-2,ns)
               d_im12 = cons(i,j,ns) - cons(i,j-1,ns)
               d_ip12 = cons(i,j+1,ns) - cons(i,j,ns)

               cons_intyl(i,j,ns) = cons(i,j-1,ns)   &
                + eps1/4.0d0*((1.0d0-kappa)*d_im32+(1.0d0+kappa)*d_im12)
               cons_intyr(i,j,ns) = cons(i,j,ns)     &
                - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

              do ns = 1, nvar2
               d_im32 = cons2(i,j-1,ns) - cons2(i,j-2,ns)
               d_im12 = cons2(i,j,ns) - cons2(i,j-1,ns)
               d_ip12 = cons2(i,j+1,ns) - cons2(i,j,ns)

               cons_intyl2(i,j,ns) = cons2(i,j-1,ns)   &
                + eps1/4.0d0*((1.0d0-kappa)*d_im32+(1.0d0+kappa)*d_im12)
               cons_intyr2(i,j,ns) = cons2(i,j,ns)     &
                - eps1/4.0d0*((1.0d0+kappa)*d_im12+(1.0d0-kappa)*d_ip12)
              enddo

             else if(space_scheme.eq.2) then
                  
               cons_intyl(i,j,:) = cons(i,j-1,:)
               cons_intyr(i,j,:) = cons(i,j,:)

               cons_intyl2(i,j,:) = cons2(i,j-1,:)
               cons_intyr2(i,j,:) = cons2(i,j,:)

             else if(space_scheme.eq.3) then

             else if(space_scheme.eq.4) then

              ! WENO5 USES PRIMITIVE VARIABLES
 
              do ns = 1, nprim
                call weno5(prim(i,j-2:j+2,ns),prim_intyl(i,j+1,ns),prim_intyr(i,j,ns))
              enddo

              do ns = 1, nvar2
                call weno5(prim2(i,j-2:j+2,ns),prim_intyl2(i,j+1,ns),prim_intyr2(i,j,ns))
              enddo

         call compute_pres(prim_intyl(i,j+1,1), prim_intyl(i,j+1,6:5+nspec), prim_intyl(i,j+1,4), prim_intyl(i,j+1,5)) 

                       if(prim_intyl(i,j+1,5).le.0.0d0.or.prim_intyl(i,j+1,1).le.0.0d0) then
                         print*, 'something wrong: p, rho (intyl): ', prim_intyl(i,j+1,5), prim_intyl(i,j+1,1), i, j
                         stop
                       endif

         call compute_pres(prim_intyr(i,j,1), prim_intyr(i,j,6:5+nspec), prim_intyr(i,j,4), prim_intyr(i,j,5)) 

                       if(prim_intyr(i,j,5).le.0.0d0.or.prim_intyr(i,j,1).le.0.0d0) then
                         print*, 'something wrong: p, rho (intyr): ', prim_intyr(i,j,5), prim_intyr(i,j,1), i, j
                         stop
                       endif

              prim_intyl(i,j+1,5+nspec) = min(max(1.0d0 - sum(prim_intyl(i,j+1,6:5+nspec-1)),0.0d0),1.0d0)
              prim_intyr(i,j,5+nspec) = min(max(1.0d0 - sum(prim_intyr(i,j,6:5+nspec-1)),0.0d0),1.0d0)


             else

               print*, 'wrong entry for space_scheme ', space_scheme

             endif

                 enddo
                enddo  


                     
                       if(space_scheme.ne.4) then

                do i = 0, imax+1 
                 do j = 0, jmax+1 

                          ! left
                 call primitive_gas(cons_intyl(i,j,1:ncons),prim_intyl(i,j,1:nprim),rhocut,nspec)  
                 call primitive_solid(cons_intyl2(i,j,1:nvar2),prim_intyl2(i,j,1:nvar2),nvar2,rho2cut)


                          ! right
                 call primitive_gas(cons_intyr(i,j,1:ncons),prim_intyr(i,j,1:nprim),rhocut,nspec)  
                 call primitive_solid(cons_intyr2(i,j,1:nvar2),prim_intyr2(i,j,1:nvar2),nvar2,rho2cut)

                 enddo
                enddo  

                     endif

                !-----------------------------------------------------
          
 
           xflux = 0.0d0
           yflux = 0.0d0
           xflux2 = 0.0d0
           yflux2 = 0.0d0

 
                 ! compute x-fluxes
           do i = 1, imax !0, imax
            do j = 1, jmax !0, jmax


                 rhol = prim_intxl(i,j,1)
                 rhor = prim_intxr(i,j,1)
                 pl = prim_intxl(i,j,5)
                 pr = prim_intxr(i,j,5)
                 ul = prim_intxl(i,j,2) 
                 ur = prim_intxr(i,j,2)
                 vl = prim_intxl(i,j,3)
                 vr = prim_intxr(i,j,3)
                 iel = prim_intxl(i,j,4)
                 ier = prim_intxr(i,j,4)
                 do ns = 1, nspec
                  ysl(ns) = prim_intxl(i,j,5+ns) 
                  ysr(ns) = prim_intxr(i,j,5+ns) 
                 enddo
 
                 call compute_sound(pl,rhol,ysl,cc) 
                 gaml = cc*cc*rhol/pl 
                 call compute_sound(pr,rhor,ysr,cc) 
                 gamr = cc*cc*rhor/pr 
                 

                   if(iel.le.0.0d0) then
                     print*, 'iel ', iel
                     stop
                   endif 
                   if(ier.le.0.0d0) then
                     print*, 'ier ', ier
                     stop
                   endif 

                 sweep = 1
                call hllc_2d(gaml, gamr, rhol, rhor, pl, pr, ul, ur, vl, vr, iel,   &
                   ier, ysl, ysr, sweep, flux)

                 xflux(i,j,1:ncons) = flux(1:ncons)



                            ! SOLID 

              spl = abs(prim_intxl2(i,j,2))  
              spr = abs(prim_intxr2(i,j,2))  
              wavespeed = max(spl,spr)
         
                 ! continuity
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)
                 termrl = prim_intxr2(i,j,1) - prim_intxl2(i,j,1)
                 xflux2(i,j,1) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! x-momentum
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,2)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,2)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,2) - prim_intxl2(i,j,1)*prim_intxl2(i,j,2)
                 xflux2(i,j,2) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! y-momentum
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,3)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,3)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,3) - prim_intxl2(i,j,1)*prim_intxl2(i,j,3)
                 xflux2(i,j,3) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! energy
                 ke = 0.5d0*prim_intxl2(i,j,1)*(prim_intxl2(i,j,2)**2.0d0 + prim_intxl2(i,j,3)**2.0d0)
                 etotl = prim_intxl2(i,j,1)*prim_intxl2(i,j,4) + ke
                 terml = etotl*prim_intxl2(i,j,2)
                 ke = 0.5d0*prim_intxr2(i,j,1)*(prim_intxr2(i,j,2)**2.0d0 + prim_intxr2(i,j,3)**2.0d0)
                 etotr = prim_intxr2(i,j,1)*prim_intxr2(i,j,4) + ke
                 termr = etotr*prim_intxr2(i,j,2)
                 termrl = etotr - etotl
                 xflux2(i,j,4) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 

                 ! fign
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,5)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,5)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,5) - prim_intxl2(i,j,1)*prim_intxl2(i,j,5)
                 xflux2(i,j,5) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! Num
                 terml = prim_intxl2(i,j,1)*prim_intxl2(i,j,2)*prim_intxl2(i,j,6)  
                 termr = prim_intxr2(i,j,1)*prim_intxr2(i,j,2)*prim_intxr2(i,j,6)  
                 termrl = prim_intxr2(i,j,1)*prim_intxr2(i,j,6) - prim_intxl2(i,j,1)*prim_intxl2(i,j,6)
                 xflux2(i,j,6) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl


            enddo
           enddo

                          !--------!


                 ! compute y-fluxes
           do i = 1, imax !0, imax
            do j = 1, jmax !0, jmax

 
                 rhol = prim_intyl(i,j,1)
                 rhor = prim_intyr(i,j,1)
                 pl = prim_intyl(i,j,5)
                 pr = prim_intyr(i,j,5)
                 ul = prim_intyl(i,j,2)
                 ur = prim_intyr(i,j,2)
                 vl = prim_intyl(i,j,3)
                 vr = prim_intyr(i,j,3)
                 iel = prim_intyl(i,j,4)
                 ier = prim_intyr(i,j,4)
                 do ns = 1, nspec
                  ysl(ns) = prim_intyl(i,j,5+ns)
                  ysr(ns) = prim_intyr(i,j,5+ns)
                 enddo

                 call compute_sound(pl,rhol,ysl,cc)
                 gaml = cc*cc*rhol/pl
                 call compute_sound(pr,rhor,ysr,cc)
                 gamr = cc*cc*rhor/pr

                   if(iel.le.0.0d0) then
                     print*, 'iel ', iel
                     stop
                   endif 
                   if(ier.le.0.0d0) then
                     print*, 'ier ', ier
                     stop
                   endif 

                 sweep = 2
                call hllc_2d(gaml, gamr, rhol, rhor, pl, pr, ul, ur, vl, vr, iel,   &
                   ier, ysl, ysr, sweep, flux)

                 yflux(i,j,1:ncons) = flux(1:ncons)



                            ! SOLID


              spl = abs(prim_intyl2(i,j,3))
              spr = abs(prim_intyr2(i,j,3)) 
              wavespeed = max(spl,spr)


                 ! continuity
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)
                 termrl = prim_intyr2(i,j,1) - prim_intyl2(i,j,1)
                 yflux2(i,j,1) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! x-momentum
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,2)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,2)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,2) - prim_intyl2(i,j,1)*prim_intyl2(i,j,2)
                 yflux2(i,j,2) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! y-momentum
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,3)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,3)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,3) - prim_intyl2(i,j,1)*prim_intyl2(i,j,3)
                 yflux2(i,j,3) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! energy
                 ke = 0.5d0*prim_intyl2(i,j,1)*(prim_intyl2(i,j,2)**2.0d0 + prim_intyl2(i,j,3)**2.0d0)
                 etotl = prim_intyl2(i,j,1)*prim_intyl2(i,j,4) + ke
                 terml = etotl*prim_intyl2(i,j,3)
                 ke = 0.5d0*prim_intyr2(i,j,1)*(prim_intyr2(i,j,2)**2.0d0 + prim_intyr2(i,j,3)**2.0d0)
                 etotr = prim_intyr2(i,j,1)*prim_intyr2(i,j,4) + ke
                 termr = etotr*prim_intyr2(i,j,3)
                 termrl = etotr - etotl
                 yflux2(i,j,4) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl

                 ! fign
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,5)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,5)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,5) - prim_intyl2(i,j,1)*prim_intyl2(i,j,5)
                 yflux2(i,j,5) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 
                
                 ! Num
                 terml = prim_intyl2(i,j,1)*prim_intyl2(i,j,3)*prim_intyl2(i,j,6)  
                 termr = prim_intyr2(i,j,1)*prim_intyr2(i,j,3)*prim_intyr2(i,j,6)  
                 termrl = prim_intyr2(i,j,1)*prim_intyr2(i,j,6) - prim_intyl2(i,j,1)*prim_intyl2(i,j,6)
                 yflux2(i,j,6) = 0.5d0*(terml+termr) - 0.5d0*wavespeed*termrl 
                

            enddo
           enddo


                !-----------------------------------------------------



                ! Runge-Kutta increment
                do i = 1, imax-1 
                 do j = 1, jmax-1 
             k_rk(i,j,1:ncons,rkstep) = - (xflux(i+1,j,:)-xflux(i,j,:))*dt/dx   &
                     - (yflux(i,j+1,:)-yflux(i,j,:))*dt/dy
             k_rk2(i,j,1:nvar2,rkstep) = - (xflux2(i+1,j,:)-xflux2(i,j,:))*dt/dx   &
                     - (yflux2(i,j+1,:)-yflux2(i,j,:))*dt/dy
                 enddo
                enddo


                ! update the conserved variables
                ! assume no source terms for now
            do i = 1, imax-1 
             do j = 1, jmax-1 
              if(rkstep.eq.1) then
                   cons(i,j,:) = cons0(i,j,:) + k_rk(i,j,:,1)/2.0d0
                   cons2(i,j,:) = cons02(i,j,:) + k_rk2(i,j,:,1)/2.0d0
              endif 
              if(rkstep.eq.2) then
                   cons(i,j,:) = cons0(i,j,:) - k_rk(i,j,:,1) + 2.0d0*k_rk(i,j,:,2)
                   cons2(i,j,:) = cons02(i,j,:) - k_rk2(i,j,:,1) + 2.0d0*k_rk2(i,j,:,2)
              endif 
             enddo
            enddo

               if(rkstep.eq.rk) then
                do i = 1, imax-1
                 do j = 1, jmax-1
         cons(i,j,:) = cons0(i,j,:) + 1.0d0/6.0d0*(1.0d0*k_rk(i,j,:,1)    &
         + 4.0d0*k_rk(i,j,:,2) + 1.0d0*k_rk(i,j,:,3))
         cons2(i,j,:) = cons02(i,j,:) + 1.0d0/6.0d0*(1.0d0*k_rk2(i,j,:,1)    &
         + 4.0d0*k_rk2(i,j,:,2) + 1.0d0*k_rk2(i,j,:,3))
                 enddo
                enddo
               endif

                
                ! recompute primitive variables
                do i = 1, imax-1
                 do j = 1, jmax-1

                  call primitive_gas(cons(i,j,1:ncons),prim(i,j,1:nprim),rhocut,nspec) 
                  call primitive_solid(cons2(i,j,1:nvar2),prim2(i,j,1:nvar2),nvar2,rho2cut)

                  prim(i,j,5+nspec) = min(max(1.0d0 - sum(prim(i,j,6:5+nspec-1)),0.0d0),1.0d0)

                 enddo ! do j             
                enddo ! do i             

 
         ! close rk loop
         enddo   


               !----------------------------------------------------------------
               ! add source terms (gravity)
             
              
             !do i = 1, imax-1 !1-nghost, imax-1+nghost
             ! do j = 1, jmax-1 !1-nghost, jmax-1+nghost
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               source(1:ncons) = 0.0d0
               source(3) = gravity*dt*cons(i,j,1)
               source(4) = gravity*dt*cons(i,j,3)

               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)
              enddo
             enddo


                       ! chemistry

             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
               source(1:ncons) = 0.0d0

                  if(cons(i,j,5).gt.0.0d0.and.cons(i,j,6).gt.0.0d0) then
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 

                     temp = prim(i,j,5)/Rgas/prim(i,j,1)

              do ns = 1, nspec
               !conc(ns) = (cons(i,j,4+ns)/1000.0d0)/MW(ns)
               conc(ns) = cons(i,j,4+ns)
              enddo


              omegadot = 0.0d0
                  ! DANGER
              !call combustion(temp,conc,omegadot)
              !omegadot = omegadot/dt

              do ns = 1, nspec
               source(4+ns) = omegadot(ns)
              enddo

                  endif

               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)*dt

                     do ns = 1, nspec
                      cons(i,j,4+ns) = max(cons(i,j,4+ns),0.0d0)
                     enddo


              enddo
             enddo
              
 
 
                        ! TWO-PHASE SOURCE TERMS
             do i = 1, imax-1
              do j = 1, jmax-1
                source(1:ncons) = 0.0d0
                source_fign = 0.0d0
                hfs = 0.0d0
                sigmadot = 0.0d0

                if(prim2(i,j,1).gt.rho2cut) then

               ! find gas temperature 
               molmix = 0.0d0
               do ns = 1, nspec
                molmix = molmix + prim(i,j,5+ns)/molwt(ns)
               enddo 
               molmix = 1.0d0/molmix
               Rgas = Runiv/molmix 
               T1 = prim(i,j,5)/Rgas/prim(i,j,1)


                 ! find particle diameter
                 num = prim2(i,j,6)
                 dia1 = (6.0d0/3.14159265d0*prim2(i,j,1)/rhosol/num)**(1.0d0/3.0d0)  
                 dia1 = max(min(dia1,dia),0.25d-6) 

                 ! find reynolds number
                 vel = sqrt((prim(i,j,2)-prim2(i,j,2))**2.0d0 + (prim(i,j,3)-prim2(i,j,3))**2.0d0)
                 mu = 4.0d-5*sqrt(T1/300.0d0) 
                 reyn = dia1*vel*prim(i,j,1)/mu   

                 ! find Mach number
                 call compute_sound(prim(i,j,5),prim(i,j,1),prim(i,j,6:5+nspec),cc) 
                 mach = vel/cc
              

                 ! find drag coeff.      
                 call drag(reyn,mach,cd)

 
                 ! find nusselt number 
                 nuss = 2.0d0 + 0.6d0*(prandtl**0.333d0)*sqrt(reyn)    


                 ! find cpmix
                 call COMPUTE_CPK (T1, CPK)
                 cpmix = 0.0d0
                 do ns = 1, nspeci
                  cpmix = cpmix + prim(i,j,5+ns)*cpk(ns)
                 enddo

                 ! find lambda
                 lambda = mu*cpmix/prandtl

                 ! T2
                 T2 = max(prim2(i,j,4)/cvsolid, 0.0d0)
                 T2 = min(T2,6000.0d0) ! DANGER: TO ENSURE SOURCE TERM DOESN'T BLOW UP  

                 ! find fdrag and qheat
                 fdrag(1) = 3.0d0/4.0d0*prim(i,j,1)/rhosol*prim2(i,j,1)/dia1*cd*vel*(prim(i,j,2)-prim2(i,j,2)) 
                 fdrag(2) = 3.0d0/4.0d0*prim(i,j,1)/rhosol*prim2(i,j,1)/dia1*cd*vel*(prim(i,j,3)-prim2(i,j,3)) 
                 qheat = 6.0d0*prim2(i,j,1)/rhosol/dia1*nuss*lambda*(T1-T2)/dia1

                 Tign = 933.0d0
                 ! mass transfer
                 if(T2.ge.Tign.and.prim2(i,j,5).ge.0.99d0.and.prim2(i,j,6).gt.0.0d0) then
                 !if(prim2(i,j,5).ge.0.99d0) then
                  Jkin = 3.1415925d0*dia1*dia1*prim2(i,j,6)*Zhyb*exp(-Ea_R/T2)
                  taup = 4.0d6*(dia*dia) !/(prim(i,j,7)**0.9d0) 
                  Jdiff = 3.0d0*prim2(i,j,1)*(1.0d0 + 0.276d0*sqrt(reyn))/taup
                  muc = 1.0d0/(1.0d0 + exp((130.0d-3-prim2(i,j,1))/(20.0d-3)))

                  sigmadot = 1.0d0/(1.0d0/Jkin + 1.0d0/Jdiff)*muc

                 else
                  sigmadot = 0.0d0  
                 endif             


                  if(sigmadot.lt.0.0d0) then
                   print*, 'bug in sigmadot: ', sigmadot, Jkin, Jdiff
                   call flush()
                   stop
                  endif

                  sigmadot = min(sigmadot,cons2(i,j,1)/dt)


                  call HTOTAL(T2,hfs)

                  ! DANGER
                  sigmadot = 0.0d0


                  source(1) = sigmadot
                  source(2) = sigmadot*prim2(i,j,2) - fdrag(1)
                  source(3) = sigmadot*prim2(i,j,3) - fdrag(2)
                  source(4) = sigmadot*(prim2(i,j,4) + 0.5d0*(prim2(i,j,2)**2.0d0 + prim2(i,j,3)**2.0d0)) &
                               - fdrag(1)*prim2(i,j,2) - fdrag(2)*prim2(i,j,3) - qheat 
                  source(5) = sigmadot
                  source(6:4+nspec) = 0.0d0

                  ! KUHL-BOIKO CURVE FIT
                  tign = 1.0d0/(6.25d10)*exp(30000/T1) 
                  source_fign = prim2(i,j,1)/tign

                else
                 source = 0.0d0 
                 source_fign = 0.0d0
                endif



               cons(i,j,1:ncons) = cons(i,j,1:ncons) + source(1:ncons)*dt
               cons(i,j,4) = cons(i,j,4) + sigmadot*dt*hfs ! ENERGY FROM MASS TRANSFER

               cons2(i,j,1:4) = cons2(i,j,1:4) - source(1:4)*dt  
               cons2(i,j,5) = cons2(i,j,5) + source_fign*dt 

               if(cons2(i,j,1).le.rho2cut) then
                 prim2(i,j,1:nvar2) = 0.0d0 
                 cons2(i,j,1:nvar2) = 0.0d0
               endif

              enddo  
             enddo

               !----------------------------------------------------------------

                    ! upper bound on fign
                 
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost
                cons2(i,j,5) = min(cons2(i,j,5),1.5d0*cons2(i,j,1))
              enddo
             enddo

               !----------------------------------------------------------------

                    ! recompute primitive variables
             do i = 1-nghost, imax-1+nghost
              do j = 1-nghost, jmax-1+nghost

                  call primitive_gas(cons(i,j,1:ncons),prim(i,j,1:nprim),rhocut,nspec) 
                  call primitive_solid(cons2(i,j,1:nvar2),prim2(i,j,1:nvar2),nvar2,rho2cut)

                  prim(i,j,5+nspec) = min(max(1.0d0 - sum(prim(i,j,6:5+nspec-1)),0.0d0),1.0d0)

             enddo ! do j             
            enddo ! do i

               !----------------------------------------------------------------
                  maxvel = 0.0d0
                  maxp = 0.0d0
                  maxrho = 0.0d0
                  minvel = 1.0d10
                  minp = 1.0d10
                  minrho = 1.0d10

                  do i = 1, imax-1 
                   do j = 1, jmax-1
                      maxp = max(maxp,prim(i,j,5))
                      maxrho = max(maxrho,prim(i,j,1))
                      maxvel = max(maxvel,abs(prim(i,j,2)),abs(prim(i,j,3)))

                      minp = min(minp,prim(i,j,5))
                      minrho = min(minrho,prim(i,j,1))
                      minvel = min(minvel,abs(prim(i,j,2)),abs(prim(i,j,3)))
                   enddo
                  enddo  
                 print*, 'max (pressure,density,velocity): ', maxp, maxrho, maxvel 
                 print*, 'min (pressure,density,velocity): ', minp, minrho, minvel 

             augmented_prim(:,:,1:nprim) = prim(:,:,1:nprim) 
             augmented_prim(:,:,nprim+1:nprim+nvar2) = prim2(:,:,1:nvar2) 


                 ! write paraview vtk file 
                 if(mod(n,ofile).eq.0) then

              print*, ' '
              print*, 'writing vtk file for paraview '
              title = 'vtk/vtk_output'       

             do i = 1, imax
              xyz(i,:,1) = x(i) + dx 
             enddo
             do j = 1, jmax
              xyz(:,j,2) = y(j) + dy
             enddo
             xyz(:,:,3) = 0.0d0 ! no z-direction

              write(filename,'("vtk/output_",I5.5,".vtk")'),n   
              call vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim+nvar2    & 
                       ,augmented_prim(1:imax,1:jmax,1:nprim+nvar2))

              write(filename,'("output/output_",I5.5,".dat")'),n
              open(34,file=filename,form='formatted') 
               write(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1 
                 call compute_temperature(prim(i,j,1), prim(i,j,5), prim(i,j,6:5+nspec), T1)
                 write(34,*) prim(i,j,1), prim(i,j,2), prim(i,j,3), T1, prim(i,j,5), prim(i,j,6), & 
                   prim(i,j,7), prim2(i,j,1), prim2(i,j,2), & 
                   prim2(i,j,3), prim2(i,j,4)/cvsolid, prim2(i,j,5), prim2(i,j,6) 
                enddo 
               enddo 
              close(34)

                 endif 


                ! close time loop (do n)
              enddo
              

 
              print*, 'done ' 


                               end program

!--------------------------------------------------------------------------------
