                 subroutine drag(reyn,Mach,cd)

          implicit none 

          real(kind=8) :: reyn, Mach, cd
          real(kind=8) :: cd_std, cd_sub, cd_sup
          real(kind=8) :: cd_M1, cd_M175, cd_M06, zeta, f1, f2, f3
          real(kind=8) :: term1, term2, term3, c1, c2, c3
          real(kind=8) :: c0, delta0
          real(kind=8) :: cd_parmar, cd_abraham


            if(reyn.ge.0.1d0) then


                 ! ABRAHAM

             c0 = 0.2924d0
             delta0 = 9.06d0
             cd_abraham = c0*(1.0d0 + delta0/sqrt(reyn))**2.0d0


                   if(reyn.le.50.0d0) then
                     cd = cd_abraham
                     return
                   endif              
 




             ! CLIFT & GAUVIN - STANDARD DRAG LAW 
             cd_std = 24.0d0/reyn*(1.0d0 + 0.15d0*reyn**0.687d0) + &
                      0.42d0/(1.0d0 + 42500.0d0/reyn**1.16d0)


             ! SUPERSONIC REGIME 

             cd_M1 = 24.0d0/reyn*(1.0d0 + 0.118d0*reyn**0.813d0) + &
                     0.69d0/(1.0d0 + 3550.0d0/reyn**0.793d0)

             cd_M175 = 24.0d0/reyn*(1.0d0 + 0.107d0*reyn**0.867d0) + &
                       0.646d0/(1.0d0 + 861.0d0/reyn**0.634d0)

             f1 = -2.963d0 + 4.392d0*Mach - 1.169d0*Mach**2.0d0 &
               - 0.027d0*Mach**3.0d0 - 0.233d0*exp((1.0d0-Mach)/0.011d0)
             f2 = -6.617d0 + 12.11d0*Mach - 6.501d0*Mach**2.0d0 &
               + 1.182d0*Mach**3.0d0 - 0.174d0*exp((1.0d0-Mach)/0.01d0)
             f3 = -5.866d0 + 11.57d0*Mach - 6.665d0*Mach**2.0d0 &
               + 1.312d0*Mach**3.0d0 - 0.35d0*exp((1.0d0-Mach)/0.012d0)

             c1 = 6.48d0
             c2 = 8.93d0
             c3 = 12.21d0

             term1 = f1*(log(reyn)-c2)/(c1-c2)*(log(reyn)-c3)/(c1-c3)
             term2 = f2*(log(reyn)-c1)/(c2-c1)*(log(reyn)-c3)/(c2-c3)
             term3 = f3*(log(reyn)-c1)/(c3-c1)*(log(reyn)-c2)/(c3-c2)

             zeta = term1 + term2 + term3

             cd_sup = cd_M1 + (cd_M175 - cd_M1)*zeta


             ! INTERMEDIATE REGIME 

             cd_M06 = 24.0d0/reyn*(1.0d0 + 0.15d0*reyn**0.684d0) + &
                     0.5130/(1.0d0 + 483.0d0/reyn**0.669d0)

             cd_M1 = 24.0d0/reyn*(1.0d0 + 0.118d0*reyn**0.813d0) + &
                     0.69d0/(1.0d0 + 3550.0d0/reyn**0.793d0)

             f1 = -1.884d0 + 8.422d0*Mach - 13.7d0*Mach**2.0d0 +   &
                   8.162d0*Mach**3.0d0 
             f2 = -2.228d0 + 10.35d0*Mach - 16.96d0*Mach**2.0d0 +  &
                   9.840d0*Mach**3.0d0
             f3 = 4.362d0 - 16.91d0*Mach + 19.84d0*Mach**2.0d0 -   &
                   6.296d0*Mach**3.0d0 

             c1 = 6.48d0
             c2 = 9.28d0
             c3 = 12.21d0

             term1 = f1*(log(reyn)-c2)/(c1-c2)*(log(reyn)-c3)/(c1-c3)
             term2 = f2*(log(reyn)-c1)/(c2-c1)*(log(reyn)-c3)/(c2-c3)
             term3 = f3*(log(reyn)-c1)/(c3-c1)*(log(reyn)-c2)/(c3-c2)

             zeta = term1 + term2 + term3

             cd_sub = cd_M06 + (cd_M1 - cd_M06)*zeta


              if(Mach.le.0.6d0) then
               cd_parmar = cd_std + (cd_M06 - cd_std)*(Mach/0.6d0)
              else if(Mach.gt.0.6d0.and.Mach.le.1.0d0) then
               cd_parmar = cd_sub
              else
               cd_parmar = cd_sup
              endif


                   if(reyn.ge.100.0d0) then
                     cd = cd_parmar
                     return
                   endif              
 

                  if(reyn.gt.50.0d0.and.reyn.lt.100.0d0) then
                    cd = ((100.0d0-reyn)*cd_abraham + (reyn-50.0d0)*cd_parmar)/(100.0d0 - 50.0d0)
                    return
                  endif



            else 

                cd = 0.0d0

            endif



             return
               end subroutine

!---------------------------------------------------------------------------------  
