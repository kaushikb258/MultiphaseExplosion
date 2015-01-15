
             subroutine combustion(T,conc,omegadot)

          USE THERMO
          real(kind=8) :: T, conc(1:nspeci)
          real(kind=8), dimension(1:nspeci) :: omegadot
          real(kind=8) :: fa
          
           
              omegadot = 0.0d0



            return 
            end subroutine
          
!--------------------------------------------------------------------------------

           SUBROUTINE COMPUTE_CPK (TEMP, CPK)
      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: CPK(1:NSPECI)
      REAL (KIND = 8) :: CPP


           CPK = 719.2D0 + 8314.0D0/(28.9D0)


      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------

      SUBROUTINE COMPUTE_HFK (TEMP, HFK)

      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: HFK(1:NSPECI)
      REAL (KIND = 8) :: HFP, CPP

      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

         HFK = (719.2D0 + 8314.0D0/(28.9D0))*TEMP


      RETURN
      END


!--------------------------------------------------------------------

      SUBROUTINE COMPUTE_SK (TEMP, SK)

      USE THERMO

      IMPLICIT NONE

      INTEGER :: NS
      REAL (KIND = 8), INTENT(IN) :: TEMP
      REAL (KIND = 8), INTENT(OUT) :: SK(1:NSPECI)
      REAL (KIND = 8) :: SP

      REAL(KIND=8), SAVE :: R5I = 1.0D0/5.0D0  
      REAL(KIND=8), SAVE :: R4I = 1.0D0/4.0D0  
      REAL(KIND=8), SAVE :: R3I = 1.0D0/3.0D0  
      REAL(KIND=8), SAVE :: R2I = 1.0D0/2.0D0  

        SK = 0.0D0

      RETURN
      END


!--------------------------------------------------------------------

        SUBROUTINE FROM_IE_TO_T(E_I, YK, TEMP)

      USE THERMO  

      IMPLICIT NONE

      REAL (KIND = 8) :: E_I, TEMP, YK(1:NSPECI)


        TEMP = E_I/(719.2D0)

        RETURN
        END
    
!--------------------------------------------------------------------


        SUBROUTINE FROM_IH_TO_T(H_I, YK, TEMP)

      USE THERMO

      IMPLICIT NONE

      REAL (KIND = 8) :: H_I, TEMP, YK(1:NSPECI)

        TEMP = H_I/(719.2D0 + 8314.0D0/(28.9D0))

        RETURN
        END  

 
!--------------------------------------------------------------------


         subroutine compute_sound(p,rho,Yk,c)

         use thermo

          real(kind=8) :: p, rho, Yk(1:nspeci), c, temp, molwt, Rgas
          real(kind=8) :: cpk(1:nspeci), cvk(1:nspeci), gam
          integer :: ns
          real(kind=8) :: cpmix, cvmix

          molwt = 0.0d0
          do ns = 1, nspeci
           molwt = molwt + Yk(ns)/MW(ns) 
          enddo  
          molwt = 1.0d0/molwt
          Rgas = runiv/molwt
          temp = p/rho/Rgas 


          call COMPUTE_CPK (TEMP, CPK)
          do ns = 1, nspeci
           cvk(ns) = cpk(ns) - runiv/MW(ns)  
          enddo 
          
          cpmix = 0.0d0
          cvmix = 0.0d0
          do ns = 1, nspeci
           cpmix = cpmix + Yk(ns)*cpk(ns)
           cvmix = cvmix + Yk(ns)*cvk(ns)
          enddo
           
          gam = cpmix/cvmix


          if(p.le.0.0d0.or.rho.le.0.0d0) then
            print*, 'sound speed not physical ', p, rho 
            stop
          endif
 
          c = sqrt(gam*p/rho)


!          if(c.lt.50.0d0.or.c.gt.1000.0d0) then
!           print*, 'something awry with c ', c, p, rho
!           stop
!          endif 

         end subroutine

!--------------------------------------------------------------

        subroutine compute_pres(rho, yk, ie, pres) 

        USE THERMO 
 
        implicit none

        real(kind=8) :: rho, yk(1:nspeci), ie, pres
        real(kind=8) :: temp, molwt, Rgas  
        integer :: ns

               if(ie.le.0.0d0) then
                print*, 'ie < 0 ', ie, rho, yk
                stop
               endif

        call FROM_IE_TO_T(ie, YK, TEMP)
       
            if(temp.le.0.0d0) then
             print*, 'temp < 0 ', temp, yk, ie
             stop 
            endif
 
        molwt = 0.0d0
        do ns = 1, nspeci
          molwt = molwt + yk(ns)/MW(ns)
        enddo 
        molwt = 1.0d0/molwt
        Rgas = runiv/molwt 

        pres = rho*Rgas*temp


        end subroutine



!--------------------------------------------------------------

        subroutine compute_temperature(rho, pres, yk, T)

        USE THERMO

        implicit none

        real(kind=8) :: rho, pres, T, yk(1:nspeci)
        real(kind=8) :: molwt, Rgas
        integer :: ns

        molwt = 0.0d0
        do ns = 1, nspeci
          molwt = molwt + yk(ns)/MW(ns)
        enddo
        molwt = 1.0d0/molwt
        Rgas = runiv/molwt

        T= pres/rho/Rgas

           T = max(T,200.0d0)


        end subroutine
  


!--------------------------------------------------------------

        SUBROUTINE HTOTAL(TEMP,HT)

      IMPLICIT NONE
      REAL(KIND=8) :: TEMP, TT, A1, A2, A3, A4, A5, A6, A7, HT, RMW
      REAL(KIND=8) :: HTOINT 

       RMW = 26.98D0
       TT = MIN(TEMP,3000.0D0)

         IF(TT.GT.600.D0) THEN
            A1 =  0.02559589D+02
            A2 = -0.10632239D-03
            A3 =  0.07202828D-06
            A4 = -0.02121105D-09
            A5 =  0.02289429D-13
            A6 =  0.03890214D+06
            A7 =  0.05234522D+02
         ELSE
            A1 =  0.02736825D+02
            A2 = -0.05912374D-02
            A3 = -0.04033937D-05
            A4 =  0.02322343D-07
            A5 = -0.01705599D-10
            A6 =  0.03886794D+06
            A7 =  0.04363879D+02
         END IF

       HTOINT =  A5 * TT**5 / 5.D0 +   &
                 A4 * TT**4 / 4.D0 +   &
                 A3 * TT**3 / 3.D0 +   &
                 A2 * TT**2 / 2.D0 +   &
                 A1 * TT           +   &
                 A6

       HT = HTOINT * 8314.15D0 / RMW

      RETURN
      END SUBROUTINE

!-------------------------------------------------------

