
         subroutine hllc_2d(gaml, gamr, rhol, rhor, pl, pr, ul, ur, vl, vr, iel,   &
                   ier, ysl, ysr, sweep, flux)

              use thermo

        real(kind=8) :: gaml, gamr, rhol, rhor, pl, pr, ul, ur, vl, vr, iel,  &
                   ier, ysl(1:2), ysr(1:2)
        real(kind=8) :: spl, spr, term1, term2, sm, pstar
        real(kind=8) :: rhostarl, rhostarr, ustarl(1:6), ustarr(1:6)
        real(kind=8) :: flux(1:6), energyl, energyr
        integer :: sweep   
        real(kind=8) :: ql, qr, qtl, qtr, omegal, omegar, cl, cr
        real(kind=8) :: rho_roe, ys_roe(1:2), u_roe, v_roe, h_roe, hl, hr
        real(kind=8) :: temp_roe, molwt, Rgas, p_roe, q_roe, c_roe  

! HLLC FLUX FOR THE GAS

               rho_roe = sqrt(rhol*rhor)
               ys_roe(1) = (sqrt(rhol)*ysl(1) + sqrt(rhor)*ysr(1))/(sqrt(rhol)+sqrt(rhor))
               ys_roe(2) = (sqrt(rhol)*ysl(2) + sqrt(rhor)*ysr(2))/(sqrt(rhol)+sqrt(rhor))
               u_roe = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol) + sqrt(rhor))   
               v_roe = (sqrt(rhol)*vl + sqrt(rhor)*vr)/(sqrt(rhol) + sqrt(rhor))   
               hl = iel + pl/rhol + 0.5d0*(ul*ul + vl*vl)
               hr = ier + pr/rhor + 0.5d0*(ur*ur + vr*vr)
               h_roe = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol) + sqrt(rhor))
               h_roe = h_roe - 0.5d0*(u_roe*u_roe + v_roe*v_roe)       
 
         call FROM_IH_TO_T(h_roe, ys_roe, temp_roe)
         
          molwt = ys_roe(1)/MW(1) + ys_roe(2)/MW(2)
          molwt = 1.0d0/molwt
          Rgas = runiv/molwt
          p_roe = rho_roe*Rgas*temp_roe
          call compute_sound(p_roe,rho_roe,ys_roe,c_roe) 
            

          call compute_sound(pl,rhol,ysl,cl) 
          call compute_sound(pr,rhor,ysr,cr) 



                   if(iel.le.0.0d0) then
                     print*, 'iel ', iel
                     stop
                   endif 
                   if(ier.le.0.0d0) then
                     print*, 'ier ', ier
                     stop
                   endif 

              if(abs(gaml-gamr).gt.1.0d-4) then
               print*, 'gamma: ', gaml, gamr
               stop
              endif


            if(sweep.eq.1) then
              ql = ul
              qr = ur
              qtl = vl
              qtr = vr
              q_roe = u_roe
            else if(sweep.eq.2) then
              ql = vl
              qr = vr
              qtl = ul
              qtr = ur
              q_roe = v_roe
            else
             print*, 'bug in sweep ', sweep
             stop
            endif

            SPL = min(ql - cl, q_roe - c_roe)
            SPR = min(qr + cr, q_roe + c_roe)
            TERM1 = RHOR*qr*(SPR-qr) - RHOL*ql*(SPL-ql)+pl-pr
            TERM2 = RHOR*(SPR-qr) - RHOL*(SPL-ql)
            SM = TERM1/TERM2

           PSTAR = RHOL*(ql-SPL)*(ql-SM) + pl

           RHOSTARL = RHOL*(SPL-ql)/(SPL-SM)
           RHOSTARR = RHOR*(SPR-qr)/(SPR-SM)
           IESTARL = PSTAR/(gaml-1.0D0)/RHOSTARL
           IESTARR = PSTAR/(gamr-1.0D0)/RHOSTARR

           USTARL(1) = RHOSTARL
           omegal = -(SM-ql)/(SPL-SM)*(ql-SPL)
           if(sweep.eq.1) then
            !USTARL(2) = RHOSTARL*SM
            !USTARL(3) = RHOSTARL*vl
            USTARL(2) = RHOSTARL*ul + RHOL*omegal
            USTARL(3) = RHOSTARL*vl  
           else if(sweep.eq.2) then
            !USTARL(2) = RHOSTARL*ul
            !USTARL(3) = RHOSTARL*SM
            USTARL(2) = RHOSTARL*ul
            USTARL(3) = RHOSTARL*vl + RHOL*omegal
           else
             print*, 'bug in sweep ', sweep
             stop
           endif 
           !USTARL(4) = RHOSTARL*(SM*SM/2.0D0 + IESTARL)
           ENERGYL = rhol*(IEL + 0.5D0*(ul*ul + vl*vl)) 
           USTARL(4) = ((SPL-ql)*ENERGYL + PSTAR*SM - pl*ql)/(SPL - SM)
           USTARL(5) = RHOSTARL*ysl(1)
           USTARL(6) = RHOSTARL*ysl(2)

           USTARR(1) = RHOSTARR
           omegar = -(SM-qr)/(SPR-SM)*(qr-SPR)
           if(sweep.eq.1) then
            !USTARR(2) = RHOSTARR*SM
            !USTARR(3) = RHOSTARR*vr
            USTARR(2) = RHOSTARR*ur + RHOR*omegar
            USTARR(3) = RHOSTARR*vr  
           else if(sweep.eq.2) then
            !USTARR(2) = RHOSTARR*ur
            !USTARR(3) = RHOSTARR*SM
            USTARR(2) = RHOSTARR*ur
            USTARR(3) = RHOSTARR*vr + RHOR*omegar
           else
             print*, 'bug in sweep ', sweep
             stop
           endif
           !USTARR(4) = RHOSTARR*(SM*SM/2.0D0 + IESTARR)
           ENERGYR = rhor*(IER + 0.5D0*(ur*ur + vr*vr))
           USTARR(4) = ((SPR-qr)*ENERGYR + PSTAR*SM - pr*qr)/(SPR - SM)
           USTARR(5) = RHOSTARR*ysr(1)
           USTARR(6) = RHOSTARR*ysr(2)

           IF(SPL.GT.0.0D0) THEN
            FLUX(1) = RHOL*ql
            if(sweep.eq.1) then
             FLUX(2) = RHOL*ql*ql + pl
             FLUX(3) = RHOL*ql*qtl
            else if(sweep.eq.2) then
             FLUX(2) = RHOL*ql*qtl 
             FLUX(3) = RHOL*ql*ql + pl
            else
             print*, 'bug in sweep ', sweep
             stop
            endif
            ENERGYL = (0.5D0*(ql*ql + qtl*qtl) + IEL)*RHOL
            FLUX(4) = (ENERGYL+pl)*ql
            FLUX(5) = RHOL*ql*ysl(1)
            FLUX(6) = RHOL*ql*ysl(2)
           ENDIF
           IF(SPL.LE.0.0D0.AND.SM.GT.0.0D0) THEN
            FLUX(1) = RHOL*ql + SPL*(USTARL(1)-RHOL)
            if(sweep.eq.1) then
             FLUX(2) = RHOL*ql*ql + pl + SPL*(USTARL(2)-RHOL*ql)
             FLUX(3) = RHOL*ql*qtl + SPL*(USTARL(3)-RHOL*qtl)
            else if(sweep.eq.2) then
             FLUX(2) = RHOL*ql*qtl + SPL*(USTARL(2)-RHOL*qtl)
             FLUX(3) = RHOL*ql*ql + pl + SPL*(USTARL(3)-RHOL*ql)
            else
             print*, 'bug in sweep ', sweep
             stop
            endif 
            ENERGYL = (0.5D0*(ql*ql + qtl*qtl) + IEL)*RHOL
            FLUX(4) = (ENERGYL+pl)*ql + SPL*(USTARL(4)-ENERGYL)
            FLUX(5) = RHOL*ql*ysl(1) + SPL*(USTARL(5)-RHOL*ysl(1))
            FLUX(6) = RHOL*ql*ysl(2) + SPL*(USTARL(6)-RHOL*ysl(2))
           ENDIF
           IF(SM.LE.0.0D0.AND.SPR.GE.0.0D0) THEN
            FLUX(1) = RHOR*qr + SPR*(USTARR(1)-RHOR)
            if(sweep.eq.1) then
             FLUX(2) = RHOR*qr*qr + pr + SPR*(USTARR(2)-RHOR*qr)
             FLUX(3) = RHOR*qr*qtr + SPR*(USTARR(3)-RHOR*qtr)
            else if(sweep.eq.2) then
             FLUX(2) = RHOR*qr*qtr + SPR*(USTARR(2)-RHOR*qtr)
             FLUX(3) = RHOR*qr*qr + pr + SPR*(USTARR(3)-RHOR*qr)
            else
             print*, 'bug in sweep ', sweep
             stop
            endif  
            ENERGYR = (0.5D0*(qr*qr + qtr*qtr) + IER)*RHOR
            FLUX(4) = (ENERGYR+pr)*qr + SPR*(USTARR(4)-ENERGYR)
            FLUX(5) = RHOR*qr*ysr(1) + SPR*(USTARR(5)-RHOR*ysr(1))
            FLUX(6) = RHOR*qr*ysr(2) + SPR*(USTARR(6)-RHOR*ysr(2))
           ENDIF
           IF(SPR.LT.0.0D0) THEN
            FLUX(1) = RHOR*qr
            if(sweep.eq.1) then
             FLUX(2) = RHOR*qr*qr + pr
             FLUX(3) = RHOR*qr*qtr 
            else if(sweep.eq.2) then
             FLUX(2) = RHOR*qr*qtr  
             FLUX(3) = RHOR*qr*qr + pr 
            else
             print*, 'bug in sweep ', sweep
             stop
            endif
            ENERGYR = (0.5D0*(qr*qr + qtr*qtr) + IER)*RHOR
            FLUX(4) = (ENERGYR+pr)*qr
            FLUX(5) = RHOR*qr*ysr(1)
            FLUX(6) = RHOR*qr*ysr(2)
           ENDIF


               return
               end subroutine
!----------------------------------------------------------------------------------


 
