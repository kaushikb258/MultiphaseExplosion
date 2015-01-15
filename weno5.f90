                 subroutine weno5(u,uiphm,uimhp)

           implicit none

           real(kind=8) :: u(1:5), uiphm, uimhp
           real(kind=8) :: d0, d1, d2, dt2, dt1, dt0
           real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, & 
                               omega0, omega1, omega2
           real(kind=8) :: alphat0, alphat1, alphat2
           real(kind=8) :: sum1, p0, p1, p2, eps
           integer :: i


                  i = 3


                  d0 = 3.0d0/10.0d0
                  d1 = 6.0d0/10.0d0
                  d2 = 1.0d0/10.0d0
                  dt2 = 3.0d0/10.0d0
                  dt1 = 6.0d0/10.0d0
                  dt0 = 1.0d0/10.0d0
                  eps = 1.0d-6


                  uiphm = 0.0d0
                  uimhp = 0.0d0


              ! compute u(i+1/2)-

           beta0 = 13.0d0/12.0d0*(u(i) - 2.0d0*u(i+1) + &
            u(i+2))**2.0d0 + 1.0d0/4.0d0*(3.0d0*u(i) - 4.0d0*u(i+1) + u(i+2))**2.0d0

           beta1 = 13.0d0/12.0d0*(u(i-1) - 2.0d0*u(i) + &
                u(i+1))**2.0d0 + 1.0d0/4.0d0*(u(i-1) - u(i+1))**2.0d0

           beta2 = 13.0d0/12.0d0*(u(i-2) - 2.0d0*u(i-1) + u(i))**2.0d0 + &
                 1.0d0/4.0d0*(u(i-2) - 4.0d0*u(i-1) + 3.0d0*u(i))**2.0d0

           alpha0 = d0/(eps + beta0)**2.0d0
           alpha1 = d1/(eps + beta1)**2.0d0
           alpha2 = d2/(eps + beta2)**2.0d0

           sum1 = alpha0 + alpha1 + alpha2

           omega0 = alpha0/sum1
           omega1 = alpha1/sum1
           omega2 = alpha2/sum1


           p0 = 1.0d0/3.0d0*u(i) + 5.0d0/6.0d0*u(i+1) - 1.0d0/6.0d0*u(i+2)
           p1 = -1.0d0/6.0d0*u(i-1) + 5.0d0/6.0d0*u(i) + 1.0d0/3.0d0*u(i+1)
           p2 = 1.0d0/3.0d0*u(i-2) - 7.0d0/6.0d0*u(i-1) + 11.0d0/6.0d0*u(i)


           uiphm = omega0*p0 + omega1*p1 + omega2*p2

                  ! compute u(i-1/2)+

           alphat0 = dt0/(eps + beta0)**2.0d0
           alphat1 = dt1/(eps + beta1)**2.0d0
           alphat2 = dt2/(eps + beta2)**2.0d0

           sum1 = alphat0 + alphat1 + alphat2

           omega0 = alphat0/sum1
           omega1 = alphat1/sum1
           omega2 = alphat2/sum1

           p0 = 11.0d0/6.0d0*u(i) - 7.0d0/6.0d0*u(i+1) + 1.0d0/3.0d0*u(i+2)
           p1 = 1.0d0/3.0d0*u(i-1) + 5.0d0/6.0d0*u(i) - 1.0d0/6.0d0*u(i+1)
           p2 = -1.0d0/6.0d0*u(i-2) + 5.0d0/6.0d0*u(i-1) + 1.0d0/3.0d0*u(i)


           uimhp = omega0*p0 + omega1*p1 + omega2*p2



           return
           end subroutine

!---------------------------------------------------------------------------


