                                     program post

                implicit none

                 integer, parameter :: imax = 501
                 integer, parameter :: jmax = 501
                 character ( len = 100 ) filename
                 integer :: n, i, j, l, ii
                 real(kind=8) :: prim(1:imax,1:jmax,1:13), Yk1(1:imax,1:jmax)
                 real(kind=8) :: dx, dy
                 logical :: file_exists
                 real(kind=8), dimension(1:imax,1:jmax) :: x, y        
                 real(kind=8) :: rc, xc, yc
                 integer, parameter :: nbins = 1000
                 real(kind=8), dimension(1:nbins) :: r, yk, entries, sumyk
                 real(kind=8) :: low, high, r0, time, delta
                 real(kind=8) :: DM, term1, term2, term3


                 dx = 1.5d0/(imax-1)
                 dy = 1.5d0/(jmax-1)

                 r0 = 0.1d0 



                 do i = 1, imax
                  do j = 1, jmax
                   x(i,j) = dble(i-1)*dx
                   y(i,j) = dble(j-1)*dy 
                  enddo
                 enddo


                 do ii = 1, nbins
                  r(ii) = sqrt(dx**2.0d0+dy**2.0d0)*dble(ii-1)
                 enddo


                 open(17,file='ml.txt',form='formatted')

                 do n = 0, 50000, 25
                  write(filename,'("output/output_",I5.5,".dat")'),n
                  INQUIRE(FILE=filename, EXIST=file_exists)

                       if(file_exists) then



                  open(34,file=filename,form='formatted')
                  read(34,*) time
               do i = 1, imax-1
                do j = 1, jmax-1
                 read(34,*) (prim(i,j,l),l=1,13)
                enddo
               enddo
                  close(34)

                  Yk1(:,:) = prim(:,:,6)


                  yk = 0.0d0
                  sumyk = 0.0d0
                  entries = 0.0d0
                  ! cylindrical average
                  do i = 1, imax
                   do j = 1, jmax
                     xc = 0.5d0*(x(i,j)+x(i+1,j)) 
                     yc = 0.5d0*(y(i,j)+y(i,j+1)) 
                     rc = sqrt(xc**2.0d0 + yc**2.0d0)
                     do ii = 1, nbins
                      if(rc.ge.r(ii).and.rc.lt.r(ii+1)) then
                        sumyk(ii) = sumyk(ii) + yk1(i,j)
                        entries(ii) = entries(ii) + 1.0d0
                      endif  
                     enddo 
                   enddo 
                  enddo 

                 
                 do ii = 1, nbins
                  if(entries(ii).gt.0.0d0) then
                   yk(ii) = sumyk(ii)/entries(ii)
                  else
                   yk(ii) = 0.0d0 
                  endif 
                 enddo   


                  low = 0.0d0
                  high = 0.0d0
                  do ii = 1, nbins-1
                   if(yk(ii).ge.0.95d0.and.yk(ii+1).le.0.95d0) then
                   !if(yk(ii).ge.0.925d0.and.yk(ii+1).le.0.925d0) then
                     low = r(ii)
                     goto 123
                   endif
                  enddo  

123               continue

                  do ii = nbins-1, 1, -1
                   if(yk(ii).ge.0.05d0.and.yk(ii+1).le.0.05d0) then
                   !if(yk(ii).ge.0.075d0.and.yk(ii+1).le.0.075d0) then
                     high = r(ii)
                     goto 124
                   endif 
                  enddo 

124               continue


                  high = high/r0
                  low = low/r0  
                  delta = high - low



                  term1 = 0.0d0
                  term2 = 0.0d0
                  term3 = 0.0d0
                  DM = 0.0d0
               do i = 1, imax-1
                do j = 1, jmax-1
                 term1 = term1 + prim(i,j,6)
                 term2 = term2 + prim(i,j,7)
                 term3 = term3 + prim(i,j,6)*prim(i,j,7)
                enddo
               enddo
               term1 = term1/dble(imax-1)/dble(jmax-1)
               term2 = term2/dble(imax-1)/dble(jmax-1)
               term3 = term3/dble(imax-1)/dble(jmax-1)

               DM = term3/(term1*term2)


                  write(17,*) time*(1.0d3), low, high, delta, DM
                  print*, time*(1.0d3), low, high, delta, DM
                  call flush()



                   endif

                 enddo 

                 close(17)
                  

                               end program
