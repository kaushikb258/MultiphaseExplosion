
        subroutine digit_to_ch ( digit, ch )

        implicit none
        character ch
        integer ( kind = 4 ) digit

        if ( 0 <= digit .and. digit <= 9 ) then
          ch = achar ( digit + 48 )
        else
         ch = '*'
         end if

        return
        end

!*****************************************************************************
       subroutine i4_to_s_left ( i4, s )

        implicit none

        character c
        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4
        integer ( kind = 4 ) idig
        integer ( kind = 4 ) ihi
        integer ( kind = 4 ) ilo
        integer ( kind = 4 ) ipos
        integer ( kind = 4 ) ival
        character ( len = * ) s

         s = ' '
         ilo = 1
         ihi = len ( s )
         if ( ihi <= 0 ) then
          return
         end if

         ival = i4

          if ( ival < 0 ) then
           if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
           end if

         ival = -ival
          s(1:1) = '-'
          ilo = 2

         end if

         ipos = ihi

        do

         idig = mod ( ival, 10 )
         ival = ival / 10

           if ( ipos < ilo ) then
       do i = 1, ihi
        s(i:i) = '*'
       end do
       return
        end if

         call digit_to_ch ( idig, c )

         s(ipos:ipos) = c
         ipos = ipos - 1

        if ( ival == 0 ) then
       exit
        end if

         end do

         s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
         s(ilo+ihi-ipos:ihi) = ' '

         return
        end

!*****************************************************************************
        subroutine vtk_write(output_unit,filename,title,imax,jmax,xyz,nprim,prim)

        USE THERMO 

        implicit none

         integer :: imax, jmax

         character ( len = 20 ) cell_size_string
         integer :: node
         character ( len = 20 ) node_num_string
         integer :: output_unit
         character ( len = 20 ) s_node_num, s_imax, s_jmax, s_cells
         character ( len = 100 ) title
         character ( len = 100 ) filename
         real(kind=8) :: xyz(1:imax,1:jmax,1:3)
         integer :: nprim 
         real(kind=8) :: prim(1:imax,1:jmax,1:nprim)
         real(kind=8) :: uvw(1:imax,1:jmax,1:3), uvw2(1:imax,1:jmax,1:3)
         integer :: i, j, k
         real(kind=8) :: ie, yk(1:2), temp(1:imax,1:jmax) 


         call i4_to_s_left ( imax*jmax, s_node_num )
         call i4_to_s_left ( (imax-1)*(jmax-1), s_cells )  ! don't use kmax for 2D
         call i4_to_s_left ( imax, s_imax )
         call i4_to_s_left ( jmax, s_jmax )

         uvw(1:imax,1:jmax,1) = prim(1:imax,1:jmax,2)
         uvw(1:imax,1:jmax,2) = prim(1:imax,1:jmax,3)
         uvw(1:imax,1:jmax,3) = 0.0d0
 
         uvw2(1:imax,1:jmax,1) = prim(1:imax,1:jmax,12)
         uvw2(1:imax,1:jmax,2) = prim(1:imax,1:jmax,13)
         uvw2(1:imax,1:jmax,3) = 0.0d0

 

         open(output_unit,file=filename,form='formatted')

         write ( output_unit, '(a)' ) '# vtk DataFile Version 2.0'
         write ( output_unit, '(a)' ) title
         write ( output_unit, '(a)' ) 'ASCII'
         write ( output_unit, '(a)' ) 'DATASET STRUCTURED_GRID'
       write ( output_unit, '(a)' ) 'DIMENSIONS ' // (s_imax) //        &
                        (s_jmax) // '1'  ! kmax = 1 for now
       write (output_unit, '(a)' ) 'POINTS ' // (s_node_num) // 'double'


         do j = 1, jmax
          do i = 1, imax
         write ( output_unit, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) xyz(i,j,1:3)
          end do
         end do

        write ( output_unit, '(a)' ) 'CELL_DATA ' // (s_cells)
        write ( output_unit, '(a)' ) 'POINT_DATA ' // (s_node_num)

        write ( output_unit, '(a)' ) 'SCALARS pressure double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,5),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'VECTORS velocity double'
        write ( output_unit, *) ((uvw(i,j,1:3),i=1,imax),j=1,jmax)

        write ( output_unit, '(a)' ) 'SCALARS density double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,1),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS internal_energy double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,4),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y1 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,6),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS Y2 double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,7),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS solid_density double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,8),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'VECTORS solid_velocity double'
        write ( output_unit, *) ((uvw2(i,j,1:3),i=1,imax),j=1,jmax)

        write ( output_unit, '(a)' ) 'SCALARS solid_energy double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,11),i=1,imax),j=1,jmax)  

        write ( output_unit, '(a)' ) 'SCALARS fign double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((prim(i,j,12),i=1,imax),j=1,jmax)  

        
        ! COMPUTE TEMPERATURE   
        do i = 1, imax
         do j = 1, jmax
          ie = prim(i,j,4)
          Yk(1:2) = prim(i,j,6:7)   
        call FROM_IE_TO_T(ie, YK, temp(i,j))
         enddo 
        enddo    

        write ( output_unit, '(a)' ) 'SCALARS temperature double '
        write ( output_unit, '(a)' ) 'LOOKUP_TABLE default '
        write ( output_unit, *) ((temp(i,j),i=1,imax),j=1,jmax)  


        call flush()
        close(output_unit)

        return
        end subroutine

!---------------------------------------------------------------------------------
