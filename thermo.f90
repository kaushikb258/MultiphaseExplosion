        MODULE   THERMO

        IMPLICIT NONE

       real(kind=8), save :: runiv = 8.314d3
       integer, parameter :: nspeci = 2   
       REAL(KIND=8), SAVE :: MW(1:NSPECI) = (/28.96d0, 28.96d0/) 

       REAL(KIND=8) :: THCOEF(1:7,1:2,1:NSPECI)
       REAL(KIND=8) :: TLOW(1:NSPECI), TUPP(1:NSPECI), TMID(1:NSPECI)
       CHARACTER (LEN = 10) :: SPNAME(1:NSPECI) 
       REAL(KIND=8) :: RGK(1:NSPECI) 


CONTAINS
         subroutine thermo_data() 

       INTEGER :: J, IOERROR, NS, NUNIT
       LOGICAL, DIMENSION(:), ALLOCATABLE :: EXISTDATA
       LOGICAL :: FLAG_EXIST
       REAL (KIND = 8) :: TLOWTEMP, TUPPTEMP, TMIDTEMP, TLOWDEF, TMIDDEF
       REAL (KIND = 8) :: TUPPDEF
       CHARACTER :: SPECIESNAME*18, DUMMY1*27


          do ns = 1, nspeci
           RGK(ns) = runiv/MW(ns)
          enddo


        return
        end subroutine 


       END MODULE

!--------------------------------------------------------------------------------


