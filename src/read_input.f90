!
! Copyright (C) 2019- 
!
! This program is free software; you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2 of the License, or 
! any later version. See the file LICENSE in the root directory of the 
! present distribution, or http://www.gnu.org/copyleft/gpl.txt ,
! or contact developers via e-mail: cono.dipaola@gmail.com 
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! Author: Cono Di Paola
!-----------------------------------------------------------------------
module input
    use prec
    use params

    implicit none
!
!   service variables for input reading purposes
    INTEGER :: IDUM
    REAL(q) :: RDUM
    COMPLEX(q)  :: CDUM   
    LOGICAL  :: LDUM
    LOGICAL :: LOPEN
    CHARACTER (LEN=5) :: FINPUT
    CHARACTER (LEN=40) :: SZNAM
    CHARACTER (LEN=1) ::  CHARAC
!
    INTEGER :: IERR,N
!
contains
!=================================
subroutine read_input
! local variables
    integer :: iov,ipos
    integer, parameter :: nlen = 100
    character (len=nlen) :: str,str1,str2
!

    FINPUT="INPUT"
    N=0
    IU1=5
    LOPEN=.TRUE.
    OPEN(UNIT=IU1,FILE=FINPUT,STATUS='OLD')
!   init calculations parameters
    SZNAM='unknown system'
    FNAME='nabla.dat'
    nbands=0.0_q
    nkpoints=0
    nfreq   = 1000
    de = 0.0_q
    omegamx = 5.0_q
    omegamn = 0.01_q
    dispers  = 0.01_q
    volume   = 0.0_q
    dsigma   = 0.0_q
    intrab   = .TRUE.
    e_f      = 0.0_q
    FROM_VASP  = .TRUE.
    nelec = 0
    FROM_CODE  = .TRUE.
!
    n_bin=.FALSE.
!
!
    write(str1,*) "(a", nlen, ")"
!
    DO
        read(IU1,str1,IOSTAT=iov) str
        if (iov > 0) THEN
            write(*,*) 'something went wrong exit > 0',iov
            exit
        else if (iov < 0) then
            write(*,*) 'end of file reached exit < 0',iov
            exit
        else    
            ipos = scan(str,"=",back=.false.)
            !write(*,*) trim(str)
            !write(*,*) str
            !write(*,*) ipos
            !write(*,*) ipos1
            !write(*,*) ipos2
            !write(*,*) ipos3
            !write(*,*) trim(str(:ipos-1))
            str2=str(:ipos-1)
            call noblanks(str2)
            IF (ANY([ 'NBANDS    ', 'NKPOINTS  ', 'NP_FREQ   ', 'N_ELECTRON'] == trim(str2))) THEN
                SELECT CASE(str2)
                    CASE ('NBANDS')
                        read (str(1+ipos:),*) nbands
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) nbands
                    CASE ('NKPOINTS')
                        read (str(1+ipos:),*) nkpoints
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) nkpoints
                    CASE ('NP_FREQ')
                        read (str(1+ipos:),*) nfreq
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) nfreq
                    CASE ('N_ELECTRON')
                        read (str(1+ipos:),*) nelec
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) nelec
                END SELECT
              !write(*,*)  'found integer',str2
            END IF
            IF (ANY([ 'DELTAE      ', 'MIN_FREQ    ', 'MAX_FREQ    ', 'GAUS_SIGMA  ', 'VOLUME      ', &
            'DEGAUSS     ','FERMI_ENERGY'] == trim(str2))) THEN
                SELECT CASE(str2)
                    CASE ('DELTAE')
                        read (str(1+ipos:),*) de
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) de
                    CASE ('MIN_FREQ')
                        read (str(1+ipos:),*) omegamn
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) omegamn
                    CASE ('MAX_FREQ')
                        read (str(1+ipos:),*) omegamx
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) omegamx
                    CASE ('GAUS_SIGMA')
                        read (str(1+ipos:),*) dispers
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) dispers
                    CASE ('VOLUME')
                        read (str(1+ipos:),*) volume
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) volume
                    CASE ('DEGAUSS')
                        read (str(1+ipos:),*) dsigma
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) dsigma
                    CASE ('FERMI_ENERGY')
                        read (str(1+ipos:),*) e_f
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) e_f
                END SELECT
            END IF
            IF (ANY([ 'SYSTEM    ', 'FILE_NABLA'] == trim(str2))) THEN
                SELECT CASE(str2)
                    CASE ('SYSTEM')
                        read (str(1+ipos:),*) SZNAM
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) SZNAM
                    CASE ('FILE_NABLA')
                        read (str(1+ipos:),*) FNAME
                        !write(*,*) 'case select',str2,'='
                        !write(*,*)  FNAME
                END SELECT   
            END IF
            IF (ANY([ 'FROM_VASP    ', 'OCC_FROM_CODE', 'INTRA_BAND   ','NABLA_BIN    '] == trim(str2))) THEN
                SELECT CASE(str2)
                    CASE ('FROM_VASP')
                        read (str(1+ipos:),*) FROM_VASP
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) FROM_VASP
                    CASE ('OCC_FROM_CODE')
                        read (str(1+ipos:),*) FROM_CODE
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) FROM_CODE
                    CASE ('INTRA_BAND')
                        read (str(1+ipos:),*) intrab
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) intrab
                    CASE ('NABLA_BIN')
                        read (str(1+ipos:),*) n_bin
                        !write(*,*) 'case select',str2,'='
                        !write(*,*) n_bin 
                END SELECT
            END IF
        ENDIF
    ENDDO
!
    CLOSE(IU1)
!
return
end subroutine read_input
!==============
subroutine write_param(IU3)
!
    integer :: ierr1,IU3
!
    ierr1=0
!    OPEN(UNIT=IU3,FILE='input_parameters')
    write(IU3,*) '===== INPUT PARAMETERS ======'
    write(IU3,*) '-----------------------------'
    write(IU3,*) 'SYSTEM = ',SZNAM
    write(IU3,*) 'FILE_NABLA = ',FNAME
    if (nbands == 0) then 
        write(IU3,*) 'NBANDS = ',nbands, '! ERROR: UNDEFINED CRUCIAL DATA'
        ierr1=ierr1+1
    else if (nbands /= 0) then 
        write(IU3,*) 'NBANDS = ',nbands
    endif
    if (nkpoints == 0) then
        write(IU3,*) 'NKPOINTS = ',nkpoints,'! ERROR: ENDEFINED CRUCIAL DATA'
        ierr1=ierr1+1
    else if (nkpoints /= 0) then
        write(IU3,*) 'NKPOINTS = ',nkpoints
    endif
    write(IU3,*) 'NP_FREQ = ', nfreq
    write(IU3,*) 'DELTAE = ', de
    write(IU3,*) 'MAX_FREQ = ', omegamx
    write(IU3,*) 'MIN_FREQ = ', omegamn
    write(IU3,*) 'GAUS_SIGMA = ', dispers
    write(IU3,*) 'VOLUME = ', volume
    write(IU3,*) 'DSIGMA = ',dsigma
    write(IU3,*) 'INTRA_BAND = ', intrab
    write(IU3,*) 'FERMI_ENERGY = ',e_f
    write(IU3,*) 'OCC_FROM_CODE = ', FROM_CODE
    write(IU3,*) 'MATRIX ELEMENTS FROM VASP = ', FROM_VASP
    write(IU3,*) 'NABLA FORMAT BINARY? = ', n_bin
!
   if (ierr1 /= 0) then
      write(*,*) 'CRUCIAL ',ierr1,' PARAMETERS MISSED, PLEASE CHECK INPUT.'
      write(*,*) 'I AM STOPPING RUN'
      stop
   endif
!
!    CLOSE(12)
!
return
end subroutine write_param
!
!================
subroutine noblanks(str1)
implicit none

character (len=100) :: str1,str2
integer :: i,ipos,ipos1

str2=''                          !in this variable will be save non-blank spaces  
!str1='   a   b   c  de   '            !Test string with blank spaces

!write(*,*)'sub1= ',len_trim(str1), str1

do i=1,len(str1)   
   if (str1(i:i).ne.' ')str2=trim(str2)//trim(str1(i:i))   
end do
ipos = scan(str2,"#",back=.false.)
ipos1= scan(str2,"!",back=.false.)

if ( ipos == 1 .or. ipos1 == 1 ) then
    write(*,*) 'COMMENT:   ',trim(str2)
    str2='NULL'
else if (ipos > 1 .or. ipos1 >1 ) then
    write(*,*) trim(str2),'  NOT A PARAMETER, CHECK INPUT'
    str2='NULL' 
endif

!write(*,*)'sub1= ',len_trim(str2), str2

str1=str2

return
end subroutine noblanks
!
end module
