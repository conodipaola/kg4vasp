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
! 
program kubo_greenwood
    use prec
    use params
    use input
    use nabla
    use sigma_part_real
    use sum_rule
    use L11_L12_exct_dc
!    use read_nabla
!    implicit none

!    real(q) ::
!    real(q) ::

!    integer ::

!    character (LEN=) ::
!
    OPEN(UNIT=12,FILE='input_parameters')
!
    write(*,*) '.....READING PARAMETERS FILE.....'
    call read_input()
    write(*,*) '.....ECHO PARAMETRS IN FILE.....'
    call write_param(12)
    write(*,*) 'JOB DONE!'
    write(*,*) '.....READING NABLA MATRIX ELEMENTS.....'
    call read_nabla()
    write(*,*) 'JOB DONE!'
!    open(unit=10, file='')
!    OPEN(UNIT=13,FILE='input_parameters_check')
!    write(13,*) '===== INPUT PARAMETERS ======'
!    write(13,*) '-----------------------------'
!    write(13,*) 'AUTOEV = ',AUTOEV
!    write(13,*) 'FILE_NABLA = ',FNAME
!    write(13,*) 'NBANDS = ',nbands
!    write(13,*) 'NKPOINTS = ',nkpoints
!    write(13,*) 'NP_FREQ = ', nfreq
!    write(13,*) 'DELTAE = ', de
!    write(13,*) 'MAX_FREQ = ', omegamx
!    write(13,*) 'MIN_FREQ = ', omegamn
!    write(13,*) 'GAUS_SIGMA = ', dispers,' (in eV)'
!    write(13,*) 'VOLUME = ',volume
!    write(13,*) 'DSIGMA = ',dsigma
!    write(13,*) 'INTRA_BAND = ', intrab
!    write(13,*) 'FERMI_ENERGY = ',e_f
!    write(13,*) 'OCC_FROM_CODE = ', FROM_CODE
!    write(13,*) 'MATRIX ELEMENTS FROM VASP = ', FROM_VASP
!    write(13,*) IU1
!
    write(*,*) '.....ELEC. CONDUCTIVITY VS. FREQUENCY: STARTING CALCULATION.....'
    call freq_binning()
!
    call real_exact_conductivity()
    write(*,*) 'JOB DONE!!'
    call sumrule_1()
    call sumrule_2()
!    write(13,*) 'sum rule (without frequency dependency) = ', sumrule1, ' over ',nbands,'bands'
!    write(13,*) 'sum rule (integral with frequency dependency) = ', sumrule2, ' over ',nbands,'bands'
!    CLOSE(13)
!
    write(12,*) 'sum rule (without frequency dependency) = ', sumrule1, ' over ',nbands,'bands'
    write(12,*) 'sum rule (integral with frequency dependency) = ', sumrule2, ' over ',nbands,'bands'
    CLOSE(12)
    call L11_L12_exct_dc_trace ()
!    CLOSE(13)
!    
stop 
end program kubo_greenwood 
