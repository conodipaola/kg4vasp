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

module sum_rule
    use prec
    use params
    use sigma_part_real, ONLY: domega,sigma_real,sigma_real1
    use nabla, ONLY : fw,en,psi_nabla_psi,psi_nabla_psi_i,wk
!
    implicit none
    real(q) :: rule1,sumrule1,sumrule2
    real(q) :: temp,denerg,dfermi
    integer :: ibands,jbands,ikpoints,icart
contains
!==================================
subroutine sumrule_1

    rule1=0.d0

    do ikpoints=1,nkpoints
!
        do ibands = 1,nbands
            do jbands = 1,nbands
!
                if (ibands == jbands) cycle ! no cycle between same bands
!
                denerg=(en(ibands,ikpoints) - en(jbands,ikpoints)) * 0.5d0
                if (dabs(denerg) <= de) cycle ! no degeneracies
!
                if (.not.FROM_VASP) then
                    dfermi = (fw(jbands,ikpoints)-fw(ibands,ikpoints))*nkpoints / denerg
                else if (FROM_VASP) then
                    dfermi = (fw(jbands,ikpoints)-fw(ibands,ikpoints))*nkpoints/2.0_q / denerg
                endif
                !
                temp = 0.d0
                do icart=1,3
                    temp =  temp &
                            + conjg(psi_nabla_psi(icart,jbands,ibands,ikpoints)) &
                            * psi_nabla_psi(icart,jbands,ibands,ikpoints)
                enddo ! end cartesian 
!
                rule1 = rule1 + dfermi * temp
!
            enddo ! end bands internal cycle
        enddo ! end bands external cycle
    enddo ! end kpoints cycle
!
    sumrule1 = rule1 / (3.d0 * nelec)
!    write(35,*) nelec,sumrule1
return
end subroutine sumrule_1
!==================================================
subroutine sumrule_2
    real(q) :: sum1, sum2
    integer :: ifreq
    
    ! 
    ! Sum rule via Simpson's integration (with delta function)
    ! 
    
    sum1 = 0.0_q
    sum2 = 0.0_q
!    write(91,*) sigma_real(:,:,0)
!    write(91,*) sigma_real(:,:,nkpoints)
!    write(91,*) nkpoints

!    do ifreq=0, nfreq
!        write(301,*) sigma_real1(1:3,1:3,ifreq)/((0.0367493**2)*(1.889725989**5))
!    enddo
 
    do ifreq=2, nfreq-2, 2
        sum1 = sum1 + (sigma_real(1,1,ifreq) + sigma_real(2,2,ifreq) + sigma_real(3,3,ifreq))
!        write(91,*) 'even: ',ikpoints
    enddo
 
    do ifreq=1, nfreq-1, 2
        sum2 = sum2 + (sigma_real(1,1,ifreq) + sigma_real(2,2,ifreq) + sigma_real(3,3,ifreq))
!        write(91,*) 'odd; ',ikpoints
    enddo

!    write(91,*) sum1,sum2,domega,volume,nelec,pi,nfreq
 
    sum1 = 2.0_q * sum1 + 4.0_q * sum2
    sum1 = sum1 + (sigma_real(1,1,0) + sigma_real(2,2,0) + sigma_real(3,3,0))
    sum1 = (sum1 +(sigma_real(1,1,nfreq)+sigma_real(2,2,nfreq)+sigma_real(3,3,nfreq)))*domega/3.0_q
 
    sumrule2 = 2.0_q * volume * sum1 / ( 3.0_q * pi * real(nelec) )
!
    sumrule2 = (sumrule2 * ELEMASS *  (ANGTOMT**3))/(ELECT*ELECT*HBAR)

!    write(90,*) volume,pi,nelec,ELEMASS,ANGTOMT,ELECT,HBAR,domega,nkpoints
!    
end subroutine sumrule_2
!====================================
end module sum_rule
