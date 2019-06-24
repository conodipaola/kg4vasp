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

module sigma_part_real
    use prec
    use params
    use nabla, ONLY : fw,en,psi_nabla_psi,psi_nabla_psi_i,wk
    IMPLICIT NONE
! private variables
    integer :: ifreq                             ! counter for Frequency
    integer :: ikpoints, ibands, jbands          ! counters for k-points and internal and external bands
    integer :: ig                                ! counter on plane waves
    integer :: icart,jcart                       ! counters for x, y and z component

    real(q) :: pref, dfermi,pref1,pref2
    real(q) :: domega,sigma

    real(q) :: calc1, denerg, denerg_1
    real(q) :: psi_grad_psi_aux(3), SIGMA1(3,3)

   integer  :: ios , ierr ! Input and memory allocation error status variables
! global variables
   real(q), allocatable :: sigma_real(:,:,:),sigma_real1(:,:,:)
!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine freq_binning()
    !
    ! Frequency grid
    ! Input units: eV    Output Units: Ha
    
        if ( nfreq == 0 ) then
            domega = 0.0_q
        else
            domega   = ( omegamx - omegamn ) / dble(nfreq)  ! Frecuency step
        endif
    
        allocate( freq_g(0:nfreq), STAT=ierr )
    !    if (ierr/=0) call errore('kgsigma1','allocating wgrid', abs(ierr))
          
    ! Frequency grid in eV
        do ifreq = 0, nfreq
            freq_g(ifreq) = omegamn + ifreq * domega
        enddo
    
    
        write(12,*) '    Frequency grid:'
        write(12,*) "       Frequency step (eV)       :", domega
        write(12,*) "       Frequency steps           :", nfreq
        write(12,'(a, f10.4,a,f10.4 )') &
            "        Frequency range (eV)      :", freq_g(0),  ' to ',freq_g(nfreq)
    
    ! Change units from eV to Ha
    !     freq_g = freq_g / AUTOEV
    !     domega    = domega    / AUTOEV
    !
    end subroutine freq_binning
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
subroutine real_exact_conductivity
! Gaussian representation of Dirac's delta: exp(-e*e)
!    sigma=dispers/AUTOEV
    sigma=dispers
    calc1=dsqrt(pi)/sigma

    allocate(sigma_real(1:3,1:3,0:nfreq),STAT=ierr)
    allocate(sigma_real1(1:3,1:3,0:nfreq),STAT=ierr)
    sigma1=0.0_q
!
    OPEN(unit=10,FILE='elec_conductivity_vs_freq.dat',STATUS='NEW')
    write(10,'(2x,11(A20,2x))')'## freq (eV)', '11', '12', '13', '21', '22', '23',&
         '31', '32', '33','TRACE/3'
!    if (ierr/=0) then error()
!! debug
    do ifreq = 0,nfreq ! starting cycle over frequency range
        sigma1=0.0_q
!
        do icart=1,3   ! starting cycle over cartesian coordinates x,y,z
            do jcart=1,3  ! starting second cycle over cartesian coordinates x,y,z
!                sigma1=0.0_q
!
                do ikpoints=1,nkpoints  ! starting cycle over kpoints
!
                    do ibands=1,nbands  ! starting cycle over bands
                        do jbands=1,nbands   ! starting second cycle over bands
!
                            if (.not.FROM_VASP) then
                                denerg=(en(ibands,ikpoints) - en(jbands,ikpoints)) * 0.5d0
                            else if (FROM_VASP) then
                                denerg=(en(ibands,ikpoints) - en(jbands,ikpoints))
                            endif
!
                            denerg_1=denerg-freq_g(ifreq)
!                            write(100,*)denerg,dispers,denerg-freq_g(ifreq)
!
                            if (dabs(denerg) <= de) then
                                if (.not.intrab) cycle
                                pref = 1.0_q / (dexp((en(jbands,ikpoints)-e_f)/dsigma)+1.d0)
                                dfermi = -wk(ikpoints)*pref*(pref-1.d0)/dsigma*2.d0  ! 2.d0 facto is for ISPIN=1 (non spin-polarised)
!                                write(90,*) 'intra',jbands,ibands,ikpoints,icart,jcart,dfermi
!                                write(91,*) jbands,ibands,ikpoints,en(jbands,ikpoints),e_f,dsigma,pref,wk(ikpoints)
                            endif
!
                            if (dabs(denerg) > de) then
                                if (.not.FROM_CODE) then  ! fermi-dirac formula to be used: e_f=mu, dsigma=Kb*T(K)
                                    pref1=1.0_q / (dexp((en(jbands,ikpoints)-e_f)/dsigma)+1.d0)
                                    pref2=1.0_q / (dexp((en(ibands,ikpoints)-e_f)/dsigma)+1.d0)
                                    dfermi = wk(ikpoints)*(pref1-pref2) / denerg * 2.d0  !  2.d0 facto is for ISPIN=1 (non spin-polarised)
                                else if (FROM_CODE) then  ! occupation output from the code: VASP or QE
                                    dfermi = (fw(jbands,ikpoints)-fw(ibands,ikpoints)) / denerg
                                endif

!                                write(90,*) 'inter',jbands,ibands,ikpoints,icart,jcart,dfermi
                            endif
                            sigma1(jcart,icart) = sigma1(jcart,icart)&
                                + dfermi &
                                * real(conjg(psi_nabla_psi(icart,jbands,ibands,ikpoints)) &
                                *            psi_nabla_psi(jcart,jbands,ibands,ikpoints)) &
                                * dexp((-denerg_1*denerg_1)/(sigma*sigma))
!                            write(101,*) jcart,icart,sigma1(jcart,icart)
!                            write(101,*) ifreq,jcart,icart,ibands,jbands,real(conjg(psi_nabla_psi_i(icart,jbands,ibands,ikpoints)&
!                                *psi_nabla_psi_i(jcart,jbands,ibands,ikpoints)))
!                            write(102,*) ifreq,dfermi
!                            write(103,*) jcart,icart,ibands,jbands,dfermi,sigma,denerg_1,dexp((-denerg_1*denerg_1)/(sigma*sigma))
!
                        enddo ! end jbands
                    enddo ! end ibands
!
                enddo ! end ikpoints
!
            enddo ! end jcart
        enddo ! end icart
!
!        sigma_real(1:3,1:3,ifreq) = ((2.0_q*(HBAR**3)*(EVTOJ**3))/((ELEMASS**2)*(ANGTOMT**5)))*(sigma1(1:3,1:3))&
        sigma_real1(1:3,1:3,ifreq)=(sigma1(1:3,1:3))* calc1 / volume
!        write(300,*) sigma_real1(1:3,1:3,ifreq)/((0.0367493**2)*(1.889725989**5))
        sigma_real(1:3,1:3,ifreq) = (((HBAR**3)*(EVTOJ**3))/((ELEMASS**2)*(ANGTOMT**5)))*(sigma1(1:3,1:3))&
            * calc1 / volume
!
!            write(100,*) calc1, volume, pi, sigma
!        write(9,'(2x,11(F20.14,2x))') freq_g(ifreq),sigma1(1,1),sigma1(1,2),sigma1(1,3)&
!            ,sigma1(2,1),sigma1(2,2),sigma1(2,3)&
!            ,sigma1(3,1),sigma1(3,2),sigma1(3,3)&
!            ,((sigma1(1,1)+sigma1(2,2)+sigma1(1,1))/3.0_q)
        write(10,'(2x,11(E20.14,2x))') freq_g(ifreq),(((sigma_real(icart,jcart,ifreq)), icart=1,3),jcart=1,3), &
        ( (sigma_real(1,1,ifreq)+sigma_real(2,2,ifreq)+sigma_real(3,3,ifreq))/3.0_q )
!
!        write(9,'(2x,11(F20.14,2x))') freq_g(ifreq)*AUTOEV &
!            ,(((sigma_real(icart,jcart,ifreq)), jcart=1,3),icart=1,3), &
!            ( (sigma_real(1,1,ifreq)+sigma_real(2,2,ifreq)+sigma_real(3,3,ifreq))/3.0_q )
    enddo ! end ifreq
!
    CLOSE(10)
!
!    write(92,*) sigma_real
return
end subroutine real_exact_conductivity
!
end module sigma_part_real
