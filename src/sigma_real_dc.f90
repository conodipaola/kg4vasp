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

module L11_L12_exct_dc
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
    real(q) :: domega,sigma,T
    real(q) :: tmp1,tmp2,tmp3,tmp4,tmp5

    real(q) :: calc1, denerg, denerg_1
    real(q) :: psi_grad_psi_aux(3), t11tmp(3,3), t12tmp(3,4), t22tmp(3,3)

   integer  :: ios , ierr ! Input and memory allocation error status variables
! global variables
   real(q) :: L11_dc_tensor(3,3),L12_dc_tensor(3,3),L22_dc_tensor(3,3)
   real(q) :: L11_dc_trace, L12_dc_trace, L22_dc_trace,S,KC
!
contains
!!!!!!!!!!!!!!!!!
!!!
subroutine L11_L12_exct_dc_trace
! Gaussian representation of Dirac's delta: exp(-e*e)
!    sigma=dispers/AUTOEV
    sigma=dispers
    calc1=dsqrt(pi)/sigma

    L11_dc_tensor=0.0_q
    L12_dc_tensor=0.0_q
    L22_dc_tensor=0.0_q
!
    t11tmp=0.0_q
    t12tmp=0.0_q
    t22tmp=0.0_q

    do icart=1,3   ! starting cycle over cartesian coordinates x,y,z
        do jcart=1,3  ! starting second cycle over cartesian coordinates x,y,z
!                sigma1=0.0_q
!
            do ikpoints=1,nkpoints  ! starting cycle over kpoints
!
                do ibands=1,nbands  ! starting cycle over bands
                    do jbands=1,nbands   ! starting second cycle over bands

                        if ( (.not.intrab) .and. (ibands==jbands) )  cycle
!
                        if (.not.FROM_VASP) then
                            denerg=(en(ibands,ikpoints) - en(jbands,ikpoints)) * 0.5d0
                        else if (FROM_VASP) then
                            denerg=(en(ibands,ikpoints) - en(jbands,ikpoints))
                        endif
!
!
                        if (dabs(denerg) <= de) then
                            if (.not.intrab) cycle
                            pref = 1.0_q / (dexp((en(jbands,ikpoints)-e_f)/dsigma)+1.d0)
                            dfermi = -wk(ikpoints)*pref*(pref-1.d0)/dsigma*2.d0  ! 2.d0 facto is for ISPIN=1 (non spin-polarised)
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
                        endif

!                        if ( ibands /=  jbands ) dfermi = 2.0_q * dfermi
!
                        tmp1=real(conjg(psi_nabla_psi(icart,jbands,ibands,ikpoints)) &
                        *  psi_nabla_psi(jcart,jbands,ibands,ikpoints))
                        tmp2=dexp((-denerg*denerg)/(sigma*sigma))
                        tmp3=en(ibands,ikpoints)-e_f
                        tmp4=en(jbands,ikpoints)-e_f
                        tmp5=dfermi &
                        * tmp1 &
                        * tmp2
                        t11tmp(jcart,icart) = t11tmp(jcart,icart)&
                            + tmp5
                        t12tmp(jcart,icart) = t12tmp(jcart,icart)&
                            + tmp5 &
                            * tmp3
                        t22tmp(jcart,icart) = t22tmp(jcart,icart)&
                            + tmp5 &
                            * tmp3 &
                            * tmp4
                    enddo ! end jbands
                enddo ! end ibands
!
            enddo ! end ikpoints
!
        enddo ! end jcart
    enddo ! end icart
!
    L11_dc_tensor(1:3,1:3)=(((HBAR**3)*(EVTOJ**3))/((ELEMASS**2)*(ANGTOMT**5)))*(t11tmp(1:3,1:3))&
        * calc1 / volume
    L12_dc_tensor(1:3,1:3)=(((HBAR**3)*(EVTOJ**3))/((ELEMASS**2)*(ANGTOMT**5)))*(t12tmp(1:3,1:3))&
        * calc1 / volume
    L22_dc_tensor(1:3,1:3)=(((HBAR**3)*(EVTOJ**3))/((ELEMASS**2)*(ANGTOMT**5)))*(t22tmp(1:3,1:3))&
        * calc1 / volume
    L11_dc_trace=(L11_dc_tensor(1,1) &
                    + L11_dc_tensor(2,2) &
                    + L11_dc_tensor(3,3)) / 3.0_q
    L12_dc_trace=(L12_dc_tensor(1,1) &
                    + L12_dc_tensor(2,2) &
                    + L12_dc_tensor(3,3)) / 3.0_q
    L22_dc_trace=(L22_dc_tensor(1,1) &
                    + L22_dc_tensor(2,2) &
                    + L22_dc_tensor(3,3)) / 3.0_q
!
    OPEN(unit=11,file='data_transport.dat',STATUS='NEW')
    write(11,*) '===== sigma real part DC exact tensor ===='
    do icart=1,3
        write(11,'(2x,3es20.7e3)')  ( L11_dc_tensor(icart,jcart), jcart=1,3 )
    enddo
    write(11,*) 'Trace average of elec conductivity DC: ', L11_dc_trace
    write(11,*) '=================================='
    write(11,*) '===== L12 DC exact tensor ===='
    do icart=1,3
        write(11,'(2x,3es20.7e3)')  ( L12_dc_tensor(icart,jcart), jcart=1,3 )
    enddo
    write(11,*) 'Trace average of L12 DC exact: ', L12_dc_trace
    write(11,*) '=================================='
    write(11,*) '===== L22 DC exact tensor ===='
    do icart=1,3
        write(11,'(2x,3es20.7e3)')  ( L22_dc_tensor(icart,jcart), jcart=1,3 )
    enddo
    write(11,*) 'Trace average of L22 DC exact: ', L22_dc_trace
    write(11,*) '=================================='
!
    S=0.0_q
    T=dsigma/KB
    S=T*L11_dc_trace
    S=-L12_dc_trace/S
    write(11,*) 'Temp (K)=',T,'Elec. Conductivity=',L11_dc_trace,' (1/(Ohm*m))'
    write(11,*) 'Temp (K)=',T,'Seebeck=',S,' (V/K)'
!
    KC=0.0_q
    KC=L12_dc_trace**2
    KC=KC/L11_dc_trace
    KC=L22_dc_trace-KC
    KC=KC/T
    write(11,*) 'Temp (K)=',T,'Thermal Conductivity=',KC,'[W/(m*K)]'
!
    CLOSE(11)

return
end subroutine L11_L12_exct_dc_trace
!
end module L11_L12_exct_dc
