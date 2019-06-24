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

module params
!
    use PREC
!
    IMPLICIT NONE
!
                                          ! Dirac delta function in eV
    integer  :: nfreq              ! Number of frequency points
    real(q) :: dsigma
    real(q) :: e_f   ! energy of the fermi level
!    real(q), parameter :: omegamn = 0.01_q              ! Minimum frequency in eV
    real(q), parameter :: pi = 3.1415926535897932384626433832795_q
    real(q), parameter :: ELECTRONVOLT_SI  = 1.602176487E-19_q
    real(q), parameter :: HARTREE_SI       = 4.35974394E-18_q
    real(q), parameter :: AUTOEV = HARTREE_SI/ELECTRONVOLT_SI
    real(q), parameter :: RYDTOEV = 13.6056980659_q
    real(q), parameter :: HBAR=6.58211928E-16_q ! in eV*sec
    real(q), parameter :: EVTOJ=1.60217733E-19_q ! conversion ev to Joule
    real(q), parameter :: ANGTOMT=1.E-10_q ! conevrsion angstroem to meter
    real(q), parameter :: ELEMASS=9.10938356E-31_q ! electron mass
    real(q), parameter :: BRTOANG=0.529177_q ! conversion bohr radius to angstroem
    real(q), parameter :: ELECT  = 1.602199E-19_q ! in coulomb=A*sec
    real(q), parameter :: KB=8.6173303E-5_q ! Boltzmann constant in eV/K
    real(q) :: omegamx,omegamn              ! Maximum frequency in eV
    real(q), allocatable :: freq_g(:)      ! Frequency grid
    real(q) :: de               ! init minimum energy equality threshold
    real(q) :: volume ! real space cell volume
    ! Gradient matrix elements ( defined and calculated in the main program )
!    complex(q), allocatable    :: psi_nabla_psi(:,:,:,:)
!
    real(q) :: dispers                    ! Gaussian width (delta function) in eV
    integer :: nbands,nkpoints
    integer :: IU1,IU2,OU1
    integer :: nelec
!
    logical :: intrab,FROM_VASP,FROM_CODE,n_bin
    CHARACTER (LEN=30) :: FNAME
!
end module
