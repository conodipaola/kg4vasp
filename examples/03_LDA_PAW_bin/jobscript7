#!/bin/bash --login
#PBS -N nabla            
#PBS -l select=1        
#PBS -l walltime=01:00:00
# Replace this with your budget code
#PBS -A e606

module swap PrgEnv-cray PrgEnv-intel
module add fftw
module add craype-hugepages2M

# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

# Set the number of threads to 1
#   This prevents any system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1

cp INCAR_SCF INCAR

  aprun -n 1 -N 1  /home/e606/e606/cdpa71/COMPILED/vasp.5.4.4_intel_mkl_Cubic_one_third_RELAX_wan90_300_nabla_new/bin/vasp_std_bin > vasp_scf

cp INCAR_nabla INCAR
mv OUTCAR OUTCAR_SCF
mv vasprun.xml vasprun.xml_SCF

aprun -n 1 -N 1  /home/e606/e606/cdpa71/COMPILED/vasp.5.4.4_intel_mkl_Cubic_one_third_RELAX_wan90_300_nabla_new/bin/vasp_std_bin > vasp_nabla



