set -x
module load fftw/3.3.10 fftw_mpi/3.3.10 gcc/base openmpi/5.0.7-ucc1.3.0-ucx1.18.0 hipfort/6.4.0

export OMP_NUM_THREADS=1
export HSA_XNACK=1

#mpirun -np 8 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.cpu.inp 2>&1 | tee 8_mpi_1_omp_libint
#wait
mpirun -np 3 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 3_mpi_1_omp_libgint_1
mpirun -np 6 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 6_mpi_1_omp_libgint_1
mpirun -np 2 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 2_mpi_1_omp_libgint_1
mpirun -np 1 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 1_mpi_1_omp_libgint_1
wait

mpirun -np 3 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 3_mpi_1_omp_libgint_2
wait

mpirun -np 6 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 6_mpi_1_omp_libgint_2
wait

mpirun -np 3 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 3_mpi_1_omp_libgint_3
wait

mpirun -np 6 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 6_mpi_1_omp_libgint_3
wait



#mpirun -np 4 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.cpu.inp 2>&1 | tee 4_mpi_1_omp_libint
#wait
#mpirun -np 4 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 4_mpi_1_omp_libgint
#wait

#mpirun -np 16 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.cpu.inp 2>&1 | tee 16_mpi_1_omp_libint
#wait
mpirun -np 16 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.gpu.inp 2>&1 | tee 16_mpi_1_omp_libgint
wait

export OMP_NUM_THREADS=2
mpirun -np 8 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.cpu.inp 2>&1 | tee 8_mpi_2_omp_libint
wait
export OMP_NUM_THREADS=4
mpirun -np 8 ../../../exe/local_hip/cp2k.psmp -i H2O-HFX-32.cpu.inp 2>&1 | tee 8_mpi_4_omp_libint
wait


export OMP_NUM_THREADS=1
set +x

