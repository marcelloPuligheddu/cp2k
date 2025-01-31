EXE=../../../exe/local_cuda/cp2k.psmp
export OMP_NUM_THREADS=4
mpirun -np 4 ${EXE} -i water_pbed3.inp | tee prepare.sh.log 
cp WATER-RESTART.wfn WATER-RESTART-GGA.wfn 
