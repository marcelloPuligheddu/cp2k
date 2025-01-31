
## Folder File MPI OMP ##
tests=(
# "H;H1.inp;1;1"
# "H;H2.inp;1;1"
# "H;H3.inp;1;1"
# "H2O;H2O-HFX-1.inp;1;1"
# "H2O;H2O-HFX-1.inp;1;2"
# "H2O;H2O-HFX-1.inp;2;1"
# "H2O;H2O-HFX-6.inp;1;1"
# "H2O;H2O-HFX-6.inp;2;2"
# "Silicon;silicon8.inp;2;2" 
# "Silicon;silicon8.inp;4;4" 
# "H2O_PBE;water_pbe0_from_pbe.inp;4;4"
# "MgO;MgO_pbe0_energy_smaller.inp;4;4"
 "Ag2;SilverDimer.inp;4;4"
 "Ethybenzene;Ethybenzene.inp;4;4"
 "KMnF3;KMnF3.inp;4;4"
)

echo $tests

#TODO: read from command line
EXE=../../../exe/local_cuda/cp2k.psmp
TOL=0.0001



for test_values in "${tests[@]}"; do

   IFS=";"

   read -r -a arr <<< "${test_values}"

   echo "${arr} ${test_values}"

   FIRST=true
   FOL="${arr[0]}"
   INP="${arr[1]}"
   MPI="${arr[2]}"
   OMP="${arr[3]}"

   CPU_INPUT=${INP}
   GPU_INPUT=${CPU_INPUT}.gpu
   LOG=${CPU_INPUT}".log"
   TMP=${CPU_INPUT}".tmp"
   echo "Testing: ${FOL}/${CPU_INPUT} "

   cd ${FOL}
   # Create the input file using libGint from the regular input file
   sed -e 's\libint\libgint\' ${CPU_INPUT} > ${GPU_INPUT}

   echo "MPI OMP INP || Tot Tot-Ref "


   for INP in ${CPU_INPUT} ${GPU_INPUT}; do
      export OMP_NUM_THREADS=${OMP}
      # Get the Total Energy from the run
      TE=$( mpirun -np ${MPI} ${EXE} -i ${INP} 2>&1 | tee ${TMP} | grep "Total energy:" | awk '{print $NF}' )
      # If this is the first run, use its value oa the ref value
      if $FIRST ; then
        REF=${TE}
        FIRST=false
      fi
      # Compute the difference between the current energy and the ref value
      DE=$( echo "(${TE} - ${REF})" | bc -l )
      # Hacky way to take the abs value of a number
      DE=${DE#-}
      # Exit if the eerngy is too far from the reference
      if (( $(echo "${DE} > ${TOL}" |bc -l) )); then
         echo "$MPI $OMP ${INP} || ${TE} ${DE} FAILED. Check ${TMP}"
         if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
            # sourced
            return 1
         else
            # not sourced
            exit 1
         fi
      fi
      # Otherwise print some info and pass
      echo "$MPI $OMP ${INP} || ${TE} ${DE} PASS "
      wait
      cat ${TMP} >> ${LOG}
      rm ${TMP}
   done

   echo ""
   rm ${GPU_INPUT}
   cd ..
done

