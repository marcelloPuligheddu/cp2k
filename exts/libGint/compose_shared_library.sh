set -x

rm obj/* libcp2kGint.a *.mod

NVCC_C_OPTS="-rdc=true -std=c++14 -gencode arch=compute_70,code=sm_70 -lcudart -Xcompiler -fPIC -lcublas"
NVCC_D_OPTS="-arch=sm_70 "

MPICPP_EX_OPTS="-std=c++14 -Wall -fPIC " 

LIS="-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -L/usr/lib/x86_64-linux-gnu/o -lmpi -lcudart -lstdc -lcudadevrt -lnvToolsExt"

#for NAME in util fgamma t_c_g0_n compute_Fm compute_VRR compute_HRR compute_SPH compute_TRA compute_KS 

# almost reasonable -c compile
nvcc -c ${NVCC_C_OPTS} -x cu src/util.cpp -o obj/util.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/fgamma.cpp -o obj/fgamma.o ${LIS}
#nvcc -c ${NVCC_C_OPTS} -x cu src/t_c_g0_n.cpp -o obj/t_c_g0_n.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_Fm.cpp -o obj/compute_Fm.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_VRR.cpp -o obj/compute_VRR.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_HRR.cpp -o obj/compute_HRR.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_SPH.cpp -o obj/compute_SPH.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_TRA.cpp -o obj/compute_TRA.o ${LIS}
nvcc -c ${NVCC_C_OPTS} -x cu src/compute_KS.cpp -o obj/compute_KS.o ${LIS}

# ? compile with -c 
nvcc -c ${NVCC_C_OPTS} -x cu src/libGint.cpp -o obj/libGint_unlinked.o ${LIS}

# ?? recompile with -dlink
nvcc -dlink ${NVCC_D_OPTS} -o obj/libGint.o obj/libGint_unlinked.o obj/compute_*.o obj/fgamma.o obj/util.o ${LIS} #obj/t_c_g0_n.o

## Single step for pure cpp files
mpic++ -c -std=c++14  -Wall -fPIC src/plan.cpp -o obj/plan.o -lcudart -lstdc
mpic++ -c -std=c++14  -Wall -fPIC src/UniqueArray.cpp -o obj/UniqueArray.o -lcudart -lstdc
mpic++ -c src/extern_functions.cpp -o obj/extern_functions.o ${MPICPP_EX_OPTS} ${LIS}

gfortran -c src/interface_libgint.F90 -o obj/interface_libgint.o -lstdc -lcudadevrt -lcudart


OBJECTS_1="obj/interface_libgint.o obj/extern_functions.o obj/plan.o obj/UniqueArray.o obj/libGint.o"
OBJECTS_2="obj/libGint_unlinked.o  obj/util.o obj/compute_*.o obj/fgamma.o"

ar -rcs libcp2kGint.a ${OBJECTS_1} ${OBJECTS_2}

#mpifort fortran_try.F90 -L. -lcp2kGint -lstdc++ -lcudart -lcublas -lblas

CP2K_LIB_DIR="/home/qjn24437/cp2k/lib/local_cuda/psmp/"
CP2K_SRC_DIR="/home/qjn24437/cp2k/src/Gint/"

cp src/*.cpp ${CP2K_SRC_DIR}
cp libcp2kGint.a ${CP2K_LIB_DIR}
cp libgint.mod ${CP2K_SRC_DIR}
touch ${CP2K_SRC_DIR}libGint_unlinked.cpp
touch ${CP2K_SRC_DIR}libgint.F

## remove the cp2k.o file
## recompile cp2k

set +x


