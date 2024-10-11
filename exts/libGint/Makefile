# Compiler and Options
NVCC = nvcc
MPICPP = mpic++
GFORTRAN = gfortran

NVCC_C_OPTS = -rdc=true -g -std=c++14 -gencode arch=compute_70,code=sm_70 -lcudart -Xcompiler -fPIC -lcublas -Xcompiler -fopenmp -Xcompiler -O3 -Xcompiler -g
NVCC_D_OPTS = -arch=sm_70 -lgomp
MPICPP_EX_OPTS = -std=c++14 -Wall -fPIC -O3
GFORTRAN_OPTS = -O3 -g -fopenmp

LIS = -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include -L/usr/lib/x86_64-linux-gnu/o -lmpi -lcudart -lstdc++ -lcudadevrt -lnvToolsExt

# Directories
SRC_DIR = src
OBJ_DIR = obj
CP2K_LIB_DIR = /home/qjn24437/cp2k/lib/local_cuda/psmp/
CP2K_PKG_DIR = /home/qjn24437/cp2k/src/Gint/

# Objects
OBJECTS_1 = $(OBJ_DIR)/interface_libgint.o $(OBJ_DIR)/extern_functions.o $(OBJ_DIR)/plan.o $(OBJ_DIR)/UniqueArray.o $(OBJ_DIR)/libGint.o
OBJECTS_2 = $(OBJ_DIR)/libGint_unlinked.o $(OBJ_DIR)/util.o $(OBJ_DIR)/fgamma.o $(OBJ_DIR)/compute_Fm.o $(OBJ_DIR)/compute_VRR.o $(OBJ_DIR)/compute_ECO.o $(OBJ_DIR)/compute_HRR.o $(OBJ_DIR)/compute_SPH.o $(OBJ_DIR)/compute_TRA.o $(OBJ_DIR)/compute_KS.o

# Targets
all: libcp2kGint.a

# Compilation rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(NVCC) -c $(NVCC_C_OPTS) -x cu $< -o $@ $(LIS)

$(OBJ_DIR)/libGint_unlinked.o: $(SRC_DIR)/libGint.cpp
	$(NVCC) -c $(NVCC_C_OPTS) -x cu $< -o $@ $(LIS)

$(OBJ_DIR)/libGint.o: $(OBJ_DIR)/libGint_unlinked.o $(OBJECTS_2)
	$(NVCC) -dlink $(NVCC_D_OPTS) -o $@ $^ $(LIS)

#$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
#	$(MPICPP) -c -g -O3 $(MPICPP_EX_OPTS) $< -o $@ $(LIS) -fopenmp

$(OBJ_DIR)/interface_libgint.o: $(SRC_DIR)/interface_libgint.F90
	$(GFORTRAN) $(GFORTRAN_OPTS) -c $< -o $@ -lstdc++ -lcudadevrt -lcudart

# Static library creation
libcp2kGint.a: $(OBJECTS_1) $(OBJECTS_2)
	ar -rcs $@ $^

# File copying and packaging
install: libcp2kGint.a
	cp $(SRC_DIR)/*.cpp $(CP2K_PKG_DIR)
	cp libcp2kGint.a $(CP2K_LIB_DIR)
	cp libgint.mod $(CP2K_PKG_DIR)
	touch $(CP2K_PKG_DIR)/libGint_unlinked.cpp
	touch $(CP2K_PKG_DIR)/libgint.F
	echo "{\n   \"description\": \"interface for the calculation of 4 centres 2 electrons integrals\",\n      \"requires\": [\"../base\"],\n}" > $(CP2K_PKG_DIR)/PACKAGE
	rm -f /home/qjn24437/cp2k/obj/local_cuda/psmp/cp2k.o
#	cd /home/qjn24437/cp2k
#	make  -j 32 ARCH=local_cuda VERSION="psmp"
#	cd -

# Clean the build
clean:
	rm -f $(OBJ_DIR)/*.o libcp2kGint.a *.mod

.PHONY: all clean install
