


# the name of the executable file
#######################################################################################################
EXECUTABLE  := swim



# simulation options
#######################################################################################################

# GPU acceleration
SIMU_OPTION += -DGPU

# use soften length in acceleration, jerk, and potential calculation
SIMU_OPTION += -DSOFTEN

# use implicit 4-th order hermite scheme
#SIMU_OPTION += -DIMPLICIT_HERMITE

# improve the GPU performance if N is a multiple of BLOCK_SIZE
SIMU_OPTION += -DN_IS_MULTIPLE_OF_BS

# use the HYBRID scheme (Harfst, 2007) for the force calculation (otherwise the RING scheme is adopted)
#SIMU_OPTION += -DHYBRID_SCHEME

# declare the accumulation arrays of acceleration, jerk, and potential as double precision
#SIMU_OPTION += -DFLOAT8_ACC

# declare all floating-point variables as double precision (must turn on FLOAT8_ACC as well)
#SIMU_OPTION += -DFLOAT8

# shared time-step (otherwise the individual time-step is adopted)
#SIMU_OPTION += -DSHARED_TIMESTEP

# GPU architecture: FERMI/KEPLER/MAXWELL/PASCAL/VOLTA/TURING
SIMU_OPTION += -DGPU_ARCH=TURING

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP



# simulation parameters
#######################################################################################################

# size of CUDA thread block and grid
ifeq      "$(filter -DGPU_ARCH=FERMI,   $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
BLOCK_SIZE = 128
GRID_SIZE  = 56
else ifeq "$(filter -DGPU_ARCH=TURING,  $(SIMU_OPTION))" "-DGPU_ARCH=TURING"
BLOCK_SIZE = 128
GRID_SIZE  = 48
else
BLOCK_SIZE = 128
GRID_SIZE  = 32
endif

BLOCK_SIZE := $(strip $(BLOCK_SIZE))
GRID_SIZE  := $(strip $(GRID_SIZE))

SIMU_PARA = -DBLOCK_SIZE=$(BLOCK_SIZE) -DGRID_SIZE=$(GRID_SIZE)



# source files
#######################################################################################################

# Cuda source files (compiled with nvcc)
CUDA_FILE   := CUAPI.cu  CUCAL_Acc_Jerk.cu  CUCAL_Acc_Jerk_Split.cu  CUCAL_Pot.cu

# C/C++ source files (compiled with mpiCC)
CC_FILE     := Main.cpp  DataTransfer.cpp  Auxiliary.cpp  Init.cpp  Evolve.cpp  CPU_Acc_Jerk.cpp  CPU_Pot.cpp \
               ExternalForce.cpp



# Rules and targets
#######################################################################################################

NVCC = $(CUDA_PATH)/bin/nvcc
CXX  = $(MPI_PATH)/bin/mpicxx


# NTU-eureka (default: openmpi-intel)
CUDA_PATH := /software/cuda/default
MPI_PATH  := /software/openmpi/default

OBJ_PATH  := Object

vpath %.cu CUDA


LIB := -limf
#LIB :=

ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
   LIB += -L$(CUDA_PATH)/lib64 -lcudart
endif

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" "OPENMP"
   OPENMP := -fopenmp
endif


INCLUDE := -I../include -I$(MPI_PATH)/include

COMMONFLAG := $(INCLUDE) $(SIMU_OPTION) $(SIMU_PARA)

CXXWARN_FLAG := -w1
#CXXWARN_FLAG := -Wextra -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses \
                -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" ""
CXXWARN_FLAG += -Wno-unknown-pragmas -diag-disable 3180
endif

CXXFLAG  := $(CXXWARN_FLAG) $(COMMONFLAG) $(OPENMP) -O3 -g

NVCCFLAG := -Xcompiler $(COMMONFLAG) -O3 #--ptxas-options=-v #-D_DEBUG #-deviceemu

ifeq      "$(filter -DGPU_ARCH=FERMI,   $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
   NVCCFLAG += -gencode arch=compute_20,code=\"compute_20,sm_20\"
else ifeq "$(filter -DGPU_ARCH=KEPLER,  $(SIMU_OPTION))" "-DGPU_ARCH=KEPLER"
   NVCCFLAG += -gencode arch=compute_30,code=\"compute_30,sm_30\"
   NVCCFLAG += -gencode arch=compute_35,code=\"compute_35,sm_35\"
   NVCCFLAG += -gencode arch=compute_37,code=\"compute_37,sm_37\"
else ifeq "$(filter -DGPU_ARCH=MAXWELL, $(SIMU_OPTION))" "-DGPU_ARCH=MAXWELL"
   NVCCFLAG += -gencode arch=compute_50,code=\"compute_50,sm_50\"
   NVCCFLAG += -gencode arch=compute_52,code=\"compute_52,sm_52\"
else ifeq "$(filter -DGPU_ARCH=PASCAL,  $(SIMU_OPTION))" "-DGPU_ARCH=PASCAL"
   NVCCFLAG += -gencode arch=compute_60,code=\"compute_60,sm_60\"
   NVCCFLAG += -gencode arch=compute_61,code=\"compute_61,sm_61\"
else ifeq "$(filter -DGPU_ARCH=VOLTA,  $(SIMU_OPTION))" "-DGPU_ARCH=VOLTA"
   NVCCFLAG += -gencode arch=compute_70,code=\"compute_70,sm_70\"
else ifeq "$(filter -DGPU_ARCH=TURING,  $(SIMU_OPTION))" "-DGPU_ARCH=TURING"
   NVCCFLAG += -gencode arch=compute_75,code=\"compute_75,sm_75\"
else ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
   $(error unknown GPU_ARCH (please set it in the Makefile))
endif

NVCCFLAG += -prec-div=false -prec-sqrt=false -ftz=true


OBJ := $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(CC_FILE))

ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
   OBJ += $(patsubst %.cu, $(OBJ_PATH)/%.o, $(CUDA_FILE))
endif


# implicit rules
$(OBJ_PATH)/%.o : %.cpp
	$(CXX) $(CXXFLAG) -o $@ -c $<

$(OBJ_PATH)/%.o : %.cu
	$(NVCC) $(NVCCFLAG) -o $@ -c $<


# link all object files
$(EXECUTABLE) : $(OBJ)
	$(CXX) -o $@ $^ $(LIB) $(OPENMP)
	cp $(EXECUTABLE) ../bin/
	rm ./*.linkinfo -f


# clean
# -------------------------------------------------------------------------------
.PHONY: clean
clean :
	@rm -f $(OBJ_PATH)/*
	@rm -f $(EXECUTABLE)
	@rm -f ./*.linkinfo
