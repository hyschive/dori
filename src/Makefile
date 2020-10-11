


# the name of the executable file
#######################################################################################################
EXECUTABLE  := Swim



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

# enable the performance optimization in Fermi GPUs
SIMU_OPTION += -DFERMI

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP



# simulation parameters
#######################################################################################################

ifeq "$(findstring FERMI, $(SIMU_OPTION))" "FERMI"
BLOCK_SIZE = 128   # size of CUDA thread block (for Fermi GPUs)
GRID_SIZE  = 56    # size of CUDA thread grid  (for Fermi GPUs)
else
BLOCK_SIZE = 128   # size of CUDA thread block
GRID_SIZE  = 32    # size of CUDA thread grid
endif

BLOCK_SIZE         := $(strip $(BLOCK_SIZE))
GRID_SIZE          := $(strip $(GRID_SIZE))

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

NVCC := nvcc
CXX  := mpiCC


# hulk (openmpi-intel)
CUDA_TOOLKIT_PATH := /opt/gpu/cuda/default
CUDA_SDK_PATH     := /opt/gpu/cuda_sdk/default/C
MPI_PATH          := /opt/mpi/openmpi/1.6.0-intel

#CUDA_TOOLKIT_PATH := /usr/local/cuda
#CUDA_SDK_PATH     := /opt/devel/nvidia/C
#MPI_PATH          := /opt/mpi/openmpi_intel

#CUDA_TOOLKIT_PATH := /usr/local/cuda
#CUDA_SDK_PATH     := /opt/devel/nvidia/C
#MPI_PATH          := /opt/mpi/openmpi

#CUDA_TOOLKIT_PATH := /usr/local/cuda
#CUDA_SDK_PATH     := /opt/cuda_sdk.3.1/C
#MPI_PATH          := /usr/local

#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/3.1
#CUDA_SDK_PATH     := /global/scratch/sd/hsiyu/cuda_sdk.3.1/C
#MPI_PATH          := /usr/common/usg/openmpi/1.4.2/gnu

OBJ_PATH          := Object

vpath %.cu CUDA


LIB := -limf
#LIB :=

ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
#  LIB := -L$(CUDA_SDK_PATH)/lib -lcutil
   LIB += -L$(CUDA_SDK_PATH)/lib -lcutil_x86_64
   LIB += -L$(CUDA_TOOLKIT_PATH)/lib64 -lcudart
endif

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" "OPENMP"
   OPENMP := -openmp
endif


INCLUDE := -I../include -I$(CUDA_SDK_PATH)/common/inc -I$(MPI_PATH)/include

COMMONFLAG := $(INCLUDE) $(SIMU_OPTION) $(SIMU_PARA)

CXXWARN_FLAG := -w1
#CXXWARN_FLAG := -Wextra -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses \
                -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" ""
CXXWARN_FLAG += -Wno-unknown-pragmas
endif

CXXFLAG  := $(CXXWARN_FLAG) $(COMMONFLAG) $(OPENMP) -O3 -xSSE4.2 -fp-model precise #-debug -g

NVCCFLAG := -Xcompiler $(COMMONFLAG) -O3 #--ptxas-options=-v #-D_DEBUG #-deviceemu

NVCCFLAG += -gencode arch=compute_20,code=\"compute_20,sm_20\"
ifeq "$(findstring FERMI, $(SIMU_OPTION))" ""
   NVCCFLAG += -gencode arch=compute_13,code=\"compute_13,sm_13\"
   ifeq "$(findstring FLOAT8_ACC, $(SIMU_OPTION))" ""
      NVCCFLAG += -gencode arch=compute_10,code=\"compute_10,sm_10\"
   endif
endif

ifeq "$(findstring FERMI, $(SIMU_OPTION))" "FERMI"
   NVCCFLAG += -prec-div=false -prec-sqrt=false -ftz=true
endif


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
	cp $(EXECUTABLE) ../bin/run/
	rm ./*.linkinfo -f

clean :
	rm -f $(OBJ_PATH)/*
	rm -f $(EXECUTABLE)
	rm ./*.linkinfo -f



