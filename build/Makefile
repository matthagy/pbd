
all: install

# # # # # # # # #
# Configuration #
# # # # # # # # #
# 1 to enable, defaults to enviromental variables when set

#Generate optimized executables
OPTIMIZE ?= 1

#Generate executables that can be profiled with gprof
PROFILE ?= 0

#Building for Mac OS X
MACOSX ?= 1

#1 for MVAPICH2 else LAM
USING_MPI_MVAPICH2 ?= 0

#Use Intel Math Kernel Library (MKL) for random number generation
MKL_RANDOM ?= 0

#Parallelize outer force evaluation loop with OpenMP
#require less efficient data structures and is really
#only efficient when compiled with icc as there is special
#support for intel hyper-threading
OMP_PARALLELIZE_FORCES ?= 0

#Parallelize outer integration loop with OpenMP
OMP_PARALLELIZE_INTEGRATION ?= 0

#Simplify force interpolation to use a single point as opposed
#to linear interpolation between two points.  Can substantially
#speed up executation at the cost of requiring larger interpolation
#table and possibly adding artifiacts to the force field
SINGLE_POINT_INTERPOLATION ?= 0

#When ran with multiple threads and have junctioned cells, we evaluate 
#internal forces while exchanging positions with other cells.
OMP_CONCURRENT_FORCE_EVALUATION ?= 0


# # # # # # # # # # # # #
# C-Preprocessor Macros #
# # # # # # # # # # # # #

MACRO_DEFINES = 

ifeq ($(MACOSX), 1)
  MACRO_DEFINES += MACOSX
endif

ifeq ($(USING_MPI_MVAPICH2), 1)
  MACRO_DEFINES += USING_MPI_MVAPICH2
else
  MACRO_DEFINES += USING_MPI_LAM
endif

ifeq ($(OMP_PARALLELIZE_FORCES), 1)
  MACRO_DEFINES += OMP_PARALLELIZE_FORCES
endif

ifeq ($(OMP_PARALLELIZE_INTEGRATION), 1)
  MACRO_DEFINES += OMP_PARALLELIZE_INTEGRATION
endif

ifeq ($(SINGLE_POINT_INTERPOLATION), 1)
  MACRO_DEFINES += SINGLE_POINT_INTERPOLATION
endif

ifeq ($(MKL_RANDOM), 1)
  MACRO_DEFINES += MKL_RANDOM
endif

ifeq ($(OMP_CONCURRENT_FORCE_EVALUATION), 1)
  MACRO_DEFINES += OMP_CONCURRENT_FORCE_EVALUATION
endif


# # # # # # #
# Compiler  #
# # # # # # #

CC = gcc #MPI C Compiler (icc internally)

COMMON_FLAGS =
COMMON_FLAGS += -Wall #Warn All
COMMON_FLAGS += -std=c99 #Use c99 standard
COMMON_FLAGS += -fstrict-aliasing #Enable optimization requiring strict pointer aliasing
#COMMON_FLAGS += -fopenmp #OpenMP and implicit parallelization

CC_EXEC = $(CC) $(COMMON_FLAGS) $(MACRO_DEFINES:%=-D%)

OPTIMIZE_FLAGS =
OPTIMIZE_FLAGS += -DNDEBUG #No debug (disablle assert(...))
OPTIMIZE_FLAGS += -fomit-frame-pointer #Free up EBP register for general operations
OPTIMIZE_FLAGS += -O2 #Intensive optimization
#OPTIMIZE_FLAGS += -march=native #Optimize specifically for this machine
OPTIMIZE_FLAGS += -ftree-vectorize #Vectorize at gimple level
#OPTIMIZE_FLAGS += -m32 #32 bit output
#OPTIMIZE_FLAGS += -msse3 -mfpmath=sse,387 #Enable SSE3 but allow 387 floating pointer operations when more efficient

ifeq ($(OPTIMIZE),1)
  CC_EXEC += $(OPTIMIZE_FLAGS)
endif

PROFILE_FLAGS = -g -p #Include symbols and profiling information for gprof

ifeq ($(PROFILE),1)
  CC_EXEC += $(PROFILE_FLAGS)
endif


# # # # # # #
# Libraries #
# # # # # # #

LIBRARIES = 
LIBRARIES += m #Standard math library
LIBRARIES += pthread #Threading needed for OpenMP
LINKING_OPTIONS = $(LIBRARIES:%=-l%)

# MKL_PATH = /opt/intel/mkl/10.0.3.020/lib/em64t/
# MKL_LIBS = 
# MKL_LIBS += $(MKL_PATH)/libmkl_intel_ilp64.so
# MKL_LIBS += $(MKL_PATH)/libmkl_intel_thread.so
# MKL_LIBS += $(MKL_PATH)/libmkl_core.so

# MKL_LINKING_OPTIONS =  -Wl,--start-group $(MKL_LIBS) -Wl,--end-group $(LIBRARIES:%=-l%)

# ifeq ($(MKL_RANDOM),1)
#    LINKING_OPTIONS += $(MKL_LINKING_OPTIONS)
# endif

# # # # # # # # # #
# Build Commands  #
# # # # # # # # # #

BUILD_OBJ = $(CC_EXEC) -c -I/home/rhmac/hagy/opt/lam-7.1.4//include -D_REENTRANT 
BUILD_ASM = $(CC_EXEC) -S
BUILD_EXC = $(CC_EXEC) $(LINKING_OPTIONS) -Wl,-u -Wl,_lam_darwin_malloc_linker_hack -Wl,-multiply_defined,suppress -Wl,-force_flat_namespace -Wl,-flat_namespace -Wl,-multiply_defined -Wl,suppress -L/home/rhmac/hagy/opt/lam-7.1.4//lib -llammpio -lmpi -llam -lutil -ldl                                                                                                  


# # # # # # # # # #
# cex Executable  #
# # # # # # # # # #

HEADERS = 
OBJECTS = 

#Platform compatibility
HEADERS += compat.h

#Numerical constants
HEADERS += constants.h

#Compiler options/optimizations
HEADERS += opt.h

#Basic debugging
HEADERS += debug.h
OBJECTS += debug.o

#Memory management utilties
HEADERS += mem.h
OBJECTS += mem.o

#3D vector arithmetic
HEADERS += vector.h

#Generic dynamically-sized arrays
HEADERS += array.h array-inline.h
OBJECTS += array.o

#Random number generation
HEADERS += random.h
OBJECTS += random.o

#Messaging between C simulation process and Python control process
HEADERS += msg.h msg-inline.h
OBJECTS += msg.o

#MPI-based communication between simulation processes
HEADERS += comm.h comm-inline.h
OBJECTS += comm.o

#Periodic bondary condition helpers
HEADERS += xperiodic.h periodic.h
OBJECTS += periodic.o

#Cell decomposition of Cartesian space
HEADERS += cells.h
OBJECTS += cells.o

#Brownian dynamics simulation algorithms
HEADERS += bd.h
OBJECTS += bd.o

#Initialization routines
HEADERS += init.h
OBJECTS += init.o

#Main executable
OBJECTS += main.o


# # # # # # # #
#Build rules  #
# # # # # # # #

#Executable
cex: $(OBJECTS)
	$(BUILD_EXC) $^ -o $@

testarray: testarray.o array.o mem.o debug.o
	$(BUILD_EXC) $^ -o $@

testmsg: testmsg.o msg.o array.o mem.o debug.o
	$(BUILD_EXC) $^ -o $@

testrandom: testrandom.o random.o array.o mem.o debug.o
	$(BUILD_EXC) $^ -o $@

#Generic object build
COMMON_DEPS = Makefile $(HEADERS:%=../src/%)

%.o: ../src/%.c $(COMMON_DEPS)
	$(BUILD_OBJ) $< -o $@

#Special rules for configurable objects
random.o: ../src/random.c ../src/random-mt19937ar.c ../src/random-mkl.c $(COMMON_DEPS)
	$(BUILD_OBJ) ../src/random.c -o $@

bd.o: ../src/bd.c ../src/eval-forces-simple.c ../src/eval-forces-openmp.c $(COMMON_DEPS)
	$(BUILD_OBJ) ../src/bd.c -o $@

#Assemblies for debugging
%.S: ../src/%.c $(COMMON_DEPS)
	$(BUILD_ASM) $< -o $@

clean:
	rm -rf *.o cex testmsg testarray testrandom

install: cex testarray testmsg testrandom
	cp $^ ../bin
