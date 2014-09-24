# build with:
# module load openmpi/intel/1.6.2

# binary names
EXE_NAME = piny-mpi-plumed
LIB_NAME = libpiny.so

# find PLUMED
PLUMED = $(shell plumed info --root)
ifeq ($(PLUMED),)
    $(error PLUMED not found.)
endif

# compilers and options
FC = mpif77
CC = mpicc
DEFINE = -DPLUMED
OPT = -O3 -no-prec-div -xHost
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC -static-intel $(DEFINE) -I $(PLUMED)/src/wrapper/
FFLAGS = $(CFLAGS) -nofor-main
LIBS = -L $(PLUMED)/src/lib -lplumed -ldl -Wl,-rpath,$(PLUMED)/src/lib -Wl,-rpath,/share/apps/intel/12.1.3.293/composer_xe_2011_sp1.9.293/compiler/lib/intel64 -Wl,-rpath,/home/om15/opt/libmatheval-1.1.11/lib

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux_parallel
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
