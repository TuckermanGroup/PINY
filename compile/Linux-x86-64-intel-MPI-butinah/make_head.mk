# build with:
# module load openmpi/intel/1.6.2

# binary names
EXE_NAME = piny-mpi
LIB_NAME = libpiny.so

# compilers and options
FC = mpif77
CC = mpicc
DEFINE =
OPT = -O3 -no-prec-div -xHost
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC -static-intel $(DEFINE)
FFLAGS = $(CFLAGS) -nofor-main
LIBS = $(LIB_PATH) $(MALLOC) -Wl,-rpath,/share/apps/intel/12.1.3.293/composer_xe_2011_sp1.9.293/compiler/lib/intel64

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux_parallel
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
