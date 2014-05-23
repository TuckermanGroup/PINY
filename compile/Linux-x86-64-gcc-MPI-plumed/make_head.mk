# binary names
EXE_NAME = piny-plumed-mpi
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
#DEFINE =
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE) -I $(PLUMED)/src/wrapper -g
FFLAGS = -fPIC $(DEFINE)
LIBS = -L $(PLUMED)src/lib -lplumed -ldl -Wl,-rpath,$(PLUMED)src/lib

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux_parallel
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
