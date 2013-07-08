#
# user specific directories
#

ROOT  = /home/om15/code/PINY/PINY
CODE  = $(ROOT)
DCODE = $(ROOT)
ECODE = $(DCODE)/runable


#
# general directories
#

INCLUDES = $(DCODE)/include/pentium
EXE      = $(ECODE)/piny-mpi
LIBPINY  = $(ECODE)/libpiny.so
CMALLOC =


#
# compilers and options
#

DEFINE = -DLINUX
FC = mpif77 $(DEFINE)
CC = mpicc $(DEFINE)
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC
FFLAGS = -fPIC
LIBS = $(LIB_PATH) $(MALLOC) -lm
