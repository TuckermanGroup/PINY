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

INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny
LIBPINY  = $(ECODE)/libpiny.so
CMALLOC =


#
# compilers and options
#

DEFINE = -DLINUX
FC = gfortran $(DEFINE)
CC = gcc $(DEFINE)
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC
FFLAGS = -fPIC
LIBS = $(LIB_PATH) $(MALLOC) -lm
