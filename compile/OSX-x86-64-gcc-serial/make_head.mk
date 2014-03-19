# binary names
EXE_NAME = piny
LIB_NAME = libpiny.so

# compilers and options
FC = gfortran
CC = gcc
DEFINE = -DLINUX -DNO_CFREE
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE)
FFLAGS = -fPIC $(DEFINE)
LIBS = -lm

INCLUDES = $(CODE)/include/pentium_nopar
CODE = $(realpath ../..)
EXE = $(CODE)/bin/$(EXE_NAME)
LIBPINY = $(CODE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
