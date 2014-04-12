# binary names
EXE_NAME = piny-plumed
LIB_NAME = libpiny.so

# compilers and options
PLUMED = /Users/andy/build/plumed-2.0.2
include $(PLUMED)/src/lib/Plumed.inc.shared
FC = gfortran
CC = gcc
DEFINE = -DNO_CFREE -DPLUMED
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE) -I $(PLUMED)/src/wrapper
FFLAGS = -fPIC $(DEFINE)
LIBS = $(PLUMED_LOAD) -ldl -Wl,-rpath,$(PLUMED)/src/lib

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
