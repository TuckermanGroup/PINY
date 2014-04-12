# binary names
EXE_NAME = piny-plumed
LIB_NAME = libpiny.so

# compilers and options
PLUMED = /Users/andy/opt/plumed-2.0.2
FC = gfortran
CC = gcc
#DEFINE = -DNO_CFREE -DPLUMED -DPLUMED_DEBUG
DEFINE = -DNO_CFREE -DPLUMED
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE) -I $(PLUMED)/include/Plumed/wrapper
FFLAGS = -fPIC $(DEFINE)
LIBS = -L $(PLUMED)/lib -lplumed -ldl

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
