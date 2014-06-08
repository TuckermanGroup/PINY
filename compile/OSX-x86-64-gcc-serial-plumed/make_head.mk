# binary names
EXE_NAME = piny-plumed
LIB_NAME = libpiny.so

# find PLUMED
PLUMED = $(shell plumed info --root)
ifeq ($(PLUMED),)
    $(error PLUMED not found.)
endif

# compilers and options
FC = gfortran
CC = gcc
DEFINE = -DNO_CFREE -DPLUMED
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE) -I $(PLUMED)src/wrapper
FFLAGS = -fPIC $(DEFINE)
LIBS = -L $(PLUMED)src/lib -lplumed -ldl

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
