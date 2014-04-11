# binary names
EXE_NAME = piny
LIB_NAME = libpiny.so

# compilers and options
FC = gfortran
CC = gcc
DEFINE =
OPT = -O2
OPT_CARE = -O2
OPT_GRP = -O2
CFLAGS = -fPIC $(DEFINE)
FFLAGS = -fPIC $(DEFINE)
LIBS = -lm

BASE = $(realpath ../..)
CODE = $(BASE)/src
INCLUDES = $(BASE)/include/linux
EXE = $(BASE)/bin/$(EXE_NAME)
LIBPINY = $(BASE)/lib/$(LIB_NAME)
SPEC_FILES = math_generic.o
