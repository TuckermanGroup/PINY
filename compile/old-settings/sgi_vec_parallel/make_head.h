#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
#
#  Machine specific compilation information
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================

#============================
# user specific directories
#============================
ROOT   = /home/dcode
CODE   = $(ROOT)/work_PARA_VER_3_JULY_03
DCODE  = $(ROOT)/work_PARA_VER_3_JULY_03
ECODE  = $(DCODE)/runable

#============================
# general directories
#============================
INCLUDES = $(DCODE)/include/sgi_vec_parallel
EXE      = $(ECODE)/piny_md_sgi_par.e

#============================
#SGI compiler:
#============================
FC       = f77
CC       = cc 
FFLAGS   = 
LIBS     = -lm -lcomplib.sgimath -lmalloc  -lmpi
#----------------------------
# Full opt
#----------------------------
OPT      = -O2 -mips4 -OPT:IEEE_arithmetic=2 \
           -TARG:platform=ip27 -TARG:processor=r10000
OPT_CARE = 
OPT_GRP  = -O2 -mips4 -OPT:IEEE_arithmetic=2 \
           -TARG:platform=ip27 -TARG:processor=r10000
CFLAGS   =
#----------------------------
# No opt
#----------------------------
#OPT      =
#OPT_CARE =
#OPT_GRP  =
#CFLAGS   = -fullwarn
#============================

