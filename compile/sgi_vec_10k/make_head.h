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
ROOT   = /home/dcode/work_PARA_VER_3_JULY_03
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(CODE)/runable

#============================
# general directories
#============================
INCLUDES = $(CODE)/include/sgi_vec
EXE      = $(ECODE)/piny_md_sgi_ser.e

#============================
#SGI compiler:
#============================
FC       = f77
CC       = cc 
FFLAGS   = 
LIBS     = -lm -lcomplib.sgimath -lmalloc
#----------------------------
# Full opt
#----------------------------
OPT      = -O2 -mips4 -64
OPT_CARE = -O2 -mips4 -64
OPT_GRP  = -O2 -mips4 -64
CFLAGS   = -O2 -mips4 -64
#----------------------------
# Full opt
#----------------------------
#OPT      = 
#OPT_CARE = 
#OPT_GRP  = 
#CFLAGS   =
#----------------------------
# No opt
#----------------------------
#OPT      =
#OPT_CARE =
#OPT_GRP  =
#CFLAGS   = -fullwarn
#============================

