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
ROOT   = /N/fs19/gmartyna/Solar
CODE   = $(ROOT)/PARA_VER_3_NOV_11_2000
DCODE  = $(ROOT)/PARA_VER_3_NOV_11_2000
ECODE  = $(DCODE)/runable

#============================
# general directories
#============================
INCLUDES = $(DCODE)/include/sun_parallel
EXE      = $(ECODE)/pi_md_sun_par.e

#============================
#SGI compiler:
#============================
FC       = mpf77
CC       = mpcc 
CFLAGS   =
FFLAGS   =  
LIBS     = -lm -lmpi -xlic_lib=sunperf_mt,mtsk,thread,c,fui,fsu,sunmath
#----------------------------
# Full opt
#----------------------------
OPT      =  -xO2
OPT_CARE =  -xO2
OPT_GRP  =  -xO2
#----------------------------
# No opt
#----------------------------
#OPT      =
#OPT_CARE =
#OPT_GRP  =
#CFLAGS   = -fullwarn
#============================

