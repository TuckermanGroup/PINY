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
MY_ROOT = /u/gmartyna/PROJECTS
CODE    = $(MY_ROOT)/work_PARA_VER_3_MAR_2002
DCODE   = $(MY_ROOT)/work_PARA_VER_3_MAR_2002
ECODE   = $(CODE)/runable

#============================
# general directories
#============================
INCLUDES = $(CODE)/include/ibm
EXE      = $(ECODE)/pi_md_ibm.e
CMALLOC  = 

#============================
# IBM compiler 
#============================
FC     = xlf 
CC     = xlc
CFLAGS =
FFLAGS = 
LIBS   = -lm  -lessl
#----------------------------
#  No Opt
#----------------------------
OPT      = -O2
OPT_CARE = -O2
OPT_GRP  = -O2
#----------------------------
#  Full Opt
#----------------------------
#OPT      = -O3 -qstrict -qarch=pwr2
#OPT_CARE = -O3 -qstrict -qarch=pwr2
#OPT_GRP  = -O3 -qstrict -qarch=pwr2
#----------------------------
#  NOTE : 
#  If you are compiling for an IBM RISC workstation, 
#  take off the -qarch=pwr2 flag, since that is 
#  special for the SP2.
#============================

