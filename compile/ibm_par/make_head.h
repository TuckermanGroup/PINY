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
ROOT   = /home/yarne/CODE
CODE   = $(ROOT)/work_PARA_VER_3
DCODE  = $(ROOT)/work_PARA_VER_3
ECODE  = $(DCODE)/runable

#============================
# general directories
#============================
INCLUDES = $(CODE)/include/ibm_par
EXE      = $(ECODE)/pi_md_ibm.e
CMALLOC  = 

#============================
# IBM compiler 
#============================
FC     = mpxlf  
CC     = mpcc
CFLAGS =
FFLAGS = 
# IF you need to allocate > 256 Mbytes of real memory you need the flags below 
#   these flags are for  2 gigabyte allocation
#   see man pages for xlf 
#CFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x10000000  
#FFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x10000000  
LIBS   = -lm  -lessl
#----------------------------
#  No Opt
#----------------------------
OPT      = -O2 -qmaxmem=-1
OPT_CARE = -O2 -qmaxmem=-1
OPT_GRP  = -O2 -qmaxmem=-1
#----------------------------
#  Full Opt
#----------------------------
#OPT      = -O3 -qstrict -qarch=pwr2 -qmaxmem=-1
#OPT_CARE = -O3 -qstrict -qarch=pwr2 -qmaxmem=-1
#OPT_GRP  = -O3 -qstrict -qarch=pwr2 -qmaxmem=-1
#----------------------------
#  NOTE : 
#  If you are compiling for an IBM RISC workstation, 
#  take off the -qarch=pwr2 flag, since that is 
#  special for the SP2.
#============================

