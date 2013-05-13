# user specific directories
#--------------------------
ROOT   = /home/dcode/work_PARA_VER_3_JULY_03
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(DCODE)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium
EXE      = $(ECODE)/pi_md_pentium_par
CMALLOC = 

# HP compiler
#--------------------------
FC  = mpif77  -fno-second-underscore
CC  = mpicc
OPT = -O2 
OPT_CARE = -O2 
OPT_GRP = -O2
CFLAGS =  
FFLAGS = 
LIBS =   $(LIB_PATH) $(MALLOC) -lm 


