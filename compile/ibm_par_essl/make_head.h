# user specific directories
#--------------------------
ROOT   = ${HOME}/code
CODE   = $(ROOT)/work_PARA_VER_3_JULY_03
DCODE  = $(ROOT)/work_PARA_VER_3_JULY_03
ECODE  = $(DCODE)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/ibm_ppc64_par
EXE      = $(ECODE)/piny_md_par_ibm_essl
CMALLOC = 

# HP compiler
#--------------------------
FC = /usr/local/mpich/1.2.6..0.94/mx/ppc64/smp/ibmcmp64/ssh/bin/mpif77  
CC = /usr/local/mpich/1.2.6..0.94/mx/ppc64/smp/ibmcmp64/ssh/bin/mpicc  
FC = mpif77 -qextname
CC = mpicc
OPT =       -O3
OPT_CARE =  -O0 
OPT_GRP =   -O3
CFLAGS = -qmaxmem=-1 -qarch=ppc970 -qtune=ppc970 -qhot=simd -DFORTRANUNDERSCORE #-DWRITE_DENSITY 
FFLAGS = -qmaxmem=-1 -qarch=ppc970 -qtune=ppc970 -qnosave 
LIBS = -lessl
LIBSP = 


