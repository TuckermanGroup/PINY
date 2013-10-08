# user specific directories
#--------------------------
ROOT   = /home/d/daniel/PINY
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(DCODE)

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny_md_pentium_normal
CMALLOC = 

# HP compiler
#--------------------------
FC = gfortran -DLINUX 
CC = gcc -DLINUX 
OPT = -O0 
OPT_CARE = -O0 
OPT_GRP = -O0
CFLAGS =  
FFLAGS =
LIBS =   $(LIB_PATH) $(MALLOC) -lm 


