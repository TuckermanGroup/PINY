# user specific directories
#--------------------------
ROOT   = /home/dcode
CODE   = $(ROOT)/work_PARA_VER_3_JULY_03
DCODE  = $(ROOT)/work_PARA_VER_3_JULY_03
ECODE  = $(DCODE)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny_md_pentium_test
CMALLOC = 

# HP compiler
#--------------------------
FC = f77  -DLINUX -fno-second-underscore
CC = gcc -DLINUX -fno-second-underscore
OPT = -O0 
OPT_CARE = -O0 
OPT_GRP = -O0
CFLAGS =  
FFLAGS = -nofor_main
LIBS =   $(LIB_PATH) $(MALLOC) -lm 


