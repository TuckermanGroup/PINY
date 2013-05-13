#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
# MAKEFILE FOR PI_MD
#
# This makefile can be used from anywhere NFSd to the directory.
# The objs defined in ../make_defs/obj_files.h, are placed in the same 
# directory as the Makefile. You may build the executable by typing ``make''.
# You may delete all of the object files by typing "make clean"
# You may delete some of the object files by typing "make clean_care"
# or "make clean_group". Compilation instructions are in 
# ../make_defs/make_*.h files and make_spec.h. Machine specific routines 
# are in spec_files.h.  Machine specific directory tree and compiler options 
# are given in make_head.h.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================


#==========================
# General unix instructions
#==========================
ECHO    = /bin/echo making
RM      = /bin/rm -f
CP      = /bin/cp
MV      = /bin/mv -f
TOUCH   = touch

#==========================
# Machine specific info
#==========================
include make_head.h

#==========================
# Compiler commands
#==========================
FOBJ       = $(FC) $(OPT) -I$(INCLUDES) -c -o $@
COBJ       = $(CC) $(OPT) $(CFLAGS) -I$(INCLUDES) -c -o $@
COBJ_FUSS  = $(CC) $(OPT) $(CFLAGS) -I$(INCLUDES) -c -o $@
COBJ_NOOPT = $(CC) $(CFLAGS)  -I$(INCLUDES) -c -o $@
COBJ_CARE  = $(CC) $(CFLAGS) $(OPT_CARE) -I$(INCLUDES) -c -o $@
COBJ_GRP   = $(CC) $(OPT_GRP) $(CFLAGS) -I$(INCLUDES) -c -o $@
LD         = $(FC) $(OPT) $(FFLAGS) 

#==========================
# Make the executable 
#==========================
include spec_files.h
include ../make_defs/obj_files.h

$(EXE): $(OBJS)
	$(ECHO) $(EXE)
	$(RM) $(EXE)
	$(LD)  $(OBJS) -o  $(EXE) $(LIBS)

#==========================
# Compilation instructions
#==========================
include ../make_defs/proto_files.h
include ../make_defs/make_analysis.h
include ../make_defs/make_analysis_cp.h
include ../make_defs/make_cp_energy.h
include ../make_defs/make_energy_inter.h
include ../make_defs/make_energy_intra.h
include ../make_defs/make_energy_intra_con.h
include ../make_defs/make_energy_pimd.h
include ../make_defs/make_friend.h
include ../make_defs/make_integrate_cp.h
include ../make_defs/make_integrate_md.h
include ../make_defs/make_integrate_pimd.h
include ../make_defs/make_integrate_cppimd.h
include ../make_defs/make_interface.h
include ../make_defs/make_interface_cp.h
include ../make_defs/make_interface_intra.h
include ../make_defs/make_interface_mol.h
include ../make_defs/make_main.h
include ../make_defs/make_main_cp.h
include ../make_defs/make_math.h
include ../make_defs/make_output.h
include ../make_defs/make_output_cp.h
include ../make_defs/make_communicate.h
include make_spec.h

#==========================
#remove all objs and libs
#==========================
clean       :
	/bin/rm -f $(OBJS) $(ULIBS) *.o

#==========================
#remove care objs and libs
#==========================
clean_care  :
	/bin/rm -f $(OBJS_CARE) $(ULIBS) *.o

#==========================
#remove group objs and libs
#==========================
clean_group :
	/bin/rm -f $(OBJS_GRP) $(ULIBS) *.o
