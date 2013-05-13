#=================================================================
#             Communicate Files
#==================================================================



#==================================================================
communicate_wrappers.o : $(STANDARD) $(DEFINES) \
                         $(COMM_WRAP)  \
                         $(CODE)/communicate/communicate_wrappers.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/communicate_wrappers.c

#-----------------------------------------------------------------
control_group_communicators.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_PAR) $(COMM_ENT) $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/control_group_communicators.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/control_group_communicators.c

#------------------------------------------------------------------
com_interface.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_ENT) $(COMM_WRAP) \
                         $(COMM_LOC) \
                         $(CODE)/communicate/com_interface.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/com_interface.c

#------------------------------------------------------------------
comm_class_info.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CP) $(TYP_CLASS) \
                         $(TYP_BND) $(COMM_WRAP) $(COMM_LOC) $(MATH) \
                         $(CODE)/communicate/comm_class_info.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_class_info.c

#------------------------------------------------------------------
comm_bond_info.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/comm_bond_info.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_bond_info.c

#------------------------------------------------------------------
comm_bond_data.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/comm_bond_data.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_bond_data.c

#------------------------------------------------------------------
comm_cp_info.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CLASS) $(TYP_CP) \
                         $(TYP_PAR) $(COMM_WRAP) $(COMM_LOC) $(MATH) \
                         $(CODE)/communicate/comm_cp_info.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_cp_info.c

#------------------------------------------------------------------
comm_class_data.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/comm_class_data.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_class_data.c

#------------------------------------------------------------------
comm_class_list.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/comm_class_list.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_class_list.c

#------------------------------------------------------------------
comm_cp_data.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) $(FRND_ENT) \
                         $(CODE)/communicate/comm_cp_data.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/comm_cp_data.c

#------------------------------------------------------------------
communicate_options.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS)  \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
                         $(CODE)/communicate/communicate_options.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/communicate_options.c

#------------------------------------------------------------------
mall_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT)\
                         $(TYP_CP) $(COMM_LOC)  $(FRND_ENT) $(PARSE_ENT) \
                         $(CODE)/communicate/mall_class.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/mall_class.c
#------------------------------------------------------------------
mall_bond.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT)\
                         $(TYP_CP) $(COMM_LOC) $(FRND_ENT)  $(PARSE_ENT) \
                         $(CODE)/communicate/mall_bond.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/mall_bond.c
#------------------------------------------------------------------
communicate_test_energy_pimd.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(TYP_CP)  $(ENR_CTRL_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/control/communicate_test_energy_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/communicate_test_energy_pimd.c
#------------------------------------------------------------------
communicate_utils_pimd.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INT_PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/integrate/pimd/communicate_utils_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/pimd/communicate_utils_pimd.c
#------------------------------------------------------------------
simpavg_pimd_communicate.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(COMM_WRAP)  \
                         $(CODE)/output/pimd/simpavg_pimd_communicate.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/pimd/simpavg_pimd_communicate.c
#------------------------------------------------------------------
communicate_output_pimd.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(OUTPUT_ENT) $(COMM_WRAP) $(FRND_ENT) \
                         $(CODE)/output/pimd/communicate_output_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/pimd/communicate_output_pimd.c
#------------------------------------------------------------------
communicate_output_cp_pimd.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/output/cppimd/communicate_output_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/cppimd/communicate_output_cp_pimd.c
#------------------------------------------------------------------
communicate_simpavg_cp_pimd.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) \
                         $(COMM_WRAP) $(FRND_ENT) $(MATH) $(OUTPUT_LOC) \
                         $(CODE)/output/cppimd/communicate_simpavg_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/cppimd/communicate_simpavg_cp_pimd.c
#------------------------------------------------------------------
pimd_trans_comm.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC) $(MATH)  $(COMM_WRAP) \
                         $(CODE)/energy/pimd/pimd_trans_comm.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/pimd_trans_comm.c
#------------------------------------------------------------------
build_communicate_groups.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) $(FRND_ENT) \
                         $(CODE)/communicate/build_communicate_groups.c
	$(ECHO) $@
	$(COBJ) $(CODE)/communicate/build_communicate_groups.c
#------------------------------------------------------------------
cp_communicate_coord_class.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) \
                         $(FRND_ENT)  $(MATH) $(PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/coords/cp_communicate_coord_class.c
	$(ECHO) $@
	$(COBJ) $(CODE)/interface/coords/cp_communicate_coord_class.c

#====================================================================
