#=================================================================
#             MAIN_CP_FILES
#=================================================================


#=================================================================
control_debug_cp.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(TYP_STAT)\
                         $(MAIN_CP_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(REAL_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_CP_ENT) $(COMM_WRAP) \
                         $(CODE)/main/debug/control_debug_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/debug/control_debug_cp.c

#------------------------------------------------------------------
control_debug_cp_pimd.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_STAT)\
                         $(MAIN_CP_LOC) $(MATH) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(INTRA_CON_ENT) $(REAL_LOC) \
                         $(INT_PIMD_LOC) $(ENR_CTRL_CP_ENT) $(INT_MD_ENT) \
                         $(INT_PIMD_ENT) $(INT_CP_LOC) $(INT_CPPIMD_LOC) \
                         $(COMM_WRAP) \
                         $(CODE)/main/debug/control_debug_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/debug/control_debug_cp_pimd.c

#------------------------------------------------------------------
control_cp_min.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_STAT) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CPMIN_ENT) \
                         $(ENR_CPCON_ENT) $(ENR_CPCON_LOC) $(ENR_CTRL_CP_ENT) \
                         $(OUTPUT_ENT) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(ANAL_MD_ENT) $(ANAL_CP_ENT) \
                         $(COMM_WRAP) $(MATH) $(FRND_ENT) \
                         $(CODE)/main/cp/control_cp_min.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/cp/control_cp_min.c

#------------------------------------------------------------------
control_cp_pimd_min.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_STAT) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CPMIN_ENT) \
                         $(ENR_CPCON_ENT) $(OUTPUT_ENT) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(ANAL_MD_ENT) $(ANAL_CP_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/main/cp/control_cp_pimd_min.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/cp/control_cp_pimd_min.c

#------------------------------------------------------------------
control_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_STAT) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CP_ENT) $(INT_CP_LOC) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_ENT) $(SMPL_CP_ENT) $(SMPL_CLASS_ENT) \
                         $(ANAL_CP_ENT) $(OUTPUT_CP_ENT) $(INT_MD_ENT) \
                         $(INT_MD_LOC) $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CPCON_ENT) $(COMM_WRAP) \
                         $(CODE)/main/cp/control_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/cp/control_cp.c

#------------------------------------------------------------------
control_cp_pimd.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(TYP_STAT) \
                         $(TYP_BND) $(TYP_PAR) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CPPIMD_ENT) \
                         $(INT_CPPIMD_LOC) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_ENT) $(ENR_CTRL_CP_ENT) $(SMPL_CP_ENT) \
                         $(SMPL_CLASS_ENT) $(ANAL_CP_ENT) $(INTRA_CON_ENT) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(INT_MD_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(PIMD_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/main/cp/control_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/cp/control_cp_pimd.c

#=================================================================




