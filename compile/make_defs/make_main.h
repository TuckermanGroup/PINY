#=================================================================
#              MAIN_FILES 
#=================================================================


#=================================================================
main.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(MAIN_ENT) $(MAIN_LOC) $(MAIN_CP_LOC) $(PARSE_ENT) \
                         $(COMM_WRAP) $(COMM_ENT) \
                         $(CODE)/main/main.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/main.c

#------------------------------------------------------------------
control_md.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(TYP_STAT) \
                         $(MAIN_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(OUTPUT_ENT) $(ANAL_MD_ENT) $(SMPL_CLASS_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(CODE)/main/md/control_md.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/md/control_md.c

#------------------------------------------------------------------
control_min.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_STAT) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_MIN_ENT) \
                         $(INT_MIN_LOC) $(OUTPUT_ENT) $(OUTPUT_LOC) $(ANAL_MD_ENT) \
                         $(SMPL_CLASS_ENT) $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(CODE)/main/md/control_min.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/md/control_min.c


#------------------------------------------------------------------
control_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(TYP_STAT) \
                         $(MAIN_LOC) $(MATH) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_ENT) $(OUTPUT_ENT) $(ANAL_MD_ENT) \
                         $(SMPL_CLASS_ENT) $(ENR_CTRL_ENT) $(PIMD_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/main/pimd/control_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/pimd/control_pimd.c

#------------------------------------------------------------------
control_debug.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(TYP_STAT)\
                         $(MAIN_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(REAL_LOC) $(COMM_WRAP) \
                         $(CODE)/main/debug/control_debug.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/debug/control_debug.c

#------------------------------------------------------------------
control_debug_pimd.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_STAT)\
                         $(MAIN_LOC) $(MATH) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_ENT) $(OUTPUT_ENT) $(OUTPUT_LOC) \
                         $(REAL_LOC) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(CODE)/main/debug/control_debug_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/debug/control_debug_pimd.c

#------------------------------------------------------------------
auto_exit.o     :        $(STANDARD) $(DEFINES) \
                         $(MAIN_LOC) \
                         $(CODE)/main/auto_exit.c
	$(ECHO) $@
	$(COBJ) $(CODE)/main/auto_exit.c

#=================================================================






