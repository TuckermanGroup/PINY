#==================================================================
#              INTEGRATE_PIMD_FILES
#==================================================================


#==================================================================
int_NVT_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(COMM_WRAP) \
                         $(INT_MD_LOC) $(INTRA_CON_ENT) $(ENR_CTRL_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) \
                         $(CODE)/integrate/pimd/int_NVT_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/pimd/int_NVT_pimd.c

#------------------------------------------------------------------
int_NPTI_pimd.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(COMM_WRAP) \
                         $(INT_MD_LOC)  $(INTRA_CON_ENT)  $(ENR_CTRL_ENT) \
                         $(PIMD_ENT)  $(PIMD_LOC)  $(MATH) \
                         $(CODE)/integrate/pimd/int_NPTI_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/pimd/int_NPTI_pimd.c

#------------------------------------------------------------------
int_NPTF_pimd.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_STAT) \
                         $(COMM_WRAP) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_LOC) $(INT_MD_ENT)  $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(PIMD_ENT)  $(PIMD_LOC) $(MATH) \
                         $(CODE)/integrate/pimd/int_NPTF_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/pimd/int_NPTF_pimd.c

#------------------------------------------------------------------
int_utils_pimd.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(COMM_WRAP) $(INT_PIMD_LOC) \
                         $(CODE)/integrate/pimd/int_utils_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/pimd/int_utils_pimd.c

#-------------------------------------------------------------------------
int_0_to_dt2_pimd.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_STAT)\
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(INT_MD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(PIMD_ENT) \
                         $(MATH)  $(COMM_WRAP) \
                         $(CODE)/integrate/int_0_to_dt2_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/int_0_to_dt2_pimd.c

#------------------------------------------------------------------
int_dt2_to_dt_pimd.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_PIMD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(CODE)/integrate/int_dt2_to_dt_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/int_dt2_to_dt_pimd.c

#==================================================================

