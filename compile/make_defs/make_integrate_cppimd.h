#==================================================================
#              INTEGRATE_CPPIMD_FILES
#==================================================================


#==================================================================
int_NVT_cp_pimd.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(INT_MD_LOC) $(INT_PIMD_LOC) $(ENR_CPCON_ENT) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) $(PIMD_ENT) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
                         $(CODE)/integrate/cppimd/int_NVT_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cppimd/int_NVT_cp_pimd.c

#------------------------------------------------------------------
int_NPTI_cp_pimd.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(INT_MD_LOC) $(INT_PIMD_LOC) $(ENR_CPCON_ENT) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) $(PIMD_ENT) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
                         $(CODE)/integrate/cppimd/int_NPTI_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cppimd/int_NPTI_cp_pimd.c

#------------------------------------------------------------------
int_NPTF_cp_pimd.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(COMM_WRAP) $(INT_MD_LOC) $(INT_PIMD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) $(MATH) $(FRND_ENT) \
                         $(COORD_LOC) \
                         $(CODE)/integrate/cppimd/int_NPTF_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cppimd/int_NPTF_cp_pimd.c

#------------------------------------------------------------------
int_utils_cp_pimd.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) \
                         $(CODE)/integrate/cppimd/int_utils_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cppimd/int_utils_cp_pimd.c

#==================================================================








