#==================================================================
#               INTEGRATE_CP_FILES
#==================================================================


#==================================================================
int_NVE_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
                         $(CODE)/integrate/cp/int_NVE_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_NVE_cp.c

#------------------------------------------------------------------
int_NVT_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_ENT) \
                         $(INT_MD_LOC)  $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
                         $(CODE)/integrate/cp/int_NVT_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_NVT_cp.c

#------------------------------------------------------------------
int_NPTI_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_CP_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) $(COORD_LOC) \
                         $(CODE)/integrate/cp/int_NPTI_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_NPTI_cp.c

#------------------------------------------------------------------
int_NPTF_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_CP_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) $(COORD_LOC) \
                         $(CODE)/integrate/cp/int_NPTF_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_NPTF_cp.c

#------------------------------------------------------------------
int_utils_cp.o    :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/integrate/cp/int_utils_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_utils_cp.c

#------------------------------------------------------------------
min_STD_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(ENR_CP_ENT) \
                         $(ENR_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/integrate/cpmin/min_STD_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cpmin/min_STD_cp.c

#------------------------------------------------------------------
min_CG_cp.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CP_LOC) $(FRND_ENT) $(COMM_WRAP) $(MATH) \
                         $(CODE)/integrate/cpmin/min_CG_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cpmin/min_CG_cp.c
#------------------------------------------------------------------
move_atm.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(MATH) \
                         $(CODE)/integrate/cpmin/move_atm.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cpmin/move_atm.c

#------------------------------------------------------------------
min_DIIS_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) \
                         $(ENR_CPCON_ENT) $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(MATH) $(FRND_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/integrate/cpmin/min_DIIS_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cpmin/min_DIIS_cp.c

#------------------------------------------------------------------
shuffle_states.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) $(ENR_CPCON_LOC) \
                         $(CODE)/integrate/cpmin/shuffle_states.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cpmin/shuffle_states.c

#------------------------------------------------------------------
cp_gauss.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) \
                         $(CODE)/integrate/cp/cp_gauss.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/cp_gauss.c

#------------------------------------------------------------------
int_0_to_dt2_cp.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) $(COORD_LOC) \
                         $(CODE)/integrate/cp/int_0_to_dt2_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_0_to_dt2_cp.c

#------------------------------------------------------------------
int_dt2_to_dt_cp.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) $(INT_CPPIMD_LOC) $(INT_PIMD_LOC) \
                         $(CODE)/integrate/cp/int_dt2_to_dt_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/cp/int_dt2_to_dt_cp.c

#==================================================================
