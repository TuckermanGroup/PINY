#==================================================================
#               INTEGRATE_MD_FILES
#==================================================================


#==================================================================
int_NVE.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT)  $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
                         $(CODE)/integrate/md/int_NVE.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NVE.c

#------------------------------------------------------------------
min_STD.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MIN_ENT)  $(INT_MIN_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/integrate/min/min_STD.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/min/min_STD.c

#------------------------------------------------------------------
min_CG.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MIN_ENT) $(INT_MIN_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(FRND_ENT) \
                         $(CODE)/integrate/min/min_CG.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/min/min_CG.c



#------------------------------------------------------------------
int_NVE_res.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
                         $(CODE)/integrate/md/int_NVE_res.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NVE_res.c

#------------------------------------------------------------------
int_NVT.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(CODE)/integrate/md/int_NVT.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NVT.c

#------------------------------------------------------------------
int_NVT_res.o    :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
                         $(CODE)/integrate/md/int_NVT_res.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NVT_res.c

#------------------------------------------------------------------
int_NPTI.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) $(COMM_WRAP) \
                         $(CODE)/integrate/md/int_NPTI.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NPTI.c

#------------------------------------------------------------------
int_NPTI_res.o   :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) \
                         $(CODE)/integrate/md/int_NPTI_res.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NPTI_res.c

#------------------------------------------------------------------
int_NPTF.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH)  $(COMM_WRAP) \
                         $(CODE)/integrate/md/int_NPTF.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NPTF.c

#------------------------------------------------------------------
int_NPTF_res.o   :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) \
                         $(CODE)/integrate/md/int_NPTF_res.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_NPTF_res.c

#------------------------------------------------------------------
int_utils.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_LOC) \
                         $(CODE)/integrate/md/int_utils.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/md/int_utils.c

#------------------------------------------------------------------
int_0_to_dt2.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(INT_PIMD_LOC) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH)  $(COMM_WRAP) \
                         $(CODE)/integrate/int_0_to_dt2.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/int_0_to_dt2.c

#------------------------------------------------------------------
int_dt2_to_dt.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_STAT)\
                         $(INT_PIMD_ENT) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(MATH) \
                         $(CODE)/integrate/int_dt2_to_dt.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/int_dt2_to_dt.c

#------------------------------------------------------------------
move_pos_box.o : $(STANDARD) $(DEFINES)  $(TYP_CLASS) $(INT_PIMD_LOC)\
                   $(TYP_GEN) $(TYP_BND) $(TYP_STAT)\
                   $(ENR_PIMD_ENT)\
                   $(INTRA_CON_LOC) $(ENR_CTRL_ENT)\
                   $(CODE)/integrate/move_pos_box.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/move_pos_box.c
#------------------------------------------------------------------
move_vel_vbox.o : $(STANDARD) $(DEFINES)  $(TYP_CLASS) $(INT_PIMD_LOC)\
                   $(TYP_GEN) $(TYP_BND) $(ENR_PIMD_ENT) $(TYP_STAT)\
                   $(INTRA_CON_LOC) $(ENR_CTRL_ENT)\
                   $(CODE)/integrate/move_vel_vbox.c
	$(ECHO) $@
	$(COBJ) $(CODE)/integrate/move_vel_vbox.c

#==================================================================
