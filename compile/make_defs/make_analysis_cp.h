#=================================================================
#             ANALYSIS_CP_FILES 
#==================================================================

#=================================================================
analysis_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ANAL_CP_ENT) $(ANAL_LOC_ENT) $(ANAL_MD_ENT) $(COMM_WRAP) $(MATH)\
                         $(CODE)/analysis/analysis_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/analysis_cp.c

#-------------------------------------------------------------------
analysis_cp_pimd.o     : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ANAL_CP_ENT) $(ANAL_LOC_ENT) \
                         $(CODE)/analysis/analysis_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/analysis_cp_pimd.c

#-------------------------------------------------------------------
