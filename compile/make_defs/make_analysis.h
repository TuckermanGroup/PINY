#=================================================================
#             ANALYSIS_FILES 
#=================================================================


#=================================================================
analysis_md.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) $(COMM_WRAP) $(ANAL_LOC_ENT)\
                         $(CODE)/analysis/analysis_md.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/analysis_md.c

#-------------------------------------------------------------------
prelim_analysis.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
                         $(CODE)/analysis/prelim_analysis.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/prelim_analysis.c

#-------------------------------------------------------------------
calcul_gr.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) $(COMM_WRAP)\
                         $(CODE)/analysis/calcul_gr.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/calcul_gr.c

#-------------------------------------------------------------------
harmonic_analysis.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) $(ANAL_LOC_ENT) $(COMM_WRAP) $(MATH)\
                         $(CODE)/analysis/harmonic_analysis.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/harmonic_analysis.c

#-------------------------------------------------------------------
calcul_ickt_iso_md.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
                         $(CODE)/analysis/calcul_ickt_iso_md.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/calcul_ickt_iso_md.c

#-------------------------------------------------------------------
calcul_iikt_iso_cmd.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
                         $(CODE)/analysis/calcul_iikt_iso_cmd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/calcul_iikt_iso_cmd.c

#-------------------------------------------------------------------
calcul_msqd_atm.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
                         $(CODE)/analysis/calcul_msqd_atm.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/calcul_msqd_atm.c

#-------------------------------------------------------------------
calcul_vovt_atm.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
                         $(CODE)/analysis/calcul_vovt_atm.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/calcul_vovt_atm.c

#-------------------------------------------------------------------
analysis_pimd.o :        $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) $(COMM_WRAP)\
                         $(CODE)/analysis/analysis_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/analysis/analysis_pimd.c


#==================================================================


