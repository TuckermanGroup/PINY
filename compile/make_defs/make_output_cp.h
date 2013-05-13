#====================================================================
#             OUTPUT_CP_FILES
#====================================================================


#=================================================================
output_cp_min.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(MATH) $(OUTPUT_LOC) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/output/cp/output_cp_min.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/cp/output_cp_min.c
#------------------------------------------------------------------
output_cp.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(FRND_ENT) $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \
                         $(CODE)/output/cp/output_cp.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/cp/output_cp.c
#------------------------------------------------------------------
simpavg_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_ENT) $(OUTPUT_CP_LOC)  \
                         $(OUTPUT_LOC) $(COMM_WRAP) \
                         $(CODE)/output/cp/simpavg_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/cp/simpavg_cp.c

#------------------------------------------------------------------
simpavg_cp_communicate.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_CP) \
                         $(OUTPUT_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/output/cp/simpavg_cp_communicate.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/cp/simpavg_cp_communicate.c

#------------------------------------------------------------------
output_cp_pimd.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(FRND_ENT) $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \
                         $(CODE)/output/cppimd/output_cp_pimd.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/cppimd/output_cp_pimd.c

#------------------------------------------------------------------
simpavg_cp_pimd.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(FRND_ENT) \
                         $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \
                         $(CODE)/output/cppimd/simpavg_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/cppimd/simpavg_cp_pimd.c

#====================================================================

