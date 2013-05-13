#====================================================================
#                OUTPUT_FILES
#====================================================================


#=================================================================
output_md.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
                         $(CODE)/output/md/output_md.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/md/output_md.c

#------------------------------------------------------------------
output_min.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
                         $(CODE)/output/md/output_min.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/md/output_min.c

#------------------------------------------------------------------
output_pimd.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(FRND_ENT) $(MATH) \
                         $(CODE)/output/pimd/output_pimd.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/pimd/output_pimd.c

#------------------------------------------------------------------
simpavg_md.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) \
                         $(CODE)/output/md/simpavg_md.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/md/simpavg_md.c

#------------------------------------------------------------------
simpavg_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(COMM_WRAP) \
                         $(CODE)/output/pimd/simpavg_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/pimd/simpavg_pimd.c

#------------------------------------------------------------------
simpavg_md_communicate.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(COMM_WRAP) \
                         $(CODE)/output/md/simpavg_md_communicate.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/md/simpavg_md_communicate.c

#------------------------------------------------------------------
get_cell.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_BND) $(TYP_GEN) $(TYP_CLASS) \
                         $(OUTPUT_LOC) \
                         $(CODE)/output/md/get_cell.c
	$(ECHO) $@
	$(COBJ) $(CODE)/output/md/get_cell.c

#------------------------------------------------------------------
write_gen_header.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(FRND_ENT)  $(OUTPUT_LOC) $(OUTPUT_CP_LOC) \
                         $(CODE)/output/write_gen_header.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/output/write_gen_header.c

#====================================================================
