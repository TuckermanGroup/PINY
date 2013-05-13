#======================================================================
#             ENERGY_PIMD_FILES 
#======================================================================


#======================================================================
energy_control_pimd.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/control/energy_control_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_pimd.c

#------------------------------------------------------------------
test_energy_pimd.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/control/test_energy_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/test_energy_pimd.c

#------------------------------------------------------------------
cp_energy_control_pimd.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_CP_LOC) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_ENT) $(ENR_CP_ENT) \
                         $(ENR_CP_LOC) $(REAL_ENT) $(REC_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_LOC) $(INTRA_CON_ENT) $(MATH) \
                         $(PIMD_ENT) $(PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/control_cp/cp_energy_control_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control_cp/cp_energy_control_pimd.c

#------------------------------------------------------------------
test_energy_cp_pimd.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CTRL_CP_LOC) $(ENR_CP_ENT) $(ENR_CP_LOC) \
                         $(INTRA_ENT) $(INTRA_CON_ENT) $(REAL_ENT) \
                         $(REC_ENT) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/control_cp/test_energy_cp_pimd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control_cp/test_energy_cp_pimd.c

#------------------------------------------------------------------
transform_cnt.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) \
                         $(CODE)/energy/pimd/transform_cnt.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/transform_cnt.c

#------------------------------------------------------------------
transform_stg.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) \
                         $(CODE)/energy/pimd/transform_stg.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/transform_stg.c

#------------------------------------------------------------------
transform_cnt_par.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/pimd/transform_cnt_par.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/transform_cnt_par.c

#------------------------------------------------------------------
transform_stg_par.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) \
                         $(CODE)/energy/pimd/transform_stg_par.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/transform_stg_par.c

#------------------------------------------------------------------
control_pimd_trans.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_ENT)  $(PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/pimd/control_pimd_trans.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/control_pimd_trans.c

#------------------------------------------------------------------
pimd_estimators.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(PIMD_LOC) $(ENR_CTRL_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/pimd/pimd_estimators.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/pimd/pimd_estimators.c

#======================================================================






