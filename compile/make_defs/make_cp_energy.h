#=====================================================================
#            CP_ENERGY_FILES 
#=====================================================================


#=====================================================================
cp_energy_control.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_ENT) $(ENR_CTRL_CP_LOC) \
                         $(INTRA_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(ENR_CTRL_CP_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/control_cp/cp_energy_control.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control_cp/cp_energy_control.c

#------------------------------------------------------------------
energy_control_elec.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INTRA_ENT) $(ENR_CP_ENT) $(ENR_CP_LOC) \
                         $(COMM_WRAP) \
                         $(CODE)/energy/control_cp/energy_control_elec.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control_cp/energy_control_elec.c

#------------------------------------------------------------------
test_energy_cp.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(COMM_WRAP) $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_ENT) $(ENR_CTRL_CP_LOC) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) \
                         $(CODE)/energy/control_cp/test_energy_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control_cp/test_energy_cp.c

#------------------------------------------------------------------
cp_ks_energy.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) $(REAL_ENT) \
                         $(CODE)/energy/cp/cp_ks_energy.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_ks_energy.c

#------------------------------------------------------------------
cp_energy_ee_rho.o :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp/cp_energy_ee_rho.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_ee_rho.c

#------------------------------------------------------------------
cp_energy_ee_rho_ke.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp/cp_energy_ee_rho_ke.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_ee_rho_ke.c

#------------------------------------------------------------------
cp_coef_force_tau_fun.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp/cp_coef_force_tau_fun.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_coef_force_tau_fun.c

#------------------------------------------------------------------
cp_energy_ee_grad_rho.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(FRND_ENT) \
                         $(CODE)/energy/cp/cp_energy_ee_grad_rho.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_ee_grad_rho.c

#------------------------------------------------------------------
cp_energy_eext.o    :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(ENR_CPCON_LOC) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp/cp_energy_eext.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_eext.c

#------------------------------------------------------------------
cp_energy_eext_nonloc.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(FRND_ENT) \
                         $(ENR_CPCON_LOC) $(MATH) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/cp/cp_energy_eext_nonloc.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_eext_nonloc.c

#------------------------------------------------------------------
cp_energy_eext_nonloc_gh.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(FRND_ENT) \
                         $(ENR_CPCON_LOC) $(MATH) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/cp/cp_energy_eext_nonloc_gh.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/cp_energy_eext_nonloc_gh.c

#------------------------------------------------------------------
xc_functionals.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) \
                         $(CODE)/energy/cp/xc_functionals.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/xc_functionals.c

#------------------------------------------------------------------
constraint_control_cp.o  : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_ENT) $(ENR_CPCON_LOC) \
                         $(CODE)/energy/cp_con/constraint_control_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/constraint_control_cp.c

#------------------------------------------------------------------
cp_con.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp_con/cp_con.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/cp_con.c

#------------------------------------------------------------------
cp_orth_rot_utils.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
                         $(CODE)/energy/cp_con/cp_orth_rot_utils.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/cp_orth_rot_utils.c

#------------------------------------------------------------------
cp_con_utils.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp_con/cp_con_utils.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/cp_con_utils.c

#------------------------------------------------------------------
orth_rot_control_cp.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/cp_con/orth_rot_control_cp.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/orth_rot_control_cp.c

#------------------------------------------------------------------
cp_transpose.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
                         $(CODE)/energy/cp_con/cp_transpose.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp_con/cp_transpose.c

#------------------------------------------------------------------
control_spread_rho.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(MATH) $(FRND_ENT) $(ENR_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/cp/control_spread_rho.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/control_spread_rho.c

#------------------------------------------------------------------
control_contract_rho.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(MATH) $(FRND_ENT) $(ENR_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/cp/control_contract_rho.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/cp/control_contract_rho.c

#===================================================================











