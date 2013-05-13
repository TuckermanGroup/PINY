#==================================================================
#           INTERFACE_CP_FILES
#==================================================================


#==================================================================
MOL_PARMS1 = $(CODE)/interface/mol_params
DMOL_PARMS1 = $(DCODE)/interface/mol_params
#==================================================================


#==================================================================
set_wave_params.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CLASS) $(TYP_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS1)/set_params/set_wave_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS1)/set_params/set_wave_params.c

#------------------------------------------------------------------
set_coef_NHC.o    :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords_cp/set_coef_NHC.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/set_coef_NHC.c

#------------------------------------------------------------------
read_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_ENT) $(COORD_CP_ENT) $(HANDLE_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords_cp/read_coef.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/read_coef.c

#------------------------------------------------------------------
mall_properties.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) \
                         $(FRND_ENT) \
                         $(CODE)/interface/coords_cp/mall_properties.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/mall_properties.c

#------------------------------------------------------------------
gen_wave.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) $(TYP_CLASS) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP) $(ENR_CPCON_LOC) \
                         $(ENR_CP_LOC) $(CPEWALD_LOC) $(COORD_CP_LOC) \
                         $(HANDLE_ENT) \
                         $(CODE)/interface/coords_cp/gen_wave.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/gen_wave.c

#------------------------------------------------------------------
mall_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) $(FRND_ENT) \
                         $(CODE)/interface/coords_cp/mall_coef.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/mall_coef.c

#------------------------------------------------------------------
control_set_cp_ewald.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(TYP_BND) $(TYP_PAR) \
                         $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/cp_ewald/control_set_cp_ewald.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/cp_ewald/control_set_cp_ewald.c

#------------------------------------------------------------------
set_cp_ewald.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(CPEWALD_LOC) $(MATH) $(FRND_ENT) \
                         $(CODE)/interface/cp_ewald/set_cp_ewald.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/cp_ewald/set_cp_ewald.c

#------------------------------------------------------------------
search_base_cp.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(SEARCH_ENT) $(VPS_LOC) $(HANDLE_ENT) \
                         $(CODE)/interface/search_base/search_base_cp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/search_base_cp.c

#------------------------------------------------------------------
proj_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(SMPL_CP_LOC) $(SMPL_CP_ENT) $(ENR_CPCON_ENT) \
                         $(CODE)/interface/vel_sampl_cp/proj_vel_cp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/proj_vel_cp.c

#------------------------------------------------------------------
set_vps_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(CODE)/interface/vps_params/set_vps_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/set_vps_dict.c

#------------------------------------------------------------------
samp_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(SMPL_CP_LOC)  $(MATH) \
                         $(CODE)/interface/vel_sampl_cp/samp_vel_cp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/samp_vel_cp.c

#------------------------------------------------------------------
control_vps_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) \
                         $(TYP_CLASS) $(TYP_BND) \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(COMM_WRAP) $(MATH) \
                         $(CODE)/interface/vps_params/control_vps_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/control_vps_params.c

#------------------------------------------------------------------
weight_node_gauss_hermite.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(CODE)/interface/vps_params/weight_node_gauss_hermite.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/weight_node_gauss_hermite.c

#------------------------------------------------------------------
control_vc_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_vc_smpl.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_vc_smpl.c

#------------------------------------------------------------------
control_vcnhc_smpl.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_vcnhc_smpl.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_vcnhc_smpl.c

#------------------------------------------------------------------
control_scale_cp.o    :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP)  \
                         $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_scale_cp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_scale_cp.c

#==================================================================



