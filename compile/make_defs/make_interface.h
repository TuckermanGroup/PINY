#==================================================================
#             INTERFACE_FILES
#==================================================================


#==================================================================
parse.o     :            $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_PAR) $(TYP_STAT) $(PARSE_ENT) $(PARSE_LOC) \
                         $(SIM_ENT) $(MOL_ENT) $(INTRA_ENT) $(COORD_ENT) \
                         $(COORD_LOC) $(CPEWALD_ENT) $(INTER_ENT) \
                         $(INTER_LOC) $(VPS_ENT) $(LISTS_ENT) $(SCRATCH_ENT) \
                         $(ENR_CPCON_ENT) $(REAL_LOC) $(SAMPL_CLASS_ENT) \
                         $(SAMPL_CP_ENT) $(SAMPL_CP_LOC) $(COORD_CP_ENT) \
                         $(COORD_CP_LOC) $(MATH) $(PIMD_ENT) $(PIMD_LOC) \
                         $(FRND_ENT) $(COMM_ENT) $(COMM_LOC) $(COMM_WRAP) \
                         $(INT_CPMIN_ENT) \
                         $(CODE)/interface/parse/parse.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/parse.c

#------------------------------------------------------------------
zero_bnd.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CP) \
                         $(TYP_CLASS) $(TYP_PAR) \
                         $(PARSE_LOC) \
                         $(CODE)/interface/parse/zero_bnd.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/zero_bnd.c

#------------------------------------------------------------------
zero_par.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_PAR) \
                         $(TYP_BND) $(TYP_CP) \
                         $(PARSE_LOC) \
                         $(CODE)/interface/parse/zero_par.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/zero_par.c

#------------------------------------------------------------------
zero_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(PARSE_LOC) \
                         $(CODE)/interface/parse/zero_class.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/zero_class.c

#------------------------------------------------------------------
zero_cp.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(PARSE_LOC) \
                         $(CODE)/interface/parse/zero_cp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/zero_cp.c

#==================================================================
interface_hand.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(HANDLE_ENT) \
                        $(CODE)/interface/handle/interface_hand.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/handle/interface_hand.c

#------------------------------------------------------------------
search_base_class.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(SEARCH_LOC) \
                         $(FRND_ENT) \
                         $(CODE)/interface/search_base/search_base_class.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/search_base_class.c

#------------------------------------------------------------------
data_base_handle.o    :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(SEARCH_ENT) $(INTER_ENT) $(INTER_LOC) \
                         $(INTRA_LOC)  $(HANDLE_ENT) $(FRND_ENT) \
                         $(CODE)/interface/search_base/data_base_handle.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/data_base_handle.c

#=========================================================================
control_sim_params.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_ENT) $(SIM_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(CODE)/interface/sim_params/control_sim_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/control_sim_params.c

#------------------------------------------------------------------
set_sim_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND)  \
                         $(TYP_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_LOC) $(FRND_ENT) \
                         $(CODE)/interface/sim_params/set_sim_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/set_sim_dict.c

#------------------------------------------------------------------
set_sim_params.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) $(TYP_CLASS) \
                         $(TYP_BND) $(SIM_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(TYP_STAT)\
                         $(CODE)/interface/sim_params/set_sim_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/set_sim_params.c

#=========================================================================
set_atm_NHC.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords/set_atm_NHC.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/set_atm_NHC.c

#------------------------------------------------------------------
read_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MATH) $(PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/coords/read_coord.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/read_coord.c

#------------------------------------------------------------------
molecule_decomp.o  :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords/molecule_decomp.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/molecule_decomp.c

#------------------------------------------------------------------
read_hmat.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(HANDLE_ENT) $(FRND_ENT) $(MATH) \
                         $(PIMD_LOC) \
                         $(CODE)/interface/coords/read_hmat.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/read_hmat.c

#------------------------------------------------------------------
mall_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_LOC) $(FRND_ENT) \
                         $(CODE)/interface/coords/mall_coord.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/mall_coord.c

#------------------------------------------------------------------
mall_pressure.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(COORD_ENT) $(FRND_ENT) \
                         $(CODE)/interface/coords/mall_pressure.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/mall_pressure.c

#=========================================================================
control_surf_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(SURF_PRMS_ENT) $(SURF_PRMS_LOC) \
                         $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/surf_params/control_surf_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/surf_params/control_surf_params.c

#------------------------------------------------------------------
set_surf_dict.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(SURF_PRMS_LOC) $(FRND_ENT) \
                         $(CODE)/interface/surf_params/set_surf_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/surf_params/set_surf_dict.c

#=========================================================================
control_inter_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTER_ENT) $(INTER_LOC) $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/inter_params/control_inter_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/control_inter_params.c

#------------------------------------------------------------------
set_inter_dict.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/set_inter_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/set_inter_dict.c

#------------------------------------------------------------------
get_clong.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/get_clong.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/get_clong.c

#------------------------------------------------------------------
spline_fit.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/spline_fit.c
	$(ECHO) $@
	$(COBJ) $(CODE)/interface/inter_params/spline_fit.c

#==================================================================
mall_scratch.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(SCRATCH_ENT) $(SCRATCH_LOC) $(FRND_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/interface/scratch/mall_scratch.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/scratch/mall_scratch.c

#==================================================================
control_vnhc_smpl.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/control_vnhc_smpl.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_vnhc_smpl.c

#------------------------------------------------------------------
control_vx_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COORD_LOC) \
                         $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/control_vx_smpl.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_vx_smpl.c

#------------------------------------------------------------------
control_scale_class.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                        $(CODE)/interface/vel_sampl_class/control_scale_class.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_scale_class.c

#------------------------------------------------------------------
proj_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/proj_vel_class.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/proj_vel_class.c

#------------------------------------------------------------------
samp_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CLASS) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(MATH) \
                         $(CODE)/interface/vel_sampl_class/samp_vel_class.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/samp_vel_class.c

#==================================================================
set_exclude.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_PAR) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) $(WEIGH_NODE) \
                         $(CODE)/interface/lists/set_exclude.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/set_exclude.c

#------------------------------------------------------------------
exl_sort.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_PAR) \
                         $(LISTS_LOC) $(FRND_ENT) \
                         $(CODE)/interface/lists/exl_sort.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/exl_sort.c

#------------------------------------------------------------------
mall_lists.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(LISTS_ENT) $(LISTS_LOC) $(REAL_LOC) $(FRND_ENT) \
                         $(COMM_WRAP) \
                         $(CODE)/interface/lists/mall_lists.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/mall_lists.c

#------------------------------------------------------------------
control_brnch_root_list.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(FRND_ENT) $(LISTS_ENT) $(LISTS_LOC) $(REAL_LOC) \
                         $(CODE)/interface/lists/control_brnch_root_list.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/control_brnch_root_list.c

#------------------------------------------------------------------
block_intra_lists.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) \
                         $(CODE)/interface/lists/block_intra_lists.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/block_intra_lists.c

#------------------------------------------------------------------
class_par_forc_lists.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(FRND_ENT) $(LISTS_ENT) $(LISTS_LOC) \
                         $(CODE)/interface/lists/class_par_forc_lists.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/class_par_forc_lists.c

#==================================================================





