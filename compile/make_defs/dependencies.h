
 
=============================================================
MAKE_ANALYSIS
=============================================================
 
---------------------------------------
analysis_md.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
---------------------------------------
anal_pimd.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ANAL_MD_ENT) \
 
=============================================================
MAKE_ANAL_CP
=============================================================
 
---------------------------------------
analysis_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ANAL_CP_ENT) \
---------------------------------------
analysis_cp_pimd.o     : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ANAL_CP_ENT) \
 
=============================================================
MAKE_COMMUNICATE
=============================================================
 
---------------------------------------
communicate_wrappers.o : $(STANDARD) $(DEFINES) \
                         $(COMM_WRAP)  \
---------------------------------------
control_group_communicators.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_PAR) $(COMM_ENT) $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
com_interface.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_ENT) $(COMM_WRAP) \
                         $(COMM_LOC) \
---------------------------------------
comm_class_info.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CP) $(TYP_CLASS) \
                         $(TYP_BND) $(COMM_WRAP) $(COMM_LOC) $(MATH) \
---------------------------------------
comm_bond_info.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
comm_bond_data.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
comm_cp_info.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CLASS) $(TYP_CP) \
                         $(TYP_PAR) $(COMM_WRAP) $(COMM_LOC) $(MATH) \
---------------------------------------
comm_class_data.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
comm_class_list.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
comm_cp_data.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) $(FRND_ENT) \
---------------------------------------
communicate_options.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) $(TYP_CLASS)  \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) \
---------------------------------------
mall_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_CP)    $(COMM_LOC)  $(FRND_ENT) $(PARSE_ENT) \
---------------------------------------
mall_bond.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_CP) $(COMM_LOC) $(FRND_ENT)  $(PARSE_ENT) \
---------------------------------------
communicate_coord.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) \
                         $(FRND_ENT)  $(MATH) $(PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
communicate_test_energy_pimd.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(TYP_CP)  $(ENR_CTRL_LOC) $(COMM_WRAP) \
---------------------------------------
communicate_utils_pimd.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INT_PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
simpavg_pimd_communicate.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(COMM_WRAP)  \
---------------------------------------
communicate_output_pimd.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(OUTPUT_ENT) $(COMM_WRAP) $(FRND_ENT) \
---------------------------------------
communicate_output_cp_pimd.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
communicate_simpavg_cp_pimd.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC)
                         $(COMM_WRAP) $(FRND_ENT) $(MATH) $(OUTPUT_LOC) \
---------------------------------------
pimd_trans_comm.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(PIMD_LOC) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
build_communicate_groups.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(TYP_CP)  $(COMM_WRAP) $(COMM_LOC) $(FRND_ENT) \
---------------------------------------
cp_communicate_coord_class.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) \
                         $(FRND_ENT)  $(MATH) $(PIMD_LOC) $(COMM_WRAP) \
 
=============================================================
MAKE_CP_ENERGY
=============================================================
 
---------------------------------------
cp_energy_control.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_ENT) $(ENR_CTRL_CP_LOC) \
                         $(INTRA_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(ENR_CTRL_CP_ENT) $(COMM_WRAP) \
---------------------------------------
energy_control_elec.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INTRA_ENT) $(ENR_CP_ENT) $(ENR_CP_LOC) \
                         $(COMM_WRAP) \
---------------------------------------
test_energy_cp.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(COMM_WRAP) $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_ENT) $(ENR_CTRL_CP_LOC) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) \
---------------------------------------
cp_ks_energy.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(ENR_CP_LOC) $(ENR_CPCON_LOC) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) \
---------------------------------------
cp_energy_ee_rho.o :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(COMM_WRAP) \
---------------------------------------
cp_energy_ee_grad_rho.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) $(FRND_ENT) \
---------------------------------------
cp_energy_eext.o    :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(ENR_CPCON_LOC) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
cp_energy_eext_nonloc.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(FRND_ENT) \
                         $(ENR_CPCON_LOC) $(MATH) \
---------------------------------------
xc_functionals.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CP_LOC) $(MATH) \
---------------------------------------
constraint_control_cp.o  : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_ENT) $(ENR_CPCON_LOC) \
---------------------------------------
cp_con.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
cp_orth_rot_utils.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
CPCON_utils.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
orth_rot_control_cp.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
cp_transpose.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
 
=============================================================
MAKE_ENERGY_INTER
=============================================================
 
---------------------------------------
energy_control.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(MATH) \
---------------------------------------
energy_control_initial.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(COMM_WRAP)    $(MATH) \
---------------------------------------
energy_control_intra.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)     $(MATH) \
---------------------------------------
energy_control_inter_real.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)  $(MATH) \
---------------------------------------
energy_control_inter_recip.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)  $(MATH) \
---------------------------------------
energy_control_final.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT)  \
                         $(REC_ENT) $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) $(MATH) \
---------------------------------------
test_energy.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
force_control.o   :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_ENT)  $(REAL_LOC) \
---------------------------------------
force_nolst.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
---------------------------------------
force_verlst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
---------------------------------------
force_lnklst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
---------------------------------------
force_npol.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_LOC) $(ENR_CTRL_ENT) $(MATH) \
---------------------------------------
period.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) \
---------------------------------------
nbr_list_control.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_ENT)  $(REAL_LOC) $(PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
verlist_control.o    :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
verlist_create.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_LOC)  $(ENR_CTRL_ENT) \
---------------------------------------
make_lnk_lst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) \
---------------------------------------
make_lnk_map.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) \
---------------------------------------
lnk_lst_dis.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(MATH) \
---------------------------------------
ewald3d.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(MATH)    $(REC_ENT) \
---------------------------------------
ewald3d_both.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(REC_ENT) $(MATH) \
---------------------------------------
ewald3d_self_bgr.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(REC_ENT) $(MATH) \
---------------------------------------
ewald3d_pme.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) \
                         $(MATH) $(REC_ENT) $(COMM_WRAP) \
 
=============================================================
MAKE_ENERGY_INTRA
=============================================================
 
---------------------------------------
bond.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
bond_both.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
bond_watts.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
bend.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
bend_bnd.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
tors.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
onefour.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
ecorr.o     :            $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
ghost_control.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_CON_ENT) \
---------------------------------------
constraint_control.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_ENT) $(INTRA_CON_LOC) $(COMM_WRAP) \
                         $(MATH) \
---------------------------------------
bond_con.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_CON_LOC) $(ENR_CTRL_ENT) \
---------------------------------------
bond_con_rolli.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_CON_LOC) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
bond_con_rollf.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_CON_LOC) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
group_bond_con_21.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
group_bond_con_23.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
group_bond_con_33.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
group_bond_con_43.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
group_bond_con_46.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
group_bond_con_rolli_21.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
group_bond_con_rolli_23.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rolli_33.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
group_bond_con_rolli_43.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \
                         $(COMM_WRAP) \
---------------------------------------
group_bond_con_rolli_46.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rollf_21.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rollf_23.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rollf_33.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rollf_43.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
group_bond_con_rollf_46.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
group_bond_con_util.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(INTRA_CON_LOC) $(FRND_ENT) $(MATH) \

=============================================================
MAKE_ENERGY_PIMD
=============================================================

---------------------------------------
energy_control_pimd.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) $(MATH) $(COMM_WRAP) \
---------------------------------------
test_ENR_pimd.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
cp_energy_control_pimd.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_CP_LOC) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_ENT) $(ENR_CP_ENT) \
                         $(ENR_CP_LOC) $(REAL_ENT) $(REC_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_LOC) $(INTRA_CON_ENT) $(MATH) \
                         $(PIMD_ENT) $(PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
test_energy_cp_pimd.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CTRL_CP_LOC) $(ENR_CP_ENT) $(ENR_CP_LOC) \
                         $(INTRA_ENT) $(INTRA_CON_ENT) $(REAL_ENT) \
                         $(REC_ENT) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
transform_cnt.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) \
---------------------------------------
transform_stg.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) \
---------------------------------------
transform_cnt_par.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC)  $(MATH) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
transform_stg_par.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) \
---------------------------------------
control_pimd_trans.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(PIMD_ENT)  $(PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
pimd_estimators.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(PIMD_LOC) $(ENR_CTRL_ENT) $(FRND_ENT) $(COMM_WRAP) \
=============================================================
MAKE_FRIEND
============================================================

---------------------------------------
cmalloc.o     :          $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
---------------------------------------
friend_lib.o   :         $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
---------------------------------------
mal_verify.o     :       $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \

=============================================================
MAKE_INT_CP
=============================================================

---------------------------------------
int_nve_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
---------------------------------------
int_nvt_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_ENT) \
                         $(INT_MD_LOC)  $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(COORD_LOC)
---------------------------------------
int_npti_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_CP_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) $(COORD_LOC) \
---------------------------------------
int_nptf_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_CP_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(MATH) \
                         $(COMM_WRAP) $(COORD_LOC) \
---------------------------------------
int_utils_cp.o    :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) $(COMM_WRAP) \
---------------------------------------
min_std_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(ENR_CP_ENT) \
                         $(ENR_CP_LOC) $(COMM_WRAP) \
---------------------------------------
min_cg_cp.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CP_LOC) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
move_atm.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(ENR_CPCON_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(MATH) \
---------------------------------------
min_diis_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) \
                         $(ENR_CPCON_ENT) $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CP_ENT) $(ENR_CP_LOC) $(MATH) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
shuffle_states.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(INT_CPMIN_ENT) $(INT_CPMIN_LOC) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) $(ENR_CPCON_LOC) \
---------------------------------------
cp_gauss.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) \
---------------------------------------
int_0_to_dt2_cp.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) $(COORD_LOC) \
---------------------------------------
int_dt2_to_dt_cp.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(INT_MD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(INTRA_CON_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \

=============================================================
MAKE_INT_CPPIMD
=============================================================

---------------------------------------
int_nvt_cp_pimd.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(INT_MD_LOC) $(INT_PIMD_LOC) $(ENR_CPCON_ENT) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) $(PIMD_ENT) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
---------------------------------------
int_npti_cp_pimd.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(INT_MD_LOC) $(INT_PIMD_LOC) $(ENR_CPCON_ENT) \
                         $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) $(PIMD_ENT) \
                         $(PIMD_LOC) $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(COORD_LOC) \
---------------------------------------
int_nptf_cp_pimd.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CPPIMD_ENT) $(INT_CPPIMD_LOC) $(INT_CP_LOC) \
                         $(COMM_WRAP) $(INT_MD_LOC) $(INT_PIMD_LOC) \
                         $(ENR_CPCON_ENT) $(ENR_CTRL_CP_ENT) $(INTRA_CON_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) $(MATH) $(FRND_ENT) \
                         $(COORD_LOC) \
---------------------------------------
int_utils_cp_pimd.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(INT_CP_LOC) \

=============================================================
MAKE_INT_MD
=============================================================

---------------------------------------
int_nve.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT)  $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
---------------------------------------
min_std.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MIN_ENT)  $(INT_MIN_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
min_cg.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MIN_ENT) $(INT_MIN_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(FRND_ENT) \
---------------------------------------
int_nve_res.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
---------------------------------------
int_nvt.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(COMM_WRAP) \
---------------------------------------
int_nvt_res.o    :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
---------------------------------------
int_nvt_isok.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(COMM_WRAP) \
---------------------------------------
int_nvt_isok_res.o    :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) \
---------------------------------------
int_npti.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) $(COMM_WRAP) \
---------------------------------------
int_npti_res.o   :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) \
---------------------------------------
int_nptf.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH)  $(COMM_WRAP) \
---------------------------------------
int_nptf_res.o   :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH) \
---------------------------------------
int_utils.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_LOC) \
---------------------------------------
int_0_to_dt2.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(MATH)  $(COMM_WRAP) \
---------------------------------------
int_dt2_to_dt.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_PIMD_ENT) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
                         $(MATH) \

=============================================================
MAKE_INT_PIMD
=============================================================

---------------------------------------
int_nvt_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(COMM_WRAP) \
                         $(INT_MD_LOC) $(INTRA_CON_ENT) $(ENR_CTRL_ENT) \
                         $(PIMD_ENT) $(PIMD_LOC) \
---------------------------------------
int_npti_pimd.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(COMM_WRAP) \
                         $(INT_MD_LOC)  $(INTRA_CON_ENT)  $(ENR_CTRL_ENT) \
                         $(PIMD_ENT)  $(PIMD_LOC)  $(MATH) \
---------------------------------------
int_nptf_pimd.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(COMM_WRAP) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_LOC) $(INT_MD_ENT)  $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(PIMD_ENT)  $(PIMD_LOC) $(MATH) \
---------------------------------------
int_utils_pimd.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(COMM_WRAP) $(INT_PIMD_LOC) \
---------------------------------------
int_0_to_dt2_pimd.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(INT_MD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(PIMD_ENT) \
                         $(MATH)  $(COMM_WRAP) \
---------------------------------------
int_dt2_to_dt_pimd.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INT_MD_ENT) $(INT_MD_LOC) $(INT_PIMD_LOC) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \

=============================================================
MAKE_INTERFACE
=============================================================

---------------------------------------
parse.o     :            $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_PAR) $(PARSE_ENT) $(PARSE_LOC) \
                         $(SIM_ENT) $(MOL_ENT) $(INTRA_ENT) $(COORD_ENT) \
                         $(COORD_LOC) $(CPEWALD_ENT) $(INTER_ENT) \
                         $(INTER_LOC) $(VPS_ENT) $(LISTS_ENT) $(SCRATCH_ENT) \
                         $(ENR_CPCON_ENT) $(REAL_LOC) $(SAMPL_CLASS_ENT) \
                         $(SAMPL_CP_ENT) $(SAMPL_CP_LOC) $(COORD_CP_ENT) \
                         $(COORD_CP_LOC) $(MATH) $(PIMD_ENT) $(PIMD_LOC) \
                         $(FRND_ENT) $(COMM_ENT) $(COMM_LOC) $(COMM_WRAP) \
                         $(INT_CPMIN_ENT) \
---------------------------------------
zero_bnd.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CP) \
                         $(TYP_CLASS) $(TYP_PAR) \
                         $(PARSE_LOC) \
---------------------------------------
zero_par.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_PAR) \
                         $(TYP_BND) $(TYP_CP) \
                         $(PARSE_LOC) \
---------------------------------------
zero_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(PARSE_LOC) \
---------------------------------------
zero_cp.o     :
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(PARSE_LOC) \
---------------------------------------
interface_hand.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(HANDLE_ENT) \
---------------------------------------
search_base_class.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(SEARCH_LOC) \
                         $(FRND_ENT) \
---------------------------------------
data_base_handle.o    :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(SEARCH_ENT) $(INTER_ENT) $(INTER_LOC) \
                         $(INTRA_LOC)  $(HANDLE_ENT) $(FRND_ENT) \
---------------------------------------
control_sim_params.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(SIM_ENT) $(SIM_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_sim_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND)  \
                         $(TYP_CP) $(TYP_PAR) \
                         $(SIM_LOC) $(FRND_ENT) \
---------------------------------------
set_sim_params.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) $(TYP_CLASS) \
                         $(TYP_BND) $(SIM_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
set_atm_nhc.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
read_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MATH) $(PIMD_LOC) $(COMM_WRAP) \
---------------------------------------
molecule_decomp.o  :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(COORD_ENT) $(COORD_LOC) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
read_hmat.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_ENT) $(HANDLE_ENT) $(FRND_ENT) $(MATH) \
                         $(PIMD_LOC) \
---------------------------------------
mall_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_BND) \
                         $(COORD_LOC) $(FRND_ENT) \
---------------------------------------
mall_pressure.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(COORD_ENT) $(FRND_ENT) \
---------------------------------------
close_intra_params.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
control_intra_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) $(PIMD_LOC) \
                         $(MATH) \
---------------------------------------
control_res_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
fetch_residue.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
---------------------------------------
fetch_resbond_prm.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) \
---------------------------------------
fetch_free_energy_index.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) \
---------------------------------------
fetch_freeze.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
init_intra_params.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
manipulate_res_bonds.o : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
replicate_mol.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
residue_bond.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
control_inter_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTER_ENT) $(INTER_LOC) $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
---------------------------------------
set_inter_dict.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
---------------------------------------
get_clong.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
---------------------------------------
spline_fit.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) \
                         $(INTER_LOC) $(FRND_ENT) \
---------------------------------------
mall_scratch.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(SCRATCH_ENT) $(SCRATCH_LOC) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
control_vnhc_smpl.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
---------------------------------------
control_vx_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COORD_LOC) \
                         $(COMM_WRAP) \
--------------------------------------
control_scale_class.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \
---------------------------------------
proj_vel_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
---------------------------------------
vel_samp_class.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_BND) $(TYP_CLASS) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(MATH) \
---------------------------------------
set_exclude.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_PAR) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) \
---------------------------------------
exl_sort.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_PAR) \
                         $(LISTS_LOC) $(FRND_ENT) \
---------------------------------------
mall_lists.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(LISTS_ENT) $(LISTS_LOC) $(REAL_LOC) $(FRND_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
control_brnch_root_list.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(FRND_ENT) $(LISTS_ENT) $(LISTS_LOC) $(REAL_LOC) \
---------------------------------------
block_intra_lists.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) \
---------------------------------------
class_par_forc_lists.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(FRND_ENT) $(LISTS_ENT) $(LISTS_LOC) \

=============================================================
MAKE_INTERFACE_CP
=============================================================

---------------------------------------
set_wave_params.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CLASS) $(TYP_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) \
---------------------------------------
set_coef_nhc.o    :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
read_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(ENR_CPCON_ENT) $(COORD_CP_ENT) $(HANDLE_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
gen_wave.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) $(TYP_CLASS) \
                         $(FRND_ENT) $(MATH) $(COMM_WRAP) $(ENR_CPCON_LOC) \
                         $(ENR_CP_LOC) $(CPEWALD_LOC) $(COORD_CP_LOC) \
---------------------------------------
mall_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_GEN) $(TYP_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) $(FRND_ENT) \
---------------------------------------
control_set_cp_ewald.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_CP) $(TYP_PAR) \
                         $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) \
                         $(MATH) $(FRND_ENT) \
---------------------------------------
search_base_cp.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(SEARCH_ENT) $(VPS_LOC) $(HANDLE_ENT) \

---------------------------------------
set_cpewald.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(CPEWALD_LOC) $(MATH) $(FRND_ENT) \
---------------------------------------
proj_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) $(TYP_BND) \
                         $(SMPL_CP_LOC) $(SMPL_CP_ENT) $(ENR_CPCON_ENT)
---------------------------------------
set_vps_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
---------------------------------------
vel_samp_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(SMPL_CP_LOC)  $(MATH) \
---------------------------------------
control_vps_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_PAR) \
                         $(TYP_CLASS) $(TYP_BND) \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
control_vc_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
---------------------------------------
control_vcnhc_smpl.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
---------------------------------------
control_scale_cp.o    :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP)  \
                         $(TYP_BND) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \

=============================================================
MAKE_INTERFACE_INTRA
=============================================================

---------------------------------------
set_atm_mask.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_atm_morph.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_atm_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_bend_bnd_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_bend_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_bond_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_intra_dict.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
set_intra_dict_pot.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
---------------------------------------
set_intra_potent.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(FRND_ENT) \
---------------------------------------
intra_coefs.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_mol_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
---------------------------------------
set_onfo_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_res_bond_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
---------------------------------------
set_res_def_params.o   : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_res_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
---------------------------------------
set_res_morph_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
---------------------------------------
set_grp_con_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_tors_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
---------------------------------------
fetch_hydrog_mass.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \

=============================================================
MAKE_INTERFACE_MOL
=============================================================

---------------------------------------
control_mol_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MOL_ENT) $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
---------------------------------------
control_set_mol_params.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
---------------------------------------
set_def_base_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(MOL_LOC) \
---------------------------------------
set_free_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(TYP_PAR) $(TYP_BND) \
                         $(MOL_LOC) $(HANDLE_ENT) \
---------------------------------------
set_mol_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CLASS) \
                         $(TYP_BND) $(TYP_CP) \
                         $(MOL_LOC) $(FRND_ENT) \
---------------------------------------
set_mol_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(MOL_LOC) $(FRND_ENT) $(HANDLE_ENT) \
---------------------------------------
set_user_base_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CLASS) \
                         $(TYP_BND) $(TYP_CP) \
                         $(MOL_LOC) \

=============================================================
MAKE_MAIN
=============================================================

---------------------------------------
main.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(TYP_CP) $(TYP_PAR) \
                         $(MAIN_ENT) $(MAIN_LOC) $(MAIN_CP_LOC) $(PARSE_ENT) \
                         $(COMM_WRAP) $(COMM_ENT) \
---------------------------------------
control_md.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(MAIN_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(OUTPUT_ENT) $(ANAL_MD_ENT) $(SMPL_CLASS_ENT) \
                         $(INTRA_CON_ENT) $(ENR_CTRL_ENT) $(COMM_WRAP) \
---------------------------------------
control_min.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_MIN_ENT) \
                         $(INT_MIN_LOC) $(OUTPUT_ENT) $(ANAL_MD_ENT) \
                         $(SMPL_CLASS_ENT) $(INTRA_CON_ENT) $(ENR_CTRL_ENT) \
---------------------------------------
control_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_PAR) \
                         $(MAIN_LOC) $(MATH) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_ENT) $(OUTPUT_ENT) $(ANAL_MD_ENT) \
                         $(SMPL_CLASS_ENT) $(ENR_CTRL_ENT) $(PIMD_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
control_debug.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(MAIN_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(REAL_LOC) $(COMM_WRAP) \
---------------------------------------
control_debug_pimd.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(MAIN_LOC) $(MATH) $(INT_PIMD_ENT) $(INT_PIMD_LOC) \
                         $(INT_MD_ENT) $(OUTPUT_ENT) $(OUTPUT_LOC) \
                         $(REAL_LOC) $(ENR_CTRL_ENT) $(COMM_WRAP) \

=============================================================
MAKE_MAIN_CP
=============================================================

---------------------------------------
control_debug_cp.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(MAIN_CP_LOC) $(MATH) $(INT_MD_ENT) $(INT_MD_LOC) \
                         $(INT_CP_ENT) $(INT_CP_LOC) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(REAL_LOC) $(INTRA_CON_ENT) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_CP_ENT) $(COMM_WRAP) \
---------------------------------------
control_debug_cp_pimd.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(MAIN_CP_LOC) $(MATH) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(INTRA_CON_ENT) $(REAL_LOC) \
                         $(INT_PIMD_LOC) $(ENR_CTRL_CP_ENT) $(INT_MD_ENT) \
                         $(INT_PIMD_ENT) $(INT_CP_LOC) $(INT_CPPIMD_LOC) \
                         $(COMM_WRAP) \
---------------------------------------
control_cp_min.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CPMIN_ENT) \
                         $(ENR_CPCON_ENT) $(OUTPUT_ENT) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(ANAL_MD_ENT) $(ANAL_CP_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
control_cp_pimd_min.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(MAIN_LOC) $(MAIN_CP_LOC) $(MATH) $(INT_CPMIN_ENT) \
                         $(ENR_CPCON_ENT) $(OUTPUT_ENT) $(OUTPUT_CP_ENT) \
                         $(OUTPUT_CP_LOC) $(ANAL_MD_ENT) $(ANAL_CP_ENT) \
                         $(COMM_WRAP) \
---------------------------------------
CONTROL_CP.O     :       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) $(TYP_BND) \
                         $(MAIN_CP_LOC) $(MATH) $(INT_CP_ENT) $(INT_CP_LOC) \
                         $(OUTPUT_CP_ENT) $(SMPL_CP_ENT) $(SMPL_CLASS_ENT) \
                         $(ANAL_CP_ENT) $(OUTPUT_CP_ENT) $(INT_MD_ENT) \
                         $(INT_MD_LOC) $(INTRA_CON_ENT) $(ENR_CTRL_CP_ENT) \
                         $(ENR_CPCON_ENT) $(COMM_WRAP) \
---------------------------------------
CONTROL_CP_PIMD.O   :    $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(TYP_BND) $(TYP_PAR) \
                         $(MAIN_CP_LOC) $(MATH) $(INT_CPPIMD_ENT) \
                         $(INT_CPPIMD_LOC) $(OUTPUT_CP_ENT)
                         $(OUTPUT_ENT) $(ENR_CTRL_CP_ENT) $(SMPL_CP_ENT) \
                         $(SMPL_CLASS_ENT) $(ANAL_CP_ENT) $(INTRA_CON_ENT) \
                         $(INT_PIMD_ENT) $(INT_PIMD_LOC) $(INT_MD_ENT) \
                         $(INT_CP_LOC) $(ENR_CPCON_ENT) $(PIMD_ENT) \
                         $(COMM_WRAP) \


=============================================================
MAKE_MATH
=============================================================

---------------------------------------
mathlib.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) \
                         $(MATH) \
---------------------------------------
fft_package.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(MATH) $(FRND_ENT) $(ENR_CP_LOC) $(COMM_WRAP) \

=============================================================
MAKE_OUTPUT
=============================================================

---------------------------------------
output_md.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
output_pimd.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(FRND_ENT) $(MATH) \
---------------------------------------
simpavg_md.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) \
---------------------------------------
simpavg_pimd.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(OUTPUT_LOC) $(COMM_WRAP) \
---------------------------------------
simpavg_md_communicate.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(OUTPUT_ENT) $(COMM_WRAP) \
---------------------------------------
get_cell.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_BND) $(TYP_GEN) $(TYP_CLASS) \
                         $(OUTPUT_LOC) \
---------------------------------------
write_gen_header.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(FRND_ENT)  $(OUTPUT_LOC) $(OUTPUT_CP_LOC) \

=============================================================
MAKE_OUTPUT_CP
=============================================================

---------------------------------------
output_cp_min.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(MATH) $(OUTPUT_LOC) $(FRND_ENT) $(COMM_WRAP) \
---------------------------------------
output_cp.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(FRND_ENT) $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \
---------------------------------------
simpavg_cp.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(OUTPUT_LOC) \
                         $(COMM_WRAP) \
---------------------------------------
simpavg_cp_communicate.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) $(TYP_CP) \
                         $(OUTPUT_CP_LOC) $(COMM_WRAP) \
---------------------------------------
output_cp_pimd.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(ENR_CPCON_ENT) \
                         $(FRND_ENT) $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \
---------------------------------------
simpavg_cp_pimd.o   :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) $(TYP_CP) \
                         $(OUTPUT_CP_ENT) $(OUTPUT_CP_LOC) $(FRND_ENT) \
                         $(MATH) $(OUTPUT_LOC) $(COMM_WRAP) \

