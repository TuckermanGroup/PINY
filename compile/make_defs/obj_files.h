#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
#
#  Object files for PI_MD.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================

ANALYSIS_FILES        = analysis_md.o analysis_pimd.o prelim_analysis.o \
                        calcul_gr.o harmonic_analysis.o \
                        calcul_ickt_iso_md.o calcul_iikt_iso_cmd.o \
                        calcul_msqd_atm.o calcul_vovt_atm.o 
ANALYSIS_CP_FILES     = analysis_cp.o analysis_cp_pimd.o
CP_ENERGY_FILES       = test_energy_cp.o \
                        cp_energy_control.o \
                        cp_ks_energy.o cp_energy_ee_rho.o cp_energy_ee_rho_ke.o \
                        cp_energy_ee_grad_rho.o cp_energy_eext.o \
                        cp_energy_eext_nonloc.o cp_energy_eext_nonloc_gh.o\
                        cp_coef_force_tau_fun.o xc_functionals.o \
                        constraint_control_cp.o cp_con.o cp_orth_rot_utils.o \
                        cp_con_utils.o energy_control_elec.o \
                        orth_rot_control_cp.o cp_transpose.o \
                        control_contract_rho.o control_spread_rho.o
ENERGY_INTER_FILES    = energy_control.o energy_control_initial.o \
                        energy_control_intra.o \
                        energy_control_inter_real.o \
                        energy_control_inter_recip.o \
                        energy_control_surf.o \
                        energy_control_final.o test_energy.o \
                        force_control.o force_nolst.o \
                        force_verlst.o force_npol.o force_lnklst.o period.o \
                        nbr_list_control.o verlist_control.o verlist_create.o \
                        make_lnk_lst.o make_lnk_map.o lnk_lst_dis.o \
                        ewald3d.o ewald3d_both.o ewald3d_self_bgr.o \
                        ewald3d_pme.o surf_pot.o
ENERGY_INTRA_FILES    = bond.o bond_watts.o bend.o bend_bnd.o tors.o \
                        onefour.o ecorr.o rbar_sigma.o
ENERGY_INTRA_CON_FILES= ghost_control.o constraint_control.o bond_con.o \
                        bond_con_rolli.o bond_both.o bond_con_rollf.o \
                        group_bond_con_21.o group_bond_con_23.o \
                        group_bond_con_33.o group_bond_con_43.o \
                        group_bond_con_46.o proj_com_min.o \
                        group_bond_con_rolli_21.o group_bond_con_rolli_23.o \
                        group_bond_con_rolli_33.o group_bond_con_rolli_43.o \
                        group_bond_con_rolli_46.o \
                        group_bond_con_rollf_21.o group_bond_con_rollf_23.o \
                        group_bond_con_rollf_33.o group_bond_con_rollf_43.o \
                        group_bond_con_rollf_46.o group_bond_con_util.o
ENERGY_PIMD_FILES     = test_energy_pimd.o energy_control_pimd.o \
                        test_energy_cp_pimd.o \
                        transform_cnt.o transform_stg.o pimd_estimators.o \
                        transform_cnt_par.o transform_stg_par.o \
                        cp_energy_control_pimd.o control_pimd_trans.o
FRIEND_FILES          = cmalloc.o friend_lib.o
INTEGRATE_CP_FILES    = int_NVE_cp.o int_NVT_cp.o int_NPTI_cp.o int_NPTF_cp.o \
                        min_STD_cp.o min_CG_cp.o min_DIIS_cp.o move_atm.o \
                        int_utils_cp.o int_0_to_dt2_cp.o int_dt2_to_dt_cp.o \
                        shuffle_states.o cp_gauss.o
INTEGRATE_PIMD_FILES  = int_NVT_pimd.o int_NPTI_pimd.o int_NPTF_pimd.o \
                        int_utils_pimd.o int_0_to_dt2_pimd.o \
                        int_dt2_to_dt_pimd.o
INTEGRATE_CPPIMD_FILES = int_NVT_cp_pimd.o int_NPTI_cp_pimd.o \
                        int_NPTF_cp_pimd.o int_utils_cp_pimd.o

INTEGRATE_MD_FILES    = int_NVE.o int_NVE_res.o int_NVT.o int_NVT_res.o \
                        int_NPTI.o int_NPTI_res.o int_NPTF.o int_NPTF_res.o \
                        int_utils.o min_STD.o min_CG.o move_pos_box.o \
                        int_0_to_dt2.o int_dt2_to_dt.o move_vel_vbox.o
INTERFACE_FILES       = parse.o zero_class.o zero_bnd.o zero_par.o zero_cp.o \
                        interface_hand.o search_base_class.o \
                        data_base_handle.o \
                        control_sim_params.o set_sim_dict.o set_sim_params.o \
                        set_atm_NHC.o read_coord.o molecule_decomp.o \
                        read_hmat.o mall_coord.o \
                        control_surf_params.o set_surf_dict.o \
                        control_inter_params.o spline_fit.o \
                        set_inter_dict.o get_clong.o mall_scratch.o \
                        control_vnhc_smpl.o control_vx_smpl.o \
                        proj_vel_class.o samp_vel_class.o set_exclude.o \
                        exl_sort.o mall_lists.o mall_pressure.o \
	                block_intra_lists.o control_scale_class.o \
                        control_brnch_root_list.o class_par_forc_lists.o
INTERFACE_CP_FILES    = set_wave_params.o set_coef_NHC.o read_coef.o \
                        mall_properties.o \
                        gen_wave.o mall_coef.o control_scale_cp.o \
                        control_set_cp_ewald.o set_cp_ewald.o \
                        search_base_cp.o proj_vel_cp.o \
                        set_vps_dict.o samp_vel_cp.o control_vps_params.o \
                        weight_node_gauss_hermite.o \
                        control_vc_smpl.o control_vcnhc_smpl.o
INTERFACE_INTRA_FILES = close_intra_params.o control_intra_params.o \
                        control_res_params.o fetch_residue.o \
                        fetch_resbond_prm.o fetch_free_energy_index.o \
                        fetch_freeze.o init_intra_params.o \
                        manipulate_res_bonds.o replicate_mol.o residue_bond.o \
                        set_atm_mask.o set_atm_morph.o set_atm_params.o \
                        set_bend_bnd_params.o set_bend_params.o \
                        set_bond_params.o set_intra_dict.o \
                        set_intra_dict_pot.o set_intra_potent.o \
                        set_mol_name_params.o set_onfo_params.o \
                        set_res_bond_params.o set_res_def_params.o \
                        set_res_name_params.o set_res_morph_params.o \
                        set_grp_con_params.o set_tors_params.o intra_coefs.o \
                        fetch_hydrog_mass.o
INTERFACE_MOL_FILES   = control_mol_params.o control_set_mol_params.o \
                        set_base_file_params.o set_free_params.o \
                        set_mol_dict.o set_mol_params.o \
                        set_surf_params.o 
MAIN_FILES            = main.o control_md.o control_pimd.o \
                        control_debug.o control_debug_pimd.o control_min.o auto_exit.o
MAIN_CP_FILES         = control_cp.o control_cp_min.o control_debug_cp.o \
                        control_cp_pimd.o control_debug_cp_pimd.o \
                        control_cp_pimd_min.o 
MATH_FILES            = mathlib.o blas_wrappers.o fft_package.o \
                        fft_create_package.o fft_generic.o
OUTPUT_FILES          = output_md.o simpavg_md.o output_pimd.o simpavg_pimd.o \
                        get_cell.o write_gen_header.o simpavg_md_communicate.o
OUTPUT_CP_FILES       = output_cp.o output_cp_min.o output_min.o simpavg_cp.o \
                        simpavg_cp_communicate.o simpavg_cp_pimd.o \
                        output_cp_pimd.o 
COMMUNICATE_FILES     = communicate_wrappers.o \
                        mall_class.o com_interface.o comm_class_info.o \
                        comm_class_data.o \
                        comm_bond_info.o comm_cp_info.o mall_bond.o \
                        comm_bond_data.o comm_class_list.o comm_cp_data.o \
                        communicate_options.o \
                        communicate_test_energy_pimd.o \
                        communicate_utils_pimd.o simpavg_pimd_communicate.o \
                        pimd_trans_comm.o build_communicate_groups.o \
                        communicate_simpavg_cp_pimd.o \
                        cp_communicate_coord_class.o \
                        communicate_output_pimd.o \
                        communicate_output_cp_pimd.o \
                        control_group_communicators.o 

OBJS = $(ANALYSIS_FILES) $(ANALYSIS_CP_FILES) $(INTEGRATE_CPPIMD_FILES)\
       $(CP_ENERGY_FILES) $(ENERGY_INTER_FILES) $(ENERGY_INTRA_FILES)\
       $(ENERGY_INTRA_CON_FILES) \
       $(ENERGY_PIMD_FILES) $(FRIEND_FILES) $(INTEGRATE_CP_FILES) \
       $(INTEGRATE_MD_FILES) $(INTEGRATE_PIMD_FILES) $(INTERFACE_FILES) \
       $(INTERFACE_CP_FILES) $(INTERFACE_INTRA_FILES) $(INTERFACE_MOL_FILES)\
       $(MAIN_FILES) $(MAIN_CP_FILES) $(MATH_FILES) $(OUTPUT_FILES) \
       $(OUTPUT_CP_FILES) $(SPEC_FILES) $(COMMUNICATE_FILES)

OBJS_CARE = make_lnk_lst.o make_lnk_map.o lnk_lst_dis.o parse.o \
            zero_bnd.o zero_par.o zero_class.o zero_cp.o interface_hand.o \
            search_base_class.o data_base_handle.o control_sim_params.o \
            set_sim_dict.o set_sim_params.o set_atm_NHC.o read_coord.o \
            read_hmat.o mall_coord.o mall_pressure.o close_intra_params.o \
            control_intra_params.o control_res_params.o fetch_residue.o \
            fetch_resbond_prm.o fetch_free_energy_index.o fetch_freeze.o \
            init_intra_params.o manipulate_res_bonds.o replicate_mol.o \
            residue_bond.o control_inter_params.o set_inter_dict.o \
            get_clong.o mall_scratch.o control_vnhc_smpl.o control_vx_smpl.o \
            control_scale_class.o proj_vel_class.o samp_vel_class.o \
            set_exclude.o exl_sort.o mall_lists.o block_intra_lists.o \
            set_wave_params.o set_coef_NHC.o read_coef.o mall_coef.o \
            control_set_cp_ewald.o set_recip.o search_base_cp.o \
            set_cp.o proj_vel_cp.o set_vps_dict.o samp_vel_cp.o \
            control_vps_params.o control_vc_smpl.o control_vcnhc_smpl.o \
            control_scale_cp.o set_atm_mask.o set_atm_morph.o \
            set_atm_params.o set_bend_bnd_params.o set_bend_params.o \
            set_bond_params.o set_intra_dict.o set_intra_dict_pot.o \
            set_intra_potent.o intra_coefs.o set_mol_name_params.o \
            set_onfo_params.o set_res_bond_params.o set_res_def_params.o \
            set_res_name_params.o set_res_morph_params.o set_grp_con_params.o \
            set_tors_params.o fetch_hydrog_mass.o control_mol_params.o \
            control_set_mol_params.o set_def_base_params.o set_free_params.o \
            set_mol_dict.o set_mol_params.o set_user_base_params.o \
            output_cp_min.o

OBJS_GRP  = group_bond_con.o group_bond_con_21.o group_bond_con_43.o \
            group_bond_con_rolli.o group_bond_con_rolli_21.o \
            group_bond_con_rolli_43.o group_bond_con_rollf.o \
            group_bond_con_rollf_21.o group_bond_con_rollf_43.o




