/*--------------------------------------------------------------*/
/* Group communication                                          */

void build_communicate_groups(COMMUNICATE *, int);

void build_cp_comm_pkg(CP *,MPI_Comm);

MPI_Comm build_grp_comm_outer(int , int , int , int , int *,int *, MPI_Comm );

MPI_Comm build_grp_comm_inner(int , int , int , int , int *,int *, MPI_Comm );

/*--------------------------------------------------------------*/
/* Interface communication                                      */

void communicate_class_info(CLASS *,GENERAL_DATA *,CLASS_PARSE *);

void communicate_class_list(CLASS *,GENERAL_DATA *,CLASS_PARSE *,int );

void comm_cell_info(CELL *,MPI_Comm);

void comm_clatoms_info(CLATOMS_INFO *,MPI_Comm);

void communicate_class_data(CLASS *,GENERAL_DATA *,CLASS_PARSE *,int);

void comm_statepoint(STATEPOINT *,MPI_Comm);

void comm_cell_data(CELL *,double *,int *,MPI_Comm);

void comm_clatoms_data(CLATOMS_INFO *,int,int,MPI_Comm);

void comm_ghost_atoms_data(GHOST_ATOMS *,MPI_Comm);

void comm_therm_data(THERM_INFO *,MPI_Comm);

void communicate_cp_info(CP *,CP_PARSE *,MPI_Comm,int);

void comm_cp_parse_info(CP_PARSE *, MPI_Comm);

void comm_cpewald(CPEWALD *,MPI_Comm);

void comm_pseudo(PSEUDO *,MPI_Comm,int);

void comm_cpopts(CPOPTS *,MPI_Comm);

void communicate_cp_list(CP *,int,MPI_Comm);

void comm_cpcoeffs_info(CPCOEFFS_INFO *,MPI_Comm);

void comm_cptherm_info(CPTHERM_INFO *,MPI_Comm);

void comm_cpewald_info(CPEWALD *,MPI_Comm);

void comm_class_parse_info(CLASS_PARSE *,MPI_Comm);

void mall_class(CLASS *,GENERAL_DATA *,CLASS_PARSE *,int);

void mall_bond(BONDED *,NULL_INTER_PARSE *);

void comm_cell_info(CELL *,MPI_Comm);

void comm_clatoms_info(CLATOMS_INFO *,MPI_Comm);

void comm_ghost_info(GHOST_ATOMS *,MPI_Comm);

void comm_atommaps_info(ATOMMAPS *,MPI_Comm);

void comm_therm_info_class(THERM_INFO *,MPI_Comm);

void comm_therm_info_bead(THERM_INFO *,MPI_Comm);

void comm_nbr_list_info(NBR_LIST *,MPI_Comm);

void comm_verlist_info(VERLIST *,MPI_Comm);

void comm_lnklist_info(LNKLIST *,MPI_Comm);

void comm_interact_info(INTERACT *,MPI_Comm);

void comm_interact_data(INTERACT *,MPI_Comm);

void comm_ewald_info(EWALD *,MPI_Comm);

void comm_ewald_data(EWALD *,MPI_Comm);

void comm_pme_info(PART_MESH *,MPI_Comm);

void comm_part_mesh_data(PART_MESH *,EWALD *,MPI_Comm);

void comm_baro(BARO *,MPI_Comm);

void comm_par_rahman(PAR_RAHMAN *,MPI_Comm);

void comm_ptens(PTENS *,MPI_Comm);

void comm_bond_info(BOND *,MPI_Comm);

void communicate_bond_data(BONDED *,NULL_INTER_PARSE *,MPI_Comm);

void comm_grp_bond_con_info(GRP_BOND_CON *,MPI_Comm);

void comm_grp_bond_con_data(GRP_BOND_CON *,MPI_Comm);

void comm_grp_bond_watts_info(GRP_BOND_WATTS *,MPI_Comm);

void comm_grp_bond_watts_data(GRP_BOND_WATTS *,MPI_Comm);

void comm_bond_free_info(BOND_FREE *,MPI_Comm);

void comm_bend_info(BEND *,MPI_Comm);

void comm_bend_free_info(BEND_FREE *,MPI_Comm);

void comm_bend_bnd_info(BEND_BND *,MPI_Comm);

void comm_tors_info(TORS *,MPI_Comm);

void comm_tors_free_info(TORS_FREE *,MPI_Comm);

void comm_onfo_info(ONFO *,MPI_Comm);

void comm_ecor_info(ECOR *,MPI_Comm);

void comm_null_inter_parse_info(NULL_INTER_PARSE *,MPI_Comm);

void communicate_bond_info(BONDED *,NULL_INTER_PARSE *,MPI_Comm);

void comm_bend_bnd_data(BEND_BND *,MPI_Comm);

void comm_bond_data(BOND *,MPI_Comm);

void comm_bond_free_data(BOND_FREE *, MPI_Comm);

void comm_bend_data(BEND *,MPI_Comm);

void comm_bend_free_data(BEND_FREE *, MPI_Comm);

void comm_nbr_data(NBR_LIST *,MPI_Comm);

void communicate_general_data(GENERAL_DATA *,CP *,MPI_Comm);

void comm_simopts(SIMOPTS *,MPI_Comm);

void communicate_class(CLASS *,CLASS_PARSE *,int );

void comm_clatoms(CLATOMS_INFO *,MPI_Comm);

void communicate_cp_data(CP *,int ,int,MPI_Comm,int );

void comm_cpcoeffs_data(CPCOEFFS_INFO *,MPI_Comm);

void comm_cpconstrnt(CPCONSTRNT *,MPI_Comm);

void communicate_cp(CP *,MPI_Comm);

void comm_tors_data(TORS *,MPI_Comm);

void comm_rbar_sig_free_info(RBAR_SIG_FREE *,MPI_Comm );

void comm_rbar_sig_free_data(RBAR_SIG_FREE *,MPI_Comm );

void comm_bond_free_data(BOND_FREE *,MPI_Comm);

void comm_bend_data(BEND *,MPI_Comm);

void comm_bend_free_data(BEND_FREE *,MPI_Comm);

void comm_bend_bnd_data(BEND_BND *,MPI_Comm);

void comm_tors_free_data(TORS_FREE *,MPI_Comm);

void comm_onfo_data(ONFO *,MPI_Comm);

void comm_ecor_data(ECOR *,MPI_Comm);

void comm_bond(BOND *,MPI_Comm);

void comm_grp_bond_con_data(GRP_BOND_CON *,MPI_Comm);

void comm_ghost_atoms(GHOST_ATOMS *,MPI_Comm);

void comm_atommaps(ATOMMAPS *,int,MPI_Comm);

void comm_ewald(EWALD *,int,MPI_Comm);

void comm_part_mesh(PART_MESH *part_mesh,int,MPI_Comm );

void comm_timeinfo(TIMEINFO *,MPI_Comm);

void comm_vel_samp_class(VEL_SAMP_CLASS *,MPI_Comm,int);

void comm_vel_samp_cp(VEL_SAMP_CP *,MPI_Comm,int);

void comm_energy_ctrl(ENERGY_CTRL *,MPI_Comm);

void comm_bend(BEND *,MPI_Comm);

void comm_bend_bnd(BEND_BND *,MPI_Comm);

void comm_tors(TORS *,MPI_Comm);

void comm_onfo(ONFO *,MPI_Comm);

void comm_ghost_atoms_data(GHOST_ATOMS *,MPI_Comm);

void comm_ghost_atoms_list(GHOST_ATOMS *,MPI_Comm);

void comm_class_parse(CLASS_PARSE *,int,MPI_Comm );

void comm_class_parse_data(CLASS_PARSE *,int,MPI_Comm );

void comm_ensopts(ENSOPTS *,MPI_Comm);

void comm_minopts(MINOPTS *,MPI_Comm);

void comm_intra_scr(INTRA_SCR *,MPI_Comm);

void comm_constrnt_info(CONSTRNT *,MPI_Comm);

void comm_cpewald(CPEWALD *,MPI_Comm);

void comm_cpscr_info(CPSCR *,MPI_Comm);

void comm_bond_list(BOND *,MPI_Comm world);

void comm_grp_bond_con_list(GRP_BOND_CON *,MPI_Comm);

void comm_grp_bond_watts_list(GRP_BOND_WATTS *,MPI_Comm);

void comm_bend_list(BEND *,MPI_Comm );

void comm_bend_bnd_list(BEND_BND *,MPI_Comm );

void comm_tors_list(TORS *,MPI_Comm );

void comm_onfo_list(ONFO *,MPI_Comm );

void comm_cpewald_list(CPEWALD *,MPI_Comm);

void comm_nbr_list_data(NBR_LIST *,MPI_Comm );

void comm_clatoms_list(CLATOMS_INFO *,MPI_Comm);

void comm_atommaps_list(ATOMMAPS *,int,MPI_Comm);

void comm_ewald_list(EWALD *,int ,MPI_Comm );

void comm_part_mesh_list(PART_MESH *,int,MPI_Comm );

void comm_class_parse_list(CLASS_PARSE *,int ,MPI_Comm);

void comm_stat_avg_info(STAT_AVG *,MPI_Comm);

void comm_null_inter_parse_list(NULL_INTER_PARSE *, MPI_Comm);

void comm_filenames(FILENAMES *,MPI_Comm);

void comm_communicate_info(COMMUNICATE *,MPI_Comm);

void comm_cpopts_data(CPOPTS *,MPI_Comm );

void comm_surface_info(SURFACE *, MPI_Comm );

void comm_surface_data(SURFACE *, MPI_Comm );

                              

                              


