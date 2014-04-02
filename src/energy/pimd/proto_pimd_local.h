/*--------------------------------------------------------------------------*/
/* get_estimators.c */

void get_estimators(CLATOMS_INFO *, CLATOMS_POS *,
                    EWD_SCR *, CELL *, STAT_AVG *,STATEPOINT *);          

void get_pimd_spread(CLATOMS_INFO *, CLATOMS_POS *, 
                     double *,COMMUNICATE *);

void prim_quantum_KE(CLATOMS_INFO *, CLATOMS_POS *,PTENS *,
                     double *);

/*--------------------------------------------------------------------------*/
/* /pimd/transform_cnt.c */

void path_integral_init(CLATOMS_INFO *,CLATOMS_POS *,CLATOMS_TRAN *,
                        GHOST_ATOMS *,SIMOPTS *,ATOMMAPS *,
                        COMMUNICATE *);

void convert_pimd_force_cent(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void convert_pimd_mode_cent(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void convert_pimd_pos_cent(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void assign_mode_bead_cent(CLATOMS_INFO *, CLATOMS_POS *);

/*--------------------------------------------------------------------------*/
/* /pimd/transform_stg.c */

void convert_pimd_force_stag(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void convert_pimd_mode_stag(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void convert_pimd_pos_stag(CLATOMS_INFO *, CLATOMS_POS *,CLATOMS_TRAN *);

void stage_1_part(double *,double *,double ,double ,
            double ,int ,int *,int *,double *,int *);

void assign_mode_bead_stag(CLATOMS_INFO *, CLATOMS_POS *);

void convert_pimd_mode_cent_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void convert_pimd_pos_cent_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void convert_pimd_force_cent_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void convert_pimd_mode_stag_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void convert_pimd_pos_stag_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void convert_pimd_force_stag_par(CLATOMS_INFO *, CLATOMS_POS *,
                                  CLATOMS_TRAN *, COMMUNICATE *);

void pimd_trans_comm_fwd(CLATOMS_INFO *,CLATOMS_POS *,
                      CLATOMS_TRAN *,COMMUNICATE *,int);

void pimd_trans_comm_bck(CLATOMS_INFO *,CLATOMS_POS *,
                      CLATOMS_TRAN *,COMMUNICATE *,int);

void get_kin_vir(CLASS *);
