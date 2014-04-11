void control_pimd_trans_force(CLASS *,GENERAL_DATA *);

void control_pimd_trans_mode(CLASS *,GENERAL_DATA *);

void control_pimd_trans_pos(CLASS *,GENERAL_DATA *);

void mode_energy_control(CLASS *,GENERAL_DATA *);

void pimd_pvten_forc_corr(CLASS *, GENERAL_DATA *,double *);

void get_vir_press(CLASS *, GENERAL_DATA *general_data, double *);

void spread_coord(CLATOMS_INFO *,CLATOMS_POS *,
                  double *,double *,double *,
                  int *,int *,double *,
                  ATOMMAPS *,int );

