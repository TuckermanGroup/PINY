void prelim_vovt_atm_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_vovt_atm_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);
void output_vovt_atm_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);

void prelim_msqd_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_msqd_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);
void output_msqd_pimd(CLASS *,GENERAL_DATA *,ANALYSIS *);

void prelim_iikt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_iikt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void output_iikt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void test_box_matrix(CLASS *,GENERAL_DATA *,ANALYSIS *);
void get_kvec_1d_iso(CLASS *,GENERAL_DATA *,ANALYSIS *,int );
void get_kvec_3d_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void test_kvec_param(CLASS *,GENERAL_DATA *,ANALYSIS *);

void prelim_ickt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_ickt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void output_ickt_iso(CLASS *,GENERAL_DATA *,ANALYSIS *);
void test_box_matrix_ickt(CLASS *,GENERAL_DATA *,ANALYSIS *);
void get_kvec_1d_iso_ickt(CLASS *,GENERAL_DATA *,ANALYSIS *,int );
void get_kvec_3d_iso_ickt(CLASS *,GENERAL_DATA *,ANALYSIS *);

void prelim_gr(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_full_gr(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_inter_gr(CLASS *,GENERAL_DATA *,ANALYSIS *);
void test_same_molecule(int, int, int);
void output_gr(CLASS *,GENERAL_DATA *,ANALYSIS *);
void get_and_print_sq_from_gr(CLASS *,GENERAL_DATA *,ANALYSIS *);

void normal_mode_water(CLASS *,GENERAL_DATA *,ANALYSIS *);
void static_charge_density(CLASS *,GENERAL_DATA *,ANALYSIS *);
void hf_stuff(CLASS *,GENERAL_DATA *,ANALYSIS *);

void calcul_freqs(CLASS *,GENERAL_DATA *,ANALYSIS *);

void prelim_sim_sc_op(CLASS *,GENERAL_DATA *,ANALYSIS *);
void get_sim_sc_op(CLASS *,GENERAL_DATA *,ANALYSIS *);
void correl_sim_sc_op(CLASS *,GENERAL_DATA *,ANALYSIS *);
void output_sim_sc_op(CLASS *,GENERAL_DATA *,ANALYSIS *);

