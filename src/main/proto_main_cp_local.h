/*---------------------------------------------------------------------*/
void control_debug_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void control_debug_cp_pimd(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void rolf_cp_vectors(CP *,int ,int );

void hurl_cp_scalars(CP *);

/*---------------------------------------------------------------------*/
/* cp/control_cp.c control_cp_min   */

void control_cp_min(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);

void control_cp_pimd_min(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);

void control_cp_pimd(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);

void control_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *,ANALYSIS *);

void prelim_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void prelim_cp_pimd(CLASS *,BONDED *,GENERAL_DATA *,CP *);

void check_coef_grad_mag(CP *,SIMOPTS *,
                         double *,double *,int *, int *, double, 
                         int,int,STAT_AVG *);

void check_atm_grad_mag(CLASS *,GENERAL_DATA *,double *,int *, int *, 
                        double);

/*---------------------------------------------------------------------*/
