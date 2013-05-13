/*--------------------------------------------------------------------------*/
/* min_STD_cp.c */

void min_STD_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int );

void min_DIIS_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int,int );

/*--------------------------------------------------------------------------*/
/* min_CG_cp.c */

void min_CG_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *,int );

/*--------------------------------------------------------------------------*/

void move_atm_std(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

void move_atm_cg(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

void move_atm_diis(CLASS *,BONDED *,GENERAL_DATA *,CP *,int);

void displace_atm(CLASS *,BONDED *,GENERAL_DATA *,CP *,double,int *,int *);

void assign_hessian(CLASS *,BONDED *,GENERAL_DATA *,CP *,double,int ,int );
/*--------------------------------------------------------------------------*/
/* shuffle_states.c */

void cp_shuffle_states(CP *,int);

