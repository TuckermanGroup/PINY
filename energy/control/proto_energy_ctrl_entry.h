/*--------------------------------------------------------------------------*/
/* Energy_control.c */

void energy_control(CLASS *, BONDED *, GENERAL_DATA *);

void energy_control_initial(CLASS *, BONDED *, GENERAL_DATA *);
void energy_control_inter_real(CLASS *, BONDED *, GENERAL_DATA *);
void energy_control_inter_recip(CLASS *, BONDED *, GENERAL_DATA *);
void energy_control_intra(CLASS *, BONDED *, GENERAL_DATA *);
void energy_control_final(CLASS *, BONDED *, GENERAL_DATA *);

void energy_control_surf(CLASS *, BONDED *, GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* Energy_control_pimd.c */

void energy_control_pimd(CLASS *, BONDED *, GENERAL_DATA *);

/*--------------------------------------------------------------------------*/

void test_energy(CLASS *, BONDED *, GENERAL_DATA *);

void test_energy_pimd(CLASS *, BONDED *, GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* Period.c */

void period(int ,double [],double [],double [],CELL *);

void period_one(int ,double *,double *,double *,CELL *);

void period_pimd(int ,double [],double [],double [],
                 double [],double [],double [],CELL *);

void constr_cell_mat(int , int , int , double *);

void chck_constr_cell_mat(int ,int ,int ,double *,double *,int );

/*--------------------------------------------------------------------------*/

