/*--------------------------------------------------------------------------*/
/* Energy_control.c */

void long_range_corr(int, double *, double ,PTENS *, FOR_SCR *,
                     INTERACT *, ENERGY_CTRL *, double);

void check_cutoff_clus(int ,double *,double *,double *, int ,double );

void mix_coul_corr(CLATOMS_INFO *,CLATOMS_POS *,INTRA_SCR *,EXCL *,INTERACT * ,
                   CELL *,PTENS *,double *,double );

/*--------------------------------------------------------------------------*/
/* energy_control_final.c                                                   */

void box_force_iso(double *,double *,double *,double *,double, double, int);

void box_force_flex(double *, double *,double *,double ,double , int , CELL *);

/*--------------------------------------------------------------------------*/
/* get_estimators.c */

void get_estimators(CLATOMS_INFO *, CLATOMS_POS *,
                    EWD_SCR *, CELL *, STAT_AVG *,STATEPOINT *);          
/*--------------------------------------------------------------------------*/
/* test_energy.c */

void test_force(CLASS *, BONDED *, GENERAL_DATA *,int , double *,double *, 
                double *);
void communicate_test_energy(double *,double *,double *,
                                      double *, MPI_Comm,int);


/*--------------------------------------------------------------------------*/
/* test_energy_pimd.c */

void test_force_pimd(CLASS *, BONDED *, GENERAL_DATA *,int , double *,double *, 
                double *);

void communicate_test_energy_pimd(double *,double *,double *,
                                              double *, MPI_Comm);

/*--------------------------------------------------------------------------*/
/* Period.c */

void check_cutoff(int ,int ,double *,double );


/*--------------------------------------------------------------------------*/
