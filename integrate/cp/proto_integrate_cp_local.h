/*--------------------------------------------------------------------------*/
/* /cp/int_NVT.c */

void apply_c_nhc(CPTHERM_INFO *,CPTHERM_POS *,CPCOEFFS_INFO *,CPCOEFFS_POS *,
                 CPSCR *,int,COMMUNICATE * );

void apply_c_nhc_massiv(CPTHERM_INFO *,CPTHERM_POS *,CPCOEFFS_INFO *,
                        CPCOEFFS_POS *, CPSCR *,int,int ,int );

/*--------------------------------------------------------------------------*/
/* half step cp integrators */

 void int_0_to_dt2_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);

 void int_dt2_to_dt_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *);

/*--------------------------------------------------------------------------*/
/* /cp/int_utils.c */

void get_cpke(CPCOEFFS_INFO *,CPCOEFFS_POS *,STAT_AVG *, int ,int);

void nhc_cp_potkin(CPTHERM_INFO *,CPTHERM_POS *,STAT_AVG *,int,MPI_Comm);

void init_cp_NHC(CPTHERM_INFO *,CPTHERM_POS *,CPCOEFFS_INFO *,CPCOEFFS_POS *,
                 CPSCR *,int, MPI_Comm ,int ,int);

void nhc_cp_potkin_massiv(CPTHERM_INFO *,CPTHERM_POS *,STAT_AVG *);

/*--------------------------------------------------------------------------*/

void add_gauss_force(CPCOEFFS_INFO *,CPCOEFFS_POS *,
                     CPSCR_OVMAT *,CPOPTS *);

void apply_lgauss(CPCOEFFS_INFO *,CPCOEFFS_POS *,
                  CPSCR_WAVE *,CPSCR_OVMAT *,CPOPTS *,double );

/*--------------------------------------------------------------------------*/
/* CP Isokinetic routines */

double HF(double , double , double );

double HFDOT(double , double , double );

void get_isok_ab_facts(double *,double *,double *,double *,
                       double *,int ,int *,int ,int ,double ,
                       int ,MPI_Comm ,double *, double *);






