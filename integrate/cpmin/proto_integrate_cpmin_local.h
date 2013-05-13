/*--------------------------------------------------------------------------*/
/* diag_ovlap */

double diag_ovlap(int ,double ,double,double *,double *,double *,double *,int );

void get_diag_cp_hess(CP *,int,CELL *,double);


void cp_shuffle_prim(double *,double *,int ,double *,double *,int *,
                     double *,double *, double *, int,
                     VEL_SAMP_CP *,
                     CP_COMM_STATE_PKG *);

void realloc_DIIS_mem(int **,double **,double **,
                      double **,double **,double **,double **,
                      int ,int ,int ,int );

void init_alloc_DIIS_mem(int **,double **,double **,
                      double **,double **,double **,double **,
                      int ,int ,int ,int );

void shift_DIIS_hist(double *,double *,double *,double *,int,int);

void setup_DIIS_eqs(double *,double *,double *,double *,
                    COMMUNICATE *,int ,int ,int ,int *,int ,int );

void init_alloc_atm_DIIS_mem(int **,double **,double **,
                             double **,double **,double **,
                             double **,double **,double **,
                             int ,int ,int );

void realloc_atm_DIIS_mem(int **,double **,double **,
                          double **,double **,double **,
                          double **,double **,double **,
                          int ,int ,int );

void shift_atm_DIIS_mem(double *,double *,double *,
                        double *,double *,double *,
                        int ,int );

void setup_atm_DIIS_eqs(double *,double *,
                        double *,double *,double *,
                        int ,int ,int );

void act_hess_inv_on_grad(double *,double *,double *,
                          double *,double *,double *,
                          double *,double *,double *,
                          double *,double *,double *,int );

void act_hess_inv_on_grad_diag(double *,double *,double *,
                               double *,double *,double *,
                               double *,double *,double *,
                               double *,double *,double *,int );


void line_min_cp(CLASS *,BONDED *,GENERAL_DATA *,CP *,
                 double *,double *,int ,double );

