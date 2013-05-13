void mall_coef(CP *,SIMOPTS *,int);

void assign_coef(CPCOEFFS_INFO *,CPCOEFFS_POS *,double *,double *,double *,
                 double *,int, MPI_Comm, int, int, int,COMMUNICATE *,int,int);

void assign_nhc_coef_vel(CPTHERM_INFO *,CPTHERM_POS *,double *,int,
                         MPI_Comm,int,int );

void bess_trans(double *,int ,double ,double *,
             double *,int ,double *,int ,double *);

void get_gpsi(double ,int ,
              double **,double **,double **,double **,
              double *,double ,double ,int );

void  splin_btrans(int *,double **,double **,
                    double **,double **,
                    double ,double ,double *,
                    double *,int *,char *,int ,
                    int ,int ,MPI_Comm  );

void  fit_spline(double *,double *,double *,
                 double *,double *,int );




