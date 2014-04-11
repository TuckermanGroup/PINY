void control_set_cp_ewald(SIMOPTS *,CELL *,CPCOEFFS_INFO *,
                          EWALD *, CPEWALD *,CP_PARSE *,
                          double *,double *,double *,
                          EWD_SCR *,int ,int ,
                          double *,int ,PART_MESH *,ECOR *,int,
                          int,int,int);


void control_fft_pkg(PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                     PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                     PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                     PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                     EWALD* ,CPEWALD *,
                     PART_MESH *, CPCOEFFS_INFO *,
                     COMMUNICATE *,int ,int,
                     double *,int ,int,int );


void get_coul_clus_corr(EWALD *,COMMUNICATE *,CELL *,
                        int ,int ,double ,double ,int ,int ,
                        double *,int ,int);


void get_coul_2D_corr(EWALD *,COMMUNICATE *,CELL *,int ,int ,int,
                      PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,double *,int );

void get_coul_1D_corr(EWALD *,COMMUNICATE *,CELL *,int ,int ,int,
                      PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,double *,int );



void init_cp_dual_pme_bw(CPSCR_DUAL_PME *,EWALD *,CPEWALD *,int ,
                         PARA_FFT_PKG3D *);

void init_cp_atom_pme_bw(CPSCR_ATOM_PME * ,EWALD * ,CPEWALD *,
                         PARA_FFT_PKG3D * );

