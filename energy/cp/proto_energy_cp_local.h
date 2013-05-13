/*===========================================================================*/
/* cp_ks_energy.c */

void cp_ks_energy_ctrl(CP *,int ,EWALD *,EWD_SCR *,CELL *, CLATOMS_INFO *,
                       CLATOMS_POS *,ATOMMAPS *,STAT_AVG *,
                       PTENS *,SIMOPTS *,FOR_SCR *);

void cp_ks_energy_hybrid(CP *,int ,EWALD *,EWD_SCR *,CELL *, CLATOMS_INFO *,
                         CLATOMS_POS *,ATOMMAPS *,STAT_AVG *,
                         PTENS *,SIMOPTS *,FOR_SCR *);

void cp_ks_energy_full_g(CP *,int ,EWALD *,EWD_SCR *,CELL *, CLATOMS_INFO *,
                         CLATOMS_POS *,ATOMMAPS *,STAT_AVG *,
                         PTENS *,SIMOPTS *,FOR_SCR *);

void cp_boundary_check(CELL *,CLATOMS_INFO *,CLATOMS_POS *,double );

void cp_dual_check(CELL *,CLATOMS_INFO *,CLATOMS_POS *,int *,double );

/*===========================================================================*/
/* cp_energy_eext.c */

void control_cp_eext_recip(CLATOMS_INFO *,CLATOMS_POS *,CPCOEFFS_INFO *,
                           CPCOEFFS_POS *,CPEWALD *,CPSCR *,
                           CPOPTS *,PSEUDO *,EWD_SCR *,
                           ATOMMAPS *,CELL *,EWALD *,PTENS *,
                           double *,double *,COMMUNICATE *,FOR_SCR *,
                           int ,PARA_FFT_PKG3D * );


void control_ewd_loc(CLATOMS_INFO *,CLATOMS_POS *,
                     CELL *, PTENS *, EWALD *, CPEWALD *, 
                     CPSCR *, PSEUDO *, EWD_SCR *,CPOPTS *,
                     ATOMMAPS *, double *,double *,COMMUNICATE *,
                     FOR_SCR *,int ,int );

void control_ewd_loc_pme(CLATOMS_INFO *, CLATOMS_POS *,
                     CELL *, PTENS *, EWALD *, CPEWALD *, 
                     CPSCR *, PSEUDO *, EWD_SCR *, CPOPTS *, 
                     ATOMMAPS *, double *, double *,COMMUNICATE *,
                     FOR_SCR *,int ,int ,PARA_FFT_PKG3D * );

void get_vpsnow(int *,int ,double ,double ,double ,
                double *,double *,double *,double *,
                double *,int *,int ,int ,int );


void ewald3d_selfbgr_cp(CLATOMS_INFO *,EWALD *, PTENS *,
                        double, double *, double *,int );

void get_ak2_sm(CPEWALD *,CELL *);

void real_hess_approx(double *,double *,double *,
                      double *,double *,double *,
                      double *,double *,double *,double *,int ,int);

void reduce_cp_hess_stuff(double *,double *,double *,
                           double *,double *,double *,
                           int ,int , MPI_Comm );

void reduce_cp_atm_forc(int ,double *,double *,double *,
                        double *,double *,double *, MPI_Comm , int );

void get_vpslong(int ,double *,double ,double *,double ,double );

void get_dvpslong(int ,double *,double ,double *,double ,double );

void get_atm_hess_recip(double ,double ,double , double *,double *,double *,
                        double *,double *, double *, double *,double *,
                        double *,double *, double *, double *, double *, 
                        double ,double , double ,double ,
                        int , int ,int , int , int ,int *);

/*===========================================================================*/
/* cp_energy_eext_nonloc.c */

void control_ewd_non_loc(CLATOMS_INFO *, CLATOMS_POS *,CPCOEFFS_INFO *,
                         CPCOEFFS_POS *,
                         CELL *, PTENS *, CPEWALD *,
                         CPSCR *, PSEUDO *, EWD_SCR *, CPOPTS *, 
                         ATOMMAPS *, COMMUNICATE *,FOR_SCR *);

void control_nonloc_gh(CLATOMS_INFO *,CLATOMS_POS *,
                       CPCOEFFS_INFO *,CPCOEFFS_POS *,
                       CELL *, PTENS *,CPEWALD *,
                       CPSCR *, PSEUDO *, EWD_SCR *,
                       CPOPTS *, ATOMMAPS *,
                       COMMUNICATE *,FOR_SCR *,
                       double );


void control_nlmat(CLATOMS_INFO *,CPCOEFFS_INFO *,
                   CPCOEFFS_POS *,CPSCR *,CPOPTS *,PSEUDO *,
                   EWD_SCR *,ATOMMAPS *,int ,int,int *,
                   double ,double ,double ,double ,int,YLM_CONS *);

void control_nlmat_gh(CLATOMS_INFO *,CPCOEFFS_INFO *,
                      CPCOEFFS_POS *,CPSCR *,CPOPTS *,PSEUDO *,
                      EWD_SCR *,ATOMMAPS *,int ,int ,
                      double ,double ,double ,double ,
                      int ,YLM_CONS *,double );



void get_nlmat(int ,int ,int ,int ,int ,int ,int ,
               int ,int ,double *,double *,
               double *,double *,
               double *,double *,
               double *,double *,double *,
               double *,double *,double *,
               double ,double ,double ,double ,double ,
               double ,double *,double *);


void get_nlmat_gh(int ,int ,int ,int ,int ,
               int ,int ,int ,int ,
               double *,double *,
               double *,double *,
               double *,double *,
               double *,double *,double *,
               double *,double *,double *,
               double ,double ,double ,double ,double ,
               double *,double *,double ,
               double ,double );


void get_nlmat_gh_hess(int ,int ,int ,int ,int ,
               int ,int ,int ,int ,
               double *,double *,
               double *,double *,
               double *,double *,
               double *,double *,double *,
               double *,double *,double *,
	       double *,double *,double *,
	       double *,double *,double *,
	       double *,double *,double *,
	       double *,double *,double *,
               double ,double ,double ,double ,double ,
               double *,double *,double ,
               double ,double );


void get_nlmat_pv(int ,int ,int ,int ,int ,int ,int ,
                  int , int ,double *,double *,
                  double *,double *,double *,double *,
                  double *,double *,
                  double *,double *,
                  double *,double *,double *,
                  double *,double *,double *,
                  double *,double *,double *,
                  double *,double *,double *,
                  double *,double *,double *,
                  double *,double ,double ,double ,double ,
                  double ,double ,double ,
                  double ,double ,double ,double ,
                  double ,double *,double *,double *);


void get_nlmat_hess(int ,int ,int ,int ,int ,
         int ,int ,int ,int ,
         double *,double *,
         double *,double *,
         double *,double *,
         double *,double *,
         double *,double *,
         double *,double *,double *,
         double *,double *,double *,
         double *,double *,double *,
         double *,double *,double *,
         double *,double *,double *,
         double *,double ,double ,double ,
         double ,double ,double ,double *,double *);


void get_nlmat0(int ,int ,
                int ,int ,int , int, int,int,
                double *,double *,double *,
                double *,double *,double );

void get_nlmat0_gh(int ,int ,int ,
                   int ,int ,int ,int ,int ,
                   double *,double *,double );


void getnl_pot_pv_fatm(CLATOMS_INFO *,CLATOMS_POS *,
                  CELL *cell,CPCOEFFS_INFO *,
                  CPSCR *,EWD_SCR *,CPOPTS *,
                  PSEUDO *,ATOMMAPS *,
                  double *, int ,double *);

void getnl_pot_pv_fatm_gh(CLATOMS_INFO *,CLATOMS_POS *,
                          CELL *,CPCOEFFS_INFO *,
                          CPSCR *,EWD_SCR *,CPOPTS *,
                          PSEUDO *,ATOMMAPS *,
                          double *,int ,double *,double *,int );

void sumnl_pot_pv_fatm_hess(int ,int ,int ,int ,int ,int ,int ,int , int ,
                            int *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,
                            double *,double *,double *,
                            double *,double *,double *,
                            double *,double *,double *,
                            double *,double *,double *,
                            int ,int ,double *, double *);


void getnl_fcoef(CLATOMS_INFO *,CLATOMS_POS *,CPCOEFFS_INFO *,
                  CPCOEFFS_POS *,CPSCR *,EWD_SCR *,CPOPTS *,
                  PSEUDO *,CPEWALD *,ATOMMAPS *,
                  CELL *, int ,double *,FOR_SCR *);

void getnl_fcoef_gh(CLATOMS_INFO *,CLATOMS_POS *,
                    CPCOEFFS_INFO *,CPCOEFFS_POS *,
                    CPSCR *,EWD_SCR *,CPOPTS *,
                    PSEUDO *,CPEWALD *,ATOMMAPS *,
                    CELL *, int ,double *,FOR_SCR *,
                    double ,double *,int );


void control_nlfcoef(CLATOMS_INFO *,
                     CPCOEFFS_INFO *,
                     CPCOEFFS_POS *,
                     CPSCR *,CPOPTS *,PSEUDO *,
                     EWD_SCR *,ATOMMAPS *,
                     int ,int ,int *,
                     double ,double ,double ,double ,double ,
                     int ,YLM_CONS *);

void control_nlfcoef_gh(CLATOMS_INFO *,
                        CPCOEFFS_INFO *,
                        CPCOEFFS_POS *,
                        CPSCR *,CPOPTS *,PSEUDO *,
                        EWD_SCR *,ATOMMAPS *,
                        int ,int ,
                        double ,double ,double ,double ,double ,
                        int ,YLM_CONS *,double ,double *,int );



void get_nlfor(int ,int ,int ,int ,int ,
               int ,int ,int ,int ,
               double *,double *,double *,double *,
               double *,double *,double *,
               double ,double ,
               double ,double );

void get_nlfor_gh(int ,int ,int ,int ,
                  int ,int ,int ,int ,int ,
                  double *,double *,double *,double *,
                  double *,double *,double *,double ); 


void get_nlhess(int ,int ,int ,int ,
                int ,int ,int ,int ,int ,
                double *,double *,double *,double *,
                double *,double *,double ,double ,
                double ,double );

void get_nlfor0(int ,int ,int ,
                int ,int ,int , int , int,
                double *,double *,double *,double *,
                double *,double *,double ,double );


void get_nlfor0_gh(int ,int ,int ,int ,int ,
                   int ,int ,int ,
                   double *,double *,
                   double  ,double ,double );


void get_nlhess0(int ,int ,int ,
                 int ,int ,int , int , int,
                 double *,double *,double *,
                 double *,double *,double ,double );


void get_vpsnorm(double *,double *,double *,int *,
                 int ,int ,int ,int ,int , int , int);

void get_ylm(double ,double ,double ,double ,
             double *,double *,double *,double *,
             double *,double *,double *,double *,
             YLM_CONS *);

void vps_atm_list(PSEUDO *, CELL *, CLATOMS_POS *, CLATOMS_INFO *, ATOMMAPS *, 
                  EWD_SCR *, FOR_SCR *,int ,int );

void control_vps_atm_list(PSEUDO *, CELL *, CLATOMS_POS *, CLATOMS_INFO *, ATOMMAPS *, 
                          EWD_SCR *, FOR_SCR *,int ,int );

void non_loc_chng_ord(CLATOMS_POS *, CLATOMS_INFO *,ATOMMAPS *, PSEUDO *, 
                      EWD_SCR *, FOR_SCR *,int );

void non_loc_restore_ord(CLATOMS_POS *, CLATOMS_INFO *,
                         ATOMMAPS *, PSEUDO *,EWD_SCR *, FOR_SCR *);

 double get_jl(double ,double ,int );
/*====================================================================*/
/* cp_energy_ee_rho.c */


void cp_rho_calc_hybrid(CPEWALD *,CPSCR *, CPCOEFFS_INFO *,EWALD *,
                        CELL *,double *, double *,
                        int ,int , double *, double *,double *,
                        double *, double *,double *, double *,double *,
                        double *,int ,int ,int ,int ,int ,COMMUNICATE *, 
                        PARA_FFT_PKG3D *, PARA_FFT_PKG3D *,
                        PARA_FFT_PKG3D *, PARA_FFT_PKG3D *, PARA_FFT_PKG3D *);

void cp_ke_dens_calc_hybrid(CPEWALD *,CPSCR *,CELL *,double *, double *,
                            int ,int ,double *,int ,int ,int ,
                            COMMUNICATE *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *);

void cp_ke_dens_calc_full_g(CPEWALD *,CPSCR *,CELL *,double *, double *,
                            int ,int ,double *,int ,int ,int ,
                            COMMUNICATE *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *,
                            PARA_FFT_PKG3D *);

void construct_elf(double *,double *,
                   double *,double *,double *,double *,double,
                   PARA_FFT_PKG3D *);

      
void cp_rho_calc_full_g(CPEWALD *,CPSCR *, CPCOEFFS_INFO *,EWALD *,
                        CELL *,double *, double *,
                        int ,int , double *, double *,double *,
                        double *, double *,double *, double *,double *,
                        double *,int ,int ,int ,int ,int ,COMMUNICATE *, 
                        PARA_FFT_PKG3D *, PARA_FFT_PKG3D *, PARA_FFT_PKG3D *);
      
void coef_force_control(CPOPTS *,CPCOEFFS_INFO *,CPCOEFFS_POS *,
                               CPSCR *,EWALD *,CPEWALD *,CELL *,STAT_AVG *,
                               char *,double *,double ,double ,int ,int ,COMMUNICATE *,
                               CP_COMM_STATE_PKG *,CP_COMM_STATE_PKG *,
                               PARA_FFT_PKG3D *, PARA_FFT_PKG3D *,
                               PARA_FFT_PKG3D *, PARA_FFT_PKG3D *,
                               PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,int );

void cp_get_vks(CPOPTS *,CPSCR *,CPEWALD *, EWALD *,
                COMMUNICATE *, CP_COMM_STATE_PKG *, CP_COMM_STATE_PKG *,
                STAT_AVG *, double *,CELL *, char *, double *,double ,double ,int ,
                int ,int ,int ,int , int , PARA_FFT_PKG3D *, PARA_FFT_PKG3D *,
                PARA_FFT_PKG3D *, PARA_FFT_PKG3D *, int );

void coef_force_calc_hybrid(CPEWALD *,int ,double *,double *, 
                             double *,double *,double *,double *,double *,double *,
                             double *,double *,double *,double *,double *,
                             double *,double *,int ,double *,
                             COMMUNICATE *,int ,int ,int ,int ,int ,
                             PARA_FFT_PKG3D *);

void coef_force_calc_full_g(CPEWALD *,int ,int ,int ,
                             double *,double *,double *,double *,
                             double *,double *,double *,double *,
                             double *,double *,double *,double *,double *,
                             double *,double *,int ,double *,
                             COMMUNICATE *,int ,int ,int ,int ,int ,
                             PARA_FFT_PKG3D *);

void coef_force_tau_fun_hybrid(CPEWALD *,int ,double *,double *, 
                               double *,double *,double *,double *,
                               double *,double *,double *,double *,
                               double *,int ,double *,
                               COMMUNICATE *,int ,int ,int ,
                               PARA_FFT_PKG3D *);

void coef_force_tau_fun_full_g(CPEWALD *,int ,double *,double *, 
                               double *,double *,double *,double *,
                               double *,double *,double *,double *,
                               double *,int ,double *,
                               COMMUNICATE *,int ,int ,int ,
                               PARA_FFT_PKG3D *);

void cp_vpsi(double *,double *,int );

void cp_pack_vks(double *,double *,int );

/*====================================================================*/
/* cp_energy_ee_rho.c */

void control_grad_rho(CPEWALD *,CPSCR *,EWALD *,
                      double *,double *,
                      double *,double *,double *,
                      double *,
                      double *,double ,int ,int ,
                      PARA_FFT_PKG3D *);

void  grad_corr_lda(CPSCR *,CPEWALD *,EWALD *,CELL *,double *,
                    double *,double ,CPOPTS *,double ,
                    double *,int ,int ,int ,int ,PARA_FFT_PKG3D *,int);

void  grad_corr_lsda(CPSCR *,CPEWALD *,EWALD *,CELL *, double *,
                     double *,double ,CPOPTS *,double ,
                     double *,int ,int ,int , int ,int ,int ,PARA_FFT_PKG3D *);

void yon_hither(CPSCR *,CPEWALD *,EWALD *,double *,double *,
                double *,double *,double *,int ,int ,PARA_FFT_PKG3D *);

void ptens_gga(double *,double ,double ,double ,double ,double ,
                double , double ,double );

void ptens_gga_lyp(double *,double ,double ,double ,double ,
                   double ,double , double ,double ,double );

void ptens_gga_lap(double *,CPSCR *,CPEWALD *,EWALD *,double *,double ,int ,
                   PARA_FFT_PKG3D *);


void  get_igdot(CPSCR *,EWALD *,CPEWALD *,double *,int ,int , int ,int );

void get_igcoef_hybrid(double *,double *, double *,double *,
                       CPEWALD *,double *,int ,int ,int );

void get_igcoef_full_g(double *,double *,double *,double *,
                       CPEWALD *,double *,int ,int ,int ,int ,int );

void  get_g2times(CPSCR *,EWALD *,double *,int *,int *,int *,int ,
                  int ,int ,int ,int);


void get_igrho(CPSCR *,int *,int *,int *,
               double *,double *,double *,int ,int ,int );

void get_igigrho(CPSCR *,EWALD *,double *,double *,double *,
                 int ,int ,int,int );
      

void get_g2rho(CPSCR *,int *,int *,int *,
               double *,double *,double *,int ,int ,int );

/*====================================================================*/
/* xc_functionals.c */

/*---------------------------------------------------------------------*/
/* GGA Exchange functionals */

void becke_gcx_lda(double ,double ,
                     double *,double *,double *,double);

void becke_gcx_lsda(double ,double ,
                     double *,double *,double *,double);

void fila_1x_lsda(double ,double ,double *,double *,double *);


void fila_2x_lsda(double ,double ,double *,double *,double *);


void fila_1x_lda(double ,double ,double *,double *,double *);


void fila_2x_lda(double ,double ,double *,double *,double *);


void pw91_gcx(double ,double ,double *,double *,double *);

void pbe_gcx(double ,double ,double *,double *,double *,double,double);

void rpbe_gcx(double ,double ,double *,double *,double *);

void brx89_lda(double , double , double , double ,double *,
                double *,double *,double *,double *);

void brx89_lsda(double , double , double , double ,double *,
                double *,double *,double *,double *);

void brx2K_lda(double , double , double , double ,double *,
                double *,double *,double *,double *);

void brx2K_lsda(double , double , double , double ,double *,
                double *,double *,double *,double *);

double newt_raph(double,int *);

double bisect(double);

void tau1_lda(double , double , double , double ,double *,
                double *,double *,double *,double *,int );

void tau1_lsda(double , double ,double ,double ,
               double ,double , double ,double ,
               double *,double *,double *,
               double *,double *,
               double *,double *,
               double *,double *,int ,int );

/*---------------------------------------------------------------------*/
/* GGA Correlation functionals */

void pw91_gcc(double ,double ,double *,double *,double *);

void pbe_gcc(double ,double ,double *,double *,double *,double ,double);

void lyp_gcc(double ,double ,double *,double *,double *);

void lypm1(double ,double ,double ,double *,double *,double *,double *);


void lyp_lsda(double ,double ,double ,double ,double ,
              double ,double ,double , double *,double *,
              double *,double *,double *,double *,double *,double *,double *);
void pbe_gcc_lsda(double ,double ,double ,double ,double ,double ,double ,
		  double ,double ,double ,double *,double *, 
		  double *,double *,double *,double ,double);

/*---------------------------------------------------------------------*/
/* Debug GGA functionals */

void debug97x(double ,double ,double *,double *,double *);

void debug_lap1(double ,double ,double ,double *,double *,double *,double *);

void debug_tau(double ,double ,double ,
               double *,double *,double *,double *);

/*---------------------------------------------------------------------*/
/* LDA/LSDA functionals */

void  excpot_pz_lda(double *,double *,double *,double *,int ,
                    int ,double ,int ,int ,double *,int );


void  excpot_pade_lda(double *,double *,double *,double *,int ,
                      int ,double ,int ,int ,double *,int );

void  excpot_pw_lda(double *,double *,double *,double *,
                    int ,int ,double , int ,int ,double *,int );

void excpot_pz_lsda(double *,double *,double *,double *,double *,
                    double *,int ,int ,double ,int ,int ,double *,int );

void excpot_pw_lsda(double *,double *,double *,double *,double *,
		    double *,int ,int ,double ,int ,int ,double *,int );

void pw_c_lsda(double ,double ,double *,double *,double *);

void pw_ec(double, double, double, double, double, double, double, 
	   double *, double *);

void excpot_pade_lsda(double *,double *,double *,double *,double *,
                      double *,int ,int ,double ,int ,int ,double *,int );


/*===========================================================================*/
/* control_spread_rho.c */

void control_spread_rho(CPSCR *,double *,CELL *,double ,int ,int ,
                        PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,int );

void sngl_pack_rho_dual_ser_pme(double *,double *,CELL *,double ,int ,int ,
                                CPSCR *,PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void sngl_pack_rho_dual_par_pme(CPSCR_WAVE *,double *,CELL *,
                                double ,int ,int ,CPSCR_DUAL_PME *,
                                PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void sngl_upack_rho_dual_ser_pme(double *,double *,CELL *,int ,CPSCR *cpscr,
                                 PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void sngl_upack_rho_dual_par_pme(CPSCR_WAVE *,double *, CELL *,
                                 int ,CPSCR_DUAL_PME *,
                                PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);

void control_contract_rho(CPSCR *,double *,CELL *,int ,int ,
                          PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,int ); 

void make_dual_pme_wghts(CPSCR_DUAL_PME *,CELL *,double ,int ,PARA_FFT_PKG3D *,
                       PARA_FFT_PKG3D *);

void make_pme_para_dual_map(int ,CPSCR_WAVE *,CPSCR_DUAL_PME *,
                            PARA_FFT_PKG3D *,PARA_FFT_PKG3D *);


/*===========================================================================*/
