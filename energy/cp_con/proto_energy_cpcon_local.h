/*--------------------------------------------------------------------------*/
/* CP constraint routines:cp_con.c */

void cp_shake(double *,double *,  int , int ,
              double *,double *,int , 
              double *,double *,int ,
              double ,double ,
              CPSCR_OVMAT *, 
              int *,double *,double *,int ,
              int *,CP_COMM_STATE_PKG *);

void cp_gram_schmidt_par(double *,double *,int ,
                         double *,int *,
                         double *, double *,
                         CP_COMM_STATE_PKG *);

void cp_shake_mass(double *,double *,int , int ,
                   double *,double *,int ,
                   double *,double *,int ,
                   double *,
                   double ,double ,
                   CPSCR_OVMAT *, 
                   int *,double *,double *,int ,
                   int *,
                   CP_COMM_STATE_PKG *);

void cp_shake_norb(double *,double *, int ,int ,
                   double *,double *,int ,
                   double *,double *,int ,
                   double ,int *,double *,
                   double *,int,
                   CPSCR_OVMAT *,
                   CP_COMM_STATE_PKG *);

void cp_rattle(double *,double *,int , int ,
               double *,double *,int ,
               CPSCR_OVMAT *,
               int *,double *,
               CP_COMM_STATE_PKG *,int );


void cp_rattle_mass(double *,double *,int , int ,
                    double *,double *,int ,
                    double *,double *,
                    double , double *,
                    CPSCR_OVMAT *, 
                    int *,double *, double *,
                    int *,
                    CP_COMM_STATE_PKG *,int);

void cp_rattle_norb(double *,double *,int , int ,
                    double *,double *,int ,
                    double *,double *,
                    int *,double *, 
                    double *,int,
                    CPSCR_OVMAT *,
                    CP_COMM_STATE_PKG *);




/*--------------------------------------------------------------------------*/
/* CP Rotation-orthogonalization routines:cp_orth_rot_utils.c */


void control_cp_gram_schmidt(double *,double *cimag,int ,
                          double *,double *,
                          double *,int *,
                          CP_COMM_STATE_PKG *);

void cp_gram_schmidt_scalar(double *,double *,double *,
                            double *,int *,int ,int );

void cp_gram_schmidt_par(double *,double *,int ,
                         double *,int *,double *, double *,
                         CP_COMM_STATE_PKG *);

void cp_normalize(double *,double *,int ,
                  double *,int *,
                  double *,double *,
                  CP_COMM_STATE_PKG *);

void cp_normalize_all(double *,double *,int ,
                  double *,double *,int , 
                  double *,int *,
                  double *,double *,
                  CP_COMM_STATE_PKG *);

void cp_add_ksmat_force(double *,double *,
                        int ,int ,
                        double *,double *,
                        int ,int ,
                        double *,double *,
                        int *,int ,int ,double *,
                        CP_COMM_STATE_PKG *);

void cp_rotate_gen_nonortho(double *,double *, 
                            int ,int *,
                            double *,int *,
                            double *,double *,
                            CP_COMM_STATE_PKG *);

void cp_rotate_coef_ortho(double *,double *,
                          int ,int *,
                          double *,double *,
                          double *,
                          double *,double *,
                          double *,int *,
                          double *,double *,
                          CPSCR_OVMAT *,
                          CP_COMM_STATE_PKG *);


void cp_rotate_ks_basis(double *,double *,
                        int , int ,
                        double *,double *,
                        int , int ,
                        double *,double *,
                        int , int ,
                        int *,double *,double *,
                        double *, double *,
                        double *,double *,
                        CPSCR_OVMAT *,
                        CP_COMM_STATE_PKG *,double *);


void cp_condiag_ksmat(double *,double *,
                      int , int ,
                      double *,double *,
                      int , int ,
                      double *,double *,
                      double *, double *,
                      double *,double *,
                      int *,
                      CP_COMM_STATE_PKG *,double *);

void cp_construct_orth_rotmat(double *,double *,
                              int ,int ,
                              double *,double *,
                              double *,
                              double *,int *,
                              double *,double *,
                              CPSCR_OVMAT *,
                              CP_COMM_STATE_PKG *);

void cp_rotate_vector(double *,double *,int , 
                      double *,int *,
                      double *, double *,
                      CP_COMM_STATE_PKG *);

void cp_rotate_all(double *,double *,int ,
                   double *,double *,int ,
                   double *,double *,int ,
                   double *,int *ioff,
                   double *, double *,
                   CP_COMM_STATE_PKG *);

void cp_rotation_prim(double *,int ,
                      double *,int ,
                      double *,double ,double ,
                      double ,int *,
                      CP_COMM_STATE_PKG *);

void cp_rotate_prim_par(double *,int ,
                        double *,int ,
                        double *,double ,double ,
                        double ,CP_COMM_STATE_PKG *);

void  rotate_occ_orth_shuffle(double *,double *,int );

void rotate_occ_shuffle(double *,double *,int ,
                        double *,double *,
                        double *,double *,
                        int *,CP_COMM_STATE_PKG *);


void occ_sort(int , double *);


/*--------------------------------------------------------------------------*/
/* Utilities:cp_con_utils.c                                                 */ 

void cp_ovlap_mat_same(double *,double *, int ,
                       double *,double *,double ,int *,
                       CP_COMM_STATE_PKG *);

void cp_par_ovlap_same(double *,int ,
                       double *,double ,
                       CP_COMM_STATE_PKG *);

void cp_ovlap_mat_diff(double *,double *,int ,
                       double *,double *,int ,
                       double *,double *,double ,int *,
                       CP_COMM_STATE_PKG *);

void cp_par_ovlap_diff(double *,int ,
                       double *,int ,
                       double *,double ,
                       CP_COMM_STATE_PKG *);

void cp_iter_mat_solve_shake(double *,double *,double *,
                             double *,double *,double *,double *,
                             double *,double *,int *,int *,
                             double );

void cp_iter_mat_solve_shake_mass(double *,double *,double *,
                                  double *,double *,double *,
                                  double *,double *,double *,
                                  int *,int *, double );

void cp_iter_mat_solve_rattle_mass(double *,double *,
                                   double *,double *,double *,
                                   double *,double *,double *,
                                   int *,int *, double );

/*--------------------------------------------------------------------------*/
/* Transpose routines                                                       */ 


void cp_transpose_fwd(double *,double *,int *,
                      double *,double *,
                      CP_COMM_STATE_PKG *);

void cp_transpose_bck(double *,double *,int *,
                      double *,double *,
                      CP_COMM_STATE_PKG *);

void cp_therm_transpose_fwd(CPTHERM_POS *,CPTHERM_INFO *,
                            double *,double *,int ,
                            CP_COMM_STATE_PKG *,int *,int);

void cp_therm_transpose_bck(CPTHERM_POS *,CPTHERM_INFO *,
                            double *,double *,int ,
                            CP_COMM_STATE_PKG *,int *,int);

void cp_therm_tran_pack(CPTHERM_POS *,CPTHERM_INFO *,
                          int ,int );

void cp_therm_tran_unpack(CPTHERM_POS *,CPTHERM_INFO *,
                          int ,int );

void cp_transpose_fwd_prim(double *,int *,double *,
                           int ,int ,int ,
                           int ,int ,
                           int ,int ,
                           int ,int ,MPI_Comm );

void cp_transpose_bck_prim(double *,int *,double *,
                           int ,int ,int ,
                           int ,int ,
                           int ,int ,
                           int ,int ,MPI_Comm );


/*--------------------------------------------------------------------------*/
