/*--------------------------------------------------------------------------*/
/* Force_controlers.c */

void force_nolst(CLATOMS_INFO *,CLATOMS_POS *,
                 FOR_SCR *,ATOMMAPS *,CELL *,PTENS *,INTERACT *,
                 ENERGY_CTRL *, NBR_LIST *, EXCL *, INTRA_SCR *, double *,
                 double *,double *,CLASS_COMM_FORC_PKG *);

void force_verlst(CLATOMS_INFO *,CLATOMS_POS *,
                  FOR_SCR *,ATOMMAPS *,CELL *,PTENS *,INTERACT *,
                  ENERGY_CTRL *,int [],list_int *,int [],INTRA_SCR *,double *,
                  double *,double *);


void force_lnklst(CLATOMS_INFO *,CLATOMS_POS *,
                  FOR_SCR *,ATOMMAPS *,CELL *,PTENS *,INTERACT *,
                  ENERGY_CTRL *,int ,int ,int ,int ,list_int [],
                  int ,double [],int [],int [],int [],int [],
                  EXCL *, INTRA_SCR *, double *,double *,double *,
                  CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* Force_npol.c */

void force_npol(CLATOMS_INFO *,CLATOMS_POS *,
       FOR_SCR *,ATOMMAPS *,
       CELL *,PTENS *,INTERACT *,
       ENERGY_CTRL *,INTRA_SCR *,double *,double *,double *, int ,int );

void vspl_fetch(int ,double [], double [],
                int [],double [],double [],
                double [],double [],
                double [],double [],double [],double []);

void vspl_fetch_new(int intact_now,double [],double [],
               int [],double [],double [],
               double [],double [],
               double [],double [],double [],double [],
               double [],double []);

void pten_sum(int ,int , double [], double [],
              double [], double [], double [],
              double [], double [], double [],
              double [], double [], double [],
              double [], double [], double [],
              double [], double [], double []);


/*--------------------------------------------------------------------------*/
/* verlist_control.c                                                        */

void verlist_control(CLATOMS_INFO *,CLATOMS_POS *,
                     NBR_LIST *,EXCL *,
                     ATOMMAPS *,CELL *,INTRA_SCR *,
                     FOR_SCR *,TIMEINFO *,
                     INTERACT *,int,CLASS_COMM_FORC_PKG *);

void nolst_ver_gen(CLATOMS_INFO *,CLATOMS_POS *,
                   FOR_SCR *,ATOMMAPS *,
                   CELL *,INTERACT *,
                   TIMEINFO *, NBR_LIST *, EXCL *,
                   INTRA_SCR *,int *,int *,int,CLASS_COMM_FORC_PKG *);

void nolst_ver_gen_root(CLATOMS_INFO *,CLATOMS_POS *,
                   FOR_SCR *,ATOMMAPS *,
                   CELL *,INTERACT *,
                   TIMEINFO *, NBR_LIST *, EXCL *,
                   INTRA_SCR *,int *,int *,
                   int ,
                   CLASS_COMM_FORC_PKG *);

void lnk_ver_gen_root(CLATOMS_INFO *,CLATOMS_POS *,
                      FOR_SCR *,ATOMMAPS *,
                      CELL *,INTERACT *,
                      TIMEINFO *,NBR_LIST *,
                      EXCL *,INTRA_SCR *, int *, 
                      int *,int ,
                      CLASS_COMM_FORC_PKG *);

void lnk_ver_gen(CLATOMS_INFO *,CLATOMS_POS *,
                 FOR_SCR *,ATOMMAPS *,
                 CELL *,INTERACT *,
                 TIMEINFO *,NBR_LIST *,EXCL *,
                 INTRA_SCR *,int *,int *,int,CLASS_COMM_FORC_PKG *);

void lnk_ver_gen_pad(CLATOMS_INFO *,CLATOMS_POS *,
                     FOR_SCR *,ATOMMAPS *,
                     CELL *,INTERACT *,
                     TIMEINFO *,NBR_LIST *,EXCL *,
                     INTRA_SCR *,int *,int *, int *,int,
                     CLASS_COMM_FORC_PKG *);

void lnk_ver_gen_pad_root(CLATOMS_INFO *,CLATOMS_POS *,
                          FOR_SCR *,ATOMMAPS *,
                          CELL *,INTERACT *,
                          TIMEINFO *,NBR_LIST *,EXCL *,
                          INTRA_SCR *,int *,int *, int *,int,
                                           CLASS_COMM_FORC_PKG *);

void ver_list_mem_check(int ,int *,int ,int ,int ,double ,int *,int);

void verpad_list_mem_check(int ,int *,int , int ,int , int *,int);

/*--------------------------------------------------------------------------*/
/* verlist_create.c                                                         */

void verlist_create(CLATOMS_INFO *,CLATOMS_POS *,
                    INTRA_SCR *,FOR_SCR *,
                    ATOMMAPS *,CELL *,NBR_LIST *,
                    TIMEINFO *, INTERACT *,
                    int ,int *,int *,int );

/*--------------------------------------------------------------------------*/
/* make_lnk_map.c                                                           */

void make_lnk_map(CLATOMS_INFO *,CLATOMS_POS *,
                  NBR_LIST *,TIMEINFO *,
                  CELL *,EXCL *);

void lnk_map_cnt(CLATOMS_INFO *,CLATOMS_POS *,
                 NBR_LIST *,TIMEINFO *,CELL *,int);


/*--------------------------------------------------------------------------*/
/* make_lnk_lst.c                                                           */

void make_lnk_lst(CLATOMS_INFO *,CLATOMS_POS *, NBR_LIST *,CELL *,
                  FOR_SCR *,INTRA_SCR *,TIMEINFO *, EXCL *,INTERACT *,int);
               
void lnkcell_sort(int , int *);

/*--------------------------------------------------------------------------*/
/* lnk_cell_dis.c                                                           */

void lnk_cell_dis(int , int , int , 
                  int , int , int , 
                  double *, double [], int );

void max_cut_excl(CLATOMS_INFO *,CLATOMS_POS *,
                  CELL *,NBR_LIST *, EXCL *);

void max_cut_part(INTERACT *interact, NBR_LIST *nbr_list);

void period_lnk(int ,double [],double [],double [],
                double [], double [], int );

void make_hmat_lnk(CLATOMS_INFO *,CLATOMS_POS *,
                   CELL *,NBR_LIST *,
                   double *);

/*--------------------------------------------------------------------------*/
/* Nbr_list_control.c */

void check_ver_list(CLATOMS_INFO *,CLATOMS_POS *,
                    NBR_LIST *,INTERACT *,int *);

void get_cut_skin(double *,double *,double *,double *,double *,double *,
                  double,double,double,int);

void find_max_brnch_root_dist(CLATOMS_POS *,int ,int *,
                              int *,double *,int);

/*--------------------------------------------------------------------------*/
/* force npol utilities                                                 */


void npol_posdata_fetch(double *,double *,double *,int *,
                        int *,int *,
                        double *,double *,double *,
                        double *,double *,double *,
                        int *,
                        double *,double *,double *,
                        double *,double *,
                        int ,int *,int ,int ,
                        int ,int *,int ,int ,
                        int );

void npol_spldata_fetch(double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        int *,int *,
                        double ,double ,
                        int ,int ,
                        int ,
                        int ,int ,int ,int ,
                        int );

void npol_dist_calc(double *,double *,double *,
                    double *,double *,double *,
                    double *,
                    int ,int ,int ,CELL *);

void npol_pot_calc(int ,double *,
                   double *,int *,
                   double *,double *,
                   double *,double *,
                   int ,double ,
                   int ,int ,int , double *,
                   double *,double *,
                   double *,double *,double *);

void npol_force_calc(int ,double *,
                        double *,int *,
                        double *,
                        double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,double *,
                        double *,
                        double *,double *,
                        double *,double *,
                        double *,
                        int ,int ,
                        int ,double ,double ,
                        int ,
                        int );

void npol_hess_calc(CLATOMS_INFO *,CLATOMS_POS *,CELL *,double *);

void npol_force_reduc(int ,int *,
                      int ,double *,double *,
                      double *,double *,
                      double *,double *,
                      int *,int *);
       
void npol_press_reduc(int,int,double *,
                      double *,
                      double *,double *,double *,
                      double *,double *,double *,
                      double *,double *,double *,
                      double *,double *,double *,
                      double *,double *,double *,
                      double *,double *,
                      double *,double *,
                      int , int ,int ,double ,double );
       
void npol_vir_reduc(int ,int ,
                    double *,double *,
                    double *,double,
                    int,int *,
                    double *,double *,
                    double *,int *,
                    int *);
       

void npol_shave_skin(int *,double *,
                    double *,double *,
                    double *,double *,
                    int *,int *,
                    int,int *,
                    double *,double *,
                    double *,double *,
                    int *,int,int,int,int);

void pten_sum_roll(int ,int , double [],double [],
                   double [],  double [],  double [],
                   double [],  double [],  double [],
                   double [], 
                   double [],  double [],  double [],
                   double [], double [], double []);

double dsum1_npol(int,double *);

void npol_vspl_fetch(int ,double [], double [],
                    int [],double [],double []);



/*--------------------------------------------------------------------------*/
/* vercreate utilities                                                      */



void vercreate_posdata_fetch(double *,double *,double *,
                             int *,int *,int *,
                             double *,double *,double *,
                             double *,double *,double *,
                             int *,
                             double *,double *,double *,
                             int ,int *,int ,int ,
                             int ,int ,int ,CELL *);

void vercreate_list_count_fill(
                          int ,int *,int ,
                          double *,double *,double *,
                          int ,int *,int *,
                          int *,
                          double *,double *,
                          double *,double *,
                          double *,
                          int *,list_int *,int *,
                          int *,list_int *,int *,
                          int );

void vercreate_list_count_fill_root(
                          int ,int *,int ,
                          double *,double *,double *,
                          int ,int *,int *,
                          int *,
                          double *,double *,
                          double *,double *,
                          double *,
                          int *,list_int *,int *,
                          int *,list_int *,int *,
                          int ,int *);



void vercreate_list_count(double *,double *,double *,
                          int ,int *,int *,
                          int *,
                          double *,double *,
                          double *,double *,
                          double *,
                          int *,list_int *,int *,
                          int *,list_int *,int *,
                          int );


void vercreate_list_count_root(double *,double *,double *,
                               int ,int *,int *,
                               int *,
                               double *,double *,
                               double *,double *,
                               double *,
                               int *,list_int *,int *,
                               int *,list_int *,int *,
                               int ,int *);


void vercreate_list_fill(double *,double *,double *,
                         int ,int *,int *,
                         int *,
                         double *,double *,
                         double *,double *,
                         double *,
                         int *,list_int *,int *,
                         int *,list_int *,int *,
                         int ,int , int *);

void vercreate_list_fill_root(double *,double *,double *,
                              int ,int *,int *,
                              int *,
                              double *,double *,
                              double *,double *,
                              double *,
                              int *,list_int *,int *,
                              int *,list_int *,int *,
                              int ,int , int *,int *);



void expand_brnch_root_cnt(int *,int *,
                        int *,
                        int ,int *,int *,
                        int ,int *,
                        int ,int *,
                        int **,int *,
                        CLASS_COMM_FORC_PKG *,int);



void expand_brnch_root_cntfll(NBR_LIST *,int ,
                              int *,
                              int ,int *,
                              int *,int *,
                              int ,int *,
                              int *,int **,
                              int ,int *,
                              int *,int *,
                              int *,
                              int ,int ,
                              int ,double ,int *,
                              CLASS_COMM_FORC_PKG *);


void expand_brnch_root_fll(NBR_LIST *,int ,
                              int *,
                              int ,int *,
                              int *,int *,
                              int ,int *,
                              int *,int **,
                              int ,int *,
                              int *,int *,
                              int *,
                              int ,int ,int ,
                              int ,double ,
                              CLASS_COMM_FORC_PKG *);

void ver_pare_down_root(CLATOMS_INFO *,CLATOMS_POS *,
                        FOR_SCR *,ATOMMAPS *,
                        CELL *,INTERACT *,
                        TIMEINFO *, NBR_LIST *, EXCL *,
                        INTRA_SCR *,int *,int *,
                        int);

