/*-----------------------------------------------------------------------*/
/* bond_con.c */

void shake_bond(BOND *,CLATOMS_INFO *,CLATOMS_POS *, 
                CELL *,INTRA_SCR *,PTENS *,double ,int );

void rattle_bond(BOND *,CLATOMS_INFO *,CLATOMS_POS *, 
                 CELL *,INTRA_SCR *,PTENS *,double ,int );

void flip_bond_con(BOND *);

/*-----------------------------------------------------------------------*/
/* bond_con_rolli.c */

void shake_bond_roll_i(BOND *,CLATOMS_INFO *,CLATOMS_POS *, 
                       CELL *,INTRA_SCR *,PTENS *,
                       BARO *,double ,int,CLASS_COMM_FORC_PKG * );

void rattle_bond_roll_i(BOND *,CLATOMS_INFO *,CLATOMS_POS *, 
                        CELL *,INTRA_SCR *,PTENS *,
                        BARO *,double ,int,CLASS_COMM_FORC_PKG * );

/*-----------------------------------------------------------------------*/
/* bond_con_fullf.c */

void shake_bond_roll_f(BOND *,CLATOMS_INFO *,CLATOMS_POS *, 
                       CELL *,INTRA_SCR *,PTENS *,
                       PAR_RAHMAN *,double ,int,CLASS_COMM_FORC_PKG * );

void rattle_bond_roll_f(BOND *bond,
                        CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                        CELL *cell,INTRA_SCR *intra_scr,PTENS *ptens,
                        PAR_RAHMAN *, double ,int,CLASS_COMM_FORC_PKG * );

/*-----------------------------------------------------------------------*/
/* Group 21 constraint routines                                          */

void shake_21(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,CLASS_COMM_FORC_PKG *);

void shake_21_rolli(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,
              BARO *,int ,CLASS_COMM_FORC_PKG *);

void shake_21_rollf(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,
              PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

void rattle_21(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,CLASS_COMM_FORC_PKG *);

void rattle_21_rolli(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,BARO *,int ,CLASS_COMM_FORC_PKG *);

void rattle_21_rollf(GRP_BOND_CON *,
              CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,
              PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

/*-----------------------------------------------------------------------*/
/* Group 23 constraint routines                                          */

void shake_23(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double ,
              double *,CLASS_COMM_FORC_PKG *);

void shake_23_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
                    double *,BARO *,int ,CLASS_COMM_FORC_PKG *);

void shake_23_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
                    double *,PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

void rattle_23(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double ,
               CLASS_COMM_FORC_PKG *);

void rattle_23_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                     double ,BARO *, int ,CLASS_COMM_FORC_PKG *);

void rattle_23_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                     double ,PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

/*-----------------------------------------------------------------------*/
/* Group 33 constraint routines                                          */

void shake_33(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double
                            ,double *,CLASS_COMM_FORC_PKG *);

void shake_33_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
              double *,BARO *,int ,CLASS_COMM_FORC_PKG *);

void shake_33_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
              double *,PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

void rattle_33(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double ,
               CLASS_COMM_FORC_PKG *);

void rattle_33_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                    double , BARO *, int , CLASS_COMM_FORC_PKG *);

void rattle_33_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                     double ,PAR_RAHMAN *,int ,CELL * , CLASS_COMM_FORC_PKG *);

/*-----------------------------------------------------------------------*/
/* Group 43 constraint routines                                          */

void shake_43(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,CLASS_COMM_FORC_PKG *);

void shake_43_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,BARO *,int ,CLASS_COMM_FORC_PKG *);

void shake_43_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,double *,
              PAR_RAHMAN *,int ,CELL *  ,CLASS_COMM_FORC_PKG *);

void rattle_43(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,CLASS_COMM_FORC_PKG *);

void rattle_43_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,BARO *,int ,CLASS_COMM_FORC_PKG *);

void rattle_43_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,
              PTENS *,double ,PAR_RAHMAN *,int ,CELL * ,
              CLASS_COMM_FORC_PKG *);

/*-----------------------------------------------------------------------*/
/* Group 46 constraint routines                                          */

void shake_46(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double ,
              double *,CLASS_COMM_FORC_PKG *);

void shake_46_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
                    double *,BARO *, int ,CLASS_COMM_FORC_PKG *);

void shake_46_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double,
                    double *,PAR_RAHMAN *,int ,CELL * ,CLASS_COMM_FORC_PKG *);

void rattle_46(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,double ,
               CLASS_COMM_FORC_PKG *);

void rattle_46_rolli(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                    double ,BARO *, int ,CLASS_COMM_FORC_PKG *);

void rattle_46_rollf(GRP_BOND_CON *,CLATOMS_INFO *,CLATOMS_POS *,PTENS *,
                     double ,PAR_RAHMAN *, int ,CELL * ,
                     CLASS_COMM_FORC_PKG *);

/*-----------------------------------------------------------------------*/
/* Constraint_control                                                    */

void recalc_pos_rollf(CLATOMS_INFO *,CLATOMS_POS *,CELL *,PTENS *,
                      PAR_RAHMAN *,double ,int,CLASS_COMM_FORC_PKG *);

void recalc_pos_rolli(CLATOMS_INFO *,CLATOMS_POS *,CELL *,PTENS *,
                      BARO *,double , int,CLASS_COMM_FORC_PKG *);

void recalc_vel_rollf(CELL *,PTENS *,PAR_RAHMAN *,double ,
                      CLASS_COMM_FORC_PKG *);

void recalc_vel_rolli(CELL *,PTENS *,BARO *,double ,CLASS_COMM_FORC_PKG *);

void all_gather_class_vec(double *,double *,double *,CLATOMS_INFO *,
                          CLASS_COMM_FORC_PKG *,double *,double *,
                          GRP_BOND_CON *);

/*--------------------------------------------------------------------------*/
/* group_bond_con_util                                                      */

void nrerror(char []);

void free_dvector(double *, int , int );

void free_dmatrix(double **, int , int , int , int );

void free_d3tensor(double ***, int , int , int , int , int , int );

double *dvector(int , int );

double **dmatrix(int , int , int , int );

double ***d3tensor(int , int , int , int , int , int );

/*--------------------------------------------------------------------------*/






