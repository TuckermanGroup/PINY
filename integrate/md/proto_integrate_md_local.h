/*--------------------------------------------------------------------------*/
/* /md/int_NVT.c */

void apply_NHC(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *,
               INT_SCR *, int);

void apply_NHC_par(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *,
               INT_SCR *, int,CLASS_COMM_FORC_PKG *);

void apply_GGMT2_par(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *, 
               INT_SCR *,int ,CLASS_COMM_FORC_PKG *);

void apply_GGMT3_par(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *, 
               INT_SCR *,int ,CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* /md/int_NPTI.c */

void apply_NHCPI(CLATOMS_INFO *,CLATOMS_POS *, THERM_INFO *,
                 THERM_POS *,BARO *, INT_SCR *);

void apply_NHCPI0(CLATOMS_INFO *,CLATOMS_POS *, THERM_INFO *,
                  THERM_POS *,BARO *, INT_SCR *);

void apply_NHCPI0_par(CLATOMS_INFO *,CLATOMS_POS *, THERM_INFO *,
                  THERM_POS *,BARO *, INT_SCR *,CLASS_COMM_FORC_PKG *);

void apply_NHCPI_par(CLATOMS_INFO *,CLATOMS_POS *, 
                 THERM_INFO *,THERM_POS *, BARO *, INT_SCR *,
                 CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* /md/int_NPTF.c */

void apply_NHCPF(CLATOMS_INFO *,CLATOMS_POS *, 
                 THERM_INFO *,THERM_POS *,
                 BARO *,PAR_RAHMAN *,CELL *, INT_SCR *,double);

void apply_NHCPF_par(CLATOMS_INFO *,CLATOMS_POS *, 
                 THERM_INFO *,THERM_POS *,
                 BARO *,PAR_RAHMAN *,CELL *, INT_SCR *,CLASS_COMM_FORC_PKG *);
                 

void apply_NHCPF0(CLATOMS_INFO *,CLATOMS_POS *, THERM_INFO *,
                  THERM_POS *,BARO *, 
                  PAR_RAHMAN *,CELL *, INT_SCR *,double );

void apply_NHCPF0_par(CLATOMS_INFO *,CLATOMS_POS *, THERM_INFO *,
                  THERM_POS *,BARO *, 
                  PAR_RAHMAN *,CELL *, INT_SCR *,CLASS_COMM_FORC_PKG *);

void move_vel_vbox0(CLATOMS_POS *, CLATOMS_INFO *, CELL *,PAR_RAHMAN *,
                    double *,double *,double *,double *,double *,int,int );

void move_vel_vbox_upper0(CLATOMS_POS *, CLATOMS_INFO *, CELL *,PAR_RAHMAN *,
                    double *,double *,double *,double *,double *,int,int,
                    double );

/*--------------------------------------------------------------------------*/
/* /md/int_utils.c */

void get_tvten(CLATOMS_INFO *, CLATOMS_POS *,STAT_AVG *,PTENS *, CELL *);


void get_pvten_inc_t(PTENS *, TIMEINFO *,int , int , int );


void get_pvten_inc_a(PTENS *, TIMEINFO *);


void set_yosh(int ,double ,double [],double [],double [],double [],double[]);

void cpysys_NPT(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *,
                BARO *,PAR_RAHMAN *,INT_SCR *,int );

void getsys_NPT(CLATOMS_INFO *,CLATOMS_POS *,THERM_INFO *,THERM_POS *,
                BARO *,PAR_RAHMAN *,INT_SCR *,int );

void nhc_vol_potkin(THERM_INFO *, THERM_POS *,BARO *, PAR_RAHMAN *,
                    STAT_AVG *,STATEPOINT *, int,int );

/*--------------------------------------------------------------------------*/

void int_0_to_dt2_nve(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_0_to_dt2_nvt(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_0_to_dt2_npti(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_0_to_dt2_nptf(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_dt2_to_dt_nve(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_dt2_to_dt_nvt(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_dt2_to_dt_npti(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_dt2_to_dt_nptf(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,double);

void int_final_class(CLASS *,BONDED *,GENERAL_DATA *,int );

void anneal_class(CLASS *,double,int,int,double,int *);

void class_int_allgather(GENERAL_DATA *,CLASS *);



/*--------------------------------------------------------------------------*/
/* move_vel_vbox.c */

void move_vel_vbox(CLATOMS_POS *, CLATOMS_INFO *,
                   CELL *,PAR_RAHMAN *, double *,double *,
                   double *,double *, double *,
                   double **,double **,
                   int *,int *,int *,int , 
                   int , double *, double *,double , int, int, int );

void move_vel_vbox_upper(CLATOMS_POS *, CLATOMS_INFO *,
                   CELL *,PAR_RAHMAN *, double *,double *,
                   double *,double *, double *,
                   double **,double **,
                   int *,int *,int *,int , 
                   int , double *, double *,double , int);




