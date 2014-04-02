/*--------------------------------------------------------------------------*/
/* /pimd/int_NVT_pimd.c */

void apply_NHC_bead(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *, 
               INT_SCR *,int );

void apply_NHC_bead_par(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *, 
               INT_SCR *,int ,CLASS_COMM_FORC_PKG *);


/*--------------------------------------------------------------------------*/
/* /pimd/int_NPTI.c */

void apply_NHCPI_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                      THERM_INFO *,THERM_POS *,
                      THERM_INFO *,THERM_POS *,
                      BARO *, INT_SCR *,int , MPI_Comm);

void apply_NHCPI0_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                       THERM_INFO *,THERM_POS *,
                       THERM_INFO *,THERM_POS *,
                       BARO *, INT_SCR *,int , MPI_Comm);

/*--------------------------------------------------------------------------*/
/* /pimd/int_NPTF.c */

void apply_NHCPF_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                      THERM_INFO *,THERM_POS *,
                      THERM_INFO *,THERM_POS *,
                      BARO *,PAR_RAHMAN *,CELL *, INT_SCR *,int,MPI_Comm,int);

void apply_NHCPF0_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                       THERM_INFO *,THERM_POS *,
                       THERM_INFO *,THERM_POS *,
                       BARO *, 
                       PAR_RAHMAN *,CELL *, INT_SCR *,int,MPI_Comm,int);

/*--------------------------------------------------------------------------*/
/* /pimd/int_utils.c */

void get_tvten_pimd(CLASS *,GENERAL_DATA *);


void nhc_vol_potkin_pimd(CLASS *,GENERAL_DATA *,int);

void anneal_bead(CLATOMS_INFO *,CLATOMS_POS *,INT_SCR *,CLASS_COMM_FORC_PKG *,
                 THERM_INFO *,THERM_POS *,
                 double ,int ,int ,int ,double *);


/*--------------------------------------------------------------------------*/


void communicate_utils_pimd(double *, COMMUNICATE *);


void int_0_to_dt2_nvt_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);

void int_0_to_dt2_npti_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);

void int_0_to_dt2_nptf_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);


void int_dt2_to_dt_nvt_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);

void int_dt2_to_dt_npti_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);

void int_dt2_to_dt_nptf_pimd(CLASS *,BONDED *,GENERAL_DATA *,int,int,int,int, 
                        double);


/*--------------------------------------------------------------------------*/
/* move_pos_box.c */

 void move_pos_box(CELL *,PAR_RAHMAN *,
                         double *,double *,
                         double *,double *,
                         double *,double *,
                         double ,int);



 void move_pos_box_upper(double *,double *,
                         double *,double *,
                         double *,double *,
                         double *, double *,double *, 
                         double ,int,int, int );



