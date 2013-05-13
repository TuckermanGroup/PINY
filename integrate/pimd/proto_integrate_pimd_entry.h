/*--------------------------------------------------------------------------*/
/* /pimd/int_NVT_res.c */

void int_NVT_pimd(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

/*--------------------------------------------------------------------------*/
/* /pimd/int_NPTI.c */

void int_NPTI_pimd(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

void init_NHCPI_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                     THERM_INFO *, THERM_POS *,
                     THERM_INFO *, THERM_POS *,
                     BARO *, INT_SCR *,int , MPI_Comm );

/*--------------------------------------------------------------------------*/
/* /pimd/int_NPTF.c */

void int_NPTF_pimd(CLASS *,BONDED *,GENERAL_DATA *,ANALYSIS *);

void init_NHCPF_pimd(CLATOMS_INFO *,CLATOMS_POS *,
                     THERM_INFO *,THERM_POS *,
                     THERM_INFO *,THERM_POS *,
                     BARO *, 
                     PAR_RAHMAN *, CELL *, INT_SCR *,
                     int , MPI_Comm);

void get_fnhc1_pimd(CLATOMS_INFO *,CLATOMS_POS *, 
                    THERM_INFO *,THERM_POS *,
                    THERM_INFO *,THERM_POS *,
                    BARO *, PAR_RAHMAN *,int);

void get_fgmatv_pimd(CLATOMS_INFO *,CLATOMS_POS *,
                PAR_RAHMAN *, CELL *,int ,MPI_Comm);

/*--------------------------------------------------------------------------*/
