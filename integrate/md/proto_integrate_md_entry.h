/*--------------------------------------------------------------------------*/
/* /md/int_NVE.c */

void int_NVE(CLASS *,BONDED *,GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* /md/int_NVE_res.c */

void int_NVE_res(CLASS *,BONDED *,GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* /md/int_NVT.c */

void int_NVT(CLASS *,BONDED *,GENERAL_DATA *);

void init_NHC(CLATOMS_INFO *,CLATOMS_POS *, 
              THERM_INFO *, THERM_POS *,
              INT_SCR *, int);

void init_NHC_par(CLATOMS_INFO *,CLATOMS_POS *, 
              THERM_INFO *, THERM_POS *,
              INT_SCR *, int, CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* /md/int_NVT_res.c */

void int_NVT_res(CLASS *,BONDED *,GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* /md/int_NPTI.c */

void int_NPTI(CLASS *,BONDED *,GENERAL_DATA *);

void init_NHCPI(CLATOMS_INFO *,CLATOMS_POS *, 
                THERM_INFO *, THERM_POS *,
                BARO *, INT_SCR *);

void init_NHCPI_par(CLATOMS_INFO *,CLATOMS_POS *, 
                THERM_INFO *, THERM_POS *,
                BARO *, INT_SCR *,CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* /md/int_NVPI_res.c */

void int_NPTI_res(CLASS *,BONDED *,GENERAL_DATA *);

/*--------------------------------------------------------------------------*/
/* /md/int_NPTF.c */

void int_NPTF(CLASS *,BONDED *,GENERAL_DATA *);

void init_NHCPF(CLATOMS_INFO *,CLATOMS_POS *, 
                THERM_INFO *,THERM_POS *,
                BARO *, PAR_RAHMAN *,CELL *, INT_SCR *);

void init_NHCPF_par(CLATOMS_INFO *,CLATOMS_POS *, 
                THERM_INFO *,THERM_POS *,
                BARO *, PAR_RAHMAN *,CELL *, INT_SCR *,CLASS_COMM_FORC_PKG *);

void get_fnhc1(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *,
               BARO *, PAR_RAHMAN *);

void get_fnhc1_par(CLATOMS_INFO *,CLATOMS_POS *, 
               THERM_INFO *,THERM_POS *,
               BARO *, PAR_RAHMAN *,CLASS_COMM_FORC_PKG *,
               INT_SCR *);

void get_fgmatv(CLATOMS_INFO *,CLATOMS_POS *, PAR_RAHMAN *, CELL *);

void get_fgmatv_par(CLATOMS_INFO *,CLATOMS_POS *, PAR_RAHMAN *, CELL *,
                    CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/
/* /md/int_NVPF_res.c */

void int_NPTF_res(CLASS *,BONDED *,GENERAL_DATA *);


