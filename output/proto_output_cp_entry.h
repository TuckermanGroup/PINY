/*----------------------------------------------------------------------*/
/* Output_cp */

void output_cp(CLASS *,GENERAL_DATA *,BONDED *,CP *);

/*----------------------------------------------------------------------*/
/* Output_cp_min */

void output_cp_min(CLASS *,GENERAL_DATA *,BONDED *,CP *,double ,int);

/*----------------------------------------------------------------------*/
/* Output_cp_pimd */

void output_cp_pimd(CLASS *,GENERAL_DATA *,BONDED *,CP *);

/*----------------------------------------------------------------------*/
/* Simpavg_cp.c */

void simpavg_cp(TIMEINFO *,STAT_AVG *,CELL *, CONSTRNT *,
                ENSOPTS *, SIMOPTS *, PTENS *, CP *cp,
                COMMUNICATE *,VERLIST *,ENERGY_CTRL *);

/*----------------------------------------------------------------------*/
/* Simpavg_cp_pimd.c */

void simpavg_cp_pimd(TIMEINFO *,STAT_AVG *,CELL *, CONSTRNT *,
                ENSOPTS *, SIMOPTS *, PTENS *, CP *,COMMUNICATE *);

/*----------------------------------------------------------------------*/


