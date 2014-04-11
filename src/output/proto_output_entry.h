/*----------------------------------------------------------------------*/
/* Output_md */

void output_md(CLASS *,GENERAL_DATA *,BONDED *);

/*----------------------------------------------------------------------*/
/* Output_pimd */

void output_pimd(CLASS *,GENERAL_DATA *,BONDED *);

/*----------------------------------------------------------------------*/
/* Output_min.c */

 void output_min(CLASS *,GENERAL_DATA *,BONDED *);

/*----------------------------------------------------------------------*/
/* Simpavg_md.c */

void simpavg_md(TIMEINFO *,STAT_AVG *,CELL *, CONSTRNT *,
                ENSOPTS *, SIMOPTS *, PTENS *,COMMUNICATE *,VERLIST *,
                ENERGY_CTRL *);

void simpavg_md_communicate(STAT_AVG *,ENSOPTS *,PTENS *, int,COMMUNICATE *,
                             VERLIST *,int,int ,int);

/*----------------------------------------------------------------------*/
/* Simpavg_pimd.c */

void simpavg_pimd(TIMEINFO *,STAT_AVG *,CELL *, CONSTRNT *,
                  ENSOPTS *, SIMOPTS *, PTENS *,COMMUNICATE *);

void simpavg_pimd_communicate(STAT_AVG *,ENSOPTS *,PTENS *, int,COMMUNICATE *);

/*----------------------------------------------------------------------*/

void communicate_output_pimd(CLASS *);


void forc_level_vel_gather(CLASS *);
