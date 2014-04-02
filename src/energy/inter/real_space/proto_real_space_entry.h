/*--------------------------------------------------------------------------*/
/* Force_control.c */

void force_control(CLATOMS_INFO *,CLATOMS_POS *,
                   FOR_SCR *,ATOMMAPS *,CELL *,PTENS *,INTERACT *,
                   ENERGY_CTRL *, NBR_LIST *,EXCL *, INTRA_SCR *, double *,
		   double *,double *,int,CLASS_COMM_FORC_PKG *);
/*--------------------------------------------------------------------------*/
/* Nbr_list_control.c */

void nbr_list_control(CLATOMS_INFO *,CLATOMS_POS *,
                      FOR_SCR *,ATOMMAPS *, CELL *,INTERACT *,
                      TIMEINFO *, NBR_LIST *,EXCL *,INTRA_SCR *,STAT_AVG *,
                      COMMUNICATE *,int,CLASS_COMM_FORC_PKG *);

/*--------------------------------------------------------------------------*/



