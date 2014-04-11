/*-------------------------------------------------------------------------*/
/* control_intra_params.c */

void control_intra_params(double *,CLATOMS_INFO *,CLATOMS_POS *,
                          GHOST_ATOMS *,
                          ATOMMAPS *,BONDED *,
			  FILENAME_PARSE *,
			  FREE_PARSE *,CLASS_PARSE *,NULL_INTER_PARSE *,
                          SIMOPTS *,COMMUNICATE *, int );

/*-------------------------------------------------------------------------*/
/* close_intra_params.c */

void close_intra_params(CLATOMS_INFO *,CLATOMS_POS *,GHOST_ATOMS *,
                   ATOMMAPS *,
                   BUILD_INTRA *,
                   BONDED *,
                   NULL_INTER_PARSE *,double *, SIMOPTS *,int);

