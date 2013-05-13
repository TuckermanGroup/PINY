/*-----------------------------------------------------------------------*/
/* ewald3d_recip.c */

void ewald3d_recip(CLATOMS_INFO *,CLATOMS_POS *,
                   CELL *,PTENS *,double ,double ,int ,
                   int [],int [],int [], int [],int [],
                   EWD_SCR *,double *, double , int,
                   CLASS_COMM_FORC_PKG *,int,
                   double *, double *);

/*-----------------------------------------------------------------------*/
/* ewald3d_recip_both.c */

void ewald3d_recip_both(CLATOMS_INFO *,CLATOMS_POS *,
                     CELL *,PTENS *,
                     double ,int , int [],int [],int [],
                     int [],int [],int [],EWD_SCR *,
                     double *,double , double , int,CLASS_COMM_FORC_PKG *,int, 
                     double *, double *);

/*-----------------------------------------------------------------------*/
/* ewald3d_self_bgr.c */

void ewald3d_selfbgr(CLATOMS_INFO *,EWALD *,PTENS *,double ,double ,
                     double *, double *,int,int,int);

/*-----------------------------------------------------------------------*/
/* ewald3d_recip_pme */

void ewald3d_recip_pme(CLATOMS_INFO *,CLATOMS_POS *,
                       CELL *,PTENS *,double ,int ,
                       int [],int [],int [],
                       EWD_SCR *,double *, double , int,PART_MESH *,
                       int , CLASS_COMM_FORC_PKG *,int,PARA_FFT_PKG3D *,
                       FOR_SCR *,double *, double *);

/*-----------------------------------------------------------------------*/







