/*==========================================================================*/

void assign_atm_pme(int ,int *,double *,double *,double *,
                    double *, double *, double *,
                    int *,list_int *,list_int *,
                    int *,int *, int ,int ,int , 
                    int ,int ,int *,int *);

void sort_pme(int , double [],int []);
void sort_pme_short(int , double [],list_int []);



void ewald3d_recip_pme_full_g(CLATOMS_INFO *,CLATOMS_POS *,
                       CELL *,PTENS *,double ,int ,
                       int [],int [],int [],
                       EWD_SCR *,double *, double , int,PART_MESH *,
                       int , CLASS_COMM_FORC_PKG *,int,PARA_FFT_PKG3D *,
                       FOR_SCR *,double *, double *);


void ewald3d_recip_pme_none(CLATOMS_INFO *,CLATOMS_POS *,
                       CELL *,PTENS *,double ,int ,
                       int [],int [],int [],
                       EWD_SCR *,double *, double , int,PART_MESH *,
                       int , CLASS_COMM_FORC_PKG *,int,PARA_FFT_PKG3D *, 
                       double *, double *);

void ewald3d_recip_pme_hybr(CLATOMS_INFO *,CLATOMS_POS *,
                       CELL *,PTENS *,double ,int ,
                       int [],int [],int [],
                       EWD_SCR *,double *, double , int,PART_MESH *,
                       int , CLASS_COMM_FORC_PKG *,int,PARA_FFT_PKG3D *, 
                       double *, double *);

/*==========================================================================*/
