/*========================================================================*/

void mall_integrator_scr(int ,int ,CLATOMS_INFO *,THERM_INFO *,
                        THERM_INFO *,INT_SCR *,double *,int ,MPI_Comm );

void mall_intra_scr(INTRA_SCR *,double *, int ,MPI_Comm );

void mall_ewald_scr(int ,int ,int ,int ,int,PART_MESH *,EWD_SCR *, 
                    double *,int ,MPI_Comm );

void mall_cp_scr(CPTHERM_INFO *,CPOPTS *,CPEWALD *,
                 CPSCR *,CPCOEFFS_INFO *,PSEUDO *,
                 PARA_FFT_PKG3D *,PARA_FFT_PKG3D *,
                 CP_COMM_STATE_PKG *,CP_COMM_STATE_PKG *,int,double *,
                 int ,int ,MPI_Comm);

void mall_atm_forc_scr(int ,FOR_SCR *,int ,int ,int , int ,
                       double *,int ,int, MPI_Comm);

/*==========================================================================*/
