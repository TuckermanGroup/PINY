/*-------------------------------------------------------------*/
/* Bond.c */

void bond(CLATOMS_INFO *,CLATOMS_POS *,
          BOND *,CELL *,INTRA_SCR *,PTENS *,double *, int, int,
          CLASS_COMM_FORC_PKG *,int );

void bond_both(CLATOMS_INFO *,CLATOMS_POS *,
          BOND *,CELL *,
          INTRA_SCR *,PTENS *, double *, int ,
          CLASS_COMM_FORC_PKG *,int );

void bond_watts_33(CLATOMS_INFO *,CLATOMS_POS *,
          GRP_BOND_WATTS *,CELL *,
          INTRA_SCR *,PTENS *, double *,
          double *, double *, int,
          CLASS_COMM_FORC_PKG *,int);

void bond_free(CLATOMS_INFO *,CLATOMS_POS *,
               BOND_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,
               ENERGY_CTRL *,int );

void bond_free_mode(CLATOMS_INFO *,CLATOMS_POS *,
               BOND_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,
               ENERGY_CTRL *,int );

/*-------------------------------------------------------------*/
/* Bend.c */

void bend(CLATOMS_INFO *,CLATOMS_POS *,
          BEND *,CELL *,
          INTRA_SCR *,PTENS *,double *, int ,
          CLASS_COMM_FORC_PKG *,int );

void bend_free(CLATOMS_INFO *,CLATOMS_POS *,
               BEND_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,
               ENERGY_CTRL *,int);

void bend_free_mode(CLATOMS_INFO *,CLATOMS_POS *,
          BEND_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,ENERGY_CTRL *,
          int );

/*-------------------------------------------------------------*/
/* Bend_bnd.c */

void bend_bnd(CLATOMS_INFO *,CLATOMS_POS *,
              BEND_BND *,CELL *,INTRA_SCR *,PTENS *,double *,
              double *,double *,int,
              CLASS_COMM_FORC_PKG *,int );

/*-------------------------------------------------------------*/
/* Tors.c */

void tors(CLATOMS_INFO *,CLATOMS_POS *,
          TORS *,CELL *,INTRA_SCR *,PTENS *,double *,int,
          CLASS_COMM_FORC_PKG *,int );

void tors_free(CLATOMS_INFO *,CLATOMS_POS *,
               TORS_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,
               ENERGY_CTRL *,int);

void tors_free_mode(CLATOMS_INFO *,CLATOMS_POS *,
   TORS_FREE *,CELL *,INTRA_SCR *,PTENS *,double *,ENERGY_CTRL *,int );

/*-------------------------------------------------------------*/
/* onefour.c */

void onfo(CLATOMS_INFO *,CLATOMS_POS *,
         ONFO *,CELL *,INTRA_SCR *,PTENS *, double *, double *, double *,int,
         CLASS_COMM_FORC_PKG *,int );

/*-------------------------------------------------------------*/
/* ecorr.c */

void ecor(CLATOMS_INFO *,CLATOMS_POS *,
          ECOR *,CELL *,INTRA_SCR *,PTENS *, double *,int,
          CLASS_COMM_FORC_PKG *,int ,int);

void ecor_both(CLATOMS_INFO *,CLATOMS_POS *,
               ECOR *,CELL *,INTRA_SCR *,PTENS *, double *, int ,
               CLASS_COMM_FORC_PKG *,int ,int );

void ecor_vspl_fetch(int ,double [], double [],int [],double []);
void ecor_pvten_sum_roll(int ,double *, double *, double *,
                         double *,double *,double *,
                         double *,double *,double *,double *,
                         double *,double *,double *,int );

/*-------------------------------------------------------------*/
/* rbar_sigma.c */

void rbar_sig_free(CLATOMS_INFO *,CLATOMS_POS *,RBAR_SIG_FREE *,CELL *,
                   INTRA_SCR *,PTENS *, double *, ENERGY_CTRL *,int );

/*-------------------------------------------------------------*/
