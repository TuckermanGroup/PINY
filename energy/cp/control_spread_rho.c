/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_spread_rho.c                         */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define REAL_CASE
#define DEBUG_PME_PAR1_OFF
#define ORIG_OFF
#define PME

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*Controller: If parallel create the maps for the packing and upacking of   */
/*            the density from the small grid to the large grid             */
/*  COMMENSUERATE GRID OPTION                                               */
/* This routine is called from parse                                        */
/*==========================================================================*/
void control_init_dual_maps(CPSCR *cpscr, CELL *cell,double dbox_rat,
                            PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                            PARA_FFT_PKG3D *para_fft_pkg3d_lg,
                            int cp_dual_grid_opt)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/* local variables  */  

  int iset_map_flag,iset_map_upack_flag;
  int  nfft2_up_proc_dens_cp_box;
  int i;

/* local pointers  */  
  double *rfft         = cpscr->cpscr_rho.rho_up;
  double *zfft         = cpscr->cpscr_wave.zfft;
  double *zfft_tmp     = cpscr->cpscr_wave.zfft_tmp;
  int *map_dual        = cpscr->cpscr_rho.map_dual;
  int *map_dual_upack  = cpscr->cpscr_rho.map_dual_upack; 

  int  nfft_up_proc_dens_cp_box  = para_fft_pkg3d_cp_box->nfft_proc; 

/*--------------------------------------------------------------------------*/
/* Zero the rfft array  */
  nfft2_up_proc_dens_cp_box = nfft_up_proc_dens_cp_box/2;

 for(i=1; i<= nfft2_up_proc_dens_cp_box;i++){
    rfft[i] = 0.0;
 }/*endfor */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*I) cp_dual_grid_opt = 1  commensuerate grids      */
/* make the map: small cp_box_grid  to large grid   */

 if(cp_dual_grid_opt == 1){
  cpscr->cpscr_rho.iset_map_flag = 0;
  iset_map_flag = cpscr->cpscr_rho.iset_map_flag;

  sngl_pack_rho_dual_par(zfft,zfft_tmp,rfft,
                         map_dual,iset_map_flag,
                         cell,dbox_rat, 
                         para_fft_pkg3d_cp_box,
                         para_fft_pkg3d_lg);

   cpscr->cpscr_rho.iset_map_flag = 1;
 }/*endif*/


/*--------------------------------------------------------------------------*/
/* make the unpack map: large grid to small cp box grid   */

 if(cp_dual_grid_opt == 1){
  cpscr->cpscr_rho.iset_map_upack_flag = 0;
  iset_map_upack_flag = cpscr->cpscr_rho.iset_map_upack_flag;

   sngl_upack_rho_dual_par(zfft,zfft_tmp,rfft,
                           map_dual_upack,
                           iset_map_upack_flag,cell,
                           para_fft_pkg3d_cp_box,
                           para_fft_pkg3d_lg);   

  cpscr->cpscr_rho.iset_map_upack_flag = 1;
 }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*Controller: If PME create the maps for the packing and upacking of   */
/*  the density from the small dense grid to the large sparse grid     */
/*  INCOMMENSUERATE GRID OPTION                                        */
/* This routine is called from parse                                   */
/*==========================================================================*/
void control_init_dual_maps_pme(double dbox_rat,int n_interp_pme_dual,
                                int cp_dual_grid_opt,
                                PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                                PARA_FFT_PKG3D *para_fft_pkg3d_lg,
                                CPSCR *cpscr, CELL *cell)
                            
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/* local variables  */  

  int iset_map_flag;
  int i;

/* local pointers  */  
  double *rfft         = cpscr->cpscr_rho.rho_up;
  double *zfft         = cpscr->cpscr_wave.zfft;

  int  num_proc  = para_fft_pkg3d_lg->num_proc;
  int  nfft_up_proc_dens_cp_box = para_fft_pkg3d_cp_box->nfft_proc;
  int  nfft2_up_proc_dens_cp_box;

/*--------------------------------------------------------------------------*/
/* Zero the rfft array  */

  nfft2_up_proc_dens_cp_box = nfft_up_proc_dens_cp_box/2;

 for(i=1; i<= nfft2_up_proc_dens_cp_box;i++){
    rfft[i] = 0.0;
 }/*endfor */

/*--------------------------------------------------------------------------*/
/*I) cp_dual_grid_opt = 2  incommensuerate grids  */

 if(cp_dual_grid_opt == 2){
  if(num_proc == 1 ){  /*SERIAL create PME map */
  cpscr->cpscr_rho.iset_map_flag = 0;
  iset_map_flag = cpscr->cpscr_rho.iset_map_flag;


#ifdef DEBUG_PME_PAR1
printf("checking pme_par_Upack routine in serial \n");
    sngl_pack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,cell,dbox_rat,
                               n_interp_pme_dual,iset_map_flag,
                               &(cpscr->cpscr_dual_pme),para_fft_pkg3d_cp_box,
                               para_fft_pkg3d_lg);

#endif


#ifdef  REAL_CASE

    sngl_pack_rho_dual_ser_pme(zfft,rfft,cell,dbox_rat,
                               n_interp_pme_dual,
                               iset_map_flag,cpscr,
                               para_fft_pkg3d_cp_box,
                               para_fft_pkg3d_lg);  
#endif

  cpscr->cpscr_rho.iset_map_flag = 1;

  }else{ /*PARALLEL create PME map */

  cpscr->cpscr_rho.iset_map_flag = 0;
  iset_map_flag = cpscr->cpscr_rho.iset_map_flag;

    sngl_pack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,cell,dbox_rat,
                               n_interp_pme_dual,iset_map_flag,
                               &(cpscr->cpscr_dual_pme),para_fft_pkg3d_cp_box,
                               para_fft_pkg3d_lg);

   cpscr->cpscr_rho.iset_map_flag = 1;
  }
 }/*endif cp_dual_grid_opt*/
/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*Controller: call correct spread routine for rho when use dual grid opt    */
/*==========================================================================*/

void control_spread_rho(CPSCR *cpscr,double *rfft,
                        CELL *cell,double dbox_rat,int np_states,
                        int n_interp_pme_dual,
                        PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                        PARA_FFT_PKG3D *para_fft_pkg3d_lg,
                        int cp_dual_grid_opt)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/* local pointers */
/* rfft contains rho*vol_cp  */                       
 
 double *zfft           =  cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =  cpscr->cpscr_wave.zfft_tmp;
 int    *map_dual       =  cpscr->cpscr_rho.map_dual;
 double cpu1,cpu2; 

 if(np_states == 1){
  switch(cp_dual_grid_opt){
   case 1:   
      sngl_pack_rho_dual_ser(zfft,rfft,cell,dbox_rat,
                             para_fft_pkg3d_cp_box,
                             para_fft_pkg3d_lg);  
   break;
   case 2:

#ifdef  ORIG
      sngl_pack_rho_dual_ser(zfft,rfft,cell,dbox_rat,
                             para_fft_pkg3d_cp_box,
                             para_fft_pkg3d_lg);  
#endif

#ifdef  PME
  cputime(&cpu1);

#ifdef DEBUG_PME_PAR1

   sngl_pack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,
                              cell,dbox_rat,
                              n_interp_pme_dual,
                              (cpscr->cpscr_rho.iset_map_flag),
                              &(cpscr->cpscr_dual_pme),
                              para_fft_pkg3d_cp_box,
                              para_fft_pkg3d_lg);

#endif

#ifdef  REAL_CASE

      sngl_pack_rho_dual_ser_pme(zfft,rfft,cell,dbox_rat,
                                 n_interp_pme_dual,
                                (cpscr->cpscr_rho.iset_map_flag),
                                 cpscr,
                                 para_fft_pkg3d_cp_box,
                                 para_fft_pkg3d_lg);  

#endif
  cputime(&cpu2);

#ifdef TIME_ME
  printf("cpu time spread %g \n",cpu2-cpu1);
#endif

#endif
   break;
  }

 }else{
  switch(cp_dual_grid_opt){
   case 1:
    sngl_pack_rho_dual_par(zfft,zfft_tmp,rfft,map_dual,
                           (cpscr->cpscr_rho.iset_map_flag),
                           cell,dbox_rat,
                           para_fft_pkg3d_cp_box,
                           para_fft_pkg3d_lg);
   break;
   case 2:
   sngl_pack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,
                              cell,dbox_rat,
                              n_interp_pme_dual,
                              (cpscr->cpscr_rho.iset_map_flag),
                              &(cpscr->cpscr_dual_pme),
                              para_fft_pkg3d_cp_box,
                              para_fft_pkg3d_lg);
   break;
  }/*end switch*/
 }/*endif np_states*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_pack_rho_dual_ser(double *zfft,double *rfft,
                            CELL *cell,double dbox_rat,
                            PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                            PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
 int i,iii;
 int nfft = para_fft_pkg3d_lg->nfft;   /* on large grid*/
 int ndata      = nfft;

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 double *hmati             = cell->hmati;
 double *cp_box_center     = cell->cp_box_center;

 int ia,ib,ic;
 int iap,ibp,icp;
 int ia_shift,ib_shift,ic_shift;

 double dia_shift,dib_shift,dic_shift;
 int index_lg_cmplx,index_sm,index_lg;

 double scale;
 double sx,sy,sz;
 double dnkf1,dnkf2,dnkf3;
 double dia,dib,dic;

/*=======================================================================*/
/* Determine useful constants                                            */

  dnkf1 = (double) nkf1;
  dnkf2 = (double) nkf2;
  dnkf3 = (double) nkf3;

   scale = dbox_rat*dbox_rat*dbox_rat;

   sx = cp_box_center[1]*hmati[1]+cp_box_center[2]*hmati[4] 
      + cp_box_center[3]*hmati[7];

   sy = cp_box_center[1]*hmati[2]+cp_box_center[2]*hmati[5]
      + cp_box_center[3]*hmati[8];

   sz = cp_box_center[1]*hmati[3]+cp_box_center[2]*hmati[6]
      + cp_box_center[3]*hmati[9];

   dia_shift = ( -sx *(double)nkf1 - 0.5*(double)nkf1_cp_box);
   dib_shift = ( -sy *(double)nkf2 - 0.5*(double)nkf2_cp_box);
   dic_shift = ( -sz *(double)nkf3 - 0.5*(double)nkf3_cp_box);

   dia_shift /= dnkf1;
   dib_shift /= dnkf2;
   dic_shift /= dnkf3;

   dia_shift -= NINT((dia_shift - .5));
   dib_shift -= NINT((dib_shift - .5));
   dic_shift -= NINT((dic_shift - .5)); 

   ia_shift = (int)(dia_shift*dnkf1);
   ib_shift = (int)(dib_shift*dnkf2);
   ic_shift = (int)(dic_shift*dnkf3);

/*=======================================================================*/
/*=======================================================================*/
/* 1.) Zero the zfft array on the large grid                             */

 for(i=1; i<=ndata; i++){
   zfft[i] = 0.0;
 }/*endfor*/

/*=======================================================================*/
/* RFFT defined on real space grid corresponding to smaller CP_BOX       */

 for(ic=1; ic <= nkf3_cp_box; ic++){
  for(ib=1; ib <= nkf2_cp_box; ib++){
   for(ia=1; ia <= nkf1_cp_box; ia++){
      index_sm = (ia) + (ib - 1)*nkf1_cp_box
               + (ic - 1)*nkf1_cp_box*nkf2_cp_box;

      iap = ia + ia_shift;
      ibp = ib + ib_shift;
      icp = ic + ic_shift;

      iap = (iap > 0 ? iap : nkf1 + iap);
      ibp = (ibp > 0 ? ibp : nkf2 + ibp);
      icp = (icp > 0 ? icp : nkf3 + icp);

      iap = (iap > nkf1 ? iap - nkf1 : iap);
      ibp = (ibp > nkf2 ? ibp - nkf2 : ibp);
      icp = (icp > nkf3 ? icp - nkf3 : icp);

   index_lg = (iap) + (ibp - 1)*nkf1
               + (icp - 1)*nkf1*nkf2; 

   index_lg_cmplx = 2*index_lg - 1; 

   zfft[index_lg_cmplx] = rfft[index_sm]*scale;  

/* rfft =  rho*vol_cp, in the next section of cp_rho_calc,*/
/* rho is divided by vol. Sinc vol = (box_rat)^3*vol_cp,  */
/* you need to multiply rfft by scale = (box_rat)^3 here  */

   }/*endfor ia*/
  }/*endfor ib*/
 }/*endfor ic*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_pack_rho_dual_ser_pme(double *zfft,double *rfft,
                                CELL *cell,double dbox_rat,
                                int n_interp_pme_dual,
                                int iset_map_flag,CPSCR *cpscr,
                                PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                                PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
 int i,iii;

 int ja,jb,jc; 
 int ka,kb,kc;

 int ioff,joff;

 double scale;
 double trap_corr;
 double temp;

 double di_nkf1_cp,di_nkf2_cp,di_nkf3_cp;
 double dnkf1,dnkf2,dnkf3;

/* Local pointers */
 int  nfft = para_fft_pkg3d_lg->nfft;         /* on large sparse grid*/
 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int **igrid_a = cpscr->cpscr_dual_pme.igrid_a; /*length nkf1_cp_box*/
 int **igrid_b = cpscr->cpscr_dual_pme.igrid_b;
 int **igrid_c = cpscr->cpscr_dual_pme.igrid_c;

 double **mn_a  = cpscr->cpscr_dual_pme.mn_a; 
 double **mn_b  = cpscr->cpscr_dual_pme.mn_b;
 double **mn_c  = cpscr->cpscr_dual_pme.mn_c;

/*=======================================================================*/
/* 1.) Zero the zfft array on the large grid                             */

 for(i=1; i<=nfft; i++){
   zfft[i] = 0.0;
 }/*endfor*/

/*=======================================================================*/
/* Determine useful constants                                            */

   scale    = dbox_rat*dbox_rat*dbox_rat;

   di_nkf1_cp = 1.0/((double) (nkf1_cp_box));
   di_nkf2_cp = 1.0/((double) (nkf2_cp_box));
   di_nkf3_cp = 1.0/((double) (nkf3_cp_box));

   dnkf1 = (double) nkf1;  /*large sparse grid*/
   dnkf2 = (double) nkf2;  /*large sparse grid*/
   dnkf3 = (double) nkf3;  /*large sparse grid*/

   trap_corr = (dnkf1*dnkf2*dnkf3*di_nkf1_cp*di_nkf2_cp*di_nkf3_cp)/
                (dbox_rat*dbox_rat*dbox_rat);

/*--------------------------------------------------------------------------*/ 
/* make PME map only needs to be done once if hmat_cp does not move         */

 if(iset_map_flag == 0){ /*create maps */
   make_dual_pme_wghts(&(cpscr->cpscr_dual_pme),cell,dbox_rat,
                       n_interp_pme_dual,
                       para_fft_pkg3d_cp_box,para_fft_pkg3d_lg);
 }/*endif set maps*/

/*--------------------------------------------------------------------------*/ 
/*--------------------------------------------------------------------------*/ 
/* C) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 


     for(jc=1;jc<=n_interp_pme_dual;jc++){
       for(jb=1;jb<=n_interp_pme_dual;jb++){

     for(kc=1; kc <= nkf3_cp_box ; kc++){
       for(kb=1; kb <= nkf2_cp_box ; kb++){

         ioff = (kb - 1)*nkf1_cp_box + (kc-1)*nkf1_cp_box*nkf2_cp_box;
         joff = igrid_b[jb][kb]+igrid_c[jc][kc];
         temp = mn_b[jb][kb]*mn_c[jc][kc]*scale*trap_corr;

/* rfft =  rho*vol_cp, in the next section of cp_rho_calc,*/
/* rho is divided by vol. Sinc vol = (box_rat)^3*vol_cp,  */
/* you need to multiply rfft by scale = (box_rat)^3 here  */
/* trap_corr is needed because we are changing the grid spacing */
/*  of the density and we still want the integral of the density */
/* to generate the number of electrons   */

         for(ja=1;ja<=n_interp_pme_dual;ja++){
         for(ka=1; ka <= nkf1_cp_box ; ka++){
          zfft[(igrid_a[ja][ka] + joff)] += 
                (mn_a[ja][ka]*temp*rfft[(ka+ioff)]);

         }}/* endfor : ka,ja */

       }}/*endfor :  kb,jb */
     }}/*endfor : kc,jc*/


/*==========================================================================*/
/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_pack_rho_dual_par(double *zfft,double *zfft_tmp,double *rfft,
                           int *map_dual,int iset_map_flag,
                           CELL *cell,double dbox_rat, 
                           PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                           PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int i;
 int nfft_proc  = para_fft_pkg3d_lg->nfft_proc;   /* on large grid*/

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 double *hmati             = cell->hmati;
 double *cp_box_center     = cell->cp_box_center;

 int ia_shift,ib_shift,ic_shift;
 double dia_shift,dib_shift,dic_shift;
 int index_lg_cmplx,index_sm,index_lg;

 double scale;

 double sx,sy,sz;
 
/*-----------------------------------------------------------------------*/
/* local pointers from fft packages                                      */

/* Large Box                                                             */

 int *skb_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->skb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->skc_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->ekb_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->ekc_fft_ka_proc_all;


/* Small CP Box                                                          */

 int   skc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekb_fft_ka_proc;

 int num_proc = para_fft_pkg3d_lg->num_proc;

 int iproc;
 int kc_sm,kb_sm,ka_sm;
 int kc_lg,kb_lg,ka_lg;
 int kb_off_sm,kc_off_sm;
 int myid             = para_fft_pkg3d_lg->myid;
 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 
 MPI_Comm comm_cp_box = para_fft_pkg3d_cp_box->comm; 
 int np_states        = para_fft_pkg3d_lg->num_proc;

 int *send_counts_row_big_small = 
                         para_fft_pkg3d_cp_box->send_counts_row_big_small;
 int *recv_counts_row_big_small =
                        para_fft_pkg3d_cp_box->recv_counts_row_big_small;

 int *sdispls_row_big_small  = para_fft_pkg3d_cp_box->sdispls_row_big_small;
 int *rdispls_row_big_small  = para_fft_pkg3d_cp_box->rdispls_row_big_small;

 int *send_counts_dual_map       = para_fft_pkg3d_cp_box->send_counts_dual_map;
 int *recv_counts_dual_map       = para_fft_pkg3d_cp_box->recv_counts_dual_map;
 int *sdispls_dual_map           = para_fft_pkg3d_cp_box->sdispls_dual_map;
 int *rdispls_dual_map           = para_fft_pkg3d_cp_box->rdispls_dual_map;


 int isend_tot,isend_tot_map,icount,icount_tot;    
 int kb_str_sm,kb_end_sm;
 int kb_str_lg,kb_end_lg;
 int kb_off,kc_off;
 int one= 1;
 int master = 0;

 int sum_count,sum_count_temp;
 int nfft_up_cp_box      = para_fft_pkg3d_cp_box->nfft;
 int nfft2_up_cp_box;
 int map_count,map_count_proc;
 int *ind_tmp;
 int ioff,ioff_map;
 int *pcount,*pcount_map;
 int ipack_on = 0;
 double dnkf1,dnkf2,dnkf3;

   dnkf1 = (double) nkf1;
   dnkf2 = (double) nkf2;
   dnkf3 = (double) nkf3;
/*=======================================================================*/
/* A local malloc that will be removed in a short while                  */

   nfft2_up_cp_box =  nfft_up_cp_box/2;
   map_count = nfft2_up_cp_box/nkf1_cp_box;   

  if(iset_map_flag == 0){
   ind_tmp = (int *) cmalloc(map_count*sizeof(int))-1;
  }/*endif*/
 
   pcount     = (int *) cmalloc(num_proc*sizeof(int))-1;
   pcount_map = (int *) cmalloc(num_proc*sizeof(int))-1;

/*=======================================================================*/
/* zero the array to hold the real space on the large grid               */

   for(i=1; i<= nfft_proc; i++){
     zfft[i]     = 0.0;
     zfft_tmp[i] = 0.0; 
   }

/*=======================================================================*/
/* Determine useful constants                                            */

   scale = dbox_rat*dbox_rat*dbox_rat;

   sx = cp_box_center[1]*hmati[1]+cp_box_center[2]*hmati[4] 
      + cp_box_center[3]*hmati[7];

   sy = cp_box_center[1]*hmati[2]+cp_box_center[2]*hmati[5]
      + cp_box_center[3]*hmati[8];

   sz = cp_box_center[1]*hmati[3]+cp_box_center[2]*hmati[6]
      + cp_box_center[3]*hmati[9];


   dia_shift = ( -sx *(double)nkf1 - 0.5*(double)nkf1_cp_box);
   dib_shift = ( -sy *(double)nkf2 - 0.5*(double)nkf2_cp_box);
   dic_shift = ( -sz *(double)nkf3 - 0.5*(double)nkf3_cp_box);

   dia_shift /= dnkf1;
   dib_shift /= dnkf2;
   dic_shift /= dnkf3;

   dia_shift -= NINT((dia_shift - .5));
   dib_shift -= NINT((dib_shift - .5));
   dic_shift -= NINT((dic_shift - .5)); 

   ia_shift = (int)(dia_shift*dnkf1);
   ib_shift = (int)(dib_shift*dnkf2);
   ic_shift = (int)(dic_shift*dnkf3);

   if( ib_shift < 0 || ic_shift < 0){
     ipack_on = 1;
   }

/*=======================================================================*/
/* The boxes are stationary  this means the map only needs to be constructed */
/*   the very first time                                                     */


 if(iset_map_flag == 0){
   for(iproc = 1; iproc<= num_proc ; iproc++){
      send_counts_row_big_small[iproc]      = 0;
      recv_counts_row_big_small[iproc]      = 0;
      send_counts_dual_map[iproc]           = 0;
      recv_counts_dual_map[iproc]           = 0;
   }
 }


   isend_tot         = 0;
   isend_tot_map     = 0;
   icount            = 0;
   for(kc_sm=skc_fft_ka_proc_cp_box;kc_sm<=ekc_fft_ka_proc_cp_box;kc_sm++){
     kb_str_sm = (kc_sm==skc_fft_ka_proc_cp_box ? skb_fft_ka_proc_cp_box : 1);
     kb_end_sm = (kc_sm==ekc_fft_ka_proc_cp_box ? ekb_fft_ka_proc_cp_box 
                                                : nkf2_cp_box);
     for(kb_sm=kb_str_sm;kb_sm<=kb_end_sm;kb_sm++){
       kc_lg = kc_sm + ic_shift;
       kc_lg = (kc_lg > 0 ? kc_lg : kc_lg + nkf3);
       kc_lg = (kc_lg > nkf3 ? kc_lg -nkf3 : kc_lg );

       for(iproc=1;iproc<=num_proc;iproc++){
         if( (kc_lg>=skc_fft_ka_proc_all_lg[iproc]) && 
             (kc_lg<=ekc_fft_ka_proc_all_lg[iproc])){
           kb_lg     = kb_sm + ib_shift;
           kb_lg     = (kb_lg > 0 ? kb_lg : kb_lg + nkf2);
           kb_lg     = (kb_lg > nkf2 ? kb_lg -nkf2 : kb_lg );
           kb_str_lg = (kc_lg==skc_fft_ka_proc_all_lg[iproc] ? 
                               skb_fft_ka_proc_all_lg[iproc] : 1);
           kb_end_lg = (kc_lg==ekc_fft_ka_proc_all_lg[iproc] ? 
                               ekb_fft_ka_proc_all_lg[iproc] : nkf2);

           if((kb_lg>=kb_str_lg) && (kb_lg<=kb_end_lg)){

             kc_off = kc_lg-skc_fft_ka_proc_all_lg[iproc];
             kb_off = kb_lg-kb_str_lg;

             index_lg  = kb_off*nkf1;

             if(kc_off>0){
	        index_lg += (kc_off-1)*nkf2*nkf1
		         +  (nkf2-skb_fft_ka_proc_all_lg[iproc]+1)*nkf1;  
             }/*endif*/

           if(ipack_on == 0 && iset_map_flag == 0){
	    ind_tmp[1+isend_tot_map] = index_lg; 
	  }

            for(ka_sm=1;ka_sm<=nkf1_cp_box;ka_sm++){
               ka_lg =  ka_sm + ia_shift;
               ka_lg = (ka_lg > 0 ? ka_lg : ka_lg + nkf1);
               ka_lg = (ka_lg > nkf1 ? ka_lg -nkf1 : ka_lg);

             if(ipack_on == 0){
    	        zfft[ka_sm+isend_tot]    = rfft[ka_sm+icount];   
	     }/*endif ipack*/

            }/*endfor*/
          if(iset_map_flag == 0){
             send_counts_row_big_small[iproc] += nkf1_cp_box;
             send_counts_dual_map[iproc] += 1;
	   }
             isend_tot                  += nkf1_cp_box;
             isend_tot_map              += 1;

	   }/*endif : kb in range*/
         }/*endif : kc in range */
       }/*endfor : iproc */
       icount  += nkf1_cp_box;
     }/*endfor : kb*/
   }/*endfor : kc */

/*=======================================================================*/
/* AllReduce the sum_counts_temp                                         */
/*  If send_counts does not equal the total number of elements of rho sm */
/*   there is an error                                                   */

 if(iset_map_flag == 0){

    sum_count_temp = 0;
  for(iproc=1; iproc<= num_proc; iproc++){
    sum_count_temp += send_counts_row_big_small[iproc];
   }

  if(np_states > 1){
    Allreduce(&(sum_count_temp),&(sum_count),
              1,MPI_INT,MPI_SUM,master,comm_lg);
  }else{
    sum_count = sum_count_temp;
  }
  printf("1 sum_count %d  nfft2_up_cp_box %d \n",sum_count,nfft2_up_cp_box);
  if( sum_count != nfft2_up_cp_box ){
   fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@ \n\n");
   fprintf(stderr,"Mismatch in the  number of values of rho to \n");
   fprintf(stderr,"be communicated in the routine that maps the\n");
   fprintf(stderr," density from the small to large grid       \n"); 
   fprintf(stderr," contact technical support\n");
   fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@\n");
   fflush(stderr);
   exit(1);
  }/*endif*/

/*=======================================================================*/
/* I) Alltoall  the send_counts to the receive counts                    */


  if(np_states > 1){
   Alltoall(&(send_counts_row_big_small[1]),one,MPI_INT,
            &(recv_counts_row_big_small[1]),one,MPI_INT,comm_cp_box);

   Alltoall(&(send_counts_dual_map[1]),one,MPI_INT,
            &(recv_counts_dual_map[1]),one,MPI_INT,comm_cp_box);
  }else{
    recv_counts_row_big_small[1] = send_counts_row_big_small[1];
    recv_counts_dual_map[1]       = send_counts_dual_map[1];
  }

/*-----------------------------------------------------------------------*/
/* I.V) Construct the displacements for the alltoallv                    */

    rdispls_row_big_small[1] = 0;
    sdispls_row_big_small[1] = 0;

    rdispls_dual_map[1] = 0;
    sdispls_dual_map[1] = 0;

  for(iproc=2; iproc<=num_proc; iproc++){
   rdispls_row_big_small[iproc] = rdispls_row_big_small[iproc-1]
                                + recv_counts_row_big_small[iproc-1];

   sdispls_row_big_small[iproc] = sdispls_row_big_small[iproc-1] 
                                + send_counts_row_big_small[iproc-1];

   rdispls_dual_map[iproc] = rdispls_dual_map[iproc-1]
                          + recv_counts_dual_map[iproc-1];

   sdispls_dual_map[iproc] = sdispls_dual_map[iproc-1] 
                          + send_counts_dual_map[iproc-1];
  }

}/*endif iset_map_flag*/


/*-----------------------------------------------------------------------*/
/* If the wrap condition applies, ie ib_shift < 0 ic shift < 0 then      */
/* the zfft array is packed in a different order for the alltoallv       */
/* where rfft elements are rearranged    this guarantees proper ordering */
 
 if(ipack_on == 1){

   for(iproc=1;iproc<= num_proc;iproc++){
      pcount[iproc] = 0;
      pcount_map[iproc] = 0;
   }

   isend_tot         = 0;
   isend_tot_map     = 0;
   icount            = 0;

   for(kc_sm=skc_fft_ka_proc_cp_box;kc_sm<=ekc_fft_ka_proc_cp_box;kc_sm++){
     kb_str_sm = (kc_sm==skc_fft_ka_proc_cp_box ? skb_fft_ka_proc_cp_box : 1);
     kb_end_sm = (kc_sm==ekc_fft_ka_proc_cp_box ? ekb_fft_ka_proc_cp_box 
                                                : nkf2_cp_box);
     for(kb_sm=kb_str_sm;kb_sm<=kb_end_sm;kb_sm++){
       kc_lg = kc_sm + ic_shift;
       kc_lg = (kc_lg > 0 ? kc_lg : kc_lg + nkf3);

       kc_off_sm = kc_sm - skc_fft_ka_proc_cp_box;
       kb_off_sm = kb_sm - kb_str_sm;

       for(iproc=1;iproc<=num_proc;iproc++){
         ioff     = sdispls_row_big_small[iproc];
         ioff_map = sdispls_dual_map[iproc];

         if( (kc_lg>=skc_fft_ka_proc_all_lg[iproc]) && 
             (kc_lg<=ekc_fft_ka_proc_all_lg[iproc])){
           kb_lg     = kb_sm + ib_shift;
           kb_lg     = (kb_lg > 0 ? kb_lg : kb_lg + nkf2);

           kb_str_lg = (kc_lg==skc_fft_ka_proc_all_lg[iproc] ? 
                               skb_fft_ka_proc_all_lg[iproc] : 1);
           kb_end_lg = (kc_lg==ekc_fft_ka_proc_all_lg[iproc] ? 
                               ekb_fft_ka_proc_all_lg[iproc] : nkf2);

           if((kb_lg>=kb_str_lg) && (kb_lg<=kb_end_lg)){

             kc_off = kc_lg-skc_fft_ka_proc_all_lg[iproc];
             kb_off = kb_lg-kb_str_lg;

             index_lg  = kb_off*nkf1;

             if(kc_off>0){
	        index_lg += (kc_off-1)*nkf2*nkf1
		         +  (nkf2-skb_fft_ka_proc_all_lg[iproc]+1)*nkf1;  
             }/*endif*/

               pcount_map[iproc]++;
           if(iset_map_flag == 0){
	       ind_tmp[(ioff_map+pcount_map[iproc])] = index_lg; 
           }

            for(ka_sm=1;ka_sm<=nkf1_cp_box;ka_sm++){
               pcount[iproc]++;

               ka_lg =  ka_sm + ia_shift;
               ka_lg = (ka_lg > 0 ? ka_lg : ka_lg + nkf1);

    	        zfft[(ioff+pcount[iproc])] = rfft[ka_sm+icount];   
             }/*endfor*/
             isend_tot                  += nkf1_cp_box;
             isend_tot_map              += 1;
	   }/*endif : kb in range*/
         }/*endif : kc in range */
       }/*endfor : iproc */
       icount  += nkf1_cp_box;
     }/*endfor : kb*/
   }/*endfor : kc */
   
}/*endif ipack_on*/  
 
/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*II) Alltoallv the dual_map which holds the index for the placement of   */
/*       the small density onto the large density in parallel            */
 if(iset_map_flag == 0){
  if(np_states > 1){
   Alltoallv(&(ind_tmp[1]),&(send_counts_dual_map[1]),
             &(sdispls_dual_map[1]),MPI_INT,
             &(map_dual[1]),&(recv_counts_dual_map[1]),
             &(rdispls_dual_map[1]),MPI_INT,comm_lg); 
  }else{
    for(i=1; i<= map_count;i++){
      map_dual[i] = ind_tmp[i];
    }/*endfor*/

  }

   cfree(&(ind_tmp[1]));

 }/*endif iset_map_flag */

       map_count_proc = 0;
     for(iproc=1; iproc<= num_proc; iproc++){
       map_count_proc += recv_counts_dual_map[iproc];
     }/*endfor*/

/*-----------------------------------------------------------------------*/
/*III) Alltoallv the real space density on the small real space grid     */
/*  store resulting vector in zfft_tmp                                   */

 if(np_states > 1 ){
   Alltoallv(&(zfft[1]),&(send_counts_row_big_small[1]),
             &(sdispls_row_big_small[1]),MPI_DOUBLE,
             &(zfft_tmp[1]),&(recv_counts_row_big_small[1]),
             &(rdispls_row_big_small[1]),MPI_DOUBLE,comm_lg);
 }else{
   for(i=1; i<= nfft2_up_cp_box;i++){
     zfft_tmp[i] = zfft[i];
   }/*endfor*/
 }

/*=======================================================================*/
/* Assign the small real space density elements to the appropriate       */
/*  on the large real space grid                                         */

  icount_tot = 0;

  for(iproc=1; iproc<= num_proc; iproc++){
    icount_tot += recv_counts_row_big_small[iproc];
  }


  for(i=1; i<= nfft_proc; i++){
    zfft[i] = 0.0;
  }
 
  
      index_sm = 1;
   for(i=1; i<= map_count_proc; i++){
    for(ka_sm=1 ;ka_sm <= nkf1_cp_box; ka_sm++){
     ka_lg = ka_sm + ia_shift;
     ka_lg = (ka_lg > 0 ? ka_lg : ka_lg + nkf1);
     index_lg_cmplx       = 2*(map_dual[i]+ka_lg) - 1;   
     zfft[index_lg_cmplx] = zfft_tmp[index_sm]*scale;
/* rfft =  rho*vol_cp, in the next section of cp_rho_calc,*/
/* rho is divided by vol. Sinc vol = (box_rat)^3*vol_cp,  */
/* you need to multiply rfft by scale = (box_rat)^3 here  */

     index_sm++;
    }/*endfor ka_sm*/
   }/*endfor*/


/*-----------------------------------------------------------------------*/
/* free locally assigned memory                                          */
 
   cfree(&(pcount[1]));
   cfree(&(pcount_map[1]));

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_pack_rho_dual_par_pme(CPSCR_WAVE *cpscr_wave,double *rfft,
                                CELL *cell,double dbox_rat,
                                int n_interp_pme_dual,
                                int iset_map_flag,
                                CPSCR_DUAL_PME *cpscr_dual_pme,
                                PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                                PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int m,i,iii;
 int ja,jb,jc; 
 int ka,kb,kc;
 int irow;
 int joff;

 int icount;
 int kb_str_sm,kb_end_sm;

 double scale,trap_corr;
 double temp;
 double di_nkf1_cp,di_nkf2_cp,di_nkf3_cp;
 double dnkf1,dnkf2,dnkf3;

/* Local pointers */

 double *zfft,*zfft_tmp;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

/* PME MAP variables */
 int **igrid_a;
 int **igrid_b;
 int **igrid_c;

 double **mn_a;
 double **mn_b;
 double **mn_c;

 int nfft_send_big_small;
 int num_rows_tot;

 int **joff_sm_lg;
 int *ioff_lg_sm;

 int *send_counts_row_big_small = 
                         para_fft_pkg3d_cp_box->send_counts_row_big_small;
 int *recv_counts_row_big_small =
                        para_fft_pkg3d_cp_box->recv_counts_row_big_small;

 int *sdispls_row_big_small  = para_fft_pkg3d_cp_box->sdispls_row_big_small;
 int *rdispls_row_big_small  = para_fft_pkg3d_cp_box->rdispls_row_big_small;


/* Small CP Box                                                          */

 int   skc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekb_fft_ka_proc;

 int nfft_proc        = para_fft_pkg3d_lg->nfft_proc;   /* on large grid*/

 int nproc            = para_fft_pkg3d_lg->num_proc;
 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 
 int myid             = para_fft_pkg3d_lg->myid;

/*=======================================================================*/
/* I) Determine useful constants                                         */

  scale = dbox_rat*dbox_rat*dbox_rat;

  di_nkf1_cp = 1.0/((double) (nkf1_cp_box));
  di_nkf2_cp = 1.0/((double) (nkf2_cp_box));
  di_nkf3_cp = 1.0/((double) (nkf3_cp_box));

  dnkf1 = (double) nkf1;  /*large sparse grid*/
  dnkf2 = (double) nkf2;  /*large sparse grid*/
  dnkf3 = (double) nkf3;  /*large sparse grid*/

  trap_corr = (dnkf1*dnkf2*dnkf3*di_nkf1_cp*di_nkf2_cp*di_nkf3_cp)
             /(dbox_rat*dbox_rat*dbox_rat);

/*=============================================================================*/
/* Make the  PME map  if cp_box is stationary only make map the very first time */

 if(iset_map_flag == 0){ /*create maps */

   make_dual_pme_wghts(cpscr_dual_pme,cell,dbox_rat,
                       n_interp_pme_dual,
                       para_fft_pkg3d_cp_box,para_fft_pkg3d_lg);

   make_pme_para_dual_map(n_interp_pme_dual,cpscr_wave,
                          cpscr_dual_pme,
                          para_fft_pkg3d_cp_box,para_fft_pkg3d_lg);

 }/*endif set maps*/

/* assign local pointers */

 igrid_a = cpscr_dual_pme->igrid_a; 
 igrid_b = cpscr_dual_pme->igrid_b;
 igrid_c = cpscr_dual_pme->igrid_c;

 mn_a  = cpscr_dual_pme->mn_a;  
 mn_b  = cpscr_dual_pme->mn_b;
 mn_c  = cpscr_dual_pme->mn_c;

 nfft_send_big_small = cpscr_dual_pme->nfft_send_big_small;
 num_rows_tot        = cpscr_dual_pme->num_rows_tot;

 joff_sm_lg   = cpscr_dual_pme->joff_sm_lg;
 ioff_lg_sm   = cpscr_dual_pme->ioff_lg_sm;

 zfft     = cpscr_wave->zfft; 
 zfft_tmp = cpscr_wave->zfft_tmp; 

/*=======================================================================*/
/* C) Interpolate the small grid to the large using the  Mn's weighting  */
/*    factors */

  for(i=1; i<=nfft_send_big_small; i++){
    zfft[i]     = 0.0;
  }/*endfor*/


  for(jc=1;jc<=n_interp_pme_dual;jc++){
  for(jb=1;jb<=n_interp_pme_dual;jb++){

    icount  = 0;
    for(kc=skc_fft_ka_proc_cp_box;kc<=ekc_fft_ka_proc_cp_box;kc++){

      kb_str_sm = (kc==skc_fft_ka_proc_cp_box ? skb_fft_ka_proc_cp_box : 1);
      kb_end_sm = (kc==ekc_fft_ka_proc_cp_box ? ekb_fft_ka_proc_cp_box 
                                              : nkf2_cp_box);

      for(kb=kb_str_sm;kb<=kb_end_sm;kb++){
       /* rfft =  rho*vol_cp, in the next section of cp_rho_calc,       */
       /* rho is divided by vol. Sinc vol = (box_rat)^3*vol_cp,         */
       /* you need to multiply rfft by scale = (box_rat)^3 here         */
       /* trap_corr is needed because we are changing the grid spacing  */
       /*  of the density and we still want the integral of the density */
       /* to generate the number of electrons   */

        temp = mn_b[jb][kb]*mn_c[jc][kc]*scale*trap_corr;

       /* All this goes to the same processor. That is, full rows of ka  */
       /* are stored on the same proc. joff says where this row should   */
       /* be stored in zfft so that the alltoallv will work properly.    */

         joff = joff_sm_lg[igrid_c[jc][kc]][igrid_b[jb][kb]];

         for(ja=1;ja<=n_interp_pme_dual;ja++){
         for(ka=1; ka <= nkf1_cp_box ; ka++){
          zfft[(igrid_a[ja][ka] + joff)] += 
	      (mn_a[ja][ka]*temp*rfft[(ka+icount)]); 

         }}/* endfor : ka,ja */
         icount += nkf1_cp_box;

      }/*endfor :  kb */
    }/*endfor :  kc */

  }}/*endfor : jb,jc*/

/*=======================================================================*/
/* Communicate zfft across processors to zfft_tmp */
 
 if(nproc > 1){
     Alltoallv(&(zfft[1]),&(send_counts_row_big_small[1]),
             &(sdispls_row_big_small[1]),MPI_DOUBLE,
             &(zfft_tmp[1]),&(recv_counts_row_big_small[1]),
             &(rdispls_row_big_small[1]),MPI_DOUBLE,comm_lg);

 }else{
   for(i=1; i<= nfft_proc; i++){zfft_tmp[i] = zfft[i];}
 }

/*=======================================================================*/
/* Unpack zfft_tmp into zfft on the large sparse grid */

  for(i=1; i<= nfft_proc; i++){zfft[i] = 0.0;}


  /* unpack the rows into zfft using the offset */
  icount = 0;
  for(irow =1; irow <= num_rows_tot; irow++){

   for(ka=1,m=1; m<= nkf1; ka+=2,m++){
     zfft[ioff_lg_sm[irow] + ka]  += zfft_tmp[m+icount];
   }/*endfor ka*/
   icount+= nkf1;

  }/*endfor*/


/*==========================================================================*/
/*-----------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Routine to create the PME weights for the dual grid option               */
/*==========================================================================*/

void make_dual_pme_wghts(CPSCR_DUAL_PME *cpscr_dual_pme,CELL *cell,
                       double dbox_rat,int n_interp_pme_dual,
                       PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                       PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
 int i,iii;

 int ia,ib,ic;
 int j,j1,j2,n;

 double sa_shift,sb_shift,sc_shift;
 double atemp,btemp,ctemp;
 double dia1,dib1,dic1;

 double di_nkf1_cp,di_nkf2_cp,di_nkf3_cp;
 double dnkf1,dnkf2,dnkf3;
 double mn_a_tmp,mn_b_tmp,mn_c_tmp;

/* Local pointers */
 int  np_states  = para_fft_pkg3d_lg->num_proc;

/* on large sparse grid*/
 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 double *hmati             = cell->hmati;
 double *cp_box_center     = cell->cp_box_center;

/*length nkf1(2,3)_cp_box*/
 double *a_pme = cpscr_dual_pme->a_pme; 
 double *b_pme = cpscr_dual_pme->b_pme;
 double *c_pme = cpscr_dual_pme->c_pme;

 int *iatemp = cpscr_dual_pme->iatemp;
 int *ibtemp = cpscr_dual_pme->ibtemp;
 int *ictemp = cpscr_dual_pme->ictemp;

/*length nkf1(2,3)_cp_box*/
 int **igrid_a = cpscr_dual_pme->igrid_a; 
 int **igrid_b = cpscr_dual_pme->igrid_b;
 int **igrid_c = cpscr_dual_pme->igrid_c;

/*length nkf1(2,3)_cp_box*/
 double *frac_a = cpscr_dual_pme->frac_a; 
 double *frac_b = cpscr_dual_pme->frac_b;
 double *frac_c = cpscr_dual_pme->frac_c;

 /*length: n_interp_pme_dual X nkf1(2,3)_cp_box*/
 double **mn_a  = cpscr_dual_pme->mn_a; 
 double **mn_b  = cpscr_dual_pme->mn_b;
 double **mn_c  = cpscr_dual_pme->mn_c;
 double **ua    = cpscr_dual_pme->ua;
 double **ub    = cpscr_dual_pme->ub;
 double **uc    = cpscr_dual_pme->uc;

/*length n_interp_pme_dual*/
 double *aj     = cpscr_dual_pme->aj;   
 double *rn     = cpscr_dual_pme->rn;
 double *rn1    = cpscr_dual_pme->rn1;

 int nproc            = para_fft_pkg3d_lg->num_proc;
 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 
 int myid             = para_fft_pkg3d_lg->myid;

 int ic_old;

/*=======================================================================*/
/* Determine useful constants                                            */

   sa_shift = cp_box_center[1]*hmati[1]+cp_box_center[2]*hmati[4] 
            + cp_box_center[3]*hmati[7];

   sb_shift = cp_box_center[1]*hmati[2]+cp_box_center[2]*hmati[5]
            + cp_box_center[3]*hmati[8];

   sc_shift = cp_box_center[1]*hmati[3]+cp_box_center[2]*hmati[6]
            + cp_box_center[3]*hmati[9];


    di_nkf1_cp = 1.0/((double) (nkf1_cp_box));
    di_nkf2_cp = 1.0/((double) (nkf2_cp_box));
    di_nkf3_cp = 1.0/((double) (nkf3_cp_box));

    dnkf1 = (double) nkf1;  /*large sparse grid*/
    dnkf2 = (double) nkf2;  /*large sparse grid*/
    dnkf3 = (double) nkf3;  /*large sparse grid*/


/*=======================================================================*/

   for(j=1;j<=n_interp_pme_dual;j++){
     aj[j] = (double) (j-1);
     rn[j] = (double) (j);
     if(j > 1){rn1[j] = 1.0/((double)(j-1));}
   }/*endfor*/
     rn1[1] = 0.0;

/*=======================================================================*/
/* I) For each small dense grid point, get its scaled coordinate         */
/*     in the big box frame */

/* Scaled coordinate of A */
    for(ia=1; ia <= nkf1_cp_box; ia++){
      dia1 = (double)(ia-1);
      atemp = (dia1*di_nkf1_cp - 0.5)/(dbox_rat) - sa_shift; 
      atemp = atemp - NINT((atemp-0.5));
      atemp = (atemp >= 1.0 ? 0.0 : atemp );
      a_pme[ia] = atemp*dnkf1; 
    }


/* Scaled coordinate of B */
    for(ib=1; ib <= nkf2_cp_box; ib++){
      dib1 = (double)(ib-1);
      btemp = (dib1*di_nkf2_cp - 0.5)/(dbox_rat) - sb_shift; 
      btemp = btemp - NINT((btemp-0.5));
      btemp = (btemp >= 1.0 ? 0.0 : btemp );
      b_pme[ib] = btemp*dnkf2; 
    }

/* Scaled coordinate of C */
    for(ic=1; ic <= nkf3_cp_box; ic++){
      dic1 = (double)(ic-1);
      ctemp = (dic1*di_nkf3_cp - 0.5)/(dbox_rat) - sc_shift; 
      ctemp = ctemp - NINT((ctemp-0.5));
      ctemp = (ctemp >= 1.0 ? 0.0 : ctemp );
      c_pme[ic] = ctemp*dnkf3; 
    }


/*==========================================================================*/
/*==========================================================================*/
/* III) Calculate the Cardinal B spline interpolation functions of the      */
/*     scaled lattice vectors                                               */ 
/*--------------------------------------------------------------------------*/ 
/* A) Using current fraction, find the grid points on which M_n is non-zero */


/*---- A ---- */

     for(i=1;i<=nkf1_cp_box;i++){
      iatemp[i] = (int) (a_pme[i]);
      frac_a[i] = a_pme[i] - (double) (iatemp[i]);
     }/*endfor*/

     for(j=1;j<=n_interp_pme_dual;j++){
     if(np_states == 1){
      for(i=1;i<= nkf1_cp_box;i++){
       ua[j][i]  = frac_a[i] + aj[j];
       j2        = j-2;
       ia        = iatemp[i] - j2;
       ia        = (ia>0 ? ia:nkf1 + ia);

#ifdef REAL_CASE
       igrid_a[j][i] = 2*ia - 1;
#else
       igrid_a[j][i] = ia;
#endif
      }/*endfor*/
     }else{
      for(i=1;i<= nkf1_cp_box;i++){
       ua[j][i]  = frac_a[i] + aj[j];
       j2        = j-2;
       ia        = iatemp[i] - j2;
       ia        = (ia>0 ? ia:nkf1 + ia);
       igrid_a[j][i] = ia;
      }/*endfor*/
      }/*endif np_states */
     }/*endfor*/

/*---- B ---- */

     for(i=1;i<=nkf2_cp_box;i++){
      ibtemp[i] = (int) (b_pme[i]);
      frac_b[i] = b_pme[i] - (double) (ibtemp[i]);
     }/*endfor*/

     for(j=1;j<=n_interp_pme_dual;j++){
     if(np_states == 1){
      for(i=1;i<= nkf2_cp_box;i++){
       ub[j][i]  = frac_b[i] + aj[j];
       j2        = j-2;
       ib        = ibtemp[i] - j2;
       ib        = (ib>0 ? ib:nkf2 + ib);
       
#ifdef REAL_CASE
       igrid_b[j][i] = 2*(ib - 1)*nkf1;
#else
       igrid_b[j][i] = ib;
#endif
      }/*endfor*/
     }else{
      for(i=1;i<= nkf2_cp_box;i++){
       ub[j][i]  = frac_b[i] + aj[j];
       j2        = j-2;
       ib        = ibtemp[i] - j2;
       ib        = (ib>0 ? ib:nkf2 + ib);
       igrid_b[j][i] = ib;
      }/*endfor*/
     }/*endif np_states*/
     }/*endfor*/

/*---- C ---- */
     for(i=1;i<=nkf3_cp_box;i++){
      ictemp[i] = (int) (c_pme[i]);
      frac_c[i] = c_pme[i] - (double) (ictemp[i]);
     }/*endfor*/

     for(j=1;j<=n_interp_pme_dual;j++){
      if(np_states == 1){
      for(i=1;i<= nkf3_cp_box;i++){
       uc[j][i]  = frac_c[i] + aj[j];
       j2        = j-2;
       ic        = ictemp[i] - j2;
       ic        = (ic>0 ? ic:nkf3 +ic);
       ic        = (ic<=nkf3 ? ic:-nkf3+ic);

#ifdef REAL_CASE
       igrid_c[j][i] = 2*(ic - 1)*nkf1*nkf2;
#else
       igrid_c[j][i] = ic;
#endif
      }/*endfor*/
      }else{
      for(i=1;i<= nkf3_cp_box;i++){
       uc[j][i]  = frac_c[i] + aj[j];
       j2        = j-2;
       ic        = ictemp[i] - j2;
       ic        = (ic>0 ? ic:nkf3 +ic);
       ic        = (ic<=nkf3 ? ic:-nkf3+ic);
       igrid_c[j][i] = ic;
      }/*endfor*/
      }/*endif np_states*/
     }/*endfor*/
 
/*==========================================================================*/
/*--------------------------------------------------------------------------*/ 
/* B) Initialize M2 and get the Mn's using the recursion relation           */ 
/*    Note: M_n is defined on 0 to n. Since frac between 0 and 1, the       */
/*          calculation is performed in an order that takes advantage of it */

/*---- A ----*/
     for(i=1;i<=nkf1_cp_box;i++){
      mn_a[1][i] = 1.0 - fabs(ua[1][i]-1.0);
      mn_a[2][i] = 1.0 - fabs(ua[2][i]-1.0);
     }/*endfor*/

     for(j=3;j<=n_interp_pme_dual;j++){
      for(i=1;i<=nkf1_cp_box;i++){
       mn_a[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp_pme_dual;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nkf1_cp_box;i++){
         mn_a_tmp  = (ua[j][i]*mn_a[j][i]+(rn[n]-ua[j][i])*mn_a[j1][i])*rn1[n];
         mn_a[j][i] = mn_a_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nkf1_cp_box;i++){
        mn_a[1][i] = ua[1][i]*mn_a[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*---- B ----*/
     for(i=1;i<=nkf2_cp_box;i++){
      mn_b[1][i] = 1.0 - fabs(ub[1][i]-1.0);
      mn_b[2][i] = 1.0 - fabs(ub[2][i]-1.0);
     }/*endfor*/

     for(j=3;j<=n_interp_pme_dual;j++){
      for(i=1;i<=nkf2_cp_box;i++){
       mn_b[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp_pme_dual;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nkf2_cp_box;i++){
         mn_b_tmp  = (ub[j][i]*mn_b[j][i]+(rn[n]-ub[j][i])*mn_b[j1][i])*rn1[n];
         mn_b[j][i] = mn_b_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nkf2_cp_box;i++){
        mn_b[1][i] = ub[1][i]*mn_b[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*---- C ----*/
     for(i=1;i<=nkf3_cp_box;i++){
      mn_c[1][i] = 1.0 - fabs(uc[1][i]-1.0);
      mn_c[2][i] = 1.0 - fabs(uc[2][i]-1.0);
     }/*endfor*/

     for(j=3;j<=n_interp_pme_dual;j++){
      for(i=1;i<=nkf3_cp_box;i++){
       mn_c[j][i]   = 0.0;
      }/*endfor*/
     }/*endfor*/

     for(n=3;n<=n_interp_pme_dual;n++){
       for(j=n;j>=2;j--){
        j1 = j-1;
        for(i=1;i<=nkf3_cp_box;i++){
         mn_c_tmp  = (uc[j][i]*mn_c[j][i]+(rn[n]-uc[j][i])*mn_c[j1][i])*rn1[n];
         mn_c[j][i] = mn_c_tmp;
        }/*end for: i*/
       }/*end for: j*/
       for(i=1;i<=nkf3_cp_box;i++){
        mn_c[1][i] = uc[1][i]*mn_c[1][i]*rn1[n];
      }/*endfor */
     }/*end for: k*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Routine to create the parallel maps for dual grid PME option             */
/*==========================================================================*/

void make_pme_para_dual_map(int n_interp_pme_dual,
                            CPSCR_WAVE *cpscr_wave,
                            CPSCR_DUAL_PME *cpscr_dual_pme,
                            PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                            PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int iproc,jproc,mem_now;
 int icase,nka_hit,iii,i;
 int kc_lg,kb_lg;
 int kc,kb,jc,jb;
 int kc_hit_mid,kc_hit_high,kc_hit_low,kb_hit;
 int nfft_recv_big_small,imax1;
 int num_rows_tot;
 int one = 1;

/*          Local pointer declarations                                  */

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;

 int   skc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekb_fft_ka_proc;

 /* remove after debug*/
 int   skc_fft_ka_proc_lg = para_fft_pkg3d_lg->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_lg = para_fft_pkg3d_lg->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_lg = para_fft_pkg3d_lg->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_lg = para_fft_pkg3d_lg->ekb_fft_ka_proc;
 /* remove after debug*/

 int *skb_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->skb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->skc_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->ekb_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all_lg  = para_fft_pkg3d_lg->ekc_fft_ka_proc_all;

 int nproc = para_fft_pkg3d_lg->num_proc;

 int *send_counts_big_small = 
                         para_fft_pkg3d_cp_box->send_counts_row_big_small;
 int *recv_counts_big_small =
                         para_fft_pkg3d_cp_box->recv_counts_row_big_small;

 int *sdispls_big_small = 
                         para_fft_pkg3d_cp_box->sdispls_row_big_small;
 int *rdispls_big_small =
                         para_fft_pkg3d_cp_box->rdispls_row_big_small;


 int *send_counts_row_big_small
                     = para_fft_pkg3d_cp_box->send_counts_ioff_big_small;
 int *recv_counts_row_big_small 
                     = para_fft_pkg3d_cp_box->recv_counts_ioff_big_small;
 int *sdispls_row_big_small
                     = para_fft_pkg3d_cp_box->sdispls_ioff_big_small;
 int *rdispls_row_big_small
                     = para_fft_pkg3d_cp_box->rdispls_ioff_big_small;

 int zfft_mall_size  =  cpscr_wave->zfft_mall_size;

/* PME MAP variables */
 int **joff_sm_lg;
 int *ioff_lg_sm;   
 int *ioff_lg_sm_tmp;

  /*length nkf1_cp_box*/
 int **igrid_b = cpscr_dual_pme->igrid_b;
 int **igrid_c = cpscr_dual_pme->igrid_c;

 int nfft_send_big_small = cpscr_dual_pme->nfft_send_big_small;
 int myid             = para_fft_pkg3d_lg->myid;
 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 

/*=============================================================================*/
/* I)Determine which kc and kb on the big grid are hit by the PME procedure   */
/*   from the portion of the small grid stored on this processor. Also,       */
/*   determine which processor these points must go to and where in the array */
/*   they must be stored so that an Alltoallv can be performed.               */
/*   The method takes advantage of the separablility of the different         */
/*   dimensions.                                                              */

  for(iproc=1; iproc<= nproc; iproc++){
    send_counts_big_small[iproc] = 0;
  }/*endfor*/

/*-------------------------------------------------------------------*/
/* For now we are looping over all possible kc_lg, kb_lg             */
/* This can be restricted later because all possible kc_lg and kb_lg */
/* cannot be reached from the portion of the small grid on this proc */
/* or indeed the small grid in general. Later, we will make the      */
/* routine more complex by looping over a restricted range of kc_lg  */
/* and kb_lg.                                                        */
/*-------------------------------------------------------------------*/

 cpscr_dual_pme->joff_sm_lg = cmall_int_mat(1,nkf3,1,nkf2);

 cpscr_dual_pme->ioff_lg_sm = (int *) 
                                 cmalloc((nkf3*nkf2)*sizeof(int))-1;

 /* assign local pointers only after you do the mallocs*/
 joff_sm_lg = cpscr_dual_pme->joff_sm_lg;
 ioff_lg_sm = cpscr_dual_pme->ioff_lg_sm;

 jproc   = 1;
 mem_now = 0;
 nka_hit  = 0;

  for(kc_lg=1; kc_lg <= nkf3; kc_lg++){

    /*---------------------------------------------------------------*/
    /* Is there a mid point value with a hit on this big grid point? */
    kc_hit_mid = 0;

    for(kc=skc_fft_ka_proc_cp_box+1;kc<=ekc_fft_ka_proc_cp_box-1;kc++){
       for(jc=1; jc <= n_interp_pme_dual; jc++){
        if( igrid_c[jc][kc] == kc_lg ){
          kc_hit_mid = 1;  /* found a kc_sm that maps */
          break;
        }
       }/*endfor jc*/
        if(kc_hit_mid == 1){break;} /*once I find one I break */
    }/*endfor*/

    /*-------------------------------------*/
    /* Does the low kc_sm value have a hit */
    kc_hit_low = 0;
    kc = skc_fft_ka_proc_cp_box;
       for(jc=1; jc <= n_interp_pme_dual; jc++){
        if( igrid_c[jc][kc] == kc_lg ){
          kc_hit_low = 1;
          break;
        }
       }/*endfor jc*/

    /*--------------------------------------*/
    /* Does the high kc_sm value have a hit */
    kc_hit_high = 0;
    kc = ekc_fft_ka_proc_cp_box;
       for(jc=1; jc <= n_interp_pme_dual; jc++){
        if( igrid_c[jc][kc] == kc_lg ){
          kc_hit_high = 1;
          break;
        }
       }/*endfor jc*/

    /*----------------------------------------------------------------------*/
    /* Determine how many kb_sm values one must loop over to determine hits */
    /* on kc_lg and kb_lg */
    if( kc_hit_mid == 0 && kc_hit_low == 0 && kc_hit_high == 0 ){ icase = 0;}
    if( kc_hit_mid > 0)                                         { icase = 1;}
    if( kc_hit_mid == 0 && kc_hit_low == 1 && kc_hit_high == 0) { icase = 2;}
    if( kc_hit_mid == 0 && kc_hit_low == 0 && kc_hit_high == 1) { icase = 3;}
    if( kc_hit_mid == 0 && kc_hit_low == 1 && kc_hit_high == 1 &&
         skc_fft_ka_proc_cp_box != ekc_fft_ka_proc_cp_box     ) { icase = 4;}
    if( kc_hit_mid == 0 && kc_hit_low == 1 && kc_hit_high == 1 &&
        skc_fft_ka_proc_cp_box == ekc_fft_ka_proc_cp_box     )  { icase = 5;}

    /*-----------------------------------------------------------------*/
    /* loop over the kb_lg and determine if there is a kb_sm that hits */
    for(kb_lg=1; kb_lg <= nkf2; kb_lg++){
      if(icase > 0 ){
         kb_hit = 0;
         switch(icase){
           case 1: 
             for(kb=1; kb <= nkf2_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                 if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; } /* You got a hit */
             }
           break;
           case 2: 
             for(kb=skb_fft_ka_proc_cp_box; kb <= nkf2_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                 if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; }/* You got a hit */
             }
           break;
           case 3: 
             for(kb=1; kb <= ekb_fft_ka_proc_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                 if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; }/* You got a hit */
             }
           break;
           case 4: 
             for(kb=skb_fft_ka_proc_cp_box; kb <= nkf2_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                 if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; }/* You got a hit */
             }
             for(kb=1; kb <= ekb_fft_ka_proc_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; }/* You got a hit */
             }
           break;
           case 5: 
             for(kb=skb_fft_ka_proc_cp_box; kb <= ekb_fft_ka_proc_cp_box; kb++){
               for(jb=1; jb <= n_interp_pme_dual; jb++){
                 if(kb_lg == igrid_b[jb][kb]){kb_hit = 1;break;}
               }
               if(kb_hit == 1){ break; }/* You got a hit */
             }
           break;
         }/*end switch*/
         /*-----------------------------------------*/
         /* if you have a hit, thanks for the memory*/
         if(kb_hit == 1){
           joff_sm_lg[kc_lg][kb_lg]    = mem_now;  

           mem_now                     += nkf1;

           send_counts_big_small[jproc] += nkf1;
           nka_hit++;
           if(kc_lg > skc_fft_ka_proc_all_lg[jproc]){
            ioff_lg_sm[nka_hit]  = ((nkf2-skb_fft_ka_proc_all_lg[jproc]+1)
                                  + (kc_lg-skc_fft_ka_proc_all_lg[jproc]-1)*nkf2
                                  + (kb_lg-1) )*2*nkf1;
	  }else{
            ioff_lg_sm[nka_hit] = (kb_lg-skb_fft_ka_proc_all_lg[jproc])*2*nkf1;
          }/*endif*/
    
         }/*endif : I have a hit */

      }/*endif icase : I can get a hit*/
  
     /*-----------------------------------------------*/
     /* check to see if you are on the next processor */
     if(kc_lg == ekc_fft_ka_proc_all_lg[jproc] &&
	  kb_lg >= ekb_fft_ka_proc_all_lg[jproc] ){jproc++;}
                   
    }/*endfor kb_lg*/
  }/*endfor kc_lg*/

  nfft_send_big_small = mem_now;
  cpscr_dual_pme->nfft_send_big_small = nfft_send_big_small;

/*=======================================================================*/
/* I) Alltoall  the send_counts to the receive counts                    */

/*-----------------------------------------------------------------------*/
/* i) Make send displs */

  sdispls_big_small[1] = 0;
  for(iproc=2; iproc<= nproc; iproc++){
      sdispls_big_small[iproc] = sdispls_big_small[iproc-1] 
                             + send_counts_big_small[iproc-1];  
  }/*endfor*/


/*Create number of send counts and sdispls for ioff_lg_sm  */

  for(iproc=1; iproc <=nproc; iproc++){
   send_counts_row_big_small[iproc] = send_counts_big_small[iproc]/nkf1;
   sdispls_row_big_small[iproc]     = sdispls_big_small[iproc]/nkf1;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
/* ii) Create receive count  */

  if(nproc > 1){
  Barrier(comm_lg);
   Alltoall(&(send_counts_big_small[1]),one,MPI_INT,
            &(recv_counts_big_small[1]),one,MPI_INT,comm_lg);
  Barrier(comm_lg);
  }else{
    recv_counts_big_small[1] = send_counts_big_small[1];
  }

   nfft_recv_big_small = 0;

   for(iproc=1; iproc <= nproc; iproc++){
      nfft_recv_big_small += recv_counts_big_small[iproc];
   }/*endfor*/

   cpscr_dual_pme->nfft_recv_big_small = nfft_recv_big_small;

/*-----------------------------------------------------------------------*/
/* iii) Create receive displs */

  rdispls_big_small[1] = 0;

  for(iproc=2; iproc<=nproc; iproc++){
   rdispls_big_small[iproc] = rdispls_big_small[iproc-1]
                            + recv_counts_big_small[iproc-1];
  }/*endfor*/


  for(iproc = 1; iproc <= nproc; iproc++){
    rdispls_row_big_small[iproc]     = rdispls_big_small[iproc]/nkf1;
    recv_counts_row_big_small[iproc] = recv_counts_big_small[iproc]/nkf1;
  }

/*-----------------------------------------------------------------------*/
/* iv) Put row offsets onto correct processors */

  if(nproc > 1){
    ioff_lg_sm_tmp = (int *) cmalloc((nkf3*nkf2)*sizeof(int))-1;
  
    for(i=1; i<= (nkf3*nkf2); i++){
      ioff_lg_sm_tmp[i] = ioff_lg_sm[i];
    }

    Alltoallv(&(ioff_lg_sm_tmp[1]),&(send_counts_row_big_small[1]),
              &(sdispls_row_big_small[1]),MPI_INT,
              &(ioff_lg_sm[1]),&(recv_counts_row_big_small[1]),
              &(rdispls_row_big_small[1]),MPI_INT,comm_lg);

    cfree(&(ioff_lg_sm_tmp[1]));

  }/*endif nproc*/

/*=========================================================================*/
/* III) Reallocate zfft and zfft_tmp if necessary                           */

 imax1 = MAX3(nfft_recv_big_small,nfft_send_big_small,zfft_mall_size); 

 if(imax1>zfft_mall_size){

     printf("Reallocating zfft for pme dual grid from %d to %d on proc %d \n",
           zfft_mall_size,imax1,myid+1);

    cpscr_wave->zfft_mall_size = imax1;
    cpscr_wave->zfft     = (double *) crealloc(&((cpscr_wave->zfft[1])),
                                              imax1*sizeof(double))-1;

    cpscr_wave->zfft_tmp = (double *) crealloc(&((cpscr_wave->zfft_tmp[1])),
                                              imax1*sizeof(double))-1; 

 }/*endif*/


  /* rows of length nkf1 sent to this proc */
  num_rows_tot  = 0;
  for(iproc=1; iproc <= nproc; iproc++){
    num_rows_tot += recv_counts_big_small[iproc];
  }
  num_rows_tot /= nkf1;  

  cpscr_dual_pme->num_rows_tot = num_rows_tot;


/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/









