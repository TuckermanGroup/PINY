/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_bond.c                                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_parse_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_bond(BONDED *bonded, NULL_INTER_PARSE *null_inter_parse)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/

   int i,ntyp_pow,nmol_typ_max,dim_k;
   int n_interp,nlen_pme,iii;

/*========================================================================*/
/* Local Pointer Sizes                                                    */

   int ngrp_21_mall          = bonded->grp_bond_con.ngrp_21_mall; 
   int ngrp_33_mall          = bonded->grp_bond_con.ngrp_33_mall;
   int ngrp_watt_33_mall     = bonded->grp_bond_watts.ngrp_33_mall;
   int ngrp_43_mall          = bonded->grp_bond_con.ngrp_43_mall;
   int ngrp_23_mall          = bonded->grp_bond_con.ngrp_23_mall;
   int ngrp_46_mall          = bonded->grp_bond_con.ngrp_46_mall; 
   int ngrp_typ_21_mall      = bonded->grp_bond_con.ngrp_typ_21_mall; 
   int ngrp_typ_33_mall      = bonded->grp_bond_con.ngrp_typ_33_mall;
   int ngrp_typ_watt_33_mall = bonded->grp_bond_watts.ngrp_typ_33_mall;
   int ngrp_typ_43_mall      = bonded->grp_bond_con.ngrp_typ_43_mall;
   int ngrp_typ_23_mall      = bonded->grp_bond_con.ngrp_typ_23_mall;
   int ngrp_typ_46_mall      = bonded->grp_bond_con.ngrp_typ_46_mall;

   int nbond_pow_mall        = bonded->bond.nbond_pow_mall;
   int nbond_typ_pow_mall    = bonded->bond.nbond_typ_pow_mall;
   int nbond_con_mall        = bonded->bond.nbond_con_mall;
   int nbond_typ_con_mall    = bonded->bond.nbond_typ_con_mall;
   int nbend_pow_mall        = bonded->bend.nbend_pow_mall;
   int nbend_typ_pow_mall    = bonded->bend.nbend_typ_pow_mall;
   int nbend_con_mall        = bonded->bend.nbend_con_mall;
   int nbend_typ_con_mall    = bonded->bend.nbend_typ_con_mall;
   int nbend_bnd_mall        = bonded->bend_bnd.nbend_bnd_mall;
   int nbend_bnd_typ_mall    = bonded->bend_bnd.nbend_bnd_typ_mall;
   int ntors_pow_mall        = bonded->tors.ntors_pow_mall;
   int ntors_typ_pow_mall    = bonded->tors.ntors_typ_pow_mall;
   int ntors_con_mall        = bonded->tors.ntors_con_mall;
   int ntors_typ_con_mall    = bonded->tors.ntors_typ_con_mall;
   int nonfo_mall            = bonded->onfo.nonfo_mall;
   int nonfo_typ_mall        = bonded->onfo.nonfo_typ_mall;

#ifdef DEBUG_MALLOC
 printf("%d %d %d %d %d %d %d %d %d %d\n",
    ngrp_21_mall,ngrp_33_mall,ngrp_43_mall,ngrp_23_mall,ngrp_46_mall,
    ngrp_typ_21_mall,ngrp_typ_33_mall,ngrp_typ_43_mall,ngrp_typ_23_mall,
    ngrp_typ_46_mall);

 printf("%d %d %d %d %d %d %d %d %d\n",nbond_pow_mall,nbond_typ_pow_mall,
    nbond_con_mall,nbond_typ_con_mall,nbend_pow_mall,nbend_typ_pow_mall,
    nbend_con_mall,nbend_typ_con_mall,nbend_bnd_mall);

 printf("%d %d %d %d %d %d %d \n",nbend_bnd_typ_mall,ntors_pow_mall,
    ntors_typ_pow_mall,ntors_con_mall,ntors_typ_con_mall,nonfo_mall,
    nonfo_typ_mall);
    scanf("%d",&iii);   
#endif

/*========================================================================*/
/*     bonded->grp_bond_con.grp_bond_con_list                             */
/*========================================================================*/
  bonded->grp_bond_con.j1_21 = (int *)cmalloc(ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_21 = (int *)cmalloc(ngrp_21_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j1_23 = (int *)cmalloc(ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j2_23 = (int *)cmalloc(ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_23 = (int *)cmalloc(ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j1_33 = (int *)cmalloc(ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_33 = (int *)cmalloc(ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j3_33 = (int *)cmalloc(ngrp_33_mall*sizeof(int))-1;  

  bonded->grp_bond_watts.j1_33 = (int *)cmalloc(ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j2_33 = (int *)cmalloc(ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j3_33 = (int *)cmalloc(ngrp_watt_33_mall*sizeof(int))-1;

  bonded->grp_bond_con.j1_43   = (int *)cmalloc(ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_43   = (int *)cmalloc(ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_43   = (int *)cmalloc(ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j4_43   = (int *)cmalloc(ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j1_46   = (int *)cmalloc(ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j2_46   = (int *)cmalloc(ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_46   = (int *)cmalloc(ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j4_46   = (int *)cmalloc(ngrp_46_mall*sizeof(int))-1;  

  bonded->grp_bond_con.jtyp_33 = (int *)cmalloc(ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.jtyp_33 = 
                                 (int *)cmalloc(ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_21 = (int *)cmalloc(ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_43 = (int *)cmalloc(ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_23 = (int *)cmalloc(ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.jtyp_46 = (int *)cmalloc(ngrp_46_mall*sizeof(int))-1;

  bonded->grp_bond_watts.cos_thet0_2   = 
                         (double *)malloc(ngrp_typ_watt_33_mall*sizeof(double))-1;
  bonded->grp_bond_watts.sin_thet0_2   = 
                         (double *)malloc(ngrp_typ_watt_33_mall*sizeof(double))-1;

/*========================================================================*/
/*     bonded->grp_bond_con.grp_bond_con_data                             */
/*========================================================================*/

  bonded->grp_bond_con.al_21   = cmall_mat(1,1,1,ngrp_21_mall);
  bonded->grp_bond_con.eq_21   = cmall_mat(1,1,1,ngrp_typ_21_mall);
  bonded->grp_bond_con.al_23   = cmall_mat(1,2,1,ngrp_23_mall);
  bonded->grp_bond_con.eq_23   = cmall_mat(1,2,1,ngrp_typ_23_mall);
  bonded->grp_bond_con.al_33   = cmall_mat(1,3,1,ngrp_33_mall);
  bonded->grp_bond_con.eq_33   = cmall_mat(1,3,1,ngrp_typ_33_mall);
  bonded->grp_bond_watts.eq_33 = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_con.al_43   = cmall_mat(1,3,1,ngrp_43_mall);
  bonded->grp_bond_con.eq_43   = cmall_mat(1,3,1,ngrp_typ_43_mall);
  bonded->grp_bond_con.al_46   = cmall_mat(1,6,1,ngrp_46_mall);
  bonded->grp_bond_con.eq_46   = cmall_mat(1,6,1,ngrp_typ_46_mall); 

  bonded->grp_bond_watts.c_0_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_1_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_2_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_3_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_4_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_5_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_6_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_0_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_1_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_2_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_3_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_4_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_5_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_6_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);

/*========================================================================*/
/*     bonded->bond.bond_list                                             */
/*========================================================================*/

  bonded->bond.j1_pow   = (int *)cmalloc(nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j2_pow   = (int *)cmalloc(nbond_pow_mall*sizeof(int))-1;
  bonded->bond.jtyp_pow = (int *)cmalloc(nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j1_con   = (int *)cmalloc(nbond_con_mall*sizeof(int))-1;
  bonded->bond.j2_con   = (int *)cmalloc(nbond_con_mall*sizeof(int))-1;
  bonded->bond.jtyp_con = (int *)cmalloc(nbond_con_mall*sizeof(int))-1;

/*========================================================================*/
/*     bonded->bond.bond_data                                             */
/*========================================================================*/
  bonded->bond.al_con = (double *)cmalloc(nbond_con_mall*sizeof(double))-1;
  bonded->bond.eq_pow = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.eq_pow_res = 
                        (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.eq_con = (double *)cmalloc(nbond_typ_con_mall*sizeof(double))-1;

  bonded->bond.c_0  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_1  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_2  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_3  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_4  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_5  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.c_6  = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;

  bonded->bond.dc_0 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_1 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_2 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_3 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_4 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_5 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;
  bonded->bond.dc_6 = (double *)cmalloc(nbond_typ_pow_mall*sizeof(double))-1;

/*========================================================================*/
/*     bonded->bond_free.bond_free_data                                   */
/*========================================================================*/

  if(bonded->bond_free.num > 0){
   bonded->bond_free.hist = 
    (double *)cmalloc(bonded->bond_free.nhist*sizeof(double))-1;
  }/*endif*/

/*========================================================================*/
/*     bonded->bend_free.bend_free_data                                   */
/*========================================================================*/

  if(bonded->bend_free.num > 0){
   bonded->bend_free.hist = 
    (double *)cmalloc(bonded->bend_free.nhist*sizeof(double))-1;
  }/*endif*/

/*========================================================================*/
/*     bonded->tors_free.tors_free_data                                   */
/*========================================================================*/

  if(bonded->tors_free.num == 1){
   bonded->tors_free.hist = 
    (double *)cmalloc(bonded->tors_free.nhist*sizeof(double))-1;
  }/*endif*/

  if(bonded->tors_free.num == 2){
   bonded->tors_free.hist_2d = 
         cmall_mat(1,bonded->tors_free.nhist,1,bonded->tors_free.nhist);
  }/*endif*/

/*========================================================================*/
/*     bonded->bend.bend_list                                             */
/*========================================================================*/

  bonded->bend.j1_pow   = (int *) cmalloc(nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j2_pow   = (int *) cmalloc(nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j3_pow   = (int *) cmalloc(nbend_pow_mall*sizeof(int))-1;
  bonded->bend.jtyp_pow = (int *) cmalloc(nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j1_con   = (int *) cmalloc(nbend_con_mall*sizeof(int))-1;
  bonded->bend.j2_con   = (int *) cmalloc(nbend_con_mall*sizeof(int))-1;
  bonded->bend.j3_con   = (int *) cmalloc(nbend_con_mall*sizeof(int))-1;
  bonded->bend.jtyp_con = (int *) cmalloc(nbend_con_mall*sizeof(int))-1;

/*========================================================================*/
/*     bonded->bend.bend_data                                             */
/*========================================================================*/

  ntyp_pow = nbend_typ_pow_mall;
 
  bonded->bend.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.eq_con = (double *)cmalloc(nbend_typ_con_mall*sizeof(double))-1;
  bonded->bend.al_con = (double *)cmalloc(nbend_con_mall*sizeof(double))-1;

/*========================================================================*/
/*     bonded->bend_bnd.bend_bnd_list                                     */
/*========================================================================*/

  bonded->bend_bnd.j1   = (int *) cmalloc(nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j2   = (int *) cmalloc(nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j3   = (int *) cmalloc(nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.jtyp = (int *) cmalloc(nbend_bnd_mall*sizeof(int))-1;

/*========================================================================*/
/*     bonded->bend_bnd.bend_bnd_data                                     */
/*========================================================================*/

  bonded->bend_bnd.eq_bond = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.eq_bend = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_0    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_1    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_2    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_3    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_4    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_5    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbond_6    = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_0   = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_1   = 
                   (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_2   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_3   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_4   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_5   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbond_6   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;

  bonded->bend_bnd.cbend_0    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_1    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_2    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_3    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_4    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_5    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.cbend_6    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_0    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_1    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_2    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_3    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_4    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_5    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.sbend_6    = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_0   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_1   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_2   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_3   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_4   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_5   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dcbend_6   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_0   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_1   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_2   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_3   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_4   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_5   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;
  bonded->bend_bnd.dsbend_6   = 
                  (double *)cmalloc(nbend_bnd_typ_mall*sizeof(double))-1;

/*========================================================================*/
/*     bonded->tors.tors_list                                             */
/*========================================================================*/

  bonded->tors.j1_pow   = (int *) cmalloc(ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j2_pow   = (int *) cmalloc(ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j3_pow   = (int *) cmalloc(ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j4_pow   = (int *) cmalloc(ntors_pow_mall*sizeof(int))-1;
  bonded->tors.jtyp_pow = (int *) cmalloc(ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j1_con   = (int *) cmalloc(ntors_con_mall*sizeof(int))-1;
  bonded->tors.j2_con   = (int *) cmalloc(ntors_con_mall*sizeof(int))-1;
  bonded->tors.j3_con   = (int *) cmalloc(ntors_con_mall*sizeof(int))-1;
  bonded->tors.j4_con   = (int *) cmalloc(ntors_con_mall*sizeof(int))-1;
  bonded->tors.jtyp_con = (int *) cmalloc(ntors_con_mall*sizeof(int))-1;

/*========================================================================*/
/*     bonded->tors.tors_data                                             */
/*========================================================================*/

  bonded->tors.al_con = (double *)cmalloc(ntors_con_mall*sizeof(double))-1;

  ntyp_pow = ntors_typ_pow_mall;

  bonded->tors.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.eq_con = (double *)cmalloc(ntyp_pow*sizeof(double))-1;

/*========================================================================*/
/*     bonded->onfo.onfo_list                                             */
/*========================================================================*/

  bonded->onfo.j1     =  (int *)   cmalloc(nonfo_mall*sizeof(int))-1;
  bonded->onfo.j2     =  (int *)   cmalloc(nonfo_mall*sizeof(int))-1;
  bonded->onfo.jtyp   =  (int *)   cmalloc(nonfo_mall*sizeof(int))-1;

/*========================================================================*/
/*     bonded->onfo.onfo_data                                             */
/*========================================================================*/

  bonded->onfo.feps   =  (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.s6     =  (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.sc     =  (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;

/*========================================================================*/
/*     null_inter_parse->null_inter_parse_list                            */
/*========================================================================*/

  null_inter_parse->jbond1_nul= (int *)
                cmalloc(null_inter_parse->nbond_nul*sizeof(int))-1;
  null_inter_parse->jbond2_nul= (int *)
                cmalloc(null_inter_parse->nbond_nul*sizeof(int))-1;
  null_inter_parse->jbend1_nul= (int *)  
                cmalloc(null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend2_nul= (int *)  
                cmalloc(null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend3_nul= (int *)  
                cmalloc(null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jtors1_nul= (int *) 
                cmalloc(null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors2_nul= (int *) 
                cmalloc(null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors3_nul= (int *) 
                cmalloc(null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors4_nul= (int *) 
                cmalloc(null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jonfo1_nul= (int *) 
                cmalloc(null_inter_parse->nonfo_nul*sizeof(int))-1;
  null_inter_parse->jonfo2_nul  = (int *) 
                cmalloc(null_inter_parse->nonfo_nul*sizeof(int))-1;

/*========================================================================*/
/*     bonded->rbar_sig_free.rbar_sig_free_data                           */
/*========================================================================*/

  if(bonded->rbar_sig_free.nfree>0){
    bonded->rbar_sig_free.hist = cmall_mat(1,bonded->rbar_sig_free.nhist_bar,
                                           1,bonded->rbar_sig_free.nhist_sig);
    bonded->rbar_sig_free.hist_rn = 
                                 cmall_mat(1,bonded->rbar_sig_free.nfree, 
                                           1,bonded->rbar_sig_free.nhist_bar);
    bonded->rbar_sig_free.j1  = (int *)
                           cmalloc(bonded->rbar_sig_free.nfree*sizeof(int))-1;
    bonded->rbar_sig_free.j2  = (int *)
                           cmalloc(bonded->rbar_sig_free.nfree*sizeof(int))-1;
  }/*endif*/
 
/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/








