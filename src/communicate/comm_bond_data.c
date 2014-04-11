/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_bond_data.c                           */
/*                                                                          */
/* Subprogram contains MPI utils and communication routines for interface   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_bond_data(BONDED *bonded,NULL_INTER_PARSE *null_inter_parse,
                                MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
   {/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int iii;

  comm_grp_bond_con_data(&(bonded->grp_bond_con),world);     Barrier(world);
  comm_grp_bond_con_list(&(bonded->grp_bond_con),world);     Barrier(world);

  comm_grp_bond_watts_data(&(bonded->grp_bond_watts),world); Barrier(world);
  comm_grp_bond_watts_list(&(bonded->grp_bond_watts),world); Barrier(world);

  comm_bond_data(&(bonded->bond),world);                     Barrier(world);
  comm_bond_list(&(bonded->bond),world);                     Barrier(world);
  comm_bond_free_data(&(bonded->bond_free),world);           Barrier(world);

  comm_bend_data(&(bonded->bend),world);                     Barrier(world);
  comm_bend_list(&(bonded->bend),world);                     Barrier(world);
  comm_bend_free_data(&(bonded->bend_free),world);           Barrier(world);

  comm_bend_bnd_data(&(bonded->bend_bnd),world);             Barrier(world);
  comm_bend_bnd_list(&(bonded->bend_bnd),world);             Barrier(world);

  comm_tors_data(&(bonded->tors),world);                     Barrier(world);
  comm_tors_list(&(bonded->tors),world);                     Barrier(world);
  comm_tors_free_data(&(bonded->tors_free),world);           Barrier(world);

  comm_onfo_data(&(bonded->onfo),world);                     Barrier(world);
  comm_onfo_list(&(bonded->onfo),world);                     Barrier(world);

  comm_ecor_data(&(bonded->ecor),world);                     Barrier(world);
  comm_null_inter_parse_list(null_inter_parse,world);        Barrier(world);
  comm_rbar_sig_free_data(&(bonded->rbar_sig_free),world);   Barrier(world);

/*------------------------------------------------------------------------*/
     }/*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_con_data(GRP_BOND_CON *grp_bond_con,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ngrp_21_mall     = grp_bond_con->ngrp_21_mall; 
  int ngrp_33_mall     = grp_bond_con->ngrp_33_mall;
  int ngrp_43_mall     = grp_bond_con->ngrp_43_mall;
  int ngrp_23_mall     = grp_bond_con->ngrp_23_mall;
  int ngrp_46_mall     = grp_bond_con->ngrp_46_mall; 
  int ngrp_typ_21_mall = grp_bond_con->ngrp_typ_21_mall; 
  int ngrp_typ_33_mall = grp_bond_con->ngrp_typ_33_mall;
  int ngrp_typ_43_mall = grp_bond_con->ngrp_typ_43_mall;
  int ngrp_typ_23_mall = grp_bond_con->ngrp_typ_23_mall;
  int ngrp_typ_46_mall = grp_bond_con->ngrp_typ_46_mall;

/*------------------------------------------------------------*/

  Barrier(world);
  Bcast(&(grp_bond_con->tol),1,MPI_DOUBLE,0,world);

 if(ngrp_typ_21_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->eq_21[1][1]),ngrp_typ_21_mall,MPI_DOUBLE,0,world);
 }/*endif*/

 if(ngrp_typ_23_mall!=0){
  ngrp_typ_23_mall *=2;
  Barrier(world);
  Bcast(&(grp_bond_con->eq_23[1][1]),ngrp_typ_23_mall,MPI_DOUBLE,0,world);
 }/*endif*/

 if(ngrp_typ_33_mall!=0){
  ngrp_typ_33_mall *=3;
  Barrier(world);
  Bcast(&(grp_bond_con->eq_33[1][1]),ngrp_typ_33_mall,MPI_DOUBLE,0,world);
 }/*endif*/

 if(ngrp_typ_43_mall!=0){
  ngrp_typ_43_mall *=3;
  Barrier(world);
  Bcast(&(grp_bond_con->eq_43[1][1]),ngrp_typ_43_mall,MPI_DOUBLE,0,world);
 }/*endif*/

 if(ngrp_typ_46_mall!=0){
  ngrp_typ_46_mall *=6;
  Barrier(world);
  Bcast(&(grp_bond_con->eq_46[1][1]),ngrp_typ_46_mall,MPI_DOUBLE,0,world);
 }/*endif*/

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_watts_data(GRP_BOND_WATTS *grp_bond_watts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ngrp_33_mall      =  grp_bond_watts->ngrp_33_mall;
  int ngrp_typ_33_mall  =  grp_bond_watts->ngrp_typ_33_mall;
  int ngrp3_typ_33_mall; 
      ngrp3_typ_33_mall = 3*ngrp_typ_33_mall;

 if(ngrp_typ_33_mall!=0){

  Barrier(world);
  Bcast(&(grp_bond_watts->eq_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->cos_thet0_2[1]),ngrp_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->sin_thet0_2[1]),ngrp_typ_33_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(grp_bond_watts->c_0_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_1_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_2_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_3_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_4_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_5_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->c_6_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(grp_bond_watts->dc_0_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_1_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_2_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_3_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_4_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_5_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);
  Bcast(&(grp_bond_watts->dc_6_33[1][1]),ngrp3_typ_33_mall,MPI_DOUBLE,0,world);

 }/*endif*/

/*------------------------------------------------------------------------*/
    } /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bond_data(BOND *bond,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nbond_con_mall     = bond->nbond_con_mall;
  int nbond_typ_pow_mall = bond->nbond_typ_pow_mall;
  int nbond_typ_con_mall = bond->nbond_typ_con_mall;

 if(nbond_typ_pow_mall!=0){

  Barrier(world);
  Bcast(&(bond->eq_pow[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->eq_pow_res[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bond->c_0[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_1[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_2[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_3[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_4[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_5[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->c_6[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bond->dc_0[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_1[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_2[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_3[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_4[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_5[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bond->dc_6[1]),nbond_typ_pow_mall,MPI_DOUBLE,0,world);

 }/*endif*/

 if(nbond_typ_con_mall!=0){
  Barrier(world);
  Bcast(&(bond->eq_con[1]),nbond_typ_con_mall,MPI_DOUBLE,0,world);
 }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bond_free_data(BOND_FREE *bond_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(bond_free->fk),1,MPI_DOUBLE,0,world);
  Bcast(&(bond_free->eq),1,MPI_DOUBLE,0,world);
  Bcast(&(bond_free->del),1,MPI_DOUBLE,0,world);
  Bcast(&(bond_free->rmin),1,MPI_DOUBLE,0,world);
  Bcast(&(bond_free->rmax),1,MPI_DOUBLE,0,world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_data(BEND *bend,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nbend_typ_pow_mall = bend->ntyp_pow;
  int nbend_typ_con_mall = bend->ntyp_con;


 if(nbend_typ_pow_mall!=0){

  Barrier(world);
  Bcast(&(bend->eq_pow[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_0[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_1[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_2[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_3[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_4[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_5[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->c_6[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);  
  Bcast(&(bend->s_0[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_1[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_2[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_3[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_4[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_5[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->s_6[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend->dc_0[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_1[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_2[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_3[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_4[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_5[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->dc_6[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend->ds_0[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_1[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_2[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_3[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_4[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_5[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend->ds_6[1]),nbend_typ_pow_mall,MPI_DOUBLE,0,world);
 }/*endif*/

 if(nbend_typ_con_mall!=0){
  Barrier(world);
  Bcast(&(bend->eq_con[1]),nbend_typ_con_mall,MPI_DOUBLE,0,world);
 }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_free_data(BEND_FREE *bend_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(bend_free->fk),1,MPI_DOUBLE,0,world);
  Bcast(&(bend_free->eq),1,MPI_DOUBLE,0,world);
  Bcast(&(bend_free->del),1,MPI_DOUBLE,0,world);


/*------------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_bnd_data(BEND_BND *bend_bnd,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nbend_bnd_typ_mall = bend_bnd->nbend_bnd_typ_mall;

 if(nbend_bnd_typ_mall!=0){

  Barrier(world);
  Bcast(&(bend_bnd->eq_bond[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->eq_bend[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->cbond_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbond_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->dcbond_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbond_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->cbend_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->cbend_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->sbend_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->sbend_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->dcbend_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dcbend_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(bend_bnd->dsbend_0[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_1[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_2[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_3[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_4[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_5[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(bend_bnd->dsbend_6[1]),nbend_bnd_typ_mall,MPI_DOUBLE,0,world);

 }/*endif*/
/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_tors_data(TORS *tors,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ntors_con_mall = tors->ntors_con_mall;
  int ntors_typ_pow_mall = tors->ntors_typ_pow_mall;

 if(ntors_typ_pow_mall!=0){

  Barrier(world);
  Bcast(&(tors->eq_pow[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->eq_con[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(tors->c_0[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_1[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_2[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_3[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_4[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_5[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->c_6[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(tors->s_0[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_1[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_2[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_3[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_4[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_5[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->s_6[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(tors->dc_0[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_1[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_2[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_3[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_4[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_5[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->dc_6[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);

  Barrier(world);
  Bcast(&(tors->ds_0[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_1[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_2[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_3[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_4[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_5[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);
  Bcast(&(tors->ds_6[1]),ntors_typ_pow_mall,MPI_DOUBLE,0,world);

 }/*endif*/

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_tors_free_data(TORS_FREE *tors_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(tors_free->fk),1,MPI_DOUBLE,0,world);    Barrier(world);
  Bcast(&(tors_free->eq[1]),2,MPI_DOUBLE,0,world); Barrier(world);
  Bcast(&(tors_free->del),1,MPI_DOUBLE,0,world);   Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_onfo_data(ONFO *onfo,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nonfo_typ_mall = onfo->nonfo_typ_mall;

 if(nonfo_typ_mall!=0){
  Barrier(world);
  Bcast(&(onfo->feps[1]),nonfo_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(onfo->s6[1]),nonfo_typ_mall,MPI_DOUBLE,0,world);
  Bcast(&(onfo->sc[1]),nonfo_typ_mall,MPI_DOUBLE,0,world);
  Barrier(world);
 }/*endif*/

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ecor_data(ECOR *ecor,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(ecor->alp_ewd),1,MPI_DOUBLE,0,world);

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bond_list(BOND *bond,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nbond_pow_mall = bond->nbond_pow_mall;
  int nbond_con_mall = bond->nbond_con_mall;

 if(nbond_pow_mall>0){
   Barrier(world);
   Bcast(&(bond->j1_pow[1]),nbond_pow_mall,MPI_INT,0,world);
   Bcast(&(bond->j2_pow[1]),nbond_pow_mall,MPI_INT,0,world);
   Bcast(&(bond->jtyp_pow[1]),nbond_pow_mall,MPI_INT,0,world);
 }/*endif*/

 if(nbond_con_mall>0){
   Barrier(world);
   Bcast(&(bond->j1_con[1]),nbond_con_mall,MPI_INT,0,world);
   Bcast(&(bond->j2_con[1]),nbond_con_mall,MPI_INT,0,world);
   Bcast(&(bond->jtyp_con[1]),nbond_con_mall,MPI_INT,0,world);
 }/*endif*/

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_con_list(GRP_BOND_CON *grp_bond_con,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ngrp21_mall = grp_bond_con->ngrp_21_mall;
  int ngrp23_mall = grp_bond_con->ngrp_23_mall;
  int ngrp33_mall = grp_bond_con->ngrp_33_mall;
  int ngrp43_mall = grp_bond_con->ngrp_43_mall;
  int ngrp46_mall = grp_bond_con->ngrp_46_mall;

if(ngrp21_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->j1_21[1]),ngrp21_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j2_21[1]),ngrp21_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->jtyp_21[1]),ngrp21_mall,MPI_INT,0,world);
}/*endif*/

if(ngrp23_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->j1_23[1]),ngrp23_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j2_23[1]),ngrp23_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j3_23[1]),ngrp23_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->jtyp_23[1]),ngrp23_mall,MPI_INT,0,world);
}/*endif*/

if(ngrp33_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->j1_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j2_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j3_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->jtyp_33[1]),ngrp33_mall,MPI_INT,0,world);
}/*endif*/

if(ngrp43_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->j1_43[1]),ngrp43_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j2_43[1]),ngrp43_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j3_43[1]),ngrp43_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j4_43[1]),ngrp43_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->jtyp_43[1]),ngrp43_mall,MPI_INT,0,world);
}/*endif*/

if(ngrp46_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_con->j1_46[1]),ngrp46_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j2_46[1]),ngrp46_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j3_46[1]),ngrp46_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->j4_46[1]),ngrp46_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_con->jtyp_46[1]),ngrp46_mall,MPI_INT,0,world);
}/*endif*/


/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_watts_list(GRP_BOND_WATTS *grp_bond_watts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ngrp33_mall = grp_bond_watts->ngrp_33_mall;

if(ngrp33_mall!=0){
  Barrier(world);
  Bcast(&(grp_bond_watts->j1_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_watts->j2_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_watts->j3_33[1]),ngrp33_mall,MPI_INT,0,world);
  Bcast(&(grp_bond_watts->jtyp_33[1]),ngrp33_mall,MPI_INT,0,world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_list(BEND *bend,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nbend_pow_mall = bend->nbend_pow_mall;
  int nbend_con_mall = bend->nbend_con_mall;

if(nbend_pow_mall!=0){
  Barrier(world);
  Bcast(&(bend->j1_pow[1]),nbend_pow_mall,MPI_INT,0,world);
  Bcast(&(bend->j2_pow[1]),nbend_pow_mall,MPI_INT,0,world);
  Bcast(&(bend->j3_pow[1]),nbend_pow_mall,MPI_INT,0,world);
  Bcast(&(bend->jtyp_pow[1]),nbend_pow_mall,MPI_INT,0,world);
}/*endif*/
if(nbend_con_mall!=0){
  Barrier(world);
  Bcast(&(bend->j1_con[1]),nbend_con_mall,MPI_INT,0,world);
  Bcast(&(bend->j2_con[1]),nbend_con_mall,MPI_INT,0,world);
  Bcast(&(bend->j3_con[1]),nbend_con_mall,MPI_INT,0,world);
  Bcast(&(bend->jtyp_con[1]),nbend_con_mall,MPI_INT,0,world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_bnd_list(BEND_BND *bend_bnd,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nbend_bnd_mall = bend_bnd->nbend_bnd_mall;

if(nbend_bnd_mall!=0){
  Barrier(world);
  Bcast(&(bend_bnd->j1[1]),nbend_bnd_mall,MPI_INT,0,world);
  Bcast(&(bend_bnd->j2[1]),nbend_bnd_mall,MPI_INT,0,world);
  Bcast(&(bend_bnd->j3[1]),nbend_bnd_mall,MPI_INT,0,world);
  Bcast(&(bend_bnd->jtyp[1]),nbend_bnd_mall,MPI_INT,0,world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_tors_list(TORS *tors,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ntors_pow_mall = tors->ntors_pow_mall;
  int ntors_con_mall = tors->ntors_con_mall;

if(ntors_pow_mall!=0){
  Barrier(world);
  Bcast(&(tors->j1_pow[1]),ntors_pow_mall,MPI_INT,0,world);
  Bcast(&(tors->j2_pow[1]),ntors_pow_mall,MPI_INT,0,world);
  Bcast(&(tors->j3_pow[1]),ntors_pow_mall,MPI_INT,0,world);
  Bcast(&(tors->j4_pow[1]),ntors_pow_mall,MPI_INT,0,world);
  Bcast(&(tors->jtyp_pow[1]),ntors_pow_mall,MPI_INT,0,world);
}/*endif*/

if(ntors_con_mall!=0){
  Barrier(world);
  Bcast(&(tors->j1_con[1]),ntors_con_mall,MPI_INT,0,world);
  Bcast(&(tors->j2_con[1]),ntors_con_mall,MPI_INT,0,world);
  Bcast(&(tors->j3_con[1]),ntors_con_mall,MPI_INT,0,world);
  Bcast(&(tors->j4_con[1]),ntors_con_mall,MPI_INT,0,world);
  Bcast(&(tors->jtyp_con[1]),ntors_con_mall,MPI_INT,0,world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_onfo_list(ONFO *onfo,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nonfo_mall = onfo->nonfo_mall;

if(nonfo_mall!=0){
  Barrier(world);
  Bcast(&(onfo->j1[1]),nonfo_mall,MPI_INT,0,world);
  Bcast(&(onfo->j2[1]),nonfo_mall,MPI_INT,0,world);
  Bcast(&(onfo->jtyp[1]),nonfo_mall,MPI_INT,0,world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_null_inter_parse_list(NULL_INTER_PARSE *null_inter_parse,
                                  MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nbond = null_inter_parse->nbond_nul;
  int nbend = null_inter_parse->nbend_nul;
  int ntors = null_inter_parse->ntors_nul;
  int nonfo = null_inter_parse->nonfo_nul;

  Barrier(world);
  if(nbond>0){
    Barrier(world);
    Bcast(&(null_inter_parse->jbond1_nul[1]),nbond,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jbond2_nul[1]),nbond,MPI_INT,0,world);
  }/*endif*/
  if(nbend>0){
    Barrier(world);
    Bcast(&(null_inter_parse->jbend1_nul[1]),nbend,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jbend2_nul[1]),nbend,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jbend3_nul[1]),nbend,MPI_INT,0,world);
  }/*endif*/
  if(ntors>0){
    Barrier(world);
    Bcast(&(null_inter_parse->jtors1_nul[1]),ntors,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jtors2_nul[1]),ntors,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jtors3_nul[1]),ntors,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jtors4_nul[1]),ntors,MPI_INT,0,world);
  }/*endif*/
  if(nonfo>0){
    Barrier(world);
    Bcast(&(null_inter_parse->jonfo1_nul[1]),nonfo+1,MPI_INT,0,world);
    Bcast(&(null_inter_parse->jonfo2_nul[1]),nonfo+1,MPI_INT,0,world);
  }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_rbar_sig_free_data(RBAR_SIG_FREE *rbar_sig_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
  {/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nfree = rbar_sig_free->nfree;

  if(nfree>0){  

    Barrier(world);
    Bcast(&(rbar_sig_free->rnfree)  ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->fk_bar)  ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->fk_sigma),1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->eq_bar)  ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->eq_sigma),1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->del_bar) ,1,MPI_DOUBLE,0,world);


    Barrier(world);
    Bcast(&(rbar_sig_free->del_sig) ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->rmin)    ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->rmax)    ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->smin)    ,1,MPI_DOUBLE,0,world);
    Bcast(&(rbar_sig_free->smax)    ,1,MPI_DOUBLE,0,world);

    Barrier(world);
    Bcast(&(rbar_sig_free->j1[0]),nfree,MPI_INT,0,world);
    Bcast(&(rbar_sig_free->j2[0]),nfree,MPI_INT,0,world);

  }/*endif*/

/*------------------------------------------------------------------------*/
    }/*end routine*/ 
/*==========================================================================*/







