/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: class_par_forc_lists.c                       */
/*                                                                          */
/*                                                                          */
/* This subprogram  constructus lists needed for classical force            */
/* parallelism                                                              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_class_par_forc_ind(CLASS *class,GENERAL_DATA *general_data,
                            BONDED *bonded)

/*=======================================================================*/
/*  Begin Routine */
  {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
  int num21_me,num23_me,num33_me,num43_me,num46_me;
  int num21_max,num23_max,num33_max,num43_max,num46_max;
  int num21_now,num23_now,num33_now,num43_now,num46_now;
  int igrp_off_21;
  int igrp_off_23;
  int igrp_off_33;
  int igrp_off_43;
  int igrp_off_46;
  int irem,ioff_fwd,ioff_bck,ioff_scr;
  int i,iproc,iii;
/*=======================================================================*/
/*          Local pointers                                               */
  int num21 = bonded->grp_bond_con.num_21;
  int num23 = bonded->grp_bond_con.num_23;
  int num33 = bonded->grp_bond_con.num_33;
  int num43 = bonded->grp_bond_con.num_43;
  int num46 = bonded->grp_bond_con.num_46;

  int *j1_21     = bonded->grp_bond_con.j1_21;
  int *j1_23     = bonded->grp_bond_con.j1_23;
  int *j1_33     = bonded->grp_bond_con.j1_33;
  int *j1_43     = bonded->grp_bond_con.j1_43;
  int *j1_46     = bonded->grp_bond_con.j1_46;
  int *j2_21     = bonded->grp_bond_con.j2_21;
  int *j2_23     = bonded->grp_bond_con.j2_23;
  int *j2_33     = bonded->grp_bond_con.j2_33;
  int *j2_43     = bonded->grp_bond_con.j2_43;
  int *j2_46     = bonded->grp_bond_con.j2_46;
  int *j3_23     = bonded->grp_bond_con.j3_23;
  int *j3_33     = bonded->grp_bond_con.j3_33;
  int *j3_43     = bonded->grp_bond_con.j3_43;
  int *j3_46     = bonded->grp_bond_con.j3_46;
  int *j4_43     = bonded->grp_bond_con.j4_43;
  int *j4_46     = bonded->grp_bond_con.j4_46;
  int ncons_gather     = bonded->grp_bond_con.ncons_gather;
  int ncons_gather_tot = bonded->grp_bond_con.ncons_gather_tot;
  int ncons_gather_max = bonded->grp_bond_con.ncons_gather_max;

  int np_forc        = class->class_comm_forc_pkg.num_proc;
  int myid_forc      = class->class_comm_forc_pkg.myid;

  int *ind_gather_fwd; /* assigned after malloc */
  int *ind_gather_bck;
  int *ind_gather_scr;
  
/*=======================================================================*/
/* 0) Counting                                                           */

  num21_me = num21/np_forc;
  irem = num21%np_forc;
  num21_max = (irem>0 ? num21_me + 1 : num21_me);
  if(irem>myid_forc){num21_me++;}
  if(myid_forc<=irem){
    igrp_off_21 = (num21_me + 1)*myid_forc;
  }else{
    igrp_off_21 = (num21_me + 1)*irem + num21_me*(myid_forc - irem);
  }/*endif*/

  num23_me = num23/np_forc;
  irem = num21%np_forc;
  num23_max = (irem>0 ? num23_me + 1 : num23_me);
  if(irem>myid_forc){num23_me++;}
  if(myid_forc<=irem){
    igrp_off_23 = (num23_me + 1)*myid_forc;
  }else{
    igrp_off_23 = (num23_me + 1)*irem + num23_me*(myid_forc - irem);
  }/*endif*/

  num33_me = num33/np_forc;
  irem = num33%np_forc;
  num33_max = (irem>0 ? num33_me + 1 : num33_me);
  if(irem>myid_forc){num33_me++;}
  if(myid_forc<=irem){
    igrp_off_33 = (num33_me + 1)*myid_forc;
  }else{
    igrp_off_33 = (num33_me + 1)*irem + num33_me*(myid_forc - irem);
  }/*endif*/

  num43_me = num43/np_forc;
  irem = num43%np_forc;
  num43_max = (irem>0 ? num43_me + 1 : num43_me);
  if(irem>myid_forc){num43_me++;}
  if(myid_forc<=irem){
    igrp_off_43 = (num43_me + 1)*myid_forc;
  }else{
    igrp_off_43 = (num43_me + 1)*irem + num43_me*(myid_forc - irem);
  }/*endif*/

  num46_me = num46/np_forc;
  irem = num46%np_forc;
  num46_max = (irem>0 ? num46_me + 1 : num46_me);
  if(irem>myid_forc){num46_me++;}
  if(myid_forc<=irem){
    igrp_off_46 = (num46_me + 1)*myid_forc;
  }else{
    igrp_off_46 = (num46_me + 1)*irem + num46_me*(myid_forc - irem);
  }/*endif*/

/*========================================================================*/
/* I) Malloc the indices                                                  */

  ncons_gather_tot = 2*num21 + 3*num23 + 3*num33 + 4*num43 + 4*num46;
  ncons_gather = 2*num21_me + 3*num23_me + 3*num33_me
                   + 4*num43_me + 4*num46_me;
  ncons_gather_max = 2*num21_max + 3*num23_max + 3*num33_max
                   + 4*num43_max + 4*num46_max;

  bonded->grp_bond_con.ind_gather_fwd =
                  (int *)cmalloc(ncons_gather_max*sizeof(int)) - 1;
  bonded->grp_bond_con.ind_gather_bck =
                  (int *)cmalloc(ncons_gather_max*np_forc*sizeof(int)) - 1;
  bonded->grp_bond_con.ind_gather_scr =
                  (int *)cmalloc(ncons_gather_max*np_forc*sizeof(int)) - 1;

  ind_gather_fwd =  bonded->grp_bond_con.ind_gather_fwd;
  ind_gather_bck =  bonded->grp_bond_con.ind_gather_bck;
  ind_gather_scr =  bonded->grp_bond_con.ind_gather_scr;

/*========================================================================*/
/* II) Fill the forward lists                                             */

  ioff_fwd = 0;

  for(i=1;i<=num21_me;i++){
    ind_gather_fwd[2*i-1+ioff_fwd] = j1_21[i+igrp_off_21];
    ind_gather_fwd[2*i+ioff_fwd]   = j2_21[i+igrp_off_21];
  }/*endfor*/
  ioff_fwd += 2*num21_me;

  for(i=1;i<=num23_me;i++){
    ind_gather_fwd[3*i-2+ioff_fwd] = j1_23[i+igrp_off_23];
    ind_gather_fwd[3*i-1+ioff_fwd] = j2_23[i+igrp_off_23];
    ind_gather_fwd[3*i+ioff_fwd]   = j3_23[i+igrp_off_23];
  }/*endfor*/
  ioff_fwd += 3*num23_me;

  for(i=1;i<=num33_me;i++){
    ind_gather_fwd[3*i-2+ioff_fwd] = j1_33[i+igrp_off_33];
    ind_gather_fwd[3*i-1+ioff_fwd] = j2_33[i+igrp_off_33];
    ind_gather_fwd[3*i+ioff_fwd]   = j3_33[i+igrp_off_33];
  }/*endfor*/
  ioff_fwd += 3*num33_me;

  for(i=1;i<=num43_me;i++){
    ind_gather_fwd[4*i-3+ioff_fwd] = j1_43[i+igrp_off_43];
    ind_gather_fwd[4*i-2+ioff_fwd] = j2_43[i+igrp_off_43];
    ind_gather_fwd[4*i-1+ioff_fwd] = j3_43[i+igrp_off_43];
    ind_gather_fwd[4*i+ioff_fwd]   = j4_43[i+igrp_off_43];
  }/*endfor*/
  ioff_fwd += 4*num43_me;

  for(i=1;i<=num46_me;i++){
    ind_gather_fwd[4*i-3+ioff_fwd] = j1_46[i+igrp_off_46];
    ind_gather_fwd[4*i-2+ioff_fwd] = j2_46[i+igrp_off_46];
    ind_gather_fwd[4*i-1+ioff_fwd] = j3_46[i+igrp_off_46];
    ind_gather_fwd[4*i+ioff_fwd]   = j4_46[i+igrp_off_46];
  }/*endfor*/
  ioff_fwd += 4*num46_me;

/*========================================================================*/
/* II) Fill the back lists                                                */

  ioff_bck = 0;
  ioff_scr = 0;

  for(iproc=0;iproc<np_forc;iproc++){

    num21_now = num21/np_forc;
    irem = num21%np_forc;
    if(iproc<=irem){
      igrp_off_21 = (num21_now + 1)*iproc;
    }else{
      igrp_off_21 = (num21_now + 1)*irem + num21_now*(iproc - irem);
    }/*endif*/
    if(irem>iproc){num21_now++;}

    for(i=1;i<=num21_now;i++){
      ind_gather_bck[2*i-1+ioff_bck] = j1_21[i+igrp_off_21];
      ind_gather_bck[2*i+ioff_bck]   = j2_21[i+igrp_off_21];
      ind_gather_scr[2*i-1+ioff_bck] = 2*i-1+ioff_scr;
      ind_gather_scr[2*i+ioff_bck]   = 2*i+ioff_scr;
    }/*endfor*/
    ioff_bck += 2*num21_now;
    ioff_scr += 2*num21_max;

    num23_now = num23/np_forc;
    irem = num23%np_forc;
    if(iproc<=irem){
      igrp_off_23 = (num23_now + 1)*iproc;
    }else{
      igrp_off_23 = (num23_now + 1)*irem + num23_now*(iproc - irem);
    }/*endif*/
    if(irem>iproc){num23_now++;}

    for(i=1;i<=num23_now;i++){
      ind_gather_bck[3*i-2+ioff_bck] = j1_23[i+igrp_off_23];
      ind_gather_bck[3*i-1+ioff_bck] = j2_23[i+igrp_off_23];
      ind_gather_bck[3*i+ioff_bck]   = j3_23[i+igrp_off_23];
      ind_gather_scr[3*i-2+ioff_bck] = 3*i-2+ioff_scr;
      ind_gather_scr[3*i-1+ioff_bck] = 3*i-1+ioff_scr;
      ind_gather_scr[3*i+ioff_bck]   = 3*i+ioff_scr;
    }/*endfor*/
    ioff_bck += 3*num23_now;
    ioff_scr += 3*num23_max;

    num33_now = num33/np_forc;
    irem = num33%np_forc;
    if(iproc<=irem){
      igrp_off_33 = (num33_now + 1)*iproc;
    }else{
      igrp_off_33 = (num33_now + 1)*irem + num33_now*(iproc - irem);
    }/*endif*/
    if(irem>iproc){num33_now++;}

    for(i=1;i<=num33_now;i++){
      ind_gather_bck[3*i-2+ioff_bck] = j1_33[i+igrp_off_33];
      ind_gather_bck[3*i-1+ioff_bck] = j2_33[i+igrp_off_33];
      ind_gather_bck[3*i+ioff_bck]   = j3_33[i+igrp_off_33];
      ind_gather_scr[3*i-2+ioff_bck] = 3*i-2+ioff_scr;
      ind_gather_scr[3*i-1+ioff_bck] = 3*i-1+ioff_scr;
      ind_gather_scr[3*i+ioff_bck]   = 3*i+ioff_scr;
    }/*endfor*/
    ioff_bck += 3*num33_now;
    ioff_scr += 3*num33_max;

    num43_now = num43/np_forc;
    irem = num43%np_forc;
    if(iproc<=irem){
      igrp_off_43 = (num43_now + 1)*iproc;
    }else{
      igrp_off_43 = (num43_now + 1)*irem + num43_now*(iproc - irem);
    }/*endif*/
    if(irem>iproc){num43_now++;}

    for(i=1;i<=num43_now;i++){
      ind_gather_bck[4*i-3+ioff_bck] = j1_43[i+igrp_off_43];
      ind_gather_bck[4*i-2+ioff_bck] = j2_43[i+igrp_off_43];
      ind_gather_bck[4*i-1+ioff_bck] = j3_43[i+igrp_off_43];
      ind_gather_bck[4*i+ioff_bck]   = j4_43[i+igrp_off_43];
      ind_gather_scr[4*i-3+ioff_bck] = 4*i-3+ioff_scr;
      ind_gather_scr[4*i-2+ioff_bck] = 4*i-2+ioff_scr;
      ind_gather_scr[4*i-1+ioff_bck] = 4*i-1+ioff_scr;
      ind_gather_scr[4*i+ioff_bck]   = 4*i+ioff_scr;
    }/*endfor*/
    ioff_bck += 4*num43_now;
    ioff_scr += 4*num43_max;

    num46_now = num46/np_forc;
    irem = num46%np_forc;
    if(iproc<=irem){
      igrp_off_46 = (num46_now + 1)*iproc;
    }else{
      igrp_off_46 = (num46_now + 1)*irem + num46_now*(iproc - irem);
    }/*endif*/
    if(irem>iproc){num46_now++;}

    for(i=1;i<=num46_now;i++){
      ind_gather_bck[4*i-3+ioff_bck] = j1_46[i+igrp_off_46];
      ind_gather_bck[4*i-2+ioff_bck] = j2_46[i+igrp_off_46];
      ind_gather_bck[4*i-1+ioff_bck] = j3_46[i+igrp_off_46];
      ind_gather_bck[4*i+ioff_bck]   = j4_46[i+igrp_off_46];
      ind_gather_scr[4*i-3+ioff_bck] = 4*i-3+ioff_scr;
      ind_gather_scr[4*i-2+ioff_bck] = 4*i-2+ioff_scr;
      ind_gather_scr[4*i-1+ioff_bck] = 4*i-1+ioff_scr;
      ind_gather_scr[4*i+ioff_bck]   = 4*i+ioff_scr;
    }/*endfor*/
    ioff_bck += 4*num46_now;
    ioff_scr += 4*num46_max;

  }/*endfor:np_proc*/

/*========================================================================*/
}/*end routine*/
/*========================================================================*/

