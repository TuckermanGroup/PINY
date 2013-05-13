/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: verlist_create.c                             */
/*                                                                          */
/* This routine generates a verlet neighbor list                            */
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void verlist_create(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    INTRA_SCR *intra_scr,FOR_SCR *for_scr,
                    ATOMMAPS *atommaps,CELL *cell,NBR_LIST *nbr_list,
                    TIMEINFO *timeinfo, INTERACT *interact,
                    int iswitch, int *nver_now,int *nver_res_now,int lst_typ)

/*========================================================================*/
    {/*Begin routine*/
/*========================================================================*/ 
/*            Local variable                                              */

  int intact_now,ifdx,natm1_typ2,iii,nstart,nstart_res;
  int isum,ibrk,jpart,num_brk;
  int i2,i3,ivdx;
  int iitemp,jjtemp,kktemp,mtemp,itemp,nver_now_loc;
  int upper,lower;
  int n2;
  double q1t,q2t,qtemp,qsum,rtemp;

/* Local Pointers */

  int *for_scr_i_index         = for_scr->i_index;
  int *for_scr_j_index         = for_scr->j_index;
  int *for_scr_i_indext        = for_scr->i_indext;
  int *for_scr_j_indext        = for_scr->j_indext;
  int *atommaps_iatm_atm_typ   = atommaps->iatm_atm_typ;
  double *clatoms_x            = clatoms_pos->x;
  double *clatoms_y            = clatoms_pos->y;
  double *clatoms_z            = clatoms_pos->z;
  double *intra_scr_x1         = intra_scr->x1;  
  double *intra_scr_y1         = intra_scr->y1;  
  double *intra_scr_z1         = intra_scr->z1;  
  double *intra_scr_x2         = intra_scr->x2;  
  double *intra_scr_y2         = intra_scr->y2;  
  double *intra_scr_z2         = intra_scr->z2;  
  double *intra_scr_dx12       = intra_scr->dx12;  
  double *intra_scr_dy12       = intra_scr->dy12;  
  double *intra_scr_dz12       = intra_scr->dz12;  
  double *intra_scr_x5         = intra_scr->x5;  
  double *intra_scr_y5         = intra_scr->y5;  
  double *intra_scr_z5         = intra_scr->z5;  
  double *intra_scr_x6         = intra_scr->x6;  
  double *intra_scr_y6         = intra_scr->y6;  
  double *intra_scr_z6         = intra_scr->z6;  
  double *intra_scr_dx56       = intra_scr->dx56;  
  double *intra_scr_dy56       = intra_scr->dy56;  
  double *intra_scr_dz56       = intra_scr->dz56;  
  double *interact_cutskin     = interact->cutskin;
  double *interact_cutskin_res = interact->cutskin_res;
  double *interact_cutskin_root     = interact->cutskin_root;
  double *interact_cutskin_root_res = interact->cutskin_root_res;
  double *xmod                 = clatoms_info->xmod;
  double *ymod                 = clatoms_info->ymod;
  double *zmod                 = clatoms_info->zmod;
  int    *nbr_list_nter        = nbr_list->verlist.nter;
  int    *nbr_list_nter_res    = nbr_list->verlist.nter_res;
  list_int  *nbr_list_jter_res    = nbr_list->verlist.jter_res;
  list_int  *nbr_list_jter        = nbr_list->verlist.jter;
  int    *for_scr_num_brk_i    = for_scr->num_brk_i;
  int    *nbrnch_of_root       = nbr_list->brnch_root.nbrnch_of_root_big;
  int    brnch_root_list_opt   = nbr_list->brnch_root_list_opt; 

  int iperd    = cell->iperd;
  int pi_beads = clatoms_info->pi_beads;
  int natm_typ = atommaps->natm_typ;

  double *intra_scr_cutoff     = intra_scr->cutoff;
  double *intra_scr_cutoff_res = intra_scr->cutoff_res;
  double *intra_scr_r          = intra_scr->r;  
  int iver_fill  = nbr_list->verlist.iver_fill;
  int iver_count = nbr_list->verlist.iver_count;
  int int_res_ter = timeinfo->int_res_ter;
  int ishave_opt,intact_use;
  int icoul;

  double *interact_cutoff         = interact->cutoff;
  double *interact_cutoff_res     = interact->cutoff_res;
  double *interact_rmin_spl       = interact->rmin_spl;
  double *interact_dri_spl        = interact->dri_spl;
  double *interact_dr_spl         = interact->dr_spl;
  double *interact_vcut_coul      = interact->vcut_coul;
  double interact_rheal_res       = interact->rheal_res;
  double interact_rheal_resi;
  double *interact_cv0            = interact->cv0;
  double *interact_cv0_c          = interact->cv0_c;
  double *interact_cdv0           = interact->cdv0;
  double *interact_cdv0_c         = interact->cdv0_c;
  int diele_opt = interact->dielectric_opt;

  double *intra_scr_q2            = intra_scr->q2;
  double *intra_scr_del_r         = intra_scr->del_r;
  double *intra_scr_swit_hard     = intra_scr->dz43;
  double *intra_scr_swit          = intra_scr->swit;
  double *intra_scr_vcut_coul     = intra_scr->vcut_coul;
  double *intra_scr_vpot          = intra_scr->vpot;
  double *intra_scr_spl_tmp       = intra_scr->spl_tmp;
  double *intra_scr_fx1           = intra_scr->fx1;
  double *intra_scr_fy1           = intra_scr->fy1;
  double *intra_scr_fz1           = intra_scr->fz1;
  double *intra_scr_fx2           = intra_scr->fx2;
  double *intra_scr_fy2           = intra_scr->fy2;
  double *intra_scr_fz2           = intra_scr->fz2;
  double *intra_scr_c_0           = intra_scr->c_0;
  double * intra_scr_p11          = intra_scr->p11;
  double * intra_scr_p22          = intra_scr->p22;
  double * intra_scr_p33          = intra_scr->p33;
  double * intra_scr_p12          = intra_scr->p12;
  double * intra_scr_p21          = intra_scr->p21;
  double * intra_scr_p13          = intra_scr->p13;
  double * intra_scr_p31          = intra_scr->p31;
  double * intra_scr_p23          = intra_scr->p23;
  double * intra_scr_p32          = intra_scr->p32;
 
  int nchrg                       = clatoms_info->nchrg;
  double *clatoms_fxt             = clatoms_pos->fxt;
  double *clatoms_fyt             = clatoms_pos->fyt;
  double *clatoms_fzt             = clatoms_pos->fzt;
  double *clatoms_fx              = clatoms_pos->fx;
  double *clatoms_fy              = clatoms_pos->fy;
  double *clatoms_fz              = clatoms_pos->fz;
  double *clatoms_q               = clatoms_info->q;
  int interact_nsplin;
  int interact_nsplin_m2;

  double wght_res;
  double wght_full;
  double wght_lnk=1.0;
  double wght_ter;
  int energy_ctrl_iget_res_inter;
  int energy_ctrl_iget_full_inter;
  int energy_ctrl_int_res_ter;
  int energy_ctrl_isep_vvdw;

  interact_nsplin         = interact->nsplin;
  interact_nsplin_m2      = interact->nsplin-2;
  interact_rheal_resi     = 1.0/interact->rheal_res;

/*========================================================================*/
/* I) Error check                                                         */

  intact_now = for_scr->intact; 
  num_brk = for_scr->num_brk;

  isum    = 0;
  for (ibrk = 1; ibrk <= num_brk; ++ibrk) {
    isum += for_scr_num_brk_i[ibrk];
  }/*endfor*/
  if (isum != for_scr->intact) {
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error:\n ");
    printf("number of interactions inconsistent with\n");
    printf("the breakpoints: in verlist_create\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
  
  if(lst_typ==3 && iver_count==1 && iver_fill==1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal error: No ver_count plus ver_fill with lnk lst \n");    
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Gather the displacements/interaction types in this set     */

  vercreate_posdata_fetch(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                        for_scr_i_index,
                        for_scr_j_index,for_scr_i_indext,
                        intra_scr_dx56,intra_scr_dy56,intra_scr_dz56,
                        clatoms_x,clatoms_y,clatoms_z,atommaps_iatm_atm_typ,
                        xmod,ymod,zmod,
                        num_brk,for_scr_num_brk_i,iperd,pi_beads,lst_typ,
                        intact_now,natm_typ,cell);
     
/*========================================================================*/
/*  III) Calculate distances, shave the interactions, count the list       */

  if(iver_fill == 0 && iver_count == 1 ){
   if(brnch_root_list_opt==0){
    vercreate_list_count(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                        intact_now,for_scr_i_index,for_scr_j_index,
                        for_scr_i_indext,
                        interact_cutskin,interact_cutskin_res,
                        intra_scr_r,intra_scr_cutoff_res,
                        intra_scr_cutoff,
                        nbr_list_nter,nbr_list_jter,nver_now,
                        nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                        int_res_ter);
  }/*endif*/ else{  
    vercreate_list_count_root(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                        intact_now,for_scr_i_index,for_scr_j_index,
                        for_scr_i_indext,
                        interact_cutskin_root,interact_cutskin_root_res,
                        intra_scr_r,intra_scr_cutoff_res,
                        intra_scr_cutoff,
                        nbr_list_nter,nbr_list_jter,nver_now,
                        nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                        int_res_ter,nbrnch_of_root);

   }/*endelse*/
  }/*endif*/

/*========================================================================*/
/*  IV) Calculate distances, shave the interactions, fill the list        */

  if(iver_fill == 1 && iver_count == 0){
   if(brnch_root_list_opt==0){
   vercreate_list_fill(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intact_now,for_scr_i_index,for_scr_j_index,
                       for_scr_i_indext,
                       interact_cutskin,interact_cutskin_res,
                       intra_scr_r,intra_scr_cutoff_res,
                       intra_scr_cutoff,
                       nbr_list_nter,nbr_list_jter,nver_now,
                       nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                       int_res_ter,num_brk,for_scr_num_brk_i);

  }/*endif*/ else{  
   vercreate_list_fill_root(intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intact_now,for_scr_i_index,for_scr_j_index,
                       for_scr_i_indext,
                       interact_cutskin_root,interact_cutskin_root_res,
                       intra_scr_r,intra_scr_cutoff_res,
                       intra_scr_cutoff,
                       nbr_list_nter,nbr_list_jter,nver_now,
                       nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                       int_res_ter,num_brk,for_scr_num_brk_i,nbrnch_of_root);
   }/*endelse*/
  }/*endif*/

/*========================================================================*/
/*  V) Calculate distances, shave the interactions, fill/count the list   */

  if(iver_fill == 1 && iver_count ==1){
   if(brnch_root_list_opt==0){
   vercreate_list_count_fill(
                       num_brk,for_scr_num_brk_i,lst_typ,
                       intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intact_now,for_scr_i_index,for_scr_j_index,
                       for_scr_i_indext,
                       interact_cutskin,interact_cutskin_res,
                       intra_scr_r,intra_scr_cutoff_res,
                       intra_scr_cutoff,
                       nbr_list_nter,nbr_list_jter,nver_now,
                       nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                       int_res_ter);

   }/*endif*/ else{  
   vercreate_list_count_fill_root(
                       num_brk,for_scr_num_brk_i,lst_typ,
                       intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,
                       intact_now,for_scr_i_index,for_scr_j_index,
                       for_scr_i_indext,
                       interact_cutskin_root,interact_cutskin_root_res,
                       intra_scr_r,intra_scr_cutoff_res,
                       intra_scr_cutoff,
                       nbr_list_nter,nbr_list_jter,nver_now,
                       nbr_list_nter_res,nbr_list_jter_res,nver_res_now,
                       int_res_ter,nbrnch_of_root);
   }/*endelse*/
  }/*endif*/

/*========================================================================*/
/*========================================================================*/


/*--------------------------------------------------------------------------*/
   }/*end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vercreate_posdata_fetch(double *dx12,double *dy12,double *dz12,
                        int *i_index,int *j_index,int *jatm_typ,
                        double *dx56,double *dy56,double *dz56,
                        double *clatoms_x,double *clatoms_y,double *clatoms_z,
                        int *clatoms_atm_typ,
                        double *xmod,double *ymod,double *zmod,
                        int num_brk,int *num_brk_i,int iperd,int pi_beads,
                        int lst_typ,int intact_now,int natm_typ,
                        CELL *cell)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i,ibrk,lower,iii;
  int jpart,ktemp,icoul;
  int upper,n2;
  double qsum,qtemp;
  double x1t,y1t,z1t;
  double x2t,y2t,z2t;
  double x5t,y5t,z5t;
  double x6t,y6t,z6t;
  double q1t,q2t;
  int natm1_typ2,i2,i3,jjtemp,iitemp,kktemp,itemp,jtemp;
  int iatm_typt,jatm_typt;

/*=======================================================================*/
/* I) Fetch the positions, get displacements and interaction number      */

   natm1_typ2 = (natm_typ+1) * 2;
  /*---------------------------------------------------------------*/
  /* For a link list */
if(lst_typ == 3){
    for (jpart = 1; jpart <= intact_now; jpart++) {
      itemp          =  i_index[jpart];
      jtemp          =  j_index[jpart];
      j_index[jpart] = MIN(itemp,jtemp);
      i_index[jpart] = itemp+jtemp-j_index[jpart];  /*MAX(itemp,jtemp)*/
      x1t         = clatoms_x[i_index[jpart]];
      y1t         = clatoms_y[i_index[jpart]];
      z1t         = clatoms_z[i_index[jpart]];
      x2t         = clatoms_x[j_index[jpart]];
      y2t         = clatoms_y[j_index[jpart]];
      z2t         = clatoms_z[j_index[jpart]];
      dx12[jpart] = x1t-x2t;
      dy12[jpart] = y1t-y2t;
      dz12[jpart] = z1t-z2t;
      jatm_typt   = clatoms_atm_typ[j_index[jpart]];
      iatm_typt   = clatoms_atm_typ[i_index[jpart]];
      i2          = iatm_typt;
      i3          = jatm_typt;
      jjtemp      = MIN(i2,i3); 
      iitemp      = i2+i3-jjtemp;  /*MAX(i2,i3)*/
      kktemp      = ((jjtemp - 1) * (natm1_typ2 - jjtemp))/2; 
      jatm_typ[jpart] = (kktemp + iitemp - jjtemp + 1); 
    }/*endfor*/
    if(iperd>0&&pi_beads>1){     
      for (jpart = 1; jpart <= intact_now; jpart++) {
        x5t         = xmod[i_index[jpart]];
        y5t         = ymod[i_index[jpart]];
        z5t         = zmod[i_index[jpart]];
        x6t         = xmod[j_index[jpart]];
        y6t         = ymod[j_index[jpart]];
        z6t         = zmod[j_index[jpart]];
        dx56[jpart] = x5t-x6t;
        dy56[jpart] = y5t-y6t;
        dz56[jpart] = z5t-z6t;
      }/*endfor*/
    }/*endif*/
}/*endif*/

  /*---------------------------------------------------------------*/
  /* For a nolist or verlet list */

if(lst_typ != 3){
   lower = 1;
   for (ibrk = 1; ibrk <= num_brk; ibrk++) {
     upper = num_brk_i[ibrk] + lower - 1;
     ktemp = i_index[lower];
#ifndef NO_PRAGMA
#pragma ivdep
#endif
     for(jpart=lower; jpart<=upper;jpart++){
      x1t         = clatoms_x[ktemp];
      y1t         = clatoms_y[ktemp];
      z1t         = clatoms_z[ktemp];
      x2t         = clatoms_x[j_index[jpart]];
      y2t         = clatoms_y[j_index[jpart]];
      z2t         = clatoms_z[j_index[jpart]];
      dx12[jpart] = x1t-x2t;
      dy12[jpart] = y1t-y2t;
      dz12[jpart] = z1t-z2t;
      jatm_typt   = clatoms_atm_typ[j_index[jpart]];
      iatm_typt   = clatoms_atm_typ[ktemp];
      i2          = iatm_typt;
      i3          = jatm_typt;
      jjtemp      = MIN(i2,i3); 
      iitemp      = i2+i3-jjtemp;  /*MAX(i2,i3)*/
      kktemp      = ((jjtemp - 1) * (natm1_typ2 - jjtemp))/2; 
      jatm_typ[jpart] = (kktemp + iitemp - jjtemp + 1); 
     }/*endfor*/
     if(iperd>0&&pi_beads>1){     
       for(jpart=lower; jpart<=upper;jpart++){
        x5t         = xmod[ktemp];
        y5t         = ymod[ktemp];
        z5t         = zmod[ktemp];
        x6t         = xmod[j_index[jpart]];
        y6t         = ymod[j_index[jpart]];
        z6t         = zmod[j_index[jpart]];
        dx56[jpart] = x5t-x6t;
        dy56[jpart] = y5t-y6t;
        dz56[jpart] = z5t-z6t;
       }/*endfor*/
     }/*endif*/
     lower += num_brk_i[ibrk];
    }/*endfor*/
}/*endif lst_typ */ 

/*=======================================================================*/
/* II) Periodically image the distances                                  */

  if(iperd>0){
#define DEBUG_OFF
#ifdef DEBUG
    period(intact_now,dx12,dy12,dz12,cell);
#endif
#ifdef DEBUG_OFF
   if(pi_beads==1){     
    period(intact_now,dx12,dy12,dz12,cell);
   }else{
    period_pimd(intact_now,dx12,dy12,dz12,dx56,dy56,dz56,cell);
   }/*endif*/
#endif
  }/*endif*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count and fill the list at the same time                                 */
/*==========================================================================*/

void vercreate_list_count_fill(
              int num_brk,int *for_scr_num_brk_i,int lst_typ,
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart,ifdx,ibrk,upper,lower,nstart,nstart_res;
  int iii;
/*==========================================================================*/
/* 0) Get distances and cutoffs                                             */

     if(int_res_ter==0){
       for (jpart = 1; jpart <= intact_now; ++jpart) {
          intra_scr_cutoff[jpart]     = 
                            interact_cutskin[(for_scr_i_indext[jpart])];
       }/*endfor*/
     }/*endif*/

     if(int_res_ter==1){
       for (jpart = 1; jpart <= intact_now; ++jpart) {
          intra_scr_cutoff[jpart]     = 
                            interact_cutskin[(for_scr_i_indext[jpart])];
          intra_scr_cutoff_res[jpart]     = 
                            interact_cutskin_res[(for_scr_i_indext[jpart])];
       }/*endfor*/
     }/*endif*/

     for (jpart = 1; jpart <= intact_now; ++jpart) {
      intra_scr_r[jpart] = dx12[jpart]*dx12[jpart]
                         + dy12[jpart]*dy12[jpart]
                         + dz12[jpart]*dz12[jpart];
     }/*endfor*/

/*==========================================================================*/
/* I) Inter Respa off*/

     if(int_res_ter==0){
  /*----------------------------------------------------------------------*/
  /* From lnk lst data                                                    */
       if(lst_typ == 3){
        for (jpart = 1; jpart <= intact_now; ++jpart) {
          if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
             ifdx =  for_scr_i_index[jpart];
             nbr_list_nter[ifdx]++;
             (*nver_now)++;
             nbr_list_jter[(*nver_now)] = for_scr_j_index[jpart];
          }/*endif*/
        }/*endfor*/
       }else{
  /*----------------------------------------------------------------------*/
  /* From no_lst data                                                    */
        lower = 1;
        for (ibrk = 1; ibrk <= num_brk; ibrk++) {
          upper = for_scr_num_brk_i[ibrk] + lower - 1;
          ifdx = for_scr_i_index[lower];
          nstart = (*nver_now);
          for(jpart=lower; jpart<=upper;jpart++){
            if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
             (*nver_now)++;
             nbr_list_jter[(*nver_now)] = for_scr_j_index[jpart];
            }/*endif*/
          }/*endfor*/
          nbr_list_nter[ifdx] += (*nver_now)-nstart;
          lower += for_scr_num_brk_i[ibrk];
        }/*endfor*/
       }/*endif*/
     }/*endif: RESPA off*/

/*==========================================================================*/
/* ii) Inter Respa on*/

     if(int_res_ter==1){
  /*----------------------------------------------------------------------*/
  /* From lnk lst data                                                    */
       if(lst_typ == 3){
        for (jpart = 1; jpart <= intact_now; ++jpart) {
         if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
             ifdx =  for_scr_i_index[jpart];
             nbr_list_nter[ifdx]++;
             (*nver_now)++;
             nbr_list_jter[(*nver_now)] = for_scr_j_index[jpart];

             if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
              (*nver_res_now)++;
              nbr_list_nter_res[ifdx]++;
              nbr_list_jter_res[(*nver_res_now)] = for_scr_j_index[jpart];
             }/*endif: add to RESPA list*/
         }/*endif: add to big list*/
       }/*endfor: jpart interaction*/
      }else{
  /*----------------------------------------------------------------------*/
  /* From no lst data                                                    */
        lower = 1;
        for (ibrk = 1; ibrk <= num_brk; ibrk++) {
          upper = for_scr_num_brk_i[ibrk] + lower - 1;
          ifdx = for_scr_i_index[lower];
          nstart     = (*nver_now);
          nstart_res = (*nver_res_now);
          for(jpart=lower; jpart<=upper;jpart++){
            if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
             (*nver_now)++;
             nbr_list_jter[(*nver_now)] = for_scr_j_index[jpart];
             if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
              (*nver_res_now)++;
              nbr_list_jter_res[(*nver_res_now)] = for_scr_j_index[jpart];
             }/*endif: add to RESPA list*/
            }/*endif: add to big list*/
          }/*endfor*/
          nbr_list_nter[ifdx]     += ((*nver_now)-nstart);
          nbr_list_nter_res[ifdx] += ((*nver_res_now)-nstart_res);
          lower += for_scr_num_brk_i[ibrk];
        }/*endfor*/
       }/*endif*/
     }/*endif:RESPA on*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Count and fill the list at the same time                                 */
/*==========================================================================*/

void vercreate_list_count_fill_root(
              int num_brk,int *for_scr_num_brk_i,int lst_typ,
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter, int *nbrnch_of_root)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart,ifdx,ibrk,upper,lower,nstart,nstart_res;
  int iii,jtemp,ntemp;
/*==========================================================================*/
/* 0) Get distances and cutoffs                                             */

  if(int_res_ter==0){
     for (jpart = 1; jpart <= intact_now; ++jpart) {
       intra_scr_cutoff[jpart]     = 
                          interact_cutskin[(for_scr_i_indext[jpart])];
     }/*endfor*/
  }/*endif*/

  if(int_res_ter==1){
     for (jpart = 1; jpart <= intact_now; ++jpart) {
        intra_scr_cutoff[jpart]     = 
                           interact_cutskin[(for_scr_i_indext[jpart])];
        intra_scr_cutoff_res[jpart]     = 
                           interact_cutskin_res[(for_scr_i_indext[jpart])];
     }/*endfor*/
  }/*endif*/

  for (jpart = 1; jpart <= intact_now; ++jpart) {
    intra_scr_r[jpart] = dx12[jpart]*dx12[jpart]
                       + dy12[jpart]*dy12[jpart]
                       + dz12[jpart]*dz12[jpart];
  }/*endfor*/

/*==========================================================================*/
/* I) Inter Respa off*/

  if(int_res_ter==0){

    lower = 1;
    for (ibrk = 1; ibrk <= num_brk; ibrk++) {
       upper = for_scr_num_brk_i[ibrk] + lower - 1;
       ifdx  = for_scr_i_index[lower];
       nstart = (*nver_now);
       for(jpart=lower; jpart<=upper;jpart++){
         if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
             nbr_list_jter[((*nver_now)+1)] = for_scr_j_index[jpart];
             (*nver_now)+= nbrnch_of_root[for_scr_j_index[jpart]];
         }/*endif*/
       }/*endfor*/
       nbr_list_nter[ifdx] += ((*nver_now)-nstart);
       lower += for_scr_num_brk_i[ibrk];
    }/*endfor*/

  }/*endif: RESPA off*/


/*==========================================================================*/
/* ii) Inter Respa on*/

  if(int_res_ter==1){

    lower = 1;
    for (ibrk = 1; ibrk <= num_brk; ibrk++) {
       upper = for_scr_num_brk_i[ibrk] + lower - 1;
       ifdx = for_scr_i_index[lower];
       nstart     = (*nver_now);
       nstart_res = (*nver_res_now);
       for(jpart=lower; jpart<=upper;jpart++){
         if(intra_scr_r[jpart] < intra_scr_cutoff[jpart]){
           nbr_list_jter[((*nver_now)+1)] = for_scr_j_index[jpart];
           ntemp = nbrnch_of_root[for_scr_j_index[jpart]];
           (*nver_now)+=ntemp;
           if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
              nbr_list_jter_res[((*nver_res_now)+1)]=for_scr_j_index[jpart];
              (*nver_res_now)+=ntemp;
           }/*endif: add to RESPA list*/
         }/*endif: add to big list*/
       }/*endfor*/
       nbr_list_nter[ifdx]     += ((*nver_now)-nstart);
       nbr_list_nter_res[ifdx] +=((*nver_res_now)-nstart_res);
       lower += for_scr_num_brk_i[ibrk];
    }/*endfor*/

  }/*endif:RESPA on*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Shave the interactions down                                              */
/*==========================================================================*/

void vercreate_list_count(
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   int jpart,ic,intact_use,ifdx,ivdx,iii;
   double cutoff,r;

/*=======================================================================*/
/* II) Shrink and count the list no respa                               */

  if(int_res_ter==0){

    ic = 0;
    for (jpart = 1; jpart <= intact_now; ++jpart) {
      cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
      r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
              + dz12[jpart]*dz12[jpart];
      if(r < cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
      }/*endif*/
    }/*endfor*/
    intact_use = ic;

    for (jpart = 1; jpart <= intact_use; ++jpart) {
      ifdx = for_scr_i_index[jpart];
      nbr_list_nter[ifdx]++;
    }/*endfor:jpart interaction*/
    (*nver_now)+=intact_use;

  }/*endif:RESPA off*/

/*=======================================================================*/
/* II) Shrink and count the list respa                                   */

  if(int_res_ter==1){

    ic = 0;
    for (jpart = 1; jpart <= intact_now; ++jpart) {
      cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
      r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
              + dz12[jpart]*dz12[jpart];
      if(r<cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
        intra_scr_cutoff_res[ic] = 
                            interact_cutskin_res[(for_scr_i_indext[jpart])];
        intra_scr_r[ic]     = r;
      }/*endif*/
    }/*endfor*/
    intact_use = ic;

    for (jpart = 1; jpart <= intact_use; ++jpart) {
      ifdx = for_scr_i_index[jpart];
      nbr_list_nter[ifdx]++;
      if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
         (*nver_res_now)++;
         nbr_list_nter_res[ifdx]++;
      }/*endif:add to RESPA list*/
    }/*endfor: jpart interaction*/
    (*nver_now)+=intact_use;

  }/*endif:RESPA on*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Shave the interactions down                                              */
/*==========================================================================*/

void vercreate_list_count_root(
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter, int *nbrnch_of_root)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   int jpart,ic,intact_use,ifdx,ivdx,ntemp,iii;
   double cutoff,r;

/*=======================================================================*/
/* II) Shrink and count the list no respa                               */

  if(int_res_ter==0){

    ic = 0;
    for (jpart = 1; jpart <= intact_now; ++jpart) {
      cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
      r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
              + dz12[jpart]*dz12[jpart];
      if(r < cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
      }/*endif*/
    }/*endfor*/
    intact_use = ic;

    ic = 0;
    for (jpart = 1; jpart <= intact_use; ++jpart) {
      ifdx = for_scr_i_index[jpart];
      ntemp = nbrnch_of_root[for_scr_j_index[jpart]];
      nbr_list_nter[ifdx] += ntemp;
      ic += ntemp;
    }/*endfor:jpart interaction*/
    (*nver_now)+=ic;
  }/*endif:RESPA off*/

/*=======================================================================*/
/* II) Shrink and count the list respa                                   */

  if(int_res_ter==1){

    ic = 0;
    for (jpart = 1; jpart <= intact_now; ++jpart) {
      cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
      r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
              + dz12[jpart]*dz12[jpart];
      if(r<cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
        intra_scr_cutoff_res[ic] = 
                            interact_cutskin_res[(for_scr_i_indext[jpart])];
        intra_scr_r[ic]     = r;
      }/*endif*/
    }/*endfor*/
    intact_use = ic;

    ic = 0;
    for (jpart = 1; jpart <= intact_use; ++jpart) {
      ifdx                = for_scr_i_index[jpart];
      ntemp               = nbrnch_of_root[for_scr_j_index[jpart]];
      nbr_list_nter[ifdx]+= ntemp;
      ic                 += ntemp;
      if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
        nbr_list_nter_res[ifdx] += ntemp;
        (*nver_res_now)         += ntemp;
      }/*endif:add to RESPA list*/
    }/*endfor: jpart interaction*/
    (*nver_now)+=ic;

  }/*endif:RESPA on*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Shave the interactions down                                              */
/*==========================================================================*/

void vercreate_list_fill(
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter,int num_brk,int *num_brk_i)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   int jpart,ic,intact_use,ifdx,ivdx;
   double cutoff,r,cutoff_res;
   int upper,lower,ibrk,iii;

/*=======================================================================*/
/* II) Shrink and fill the list no respa                                 */

  if(int_res_ter==0){

     ic = 0;
     for (jpart = 1; jpart <= intact_now; ++jpart) {
       cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
       r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
               + dz12[jpart]*dz12[jpart];
       if(r < cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
       }/*endif*/
     }/*endfor*/
     intact_use = ic;

     for (jpart = 1; jpart <= intact_use; ++jpart) {
       ifdx                   = for_scr_i_index[jpart];
       nbr_list_nter[ifdx]++;
       for_scr_i_index[jpart] = nbr_list_nter[ifdx];
     }/*endfor: jpart interactions*/

#ifndef NO_PRAGMA
#pragma ivdep
#endif

     for (jpart = 1; jpart <= intact_use; ++jpart) {
       nbr_list_jter[for_scr_i_index[jpart]] = for_scr_j_index[jpart];
     }/*endfor: jpart interactions*/
     (*nver_now)+=intact_use;

  }/*endif*/

/*=======================================================================*/
/* II) Shrink and fill the list respa                                    */

  if(int_res_ter==1){

     ic = 0;
     for (jpart = 1; jpart <= intact_now; ++jpart) {
       cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
       r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
               + dz12[jpart]*dz12[jpart];
       if(r<cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
        intra_scr_cutoff_res[ic] = 
                            interact_cutskin_res[(for_scr_i_indext[jpart])];
        intra_scr_r[ic]     = r;
       }/*endif*/
     }/*endfor*/
     intact_use = ic;

     for (jpart = 1; jpart <= intact_use; ++jpart) {
       ifdx = for_scr_i_index[jpart];
       nbr_list_nter[ifdx]++;
       ivdx =  nbr_list_nter[ifdx];
       nbr_list_jter[ivdx] = for_scr_j_index[jpart];
       if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
        nbr_list_nter_res[ifdx]++;
        (*nver_res_now)++;
        nbr_list_jter_res[(nbr_list_nter_res[ifdx])] = for_scr_j_index[jpart];
       }/*endif: add to RESPA list*/
     }/*endfor: jpart interaction*/
     (*nver_now)+=intact_use;

  }/*endif*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Shave the interactions down                                              */
/*==========================================================================*/

void vercreate_list_fill_root(
              double *dx12,double *dy12,double *dz12,
              int intact_now,int *for_scr_i_index,int *for_scr_j_index,
              int *for_scr_i_indext,
              double *interact_cutskin,double *interact_cutskin_res,
              double *intra_scr_r,double *intra_scr_cutoff_res,
              double *intra_scr_cutoff,
              int *nbr_list_nter,list_int *nbr_list_jter,int *nver_now,
              int *nbr_list_nter_res,list_int *nbr_list_jter_res,
              int *nver_res_now,
              int int_res_ter,int num_brk,int *num_brk_i,
              int *nbrnch_of_root)

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

   int jpart,ic,intact_use,ifdx,ivdx;
   double cutoff,r,cutoff_res;
   int upper,lower,ibrk;
   int ntemp;

/*=======================================================================*/
/* II) Shrink and fill the list no respa                                 */

  if(int_res_ter==0){

     ic = 0;
     for (jpart = 1; jpart <= intact_now; ++jpart) {
       cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
       r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
               + dz12[jpart]*dz12[jpart];
       if(r < cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
       }/*endif*/
     }/*endfor*/
     intact_use = ic;

     ic = 0;
     for (jpart = 1; jpart <= intact_use; ++jpart) {
       ifdx                   = for_scr_i_index[jpart];
       ntemp                  = nbrnch_of_root[for_scr_j_index[jpart]];
       for_scr_i_index[jpart] = (nbr_list_nter[ifdx]+1);
       nbr_list_nter[ifdx]   += ntemp;
       ic                    += ntemp;
     }/*endfor: jpart interactions*/
     (*nver_now)+=ic;

#ifndef NO_PRAGMA
#pragma ivdep
#endif
     for (jpart = 1; jpart <= intact_use; ++jpart) {
       nbr_list_jter[for_scr_i_index[jpart]] = for_scr_j_index[jpart];
     }/*endfor: jpart interactions*/

  }/*endif:RESPA off*/

/*=======================================================================*/
/* II) Shrink and fill the list respa                                    */

  if(int_res_ter==1){

     ic = 0;
     for (jpart = 1; jpart <= intact_now; ++jpart) {
       cutoff  =  interact_cutskin[(for_scr_i_indext[jpart])];
       r       = dx12[jpart]*dx12[jpart]+dy12[jpart]*dy12[jpart]
               + dz12[jpart]*dz12[jpart];
       if(r<cutoff){
        ic++;
        for_scr_i_index[ic] = for_scr_i_index[jpart];
        for_scr_j_index[ic] = for_scr_j_index[jpart];
        intra_scr_cutoff_res[ic] = 
                            interact_cutskin_res[(for_scr_i_indext[jpart])];
        intra_scr_r[ic]     = r;
       }/*endif*/
     }/*endfor*/
     intact_use = ic;

     ic = 0;
     for (jpart = 1; jpart <= intact_use; ++jpart) {
       ifdx                 = for_scr_i_index[jpart];
       ntemp                = nbrnch_of_root[for_scr_j_index[jpart]];
       ivdx                 = (nbr_list_nter[ifdx]+1);
       nbr_list_jter[ivdx]  = for_scr_j_index[jpart];
       nbr_list_nter[ifdx] += ntemp;
       ic                  += ntemp;
       if(intra_scr_r[jpart] < intra_scr_cutoff_res[jpart]){
        ivdx                    = (nbr_list_nter_res[ifdx]+1);
        nbr_list_jter_res[ivdx] =  for_scr_j_index[jpart];
        nbr_list_nter_res[ifdx]+=ntemp;
        (*nver_res_now)        +=ntemp;
       }/*endif: add to RESPA list*/
     }/*endfor: jpart interaction*/
     (*nver_now)+=ic;

  }/*endif*/


/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/

