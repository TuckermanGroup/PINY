/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_pot_nonlocal.c                       */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

#define HESS_OFF
#define DEBUG_GJM_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_ewd_non_loc(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                      CELL *cell, PTENS *ptens, CPEWALD *cpewald,
                      CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,  
                      CPOPTS *cpopts, ATOMMAPS *atommaps, 
                      COMMUNICATE *communicate,FOR_SCR *for_scr)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/

/*=======================================================================*/
/*         Local Variable declarations                                   */
     
#include "../typ_defs/typ_mask.h"
  double arg;
  double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
  double xtemp,ytemp,ztemp;
  double g2,tpi,pi,g,fpi;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;

  double temp;
  double xtrans,ytrans,ztrans;
  int nl_chan_max,irad,jrad;
  int ipart,iii,i,iatm;
  int nl_max,np_nlmax;
  int icount;
  YLM_CONS ylm_cons;
  double ylmr[2];

/*         Local Pointer declarations                                   */
  double *creal_up   = cpcoeffs_pos->cre_up;
  double *cimag_up   = cpcoeffs_pos->cim_up;
  double *creal_dn   = cpcoeffs_pos->cre_dn;
  double *cimag_dn   = cpcoeffs_pos->cim_dn;
  int icoef_orth_up  = cpcoeffs_pos->icoef_orth_up;
  int icoef_form_up  = cpcoeffs_pos->icoef_form_up;
  int ifcoef_form_up = cpcoeffs_pos->ifcoef_form_up;
  int icoef_orth_dn  = cpcoeffs_pos->icoef_orth_dn;
  int icoef_form_dn  = cpcoeffs_pos->icoef_form_dn;
  int ifcoef_form_dn = cpcoeffs_pos->ifcoef_form_dn;
  int nstate_up      = cpcoeffs_info->nstate_up_proc;
  int nstate_dn      = cpcoeffs_info->nstate_dn_proc;
  int ncoef          = cpcoeffs_info->ncoef;
  int cp_hess_calc   = cpopts->cp_hess_calc;
  int cp_lsda        = cpopts->cp_lsda;

  int npart          = clatoms_info->natm_tot;
  double *x          = clatoms_pos->x;
  double *y          = clatoms_pos->y;
  double *z          = clatoms_pos->z;
  double *hmat_cp    = cell->hmat_cp;
  double *hmati_cp   = cell->hmati_cp;
  double *hmat_big  = cell->hmat;
  double *hmati_big  = cell->hmati;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;
  int natm_typ       = atommaps->natm_typ;
  int *iatm_typ      = atommaps->iatm_atm_typ;
  int *iatm_typ_nl   = atommaps->iatm_atm_typ_nl;
  int *index_atm     = for_scr->iexcl;

  int nktot_sm       = cpewald->nktot_sm;
  int *kastore_sm    = cpewald->kastr_sm;
  int *kbstore_sm    = cpewald->kbstr_sm;
  int *kcstore_sm    = cpewald->kcstr_sm;
  double *ak2_sm     = cpewald->ak2_sm;
  int *ibreak1_sm    = cpewald->ibrk1_sm;
  int *ibreak2_sm    = cpewald->ibrk2_sm;
  double *cossc      = ewd_scr->cossc;  
  double *sinsc      = ewd_scr->sinsc;
  double *helr       = ewd_scr->helr;   
  double *heli       = ewd_scr->heli;
  double *vtemp      = ewd_scr->fx2;
  double *vtemp_now  = ewd_scr->vtemp_now;
  double *dvtemp     = ewd_scr->q;
  double *ewd_scr_x  = ewd_scr->x;
  double *ewd_scr_y  = ewd_scr->y;
  double *ewd_scr_z  = ewd_scr->z;
  int n_ang_max      = pseudo->n_ang_max;
  int n_ang_max_kb   = pseudo->n_ang_max_kb;
  int n_rad_max      = pseudo->n_rad_max;
  int *ip_nl         = pseudo->ip_nl;
  int *np_nl         = pseudo->np_nl;
  double *gzvps      = pseudo->gzvps;
  double *gzvps0     = pseudo->gzvps0;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  int *nrad_max_l       = pseudo->nrad_max_l;
  int **np_nl_rad_str   = pseudo->np_nl_rad_str;
  double *vnlreal_up = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn = cpscr->cpscr_nonloc.vnlim_dn;

  int myid_state     = communicate->myid_state;
  int np_states      = communicate->np_states;
  

/*======================================================================*/
/* 0) Check the forms                                                   */

  if(icoef_orth_up!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The UP coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_ewd_non_loc   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
  }/*endif*/
  if(cp_lsda==1){
   if(icoef_orth_dn!=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Dn coefficients must be in orthogonal form    \n");
    printf("on state processor %d in control_ewd_non_loc   \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
  }/*endif*/

  if(np_states>1){
   if((icoef_form_up+ifcoef_form_up)!=0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("The Up coefs and coef forces must not be in transposed form \n");
    printf("on state processor %d in control_ewd_non_loc  \n",myid_state);
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);exit(1);
   }/*endif*/
   if(cp_lsda==1){
    if((icoef_form_up+ifcoef_form_up)!=0){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The Up coefs and coef forces must not be in transposed form \n");
     printf("on state processor %d in control_ewd_non_loc  \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
    }/*endif*/
   }/*endif*/
  }/*endif*/

/*======================================================================*/
/* I) Get some useful constants                                         */

  tpi                 = 2.0*M_PI;
  fpi                 = 4.0*M_PI;
  ylm_cons.rt_fpi     = 1.0/sqrt(fpi);
  ylm_cons.rt_thrfpi  = sqrt(3.0/fpi);
  ylm_cons.rt_threpi  = sqrt(1.50/fpi);
  ylm_cons.hrt_fivfpi = 0.50*sqrt(5.0/fpi);
  ylm_cons.rt_fiftepi = sqrt(7.50/fpi);
  ylm_cons.hrt_sevfpi = 0.50*sqrt(7.0/fpi);
  ylm_cons.hrt_toepi  = 0.50*sqrt(10.50/fpi)/sqrt(2.0);
  ylm_cons.hrt_ohffpi = 0.50*sqrt(105.0/fpi)/sqrt(2.0);
  ylm_cons.hrt_tfepi  = 0.50*sqrt(17.50/fpi)/sqrt(2.0);

  

/*======================================================================*/
/* II) Determine the maximum open non-local angular momentum channel     */


  nl_max = -1;


  for(i=1;i<=(n_ang_max_kb +1);i++){
    if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/

  nl_chan_max = (nl_max + 1)*(nl_max + 1);

/*======================================================================*/
/* III) Determine the maximum number of atoms in any                    */
/*       open angular momentum channel                                  */

  np_nlmax = 1;
  for(i = 1;i<=(nl_max+1);i++){
    np_nlmax = MAX(np_nlmax,np_nl[i]);
  }/*endfor*/


/*======================================================================*/
/* IV) Find cos and sin of sc components of the particles               */
/*    ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  for(ipart=1;ipart<= np_nonloc_cp_box_kb;ipart++){
    iatm = ip_nl[ipart];
    dx  = x[iatm] - cp_box_center[1];
    dy  = y[iatm] - cp_box_center[2];
    dz  = z[iatm] - cp_box_center[3];
    asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
    asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
    asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
    sx  = asx - NINT(asx);
    sy  = asy - NINT(asy);
    sz  = asz - NINT(asz);
    dx  = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
    dy  = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
    dz  = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];
    xtemp = dx + cp_box_center_rel[1];
    ytemp = dy + cp_box_center_rel[2];
    ztemp = dz + cp_box_center_rel[3];
    ewd_scr_x[ipart] = xtemp*hmati_cp[1]
                     + ytemp*hmati_cp[4]
                     + ztemp*hmati_cp[7];
    ewd_scr_y[ipart] = xtemp*hmati_cp[2]
                     + ytemp*hmati_cp[5]
                     + ztemp*hmati_cp[8];
    ewd_scr_z[ipart] = xtemp*hmati_cp[3]
                     + ytemp*hmati_cp[6]
                     + ztemp*hmati_cp[9];
    ctemp = ewd_scr_z[ipart]*tpi;
    cossc[ipart] = cos(ctemp);
    sinsc[ipart] = sin(ctemp);
  }/*endfor*/

/*======================================================================*/
/* V) Perform the ewald sum/ non-local potential calculation            */

 for(icount=1;icount<=nktot_sm;icount++){
  
/*----------------------------------------------------------------------*/
/* i) Get the k vectors                                                 */

   aka = (double)(kastore_sm[icount]);
   akb = (double)(kbstore_sm[icount]);
   akc = (double)(kcstore_sm[icount]);
   xk = (aka*hmati_cp[1]+akb*hmati_cp[2]+akc*hmati_cp[3])*tpi;
   yk = (aka*hmati_cp[4]+akb*hmati_cp[5]+akc*hmati_cp[6])*tpi;
   zk = (aka*hmati_cp[7]+akb*hmati_cp[8]+akc*hmati_cp[9])*tpi;
   g2 = xk*xk+yk*yk+zk*zk;
   g  = sqrt(g2);
   ak2_sm[icount] = g2;

/*----------------------------------------------------------------------*/
/* ii) If break point number one calculate the helpful vectors          */
 
   if(ibreak1_sm[icount]==1){
    for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
     atemp = ewd_scr_x[ipart];
     btemp = ewd_scr_y[ipart];
     ctemp = ewd_scr_z[ipart];
     arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
     helr[ipart] = cos(arg);
     heli[ipart] = sin(arg);
    }/*endfor*/
   }/*endif*/

/*----------------------------------------------------------------------*/
/* iii) nonlocal matrix                                                 */
/*     use structure factor and wavefunction to compute                 */
/*     the array znl and its derivative dznl for the                    */
/*     calculation of the nonlocal forces and energy                    */
/*     (done separately for each angular momentum component)            */


   control_nlmat(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
                 cpscr,cpopts,pseudo,ewd_scr,atommaps,
                 np_nlmax,nl_max,index_atm,g,xk,yk,zk,
                 icount,&ylm_cons);

/*----------------------------------------------------------------------*/
/* iv) If break point two, increment the helpful vectors                */

   if(ibreak2_sm[icount]==1){
     for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
     temp = helr[ipart];
     helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
     heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
    }/*endfor*/
   }/*endif*/

 }/*endfor:icount loop over k vectors */

/*======================================================================*/
/* VI) g=0 term                                                         */

  ak2_sm[ncoef] = 0.0;
  if(np_nl[1]>0){

    ylmr[1] = 1.0/sqrt(fpi);
    for(irad=1;irad<=nrad_max_l[1];irad++){

      for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
        iii = n_rad_max*(iatm_typ_nl[ipart]-1) + irad;
        vtemp[ipart] = gzvps0[iii];
      }/*endfor*/

      for(ipart=np_nl_rad_str[1][irad];ipart<=np_nl[1];ipart++){
        vtemp_now[ipart] = vtemp[ipart]*ylmr[1];
      }/*endfor*/

      get_nlmat0(ncoef,ncoef,nstate_up,np_nlmax,
                 np_nl_rad_str[1][irad],np_nl[1],nl_chan_max,irad,
                 creal_up,cimag_up,vtemp_now,vnlreal_up,vnlimag_up,ylmr[1]); 
      if(cp_lsda==1){
        get_nlmat0(ncoef,ncoef,nstate_dn,np_nlmax,
                  np_nl_rad_str[1][irad],np_nl[1],nl_chan_max,irad,
                  creal_dn,cimag_dn,vtemp_now,vnlreal_dn,vnlimag_dn,ylmr[1]);
      }/*endif*/

    }/*endfor*/

  }/*endif: l=0 nonlocal*/

/*======================================================================*/
   }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_nlmat(CLATOMS_INFO *clatoms_info,
                   CPCOEFFS_INFO *cpcoeffs_info,
                   CPCOEFFS_POS *cpcoeffs_pos,
                   CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                   EWD_SCR *ewd_scr,ATOMMAPS *atommaps,
                   int np_nlmax,int nl_max,int *index_atm,
                   double g,double xk,double yk,double zk,
                   int ismcount,YLM_CONS *ylm_cons)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/

/*=======================================================================*/
/*         Local Variable declarations                                   */

  int ind_lm,ipart,lp1,irad,jrad,ipart_nl;
  double sgn_l;
  double ylmr[21], ylmi[21];
  double dylmr_gx[21], dylmi_gx[21];
  double dylmr_gy[21], dylmi_gy[21];
  double dylmr_gz[21], dylmi_gz[21];
  int itype,itype_nl,l,m,i_shift,iii;
  int ktemp,ltemp;
  int nl_chan_max;
/* Local pointers */     
  int npart                = clatoms_info->natm_tot;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  int *iatm_typ_nl_rev     = atommaps->iatm_atm_typ_nl_rev;
  int *imap_atm_typ_nl     = atommaps->imap_atm_typ_nl;
  int  natm_typ_nl         = pseudo->natm_typ_nl;

  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;
  int ncoef                = cpcoeffs_info->ncoef;
  double *creal_up         = cpcoeffs_pos->cre_up;  
  double *cimag_up         = cpcoeffs_pos->cim_up;
  double *creal_dn         = cpcoeffs_pos->cre_dn;
  double *cimag_dn         = cpcoeffs_pos->cim_dn;
  int cp_lsda              = cpopts->cp_lsda;
  int cp_hess_calc         = cpopts->cp_hess_calc;
  int cp_ptens             = cpopts->cp_ptens_calc;
  int atm_hess_calc        = clatoms_info->hess_calc;

  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int nsplin_g             = pseudo->nsplin_g;
  double dg_spl            = pseudo->dg_spl;
  double gmin_spl          = pseudo->gmin_spl;
  double *vps0             = pseudo->vps0;
  double *vps1             = pseudo->vps1;
  double *vps2             = pseudo->vps2;
  double *vps3             = pseudo->vps3;
  double *dvps0            = pseudo->dvps0;
  double *dvps1            = pseudo->dvps1;
  double *dvps2            = pseudo->dvps2;
  double *dvps3            = pseudo->dvps3;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;
 
  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *helr             = ewd_scr->helr;   
  double *heli             = ewd_scr->heli;
  double *helr_now         = ewd_scr->helr_now;
  double *heli_now         = ewd_scr->heli_now;
  double *dheli_now        = ewd_scr->temp;     
  double *dhelr_now        = ewd_scr->vtemp_now; 
  double *vtemp            = ewd_scr->fx2;
  double *dvtemp           = ewd_scr->q;
  double *scr1             = ewd_scr->fx2; /* yes the same as vtemp, its ok */
  double *scr2             = ewd_scr->q;   /* yes the same as dvtemp,its ok */
  double *scr3             = ewd_scr->fz2;  

/*======================================================================*/
/* I) Get the ylm(g)                                                   */

  get_ylm(xk,yk,zk,g,ylmr,ylmi,dylmr_gx,dylmi_gx,dylmr_gy,dylmi_gy,
          dylmr_gz,dylmi_gz,ylm_cons);

  nl_chan_max = (nl_max + 1)*(nl_max + 1);
/*======================================================================*/
/* II) Calculate the nl-pseudoponential matrix elements by looping over  */
/* the channels, l, and then the 2l+1 components of the channel         */

  for(l=0;l<=nl_max;l++){

    lp1 = l+1;
    if(np_nl[lp1]>0){

     for(irad=1;irad<=nrad_max_l[lp1];irad++){
/*---------------------------------------------------------------------*/
/* i) Get the bessel transform of the pseudopotential at this g vector */
/*    and scale the structure factor appropriately                     */

      for(itype=1;itype<=natm_typ_nl;itype++){
        itype_nl = imap_atm_typ_nl[itype];
        index_atm[itype] =  (itype_nl-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                           + l*nsplin_g*n_rad_max +  (irad-1)*nsplin_g;
      }/*endfor*/

      get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                 vps0,vps1,vps2,vps3,vtemp,iatm_typ_nl_rev,
                 natm_typ_nl,np_nonloc_cp_box_kb,
                 np_nl_rad_str[lp1][irad]);

      i_shift  = l*npart;

      if(cp_ptens==0) {

        for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
          ktemp             =  ipart+i_shift;
          ltemp             =  ip_nl_rev[ktemp];
          helr_now[ipart]   =  helr[ltemp]*vtemp[ltemp];
          heli_now[ipart]   =  heli[ltemp]*vtemp[ltemp];
        }/*endfor*/
      }else{
	/*PRESSURE TENSOR CALC*/
        get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                   dvps0,dvps1,dvps2,dvps3,dvtemp,iatm_typ,
                   natm_typ,npart,np_nl_rad_str[lp1][irad]);
        for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
          ktemp              = ipart+i_shift;
          ltemp             =  ip_nl_rev[ktemp];
          helr_now[ipart]    = helr[ltemp]*vtemp[ltemp];
          heli_now[ipart]    = heli[ltemp]*vtemp[ltemp];
          dhelr_now[ipart]   = helr[ltemp]*dvtemp[ltemp];
          dheli_now[ipart]   = heli[ltemp]*dvtemp[ltemp];
        }/*endfor*/
      }/*endif*/
    
/*---------------------------------------------------------------------*/
/* ii) Loop over m components of the channel and get the matrix elements */ 
/*     vtemp,dvtemp, */
      sgn_l = 1.0;
      if((l%2)==1)sgn_l = -1.0;
      for(m=1;m<=(2*l+1);m++){
        ind_lm = m + l*l;
        if(cp_ptens==1) {
          get_nlmat_pv(ncoef,ismcount,nstate_up,ind_lm,irad,
            np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
            helr_now,heli_now,
            creal_up,cimag_up,dhelr_now,dheli_now,
            vnlreal_up,vnlimag_up,dvnlreal_x_up,dvnlreal_y_up,
            dvnlreal_z_up,dvnlimag_x_up,dvnlimag_y_up,
            dvnlimag_z_up,dvnlreal_gxgx_up,dvnlreal_gxgy_up,
            dvnlreal_gxgz_up,dvnlreal_gygy_up,dvnlreal_gygz_up,
            dvnlreal_gzgz_up,
            dvnlimag_gxgx_up,dvnlimag_gxgy_up,dvnlimag_gxgz_up,
            dvnlimag_gygy_up,dvnlimag_gygz_up,dvnlimag_gzgz_up,
            xk,yk,zk,
            ylmr[ind_lm],ylmi[ind_lm],
            dylmr_gx[ind_lm],dylmi_gx[ind_lm],
            dylmr_gy[ind_lm],dylmi_gy[ind_lm],
            dylmr_gz[ind_lm],dylmi_gz[ind_lm],
            sgn_l,scr1,scr2,scr3);
        }/* endif */
        if(atm_hess_calc == 3){
          get_nlmat_hess(ncoef,ismcount,nstate_up,ind_lm,irad,
            np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
            helr_now,heli_now,
            creal_up,cimag_up,dhelr_now,dheli_now,
            vnlreal_up,vnlimag_up,dvnlreal_x_up,dvnlreal_y_up,
            dvnlreal_z_up,dvnlimag_x_up,dvnlimag_y_up,
            dvnlimag_z_up,dvnlreal_gxgx_up,dvnlreal_gxgy_up,
            dvnlreal_gxgz_up,dvnlreal_gygy_up,dvnlreal_gygz_up,
            dvnlreal_gzgz_up,
            dvnlimag_gxgx_up,dvnlimag_gxgy_up,dvnlimag_gxgz_up,
            dvnlimag_gygy_up,dvnlimag_gygz_up,dvnlimag_gzgz_up,
            xk,yk,zk,
            ylmr[ind_lm],ylmi[ind_lm],
            sgn_l,scr1,scr2);
        }/* endif */
        if(atm_hess_calc != 3 && cp_ptens == 0){
          get_nlmat(ncoef,ismcount,nstate_up,ind_lm,irad,
            np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
            helr_now,heli_now,
            creal_up,cimag_up,
            vnlreal_up,vnlimag_up,
            dvnlreal_x_up,dvnlreal_y_up,dvnlreal_z_up,
            dvnlimag_x_up,dvnlimag_y_up,dvnlimag_z_up,
            xk,yk,zk,
            ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
        }/*endif*/
      
        if(cp_lsda==1) {
          if(cp_ptens==1) {
            get_nlmat_pv(ncoef,ismcount,nstate_dn,ind_lm,irad,
              np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
              helr_now,heli_now,
              creal_dn,cimag_dn,dhelr_now,dheli_now,
              vnlreal_dn,vnlimag_dn,
              dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
              dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
              dvnlreal_gxgx_dn,dvnlreal_gxgy_dn,dvnlreal_gxgz_dn,
              dvnlreal_gygy_dn,dvnlreal_gygz_dn,dvnlreal_gzgz_dn,
              dvnlimag_gxgx_dn,dvnlimag_gxgy_dn,dvnlimag_gxgz_dn,
              dvnlimag_gygy_dn,dvnlimag_gygz_dn,dvnlimag_gzgz_dn,
              xk,yk,zk,
              ylmr[ind_lm],ylmi[ind_lm],
              dylmr_gx[ind_lm],dylmi_gx[ind_lm],
              dylmr_gy[ind_lm],dylmi_gy[ind_lm],
              dylmr_gz[ind_lm],dylmi_gz[ind_lm],sgn_l,scr1,scr2,scr3);
          }/* endif */
          if(atm_hess_calc == 3){
            get_nlmat_hess(ncoef,ismcount,nstate_dn,ind_lm,irad,
              np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
              helr_now,heli_now,
              creal_dn,cimag_dn,dhelr_now,dheli_now,
              vnlreal_dn,vnlimag_dn,
              dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
              dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
              dvnlreal_gxgx_dn,dvnlreal_gxgy_dn,dvnlreal_gxgz_dn,
              dvnlreal_gygy_dn,dvnlreal_gygz_dn,dvnlreal_gzgz_dn,
              dvnlimag_gxgx_dn,dvnlimag_gxgy_dn,dvnlimag_gxgz_dn,
              dvnlimag_gygy_dn,dvnlimag_gygz_dn,dvnlimag_gzgz_dn,
              xk,yk,zk,
              ylmr[ind_lm],ylmi[ind_lm],
              sgn_l,scr1,scr2);
          }/* endif */
          if(atm_hess_calc != 3 && cp_ptens == 0){
             get_nlmat(ncoef,ismcount,nstate_dn,ind_lm,irad,
               np_nlmax,np_nl[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
               helr_now,heli_now,
               creal_dn,cimag_dn,
               vnlreal_dn,vnlimag_dn,
               dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
               dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,xk,yk,zk,
               ylmr[ind_lm],ylmi[ind_lm],sgn_l,scr1,scr2);
          }/*endif:pvten*/
        }/*endif:lsda*/
      }/*endfor:m quantum number loop*/
     }/*endfor: irad quantum number loop */
    }/*endif: this channel has atoms in it*/
  }/*endfor:l quantum number loop*/   

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlmat(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
               int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
               double *helr,double *heli,
               double *creal,double *cimag,
               double *vnlreal,double *vnlimag,
               double *dvnlreal_x,double *dvnlreal_y,double *dvnlreal_z,
               double *dvnlimag_x,double *dvnlimag_y,double *dvnlimag_z,
               double xk,double yk,double zk,double ylmr,double ylmi,
               double sgn_l,double *cbyhelr,double *cbyheli)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

 
  int is,i,ipart,iii;
  int ind_loc_c,ioff_v,ind_loc_v;
  int isgn_l;
  double cre_ind,cim_ind;
  double tylmr, tylmi;
  double tmp1,tmp2,tmp3,tmp4;
/*==========================================================================*/
/* I) Loop over states and particles:                                       */
/*     Here helr,heli have vtemp mutliplied in them                         */
/*     See control routine                                                  */
/* storage : atm,state,l,m,n */

  isgn_l = (int) (0.5*(sgn_l+1.0)+1.0);
  tylmr = 2.0*ylmr;
  tylmi = 2.0*ylmi;

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind   = creal[ind_loc_c];
    cim_ind   = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    switch(isgn_l){
       case 1: 
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           cbyhelr[ipart]     =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
           cbyheli[ipart]     =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
         }/* endfor : loop over particles */ 
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        = -cbyheli[i]*tylmi;
           tmp2        =  cbyheli[i]*tylmr;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
         }/*endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp3        = -cbyhelr[i]*tylmi;
            dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
            dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
            dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp4        =  cbyhelr[i]*tylmr;
            dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
            dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
            dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
         }/* endfor */
       break;
       case 2: 
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           cbyhelr[ipart]     =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
           cbyheli[ipart]     =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
         }/* endfor : loop over particles */ 
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        =  cbyhelr[i]*tylmr;
           tmp2        =  cbyhelr[i]*tylmi;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
         }/*endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp3        = -cbyheli[i]*tylmr;
           dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
           dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
           dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp4        = -cbyheli[i]*tylmi;
           dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
           dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
           dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
         }/* endfor */
    }/*end switch: sgn of m in Ylm*/

  }/* endfor : loop over states */
   
/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlmat_pv(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
         int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
         double *helr,double *heli,
         double *creal,double *cimag,
         double *dhelr,double *dheli,
         double *vnlreal,double *vnlimag,
         double *dvnlreal_x,double *dvnlreal_y,
         double *dvnlreal_z,double *dvnlimag_x,double *dvnlimag_y,
         double *dvnlimag_z,double *dvnlreal_gxgx,double *dvnlreal_gxgy,
         double *dvnlreal_gxgz,double *dvnlreal_gygy,double *dvnlreal_gygz,
         double *dvnlreal_gzgz,double *dvnlimag_gxgx,double *dvnlimag_gxgy,
         double *dvnlimag_gxgz,double *dvnlimag_gygy,double *dvnlimag_gygz,
         double *dvnlimag_gzgz,double xk,double yk,double zk,
         double ylmr,double ylmi,double dylmr_gx,double dylmi_gx,
         double dylmr_gy,double dylmi_gy,double dylmr_gz,double dylmi_gz,
         double sgn_l,double *cbyhelr,double *cbyheli,double *dcbyhel)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  int is,i,ipart,iii;
  int ind_loc_c,ioff_v,ind_loc_v;
  int isgn_l ;
  double cre_ind,cim_ind;
  double tmp1,tmp2,tmp3,tmp4;
  double dtmp1,dtmp2;
  double tylmr,tylmi;
  double tdylmr_gx,tdylmr_gy,tdylmr_gz;
  double tdylmi_gx,tdylmi_gy,tdylmi_gz,dtmp3;

  double *dcbyhelr = dcbyhel;  /* These are not used together */
  double *dcbyheli = dcbyhel;

/*==========================================================================*/
/* I) Loop over states and particles                                        */
/*     Here helr,heli have vtemp mutliplied in them                         */
/*     Here dhelr,dheli have dvtemp mutliplied in them                      */
/*     See control routine                                                  */

  isgn_l    = (int) (0.5*(sgn_l+1.0)+1.0);
  tylmr     = 2.0*ylmr;
  tylmi     = 2.0*ylmi;
  tdylmr_gx = 2.0*dylmr_gx;
  tdylmr_gy = 2.0*dylmr_gy;
  tdylmr_gz = 2.0*dylmr_gz;
  tdylmi_gx = 2.0*dylmi_gx;
  tdylmi_gy = 2.0*dylmi_gy;
  tdylmi_gz = 2.0*dylmi_gz;

/* storage : atm,state,l,m,n */

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind = creal[ind_loc_c];
    cim_ind = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;

    switch(isgn_l){
      case 1: 
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]    =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]    =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
          dcbyheli[ipart]   =  cre_ind*dheli[ipart] + cim_ind*dhelr[ipart];
        }/* endfor : loop over particles */ 
        for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        = -cbyheli[i]*tylmi;
           tmp2        =  cbyheli[i]*tylmr;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
        }/*endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp3       = -cbyhelr[i]*tylmi;
            dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
            dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
            dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
        }/* endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp4       =  cbyhelr[i]*tylmr;
            dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
            dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
            dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
        }/* endfor */

        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = (-dcbyheli[i]*tylmi*xk - cbyheli[i]*tdylmi_gx);
           dvnlreal_gxgx[(i+ioff_v)] += (dtmp3*xk);
           dvnlreal_gxgy[(i+ioff_v)] += (dtmp3*yk);
           dvnlreal_gxgz[(i+ioff_v)] += (dtmp3*zk);
        }/*endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp1     = -dcbyheli[i]*tylmi;
           dtmp3     = (dtmp1*yk - cbyheli[i]*tdylmi_gy);
           dvnlreal_gygy[(i+ioff_v)] += (dtmp3*yk);
           dvnlreal_gygz[(i+ioff_v)] += (dtmp3*zk);
           dvnlreal_gzgz[(i+ioff_v)] += ((dtmp1*zk - cbyheli[i]*tdylmi_gz)*zk);
        }/* endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = (dcbyheli[i]*tylmr*xk + cbyheli[i]*tdylmr_gx);
           dvnlimag_gxgx[(i+ioff_v)] += (dtmp3*xk);
           dvnlimag_gxgy[(i+ioff_v)] += (dtmp3*yk);
           dvnlimag_gxgz[(i+ioff_v)] += (dtmp3*zk);
        }/* endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp2     =  dcbyheli[i]*tylmr;
           dtmp3     = (dtmp2*yk + cbyheli[i]*tdylmr_gy);
           dvnlimag_gygy[(i+ioff_v)] += (dtmp3*yk);
           dvnlimag_gygz[(i+ioff_v)] += (dtmp3*zk);
           dvnlimag_gzgz[(i+ioff_v)] += ((dtmp2*zk + cbyheli[i]*tdylmr_gz)*zk);
        }/* endfor */
      break;
      case 2: 
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]    =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]    =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
          dcbyhelr[ipart]   =  cre_ind*dhelr[ipart] - cim_ind*dheli[ipart];
        }/* endfor : loop over particles */ 
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        =  cbyhelr[i]*tylmr;
           tmp2        =  cbyhelr[i]*tylmi;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
         }/*endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp3        = -cbyheli[i]*tylmr;
           dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
           dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
           dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp4       = -cbyheli[i]*tylmi;
            dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
            dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
            dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
         }/* endfor */
         
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = (dcbyhelr[i]*tylmr*xk + cbyhelr[i]*tdylmr_gx);
           dvnlreal_gxgx[(i+ioff_v)] += (dtmp3*xk);
           dvnlreal_gxgy[(i+ioff_v)] += (dtmp3*yk);
           dvnlreal_gxgz[(i+ioff_v)] += (dtmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp1     =  dcbyhelr[i]*tylmr;
           dtmp3     = (dtmp1*yk + cbyhelr[i]*tdylmr_gy);
           dvnlreal_gygy[(i+ioff_v)] += (dtmp3*yk);
           dvnlreal_gygz[(i+ioff_v)] += (dtmp3*zk);
           dvnlreal_gzgz[(i+ioff_v)] += ((dtmp1*zk + cbyhelr[i]*tdylmr_gz)*zk);
         }/*endfor */
 
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = (dcbyhelr[i]*tylmi*xk + cbyhelr[i]*tdylmi_gx);
           dvnlimag_gxgx[(i+ioff_v)] += (dtmp3*xk);
           dvnlimag_gxgy[(i+ioff_v)] += (dtmp3*yk);
           dvnlimag_gxgz[(i+ioff_v)] += (dtmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp2     =  dcbyhelr[i]*tylmi;
           dtmp3     = (dtmp2*yk + cbyhelr[i]*tdylmi_gy);
           dvnlimag_gygy[(i+ioff_v)] += (dtmp3*yk);
           dvnlimag_gygz[(i+ioff_v)] += (dtmp3*zk);
           dvnlimag_gzgz[(i+ioff_v)] += ((dtmp2*zk + cbyhelr[i]*tdylmi_gz)*zk);
         }/* endfor */
    }/*end switch : sign of m in Ylm*/

  }/* endfor : loop over states */


/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlmat_hess(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
         int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
         double *helr,double *heli,
         double *creal,double *cimag,
         double *dhelr,double *dheli,
         double *vnlreal,double *vnlimag,
         double *dvnlreal_x,double *dvnlreal_y,
         double *dvnlreal_z,double *dvnlimag_x,double *dvnlimag_y,
         double *dvnlimag_z,double *dvnlreal_gxgx,double *dvnlreal_gxgy,
         double *dvnlreal_gxgz,double *dvnlreal_gygy,double *dvnlreal_gygz,
         double *dvnlreal_gzgz,double *dvnlimag_gxgx,double *dvnlimag_gxgy,
         double *dvnlimag_gxgz,double *dvnlimag_gygy,double *dvnlimag_gygz,
         double *dvnlimag_gzgz,double xk,double yk,double zk,
         double ylmr,double ylmi,double sgn_l,double *cbyhelr,double *cbyheli)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  int is,i,ipart,iii;
  int ind_loc_c,ioff_v,ind_loc_v;
  int isgn_l ;
  double cre_ind,cim_ind;
  double tmp1,tmp2,tmp3,tmp4;
  double dtmp1,dtmp2;
  double tylmr,tylmi;
  double dtmp3;


/*==========================================================================*/
/* I) Loop over states and particles                                        */
/*     Here helr,heli have vtemp mutliplied in them                         */
/*     Here dhelr,dheli have dvtemp mutliplied in them                      */
/*     See control routine                                                  */

  isgn_l    = (int) (0.5*(sgn_l+1.0)+1.0);
  tylmr     = 2.0*ylmr;
  tylmi     = 2.0*ylmi;

/* storage : atm,state,l,m,n */

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind = creal[ind_loc_c];
    cim_ind = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;

    switch(isgn_l){
      case 1: 
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]    =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]    =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
        }/* endfor : loop over particles */ 
        for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        = -cbyheli[i]*tylmi;
           tmp2        =  cbyheli[i]*tylmr;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
        }/*endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp3       = -cbyhelr[i]*tylmi;
            dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
            dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
            dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
        }/* endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp4       =  cbyhelr[i]*tylmr;
            dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
            dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
            dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
        }/* endfor */

        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = cbyheli[i]*tylmi;
           dvnlreal_gxgx[(i+ioff_v)] += (dtmp3*xk*xk);
           dvnlreal_gxgy[(i+ioff_v)] += (dtmp3*xk*yk);
           dvnlreal_gxgz[(i+ioff_v)] += (dtmp3*xk*zk);
           dvnlreal_gygy[(i+ioff_v)] += (dtmp3*yk*yk);
           dvnlreal_gygz[(i+ioff_v)] += (dtmp3*yk*zk);
           dvnlreal_gzgz[(i+ioff_v)] += (dtmp3*zk*zk);
        }/*endfor */
        for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = -cbyheli[i]*tylmr;
           dvnlimag_gxgx[(i+ioff_v)] += (dtmp3*xk*xk);
           dvnlimag_gxgy[(i+ioff_v)] += (dtmp3*xk*yk);
           dvnlimag_gxgz[(i+ioff_v)] += (dtmp3*xk*zk);
           dvnlimag_gygy[(i+ioff_v)] += (dtmp3*yk*yk);
           dvnlimag_gygz[(i+ioff_v)] += (dtmp3*yk*zk);
           dvnlimag_gzgz[(i+ioff_v)] += (dtmp3*zk*zk);
        }/* endfor */
      break;
      case 2: 
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          cbyhelr[ipart]    =  cre_ind*helr[ipart] - cim_ind*heli[ipart];
          cbyheli[ipart]    =  cre_ind*heli[ipart] + cim_ind*helr[ipart];
        }/* endfor : loop over particles */ 
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp1        =  cbyhelr[i]*tylmr;
           tmp2        =  cbyhelr[i]*tylmi;
           vnlreal[(i+ioff_v)]    += tmp1;
           vnlimag[(i+ioff_v)]    += tmp2;
         }/*endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           tmp3        = -cbyheli[i]*tylmr;
           dvnlreal_x[(i+ioff_v)] += (tmp3*xk);
           dvnlreal_y[(i+ioff_v)] += (tmp3*yk);
           dvnlreal_z[(i+ioff_v)] += (tmp3*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
            tmp4       = -cbyheli[i]*tylmi;
            dvnlimag_x[(i+ioff_v)] += (tmp4*xk);
            dvnlimag_y[(i+ioff_v)] += (tmp4*yk);
            dvnlimag_z[(i+ioff_v)] += (tmp4*zk);
         }/* endfor */
         
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = -cbyhelr[i]*tylmr;
           dvnlreal_gxgx[(i+ioff_v)] += (dtmp3*xk*xk);
           dvnlreal_gxgy[(i+ioff_v)] += (dtmp3*xk*yk);
           dvnlreal_gxgz[(i+ioff_v)] += (dtmp3*xk*zk);
           dvnlreal_gygy[(i+ioff_v)] += (dtmp3*yk*yk);
           dvnlreal_gygz[(i+ioff_v)] += (dtmp3*yk*zk);
           dvnlreal_gzgz[(i+ioff_v)] += (dtmp3*zk*zk);
         }/* endfor */
         for(i=np_nl_rad_str;i<=np_nl;i++){
           dtmp3     = -cbyhelr[i]*tylmi;
           dvnlimag_gxgx[(i+ioff_v)] += (dtmp3*xk*xk);
           dvnlimag_gxgy[(i+ioff_v)] += (dtmp3*xk*yk);
           dvnlimag_gxgz[(i+ioff_v)] += (dtmp3*xk*zk);
           dvnlimag_gygy[(i+ioff_v)] += (dtmp3*yk*yk);
           dvnlimag_gygz[(i+ioff_v)] += (dtmp3*yk*zk);
           dvnlimag_gzgz[(i+ioff_v)] += (dtmp3*zk*zk);
         }/* endfor */
    }/*end switch : sign of m in Ylm*/

  }/* endfor : loop over states */


/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void get_nlmat0(int ncoef,int ismcount,int nstate,
                int np_nlmax,int np_nl_rad_str,int np_nl,int nl_chan_max,
                int irad,
                double *creal,double *cimag,double *vtemp,
                double *vnlreal,double *vnlimag,
                double ylmr)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations                                     */

 int is,ipart,ioff_v,ind_loc_c,ind_loc_v;
 double cre_ind,term;

/*========================================================================*/
/* Loop over the number of states and the number of particles             */

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind   = creal[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
      term  = cre_ind*vtemp[ipart];
      vnlreal[(ipart+ioff_v)] += term;
    }/*endfor*/
  }/* endfor : loop over the number of states */

/*------------------------------------------------------------------------*/
   }/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_pot_pv_fatm(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
                  EWD_SCR *ewd_scr,CPOPTS *cpopts,
                  PSEUDO *pseudo,ATOMMAPS *atommaps,
                  double *cp_enl_ret, int np_nlmax,double *pvten)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm;
  int ind_loc,nl_max,irad,jrad;
  int nl_chan_max;
  double rvol_cp,vol_cp,cp_enl;
  double p11,p22,p33,p12,p13,p23;

/* Local pointers */
  int npart                = clatoms_info->natm_tot;
  double *fx               = clatoms_pos->fx;     
  double *fy               = clatoms_pos->fy;
  double *fz               = clatoms_pos->fz;
  double *hess_xx          = clatoms_pos->hess_xx;
  double *hess_xy          = clatoms_pos->hess_xy;
  double *hess_xz          = clatoms_pos->hess_xz;
  double *hess_yy          = clatoms_pos->hess_yy;
  double *hess_yz          = clatoms_pos->hess_yz;
  double *hess_zz          = clatoms_pos->hess_zz;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;
  double *hmat_cp          = cell->hmat_cp;

  int cp_ptens             = cpopts->cp_ptens_calc;
  int cp_lsda              = cpopts->cp_lsda;
  int atm_hess_calc        = clatoms_info->hess_calc;
  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  double *vpsnorm          = pseudo->vpsnorm; 
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *loc_opt             = pseudo->loc_opt;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int *np_nl               = pseudo->np_nl;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;

  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *dvnlreal_x_up    = cpscr->cpscr_nonloc.dvnlre_x_up;
  double *dvnlreal_y_up    = cpscr->cpscr_nonloc.dvnlre_y_up;
  double *dvnlreal_z_up    = cpscr->cpscr_nonloc.dvnlre_z_up;
  double *dvnlimag_x_up    = cpscr->cpscr_nonloc.dvnlim_x_up;
  double *dvnlimag_y_up    = cpscr->cpscr_nonloc.dvnlim_y_up;
  double *dvnlimag_z_up    = cpscr->cpscr_nonloc.dvnlim_z_up;
  double *dvnlreal_x_dn    = cpscr->cpscr_nonloc.dvnlre_x_dn;
  double *dvnlreal_y_dn    = cpscr->cpscr_nonloc.dvnlre_y_dn;
  double *dvnlreal_z_dn    = cpscr->cpscr_nonloc.dvnlre_z_dn;
  double *dvnlimag_x_dn    = cpscr->cpscr_nonloc.dvnlim_x_dn;
  double *dvnlimag_y_dn    = cpscr->cpscr_nonloc.dvnlim_y_dn;
  double *dvnlimag_z_dn    = cpscr->cpscr_nonloc.dvnlim_z_dn;
  double *dvnlreal_gxgx_up = cpscr->cpscr_nonloc.dvnlre_gxgx_up;
  double *dvnlimag_gxgx_up = cpscr->cpscr_nonloc.dvnlim_gxgx_up;
  double *dvnlreal_gygy_up = cpscr->cpscr_nonloc.dvnlre_gygy_up;
  double *dvnlimag_gygy_up = cpscr->cpscr_nonloc.dvnlim_gygy_up;
  double *dvnlreal_gzgz_up = cpscr->cpscr_nonloc.dvnlre_gzgz_up;
  double *dvnlimag_gzgz_up = cpscr->cpscr_nonloc.dvnlim_gzgz_up;
  double *dvnlreal_gxgy_up = cpscr->cpscr_nonloc.dvnlre_gxgy_up;
  double *dvnlimag_gxgy_up = cpscr->cpscr_nonloc.dvnlim_gxgy_up;
  double *dvnlreal_gygz_up = cpscr->cpscr_nonloc.dvnlre_gygz_up;
  double *dvnlimag_gygz_up = cpscr->cpscr_nonloc.dvnlim_gygz_up;
  double *dvnlreal_gxgz_up = cpscr->cpscr_nonloc.dvnlre_gxgz_up;
  double *dvnlimag_gxgz_up = cpscr->cpscr_nonloc.dvnlim_gxgz_up;
  double *dvnlreal_gxgx_dn = cpscr->cpscr_nonloc.dvnlre_gxgx_dn;
  double *dvnlimag_gxgx_dn = cpscr->cpscr_nonloc.dvnlim_gxgx_dn;
  double *dvnlreal_gygy_dn = cpscr->cpscr_nonloc.dvnlre_gygy_dn;
  double *dvnlimag_gygy_dn = cpscr->cpscr_nonloc.dvnlim_gygy_dn;
  double *dvnlreal_gzgz_dn = cpscr->cpscr_nonloc.dvnlre_gzgz_dn;
  double *dvnlimag_gzgz_dn = cpscr->cpscr_nonloc.dvnlim_gzgz_dn;
  double *dvnlreal_gxgy_dn = cpscr->cpscr_nonloc.dvnlre_gxgy_dn;
  double *dvnlimag_gxgy_dn = cpscr->cpscr_nonloc.dvnlim_gxgy_dn;
  double *dvnlreal_gygz_dn = cpscr->cpscr_nonloc.dvnlre_gygz_dn;
  double *dvnlimag_gygz_dn = cpscr->cpscr_nonloc.dvnlim_gygz_dn;
  double *dvnlreal_gxgz_dn = cpscr->cpscr_nonloc.dvnlre_gxgz_dn;
  double *dvnlimag_gxgz_dn = cpscr->cpscr_nonloc.dvnlim_gxgz_dn;

  double *fxtemp           = ewd_scr->fx;
  double *fytemp           = ewd_scr->fy;
  double *fztemp           = ewd_scr->fz; 
  double *vscr             = ewd_scr->fx2;
  double *vnorm            = ewd_scr->fy2;
  double *vnorm_now        = ewd_scr->fz2;

/*======================================================================*/
/* I) Useful constants                                                  */

  vol_cp   = getdeth(hmat_cp);
  rvol_cp  = 1.0/vol_cp;
  cp_enl   = 0.0;
  nl_max   = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/
  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Loop over the open channels, the states and get the nl potent,   */
/*     pvten and  particle forces                                       */

  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl[lp1]>0){
     for(irad=1;irad<=nrad_max_l[lp1];irad++){
     for(jrad=irad;jrad<=nrad_max_l[lp1];jrad++){
       
/*-----------------------------------------------------------------------*/
/* i) Get the normalization scaled by the volume                        */

     get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,l,
               n_ang_max,irad,jrad,n_rad_max);


      i_shift = l*npart;
      for(ipart=np_nl_rad_str[lp1][jrad];ipart<=np_nl[lp1];ipart++){
        vnorm_now[ipart] = vnorm[ip_nl_rev[(ipart+i_shift)]]*rvol_cp; 
      }/*endfor*/

/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */
    
       sumnl_pot_pv_fatm_hess(npart,nstate_up,np_nlmax,nl_chan_max,
                              np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
                              irad,jrad,
                              ip_nl,vnorm_now,vnlreal_up,vnlimag_up,
                              dvnlreal_gxgx_up,dvnlimag_gxgx_up,
                              dvnlreal_gygy_up,dvnlimag_gygy_up,
                              dvnlreal_gzgz_up,dvnlimag_gzgz_up,
                              dvnlreal_gxgy_up,dvnlimag_gxgy_up,
                              dvnlreal_gxgz_up,dvnlimag_gxgz_up,
                              dvnlreal_gygz_up,dvnlimag_gygz_up,
                              dvnlreal_x_up,dvnlimag_x_up,
                              dvnlreal_y_up,dvnlimag_y_up,
                              dvnlreal_z_up,dvnlimag_z_up,
                              fx,fy,fz,fxtemp,fytemp,fztemp,
                              hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                              atm_hess_calc,cp_ptens,pvten,&cp_enl);
       if(cp_lsda==1){
         sumnl_pot_pv_fatm_hess(npart,nstate_dn,np_nlmax,nl_chan_max,
                                np_nl[lp1],l,np_nl_rad_str[lp1][jrad],
                                irad,jrad,
                                ip_nl,vnorm_now,vnlreal_dn,vnlimag_dn,
                                dvnlreal_gxgx_dn,dvnlimag_gxgx_dn,
                                dvnlreal_gygy_dn,dvnlimag_gygy_dn,
                                dvnlreal_gzgz_dn,dvnlimag_gzgz_dn,
                                dvnlreal_gxgy_dn,dvnlimag_gxgy_dn,
                                dvnlreal_gxgz_dn,dvnlimag_gxgz_dn,
                                dvnlreal_gygz_dn,dvnlimag_gygz_dn,
                                dvnlreal_x_dn,dvnlimag_x_dn,
                                dvnlreal_y_dn,dvnlimag_y_dn,
                                dvnlreal_z_dn,dvnlimag_z_dn,
                                fx,fy,fz,fxtemp,fytemp,fztemp,
                                hess_xx,hess_xy,hess_xz,hess_yy,hess_yz,hess_zz,
                                atm_hess_calc,cp_ptens,pvten,&cp_enl);
       }/*endif*/
     }}/*endfor: radial channels */
    }/*endif: l channel open */
  }/*endfor: l channels     */

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;
  if(cp_ptens==1){
    pvten[1] += cp_enl;
    pvten[5] += cp_enl;
    pvten[9] += cp_enl;
  }/*endif*/

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sumnl_pot_pv_fatm_hess(int npart,int nstate,int np_nlmax,
                            int nl_chan_max,int np_nl,int l,int np_nl_rad_str,
                            int irad, int jrad,
                            int *ip_nl,double *vnorm_now,
                            double *vnlreal,double *vnlimag,
                            double *dvnlreal_gxgx,double *dvnlimag_gxgx,
                            double *dvnlreal_gygy,double *dvnlimag_gygy,
                            double *dvnlreal_gzgz,double *dvnlimag_gzgz,
                            double *dvnlreal_gxgy,double *dvnlimag_gxgy,
                            double *dvnlreal_gxgz,double *dvnlimag_gxgz,
                            double *dvnlreal_gygz,double *dvnlimag_gygz,
                            double *dvnlreal_x,double *dvnlimag_x,
                            double *dvnlreal_y,double *dvnlimag_y,
                            double *dvnlreal_z,double *dvnlimag_z,
                            double *fx,double *fy,double *fz,
                            double *fxtemp,double *fytemp,double *fztemp,
                            double *hess_xx,double *hess_xy,double *hess_xz,
                            double *hess_yy,double *hess_yz,double *hess_zz,
                            int atm_hess_calc,
                            int cp_ptens,double *pvten, double *cp_enl_ret)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,m,i_shift,ipart,is,iii,ioff,ltemp;
  int ind_loc,ind_lm;
  int ioff_i,ioff_j,ind_loc_i,ind_loc_j;
  int hess_ind;
  double p11,p22,p33,p12,p13,p23,cp_enl_now;

  /* Local pointers */
  double cp_enl = *cp_enl_ret;

/*==========================================================================*/
/* I) Loop over the 2*l+1 directions and sum the nl contributions           */

    for(m = 1;m<=(2*l+1);m++){
      ind_lm = m + l*l;
      for(is=1;is<=nstate;is++){

/*----------------------------------------------------------------------*/
/*  i) Get the contrib to the non-local energy                          */

        cp_enl_now = 0.0;
        ioff_i = (is-1)*np_nlmax + 
               +(ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
        ioff_j = (is-1)*np_nlmax + 
               +(ind_lm-1)*nstate*np_nlmax
               +(jrad-1)*nl_chan_max*nstate*np_nlmax;
        if(irad==jrad){
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ind_loc = ipart + ioff_i;
           cp_enl_now = (vnlreal[ind_loc]*vnlreal[ind_loc]
                      +vnlimag[ind_loc]*vnlimag[ind_loc])*vnorm_now[ipart];

           cp_enl += cp_enl_now;
           
         }/*endfor*/

       }else{
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ind_loc_i = ipart + ioff_i;
           ind_loc_j = ipart + ioff_j;
           cp_enl += 2.0*(vnlreal[ind_loc_i]*vnlreal[ind_loc_j]
                         +vnlimag[ind_loc_i]*vnlimag[ind_loc_j])
                         *vnorm_now[ipart];

         }/*endfor*/
       }/*endif*/

/*----------------------------------------------------------------------*/
/* ii) Get the contrib to non-local piece of the pressure tensor        */

        if(cp_ptens==1){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                      +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }/*endfor:ipart*/
         }else{
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            p11 = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p22 = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            p33 = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p12 = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p13 = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i])
                      *vnorm_now[ipart];
            p23 = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                      +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                      +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                      +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i])
                     *vnorm_now[ipart];
            pvten[1] += p11;
            pvten[5] += p22;
            pvten[9] += p33;
            pvten[4] += p12;
            pvten[2] += p12;
            pvten[7] += p13;
            pvten[3] += p13;
            pvten[8] += p23;
            pvten[6] += p23;
          }/*endfor:ipart*/
	 }/*endif*/
        }/*endif:cp_ptens on*/

/*----------------------------------------------------------------------*/
/* iii) Sum the non-local particle force                                  */

        if(irad==jrad){
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          ind_loc = ipart + ioff_i;
          fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_x[ind_loc]*vnlimag[ind_loc])
                              *vnorm_now[ipart];
          fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_y[ind_loc]*vnlimag[ind_loc])
                              *vnorm_now[ipart];
          fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc]*vnlreal[ind_loc]
                               +dvnlimag_z[ind_loc]*vnlimag[ind_loc])
                              *vnorm_now[ipart];
         }/*endfor*/
       }else{
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          ind_loc_i = ipart + ioff_i;
          ind_loc_j = ipart + ioff_j;
          fxtemp[ipart] = -2.0*(dvnlreal_x[ind_loc_i]*vnlreal[ind_loc_j]
                               +dvnlreal_x[ind_loc_j]*vnlreal[ind_loc_i]
                               +dvnlimag_x[ind_loc_i]*vnlimag[ind_loc_j]
                               +dvnlimag_x[ind_loc_j]*vnlimag[ind_loc_i])
                              *vnorm_now[ipart];
          fytemp[ipart] = -2.0*(dvnlreal_y[ind_loc_i]*vnlreal[ind_loc_j]
                               +dvnlreal_y[ind_loc_j]*vnlreal[ind_loc_i]
                               +dvnlimag_y[ind_loc_i]*vnlimag[ind_loc_j]
                               +dvnlimag_y[ind_loc_j]*vnlimag[ind_loc_i])
                              *vnorm_now[ipart];
          fztemp[ipart] = -2.0*(dvnlreal_z[ind_loc_i]*vnlreal[ind_loc_j]
                               +dvnlreal_z[ind_loc_j]*vnlreal[ind_loc_i]
                               +dvnlimag_z[ind_loc_i]*vnlimag[ind_loc_j]
                               +dvnlimag_z[ind_loc_j]*vnlimag[ind_loc_i])
                              *vnorm_now[ipart];
         }/*endfor*/
        }/*endif*/
        i_shift = l*npart;
        for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
          ltemp = ip_nl[(ipart+i_shift)];
          fx[ltemp] += fxtemp[ipart];
          fy[ltemp] += fytemp[ipart];
          fz[ltemp] += fztemp[ipart];
        }/*endfor:atomic forces*/

/*----------------------------------------------------------------------*/
/* iv) Sum the non-local atomic hessian                                  */

        if(atm_hess_calc == 3){
         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgx[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_x[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_x[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gxgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_x[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_x[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }/* endfor */
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gxgx[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgx[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgx[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgx[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_x[ind_loc_j]
	         		     +dvnlimag_x[ind_loc_i]*dvnlimag_x[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gxgy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_y[ind_loc_j]
				     +dvnlimag_x[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gxgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gxgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gxgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gxgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_x[ind_loc_i]*dvnlreal_z[ind_loc_j]
				     +dvnlimag_x[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }/* endfor */
         }/* endif */
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_xx[hess_ind] += fxtemp[ipart];
           hess_xy[hess_ind] += fytemp[ipart];
           hess_xz[hess_ind] += fztemp[ipart];
         }/*endfor:atomic forces*/

         if(irad==jrad){
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc = ipart + ioff_i;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygy[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_y[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_y[ind_loc])
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gygz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_y[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_y[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc]*vnlreal[ind_loc]
                                +dvnlimag_gzgz[ind_loc]*vnlimag[ind_loc]
                                +dvnlreal_z[ind_loc]*dvnlreal_z[ind_loc]
                                +dvnlimag_z[ind_loc]*dvnlimag_z[ind_loc])
                                *vnorm_now[ipart];
          }/* endfor */
         } else {
          for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
            ind_loc_i = ipart + ioff_i;
            ind_loc_j = ipart + ioff_j;
            fxtemp[ipart] = 2.0*(dvnlreal_gygy[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygy[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygy[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygy[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_y[ind_loc_j]
				     +dvnlimag_y[ind_loc_i]*dvnlimag_y[ind_loc_j]))
                                *vnorm_now[ipart];
            fytemp[ipart] = 2.0*(dvnlreal_gygz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gygz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gygz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gygz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_y[ind_loc_i]*dvnlreal_z[ind_loc_j]
				     +dvnlimag_y[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
            fztemp[ipart] = 2.0*(dvnlreal_gzgz[ind_loc_i]*vnlreal[ind_loc_j]
                                +dvnlreal_gzgz[ind_loc_j]*vnlreal[ind_loc_i]
                                +dvnlimag_gzgz[ind_loc_i]*vnlimag[ind_loc_j]
                                +dvnlimag_gzgz[ind_loc_j]*vnlimag[ind_loc_i]
                                +2.0*(dvnlreal_z[ind_loc_i]*dvnlreal_z[ind_loc_j]
				     +dvnlimag_z[ind_loc_i]*dvnlimag_z[ind_loc_j]))
                                *vnorm_now[ipart];
          }/* endfor */
         }/* endif */
         i_shift = l*npart;
         for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
           ltemp = ip_nl[(ipart+i_shift)];
           hess_ind = (ltemp-1)*npart + ltemp;
           hess_yy[hess_ind] += fxtemp[ipart];
           hess_yz[hess_ind] += fytemp[ipart];
           hess_zz[hess_ind] += fztemp[ipart];
         }/*endfor:atomic forces*/

        }/* endif: atm hess calc */

      }/*endfor:loop over states*/
    }/*endfor: loop over the m channels */

/*==========================================================================*/
/* II) Set the return values                                               */

 *cp_enl_ret = cp_enl;

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_fcoef(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                  CPSCR *cpscr,EWD_SCR *ewd_scr,CPOPTS *cpopts,
                  PSEUDO *pseudo,CPEWALD *cpewald,ATOMMAPS *atommaps,
                  CELL *cell, int np_nlmax,double *pvten,FOR_SCR *for_scr)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ipart,nl_max,iii,i,l,istr,irad,jrad,iatm;
  int icount,ismcount,nl_chan_max;
  double vol_cp;
  double tpi,fpi;
  double atemp,btemp,ctemp,aka,akb,akc,arg,xk,yk,zk;
  double xtemp,ytemp,ztemp;
  double g,g2,temp;
  double ylmr[2];
  YLM_CONS ylm_cons;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;

/* Local pointers */
  int npart          = clatoms_info->natm_tot;
  double *x          = clatoms_pos->x;           
  double *y          = clatoms_pos->y;
  double *z          = clatoms_pos->z;
  int natm_typ       = atommaps->natm_typ;
  int *iatm_typ      = atommaps->iatm_atm_typ;
  int *iatm_typ_nl   = atommaps->iatm_atm_typ_nl;
  int *index_atm     = for_scr->iexcl;
  double *cp_box_center     = cell->cp_box_center;
  double *cp_box_center_rel = cell->cp_box_center_rel;

  double *hmat_cp    = cell->hmat_cp;
  double *hmati_cp   = cell->hmati_cp;

  double *cossc      = ewd_scr->cossc;  
  double *sinsc      = ewd_scr->sinsc;
  double *helr       = ewd_scr->helr;   
  double *heli       = ewd_scr->heli;
  double *ewd_scr_x  = ewd_scr->x;
  double *ewd_scr_y  = ewd_scr->y;
  double *ewd_scr_z  = ewd_scr->z;
  double *vtemp      = ewd_scr->q;
  double *vtemp_now  = ewd_scr->vtemp_now;
  double *vscr       = ewd_scr->fx2;
  double *vnorm      = ewd_scr->fy2;
  double *vnorm_now  = ewd_scr->fz2;

  int *kastore_sm    = cpewald->kastr_sm;
  int *kbstore_sm    = cpewald->kbstr_sm;
  int *kcstore_sm    = cpewald->kcstr_sm;
  int *ibreak1_sm    = cpewald->ibrk1_sm;
  int *ibreak2_sm    = cpewald->ibrk2_sm;
  int nktot_sm       = cpewald->nktot_sm;

  int n_ang_max      = pseudo->n_ang_max;
  int n_ang_max_kb   = pseudo->n_ang_max_kb;
  int *ip_nl         = pseudo->ip_nl;
  int *ip_nl_rev        = pseudo->ip_nl_rev;
  int  np_nonloc_cp_box_kb = pseudo->np_nonloc_cp_box_kb;
  int *np_nl         = pseudo->np_nl;
  double *gzvps0     = pseudo->gzvps0;
  double *vpsnorm    = pseudo->vpsnorm;
  int n_rad_max      = pseudo->n_rad_max;
  int *nrad_max_l    = pseudo->nrad_max_l;
  int **np_nl_rad_str= pseudo->np_nl_rad_str;

  int ncoef             = cpcoeffs_info->ncoef;
  int cp_lsda           = cpopts->cp_lsda;
  int cp_hess_calc      = cpopts->cp_hess_calc;
  int nstate_up         = cpcoeffs_info->nstate_up_proc;
  int nstate_dn         = cpcoeffs_info->nstate_dn_proc;
  double *fcreal_up     = cpcoeffs_pos->fcre_up;  
  double *fcimag_up     = cpcoeffs_pos->fcim_up;
  double *fcreal_dn     = cpcoeffs_pos->fcre_dn;
  double *fcimag_dn     = cpcoeffs_pos->fcim_dn;
  double *cp_hess_re_up = cpcoeffs_pos->cp_hess_re_up;
  double *cp_hess_im_up = cpcoeffs_pos->cp_hess_im_up;
  double *cp_hess_re_dn = cpcoeffs_pos->cp_hess_re_dn;
  double *cp_hess_im_dn = cpcoeffs_pos->cp_hess_im_dn;
  double *vnlreal_up    = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up    = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn    = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn    = cpscr->cpscr_nonloc.vnlim_dn;


  double *hmat_big          = cell->hmat;
  double *hmati_big         = cell->hmati;

/*======================================================================*/
/* 0) Get some useful constants                                         */

  vol_cp              = getdeth(hmat_cp);
  tpi                 = 2.0*M_PI;

  fpi                 = 4.0*M_PI;
  ylm_cons.rt_fpi     = 1.0/sqrt(fpi);
  ylm_cons.rt_thrfpi  = sqrt(3.0/fpi);
  ylm_cons.rt_threpi  = sqrt(1.50/fpi);
  ylm_cons.hrt_fivfpi = 0.50*sqrt(5.0/fpi);
  ylm_cons.rt_fiftepi = sqrt(7.50/fpi);
  ylm_cons.hrt_sevfpi = 0.50*sqrt(7.0/fpi);
  ylm_cons.hrt_toepi  = 0.50*sqrt(10.50/fpi)/sqrt(2.0);
  ylm_cons.hrt_ohffpi = 0.50*sqrt(105.0/fpi)/sqrt(2.0);
  ylm_cons.hrt_tfepi  = 0.50*sqrt(17.50/fpi)/sqrt(2.0);



/*======================================================================*/
/* I) Find the open channels                                            */

  nl_max = -1;
  for(i=1;i<=(n_ang_max_kb+1);i++){
   if(np_nl[i]>0){nl_max=i-1;}
  }/*endfor*/

  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Find cos and sin of sc components of charged particles           */
/*  ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                         */

  for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
    iatm = ip_nl[ipart];
    dx  = x[iatm] - cp_box_center[1];
    dy  = y[iatm] - cp_box_center[2];
    dz  = z[iatm] - cp_box_center[3];
    asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
    asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
    asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
    sx  = asx - NINT(asx);
    sy  = asy - NINT(asy);
    sz  = asz - NINT(asz);
    dx  = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
    dy  = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
    dz  = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];
    xtemp = dx - cp_box_center_rel[1];
    ytemp = dy - cp_box_center_rel[2];
    ztemp = dz - cp_box_center_rel[3];

    ewd_scr_x[ipart] = xtemp*hmati_cp[1]
                     + ytemp*hmati_cp[4]
                     + ztemp*hmati_cp[7];
    ewd_scr_y[ipart] = xtemp*hmati_cp[2]
                     + ytemp*hmati_cp[5]
                     + ztemp*hmati_cp[8];
    ewd_scr_z[ipart] = xtemp*hmati_cp[3]
                     + ytemp*hmati_cp[6]
                     + ztemp*hmati_cp[9];
    ctemp = ewd_scr_z[ipart]*tpi;
    cossc[ipart] = cos(ctemp);
    sinsc[ipart] = sin(ctemp);
  }/*endfor*/

/*======================================================================*/
/* III) Perform the recip space sum                                     */

  for(icount=1;icount<=nktot_sm;icount++){

/*----------------------------------------------------------------------*/
/* i) Get the k vectors                                                 */

    aka = (double)(kastore_sm[icount]);
    akb = (double)(kbstore_sm[icount]);
    akc = (double)(kcstore_sm[icount]);
    xk = (aka*hmati_cp[1]+akb*hmati_cp[2]+akc*hmati_cp[3])*tpi;
    yk = (aka*hmati_cp[4]+akb*hmati_cp[5]+akc*hmati_cp[6])*tpi;
    zk = (aka*hmati_cp[7]+akb*hmati_cp[8]+akc*hmati_cp[9])*tpi;
    g2 = xk*xk+yk*yk+zk*zk;
    g  = sqrt(g2);

/*----------------------------------------------------------------------*/
/* ii) If break point number one calculate the initial helpful vectors   */

    if(ibreak1_sm[icount]==1){
      for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
        atemp = ewd_scr_x[ipart];
        btemp = ewd_scr_y[ipart];
        ctemp = ewd_scr_z[ipart];
        arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
        helr[ipart] = cos(arg);
        heli[ipart] = sin(arg);
      }/*endfor*/
    }/*endif*/

/*----------------------------------------------------------------------*/
/* iii) Get the external non-local force on the coefs                   */
/*      use structure factor and wavefunction to compute                */
/*      the array znl and its derivative dznl for the                   */
/*      calculation of the nonlocal forces on the coefs                 */
/*      (done separately for each angular momentum component)           */

    ismcount = icount;
    control_nlfcoef(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
                    cpscr,cpopts,pseudo,ewd_scr,atommaps,
                    np_nlmax,nl_max,index_atm,
                    g,xk,yk,zk,vol_cp,ismcount,&ylm_cons);

/*----------------------------------------------------------------------*/
/* iv) If break point two, increment the helpful vectors                */

    if(ibreak2_sm[icount]==1){
      for(ipart=1;ipart<=np_nonloc_cp_box_kb;ipart++){
        temp = helr[ipart];
        helr[ipart] = helr[ipart]*cossc[ipart] 
                    - heli[ipart]*sinsc[ipart];
        heli[ipart] = heli[ipart]*cossc[ipart] 
                    + temp*sinsc[ipart];
      }/*endfor:ipart*/

    }/*endif:ibreak2_sm==1*/
  }/*endfor:icount*/

/*======================================================================*/
/* IV) g=0 term (non-local potential,l=0 only!!!!)                      */

  if(np_nl[1]>0){

    l = 0;
    ylmr[1] = 1.0/sqrt(fpi);
    ismcount = nktot_sm + 1;
    for(irad=1;irad<=nrad_max_l[1];irad++){
      for(ipart=1;ipart<=npart;ipart++){
        iii = n_rad_max*(iatm_typ[ipart]-1) + irad;
        vtemp[ipart] = gzvps0[iii];
      }/*endfor*/
      for(ipart=np_nl_rad_str[1][irad];ipart<=np_nl[1];ipart++){
        vtemp_now[ipart] = vtemp[ip_nl[ipart]];
      }/*endfor*/
      for(jrad=1;jrad<=nrad_max_l[1];jrad++){
        get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ,natm_typ,npart,l,
                    n_ang_max,irad,jrad,n_rad_max);
        istr = MAX(np_nl_rad_str[1][irad],np_nl_rad_str[1][jrad]);
        for(ipart=istr;ipart<=np_nl[1];ipart++){
          vnorm_now[ipart] = vnorm[ip_nl[ipart]];
	}/*endfor*/
        get_nlfor0(ncoef,ismcount,nstate_up,np_nlmax,np_nl[1],
                   istr,nl_chan_max,jrad,vtemp_now,vnorm_now,
                   fcreal_up,fcimag_up,vnlreal_up,vnlimag_up,ylmr[1],vol_cp);
#ifdef HESS
        if(cp_hess_calc == 1){
          get_nlhess0(ncoef,ismcount,nstate_up,np_nlmax,np_nl[1],
                      istr,nl_chan_max,jrad,vtemp_now,vnorm_now,cp_hess_re_up,
                      vnlreal_up,vnlimag_up,ylmr[1],vol_cp);
        }/*endif*/
#endif
        if(cp_lsda==1){
          get_nlfor0(ncoef,ismcount,nstate_dn,np_nlmax,np_nl[1],
                     istr,nl_chan_max,jrad,vtemp_now,vnorm_now,
                     fcreal_dn,fcimag_dn,vnlreal_dn,vnlimag_dn,ylmr[1],vol_cp);
#ifdef HESS
          if(cp_hess_calc == 1){
           get_nlhess0(ncoef,ismcount,nstate_dn,np_nlmax,np_nl[1],
                       istr,nl_chan_max,jrad,vtemp_now,vnorm_now,cp_hess_re_dn,
                       vnlreal_dn,vnlimag_dn,ylmr[1],vol_cp);
          }/*endif*/
#endif
        }/*endif:lsda*/
      }/*endfor*/
    }/*endfor*/
  }/*endif:l=0 channel open*/

/*======================================================================*/
   }/*end routine*/
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_nlfcoef(CLATOMS_INFO *clatoms_info,
                     CPCOEFFS_INFO *cpcoeffs_info,
                     CPCOEFFS_POS *cpcoeffs_pos,
                     CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                     EWD_SCR *ewd_scr,ATOMMAPS *atommaps,
                     int np_nlmax,int nl_max,int *index_atm,
                     double g,double xk,double yk,double zk,double vol,
                     int ismcount,YLM_CONS *ylm_cons)

/*==========================================================================*/
    {/*Begin Routine */
/*==========================================================================*/
/*     Local Variables */

  int l,itype,ind_lm,m,ipart,i_shift,lp1,ltemp,ltemp2,ktemp;
  int nl_chan_max,istr,irad,jrad;
  double sgn_l,rvol4,tmp1;
  double ylmr[21], ylmi[21];
  double dylmr_gx[21], dylmi_gx[21];
  double dylmr_gy[21], dylmi_gy[21];
  double dylmr_gz[21], dylmi_gz[21];

/* Local pointers */     
  int npart                = clatoms_info->natm_tot;
  int natm_typ             = atommaps->natm_typ;
  int *iatm_typ            = atommaps->iatm_atm_typ;
  int *iatm_typ_nl         = atommaps->iatm_atm_typ_nl;

  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;
  int ncoef                = cpcoeffs_info->ncoef;
  double *fcreal_up        = cpcoeffs_pos->fcre_up;  
  double *fcimag_up        = cpcoeffs_pos->fcim_up;
  double *fcreal_dn        = cpcoeffs_pos->fcre_dn;
  double *fcimag_dn        = cpcoeffs_pos->fcim_dn;
  int cp_lsda              = cpopts->cp_lsda;
  int cp_hess_calc         = cpopts->cp_hess_calc;

  int *np_nl               = pseudo->np_nl;
  int *ip_nl               = pseudo->ip_nl;
  int *ip_nl_rev           = pseudo->ip_nl_rev;
  int  np_nonloc_cp_box_kb = pseudo->np_nonloc_cp_box_kb;
  int nsplin_g             = pseudo->nsplin_g;
  double dg_spl            = pseudo->dg_spl;
  double gmin_spl          = pseudo->gmin_spl;
  double *vps0             = pseudo->vps0;
  double *vps1             = pseudo->vps1;
  double *vps2             = pseudo->vps2;
  double *vps3             = pseudo->vps3;
  double *vpsnorm          = pseudo->vpsnorm;
  int n_ang_max            = pseudo->n_ang_max;
  int n_ang_max_kb         = pseudo->n_ang_max_kb;
  int n_rad_max            = pseudo->n_rad_max;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;


  double *vnlreal_up          = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up          = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn          = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn          = cpscr->cpscr_nonloc.vnlim_dn;
  double *cp_hess_re_up       = cpcoeffs_pos->cp_hess_re_up;
  double *cp_hess_im_up       = cpcoeffs_pos->cp_hess_im_up;
  double *cp_hess_re_dn       = cpcoeffs_pos->cp_hess_re_dn;
  double *cp_hess_im_dn       = cpcoeffs_pos->cp_hess_im_dn;

  double *helr             = ewd_scr->helr;   
  double *heli             = ewd_scr->heli;
  double *helr_now         = ewd_scr->helr_now;
  double *heli_now         = ewd_scr->heli_now;
  double *vtemp            = ewd_scr->q;
  double *vtemp_now        = ewd_scr->vtemp_now;
  double *vnorm            = ewd_scr->fx;
  double *vnorm_now        = ewd_scr->fy;
  double *vscr             = ewd_scr->fz; 
  double *vfactr           = ewd_scr->fx2; 
  double *vfacti           = ewd_scr->fy2; 

/*==========================================================================*/
/* I) Get the ylm and other constants                                       */

  get_ylm(xk,yk,zk,g,ylmr,ylmi,dylmr_gx,dylmi_gx,dylmr_gy,dylmi_gy,
          dylmr_gz,dylmi_gz,ylm_cons);
  rvol4 = 4.0/vol;

/*==========================================================================*/
/* II) Loop over the channels                                               */

    nl_chan_max = (nl_max+1)*(nl_max+1);
    for(l=0;l<=nl_max;l++){
      lp1 = l+1;
      sgn_l = 1.0;
      if(l%2==1)sgn_l = -1.0;
      if(np_nl[lp1]>0){
       for(irad=1;irad<=nrad_max_l[lp1];irad++){
#ifdef DEBUG_GJM_INNER
         printf("Hi! Lang =%d Nrad =%d in nlcoef\n",l,irad);
#endif
/*-----------------------------------------------------------------------*/
/* i) Get the norm and bessel transform of the projector                 */
/*    and scale the particle structure factor                            */

        for(itype=1;itype<=natm_typ;itype++){
          index_atm[itype] =  (itype-1)*nsplin_g*(n_ang_max+1)*n_rad_max
                            + l*nsplin_g*n_rad_max +  (irad-1)*nsplin_g;
        }/*endfor*/

        get_vpsnow(index_atm,nsplin_g,gmin_spl,dg_spl,g,
                   vps0,vps1,vps2,vps3,vtemp,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,
                   np_nl_rad_str[lp1][irad]);

        i_shift = l*npart;
        for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
          ktemp            = ipart+i_shift;
          ltemp            = ip_nl_rev[ktemp];
          vtemp_now[ipart] = vtemp[ltemp];
          helr_now[ipart]  = helr[ltemp];
          heli_now[ipart]  = heli[ltemp];
          tmp1             = vtemp_now[ipart]*rvol4;
          helr_now[ipart]  *=  tmp1;
          heli_now[ipart]  *=  tmp1;
        }/*endfor*/
     
/*-----------------------------------------------------------------------*/
/* i) Loop over the channels and compute the forces                      */

       for(m=1;m<=(2*l+1);m++){
         ind_lm = m + l*l;

         for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl[lp1];ipart++){
           vfactr[ipart]   =  (helr_now[ipart]*ylmr[ind_lm]
                              -heli_now[ipart]*ylmi[ind_lm]);
           vfacti[ipart]   = -(heli_now[ipart]*ylmr[ind_lm]
                              +helr_now[ipart]*ylmi[ind_lm]);
         }/*endfor*/
         for(jrad=1;jrad<=nrad_max_l[lp1];jrad++){
          get_vpsnorm(vscr,vpsnorm,vnorm,iatm_typ_nl,natm_typ,np_nonloc_cp_box_kb,l,
                       n_ang_max,irad,jrad,n_rad_max);

           istr = MAX(np_nl_rad_str[lp1][irad],np_nl_rad_str[lp1][jrad]);

           for(ipart=istr;ipart<=np_nl[lp1];ipart++){
             ltemp            = ip_nl_rev[(ipart+i_shift)];
             vnorm_now[ipart] = vnorm[ltemp];
           }/*endfor*/

           get_nlfor(ncoef,ismcount,nstate_up,ind_lm,
                     np_nlmax,np_nl[lp1],
                     istr,nl_chan_max,jrad,
                     vfactr,vfacti,fcreal_up,fcimag_up,
                     vnlreal_up,vnlimag_up,vnorm_now,
                     ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            if(cp_hess_calc == 1){
              get_nlhess(ncoef,ismcount,nstate_up,ind_lm,
                         np_nlmax,np_nl[lp1],
                         istr,nl_chan_max,jrad,
                         vtemp_now,vnorm_now,cp_hess_re_up,cp_hess_im_up,
                         helr,heli,ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            }/*endif*/
            if(cp_lsda==1){
              get_nlfor(ncoef,ismcount,nstate_dn,ind_lm,
                        np_nlmax,np_nl[lp1],
                        istr,nl_chan_max,jrad,
                        vfactr,vfacti,fcreal_dn,fcimag_dn,
                        vnlreal_dn,vnlimag_dn,vnorm_now,
                        ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            if(cp_hess_calc == 1){
              get_nlhess(ncoef,ismcount,nstate_dn,ind_lm,
                         np_nlmax,np_nl[lp1],
                         istr,nl_chan_max,jrad,
                         vtemp_now,vnorm_now,cp_hess_re_dn,cp_hess_im_dn,
                         helr,heli,ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            }/*endif*/

          }/*endif:lsda on*/
         }/*endfor: loop over radial quantum number */
        }/*endfor:loop over m quantum number*/
       }/*endfor: loop over radial quantum number */
      }/*endif:np_nl>0*/
    }/*endfor:loop over l quantum number*/   

/*======================================================================*/
   }/*end routine*/
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlfor(int ncoef,int ismcount,int nstate,int ind_lm,int np_nlmax,
               int np_nl,int np_nl_rad_str,int nl_chan_max,int jrad,
               double *vfactr,double *vfacti,double *fcreal,double *fcimag,
               double *vnlreal,double *vnlimag,double *vnorm_now,
               double ylmr,double ylmi,
               double sgn_l,double vol)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  int is,ipart,ind_loc_c,ind_loc_v,ioff_v;
  double tmp1,tmp2;

/*==========================================================================*/
/* Loop over the number of states and the number of particles               */
/*     Here helr,heli have vtemp*vnorm*4.0/vol mutliplied in them           */
/*     See control routine */

    for(is=1;is<=nstate;is++){
      ioff_v = (is-1)*np_nlmax + 
               +(ind_lm-1)*nstate*np_nlmax
               +(jrad-1)*nl_chan_max*nstate*np_nlmax;
      ind_loc_c  = ismcount + ncoef*(is-1);

      for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
        ind_loc_v = ipart + ioff_v;
        tmp1      = -(vfactr[ipart]*vnlreal[ind_loc_v]
                     -vfacti[ipart]*vnlimag[ind_loc_v])*vnorm_now[ipart];
        tmp2      = -(vfacti[ipart]*vnlreal[ind_loc_v]
                     +vfactr[ipart]*vnlimag[ind_loc_v])*vnorm_now[ipart];
        fcreal[ind_loc_c] += tmp1;
        fcimag[ind_loc_c] += tmp2;
      }/* endfor ipart */

    }/* endfor : loop over states*/


/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlhess(int ncoef,int ismcount,int nstate,int ind_lm,
                int np_nlmax,int np_nl,int istr,int nl_chan_max,int jrad,
                double *vtemp,double *vnorm,double *cp_hess_re,double *cp_hess_im,
                double *helr,double *heli,double ylmr,double ylmi,
                double sgn_l,double vol)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  double rvol;
  int is,ipart,ind_loc_c,ind_loc_v,ioff_v;
  double hfactrr,hfactir;

/*==========================================================================*/
/* Loop over the number of states and the number of particles               */

  rvol = 1.0/vol;


    for(ipart=istr;ipart<=np_nl;ipart++){
      hfactrr = (ylmr*helr[ipart] - ylmi*heli[ipart])*vtemp[ipart];
      hfactir = (ylmr*heli[ipart] + ylmr*helr[ipart])*vtemp[ipart];
      cp_hess_re[ismcount] += 4.0*(hfactrr*hfactrr + hfactir*hfactir)
	                    *vnorm[ipart]*rvol*(1.0+sgn_l);
      cp_hess_im[ismcount] += 4.0*(hfactrr*hfactrr + hfactir*hfactir)
	                    *vnorm[ipart]*rvol*(1.0+sgn_l);
    }/* endfor : loop over particles */ 

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlfor0(int ncoef,int ismcount,int nstate,
                int np_nlmax,int np_nl,
                int istr,int nl_chan_max,int jrad,
                double *vtemp,double *vnorm,double *fcreal,double *fcimag,
                double *vnlreal,double *vnlimag,double ylmr,double vol)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations                                     */

   int is,ipart,ind_loc_c,ind_loc_v,ioff_v;
   double vfactr;
   double rvol;

/*========================================================================*/
/* Loop over the number of states and the number of particles             */

  rvol = 1.0/vol;
  for(is=1;is<=nstate;is++){

    ind_loc_c = ismcount + (is-1)*ncoef;
    ioff_v = (is-1)*np_nlmax  
            +(jrad-1)*nl_chan_max*nstate*np_nlmax;
    for(ipart=istr;ipart<=np_nl;ipart++){
      ind_loc_v = ipart + ioff_v;
      vfactr  =  ylmr*vtemp[ipart]*vnorm[ipart]*rvol;
      fcreal[ind_loc_c] -= 2.0*vfactr*vnlreal[ind_loc_v];
    }/* endfor : loop over the number of particles*/

  }/* endfor : loop over the number of states */

/*------------------------------------------------------------------------*/
  }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlhess0(int ncoef,int ismcount,int nstate,
                int np_nlmax,int np_nl,
                int istr,int nl_chan_max,int jrad,
                double *vtemp,double *vnorm,double *cp_hess,
                double *vnlreal,double *vnlimag,double ylmr,double vol)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations                                     */

   int is,ipart,ind_loc_c,ind_loc_v;
   double hfactrr,rvol;

/*========================================================================*/
/* Loop over the number of states and the number of particles             */

  rvol = 1.0/vol;
  for(is=1;is<=nstate;is++){

    for(ipart=istr;ipart<=np_nl;ipart++){
      hfactrr = ylmr*vtemp[ipart]*vnorm[ipart]*rvol;
      cp_hess[ismcount] += 2.0*hfactrr*hfactrr;
    }/* endfor : loop over particles */ 

  }/* endfor : loop over the number of states */

/*------------------------------------------------------------------------*/
   }/*end routine*/
/*========================================================================*/




/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void get_vpsnorm(double *vnorm_typ,double *vpsnorm,double *vnorm,int *iatm_typ,
                 int natm_typ,int npart,int l,int n_ang_max,
                 int irad, int jrad, int n_rad_max)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*        Local variable declarations                                     */

  int ipart,itype,ind_ityp,n_rad_max_sq;

/*========================================================================*/
/* Find norm of pseudopotential                                           */

  n_rad_max_sq = n_rad_max*n_rad_max;
  for(itype=1;itype<=natm_typ;itype++){
    ind_ityp = (itype-1)*(n_ang_max+1)*n_rad_max_sq
             + l*n_rad_max_sq + (irad-1)*n_rad_max + jrad;
    vnorm_typ[itype] = vpsnorm[ind_ityp];
  }/*endfor*/

  for(ipart=1;ipart<=npart;ipart++){
    vnorm[ipart] = vnorm_typ[iatm_typ[ipart]];
  }/*endfor*/

/*------------------------------------------------------------------------*/
  }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_ylm(double xk,double yk,double zk,double g,
             double *ylmr,double *ylmi,double *dylmr_x,double *dylmi_x,
             double *dylmr_y,double *dylmi_y,double *dylmr_z,double *dylmi_z,
             YLM_CONS *ylm_cons)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  double y00r,y00i,y10r,y10i,y11r,y11i,y20r,y20i,y21r,y21i;
  double y22r,y22i,y30r,y30i,y31r,y31i,y32r,y32i,y33r,y33i;
  double dy00r_phi,dy00i_phi,dy10r_phi,dy10i_phi,dy11r_phi,dy11i_phi;
  double dy20r_phi,dy20i_phi,dy21r_phi,dy21i_phi,dy22r_phi,dy22i_phi;
  double dy30r_phi,dy30i_phi,dy31r_phi,dy31i_phi,dy32r_phi,dy32i_phi;
  double dy33r_phi,dy33i_phi;
  double dy00r_theta,dy00i_theta,dy10r_theta,dy10i_theta;
  double dy11r_theta,dy11i_theta,dy20r_theta,dy20i_theta;
  double dy21r_theta,dy21i_theta,dy22r_theta,dy22i_theta;
  double dy30r_theta,dy30i_theta,dy31r_theta,dy31i_theta;
  double dy32r_theta,dy32i_theta,dy33r_theta,dy33i_theta;
  double ctheta,stheta,cphi,sphi,c2phi,s2phi,c3phi,s3phi;
  double xydist;
  double pre_phi_x,pre_phi_y,pre_theta_x,pre_theta_y,pre_theta_z;
  int iii;

 /* Local pointers */
  double rt_fpi     = ylm_cons->rt_fpi;
  double rt_thrfpi  = ylm_cons->rt_thrfpi;
  double rt_threpi  = ylm_cons->rt_threpi;
  double hrt_fivfpi = ylm_cons->hrt_fivfpi;
  double rt_fiftepi = ylm_cons->rt_fiftepi;
  double hrt_sevfpi = ylm_cons->hrt_sevfpi;
  double hrt_toepi  = ylm_cons->hrt_toepi;
  double hrt_ohffpi = ylm_cons->hrt_ohffpi;
  double hrt_tfepi  = ylm_cons->hrt_tfepi;

/*==========================================================================*/
/* I) Calculate polar angles of the vector g                                */

  ctheta = zk/g;
  stheta = sqrt(xk*xk + yk*yk)/g;
  xydist = sqrt(xk*xk + yk*yk);
  if(xydist==0){
    cphi = 1.0;
    sphi = 0.0;
  }else{
    cphi = xk/xydist;
    sphi = yk/xydist;
  }/*endif*/
  c2phi = cphi*cphi - sphi*sphi;
  s2phi = 2.0*cphi*sphi;
  c3phi = cphi*c2phi - sphi*s2phi;
  s3phi = cphi*s2phi + c2phi*sphi;

/*==========================================================================*/
/* I.i) l=0                                                                 */

  y00r    = rt_fpi;
  y00i    = 0.0;
  ylmr[1] = y00r;
  ylmi[1] = y00i;

/*==========================================================================*/
/* I.ii) l=1 (phi derivatives have stheta divided out)                      */

  y10r        =  rt_thrfpi*ctheta;
  dy10r_theta = -rt_thrfpi*stheta;
  y10i        = 0.0;
  y11r        =  rt_threpi*stheta*cphi;
  dy11r_theta =  rt_threpi*ctheta*cphi;
  dy11r_phi   = -rt_threpi*sphi;
  y11i        =  rt_threpi*stheta*sphi;
  dy11i_theta =  rt_threpi*ctheta*sphi;
  dy11i_phi   =  rt_threpi*cphi;
  ylmr[2] = y10r;
  ylmi[2] = y10i;
  ylmr[3] = y11r;
  ylmi[3] = y11i;
  ylmr[4] = y11r;
  ylmi[4] = -y11i;

/*==========================================================================*/
/* I.iii) l=2 (phi derivatives have stheta divided out)                     */

  y20r    = hrt_fivfpi*(3.0*ctheta*ctheta - 1.0);
  dy20r_theta = -hrt_fivfpi*6.0*ctheta*stheta;
  y20i    = 0.0;
  y21r    = rt_fiftepi*stheta*ctheta*cphi;
  dy21r_theta = rt_fiftepi*cphi*(ctheta*ctheta-stheta*stheta);
  dy21r_phi   = -rt_fiftepi*ctheta*sphi;
  y21i    = rt_fiftepi*stheta*ctheta*sphi;
  dy21i_theta = rt_fiftepi*sphi*(ctheta*ctheta-stheta*stheta);
  dy21i_phi   = rt_fiftepi*ctheta*cphi;
  y22r    = 0.50*rt_fiftepi*stheta*stheta*c2phi;
  dy22r_theta = rt_fiftepi*ctheta*stheta*c2phi;
  dy22r_phi   = -rt_fiftepi*stheta*s2phi;
  y22i    = 0.50*rt_fiftepi*stheta*stheta*s2phi;
  dy22i_theta = rt_fiftepi*ctheta*stheta*s2phi;
  dy22i_phi   = rt_fiftepi*stheta*c2phi;
  ylmr[5] = y20r;
  ylmi[5] = y20i;
  ylmr[6] = y21r;
  ylmi[6] = y21i;
  ylmr[7] = y21r;
  ylmi[7] = -y21i;
  ylmr[8] = y22r;
  ylmi[8] = y22i;
  ylmr[9] = y22r;
  ylmi[9] = -y22i;

/*==========================================================================*/
/* I.iv) l=3 (phi derivatives have stheta divided out)                      */

  y30r        = hrt_sevfpi*(5.0*ctheta*ctheta*ctheta - 3.0*ctheta);
  dy30r_theta = -hrt_sevfpi*stheta*(15.0*ctheta*ctheta-3.0);
  y30i        = 0.0;

  y31r        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*cphi;
  dy31r_theta = hrt_toepi*(ctheta*(5.0*ctheta*ctheta - 1.0)
                - stheta*stheta*10.0*ctheta)*cphi;
  dy31r_phi   = -hrt_toepi*(5.0*ctheta*ctheta - 1.0)*sphi;
  y31i        = hrt_toepi*stheta*(5.0*ctheta*ctheta - 1.0)*sphi;
  dy31i_theta = hrt_toepi*(ctheta*(5.0*ctheta*ctheta - 1.0)
                - stheta*stheta*10.0*ctheta)*sphi;
  dy31i_phi   = hrt_toepi*(5.0*ctheta*ctheta - 1.0)*cphi;
  y32r        = hrt_ohffpi*stheta*stheta*ctheta*c2phi;
  dy32r_theta = hrt_ohffpi*(2.0*ctheta*stheta*ctheta
                - stheta*stheta*stheta)*c2phi;
  dy32r_phi   = -2.0*hrt_ohffpi*stheta*ctheta*s2phi;
  y32i        = hrt_ohffpi*stheta*stheta*ctheta*s2phi;
  dy32i_theta = hrt_ohffpi*(2.0*ctheta*stheta*ctheta
                - stheta*stheta*stheta)*s2phi;
  dy32i_phi   = 2.0*hrt_ohffpi*ctheta*stheta*c2phi;
  y33r        = hrt_tfepi*stheta*stheta*stheta*c3phi;
  dy33r_theta = 3.0*hrt_tfepi*ctheta*stheta*stheta*c3phi;
  dy33r_phi   = -3.0*hrt_tfepi*stheta*stheta*s3phi;
  y33i        = hrt_tfepi*stheta*stheta*stheta*s3phi;
  dy33i_theta = 3.0*hrt_tfepi*ctheta*stheta*stheta*s3phi;
  dy33i_phi   = 3.0*hrt_tfepi*stheta*stheta*c3phi;

  ylmr[10] = y30r;
  ylmi[10] = y30i;
  ylmr[11] = y31r;
  ylmi[11] = y31i;
  ylmr[12] = y31r;
  ylmi[12] = -y31i;
  ylmr[13] = y32r;
  ylmi[13] = y32i;
  ylmr[14] = y32r;
  ylmi[14] = -y32i;
  ylmr[15] = y33r;
  ylmi[15] = y33i;
  ylmr[16] = y33r;
  ylmi[16] = -y33i;

/*==========================================================================*/
/* II) Metrix tensor factors to get gradients of ylm's                      */

  pre_phi_x   = -sphi/g;
  pre_phi_y   =  cphi/g;
  pre_theta_x =  ctheta*cphi/g;
  pre_theta_y =  ctheta*sphi/g;
  pre_theta_z =  -stheta/g;

/*==========================================================================*/
/* II.i) grad l=0                                                           */

  dylmr_x[1] = 0.00;
  dylmi_x[1] = 0.00;
  dylmr_y[1] = 0.00;
  dylmi_y[1] = 0.00;
  dylmr_z[1] = 0.00;
  dylmi_z[1] = 0.00;

/*==========================================================================*/
/* II.ii) grad l=1                                                          */


/*--------------------------------------------------------------------------*/
/* m=0                                                                      */

  dylmr_x[2] = dy10r_theta*pre_theta_x;
  dylmi_x[2] = 0.0;
  dylmr_y[2] = dy10r_theta*pre_theta_y;
  dylmi_y[2] = 0.0;
  dylmr_z[2] = dy10r_theta*pre_theta_z;
  dylmi_z[2] = 0.0;

/*--------------------------------------------------------------------------*/
/* m=1,-1                                                                   */

  dylmr_x[3] = dy11r_phi*pre_phi_x + dy11r_theta*pre_theta_x;
  dylmi_x[3] = dy11i_phi*pre_phi_x + dy11i_theta*pre_theta_x;
  dylmr_y[3] = dy11r_phi*pre_phi_y + dy11r_theta*pre_theta_y;
  dylmi_y[3] = dy11i_phi*pre_phi_y + dy11i_theta*pre_theta_y;
  dylmr_z[3] =         + dy11r_theta*pre_theta_z;
  dylmi_z[3] =         + dy11i_theta*pre_theta_z;
  dylmr_x[4] =  dylmr_x[3];
  dylmi_x[4] = -dylmi_x[3];
  dylmr_y[4] =  dylmr_y[3];
  dylmi_y[4] = -dylmi_y[3];
  dylmr_z[4] =  dylmr_z[3];
  dylmi_z[4] = -dylmi_z[3];

/*==========================================================================*/
/* II.iii) grad l=2                                                         */

/*--------------------------------------------------------------------------*/
/* m=0                                                                      */

  dylmr_x[5] = dy20r_theta*pre_theta_x;
  dylmi_x[5] = 0.0;
  dylmr_y[5] = dy20r_theta*pre_theta_y;
  dylmi_y[5] = 0.0;
  dylmr_z[5] = dy20r_theta*pre_theta_z;
  dylmi_z[5] = 0.0;

/*--------------------------------------------------------------------------*/
/* m=1,-1                                                                   */

  dylmr_x[6] = dy21r_phi*pre_phi_x + dy21r_theta*pre_theta_x;
  dylmi_x[6] = dy21i_phi*pre_phi_x + dy21i_theta*pre_theta_x;
  dylmr_y[6] = dy21r_phi*pre_phi_y + dy21r_theta*pre_theta_y;
  dylmi_y[6] = dy21i_phi*pre_phi_y + dy21i_theta*pre_theta_y;
  dylmr_z[6] =                     + dy21r_theta*pre_theta_z;
  dylmi_z[6] =                     + dy21i_theta*pre_theta_z;
  dylmr_x[7] =  dylmr_x[6];
  dylmi_x[7] = -dylmi_x[6];
  dylmr_y[7] =  dylmr_y[6];
  dylmi_y[7] = -dylmi_y[6];
  dylmr_z[7] =  dylmr_z[6];
  dylmi_z[7] = -dylmi_z[6];

/*--------------------------------------------------------------------------*/
/* m=2,-2                                                                   */

  dylmr_x[8] = dy22r_phi*pre_phi_x + dy22r_theta*pre_theta_x;
  dylmi_x[8] = dy22i_phi*pre_phi_x + dy22i_theta*pre_theta_x;
  dylmr_y[8] = dy22r_phi*pre_phi_y + dy22r_theta*pre_theta_y;
  dylmi_y[8] = dy22i_phi*pre_phi_y + dy22i_theta*pre_theta_y;
  dylmr_z[8] =                     + dy22r_theta*pre_theta_z;
  dylmi_z[8] =                     + dy22i_theta*pre_theta_z;
  dylmr_x[9] =  dylmr_x[8];
  dylmi_x[9] = -dylmi_x[8];
  dylmr_y[9] =  dylmr_y[8];
  dylmi_y[9] = -dylmi_y[8];
  dylmr_z[9] =  dylmr_z[8];
  dylmi_z[9] = -dylmi_z[8];

/*==========================================================================*/
/* II.iv) grad l=3                                                          */

/*--------------------------------------------------------------------------*/
/* m=0                                                                      */

  dylmr_x[10] = dy30r_theta*pre_theta_x;
  dylmi_x[10] = 0.0;
  dylmr_y[10] = dy30r_theta*pre_theta_y;
  dylmi_y[10] = 0.0;
  dylmr_z[10] = dy30r_theta*pre_theta_z;
  dylmi_z[10] = 0.0;

/*--------------------------------------------------------------------------*/
/* m=1,-1                                                                   */

  dylmr_x[11] = dy31r_phi*pre_phi_x + dy31r_theta*pre_theta_x;
  dylmi_x[11] = dy31i_phi*pre_phi_x + dy31i_theta*pre_theta_x;
  dylmr_y[11] = dy31r_phi*pre_phi_y + dy31r_theta*pre_theta_y;
  dylmi_y[11] = dy31i_phi*pre_phi_y + dy31i_theta*pre_theta_y;
  dylmr_z[11] =         + dy31r_theta*pre_theta_z;
  dylmi_z[11] =         + dy31i_theta*pre_theta_z;
  dylmr_x[12] =  dylmr_x[11];
  dylmi_x[12] = -dylmi_x[11];
  dylmr_y[12] =  dylmr_y[11];
  dylmi_y[12] = -dylmi_y[11];
  dylmr_z[12] =  dylmr_z[11];
  dylmi_z[12] = -dylmi_z[11];

/*--------------------------------------------------------------------------*/
/* m=2,-2                                                                   */

  dylmr_x[13] = dy32r_phi*pre_phi_x + dy32r_theta*pre_theta_x;
  dylmi_x[13] = dy32i_phi*pre_phi_x + dy32i_theta*pre_theta_x;
  dylmr_y[13] = dy32r_phi*pre_phi_y + dy32r_theta*pre_theta_y;
  dylmi_y[13] = dy32i_phi*pre_phi_y + dy32i_theta*pre_theta_y;
  dylmr_z[13] =                     + dy32r_theta*pre_theta_z;
  dylmi_z[13] =                     + dy32i_theta*pre_theta_z;
  dylmr_x[14] =  dylmr_x[13];
  dylmi_x[14] = -dylmi_x[13];
  dylmr_y[14] =  dylmr_y[13];
  dylmi_y[14] = -dylmi_y[13];
  dylmr_z[14] =  dylmr_z[13];
  dylmi_z[14] = -dylmi_z[13];

/*--------------------------------------------------------------------------*/
/* m=3,-3                                                                   */

  dylmr_x[15] = dy33r_phi*pre_phi_x + dy33r_theta*pre_theta_x;
  dylmi_x[15] = dy33i_phi*pre_phi_x + dy33i_theta*pre_theta_x;
  dylmr_y[15] = dy33r_phi*pre_phi_y + dy33r_theta*pre_theta_y;
  dylmi_y[15] = dy33i_phi*pre_phi_y + dy33i_theta*pre_theta_y;
  dylmr_z[15] =         + dy33r_theta*pre_theta_z;
  dylmi_z[15] =         + dy33i_theta*pre_theta_z;
  dylmr_x[16] =  dylmr_x[15];
  dylmi_x[16] = -dylmi_x[15];
  dylmr_y[16] =  dylmr_y[15];
  dylmi_y[16] = -dylmi_y[15];
  dylmr_z[16] =  dylmr_z[15];
  dylmi_z[16] = -dylmi_z[15];

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_vps_atm_list(PSEUDO *pseudo, CELL *cell, CLATOMS_POS *clatoms_pos, 
                          CLATOMS_INFO *clatoms_info, ATOMMAPS *atommaps, 
                          EWD_SCR *ewd_scr, FOR_SCR *for_scr,
                          int cp_dual_grid_opt,int itime)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int i;
  int natm_tot          = clatoms_info->natm_tot;
  double *x             = clatoms_pos->x;
  double *y             = clatoms_pos->y;
  double *z             = clatoms_pos->z;
  double *x0            = pseudo->x0w;
  double *y0            = pseudo->y0w;
  double *z0            = pseudo->z0w;
  double tol_update     = 0.5*(pseudo->nlvps_skin);

  double dx,dy,dz,rsq,tol_now;

/*=======================================================================*/ 
  if(itime == 0){
    vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                atommaps,ewd_scr,for_scr,cp_dual_grid_opt,itime);
  }else{
    tol_now = 0.0;
   for(i=1; i<= natm_tot; i++){
     dx = x[i] - x0[i];
     dy = y[i] - y0[i];
     dz = z[i] - z0[i];
     rsq = dx*dx + dy*dy + dz*dz;
     tol_now = MAX(rsq,tol_now);
   }
    tol_now = sqrt(tol_now);
   if(tol_now > tol_update){
    vps_atm_list(pseudo,cell,clatoms_pos,clatoms_info,
                atommaps,ewd_scr,for_scr,cp_dual_grid_opt,itime);
   }/*endif*/
  }/*endif*/



/*-----------------------------------------------------------------------*/
   }/* end routine */
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vps_atm_list(PSEUDO *pseudo, CELL *cell, CLATOMS_POS *clatoms_pos, 
                  CLATOMS_INFO *clatoms_info, ATOMMAPS *atommaps, 
                  EWD_SCR *ewd_scr, FOR_SCR *for_scr,
                  int cp_dual_grid_opt,int itime)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

  int iatm,iang,ishift,loc_now,iii,i,j,iii_old;
  int natm_typ_nl,natm_typ_nl_gh,isame,icount,ioff;
  int nmall,ishift2;
  int np_loc_cp_box;
  int np_nonloc_cp_box;
  int np_nonloc_cp_box_kb;
  int np_nonloc_cp_box_gh;
  int count,count_np_nonloc_cp_box_kb,count_np_nonloc_cp_box_gh;
  int *inonloc_index;
  static int ifirst = 0;
  double dx,dy,dz;
  double sx,sy,sz;
  double asx,asy,asz;

 /* Local pointers */
 /*----------------*/

  int natm_tot          = clatoms_info->natm_tot;
  int nab_initio        = clatoms_info->nab_initio;
  int *iatm_typ         = atommaps->iatm_atm_typ;
  int *iatm_typ_nl,*iatm_typ_nl_rev;
  int *imap_atm_typ_nl,*imap_atm_typ_nl_rev;
  int *imap_atm_typ_nl_gh;

  double *x             = clatoms_pos->x;
  double *y             = clatoms_pos->y;
  double *z             = clatoms_pos->z;

  double *x0            = pseudo->x0w;
  double *y0            = pseudo->y0w;
  double *z0            = pseudo->z0w;

  int *ip_nl            = pseudo->ip_nl;
  int *ip_nl_gh         = pseudo->ip_nl_gh;

  int *np_nl            = pseudo->np_nl;
  int *np_nl_gh         = pseudo->np_nl_gh;

  int **np_nl_rad_str   = pseudo->np_nl_rad_str;
  int **np_nl_rad_end   = pseudo->np_nl_rad_end;
  int *ip_loc_cp_box    = pseudo->ip_loc_cp_box;

  int *ip_nl_rev        = pseudo->ip_nl_rev;     /*length natm*nang_max*/
  int *ip_nl_rev_gh     = pseudo->ip_nl_rev_gh;  /*length natm*nang_max*/

  /* list: index from full atom array into nonlocal cp atom array*/
  int *ivps_label       = pseudo->ivps_label;
  int *n_ang            = pseudo->n_ang;
  int *loc_opt          = pseudo->loc_opt;
  int n_ang_max         = pseudo->n_ang_max;
  int n_rad_max         = pseudo->n_rad_max;
  int *nrad;
  int *nrad_0           = pseudo->nrad_0;
  int *nrad_1           = pseudo->nrad_1;
  int *nrad_2           = pseudo->nrad_2;
  int *nrad_3           = pseudo->nrad_3;

  double *hmat_cp       = cell->hmat_cp;
  double *cp_box_center = cell->cp_box_center;
  double *hmat_big      = cell->hmat;
  double *hmati_big     = cell->hmati;

/*=======================================================================*/

  for(iatm=1;iatm<=natm_tot;iatm++){
    x0[iatm] = x[iatm];
    y0[iatm] = y[iatm];
    z0[iatm] = z[iatm];
  }/*endfor*/

/*=======================================================================*/
/* 0) Reorder the system                                                 */

   if(n_rad_max>1){
     non_loc_chng_ord(clatoms_pos,clatoms_info,atommaps,pseudo,
                      ewd_scr,for_scr,itime); 
   }/*endif*/

/*=======================================================================*/
/* I) Make a list of the local/non-local atoms on the small grid         */
 
  np_loc_cp_box = 0;
  np_nonloc_cp_box = 0;
  np_nonloc_cp_box_kb = 0;
  np_nonloc_cp_box_gh = 0;

 /* A. Count up the total number of atoms with non-local pseudo potentials */
  for(iatm=1;iatm<=natm_tot;iatm++){
   if(ivps_label[iatm_typ[iatm]]==1 || 
      ivps_label[iatm_typ[iatm]]==2 || 
      ivps_label[iatm_typ[iatm]]==3 || 
      ivps_label[iatm_typ[iatm]]==5 ){
        loc_now = loc_opt[iatm_typ[iatm]];
        if(loc_now > 0){  
         np_nonloc_cp_box++;
        }/*endif*/
   }/*endif*/
  }/*endfor*/

  inonloc_index = (int *) cmalloc(np_nonloc_cp_box*sizeof(int))-1;

  /* Make a list of all the nonlocal atoms EXCEPT gauss-hermite */
    icount = 0;
  for(iatm=1;iatm<=natm_tot;iatm++){
   if(ivps_label[iatm_typ[iatm]]==1 ||  /* kleinman-bylander */
      ivps_label[iatm_typ[iatm]]==3 ||  /* vanderbilt */
      ivps_label[iatm_typ[iatm]]==5 ){  /* goedecker */
        loc_now = loc_opt[iatm_typ[iatm]];
        if(loc_now > 0){  
         icount++;
         np_nonloc_cp_box_kb++;
         inonloc_index[icount] = iatm;         
        }/*endif*/
   }/*endif*/
  }/*endfor*/

  /* Add the GAUSS-HERMITE atoms at the  bottom of the nonlocal list */

  for(iatm=1;iatm<=natm_tot;iatm++){
   if(ivps_label[iatm_typ[iatm]]== 2  ){  
        loc_now = loc_opt[iatm_typ[iatm]];
        if(loc_now > 0){  
         icount++;
         np_nonloc_cp_box_gh++;
         inonloc_index[icount] = iatm;         
        }/*endif*/
   }/*endif*/
  }/*endfor*/

  pseudo->np_nonloc_cp_box    = np_nonloc_cp_box;
  pseudo->np_nonloc_cp_box_kb = np_nonloc_cp_box_kb;
  pseudo->np_nonloc_cp_box_gh = np_nonloc_cp_box_gh;

/*==================================================================*/
/* Assign number of KB and GH atoms in each angular momentum channel*/

 /*i) Zero the arrays */

 for(iang=1;iang<=(n_ang_max+1);iang++){np_nl[iang] = 0; np_nl_gh[iang]= 0;}

  count_np_nonloc_cp_box_kb = 0; 
  count_np_nonloc_cp_box_gh = 0; 

  for(iatm=1;iatm<=natm_tot;iatm++){
    dx  = x[iatm] - cp_box_center[1];
    dy  = y[iatm] - cp_box_center[2];
    dz  = z[iatm] - cp_box_center[3]; 

    asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
    asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
    asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];
    sx  = asx - NINT(asx);
    sy  = asy - NINT(asy);
    sz  = asz - NINT(asz);
    dx  = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7];
    dy  = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8];
    dz  = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9];

    if(((fabs(dx) <= 0.5*hmat_cp[1]) &&
        (fabs(dy) <= 0.5*hmat_cp[5]) &&
        (fabs(dz) <= 0.5*hmat_cp[9]))|| (cp_dual_grid_opt==0)){
      if(cp_dual_grid_opt>=1){
        np_loc_cp_box++;
        ip_loc_cp_box[np_loc_cp_box] = iatm;
      }/*endif*/

      /*  KLEINMAN-BYLANDER and GOEDECKER */
      if(ivps_label[iatm_typ[iatm]]==1 || /*KB*/
         ivps_label[iatm_typ[iatm]]==3 || /*VdB*/
         ivps_label[iatm_typ[iatm]]==5 ){ /*Goedecker*/

        loc_now = loc_opt[iatm_typ[iatm]];

        if(loc_now > 0 ){count_np_nonloc_cp_box_kb++;}

        for(iang=1;iang<=loc_now;iang++){
          np_nl[iang]  += 1;    
          ishift        = (iang - 1)*natm_tot + np_nl[iang];
          ip_nl[ishift] = iatm;
          ip_nl_rev[ishift] = count_np_nonloc_cp_box_kb;
        }/*endfor*/

        for(iang=loc_now+2;iang<=(n_ang[iatm_typ[iatm]]+1);iang++){
          np_nl[iang]  += 1;    
          ishift        = (iang - 1)*natm_tot + np_nl[iang];
          ip_nl[ishift] = iatm;
          ip_nl_rev[ishift] = count_np_nonloc_cp_box_kb;
        }/*endfor*/
      }/*endif: atm has a KB GOEDECKER pseudopotential*/


      /*  GAUSS-HERMITE */
      if(ivps_label[iatm_typ[iatm]]==2 ){ /*GAUSS-HERMITE*/
        loc_now = loc_opt[iatm_typ[iatm]];

        if(loc_now > 0 ){count_np_nonloc_cp_box_gh++;}

        for(iang=1;iang<=loc_now;iang++){
          np_nl_gh[iang]      += 1;  
          ishift               = (iang - 1)*natm_tot + np_nl_gh[iang];
          ip_nl_gh[ishift]     = iatm;
          ip_nl_rev_gh[ishift] = count_np_nonloc_cp_box_gh;
        }/*endfor*/

        for(iang=loc_now+2;iang<=(n_ang[iatm_typ[iatm]]+1);iang++){
          np_nl_gh[iang]      += 1;  
          ishift               = (iang - 1)*natm_tot + np_nl_gh[iang];
          ip_nl_gh[ishift]     = iatm;
          ip_nl_rev_gh[ishift] = count_np_nonloc_cp_box_gh;
        }/*endfor*/
      }/*endif: atm has a GAUSS-HERMITE pseudopotential*/


    }/*endif particle in cp_box or the cp_box is the box*/
  }/*endfor*/

  pseudo->np_loc_cp_box    = np_loc_cp_box;

/*-----------------------------------------------------------------*/
/* First time only: malloc and assign array that holds atom types  */
/*   for atoms with non-local pseudopotentials                     */
/*-----------------------------------------------------------------*/

  if(ifirst == 0){

   atommaps->iatm_atm_typ_nl     = (int *) cmalloc(np_nonloc_cp_box*sizeof(int))-1;
   atommaps->iatm_atm_typ_nl_rev = (int *) cmalloc(np_nonloc_cp_box*sizeof(int))-1;

   iatm_typ_nl     = atommaps->iatm_atm_typ_nl;
   iatm_typ_nl_rev = atommaps->iatm_atm_typ_nl_rev;

/*Nonlocal atom's  global atom type  */

   for(i=1; i<= np_nonloc_cp_box; i++){
     iatm = inonloc_index[i]; 
     iatm_typ_nl[i] = iatm_typ[iatm];
   }/*endfor*/

/*Count how many non-local atom types there are*/
    natm_typ_nl = 0;

    for(i=1; i<= np_nonloc_cp_box; i++){
      isame = 0;
     for(j=(i-1); j >= 1; j--){
      if(iatm_typ_nl[j] == iatm_typ_nl[i]){
        isame = 1;
        break;
      }/*endif*/
     }/*endfor*/
      if(isame == 0){natm_typ_nl++;}/*endif*/
    }/*endfor*/

    pseudo->natm_typ_nl = natm_typ_nl;

/* Assign atom's nonlocal atom type */

   if( np_nonloc_cp_box > 0){
    iatm_typ_nl_rev[1] = 1;
    icount = 1;
   }/*endif*/
   for(i=2; i<= np_nonloc_cp_box; i++){
    for(j=(i-1); j>= 1;j--){
       isame = 0;
      if(iatm_typ_nl[i] == iatm_typ_nl[j]){
        isame = 1;
        break;
      }/*endif*/
    }/*endfor*/
    if(isame==0){/* not the same*/
      icount++;
      iatm_typ_nl_rev[i] = icount;
    }else{/*the same*/
      iatm_typ_nl_rev[i] = iatm_typ_nl_rev[j];
    }

   }/*endfor*/


/* Create a map: nonlocal atom type j is of global atom type k */

   atommaps->imap_atm_typ_nl = (int *) cmalloc(natm_typ_nl*sizeof(int))-1;
   imap_atm_typ_nl  = atommaps->imap_atm_typ_nl;

   if( np_nonloc_cp_box > 0){
     icount = 1;
     imap_atm_typ_nl[1] = iatm_typ_nl[1];
   }

    for(i=2; i<= np_nonloc_cp_box; i++){
      isame = 0;
     for(j=(i-1); j >= 1; j--){
      if(iatm_typ_nl_rev[i] == iatm_typ_nl_rev[j]){
        isame = 1;
        break;
      }/*endif*/
     }/*endfor*/
      if(isame == 0){
	icount++;
       imap_atm_typ_nl[icount] = iatm_typ_nl[i];
      }/*endif*/
    }/*endfor*/

/*--------------------------------------------------------------*/
/* Assign GAUSS-HERMITE atom's GAUSS-HERMITE nonlocal atom type */
/*--------------------------------------------------------------*/
  if(np_nonloc_cp_box_gh > 0){
   atommaps->imap_atm_typ_nl_gh = (int *) cmalloc(np_nonloc_cp_box_gh*sizeof(int))-1;
   imap_atm_typ_nl_gh  = atommaps->imap_atm_typ_nl_gh;

             count = 1;
    natm_typ_nl_gh = 1;
    imap_atm_typ_nl_gh[count] = natm_typ_nl_gh;
    
    ioff = (np_nonloc_cp_box-np_nonloc_cp_box_gh);

   for(i=(ioff+2); i<= np_nonloc_cp_box; i++){ 
     if(ivps_label[(iatm_typ_nl[i])] == 2){
      isame = 0;
      count++;

     for(j=(i-1); j >= (ioff+1); j--){
      if((ivps_label[(iatm_typ_nl[j])] == 2) /*GAUSS-HERMITE*/
       &&(iatm_typ_nl[i] == iatm_typ_nl[j])){
        isame = j;
        break;
      }/*endif*/
     }/*endfor*/

       if(isame == 0){ 
        natm_typ_nl_gh++;
        imap_atm_typ_nl_gh[count] = natm_typ_nl_gh;
       }else{ 
        imap_atm_typ_nl_gh[count] = imap_atm_typ_nl_gh[j-ioff];
       }

     }/*end if atom1 GAUSS-HERMITE*/
    }/*endfor*/

  }/*endif there are gauss-hermite nonlocal atoms */

     ifirst = 1;

  }/*endif ifirst*/

/*=======================================================================*/
/* Find the radial channel divisions in the list                         */

  if(n_rad_max>1){

    for(iang=1;iang<=(n_ang_max+1);iang++){
      for(i=1;i<=n_rad_max;i++){np_nl_rad_end[iang][i]=0;}
      if(iang==1){nrad=nrad_0;}
      if(iang==2){nrad=nrad_1;}
      if(iang==3){nrad=nrad_2;}
      if(iang==4){nrad=nrad_3;}
      ishift = (iang - 1)*natm_tot;
      iii_old = 0;
      for(i=1;i<=np_nl[iang];i++){
        iii = nrad[iatm_typ[ip_nl[(ishift+i)]]];
        np_nl_rad_end[iang][iii]++;
#ifdef DEBUG_GJM
        if(iii<iii_old){printf("Bogus dude \n");}
        iii_old = iii;
#endif
      }/*endfor*/
      np_nl_rad_str[iang][1]=1;
      for(i=2;i<=n_rad_max;i++){
        np_nl_rad_str[iang][i] =np_nl_rad_end[iang][(i-1)]+1;
        np_nl_rad_end[iang][i]+=np_nl_rad_end[iang][(i-1)];
      }/*endfor*/
#ifdef DEBUG_GJM
      printf("L %d num %d\n",iang-1,np_nl[iang]);
      for(i=1;i<=n_rad_max;i++){
        printf("L %d N %d str %d end %d\n",iang-1,i,
                         np_nl_rad_str[iang][i],np_nl_rad_end[iang][i]);
      }/*endfor*/
      printf("scanf : ");scanf("%d",&i);
#endif
    }/*endfor*/

  }else{

    for(iang=1;iang<=(n_ang_max+1);iang++){
      np_nl_rad_str[iang][1]=1;
      np_nl_rad_end[iang][1]=np_nl[iang];
    }/*endfor*/

  }/*endif*/

/*=======================================================================*/
/* III) Restore order                                                    */

   if(n_rad_max>1){
     non_loc_restore_ord(clatoms_pos,clatoms_info,atommaps,pseudo,
                         ewd_scr,for_scr); 
   }/*endif*/

   cfree(&(inonloc_index[1]));

/*-----------------------------------------------------------------------*/
   }/* end routine */
/*=======================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void non_loc_restore_ord(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                         ATOMMAPS *atommaps, PSEUDO *pseudo, 
                         EWD_SCR *ewd_scr, FOR_SCR *for_scr) 

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*==========================================================================*/
/*         Local Variable declarations                                      */

  int i;

  int natm_tot      = clatoms_info->natm_tot;
  int *iatm_typ     = atommaps->iatm_atm_typ;
  int *i_tmp        = for_scr->iexcl;
  int *map          = pseudo->map_nl;

  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;

  double *q         = clatoms_info->q;

  double *fx        = clatoms_pos->fx;
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;

  double *x_tmp     = ewd_scr->fx;
  double *y_tmp     = ewd_scr->fy;
  double *z_tmp     = ewd_scr->fz;

/*==========================================================================*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[map[i]] = x[i];
    y_tmp[map[i]] = y[i];
    z_tmp[map[i]] = z[i];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    x[i] = x_tmp[i];
    y[i] = y_tmp[i];
    z[i] = z_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[map[i]] = q[i];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    q[i] = x_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    i_tmp[map[i]] = iatm_typ[i];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    iatm_typ[i] = i_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[map[i]] = fx[i];
    y_tmp[map[i]] = fy[i];
    z_tmp[map[i]] = fz[i];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    fx[i] = x_tmp[i];
    fy[i] = y_tmp[i];
    fz[i] = z_tmp[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void non_loc_chng_ord(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                      ATOMMAPS *atommaps, PSEUDO *pseudo, 
                      EWD_SCR *ewd_scr, FOR_SCR *for_scr,int itime) 

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*==========================================================================*/
/*         Local Variable declarations                                      */

  int i;

  int natm_tot      = clatoms_info->natm_tot;
  int *iatm_typ     = atommaps->iatm_atm_typ;
  int *i_tmp        = for_scr->iexcl;

  int *map_nl       = pseudo->map_nl;
  int *nrad_0       = pseudo->nrad_0;
  int *nrad_1       = pseudo->nrad_1;
  int *nrad_2       = pseudo->nrad_2;
  int *nrad_3       = pseudo->nrad_3;
  int *ivps_label   = pseudo->ivps_label;
  int *rank;

  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;

  double *q         = clatoms_info->q;

  double *fx        = clatoms_pos->fx;
  double *fy        = clatoms_pos->fy;
  double *fz        = clatoms_pos->fz;

  double *x_tmp     = ewd_scr->fx;
  double *y_tmp     = ewd_scr->fy;
  double *z_tmp     = ewd_scr->fz;

/*==========================================================================*/
/* Sort the atoms and make a map so that they appear in */
/* the order of fewest to most channels                 */

  if(itime==0){

    rank = (int *)cmalloc(natm_tot*sizeof(int))-1;
 
    for(i=1;i<=natm_tot;i++){
      if(ivps_label[iatm_typ[i]]==1 || 
         ivps_label[iatm_typ[i]]==2 ||
         ivps_label[iatm_typ[i]]==3 ||
         ivps_label[iatm_typ[i]]==5){
        rank[i] = 1000*nrad_0[iatm_typ[i]]
                +  100*nrad_1[iatm_typ[i]]
                +   10*nrad_2[iatm_typ[i]]
                +      nrad_3[iatm_typ[i]];
      }else{
        rank[i] = 0;
      }/*endif*/
      map_nl[i]  = i;
    }/*endfor*/

    if(natm_tot>1){sort_commence(natm_tot,rank,map_nl);}

#ifdef DEBUG_GJM
    for(i=1;i<=natm_tot;i++){
      printf("%d rank %d map_nl %d\n",i,rank[i],map_nl[i]);
    }/*endfor*/
    printf("scanf : ");scanf("%d",&i);
#endif

    cfree(&rank[1]);

  }/*endif*/

/*==========================================================================*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[i] = x[map_nl[i]];
    y_tmp[i] = y[map_nl[i]];
    z_tmp[i] = z[map_nl[i]];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    x[i] = x_tmp[i];
    y[i] = y_tmp[i];
    z[i] = z_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[i] = q[map_nl[i]];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    q[i] = x_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    i_tmp[i] = iatm_typ[map_nl[i]];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    iatm_typ[i] = i_tmp[i];
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){
    x_tmp[i] = fx[map_nl[i]];
    y_tmp[i] = fy[map_nl[i]];
    z_tmp[i] = fz[map_nl[i]];
  }/*endfor*/
  for(i=1;i<=natm_tot;i++){
    fx[i] = x_tmp[i];
    fy[i] = y_tmp[i];
    fz[i] = z_tmp[i];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/



 
