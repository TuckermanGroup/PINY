/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_eext_nonlocal_gh.c                   */
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

#define JUERG_FACTOR 0.72

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_nonloc_gh(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                      CELL *cell, PTENS *ptens,CPEWALD *cpewald,
                      CPSCR *cpscr, PSEUDO *pseudo, EWD_SCR *ewd_scr,  
                      CPOPTS *cpopts, ATOMMAPS *atommaps, 
                      COMMUNICATE *communicate,FOR_SCR *for_scr,
                      double rgh)
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
  int ipart,iii,i,iatm,ioff;
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

  int n_ang_max_gh   = pseudo->n_ang_max_gh;
  int n_rad_max      = pseudo->n_rad_max;

  int *ip_nl_gh      = pseudo->ip_nl_gh;
  int *np_nl_gh      = pseudo->np_nl_gh;
  double *gzvps      = pseudo->gzvps;
  double *gzvps0     = pseudo->gzvps0;

  int np_nonloc_cp_box_kb  = pseudo->np_nonloc_cp_box_kb;
  int np_nonloc_cp_box_gh  = pseudo->np_nonloc_cp_box_gh;
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

  for(i=1;i<=(n_ang_max_gh+1);i++){
    if(np_nl_gh[i]>0){nl_max=i-1; }
  }/*endfor*/

  nl_chan_max = (nl_max + 1)*(nl_max + 1);

/*======================================================================*/
/* III) Determine the maximum number of atoms in any                    */
/*       open angular momentum channel                                  */

  np_nlmax = 1;
  for(i = 1;i<=(nl_max+1);i++){
    np_nlmax = MAX(np_nlmax,np_nl_gh[i]);
  }/*endfor*/

/*======================================================================*/
/* IV) Find cos and sin of sc components of the particles               */
/*    ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                       */

  for(ipart=1;ipart<= np_nonloc_cp_box_gh;ipart++){
    iatm = ip_nl_gh[ipart];
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

/*----------------------------------------------------------------------*/
/* ii) If break point number one calculate the helpful vectors          */
 
 
   if(ibreak1_sm[icount]==1 ){
     for(ipart=1;ipart<=np_nonloc_cp_box_gh;ipart++){
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

   control_nlmat_gh(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
                    cpscr,cpopts,pseudo,ewd_scr,atommaps,
                    np_nlmax,nl_max,g,xk,yk,zk,
                    icount,&ylm_cons,rgh);

/*----------------------------------------------------------------------*/
/* iv) If break point two, increment the helpful vectors                */

   if(ibreak2_sm[icount]==1 ){
     for(ipart=1;ipart<=np_nonloc_cp_box_gh;ipart++){
     temp = helr[ipart];
     helr[ipart] = helr[ipart]*cossc[ipart] - heli[ipart]*sinsc[ipart];
     heli[ipart] = heli[ipart]*cossc[ipart] + temp*sinsc[ipart];
    }/*endfor*/
   }/*endif*/

 }/*endfor:icount loop over k vectors */


/*======================================================================*/
/* VI) g=0 term                                                         */


  ak2_sm[ncoef] = 0.0;

  if(np_nl_gh[1]>0){

    ylmr[1] = 1.0/sqrt(fpi);
    for(irad=1;irad<=nrad_max_l[1];irad++){

      get_nlmat0_gh(ncoef,ncoef,nstate_up,np_nlmax,
                    np_nl_rad_str[1][irad],np_nl_gh[1],nl_chan_max,irad,
                    creal_up,vnlreal_up,ylmr[1]);

      if(cp_lsda==1){
       get_nlmat0_gh(ncoef,ncoef,nstate_dn,np_nlmax,
                  np_nl_rad_str[1][irad],np_nl_gh[1],nl_chan_max,irad,
                  creal_dn,vnlreal_dn,ylmr[1]);
      }/*endif lsda*/

    }/*endfor*/

  }/*endif: l=0 nonlocal*/


/*======================================================================*/
   }/*end routine*/
/*======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_nlmat_gh(CLATOMS_INFO *clatoms_info,
                   CPCOEFFS_INFO *cpcoeffs_info,
                   CPCOEFFS_POS *cpcoeffs_pos,
                   CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                   EWD_SCR *ewd_scr,ATOMMAPS *atommaps,
                   int np_nlmax,int nl_max,
                   double g,double xk,double yk,double zk,
                   int ismcount,YLM_CONS *ylm_cons,double rgh)

/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/

/*=======================================================================*/
/*         Local Variable declarations                                   */

  int ind_lm,ipart,lp1,irad,jrad,ipart_nl;
  double sgn_l;
  double jl;
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

  int *ip_nl_rev_gh        = pseudo->ip_nl_rev_gh;
  int *np_nl_gh            = pseudo->np_nl_gh;

  int n_rad_max            = pseudo->n_rad_max;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;
 
  int np_nonloc_cp_box_gh  = pseudo->np_nonloc_cp_box_gh;
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
    jl = get_jl(g,rgh,l); /*Calculate jl */

    lp1 = l+1;
    if(np_nl_gh[lp1]>0){

     for(irad=1;irad<=nrad_max_l[lp1];irad++){
/*---------------------------------------------------------------------*/

      i_shift  = l*npart;

      if(cp_ptens==0) {
        for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl_gh[lp1];ipart++){
          ktemp             =  ipart+i_shift;
          ltemp             =  ip_nl_rev_gh[ktemp];
          helr_now[ipart]   =  helr[ltemp];
          heli_now[ipart]   =  heli[ltemp];
        }/*endfor*/
      }else{
        printf("@@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@@@\n");
        printf("PRESSURE TENSOR FOR GAUSS-HERMITE IS NOT IMPLEMENTED \n");
        printf("@@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@@@\n");
        exit(1);
      }/*endif*/
    
/*---------------------------------------------------------------------*/
/*ii) Loop over m components of the channel and get the matrix elements*/ 
/*     vtemp,dvtemp, */

      sgn_l = 1.0;
      if((l%2)==1)sgn_l = -1.0;

      for(m=1;m<=(2*l+1);m++){
        ind_lm = m + l*l;
          if(atm_hess_calc == 3){
             get_nlmat_gh_hess(ncoef,ismcount,nstate_up,ind_lm,irad,
                np_nlmax,np_nl_gh[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
                helr_now,heli_now,
                creal_up,cimag_up,
                vnlreal_up,vnlimag_up,
                dvnlreal_x_up,dvnlreal_y_up,dvnlreal_z_up,
                dvnlimag_x_up,dvnlimag_y_up,dvnlimag_z_up,
                dvnlreal_gxgx_up,dvnlreal_gxgy_up,
                dvnlreal_gxgz_up,dvnlreal_gygy_up,dvnlreal_gygz_up,
                dvnlreal_gzgz_up,
                dvnlimag_gxgx_up,dvnlimag_gxgy_up,dvnlimag_gxgz_up,
                dvnlimag_gygy_up,dvnlimag_gygz_up,dvnlimag_gzgz_up,
                xk,yk,zk,
                ylmr[ind_lm],ylmi[ind_lm],scr1,scr2,sgn_l,rgh,jl);
	  }/* endif */
          if(atm_hess_calc != 3){
             get_nlmat_gh(ncoef,ismcount,nstate_up,ind_lm,irad,
                np_nlmax,np_nl_gh[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
                helr_now,heli_now,
                creal_up,cimag_up,
                vnlreal_up,vnlimag_up,
                dvnlreal_x_up,dvnlreal_y_up,dvnlreal_z_up,
                dvnlimag_x_up,dvnlimag_y_up,dvnlimag_z_up,
                xk,yk,zk,
                ylmr[ind_lm],ylmi[ind_lm],scr1,scr2,sgn_l,rgh,jl);
          }/* endif */

        if(cp_lsda==1) {
          if(atm_hess_calc == 3){
             get_nlmat_gh_hess(ncoef,ismcount,nstate_dn,ind_lm,irad,
               np_nlmax,np_nl_gh[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
               helr_now,heli_now,
               creal_dn,cimag_dn,
               vnlreal_dn,vnlimag_dn,
               dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
               dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
               dvnlreal_gxgx_dn,dvnlreal_gxgy_dn,dvnlreal_gxgz_dn,
               dvnlreal_gygy_dn,dvnlreal_gygz_dn,dvnlreal_gzgz_dn,
               dvnlimag_gxgx_dn,dvnlimag_gxgy_dn,dvnlimag_gxgz_dn,
               dvnlimag_gygy_dn,dvnlimag_gygz_dn,dvnlimag_gzgz_dn,
               xk,yk,zk,
               ylmr[ind_lm],ylmi[ind_lm],scr1,scr2,sgn_l,rgh,jl);
          }/* endif */

          if(atm_hess_calc != 3){
             get_nlmat_gh(ncoef,ismcount,nstate_dn,ind_lm,irad,
               np_nlmax,np_nl_gh[lp1],np_nl_rad_str[lp1][irad],nl_chan_max,
               helr_now,heli_now,
               creal_dn,cimag_dn,
               vnlreal_dn,vnlimag_dn,
               dvnlreal_x_dn,dvnlreal_y_dn,dvnlreal_z_dn,
               dvnlimag_x_dn,dvnlimag_y_dn,dvnlimag_z_dn,
               xk,yk,zk,
               ylmr[ind_lm],ylmi[ind_lm],scr1,scr2,sgn_l,rgh,jl);
          }/* endif */

        }/*endif:lsda*/

      }/*endfor:m quantum number loop*/
     }/*endfor: irad quantum number loop */
    }/*endif: this channel has atoms in it*/
  }/*endfor:l quantum number loop*/   


/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*==========================================================================*/
  double get_jl(double g,double rgh,int l)
/*==========================================================================*/
   {/*start routine*/
/*==========================================================================*/
  
   double arg,arg2;
   double jl;

/*==========================================================================*/
/* II) Get some constants                                                    */

   arg  = g*rgh*JUERG_FACTOR;

   switch(l){
     case 0:
        jl = sin(arg)/arg;
      break;
     case 1:
        jl = sin(arg)/(arg*arg) - cos(arg)/arg; 
      break;
     case 2:
       arg2 = arg*arg;
       jl = (3.0/(arg2*arg) - 1.0/arg)*sin(arg) - (3.0/(arg2))*cos(arg);
/*     jl = (3.0/(arg*arg*arg) - 1.0/arg)*sin(arg) - (3.0/(arg*arg))*cos(arg); */
      break;
     case 3:
       arg2 = arg*arg;
       jl = (15.0/(arg2*arg2) - 6.0/arg2)*sin(arg) 
          - (15.0/(arg2*arg)  - 1.0/arg )*cos(arg);
      break;


   }/*end switch*/

    return jl;
/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_nlmat_gh(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
               int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
               double *helr,double *heli,
               double *creal,double *cimag,
               double *vnlreal,double *vnlimag,
               double *dvnlreal_x,double *dvnlreal_y,double *dvnlreal_z,
               double *dvnlimag_x,double *dvnlimag_y,double *dvnlimag_z,
               double xk,double yk,double zk,double ylmr,double ylmi,
               double *cbyhelr,double *cbyheli,double sgn_l,
               double rgh,double jl)

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
/* Here helr,heli DO NOT have vtemp mutliplied in them                      */
/*   See control routine                                                    */
/* storage : atm,state,l,m,n */

  isgn_l = (int) (0.5*(sgn_l+1.0)+1.0);

  tylmr = 2.0*ylmr*jl;
  tylmi = 2.0*ylmi*jl;

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind   = creal[ind_loc_c];
    cim_ind   = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    switch(isgn_l){
       case 1:    /*ODD*/
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
       case 2:   /*EVEN*/
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

void get_nlmat_gh_hess(int ncoef,int ismcount,int nstate,int ind_lm,int irad,
               int np_nlmax,int np_nl,int np_nl_rad_str,int nl_chan_max,
               double *helr,double *heli,
               double *creal,double *cimag,
               double *vnlreal,double *vnlimag,
               double *dvnlreal_x,double *dvnlreal_y,double *dvnlreal_z,
               double *dvnlimag_x,double *dvnlimag_y,double *dvnlimag_z,
	       double *dvnlreal_gxgx,double *dvnlreal_gxgy,
               double *dvnlreal_gxgz,double *dvnlreal_gygy,double *dvnlreal_gygz,
               double *dvnlreal_gzgz,double *dvnlimag_gxgx,double *dvnlimag_gxgy,
               double *dvnlimag_gxgz,double *dvnlimag_gygy,double *dvnlimag_gygz,
               double *dvnlimag_gzgz,
               double xk,double yk,double zk,double ylmr,double ylmi,
               double *cbyhelr,double *cbyheli,double sgn_l,
               double rgh,double jl)

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
  double dtmp3;

/*==========================================================================*/
/* I) Loop over states and particles:                                       */
/* Here helr,heli DO NOT have vtemp mutliplied in them                      */
/*   See control routine                                                    */
/* storage : atm,state,l,m,n */

  isgn_l = (int) (0.5*(sgn_l+1.0)+1.0);

  tylmr = 2.0*ylmr*jl;
  tylmi = 2.0*ylmi*jl;

  for(is=1;is<=nstate;is++){
    ind_loc_c = ismcount + ncoef*(is-1);
    cre_ind   = creal[ind_loc_c];
    cim_ind   = cimag[ind_loc_c];
    ioff_v    = (is-1)*np_nlmax + 
                (ind_lm-1)*nstate*np_nlmax
               +(irad-1)*nl_chan_max*nstate*np_nlmax;
    switch(isgn_l){
       case 1:    /*ODD*/
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
       case 2:   /*EVEN*/
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
    }/*end switch: sgn of m in Ylm*/

  }/* endfor : loop over states */

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void get_nlmat0_gh(int ncoef,int ismcount,int nstate,
                   int np_nlmax,int np_nl_rad_str,
                   int np_nl,int nl_chan_max,int irad,
                   double *creal,double *vnlreal,double ylmr)

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
      term  = cre_ind;

      vnlreal[(ipart+ioff_v)] += term*ylmr;
    }/*endfor*/

  }/* endfor : loop over the number of states */

/*------------------------------------------------------------------------*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_pot_pv_fatm_gh(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  CELL *cell,CPCOEFFS_INFO *cpcoeffs_info,CPSCR *cpscr,
                  EWD_SCR *ewd_scr,CPOPTS *cpopts,
                  PSEUDO *pseudo,ATOMMAPS *atommaps,
                  double *cp_enl_ret, int np_nlmax,double *pvten,
                  double *wgh,int igh)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,l,m,i_shift,ipart,is,iii,ioff,lp1;
  int iatm,index_gh;
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

  double *hmat_cp          = cell->hmat_cp;

  int cp_ptens             = cpopts->cp_ptens_calc;
  int cp_lsda              = cpopts->cp_lsda;
  int atm_hess_calc        = clatoms_info->hess_calc;
  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;
  int *imap_atm_typ_nl_gh  = atommaps->imap_atm_typ_nl_gh;
  int *ip_nl_rev_gh        = pseudo->ip_nl_rev_gh;
  double *vpsnorm          = pseudo->vpsnorm; 
  int    ngh               = pseudo->ngh; 
  int n_ang_max_gh         = pseudo->n_ang_max_gh;
  int n_rad_max            = pseudo->n_rad_max;
  int *loc_opt             = pseudo->loc_opt;

  int *ip_nl_gh            = pseudo->ip_nl_gh;

  int *np_nl_gh            = pseudo->np_nl_gh;
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
  cp_enl   = 0.0000;
  nl_max   = -1;

  for(i=1;i<=(n_ang_max_gh+1);i++){
   if(np_nl_gh[i]>0){nl_max=i-1;}
  }/*endfor*/

  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Loop over the open channels, the states and get the nl potent,   */
/*     pvten and  particle forces                                       */


  for(l=0;l<=nl_max;l++){
    lp1 = l+1;
    if(np_nl_gh[lp1]>0){
     for(irad=1;irad<=nrad_max_l[lp1];irad++){
     for(jrad=irad;jrad<=nrad_max_l[lp1];jrad++){
       
/*-----------------------------------------------------------------------*/
/* i) Get the Gauss-Hermite weight scaled by the volume     */

      i_shift = l*npart;
      for(ipart=np_nl_rad_str[lp1][jrad];ipart<=np_nl_gh[lp1];ipart++){
        index_gh = ip_nl_rev_gh[i_shift+ipart];
        ioff = (imap_atm_typ_nl_gh[index_gh]-1)*(n_ang_max_gh)*ngh
             +  l*ngh; 
         vnorm_now[ipart] = wgh[ioff + igh]*rvol_cp; 
      }/*endfor*/


/*-----------------------------------------------------------------------*/
/* ii) Sum the contributions over the 2l+1 directions and the states    */


       sumnl_pot_pv_fatm_hess(npart,nstate_up,np_nlmax,nl_chan_max,
                              np_nl_gh[lp1],l,np_nl_rad_str[lp1][jrad],
                              irad,jrad,
                              ip_nl_gh,vnorm_now,vnlreal_up,vnlimag_up,
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
                              hess_xx,hess_xy,hess_xz,
                              hess_yy,hess_yz,hess_zz,
                              atm_hess_calc,cp_ptens,pvten,&cp_enl);

       if(cp_lsda==1){
         sumnl_pot_pv_fatm_hess(npart,nstate_dn,np_nlmax,nl_chan_max,
                                np_nl_gh[lp1],l,np_nl_rad_str[lp1][jrad],
                                irad,jrad,
                                ip_nl_gh,vnorm_now,vnlreal_dn,vnlimag_dn,
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
                                hess_xx,hess_xy,hess_xz,
                                hess_yy,hess_yz,hess_zz,
                                atm_hess_calc,cp_ptens,pvten,&cp_enl);
       }/*endif*/
     }}/*endfor: radial channels */
    }/*endif: l channel open */
  }/*endfor: l channels     */

/*======================================================================*/
/* III) Assign the non-local energy  and add it to the pvten            */

  *cp_enl_ret = cp_enl;
  if(cp_ptens==1){
     printf("@@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@@@\n");
     printf("PRESSURE TENSOR FOR GAUSS-HERMITE IS NOT IMPLEMENTED \n");
     printf("@@@@@@@@@@@@@@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@@@@@@@@@\n");
     exit(1);

    /*
    pvten[1] += cp_enl;
    pvten[5] += cp_enl;
    pvten[9] += cp_enl;
    */
  }/*endif*/

/*======================================================================*/
  }/*end routine*/
/*======================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void getnl_fcoef_gh(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                    CPSCR *cpscr,EWD_SCR *ewd_scr,CPOPTS *cpopts,
                    PSEUDO *pseudo,CPEWALD *cpewald,ATOMMAPS *atommaps,
                    CELL *cell, int np_nlmax,double *pvten,FOR_SCR *for_scr,
                    double rgh,double *wgh,int igh)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int ipart,nl_max,iii,i,l,istr,irad,jrad,iatm,ioff;
  int icount,ismcount,nl_chan_max;
  double wgh0;
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

  int n_ang_max_gh   = pseudo->n_ang_max_gh;

  int *ip_nl_gh      = pseudo->ip_nl_gh;
  int *ip_nl_rev_gh  = pseudo->ip_nl_rev_gh;
  int *imap_atm_typ_nl_gh = atommaps->imap_atm_typ_nl_gh;
  int  np_nonloc_cp_box_gh = pseudo->np_nonloc_cp_box_gh;

  int *np_nl_gh      = pseudo->np_nl_gh;
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
  for(i=1;i<=(n_ang_max_gh+1);i++){
   if(np_nl_gh[i]>0){nl_max=i-1;}
  }/*endfor*/

  nl_chan_max = (nl_max+1)*(nl_max+1);

/*======================================================================*/
/* II) Find cos and sin of sc components of charged particles           */
/*  ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                         */

  for(ipart=1;ipart<= np_nonloc_cp_box_gh;ipart++){
    iatm = ip_nl_gh[ipart];
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

    if(ibreak1_sm[icount]==1 ){
      for(ipart=1;ipart<=np_nonloc_cp_box_gh;ipart++){
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
    control_nlfcoef_gh(clatoms_info,cpcoeffs_info,cpcoeffs_pos,
                       cpscr,cpopts,pseudo,ewd_scr,atommaps,
                       np_nlmax,nl_max,
                       g,xk,yk,zk,vol_cp,ismcount,&ylm_cons,
                       rgh,wgh,igh);

/*----------------------------------------------------------------------*/
/* iv) If break point two, increment the helpful vectors                */

    if(ibreak2_sm[icount]==1 ){
      for(ipart=1;ipart<=np_nonloc_cp_box_gh;ipart++){
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

  if(np_nl_gh[1]>0){
   
    l = 0;
    ylmr[1] = 1.0/sqrt(fpi);
    ismcount = nktot_sm + 1;

    wgh0 = 0.0;
    for(i=1; i<= np_nonloc_cp_box_gh; i++){
      ioff = (imap_atm_typ_nl_gh[i]-1)*(pseudo->n_ang_max_gh)*pseudo->ngh
           +  l*pseudo->ngh; 
      wgh0 += wgh[ioff+igh];
    }  


    for(irad=1;irad<=nrad_max_l[1];irad++){
      for(jrad=1;jrad<=nrad_max_l[1];jrad++){

        istr = MAX(np_nl_rad_str[1][irad],np_nl_rad_str[1][jrad]);

        get_nlfor0_gh(ncoef,ismcount,nstate_up,np_nlmax,np_nl_gh[1],
                      istr,nl_chan_max,jrad,
                      fcreal_up,vnlreal_up,ylmr[1],vol_cp,wgh0);
 /*note: wgh 1D array of length n_ang*ngh - elements 1-ngh are l=0 channel*/

#ifdef HESS_GH
        if(cp_hess_calc == 1){
          get_nlhess0(ncoef,ismcount,nstate_up,np_nlmax,np_nl_gh[1],
                      istr,nl_chan_max,jrad,vtemp_now,vnorm_now,cp_hess_re_up,
                      vnlreal_up,vnlimag_up,ylmr[1],vol_cp);
        }/*endif*/
#endif

        if(cp_lsda==1){
          get_nlfor0_gh(ncoef,ismcount,nstate_dn,np_nlmax,np_nl_gh[1],
                        istr,nl_chan_max,jrad,
                        fcreal_dn,vnlreal_dn,ylmr[1],vol_cp,wgh0);

#ifdef HESS_GH
          if(cp_hess_calc == 1){
           get_nlhess0(ncoef,ismcount,nstate_dn,np_nlmax,np_nl_gh[1],
                       istr,nl_chan_max,jrad,vtemp_now,vnorm_now,cp_hess_re_dn,
                       vnlreal_dn,vnlimag_dn,ylmr[1],vol_cp);
          }/*endif*/
#endif
        }/*endif:lsda*/

      }/*endfor radial channel*/
    }/*endfor radial channel*/
  }/*endif:l=0 channel open*/

/*======================================================================*/
   }/*end routine*/
/*======================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_nlfcoef_gh(CLATOMS_INFO *clatoms_info,
                        CPCOEFFS_INFO *cpcoeffs_info,
                        CPCOEFFS_POS *cpcoeffs_pos,
                        CPSCR *cpscr,CPOPTS *cpopts,PSEUDO *pseudo,
                        EWD_SCR *ewd_scr,ATOMMAPS *atommaps,
                        int np_nlmax,int nl_max,
                        double g,double xk,double yk,double zk,double vol,
                        int ismcount,YLM_CONS *ylm_cons,
                        double rgh,double *wgh,int igh)

/*==========================================================================*/
    {/*Begin Routine */
/*==========================================================================*/
/*     Local Variables */

  int l,itype,ind_lm,m,ipart,i_shift,lp1,ltemp,ltemp2,ktemp;
  int nl_chan_max,istr,irad,jrad,ioff,index_gh;
  double sgn_l,rvol4,tmp1,jl;
  double ylmr[21], ylmi[21];
  double dylmr_gx[21], dylmi_gx[21];
  double dylmr_gy[21], dylmi_gy[21];
  double dylmr_gz[21], dylmi_gz[21];


/* Local pointers */     
  int npart                = clatoms_info->natm_tot;
  int natm_typ             = atommaps->natm_typ;

  int nstate_up            = cpcoeffs_info->nstate_up_proc;
  int nstate_dn            = cpcoeffs_info->nstate_dn_proc;
  int ncoef                = cpcoeffs_info->ncoef;
  double *fcreal_up        = cpcoeffs_pos->fcre_up;  
  double *fcimag_up        = cpcoeffs_pos->fcim_up;
  double *fcreal_dn        = cpcoeffs_pos->fcre_dn;
  double *fcimag_dn        = cpcoeffs_pos->fcim_dn;
  int cp_lsda              = cpopts->cp_lsda;
  int cp_hess_calc         = cpopts->cp_hess_calc;
  int *imap_atm_typ_nl_gh  = atommaps->imap_atm_typ_nl_gh;
  
  int  n_ang_max_gh        = pseudo->n_ang_max_gh;
  int  ngh                 = pseudo->ngh;
  int *np_nl_gh            = pseudo->np_nl_gh;

  int *ip_nl_rev_gh          = pseudo->ip_nl_rev_gh;
  int  np_nonloc_cp_box_kb = pseudo->np_nonloc_cp_box_kb;
  int nsplin_g             = pseudo->nsplin_g;
  double dg_spl            = pseudo->dg_spl;
  double gmin_spl          = pseudo->gmin_spl;
  double *vps0             = pseudo->vps0;
  double *vps1             = pseudo->vps1;
  double *vps2             = pseudo->vps2;
  double *vps3             = pseudo->vps3;
  double *vpsnorm          = pseudo->vpsnorm;
  int n_rad_max            = pseudo->n_rad_max;
  int *nrad_max_l          = pseudo->nrad_max_l;
  int **np_nl_rad_str      = pseudo->np_nl_rad_str;


  double *vnlreal_up       = cpscr->cpscr_nonloc.vnlre_up;
  double *vnlimag_up       = cpscr->cpscr_nonloc.vnlim_up;
  double *vnlreal_dn       = cpscr->cpscr_nonloc.vnlre_dn;
  double *vnlimag_dn       = cpscr->cpscr_nonloc.vnlim_dn;
  double *cp_hess_re_up    = cpcoeffs_pos->cp_hess_re_up;
  double *cp_hess_im_up    = cpcoeffs_pos->cp_hess_im_up;
  double *cp_hess_re_dn    = cpcoeffs_pos->cp_hess_re_dn;
  double *cp_hess_im_dn    = cpcoeffs_pos->cp_hess_im_dn;

  double *helr             = ewd_scr->helr;   
  double *heli             = ewd_scr->heli;
  double *helr_now         = ewd_scr->helr_now;
  double *heli_now         = ewd_scr->heli_now;
  double *vtemp            = ewd_scr->q;
  double *vtemp_now        = ewd_scr->vtemp_now;


  double *vscr             = ewd_scr->fz; 
  double *vfactr           = ewd_scr->fx2; 
  double *vfacti           = ewd_scr->fy2; 
  double *vnorm            = ewd_scr->fx;

/*==========================================================================*/
/* I) Get the ylm and other constants                                       */

  get_ylm(xk,yk,zk,g,ylmr,ylmi,dylmr_gx,dylmi_gx,dylmr_gy,dylmi_gy,
          dylmr_gz,dylmi_gz,ylm_cons);

  rvol4 = 4.0/vol;

/*==========================================================================*/
/* II) Loop over the channels                                               */

    nl_chan_max = (nl_max+1)*(nl_max+1);
    for(l=0;l<=nl_max;l++){

      jl = get_jl(g,rgh,l); /*Calculate jl */

      lp1 = l+1;

      sgn_l = 1.0;
      if(l%2==1)sgn_l = -1.0;

      if(np_nl_gh[lp1]>0){
       for(irad=1;irad<=nrad_max_l[lp1];irad++){

/*-----------------------------------------------------------------------*/
/* i) Create vnorm = wgh*4*rvol*jl(gR)                                   */
        i_shift = l*npart;
        for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl_gh[lp1];ipart++){
          ktemp            = ipart+i_shift;
          ltemp            = ip_nl_rev_gh[ktemp];
          helr_now[ipart]  = helr[ltemp];
          heli_now[ipart]  = heli[ltemp];
          index_gh = ip_nl_rev_gh[i_shift+ipart];
          ioff = (imap_atm_typ_nl_gh[index_gh]-1)*(n_ang_max_gh)*ngh
               +  l*ngh; 

          vnorm[ipart] = wgh[ioff+igh]*jl*rvol4;
        }/*endfor*/
     
/*-----------------------------------------------------------------------*/
/* i) Loop over the channels and compute the forces                      */

       for(m=1;m<=(2*l+1);m++){
         ind_lm = m + l*l;

         for(ipart=np_nl_rad_str[lp1][irad];ipart<=np_nl_gh[lp1];ipart++){
           vfactr[ipart]   =  (helr_now[ipart]*ylmr[ind_lm]
                              -heli_now[ipart]*ylmi[ind_lm]);
           vfacti[ipart]   = -(heli_now[ipart]*ylmr[ind_lm]
                              +helr_now[ipart]*ylmi[ind_lm]);
         }/*endfor*/

         for(jrad=1;jrad<=nrad_max_l[lp1];jrad++){

           istr = MAX(np_nl_rad_str[lp1][irad],np_nl_rad_str[lp1][jrad]);

           get_nlfor_gh(ncoef,ismcount,nstate_up,ind_lm,
                     np_nlmax,np_nl_gh[lp1],
                     istr,nl_chan_max,jrad,
                     vfactr,vfacti,fcreal_up,fcimag_up,
                     vnlreal_up,vnlimag_up,vnorm,sgn_l);
                     
#ifdef HESS_GH
/* CHECK GET_NLHESS for GAUSS-HERMITE*/
            if(cp_hess_calc == 1){   /* MODIFY THIS?!*/
printf("am I going to get_nlhess?  exit \n"); exit(1);
              get_nlhess(ncoef,ismcount,nstate_up,ind_lm,
                         np_nlmax,np_nl_gh[lp1],
                         istr,nl_chan_max,jrad,
                         vtemp_now,vnorm_now,cp_hess_re_up,
                         vnlreal_up,vnlimag_up,
                         ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            }/*endif*/
#endif
            if(cp_lsda==1){
              get_nlfor_gh(ncoef,ismcount,nstate_dn,ind_lm,
                        np_nlmax,np_nl_gh[lp1],
                        istr,nl_chan_max,jrad,
                        vfactr,vfacti,fcreal_dn,fcimag_dn,
                        vnlreal_dn,vnlimag_dn,vnorm,sgn_l);
                        
#ifdef HESS_GH
	      /* MODIFY THIS?!*/
            if(cp_hess_calc == 1){
              get_nlhess(ncoef,ismcount,nstate_up,ind_lm,
                           np_nlmax,np_nl_gh[lp1],
                           istr,nl_chan_max,jrad,
                           vtemp_now,vnorm_now,cp_hess_re_dn,
                           vnlreal_dn,vnlimag_dn,
                           ylmr[ind_lm],ylmi[ind_lm],sgn_l,vol);
            }/*endif*/
#endif
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

void get_nlfor_gh(int ncoef,int ismcount,int nstate,int ind_lm,int np_nlmax,
               int np_nl,int np_nl_rad_str,int nl_chan_max,int jrad,
               double *vfactr,double *vfacti,double *fcreal,double *fcimag,
               double *vnlreal,double *vnlimag,double *vnorm,double sgn_l)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*        Local variable declarations            */

  int iii,is,ipart,ind_loc_c,ind_loc_v,ioff_v;
  double tmp1,tmp2;

/*==========================================================================*/
/* Loop over the number of states and the number of particles               */
/*     Here helr,heli are cos(gR) and sin(gR)                               */
/*       vnorm = 4.0*jl(gR)*wgh                                             */

    for(is=1;is<=nstate;is++){
      ioff_v = (is-1)*np_nlmax + 
               +(ind_lm-1)*nstate*np_nlmax
               +(jrad-1)*nl_chan_max*nstate*np_nlmax;
      ind_loc_c  = ismcount + ncoef*(is-1);

      for(ipart=np_nl_rad_str;ipart<=np_nl;ipart++){
        ind_loc_v = ipart + ioff_v;
        tmp1      = -(vfactr[ipart]*vnlreal[ind_loc_v]
                     -vfacti[ipart]*vnlimag[ind_loc_v])*vnorm[ipart];
        tmp2      = -(vfacti[ipart]*vnlreal[ind_loc_v]
                     +vfactr[ipart]*vnlimag[ind_loc_v])*vnorm[ipart];
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

void get_nlfor0_gh(int ncoef,int ismcount,int nstate,
                   int np_nlmax,int np_nl,
                   int istr,int nl_chan_max,int jrad,
                   double *fcreal,double *vnlreal,
                   double ylmr,double vol,double wgh)

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

      vfactr  =  wgh*ylmr*rvol;
      fcreal[ind_loc_c] -= 2.0*vfactr*vnlreal[ind_loc_v];
    }/* endfor : loop over the number of particles*/

  }/* endfor : loop over the number of states */

/*------------------------------------------------------------------------*/
  }/*end routine*/
/*========================================================================*/

