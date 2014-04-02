/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: min_DIIS_cp                                  */
/*                                                                          */
/* This subprogram minimizes the CP wavefunction using DIIS                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_integrate_cpmin_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void min_DIIS_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                 CP *cp,int iatm_count,int ip_now)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int i,ip=1,ipart,icoef,is;
    double dt,wght1,wght2;
    int iist,iipst,iii;
    int ncoef_up_tot;
    int ncoef_dn_tot;
    int maxcoef;
    int ipvt_flag;
    double ovlap_tmp;
    double gamma=1.0;


/*----------------*/
/* Local Pointers */

    int  natm_tot               = class->clatoms_info.natm_tot;
    double *class_clatoms_fx   = class->clatoms_pos[ip_now].fx;
    double *class_clatoms_fy   = class->clatoms_pos[ip_now].fy;
    double *class_clatoms_fz   = class->clatoms_pos[ip_now].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;
    double *ptens_pvten         = general_data->ptens.pvten;
    double *ptens_pvten_tot     =  general_data->ptens.pvten_tot;

    double *cre_up     = cp->cpcoeffs_pos[ip_now].cre_up;
    double *cim_up     = cp->cpcoeffs_pos[ip_now].cim_up;
    double *cre_dn     = cp->cpcoeffs_pos[ip_now].cre_dn;
    double *cim_dn     = cp->cpcoeffs_pos[ip_now].cim_dn;
    double *fcre_up    = cp->cpcoeffs_pos[ip_now].fcre_up;
    double *fcim_up    = cp->cpcoeffs_pos[ip_now].fcim_up;
    double *fcre_dn    = cp->cpcoeffs_pos[ip_now].fcre_dn;
    double *fcim_dn    = cp->cpcoeffs_pos[ip_now].fcim_dn;
    double *cmass      = cp->cpcoeffs_info.cmass;
    double *cp_hess_re_up = cp->cpcoeffs_pos[ip_now].cp_hess_re_up;
    double *cp_hess_im_up = cp->cpcoeffs_pos[ip_now].cp_hess_im_up;
    double *cp_hess_re_dn = cp->cpcoeffs_pos[ip_now].cp_hess_re_dn;
    double *cp_hess_im_dn = cp->cpcoeffs_pos[ip_now].cp_hess_im_dn;
    int ncoef          = cp->cpcoeffs_info.ncoef;
    int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
    int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
    int icoef_form_up        = cp->cpcoeffs_pos[ip_now].icoef_form_up;
    int icoef_orth_up        = cp->cpcoeffs_pos[ip_now].icoef_orth_up;
    int ifcoef_form_up        = cp->cpcoeffs_pos[ip_now].ifcoef_form_up;
    int ifcoef_orth_up        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_up;
    int icoef_form_dn        = cp->cpcoeffs_pos[ip_now].icoef_form_dn;
    int icoef_orth_dn        = cp->cpcoeffs_pos[ip_now].icoef_orth_dn;
    int ifcoef_form_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_form_dn;
    int ifcoef_orth_dn        = cp->cpcoeffs_pos[ip_now].ifcoef_orth_dn;
    int *ioff_up       = cp->cpcoeffs_info.ioff_upt;
    int *ioff_dn       = cp->cpcoeffs_info.ioff_dnt;

    int nstate_up               = cp->cpcoeffs_info.nstate_up;
    int nstate_dn               = cp->cpcoeffs_info.nstate_dn;
    int cp_lsda                 = cp->cpopts.cp_lsda;
    int cg_reset_flag           = cp->cpcoeffs_info.cg_reset_flag;
    int cp_norb                 = cp->cpopts.cp_norb;
    int ireset                  = cp->cpcoeffs_info.diis_reset_flag;
    int np_states            = cp->communicate.np_states;
    int myid_state           = cp->communicate.myid_state;
    int cp_para_opt          = cp->cpopts.cp_para_opt;
    int ioff_hyb;
    int ncoef_up,ncoef_dn,ncoef_up_max,ncoef_dn_max;
    int MVEC_MAX = general_data->minopts.cp_diis_hist_len;
    MPI_Comm comm_states = class->communicate.comm_states;

/*-----------*/
/* DIIS defs */

    int job;
    int j,k,ist,ipst,jst,jpst;
    int mvec,mp1;  
    int info;
    static int itcount;

    double ovlap;  
    double vtemp,mag;
    double sum_check;

/* DIIS local pointers */
    static int *ipvt;
    static double *pmat_up;
    static double *vect_up;
    static double *cstre_up;
    static double *cstim_up;
    static double *estre_up;
    static double *estim_up;
    static double *pmat_dn;
    static double *vect_dn;
    static double *cstre_dn;
    static double *cstim_dn;
    static double *estre_dn;
    static double *estim_dn;


/*==========================================================================*/
/* 0) Checks                                                                */

  if(cp_norb>0){
    if((icoef_orth_up+ifcoef_orth_up)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are in orthonormal form \n");
      printf("on state processor %d in min_STD_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_orth_dn+ifcoef_orth_dn)!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn CP vectors are in orthonormal form \n");
       printf("on state processor %d in min_STD_cp \n",myid_state);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/

  if(np_states>1){
    if((icoef_form_up+ifcoef_form_up)!=2){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up CP vectors are not in transposed form \n");
     printf("on state processor %d in min_STD_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if((icoef_form_dn+ifcoef_form_dn)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up CP vectors are not in transposed form \n");
      printf("on state processor %d in min_STD_cp \n",myid_state);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
    }/*endif*/
  }/*endif*/


/*==========================================================================*/
/* 0) Useful constants                                                      */

    general_data->timeinfo.int_res_tra    = 0;
    general_data->timeinfo.int_res_ter    = 0;
    dt  = general_data->timeinfo.dt;
    general_data->stat_avg.iter_shake     = 0; 
    general_data->stat_avg.iter_ratl     = 0; 

    if(np_states==1){
     ncoef_up     = cp->cpcoeffs_info.ncoef;
     ncoef_dn     = cp->cpcoeffs_info.ncoef;
     ncoef_up_max = cp->cpcoeffs_info.ncoef;
     ncoef_dn_max  = cp->cpcoeffs_info.ncoef;
    }else{
     ncoef_up     = cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
     ncoef_dn     = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
     ncoef_up_max = cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
     ncoef_dn_max = cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
    }/*endif*/

    ncoef_up_tot = ncoef_up_max*nstate_up;
    ncoef_dn_tot = ncoef_dn_max*nstate_dn; 

    if(ireset == 1) itcount=1;
    mvec = ( itcount < MVEC_MAX ? itcount:MVEC_MAX);
    mp1 = mvec + 1;

/*==========================================================================*/
/* I.a) Setup initial DIIS memory                                           */
    
    if(itcount == 1){
      if(iatm_count > 0){
        ipvt_flag=1;
        realloc_DIIS_mem(&ipvt,&pmat_up,&vect_up,&cstre_up,&cstim_up,
                         &estre_up,&estim_up,mp1,ncoef_up_tot,mvec,ipvt_flag);
      } else {
       ipvt_flag=1;
       init_alloc_DIIS_mem(&ipvt,&pmat_up,&vect_up,&cstre_up,&cstim_up,
                           &estre_up,&estim_up,mp1,ncoef_up_tot,mvec,ipvt_flag);
      } /* endif iatm_count */

      if(cp_lsda == 1){
       if(iatm_count > 0){
        ipvt_flag=0;
        realloc_DIIS_mem(&ipvt,&pmat_dn,&vect_dn,&cstre_dn,&cstim_dn,
                         &estre_dn,&estim_dn,mp1,ncoef_dn_tot,mvec,ipvt_flag);
       } else {
        ipvt_flag=0;
        init_alloc_DIIS_mem(&ipvt,&pmat_dn,&vect_dn,&cstre_dn,&cstim_dn,
                           &estre_dn,&estim_dn,mp1,ncoef_dn_tot,mvec,ipvt_flag);
       } /* endif iatm_count */
      } /*endif lsda */

    } /*endif I.a) */

 if ( itcount >= 1 && itcount <= MVEC_MAX){

/*==========================================================================*/
/* II.a) Reallocate DIIS memory                                             */


    ipvt_flag=1;
    realloc_DIIS_mem(&ipvt,&pmat_up,&vect_up,&cstre_up,&cstim_up,
                     &estre_up,&estim_up,mp1,ncoef_up_tot,mvec,ipvt_flag);

    if(cp_lsda == 1){
     ipvt_flag=0;      
     realloc_DIIS_mem(&ipvt,&pmat_dn,&vect_dn,&cstre_dn,&cstim_dn,
                      &estre_dn,&estim_dn,mp1,ncoef_dn_tot,mvec,ipvt_flag);
    } /*endif lsda */

 } else {

/*==========================================================================*/
/* II.b) Maximum length of history vector reached, so shift coef histories */

  shift_DIIS_hist(cstre_up,cstim_up,estre_up,estim_up,mvec,ncoef_up_tot);


  if(cp_lsda == 1){
   shift_DIIS_hist(cstre_dn,cstim_dn,estre_dn,estim_dn,mvec,ncoef_dn_tot);
  } /*endif lsda */

 }/*endif (II) */

/*==========================================================================*/
/* III) Get forces                                                          */

   for(i=1;i<= class->clatoms_info.natm_tot;i++){
     class_clatoms_fx[i]  = 0.0;
     class_clatoms_fy[i]  = 0.0;
     class_clatoms_fz[i]  = 0.0;
   }/*endfor*/

   for(i=1;i<=ncoef_up_tot;i++){
     fcre_up[i] = 0.0;
     fcim_up[i] = 0.0;
   }/*endfor*/
   if( (cp_lsda == 1) && (nstate_dn != 0) ){
     for(i=1;i<=ncoef_dn_tot; i++){
      fcre_dn[i] = 0.0;
      fcim_dn[i] = 0.0;
     }/*endfor*/
   }/*endif*/


   for(i=1;i<=9;i++){
     ptens_pvten[i]     = 0.0;
     ptens_pvten_tot[i] = 0.0;
   }/*endfor*/
   general_data->stat_avg.cp_ehart = 0.0;
   general_data->stat_avg.cp_exc   = 0.0;
   general_data->stat_avg.cp_muxc  = 0.0;
   general_data->stat_avg.cp_eext  = 0.0;
   general_data->stat_avg.cp_enl   = 0.0;
   general_data->stat_avg.cp_eke   = 0.0;
   general_data->stat_avg.vrecip   = 0.0;

   for(i=1;i<=ncoef;i++){
     cp_hess_re_up[i] = 0.0;
     cp_hess_im_up[i] = 0.0;
   }
   if(cp->cpopts.cp_lsda == 1 && nstate_dn != 0){
    for(i=1;i<=ncoef;i++){
      cp_hess_re_dn[i] = 0.0;
      cp_hess_im_dn[i] = 0.0;
    }/* endfor */
   }/* endfor */
   cp_ks_energy_ctrl(cp,ip_now,&(general_data->ewald),&(class->ewd_scr),
                       &(general_data->cell),
                       &(class->clatoms_info),
                       &(class->clatoms_pos[ip_now]),
                       &(class->atommaps),&(general_data->stat_avg),
                       &(general_data->ptens),
                       &(general_data->simopts),
                       &(class->for_scr));

   get_diag_cp_hess(cp,ip_now,&(general_data->cell),gamma);


/*==========================================================================*/
/* IV) Get approximate gradients, evolve coeffs and save new histories      */ 

    ist=(mvec-1)*ncoef_up_tot;
    ioff_hyb = (cp_para_opt == 0 ? myid_state*ncoef_up_max : 0);
    for(is=1;is<=nstate_up;is++) {
      for(i=1;i<=ncoef_up;i++) {
        icoef = i+ioff_up[is];
        iist  = ist+icoef;
        estre_up[iist] = fcre_up[icoef]/cp_hess_re_up[ioff_hyb+i];
        estim_up[iist] = fcim_up[icoef]/cp_hess_im_up[ioff_hyb+i];
      }/* endfor */
    }/* endfor */

    for(is=1;is <= nstate_up; is++){
      for(i=1;i<=ncoef_up;i++){
        icoef = i+ioff_up[is];
        iist  = ist+icoef;
        cre_up[icoef] += estre_up[iist];
        cim_up[icoef] += estim_up[iist];
      }
    }

    for(is=1;is<=nstate_up;is++){
      for(i=1;i<=ncoef_up;i++){
        icoef = i+ioff_up[is];
        iist = ist+icoef;
        cstre_up[iist] = cre_up[icoef];
        cstim_up[iist] = cim_up[icoef];
       }/*endfor i */
    }/*endfor is */

    if(cp_lsda == 1){
      ioff_hyb = (cp_para_opt == 0 ? myid_state*ncoef_dn_max : 0);
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+ioff_dn[is];
          iist  = ist+icoef;
          estre_dn[iist] = fcre_dn[icoef]/cp_hess_re_dn[ioff_hyb+i];
          estim_dn[iist] = fcim_dn[icoef]/cp_hess_im_dn[ioff_hyb+i];
        }/* endfor */
      }/* endfor */

      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_dn;i++) {
          icoef = i+ioff_dn[is];
          iist  = ist+icoef;
          cre_dn[icoef] += estre_dn[iist];
          cim_dn[icoef] += estim_dn[iist];
        }/* endfor */
      }/* endfor */

      for(is=1;is<=nstate_dn;is++){
        for(i=1;i<=ncoef_dn;i++){
          icoef = i+ioff_dn[is];
          iist = ist+icoef;
          cstre_dn[iist] = cre_dn[icoef];
          cstim_dn[iist] = cim_dn[icoef];
         }/*endfor i */
      }/*endfor is */
    }/* endif lsda */

    if(itcount >= 2){

/*==========================================================================*/
/* VI.a) Construct DIIS linear equations                                    */ 
 
    setup_DIIS_eqs(pmat_up,vect_up,estre_up,estim_up,&(cp->communicate),
                   nstate_up,ncoef_up,ncoef_up_tot,ioff_up,mp1,mvec);

  if(cp_lsda == 1){
    setup_DIIS_eqs(pmat_dn,vect_dn,estre_dn,estim_dn,&(cp->communicate),
                   nstate_dn,ncoef_dn,ncoef_dn_tot,ioff_dn,mp1,mvec);
  }/* endif */
   
/*==========================================================================*/
/* VI.b) Solve DIIS equations to get DIIS vectors                           */

#ifdef IBM_ESSL
      dgef(&(pmat_up[1]),&mp1,&mp1,&(ipvt[1]));
      job = 0;
      dges(&(pmat_up[1]),&mp1,&mp1,&(ipvt[1]),&(vect_up[1]),&job);
 
      if(cp_lsda == 1){
        dgef(&(pmat_dn[1]),&mp1,&mp1,&(ipvt[1]));
        job = 0;
        dges(&(pmat_dn[1]),&mp1,&mp1,&(ipvt[1]),&(vect_dn[1]),&job);
      }/* endif lsda */ 
#else


      DGEFA(&(pmat_up[1]),&mp1,&mp1,&(ipvt[1]),&info);
      job=0;
      DGESL(&(pmat_up[1]),&mp1,&mp1,&(ipvt[1]),&(vect_up[1]),&job);
 
      if(cp_lsda == 1){
        DGEFA(&(pmat_dn[1]),&mp1,&mp1,&(ipvt[1]),&info);
        job = 0;
        DGESL(&(pmat_dn[1]),&mp1,&mp1,&(ipvt[1]),&(vect_dn[1]),&job);
      }/* endif lsda */ 


#endif

/*==========================================================================*/
/* VI.c) Compute new coefs from DIIS vectors                                */

    for(is=1;is<=ncoef_up_tot;is++){
      cre_up[is] = 0.0;
      cim_up[is] = 0.0;
    }/*endfor*/
    for(i=1;i<=mvec;i++){
      ist=(i-1)*ncoef_up_tot;
      for(is=1;is<=ncoef_up_tot;is++) {
        iist = ist+is;
        cre_up[is] += vect_up[i]*cstre_up[iist];
        cim_up[is] += vect_up[i]*cstim_up[iist];
      }/*endfor is */
    }/*endfor i */

    if(cp_lsda == 1){
      for(is=1;is<=ncoef_dn_tot;is++){
        cre_dn[is] = 0.0;
        cim_dn[is] = 0.0;
      }/*endfor*/
      for(i=1;i<=mvec;i++){
        ist=(i-1)*ncoef_dn_tot;
        for(is=1;is<=ncoef_dn_tot;is++) {
          iist = ist+is;
          cre_dn[is] += vect_dn[i]*cstre_dn[iist];
          cim_dn[is] += vect_dn[i]*cstim_dn[iist];
        }/*endfor is */
      }/*endfor i */
    }/* endif lsda */ 

   }/*endif DIIS lin equations */

/*==========================================================================*/
/* V) Orthogonalize wave functions                                          */

   orthog_control_cp(cp,ip_now);
   ++itcount;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void realloc_DIIS_mem(int **ipvt,double **pmat,double **vect,
                 double **cstre,double **cstim,double **estre,double **estim,
                 int mp1,int ncoef_tot,int mvec,int ipvt_flag)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* I) Reallocate DIIS arrays */

 if(ipvt_flag == 1){
  (*ipvt) = (int *) crealloc(&((*ipvt)[1]),mp1*sizeof(int))-1;
 }
 (*pmat) = (double *) crealloc(&((*pmat)[1]),mp1*mp1*sizeof(double))-1;
 (*vect) = (double *) crealloc(&((*vect)[1]),mp1*sizeof(double))-1;
 (*cstre) = (double *) crealloc(&((*cstre)[1]),ncoef_tot*mvec*sizeof(double))-1;
 (*cstim) = (double *) crealloc(&((*cstim)[1]),ncoef_tot*mvec*sizeof(double))-1;
 (*estre) = (double *) crealloc(&((*estre)[1]),ncoef_tot*mvec*sizeof(double))-1;
 (*estim) = (double *) crealloc(&((*estim)[1]),ncoef_tot*mvec*sizeof(double))-1;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_alloc_DIIS_mem(int **ipvt,double **pmat,double **vect,
                 double **cstre,double **cstim,double **estre,double **estim,
                 int mp1,int ncoef_tot,int mvec,int ipvt_flag)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/* I) Allocate DIIS arrays */

  (*ipvt) =  (int *) cmalloc(mp1*sizeof(int))-1;
  (*pmat) =  (double *) cmalloc(mp1*mp1*sizeof(double))-1;
  (*vect) =  (double *) cmalloc(mp1*sizeof(double))-1;
  (*cstre) = (double *) cmalloc(ncoef_tot*mvec*sizeof(double))-1;
  (*cstim) = (double *) cmalloc(ncoef_tot*mvec*sizeof(double))-1;
  (*estre) = (double *) cmalloc(ncoef_tot*mvec*sizeof(double))-1;
  (*estim) = (double *) cmalloc(ncoef_tot*mvec*sizeof(double))-1;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shift_DIIS_hist(double *cstre,double *cstim,double *estre,double *estim,
                     int mvec,int ncoef_tot)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

 int ist,ipst,iist,iipst,i,is;

/*========================================================================*/
/* I) Revise history                                                      */

  for(i=1;i<=mvec-1;i++){
    ist=(i-1)*ncoef_tot;
    ipst=i*ncoef_tot;
    for(is=1;is<=ncoef_tot;is++) {
      iist = ist+is;
      iipst = ipst+is;
      cstre[iist] = cstre[iipst];
      cstim[iist] = cstim[iipst];
      estre[iist] = estre[iipst];
      estim[iist] = estim[iipst];
    }/*endfor is */
  }/*endfor i */


/*========================================================================*/
   }/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void setup_DIIS_eqs(double *pmat,double *vect,
                    double *estre,double *estim,
                    COMMUNICATE *communicate,int nstate,int ncoef,
                    int ncoef_tot,int *ioff,int mp1,int mvec)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*  Local variable declarations                                           */

#include "../typ_defs/typ_mask.h"
  int k,i,j,l,m,n,ist,jst,is;
  int np_states        = communicate->np_states;
  int myid_state       = communicate->myid_state;
  MPI_Comm comm_states = communicate->comm_states;
  double ovlap,ovlap_tmp;
  double wght1,wght2;

/*========================================================================*/
/* I) Compute overlaps and fill DIIS matrix                               */


  k=0;
  wght1 = 2.0;
  wght2 = 2.0;
  if(myid_state+1==np_states){wght1=1.0;wght2=0.0;}
  for(i=1;i<=mp1;i++){
    for(j=1;j<=mp1;j++){
      k++;
      if(i == mp1 || j == mp1)  {
        pmat[k] = -1.0;
      } else {
        ist=(i-1)*ncoef_tot;
        jst=(j-1)*ncoef_tot;           
        ovlap = 0.0;
        for(is=1;is<=nstate;is++){
#define ALT_OVERLAP
#ifndef ALT_OVERLAP
            ovlap += diag_ovlap(ncoef,wght1,wght2,
                                &(estre[ist]),&(estim[ist]),
                                &(estre[jst]),&(estim[jst]),  
                                ioff[is]);
#endif
#ifdef ALT_OVERLAP
	    for(l=1;l<=ncoef-1;l++){
	      m = ist + ioff[is] + l;
              n = jst + ioff[is] + l;
	      ovlap += 2.0*(estre[m]*estre[n] + estim[m]*estim[n]);
            }
            m = ist + ioff[is] + ncoef;
            n = jst + ioff[is] + ncoef;
            ovlap += estre[m]*estre[n];
#endif
        }/*endfor*/
        if(np_states>1){
         ovlap_tmp = ovlap;
         Allreduce(&(ovlap_tmp),&(ovlap),1,MPI_DOUBLE,MPI_SUM,0,
                   comm_states);
        }/*endif*/
       pmat[k] = ovlap;
      } /* endif */
    } /* endfor j */
  } /* endfor i */
  pmat[(mp1*mp1)] = 0.0;
  for(i=1;i<=mvec;i++){
    vect[i] = 0.0;
  }/*endfor*/
  vect[mp1] = -1.0;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

