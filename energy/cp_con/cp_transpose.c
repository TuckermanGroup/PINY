/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: cp_con_utils.c                               */
/*                                                                          */
/* Contains utilities for wave function constraints and norb                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_OFF


/*==========================================================================*/
/* High level Fwd Transpose controller:nstate_procXncoef->nstateXncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_coef_transpose_fwd(CP *cp, int iopt)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  int ip,iii;
  int np_states         = cp->communicate.np_states;
  double *cpscr_cre_up  = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up  = cp->cpscr.cpscr_wave.cim_up;
  int pi_beads_proc     = cp->cpcoeffs_info.pi_beads_proc;
  int cp_lsda           = cp->cpopts.cp_lsda;
  int nstate_dn         = cp->cpcoeffs_info.nstate_dn;
  int iup=1;int idn=0;

/*======================================================================*/

 if(np_states>1){

   for(ip=1;ip<=pi_beads_proc;ip++){
     cp_transpose_fwd((cp->cpcoeffs_pos[ip].cre_up),
                      (cp->cpcoeffs_pos[ip].cim_up),
                      &(cp->cpcoeffs_pos[ip].icoef_form_up),
                      cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                      &(cp->cp_comm_state_pkg_up));
     if(cp_lsda==1 && nstate_dn != 0){
       cp_transpose_fwd((cp->cpcoeffs_pos[ip].cre_dn),
                        (cp->cpcoeffs_pos[ip].cim_dn),
                        &(cp->cpcoeffs_pos[ip].icoef_form_dn),
                        cp->cpscr.cpscr_wave.cre_dn,
                        cp->cpscr.cpscr_wave.cim_dn,
                        &(cp->cp_comm_state_pkg_dn));
     }/*endif*/
   }/*endfor*/
   if(iopt>=2){
    for(ip=1;ip<=pi_beads_proc;ip++){
     cp_transpose_fwd((cp->cpcoeffs_pos[ip].vcre_up),
                      (cp->cpcoeffs_pos[ip].vcim_up),
                      &(cp->cpcoeffs_pos[ip].ivcoef_form_up),
                      cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                      &(cp->cp_comm_state_pkg_up));
     if(cp_lsda==1 && nstate_dn != 0){
       cp_transpose_fwd((cp->cpcoeffs_pos[ip].vcre_dn),
                        (cp->cpcoeffs_pos[ip].vcim_dn),
                        &(cp->cpcoeffs_pos[ip].ivcoef_form_dn),
                        cp->cpscr.cpscr_wave.cre_dn,
                        cp->cpscr.cpscr_wave.cim_dn,
                        &(cp->cp_comm_state_pkg_dn));
      }/*endif*/
     }/*endfor*/
    }/* endif iopt */
    if(iopt==3 && cp->cptherm_info.massiv_flag == 1){
      cp_therm_transpose_fwd(cp->cptherm_pos,&(cp->cptherm_info),
                             cp->cpscr.cpscr_wave.cre_up,
                             cp->cpscr.cpscr_wave.cim_up,pi_beads_proc,
                             &(cp->cp_comm_state_pkg_up),
                             cp->cpcoeffs_info.ioff_upt,iup);
      if(cp->cpopts.cp_lsda == 1 && nstate_dn != 0){
       cp_therm_transpose_fwd(cp->cptherm_pos,&(cp->cptherm_info),
                              cp->cpscr.cpscr_wave.cre_dn,
                              cp->cpscr.cpscr_wave.cim_dn,pi_beads_proc,
                              &(cp->cp_comm_state_pkg_dn),
                              cp->cpcoeffs_info.ioff_dnt,idn);
      }/* endif lsda */
      cp_therm_tran_unpack(cp->cptherm_pos,&(cp->cptherm_info),
                           cp->cp_comm_state_pkg_up.ioff_therm_tran[0],
                           pi_beads_proc);
   }/*endif iopt*/
 }/*endif*/

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/* High level Bck Transpose controller:nstate_procXncoef->nstateXncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_coef_transpose_bck(CP *cp,int iopt)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  int ip,iii;
  int np_states         = cp->communicate.np_states;
  double *cpscr_cre_up  = cp->cpscr.cpscr_wave.cre_up;
  double *cpscr_cim_up  = cp->cpscr.cpscr_wave.cim_up;
  int pi_beads_proc     = cp->cpcoeffs_info.pi_beads_proc;
  int cp_lsda           = cp->cpopts.cp_lsda;
  int nstate_dn         = cp->cpcoeffs_info.nstate_dn;
  int iup=1;int idn=0;

/*======================================================================*/

 if(np_states>1){

   for(ip=1;ip<=pi_beads_proc;ip++){
     cp_transpose_bck((cp->cpcoeffs_pos[ip].cre_up),
                      (cp->cpcoeffs_pos[ip].cim_up),
                      &(cp->cpcoeffs_pos[ip].icoef_form_up),
                      cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                      &(cp->cp_comm_state_pkg_up));
     if(cp_lsda==1){
       cp_transpose_bck((cp->cpcoeffs_pos[ip].cre_dn),
                        (cp->cpcoeffs_pos[ip].cim_dn),
                        &(cp->cpcoeffs_pos[ip].icoef_form_dn),
                        cp->cpscr.cpscr_wave.cre_dn,
                        cp->cpscr.cpscr_wave.cim_dn,
                        &(cp->cp_comm_state_pkg_dn));
     }/*endif*/
   }/*endfor*/
   if(iopt>=2){
    for(ip=1;ip<=pi_beads_proc;ip++){
     cp_transpose_bck((cp->cpcoeffs_pos[ip].vcre_up),
                      (cp->cpcoeffs_pos[ip].vcim_up),
                      &(cp->cpcoeffs_pos[ip].ivcoef_form_up),
                      cp->cpscr.cpscr_wave.cre_up,cp->cpscr.cpscr_wave.cim_up,
                      &(cp->cp_comm_state_pkg_up));
     if(cp_lsda==1){
       cp_transpose_bck((cp->cpcoeffs_pos[ip].vcre_dn),
                        (cp->cpcoeffs_pos[ip].vcim_dn),
                        &(cp->cpcoeffs_pos[ip].ivcoef_form_dn),
                        cp->cpscr.cpscr_wave.cre_dn,
                        cp->cpscr.cpscr_wave.cim_dn,
                        &(cp->cp_comm_state_pkg_dn));
     }/*endif*/
    }/*endfor*/
   }/*endif*/
   if(iopt==3 && cp->cptherm_info.massiv_flag == 1){
      cp_therm_transpose_bck(cp->cptherm_pos,&(cp->cptherm_info),
                             cp->cpscr.cpscr_wave.cre_up,
                             cp->cpscr.cpscr_wave.cim_up,pi_beads_proc,
                             &(cp->cp_comm_state_pkg_up),
                             cp->cpcoeffs_info.ioff_upt,iup);
      if(cp->cpopts.cp_lsda == 1 && nstate_dn != 0){
       cp_therm_transpose_bck(cp->cptherm_pos,&(cp->cptherm_info),
                              cp->cpscr.cpscr_wave.cre_dn,
                              cp->cpscr.cpscr_wave.cim_dn,pi_beads_proc,
                              &(cp->cp_comm_state_pkg_dn),
                              cp->cpcoeffs_info.ioff_dnt,idn);
      }/* endif lsda */
        cp_therm_tran_pack(cp->cptherm_pos,&(cp->cptherm_info),
                           cp->cp_comm_state_pkg_up.ioff_therm_tran[0],
                           pi_beads_proc);
   }/*endif iopt*/
 
 }/*endif*/

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/* Fwd Transpose controller: nstate_proc x ncoef ->  nstate x ncoef_proc    */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_fwd(double *cre,double *cim,int *icoef_form,
                      double *c1_temp,double *c2_temp,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 int iii;
 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int icoef_form_tmp;

/*========================================================================*/
/* I) Check the form                                                      */
 
  if(num_proc==1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the vectors on 1 processor\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/*endif*/
  if((*icoef_form==1)){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in normal form\n");
      printf("on state processor %d in cp_transpose_fwd \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Transpose the vector and flip the flag */

 icoef_form_tmp = 0;
 cp_transpose_fwd_prim(cre,&icoef_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                       nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                       nstate_ncoef_proc,num_proc,myid,comm);
 icoef_form_tmp = 0;
 cp_transpose_fwd_prim(cim,&icoef_form_tmp,c2_temp,nstate,nstate_max,ncoef,
                       nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                       nstate_ncoef_proc,num_proc,myid,comm);
 (*icoef_form) = 1;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/* Bck Transpose controller: nstate_proc x ncoef ->  nstate x ncoef_proc    */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_bck(double *cre,double *cim,int *icoef_form,
                      double *c1_temp,double *c2_temp,
                      CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int icoef_form_tmp;

/*========================================================================*/
/* I) Check the form                                                      */
 
  if(num_proc==1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the vectors on 1 processor\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/*endif*/
  if((*icoef_form==0)){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in normal form\n");
      printf("on state processor %d in cp_transpose_bck \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Transpose the vector and flip the flag                             */

 icoef_form_tmp = 1;
 cp_transpose_bck_prim(cre,&icoef_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                       nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                       nstate_ncoef_proc,num_proc,myid,comm);
 icoef_form_tmp = 1;
 cp_transpose_bck_prim(cim,&icoef_form_tmp,c2_temp,nstate,nstate_max,ncoef,
                       nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                       nstate_ncoef_proc,num_proc,myid,comm);
 (*icoef_form) = 0;

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/* Fwd Therm Transpose controller: nstate_proc x ncoef  ->  nstate x ncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_therm_transpose_fwd(CPTHERM_POS *cptherm_pos,CPTHERM_INFO *cptherm_info,
                            double *c1_temp,double *c2_temp,int pi_beads,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg,
                            int *ioff,int iup)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 double kinet,kinet_tmp,trash;
 int iii,ip,ich,igo,k,icf,is,itest,koff;

 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 int *ioff_therm_norm      = cp_comm_state_pkg->ioff_therm_norm;
 int *ioff_therm_tran      = cp_comm_state_pkg->ioff_therm_tran;

 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int ncoef_tot;
 int itherm_form_tmp;
 int len_c_nhc             = cptherm_info->len_c_nhc;

/*========================================================================*/
/* I) Check the form                                                      */
 
  if(cptherm_info->massiv_flag != 1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the thermostats \n");
     printf("if you are not using the massive coef \n");
     printf("thermostatting option\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }
  if(num_proc==1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the vectors on 1 processor\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Transpose the vector and flip the flag */

  ncoef_tot = ncoef*nstate_proc;
  kinet = 0.0;
  koff= 0;
  for(ip=1;ip<=pi_beads;ip++){
    for(ich=1;ich<=len_c_nhc;ich++){

     itest = (iup == 1 ? cptherm_pos[ip].itherm_form_up :
                         cptherm_pos[ip].itherm_form_dn); 
     if(itest != 0){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("The thermostats must be in normal form\n");
         printf("on state processor %d in cp_therm_transpose_fwd \n",myid);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);exit(1);
     }/*endif*/
   
     for(igo=1;igo<=2;igo++){
      koff = ioff_therm_norm[igo]-ioff_therm_norm[0];
      for(k=1;k<=ncoef_tot;k++){
        c2_temp[k] = cptherm_pos[ip].vc_nhc[ich][(k+koff)];
      }
       itherm_form_tmp = 0;
       cp_transpose_fwd_prim(c2_temp,&itherm_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                            nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                            nstate_ncoef_proc,num_proc,myid,comm);

      koff = ioff_therm_tran[igo];
      for(is=1;is <= nstate; is++){
       for(k=1;k<=nstate_ncoef_proc;k++){
         icf = k + ioff[is];
         cptherm_pos[ip].vc_nhc[ich][(k+koff)] = c2_temp[icf];
       }/* endfor  k*/
       koff += nstate_ncoef_proc;
      }/* endfor is*/
     }/* endfor igo */
    }/* endfor ich */
    if(iup==1){cptherm_pos[ip].itherm_form_up = 1;}
    if(iup==0){cptherm_pos[ip].itherm_form_dn = 1;}
   }/* endfor ip */
   

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/* Bck Therm Transpose controller: nstate_proc x ncoef  ->  nstate x ncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_therm_transpose_bck(CPTHERM_POS *cptherm_pos,CPTHERM_INFO *cptherm_info,
                            double *c1_temp,double *c2_temp,int pi_beads,
                            CP_COMM_STATE_PKG *cp_comm_state_pkg,int *ioff,int iup)

/*======================================================================*/
/*                Begin Routine */
     {/*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */
#include "../typ_defs/typ_mask.h"

 int iii,ip,ich,igo,k,icf,is,itest,koff;
 int nstate                = cp_comm_state_pkg->nstate;
 int nstate_max            = cp_comm_state_pkg->nstate_max;
 int ncoef                 = cp_comm_state_pkg->ncoef;
 int nstate_proc           = cp_comm_state_pkg->nstate_proc;
 int nstate_proc_max       = cp_comm_state_pkg->nstate_proc_max;
 int nstate_ncoef_proc_max = cp_comm_state_pkg->nstate_ncoef_proc_max;
 int nstate_ncoef_proc     = cp_comm_state_pkg->nstate_ncoef_proc;
 int num_proc              = cp_comm_state_pkg->num_proc;
 int myid                  = cp_comm_state_pkg->myid;
 int *ioff_therm_norm      = cp_comm_state_pkg->ioff_therm_norm;
 int *ioff_therm_tran      = cp_comm_state_pkg->ioff_therm_tran;
 MPI_Comm comm  = cp_comm_state_pkg->comm;
 int ncoef_tot;
 int itherm_form_tmp;
 int len_c_nhc             = cptherm_info->len_c_nhc;

/*========================================================================*/
/* I) Check the form                                                      */
 
  if(cptherm_info->massiv_flag != 1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the thermostats \n");
     printf("if you are not using the massive coef \n");
     printf("thermostatting option\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }
  if(num_proc==1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There is no need to transpose the vectors on 1 processor\n");
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);exit(1);
  }/*endif*/

/*========================================================================*/
/* II) Transpose the vector and flip the flag */

  ncoef_tot = ncoef*nstate_proc;
  for(ip=1;ip<=pi_beads;ip++){
    for(ich=1;ich<=len_c_nhc;ich++){

     itest = (iup == 1 ? cptherm_pos[ip].itherm_form_up :
                         cptherm_pos[ip].itherm_form_dn); 

     if(itest != 1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("The thermostats must be in transposed form\n");
         printf("on state processor %d in cp_therm_transpose_bck \n",myid);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);exit(1);
     }/*endif*/

      
     for(igo=1;igo<=2;igo++){
      koff = ioff_therm_tran[igo]-ioff_therm_tran[0];
      for(is=1;is <= nstate; is++){
       for(k=1;k<=nstate_ncoef_proc;k++){
         icf = k + ioff[is];
         c2_temp[icf] = cptherm_pos[ip].vc_nhc[ich][(k+koff)];
       }/* endfor  k*/
       koff += nstate_ncoef_proc;
       if(nstate_ncoef_proc < nstate_ncoef_proc_max) c2_temp[(k+ioff[is])] = 0;
      }/* endfor is*/
      itherm_form_tmp = 1;
      cp_transpose_bck_prim(c2_temp,&itherm_form_tmp,c1_temp,nstate,nstate_max,ncoef,
                            nstate_proc,nstate_proc_max,nstate_ncoef_proc_max,
                            nstate_ncoef_proc,num_proc,myid,comm);

      koff = ioff_therm_norm[igo];
      for(k=1;k<=ncoef_tot;k++){
        cptherm_pos[ip].vc_nhc[ich][(k+koff)] = c2_temp[k];
      }/* endfor k */
     }/* endfor igo */
    }/* endfor ich */
    if(iup==1){cptherm_pos[ip].itherm_form_up = 0;}
    if(iup==0){cptherm_pos[ip].itherm_form_dn = 0;}
   }/* endfor ip */
   

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/


/*==========================================================================*/
/* Perform the forward transpose: nstate_proc x ncoef-> nstate x ncoef_proc */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_transpose_fwd_prim(double *c,int *icoef_form,double *c_temp,
                          int nstate,int nstate_max,int ncoef,
                          int nstate_proc,int nstate_proc_max,
                          int nstate_ncoef_proc_max,int nstate_ncoef_proc,
                          int num_proc,int myid,MPI_Comm comm)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
  
#include "../typ_defs/typ_mask.h"

  double *c_temp_pt,*c_pt;
  int nfull,nblock;
  int irem,ig,is,ioff;
  int i,ioff_temp;
  int ioff_c;
  int joff;
  int iproc,itemp,iii;
  int proc_rem;
  int nstate_ncoef_proc_now;
  int sendcounts,recvcounts;

/*========================================================================*/
/*              Incoming variable declarations                            */

/* nstate            = Total number of states                             */
/* nstate_proc       = number of states on this processor                 */
/* nstate_proc_max   = maximum number of states on any processor          */
/* ncoef             = Total number of coefficients in a state            */
/* nstate_ncoef_proc = Number of coefficients in a state on this processor*/
/*                     in the transposed data                             */
/* nstate_ncoef_proc_max = Maximum number of coefficients in a state on   */
/*                          any processesor in the transposed data        */
/* nstate_max        = nstate_proc_max*num_proc                           */
/* c                 = nstate_proc x ncoef array of coefficients          */
/* c_temp            = transposed data stored as nstate x ncoef_proc_max  */
/* ct_temp           = scratch space to help make transposed data         */
/* nscr_size         = size of scratch nstate_ncoef_proc_max*nstate_max   */

/*========================================================================*/
/* 0) Check the forms */

    if((*icoef_form) !=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in normal form\n");
      printf("on state processor %d in cp_transpose_fwd_prim \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
    *icoef_form = 1;

/*========================================================================*/
/* I) Rearrange the coefficient data in c                                 */

  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c_temp[i]=0.0;}
  proc_rem = (ncoef % num_proc);
  for(is=1;is<=nstate_proc;is++){
    ioff = 0;
    ioff_c = (is-1)*ncoef;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (is-1)*nstate_ncoef_proc_max
                + (iproc-1)*(nstate_proc_max*nstate_ncoef_proc_max);
      nstate_ncoef_proc_now = nstate_ncoef_proc_max;
      if((iproc > proc_rem) && (proc_rem >0)) nstate_ncoef_proc_now--;
      for(ig=1;ig<=nstate_ncoef_proc_now;ig++){
        itemp = ig+ioff_temp;
        i     = ig+ioff;
        c_temp[itemp] = c[(i+ioff_c)];
      }/*endfor*/   
      ioff += nstate_ncoef_proc_now;
    }/*endfor*/   
  }/*endfor*/   
  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=c_temp[i];}

#ifdef DEBUG
  if(myid==0){
    printf("This is c in trans_fwd %d %d %d %d %d %d\n",
           nstate_proc,ncoef,nstate_ncoef_proc_max,nstate_proc_max,
           proc_rem,num_proc);
  }/*endif*/
  scanf("%d",&iii);
  Dbx_Barrier(comm);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
      printf("fwd %d %d %g\n",iproc,ig,c[ig]);
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Dbx_Barrier(comm);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Transpose the coef data                                            */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  c_pt = c+1;
  c_temp_pt  = c_temp+1;
  Alltoall(c_pt,sendcounts,MPI_DOUBLE,c_temp_pt,recvcounts,MPI_DOUBLE,comm);

#ifdef DEBUG
  if(myid==0){
    printf("This is c_temp in trans_fwd\n");
  }/*endif*/
  Dbx_Barrier(comm);
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
         printf("fwd2 %d %d %g\n",ig,is,c_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Dbx_Barrier(comm);
  }/*endfor*/
#endif
  

/*========================================================================*/
/* V) Internal rearrangement of coeff data                                */

  irem = (nstate % num_proc);

  if(irem != 0){

/*------------------------------------------------------------------------*/
/* A) copy full blocks                                                    */
    nfull = nstate_proc_max*nstate_ncoef_proc_max*irem;
    for(ig=1;ig<=nfull;ig++){c[ig] = c_temp[ig];}
/*------------------------------------------------------------------------*/
/* B) copy partial blocks                                                 */
    nblock = (nstate_proc_max-1)*nstate_ncoef_proc_max;
    ioff = nfull; joff=nfull;
    for(iproc=irem+1;iproc<=num_proc;iproc++){
     for(ig=1;ig<=nblock;ig++){c[(ig+ioff)] = c_temp[(ig+joff)];}
     ioff += nblock;joff += (nblock+nstate_ncoef_proc_max);
    }/*endfor*/

  }else{

/*------------------------------------------------------------------------*/
/* A) copy all blocks                                                     */
    nfull = nstate_ncoef_proc_max*nstate_max;
    for(ig=1;ig<=nfull;ig++){c[ig] = c_temp[ig];}

  }/*endif: remainder */
    
#ifdef DEBUG
  if(myid==0){
    printf("This is c trans_fwd\n");
  }/*endif*/
  Dbx_Barrier(comm);
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
         printf("fwd3 %d %d %g\n",ig,is,c[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Dbx_Barrier(comm);
  }/*endfor*/
    if(myid==0){printf("Completed tranform fwd\n");}
    scanf("%d",&ig);
    Dbx_Barrier(comm);
#endif
  

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/* Perform the back transpose: nstate x ncoef_proc -> nstate_proc x ncoef   */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void cp_transpose_bck_prim(double *c,int *icoef_form,double *c_temp,
                           int nstate,int nstate_max,int ncoef,
                           int nstate_proc,int nstate_proc_max,
                           int nstate_ncoef_proc_max,int nstate_ncoef_proc,
                           int num_proc,int myid,MPI_Comm comm)


/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
  
#include "../typ_defs/typ_mask.h"

  double *c_temp_pt,*c_pt,sum;
  int nfull,nblock,isum,nread,istop;
  int irem,ig,is,ioff;
  int i,ioff_temp;
  int ioff_c;
  int joff;
  int iproc,itemp,iii;
  int proc_rem;
  int nstate_ncoef_proc_now;
  int sendcounts,recvcounts;

/*========================================================================*/
/*              Incoming variable declarations                            */

/* nstate            = Total number of states                             */
/* nstate_proc       = number of states on this processor                 */
/* nstate_proc_max   = maximum number of states on any processor          */
/* ncoef             = Total number of coefficients in a state            */
/* nstate_ncoef_proc = Number of coefficients in a state on this processor*/
/*                     in the transposed data                             */
/* nstate_ncoef_proc_max = Maximum number of coefficients in a state on   */
/*                          any processesor in the transposed data        */
/* nstate_max        = nstate_proc_max*num_proc                           */
/* c                 = nstate_proc x ncoef array of coefficients          */
/* c_temp            = transposed data stored as nstate x ncoef_proc_max  */
/* ct_temp           = scratch space to help make transposed data         */
/* nscr_size         = size of scratch nstate_ncoef_proc_max*nstate_max   */

/*========================================================================*/
/* 0) Check the forms */

    if((*icoef_form) !=1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The coefficients must be in transposed form\n");
      printf("on state processor %d in cp_transpose_bck_prim \n",myid);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);exit(1);
    }/*endif*/
    *icoef_form = 0;

/*========================================================================*/
/* I) Internal rearrangement of coeff data                                */

  irem = (nstate % num_proc);

  if(irem != 0){

/*------------------------------------------------------------------------*/
/* A) copy full blocks                                                    */
    nfull = nstate_proc_max*nstate_ncoef_proc_max*irem;
    for(ig=1;ig<=nfull;ig++){c_temp[ig] = c[ig];}
/*------------------------------------------------------------------------*/
/* B) copy partial blocks                                                 */
    nblock = (nstate_proc_max-1)*nstate_ncoef_proc_max;
    ioff = nfull; joff=nfull;
    for(iproc=irem+1;iproc<=num_proc;iproc++){
     for(ig=1;ig<=nblock;ig++){c_temp[(ig+joff)] = c[(ig+ioff)];}
     ioff += nblock;joff += (nblock+nstate_ncoef_proc_max);
    }/*endfor*/
/*------------------------------------------------------------------------*/
/* C) copy back into c                                                    */
    nfull = nstate_ncoef_proc_max*nstate_max;
    for(ig=1;ig<=nfull;ig++){c[ig] = c_temp[ig];}

  }/*endif: remainder */

#ifdef DEBUG
  if(myid==0){
    printf("This is c in trans_bck %d %d %d %d %d %d\n",
           nstate_proc,ncoef,nstate_ncoef_proc_max,nstate_proc_max,
           proc_rem,num_proc);
  }/*endif*/
  scanf("%d",&iii);
  Dbx_Barrier(comm);
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate;ig++){
         printf("bck1 %d %d %g\n",ig,is,c[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Dbx_Barrier(comm);
  }/*endfor*/
#endif
  
/*========================================================================*/
/* II) Send the transformed position data                               */

  sendcounts = nstate_ncoef_proc_max*nstate_proc_max;
  recvcounts = nstate_ncoef_proc_max*nstate_proc_max;

  c_pt = c+1;
  c_temp_pt = c_temp+1;
  Alltoall(c_pt,sendcounts,MPI_DOUBLE,c_temp_pt,recvcounts,
                MPI_DOUBLE,comm);

#ifdef DEBUG
  if(myid==0){
    printf("This is c_temp in trans bck\n");
  }/*endif*/
  Dbx_Barrier(comm);
  for(is=1;is<=num_proc;is++){
    if(myid==is-1){
      for(ig=1;ig<=nstate_ncoef_proc_max*nstate_max;ig++){
         printf("bck2 %d %d %g\n",ig,is,c_temp[ig]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&ig);
    Dbx_Barrier(comm);
  }/*endfor*/
#endif


/*========================================================================*/
/* III) Extract the transformed position data                               */

  for(i=1;i<=nstate_ncoef_proc_max*nstate_max;i++){c[i]=0.0;}

  proc_rem = ncoef % num_proc;
  for(is=1;is<=nstate_proc;is++){
    ioff = 0;
    ioff_c = (is-1)*ncoef;
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (is-1)*nstate_ncoef_proc_max
                + (iproc-1)*(nstate_proc_max*nstate_ncoef_proc_max);
      nstate_ncoef_proc_now = nstate_ncoef_proc_max;
      if((iproc>proc_rem)&&(proc_rem>0)) nstate_ncoef_proc_now--;
      for(ig=1;ig<=nstate_ncoef_proc_now;ig++){
        itemp = ig+ioff_temp;
        i     = ig+ioff;
        c[(i+ioff_c)] = c_temp[itemp];
      }/*endfor*/   
      ioff += nstate_ncoef_proc_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is c in trans bck\n");
  }/*endif*/
  Dbx_Barrier(comm);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(is=1;is<=nstate_proc;is++){
       ioff_c = (is-1)*ncoef;
       for(ig=1;ig<=ncoef;ig++){
        printf("bck3 %d %d %d %g\n",iproc,is,ig,c[(ig+ioff_c)]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    Dbx_Barrier(comm);
    if(myid==0){printf("completed back transform\n");
    scanf("%d",&iii);
    Dbx_Barrier(comm);
  }/*endfor*/
  if(myid==0){printf("completed back transform\n");
    scanf("%d",&iii);
  Dbx_Barrier(comm);
#endif

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_therm_tran_pack(CPTHERM_POS *cptherm_pos,CPTHERM_INFO *cptherm_info,
                          int ioff_therm_safe,int pi_beads)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

 int ip,ich,k;
 int len_c_nhc      = cptherm_info->len_c_nhc;
 int num_c_nhc_norm = cptherm_info->num_c_nhc_norm;

/*======================================================================*/

 for(ip=1;ip<=pi_beads;ip++){
   for(ich=1;ich <= len_c_nhc;ich++){
     for(k=num_c_nhc_norm;k>=1;k--){
       cptherm_pos[ip].vc_nhc[ich][k] = 
                cptherm_pos[ip].vc_nhc[ich][(k+ioff_therm_safe)];
     }/* endfor */
   }/* endfor */
  }/* endfor */

/*========================================================================*/
 }/*end routine*/
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_therm_tran_unpack(CPTHERM_POS *cptherm_pos,CPTHERM_INFO *cptherm_info,
                          int ioff_therm_safe,int pi_beads)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

 int ip,ich,k;
 int num_c_nhc_proc = cptherm_info->num_c_nhc_proc;
 int len_c_nhc      = cptherm_info->len_c_nhc;

/*======================================================================*/

 for(ip=1;ip<=pi_beads;ip++){
   for(ich=1;ich <= len_c_nhc;ich++){
     for(k=num_c_nhc_proc;k>=1;k--){
       cptherm_pos[ip].vc_nhc[ich][k] = 
                cptherm_pos[ip].vc_nhc[ich][(k+ioff_therm_safe)];
     }/* endfor */
   }/* endfor */
  }/* endfor */

/*========================================================================*/
 }/*end routine*/
/*========================================================================*/



