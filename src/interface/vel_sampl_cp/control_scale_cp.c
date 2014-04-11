/*==================================  ===============================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                Module: control_vx_smpl.c                          */
/*                                                                   */
/* Control routine for velocity resampling                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#define DEBUG_OFF



/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vc_scale(CP *cp)

/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*  Local Variables */
#include "../typ_defs/typ_mask.h"

  int i,ip,is,icoef,iii;
  double ake_tmp,ake,sc,temp_now;

/* Local pointers */

  double *vcre_up,*vcim_up,*vcre_dn,*vcim_dn;
  int nscale           = cp->cpcoeffs_info.cp_nfree;
  int nscale_up        = cp->cpcoeffs_info.cp_nfree_up;
  int nscale_dn        = cp->cpcoeffs_info.cp_nfree_dn;
  int cp_lsda          = cp->cpopts.cp_lsda;
  int  pi_beads_proc   = cp->cpcoeffs_info.pi_beads_proc;
  int nstate_up        = cp->cpcoeffs_info.nstate_up;
  int nstate_dn        = cp->cpcoeffs_info.nstate_dn;
  int ncoef_tot        = cp->cpcoeffs_info.ncoef;
  int ncoef_use        = cp->cpcoeffs_info.ncoef;
  int icmoff_up        = cp->cpcoeffs_info.icoef_start_up-1;
  int icmoff_dn        = cp->cpcoeffs_info.icoef_start_dn-1;
  double te_ext        = cp->cpopts.te_ext;
  double *cmass        = cp->cpcoeffs_info.cmass;
  int *ioff_upt        = cp->cpcoeffs_info.ioff_upt;
  int *ioff_dnt        = cp->cpcoeffs_info.ioff_dnt;
  int myid_state       = cp->communicate.myid_state;
  int myid             = cp->communicate.myid;
  int np_states        = cp->communicate.np_states;
  MPI_Comm comm_states = cp->communicate.comm_states;
  

/*===================================================================*/
/* 0) Write to screen                                                */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Scaling coefficient velocities\n");
    PRINT_LINE_DASH;printf("\n");
  }/*endif : error_check_on*/

/*====================================================================*/
/* I) Scale up state coef velocities                                  */

  if(np_states>1){
   for(ip=1;ip<=pi_beads_proc;ip++){
    if(cp->cpcoeffs_pos[ip].ivcoef_form_up!=1){
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Up Coef velocities are not in transposed form \n");
     printf("on state processor %d in control_scale_cp \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/
   }/*endfor*/
  }/*endif*/

  if(np_states>1){ncoef_use= cp->cp_comm_state_pkg_up.nstate_ncoef_proc;}
  for(ip=1;ip<=pi_beads_proc;ip++){

    vcre_up = cp->cpcoeffs_pos[ip].vcre_up;
    vcim_up = cp->cpcoeffs_pos[ip].vcim_up;
    ake_tmp = 0.0;
    for(is=1;is<=nstate_up;is++) {
      for(i=1;i<=ncoef_use;i++) {
        icoef = i+ioff_upt[is];
        ake_tmp += (vcre_up[icoef]*vcre_up[icoef]*cmass[(i+icmoff_up)]);
        ake_tmp += (vcim_up[icoef]*vcim_up[icoef]*cmass[(i+icmoff_up)]);
      }/*endfor:coeffs*/
    }/*endfor:states*/
    ake = ake_tmp;
    if(np_states>1){
     ake = 0.0;
     Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
    }/*endif*/
    temp_now = (BOLTZ*ake)/((double)(nscale_up));
    sc = sqrt(te_ext/temp_now);

    ake_tmp = 0.0;
    for(is=1;is<=nstate_up;is++) {
      for(i=1;i<=ncoef_use;i++) {
        icoef = i+ioff_upt[is];
        vcre_up[icoef] *= sc;
        vcim_up[icoef] *= sc;
        ake_tmp += (vcre_up[icoef]*vcre_up[icoef]*cmass[(i+icmoff_up)]);
        ake_tmp += (vcim_up[icoef]*vcim_up[icoef]*cmass[(i+icmoff_up)]);
      }/*endfor:coeffs*/
    }/*endfor:states*/

#ifdef DEBUG
    ake = ake_tmp;
    if(np_states>1){
     ake = 0.0;
     Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
    }
    temp_now = (BOLTZ*ake)/((double)(2*nscale_up));
    printf("temp_now %d %g\n",ip,temp_now);
#endif

  }/*endfor:beads*/

/*====================================================================*/
/* I) Scale dn state cp velocities                                    */

  if(cp_lsda==1 && nstate_dn != 0){


    if(np_states>1){
     for(ip=1;ip<=pi_beads_proc;ip++){
       if(cp->cpcoeffs_pos[ip].ivcoef_form_dn!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Dn Coef velocities are not in transposed form \n");
        printf("on state processor %d in control_scale_cp \n",myid_state);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
       }/*endif*/
     }/*endfor*/
    }/*endif*/

    if(np_states>1){ncoef_use= cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;}
    for(ip=1;ip<=pi_beads_proc;ip++){

      vcre_dn = cp->cpcoeffs_pos[ip].vcre_dn;
      vcim_dn = cp->cpcoeffs_pos[ip].vcim_dn;
      ake_tmp = 0.0;
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_use;i++) {
          icoef = i+ioff_dnt[is];
          ake_tmp += (vcre_dn[icoef]*vcre_dn[icoef]*cmass[(i+icmoff_dn)]);
          ake_tmp += (vcim_dn[icoef]*vcim_dn[icoef]*cmass[(i+icmoff_dn)]);
        }/*endfor:coeffs*/
      }/*endfor:states*/
      if(np_states>1){
       Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      } else {
       ake = ake_tmp;
      }
      temp_now = (BOLTZ*ake)/((double)(nscale_dn));
      sc = sqrt(te_ext/temp_now);

      ake_tmp = 0.0;
      for(is=1;is<=nstate_dn;is++) {
        for(i=1;i<=ncoef_use;i++) {
          icoef = i+ioff_dnt[is];
          vcre_dn[icoef] *= sc;
          vcim_dn[icoef] *= sc;
          ake_tmp += (vcre_dn[icoef]*vcre_dn[icoef]*cmass[(i+icmoff_dn)]);
          ake_tmp += (vcim_dn[icoef]*vcim_dn[icoef]*cmass[(i+icmoff_dn)]);
        }/*endfor:coeffs*/
      }/*endfor:states*/

#ifdef DEBUG
      ake = ake_tmp;
      if(np_states>1){
       ake = 0.0;
       Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      }/*endif*/
      temp_now = (BOLTZ*ake)/((double)(2*nscale_dn));
      printf("temp_now %d %g\n",ip,temp_now);
#endif

    }/*endfor:ip*/

  }/*endif:cp_lsda*/


/*====================================================================*/
/* III) Write to screen                                               */

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Coefficient velocity scaling complete\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif : myid*/

/*====================================================================*/
  }/*end routine*/
/*====================================================================*/





/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vcnhc_scale(CP *cp)

/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*  Local Variables */
#include "../typ_defs/typ_mask.h"

 int ip,inhc,ichain,iii;
 double ake_tmp,ake,sc,temp_now;
 double **cmass_nhc;
 double **vc_nhc;
 double cmass_nhc_massiv;

/*  Local pointers */
 int pi_beads      = cp->cpcoeffs_info.pi_beads;
 int pi_beads_proc = cp->cpcoeffs_info.pi_beads_proc;
 int num_c_nhc     = cp->cptherm_info.num_c_nhc_proc;
 int num_c_nhc_tot = cp->cptherm_info.num_c_nhc;
 int len_c_nhc     = cp->cptherm_info.len_c_nhc;
 int massiv_flag   = cp->cptherm_info.massiv_flag;
 double *vc_nhc_tmp= cp->cpscr.cpscr_therm.coef_kin;
 
 int myid_state       = cp->communicate.myid_state;
 int myid             = cp->communicate.myid;
 int np_states        = cp->communicate.np_states;
 MPI_Comm comm_states = cp->communicate.comm_states;
 MPI_Comm world       = cp->communicate.world;
 double te_ext        = cp->cpopts.te_ext;
 double heat_fact     = cp->cptherm_info.cp_therm_heat_fact;

/*===================================================================*/
/* 0) Write to screen                                                */

 if(myid==0){
  PRINT_LINE_STAR;
  printf("Scaling coefficient NHC velocities\n");
  PRINT_LINE_DASH;printf("\n");
 }/*endif*/

/*====================================================================*/
/* I) Scale cp NHC velocities under massiv                            */

 if(massiv_flag==1){

  for(ip=1;ip<=pi_beads_proc;ip++){
      vc_nhc = cp->cptherm_pos[ip].vc_nhc;
      ake_tmp = 0.0;
      cmass_nhc_massiv = cp->cptherm_info.cmass_nhc_massiv;
      for(inhc=1;inhc<=num_c_nhc;inhc++){
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        ake_tmp += cmass_nhc_massiv*vc_nhc[ichain][inhc]
                               *vc_nhc[ichain][inhc];     
      }}/*endfor : ichain,inhc*/
      if(np_states>1){
       Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      } else {
       ake = ake_tmp;
      }
      temp_now = (BOLTZ*ake)/((double)(len_c_nhc*num_c_nhc_tot));
      sc = sqrt(te_ext*heat_fact/temp_now);
      ake_tmp = 0.0;
      for(inhc=1;inhc<=num_c_nhc;inhc++){
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        vc_nhc[ichain][inhc] *= sc;
        ake_tmp += cmass_nhc_massiv*vc_nhc[ichain][inhc]
                               *vc_nhc[ichain][inhc];     
      }}/*endfor : ichain,inhc*/
      if(np_states>1){
       Allreduce(&(ake_tmp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_states);
      } else {
       ake = ake_tmp;
      }
      temp_now = (BOLTZ*ake)/((double)(len_c_nhc*num_c_nhc_tot));
#ifdef DEBUG
      printf("temp_now %d %g\n",ip,temp_now);
#endif
    }/*endfor:beads*/

 }/*endif:massiv*/

/*====================================================================*/
/* II) All thermostats on all processors                              */

 if(massiv_flag==0){

  for(ip=1;ip<=pi_beads_proc;ip++){
      vc_nhc = cp->cptherm_pos[ip].vc_nhc;
  /*----------------------*/ 
  /* i)Scale on proc zero */
    if(myid_state==0){
      ake_tmp = 0.0;
      cmass_nhc = cp->cptherm_info.cmass_nhc;
      for(inhc=1;inhc<=num_c_nhc;inhc++){
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        ake_tmp += cmass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                          *vc_nhc[ichain][inhc];     
      }}/*endfor : ichain,inhc*/
      ake = ake_tmp;
      temp_now = (BOLTZ*ake)/((double)(len_c_nhc*num_c_nhc_tot));
      sc = sqrt(te_ext/temp_now);
      ake_tmp = 0.0;
      for(inhc=1;inhc<=num_c_nhc;inhc++){
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        vc_nhc[ichain][inhc] *= sc;
        ake_tmp += cmass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                           *vc_nhc[ichain][inhc];     
      }}/*endfor : ichain,inhc*/
    }/*endif:myid_state*/
  /*--------------------------*/
  /* ii) Bcast to other procs */
    if(np_states>1){
      for(inhc=1;inhc<=num_c_nhc;inhc++){
        if(myid_state==0){
          for(ichain=1;ichain<=len_c_nhc;ichain++){
             vc_nhc_tmp[ichain] = vc_nhc[ichain][inhc];
          }/*endfor : ichain*/
        }/*endif*/
        Barrier(comm_states);
        Bcast(&vc_nhc_tmp[1],len_c_nhc,MPI_DOUBLE,0,world);
        for(ichain=1;ichain<=len_c_nhc;ichain++){
          vc_nhc[ichain][inhc] = vc_nhc_tmp[ichain];
        }/*endfor : ichain*/
      }/*endfor : inhc*/
    }/*endif:np_states*/
  }/*endfor:beads*/ 
 
 }/*endif : massiv_flag*/

/*====================================================================*/
/* III) Write to screen                                               */

 if(myid==0){
   PRINT_LINE_DASH;
   printf("Coefficient velocity scaling complete\n");
   PRINT_LINE_STAR;printf("\n");
 }/*endif : error_check_on*/

/*====================================================================*/
}/*end routine*/
/*====================================================================*/
