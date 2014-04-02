/*==================================  ===============================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                Module: control_vcnhc_smpl.c                       */
/*                                                                   */
/* Control routine for velocity resampling                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
void control_vcnhc_smpl(CP *cp)
/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*        Local Variables */
#include "../typ_defs/typ_mask.h"

 int iii,i,is,icoef,inhc,ichain,ip;
 double **vc_nhc;

/*        Local pointers */
 int pi_beads       = cp->cpcoeffs_info.pi_beads_proc;
 int num_c_nhc      = cp->cptherm_info.num_c_nhc;
 int len_c_nhc      = cp->cptherm_info.len_c_nhc;
 int massiv_flag    = cp->cptherm_info.massiv_flag;
 double *vc_nhc_tmp = cp->cpscr.cpscr_therm.coef_kin;

 int myid             = cp->communicate.myid;
 MPI_Comm world       = cp->communicate.world;
 int np_states        = cp->communicate.np_states;
 int myid_state       = cp->communicate.myid_state;
 MPI_Comm comm_states = cp->communicate.comm_states;

/*===================================================================*/
/* 0) Write to screen                                                */

 if(myid==0){
   PRINT_LINE_STAR;
   printf("Sampling coefficient NHC velocities\n");
   PRINT_LINE_DASH;printf("\n");
 }/*endif*/

/*===================================================================*/
/* I) Sample coefficient extended system velocities under massiv     */

 if(massiv_flag==1){

   sampl_vcnhc(&(cp->cptherm_info),(cp->cptherm_pos),
               &(cp->cpscr),
               &(cp->cpopts),pi_beads,
               &(cp->vel_samp_cp.iseed),
               &(cp->vel_samp_cp.iseed2),&(cp->vel_samp_cp.qseed),
               &(cp->cpcoeffs_info),cp->communicate.np_states);

 }/*endif:massiv*/

/*===================================================================*/
/* I) Sample coefficient extended system velocities: no massiv       */

 if(massiv_flag==0){

   if(myid_state==0){
     sampl_vcnhc(&(cp->cptherm_info),(cp->cptherm_pos),
                 &(cp->cpscr),
                 &(cp->cpopts),pi_beads,
                 &(cp->vel_samp_cp.iseed),
                 &(cp->vel_samp_cp.iseed2),&(cp->vel_samp_cp.qseed),
                 &(cp->cpcoeffs_info),cp->communicate.np_states);
   }/*endif*/

   if(np_states>1){
     for(ip=1;ip<=pi_beads;ip++){
       vc_nhc = cp->cptherm_pos[ip].vc_nhc;
       for(ichain=1;ichain<=len_c_nhc;ichain++){
         if(myid_state==0){
           for(inhc=1;inhc<=num_c_nhc;inhc++){
             vc_nhc_tmp[inhc] = vc_nhc[ichain][inhc];
           }/*endfor : ichain*/
         }/*endif*/
         Barrier(comm_states);
         Bcast(&vc_nhc_tmp[1],num_c_nhc,MPI_DOUBLE,0,world);
          for(inhc=1;inhc<=num_c_nhc;inhc++){
           vc_nhc[ichain][inhc] = vc_nhc_tmp[inhc];
         }/*endfor : ichain*/
       }/*endfor : inhc*/
     }/*endif:ip*/
   }/*endif:npstates*/

 }/*endif:massiv*/

/*====================================================================*/
/* II) Write to screen                                               */

 if(myid==0){
   PRINT_LINE_DASH;
   printf("Coefficient NHC velocity sampling complete\n");
   PRINT_LINE_STAR;printf("\n");
 }/*endif*/

/*====================================================================*/
    }/*end routine*/
/*====================================================================*/


















