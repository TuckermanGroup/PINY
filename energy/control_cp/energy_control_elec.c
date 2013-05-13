/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: cp_energy_control.c                            */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void energy_control_elec(CLASS *class,BONDED *bonded,
                         GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

 double vrecip,vrecip_temp,vself,vbgr,vecor,vecor_temp;
 double *pvten_temp = general_data->ptens.pvten_tmp;
 int ip,iver_get,iii,i;
 int myid_state = cp->communicate.myid_state;
 int myid       = cp->communicate.myid;
 int pi_beads_proc = cp->cpcoeffs_info.pi_beads_proc;
 int pi_beads      = cp->cpcoeffs_info.pi_beads;
 int iperd         = general_data->cell.iperd;
 int np_beads      = class->communicate.np_beads;
 int np_forc       = class->communicate.np_forc;
 int myid_bead     = class->communicate.myid_bead;
  
/*======================================================================*/
/* 0)   Initialize local variables                                      */

  iver_get = 0;
  if( ((class->energy_ctrl.iget_full_inter)==1)&&(pi_beads>1)){
   iver_get = 1;
  }/*endif*/


  vself       = 0.0;
  vbgr        = 0.0;   
  vecor       = 0.0;
  vecor_temp  = 0.0;
  vrecip      = 0.0;
  general_data->stat_avg.vrecip     = 0.0;

/*======================================================================*/
/* I)   Get electronic energy and forces                                */

 for(ip=1;ip<=pi_beads_proc;ip++){
  cp_ks_energy_ctrl(cp,ip,&(general_data->ewald),&(class->ewd_scr),
                      &(general_data->cell),
                      &(class->clatoms_info),
                      &(class->clatoms_pos[ip]),
                      &(class->atommaps),&(general_data->stat_avg),
                      &(general_data->ptens),
                      &(general_data->simopts),
                      &(class->for_scr));
 }/* endif ip */


 general_data->stat_avg.cp_eke   /= pi_beads;
 general_data->stat_avg.cp_enl   /= pi_beads;
 general_data->stat_avg.cp_exc   /= pi_beads;
 general_data->stat_avg.cp_muxc  /= pi_beads;
 general_data->stat_avg.cp_eext  /= pi_beads;
 general_data->stat_avg.cp_ehart /= pi_beads;
 general_data->stat_avg.vrecip   /= pi_beads;


/*======================================================================*/
/* II)   Ewald self and background terms                                */

 if(myid == 0 && iperd > 0){
   ewald3d_selfbgr_cp(&(class->clatoms_info), &(general_data->ewald),
                      &(general_data->ptens), (general_data->cell.vol),
                      &vself,&vbgr,iperd);
 }/* endif myid_state */

/*======================================================================*/
/* III) Ewald corrections          */

 if( (myid_state == 0) || (np_forc > 1) ){
  if(((class->energy_ctrl.iget_full_intra)==1)||
     ((class->energy_ctrl.iget_res_intra)==1)){
    if((bonded->ecor.num)!=0){
      for(ip=1;ip<=pi_beads_proc;ip++){
       ecor(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->ecor),&(general_data->cell),
            &(bonded->intra_scr),&(general_data->ptens),&vecor,iver_get,
            &(class->class_comm_forc_pkg),
            class->energy_ctrl.iget_pe_real_inter,
            class->energy_ctrl.iget_pv_real_inter);
      }/*endfor*/
      vecor/=pi_beads;
    }/*endif*/
  }/*endif*/
 }/*endif myid_state*/


/*======================================================================*/
/* IV) Collect and store PE                                           */

 if(class->communicate.np_forc > 1){
  Reduce(&vecor, &vecor_temp,1,MPI_DOUBLE,MPI_SUM,0,
         class->communicate.comm_forc);
  vecor = vecor_temp;
 }

 if((class->communicate.np_states>1)||(class->communicate.np_beads>1)
        && iperd > 0){
   Allreduce(&(general_data->stat_avg.vrecip), &(vrecip_temp),1,MPI_DOUBLE,
               MPI_SUM,0,class->communicate.world);
   general_data->stat_avg.vrecip = vrecip_temp;
 }/*endif*/

 if( myid_state==0 ){
  if( ((class->energy_ctrl.iget_full_intra)==1)||
     ((class->energy_ctrl.iget_full_inter)==1)){

    vrecip                      = general_data->stat_avg.vrecip;
    vrecip                      += vself+vbgr+vecor;
    general_data->stat_avg.vcoul      += vrecip;
    general_data->stat_avg.vintert   += vrecip;
    general_data->stat_avg.vrecip     = vrecip;

  }/*endif*/
 }/*endif myid*/

/*======================================================================*/
}/* end routine */
/*======================================================================*/









