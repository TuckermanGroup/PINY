/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/
/*                                                                        */
/*                       PI_MD:                                           */
/*           The future of simulation technology                          */
/*           ------------------------------------                         */
/*                 Module: energy_control.c                               */
/*                                                                        */
/* This routine calls the required force and PE routines                  */
/*                                                                        */
/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_math.h"

/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

void energy_control_intra(CLASS *class, BONDED *bonded,
                          GENERAL_DATA *general_data)

/*========================================================================*/
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  int i,iii,iver_get,ires_bond,ip;
  double vbond,vbend,vbend_bnd,vtors,vonfo,vonfo_vdw,vonfo_coul;
  double vbend_bnd_bond, vbend_bnd_bend;
  double vbond_free,vbend_free,vtors_free,vbar_free,pold;
  double vwatts_bond,vwatts_bend,vwatts_tot;
  double press_intra;

  int bond_npow         = bonded->bond.npow;
  int bond_nfree        = bonded->bond_free.num;
  int bond_nwatts_33    = bonded->grp_bond_watts.num_33;
  int bend_npow         = bonded->bend.npow;
  int bend_nfree        = bonded->bend_free.num;
  int bend_bnd_npow     = bonded->bend_bnd.num;
  int tors_npow         = bonded->tors.npow;
  int tors_nfree        = bonded->tors_free.num;
  int onfo_num          = bonded->onfo.num;
  int rbar_sig_nfree    = bonded->rbar_sig_free.nfree;

  double cell_vol       = general_data->cell.vol;

  int nres_tra          = general_data->timeinfo.nres_tra;
  int nres_tor          = general_data->timeinfo.nres_tor;
  int simopts_md        = general_data->simopts.md;
  int simopts_pimd      = general_data->simopts.pimd;
  double *pvten_tot     = general_data->ptens.pvten_tot;

  int iget_full_inter   = class->energy_ctrl.iget_full_inter;
  int iget_res_inter    = class->energy_ctrl.iget_res_inter;
  int iget_full_intra   = class->energy_ctrl.iget_full_intra;
  int iget_res_intra    = class->energy_ctrl.iget_res_intra;
  int iget_pe_real_inter= class->energy_ctrl.iget_pe_real_inter;
  int iget_pv_real_inter= class->energy_ctrl.iget_pv_real_inter;

  int np_forc           = class->communicate.np_forc;

  int pi_beads          = class->clatoms_info.pi_beads;  
  int pi_beads_proc     = class->clatoms_info.pi_beads_proc;  
  int *ip_lab           = class->clatoms_info.ip_lab;

  int pimd_on;
  pimd_on               =  general_data->simopts.pimd
                          +general_data->simopts.debug_pimd
                          +general_data->simopts.cp_pimd
                          +general_data->simopts.cp_wave_pimd
                          +general_data->simopts.cp_wave_min_pimd
                          +general_data->simopts.debug_cp_pimd;

/*======================================================================*/
/* I) Zero stuff and set some constants */

  vbond           = 0.0;
  vbend           = 0.0;
  vbar_free       = 0.0;
  vbend_bnd       = 0.0;
  vbend_bnd_bond  = 0.0;
  vbend_bnd_bend  = 0.0;
  vtors           = 0.0;
  vonfo           = 0.0;
  vonfo_vdw       = 0.0;
  vonfo_coul      = 0.0;
  vbond_free      = 0.0;
  vbend_free      = 0.0;
  vtors_free      = 0.0;
  vwatts_bond     = 0.0;
  vwatts_bend     = 0.0;
  vwatts_tot      = 0.0;

  iver_get = 0;
  if((pimd_on==1)&&(iget_full_inter==1)){iver_get = 1;}

  pold        = pvten_tot[1] + pvten_tot[5] + pvten_tot[9];
  if(pimd_on==1){get_vir_press(class,general_data,&pold);}

  ires_bond = 0;  
  if( (simopts_md==1) || (simopts_pimd==1) ){
    if( (nres_tra*nres_tor) > 1 ){ires_bond = 1;}
  }/*endif*/

/*======================================================================*/
/* II) Get intramolecular bond force and  PE         */

  if( (iget_full_intra==1) || (iget_res_intra==1) ){

    if(bond_npow!=0){
      if( (iget_full_inter==1) || (iget_res_inter==1) ){
        for(ip=1;ip<=pi_beads_proc;ip++){
            bond_both(&(class->clatoms_info),&(class->clatoms_pos[ip]),
              &(bonded->bond),&(general_data->cell),
              &(bonded->intra_scr),&(general_data->ptens),&vbond,iver_get,
              &(class->class_comm_forc_pkg),iget_pv_real_inter);
	}/*endfor*/
      }else{
        for(ip=1;ip<=pi_beads_proc;ip++){
            bond(&(class->clatoms_info),&(class->clatoms_pos[ip]),
              &(bonded->bond),&(general_data->cell),
              &(bonded->intra_scr),&(general_data->ptens),&vbond,iver_get,
              ires_bond,&(class->class_comm_forc_pkg),iget_pv_real_inter);
        }/*endfor*/
      }/*endif*/
      vbond/=pi_beads;
    }/*endif*/


    if( (bond_nfree!=0) && (pimd_on==0) ){
      bond_free(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->bond_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vbond_free,
                &(class->energy_ctrl),np_forc);
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/* III) Get mean-std umbrella sampling potential                        */

  if( (iget_full_intra==1) || (iget_res_intra==1) ){

    if( (rbar_sig_nfree!=0) && (pimd_on==0) ){
      rbar_sig_free(&(class->clatoms_info),&(class->clatoms_pos[1]),
                    &(bonded->rbar_sig_free),&(general_data->cell),
                    &(bonded->intra_scr),&(general_data->ptens),&vbar_free, 
                    &(class->energy_ctrl),np_forc);
    }/*endif*/

  }/*endif*/

/*======================================================================*/
/* IV) Get intramolecular group bond Watts                               */

  if( (iget_full_intra==1) || (iget_res_intra==1) ){

    if(bond_nwatts_33!=0){
        for(ip=1;ip<=pi_beads_proc;ip++){
            bond_watts_33(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                &(bonded->grp_bond_watts),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vwatts_bend,
                &vwatts_bond,&vwatts_tot,iver_get,
                &(class->class_comm_forc_pkg),iget_pv_real_inter);
        }/*endfor*/
      vwatts_bend/=pi_beads;
      vwatts_bond/=pi_beads;
      vwatts_tot /=pi_beads;
    }/*endif*/

  }/*endif*/


/*=======================================================================*/
/* V) Get intramolecular bend force and  PE          */

  if( (iget_full_intra==1) || (iget_res_intra==1) ){  

    if(bend_npow!=0){
     for(ip=1;ip<=pi_beads_proc;ip++){
         bend(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->bend),&(general_data->cell),
            &(bonded->intra_scr),&(general_data->ptens),&vbend,iver_get,
            &(class->class_comm_forc_pkg),iget_pv_real_inter);
     }/*endfor*/
    }/*endif*/
    vbend/=pi_beads;

    if( (bend_nfree!=0) && (pimd_on==0) ){
          bend_free(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->bend_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vbend_free,
                &(class->energy_ctrl),np_forc);
   }/*endif*/

  }/*endif*/

/*=======================================================================*/
/* VI) Get intramolecular Uri-Bradley force and  PE                      */

  if( (iget_full_intra==1) || (iget_res_intra==1) ){  

    if(bend_bnd_npow!=0){
     for(ip=1;ip<=pi_beads_proc;ip++){
         bend_bnd(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                &(bonded->bend_bnd),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vbend_bnd,
                &vbend_bnd_bend,&vbend_bnd_bond,iver_get,
                &(class->class_comm_forc_pkg),iget_pv_real_inter);
     }/*endfor*/
     vbend_bnd/=pi_beads;
     vbend_bnd_bend/=pi_beads;
     vbend_bnd_bond/=pi_beads;
    }/*endif*/

  }/*endif*/

/*=====================================================================*/
/* VII) Get intramolecular tors force and  PE                           */

  if(iget_full_intra==1){

    if(tors_npow!=0){
     for(ip=1;ip<=pi_beads_proc;ip++){
         tors(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->tors),&(general_data->cell),
            &(bonded->intra_scr),&(general_data->ptens),&vtors,iver_get,
            &(class->class_comm_forc_pkg),iget_pv_real_inter);
     }/*endfor*/
     vtors/=pi_beads;
    }/*endif*/

    if( (tors_nfree!=0) && (pimd_on==0)){
      tors_free(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(bonded->tors_free),&(general_data->cell),
                &(bonded->intra_scr),&(general_data->ptens),&vtors_free,
                &(class->energy_ctrl),np_forc);
    }/*endif*/

  }/*endif*/


/*======================================================================*/
/* VIII) Get intramolecular onefour force and  PE          */
  

  if( (iget_full_intra==1) || (iget_res_intra==1) ){  

    if(onfo_num!=0){
     for(ip=1;ip<=pi_beads_proc;ip++){
         onfo(&(class->clatoms_info),&(class->clatoms_pos[ip]),
            &(bonded->onfo),&(general_data->cell),
            &(bonded->intra_scr),&(general_data->ptens),&vonfo,
            &vonfo_vdw,&vonfo_coul,iver_get,&(class->class_comm_forc_pkg),
            iget_pv_real_inter);
     }/*endfor*/
     vonfo/=pi_beads;
     vonfo_coul /=(pi_beads);
     vonfo_vdw  /=(pi_beads);
    }/*endif*/

  }/*endif*/


/*======================================================================*/
/* IX) Collect and store PE                                           */
  
  if( (iget_full_intra==1) || (iget_full_inter==1) ){

    if(iget_pe_real_inter==1){
      general_data->stat_avg.vvdw       += vonfo_vdw;
      general_data->stat_avg.vcoul      += vonfo_coul;
      general_data->stat_avg.vintert    += (vonfo_vdw+vonfo_coul);
    }/*endif*/

    general_data->stat_avg.vbondt         = vbond;
    general_data->stat_avg.vbendt         = vbend;
    general_data->stat_avg.vbondt_watts   = vwatts_bond;
    general_data->stat_avg.vbendt_watts   = vwatts_bend;
    general_data->stat_avg.vtot_watts     = vwatts_tot;
    general_data->stat_avg.vbend_bndt     = vbend_bnd;
    general_data->stat_avg.vbend_bnd_bend = vbend_bnd_bend;
    general_data->stat_avg.vbend_bnd_bond = vbend_bnd_bond;
    general_data->stat_avg.vtorst         = vtors;
    general_data->stat_avg.vbond_free     = vbond_free;
    general_data->stat_avg.vbend_free     = vbend_free;
    general_data->stat_avg.vtors_free     = vtors_free;
    general_data->stat_avg.vonfot         = vonfo;
    general_data->stat_avg.vbar_free      = vbar_free;
    general_data->stat_avg.vintrat        = (vbond+vbend+vbend_bnd+vtors
                                            +vbond_free+vbend_free+vtors_free
                                            +vbar_free+vwatts_tot);
  }/*endif*/

  /*----------------------------------------------------------------------*/
  /*   Get intramolecular contribution to the pressure                    */
  
  if(iget_full_inter==1){

    if(pimd_on==0){
     if(iget_pv_real_inter==1){
      general_data->stat_avg.press_intra = 
              (pvten_tot[1]+ pvten_tot[5]+pvten_tot[9]-pold)/(3.0*cell_vol);
     }/*endif*/
    }else{
      get_vir_press(class,general_data,&press_intra);
      general_data->stat_avg.press_intra =  (press_intra-pold)/(cell_vol);
    }/*endif*/

  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/
