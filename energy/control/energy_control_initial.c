/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: energy_control.c                               */
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
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void energy_control_initial(CLASS *class, BONDED *bonded, 
                            GENERAL_DATA *general_data)
/*==========================================================================*/
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
#include "../typ_defs/typ_mask.h"

  double vol,area;
  int i,iii,iver_get,ip;
  int energy_ctrl_iget_pe_real_inter;
  int energy_ctrl_iget_pv_real_inter;

/*         Local Pointers                                                */

  int nres_tra                    = general_data->timeinfo.nres_tra;
  int nres_tor                    = general_data->timeinfo.nres_tor;
  int nres_ter                    = general_data->timeinfo.nres_ter;
  int nres_pimd                   = general_data->timeinfo.nres_pimd;

  int energy_ctrl_iget_full_inter = class->energy_ctrl.iget_full_inter;
  int energy_ctrl_iget_res_inter  = class->energy_ctrl.iget_res_inter;
  int energy_ctrl_int_res_ter     = class->energy_ctrl.int_res_ter;
  int energy_ctrl_itime           = class->energy_ctrl.itime;
  int energy_ctrl_iget_pe_real_inter_freq
                                  =class->energy_ctrl.iget_pe_real_inter_freq;

  int pi_beads                    = class->clatoms_info.pi_beads;  
  int pi_beads_proc               = class->clatoms_info.pi_beads_proc;  
  int natm_tot                    = class->clatoms_info.natm_tot;
  int atm_hess_calc               = class->clatoms_info.hess_calc;
  double *xtemp                   = class->ewd_scr.x;
  double *ytemp                   = class->ewd_scr.y;
  double *ztemp                   = class->ewd_scr.z;
  double *x,*y,*z;
  double *fx,*fy,*fz;
  double *fxt,*fyt,*fzt;
  double *hess_xx,*hess_xy,*hess_xz,*hess_yy,*hess_yz,*hess_zz;

  int iperd                       = general_data->cell.iperd;
  double *hmat                    = general_data->cell.hmat;
  double *hmati                   = general_data->cell.hmati;
  double *pvten                   = general_data->ptens.pvten;
  double *pvten_tot               = general_data->ptens.pvten_tot;

  int np_forc                     = class->communicate.np_forc;  
  int *recv_count_atm             = class->class_comm_forc_pkg.recv_count_atm;
  int *displs_atm                 = class->class_comm_forc_pkg.displs_atm;
  int myid_forc                   = class->class_comm_forc_pkg.myid;
  int myatm_start                 = class->clatoms_info.myatm_start;
  MPI_Comm comm_forc              = class->class_comm_forc_pkg.comm;

  int simopts_md                  = general_data->simopts.md;

  int cp_on,pimd_on;

  cp_on   = general_data->simopts.cp_min  +general_data->simopts.cp_wave_min
           +general_data->simopts.cp      +general_data->simopts.cp_wave
           +general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd
           +general_data->simopts.debug_cp+general_data->simopts.debug_cp_pimd
           +general_data->simopts.cp_wave_min_pimd;
  pimd_on = general_data->simopts.pimd
           +general_data->simopts.debug_pimd
           +general_data->simopts.cp_pimd
           +general_data->simopts.cp_wave_pimd
           +general_data->simopts.cp_wave_min_pimd
           +general_data->simopts.debug_cp_pimd;

/*======================================================================*/
/* 0) Calculate flags and volume stuff */

  iver_get = 0;
  if(pimd_on==1){if(energy_ctrl_iget_full_inter==1){iver_get = 1;}}

  energy_ctrl_iget_pe_real_inter = 0;
  energy_ctrl_iget_pv_real_inter = 0;
  if(energy_ctrl_iget_full_inter==1){
    energy_ctrl_iget_pe_real_inter = 1;
    energy_ctrl_iget_pv_real_inter = 1;
  }/*endif*/
  if(energy_ctrl_itime > 1 && cp_on == 0){
   if((energy_ctrl_itime % energy_ctrl_iget_pe_real_inter_freq)!=0){
    energy_ctrl_iget_pe_real_inter = 0;
    if(energy_ctrl_int_res_ter ==1){energy_ctrl_iget_pv_real_inter = 0;}
   }/*endif*/
  }/*endif*/

  gethinv(hmat,hmati,&vol,iperd);

  class->energy_ctrl.iget_pe_real_inter = energy_ctrl_iget_pe_real_inter;
  class->energy_ctrl.iget_pv_real_inter = energy_ctrl_iget_pv_real_inter;
  general_data->cell.vol                = vol;

  area               = hmat[1]*hmat[5] - hmat[2]*hmat[4];
  general_data->cell.area               = area;
  general_data->baro.area               = area;
  general_data->par_rahman.area         = area;

/*======================================================================*/
/* I) Zero stuff */


  general_data->stat_avg.vbondt         = 0.0;
  general_data->stat_avg.vbendt         = 0.0;
  general_data->stat_avg.vbend_bndt     = 0.0;
  general_data->stat_avg.vbend_bnd_bond = 0.0;
  general_data->stat_avg.vbend_bnd_bend = 0.0;
  general_data->stat_avg.vtorst         = 0.0;
  general_data->stat_avg.vonfot         = 0.0;
  general_data->stat_avg.vintrat        = 0.0;
  general_data->stat_avg.vbond_free     = 0.0;
  general_data->stat_avg.vbend_free     = 0.0;
  general_data->stat_avg.vtors_free     = 0.0;
  general_data->stat_avg.vrecip         = 0.0;
  general_data->stat_avg.vlong          = 0.0;
  if(energy_ctrl_iget_pe_real_inter==1){
   general_data->stat_avg.vintert       = 0.0;
   general_data->stat_avg.vvdw          = 0.0;
   general_data->stat_avg.vcoul         = 0.0;
   general_data->stat_avg.vreal         = 0.0;
  }/*endif*/

  if(cp_on==1){
    general_data->stat_avg.cp_eke   = 0.0;
    general_data->stat_avg.cp_enl   = 0.0;
    general_data->stat_avg.cp_exc   = 0.0;
    general_data->stat_avg.cp_muxc  = 0.0;
    general_data->stat_avg.cp_eext  = 0.0;
    general_data->stat_avg.cp_ehart = 0.0;
    general_data->stat_avg.vrecip   = 0.0;
  }/*endif*/

  for(ip=1;ip<= pi_beads_proc;ip++){
    fx = class->clatoms_pos[ip].fx;
    fy = class->clatoms_pos[ip].fy;
    fz = class->clatoms_pos[ip].fz;
    for(i=1;i<= natm_tot;i++){
      fx[i] = 0.0;
      fy[i] = 0.0;
      fz[i] = 0.0;
    }/*endfor*/
    if( (iver_get==1) && (pi_beads>1)){
      fxt = class->clatoms_pos[ip].fxt;
      fyt = class->clatoms_pos[ip].fyt;
      fzt = class->clatoms_pos[ip].fzt;
      for(i=1;i<= natm_tot;i++){
        fxt[i] = 0.0;
        fyt[i] = 0.0;
        fzt[i] = 0.0;
      }/*endfor*/
    }/*endif*/
    if(atm_hess_calc == 3){
      hess_xx = class->clatoms_pos[ip].hess_xx;
      hess_xy = class->clatoms_pos[ip].hess_xy;
      hess_xz = class->clatoms_pos[ip].hess_xz;
      hess_yy = class->clatoms_pos[ip].hess_yy;
      hess_yz = class->clatoms_pos[ip].hess_yz;
      hess_zz = class->clatoms_pos[ip].hess_zz;
      for(i=1;i<=natm_tot*natm_tot;i++){
        hess_xx[i] = 0.0;
        hess_xy[i] = 0.0;
        hess_xz[i] = 0.0;
        hess_yy[i] = 0.0;
        hess_yz[i] = 0.0;
        hess_zz[i] = 0.0;
      }/* endfor */
    }/* endif */
   }/*endfor*/

   for(i=1;i<=9;i++){pvten[i]=0.0;}
   if(energy_ctrl_iget_pv_real_inter==1){
     for(i=1;i<=9;i++){pvten_tot[i] = 0.0;}
   }/*endif*/

  
/*======================================================================*/
/* II) Initialize weights */

  class->for_scr.wght_ter         = 1.0;
  class->for_scr.wght_ter_res     = 1.0;
  bonded->intra_scr.wght_tra      = 1.0;
  bonded->intra_scr.wght_tra_res  = 1.0;
  bonded->intra_scr.int_res_inter = 0;

  if(simopts_md==1){
    class->for_scr.wght_ter        = (double)(nres_tra*nres_tor*nres_ter);
    class->for_scr.wght_ter_res    = (double)(nres_tra*nres_tor);
    bonded->intra_scr.wght_tra     = (double)(nres_tra);
    bonded->intra_scr.wght_tra_res = 1.0;
    if(nres_ter>1){bonded->intra_scr.int_res_inter = 1;}
  }/*endif*/

  if(pimd_on==1){
    bonded->intra_scr.int_res_inter = 0;
    if(nres_ter>1){bonded->intra_scr.int_res_inter = 1;}
  
    class->for_scr.wght_ter   = (double)(nres_pimd*nres_tra*nres_tor*nres_ter);
    class->for_scr.wght_ter_res    = (double)(nres_pimd*nres_tra*nres_tor);
    bonded->intra_scr.wght_tra     = (double)(nres_pimd*nres_tra);
    bonded->intra_scr.wght_tra_res = (double)(nres_pimd);
    class->clatoms_info.wght_pimd  = 1.0;
    class->for_scr.wght_pimd       = 1.0;
  }/*endif*/

  class->for_scr.wght_tra        = bonded->intra_scr.wght_tra;
  class->for_scr.wght_tra_res    = bonded->intra_scr.wght_tra_res;
  bonded->intra_scr.wght_ter     = class->for_scr.wght_ter;
  bonded->intra_scr.wght_ter_res = class->for_scr.wght_ter_res;

/*======================================================================*/
/* III) Allgatherv positions for force level parallel      */

  if(np_forc > 1){

    if((energy_ctrl_iget_full_inter==1)||(energy_ctrl_iget_res_inter ==1)){
      for(ip=1;ip<=pi_beads_proc;ip++){
        x = class->clatoms_pos[ip].x;
        y = class->clatoms_pos[ip].y;
        z = class->clatoms_pos[ip].z;
        Barrier(comm_forc); /* Needed by T3E */
        Allgatherv(&(x[myatm_start]),recv_count_atm[(myid_forc+1)],
                   MPI_DOUBLE,&(xtemp[1]),&recv_count_atm[1],&displs_atm[1],
                   MPI_DOUBLE,0,comm_forc);
        Allgatherv(&(y[myatm_start]),recv_count_atm[(myid_forc+1)],
                   MPI_DOUBLE,&(ytemp[1]),&recv_count_atm[1],&displs_atm[1],
                   MPI_DOUBLE,0,comm_forc);
        Allgatherv(&(z[myatm_start]),recv_count_atm[(myid_forc+1)],
                   MPI_DOUBLE,&(ztemp[1]),&recv_count_atm[1],&displs_atm[1],
                   MPI_DOUBLE,0,comm_forc);
        for(i=1;i<=natm_tot;i++){
          x[i] = xtemp[i];
          y[i] = ytemp[i];
          z[i] = ztemp[i];
        }/*endfor*/
      }/*endfor : pi_beads_proc*/ 
    }/*endif*/

  }/*endif*/

/*-----------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/
