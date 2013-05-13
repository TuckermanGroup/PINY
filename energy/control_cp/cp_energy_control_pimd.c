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
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_local.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_energy_ctrl_cp_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/

void cp_energy_control_pimd(CLASS *class,BONDED *bonded,
                             GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
/*        Begin Routine                                                  */
{/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
#include "../typ_defs/typ_mask.h"

  int i,ip,iii;
  double *fx,*fy,*fz,*fxm,*fym,*fzm,*fxt,*fyt,*fzt;
  double vbond_free,vbend_free,vtors_free;
  double akvirx;

  int natm_tot        = class->clatoms_info.natm_tot; 
  int pi_beads        = class->clatoms_info.pi_beads; 
  int pi_beads_proc   = class->clatoms_info.pi_beads_proc; 
  int myid            = class->communicate.myid; 
  int myid_bead       = class->communicate.myid_bead; 
  int myid_state      = class->communicate.myid_state; 
  double *pvten       = general_data->ptens.pvten;
  double *pvten_tot   = general_data->ptens.pvten_tot;
  int nproc           = class->communicate.np;
  MPI_Comm world      = class->communicate.world;
  MPI_Comm comm_beads = class->communicate.comm_beads;
  int simopts_cp      = general_data->simopts.cp;
  int simopts_cp_min  = general_data->simopts.cp_min;
  int simopts_cp_pimd = general_data->simopts.cp_pimd;
  int iget_full_inter = class->energy_ctrl.iget_full_inter;
  int iget_res_inter  = class->energy_ctrl.iget_res_inter;
  int iget_full_intra = class->energy_ctrl.iget_full_intra;
  int iget_res_intra  = class->energy_ctrl.iget_res_intra;
  double t_ext        = general_data->statepoint.t_ext;
  double volume       = general_data->cell.vol; 
  int bond_free_num      = bonded->bond_free.num;
  int bend_free_num      = bonded->bend_free.num;
  int tors_free_num      = bonded->tors_free.num;
  
/*======================================================================*/
/* I) Zero stuff */

  akvirx     = 0.0;
  vbond_free  = 0.0;
  vbend_free  = 0.0;
  vtors_free  = 0.0;

/*======================================================================*/
/* Initialize */

  energy_control_initial(class,bonded,general_data);


/*======================================================================*/
/* III) Get electron contribution to the atm force and total E          */
/*      Get forces on basis set parameters                              */

  energy_control_elec(class,bonded,general_data,cp);

/*======================================================================*/
/*======================================================================*/
/*  ONLY GET CLASSICAL FORCES AND ENERGY IF NECESSARY                   */

  if( (simopts_cp==1) || (simopts_cp_min==1) || (simopts_cp_pimd==1)){
  
    if(myid_state==0){
/*======================================================================*/
/* IV) Get intermolecular real space force and  PE   */
  
      energy_control_inter_real(class,bonded,general_data);

/*======================================================================*/
/* IV) Get surface PE                                                   */

      energy_control_surf(class,bonded,general_data);

/*======================================================================*/
/* VI) Get intramolecular bond force and  PE         */

      energy_control_intra(class,bonded,general_data);

/*==========================================================================*/
/* XIII) Construct staging/centroid forces                                  */

      for(ip=1;ip<= pi_beads_proc;ip++){
        fx = class->clatoms_pos[ip].fx;
        fy = class->clatoms_pos[ip].fy;
        fz = class->clatoms_pos[ip].fz;
        for(i=1;i<=natm_tot;i++){
          fx[i] /= pi_beads;
          fy[i] /= pi_beads;
          fz[i] /= pi_beads;
        }/*endfor*/
      }/*endfor*/

      if( (iget_full_inter==1)){
        if(pi_beads>1){
          for(ip=1;ip<= pi_beads_proc;ip++){
            fxt = class->clatoms_pos[ip].fxt;
            fyt = class->clatoms_pos[ip].fyt;
            fzt = class->clatoms_pos[ip].fzt;
            for(i=1;i<= class->clatoms_info.natm_tot;i++){
              fxt[i] /= pi_beads;
              fyt[i] /= pi_beads;
              fzt[i] /= pi_beads;
            }/*endfor*/
          }/*endfor*/
        }/*endif*/
      }/*endif*/

    }/*endif : myid_state == 0*/

    for(i=1;i<=9;i++){
      pvten[i]     /= pi_beads;
      pvten_tot[i] /= pi_beads;
    }/*endfor*/

  /*----------------------------------------------------------------------*/
  /*  Add in mode forces and correct pressure tensor                      */

    if(myid_state==0){
      if(pi_beads>1){pimd_pvten_forc_corr(class,general_data,&akvirx);}

      control_pimd_trans_force(class,general_data);

      for(ip=1;ip<=pi_beads_proc;ip++){
        fx = class->clatoms_pos[ip].fx;
        fy = class->clatoms_pos[ip].fy;
        fz = class->clatoms_pos[ip].fz;
        fxm = class->clatoms_pos[ip].fxm;
        fym = class->clatoms_pos[ip].fym;
        fzm = class->clatoms_pos[ip].fzm;
        for(i=1;i<=natm_tot;i++){
          fx[i] += fxm[i];
          fy[i] += fym[i];
          fz[i] += fzm[i];
        }/*endfor*/
      }/*endfor*/

  /*----------------------------------------------------------------------*/
  /*   Get kinetic contribution to the pressure                    */

     general_data->stat_avg.press_kin =  (natm_tot*t_ext)/(volume*BOLTZ);

    }/*endif*/

/*==========================================================================*/
/* XIII) Mode free energy                                                   */

    if(myid_state==0){
    if(myid_bead==0){
     if( (iget_full_intra==1) || (iget_res_intra==1) ){
      if(bond_free_num!=0){
        bond_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                       &(bonded->bond_free),&(general_data->cell),
                       &(bonded->intra_scr),&(general_data->ptens),&vbond_free,
                       &(class->energy_ctrl),class->communicate.np_forc);
      }/*endif*/
      if(bend_free_num!=0){
        bend_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                       &(bonded->bend_free),&(general_data->cell),
                       &(bonded->intra_scr),&(general_data->ptens),&vbend_free,
                       &(class->energy_ctrl),class->communicate.np_forc);
      }/*endif*/
      if(tors_free_num!=0){
        tors_free_mode(&(class->clatoms_info),&(class->clatoms_pos[1]),
                       &(bonded->tors_free),&(general_data->cell),
                       &(bonded->intra_scr),&(general_data->ptens),&vtors_free,
                       &(class->energy_ctrl),class->communicate.np_forc);
      }/*endif*/
      if((bond_free_num+bend_free_num+tors_free_num)!=0){
        distrib_ghost_force_mode(&(class->clatoms_info),
                                 &(class->clatoms_pos[1]),
                                 &(class->ghost_atoms));
      }/*endif*/
     }/*endif*/

    }/*endif : myid_bead=0*/
    }/*endif : myid_state==0*/

/*======================================================================*/
/* VI) Finish energy routine                                            */

    energy_control_final(class,bonded,general_data);

  }/*endif: Classical stuff */
  if(nproc>1){Barrier(world);}
/*======================================================================*/
/*======================================================================*/
/* XIII) Collect and store PE                                           */
  
  if( ((class->energy_ctrl.iget_full_inter)==1)){
    general_data->stat_avg.vbond_free = vbond_free;
    general_data->stat_avg.vbend_free = vbend_free;
    general_data->stat_avg.vtors_free = vtors_free;
    general_data->stat_avg.vintrat += vbond_free+vbend_free+vtors_free;
    general_data->stat_avg.pi_ke_vir  = akvirx; 
  }/*endif*/

/*-----------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/













