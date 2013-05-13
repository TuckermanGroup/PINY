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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

void energy_control_final(CLASS *class, BONDED *bonded,
                          GENERAL_DATA *general_data)

/*========================================================================*/
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  
  int i,iii,igloc,iver_get,ip,iup,j,igo;
  double *fx,*fy,*fz,temp;

  double *fxtemp      = class->ewd_scr.x;
  double *fytemp      = class->ewd_scr.y;
  double *fztemp      = class->ewd_scr.z;

  int nfreeze         = class->atommaps.nfreeze;
  int *freeze_map     = class->atommaps.freeze_map;
  int pimd_freez_typ  = class->atommaps.pimd_freez_typ;

  int nfree           = class->clatoms_info.nfree;
  int pi_beads_proc   = class->clatoms_info.pi_beads_proc;  
  int pi_beads_proc_st= class->clatoms_info.pi_beads_proc_st;  
  int *recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
  int myatm_start     = class->clatoms_info.myatm_start;
  int myatm_end       = class->clatoms_info.myatm_end;

  int np_beads        = class->communicate.np_beads;
  int np_states       = class->communicate.np_states;
  int np_tot          = class->communicate.np;
  int np_forc         = class->communicate.np_forc;
  int myid            = class->communicate.myid_forc;
  int myid_state      = class->communicate.myid_state;
  int myid_bead       = class->communicate.myid_bead;
  MPI_Comm comm_forc  = class->communicate.comm_forc;
  MPI_Comm world      = class->communicate.world;
  MPI_Comm comm_beads = class->communicate.comm_beads;

  int iget_full_inter = class->energy_ctrl.iget_full_inter;
  int iget_res_inter  = class->energy_ctrl.iget_res_inter;

  double *pvten_temp  = general_data->ptens.pvten_tmp;
  double *pvten_tot   = general_data->ptens.pvten_tot;
  double *pvten       = general_data->ptens.pvten;
  double t_ext        = general_data->statepoint.t_ext;
  double volume       = general_data->cell.vol;
  int hmat_int_typ    = general_data->cell.hmat_int_typ;

  int npt_i           = general_data->ensopts.npt_i;
  int npt_f           = general_data->ensopts.npt_f;
  int pimd_on;
  pimd_on             =  general_data->simopts.pimd
                        +general_data->simopts.debug_pimd
                        +general_data->simopts.cp_pimd
                        +general_data->simopts.cp_wave_pimd
                        +general_data->simopts.cp_wave_min_pimd
                        +general_data->simopts.debug_cp_pimd;
  

/*======================================================================*/
/* I) Set virial estimator flag stuff                                   */

  iver_get = 0;
  if( (pimd_on==1) && (iget_full_inter==1) ){iver_get = 1;}

/*======================================================================*/
/* II)  Get kinetic contribution to the pressure                        */
  
  if(iget_full_inter==1){
    general_data->stat_avg.press_kin = (nfree*t_ext)/(3.0*BOLTZ*volume);
  }/*endif*/

/*======================================================================*/
/* III) Get baro/par_rahman forces                                      */

/*----------------------------------------------------------------------*/
/*  A) Reduce the pressure tensor */

  if(np_tot > 1){

    if((np_beads>1)||((pimd_on==1)&&(np_states>1||np_forc>1))){
      for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
      Allreduce(&(pvten_tot[1]),&(pvten_temp[1]),9,MPI_DOUBLE,MPI_SUM,0,world);
      for(i=1;i<=9;i++){pvten_tot[i] = pvten_temp[i];}

      if(np_states>1||np_forc>1){
        for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
        Allreduce(&(pvten[1]),&(pvten_temp[1]),9,MPI_DOUBLE,MPI_SUM,0,world);
        for(i=1;i<=9;i++){pvten[i] = pvten_temp[i];}
      }else{
        if(myid_state==0){
          for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
          Allreduce(&(pvten[1]),&(pvten_temp[1]),9,MPI_DOUBLE,MPI_SUM,0,
                                                              comm_beads);
          for(i=1;i<=9;i++){pvten[i] = pvten_temp[i];}
        }/*endif*/
      }/*endif : np_states>1*/

    }else{

      if( (np_forc>1) || (np_states>1) ){
        for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
        Allreduce(&(pvten[1]),&(pvten_temp[1]),9,MPI_DOUBLE,MPI_SUM,0,world);
        for(i=1;i<=9;i++){pvten[i] = pvten_temp[i];}
      }/*endif : np_forc*/

    }/*endif : np_beads*/

  }/* endif np_tot */

/*----------------------------------------------------------------------*/
/* B) Get forces on the h-matrix                                        */

  if(myid_state == 0 || np_forc > 0){

    if(npt_i==1){
     box_force_iso(general_data->par_rahman.fgmat_p,
                 &(general_data->baro.f_lnv_p),
                   general_data->ptens.pvten,general_data->ptens.pvten_inc, 
                   general_data->statepoint.pext,general_data->cell.vol, 
                   bonded->constrnt.iconstrnt);
    }/*endif : isotropic box*/

    if(npt_f==1){
      box_force_flex(general_data->par_rahman.fgmat_p,
                     general_data->ptens.pvten,general_data->ptens.pvten_inc, 
                     general_data->statepoint.pext, 
                     general_data->statepoint.stens_ext, 
                     bonded->constrnt.iconstrnt,&(general_data->cell));

    }/*endif : flexible box*/

  }/*endif*/ 

/*========================================================================*/
/* IV) Distribute ghost forces into atom forces                          */

  if(myid_state == 0 || np_forc > 0){

    for(ip=1;ip<=pi_beads_proc;ip++){
      distrib_ghost_force(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                          &(class->ghost_atoms),iver_get);
    }/*endfor*/

  }/*endif*/

/*========================================================================*/
/* V) Zero forces of frozen atoms if any                                */

  if(myid_state == 0 || np_forc > 1){
    if(nfreeze > 0) {
      iup = 1; igo=0;
      if(pimd_freez_typ==2){iup=pi_beads_proc;igo=1;}
      if(pimd_freez_typ==1 && pi_beads_proc_st==1){igo=1;}
      if(igo==1){
        for(ip=1;ip<=iup;ip++){
          fx = class->clatoms_pos[ip].fx;
          fy = class->clatoms_pos[ip].fy;
          fz = class->clatoms_pos[ip].fz;
          for(i=1;i <= nfreeze;i++){
            igloc     = freeze_map[i];
            fx[igloc] = 0.0;
            fy[igloc] = 0.0;
            fz[igloc] = 0.0;
          }/*endfor*/
        }/*endfor*/
      }/*endif*/
    }/*endif*/

  }/* endif */


/*========================================================================*/
/* VI) Allreduce forces if force level par is on                          */

  if( (iget_full_inter+iget_res_inter)>=1){

    if(np_forc>1&&pimd_on==0){
      for(ip=1;ip<=pi_beads_proc;ip++){
        fx = class->clatoms_pos[ip].fx;
        fy = class->clatoms_pos[ip].fy;
        fz = class->clatoms_pos[ip].fz;
        Reduce_scatter(&fx[1],&fxtemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        Reduce_scatter(&fy[1],&fytemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        Reduce_scatter(&fz[1],&fztemp[myatm_start],&recv_count_atm[1],
                       MPI_DOUBLE,MPI_SUM,comm_forc);
        for(i=myatm_start;i<=myatm_end;i++){
          fx[i] = fxtemp[i];
          fy[i] = fytemp[i];
          fz[i] = fztemp[i];
        }/*endfor : i*/
      }/*endfor : ip*/
    }/*endif : np_forc>1*/

  }/*endif : inter RESPA step*/

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*=======================================================================*/




/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/


void box_force_iso(double *fgmat_p,double *f_lnv_p,double *pvten,
                        double *pvten_inc, 
                        double pext, double vol, int iconstrnt)

/*========================================================================*/
   {/*Begin Routine*/
/*=======================================================================*/

int iii;

/*=======================================================================*/

   *f_lnv_p = (pvten[1]+pvten[5]+pvten[9])-3.0*pext*vol;
   if(iconstrnt==1){
      *f_lnv_p += (pvten_inc[1] + pvten_inc[5] + pvten_inc[9]);
   }/*endif*/

/*-----------------------------------------------------------------------*/
}/*end routine */
/*=======================================================================*/



/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

void box_force_flex(double *fgmat_p, double *pvten,double *pvten_inc,
                    double pext, double stens_ext, int iconstrnt, CELL *cell)

/*========================================================================*/
/*  Begin Routine                                                         */
   {/*Begin Routine*/
/*========================================================================*/
/*         Local Variable declarations                                    */

  int i;

  double vol        = cell->vol;
  double area       = cell->area;
  double *hmat      = cell->hmat;
  int iperd         = cell->iperd;
  int hmat_cons_typ = cell->hmat_cons_typ;
  int hmat_int_typ  = cell->hmat_int_typ;

/*=======================================================================*/

  area = hmat[1]*hmat[5] - hmat[2]*hmat[4];

  for(i=1;i<=9;i++){fgmat_p[i]=pvten[i];}
  if(iconstrnt==1){for(i=1;i<=9;i++){fgmat_p[i]+=pvten_inc[i];}}

  fgmat_p[1]-= (pext*vol - stens_ext*area);
  fgmat_p[5]-= (pext*vol - stens_ext*area);
  fgmat_p[9]-= (pext*vol - stens_ext*area);

  constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,fgmat_p);

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*=======================================================================*/




