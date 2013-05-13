/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_dt2_to_dt_pimd                           */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nvt_pimd(CLASS *class,BONDED *bonded,
                            GENERAL_DATA *general_data,
                            int ir_tra,int ir_tor,int ir_ter,int ir_pimd,
                            double dt) 

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int ipart,iflag,ip;
   int ix_now;
   int iii,iproc;
   int iflag_mass;
   double dt2;
   double *class_clatoms_mass;
   double *class_clatoms_x;
   double *class_clatoms_y;
   double *class_clatoms_z;
   double *class_clatoms_vx;
   double *class_clatoms_vy;
   double *class_clatoms_vz;
   double *class_clatoms_fx;
   double *class_clatoms_fy;
   double *class_clatoms_fz;
   double *class_clatoms_fxm;
   double *class_clatoms_fym;
   double *class_clatoms_fzm;

   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int nres_tra         = (general_data->timeinfo.nres_tra);
   int nres_tor         = (general_data->timeinfo.nres_tor);
   int nres_ter         = (general_data->timeinfo.nres_ter);
   int nres_pimd        = (general_data->timeinfo.nres_pimd);
   int ix_respa         = general_data->timeinfo.ix_respa;
   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

   class->clatoms_info.wght_pimd = 1.0;
   dt2 = dt/2.0;

/*==========================================================================*/
/* 1) Evolve velocities                                                     */
 
   for(ip=1;ip<=pi_beads_proc;ip++){
     class_clatoms_vx   = class->clatoms_pos[ip].vx;
     class_clatoms_vy   = class->clatoms_pos[ip].vy;
     class_clatoms_vz   = class->clatoms_pos[ip].vz;
     class_clatoms_fx   = class->clatoms_pos[ip].fx;
     class_clatoms_fy   = class->clatoms_pos[ip].fy;
     class_clatoms_fz   = class->clatoms_pos[ip].fz;
     class_clatoms_mass = class->clatoms_pos[ip].mass;
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                  /class_clatoms_mass[ipart];
       class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                  /class_clatoms_mass[ipart];
       class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                  /class_clatoms_mass[ipart];
     }/*endfor*/
   }/*endfor*/

/*==========================================================================*/
/* 2) Evolve NHCs and velocities                                            */

   ix_now = 5;
   if((ir_pimd==nres_pimd)){ix_now=4;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)){ix_now=3;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)&&(ir_tor==nres_tor))
                            {ix_now=2;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)&&(ir_tor==nres_tor)
                                 &&(ir_ter==nres_ter)){ix_now=1;}

   if(ix_respa>=ix_now){
     if(pi_beads_proc_st == 1){
       iflag_mass = 1;
       apply_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                     &(class->therm_info_class),&(class->therm_class),
                     &(class->int_scr),iflag_mass,
                     &(class->class_comm_forc_pkg));
       iflag_mass = 2;
       for(ip=2;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
         apply_NHC_bead_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
#endif

         apply_NHC_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
       }/*endfor*/
     }else{
       iflag_mass = 2;
       for(ip=1;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
         apply_NHC_bead_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
#endif

         apply_NHC_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
       }/*endfor*/
      }/*endif : rank==0*/
   }/*endif respa*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_npti_pimd(CLASS *class,BONDED *bonded,
                             GENERAL_DATA *general_data,
                             int ir_tra,int ir_tor,int ir_ter,int ir_pimd,
                             double dt) 

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

   int ipart,iflag,ip;
   int ix_now;
   int iii,iproc;
   int iflag_mass;
   double dt2;
   double aa,aa2,arg2,poly,bb,dlen;
   double e2,e4,e6,e8;

   double *class_clatoms_mass;
   double *class_clatoms_x;
   double *class_clatoms_y;
   double *class_clatoms_z;
   double *class_clatoms_vx;
   double *class_clatoms_vy;
   double *class_clatoms_vz;
   double *class_clatoms_fx;
   double *class_clatoms_fy;
   double *class_clatoms_fz;
   double *class_clatoms_fxm;
   double *class_clatoms_fym;
   double *class_clatoms_fzm;

   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int nres_tra         = (general_data->timeinfo.nres_tra);
   int nres_tor         = (general_data->timeinfo.nres_tor);
   int nres_ter         = (general_data->timeinfo.nres_ter);
   int nres_pimd        = (general_data->timeinfo.nres_pimd);
   int ix_respa         = general_data->timeinfo.ix_respa;
   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;
   e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
   e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* 1) Evolve velocities: v_lnv_g = guess to v_lnv(dt)                       */

   for(ip=1;ip<=pi_beads_proc;ip++){
     class_clatoms_vx   = class->clatoms_pos[ip].vx;
     class_clatoms_vy   = class->clatoms_pos[ip].vy;
     class_clatoms_vz   = class->clatoms_pos[ip].vz;
     class_clatoms_fx   = class->clatoms_pos[ip].fx;
     class_clatoms_fy   = class->clatoms_pos[ip].fy;
     class_clatoms_fz   = class->clatoms_pos[ip].fz;
     class_clatoms_mass = class->clatoms_pos[ip].mass;
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                 /class_clatoms_mass[ipart];
     }/*endfor*/
   }/*endfor*/

/*==========================================================================*/
/* 2) Evolve NHCs and velocities                                            */

    ix_now = 5;
    if((ir_pimd==nres_pimd)){ix_now=4;}
    if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)){ix_now=3;}
    if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)&&(ir_tor==nres_tor))
                                   {ix_now=2;}
    if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)&&(ir_tor==nres_tor)
                                 &&(ir_ter==nres_ter)){ix_now=1;}
 
    if(ix_respa>=ix_now){
      if(pi_beads_proc_st == 1){
        iflag_mass = 1;
        apply_NHCPI_par(&(class->clatoms_info),
                        &(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),&(class->int_scr),
                        &(class->class_comm_forc_pkg));
        iflag_mass = 2;
        for(ip=2;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
          apply_NHC_bead_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
#endif
          apply_NHC_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
        }/*endfor*/
      }else{
        iflag_mass = 2;
        for(ip=1;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
          apply_NHC_bead_par(&(class->clatoms_info),
                             &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                             &(class->int_scr),iflag_mass,
                             &(class->class_comm_forc_pkg));
#endif

          apply_NHC_par(&(class->clatoms_info),
                             &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                             &(class->int_scr),iflag_mass,
                             &(class->class_comm_forc_pkg));
         }/*endfor*/
      }/*endif : rank*/
    }else{
      if(pi_beads_proc_st == 1){
         apply_NHCPI0_par(&(class->clatoms_info),
                          &(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(general_data->baro),&(class->int_scr),
                         &(class->class_comm_forc_pkg));
      }/*endif*/
    }/*endif respa*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nptf_pimd(CLASS *class,BONDED *bonded,
                             GENERAL_DATA *general_data,
                             int ir_tra,int ir_tor,int ir_ter,int ir_pimd,
                             double dt) 

/*========================================================================*/
      {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
   int ipart,iflag,ip;
   int ix_now;
   int iii,iproc;
   int iflag_mass;
   double dt2;

   double *class_clatoms_mass;
   double *class_clatoms_x;
   double *class_clatoms_y;
   double *class_clatoms_z;
   double *class_clatoms_vx;
   double *class_clatoms_vy;
   double *class_clatoms_vz;
   double *class_clatoms_fx;
   double *class_clatoms_fy;
   double *class_clatoms_fz;
   double *class_clatoms_fxm;
   double *class_clatoms_fym;
   double *class_clatoms_fzm;

   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int nres_tra         = (general_data->timeinfo.nres_tra);
   int nres_tor         = (general_data->timeinfo.nres_tor);
   int nres_ter         = (general_data->timeinfo.nres_ter);
   int nres_pimd        = (general_data->timeinfo.nres_pimd);
   int ix_respa         = general_data->timeinfo.ix_respa;
   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;
 
/*==========================================================================*/
/* 1) Evolve velocities                                                     */

   for(ip=1;ip<=pi_beads_proc;ip++){
     class_clatoms_vx   = class->clatoms_pos[ip].vx;
     class_clatoms_vy   = class->clatoms_pos[ip].vy;
     class_clatoms_vz   = class->clatoms_pos[ip].vz;
     class_clatoms_fx   = class->clatoms_pos[ip].fx;
     class_clatoms_fy   = class->clatoms_pos[ip].fy;
     class_clatoms_fz   = class->clatoms_pos[ip].fz;
     class_clatoms_mass = class->clatoms_pos[ip].mass;
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                 /class_clatoms_mass[ipart];
     }/*endfor*/
   }/*endfor*/
  
/*==========================================================================*/
/* 2) Evolve NHCs and velocities                                            */

   ix_now = 5;
   if((ir_pimd==nres_pimd)){ix_now=4;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)){ix_now=3;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)
                                   &&(ir_tor==nres_tor)){ix_now=2;}
   if((ir_pimd==nres_pimd)&&(ir_tra==nres_tra)&&(ir_tor==nres_tor)
                                 &&(ir_ter==nres_ter)){ix_now=1;}

   if(ix_respa>=ix_now){
     if(pi_beads_proc_st == 1){
       iflag_mass = 1;
       apply_NHCPF_par(&(class->clatoms_info),
                       &(class->clatoms_pos[1]),
                       &(class->therm_info_class),&(class->therm_class),
                       &(general_data->baro),
                       &(general_data->par_rahman),&(general_data->cell),
                       &(class->int_scr),&(class->class_comm_forc_pkg));
       iflag_mass = 2;
       for(ip=2;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
         apply_NHC_bead_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
#endif
         apply_NHC_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
       }/*endfor*/
     }else{
       iflag_mass = 2;
       for(ip=1;ip<=pi_beads_proc;ip++){

#ifdef FIX_OR_NUKE_ME
         apply_NHC_bead_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
#endif
         apply_NHC_par(&(class->clatoms_info),
                            &(class->clatoms_pos[ip]),
                            &(class->therm_info_bead),&(class->therm_bead[ip]),
                            &(class->int_scr),iflag_mass,
                            &(class->class_comm_forc_pkg));
       }/*endfor*/
     }/*endif : rank*/
   }else{
     if(pi_beads_proc_st == 1){
       apply_NHCPF0_par(&(class->clatoms_info),
                        &(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),
                        &(general_data->par_rahman),&(general_data->cell),
                        &(class->int_scr),&(class->class_comm_forc_pkg));
     }/*endif : rank==0*/
   }/*endif*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void anneal_bead(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                 INT_SCR *int_scr,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                 THERM_INFO *therm_info,THERM_POS *therm_pos,
                 double ann_rate,int iflag,int iflag_mass,int ip,
                 double *kinet_tmp)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/  

  int ipart,ichain,inhc;
  int natm_tot              = clatoms_info->natm_tot;
  int len_nhc               = therm_info->len_nhc;
  int num_nhc               = therm_info->num_nhc;
  int num_nhc_proc          = therm_info->num_nhc_proc;
  int mytherm_start         = therm_info->mytherm_start;
  int mytherm_end           = therm_info->mytherm_end;
  double v2;
  double *class_clatoms_vx  = clatoms_pos->vx;
  double *class_clatoms_vy  = clatoms_pos->vy;
  double *class_clatoms_vz  = clatoms_pos->vz;
  double *class_mass        = clatoms_pos->mass;
  double **therm_v_nhc      = therm_pos->v_nhc;
  double **therm_f_nhc      = therm_pos->f_nhc;
  double **therm_gkt        = therm_info->gkt;
  double **therm_mass_nhc   = therm_info->mass_nhc;
  int myatm_start           = clatoms_info->myatm_start; 
  int myatm_end             = clatoms_info->myatm_end; 

/*========================================================================*/
/* I) Particle annealing                                                  */


  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    class_clatoms_vx[ipart] *= ann_rate;
    class_clatoms_vy[ipart] *= ann_rate;
    class_clatoms_vz[ipart] *= ann_rate;
  }/*endfor*/

/*========================================================================*/
/* II) Nose-Hoover chain annaling                                         */


  if(iflag>=0){

    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        therm_v_nhc[ichain][inhc]    *= ann_rate;
        therm_f_nhc[ichain][inhc]    *= ann_rate*ann_rate;
        therm_gkt[ichain][inhc]      *= (ann_rate*ann_rate);
      }/*endfor*/
    }/*endfor*/

#ifdef JUNK
    init_NHC_par(clatoms_info,clatoms_pos,therm_info,therm_pos, 
                 int_scr,iflag_mass,class_comm_forc_pkg);
#endif

  }/* endif */

/*========================================================================*/
/* III) Check to see if target temperature has been reached               */
/*      If so, set an exit flag to be passed out                          */

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    v2 = class_clatoms_vx[ipart]*class_clatoms_vx[ipart]
       + class_clatoms_vy[ipart]*class_clatoms_vy[ipart]
       + class_clatoms_vz[ipart]*class_clatoms_vz[ipart];
    *kinet_tmp += class_mass[ipart]*v2;
  }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/




