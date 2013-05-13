/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_0_to_dt2_pimd                            */
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
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_integrate_pimd_entry.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nvt_pimd(CLASS *class,BONDED *bonded,
                       GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,int ir_pimd,double dt) 

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int ipart,iflag,ip;
   int ix_now;
   int iii;
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

   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
   int ix_respa         = general_data->timeinfo.ix_respa;
   int nres_pimd        = general_data->timeinfo.nres_pimd;
   int myid_bead        = class->communicate.myid_bead;
   MPI_Comm comm_bead   = class->communicate.comm_beads;

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;

/*==========================================================================*/
/* 1) Evolve NHCs and velocities                                            */
          
   ix_now = 5;
   if(ir_pimd==1){ix_now=4;}
   if((ir_pimd==1)&&(ir_tra==1)){ix_now=3;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)){ix_now=2;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}

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
     }/*endif : pi_beads_st*/
   }/*endif respa*/
/*==========================================================================*/
/* 2) Evolve velocities                                                     */

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
/* 3) Evolve positions                                                      */

   if(ir_pimd==1){         
     control_pimd_trans_mode(class,general_data);
   }/*endif*/
   for(ip=1;ip<=pi_beads_proc;ip++){
     class_clatoms_x    = class->clatoms_pos[ip].x;
     class_clatoms_y    = class->clatoms_pos[ip].y;
     class_clatoms_z    = class->clatoms_pos[ip].z;
     class_clatoms_vx    = class->clatoms_pos[ip].vx;
     class_clatoms_vy    = class->clatoms_pos[ip].vy;
     class_clatoms_vz    = class->clatoms_pos[ip].vz;
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_x[ipart] += class_clatoms_vx[ipart]*dt;
       class_clatoms_y[ipart] += class_clatoms_vy[ipart]*dt;
       class_clatoms_z[ipart] += class_clatoms_vz[ipart]*dt;
     }/*endfor*/
   }/*endfor*/

   mode_energy_control(class,general_data);
         
   if(ir_pimd==nres_pimd){
     control_pimd_trans_pos(class,general_data);
   }/*endif*/

/*==========================================================================*/
/* 4) Recalculate positions of ghost atoms                                  */

   for(ip=1;ip<=pi_beads_proc;ip++){
     get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                   &(class->ghost_atoms));
   }/*endfor*/


/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_npti_pimd(CLASS *class,BONDED *bonded,
                       GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,int ir_pimd,
                       double dt) 

/*========================================================================*/
{/*begin routine*/
#include "../typ_defs/typ_mask.h"

/*========================================================================*/
/*             Local variable declarations                                */

   int i,ipart,ip,iflag,iflag_mass;
   int ix_now;
   int iii;
   double aa,aa2,arg2,poly,bb,dlen,vol;
   double e2,e4,e6,e8;
   double dt2,v_lnv,x_lnv,x_lnv_o;

   double *class_clatoms_x;   
   double *class_clatoms_y;   
   double *class_clatoms_z;   
   double *class_clatoms_vx;  
   double *class_clatoms_vy;  
   double *class_clatoms_vz;  
   double *class_clatoms_fx;  
   double *class_clatoms_fy;  
   double *class_clatoms_fz;  
   double *class_clatoms_mass;

   double *xmod = class->clatoms_info.xmod;
   double *ymod = class->clatoms_info.ymod;
   double *zmod = class->clatoms_info.zmod;

   double *hmat = general_data->cell.hmat;
   double *hmati = general_data->cell.hmati;
   int    iperd  = general_data->cell.iperd;

   int np_forc          = class->communicate.np_forc;
   int np_beads         = class->communicate.np_beads;
   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
   MPI_Comm comm_beads  = class->communicate.comm_beads;
   MPI_Comm comm_beads_forc   = class->communicate.comm_beads_forc;
   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int ix_respa         = general_data->timeinfo.ix_respa;
   int nres_pimd        = general_data->timeinfo.nres_pimd;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;

   e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
   e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* 1) Evolve NHCs and velocities                                            */

   ix_now = 5;
   if(ir_pimd==1){ix_now=4;}
   if((ir_pimd==1)&&(ir_tra==1)){ix_now=3;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)){ix_now=2;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}

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

#ifdef FIX_OR_NUK_ME
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

/*==========================================================================*/
/* 2) Evolve velocities                                                     */

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
/* 3) Evolve positions                                                      */
/*--------------------------------------------------------------------------*/
/*    Send v_lnv to all processors                                          */

   if(np_forc>1){
      Bcast(&(general_data->baro.v_lnv),1,MPI_DOUBLE,0,comm_beads_forc);
   }/*endif*/
   v_lnv = general_data->baro.v_lnv;

   if(ir_pimd==1){control_pimd_trans_mode(class,general_data);}

   if(pi_beads_proc_st==1){
       aa = exp(dt2*(general_data->baro.v_lnv));
       aa2 = aa*aa;
       arg2 = (v_lnv*dt2)*(v_lnv*dt2); 
       poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
       bb   = aa*poly;
   }/*endif*/
   for(ip=1;ip<=pi_beads_proc;ip++){
     class_clatoms_x    = class->clatoms_pos[ip].x;
     class_clatoms_y    = class->clatoms_pos[ip].y;
     class_clatoms_z    = class->clatoms_pos[ip].z;
     class_clatoms_vx   = class->clatoms_pos[ip].vx;
     class_clatoms_vy   = class->clatoms_pos[ip].vy;
     class_clatoms_vz   = class->clatoms_pos[ip].vz;
     if(pi_beads_proc_st==1&&ip==1){
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         class_clatoms_x[ipart] = class_clatoms_x[ipart]*aa2
                                + class_clatoms_vx[ipart]*bb*dt;
         class_clatoms_y[ipart] = class_clatoms_y[ipart]*aa2 
                                + class_clatoms_vy[ipart]*bb*dt;
         class_clatoms_z[ipart] = class_clatoms_z[ipart]*aa2
                                 + class_clatoms_vz[ipart]*bb*dt;
       }/*endfor*/
     }else{
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         class_clatoms_x[ipart] += class_clatoms_vx[ipart]*dt;
         class_clatoms_y[ipart] += class_clatoms_vy[ipart]*dt;
         class_clatoms_z[ipart] += class_clatoms_vz[ipart]*dt;
       }/*endfor*/
     }/*endif*/
   }/*endfor*/

   mode_energy_control(class,general_data);

   if(ir_pimd==nres_pimd){control_pimd_trans_pos(class,general_data);}
 
/*==========================================================================*/
/* 4) Evolve cell                                                           */

   if(np_beads>1){
     Bcast(&(general_data->baro.v_lnv),1,MPI_DOUBLE,0,comm_beads);
   }/*endif*/
   v_lnv   = general_data->baro.v_lnv;
   x_lnv   = general_data->baro.x_lnv;
   x_lnv_o = general_data->baro.x_lnv;

   x_lnv   = x_lnv_o + v_lnv*dt;
   dlen    = exp(x_lnv-x_lnv_o);
   for(i=1;i<=9;i++){hmat[i]*=dlen;}
   gethinv(hmat,hmati,&(vol),iperd);

   general_data->cell.vol     = vol;
   general_data->baro.vol     = vol;
   general_data->baro.x_lnv   = x_lnv;
   general_data->baro.x_lnv_o = x_lnv_o;

/*==========================================================================*/
/*  5) Recalculate positions of ghost atoms                                 */

   for(ip=1;ip<=pi_beads_proc;ip++){
      get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                    &(class->ghost_atoms));
   }/*endfor*/
 
/*--------------------------------------------------------------------------*/
      }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nptf_pimd(CLASS *class,BONDED *bonded,
                            GENERAL_DATA *general_data,
                            int ir_tra,int ir_tor,int ir_ter,int ir_pimd,
                            double dt) 

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   int i,j,ip,ip_start,joff,n=3,ipart,iflag;
   int ix_now,iflag_mass;
   int iii;
   double aa,arg2,poly,dt_neg;
   double tempx,tempy,tempz;
   double tempvx,tempvy,tempvz;
   double e2,e4,e6,e8;
   double dt2;

   double *class_clatoms_x;  
   double *class_clatoms_y;  
   double *class_clatoms_z;  
   double *class_clatoms_vx; 
   double *class_clatoms_vy; 
   double *class_clatoms_vz; 
   double *class_clatoms_fx; 
   double *class_clatoms_fy; 
   double *class_clatoms_fz; 
   double *class_clatoms_mass;

   double *hmat         = general_data->cell.hmat;
   double *hmato        = general_data->par_rahman.hmato;

   int pi_beads         = class->clatoms_info.pi_beads;
   int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
   int np_beads         = class->communicate.np_beads;
   MPI_Comm comm_beads  = class->communicate.comm_beads;

   int myatm_start      = class->clatoms_info.myatm_start;
   int myatm_end        = class->clatoms_info.myatm_end;
   int nres_pimd        = general_data->timeinfo.nres_pimd;
   int ix_respa         = general_data->timeinfo.ix_respa;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   ip_start = (pi_beads_proc_st == 1 ? 2:1);

   dt2 = dt/2.0;
   e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
   e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* 1)Save class cell                                                       */

   for(i=1;i<=9;i++){hmato[i]=hmat[i];}

/*==========================================================================*/
/* 2) Evolve NHCs and velocities                                            */
          
   ix_now = 5;
   if(ir_pimd==1){ix_now=4;}
   if((ir_pimd==1)&&(ir_tra==1)){ix_now=3;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)){ix_now=2;}
   if((ir_pimd==1)&&(ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}

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
   }/*endif : respa*/

/*==========================================================================*/
/* 3) Evolve velocities                                                     */

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
/* 4) Evolve the h-matrix and all positions except the first bead           */

   if(ir_pimd==1){control_pimd_trans_mode(class,general_data);}


   for(ip=ip_start;ip<=pi_beads_proc;ip++){
     class_clatoms_x    = class->clatoms_pos[ip].x;
     class_clatoms_y    = class->clatoms_pos[ip].y;
     class_clatoms_z    = class->clatoms_pos[ip].z;
     class_clatoms_vx   = class->clatoms_pos[ip].vx;
     class_clatoms_vy   = class->clatoms_pos[ip].vy;
     class_clatoms_vz   = class->clatoms_pos[ip].vz;
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_x[ipart] += class_clatoms_vx[ipart]*dt;
       class_clatoms_y[ipart] += class_clatoms_vy[ipart]*dt;
       class_clatoms_z[ipart] += class_clatoms_vz[ipart]*dt;
     }/*endfor*/
   }/*endfor*/

   if(pi_beads_proc_st == 1){
     class_clatoms_x    = class->clatoms_pos[1].x;
     class_clatoms_y    = class->clatoms_pos[1].y;
     class_clatoms_z    = class->clatoms_pos[1].z;
     class_clatoms_vx   = class->clatoms_pos[1].vx;
     class_clatoms_vy   = class->clatoms_pos[1].vy;
     class_clatoms_vz   = class->clatoms_pos[1].vz;
     move_pos_box_upper(class_clatoms_x,class_clatoms_y,class_clatoms_z,
                        class_clatoms_vx,class_clatoms_vy,class_clatoms_vz,
                        general_data->cell.hmat,general_data->cell.hmati, 
                        general_data->par_rahman.vgmat,dt,
                        class->clatoms_info.natm_tot,myatm_start,myatm_end);
   }/*endif*/

   mode_energy_control(class,general_data);

   if(ir_pimd==nres_pimd){control_pimd_trans_pos(class,general_data);}

/*==========================================================================*/
/* 5) Communicate the cell information                                      */

  if(np_beads>1&&ir_pimd==nres_pimd){
      Bcast(&(general_data->cell.hmat[1]),9,MPI_DOUBLE,0,
            class->communicate.comm_beads);
  }/*endif*/
   
  gethinv(general_data->cell.hmat,general_data->cell.hmati,
          &(general_data->cell.vol),general_data->cell.iperd);


  general_data->par_rahman.vol = general_data->cell.vol;
  general_data->baro.vol       = general_data->cell.vol;
      
/*==========================================================================*/
/* 5) Recalculate positions of ghost atoms                                  */

  if(ir_pimd==nres_pimd){         
     for(ip=1;ip<=pi_beads_proc;ip++){
       get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                     &(class->ghost_atoms));
     }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/







