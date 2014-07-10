/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVT                                      */
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
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nve(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
      {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,ifirst=1,iproc;
  double dt2,tol_glob;
  int iii;

  double *xold = class->clatoms_info.xold;
  double *yold = class->clatoms_info.yold;
  double *zold = class->clatoms_info.zold;
  double *x    = class->clatoms_pos[1].x;
  double *y    = class->clatoms_pos[1].y;
  double *z    = class->clatoms_pos[1].z;
  double *vx   = class->clatoms_pos[1].vx;
  double *vy   = class->clatoms_pos[1].vy;
  double *vz   = class->clatoms_pos[1].vz;
  double *fx   = class->clatoms_pos[1].fx;
  double *fy   = class->clatoms_pos[1].fy;
  double *fz   = class->clatoms_pos[1].fz;
  double *mass = class->clatoms_info.mass;

  int int_res_tra    = general_data->timeinfo.int_res_tra;
  int int_res_ter    = general_data->timeinfo.int_res_ter;
  int iconstrnt      = bonded->constrnt.iconstrnt;
  int np_forc        = class->communicate.np_forc;
  int myid_forc      = class->communicate.myid_forc;
  int myatm_start    = class->clatoms_info.myatm_start;
  int myatm_end      = class->clatoms_info.myatm_end;
  MPI_Comm comm_forc = class->communicate.comm_forc;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;

/*==========================================================================*/
/* VI) Evolve velocities                                                    */

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     vx[ipart] += fx[ipart]*dt2/mass[ipart];
     vy[ipart] += fy[ipart]*dt2/mass[ipart];
     vz[ipart] += fz[ipart]*dt2/mass[ipart];
   }/*endfor*/

/*==========================================================================*/
/* VII) Rattle if necessary (roll not)                                      */

   if( (int_res_tra==0) && (int_res_ter==0) ){
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 0;
       rattle_control(bonded,
                      &(class->clatoms_info),&(class->clatoms_pos[1]), 
                      &(general_data->cell),&(general_data->ptens),
                      &(general_data->statepoint),
                      &(general_data->baro),&(general_data->par_rahman),
                      &(general_data->stat_avg),dt,&tol_glob,ifirst,
                      &(class->class_comm_forc_pkg),&(class->ewd_scr));
     }/*endif*/
   }else{
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 0;
       rattle_control(bonded,
                      &(class->clatoms_info),&(class->clatoms_pos[1]), 
                      &(general_data->cell),&(general_data->ptens),
                      &(general_data->statepoint),
                      &(general_data->baro),&(general_data->par_rahman),
                      &(general_data->stat_avg),dt,&tol_glob,ifirst,
                      &(class->class_comm_forc_pkg),&(class->ewd_scr));
         get_pvten_inc_t(&(general_data->ptens),&(general_data->timeinfo),
                         ir_tra,ir_tor,ir_ter);
     }/*endif*/
  }/*endelse*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nvt(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,ipart,ifirst=1;
   int iflag_mass = 1;
   double dt2,tol_glob;
   int iii;
   int ix_now;

   double *class_clatoms_xold = class->clatoms_info.xold;
   double *class_clatoms_yold = class->clatoms_info.yold;
   double *class_clatoms_zold = class->clatoms_info.zold;
   double *class_clatoms_x    = class->clatoms_pos[1].x;
   double *class_clatoms_y    = class->clatoms_pos[1].y;
   double *class_clatoms_z    = class->clatoms_pos[1].z;
   double *class_clatoms_vx   = class->clatoms_pos[1].vx;
   double *class_clatoms_vy   = class->clatoms_pos[1].vy;
   double *class_clatoms_vz   = class->clatoms_pos[1].vz;
   double *class_clatoms_fx   = class->clatoms_pos[1].fx;
   double *class_clatoms_fy   = class->clatoms_pos[1].fy;
   double *class_clatoms_fz   = class->clatoms_pos[1].fz;
   double *class_clatoms_mass = class->clatoms_info.mass;

   int ix_respa               = general_data->timeinfo.ix_respa;
   int int_res_tra            = general_data->timeinfo.int_res_tra;
   int int_res_ter            = general_data->timeinfo.int_res_ter;
   int iconstrnt              = bonded->constrnt.iconstrnt;
   int myatm_start            = class->clatoms_info.myatm_start;
   int myatm_end              = class->clatoms_info.myatm_end;

   int nres_tra               = general_data->timeinfo.nres_tra;
   int nres_tor               = general_data->timeinfo.nres_tor;
   int nres_ter               = general_data->timeinfo.nres_ter;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;

/*==========================================================================*/
/* VII) Evolve velocities                                                   */

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                 /class_clatoms_mass[ipart];
      class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                 /class_clatoms_mass[ipart];
      class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                 /class_clatoms_mass[ipart];
   }/*endfor*/

/*==========================================================================*/
/* VIII) Rattle if necessary (roll not)                                     */

   if( (int_res_tra==0) && (int_res_ter==0) ){
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 0;
        rattle_control(bonded,
                       &(class->clatoms_info),&(class->clatoms_pos[1]), 
                       &(general_data->cell),&(general_data->ptens),
                       &(general_data->statepoint),
                       &(general_data->baro),&(general_data->par_rahman),
                       &(general_data->stat_avg),dt,&tol_glob,ifirst,
                       &(class->class_comm_forc_pkg),&(class->ewd_scr));
     }/*endif*/
   }else{
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 0;
       rattle_control(bonded,
                      &(class->clatoms_info),&(class->clatoms_pos[1]), 
                      &(general_data->cell),&(general_data->ptens),
                      &(general_data->statepoint),
                      &(general_data->baro),&(general_data->par_rahman),
                      &(general_data->stat_avg),dt,&tol_glob,ifirst,
                      &(class->class_comm_forc_pkg),&(class->ewd_scr));
       get_pvten_inc_t(&(general_data->ptens),&(general_data->timeinfo),
                       ir_tra,ir_tor,ir_ter);
     }/*endif*/  
   }/*endelse*/

/*==========================================================================*/
/* II) Evolve NHCs and velocities                                           */

   if( (int_res_tra==0) && (int_res_ter==0) ){
     if(class->therm_info_class.therm_typ == 1) {
       apply_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
     if(class->therm_info_class.therm_typ == 2 
	&& class->therm_info_class.len_nhc == 2) {
       apply_GGMT2_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
     if(class->therm_info_class.therm_typ == 2 
	&& class->therm_info_class.len_nhc == 3) {
       apply_GGMT3_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
   }else{
     ix_now = 4;
     if(ir_tra==nres_tra) {ix_now=3;}
     if((ir_tra==nres_tra) && (ir_tor==nres_tor)) {ix_now=2;}
     if((ir_tra==nres_tra) && (ir_tor==nres_tor) && (ir_ter==nres_ter)) {
         ix_now=1;
     }
     if(ix_respa>=ix_now){
      if(class->therm_info_class.therm_typ == 1) {
       apply_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(class->therm_info_class),&(class->therm_class),
                   &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
      }/*endif*/
      if(class->therm_info_class.therm_typ == 2 
	&& class->therm_info_class.len_nhc == 2) {
        apply_GGMT2_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
      }/*endif*/
      if(class->therm_info_class.therm_typ == 2 
	&& class->therm_info_class.len_nhc == 3) {
       apply_GGMT3_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
      }/*endif*/
     }/*endif*/

   }/*endelse*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nvt_isok(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                            int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
   int *inhc_x	   	          = class->therm_info_class.inhc_x;
   int *inhc_y		            = class->therm_info_class.inhc_y;
   int *inhc_z		            = class->therm_info_class.inhc_z;
   int i,ipart,ichain;
   int ifirst                 = 1;
   int iflag_mass             = 1;
   int tol_glob               = 1;
   int iii;
   int ix_now;
   double arg,s,sdot,a,b,rb;
   double *clatoms_vx         = class->clatoms_pos[1].vx;
   double *clatoms_vy         = class->clatoms_pos[1].vy;
   double *clatoms_vz         = class->clatoms_pos[1].vz;
   double *clatoms_fx         = class->clatoms_pos[1].fx;
   double *clatoms_fy         = class->clatoms_pos[1].fy;
   double *clatoms_fz         = class->clatoms_pos[1].fz;
   double *clatoms_mass       = class->clatoms_info.mass;
   double **therm_v1          = class->therm_class.x_nhc;
   double **therm_gkt	      = class->therm_info_class.gkt;
   int ix_respa               = general_data->timeinfo.ix_respa;
   int int_res_tra            = general_data->timeinfo.int_res_tra;
   int int_res_ter            = general_data->timeinfo.int_res_ter;
   int iconstrnt              = bonded->constrnt.iconstrnt;
   int myatm_start            = class->clatoms_info.myatm_start;
   int myatm_end              = class->clatoms_info.myatm_end;
   int mytherm_start          = class->therm_info_class.mytherm_start;
   int mytherm_end            = class->therm_info_class.mytherm_end;
   int nres_tra               = general_data->timeinfo.nres_tra;
   int nres_tor               = general_data->timeinfo.nres_tor;
   int nres_ter               = general_data->timeinfo.nres_ter;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   int len_nhc            = class->therm_info_class.len_nhc;    
   double lennhc          = (double)len_nhc;      
   double LkT             = lennhc*therm_gkt[1][inhc_x[1]];
   double kinetic         = 0.0;
   double dt2             = 0.5*dt;

   /*==========================================================================*/
/* VI) Isokinetic Constraint, Force-dependent portion                      */
	  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
		  a=clatoms_fx[ipart]*clatoms_vx[ipart]/LkT;
		  b=clatoms_fx[ipart]*clatoms_fx[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vx[ipart]=(clatoms_vx[ipart] + (clatoms_fx[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_x[ipart]]/=sdot;
		  }

		  a=clatoms_fy[ipart]*clatoms_vy[ipart]/LkT;
		  b=clatoms_fy[ipart]*clatoms_fy[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vy[ipart]=(clatoms_vy[ipart] + (clatoms_fy[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_y[ipart]]/=sdot;
		  }

		  a=clatoms_fz[ipart]*clatoms_vz[ipart]/LkT;
		  b=clatoms_fz[ipart]*clatoms_fz[ipart]/(clatoms_mass[ipart]*LkT);
		  rb=sqrt(b);
		  arg=dt2*rb;
                  if(arg > 0.00001) {
                     s = (1.0/rb)*sinh(arg) + (a/b)*(cosh(arg)-1.0);
                     sdot = cosh(arg) + (a/rb)*sinh(arg);
                  } else {
                     s = ((((b*a/24.0)*dt2 + b/6.0)*dt2 + 0.5*a)*dt2 + 1.0)*dt2;
                     sdot = (((b*a/6.0)*dt2 + 0.5*b)*dt2 + a)*dt2 + 1.0;
                  }
		  clatoms_vz[ipart]=(clatoms_vz[ipart] + (clatoms_fz[ipart]*s/clatoms_mass[ipart]))/sdot;
		  for(ichain=1;ichain<=len_nhc;ichain++){
			  therm_v1[ichain][inhc_z[ipart]]/=sdot;
		  }
	  }
/*==========================================================================*/

/*==========================================================================*/

/* VII) Isokinetic Constraint, Nose-like portion                             */

   if( (int_res_tra==0) && (int_res_ter==0) ){
       apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
      }else{
       ix_now = 4;
       if(ir_tra==1) {ix_now=3;}
       if((ir_tra==1) && (ir_tor==1)) {ix_now=2;}
       if((ir_tra==1) && (ir_tor==1) && (ir_ter==1)) {ix_now=1;}
       if(ix_respa>=ix_now){
           apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                              &(class->therm_info_class),&(class->therm_class),
                              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
      }
   }

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_npti(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                        int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,iflag,ifirst,ix_now;
  double dt2,tol_glob;
  double temp,temp_now;
  int iii;
  double aa,aa2,arg2,bb,poly;
  double e2,e4,e6,e8; 

  double *class_clatoms_vx   = class->clatoms_pos[1].vx;
  double *class_clatoms_vy   = class->clatoms_pos[1].vy;
  double *class_clatoms_vz   = class->clatoms_pos[1].vz;
  double *class_clatoms_fx   = class->clatoms_pos[1].fx;
  double *class_clatoms_fy   = class->clatoms_pos[1].fy;
  double *class_clatoms_fz   = class->clatoms_pos[1].fz;
  double *class_clatoms_mass = class->clatoms_info.mass;

  int int_res_tra            = general_data->timeinfo.int_res_tra;
  int int_res_ter            = general_data->timeinfo.int_res_ter;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  double tolratl             = bonded->constrnt.tolratl;
  int ix_respa               = general_data->timeinfo.ix_respa;
  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  double mass_lnv            = general_data->baro.mass_lnv;

  int nres_tra               = general_data->timeinfo.nres_tra;
  int nres_tor               = general_data->timeinfo.nres_tor;
  int nres_ter               = general_data->timeinfo.nres_ter;

        /* Careful because these values change */
  double v_lnv_g;
  double v_lnv               = general_data->baro.v_lnv;
  double f_lnv_v             = general_data->baro.f_lnv_v;
  double f_lnv_p             = general_data->baro.f_lnv_p;
  double roll_scg            = general_data->baro.roll_scg;
  double roll_scg0           = general_data->baro.roll_scg0;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);
  dt2 = dt/2.0;

/*==========================================================================*/
/* V.V) Approx v_lnv(dt)                                                    */

  // for parallel running, the following might be correct. To check!!
  // v_lnv_g = v_lnv*roll_scg0 + (f_lnv_p+f_lnv_v)*dt2*roll_scg/mass_lnv;
  v_lnv_g = v_lnv*roll_scg + (f_lnv_p+f_lnv_v)*dt2*roll_scg/mass_lnv;
  general_data->baro.v_lnv_g  = v_lnv_g;

/*==========================================================================*/
/*  RATTLE/ROLL CONVERGENCE LOOP                                            */ 

  if(iconstrnt==1){
    iflag = 1;
    cpysys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
               &(class->therm_info_class),&(class->therm_class),
               &(general_data->baro),
               &(general_data->par_rahman),&(class->int_scr),iflag);
  }/*endif*/
  tol_glob = tolratl+1.0;
  ifirst = 1;
  while(tol_glob>=tolratl){
     if(iconstrnt==1){
       getsys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                  &(class->therm_info_class),&(class->therm_class),
                  &(general_data->baro),
                  &(general_data->par_rahman),&(class->int_scr),iflag);
     }/*endif*/

/*==========================================================================*/
/* VII) Evolve velocities:                                                  */
       
     v_lnv = general_data->baro.v_lnv;
     aa   = exp(-0.5*dt2*v_lnv*(general_data->baro.c2_lnv));
     aa2  = aa*aa;
     arg2 = (0.5*dt2*v_lnv*(general_data->baro.c2_lnv))
          * (0.5*dt2*v_lnv*(general_data->baro.c2_lnv));
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     bb   = poly*aa;
     general_data->baro.roll_scvv = aa2;
     general_data->baro.roll_scf = bb;

     for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
        class_clatoms_vx[ipart] = class_clatoms_vx[ipart]*aa2
                                + bb*class_clatoms_fx[ipart]*dt2
                                / class_clatoms_mass[ipart];
        class_clatoms_vy[ipart] = class_clatoms_vy[ipart]*aa2
                                + bb*class_clatoms_fy[ipart]*dt2
                                / class_clatoms_mass[ipart];
       class_clatoms_vz[ipart] = class_clatoms_vz[ipart]*aa2
                                + bb*class_clatoms_fz[ipart]*dt2
                                / class_clatoms_mass[ipart];
     }/*endfor*/
     general_data->baro.roll_scvv = aa2;
     general_data->baro.roll_scf = bb;
    
/*==========================================================================*/
/* VIII) Rattle/Roll if necessary :                                         */
       
     tol_glob = 0.0;
     if(iconstrnt==1){
        (bonded->constrnt.iroll) = 1;
        rattle_control(bonded,
                       &(class->clatoms_info),&(class->clatoms_pos[1]), 
                       &(general_data->cell),&(general_data->ptens),
                       &(general_data->statepoint),
                       &(general_data->baro),&(general_data->par_rahman),
                       &(general_data->stat_avg),dt,&tol_glob,ifirst,
                       &(class->class_comm_forc_pkg),&(class->ewd_scr));
     }/*endif*/  
     ifirst = 0;

/*==========================================================================*/
/* VII) Evolve velocities of "log volumn" (epsilon)                         */

      temp = 0.0;
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        temp += class_clatoms_mass[ipart]
             *  class_clatoms_vx[ipart]*class_clatoms_vx[ipart];
        temp += class_clatoms_mass[ipart]
             *  class_clatoms_vy[ipart]*class_clatoms_vy[ipart];
        temp += class_clatoms_mass[ipart]
             *  class_clatoms_vz[ipart]*class_clatoms_vz[ipart];
      }/*endfor*/
      (general_data->baro.f_lnv_v)  = (general_data->baro.c2_lnv*temp);
      (general_data->baro.v_lnv) += 
                (general_data->baro.f_lnv_v+general_data->baro.f_lnv_p)
                                 *  dt2/(general_data->baro.mass_lnv);
     
/*==========================================================================*/
/* IX) Evolve NHCS and velocities                                           */

     if( (int_res_ter==0) && (int_res_tra==0) ){
        apply_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),
                        &(class->int_scr),&(class->class_comm_forc_pkg));
     }else{
        ix_now = 4;
        if(ir_tra==nres_tra) {ix_now=3;}
        if((ir_tra==nres_tra) && (ir_tor==nres_tor)) {ix_now=2;}
        if((ir_tra==nres_tra) && (ir_tor==nres_tor) && (ir_ter==nres_ter)) {
            ix_now=1;
        }
        if(ix_respa>=ix_now){
          apply_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(general_data->baro),
                          &(class->int_scr),&(class->class_comm_forc_pkg));
        }else{
          printf("respa for NPT: not ready for this part!\n");
          exit(1);
        }/*endif*/
     }/*endelse*/

  }/*endwhile:rattle roll*/

/*========================================================================*/
/* Respa pressure tensor update */

  if( (int_res_ter!=0) || (int_res_tra!=0) ){

    if(iconstrnt==1){
      get_pvten_inc_t(&(general_data->ptens),&(general_data->timeinfo),
                      ir_tra,ir_tor,ir_ter);
    }/*endif*/

  }/*endif*/

/*--------------------------------------------------------------------------*/
     }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_dt2_to_dt_nptf(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                        int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,j,joff,n,ipart,iflag,ifirst;
   int iii;
   double dt2,tol_glob;
   double aa,arg2,poly;
   double tempx,tempy,tempz;
   double tempvx,tempvy,tempvz;
   double tempfx,tempfy,tempfz;
   int ix_now;

   double *class_clatoms_vx   = class->clatoms_pos[1].vx;
   double *class_clatoms_vy   = class->clatoms_pos[1].vy;
   double *class_clatoms_vz   = class->clatoms_pos[1].vz;
   double *class_clatoms_fx   = class->clatoms_pos[1].fx;
   double *class_clatoms_fy   = class->clatoms_pos[1].fy;
   double *class_clatoms_fz   = class->clatoms_pos[1].fz;
   double *class_clatoms_mass = class->clatoms_info.mass;

   int myatm_start            = class->clatoms_info.myatm_start;
   int myatm_end              = class->clatoms_info.myatm_end;
   int hmat_int_typ           = general_data->cell.hmat_int_typ;

   double *vgmat_g            = general_data->par_rahman.vgmat_g;
   double *vgmat              = general_data->par_rahman.vgmat;
   double *hmat               = general_data->cell.hmat;
   double *fgmat_p            = general_data->par_rahman.fgmat_p;
   double *fgmat_v            = general_data->par_rahman.fgmat_v;
   double mass_hm             = general_data->par_rahman.mass_hm;
   int int_res_tra            = general_data->timeinfo.int_res_tra;
   int int_res_ter            = general_data->timeinfo.int_res_ter;
   int iconstrnt              = bonded->constrnt.iconstrnt;
   double tolratl             = bonded->constrnt.tolratl;
   int ix_respa               = general_data->timeinfo.ix_respa;

   int nres_tra               = general_data->timeinfo.nres_tra;
   int nres_tor               = general_data->timeinfo.nres_tor;
   int nres_ter               = general_data->timeinfo.nres_ter;

   double *roll_mtvv           = general_data->par_rahman.roll_mtvv;
   double *roll_mtf            = general_data->par_rahman.roll_mtf;
      /* careful because values are changed */
   double roll_scg            = general_data->par_rahman.roll_scg;
   double roll_scg0           = general_data->par_rahman.roll_scg0;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;

/*==========================================================================*/
/* VII) vgmat_g = guess to vgmat(dt)                     */

    for(i=1;i<=9;i++){
      vgmat_g[i] = vgmat[i]*roll_scg
                 + (fgmat_p[i]+fgmat_v[i])*dt2*roll_scg/mass_hm;
    }/*endfor*/

//  parallel running: the following might be important. To check!
//   for(i=1;i<=9;i++){
//      vgmat_g[i] = vgmat[i]*roll_scg0  
//                 + (fgmat_p[i]+fgmat_v[i])*dt2*roll_scg/mass_hm;
//   }/*endfor*/

/*==========================================================================*/
/*  RATTLE/ROLL CONVERGENCE LOOP                                            */ 

   if(iconstrnt==1){
       iflag = 2;
       cpysys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                  &(class->therm_info_class),&(class->therm_class),
                  &(general_data->baro),
                  &(general_data->par_rahman),&(class->int_scr),iflag);
   }/*endfor*/
   tol_glob = tolratl+1.0;
   ifirst = 1;
   while(tol_glob>=tolratl){
      if(iconstrnt==1){
         getsys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                    &(class->therm_info_class),&(class->therm_class),
                    &(general_data->baro),
                    &(general_data->par_rahman),&(class->int_scr),iflag);
      }/*endif*/      

/*==========================================================================*/
/* VII) Evolve velocities                                                   */

   if(hmat_int_typ==0){
      move_vel_vbox(&(class->clatoms_pos[1]),&(class->clatoms_info),
                    &(general_data->cell),&(general_data->par_rahman),dt);
   }else{
      move_vel_vbox_upper(&(class->clatoms_pos[1]),&(class->clatoms_info),
                    &(general_data->cell),&(general_data->par_rahman),dt);
   }

/*==========================================================================*/
/* VIII) Rattle/Roll if necessary :                                         */
      
      tol_glob = 0.0;
      if(iconstrnt==1){
        (bonded->constrnt.iroll) = 2;
        rattle_control(bonded,&(class->clatoms_info),
                       &(class->clatoms_pos[1]),
                       &(general_data->cell),&(general_data->ptens),
                       &(general_data->statepoint),
                       &(general_data->baro),&(general_data->par_rahman),
                       &(general_data->stat_avg),dt,&tol_glob,ifirst,
                       &(class->class_comm_forc_pkg),&(class->ewd_scr));
      }/*endif*/
      ifirst = 0;
      
/*==========================================================================*/
/* VIII) Evolve h-matrix velocities                                          */

    get_fgmatv_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(general_data->par_rahman),&(general_data->cell),
                   &(class->class_comm_forc_pkg));

    for(i=1;i<=9;i++){
      vgmat[i] += (fgmat_v[i] + fgmat_p[i])
                * dt2/mass_hm;
    }/*endfor*/

/*==========================================================================*/
/* IX) Evolve NHCS and velocities                                           */

      if((int_res_ter==0) && (int_res_tra==0) ){
        apply_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),
                        &(general_data->par_rahman),&(general_data->cell),
                       &(class->int_scr),&(class->class_comm_forc_pkg));
      }else{
        ix_now = 4;
        if(ir_tra==nres_tra) {ix_now=3;}
        if((ir_tra==nres_tra) && (ir_tor==nres_tor)) {ix_now=2;}
        if((ir_tra==nres_tra) && (ir_tor==nres_tor) && (ir_ter==nres_ter)) {
            ix_now=1;
        }
        if(ix_respa>=ix_now){
          apply_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(general_data->baro),
                          &(general_data->par_rahman),&(general_data->cell),
                          &(class->int_scr),&(class->class_comm_forc_pkg));
        }else{
          printf("respa for NPT: this part is not ready!\n");
          exit(1);
        }/*endif*/
      }/*endelse*/

   }/*endwhile:rattle roll*/

/*=========================================================================*/
/* Respa pressure tensor update */

  if( (int_res_ter!=0) || (int_res_tra!=0) ){
    if(iconstrnt==1){
       get_pvten_inc_t(&(general_data->ptens),&(general_data->timeinfo),
                       ir_tra,ir_tor,ir_ter);
    }/*endif*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_final_class(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                     int iflag)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/  

  int i,iii;
  double wght;

  int np_forc            = class->communicate.np_forc;
  int myid_forc          = class->communicate.myid_forc;
  int myid               = class->communicate.myid;
  int iconstrnt          = bonded->constrnt.iconstrnt;
  int iget_pv_real_inter = class->energy_ctrl.iget_pv_real_inter;
  int int_res_ter        = general_data->timeinfo.int_res_ter;
  int int_res_tra        = general_data->timeinfo.int_res_tra;
  double *pvten_tot      = general_data->ptens.pvten_tot;
  double *pvten_inc_a    = general_data->ptens.pvten_inc_a;
  double *pvten_inc      = general_data->ptens.pvten_inc;
  double *vgmat          = general_data->par_rahman.vgmat;

  int iperd              = general_data->cell.iperd;
  int hmat_cons_typ      = general_data->cell.hmat_cons_typ;
  int hmat_int_typ       = general_data->cell.hmat_int_typ;
  double *hmat           = general_data->cell.hmat;
 
  wght = 1.0/(double)np_forc;

/*==========================================================================*/
/* 0) Allgather the volume                         */

  if(np_forc > 1){class_int_allgather(general_data,class);}

/*==========================================================================*/
/* I) Get NHC contribution to energy                                        */

  if(iflag >=0){
    if(general_data->ensopts.nvt_isok==1){
    	nhc_vol_potkin_isok(&(class->clatoms_info),&(class->clatoms_pos[1]),&(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),&(general_data->par_rahman),
                &(general_data->stat_avg),&(general_data->statepoint),
                iflag,myid_forc,general_data->timeinfo.itime);
    }else{
        nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                       &(general_data->baro),&(general_data->par_rahman),
                       &(general_data->stat_avg),&(general_data->statepoint),
                       iflag,myid_forc);
    }
  }/* endif */

/*==========================================================================*/
/* II) Get Kinetic energy and kinetic energy tensor                       */

  get_tvten(&(class->clatoms_info),&(class->clatoms_pos[1]),
             &(general_data->stat_avg),&(general_data->ptens),
             &(general_data->cell));


/*==========================================================================*/
/* III) Add constraint contribution to the total Pressure Tensor             */

  if( (iconstrnt==1)&& (iget_pv_real_inter==1) ){

    if( (int_res_ter!=0) ||(int_res_tra!=0) ){
      get_pvten_inc_a(&(general_data->ptens),&(general_data->timeinfo));
      for(i=1;i<=9;i++){pvten_tot[i] += pvten_inc_a[i]*wght;}
    }else{
      for(i=1;i<=9;i++){pvten_tot[i] += pvten_inc[i]*wght;}
    }/*endif*/ 

  }/*endif*/ 

/*==========================================================================*/
/* IV) Cell velocity matrix symmetry checker                                */

  if(iflag == 2){
    chck_constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,vgmat,hmat,myid);
  }/*endif */

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void anneal_class(CLASS *class,double ann_rate,int iflag,int iflag_mass,
                  double target_temp,int *exit_flag)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/  
#include "../typ_defs/typ_mask.h"

  int ipart,ichain,inhc;
  int natm_tot              = class->clatoms_info.natm_tot;
  int len_nhc               = class->therm_info_class.len_nhc;
  int num_nhc               = class->therm_info_class.num_nhc;
  int num_proc              = class->class_comm_forc_pkg.num_proc;
  int myid                  = class->class_comm_forc_pkg.myid;
  MPI_Comm forc_comm        = class->class_comm_forc_pkg.comm;
  double *class_clatoms_vx  = class->clatoms_pos[1].vx;
  double *class_clatoms_vy  = class->clatoms_pos[1].vy;
  double *class_clatoms_vz  = class->clatoms_pos[1].vz;
  double *class_mass        = class->clatoms_info.mass;
  double **therm_v_nhc      = class->therm_class.v_nhc;
  double **therm_f_nhc      = class->therm_class.f_nhc;
  double **therm_gkt        = class->therm_info_class.gkt;
  double **therm_mass_nhc   = class->therm_info_class.mass_nhc;
  double v2,kinet_tmp,kinet,temp;
  int myatm_start           = class->clatoms_info.myatm_start; 
  int myatm_end             = class->clatoms_info.myatm_end; 
  int nfree                 = class->clatoms_info.nfree;

/*========================================================================*/
/* I) Particle annealing                                                  */

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    class_clatoms_vx[ipart] *= ann_rate;
    class_clatoms_vy[ipart] *= ann_rate;
    class_clatoms_vz[ipart] *= ann_rate;
  }/*endfor*/

/*========================================================================*/
/* II) Nose-Hoover) chain ennaling                                        */

  if(iflag>=0){

    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=1;inhc<=num_nhc;inhc++){
        therm_v_nhc[ichain][inhc]    *= ann_rate;
        therm_f_nhc[ichain][inhc]    *= ann_rate*ann_rate;
        therm_gkt[ichain][inhc]      *= (ann_rate*ann_rate);
      }/*endfor*/
    }/*endfor*/

#ifdef JUNK
    init_NHC_par(&class->clatoms_info,&class->clatoms_pos[1], 
                 &class->therm_info_class,&class->therm_class, 
                 &class->int_scr,iflag_mass,&(class->class_comm_forc_pkg));
#endif

  }/* endif */

/*========================================================================*/
/* III) Check to see if target temperature has been reached               */
/*      If so, set an exit flag to be passed out                          */

  kinet_tmp = 0.0;
  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    v2 = class_clatoms_vx[ipart]*class_clatoms_vx[ipart]
       + class_clatoms_vy[ipart]*class_clatoms_vy[ipart]
       + class_clatoms_vz[ipart]*class_clatoms_vz[ipart];
    kinet_tmp += class_mass[ipart]*v2;
  }/*endfor*/

  if(num_proc > 1){
    Allreduce(&kinet_tmp,&kinet,1,MPI_DOUBLE,MPI_SUM,0,forc_comm);
  } else {    kinet = kinet_tmp;
  }/* endif */

  temp = kinet*BOLTZ/((double) (nfree));

  if(ann_rate > 1.0 && temp > target_temp) {
     *exit_flag = 1;
     if(myid == 0) {
       printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Target temperature reached in annealing run\n");
       printf("Performing last step and preparing to exit\n");
       printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     }/* endif */
  }/* endif */
  if(ann_rate < 1.0 && temp < target_temp) {
     *exit_flag = 1;
     if(myid == 0){
       printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("Target temperature reached in annealing run\n");
       printf("Performing to last step and preparing to exit\n");
       printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     }/* endif */
  }/* endif */
  

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void class_int_allgather(GENERAL_DATA *general_data,CLASS *class)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/

#include "../typ_defs/typ_mask.h"

  int iii,i,j;
  int nvt              = general_data->ensopts.nvt;
  int npti             = general_data->ensopts.npt_i;
  int nptf             = general_data->ensopts.npt_f;
  int natm_tot         = class->clatoms_info.natm_tot;
  int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
  int np_forc          = class->communicate.np_forc;
  int myid_forc        = class->communicate.myid_forc;
  int myid_bead_prime  = class->communicate.myid_bead_prime;
  MPI_Comm comm_forc   = class->communicate.comm_forc;
  MPI_Comm comm_states = class->communicate.comm_states;
  double **v_nhc = class->therm_class.v_nhc;
  double *int_scr_sc = class->int_scr.sc;  
  double *int_scr_sc_temp = class->int_scr.sc_temp;  
  int *map_share    = class->therm_info_class.map_share;
  int num_nhc_share = class->therm_info_class.num_nhc_share;
  int len_nhc       = class->therm_info_class.len_nhc;
  int *displs_therm = class->class_comm_forc_pkg.displs_therm;
  int *recv_count_therm = class->class_comm_forc_pkg.recv_count_therm;

/*==========================================================================*/
/* I) Bcast the thermostats */
  
  if(nvt + npti + nptf > 0){

     for(i=1;i<=len_nhc;i++){

      for(j=1;j<=num_nhc_share+1;j++){
        int_scr_sc[j] = 0.0;
      }/*endfor*/

      for(j=1;j<=num_nhc_share;j++){
        int_scr_sc[j] = v_nhc[i][map_share[j]];      
      }/*endfor*/

      if(num_nhc_share > 0){
       Allgatherv(&(int_scr_sc[displs_therm[myid_forc+1]+1]),
                   recv_count_therm[myid_forc+1],MPI_DOUBLE,
                   &(int_scr_sc_temp[1]),&recv_count_therm[1],&displs_therm[1],
                   MPI_DOUBLE,0,comm_forc);
      }/* endif */


      for(j=1;j<=num_nhc_share;j++){
         v_nhc[i][map_share[j]] = int_scr_sc_temp[j];      
      }/*endfor*/

     }/*endfor*/

  }/*endif*/

/*==========================================================================*/
/* I) Bcast the barostats/cell information */

  if(npti==1){
    Bcast(&(general_data->cell.hmat[0]),10,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.x_lnv),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.x_lnv_o),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.v_lnv),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.v_lnv_g),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.f_lnv_p),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.f_lnv_v),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.roll_scg0),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->baro.roll_scg),1,MPI_DOUBLE,0,comm_forc);
    gethinv((general_data->cell.hmat),(general_data->cell.hmati),
           &(general_data->cell.vol),general_data->cell.iperd);
    general_data->baro.vol        = general_data->cell.vol;
    general_data->par_rahman.vol  = general_data->cell.vol;
  }/*endif*/

  if(nptf==1){
    Bcast(&(general_data->cell.hmat[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->baro.hmato[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->par_rahman.vgmat[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->par_rahman.vgmat_g[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->par_rahman.fgmat_p[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->par_rahman.fgmat_v[0]),10,MPI_DOUBLE,0,comm_states);
    Bcast(&(general_data->par_rahman.roll_scg),1,MPI_DOUBLE,0,comm_forc);
    Bcast(&(general_data->par_rahman.roll_scg0),1,MPI_DOUBLE,0,comm_forc);
    gethinv((general_data->cell.hmat),(general_data->cell.hmati),
           &(general_data->cell.vol),general_data->cell.iperd);
    general_data->baro.vol        = general_data->cell.vol;
    general_data->par_rahman.vol  = general_data->cell.vol;
  }/*endif*/
 
/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/








