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
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_pimd_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define JUNK_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nve(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                      int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,ifirst=1;
  double dt2,tol_glob=0.0;
  int iii;

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

  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  int myid_forc              = class->communicate.myid_forc;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  dt2 = dt/2.0;

/*==========================================================================*/
/* I)Save positions                                                         */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_xold[ipart] = class_clatoms_x[ipart];
    class_clatoms_yold[ipart] = class_clatoms_y[ipart];
    class_clatoms_zold[ipart] = class_clatoms_z[ipart];
  }/*endfor*/


/*==========================================================================*/
/* II) Evolve velocities                                                    */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                              /class_clatoms_mass[ipart];
    class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                              /class_clatoms_mass[ipart];
    class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                              /class_clatoms_mass[ipart];
  }/*endfor*/

/*==========================================================================*/
/* III) Evolve positions                                                    */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_x[ipart] += class_clatoms_vx[ipart]*dt;
    class_clatoms_y[ipart] += class_clatoms_vy[ipart]*dt;
    class_clatoms_z[ipart] += class_clatoms_vz[ipart]*dt;
  }/*endfor*/

/*==========================================================================*/
/* IV) Shake if necessary (roll not)                                        */

  if(iconstrnt==1){
    (bonded->constrnt.iroll) = 0;
    shake_control(bonded,
                 &(class->clatoms_info),&(class->clatoms_pos[1]), 
                 &(general_data->cell),&(general_data->ptens),
                 &(general_data->statepoint),
                 &(general_data->baro),&(general_data->par_rahman),
                 &(general_data->stat_avg),dt,&tol_glob,ifirst,
                 &(class->class_comm_forc_pkg),&(class->ewd_scr));
  }/*endif*/

/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

  get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->ghost_atoms));

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nvt(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                      int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,ifirst=1;
  int iflag_mass = 1;
  double dt2,tol_glob=0.0;
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

  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  int myid_forc              = class->communicate.myid_forc;
  int int_res_tra            = general_data->timeinfo.int_res_tra;
  int int_res_ter            = general_data->timeinfo.int_res_ter;
  int ix_respa               = general_data->timeinfo.ix_respa;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  dt2 = dt/2.0;

/*==========================================================================*/
/* I)Save positions                                                         */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_xold[ipart] = class_clatoms_x[ipart];
    class_clatoms_yold[ipart] = class_clatoms_y[ipart];
    class_clatoms_zold[ipart] = class_clatoms_z[ipart];
  }/*endfor*/

/*==========================================================================*/
/* II) Evolve NHCs and velocities                                           */

  if( (int_res_tra==0) && (int_res_ter==0) ){
    if(class->therm_info_class.therm_typ == 1){
      apply_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->therm_info_class),&(class->therm_class),
                 &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
    }/*endif*/
    if(class->therm_info_class.therm_typ == 2 
       && class->therm_info_class.len_nhc == 2){
      apply_GGMT2_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->therm_info_class),&(class->therm_class),
                 &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
    }/*endif*/
    if(class->therm_info_class.therm_typ == 2 
       && class->therm_info_class.len_nhc == 3){
      apply_GGMT3_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->therm_info_class),&(class->therm_class),
                 &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
    }/*endif*/
  }else{
    ix_now = 4;
    if((ir_tra==1)){ix_now=3;}
    if((ir_tra==1)&&(ir_tor==1)){ix_now=2;}
    if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}
    if(ix_respa>=ix_now){
     if(class->therm_info_class.therm_typ == 1){     
      apply_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(class->therm_info_class),&(class->therm_class),
                   &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
     if(class->therm_info_class.therm_typ == 2 
       && class->therm_info_class.len_nhc == 2){
       apply_GGMT2_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(class->therm_info_class),&(class->therm_class),
                   &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
     if(class->therm_info_class.therm_typ == 2 
       && class->therm_info_class.len_nhc == 3){
       apply_GGMT3_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(class->therm_info_class),&(class->therm_class),
                   &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
     }/*endif*/
    }/*endif*/

  }/*endelse*/

/*==========================================================================*/
/* III) Evolve velocities                                                   */

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                              /class_clatoms_mass[ipart];
    class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                              /class_clatoms_mass[ipart];
    class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                              /class_clatoms_mass[ipart];
  }/*endfor*/

/*==========================================================================*/
/* IV) Evolve positions                                                     */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_x[ipart] += class_clatoms_vx[ipart]*dt;
    class_clatoms_y[ipart] += class_clatoms_vy[ipart]*dt;
    class_clatoms_z[ipart] += class_clatoms_vz[ipart]*dt;
  }/*endfor*/

/*==========================================================================*/
/* V) Shake if necessary (roll not)                                         */

  if(iconstrnt==1){
    (bonded->constrnt.iroll) = 0;
    shake_control(bonded,
                  &(class->clatoms_info),&(class->clatoms_pos[1]), 
                  &(general_data->cell),&(general_data->ptens),
                  &(general_data->statepoint),
                  &(general_data->baro),&(general_data->par_rahman),
                  &(general_data->stat_avg),dt,&tol_glob,ifirst,
                  &(class->class_comm_forc_pkg),&(class->ewd_scr));
  }/*endif*/

/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

  get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->ghost_atoms));

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nvt_isok(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                           int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,ipart,ichain,ifirst=1;
  int iflag_mass = 1;
  double dt2,tol_glob=0.0;
  int iii;
  int ix_now;
  double *random;
  double cst1, cst2, cst3;
  double a,b,rb,s,sdot,arg;   /* Isokinetic Constraint Variables */
  int *inhc_x			           = class->therm_info_class.inhc_x;
  int *inhc_y			           = class->therm_info_class.inhc_y;
  int *inhc_z			           = class->therm_info_class.inhc_z;
  double gamma				       = 1.0/(200.0*dt);  /* Stochastic Term */
  int *iseed			           = &class->vel_samp_class.iseed;
  int *iseed2				         = &class->vel_samp_class.iseed2;
  double *qseed				       = &class->vel_samp_class.qseed;
  double *clatoms_xold       = class->clatoms_info.xold;
  double *clatoms_yold       = class->clatoms_info.yold;
  double *clatoms_zold       = class->clatoms_info.zold;
  double *clatoms_x          = class->clatoms_pos[1].x;
  double *clatoms_y          = class->clatoms_pos[1].y;
  double *clatoms_z          = class->clatoms_pos[1].z;
  double *clatoms_vx         = class->clatoms_pos[1].vx;
  double *clatoms_vy         = class->clatoms_pos[1].vy;
  double *clatoms_vz         = class->clatoms_pos[1].vz;
  double *clatoms_fx         = class->clatoms_pos[1].fx;
  double *clatoms_fy         = class->clatoms_pos[1].fy;
  double *clatoms_fz         = class->clatoms_pos[1].fz;
  double *clatoms_mass       = class->clatoms_info.mass;
  double **therm_mass_nhc 	 = class->therm_info_class.mass_nhc;
  double **therm_v1          = class->therm_class.x_nhc;
  double **therm_v2          = class->therm_class.v_nhc;
  double **therm_gkt		     = class->therm_info_class.gkt;
  double *int_scr_atm_kin    = class->int_scr.atm_kin;			  

  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  int mytherm_start          = class->therm_info_class.mytherm_start;
  int mytherm_end            = class->therm_info_class.mytherm_end;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  int myid_forc              = class->communicate.myid_forc;
  int int_res_tra            = general_data->timeinfo.int_res_tra;
  int int_res_ter            = general_data->timeinfo.int_res_ter;
  int ix_respa               = general_data->timeinfo.ix_respa;
  int len_nhc 				       = class->therm_info_class.len_nhc;	  
  double lennhc		  	       = (double)len_nhc;			
  double LkT                 = lennhc*therm_gkt[1][inhc_x[1]];
  double stoch_e             = exp(-gamma*dt);
  double mass_factor         = 1.0;
  double sigma 		     = sqrt(therm_gkt[1][inhc_x[1]]*(1.0-exp(-2.0*gamma*dt))
		  	  	  	  	    /(mass_factor*therm_mass_nhc[1][inhc_x[1]]));
  double kinetic             = 0.0;




/*==========================================================================*/
/* 0) Useful constants                                                      */


  dt2 = dt/2.0;

/*==========================================================================*/
/* I)Save positions                                                         */

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    clatoms_xold[ipart] = clatoms_x[ipart];
    clatoms_yold[ipart] = clatoms_y[ipart];
    clatoms_zold[ipart] = clatoms_z[ipart];
  }/*endfor*/

/*==========================================================================*/
/* II) Isokinetic Constraint, Nose-like portion                             */
  if( (int_res_tra==0) && (int_res_ter==0) ){
      apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                         &(class->therm_info_class),&(class->therm_class),
                         &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
           /* kinetic=0.0 ;
           for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        }
        printf("after nhc update,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
  }else{
      ix_now = 4;
      if((ir_tra==1)){ix_now=3;}
      if((ir_tra==1)&&(ir_tor==1)){ix_now=2;}
      if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}
      if(ix_respa>=ix_now){
          apply_NH_ISOK_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                             &(class->therm_info_class),&(class->therm_class),
                             &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
           /* kinetic=0.0;
            for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        }
        printf("after nhc update,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
      }
  }



/*==========================================================================*/
/* III) Random Force Generation and Application to {v_2j}                   */

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_x[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_x[ipart]]+=(sigma*int_scr_atm_kin[ichain]);
	  }
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_y[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_y[ipart]]+=(sigma*int_scr_atm_kin[ichain]);

	  }
	  gaussran(len_nhc,iseed,iseed2,qseed,int_scr_atm_kin);
	  for(ichain=1;ichain<=len_nhc;ichain++){
		  therm_v2[ichain][inhc_z[ipart]]*=stoch_e;
		  therm_v2[ichain][inhc_z[ipart]]+=(sigma*int_scr_atm_kin[ichain]);

	  }
  }


/*==========================================================================*/
/* IV) Isokinetic Constraint, Force-dependent portion                       */

    	  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    		  a=clatoms_fx[ipart]*clatoms_vx[ipart]/(LkT);
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

    		  a=clatoms_fy[ipart]*clatoms_vy[ipart]/(LkT);
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

    		  a=clatoms_fz[ipart]*clatoms_vz[ipart]/(LkT);
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
   	  /* kinetic=0.0;
	      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
	    	  kinetic+=clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
	      }
        printf("after first v-v1 forces,dti4 %g dti2 %g dti %g ATM KE %g \n",class->therm_info_class.wdti4[1],
                class->therm_info_class.wdti2[1],class->therm_info_class.dti_nhc,
                kinetic);*/
/*==========================================================================*/
/* V) Evolve positions                                                     */


  	  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
  		  clatoms_x[ipart] += clatoms_vx[ipart]*dt;
  		  clatoms_y[ipart] += clatoms_vy[ipart]*dt;
  		  clatoms_z[ipart] += clatoms_vz[ipart]*dt;
  	  }/*endfor*/


/*==========================================================================*/


/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

  get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->ghost_atoms));


/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_npti(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int i,ipart,iflag,ifirst,ix_now;
  double dt2,tol_glob,v_lnv,x_lnv,x_lnv_o;
  double aa,aa2,arg2,poly,bb,dlen;
  double e2,e4,e6,e8;
  int iii;

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
  
  double *hmato              = general_data->baro.hmato;
  double *hmat               = general_data->cell.hmat;
  double *hmati              = general_data->cell.hmati;
  int    iperd               = general_data->cell.iperd;
  double tolshake            = bonded->constrnt.tolshake;

  MPI_Comm comm_forc         = class->communicate.comm_forc;
  int myid_forc              = class->communicate.myid_forc;
  int np_forc                = class->communicate.np_forc;
  int myatm_start            = class->clatoms_info.myatm_start;
  int myatm_end              = class->clatoms_info.myatm_end;
  int iconstrnt              = bonded->constrnt.iconstrnt;
  int int_res_ter            = general_data->timeinfo.int_res_ter;
  int int_res_tra            = general_data->timeinfo.int_res_tra;
  int ix_respa               = general_data->timeinfo.ix_respa;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);
  dt2 = dt/2.0;

/*==========================================================================*/
/* 1)Save positions and class cell/barostat                                */

  general_data->baro.x_lnv_o = general_data->baro.x_lnv;
  for(i=1;i<=9;i++){hmato[i] = hmat[i];}

  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    class_clatoms_xold[ipart] = class_clatoms_x[ipart];
    class_clatoms_yold[ipart] = class_clatoms_y[ipart];
    class_clatoms_zold[ipart] = class_clatoms_z[ipart];
  }/*endfor*/

/*==========================================================================*/
/*  SHAKE/ROLL CONVERGENCE LOOP                                             */ 

  if(iconstrnt==1){
     iflag = 1;
     cpysys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(general_data->par_rahman),&(class->int_scr),iflag);
  }/*endif*/

  tol_glob = tolshake+1.0;
  ifirst   = 1;

  while(tol_glob>=tolshake){
    if(iconstrnt==1){
      getsys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->therm_info_class),&(class->therm_class),
                 &(general_data->baro),
                 &(general_data->par_rahman),&(class->int_scr),iflag);
    }/*endif*/
  
/*==========================================================================*/
/* II) Evolve NHCs and velocities                                           */
   
    if( (int_res_ter==0) && (int_res_tra==0) ){
      apply_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                      &(class->therm_info_class),&(class->therm_class),
                      &(general_data->baro),
                      &(class->int_scr),&(class->class_comm_forc_pkg));
    }else{
      ix_now = 4;
      if((ir_tra==1)){ix_now=3;}
      if((ir_tra==1)&&(ir_tor==1)){ix_now=2;}
      if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}
      if(ix_respa>=ix_now){
        apply_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),
                        &(class->int_scr),&(class->class_comm_forc_pkg));
      }else{
        apply_NHCPI0_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                         &(class->therm_info_class),&(class->therm_class),
                         &(general_data->baro),
                         &(class->int_scr),&(class->class_comm_forc_pkg));
      }/*endif*/
    }/*endif*/
 
/*==========================================================================*/
/* III) Evolve velocities                                                   */

    for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
      class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                /class_clatoms_mass[ipart];
      class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                /class_clatoms_mass[ipart];
      class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                /class_clatoms_mass[ipart];
    }/*endfor*/
  
/*==========================================================================*/
/* IV) Evolve positions                                                     */

    v_lnv = general_data->baro.v_lnv;

    aa   = exp(dt2*v_lnv);
    aa2  = aa*aa;
    arg2 = (v_lnv*dt2)*(v_lnv*dt2); 
    poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
    bb   = aa*poly;
    for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
      class_clatoms_x[ipart] = class_clatoms_xold[ipart]*aa2
                             + class_clatoms_vx[ipart]*bb*dt;
      class_clatoms_y[ipart] = class_clatoms_yold[ipart]*aa2 
                             + class_clatoms_vy[ipart]*bb*dt;
      class_clatoms_z[ipart] = class_clatoms_zold[ipart]*aa2
                             + class_clatoms_vz[ipart]*bb*dt;
     }/*endfor*/

     general_data->baro.roll_scv = bb;
  
/*==========================================================================*/
/* IV) Evolve the cell                                                      */

     x_lnv_o                  = general_data->baro.x_lnv_o;
     x_lnv                    = x_lnv_o + v_lnv*dt;
     general_data->baro.x_lnv = x_lnv;

     dlen = exp(x_lnv-x_lnv_o);
     for(i=1;i<=9;i++){hmat[i]*=dlen;}
     gethinv(hmat,hmati,&(general_data->cell.vol),iperd);

     (general_data->baro.vol)       = (general_data->cell.vol);
     (general_data->par_rahman.vol) = (general_data->cell.vol);

  
/*==========================================================================*/
/* V) Shake/ROLL if necessary                                               */

     tol_glob = 0.0;
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 1;
       shake_control(bonded,
                    &(class->clatoms_info),&(class->clatoms_pos[1]), 
                    &(general_data->cell),&(general_data->ptens),
                    &(general_data->statepoint),
                    &(general_data->baro),&(general_data->par_rahman),
                    &(general_data->stat_avg),dt,&tol_glob,ifirst,
                    &(class->class_comm_forc_pkg),&(class->ewd_scr));
     }/*endif*/ 
     ifirst = 0;
   }/*endwhile:shake roll*/

/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void int_0_to_dt2_nptf(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                       int ir_tra,int ir_tor,int ir_ter,double dt)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,ipart,iflag,ifirst,iproc;
   int iii,jjj=0;
   double dt2,tol_glob;
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

   int iconstrnt      = bonded->constrnt.iconstrnt;
   int myid_forc      = class->communicate.myid_forc;
   int np_forc        = class->communicate.np_forc;
   int myatm_start    = class->clatoms_info.myatm_start;
   int myatm_end      = class->clatoms_info.myatm_end;
   MPI_Comm comm_forc = class->communicate.comm_forc;

   int int_res_ter    = general_data->timeinfo.int_res_ter;
   int int_res_tra    = general_data->timeinfo.int_res_tra;
   int ix_respa       = general_data->timeinfo.ix_respa;

   double *hmato      = general_data->par_rahman.hmato;
   double *hmat       = general_data->cell.hmat;
   double tolshake    = bonded->constrnt.tolshake;
   int hmat_int_typ   = general_data->cell.hmat_int_typ;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  dt2 = dt/2.0;

/*==========================================================================*/
/* I)Save positions and class cell                                         */

  for(i=1;i<=9;i++){hmato[i] = hmat[i];}
  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    class_clatoms_xold[ipart] = class_clatoms_x[ipart];
    class_clatoms_yold[ipart] = class_clatoms_y[ipart];
    class_clatoms_zold[ipart] = class_clatoms_z[ipart];
  }/*endfor*/

/*==========================================================================*/
/*  SHAKE/ROLL CONVERGENCE LOOP                                             */ 

  if(iconstrnt==1){
    iflag = 2;
    cpysys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
               &(class->therm_info_class),&(class->therm_class),
               &(general_data->baro),
               &(general_data->par_rahman),&(class->int_scr),iflag);
  }/*endif*/
  tol_glob = tolshake+1.0;
  ifirst  = 1;
  while(tol_glob>=tolshake){
     if(iconstrnt==1){
        getsys_NPT(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),
                   &(general_data->par_rahman),&(class->int_scr),iflag);
     }/*endif*/

/*==========================================================================*/
/* II) Evolve NHCs and velocities                                           */


#ifdef JUNK
      printf("before apply %d %d %d\n",ir_ter,ir_tor,ir_tra);
      jjj = jjj + 1;
      printf("roll iteration: %d\n",jjj);
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif

    if( (int_res_ter==0) && (int_res_tra==0) ){
        apply_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                        &(class->therm_info_class),&(class->therm_class),
                        &(general_data->baro),
                        &(general_data->par_rahman),&(general_data->cell),
                        &(class->int_scr),&(class->class_comm_forc_pkg));
     }else{
        ix_now = 4;
        if((ir_tra==1)){ix_now=3;}
        if((ir_tra==1)&&(ir_tor==1)){ix_now=2;}
        if((ir_tra==1)&&(ir_tor==1)&&(ir_ter==1)){ix_now=1;}
        if(ix_respa>=ix_now){
          apply_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                          &(class->therm_info_class),&(class->therm_class),
                          &(general_data->baro),
                          &(general_data->par_rahman),
                          &(general_data->cell),&(class->int_scr),
                          &(class->class_comm_forc_pkg));
        }else{
          apply_NHCPF0_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                           &(class->therm_info_class),&(class->therm_class),
                           &(general_data->baro),
                           &(general_data->par_rahman),
                           &(general_data->cell),&(class->int_scr),
                           &(class->class_comm_forc_pkg));
        }/*endif*/
     }/*endelse*/

#ifdef JUNK
      printf("after apply %d %d %d\n",ir_ter,ir_tor,ir_tra);
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                            class->clatoms_pos[1].vx[1],
                                            class->clatoms_pos[1].vy[1],
                                            class->clatoms_pos[1].vz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                  class->therm_class.v_nhc[1][1],
                                  class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                  class->therm_class.x_nhc[1][1],
                                  class->therm_class.x_nhc[2][1]); 
      printf("v_vol_nhc[1],v_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.v_vol_nhc[1],
                                   general_data->baro.v_vol_nhc[2]); 
      printf("x_vol_nhc[1],x_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.x_vol_nhc[1],
                                   general_data->baro.x_vol_nhc[2]); 
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
                                 general_data->par_rahman.vgmat[2],
                                 general_data->par_rahman.vgmat[3]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[4],
                                 general_data->par_rahman.vgmat[5],
                                 general_data->par_rahman.vgmat[6]);
      printf("vgmat  %.13g %.13g %.13g\n",general_data->par_rahman.vgmat[7],
                                 general_data->par_rahman.vgmat[8],
                                 general_data->par_rahman.vgmat[9]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[1],
                                        general_data->cell.hmat[2],
                                        general_data->cell.hmat[3]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[4],
                                        general_data->cell.hmat[5],
                                        general_data->cell.hmat[6]);
      printf("hmat %.13g %.13g %.13g\n",general_data->cell.hmat[7],
                                        general_data->cell.hmat[8],
                                        general_data->cell.hmat[9]);
      scanf("%d",&iii);
#endif
/*==========================================================================*/
/* III) Evolve velocities                                                   */

     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       class_clatoms_vx[ipart] += class_clatoms_fx[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vy[ipart] += class_clatoms_fy[ipart]*dt2
                                 /class_clatoms_mass[ipart];
       class_clatoms_vz[ipart] += class_clatoms_fz[ipart]*dt2
                                 /class_clatoms_mass[ipart];
     }/*endfor*/

/*==========================================================================*/
/* IV) Evolve positions and h-matrix                                        */

     if(hmat_int_typ==0){
        move_pos_box(&(general_data->cell),&(general_data->par_rahman),
                     class_clatoms_x,class_clatoms_y,class_clatoms_z,
                     class_clatoms_vx,class_clatoms_vy,class_clatoms_vz,
                     dt,class->clatoms_info.natm_tot);
     }else{
        move_pos_box_upper(class_clatoms_x,class_clatoms_y,class_clatoms_z,
                           class_clatoms_vx,class_clatoms_vy,class_clatoms_vz,
                           general_data->cell.hmat,general_data->cell.hmati, 
                           general_data->par_rahman.vgmat, 
                           dt,class->clatoms_info.natm_tot,myatm_start,myatm_end);
     }/*endif : hmat_int_typ*/

   
/*==========================================================================*/
/* V) Shake/ROLL if necessary                                               */

     tol_glob = 0.0;
     if(iconstrnt==1){
       (bonded->constrnt.iroll) = 2;
       shake_control(bonded,
                     &(class->clatoms_info),
                     &(class->clatoms_pos[1]), 
                     &(general_data->cell),&(general_data->ptens),
                     &(general_data->statepoint),
                     &(general_data->baro),&(general_data->par_rahman),
                     &(general_data->stat_avg),dt,&tol_glob,ifirst,
                     &(class->class_comm_forc_pkg),&(class->ewd_scr));
     }/*endif*/ 
     ifirst = 0;

   }/*endwhile:shake roll*/

/*==========================================================================*/
/* i) Recalculate positions of ghost atoms                                  */

   get_ghost_pos(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->ghost_atoms));

/*==========================================================================*/
#ifdef DEBUG
   for(iproc=0;iproc<np_forc;iproc++){
      Barrier(comm_forc);
      if(myid_forc==iproc){
	for(i=1;i<=9;i++){printf("%d %.12g %.12g\n",i,
                    general_data->cell.hmat[i],
                    general_data->par_rahman.vgmat[i]);}
      }
   }
   if(myid_forc==0){scanf("%d",&iii);}
   Barrier(comm_forc);
#endif
/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


















