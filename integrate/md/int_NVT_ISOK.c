/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVT_ISOK                                 */
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


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVT_ISOK(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */ 
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
	  int iii;
    int i,ipart,iflag;
    int ifirst                 = 1;
    int exit_flag              = 0;
    double dt,dti,dt2,tol_glob;
    int natm_tot,ichain,inhc;
    int iflag_mass             = 1;
    double ann_rate            = general_data->simopts.ann_rate;
    int num_nhc                = class->therm_info_class.num_nhc;
    int len_nhc                = class->therm_info_class.len_nhc;
    int ir_tra                 = 1;
    int ir_tor                 = 1;
    int ir_ter                 = 1;
    int anneal_opt             = general_data->simopts.anneal_opt;
    double anneal_target_temp  = general_data->simopts.ann_target_temp;
    double **therm_gkt         = class->therm_info_class.gkt;
    double **therm_mass        = class->therm_info_class.mass_nhc;
    double *class_clatoms_vx   = class->clatoms_pos[1].vx;
    double *class_clatoms_vy   = class->clatoms_pos[1].vy;
    double *class_clatoms_vz   = class->clatoms_pos[1].vz;
    double *class_clatoms_fx   = class->clatoms_pos[1].fx;
    double *class_clatoms_fy   = class->clatoms_pos[1].fy;
    double *class_clatoms_fz   = class->clatoms_pos[1].fz;
    double *class_clatoms_mass = class->clatoms_info.mass;
    double cst1,cst2,cst3;
    double **therm_v2          = class->therm_class.v_nhc;
    double **therm_v1          = class->therm_class.x_nhc;
    int *inhc_x                = class->therm_info_class.inhc_x;
    int *inhc_y                = class->therm_info_class.inhc_y;
    int *inhc_z                = class->therm_info_class.inhc_z;
    double lennhc	             = (double)len_nhc;  

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    (general_data->timeinfo.int_res_tra)    = 0;
    (general_data->timeinfo.int_res_ter)    = 0;
    dt  = (general_data->timeinfo.dt);
    dti = dt;
    dt2 = dt/2.0;
    class->therm_info_class.dt_nhc  = dt;
    class->therm_info_class.dti_nhc = dt/( (double)(
                             class->therm_info_class.nres_nhc) );
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);
    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    energy_control(class,bonded,general_data);


/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

    int_dt2_to_dt_nvt_isok(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 4) Scale by annealing factor                                          */

  iflag=0;
  if(anneal_opt == 1){
    anneal_class(class,ann_rate,iflag,iflag_mass,anneal_target_temp,&exit_flag);
    general_data->timeinfo.exit_flag = exit_flag;
  }/*endif*/

/*==========================================================================*/
/* 5) Finalize                                                              */

   int_final_class(class,bonded,general_data,iflag);

/*==========================================================================*/

#ifdef DEBUG_GLENN
    for(iproc=0;iproc<np_forc;iproc++){
      Barrier(comm_forc);
      if(myid_forc==iproc){
       printf("1 %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[1],y[1],z[1],vx[1],vy[1],vz[1],fx[1],fy[1],fz[1]);
       printf("n %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n",
                  x[n],y[n],z[n],vx[n],vy[n],vz[n],fx[n],fy[n],fz[n]);
      }
    }
    if(myid_forc==0){scanf("%d",&iii);}
    Barrier(comm_forc);
#endif
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NH_ISOK_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                       THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                       INT_SCR *int_scr,int iflag_mass,
                       CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    double arg,arg2,aa;               
    int iii;
    double cst1, cst2, cst3;

/* Define local pointers                                          */
      double *int_scr_atm_kin   = int_scr->atm_kin;
      double *int_scr_sc        = int_scr->sc;
      double *int_scr_sc_temp   = int_scr->sc_temp;
      int *inhc_x               = therm_info_class->inhc_x;
      int *inhc_y               = therm_info_class->inhc_y;
      int *inhc_z               = therm_info_class->inhc_z;
      double **therm_f          = therm_class->f_nhc;
      double **therm_v2         = therm_class->v_nhc;
      double **therm_v1         = therm_class->x_nhc;
      double **therm_gkt      	= therm_info_class->gkt;
      double **therm_mass 	    = therm_info_class->mass_nhc;
      double *clatoms_vx        = clatoms_pos->vx;
      double *clatoms_vy        = clatoms_pos->vy;
      double *clatoms_vz        = clatoms_pos->vz;
      double *therm_wdti2    	  = therm_info_class->wdti2;
      double *therm_wdti4    	  = therm_info_class->wdti4;
      double *therm_wdti8     	= therm_info_class->wdti8;
      int len_nhc				        = therm_info_class->len_nhc; 
      double lennhc	          	= (double)len_nhc;	     
      double *clatoms_mass		  = clatoms_info->mass;
      int natm_tot              = clatoms_info->natm_tot;;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
      MPI_Comm world            = class_comm_forc_pkg->world;
      double len_fac            = ((double)(len_nhc))/(((double)(len_nhc))+1.0);
      double kinetic            = 0.0;
      double k_BT               = therm_gkt[1][inhc_x[1]];
      double LkT                = lennhc*k_BT;
      double mass_factor        = 1.0;      

/* I) Apply the nhc evolution operator using RESPA                         */


   for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
    for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Calculate forces {G(v_1j)} and apply to {v_2j}                      */
      for(ipart=1;ipart<=natm_tot;ipart++){
    	  for(ichain=1;ichain<=len_nhc;ichain++){
              therm_f[ichain][inhc_x[ipart]] =
            		  therm_mass[ichain][inhc_x[ipart]]*therm_v1[ichain][inhc_x[ipart]]*
            		  therm_v1[ichain][inhc_x[ipart]]-k_BT;
              therm_v2[ichain][inhc_x[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f[ichain][inhc_x[ipart]]/(mass_factor*therm_mass[ichain][inhc_x[ipart]]);

              therm_f[ichain][inhc_y[ipart]] =
            		  therm_mass[ichain][inhc_y[ipart]]*therm_v1[ichain][inhc_y[ipart]]*
            		  therm_v1[ichain][inhc_y[ipart]]-k_BT;
              therm_v2[ichain][inhc_y[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f[ichain][inhc_y[ipart]]/(mass_factor*therm_mass[ichain][inhc_y[ipart]]);

              therm_f[ichain][inhc_z[ipart]] =
            		  therm_mass[ichain][inhc_z[ipart]]*therm_v1[ichain][inhc_z[ipart]]*
            		  therm_v1[ichain][inhc_z[ipart]]-k_BT;
              therm_v2[ichain][inhc_z[ipart]]+=
            		  therm_wdti4[iyosh]*therm_f[ichain][inhc_z[ipart]]/(mass_factor*therm_mass[ichain][inhc_z[ipart]]);

    	  /*endfor chains and endfor atoms*/}   }

/*--------------------------------------------------------------------------*/
/*  2) Calculate and apply Nose-like portion of isokinetic constraint       */

      for(ipart=1;ipart<=natm_tot;ipart++){
    	  aa = 0.0;
    	  arg=(clatoms_vx[ipart]*clatoms_vx[ipart]*clatoms_mass[ipart]);
       	  for(ichain=1;ichain<=len_nhc;ichain++){
       	  aa+=(len_fac*therm_mass[ichain][inhc_x[ipart]]*therm_v1[ichain][inhc_x[ipart]]
    				 *therm_v1[ichain][inhc_x[ipart]]
    				 *exp(-2.0*therm_wdti2[iyosh]*therm_v2[ichain][inhc_x[ipart]]));
       	  }
       	  arg +=aa;
       	  arg2=sqrt(LkT/arg);
    	  clatoms_vx[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1[ichain][inhc_x[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2[ichain][inhc_x[ipart]]);
    	  }

    	  aa = 0.0;
          arg=(clatoms_vy[ipart]*clatoms_vy[ipart]*clatoms_mass[ipart]);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
          aa+=(len_fac*therm_mass[ichain][inhc_y[ipart]]*therm_v1[ichain][inhc_y[ipart]]
    		    	  *therm_v1[ichain][inhc_y[ipart]]
    		   		  *exp(-2.0*therm_wdti2[iyosh]*therm_v2[ichain][inhc_y[ipart]]));
    	  }

    	  arg+=aa;
          arg2=sqrt(LkT/arg);
    	  clatoms_vy[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1[ichain][inhc_y[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2[ichain][inhc_y[ipart]]);
    	  }

    	  aa = 0.0;
          arg=(clatoms_vz[ipart]*clatoms_vz[ipart]*clatoms_mass[ipart]);
    	  for(ichain=1;ichain<=len_nhc;ichain++){
          aa+=(len_fac*therm_mass[ichain][inhc_z[ipart]]*therm_v1[ichain][inhc_z[ipart]]
    		    	  *therm_v1[ichain][inhc_z[ipart]]
    		   		  *exp(-2.0*therm_wdti2[iyosh]*therm_v2[ichain][inhc_z[ipart]]));
    	  }

    	  arg+=aa;
          arg2=sqrt(LkT/arg);
    	  clatoms_vz[ipart]*=arg2;
    	  for(ichain=1;ichain<=len_nhc;ichain++){
    	  therm_v1[ichain][inhc_z[ipart]]*=arg2
    			  *exp(-therm_wdti2[iyosh]*therm_v2[ichain][inhc_z[ipart]]);
    	  }
      }
/*--------------------------------------------------------------------------*/
/*  3) Calculate forces {G(v_1j)} and apply to {v_2j}                      */


      for(ipart=1;ipart<=natm_tot;ipart++){
    	  for(ichain=1;ichain<=len_nhc;ichain++){
              therm_f[ichain][inhc_x[ipart]] =
            		  therm_mass[ichain][inhc_x[ipart]]*therm_v1[ichain][inhc_x[ipart]]*
            		  therm_v1[ichain][inhc_x[ipart]]-k_BT;
              therm_v2[ichain][inhc_x[ipart]]+=
			therm_wdti4[iyosh]*therm_f[ichain][inhc_x[ipart]]/(mass_factor*therm_mass[ichain][inhc_x[ipart]]);

              therm_f[ichain][inhc_y[ipart]] =
            		  therm_mass[ichain][inhc_y[ipart]]*therm_v1[ichain][inhc_y[ipart]]*
            		  therm_v1[ichain][inhc_y[ipart]]-k_BT;
              therm_v2[ichain][inhc_y[ipart]]+=
			therm_wdti4[iyosh]*therm_f[ichain][inhc_y[ipart]]/(mass_factor*therm_mass[ichain][inhc_y[ipart]]);

              therm_f[ichain][inhc_z[ipart]] =
            		  therm_mass[ichain][inhc_z[ipart]]*therm_v1[ichain][inhc_z[ipart]]*
            		  therm_v1[ichain][inhc_z[ipart]]-k_BT;
              therm_v2[ichain][inhc_z[ipart]]+=
			therm_wdti4[iyosh]*therm_f[ichain][inhc_z[ipart]]/(mass_factor*therm_mass[ichain][inhc_z[ipart]]);

    	  /*endfor chains and endfor atoms*/}}

    }} /* endfor iyosh, iresn */

   /*end routine*/}




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NH_ISOK_par(GENERAL_DATA *general_data,CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                      INT_SCR *int_scr, int iflag_mass,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg, VEL_SAMP_CLASS *vel_samp_class)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/









