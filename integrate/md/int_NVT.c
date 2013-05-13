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
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVT(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */ 
    int iii;
    int i,ipart,iflag,ifirst=1;
    int exit_flag=0;
    double dt,dt2,tol_glob;
    int natm_tot,ichain,inhc;
    int iflag_mass = 1;
    double ann_rate = general_data->simopts.ann_rate;
    int     num_nhc = class->therm_info_class.num_nhc;
    int     len_nhc = class->therm_info_class.len_nhc;
    int ir_tra = 1;
    int ir_tor = 1;
    int ir_ter = 1;
    int anneal_opt = general_data->simopts.anneal_opt;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    double **therm_v_nhc    = class->therm_class.v_nhc;
    double **therm_gkt      = class->therm_info_class.gkt;
    double **therm_mass_nhc = class->therm_info_class.mass_nhc;
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

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    (general_data->timeinfo.int_res_tra)    = 0;
    (general_data->timeinfo.int_res_ter)    = 0;
    dt  = (general_data->timeinfo.dt);
    dt2 = (general_data->timeinfo.dt)/2.0;
    class->therm_info_class.dt_nhc  = dt;
    class->therm_info_class.dti_nhc = dt/( (double)(
                             class->therm_info_class.nres_nhc) );
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc ,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);

    zero_constrt_iters(&(general_data->stat_avg));
    general_data->timeinfo.exit_flag = 0;

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_nvt(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);


/*==========================================================================*/
/* 2) Get the new energy/force                                              */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    energy_control(class,bonded,general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

    int_dt2_to_dt_nvt(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

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

void apply_NHC_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
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
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa;                      /* Num: scalar temps      */
    int iii;
    int natm_tot,num_nhc;

/* Define local pointers                                          */
      double *int_scr_atm_kin = int_scr->atm_kin;
      double *int_scr_sc      = int_scr->sc;
      double *int_scr_sc_temp = int_scr->sc_temp;
      int *therm_inhc_x       = therm_info_class->inhc_x;
      int *therm_inhc_y       = therm_info_class->inhc_y;
      int *therm_inhc_z       = therm_info_class->inhc_z;
      double *clatoms_vx      = clatoms_pos->vx;
      double *clatoms_vy      = clatoms_pos->vy;
      double *clatoms_vz      = clatoms_pos->vz;
      double **therm_f_nhc    = therm_class->f_nhc;
      double **therm_v_nhc    = therm_class->v_nhc;
      double **therm_x_nhc    = therm_class->x_nhc;
      double **therm_gkt      = therm_info_class->gkt;
      double **therm_mass_nhc = therm_info_class->mass_nhc;
      double *therm_wdti2     = therm_info_class->wdti2;
      double *therm_wdti4     = therm_info_class->wdti4;
      double *therm_wdti8     = therm_info_class->wdti8;
      double *clatoms_mass;
      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
      MPI_Comm world            = class_comm_forc_pkg->world;

      if(iflag_mass==1){
       clatoms_mass    = clatoms_info->mass;
     }else{
       clatoms_mass    = clatoms_pos->mass;
     }/*endif*/


/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*     and assemble the total particle ke associated                        */
/*     with each chain. The num_nhc+1 thermo is the null thermo             */

    natm_tot = clatoms_info->natm_tot;
    num_nhc = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]       = 0.0;
    /*endfor*/}

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_x[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_y[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_z[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
    /*endfor*/}

 
    if(np_forc > 1 && num_nhc_share > 0){
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
     }/*endfor*/
     Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
     }/*endfor*/
    }/*endif*/

     for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_sc[inhc]       = 1.0;
     /*endfor*/}

  
/*==========================================================================*/
/* III) Get the force on the first NHC in each chain                        */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                             /therm_mass_nhc[1][inhc];
    }/*endfor*/

    for(ichain=1;ichain<=len_nhcm1;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_f_nhc[ichain+1][inhc] = 
               (therm_mass_nhc[ichain][inhc]*
                therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]-
                therm_gkt[ichain+1][inhc])/therm_mass_nhc[ichain+1][inhc];
       } /*endfor*/
    }/* endfor */

/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                              therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            arg = -therm_wdti8[iyosh]*
                   therm_v_nhc[len_nhcp1-ichain][inhc];
            aa = exp(arg);
            therm_v_nhc[len_nhc-ichain][inhc] = 
                  therm_v_nhc[len_nhc-ichain][inhc]*aa*aa
                + therm_wdti4[iyosh]*therm_f_nhc[len_nhc-ichain][inhc]*aa;
          /*endfor*/}
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  3) Evolve the particle velocities (by adding to the scaling factor)     */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          arg = -therm_wdti2[iyosh]*therm_v_nhc[1][inhc];
          aa = exp(arg);
          int_scr_sc[inhc]      *= aa;
          int_scr_atm_kin[inhc] *= aa*aa;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  4) Evolve the therm positions                                           */
        for(ichain=1;ichain<=therm_info_class->len_nhc;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_x_nhc[ichain][inhc] += 
                       therm_v_nhc[ichain][inhc]*therm_wdti2[iyosh];
          /*endfor*/}
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  5) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_f_nhc[1][inhc] = 
                   (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                   /therm_mass_nhc[1][inhc];
        /*endfor*/}
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            arg = -therm_wdti8[iyosh]*therm_v_nhc[ichain+1][inhc];
            aa = exp(arg);
            therm_v_nhc[ichain][inhc] = therm_v_nhc[ichain][inhc]*aa*aa
                      + therm_wdti4[iyosh]*therm_f_nhc[ichain][inhc]*aa;
          /*endfor*/}
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_f_nhc[ichain+1][inhc] = 
                 (therm_mass_nhc[ichain][inhc]*
                  therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]-
                  therm_gkt[ichain+1][inhc])/therm_mass_nhc[ichain+1][inhc];

          /*endfor*/}
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  6) Evolve the last therm velocotiy in each chain                        */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                       therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vx[ipart] *= int_scr_sc[therm_inhc_x[ipart]];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vy[ipart] *= int_scr_sc[therm_inhc_y[ipart]];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vz[ipart] *= int_scr_sc[therm_inhc_z[ipart]];

    /*endfor*/}


/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_GGMT2_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
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
    int len_nhc;                        /* Num: length of chains  */
    double arg,aa;                      /* Num: scalar temps      */
    int iii;
    int natm_tot,num_nhc;
    double cst1,cst2;
    double ka,sn;
    double text;

/* Define local pointers                                          */
      double *int_scr_atm_kin = int_scr->atm_kin;
      double *int_scr_sc      = int_scr->sc;
      double *int_scr_sc_temp = int_scr->sc_temp;
      int *therm_inhc_x       = therm_info_class->inhc_x;
      int *therm_inhc_y       = therm_info_class->inhc_y;
      int *therm_inhc_z       = therm_info_class->inhc_z;
      double *clatoms_vx      = clatoms_pos->vx;
      double *clatoms_vy      = clatoms_pos->vy;
      double *clatoms_vz      = clatoms_pos->vz;
      double **therm_f_nhc    = therm_class->f_nhc;
      double **therm_v_nhc    = therm_class->v_nhc;
      double **therm_x_nhc    = therm_class->x_nhc;
      double **therm_gkt      = therm_info_class->gkt;
      double **therm_mass_nhc = therm_info_class->mass_nhc;
      double *text_nhc        = therm_info_class->text_nhc;
      double *therm_wdti2     = therm_info_class->wdti2;
      double *therm_wdti4     = therm_info_class->wdti4;
      double *therm_wdti8     = therm_info_class->wdti8;
      double *clatoms_mass;

      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;

      if(iflag_mass==1){
       clatoms_mass    = clatoms_info->mass;
     }else{
       clatoms_mass    = clatoms_pos->mass;
     }/*endif*/


/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate GGMT chain             */
/*    and assemble the total particle ke associated                         */
/*    with each chain. The num_nhc+1 thermo is the null thermo              */

    natm_tot = clatoms_info->natm_tot;
    num_nhc = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]      = 0.0;
    /*endfor*/}

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_x[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_y[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_z[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
    /*endfor*/}

 
    if(np_forc > 1 && num_nhc_share > 0){
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
     }/*endfor*/
     Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
     }/*endfor*/
    }/*endif*/

/*==========================================================================*/
/* III) Get the force on each GGMT chain                                   */

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      text = text_nhc[inhc]/BOLTZ;
      cst1 = therm_gkt[1][inhc]/text+2.0;

      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                 /therm_mass_nhc[1][inhc];
      therm_f_nhc[2][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                 /cst1-therm_gkt[1][inhc]*text)/therm_mass_nhc[2][inhc];
    }/*endfor*/

/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
     for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
     
/*==========================================================================*/
/*  1) Evolve the therm velocity in each chain                              */
 
       for (ichain=1;ichain<=len_nhc;ichain++){
         for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_v_nhc[ichain][inhc] +=
                           therm_f_nhc[ichain][inhc]*therm_wdti4[iyosh];
           } /*endfor*/
        } /*endfor*/
     
/*--------------------------------------------------------------------------*/
/*  2) Evolve the particle velocities -- step 1                             */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/
     
/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 2                                 */
     
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = therm_v_nhc[2][therm_inhc_x[ipart]]/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = therm_v_nhc[2][therm_inhc_y[ipart]]/cst1;
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = therm_v_nhc[2][therm_inhc_z[ipart]]/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/

/*======================================================================*/
/* Evolve the particle velocities -- step 3 */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
              +text*therm_v_nhc[2][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/  

/*==========================================================================*/
/*  Update S                                                                */

        for(inhc=1;inhc<=num_nhc+1;inhc++){
          int_scr_atm_kin[inhc] = 0.0;
          int_scr_sc[inhc]      = 0.0;
          /*endfor*/}

        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
          /*endfor*/}
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
          /*endfor*/}
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
          /*endfor*/}

        if(np_forc > 1 && num_nhc_share > 0){
          for(inhc = 1;inhc<=num_nhc_share;inhc++){
            int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
          Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
          for(inhc = 1;inhc<=num_nhc_share;inhc++){
            int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
        }/*endfor*/

/*==========================================================================*/
/*  3) Evolve the therm positions                                           */

       for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         text = text_nhc[inhc]/BOLTZ;
         therm_x_nhc[1][inhc] +=
                    therm_v_nhc[1][inhc]*therm_wdti2[iyosh];
         therm_x_nhc[2][inhc] +=
                   (int_scr_atm_kin[inhc]*text/therm_gkt[1][inhc]+text)
                   *therm_v_nhc[2][inhc]*therm_wdti2[iyosh];
       }/*endfor*/

/*==========================================================================*/
/*  4) Evolve the particle velocities, again, step 1                        */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/

/*========================================================================*/
/*     Obtained the updated S values                                      */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*======================================================================*/
/* Evolve the particle velocities -- step 2                             */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = therm_v_nhc[2][therm_inhc_x[ipart]]/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = therm_v_nhc[2][therm_inhc_y[ipart]]/cst1;
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = therm_v_nhc[2][therm_inhc_z[ipart]]/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti2[iyosh]*ka*sn);
       }/*endfor*/

/*======================================================================*/
/* Evolve the particle velocities -- step 3 */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/ 
     
/*==========================================================================*/
/* Obtain the updated S value again                                         */

      for(inhc=1;inhc<=num_nhc+1;inhc++){
        int_scr_atm_kin[inhc] = 0.0;
        int_scr_sc[inhc]      = 0.0;
        /*endfor*/}

      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
        /*endfor*/}
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
        /*endfor*/}
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
        /*endfor*/}

        if(np_forc > 1 && num_nhc_share > 0){
          for(inhc = 1;inhc<=num_nhc_share;inhc++){
            int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
          Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
          for(inhc = 1;inhc<=num_nhc_share;inhc++){
            int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
        }/*endfor*/

/*==========================================================================*/
/*  5) Get the force for each thermostats                                 */

      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        text = text_nhc[inhc]/BOLTZ;
        cst1 = therm_gkt[1][inhc]/text+2.0;
        therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                              /therm_mass_nhc[1][inhc];
        therm_f_nhc[2][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                              /cst1-therm_gkt[1][inhc]*text)
                              /therm_mass_nhc[2][inhc];
      }/*endfor*/

/*==========================================================================*/
/*  6) Evolve the therm velocity in each chain                              */

       for (ichain=1;ichain<=len_nhc;ichain++){
         for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_v_nhc[ichain][inhc] +=
                           therm_f_nhc[ichain][inhc]*therm_wdti4[iyosh];
           } /*endfor*/
        } /*endfor*/
     
/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
     }/*endfor: iyosh */
    }/*endfor:iresn */

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_GGMT3_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
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
    int len_nhc;                        /* Num: length of chains  */
    double arg,aa;                      /* Num: scalar temps      */
    int iii;
    int natm_tot,num_nhc;
    double cst1,cst2;
    double ka,kaa,sn,sq;
    double text;

/* Define local pointers                                          */
      double *int_scr_atm_kin = int_scr->atm_kin;
      double *int_scr_sc      = int_scr->sc;
      double *int_scr_sc_temp = int_scr->sc_temp;
      int *therm_inhc_x       = therm_info_class->inhc_x;
      int *therm_inhc_y       = therm_info_class->inhc_y;
      int *therm_inhc_z       = therm_info_class->inhc_z;
      double *clatoms_vx      = clatoms_pos->vx;
      double *clatoms_vy      = clatoms_pos->vy;
      double *clatoms_vz      = clatoms_pos->vz;
      double **therm_f_nhc    = therm_class->f_nhc;
      double **therm_v_nhc    = therm_class->v_nhc;
      double **therm_x_nhc    = therm_class->x_nhc;
      double **therm_gkt      = therm_info_class->gkt;
      double **therm_mass_nhc = therm_info_class->mass_nhc;
      double *text_nhc        = therm_info_class->text_nhc;
      double *therm_wdti      = therm_info_class->wdti;
      double *therm_wdti2     = therm_info_class->wdti2;
      double *therm_wdti4     = therm_info_class->wdti4;
      double *therm_wdti8     = therm_info_class->wdti8;
      double *clatoms_mass;

      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;

      if(iflag_mass==1){
       clatoms_mass    = clatoms_info->mass;
     }else{
       clatoms_mass    = clatoms_pos->mass;
     }/*endif*/


/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate GGMT chain             */
/*    and assemble the total particle ke associated                         */
/*    with each chain. The num_nhc+1 thermo is the null thermo              */

    natm_tot = clatoms_info->natm_tot;
    num_nhc = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]      = 0.0;
    /*endfor*/}

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_x[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_y[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      int_scr_atm_kin[therm_inhc_z[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
    /*endfor*/}

 
    if(np_forc > 1 && num_nhc_share > 0){
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
     }/*endfor*/
     Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
     }/*endfor*/
    }/*endif*/

/*==========================================================================*/
/* III) Get the force on each GGMT chain                                   */

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      text = text_nhc[inhc]/BOLTZ;
      cst1 = therm_gkt[1][inhc]/text+2.0;
      cst2 = cst1*(therm_gkt[1][inhc]/text+4.0);

      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                 /therm_mass_nhc[1][inhc];
      therm_f_nhc[2][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                 /cst1-therm_gkt[1][inhc]*text)/therm_mass_nhc[2][inhc];
      therm_f_nhc[3][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                            *int_scr_atm_kin[inhc]/cst2-therm_gkt[1][inhc]
                            *text*text)/therm_mass_nhc[3][inhc];
    }/*endfor*/

/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
     for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
     
/*==========================================================================*/
/*  1) Evolve the therm velocity in each chain                              */
 
       for (ichain=1;ichain<=len_nhc;ichain++){
         for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_v_nhc[ichain][inhc] +=
                           therm_f_nhc[ichain][inhc]*therm_wdti4[iyosh];
           } /*endfor*/
        } /*endfor*/
     
/*--------------------------------------------------------------------------*/
/*  2) Evolve the particle velocities -- step 1                             */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 2                                 */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_x[ipart]]
             +text*therm_v_nhc[3][therm_inhc_x[ipart]])/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_y[ipart]]
             +text*therm_v_nhc[3][therm_inhc_y[ipart]])/cst1; 
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_z[ipart]]
             +text*therm_v_nhc[3][therm_inhc_z[ipart]])/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
     
/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 3                                 */
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_x[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = therm_v_nhc[3][therm_inhc_x[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vx[ipart] = clatoms_vx[ipart]/sqrt(sq);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_y[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = therm_v_nhc[3][therm_inhc_y[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vy[ipart] = clatoms_vy[ipart]/sqrt(sq);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_z[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = therm_v_nhc[3][therm_inhc_z[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vz[ipart] = clatoms_vz[ipart]/sqrt(sq);
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;

       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  Evolve the particle velocities -- step 2 again                          */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_x[ipart]]
             +text*therm_v_nhc[3][therm_inhc_x[ipart]])/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_y[ipart]]
             +text*therm_v_nhc[3][therm_inhc_y[ipart]])/cst1; 
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_z[ipart]]
             +text*therm_v_nhc[3][therm_inhc_z[ipart]])/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
     
/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 1 again                           */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  3) Evolve the therm positions                                           */

       for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         text = text_nhc[inhc]/BOLTZ;
         cst1 = therm_gkt[1][inhc]/text+2.0;
         therm_x_nhc[1][inhc] +=
                    therm_v_nhc[1][inhc]*therm_wdti2[iyosh];
         therm_x_nhc[2][inhc] +=
                   (int_scr_atm_kin[inhc]*text/therm_gkt[1][inhc]+text)
                   *therm_v_nhc[2][inhc]*therm_wdti2[iyosh];
         therm_x_nhc[3][inhc] +=
                   (int_scr_atm_kin[inhc]*text*text/therm_gkt[1][inhc]
                   +int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]*text
                   /(cst1*therm_gkt[1][inhc])+text*text)
                   *therm_v_nhc[3][inhc]*therm_wdti2[iyosh];
       }/*endfor*/

/*==========================================================================*/
/*  4) Evolve the particle velocities -- step 1                             */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 2                                 */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_x[ipart]]
             +text*therm_v_nhc[3][therm_inhc_x[ipart]])/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_y[ipart]]
             +text*therm_v_nhc[3][therm_inhc_y[ipart]])/cst1; 
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_z[ipart]]
             +text*therm_v_nhc[3][therm_inhc_z[ipart]])/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
     
/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 3                                 */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_x[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = therm_v_nhc[3][therm_inhc_x[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vx[ipart] = clatoms_vx[ipart]/sqrt(sq);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_y[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = therm_v_nhc[3][therm_inhc_y[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vy[ipart] = clatoms_vy[ipart]/sqrt(sq);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;
         cst2 = cst1*(therm_gkt[1][therm_inhc_z[ipart]]/text+4.0); 
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = therm_v_nhc[3][therm_inhc_z[ipart]]/cst2;
         sq = sqrt(1.0+therm_wdti[iyosh]*ka*sn*sn);
         clatoms_vz[ipart] = clatoms_vz[ipart]/sqrt(sq);
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;


       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  Evolve the particle velocities -- step 2 again                          */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_x[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_x[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_x[ipart]]
             +text*therm_v_nhc[3][therm_inhc_x[ipart]])/cst1; 
         clatoms_vx[ipart] = clatoms_vx[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_y[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_y[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_y[ipart]]
             +text*therm_v_nhc[3][therm_inhc_y[ipart]])/cst1; 
         clatoms_vy[ipart] = clatoms_vy[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         cst1 = therm_gkt[1][therm_inhc_z[ipart]]/text+2.0;       
         sn = int_scr_atm_kin[therm_inhc_z[ipart]];
         ka = (therm_v_nhc[2][therm_inhc_z[ipart]]
             +text*therm_v_nhc[3][therm_inhc_z[ipart]])/cst1; 
         clatoms_vz[ipart] = clatoms_vz[ipart]/
                             sqrt(1.0+therm_wdti4[iyosh]*ka*sn);
       }/*endfor*/
     
/*--------------------------------------------------------------------------*/
/* Evolve the particle velocities -- step 1 again                           */

       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_x[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_x[ipart]]
                +text*therm_v_nhc[2][therm_inhc_x[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_x[ipart]]);
         aa   = exp(arg);
         clatoms_vx[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_y[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_y[ipart]]
                +text*therm_v_nhc[2][therm_inhc_y[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_y[ipart]]);
         aa   = exp(arg);
         clatoms_vy[ipart] *= aa;
       }/*endfor*/
       for(ipart=1;ipart<=natm_tot;ipart++){
         text = text_nhc[therm_inhc_z[ipart]]/BOLTZ;
         arg  = -therm_wdti8[iyosh]*(therm_v_nhc[1][therm_inhc_z[ipart]]
                +text*therm_v_nhc[2][therm_inhc_z[ipart]]
                +text*text*therm_v_nhc[3][therm_inhc_z[ipart]]);
         aa   = exp(arg);
         clatoms_vz[ipart] *= aa;
       }/*endfor*/

/*--------------------------------------------------------------------------*/
/*     Obtained the updated S values                                        */

       for(inhc=1;inhc<=num_nhc+1;inhc++){
         int_scr_atm_kin[inhc] = 0.0;
         int_scr_sc[inhc]      = 0.0;
       } /*endfor*/

       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_x[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_y[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
         /*endfor*/}
       for(ipart=myatm_start;ipart<=myatm_end;ipart++){
         int_scr_atm_kin[therm_inhc_z[ipart]] +=
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
         /*endfor*/}

      if(np_forc > 1 && num_nhc_share > 0){
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
          }/*endfor*/
        Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
        for(inhc = 1;inhc<=num_nhc_share;inhc++){
          int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
          }/*endfor*/
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  5) Get the force for each thermostats                                   */
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        text = text_nhc[inhc]/BOLTZ;
        cst1 = therm_gkt[1][inhc]/text+2.0;
        cst2 = cst1*(therm_gkt[1][inhc]/text+4.0);

      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                 /therm_mass_nhc[1][inhc];
      therm_f_nhc[2][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                 /cst1-therm_gkt[1][inhc]*text)/therm_mass_nhc[2][inhc];
      therm_f_nhc[3][inhc] = (int_scr_atm_kin[inhc]*int_scr_atm_kin[inhc]
                            *int_scr_atm_kin[inhc]/cst2-therm_gkt[1][inhc]
                            *text*text)/therm_mass_nhc[3][inhc];
     }/*endfor*/

/*--------------------------------------------------------------------------*/
/* 6) Evolve the therm velocity in each chain                              */
        for (ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_v_nhc[ichain][inhc] +=
                              therm_f_nhc[ichain][inhc]*therm_wdti4[iyosh];
            }/*endfor*/
          }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  7) End Respa                                                            */
  }/*endfor iyosh*/
}/*endfor iresn*/

/*--------------------------------------------------------------------------*/
}/*endof routine*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NHC_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
              THERM_INFO *therm_info_class,THERM_POS *therm_class, 
              INT_SCR *int_scr, int iflag_mass,
              CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,inhc,ichain;              /* Num: for loop counters */
    int len_nhcm1;    /* Num: length of chains  */
    int natm_tot;
    int num_nhc,iii;
    double cst1,cst2;
    double text;
    double *mass;
    int myatm_start = clatoms_info->myatm_start;
    int myatm_end = clatoms_info->myatm_end;
    int mytherm_start = therm_info_class->mytherm_start;
    int mytherm_end = therm_info_class->mytherm_end;
    int num_nhc_share = therm_info_class->num_nhc_share;
    int *map_share    = therm_info_class->map_share;
    int np_forc       = class_comm_forc_pkg->num_proc;
    MPI_Comm comm_forc   = class_comm_forc_pkg->comm;

    if(iflag_mass==1){
     mass = clatoms_info->mass;
    }else{
     mass = clatoms_pos->mass;
   }/*endif*/

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain     */
/*     and assemble the total particle ke associated                        */
/*     with each chain. The num_nhc+1 thermo is the null thermo             */

    natm_tot = clatoms_info->natm_tot;
    num_nhc = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr->atm_kin[inhc] = 0.0;
    /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_x[ipart]] += 
                  mass[ipart]*
                  clatoms_pos->vx[ipart]*clatoms_pos->vx[ipart];
     /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_y[ipart]] += 
                  mass[ipart]*
                  clatoms_pos->vy[ipart]*clatoms_pos->vy[ipart];
     /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_z[ipart]] += 
                  mass[ipart]*
                  clatoms_pos->vz[ipart]*clatoms_pos->vz[ipart];
     /*endfor*/}

    if(np_forc > 1){
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr->sc[inhc] = int_scr->atm_kin[map_share[inhc]];
     }/*endfor*/
     Allreduce(&(int_scr->sc[1]), &(int_scr->sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
      int_scr->atm_kin[map_share[inhc]] = int_scr->sc_temp[inhc];
     }/*endfor*/
    }/*endif*/
  
/*==========================================================================*/
/* II) Get the force on the first NHC chain                         */

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      therm_class->f_nhc[1][inhc] = (int_scr->atm_kin[inhc]
                             -therm_info_class->gkt[1][inhc])
                             /therm_info_class->mass_nhc[1][inhc];
    /*endfor*/}
  
    if(therm_info_class->therm_typ == 1) {
/*==========================================================================*/
/* III) Get the force on the rest of the NHC in each chain                  */

     len_nhcm1 = (therm_info_class->len_nhc)-1;        
     for(ichain=1;ichain<=len_nhcm1;ichain++){
       for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         therm_class->f_nhc[ichain+1][inhc] = 
                  (therm_info_class->mass_nhc[ichain][inhc]*
                  therm_class->v_nhc[ichain][inhc]*
                  therm_class->v_nhc[ichain][inhc]-
                  therm_info_class->gkt[ichain+1][inhc])/
                  therm_info_class->mass_nhc[ichain+1][inhc];
       /*endfor*/}
     /*endfor*/}

/*==========================================================================*/
/* III) Get the force on the rest of the GGMT in each chain                  */
    } else {
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        text = therm_info_class->text_nhc[inhc]/BOLTZ;
        cst1 = therm_info_class->gkt[1][inhc]/text+2.0;        
        therm_class->f_nhc[2][inhc] = (int_scr->atm_kin[inhc]
           *int_scr->atm_kin[inhc]/cst1-therm_info_class->gkt[1][inhc]
           *text)/therm_info_class->mass_nhc[2][inhc];
    
        if(therm_info_class->len_nhc == 3) {
          cst2 = cst1*(therm_info_class->gkt[1][inhc]/text+4.0);
          therm_class->f_nhc[3][inhc] = (int_scr->atm_kin[inhc]
                *int_scr->atm_kin[inhc]*int_scr->atm_kin[inhc]
                /cst2-therm_info_class->gkt[1][inhc]*text*text)
                /therm_info_class->mass_nhc[3][inhc];
	}/*endif*/
      }/*endfor*/
    }/*endelse*/
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/









