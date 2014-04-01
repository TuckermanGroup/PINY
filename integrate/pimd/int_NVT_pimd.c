/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVT_pimd                                 */
/*                                                                          */
/*     This subprogram integrates the system using Vel Verlet RESPA         */
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
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NVT_pimd(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                  ANALYSIS *analysis)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,iflag,ip;
    double dt,dti2,dti;
    int anneal_opt  = general_data->simopts.anneal_opt;
    double ann_rate = general_data->simopts.ann_rate;
    double anneal_target_temp = general_data->simopts.ann_target_temp;
    int ir_tra,ir_tor,ir_ter,ir_pimd;
    int ix_now;
    int nres_tra,nres_tor,nres_ter,nres_pimd;
    int int_res_ter,iii,iproc;
    int natm_tot,pi_beads,pi_beads_proc,pi_beads_proc_st,pi_beads_proc_end,rank;
    int iflag_mass;
    int myid = class->communicate.myid;
    int num_proc = class->communicate.np;
    double pkin;
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
    int nfree = class->clatoms_info.nfree;
    MPI_Comm comm_beads = class->communicate.comm_beads;
    MPI_Comm world      = class->communicate.world;
    double kinet_tmp,kinet,temp;

    TIMER_START("integrate NVT PIMD");

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

    general_data->timeinfo.exit_flag = 0;
    natm_tot      = class->clatoms_info.natm_tot;
    int_res_ter   = general_data->timeinfo.int_res_ter;
    nres_tra      = (general_data->timeinfo.nres_tra);
    nres_tor      = (general_data->timeinfo.nres_tor);
    nres_ter      = (general_data->timeinfo.nres_ter);
    nres_pimd     = (general_data->timeinfo.nres_pimd);
    class->clatoms_info.wght_pimd = 1.0;
    pi_beads      = class->clatoms_info.pi_beads;
    pi_beads_proc_st  = class->clatoms_info.pi_beads_proc_st;
    pi_beads_proc_end = class->clatoms_info.pi_beads_proc_end;
    pi_beads_proc    = class->clatoms_info.pi_beads_proc;
    rank          = class->communicate.myid;
    dt            = (general_data->timeinfo.dt);
    dti           = dt/((double)(nres_ter*nres_tor*nres_tra*nres_pimd));
    dti2          = dti/2.0;
    if(general_data->timeinfo.ix_respa==1){
      class->therm_info_class.dt_nhc   = dt;
      class->therm_info_bead.dt_nhc    = dt;
    }/*endif*/
    if(general_data->timeinfo.ix_respa==2){
      class->therm_info_class.dt_nhc  = dt/((double)(nres_ter));
      class->therm_info_bead.dt_nhc  = dt/((double)(nres_ter));
    }/*endif*/
    if(general_data->timeinfo.ix_respa==3){
      class->therm_info_class.dt_nhc  = dt/((double)(nres_ter*nres_tor));
      class->therm_info_bead.dt_nhc  = dt/((double)(nres_ter*nres_tor));
    }/*endif*/
    if(general_data->timeinfo.ix_respa==4){
      class->therm_info_class.dt_nhc  = dt/((double)(
                                         nres_ter*nres_tor*nres_tra));
      class->therm_info_bead.dt_nhc  = dt/((double)(
                                         nres_ter*nres_tor*nres_tra));
    }/*endif*/
    if(general_data->timeinfo.ix_respa==5){
      class->therm_info_class.dt_nhc  = dt/((double)(
                             nres_ter*nres_tor*nres_tra*nres_pimd));
      class->therm_info_bead.dt_nhc  = dt/((double)(
                             nres_ter*nres_tor*nres_tra*nres_pimd));
    }/*endif*/
    class->therm_info_class.dti_nhc  = (class->therm_info_class.dt_nhc)/
                              ( (double)(class->therm_info_class.nres_nhc) );
    class->therm_info_bead.dti_nhc  = (class->therm_info_bead.dt_nhc)/
                              ( (double)(class->therm_info_bead.nres_nhc) );

    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
             class->therm_info_class.wdti2,class->therm_info_class.wdti4,
             class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);
    set_yosh(class->therm_info_bead.nyosh_nhc,
             class->therm_info_bead.dti_nhc,class->therm_info_bead.wdti,
             class->therm_info_bead.wdti2,class->therm_info_bead.wdti4,
             class->therm_info_bead.wdti8,
             class->therm_info_bead.wdti16);

    zero_constrt_iters(&(general_data->stat_avg));

/*==========================================================================*/
/*==========================================================================*/
/* I) Loop over inter RESPA                                                 */

  for(ir_ter=1;ir_ter<=nres_ter;ir_ter++){

/*==========================================================================*/
/*==========================================================================*/
/* II) Loop over tors RESPA                                                 */

    for(ir_tor=1;ir_tor<=nres_tor;ir_tor++){

/*==========================================================================*/
/*==========================================================================*/
/* III) Loop over intra RESPA                                               */

      for(ir_tra=1;ir_tra<=nres_tra;ir_tra++){

/*==========================================================================*/
/*==========================================================================*/
/* IV) Loop over bead RESPA                                               */

        for(ir_pimd=1;ir_pimd<=nres_pimd;ir_pimd++){

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

           int_0_to_dt2_nvt_pimd(class,bonded,general_data,ir_tra,ir_tor,
                               ir_ter,ir_pimd,dti);

/*==========================================================================*/
/* 2) Get the new energy/force                                              */

      if(ir_pimd==nres_pimd){
          (class->energy_ctrl.iget_full_inter) = 0;
          (class->energy_ctrl.iget_res_inter) = 0;
          (class->energy_ctrl.iget_full_intra) = 0;
          (class->energy_ctrl.iget_res_intra)  = 0;
          (class->energy_ctrl.iget_res_pimd)  = 1;
          if((ir_ter==nres_ter)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra)
                               &&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_full_inter) = 1;}
          if((int_res_ter==1)&&(ir_tor==nres_tor)&&(ir_tra==nres_tra)
                             &&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_res_inter) = 1;}
          if((ir_tra==nres_tra)&&(ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_full_intra) = 1;}
          if(ir_pimd == nres_pimd)
             {(class->energy_ctrl.iget_res_intra) = 1;}

          energy_control_pimd(class,bonded,general_data);

      }else{
        for(ip=1;ip<=pi_beads_proc;ip++){
          class_clatoms_fx   = class->clatoms_pos[ip].fx;
          class_clatoms_fy   = class->clatoms_pos[ip].fy;
          class_clatoms_fz   = class->clatoms_pos[ip].fz;
          class_clatoms_fxm   = class->clatoms_pos[ip].fxm;
          class_clatoms_fym   = class->clatoms_pos[ip].fym;
          class_clatoms_fzm   = class->clatoms_pos[ip].fzm;
          for(ipart=1;ipart<=(natm_tot);ipart++){
            class_clatoms_fx[ipart] = class_clatoms_fxm[ipart];
            class_clatoms_fy[ipart] = class_clatoms_fym[ipart];
            class_clatoms_fz[ipart] = class_clatoms_fzm[ipart];
          }/*endfor*/
        }/*endfor*/
      }/*endif*/
/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

           int_dt2_to_dt_nvt_pimd(class,bonded,general_data,ir_tra,ir_tor,
                               ir_ter,ir_pimd,dti); 

/*==========================================================================*/
    /*endfor:ir_pimd*/}

   /*endfor:ir_tra*/}
  /*endfor:ir_tor*/}
 /*endfor:ir_ter*/}

/*==========================================================================*/
/*==========================================================================*/
/* IV) Scale by annealing factor and check for exit condition               */

     iflag = 0;
     if(anneal_opt == 1){
       kinet_tmp = 0.0;
       if(pi_beads_proc_st == 1){
         iflag_mass=1;
         anneal_bead(&(class->clatoms_info),&(class->clatoms_pos[1]),
                     &(class->int_scr),&(class->class_comm_forc_pkg),
                     &(class->therm_info_class),&(class->therm_class),
                     ann_rate,iflag,iflag_mass,1,&kinet_tmp);
         iflag_mass=2;
         for(ip=2;ip<=pi_beads_proc;ip++){
           anneal_bead(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                       &(class->int_scr),&(class->class_comm_forc_pkg),
                       &(class->therm_info_bead),&(class->therm_bead[ip]),
                       ann_rate,iflag,iflag_mass,ip,&kinet_tmp);
         }/* endfor */
       } else {
         iflag_mass=2;
         for(ip=1;ip<=pi_beads_proc;ip++){
           anneal_bead(&(class->clatoms_info),&(class->clatoms_pos[ip]),
                       &(class->int_scr),&(class->class_comm_forc_pkg),
                       &(class->therm_info_bead),&(class->therm_bead[ip]),
                       ann_rate,iflag,iflag_mass,ip,&kinet_tmp);
         }/* endfor */
       }/* endif */

      if(num_proc > 1){
          Allreduce(&kinet_tmp,&kinet,1,MPI_DOUBLE,MPI_SUM,0,world);
       } else {
          kinet = kinet_tmp;
       }/* endif */

      temp = kinet*BOLTZ/((double) (nfree*pi_beads));
      if(num_proc > 1) Bcast(&temp,1,MPI_DOUBLE,0,world);

      if(ann_rate > 1.0 && temp > anneal_target_temp) {
          general_data->timeinfo.exit_flag = 1;
          if(myid == 0){
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf("Target temperature reached in annealing run\n");
            printf("Performing to last step and preparing to exit\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          }/* endif */
      }/* endif */
      if(ann_rate < 1.0 && temp < anneal_target_temp) {
          general_data->timeinfo.exit_flag = 1;
          if(myid == 0){
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
            printf("Target temperature reached in annealing run\n");
            printf("Performing to last step and preparing to exit\n");
            printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          }/* endif */
      }/* endif */
 
     }/* endif annealing */


/*==========================================================================*/
/* V) Get Kinetic energy and kinetic energy tensor                         */

       get_tvten_pimd(class,general_data);


/*==========================================================================*/
/* VI) Get NHC contribution to energy                                       */

    nhc_vol_potkin_pimd(class,general_data,iflag);

/*==========================================================================*/
/*
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].x[1],
                                         class->clatoms_pos[1].y[1],
                                         class->clatoms_pos[1].z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].vx[1],
                                         class->clatoms_pos[1].vy[1],
                                         class->clatoms_pos[1].vz[1]); 
      printf("fx(1),fy(1),fz(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[1].fx[1]/
                                         class->clatoms_pos[1].mass[1],
                                         class->clatoms_pos[1].fy[1]/
                                         class->clatoms_pos[1].mass[1],
                                         class->clatoms_pos[1].fz[1]/ 
                                         class->clatoms_pos[1].mass[1]);
      printf("fx(1),fy(1),fz(1) %.13g %.13g %.13g\n",
                                         class->clatoms_pos[2].fx[1]/
                                         class->clatoms_pos[2].mass[1],
                                         class->clatoms_pos[2].fy[1]/
                                         class->clatoms_pos[2].mass[1],
                                         class->clatoms_pos[2].fz[1]/ 
                                         class->clatoms_pos[2].mass[1]);
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                         class->therm_class.v_nhc[1][1],
                                         class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                         class->therm_class.x_nhc[1][1],
                                         class->therm_class.x_nhc[2][1]); 
      printf("mass_nhc[1][1],mass_nhc[2][1] %.13g %.13g\n",
                          class->therm_info_class.mass_nhc[1][1],
                          class->therm_info_class.mass_nhc[2][1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                          class->therm_bead[2].v_nhc[1][1],
                          class->therm_bead[2].v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                          class->therm_bead[2].x_nhc[1][1],
                          class->therm_bead[2].x_nhc[2][1]); 
      printf("mass_nhc[1][1],mass_nhc[2][1] %.13g %.13g\n",
                          class->therm_info_bead.mass_nhc[1][1],
                          class->therm_info_bead.mass_nhc[2][1]); 
     pkin = 14.0*32.0*180.0/BOLTZ;
     printf("pressure %g %g %g\n",
    (general_data->ptens.pvten_tot[1]+pkin)/(general_data->cell.vol*PCONV),
     general_data->ptens.pvten_tot[2]/(general_data->cell.vol*PCONV),
     general_data->ptens.pvten_tot[3]/(general_data->cell.vol*PCONV));
    printf("pressure %g %g %g\n",
       general_data->ptens.pvten_tot[4]/(general_data->cell.vol*PCONV),
      (general_data->ptens.pvten_tot[5]+pkin)/(general_data->cell.vol*PCONV),
       general_data->ptens.pvten_tot[6]/(general_data->cell.vol*PCONV));
     printf("pressure %g %g %g\n",
     general_data->ptens.pvten_tot[7]/(general_data->cell.vol*PCONV),
     general_data->ptens.pvten_tot[8]/(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[9]+pkin)/(general_data->cell.vol*PCONV));
     printf("pressure %g %g %g\n",
    (general_data->ptens.pvten_tot[1]+general_data->ptens.tvten[1])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[2]+general_data->ptens.tvten[2])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[3]+general_data->ptens.tvten[3])
       /(general_data->cell.vol*PCONV));
     printf("pressure %g %g %g\n",
    (general_data->ptens.pvten_tot[4]+general_data->ptens.tvten[4])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[5]+general_data->ptens.tvten[5])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[6]+general_data->ptens.tvten[6])
       /(general_data->cell.vol*PCONV));
     printf("pressure %g %g %g\n",
    (general_data->ptens.pvten_tot[7]+general_data->ptens.tvten[7])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[8]+general_data->ptens.tvten[8])
       /(general_data->cell.vol*PCONV),
    (general_data->ptens.pvten_tot[9]+general_data->ptens.tvten[9])
       /(general_data->cell.vol*PCONV));
     printf("tvten %g %g %g %g",general_data->ptens.tvten[1],
                                general_data->ptens.tvten[5],
                                general_data->ptens.tvten[9],pkin);
     scanf("%d",&iii);
*/

  TIMER_STOP("integrate NVT PIMD");

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHC_bead(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
               THERM_INFO *therm_info_class,THERM_POS *therm_class, 
               INT_SCR *int_scr,int iflag_mass)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */


    int ipart,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa,junk;                      /* Num: scalar temps      */
    double arg16;
    int iii,i;
    int natm_tot,num_nhc;

/* Define local pointers                                          */
      double *int_scr_atm_kin = int_scr->atm_kin;
      double *int_scr_sc      = int_scr->sc;
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
      double *therm_wdti16    = therm_info_class->wdti16;
      double *clatoms_mass,wdti2,wdti4,wdti8,wdti16;
      double x_nhc_tot = therm_class->x_nhc_tot;
      double *x_loc;
      double f_nhc_temp = 0.0;
      double *v_nhc_temp;
      double mass_fact,gkt,mass_nhc;
      double *v_nhc;
      double *kinet_nhc_temp;

      TIMER_START("NHC beads");

      if(iflag_mass==1){
       clatoms_mass    = clatoms_info->mass;
     }else{
       clatoms_mass    = clatoms_pos->mass;
     }/*endif*/
     x_nhc_tot = 0.0;


/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*     and assemble the total particle ke associated                        */
/*     with each chain. The num_nhc+1 thermo is the null thermo             */

    natm_tot = clatoms_info->natm_tot;
    num_nhc = therm_info_class->num_nhc;

    kinet_nhc_temp =  (double *)cmalloc(num_nhc*sizeof(double))-1;
    v_nhc_temp     =  (double *)cmalloc(num_nhc*sizeof(double))-1;

    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]       = 1.0;
    /*endfor*/}
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[therm_inhc_x[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
    /*endfor*/}
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[therm_inhc_y[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
    /*endfor*/}
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[therm_inhc_z[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
    /*endfor*/}

  
/*==========================================================================*/
/* III) Get the force on the first NHC in each chain                        */

    for(inhc=1;inhc<=num_nhc;inhc++){
      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                             /therm_mass_nhc[1][inhc];
    /*endfor*/}
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    mass_fact = 1.0/therm_mass_nhc[1][1];  
    mass_nhc  = therm_mass_nhc[1][1];  
    gkt       = therm_gkt[1][1];
  
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
       wdti2  = therm_wdti2[iyosh];
       wdti4  = therm_wdti4[iyosh];
       wdti8  = therm_wdti8[iyosh];
       wdti16 = therm_wdti16[iyosh];


/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */

       for(inhc=1;inhc<=num_nhc;inhc++){
          f_nhc_temp = (mass_nhc*therm_v_nhc[len_nhc-1][inhc]
                                *therm_v_nhc[len_nhc-1][inhc]-gkt)*mass_fact;
          therm_v_nhc[len_nhc][inhc] += 
                              f_nhc_temp*wdti4;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */

        for(ichain=1;ichain<=len_nhcm1;ichain++){
           v_nhc = &(therm_v_nhc[len_nhc-ichain][1]);
          if(ichain>1){
            for(inhc=1;inhc<=num_nhc;inhc++){
              kinet_nhc_temp[inhc] = 
                  mass_nhc*v_nhc[inhc-1]*
                           v_nhc[inhc-1];
              v_nhc_temp[inhc] =  therm_v_nhc[ichain][inhc];
            }
	  }else{
           for(inhc=1;inhc<=num_nhc;inhc++){
              kinet_nhc_temp[inhc] = int_scr_atm_kin[inhc];
            }
	  }
          for(inhc=1;inhc<=num_nhc;inhc++){
            f_nhc_temp = (kinet_nhc_temp[inhc]-gkt)*mass_fact;
            arg = -wdti16*therm_v_nhc[len_nhcp1-ichain][inhc];
            aa =  (1.0 + arg)/(1.0-arg);
            v_nhc_temp[inhc] = 
                  v_nhc[inhc-1]*aa*aa
                + wdti4*f_nhc_temp*aa;
          /*endfor*/}
          for(inhc=1;inhc<=num_nhc;inhc++){
           therm_v_nhc[ichain][inhc] = v_nhc_temp[inhc];
          }          
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  3) Evolve the particle velocities (by adding to the scaling factor)     */

        for(inhc=1;inhc<=num_nhc;inhc++){
          arg = -wdti4*v_nhc_temp[inhc];
          aa =  (1.0 + arg)/(1.0-arg);
          int_scr_sc[inhc]      *= aa;
          int_scr_atm_kin[inhc] *= aa*aa;
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  4) Evolve the therm positions                                           */


        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            x_nhc_tot += 
                       therm_v_nhc[ichain][inhc]*wdti2;
          /*endfor*/}
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  5) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */

        for(ichain=1;ichain<=len_nhcm1;ichain++){
          if(ichain>1){
            for(inhc=1;inhc<=num_nhc;inhc++){
              kinet_nhc_temp[inhc] = 
                  mass_nhc*v_nhc[inhc-1]
                          *v_nhc[inhc-1];
            }
	  }/*endif*/
          v_nhc = &(therm_v_nhc[ichain][1]);
          if(ichain==1){
            for(inhc=1;inhc<=num_nhc;inhc++){
              kinet_nhc_temp[inhc] = int_scr_atm_kin[inhc];
            }
	  }
          for(inhc=1;inhc<=num_nhc;inhc++){
            f_nhc_temp = (kinet_nhc_temp[inhc]-gkt)*mass_fact;
            arg = -wdti16*therm_v_nhc[ichain+1][inhc];
            aa =  (1.0 + arg)/(1.0-arg);
            v_nhc[inhc-1] = v_nhc[inhc-1]*aa*aa
                      + wdti4*f_nhc_temp*aa;
          /*endfor*/}
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  6) Evolve the last therm velocotiy in each chain                        */


        for(inhc=1;inhc<=num_nhc;inhc++){
          f_nhc_temp = (mass_nhc*therm_v_nhc[len_nhc-1][inhc]
                                *therm_v_nhc[len_nhc-1][inhc]-gkt)*mass_fact;
          therm_v_nhc[len_nhc][inhc] += 
                       f_nhc_temp*wdti4;
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
      /*endfor: iyosh */}
    /*endfor: iresn*/}
    therm_class->x_nhc_tot += x_nhc_tot;

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */

    for(ipart=1;ipart<=natm_tot;ipart++){
      clatoms_vx[ipart] *= int_scr_sc[therm_inhc_x[ipart]];
    /*endfor*/}
    for(ipart=1;ipart<=natm_tot;ipart++){
      clatoms_vy[ipart] *= int_scr_sc[therm_inhc_y[ipart]];
    /*endfor*/}
    for(ipart=1;ipart<=natm_tot;ipart++){
      clatoms_vz[ipart] *= int_scr_sc[therm_inhc_z[ipart]];
    /*endfor*/}

    cfree(&(v_nhc_temp[1]));
    cfree(&(kinet_nhc_temp[1]));

    TIMER_STOP("NHC beads");

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHC_bead_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
               THERM_INFO *therm_info_class,THERM_POS *therm_class, 
               INT_SCR *int_scr,int iflag_mass,
               CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

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
      double *therm_wdti16    = therm_info_class->wdti16;
      double *clatoms_mass;
      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int *map_share            = therm_info_class->map_share;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
      double x_nhc_tot = therm_class->x_nhc_tot;
      double wdti2,wdti4,wdti8,wdti16;
      double f_nhc_temp = 0.0;
      double *v_nhc_temp;
      double mass_fact,gkt,mass_nhc;
      double *v_nhc;
      double *kinet_nhc_temp;

      x_nhc_tot = 0.0;

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

    kinet_nhc_temp =  (double *)cmalloc(num_nhc*sizeof(double))-1;
    v_nhc_temp     =  (double *)cmalloc(num_nhc*sizeof(double))-1;

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

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                             /therm_mass_nhc[1][inhc];
    /*endfor*/}
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    mass_fact = 1.0/therm_mass_nhc[1][1];  
    mass_nhc  = therm_mass_nhc[1][1];  
    gkt       = therm_gkt[1][1];
  
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){
       wdti2  = therm_wdti2[iyosh];
       wdti4  = therm_wdti4[iyosh];
       wdti8  = therm_wdti8[iyosh];
       wdti16 = therm_wdti16[iyosh];


/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */

       for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          f_nhc_temp = (mass_nhc*therm_v_nhc[len_nhc-1][inhc]
                                *therm_v_nhc[len_nhc-1][inhc]-gkt)*mass_fact;
          therm_v_nhc[len_nhc][inhc] += 
                              f_nhc_temp*wdti4;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */

        for(ichain=1;ichain<=len_nhcm1;ichain++){
           v_nhc = &(therm_v_nhc[len_nhc-ichain][1]);
          if(ichain>1){
            for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
              kinet_nhc_temp[inhc] = 
                  mass_nhc*v_nhc[inhc-1]*
                           v_nhc[inhc-1];
              v_nhc_temp[inhc] =  therm_v_nhc[ichain][inhc];
            }
	  }else{
           for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
              kinet_nhc_temp[inhc] = int_scr_atm_kin[inhc];
            }
	  }
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            f_nhc_temp = (kinet_nhc_temp[inhc]-gkt)*mass_fact;
            arg = -wdti16*therm_v_nhc[len_nhcp1-ichain][inhc];
            aa =  (1.0 + arg)/(1.0-arg);
            v_nhc_temp[inhc] = 
                  v_nhc[inhc-1]*aa*aa
                + wdti4*f_nhc_temp*aa;
          /*endfor*/}
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_v_nhc[ichain][inhc] = v_nhc_temp[inhc];
          }          
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  3) Evolve the particle velocities (by adding to the scaling factor)     */

        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          arg = -wdti4*v_nhc_temp[inhc];
          aa =  (1.0 + arg)/(1.0-arg);
          int_scr_sc[inhc]      *= aa;
          int_scr_atm_kin[inhc] *= aa*aa;
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  4) Evolve the therm positions                                           */

        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            x_nhc_tot += 
                       therm_v_nhc[ichain][inhc]*wdti2;
          /*endfor*/}
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  5) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */

        for(ichain=1;ichain<=len_nhcm1;ichain++){
          if(ichain>1){
            for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
              kinet_nhc_temp[inhc] = 
                  mass_nhc*v_nhc[inhc-1]
                          *v_nhc[inhc-1];
            }
	  }/*endif*/
          v_nhc = &(therm_v_nhc[ichain][1]);
          if(ichain==1){
            for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
              kinet_nhc_temp[inhc] = int_scr_atm_kin[inhc];
            }
	  }
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            f_nhc_temp = (kinet_nhc_temp[inhc]-gkt)*mass_fact;
            arg = -wdti16*therm_v_nhc[ichain+1][inhc];
            aa =  (1.0 + arg)/(1.0-arg);
            v_nhc[inhc-1] = v_nhc[inhc-1]*aa*aa
                      + wdti4*f_nhc_temp*aa;
          /*endfor*/}
        /*endfor*/}

/*--------------------------------------------------------------------------*/
/*  6) Evolve the last therm velocotiy in each chain                        */

        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          f_nhc_temp = (mass_nhc*therm_v_nhc[len_nhc-1][inhc]
                                *therm_v_nhc[len_nhc-1][inhc]-gkt)*mass_fact;
          therm_v_nhc[len_nhc][inhc] += 
                       f_nhc_temp*wdti4;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/* 7) End Respa                                                             */
      /*endfor: iyosh */}
    /*endfor: iresn*/}
    therm_class->x_nhc_tot += x_nhc_tot;

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

    cfree(&(v_nhc_temp[1]));
    cfree(&(kinet_nhc_temp[1]));


/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





