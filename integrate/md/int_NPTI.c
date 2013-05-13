/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NVE                                      */
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
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NPTI(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int i,ipart,iflag,ifirst;
    double dt,dt2,tol_glob;
    double aa,aa2,arg2,poly,bb,dlen;
    double e2,e4,e6,e8;
    int iii;
    int natm_tot;
    int ir_tra = 1;
    int ir_tor = 1;
    int ir_ter = 1;
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
    int myid_forc              = class->communicate.myid_forc;
    e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
    e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    (general_data->timeinfo.int_res_tra)    = 0;
    (general_data->timeinfo.int_res_ter)    = 0;
    dt  = (general_data->timeinfo.dt);
    dt2 = (general_data->timeinfo.dt)/2.0;
    class->therm_info_class.wght    = 1.0;
    class->therm_info_class.dt_nhc  = dt;
    class->therm_info_class.dti_nhc = dt/
                              ( (double)(class->therm_info_class.nres_nhc) );
    set_yosh(class->therm_info_class.nyosh_nhc,
             class->therm_info_class.dti_nhc,
             class->therm_info_class.wdti,class->therm_info_class.wdti2,
             class->therm_info_class.wdti4,class->therm_info_class.wdti8,
             class->therm_info_class.wdti16);

    zero_constrt_iters(&(general_data->stat_avg));

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

    int_0_to_dt2_npti(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/*==========================================================================*/
/* 2) Get the new energy/force                                             */

    (class->energy_ctrl.iget_full_inter)= 1;
    (class->energy_ctrl.iget_res_inter) = 0;
    (class->energy_ctrl.iget_full_intra)= 1;
    (class->energy_ctrl.iget_res_intra) = 0;

    energy_control(class,bonded,general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

    int_dt2_to_dt_npti(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 4) Finalize                                                              */

   iflag=1;
   int_final_class(class,bonded,general_data,iflag);

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
      printf("v_lnv f_lnv_v f_lnv_p %.13g %.13g %.13g\n",
                      general_data->baro.v_lnv,
                      general_data->baro.f_lnv_v/general_data->baro.mass_lnv,
                      general_data->baro.f_lnv_p/general_data->baro.mass_lnv); 
      printf("x_lnv,diff %.13g %.13g \n",general_data->baro.x_lnv,
                       general_data->baro.x_lnv-general_data->baro.x_lnv_o);

      scanf("%d",&iii);
*/
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPI_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                 THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                 BARO *baro, INT_SCR *int_scr,
                 CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa,aa2;                  /* Num: scalar temps      */
    double temp,temp_now;
    double wght;
    int iii,ktemp,i;
    int natm_tot,num_nhc;

/* Define local pointers                                                */
      double *int_scr_atm_kin   = int_scr->atm_kin;
      double *int_scr_sc        = int_scr->sc;
      double *ditherm_nshare_i = therm_info_class->ditherm_nshare_i;
      int *therm_inhc_x         = therm_info_class->inhc_x;
      int *therm_inhc_y         = therm_info_class->inhc_y;
      int *therm_inhc_z         = therm_info_class->inhc_z;
      double *clatoms_mass      = clatoms_info->mass;
      double *clatoms_vx        = clatoms_pos->vx;
      double *clatoms_vy        = clatoms_pos->vy;
      double *clatoms_vz        = clatoms_pos->vz;
      double *clatoms_roll_sc   = clatoms_info->roll_sc;
      double **therm_f_nhc      = therm_class->f_nhc;
      double **therm_v_nhc      = therm_class->v_nhc;
      double **therm_x_nhc      = therm_class->x_nhc;
      double **therm_gkt        = therm_info_class->gkt;
      double **therm_mass_nhc   = therm_info_class->mass_nhc;
      double *therm_wdti2       = therm_info_class->wdti2;
      double *therm_wdti4       = therm_info_class->wdti4;
      double *therm_wdti8       = therm_info_class->wdti8;
      double *baro_v_vol_nhc    = baro->v_vol_nhc;
      double *baro_f_vol_nhc    = baro->f_vol_nhc;
      double *baro_x_vol_nhc    = baro->x_vol_nhc;
      double *baro_mass_vol_nhc = baro->mass_vol_nhc;
      double *baro_gkt_vol      = baro->gkt_vol;
      int mytherm_start         = therm_info_class->mytherm_start;
      int mytherm_end           = therm_info_class->mytherm_end;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int np_forc               = class_comm_forc_pkg->num_proc;
      int num_nhc_share         = therm_info_class->num_nhc_share;
      int *map_share            = therm_info_class->map_share;
      double *int_scr_sc_temp      = int_scr->sc_temp;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain     */
/*    and assemble the total particle ke associated                        */
/*    with each chain. The num_nhc+1 thermo is the null thermo             */
/*    Initialize the v_lnv roll scaling factors                            */

    natm_tot = clatoms_info->natm_tot;
    num_nhc  = therm_info_class->num_nhc;
    baro->roll_scg  = 0.0;
    baro->roll_scg0 = 1.0;
    wght            = therm_info_class->wght;
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
/* III) Get the force on the first NHC in each chain and on the baro stat   */
     temp = 0.0;
    for(i=mytherm_start;i<=mytherm_end;i++){
     temp += int_scr_atm_kin[i]*ditherm_nshare_i[i];
    }/*endfor*/
    if(np_forc > 1){
     temp_now = temp;
     Allreduce(&temp_now, &temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
    }/*endif*/

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                             /therm_mass_nhc[1][inhc];
    /*endfor*/}
    baro_f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -baro_gkt_vol[1])/baro_mass_vol_nhc[1];
    baro->f_lnv_v  = (baro->c2_lnv*temp);

  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        /*endfor*/}
        baro_v_vol_nhc[len_nhc] += 
                            baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            arg = -therm_wdti8[iyosh]*wght*
                   therm_v_nhc[len_nhcp1-ichain][inhc];
            aa = exp(arg);
            therm_v_nhc[len_nhc-ichain][inhc] = 
                therm_v_nhc[len_nhc-ichain][inhc]*aa*aa
              + therm_wdti4[iyosh]*therm_f_nhc[len_nhc-ichain][inhc]*aa*wght;
          /*endfor*/}
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          baro_v_vol_nhc[len_nhc-ichain] = 
                     baro_v_vol_nhc[len_nhc-ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[len_nhc-ichain]*aa;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/* 3) Evolve the dlog(v)/dt                                                 */
        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->roll_scg  = baro->roll_scg*aa2 + therm_wdti4[iyosh]*aa;
        baro->roll_scg0 = baro->roll_scg0*aa2; 
        baro->v_lnv = baro->v_lnv*aa2 + 
                   therm_wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/(baro->mass_lnv);
/*--------------------------------------------------------------------------*/
/*  4) Evolve the particle velocities (by adding to the scaling factor)     */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          arg = -therm_wdti2[iyosh]*
                 (therm_v_nhc[1][inhc]*wght+(baro->v_lnv)*(baro->c2_lnv));
          aa  = exp(arg);
          int_scr_sc[inhc]      *= aa;
          int_scr_atm_kin[inhc] *= aa*aa;
        /*endfor*/}

         temp = 0.0;
        for(i=mytherm_start;i<=mytherm_end;i++){
         temp += int_scr_atm_kin[i]*ditherm_nshare_i[i];
        }/*endfor*/
        if(np_forc > 1){
         temp_now = temp;
         Allreduce(&temp_now, &temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
        }/*endif*/
/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_x_nhc[ichain][inhc] += 
                       therm_v_nhc[ichain][inhc]*therm_wdti2[iyosh]*wght;
          /*endfor*/}
          baro_x_vol_nhc[ichain] += 
                       baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  6) Evolve the dlog(v)/dt                                                */

        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->roll_scg  = baro->roll_scg*aa2 + aa*therm_wdti4[iyosh];
        baro->roll_scg0 = baro->roll_scg0*aa2; 
        baro->f_lnv_v  = (baro->c2_lnv*temp);
        baro->v_lnv = baro->v_lnv*aa2 + 
                   therm_wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/baro->mass_lnv;
/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
           therm_f_nhc[1][inhc] = 
                   (int_scr_atm_kin[inhc]-therm_gkt[1][inhc])
                   /therm_mass_nhc[1][inhc];
        /*endfor*/}
        baro_f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                             -baro_gkt_vol[1])/baro_mass_vol_nhc[1];
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            arg = -therm_wdti8[iyosh]*therm_v_nhc[ichain+1][inhc]*wght;
            aa = exp(arg);
            therm_v_nhc[ichain][inhc] = therm_v_nhc[ichain][inhc]*aa*aa
              + therm_wdti4[iyosh]*therm_f_nhc[ichain][inhc]*aa*wght;
          /*endfor*/}
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_f_nhc[ichain+1][inhc] = 
                 (therm_mass_nhc[ichain][inhc]*
                  therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]-
                  therm_gkt[ichain+1][inhc])/therm_mass_nhc[ichain+1][inhc];

          /*endfor*/}
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                       therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        /*endfor*/}
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */
/*    Store the atm velocity roll scaling factor                            */

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      ktemp = therm_inhc_x[ipart];
      clatoms_vx[ipart] *= int_scr_sc[ktemp];
      clatoms_roll_sc[ipart] = int_scr_sc[ktemp];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vy[ipart] *= int_scr_sc[therm_inhc_y[ipart]];
    /*endfor*/}
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vz[ipart] *= int_scr_sc[therm_inhc_z[ipart]];
    /*endfor*/}
/*==========================================================================*/
/* VI) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    baro->v_lnv_g   = baro->v_lnv;       
    baro->roll_scg *= (2.0/(therm_info_class->dt_nhc));

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPI0_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                  THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                  BARO *baro, INT_SCR *int_scr,
                  CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa,aa2;                  /* Num: scalar temps      */
    double atm_kin,atm_kin_temp,atm_sc;
    int natm_tot,num_nhc;
/* Define local pointers                                                */
      double *clatoms_mass      = clatoms_info->mass;
      double *clatoms_vx        = clatoms_pos->vx;
      double *clatoms_vy        = clatoms_pos->vy;
      double *clatoms_vz        = clatoms_pos->vz;
      double *clatoms_roll_sc   = clatoms_info->roll_sc;
      double *therm_wdti2       = therm_info_class->wdti2;
      double *therm_wdti4       = therm_info_class->wdti4;
      double *therm_wdti8       = therm_info_class->wdti8;
      double *baro_v_vol_nhc    = baro->v_vol_nhc;
      double *baro_f_vol_nhc    = baro->f_vol_nhc;
      double *baro_x_vol_nhc    = baro->x_vol_nhc;
      double *baro_mass_vol_nhc = baro->mass_vol_nhc;
      double *baro_gkt_vol      = baro->gkt_vol;
      int myatm_start           = clatoms_info->myatm_start;
      int myatm_end             = clatoms_info->myatm_end;
      int np_forc               = class_comm_forc_pkg->num_proc;
      MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain     */
/*    and assemble the total particle ke associated                        */
/*    with each chain. The num_nhc+1 thermo is the null thermo             */
/*    Initialize the v_lnv roll scaling factors                            */

    natm_tot = clatoms_info->natm_tot;
    num_nhc  = therm_info_class->num_nhc;
    baro->roll_scg  = 0.0;
    baro->roll_scg0 = 1.0;
    atm_sc          = 1.0;
    atm_kin         = 0.0;
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      atm_kin += 
                (clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart]
                +clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart]
                +clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart]);
    /*endfor*/}

    if(np_forc > 1){
     atm_kin_temp = atm_kin;
     Allreduce(&atm_kin_temp, &atm_kin,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
    }/*endif*/

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain and on the baro stat   */

    baro_f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -baro_gkt_vol[1])/baro_mass_vol_nhc[1];
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        baro_v_vol_nhc[len_nhc] += 
                              baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          baro_v_vol_nhc[len_nhc-ichain] = 
                     baro_v_vol_nhc[len_nhc-ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[len_nhc-ichain]*aa;
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/* 3) Evolve the dlog(v)/dt                                                 */
        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->roll_scg  = baro->roll_scg*aa2 + therm_wdti4[iyosh]*aa;
        baro->roll_scg0 = baro->roll_scg0*aa2; 
        baro->f_lnv_v  = (baro->c2_lnv*atm_kin);
        baro->v_lnv = baro->v_lnv*aa2 + 
                   therm_wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/(baro->mass_lnv);
/*--------------------------------------------------------------------------*/
/*  4) Evolve the particle velocities (by adding to the scaling factor)     */
         arg = -(therm_wdti2[iyosh])*(baro->v_lnv)*(baro->c2_lnv);
         aa  = exp(arg);
         atm_sc  *= aa;
         atm_kin *= aa*aa;
/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          baro_x_vol_nhc[ichain] += 
                       baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  6) Evolve the dlog(v)/dt                                                */
        baro->f_lnv_v  = (baro->c2_lnv*atm_kin);
        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->roll_scg  = baro->roll_scg*aa2 + aa*therm_wdti4[iyosh];
        baro->roll_scg0 = baro->roll_scg0*aa2; 
        baro->v_lnv = baro->v_lnv*aa2 + 
                   therm_wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/baro->mass_lnv;
/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        baro_f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                             -baro_gkt_vol[1])/baro_mass_vol_nhc[1];
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        /*endfor*/}
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */
/*    Store the atm velocity roll scaling factor                            */

    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
      clatoms_vx[ipart] *= atm_sc;
      clatoms_vy[ipart] *= atm_sc;
      clatoms_vz[ipart] *= atm_sc;
      clatoms_roll_sc[ipart] = atm_sc;
    /*endfor*/}

/*==========================================================================*/
/* VI) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    baro->v_lnv_g   = baro->v_lnv;       
    baro->roll_scg *= (2.0/(therm_info_class->dt_nhc));

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NHCPI_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                BARO *baro, INT_SCR *int_scr,
                CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

    int ipart,inhc,ichain,i;              /* Num: for loop counters */
    int len_nhcm1;    /* Num: length of chains  */
    double temp,temp_now;
    int natm_tot,num_nhc;
    int mytherm_start         = therm_info_class->mytherm_start;
    int mytherm_end           = therm_info_class->mytherm_end;
    int myatm_start           = clatoms_info->myatm_start;
    int myatm_end             = clatoms_info->myatm_end;
    int np_forc               = class_comm_forc_pkg->num_proc;
    int num_nhc_share         = therm_info_class->num_nhc_share;
    int *map_share            = therm_info_class->map_share;
    double *ditherm_nshare_i     = therm_info_class->ditherm_nshare_i;
    MPI_Comm comm_forc        = class_comm_forc_pkg->comm;
/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*    and assemble the total particle ke associated                         */
/*    with each chain. The num_nhc+1 thermo is the null thermo              */

    natm_tot = clatoms_info->natm_tot;
    num_nhc  = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr->atm_kin[inhc] = 0.0;
    /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_x[ipart]] += 
                  clatoms_info->mass[ipart]*
                  clatoms_pos->vx[ipart]*clatoms_pos->vx[ipart];
     /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_y[ipart]] += 
                  clatoms_info->mass[ipart]*
                  clatoms_pos->vy[ipart]*clatoms_pos->vy[ipart];
     /*endfor*/}
     for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       int_scr->atm_kin[therm_info_class->inhc_z[ipart]] += 
                  clatoms_info->mass[ipart]*
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
/* II) Get the force on the first NHC in each chain                         */
/*     and on the baro stat                                                 */

    for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
      therm_class->f_nhc[1][inhc] = (int_scr->atm_kin[inhc]
                             -therm_info_class->gkt[1][inhc])
                             /therm_info_class->mass_nhc[1][inhc];
    /*endfor*/}
    baro->f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -baro->gkt_vol[1])/baro->mass_vol_nhc[1];

    temp = 0.0;
    for(i=mytherm_start;i<=mytherm_end;i++){
      temp += int_scr->atm_kin[i]*ditherm_nshare_i[i];
    }/*endfor*/
    if(np_forc > 1){
     temp_now = temp;
     Allreduce(&temp_now, &temp,1,MPI_DOUBLE,
                  MPI_SUM,0,comm_forc);
    }/*endif*/

    baro->f_lnv_v  = (baro->c2_lnv*temp);
  
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
      baro->f_vol_nhc[ichain+1] = 
                   (baro->mass_vol_nhc[ichain]
                   *baro->v_vol_nhc[ichain]*baro->v_vol_nhc[ichain]
                   -baro->gkt_vol[ichain+1])/baro->mass_vol_nhc[ichain+1];
    /*endfor*/}

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/










