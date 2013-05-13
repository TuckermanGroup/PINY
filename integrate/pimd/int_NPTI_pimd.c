/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NPTI_res_pimd                            */
/*                                                                          */
/* This subprogram integrates the system using Vel Verlet RESPA             */
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
#include "../proto_defs/proto_math.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void int_NPTI_pimd(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                   ANALYSIS *analysis)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int i,ipart,ip,iflag,iflag_mass;
    double dt,dti2,dti,temp,temp_all;
    int ir_tra,ir_tor,ir_ter,ir_pimd;
    int ix_now;
    int nres_tra,nres_tor,nres_ter,nres_pimd;
    int int_res_ter,iii;
    double aa,aa2,arg2,poly,bb,dlen;
    double e2,e4,e6,e8;
    int natm_tot,pi_beads;
    int pi_beads_proc,rank;
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
    int myatm_start = class->clatoms_info.myatm_start;
    int myatm_end = class->clatoms_info.myatm_end;
    int mytherm_start         = class->therm_info_class.mytherm_start;
    int mytherm_end           = class->therm_info_class.mytherm_end;
    int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
    int np_forc = class->communicate.np_forc;
    int np = class->communicate.np;
    MPI_Comm comm_beads = class->communicate.comm_beads;
    e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
    e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* 0) Useful constants                                                      */

    natm_tot = class->clatoms_info.natm_tot;
    int_res_ter = general_data->timeinfo.int_res_ter;
    nres_tra = (general_data->timeinfo.nres_tra);
    nres_tor = (general_data->timeinfo.nres_tor);
    nres_ter = (general_data->timeinfo.nres_ter);
    nres_pimd = (general_data->timeinfo.nres_pimd);
    pi_beads = class->clatoms_info.pi_beads;
    pi_beads_proc = class->clatoms_info.pi_beads_proc;
    rank = class->communicate.myid;
    dt   = (general_data->timeinfo.dt);
    dti  = dt/((double)(nres_ter*nres_tor*nres_tra*nres_pimd));
    dti2 = dti/2.0;
    if(general_data->timeinfo.ix_respa==1){
      class->therm_info_class.wght  = 
                  (double)(nres_pimd*nres_tra*nres_tor*nres_ter);
      class->therm_info_bead.wght  = 
                  (double)(nres_pimd*nres_tra*nres_tor*nres_ter);
    }/*endif*/
    if(general_data->timeinfo.ix_respa==2){
      class->therm_info_class.wght  = (double)(nres_pimd*nres_tra*nres_tor);
      class->therm_info_bead.wght  = (double)(nres_pimd*nres_tra*nres_tor);
    }/*endif*/
    if(general_data->timeinfo.ix_respa==3){
      class->therm_info_class.wght  =  (double)(nres_pimd*nres_tra);
      class->therm_info_bead.wght  =  (double)(nres_pimd*nres_tra);
    }/*endif*/
    if(general_data->timeinfo.ix_respa==4){
      class->therm_info_class.wght  =  (double)(nres_pimd);
      class->therm_info_bead.wght  =  (double)(nres_pimd);
    }/*endif*/
    if(general_data->timeinfo.ix_respa==5){
      class->therm_info_class.wght  =  1.0;
      class->therm_info_bead.wght  =  1.0;
    }/*endif*/
    class->therm_info_class.dt_nhc   =  dti;
    class->therm_info_bead.dt_nhc   =  dti;
    class->therm_info_class.dti_nhc  = (class->therm_info_class.dt_nhc)/
                             ( (double)(class->therm_info_class.nres_nhc) );
    class->therm_info_bead.dti_nhc  = (class->therm_info_bead.dt_nhc)/
                             ( (double)(class->therm_info_bead.nres_nhc) );
    if(pi_beads_proc_st==1){
      set_yosh(class->therm_info_class.nyosh_nhc,
               class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
               class->therm_info_class.wdti2,class->therm_info_class.wdti4,
               class->therm_info_class.wdti8,
               class->therm_info_class.wdti16);
    }/*endif*/
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

	          int_0_to_dt2_npti_pimd(class,bonded,general_data,
                                  ir_tra,ir_tor,ir_ter,ir_pimd,dti);
				  

/*==========================================================================*/
/*==========================================================================*/
/* 2) Get the new energy/force                                              */


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
          if((ir_pimd==nres_pimd))
             {(class->energy_ctrl.iget_res_intra) = 1;}
          energy_control_pimd(class,bonded,general_data);


/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */


          int_dt2_to_dt_npti_pimd(class,bonded,general_data,
                                  ir_tra,ir_tor,ir_ter,ir_pimd,dti);

/*==========================================================================*/
/*==========================================================================*/
 
         /*endfor:ir_pimd*/}
        /*endfor:ir_tra*/}
      /*endfor:ir_tor*/}
    /*endfor:ir_ter*/}
/*==========================================================================*/
/*==========================================================================*/
/* IV) Get Kinetic energy and kinetic energy tensor                         */


    get_tvten_pimd(class,general_data);


/*==========================================================================*/
/*==========================================================================*/
/* VI) Get NHC contribution to energy                                       */


    iflag = 1;
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
                                         class->clatoms_pos[1].fx[1],
                                         class->clatoms_pos[1].fy[1],
                                         class->clatoms_pos[1].fz[1]); 
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                            class->therm_class.v_nhc[1][1],
                                            class->therm_class.v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                            class->therm_class.x_nhc[1][1],
                                            class->therm_class.x_nhc[2][1]); 
     
      printf("v_nhc[1][1],v_nhc[2,1] %.13g %.13g\n",
                                            class->therm_bead[2].v_nhc[1][1],
                                            class->therm_bead[2].v_nhc[2][1]); 
      printf("x_nhc[1][1],x_nhc[2][1] %.13g %.13g\n",
                                            class->therm_bead[2].x_nhc[1][1],
                                            class->therm_bead[2].x_nhc[2][1]); 
     
      printf("mass_vol_nhc[1],mass_vol_nhc(2) %.13g %.13g\n",
                                   general_data->baro.mass_vol_nhc[1],
                                   general_data->baro.mass_vol_nhc[2]); 
      printf("mass_lnv %.13g\n",  general_data->baro.mass_lnv);
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

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPI_pimd(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                      THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                      THERM_INFO *therm_info_bead,THERM_POS *therm_bead, 
                      BARO *baro, INT_SCR *int_scr,int rank,MPI_Comm world)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int ipart,ip,inhc,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa,aa2;                  /* Num: scalar temps      */
    double temp,temp_loc;
    double wght;
    double v_prod;
    int iii,ktemp;
    int natm_tot,num_nhc,num_nhc_bead;

/* Define local pointers                                                */
      double *int_scr_atm_kin;
      double *int_scr_sc;     
      int *inhc_x;      
      int *inhc_y;      
      int *inhc_z;      
      double *mass;   
      double *vx;     
      double *vy;     
      double *vz;     
      double **f_nhc;   
      double **v_nhc;   
      double **x_nhc;   
      double **gkt ;    
      double **mass_nhc;
      double *wdti2        = therm_info_class->wdti2;
      double *wdti4        = therm_info_class->wdti4;
      double *wdti8        = therm_info_class->wdti8;
      double *v_vol_nhc    = baro->v_vol_nhc;
      double *f_vol_nhc    = baro->f_vol_nhc;
      double *x_vol_nhc    = baro->x_vol_nhc;
      double *mass_vol_nhc = baro->mass_vol_nhc;
      double *gkt_vol      = baro->gkt_vol;
      int pi_beads         = clatoms_info->pi_beads;
      int pi_beads_proc    = clatoms_info->pi_beads_proc;
      int start;

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain     */
/*    and assemble the total particle ke associated                        */
/*    with each chain. The num_nhc+1 thermo is the null thermo             */
/*--------------------------------------------------------------------------*/

    start = (rank == 0 ? 2 : 1);
    natm_tot = clatoms_info->natm_tot;
    num_nhc  = therm_info_class->num_nhc;
    num_nhc_bead  = therm_info_bead->num_nhc;
    wght            = therm_info_class->wght;
/*--------------------------------------------------------------------------*/
/*        Centroid bead                                                     */
 if(rank == 0){
    int_scr_atm_kin = int_scr->therm_scr[1].atm_kin;
    int_scr_sc      = int_scr->therm_scr[1].sc;
    inhc_x    = therm_info_class->inhc_x;
    inhc_y    = therm_info_class->inhc_y;
    inhc_z    = therm_info_class->inhc_z;
    mass    = clatoms_pos[1].mass;
    vx      = clatoms_pos[1].vx;
    vy      = clatoms_pos[1].vy;
    vz      = clatoms_pos[1].vz;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]       = 1.0;
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_x[ipart]] += 
                 mass[ipart]*vx[ipart]*vx[ipart];
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_y[ipart]] += 
                 mass[ipart]*vy[ipart]*vy[ipart];
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_z[ipart]] += 
                 mass[ipart]*vz[ipart]*vz[ipart];
    }/*endfor*/
  }/*endif*/
/*--------------------------------------------------------------------------*/
/*        Other beads                                                      */

    inhc_x    = therm_info_bead->inhc_x;
    inhc_y    = therm_info_bead->inhc_y;
    inhc_z    = therm_info_bead->inhc_z;

/* start = 2 serial or process 0 else start = 1  */
    for(ip=start ;ip<=pi_beads_proc;ip++){
     int_scr_atm_kin = int_scr->therm_scr[ip].atm_kin;
     int_scr_sc      = int_scr->therm_scr[ip].sc;
     mass    = clatoms_pos[ip].mass;
     vx      = clatoms_pos[ip].vx;
     vy      = clatoms_pos[ip].vy;
     vz      = clatoms_pos[ip].vz;
     for(inhc=1;inhc<=num_nhc_bead+1;inhc++){
      int_scr_atm_kin[inhc] = 0.0;
      int_scr_sc[inhc]       = 1.0;
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_x[ipart]] += 
                 mass[ipart]*vx[ipart]*vx[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_y[ipart]] += 
                 mass[ipart]*vy[ipart]*vy[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
      int_scr_atm_kin[inhc_z[ipart]] += 
                 mass[ipart]*vz[ipart]*vz[ipart];
     }/*endfor*/
    }/*endfor*/

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain and on the baro stat   */

  if( rank == 0){
    int_scr_atm_kin = int_scr->therm_scr[1].atm_kin;
    gkt       = therm_info_class->gkt;
    mass_nhc  = therm_info_class->mass_nhc;
    f_nhc     = therm_class->f_nhc;
    for(inhc=1;inhc<=num_nhc;inhc++){
      f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-gkt[1][inhc])
                             /mass_nhc[1][inhc];
    }/*endfor*/
  }/*endif*/

    gkt       = therm_info_bead->gkt;
    mass_nhc  = therm_info_bead->mass_nhc;
    for(ip=start;ip<=pi_beads_proc;ip++){
     int_scr_atm_kin = int_scr->therm_scr[ip].atm_kin;
     f_nhc     = therm_bead[ip].f_nhc;
     for(inhc=1;inhc<=num_nhc_bead;inhc++){
       f_nhc[1][inhc] = (int_scr_atm_kin[inhc]-gkt[1][inhc])
                             /mass_nhc[1][inhc];
     }/*endfor*/
    }/*endfor*/

  if( rank == 0){
    f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -gkt_vol[1])/mass_vol_nhc[1];
    temp = dsum1(num_nhc,int_scr->therm_scr[1].atm_kin,1);
  }else{ temp = 0.0 ; }

    for(ip=start;ip<=pi_beads_proc;ip++){
     temp += dsum1(num_nhc_bead,int_scr->therm_scr[ip].atm_kin,1);
   }/*endfor*/
    temp_loc = temp;
    Reduce(&temp_loc,&temp,1,MPI_DOUBLE,MPI_SUM,0,world);

  if(rank == 0){
    baro->f_lnv_v  = (baro->c2_lnv*temp);
  }/*endif*/  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
    if( rank == 0){
        v_nhc = therm_class->v_nhc;
        f_nhc = therm_class->f_nhc;
        for(inhc=1;inhc<=num_nhc;inhc++){
          v_nhc[len_nhc][inhc] += 
                f_nhc[len_nhc][inhc]*wdti4[iyosh]*wght;
        }/*endfor*/
      }/*endif*/

        for(ip=start;ip<=pi_beads_proc;ip++){
          v_nhc = therm_bead[ip].v_nhc;
          f_nhc = therm_bead[ip].f_nhc;
          for(inhc=1;inhc<=num_nhc_bead;inhc++){
            v_nhc[len_nhc][inhc] += 
                  f_nhc[len_nhc][inhc]*wdti4[iyosh]*wght;
          }/*endfor*/
        }/*endfor*/
    if(rank == 0){
        v_vol_nhc[len_nhc] += 
                            f_vol_nhc[len_nhc]*wdti4[iyosh];
      }/*endif*/
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
      
    if(rank == 0){    
        v_nhc = therm_class->v_nhc;
        f_nhc = therm_class->f_nhc;
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            arg = -wdti8[iyosh]*wght*
                   v_nhc[len_nhcp1-ichain][inhc];
            aa = exp(arg);
            v_nhc[len_nhc-ichain][inhc] = 
                v_nhc[len_nhc-ichain][inhc]*aa*aa
              + wdti4[iyosh]*f_nhc[len_nhc-ichain][inhc]*aa*wght;
          }/*endfor*/
          arg = -wdti8[iyosh]*v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          v_vol_nhc[len_nhc-ichain] = 
                     v_vol_nhc[len_nhc-ichain]*aa*aa
                   + wdti4[iyosh]*f_vol_nhc[len_nhc-ichain]*aa;
        }/*endfor*/
      }/*endif*/

       for(ip=start;ip<=pi_beads_proc;ip++){
        v_nhc = therm_bead[ip].v_nhc;
        f_nhc = therm_bead[ip].f_nhc;
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=1;inhc<=num_nhc_bead;inhc++){
            arg = -wdti8[iyosh]*wght*
                   v_nhc[len_nhcp1-ichain][inhc];
            aa = exp(arg);
            v_nhc[len_nhc-ichain][inhc] = 
                v_nhc[len_nhc-ichain][inhc]*aa*aa
              + wdti4[iyosh]*f_nhc[len_nhc-ichain][inhc]*aa*wght;
          }/*endfor*/
        }/*endfor*/
       }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the dlog(v)/dt                                                 */
 if( rank == 0){
        arg = -wdti8[iyosh]*v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->v_lnv = baro->v_lnv*aa2 + 
                   wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/(baro->mass_lnv);
      }/*endif*/
/*--------------------------------------------------------------------------*/
/*  4) Evolve the particle velocities (by adding to the scaling factor)     */
     
    if ( rank == 0){     
        v_nhc = therm_class->v_nhc;
        v_prod = (baro->v_lnv)*(baro->c2_lnv);
        int_scr_atm_kin = int_scr->therm_scr[1].atm_kin;
        int_scr_sc      = int_scr->therm_scr[1].sc;
        for(inhc=1;inhc<=num_nhc;inhc++){
          arg = -wdti2[iyosh]*
                 (v_nhc[1][inhc]*wght+v_prod);
          aa  = exp(arg);
          int_scr_sc[inhc]      *= aa;
          int_scr_atm_kin[inhc] *= aa*aa;
        }/*endfor*/
      }/*endif*/

#ifdef PARALLEL
   Bcast(&v_prod,1,MPI_DOUBLE,0,world);
#endif
        for(ip=start;ip<=pi_beads_proc;ip++){
          v_nhc = therm_bead[ip].v_nhc;
          int_scr_atm_kin = int_scr->therm_scr[ip].atm_kin;
          int_scr_sc      = int_scr->therm_scr[ip].sc;
          for(inhc=1;inhc<=num_nhc_bead;inhc++){
            arg = -wdti2[iyosh]*
                   (v_nhc[1][inhc]*wght+v_prod);
            aa  = exp(arg);
            int_scr_sc[inhc]      *= aa;
            int_scr_atm_kin[inhc] *= aa*aa;
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
    if( rank == 0){
        x_nhc = therm_class->x_nhc;
        v_nhc = therm_class->v_nhc;
        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            x_nhc[ichain][inhc] += 
                       v_nhc[ichain][inhc]*wdti2[iyosh]*wght;
          }/*endfor*/
          x_vol_nhc[ichain] += 
                       v_vol_nhc[ichain]*wdti2[iyosh];
        }/*endfor*/
      }/*endif*/

        for(ip=start;ip<=pi_beads_proc;ip++){
          x_nhc = therm_bead[ip].x_nhc;
          v_nhc = therm_bead[ip].v_nhc;
          for(ichain=1;ichain<=len_nhc;ichain++){
            for(inhc=1;inhc<=num_nhc_bead;inhc++){
              x_nhc[ichain][inhc] += 
                         v_nhc[ichain][inhc]*wdti2[iyosh]*wght;
            }/*endfor*/
          }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the dlog(v)/dt                                                */
   if( rank == 0){
        temp = dsum1(num_nhc,int_scr->therm_scr[1].atm_kin,1);
   }else{ temp = 0.0; }

        for(ip=start;ip<=pi_beads_proc;ip++){
         temp += dsum1(num_nhc_bead,int_scr->therm_scr[ip].atm_kin,1);
        }/*endfor*/
         temp_loc = temp;
       Reduce(&temp_loc,&temp,1,MPI_DOUBLE,MPI_SUM,0,world);

    if( rank == 0){
        baro->f_lnv_v  = (baro->c2_lnv*temp);
        arg = -wdti8[iyosh]*v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->v_lnv = baro->v_lnv*aa2 + 
                   wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/baro->mass_lnv;
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
  if( rank == 0){
        gkt = therm_info_class->gkt;
        mass_nhc = therm_info_class->mass_nhc;
        x_nhc = therm_class->x_nhc;
        f_nhc = therm_class->f_nhc;
        int_scr_atm_kin = int_scr->therm_scr[1].atm_kin;
        for(inhc=1;inhc<=num_nhc;inhc++){
           f_nhc[1][inhc] = 
                   (int_scr_atm_kin[inhc]-gkt[1][inhc])
                   /mass_nhc[1][inhc];
        }/*endfor*/
        f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                             -gkt_vol[1])/mass_vol_nhc[1];
      }/*endif*/

        gkt = therm_info_bead->gkt;
        mass_nhc = therm_info_bead->mass_nhc;
        for(ip=start;ip<=pi_beads_proc;ip++){
         x_nhc = therm_bead[ip].x_nhc;
         f_nhc = therm_bead[ip].f_nhc;
         int_scr_atm_kin = int_scr->therm_scr[ip].atm_kin;
         for(inhc=1;inhc<=num_nhc_bead;inhc++){
            f_nhc[1][inhc] = 
                    (int_scr_atm_kin[inhc]-gkt[1][inhc])
                    /mass_nhc[1][inhc];
         }/*endfor*/
        }/*endfor*/

   if( rank == 0){
        v_nhc = therm_class->v_nhc;
        f_nhc = therm_class->f_nhc;
        mass_nhc = therm_info_class->mass_nhc;
        gkt = therm_info_class->gkt;
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            arg = -wdti8[iyosh]*v_nhc[ichain+1][inhc]*wght;
            aa = exp(arg);
            v_nhc[ichain][inhc] = v_nhc[ichain][inhc]*aa*aa
              + wdti4[iyosh]*f_nhc[ichain][inhc]*aa*wght;
          }/*endfor*/
          arg = -wdti8[iyosh]*v_vol_nhc[ichain+1];
          aa = exp(arg);
          v_vol_nhc[ichain] = v_vol_nhc[ichain]*aa*aa
                   + wdti4[iyosh]*f_vol_nhc[ichain]*aa;
          for(inhc=1;inhc<=num_nhc;inhc++){
            f_nhc[ichain+1][inhc] = 
                 (mass_nhc[ichain][inhc]*
                  v_nhc[ichain][inhc]*v_nhc[ichain][inhc]-
                  gkt[ichain+1][inhc])/mass_nhc[ichain+1][inhc];

          }/*endfor*/
          f_vol_nhc[ichain+1] = 
                   (mass_vol_nhc[ichain]
                   *v_vol_nhc[ichain]*v_vol_nhc[ichain]
                   -gkt_vol[ichain+1])/mass_vol_nhc[ichain+1];
        }/*endfor*/
      }/*endif*/

        mass_nhc = therm_info_bead->mass_nhc;
        gkt = therm_info_bead->gkt;
        for(ip=start;ip<=pi_beads_proc;ip++){
         v_nhc = therm_bead[ip].v_nhc;
         f_nhc = therm_bead[ip].f_nhc;
         for(ichain=1;ichain<=len_nhcm1;ichain++){
           for(inhc=1;inhc<=num_nhc_bead;inhc++){
             arg = -wdti8[iyosh]*v_nhc[ichain+1][inhc]*wght;
             aa = exp(arg);
             v_nhc[ichain][inhc] = v_nhc[ichain][inhc]*aa*aa
               + wdti4[iyosh]*f_nhc[ichain][inhc]*aa*wght;
           }/*endfor*/
           for(inhc=1;inhc<=num_nhc_bead;inhc++){
             f_nhc[ichain+1][inhc] = 
                  (mass_nhc[ichain][inhc]*
                   v_nhc[ichain][inhc]*v_nhc[ichain][inhc]-
                   gkt[ichain+1][inhc])/mass_nhc[ichain+1][inhc];

           }/*endfor*/
         }/*endfor*/
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
     if(rank == 0){
        v_nhc = therm_class->v_nhc;
        f_nhc = therm_class->f_nhc;
        for(inhc=1;inhc<=num_nhc;inhc++){
          v_nhc[len_nhc][inhc] += 
                f_nhc[len_nhc][inhc]*wdti4[iyosh]*wght;
        }/*endfor*/
      }/*endif*/

        for(ip=start;ip<=pi_beads_proc;ip++){
          v_nhc = therm_bead[ip].v_nhc;
          f_nhc = therm_bead[ip].f_nhc;
          for(inhc=1;inhc<=num_nhc_bead;inhc++){
            v_nhc[len_nhc][inhc] += 
                  f_nhc[len_nhc][inhc]*wdti4[iyosh]*wght;
          }/*endfor*/
        }/*endfor*/
     if(rank == 0){
        v_vol_nhc[len_nhc] += 
                            f_vol_nhc[len_nhc]*wdti4[iyosh];
      }/*endif*/
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */
    for(ip=1;ip<=pi_beads_proc;ip++){
     int_scr_sc      = int_scr->therm_scr[ip].sc;
     if(ip==1&&rank==0){
      inhc_x    = therm_info_class->inhc_x;
      inhc_y    = therm_info_class->inhc_y;
      inhc_z    = therm_info_class->inhc_z;
     }else{
      inhc_x    = therm_info_bead->inhc_x;
      inhc_y    = therm_info_bead->inhc_y;
      inhc_z    = therm_info_bead->inhc_z;
    }/*endif*/
     vx      = clatoms_pos[ip].vx;
     vy      = clatoms_pos[ip].vy;
     vz      = clatoms_pos[ip].vz;
     for(ipart=1;ipart<=natm_tot;ipart++){
       ktemp = inhc_x[ipart];
       vx[ipart] *= int_scr_sc[ktemp];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       vy[ipart] *= int_scr_sc[inhc_y[ipart]];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       vz[ipart] *= int_scr_sc[inhc_z[ipart]];
     }/*endfor*/
    }/*endfor*/
/*--------------------------------------------------------------------------*/
/* broadcast volume and pressure information                                */
#ifdef PARALLEL
   Bcast(&baro->v_lnv,1,MPI_DOUBLE,0,world);
#endif
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPI0_pimd(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                       THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                       THERM_INFO *therm_info_bead,THERM_POS *therm_bead, 
                       BARO *baro, INT_SCR *int_scr,int rank,MPI_Comm world)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int ipart,ip,ichain;              /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhc,len_nhcm1,len_nhcp1;    /* Num: length of chains  */
    double arg,aa,aa2;                  /* Num: scalar temps      */
    double atm_kin,atm_sc;
    int natm_tot,num_nhc;
/* Define local pointers                                                */
      double *mass; 
      double *vx;   
      double *vy;   
      double *vz;   
      double *wdti2        = therm_info_class->wdti2;
      double *wdti4        = therm_info_class->wdti4;
      double *wdti8        = therm_info_class->wdti8;
      double *v_vol_nhc    = baro->v_vol_nhc;
      double *f_vol_nhc    = baro->f_vol_nhc;
      double *x_vol_nhc    = baro->x_vol_nhc;
      double *mass_vol_nhc = baro->mass_vol_nhc;
      double *gkt_vol      = baro->gkt_vol;
      int pi_beads         = clatoms_info->pi_beads;
      int pi_beads_proc    = clatoms_info->pi_beads_proc;

/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain     */
/*    and assemble the total particle ke associated                        */
/*    with each chain. The num_nhc+1 thermo is the null thermo             */

    natm_tot = clatoms_info->natm_tot;
    num_nhc  = therm_info_class->num_nhc;
    atm_sc          = 1.0;
    atm_kin         = 0.0;
    for(ip=1;ip<=pi_beads_proc;ip++){
     mass = clatoms_pos[ip].mass;
     vx = clatoms_pos[ip].vx;
     vy = clatoms_pos[ip].vy;
     vz = clatoms_pos[ip].vz;
     for(ipart=1;ipart<=natm_tot;ipart++){
      atm_kin += 
                (mass[ipart]*vx[ipart]*vx[ipart]
                +mass[ipart]*vy[ipart]*vy[ipart]
                +mass[ipart]*vz[ipart]*vz[ipart]);
     }/*endfor*/
    }/*endfor*/

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain and on the baro stat   */

  if( rank == 0){
    f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -gkt_vol[1])/mass_vol_nhc[1];
    baro->f_lnv_v  = (baro->c2_lnv*atm_kin);
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhc   = (therm_info_class->len_nhc);
    len_nhcm1 = (therm_info_class->len_nhc)-1;
    len_nhcp1 = (therm_info_class->len_nhc)+1;
    for(iresn=1;iresn<=therm_info_class->nres_nhc;iresn++){
      for(iyosh=1;iyosh<=therm_info_class->nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        v_vol_nhc[len_nhc] += 
                              f_vol_nhc[len_nhc]*wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -wdti8[iyosh]*v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          v_vol_nhc[len_nhc-ichain] = 
                     v_vol_nhc[len_nhc-ichain]*aa*aa
                   + wdti4[iyosh]*f_vol_nhc[len_nhc-ichain]*aa;
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the dlog(v)/dt                                                 */
        arg = -wdti8[iyosh]*v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->v_lnv = baro->v_lnv*aa2 + 
                   wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/(baro->mass_lnv);
/*--------------------------------------------------------------------------*/
/*  4) Evolve the particle velocities (by adding to the scaling factor)     */
         arg = -(wdti2[iyosh])*(baro->v_lnv)*(baro->c2_lnv);
         aa  = exp(arg);
         atm_sc  *= aa;
         atm_kin *= aa*aa;
/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          x_vol_nhc[ichain] += 
                       v_vol_nhc[ichain]*wdti2[iyosh];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the dlog(v)/dt                                                */
        baro->f_lnv_v  = (baro->c2_lnv*atm_kin);
        arg = -wdti8[iyosh]*v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        baro->v_lnv = baro->v_lnv*aa2 + 
                   wdti4[iyosh]*
                   (baro->f_lnv_v+baro->f_lnv_p)*aa/baro->mass_lnv;
/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                             -gkt_vol[1])/mass_vol_nhc[1];
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -wdti8[iyosh]*v_vol_nhc[ichain+1];
          aa = exp(arg);
          v_vol_nhc[ichain] = v_vol_nhc[ichain]*aa*aa
                   + wdti4[iyosh]*f_vol_nhc[ichain]*aa;
          f_vol_nhc[ichain+1] = 
                   (mass_vol_nhc[ichain]
                   *v_vol_nhc[ichain]*v_vol_nhc[ichain]
                   -gkt_vol[ichain+1])/mass_vol_nhc[ichain+1];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        v_vol_nhc[len_nhc] += 
                f_vol_nhc[len_nhc]*wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}
  }/*endif for rank*/
#ifdef PARALLEL
  Bcast(&atm_sc,1,MPI_DOUBLE,0,world);
#endif
/*==========================================================================*/
/* V) Apply the accumulated scaling factor to the velocities                */

    for(ip=1;ip<=pi_beads_proc;ip++){
     vx = clatoms_pos[ip].vx;
     vy = clatoms_pos[ip].vy;
     vz = clatoms_pos[ip].vz;
     for(ipart=1;ipart<=natm_tot;ipart++){
      vx[ipart] *= atm_sc;
      vy[ipart] *= atm_sc;
      vz[ipart] *= atm_sc;
     }/*endfor*/
    }/*endfor*/

#ifdef PARALLEL
     Bcast(&baro->v_lnv,1,MPI_DOUBLE,0,world);
#endif
 
/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NHCPI_pimd(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                     THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                     THERM_INFO *therm_info_bead,THERM_POS *therm_bead, 
                     BARO *baro, INT_SCR *int_scr,int rank,
                     MPI_Comm world)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

    int ip,ipart,inhc,ichain;              /* Num: for loop counters */
    int len_nhcm1;    /* Num: length of chains  */
    double temp,temp_loc;
    int natm_tot,num_nhc,num_nhc_bead;
    int pi_beads = clatoms_info->pi_beads;
    int pi_beads_proc  = clatoms_info->pi_beads_proc;
    int start;
        start = (rank == 0 ? 2 : 1);


/*==========================================================================*/
/* I) Map the ke of each particle to the appropriate nose-hoover chain      */
/*    and assemble the total particle ke associated                         */
/*    with each chain. The num_nhc+1 thermo is the null thermo              */

    natm_tot = clatoms_info->natm_tot;
    num_nhc_bead  = therm_info_bead->num_nhc;

  if(rank == 0){
    num_nhc       = therm_info_class->num_nhc;
    for(inhc=1;inhc<=num_nhc+1;inhc++){
      int_scr->therm_scr[1].atm_kin[inhc] = 0.0;
    }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[1].atm_kin[therm_info_class->inhc_x[ipart]] += 
                  clatoms_pos[1].mass[ipart]*
                  clatoms_pos[1].vx[ipart]*clatoms_pos[1].vx[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[1].atm_kin[therm_info_class->inhc_y[ipart]] += 
                  clatoms_pos[1].mass[ipart]*
                  clatoms_pos[1].vy[ipart]*clatoms_pos[1].vy[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[1].atm_kin[therm_info_class->inhc_z[ipart]] += 
                  clatoms_pos[1].mass[ipart]*
                  clatoms_pos[1].vz[ipart]*clatoms_pos[1].vz[ipart];
     }/*endfor*/
  }/*endif*/

    for(ip=start;ip<=pi_beads_proc;ip++){
     for(inhc=1;inhc<=num_nhc_bead+1;inhc++){
       int_scr->therm_scr[ip].atm_kin[inhc] = 0.0;
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[ip].atm_kin[therm_info_bead->inhc_x[ipart]] += 
                  clatoms_pos[ip].mass[ipart]*
                  clatoms_pos[ip].vx[ipart]*clatoms_pos[ip].vx[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[ip].atm_kin[therm_info_bead->inhc_y[ipart]] += 
                  clatoms_pos[ip].mass[ipart]*
                  clatoms_pos[ip].vy[ipart]*clatoms_pos[ip].vy[ipart];
     }/*endfor*/
     for(ipart=1;ipart<=natm_tot;ipart++){
       int_scr->therm_scr[ip].atm_kin[therm_info_bead->inhc_z[ipart]] += 
                  clatoms_pos[ip].mass[ipart]*
                  clatoms_pos[ip].vz[ipart]*clatoms_pos[ip].vz[ipart];
     }/*endfor*/
    }/*endfor*/
/*==========================================================================*/
/* II) Get the force on the first NHC in each chain                         */
/*     and on the baro stat                                                 */

   if( rank == 0){
    for(inhc=1;inhc<=num_nhc;inhc++){
      therm_class->f_nhc[1][inhc] = (int_scr->therm_scr[1].atm_kin[inhc]
                             -therm_info_class->gkt[1][inhc])
                             /therm_info_class->mass_nhc[1][inhc];
    }/*endfor*/
    baro->f_vol_nhc[1] = (baro->mass_lnv*baro->v_lnv*baro->v_lnv
                         -baro->gkt_vol[1])/baro->mass_vol_nhc[1];
    temp = dsum1(num_nhc,int_scr->therm_scr[1].atm_kin,1);
  }else{ temp = 0.0; }

    for(ip=start;ip<=pi_beads_proc;ip++){
     temp += dsum1(num_nhc_bead,int_scr->therm_scr[ip].atm_kin,1);
    }/*endfor*/
       temp_loc = temp;
       Reduce(&temp_loc,&temp,1,MPI_DOUBLE,MPI_SUM,0,world);

       if(rank  == 0){
        baro->f_lnv_v  = (baro->c2_lnv*temp);
      }/*endif*/

    for(ip=start;ip<=pi_beads_proc;ip++){
     for(inhc=1;inhc<=num_nhc_bead;inhc++){
       therm_bead[ip].f_nhc[1][inhc] = (int_scr->therm_scr[ip].atm_kin[inhc]
                              -therm_info_bead->gkt[1][inhc])
                              /therm_info_bead->mass_nhc[1][inhc];
     }/*endfor*/
    }/*endfor*/
  
/*==========================================================================*/
/* III) Get the force on the rest of the NHC in each chain                  */

  if(rank == 0){
    len_nhcm1 = (therm_info_class->len_nhc)-1;        
    for(ichain=1;ichain<=len_nhcm1;ichain++){
      for(inhc=1;inhc<=num_nhc;inhc++){
        therm_class->f_nhc[ichain+1][inhc] = 
                  (therm_info_class->mass_nhc[ichain][inhc]*
                  therm_class->v_nhc[ichain][inhc]*
                  therm_class->v_nhc[ichain][inhc]-
                  therm_info_class->gkt[ichain+1][inhc])/
                  therm_info_class->mass_nhc[ichain+1][inhc];
      }/*endfor*/
      baro->f_vol_nhc[ichain+1] = 
                   (baro->mass_vol_nhc[ichain]
                   *baro->v_vol_nhc[ichain]*baro->v_vol_nhc[ichain]
                   -baro->gkt_vol[ichain+1])/baro->mass_vol_nhc[ichain+1];
    }/*endfor*/
  }/*endify*/

    for(ip=start;ip<=pi_beads_proc;ip++){
     for(ichain=1;ichain<=len_nhcm1;ichain++){
      for(inhc=1;inhc<=num_nhc_bead;inhc++){
         therm_bead[ip].f_nhc[ichain+1][inhc] = 
                   (therm_info_bead->mass_nhc[ichain][inhc]*
                   therm_bead[ip].v_nhc[ichain][inhc]*
                   therm_bead[ip].v_nhc[ichain][inhc]-
                   therm_info_bead->gkt[ichain+1][inhc])/
                   therm_info_bead->mass_nhc[ichain+1][inhc];
      }/*endfor*/
     }/*endfor*/
    }/*endfor*/

/*--------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/

















