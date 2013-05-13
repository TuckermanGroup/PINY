/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_NPTF                                     */
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

void int_NPTF(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*========================================================================*/
  {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int iii;

   int iflag;
   int ir_tra = 1;
   int ir_tor = 1;
   int ir_ter = 1;

   double dt       = general_data->timeinfo.dt;
   int nres_nhc    = class->therm_info_class.nres_nhc;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   (general_data->timeinfo.int_res_tra)    = 0;
   (general_data->timeinfo.int_res_ter)    = 0;

   class->therm_info_class.wght            = 1.0;
   class->therm_info_class.dt_nhc          = dt;
   class->therm_info_class.dti_nhc         = dt/( (double)(nres_nhc) );

   set_yosh(class->therm_info_class.nyosh_nhc,
            class->therm_info_class.dti_nhc,class->therm_info_class.wdti,
            class->therm_info_class.wdti2,class->therm_info_class.wdti4,
            class->therm_info_class.wdti8,
            class->therm_info_class.wdti16);

   zero_constrt_iters(&(general_data->stat_avg));

/*==========================================================================*/
/* 1) Evolve system from t=0 to dt/2                                        */

   int_0_to_dt2_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 2) Get the new energy/force                                             */

   (class->energy_ctrl.iget_full_inter)= 1;
   (class->energy_ctrl.iget_res_inter) = 0;
   (class->energy_ctrl.iget_full_intra)= 1;
   (class->energy_ctrl.iget_res_intra) = 0;

   energy_control(class,bonded,general_data);

/*==========================================================================*/
/* 3) Evolve system from dt/2 to dt                                         */

   int_dt2_to_dt_nptf(class,bonded,general_data,ir_tra,ir_tor,ir_ter,dt);

/*==========================================================================*/
/* 4) Finalize                                                              */

   iflag=2;
   int_final_class(class,bonded,general_data,iflag);

/*==========================================================================*/

#ifdef DEBUG
      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",class->clatoms_pos[1].x[1],
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

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPF(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                 THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                 BARO *baro, 
                 PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr,
                 double dt)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int i,j,joff,n,ipart,inhc,ichain;   /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhcm1,len_nhcp1;            /* Num: length of chains  */
    int ktemp;
    double arg,aa,aa2;                  /* Num: scalar temps      */
    double tempvx,tempvy,tempvz;  
    int iii;

/* Define local pointers                                          */
    int natm_tot                 = clatoms_info->natm_tot;
    double *clatoms_roll_sc      = clatoms_info->roll_sc;

    double *clatoms_vx           = clatoms_pos->vx;
    double *clatoms_vy           = clatoms_pos->vy;
    double *clatoms_vz           = clatoms_pos->vz;

    int num_nhc                  = therm_info_class->num_nhc; 
    int len_nhc                  = therm_info_class->len_nhc;
    int nres_nhc                 = therm_info_class->nres_nhc;
    int nyosh_nhc                = therm_info_class->nyosh_nhc;
    double dt_nhc                = therm_info_class->dt_nhc;
    double     wght              = therm_info_class->wght;
    int *therm_inhc_x            = therm_info_class->inhc_x;
    int *therm_inhc_y            = therm_info_class->inhc_y;
    int *therm_inhc_z            = therm_info_class->inhc_z;
    double **therm_gkt           = therm_info_class->gkt;
    double **therm_mass_nhc      = therm_info_class->mass_nhc;
    double *therm_wdti2          = therm_info_class->wdti2;
    double *therm_wdti4          = therm_info_class->wdti4;
    double *therm_wdti8          = therm_info_class->wdti8;

    double *int_scr_sc           = int_scr->sc;

    double **therm_f_nhc         = therm_class->f_nhc;
    double **therm_v_nhc         = therm_class->v_nhc;
    double **therm_x_nhc         = therm_class->x_nhc;

    double *baro_v_vol_nhc       = baro->v_vol_nhc;
    double *baro_f_vol_nhc       = baro->f_vol_nhc;
    double *baro_x_vol_nhc       = baro->x_vol_nhc;
    double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
    double *baro_gkt_vol         = baro->gkt_vol;

    double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
    double *par_rahman_vgmat     = par_rahman->vgmat;
    double *par_rahman_vgmat_g   = par_rahman->vgmat_g;
    double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
    double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
    double *par_rahman_vtemps    = par_rahman->vtemps;
    double *par_rahman_vtempv    = par_rahman->vtempv;
    double *par_rahman_veig      = par_rahman->veig;
    double *par_rahman_veigv     = par_rahman->veigv;
    double *par_rahman_vexpdt    = par_rahman->vexpdt;
    double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
    double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);
    double mass_hm               = par_rahman->mass_hm;

    int hmat_int_typ             = cell->hmat_int_typ;
    int hmat_cons_typ            = cell->hmat_cons_typ;
    int iperd                    = cell->iperd;

/*==========================================================================*/
/* I) Initialize the roll scaling factors                                   */

    for(i=1;i<=9;i++){par_rahman_roll_mtvv[i] = 0.0;}
    par_rahman_roll_mtvv[1]=1.0;
    par_rahman_roll_mtvv[5]=1.0;
    par_rahman_roll_mtvv[9]=1.0;

    for(ipart=1;ipart<=natm_tot;ipart++){clatoms_roll_sc[ipart] = 1.0;}
    int_scr_sc[num_nhc+1]   = 1.0;
    (*par_rahman_roll_scg)  = 0.0;
    (*par_rahman_roll_scg0) = 1.0;

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain and on the baro stat   */

    get_fnhc1(clatoms_info,clatoms_pos,therm_info_class,therm_class,
              baro,par_rahman);
    get_fgmatv(clatoms_info,clatoms_pos,par_rahman,cell);
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */
 

    len_nhcm1 = len_nhc-1;
    len_nhcp1 = len_nhc+1;

    for(iresn=1;iresn<=nres_nhc;iresn++){
      for(iyosh=1;iyosh<=nyosh_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        for(inhc=1;inhc<=num_nhc;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                        therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        }/*endfor*/
        baro_v_vol_nhc[len_nhc] += 
                            baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocitiy in each chain        */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            arg = -therm_wdti8[iyosh]*wght*
                   therm_v_nhc[len_nhcp1-ichain][inhc];
            aa = exp(arg);
            therm_v_nhc[len_nhc-ichain][inhc] = 
                therm_v_nhc[len_nhc-ichain][inhc]*aa*aa
              + therm_wdti4[iyosh]*therm_f_nhc[len_nhc-ichain][inhc]*aa*wght;
          }/*endfor*/
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          baro_v_vol_nhc[len_nhc-ichain] = 
                     baro_v_vol_nhc[len_nhc-ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[len_nhc-ichain]*aa;
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the vgmat                                                      */

        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = 
              (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 

        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                             therm_wdti4[iyosh]*
                             (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                             /mass_hm;
        }/*endfor*/


/*--------------------------------------------------------------------------*/
/* 4) Evolve the classical velocities                                       */

       if(hmat_int_typ==0){
         move_vel_vbox(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4, therm_wdti8,
                   therm_x_nhc,therm_v_nhc,therm_inhc_x,therm_inhc_y,
                   therm_inhc_z,iyosh,num_nhc, int_scr_sc, 
                   clatoms_roll_sc,wght,len_nhc,1,num_nhc);
       }else{
         move_vel_vbox_upper(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4, therm_wdti8,
                   therm_x_nhc,therm_v_nhc,therm_inhc_x,therm_inhc_y,
                   therm_inhc_z,iyosh,num_nhc, int_scr_sc, 
                   clatoms_roll_sc,wght,len_nhc);

       }/*endif : hmat_int_typ*/

/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            therm_x_nhc[ichain][inhc] += 
                       therm_v_nhc[ichain][inhc]*therm_wdti2[iyosh]*wght;
          }/*endfor*/
          baro_x_vol_nhc[ichain] += 
                       baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the vgmat                                                     */


        get_fgmatv(clatoms_info,clatoms_pos,par_rahman,cell);

        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = 
                  (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 
 

        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                             therm_wdti4[iyosh]*
                             (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                             /(mass_hm);
        }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        get_fnhc1(clatoms_info,clatoms_pos,therm_info_class,therm_class,
                  baro,par_rahman);
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            arg = -therm_wdti8[iyosh]*therm_v_nhc[ichain+1][inhc]*wght;
            aa = exp(arg);
            therm_v_nhc[ichain][inhc] = therm_v_nhc[ichain][inhc]*aa*aa
              + therm_wdti4[iyosh]*therm_f_nhc[ichain][inhc]*aa*wght;
          }/*endfor*/
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          for(inhc=1;inhc<=num_nhc;inhc++){
            therm_f_nhc[ichain+1][inhc] = 
                 (therm_mass_nhc[ichain][inhc]*
                  therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]-
                  therm_gkt[ichain+1][inhc])/therm_mass_nhc[ichain+1][inhc];

          }/*endfor*/
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        for(inhc=1;inhc<=num_nhc;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                       therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        }/*endfor*/
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* VI) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    for(i=1;i<=9;i++){par_rahman_vgmat_g[i]=par_rahman_vgmat[i];}
    (*par_rahman_roll_scg) *= (2.0/dt_nhc);

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPF_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                 THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                 BARO *baro, 
                 PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr,
                 CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
      {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int iii;
    int i,j,ipart,inhc,ichain;          /* Num: for loop counters */
    int iresn,iyosh;                    /* Num: for loop counters */
    int len_nhcm1,len_nhcp1;
    double arg,aa,aa2;                  /* Num: scalar temps      */

/* Define local pointers                                          */

    int natm_tot                 = clatoms_info->natm_tot;
    int myatm_start              = clatoms_info->myatm_start;
    int myatm_end                = clatoms_info->myatm_end;
    double *clatoms_roll_sc      = clatoms_info->roll_sc;

    double *clatoms_vx           = clatoms_pos->vx;
    double *clatoms_vy           = clatoms_pos->vy;
    double *clatoms_vz           = clatoms_pos->vz;

    int num_nhc                  = therm_info_class->num_nhc; 
    int len_nhc                  = therm_info_class->len_nhc;
    int nres_nhc                 = therm_info_class->nres_nhc;
    int nyosh_nhc                = therm_info_class->nyosh_nhc;

    int mytherm_start            = therm_info_class->mytherm_start;
    int mytherm_end              = therm_info_class->mytherm_end;
    int *therm_inhc_x            = therm_info_class->inhc_x;
    int *therm_inhc_y            = therm_info_class->inhc_y;
    int *therm_inhc_z            = therm_info_class->inhc_z;
    double **therm_gkt           = therm_info_class->gkt;
    double **therm_mass_nhc      = therm_info_class->mass_nhc;
    double *therm_wdti2          = therm_info_class->wdti2;
    double *therm_wdti4          = therm_info_class->wdti4;
    double *therm_wdti8          = therm_info_class->wdti8;
    double  wght                 = therm_info_class->wght;
    double dt_nhc                = therm_info_class->dt_nhc;

    double **therm_f_nhc         = therm_class->f_nhc;
    double **therm_v_nhc         = therm_class->v_nhc;
    double **therm_x_nhc         = therm_class->x_nhc;

    double *int_scr_sc           = int_scr->sc;

    int hmat_int_typ             = cell->hmat_int_typ;


    double *baro_v_vol_nhc       = baro->v_vol_nhc;
    double *baro_f_vol_nhc       = baro->f_vol_nhc;
    double *baro_x_vol_nhc       = baro->x_vol_nhc;
    double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
    double *baro_gkt_vol         = baro->gkt_vol;

    double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
    double *par_rahman_vgmat     = par_rahman->vgmat;
    double *par_rahman_vgmat_g   = par_rahman->vgmat_g;
    double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
    double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
    double *par_rahman_vtemps    = par_rahman->vtemps;
    double *par_rahman_vtempv    = par_rahman->vtempv;
    double *par_rahman_veig      = par_rahman->veig;
    double *par_rahman_veigv     = par_rahman->veigv;
    double *par_rahman_vexpdt    = par_rahman->vexpdt;
    double par_rahman_mass_hm    = par_rahman->mass_hm;
    double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
    double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);

    len_nhcm1    = len_nhc-1;
    len_nhcp1    = len_nhc+1;

/*==========================================================================*/
/* I) Initialize the roll scaling factors                                   */

    for(i=1;i<=9;i++){par_rahman_roll_mtvv[i] = 0.0;}
    par_rahman_roll_mtvv[1]=1.0;
    par_rahman_roll_mtvv[5]=1.0;
    par_rahman_roll_mtvv[9]=1.0;

    for(ipart=1;ipart<=natm_tot;ipart++){
      clatoms_roll_sc[ipart] = 1.0;
    }/*endfor*/
    (*par_rahman_roll_scg)  = 0.0;
    (*par_rahman_roll_scg0) = 1.0;

    int_scr_sc[num_nhc+1] = 1.0;

/*==========================================================================*/
/* III) Get the force on the first NHC in each chain and on the baro stat   */

    baro_f_vol_nhc[1] = -baro_gkt_vol[1];
    for(i=1;i<=9;i++){
      baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                           par_rahman_vgmat[i]*par_rahman_vgmat[i];
    }/*endfor*/
    baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];


    get_fnhc1_par(clatoms_info,clatoms_pos,
                  therm_info_class,therm_class,
                  baro,par_rahman,class_comm_forc_pkg,int_scr);
    get_fgmatv_par(clatoms_info,clatoms_pos,par_rahman,cell,
                   class_comm_forc_pkg);
  

/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */
 
    for(iresn=1;iresn<=nres_nhc;iresn++){
      for(iyosh=1;iyosh<=nyosh_nhc;iyosh++){
/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        }/*endfor*/
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
          }/*endfor*/
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          baro_v_vol_nhc[len_nhc-ichain] = 
                     baro_v_vol_nhc[len_nhc-ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[len_nhc-ichain]*aa;
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the vgmat                                                      */

        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = 
              (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 

        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                               therm_wdti4[iyosh]*
                              (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                              /(par_rahman_mass_hm);
        }/*endfor*/

/*--------------------------------------------------------------------------*/
/* 4) Evolve the classical velocities                                       */

       if(hmat_int_typ==0){
         move_vel_vbox(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4, therm_wdti8,
                   therm_x_nhc,therm_v_nhc,therm_inhc_x,therm_inhc_y,
                   therm_inhc_z,iyosh,num_nhc, int_scr_sc, 
                   clatoms_roll_sc,wght,len_nhc,mytherm_start,mytherm_end);
       }else{
         move_vel_vbox_upper(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4, therm_wdti8,
                   therm_x_nhc,therm_v_nhc,therm_inhc_x,therm_inhc_y,
                   therm_inhc_z,iyosh,num_nhc, int_scr_sc, 
                   clatoms_roll_sc,wght,len_nhc);

       }/*endif : hmat_int_typ*/

/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          for(inhc=1;inhc<=num_nhc;inhc++){
            therm_x_nhc[ichain][inhc] += 
                       therm_v_nhc[ichain][inhc]*therm_wdti2[iyosh]*wght;
          }/*endfor*/
          baro_x_vol_nhc[ichain] += 
                       baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the vgmat                                                     */


        get_fgmatv_par(clatoms_info,clatoms_pos,par_rahman,cell,
                   class_comm_forc_pkg);

        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = (*par_rahman_roll_scg)*aa2
                                 +therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 
 
        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                                therm_wdti4[iyosh]*
                             (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                             /(par_rahman_mass_hm);
        }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */

        get_fnhc1_par(clatoms_info,clatoms_pos,
                      therm_info_class,therm_class,
                      baro,par_rahman,class_comm_forc_pkg,int_scr);

        baro_f_vol_nhc[1] = -baro_gkt_vol[1];
        for(i=1;i<=9;i++){
         baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                           par_rahman_vgmat[i]*par_rahman_vgmat[i];
        }/*endfor*/
        baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];

        for(ichain=1;ichain<=len_nhcm1;ichain++){
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            arg = -therm_wdti8[iyosh]*therm_v_nhc[ichain+1][inhc]*wght;
            aa = exp(arg);
            therm_v_nhc[ichain][inhc] = therm_v_nhc[ichain][inhc]*aa*aa
              + therm_wdti4[iyosh]*therm_f_nhc[ichain][inhc]*aa*wght;
          }/*endfor*/
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
            therm_f_nhc[ichain+1][inhc] = 
                 (therm_mass_nhc[ichain][inhc]*
                  therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]-
                  therm_gkt[ichain+1][inhc])/therm_mass_nhc[ichain+1][inhc];

          }/*endfor*/
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
          therm_v_nhc[len_nhc][inhc] += 
                       therm_f_nhc[len_nhc][inhc]*therm_wdti4[iyosh]*wght;
        }/*endfor*/
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      }/*endfor: iyosh */
    }/*endfor: iresn*/

/*==========================================================================*/
/* VI) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    for(i=1;i<=9;i++){par_rahman_vgmat_g[i]=par_rahman_vgmat[i];}
    (*par_rahman_roll_scg) *= (2.0/dt_nhc);

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPF0(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                  THERM_INFO *therm_info_class, 
                  THERM_POS *therm_class, 
                  BARO *baro, 
                  PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr,
                  double dt)

/*========================================================================*/
    {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int i,j,joff,n,ipart,ichain;   /* Num: for loop counters */
    int iresn,iyosh;               /* Num: for loop counters */
    int len_nhcm1,len_nhcp1;       /* Num: length of chains  */
    double arg,aa,aa2;             /* Num: scalar temps      */
    double tempvx,tempvy,tempvz;  

/* Define local pointers                                                */

    int natm_tot                 = clatoms_info->natm_tot;
    double *clatoms_roll_sc      = clatoms_info->roll_sc;

    double *clatoms_vx           = clatoms_pos->vx;
    double *clatoms_vy           = clatoms_pos->vy;
    double *clatoms_vz           = clatoms_pos->vz;

    double *int_scr_sc           = int_scr->sc;

    int num_nhc                  = therm_info_class->num_nhc; 
    int len_nhc                  = therm_info_class->len_nhc;
    int nres_nhc                 = therm_info_class->nres_nhc;
    int nyosh_nhc                = therm_info_class->nyosh_nhc;
    double *therm_wdti2          = therm_info_class->wdti2;
    double *therm_wdti4          = therm_info_class->wdti4;
    double *therm_wdti8          = therm_info_class->wdti8;
    double wght                  = therm_info_class->wght;
    double therm_info_class_dt_nhc = therm_info_class->dt_nhc;

    int hmat_int_typ             = cell->hmat_int_typ;

    double *baro_v_vol_nhc       = baro->v_vol_nhc;
    double *baro_f_vol_nhc       = baro->f_vol_nhc;
    double *baro_x_vol_nhc       = baro->x_vol_nhc;
    double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
    double *baro_gkt_vol         = baro->gkt_vol;

    double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
    double *par_rahman_vgmat     = par_rahman->vgmat;
    double *par_rahman_vgmat_g   = par_rahman->vgmat_g;
    double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
    double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
    double *par_rahman_vtemps    = par_rahman->vtemps;
    double *par_rahman_vtempv    = par_rahman->vtempv;
    double *par_rahman_veig      = par_rahman->veig;
    double *par_rahman_veigv     = par_rahman->veigv;
    double *par_rahman_vexpdt    = par_rahman->vexpdt;
    double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
    double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);  
    double c1_hm                 = par_rahman->c1_hm;
    double par_rahman_mass_hm    = par_rahman->mass_hm;

/*==========================================================================*/
/* I) Initialize the roll scaling factors                                   */

    for(i=1;i<=9;i++){par_rahman_roll_mtvv[i] = 0.0;}
    par_rahman_roll_mtvv[1]=1.0;par_rahman_roll_mtvv[5]=1.0;
    par_rahman_roll_mtvv[9]=1.0;
    for(ipart=1;ipart<=natm_tot;ipart++){
       clatoms_roll_sc[ipart] = 1.0;
    }/*endfor*/
    (*par_rahman_roll_scg)  = 0.0;
    (*par_rahman_roll_scg0) = 1.0;
    int_scr_sc[num_nhc+1] = 1.0;

/*==========================================================================*/
/* III) Get first NHC on the baro stat                                      */

    baro_f_vol_nhc[1] = -baro_gkt_vol[1];
    for(i=1;i<=9;i++){
      baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                           par_rahman_vgmat[i]*par_rahman_vgmat[i];
    }/*endfor*/
    baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhcm1 = (len_nhc)-1;
    len_nhcp1 = (len_nhc)+1;
    for(iresn=1;iresn<=nres_nhc;iresn++){
      for(iyosh=1;iyosh<=nyosh_nhc;iyosh++){

/*--------------------------------------------------------------------------*/
/*  1) Evolve the last therm velocity in each chain                         */
        baro_v_vol_nhc[len_nhc] += 
                            baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/*  2) Evolve the last-1 to the first thermo velocity in each chain         */
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[len_nhcp1-ichain];
          aa = exp(arg);
          baro_v_vol_nhc[len_nhc-ichain] = 
                     baro_v_vol_nhc[len_nhc-ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[len_nhc-ichain]*aa;
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the vgmat,classical velocities, and therm positions            */

       if(hmat_int_typ==0){
         move_vel_vbox0(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4,therm_wdti8,iyosh,len_nhc);
       }else{
         move_vel_vbox_upper0(clatoms_pos,clatoms_info,cell,par_rahman, 
                   baro_x_vol_nhc,baro_v_vol_nhc,
                   therm_wdti2,therm_wdti4,therm_wdti8,iyosh,len_nhc,dt);
       }/*endif : hmat_int_typ*/

/*--------------------------------------------------------------------------*/
/*  4) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        baro_f_vol_nhc[1] = -baro_gkt_vol[1];
        for(i=1;i<=9;i++){
          baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                               par_rahman_vgmat[i]*par_rahman_vgmat[i];
        }/*endfor*/
        baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  5) Evolve the last therm velocotiy in each chain                        */
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 6) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* 7) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    for(i=1;i<=9;i++){par_rahman_vgmat_g[i]=par_rahman_vgmat[i];}
    (*par_rahman_roll_scg) *= (2.0/(therm_info_class_dt_nhc));

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void apply_NHCPF0_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                  THERM_INFO *therm_info_class, 
                  THERM_POS *therm_class, 
                  BARO *baro, 
                  PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr,
                  CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
        {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int i,j,joff,n,ipart,ichain;     /* Num: for loop counters */
    int iresn,iyosh;                 /* Num: for loop counters */
    int len_nhcm1,len_nhcp1;         /* Num: length of chains  */
    double arg,aa,aa2;               /* Num: scalar temps      */
    double trace3;
    double tempvx,tempvy,tempvz;  

/* Define local pointers                                                */

    int natm_tot                 = clatoms_info->natm_tot;
    int myatm_start              = clatoms_info->myatm_start;
    int myatm_end                = clatoms_info->myatm_end;
    double *clatoms_roll_sc      = clatoms_info->roll_sc;

    double *int_scr_sc           = int_scr->sc;
 
    double *clatoms_vx           = clatoms_pos->vx;
    double *clatoms_vy           = clatoms_pos->vy;
    double *clatoms_vz           = clatoms_pos->vz;

    int num_nhc                  = therm_info_class->num_nhc; 
    int len_nhc                  = therm_info_class->len_nhc;
    int nres_nhc                 = therm_info_class->nres_nhc;
    int nyosh_nhc                = therm_info_class->nyosh_nhc;
    double *therm_wdti2          = therm_info_class->wdti2;
    double *therm_wdti4          = therm_info_class->wdti4;
    double *therm_wdti8          = therm_info_class->wdti8;
    double wght                  = therm_info_class->wght;
    double therm_info_class_dt_nhc = therm_info_class->dt_nhc;

    double *baro_v_vol_nhc       = baro->v_vol_nhc;
    double *baro_f_vol_nhc       = baro->f_vol_nhc;
    double *baro_x_vol_nhc       = baro->x_vol_nhc;
    double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
    double *baro_gkt_vol         = baro->gkt_vol;

    double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
    double *par_rahman_vgmat     = par_rahman->vgmat;
    double *par_rahman_vgmat_g   = par_rahman->vgmat_g;
    double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
    double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
    double *par_rahman_vtemps    = par_rahman->vtemps;
    double *par_rahman_vtempv    = par_rahman->vtempv;
    double *par_rahman_veig      = par_rahman->veig;
    double *par_rahman_veigv     = par_rahman->veigv;
    double *par_rahman_vexpdt    = par_rahman->vexpdt;
    double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
    double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);  
    double par_rahman_mass_hm    = par_rahman->mass_hm;
    double c1_hm                 = par_rahman->c1_hm;
  
/*==========================================================================*/
/* I) Initialize the roll scaling factors                                   */

    for(i=1;i<=9;i++){par_rahman_roll_mtvv[i] = 0.0;}
    par_rahman_roll_mtvv[1]=1.0;par_rahman_roll_mtvv[5]=1.0;
    par_rahman_roll_mtvv[9]=1.0;
    for(ipart=1;ipart<=natm_tot;ipart++){
       clatoms_roll_sc[ipart] = 1.0;
    }/*endfor*/
    (*par_rahman_roll_scg)  = 0.0;
    (*par_rahman_roll_scg0) = 1.0;
    int_scr_sc[num_nhc+1]   = 1.0;

/*==========================================================================*/
/* III) Get first NHC on the baro stat                                      */

    baro_f_vol_nhc[1] = -baro_gkt_vol[1];
    for(i=1;i<=9;i++){
      baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                          par_rahman_vgmat[i]*par_rahman_vgmat[i];
    }/*endfor*/
    baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];
  
/*==========================================================================*/
/* IV) Apply the nhc evolution operator using RESPA                         */

    len_nhcm1 = (len_nhc)-1;
    len_nhcp1 = (len_nhc)+1;
    for(iresn=1;iresn<=nres_nhc;iresn++){
      for(iyosh=1;iyosh<=nyosh_nhc;iyosh++){

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
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/* 3) Evolve the vgmat                                                      */
        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = (*par_rahman_roll_scg)*aa2
                                  +therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 
        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                               therm_wdti4[iyosh]*
                             (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                             /(par_rahman_mass_hm);
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  4.2) Evolve the particle velocities with vgmat                          */
        for(i=1;i<=9;i++){par_rahman_vtemps[i] = par_rahman_vgmat[i];}
        trace3 = c1_hm*(par_rahman_vgmat[1]+par_rahman_vgmat[5]
                       +par_rahman_vgmat[9]);
        par_rahman_vtemps[1] += trace3;    
        par_rahman_vtemps[5] += trace3;
        par_rahman_vtemps[9] += trace3;
        diag33(par_rahman->vtemps,par_rahman->veig,
               par_rahman->veigv,par_rahman->fv1,
               par_rahman->fv2);
        for(i=1;i<=3;i++){
          par_rahman_vexpdt[i] = exp(-par_rahman_veig[i]*therm_wdti2[iyosh]);
        }/*endfor*/
        for(i=1;i<=3;i++){
          joff = (i-1)*3 ;
          for(j=1;j<=3;j++){ 
            (par_rahman_vtemps)[j+joff] = (par_rahman_veigv)[j+joff]*
                                           (par_rahman_vexpdt)[i];
          }/*endfor*/
        }/*endfor*/
        n = 3;
        matmul_t2(par_rahman_veigv,par_rahman_vtemps,
                  par_rahman_vtempv,n);
        matmul_t2(par_rahman_vtempv,par_rahman_roll_mtvv,
                  par_rahman_vtemps,n);
        for(i=1;i<=9;i++){par_rahman_roll_mtvv[i]=par_rahman_vtemps[i];}
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          tempvx =  clatoms_vx[ipart]*par_rahman_vtempv[1]
                 +  clatoms_vy[ipart]*par_rahman_vtempv[2]
                 +  clatoms_vz[ipart]*par_rahman_vtempv[3];
          tempvy =  clatoms_vx[ipart]*par_rahman_vtempv[4]
                 +  clatoms_vy[ipart]*par_rahman_vtempv[5]
                 +  clatoms_vz[ipart]*par_rahman_vtempv[6];
          tempvz =  clatoms_vx[ipart]*par_rahman_vtempv[7]
                 +  clatoms_vy[ipart]*par_rahman_vtempv[8]
                 +  clatoms_vz[ipart]*par_rahman_vtempv[9];
          clatoms_vx[ipart] = tempvx;
          clatoms_vy[ipart] = tempvy;
          clatoms_vz[ipart] = tempvz;
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */
        for(ichain=1;ichain<=len_nhc;ichain++){
          baro_x_vol_nhc[ichain] += 
                       baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  6) Evolve the vgmat                                                     */
        get_fgmatv_par(clatoms_info,clatoms_pos,par_rahman,cell,
                       class_comm_forc_pkg);
        arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
        aa = exp(arg);
        aa2 = aa*aa;
        (*par_rahman_roll_scg) = (*par_rahman_roll_scg)*aa2
                                  +therm_wdti4[iyosh]*aa;
        (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 
        for(i=1;i<=9;i++){       
          par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                             therm_wdti4[iyosh]*
                             (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                             /(par_rahman_mass_hm);
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  7) Evolve the 1 to last-1 therm velocity in each chain                  */
/*     calculting therm forces as you go along                              */
        baro_f_vol_nhc[1] = -baro_gkt_vol[1];
        for(i=1;i<=9;i++){
          baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                               par_rahman_vgmat[i]*par_rahman_vgmat[i];
        }/*endfor*/
        baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];
        for(ichain=1;ichain<=len_nhcm1;ichain++){
          arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[ichain+1];
          aa = exp(arg);
          baro_v_vol_nhc[ichain] = baro_v_vol_nhc[ichain]*aa*aa
                   + therm_wdti4[iyosh]*baro_f_vol_nhc[ichain]*aa;
          baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
        }/*endfor*/
/*--------------------------------------------------------------------------*/
/*  8) Evolve the last therm velocotiy in each chain                        */
        baro_v_vol_nhc[len_nhc] += 
                baro_f_vol_nhc[len_nhc]*therm_wdti4[iyosh];
/*--------------------------------------------------------------------------*/
/* 9) End Respa                                                             */

      /*endfor: iyosh */}
    /*endfor: iresn*/}

/*==========================================================================*/
/* VI) Save the present value of the barostat velocity                      */
/*     Scale barostat roll_scg variable                                     */

    for(i=1;i<=9;i++){par_rahman_vgmat_g[i]=par_rahman_vgmat[i];}
    (*par_rahman_roll_scg) *= (2.0/(therm_info_class_dt_nhc));

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NHCPF(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                BARO *baro, 
                PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int inhc,ichain;       /* Num: for loop counters */
   int len_nhcm1;        /* Num: length of chains -1*/

   int num_nhc = therm_info_class->num_nhc; 
   int len_nhc = therm_info_class->len_nhc;

   double **therm_class_f_nhc         = therm_class->f_nhc;
   double **therm_info_class_mass_nhc = therm_info_class->mass_nhc;
   double **therm_class_v_nhc         = therm_class->v_nhc;
   double **therm_info_class_gkt      = therm_info_class->gkt;

   double *baro_f_vol_nhc             = baro->f_vol_nhc;
   double *baro_mass_vol_nhc          = baro->mass_vol_nhc;
   double *baro_v_vol_nhc             = baro->v_vol_nhc;
   double *baro_gkt_vol               = baro->gkt_vol;

/*==========================================================================*/
/* I) Get the force on the first NHC in each chain and on the baro stat     */

    get_fnhc1(clatoms_info,clatoms_pos,therm_info_class,therm_class,
              baro,par_rahman);
    get_fgmatv(clatoms_info,clatoms_pos,par_rahman,cell);

/*==========================================================================*/
/* II) Get the force on the rest of the NHC in each chain                   */

    len_nhcm1 = len_nhc-1;        
    for(ichain=1;ichain<=len_nhcm1;ichain++){
      for(inhc=1;inhc<=num_nhc;inhc++){
        therm_class_f_nhc[ichain+1][inhc] = 
                     (therm_info_class_mass_nhc[ichain][inhc]*
                     therm_class_v_nhc[ichain][inhc]*
                     therm_class_v_nhc[ichain][inhc]-
                     therm_info_class_gkt[ichain+1][inhc])/
                     therm_info_class_mass_nhc[ichain+1][inhc];
      }/*endfor*/
      baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
    }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_NHCPF_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                THERM_INFO *therm_info_class,THERM_POS *therm_class, 
                BARO *baro, 
                PAR_RAHMAN *par_rahman, CELL *cell, INT_SCR *int_scr,
                CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
     {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int inhc,ichain;       /* Num: for loop counters */
   int len_nhcm1;        /* Num: length of chains -1*/
   int i;

/*-----------------------------------------------------------------------*/
/*             Local pointer declarations                                */

   int mytherm_start             = therm_info_class->mytherm_start;
   int mytherm_end               = therm_info_class->mytherm_end;
   int num_nhc                   = therm_info_class->num_nhc; 
   int len_nhc                   = therm_info_class->len_nhc;
   double **therm_info_class_gkt      = therm_info_class->gkt;
   double **therm_info_class_mass_nhc = therm_info_class->mass_nhc;

   double *baro_f_vol_nhc        = baro->f_vol_nhc;
   double *baro_mass_vol_nhc     = baro->mass_vol_nhc;
   double *baro_v_vol_nhc        = baro->v_vol_nhc;
   double *baro_gkt_vol          = baro->gkt_vol;

   double par_rahman_mass_hm     = par_rahman->mass_hm;
   double *par_rahman_vgmat      = par_rahman->vgmat;

   double **therm_class_f_nhc    = therm_class->f_nhc;
   double **therm_class_v_nhc    = therm_class->v_nhc;

/*==========================================================================*/
/* I) Get the force on the first NHC in each chain and on the baro stat     */

   baro_f_vol_nhc[1] = -baro_gkt_vol[1];
   for(i=1;i<=9;i++){
     baro_f_vol_nhc[1]+= par_rahman_mass_hm*
                         par_rahman_vgmat[i]*par_rahman_vgmat[i];
   }/*endfor*/
   baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];

   get_fnhc1_par(clatoms_info,clatoms_pos,therm_info_class,therm_class,
                 baro,par_rahman,class_comm_forc_pkg,int_scr);
   get_fgmatv_par(clatoms_info,clatoms_pos,par_rahman,cell,
                  class_comm_forc_pkg);

/*==========================================================================*/
/* II) Get the force on the rest of the NHC in each chain                   */

   len_nhcm1 = len_nhc-1;        
   for(ichain=1;ichain<=len_nhcm1;ichain++){
     for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        therm_class_f_nhc[ichain+1][inhc] = 
                     (therm_info_class_mass_nhc[ichain][inhc]*
                     therm_class_v_nhc[ichain][inhc]*
                     therm_class_v_nhc[ichain][inhc]-
                     therm_info_class_gkt[ichain+1][inhc])/
                     therm_info_class_mass_nhc[ichain+1][inhc];
     }/*endfor*/
     baro_f_vol_nhc[ichain+1] = 
                   (baro_mass_vol_nhc[ichain]
                   *baro_v_vol_nhc[ichain]*baro_v_vol_nhc[ichain]
                   -baro_gkt_vol[ichain+1])/baro_mass_vol_nhc[ichain+1];
   }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_fnhc1(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               THERM_INFO *therm_info_class,THERM_POS *therm_class, 
               BARO *baro, 
               PAR_RAHMAN *par_rahman)

/*========================================================================*/
     {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int ipart,inhc,i;

/* Define local pointers                                                  */

   int num_nhc                  = therm_info_class->num_nhc;
   int *therm_inhc_x            = therm_info_class->inhc_x;
   int *therm_inhc_y            = therm_info_class->inhc_y;
   int *therm_inhc_z            = therm_info_class->inhc_z;
   double **therm_gkt           = therm_info_class->gkt;
   double **therm_mass_nhc      = therm_info_class->mass_nhc;

   double **therm_f_nhc         = therm_class->f_nhc;

   int natm_tot                 = clatoms_info->natm_tot;
   double *clatoms_mass         = clatoms_info->mass;

   double *clatoms_vx           = clatoms_pos->vx;
   double *clatoms_vy           = clatoms_pos->vy;
   double *clatoms_vz           = clatoms_pos->vz;

   double *par_rahman_vgmat     = par_rahman->vgmat;

   double *baro_f_vol_nhc       = baro->f_vol_nhc;
   double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
   double *baro_gkt_vol         = baro->gkt_vol;

/*========================================================================*/
/* I) f_nhc[1][inhc]                                                      */

    for(inhc=1;inhc<=num_nhc;inhc++){
      therm_f_nhc[1][inhc] = -therm_gkt[1][inhc];
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      therm_f_nhc[1][therm_inhc_x[ipart]] +=
                clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      therm_f_nhc[1][therm_inhc_y[ipart]] +=
                clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
    }/*endfor*/
    for(ipart=1;ipart<=natm_tot;ipart++){
      therm_f_nhc[1][therm_inhc_z[ipart]] +=
                clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
    }/*endfor*/
    for(inhc=1;inhc<=num_nhc;inhc++){
      therm_f_nhc[1][inhc] /= therm_mass_nhc[1][inhc];
    }/*endfor*/

/*========================================================================*/
/* III) f_vol_nhc[1][inhc]                                                  */

    baro_f_vol_nhc[1] = -baro_gkt_vol[1];
    for(i=1;i<=9;i++){
      baro_f_vol_nhc[1]+= par_rahman->mass_hm*
                           par_rahman_vgmat[i]*par_rahman_vgmat[i];
    }/*endfor*/
    baro_f_vol_nhc[1] /= baro_mass_vol_nhc[1];

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_fnhc1_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
               THERM_INFO *therm_info_class,THERM_POS *therm_class, 
               BARO *baro,PAR_RAHMAN *par_rahman,
               CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
               INT_SCR *int_scr)

/*========================================================================*/
      {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

int ipart,inhc,i;

/* Define local pointers                                                  */

   int num_nhc                  = therm_info_class->num_nhc;
   int *therm_inhc_x            = therm_info_class->inhc_x;
   int *therm_inhc_y            = therm_info_class->inhc_y;
   int *therm_inhc_z            = therm_info_class->inhc_z;
   int mytherm_start            = therm_info_class->mytherm_start;
   int mytherm_end              = therm_info_class->mytherm_end;
   double **therm_gkt           = therm_info_class->gkt;
   double **therm_mass_nhc      = therm_info_class->mass_nhc;
   int num_nhc_share            = therm_info_class->num_nhc_share;
   int *map_share               = therm_info_class->map_share;

   double **therm_f_nhc         = therm_class->f_nhc;

   int natm_tot                 = clatoms_info->natm_tot;
   int myatm_start              = clatoms_info->myatm_start;
   int myatm_end                = clatoms_info->myatm_end;
   double *clatoms_mass         = clatoms_info->mass;

   double *clatoms_vx           = clatoms_pos->vx;
   double *clatoms_vy           = clatoms_pos->vy;
   double *clatoms_vz           = clatoms_pos->vz;

   double *par_rahman_vgmat     = par_rahman->vgmat;

   double *baro_f_vol_nhc       = baro->f_vol_nhc;
   double *baro_mass_vol_nhc    = baro->mass_vol_nhc;
   double *baro_gkt_vol         = baro->gkt_vol;

   int np_forc                  = class_comm_forc_pkg->num_proc;
   MPI_Comm comm_forc           = class_comm_forc_pkg->comm;
 
   double *int_scr_atm_kin      = int_scr->atm_kin;
   double *int_scr_sc           = int_scr->sc;
   double *int_scr_sc_temp      = int_scr->sc_temp;

/*========================================================================*/
/* I) f_nhc[1][inhc]                                                      */


  for(inhc=1;inhc<=num_nhc+1;inhc++){
    int_scr_atm_kin[inhc] = 0.0;
    int_scr_sc[inhc]       = 0.0;
  }/*endfor*/

  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    int_scr_atm_kin[therm_inhc_x[ipart]] += 
                 clatoms_mass[ipart]*clatoms_vx[ipart]*clatoms_vx[ipart];
  }/*endfor*/
  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    int_scr_atm_kin[therm_inhc_y[ipart]] += 
                clatoms_mass[ipart]*clatoms_vy[ipart]*clatoms_vy[ipart];
  }/*endfor*/
  for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    int_scr_atm_kin[therm_inhc_z[ipart]] += 
               clatoms_mass[ipart]*clatoms_vz[ipart]*clatoms_vz[ipart];
  }/*endfor*/

  if(np_forc > 1){

     for(inhc = 1;inhc<=num_nhc_share;inhc++){
       int_scr_sc[inhc] = int_scr_atm_kin[map_share[inhc]];
     }/*endfor*/
     Allreduce(&(int_scr_sc[1]), &(int_scr_sc_temp[1]),num_nhc_share,
                 MPI_DOUBLE,MPI_SUM,0,comm_forc);
     for(inhc = 1;inhc<=num_nhc_share;inhc++){
       int_scr_atm_kin[map_share[inhc]] = int_scr_sc_temp[inhc];
     }/*endfor*/

  }/*endif*/
  
  for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
     therm_f_nhc[1][inhc] = (int_scr_atm_kin[inhc] - therm_gkt[1][inhc])/
                             therm_mass_nhc[1][inhc];
  }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_fgmatv(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                PAR_RAHMAN *par_rahman, CELL *cell)

/*========================================================================*/
       {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

   int i,ipart,iii;
   double trace3;

/* Define local pointers                                                  */

   int natm_tot                 = clatoms_info->natm_tot;
   double *clatoms_mass         = clatoms_info->mass;

   double *clatoms_vx           = clatoms_pos->vx;
   double *clatoms_vy           = clatoms_pos->vy;
   double *clatoms_vz           = clatoms_pos->vz;

   double *par_rahman_fgmat_v   = par_rahman->fgmat_v;

   int iperd                    = cell->iperd;
   int hmat_cons_typ            = cell->hmat_cons_typ;
   int hmat_int_typ             = cell->hmat_int_typ;

/*========================================================================*/
/* II) fgmat_v                                                            */

   for(i=1;i<=9;i++){par_rahman_fgmat_v[i] = 0.0;}
   for(ipart=1;ipart<=natm_tot;ipart++){
     par_rahman_fgmat_v[1] += clatoms_mass[ipart]*
                              clatoms_vx[ipart]*clatoms_vx[ipart];
     par_rahman_fgmat_v[5] += clatoms_mass[ipart]*
                              clatoms_vy[ipart]*clatoms_vy[ipart];
     par_rahman_fgmat_v[9] += clatoms_mass[ipart]*
                              clatoms_vz[ipart]*clatoms_vz[ipart];
     par_rahman_fgmat_v[2] += clatoms_mass[ipart]*
                              clatoms_vx[ipart]*clatoms_vy[ipart];
   }/*endfor*/
   trace3 = (par_rahman_fgmat_v[1]+par_rahman_fgmat_v[5]
            +par_rahman_fgmat_v[9]);
   par_rahman_fgmat_v[4] = par_rahman_fgmat_v[2];
   par_rahman_fgmat_v[5] += trace3;
   par_rahman_fgmat_v[9] += trace3;

   if(iperd==3){
     for(ipart=1;ipart<=natm_tot;ipart++){
       par_rahman_fgmat_v[3] += clatoms_mass[ipart]*
                                clatoms_vx[ipart]*clatoms_vz[ipart];
       par_rahman_fgmat_v[6] += clatoms_mass[ipart]*
                                clatoms_vy[ipart]*clatoms_vz[ipart];
     }/*endfor*/
     par_rahman_fgmat_v[7] = par_rahman_fgmat_v[3];
     par_rahman_fgmat_v[8] = par_rahman_fgmat_v[6];
   }/*endif*/

   constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,par_rahman_fgmat_v);

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void get_fgmatv_par(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                PAR_RAHMAN *par_rahman, CELL *cell,
                CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*========================================================================*/
     {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   int i,ipart;
   double trace3;

/* Define local pointers                                                  */

   int myatm_start                 = clatoms_info->myatm_start;
   int myatm_end                   = clatoms_info->myatm_end;
   double *clatoms_mass            = clatoms_info->mass;

   double *clatoms_vx              = clatoms_pos->vx;
   double *clatoms_vy              = clatoms_pos->vy;
   double *clatoms_vz              = clatoms_pos->vz;

   int np_forc                     = class_comm_forc_pkg->num_proc;
   MPI_Comm comm_forc              = class_comm_forc_pkg->comm;

   double *par_rahman_fgmat_v      = par_rahman->fgmat_v;
   double *par_rahman_fgmat_v_temp = par_rahman->vtemps;
   double c1_hm                    = par_rahman->c1_hm;

   int iperd                       = cell->iperd;
   int hmat_cons_typ               = cell->hmat_cons_typ;
   int hmat_int_typ                = cell->hmat_int_typ;

/*========================================================================*/
/* II) fgmat_v                                                            */

   for(i=1;i<=9;i++){par_rahman_fgmat_v[i] = 0.0;}
   for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
     par_rahman_fgmat_v[1] += clatoms_mass[ipart]*
                              clatoms_vx[ipart]*clatoms_vx[ipart];
     par_rahman_fgmat_v[5] += clatoms_mass[ipart]*
                              clatoms_vy[ipart]*clatoms_vy[ipart];
     par_rahman_fgmat_v[9] += clatoms_mass[ipart]*
                              clatoms_vz[ipart]*clatoms_vz[ipart];
     par_rahman_fgmat_v[2] += clatoms_mass[ipart]*
                              clatoms_vx[ipart]*clatoms_vy[ipart];
    }/*endfor*/
    trace3 = c1_hm*(par_rahman_fgmat_v[1]+par_rahman_fgmat_v[5]
             +par_rahman_fgmat_v[9]);
    par_rahman_fgmat_v[1] += trace3;
    par_rahman_fgmat_v[5] += trace3;
    par_rahman_fgmat_v[9] += trace3;
    par_rahman_fgmat_v[4] = par_rahman_fgmat_v[2];

    if(iperd==3){
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        par_rahman_fgmat_v[3] += clatoms_mass[ipart]*
                                 clatoms_vx[ipart]*clatoms_vz[ipart];
        par_rahman_fgmat_v[6] += clatoms_mass[ipart]*
                                 clatoms_vy[ipart]*clatoms_vz[ipart];
      }/*endfor*/
      par_rahman_fgmat_v[7] = par_rahman_fgmat_v[3];
      par_rahman_fgmat_v[8] = par_rahman_fgmat_v[6];
    }/*endif*/

    constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,par_rahman_fgmat_v);

    if(np_forc > 1){

      for(i=1;i<=9;i++){
        par_rahman_fgmat_v_temp[i] = par_rahman_fgmat_v[i];
      }/*endfor*/    

      Allreduce(&(par_rahman_fgmat_v_temp[1]),&(par_rahman_fgmat_v[1]),9,
                  MPI_DOUBLE,MPI_SUM,0,comm_forc);
    }/*endif*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_vel_vbox0(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                    CELL *cell,PAR_RAHMAN *par_rahman, 
                    double *baro_x_vol_nhc,double *baro_v_vol_nhc,
                    double *therm_wdti2,double *therm_wdti4,
                    double *therm_wdti8,int iyosh,int len_nhc)

/*==========================================================================*/
      {/*Begin Routine*/
/*==========================================================================*/
/*             Local variable declarations                                  */

   int i,j,joff,ipart,n,ichain;
   double arg,aa,aa2,trace3;
   double tempvx,tempvy,tempvz;

   double c1_hm                 = par_rahman->c1_hm;
   double *par_rahman_vgmat     = par_rahman->vgmat;
   double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
   double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
   double *par_rahman_vtemps    = par_rahman->vtemps;
   double *par_rahman_vtempv    = par_rahman->vtempv;
   double *par_rahman_vexpdt    = par_rahman->vexpdt;
   double *par_rahman_veig      = par_rahman->veig;
   double *par_rahman_veigv     = par_rahman->veigv;
   double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
   double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
   double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);
   double par_rahman_mass_hm    = par_rahman->mass_hm;


   double *clatoms_vx           = clatoms_pos->vx;
   double *clatoms_vy           = clatoms_pos->vy;
   double *clatoms_vz           = clatoms_pos->vz;
   int natm_tot                 = clatoms_info->natm_tot;
   int myatm_start              = clatoms_info->myatm_start;
   int myatm_end                = clatoms_info->myatm_end;

/*--------------------------------------------------------------------------*/
/* 3) Evolve the vgmat                                                      */

   arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
   aa = exp(arg);
   aa2 = aa*aa;
   (*par_rahman_roll_scg)  = (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
   (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 

   for(i=1;i<=9;i++){       
     par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                           therm_wdti4[iyosh]*
                           (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                          /(par_rahman_mass_hm);
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  4.2) Evolve the particle velocities with vgmat                          */

   for(i=1;i<=9;i++){par_rahman_vtemps[i] = par_rahman_vgmat[i];}
   trace3 = c1_hm*(par_rahman_vgmat[1]+par_rahman_vgmat[5]
                  +par_rahman_vgmat[9]);
   par_rahman_vtemps[1] += trace3; 
   par_rahman_vtemps[5] += trace3;
   par_rahman_vtemps[9] += trace3;

   diag33(par_rahman->vtemps,par_rahman->veig,
          par_rahman->veigv,par_rahman->fv1,par_rahman->fv2);
   for(i=1;i<=3;i++){
     par_rahman_vexpdt[i] = exp(-par_rahman_veig[i]*therm_wdti2[iyosh]);
   }/*endfor*/

   for(i=1;i<=3;i++){
     joff = (i-1)*3 ;
     for(j=1;j<=3;j++){ 
       par_rahman_vtemps[j+joff] = par_rahman_veigv[j+joff]*
                                   par_rahman_vexpdt[i];
     }/*endfor*/
   }/*endfor*/

   n = 3;
   matmul_t2(par_rahman_veigv,par_rahman_vtemps,par_rahman_vtempv,n);
   matmul_t2(par_rahman_vtempv,par_rahman_roll_mtvv,par_rahman_vtemps,n);
   for(i=1;i<=9;i++){par_rahman_roll_mtvv[i]=par_rahman_vtemps[i];}

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     tempvx =  clatoms_vx[ipart]*par_rahman_vtempv[1]
            +  clatoms_vy[ipart]*par_rahman_vtempv[2]
            +  clatoms_vz[ipart]*par_rahman_vtempv[3];
     tempvy =  clatoms_vx[ipart]*par_rahman_vtempv[4]
            +  clatoms_vy[ipart]*par_rahman_vtempv[5]
            +  clatoms_vz[ipart]*par_rahman_vtempv[6];
     tempvz =  clatoms_vx[ipart]*par_rahman_vtempv[7]
            +  clatoms_vy[ipart]*par_rahman_vtempv[8]
            +  clatoms_vz[ipart]*par_rahman_vtempv[9];
     clatoms_vx[ipart] = tempvx;
     clatoms_vy[ipart] = tempvy;
     clatoms_vz[ipart] = tempvz;
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  5) Evolve the therm positions                                           */

   for(ichain=1;ichain<=len_nhc;ichain++){
     baro_x_vol_nhc[ichain] += 
                              baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  6) Evolve the vgmat                                                     */

   get_fgmatv(clatoms_info,clatoms_pos,par_rahman,cell);
   arg  = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
   aa   = exp(arg);
   aa2  = aa*aa;
   (*par_rahman_roll_scg)  = (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
   (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 

   for(i=1;i<=9;i++){       
     par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                           therm_wdti4[iyosh]*
                           (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                           /(par_rahman_mass_hm);
   }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_vel_vbox_upper0(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                          CELL *cell,PAR_RAHMAN *par_rahman, 
                          double *baro_x_vol_nhc,double *baro_v_vol_nhc,
                          double *therm_wdti2,double *therm_wdti4,
                          double *therm_wdti8,int iyosh,int len_nhc, double dt)
/*==========================================================================*/
    {/*Begin Routine*/
/*==========================================================================*/
/*             Local variable declarations                                  */

   int i,ipart,iii,ichain;
   double arg,aa,aa2,temp;
   double ax,ay,az;
   double dt2,dt24;
   double gammax,gammay,gammaz;
   double vg12_vy0,vg13_vy0,vg13_vz0,vg23_vz0,tr_vgmat;

   double c1_hm                 = par_rahman->c1_hm;
   double *par_rahman_vgmat     = par_rahman->vgmat;
   double *par_rahman_fgmat_p   = par_rahman->fgmat_p;
   double *par_rahman_fgmat_v   = par_rahman->fgmat_v;
   double *par_rahman_vtemps    = par_rahman->vtemps;
   double *par_rahman_vtempv    = par_rahman->vtempv;
   double *par_rahman_vexpdt    = par_rahman->vexpdt;
   double *par_rahman_veig      = par_rahman->veig;
   double *par_rahman_veigv     = par_rahman->veigv;
   double *par_rahman_roll_mtvv = par_rahman->roll_mtvv;
   double *par_rahman_roll_scg  = &(par_rahman->roll_scg);
   double *par_rahman_roll_scg0 = &(par_rahman->roll_scg0);
   double par_rahman_mass_hm    = par_rahman->mass_hm;

   double *clatoms_vx           = clatoms_pos->vx;
   double *clatoms_vy           = clatoms_pos->vy;
   double *clatoms_vz           = clatoms_pos->vz;
   double *clatoms_mass         = clatoms_pos->mass;
   int natm_tot                 = clatoms_info->natm_tot;
   int myatm_start              = clatoms_info->myatm_start;
   int myatm_end                = clatoms_info->myatm_end;

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

   dt2 = dt/2.0;      
   dt24 = dt2*dt2;

/*--------------------------------------------------------------------------*/
/* 1) Evolve the vgmat                                                      */

   arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
   aa = exp(arg);
   aa2 = aa*aa;
   (*par_rahman_roll_scg)  = (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
   (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 

   for(i=1;i<=9;i++){       
     par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                           therm_wdti4[iyosh]*
                           (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                           /(par_rahman_mass_hm);
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  2) Evolve the particle velocities                                       */

   tr_vgmat = (par_rahman_vgmat[1]+par_rahman_vgmat[5]+
               par_rahman_vgmat[9])*c1_hm;       

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){

     gammax = -(tr_vgmat+par_rahman_vgmat[1]);
     gammay = -(tr_vgmat+par_rahman_vgmat[5]);
     gammaz = -(tr_vgmat+par_rahman_vgmat[9]);

     ax     =  exp(gammax*therm_wdti2[iyosh]);
     ay     =  exp(gammay*therm_wdti2[iyosh]);
     az     =  exp(gammaz*therm_wdti2[iyosh]);

     vg12_vy0 = par_rahman_vgmat[4]*clatoms_vy[ipart];
     vg13_vy0 = par_rahman_vgmat[7]*clatoms_vy[ipart];
     vg13_vz0 = par_rahman_vgmat[7]*clatoms_vz[ipart];
     vg23_vz0 = par_rahman_vgmat[8]*clatoms_vz[ipart];

     clatoms_vx[ipart] = clatoms_vx[ipart]*ax-vg13_vz0*dt2*(ax+az)
                       - vg12_vy0*dt2*(ax+ay) 
                       + vg23_vz0*par_rahman_vgmat[4]*dt24*(0.5*(ax+az) + ay);
     clatoms_vy[ipart] = clatoms_vy[ipart]*ay-vg23_vz0*dt2*(az+ay);
     clatoms_vz[ipart] = clatoms_vz[ipart]*az;

   }/*endfor : ipart*/

/*--------------------------------------------------------------------------*/
/*  3) Evolve the therm positions                                           */

   for(ichain=1;ichain<=len_nhc;ichain++){
     baro_x_vol_nhc[ichain] += baro_v_vol_nhc[ichain]*therm_wdti2[iyosh];
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  4) Evolve the vgmat                                                     */

   get_fgmatv(clatoms_info,clatoms_pos,par_rahman,cell);
   arg = -therm_wdti8[iyosh]*baro_v_vol_nhc[1];
   aa  = exp(arg);
   aa2 = aa*aa;
   (*par_rahman_roll_scg)  = (*par_rahman_roll_scg)*aa2+therm_wdti4[iyosh]*aa;
   (*par_rahman_roll_scg0) = (*par_rahman_roll_scg0)*aa2; 
   for(i=1;i<=9;i++){       
     par_rahman_vgmat[i] = par_rahman_vgmat[i]*aa2 + 
                           therm_wdti4[iyosh]*
                           (par_rahman_fgmat_v[i]+par_rahman_fgmat_p[i])*aa
                           /(par_rahman_mass_hm);
   }/*endfor*/

/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/


