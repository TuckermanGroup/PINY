/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: move_pos_box                                 */
/*                                                                          */
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

void move_vel_vbox(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                   CELL *cell,PAR_RAHMAN *par_rahman, 
                   double *baro_x_vol_nhc,double *baro_v_vol_nhc,
                   double *therm_wdti2,double *therm_wdti4,
                   double *therm_wdti8,
                   double **therm_x_nhc,double **therm_v_nhc,
                   int *therm_inhc_x,int *therm_inhc_y,
                   int *therm_inhc_z,int iyosh, 
                   int num_nhc, double *int_scr_sc, double *clatoms_roll_sc,
                   double wght, int len_nhc, int mytherm_start,
                   int mytherm_end)

/*==========================================================================*/
    {/*Begin Routine*/
/*==========================================================================*/
/*             Local variable declarations                                  */

   int i,j,joff,inhc,ipart,ktemp,ichain,n;
   double arg,aa,aa2,trace3;
   double tempvx,tempvy,tempvz;

   double *vgmat      = par_rahman->vgmat;
   double *fgmat_p    = par_rahman->fgmat_p;
   double *fgmat_v    = par_rahman->fgmat_v;
   double *vtemps     = par_rahman->vtemps;
   double *vtempv     = par_rahman->vtempv;
   double *vexpdt     = par_rahman->vexpdt;
   double *veig       = par_rahman->veig;
   double *veigv      = par_rahman->veigv;
   double *roll_mtvv  = par_rahman->roll_mtvv;
   double *fv1        = par_rahman->fv1;
   double *fv2        = par_rahman->fv2;
   double c1_hm       = par_rahman->c1_hm;

   double *clatoms_vx = clatoms_pos->vx;
   double *clatoms_vy = clatoms_pos->vy;
   double *clatoms_vz = clatoms_pos->vz;

   int natm_tot       = clatoms_info->natm_tot;
   int myatm_start    = clatoms_info->myatm_start;
   int myatm_end      = clatoms_info->myatm_end;
 
/*--------------------------------------------------------------------------*/
/*  2.1) Evolve the particle velocities with NHC                            */

   for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
     arg = -therm_wdti4[iyosh]*therm_v_nhc[1][inhc]*wght;
     int_scr_sc[inhc] = exp(arg);
   }/*endfor*/

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     ktemp = therm_inhc_x[ipart];
     clatoms_vx[ipart]      *= int_scr_sc[ktemp];
     clatoms_roll_sc[ipart] *= int_scr_sc[ktemp];
   }/*endfor*/
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     clatoms_vy[ipart] *= int_scr_sc[therm_inhc_y[ipart]];
   }/*endfor*/
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     clatoms_vz[ipart] *= int_scr_sc[therm_inhc_z[ipart]];
   }/*endfor*/
 
/*--------------------------------------------------------------------------*/
/*  2.2) Evolve the particle velocities with vgmat                          */

   for(i=1;i<=9;i++){vtemps[i] = vgmat[i];}
   trace3     = c1_hm*(vgmat[1]+vgmat[5]+vgmat[9]);
   vtemps[1] += trace3; 
   vtemps[5] += trace3;
   vtemps[9] += trace3;

   diag33(vtemps,veig,veigv,fv1,fv2);
   for(i=1;i<=3;i++){vexpdt[i] = exp(-veig[i]*therm_wdti2[iyosh]);}
   for(i=1;i<=3;i++){
     joff = (i-1)*3 ;
     for(j=1;j<=3;j++){ 
        vtemps[j+joff] = veigv[j+joff]*vexpdt[i];
     }/*endfor*/
   }/*endfor*/
   n = 3;
   matmul_t2(veigv,vtemps,vtempv,n);
   matmul_t2(vtempv,roll_mtvv,vtemps,n);
   for(i=1;i<=9;i++){roll_mtvv[i]=vtemps[i];}
		   
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     tempvx =  clatoms_vx[ipart]*vtempv[1]
            +  clatoms_vy[ipart]*vtempv[2]
            +  clatoms_vz[ipart]*vtempv[3];
     tempvy =  clatoms_vx[ipart]*vtempv[4]
            +  clatoms_vy[ipart]*vtempv[5]
            +  clatoms_vz[ipart]*vtempv[6];
     tempvz =  clatoms_vx[ipart]*vtempv[7]
            +  clatoms_vy[ipart]*vtempv[8]
            +  clatoms_vz[ipart]*vtempv[9];
     clatoms_vx[ipart] = tempvx;
     clatoms_vy[ipart] = tempvy;
     clatoms_vz[ipart] = tempvz;
   }/*endfor*/

/*--------------------------------------------------------------------------*/
/*  2.3) Evolve the particle velocities with NHC                            */

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     ktemp = therm_inhc_x[ipart];
     clatoms_vx[ipart]      *= int_scr_sc[ktemp];
     clatoms_roll_sc[ipart] *= int_scr_sc[ktemp];
   }/*endfor*/
   for(ipart=1;ipart<=natm_tot;ipart++){
     clatoms_vy[ipart] *= int_scr_sc[therm_inhc_y[ipart]];
   }/*endfor*/
   for(ipart=1;ipart<=natm_tot;ipart++){
     clatoms_vz[ipart] *= int_scr_sc[therm_inhc_z[ipart]];
   }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/








/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_vel_vbox_upper(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                   CELL *cell,PAR_RAHMAN *par_rahman, 
                   double *baro_x_vol_nhc,double *baro_v_vol_nhc,
                   double *therm_wdti2,double *therm_wdti4,
                   double *therm_wdti8,
                   double **therm_x_nhc,double **therm_v_nhc,
                   int *therm_inhc_x,int *therm_inhc_y,
                   int *therm_inhc_z,int iyosh, 
                   int num_nhc, double *int_scr_sc, double *clatoms_roll_sc,
                   double wght, int len_nhc)

/*==========================================================================*/
/*Begin Routine*/{
/*==========================================================================*/
/*             Local variable declarations                                  */

  int i,inhc,ipart,ktemp,ichain,iii;
  double arg,aa,aa2,temp,dt,dtl;
  double ax,ay,az;
  double v_nhc_x,v_nhc_y,v_nhc_z;
  double dt2,dt24,nf;
  double gammax,gammay,gammaz;
  double sindt;
  double e2,e4,e6,e8;
  double vg12_vy0,vg13_vy0,vg13_vz0,vg23_vz0,tr_vgmat;
  double vx_temp,vy_temp,vz_temp;

  double c1_hm      = par_rahman->c1_hm;
  double *vgmat     = par_rahman->vgmat;
  double *fgmat_p   = par_rahman->fgmat_p;
  double *fgmat_v   = par_rahman->fgmat_v;
  double *vtemps    = par_rahman->vtemps;
  double *vtempv    = par_rahman->vtempv;
  double *vexpdt    = par_rahman->vexpdt;
  double *veig      = par_rahman->veig;
  double *veigv     = par_rahman->veigv;
  double *roll_mtvv = par_rahman->roll_mtvv;

  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  double *clatoms_mass         = clatoms_pos->mass;

  int natm_tot                 = clatoms_info->natm_tot;
  int nfree                    = clatoms_info->nfree;
  int myatm_start              = clatoms_info->myatm_start;
  int myatm_end                = clatoms_info->myatm_end;

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

  dt   = therm_wdti2[iyosh];
  dt2  = dt/2.0;      
  dt24 = dt2*dt2;

/*--------------------------------------------------------------------------*/
/* 2) Evolve the particle velocities with vgmat                             */


   tr_vgmat = (vgmat[1]+vgmat[5]+vgmat[9])*c1_hm;

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){

     v_nhc_x = therm_v_nhc[1][therm_inhc_x[ipart]];
     v_nhc_y = therm_v_nhc[1][therm_inhc_y[ipart]];
     v_nhc_z = therm_v_nhc[1][therm_inhc_z[ipart]];

     gammax = -(wght*v_nhc_x+tr_vgmat+vgmat[1]);
     gammay = -(wght*v_nhc_y+tr_vgmat+vgmat[5]);
     gammaz = -(wght*v_nhc_z+tr_vgmat+vgmat[9]);

     ax     =  exp(gammax*therm_wdti2[iyosh]);
     ay     =  exp(gammay*therm_wdti2[iyosh]);
     az     =  exp(gammaz*therm_wdti2[iyosh]);

     vg12_vy0 = vgmat[4]*clatoms_vy[ipart];
     vg13_vy0 = vgmat[7]*clatoms_vy[ipart];
     vg13_vz0 = vgmat[7]*clatoms_vz[ipart];
     vg23_vz0 = vgmat[8]*clatoms_vz[ipart];

     clatoms_vx[ipart] = clatoms_vx[ipart]*ax-vg13_vz0*dt2*(ax+az)
                       - vg12_vy0*dt2*(ax+ay) 
                       + vg23_vz0*vgmat[4]*dt24*(0.5*(ax+az) + ay);
     clatoms_vy[ipart] = clatoms_vy[ipart]*ay
                       - vg23_vz0*dt2*(az+ay);
     clatoms_vz[ipart] = clatoms_vz[ipart]*az;

   }/*endfor : ipart*/

/*---------------------------------------------------------------------------*/
  }/*end routine*/
/*===========================================================================*/







