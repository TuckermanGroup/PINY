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
                   CELL *cell,PAR_RAHMAN *par_rahman,double dt)

/*==========================================================================*/
    {/*Begin Routine*/
/*==========================================================================*/
/*             Local variable declarations                                  */

   int i,j,joff,ipart,ktemp,n;
   double arg,arg2,aa,aa2,poly,trace3;
   double tempvx,tempvy,tempvz;
   double tempfx,tempfy,tempfz;
   double dt2,dt4;
   double e2,e4,e6,e8;

   double *roll_mtf   = par_rahman->roll_mtf;
   double *roll_mtvv  = par_rahman->roll_mtvv;
   double *vgmat      = par_rahman->vgmat;
   double *vgmat_g    = par_rahman->vgmat_g;
   double *fgmat_p    = par_rahman->fgmat_p;
   double *fgmat_v    = par_rahman->fgmat_v;
   double *vtemps     = par_rahman->vtemps;
   double *vtempv     = par_rahman->vtempv;
   double *vtempf     = par_rahman->vtempf;
   double *vexpdt     = par_rahman->vexpdt;
   double *vsindt     = par_rahman->vsindt;
   double *veig       = par_rahman->veig;
   double *veigv      = par_rahman->veigv;
   double *fv1        = par_rahman->fv1;
   double *fv2        = par_rahman->fv2;
   double c1_hm       = par_rahman->c1_hm;

   double *clatoms_vx = clatoms_pos->vx;
   double *clatoms_vy = clatoms_pos->vy;
   double *clatoms_vz = clatoms_pos->vz;
   double *clatoms_fx = clatoms_pos->fx;
   double *clatoms_fy = clatoms_pos->fy;
   double *clatoms_fz = clatoms_pos->fz;
   double *clatoms_mass = clatoms_info->mass;

   int natm_tot       = clatoms_info->natm_tot;
   int myatm_start    = clatoms_info->myatm_start;
   int myatm_end      = clatoms_info->myatm_end;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

  dt2  =  dt/2.0;
  dt4  =  dt2/2.0;

/*--------------------------------------------------------------------------*/
/*  Evolve the particle velocities with vgmat                               */

   for(i=1;i<=9;i++){vtemps[i] = vgmat[i];}
   trace3     = c1_hm*(vgmat[1]+vgmat[5]+vgmat[9]);
   vtemps[1] += trace3;
   vtemps[5] += trace3;
   vtemps[9] += trace3;

   diag33(vtemps,veig,veigv,fv1,fv2);

   for(i=1;i<=3;i++){
     aa = exp(-dt4*(veig)[i]);
     vexpdt[i] = aa*aa;
     arg2 = ((veig)[i]*dt4)*((veig)[i]*dt4);
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     vsindt[i] = aa*poly;
   }/*endfor*/

   for(i=1;i<=3;i++){
     joff = (i-1)*3 ;
     for(j=1;j<=3;j++){ 
     	  vtempf[j+joff] = veigv[j+joff]*vsindt[i];     	  
     	  vtempv[j+joff] = veigv[j+joff]*vexpdt[i];
     }/*endfor*/
   }/*endfor*/
   n = 3;
   matmul_t2(veigv,vtempf,roll_mtf,n); 
   matmul_t2(veigv,vtempv,roll_mtvv,n);

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     tempvx =  clatoms_vx[ipart]*roll_mtvv[1]
            +  clatoms_vy[ipart]*roll_mtvv[2]
            +  clatoms_vz[ipart]*roll_mtvv[3];
     tempvy =  clatoms_vx[ipart]*roll_mtvv[4]
            +  clatoms_vy[ipart]*roll_mtvv[5]
            +  clatoms_vz[ipart]*roll_mtvv[6];
     tempvz =  clatoms_vx[ipart]*roll_mtvv[7]
            +  clatoms_vy[ipart]*roll_mtvv[8]
            +  clatoms_vz[ipart]*roll_mtvv[9];
     tempfx =  clatoms_fx[ipart]*roll_mtf[1]
            +  clatoms_fy[ipart]*roll_mtf[2]
            +  clatoms_fz[ipart]*roll_mtf[3];
     tempfy =  clatoms_fx[ipart]*roll_mtf[4]
            +  clatoms_fy[ipart]*roll_mtf[5]
            +  clatoms_fz[ipart]*roll_mtf[6];
     tempfz =  clatoms_fx[ipart]*roll_mtf[7]
            +  clatoms_fy[ipart]*roll_mtf[8]
            +  clatoms_fz[ipart]*roll_mtf[9];
     clatoms_vx[ipart] = tempvx+tempfx*dt2/clatoms_mass[ipart];
     clatoms_vy[ipart] = tempvy+tempfy*dt2/clatoms_mass[ipart];
     clatoms_vz[ipart] = tempvz+tempfz*dt2/clatoms_mass[ipart];
   }/*endfor*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_vel_vbox_upper(CLATOMS_POS *clatoms_pos, CLATOMS_INFO *clatoms_info,
                         CELL *cell,PAR_RAHMAN *par_rahman,double dt)

/*==========================================================================*/
/*Begin Routine*/{
/*==========================================================================*/
/*             Local variable declarations                                  */

  int i,ipart,ktemp,iii;
  double arg2,aa,poly,tempv,tempf;
  double dt2,dt4;
  double vg1,vg4,vg7,vg5,vg8,vg9;
  double vxo,vyo,vzo;
  double vsindt,vexpdt;
  double e2,e4,e6,e8;
  double tr_vgmat;

  double c1_hm      = par_rahman->c1_hm;
  double *vgmat     = par_rahman->vgmat;
  double *fgmat_p   = par_rahman->fgmat_p;
  double *fgmat_v   = par_rahman->fgmat_v;

  double *clatoms_vx           = clatoms_pos->vx;
  double *clatoms_vy           = clatoms_pos->vy;
  double *clatoms_vz           = clatoms_pos->vz;
  double *clatoms_fx           = clatoms_pos->fx;
  double *clatoms_fy           = clatoms_pos->fy;
  double *clatoms_fz           = clatoms_pos->fz;
  double *clatoms_mass         = clatoms_info->mass;

  int natm_tot                 = clatoms_info->natm_tot;
  int nfree                    = clatoms_info->nfree;
  int myatm_start              = clatoms_info->myatm_start;
  int myatm_end                = clatoms_info->myatm_end;

/*==========================================================================*/
/*==========================================================================*/
/* 0) Useful constants                                                      */

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

  dt2  =  dt/2.0;
  dt4  =  dt2/2.0;
  
/*--------------------------------------------------------------------------*/
/* 2) Evolve the particle velocities with vgmat                             */

   tr_vgmat = (vgmat[1]+vgmat[5]+vgmat[9])*c1_hm;
   vg1 = (tr_vgmat+vgmat[1]);
   vg5 = (tr_vgmat+vgmat[5]);
   vg9 = (tr_vgmat+vgmat[9]);
   vg4 = vgmat[4];
   vg7 = vgmat[7];
   vg8 = vgmat[8];

   for(ipart=myatm_start;ipart<=myatm_end;ipart++){

     vzo = clatoms_vz[ipart];
     vyo = clatoms_vy[ipart];
     vxo = clatoms_vx[ipart];

     aa = exp(-dt4*vg9);
     vexpdt = aa*aa;
     arg2 = (vg9*dt4)*(vg9*dt4);
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     vsindt = aa*poly;
     tempv = clatoms_vz[ipart]*vexpdt;
     tempf = clatoms_fz[ipart]*vsindt;
     clatoms_vz[ipart] = tempv + tempf*dt2/clatoms_mass[ipart]; 
     
     aa = exp(-dt4*vg5);
     vexpdt = aa*aa;
     arg2 = (vg5*dt4)*(vg5*dt4);
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     vsindt = aa*poly;
     tempv = clatoms_vy[ipart]*vexpdt;
     tempf = clatoms_fy[ipart]*vsindt;
     clatoms_vy[ipart] = tempv + tempf*dt2/clatoms_mass[ipart]
                        -vg8*vzo*vsindt*dt2;

     aa = exp(-dt4*vg1);
     vexpdt = aa*aa;
     arg2 = (vg1*dt4)*(vg1*dt4);
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     vsindt = aa*poly;
     tempv = clatoms_vx[ipart]*vexpdt;
     tempf = clatoms_fx[ipart]*vsindt;
     clatoms_vx[ipart] 
         = tempv + tempf*dt2/clatoms_mass[ipart]
          -(vg4*vyo + vg7*vzo)*vsindt*dt2;

   }/*endfor : ipart*/

/*---------------------------------------------------------------------------*/
  }/*end routine*/
/*===========================================================================*/







