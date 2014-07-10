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

void move_pos_box(CELL *cell,PAR_RAHMAN *par_rahman,
                  double *clatoms_x,double *clatoms_y,
                  double *clatoms_z,double *clatoms_xold,
                  double *clatoms_yold,double *clatoms_zold,
                  double *clatoms_vx,
                  double *clatoms_vy,double *clatoms_vz,
                  double dt,int natm_tot)

/*==========================================================================*/
    {/*Begin Routine*/
/*==========================================================================*/
/*             Local variable declarations                                  */

#include "../typ_defs/typ_mask.h"

  int i,j,ip,ipart,joff,n;

  double e2,e4,e6,e8;
  double tempx,tempy,tempz;
  double tempvx,tempvy,tempvz;
  double dt2,aa,arg2,poly;

  double *roll_mtx = par_rahman->roll_mtx;
  double *roll_mtv = par_rahman->roll_mtv;
  double *vtemps   = par_rahman->vtemps;
  double *vgmat    = par_rahman->vgmat;
  double *vtempx   = par_rahman->vtempx;
  double *vtempv   = par_rahman->vtempv;
  double *veig     = par_rahman->veig;
  double *veigv    = par_rahman->veigv;
  double *fv1      = par_rahman->fv1;
  double *fv2      = par_rahman->fv2;
  double *vexpdt   = par_rahman->vexpdt;
  double *vsindt   = par_rahman->vsindt;
  double *hmat_t   =  par_rahman->hmat_t;
  double *hmato    =  par_rahman->hmato;

  double *hmat     = cell->hmat;
  double *hmati    = cell->hmati;
  int iperd        = cell->iperd;

/*==========================================================================*/
/* 0) Useful constants                                                      */

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
  e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

  dt2  =  dt/2.0;

/*==========================================================================*/

  for(i=1;i<=9;i++){vtemps[i]=vgmat[i];}
  diag33(vtemps,veig,veigv,fv1,fv2);
  for(i=1;i<=3;i++){
    aa         = exp(dt2*(veig)[i]);
    vexpdt[i]  = aa*aa;
    arg2       = ((veig)[i]*dt2)*(veig[i]*dt2);
    poly       = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
    vsindt[i]  = aa*poly;
  }/*endfor*/

  for(i=1;i<=3;i++){
    joff = (i-1)*3 ;
    for(j=1;j<=3;j++){ 
      vtempx[j+joff] = veigv[j+joff]*vexpdt[i];
      vtempv[j+joff] = veigv[j+joff]*vsindt[i];
    }/*endfor*/
  }/*endfor*/

  n = 3;
  matmul_t2(veigv,vtempx,roll_mtx,n); 
  matmul_t2(veigv,vtempv,roll_mtv,n);
     	
  for(ipart=1;ipart<=(natm_tot); ++ipart) {
    tempx  =  clatoms_xold[ipart]*roll_mtx[1]
             +clatoms_yold[ipart]*roll_mtx[2]
             +clatoms_zold[ipart]*roll_mtx[3];
    tempy  =  clatoms_xold[ipart]*roll_mtx[4]
             +clatoms_yold[ipart]*roll_mtx[5]
             +clatoms_zold[ipart]*roll_mtx[6];
    tempz  =  clatoms_xold[ipart]*roll_mtx[7]
             +clatoms_yold[ipart]*roll_mtx[8]
             +clatoms_zold[ipart]*roll_mtx[9];
    tempvx = clatoms_vx[ipart]*roll_mtv[1]
            +clatoms_vy[ipart]*roll_mtv[2]
            +clatoms_vz[ipart]*roll_mtv[3];
    tempvy = clatoms_vx[ipart]*roll_mtv[4]
            +clatoms_vy[ipart]*roll_mtv[5]
            +clatoms_vz[ipart]*roll_mtv[6];
    tempvz = clatoms_vx[ipart]*roll_mtv[7]
            +clatoms_vy[ipart]*roll_mtv[8]
            +clatoms_vz[ipart]*roll_mtv[9];
    clatoms_x[ipart] = tempx+tempvx*dt;
    clatoms_y[ipart] = tempy+tempvy*dt;
    clatoms_z[ipart] = tempz+tempvz*dt;
  }/*endfor*/

/*==========================================================================*/
/* 5) Evolve the cell                                              */

  n = 3;
  matmul_t(hmato,veigv,hmat_t,n);
  for(i=1;i<=3;i++){
    joff = (i-1)*3 ;
    for(j=1;j<=3;j++){ 
      hmat_t[j+joff] *=  vexpdt[j];
    }/*endfor*/
  }/*endfor*/
  matmul_tt(hmat_t,veigv,hmat,n);
  gethinv(hmat,hmati,&(cell->vol),iperd);
  par_rahman->vol = cell->vol;

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void move_pos_box_upper(CELL *cell,PAR_RAHMAN *par_rahman,
                        double *clatoms_x,double *clatoms_y,
                        double *clatoms_z,double *clatoms_vx,
                        double *clatoms_vy,double *clatoms_vz,
                        double dt,int natm_tot,int myatm_start, int myatm_end)

/*==========================================================================*/
/*Begin Routine*/{
/*==========================================================================*/
/*             Local variable declarations                                  */

   int i,ip,iatm,iii,ires,nres;
   double dti24,vg48,dt_res,dt2_res,dt24_res;
   double vg15,vg19,vg59,dt2,dt24;
   double gamma11,gamma22,gamma33;
   double sinh11,sinh22,sinh33;
   double sindt,arg2;
   double e2,e4,e6,e8;
   double ax,ay,az,bx,by,bz;
   double vol;
   double hmat_old[10],hmat_temp[10];

   double *vgmat    = par_rahman->vgmat;
   double *hmat     = cell->hmat;
   double *hmati    = cell->hmati;
   int iperd        = cell->iperd;

/*==========================================================================*/
/* 0) Useful constants                                                      */

   e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);
   e6= e4/(6.0*7.0);e8=e6/(8.0*9.0);

   dt2     = dt/2.0;
   dt24    = dt2*dt2;
   vg48    = vgmat[4]*vgmat[8];

   nres     = 1;
   dt_res   = dt  /((double)nres);
   dt2_res  = dt2 /((double)nres);
   dt24_res = dt24/((double)(nres*nres));

   gamma11  = (vgmat[1]*dt_res)/2.0;
   arg2     = gamma11*gamma11;
   sinh11   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

   gamma22  = (vgmat[5]*dt_res)/2.0;
   arg2     = gamma22*gamma22;
   sinh22   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

   gamma33  = (vgmat[9]*dt_res)/2.0;
   arg2     = gamma33*gamma33;
   sinh33   = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

   ax       = dt_res*sinh11*exp(vgmat[1]*dt2_res);
   ay       = dt_res*sinh22*exp(vgmat[5]*dt2_res);
   az       = dt_res*sinh33*exp(vgmat[9]*dt2_res);

   bx       = exp(vgmat[1]*dt_res);
   by       = exp(vgmat[5]*dt_res);
   bz       = exp(vgmat[9]*dt_res);

/*==========================================================================*/
/* I) Evolve the particles from 0 to dt                                     */

   for(ires=1;ires<=nres;ires++){
     for(iatm=myatm_start;iatm<=myatm_end;iatm++){
  /*--------------------------------------------------------------------*/
  /* 1) Evolve the postitions from 0 to dt with exp(iL_corr*dt/2)      */
       clatoms_x[iatm] += (vgmat[4]*clatoms_y[iatm]*dt2_res
                           +0.5*vg48*dt24_res*clatoms_z[iatm] 
                           +vgmat[7]*clatoms_z[iatm]*dt2_res);
       clatoms_y[iatm] +=  (vgmat[8]*clatoms_z[iatm]*dt2_res);

  /*--------------------------------------------------------------------*/
  /* 2) Evolve the postitions from 0 to dt with exp(iL_reference*dt)    */

       clatoms_x[iatm] = ax*clatoms_vx[iatm] + bx*clatoms_x[iatm];
       clatoms_y[iatm] = ay*clatoms_vy[iatm] + by*clatoms_y[iatm];
       clatoms_z[iatm] = az*clatoms_vz[iatm] + bz*clatoms_z[iatm];

  /*--------------------------------------------------------------------*/
  /* 3) Evolve the postitions from dt2 to dt with exp(iL_correction*dt2)*/

       clatoms_x[iatm] += (vgmat[4]*clatoms_y[iatm]*dt2_res
                           +0.5*vg48*dt24_res*clatoms_z[iatm] 
                           +vgmat[7]*clatoms_z[iatm]*dt2_res);
       clatoms_y[iatm] += (vgmat[8]*clatoms_z[iatm]*dt2_res);
     }/*endfor : iatm*/
   }/*endfor : ires*/

/*==========================================================================*/
/* II) Evolve the cell matrix                                               */

   nres     = 5;
   dt_res   = dt/((double)nres);
   dt24_res = dt24/((double)(nres*nres));

   for(ires=1;ires<=nres;ires++){
  /*--------------------------------------------------------------------*/
  /* 1) Evolve the matrix from 0 to dt2 with exp(iL_correction*dt2)     */

     hmat[7] += (vgmat[4]*hmat[8]*dt_res);
  /*--------------------------------------------------------------------*/
  /* 2) Save the old h-matrix                                           */

     for(i=1;i<=9;i++){hmat_old[i] = hmat[i];}

  /*--------------------------------------------------------------------*/
  /* 3) Evolve the matrix from 0 to dt with exp(iL_reference*dt)        */

     hmat[1] = hmat_old[1]*exp(vgmat[1]*dt_res);
     hmat[5] = hmat_old[5]*exp(vgmat[5]*dt_res);
     hmat[9] = hmat_old[9]*exp(vgmat[9]*dt_res);

     vg15  = (vgmat[1]-vgmat[5])*dt_res/2.0;
     arg2  = vg15*vg15;
     sindt = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

     hmat[4] =  hmat_old[4]*exp(vgmat[1]*dt_res) 
              + hmat_old[5]*vgmat[4]*dt_res
                *exp((vgmat[1]+vgmat[5])*dt_res/2.0)*sindt;

     vg59  = (vgmat[5]-vgmat[9])*dt_res/2.0;
     arg2  = vg59*vg59;
     sindt = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

     hmat[8] =  hmat_old[8]*exp(vgmat[5]*dt_res) 
              + hmat_old[9]*vgmat[8]*dt_res
                *exp((vgmat[5]+vgmat[9])*dt_res/2.0)*sindt;

     vg19  = (vgmat[1]-vgmat[9])*dt_res;
     arg2  = vg19*vg19;
     sindt = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;

     hmat[7] =  hmat_old[7]*exp(vgmat[1]*dt_res) 
              + hmat_old[9]*vgmat[7]*dt_res
                *exp((vgmat[1]+vgmat[9])*dt_res)*sindt;

  /*--------------------------------------------------------------------*/
  /* 4) Evolve the matrix from dt2 to dt with exp(iL_correction*dt)     */

     hmat[7] += vgmat[4]*hmat[8]*dt_res;

   }/*endfor : ires*/
   
  /*--------------------------------------------------------------------*/
  /* 5) recalculate the inverse of h matrix and cell volume             */

   cell->vol = 0.0;
   for(i=1;i<=9;i++){hmati[i]=0.0;}
   if (iperd == 3) {
     vol = hmat[1] * hmat[5] * hmat[9] ;
     cell->vol = vol;
     hmati[1] = (hmat[5] * hmat[9]) / vol;
     hmati[5] = (hmat[1] * hmat[9]) / vol;
     hmati[9] = (hmat[1] * hmat[5]) / vol;
     hmati[4] = (- hmat[4] * hmat[9]) / vol;
     hmati[2] = 0.0;
     hmati[7] = (hmat[4] * hmat[8] - hmat[7] * hmat[5]) / vol;
     hmati[3] = 0.0;
     hmati[8] = (- hmat[8] * hmat[1]) / vol;
     hmati[6] = 0.0;
   }/*endif*/
   par_rahman->vol = cell->vol;

/*===============================================================*/
/* Errors */

  if((cell->vol)==0.0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");    
    printf("The present volume is zero.                 \n");
    printf("If this is not an error in your input data, \n");
    printf("contact technical support                  \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  } /*endif*/

/*--------------------------------------------------------------------------*/
     }/*end routine*/
/*==========================================================================*/






