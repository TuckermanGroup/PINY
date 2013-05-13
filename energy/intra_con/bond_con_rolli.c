/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: constraint_control.c                           */
/*                                                                          */
/* This routine controls the atom based constraint routines                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_bond_roll_i(BOND *bond,
                       CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                       CELL *cell,
                       INTRA_SCR *intra_scr,PTENS *ptens,BARO *baro,
                       double dt,int iter,
                       CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/* WARNING DANGER dt may not be timeinfo.dt so let it come in    */
/* as a parameter                                                */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

int ibond,nbond,i,iii;
int ibig,ibig1,ioff,iend,nnow;
int ktemp,iboff;
  int ntot_use;
  int nlen_use,nlen_now;
  int igo,ntot;
  int block_on            = bond->block_con_on;
  int iloop=0;

double dt2,dt22,fdot1,fdot2,dt22i;
double r122,r12,alam,vcons,r12i;
double fxo1,fyo1,fzo1,f_lnv_inc;
double fxo2,fyo2,fzo2;
double fxo1r,fyo1r,fzo1r;
double fxo2r,fyo2r,fzo2r;
double fxn1,fyn1,fzn1;
double fxn2,fyn2,fzn2;
double aa,aa2,arg2,poly,bb,dlen;
double e2,e4,e6,e8;
double baro_roll_scv;
double temp1,temp2; 
int ipart;
  /* Define local pointers                                                */
  double *intra_scr_x1    = intra_scr->x1;
  double *intra_scr_y1    = intra_scr->y1;
  double *intra_scr_z1    = intra_scr->z1;
  double *intra_scr_x2    = intra_scr->x2;
  double *intra_scr_y2    = intra_scr->y2;
  double *intra_scr_z2    = intra_scr->z2;
  double *intra_scr_x3    = intra_scr->x3;
  double *intra_scr_y3    = intra_scr->y3;
  double *intra_scr_z3    = intra_scr->z3;
  double *intra_scr_x4    = intra_scr->x4;
  double *intra_scr_y4    = intra_scr->y4;
  double *intra_scr_z4    = intra_scr->z4;
  double *intra_scr_m1    = intra_scr->m1;
  double *intra_scr_m2    = intra_scr->m2;
  double *intra_scr_eq    = intra_scr->eq;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p23   = intra_scr->p23;
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_fx2   = intra_scr->fx2;
  double *intra_scr_fy2   = intra_scr->fy2;
  double *intra_scr_fz2   = intra_scr->fz2;
  double *intra_scr_fx3   = intra_scr->fx3;
  double *intra_scr_fy3   = intra_scr->fy3;
  double *intra_scr_fz3   = intra_scr->fz3;
  double *intra_scr_fx4   = intra_scr->fx4;
  double *intra_scr_fy4   = intra_scr->fy4;
  double *intra_scr_fz4   = intra_scr->fz4;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *intra_scr_dx43  = intra_scr->dx43;
  double *intra_scr_dy43  = intra_scr->dy43;
  double *intra_scr_dz43  = intra_scr->dz43;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_mass    = clatoms_info->mass;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *clatoms_xold    = clatoms_info->xold;
  double *clatoms_yold    = clatoms_info->yold;
  double *clatoms_zold    = clatoms_info->zold;
  double *bond_eq_con     = bond->eq_con;
  double *bond_al_con     = bond->al_con;
  double *ptens_pvten_inc = ptens->pvten_inc;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *ptens_pvten_tmp2 = ptens->pvten_tmp_res;
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int *bond_jtyp_con      = bond->jtyp_con;
  int *iblock_size        = bond->iblock_con_size;
  int *iblock_conflict_1  = bond->iblock_con_conflict_1;
  int *iblock_conflict_2  = bond->iblock_con_conflict_2;
  int np_forc = class_comm_forc_pkg->num_proc;
  int myid_forc = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc = class_comm_forc_pkg->comm;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
/*==========================================================================*/
/* I) Flip the order of the constraints on even iteration numbers            */

  nbond = bond->ncon;
  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* II) Loop over bonds in strips of length nlen                             */

  dt2  = dt/2.0;
  dt22 = dt2*dt;
  dt22i = 1.0/dt22; 
  nlen_use = (intra_scr->nlen);
  ntot = (bond->ncon);
  ntot_use = bond->ncon_max;
  for(ibig=1;ibig <= ntot_use;ibig += nlen_use) {

/*----------------------------------------------------------------------*/
/*  A) Calculate Offsets                                                */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;
    iloop++;
/*----------------------------------------------------------------------*/
/*  B) Gather positions  masses, multipliers and req                    */
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      intra_scr_x1[iboff]   = clatoms_x[ktemp];
      intra_scr_y1[iboff]   = clatoms_y[ktemp];
      intra_scr_z1[iboff]   = clatoms_z[ktemp];
      intra_scr_x3[iboff]   = clatoms_xold[ktemp];
      intra_scr_y3[iboff]   = clatoms_yold[ktemp];
      intra_scr_z3[iboff]   = clatoms_zold[ktemp];
      intra_scr_m1[iboff]   = clatoms_mass[ktemp];
    }
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      intra_scr_x2[iboff]   = clatoms_x[ktemp];
      intra_scr_y2[iboff]   = clatoms_y[ktemp];
      intra_scr_z2[iboff]   = clatoms_z[ktemp];
      intra_scr_x4[iboff]   = clatoms_xold[ktemp];
      intra_scr_y4[iboff]   = clatoms_yold[ktemp];
      intra_scr_z4[iboff]   = clatoms_zold[ktemp];
      intra_scr_m2[iboff]   = clatoms_mass[ktemp];
    }/*endfor*/
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      intra_scr_eq[iboff] = bond_eq_con[bond_jtyp_con[ibond]];
    }/*endfor*/
    for(ibond=1;ibond<=nnow;ibond++){
      intra_scr_m1[ibond] = 1.0/intra_scr_m1[ibond];
      intra_scr_m2[ibond] = 1.0/intra_scr_m2[ibond];
    }/*endfor*/
/*----------------------------------------------------------------------*/
/*  C) Displacements                                                    */
    for(ibond=1;ibond <= nnow; ++ibond) {       
      intra_scr_dx12[ibond] = intra_scr_x1[ibond]-intra_scr_x2[ibond];
      intra_scr_dy12[ibond] = intra_scr_y1[ibond]-intra_scr_y2[ibond];
      intra_scr_dz12[ibond] = intra_scr_z1[ibond]-intra_scr_z2[ibond];
      intra_scr_dx43[ibond] = intra_scr_x3[ibond]-intra_scr_x4[ibond];
      intra_scr_dy43[ibond] = intra_scr_y3[ibond]-intra_scr_y4[ibond];
      intra_scr_dz43[ibond] = intra_scr_z3[ibond]-intra_scr_z4[ibond];
    }/*endfor*/
    if(cell->intra_perds == 1) {
     period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
     period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
    }
/*----------------------------------------------------------------------*/
/*  D) First time get force with old multiplier                         */
    if(iter==1){
      baro_roll_scv = baro->roll_scv;
      for(ibond=1;ibond <= nnow; ++ibond) {           
/* i) calculate the basis vectors (r12) of old positions and get force */
        r122 = intra_scr_dx43[ibond]*intra_scr_dx43[ibond]
              +intra_scr_dy43[ibond]*intra_scr_dy43[ibond]
              +intra_scr_dz43[ibond]*intra_scr_dz43[ibond];
        r12  = sqrt(r122);
        r12i = 1.0/r12;

        fxo1  = -intra_scr_dx43[ibond]*r12i;
        fyo1  = -intra_scr_dy43[ibond]*r12i;
        fzo1  = -intra_scr_dz43[ibond]*r12i;
        fxo2  = -fxo1;
        fyo2  = -fyo1;
        fzo2  = -fzo1;

        fxo1r =  fxo1*baro_roll_scv;
        fyo1r =  fyo1*baro_roll_scv;
        fzo1r =  fzo1*baro_roll_scv;
        fxo2r =  fxo2*baro_roll_scv;
        fyo2r =  fyo2*baro_roll_scv;
        fzo2r =  fzo2*baro_roll_scv;

/* ii) get the multiplier                                 */

        alam  = bond_al_con[(ibond+ibig1)];

/* iii) get the forces                                    */
        intra_scr_fx1[ibond] = alam*fxo1;
        intra_scr_fy1[ibond] = alam*fyo1;
        intra_scr_fz1[ibond] = alam*fzo1;
        intra_scr_fx2[ibond] = alam*fxo2;
        intra_scr_fy2[ibond] = alam*fyo2;
        intra_scr_fz2[ibond] = alam*fzo2;
        intra_scr_fx3[ibond] = alam*fxo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fy3[ibond] = alam*fyo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fz3[ibond] = alam*fzo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fx4[ibond] = alam*fxo2r*dt22*intra_scr_m2[ibond];
        intra_scr_fy4[ibond] = alam*fyo2r*dt22*intra_scr_m2[ibond];
        intra_scr_fz4[ibond] = alam*fzo2r*dt22*intra_scr_m2[ibond];
      /*endfor*/}
    /*endif*/}
/*----------------------------------------------------------------------*/
/*  E) Thereafter refine the multipliers  and get the force             */
    if(iter!=1){
      baro_roll_scv = baro->roll_scv;
      for(ibond=1;ibond <= nnow; ++ibond) {           
/* i) calculate the basis vectors (r12) of old positions and get force */
        r122  = intra_scr_dx43[ibond]*intra_scr_dx43[ibond]
               +intra_scr_dy43[ibond]*intra_scr_dy43[ibond]
               +intra_scr_dz43[ibond]*intra_scr_dz43[ibond];
        r12   =  sqrt(r122);
        r12i  = 1.0/r12;

        fxo1  = -intra_scr_dx43[ibond]*r12i;
        fyo1  = -intra_scr_dy43[ibond]*r12i;
        fzo1  = -intra_scr_dz43[ibond]*r12i;
        fxo2  = -fxo1;
        fyo2  = -fyo1;
        fzo2  = -fzo1;

        fxo1r =  fxo1*baro_roll_scv;
        fyo1r =  fyo1*baro_roll_scv;
        fzo1r =  fzo1*baro_roll_scv;
        fxo2r =  fxo2*baro_roll_scv;
        fyo2r =  fyo2*baro_roll_scv;
        fzo2r =  fzo2*baro_roll_scv;
/*ii) calculate the basis vectors (r12) of new positions and get force*/
        r122  = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
               +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
               +intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
        r12   =  sqrt(r122); 
        r12i  = 1.0/r12;

        fxn1  = -intra_scr_dx12[ibond]*r12i;
        fyn1  = -intra_scr_dy12[ibond]*r12i;
        fzn1  = -intra_scr_dz12[ibond]*r12i;
        fxn2  = -fxn1;
        fyn2  = -fyn1;
        fzn2  = -fzn1;
/*iii) get the present value of the constraint */

        vcons  = (r12-intra_scr_eq[ibond]);

/*iv) get the increment to the multiplier      */
        fdot1 = (fxn1*fxo1r+fyn1*fyo1r+fzn1*fzo1r)*intra_scr_m1[ibond];
        fdot2 = (fxn2*fxo2r+fyn2*fyo2r+fzn2*fzo2r)*intra_scr_m2[ibond];
        alam = dt22i*vcons/(fdot1+fdot2);
        bond_al_con[(ibond+ibig1)] += alam;
/*v) get the force on the atoms                */
        intra_scr_fx1[ibond] = alam*fxo1;
        intra_scr_fy1[ibond] = alam*fyo1;
        intra_scr_fz1[ibond] = alam*fzo1;
        intra_scr_fx2[ibond] = alam*fxo2;
        intra_scr_fy2[ibond] = alam*fyo2;
        intra_scr_fz2[ibond] = alam*fzo2;
        intra_scr_fx3[ibond] = alam*fxo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fy3[ibond] = alam*fyo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fz3[ibond] = alam*fzo1r*dt22*intra_scr_m1[ibond];
        intra_scr_fx4[ibond] = alam*fxo2r*dt22*intra_scr_m2[ibond];
        intra_scr_fy4[ibond] = alam*fyo2r*dt22*intra_scr_m2[ibond];
        intra_scr_fz4[ibond] = alam*fzo2r*dt22*intra_scr_m2[ibond];
      /*endfor*/}
    /*endif*/}    
/*----------------------------------------------------------------------*/
/*  F) Get Pressure tensor                                              */
    for(i=1;i<=9;i++){ptens_pvten_tmp[i]=0;}
    if((cell->iperd==2)||(cell->iperd==3)){
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy43[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz43[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx43[ibond]*intra_scr_fy1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[1] += intra_scr_p11[ibond];
        ptens_pvten_tmp[5] += intra_scr_p22[ibond];
        ptens_pvten_tmp[9] += intra_scr_p33[ibond];
        ptens_pvten_tmp[2] += intra_scr_p12[ibond];
      /*endfor*/}
      ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
      if(cell->iperd==3){
        for(ibond=1;ibond <= nnow; ++ibond) {           
          intra_scr_p13[ibond] = intra_scr_dx43[ibond]*intra_scr_fz1[ibond];
          intra_scr_p23[ibond] = intra_scr_dy43[ibond]*intra_scr_fz1[ibond];
        /*endfor*/}
        for(ibond=1;ibond <= nnow; ++ibond) {
          ptens_pvten_tmp[3] += intra_scr_p13[ibond];
          ptens_pvten_tmp[6] += intra_scr_p23[ibond];
        /*endfor*/}
        ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
        ptens_pvten_tmp[8] = ptens_pvten_tmp[6];
      /*endif*/}
     }else{
      for(i=1;i<=9;i++){ptens_pvten_tmp[i]=0;}
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[1] += intra_scr_p11[ibond];}
    /*endif*/}

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    ptens_pvten_tmp2[i] = ptens_pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(ptens_pvten_tmp2[1]), &(ptens_pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/


    for(i=1;i<=9;i++){ptens_pvten_inc[i]+=ptens_pvten_tmp[i];}    
/*----------------------------------------------------------------------*/
/*  G) Add the force to the velocities then the positions               */
    for(ibond=1;ibond <= nnow; ++ibond) {
      temp1 = (dt2*intra_scr_m1[ibond]);
      temp2 = (dt2*intra_scr_m2[ibond]);
      intra_scr_fx1[ibond] *= temp1;
      intra_scr_fy1[ibond] *= temp1;
      intra_scr_fz1[ibond] *= temp1;
      intra_scr_fx2[ibond] *= temp2;
      intra_scr_fy2[ibond] *= temp2;
      intra_scr_fz2[ibond] *= temp2;
    /*endfor*/}
    if(iter!=1){
      f_lnv_inc   = (ptens_pvten_tmp[1]+ptens_pvten_tmp[5]
                    +ptens_pvten_tmp[9]);
      baro->f_lnv_p += f_lnv_inc;
      baro->v_lnv   += f_lnv_inc*(baro->roll_scg)*dt2/(baro->mass_lnv);
    /*endif*/}     

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_2[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx3[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy3[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz3[iboff];
     }/*endfor*/
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx3[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy3[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz3[iboff];
     }/*endfor*/
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_2[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx4[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy4[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz4[iboff];
     }/*endfor*/
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx4[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy4[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz4[iboff];
     }/*endfor*/
    }/*endif*/

 }/* endfor ibig */
/*------------------------------------------------------------------------*/
/* H) Initial correction: must be done outside ibig loop                  */
/*                        Needed in Shake only if multipliers interpolated*/
  if(iter==1){
    f_lnv_inc   = (ptens_pvten_inc[1]+ptens_pvten_inc[5]+ptens_pvten_inc[9])
                - (ptens->pvten_inc_old[1]+ptens->pvten_inc_old[5]
                  +ptens->pvten_inc_old[9]);
    baro->f_lnv_p += f_lnv_inc;
    baro->v_lnv   += f_lnv_inc*(baro->roll_scg)*dt2/(baro->mass_lnv);
  }/*endif*/ 

/*==========================================================================*/
/* III) Flip back                                                           */

  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* IV) consistently re-calculate everything                                 */
  if(iter!=1){
     e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);e6=e4/(6.0*7.0);e8=e6/(8.0*9.0);
/*----------------------------------------------------------------------*/
/* A) new atom positions */
     aa = exp( dt2*(baro->v_lnv) );
     aa2 = aa*aa;
     arg2 = ((baro->v_lnv)*dt2)*((baro->v_lnv)*dt2);
     poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
     bb   = aa*poly;
     for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
       clatoms_x[ipart] = aa2*(clatoms_xold)[ipart]
                           +bb*(clatoms_vx)[ipart]*dt;
       clatoms_y[ipart] = aa2*(clatoms_yold)[ipart]
                           +bb*(clatoms_vy)[ipart]*dt;
       clatoms_z[ipart] = aa2*(clatoms_zold)[ipart]
                           +bb*(clatoms_vz)[ipart]*dt;
     /*endfor*/}
     baro->roll_scv = bb;
/*----------------------------------------------------------------------*/
/* B) new log(vol) */
     (baro->x_lnv) = (baro->x_lnv_o) + (baro->v_lnv)*dt;
     dlen = exp( (baro->x_lnv)-(baro->x_lnv_o) );
/*----------------------------------------------------------------------*/
/* C) get the new matrix of cell parameters and their inverse */
     for(i=1;i<=9;i++){cell->hmat[i] = baro->hmato[i]*dlen;}
     gethinv(cell->hmat,cell->hmati,&(cell->vol),cell->iperd);
     baro->vol = cell->vol;
  /*endif*/}

/*--------------------------------------------------------------------------*/
/*end routine */}
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_bond_roll_i(BOND *bond,
                        CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                        CELL *cell,
                        INTRA_SCR *intra_scr,PTENS *ptens,BARO *baro,
                        double dt,int iter,
                        CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/* WARNING DANGER dt may not be timeinfo.dt so let it come in    */
/* as a parameter                                                */

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

int i,ibond,nbond,iii;
int ibig,ibig1,ioff,iend,nnow;
int ktemp,iboff;
  int nlen_use,nlen_now;
  int igo,ntot;
  int block_on            = bond->block_con_on;
  int iloop=0;
  int ntot_use;

double dt2,fdot1,fdot2,dt2i;
double r122,r12,alam,vcons,r12i;
double fxo1,fyo1,fzo1;
double fxo2,fyo2,fzo2;
double fxo1r,fyo1r,fzo1r;
double fxo2r,fyo2r,fzo2r;
double vx1r,vy1r,vz1r;
double vx2r,vy2r,vz2r;
double vxg1r,vyg1r,vzg1r;
double vxg2r,vyg2r,vzg2r;
double f_lnv_inc;
double baro_v_lnv_g;
double temp1,temp2;
  /* Define local pointers                                                */
  double *intra_scr_x1    = intra_scr->x1;
  double *intra_scr_y1    = intra_scr->y1;
  double *intra_scr_z1    = intra_scr->z1;
  double *intra_scr_x2    = intra_scr->x2;
  double *intra_scr_y2    = intra_scr->y2;
  double *intra_scr_z2    = intra_scr->z2;
  double *intra_scr_vx1   = intra_scr->vx1;
  double *intra_scr_vy1   = intra_scr->vy1;
  double *intra_scr_vz1   = intra_scr->vz1;
  double *intra_scr_vx2   = intra_scr->vx2;
  double *intra_scr_vy2   = intra_scr->vy2;
  double *intra_scr_vz2   = intra_scr->vz2;
  double *intra_scr_sc_1   = intra_scr->sc_1;
  double *intra_scr_sc_2   = intra_scr->sc_2;
  double *intra_scr_m1    = intra_scr->m1;
  double *intra_scr_m2    = intra_scr->m2;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p23   = intra_scr->p23;
  double *intra_scr_fx1   = intra_scr->fx1;
  double *intra_scr_fy1   = intra_scr->fy1;
  double *intra_scr_fz1   = intra_scr->fz1;
  double *intra_scr_fx2   = intra_scr->fx2;
  double *intra_scr_fy2   = intra_scr->fy2;
  double *intra_scr_fz2   = intra_scr->fz2;
  double *intra_scr_dx12  = intra_scr->dx12;
  double *intra_scr_dy12  = intra_scr->dy12;
  double *intra_scr_dz12  = intra_scr->dz12;
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *clatoms_mass    = clatoms_info->mass;
  double *clatoms_roll_sc = clatoms_info->roll_sc;
  double *bond_al_con     = bond->al_con;
  double *ptens_pvten_inc = ptens->pvten_inc;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  double *ptens_pvten_tmp2 = ptens->pvten_tmp_res;
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int *iblock_size        = bond->iblock_con_size;
  int *iblock_conflict_1  = bond->iblock_con_conflict_1;
  int *iblock_conflict_2  = bond->iblock_con_conflict_2;
  int np_forc = class_comm_forc_pkg->num_proc;
  int myid_forc = class_comm_forc_pkg->myid;
  MPI_Comm comm_forc = class_comm_forc_pkg->comm;
/*==========================================================================*/
/* I) Flip the order of the constraints on even iteration numbers            */

  nbond = bond->ncon;
  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* II) Loop over bonds in strips of length nlen                             */

  dt2  = dt/2.0;
  dt2i = 1.0/dt2;
  nlen_use = (intra_scr->nlen);
  ntot = (bond->ncon);

  ntot_use = bond->ncon_max;
  for(ibig=1;ibig <= ntot_use;ibig += nlen_use) {
/*----------------------------------------------------------------------*/
/*  A) Calculate Offsets                                                */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;
    iloop++;
/*----------------------------------------------------------------------*/
/*  B) Gather positions  masses, multipliers and req                    */
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      intra_scr_x1[iboff]   = clatoms_x[ktemp];
      intra_scr_y1[iboff]   = clatoms_y[ktemp];
      intra_scr_z1[iboff]   = clatoms_z[ktemp];
      intra_scr_vx1[iboff]  = clatoms_vx[ktemp];
      intra_scr_vy1[iboff]  = clatoms_vy[ktemp];
      intra_scr_vz1[iboff]  = clatoms_vz[ktemp];
      intra_scr_m1[iboff]   = clatoms_mass[ktemp];
      intra_scr_sc_1[iboff] = clatoms_roll_sc[ktemp];
    }
    for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      intra_scr_x2[iboff]   = clatoms_x[ktemp];
      intra_scr_y2[iboff]   = clatoms_y[ktemp];
      intra_scr_z2[iboff]   = clatoms_z[ktemp];
      intra_scr_vx2[iboff]  = clatoms_vx[ktemp];
      intra_scr_vy2[iboff]  = clatoms_vy[ktemp];
      intra_scr_vz2[iboff]  = clatoms_vz[ktemp];
      intra_scr_m2[iboff]   = clatoms_mass[ktemp];
      intra_scr_sc_2[iboff] = clatoms_roll_sc[ktemp];
    }/*endfor*/
    for(ibond=1;ibond<=nnow;ibond++){
      intra_scr_m1[ibond] = 1.0/intra_scr_m1[ibond];
      intra_scr_m2[ibond] = 1.0/intra_scr_m2[ibond];
    }/*endfor*/
/*----------------------------------------------------------------------*/
/*  C) Displacements                                                    */
    for(ibond=1;ibond <= nnow; ++ibond) {       
      intra_scr_dx12[ibond] = intra_scr_x1[ibond]-intra_scr_x2[ibond];
      intra_scr_dy12[ibond] = intra_scr_y1[ibond]-intra_scr_y2[ibond];
      intra_scr_dz12[ibond] = intra_scr_z1[ibond]-intra_scr_z2[ibond];
    }/*endfor*/
    if(cell->intra_perds == 1) {
     period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
    }
/*----------------------------------------------------------------------*/
/*  D) First time get the force with the old multiplier                 */
    if(iter==1){
      for(ibond=1;ibond <= nnow; ++ibond) {           
/* i) calculate the basis vectors (r12) of old positions and get force */
         r122 = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
               +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
               +intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
         r12  = sqrt(r122);
         r12i = 1.0/r12;
 
         fxo1  = -intra_scr_dx12[ibond]*r12i;
         fyo1  = -intra_scr_dy12[ibond]*r12i;
         fzo1  = -intra_scr_dz12[ibond]*r12i;
         fxo2  = -fxo1;
         fyo2  = -fyo1;
         fzo2  = -fzo1;
/* ii) get the multiplier                                 */
         alam  = bond_al_con[(ibond+ibig1)];
/* iii) get the forces forces                             */
         intra_scr_fx1[ibond]  = alam*fxo1;
         intra_scr_fy1[ibond]  = alam*fyo1;
         intra_scr_fz1[ibond]  = alam*fzo1;
         intra_scr_fx2[ibond]  = alam*fxo2;
         intra_scr_fy2[ibond]  = alam*fyo2;
         intra_scr_fz2[ibond]  = alam*fzo2;
      /*endfor*/}
    /*endif*/}
/*----------------------------------------------------------------------*/
/*  E) Thereafter refine the multipliers                                */
    if(iter!=1){
      baro_v_lnv_g = baro->v_lnv_g;
      for(ibond=1;ibond <= nnow; ++ibond) {           
/* i) calculate the basis vectors (r12) of old positions and get force */
         r122  = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
                +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
                +intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
         r12   =  sqrt(r122);
         r12i = 1.0/r12;

         fxo1  = -intra_scr_dx12[ibond]*r12i;
         fyo1  = -intra_scr_dy12[ibond]*r12i;
         fzo1  = -intra_scr_dz12[ibond]*r12i;
         fxo2  = -fxo1;
         fyo2  = -fyo1;
         fzo2  = -fzo1;
/* ii) scaled forces */
         fxo1r =  fxo1*intra_scr_sc_1[ibond];
         fyo1r =  fyo1*intra_scr_sc_1[ibond];
         fzo1r =  fzo1*intra_scr_sc_1[ibond];
         fxo2r =  fxo2*intra_scr_sc_2[ibond];
         fyo2r =  fyo2*intra_scr_sc_2[ibond];
         fzo2r =  fzo2*intra_scr_sc_2[ibond];
/* iii) scaled velocities */
         vx1r  =   intra_scr_vx1[ibond]*intra_scr_sc_1[ibond];
         vy1r  =   intra_scr_vy1[ibond]*intra_scr_sc_1[ibond];
         vz1r  =   intra_scr_vz1[ibond]*intra_scr_sc_1[ibond];
         vx2r  =   intra_scr_vx2[ibond]*intra_scr_sc_2[ibond];
         vy2r  =   intra_scr_vy2[ibond]*intra_scr_sc_2[ibond];
         vz2r  =   intra_scr_vz2[ibond]*intra_scr_sc_2[ibond];
/* iv) box contribution to dr/dt */
         vxg1r  =  intra_scr_x1[ibond]*baro_v_lnv_g;
         vyg1r  =  intra_scr_y1[ibond]*baro_v_lnv_g;
         vzg1r  =  intra_scr_z1[ibond]*baro_v_lnv_g;
         vxg2r  =  intra_scr_x2[ibond]*baro_v_lnv_g;
         vyg2r  =  intra_scr_y2[ibond]*baro_v_lnv_g;
         vzg2r  =  intra_scr_z2[ibond]*baro_v_lnv_g;
/* v) combine */
         vx1r  = vx1r + vxg1r;
         vy1r  = vy1r + vyg1r;
         vz1r  = vz1r + vzg1r;
         vx2r  = vx2r + vxg2r;
         vy2r  = vy2r + vyg2r;
         vz2r  = vz2r + vzg2r;
/* vi) get the present value of the constraint */
         vcons = -(vx1r*fxo1+vy1r*fyo1+vz1r*fzo1)
                 -(vx2r*fxo2+vy2r*fyo2+vz2r*fzo2);
         fdot1 = (fxo1r*fxo1+fyo1r*fyo1+fzo1r*fzo1)*intra_scr_m1[ibond];
         fdot2 = (fxo2r*fxo2+fyo2r*fyo2+fzo2r*fzo2)*intra_scr_m2[ibond];
         alam  = dt2i*vcons/(fdot1+fdot2);
         bond_al_con[(ibond+ibig1)] += alam;
/*v) get the force on the atoms                */
         intra_scr_fx1[ibond]  = alam*fxo1;
         intra_scr_fy1[ibond]  = alam*fyo1;
         intra_scr_fz1[ibond]  = alam*fzo1;
         intra_scr_fx2[ibond]  = alam*fxo2;
         intra_scr_fy2[ibond]  = alam*fyo2;
         intra_scr_fz2[ibond]  = alam*fzo2;
      /*endfor*/}
    /*endif*/}    
/*----------------------------------------------------------------------*/
/*  F) Get Pressure tensor                                              */
    for(i=1;i<=9;i++){ptens_pvten_tmp[i]=0.0;}
    if((cell->iperd==2)||(cell->iperd==3)){
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx12[ibond]*intra_scr_fy1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[1] += intra_scr_p11[ibond];
        ptens_pvten_tmp[5] += intra_scr_p22[ibond];
        ptens_pvten_tmp[9] += intra_scr_p33[ibond];
        ptens_pvten_tmp[2] += intra_scr_p12[ibond];
      /*endfor*/}
      ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
      if(cell->iperd==3){
        for(ibond=1;ibond <= nnow; ++ibond) {           
          intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
          intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
        /*endfor*/}
        for(ibond=1;ibond <= nnow; ++ibond) {
          ptens_pvten_tmp[3] += intra_scr_p13[ibond];
          ptens_pvten_tmp[6] += intra_scr_p23[ibond];
        /*endfor*/}
        ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
        ptens_pvten_tmp[8] = ptens_pvten_tmp[6];
      /*endif*/}
    }else{
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_tmp[1] += intra_scr_p11[ibond];}
    }/*endif*/

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */


  if(np_forc > 1 ){
   for(i=1;i<=9;i++){
    ptens_pvten_tmp2[i] = ptens_pvten_tmp[i];
   }/*endfor*/
   Allreduce(&(ptens_pvten_tmp2[1]), &(ptens_pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
  }/*endif*/

    for(i=1;i<=9;i++){ptens_pvten_inc[i]+=ptens_pvten_tmp[i];}     
/*----------------------------------------------------------------------*/
/*  G) Add the force to the velocities                                  */
    for(ibond=1;ibond <= nnow; ++ibond) {
      temp1 = (dt2*intra_scr_m1[ibond]);
      temp2 = (dt2*intra_scr_m2[ibond]);
      intra_scr_fx1[ibond] *=temp1;
      intra_scr_fy1[ibond] *=temp1;
      intra_scr_fz1[ibond] *=temp1;
      intra_scr_fx2[ibond] *=temp2;
      intra_scr_fy2[ibond] *=temp2;
      intra_scr_fz2[ibond] *=temp2;
    /*endfor*/}

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_2[iloop]==0){igo=1;}}
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }/*endif*/

    if(iter!=1){
      f_lnv_inc   = (ptens_pvten_tmp[1]+ptens_pvten_tmp[5]
                    +ptens_pvten_tmp[9]);
      baro->f_lnv_p += f_lnv_inc;
      baro->v_lnv_g += f_lnv_inc*(baro->roll_scg)*dt2/(baro->mass_lnv);
    /*endif*/}     
  /* endfor: ibig */}
/*------------------------------------------------------------------------*/
/* H) Initial correction: must be done outside ibig loop                  */
/*                        Necessary in Rattle                             */
  if(iter==1){
    f_lnv_inc   = (ptens_pvten_inc[1]+ptens_pvten_inc[5]+ptens_pvten_inc[9])
                - (ptens->pvten_inc_old[1]+ptens->pvten_inc_old[5]
                  +ptens->pvten_inc_old[9]);
    baro->f_lnv_p += f_lnv_inc;
    baro->v_lnv_g += f_lnv_inc*(baro->roll_scg)*dt2/(baro->mass_lnv);
  }/*endif*/ 

/*==========================================================================*/
/* III) Flip back                                                            */

  if((iter % 2)==0){flip_bond_con(bond);}

/*--------------------------------------------------------------------------*/
/*end routine */}
/*==========================================================================*/
