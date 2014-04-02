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



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_bond(BOND *bond,
                CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                CELL *cell,INTRA_SCR *intra_scr,
                PTENS *ptens,double dt,int iter)

/* WARNING DANGER dt may not be timeinfo.dt so let it come in    */
/* as a parameter                                                */

{/*Begin Routine*/
  /*=======================================================================*/
  /*         Local Variable declarations                                   */
  
  int ibond,nbond;
  int ibig,ibig1,ioff,iend,nnow;
  int ktemp,iboff;
  double dt2,dt22,fdot1,fdot2,dt22i;
  double r122,r12,alam,vcons,r12i;
  double fxo1,fyo1,fzo1;
  double fxo2,fyo2,fzo2;
  double fxn1,fyn1,fzn1;
  double fxn2,fyn2,fzn2;
  double temp1,temp2;
  int nlen_use,nlen_now;
  int igo,ntot;
  int block_on            = bond->block_con_on;
  int iloop;

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
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p31   = intra_scr->p31;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p32   = intra_scr->p32;
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
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int *bond_jtyp_con      = bond->jtyp_con;
  int *iblock_size        = bond->iblock_con_size;
  int *iblock_conflict_1  = bond->iblock_con_conflict_1;
  int *iblock_conflict_2  = bond->iblock_con_conflict_2;

  /*=======================================================================*/
  /* I) Flip the order of the constraints on even iteration numbers        */

  nbond = bond->ncon;
  if((iter % 2)==0){flip_bond_con(bond);}

  /*=======================================================================*/
  /* II) Loop over bonds in strips of length nlen                          */

  dt2  = dt/2.0;
  dt22 = dt2*dt;
  dt22i = 1.0/dt22;

  nlen_use = (intra_scr->nlen);
  ntot = (bond->ncon);
  if(block_on==1){
   nlen_use  =  bond->nblock_size_con;
   ntot      =  (bond->nblock_size_con)*bond->nblock_con;
  }/*endif*/

  iloop = 0;
  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
    iloop++;

    /*----------------------------------------------------------------------*/
    /*  A) Calculate Offsets                                                */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    if(block_on==1){nlen_now = iblock_size[iloop];}
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;

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
      for(ibond=1;ibond <= nnow; ++ibond) {           

        /* i) calculate the basis vectors of old positions and get force */
        r122 = (intra_scr_dx43[ibond]*intra_scr_dx43[ibond]
               +intra_scr_dy43[ibond]*intra_scr_dy43[ibond]
               +intra_scr_dz43[ibond]*intra_scr_dz43[ibond]);
        r12  = sqrt(r122);
        r12i = 1.0/r12;
        
        fxo1  = -intra_scr_dx43[ibond]*r12i;
        fyo1  = -intra_scr_dy43[ibond]*r12i;
        fzo1  = -intra_scr_dz43[ibond]*r12i;
        fxo2  = -fxo1;
        fyo2  = -fyo1;
        fzo2  = -fzo1;

        /* ii) get the multiplier                                 */

        alam  = bond_al_con[(ibond+ibig1)];

        /* iii) get the forces                             */
        intra_scr_fx1[ibond] = alam*fxo1;
        intra_scr_fy1[ibond] = alam*fyo1;
        intra_scr_fz1[ibond] = alam*fzo1;
        intra_scr_fx2[ibond] = alam*fxo2;
        intra_scr_fy2[ibond] = alam*fyo2;
        intra_scr_fz2[ibond] = alam*fzo2;
      }/*endfor*/
    }/*endif*/
    /*----------------------------------------------------------------------*/
    /*  E) Thereafter refine the multipliers  and get the force             */
    if(iter!=1){
      for(ibond=1;ibond <= nnow; ++ibond) {           
        /* i) calculate the basis vectors of old positions and get force */
        r122  = ( intra_scr_dx43[ibond]*intra_scr_dx43[ibond]
                 +intra_scr_dy43[ibond]*intra_scr_dy43[ibond]
                 +intra_scr_dz43[ibond]*intra_scr_dz43[ibond]);
        r12   =  sqrt(r122);
        r12i  = 1.0/r12;

        fxo1  = -intra_scr_dx43[ibond]*r12i;
        fyo1  = -intra_scr_dy43[ibond]*r12i;
        fzo1  = -intra_scr_dz43[ibond]*r12i;
        fxo2  = -fxo1;
        fyo2  = -fyo1;
        fzo2  = -fzo1;
        /*ii) calculate the basis vectors of new positions and get force*/
        r122  = (intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
                 +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
                 +intra_scr_dz12[ibond]*intra_scr_dz12[ibond]);
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
        fdot1 = (fxn1*fxo1+fyn1*fyo1+fzn1*fzo1)*intra_scr_m1[ibond];
        fdot2 = (fxn2*fxo2+fyn2*fyo2+fzn2*fzo2)*intra_scr_m2[ibond];
        alam = dt22i*vcons/(fdot1+fdot2);
        bond_al_con[(ibond+ibig1)] += alam;
        /*v) get the force on the atoms                */
        intra_scr_fx1[ibond] = alam*fxo1;
        intra_scr_fy1[ibond] = alam*fyo1;
        intra_scr_fz1[ibond] = alam*fzo1;
        intra_scr_fx2[ibond] = alam*fxo2;
        intra_scr_fy2[ibond] = alam*fyo2;
        intra_scr_fz2[ibond] = alam*fzo2;
      }/*endfor*/
    }/*endif*/
    /*----------------------------------------------------------------------*/
    /*  F) Get Pressure tensor                                              */
    if((cell->iperd==2)||(cell->iperd==3)){

      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy43[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz43[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx43[ibond]*intra_scr_fy1[ibond];
        intra_scr_p21[ibond] = intra_scr_dy43[ibond]*intra_scr_fx1[ibond];
      }/*endfor*/

      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_inc[1] += intra_scr_p11[ibond];
        ptens_pvten_inc[5] += intra_scr_p22[ibond];
        ptens_pvten_inc[9] += intra_scr_p33[ibond];
        ptens_pvten_inc[2] += intra_scr_p12[ibond];
        ptens_pvten_inc[4] += intra_scr_p21[ibond];
      }/*endfor*/

      if(cell->iperd==3){
        for(ibond=1;ibond <= nnow; ++ibond) {           
          intra_scr_p13[ibond] = intra_scr_dx43[ibond]*intra_scr_fz1[ibond];
          intra_scr_p31[ibond] = intra_scr_dz43[ibond]*intra_scr_fx1[ibond];
          intra_scr_p23[ibond] = intra_scr_dy43[ibond]*intra_scr_fz1[ibond];
          intra_scr_p32[ibond] = intra_scr_dz43[ibond]*intra_scr_fy1[ibond];
        }/*endfor*/

        for(ibond=1;ibond <= nnow; ++ibond) {
          ptens_pvten_inc[3] += intra_scr_p13[ibond];
          ptens_pvten_inc[7] += intra_scr_p31[ibond];
          ptens_pvten_inc[6] += intra_scr_p23[ibond];
          ptens_pvten_inc[8] += intra_scr_p32[ibond];
        }/*endfor*/
      }/*endif*/
     }else{
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_inc[1] += intra_scr_p11[ibond];}
     }/*endif*/

    /*--------------------------------------------------------------------*/
    /*  G) Add the force to the velocities then the positions             */

    for(ibond=1;ibond <= nnow; ++ibond) {
      temp1 = (dt2*intra_scr_m1[ibond]);
      temp2 = (dt2*intra_scr_m2[ibond]);
      intra_scr_fx1[ibond] *= temp1;
      intra_scr_fy1[ibond] *= temp1;
      intra_scr_fz1[ibond] *= temp1;
      intra_scr_fx2[ibond] *= temp2;
      intra_scr_fy2[ibond] *= temp2;
      intra_scr_fz2[ibond] *= temp2;
    }/*endfor*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }/*endfor*/
    }else{
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }/*endfor*/
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_2[iloop]==0){igo=1;}}
    if(igo==0){
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }else{
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }/*endif*/

    for(ibond=1;ibond <= nnow; ++ibond) {
      intra_scr_fx1[ibond] *= dt;
      intra_scr_fy1[ibond] *= dt;
      intra_scr_fz1[ibond] *= dt;
      intra_scr_fx2[ibond] *= dt;
      intra_scr_fy2[ibond] *= dt;
      intra_scr_fz2[ibond] *= dt;
    }/*endfor*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz1[iboff];
     }/*endfor*/
    }else{
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz1[iboff];
     }/*endfor*/
    }/*endif*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_2[iloop]==0){igo=1;}}
    if(igo==0){
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }else{
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }/*endif*/

  }/* endfor ibig */
  /*========================================================================*/
  /* III) Flip back                                                         */

  if((iter % 2)==0){flip_bond_con(bond);}
  
  /*------------------------------------------------------------------------*/
}/*end routine */
  /*=======================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_bond(BOND *bond,
                 CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                 CELL *cell,INTRA_SCR *intra_scr,
                 PTENS *ptens,double dt,int iter)

/* WARNING DANGER dt may not be timeinfo.dt so let it come in    */
/* as a parameter                                                */

{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

  int ibond,nbond;
  int ibig,ibig1,ioff,iend,nnow;
  int ktemp,iboff;
  int nlen_use,nlen_now;
  int igo,ntot;
  int block_on            = bond->block_con_on;
  int iloop;
  double dt2,fdot1,fdot2,dt2i;
  double r122,r12,alam,vcons,r12i;
  double fxo1,fyo1,fzo1;
  double fxo2,fyo2,fzo2;
  double temp1,temp2;

  /* Define local pointers                                                */
  double *intra_scr_x1    = intra_scr->x1;
  double *intra_scr_y1    = intra_scr->y1;
  double *intra_scr_z1    = intra_scr->z1;
  double *intra_scr_x2    = intra_scr->x2;
  double *intra_scr_y2    = intra_scr->y2;
  double *intra_scr_z2    = intra_scr->z2;
  double *intra_scr_vx1    = intra_scr->vx1;
  double *intra_scr_vy1    = intra_scr->vy1;
  double *intra_scr_vz1    = intra_scr->vz1;
  double *intra_scr_vx2    = intra_scr->vx2;
  double *intra_scr_vy2    = intra_scr->vy2;
  double *intra_scr_vz2    = intra_scr->vz2;
  double *intra_scr_m1    = intra_scr->m1;
  double *intra_scr_m2    = intra_scr->m2;
  double *intra_scr_p11   = intra_scr->p11;
  double *intra_scr_p22   = intra_scr->p22;
  double *intra_scr_p33   = intra_scr->p33;
  double *intra_scr_p21   = intra_scr->p21;
  double *intra_scr_p12   = intra_scr->p12;
  double *intra_scr_p31   = intra_scr->p31;
  double *intra_scr_p13   = intra_scr->p13;
  double *intra_scr_p32   = intra_scr->p32;
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
  double *clatoms_mass    = clatoms_info->mass;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *bond_al_con     = bond->al_con;
  double *ptens_pvten_inc = ptens->pvten_inc;
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int *iblock_size        = bond->iblock_con_size;
  int *iblock_conflict_1  = bond->iblock_con_conflict_1;
  int *iblock_conflict_2  = bond->iblock_con_conflict_2;

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
  if(block_on==1){
   nlen_use  =  bond->nblock_size_con;
   ntot      =  (bond->nblock_size_con)*bond->nblock_con;
  }/*endif*/

  iloop = 0;
  for(ibig=1;ibig <= ntot;ibig += nlen_use) {
   iloop++;

/*----------------------------------------------------------------------*/
/*  A) Calculate Offsets                                                */
    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
    if(block_on==1){nlen_now = iblock_size[iloop];}
    iend = MIN(ntot,ibig1+nlen_now);
    nnow = iend-ibig1;
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
      }/*endfor*/
    }/*endif*/
/*----------------------------------------------------------------------*/
/*  E) Thereafter refine the multipliers                                */
    if(iter!=1){
      for(ibond=1;ibond <= nnow; ++ibond) {           
/* i) calculate the basis vectors (r12) of old positions and get force */
         r122  = intra_scr_dx12[ibond]*intra_scr_dx12[ibond]
                +intra_scr_dy12[ibond]*intra_scr_dy12[ibond]
                +intra_scr_dz12[ibond]*intra_scr_dz12[ibond];
         r12   =  sqrt(r122);
         r12i  =  1.0/r12;

         fxo1  = -intra_scr_dx12[ibond]*r12i;
         fyo1  = -intra_scr_dy12[ibond]*r12i;
         fzo1  = -intra_scr_dz12[ibond]*r12i;
         fxo2  = -fxo1;
         fyo2  = -fyo1;
         fzo2  = -fzo1;
         /*iii) get the present value of the constraint */

         vcons = -( intra_scr_vx1[ibond]*fxo1 +intra_scr_vy1[ibond]*fyo1
                                              +intra_scr_vz1[ibond]*fzo1 )
                 -( intra_scr_vx2[ibond]*fxo2 +intra_scr_vy2[ibond]*fyo2
                                              +intra_scr_vz2[ibond]*fzo2 );

         fdot1 = (fxo1*fxo1+fyo1*fyo1+fzo1*fzo1)*intra_scr_m1[ibond];
         fdot2 = (fxo2*fxo2+fyo2*fyo2+fzo2*fzo2)*intra_scr_m2[ibond];

         alam = dt2i*vcons/(fdot1+fdot2);

         bond_al_con[(ibond+ibig1)] += alam;
/*v) get the force on the atoms                */
         intra_scr_fx1[ibond]  = alam*fxo1;
         intra_scr_fy1[ibond]  = alam*fyo1;
         intra_scr_fz1[ibond]  = alam*fzo1;
         intra_scr_fx2[ibond]  = alam*fxo2;
         intra_scr_fy2[ibond]  = alam*fyo2;
         intra_scr_fz2[ibond]  = alam*fzo2;
      }/*endfor*/
    }/*endif*/
/*----------------------------------------------------------------------*/
/*  F) Get Pressure tensor                                              */
    if((cell->iperd==2)||(cell->iperd==3)){
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p21[ibond] = intra_scr_dy12[ibond]*intra_scr_fx1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_inc[1] += intra_scr_p11[ibond];
        ptens_pvten_inc[5] += intra_scr_p22[ibond];
        ptens_pvten_inc[9] += intra_scr_p33[ibond];
        ptens_pvten_inc[2] += intra_scr_p12[ibond];
        ptens_pvten_inc[4] += intra_scr_p21[ibond];
      }/*endfor*/
      if(cell->iperd==3){
       for(ibond=1;ibond <= nnow; ++ibond) {           
         intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
         intra_scr_p31[ibond] = intra_scr_dz12[ibond]*intra_scr_fx1[ibond];
         intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
         intra_scr_p32[ibond] = intra_scr_dz12[ibond]*intra_scr_fy1[ibond];
       }/*endfor*/
       for(ibond=1;ibond <= nnow; ++ibond) {
         ptens_pvten_inc[3] += intra_scr_p13[ibond];
         ptens_pvten_inc[7] += intra_scr_p31[ibond];
         ptens_pvten_inc[6] += intra_scr_p23[ibond];
         ptens_pvten_inc[8] += intra_scr_p32[ibond];
       }/*endfor*/
     }/*endif*/
   }else{
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
      /*endfor*/}
      for(ibond=1;ibond <= nnow; ++ibond) {
        ptens_pvten_inc[1] += intra_scr_p11[ibond];}
    }/*endif*/

    /*----------------------------------------------------------------------*/
    /*  G) Add the force to the velocities                                  */

    for(ibond=1;ibond <= nnow; ++ibond) {
      temp1 = (dt2*intra_scr_m1[ibond]);
      temp2 = (dt2*intra_scr_m2[ibond]);
      intra_scr_fx1[ibond] *= temp1;
      intra_scr_fy1[ibond] *= temp1;
      intra_scr_fz1[ibond] *= temp1;
      intra_scr_fx2[ibond] *= temp2;
      intra_scr_fy2[ibond] *= temp2;
      intra_scr_fz2[ibond] *= temp2;
    }/*endfor*/

    igo = 0;
    if(block_on==1){if(iblock_conflict_1[iloop]==0){igo=1;}}
    if(igo==0){
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx1[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy1[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz1[iboff];
     }
    }else{
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
#ifndef NO_PRAGMA
#pragma NOVECTOR
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }else{
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j2_con[ibond];
      clatoms_vx[ktemp] +=  intra_scr_fx2[iboff];
      clatoms_vy[ktemp] +=  intra_scr_fy2[iboff];
      clatoms_vz[ktemp] +=  intra_scr_fz2[iboff];
     }/*endfor*/
    }/*endif*/

  }/* endfor ibig */
/*==========================================================================*/
/* III) Flip back                                                            */

  if((iter % 2)==0){flip_bond_con(bond);}

/*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void flip_bond_con(BOND *bond) 

{/*Begin Routine*/
  /*========================================================================*/
  /*          Local Variable declarations                                   */

  int ibond,nbond2,nbond,kbond;
  int itemp;
  double temper;

  /* Define local pointers                                                 */
  int *bond_j1_con    = bond->j1_con;
  int *bond_j2_con    = bond->j2_con;
  int *bond_jtyp_con  = bond->jtyp_con;
  double *bond_al_con = bond->al_con;

  /*========================================================================*/

  nbond  = bond->ncon;
  nbond2 = (bond->ncon)/2;
  for(ibond=1;ibond<=nbond2;ibond++){
    kbond                 = nbond-ibond+1;
    itemp                 = bond_j1_con[ibond];
    bond_j1_con[ibond]    = bond_j1_con[kbond];
    bond_j1_con[kbond]    = itemp;
    itemp                 = bond_j2_con[ibond];
    bond_j2_con[ibond]    = bond_j2_con[kbond];
    bond_j2_con[kbond]    = itemp;
    itemp                 = bond_jtyp_con[ibond];
    bond_jtyp_con[ibond]  = bond_jtyp_con[kbond];
    bond_jtyp_con[kbond]  = itemp;
    temper                = bond_al_con[ibond];
    bond_al_con[ibond]    = bond_al_con[kbond];
    bond_al_con[kbond]    = temper;
  }/*endfor*/

}/*end routine */
/*==========================================================================*/














