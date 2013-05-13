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

void shake_bond_roll_f(BOND *bond,
                       CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                       CELL *cell,
                       INTRA_SCR *intra_scr,PTENS *ptens,
                       PAR_RAHMAN *par_rahman,double dt,int iter,
                       CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*   WARNING DANGER dt may not be timeinfo.dt so let it come in             */
/*   as a parameter                                                         */
/*==========================================================================*/
    {/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

  int ibond,i,iii;
  int ktemp,iboff;
  int ibig,ibig1,ioff,iend,nnow;
  int nlen_now;
  int igo;
  int ipart;
  int j,joff,n;

  double dt2,dt22,fdot1,fdot2,aa,dt22i;
  double r122,r12,alam,vcons,r12i;
  double fxo1,fyo1,fzo1;
  double fxo2,fyo2,fzo2;
  double fxo1r,fyo1r,fzo1r;
  double fxo2r,fyo2r,fzo2r;
  double fxn1,fyn1,fzn1;
  double fxn2,fyn2,fzn2;
  double arg2,poly;
  double e2,e4,e6,e8;
  double tempx,tempy,tempz;
  double tempvx,tempvy,tempvz;
  double temp1,temp2;

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
  int nlen_use            = intra_scr->nlen;

  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *clatoms_mass    = clatoms_info->mass;
  double *clatoms_xold    = clatoms_info->xold;
  double *clatoms_yold    = clatoms_info->yold;
  double *clatoms_zold    = clatoms_info->zold;
  int myatm_start         = clatoms_info->myatm_start;
  int myatm_end           = clatoms_info->myatm_end;

  double *bond_eq_con     = bond->eq_con;
  double *bond_al_con     = bond->al_con;
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int *bond_jtyp_con      = bond->jtyp_con;
  int ntot                = bond->ncon;
  int ntot_use            = bond->ncon_max;

  double *pvten_inc       = ptens->pvten_inc;
  double *pvten_tmp       = ptens->pvten_tmp;
  double *pvten_tmp2      = ptens->pvten_tmp_res;
  double *pvten_inc_old   = ptens->pvten_inc_old;

  int intra_perds         = cell->intra_perds;
  int iperd               = cell->iperd;
  int hmat_cons_typ       = cell->hmat_cons_typ;
  int hmat_int_typ        = cell->hmat_int_typ;
  double *hmat            = cell->hmat;
  double *hmati           = cell->hmati;

  int np_forc             = class_comm_forc_pkg->num_proc;
  MPI_Comm comm_forc      = class_comm_forc_pkg->comm;

  double *roll_mtv        = par_rahman->roll_mtv;
  double *roll_mtx        = par_rahman->roll_mtx;
  double *fgmat_p         = par_rahman->fgmat_p;
  double *vgmat           = par_rahman->vgmat;
  double roll_scg         = par_rahman->roll_scg;
  double mass_hm          = par_rahman->mass_hm;
  double *vtemps          = par_rahman->vtemps;
  double *vtempv          = par_rahman->vtempv;
  double *vtempx          = par_rahman->vtempx;
  double *veig            = par_rahman->veig;
  double *veigv           = par_rahman->veigv;
  double *fv1             = par_rahman->fv1;
  double *fv2             = par_rahman->fv2;
  double *vsindt          = par_rahman->vsindt;
  double *vexpdt          = par_rahman->vexpdt;
  double *hmato           = par_rahman->hmato;
  double *hmat_t          = par_rahman->hmat_t;

  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);e6=e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* I) Flip the order of the constraints on even iteration numbers            */

  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* II) Loop over bonds in strips of length nlen                             */

  dt2  = dt/2.0;
  dt22 = dt2*dt;
  dt22i = 1.0/dt22;

  for(ibig=1;ibig <= ntot_use;ibig += nlen_use) {

/*----------------------------------------------------------------------*/
/*  A) Calculate Offsets                                                */

    ibig1    = ibig-1;
    ioff     = -ibig1;
    nlen_now = nlen_use;
    iend     = MIN(ntot,ibig1+nlen_now);
    nnow     = iend-ibig1;

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
    if(intra_perds == 1) {
      period(nnow,intra_scr_dx12,intra_scr_dy12,intra_scr_dz12,cell);
      period(nnow,intra_scr_dx43,intra_scr_dy43,intra_scr_dz43,cell);
    }
/*----------------------------------------------------------------------*/
/*  D) First time get force with old multiplier                         */
    if(iter==1){
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
/* ii) rolled forces*/
         fxo1r =  fxo1*roll_mtv[1]+fyo1*roll_mtv[2]
               +  fzo1*roll_mtv[3];
         fyo1r =  fxo1*roll_mtv[4]+fyo1*roll_mtv[5]
               +  fzo1*roll_mtv[6];
         fzo1r =  fxo1*roll_mtv[7]+fyo1*roll_mtv[8]
               +  fzo1*roll_mtv[9];
         fxo2r =  fxo2*roll_mtv[1]+fyo2*roll_mtv[2]
               +  fzo2*roll_mtv[3];
         fyo2r =  fxo2*roll_mtv[4]+fyo2*roll_mtv[5]
               +  fzo2*roll_mtv[6];
         fzo2r =  fxo2*roll_mtv[7]+fyo2*roll_mtv[8]
               +  fzo2*roll_mtv[9];
/* ii) get the multiplier                                 */
         alam  = bond_al_con[(ibond+ibig1)];
/* iii) get the forces forces                             */
         intra_scr_fx1[ibond]  = alam*fxo1;
         intra_scr_fy1[ibond]  = alam*fyo1;
         intra_scr_fz1[ibond]  = alam*fzo1;
         intra_scr_fx2[ibond]  = alam*fxo2;
         intra_scr_fy2[ibond]  = alam*fyo2;
         intra_scr_fz2[ibond]  = alam*fzo2;
         intra_scr_fx3[ibond]  = alam*fxo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fy3[ibond]  = alam*fyo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fz3[ibond]  = alam*fzo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fx4[ibond]  = alam*fxo2r*dt22*intra_scr_m2[ibond];
         intra_scr_fy4[ibond]  = alam*fyo2r*dt22*intra_scr_m2[ibond];
         intra_scr_fz4[ibond]  = alam*fzo2r*dt22*intra_scr_m2[ibond];
      }/*endfor*/
    }/*endif*/
/*----------------------------------------------------------------------*/
/*  E) Thereafter refine the multipliers  and get the force             */
    if(iter!=1){
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
/* ii) rolled forces*/
         fxo1r =  fxo1*roll_mtv[1]+fyo1*roll_mtv[2]
               +  fzo1*roll_mtv[3];
         fyo1r =  fxo1*roll_mtv[4]+fyo1*roll_mtv[5]
               +  fzo1*roll_mtv[6];
         fzo1r =  fxo1*roll_mtv[7]+fyo1*roll_mtv[8]
               +  fzo1*roll_mtv[9];
         fxo2r =  fxo2*roll_mtv[1]+fyo2*roll_mtv[2]
               +  fzo2*roll_mtv[3];
         fyo2r =  fxo2*roll_mtv[4]+fyo2*roll_mtv[5]
               +  fzo2*roll_mtv[6];
         fzo2r =  fxo2*roll_mtv[7]+fyo2*roll_mtv[8]
               +  fzo2*roll_mtv[9];
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
         intra_scr_fx1[ibond]  = alam*fxo1;
         intra_scr_fy1[ibond]  = alam*fyo1;
         intra_scr_fz1[ibond]  = alam*fzo1;
         intra_scr_fx2[ibond]  = alam*fxo2;
         intra_scr_fy2[ibond]  = alam*fyo2;
         intra_scr_fz2[ibond]  = alam*fzo2;
         intra_scr_fx3[ibond]  = alam*fxo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fy3[ibond]  = alam*fyo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fz3[ibond]  = alam*fzo1r*dt22*intra_scr_m1[ibond];
         intra_scr_fx4[ibond]  = alam*fxo2r*dt22*intra_scr_m2[ibond];
         intra_scr_fy4[ibond]  = alam*fyo2r*dt22*intra_scr_m2[ibond];
         intra_scr_fz4[ibond]  = alam*fzo2r*dt22*intra_scr_m2[ibond];
      }/*endfor*/
    }/*endif*/    
/*----------------------------------------------------------------------*/
/*  F) Get Pressure tensor                                              */
    for(i=1;i<=9;i++){pvten_tmp[i]=0;}
    if((iperd==2)||(iperd==3)){
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy43[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz43[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx43[ibond]*intra_scr_fy1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        pvten_tmp[1] += intra_scr_p11[ibond];
        pvten_tmp[5] += intra_scr_p22[ibond];
        pvten_tmp[9] += intra_scr_p33[ibond];
        pvten_tmp[2] += intra_scr_p12[ibond];
      }/*endfor*/
      pvten_tmp[4] = pvten_tmp[2];
      if(iperd==3){
        for(ibond=1;ibond <= nnow; ++ibond) {           
          intra_scr_p13[ibond] = intra_scr_dx43[ibond]*intra_scr_fz1[ibond];
          intra_scr_p23[ibond] = intra_scr_dy43[ibond]*intra_scr_fz1[ibond];
        }/*endfor*/
        for(ibond=1;ibond <= nnow; ++ibond) {
          pvten_tmp[3] += intra_scr_p13[ibond];
          pvten_tmp[6] += intra_scr_p23[ibond];
        }/*endfor*/
        pvten_tmp[7] = pvten_tmp[3];
        pvten_tmp[8] = pvten_tmp[6];
      }/*endif*/     

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

      if(np_forc > 1 ){
        for(i=1;i<=9;i++){pvten_tmp2[i] = pvten_tmp[i];}
        Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
      }/*endif*/

      for(i=1;i<=9;i++){pvten_inc[i]+=pvten_tmp[i];}    
      constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
    }else{
      for(i=1;i<=9;i++){pvten_tmp[i]=0;}
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx43[ibond]*intra_scr_fx1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        pvten_tmp[1] += intra_scr_p11[ibond];}
      for(i=1;i<=9;i++){pvten_inc[i]+=pvten_tmp[i];}    
    }/*endif*/
/*----------------------------------------------------------------------*/
/*  G) Add the force to the velocities then the positions               */
    for(ibond=1;ibond <= nnow; ++ibond) {
      temp1 = (dt2*intra_scr_m1[ibond]);
      temp2 = (dt2*intra_scr_m2[ibond]);
      intra_scr_fx1[ibond] *=temp1;
      intra_scr_fy1[ibond] *=temp1;
      intra_scr_fz1[ibond] *=temp1;
      intra_scr_fx2[ibond] *=temp2;
      intra_scr_fy2[ibond] *=temp2;
      intra_scr_fz2[ibond] *=temp2;
    }/*endfor*/

    igo = 0;
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
       for(i=1;i<=9;i++){
         fgmat_p[i]+= pvten_tmp[i];    
         vgmat[i]  += (pvten_tmp[i]*roll_scg*dt2/mass_hm);
       }/*endfor*/
    }/*endif*/

    igo = 0;
    if(igo==0){
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx3[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy3[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz3[iboff];
     }
    }else{
#ifndef NO_PRAGMA
#pragma IVDEP
#endif
     for(iboff=1,ibond=ibig;ibond <= iend; ++ibond,++iboff) {       
      ktemp = bond_j1_con[ibond];
      clatoms_x[ktemp] +=  intra_scr_fx3[iboff];
      clatoms_y[ktemp] +=  intra_scr_fy3[iboff];
      clatoms_z[ktemp] +=  intra_scr_fz3[iboff];
     }
    }/*endif*/

    igo = 0;
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
    for(i=1;i<=9;i++){pvten_tmp[i]=pvten_inc[i];}
    constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
    for(i=1;i<=9;i++){pvten_tmp[i]-=pvten_inc_old[i];}
    for(i=1;i<=9;i++){      
      fgmat_p[i] += pvten_tmp[i];
      vgmat[i]   += (pvten_tmp[i]*roll_scg*dt2/mass_hm);
    }/*endfor*/
  }/*endif*/ 

/*==========================================================================*/
/* III) Flip back                                                           */

  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* IV) consistently re-calculate everything                                 */
  if(iter!=1){
/*----------------------------------------------------------------------*/
/* A) new atom positions */
    for(i=1;i<=9;i++){vtemps[i] = vgmat[i];}
    diag33(vtemps,veig,veigv,fv1,fv2);
    for(i=1;i<=3;i++){
      aa        = exp(dt2*(veig)[i]);
      vexpdt[i] = aa*aa;
      arg2      = (veig[i]*dt2)*(veig[i]*dt2);
      poly      = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
      vsindt[i] = aa*poly;
    }/*endfor*/
    for(i=1;i<=3;i++){
      joff = (i-1)*3 ;
      for(j=1;j<=3;j++){ 
        vtempx[j+joff] = veigv[j+joff]*(vexpdt)[i];
        vtempv[j+joff] = veigv[j+joff]*(vsindt)[i];
      }/*endfor*/
    }/*endfor*/
    n = 3;
    matmul_t2(veigv,vtempx,roll_mtx,n);
    matmul_t2(veigv,vtempv,roll_mtv,n);
    for(ipart=myatm_start;ipart<=(myatm_end); ++ipart) {
      tempx  =  clatoms_xold[ipart]*roll_mtx[1]
               +clatoms_yold[ipart]*roll_mtx[2]
               +clatoms_zold[ipart]*roll_mtx[3];
      tempy  =  clatoms_xold[ipart]*roll_mtx[4]
               +clatoms_yold[ipart]*roll_mtx[5]
               +clatoms_zold[ipart]*roll_mtx[6];
      tempz  = clatoms_xold[ipart]*roll_mtx[7]
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
/*----------------------------------------------------------------------*/
/* B) get the new matrix of cell parameters and their inverse */

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

  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_bond_roll_f(BOND *bond,
                        CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos, 
                        CELL *cell,
                        INTRA_SCR *intra_scr,PTENS *ptens,
                        PAR_RAHMAN *par_rahman, double dt,int iter,
                        CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*   WARNING DANGER dt may not be timeinfo.dt so let it come in             */
/*   as a parameter                                                         */
/*==========================================================================*/
     {/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"

  int i,ibond,iii;
  int ktemp,iboff;
  int ibig,ibig1,ioff,iend,nnow;
  int nlen_now;
  int igo;

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
  double *intra_scr_sc_1  = intra_scr->sc_1;
  double *intra_scr_sc_2  = intra_scr->sc_2;
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
  int nlen_use            = intra_scr->nlen;

  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *clatoms_mass    = clatoms_info->mass;
  double *clatoms_roll_sc = clatoms_info->roll_sc;

  double *bond_al_con     = bond->al_con;
  int *bond_j1_con        = bond->j1_con;
  int *bond_j2_con        = bond->j2_con;
  int ntot_use            = bond->ncon_max;
  int ntot                = bond->ncon;

  int intra_perds         = cell->intra_perds;
  int iperd               = cell->iperd;
  int hmat_cons_typ       = cell->hmat_cons_typ;
  int hmat_int_typ        = cell->hmat_int_typ;

  double *pvten_inc       = ptens->pvten_inc;
  double *pvten_inc_old   = ptens->pvten_inc_old;
  double *pvten_tmp       = ptens->pvten_tmp;
  double *pvten_tmp2      = ptens->pvten_tmp_res;

  double *roll_mtvv       = par_rahman->roll_mtvv;
  double *vgmat_g         = par_rahman->vgmat_g;
  double *fgmat_p         = par_rahman->fgmat_p;
  double roll_scg         = par_rahman->roll_scg;
  double mass_hm          = par_rahman->mass_hm;

  int np_forc             = class_comm_forc_pkg->num_proc;
  MPI_Comm comm_forc      = class_comm_forc_pkg->comm;


/*==========================================================================*/
/* I) Flip the order of the constraints on even iteration numbers            */

  if((iter % 2)==0){flip_bond_con(bond);}

/*==========================================================================*/
/* II) Loop over bonds in strips of length nlen                             */

  dt2  = dt/2.0;
  dt2i = 1.0/dt2; 
  for(ibig=1;ibig <= ntot_use;ibig += nlen_use) {

/*----------------------------------------------------------------------*/
/*  A) Calculate Offsets                                                */

    ibig1 = ibig-1;
    ioff = -ibig1;
    nlen_now = nlen_use;
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
   if(intra_perds == 1) {
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
/* ii) rolled forces */
         fxo1r =  (fxo1*roll_mtvv[1]+fyo1*roll_mtvv[2]
               +   fzo1*roll_mtvv[3])*intra_scr_sc_1[ibond];
         fyo1r =  (fxo1*roll_mtvv[4]+fyo1*roll_mtvv[5]
               +  fzo1*roll_mtvv[6])*intra_scr_sc_1[ibond];
         fzo1r =  (fxo1*roll_mtvv[7]+fyo1*roll_mtvv[8]
               +  fzo1*roll_mtvv[9])*intra_scr_sc_1[ibond];
         fxo2r =  (fxo2*roll_mtvv[1]+fyo2*roll_mtvv[2]
               +  fzo2*roll_mtvv[3])*intra_scr_sc_2[ibond];
         fyo2r =  (fxo2*roll_mtvv[4]+fyo2*roll_mtvv[5]
               +  fzo2*roll_mtvv[6])*intra_scr_sc_2[ibond];
         fzo2r =  (fxo2*roll_mtvv[7]+fyo2*roll_mtvv[8]
               +  fzo2*roll_mtvv[9])*intra_scr_sc_2[ibond];
/*rolled velocities */
         vx1r  =(intra_scr_vx1[ibond]*roll_mtvv[1]
              + intra_scr_vy1[ibond]*roll_mtvv[2]
              + intra_scr_vz1[ibond]*roll_mtvv[3])*
                intra_scr_sc_1[ibond];
         vy1r  =(intra_scr_vx1[ibond]*roll_mtvv[4]
              + intra_scr_vy1[ibond]*roll_mtvv[5]
              + intra_scr_vz1[ibond]*roll_mtvv[6])*
                intra_scr_sc_1[ibond];
         vz1r  =(intra_scr_vx1[ibond]*roll_mtvv[7]
              + intra_scr_vy1[ibond]*roll_mtvv[8]
              + intra_scr_vz1[ibond]*roll_mtvv[9])*
                intra_scr_sc_1[ibond];
         vx2r  =(intra_scr_vx2[ibond]*roll_mtvv[1]
              + intra_scr_vy2[ibond]*roll_mtvv[2]
              + intra_scr_vz2[ibond]*roll_mtvv[3])*
                intra_scr_sc_2[ibond];
         vy2r  =(intra_scr_vx2[ibond]*roll_mtvv[4]
              + intra_scr_vy2[ibond]*roll_mtvv[5]
              + intra_scr_vz2[ibond]*roll_mtvv[6])*
                intra_scr_sc_2[ibond];
         vz2r  =(intra_scr_vx2[ibond]*roll_mtvv[7]
              + intra_scr_vy2[ibond]*roll_mtvv[8]
              + intra_scr_vz2[ibond]*roll_mtvv[9])*
                intra_scr_sc_2[ibond];
/* iv) box contribution to dr/dt */
         vxg1r  = intra_scr_x1[ibond]*vgmat_g[1]
                + intra_scr_y1[ibond]*vgmat_g[2]
                + intra_scr_z1[ibond]*vgmat_g[3];
         vyg1r  = intra_scr_x1[ibond]*vgmat_g[4]
                + intra_scr_y1[ibond]*vgmat_g[5]
                + intra_scr_z1[ibond]*vgmat_g[6];
         vzg1r  = intra_scr_x1[ibond]*vgmat_g[7]
                + intra_scr_y1[ibond]*vgmat_g[8]
                + intra_scr_z1[ibond]*vgmat_g[9];
         vxg2r  = intra_scr_x2[ibond]*vgmat_g[1]
                + intra_scr_y2[ibond]*vgmat_g[2]
                + intra_scr_z2[ibond]*vgmat_g[3];
         vyg2r  = intra_scr_x2[ibond]*vgmat_g[4]
                + intra_scr_y2[ibond]*vgmat_g[5]
                + intra_scr_z2[ibond]*vgmat_g[6];
         vzg2r  = intra_scr_x2[ibond]*vgmat_g[7]
                + intra_scr_y2[ibond]*vgmat_g[8]
                + intra_scr_z2[ibond]*vgmat_g[9];
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
      }/*endfor*/
    }/*endif*/    
/*----------------------------------------------------------------------*/
/*  F) Get Pressure tensor                                              */
    if((iperd==2)||(iperd==3)){
      for(i=1;i<=9;i++){pvten_tmp[i]=0.0;}
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
        intra_scr_p22[ibond] = intra_scr_dy12[ibond]*intra_scr_fy1[ibond];
        intra_scr_p33[ibond] = intra_scr_dz12[ibond]*intra_scr_fz1[ibond];
        intra_scr_p12[ibond] = intra_scr_dx12[ibond]*intra_scr_fy1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        pvten_tmp[1] += intra_scr_p11[ibond];
        pvten_tmp[5] += intra_scr_p22[ibond];
        pvten_tmp[9] += intra_scr_p33[ibond];
        pvten_tmp[2] += intra_scr_p12[ibond];
      }/*endfor*/
      pvten_tmp[4] = pvten_tmp[2];
      if(iperd==3){
        for(ibond=1;ibond <= nnow; ++ibond) {           
          intra_scr_p13[ibond] = intra_scr_dx12[ibond]*intra_scr_fz1[ibond];
          intra_scr_p23[ibond] = intra_scr_dy12[ibond]*intra_scr_fz1[ibond];
        }/*endfor*/
        for(ibond=1;ibond <= nnow; ++ibond) {
          pvten_tmp[3] += intra_scr_p13[ibond];
          pvten_tmp[6] += intra_scr_p23[ibond];
        }/*endfor*/
        pvten_tmp[7] = pvten_tmp[3];
        pvten_tmp[8] = pvten_tmp[6];
      }/*endif*/

/*=======================================================================*/
/*  IV)Allreduce pvten_tmp     */

      if(np_forc > 1 ){
       for(i=1;i<=9;i++){pvten_tmp2[i] = pvten_tmp[i];}
       Allreduce(&(pvten_tmp2[1]), &(pvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
      }/*endif*/
      for(i=1;i<=9;i++){pvten_inc[i]+=pvten_tmp[i];}     
      constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
    }else{
      for(i=1;i<=9;i++){pvten_tmp[i]=0;}
      for(ibond=1;ibond <= nnow; ++ibond) {           
        intra_scr_p11[ibond] = intra_scr_dx12[ibond]*intra_scr_fx1[ibond];
      }/*endfor*/
      for(ibond=1;ibond <= nnow; ++ibond) {
        pvten_tmp[1] += intra_scr_p11[ibond];
      }/*endfor*/
      for(i=1;i<=9;i++){pvten_inc[i]+=pvten_tmp[i];}     
    }/*endif*/
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
    }/*endfor*/

    igo = 0;
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
      for(i=1;i<=9;i++){      
        fgmat_p[i] += pvten_tmp[i];
        vgmat_g[i] += (pvten_tmp[i]*roll_scg*dt2/mass_hm);
      }/*endfor*/
    }/*endif*/    
  }/* endfor ibig */
/*------------------------------------------------------------------------*/
/* H) Initial correction: must be done outside ibig loop                  */
/*                        Necessary in Rattle                             */
  if(iter==1){
    for(i=1;i<=9;i++){pvten_tmp[i]=pvten_inc[i];}
    constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);
    for(i=1;i<=9;i++){pvten_tmp[i]-=pvten_inc_old[i];}
    for(i=1;i<=9;i++){      
      fgmat_p[i] += pvten_tmp[i];
      vgmat_g[i] += (pvten_tmp[i]*roll_scg*dt2/mass_hm);
    }/*endfor*/
  }/*endif*/ 

/*==========================================================================*/
/* III) Flip back                                                            */

  if((iter % 2)==0){flip_bond_con(bond);}

/*--------------------------------------------------------------------------*/
    }/*end routine */
/*==========================================================================*/





