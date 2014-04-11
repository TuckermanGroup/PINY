/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: ewald3d_both.c                                 */
/*                                                                          */
/* Performs the three-dimensional Ewald sum with RESPA                      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_recip_both(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                        CELL *cell,PTENS *ptens,
                        double alp_ewd,int nktot,
                        int kastr[],int kbstr[],int kcstr[],
                        int ibrk1[],int ibrk2[],int ibrk3[],
		        EWD_SCR *ewd_scr,double *vrecip,
		        double wght_ter,
			double wght_ter_res,int iver_get,
                        CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                        int iget_pv_real_inter, 
                        double *clus_corr_r, double *dclus_corr_r)

/*==========================================================================*/
{/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

      double falp2,vol,pivol,arg,rvol;
      double aka,akb,akc,xk,yk,zk,atemp,btemp,ctemp;
      double sumr,sumi,g2,preg,prep,tpi;
      double srx,sry,srz,six,siy,siz,temp,smag,q_sum1;
      double vrecip_long,vrecip_short;
      int i,ipart,icount,ktemp;
      int istart,iend,irem,ngo,icoff;

/* Define local pointers                                                  */
  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_q       = clatoms_info->q;
  double *clatoms_fx      = clatoms_pos->fx;
  double *clatoms_fy      = clatoms_pos->fy;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fxt     = clatoms_pos->fxt;
  double *clatoms_fyt     = clatoms_pos->fyt;
  double *clatoms_fzt     = clatoms_pos->fzt;
  double *cell_hmati      = cell->hmati;
  double *ptens_pvten     = ptens->pvten;
  double *ptens_pvten_tot = ptens->pvten_tot;
  double *ewd_scr_x       = ewd_scr->x;
  double *ewd_scr_y       = ewd_scr->y;
  double *ewd_scr_z       = ewd_scr->z;
  double *ewd_scr_q       = ewd_scr->q;
  double *ewd_scr_fx      = ewd_scr->fx;
  double *ewd_scr_fy      = ewd_scr->fy;
  double *ewd_scr_fz      = ewd_scr->fz;
  double *ewd_scr_fx2     = ewd_scr->fx2;
  double *ewd_scr_fy2     = ewd_scr->fy2;
  double *ewd_scr_fz2     = ewd_scr->fz2;
  double *ewd_scr_heli    = ewd_scr->heli;
  double *ewd_scr_helr    = ewd_scr->helr;
  double *ewd_scr_cossc   = ewd_scr->cossc;
  double *ewd_scr_sinsc   = ewd_scr->sinsc;
  int *clatoms_ichrg      = clatoms_info->ichrg;
  int clatoms_nchrg       = clatoms_info->nchrg;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  int iperd               = cell->iperd; 

/*========================================================================*/
/* I) Get some useful constants                                           */

      falp2 = 4.0*alp_ewd*alp_ewd;
      vol = getdeth(cell->hmat);
      tpi = 2.0*M_PI;
      pivol = vol/(4.0*M_PI);
      rvol  = 1.0/vol; 

/*------------------------------------------------------------------------*/
/* II) Get the charged particles                                          */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       ktemp = clatoms_ichrg[ipart];
       ewd_scr_x[ipart] = clatoms_x[ktemp];
       ewd_scr_y[ipart] = clatoms_y[ktemp];
       ewd_scr_z[ipart] = clatoms_z[ktemp];
       ewd_scr_q[ipart] = clatoms_q[ktemp];
    }/*endfor*/

/*------------------------------------------------------------------------*/
/* III) Zero the forces                                                   */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       ewd_scr_fx[ipart]  = 0.0;
       ewd_scr_fy[ipart]  = 0.0;
       ewd_scr_fz[ipart]  = 0.0;
       ewd_scr_fx2[ipart] = 0.0;
       ewd_scr_fy2[ipart] = 0.0;
       ewd_scr_fz2[ipart] = 0.0;
    }/*endfor*/

/*------------------------------------------------------------------------*/
/* IV) Find cos and sin of sc components of charged particles             */
/*  ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                           */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       ctemp = (ewd_scr_x[ipart]*cell_hmati[3]
	     +  ewd_scr_y[ipart]*cell_hmati[6]
             +  ewd_scr_z[ipart]*cell_hmati[9])*tpi;
       ewd_scr_cossc[ipart] = cos(ctemp);
       ewd_scr_sinsc[ipart] = sin(ctemp);
    }/*endfor*/

/*------------------------------------------------------------------------*/
/* V) Perform the ewald sum                                               */

    vrecip_long     = 0.0;
    vrecip_short = 0.0;

    ngo      = nktot/np_forc;
    irem     = nktot % np_forc;
    if(myid_forc<=irem){
      istart = (ngo + 1)*myid_forc + 1;
    }else{
      istart = (ngo + 1)*irem + ngo*(myid_forc - irem) + 1;
    }/*endif*/
    if(myid_forc<irem){ngo++;}
    iend = istart + ngo - 1;
    icoff = -istart+1;

    for(icount=istart;icount <= iend; ++icount) {

/*------------------------------------------------------------------------*/
/*  A) Get the k vectors                                                  */

     aka = (double) kastr[icount];
     akb = (double) kbstr[icount];
     akc = (double) kcstr[icount];
     xk = (aka*cell_hmati[1] + akb*cell_hmati[2] + akc*cell_hmati[3])*tpi;
     yk = (aka*cell_hmati[4] + akb*cell_hmati[5] + akc*cell_hmati[6])*tpi;
     zk = (aka*cell_hmati[7] + akb*cell_hmati[8] + akc*cell_hmati[9])*tpi;
     g2 = xk*xk + yk*yk + zk*zk;
     preg = exp((-g2/falp2))/(g2*pivol);
     prep = -2.0*preg*((g2/falp2)+1.0)/g2;

     if(iperd != 3){
       preg += clus_corr_r[icount+icoff]*rvol;
       prep += dclus_corr_r[icount+icoff]*rvol;
     }/*endif*/

/*------------------------------------------------------------------------*/
/*  B) If break point number one calculate the initial helpful vectors    */

     if((ibrk1[icount] == 1) || (istart==icount)) {
      for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
         atemp = ewd_scr_x[ipart]*cell_hmati[1]
	       + ewd_scr_y[ipart]*cell_hmati[4]
               + ewd_scr_z[ipart]*cell_hmati[7];
         btemp = ewd_scr_x[ipart]*cell_hmati[2]
	       + ewd_scr_y[ipart]*cell_hmati[5]
               + ewd_scr_z[ipart]*cell_hmati[8];
         ctemp = ewd_scr_x[ipart]*cell_hmati[3]
	       + ewd_scr_y[ipart]*cell_hmati[6]
               + ewd_scr_z[ipart]*cell_hmati[9];
         arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;
         ewd_scr_helr[ipart] = cos(arg);
         ewd_scr_heli[ipart] = sin(arg);
       }/* endfor */
      }/* endif*/
/*------------------------------------------------------------------------*/
/*  C) Get the structure factor                                           */

       sumr = 0.0;
       sumi = 0.0;
       for(i=1;i <= clatoms_nchrg;i++){
        sumr += ewd_scr_helr[i]*ewd_scr_q[i];
        sumi += ewd_scr_heli[i]*ewd_scr_q[i];
       }
       smag = sumr*sumr + sumi*sumi;

/*------------------------------------------------------------------------*/
/*  D) Use the stucture factor to get the potential energy and            */
/*     the pressure*volume tensor                                         */

       prep   *= smag;
       if(ibrk3[icount]==0) {
        ptens_pvten[1] += prep*xk*xk*wght_ter;
        ptens_pvten[5] += prep*yk*yk*wght_ter;
        ptens_pvten[9] += prep*zk*zk*wght_ter;
        ptens_pvten[2] += prep*xk*yk*wght_ter;
        ptens_pvten[4] += prep*xk*yk*wght_ter;
        ptens_pvten[3] += prep*xk*zk*wght_ter;
        ptens_pvten[7] += prep*xk*zk*wght_ter;
        ptens_pvten[6] += prep*yk*zk*wght_ter;
        ptens_pvten[8] += prep*yk*zk*wght_ter;

        vrecip_long += smag*preg;
       } else {
        ptens_pvten[1] += (prep*xk*xk)*wght_ter_res;
        ptens_pvten[5] += (prep*yk*yk)*wght_ter_res;
        ptens_pvten[9] += (prep*zk*zk)*wght_ter_res;
        ptens_pvten[2] += (prep*xk*yk)*wght_ter_res;
        ptens_pvten[4] += (prep*xk*yk)*wght_ter_res;
        ptens_pvten[3] += (prep*xk*zk)*wght_ter_res;
        ptens_pvten[7] += (prep*xk*zk)*wght_ter_res;
        ptens_pvten[6] += (prep*yk*zk)*wght_ter_res;
        ptens_pvten[8] += (prep*yk*zk)*wght_ter_res;
        vrecip_short += smag*preg;
      }/*endif*/
      if(iget_pv_real_inter==1){
        ptens_pvten_tot[1] += prep*xk*xk;
        ptens_pvten_tot[5] += prep*yk*yk;
        ptens_pvten_tot[9] += prep*zk*zk;
        ptens_pvten_tot[2] += prep*xk*yk;
        ptens_pvten_tot[4] += prep*xk*yk;
        ptens_pvten_tot[3] += prep*xk*zk;
        ptens_pvten_tot[7] += prep*xk*zk;
        ptens_pvten_tot[6] += prep*yk*zk;
        ptens_pvten_tot[8] += prep*yk*zk;
      }/*endif*/
/*------------------------------------------------------------------------*/
/*  E) Get the force on the particles                                     */

   sumr *= preg*2.0;
   sumi *= preg*2.0;
   if(ibrk3[icount]==0) {
    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
     srx = xk*sumr*ewd_scr_q[ipart];
     sry = yk*sumr*ewd_scr_q[ipart];
     srz = zk*sumr*ewd_scr_q[ipart];
     six = xk*sumi*ewd_scr_q[ipart];
     siy = yk*sumi*ewd_scr_q[ipart];
     siz = zk*sumi*ewd_scr_q[ipart];
     ewd_scr_fx[ipart] += srx*ewd_scr_heli[ipart] - six*ewd_scr_helr[ipart];
     ewd_scr_fy[ipart] += sry*ewd_scr_heli[ipart] - siy*ewd_scr_helr[ipart];
     ewd_scr_fz[ipart] += srz*ewd_scr_heli[ipart] - siz*ewd_scr_helr[ipart];
   }/*endfor*/
  } else {
    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
     srx = xk*sumr*ewd_scr_q[ipart];
     sry = yk*sumr*ewd_scr_q[ipart];
     srz = zk*sumr*ewd_scr_q[ipart];
     six = xk*sumi*ewd_scr_q[ipart];
     siy = yk*sumi*ewd_scr_q[ipart];
     siz = zk*sumi*ewd_scr_q[ipart];
     ewd_scr_fx2[ipart] +=srx*ewd_scr_heli[ipart] - six*ewd_scr_helr[ipart];
     ewd_scr_fy2[ipart] +=sry*ewd_scr_heli[ipart] - siy*ewd_scr_helr[ipart];
     ewd_scr_fz2[ipart] +=srz*ewd_scr_heli[ipart] - siz*ewd_scr_helr[ipart];
   }/*endfor*/
  }/*endif*/
   if(ibrk2[icount]==1) {
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      temp = ewd_scr_helr[ipart];
      ewd_scr_helr[ipart] = ewd_scr_helr[ipart]*ewd_scr_cossc[ipart] 
                           - ewd_scr_heli[ipart]*ewd_scr_sinsc[ipart];
      ewd_scr_heli[ipart] = ewd_scr_heli[ipart]*ewd_scr_cossc[ipart] 
                           + temp*ewd_scr_sinsc[ipart];
      }/*endfor*/
     }/*endif*/
/*------------------------------------------------------------------------*/
  } /* endfor icount */

   if(iperd>0 && iperd!=3) {
     q_sum1 = dsum1(clatoms_nchrg,ewd_scr_q,1);
     vrecip_short += (0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol);
   }/*endif*/

/*------------------------------------------------------------------------*/
/* VI) Put the force back                                                 */
  if(iver_get==1){
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fxt[ktemp] += ewd_scr_fx[ipart];
      clatoms_fyt[ktemp] += ewd_scr_fy[ipart];
      clatoms_fzt[ktemp] += ewd_scr_fz[ipart];
    }
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fxt[ktemp] += ewd_scr_fx2[ipart];
      clatoms_fyt[ktemp] += ewd_scr_fy2[ipart];
      clatoms_fzt[ktemp] += ewd_scr_fz2[ipart];
    } /*endfor*/
   } /*endif*/

     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ewd_scr_fx[ipart] *= wght_ter;
      ewd_scr_fy[ipart] *= wght_ter;
      ewd_scr_fz[ipart] *= wght_ter;

      ewd_scr_fx2[ipart] *= wght_ter_res;
      ewd_scr_fy2[ipart] *= wght_ter_res;
      ewd_scr_fz2[ipart] *= wght_ter_res;
    }/*endfor*/

     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fx[ktemp] += ewd_scr_fx[ipart];
      clatoms_fy[ktemp] += ewd_scr_fy[ipart];
      clatoms_fz[ktemp] += ewd_scr_fz[ipart];
    }
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fx[ktemp] += ewd_scr_fx2[ipart];
      clatoms_fy[ktemp] += ewd_scr_fy2[ipart];
      clatoms_fz[ktemp] += ewd_scr_fz2[ipart];
    } /*endfor*/
/*------------------------------------------------------------------------*/
/* VII) Get the diagonal piece of the pressure tensor                     */

    ptens_pvten[1] += (vrecip_long)*wght_ter + (vrecip_short)*wght_ter_res;
    ptens_pvten[5] += (vrecip_long)*wght_ter + (vrecip_short)*wght_ter_res;
    ptens_pvten[9] += (vrecip_long)*wght_ter + (vrecip_short)*wght_ter_res;

    ptens_pvten_tot[1] += vrecip_long + vrecip_short;
    ptens_pvten_tot[5] += vrecip_long + vrecip_short;
    ptens_pvten_tot[9] += vrecip_long + vrecip_short;

    *vrecip = vrecip_long + vrecip_short + *vrecip;

/*------------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/

