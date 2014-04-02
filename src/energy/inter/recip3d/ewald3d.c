/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: ewald3d.c                                      */
/*                                                                          */
/* Performs the three-dimensional Ewald sum                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_recip3d_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_recip(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   CELL *cell,PTENS *ptens,
                   double alp_ewd,double alp_clus,int nktot,
                   int kastr[],int kbstr[],int kcstr[],
                   int ibrk1[],int ibrk2[],
                   EWD_SCR *ewd_scr,double *vrecip,
                   double wght_now, int iver_get,
                   CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                   int iget_pv_real_inter,
                   double *clus_corr_r, double *dclus_corr_r)

/*==========================================================================*/
{/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */
  double falp2,falp_clus2,vol,pivol,arg,rvol;
  double aka,akb,akc,xk,yk,zk;
  double sumr,sumi,g2,preg,prep,tpi;
  double smag,vnow;
  double q_sum1;
  double qi,helr,heli,cossc,sinsc;
  double fxtemp,fytemp,fztemp;
  double xtemp,ytemp,ztemp,fdiff;
  double vtest,fz;
  int i,ipart,jpart,icount,ktemp,ltemp,iii;
  int istart,iend,irem,ngo,icoff;

  double dz,argp,fargp,argm,fargm,area;
  double phase;

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
  double *cell_hmat       = cell->hmat;
  double *ptens_pvten     = ptens->pvten;
  double *ptens_pvten_tot = ptens->pvten_tot;
  double *ewd_scr_x       = ewd_scr->x;
  double *ewd_scr_y       = ewd_scr->y;
  double *ewd_scr_z       = ewd_scr->z;
  double *ewd_scr_q       = ewd_scr->q;
  double *ewd_scr_fx      = ewd_scr->fx;
  double *ewd_scr_fy      = ewd_scr->fy;
  double *ewd_scr_fz      = ewd_scr->fz; 
  double *ewd_scr_heli    = ewd_scr->heli;
  double *ewd_scr_helr    = ewd_scr->helr;
  double *ewd_scr_cossc   = ewd_scr->cossc;
  double *ewd_scr_sinsc   = ewd_scr->sinsc;
  double *ptens_pvten_tmp = ptens->pvten_tmp;
  int *clatoms_ichrg      = clatoms_info->ichrg;
  int clatoms_nchrg       = clatoms_info->nchrg;
  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  int iperd               = cell->iperd; 
 
/*========================================================================*/
/* I) Get some useful constants                                           */

      falp2      = 4.0*alp_ewd*alp_ewd;
      falp_clus2 = 4.0*alp_clus*alp_clus;
      vol = getdeth(cell->hmat);
      tpi = 2.0*M_PI;
      pivol = vol/(4.0*M_PI);
      rvol  = 1.0/vol; 
 
/*========================================================================*/
/* II) Get the charged particles                                          */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       ktemp = clatoms_ichrg[ipart];
       ewd_scr_x[ipart] = clatoms_x[ktemp];
       ewd_scr_y[ipart] = clatoms_y[ktemp];
       ewd_scr_z[ipart] = clatoms_z[ktemp];
       ewd_scr_q[ipart] = clatoms_q[ktemp];
    }/*endfor*/

/*========================================================================*/
/* III) Zero the forces                                                   */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       ewd_scr_fx[ipart] = 0.0;
       ewd_scr_fy[ipart] = 0.0;
       ewd_scr_fz[ipart] = 0.0;
    }/*endfor*/
    for(i=1;i<=9;i++){ptens_pvten_tmp[i]=0.0;}
    vnow = 0.0;

/*========================================================================*/
/* IV) Find cos and sin of sc components of charged particles             */
/*  ( hmati rvec = svec   r=(x,y,z) s=(a,b,c) )                           */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       xtemp = ewd_scr_x[ipart];
       ytemp = ewd_scr_y[ipart];
       ztemp = ewd_scr_z[ipart];
       ewd_scr_x[ipart] = xtemp*cell_hmati[1]
                        + ytemp*cell_hmati[4]
                        + ztemp*cell_hmati[7];
       ewd_scr_y[ipart] = xtemp*cell_hmati[2]
                        + ytemp*cell_hmati[5]
                        + ztemp*cell_hmati[8];
       ewd_scr_z[ipart] = xtemp*cell_hmati[3]
                        + ytemp*cell_hmati[6]
                        + ztemp*cell_hmati[9];
    }/*endfor*/
    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
       arg = ewd_scr_z[ipart]*tpi;
       ewd_scr_cossc[ipart] = cos(arg);
       ewd_scr_sinsc[ipart] = sin(arg);
    }/*endfor*/

/*========================================================================*/
/* V) Perform the ewald sum                                               */

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

     aka = tpi*( (double) kastr[icount] );
     akb = tpi*( (double) kbstr[icount] );
     akc = tpi*( (double) kcstr[icount] );
     xk = (aka*cell_hmati[1] + akb*cell_hmati[2] + akc*cell_hmati[3]);
     yk = (aka*cell_hmati[4] + akb*cell_hmati[5] + akc*cell_hmati[6]);
     zk = (aka*cell_hmati[7] + akb*cell_hmati[8] + akc*cell_hmati[9]);
     g2 = xk*xk + yk*yk + zk*zk;
     preg = exp((-g2/falp2))/(g2*pivol);
     if(iperd == 2 && ((kastr[icount] == 0) && (kbstr[icount] == 0))){
       phase = cos(0.5*zk*cell_hmat[9]);
       preg  += phase*(1.0-exp(-g2/falp_clus2))/(g2*pivol);
     }
     prep = -2.0*preg*((g2/falp2)+1.0)/g2;
     if(iperd == 2 && ((kastr[icount] == 0) && (kbstr[icount] == 0))){
       phase = cos(0.5*zk*cell_hmat[9]);
       prep += 2.0*phase/(g2*pivol*falp_clus2);
     }
     
     if(iperd != 3){
       preg += clus_corr_r[(icount+icoff)]*rvol;
       prep += dclus_corr_r[(icount+icoff)]*rvol;
     }/*endif*/

/*------------------------------------------------------------------------*/
/*  B) If break point number one calculate the initial helpful vectors    */

     if((ibrk1[icount] == 1) || (icount==istart)) {
      for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
         arg = aka*ewd_scr_x[ipart]+akb*ewd_scr_y[ipart]+akc*ewd_scr_z[ipart];
         ewd_scr_helr[ipart] = cos(arg);
         ewd_scr_heli[ipart] = sin(arg); 
       }  /* endfor */
      }/* endif*/
/*------------------------------------------------------------------------*/
/*  C) Get the structure factor                                           */

       sumr = 0.0;
       sumi = 0.0;
       for(i=1;i <= clatoms_nchrg;i++){
        qi = ewd_scr_q[i];
        sumr += ewd_scr_helr[i]*qi;
        sumi += ewd_scr_heli[i]*qi;
       }
       smag = sumr*sumr + sumi*sumi;

/*------------------------------------------------------------------------*/
/*  D) Use the stucture factor to get the potential energy and            */
/*     the pressure*volume tensor                                         */

       vnow += smag*preg;
       prep    *= smag;

       ptens_pvten_tmp[1] += prep*xk*xk;
       ptens_pvten_tmp[5] += prep*yk*yk;
       ptens_pvten_tmp[9] += prep*zk*zk;
       ptens_pvten_tmp[2] += prep*xk*yk;
       ptens_pvten_tmp[3] += prep*xk*zk;
       ptens_pvten_tmp[6] += prep*yk*zk;

/*------------------------------------------------------------------------*/
/*  E) Get the force on the particles                                     */

    sumr *= (preg*2.0);
    sumi *= (preg*2.0);
    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
     fdiff = ewd_scr_heli[ipart]*sumr-ewd_scr_helr[ipart]*sumi;
     fxtemp = ewd_scr_fx[ipart] + xk*fdiff;
     fytemp = ewd_scr_fy[ipart] + yk*fdiff;
     fztemp = ewd_scr_fz[ipart] + zk*fdiff;
     ewd_scr_fx[ipart] = fxtemp;
     ewd_scr_fy[ipart] = fytemp;
     ewd_scr_fz[ipart] = fztemp; 
    } 
/*------------------------------------------------------------------------*/
/*  F) If break point two, increment the helpful vectors                  */

    if(ibrk2[icount] == 1) {
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      helr = ewd_scr_helr[ipart];
      heli = ewd_scr_heli[ipart];
      cossc = ewd_scr_cossc[ipart];
      sinsc = ewd_scr_sinsc[ipart];
      ewd_scr_helr[ipart] = helr*cossc - heli*sinsc;
      ewd_scr_heli[ipart] = heli*cossc + helr*sinsc;
     }  /*endfor*/
    }/*endif*/
/*------------------------------------------------------------------------*/
  } /* endfor icount */

   if(iperd>0 && iperd!=3) {
     q_sum1 = dsum1(clatoms_nchrg,ewd_scr_q,1);
     if(iperd == 2){
       vnow += (0.5*q_sum1*q_sum1*(clus_corr_r[(ngo+1)]
             + M_PI/(alp_clus*alp_clus))*rvol);
     } else {
       vnow += (0.5*q_sum1*q_sum1*clus_corr_r[(ngo+1)]*rvol);
     }
   }/*endif*/

/*========================================================================*/
/* VI) Put the force back                                                 */

    for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ewd_scr_fx[ipart] *= ewd_scr_q[ipart];
      ewd_scr_fy[ipart] *= ewd_scr_q[ipart];
      ewd_scr_fz[ipart] *= ewd_scr_q[ipart];
    }/*endfor*/

   if(iver_get==1){ 
    if(clatoms_info->pi_beads>1){
     for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fxt[ktemp] += ewd_scr_fx[ipart];
      clatoms_fyt[ktemp] += ewd_scr_fy[ipart];
      clatoms_fzt[ktemp] += ewd_scr_fz[ipart];
     } /*endfor*/
    } /*endif*/
   } /*endif*/

   for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ewd_scr_fx[ipart] *= wght_now;
      ewd_scr_fy[ipart] *= wght_now;
      ewd_scr_fz[ipart] *= wght_now;
   }/*endfor*/

   for(ipart = 1;ipart <= clatoms_nchrg;++ipart){
      ktemp = clatoms_ichrg[ipart];
      clatoms_fx[ktemp] += ewd_scr_fx[ipart];
      clatoms_fy[ktemp] += ewd_scr_fy[ipart];
      clatoms_fz[ktemp] += ewd_scr_fz[ipart];
   } /*endfor*/

/*========================================================================*/
/* VII) Get the diagonal piece of the pressure tensor                     */

    ptens_pvten_tmp[4] = ptens_pvten_tmp[2];
    ptens_pvten_tmp[7] = ptens_pvten_tmp[3];
    ptens_pvten_tmp[8] = ptens_pvten_tmp[6];

    for(i=1;i<=9;i++){
     ptens_pvten[i]     += ptens_pvten_tmp[i]*wght_now;
    }/*endfor*/
    ptens_pvten[1] += vnow*wght_now;
    ptens_pvten[5] += vnow*wght_now;
    ptens_pvten[9] += vnow*wght_now;

    if(iget_pv_real_inter==1){
     for(i=1;i<=9;i++){
      ptens_pvten_tot[i] += ptens_pvten_tmp[i];
     }/*endfor*/
     ptens_pvten_tot[1] += vnow;
     ptens_pvten_tot[5] += vnow;
     ptens_pvten_tot[9] += vnow;
    }/*endif*/

/*======================================================================*/
/* XII) Add in the surface term the dumb way */

  if(iperd==2){
    vtest = 0.0;

    area = cell_hmat[1]*cell_hmat[5]-cell_hmat[2]*cell_hmat[4];
    for(ipart=1;ipart<=clatoms_nchrg;ipart++){
      ktemp = clatoms_ichrg[ipart];
      for(jpart=1;jpart<=clatoms_nchrg;jpart++){
        ltemp = clatoms_ichrg[jpart];
          dz = clatoms_z[ktemp]-clatoms_z[ltemp];
          dz -= cell_hmat[9]*NINT(dz/cell_hmat[9]);
          arg = alp_clus*fabs(dz-0.5*cell_hmat[9]);
          vnow  += 0.5*ewd_scr_q[ipart]*ewd_scr_q[jpart]*surf_corr(arg)/(area*alp_clus);
          vtest += 0.5*ewd_scr_q[ipart]*ewd_scr_q[jpart]*surf_corr(arg)/(area*alp_clus);
      }/* endfor jpart */
    } /* endfor ipart */


    for(ipart=1;ipart<=clatoms_nchrg;ipart++){
      ktemp = clatoms_ichrg[ipart];
      for(jpart=1;jpart<=clatoms_nchrg;jpart++){
        ltemp = clatoms_ichrg[jpart];
        if(ipart != jpart){
          dz = clatoms_z[ktemp]-clatoms_z[ltemp];
          dz -= cell_hmat[9]*NINT(dz/cell_hmat[9]);
          argm = alp_clus*(dz - 0.5*cell_hmat[9]);
          fargm = fabs(argm);
          argp = alp_clus*(dz + 0.5*cell_hmat[9]);
          fargp = fabs(argp);
          clatoms_fz[ktemp] -= 0.5*ewd_scr_q[ipart]*ewd_scr_q[jpart]*alp_clus
                             *(dsurf_corr(fargm)*(argm/fargm) 
                             + dsurf_corr(fargp)*(argp/fargp))/(area*alp_clus);
	}/* endif */
      }/* endfor */
    }/* endfor */


  }/* endif iperd */

  *vrecip += vnow;
/*------------------------------------------------------------------------*/
   }/*end routine*/
/*========================================================================*/









