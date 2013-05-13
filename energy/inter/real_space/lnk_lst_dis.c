/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: lnk_cell_dis.c                               */

/* ROUTINE to get the distance between a cell at the origen */
/* and a vector displacement (i.e. a shift) */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void lnk_cell_dis(int ia, int ib, int ic, 
		  int ncell_a, int ncell_b, int ncell_c, 
		  double *rcell, double hmat_lnk[], int iperd)

/*==========================================================================*/
{/*Begin Routine */
/*==========================================================================*/
/*                 Local variables                                         */
  double acelli, bcelli, ccelli, dx, dy, dz, dsx, dsy, dsz;
  double dsxt, dsyt, dszt;

/*==========================================================================*/
/*  ia,ib,ic are the shifts along a,b,c                                     */
/*  ncell_a,ncell_b,ncell_c are the number of divisions along a,b,c         */
/*  hmat matrix of cell parameters                                          */
/*  rcell is return with cartesian distance associated with shifts          */
/* scalar temp                                                              */
/* get inverse number of divisions along cell axis                          */
/*==========================================================================*/

  acelli = 1. / (double) ncell_a;
  bcelli = 1. / (double) ncell_b;
  ccelli = 1. / (double) ncell_c;

  /* convert shifts to a displacement along fractional coord. (the s) */

  dsx = (double) ia * acelli;
  dsy = (double) ib * bcelli;
  dsz = (double) ic * ccelli;

  /* periodic image */
/*  if (iperd >= 1) {  dsx -= NINT(dsx);}*/
/*  if (iperd >= 2) {  dsy -= NINT(dsy);}*/
/*  if (iperd >= 3) {  dsz -= NINT(dsz);}*/
    dsx -= NINT(dsx);
    dsy -= NINT(dsy);
    dsz -= NINT(dsz);

  /* subtract out one cell to account that cell that share edges can have */
  /* particles separteded by one less cell distances */

/*  dsxt = copysign(acelli,dsx); */
/*  dsyt = copysign(bcelli,dsy); */
/*  dszt = copysign(ccelli,dsz); */
  if(dsx>=0){dsxt=acelli;}else{dsxt=-acelli;}
  if(dsy>=0){dsyt=bcelli;}else{dsyt=-bcelli;}
  if(dsz>=0){dszt=ccelli;}else{dszt=-ccelli;}
  if (dsx != 0.) {dsx -= dsxt;}
  if (dsy != 0.) {dsy -= dsyt;}
  if (dsz != 0.) {dsz -= dszt;}

  /* convert to cartesian displacements */
  dx = dsx * hmat_lnk[1] + dsy * hmat_lnk[4] + dsz * hmat_lnk[7];
  dy = dsx * hmat_lnk[2] + dsy * hmat_lnk[5] + dsz * hmat_lnk[8];
  dz = dsx * hmat_lnk[3] + dsy * hmat_lnk[6] + dsz * hmat_lnk[9];

  /* get distance */
  *rcell = sqrt(dx * dx + dy * dy + dz * dz);
} /* lnk_cell_dis */
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void max_cut_part(INTERACT *interact, NBR_LIST *nbr_list)
/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*              Local Variables                                             */

  int i;
  double d1, d2;

/*==========================================================================*/
/* I) get maximum cutoff distances                                          */

  nbr_list->lnklist.rcut_max     = 0.;
  nbr_list->lnklist.rcut_max_res = 0.;
  for (i = 1; i <= interact->nter_typ; ++i) {
    d1 = nbr_list->lnklist.rcut_max;      d2 = interact->cutoff[i];
    nbr_list->lnklist.rcut_max = MAX(d1,d2);

    d1 = nbr_list->lnklist.rcut_max_res;  d2 = interact->cutoff_res[i];
    nbr_list->lnklist.rcut_max_res = MAX(d1,d2);
  }/*endfor*/

  nbr_list->lnklist.rcut_max_res += interact->spread_now*1.1;
  nbr_list->lnklist.rcut_max     += interact->spread_now*1.1;
  interact->spread_lnk = interact->spread_now*1.1;

  nbr_list->lnklist.rcut_max     += interact->brnch_root_skin;
  nbr_list->lnklist.rcut_max_res += interact->brnch_root_skin;
  interact->brnch_root_skin_lnk   = interact->brnch_root_skin;

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void max_cut_excl(CLATOMS_INFO *clatoms_info,
                  CLATOMS_POS *clatoms_pos,
                  CELL *cell,NBR_LIST *nbr_list,EXCL *excl)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*                 Local variables                                  */
  int i1, i2;
  int j, ipart;
  double r2, dx, dy, dz, sx, sy, sz,tmp;
  int ktemp,ltemp;

/*==========================================================================*/
/* Initialize the distance */

  nbr_list->lnklist.rexl_max = 0.0;

/*------------------------------------------------------------------*/
/* iperd=0                                                          */

  if (cell->iperd == 0) {
    i1 = clatoms_info->natm_tot;
    for (ipart = 1; ipart <= i1; ++ipart) {
      i2 = excl->num[ipart];
      for (j = 1; j <= i2; ++j) {
        ltemp = (j + excl->j_off[ipart]);
        ktemp = excl->j[ltemp];
	dx = clatoms_pos->x[ipart] - clatoms_pos->x[ktemp];
	dy = clatoms_pos->y[ipart] - clatoms_pos->y[ktemp];
	dz = clatoms_pos->z[ipart] - clatoms_pos->z[ktemp];
	r2 = dx * dx + dy * dy + dz * dz;
        tmp = nbr_list->lnklist.rexl_max ;
	nbr_list->lnklist.rexl_max = MAX(tmp,r2);
      }
    }
  }
/*------------------------------------------------------------------*/
/* iperd=1                                                          */

  if (cell->iperd == 2) {
    i1 = clatoms_info->natm_tot;
    for (ipart = 1; ipart <= i1; ++ipart) {
      i2 = excl->num[ipart];
      for (j = 1; j <= i2; ++j) {
        ltemp = (j + excl->j_off[ipart]);
        ktemp = excl->j[ltemp];
	dx = clatoms_pos->x[ipart] - clatoms_pos->x[ktemp];
	dy = clatoms_pos->y[ipart] - clatoms_pos->y[ktemp];
	dz = clatoms_pos->z[ipart] - clatoms_pos->z[ktemp];
	sx = dx * cell->hmati[4];
	sx -= NINT(sx);
	dx = sx * cell->hmat[1];
	r2 = dx * dx + dy * dy + dz * dz;
        tmp = nbr_list->lnklist.rexl_max ;
	nbr_list->lnklist.rexl_max = MAX(tmp,r2);
      }
    }
  }
/*------------------------------------------------------------------*/
/* iperd=2                                                          */

  if (cell->iperd == 2) {
    i1 = clatoms_info->natm_tot;
    for (ipart = 1; ipart <= i1; ++ipart) {
      i2 = excl->num[ipart];
      for (j = 1; j <= i2; ++j) {
        ltemp = (j + excl->j_off[ipart]);
        ktemp = excl->j[ltemp];
	dx = clatoms_pos->x[ipart] - clatoms_pos->x[ktemp];
	dy = clatoms_pos->y[ipart] - clatoms_pos->y[ktemp];
	dz = clatoms_pos->z[ipart] - clatoms_pos->z[ktemp];
	sx = dx * cell->hmati[1] + dy * cell->hmati[4];
	sy = dx * cell->hmati[2] + dy * cell->hmati[5];
	sx -= NINT(sx);
	sy -= NINT(sy);
	dx = sx * cell->hmat[1] + sy * cell->hmat[4];
	dy = sx * cell->hmat[2] + sy * cell->hmat[5];
	r2 = dx * dx + dy * dy + dz * dz;
        tmp = nbr_list->lnklist.rexl_max ;
	nbr_list->lnklist.rexl_max = MAX(tmp,r2);
      }
    }
  }
/*------------------------------------------------------------------*/
/* iperd=3                                                         */
  if (cell->iperd == 3) {
    i1 = clatoms_info->natm_tot;
    for (ipart = 1; ipart <= i1; ++ipart) {
      i2 = excl->num[ipart];
      for (j = 1; j <= i2; ++j) {
        ltemp = (j + excl->j_off[ipart]);
        ktemp = excl->j[ltemp];
	dx = clatoms_pos->x[ipart] - clatoms_pos->x[ktemp];
	dy = clatoms_pos->y[ipart] - clatoms_pos->y[ktemp];
	dz = clatoms_pos->z[ipart] - clatoms_pos->z[ktemp];
	sx = dx * cell->hmati[1] + dy * cell->hmati[4] + dz * cell->hmati[7];
	sy = dx * cell->hmati[2] + dy * cell->hmati[5] + dz * cell->hmati[8];
	sz = dx * cell->hmati[3] + dy * cell->hmati[6] + dz * cell->hmati[9];
	sx -= NINT(sx);
	sy -= NINT(sy);
	sz -= NINT(sz);
	dx = sx * cell->hmat[1] + sy * cell->hmat[4] + sz * cell->hmat[7];
	dy = sx * cell->hmat[2] + sy * cell->hmat[5] + sz * cell->hmat[8];
	dz = sx * cell->hmat[3] + sy * cell->hmat[6] + sz * cell->hmat[9];
	r2 = dx * dx + dy * dy + dz * dz;
        tmp = nbr_list->lnklist.rexl_max ;
	nbr_list->lnklist.rexl_max = MAX(tmp,r2);
      }
    }
  }
/*------------------------------------------------------------------*/
/* make a distance */ 

  tmp = nbr_list->lnklist.rexl_max ;
  nbr_list->lnklist.rexl_max = sqrt(tmp);

/*------------------------------------------------------------------*/
} /* max_cut_excl */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void period_lnk(int intact,double dxs[],double dys[],double dzs[],
                double hmat[], double hmati[], int iperd)

{/*begin routine*/
/*=======================================================================*/
/*            Begin subprogram:                                          */
/*=======================================================================*/
/*            Local variable declarations                                */

  int jpart;
  double sx,sy,sz;

/*=======================================================================*/
/* One dimensions */

  if((iperd) == 1) { 
    for(jpart=1;jpart <= intact; ++jpart) {
     sx   =dxs[jpart]*(hmati)[1];
     sx  -= NINT(sx);
     dxs[jpart] = sx*(hmat)[1];
   }/*endfor*/
  }/* endif */

/*=======================================================================*/
/* Two dimensions */

  if((iperd) == 2) { 
    for(jpart=1;jpart <= intact; ++jpart) {
      sx =dxs[jpart]*(hmati)[1]+dys[jpart]*(hmati)[2];
      sy =dxs[jpart]*(hmati)[3]+dys[jpart]*(hmati)[4];
      sx -= NINT(sx);
      sy -= NINT(sy);
      dxs[jpart] = sx*(hmat)[1]+sy*(hmat)[2];
      dys[jpart] = sx*(hmat)[3]+sy*(hmat)[4];
    }/* endfor */
  }/* endif */

/*======================================================================*/
/* Three dimensions */

  if((iperd) == 3) {
    for(jpart=1;jpart <= intact; ++jpart) {
      sx =dxs[jpart]*(hmati)[1] +dys[jpart]*(hmati)[4]
	+dzs[jpart]*(hmati)[7];
      sy =dxs[jpart]*(hmati)[2] +dys[jpart]*(hmati)[5]
	+dzs[jpart]*(hmati)[8];
      sz =dxs[jpart]*(hmati)[3] +dys[jpart]*(hmati)[6]
	+dzs[jpart]*(hmati)[9];
      sx -= NINT(sx);
      sy -= NINT(sy);
      sz -= NINT(sz);
      dxs[jpart] =sx*(hmat)[1]+sy*(hmat)[4]+sz*(hmat)[7];
      dys[jpart] =sx*(hmat)[2]+sy*(hmat)[5]+sz*(hmat)[8];
      dzs[jpart] =sx*(hmat)[3]+sy*(hmat)[6]+sz*(hmat)[9];
    }/* endfor */
  }/* endif */ 
  /*========================================================================*/
}/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_hmat_lnk(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   CELL *cell,NBR_LIST *nbr_list,
		   double *vol_lnk)

/*==========================================================================*/
{/*Begin Routine */
/*==========================================================================*/
/*                 Local variables                                         */
  int i1,i, iperd_now;
  double xmax, ymax, zmax;
  double dx2, dy2, dz2, xcm, ycm, zcm;

/*==========================================================================*/
/* 3D periodicity */

  for (i = 1; i <= 9; ++i) {nbr_list->lnklist.hmat_lnk[i] = cell->hmat[i];}

/*==========================================================================*/
/* Not 3D periodicity */

  if (cell->iperd < 3) {
    xcm = 0.; ycm = 0.; zcm = 0.;
    i1 = clatoms_info->natm_tot;
    for (i = 1; i <= i1; ++i) {
      xcm += clatoms_pos->x[i];
      ycm += clatoms_pos->y[i];
      zcm += clatoms_pos->z[i];
    }
    xcm /= (double) (clatoms_info->natm_tot);
    ycm /= (double) (clatoms_info->natm_tot);
    zcm /= (double) (clatoms_info->natm_tot);
    xmax = ymax = zmax = 0.;
    i1 = clatoms_info->natm_tot;
    for (i = 1; i <= i1; ++i) {
      dx2 = (clatoms_pos->x[i] - xcm) * (clatoms_pos->x[i] - xcm);
      dy2 = (clatoms_pos->y[i] - ycm) * (clatoms_pos->y[i] - ycm);
      dz2 = (clatoms_pos->z[i] - zcm) * (clatoms_pos->z[i] - ycm);
      xmax = MAX(dx2,xmax);
      ymax = MAX(dy2,ymax);
      zmax = MAX(dz2,zmax);
    }
    xmax = sqrt(xmax);
    ymax = sqrt(ymax);
    zmax = sqrt(zmax);
    if (cell->iperd <= 2) {
      nbr_list->lnklist.hmat_lnk[3] = 0.;
      nbr_list->lnklist.hmat_lnk[6] = 0.;
      nbr_list->lnklist.hmat_lnk[9] = zmax * 1.2;
    }
    if (cell->iperd <= 1) {
      nbr_list->lnklist.hmat_lnk[2] = 0.;
      nbr_list->lnklist.hmat_lnk[5] = ymax * 1.2;
      nbr_list->lnklist.hmat_lnk[8] = 0.;
    }
    if (cell->iperd <= 0) {
      nbr_list->lnklist.hmat_lnk[1] = xmax * 1.2;
      nbr_list->lnklist.hmat_lnk[4] = 0.;
      nbr_list->lnklist.hmat_lnk[7] = 0.;
    }
  }

/*==========================================================================*/
/* Get inverse matrix: used as if 3d periodicity exists */

  iperd_now = 3;
  gethinv(nbr_list->lnklist.hmat_lnk,nbr_list->lnklist.hmati_lnk, 
	  vol_lnk, iperd_now);

} /* make_hmat_lnk */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/









