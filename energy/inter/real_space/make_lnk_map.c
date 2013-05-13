/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: make_map.c                                   */

/* these subroutines and force control run of the linklist stuff */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/* Map the list to the cell shape */

void make_lnk_map(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  NBR_LIST *nbr_list,TIMEINFO *timeinfo,
                  CELL *cell,EXCL *excl)

/*==========================================================================*/
{/*Begin Routine */
/*==========================================================================*/
/*                 Local variables                                         */
  int i1, i2, i3,iii;
  int ia, ib, ic;
  int na12, nb12, nc12;
  double rcell; 

/*--------------------------------------------------------------------------*/
/*I) GET THE allowed SHIFTS based on the cutoff and the                     */
/*    cell division selected in lnk_cell_cnt.                               */


/*--------------------------------------------------------------------------*/
/*  1) zero total number of shifts                                          */

  nbr_list->lnklist.nshft_lnk = 1;
  nbr_list->lnklist.ishft_a[1] = nbr_list->lnklist.ishft_b[1] = 
                                 nbr_list->lnklist.ishft_c[1] = 0;
  nbr_list->lnklist.shft_wght[1] = 1.0;
  if (excl->nlst != 0) {
    nbr_list->lnklist.iexl_chk[1] = 1;
  } else {
    nbr_list->lnklist.iexl_chk[1] = 0;
  }/*endif*/

  if(timeinfo->int_res_ter==1){
    nbr_list->lnklist.nshft_lnk_res = 1;
    nbr_list->lnklist.ishft_a_res[1] = nbr_list->lnklist.ishft_b_res[1] = 
    nbr_list->lnklist.ishft_c_res[1] = 0;
    nbr_list->lnklist.shft_wght_res[1] = 1.0;
    if (excl->nlst != 0) {
      nbr_list->lnklist.iexl_chk_res[1] = 1;
    } else {
      nbr_list->lnklist.iexl_chk_res[1] = 0;
    }/*endif*/
  /*endif*/}


/*--------------------------------------------------------------------------*/
/*  2) Number of shifts along an axis when it is an outer loop              */

  na12 = (nbr_list->lnklist.ncell_a - 1) / 2;
  nb12 = (nbr_list->lnklist.ncell_b - 1) / 2;
  nc12 = (nbr_list->lnklist.ncell_c - 1) / 2;

/*--------------------------------------------------------------------------*/
/*  3) Shift along a, no shift along c or along b (ic=ib=0)                 */

  ib = ic = 0;
  i1 = na12;
  for (ia = 1; ia <= i1; ++ia) {
    lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,nbr_list->lnklist.ncell_b,
     nbr_list->lnklist.ncell_c,&rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);

    if (rcell <= nbr_list->lnklist.rcut_max) {
      ++(nbr_list->lnklist.nshft_lnk);
      nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
      nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
      nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
      nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = 1.;

      if ((rcell < nbr_list->lnklist.rexl_max)&&( excl->nlst > 0)) {
        nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
      } else {
        nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
      }
    }

/* RESPA */
    if(timeinfo->int_res_ter==1){
    if (rcell <= nbr_list->lnklist.rcut_max_res) {
      ++(nbr_list->lnklist.nshft_lnk_res);
      nbr_list->lnklist.ishft_a_res[(nbr_list->lnklist.nshft_lnk_res)] = ia;
      nbr_list->lnklist.ishft_b_res[(nbr_list->lnklist.nshft_lnk_res)] = ib;
      nbr_list->lnklist.ishft_c_res[(nbr_list->lnklist.nshft_lnk_res)] = ic;
      nbr_list->lnklist.shft_wght_res[(nbr_list->lnklist.nshft_lnk_res)] = 1.;
      if (rcell < nbr_list->lnklist.rexl_max && excl->nlst > 0) {
        nbr_list->lnklist.iexl_chk_res[(nbr_list->lnklist.nshft_lnk_res)] = 1;
      } else {
        nbr_list->lnklist.iexl_chk_res[(nbr_list->lnklist.nshft_lnk_res)] = 0;
      }
    }
    }
  }


  if (nbr_list->lnklist.ncell_a % 2 == 0) {
    ia = na12 + 1;
    lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,nbr_list->lnklist.ncell_b,
     nbr_list->lnklist.ncell_c,&rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
    if (rcell <= nbr_list->lnklist.rcut_max) {
      ++(nbr_list->lnklist.nshft_lnk);
      nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
      nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
      nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
      nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = .5;
      if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
        nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
      } else {
        nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
      }
    }

    /* RESPA */
    if(timeinfo->int_res_ter==1){
    if (rcell <= nbr_list->lnklist.rcut_max_res) {
      ++(nbr_list->lnklist.nshft_lnk_res);
      nbr_list->lnklist.ishft_a_res[(nbr_list->lnklist.nshft_lnk_res)] = ia;
      nbr_list->lnklist.ishft_b_res[(nbr_list->lnklist.nshft_lnk_res)] = ib;
      nbr_list->lnklist.ishft_c_res[(nbr_list->lnklist.nshft_lnk_res)] = ic;
      nbr_list->lnklist.shft_wght_res[(nbr_list->lnklist.nshft_lnk_res)] = .5;
      if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
        nbr_list->lnklist.iexl_chk_res[(nbr_list->lnklist.nshft_lnk_res)] = 1;
      } else {
        nbr_list->lnklist.iexl_chk_res[(nbr_list->lnklist.nshft_lnk_res)] = 0;
      }
    }
    }
  }
  /*  4) SHIFT ALONG A, ALONG B, NO SHIFT ALONG C (IC=0) */

  ic = 0;
  i1 = nb12;
  for (ib = 1; ib <= i1; ++ib) {
    i2 = nbr_list->lnklist.ncell_a;
    for (ia = 1; ia <= i2; ++ia) {
      lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                            nbr_list->lnklist.ncell_b,
                            nbr_list->lnklist.ncell_c,
                            &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
      if (rcell <= nbr_list->lnklist.rcut_max) {
        ++(nbr_list->lnklist.nshft_lnk);
        nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
        nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
        nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
        nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = 1.;
        if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
          nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
        } else {
          nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
        }
      }
      /* RESPA */
      if(timeinfo->int_res_ter==1){
      if (rcell <= nbr_list->lnklist.rcut_max_res) {
        ++(nbr_list->lnklist.nshft_lnk_res);
        nbr_list->lnklist.ishft_a_res[(nbr_list->lnklist.nshft_lnk_res)] = ia;
        nbr_list->lnklist.ishft_b_res[(nbr_list->lnklist.nshft_lnk_res)] = ib;
        nbr_list->lnklist.ishft_c_res[(nbr_list->lnklist.nshft_lnk_res)] = ic;
        nbr_list->lnklist.shft_wght_res[
                                (nbr_list->lnklist.nshft_lnk_res)] = 1.;
        if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
          nbr_list->lnklist.iexl_chk_res[
                                (nbr_list->lnklist.nshft_lnk_res)] = 1;
        } else {
          nbr_list->lnklist.iexl_chk_res[
                                (nbr_list->lnklist.nshft_lnk_res)] = 0;
        }
      }
      } 
    }
  }


  if (nbr_list->lnklist.ncell_b % 2 == 0) {
    ib = nb12 + 1;
    i1 = nbr_list->lnklist.ncell_a;
    for (ia = 1; ia <= i1; ++ia) {
      lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                            nbr_list->lnklist.ncell_b,
                            nbr_list->lnklist.ncell_c,&rcell,
                            nbr_list->lnklist.hmat_lnk,cell->iperd);
      if (rcell <= nbr_list->lnklist.rcut_max) {
        ++(nbr_list->lnklist.nshft_lnk);
        nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
        nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
        nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
        nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = .5;
        if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
          nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
        } else {
          nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
        }
      }
      /* RESPA */
      if (rcell <= nbr_list->lnklist.rcut_max_res) {
        ++(nbr_list->lnklist.nshft_lnk_res);
        nbr_list->lnklist.ishft_a_res[(nbr_list->lnklist.nshft_lnk_res)] = ia;
        nbr_list->lnklist.ishft_b_res[(nbr_list->lnklist.nshft_lnk_res)] = ib;
        nbr_list->lnklist.ishft_c_res[(nbr_list->lnklist.nshft_lnk_res)] = ic;
        nbr_list->lnklist.shft_wght_res[
                          (nbr_list->lnklist.nshft_lnk_res)] = .5;
        if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
          nbr_list->lnklist.iexl_chk_res[
                          (nbr_list->lnklist.nshft_lnk_res)] = 1;
        } else {
          nbr_list->lnklist.iexl_chk_res[
                          (nbr_list->lnklist.nshft_lnk_res)] = 0;
        }
      }
    }
  }
  /*  5) SHIFT ALONG A, ALONG B and ALONG C */
  i1 = nc12;
  for (ic = 1; ic <= i1; ++ic) {
    i2 = nbr_list->lnklist.ncell_b;
    for (ib = 1; ib <= i2; ++ib) {
      i3 = nbr_list->lnklist.ncell_a;
      for (ia = 1; ia <= i3; ++ia) {
        lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                     nbr_list->lnklist.ncell_b,
                     nbr_list->lnklist.ncell_c,
                     &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
        if (rcell <= nbr_list->lnklist.rcut_max) {
          ++(nbr_list->lnklist.nshft_lnk);
          nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
          nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
          nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
          nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = 1.;
          if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
            nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
          } else {
            nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
          }
        }
        /*     RESPA */
        if(timeinfo->int_res_ter==1){
        if (rcell <= nbr_list->lnklist.rcut_max_res) {
          ++(nbr_list->lnklist.nshft_lnk_res);
          nbr_list->lnklist.ishft_a_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = ia;
          nbr_list->lnklist.ishft_b_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = ib;
          nbr_list->lnklist.ishft_c_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = ic;
          nbr_list->lnklist.shft_wght_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = 1.;
          if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
            nbr_list->lnklist.iexl_chk_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = 1;
          } else {
            nbr_list->lnklist.iexl_chk_res[
                        (nbr_list->lnklist.nshft_lnk_res)] = 0;
          }
        }
        }
      }
    }
  }


  if (nbr_list->lnklist.ncell_c % 2 == 0) {
    ic = nc12 + 1;
    i1 = nbr_list->lnklist.ncell_b;
    for (ib = 1; ib <= i1; ++ib) {
      i2 = nbr_list->lnklist.ncell_a;
      for (ia = 1; ia <= i2; ++ia) {
        lnk_cell_dis(ia,ib,ic, nbr_list->lnklist.ncell_a, 
                               nbr_list->lnklist.ncell_b, 
                               nbr_list->lnklist.ncell_c,
                     &rcell, nbr_list->lnklist.hmat_lnk,cell->iperd);
        if (rcell <= nbr_list->lnklist.rcut_max) {
          ++(nbr_list->lnklist.nshft_lnk);
          nbr_list->lnklist.ishft_a[(nbr_list->lnklist.nshft_lnk)] = ia;
          nbr_list->lnklist.ishft_b[(nbr_list->lnklist.nshft_lnk)] = ib;
          nbr_list->lnklist.ishft_c[(nbr_list->lnklist.nshft_lnk)] = ic;
          nbr_list->lnklist.shft_wght[(nbr_list->lnklist.nshft_lnk)] = .5;
          if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
            nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 1;
          } else {
            nbr_list->lnklist.iexl_chk[(nbr_list->lnklist.nshft_lnk)] = 0;
          }
        }
        /*     RESPA */
        if(timeinfo->int_res_ter==1){
        if (rcell <= nbr_list->lnklist.rcut_max_res) {
          ++(nbr_list->lnklist.nshft_lnk_res);
          nbr_list->lnklist.ishft_a_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = ia;
          nbr_list->lnklist.ishft_b_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = ib;
          nbr_list->lnklist.ishft_c_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = ic;
          nbr_list->lnklist.shft_wght_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = .5;
          if ((rcell < nbr_list->lnklist.rexl_max) && (excl->nlst > 0)) {
            nbr_list->lnklist.iexl_chk_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = 1;
          } else {
            nbr_list->lnklist.iexl_chk_res[
                      (nbr_list->lnklist.nshft_lnk_res)] = 0;
          }
        }
        }
      }
    }
  }

} /* make_map */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_map_cnt(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                 NBR_LIST *nbr_list,TIMEINFO *timeinfo,
                 CELL *cell,int error_check_on)

/*==========================================================================*/
{/*Begin Routine */
/*==========================================================================*/
/*                 Local variables                                         */
  int i1, i2, i3,i;
  int ia, ib, ic;
  int na12, nb12, nc12;
  double d1, d2, d3;
  double al, bl, cl;
  double ccell,rcell;
  double vol_lnk_map,fcell;
  int nshft_lnk_old,nshft_lnk_res_old,iii;
  int ncell_tot_old,ncell_tot,nshft_res,nshft;
  int ncell_tot_mall,nshft_mall;

/*==========================================================================*/
/* I)  SET UP CELL                                                          */

/*--------------------------------------------------------------------------*/
/*  0) get cell shape                                                       */

  make_hmat_lnk(clatoms_info,clatoms_pos,
                cell,nbr_list,&vol_lnk_map);
  nbr_list->lnklist.vol_lnk_map = vol_lnk_map;
/*--------------------------------------------------------------------------*/
/*  1)get cell edges                                                        */
  d1 = nbr_list->lnklist.hmat_lnk[1]; 
  d2 = nbr_list->lnklist.hmat_lnk[4]; 
  d3 = nbr_list->lnklist.hmat_lnk[7];
  al = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  d1 = nbr_list->lnklist.hmat_lnk[2]; 
  d2 = nbr_list->lnklist.hmat_lnk[5]; 
  d3 = nbr_list->lnklist.hmat_lnk[8];
  bl = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
  d1 = nbr_list->lnklist.hmat_lnk[3]; 
  d2 = nbr_list->lnklist.hmat_lnk[6]; 
  d3 = nbr_list->lnklist.hmat_lnk[9];
  cl = sqrt(d1 * d1 + d2 * d2 + d3 * d3);

  ccell = pow(nbr_list->lnklist.vol_lnk_map, (double)(1./3.));

/*--------------------------------------------------------------------------*/
/*  3) calculate the ncell_a,ncell_b,ncell_c                                */
/*     in each direction to be used. Save old stuff                         */

  if(nbr_list->lnklist.ilnk_init==0){
    ncell_tot_old =(nbr_list->lnklist.ncell_a)*
                   (nbr_list->lnklist.ncell_b)*
                   (nbr_list->lnklist.ncell_c);
    nshft_lnk_old =(nbr_list->lnklist.nshft_lnk);
    nshft_lnk_res_old =(nbr_list->lnklist.nshft_lnk_res);
  /*endif*/}
  fcell = (double) (nbr_list->lnklist.ncell_div_avg);
  i1 =  (int) (fcell *  (al / ccell));
  nbr_list->lnklist.ncell_a = MAX(i1,1);
  i1 = (int) (fcell *  (bl / ccell));
  nbr_list->lnklist.ncell_b = MAX(i1,1);
  i1 = (int) (fcell *  (cl / ccell));
  nbr_list->lnklist.ncell_c = MAX(i1,1);
  if (nbr_list->lnklist.lnk_for_odd == 1) {
    nbr_list->lnklist.ncell_a += ((nbr_list->lnklist.ncell_a - 1) % 2);
    nbr_list->lnklist.ncell_b += ((nbr_list->lnklist.ncell_b - 1) % 2);
    nbr_list->lnklist.ncell_c += ((nbr_list->lnklist.ncell_c - 1) % 2);
  /*endif*/}

/*--------------------------------------------------------------------------*/
/*II) GET THE allowed SHIFTS based on the cutoff and the                    */
/*    cell division selected                                                */

/*--------------------------------------------------------------------------*/
/*  1) zero total number of shifts                                          */
  nbr_list->lnklist.nshft_lnk = 1;
  nbr_list->lnklist.nshft_lnk_res = 1;

/*--------------------------------------------------------------------------*/
/*  2) number of shifts along an axis when it is an outer loop              */

  na12 = (nbr_list->lnklist.ncell_a - 1) / 2;
  nb12 = (nbr_list->lnklist.ncell_b - 1) / 2;
  nc12 = (nbr_list->lnklist.ncell_c - 1) / 2;
/*--------------------------------------------------------------------------*/
/*  3) shift along a, no shift along c or along b (ic=ib=0)                 */

  ib = 0;
  ic = 0;
  i1 = na12;
  for (ia = 1; ia <= i1; ++ia) {
    lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                 nbr_list->lnklist.ncell_b,
		 nbr_list->lnklist.ncell_c,&rcell,
                 nbr_list->lnklist.hmat_lnk,cell->iperd);
    if (rcell <= nbr_list->lnklist.rcut_max) {
      ++(nbr_list->lnklist.nshft_lnk);
    }
    if (rcell <= nbr_list->lnklist.rcut_max_res) {
      ++(nbr_list->lnklist.nshft_lnk_res);
    }
  }
  if (nbr_list->lnklist.ncell_a % 2 == 0) {
    ia = na12 + 1;
    lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                          nbr_list->lnklist.ncell_b,
                          nbr_list->lnklist.ncell_c,
                        &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
    if (rcell <= nbr_list->lnklist.rcut_max) {
      ++(nbr_list->lnklist.nshft_lnk);
    }
    if (rcell <= nbr_list->lnklist.rcut_max_res) {
      ++(nbr_list->lnklist.nshft_lnk_res);
    }
  }
/*--------------------------------------------------------------------------*/
/*  4) shift along a, along b, no shift along c (ic=0)                      */

  ic = 0;
  i1 = nb12;
  for (ib = 1; ib <= i1; ++ib) {
    i2 = nbr_list->lnklist.ncell_a;
    for (ia = 1; ia <= i2; ++ia) {
      lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                   nbr_list->lnklist.ncell_b,
		   nbr_list->lnklist.ncell_c,
                   &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
      if (rcell <= nbr_list->lnklist.rcut_max) {
	++(nbr_list->lnklist.nshft_lnk);
      }
      if (rcell <= nbr_list->lnklist.rcut_max_res) {
	++(nbr_list->lnklist.nshft_lnk_res);
      }
    }
  }
  if (nbr_list->lnklist.ncell_b % 2 == 0) {
    ib = nb12 + 1;
    i1 = nbr_list->lnklist.ncell_a;
    for (ia = 1; ia <= i1; ++ia) {
      lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                   nbr_list->lnklist.ncell_b,
		   nbr_list->lnklist.ncell_c,
                   &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
      if (rcell <= nbr_list->lnklist.rcut_max) {
	++(nbr_list->lnklist.nshft_lnk);
      }
      if (rcell <= nbr_list->lnklist.rcut_max_res) {
	++(nbr_list->lnklist.nshft_lnk_res);
      }
    }
  }
/*--------------------------------------------------------------------------*/
/*  5) shift along a, along b and along c                                   */

  i1 = nc12;
  for (ic = 1; ic <= i1; ++ic) {
    i2 = nbr_list->lnklist.ncell_b;
    for (ib = 1; ib <= i2; ++ib) {
      i3 = nbr_list->lnklist.ncell_a;
      for (ia = 1; ia <= i3; ++ia) {
	lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                              nbr_list->lnklist.ncell_b,
                              nbr_list->lnklist.ncell_c,
                   &rcell,nbr_list->lnklist.hmat_lnk,cell->iperd);
	if (rcell <= nbr_list->lnklist.rcut_max) {
	  ++(nbr_list->lnklist.nshft_lnk);
	}
	if (rcell <= nbr_list->lnklist.rcut_max_res) {
	  ++(nbr_list->lnklist.nshft_lnk_res);
	}
      }
    }
  }
  if (nbr_list->lnklist.ncell_c % 2 == 0) {
    ic = nc12 + 1;
    i1 = nbr_list->lnklist.ncell_b;
    for (ib = 1; ib <= i1; ++ib) {
      i2 = nbr_list->lnklist.ncell_a;
      for (ia = 1; ia <= i2; ++ia) {
	lnk_cell_dis(ia,ib,ic,nbr_list->lnklist.ncell_a,
                     nbr_list->lnklist.ncell_b,
		     nbr_list->lnklist.ncell_c,&rcell,
                     nbr_list->lnklist.hmat_lnk,cell->iperd);
	if (rcell <= nbr_list->lnklist.rcut_max) {
	  ++(nbr_list->lnklist.nshft_lnk);
	}
	if (rcell <= nbr_list->lnklist.rcut_max_res) {
	  ++(nbr_list->lnklist.nshft_lnk_res);
	}
      }
    }
  }
  if (timeinfo->int_res_ter == 0) { nbr_list->lnklist.nshft_lnk_res = 0; }

/*==========================================================================*/
/* Malloc the shift lists                                                   */

  ncell_tot =(nbr_list->lnklist.ncell_a)*
             (nbr_list->lnklist.ncell_b)*(nbr_list->lnklist.ncell_c);
  if(nbr_list->lnklist.ilnk_init==1){
    nshft = nbr_list->lnklist.nshft_lnk;
    nshft_mall = nshft;
    ncell_tot_mall = ncell_tot;
    if((ncell_tot_mall % 2)==0){ncell_tot_mall+=1;}
    if((nshft_mall % 2)==0){nshft_mall+=1;}
    nbr_list->lnklist.natm_cell= (int *)cmalloc(ncell_tot_mall*sizeof(int))-1;
    nbr_list->lnklist.ishft_a  = (int *)cmalloc(nshft_mall*sizeof(int))-1;
    nbr_list->lnklist.ishft_b  = (int *)cmalloc(nshft_mall*sizeof(int))-1;
    nbr_list->lnklist.ishft_c  = (int *)cmalloc(nshft_mall*sizeof(int))-1;
    nbr_list->lnklist.iexl_chk = (int *)cmalloc(nshft_mall*sizeof(int))-1;
    nbr_list->lnklist.shft_wght = 
                           (double *)cmalloc(nshft_mall*sizeof(double))-1;
    if(timeinfo->int_res_ter==1){
     nshft_res = nbr_list->lnklist.nshft_lnk_res;
     nshft_mall = nshft_res;
     if((nshft_mall % 2)==0){nshft_mall+=1;}
      nbr_list->lnklist.ishft_a_res  = 
                    (int *)cmalloc(nshft_mall*sizeof(int))-1;
      nbr_list->lnklist.ishft_b_res  = 
                    (int *)cmalloc(nshft_mall*sizeof(int))-1;
      nbr_list->lnklist.ishft_c_res  = 
                    (int *)cmalloc(nshft_mall*sizeof(int))-1;
      nbr_list->lnklist.iexl_chk_res = 
                    (int *)cmalloc(nshft_mall*sizeof(int))-1;
      nbr_list->lnklist.shft_wght_res = 
                    (double *) cmalloc(nshft_mall*sizeof(double))-1;
    /*endif: RESPA malloc*/}
  /*endif: the first full malloc*/}

/*==========================================================================*/
/* Realloc the shift lists                                                  */

  if(nbr_list->lnklist.ilnk_init==0){
    if(ncell_tot>ncell_tot_old){
if(error_check_on==1){
       printf("\n=============================================\n"); 
       printf("Allocating more link list memory: old %d new  %d\n",
                           ncell_tot,ncell_tot_old);
}/*endif*/
       ncell_tot_mall = ncell_tot;
       if((ncell_tot_mall % 2)==0){ncell_tot_mall+=1;}
       nbr_list->lnklist.natm_cell = 
                 (int *)crealloc(&(nbr_list->lnklist.natm_cell)[1],
                                   ncell_tot_mall*sizeof(int))-1;
if(error_check_on==1){
       printf("=============================================\n\n"); 
}/*endif*/
    /*endif:Realloc neede*/}
    if(nbr_list->lnklist.nshft_lnk>nshft_lnk_old){
if(error_check_on==1){
       printf("\n=============================================\n"); 
       printf("Allocating more link list memory: old %d new  %d\n",
                         nbr_list->lnklist.nshft_lnk,nshft_lnk_old);
}/*endif*/
       nshft = nbr_list->lnklist.nshft_lnk;
       nshft_mall = nshft;
       if((nshft_mall % 2)==0){nshft_mall+=1;}
       nbr_list->lnklist.ishft_a = 
                   (int *)crealloc(&(nbr_list->lnklist.ishft_a)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.ishft_b = 
                   (int *)crealloc(&(nbr_list->lnklist.ishft_b)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.ishft_c = 
                   (int *)crealloc(&(nbr_list->lnklist.ishft_c)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.iexl_chk = 
                   (int *)crealloc(&(nbr_list->lnklist.iexl_chk)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.shft_wght = 
                   (double *)crealloc(&(nbr_list->lnklist.shft_wght)[1],
                                              nshft_mall*sizeof(double))-1;

if(error_check_on==1){
       printf("=============================================\n\n"); 
}/*endif*/
    /*endif:Realloc needed*/}
    if(timeinfo->int_res_ter==1){
    if(nbr_list->lnklist.nshft_lnk_res>nshft_lnk_res_old){   
if(error_check_on==1){
       printf("\n=============================================\n"); 
       printf("Allocating more RESPA link list memory: old %d new  %d\n",
                      nshft_lnk_res_old,nbr_list->lnklist.nshft_lnk_res);
}/*endif*/
       nshft_res = nbr_list->lnklist.nshft_lnk_res;
       nshft_mall = nshft_res;
       if((nshft_mall % 2)==0){nshft_mall+=1;}
if(error_check_on==1){
       printf("Shifts_res: old mem=%d,new mem=%d\n",
                nbr_list->lnklist.nshft_lnk_res,nshft_lnk_res_old);
}/*endif*/
       nbr_list->lnklist.ishft_a_res  = 
                 (int *) crealloc(&(nbr_list->lnklist.ishft_a_res)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.ishft_b_res  = 
                 (int *) crealloc(&(nbr_list->lnklist.ishft_b_res)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.ishft_c_res  = 
                 (int *) crealloc(&(nbr_list->lnklist.ishft_c_res)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.iexl_chk_res  = 
                 (int *) crealloc(&(nbr_list->lnklist.iexl_chk_res)[1],
                                              nshft_mall*sizeof(int))-1;
       nbr_list->lnklist.shft_wght_res =(double *) 
          crealloc(&(nbr_list->lnklist.shft_wght_res)[1],
          nshft_mall*sizeof(double))-1;
if(error_check_on==1){
       printf("=============================================\n\n"); 
}/*endif*/
    /*endif:RESPA realloc needed*/}
    /*endif:RESPA on*/}
  /*endif:not the first time*/}

/*--------------------------------------------------------------------------*/
/* lnk_map_cnt */}
/*==========================================================================*/







