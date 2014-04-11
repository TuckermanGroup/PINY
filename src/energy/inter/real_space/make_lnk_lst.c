/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: make_link.c                                  */

/* Construct the link list                                                  */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

void lnk_assign_atm_cell_full(double *,double *,
                              double *,int,int,
                              int,double *,
                              list_int *, list_int *,
                              int);

void lnk_assign_atm_cell_root(double *,double *,
                              double *,int ,int ,
                              int ,double *,
                              list_int *, list_int *,
                              int , int *,int );

/*==========================================================================*/

void make_lnk_lst(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                NBR_LIST *nbr_list,CELL *cell, FOR_SCR *for_scr,
                INTRA_SCR *intra_scr,TIMEINFO *timeinfo,EXCL *excl,
                INTERACT *interact,int error_check_on)
               

/*==========================================================================*/
{/*Begin Routine */
/*==========================================================================*/
/*                 Local variables                                         */

  int iii,iflag,iflag2,iflag3;
  int i1,i2, i3,i4;
  int lnk_list_dim_now, i;
  int idx, idy, idz;
  int ncell_tot;
  int nmax,natm_mall,natm_use;
  double dx, dy, dz,rat_vol,dskin;
  double sx, sy, sz;
  double xcm, ycm, zcm, vol_lnk;
  double rexl_max_old;
  double dspread;
  int natm_cell_max;
  int ncell_a,ncell_b,ncell_c;
  int ncell_ap,ncell_ab,ncell_abp,ncell_abc,ncell_abci;
  int ia,ib,ic,iatm,ioff;
/* Local pointers */
  int *nbr_list_natm_cell; /* assigned below after it is allocated */
  list_int *nbr_list_lnk_list;  /* assigned below after it is allocated */
  list_int *for_scr_i_lnk;      /* assigned below after it is allocated */
  list_int *for_scr_j_lnk;      /* assigned below after it is allocated */
  double *nbr_list_hmati_lnk;/* assigned below after it is allocated */
  int *for_scr_iexcl         = for_scr->iexcl;
  double *clatoms_x          = clatoms_pos->x;
  double *clatoms_y          = clatoms_pos->y;
  double *clatoms_z          = clatoms_pos->z;
  int natm_tot = clatoms_info->natm_tot;
  int brnch_root_list_opt = nbr_list->brnch_root_list_opt;
  int nroot_tot = nbr_list->brnch_root.nroot_tot;
  int *root_atm_list = nbr_list->brnch_root.root_atm_list;

/*==========================================================================*/
/* I) Get the cell shape and the tot number of cells                        */
/*    and the max exclusion distance                                        */


  make_hmat_lnk(clatoms_info,clatoms_pos,
                cell,nbr_list,&vol_lnk);
  iflag = 0;
  if((excl->nlst)!=0){
      rexl_max_old = (nbr_list->lnklist.rexl_max);
      max_cut_excl(clatoms_info,clatoms_pos,
                   cell,nbr_list,excl);
      if(rexl_max_old<nbr_list->lnklist.rexl_max){
        iflag=1;nbr_list->lnklist.rexl_max+=nbr_list->lnklist.lnk_excl_safe;
      }else{
        iflag=0;nbr_list->lnklist.rexl_max = rexl_max_old;
      /*endif*/}
  /*endif*/}
  iflag2 = 0;
  if(clatoms_info->pi_beads>1){
   if(interact->spread_now>interact->spread_lnk){
    dspread = 1.1*interact->spread_now - interact->spread_lnk;
    interact->spread_lnk = 1.1*interact->spread_now;
    nbr_list->lnklist.rcut_max_res += dspread;
    nbr_list->lnklist.rcut_max     += dspread;
    iflag2 = 1;
   }
  }
  iflag3 = 0;
  if(interact->brnch_root_skin>interact->brnch_root_skin_lnk){
   dskin = interact->brnch_root_skin-interact->brnch_root_skin_lnk;
   interact->brnch_root_skin_lnk   = interact->brnch_root_skin;
   nbr_list->lnklist.rcut_max_res += dskin;
   nbr_list->lnklist.rcut_max     += dskin;
   iflag3 = 1;
  }

/*==========================================================================*/
/* II)The first time make the map and malloc up some memory                 */
/*    Therafter, if volume changes more then 5%                             */
/*    or if the maximum exlcusion distances changes then redo map           */
/*    and realloc any extra memory necessary.                               */

  rat_vol = 1.0;
  if(nbr_list->lnklist.ilnk_init==0){
   rat_vol = (vol_lnk/nbr_list->lnklist.vol_lnk_map);
   if(rat_vol<1.0){rat_vol=(nbr_list->lnklist.vol_lnk_map/vol_lnk);} 
  }/*endif*/
  if((nbr_list->lnklist.ilnk_init==1) || (iflag==1) || (iflag2==1)||
     (iflag3==1) || (nbr_list->lnklist.lnk_vol_safe > rat_vol  )){
    if((nbr_list->lnklist.ilnk_init==0)&&
       (nbr_list->lnklist.vol_lnk_map*nbr_list->lnklist.lnk_vol_safe 
                        < vol_lnk )){
     if(error_check_on==1){
       printf("\n$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Recalculating link list shifts due to volume change\n");
       printf("old vol=%g,present vol=%g\n",
                                  nbr_list->lnklist.vol_lnk_map,vol_lnk); 
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n\n");
     }/*endif*/
    }/*endif*/
    if((nbr_list->lnklist.ilnk_init==0)&&(iflag==1)){
     if(error_check_on==1){
       printf("\n$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Recalculating link list shifts due to a change in \n");
       printf("maximum exclusion distance %g vs %g\n",
             (nbr_list->lnklist.rexl_max-1.0),(rexl_max_old-1.0)); 
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n\n");
     }/*endif*/
    }/*endif*/
    if((nbr_list->lnklist.ilnk_init==0)&&(iflag2==1)){
     if(error_check_on==1){
       printf("\n$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Recalculating link list shifts due to a change in \n");
       printf("maximum path integral bead spread %g\n",dspread   );
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n\n");
     }/*endif*/
    }/*endif*/
    if((nbr_list->lnklist.ilnk_init==0)&&(iflag3==1)){
     if(error_check_on==1){
       printf("\n$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Recalculating link list shifts due to a change in \n");
       printf("maximum brnch root skin of %g\n",dskin   );
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n\n");
     }/*endif*/
    }/*endif*/
    nbr_list->lnklist.vol_lnk_map = vol_lnk;
    lnk_map_cnt(clatoms_info,clatoms_pos,
                nbr_list,timeinfo,cell,error_check_on);
    make_lnk_map(clatoms_info,clatoms_pos,
                 nbr_list,timeinfo,cell,excl);
  /*endif*/}


  if(nbr_list->lnklist.ilnk_init==1){
     natm_mall = natm_tot;
     if((natm_mall % 2)==0){natm_mall+=1;} 
     for_scr->i_lnk     = (list_int *) cmalloc(natm_mall*sizeof(list_int))-1;
     for_scr->j_lnk     = (list_int *) cmalloc(natm_mall*sizeof(list_int))-1;
  }/*endif*/
  for_scr_i_lnk = for_scr->i_lnk;
  for_scr_j_lnk = for_scr->j_lnk;
  nbr_list_natm_cell = nbr_list->lnklist.natm_cell;
  ncell_tot = (nbr_list->lnklist.ncell_a)*(nbr_list->lnklist.ncell_b)*
                                          (nbr_list->lnklist.ncell_c);
  nbr_list_hmati_lnk = nbr_list->lnklist.hmati_lnk;

/*==========================================================================*/
/* III)put particles in cell                                                */

  ncell_a  = nbr_list->lnklist.ncell_a;
  ncell_b  = nbr_list->lnklist.ncell_b;
  ncell_c  = nbr_list->lnklist.ncell_c;

 if(brnch_root_list_opt == 0){
  lnk_assign_atm_cell_full(clatoms_x,clatoms_y,clatoms_z,ncell_a,ncell_b,
                           ncell_c,nbr_list_hmati_lnk,for_scr_i_lnk,
                           for_scr_j_lnk,natm_tot);

  natm_use = natm_tot;
 }/*endif*/

 if(brnch_root_list_opt > 0){
  lnk_assign_atm_cell_root(clatoms_x,clatoms_y,clatoms_z,ncell_a,ncell_b,
                           ncell_c,nbr_list_hmati_lnk,for_scr_i_lnk,
                           for_scr_j_lnk,nroot_tot,root_atm_list,natm_tot);

  natm_use = nroot_tot;
 }/*endif*/


/*==========================================================================*/
/* III)Count up the particles in each cell                                  */

  for (i = 1; i <= ncell_tot; ++i) {
    nbr_list_natm_cell[i] = 0;
  /*endfor*/}
  for (i = 1; i <= natm_use; ++i) {
    ++nbr_list_natm_cell[(for_scr_i_lnk[i])];
  /*endfor*/}
  
/*==========================================================================*/
/* IV) Find the maximum number of particles in any one cell                 */

  natm_cell_max = 0;
  for (i = 1; i <= ncell_tot; ++i) {
    i2 = natm_cell_max;
    i3 = nbr_list_natm_cell[i];
    natm_cell_max = MAX(i2,i3);
  /*endfor*/}
  nbr_list->lnklist.natm_cell_max = natm_cell_max;

/*==========================================================================*/
/* V) Using the maximum number of particles in any one cell                 */
/*    determine the size of the new link list.                              */
/*    If it is the first time through, malloc the memory                    */
/*    If the new link list requires too much memory get some more.          */

  lnk_list_dim_now = (nbr_list->lnklist.natm_cell_max)*ncell_tot*16;
  if(nbr_list->lnklist.ilnk_init==1){
     nmax = (nbr_list->lnklist.natm_cell_max)*(nbr_list->verlist.mem_safe);
     if(nmax==nbr_list->lnklist.natm_cell_max){nmax+=1;} 
     lnk_list_dim_now   = nmax*ncell_tot*16;
     if((lnk_list_dim_now % 2)==0){lnk_list_dim_now+=1;}
     nbr_list->lnklist.lnk_list_dim  = lnk_list_dim_now;
     nbr_list->lnklist.lnk_list = 
                    (list_int *) cmalloc(lnk_list_dim_now*sizeof(list_int))-1;
    /*endif*/}
    if(lnk_list_dim_now > nbr_list->lnklist.lnk_list_dim){
       nmax = nbr_list->lnklist.natm_cell_max*nbr_list->verlist.mem_safe;
       if(nmax==nbr_list->lnklist.natm_cell_max){nmax+=1;} 
     if(error_check_on==1){
       printf("\n=============================================\n"); 
       printf("Allocating more link list memory: old %d,new %d \n",
                    nbr_list->lnklist.lnk_list_dim,
                                          nmax*ncell_tot);
     }/*endif*/
       lnk_list_dim_now   = nmax*ncell_tot*16;
       if((lnk_list_dim_now % 2)==0){lnk_list_dim_now+=1;}
       nbr_list->lnklist.lnk_list_dim  = lnk_list_dim_now;
       nbr_list->lnklist.lnk_list = 
          (list_int *)crealloc(&(nbr_list->lnklist.lnk_list)[1],
                           (nbr_list->lnklist.lnk_list_dim)*sizeof(list_int))-1;
     if(error_check_on==1){
       printf("=============================================\n\n"); 
     }/*endif*/
  /*endif*/}
  nbr_list_lnk_list = nbr_list->lnklist.lnk_list;

/*==========================================================================*/
/* VI) rezero the # natoms in cell list and fill the link list              */

  for (i = 1; i <= ncell_tot; ++i) {
    nbr_list_natm_cell[i] = 0;
  /*endfor*/}
  lnk_list_dim_now = nbr_list->lnklist.natm_cell_max*ncell_tot*16;
  for (i = 1; i <= lnk_list_dim_now; ++i) {
    nbr_list_lnk_list[i]   = 0;
  /*endfor*/}

/*--------------------------------------------------------------------------*/
/* Filling the link list:                                                   */
/*   The list is filled so that a full 4 index beast is created             */
/*   ncell_a,ncell_b,ncell_c,natm_cell_max.                                 */
/*   Zero's are placed where no are particles are present, i.e.             */
/*   natm_cell(i) < natm_cell_max. The zeros allow for efficient shifts     */
/*   to be made along the 4th dimension of the list i.e. shifting           */
/*   the atoms in each cell around.                                         */
/*--------------------------------------------------------------------------*/

 if(brnch_root_list_opt == 0){
  for (i = 1; i <= natm_tot; ++i) {
    i3 = for_scr_i_lnk[i];
    nbr_list_natm_cell[i3]++;
    i4 = for_scr_j_lnk[i] + (nbr_list_natm_cell[i3] - 1)*8*ncell_tot;
    for_scr_iexcl[i] = i4;
    nbr_list_lnk_list[i4] = i;
  /*endfor*/}
  if(natm_tot>1){lnkcell_sort(natm_tot,for_scr->iexcl);}
 }/*endif*/

 if(brnch_root_list_opt > 0){
  for (i = 1; i <= nroot_tot; ++i) {
    i3 = for_scr_i_lnk[i];
    nbr_list_natm_cell[i3]++;
    i4 = for_scr_j_lnk[i] + (nbr_list_natm_cell[i3] - 1)*8*ncell_tot;
    for_scr_iexcl[i] = i4;
    nbr_list_lnk_list[i4] = root_atm_list[i];
  /*endfor*/}
  if(natm_use>1){lnkcell_sort(natm_use,for_scr->iexcl);}
 }/*endif*/

/*--------------------------------------------------------------------------*/
/* Perform 16 full shifts: All possible shifts are now stored    */
/*--------------------------------------------------------------------------*/

 if(brnch_root_list_opt == 0){
  ncell_ab   = 2*ncell_a*ncell_b;
  ncell_abc  = 4*ncell_a*ncell_b*ncell_c;
  ncell_abci = 8*ncell_a*ncell_b*ncell_c*natm_cell_max;
  for(ia=0;ia<=1;ia++){
  for(ib=0;ib<=1;ib++){
  for(ic=0;ic<=1;ic++){
  for(iatm=0;iatm<=1;iatm++){
  ioff = ia*ncell_a+ib*ncell_ab + ic*ncell_abc + iatm*ncell_abci;
  for(i=1;i<=natm_tot;i++){
    i1 = for_scr_iexcl[i];
    i2 = i1+ioff;
    nbr_list_lnk_list[i2] = nbr_list_lnk_list[i1];
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
 }/*endif*/

 if(brnch_root_list_opt > 0){
  ncell_ab   = 2*ncell_a*ncell_b;
  ncell_abc  = 4*ncell_a*ncell_b*ncell_c;
  ncell_abci = 8*ncell_a*ncell_b*ncell_c*natm_cell_max;
  for(ia=0;ia<=1;ia++){
  for(ib=0;ib<=1;ib++){
  for(ic=0;ic<=1;ic++){
  for(iatm=0;iatm<=1;iatm++){
  ioff = ia*ncell_a+ib*ncell_ab + ic*ncell_abc + iatm*ncell_abci;
  for(i=1;i<=nroot_tot;i++){
    i1 = for_scr_iexcl[i];
    i2 = i1+ioff;
    nbr_list_lnk_list[i2] = nbr_list_lnk_list[i1];
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
  }/*endfor*/
 }/*endif*/

/*--------------------------------------------------------------------------*/
}/* make_link */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnkcell_sort(int n, int *index)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m-=1;
      rindex = index[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      index[ir]=index[1];
      ir-=1;
      if(ir==1){
	index[1]=rindex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)]))j+=1;
      /*    b)demote */
      if(rindex<index[j]){
	index[i]=index[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_assign_atm_cell_full(double *clatoms_x,double *clatoms_y,
                              double *clatoms_z,int ncell_a,int ncell_b,
                              int ncell_c,double *nbr_list_hmati_lnk,
                              list_int *for_scr_i_lnk,list_int *for_scr_j_lnk,
                              int natm_tot)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,ncell_ab,ncell_ap,ncell_abp;
  int idx,idy,idz;
  double xcm,ycm,zcm,dx,dy,dz,sx,sy,sz;

  xcm = 0.; ycm = 0.; zcm = 0.;
  for (i = 1; i <= natm_tot; ++i) {
    xcm += clatoms_x[i];
    ycm += clatoms_y[i];
    zcm += clatoms_z[i];
  /*endfor*/}
  xcm /= (double) (natm_tot);
  ycm /= (double) (natm_tot);
  zcm /= (double) (natm_tot);
  ncell_ab  = ncell_a*ncell_b;
  ncell_ap  = ncell_a*2;
  ncell_abp = ncell_ab*4;
  for (i = 1; i <= natm_tot; ++i) {
    dx = clatoms_x[i] - xcm;
    dy = clatoms_y[i] - ycm;
    dz = clatoms_z[i] - zcm;
    sx = dx*nbr_list_hmati_lnk[1]
       + dy*nbr_list_hmati_lnk[4]
       + dz*nbr_list_hmati_lnk[7];
    sy = dx*nbr_list_hmati_lnk[2]
       + dy*nbr_list_hmati_lnk[5]
       + dz*nbr_list_hmati_lnk[8];
    sz = dx*nbr_list_hmati_lnk[3]
       + dy*nbr_list_hmati_lnk[6]
       + dz*nbr_list_hmati_lnk[9];
/*--------------------------------------------------------------------------*/
/* Put particles in the box (from zero to 1)                                */
/* Hmat_lnk is constructed such that this operation is OK for all iperds    */
    sx = sx - NINT(sx-0.5);
    sy = sy - NINT(sy-0.5);
    sz = sz - NINT(sz-0.5);
/*--------------------------------------------------------------------------*/
/*  Find index of cell in the box where the particle is located             */ 
/*  MOD used in case particle is exactly at the boundry                     */

    idx = ((int)(sx*ncell_a)+ncell_a) % ncell_a;
    idy = ((int)(sy*ncell_b)+ncell_b) % ncell_b;
    idz = ((int)(sz*ncell_c)+ncell_c) % ncell_c;
    for_scr_i_lnk[i] = idx + 1 + idy*ncell_a  + idz*ncell_ab;
    for_scr_j_lnk[i] = idx + 1 + idy*ncell_ap + idz*ncell_abp;
  /*endfor*/}



/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void lnk_assign_atm_cell_root(double *clatoms_x,double *clatoms_y,
                              double *clatoms_z,int ncell_a,int ncell_b,
                              int ncell_c,double *nbr_list_hmati_lnk,
                              list_int *for_scr_i_lnk, list_int *for_scr_j_lnk,
                              int nroot_tot, int *root_atm_list,int 
                              natm_tot)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,ncell_ab,ncell_ap,ncell_abp;
  int idx,idy,idz,ktemp;
  double xcm,ycm,zcm,dx,dy,dz,sx,sy,sz;

  xcm = 0.; ycm = 0.; zcm = 0.;
  for (i = 1; i <= natm_tot; ++i) {
    xcm += clatoms_x[i];
    ycm += clatoms_y[i];
    zcm += clatoms_z[i];
  /*endfor*/}
  xcm /= (double) (natm_tot);
  ycm /= (double) (natm_tot);
  zcm /= (double) (natm_tot);
  ncell_ab  = ncell_a*ncell_b;
  ncell_ap  = ncell_a*2;
  ncell_abp = ncell_ab*4;
  for (i = 1; i <= nroot_tot; ++i) {
    ktemp = root_atm_list[i];
    dx = clatoms_x[ktemp] - xcm;
    dy = clatoms_y[ktemp] - ycm;
    dz = clatoms_z[ktemp] - zcm;
    sx = dx*nbr_list_hmati_lnk[1]
       + dy*nbr_list_hmati_lnk[4]
       + dz*nbr_list_hmati_lnk[7];
    sy = dx*nbr_list_hmati_lnk[2]
       + dy*nbr_list_hmati_lnk[5]
       + dz*nbr_list_hmati_lnk[8];
    sz = dx*nbr_list_hmati_lnk[3]
       + dy*nbr_list_hmati_lnk[6]
       + dz*nbr_list_hmati_lnk[9];
/*--------------------------------------------------------------------------*/
/* Put particles in the box (from zero to 1)                                */
/* Hmat_lnk is constructed such that this operation is OK for all iperds    */
    sx = sx - NINT(sx-0.5);
    sy = sy - NINT(sy-0.5);
    sz = sz - NINT(sz-0.5);
/*--------------------------------------------------------------------------*/
/*  Find index of cell in the box where the particle is located             */ 
/*  MOD used in case particle is exactly at the boundry                     */

    idx = ((int)(sx*ncell_a)+ncell_a) % ncell_a;
    idy = ((int)(sy*ncell_b)+ncell_b) % ncell_b;
    idz = ((int)(sz*ncell_c)+ncell_c) % ncell_c;
    for_scr_i_lnk[i] = idx + 1 + idy*ncell_a  + idz*ncell_ab;
    for_scr_j_lnk[i] = idx + 1 + idy*ncell_ap + idz*ncell_abp;
  /*endfor*/}



/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


