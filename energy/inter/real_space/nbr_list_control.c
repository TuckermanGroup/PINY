/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: nbr_list_control.c                             */
/*                                                                          */
/* This routine controls neighbor list generation                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

void check_root_safety(CLATOMS_POS *,NBR_LIST *,
                       INTERACT *,int );

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void nbr_list_control(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      FOR_SCR *for_scr,ATOMMAPS *atommaps,
                      CELL *cell,INTERACT *interact,
                      TIMEINFO *timeinfo, NBR_LIST *nbr_list,EXCL *excl,
                      INTRA_SCR *intra_scr,STAT_AVG *stat_avg,COMMUNICATE 
                      *communicate,int error_check_on,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
 double brnch_root_dist_max,brnch_root_dist_max_now;
 double brnch_root_skin_br,brnch_root_skin_bb,brnch_root_skin;
 double safety_br,safety_bb;
 int iflag,iflag_tmp;
 int ninter,i;
 int numprocs = communicate->np_beads;
 int brnch_root_list_opt = nbr_list->brnch_root_list_opt;
 int pi_beads=clatoms_info->pi_beads;

  ninter = (atommaps->natm_typ*(atommaps->natm_typ+1))/2;
/*==========================================================================*/
/* 0) Regenerate cut skin if you have more than one PI bead                 */
 
  if(brnch_root_list_opt>0){
    check_root_safety(clatoms_pos,nbr_list,interact,pi_beads);
  }/*endif*/

  if(clatoms_info->pi_beads>1){
   get_pimd_spread(clatoms_info,clatoms_pos,&(interact->spread_now),
                   communicate);
   get_cut_skin(interact->cutskin,
         interact->cutskin_root,
         interact->cutoff,
         interact->cutskin_res,
         interact->cutskin_root_res,
         interact->cutoff_res,
         interact->skin,
         interact->spread_now,
         interact->brnch_root_skin,
         interact->nter_typ);

  }/*endif*/

/*==========================================================================*/
/* I) Link list generation: Generate a link list.                           */

  if((nbr_list->ilnk)==1){
    make_lnk_lst(clatoms_info,&clatoms_pos[1],
                 nbr_list,cell,for_scr,intra_scr,timeinfo,excl,interact,
                 error_check_on); 
  /*endif*/}

/*==========================================================================*/
/* II) Verlet List generation: if the atoms have moved too                  */
/*                            much make a new list                          */

  if((nbr_list->iver)==1){
      check_ver_list(clatoms_info,&clatoms_pos[1],
                     nbr_list,interact,&iflag);
  
      if(numprocs > 1){
        iflag_tmp = iflag;
        Allreduce(&iflag_tmp,&iflag,1,MPI_INT,MPI_MAX,0,
                                           communicate->comm_beads);
      }/*endif*/

      if(iflag == 1){
        (stat_avg->updates) += 1.0;
        (stat_avg->itime_update) = timeinfo->itime;
         verlist_control(clatoms_info,&clatoms_pos[1],
                         nbr_list,excl,atommaps,cell,
                         intra_scr,for_scr,timeinfo,interact,error_check_on,
                         class_comm_forc_pkg);
      /*endif*/}
  /*endif*/}

/*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_root_safety(CLATOMS_POS *clatoms_pos,NBR_LIST *nbr_list,
                       INTERACT *interact,int pi_beads)

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */
 double brnch_root_dist_max,brnch_root_dist_max_now;
 double brnch_root_skin_br,brnch_root_skin_bb;
 double safety_br,safety_bb;
 double brnch_root_skin=interact->brnch_root_skin;
 double cutdif_min_brnch_brnch=nbr_list->brnch_root.cutdif_min_brnch_brnch;
 double cutdif_min_brnch_root=nbr_list->brnch_root.cutdif_min_brnch_root;

     brnch_root_dist_max = nbr_list->brnch_root.brnch_root_dist_max;
     find_max_brnch_root_dist(clatoms_pos,
                              nbr_list->brnch_root.nbrnch_tot,
                              nbr_list->brnch_root.brnch_atm_list,
                              nbr_list->brnch_root.brnch_atm_root,
                              &brnch_root_dist_max_now,
                              pi_beads);

     brnch_root_dist_max = MAX(brnch_root_dist_max,brnch_root_dist_max_now);
     nbr_list->brnch_root.brnch_root_dist_max = brnch_root_dist_max;

     safety_br = cutdif_min_brnch_root -  brnch_root_dist_max+brnch_root_skin;
     safety_bb = cutdif_min_brnch_brnch-2*brnch_root_dist_max+brnch_root_skin;
     if(safety_br<0.0 || safety_bb<0.0){
       brnch_root_skin          += (-MIN(safety_br,safety_bb) + 0.1);
       interact->brnch_root_skin = brnch_root_skin;
       get_cut_skin(interact->cutskin,
             interact->cutskin_root,
             interact->cutoff,
             interact->cutskin_res,
             interact->cutskin_root_res,
             interact->cutoff_res,
             interact->skin,
             interact->spread_now,
             interact->brnch_root_skin,
             interact->nter_typ);
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Your branch root safety is less than zero, %gA or %gA < 0\n",
               safety_br*BOHR,safety_bb*BOHR);
       printf("Increasing the branch root skin to %g A\n",brnch_root_skin*BOHR);
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");

     }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_ver_list(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                    NBR_LIST *nbr_list,INTERACT *interact,
                    int *iflag)

/*==========================================================================*/
{/*Begin Routine*/
/*==========================================================================*/
/*            Local Variable declarations                                   */

  int i,iii;
  double dx,dy,dz,dspread;
  double r2,rmax,skin2,rtmp;

/*            Local Pointer declarations                                   */

  int pi_beads      = clatoms_info->pi_beads;
  int natm_tot      = clatoms_info->natm_tot;
  double *x         = clatoms_pos->x;
  double *y         = clatoms_pos->y;
  double *z         = clatoms_pos->z;
  double *x0        = nbr_list->x0;
  double *y0        = nbr_list->y0;
  double *z0        = nbr_list->z0;
  double spread_now = interact->spread_now;
  double spread     = interact->spread;
  double skin       = interact->skin;

/*==========================================================================*/
/* I) Find the max distance between the particle positions when the list    */
/*     was calculated and now                                               */

      rmax = 0.0;
      for(i=1;i<=natm_tot;i++){
         dx = x[i]-x0[i];
         dy = y[i]-y0[i];
         dz = z[i]-z0[i];
         r2 = dx*dx+dy*dy+dz*dz;
         rtmp=rmax;
         rmax = MAX(rtmp,r2);
      }/*endfor*/
      rmax=sqrt(rmax);

/*==========================================================================*/
/* II) If the max distance is too large: a list recalculation is needed     */

      dspread = 0.0;
      if(pi_beads>1){dspread = MAX(0.0,spread_now - spread);}
      skin2    = (skin)/2.0 + dspread;

      (*iflag) = 0;
      if(rmax>skin2){ 
         (*iflag)=1;
         if(pi_beads>1){interact->spread=spread_now;}
      }/*endif*/

/*--------------------------------------------------------------------------*/
/*end routine */}
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* get_cut_skin:This subroutine gets the cut skin                           */
/*--------------------------------------------------------------------------*/
void get_cut_skin(double *cutskin,double *cutskin_root,double *cutoff,
                  double *cutskin_res,double *cutskin_root_res,
                  double *cutoff_res,double skin,double spread,
                  double brnch_root_skin,int ninter)


/*=======================================================================*/
/*            Begin subprogram:                                          */
  {/*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

  for(i=1;i<=ninter;i++){

      cutskin[i]          = cutoff[i]     + skin + spread;
      cutskin_res[i]      = cutoff_res[i] + skin + spread;
      cutskin_res[i]      = (cutskin_res[i] < cutskin[i] ? 
                             cutskin_res[i]:cutskin[i]);

      cutskin_root[i]     = cutoff[i]     + skin + spread + brnch_root_skin;
      cutskin_root_res[i] = cutoff_res[i] + skin + spread + brnch_root_skin;
      cutskin_root_res[i] = (cutskin_root_res[i] < cutskin[i] ? 
                             cutskin_root_res[i]:cutskin[i]);

      cutskin[i]          = cutskin[i]*cutskin[i];
      cutskin_res[i]      = cutskin_res[i]*cutskin_res[i];
      cutskin_root[i]     = cutskin_root[i]*cutskin_root[i];
      cutskin_root_res[i] = cutskin_root_res[i]*cutskin_root_res[i];
    }/*endfor*/

/*--------------------------------------------------------------------*/
}/* end routine */
/*==========================================================================*/









/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void find_max_brnch_root_dist(CLATOMS_POS *clatoms_pos,
                              int nbrnch_tot,int *brnch_atm_list,
                              int *brnch_atm_root,double *brnch_root_dist_max,
                              int pi_beads)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,ip;
  double dist_max;
  double *x,*y,*z;
  double x_1,x_2,y_1,y_2,z_1,z_2,r2;

  dist_max = 0.0;
  for(ip=1;ip<=pi_beads;ip++){
   x = clatoms_pos[ip].x;
   y = clatoms_pos[ip].y;
   z = clatoms_pos[ip].z;
   for(i=1;i<=nbrnch_tot;i++){
    x_1 = x[brnch_atm_list[i]];
    y_1 = y[brnch_atm_list[i]];
    z_1 = z[brnch_atm_list[i]];
    x_2 = x[brnch_atm_root[i]];
    y_2 = y[brnch_atm_root[i]];
    z_2 = z[brnch_atm_root[i]];
    r2  = (x_2 - x_1)*(x_2 - x_1) + (y_2 - y_1)*(y_2 - y_1) +
          (z_2 - z_1)*(z_2 - z_1); 
    dist_max = MAX(dist_max,r2);
   }/*endfor*/
  }/*endfor*/
  dist_max = sqrt(dist_max);
  *brnch_root_dist_max = dist_max;

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/
