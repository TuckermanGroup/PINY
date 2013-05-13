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
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_intra_con_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void shake_control(BONDED *bonded,
                   CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                   CELL *cell,PTENS *ptens,
                   STATEPOINT *statepoint,BARO *baro,PAR_RAHMAN *par_rahman,
                   STAT_AVG *stat_avg,double dt, double *tol_glob,
                   int ifirst,CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                   EWD_SCR *ewd_scr)

/*=========================================================================*/
/* WARNING DANGER dt may not be timeinfo.dt so let it come in              */
/* as a parameter                                                          */
/*=========================================================================*/
    {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */

#include "../typ_defs/typ_mask.h"
  
  int i,iii,ifirst2,init=0;
  int iter_std,iter_whol;
  int ncon_use;
  double tol_whol,tol_std,pinc_now,pinc_max;
  double aiter_21_now;
  double aiter_23_now;
  double aiter_33_now;
  double aiter_46_now;
  double aiter_43_now;

  double beta            = BOLTZ/(statepoint->t_ext); 
  MPI_Comm comm_forc     = class_comm_forc_pkg->comm;
  int np_forc            = class_comm_forc_pkg->num_proc;
  int myid_forc          = class_comm_forc_pkg->myid;

  double tolshake        = bonded->constrnt.tolshake;
  int max_iter           = bonded->constrnt.max_iter;
  int iroll              = bonded->constrnt.iroll;

  int bond_ncon          = bonded->bond.ncon;
  int bond_ncon_tot      = bonded->bond.ncon_tot;
  int num_21             = bonded->grp_bond_con.num_21;
  int num_21_tot         = bonded->grp_bond_con.num_21_tot;
  int num_23             = bonded->grp_bond_con.num_23;
  int num_23_tot         = bonded->grp_bond_con.num_23_tot;
  int num_33             = bonded->grp_bond_con.num_33;
  int num_33_tot         = bonded->grp_bond_con.num_33_tot;
  int num_43             = bonded->grp_bond_con.num_43;
  int num_43_tot         = bonded->grp_bond_con.num_43_tot;
  int num_46             = bonded->grp_bond_con.num_46;
  int num_46_tot         = bonded->grp_bond_con.num_46_tot;

  double *pvten_inc_glob = ptens->pvten_inc_glob;
  double *pvten_inc      = ptens->pvten_inc;
  double *pvten_tmp      = ptens->pvten_tmp;
  double *pvten_inc_old  = ptens->pvten_inc_old;
  double *pvten_inc_whol = ptens->pvten_inc_whol;
  double *pvten_inc_std  = ptens->pvten_inc_std;

  double *vgmat          = par_rahman->vgmat;
  double *hmat           = cell->hmat;

  int iperd              = cell->iperd;
  int hmat_cons_typ      = cell->hmat_cons_typ;
  int hmat_int_typ       = cell->hmat_int_typ;

/*=======================================================================*/
/* I) Save old pressure tensor                                          */

  for(i=1;i<=9;i++){pvten_inc_glob[i] = pvten_inc[i];}
  for(i=1;i<=9;i++){pvten_inc_old[i]  = pvten_inc[i];}

  constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_inc_old);

  for(i=1;i<=9;i++){pvten_inc[i]=0;}

/*=======================================================================*/
/*=======================================================================*/
/* Big while convergence loop                                            */

  tol_whol = tolshake+1.0;
  iter_std = 0;
  iter_whol = 0;
  if(ifirst==2){init=1;ifirst=1;}
  ifirst2 =ifirst+1;

  aiter_21_now = 0.0;
  aiter_23_now = 0.0;
  aiter_33_now = 0.0;
  aiter_46_now = 0.0;
  aiter_43_now = 0.0;

  while((tol_whol>tolshake)&&
        (iter_std<=max_iter)&&
        (iter_whol<=max_iter)){
    iter_whol++;
    for(i=1;i<=9;i++){pvten_inc_whol[i]=pvten_inc[i];}

/*=======================================================================*/
/* I) Group constraints                                                  */

    switch(iroll){
     case 0:
      if(num_21 > 0){
        shake_21(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,&aiter_21_now,class_comm_forc_pkg);}
      if(num_23 > 0){
        shake_23(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,&aiter_23_now,class_comm_forc_pkg);}
      if(num_33 > 0){
        shake_33(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,&aiter_33_now,class_comm_forc_pkg);}
      if(num_43 > 0){
       shake_43(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_43_now,class_comm_forc_pkg);}
      if(num_46 > 0){
       shake_46(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_46_now,class_comm_forc_pkg);}
      break;
     case 1:
      if(num_21_tot > 0 ){
        shake_21_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_21_now,baro,ifirst2,class_comm_forc_pkg);}
      if(num_23_tot > 0){
        shake_23_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_23_now,baro,ifirst2,class_comm_forc_pkg);}
      if(num_33_tot > 0){
        shake_33_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_33_now,baro,ifirst2,class_comm_forc_pkg);}
      if(num_46_tot > 0){
        shake_46_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_46_now,baro,ifirst2,class_comm_forc_pkg);}
      if(num_43_tot > 0){
        shake_43_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_43_now,baro,ifirst2,class_comm_forc_pkg);}
      break;
     case 2:
      if(num_21_tot > 0){
        shake_21_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_21_now,par_rahman,ifirst2,cell,
                  class_comm_forc_pkg);}
      if(num_23_tot > 0){
        shake_23_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_23_now,par_rahman,ifirst2,cell,
                  class_comm_forc_pkg);}
      if(num_33_tot > 0){
        shake_33_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_33_now,par_rahman,ifirst2,cell,
                  class_comm_forc_pkg);}
      if(num_46_tot > 0){
        shake_46_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_46_now,par_rahman,ifirst2,cell,
                  class_comm_forc_pkg);}
      if(num_43_tot > 0){
        shake_43_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,&aiter_43_now,par_rahman,ifirst2,cell,
                  class_comm_forc_pkg);}
      break;
    }/*endswitch*/
    stat_avg->iter_21 += aiter_21_now;
    stat_avg->iter_23 += aiter_23_now;
    stat_avg->iter_33 += aiter_33_now;
    stat_avg->iter_46 += aiter_46_now;
    stat_avg->iter_43 += aiter_43_now;

/*=======================================================================*/
/* II) Standard constraints                                              */

    if(iter_whol==1){ifirst2 = 1;}
    tol_std = tolshake+1.0;
    ncon_use = bond_ncon;
    if(iroll > 0){ncon_use = bond_ncon_tot;}   
    while((tol_std>tolshake)&&
          (iter_std<=max_iter)&&
          (ncon_use>0)){
      iter_std++;
      for(i=1;i<=9;i++){pvten_inc_std[i]=pvten_inc[i];}

      switch(iroll){
       case 0:
         shake_bond(&(bonded->bond),clatoms_info,clatoms_pos,
                    cell,&(bonded->intra_scr),ptens,dt,iter_std); 
         break;
       case 1:
         shake_bond_roll_i(&(bonded->bond),clatoms_info,clatoms_pos,
                    cell,&(bonded->intra_scr),ptens,baro,dt,iter_std,
                    class_comm_forc_pkg); 
         break;
       case 2:
         shake_bond_roll_f(&(bonded->bond),clatoms_info,clatoms_pos,
                    cell,&(bonded->intra_scr),ptens,par_rahman,dt,iter_std,
                    class_comm_forc_pkg);
         break;
      }/*endswitch*/
      if(iter_std==1){tol_std = tolshake+1;}
      if(iter_std>1){
        tol_std = 0.0;
        pinc_max = 1.0;
        for(i=1;i<=9;i++){
          tol_std +=beta*fabs(pvten_inc_std[i]-pvten_inc[i]);
          pinc_now = beta*fabs(pvten_inc[i]);
          pinc_max = MAX(pinc_max,pinc_now);
        }/*endfor*/
        if(init==1){tol_std /= pinc_max;}
      }/*endif*/
      ifirst2 = 0;
    }/*endwhile: standard shake not converged*/

/*=======================================================================*/
/* III) Recalc and check tolerence                                       */

    if(bond_ncon==0){
     switch(iroll){
      case 1: recalc_pos_rolli(clatoms_info,clatoms_pos,
                   cell,ptens,baro,dt,ifirst2,class_comm_forc_pkg); 
              break;
      case 2: recalc_pos_rollf(clatoms_info,clatoms_pos,
                   cell,ptens,par_rahman,dt,ifirst2,class_comm_forc_pkg);
              break;
     }/*switch*/
    }/*endif*/
    ifirst2 = 0;
    tol_whol = 0.0;
    for(i=1;i<=9;i++){
       tol_whol+=beta*fabs(pvten_inc_whol[i]-pvten_inc[i]);
    }/*endfor*/
    if(iroll==0){break;}

  }/*endwhile: standard+group not converged */
  stat_avg->iter_shake += iter_std;

/*=======================================================================*/
/*=======================================================================*/
/*  III)Error message                                                    */

  if((iter_std>=max_iter) || 
     (iter_whol>=max_iter)){
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Shake not converged after %d iterations.\n",
            max_iter);
    printf("The present tolerance is %g %g\n",tol_std,tol_whol);
    printf("The desired tolerance is %g \n",tolshake);
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

/*=======================================================================*/
/*  IV)Calculate the global tolernce for main program                    */

  *tol_glob = 0.0;
  for(i=1;i<=9;i++){
    *tol_glob+=beta*fabs(pvten_inc_glob[i]-pvten_inc[i]);
  }/*endfor*/

/*=======================================================================*/
/*  IV)Allreduce pvten_inc     */

  if(np_forc > 1 && iroll==0){
   Allreduce(&(pvten_inc[1]),&(pvten_tmp[1]),9,MPI_DOUBLE,MPI_SUM,0,comm_forc);
   for(i=1;i<=9;i++){
    pvten_inc[i] = pvten_tmp[i];
   }/*endfor*/
  }/*endif*/
 
/*=======================================================================*/
/*  IV)DEBUG  */

#ifdef DEBUG
  for(i=1;i<=9;i+=3){
     printf("SHAKE Pkj %.10g %.10g %.10g\n",
              pvten_inc[i],pvten_inc[i+1],pvten_inc[i+2]);
  }/*endfor*/
  printf("iter_whol=%d,iter_std=%d\n",iter_whol,iter_std);
  scanf("%d",&iii);

      printf("x(1),y(1),z(1) %.13g %.13g %.13g\n",
                                         clatoms_pos->x[1],
                                         clatoms_pos->y[1],
                                         clatoms_pos->z[1]); 
      printf("vx(1),vy(1),vz(1) %.13g %.13g %.13g\n",
                                         clatoms_pos->vx[1],
                                         clatoms_pos->vy[1],
                                         clatoms_pos->vz[1]); 
      printf("fx(1),fy(1),fz(1) %.13g %.13g %.13g\n",
                                         clatoms_pos->fx[1],
                                         clatoms_pos->fy[1],
                                         clatoms_pos->fz[1]); 
      scanf("%d",&iii);
#endif

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rattle_control(BONDED *bonded,
                    CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos, 
                    CELL *cell,PTENS *ptens,
                    STATEPOINT *statepoint,BARO *baro,PAR_RAHMAN *par_rahman,
                    STAT_AVG *stat_avg,double dt,double *tol_glob,int ifirst,
                    CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
                    EWD_SCR *ewd_scr)

/*==========================================================================*/
/* WARNING DANGER dt may not be timeinfo.dt so let it come in    */
/* as a parameter                                                */
/*==========================================================================*/
/*           Begin Routine  */
    {/*Begin Routine*/
/*=======================================================================*/
/*           Local Variable declarations                                   */
  
#include "../typ_defs/typ_mask.h"

  int i,iii,ifirst2;
  int iter_std,iter_whol;
  int ncon_use;
  double tol_whol,tol_std;

  double beta            = BOLTZ/(statepoint->t_ext); 
  int np_forc            = class_comm_forc_pkg->num_proc;
  MPI_Comm comm_forc     = class_comm_forc_pkg->comm;

  double tolratl         = bonded->constrnt.tolratl;
  int max_iter           = bonded->constrnt.max_iter;
  int iroll              = bonded->constrnt.iroll;

  int bond_ncon          = bonded->bond.ncon;
  int bond_ncon_tot      = bonded->bond.ncon_tot;
  int num_21             = bonded->grp_bond_con.num_21;
  int num_21_tot         = bonded->grp_bond_con.num_21_tot;
  int num_23             = bonded->grp_bond_con.num_23;
  int num_23_tot         = bonded->grp_bond_con.num_23_tot;
  int num_33             = bonded->grp_bond_con.num_33;
  int num_33_tot         = bonded->grp_bond_con.num_33_tot;
  int num_43             = bonded->grp_bond_con.num_43;
  int num_43_tot         = bonded->grp_bond_con.num_43_tot;
  int num_46             = bonded->grp_bond_con.num_46;
  int num_46_tot         = bonded->grp_bond_con.num_46_tot;

  double *pvten_inc_glob = ptens->pvten_inc_glob;
  double *pvten_inc      = ptens->pvten_inc;
  double *pvten_tmp      = ptens->pvten_tmp;
  double *pvten_inc_old  = ptens->pvten_inc_old;
  double *pvten_inc_whol = ptens->pvten_inc_whol;
  double *pvten_inc_std  = ptens->pvten_inc_std;

  int iperd              = cell->iperd;
  int hmat_cons_typ      = cell->hmat_cons_typ;
  int hmat_int_typ       = cell->hmat_int_typ;

/*=======================================================================*/
/* I) Save old pressure tensor                                          */

  for(i=1;i<=9;i++){pvten_inc_glob[i] = pvten_inc[i];}
  for(i=1;i<=9;i++){pvten_inc_old[i]  = pvten_inc[i];}

  constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_inc_old);

  for(i=1;i<=9;i++){pvten_inc[i]=0;}

/*=======================================================================*/
/*=======================================================================*/
/* Big while convergence loop                                            */

  tol_whol = tolratl+1.0;
  iter_std = 0;
  iter_whol = 0;
  ifirst2 =ifirst + 1;
  while((tol_whol>tolratl)&&
        (iter_std<=max_iter)&&
        (iter_whol<=max_iter)){
    iter_whol++;
    for(i=1;i<=9;i++){pvten_inc_whol[i]=pvten_inc[i];}

/*=======================================================================*/
/* I) Group constraints                                                  */

    switch(iroll){
     case 0:
      if(num_23 > 0){
        rattle_23(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,class_comm_forc_pkg);}
      if(num_21 > 0){
        rattle_21(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,class_comm_forc_pkg);}
      if(num_33 > 0){
       rattle_33(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                 ptens,dt,class_comm_forc_pkg);}
      if(num_46 > 0){
        rattle_46(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,class_comm_forc_pkg);}
      if(num_43 > 0){
        rattle_43(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                  ptens,dt,class_comm_forc_pkg);}
      break;
     case 1:
      if(num_21_tot > 0){
       rattle_21_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,baro,ifirst2,class_comm_forc_pkg);}
      if(num_23_tot > 0){
       rattle_23_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,baro,ifirst2,class_comm_forc_pkg);}
      if(num_33_tot > 0){
       rattle_33_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,baro,ifirst2,class_comm_forc_pkg);}
      if(num_46_tot > 0){
       rattle_46_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,baro,ifirst2,class_comm_forc_pkg);}
      if(num_43_tot > 0){
       rattle_43_rolli(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                   ptens,dt,baro,ifirst2,class_comm_forc_pkg);}
      break;
     case 2:
      if(num_21_tot > 0){
        rattle_21_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                    ptens,dt,par_rahman,ifirst2,cell,class_comm_forc_pkg);}
      if(num_23_tot > 0){
        rattle_23_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                    ptens,dt,par_rahman,ifirst2,cell,class_comm_forc_pkg);}
      if(num_33_tot > 0){
        rattle_33_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                    ptens,dt,par_rahman,ifirst2,cell,class_comm_forc_pkg);}
      if(num_46_tot > 0){
        rattle_46_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                    ptens,dt,par_rahman,ifirst2,cell,class_comm_forc_pkg);}
      if(num_43_tot > 0){
        rattle_43_rollf(&(bonded->grp_bond_con),clatoms_info,clatoms_pos,
                    ptens,dt,par_rahman,ifirst2,cell,class_comm_forc_pkg);}
      break;
    }/*endswitch*/

/*=======================================================================*/
/* II) Standard constraints                                              */

    if(iter_whol==1){ifirst2 = 1;}
    tol_std = tolratl+1.0;
    ncon_use = bond_ncon;
    if(iroll > 0){ncon_use = bond_ncon_tot;}   
    while((tol_std>tolratl)&&
          (iter_std<=max_iter)&&
          (ncon_use>0)){
      iter_std++;
      for(i=1;i<=9;i++){pvten_inc_std[i]=pvten_inc[i];}
      switch(iroll){
       case 0:
         rattle_bond(&(bonded->bond),clatoms_info,clatoms_pos,
                     cell,&(bonded->intra_scr),ptens,dt,iter_std);
         break;
       case 1:
         rattle_bond_roll_i(&(bonded->bond),clatoms_info,clatoms_pos,
                     cell,&(bonded->intra_scr),ptens,baro,dt,iter_std,
                     class_comm_forc_pkg);
         break;
       case 2:
         rattle_bond_roll_f(&(bonded->bond),clatoms_info,clatoms_pos,
                      cell,&(bonded->intra_scr),ptens,par_rahman,dt,iter_std,
                      class_comm_forc_pkg);
         break;
      }/*endswitch*/
      if(iter_std==1){tol_std = tolratl+1;}
      if(iter_std>1){
        tol_std = 0.0;
        for(i=1;i<=9;i++){
          tol_std+=beta*fabs(pvten_inc_std[i]-pvten_inc[i]);
        }/*endfor*/
      }/*endif*/
      ifirst2 = 0;
    }/*endwhile: standard ratl not converged*/

/*=======================================================================*/
/* III) Recalc and check tolerence                                       */

    if(ifirst2==1){
     switch(iroll){
      case 1: recalc_vel_rolli(cell,ptens,baro,dt,class_comm_forc_pkg); 
              break;
      case 2: recalc_vel_rollf(cell,ptens,par_rahman,dt,class_comm_forc_pkg);
              break;
     }/*switch*/
    }/*endif*/
    ifirst2 = 0;
    tol_whol = 0.0;
    for(i=1;i<=9;i++){
       tol_whol+=beta*fabs(pvten_inc_whol[i]-pvten_inc[i]);
    }/*endfor*/
    if(iroll==0){break;}

  }/*endwhile: standard+group not converged */

  if(num_23 > 0){stat_avg->iter_23r += iter_whol;}
  if(num_21 > 0){stat_avg->iter_21r += iter_whol;}
  if(num_33 > 0){stat_avg->iter_33r += iter_whol;}
  if(num_46 > 0){stat_avg->iter_46r += iter_whol;}
  if(num_43 > 0){stat_avg->iter_43r += iter_whol;}
  stat_avg->iter_ratl += iter_std;

/*=======================================================================*/
/*=======================================================================*/
/*  III)Error message                                                    */

  if((iter_std>=max_iter) || (iter_whol>=max_iter)){
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Rattle not converged after %d iterations.\n",max_iter);
    printf("The present tolerance is %g %g\n",tol_std,tol_whol);
    printf("The desired tolerance is %g \n",tolratl);
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

/*=======================================================================*/
/*  IV)Calculate the global tolernce for main program                    */

  *tol_glob = 0.0;
  for(i=1;i<=9;i++){
    *tol_glob+=beta*fabs(pvten_inc_glob[i]-pvten_inc[i]);
  }/*endfor*/

/*=======================================================================*/
/*  IV)Allreduce pvten_inc     */

  if(np_forc > 1 && iroll==0){
   Allreduce(&(pvten_inc[1]),&(pvten_tmp[1]),9,MPI_DOUBLE,MPI_SUM,0,comm_forc);
   for(i=1;i<=9;i++){
     pvten_inc[i] = pvten_tmp[i];
   }/*endfor*/
  }/*endif*/

#ifdef DEBUG
  for(i=1;i<=9;i+=3){
     printf("RATTLE Pkj %.10g %.10g %.10g\n",
              pvten_inc[i],
              pvten_inc[i+1],
              pvten_inc[i+2]);
  }/*endfor*/
  printf("iter_whol=%d,iter_std=%d\n",iter_whol,iter_std);
  scanf("%d",&iii);
#endif

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_constraint(BONDED *bonded,PTENS *ptens)

/*==========================================================================*/
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  
  int i;

  int nbond           = bonded->bond.ncon;
  int nbend           = bonded->bend.ncon;
  int ntors           = bonded->tors.ncon;

  double *bond_al_con = bonded->bond.al_con;
  double *bend_al_con = bonded->bend.al_con;
  double *tors_al_con = bonded->tors.al_con;
 
  double *pvten_inc   = ptens->pvten_inc;

/*=======================================================================*/
/* I) Initialize the multipliers */

  for(i=1;i<=nbond;i++){bond_al_con[i]=0.0;}
  for(i=1;i<=nbend;i++){bend_al_con[i]=0.0;}
  for(i=1;i<=ntors;i++){tors_al_con[i]=0.0;}

/*=======================================================================*/
/* II) Initialize the pressure tensor */

  for(i=1;i<=9;i++){pvten_inc[i]=0;}

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void recalc_pos_rolli(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CELL *cell,PTENS *ptens,
                      BARO *baro,double dt,int ifirst,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*           Begin Routine  */
{/*Begin Routine*/
/*=======================================================================*/
/*           Local Variable declarations                                   */

  int ipart,i;
  double ftemp;
  double e2,e4,e6,e8;
  double aa,aa2,arg2,poly,bb,dlen;

  double *clatoms_x       = clatoms_pos->x;
  double *clatoms_y       = clatoms_pos->y;
  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_vx      = clatoms_pos->vx;
  double *clatoms_vy      = clatoms_pos->vy;
  double *clatoms_vz      = clatoms_pos->vz;
  double *clatoms_xold    = clatoms_info->xold;
  double *clatoms_yold    = clatoms_info->yold;
  double *clatoms_zold    = clatoms_info->zold;
  int myatm_start         = clatoms_info->myatm_start;
  int myatm_end           = clatoms_info->myatm_end;

  double *x_lnv           = &(baro->x_lnv);
  double *x_lnv_o         = &(baro->x_lnv_o);
  double *v_lnv           = &(baro->v_lnv);
  double *f_lnv_p         = &(baro->f_lnv_p);
  double mass_lnv         = baro->mass_lnv;
  double roll_scg         = baro->roll_scg;
  double *hmato           = baro->hmato;

  double *pvten_inc       = ptens->pvten_inc;
  double *pvten_inc_old   = ptens->pvten_inc_old;
  double *hmat            = cell->hmat;
  double *hmati           = cell->hmati;
  int    iperd            = cell->iperd;

  double dt2 = dt*0.5;
  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);e6=e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*==========================================================================*/
/* I) Evolve the volume velocity            */

  if(ifirst != 0){
    ftemp   = (pvten_inc[1]+pvten_inc[5]+pvten_inc[9])
            - (pvten_inc_old[1]+pvten_inc_old[5]+pvten_inc_old[9]);
    (*f_lnv_p) += ftemp;
    (*v_lnv)   += 0.5*ftemp*(roll_scg)*dt/(mass_lnv);
  }/*endif*/

/*==========================================================================*/
/* II) Evolve atom positions */

  aa   = exp( dt2*((*v_lnv)) );
  aa2  = aa*aa;
  arg2 = (((*v_lnv))*dt2)*(((*v_lnv))*dt2);
  poly = (((e8*arg2+e6)*arg2+e4)*arg2+e2)*arg2+1.0;
  bb   = aa*poly;
  for(ipart=myatm_start;ipart<=(myatm_end);ipart++){
    clatoms_x[ipart] = aa2*clatoms_xold[ipart] + bb*clatoms_vx[ipart]*dt;
    clatoms_y[ipart] = aa2*clatoms_yold[ipart] + bb*clatoms_vy[ipart]*dt;
    clatoms_z[ipart] = aa2*clatoms_zold[ipart] + bb*clatoms_vz[ipart]*dt;
  }/*endfor*/
  baro->roll_scv = bb;

/*==========================================================================*/
/* III) Evolve volume positions */

  (*x_lnv) = (*x_lnv_o) + ((*v_lnv))*dt;
  dlen     = exp( (*x_lnv)-(*x_lnv_o) );

  for(i=1;i<=9;i++){hmat[i] = hmato[i]*dlen;}
  gethinv(hmat,hmati,&(cell->vol),iperd);
  baro->vol = cell->vol;

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void recalc_pos_rollf(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CELL *cell,PTENS *ptens,
                      PAR_RAHMAN *par_rahman,double dt,int ifirst,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*           Begin Routine  */
   {/*Begin Routine*/
/*=======================================================================*/
/*           Local Variable declarations                                   */

  int ipart,i,j,n,joff;
  double ftemp;
  double e2,e4,e6,e8;
  double aa,arg2,poly;
  double tempx,tempy,tempz;
  double tempvx,tempvy,tempvz;

  double *clatoms_x     = clatoms_pos->x;
  double *clatoms_y     = clatoms_pos->y;
  double *clatoms_z     = clatoms_pos->z;
  double *clatoms_vx    = clatoms_pos->vx;
  double *clatoms_vy    = clatoms_pos->vy;
  double *clatoms_vz    = clatoms_pos->vz;
  double *clatoms_xold  = clatoms_info->xold;
  double *clatoms_yold  = clatoms_info->yold;
  double *clatoms_zold  = clatoms_info->zold;
  int myatm_start       = clatoms_info->myatm_start;
  int myatm_end         = clatoms_info->myatm_end;

  double *roll_mtv      = par_rahman->roll_mtv;
  double *roll_mtx      = par_rahman->roll_mtx;
  double roll_scg       = par_rahman->roll_scg;
  double mass_hm        = par_rahman->mass_hm;
  double *fgmat_p       = par_rahman->fgmat_p;
  double *vgmat         = par_rahman->vgmat;
  double *vtemps        = par_rahman->vtemps;
  double *vtempv        = par_rahman->vtempv;
  double *veig          = par_rahman->veig;
  double *veigv         = par_rahman->veigv;
  double *vexpdt        = par_rahman->vexpdt;
  double *vsindt        = par_rahman->vsindt;
  double *vtempx        = par_rahman->vtempx;
  double *fv1           = par_rahman->fv1;
  double *fv2           = par_rahman->fv2;
  double *hmato         = par_rahman->hmato;
  double *hmat_t        = par_rahman->hmat_t;

  double *hmat          = cell->hmat;
  double *hmati         = cell->hmati;
  int iperd             = cell->iperd;
  int hmat_cons_typ     = cell->hmat_cons_typ;
  int hmat_int_typ      = cell->hmat_int_typ;

  double *pvten_tmp     = ptens->pvten_tmp;
  double *pvten_inc     = ptens->pvten_inc;
  double *pvten_inc_old = ptens->pvten_inc_old;

  double dt2 = dt*0.5;
  e2=1.0/(2.0*3.0);e4=e2/(4.0*5.0);e6=e4/(6.0*7.0);e8=e6/(8.0*9.0);

/*=======================================================================*/
/* I) Update box variables */

  if(ifirst==1){
    for(i=1;i<=9;i++){pvten_tmp[i]=pvten_inc[i];}
    constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);    
    for(i=1;i<=9;i++){      
       ftemp       = pvten_tmp[i]-pvten_inc_old[i];
       fgmat_p[i] += ftemp;
       vgmat[i]   += (0.50*ftemp*roll_scg*dt/mass_hm);
    }/*endfor*/

  }/*endfor*/

/*=======================================================================*/
/* II) Update the particles  */

  for(i=1;i<=9;i++){vtemps[i]=vgmat[i];}
  diag33(vtemps,veig,veigv,fv1,fv2);
  for(i=1;i<=3;i++){
    aa         = exp(dt2*veig[i]);
    vexpdt[i]  = aa*aa;
    arg2       = (veig[i]*dt2)*(veig[i]*dt2);
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
  for(ipart=myatm_start;ipart<=myatm_end; ++ipart) {
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

/*----------------------------------------------------------------------*/
/* B) get the new matrix of cell parameters and their inverse */

   n = 3;
   matmul_t(hmato,veigv,hmat_t,n);
   for(i=1;i<=3;i++){
     joff = (i-1)*3 ;
     for(j=1;j<=3;j++){ 
       hmat_t[j+joff] *=  vexpdt[j];
     }/*endfor*/
   }/*endfor*/
   matmul_tt(hmat_t,veigv,cell->hmat,n);
   gethinv(hmat,hmati,&(cell->vol),iperd);
   par_rahman->vol = cell->vol;

/*--------------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void recalc_vel_rolli(CELL *cell,PTENS *ptens,BARO *baro,double dt,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)


/*==========================================================================*/
/*           Begin Routine  */
  {/*Begin Routine*/
/*=======================================================================*/
/*           Local Variable declarations                                   */

  double ftemp;
  double dt2 = dt*0.5;

  double *pvten_inc     = ptens->pvten_inc;
  double *pvten_inc_old = ptens->pvten_inc_old;
  double *f_lnv_p       = &(baro->f_lnv_p);
  double *v_lnv_g       = &(baro->v_lnv_g);
  double roll_scg       = baro->roll_scg;
  double mass_lnv       = baro->mass_lnv;

/*==========================================================================*/
/* I) Evolve the volume velocity            */

  ftemp      = (pvten_inc[1]+pvten_inc[5]+pvten_inc[9])
             - (pvten_inc_old[1]+pvten_inc_old[5]+pvten_inc_old[9]);
  (*f_lnv_p) += ftemp;
  (*v_lnv_g) += (ftemp*roll_scg*dt2/mass_lnv);

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void recalc_vel_rollf(CELL *cell,PTENS *ptens,PAR_RAHMAN *par_rahman,double dt,
                      CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
/*           Begin Routine  */
{/*Begin Routine*/
/*=======================================================================*/
/*           Local Variable declarations                                   */

  int i;
  double dt2 = dt*0.5;

  double *pvten_inc     = ptens->pvten_inc;
  double *pvten_inc_old = ptens->pvten_inc_old;
  double *pvten_tmp     = ptens->pvten_tmp;

  double *fgmat_p       = par_rahman->fgmat_p;
  double *vgmat_g       = par_rahman->vgmat_g;
  double mass_hm        = par_rahman->mass_hm;
  double roll_scg       = par_rahman->roll_scg;

  int iperd             = cell->iperd;
  int hmat_cons_typ     = cell->hmat_cons_typ;
  int hmat_int_typ      = cell->hmat_int_typ;

/*=======================================================================*/

  for(i=1;i<=9;i++){pvten_tmp[i]=pvten_inc[i];}
  constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,pvten_tmp);    
  for(i=1;i<=9;i++){pvten_tmp[i]-=pvten_inc_old[i];}

  for(i=1;i<=9;i++){      
    fgmat_p[i] += pvten_tmp[i];
    vgmat_g[i] += (pvten_tmp[i]*roll_scg*dt2/mass_hm);
  }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void zero_constrt_iters(STAT_AVG *stat_avg)

/*==========================================================================*/
    {/*begin routine*/
/*========================================================================*/

    stat_avg->iter_shake    = 0; 
    stat_avg->iter_ratl     = 0; 
    stat_avg->iter_21r      = 0.0;
    stat_avg->iter_23r      = 0.0;
    stat_avg->iter_33r      = 0.0;
    stat_avg->iter_43r      = 0.0;
    stat_avg->iter_46r      = 0.0;
    stat_avg->iter_21       = 0.0;
    stat_avg->iter_23       = 0.0;
    stat_avg->iter_33       = 0.0;
    stat_avg->iter_43       = 0.0;
    stat_avg->iter_46       = 0.0;
    stat_avg->iter_shake_cp = 0; 
    stat_avg->iter_ratl_cp  = 0; 

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/
