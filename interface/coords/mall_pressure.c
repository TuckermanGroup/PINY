/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void mall_pressure(CLASS *class,GENERAL_DATA *general_data)

/*======================================================================*/
     {/*begin routine*/
/*======================================================================*/
/*          Local variable declarations                                */

  double now_memory;
  int myid = class->communicate.myid;

/*======================================================================*/
/*         Output                                                       */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Setting up cell/pressure related quantities\n");
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/ 

  now_memory   = (sizeof(int)*0 + sizeof(double)*29*9)*1.e-06;
  class->tot_memory  += now_memory;

  if(myid==0){
    printf("Pressure allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,class->tot_memory);
  }/*endif*/

  /*=======================================================================*/
  /* 0) Average pressure Memory */

  general_data->stat_avg.apten_out = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->stat_avg.apten     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->stat_avg.aipten    = (double *)cmalloc((size_t)9*sizeof(double))-1;

  general_data->ptens.tvten        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tot    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tmp    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_tmp_res = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t1 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t2 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t3 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_t4 = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.count_inc    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_a  = (double *)cmalloc((size_t)9*sizeof(double))-1;    
  general_data->ptens.pvten_inc_old = (double *)cmalloc((size_t)9*sizeof(double))-1;    
  general_data->ptens.pvten_inc_glob = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_whol = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_inc_std  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_mode     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->ptens.pvten_mode_tot = (double *)cmalloc((size_t)9*sizeof(double))-1;

  /*=======================================================================*/
  /* 0) Baro/cell/par_rahman Memory */

  general_data->par_rahman.vgmat    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vgmat_glob= (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vgmat_g  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fgmat_p  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fgmat_v  = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtv = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtx = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.roll_mtvv = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtemps   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.veigv    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.veig     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vexpdt   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vsindt   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtempx   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.vtempv   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.hmat_t   = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.hmato    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fv1    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->par_rahman.fv2    = (double *)cmalloc((size_t)9*sizeof(double))-1;

  general_data->baro.hmato       = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat        = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmati       = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_ewd    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_ewd_cp = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmat_cp     = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.hmati_cp    = (double *)cmalloc((size_t)9*sizeof(double))-1;
  general_data->cell.cp_box_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_box_center_rel = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_vbox_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;
  general_data->cell.cp_fbox_center = (double *) cmalloc((size_t)3*sizeof(double ))-1;

  /*=======================================================================*/
  /* 0) More output */

  if(myid==0){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed cell/pressure set up\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/







