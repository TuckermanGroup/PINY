/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: simpavg_md_communicate                       */
/*                                                                          */
/* This subprogram provides communication for the                           */ 
/* instantaneous averages of PI_MD                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*               Header:                                                    */

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/

void simpavg_pimd_communicate(STAT_AVG *stat_avg,ENSOPTS *ensopts,
                    PTENS *ptens, int iconstrnt,COMMUNICATE *communicate)

/*==========================================================================*/
    {/*begin routine*/

#include "../typ_defs/typ_mask.h"

    int i,iii;
    double tvten_temp = 0.0;
    double iter_shake_temp=0.0;
    double iter_ratl_temp=0.0;
    double iter_23_temp=0.0;
    double iter_33_temp=0.0;
    double iter_46_temp=0.0;
    double iter_23r_temp=0.0;
    double iter_33r_temp=0.0;
    double iter_46r_temp=0.0;
    double kinet_temp=0.0;
    double vintert_temp=0.0;
    double vintrat_temp=0.0;
    double kinet_nhc_temp=0.0;
    double vpotnhc_temp=0.0;
    double vlong_temp=0.0;
    double kinet_v_temp=0.0;
    double vpot_v_temp=0.0;
    double kinet_nhc_bead_temp=0.0;
    double kin_harm_temp=0.0;
    double pi_ke_prim=0.0;
    double pi_ke_vir=0.0;
    double press_inter_temp=0.0;
    double press_intra_temp=0.0;
    int np_beads = communicate->np_beads;
    int np_forc  = communicate->np_forc;
    MPI_Comm comm_beads = communicate->comm_beads;
    MPI_Comm comm_forc = communicate->comm_forc;
    MPI_Comm world      = communicate->world;


  if(iconstrnt==1){
   Reduce(&(stat_avg->iter_shake), &iter_shake_temp,1,MPI_INT,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_ratl), &iter_ratl_temp,1,MPI_INT,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_23), &iter_23_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_33), &iter_33_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_46), &iter_46_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_23r), &iter_23r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_33r), &iter_33r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
   Reduce(&(stat_avg->iter_46r), &iter_46r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,world);
    stat_avg->iter_shake =  iter_shake_temp;
    stat_avg->iter_ratl  =  iter_ratl_temp;
    stat_avg->iter_23    =  iter_23_temp;
    stat_avg->iter_33    =  iter_33_temp;
    stat_avg->iter_46    =  iter_46_temp;
    stat_avg->iter_23r   =  iter_23r_temp;
    stat_avg->iter_33r   =  iter_33r_temp;
    stat_avg->iter_46r   =  iter_46r_temp;
  }/*endif*/

  /*=======================================================================*/
  /* III) Construct the Conserved Quantity */

   Reduce(&(stat_avg->vintert), &vintert_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vintrat), &vintrat_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->kin_harm), &kin_harm_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->pi_ke_vir), &pi_ke_vir,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->kinet), &kinet_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->vintert = vintert_temp;
   stat_avg->vintrat = vintrat_temp;
   stat_avg->kin_harm = kin_harm_temp;
   stat_avg->pi_ke_vir = pi_ke_vir;
   stat_avg->kinet = kinet_temp;

   Reduce(&(stat_avg->kinet_nhc_bead), &kinet_nhc_bead_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->kinet_nhc_bead = kinet_nhc_bead_temp;
  if(np_beads>1){   
   Reduce(&(stat_avg->pi_ke_prim), &pi_ke_prim,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   stat_avg->pi_ke_prim = pi_ke_prim;
  }

  if((ensopts->nvt)==1){   
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
  }/*endif for nvt*/

  if((ensopts->npt_i)==1){
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(stat_avg->vlong), &vlong_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->vlong = vlong_temp;

   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
 }/*endif npt_i*/

  if((ensopts->npt_f)==1) {
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vlong), &vlong_temp,1,MPI_DOUBLE,
                  MPI_SUM,0,world);
   stat_avg->vlong = vlong_temp;

   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
  }/*endif for npt_f*/

  /*=======================================================================*/
  /* IV) Construct the Pressure Tensor Quantities                          */

   Reduce(&(stat_avg->press_inter),&press_inter_temp,1,MPI_DOUBLE,
                 MPI_SUM,0,world);
   Reduce(&(stat_avg->press_intra),&press_intra_temp,1,MPI_DOUBLE,
                 MPI_SUM,0,world);
   stat_avg->press_inter = press_inter_temp;
   stat_avg->press_intra = press_intra_temp;
   for ( i = 1;i<=9;i++){
      Reduce(&(ptens->tvten[i]), &tvten_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
      ptens->tvten[i]     = tvten_temp;
   }
  }/*end routine*/











