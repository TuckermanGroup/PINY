/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*             Module: communicate_simpavg_cp_pimd                          */
/*                                                                          */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_output_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_simpavg_cp_pimd(STAT_AVG *stat_avg, PTENS *ptens, 
                                COMMUNICATE *communicate, SIMOPTS *simopts,
                                ENSOPTS *ensopts,int iconstrnt,int num_c_nhc)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
    int i,iii;
    double pvten_temp[10];
    double tvten_temp = 0.0;
    double cp_ehart=0.0;
    double cp_enl  =0.0;
    double cp_eke  =0.0;
    double cp_exc  =0.0;
    double cp_muxc  =0.0;
    double cp_eext  =0.0;
    double kinet_cp=0.0;
    double kinet_nhc_cp=0.0;
    double vpotnhc_cp=0.0;
    double iter_shake_temp=0.0;
    double iter_ratl_temp=0.0;
    double iter_23_temp=0.0;
    double iter_33_temp=0.0;
    double iter_46_temp=0.0;
    double iter_23r_temp=0.0;
    double iter_33r_temp=0.0;
    double iter_46r_temp=0.0;
    double kinet_temp=0.0;
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
    int np_states       = communicate->np_states;
    int myid_state      = communicate->myid_state;
    MPI_Comm world      = communicate->world;
    MPI_Comm comm_beads = communicate->comm_beads;
    MPI_Comm comm_states = communicate->comm_states;

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

 if(myid_state==0){
   Reduce(&(stat_avg->vintrat), &vintrat_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   Reduce(&(stat_avg->kin_harm), &kin_harm_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   Reduce(&(stat_avg->pi_ke_prim), &pi_ke_prim,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   Reduce(&(stat_avg->pi_ke_vir), &pi_ke_vir,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   Reduce(&(stat_avg->kinet), &kinet_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_beads);
   stat_avg->vintrat = vintrat_temp;
   stat_avg->kin_harm = kin_harm_temp;
   stat_avg->pi_ke_prim = pi_ke_prim;
   stat_avg->pi_ke_vir = pi_ke_vir;
   stat_avg->kinet = kinet_temp;
 }/*endif : myid_state==0*/

if(simopts->cp == 1 || simopts->cp_pimd==1) {
  if((ensopts->nvt)==1){   
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->kinet_nhc_bead), &kinet_nhc_bead_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
   stat_avg->kinet_nhc_bead = kinet_nhc_bead_temp;
  }/*endif for nvt*/

  if((ensopts->npt_i)==1){
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpot_v), &vpot_v_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->kinet_nhc_bead), &kinet_nhc_bead_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);

    Reduce(&(stat_avg->vlong), &vlong_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   stat_avg->vlong = vlong_temp;

   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
   stat_avg->vpot_v = vpot_v_temp;
   stat_avg->kinet_nhc_bead = kinet_nhc_bead_temp;
 }/*endif npt_i*/

  if((ensopts->npt_f)==1) {
   Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->vpot_v), &vpot_v_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
   Reduce(&(stat_avg->kinet_nhc_bead), &kinet_nhc_bead_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    Reduce(&(stat_avg->vlong), &vlong_temp,1,MPI_DOUBLE,
                   MPI_SUM,0,world);
    stat_avg->vlong = vlong_temp;

   stat_avg->kinet_nhc = kinet_nhc_temp;
   stat_avg->vpotnhc = vpotnhc_temp;
   stat_avg->vpot_v = vpot_v_temp;
   stat_avg->kinet_nhc_bead = kinet_nhc_bead_temp;
  }/*endif for npt_f*/
}/*endif : full cp_on*/
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

/*=======================================================================*/
/* III) Conserved Quantity Components                                    */
  
   if(simopts->cp_wave_min_pimd == 1 && np_states > 1){
     Reduce(&(stat_avg->cp_ehart), &cp_ehart,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);
     Reduce(&(stat_avg->cp_eext), &cp_eext,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);
     Reduce(&(stat_avg->cp_exc), &cp_exc,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);
     Reduce(&(stat_avg->cp_muxc), &cp_muxc,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);
     Reduce(&(stat_avg->cp_enl), &cp_enl,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);
     Reduce(&(stat_avg->cp_eke), &cp_eke,1,MPI_DOUBLE,
            MPI_SUM,0,comm_states);


     stat_avg->cp_ehart = cp_ehart;
     stat_avg->cp_eext  = cp_eext;
     stat_avg->cp_exc   = cp_exc;
     stat_avg->cp_muxc  = cp_muxc;
     stat_avg->cp_enl   = cp_enl;
     stat_avg->cp_eke   = cp_eke;

   }


   if(simopts->cp_wave_min_pimd != 1){
      Reduce(&stat_avg->kinet_cp,&kinet_cp,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_ehart,&cp_ehart,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_eext,&cp_eext,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_exc,&cp_exc,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_muxc,&cp_muxc,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_eke,&cp_eke,1,MPI_DOUBLE,MPI_SUM,
                0,world);
      Reduce(&stat_avg->cp_enl,&cp_enl,1,MPI_DOUBLE,MPI_SUM,
                0,world);


      stat_avg->kinet_cp = kinet_cp;
      stat_avg->cp_ehart = cp_ehart;
      stat_avg->cp_eext = cp_eext;
      stat_avg->cp_exc = cp_exc;
      stat_avg->cp_muxc = cp_muxc;
      stat_avg->cp_eke = cp_eke;
      stat_avg->cp_enl = cp_enl;

   }/* endif cp_save_min_pimd */

  if(num_c_nhc > 0) {
    Reduce(&stat_avg->kinet_nhc_cp,&kinet_nhc_cp,1,MPI_DOUBLE,MPI_SUM,
              0,world);
    Reduce(&stat_avg->vpotnhc_cp,&vpotnhc_cp,1,MPI_DOUBLE,MPI_SUM,
              0,world);
    stat_avg->kinet_nhc_cp = kinet_nhc_cp; 
    stat_avg->vpotnhc_cp = vpotnhc_cp; 
  }
   

/*======================================================================*/
}/*end routine*/
/*==========================================================================*/




