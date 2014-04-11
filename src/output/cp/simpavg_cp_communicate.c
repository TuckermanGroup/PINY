#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/

void simpavg_cp_communicate(CP *cp,STAT_AVG *stat_avg,COMMUNICATE *communicate,PTENS *ptens)

/*==========================================================================*/
    {/*begin routine*/
#include "../typ_defs/typ_mask.h"


/*==========================================================================*/
    int i,iii;
    double cp_ehart_tmp=0.0;
    double cp_eext_tmp=0.0;
    double cp_exc_tmp=0.0;
    double cp_muxc_tmp=0.0;
    double cp_enl_tmp=0.0;
    double cp_eke_tmp=0.0;
    double kinet_cp_tmp=0.0;
    double kinet_cp_nhc_tmp=0.0;
    double vpotnhc_cp_tmp=0.0;
    double pvten_temp[10];
    MPI_Comm comm_states = communicate->world;

  /*=======================================================================*/
  /* II)Electron DFT Energies                                              */
   Reduce(&(stat_avg->cp_ehart), &cp_ehart_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->cp_eext), &cp_eext_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->cp_exc), &cp_exc_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->cp_muxc), &cp_muxc_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->cp_enl), &cp_enl_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->cp_eke), &cp_eke_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);

   for(i=1;i<=9;i++){pvten_temp[i]=0.0;}
   Reduce(&(ptens->pvten_tot[1]), &(pvten_temp[1]),9,MPI_DOUBLE,
             MPI_SUM,0,comm_states);

   stat_avg->cp_ehart = cp_ehart_tmp;
   stat_avg->cp_eext  = cp_eext_tmp;
   stat_avg->cp_exc   = cp_exc_tmp;
   stat_avg->cp_muxc  = cp_muxc_tmp;
   stat_avg->cp_enl   = cp_enl_tmp;
   stat_avg->cp_eke   = cp_eke_tmp;

   for(i=1;i<=9;i++){ptens->pvten_tot[i] =  pvten_temp[i];}

  /*=======================================================================*/
  /* III)Electron Nose-Hoover quantities and fictitious electon KE  */

   Reduce(&(stat_avg->kinet_cp), &kinet_cp_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   stat_avg->kinet_cp  = kinet_cp_tmp;

  if(cp->cptherm_info.num_c_nhc > 0){
   Reduce(&(stat_avg->kinet_nhc_cp), &kinet_cp_nhc_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   Reduce(&(stat_avg->vpotnhc_cp), &vpotnhc_cp_tmp,1,MPI_DOUBLE,
                   MPI_SUM,0,comm_states);
   stat_avg->kinet_nhc_cp  = kinet_cp_nhc_tmp;
   stat_avg->vpotnhc_cp  = vpotnhc_cp_tmp;
  }/* endif */


  /*=======================================================================*/
  }/*end routine*/
  /*=======================================================================*/







