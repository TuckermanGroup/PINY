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
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void simpavg_md_communicate(STAT_AVG *stat_avg,ENSOPTS *ensopts,
                    PTENS *ptens, int iconstrnt,COMMUNICATE *communicate,
                    VERLIST *verlist,int int_res_ter,
                    int iget_pv_real_inter, int iget_pe_real_inter) 

/*==========================================================================*/
    {/*begin routine*/
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

    int i,iii;
    int npairs_tmp=0,npairs_res_tmp=0;
    double tvten_temp,tv1,tv2;
    int    iter_shake_temp=0.0;
    int    iter_ratl_temp=0.0;
    double iter_23_temp=0.0;
    double iter_33_temp=0.0;
    double iter_46_temp=0.0;
    double iter_23r_temp=0.0;
    double iter_33r_temp=0.0;
    double iter_46r_temp=0.0;
    double vintert_temp=0.0;
    double vintrat_temp=0.0;
    double vpotnhc_temp=0.0;
    double vlong_temp=0.0;
    double vpot_v_temp=0.0;
    double kin_harm_temp=0.0;
    double pi_ke_prim=0.0;
    double pi_ke_vir=0.0;
    double press_inter_temp=0.0;
    double press_intra_temp=0.0;
    double pvten_temp[10];
    double vcoul_temp = 0.0;
    double vvdw_temp = 0.0;
    double vbondt_watts_temp = 0.0; 
    double vbendt_watts_temp = 0.0; 
    double vtot_watts_temp = 0.0; 
    double vbar_free_temp = 0.0;
    double vsurf_temp = 0.0;

    double vbondt_temp = 0.0;
    double vbendt_temp = 0.0;
    double vbend_bndt_temp = 0.0;   
    double vbend_bnd_bend_temp = 0.0;
    double vbend_bnd_bond_temp = 0.0;
    double vtorst_temp = 0.0;
    double vbond_free_temp = 0.0;
    double vbend_free_temp = 0.0;
    double vtors_free_temp = 0.0;
    double vonfot_temp = 0.0;
    double kinet_temp = 0.0;
    double kinet_nhc_temp = 0.0;
    double vpot_nhc_temp = 0.0;
    double *ptens_tvten = ptens->tvten;
    double *ptens_tvten_tmp = ptens->pvten_tmp;

    int np_forc        = communicate->np_forc;
    MPI_Comm comm_forc = communicate->comm_forc;
    MPI_Comm world     = communicate->world;

/*=======================================================================*/
/* I) List stuff */

    Reduce(&verlist->nver_lst_now,&npairs_tmp,1,MPI_INT,
               MPI_SUM,0,comm_forc);
    verlist->nver_lst_now      = npairs_tmp;
    if(int_res_ter==1){
     Reduce(&verlist->nver_lst_now_res,&npairs_res_tmp,1,
                                   MPI_INT,MPI_SUM,0,comm_forc);
     verlist->nver_lst_now_res  = npairs_res_tmp;
    }/*endif*/

/*=======================================================================*/
/* II) Energies */
  
   if(iget_pe_real_inter==1){
    Reduce(&(stat_avg->vintert), &vintert_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vcoul), &vcoul_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vvdw), &vvdw_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    stat_avg->vintert = vintert_temp;
    stat_avg->vcoul   = vcoul_temp;
    stat_avg->vvdw    = vvdw_temp;
   }/*endif*/


/*=======================================================================*/
/* II) Intra energies     */

    Reduce(&(stat_avg->vintrat), &vintrat_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbondt), &vbondt_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbendt), &vbendt_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbend_bndt), &vbend_bndt_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbend_bnd_bend), &vbend_bnd_bend_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbend_bnd_bond), &vbend_bnd_bond_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vtorst), &vtorst_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbond_free), &vbond_free_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbend_free), &vbend_free_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vtors_free), &vtors_free_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vonfot), &vonfot_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->kinet), &kinet_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbar_free), &vbar_free_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vsurft), &vsurf_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);

    Reduce(&(stat_avg->vbondt_watts), &vbondt_watts_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vbendt_watts), &vbendt_watts_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vtot_watts), &vtot_watts_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);

if((ensopts->nvt + ensopts->npt_i+ensopts->npt_f) > 0){
    Reduce(&(stat_avg->kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
    Reduce(&(stat_avg->vpotnhc), &vpotnhc_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
}/*endif*/

    stat_avg->vintrat = vintrat_temp;
    stat_avg->vbondt  = vbondt_temp;      
    stat_avg->vbendt  = vbendt_temp;      
    stat_avg->vbend_bndt     = vbend_bndt_temp;   
    stat_avg->vbend_bnd_bend = vbend_bnd_bend_temp;
    stat_avg->vbend_bnd_bond = vbend_bnd_bond_temp;
    stat_avg->vbondt_watts   = vbondt_watts_temp; 
    stat_avg->vbendt_watts   = vbendt_watts_temp; 
    stat_avg->vtot_watts     = vtot_watts_temp; 
    stat_avg->vtorst         = vtorst_temp;
    stat_avg->vsurft         = vsurf_temp;      
    stat_avg->vbond_free     = vbond_free_temp;
    stat_avg->vbend_free     = vbend_free_temp;
    stat_avg->vtors_free     = vtors_free_temp;
    stat_avg->vbar_free      = vbar_free_temp;
    stat_avg->vonfot         = vonfot_temp;
    stat_avg->kinet          = kinet_temp;
    stat_avg->vpotnhc = vpotnhc_temp;
    stat_avg->kinet_nhc = kinet_nhc_temp;

  if(iconstrnt==1){
   Reduce(&(stat_avg->iter_shake), &iter_shake_temp,1,MPI_INT,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_ratl), &iter_ratl_temp,1,MPI_INT,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_23), &iter_23_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_33), &iter_33_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_46), &iter_46_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_23r), &iter_23r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_33r), &iter_33r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
   Reduce(&(stat_avg->iter_46r), &iter_46r_temp,1,MPI_DOUBLE,
                   MPI_MAX,0,comm_forc);
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
/* IV) All reduce tvten     */

  if(np_forc > 1){
   Allreduce(&(ptens_tvten[1]), &(ptens_tvten_tmp[1]),9,MPI_DOUBLE,
                   MPI_SUM,0,comm_forc);
   for(i=1;i<=9;i++){
    ptens_tvten[i] = ptens_tvten_tmp[i];
   }/*endfor*/
  }/*endif*/

/*=======================================================================*/
/* IV) Construct the Pressure Tensor Quantities                          */

  if(iget_pv_real_inter==1){
    Allreduce(&(ptens->pvten_tot[1]), &(pvten_temp[1]),9,
              MPI_DOUBLE, MPI_SUM,0,world);
    for ( i = 1;i<=9;i++){
      ptens->pvten_tot[i] = pvten_temp[i];
    }

    Reduce(&(stat_avg->press_inter),&press_inter_temp,1,MPI_DOUBLE,
                 MPI_SUM,0,comm_forc);
    stat_avg->press_inter = press_inter_temp;

    Reduce(&(stat_avg->press_intra),&press_intra_temp,1,MPI_DOUBLE,
                 MPI_SUM,0,comm_forc);
    stat_avg->press_intra = press_intra_temp;
  }/*endif*/

/*=======================================================================*/
  }/*end routine*/
/*=======================================================================*/










