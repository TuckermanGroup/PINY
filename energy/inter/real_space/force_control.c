/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: force_control.c                                */
/*                                                                          */
/* These routine bundle the interactions and send them to force nopol       */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_real_space_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void force_control(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   FOR_SCR *for_scr,ATOMMAPS *atommaps,
                   CELL *cell,PTENS *ptens,INTERACT *interact,
                   ENERGY_CTRL *energy_ctrl, NBR_LIST *nbr_list,EXCL *excl,
                   INTRA_SCR *intra_scr,double *vreal,double *vvdw,
                   double *vcoul,int error_check_on,
                   CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

/*==========================================================================*/
{/*Begin Routine*/
/*=======================================================================*/
/* I) Nolist option */

  
    
  if((nbr_list->nolst)==1) {
    force_nolst(clatoms_info,clatoms_pos,
                for_scr,atommaps,cell,
              ptens,interact,energy_ctrl,nbr_list,excl,intra_scr,vreal,
              vvdw,vcoul,class_comm_forc_pkg);
  }/*endif*/

/*======================================================================*/
/* II) Verlist option */

  if(((nbr_list->iver)==1)&&((energy_ctrl->iget_full_inter)==1)){
    force_verlst(clatoms_info,clatoms_pos,
                 for_scr,atommaps,cell,
               ptens,interact,energy_ctrl,
               nbr_list->verlist.nter,nbr_list->verlist.jter,
               nbr_list->verlist.jver_off,intra_scr,vreal,vvdw,vcoul); 
  }/*endif*/
  if(((nbr_list->iver)==1)&&((energy_ctrl->iget_full_inter)==0)
     &&((energy_ctrl->iget_res_inter)==1)){
    force_verlst(clatoms_info,clatoms_pos,
                 for_scr,atommaps,cell,
               ptens,interact,energy_ctrl,
               nbr_list->verlist.nter_res,nbr_list->verlist.jter_res,
               nbr_list->verlist.jver_off_res,intra_scr,vreal,vvdw,vcoul);
  }/*endif*/

/*=====================================================================*/
/* III) Link list option */

  if((nbr_list->ilnk==1)&&(energy_ctrl->iget_full_inter==1)){
    force_lnklst(clatoms_info,clatoms_pos,
              for_scr,atommaps,cell,ptens,interact,energy_ctrl,
              (nbr_list->lnklist.ncell_a),(nbr_list->lnklist.ncell_b), 
              (nbr_list->lnklist.ncell_c),
              (nbr_list->lnklist.natm_cell_max),
              (nbr_list->lnklist.lnk_list),
             (nbr_list->lnklist.nshft_lnk),(nbr_list->lnklist.shft_wght),
             (nbr_list->lnklist.ishft_a),(nbr_list->lnklist.ishft_b),
              (nbr_list->lnklist.ishft_c),
             (nbr_list->lnklist.iexl_chk),excl,intra_scr,vreal,vvdw,vcoul,
             class_comm_forc_pkg);
  }/*endif*/
  if( ((nbr_list->ilnk)==1)&&((energy_ctrl->iget_full_inter)==0)
     &&((energy_ctrl->iget_res_inter)==1)){
     force_lnklst(clatoms_info,clatoms_pos,
       for_scr,atommaps,cell,ptens,interact,energy_ctrl,
       (nbr_list->lnklist.ncell_a),(nbr_list->lnklist.ncell_b),
       (nbr_list->lnklist.ncell_c),
       (nbr_list->lnklist.natm_cell_max),(nbr_list->lnklist.lnk_list),
       (nbr_list->lnklist.nshft_lnk_res),(nbr_list->lnklist.shft_wght_res),
       (nbr_list->lnklist.ishft_a_res),(nbr_list->lnklist.ishft_b_res),
       (nbr_list->lnklist.ishft_c_res),(nbr_list->lnklist.iexl_chk_res),
       excl,intra_scr,vreal,vvdw,vcoul,class_comm_forc_pkg);
  }/*endif*/

/*---------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/






