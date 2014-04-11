 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Routine: zero_par.c                                    */
 /*                                                                          */
 /*  routine to zero elements in parse structures                            */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_parse_local.h"

 void zero_par(CLASS_PARSE *class_parse,NULL_INTER_PARSE *null_inter_parse)
{/*begin routine*/
              

/* zero class_parse->class_parse_info */

  class_parse->ivx_smpl           = 0;
  class_parse->ivnhc_smpl         = 0; 
  class_parse->kmax_ewd           = 0;
  class_parse->kmax_res           = 0;   
  class_parse->istart             = 0;              
  class_parse->ishift_pot         = 0;          
  class_parse->initial_spread_opt = 0;  
  class_parse->zero_com_vel       = 0;        



/* zero resbond_parse->resbond_parse_info  */
  /*  resbond_parse->ionfo         = 0;               
  resbond_parse->iconv         = 0;               
  resbond_parse->nresidue      = 0;
  resbond_parse->nresidue_max  = 0; 
  resbond_parse->nres_bond     = 0;
  resbond_parse->nres_bond_max = 0;*/

/* zero null_inter_parse->null_inter_parse_info */
  null_inter_parse->nbond_nul = 0;
  null_inter_parse->nbend_nul = 0;  
  null_inter_parse->ntors_nul = 0;
  null_inter_parse->nonfo_nul = 0;  


/* zero build_intra->build_intra_info */
  /*      build_intra->build_intra_info.natm_1res_pure_now  = 0;
      build_intra->natm_tot_max        = 0;
      build_intra->nfreeze_max         = 0;
      build_intra->nfreeze_now         = 0;
      build_intra->natm_1res_now       = 0;
      build_intra->natm_1res_max       = 0;
      build_intra->natmind_1res_now    = 0;
      build_intra->natmind_1res_max    = 0;
      build_intra->nghost_now          = 0;
      build_intra->nghost_tot_max      = 0;      
      build_intra->nbond_pow_max       = 0;       
      build_intra->nbond_con_max       = 0; 
      build_intra->nbond_nul_max       = 0;
      build_intra->nbond_typ_pow_max   = 0; 
      build_intra->nbond_typ_con_max   = 0;
      build_intra->nbend_pow_max       = 0;       
      build_intra->nbend_con_max       = 0; 
      build_intra->nbend_nul_max       = 0;
      build_intra->nbend_typ_pow_max   = 0; 
      build_intra->nbend_typ_con_max   = 0;
      build_intra->ntors_pow_max       = 0;       
      build_intra->ntors_con_max       = 0; 
      build_intra->ntors_nul_max       = 0;
      build_intra->ntors_typ_pow_max   = 0; 
      build_intra->ntors_typ_con_max   = 0;
      build_intra->nonfo_max           = 0;           
      build_intra->nonfo_nul_max       = 0;
      build_intra->nonfo_typ_max       = 0;
      build_intra->nbend_bnd_max       = 0;       
      build_intra->nbend_bnd_typ_max   = 0; 
      build_intra->ngrp_33_max         = 0;      
      build_intra->ngrp_21_max         = 0; 
      build_intra->ngrp_43_max         = 0; 
      build_intra->ngrp_23_max         = 0; 
      build_intra->ngrp_46_max         = 0;
      build_intra->ngrp_typ_21_max     = 0; 
      build_intra->ngrp_typ_43_max     = 0; 
      build_intra->ngrp_typ_33_max     = 0; 
      build_intra->ngrp_typ_23_max     = 0;
      build_intra->ngrp_typ_46_max     = 0;
      build_intra->natm_typ_max        = 0;
      build_intra->nres_typ_max        = 0;

      */
}/*end routine*/

