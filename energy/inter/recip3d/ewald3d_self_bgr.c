/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: ewald_selfbgr.c                                */
/*                                                                          */
/* Gets the self and background contributions from the Ewald sum            */
/* and the corresponding contributions to the pressure tensor               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_math.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void ewald3d_selfbgr(CLATOMS_INFO *clatoms_info,
                     EWALD *ewald, PTENS *ptens,double vol,
                     double wght_tra_res, double *vself,double *vbgr,
                     int np_forc,int iget_pv_real_inter,int iperd)

/*========================================================================*/
/*             Begin Routine                                              */
     {/*Begin Routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int step=1,iii;
  double self,bgr;

/* Define local pointers                                                  */

  int pi_beads            = clatoms_info->pi_beads;
  int natm_tot            = clatoms_info->natm_tot;
  double *clatoms_q       = clatoms_info->q;
  double *ptens_pvten_tot = ptens->pvten_tot;
  double *ptens_pvten     = ptens->pvten;
  double alp_ewd          = ewald->alp_ewd;
  double self_erf         = ewald->self_erf;

/*========================================================================*/
/* I) Self and background */

  self     = ddot1(natm_tot,clatoms_q,step,clatoms_q,step);
  self     = -self*(alp_ewd)*(self_erf)/sqrt(M_PI);
  (*vself) = self/((double)np_forc);

  bgr      = dsum1(natm_tot,clatoms_q,step);
  bgr      = -0.5*bgr*bgr*M_PI/(alp_ewd*alp_ewd*vol);
  (*vbgr)  = bgr/((double)np_forc);

/*========================================================================*/
/* II) Pressure */

  if( iperd <= 3 ){

    if(iget_pv_real_inter==1){
      ptens_pvten_tot[1] += bgr*pi_beads/((double)np_forc);
      ptens_pvten_tot[5] += bgr*pi_beads/((double)np_forc);
      ptens_pvten_tot[9] += bgr*pi_beads/((double)np_forc);
    }/*endif*/

    ptens_pvten[1] += bgr*wght_tra_res*pi_beads;
    ptens_pvten[5] += bgr*wght_tra_res*pi_beads;
    ptens_pvten[9] += bgr*wght_tra_res*pi_beads;

  }/*endif*/

/*------------------------------------------------------------------------*/
   }/*end routine */
/*========================================================================*/
