/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_utilities                                */
/*                                                                          */
/* This subprogram provides some integrator utility routines                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cp_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void get_cpke_pimd(CPCOEFFS_INFO *cpcoeffs_info,CPCOEFFS_POS *cpcoeffs_pos,
                      STAT_AVG *stat_avg, int cp_lsda, int pi_beads_proc,
                         int np_states)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int is,icoef,nstate_up,ncoef,nstate_dn,i,iii;

    double kinet_cp,kinet_cp_up,kinet_cp_dn;
    double *cpcoeffs_vcre_up;
    double *cpcoeffs_vcim_up;
    double *cpcoeffs_vcre_dn;
    double *cpcoeffs_vcim_dn;
    double *cpcoeffs_cmass    = cpcoeffs_info->cmass;
    int *cpcoeffs_ioff_up     = cpcoeffs_info->ioff_upt;
    int *cpcoeffs_ioff_dn     = cpcoeffs_info->ioff_dnt;
    int pi_beads              = cpcoeffs_info->pi_beads;
    int ip,ncoef_up,ncoef_dn;    
    int icmoff_up        = cpcoeffs_info->icoef_start_up-1;
    int icmoff_dn        = cpcoeffs_info->icoef_start_dn-1;

/*========================================================================*/
/* I) Get the CP kinetic energy                                           */


  if(np_states==1){
     ncoef_up     = cpcoeffs_info->ncoef;
     ncoef_dn     = cpcoeffs_info->ncoef;
  }else{
     ncoef_up     = cpcoeffs_info->nstate_ncoef_proc_up;
     ncoef_dn     = cpcoeffs_info->nstate_ncoef_proc_dn;
  }/*endif*/

   nstate_up = cpcoeffs_info->nstate_up;
   nstate_dn = cpcoeffs_info->nstate_dn;
   kinet_cp_up   = 0.0;
   for(ip=1;ip<=pi_beads_proc;ip++){
     cpcoeffs_vcre_up  = cpcoeffs_pos[ip].vcre_up;
     cpcoeffs_vcim_up  = cpcoeffs_pos[ip].vcim_up;
     for(is=1;is<=nstate_up;is++) {
       for(i=1;i<=ncoef_up;i++) {
         icoef = i+cpcoeffs_ioff_up[is];
         kinet_cp_up += (cpcoeffs_vcre_up[icoef]*cpcoeffs_vcre_up[icoef]
                        *cpcoeffs_cmass[(i+icmoff_up)]);
         kinet_cp_up += (cpcoeffs_vcim_up[icoef]*cpcoeffs_vcim_up[icoef]
                        *cpcoeffs_cmass[(i+icmoff_up)]);
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/

   kinet_cp_dn = 0.0;
   if( (cp_lsda == 1) &&(nstate_dn != 0) ){
    for(ip=1;ip<=pi_beads_proc;ip++){
     cpcoeffs_vcre_dn  = cpcoeffs_pos[ip].vcre_dn;
     cpcoeffs_vcim_dn  = cpcoeffs_pos[ip].vcim_dn;
     for(is=1;is<=nstate_dn;is++) {
       for(i=1;i<=ncoef_dn;i++) {
         icoef = i+cpcoeffs_ioff_dn[is];
         kinet_cp_dn += (cpcoeffs_vcre_dn[icoef]*cpcoeffs_vcre_dn[icoef]
                        *cpcoeffs_cmass[i+icmoff_dn]);
         kinet_cp_dn += (cpcoeffs_vcim_dn[icoef]*cpcoeffs_vcim_dn[icoef]
                        *cpcoeffs_cmass[i]);
       }/*endfor*/
      }/*endfor*/
     }/*endfor*/
   }/* endif */

    kinet_cp_up /= 2.0;
    kinet_cp_dn /= 2.0;
    stat_avg->kinet_cp_up   = kinet_cp_up/(double)pi_beads;
    stat_avg->kinet_cp_dn   = kinet_cp_dn/(double)pi_beads;
    stat_avg->kinet_cp      = (kinet_cp_up + kinet_cp_dn)/(double)pi_beads;

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void nhc_cp_potkin_pimd(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                         STAT_AVG *stat_avg, int pi_beads)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

int inhc,ichain,ip,iii;
int len_c_nhc,num_c_nhc;
double kinet_nhc,vpotnhc;
double **cptherm_c_nhc;
double **cptherm_vc_nhc;
double **cptherm_cmass_nhc = cptherm_info->cmass_nhc;
double **cptherm_c_gkt     = cptherm_info->c_gkt;

/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc     = 0.0;
    kinet_nhc   = 0.0;
    len_c_nhc   = cptherm_info->len_c_nhc;
    num_c_nhc   = cptherm_info->num_c_nhc;
    for(ip=1;ip<=pi_beads;ip++){
      cptherm_c_nhc     = cptherm_pos[ip].c_nhc;
      cptherm_vc_nhc    = cptherm_pos[ip].vc_nhc;
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          kinet_nhc += (cptherm_cmass_nhc[ichain][inhc]
                    *cptherm_vc_nhc[ichain][inhc]*cptherm_vc_nhc[ichain][inhc]);
          vpotnhc += (cptherm_c_gkt[ichain][inhc]*
                               cptherm_c_nhc[ichain][inhc]);   
        /*endfor*/}
      /*endfor*/}
    /*endfor*/}
    kinet_nhc /= 2.0;

    stat_avg->vpotnhc_cp   = vpotnhc;
    stat_avg->kinet_nhc_cp = kinet_nhc;

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 
void nhc_cp_potkin_massiv_pimd(CPTHERM_INFO *cptherm_info,CPTHERM_POS *cptherm_pos,
                               STAT_AVG *stat_avg, int pi_beads)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

int inhc,ichain,ip,iii;
int len_c_nhc,num_c_nhc;
double kinet_nhc,vpotnhc;
double cptherm_c_nhc;
double **cptherm_vc_nhc;
double cptherm_cmass_nhc = cptherm_info->cmass_nhc_massiv;
double cptherm_c_gkt     = cptherm_info->c_gkt_massiv;

/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc     = 0.0;
    kinet_nhc   = 0.0;
    len_c_nhc   = cptherm_info->len_c_nhc;
    num_c_nhc   = cptherm_info->num_c_nhc;
    for(ip=1;ip<=pi_beads;ip++){
      cptherm_c_nhc     = cptherm_pos[ip].c_nhc_massiv;
      cptherm_vc_nhc    = cptherm_pos[ip].vc_nhc;
      for(ichain=1;ichain<=len_c_nhc;ichain++){
        for(inhc=1;inhc<=num_c_nhc;inhc++){
          kinet_nhc += (cptherm_cmass_nhc
                    *cptherm_vc_nhc[ichain][inhc]*cptherm_vc_nhc[ichain][inhc]);
        }/*endfor*/
        vpotnhc += (cptherm_c_gkt*cptherm_c_nhc);   
      }/*endfor*/
    }/*endfor*/
    kinet_nhc /= 2.0;
    stat_avg->vpotnhc_cp   = vpotnhc;
    stat_avg->kinet_nhc_cp = kinet_nhc;

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/






