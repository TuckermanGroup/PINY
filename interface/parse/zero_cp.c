 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Routine: zero_cp.c                                     */
 /*                                                                          */
 /* The routine to zero all the integer elements of structure cp             */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_parse_local.h"

 void zero_cp(CP *cp)
{/*begin routine*/

/* zero elements of  cpcoeffs_info  */
	
  cp->cpcoeffs_info.ncoef                = 0;           
  cp->cpcoeffs_info.ncoef_l              = 0;         
  cp->cpcoeffs_info.pi_beads             = 0;        
  cp->cpcoeffs_info.nstate_up            = 0;
  cp->cpcoeffs_info.nstate_dn            = 0;
  cp->cpcoeffs_info.icmass_unif          = 0;         
  cp->cpcoeffs_info.ks_rot_on            = 0;           
  cp->cpcoeffs_info.n_ks_rot             = 0;               

/* zero cptherm  */
  cp->cptherm_info.num_c_nhc       = 0;          
  cp->cptherm_info.len_c_nhc       = 0;          
  cp->cptherm_info.nres_c_nhc      = 0;         
  cp->cptherm_info.nyosh_c_nhc     = 0;        

/* zero cpewald->cpewald_info */
  cp->cpewald.nktot_sm         = 0;   
  cp->cpewald.nktot_cp_mall    = 0;   
  cp->cpewald.nktot_cp_sm_mall = 0;


/* zero pseudo->pseudo_info  */
  cp->pseudo.n_ang_max     = 0;       
  cp->pseudo.nsplin_g      = 0;        
  cp->pseudo.nsplin_g_tot  = 0;    
  cp->pseudo.num_nl_lst    = 0;      
  cp->pseudo.nl_cut_on     = 0;       
  cp->pseudo.n_ang_mall    = 0;      
  cp->pseudo.nsplin_g_mall = 0;   
  cp->pseudo.nvpsnorm_mall = 0;       
  cp->pseudo.nlist_mall    = 0;       



/* zero cpscr->cpscr_info */
  cp->cpscr.cpscr_nonloc.natm_nls_max = 0;  
  cp->cpscr.cpscr_nonloc.nlscr_up     = 0;      
  cp->cpscr.cpscr_nonloc.nlscr_dn     = 0;      

}/*endroutine*/

