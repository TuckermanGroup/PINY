/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_properties                              */
/*                                                                          */
/* This subprogram allocates memory for the calculation of electronic       */ 
/* properties (ELF, Wannier orbitals, etc.)                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_properties(CP *cp)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  int nfft2_up,nfft2_dn;
  int nfft2_mall_up,nfft2_mall_dn;
  int nfft2_mall_up_proc,nfft2_mall_dn_proc;
  int nfft_dn,nfft2_dn_proc;
  int nfft2_up_proc;
  int nfft_up_proc_dens_cp_box,nfft_up_dens_cp_box;
  int nfft2_up_proc_dens_cp_box,nfft2_up_dens_cp_box;
  int nfft_dn_proc_dens_cp_box,nfft_dn_dens_cp_box;
  int nfft2_dn_proc_dens_cp_box,nfft2_dn_dens_cp_box;
  
  int cp_lsda             = cp->cpopts.cp_lsda;
  int nfft_up_proc        = cp->cp_para_fft_pkg3d_lg.nfft_proc;
  int nfft_up             = cp->cp_para_fft_pkg3d_lg.nfft;
  int cp_dual_grid_opt_on = cp->cpopts.cp_dual_grid_opt;

/*=========================================================================*/
/* I) Malloc size calculation */

 /*-------------------------------------------------------------------------*/
 /* i) Dual grid CP : Define the small dense grid sizes */


  if(cp_dual_grid_opt_on >= 1){

    nfft_up_proc_dens_cp_box   = cp->cp_para_fft_pkg3d_dens_cp_box.nfft_proc;
    nfft_up_dens_cp_box        = cp->cp_para_fft_pkg3d_dens_cp_box.nfft;
    nfft2_up_dens_cp_box       = nfft_up_dens_cp_box/2;

    nfft_dn_proc_dens_cp_box   = (cp_lsda == 1 ? nfft_up_proc_dens_cp_box : 0);
    nfft_dn_dens_cp_box        = (cp_lsda == 1 ? nfft_up_dens_cp_box : 0);
    nfft2_dn_dens_cp_box       = nfft_dn_dens_cp_box/2;

  }/*endif cp_dual_grid_opt_on */

 /*-------------------------------------------------------------------------*/
 /* ii) Normal CP : Define the grid size              */
 /*     Dual   CP : Define the large sparse grid size */

 nfft2_up      = nfft_up/2;
 nfft2_up_proc = nfft_up_proc/2;

 nfft_dn       = (cp_lsda == 1 ? nfft_up : 0);
 nfft2_dn      = (cp_lsda == 1 ? nfft2_up : 0);
 nfft2_dn_proc = (cp_lsda == 1 ? nfft2_up_proc : 0);

 /*-------------------------------------------------------------------------*/
 /* iii) Choose the correct size for your application                       */
 /*      This is always the small dense grid                                */

  nfft2_mall_up      = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_up_dens_cp_box/2 : nfft2_up);
  nfft2_mall_up_proc = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_up_proc_dens_cp_box/2 : nfft2_up_proc);

  nfft2_mall_dn      = (cp_dual_grid_opt_on >= 1 ?
                        nfft_dn_dens_cp_box/2 : nfft2_dn);
  nfft2_mall_dn_proc = (cp_dual_grid_opt_on >= 1 ? 
                        nfft_dn_proc_dens_cp_box/2 : nfft2_dn_proc);

 /*-------------------------------------------------------------------------*/
 /* iv) Allocate the electron localization function                         */

  cp->electronic_properties.cp_elf_up = 
               (double *) cmalloc(nfft2_mall_up_proc*sizeof(double))-1;

  cp->electronic_properties.cp_elf_dn = 
               (double *) cmalloc(nfft2_mall_dn_proc*sizeof(double))-1;

/*========================================================================*/
} /* end routine */
/*==========================================================================*/




