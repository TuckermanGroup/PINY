 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Routine :   zero_sys.c                                 */
 /*                                                                          */
 /* The routine to zero all the integer elements of structure class         */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_parse_local.h"

  void zero_sys(CLASS *class,GENERAL_DATA *general_data)

{/*begin routine*/

/* zero cell->cell_info  */
  general_data->cell.iperd       = 0;     
  general_data->cell.intra_perds = 0;

/* zero general_data->stat_avg  */

  general_data->stat_avg.vintert = 0;
  general_data->stat_avg.aivintert = 0;
  general_data->stat_avg.avintert = 0;
  general_data->stat_avg.vintrat = 0;
  general_data->stat_avg.aivintrat = 0;
  general_data->stat_avg.avintrat = 0;
  general_data->stat_avg.vbondt = 0;
  general_data->stat_avg.vbendt = 0;
  general_data->stat_avg.vtorst = 0;
  general_data->stat_avg.vonfot = 0;
  general_data->stat_avg.vbend_bndt = 0;
  general_data->stat_avg.vbend_bnd_bond = 0;
  general_data->stat_avg.vbend_bnd_bend = 0;
  general_data->stat_avg.vreal = 0;
  general_data->stat_avg.vrecip = 0;
  general_data->stat_avg.vvdw = 0;
  general_data->stat_avg.vcoul = 0;
  general_data->stat_avg.vlong = 0;
  general_data->stat_avg.vbond_free = 0;
  general_data->stat_avg.vbend_free = 0;
  general_data->stat_avg.vtors_free = 0;
  general_data->stat_avg.kinet = 0;
  general_data->stat_avg.aikinet = 0;
  general_data->stat_avg.akinet = 0;
  general_data->stat_avg.kinet_v = 0;
  general_data->stat_avg.aikinet_v = 0;
  general_data->stat_avg.akinet_v = 0;
  general_data->stat_avg.vol = 0;
  general_data->stat_avg.aivol = 0;
  general_data->stat_avg.avol = 0;
  general_data->stat_avg.kinet_nhc = 0;
  general_data->stat_avg.aikinet_nhc = 0;
  general_data->stat_avg.akinet_nhc = 0;
  general_data->stat_avg.kinet_nhc_bead = 0;
  general_data->stat_avg.aikinet_nhc_bead = 0;
  general_data->stat_avg.akinet_nhc_bead = 0;
  general_data->stat_avg.vpot_v = 0;
  general_data->stat_avg.vpotnhc = 0;
  general_data->stat_avg.aiter_shake = 0;
  general_data->stat_avg.aiter_ratl = 0;
  general_data->stat_avg.iter_23 = 0;
  general_data->stat_avg.iter_33 = 0;
  general_data->stat_avg.iter_46 = 0;
  general_data->stat_avg.iter_43 = 0;
  general_data->stat_avg.iter_21 = 0;
  general_data->stat_avg.aiter_23 = 0;
  general_data->stat_avg.aiter_33 = 0;
  general_data->stat_avg.aiter_46 = 0;
  general_data->stat_avg.aiter_21 = 0;
  general_data->stat_avg.aiter_43 = 0;
  general_data->stat_avg.iter_23r = 0;
  general_data->stat_avg.iter_33r = 0;
  general_data->stat_avg.iter_46r = 0;
  general_data->stat_avg.iter_43r = 0;
  general_data->stat_avg.iter_21r = 0;
  general_data->stat_avg.aiter_23r = 0;
  general_data->stat_avg.aiter_33r = 0;
  general_data->stat_avg.aiter_46r = 0;
  general_data->stat_avg.aiter_21r = 0;
  general_data->stat_avg.aiter_43r = 0;
  general_data->stat_avg.acella = 0;
  general_data->stat_avg.acellb = 0;
  general_data->stat_avg.acellc = 0;
  general_data->stat_avg.aicella = 0;
  general_data->stat_avg.aicellb = 0;
  general_data->stat_avg.aicellc = 0;
  general_data->stat_avg.acellab = 0;
  general_data->stat_avg.acellbc = 0;
  general_data->stat_avg.acellac = 0;
  general_data->stat_avg.aicellab = 0;
  general_data->stat_avg.aicellbc = 0;
  general_data->stat_avg.aicellac = 0;
  general_data->stat_avg.apress = 0;
  general_data->stat_avg.aipress = 0;
  general_data->stat_avg.press_inter = 0;
  general_data->stat_avg.press_intra = 0;
  general_data->stat_avg.apress_inter = 0;
  general_data->stat_avg.aipress_inter = 0;
  general_data->stat_avg.apress_intra = 0;
  general_data->stat_avg.aipress_intra = 0;
  general_data->stat_avg.press_kin = 0;
  general_data->stat_avg.apress_kin = 0;
  general_data->stat_avg.aipress_kin = 0;
  general_data->stat_avg.econv0 = 0;
  general_data->stat_avg.econv = 0;
  general_data->stat_avg.cpu1 = 0;
  general_data->stat_avg.cpu2 = 0;
  general_data->stat_avg.acpu = 0;
  general_data->stat_avg.cpu_now = 0;
  general_data->stat_avg.updates = 0;
  general_data->stat_avg.updates_w = 0;
  general_data->stat_avg.kinet_cp = 0;
  general_data->stat_avg.aikinet_cp = 0;
  general_data->stat_avg.akinet_cp = 0;
  general_data->stat_avg.kinet_nhc_cp = 0;
  general_data->stat_avg.aikinet_nhc_cp = 0;
  general_data->stat_avg.akinet_nhc_cp = 0;
  general_data->stat_avg.vpotnhc_cp = 0;
  general_data->stat_avg.cp_ehart = 0;
  general_data->stat_avg.aicp_ehart = 0;
  general_data->stat_avg.acp_ehart = 0;
  general_data->stat_avg.cp_eext = 0;
  general_data->stat_avg.aicp_eext = 0;
  general_data->stat_avg.acp_eext = 0;
  general_data->stat_avg.cp_exc = 0;
  general_data->stat_avg.cp_muxc = 0;
  general_data->stat_avg.aicp_exc = 0;
  general_data->stat_avg.acp_exc = 0;
  general_data->stat_avg.cp_eke = 0;
  general_data->stat_avg.aicp_eke = 0;
  general_data->stat_avg.acp_eke = 0;
  general_data->stat_avg.cp_enl = 0;
  general_data->stat_avg.aicp_enl = 0;
  general_data->stat_avg.acp_enl = 0;
  general_data->stat_avg.aiter_shake_cp = 0;
  general_data->stat_avg.aiter_ratl_cp = 0;
  general_data->stat_avg.maxfc = 0;
  general_data->stat_avg.maxf = 0;
  general_data->stat_avg.pi_ke_prim = 0;
  general_data->stat_avg.pi_ke_vir = 0;
  general_data->stat_avg.api_ke_prim = 0;
  general_data->stat_avg.api_ke_vir = 0;
  general_data->stat_avg.aipi_ke_prim = 0;
  general_data->stat_avg.aipi_ke_vir = 0;
  general_data->stat_avg.kin_harm = 0;
  general_data->stat_avg.akin_harm = 0;
  general_data->stat_avg.aikin_harm = 0;
  general_data->stat_avg.iter_shake = 0;


/* zero clatoms_info */
  class->clatoms_info.natm_tot           = 0;        
  class->clatoms_info.pi_beads           = 0;        
  class->clatoms_info.pi_beads_proc      = 0;   
  class->clatoms_info.nfree              = 0;              
  class->clatoms_info.nchrg              = 0;              
  class->clatoms_info.pi_beads_full_ter  = 0;  
  class->clatoms_info.pi_beads_res_ter   = 0;   
  class->clatoms_info.pi_beads_full_tra  = 0;  
  class->clatoms_info.pi_beads_res_tra   = 0;   
  class->clatoms_info.natm_mall          = 0;    
  class->clatoms_info.nchrg_mall         = 0;


/* zero element to ghost_atoms.ghost_atoms_info  */
    class->ghost_atoms.nghost_tot    = 0;        
    class->ghost_atoms.natm_comp_max = 0;     
    class->ghost_atoms.nghost_mall   = 0;
    class->ghost_atoms.nghost_old    = 0;
    class->ghost_atoms.ncomp_mall    = 0;  
    class->ghost_atoms.ncomp_old     = 0;  


/* zero elements of atommaps.atommaps_info  */
    class->atommaps.nmol_typ      = 0;            
    class->atommaps.nres_typ      = 0;            
    class->atommaps.natm_typ      = 0;            
    class->atommaps.nfreeze       = 0;             
    class->atommaps.nfreeze_mall  = 0;
    class->atommaps.nres_max      = 0;    
    class->atommaps.nres_typ_max  = 0;
    class->atommaps.nres_tot      = 0;
    class->atommaps.nres_sum      = 0;
    class->atommaps.natm_typ_mall = 0;

/* zero therm_info */
  class->therm_info_class.num_nhc   = 0;    
  class->therm_info_class.len_nhc   = 0;    
  class->therm_info_class.nres_nhc  = 0;   
  class->therm_info_class.nyosh_nhc = 0;  

  class->therm_info_bead.num_nhc   = 0;    
  class->therm_info_bead.len_nhc   = 0;    
  class->therm_info_bead.nres_nhc  = 0;   
  class->therm_info_bead.nyosh_nhc = 0;  

/* zero nbr_list.nbr_list_info */
    class->nbr_list.nolst            = 0;
    class->nbr_list.iver             = 0;                
    class->nbr_list.verlist.iver_init        = 0;           
    class->nbr_list.verlist.iver_fill        = 0;  
    class->nbr_list.verlist.iver_count       = 0;
    class->nbr_list.verlist.nolst_ver_update = 0;
    class->nbr_list.verlist.lnk_ver_update   = 0; 
    class->nbr_list.verlist.nver_lst         = 0;           
    class->nbr_list.verlist.nver_lst_res     = 0;       
    class->nbr_list.verlist.jver_pad         = 0;           
    class->nbr_list.verlist.nmem_min_lst     = 0;       
    class->nbr_list.lnklist.lnk_excl_safe    = 0;      
    class->nbr_list.lnklist.lnk_vol_safe     = 0;       
    class->nbr_list.ilnk             = 0;               
    class->nbr_list.lnklist.ilnk_init        = 0;          
    class->nbr_list.lnklist.lnk_for_odd      = 0;        
    class->nbr_list.lnklist.ncell_div_avg    = 0;     
    class->nbr_list.lnklist.ncell_a          = 0;
    class->nbr_list.lnklist.ncell_b          = 0;
    class->nbr_list.lnklist.ncell_c          = 0;
    class->nbr_list.lnklist.natm_cell_max    = 0;          
    class->nbr_list.lnklist.lnk_list_dim     = 0;           
    class->nbr_list.lnklist.nshft_lnk        = 0;              
    class->nbr_list.lnklist.nshft_lnk_res    = 0;          

/* zero interact.interact_info  */
    class->interact.nter_typ        = 0;   
    class->interact.nsplin          = 0;     
    class->interact.nsplin_tot      = 0; 
    class->interact.iswit_vdw       = 0;  
    class->interact.dielectric_opt  = 0; 
    class->interact.ninter_mall     = 0;
    class->interact.nsplin_mall_tot = 0;
  
/* zero ewald.ewald_info  */
    general_data->ewald.nktot          = 0; 
    general_data->ewald.nktot_res      = 0; 
    general_data->ewald.nktot_mall     = 0;      
    general_data->ewald.nktot_res_mall = 0;  

/* zero part_mesh.part_mesh_info  */
    class->part_mesh.pme_on           = 0;                   
    class->part_mesh.kmax_pme         = 0;                 
    class->part_mesh.n_interp         = 0;                 
    class->part_mesh.nktot_pme        = 0;                
    class->part_mesh.ngrid_a          = 0;
    class->part_mesh.ngrid_b          = 0;
    class->part_mesh.ngrid_c          = 0;
    class->part_mesh.pme_res_on       = 0;
    class->part_mesh.kmax_pme_res     = 0;                 
    class->part_mesh.n_interp_res     = 0;                  
    class->part_mesh.nktot_pme_res    = 0;            
    class->part_mesh.ngrid_a_res      = 0;
    class->part_mesh.ngrid_b_res      = 0;
    class->part_mesh.ngrid_c_res      = 0;
    class->part_mesh.nlen_pme         = 0;                 
    class->part_mesh.ninterp_mall     = 0;    
    class->part_mesh.ninterp_tot_mall = 0;   

}/*end routine*/










