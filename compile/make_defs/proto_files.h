#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
#
#  Abreviations for include files  PI_MD.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================

#-------------------------------------------------------
#              STANDARD INCLUDE FILE

STANDARD      = $(INCLUDES)/standard_include.h
DEFINES       = $(CODE)/include/typ_defs/defines.h 

#-------------------------------------------------------
#              TYPEDEFS 

TYP_GEN       = $(CODE)/include/typ_defs/typedefs_gen.h
TYP_CLASS     = $(CODE)/include/typ_defs/typedefs_class.h
TYP_PAR       = $(CODE)/include/typ_defs/typedefs_par.h
TYP_BND       = $(CODE)/include/typ_defs/typedefs_bnd.h
TYP_CP        = $(CODE)/include/typ_defs/typedefs_cp.h
TYP_STAT      = $(CODE)/include/typ_defs/typedefs_stat.h

#-------------------------------------------------------
#               PROTO_TYPES

ANAL_CP_ENT   = $(CODE)/analysis/proto_analysis_cp_entry.h
ANAL_MD_ENT   = $(CODE)/analysis/proto_analysis_md_entry.h
ANAL_LOC_ENT  = $(CODE)/analysis/proto_analysis_local_entry.h
COORD_ENT     = $(CODE)/interface/coords/proto_coords_entry.h
COORD_LOC     = $(CODE)/interface/coords/proto_coords_local.h
COORD_CP_ENT  = $(CODE)/interface/coords_cp/proto_coords_cp_entry.h
COORD_CP_LOC  = $(CODE)/interface/coords_cp/proto_coords_cp_local.h
CPEWALD_ENT   = $(CODE)/interface/cp_ewald/proto_cp_ewald_entry.h
CPEWALD_LOC   = $(CODE)/interface/cp_ewald/proto_cp_ewald_local.h
ENR_CP_LOC    = $(CODE)/energy/cp/proto_energy_cp_local.h
ENR_CP_ENT    = $(CODE)/energy/cp/proto_energy_cp_entry.h
ENR_CTRL_LOC  = $(CODE)/energy/control/proto_energy_ctrl_local.h
ENR_CTRL_ENT  = $(CODE)/energy/control/proto_energy_ctrl_entry.h
ENR_CTRL_CP_LOC = $(CODE)/energy/control_cp/proto_energy_ctrl_cp_local.h
ENR_CTRL_CP_ENT = $(CODE)/energy/control_cp/proto_energy_ctrl_cp_entry.h
ENR_CPCON_ENT = $(CODE)/energy/cp_con/proto_energy_cpcon_entry.h
ENR_CPCON_LOC = $(CODE)/energy/cp_con/proto_energy_cpcon_local.h
ENR_PIMD_LOC  = $(CODE)/energy/pimd/proto_pimd_local.h
ENR_PIMD_ENT  = $(CODE)/energy/pimd/proto_pimd_entry.h
FRND_ENT      = $(CODE)/friend_lib/proto_friend_lib_entry.h
HANDLE_ENT    = $(CODE)/interface/handle/proto_handle_entry.h
HANDLE_LOC    = $(CODE)/interface/handle/proto_handle_local.h
INTRA_CON_ENT = $(CODE)/energy/intra_con/proto_intra_con_entry.h
INTRA_CON_LOC = $(CODE)/energy/intra_con/proto_intra_con_local.h
INTRA_ENT     = $(CODE)/energy/intra/proto_intra_entry.h
INT_CP_ENT    = $(CODE)/integrate/cp/proto_integrate_cp_entry.h
INT_CP_LOC    = $(CODE)/integrate/cp/proto_integrate_cp_local.h
INT_CPPIMD_ENT= $(CODE)/integrate/cppimd/proto_integrate_cppimd_entry.h
INT_CPPIMD_LOC= $(CODE)/integrate/cppimd/proto_integrate_cppimd_local.h
INT_CPMIN_ENT = $(CODE)/integrate/cpmin/proto_integrate_cpmin_entry.h
INT_CPMIN_LOC = $(CODE)/integrate/cpmin/proto_integrate_cpmin_local.h
INT_MD_ENT    = $(CODE)/integrate/md/proto_integrate_md_entry.h
INT_MD_LOC    = $(CODE)/integrate/md/proto_integrate_md_local.h
INT_MIN_ENT   = $(CODE)/integrate/min/proto_integrate_min_entry.h
INT_MIN_LOC   = $(CODE)/integrate/min/proto_integrate_min_local.h
INT_PIMD_ENT  = $(CODE)/integrate/pimd/proto_integrate_pimd_entry.h
INT_PIMD_LOC  = $(CODE)/integrate/pimd/proto_integrate_pimd_local.h
INTER_ENT     = $(CODE)/interface/inter_params/proto_inter_params_entry.h
INTER_LOC     = $(CODE)/interface/inter_params/proto_inter_params_local.h
INTRA_LOC     = $(CODE)/interface/intra_params/proto_intra_params_local.h
INTRA_ENT     = $(CODE)/interface/intra_params/proto_intra_params_entry.h
LISTS_ENT     = $(CODE)/interface/lists/proto_lists_entry.h
LISTS_LOC     = $(CODE)/interface/lists/proto_lists_local.h
MOL_ENT       = $(CODE)/interface/mol_params/proto_mol_params_entry.h 
MOL_LOC       = $(CODE)/interface/mol_params/proto_mol_params_local.h
MAIN_ENT      = $(CODE)/main/proto_main_entry.h
MAIN_LOC      = $(CODE)/main/proto_main_local.h
MAIN_CP_LOC   = $(CODE)/main/proto_main_cp_local.h
MATH          = $(CODE)/mathlib/proto_math.h
OUTPUT_ENT    = $(CODE)/output/proto_output_entry.h
OUTPUT_LOC    = $(CODE)/output/proto_output_local.h
OUTPUT_CP_ENT = $(CODE)/output/proto_output_cp_entry.h
OUTPUT_CP_LOC = $(CODE)/output/proto_output_cp_local.h
PARSE_ENT     = $(CODE)/interface/parse/proto_parse_entry.h
PARSE_LOC     = $(CODE)/interface/parse/proto_parse_local.h
REC_ENT       = $(CODE)/energy/inter/recip3d/proto_recip3d_entry.h
REC_LOC       = $(CODE)/energy/inter/recip3d/proto_recip3d_local.h
REAL_ENT      = $(CODE)/energy/inter/real_space/proto_real_space_entry.h
REAL_LOC      = $(CODE)/energy/inter/real_space/proto_real_space_local.h
SCRATCH_ENT   = $(CODE)/interface/scratch/proto_scratch_entry.h
SCRATCH_LOC   = $(CODE)/interface/scratch/proto_scratch_local.h
SEARCH_ENT    = $(CODE)/interface/search_base/proto_search_entry.h
SEARCH_LOC    = $(CODE)/interface/search_base/proto_search_local.h
SIM_LOC       = $(CODE)/interface/sim_params/proto_sim_params_local.h
SIM_ENT       = $(CODE)/interface/sim_params/proto_sim_params_entry.h
SMPL_CLAS_ENT = $(CODE)/interface/vel_sampl_class/proto_vel_sampl_class_entry.h
SMPL_CLAS_LOC = $(CODE)/interface/vel_sampl_class/proto_vel_sampl_class_local.h
SMPL_CP_ENT   = $(CODE)/interface/vel_sampl_cp/proto_vel_sampl_cp_entry.h
SMPL_CP_LOC   = $(CODE)/interface/vel_sampl_cp/proto_vel_sampl_cp_local.h
SURF_ENT      = $(CODE)/energy/surface/proto_surf_entry.h
SURF_PRMS_ENT = $(CODE)/interface/surf_params/proto_surf_params_entry.h
SURF_PRMS_LOC = $(CODE)/interface/surf_params/proto_surf_params_local.h
VPS_ENT       = $(CODE)/interface/vps_params/proto_vps_params_entry.h
VPS_LOC       = $(CODE)/interface/vps_params/proto_vps_params_local.h
COMM_ENT      = $(CODE)/communicate/proto_communicate_entry.h
COMM_WRAP     = $(CODE)/communicate/proto_communicate_wrappers.h
COMM_LOC      = $(CODE)/communicate/proto_communicate_local.h
WEIGH_NODE    = $(CODE)/interface/lists/weights_nodes_128.h
#-------------------------------------------------------
