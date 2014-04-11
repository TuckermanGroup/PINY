/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_class_list.c                            */
/*                                                                          */
/* Subprogram contains MPI utils and communication routines for interface   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_class_list(CLASS *class,GENERAL_DATA *general_data,
                            CLASS_PARSE *class_parse,int pimd_on)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int i,iii;
  int natm_mall  = class->clatoms_info.natm_mall;  
  int nktot      = general_data->ewald.nktot;
  int nmol_typ   = class->atommaps.nmol_typ;
  MPI_Comm world = class->communicate.world;

  comm_clatoms_list(&(class->clatoms_info),world);        Barrier(world);
  comm_ghost_atoms_list(&(class->ghost_atoms),world);     Barrier(world);
  comm_atommaps_list(&(class->atommaps),natm_mall,world); Barrier(world);
  comm_class_parse_list(class_parse,nmol_typ,world);      Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_clatoms_list(CLATOMS_INFO *clatoms_info,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int pi_beads   = clatoms_info->pi_beads;
  int nchrg_mall = clatoms_info->nchrg_mall;
  int natm_mall  = clatoms_info->natm_mall;
  int iii,i;

  Barrier(world);
  if(nchrg_mall>0){
    Bcast(&(clatoms_info->ichrg[1]),nchrg_mall,MPI_INT,0,world);
  }
  Bcast(&(clatoms_info->ip_lab[1]), pi_beads,MPI_INT,0,world);
  Barrier(world);

  Bcast(&(clatoms_info->cp_vlnc_up[1]), natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(clatoms_info->cp_vlnc_dn[1]), natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(clatoms_info->cp_atm_flag[1]), natm_mall,MPI_INT,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ghost_atoms_list(GHOST_ATOMS *ghost_atoms,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nghost_mall =  ghost_atoms->nghost_mall;
  int ncomp_mall  =  ghost_atoms->ncomp_mall; 

  if(nghost_mall>0){
    Bcast(&(ghost_atoms->ighost_map[1]),nghost_mall,MPI_INT,0,world);
    Barrier(world);
    Bcast(&(ghost_atoms->natm_comp[1]),nghost_mall,MPI_INT,0,world);
    Barrier(world);
    if(ncomp_mall>0){
      Bcast(&(ghost_atoms->iatm_comp[1][1]),(ncomp_mall*nghost_mall),
                                            MPI_INT,0,world);
      Barrier(world);
    }/*endif*/
  }/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_atommaps_list(ATOMMAPS *atommaps,int natm_mall,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nres_typ_max  = atommaps->nres_typ_max;
  int natm_typ_mall = atommaps->natm_typ_mall;
  int nfreeze_mall  =  atommaps->nfreeze_mall; 
  int nres_tot      = atommaps->nres_tot;
  int nmol_typ      = atommaps->nmol_typ;
  int nres_sum      = atommaps->nres_sum;
  int iii,i;

  Bcast(&(atommaps->nmol_jmol_typ[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->nres_1mol_jmol_typ[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->jatm_jmol_typ_strt[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->jres_jmol_typ_strt[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  if(nres_tot>0){
    Bcast(&(atommaps->ires_typ_jres_jmol_typ[1]),nres_tot,MPI_INT,0,world);
    Barrier(world);
    Bcast(&(atommaps->jatm_jres_1mol_jmol_typ_strt[1]),nres_tot,MPI_INT,0,world);
    Barrier(world);
  }/*endif*/
  Bcast(&(atommaps->natm_1mol_jmol_typ[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  if(nres_tot>0){
    Bcast(&(atommaps->natm_jres_jmol_typ[1]),nres_tot,MPI_INT,0,world);
    Barrier(world);
  }/*endif*/
  Bcast(&(atommaps->nfree_1mol_jmol_typ[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  if(nres_tot>0){
    Bcast(&(atommaps->nfree_jres_jmol_typ[1]),nres_tot,MPI_INT,0,world);
    Barrier(world);
  }/*endif*/
  Bcast(&(atommaps->icons_jmol_typ[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
  if(nres_tot>0){
    Bcast(&(atommaps->icons_jres_jmol_typ[1]),nres_tot,MPI_INT,0,world);
    Barrier(world);
  }/*endif*/
  Bcast(&(atommaps->iatm_mol_typ[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->iatm_res_typ[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->iatm_atm_typ[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->iatm_mol_num[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->iatm_res_num[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->ighost_flag[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  Bcast(&(atommaps->freeze_flag[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);
  if(nfreeze_mall>0){
    Bcast(&(atommaps->freeze_map[1]),nfreeze_mall,MPI_INT,0,world);
    Barrier(world);
  }
  Bcast(&(atommaps->atom_label[1]),natm_mall,MPI_INT,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_class_parse_list(CLASS_PARSE *class_parse,int nmol_typ,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(class_parse->imol_nhc_opt[1]),nmol_typ,MPI_INT,0,world);
  Bcast(&(class_parse->ionfo_opt[1]),nmol_typ,MPI_INT,0,world);
  Bcast(&(class_parse->ires_bond_conv[1]),nmol_typ,MPI_INT,0,world);
  Bcast(&(class_parse->mol_freeze_opt[1]),nmol_typ,MPI_INT,0,world);
  Bcast(&(class_parse->mol_hydrog_mass_opt[1]),nmol_typ,MPI_INT,0,world);
  Bcast(&(class_parse->mol_hydrog_con_opt[1]),nmol_typ,MPI_INT,0,world);
  Barrier(world);
 
/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/






