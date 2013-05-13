/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_class_data.c                            */
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

void communicate_class_data(CLASS *class,GENERAL_DATA *general_data,
                            CLASS_PARSE *class_parse,int pimd_on)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nmol_typ   = class->atommaps.nmol_typ;
  MPI_Comm world = class->communicate.world;

  comm_statepoint(&(general_data->statepoint),world);   Barrier(world);
  comm_clatoms_data(&(class->clatoms_info),nmol_typ,   
                                     pimd_on,world);    Barrier(world);
  comm_nbr_list_data(&(class->nbr_list),world);         Barrier(world);
  comm_interact_data(&(class->interact),world);         Barrier(world);
  comm_ghost_atoms_data(&(class->ghost_atoms),world);   Barrier(world);
  comm_class_parse_data(class_parse,nmol_typ,world);    Barrier(world);
  comm_surface_data(&(class->surface),world);           Barrier(world);

  Bcast(&(class->tot_memory),1,MPI_DOUBLE,0,world);     Barrier(world);

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_statepoint(STATEPOINT *statepoint,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  MPI_Datatype statepoint_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(statepoint->pext),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = 3;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&statepoint_comm);
  Barrier(world);
  Type_commit(&statepoint_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,statepoint_comm,0,world);
  Barrier(world);
  Type_free(&statepoint_comm);



/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cell_data(CELL *cell,double *dbox_rat,int *box_rat,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii,i;

  Barrier(world);
  Bcast(&(cell->vol),1,MPI_DOUBLE,0,world);
  Bcast(&(cell->vol_cp),1,MPI_DOUBLE,0,world);
  Bcast(&(cell->vol0),1,MPI_DOUBLE,0,world);
  Bcast(&(cell->cubic_box_flag),1,MPI_INT,0,world);
  Bcast(&(cell->hmat[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->hmati[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->hmat_ewd[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->hmat_ewd_cp[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->hmat_cp[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->hmati_cp[1]),9,MPI_DOUBLE,0,world);
  Bcast(&(cell->cp_box_center[1]),3,MPI_DOUBLE,0,world);
  Bcast(&(cell->cp_box_center_rel[1]),3,MPI_DOUBLE,0,world);
  Bcast(&(cell->cp_vbox_center[1]),3,MPI_DOUBLE,0,world);
  Bcast(&(cell->cp_fbox_center[1]),3,MPI_DOUBLE,0,world);

  Bcast(dbox_rat,1,MPI_DOUBLE,0,world);
  Bcast(box_rat,1,MPI_INT,0,world);

  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_clatoms_data(CLATOMS_INFO *clatoms_info,int nmol_typ,
                            int pimd_on,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  /* This one still needs some work! */

  int iii;
  int ndata_scal = 9;
  MPI_Datatype clatoms_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int natm_mall = clatoms_info->natm_mall;
  int pi_beads  = clatoms_info->pi_beads;

  Address(&(clatoms_info->gamma_adb),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = ndata_scal;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&clatoms_info_comm);
  Barrier(world);
  Type_commit(&clatoms_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,clatoms_info_comm,0,world);
  Type_free(&clatoms_info_comm);

  Barrier(world);
  Bcast(&(clatoms_info->mass[1]),natm_mall,MPI_DOUBLE,0,world);
  Bcast(&(clatoms_info->q[1]),natm_mall,MPI_DOUBLE,0,world);
  Bcast(&(clatoms_info->alp_pol[1]),natm_mall,MPI_DOUBLE,0,world);
  Bcast(&(clatoms_info->b_neut[1]),natm_mall,MPI_DOUBLE,0,world);
  Bcast(&(clatoms_info->text_atm[1]),natm_mall,MPI_DOUBLE,0,world);
  Bcast(&(clatoms_info->text_mol[1]),nmol_typ ,MPI_DOUBLE,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_nbr_list_data(NBR_LIST *nbr_list,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int nnbr_list_scal = 1;
  MPI_Datatype nbr_list_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(nbr_list->verlist.mem_safe),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = nnbr_list_scal;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&nbr_list_comm);
  Barrier(world);
  Type_commit(&nbr_list_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,nbr_list_comm,0,world);
  Type_free(&nbr_list_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_interact_data(INTERACT *interact,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ninteract_scal  = 9;
  int ninter_mall     = interact->ninter_mall;
  int nsplin_mall_tot = interact->nsplin_mall_tot;

  MPI_Datatype interact_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(interact->dielectric_rheal),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = ninteract_scal;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&interact_comm);
  Barrier(world);
  Type_commit(&interact_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,interact_comm,0,world);
  Type_free(&interact_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ghost_atoms_data(GHOST_ATOMS *ghost_atoms,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;
  int ncomp_mall  = ghost_atoms->ncomp_mall;
  int nghost_mall = ghost_atoms->nghost_mall;

  Barrier(world);
  if((ncomp_mall>0)&&(nghost_mall>0)){
    Bcast(&(ghost_atoms->coef[1][1]),
                  ncomp_mall*nghost_mall,MPI_DOUBLE,0,world);
  }
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_class_parse_data(CLASS_PARSE *class_parse,int nmol_typ,
                         MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int iii;

  Bcast(&(class_parse->tau_vol),1,MPI_DOUBLE,0,world);
  Bcast(&(class_parse->tau_vol_nhc),1,MPI_DOUBLE,0,world);
  Bcast(&(class_parse->tau_nhc_def),1,MPI_DOUBLE,0,world);
  Barrier(world);

if(nmol_typ > 0){
  Bcast(&(class_parse->tau_nhc_mol[1]),nmol_typ,MPI_DOUBLE,0,world);
  Bcast(&(class_parse->text_nhc_mol[1]),nmol_typ,MPI_DOUBLE,0,world);
  Bcast(&(class_parse->mol_hydrog_mass_val[1]),nmol_typ,MPI_DOUBLE,0,world);
  Barrier(world);
}/*endif*/

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_surface_data(SURFACE *surface, MPI_Comm world)

/*=======================================================================*/
/*      Begin routine                                                    */
     {/*begin routine */
/*=======================================================================*/
/*      Local variable declarations                                      */

#include "../typ_defs/typ_mask.h"

  MPI_Datatype surface_comm;
  MPI_Datatype types[1];
  MPI_Aint     displs[1];
  int          blockcounts[1];

  Address(&(surface->zheal),&displs[0]);
  types[0]       = MPI_DOUBLE;
  blockcounts[0] = 2;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&surface_comm);
  Barrier(world);
  Type_commit(&surface_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,surface_comm,0,world);
  Barrier(world);
  Type_free(&surface_comm);

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/



