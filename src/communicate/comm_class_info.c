/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_class_info.c                            */
/*                                                                          */
/* Subprogram contains MPI utils and communication routines for interface   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_class_info(CLASS *class,GENERAL_DATA *general_data,
                                               CLASS_PARSE *class_parse)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  int myid       = class->communicate.myid;
  MPI_Comm world = class->communicate.world;


  comm_cell_info(&(general_data->cell),world);              Barrier(world);
  comm_communicate_info(&(class->communicate),world);       Barrier(world);
  comm_clatoms_info(&(class->clatoms_info),world);          Barrier(world);
  comm_ghost_info(&(class->ghost_atoms),world);             Barrier(world);
  comm_atommaps_info(&(class->atommaps),world);             Barrier(world);
  comm_therm_info_class(&(class->therm_info_class),world);  Barrier(world);
  comm_therm_info_bead(&(class->therm_info_bead),world);    Barrier(world);
  comm_nbr_list_info(&(class->nbr_list),world);             Barrier(world);
  comm_verlist_info(&(class->nbr_list.verlist),world);      Barrier(world);
  comm_lnklist_info(&(class->nbr_list.lnklist),world);      Barrier(world);
  comm_interact_info(&(class->interact),world);             Barrier(world);
  comm_ewald_info(&(general_data->ewald),world);            Barrier(world);
  comm_pme_info(&(class->part_mesh),world);                 Barrier(world);
  comm_timeinfo(&(general_data->timeinfo),world);           Barrier(world);
  comm_vel_samp_class(&(class->vel_samp_class),world,myid); Barrier(world);
  comm_energy_ctrl(&(class->energy_ctrl),world);            Barrier(world);
  comm_class_parse_info(class_parse,world);                 Barrier(world);
  comm_stat_avg_info(&(general_data->stat_avg),world);      Barrier(world);
  comm_surface_info(&(class->surface),world);               Barrier(world);

  Bcast(&(general_data->baro.len_nhc),1,MPI_DOUBLE,0,world);
  Bcast(&(general_data->pme_fft_pkg.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(general_data->pme_res_fft_pkg.igeneric_opt),1,MPI_INT,0,world);
  Barrier(world);   

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cell_info(CELL *cell,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  MPI_Datatype cell_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int ncell = 6;

  Address(&(cell->iperd),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ncell;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cell_info_comm);
  Barrier(world);
  Type_commit(&cell_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cell_info_comm,0,world);
  Barrier(world);
  Type_free(&cell_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_communicate_info(COMMUNICATE *communicate,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  MPI_Datatype communicate_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int ncommunicate = 5;

  Address(&(communicate->np_beads),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ncommunicate;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&communicate_info_comm);
  Barrier(world);
  Type_commit(&communicate_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,communicate_info_comm,0,world);
  Barrier(world);
  Type_free(&communicate_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_clatoms_info(CLATOMS_INFO *clatoms_info,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ninfo = 17;
  MPI_Datatype clatoms_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(clatoms_info->natm_tot),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ninfo;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&clatoms_info_comm);
  Barrier(world);
  Type_commit(&clatoms_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,clatoms_info_comm,0,world);
  Barrier(world);
  Type_free(&clatoms_info_comm);
  Barrier(world);
  Bcast(&clatoms_info->mass_sc_fact,1,MPI_DOUBLE,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ghost_info(GHOST_ATOMS *ghost_atoms,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nghost = 6;
  MPI_Datatype ghost_atoms_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(ghost_atoms->nghost_tot),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nghost;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&ghost_atoms_info_comm);
  Barrier(world);
  Type_commit(&ghost_atoms_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,ghost_atoms_info_comm,0,world);
  Barrier(world);
  Type_free(&ghost_atoms_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_atommaps_info(ATOMMAPS *atommaps,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int natommaps = 11;
  MPI_Datatype atommaps_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(atommaps->nmol_typ),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = natommaps;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&atommaps_info_comm);
  Barrier(world);
  Type_commit(&atommaps_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,atommaps_info_comm,0,world);
  Barrier(world);
  Type_free(&atommaps_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_therm_info_class(THERM_INFO *therm_info_class,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ntherm_class = 4;
  MPI_Datatype therm_info_class_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(therm_info_class->len_nhc),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ntherm_class;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&therm_info_class_comm);
  Barrier(world);
  Type_commit(&therm_info_class_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,therm_info_class_comm,0,world);
  Barrier(world);
  Type_free(&therm_info_class_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_therm_info_bead(THERM_INFO *therm_info_bead,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ntherm_bead = 3;
  int iii;
  MPI_Datatype therm_info_bead_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(therm_info_bead->len_nhc),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ntherm_bead;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&therm_info_bead_comm);
  Barrier(world);
  Type_commit(&therm_info_bead_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,therm_info_bead_comm,0,world);
  Barrier(world);
  Type_free(&therm_info_bead_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_nbr_list_info(NBR_LIST *nbr_list,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nlist_int = 4;
  MPI_Datatype nbr_list_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(nbr_list->iver),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nlist_int;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&nbr_list_info_comm);
  Barrier(world);
  Type_commit(&nbr_list_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,nbr_list_info_comm,0,world);
  Barrier(world);
  Type_free(&nbr_list_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_verlist_info(VERLIST *verlist,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nlist_int = 7;
  MPI_Datatype verlist_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(verlist->iver_init),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nlist_int;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&verlist_info_comm);
  Barrier(world);
  Type_commit(&verlist_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,verlist_info_comm,0,world);
  Barrier(world);
  Type_free(&verlist_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_lnklist_info(LNKLIST *lnklist,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nlist_int = 5;
  MPI_Datatype lnklist_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(lnklist->ilnk_init),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nlist_int;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&lnklist_info_comm);
  Barrier(world);
  Type_commit(&lnklist_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,lnklist_info_comm,0,world);
  Barrier(world);
  Type_free(&lnklist_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_interact_info(INTERACT *interact,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nlist = 4;
  MPI_Datatype interact_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(interact->nsplin),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nlist;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&interact_info_comm);
  Barrier(world);
  Type_commit(&interact_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,interact_info_comm,0,world);
  Barrier(world);
  Type_free(&interact_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ewald_info(EWALD *ewald,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  int newald = 1;
  MPI_Datatype ewald_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(ewald->nsplin_g),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = newald;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&ewald_info_comm);
  Barrier(world);
  Type_commit(&ewald_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,ewald_info_comm,0,world);
  Barrier(world);
  Type_free(&ewald_info_comm);
  Barrier(world);
  Bcast(&ewald->alp_ewd,1,MPI_DOUBLE,0,world);
  Bcast(&(ewald->ecut_clus),1,MPI_DOUBLE,0,world);
  Bcast(&(ewald->alp_clus),1,MPI_DOUBLE,0,world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_pme_info(PART_MESH *part_mesh,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int npart_mesh = 8;
  MPI_Datatype part_mesh_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(part_mesh->pme_on),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = npart_mesh;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&part_mesh_info_comm);
  Barrier(world);
  Type_commit(&part_mesh_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,part_mesh_info_comm,0,world);
  Barrier(world);
  Type_free(&part_mesh_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_timeinfo(TIMEINFO *timeinfo,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nint = 12;
  MPI_Datatype timeinfo_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(timeinfo->ntime),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nint;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&timeinfo_comm);
  Barrier(world);
  Type_commit(&timeinfo_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,timeinfo_comm,0,world);
  Barrier(world);
  Type_free(&timeinfo_comm);
  Barrier(world);
  Bcast(&(timeinfo->dt),1,MPI_DOUBLE,0,world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_vel_samp_class(VEL_SAMP_CLASS *vel_samp_class,MPI_Comm world,
                         int myid)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int i,itemp;
  double temp,qseed;
  int nint = 7;
  MPI_Datatype vel_samp_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(vel_samp_class->ivel_smpl_on),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nint;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&vel_samp_comm);
  Barrier(world);
  Type_commit(&vel_samp_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,vel_samp_comm,0,world);
  Barrier(world);
  Type_free(&vel_samp_comm);
  Barrier(world);
  Bcast(&(vel_samp_class->qseed),1,MPI_DOUBLE,0,world);

/*Randomize random seed */
  qseed = vel_samp_class->qseed;
  for(i=1;i<=myid;i++){
    temp=10000.0*ran_essl(&qseed);
    itemp = temp;
  }/*endfor*/
  if(myid>0){vel_samp_class->qseed = (double)itemp;}

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_energy_ctrl(ENERGY_CTRL *energy_ctrl,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nint = 11;
  MPI_Datatype energy_ctrl_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(energy_ctrl->pme_on),&displs[0]);
  types[0]       = MPI_INT;
  blockcounts[0] = nint;                                    
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&energy_ctrl_comm);
  Barrier(world);
  Type_commit(&energy_ctrl_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,energy_ctrl_comm,0,world);
  Barrier(world);
  Type_free(&energy_ctrl_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_class_parse_info(CLASS_PARSE *class_parse,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  MPI_Datatype class_parse_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int nclass_parse = 10;

  Address(&(class_parse->ivx_smpl),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nclass_parse;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&class_parse_info_comm);
  Barrier(world);
  Type_commit(&class_parse_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,class_parse_info_comm,0,world);
  Barrier(world);
  Type_free(&class_parse_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_stat_avg_info(STAT_AVG *stat_avg,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  Barrier(world);
  Bcast(&(stat_avg->iswit_vdw),1,MPI_INT,0,world);


/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_surface_info(SURFACE *surface,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  MPI_Datatype surf_info_comm;
  MPI_Datatype types[1];
  MPI_Aint     displs[1];
  int          blockcounts[1];
  int          nsurf = 4;

  Address(&(surface->isurf_on),&displs[0]);
  types[0]       = MPI_INT;
  blockcounts[0] = nsurf;
  Barrier(world);

  Type_struct(1,blockcounts,displs,types,&surf_info_comm);
  Barrier(world);

  Type_commit(&surf_info_comm);
  Barrier(world);

  Bcast(MPI_BOTTOM,1,surf_info_comm,0,world);
  Barrier(world);

  Type_free(&surf_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/
