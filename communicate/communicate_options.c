/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_options.c                           */
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

void communicate_general_data(GENERAL_DATA *general_data,CP *cp, 
                              MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  comm_simopts(&(general_data->simopts),world);
  comm_ensopts(&(general_data->ensopts),world);
  comm_minopts(&(general_data->minopts),world);
  comm_cpopts(&(cp->cpopts),world);
  comm_filenames(&(general_data->filenames),world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_simopts(SIMOPTS *simopts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  MPI_Datatype simopts_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(simopts->md),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = 18;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&simopts_comm);
  Barrier(world);
  Type_commit(&simopts_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,simopts_comm,0,world);
  Barrier(world);
  Type_free(&simopts_comm);
  Barrier(world);
  Bcast(&(simopts->ann_rate),1,MPI_DOUBLE,0,world);
  Bcast(&(simopts->ann_start_temp),1,MPI_DOUBLE,0,world);
  Bcast(&(simopts->ann_target_temp),1,MPI_DOUBLE,0,world);



/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_minopts(MINOPTS *minopts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int nmin_int    = 10;
  MPI_Datatype minopts_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(minopts->min_std),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nmin_int;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&minopts_comm);
  Barrier(world);
  Type_commit(&minopts_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,minopts_comm,0,world);
  Barrier(world);
  Type_free(&minopts_comm);
  Barrier(world);
  Bcast(&(minopts->tol_coef),1,MPI_DOUBLE,0,world);
  Barrier(world);
  Bcast(&(minopts->tol_atom),1,MPI_DOUBLE,0,world);



/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ensopts(ENSOPTS *ensopts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  int nensopts = 5;
  MPI_Datatype ensopts_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(ensopts->nve),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nensopts;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&ensopts_comm);
  Barrier(world);
  Type_commit(&ensopts_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,ensopts_comm,0,world);
  Barrier(world);
  Type_free(&ensopts_comm);



/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cpopts(CPOPTS *cpopts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  int cpopts_num = 39;
  MPI_Datatype cpopts_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(cpopts->cp_lda),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = cpopts_num;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cpopts_comm);
  Barrier(world);
  Type_commit(&cpopts_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cpopts_comm,0,world);
  Barrier(world);
  Type_free(&cpopts_comm);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_filenames(FILENAMES *filenames,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"
  int filenames_num = 14;
  MPI_Datatype filenames_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(filenames->iwrite_screen),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = filenames_num;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&filenames_comm);
  Barrier(world);
  Type_commit(&filenames_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,filenames_comm,0,world);
  Barrier(world);
  Type_free(&filenames_comm);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/












