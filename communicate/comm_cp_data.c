/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_cp_data.c                               */
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
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_cp_data(CP *cp,int natm_mall,int natm_typ_mall,
                                     MPI_Comm world,int myid)

/*=======================================================================*/
/*             Begin routine                                              */
   {/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   comm_cpconstrnt(&(cp->cpconstrnt),world);   Barrier(world);
   comm_pseudo(&(cp->pseudo),world,myid);      Barrier(world);
   comm_cpopts_data(&(cp->cpopts),world);      Barrier(world);

/*------------------------------------------------------------------------*/
    } /*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cpconstrnt(CPCONSTRNT *cpconstrnt,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nscal_cpconstrnt = 3;
  MPI_Datatype cpconstrnt_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(cpconstrnt->c_tolshake),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = nscal_cpconstrnt;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cpconstrnt_comm);
  Barrier(world);
  Type_commit(&cpconstrnt_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cpconstrnt_comm,0,world);
  Barrier(world);
  Type_free(&cpconstrnt_comm);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_pseudo(PSEUDO *pseudo,MPI_Comm world,int myid)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"


    if(myid!=0){
      pseudo->vxc_typ  = (char *)cmalloc(MAXWORD*sizeof(char));   
      pseudo->ggax_typ = (char *)cmalloc(MAXWORD*sizeof(char));   
      pseudo->ggac_typ = (char *)cmalloc(MAXWORD*sizeof(char));   
    }/*endif*/
    Barrier(world);
    Bcast(&(pseudo->vxc_typ[0]),MAXWORD,MPI_CHAR,0,world); /* must be from 0*/
    Bcast(&(pseudo->ggax_typ[0]),MAXWORD,MPI_CHAR,0,world);
    Bcast(&(pseudo->ggac_typ[0]),MAXWORD,MPI_CHAR,0,world);
    Bcast(&(pseudo->gga_cut),1,MPI_DOUBLE,0,world);
    Bcast(&(pseudo->alpha_conv_dual),1,MPI_DOUBLE,0,world);
    Bcast(&(pseudo->n_interp_pme_dual),1,MPI_INT,0,world);
    Bcast(&(pseudo->nsplin_g),1,MPI_INT,0,world);
    Bcast(&(pseudo->nl_cut_on),1,MPI_INT,0,world);
    Bcast(&(pseudo->nlvps_skin),1,MPI_DOUBLE,0,world);

    Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cpopts_data(CPOPTS *cpopts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int nscal_cpopts = 3;
  MPI_Datatype cpopts_data_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int count;

  Address(&(cpopts->te_ext),&displs[0]);
  types[0] = MPI_DOUBLE;
  blockcounts[0] = nscal_cpopts;
  
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cpopts_data_comm);
  Barrier(world);
  Type_commit(&cpopts_data_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cpopts_data_comm,0,world);
  Barrier(world);
  Type_free(&cpopts_data_comm);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/











