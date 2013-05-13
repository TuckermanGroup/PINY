/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_bond_info.c                             */
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

void communicate_bond_info(BONDED *bonded,NULL_INTER_PARSE *null_inter_parse,
                           MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */

#include "../typ_defs/typ_mask.h"

/*=======================================================================*/
/*             Local variable declarations                                */

  comm_bond_info(&(bonded->bond),world);                      Barrier(world);
  comm_grp_bond_con_info(&(bonded->grp_bond_con),world);      Barrier(world);
  comm_grp_bond_watts_info(&(bonded->grp_bond_watts),world);  Barrier(world);
  comm_bond_free_info(&(bonded->bond_free),world);            Barrier(world);
  comm_bend_info(&(bonded->bend),world);                      Barrier(world);
  comm_bend_free_info(&(bonded->bend_free),world);            Barrier(world);
  comm_bend_bnd_info(&(bonded->bend_bnd),world);              Barrier(world);
  comm_tors_info(&(bonded->tors),world);                      Barrier(world);
  comm_tors_free_info(&(bonded->tors_free),world);            Barrier(world);
  comm_onfo_info(&(bonded->onfo),world);                      Barrier(world);
  comm_ecor_info(&(bonded->ecor),world);                      Barrier(world);
  comm_null_inter_parse_info(null_inter_parse,world);         Barrier(world);
  comm_intra_scr(&(bonded->intra_scr),world);                 Barrier(world);
  comm_constrnt_info(&(bonded->constrnt),world);              Barrier(world);
  comm_rbar_sig_free_info(&(bonded->rbar_sig_free),world);    Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bond_info(BOND *bond,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  MPI_Datatype bond_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];
  int nbond = 10;

  Address(&(bond->npow),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nbond;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&bond_info_comm);
  Barrier(world);
  Type_commit(&bond_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,bond_info_comm,0,world);
  Barrier(world);
  Type_free(&bond_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_con_info(GRP_BOND_CON *grp_bond_con,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int ngrp_bond_con = 21;
  MPI_Datatype grp_bond_con_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(grp_bond_con->num_33),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ngrp_bond_con;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&grp_bond_con_info_comm);
  Barrier(world);
  Type_commit(&grp_bond_con_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,grp_bond_con_info_comm,0,world);
  Barrier(world);
  Type_free(&grp_bond_con_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_grp_bond_watts_info(GRP_BOND_WATTS *grp_bond_watts,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int ngrp_bond_watts = 5;
  MPI_Datatype grp_bond_watts_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(grp_bond_watts->num_33),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ngrp_bond_watts;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&grp_bond_watts_info_comm);
  Barrier(world);
  Type_commit(&grp_bond_watts_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,grp_bond_watts_info_comm,0,world);
  Barrier(world);
  Type_free(&grp_bond_watts_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bond_free_info(BOND_FREE *bond_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nbond_free = 5;
  MPI_Datatype bond_free_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1],iii;

  Address(&(bond_free->num),&displs[0]);
  types[0] = MPI_INT;

  blockcounts[0] = nbond_free;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&bond_free_info_comm);
  Barrier(world);
  Type_commit(&bond_free_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,bond_free_info_comm,0,world);
  Barrier(world);
  Type_free(&bond_free_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_info(BEND *bend,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nbend = 10;
  int iii;
  MPI_Datatype bend_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Barrier(world);
  Address(&(bend->npow),&displs[0]);

  types[0]       = MPI_INT;
  blockcounts[0] = nbend;

  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&bend_info_comm);

  Barrier(world);
  Type_commit(&bend_info_comm);

  Barrier(world);
  Bcast(MPI_BOTTOM,1,bend_info_comm,0,world);

  Barrier(world);
  Type_free(&bend_info_comm);

  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_free_info(BEND_FREE *bend_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nbend_free = 6;
  MPI_Datatype bend_free_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(bend_free->num),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nbend_free;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&bend_free_info_comm);
  Barrier(world);
  Type_commit(&bend_free_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,bend_free_info_comm,0,world);
  Barrier(world);
  Type_free(&bend_free_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/

/* Fixed until here !!! */

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_bend_bnd_info(BEND_BND *bend_bnd,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nbend_bnd = 5;
  MPI_Datatype bend_bnd_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(bend_bnd->num),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nbend_bnd;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&bend_bnd_info_comm);
  Barrier(world);
  Type_commit(&bend_bnd_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,bend_bnd_info_comm,0,world);
  Barrier(world);
  Type_free(&bend_bnd_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_tors_info(TORS *tors,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int ntors = 12;
  MPI_Datatype tors_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(tors->npow),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ntors;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&tors_info_comm);
  Barrier(world);
  Type_commit(&tors_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,tors_info_comm,0,world);
  Barrier(world);
  Type_free(&tors_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_tors_free_info(TORS_FREE *tors_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int ntors_free = 3;
  MPI_Datatype tors_free_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(tors_free->num),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ntors_free;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&tors_free_info_comm);
  Barrier(world);
  Type_commit(&tors_free_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,tors_free_info_comm,0,world);
  Barrier(world);
  Type_free(&tors_free_info_comm);
  Barrier(world);


  Bcast(&(tors_free->j1[1]),2,MPI_INT,0,world);   Barrier(world);
  Bcast(&(tors_free->j2[1]),2,MPI_INT,0,world);   Barrier(world);
  Bcast(&(tors_free->j3[1]),2,MPI_INT,0,world);   Barrier(world);
  Bcast(&(tors_free->j4[1]),2,MPI_INT,0,world);   Barrier(world);


/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_onfo_info(ONFO *onfo,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nonfo = 5;
  MPI_Datatype onfo_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(onfo->num),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nonfo;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&onfo_info_comm);
  Barrier(world);
  Type_commit(&onfo_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,onfo_info_comm,0,world);
  Barrier(world);
  Type_free(&onfo_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_ecor_info(ECOR *ecor,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int necor = 3;
  MPI_Datatype ecor_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(ecor->nsplin),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = necor;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&ecor_info_comm);
  Barrier(world);
  Type_commit(&ecor_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,ecor_info_comm,0,world);
  Barrier(world);
  Type_free(&ecor_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_null_inter_parse_info(NULL_INTER_PARSE *null_inter_parse,
                                   MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int iii;

  int ninfo_null_inter_parse = 4;
  MPI_Datatype null_inter_parse_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(null_inter_parse->nbond_nul),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ninfo_null_inter_parse;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&null_inter_parse_info_comm);
  Barrier(world);
  Type_commit(&null_inter_parse_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,null_inter_parse_info_comm,0,world);
  Barrier(world);
  Type_free(&null_inter_parse_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_intra_scr(INTRA_SCR *intra_scr,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  int nintra_scr = 1;
  MPI_Datatype intra_scr_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(intra_scr->nlen),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = nintra_scr;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&intra_scr_info_comm);
  Barrier(world);
  Type_commit(&intra_scr_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,intra_scr_info_comm,0,world);
  Barrier(world);
  Type_free(&intra_scr_info_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_constrnt_info(CONSTRNT *constrnt,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */


#include "../typ_defs/typ_mask.h"

  MPI_Datatype constrnt_info_comm;
  MPI_Datatype constrnt_data_comm;
  MPI_Datatype types[1];
  MPI_Datatype types1[1];
  MPI_Aint displs[1];
  MPI_Aint displs1[1];
  int blockcounts[1];
  int blockcounts1[1];

  Address(&(constrnt->iconstrnt),&displs[0]);
  Address(&(constrnt->tolshake),&displs1[0]);
  types[0] = MPI_INT;
  types1[0] = MPI_DOUBLE;
  blockcounts[0] = 3;
  blockcounts1[0] = 2;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&constrnt_info_comm);
  Barrier(world);
  Type_struct(1,blockcounts1,displs1,types1,&constrnt_data_comm);
  Barrier(world);
  Type_commit(&constrnt_info_comm);
  Barrier(world);
  Type_commit(&constrnt_data_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,constrnt_info_comm,0,world);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,constrnt_data_comm,0,world);
  Barrier(world);
  Type_free(&constrnt_info_comm);
  Barrier(world);
  Type_free(&constrnt_data_comm);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_rbar_sig_free_info(RBAR_SIG_FREE *rbar_sig_free,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                             */
   {/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                               */


#include "../typ_defs/typ_mask.h"

   int nbar_free = 4;
   MPI_Datatype rbar_free_info_comm;
   MPI_Datatype types[1];
   MPI_Aint displs[1];
   int blockcounts[1];
 
   Address(&(rbar_sig_free->nfree),&displs[0]);
   types[0]       = MPI_INT;
   blockcounts[0] = nbar_free;
   Barrier(world);
   Type_struct(1,blockcounts,displs,types,&rbar_free_info_comm);
   Barrier(world);
   Type_commit(&rbar_free_info_comm);
   Barrier(world);
   Bcast(MPI_BOTTOM,1,rbar_free_info_comm,0,world);
   Barrier(world);
   Type_free(&rbar_free_info_comm);
   Barrier(world);
 
 
/*------------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/

