/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_cp_info.c                               */
/*                                                                          */
/* Subprogram contains MPI utils and communication routines for interface   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_math.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_cp_info(CP *cp,CP_PARSE *cp_parse,MPI_Comm world,int myid)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */

#include "../typ_defs/typ_mask.h"

/*=======================================================================*/
/*             Local variable declarations                                */

  comm_cpcoeffs_info(&(cp->cpcoeffs_info),world);  Barrier(world);
  comm_cptherm_info(&(cp->cptherm_info),world);    Barrier(world);
  comm_cp_parse_info(cp_parse,world);              Barrier(world);
  comm_vel_samp_cp(&(cp->vel_samp_cp),world,myid); Barrier(world);

  Bcast(&(cp->cp_sclr_fft_pkg3d_lg.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(cp->cp_para_fft_pkg3d_lg.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(cp->cp_sclr_fft_pkg3d_sm.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(cp->cp_para_fft_pkg3d_sm.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(cp->cp_sclr_fft_pkg3d_dens_cp_box.igeneric_opt),1,MPI_INT,0,world);
  Bcast(&(cp->cp_para_fft_pkg3d_dens_cp_box.igeneric_opt),1,MPI_INT,0,world);

  Bcast(&(cp->cpscr.cpscr_atom_pme.n_interp),1,MPI_INT,0,world);
  Bcast(&(cp->cpscr.cpscr_atom_pme.pme_on),1,MPI_INT,0,world);
  Bcast(&(cp->cpscr.cpscr_atom_pme.nlen_pme),1,MPI_INT,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cpcoeffs_info(CPCOEFFS_INFO *cpcoeffs_info,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  int ninfo = 14;
  MPI_Datatype cpcoeffs_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(cpcoeffs_info->pi_beads),&displs[0]);
  types[0] = MPI_INT;
  blockcounts[0] = ninfo;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cpcoeffs_info_comm);
  Barrier(world);
  Type_commit(&cpcoeffs_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cpcoeffs_info_comm,0,world);
  Barrier(world);
  Type_free(&cpcoeffs_info_comm);
  Barrier(world);
  Bcast(&(cpcoeffs_info->ecut),1,MPI_DOUBLE,0,world);
  Bcast(&(cpcoeffs_info->cp_hess_cut),1,MPI_DOUBLE,0,world);

/*------------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cptherm_info(CPTHERM_INFO *cptherm_info,MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  
#include "../typ_defs/typ_mask.h"

  int ninfo = 5;
  MPI_Datatype cptherm_info_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(cptherm_info->len_c_nhc),&displs[0]);
  types[0]       = MPI_INT;
  blockcounts[0] = ninfo;
  Barrier(world);
  Type_struct(1,blockcounts,displs,types,&cptherm_info_comm);
  Barrier(world);
  Type_commit(&cptherm_info_comm);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cptherm_info_comm,0,world);
  Barrier(world);
  Type_free(&cptherm_info_comm);
  Barrier(world);

  Bcast(&(cptherm_info->cp_therm_heat_fact),1,MPI_DOUBLE,0,world);
  Barrier(world);

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_cp_parse_info(CP_PARSE *cp_parse, MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  
#include "../typ_defs/typ_mask.h"

  int ncp_parse_int    = 6;
  int ncp_parse_double = 10;
  MPI_Datatype cp_parse_info_comm_int;
  MPI_Datatype cp_parse_info_comm_double;
  MPI_Datatype types_int[1];
  MPI_Datatype types_double[1];
  MPI_Aint displs_int[1];
  MPI_Aint displs_double[1];
  int blockcounts_int[1];
  int blockcounts_double[1];

  Address(&(cp_parse->istart_cp),&displs_int[0]);
  Address(&(cp_parse->cp_mass_tau_def),&displs_double[0]);
  types_int[0] = MPI_INT;
  types_double[0] = MPI_DOUBLE;
  blockcounts_int[0] = ncp_parse_int;
  blockcounts_double[0] = ncp_parse_double;
  Barrier(world);
  Type_struct(1,blockcounts_int,displs_int,types_int,
                &cp_parse_info_comm_int);
  Barrier(world);
  Type_struct(1,blockcounts_double,displs_double,types_double,
                &cp_parse_info_comm_double);
  Barrier(world);
  Type_commit(&cp_parse_info_comm_int);
  Barrier(world);
  Type_commit(&cp_parse_info_comm_double);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cp_parse_info_comm_int,0,world);
  Barrier(world);
  Bcast(MPI_BOTTOM,1,cp_parse_info_comm_double,0,world);
  Barrier(world);
  Type_free(&cp_parse_info_comm_int);
  Barrier(world);
  Type_free(&cp_parse_info_comm_double);
  Barrier(world);


/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_vel_samp_cp(VEL_SAMP_CP *vel_samp_cp,MPI_Comm world,int myid)

/*=======================================================================*/
/*             Begin routine                                             */
{           /* begin routine */
/*=======================================================================*/
/*             Local variable declarations                               */

#include "../typ_defs/typ_mask.h"

  int i,itemp;
  double temp,qseed;
  int nint = 8;
  MPI_Datatype vel_samp_comm;
  MPI_Datatype types[1];
  MPI_Aint displs[1];
  int blockcounts[1];

  Address(&(vel_samp_cp->ivelc_smpl_on),&displs[0]);
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
  Bcast(&(vel_samp_cp->qseed),1,MPI_DOUBLE,0,world);
  Bcast(&(vel_samp_cp->vc_scal_tol),1,MPI_DOUBLE,0,world);

/*Randomize random seed */
  qseed = vel_samp_cp->qseed;
  for(i=1;i<=myid;i++){
    temp=10000.0*ran_essl(&qseed);
    itemp = temp;
  }/*endfor*/
  if(myid>0){vel_samp_cp->qseed = (double)itemp;}

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/












