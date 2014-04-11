/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: analysis_md                                  */
/*                                                                          */
/* This subprogram performs on the fly analysis of MD data                  */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_analysis_md_entry.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define ALL
#define MSQD_ON
#define IIKT_ON
#define ICKT_OFF
#define  DIP_OFF
#define MISC_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void analysis_md(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
		 ANALYSIS *analysis)

/*==========================================================================*/
{ /* begin routine */
 /*=========================================================================*/
 /*              Local variable declarations                                */
 /*=========================================================================*/

#include "../typ_defs/typ_mask.h"

#ifdef ALL
  unsigned int i;

  int num_proc       = class->communicate.np;
  int npforc          = class->class_comm_forc_pkg.num_proc;
  int myid            = class->communicate.myid;
  int myid_forc       = class->class_comm_forc_pkg.myid;
  int myatm_start     = class->clatoms_info.myatm_start;
  int natm_tot        = class->clatoms_info.natm_tot;

  int *recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
  int *displs_atm     = class->class_comm_forc_pkg.displs_atm;

  double *x           = class->clatoms_pos[1].x;
  double *y           = class->clatoms_pos[1].y;
  double *z           = class->clatoms_pos[1].z;
  double *xtemp       = class->ewd_scr.x;
  double *ytemp       = class->ewd_scr.y;
  double *ztemp       = class->ewd_scr.z;
  double *vx          = class->clatoms_pos[1].vx;
  double *vy          = class->clatoms_pos[1].vy;
  double *vz          = class->clatoms_pos[1].vz;
  double *vxtemp      = class->ewd_scr.fx;
  double *vytemp      = class->ewd_scr.fy;
  double *vztemp      = class->ewd_scr.fz;

  MPI_Comm comm_forc  = class->class_comm_forc_pkg.comm;
  MPI_Comm world = class->communicate.world;

 /*=========================================================================*/
 /* I- Gather positions and velocities on processor 0                       */
 /*=========================================================================*/
    if (npforc>1)
    {
     Gatherv(&(x[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(xtemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
     Gatherv(&(y[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(ytemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
     Gatherv(&(z[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(ztemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);

     Gatherv(&(vx[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vxtemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
     Gatherv(&(vy[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vytemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
     Gatherv(&(vz[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vztemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
     if(myid==0)
     {
      for(i=1;i<=natm_tot;i++)
      {
       x[i]  =  xtemp[i];  y[i] =  ytemp[i];  z[i] =  ztemp[i];
       vx[i] = vxtemp[i]; vy[i] = vytemp[i]; vz[i] = vztemp[i];
      }; /* endfor */
     }; /* endif */
    }; /* endif */

 /*=========================================================================*/
 /* II- Do some preliminary stuff on processor 0                            */
 /*=========================================================================*/

  if(myid==0)
  {
   if (general_data->timeinfo.itime==1)
   {
    prelim_analysis(class,general_data,analysis); 
   }; /* endif */ 
  }; /* endif */

 /*=========================================================================*/
 /* III- Calculate the radial distribution functions on processor 0         */
 /*=========================================================================*/

  if(myid==0)
  {
   if (analysis->rdf.calcul_gr_on==1)
   {
    calcul_gr(class,general_data,analysis); 
   }; /* endif */ 
  }; /* endif */

 /*=========================================================================*/
 /* IV- Calculate the velocity correlation functions on processor 0         */
 /*=========================================================================*/

  if(myid==0)
  {
   if (analysis->velocorel.calcul_atm_on==1)
   {
    calcul_vovt_atm(class,general_data,analysis);
   }; /* endif */
  }; /* endif */

 /*=========================================================================*/
 /* V- Calculate the mean square displacement on processor 0                */
 /*=========================================================================*/

#ifdef MSQD_ON
  if(myid==0)
  {
   if (analysis->msqdcorel.calcul_atm_on==1) 
   {
    calcul_msqd_atm(class,general_data,analysis);
   }; /* endif */
  }; /* endif */
#endif

 /*=========================================================================*/
 /* VI- Calculate the incoherent scattering functions on processor 0        */
 /*=========================================================================*/

#ifdef IIKT_ON
  if(myid==0)
  {
   if (analysis->iikt_iso_corel.calcul_on==1) 
   {
    calcul_iikt_iso_cmd(class,general_data,analysis); 
   }; /* endif */
  }; /* endif */
#endif

 /*=========================================================================*/
 /* VII- Calculate the coherent scattering functions on processor 0         */
 /*=========================================================================*/

#ifdef ICKT_ON
  if(myid==0)
  {
   if (analysis->ickt_iso_corel.calcul_on==1) 
   {
    calcul_ickt_iso_md(class,general_data,analysis); 
   }; /* endif */
  }; /* endif */
#endif

 /*=========================================================================*/
 /* VIII- Calculatethe dipole-dipole correlation functions on processor 0   */
 /*=========================================================================*/

#ifdef DIP_ON
  if(myid==0)
  {
   calcul_momt(class,general_data,analysis); */
  }; /* endif */
#endif

 /*=========================================================================*/
 /* IX- Calculate miscellaneous stuff on processor 0                        */
 /*=========================================================================*/

#ifdef MISC_ON
  if(myid==0)
  {
   calcul_miscellaneous(class,general_data,analysis); */
  }; /* endif */
#endif

 /*=========================================================================*/
 /* X- Waiting each others.....                                            */
 /*=========================================================================*/

  if(num_proc>1){  Barrier(world); }

#endif

/*=======================================================================*/
} /* end routine */
/*==========================================================================*/

