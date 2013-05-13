/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_pimd_trans_force                     */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_pimd_trans_force(CLASS *class,GENERAL_DATA *general_data)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

 if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==1)){
   convert_pimd_force_stag(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
 }/*endif*/


 if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==1)){
   convert_pimd_force_stag_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
 }/*endif*/



 if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==2)){
   convert_pimd_force_cent(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
 }/*endif*/



 if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==2)){
   convert_pimd_force_cent_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
 }/*endif*/


/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_pimd_trans_mode(CLASS *class,GENERAL_DATA *general_data)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */


   int iii;

   if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==1)){
     convert_pimd_mode_stag(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
   }/*endif*/


   if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==1)){
     convert_pimd_mode_stag_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
   }/*endif*/



   if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==2)){
     convert_pimd_mode_cent(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
   }/*endif*/



   if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==2)){
     convert_pimd_mode_cent_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
   }/*endif*/


/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_pimd_trans_pos(CLASS *class,GENERAL_DATA *general_data)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

   int i,iii;
   double *xmod = class->clatoms_info.xmod,sum=0.0;
   double *ymod = class->clatoms_info.ymod;
   double *zmod = class->clatoms_info.zmod;
   double *x    = class->clatoms_pos[1].x;
   double *y    = class->clatoms_pos[1].y;
   double *z    = class->clatoms_pos[1].z;
   double *xtemp = class->ewd_scr.x;
   double *ytemp = class->ewd_scr.y;
   double *ztemp = class->ewd_scr.z;
   int natm_tot = class->clatoms_info.natm_tot;
   int myid = class->communicate.myid_bead;
   int myid_forc = class->communicate.myid_forc;
   int numproc = class->communicate.np_beads;
   int np_forc = class->communicate.np_forc;
   int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
   int pi_beads         = class->clatoms_info.pi_beads;
   int myatm_start = class->clatoms_info.myatm_start;
   int myatm_end = class->clatoms_info.myatm_end;
   int *recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
   int *displs_atm = class->class_comm_forc_pkg.displs_atm;
   MPI_Comm comm_beads = class->communicate.comm_beads;
   MPI_Comm comm_beads_forc = class->communicate.comm_beads_forc;
   MPI_Comm comm_forc = class->communicate.comm_forc;

   if(pi_beads_proc_st==1){
     for(i=myatm_start;i<=myatm_end;i++){
       xmod[i] = x[i];
       ymod[i] = y[i];
       zmod[i] = z[i];
     }/*endfor*/
   }/*endif*/

   if(numproc>1){
    Bcast(&(xmod[0]),natm_tot+1,MPI_DOUBLE,0,comm_beads_forc);
    Bcast(&(ymod[0]),natm_tot+1,MPI_DOUBLE,0,comm_beads_forc);
    Bcast(&(zmod[0]),natm_tot+1,MPI_DOUBLE,0,comm_beads_forc);
   }/*endif*/

   if(np_forc>1){
     Allgatherv(&(xmod[myatm_start]),recv_count_atm[(myid_forc+1)],
                 MPI_DOUBLE,&(xtemp[1]),&recv_count_atm[1],&displs_atm[1],
                 MPI_DOUBLE,0,comm_forc);
     Allgatherv(&(ymod[myatm_start]),recv_count_atm[(myid_forc+1)],
                 MPI_DOUBLE,&(ytemp[1]),&recv_count_atm[1],&displs_atm[1],
                 MPI_DOUBLE,0,comm_forc);
     Allgatherv(&(zmod[myatm_start]),recv_count_atm[(myid_forc+1)],
                 MPI_DOUBLE,&(ztemp[1]),&recv_count_atm[1],&displs_atm[1],
                 MPI_DOUBLE,0,comm_forc);

     for(i=1;i<=natm_tot;i++){
       xmod[i] = xtemp[i];
       ymod[i] = ytemp[i];
       zmod[i] = ztemp[i];
     }/*endfor*/
   }/*endif*/

   if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==1)){
     convert_pimd_pos_stag(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
   }/*endif*/


   if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==1)){
     convert_pimd_pos_stag_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
   }/*endif*/



   if((class->communicate.np_beads==1)&&(general_data->simopts.pi_md_typ==2)){
     convert_pimd_pos_cent(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran));
   }/*endif*/


   if((class->communicate.np_beads>1)&&(general_data->simopts.pi_md_typ==2)){
     convert_pimd_pos_cent_par(&(class->clatoms_info),class->clatoms_pos,
              &(class->clatoms_tran),&(class->communicate));
   }/*endif*/

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/






