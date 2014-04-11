/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: communicate_output_cp_pimd                   */
/*                                                                          */
/* Routine does the parallel/serial assignment of classical coordinates     */ 
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_coeff_output(CP *cp, double *cre_temp,double *cim_temp,
                              int ip, int pflag,int source, int source_old,
                              COMMUNICATE *communicate,int pi_beads_proc)

/*======================================================================*/
/*                Begin Routine */
  {/*begin routine */
/*======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
   int nread_up,nread_dn,i,ii,tag,ip_loc,iii;
   int rank       = communicate->myid;
   int numprocs   = communicate->np;
   MPI_Comm world = communicate->world;

/*======================================================================*/

 nread_up = (cp->cpcoeffs_info.ncoef)*(cp->cpcoeffs_info.nstate_up);
 nread_dn = (cp->cpcoeffs_info.ncoef)*(cp->cpcoeffs_info.nstate_dn);

if(numprocs>1){

 tag=ip;
 ip_loc = ip - pi_beads_proc*(source_old+1);

 switch(pflag){
  case 1:
   if(rank==source){
     ii = 0;
     for(i=1;i<=nread_up;i++){
       ii++;
       cre_temp[ii] = cp->cpcoeffs_pos[ip_loc].cre_up[i];     
       cim_temp[ii] = cp->cpcoeffs_pos[ip_loc].cim_up[i];
     }/*endfor coeffs*/
   }/* endif : rank=source */
   Barrier(world);
    break;

  case 2:
   if(rank==source){
     ii = 0;
     for(i=1;i<=nread_up;i++){
       ii++;
       cre_temp[ii] = cp->cpcoeffs_pos[ip_loc].cre_dn[i];     
       cim_temp[ii] = cp->cpcoeffs_pos[ip_loc].cim_dn[i];
     }/*endfor coeffs*/
   }/* endif : rank=source */
   Barrier(world);
    break;

  case 3:
   if(rank==source){
     ii = 0;
     for(i=1;i<=nread_up;i++){
       ii++;
       cre_temp[ii] = cp->cpcoeffs_pos[ip_loc].vcre_up[i];     
       cim_temp[ii] = cp->cpcoeffs_pos[ip_loc].vcim_up[i];
     }/*endfor coeffs*/
   }/* endif : rank=source */
   Barrier(world);
    break;

  case 4:
   if(rank==source){
     ii = 0;
     for(i=1;i<=nread_up;i++){
       ii++;
       cre_temp[ii] = cp->cpcoeffs_pos[ip_loc].vcre_dn[i];     
       cim_temp[ii] = cp->cpcoeffs_pos[ip_loc].vcim_dn[i];
     }/*endfor coeffs*/
   }/* endif : rank=source */
   Barrier(world);
    break;
 }/*endif : switch*/

 if(rank==source){
   Send(&(cre_temp[0]),nread_up+1,MPI_DOUBLE,0,tag,world);
   Send(&(cim_temp[0]),nread_up+1,MPI_DOUBLE,0,tag,world);
 }/* endif : rank=source */


 if(rank==0){
    Recv(&(cre_temp[0]),nread_up+1,MPI_DOUBLE,source,tag,world);
     Recv(&(cim_temp[0]),nread_up+1,MPI_DOUBLE,source,tag,world);
 }/* endif : rank=0 */

}/*endif*/

/*==========================================================================*/
} /* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_coeff_nhc_output(CP *cp, double *nhc_temp,int ip,
                                 int source,int source_old,
                                 COMMUNICATE *communicate,int pi_beads_proc)

/*======================================================================*/
/*                Begin Routine                                         */
   {/*begin routine                                                     */
/*======================================================================*/
/*            Local variable declarations                               */
#include "../typ_defs/typ_mask.h"
  
  int ii,i,j,tag,ip_loc;
  int len_nhc    = cp->cptherm_info.len_c_nhc;
  int num_nhc    = cp->cptherm_info.num_c_nhc;
  int rank       = communicate->myid;
  int rank_state = communicate->myid_state;
  int np_states  = communicate->np_states;
  int num_proc   = communicate->np;
  MPI_Comm world = communicate->world;

/*======================================================================*/

Barrier(world);
if(ip<=pi_beads_proc){
  if(rank==0){
    ii = 0;
    for(j=1;j<=len_nhc;j++){
      for(i=1;i<=num_nhc;i++){
        ii++;
        nhc_temp[ii] = cp->cptherm_pos[ip].vc_nhc[j][i];
      }/*endfor*/
    }/*endfor*/
  }/*endif : rank==0*/
}/*endif : ip<=pi_beads_proc*/
Barrier(world);

if(num_proc>1){

 if(ip > pi_beads_proc){
  ip_loc = ip - pi_beads_proc*(source_old+1);
  tag=ip;
  Barrier(world);
  if(rank==source*np_states&&rank_state==0){
    if((len_nhc>0)&&(num_nhc>0)){
      ii = 0;
     printf("commd %d\n",ip);fflush(stdout);
      for(j=1;j<=len_nhc;j++){
        for(i=1;i<=num_nhc;i++){
          ii++;
          nhc_temp[ii] = cp->cptherm_pos[ip_loc].vc_nhc[j][i];
        }/*endfor*/
      }/*endfor*/
      Send(&(nhc_temp[0]),len_nhc*num_nhc+1,MPI_DOUBLE,0,tag,world);
    }/*endif : thermostatting on*/
  }/*endif : rank=source*/
  Barrier(world);

  if(rank==0){
    if((len_nhc>0)&&(num_nhc>0)){
      Recv(&(nhc_temp[0]),len_nhc*num_nhc+1,MPI_DOUBLE,source,tag,world);
    }/*endif : thermostatting on*/
  }/*endif : rank==0*/
  Barrier(world);

 }/*endif : ip>=pi_beads_proc*/

}/*endif*/


/*==========================================================================*/
} /* end routine */
/*==========================================================================*/











