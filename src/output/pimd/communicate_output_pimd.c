/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: communicate_output_pimd.c                    */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_friend_lib_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void communicate_output_pimd(CLASS *class)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

#include "../typ_defs/typ_mask.h"

  double *scr_x,*scr_temp_x;
  double *scr_y,*scr_temp_y;
  double *scr_z,*scr_temp_z;
  double *scr_nhc_x,*scr_temp_nhc_x;
  int ip,ipart,ii,num,i,j,k;
  int iproc,iii;
  int len_nhc       = class->therm_info_bead.len_nhc;
  int num_nhc       = class->therm_info_bead.num_nhc;
  int natm_tot      = class->clatoms_info.natm_tot;
  int pi_beads      = class->clatoms_info.pi_beads;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int num_proc      = class->communicate.np_beads;
  int myid          = class->communicate.myid_bead;
  int n2,n3;
  MPI_Comm world    = class->communicate.comm_beads;

/*========================================================================*/
/* I) Malloc the memory                                                   */

  if(myid==0){
    scr_x     = (double *)malloc(pi_beads*natm_tot*sizeof(double)) -1;
    scr_y     = (double *)malloc(pi_beads*natm_tot*sizeof(double)) -1;
    scr_z     = (double *)malloc(pi_beads*natm_tot*sizeof(double)) -1;
    num       = pi_beads*len_nhc*num_nhc;
    scr_nhc_x = (double *)cmalloc(num*sizeof(double)) -1;
  }/*endif*/
  if(myid>0){
    scr_temp_x = (double *)malloc(pi_beads_proc*natm_tot*sizeof(double)) -1;
    scr_temp_y = (double *)malloc(pi_beads_proc*natm_tot*sizeof(double)) -1;
    scr_temp_z = (double *)malloc(pi_beads_proc*natm_tot*sizeof(double)) -1;
    num       = pi_beads_proc*len_nhc*num_nhc;
    scr_temp_nhc_x = (double *)cmalloc(num*sizeof(double)) -1;
  }/*endif*/

/*========================================================================*/
/* III) Send the positions, in order, to processor 0                      */

  if(myid>0){
    ii = 0;
    for(ip=1;ip<=pi_beads_proc;ip++){
      for(ipart=1;ipart<=natm_tot;ipart++){
        ii++;
        scr_temp_x[ii] = class->clatoms_pos[ip].x[ipart]; 
        scr_temp_y[ii] = class->clatoms_pos[ip].y[ipart]; 
        scr_temp_z[ii] = class->clatoms_pos[ip].z[ipart]; 
      }/*endfor : number of atm's */
    }/*endfor : number of beads */
  }/* endif : myid>0 */


  if(myid>0){
    Send(&scr_temp_x[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
    Send(&scr_temp_y[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
    Send(&scr_temp_z[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
  }/*endif*/


  if(myid==0){
    for(iproc=1;iproc<num_proc;iproc++){
        Recv(&(scr_x[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
        Recv(&(scr_y[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
        Recv(&(scr_z[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
    }/*endfor*/
  }/*endif*/    

  if(myid==0){
    ii = 0;
    for(ip=pi_beads_proc+1;ip<=pi_beads;ip++){
      for(ipart=1;ipart<=natm_tot;ipart++){
        ii++;
        class->clatoms_pos[ip].x[ipart] = scr_x[ii];
        class->clatoms_pos[ip].y[ipart] = scr_y[ii];
        class->clatoms_pos[ip].z[ipart] = scr_z[ii];
      }/*endfor : number of atm's */
    }/*endfor : number of beads */
  }/* endif : myid=0 */

/*========================================================================*/
/* IV ) Send the velocities, in order, to processor 0                     */

  ii = 0;
  if(myid>0){
    for(ip=1;ip<=pi_beads_proc;ip++){
      for(ipart=1;ipart<=natm_tot;ipart++){
        ii++;
        scr_temp_x[ii] = class->clatoms_pos[ip].vx[ipart]; 
        scr_temp_y[ii] = class->clatoms_pos[ip].vy[ipart]; 
        scr_temp_z[ii] = class->clatoms_pos[ip].vz[ipart]; 
      }/*endfor : number of atm's */
    }/*endfor : number of beads */
  }/* endif : myid>0 */


  if(myid>0){
    Send(&scr_temp_x[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
    Send(&scr_temp_y[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
    Send(&scr_temp_z[1],(pi_beads_proc)*(natm_tot),
                      MPI_DOUBLE,0,myid,world);
  }/*endif*/


  if(myid==0){
    for(iproc=1;iproc<num_proc;iproc++){
        Recv(&(scr_x[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
        Recv(&(scr_y[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
        Recv(&(scr_z[(iproc-1)*pi_beads_proc*natm_tot+1]),
           (pi_beads_proc)*(natm_tot),MPI_DOUBLE,
                    iproc,iproc,world);
    }/*endfor*/
  }/*endif*/    
  if(myid==0){
    ii = 0;
    for(ip=pi_beads_proc+1;ip<=pi_beads;ip++){
      for(ipart=1;ipart<=natm_tot;ipart++){
        ii++;
        class->clatoms_pos[ip].vx[ipart] = scr_x[ii];
        class->clatoms_pos[ip].vy[ipart] = scr_y[ii];
        class->clatoms_pos[ip].vz[ipart] = scr_z[ii];
      }/*endfor : number of atm's */
    }/*endfor : number of beads */
  }/* endif : myid=0 */

/*========================================================================*/
/* IV ) Send the NHC velocities, in order, to processor 0                 */

  n2 = len_nhc;
  n3 = num_nhc;
  if(myid>0){
    for(k=1;k<= pi_beads_proc; k++){
       for(j=1;j<=(len_nhc);j++){
          for(i=1;i<=(num_nhc);i++){
             ii = (i-1) + (j-1)*n3 + (k-1)*n2*n3 + 1;
              scr_temp_nhc_x[ii] = class->therm_bead[k].v_nhc[j][i];
          }/*endfor : nhc length */
       }/*endfor : number of nhc's */
    }/*endfor : number of beads */
  }/* endif : myid>0 */

  if(myid>0){
    Send(&scr_temp_nhc_x[1],pi_beads_proc*len_nhc*num_nhc,
                      MPI_DOUBLE,0,myid,world);
  }/*endif*/


  if(myid==0){
    for(iproc=1;iproc<num_proc;iproc++){
        Recv(&(scr_nhc_x[(iproc-1)*pi_beads_proc*len_nhc*num_nhc+1]),
          pi_beads_proc*len_nhc*num_nhc,MPI_DOUBLE,
                    iproc,iproc,world);
    }/*endfor*/
  }/*endif*/    

  if(myid==0){
    for(k=pi_beads_proc+1;k<= pi_beads; k++){
       for(j=1;j<=(len_nhc);j++){
          for(i=1;i<=(num_nhc);i++){
            ii = (i-1) + (j-1)*n3 + (k-pi_beads_proc-1)*n2*n3 + 1;
            class->therm_bead[k].v_nhc[j][i] = scr_nhc_x[ii];
         }/*endfor : nhc length */
      }/*endfor : number of nhc's */
    }/*endfor : number of beads */
  }/*endif*/

 if(myid==0){
     cfree(&(scr_x[1]));
     cfree(&(scr_y[1]));
     cfree(&(scr_z[1]));
     cfree(&(scr_nhc_x[1]));
 }else{
     cfree(&(scr_temp_x[1]));
     cfree(&(scr_temp_y[1]));
     cfree(&(scr_temp_z[1]));
     cfree(&(scr_temp_nhc_x[1]));
 }/*endif*/

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/









