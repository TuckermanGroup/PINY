/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: transform                                    */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pimd_trans_comm_fwd(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CLATOMS_TRAN *clatoms_tran,COMMUNICATE *communicate,
                      int iopt)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"
  
  int iproc,itemp,iii;
  int iatm,ip,ioff,ioff_t;
  int i,ioff_temp;
  int pi_atm_now;
/*    Local Pointers  */
  double *clatoms_pos_x;
  double *clatoms_pos_y;
  double *clatoms_pos_z;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  double *x_temp_pt,*y_temp_pt,*z_temp_pt;
  double *xt_temp = clatoms_tran->xt_temp;
  double *yt_temp = clatoms_tran->yt_temp;
  double *zt_temp = clatoms_tran->zt_temp;
  double *xt_temp_pt,*yt_temp_pt,*zt_temp_pt;
  int natm_tot      = clatoms_info->natm_tot;
  int pi_beads      = clatoms_info->pi_beads;
  int pi_atm_proc   = clatoms_tran->pi_atm_proc;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_proc_rem   = clatoms_tran->pi_proc_rem;
  int *sendcounts   = clatoms_tran->sendcounts;
  int *senddspls    = clatoms_tran->senddspls;
  int *recvcounts   = clatoms_tran->recvcounts;
  int *recvdspls    = clatoms_tran->recvdspls;
  int num_proc      = communicate->np_beads;
  int myid          = communicate->myid_bead;
  int myatm_start   = clatoms_info->myatm_start;
  int myatm_end     = clatoms_info->myatm_end;
  MPI_Comm comm_beads = communicate->comm_beads;
  MPI_Comm comm_beads_forc = communicate->comm_beads_forc;

/*========================================================================*/
/* I) Extract the position data from clatoms_pos                        */

  for(i=1;i<=pi_atm_proc*pi_beads;i++){x_temp[i]=0.0;}
  for(i=1;i<=pi_atm_proc*pi_beads;i++){y_temp[i]=0.0;}
  for(i=1;i<=pi_atm_proc*pi_beads;i++){z_temp[i]=0.0;}
  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = myatm_start-1;
    switch(iopt) {
      case 1:
        clatoms_pos_x = clatoms_pos[ip].x;
        clatoms_pos_y = clatoms_pos[ip].y;
        clatoms_pos_z = clatoms_pos[ip].z;break;
      case 2:
        clatoms_pos_x = clatoms_pos[ip].fx;
        clatoms_pos_y = clatoms_pos[ip].fy;
        clatoms_pos_z = clatoms_pos[ip].fz;break;
    }/*end switch*/
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>pi_proc_rem)&&(pi_proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        x_temp[itemp] = clatoms_pos_x[i];
        y_temp[itemp] = clatoms_pos_y[i];
        z_temp[itemp] = clatoms_pos_z[i];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   
 
#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
      printf("%d %d %g\n",iproc,iatm,x_temp[iatm]);
      printf("%d %d %g\n",iproc,iatm,y_temp[iatm]);
      printf("%d %d %g\n",iproc,iatm,z_temp[iatm]);
    }/*endfor*/
    }/*endif*/
    Dbx_Barrier(comm_beads_forc);
    scanf("%d",&iii);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Send the position data                                             */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  yt_temp_pt = yt_temp+1;
  zt_temp_pt = zt_temp+1;
  x_temp_pt = x_temp+1;
  y_temp_pt = y_temp+1;
  z_temp_pt = z_temp+1;
  Alltoallv(x_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,xt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(y_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,yt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(z_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,zt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,yt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,zt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
/* III) Rearrange the position data                                         */

  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      x_temp[ioff] = xt_temp[i];
      y_temp[ioff] = yt_temp[i];
      z_temp[ioff] = zt_temp[i];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,y_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,z_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pimd_trans_comm_bck(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CLATOMS_TRAN *clatoms_tran,COMMUNICATE *communicate,
                      int iopt)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"
  
  int iproc,itemp,iii;
  int iatm,ip,ioff,ioff_t;
  int i,ioff_temp;
  int pi_atm_now;
/*    Local Pointers  */
  double *clatoms_pos_x;
  double *clatoms_pos_y;
  double *clatoms_pos_z;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  double *x_temp_pt,*y_temp_pt,*z_temp_pt;
  double *xt_temp = clatoms_tran->xt_temp;
  double *yt_temp = clatoms_tran->yt_temp;
  double *zt_temp = clatoms_tran->zt_temp;
  double *xt_temp_pt,*yt_temp_pt,*zt_temp_pt;
  int natm_tot      = clatoms_info->natm_tot;
  int pi_beads      = clatoms_info->pi_beads;
  int pi_atm_proc   = clatoms_tran->pi_atm_proc;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_proc_rem   = clatoms_tran->pi_proc_rem;
  int *sendcounts   = clatoms_tran->sendcounts;
  int *senddspls    = clatoms_tran->senddspls;
  int *recvcounts   = clatoms_tran->recvcounts;
  int *recvdspls    = clatoms_tran->recvdspls;
  int num_proc      = communicate->np_beads;
  int myid          = communicate->myid_bead;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end   = clatoms_info->myatm_end;
  MPI_Comm comm_beads = communicate->comm_beads;
  MPI_Comm comm_beads_forc = communicate->comm_beads_forc;

/*========================================================================*/
/* I) Rearrange the transformed position data                             */


  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      xt_temp[i] = x_temp[ioff];
      yt_temp[i] = y_temp[ioff];
      zt_temp[i] = z_temp[ioff];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,yt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,zt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Send the transformed position data                               */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  yt_temp_pt = yt_temp+1;
  zt_temp_pt = zt_temp+1;
  x_temp_pt = x_temp+1;
  y_temp_pt = y_temp+1;
  z_temp_pt = z_temp+1;

  Alltoallv(xt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,x_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(yt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,y_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(zt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,z_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,y_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,z_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif


/*========================================================================*/
/* III) Extract the transformed position data                               */

  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = myatm_start-1;
    switch(iopt) {
      case 1:
        clatoms_pos_x = clatoms_pos[ip].x;
        clatoms_pos_y = clatoms_pos[ip].y;
        clatoms_pos_z = clatoms_pos[ip].z;break;
      case 2:
        clatoms_pos_x = clatoms_pos[ip].fx;
        clatoms_pos_y = clatoms_pos[ip].fy;
        clatoms_pos_z = clatoms_pos[ip].fz;break;
    }/*end switch*/
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>pi_proc_rem)&&(pi_proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        clatoms_pos_x[i] = x_temp[itemp];
        clatoms_pos_y[i] = y_temp[itemp];
        clatoms_pos_z[i] = z_temp[itemp];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is clatoms_pos_x\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ip=1;ip<=pi_beads_proc;ip++){
       clatoms_pos_x = clatoms_pos[ip].x;
       clatoms_pos_y = clatoms_pos[ip].y;
       clatoms_pos_z = clatoms_pos[ip].z;
       for(iatm=1;iatm<=natm_tot;iatm++){
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_x[iatm]);
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_y[iatm]);
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_z[iatm]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pimd_trans_comm_fwd_full(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CLATOMS_TRAN *clatoms_tran,COMMUNICATE *communicate,
                      int iopt)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"
  
  int iproc,itemp,iii;
  int iatm,ip,ioff,ioff_t;
  int i,ioff_temp;
  int pi_atm_now;
/*    Local Pointers  */
  double *clatoms_pos_x;
  double *clatoms_pos_y;
  double *clatoms_pos_z;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  double *x_temp_pt,*y_temp_pt,*z_temp_pt;
  double *xt_temp = clatoms_tran->xt_temp;
  double *yt_temp = clatoms_tran->yt_temp;
  double *zt_temp = clatoms_tran->zt_temp;
  double *xt_temp_pt,*yt_temp_pt,*zt_temp_pt;
  int natm_tot      = clatoms_info->natm_tot;
  int pi_beads      = clatoms_info->pi_beads;
  int pi_atm_proc   = clatoms_tran->pi_atm_proc;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_proc_rem   = clatoms_tran->pi_proc_rem;
  int *sendcounts   = clatoms_tran->sendcounts;
  int *senddspls    = clatoms_tran->senddspls;
  int *recvcounts   = clatoms_tran->recvcounts;
  int *recvdspls    = clatoms_tran->recvdspls;
  int num_proc      = communicate->np_beads;
  int myid          = communicate->myid_bead;
  int myatm_start   = clatoms_info->myatm_start;
  int myatm_end     = clatoms_info->myatm_end;
  MPI_Comm comm_beads = communicate->comm_beads;
  MPI_Comm comm_beads_forc = communicate->comm_beads_forc;

/*========================================================================*/
/* I) Extract the position data from clatoms_pos                        */


  for(i=1;i<=pi_atm_proc*pi_beads;i++){x_temp[i]=0.0;}
  for(i=1;i<=pi_atm_proc*pi_beads;i++){y_temp[i]=0.0;}
  for(i=1;i<=pi_atm_proc*pi_beads;i++){z_temp[i]=0.0;}
  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = 0;
    switch(iopt) {
      case 1:
        clatoms_pos_x = clatoms_pos[ip].x;
        clatoms_pos_y = clatoms_pos[ip].y;
        clatoms_pos_z = clatoms_pos[ip].z;break;
      case 2:
        clatoms_pos_x = clatoms_pos[ip].fx;
        clatoms_pos_y = clatoms_pos[ip].fy;
        clatoms_pos_z = clatoms_pos[ip].fz;break;
    }/*end switch*/
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>pi_proc_rem)&&(pi_proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        x_temp[itemp] = clatoms_pos_x[i];
        y_temp[itemp] = clatoms_pos_y[i];
        z_temp[itemp] = clatoms_pos_z[i];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   
 
#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
      printf("%d %d %g\n",iproc,iatm,x_temp[iatm]);
      printf("%d %d %g\n",iproc,iatm,y_temp[iatm]);
      printf("%d %d %g\n",iproc,iatm,z_temp[iatm]);
    }/*endfor*/
    }/*endif*/
    Dbx_Barrier(comm_beads_forc);
    scanf("%d",&iii);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Send the position data                                             */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  yt_temp_pt = yt_temp+1;
  zt_temp_pt = zt_temp+1;
  x_temp_pt = x_temp+1;
  y_temp_pt = y_temp+1;
  z_temp_pt = z_temp+1;
  Alltoallv(x_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,xt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(y_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,yt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(z_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,zt_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,yt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,zt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
/* III) Rearrange the position data                                         */

  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      x_temp[ioff] = xt_temp[i];
      y_temp[ioff] = yt_temp[i];
      z_temp[ioff] = zt_temp[i];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,y_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,z_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void pimd_trans_comm_bck_full(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                      CLATOMS_TRAN *clatoms_tran,COMMUNICATE *communicate,
                      int iopt)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"
  
  int iproc,itemp,iii;
  int iatm,ip,ioff,ioff_t;
  int i,ioff_temp;
  int pi_atm_now;
/*    Local Pointers  */
  double *clatoms_pos_x;
  double *clatoms_pos_y;
  double *clatoms_pos_z;
  double *x_temp = clatoms_tran->x_temp;
  double *y_temp = clatoms_tran->y_temp;
  double *z_temp = clatoms_tran->z_temp;
  double *x_temp_pt,*y_temp_pt,*z_temp_pt;
  double *xt_temp = clatoms_tran->xt_temp;
  double *yt_temp = clatoms_tran->yt_temp;
  double *zt_temp = clatoms_tran->zt_temp;
  double *xt_temp_pt,*yt_temp_pt,*zt_temp_pt;
  int natm_tot      = clatoms_info->natm_tot;
  int pi_beads      = clatoms_info->pi_beads;
  int pi_atm_proc   = clatoms_tran->pi_atm_proc;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_proc_rem   = clatoms_tran->pi_proc_rem;
  int *sendcounts   = clatoms_tran->sendcounts;
  int *senddspls    = clatoms_tran->senddspls;
  int *recvcounts   = clatoms_tran->recvcounts;
  int *recvdspls    = clatoms_tran->recvdspls;
  int num_proc      = communicate->np_beads;
  int myid          = communicate->myid_bead;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end   = clatoms_info->myatm_end;
  MPI_Comm comm_beads = communicate->comm_beads;
  MPI_Comm comm_beads_forc = communicate->comm_beads_forc;

/*========================================================================*/
/* I) Rearrange the transformed position data                             */


  for(ip=1;ip<=pi_beads;ip++){
    ioff_t = (ip-1)*pi_atm_proc;
    for(iatm=1;iatm<=pi_atm_proc;iatm++){
      ioff = (iatm-1)*pi_beads+ip;
      i    = iatm+ioff_t;
      xt_temp[i] = x_temp[ioff];
      yt_temp[i] = y_temp[ioff];
      zt_temp[i] = z_temp[ioff];
    }/*endfor*/
  }/*endfor*/

#ifdef DEBUG
  if(myid==0){
    printf("This is xt_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,xt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,yt_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,zt_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
/* II) Send the transformed position data                               */

  for(ip=1;ip<=num_proc;ip++){
    sendcounts[ip] = pi_atm_proc*pi_beads_proc;
    senddspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
    recvcounts[ip] = pi_atm_proc*pi_beads_proc;
    recvdspls[ip]  = (ip-1)*pi_atm_proc*pi_beads_proc;
  }/*endfor*/

  xt_temp_pt = xt_temp+1;
  yt_temp_pt = yt_temp+1;
  zt_temp_pt = zt_temp+1;
  x_temp_pt = x_temp+1;
  y_temp_pt = y_temp+1;
  z_temp_pt = z_temp+1;

  Alltoallv(xt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,x_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(yt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,y_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);
  Alltoallv(zt_temp_pt,&(sendcounts[1]),&(senddspls[1]),
                MPI_DOUBLE,z_temp_pt,&(recvcounts[1]),
                &(recvdspls[1]), MPI_DOUBLE,comm_beads_forc);

#ifdef DEBUG
  if(myid==0){
    printf("This is x_temp\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(ip=1;ip<=num_proc;ip++){
    if(myid==ip-1){
      for(iatm=1;iatm<=pi_atm_proc*pi_beads;iatm++){
         printf("%d %d %g\n",iatm,ip,x_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,y_temp[iatm]);
         printf("%d %d %g\n",iatm,ip,z_temp[iatm]);
      }/*endfor*/
    }/*endif*/
    scanf("%d",&iatm);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif


/*========================================================================*/
/* III) Extract the transformed position data                               */

  for(ip=1;ip<=pi_beads_proc;ip++){
    ioff = myatm_start-1;
    switch(iopt) {
      case 1:
        clatoms_pos_x = clatoms_pos[ip].x;
        clatoms_pos_y = clatoms_pos[ip].y;
        clatoms_pos_z = clatoms_pos[ip].z;break;
      case 2:
        clatoms_pos_x = clatoms_pos[ip].fx;
        clatoms_pos_y = clatoms_pos[ip].fy;
        clatoms_pos_z = clatoms_pos[ip].fz;break;
    }/*end switch*/
    for(iproc=1;iproc<=num_proc;iproc++){
      ioff_temp = (ip-1)*pi_atm_proc+(iproc-1)*(pi_beads_proc*pi_atm_proc);
      pi_atm_now = pi_atm_proc;
      if((iproc>pi_proc_rem)&&(pi_proc_rem>0)){pi_atm_now = pi_atm_proc-1;}
      for(iatm=1;iatm<=pi_atm_now;iatm++){
        itemp = iatm+ioff_temp;
        i     = iatm+ioff;
        clatoms_pos_x[i] = x_temp[itemp];
        clatoms_pos_y[i] = y_temp[itemp];
        clatoms_pos_z[i] = z_temp[itemp];
      }/*endfor*/   
      ioff += pi_atm_now;
    }/*endfor*/   
  }/*endfor*/   

#ifdef DEBUG
  if(myid==0){
    printf("This is clatoms_pos_x\n");
  }/*endif*/
  Dbx_Barrier(comm_beads_forc);
  for(iproc=1;iproc<=num_proc;iproc++){
    if(myid==iproc-1){
     for(ip=1;ip<=pi_beads_proc;ip++){
       clatoms_pos_x = clatoms_pos[ip].x;
       clatoms_pos_y = clatoms_pos[ip].y;
       clatoms_pos_z = clatoms_pos[ip].z;
       for(iatm=1;iatm<=natm_tot;iatm++){
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_x[iatm]);
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_y[iatm]);
        printf("%d %d %d %g\n",iproc,ip,iatm,clatoms_pos_z[iatm]);
       }/*endfor*/
     }/*endfor*/
    }/*endif*/
    scanf("%d",&iii);
    Dbx_Barrier(comm_beads_forc);
  }/*endfor*/
#endif

/*========================================================================*/
   }/*end routine*/
/*========================================================================*/






