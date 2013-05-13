/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_coord_class_state                       */
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
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void comm_coord_class_state(CLASS *class, GENERAL_DATA *general_data, 
                            int ivx_flag, int ivnhc_flag,int ix_flag,int pimd_on)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

#include "../typ_defs/typ_mask.h"

 int ip,iproc,dest,tag,source,ipart,iii,np_use,i,j,ichain,ip_start;
 int len_nhc,nsec,iend,nsend,num_mall;
 double *x_scr;
 int natm_tot         = class->clatoms_info.natm_tot;
 int pi_beads         = class->clatoms_info.pi_beads;
 int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
 int myid             = class->communicate.myid;
 int myid_state       = class->communicate.myid_state;
 int myid_bead        = class->communicate.myid_bead;
 int myid_bead_prime  = class->communicate.myid_bead_prime;
 int myid_forc        = class->communicate.myid_forc;
 int np_states        = class->communicate.np_states;
 int np_forc          = class->communicate.np_forc;
 int num_proc         = class->communicate.np;
 int num_nhc          = class->therm_info_class.num_nhc;
 int num_nhc_bead     = class->therm_info_bead.num_nhc;
 MPI_Comm world       = class->communicate.world;
 MPI_Comm comm_states = class->communicate.comm_states;
 MPI_Comm comm_forc   = class->communicate.comm_forc;

/*======================================================================*/


 if(pi_beads==1){num_nhc_bead=0;}
 if(num_nhc>0||num_nhc_bead>0){
   num_mall = MAX(num_nhc,num_nhc_bead);
   x_scr   = (double *)cmalloc(num_mall*sizeof(double))-1;
 }/*endif*/

 if(np_forc>1){np_use=np_forc;}
 if(np_states>1){np_use=np_states;}


if(num_proc>1 ){
if(ix_flag==0){
 if(myid_state==0&&myid_forc==0){
   for(iproc=1;iproc<np_use;iproc++){
     dest = myid + iproc;
     for(ip=1;ip<=pi_beads_proc;ip++){
       Send(&(class->clatoms_pos[ip].x[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                     world);
       Send(&(class->clatoms_pos[ip].y[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                     world);
       Send(&(class->clatoms_pos[ip].z[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                     world);
     }/*endfor*/
   }/*endfor*/
 }/*endif*/
 if(myid_state!=0||myid_forc!=0){
   for(iproc=1;iproc<np_use;iproc++){
     if(myid_state==iproc||myid_forc==iproc){
       source = myid - iproc;
       for(ip=1;ip<=pi_beads_proc;ip++){
         Recv(&(class->clatoms_pos[ip].x[0]),natm_tot+1,MPI_DOUBLE,source,ip,
                                                                       world);
         Recv(&(class->clatoms_pos[ip].y[0]),natm_tot+1,MPI_DOUBLE,source,ip,
                                                                       world);
         Recv(&(class->clatoms_pos[ip].z[0]),natm_tot+1,MPI_DOUBLE,source,ip,
                                                                       world);
       }/*endfor*/
     }/*endif:myid_state*/
   }/*endfor*/
  }/*endif*/
 }/*endif*/

if(ivx_flag == 0){
  if(np_forc>1){
    if(myid_forc==0){
      for(iproc=1;iproc<np_use;iproc++){
        dest = myid + iproc;
        for(ip=1;ip<=pi_beads_proc;ip++){
          Send(&(class->clatoms_pos[ip].vx[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                       world);
          Send(&(class->clatoms_pos[ip].vy[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                       world);
          Send(&(class->clatoms_pos[ip].vz[0]),natm_tot+1,MPI_DOUBLE,dest,ip,
                                                                       world);
        }/*endfor*/
      }/*endfor*/
    }/*endif*/
    if(myid_forc!=0){
      for(iproc=1;iproc<np_use;iproc++){
        if(myid_forc==iproc){
          source = myid - iproc;
          for(ip=1;ip<=pi_beads_proc;ip++){
            Recv(&(class->clatoms_pos[ip].vx[0]),natm_tot+1,MPI_DOUBLE,
                                                           source,ip,world);
            Recv(&(class->clatoms_pos[ip].vy[0]),natm_tot+1,MPI_DOUBLE,
                                                           source,ip,world);
            Recv(&(class->clatoms_pos[ip].vz[0]),natm_tot+1,MPI_DOUBLE,
                                                           source,ip,world);
          }/*endfor*/
        }/*endif:myid_forc*/
      }/*endfor*/
    }/*endif : myid_forc*/
  }/*endif : np_forc>1*/
}/*endif : no classical vel sampling*/


if(class->therm_info_class.num_nhc>0){

if(np_forc>1){
  if(myid_forc==0){
    for(iproc=1;iproc<np_use;iproc++){
        dest = myid + iproc;
        num_nhc = class->therm_info_class.num_nhc;
        len_nhc = class->therm_info_class.len_nhc;
       for(i=1;i<=len_nhc;i++){
        for(j=1;j<=num_nhc;j++){
         x_scr[j] = class->therm_class.x_nhc[i][j];
	}/*endfor*/
         Send(&(x_scr[1]),num_nhc,MPI_DOUBLE,dest,0,world);
       }/*endfor*/

        Send(&(general_data->baro.x_vol_nhc[1]),len_nhc, MPI_DOUBLE,dest,0,
                                                                       world);
        Send(&(general_data->baro.x_lnv),1, MPI_DOUBLE,dest,0,world);
        num_nhc = class->therm_info_bead.num_nhc;
        len_nhc = class->therm_info_bead.len_nhc;
        if(pimd_on==1){
          for(ip=1;ip<=pi_beads_proc;ip++){
            for(i=1;i<=len_nhc;i++){
             for(j=1;j<=num_nhc;j++){
              x_scr[j] = class->therm_bead[ip].x_nhc[i][j];
             }/*endfor*/
              Send(&(x_scr[1]),num_nhc,MPI_DOUBLE,dest,0,world);
            }/*endfor*/
          }/*endfor*/
	}/*endif : pimd_on*/
    }/*endfor*/
  }/*endif*/

  if(myid_forc!=0){
    for(iproc=1;iproc<np_use;iproc++){
      if(myid_forc==iproc){
        source = myid - iproc;
        num_nhc = class->therm_info_class.num_nhc;
        len_nhc = class->therm_info_class.len_nhc;
       for(i=1;i<=len_nhc;i++){
        Recv(&(x_scr[1]),num_nhc,MPI_DOUBLE,source,0,world);
        for(j=1;j<=num_nhc;j++){
         class->therm_class.x_nhc[i][j] = x_scr[j];
	}/*endfor*/
       }/*endfor*/
        Recv(&(general_data->baro.x_vol_nhc[1]),len_nhc,MPI_DOUBLE,source,0,
                                                                      world);
        Recv(&(general_data->baro.x_lnv),1,MPI_DOUBLE,source,0,world);
        num_nhc = class->therm_info_bead.num_nhc;
        len_nhc = class->therm_info_bead.len_nhc;
        if(pimd_on==1){
          for(ip=1;ip<=pi_beads_proc;ip++){
           for(i=1;i<=len_nhc;i++){
             Recv(&(x_scr[1]),num_nhc,MPI_DOUBLE,source,0,world);
            for(j=1;j<=num_nhc;j++){
             class->therm_bead[ip].x_nhc[i][j] = x_scr[j];
            }/*endfor*/
           }/*endfor*/

          }/*endfor*/
	}/*endif : pimd_on*/
      }/*endif:myid_forc*/
    }/*endfor*/
  }/*endif : myid_forc*/
}/*endif : np_forc>1*/


if(ivnhc_flag == 0){
  if(np_forc>1){
    ip_start=1;
    if(myid_bead_prime==0){
      ip_start=2;
      num_nhc = class->therm_info_class.num_nhc;
      len_nhc = class->therm_info_class.len_nhc;
      for(ichain=1;ichain<=len_nhc;ichain++){
          if(myid_forc==0){
            for(j=1;j<=num_nhc;j++){
              x_scr[j] = class->therm_class.v_nhc[ichain][j];
            }/*endfor*/
	  }/*endif*/
          Bcast(&(x_scr[1]),num_nhc,MPI_DOUBLE,0,comm_forc);
          if(myid_forc!=0){
                for(j=1;j<=num_nhc;j++){
                  class->therm_class.v_nhc[ichain][j]=x_scr[j];
                }/*endfor*/
	   }/*endif*/
      }/*endfor : ichain*/
      Bcast(&(general_data->baro.v_vol_nhc[1]),len_nhc, MPI_DOUBLE,
                                                              0,comm_forc);
      Bcast(&(general_data->baro.v_lnv),1, MPI_DOUBLE,0,comm_forc);
      Bcast(&(general_data->par_rahman.vgmat[1]),9, MPI_DOUBLE,0,comm_forc);
   }/*endif : myid_bead_prime*/


          if(pimd_on==1){
  
            num_nhc = class->therm_info_bead.num_nhc;
            len_nhc = class->therm_info_bead.len_nhc;
            for(ip=ip_start;ip<=pi_beads_proc;ip++){
              for(ichain=1;ichain<=len_nhc;ichain++){
                   if(myid_forc==0){
                     for(j=1;j<=num_nhc;j++){
                        x_scr[j] = 
                            class->therm_bead[ip].v_nhc[ichain][j];
                     }/*endfor*/
	           }/*endif*/
                   Bcast(&(x_scr[1]),num_nhc,MPI_DOUBLE,0,comm_forc);
                   if(myid_forc!=0){
                     for(j=1;j<=num_nhc;j++){
                        class->therm_bead[ip].v_nhc[ichain][j]=
                               x_scr[j];
                     }/*endfor*/
	           }/*endif : myid_forc*/
               }/*endfor : ichain*/
             }/*endfor : ip*/
	  }/*endif : pimd_on*/
    }/*endif : np_forc >1*/

}/*endif : no NHC vel sampling*/

}/*endif : num_nhc*/

}/*endif : num_proc > 1 */


 if(np_forc>1){
  Bcast(&(class->interact.spread),1,MPI_DOUBLE,0,comm_forc);
  Bcast(&(class->interact.spread_now),1,MPI_DOUBLE,0,comm_forc);
 }/*endif*/

 if(num_nhc>0||num_nhc_bead>0){
   cfree(&(x_scr[1]));
 }/*endif*/

/*========================================================================*/
} /* end routine */
/*==========================================================================*/








