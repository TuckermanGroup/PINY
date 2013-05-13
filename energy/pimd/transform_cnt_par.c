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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_mode_cent_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2,ioff,ipart_ind,iproc;
  int pi_beads      = clatoms_info->pi_beads;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int natm_tot = clatoms_info->natm_tot;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int ifound;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init,iii;
  double sfact;
  double *x_c,*y_c,*z_c,anorm;
  double *x_temp,*y_temp,*z_temp;
  int ip_ind,ip_ind1,ip_ind2;
  int sign,lda;
  int num_proc        = communicate->np_beads;
  int myid            = communicate->myid_bead;
  MPI_Comm comm_beads = communicate->comm_beads;

/*==========================================================================*/
/* I) Assign local pointers */

  np2 = pi_beads/2;
  x_c = clatoms_tran->x_trans;
  y_c = clatoms_tran->y_trans;
  z_c = clatoms_tran->z_trans;
  x_temp = clatoms_tran->x_temp;
  y_temp = clatoms_tran->y_temp;
  z_temp = clatoms_tran->z_temp;

/*==========================================================================*/
/* II) Calculate the normal modes */

  
  if(pi_beads!=1){
   iopt = 1;
   pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,communicate,iopt);
   ioff = 0;
   for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
    for(ip=1;ip<=pi_beads;ip++){
     ip_ind = 2*ip-1;
     ipart_ind = ip+ioff;
     x_c[ip_ind]   = x_temp[ipart_ind];
     y_c[ip_ind]   = y_temp[ipart_ind];
     z_c[ip_ind]   = z_temp[ipart_ind];
    }/*endfor*/
    for(ip=1;ip<=pi_beads;ip++){
     ip_ind = 2*ip;
     x_c[ip_ind] = 0.0;
     y_c[ip_ind] = 0.0;
     z_c[ip_ind] = 0.0;
    }/*endfor*/
/* i) calculate normal modes*/ 
    ifound = 0;   
#ifdef HP_VECLIB
    ifound = 1;   
    iopt = -1;
    ierr  = 0;
    npp = pi_beads;
    incl = 1;
    n = 1;
    incn = 1;
    ZFFTS(&x_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
    ZFFTS(&y_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
    ZFFTS(&z_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
#endif
#ifdef SGI_COMPLIB
    ifound = 1;   
    sign   = 1;
    npp = pi_beads;
    n = 1;
    incl = 1;
    lda = 1;
    ZFFT1D(&sign,&npp,&(x_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(y_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(z_c[1]),&incl,&(clatoms_tran->work[1]));
    anorm = 1.0/((double) (pi_beads));
    for(ip=1;ip<=2*pi_beads;ip++){
      x_c[ip] *= anorm;
      y_c[ip] *= anorm;
      z_c[ip] *= anorm;
    }/*endfor*/
#endif
#ifdef IBM_ESSL  
    ifound=1;
    init=0;
    one1=one2=one3=one4=one5=one6=1;
    lwork1=lwork2=clatoms_tran->nwork;
    anorm = 1.0;
    sign = -1;
    DCFT(&init,&(x_c[1]),&one1,&one2,&(x_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
    DCFT(&init,&(y_c[1]),&one1,&one2,&(y_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
    DCFT(&init,&(z_c[1]),&one1,&one2,&(z_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
    anorm = 1.0/((double) (pi_beads));
    for(ip=1;ip<=2*pi_beads;ip++){
      x_c[ip] *= anorm;
      y_c[ip] *= anorm;
      z_c[ip] *= anorm;
    }/*endfor*/
#endif
#ifdef C90
    ifound = 1;   
    sign   = 1;
    npp = pi_beads;
    n = 1;
    incl = 1;
    lda = 1;
    CFFT99(&x_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&y_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&z_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    anorm = 1.0/((double) (pi_beads));
    for(ip=1;ip<=2*pi_beads;ip++){
      x_c[ip] *= anorm;
      y_c[ip] *= anorm;
      z_c[ip] *= anorm;
    }/*endfor*/
#endif
    if(ifound==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Centroid transformation not implemented on your platform\n");
      printf("Contact Technical support: \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
/* ii) unload normal modes */
    ipart_ind = ioff+1;
    x_temp[ipart_ind]  = x_c[1];
    y_temp[ipart_ind]  = y_c[1];
    z_temp[ipart_ind]  = z_c[1];
    ipart_ind = ioff+pi_beads;
    x_temp[ipart_ind]  = x_c[(pi_beads+1)];
    y_temp[ipart_ind]  = y_c[(pi_beads+1)];
    z_temp[ipart_ind]  = z_c[(pi_beads+1)];
    for(ip=2;ip<=np2;ip++){
     ip_ind1 = 2*ip-1;
     ipart_ind = 2*ip-2+ioff;
     x_temp[ipart_ind] = x_c[ip_ind1];
     y_temp[ipart_ind] = y_c[ip_ind1];
     z_temp[ipart_ind] = z_c[ip_ind1];
    }/*endfor*/
    for(ip=2;ip<=np2;ip++){
     ipart_ind = 2*ip-1+ioff;
     ip_ind2 = 2*ip;
     x_temp[ipart_ind] = x_c[ip_ind2];
     y_temp[ipart_ind] = y_c[ip_ind2];
     z_temp[ipart_ind] = z_c[ip_ind2];
    }/*endfor*/
    ioff += pi_beads;
   }/*endfor:ipart*/
   iopt = 1;
   pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,communicate,iopt);
  }/*endif*/
#define DEBUG_MODE_OFF
#ifdef DEBUG_MODE
  if(myid==0){  printf("DEBUGGING MODE TRANSFORM\n");}
  Dbx_Barrier(comm_beads);
for(iproc=0;iproc<num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[1]=%g\n",ip,clatoms_pos[ip].x[1]);
      printf("pos[%d].y[1]=%g\n",ip,clatoms_pos[ip].y[1]);
      printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
for(iproc=0;iproc<num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].x[natm_tot]);
      printf("pos[%d].y[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].y[natm_tot]);
      printf("pos[%d].z[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].z[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_pos_cent_par(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2,iii,ioff,ipart_ind;
  int pi_beads = clatoms_info->pi_beads;
  int pi_beads_proc = clatoms_info->pi_beads_proc;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int natm_tot = clatoms_info->natm_tot;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init;
  double sfact;
  double *x_c,*y_c,*z_c;
  double *x_temp,*y_temp,*z_temp;
  double anorm;
  int ifound;
  int ip_ind,ip_ind2,ip_ind1;
  int ip_indm2,ip_indm1;
  int sign,lda;
  int num_proc        = communicate->np_beads;
  MPI_Comm comm_beads = communicate->comm_beads;

/*==========================================================================*/
/* I) Assign local pointers */

  np2 = pi_beads/2;

  x_c = clatoms_tran->x_trans;
  y_c = clatoms_tran->y_trans;
  z_c = clatoms_tran->z_trans;
  x_temp = clatoms_tran->x_temp;
  y_temp = clatoms_tran->y_temp;
  z_temp = clatoms_tran->z_temp;

/*==========================================================================*/
/* II) Get cartesian positions */

  if(pi_beads!=1){
   iopt = 1;
   pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,communicate,iopt);
   ioff=0;
   for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
     ipart_ind = ioff+1;
     x_c[1]   = x_temp[ipart_ind];
     y_c[1]   = y_temp[ipart_ind];
     z_c[1]   = z_temp[ipart_ind];
     x_c[2] = 0.0;
     y_c[2] = 0.0;
     z_c[2] = 0.0;
     ipart_ind = ioff+pi_beads;
     x_c[(pi_beads+1)] = x_temp[ipart_ind];
     y_c[(pi_beads+1)] = y_temp[ipart_ind];
     z_c[(pi_beads+1)] = z_temp[ipart_ind];
     x_c[(pi_beads+2)] = 0.0;
     y_c[(pi_beads+2)] = 0.0;
     z_c[(pi_beads+2)] = 0.0;
		  
     for(ip=2;ip<=np2;ip++){
      ip_ind1 = 2*ip-1;
      ipart_ind = 2*ip-2+ioff;
      x_c[ip_ind1] = x_temp[ipart_ind];
      y_c[ip_ind1] = y_temp[ipart_ind];
      z_c[ip_ind1] = z_temp[ipart_ind];
     }/*endfor*/
     for(ip=2;ip<=np2;ip++){
      ipart_ind = 2*ip-1+ioff;
      ip_ind2 = 2*ip;
      x_c[ip_ind2] = x_temp[ipart_ind]; 
      y_c[ip_ind2] = y_temp[ipart_ind]; 
      z_c[ip_ind2] = z_temp[ipart_ind]; 
     }/*endfor*/
     for(ip=2;ip<=np2;ip++){
      ip_ind1  = 2*ip-1;
      ip_ind2  = 2*ip;
      ip_indm1 = 2*pi_beads-2*ip+3;
      ip_indm2 = 2*pi_beads-2*ip+4;
      x_c[ip_indm1] =  x_c[ip_ind1];
      y_c[ip_indm1] =  y_c[ip_ind1];
      z_c[ip_indm1] =  z_c[ip_ind1];
      x_c[ip_indm2] = -x_c[ip_ind2];
      y_c[ip_indm2] = -y_c[ip_ind2];
      z_c[ip_indm2] = -z_c[ip_ind2];
     }/*endfor*/
/* calculate cartesian coordinates */
     ifound = 0;
#ifdef HP_VECLIB
     ifound = 1;
     iopt = 1;
     ierr  = 0;
     npp = pi_beads;
     incl = 1;
     n = 1;
     incn = 1;
     ZFFTS(&x_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
     ZFFTS(&y_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
     ZFFTS(&z_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
#endif
#ifdef SGI_COMPLIB
    ifound = 1;   
    sign   = -1;
    npp    = pi_beads;
    n      = 1;
    incl   = 1;
    lda    = 1;
    ZFFT1D(&sign,&npp,&(x_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(y_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(z_c[1]),&incl,&(clatoms_tran->work[1]));
#endif
#ifdef IBM_ESSL
    ifound=1;
    init=0;
    one1=one2=one3=one4=one5=one6=1;
    lwork1=lwork2=clatoms_tran->nwork;
    anorm = 1.0;
    sign = 1;
    DCFT(&init,&(x_c[1]),&one1,&one2,&(x_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work[1]),
         &lwork1,&(clatoms_tran->work2[1]),&lwork2);
    DCFT(&init,&(y_c[1]),&one1,&one2,&(y_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work[1]),
         &lwork1,&(clatoms_tran->work2[1]),&lwork2);
    DCFT(&init,&(z_c[1]),&one1,&one2,&(z_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work[1]),
         &lwork1,&(clatoms_tran->work2[1]),&lwork2);
#endif
#ifdef C90
    ifound = 1;   
    sign   = -1;
    npp = pi_beads;
    n = 1;
    incl = 1;
    lda = 1;
    CFFT99(&x_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&y_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&z_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
#endif
    if(ifound==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Centroid transformation not implemented on your platform\n");
      printf("Contact Technical support: \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
/* unload cartesian coordinates */
     for(ip=1;ip<=pi_beads;ip++){
       ip_ind = 2*ip-1;
       ipart_ind = ip+ioff;
       x_temp[ipart_ind] = x_c[ip_ind];
       y_temp[ipart_ind] = y_c[ip_ind];
       z_temp[ipart_ind] = z_c[ip_ind];
     }/*endfor*/
     ioff += pi_beads;
   }/*endfor:ipart*/
   iopt = 1;
   pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,communicate,iopt);
  }/*endif*/
#ifdef DEBUG
  printf("DEBUGGING POS TRANSFORM\n");
for(iproc=0;iproc<=num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[1]=%g\n",ip,clatoms_pos[ip].x[1]);
      printf("pos[%d].y[1]=%g\n",ip,clatoms_pos[ip].y[1]);
      printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);
    }/*endfor*/
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].x[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].x[natm_tot]);
      printf("pos[%d].y[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].y[natm_tot]);
      printf("pos[%d].z[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].z[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_force_cent_par(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran, COMMUNICATE *communicate)

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2,ioff,ipart_ind,iproc,iii;
  int np2p1,ip2m1, ip2m2;
  int pi_beads = clatoms_info->pi_beads;
  int pi_atm_proc_use = clatoms_tran->pi_atm_proc_use;
  int natm_tot = clatoms_info->natm_tot;
  int ifound;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init;
  double sfact;
  double *fx_c,*fy_c,*fz_c,anorm;
  double *fx_temp,*fy_temp,*fz_temp;
  int ip_ind,ip_ind1,ip_ind2;
  int sign,lda;
  int myid = communicate->myid_bead;
  int num_proc = communicate->np;
  MPI_Comm comm_beads = communicate->comm_beads;

/*==========================================================================*/
/* I) Assign local pointers */

  np2  = pi_beads/2;
  fx_c = clatoms_tran->x_trans;
  fy_c = clatoms_tran->y_trans;
  fz_c = clatoms_tran->z_trans;
  fx_temp = clatoms_tran->x_temp;
  fy_temp = clatoms_tran->y_temp;
  fz_temp = clatoms_tran->z_temp;

/*==========================================================================*/
/* II) Get mode forces */

  if(pi_beads!=1){
   iopt = 2;
   pimd_trans_comm_fwd(clatoms_info,clatoms_pos,clatoms_tran,
                              communicate,iopt);
   ioff = 0;
   for(ipart=1;ipart<=pi_atm_proc_use;ipart++){
    for(ip=1;ip<=pi_beads;ip++){
     ip_ind = 2*ip-1;
     ipart_ind = ip+ioff;
     fx_c[ip_ind]   = fx_temp[ipart_ind];
     fy_c[ip_ind]   = fy_temp[ipart_ind];
     fz_c[ip_ind]   = fz_temp[ipart_ind];
    }/*endfor*/
   for(ip=1;ip<=pi_beads;ip++){
    ip_ind = 2*ip;
    fx_c[ip_ind] = 0.0;
    fy_c[ip_ind] = 0.0;
    fz_c[ip_ind] = 0.0;
   }/*endfor*/
/* i) transform to normal mode forces */
   ifound = 0;    
#ifdef HP_VECLIB
   ifound = 1;    
   iopt = -1;
   ierr  = 0;
   npp = pi_beads;
   incl = 1;
   n = 1;
   incn = 1;
   ZFFTS(&fx_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   ZFFTS(&fy_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   ZFFTS(&fz_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   anorm = (double)(pi_beads);
   for(ip=1;ip<=2*pi_beads;ip++){
     fx_c[ip] *= anorm;
     fy_c[ip] *= anorm;
     fz_c[ip] *= anorm;
   }/*endfor*/
#endif
#ifdef SGI_COMPLIB
    ifound = 1;   
    sign   = 1;
    npp    = pi_beads;
    n      = 1;
    incl   = 1;
    lda    = 1;
    ZFFT1D(&sign,&npp,&(fx_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(fy_c[1]),&incl,&(clatoms_tran->work[1]));
    ZFFT1D(&sign,&npp,&(fz_c[1]),&incl,&(clatoms_tran->work[1]));
#endif
#ifdef IBM_ESSL
    ifound=1;
    init=0;
    one1=one2=one3=one4=one5=one6=1;
    lwork1=lwork2=clatoms_tran->nwork;
    anorm = 1.0;
    sign = -1;
    DCFT(&init,&(fx_c[1]),&one1,&one2,&(fx_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
    DCFT(&init,&(fy_c[1]),&one1,&one2,&(fy_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
    DCFT(&init,&(fz_c[1]),&one1,&one2,&(fz_c[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&anorm,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
#endif
#ifdef C90
    ifound = 1;   
    sign   = 1;
    npp = pi_beads;
    n = 1;
    incl = 1;
    lda = 1;
    CFFT99(&fx_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&fy_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
    CFFT99(&fz_c[1],&(clatoms_tran->work[1]),&(clatoms_tran->work2[1]),
           &(clatoms_tran->ifax[1]),&incl,&lda,&pi_beads,&n,&sign);
#endif
    if(ifound==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Centroid transformation not implemented on your platform\n");
      printf("Contact Technical support: \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
/* ii) put the forces in the correct vectors */
     ipart_ind = ioff+1;
     fx_temp[ipart_ind] = fx_c[1];
     fy_temp[ipart_ind] = fy_c[1];
     fz_temp[ipart_ind] = fz_c[1];
     ipart_ind = ioff+pi_beads;
     fx_temp[ipart_ind] = fx_c[(pi_beads+1)];
     fy_temp[ipart_ind] = fy_c[(pi_beads+1)];
     fz_temp[ipart_ind] = fz_c[(pi_beads+1)];
     for(ip=2;ip<=np2;ip++){
      ipart_ind = 2*ip - 2 + ioff;
      ip_ind1 = 2*ip - 1;
      fx_temp[ipart_ind] = fx_c[ip_ind1]*2.0;
      fy_temp[ipart_ind] = fy_c[ip_ind1]*2.0;
      fz_temp[ipart_ind] = fz_c[ip_ind1]*2.0;
     }/*endfor*/
     for(ip=2;ip<=np2;ip++){
      ipart_ind = 2*ip - 1 + ioff;
      ip_ind2 = 2*ip;
      fx_temp[ipart_ind] = fx_c[ip_ind2]*2.0;
      fy_temp[ipart_ind] = fy_c[ip_ind2]*2.0;
      fz_temp[ipart_ind] = fz_c[ip_ind2]*2.0;
     }/*endfor*/
     ioff += pi_beads;
    }/*endfor:ipart*/
    iopt = 2;
    pimd_trans_comm_bck(clatoms_info,clatoms_pos,clatoms_tran,
                              communicate,iopt);
  }/*endif*/
#ifdef DEBUG
  printf("DEBUGGING FORCE TRANSFORM\n");
for(iproc=0;iproc<=num_proc;iproc++){
  if(myid==iproc){
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].fx[1]=%g\n",ip,clatoms_pos[ip].fx[1]);
      printf("pos[%d].fy[1]=%g\n",ip,clatoms_pos[ip].fy[1]);
      printf("pos[%d].fz[1]=%g\n",ip,clatoms_pos[ip].fz[1]);
    }/*endfor*/
    for(ip=1;ip<=pi_beads_proc;ip++){
      printf("pos[%d].fx[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fx[natm_tot]);
      printf("pos[%d].fy[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fy[natm_tot]);
      printf("pos[%d].fz[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fz[natm_tot]);
    }/*endfor*/
    printf("Enter an integer: ");
  }/*endif*/
  scanf("%d",&ip);
  Dbx_Barrier(comm_beads);  
}/*endfor*/
#endif

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




