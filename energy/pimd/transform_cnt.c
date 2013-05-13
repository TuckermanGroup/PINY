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

#define DEBUG_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void path_integral_init(CLATOMS_INFO *clatoms_info,
                        CLATOMS_POS *clatoms_pos,CLATOMS_TRAN *clatoms_tran, 
                        GHOST_ATOMS *ghost_atoms, SIMOPTS *simopts,
                        ATOMMAPS *atommaps, COMMUNICATE *communicate)

/*==========================================================================*/
/*            Begin Routine */
 {/*begin routine */
/*==========================================================================*/
/*           Local variable declarations                                */

   int pi_beads,pi_beads2,ip,ipart,i,igloc,natm_tot,ighost,nghost;
   int pi_beads_proc,natm_use;
   int init,iii;
   int one1,one2,one3,one4,one5,one6;
   int nwork;
   int lwork1,lwork2;
   int sign,ioff;
   int pi_atm_proc,pi_proc_rem;
   double sfact;
   double dpi_beads,*mass;
   double beta,tau,*prekf,temp,pre,tpip,gamma;
   double *bead_mass,*veig,arg;   
   double *rat1,*rat2;

   int pi_md_typ        = simopts->pi_md_typ;
   int *ighost_map      = ghost_atoms->ighost_map;

   int myid             = communicate->myid_bead;
   int num_proc         = communicate->np_beads;
   MPI_Comm comm_beads  = communicate->comm_beads;

   int pi_beads_proc_st = clatoms_info->pi_beads_proc_st;
   int myatm_start      = clatoms_info->myatm_start;
   int myatm_end        = clatoms_info->myatm_end;

   int nfreeze          = atommaps->nfreeze;
   int *freeze_map      = atommaps->freeze_map;
   int pimd_freez_typ   = atommaps->pimd_freez_typ;

   double pi_temperature = clatoms_info->pi_temperature;

/*==========================================================================*/
/* 0) Assign local pointers */

   clatoms_tran->path_eig = 
              (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
   natm_tot      = clatoms_info->natm_tot;
   nghost        = ghost_atoms->nghost_tot;
   pi_beads      = clatoms_info->pi_beads;
   pi_beads_proc = clatoms_info->pi_beads_proc;
   mass          = clatoms_info->mass;
   gamma         = clatoms_info->gamma_adb;
   prekf         = clatoms_info->prekf;
   veig          = clatoms_tran->path_eig;
   clatoms_tran->nwork  = 10000;
   nwork         = clatoms_tran->nwork;

   clatoms_tran->ifax    = (int *)cmalloc(13*sizeof(int))-1;
   clatoms_tran->work    =  (double *)cmalloc(nwork*sizeof(double))-1;
   clatoms_tran->work2   = (double *)cmalloc(nwork*sizeof(double))-1;
   clatoms_tran->work3   =  (double *)cmalloc(nwork*sizeof(double))-1;
   clatoms_tran->work4   = (double *)cmalloc(nwork*sizeof(double))-1;
   clatoms_tran->x_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;
   clatoms_tran->y_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;
   clatoms_tran->z_trans = (double *)cmalloc(2*pi_beads*sizeof(double))-1;

   dpi_beads = (double)pi_beads;

/*==========================================================================*/
/* I.V) FFT initialization                                                  */

 if(pi_md_typ == 2){

#ifdef SGI_COMPLIB
   ZFFT1DI(&pi_beads,&(clatoms_tran->work[1]));
#endif
#ifdef IBM_ESSL
    init=1;
    one1=one2=one3=one4=one5=one6=1;
    sign=1;
    lwork1=lwork2=nwork;
    sfact = 1.0;
    DCFT(&init,&(clatoms_tran->x_trans[1]),&one1,&one2,
         &(clatoms_tran->x_trans[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&sfact,&(clatoms_tran->work[1]),
         &lwork1,&(clatoms_tran->work2[1]),&lwork2);
    init=1;
    one1=one2=one3=one4=one5=one6=1;
    lwork1=lwork2=nwork;
    sfact = 1.0;
    sign = -1;
    DCFT(&init,&(clatoms_tran->x_trans[1]),&one1,&one2,
         &(clatoms_tran->x_trans[1]),
         &one3,&one4,&pi_beads,&one5,&sign,&sfact,&(clatoms_tran->work3[1]),
         &lwork1,&(clatoms_tran->work4[1]),&lwork2);
#endif
#ifdef C90
   CFTFAX(&pi_beads,&clatoms_tran->ifax[1],&clatoms_tran->work2[1]);
#endif

 }/*endif : centroid*/

/*==========================================================================*/
/* I) Malloc parallel arrays */

   if(num_proc>1){
     natm_use = (myatm_end-myatm_start+1);
     pi_atm_proc = (natm_use)/num_proc;

     if( (natm_use % num_proc) !=0){pi_atm_proc++;}
     clatoms_tran->pi_atm_proc = pi_atm_proc;
     pi_proc_rem = (natm_use%num_proc);
  
     clatoms_tran->pi_atm_proc_use = pi_atm_proc;
     if(myid>=pi_proc_rem&&pi_proc_rem!=0){
       clatoms_tran->pi_atm_proc_use = pi_atm_proc-1;  
     }/*endif*/
     clatoms_tran->pi_proc_rem = pi_proc_rem;
     clatoms_tran->x_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
     clatoms_tran->y_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
     clatoms_tran->z_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
     clatoms_tran->xt_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
     clatoms_tran->yt_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;
     clatoms_tran->zt_temp = 
                (double *)malloc((pi_atm_proc*pi_beads)*sizeof(double))-1;

     clatoms_tran->sendcounts    = (int *)malloc(num_proc*sizeof(int))-1;
     clatoms_tran->senddspls    = (int *)malloc(num_proc*sizeof(int))-1;
     clatoms_tran->recvcounts    = (int *)malloc(num_proc*sizeof(int))-1;
     clatoms_tran->recvdspls    = (int *)malloc(num_proc*sizeof(int))-1;
   }/*endif*/

/*==========================================================================*/
/* II) Set pre kinetic energy constants */

  for(i=1;i<=natm_tot;i++){
    beta = BOLTZ/pi_temperature;
    tau = BOLTZ/(pi_temperature*dpi_beads);
    prekf[i]  = mass[i]/(2.0*tau*beta);
  }/*endfor*/
   
  for(ighost=1;ighost <= nghost;ighost++){
    igloc = ighost_map[ighost];
    prekf[igloc] = 0.0;
  }/*endfor*/

  if(pimd_freez_typ==2){
    for(i=1;i <= nfreeze;i++){
      igloc = freeze_map[i];
      prekf[igloc] = 0.0;      
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
/* III) Set centroid eigenvalues and use them to define the masses          */

   if(pi_md_typ == 2){
/* setup eigenvalues in general form    */
     pi_beads2 = pi_beads/2;
     tpip = 2.0*M_PI/dpi_beads;
     pre  = 4.0*dpi_beads;
     veig[1]  = 1.0;      /*True eigenvalue is 0*/
     if(pi_beads>1){veig[pi_beads] = pre;}
     for(ip=2;ip<=pi_beads2;ip++){
      arg = tpip*(double)(ip-1);
      temp = (1.0-cos(arg))*pre;
      veig[2*ip-2] = temp;
      veig[2*ip-1] = temp;
     }/*endfor*/
     for(ip=1;ip<=pi_beads_proc;ip++){
      ioff = ip + pi_beads_proc_st-1;
      bead_mass = clatoms_pos[ip].mass;
      for(ipart=1;ipart<=natm_tot;ipart++){
        if(ip==1&&pi_beads_proc_st==1){
          bead_mass[ipart] = mass[ipart]*veig[ioff];
	}else{
          bead_mass[ipart] = mass[ipart]*veig[ioff]/(gamma*gamma);
	}/*endif*/
      }/*endfor*/
     }/*endfor*/
     veig[1]  = 0.0;      /*True eigenvalue is 0*/
   }/*endif*/

/*==========================================================================*/
/* III) Set staging eigenvalues and use them to define the masses          */

   if(pi_md_typ == 1){
     clatoms_tran->rat1_stag = 
              (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
     clatoms_tran->rat2_stag = 
              (double *)cmalloc(clatoms_info->pi_beads*sizeof(double))-1;
     rat1 = clatoms_tran->rat1_stag;
     rat2 = clatoms_tran->rat2_stag;
     for(ip=2;ip<=pi_beads;ip++){
       rat1[ip] = ((double)(ip-1))/((double)(ip));
       rat2[ip] = 1.0/((double)(ip));
     }/*endfor*/
     veig[1]  = 1.0; 
     for(ip=2;ip<=pi_beads;ip++){
      veig[ip] = ((double)(ip))/((double)(ip-1));
     }/*endfor*/
     for(ip=1;ip<=pi_beads_proc;ip++){
       ioff = ip + pi_beads_proc_st - 1;
       bead_mass = clatoms_pos[ip].mass;
       for(ipart=1;ipart<=natm_tot;ipart++){
        if(ip==1&&pi_beads_proc_st==1){
          bead_mass[ipart] = mass[ipart]*veig[ioff];
	}else{
          bead_mass[ipart] = mass[ipart]*veig[ioff]/(gamma*gamma);
	}/*endif*/
       }/*endfor*/
      }/*endfor*/
      veig[1]  = 0.0;      /*True eigenvalue is 0*/
   }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_mode_cent(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran )

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  int ifound;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init;
  double sfact;
  double *x_c,*y_c,*z_c,anorm;
  int ip_ind,ip_ind1,ip_ind2;
  int sign,lda;

/*==========================================================================*/
/* I) Assign local pointers */

  np2 = pi_beads/2;
  x_c = clatoms_tran->x_trans;
  y_c = clatoms_tran->y_trans;
  z_c = clatoms_tran->z_trans;

/*==========================================================================*/
/* II) Calculate the normal modes */

  if(pi_beads!=1){
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    for(ip=1;ip<=pi_beads;ip++){
     ip_ind = 2*ip-1;
     x_c[ip_ind]   = clatoms_pos[ip].x[ipart];
     y_c[ip_ind]   = clatoms_pos[ip].y[ipart];
     z_c[ip_ind]   = clatoms_pos[ip].z[ipart];
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
    clatoms_pos[1].x[ipart]  = x_c[1];
    clatoms_pos[1].y[ipart]  = y_c[1];
    clatoms_pos[1].z[ipart]  = z_c[1];
    clatoms_pos[pi_beads].x[ipart]  = x_c[(pi_beads+1)];
    clatoms_pos[pi_beads].y[ipart]  = y_c[(pi_beads+1)];
    clatoms_pos[pi_beads].z[ipart]  = z_c[(pi_beads+1)];
    for(ip=2;ip<=np2;ip++){
     ip_ind1 = 2*ip-1;
     ip_ind2 = 2*ip-2;
     clatoms_pos[ip_ind2].x[ipart] = x_c[ip_ind1];
     clatoms_pos[ip_ind2].y[ipart] = y_c[ip_ind1];
     clatoms_pos[ip_ind2].z[ipart] = z_c[ip_ind1];
    }/*endfor*/
    for(ip=2;ip<=np2;ip++){
     ip_ind1 = 2*ip-1;
     ip_ind2 = 2*ip;
     clatoms_pos[ip_ind1].x[ipart] = x_c[ip_ind2];
     clatoms_pos[ip_ind1].y[ipart] = y_c[ip_ind2];
     clatoms_pos[ip_ind1].z[ipart] = z_c[ip_ind2];
    }/*endfor*/
   }/*endfor*/
  }/*endif*/

#ifdef DEBUG
  printf("DEBUGGING MODE TRANSFORM\n");
  for(ip=1;ip<=pi_beads;ip++){
/*    printf("pos[%d].x[1]=%g \n",ip,clatoms_pos[ip].x[1]);*/
    printf("pos[%d].y[1]=%g .vy[1]=%g.fy[1]=%g\n",ip,
                                  clatoms_pos[ip].y[1],
                                  clatoms_pos[ip].vy[1],
                                  clatoms_pos[ip].fy[1]);
/*    printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);*/
  }/*endfor*/
  printf("Enter an integer: ");scanf("%d",&ip);
#endif
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_pos_cent(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran )

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2,iii;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init;
  double sfact;
  double *x_c,*y_c,*z_c;
  double anorm;
  int ifound;
  int ip_ind,ip_ind2,ip_ind1;
  int ip_indm2,ip_indm1;
  int sign,lda;

/*==========================================================================*/
/* I) Assign local pointers */

  np2 = pi_beads/2;

  x_c = clatoms_tran->x_trans;
  y_c = clatoms_tran->y_trans;
  z_c = clatoms_tran->z_trans;

/*==========================================================================*/
/* II) Get cartesian positions */

  if(pi_beads!=1){
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
     x_c[1]   = clatoms_pos[1].x[ipart];
     y_c[1]   = clatoms_pos[1].y[ipart];
     z_c[1]   = clatoms_pos[1].z[ipart];
     x_c[2] = 0.0;
     y_c[2] = 0.0;
     z_c[2] = 0.0;
     x_c[(pi_beads+1)] = clatoms_pos[pi_beads].x[ipart];
     y_c[(pi_beads+1)] = clatoms_pos[pi_beads].y[ipart];
     z_c[(pi_beads+1)] = clatoms_pos[pi_beads].z[ipart];
     x_c[(pi_beads+2)] = 0.0;
     y_c[(pi_beads+2)] = 0.0;
     z_c[(pi_beads+2)] = 0.0;
		  
     for(ip=2;ip<=np2;ip++){
      ip_ind1 = 2*ip-1;
      ip_ind2 = 2*ip-2;
      x_c[ip_ind1] = clatoms_pos[ip_ind2].x[ipart];
      y_c[ip_ind1] = clatoms_pos[ip_ind2].y[ipart];
      z_c[ip_ind1] = clatoms_pos[ip_ind2].z[ipart];
     }/*endfor*/
     for(ip=2;ip<=np2;ip++){
      ip_ind1 = 2*ip-1;
      ip_ind2 = 2*ip;
      x_c[ip_ind2] = clatoms_pos[ip_ind1].x[ipart]; 
      y_c[ip_ind2] = clatoms_pos[ip_ind1].y[ipart];
      z_c[ip_ind2] = clatoms_pos[ip_ind1].z[ipart]; 
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
       clatoms_pos[ip].x[ipart] = x_c[ip_ind];
       clatoms_pos[ip].y[ipart] = y_c[ip_ind];
       clatoms_pos[ip].z[ipart] = z_c[ip_ind];
     }/*endfor*/
   }/*endfor:ipart*/
  }/*endif*/
#ifdef DEBUG
  printf("DEBUGGING POS TRANSFORM\n");
  for(ip=4;ip<=pi_beads;ip+=4){
    printf("pos[%d].x[1]=%g\n",ip,clatoms_pos[ip].x[1]);
    printf("pos[%d].y[1]=%g\n",ip,clatoms_pos[ip].y[1]);
    printf("pos[%d].z[1]=%g\n",ip,clatoms_pos[ip].z[1]);
  }/*endfor*/
  printf("Enter an integer: ");scanf("%d",&ip);
#endif


/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_force_cent(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran )

/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,iopt,ierr,npp,incl,n,incn,np2,iii;
  int np2p1,ip2m1, ip2m2;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  int ifound;
  int one1,one2,one3,one4,one5,one6;
  int nwork;
  int lwork1,lwork2;
  int init;
  double sfact;
  double *x_c,*y_c,*z_c,anorm;
  int ip_ind,ip_ind1,ip_ind2;
  int sign,lda;

/*==========================================================================*/
/* I) Assign local pointers */

  np2  = pi_beads/2;
  x_c = clatoms_tran->x_trans;
  y_c = clatoms_tran->y_trans;
  z_c = clatoms_tran->z_trans;

/*==========================================================================*/
/* II) Get mode forces */

  if(pi_beads!=1){
   for(ipart=myatm_start;ipart<=myatm_end;ipart++){
    for(ip=1;ip<=pi_beads;ip++){
     ip_ind = 2*ip-1;
     x_c[ip_ind]   = clatoms_pos[ip].fx[ipart];
     y_c[ip_ind]   = clatoms_pos[ip].fy[ipart];
     z_c[ip_ind]   = clatoms_pos[ip].fz[ipart];
    }/*endfor*/
   for(ip=1;ip<=pi_beads;ip++){
    ip_ind = 2*ip;
    x_c[ip_ind] = 0.0;
    y_c[ip_ind] = 0.0;
    z_c[ip_ind] = 0.0;
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
   ZFFTS(&x_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   ZFFTS(&y_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   ZFFTS(&z_c[1],&npp,&incl,&n,&incn,&iopt,&ierr);
   anorm = (double)(pi_beads);
   for(ip=1;ip<=2*pi_beads;ip++){
     x_c[ip] *= anorm;
     y_c[ip] *= anorm;
     z_c[ip] *= anorm;
   }/*endfor*/
#endif
#ifdef SGI_COMPLIB
    ifound = 1;   
    sign   = 1;
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
     clatoms_pos[1].fx[ipart] = x_c[1];
     clatoms_pos[1].fy[ipart] = y_c[1];
     clatoms_pos[1].fz[ipart] = z_c[1];
     clatoms_pos[pi_beads].fx[ipart] = x_c[(pi_beads+1)];
     clatoms_pos[pi_beads].fy[ipart] = y_c[(pi_beads+1)];
     clatoms_pos[pi_beads].fz[ipart] = z_c[(pi_beads+1)];
     for(ip=2;ip<=np2;ip++){
      ip_ind2 = 2*ip - 2;
      ip_ind1 = 2*ip - 1;
      clatoms_pos[ip_ind2].fx[ipart] = x_c[ip_ind1]*2.0;
      clatoms_pos[ip_ind2].fy[ipart] = y_c[ip_ind1]*2.0;
      clatoms_pos[ip_ind2].fz[ipart] = z_c[ip_ind1]*2.0;
     }/*endfor*/
     for(ip=2;ip<=np2;ip++){
      ip_ind1 = 2*ip - 1;
      ip_ind2 = 2*ip;
      clatoms_pos[ip_ind1].fx[ipart] = x_c[ip_ind2]*2.0;
      clatoms_pos[ip_ind1].fy[ipart] = y_c[ip_ind2]*2.0;
      clatoms_pos[ip_ind1].fz[ipart] = z_c[ip_ind2]*2.0;
     }/*endfor*/
    }/*endfor:ipart*/
  }/*endif*/
#ifdef DEBUG
  printf("DEBUGGING FORCE TRANSFORM\n");
  for(ip=4;ip<=pi_beads;ip+=4){
    printf("pos[%d].fx[1]=%g\n",ip,clatoms_pos[ip].fx[myatm_start]);
    printf("pos[%d].fy[1]=%g\n",ip,clatoms_pos[ip].fy[myatm_start]);
    printf("pos[%d].fz[1]=%g\n",ip,clatoms_pos[ip].fz[myatm_start]);
  }/*endfor*/
  for(ip=4;ip<=pi_beads;ip+=4){
    printf("pos[%d].fx[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fx[myatm_end]);
    printf("pos[%d].fy[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fy[myatm_end]);
    printf("pos[%d].fz[%d]=%g\n",ip,natm_tot,clatoms_pos[ip].fz[myatm_end]);
  }/*endfor*/
  printf("Enter an integer: ");scanf("%d",&ip);
#endif

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/




