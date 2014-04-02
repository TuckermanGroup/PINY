/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: fft_package.c                                */
/*                                                                          */
/* 3D FFT modules                                                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"


#define BLOCK_PACK
#define DALLTOALL
#define SALLTOALL
#define POST_MAP
#define FFT
#define TRANS


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void para_fft_gen3d_fwd_to_r(double *zfft, double *zfft_tmp,
                             PARA_FFT_PKG3D *para_fft_pkg3d)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int iopt = 1,incl=1;
 int ier  = 0;
 int incn;
 int ndata;
 int i,j,k,m,ioff,joff,koff,kap,kbp,kcp;
 int iproc,iii,ka,ka_str,ka_end;
 int jb;
 int kb,kb_str,kb_end;
 int kc,kc_str,kc_end;
 int kstr,kend;
 int skc_fft_ka_now,ekc_fft_ka_now;
 int skb_fft_ka_now,ekb_fft_ka_now;
 int skc_fft_kb_now,ekc_fft_kb_now;

/*-----------------------------------------------------------------------*/
/*          Local pointers                                               */


/* Communicator */
 MPI_Comm comm = para_fft_pkg3d->comm;
 int num_proc  = para_fft_pkg3d->num_proc;
 int myid      = para_fft_pkg3d->myid; 
 int myidp1    = para_fft_pkg3d->myid+1; 

/* General Info */
 int ncoef_proc   = para_fft_pkg3d->ncoef_proc; 
 int ncoef_use    = para_fft_pkg3d->ncoef_use; 
 int igeneric_opt = para_fft_pkg3d->igeneric_opt;
 int nwork1       = para_fft_pkg3d->nwork1;
 int nwork2       = para_fft_pkg3d->nwork2;

/* FFT along kc  */
 int nkf3                = para_fft_pkg3d->nkf3; 
 int *sendcounts_fft_kc  = para_fft_pkg3d->sendcounts_fft_kc;
 int *senddspls_fft_kc   = para_fft_pkg3d->senddspls_fft_kc;
 int *recvcounts_fft_kc  = para_fft_pkg3d->recvcounts_fft_kc;
 int *recvdspls_fft_kc   = para_fft_pkg3d->recvdspls_fft_kc;
 int nfft_kc_proc        = para_fft_pkg3d->nfft_kc_proc;
 int nka_fft_kc          = para_fft_pkg3d->nka_fft_kc;
 int str_fft_kc_proc     = para_fft_pkg3d->str_fft_kc_proc;
 int end_fft_kc_proc     = para_fft_pkg3d->end_fft_kc_proc;
 int *sum_fft_kc_proc    = para_fft_pkg3d->sum_fft_kc_proc;
 int *ka_fft_kc          = para_fft_pkg3d->ka_fft_kc;
 int *kb_fft_kc          = para_fft_pkg3d->kb_fft_kc;
 int ska_fft_kc_proc     = para_fft_pkg3d->ska_fft_kc_proc;
 int eka_fft_kc_proc     = para_fft_pkg3d->eka_fft_kc_proc;
 int skb_fft_kc_proc     = para_fft_pkg3d->skb_fft_kc_proc;
 int ekb_fft_kc_proc     = para_fft_pkg3d->ekb_fft_kc_proc;
 int *nkb_fft_kc         = para_fft_pkg3d->nkb_fft_kc;
 int *ka_fft_kc_red      = para_fft_pkg3d->ka_fft_kc_red;
 int *ifax_c_f           = para_fft_pkg3d->ifax_c_f;
 double *work_1c_f       = para_fft_pkg3d->work_1c_f;
 double *work_2c_f       = para_fft_pkg3d->work_2c_f;

/* FFT along kb  */
 int nkf2                 = para_fft_pkg3d->nkf2; 
 int *sendcounts_fft_kb   = para_fft_pkg3d->sendcounts_fft_kb;
 int *senddspls_fft_kb    = para_fft_pkg3d->senddspls_fft_kb;
 int *recvcounts_fft_kb   = para_fft_pkg3d->recvcounts_fft_kb;
 int *recvdspls_fft_kb    = para_fft_pkg3d->recvdspls_fft_kb;
 int nfft_kb_proc         = para_fft_pkg3d->nfft_kb_proc;
 int str_fft_kb_proc      = para_fft_pkg3d->str_fft_kb_proc;
 int end_fft_kb_proc      = para_fft_pkg3d->end_fft_kb_proc;
 int ska_fft_kb_proc      = para_fft_pkg3d->ska_fft_kb_proc;
 int eka_fft_kb_proc      = para_fft_pkg3d->eka_fft_kb_proc;
 int skc_fft_kb_proc      = para_fft_pkg3d->skc_fft_kb_proc;
 int ekc_fft_kb_proc      = para_fft_pkg3d->ekc_fft_kb_proc;
 int *ska_fft_kb_proc_all = para_fft_pkg3d->ska_fft_kb_proc_all;
 int *eka_fft_kb_proc_all = para_fft_pkg3d->eka_fft_kb_proc_all;
 int *skc_fft_kb_proc_all = para_fft_pkg3d->skc_fft_kb_proc_all;
 int *ekc_fft_kb_proc_all = para_fft_pkg3d->ekc_fft_kb_proc_all;
 int *sum_fft_kb_proc     = para_fft_pkg3d->sum_fft_kb_proc;
 int *ifax_b_f            = para_fft_pkg3d->ifax_b_f;
 double *work_1b_f       = para_fft_pkg3d->work_1b_f;
 double *work_2b_f       = para_fft_pkg3d->work_2b_f;

/* FFT along ka  */
 int nkf1                 = para_fft_pkg3d->nkf1; 
 int *sendcounts_fft_ka   = para_fft_pkg3d->sendcounts_fft_ka;
 int *senddspls_fft_ka    = para_fft_pkg3d->senddspls_fft_ka;
 int *recvcounts_fft_ka   = para_fft_pkg3d->recvcounts_fft_ka;
 int *recvdspls_fft_ka    = para_fft_pkg3d->recvdspls_fft_ka;
 int nfft_ka_proc         = para_fft_pkg3d->nfft_ka_proc;
 int str_fft_ka_proc      = para_fft_pkg3d->str_fft_ka_proc;
 int end_fft_ka_proc      = para_fft_pkg3d->end_fft_ka_proc;
 int skb_fft_ka_proc      = para_fft_pkg3d->skb_fft_ka_proc;
 int ekb_fft_ka_proc      = para_fft_pkg3d->ekb_fft_ka_proc;
 int skc_fft_ka_proc      = para_fft_pkg3d->skc_fft_ka_proc;
 int ekc_fft_ka_proc      = para_fft_pkg3d->ekc_fft_ka_proc;
 int *skb_fft_ka_proc_all = para_fft_pkg3d->skb_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all = para_fft_pkg3d->ekb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all = para_fft_pkg3d->skc_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all = para_fft_pkg3d->ekc_fft_ka_proc_all;
 int *ifax_a_f            = para_fft_pkg3d->ifax_a_f;
 double *work_1a_f       = para_fft_pkg3d->work_1a_f;
 double *work_2a_f       = para_fft_pkg3d->work_2a_f;

/* MAPS */
 int *map_proc_post       = para_fft_pkg3d->map_proc_post;

/*==========================================================================*/
/* Transpose I  */

  if(num_proc>1){

    /*-----------------------------*/
    /* Communicate blocks          */
    /*-----------------------------*/

#ifdef SALLTOALL
    Alltoallv(&zfft[1],&(sendcounts_fft_kc[1]),&(senddspls_fft_kc[1]),
              MPI_DOUBLE,&zfft_tmp[1],&(recvcounts_fft_kc[1]),
              &(recvdspls_fft_kc[1]),MPI_DOUBLE,comm);
#endif

    /*-----------------------------*/
    /* Index Transpose plus spread */
    /*-----------------------------*/

    ndata = 2*nfft_kc_proc*nkf3; 
    for(i=1;i<=ndata;i++){zfft[i]=0.0;}
#ifdef POST_MAP
    ndata = recvdspls_fft_kc[num_proc]+recvcounts_fft_kc[num_proc];
    for(i=1;i<=ndata;i++){
     zfft[map_proc_post[i]]=zfft_tmp[i];
    }/*endfor*/
#endif

  }/*endif*/

/*==========================================================================*/
/* Perform FFT along ``c'' using zfft */

  incn = nkf3;
#ifdef FFT
  fft_gen1d(zfft,nkf3,incl,nfft_kc_proc,incn,iopt,&ier,
            work_1c_f,nwork1,work_2c_f,nwork2,ifax_c_f,igeneric_opt);
#endif

/*==========================================================================*/
/* Transpose II */


  if(num_proc>1){

   /*---------------------------------------------*/
   /* In place-block transpose                     */
   /*---------------------------------------------*/
   /*   i)Loop over procs=iproc                   */
   /*  ii)Loop over the ka range to send to iproc */
   /* iii)Loop over the kb range you have         */
   /*  iv)Loop over the kc range to send to iproc */
   /*   v)Store allowed/in range data             */

#ifdef BLOCK_PACK
   j=0;
   for(iproc=1;iproc<=num_proc;iproc++){
     ka_str = MAX(ska_fft_kc_proc,ska_fft_kb_proc_all[iproc]);
     ka_end = MIN(eka_fft_kc_proc,eka_fft_kb_proc_all[iproc]);
     for(ka=ka_str;ka<=ka_end;ka++){
       kb_str = (ka==ska_fft_kc_proc ? skb_fft_kc_proc : 1);
       kb_end = (ka==eka_fft_kc_proc ? ekb_fft_kc_proc : nkb_fft_kc[ka]);
       i = 2*sum_fft_kc_proc[ka]*nkf3;
       kc_str = (ka==ska_fft_kb_proc_all[iproc] ? 
                     skc_fft_kb_proc_all[iproc]  : 1);
       kc_end = (ka==eka_fft_kb_proc_all[iproc] ?  
                     ekc_fft_kb_proc_all[iproc]  : nkf3);
       for(kb=kb_str;kb<=kb_end;kb++){
         ioff   = i+2*(kc_str-1);
         for(kc=1;kc<=2*(kc_end-kc_str+1);kc++){
           zfft_tmp[(j+kc)] = zfft[(kc+ioff)];
         }/*endfor*/
         iii = 2*MAX(kc_end-kc_str+1,0);
         j+=iii;i+= 2*nkf3;
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/

#endif
  /*-----------------------------*/
  /* Communicate Blocks          */
  /*-----------------------------*/

#ifdef SALLTOALL
    Alltoallv(&zfft_tmp[1],&(sendcounts_fft_kb[1]),&(senddspls_fft_kb[1]),
              MPI_DOUBLE,&zfft[1],&(recvcounts_fft_kb[1]),
              &(recvdspls_fft_kb[1]),MPI_DOUBLE,comm);
#endif
  }/*endif*/

  /*-----------------------------*/
  /* Index Transpose plus spread */
  /*-----------------------------*/

   ndata = 2*nfft_kb_proc*nkf2;
   for(i=1;i<=ndata;i++){zfft_tmp[i]=0.0;}

#ifdef TRANS
   jb=0;for(i=1;i<=ska_fft_kb_proc-1;i++){jb+=nkb_fft_kc[i];}
   i=0;koff = 0;
   for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
     kc_str = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
     kc_end = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
     for(kb=1;kb<=nkb_fft_kc[ka];kb++){
      jb++;ioff = 2*kb_fft_kc[jb] - 1 + koff;
      for(kc=1,m=0;kc<=2*(kc_end-kc_str+1);kc+=2,m+=2*nkf2){
        zfft_tmp[(m+ioff)] = zfft[(i+kc)];
        zfft_tmp[(m+ioff+1)] = zfft[(i+kc+1)];
      }/*endfor*/
      i+=2*(kc_end-kc_str+1);
    }/*endfor*/
     koff+= 2*(kc_end-kc_str+1)*nkf2;
   }/*endfor*/
#endif

/*==========================================================================*/
/* Perform FFT along ``b'' using zfft_tmp */

  incn = nkf2;
#ifdef FFT
  fft_gen1d(zfft_tmp,nkf2,incl,nfft_kb_proc,incn,iopt,&ier,
            work_1b_f,nwork1,work_2b_f,nwork2,ifax_b_f,igeneric_opt);
#endif

/*==========================================================================*/
/* Transpose III */


 if(num_proc>1){

#ifdef BLOCK_PACK
   /*---------------------------------------------*/
   /* In place-block transpose                     */
   /*---------------------------------------------*/
   /*   i)Loop over ka's you have                 */
   /*  ii)Loop over procs=iproc                   */
   /* iii)Loop over the kc range to send to iproc */
   /*  iv)Loop over the kb range to send to iproc */
   /*   v)Store allowed/in range data            */

   j=0;
   for(iproc=1;iproc<=num_proc;iproc++){
     for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
       skc_fft_kb_now = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
       ekc_fft_kb_now = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
       kc_str = MAX(skc_fft_kb_now,skc_fft_ka_proc_all[iproc]);
       kc_end = MIN(ekc_fft_kb_now,ekc_fft_ka_proc_all[iproc]);
       i      = (sum_fft_kb_proc[ka] + kc_str-skc_fft_kb_now)*nkf2*2;
       for(kc=kc_str;kc<=kc_end;kc++){
         skb_fft_ka_now = (kc==skc_fft_ka_proc_all[iproc] ? 
                               skb_fft_ka_proc_all[iproc] : 1);
         ekb_fft_ka_now = (kc==ekc_fft_ka_proc_all[iproc] ? 
                               ekb_fft_ka_proc_all[iproc] : nkf2);
         kb_str = MAX(1,skb_fft_ka_now);
         kb_end = MIN(nkf2,ekb_fft_ka_now);
         ioff   = i+2*(kb_str-1);
         for(kb=1;kb<=2*(kb_end-kb_str+1);kb++){
           zfft[(j+kb)] = zfft_tmp[(ioff+kb)];
         }/*endfor*/
         iii = 2*MAX(kb_end-kb_str+1,0);
         j+=iii;i+= 2*nkf2;
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/
#endif

    /*-----------------------------*/
    /* Communicate Blocks          */
    /*-----------------------------*/

#ifdef DALLTOALL
    Alltoallv(&zfft[1],&(sendcounts_fft_ka[1]),&(senddspls_fft_ka[1]),
              MPI_DOUBLE,&zfft_tmp[1],&(recvcounts_fft_ka[1]),
              &(recvdspls_fft_ka[1]),MPI_DOUBLE,comm);
#endif

 }/*endif*/

  /*-----------------------------*/
  /* Index Transpose plus spread */
  /*-----------------------------*/

   ndata = 2*nfft_ka_proc*nkf1;
   for(i=1;i<=ndata;i++){zfft[i]=0.0;}

#ifdef TRANS
   i = 0;
   for(ka=1;ka<=nka_fft_kc;ka++){
    ioff = 2*ka_fft_kc_red[ka]-1;
    for(k=1,m=0;k<=2*nfft_ka_proc;k+=2,m+=2*nkf1){
      zfft[(m+ioff)]   = zfft_tmp[(i+k)];
      zfft[(m+ioff+1)] = zfft_tmp[(i+k+1)];
    }/*endfor*/
    i+=2*nfft_ka_proc;
   }/*endfor*/
#endif
/*==========================================================================*/
/* Perform FFT along ``a'' using zfft                                       */

  incn = nkf1;
#ifdef FFT
  fft_gen1d(zfft,nkf1,incl,nfft_ka_proc,incn,iopt,&ier,
            work_1a_f,nwork1,work_2a_f,nwork2,ifax_a_f,igeneric_opt);
#endif

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sngl_pack_coef(int ncoef_proc,int ncoef_use,int ndata,int num_proc,
                    double *cre,double *cim,double *zfft,
                    int *map_proc,int *map_c_proc)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

   int i;

/*=========================================================================*/
/* Pack the data up: top and bottom half of k-space : zero fill in scalar */

  if(num_proc==1){for(i=1;i<=ndata;i++){zfft[i]=0.0;} }

  for(i=1;i<=ncoef_use;i++){
    zfft[(map_proc[i])]     =  cre[i];
    zfft[(map_proc[i]+1)]   =  cim[i];
    zfft[(map_c_proc[i])]   =  cre[i]; 
    zfft[(map_c_proc[i]+1)] = -cim[i];
  }/*endfor*/

  if(ncoef_proc > ncoef_use){
    i = ncoef_proc;
    zfft[(map_proc[i])]     = cre[i];
    zfft[(map_proc[i]+1)]   = cim[i];
  }/*endif*/


/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dble_pack_coef(int ncoef_proc,int ncoef_use,int ndata,int num_proc,
               double *c1re, double *c1im,
               double *c2re, double *c2im,
               double *zfft,int *map_proc,int *map_c_proc)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
   int i;

/*=========================================================================*/
/* Pack the data up: top and bottom half of k-space : zero fill in scalar */

  if(num_proc==1){for(i=1;i<=ndata;i++){zfft[i]=0.0;} }

  for(i=1;i<=ncoef_use;i++){
    zfft[(map_proc[i])]     =   c1re[i] - c2im[i];
    zfft[(map_proc[i]+1)]   =   c1im[i] + c2re[i];
    zfft[(map_c_proc[i])]   =   c1re[i] + c2im[i]; 
    zfft[(map_c_proc[i]+1)] =  -c1im[i] + c2re[i];
  }/*endfor*/

  if(ncoef_proc > ncoef_use){
    i = ncoef_proc;
    zfft[(map_proc[i])]     = c1re[i];
    zfft[(map_proc[i]+1)]   = c2re[i];
  }/*endif*/


/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void para_fft_gen3d_bck_to_g(double *zfft, double *zfft_tmp,
                             PARA_FFT_PKG3D *para_fft_pkg3d)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int iopt = -1,incl=1;
 int ier  = 0;
 int incn;
 int i,j,k,m,ka,kb,kc,ioff,joff,koff,iii;
 int kap,kbp,kcp,kstr,kend,jb;
 int skc_fft_kb_now,ekc_fft_kb_now;
 int skb_fft_kb_now,ekb_fft_kb_now;
 int skb_fft_ka_now,ekb_fft_ka_now;
 int ndata,iproc;
 int kc_str,kc_end,kb_str,kb_end,ka_str,ka_end;
 double scale;
 int nfft;

/*-----------------------------------------------------------------------*/
/*          Local pointers                                               */

/* Communicator */
 MPI_Comm comm = para_fft_pkg3d->comm;
 int num_proc  = para_fft_pkg3d->num_proc;
 int myid      = para_fft_pkg3d->myid; 
 int myidp1    = para_fft_pkg3d->myid+1; 

/* General Info */
 int nktot        = para_fft_pkg3d->nktot;
 int ncoef        = para_fft_pkg3d->ncoef;
 int ncoef_proc   = para_fft_pkg3d->ncoef_proc;
 int ncoef_use    = para_fft_pkg3d->ncoef_use;
 int icoef_off    = para_fft_pkg3d->icoef_off;
 int icoef_strt   = para_fft_pkg3d->icoef_strt;
 int igeneric_opt = para_fft_pkg3d->igeneric_opt;
 int scale_opt    = para_fft_pkg3d->scale_opt; 
 int nwork1       = para_fft_pkg3d->nwork1;
 int nwork2       = para_fft_pkg3d->nwork2;

/* Map */
 int *map_proc_post = para_fft_pkg3d->map_proc_post; 

/* FFT along kc  */
 int nkf3                = para_fft_pkg3d->nkf3; 
 int *sendcounts_fft_kc  = para_fft_pkg3d->sendcounts_fft_kc;
 int *senddspls_fft_kc   = para_fft_pkg3d->senddspls_fft_kc;
 int *recvcounts_fft_kc  = para_fft_pkg3d->recvcounts_fft_kc;
 int *recvdspls_fft_kc   = para_fft_pkg3d->recvdspls_fft_kc;
 int nfft_kc_proc        = para_fft_pkg3d->nfft_kc_proc;
 int str_fft_kc_proc     = para_fft_pkg3d->str_fft_kc_proc;
 int end_fft_kc_proc     = para_fft_pkg3d->end_fft_kc_proc;
 int nka_fft_kc          = para_fft_pkg3d->nka_fft_kc;
 int *sum_fft_kc_proc    = para_fft_pkg3d->sum_fft_kc_proc;
 int *nkb_fft_kc         = para_fft_pkg3d->nkb_fft_kc;
 int *ka_fft_kc          = para_fft_pkg3d->ka_fft_kc;
 int *ka_fft_kc_red      = para_fft_pkg3d->ka_fft_kc_red;
 int *kb_fft_kc          = para_fft_pkg3d->kb_fft_kc;
 int ska_fft_kc_proc     = para_fft_pkg3d->ska_fft_kc_proc;
 int eka_fft_kc_proc     = para_fft_pkg3d->eka_fft_kc_proc;
 int skb_fft_kc_proc     = para_fft_pkg3d->skb_fft_kc_proc;
 int ekb_fft_kc_proc     = para_fft_pkg3d->ekb_fft_kc_proc;
 int *ifax_c_r           = para_fft_pkg3d->ifax_c_r;
 double *work_1c_r       = para_fft_pkg3d->work_1c_r;
 double *work_2c_r       = para_fft_pkg3d->work_2c_r;

/* FFT along kb  */
 int nkf2                 = para_fft_pkg3d->nkf2; 
 int *sendcounts_fft_kb   = para_fft_pkg3d->sendcounts_fft_kb;
 int *senddspls_fft_kb    = para_fft_pkg3d->senddspls_fft_kb;
 int *recvcounts_fft_kb   = para_fft_pkg3d->recvcounts_fft_kb;
 int *recvdspls_fft_kb    = para_fft_pkg3d->recvdspls_fft_kb;
 int nfft_kb_proc         = para_fft_pkg3d->nfft_kb_proc;
 int str_fft_kb_proc      = para_fft_pkg3d->str_fft_kb_proc;
 int end_fft_kb_proc      = para_fft_pkg3d->end_fft_kb_proc;
 int ska_fft_kb_proc      = para_fft_pkg3d->ska_fft_kb_proc;
 int eka_fft_kb_proc      = para_fft_pkg3d->eka_fft_kb_proc;
 int skc_fft_kb_proc      = para_fft_pkg3d->skc_fft_kb_proc;
 int ekc_fft_kb_proc      = para_fft_pkg3d->ekc_fft_kb_proc;
 int *ska_fft_kb_proc_all = para_fft_pkg3d->ska_fft_kb_proc_all;
 int *eka_fft_kb_proc_all = para_fft_pkg3d->eka_fft_kb_proc_all;
 int *skc_fft_kb_proc_all = para_fft_pkg3d->skc_fft_kb_proc_all;
 int *ekc_fft_kb_proc_all = para_fft_pkg3d->ekc_fft_kb_proc_all;
 int *sum_fft_kb_proc     = para_fft_pkg3d->sum_fft_kb_proc;
 int *ifax_b_r            = para_fft_pkg3d->ifax_b_r;
 double *work_1b_r        = para_fft_pkg3d->work_1b_r;
 double *work_2b_r        = para_fft_pkg3d->work_2b_r;

/* FFT along ka  */
 int nkf1                 = para_fft_pkg3d->nkf1; 
 int *sendcounts_fft_ka   = para_fft_pkg3d->sendcounts_fft_ka;
 int *senddspls_fft_ka    = para_fft_pkg3d->senddspls_fft_ka;
 int *recvcounts_fft_ka   = para_fft_pkg3d->recvcounts_fft_ka;
 int *recvdspls_fft_ka    = para_fft_pkg3d->recvdspls_fft_ka;
 int nfft_ka_proc         = para_fft_pkg3d->nfft_ka_proc;
 int str_fft_ka_proc      = para_fft_pkg3d->str_fft_ka_proc;
 int end_fft_ka_proc      = para_fft_pkg3d->end_fft_ka_proc;
 int skb_fft_ka_proc      = para_fft_pkg3d->skb_fft_ka_proc;
 int ekb_fft_ka_proc      = para_fft_pkg3d->ekb_fft_ka_proc;
 int skc_fft_ka_proc      = para_fft_pkg3d->skc_fft_ka_proc;
 int ekc_fft_ka_proc      = para_fft_pkg3d->ekc_fft_ka_proc;
 int *skb_fft_ka_proc_all = para_fft_pkg3d->skb_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all = para_fft_pkg3d->ekb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all = para_fft_pkg3d->skc_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all = para_fft_pkg3d->ekc_fft_ka_proc_all;
 int *ifax_a_r            = para_fft_pkg3d->ifax_a_r;
 double *work_1a_r        = para_fft_pkg3d->work_1a_r;
 double *work_2a_r        = para_fft_pkg3d->work_2a_r;

/*==========================================================================*/
/* Perform Inverse FFT along ``a'' using zfft                               */

  incn = nkf1;
#ifdef FFT
  fft_gen1d(zfft,nkf1,incl,nfft_ka_proc,incn,iopt,&ier,
            work_1a_r,nwork1,work_2a_r,nwork2,ifax_a_r,igeneric_opt);
#endif

/*==========================================================================*/
/*  Inverse Transpose III: Back to the coefs                                */

#ifdef TRANS
  /*---------------------------------------*/
  /* Index inverse transpose plus contract */
  /*---------------------------------------*/
   i = 0;
   for(ka=1;ka<=nka_fft_kc;ka++){
    ioff = 2*ka_fft_kc_red[ka]-1;
    for(k=1,m=0;k<=2*nfft_ka_proc;k+=2,m+=2*nkf1){
      zfft_tmp[(i+k)]   = zfft[(m+ioff)];
      zfft_tmp[(i+k+1)] = zfft[(m+ioff+1)];
    }/*endfor*/
    i+=2*nfft_ka_proc;
   }/*endfor*/
#else
   ndata = 2*nfft_ka_proc*nkf1;
   for(i=1;i<=ndata;i++){zfft_tmp[i]=0.0;}
#endif

 if(num_proc>1){

    /*-----------------------------*/
    /* Communicate Blocks          */
    /*-----------------------------*/

#ifdef DALLTOALL
    Alltoallv(&zfft_tmp[1],&(recvcounts_fft_ka[1]),&(recvdspls_fft_ka[1]),
              MPI_DOUBLE,&zfft[1],&(sendcounts_fft_ka[1]),
              &(senddspls_fft_ka[1]),MPI_DOUBLE,comm);
#endif

#ifdef BLOCK_PACK

   /*---------------------------------------------*/
   /* In place-block inverse transpose            */
   /*---------------------------------------------*/
   /*   i)Loop over ka's you have                 */
   /*  ii)Loop over procs=iproc                   */
   /* iii)Loop over the kc range to send to iproc */
   /*  iv)Loop over the kb range to send to iproc */
   /*   v)Store allowed/in range data            */
 
   j=0;ndata=0;
   for(iproc=1;iproc<=num_proc;iproc++){
     for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
       skc_fft_kb_now = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
       ekc_fft_kb_now = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
       kc_str = MAX(skc_fft_kb_now,skc_fft_ka_proc_all[iproc]);
       kc_end = MIN(ekc_fft_kb_now,ekc_fft_ka_proc_all[iproc]);
       i      = (sum_fft_kb_proc[ka] + kc_str-skc_fft_kb_now)*nkf2*2;
       for(kc=kc_str;kc<=kc_end;kc++){
         skb_fft_ka_now = (kc==skc_fft_ka_proc_all[iproc] ? 
                               skb_fft_ka_proc_all[iproc] : 1);
         ekb_fft_ka_now = (kc==ekc_fft_ka_proc_all[iproc] ? 
                               ekb_fft_ka_proc_all[iproc] : nkf2);
         kb_str = MAX(1,skb_fft_ka_now);
         kb_end = MIN(nkf2,ekb_fft_ka_now);
         ioff = i + 2*(kb_str-1);
         for(kb=1;kb<=2*(kb_end-kb_str+1);kb++){
           ndata++;
           zfft_tmp[(ioff+kb)] = zfft[(j+kb)];
         }/*endfor*/
         iii = 2*MAX(kb_end-kb_str+1,0);
         j+=iii;i+= 2*nkf2;
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/
#endif

 }/*endif*/

/*==========================================================================*/
/* Perform  Inverse FFT along ``b'' using zfft_tmp                          */

  incn = nkf2;
#ifdef FFT
  fft_gen1d(zfft_tmp,nkf2,incl,nfft_kb_proc,incn,iopt,&ier,
            work_1b_r,nwork1,work_2b_r,nwork2,ifax_b_r,igeneric_opt);
#endif

/*==========================================================================*/
/* Inverse Transpose II */

  /*---------------------------------------*/
  /* Inverse Index Transpose plus contract */
  /*---------------------------------------*/

#ifdef TRANS
   jb=0;for(i=1;i<=(ska_fft_kb_proc-1);i++){jb+=nkb_fft_kc[i];}
   i=0;koff = 0;
   for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
     kc_str = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
     kc_end = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
     for(kb=1;kb<=nkb_fft_kc[ka];kb++){
      jb++;ioff = 2*kb_fft_kc[jb]-1 + koff;
      for(kc=1,m=0;kc<=2*(kc_end-kc_str+1);kc+=2,m+=2*nkf2){
        zfft[(i+kc)]   = zfft_tmp[(m+ioff)];
        zfft[(i+kc+1)] = zfft_tmp[(m+ioff+1)];
      }/*endfor*/
      i+=2*(kc_end-kc_str+1);
    }/*endfor*/
     koff+= (kc_end-kc_str+1)*nkf2*2;
   }/*endfor*/
#else
   ndata = 2*nfft_kb_proc*nkf2;
   for(i=1;i<=ndata;i++){zfft[i]=0.0;}
#endif

   if(num_proc>1){

    /*-----------------------------*/
    /* Communicate Blocks          */
    /*-----------------------------*/

#ifdef SALLTOALL
    Alltoallv(&zfft[1],&(recvcounts_fft_kb[1]),&(recvdspls_fft_kb[1]),
              MPI_DOUBLE,&zfft_tmp[1],&(sendcounts_fft_kb[1]),
              &(senddspls_fft_kb[1]),MPI_DOUBLE,comm);
#endif

#ifdef BLOCK_PACK
    /*---------------------------------------------*/
    /* Inverse in  place-block transpose           */
    /*---------------------------------------------*/
    /*   i)Loop over procs=iproc                   */
    /*  ii)Loop over the ka range to send to iproc */
    /* iii)Loop over the kb range you have         */
    /*  iv)Loop over the kc range to send to iproc */
    /*   v)Store allowed/in range data             */

    j=0;
    for(iproc=1;iproc<=num_proc;iproc++){
     ka_str = MAX(ska_fft_kc_proc,ska_fft_kb_proc_all[iproc]);
     ka_end = MIN(eka_fft_kc_proc,eka_fft_kb_proc_all[iproc]);
     for(ka=ka_str;ka<=ka_end;ka++){
       kb_str = (ka==ska_fft_kc_proc ? skb_fft_kc_proc : 1);
       kb_end = (ka==eka_fft_kc_proc ? ekb_fft_kc_proc : nkb_fft_kc[ka]);
       i = sum_fft_kc_proc[ka]*nkf3*2;
       kc_str = (ka==ska_fft_kb_proc_all[iproc] ? 
                     skc_fft_kb_proc_all[iproc]  : 1);
       kc_end = (ka==eka_fft_kb_proc_all[iproc] ?  
                     ekc_fft_kb_proc_all[iproc]  : nkf3);
       for(kb=kb_str;kb<=kb_end;kb++){
         ioff = i + 2*(kc_str-1);
         for(kc=1;kc<=2*(kc_end-kc_str+1);kc++){
           zfft[(ioff+kc)] = zfft_tmp[(j+kc)];
         }/*endfor*/
         iii = 2*MAX(kc_end-kc_str+1,0);
         j+=iii;i+= 2*nkf3;
       }/*endfor*/
     }/*endfor*/
    }/*endfor*/
#endif

   }/*endif:num_proc > 1*/

/*==========================================================================*/
/* Perform Inverse FFT along ``c'' using zfft                              */

  incn = nkf3;
#ifdef FFT
  fft_gen1d(zfft,nkf3,incl,nfft_kc_proc,incn,iopt,&ier,
            work_1c_r,nwork1,work_2c_r,nwork2,ifax_c_r,igeneric_opt);
#endif

/*==========================================================================*/
/*  Inverse Tranpose I:                                                     */

  if(num_proc>1){

    /*---------------------------------------*/
    /* Inverse Index Transpose plus contract */
    /*---------------------------------------*/
#ifdef POST_MAP
    ndata = recvdspls_fft_kc[num_proc]+recvcounts_fft_kc[num_proc];
    for(i=1;i<=ndata;i++){
      zfft_tmp[i]=zfft[map_proc_post[i]];
    }/*endfor*/
#endif

    /*-----------------------------*/
    /* Communicate blocks          */
    /*-----------------------------*/

#ifdef SALLTOALL
    Alltoallv(&zfft_tmp[1],&(recvcounts_fft_kc[1]),&(recvdspls_fft_kc[1]),
              MPI_DOUBLE,&zfft[1],&(sendcounts_fft_kc[1]),
              &(senddspls_fft_kc[1]),MPI_DOUBLE,comm);
#endif

  }/*endif*/
   
/*==========================================================================*/
/*   IV) scale zfft if necessary                                           */

  if(scale_opt == 1){
    ndata = 2*nfft_kc_proc*nkf3;
    nfft  = nkf1*nkf2*nkf3; scale = 1.0/((double) nfft);
    for(i=1;i<=ndata;i++){zfft[i] *= scale;}
  }/* endif */

  if(scale_opt == -1){
    ndata =  2*nfft_kc_proc*nkf3;
    nfft  = nkf1*nkf2*nkf3; scale = ((double) nfft);
    for(i=1;i<=ndata;i++){zfft[i] *= scale;}
  }/* endif */

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sngl_upack_coef(int ncoef_proc,double *cre,double *cim,double *zfft,
               int *map_proc)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int i;

/*=======================================================================*/
/*  Unpack the data : Top half of k space only */

  for(i=1;i<=ncoef_proc;i++){
    cre[i]=zfft[map_proc[i]];
    cim[i]=zfft[(map_proc[i]+1)];
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sngl_upack_rho(int ndata,double *zfft,double *rfft)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int m,i;

/*=======================================================================*/
/*  Unpack the data : Top half of k space only */

  for(i=1,m=1;i<=ndata;i++,m+=2){
    rfft[i] = zfft[m];
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sngl_pack_rho(int ndata,double *zfft,double *rfft)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int m,i;

/*=======================================================================*/
/*  Unpack the data : Top half of k space only */

  for(i=1,m=1;i<=ndata;i++,m+=2){
    zfft[m] = rfft[i];
    zfft[(m+1)] = 0.0;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dble_upack_rho(int ndata,double *zfft,double *rfft1,double *rfft2)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int i;

/*=======================================================================*/
/*  Unpack the data : */

   for(i=1;i<=ndata;i++){rfft1[i] = zfft[i];}
   for(i=1;i<=ndata;i++){rfft2[i] = zfft[(i+ndata)];}

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dble_pack_rho(int ndata,double *zfft,double *rfft1,double *rfft2)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int m,i;

/*=======================================================================*/
/*  Unpack the data : Top half of k space only */

  for(i=1,m=1;i<=ndata;i++,m+=2){
    zfft[m] = rfft1[i];
    zfft[(m+1)] = rfft2[i];
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dble_upack_coef(int ncoef_proc,
                     double *c1re,double *c1im,double *c2re,double *c2im,
                     double *zfft,int *map_proc,int *map_c_proc)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int i;
double tempr ,tempi;
double temprc,tempic;

/*=======================================================================*/
/*  Unpack the data :  */

  for(i=1;i<=ncoef_proc;i++){
 
    tempr  = zfft[(map_proc[i])];
    tempi  = zfft[(map_proc[i]+1)];
    temprc = zfft[(map_c_proc[i])];
    tempic = zfft[(map_c_proc[i]+1)];

    c2im[i] = 0.5*(-tempr + temprc);
    c1re[i] = 0.5*( tempr + temprc);

    c1im[i] = 0.5*( tempi - tempic);
    c2re[i] = 0.5*( tempi + tempic);

  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sum_rho(int ndata,double *zfft,double *rfft)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int m,i;

/*=======================================================================*/
/*  Unpack the data : Top half of k space only */

  for(i=1,m=1;i<=ndata;i++,m+=2){
    rfft[i] += zfft[m]*zfft[m] + zfft[(m+1)]*zfft[(m+1)];
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void create_para_fft_pkg3d(PARA_FFT_PKG3D *para_fft_pkg3d,
                           int *kastr,int *kbstr,int *kcstr)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

  int i,ic,ic_c,i_00,i_0,ic_00,ic_0,m;
  int ioff,ifft,iii;
  int nfft_kc_half,nfft_kc_half_c;
  int ka,kb,kc,kstr,kend;
  int ka_old,kb_old,kc_old;
  int ka_old_c,kb_old_c,kc_old_c;
  int ka_str,ka_end,kb_str,kb_end,kc_str,kc_end;
  int skc_fft_kb_now,ekc_fft_kb_now;
  int skc_fft_ka_now,ekc_fft_ka_now;
  int skb_fft_ka_now,ekb_fft_ka_now;
  int kcp,kbp,kap;
  int idiv,irem,jjj,kkk,j;
  int iproc,iprocm1,isum,nfft_size;
  int ncoef2_proc;
  int *recv_counts_rho,*displs_rho,nfft2_proc;
  int *recv_counts_coef;

/*-----------------------------------------------------------------------*/
/*          Local pointers                                               */

  /* Pack maps */
  int *map,*map_c,*map_inv,*map_c_inv;
  int *indx,*indx_c;
  int *map_proc,*map_c_proc,*map_proc_post;

  /* FFT along kc  */
   int nfft_kc,nfft_kc_proc; 
   int str_fft_kc_proc,end_fft_kc_proc;
   int ska_fft_kc_proc,eka_fft_kc_proc;
   int skb_fft_kc_proc,ekb_fft_kc_proc;
   int nka_fft_kc, *nkb_fft_kc;
   int *ka_fft_kc,*kb_fft_kc,*ka_fft_kc_red;
   int *ifft_strt_kc,*ifft_strt_kc_c;
   int *sum_fft_kc_proc;
   int *nkb_fft_kc_low;
   int *nkb_fft_kc_high;
   int *nfft_kc_proc_all;
   int *str_fft_kc_proc_all,*end_fft_kc_proc_all;
   int *sendcounts_fft_kc,*senddspls_fft_kc;
   int *recvcounts_fft_kc,*recvdspls_fft_kc;


  /* FFT along kb  */
   int nfft_kb,nfft_kb_proc; 
   int str_fft_kb_proc,end_fft_kb_proc;
   int ska_fft_kb_proc,eka_fft_kb_proc;
   int skc_fft_kb_proc,ekc_fft_kb_proc;
   int *sum_fft_kb_proc;
   int *nfft_kb_proc_all;
   int *ska_fft_kb_proc_all,*eka_fft_kb_proc_all;
   int *skc_fft_kb_proc_all,*ekc_fft_kb_proc_all;
   int *str_fft_kb_proc_all,*end_fft_kb_proc_all;
   int *sendcounts_fft_kb,*senddspls_fft_kb;
   int *recvcounts_fft_kb,*recvdspls_fft_kb;

  /* FFT along ka  */
   int nfft_ka,nfft_ka_proc; 
   int str_fft_ka_proc,end_fft_ka_proc;
   int skb_fft_ka_proc,ekb_fft_ka_proc;
   int skc_fft_ka_proc,ekc_fft_ka_proc;
   int *sum_fft_ka_proc;
   int *nfft_ka_proc_all;
   int *skb_fft_ka_proc_all,*ekb_fft_ka_proc_all;
   int *skc_fft_ka_proc_all,*ekc_fft_ka_proc_all;
   int *str_fft_ka_proc_all,*end_fft_ka_proc_all;
   int *sendcounts_fft_ka,*senddspls_fft_ka;
   int *recvcounts_fft_ka,*recvdspls_fft_ka;

/* Local pointers  */
 MPI_Comm comm = para_fft_pkg3d->comm;
 int num_proc  = para_fft_pkg3d->num_proc;
 int myid      = para_fft_pkg3d->myid; 
 int myidp1    = para_fft_pkg3d->myid+1; 

 int nktot      = para_fft_pkg3d->nktot; 
 int ncoef      = para_fft_pkg3d->ncoef;
 int ncoef_proc = para_fft_pkg3d->ncoef_proc;
 int ncoef_use  = para_fft_pkg3d->ncoef_use;
 int icoef_off  = para_fft_pkg3d->icoef_off;
 int icoef_strt = para_fft_pkg3d->icoef_strt;
 int nkf3       = para_fft_pkg3d->nkf3; 
 int nkf2       = para_fft_pkg3d->nkf2; 
 int nkf1       = para_fft_pkg3d->nkf1; 

/*=======================================================================*/
/* I) Malloc some memory */

   map          = (int *) cmalloc(ncoef*sizeof(int))-1;
   map_c        = (int *) cmalloc(ncoef*sizeof(int))-1;
   map_inv      = (int *) cmalloc(4*ncoef*sizeof(int))-1;
   map_c_inv    = (int *) cmalloc(2*ncoef*sizeof(int))-1;
   indx         = (int *) cmalloc(2*ncoef*sizeof(int))-1;
   indx_c       = (int *) cmalloc(2*ncoef*sizeof(int))-1;

   map_proc     = (int *) cmalloc(ncoef*sizeof(int))-1;
   map_c_proc   = (int *) cmalloc(ncoef*sizeof(int))-1;
   map_proc_post= (int *) cmalloc(4*ncoef*sizeof(int))-1;

/*========================================================================*/
/* I) Find FFT indices, sort them to make map_inv and map_c_inv           */
/*      then create standard map which are easier to handle               */
/*      and include the c-FFT offset filling                              */

/*------------------------------------------------------------------------*/
/*  i) Find and sort indicies */

  setfft_indx(nkf1,nkf2,nkf3,ncoef,kastr,kbstr,kcstr,
              indx,indx_c);
  for(i=1;i<=ncoef;i++){map_inv[i]=i;map_c_inv[i]=i;}
  sort_commence(ncoef,indx,map_inv);
  sort_commence(nktot,indx_c,map_c_inv);

/*------------------------------------------------------------------------*/
/* ii) Find the starting indicies of FFT's along kc                       */
/*     and store them                                                     */
   
  nfft_kc_half   = 0;
  nfft_kc_half_c = 0;
  ka_old   = kastr[map_inv[1]]-1;  
  kb_old   = kbstr[map_inv[1]]-1;
  ka_old_c = kastr[map_c_inv[1]]-1;
  kb_old_c = kbstr[map_c_inv[1]]-1;
  for(i=1;i<=nktot;i++){
    if(kastr[map_inv[i]]!=ka_old||kbstr[map_inv[i]]!=kb_old){
      nfft_kc_half++;
    }/*endif*/
    if(kastr[map_c_inv[i]]!=ka_old_c||kbstr[map_c_inv[i]]!=kb_old_c){
      nfft_kc_half_c++;
    }/*endif*/
    ka_old   = kastr[map_inv[i]];  
    kb_old   = kbstr[map_inv[i]];
    ka_old_c = kastr[map_c_inv[i]];
    kb_old_c = kbstr[map_c_inv[i]];
  }/*endfor*/
  i = ncoef;
  if(kastr[map_inv[i]]!=ka_old||kbstr[map_inv[i]]!=kb_old){
      nfft_kc_half++;
  }/*endif*/

  ifft_strt_kc     = (int *) cmalloc((nfft_kc_half+1)*sizeof(int))-1;
  ifft_strt_kc_c   = (int *) cmalloc((nfft_kc_half_c+1)*sizeof(int))-1;
  ka_fft_kc_red    = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;
  ka_fft_kc        = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;
  kb_fft_kc        = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;
  nkb_fft_kc       = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;
  nkb_fft_kc_low   = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;
  nkb_fft_kc_high  = (int *) cmalloc((2*nfft_kc_half-1)*sizeof(int))-1;

  ic   = 0;
  ic_c = 0;
  ka_old   = kastr[map_inv[1]]-1;  
  kb_old   = kbstr[map_inv[1]]-1;
  ka_old_c = kastr[map_c_inv[1]]-1;
  kb_old_c = kbstr[map_c_inv[1]]-1;
  for(i=1;i<=nktot;i++){
    if(kastr[map_inv[i]]!=ka_old||kbstr[map_inv[i]]!=kb_old){
     ic++;ifft_strt_kc[ic] = i;
    }/*endif*/
    if(kastr[map_c_inv[i]]!=ka_old_c||kbstr[map_c_inv[i]]!=kb_old_c){
     ic_c++; ifft_strt_kc_c[ic_c] = i;
    }/*endif*/
    ka_old   = kastr[map_inv[i]];  
    kb_old   = kbstr[map_inv[i]];
    ka_old_c = kastr[map_c_inv[i]];
    kb_old_c = kbstr[map_c_inv[i]];
  }/*endfor*/
  i = ncoef;
  if(kastr[map_inv[i]]!=ka_old||kbstr[map_inv[i]]!=kb_old){
     ic++;ifft_strt_kc[ic] = i;
  }/*endif*/
  ifft_strt_kc[(ic+1)]     = ncoef+1;
  ifft_strt_kc_c[(ic_c+1)] = ncoef+1;

/*------------------------------------------------------------------------*/
/*  iii) Find end points of FFTs with zero for kx or kx and ky            */

  for(i=1;i<=nktot;i++){
    if(kastr[map_inv[i]]==0  &&kbstr[map_inv[i]]==0)  {i_00=i;}
    if(kastr[map_inv[i]]==0)  {i_0=i;}
    if(kastr[map_c_inv[i]]==0&&kbstr[map_c_inv[i]]==0){ic_00=i;}
    if(kastr[map_c_inv[i]]==0){ic_0=i;}
  }/*endfor*/
  i = ncoef;
  if(kastr[map_inv[i]]==0  &&kbstr[map_inv[i]]==0)  {i_00=i;}
  if(kastr[map_inv[i]]==0)  {i_0=i;}

/*------------------------------------------------------------------------*/
/*  iv) Create maps with FFT offset on c                                 */

  ioff = 0; ic = 2; ic_c = 2; ifft=1; nka_fft_kc=1;
  for(i=1;i<=2*nfft_kc_half-1;i++){nkb_fft_kc[i]=0;}
  for(i=1;i<=2*nfft_kc_half-1;i++){nkb_fft_kc_low[i]=0;}
  for(i=1;i<=2*nfft_kc_half-1;i++){nkb_fft_kc_high[i]=0;}

  /*-----------------------------------------------------------------------*/
  /*  a)FFT with 0,0,kz -- negative kc values get mapped to FFT order        */

  ka_old          =0;
  ka_fft_kc_red[1]=1;
  ka_fft_kc[1]    =1;
  kb_fft_kc[1]    =1;
  nkb_fft_kc[1]++;
  nkb_fft_kc_low[1]++;
  for(i=1;      i<=i_00; i++     ){
    kcp = (kcstr[map_inv[i]] < 0 ? kcstr[map_inv[i]]+nkf3+1 : 
                                      kcstr[map_inv[i]]+1);
    map[map_inv[i]]     = kcp;       
  }/*endfor*/
  for(i=1;      i<=ic_00;i++     ){
    kcp = (-kcstr[map_c_inv[i]] < 0 ? -kcstr[map_c_inv[i]]+nkf3+1 : 
                                         -kcstr[map_c_inv[i]]+1);
    map_c[map_c_inv[i]] = kcp;       
  }/*endfor*/

  /*--------------------------------------------*/
  /* b)FFTs with  0,ky,kz -- FFT packing order */

  for(i=i_00+1; i<=i_0;  i++     ){
    if(ifft_strt_kc[ic]==i){
      nkb_fft_kc[1]++;
      nkb_fft_kc_low[1]++;
      ic++;ioff+=nkf3;ifft++;
      ka_fft_kc[ifft]=1;
      kb_fft_kc[ifft]=kbstr[map_inv[i]]+1;
    }/*endif*/
    kcp = (kcstr[map_inv[i]] < 0 ? kcstr[map_inv[i]]+nkf3+1 : 
                                      kcstr[map_inv[i]]+1);
    map[map_inv[i]]     = kcp+ioff; 
  }/*endfor*/
  for(i=ic_00+1;i<=ic_0; i++     ){
    if(ifft_strt_kc_c[ic_c]==i){
      nkb_fft_kc[1]++;
      nkb_fft_kc_high[1]++;
      ic_c++;ioff+=nkf3;ifft++;
      ka_fft_kc[ifft]=1;
      kb_fft_kc[ifft]=-kbstr[map_c_inv[i]]+nkf2+1;
    }/*endif*/
    kcp = (-kcstr[map_c_inv[i]] < 0 ? -kcstr[map_c_inv[i]]+nkf3+1 : 
                                         -kcstr[map_c_inv[i]]+1);
    map_c[map_c_inv[i]] = kcp+ioff; 
  }/*endfor*/

  /*----------------------------------------*/
  /*   c)FFTs kx,ky,kz -- FFT packing order */
  for(i=i_0+1;  i<=ncoef;i++     ){
    if(ifft_strt_kc[ic]==i){
      ic++;ioff+=nkf3;ifft++;
      if(ka_old!=kastr[map_inv[i]]){
        nka_fft_kc++;
        ka_old             = kastr[map_inv[i]];
        ka_fft_kc_red[nka_fft_kc] = kastr[map_inv[i]]+1;
      }/*endif*/
      nkb_fft_kc[nka_fft_kc]++;
      ka_fft_kc[ifft]=kastr[map_inv[i]]+1;
      kbp = (kbstr[map_inv[i]] < 0 ? kbstr[map_inv[i]]+nkf2+1 :
                                        kbstr[map_inv[i]]+1);
      kb_fft_kc[ifft]=kbp;
      if(kbp<=nkf2/2+1){nkb_fft_kc_low[nka_fft_kc]++;}
      if(kbp>nkf2/2+1){nkb_fft_kc_high[nka_fft_kc]++;}
    }/*endif*/
    kcp = (kcstr[map_inv[i]] < 0 ? kcstr[map_inv[i]]+nkf3+1 : 
                                      kcstr[map_inv[i]]+1);
    map[map_inv[i]]     = kcp+ioff; 
  }/*endfor*/
  for(i=ic_0+1;  i<=nktot;i++     ){
    if(ifft_strt_kc_c[ic_c]==i){
     ic_c++;ioff+=nkf3;ifft++;
     if(ka_old!=-kastr[map_c_inv[i]]){
       nka_fft_kc++;
       ka_old             = -kastr[map_c_inv[i]];
       ka_fft_kc_red[nka_fft_kc] = -kastr[map_c_inv[i]]+nkf1+1;
     }/*endif*/
     nkb_fft_kc[nka_fft_kc]++;
     ka_fft_kc[ifft]=-kastr[map_c_inv[i]]+nkf1+1;
     kbp = (-kbstr[map_c_inv[i]] < 0 ? -kbstr[map_c_inv[i]]+nkf2+1 :
                                          -kbstr[map_c_inv[i]]+1);
     kb_fft_kc[ifft]=kbp;
     if(kbp<=nkf2/2+1){nkb_fft_kc_low[nka_fft_kc]++;}
     if(kbp>nkf2/2+1){nkb_fft_kc_high[nka_fft_kc]++;}
    }/*endif*/
    kcp = (-kcstr[map_c_inv[i]] < 0 ? -kcstr[map_c_inv[i]]+nkf3+1 : 
                                         -kcstr[map_c_inv[i]]+1);
    map_c[map_c_inv[i]] = kcp+ioff; 
  }/*endfor*/

/*==========================================================================*/
/* II) Transpose I maps: Create maps for on-proc shuffle                    */
/*                       Find displs and send/recvs for alltoall            */
/*                       Complicated in parallel due to half space condition*/

/*--------------------------------------------------------------------------*/
/*  i) Create proc maps for in place rearrangement                          */

  if(num_proc>1){
    for(i=1;i<=ncoef_proc+ncoef_use;i++){map_inv[i]   = i;}
    for(i=1;i<=ncoef_proc;i++){indx[i]                = map[i+icoef_off];}
    for(i=1;i<=ncoef_use ;i++){indx[(i+ncoef_proc)]   = map_c[i+icoef_off];}
    iii = (ncoef_proc+ncoef_use);
    sort_commence(iii,indx,map_inv);
    for(i=1;i<=ncoef_proc+ncoef_use;i++){
     if(map_inv[i]<=ncoef_proc){map_proc[map_inv[i]]   = i;}
     if(map_inv[i]> ncoef_proc){map_c_proc[map_inv[i]-ncoef_proc] = i;}
    }/*endfor*/
  }else{
    for(i=1;i<=ncoef;i++){map_proc[i]   = map[i];}
    for(i=1;i<=nktot;i++){map_c_proc[i] = map_c[i];}
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* ii) Determine the number of kc FFT's done on each proc                */

  nfft_kc           = (2*nfft_kc_half-1);
  idiv              = nfft_kc / num_proc;
  irem              = nfft_kc % num_proc;
  nfft_kc_proc      = (myid < irem ? idiv+1 : idiv);
  str_fft_kc_proc   = (myid <= irem ? myid*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(myid-irem)*idiv);
  end_fft_kc_proc   = str_fft_kc_proc + nfft_kc_proc - 1;

  nfft_kc_proc_all                    = (int *) cmalloc(num_proc*sizeof(int))-1;
  str_fft_kc_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  end_fft_kc_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  sum_fft_kc_proc     = (int *) cmalloc(nkf1*sizeof(int))-1;

  for(iproc=1;iproc<=num_proc;iproc++){
    iprocm1 = iproc-1;
    nfft_kc_proc_all[iproc]     = (iprocm1 < irem ? idiv+1 : idiv);
    str_fft_kc_proc_all[iproc]  = (iprocm1 <= irem ? iprocm1*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(iprocm1-irem)*idiv);
    end_fft_kc_proc_all[iproc]  = str_fft_kc_proc_all[iproc]
                                + nfft_kc_proc_all[iproc]-1;
  }/*endfor*/

/*------------------------------------------------------------------------*/
/* iii) Determine the ka and kb ranges for this processor depending on    */
/*where the kc FFTs start and end.  Part of spherical truncation condition */

  isum = 0; 
  for(ka=1;ka<=nka_fft_kc;ka++){
   for(kb=1;kb<=nkb_fft_kc[ka];kb++){
     isum++;
     if(isum==str_fft_kc_proc){ska_fft_kc_proc=ka;skb_fft_kc_proc=kb;}
     if(isum==end_fft_kc_proc){eka_fft_kc_proc=ka;ekb_fft_kc_proc=kb;}
   }/*endfor*/
  }/*endfor*/

/*-----------------------------------------------------------------------*/
/* iv) Get total number of kc FFTs on this processor for each allowed ka */
/*     by spherical truncation                                           */

  isum = 0;
  for(ka=ska_fft_kc_proc;ka<=eka_fft_kc_proc;ka++){
   sum_fft_kc_proc[ka] = isum;
   kstr = (ka==ska_fft_kc_proc ? skb_fft_kc_proc : 1);
   kend = (ka==eka_fft_kc_proc ? ekb_fft_kc_proc : nkb_fft_kc[ka]);
   isum += (kend-kstr+1);
  }/*endfor*/

/*----------------------------------------------------------------------*/
/* v) In place-block transpose/Alltoallv info                           */

  if(num_proc>1){

    sendcounts_fft_kc = (int *) cmalloc(num_proc*sizeof(int))-1;
    recvdspls_fft_kc  = (int *) cmalloc(num_proc*sizeof(int))-1; 
    senddspls_fft_kc  = (int *) cmalloc(num_proc*sizeof(int))-1;
    recvcounts_fft_kc = (int *) cmalloc(num_proc*sizeof(int))-1;

    for(i=1;i<=num_proc;i++){sendcounts_fft_kc[i]=0;}
    iproc=1;
    for(i=1;i<=ncoef_proc+ncoef_use;i++){
      iii = (indx[i]-1)/nkf3+1;
      while(iii > end_fft_kc_proc_all[iproc]&&iproc<num_proc){iproc++;}
      sendcounts_fft_kc[iproc]++;
    }/*endfor*/
  
    Alltoall(&sendcounts_fft_kc[1],1,MPI_INT,
             &recvcounts_fft_kc[1],1,MPI_INT,comm); 

    recvdspls_fft_kc[1] = 0;
    senddspls_fft_kc[1] = 0;
    for(i=1;i<=num_proc-1;i++){
      recvdspls_fft_kc[(i+1)] = recvdspls_fft_kc[i]+recvcounts_fft_kc[i];
      senddspls_fft_kc[(i+1)] = senddspls_fft_kc[i] +sendcounts_fft_kc[i];
    }/*endfor*/

/*-----------------------------------------------------------------------*/
/* vi) Create post alltoall map to complete FFT ordering.                */
/*     kc is zero filled, kb and ka are dense                            */

    for(i=1;i<=iii;i++){map_proc_post[i]=0;}

    Alltoallv(&indx[1],&(sendcounts_fft_kc[1]),&(senddspls_fft_kc[1]),
              MPI_INT,&map_proc_post[1],&(recvcounts_fft_kc[1]),
              &(recvdspls_fft_kc[1]),MPI_INT,comm);

    iii = recvdspls_fft_kc[num_proc]+recvcounts_fft_kc[num_proc];
    jjj = (str_fft_kc_proc-1)*nkf3;
    for(i=1;i<=iii;i++){
      map_proc_post[i]-=jjj;
      if(map_proc_post[i] <=0 || map_proc_post[i] > nfft_kc_proc*nkf3){
        printf("myid %d size %d  map %d\n",myid,nfft_kc_proc*nkf3,
                map_proc_post[i]);
      }/* endif */
    }/* endfor */

  }/*endif:num_proc>1*/


/*==========================================================================*/
/* III) Transpose II  :                                                     */

/*---------------------------------------------------------------------------*/
/* i) Determine the number of kb FFT's done on each proc                     */

  nfft_kb           = nka_fft_kc*nkf3;
  idiv              = nfft_kb / num_proc;
  irem              = nfft_kb % num_proc;
  nfft_kb_proc      = (myid < irem ? idiv+1 : idiv);
  str_fft_kb_proc   = (myid <= irem ? myid*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(myid-irem)*idiv);
  end_fft_kb_proc   = str_fft_kb_proc + nfft_kb_proc - 1;

  nfft_kb_proc_all    = (int *) cmalloc(num_proc*sizeof(int))-1;
  str_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  end_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  ska_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  skc_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  eka_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  ekc_fft_kb_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  sum_fft_kb_proc     = (int *) cmalloc(nkf1*sizeof(int))-1;

  for(iproc=1;iproc<=num_proc;iproc++){
      iprocm1 = iproc-1;
      nfft_kb_proc_all[iproc]    = (iprocm1 < irem ? idiv+1 : idiv);
      str_fft_kb_proc_all[iproc] = (iprocm1 <= irem ? iprocm1*(idiv+1)+1 : 
                                    irem*(idiv+1)+1+(iprocm1-irem)*idiv);
      end_fft_kb_proc_all[iproc] = str_fft_kb_proc_all[iproc] 
                                 + nfft_kb_proc_all[iproc] - 1;
  }/*endfor*/

/*----------------------------------------------------------------------------*/
/* ii) Determine kc and ka ranges for this proc depending on FFT's to be done */
/*     on this processor */

  isum = 0;
  for(ka=1;ka<=nka_fft_kc;ka++){
    for(kc=1;kc<=nkf3;kc++){
      isum++;
      for(iproc=1;iproc<=num_proc;iproc++){
        if(isum==str_fft_kb_proc_all[iproc]){
          ska_fft_kb_proc_all[iproc]=ka;
          skc_fft_kb_proc_all[iproc]=kc;
        }/*endif*/
        if(isum==end_fft_kb_proc_all[iproc]){
          eka_fft_kb_proc_all[iproc]=ka;
          ekc_fft_kb_proc_all[iproc]=kc;
        }/*endif*/
      }/*endfor*/
    }/*endfor*/
  }/*endfor*/

  ska_fft_kb_proc = ska_fft_kb_proc_all[myidp1];
  skc_fft_kb_proc = skc_fft_kb_proc_all[myidp1];
  eka_fft_kb_proc = eka_fft_kb_proc_all[myidp1];
  ekc_fft_kb_proc = ekc_fft_kb_proc_all[myidp1];
  
  isum = 0;  
  for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
    sum_fft_kb_proc[ka] = isum;
    kstr = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
    kend = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
    for(kc=kstr;kc<=kend;kc++){isum++;}
  }/*endfor*/


/*----------------------------------------------------------------------------*/
/* iii) In place-block transpose                                              */

  if(num_proc>1){

   /*---------------------------------------------*/
   /*   i)Loop over procs=iproc                   */
   /*  ii)Loop over the ka range to send to iproc */
   /* iii)Loop over the kb range you have         */
   /*  iv)Loop over the kc range to send to iproc */
   /*   v)Count                                   */

    sendcounts_fft_kb = (int *) cmalloc(num_proc*sizeof(int))-1;
    recvdspls_fft_kb  = (int *) cmalloc(num_proc*sizeof(int))-1; 
    senddspls_fft_kb  = (int *) cmalloc(num_proc*sizeof(int))-1;
    recvcounts_fft_kb = (int *) cmalloc(num_proc*sizeof(int))-1;

    for(i=1;i<=num_proc;i++){sendcounts_fft_kb[i]=0;}
    for(iproc=1;iproc<=num_proc;iproc++){
      ka_str = MAX(ska_fft_kc_proc,ska_fft_kb_proc_all[iproc]);
      ka_end = MIN(eka_fft_kc_proc,eka_fft_kb_proc_all[iproc]);
      for(ka=ka_str;ka<=ka_end;ka++){
        kb_str = (ka==ska_fft_kc_proc ? skb_fft_kc_proc : 1);
        kb_end = (ka==eka_fft_kc_proc ? ekb_fft_kc_proc : nkb_fft_kc[ka]);
        kc_str = (ka==ska_fft_kb_proc_all[iproc] ? 
                      skc_fft_kb_proc_all[iproc]  : 1);
        kc_end = (ka==eka_fft_kb_proc_all[iproc] ?  
                      ekc_fft_kb_proc_all[iproc]  : nkf3);
        for(kb=kb_str;kb<=kb_end;kb++){
          for(kc=kc_str;kc<=kc_end;kc++){sendcounts_fft_kb[iproc]++;}
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/

    Alltoall(&sendcounts_fft_kb[1],1,MPI_INT,
             &recvcounts_fft_kb[1],1,MPI_INT,comm); 

    recvdspls_fft_kb[1] = 0;
    senddspls_fft_kb[1] = 0;
    for(i=1;i<=num_proc-1;i++){
      recvdspls_fft_kb[(i+1)] = recvdspls_fft_kb[i]+recvcounts_fft_kb[i];
      senddspls_fft_kb[(i+1)] = senddspls_fft_kb[i] +sendcounts_fft_kb[i];
    }/*endfor*/

 }/*endif*/


/*==========================================================================*/
/* III) Transpose III                                                       */


/*--------------------------------------------------------------------------*/
/* i) Determine the number of ka FFT's done on each proc                    */

  nfft_ka           = nkf3*nkf2;
  idiv              = nfft_ka / num_proc;
  irem              = nfft_ka % num_proc;
  nfft_ka_proc      = (myid < irem ? idiv+1 : idiv);
  str_fft_ka_proc   = (myid <= irem ? myid*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(myid-irem)*idiv);
  end_fft_ka_proc   = str_fft_ka_proc + nfft_ka_proc - 1;


  nfft_ka_proc_all    = (int *) cmalloc(num_proc*sizeof(int))-1;
  str_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  end_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  skb_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  skc_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  ekb_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;
  ekc_fft_ka_proc_all = (int *) cmalloc(num_proc*sizeof(int))-1;

  for(iproc=1;iproc<=num_proc;iproc++){
    iprocm1 = iproc-1;
    nfft_ka_proc_all[iproc]    = (iprocm1 < irem ? idiv+1 : idiv);
    str_fft_ka_proc_all[iproc] = (iprocm1 <= irem ? iprocm1*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(iprocm1-irem)*idiv);
    end_fft_ka_proc_all[iproc] = str_fft_ka_proc_all[iproc] 
                               + nfft_ka_proc_all[iproc] - 1;
  }/*endfor*/

/*---------------------------------------------------------------------------*/
/* ii) Determine kc and kb ranges for this proc depending on FFT's to be done */
/*     on this processor                                                      */

  isum = 0;
  for(kc=1;kc<=nkf3;kc++){
   for(kb=1;kb<=nkf2;kb++){
     isum++;
     for(iproc=1;iproc<=num_proc;iproc++){
      if(isum==str_fft_ka_proc_all[iproc]){
        skb_fft_ka_proc_all[iproc]=kb;
        skc_fft_ka_proc_all[iproc]=kc;
      }/*endif*/
      if(isum==end_fft_ka_proc_all[iproc]){
        ekb_fft_ka_proc_all[iproc]=kb;
        ekc_fft_ka_proc_all[iproc]=kc;
      }/*endif*/
     }/*endfor*/
   }/*endfor*/
  }/*endfor*/

  skb_fft_ka_proc = skb_fft_ka_proc_all[myidp1];
  skc_fft_ka_proc = skc_fft_ka_proc_all[myidp1];
  ekb_fft_ka_proc = ekb_fft_ka_proc_all[myidp1];
  ekc_fft_ka_proc = ekc_fft_ka_proc_all[myidp1];
  

/*---------------------------------------------------------------------------*/
/* iii) In place-block transpose                                             */

 if(num_proc>1){

   /*---------------------------------------------*/
   /*   i)Loop over ka's you have                 */
   /*  ii)Loop over procs=iproc                   */
   /* iii)Loop over the kc range to send to iproc */
   /*  iv)Loop over the kb range to send to iproc */
   /*   v)Count                                   */

   sendcounts_fft_ka = (int *) cmalloc(num_proc*sizeof(int))-1;
   recvdspls_fft_ka  = (int *) cmalloc(num_proc*sizeof(int))-1; 
   senddspls_fft_ka  = (int *) cmalloc(num_proc*sizeof(int))-1;
   recvcounts_fft_ka = (int *) cmalloc(num_proc*sizeof(int))-1;


   kkk = 0;
   j=0;for(i=1;i<=num_proc;i++){sendcounts_fft_ka[i]=0;}
   for(ka=ska_fft_kb_proc;ka<=eka_fft_kb_proc;ka++){
     for(iproc=1;iproc<=num_proc;iproc++){
       skc_fft_kb_now = (ka==ska_fft_kb_proc ? skc_fft_kb_proc : 1);
       ekc_fft_kb_now = (ka==eka_fft_kb_proc ? ekc_fft_kb_proc : nkf3);
       kc_str = MAX(skc_fft_kb_now,skc_fft_ka_proc_all[iproc]);
       kc_end = MIN(ekc_fft_kb_now,ekc_fft_ka_proc_all[iproc]);
       for(kc=kc_str;kc<=kc_end;kc++){
         skb_fft_ka_now = (kc==skc_fft_ka_proc_all[iproc] ? 
                               skb_fft_ka_proc_all[iproc] : 1);
         ekb_fft_ka_now = (kc==ekc_fft_ka_proc_all[iproc] ? 
                               ekb_fft_ka_proc_all[iproc] : nkf2);
         kb_str = MAX(1,skb_fft_ka_now);
         kb_end = MIN(nkf2,ekb_fft_ka_now);
         for(kb=kb_str;kb<=kb_end;kb++){
           sendcounts_fft_ka[iproc]++;kkk++;
         }/*endfor*/
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/
   iii = nfft_kb_proc*nkf2;
   
   Alltoall(&sendcounts_fft_ka[1],1,MPI_INT,
            &recvcounts_fft_ka[1],1,MPI_INT,comm); 

   recvdspls_fft_ka[1] = 0;
   senddspls_fft_ka[1] = 0;
   for(i=1;i<=num_proc-1;i++){
     recvdspls_fft_ka[(i+1)] = recvdspls_fft_ka[i]+recvcounts_fft_ka[i];
     senddspls_fft_ka[(i+1)] = senddspls_fft_ka[i]+sendcounts_fft_ka[i];
   }/*endfor*/

 }/*endif*/

/*==========================================================================*/
/* IV) Find the Displs and Recv counts for rho allgather                    */


 if(num_proc>1){

    recv_counts_rho = (int *) cmalloc(num_proc*sizeof(int))-1;
    displs_rho      = (int *) cmalloc(num_proc*sizeof(int))-1; 

    nfft2_proc   = nfft_ka_proc*nkf1;
    Allgather(&nfft2_proc,1,MPI_INT,&(rev_counts_rho[1]),num_proc,MPI_INT,0,comm);

    displs_rho[1] = 0;
    for(i=2;i<=num_proc;i++){
     displs_rho[i] = recv_counts_rho[(i-1)]+displs_rho[(i-1)];
    }/*endfor*/

    recv_counts_coef = (int *) cmalloc(num_proc*sizeof(int))-1;
    ncoef2_proc = 2*nfft_kc_proc*nkf3;
    Allgather(&ncoef2_proc,1,MPI_INT,&(rev_counts_coef[1]),num_proc,MPI_INT,0,comm);

 }/*endif*/

/*==========================================================================*/
/* V) Package the useful information*/

/*-----------------------------------------------------------------------*/
/* Required Grid size  */

  nfft_size  = MAX3((2*nfft_kc_proc*nkf3),
                    (2*nfft_kb_proc*nkf2),
                    (2*nfft_ka_proc*nkf1));
  if(nfft_size % 2 == 0)nfft_size++;
  para_fft_pkg3d->nfft_size   = nfft_size;
  para_fft_pkg3d->nfft        = 2*nfft_ka*nkf1;
  para_fft_pkg3d->nfft_proc   = 2*nfft_ka_proc*nkf1;
  para_fft_pkg3d->ndata_ka    = 2*nfft_ka_proc*nkf1;
  para_fft_pkg3d->ndata_kb    = 2*nfft_kb_proc*nkf2;
  para_fft_pkg3d->ndata_kc    = 2*nfft_kc_proc*nkf3;

/*-----------------------------------------------------------------------*/
/* Refurbish/repack the Maps         */

  for(i=1;i<=ncoef_proc;i++){map_proc[i]      = 2*map_proc[i]-1;}
  for(i=1;i<=ncoef_use;i++) {map_c_proc[i]    = 2*map_c_proc[i]-1;}
  if(num_proc>1){
    iii = (recvcounts_fft_kc[num_proc]+recvdspls_fft_kc[num_proc]);
    for(i=1,m=1;i<=iii;i++,m+=2){
      map_inv[m]     = 2*map_proc_post[i]-1;
      map_inv[(m+1)] = 2*map_proc_post[i];
    }/*endfor*/
    for(i=1;i<=(2*iii);i++){map_proc_post[i] = map_inv[i];}
  }/*endif*/

  para_fft_pkg3d->map_proc      = map_proc;    
  para_fft_pkg3d->map_c_proc    = map_c_proc;  
  para_fft_pkg3d->map_proc_post = map_proc_post;

/*-----------------------------------------------------------------------*/
/* Rho receive and send counts                                           */

  if(num_proc>1){
    para_fft_pkg3d->recv_counts_rho = recv_counts_rho;
    para_fft_pkg3d->displs_rho      = displs_rho;
    para_fft_pkg3d->recv_counts_coef = recv_counts_coef;
  }/*endif*/

/*-----------------------------------------------------------------------*/
/* FFT along kc  */

  if(num_proc>1){
    for(i=1;i<=num_proc;i++){sendcounts_fft_kc[i] *= 2;}
    for(i=1;i<=num_proc;i++){recvcounts_fft_kc[i] *= 2;}
    for(i=1;i<=num_proc;i++){senddspls_fft_kc[i]  *= 2;}
    for(i=1;i<=num_proc;i++){recvdspls_fft_kc[i]  *= 2;}
  }/*endif*/
  para_fft_pkg3d->nfft_kc           = nfft_kc;
  para_fft_pkg3d->nfft_kc_proc      = nfft_kc_proc;
  para_fft_pkg3d->str_fft_kc_proc   = str_fft_kc_proc;
  para_fft_pkg3d->end_fft_kc_proc   = end_fft_kc_proc;
  para_fft_pkg3d->ska_fft_kc_proc   = ska_fft_kc_proc;
  para_fft_pkg3d->eka_fft_kc_proc   = eka_fft_kc_proc;
  para_fft_pkg3d->skb_fft_kc_proc   = skb_fft_kc_proc;
  para_fft_pkg3d->ekb_fft_kc_proc   = ekb_fft_kc_proc;
  para_fft_pkg3d->nka_fft_kc        = nka_fft_kc;
  para_fft_pkg3d->sendcounts_fft_kc = sendcounts_fft_kc;
  para_fft_pkg3d->senddspls_fft_kc  = senddspls_fft_kc;
  para_fft_pkg3d->recvcounts_fft_kc = recvcounts_fft_kc;
  para_fft_pkg3d->recvdspls_fft_kc  = recvdspls_fft_kc;
  para_fft_pkg3d->sum_fft_kc_proc   = sum_fft_kc_proc;
  para_fft_pkg3d->nkb_fft_kc        = nkb_fft_kc;
  para_fft_pkg3d->ka_fft_kc         = ka_fft_kc;
  para_fft_pkg3d->kb_fft_kc         = kb_fft_kc;
  para_fft_pkg3d->ka_fft_kc_red     = ka_fft_kc_red;

/*-----------------------------------------------------------------------*/
/* FFT along kb  */

  if(num_proc>1){
    for(i=1;i<=num_proc;i++){sendcounts_fft_kb[i] *= 2;}
    for(i=1;i<=num_proc;i++){recvcounts_fft_kb[i] *= 2;}
    for(i=1;i<=num_proc;i++){senddspls_fft_kb[i]  *= 2;}
    for(i=1;i<=num_proc;i++){recvdspls_fft_kb[i]  *= 2;}
  }/*endif*/

  para_fft_pkg3d->nfft_kb             =  nfft_kb;
  para_fft_pkg3d->nfft_kb_proc        =  nfft_kb_proc;
  para_fft_pkg3d->str_fft_kb_proc     =  str_fft_kb_proc;
  para_fft_pkg3d->end_fft_kb_proc     =  end_fft_kb_proc;
  para_fft_pkg3d->ska_fft_kb_proc     =  ska_fft_kb_proc;
  para_fft_pkg3d->eka_fft_kb_proc     =  eka_fft_kb_proc;
  para_fft_pkg3d->skc_fft_kb_proc     =  skc_fft_kb_proc;
  para_fft_pkg3d->ekc_fft_kb_proc     =  ekc_fft_kb_proc;
  para_fft_pkg3d->ska_fft_kb_proc_all = ska_fft_kb_proc_all;
  para_fft_pkg3d->eka_fft_kb_proc_all = eka_fft_kb_proc_all;
  para_fft_pkg3d->skc_fft_kb_proc_all = skc_fft_kb_proc_all;
  para_fft_pkg3d->ekc_fft_kb_proc_all = ekc_fft_kb_proc_all;
  para_fft_pkg3d->sendcounts_fft_kb   = sendcounts_fft_kb;
  para_fft_pkg3d->senddspls_fft_kb    = senddspls_fft_kb;
  para_fft_pkg3d->recvcounts_fft_kb   = recvcounts_fft_kb;
  para_fft_pkg3d->recvdspls_fft_kb    = recvdspls_fft_kb;
  para_fft_pkg3d->sum_fft_kb_proc     = sum_fft_kb_proc;

/*-----------------------------------------------------------------------*/
/* FFT along ka */

  if(num_proc>1){
    for(i=1;i<=num_proc;i++){sendcounts_fft_ka[i] *= 2;}
    for(i=1;i<=num_proc;i++){recvcounts_fft_ka[i] *= 2;}
    for(i=1;i<=num_proc;i++){senddspls_fft_ka[i]  *= 2;}
    for(i=1;i<=num_proc;i++){recvdspls_fft_ka[i]  *= 2;}
  }/*endif*/ 

  para_fft_pkg3d->nfft_ka             =  nfft_ka;
  para_fft_pkg3d->nfft_ka_proc        =  nfft_ka_proc;
  para_fft_pkg3d->str_fft_ka_proc     =  str_fft_ka_proc;
  para_fft_pkg3d->end_fft_ka_proc     =  end_fft_ka_proc;
  para_fft_pkg3d->skb_fft_ka_proc     =  skb_fft_ka_proc;
  para_fft_pkg3d->ekb_fft_ka_proc     =  ekb_fft_ka_proc;
  para_fft_pkg3d->skc_fft_ka_proc     =  skc_fft_ka_proc;
  para_fft_pkg3d->ekc_fft_ka_proc     =  ekc_fft_ka_proc;
  para_fft_pkg3d->skb_fft_ka_proc_all = skb_fft_ka_proc_all;
  para_fft_pkg3d->ekb_fft_ka_proc_all = ekb_fft_ka_proc_all;
  para_fft_pkg3d->skc_fft_ka_proc_all = skc_fft_ka_proc_all;
  para_fft_pkg3d->ekc_fft_ka_proc_all = ekc_fft_ka_proc_all;
  para_fft_pkg3d->sendcounts_fft_ka   = sendcounts_fft_ka;
  para_fft_pkg3d->senddspls_fft_ka    = senddspls_fft_ka;
  para_fft_pkg3d->recvcounts_fft_ka   = recvcounts_fft_ka;
  para_fft_pkg3d->recvdspls_fft_ka    = recvdspls_fft_ka;

/*==========================================================================*/
/*  VI) Initialize the FFTs                                                 */

  para_fft_gen3d_init(para_fft_pkg3d);

/*==========================================================================*/
/* VII) Free the excess memory                                              */

  free(&map[1]);
  free(&map_c[1]);
  free(&map_c_inv[1]);
  free(&indx[1]);
  free(&indx_c[1]);
  free(&ifft_strt_kc[1]);
  free(&nkb_fft_kc_low[1]);
  free(&nkb_fft_kc_high[1]);
  free(&nfft_kc_proc_all[1]);
  free(&nfft_kb_proc_all[1]);
  free(&nfft_ka_proc_all[1]);
  free(&map_inv[1]);

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/









/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void setfft_indx(int nkf1,int nkf2,int nkf3,int ncoef,
                 int *kastr_sm,int *kbstr_sm,int *kcstr_sm,
                 int *indx,int *indx_c)

/*==========================================================================*/
/*       Begin routine */
   {/*begin routine */
/*==========================================================================*/
/* Local variables */

  int i,ka,kb,kc,kap,kbp,kcp;

/*==========================================================================*/

  for(i=1;i<=ncoef;i++){

    ka = kastr_sm[i];
    kb = kbstr_sm[i];
    kc = kcstr_sm[i];

    kap = (ka < 0 ? ka + nkf1 + 1 : ka + 1);
    kbp = (kb < 0 ? kb + nkf2 + 1 : kb + 1);
    kcp = (kc < 0 ? kc + nkf3 + 1 : kc + 1);
    indx[i] = kcp*2 - 1 + (kbp - 1)*2*nkf3 
            + (kap - 1)*2*nkf3*nkf2;

    kap = (-ka < 0 ? -ka + nkf1 + 1 : -ka + 1);
    kbp = (-kb < 0 ? -kb + nkf2 + 1 : -kb + 1);
    kcp = (-kc < 0 ? -kc + nkf3 + 1 : -kc + 1);
    indx_c[i] = kcp*2 - 1 + (kbp - 1)*2*nkf3 
              + (kap - 1)*2*nkf3*nkf2;

  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* setkvec3d_sm */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void sort_commence(int n, int index[],int jndex[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
       jndex[1]=rjndex;
       break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
       index[i]=index[j];
       jndex[i]=jndex[j];
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void para_fft_gen3d_init(PARA_FFT_PKG3D *para_fft_pkg3d)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*   Local Variables */

  int iopt;
  int ier  = 0;
  int incl = 1;
  int incn,num;

/*   Local pointers */
  int nfft_kc_proc = para_fft_pkg3d->nfft_kc_proc;
  int nkf_c        = para_fft_pkg3d->nkf3;
  int nfft_kb_proc = para_fft_pkg3d->nfft_kb_proc;
  int nkf_b        = para_fft_pkg3d->nkf2;
  int nfft_ka_proc = para_fft_pkg3d->nfft_ka_proc;
  int nkf_a        = para_fft_pkg3d->nkf1;
  int igeneric_opt = para_fft_pkg3d->igeneric_opt;
  int scale_opt;

  int nwork1;
  int nwork2;
  int len_ifax;

  int *ifax_a_f,*ifax_a_r;
  int *ifax_b_f,*ifax_b_r;
  int *ifax_c_f,*ifax_c_r;

  double *work_1a_f,*work_2a_f;
  double *work_1a_r,*work_2a_r;
  double *work_1b_f,*work_2b_f;
  double *work_1b_r,*work_2b_r;
  double *work_1c_f,*work_2c_f;
  double *work_1c_r,*work_2c_r;

/*==========================================================================*/
/* I) Hard-wired sizes of scratch arrays                                    */

  para_fft_pkg3d->nwork1   = 60000; 
  para_fft_pkg3d->nwork2   = 60000; 
  para_fft_pkg3d->len_ifax = 13;
  nwork1                   = para_fft_pkg3d->nwork1;
  nwork2                   = para_fft_pkg3d->nwork2;
  len_ifax                 = para_fft_pkg3d->len_ifax;

/*==========================================================================*/
/* II) Allocate work and integer scratch arrays                              */

  para_fft_pkg3d->ifax_a_f    = (int *) cmalloc(len_ifax*sizeof(int))-1;
  para_fft_pkg3d->ifax_b_f    = (int *) cmalloc(len_ifax*sizeof(int))-1;
  para_fft_pkg3d->ifax_c_f    = (int *) cmalloc(len_ifax*sizeof(int))-1;

  para_fft_pkg3d->ifax_a_r    = (int *) cmalloc(len_ifax*sizeof(int))-1;
  para_fft_pkg3d->ifax_b_r    = (int *) cmalloc(len_ifax*sizeof(int))-1;
  para_fft_pkg3d->ifax_c_r    = (int *) cmalloc(len_ifax*sizeof(int))-1;

  para_fft_pkg3d->work_1a_f = (double *) cmalloc(nwork1*sizeof(double))-1;
  para_fft_pkg3d->work_1b_f = (double *) cmalloc(nwork1*sizeof(double))-1;
  para_fft_pkg3d->work_1c_f = (double *) cmalloc(nwork1*sizeof(double))-1;

  para_fft_pkg3d->work_2a_f = (double *) cmalloc(nwork2*sizeof(double))-1;
  para_fft_pkg3d->work_2b_f = (double *) cmalloc(nwork2*sizeof(double))-1;
  para_fft_pkg3d->work_2c_f = (double *) cmalloc(nwork2*sizeof(double))-1;

  para_fft_pkg3d->work_1a_r = (double *) cmalloc(nwork1*sizeof(double))-1;
  para_fft_pkg3d->work_1b_r = (double *) cmalloc(nwork1*sizeof(double))-1;
  para_fft_pkg3d->work_1c_r = (double *) cmalloc(nwork1*sizeof(double))-1;

  para_fft_pkg3d->work_2a_r = (double *) cmalloc(nwork2*sizeof(double))-1;
  para_fft_pkg3d->work_2b_r = (double *) cmalloc(nwork2*sizeof(double))-1;
  para_fft_pkg3d->work_2c_r = (double *) cmalloc(nwork2*sizeof(double))-1;

/*==========================================================================*/
/* III) Assign local pointers                                               */

  ifax_a_f = para_fft_pkg3d->ifax_a_f;  
  ifax_b_f = para_fft_pkg3d->ifax_b_f;  
  ifax_c_f = para_fft_pkg3d->ifax_c_f;  

  ifax_a_r = para_fft_pkg3d->ifax_a_r;  
  ifax_b_r = para_fft_pkg3d->ifax_b_r;  
  ifax_c_r = para_fft_pkg3d->ifax_c_r;  

  work_1a_f = para_fft_pkg3d->work_1a_f;
  work_1b_f = para_fft_pkg3d->work_1b_f;
  work_1c_f = para_fft_pkg3d->work_1c_f;
  work_2a_f = para_fft_pkg3d->work_2a_f;
  work_2b_f = para_fft_pkg3d->work_2b_f;
  work_2c_f = para_fft_pkg3d->work_2c_f;

  work_1a_r = para_fft_pkg3d->work_1a_r;
  work_1b_r = para_fft_pkg3d->work_1b_r;
  work_1c_r = para_fft_pkg3d->work_1c_r;
  work_2a_r = para_fft_pkg3d->work_2a_r;
  work_2b_r = para_fft_pkg3d->work_2b_r;
  work_2c_r = para_fft_pkg3d->work_2c_r;

/*==========================================================================*/
/* IV) Initialize fwd FFT arrays                                            */

  iopt = 1;  incn = nkf_c; num = nfft_kc_proc;
  fft_gen1d_init(nkf_c,incl,num,incn,iopt,&ier,
                 work_1c_f,nwork1,work_2c_f,nwork2,ifax_c_f,&scale_opt,
                 igeneric_opt);

  iopt = 1; incn = nkf_b; num = nfft_kb_proc;
  fft_gen1d_init(nkf_b,incl,num,incn,iopt,&ier,
                 work_1b_f,nwork1,work_2b_f,nwork2,ifax_b_f,&scale_opt,
                 igeneric_opt);

  iopt = 1; incn = nkf_a; num = nfft_ka_proc;
  fft_gen1d_init(nkf_a,incl,num,incn,iopt,&ier,
                 work_1a_f,nwork1,work_2a_f,nwork2,ifax_a_f,&scale_opt,
                 igeneric_opt);

/*==========================================================================*/
/* V) Initialize bck FFT arrays                                            */

  iopt = -1; incn = nkf_a; num = nfft_ka_proc;
  fft_gen1d_init(nkf_a,incl,num,incn,iopt,&ier,
                 work_1a_r,nwork1,work_2a_r,nwork2,ifax_a_r,&scale_opt,
                 igeneric_opt);

  iopt = -1; incn = nkf_b; num = nfft_kb_proc;
  fft_gen1d_init(nkf_b,incl,num,incn,iopt,&ier,
                 work_1b_r,nwork1,work_2b_r,nwork2,ifax_b_r,&scale_opt,
                 igeneric_opt);

  iopt = -1;  incn = nkf_c; num = nfft_kc_proc;
  fft_gen1d_init(nkf_c,incl,num,incn,iopt,&ier,
                 work_1c_r,nwork1,work_2c_r,nwork2,ifax_c_r,&scale_opt,
                 igeneric_opt);

/*==========================================================================*/
/* VI) Tuck away scale option                                               */

  para_fft_pkg3d->scale_opt = scale_opt;

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* fft_gen1d_init : Fortran wrappers for CP FFT initialization              */
/*==========================================================================*/

void fft_gen1d_init(int nfft,int incl,int num,int incn,
                    int iopt,int *ier,double *work_1,int nwork_1,
                    double *work_2,int nwork_2,int *ifax,int *scale_opt,
                    int igeneric_opt)

/*==========================================================================*/
 {/*begin routine */
/*==========================================================================*/
/*   Local Variables */

 int init         = 1;
 int nfft_for     = nfft;
 int incl_for_x   = incl;
 int incn_for_x   = incn;
 int incl_for_y   = incl;
 int incn_for_y   = incn;
 int num_for      = num;
 int nwork_1_for  = nwork_1;
 int nwork_2_for  = nwork_2;
 int iopt_for     = iopt;
 double scale     = 1.0;
 double dum_x,dum_y;


/*==========================================================================*/
/* Machine Specific FFT */

 if(igeneric_opt==0){

#ifdef IBM_ESSL
   dcft(&init,&dum_x,&incl_for_x,&incn_for_x,&dum_y,&incl_for_y,&incn_for_y,
        &nfft_for,&num_for,&iopt_for,&scale,&(work_1[1]),&nwork_1_for,
        &(work_2[1]),&nwork_2_for);
   *scale_opt = 1;
#endif

#ifdef SGI_COMPLIB
   zfft1di_(&nfft_for,&(work_2[1]));
   *scale_opt = 1;
#endif

#ifdef HP_VECLIB
   *scale_opt = 0;
#endif

#ifdef C90
   CFTFAX(&nfft_for,&(ifax[1]),&(work_2[1]));
   *scale_opt = 1;
#endif

#ifdef T3E_SCILIB
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Library FFT's not Implemented for T3E\n");
       printf("Use generic_fft_option in input file\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
#endif

#ifdef DEC_ALPHA
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Library FFT's not Implemented for DEC_ALPHA\n");
       printf("Use generic_fft_option in input file\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
#endif


/*==========================================================================*/
/* Generic FFT */

 }else{
#ifdef T3E_SCILIB
   *scale_opt = 1;
   CFFTI(&nfft_for,&(work_2[1]));
#else
#ifdef JUNK
   *scale_opt = 1;
   DCFFTI(&nfft_for,&(work_2[1]));
#endif
#endif

 }/*endif*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* zffts_gen : Fortran wrappers for 1D CP FFT                              */
/*==========================================================================*/

void fft_gen1d(double *zfft,int nfft,int incl,int num,int incn,
               int iopt,int *ier,double *work_1,int nwork_1,
               double *work_2,int nwork_2,int *ifax,int igeneric_opt)

/*==========================================================================*/
 {/*begin routine */
/*==========================================================================*/
/*   Local Variables */

 int init        = 0;
 int nfft_for    = nfft;
 int incl_for_x  = incl;
 int incn_for_x  = incn;
 int incl_for_y  = incl;
 int incn_for_y  = incn;
 int num_for     = num;
 int nwork_1_for = nwork_1;
 int nwork_2_for = nwork_2;
 double scale    = 1.0;
 int iopt_for,i,istart;
 int ier_for;


/*==========================================================================*/
/* Machine Specific FFT */

 if(igeneric_opt==0){

#ifdef IBM_ESSL
    iopt_for = iopt;
    dcft(&init,&(zfft[1]),&incl_for_x,&incn_for_x,&(zfft[1]),
         &incl_for_y,&incn_for_y,&nfft_for,&num_for,&iopt_for,&scale,
         &(work_1[1]),&nwork_1_for,&(work_2[1]),&nwork_2_for);
#endif

#ifdef SGI_COMPLIB
    iopt_for = -iopt;
    for(i=1;i<=num;i++){
      istart = 1 + (i-1)*2*incn;
      zfft1d_(&iopt_for,&nfft_for,&(zfft[istart]),&incl_for_x,&(work_2[1]));
    }/*endfor*/
#endif

#ifdef HP_VECLIB
    iopt_for = iopt;
    zffts(&(zfft[1]),&nfft_for,&incl_for_x,&num_for,&incn_for_x,&iopt_for,
          &ier_for);
    *ier = ier_for;
    if(*ier!=0){
     printf("The HP is nice and dies if the radix conditions is wrong!!!\n");
     exit(1);
    }/*endif*/
#endif

#ifdef C90
   iopt_for = -iopt;
   CFFT99(&(zfft[1]),&(work_1[1]),&(work_2[1]),&(ifax[1]),&incl_for_x,
          &incn_for_x,&nfft_for,&num_for,&iopt_for);
#endif

#ifdef T3E_SCILIB
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Library FFT's not Implemented for T3E\n");
       printf("Use generic_fft_option in input file\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
#endif

#ifdef DEC_ALPHA
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Library FFT's not Implemented for DEC_ALPHA\n");
       printf("Use generic_fft_option in input file\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
#endif

/*==========================================================================*/
/* Generic FFT */

 }else{

   if(iopt == 1){

     for(i=1;i<=num;i++){
#ifdef T3E_SCILIB
       istart = 1 + (i-1)*2*incn;
       CFFTF(&nfft_for,&(zfft[istart]),&(work_2[1])); 
#else
#ifdef JUNK
       istart = 1 + (i-1)*2*incn;
       DCFFTF(&nfft_for,&(zfft[istart]),&(work_2[1])); 
#endif
#endif
     }/*endfor*/

   } else {

     for(i=1;i<=num;i++){
#ifdef T3E_SCILIB
       istart = 1 + (i-1)*2*incn;
       CFFTB(&nfft_for,&(zfft[istart]),&(work_2[1])); 
#else
#ifdef JUNK
       istart = 1 + (i-1)*2*incn;
       DCFFTB(&nfft_for,&(zfft[istart]),&(work_2[1])); 
#endif
#endif
     }/*endif*/

   }/* endif */

 }/* endif */

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/





