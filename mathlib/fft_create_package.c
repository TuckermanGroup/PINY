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

#define DEBUG_PME_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Create para fft package */
/*==========================================================================*/

void create_para_fft_pkg3d(PARA_FFT_PKG3D *para_fft_pkg3d,
                           int *kastr,int *kbstr,int *kcstr,
                           int cp_dual_grid_opt)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

  int i,ic,ic_c,i_00,i_0,ic_00,ic_0,m,iglenn;
  int ioff,ifft,iii,imall1,imall2;
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
 int nkf3       = para_fft_pkg3d->nkf3; 
 int nkf2       = para_fft_pkg3d->nkf2; 
 int nkf1       = para_fft_pkg3d->nkf1; 
 int ncoef_proc;
 int ncoef_use;
 int icoef_off;
 int icoef_strt;

/*=======================================================================*/
/* 0) Determine useful constants                                         */

   idiv                    = ncoef / num_proc;
   irem                    = ncoef % num_proc;
   ncoef_proc              = (myid <  irem ? idiv+1 : idiv);
   icoef_strt              = (myid <= irem ? myid*(idiv+1)+1 : 
                                irem*(idiv+1)+1+(myid-irem)*idiv);
   ncoef_use               = (myidp1==num_proc ? ncoef_proc-1 : ncoef_proc);
   icoef_off               = icoef_strt-1;

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
   
/* maps for dualed option */

 if(cp_dual_grid_opt >= 1){
  para_fft_pkg3d->send_counts_row_big_small = (int *) 
                                              cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->recv_counts_row_big_small = (int *)
                                               cmalloc(num_proc*sizeof(int))-1;

  para_fft_pkg3d->sdispls_row_big_small = (int *) 
                                           cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->rdispls_row_big_small = (int *) 
                                           cmalloc(num_proc*sizeof(int))-1;

  para_fft_pkg3d->send_counts_row_big_small_upack = (int *) 
                                             cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->recv_counts_row_big_small_upack = (int *)
                                            cmalloc(num_proc*sizeof(int))-1;

  para_fft_pkg3d->sdispls_row_big_small_upack = (int *)
                                               cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->rdispls_row_big_small_upack = (int *)
                                               cmalloc(num_proc*sizeof(int))-1;


 para_fft_pkg3d->send_counts_dual_map = (int *)cmalloc(num_proc*sizeof(int))-1;
 para_fft_pkg3d->recv_counts_dual_map = (int *)cmalloc(num_proc*sizeof(int))-1;
 para_fft_pkg3d->sdispls_dual_map     = (int *)cmalloc(num_proc*sizeof(int))-1;
 para_fft_pkg3d->rdispls_dual_map     = (int *)cmalloc(num_proc*sizeof(int))-1;

  para_fft_pkg3d->send_counts_dual_map_upack = (int *)  
                                            cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->recv_counts_dual_map_upack = (int *) 
                                           cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->sdispls_dual_map_upack     = (int *) 
                                           cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->rdispls_dual_map_upack     = (int *) 
                                           cmalloc(num_proc*sizeof(int))-1;


 }/*endif cp_dual_grid_opt*/

 if(cp_dual_grid_opt == 2){

  para_fft_pkg3d->send_counts_ioff_big_small  = (int *) 
                                          cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->recv_counts_ioff_big_small  = (int *) 
                                          cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->sdispls_ioff_big_small  = (int *) 
                                          cmalloc(num_proc*sizeof(int))-1;
  para_fft_pkg3d->rdispls_ioff_big_small  = (int *) 
                                          cmalloc(num_proc*sizeof(int))-1;
 }/*endif cp_dual_grid_opt*/

/*========================================================================*/
/* I) Find FFT indices, sort them to make map_inv and map_c_inv           */
/*      then create standard map which are easier to handle               */
/*      and include the c-FFT offset filling                              */
/*------------------------------------------------------------------------*/
/*                  A Little Jingle:                                      */
/* Mapmaker mapmaker make me a map. Playing with indices a code           */
/* could get burned and it could get stuck for good.                      */
/* For mama make it vector for a papa a complex para_fft                  */
/* For me... Well i wouldn't holler if it were faster than anything.      */
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

  imall1 = MAX((nfft_kc_half+1),1);
  imall2 = MAX((2*nfft_kc_half-1),1);

  ifft_strt_kc     = (int *) cmalloc(imall1*sizeof(int))-1;
  ifft_strt_kc_c   = (int *) cmalloc(imall1*sizeof(int))-1;
  ka_fft_kc_red    = (int *) cmalloc(imall2*sizeof(int))-1;
  ka_fft_kc        = (int *) cmalloc(imall2*sizeof(int))-1;
  kb_fft_kc        = (int *) cmalloc(imall2*sizeof(int))-1;
  nkb_fft_kc       = (int *) cmalloc(imall2*sizeof(int))-1;
  nkb_fft_kc_low   = (int *) cmalloc(imall2*sizeof(int))-1;
  nkb_fft_kc_high  = (int *) cmalloc(imall2*sizeof(int))-1;

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
  /*  a)FFT with 0,0,kz -- negative kc values get mapped to FFT order      */

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
/*  i) Create maps for in place rearrangement on each proc                  */

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

  nfft_kc_proc_all    = (int *) cmalloc(num_proc*sizeof(int))-1;
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
/*      where the kc FFTs start/end.  Part of spherical truncation        */

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
    }/* endfor */

  }/*endif:num_proc>1*/

   
/*==========================================================================*/
/* III) Transpose II  :                                                     */

/*-------------------------------------------------------------------------*/
/* i) Determine the number of kb FFT's done on each proc                   */

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
  sum_fft_kb_proc     = (int *) cmalloc(nkf3*sizeof(int))-1;

  for(iproc=1;iproc<=num_proc;iproc++){
      iprocm1 = iproc-1;
      nfft_kb_proc_all[iproc]    = (iprocm1 < irem ? idiv+1 : idiv);
      str_fft_kb_proc_all[iproc] = (iprocm1 <= irem ? iprocm1*(idiv+1)+1 : 
                                    irem*(idiv+1)+1+(iprocm1-irem)*idiv);
      end_fft_kb_proc_all[iproc] = str_fft_kb_proc_all[iproc] 
                                 + nfft_kb_proc_all[iproc] - 1;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
/* ii) Determine kc,ka ranges for this proc depending on FFT's to be done */
/*     on this processor */

  isum = 0;
  for(kc=1;kc<=nkf3;kc++){
   for(ka=1;ka<=nka_fft_kc;ka++){
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
  for(kc=skc_fft_kb_proc;kc<=ekc_fft_kb_proc;kc++){
    sum_fft_kb_proc[kc] = isum;
    kstr = (kc==skc_fft_kb_proc ? ska_fft_kb_proc : 1);
    kend = (kc==ekc_fft_kb_proc ? eka_fft_kb_proc : nka_fft_kc);
    for(ka=kstr;ka<=kend;ka++){isum++;}
  }/*endfor*/

/*-------------------------------------------------------------------------*/
/* iii) In place-block transpose                                           */

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
      kc_str = MAX(1,   skc_fft_kb_proc_all[iproc]);
      kc_end = MIN(nkf3,ekc_fft_kb_proc_all[iproc]);
      for(kc=kc_str;kc<=kc_end;kc++){
        ka_str = (kc==skc_fft_kb_proc_all[iproc] ? 
                      ska_fft_kb_proc_all[iproc] : 1);
        ka_end = (kc==ekc_fft_kb_proc_all[iproc] ? 
                      eka_fft_kb_proc_all[iproc] : nka_fft_kc);
        ka_str = MAX(ska_fft_kc_proc,ka_str);
        ka_end = MIN(eka_fft_kc_proc,ka_end);
        for(ka=ka_str;ka<=ka_end;ka++){
          kb_str = (ka==ska_fft_kc_proc ? skb_fft_kc_proc : 1);
          kb_end = (ka==eka_fft_kc_proc ? ekb_fft_kc_proc : nkb_fft_kc[ka]);
          for(kb=kb_str;kb<=kb_end;kb++){sendcounts_fft_kb[iproc]++;}
        }/*endfor:ka*/
      }/*endfor:kc*/
    }/*endfor:iproc*/

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

  nfft_ka_proc_all    = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  str_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  end_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  skb_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  skc_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  ekb_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;
  ekc_fft_ka_proc_all = (int *) cmalloc((num_proc+1)*sizeof(int))-1;

  for(iproc=1;iproc<=num_proc;iproc++){
    iprocm1 = iproc-1;
    nfft_ka_proc_all[iproc]    = (iprocm1 < irem ? idiv+1 : idiv);
    str_fft_ka_proc_all[iproc] = (iprocm1 <= irem ? iprocm1*(idiv+1)+1 : 
                                      irem*(idiv+1)+1+(iprocm1-irem)*idiv);
    end_fft_ka_proc_all[iproc] = str_fft_ka_proc_all[iproc] 
                               + nfft_ka_proc_all[iproc] - 1;
  }/*endfor*/

/*-------------------------------------------------------------------------*/
/* ii) Determine kc,kb ranges for this proc depending on FFT's to be done */
/*     on this processor                                                 */

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
  
/*--------------------------------------------------------------------------*/
/* iii) In place-block transpose                                            */

 if(num_proc>1){

   /*---------------------------------------------*/
   /*   i)Loop over procs=iproc                   */
   /*  ii)Loop over the kc range to send to iproc */
   /* iii)Loop over the kb range to send to iproc */
   /*  iv)Loop over ka's you have                 */
   /*   v)Count                                   */

   sendcounts_fft_ka = (int *) cmalloc(num_proc*sizeof(int))-1;
   recvdspls_fft_ka  = (int *) cmalloc(num_proc*sizeof(int))-1; 
   senddspls_fft_ka  = (int *) cmalloc(num_proc*sizeof(int))-1;
   recvcounts_fft_ka = (int *) cmalloc(num_proc*sizeof(int))-1;


   for(i=1;i<=num_proc;i++){sendcounts_fft_ka[i]=0;}
   for(iproc=1;iproc<=num_proc;iproc++){
     kc_str = MAX(skc_fft_kb_proc,skc_fft_ka_proc_all[iproc]);
     kc_end = MIN(ekc_fft_kb_proc,ekc_fft_ka_proc_all[iproc]);
     for(kc=kc_str;kc<=kc_end;kc++){
       ka_str = (kc==skc_fft_kb_proc ? ska_fft_kb_proc : 1);
       ka_end = (kc==ekc_fft_kb_proc ? eka_fft_kb_proc : nka_fft_kc);
       kb_str = (kc==skc_fft_ka_proc_all[iproc] ? 
                     skb_fft_ka_proc_all[iproc] : 1);
       kb_end = (kc==ekc_fft_ka_proc_all[iproc] ? 
                     ekb_fft_ka_proc_all[iproc] : nkf2);
       for(kb=kb_str;kb<=kb_end;kb++){
         for(ka=ka_str;ka<=ka_end;ka++){sendcounts_fft_ka[iproc]++;}
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

    Allgather(&nfft2_proc,1,MPI_INT,&(recv_counts_rho[1]),1,MPI_INT,0,comm);

    displs_rho[1] = 0;
    for(i=2;i<=num_proc;i++){
     displs_rho[i] = recv_counts_rho[(i-1)]+displs_rho[(i-1)];
    }/*endfor*/

    recv_counts_coef = (int *) cmalloc(num_proc*sizeof(int))-1;
    ncoef2_proc = ncoef_proc;

    Allgather(&ncoef2_proc,1,MPI_INT,&(recv_counts_coef[1]),1,MPI_INT,0,comm);

 }/*endif*/

/*==========================================================================*/
/*  V) Package the useful information                                       */

/*-----------------------------------------------------------------------*/
/* Required size --  Grid and otherwise  */

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
  para_fft_pkg3d->ncoef_proc          = ncoef_proc;
  para_fft_pkg3d->icoef_strt          = icoef_strt;
  para_fft_pkg3d->icoef_off           = icoef_off;
  para_fft_pkg3d->ncoef_use           = ncoef_use;

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
    para_fft_pkg3d->recv_counts_rho  = recv_counts_rho;
    para_fft_pkg3d->displs_rho       = displs_rho;
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
/*  VI) Initialize the FFTs                                                */

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
/*  Set the fft index */
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
/*  Sort one vector commensorate with another  */
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


/*=========================================================================*/
/* Create the pme/send recv information for hybrid PME                     */
/*=========================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=========================================================================*/

void create_pme_comm_hybr(int np_forc,int myid_forc,MPI_Comm comm,
                          PARA_FFT_PKG3D *pme_pkg)

/*=========================================================================*/
   {/*begin routine*/
/*=========================================================================*/
/*       Local variables                                                   */
#include "../typ_defs/typ_mask.h"

  int idiv,irem,ncoef_proc,icoef_strt,ncoef_use,icoef_off;

  int ncoef = pme_pkg->ncoef;
  int *recv_counts_coef,*recv_dspls_coef;

/*=========================================================================*/
/* I) Get the parallelization  BUT                                         */
/*    DONT OVER WRITE THE NCOEF_PROCS IN THE PACKAGE                       */

  idiv        = ncoef / np_forc;
  irem        = ncoef % np_forc;
  ncoef_proc  = (myid_forc <  irem ? idiv+1 : idiv);
  icoef_strt  = (myid_forc <= irem ? myid_forc*(idiv+1)+1 : 
                     irem*(idiv+1)+1+(myid_forc-irem)*idiv);
  ncoef_use   = (myid_forc+1==np_forc ? ncoef_proc-1 : ncoef_proc);
  icoef_off   = icoef_strt-1;

/*=========================================================================*/
/* II) Allgather the send counts and displacements                         */

  recv_counts_coef = (int *) cmalloc(np_forc*sizeof(int))-1;
  recv_dspls_coef  = (int *) cmalloc(np_forc*sizeof(int))-1;

  Allgather(&ncoef_proc,1,MPI_INT,&(recv_counts_coef[1]),1,MPI_INT,0,comm);
  Allgather(&icoef_off,1,MPI_INT,&(recv_dspls_coef[1]),1,MPI_INT,0,comm);

/*=========================================================================*/
/* III) Tuck the stuff away */

  pme_pkg->recv_counts_coef = recv_counts_coef;
  pme_pkg->recv_dspls_coef  = recv_dspls_coef;

/*-------------------------------------------------------------------------*/
   }/*end routine*/ 
/*=========================================================================*/



/*=========================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=========================================================================*/
/* Create the pme/send recv information to account for particle smear.     */
/*   Mirror mirror on the wall, what code has the most heinous             */
/*   index logic of all.                                                   */
/*=========================================================================*/

void create_pme_comm_full_g(int n_interp,PARA_FFT_PKG3D *pme_pkg)

/*=========================================================================*/
   {/*begin routine*/
/*=========================================================================*/
/*       Local variables                                                   */
#include "../typ_defs/typ_mask.h"

  int jjj,itot,iproc,jproc,iii,istart_off;
  int inow,n_ovlap,istart,index,istart_now,iprod;
  int nrecv_up,nrecv_dn,ndata,nrecv_up_max;
  int index_now;
  int ia,ib,ic;
  int nfft_size;
  int skc_now,ekc_now;
  int skb_now,ekb_now;
  int skc_use,ekc_use;
  int skb_use,ekb_use;

  /* Local pointers */
  int nkf1                 = pme_pkg->nkf1;
  int nkf2                 = pme_pkg->nkf2;
  int nkf3                 = pme_pkg->nkf3;
  int nfft_ka_proc         = pme_pkg->nfft_ka_proc;
  int *skc_fft_all         = pme_pkg->skc_fft_ka_proc_all;
  int *ekc_fft_all         = pme_pkg->ekc_fft_ka_proc_all;
  int *skb_fft_all         = pme_pkg->skb_fft_ka_proc_all;
  int *ekb_fft_all         = pme_pkg->ekb_fft_ka_proc_all;
  int skc_fft              = pme_pkg->skc_fft_ka_proc;
  int ekc_fft              = pme_pkg->ekc_fft_ka_proc;
  int skb_fft              = pme_pkg->skb_fft_ka_proc;
  int ekb_fft              = pme_pkg->ekb_fft_ka_proc;
  int num_proc             = pme_pkg->num_proc;
  int myid                 = pme_pkg->myid;
  int myidp1               = pme_pkg->myidp1;
  MPI_Comm comm            = pme_pkg->comm;

  /* malloc and tuck later */
  int *map_pme_up;
  int *recvcounts_pme_up;
  int *recvdspls_pme_up;
  int *sendcounts_pme_up;
  int *senddspls_pme_up;

  int *map_pme_dn;
  int *recvcounts_pme_dn;
  int *recvdspls_pme_dn;
  int *sendcounts_pme_dn;
  int *senddspls_pme_dn;

  int *senddspls_tmp;
  int *ndata_tot;


/*=========================================================================*/
/* 0) Print out the data */

#ifdef DEBUG_PME
  if(myid==0){
    printf("====================================\n");
    printf("Grid partitioning\n");
    printf("------------------------------------\n");
  }/*endif*/
  for(iproc=1;iproc<=num_proc;iproc++){
   if(myidp1==iproc){
      printf("proc %d skc %d ekc %d skb %d ekb %d \n",
                   iproc,skc_fft,ekc_fft,skb_fft,ekb_fft);
   }/*endif*/
   Barrier(comm);
  }/*endfor*/
  if(myid==0){
    scanf("%d",&iii);
    printf("====================================\n\n");
  }/*endif*/
  Barrier(comm);
#endif

/*=========================================================================*/
/* I) Set constants and Malloc                                             */

  n_ovlap     = 2*nkf1*( (n_interp-1)*nkf2 + skb_fft - 1 );
  ndata       = n_ovlap + nfft_ka_proc*2*nkf1;
  iprod = 2*nkf1*nkf2;

  recvcounts_pme_dn = (int *) cmalloc(num_proc*sizeof(int))-1;
  recvdspls_pme_dn  = (int *) cmalloc(num_proc*sizeof(int))-1;
  sendcounts_pme_dn = (int *) cmalloc(num_proc*sizeof(int))-1;
  senddspls_pme_dn  = (int *) cmalloc(num_proc*sizeof(int))-1;
  map_pme_dn        = (int *) cmalloc(num_proc*sizeof(int))-1;

  recvcounts_pme_up = (int *) cmalloc(num_proc*sizeof(int))-1;
  recvdspls_pme_up  = (int *) cmalloc(num_proc*sizeof(int))-1;
  sendcounts_pme_up = (int *) cmalloc(num_proc*sizeof(int))-1;
  senddspls_pme_up  = (int *) cmalloc(num_proc*sizeof(int))-1;
  map_pme_up        = (int *) cmalloc(num_proc*sizeof(int))-1;

  senddspls_tmp     = (int *) cmalloc(num_proc*sizeof(int))-1;
  ndata_tot         = (int *) cmalloc(num_proc*sizeof(int))-1;

  Allgather(&ndata,1,MPI_INT,&(ndata_tot[1]),1,MPI_INT,0,comm);

/*=========================================================================*/
/* II) Find lower pme grid information required by each processor          */
/*     and tabulate the send/recv counts and displs. Also tabulate         */
/*     where to stick the information that has been gathered.              */

  itot = 0;
  for(iproc=1;iproc<=num_proc;iproc++){

    jjj     = (myidp1 + iproc); /* recv from upstairs */
    jproc   = (jjj > num_proc ?  jjj-num_proc : jjj);
    recvdspls_pme_dn[jproc]  = itot;
    inow                  = 0;
    istart                = -1;
    index                 = -1;
    jjj                   = skc_fft_all[jproc]-n_interp+1;
   /*--------------------------------------------------*/
   /* i) collect jproc's lower wrapped overlap region: */     
    if(jjj < 1){
      skc_now = jjj+nkf3;
      ekc_now = nkf3;
      skc_use = MAX(skc_now,skc_fft);
      ekc_use = MIN(ekc_now,ekc_fft);
      for(ic=skc_use;ic<=ekc_use;ic++){
        skb_use = (ic == skc_fft ?  skb_fft : 1);
        ekb_use = (ic == ekc_fft ?  ekb_fft : nkf2);
        if(ic==skc_use){ 
          istart = (skc_use-skc_now)*iprod + (skb_use-1)*2*nkf1;
          if(skc_use==skc_fft){
            index  = (skb_use-skb_fft)*nkf1*2 + n_ovlap;
          }else{
            index  = (nkf2-skb_fft+1)*nkf1*2 + (skb_use-1)*nkf1*2
                    + (skc_use-skc_fft-1)*iprod + n_ovlap;
          }/*endif*/
        }/*endif*/
        for(ib=skb_use;ib<=ekb_use;ib++){
          itot += 2*nkf1;
          inow += 2*nkf1;
        }/*endfor*/
      }/*endfor*/
    }/*endif: lower*/
   /*-----------------------------------------------------*/
   /* ii) collect jproc's lower unwrapped overlap region: */     
    if(jproc !=myidp1){
      istart_off = (jjj < 1 ? (-jjj+1)*iprod : 0);
      skc_now = (jjj < 1 ?  1 : jjj);
      ekc_now = (skb_fft_all[jproc]==1 ? skc_fft_all[jproc]-1 : 
                                         skc_fft_all[jproc]);
      skc_use = MAX(skc_now,skc_fft);
      ekc_use = MIN(ekc_now,ekc_fft);
      for(ic=skc_use;ic<=ekc_use;ic++){
        skb_use = (ic == skc_fft ?  skb_fft : 1);
        ekb_use = (ic == ekc_fft ?  ekb_fft : nkf2);
        if(ic==skc_use){
          istart_now = (skc_use-skc_now)*iprod + (skb_use-1)*nkf1*2
                     + istart_off;
          istart     = (istart > -1 ? istart : istart_now);
          if(skc_use==skc_fft){
            index_now  = (skb_use-skb_fft)*nkf1*2 + n_ovlap;
          }else{
            index_now  = (nkf2-skb_fft+1)*nkf1*2 + (skb_use-1)*nkf1*2
                       + (skc_use-skc_fft-1)*iprod + n_ovlap;
          }/*endif*/
          index     = (index > -1 ? index : index_now);
        }/*endif*/
        for(ib=skb_use;ib<=ekb_use;ib++){
          itot  += 2*nkf1;
          inow  += 2*nkf1;
        }/*endfor*/
      }/*endfor*/
    }/*endif:upper*/
   /*--------------------------------------------*/
   /* iii) Store the information:                */     
    senddspls_tmp[jproc]     = (istart > -1 ? istart : 0);
    map_pme_dn[jproc]        = (index  > -1 ? index  : 0);
    recvcounts_pme_dn[jproc] = inow;

  }/*endfor: processors*/
  nrecv_dn = itot;

/*=========================================================================*/
/* III) Send off the lower counts and displs                               */

  Alltoall(&recvcounts_pme_dn[1],1,MPI_INT,
           &sendcounts_pme_dn[1],1,MPI_INT,comm); 
  Alltoall(&senddspls_tmp[1] ,1,MPI_INT,
           &senddspls_pme_dn[1] ,1,MPI_INT,comm); 

/*=========================================================================*/
/* IV) Find the upper pme grid information required by each processor      */
/*     and tabulate the send/recv counts and displs. Also tabulate         */
/*     where to stick the information that has been gathered.              */

  itot = 0;
  for(iproc=1;iproc<=num_proc;iproc++){
    recvdspls_pme_up[iproc]  = itot;
    inow                     = 0;
    istart                   = -1;
    index                    = -1;
    skb_now                  = ekb_fft_all[iproc]-n_interp+1;
    if((iproc!=myidp1) && (skc_fft == ekc_fft_all[iproc]) &&(skb_now<1)){
      skb_now += nkf2;
      ekb_now  = nkf2;
      skb_use  = MAX(skb_now,skb_fft);
      ekb_use  = ekb_now;
      if(skc_fft==ekc_fft){ekb_use  = MIN(ekb_now,ekb_fft);}
      for(ib=skb_use;ib<=ekb_use;ib++){
        if(ib==skb_use){
         index  = (skb_use-skb_fft)*2*nkf1 + n_ovlap;
         istart = ndata_tot[iproc] + (skb_use-ekb_fft_all[iproc]-1)*2*nkf1;
        }/*endif*/
        itot  += 2*nkf1;
        inow  += 2*nkf1;
      }/*endfor*/
    }/*endif*/
   /* iii) Store the information:                */     
    senddspls_tmp[iproc]     = (istart > -1 ? istart : 0);
    map_pme_up[iproc]        = (index  > -1 ? index  : 0);
    recvcounts_pme_up[iproc] = inow;
  }/*endfor: processors*/
  nrecv_up = itot;

/*=========================================================================*/
/* V) Send off the upper counts and displs                               */

  Alltoall(&recvcounts_pme_up[1],1,MPI_INT,
           &sendcounts_pme_up[1],1,MPI_INT,comm); 
  Alltoall(&senddspls_tmp[1] ,1,MPI_INT,
           &senddspls_pme_up[1] ,1,MPI_INT,comm); 

  Allreduce(&(nrecv_up),&(nrecv_up_max),1,MPI_INT,MPI_MAX,0,comm);

/*=========================================================================*/
/* VI) Tuck some stuff away and let others go free                         */

  pme_pkg->nrecv_dn           = nrecv_dn;
  pme_pkg->map_pme_dn         = map_pme_dn;
  pme_pkg->recvcounts_pme_dn  = recvcounts_pme_dn;
  pme_pkg->recvdspls_pme_dn   = recvdspls_pme_dn;
  pme_pkg->sendcounts_pme_dn  = sendcounts_pme_dn;
  pme_pkg->senddspls_pme_dn   = senddspls_pme_dn;

  pme_pkg->nrecv_up           = nrecv_up;
  pme_pkg->nrecv_up_max       = nrecv_up_max;
  pme_pkg->map_pme_up         = map_pme_up;
  pme_pkg->recvcounts_pme_up  = recvcounts_pme_up;
  pme_pkg->recvdspls_pme_up   = recvdspls_pme_up;
  pme_pkg->sendcounts_pme_up  = sendcounts_pme_up;
  pme_pkg->senddspls_pme_up   = senddspls_pme_up;

  nfft_size          = pme_pkg->nfft_size;
  nfft_size          = MAX3(nfft_size,nrecv_dn,nrecv_up);
  pme_pkg->nfft_size = nfft_size;

  cfree(&senddspls_tmp[1]);
  cfree(&ndata_tot[1]);

/*-------------------------------------------------------------------------*/
   }/*end routine*/
/*=========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void destroy_para_fft_pkg3d(PARA_FFT_PKG3D *para_fft_pkg3d)

/*==========================================================================*/
/*       Begin routine */
   {/*begin routine */
/*==========================================================================*/
/* Local variables */

  int num_proc  = para_fft_pkg3d->num_proc;
  int *ifax_a_f = para_fft_pkg3d->ifax_a_f;
  int *ifax_b_f = para_fft_pkg3d->ifax_b_f;
  int *ifax_c_f = para_fft_pkg3d->ifax_c_f;

  int *ifax_a_r = para_fft_pkg3d->ifax_a_r;
  int *ifax_b_r = para_fft_pkg3d->ifax_b_r;
  int *ifax_c_r = para_fft_pkg3d->ifax_c_r;

  double *work_1a_f = para_fft_pkg3d->work_1a_f;
  double *work_1b_f = para_fft_pkg3d->work_1b_f;
  double *work_1c_f = para_fft_pkg3d->work_1c_f;

  double *work_1a_r = para_fft_pkg3d->work_1a_r;
  double *work_1b_r = para_fft_pkg3d->work_1b_r;
  double *work_1c_r = para_fft_pkg3d->work_1c_r;

  double *work_2a_f = para_fft_pkg3d->work_2a_f;
  double *work_2b_f = para_fft_pkg3d->work_2b_f;
  double *work_2c_f = para_fft_pkg3d->work_2c_f;

  double *work_2a_r = para_fft_pkg3d->work_2a_r;
  double *work_2b_r = para_fft_pkg3d->work_2b_r;
  double *work_2c_r = para_fft_pkg3d->work_2c_r;

  int *ka_fft_kc_red = para_fft_pkg3d->ka_fft_kc_red;
  int *ka_fft_kc     = para_fft_pkg3d->ka_fft_kc;
  int *kb_fft_kc     = para_fft_pkg3d->kb_fft_kc;
  int *nkb_fft_kc    = para_fft_pkg3d->nkb_fft_kc;

  int *sendcounts_fft_ka = para_fft_pkg3d->sendcounts_fft_ka;
  int *senddspls_fft_ka  = para_fft_pkg3d->senddspls_fft_ka;
  int *recvcounts_fft_ka = para_fft_pkg3d->recvcounts_fft_ka;
  int *recvdspls_fft_ka  = para_fft_pkg3d->recvdspls_fft_ka;

  int *ska_fft_kb_proc_all = para_fft_pkg3d->ska_fft_kb_proc_all;
  int *eka_fft_kb_proc_all = para_fft_pkg3d->eka_fft_kb_proc_all;
  int *skc_fft_kb_proc_all = para_fft_pkg3d->skc_fft_kb_proc_all;
  int *ekc_fft_kb_proc_all = para_fft_pkg3d->ekc_fft_kb_proc_all;
  int *sum_fft_kb_proc     = para_fft_pkg3d->sum_fft_kb_proc;

  int *sendcounts_fft_kb = para_fft_pkg3d->sendcounts_fft_kb;
  int *senddspls_fft_kb  = para_fft_pkg3d->senddspls_fft_kb;
  int *recvcounts_fft_kb = para_fft_pkg3d->recvcounts_fft_kb;
  int *recvdspls_fft_kb  = para_fft_pkg3d->recvdspls_fft_kb;

  int *skb_fft_ka_proc_all = para_fft_pkg3d->skb_fft_ka_proc_all;
  int *ekb_fft_ka_proc_all = para_fft_pkg3d->ekb_fft_ka_proc_all;
  int *skc_fft_ka_proc_all = para_fft_pkg3d->skc_fft_ka_proc_all;
  int *ekc_fft_ka_proc_all = para_fft_pkg3d->ekc_fft_ka_proc_all;

  int *recv_counts_rho  = para_fft_pkg3d->recv_counts_rho;
  int *displs_rho       = para_fft_pkg3d->displs_rho;
  int *recv_counts_coef = para_fft_pkg3d->recv_counts_coef;

  int *map_proc         = para_fft_pkg3d->map_proc;
  int *map_c_proc       = para_fft_pkg3d->map_c_proc;
  int *map_proc_post    = para_fft_pkg3d->map_proc_post;

/*==========================================================================*/

  cfree(&(ifax_a_f[1]));
  cfree(&(ifax_b_f[1]));
  cfree(&(ifax_c_f[1]));

  cfree(&(ifax_a_r[1]));
  cfree(&(ifax_b_r[1]));
  cfree(&(ifax_c_r[1]));

  cfree(&(work_1a_f[1]));
  cfree(&(work_1b_f[1]));
  cfree(&(work_1c_f[1]));

  cfree(&(work_1a_r[1]));
  cfree(&(work_1b_r[1]));
  cfree(&(work_1c_r[1]));

  cfree(&(work_2a_f[1]));
  cfree(&(work_2b_f[1]));
  cfree(&(work_2c_f[1]));

  cfree(&(work_2a_r[1]));
  cfree(&(work_2b_r[1]));
  cfree(&(work_2c_r[1]));

  cfree(&(ka_fft_kc_red[1]));
  cfree(&(ka_fft_kc[1]));
  cfree(&(kb_fft_kc[1]));
  cfree(&(nkb_fft_kc[1]));

  cfree(&(ska_fft_kb_proc_all[1]));
  cfree(&(eka_fft_kb_proc_all[1]));
  cfree(&(skc_fft_kb_proc_all[1]));
  cfree(&(ekc_fft_kb_proc_all[1]));
  cfree(&(sum_fft_kb_proc[1]));

  if(num_proc>1) {
   cfree(&(sendcounts_fft_kb[1]));
   cfree(&(senddspls_fft_kb[1]));
   cfree(&(recvcounts_fft_kb[1]));
   cfree(&(recvdspls_fft_kb[1]));
  }/* endif */

  cfree(&(skb_fft_ka_proc_all[1]));
  cfree(&(ekb_fft_ka_proc_all[1]));
  cfree(&(skc_fft_ka_proc_all[1]));
  cfree(&(ekc_fft_ka_proc_all[1]));

  if(num_proc>1){
   cfree(&(sendcounts_fft_ka[1]));
   cfree(&(senddspls_fft_ka[1]));
   cfree(&(recvcounts_fft_ka[1]));
   cfree(&(recvdspls_fft_ka[1]));
  }/* endif */


  if(num_proc>1){
    cfree(&(recv_counts_rho[1]));
    cfree(&(displs_rho[1]));
    cfree(&(recv_counts_coef[1]));
  }/*endif*/

  cfree(&(map_proc[1]));
  cfree(&(map_c_proc[1]));
  cfree(&(map_proc_post[1]));

/*-----------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Initial the 3d FFT */
/*==========================================================================*/

void para_fft_gen3d_init(PARA_FFT_PKG3D *para_fft_pkg3d)

/*==========================================================================*/
{/*begin routine */
/*==========================================================================*/
/*   Local Variables */

  int iopt,iii;
  int ier  = 0;
  int incl = 1;
  int incn,num;
  size_t dim_k;

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
/* II) Allocate work and integer scratch arrays                             */

  dim_k = (size_t) len_ifax;
  para_fft_pkg3d->ifax_a_f    = (int *) cmalloc(dim_k*sizeof(int))-1;
  para_fft_pkg3d->ifax_b_f    = (int *) cmalloc(dim_k*sizeof(int))-1;
  para_fft_pkg3d->ifax_c_f    = (int *) cmalloc(dim_k*sizeof(int))-1;

  para_fft_pkg3d->ifax_a_r    = (int *) cmalloc(dim_k*sizeof(int))-1;
  para_fft_pkg3d->ifax_b_r    = (int *) cmalloc(dim_k*sizeof(int))-1;
  para_fft_pkg3d->ifax_c_r    = (int *) cmalloc(dim_k*sizeof(int))-1;

  dim_k = (size_t) nwork1;
  para_fft_pkg3d->work_1a_f = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_1b_f = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_1c_f = (double *) cmalloc(dim_k*sizeof(double))-1;

  para_fft_pkg3d->work_1a_r = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_1b_r = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_1c_r = (double *) cmalloc(dim_k*sizeof(double))-1;

  dim_k = (size_t) nwork2;
  para_fft_pkg3d->work_2a_f = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_2b_f = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_2c_f = (double *) cmalloc(dim_k*sizeof(double))-1;

  para_fft_pkg3d->work_2a_r = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_2b_r = (double *) cmalloc(dim_k*sizeof(double))-1;
  para_fft_pkg3d->work_2c_r = (double *) cmalloc(dim_k*sizeof(double))-1;

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

 int ifound;
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
 int isign,isys; /* for T3E*/


/*==========================================================================*/
/* Machine Specific FFT */

 if(igeneric_opt==0){

   ifound = 0;

#ifdef IBM_ESSL
   ifound++;
   dcft(&init,&dum_x,&incl_for_x,&incn_for_x,&dum_y,&incl_for_y,&incn_for_y,
        &nfft_for,&num_for,&iopt_for,&scale,&(work_1[1]),&nwork_1_for,
        &(work_2[1]),&nwork_2_for);
   *scale_opt = 1;
#endif

#ifdef SGI_COMPLIB
   ifound++;
   zfft1di_(&nfft_for,&(work_2[1]));
   *scale_opt = 1;
#endif

#ifdef SUN_COMPLIB
   ifound++;
   zffti(nfft_for,&(work_2[1]));
   *scale_opt = 1;
#endif

#ifdef HP_VECLIB
   ifound++;
   *scale_opt = 0;
#endif

#ifdef T3E_SCILIB
   ifound++;
   isign = 0;
   isys  = 0;
  CCFFT(&isign,&nfft_for,&scale,&dum_x,&dum_y,&(work_1[1]),&(work_2[1]),&isys);
   *scale_opt = 1;
#endif

#ifdef C90
   ifound++;
   CFTFAX(&nfft_for,&(ifax[1]),&(work_2[1]));
   *scale_opt = 1;
#endif

   if(ifound!=1){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Library FFT's not Implemented on this machine\n");
       printf("Use generic_fft_option in input file\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/

/*==========================================================================*/
/* Generic FFT */

 }else{
#ifdef T3E_SCILIB
   *scale_opt = 1;
   CFFTI_GENERIC(&nfft_for,&(work_2[1]));
#else
   *scale_opt = 1;
   DCFFTI_GENERIC(&nfft_for,&(work_2[1]));
#endif

 }/*endif*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/

