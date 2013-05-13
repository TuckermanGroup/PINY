/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*              Module: control_contract_rho.c                              */
/*                                                                          */
/*                                                                          */
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

#define REAL_CASE_CON
#define DEBUG_PME_CON_PAR1_OFF
#define ORIG_OFF
#define PME

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*Controller: call correct spread routine for rho when use dual grid opt    */
/*==========================================================================*/

void control_contract_rho(CPSCR *cpscr,double *rfft,
                          CELL *cell,int np_states,int n_interp_pme_dual,
                          PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                          PARA_FFT_PKG3D *para_fft_pkg3d_lg,
                          int cp_dual_grid_opt)
/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/* local pointers */
    
 double *zfft           =  cpscr->cpscr_wave.zfft;
 double *zfft_tmp       =  cpscr->cpscr_wave.zfft_tmp;
 int    *map_dual_upack =  cpscr->cpscr_rho.map_dual_upack;
 double cpu1,cpu2; 

 if(np_states == 1){
  switch(cp_dual_grid_opt){
   case 1:
     sngl_upack_rho_dual_ser(zfft,rfft,cell,para_fft_pkg3d_cp_box,
                            para_fft_pkg3d_lg);   
   break;
   case 2:

#ifdef  ORIG
     sngl_upack_rho_dual_ser(zfft,rfft,cell,para_fft_pkg3d_cp_box,
                             para_fft_pkg3d_lg);   
#endif

#ifdef  PME

 cputime(&cpu1);

#ifdef REAL_CASE_CON
     sngl_upack_rho_dual_ser_pme(zfft,rfft,cell,n_interp_pme_dual,
                                 cpscr,para_fft_pkg3d_cp_box,
                                 para_fft_pkg3d_lg);
#endif

#ifdef DEBUG_PME_CON_PAR1
printf("checking pme_par spread  routine in serial \n");

    sngl_upack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,
                                cell,n_interp_pme_dual,
                                &(cpscr->cpscr_dual_pme),
                                para_fft_pkg3d_cp_box,
                                para_fft_pkg3d_lg);
#endif

 cputime(&cpu2);

#ifdef TIME_ME
 printf("cpu time contract %g \n",cpu2-cpu1);
#endif

#endif
   break;
  }/*end switch*/

 }else{
  switch(cp_dual_grid_opt){
   case 1:
   sngl_upack_rho_dual_par(zfft,zfft_tmp,rfft,
                           map_dual_upack,
                           (cpscr->cpscr_rho.iset_map_upack_flag),cell,
                           para_fft_pkg3d_cp_box,
                           para_fft_pkg3d_lg);   
    break;
    case 2:
    sngl_upack_rho_dual_par_pme(&(cpscr->cpscr_wave),rfft,
                                cell,n_interp_pme_dual,
                                &(cpscr->cpscr_dual_pme),
                                para_fft_pkg3d_cp_box,
                                para_fft_pkg3d_lg);
     break;
  }/*end switch*/
 }/*endif np_states*/

/*--------------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl unpack the density for dualed system full_g parallel                */
/*==========================================================================*/

void sngl_upack_rho_dual_par(double *zfft,double *zfft_tmp,double *rfft,
                            int *map_dual_upack,int iset_map_upack_flag,
                            CELL *cell,
                            PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                            PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int i,iii;
 int nfft_proc  = para_fft_pkg3d_lg->nfft_proc;   /* on large grid*/

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 double *hmati             = cell->hmati;
 double *cp_box_center     = cell->cp_box_center;

 int ia_shift,ib_shift,ic_shift;
 double dia_shift,dib_shift,dic_shift;
 int index_lg_cmplx,index_sm,index_lg;
 

 double sx,sy,sz;

/*-----------------------------------------------------------------------*/
/* local pointers from fft packages                                      */

/* Large Box                                                             */

 int   skc_fft_ka_proc_lg = para_fft_pkg3d_lg->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_lg = para_fft_pkg3d_lg->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_lg = para_fft_pkg3d_lg->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_lg = para_fft_pkg3d_lg->ekb_fft_ka_proc;

/* Small CP Box                                                          */

 int   skc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skc_fft_ka_proc;
 int   skb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skb_fft_ka_proc;

 int *skb_fft_ka_proc_all_cp_box  = para_fft_pkg3d_cp_box->skb_fft_ka_proc_all;
 int *skc_fft_ka_proc_all_cp_box  = para_fft_pkg3d_cp_box->skc_fft_ka_proc_all;
 int *ekb_fft_ka_proc_all_cp_box  = para_fft_pkg3d_cp_box->ekb_fft_ka_proc_all;
 int *ekc_fft_ka_proc_all_cp_box  = para_fft_pkg3d_cp_box->ekc_fft_ka_proc_all;

 int num_proc = para_fft_pkg3d_lg->num_proc;

 int iproc;
 int kc_sm,kb_sm,ka_sm;
 int kb_str_sm_ind;
 int kc_lg,kb_lg,ka_lg;

 int myid             = para_fft_pkg3d_lg->myid;
 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 
 MPI_Comm comm_cp_box = para_fft_pkg3d_cp_box->comm; 

/*NOTE the send and recv counts and displs are exactly opposite for the pack */
/*  and the upack option so rather than remake them, just assign  the send   */
/*  to the recv                                                              */
#ifdef OLD
 int *send_counts_upack = para_fft_pkg3d_cp_box->recv_counts_row_big_small;
 int *recv_counts_upack = para_fft_pkg3d_cp_box->send_counts_row_big_small;
 int *sdispls_upack     = para_fft_pkg3d_cp_box->rdispls_row_big_small;
 int *rdispls_upack     = para_fft_pkg3d_cp_box->sdispls_row_big_small;

 int *send_counts_dual_map_upack = para_fft_pkg3d_cp_box->recv_counts_dual_map;
 int *recv_counts_dual_map_upack = para_fft_pkg3d_cp_box->send_counts_dual_map;
 int *sdispls_dual_map_upack     = para_fft_pkg3d_cp_box->rdispls_dual_map;
 int *rdispls_dual_map_upack     = para_fft_pkg3d_cp_box->sdispls_dual_map; 
#endif

 int *send_counts_upack =
                    para_fft_pkg3d_cp_box->send_counts_row_big_small_upack;
 int *recv_counts_upack = 
                    para_fft_pkg3d_cp_box->recv_counts_row_big_small_upack;
 int *sdispls_upack     = para_fft_pkg3d_cp_box->sdispls_row_big_small_upack;
 int *rdispls_upack     = para_fft_pkg3d_cp_box->rdispls_row_big_small_upack;

 int *send_counts_dual_map_upack = para_fft_pkg3d_cp_box->recv_counts_dual_map;
 int *recv_counts_dual_map_upack = para_fft_pkg3d_cp_box->send_counts_dual_map;
 int *sdispls_dual_map_upack     = para_fft_pkg3d_cp_box->rdispls_dual_map;
 int *rdispls_dual_map_upack     = para_fft_pkg3d_cp_box->sdispls_dual_map; 


 int isend_tot,icount;    
 int isend_tot_map;
 int kc_str_sm,kc_end_sm;
 int kb_str_sm,kb_end_sm;
 int kb_str_lg,kb_end_lg;
 int kb_off_sm,kc_off_sm;
 int kb_off,kc_off;
 int one= 1;
 int master = 0;

 int sum_count,sum_count_temp;
 int icount_tot_lg;
 int nfft_up_cp_box      = para_fft_pkg3d_cp_box->nfft;
 int nfft2_up_cp_box;
 int np_states = para_fft_pkg3d_lg->num_proc;

 int map_count;
 int map_count_proc;
 int *ind_tmp;
 double dnkf1,dnkf2,dnkf3;

   dnkf1 = (double) nkf1;
   dnkf2 = (double) nkf2;
   dnkf3 = (double) nkf3;

/*=======================================================================*/
   nfft2_up_cp_box =  nfft_up_cp_box/2;

   map_count        =  nfft_up_cp_box/nkf1_cp_box;   
 if(iset_map_upack_flag == 0){
   ind_tmp = (int *) cmalloc(map_count*sizeof(int))-1;
 }
/*=======================================================================*/
/* Determine useful constants                                            */

   sx = cp_box_center[1]*hmati[1]+cp_box_center[2]*hmati[4] 
      + cp_box_center[3]*hmati[7];

   sy = cp_box_center[1]*hmati[2]+cp_box_center[2]*hmati[5]
      + cp_box_center[3]*hmati[8];

   sz = cp_box_center[1]*hmati[3]+cp_box_center[2]*hmati[6]
      + cp_box_center[3]*hmati[9];

   dia_shift = ( -sx *(double)nkf1 - 0.5*(double)nkf1_cp_box);
   dib_shift = ( -sy *(double)nkf2 - 0.5*(double)nkf2_cp_box);
   dic_shift = ( -sz *(double)nkf3 - 0.5*(double)nkf3_cp_box);

   dia_shift /= dnkf1;
   dib_shift /= dnkf2;
   dic_shift /= dnkf3;

   dia_shift -= NINT((dia_shift - .5));
   dib_shift -= NINT((dib_shift - .5));
   dic_shift -= NINT((dic_shift - .5)); 

   ia_shift = (int)(dia_shift*dnkf1);
   ib_shift = (int)(dib_shift*dnkf2);
   ic_shift = (int)(dic_shift*dnkf3);

/*=======================================================================*/
/*=======================================================================*/
/* The boxes are stationary  this means the map only needs to be constructed */
/*   the very first time                                                     */

   for(i=1; i<= nfft_proc; i++){
     zfft_tmp[i] = 0.0; 
   } 


 if(iset_map_upack_flag == 0){
   for(iproc = 1; iproc<= num_proc ; iproc++){
      send_counts_upack[iproc]              = 0;
      recv_counts_upack[iproc]              = 0;
      send_counts_dual_map_upack[iproc]      = 0;
      recv_counts_dual_map_upack[iproc]      = 0;
   }
 }

   isend_tot         = 0;
   isend_tot_map     = 0;
   icount            = 0;

   for(iproc=1; iproc <= num_proc; iproc++){
      kc_str_sm = skc_fft_ka_proc_all_cp_box[iproc];
      kc_end_sm = ekc_fft_ka_proc_all_cp_box[iproc];
   for(kc_sm=kc_str_sm;kc_sm<=kc_end_sm;kc_sm++){
     kb_str_sm = (kc_sm==skc_fft_ka_proc_all_cp_box[iproc] ?
                         skb_fft_ka_proc_all_cp_box[iproc] : 1);
     kb_end_sm = (kc_sm==ekc_fft_ka_proc_all_cp_box[iproc] ?
                         ekb_fft_ka_proc_all_cp_box[iproc]: nkf2_cp_box);
                         
     for(kb_sm=kb_str_sm;kb_sm<=kb_end_sm;kb_sm++){
       kc_lg = kc_sm + ic_shift;
       kc_lg = (kc_lg > 0 ? kc_lg : kc_lg + nkf3);
       kc_lg = (kc_lg > nkf3 ? kc_lg - nkf3 : kc_lg );

         if( (kc_lg>=skc_fft_ka_proc_lg) && 
             (kc_lg<=ekc_fft_ka_proc_lg)){

           kb_lg     = kb_sm + ib_shift;
           kb_lg     = (kb_lg > 0 ? kb_lg : kb_lg + nkf2);
           kb_lg     = (kb_lg > nkf2 ? kb_lg -nkf2 : kb_lg );
           kb_str_lg = (kc_lg==skc_fft_ka_proc_lg ? 
                               skb_fft_ka_proc_lg : 1);
           kb_end_lg = (kc_lg==ekc_fft_ka_proc_lg ? 
                               ekb_fft_ka_proc_lg : nkf2);

           kb_str_sm_ind = (kc_sm==skc_fft_ka_proc_cp_box ? 
                                   skb_fft_ka_proc_cp_box : 1);

           if((kb_lg>=kb_str_lg) && (kb_lg<=kb_end_lg)){

             kc_off = kc_lg-skc_fft_ka_proc_lg;
             kb_off = kb_lg-kb_str_lg;

             index_lg  = kb_off*nkf1;

             kc_off_sm = kc_sm - kc_str_sm;
             kb_off_sm = kb_sm - kb_str_sm;

             index_sm = kb_off_sm*nkf1_cp_box;

             if(kc_off>0){
	        index_lg += (kc_off-1)*nkf2*nkf1
		         +  (nkf2-skb_fft_ka_proc_lg+1)*nkf1;  
             }/*endif*/

            if(kc_off_sm>0){
              index_sm += (kc_off_sm-1)*nkf2_cp_box*nkf1_cp_box
	               +  (nkf2_cp_box-skb_fft_ka_proc_all_cp_box[iproc]+1)
                         *nkf1_cp_box;   
             }/*endif*/

             if(iset_map_upack_flag==0){
                ind_tmp[isend_tot_map+1] = index_sm; 
	      }

             for(ka_sm=1;ka_sm<=nkf1_cp_box;ka_sm++){
               ka_lg =  ka_sm + ia_shift;
               ka_lg = (ka_lg > 0 ? ka_lg : ka_lg + nkf1);
               ka_lg = (ka_lg > nkf1 ? ka_lg - nkf1  : ka_lg );

               index_lg_cmplx = 2*(index_lg + ka_lg) - 1;

    	       zfft_tmp[isend_tot+ka_sm]     = zfft[index_lg_cmplx]; 

             }/*endfor*/
          if(iset_map_upack_flag == 0){
             send_counts_upack[iproc] += nkf1_cp_box;
             send_counts_dual_map_upack[iproc] += 1;
	   }
             isend_tot                += nkf1_cp_box;
             isend_tot_map            += 1;

	   }/*endif : kb in range*/
         }/*endif : kc in range */

        icount  += nkf1_cp_box;

     }/*endfor : kb*/
   }/*endfor : kc */
  }/*endfor : iproc */

/*=======================================================================*/
/* AllReduce the sum_counts_temp                                         */
/*  If send_counts does not equal the total number of elements of rho sm */
/*   there is an error                                                   */
if(iset_map_upack_flag == 0){
    sum_count_temp = 0;
  for(iproc=1; iproc<= num_proc; iproc++){
    sum_count_temp += send_counts_upack[iproc];
   }

  if(num_proc > 1){
    Allreduce(&(sum_count_temp),&(sum_count),
              1,MPI_INT,MPI_SUM,master,comm_lg);
  }else{
    sum_count = sum_count_temp;
  }
  
printf("sum_count %d nfft2_up_cp_box %d \n",sum_count,nfft2_up_cp_box);
  if( sum_count != nfft2_up_cp_box ){
   fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@ \n\n");
   fprintf(stderr,"Mismatch in the  number of values of rho to \n");
   fprintf(stderr,"be communicated in the routine that maps the\n");
   fprintf(stderr," density from the large grid  to the small   \n"); 
   fprintf(stderr," contact technical support\n");
   fprintf(stderr,"@@@@@@@@@@@@@@@@@@-error-@@@@@@@@@@@@@@@@@@@@@@@\n");
   fflush(stderr);
   exit(1);
  }/*endif*/

/*=======================================================================*/
/* I) Alltoall  the send_counts to the receive counts                    */
  if(num_proc > 1){
   Alltoall(&(send_counts_upack[1]),one,MPI_INT,
            &(recv_counts_upack[1]),one,MPI_INT,comm_cp_box);
   Alltoall(&(send_counts_dual_map_upack[1]),one,MPI_INT,
            &(recv_counts_dual_map_upack[1]),one,MPI_INT,comm_cp_box);
  }else{
    recv_counts_upack[1]         = send_counts_upack[1];
    recv_counts_dual_map_upack[1] = send_counts_dual_map_upack[1];
  }

/*-----------------------------------------------------------------------*/
/* I.V) Construct the displacements for the alltoallv                    */

    rdispls_upack[1] = 0;
    sdispls_upack[1] = 0;

    rdispls_dual_map_upack[1] = 0;
    sdispls_dual_map_upack[1] = 0;

  for(iproc=2; iproc<=num_proc; iproc++){
   rdispls_upack[iproc] = rdispls_upack[iproc-1]
                        + recv_counts_upack[iproc-1];

   sdispls_upack[iproc] = sdispls_upack[iproc-1] 
                        + send_counts_upack[iproc-1];

   rdispls_dual_map_upack[iproc] = rdispls_dual_map_upack[iproc-1]
                                + recv_counts_dual_map_upack[iproc-1];

   sdispls_dual_map_upack[iproc] = sdispls_dual_map_upack[iproc-1] 
                                + send_counts_dual_map_upack[iproc-1];
  }
 }/*endif iset_map_upack_flag */

/*-----------------------------------------------------------------------*/
/* Determine the total number of elements the individual processors have */

       icount_tot_lg = 0;
     for(iproc=1; iproc<= num_proc; iproc++){
       icount_tot_lg += recv_counts_upack[iproc];
     }

  map_count_proc = 0;
 for(iproc=1; iproc<= num_proc; iproc++){
   map_count_proc += recv_counts_dual_map_upack[iproc];
 }

/*-----------------------------------------------------------------------*/
/*III) Alltoallv the upack map                                           */

if(iset_map_upack_flag == 0){
 if(np_states>1){
   Alltoallv(&(ind_tmp[1]),&(send_counts_dual_map_upack[1]),
             &(sdispls_dual_map_upack[1]),MPI_INT,
             &(map_dual_upack[1]),&(recv_counts_dual_map_upack[1]),
             &(rdispls_dual_map_upack[1]),MPI_INT,comm_lg); 
 }else{
   for(i=1;i<= map_count_proc;i++){
      map_dual_upack[i] = ind_tmp[i];
   }

 }
    cfree(&(ind_tmp[1]));

}/*endif iset_map_upack_flag*/


/*-----------------------------------------------------------------------*/
/*III.V) Alltoallv the real space density on the small real space grid   */
/*  store resulting vector in zfft_tmp                                   */


 if(np_states>1 ){
   Alltoallv(&(zfft_tmp[1]),&(send_counts_upack[1]),
             &(sdispls_upack[1]),MPI_DOUBLE,
             &(zfft[1]),&(recv_counts_upack[1]),
             &(rdispls_upack[1]),MPI_DOUBLE,comm_lg);
 }else{
   for(i=1;i<= icount_tot_lg;i++){
      zfft[i] = zfft_tmp[i];
    }
 }

/*=======================================================================*/
/* Assign the small real space density elements to the appropriate       */
/*  on the large real space grid                                         */

   index_sm = 1;
 for(i=1; i<= map_count_proc; i++){
   for(ka_sm = 1; ka_sm <= nkf1_cp_box; ka_sm++){
     rfft[(map_dual_upack[i]+ka_sm)] = zfft[index_sm];
     index_sm++;
   }
 }

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl unpack the density for dualed system */
/*==========================================================================*/

void sngl_upack_rho_dual_ser(double *zfft,double *rfft,CELL *cell,
                            PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                            PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
 int iii;
 int ia,ib,ic;
 int iap,ibp,icp; 
  
 double dia_shift,dib_shift,dic_shift;
 int ia_shift,ib_shift,ic_shift;
 int index_lg_cmplx,index_sm,index_lg;

 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nkf2 = para_fft_pkg3d_lg->nkf2;
 int  nkf3 = para_fft_pkg3d_lg->nkf3;

 double *hmati             = cell->hmati;
 double *cp_box_center     = cell->cp_box_center;
 double dbox_rat;
 double sx,sy,sz;
 double dia,dib,dic;
 double dnkf1,dnkf2,dnkf3;
/*=======================================================================*/
/* Determine useful constants                                            */

  dnkf1 = (double) nkf1;
  dnkf2 = (double) nkf2;
  dnkf3 = (double) nkf3;

   dbox_rat = cell->hmat[1]/cell->hmat_cp[1];

   sx = cp_box_center[1]*hmati[1]+cp_box_center[2]*hmati[4] 
      + cp_box_center[3]*hmati[7];

   sy = cp_box_center[1]*hmati[2]+cp_box_center[2]*hmati[5]
      + cp_box_center[3]*hmati[8];

   sz = cp_box_center[1]*hmati[3]+cp_box_center[2]*hmati[6]
      + cp_box_center[3]*hmati[9];

   dia_shift = ( -sx *(double)nkf1 - 0.5*(double)nkf1_cp_box);
   dib_shift = ( -sy *(double)nkf2 - 0.5*(double)nkf2_cp_box);
   dic_shift = ( -sz *(double)nkf3 - 0.5*(double)nkf3_cp_box);

   dia_shift /= dnkf1;
   dib_shift /= dnkf2;
   dic_shift /= dnkf3;

   dia_shift -= NINT((dia_shift - .5));
   dib_shift -= NINT((dib_shift - .5));
   dic_shift -= NINT((dic_shift - .5)); 

   ia_shift = (int)(dia_shift*dnkf1);
   ib_shift = (int)(dib_shift*dnkf2);
   ic_shift = (int)(dic_shift*dnkf3);

/*=======================================================================*/
/* MAP density in real space from large grid to smaller cp grid          */

 for(ic=1; ic <= nkf3_cp_box; ic++){
  for(ib=1; ib <= nkf2_cp_box; ib++){
   for(ia=1; ia <= nkf1_cp_box; ia++){
      index_sm = (ia) + (ib - 1)*nkf1_cp_box
               + (ic - 1)*nkf1_cp_box*nkf2_cp_box;

      iap = ia + ia_shift;
      ibp = ib + ib_shift;
      icp = ic + ic_shift;

      iap = (iap > 0 ? iap : nkf1 + iap);
      ibp = (ibp > 0 ? ibp : nkf2 + ibp);
      icp = (icp > 0 ? icp : nkf3 + icp);

      iap = (iap > nkf1 ? iap - nkf1 : iap);
      ibp = (ibp > nkf2 ? ibp - nkf2 : ibp);
      icp = (icp > nkf3 ? icp - nkf3 : icp);

   index_lg = (iap) + (ibp - 1)*nkf1
            + (icp - 1)*nkf1*nkf2; 

   index_lg_cmplx = 2*index_lg - 1; 

   rfft[index_sm] = zfft[index_lg_cmplx] ;    

   }/*endfor ia*/
  }/*endfor ib*/
 }/*endfor ic*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_upack_rho_dual_ser_pme(double *zfft,double *rfft,
                                CELL *cell,int n_interp_pme_dual,
                                CPSCR *cpscr,
                                PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                                PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
 int i,iii,n2;

 int ja,jb,jc; 
 int ka,kb,kc;

 int ioff,joff;
 int nfft_cp_box2;

 double temp;

/* Local pointers */

      /* on small dense grid*/
 int  nfft_cp_box = para_fft_pkg3d_cp_box->nfft;  
 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;
 int  nkf3_cp_box = para_fft_pkg3d_cp_box->nkf3;

 int **igrid_a = cpscr->cpscr_dual_pme.igrid_a; /*length nkf1_cp_box*/
 int **igrid_b = cpscr->cpscr_dual_pme.igrid_b;
 int **igrid_c = cpscr->cpscr_dual_pme.igrid_c;

 double **mn_a  = cpscr->cpscr_dual_pme.mn_a;
                                   /*length: n_interp_pme_dual X nkf1_cp_box*/
 double **mn_b  = cpscr->cpscr_dual_pme.mn_b;
 double **mn_c  = cpscr->cpscr_dual_pme.mn_c;

/*==========================================================================*/
/* Zero rfft on the small grid  */

 nfft_cp_box2 = nfft_cp_box/2;
 for(i=1; i<= nfft_cp_box2; i++){
   rfft[i] = 0.0; 
 }/*endfor*/

/*--------------------------------------------------------------------------*/ 
/* A) Stuff Mn's on the grid:Sum is ordered to remove recursive dependencies*/ 
/*                           in qgrid. igrid_now is unique for fixed i      */ 

/* zfft is on large sparse grid in complex storage form*/
/* rfft holds v_ks_up(rd) on small dense grid in real form */

     for(jc=1;jc<=n_interp_pme_dual;jc++){
     for(kc=1; kc <= nkf3_cp_box ; kc++){

       for(jb=1;jb<=n_interp_pme_dual;jb++){
       for(kb=1; kb <= nkf2_cp_box ; kb++){

         ioff = (kb - 1)*nkf1_cp_box + (kc-1)*nkf1_cp_box*nkf2_cp_box;
         joff = igrid_b[jb][kb]+igrid_c[jc][kc];
         temp = mn_b[jb][kb]*mn_c[jc][kc];

         n2 = n_interp_pme_dual-3;
         for(ja=1; ja <=n2;ja+=4){
         for(ka=1; ka <= nkf1_cp_box ; ka++){
           rfft[(ka + ioff)]   += (
               (mn_a[ja    ][ka]*temp*zfft[(igrid_a[ja][ka] + joff)])
              +(mn_a[(ja+1)][ka]*temp*zfft[(igrid_a[(ja+1)][ka] + joff)])
              +(mn_a[(ja+2)][ka]*temp*zfft[(igrid_a[(ja+2)][ka] + joff)])
              +(mn_a[(ja+3)][ka]*temp*zfft[(igrid_a[(ja+3)][ka] + joff)]) 
                                   );
         }}/* endfor : ka,ja */

         n2=ja;
         for(ja=n2; ja <=n_interp_pme_dual;ja++){
         for(ka=1; ka <= nkf1_cp_box ; ka++){
           rfft[(ka + ioff)]   += 
               (mn_a[ja    ][ka]*temp*zfft[(igrid_a[ja][ka] + joff)]);
         }}/* endfor : ka,ja */

       }}/*endfor :  kb,jb */
     }}/*endfor : kc,jc*/


/*==========================================================================*/
/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Sngl pack the density from cp grid to large real space grid             */
/*==========================================================================*/

void sngl_upack_rho_dual_par_pme(CPSCR_WAVE *cpscr_wave,double *rfft,
                                CELL *cell,int n_interp_pme_dual,
                                CPSCR_DUAL_PME *cpscr_dual_pme,
                                PARA_FFT_PKG3D *para_fft_pkg3d_cp_box,
                                PARA_FFT_PKG3D *para_fft_pkg3d_lg)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */
#include "../typ_defs/typ_mask.h"

 int m,i,iii;
 int ja,jb,jc; 
 int ka,kb,kc;
 int irow;
 int joff;
 int nfft_proc_cp_box2;
 int icount;
 int kb_str_sm,kb_end_sm;
 double temp;

/* Local pointers */

 double *zfft      = cpscr_wave->zfft;
 double *zfft_tmp  = cpscr_wave->zfft_tmp;

 int  nkf1 = para_fft_pkg3d_lg->nkf1;
 int  nfft_proc_cp_box = para_fft_pkg3d_cp_box->nfft_proc; 

 /* on small dense grid*/
 int  nkf1_cp_box = para_fft_pkg3d_cp_box->nkf1;
 int  nkf2_cp_box = para_fft_pkg3d_cp_box->nkf2;

 /* PME MAP variables */
 int **igrid_a = cpscr_dual_pme->igrid_a; /*length nkf1_cp_box*/
 int **igrid_b = cpscr_dual_pme->igrid_b;
 int **igrid_c = cpscr_dual_pme->igrid_c;

 double **mn_a  = cpscr_dual_pme->mn_a;  
 double **mn_b  = cpscr_dual_pme->mn_b;
 double **mn_c  = cpscr_dual_pme->mn_c;

 int num_rows_tot           = cpscr_dual_pme->num_rows_tot;
 int **joff_sm_lg           = cpscr_dual_pme->joff_sm_lg;
 int *ioff_lg_sm            = cpscr_dual_pme->ioff_lg_sm;
 int nfft_recv_big_small    = cpscr_dual_pme->nfft_recv_big_small;
 int *send_counts_big_small = 
                         para_fft_pkg3d_cp_box->send_counts_row_big_small;
 int *recv_counts_big_small =
                        para_fft_pkg3d_cp_box->recv_counts_row_big_small;
 int *sdispls_big_small     = para_fft_pkg3d_cp_box->sdispls_row_big_small;
 int *rdispls_big_small     = para_fft_pkg3d_cp_box->rdispls_row_big_small;

 /* Small CP Box */
 int   skc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skc_fft_ka_proc;
 int   ekc_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekc_fft_ka_proc;
 int   skb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->skb_fft_ka_proc;
 int   ekb_fft_ka_proc_cp_box = para_fft_pkg3d_cp_box->ekb_fft_ka_proc;

 /*Large grid*/
 int nfft_proc        = para_fft_pkg3d_lg->nfft_proc;
 int nproc            = para_fft_pkg3d_lg->num_proc;

 MPI_Comm comm_lg     = para_fft_pkg3d_lg->comm; 
 int      myid        = para_fft_pkg3d_lg->myid;

/*=======================================================================*/
/* pack zfft into zfft_tmp  on the large sparse grid */


 for(i=1; i<= nfft_recv_big_small; i++){
   zfft_tmp[i] = 0.0;
 }

  icount = 0;  /*should end up being the receive counts */
  for(irow =1; irow <= num_rows_tot; irow++){
   for(ka=1,m=1; m<= nkf1; ka+=2,m++){
     zfft_tmp[m+icount] = zfft[ioff_lg_sm[irow] + ka];
   }/*endfor ka*/
   icount+= nkf1;
  }/*endfor*/


/*=======================================================================*/
/* Communicate zfft across processors to zfft_tmp */

  if(nproc > 1){
    Alltoallv(&(zfft_tmp[1]),&(recv_counts_big_small[1]),
              &(rdispls_big_small[1]),MPI_DOUBLE,
              &(zfft[1]),&(send_counts_big_small[1]),
              &(sdispls_big_small[1]),MPI_DOUBLE,comm_lg);
  }else{
    for(i=1; i<= nfft_proc; i++){zfft[i] = zfft_tmp[i];}
  }/*endif*/
   
/*=======================================================================*/
/* C) Interpolate the large grid to the small using the  Mn's weighting  */
/*    factors */

  nfft_proc_cp_box2 = nfft_proc_cp_box/2; 
  for(i=1; i<= nfft_proc_cp_box2; i++){
    rfft[i] = 0.0; 
  }/*endfor*/


  for(jc=1;jc<=n_interp_pme_dual;jc++){
  for(jb=1;jb<=n_interp_pme_dual;jb++){

    icount  = 0;
    for(kc=skc_fft_ka_proc_cp_box;kc<=ekc_fft_ka_proc_cp_box;kc++){
      kb_str_sm = (kc==skc_fft_ka_proc_cp_box ? skb_fft_ka_proc_cp_box : 1);
      kb_end_sm = (kc==ekc_fft_ka_proc_cp_box ? ekb_fft_ka_proc_cp_box 
                                              : nkf2_cp_box);

      for(kb=kb_str_sm;kb<=kb_end_sm;kb++){

       /* rfft =  rho*vol_cp, in the next section of cp_rho_calc,       */
       /* rho is divided by vol. Sinc vol = (box_rat)^3*vol_cp,         */
       /* you need to multiply rfft by scale = (box_rat)^3 here         */
       /* trap_corr is needed because we are changing the grid spacing  */
       /*  of the density and we still want the integral of the density */
       /* to generate the number of electrons   */

        temp = mn_b[jb][kb]*mn_c[jc][kc];

       /* All this goes to the same processor. That is, full rows of ka  */
       /* are stored on the same proc. joff says where this row should   */
       /* be stored in zfft so that the alltoallv will work properly.    */

         joff = joff_sm_lg[igrid_c[jc][kc]][igrid_b[jb][kb]];

         for(ja=1;ja<=n_interp_pme_dual;ja++){
         for(ka=1; ka <= nkf1_cp_box ; ka++){
          rfft[(ka+icount)] += 
                (mn_a[ja][ka]*temp*zfft[(igrid_a[ja][ka] + joff)]);
         }}/* endfor : ka,ja */
         icount += nkf1_cp_box;

      }/*endfor :  kb */
    }/*endfor :  kc */

  }}/*endfor : jb,jc*/

/*==========================================================================*/
/*-----------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/






