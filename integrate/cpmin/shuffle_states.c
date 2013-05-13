
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: shuffle_states                               */
/*                                                                          */
/* This subprogram shuffles the Kohn-Sham states during minimization        */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_integrate_cpmin_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cpcon_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void cp_shuffle_states(CP *cp,int ip)

/*==========================================================================*/
{/* Begin routine */
/*========================================================================*/
/*             Local variable declarations                                */

 int cp_lsda                  = cp->cpopts.cp_lsda;
 int nstate_up                = cp->cpcoeffs_info.nstate_up;
 int *ioff_upt                = cp->cpcoeffs_info.ioff_upt;
 int nstate_dn                = cp->cpcoeffs_info.nstate_dn;
 int *ioff_dnt                = cp->cpcoeffs_info.ioff_dnt;

 int    icoef_form_up         = cp->cpcoeffs_pos[ip].icoef_form_up;
 double *cpcoeffs_cre_up      = cp->cpcoeffs_pos[ip].cre_up;
 double *cpcoeffs_cim_up      = cp->cpcoeffs_pos[ip].cim_up;
 int    icoef_form_dn         = cp->cpcoeffs_pos[ip].icoef_form_dn;
 double *cpcoeffs_cre_dn      = cp->cpcoeffs_pos[ip].cre_dn;
 double *cpcoeffs_cim_dn      = cp->cpcoeffs_pos[ip].cim_dn;
 double *cpscr_cre_up         = cp->cpscr.cpscr_wave.cre_up;
 double *cpscr_cim_up         = cp->cpscr.cpscr_wave.cim_up;
 double *cpscr_cre_dn         = cp->cpscr.cpscr_wave.cre_dn;
 double *cpscr_cim_dn         = cp->cpscr.cpscr_wave.cim_dn;

 int np_states                = cp->communicate.np_states;
 int myid_state               = cp->communicate.myid_state;

 double *occ_up               = cp->cpopts.occ_up;
 double *occ_dn               = cp->cpopts.occ_dn;
 double *rocc_sum_up          = cp->cpopts.rocc_sum_up;
 double *rocc_sum_dn          = cp->cpopts.rocc_sum_dn;

/*==========================================================================*/
/* 0) Check  */

 if(np_states>1){
  if(icoef_form_up !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Coefficients are not in transposed form \n");
     printf("on state processor %d in cp_shuffle_states \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/
 }/*endif*/

 if(np_states>1 && cp_lsda==1){
  if(icoef_form_dn !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Coefficients are not in transposed form \n");
     printf("on state processor %d in cp_shuffle_states \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/
 }/*endif*/

/*==========================================================================*/
/* I) UP STATES */

  if( nstate_up > 1 ){
   cp_shuffle_prim(cpcoeffs_cre_up,cpcoeffs_cim_up,icoef_form_up,
                   cpscr_cre_up,cpscr_cim_up,ioff_upt,
                   occ_up,rocc_sum_up,occ_dn,cp_lsda,
                   &(cp->vel_samp_cp),&(cp->cp_comm_state_pkg_up));
  }/*endif*/

/*==========================================================================*/
/* II) DN STATES */

  if( nstate_dn > 1 && cp_lsda == 1){
   cp_shuffle_prim(cpcoeffs_cre_dn,cpcoeffs_cim_dn,icoef_form_dn,
                   cpscr_cre_dn,cpscr_cim_dn,ioff_dnt,
                   occ_dn,rocc_sum_dn,occ_up,cp_lsda,
                   &(cp->vel_samp_cp),&(cp->cp_comm_state_pkg_dn));
  }/*endif*/

/*==========================================================================*/
    }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void cp_shuffle_prim(double *cpcoefs_cre,double *cpcoefs_cim,int icoef_form,
                     double *cpscr_cre,double *cpscr_cim,int *iofft,
                     double *occ,double *rocc_sum,double *occ_dn, int cp_lsda,
                     VEL_SAMP_CP *vel_samp_cp,
                     CP_COMM_STATE_PKG *cp_comm_state_pkg)

/*==========================================================================*/
{/* Begin routine */
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int rnum,count,i,is,ncoef_tot,got_it,*map;
 int icoef1,icoef2,iii;
 int iocc,j; 
 double *temp;

/*           Local Pointers */
 int np_states             = cp_comm_state_pkg->num_proc;
 int myid_state            = cp_comm_state_pkg->myid;
 int nstate                = cp_comm_state_pkg->nstate;
 int ncoef                 = cp_comm_state_pkg->nstate_ncoef_proc_max;
 MPI_Comm comm_state       = cp_comm_state_pkg->comm;


/*==========================================================================*/
/* 0) Parallel Checks */

 if(np_states>1){
  if(icoef_form !=1){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Coefficients are not in transposed form \n");
     printf("on state processor %d in cp_shuffle_prim \n",myid_state);
     printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/
 }/*endif*/

/*==========================================================================*/
/* I) malloc the map */

 map = (int *) cmalloc(nstate*sizeof(int))-1;
 temp = (double *) cmalloc(nstate*sizeof(double ))-1;

/*==========================================================================*/
/* II) Processor 0 makes the map and broadcasts it                           */

 if(myid_state == 0){ 

    rnum = (int) (((double) nstate)*ran_essl(&(vel_samp_cp->qseed)))+1;
    map[1] = rnum;

   count=1;
   do {
    do {
      rnum = (int) (((double) nstate)*ran_essl(&(vel_samp_cp->qseed)))+1;
      for(i=1; i<= count;i++){
       if(rnum == map[i]) {got_it=0;break;} else got_it=1;
      } 
    } while(!got_it);
   ++count;
   map[count] = rnum;
   } while(count < nstate);

 }/* endif myid == 0 */

 if(np_states>1){
   Bcast(&(map[1]),nstate,MPI_INT,0,comm_state);
 }/*endif*/

/*==========================================================================*/
/* III) Everybody shuffle the states according to map                       */

   ncoef_tot = ncoef*nstate;
   for(i=1;i<=ncoef_tot;i++){
    cpscr_cre[i] = cpcoefs_cre[i];
    cpscr_cim[i] = cpcoefs_cim[i];
   }/*endfor*/

   for(is=1;is<=nstate;is++) {
     for(i=1;i<=ncoef;i++) {
       icoef1 = iofft[is]  + i;
       icoef2 = iofft[map[is]] + i;
       cpcoefs_cre[icoef1] = cpscr_cre[icoef2];
       cpcoefs_cim[icoef1] = cpscr_cim[icoef2];
     }/*endfor*/
   }/*endfor*/


/* Shuffle the occupation numbers  according to map                           */
  for(i=1; i<= nstate; i++){
     temp[i] =  occ[i];
  }/*endfor*/

  for(i=1; i<= nstate; i++){
     occ[i] = temp[map[i]];
  }/*endfor*/

  if(cp_lsda == 0){

    for(i=1; i<= nstate; i++){
      temp[i] =  occ_dn[i];
    }/*endfor*/

    for(i=1; i<= nstate; i++){
     occ_dn[i] = temp[map[i]];
    }/*endfor*/

  }/*endif*/
 
/*  recalculate rocc_sum    */
  iocc=0;
  for(i=1;i<= nstate ;i++){
   for(j=1;j<= nstate ;j++){
     iocc++;
    rocc_sum[iocc] = 1.0/(occ[i]+occ[j]);
   }/*endfor i*/
  }/* endfor j*/



/*==========================================================================*/
/* IV) Free the map                                                         */

   cfree(&(map[1]));
   cfree(&(temp[1]));

/*==========================================================================*/
    }/* end routine */
/*==========================================================================*/








