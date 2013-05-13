/*==================================  ===============================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                Module: control_vc_smpl.c                          */
/*                                                                   */
/* Control routine for velocity resampling                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vc_smpl(GENERAL_DATA *general_data,CP *cp)

/*===================================================================*/
      {/*begin routine*/
/*===================================================================*/
/*        Local Variables */
#include "../typ_defs/typ_mask.h"

  int debug_on,iii,is,icoef,i,ip;
  int iproj_vel=1;
  int pi_beads = cp->cpcoeffs_info.pi_beads_proc;
  int myid    = cp->communicate.myid;
  int myid_st = cp->communicate.myid_state;
  int np_states = cp->communicate.np_states;
  int nproc = cp->communicate.np;
  int cp_norb      = cp->cpopts.cp_norb;
  int cp_lsda      = cp->cpopts.cp_lsda;
  MPI_Comm world   = cp->communicate.world;
  MPI_Comm comm_states = cp->communicate.comm_states;

/*===================================================================*/
/* I) Notify screen what you are doing                               */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Sampling coefficient velocities\n");
    PRINT_LINE_DASH;putchar('\n');
  }/*endif*/
  if(np_states>1){  Barrier(comm_states);}

/*==================================================================*/
/* II) Checks */

  if(np_states>1){

    for(ip=1;ip<=pi_beads;ip++){
     if((cp->cpcoeffs_pos[ip].ivcoef_form_up
        +cp->cpcoeffs_pos[ip].icoef_form_up)!=2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coef velocities are not in transposed form \n");
      printf("on state processor %d in samp_vel_cp \n",myid_st);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
     if(cp_lsda==1){
      if((cp->cpcoeffs_pos[ip].ivcoef_form_up
         +cp->cpcoeffs_pos[ip].icoef_form_up)!=2){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn Coef velocities are not in transposed form \n");
       printf("on state processor %d in samp_vel_cp \n",myid_st);
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
      }/*endif:form check*/
     }/*endif:lsda*/
    }/*endfor:ip*/

  }/*endif:np_states*/

  if(cp_norb>0){

   for(ip=1;ip<=pi_beads;ip++){
    if(cp->cpcoeffs_pos[ip].icoef_orth_up!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Up Coefs must be in nonorthonormal form under norb \n");
      printf("in samp_vel_cp \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(cp_lsda==1){
     if(cp->cpcoeffs_pos[ip].icoef_orth_dn!=0){
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Dn Coefs must be in nonorthonormal form under norb \n");
       printf("in samp_vel_cp \n");
       printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
     }/*endif:form check*/
    }/*endif:lsda*/
   }/*endfor:ip*/

  }/*endif:cp_norb*/

/*===================================================================*/
/* III) Sample velocities                                             */

  printf("1\n");

  sampl_vc(&(cp->cpopts),
           &(cp->cpcoeffs_info),(cp->cpcoeffs_pos),
	   &(general_data->statepoint),&(cp->vel_samp_cp.iseed),
	   &(cp->vel_samp_cp.iseed2),&(cp->vel_samp_cp.qseed), 
           &(cp->communicate),&(cp->cp_comm_state_pkg_up),
           &(cp->cp_comm_state_pkg_dn));

  printf("2\n");
/*===================================================================*/
/* IV) Project onto surface of constraint                           */

  debug_on = general_data->simopts.debug + general_data->simopts.debug_pimd
           + general_data->simopts.debug_cp 
           + general_data->simopts.debug_cp_pimd;

  if(cp->cpopts.cp_norb<3){

    if(myid==0){
     if(debug_on == 1&&myid==0) {
      printf("Do you wish to project coefficient velocities? (1 or 0)\n");
      scanf("%d",&iproj_vel);
     }/*endif*/
    }/*endif*/
    Barrier(comm_states);
    if(nproc>1){Bcast(&iproj_vel,1,MPI_INT,0,world);}

    if((debug_on == 1 && iproj_vel == 1) ||
       (debug_on != 1) && (cp->cpopts.zero_cp_vel < 1)){
      for(ip=1;ip<=pi_beads;ip++){ 
        proj_velc(cp,ip);
      }/*endfor*/
    }/*endif*/

  }/*endif*/

  printf("3\n");
/*===================================================================*/
/* IV) Reset zero_cp_vel flag if necessary                           */

 cp->cpopts.zero_cp_vel = (cp->cpopts.zero_cp_vel <= 1 ? 0:2);

/*===================================================================*/
/* V) Project onto screen and bail                                  */

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Coefficient velocity sampling complete\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/
  Barrier(comm_states);

/*===================================================================*/
    }/*end routine*/
/*===================================================================*/





