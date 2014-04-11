/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: gen_wave.c                                   */
/*                                                                          */
/* Construct an initial wave function using the radial states               */
/* of the isolated atoms                                                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_cp_ewald_local.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_handle_entry.h"

#define DEBUG_GW_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


void gen_wave(CLASS *class,GENERAL_DATA *general_data,CP *cp,
              CP_PARSE *cp_parse,NAME *vps_name)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

/* misc                                                             */
  double conv;
  int nsplin = cp->pseudo.nsplin_g;

  int i,j,k,iatm,ipart,is,iii;
  int ind,ind1,ind2,ibinary,n;

/*----------------------------------------------------------------------*/
/* particalus positionus */
  double *x,*y,*z;

  int    nab_initio = class->clatoms_info.nab_initio;

  int   *cp_atm_lst   = class->atommaps.cp_atm_lst;
  int   *cp_vlnc_up   = class->clatoms_info.cp_vlnc_up;
  int   *cp_vlnc_dn   = class->clatoms_info.cp_vlnc_dn;

/*----------------------------------------------------------------------*/
/* dual grid opt varialbes */

  int cp_dual_grid_opt_on   = cp->cpopts.cp_dual_grid_opt;
  double *cp_box_center     =  general_data->cell.cp_box_center;
  double *cp_box_center_rel =  general_data->cell.cp_box_center_rel;
  double dx,dy,dz;
  double xtrans,ytrans,ztrans;
  double asx,asy,asz;
  double sx,sy,sz;

/*----------------------------------------------------------------------*/
/* hmatrikus */

  double *hmat     = general_data->cell.hmat_cp;
  double *hmati    = general_data->cell.hmati_cp;
  double *hmat_big = general_data->cell.hmat;
  double *hmati_big= general_data->cell.hmati;
  double  vol      = general_data->cell.vol_cp;

/*----------------------------------------------------------------------*/
/* oribitalus coefficientus                                             */

  double *creal_up   = cp->cpcoeffs_pos[1].cre_up;
  double *cimag_up   = cp->cpcoeffs_pos[1].cim_up;
  double *creal_dn   = cp->cpcoeffs_pos[1].cre_dn;
  double *cimag_dn   = cp->cpcoeffs_pos[1].cim_dn;
  double *creal_tmp  = cp->cpscr.cpscr_wave.cre_up;
  double *cimag_tmp  = cp->cpscr.cpscr_wave.cim_up;

  double p_a,p_b,p_c;

/*----------------------------------------------------------------------*/
/* Up and Dn state variables  */
  int istate_up,istate_dn,istate;
  int istate_up_proc,istate_dn_proc;
  int nstate_up_gw,nstate_dn_gw;
  int  *nstate_up_atm,*nstate_dn_atm;
  int nstate_up       = cp->cpcoeffs_info.nstate_up;
  int nstate_dn       = cp->cpcoeffs_info.nstate_dn; 
  int nstate_up_proc  = cp->cp_comm_state_pkg_up.nstate_proc;
  int nstate_dn_proc  = cp->cp_comm_state_pkg_dn.nstate_proc;
  int istate_up_st    = cp->cpcoeffs_info.istate_up_st;
  int istate_up_end   = cp->cpcoeffs_info.istate_up_end;
  int istate_dn_st    = cp->cpcoeffs_info.istate_dn_st;
  int istate_dn_end   = cp->cpcoeffs_info.istate_dn_end;

/* number of coeff,maxcoef,nstate_up,nstate_dn                    */
   int  ncoef    =  cp->cp_comm_state_pkg_up.ncoef;
   int  ncoef1;

/*----------------------------------------------------------------------*/
/*spline of radial wave functions                                  */

  double gmin  = cp->cpewald.gw_gmin;
  double gmax  = cp->cpewald.gw_gmax;
  double ***gpsi0,***gpsi1, ***gpsi2, ***gpsi3;
  double **gpsi_now;                             /*(maxatmtyp,maxatmtyp)*/
  double *gpsi00;                                /*(maxatmtyp)*/
  double *psi_r,*psi_i; /*length(20)*/

/*----------------------------------------------------------------------*/
/* Gram-Schmidt and Kinetic energy Variables*/

  int icoef_form_up,icoef_form_dn;
  double *ovlap       = cp->cpcoeffs_pos[1].ovmat_eigv_up;
  double *ovlap_tmp   = cp->cpcoeffs_pos[1].ovmat_eigv_dn;
  int *ioff_upt  =  cp->cpcoeffs_info.ioff_upt;
  int *ioff_dnt  =  cp->cpcoeffs_info.ioff_dnt;

  double o,eke,eke_tot;
  double eke_tot_proc;


/* Variables to create g vectors */
  double dg,g,volrt;
  double pi,tpi,aka,akb,akc,fpi;
  double xk,yk,zk,g2;
  double rad2;
  double atemp,btemp,ctemp,arg,helr,heli;
  int itemp;
  double qseed=103481.0;

  double kappa = 1.0;

/* map from spherically cutoff half space to double full space     */
   int *kastore   = cp->cpewald.kastr_sm;
   int *kbstore   = cp->cpewald.kbstr_sm;
   int *kcstore   = cp->cpewald.kcstr_sm; 
   int nktot_sm   = cp->cpewald.nktot_sm;  /*length nktot for small grid */

/*----------------------------------------------------------------------*/
/* ylm constants                                                        */
   YLM_CONS ylm_cons;
   double *ylmr, *ylmi;
   double *dylmr_x,*dylmi_x;
   double *dylmr_y,*dylmi_y;
   double *dylmr_z,*dylmi_z;
   int *n_ang; /* length natm_typ_cp*/
/*----------------------------------------------------------------------*/

  FILE *fp_wave_in;
  NAME *fname_ps;

/*----------------------------------------------------------------------*/
/* Communicate Variables */
  int master = 0;
  int myid              = cp->communicate.myid;
  int nproc             = cp->communicate.np;
  MPI_Comm comm_states  = cp->communicate.comm_states;


  int cp_lsda         = cp->cpopts.cp_lsda;
  int cp_lda          = cp->cpopts.cp_lda;

/*----------------------------------------------------------------------*/
/* Occupation Number Variables */

  double *occ_dn      = cp->cpopts.occ_dn;
  double *occ_up      = cp->cpopts.occ_up;
  double *rocc_sum_up = cp->cpopts.rocc_sum_up;
  double *rocc_sum_dn = cp->cpopts.rocc_sum_dn;
  int iocc;

/*----------------------------------------------------------------------*/
/* Atom label Variables */

  int natm_typ_cp;
  int iflag;
  int tag;

  int *iatm_mol_typ   = class->atommaps.iatm_mol_typ;
  int *iatm_atm_typ   = class->atommaps.iatm_atm_typ;
  int natm_typ        = class->atommaps.natm_typ;

  int natm_typ_gw;
  int *iatm_atm_typ_cp;
  
/*====================================================================*/

   if(myid == master ){
     printf("\n-----------------------------------------\n");
     printf("Constructing Wave Function \n");
     printf("-----------------------------------------\n");
   }/*endif*/

/*=======================================================================*/
/*Check Sum of cp_valences is equal to the number of states in set file  */
/*=======================================================================*/


   if( nab_initio == 0){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There are no ab initio atoms \n");
     printf("Please check whether you should have assigned some atoms \n");
     printf("to be ab initio.  To proceed would be pointless.\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);
  }/*endif*/


   nstate_up_gw = 0;
   nstate_dn_gw = 0;
  for(i=1; i<= nab_initio; i++){
   iatm = cp_atm_lst[i];
   nstate_up_gw += cp_vlnc_up[iatm];
   nstate_dn_gw += cp_vlnc_dn[iatm];
  }/*endfor*/


  if( nstate_up_gw != nstate_up ){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The number of states up not equal to what is in the set file\n");
     printf("%d here and %d from the set file\n",nstate_up_gw, nstate_up);
     printf("Please check the values of cp_valence_up in the parm files\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);
  }/*endif*/

  if(  nstate_dn_gw != nstate_dn  && cp_lsda==1){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The number of states dn not equal to what is in the set file\n");
     printf("%d here and %d from the set file\n",nstate_dn_gw, nstate_dn);
     printf("Please check the values of cp_valence_dn in the parm files\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);
  }/*endif*/

  if(  nstate_dn_gw != nstate_dn  && cp_lda==1){
    if(myid == master){
     printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
     printf("The number of states dn not equal to value in the set file.\n");
     printf("This is OK only if you have some unit occupation numbers. \n");
     printf("Dawn will liberate occupation numbers from states soon!\n");
     printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif myid*/
  }/*endif*/

  if(  nstate_dn_gw > nstate_dn  && cp_lda==1){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("The number of dn nstates is assumed <= up nstates \n");
     printf("in gen_wave under LDA. Since you already have the warning\n" );
     printf("about occupation numbers and states what you have to do is\n" );
     printf("check the values of cp_valence_(dn/up) in the parm files\n");
     printf("to ensure that the up states get the extra occupation \n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);
  }/*endif*/


  /*check that the value assigned for cp_valence_up is same for */
  /* for all atoms of the same atm_typ                          */

  for(i=1; i <= nab_initio; i++){
    ind1 = cp_atm_lst[i];
    for(j=i+1; j <= nab_initio;j++){
     ind2 = cp_atm_lst[j];
     if( iatm_atm_typ[ind2] == iatm_atm_typ[ind1] &&
         cp_vlnc_up[ind2] != cp_vlnc_up[ind1]){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There are two different assigned values for    \n");
     printf("cp_valence_up for the same atom type           \n");
     printf("atom numbers %d and %d \n",ind1,ind2);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);

     }/*endif*/

     if( cp_lsda == 1 && iatm_atm_typ[ind2] == iatm_atm_typ[ind1] &&
         cp_vlnc_dn[ind2] != cp_vlnc_dn[ind1]){
    if(myid == master){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("There are two different assigned values for    \n");
     printf("cp_valence_dn for the same atom type           \n");
     printf("atom numbers %d and %d \n",ind1,ind2);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif myid*/
    if(nproc>1){
     Barrier(comm_states);
     Finalize();
    }
     exit(1);

     }/*endif*/
    }/*endfor*/
  }/*endfor*/


/*=======================================================================*/
/* Create number of ab initio atom types natm_typ_cp and make list       */
/*=======================================================================*/

    natm_typ_cp = 1;

   for(i=2; i<= nab_initio; i++){
     iflag = 0;
     ind1 = cp_atm_lst[i];
    for(j=i-1; j >= 1; j--){
     ind2 = cp_atm_lst[j];
     if( iatm_atm_typ[ind2] == iatm_atm_typ[ind1] ) {
       iflag = 0;break;
     }else{
      iflag = 1;
     }/*endif*/
    }/*endfor*/
     natm_typ_cp += iflag;
   }/*endfor*/

/*-----------------------------------------------------------------------*/
/* malloc list iatm_atm_typ_cp holds index for what cp_atm_typ this is   */
/* malloc fname_ps length natm_typ_cp  holds names of pseud files        */

   iatm_atm_typ_cp = (int *) cmalloc(nab_initio*sizeof(int)) -1;

   if(myid == 0 ){
    fname_ps = (NAME *) cmalloc(natm_typ_cp*sizeof(NAME))-1;
   }

   nstate_up_atm = (int *) cmalloc(natm_typ_cp*sizeof(int )) -1;
   nstate_dn_atm = (int *) cmalloc(natm_typ_cp*sizeof(int )) -1;

/*-----------------------------------------------------------------------*/
/* Assign cp_atm_types  and values for nstate_up_atm nstate_dn_atm       */

   /*Assign first element in arrays */
  iatm_atm_typ_cp[1] = 1;
  if(myid == 0){
   strcpy(fname_ps[1],vps_name[iatm_atm_typ[cp_atm_lst[1]]]);
  }/*endif myid*/
  nstate_up_atm[1] = cp_vlnc_up[cp_atm_lst[1]];
  nstate_dn_atm[1] = cp_vlnc_dn[cp_atm_lst[1]];

  tag = 1;  /* used as counter for unique cp atom types */

  for(i=2; i<= nab_initio; i++){
    iflag = 0;
    ind1 = cp_atm_lst[i];
   for(j=i-1; j >= 1; j--){
    ind2 = cp_atm_lst[j];
      if( iatm_atm_typ[ind1] == iatm_atm_typ[ind2] ){
        iatm_atm_typ_cp[i] = iatm_atm_typ_cp[j];
        iflag = 0;
        break;
      }else{
        iflag = 1;
      }/*endif*/
   }/*endfor*/
    tag += iflag;
   if(iflag == 1){
      iatm_atm_typ_cp[i] = tag;
     if(myid == 0){
      strcpy(fname_ps[tag],vps_name[iatm_atm_typ[ind1]]);
     }/*endif myid*/
      nstate_up_atm[tag] = cp_vlnc_up[ind1];
      nstate_dn_atm[tag] = cp_vlnc_dn[ind1];
   }/*endif*/
  }/*endfor*/

#ifdef DEBUG_GW
 if(myid == 0){
  for(i=1; i<= nab_initio; i++){
    ind1 = cp_atm_lst[i];
    printf("i %d atm_typ %d atm_typ_cp %d \n",i,iatm_atm_typ[ind1],
            iatm_atm_typ_cp[i]);
  }

  for(i=1; i<= natm_typ_cp; i++){
    printf("i %d name %s vlnc up %d dn %d \n",i,fname_ps[i],
             nstate_up_atm[i],nstate_dn_atm[i]);
  }

 }
   Dbx_Barrier(class->communicate.world);
#endif

/*=======================================================================*/
/* Assign positions array length of only number of ab initio atoms       */
/*=======================================================================*/

   x = (double *) cmalloc(nab_initio*sizeof(double )) -1;
   y = (double *) cmalloc(nab_initio*sizeof(double )) -1;
   z = (double *) cmalloc(nab_initio*sizeof(double )) -1;

  for(i=1; i<=nab_initio; i++){
   iatm = cp_atm_lst[i];
   x[i] = class->clatoms_pos[1].x[iatm];
   y[i] = class->clatoms_pos[1].y[iatm];
   z[i] = class->clatoms_pos[1].z[iatm];
  }/*endfor*/

  if(nproc > 1){
   Bcast(&(x[1]),nab_initio,MPI_DOUBLE,master,comm_states);
   Bcast(&(y[1]),nab_initio,MPI_DOUBLE,master,comm_states);
   Bcast(&(z[1]),nab_initio,MPI_DOUBLE,master,comm_states);
  }
 
#ifdef DEBUG_GW
  for(i=1; i<= nab_initio; i++){
    printf("myid %d x y z %lg %lg %lg cp_atm_typ %d \n",myid,
            x[i],y[i],z[i],iatm_atm_typ_cp[i]);
  }
 
   Dbx_Barrier(class->communicate.world);  

  if(myid == 0){
   for(i=1; i<= nab_initio; i++){
      iii = iatm_atm_typ[cp_atm_lst[i]];
      printf("i %d file name %s %s \n",i,vps_name[iii],fname_ps[iatm_atm_typ_cp[i]]);
   }   
  }
   Dbx_Barrier(class->communicate.world);

#endif

   xtrans = cp_box_center[1] - cp_box_center_rel[1];
   ytrans = cp_box_center[2] - cp_box_center_rel[2];
   ztrans = cp_box_center[3] - cp_box_center_rel[3];


/*=========================================================================*/
/*  malloc some local memory using parameters given in initial input file  */

   n_ang         = (int *) cmalloc(natm_typ_cp*sizeof(int )) -1;

   ylmr  = (double *) cmalloc(20*sizeof(double)) -1;
   ylmi  = (double *) cmalloc(20*sizeof(double)) -1;

   dylmr_x  = (double *) cmalloc(16*sizeof(double)) -1;
   dylmi_x  = (double *) cmalloc(16*sizeof(double)) -1;
   dylmr_y  = (double *) cmalloc(16*sizeof(double)) -1;
   dylmi_y  = (double *) cmalloc(16*sizeof(double)) -1;
   dylmr_z  = (double *) cmalloc(16*sizeof(double)) -1;
   dylmi_z  = (double *) cmalloc(16*sizeof(double)) -1;

   psi_r    = (double *) cmalloc(20*sizeof(double )) -1;
   psi_i    = (double *) cmalloc(20*sizeof(double )) -1;

/*===========================================================================*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*        do NOT change this 3d malloc to the cmall_tens3                    */
/*           UNLESS you want to BREAK this code                              */
/*!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*         nsplin rows 3 columns natm_typ_cp dimensions                      */

   gpsi0 = (double ***) cmalloc(natm_typ_cp*sizeof(double **))-1;
   gpsi1 = (double ***) cmalloc(natm_typ_cp*sizeof(double **))-1;
   gpsi2 = (double ***) cmalloc(natm_typ_cp*sizeof(double **))-1;
   gpsi3 = (double ***) cmalloc(natm_typ_cp*sizeof(double **))-1;

 for(i=1; i<= natm_typ_cp; i++){
   gpsi0[i] = (double **) cmalloc(4*sizeof(double *))-1;
   gpsi1[i] = (double **) cmalloc(4*sizeof(double *))-1;
   gpsi2[i] = (double **) cmalloc(4*sizeof(double *))-1;
   gpsi3[i] = (double **) cmalloc(4*sizeof(double *))-1;
   for(j=1; j<=4; j++){
      gpsi0[i][j] = (double *) cmalloc(nsplin*sizeof(double ))-1;
      gpsi1[i][j] = (double *) cmalloc(nsplin*sizeof(double ))-1;
      gpsi2[i][j] = (double *) cmalloc(nsplin*sizeof(double ))-1;
      gpsi3[i][j] = (double *) cmalloc(nsplin*sizeof(double ))-1;
   }/*endfor*/
 }/*endfor*/

   gpsi_now = cmall_mat(1,natm_typ_cp,1,4);
   gpsi00   = (double *) cmalloc(natm_typ_cp*sizeof(double)) -1;

/*===========================================================================*/
/* assign the occupation numbers  and rocc_sum matrices                      */
/*===========================================================================*/

   for(i=1; i<= nstate_up; i++){
     occ_up[i] = 0.0;
     occ_dn[i] = 0.0;
   }/*endfor*/

   for(i=1; i<= nstate_up_gw; i++){
     occ_up[i] = 1.0;
   }

   for(i=1; i<= nstate_dn_gw; i++){
     occ_dn[i] = 1.0;
   }

   if(cp_lda==1){
    for(i=1; i<= nstate_up; i++){
     occ_up[i] += occ_dn[i];
    }/*endfor*/
   }/*endif cp_lda*/

  if(nproc>1){Barrier(comm_states);}

   iocc=0;
   for(i=1;i<=nstate_up;i++){
    for(j=1;j<=nstate_up;j++){
     iocc++;
     rocc_sum_up[iocc] = 1.0/(occ_up[i]+occ_up[j]);
    }/*endfor i*/
   }/* endfor j*/

   if(cp_lsda==1){
    iocc=0;
    for(i=1;i<=nstate_dn;i++){
     for(j=1;j<=nstate_dn;j++){
      iocc++;
      rocc_sum_dn[iocc] = 1.0/(occ_dn[i]+occ_dn[j]);
     }/*endfor i*/
    }/* endfor j*/
   }/*endif*/

/* =========================================================================*/
/* Assign useful constants                                                  */
/* =========================================================================*/

   volrt = sqrt(vol);

   pi         = M_PI;
   tpi        = 2.0*pi;
   fpi        = 4.0*pi;
   ylm_cons.rt_fpi     = 1.0/sqrt(fpi);
   ylm_cons.rt_thrfpi  = sqrt(3.0/fpi);
   ylm_cons.rt_threpi  = sqrt(1.50/fpi);
   ylm_cons.hrt_fivfpi = 0.50*sqrt(5.0/fpi);
   ylm_cons.rt_fiftepi = sqrt(7.50/fpi);
   ylm_cons.hrt_sevfpi = 0.50*sqrt(7.0/fpi);
   ylm_cons.hrt_toepi  = 0.50*sqrt(10.50/fpi);
   ylm_cons.hrt_ohffpi = 0.50*sqrt(105.0/fpi);
   ylm_cons.hrt_tfepi  = 0.50*sqrt(17.50/fpi);
   rad2       = sqrt(2.0);

/*=========================================================================*/
/* bessel transform the radial wave functions and spline the result        */
/* nsplin rows 3 columns natm_typ dimensions  [d][c][r]                    */
/*=========================================================================*/

 for(i=1; i<= natm_typ_cp; i++){
   splin_btrans(&nsplin,gpsi0[i],gpsi1[i],
                gpsi2[i],gpsi3[i],gmin,gmax,&dg,
                &(gpsi00[i]),&(n_ang[i]),fname_ps[i],i,
                myid,nproc,comm_states);
 }/*endfor*/

/*=========================================================================*/
/* get the wave functions in g space                                       */
/*=========================================================================*/

 for(i=1; i<=  ncoef-1; i++){

/*-------------------------------------------------------------------------*/
/* get g vectors                                                           */
   aka = -(double)(kastore[i]);
   akb = -(double)(kbstore[i]);
   akc = -(double)(kcstore[i]);

   xk = (aka*hmati[1] + akb*hmati[2] +akc*hmati[3])*tpi;
   yk = (aka*hmati[4] + akb*hmati[5] +akc*hmati[6])*tpi;
   zk = (aka*hmati[7] + akb*hmati[8] +akc*hmati[9])*tpi;
   g2 = (xk*xk+yk*yk+zk*zk);
   g  = sqrt(g2);
  
/*-------------------------------------------------------------------------*/
/*  get ylm                                                                */

 get_ylm(xk,yk,zk,g,ylmr,ylmi,
         dylmr_x,dylmi_x,dylmr_y,dylmi_y,dylmr_z,dylmi_z,&ylm_cons);

/*-------------------------------------------------------------------------*/
/*  get the spherical bessel tranform of the radial wave functions         */
/*  at this g space point using the spline.                                */

   for(iatm=1; iatm <= natm_typ_cp; iatm++){
    get_gpsi(g,nsplin,
             gpsi0[iatm],gpsi1[iatm],gpsi2[iatm],gpsi3[iatm],
             gpsi_now[iatm],gmin,dg,n_ang[iatm]);
   }/*endfor*/

  istate_up = 0;
  istate_dn = 0;

  istate_up_proc = 0;
  istate_dn_proc = 0;

  for(ipart = 1; ipart <= nab_initio; ipart++){
/*  S STATE                                                                 */
     psi_r[1] = ylmr[1]*gpsi_now[iatm_atm_typ_cp[ipart]][1]/volrt; 
    /*     psi_r[1] = exp(-g2/(2.0*kappa*kappa))/volrt; */
     psi_i[1] = 0.0;
/* SPHERICALIZED P BAND                                                    */
    if(n_ang[iatm_atm_typ_cp[ipart]] >= 1){
      itemp = (int) (3.0*ran_essl(&qseed));
      itemp = MAX(itemp,0);
      itemp = MIN(itemp,2);
      switch(itemp){
       case 0:
        p_c = -ylmr[2]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_b = -rad2*ylmr[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_a = -rad2*ylmi[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
       break;
       case 1:
        p_a = -ylmr[2]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_c = -rad2*ylmr[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_b = -rad2*ylmi[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
       break;
       case 2:
        p_b = -ylmr[2]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_a = -rad2*ylmr[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
        p_c = -rad2*ylmi[3]*gpsi_now[iatm_atm_typ_cp[ipart]][2]/volrt;
       break;
      }/*end switch*/
      psi_r[2] = 0.0;
      psi_i[2] = (p_c + (p_b+p_a)/rad2)/rad2;
      psi_r[3] = 0.0;
      psi_i[3] = (p_c - (p_b+p_a)/rad2)/rad2;
      psi_r[4] = 0.0;
      psi_i[4] = (p_b-p_a)/rad2;
    }/*endif*/
/*  D BAND                                                                   */
   if(n_ang[iatm_atm_typ_cp[ipart]] >= 2){
     psi_r[5] = -ylmr[5]*gpsi_now[iatm_atm_typ_cp[ipart]][3]/volrt;
     psi_i[5] = 0.0;
     psi_r[6] = -rad2*ylmr[6]*gpsi_now[iatm_atm_typ_cp[ipart]][3]/volrt;
     psi_i[6] = 0.0;
     psi_r[7] = -rad2*ylmi[6]*gpsi_now[iatm_atm_typ_cp[ipart]][3]/volrt;
     psi_i[7] = 0.0;
     psi_r[8] = -rad2*ylmr[8]*gpsi_now[iatm_atm_typ_cp[ipart]][3]/volrt;
     psi_i[8] = 0.0;
     psi_r[9] = -rad2*ylmi[8]*gpsi_now[iatm_atm_typ_cp[ipart]][3]/volrt;
     psi_i[9] = 0.0;
   }/*endif*/

/*---------------------------------------------------------------------------*/
/*  structure factor                                                         */

    dx  = x[ipart] - cp_box_center[1];
    dy  = y[ipart] - cp_box_center[2];
    dz  = z[ipart] - cp_box_center[3];
    asx = dx*hmati_big[1]+dy*hmati_big[4]+dz*hmati_big[7];
    asy = dx*hmati_big[2]+dy*hmati_big[5]+dz*hmati_big[8];
    asz = dx*hmati_big[3]+dy*hmati_big[6]+dz*hmati_big[9];

    sx  = asx - NINT(asx);
    sy  = asy - NINT(asy);
    sz  = asz - NINT(asz);
    dx  = sx*hmat_big[1]+sy*hmat_big[4]+sz*hmat_big[7] + cp_box_center_rel[1];
    dy  = sx*hmat_big[2]+sy*hmat_big[5]+sz*hmat_big[8] + cp_box_center_rel[2];
    dz  = sx*hmat_big[3]+sy*hmat_big[6]+sz*hmat_big[9] + cp_box_center_rel[3];

    atemp = dx*hmati[1] + dy*hmati[4] + dz*hmati[7];
    btemp = dx*hmati[2] + dy*hmati[5] + dz*hmati[8];
    ctemp = dx*hmati[3] + dy*hmati[6] + dz*hmati[9];

    arg = (aka*atemp + akb*btemp + akc*ctemp)*tpi;

    helr = cos(arg);
    heli = sin(arg);
 
/*=========================================================================*/
/*  construct the coeff                                                    */

   for(is=1; is <= nstate_up_atm[iatm_atm_typ_cp[ipart]] ; is++){
   if( ((istate_up + is ) >= istate_up_st) &&
       ((istate_up + is) <= istate_up_end)){  
      ind = i + (istate_up - istate_up_st + is  )*ncoef;
     creal_up[ind]  =  helr*psi_r[is] - heli*psi_i[is];
     cimag_up[ind]  =  heli*psi_r[is] + helr*psi_i[is];
    }/*endif*/
   }/*endfor*/

   istate_up +=  nstate_up_atm[iatm_atm_typ_cp[ipart]];

  if((cp_lsda == 1) && (nstate_dn > 0)){
   for(is=1 ; is <= nstate_dn_atm[iatm_atm_typ_cp[ipart]] ; is++){
   if( ((istate_dn + is ) >= istate_dn_st) && 
       ((istate_dn + is) <= istate_dn_end)){  
     ind = i + (istate_dn - istate_dn_st + is  )*ncoef;
     creal_dn[ind] =  helr*psi_r[is] - heli*psi_i[is];
     cimag_dn[ind] =  heli*psi_r[is] + helr*psi_i[is];
   }/*endif*/
   }/*endfor*/
  istate_dn +=  nstate_dn_atm[iatm_atm_typ_cp[ipart]];
 }/*endif*/

 }/*endfor*/
 }/*endfor:ncoef*/

/*=======================================================================*/
/*  g=0                                                                  */
    i = ncoef;
    istate_up = 0;
    istate_dn = 0;
  for(ipart=1; ipart <= nab_initio; ipart++){
   if( ((istate_up + 1 ) >= istate_up_st) && 
       ((istate_up + 1) <= istate_up_end)){  
     ind = i + ((1+istate_up) - istate_up_st )*ncoef;
     creal_up[ind] = ylmr[1]*gpsi00[iatm_atm_typ_cp[ipart]]/volrt;
     cimag_up[ind] = 0.0;
   }/*endif*/
    for(is = 2; is <= nstate_up_atm[iatm_atm_typ_cp[ipart]]; is++){
   if( ((istate_up + is ) >= istate_up_st) && 
       ((istate_up + is) <= istate_up_end)){  
     ind = i + (istate_up - istate_up_st + is  )*ncoef;
      creal_up[ind] = 0.0;
      cimag_up[ind] = 0.0;
   }/*endif*/
  }/*endfor*/

   istate_up += nstate_up_atm[iatm_atm_typ_cp[ipart]];

   if((cp_lsda == 1) &&(nstate_dn > 0)){
   if( ((istate_dn + 1 ) >= istate_dn_st) && 
       ((istate_dn + 1) <= istate_dn_end)){  
     ind = i + ((1+istate_dn) - istate_dn_st )*ncoef; 
     creal_dn[ind] = ylmr[1]*gpsi00[iatm_atm_typ_cp[ipart]]/volrt;
     cimag_dn[ind] = 0.0;
   }/*endif*/
    for(is = 2; is<=nstate_dn_atm[iatm_atm_typ_cp[ipart]]; is++){
   if( ((istate_dn + is ) >= istate_dn_st) && 
       ((istate_dn + is) <= istate_dn_end)){  
     ind = i + (istate_dn - istate_dn_st + is  )*ncoef;
       creal_dn[ind] = 0.0;
       cimag_dn[ind] = 0.0;
   }/*endif*/
    }/*endfor*/
      istate_dn += nstate_dn_atm[iatm_atm_typ_cp[ipart]];
   }/*endif*/

  }/*endfor:ipart*/


/*=========================================================================*/
/* check overlap before  gram-schmidt                                      */
#ifdef DEBUG_GW
  if(myid == master){ fprintf(stdout,"diag overlap before gram-schmidt \n");} 

   ncoef1 = ncoef -1;
  for(i= 1; i<= nstate_up_proc ; i++){
    o = 0.0;
   for(j=1; j<= ncoef1; j++){
      ind = (i-1)*ncoef + j;
      o += (creal_up[ind]*creal_up[ind] + cimag_up[ind]*cimag_up[ind]);
    }/*endfor*/
    j= ncoef;
    ind = (i-1)*ncoef + j;
    o *= 2.0;
    o  += (creal_up[ind]*creal_up[ind]);
    fprintf(stdout,"myid %d up  %d %.12lg \n",myid,(i-1 + istate_up_st),o);
  }/*endfor*/
   
  Dbx_Barrier(comm_states);}

/* check  LSDA */
   if((cp_lsda == 1) &&(nstate_dn > 0)){
    ncoef1 = ncoef -1;
   for(i=1; i<= nstate_dn_proc ; i++){
     o = 0.0;
    for(j=1; j<= ncoef1; j++){
      ind = (i-1)*ncoef + j;
      o += (creal_dn[ind]*creal_dn[ind] + cimag_dn[ind]*cimag_dn[ind]);
    }/*endfor*/
     j= ncoef;
     ind = (i-1)*ncoef + j;
     o = (creal_dn[ind]*creal_dn[ind] + 2.0*o);
     fprintf(stdout,"myid %d dn %d  %lg \n",myid,(i-1 + istate_dn_st),o);
   }/*endfor*/
  }/*endif*/
#endif

/*=========================================================================*/
/* create an orthonormal set with gram-schmidt                             */
/* GRAM-SCHMIDT ORTHOGONALIZE                                              */

  if(nproc>1){ Barrier(comm_states);}

  icoef_form_up = 0; /*(0) normal (1) transposed */
  icoef_form_dn = 0; /*(0) normal (1) transposed */

/* before can gram-schmidt, if nproc > 1 need to transpose data         */
   if( nproc > 1 ){
     cp_transpose_fwd(creal_up,cimag_up,&icoef_form_up,
                      creal_tmp,cimag_tmp,&(cp->cp_comm_state_pkg_up));
    if(cp_lsda==1 && nstate_dn != 0){
      cp_transpose_fwd(creal_dn,cimag_dn,&icoef_form_dn,
                      creal_tmp,cimag_tmp,&(cp->cp_comm_state_pkg_dn));
     }/*endif lsda*/
   }/*endif nproc*/



  if( nproc > 1 ){
    cp_gram_schmidt_par(creal_up,cimag_up,icoef_form_up,
                        occ_up,ioff_upt,
                        ovlap,ovlap_tmp,&(cp->cp_comm_state_pkg_up));
  }else{

   cp_gram_schmidt_scalar(creal_up,cimag_up,occ_up,
                          ovlap,ioff_upt,nstate_up,ncoef); 
  }/*endif nproc*/

 if( (cp_lsda==1) && (nstate_dn > 0)){
  if( nproc > 1 ){
    cp_gram_schmidt_par(creal_dn,cimag_dn,icoef_form_dn,
                        occ_dn,ioff_dnt,
                        ovlap,ovlap_tmp,&(cp->cp_comm_state_pkg_dn));
  }else{
    cp_gram_schmidt_scalar(creal_dn,cimag_dn,occ_dn,
                           ovlap,ioff_dnt,nstate_dn,ncoef);
  }/*endif nproc*/
 }/*endif lsda*/

 if(nproc>1){Barrier(comm_states);}

   printf("17\n");
/*=========================================================================*/
/* transpose data back if necessary                                        */
/* if nproc > 1 need to transpose data                                  */

  if(nproc > 1){
     cp_transpose_bck(creal_up,cimag_up,&icoef_form_up,
                      creal_tmp,cimag_tmp,&(cp->cp_comm_state_pkg_up));
   if(cp_lsda==1){
     cp_transpose_bck(creal_dn,cimag_dn,&icoef_form_dn,
                      creal_tmp,cimag_tmp,&(cp->cp_comm_state_pkg_dn));
     }/*endif lsda*/
   }/*endif nproc*/

  if(nproc>1){Barrier(comm_states);}

   printf("18\n");
/*=========================================================================*/
/*  check the overlaps after gram-schmidt                                  */

#ifdef DEBUG_GW
   if(myid == master){fprintf(stdout,"diag overlap after gram-schmidt \n");}

  ncoef1 = ncoef -1;

   for(i=1; i<= nstate_up_proc ; i++){
     o = 0.0;
    for(j=1; j<= ncoef1; j++){
     ind = j + (i-1)*ncoef;
     o += (creal_up[ind]*creal_up[ind] + cimag_up[ind]*cimag_up[ind]);
    }/*endfor*/
     j= ncoef;
     ind = j + (i-1)*ncoef;
     o *= 2.0;
     o += (creal_up[ind]*creal_up[ind]);
/*  fprintf(stdout,"myid %d up after gram_schmidt %d %lg \n",
                   myid,(i-1 + istate_up_st),o); */
    fprintf(stdout," up after gram_schmidt %d %lg \n",(i-1 + istate_up_st),o);
   }/*endfor*/

   if(nproc>1){Barrier(comm_states);}

   if(cp_lsda ==  1){
    for(i=1; i<= nstate_dn_proc; i++){
      o = 0.0;
    for(j=1; j<= ncoef1; j++){
      ind = j + (i-1)*ncoef;
      o += (creal_dn[ind]*creal_dn[ind] + cimag_dn[ind]*cimag_dn[ind]);
    }/*endfor*/
     j= ncoef;
     ind = j + (i-1)*ncoef;
     o *= 2.0;
     o += (creal_dn[ind]*creal_dn[ind]);
    fprintf(stdout,"dn %d %lg \n",(i-1+istate_dn_st),o);
   }/*endfor*/
  }/*endif*/
#endif

   printf("19\n");
/*=========================================================================*/
/*=========================================================================*/
/* check the initial kinetic energy                                        */
#ifdef DEBUG_GW
   ncoef1 = ncoef -1;
   eke_tot_proc = 0.0;
   for(j=1; j<= nstate_up_proc ; j++){
     eke = 0.0;
     for(i=1; i<= ncoef1; i++){
       aka = (double)(kastore[i]);
       akb = (double)(kbstore[i]);
       akc = (double)(kcstore[i]);
       xk = (aka*hmati[1] +akb*hmati[4] + akc*hmati[7])*tpi;
       yk = (aka*hmati[2] +akb*hmati[5] + akc*hmati[8])*tpi;
       zk = (aka*hmati[3] +akb*hmati[6] + akc*hmati[9])*tpi;
       g2 = (xk*xk+yk*yk+zk*zk);
       ind = i + (j-1)*ncoef; /*(i,j)*/
       eke +=  2.0*(creal_up[ind]*creal_up[ind] 
                  + cimag_up[ind]*cimag_up[ind])*g2;
     }/*endfor*/
     fprintf(stdout,"kinetic energy jth up state %d %d %lg \n",
                           myid,(j-1+istate_up_st),0.5*eke);
      eke_tot_proc +=  eke;
   }/*endfor*/

   if(nproc > 1){
     Reduce(&eke_tot_proc,&eke_tot,1,MPI_DOUBLE,MPI_SUM,master,comm_states);
   }else{
      eke_tot = eke_tot_proc;
   }/*endif*/
    
  if(myid == master ){
    fprintf(stdout,"total kinetic energy up states %lg \n",0.5*eke_tot);
  }/*endif myid*/
  if(nproc>1){ Barrier(comm_states);}

   if(cp_lsda == 1){
     eke_tot_proc = 0.0;
     ncoef1 = ncoef -1;
    for(j=1; j<= nstate_dn_proc; j++){
      eke = 0.0;
     for(i=1; i<= ncoef1; i++){
       aka = (double)(kastore[i]);
       akb = (double)(kbstore[i]);
       akc = (double)(kcstore[i]);
       xk = (aka*hmati[1] +akb*hmati[4] +akc*hmati[7])*tpi;
       yk = (aka*hmati[2] +akb*hmati[5] +akc*hmati[8])*tpi;
       zk = (aka*hmati[3] +akb*hmati[6] +akc*hmati[9])*tpi;
       g2 = (xk*xk+yk*yk+zk*zk);
       ind = i + (j-1)*ncoef;
       eke +=  2.0*(creal_dn[ind]*creal_dn[ind] 
                  + cimag_dn[ind]*cimag_dn[ind])*g2;
     }/*endfor*/
     fprintf(stdout,"kinetic energy jth down state %d %lg \n",myid+j,0.50*eke);
       eke_tot_proc += eke;
    }/*endfor*/

   if(nproc > 1){
     Reduce(&eke_tot_proc,&eke_tot,1,MPI_DOUBLE,MPI_SUM,master,comm_states);
   }else{
      eke_tot = eke_tot_proc;
   }/*endif*/
  
   if(myid == master){
    fprintf(stdout,"total kinetic energy dn states %lg \n",0.50*eke_tot);
   }/*endif myid*/
  }/*endif lsda*/
#endif

/*===========================================================================*/
/* free locally assigned memory */

   cfree(&(iatm_atm_typ_cp[1]));
  if(myid == 0){ cfree(&(fname_ps[1])); }
   cfree(&(nstate_up_atm[1]));
   cfree(&(nstate_dn_atm[1]));

   cfree(&(x[1]));
   cfree(&(y[1]));
   cfree(&(z[1]));

   cfree(&(n_ang[1]));
   cfree(&(ylmr[1]));
   cfree(&(ylmi[1]));

   cfree(&(dylmr_x[1]));
   cfree(&(dylmi_x[1]));
   cfree(&(dylmr_y[1]));
   cfree(&(dylmi_y[1]));
   cfree(&(dylmr_z[1]));
   cfree(&(dylmi_z[1]));

   cfree(&psi_r[1]);
   cfree(&psi_i[1]);

 for(i=1; i<= natm_typ_cp; i++){

  for(j=1; j<=3; j++){
    cfree(&(gpsi0[i][j][1]));
    cfree(&(gpsi1[i][j][1]));
    cfree(&(gpsi2[i][j][1]));
    cfree(&(gpsi3[i][j][1]));
   }/*endfor*/

   cfree(&(gpsi0[i][1]));
   cfree(&(gpsi1[i][1]));
   cfree(&(gpsi2[i][1]));
   cfree(&(gpsi3[i][1]));
 }/*endfor*/

   cfree(&(gpsi0[1]));
   cfree(&(gpsi1[1]));
   cfree(&(gpsi2[1]));
   cfree(&(gpsi3[1]));

  cfree_mat(gpsi_now,1,natm_typ_cp,1,3);
  cfree(&(gpsi00[1]));

   printf("20\n");
/*===========================================================================*/

  if(myid == master ){
     printf("\n-----------------------------------------\n");
     printf("Completed Constructing The Wave Function \n");
     printf("-----------------------------------------\n");
  }/*endif*/

/*==========================================================================*/
    }/*  end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  splin_btrans(int *nsplin,double **gpsi0,double **gpsi1,
                    double **gpsi2,double **gpsi3,
                    double gmin,double gmax,double *dg,
                    double *gpsi00,int *pn_ang,char *fname_ps,int atm,
                    int myid,int nproc,MPI_Comm comm_states)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

#include "../typ_defs/typ_mask.h"

   double *g,*r,*rphi; /* length nsplin */

   int iii;
   int i,ir,iang_now;
   int iang,nr;
   int n_ang;
   double rmax,xx,dr;
   int nsplin_loc = *nsplin;
   int master = 0;
   int n_ang1;
   FILE *fp_name_ps;

/*========================================================================*/
/* Set up the g's */

   g   = (double *) cmalloc(nsplin_loc*sizeof(double)) -1;
   *dg = (gmax-gmin)/(double)(nsplin_loc);

   for(i=1; i<= nsplin_loc; i++){
     g[i] = (*dg)*(double)(i-1) + gmin;
   }/*endfor*/

/*========================================================================*/
/* Set up the r's */

  if(myid == 0){
   fp_name_ps = cfopen(fname_ps,"r");
   fscanf(fp_name_ps,"%d %lg %d ",&nr,&rmax,&n_ang);
   readtoendofline(fp_name_ps);
   readtoendofline(fp_name_ps);
  }

  if( nproc > 1){
     Bcast(&nr,1,MPI_INT,master,comm_states);
     Bcast(&n_ang,1,MPI_INT,master,comm_states);
     Bcast(&rmax,1,MPI_DOUBLE,master,comm_states);
  }  

     *pn_ang = n_ang;

  r    = (double *) cmalloc(nr*sizeof(double)) -1;
  rphi = (double *) cmalloc(nr*sizeof(double)) -1;

  dr = rmax/(double)(nr);

  for(ir=1; ir<= nr; ir++){
    r[ir] = (double)(ir-1)*dr;
  }/*endfor*/
  
/*========================================================================*/
/* Set up the gpsi's */
   
  n_ang1 = (*pn_ang) + 1;

  for(iang=1; iang <= n_ang1; iang++){
  
   if(myid == 0){
    for(ir=1; ir<= nr; ir++){
      fscanf(fp_name_ps,"%lg %lg ",&xx,&(rphi[ir]));
    }/*endfor*/
   }/*endif myid*/

   if(nproc>1){
    Bcast(&(rphi[1]),nr,MPI_DOUBLE,master,comm_states);
   }

    iang_now = iang-1;

    bess_trans(rphi,nr,dr,r,gpsi0[iang],*nsplin,g,iang_now,gpsi00);
    fit_spline(gpsi0[iang],gpsi1[iang],gpsi2[iang],gpsi3[iang],g,*nsplin);

  }/*endfor*/

  if( myid == master){
    fclose(fp_name_ps);
  }/*endif myid*/

/*========================================================================*/
/* free locally assigned memory                                           */

  cfree(&(g[1]));
  cfree(&(r[1]));
  cfree(&(rphi[1]));

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bess_trans(double *v_rphi,int nr,double dr,double *r,
             double *fv_rphi,int nsplin_g,double *g,
             int iang,double *gzero)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/


/*  local variables                                                         */
   double fpidr,rj0,rj1,rj2,rj3,arg,fpi,pi,tpi;
   int ig,ir;
   int iii;

/* -------------------------------------------------------------------------*/
/*  slow spherical bessel fourier transform                                 */

   pi  = M_PI;
   tpi = 2.0*pi;
   fpi = 4.0*pi;
   fpidr = fpi*dr;


  for(ig=1; ig <= nsplin_g; ig++){
    fv_rphi[ig] = 0.0;
  }/*endfor*/

  if(iang == 0){
/* l=0 or local ,g=0 ----------------------------------------------------- */
   *gzero = 0.0;
  for(ir=1; ir <= nr; ir++){
    *gzero += fpidr*r[ir]*v_rphi[ir];
  }/*endfor*/

/* l=0 or local, g ne 0 ----------------------------------------------------*/
   for(ig=1; ig <= nsplin_g; ig++){
     for(ir=2; ir <= nr; ir++){
        arg = r[ir]*g[ig];
        rj0 = sin(arg)/arg*r[ir];
        fv_rphi[ig] += fpidr*rj0*v_rphi[ir];
      }/*endfor*/
   }/*endfor*/
 }/*endif*/

   if(iang == 1){
/* l=1, g ne 0 --------------------------------------------------------------*/
   for(ig=1; ig <= nsplin_g; ig++){
     for(ir=2 ; ir <= nr; ir++){
       arg = r[ir]*g[ig];
       rj1 = (sin(arg)/arg - cos(arg))/arg*r[ir];
       fv_rphi[ig] +=  fpidr*rj1*v_rphi[ir];
     }/*endfor*/
   }/*endfor*/
  }/*endif*/

   if(iang == 2){
/*  l=2, g ne 0 ------------------------------------------------------------*/
   for(ig=1; ig <= nsplin_g; ig++){
     for(ir=2; ir <= nr; ir++){  
        arg = r[ir]*g[ig];
        rj2 = ((3.0/(arg*arg)-1.0)*sin(arg)-3.0*cos(arg)/arg)/arg*r[ir];
        fv_rphi[ig] += fpidr*rj2*v_rphi[ir];
     }/*endfor*/
   }/*endfor*/
  }/*endif*/

  if(iang == 3){
/*  l=3, g ne 0 ---------------------------------------------------------- */
   for(ig=1; ig <= nsplin_g; ig++){
     for(ir=2; ir <= nr; ir++){
       arg = r[ir]*g[ig];
       rj3 = ((15.0/(arg*arg) - 6.0)*sin(arg)/arg + 
             (1.0 - 15.0/(arg*arg))*cos(arg))/arg*r[ir];
      fv_rphi[ig] =+  fpidr*rj3*v_rphi[ir];
     }/*endfor*/
   }/*endfor*/
  }/*endif*/


/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void get_gpsi(double g,int nsplin,
               double **gpsi0,double **gpsi1,double **gpsi2,double **gpsi3,
               double *gpsi_now,double gmin,double dg,int n_ang)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/


   double partem1,partem2,partem3,partem4,h,h0;
   int iii,iang,i;

   iii = (g-gmin)/dg + 1;
   iii = MIN(iii,nsplin);
   iii = MAX(iii,1);

   h0  = (double)(iii-1)*dg+gmin;
   h = g-h0;

   for(iang =1; iang <= (n_ang+1); iang++){
     partem1 = gpsi0[iang][iii];
     partem2 = gpsi1[iang][iii];
     partem3 = gpsi2[iang][iii];
     partem4 = gpsi3[iang][iii];

     gpsi_now[iang] = ((partem4*h+partem3)*h+partem2)*h+partem1;
   }/*endfor*/

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  fit_spline(double *c0i,double *c1i,double *c2i,
                 double *c3i,double *xi,int nsplin)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
/*  positions                                                              */
/* spline coefficients, c(i,1)=function values                             */

   int iii;
/* temporary vectors                                                       */
   double *d,*diag; /*length nsplin*/
/* temporary scalars                                                       */
   double c0,c1,c2,c3,g,divdf1,divdf3,dx;
/* temporary integers                                                      */
   int m,mm1,mp1,i,ip1,n;

/*malloc local memory */
    d    = (double *) cmalloc(nsplin*sizeof(double )) -1;
    diag = (double *) cmalloc(nsplin*sizeof(double )) -1;


/* fit the spline                                                          */

/*  1st approximate initial and final derivatives                          */
   n = nsplin-1;
   c1i[1]   = (c0i[2]-c0i[1])/(xi[2]-xi[1]);
   c1i[n+1] = (c0i[n+1]-c0i[n])/(xi[n+1]-xi[n]);
   c0 = 0.0;
   c1 = 1.0;
   c2 = 2.0;
   c3 = 3.0;
   diag[1] = c1;
   d[1]    = c0;
  for(m=2; m<= n+1; m++){
    mm1 = m-1;
    d[m] = xi[m]-xi[mm1];
    diag[m] = (c0i[m]-c0i[mm1])/d[m];
  }/*endfor*/
  for(m=2; m<=n; m++){
     mp1 = m+1;
     c1i[m] = c3*(d[m]*diag[mp1]+d[mp1]*diag[m]);
     diag[m] = c2*(d[m]+d[mp1]);
  }/*endfor*/
  for(m=2; m<=n; m++){
    mp1 = m+1;
    mm1 = m-1;
    g = -d[mp1]/diag[mm1];
    diag[m] = diag[m]+g*d[mm1];
    c1i[m] = c1i[m]+g*c1i[mm1];
  }/*endfor*/
  for(m=n; m>=2; m--){
    mp1 = m+1;
    c1i[m] = (c1i[m]-d[m]*c1i[mp1])/diag[m];
  }/*endfor*/

/* calculate all other coefficients                                       */
  for(i=1; i<=n; i++){
    ip1 = i+1;
    dx = xi[ip1]-xi[i];
    divdf1 = (c0i[ip1]-c0i[i])/dx;
    divdf3 = c1i[i]+c1i[ip1]-c2*divdf1;
    c2i[i] = (divdf1-c1i[i]-divdf3)/dx;
    c3i[i] = divdf3/(dx*dx);
  }/*endfor*/


/* free locally assigned memory */
    cfree(&(d[1]));
    cfree(&(diag[1]));

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/












