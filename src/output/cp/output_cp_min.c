/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: output_md                                    */
/*                                                                          */
/* This subprogram provides output for a MD on a                            */ 
/* classical potential energy surface (PES)                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void output_cp_min(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                   CP *cp,double Delta_E,int iexit)

/*==========================================================================*/
{/*begin routine*/

/*=======================================================================*/
/*            Local variable declarations                                */
  int iii;
  double etot;
  double vtot;
  double a,b,c,tab,tac,tbc; 
  double deth;
  double c_div;

/*            Local Pointer */
  int myid              = class->communicate.myid;
  int nproc             = class->communicate.np;
  int pi_beads          = class->clatoms_info.pi_beads;
  MPI_Comm world        = class->communicate.world;
  int np_states         = cp->communicate.np_states;
  int exit_flag         = general_data->timeinfo.exit_flag;


/*=====================================================================*/
/* 0) Adjust form of CP vectors */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
    control_coef_transpose_bck(cp,1);
 }/*endif*/

/*=====================================================================*/
/*  I)Open File option :                                               */

  if(myid==0){
    if(((general_data->timeinfo.itime)==0)&&
       ((general_data->filenames.ifile_open))==1){
     initial_fopen_cp(class,general_data,bonded,cp);
    }/*endif*/
  }/*endif*/

  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  II) Write Initial Energies to screen                                 */

  if(myid==0){
   if(general_data->timeinfo.itime==0 && general_data->simopts.debug_cp 
      && general_data->simopts.cp_wave_min == 0){
    initial_output_cp_min(class,general_data,bonded,cp);
   }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/* II) Calculate some dinky quantities                                   */

  if(myid==0){
   if((general_data->timeinfo.itime)!=0){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_screen) == 0 || 
       (general_data->timeinfo.itime % 
        general_data->filenames.iwrite_inst)==0){

      get_cell(general_data->cell.hmat,&a,&b,&c,&tab,&tbc,&tac);
      dink_quant_calc_cp_min(class, general_data,cp,&etot,&vtot,&deth,&c_div);

    } /*endif*/
   }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  III) Write to the output to screen                                   */

  if(myid==0){
   if(iexit==0){
    if((general_data->timeinfo.itime) != 0){
     if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_screen) == 0){ 

      screen_write_cp_min(class,general_data,bonded,cp,
			  etot,vtot,deth,Delta_E,
                          a,b,c,tab,tbc,tac,c_div);
     }/*endif*/
    }/*endif*/
   }/*endif */
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* IV) Write to the output to dump file                                 */ 

   if((general_data->timeinfo.itime!=0)){
    if((general_data->timeinfo.itime%
        general_data->filenames.iwrite_dump)==0 || exit_flag == 1){
      if(pi_beads==1){
        write_dump_file_cp(class,bonded,general_data,cp);
       }else{
       if(myid == 0){
        write_dump_file_cp_pimd_class(class,bonded,general_data,cp);
       }
        if(nproc>1){Barrier(world);}
      }/* endif */
    }/*endif*/
   }/*endif*/

/*====================================================================*/
/* VI) Write to the config files                          */

  if((general_data->timeinfo.itime)!=0){
    write_config_files_cp(class,bonded,general_data,cp);
  }/*endif*/

/*=====================================================================*/
/* 0) Adjust form of CP vectors */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
    control_coef_transpose_fwd(cp,1);
 }/*endif*/

/*==========================================================================*/
   }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_output_cp_min(CLASS *class,GENERAL_DATA *general_data,
                           BONDED *bonded,CP *cp)

/*==========================================================================*/

{/* begin routine */

  int ncons,np_tot,npairs,nsh_tot;
  double c_div;

/*==========================================================================*/
/*  Calculate c_div                                                         */

  if(cp->cpopts.cp_lsda == 0) {
    if(cp->cpopts.cp_norb == 0) {
        ncons = (cp->cpcoeffs_info.nstate_up)*
               ((cp->cpcoeffs_info.nstate_up)+1)/2;
    }else{ ncons = cp->cpcoeffs_info.nstate_up; }
    c_div = (double) (2*(cp->cpcoeffs_info.ncoef)-ncons);
  }/* endif */

  if(cp->cpopts.cp_lsda == 1) {
    if(cp->cpopts.cp_norb == 0) {
     ncons = (cp->cpcoeffs_info.nstate_up)*((cp->cpcoeffs_info.nstate_up)+1)/2+ 
               (cp->cpcoeffs_info.nstate_dn)*
              ((cp->cpcoeffs_info.nstate_dn)+1)/2;
    } else { ncons = cp->cpcoeffs_info.nstate_up+cp->cpcoeffs_info.nstate_dn; }
    c_div = (double) (4*(cp->cpcoeffs_info.ncoef)-ncons);
  } /* endif */

/*==========================================================================*/
/* Output stuff */

  if(general_data->simopts.cp_min==1){
    printf("Initial intra      energy  %.10g\n",general_data->stat_avg.vintrat);
    printf("Initial inter      energy  %.10g\n",general_data->stat_avg.vintert);
    if(class->energy_ctrl.isep_vvdw == 1) {
      printf("Initial VDW        energy  %.10g\n",general_data->stat_avg.vvdw);
      printf("Initial Coul       energy  %.10g\n",general_data->stat_avg.vcoul);
    }/*endif*/
    printf("Initial Bond       energy  %.10g\n",general_data->stat_avg.vbondt);
    printf("Initial Bend       energy  %.10g\n",general_data->stat_avg.vbendt);
    printf("Initial Bendbnd    energy  %.10g\n",general_data->stat_avg.vbend_bndt);
    printf("Initial Torsion    energy  %.10g\n",general_data->stat_avg.vtorst);
    printf("Initial onefour    energy  %.10g\n",general_data->stat_avg.vonfot);
  }/*endif*/

    printf("Initial e-Hart+XC  energy  %.10g\n",general_data->stat_avg.cp_ehart
				           + general_data->stat_avg.cp_exc);
    printf("Initial e-Hart+muXC-XC  energy  %.10g\n",general_data->stat_avg.cp_ehart
				           + general_data->stat_avg.cp_muxc
				           - general_data->stat_avg.cp_exc);
    printf("Initial e-External energy  %.10g\n",general_data->stat_avg.cp_eext);
    printf("Initial e-Nonlocal energy  %.10g\n",general_data->stat_avg.cp_enl);
    printf("Initial e-Kinetic  energy  %.10g\n",general_data->stat_avg.cp_eke);
    printf("Initial volume             %.10g\n",general_data->stat_avg.vol);

/*==========================================================================*/
  /* Neighbor list quantities                                              */

  if(general_data->simopts.cp_min==1){
    if((class->nbr_list.iver)==1){
      npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
	+ class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
      np_tot = (class->clatoms_info.natm_tot)*
               (class->clatoms_info.natm_tot-1)/2 
	- (bonded->excl.nlst);
      printf("Number of pairs = %d out of %d\n",npairs,np_tot);
      if(general_data->timeinfo.int_res_ter==1){
	npairs = class->nbr_list.verlist.nter_res[(class->clatoms_info.natm_tot)] 
	  + class->nbr_list.verlist.jver_off_res[(class->clatoms_info.natm_tot)];
	printf("Number of pairs = %d out of %d\n",npairs,np_tot);
      }/*endif*/
    }/*endif*/

    if((class->nbr_list.ilnk)==1){
      nsh_tot = ((class->nbr_list.lnklist.ncell_a)*(class->nbr_list.lnklist.ncell_b)*
		 (class->nbr_list.lnklist.ncell_c)-1)/2 +
		   (((class->nbr_list.lnklist.ncell_a)*(class->nbr_list.lnklist.ncell_b)*
		     (class->nbr_list.lnklist.ncell_c)-1) % 2) + 1;
      printf("ncell_a = %d, ncell_b = %d, ncell_c = %d natm_cell_max %d \n",
	     class->nbr_list.lnklist.ncell_a,class->nbr_list.lnklist.ncell_b,
	     class->nbr_list.lnklist.ncell_c,class->nbr_list.lnklist.natm_cell_max);
      printf("# of shifts = %d out of %d\n",
	     class->nbr_list.lnklist.nshft_lnk,nsh_tot);
      if(general_data->timeinfo.int_res_ter==1){
	printf("# of respa shifts = %d out of %d\n",
	       class->nbr_list.lnklist.nshft_lnk_res,nsh_tot);
      }/*endif*/
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dink_quant_calc_cp_min(CLASS *class, GENERAL_DATA *general_data,CP *cp,
                            double *etot,double *vtot,double *deth,
                            double *c_div)

/*==========================================================================*/

{/* begin routine */

  int i,iii,ip,ncons,ncoef_tot;
  int write_now = general_data->stat_avg.write_cp_atm_flag;
  int pi_beads = class->clatoms_info.pi_beads;

/*==========================================================================*/
 /*  Energy */

     (*etot) = general_data->stat_avg.cp_ehart + 
	       general_data->stat_avg.cp_eext + general_data->stat_avg.cp_exc + 
               general_data->stat_avg.cp_eke + general_data->stat_avg.cp_enl;
      if(write_now == 1) {
        (*etot) += (general_data->stat_avg.vintert) + (general_data->stat_avg.vintrat);
      }/* endif  cp_min on*/

 /*==========================================================================*/
 /*  Volume                                                                  */

    if(write_now == 1) {
      (*vtot) = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;
      (*deth) = getdeth(general_data->cell.hmat);
    }/*endif*/

/*==========================================================================*/
/*  Calculate c_div and c_nhc_div                                           */

  if(cp->cpopts.cp_lsda == 0) {
    if(cp->cpopts.cp_norb == 0) {
        ncons = (cp->cpcoeffs_info.nstate_up)*
               ((cp->cpcoeffs_info.nstate_up)+1)/2;
    }else{ ncons = cp->cpcoeffs_info.nstate_up; }
    (*c_div) = (double) (2*(cp->cpcoeffs_info.ncoef)-ncons);
  }/* endif */
  if(cp->cpopts.cp_lsda == 1) {
    if(cp->cpopts.cp_norb == 0) {
     ncons = (cp->cpcoeffs_info.nstate_up)*
              ((cp->cpcoeffs_info.nstate_up)+1)/2 + 
               (cp->cpcoeffs_info.nstate_dn)*
              ((cp->cpcoeffs_info.nstate_dn)+1)/2;
    } else { ncons = cp->cpcoeffs_info.nstate_up+cp->cpcoeffs_info.nstate_dn; }
    (*c_div) = (double) (4*(cp->cpcoeffs_info.ncoef)-ncons);
  } /* endif */

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void screen_write_cp_min(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,CP *cp,
                         double etot,double vtot,double deth,double Delta_E,
                         double a,double b,double c,
                         double tab,double tbc,double tac,double c_div)

/*==========================================================================*/
{/* begin routine */

 int npairs,ip;
 double atime,atm_div;
 double updates_t,updates_true,updates_now; 
 int pi_beads = class->clatoms_info.pi_beads;
 int write_now = general_data->stat_avg.write_cp_atm_flag;
 double fc_mag_up = general_data->stat_avg.fc_mag_up;
 double fc_max_up = general_data->stat_avg.fc_max_up;
 double fc_mag_dn = general_data->stat_avg.fc_mag_dn;
 double fc_max_dn = general_data->stat_avg.fc_max_dn;
 double fatm_mag  = general_data->stat_avg.fatm_mag;
 double fatm_max  = general_data->stat_avg.fatm_max;

/*==========================================================================*/
/* Write to screen                                                          */

     atime = (double)(general_data->timeinfo.itime);

/*==========================================================================*/
/*     A) Standard                                                   */
      printf("\n");
      printf("********************************************************\n");
      printf("QUANTITY          =  MAGNITUDE        MAXIMUM           \n");
      printf("--------------------------------------------------------\n");
     if(write_now == 1 || general_data->simopts.minimize == 1){
      if(general_data->minopts.min_std ==1)  
        printf("Atm minimization      = STD \n");
      if(general_data->minopts.min_cg  ==1)  
        printf("Atm minimization      = CG  \n");
      if(general_data->minopts.min_diis==1)  
        printf("Atm minimization      = DIIS \n");
     }/*endif*/
     if(general_data->simopts.cp_wave_min == 1 || 
        general_data->simopts.cp_wave_min_pimd ==1){
      if(general_data->minopts.cp_min_std ==1)  
        printf("Coef minimization      = STD \n");
      if(general_data->minopts.cp_min_cg  ==1)  
        printf("Coef minimization      = CG  \n");
      if(general_data->minopts.cp_min_diis==1)  
        printf("Coef minimization      = DIIS \n");
     }/*endif*/
     printf("Time step         = %g\n",(double)(general_data->timeinfo.itime));
     printf("----------------- \n");
     printf("Coef force up    = %g          %g\n",fc_mag_up,
                                                         fc_max_up);
     if(cp->cpopts.cp_lsda == 1) {
       printf("Coef force dn    = %g          %g\n",fc_mag_dn,
                                                          fc_max_dn);
     }/*endif*/
     if(write_now==1){
      printf("Atm force        = %g           %g\n",fatm_mag,fatm_max);
     }/*endif*/
     if(write_now==1){
      printf("----------------- \n");
      printf("Total Energy      = %.10g \n",(vtot) +
                                (general_data->stat_avg.cp_ehart
			       + general_data->stat_avg.cp_exc
			       + general_data->stat_avg.cp_eext
			       + general_data->stat_avg.cp_eke
			       + general_data->stat_avg.cp_enl));
      printf("----------------- \n");
      printf("Atm Energy        = %.10g \n",(vtot));
      printf("Intermol Atm PE   = %.10g \n",(general_data->stat_avg.vintert));
      printf("Intramol Atm PE   = %.10g \n",(general_data->stat_avg.vintrat));
      printf("----------------- \n");
      atm_div = (double)(class->clatoms_info.nfree);
      printf("Atm Deg. Free     = %g\n",atm_div); 
    }/* endif */
      printf("----------------- \n");
      printf("Delta E elec     = %g\n",Delta_E);
      printf("e-Energy          = %.10g \n",
                                (general_data->stat_avg.cp_ehart
			       + general_data->stat_avg.cp_exc
			       + general_data->stat_avg.cp_eext
			       + general_data->stat_avg.cp_eke
			       + general_data->stat_avg.cp_enl));
      printf("e-Hartree + XC    = %.10g \n",(general_data->stat_avg.cp_ehart
				 	  + general_data->stat_avg.cp_exc));
      printf("e-External PE     = %.10g \n",(general_data->stat_avg.cp_eext));
      printf("e-Nonlocal PE     = %.10g \n",(general_data->stat_avg.cp_enl));
      printf("e-Kinetic         = %.10g \n",(general_data->stat_avg.cp_eke));
/*==========================================================================*/
/*     D)Constraint                                                  */

      if(bonded->constrnt.iconstrnt == 1) {
	if(write_now==1){
	  printf("Shake iter        = %g %g\n",
		 (double)(general_data->stat_avg.iter_shake),
		 (general_data->stat_avg.aiter_shake/atime));
	  printf("Rattle iter       = %g %g\n",
		 (double)(general_data->stat_avg.iter_ratl),
		 (general_data->stat_avg.aiter_ratl/atime));
	  printf("----------------- \n");
	}/*endif*/
      }/*endif*/

/*==========================================================================*/
/*     D)Misc                                                         */

      if(class->nbr_list.iver == 1 && write_now == 1){
	  updates_t = general_data->stat_avg.updates;
	  if(updates_t == 0) updates_t = 1;
	  updates_now = (double ) (general_data->timeinfo.itime-
	 			 general_data->stat_avg.itime_update);
          updates_true = updates_t;
          npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
	       + class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
          printf("Inst steps/update = %g\n",updates_now);
          printf("Avg. steps/update = %g\n",(atime/updates_t));
          printf("Total list updates= %g\n",updates_true);
          printf("Number of pairs   = %d \n",npairs);
	  printf("----------------- \n");
      }/*endif*/
      printf("Cpu time          = %g %g\n",(general_data->stat_avg.cpu_now),
	     (general_data->stat_avg.acpu/atime));
      printf("--------------------------------------------------------\n");
      printf("\n"); 
      printf("********************************************************\n");
      fflush(stdout);

/*==========================================================================*/
  }/* end routine */
/*==========================================================================*/



