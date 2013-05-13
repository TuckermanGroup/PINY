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
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"


void write_dump_header_cp_pimd(FILE *,CP *,GENERAL_DATA *,int ,int );
void write_dump_occ_cp_pimd(FILE *,CP *,int ,int );
void write_dump_coef_cp_pimd(FILE *,CP *,CLASS *,int );
void write_dump_vcoef_cp_pimd(FILE *,CP *,CLASS *,GENERAL_DATA *,int );
void write_dump_extended_cp_pimd(FILE *,CP *,CLASS *,GENERAL_DATA *,int );

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void output_cp_pimd(CLASS *class,GENERAL_DATA *general_data,
                    BONDED *bonded,CP *cp)

/*==========================================================================*/
{/*begin routine*/

/*=======================================================================*/
/*            Local variable declarations                                */

  int iii;
  double etot,econv_now; 
  double vtot,avtot,pnow; 
  double apnow;
  double nhc_div,nhc_div_bead;
  double a,b,c,tab,tac,tbc; 
  double deth;
  double c_nhc_div;
  double c_div = cp->cpcoeffs_info.cp_nfree;
  int exit_flag  = general_data->timeinfo.exit_flag;
  int myid       = class->communicate.myid;
  MPI_Comm world = class->communicate.world;
  int np_states  = class->communicate.np_states;
  int nproc      = class->communicate.np;

/*=====================================================================*/
/* 0) Adjust form of CP vectors */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
   if((general_data->timeinfo.itime % 
       general_data->filenames.iwrite_dump)==0 || exit_flag == 1){
    control_coef_transpose_bck(cp,3);
   }/*endif*/

   if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_confc) == 0) {
    control_coef_transpose_bck(cp,1);
   }/*endif*/
 }/*endif*/

/*=====================================================================*/
/*  I)Open File option :                                               */

  if(myid==0){
    if(((general_data->timeinfo.itime)==0)&&
       ((general_data->filenames.ifile_open))==1){
       initial_fopen_cp_pimd(class,general_data,bonded,cp);
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  II) Write Initial Energies to screen                                 */

  if(myid==0){
    if(general_data->timeinfo.itime==0 && 
       general_data->simopts.debug_cp_pimd == 0){
      initial_output_cp_pimd(class,general_data,bonded,cp);
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/* II) Calculate some dinky quantities                                   */

  if(myid==0){
    if((general_data->timeinfo.itime)!=0){
      if((general_data->timeinfo.itime
         %general_data->filenames.iwrite_screen) == 0 || 
        (general_data->timeinfo.itime
        %general_data->filenames.iwrite_inst)==0 || 
        (exit_flag == 1)){
        get_cell(general_data->cell.hmat,&a,&b,&c,&tab,&tbc,&tac);
        dink_quant_calc_cp_pimd(class, general_data,cp,&etot,&econv_now,&vtot,
                              &avtot,&deth,&pnow, &apnow,&nhc_div,
                              &nhc_div_bead,&c_nhc_div);
      } /*endif*/
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  III) Write  the output to screen                                     */

  if(myid==0){
    if((general_data->timeinfo.itime) != 0){
      if((general_data->timeinfo.itime 
        % general_data->filenames.iwrite_screen) == 0 ||
          (exit_flag == 1)){ 
        screen_write_cp_pimd(class,general_data,bonded,cp,
                          etot,econv_now,vtot,avtot,deth,pnow, apnow,
                          nhc_div,nhc_div_bead,a,b,c,tab,tbc,tac,c_div,
                          c_nhc_div);
      }/*endif*/
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* IV) Write to the output to classical dump file                       */ 

  if(myid==0){
    if((general_data->timeinfo.itime!=0)){
      if((general_data->timeinfo.itime 
        % general_data->filenames.iwrite_dump)==0 ||
          (exit_flag == 1)){
        write_dump_file_cp_pimd_class(class,bonded,general_data,cp);
      }/*endif*/
     }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* V) Write to the output to free energy                                */ 

  if(myid==0){
    if((general_data->timeinfo.itime)!=0&&
      ((general_data->simopts.cp+general_data->simopts.cp_pimd)==1)){
      if((general_data->timeinfo.itime 
        % general_data->filenames.iwrite_dump) == 0){
        write_free_energy_file_pimd(class,bonded,general_data);
      }/*endif*/
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*====================================================================*/
/* VI) Write to the config classical files                            */

  if(myid==0){
    if((general_data->timeinfo.itime)!=0){
      write_config_files_cp_pimd_class(class,bonded,general_data,cp);
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*====================================================================*/
/* VI) Write to the CP dump and config files                          */


 if((general_data->timeinfo.itime)!=0){
  if((general_data->timeinfo.itime %
       general_data->filenames.iwrite_dump) == 0){
    write_dump_file_cp_pimd(class,bonded,general_data,cp);
  }/*endif*/
 }/*endif*/

/*======================================================================*/
/* V) Write to the kseigs                                               */ 

  if((general_data->timeinfo.itime)!=0&&(general_data->simopts.cp==1)){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_kseigs) == 0){
        write_kseigs_file_cp(general_data,cp);
    }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* VI) Write to the ELF file                                             */ 

  if(cp->cpcoeffs_info.cp_elf_calc_frq > 0){
   if((general_data->timeinfo.itime)!=0 && 
    (general_data->simopts.cp==1 || general_data->simopts.cp_wave==1 )){
    if((general_data->timeinfo.itime % 
        cp->cpcoeffs_info.cp_elf_calc_frq) == 0){
         write_elf_file_cp(general_data,cp);
     }/*endif*/
    }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* VII) Write to the inst avgs to inst file                            */

  if(myid==0){
    if((general_data->timeinfo.itime)!=0){
      if((general_data->timeinfo.itime 
        % general_data->filenames.iwrite_inst) == 0 ){
        write_inst_file_cp_pimd(class,general_data,cp,etot,a,b,c,tac,tab,tbc);
      }/*endif*/
    }/*endif*/
  }/* endif : myid=0 */
  if(nproc>1){Barrier(world);}

/*=====================================================================*/
/* 0) Adjust form of CP vectors */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
   if((general_data->timeinfo.itime % 
       general_data->filenames.iwrite_dump)==0 || exit_flag == 1){
    control_coef_transpose_fwd(cp,2);
   }/*endif*/

   if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_confc) == 0) {
    control_coef_transpose_fwd(cp,1);
   }/*endif*/
 }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_output_cp_pimd(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                              CP *cp)

/*==========================================================================*/

{/* begin routine */

  int pi_beads = class->clatoms_info.pi_beads,iii;
  int ncons,np_tot,npairs,nsh_tot;
  double nhc_div,c_nhc_div;
  double c_div = cp->cpcoeffs_info.cp_nfree;

  double eu_conv=1.0;

  if(general_data->filenames.iwrite_units==0){eu_conv=1.0;}
  if(general_data->filenames.iwrite_units==1){eu_conv=KCAL;}
  if(general_data->filenames.iwrite_units==2){eu_conv=BOLTZ;}

/*==========================================================================*/
/* Output stuff */

  if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){
    printf("Initial intra      energy  %g\n",general_data->stat_avg.vintrat*eu_conv);
    printf("Initial inter      energy  %g\n",general_data->stat_avg.vintert*eu_conv);
    printf("Initial kinetic    energy  %g\n",general_data->stat_avg.kinet*eu_conv);
    printf("Initial prim kinetic      energy  %.10g\n",
                           general_data->stat_avg.pi_ke_prim*eu_conv);
    printf("Initial vir kinetic      energy  %.10g\n",
                           general_data->stat_avg.pi_ke_vir*eu_conv);
    printf("Initial harm kinetic energy  %.10g\n",
                           general_data->stat_avg.kin_harm*eu_conv);
    if(class->energy_ctrl.isep_vvdw == 1) {
      printf("Initial VDW        energy  %g\n",general_data->stat_avg.vvdw*eu_conv);
      printf("Initial Coul       energy  %g\n",general_data->stat_avg.vcoul*eu_conv);
    }/*endif*/
  printf("Initial Bond         energy  %.10g\n",
                           general_data->stat_avg.vbondt*eu_conv);
  printf("Initial Bendbnd_bond energy  %.10g\n",
                           general_data->stat_avg.vbend_bnd_bond*eu_conv);
  printf("Initial Watts_bond energy  %.10g\n",
                           general_data->stat_avg.vbondt_watts*eu_conv);
  printf("Initial Total Bond   energy  %.10g\n",
                           (general_data->stat_avg.vbend_bnd_bond+
                            general_data->stat_avg.vbondt+
                            general_data->stat_avg.vbondt_watts)*eu_conv);
  printf("Initial Bend         energy  %.10g\n",
                           general_data->stat_avg.vbendt*eu_conv);
  printf("Initial Bendbnd_bend energy  %.10g\n",
                           general_data->stat_avg.vbend_bnd_bend*eu_conv);
  printf("Initial Watts_bend energy  %.10g\n",
                           general_data->stat_avg.vbendt_watts*eu_conv);
  printf("Initial Total Bend   energy  %.10g\n",
                           (general_data->stat_avg.vbend_bnd_bend+
                            general_data->stat_avg.vbendt+
                            general_data->stat_avg.vbendt_watts)*eu_conv);
    printf("Initial Torsion    energy  %g\n",general_data->stat_avg.vtorst*eu_conv);
    printf("Initial onefour    energy  %g\n",general_data->stat_avg.vonfot*eu_conv);
    printf("Initial surface      energy  %.10g\n",
                                     general_data->stat_avg.vsurft*eu_conv);
  }/*endif*/

    printf("Initial e-Hart+XC  energy  %g\n",(general_data->stat_avg.cp_ehart
				           + general_data->stat_avg.cp_exc)*eu_conv);
    printf("Initial e-Hart+muXC-XC  energy  %g\n",(general_data->stat_avg.cp_ehart
				           + general_data->stat_avg.cp_muxc
				           - general_data->stat_avg.cp_exc)*eu_conv);
    printf("Initial e-External energy  %g\n",general_data->stat_avg.cp_eext*eu_conv);
    printf("Initial e-Nonlocal energy  %g\n",general_data->stat_avg.cp_enl*eu_conv);
    printf("Initial e-Kinetic  energy  %g\n",general_data->stat_avg.cp_eke*eu_conv);
    printf("Initial volume             %g\n",general_data->stat_avg.vol
                                                       *BOHR*BOHR*BOHR);
    if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){

     printf("Initial Atm temperature    %.10g\n",
                   ((2.0*general_data->stat_avg.kinet*BOLTZ)/
	           ((double)(class->clatoms_info.nfree_pimd))));
    }/*endif*/
    printf("Initial e-Temperature      %g\n",general_data->stat_avg.kinet_cp
                                             *2.0*BOLTZ/c_div);
/*==========================================================================*/
  /* NHC quantities                                                        */

    if(general_data->ensopts.nvt==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc) );
      printf("Initial NHC temperature %.10g\n", 
	     ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      if(class->clatoms_info.pi_beads>1){
        nhc_div = (double) ( class->therm_info_bead.len_nhc*
                           (class->therm_info_bead.num_nhc)*
                           (class->clatoms_info.pi_beads-1) );
        printf("Initial bead NHC temperature %.10g\n", 
	       ((2.0*general_data->stat_avg.kinet_nhc_bead*BOLTZ)/nhc_div));
      }/*endif*/
    }/*endif*/

    if(general_data->ensopts.npt_i==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature %.10g\n", 
	     ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      if(class->clatoms_info.pi_beads>1){
       nhc_div = (double) ( class->therm_info_bead.len_nhc*
                           (class->therm_info_bead.num_nhc)*
                           (class->clatoms_info.pi_beads-1) );
       printf("Initial bead NHC temperature %.10g\n", 
	      ((2.0*general_data->stat_avg.kinet_nhc_bead*BOLTZ)/nhc_div));
      }
      printf("Initial Volume temperature %.10g\n",
	     (2.0*general_data->stat_avg.kinet_v*BOLTZ)); 
    }/*endif*/

    if(general_data->ensopts.npt_f==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature %.10g\n", 
	     ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      if(class->clatoms_info.pi_beads>1){
       nhc_div = (double) ( class->therm_info_bead.len_nhc*
                           (class->therm_info_bead.num_nhc)*
                           (class->clatoms_info.pi_beads-1) );
       printf("Initial bead NHC temperature %.10g\n", 
	      ((2.0*general_data->stat_avg.kinet_nhc_bead*BOLTZ)/nhc_div));
      }
      printf("Initial Volume temperature %.10g\n",      
	     ((2.0*general_data->stat_avg.kinet_v*BOLTZ)/6.0));
    }/*endif*/

   if(cp->cptherm_info.num_c_nhc != 0) {
    c_nhc_div = (double) ((cp->cptherm_info.len_c_nhc)*
                          (cp->cptherm_info.num_c_nhc));
    printf("%g\n",c_nhc_div);
    printf("CP NHC Temperature  = %g\n",general_data->stat_avg.kinet_nhc_cp
                                 *BOLTZ*2.0/c_nhc_div);
    } /* endif */

/*==========================================================================*/
  /* Neighbor list quantities                                              */

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
	printf("Number RESPA pairs = %d out of %d\n",npairs,np_tot);
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

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_fopen_cp_pimd(CLASS *class,GENERAL_DATA *general_data,
                           BONDED *bonded, CP *cp)

/*==========================================================================*/

{/* begin routine */

/*==========================================================================*/
/*   Local Variables                                                        */

  int n=1,iii,ibinary,iwrite_now;
  NAME file_typ;
  FILE *fp_bond_free,*fp_bend_free,*fp_tors_free,*fp_ccname;
  FILE *fp_iname, *fp_cpname, *fp_cvname,*fp_dname,*fp_centname;
  FILE *fp_kseigs;

/*==========================================================================*/
/*     A) Open dump file                                          */

  fp_dname = cfopen(general_data->filenames.dname,"w");

  if(fp_dname != NULL){fclose(fp_dname);}

/*==========================================================================*/
/*     A) Open inst avg file                                          */

  ibinary    = 0; /* never binary */
  iwrite_now = general_data->filenames.iwrite_inst;
  strcpy(file_typ,"ins_file");
  fp_iname   = cfopen(general_data->filenames.iname,"w");
  write_gen_header(class,general_data,fp_iname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_iname);

/*==========================================================================*/
/*     B) Open atm vel conf file                                       */

  ibinary    = general_data->filenames.iwrite_conf_binary;
  iwrite_now = general_data->filenames.iwrite_confv;
  strcpy(file_typ,"vel_file");
  fp_cvname  = cfopen(general_data->filenames.cvname,"w");
  write_gen_header(class,general_data,fp_cvname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_cvname); 

/*==========================================================================*/
/*     C) Open pos conf file                                           */

  ibinary    = general_data->filenames.iwrite_conf_binary;
  iwrite_now = general_data->filenames.iwrite_confp;
  strcpy(file_typ,"pos_file");
  fp_cpname  = cfopen(general_data->filenames.cpname,"w");
  write_gen_header(class,general_data,fp_cpname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_cpname); 

/*======================================================================*/
/*     C) Open partial pos conf file                                    */

 if((general_data->filenames.low_lim_par<=general_data->filenames.high_lim_par)){
   ibinary    = general_data->filenames.iwrite_conf_binary;
   iwrite_now = general_data->filenames.iwrite_par_confp;
   strcpy(file_typ,"par_file");
   fp_cpname = cfopen(general_data->filenames.cpparname,"w"); 
   write_gen_header(class,general_data,fp_cpname,ibinary,
                   iwrite_now,file_typ);
   fclose(fp_cpname); 
 }/*endif*/

/*==========================================================================*/
/*     C) Open centroid conf file                                           */

  ibinary    = general_data->filenames.iwrite_conf_binary;
  iwrite_now = general_data->filenames.iwrite_path_cent;
  strcpy(file_typ,"cen_file");
  fp_centname  = cfopen(general_data->filenames.centname,"w");
  write_gen_header(class,general_data,fp_centname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_centname); 

/*==========================================================================*/
/*     E) Open bond free energy file                                   */

    if(bonded->bond_free.num>0){
      fp_bond_free  = cfopen(bonded->bond_free.file,"w");
      fclose(fp_bond_free);
    }/*endif*/

/*==========================================================================*/
/*     F) Open bend free energy file                                    */

    if(bonded->bend_free.num>0){
      fp_bend_free  = cfopen(bonded->bend_free.file,"w");
      fclose(fp_bend_free);
    }/*endif*/

/*==========================================================================*/
/*     G) Open tors free energy file                                    */

    if(bonded->tors_free.num>0){
      fp_tors_free  = cfopen(bonded->tors_free.file,"w");
      fclose(fp_tors_free);
    }/*endif*/
    
/*==========================================================================*/
/*     H) Open CP conf file (Must be unformatted due to hugeness of         */
/*          wave functions)                                                 */

   ibinary    = 1;
   iwrite_now = general_data->filenames.iwrite_confc;
   strcpy(file_typ,"cof_file");
   fp_ccname = cfopen(general_data->filenames.ccname,"w");
   write_gen_header_cp(class,general_data,cp,fp_ccname,ibinary,
                       iwrite_now,file_typ);
   fclose(fp_ccname); 
    
/*==========================================================================*/
/*  I) Open CP KS eigenvalues file                                          */

  if(cp->cpcoeffs_info.ks_rot_on == 1){
    fp_kseigs = cfopen(general_data->filenames.ksname,"w");
    fprintf(fp_kseigs,
            "Nstate up,Nstate dn, LDA, LSDA, Nbeads,Dump freq,Ntime\n");
    fprintf(fp_kseigs,
            "%d  %d  %d  %d  %d  %d  %d\n",cp->cpcoeffs_info.nstate_up,
                                           cp->cpcoeffs_info.nstate_dn,
                                           cp->cpopts.cp_lda,
                                           cp->cpopts.cp_lsda,
                                           cp->cpcoeffs_info.pi_beads,
                                           cp->cpcoeffs_info.n_ks_rot,
                                           general_data->timeinfo.ntime);
    fclose(fp_kseigs);
  }/* endif */
    
/*==========================================================================*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dink_quant_calc_cp_pimd(CLASS *class, GENERAL_DATA *general_data,CP *cp,
                        double *etot,
                        double *econv_now,double *vtot, double *avtot,
                        double *deth,double *pnow, double *apnow,
                        double *nhc_div,double *nhc_div_bead,
                        double *c_nhc_div)
/*==========================================================================*/

{/* begin routine */

  int i,ncons,iii;

/*==========================================================================*/
 /*  Energy */

  *etot = 0.0;

     (*etot) = general_data->stat_avg.kinet_cp + general_data->stat_avg.cp_ehart + 
	       general_data->stat_avg.cp_eext + general_data->stat_avg.cp_exc + 
               general_data->stat_avg.cp_eke + general_data->stat_avg.cp_enl;
      if(cp->cptherm_info.num_c_nhc > 0) 
	(*etot) += general_data->stat_avg.kinet_nhc_cp 
	        + general_data->stat_avg.vpotnhc_cp;
      if(general_data->simopts.cp == 1||general_data->simopts.cp_pimd == 1) {
        (*etot) += (general_data->stat_avg.kinet) + (general_data->stat_avg.vintert) 
                +  (general_data->stat_avg.vintrat) + (general_data->stat_avg.kin_harm);
        if((general_data->ensopts.nvt)==1)  
 	 {(*etot)     += (general_data->stat_avg.kinet_nhc) 
                      + (general_data->stat_avg.kinet_nhc_bead) 
 	              +  (general_data->stat_avg.vpotnhc);}
        if((general_data->ensopts.npt_i)==1)
 	 {(*etot)     += (general_data->stat_avg.kinet_nhc) 
                      + (general_data->stat_avg.kinet_nhc_bead) 
 	              +  (general_data->stat_avg.vpotnhc)
 	              +  (general_data->stat_avg.kinet_v)   
 	              +  (general_data->stat_avg.vpot_v);}
       if((general_data->ensopts.npt_f)==1)
 	 {(*etot)+= (general_data->stat_avg.kinet_nhc) 
                 + (general_data->stat_avg.kinet_nhc_bead) 
 	         +  (general_data->stat_avg.vpotnhc)
 	         +  (general_data->stat_avg.kinet_v)   
 	         +  (general_data->stat_avg.vpot_v);}
      }/* endif  cp on*/


      if(((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==1)&&
         (general_data->stat_avg.iswit_vdw<=0))
         {(*etot)+= (general_data->stat_avg.vlong);}
      (*econv_now) = fabs((*etot)-general_data->stat_avg.econv0)/
         fabs(general_data->stat_avg.econv0);

 /*==========================================================================*/
 /*  Volume and pressure                                                    */

    if(general_data->simopts.cp == 1||general_data->simopts.cp_pimd == 1) {

      (*vtot) = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;
      (*avtot) = general_data->stat_avg.avintert + general_data->stat_avg.avintrat;
      (*deth) = getdeth(general_data->cell.hmat);
      for(i=1;i<=9;i++){
	general_data->stat_avg.apten_out[i]=(general_data->ptens.pvten_tot[i]
		         	    +general_data->ptens.tvten[i])/((*deth)*PCONV);}
      (*pnow)  = (general_data->stat_avg.apten_out[1]
	       +  general_data->stat_avg.apten_out[5]
	       +  general_data->stat_avg.apten_out[9])/(3.0);
      (*apnow) =  general_data->stat_avg.apress
	/(PCONV*((double)(general_data->timeinfo.itime)));
    }/*endif*/

 /*==========================================================================*/
 /*  NHC degrees of freedom                                                  */

      (*nhc_div) = (double)((class->therm_info_class.len_nhc)*
                   (class->therm_info_class.num_nhc));
      if(general_data->ensopts.npt_f == 1 || general_data->ensopts.npt_i == 1){
	(*nhc_div) = (double) (class->therm_info_class.len_nhc*
			      (class->therm_info_class.num_nhc +1) );
      }/*endif*/
      if(class->clatoms_info.pi_beads>1){
        *nhc_div_bead = (double) ( class->therm_info_bead.len_nhc*
                           (class->therm_info_bead.num_nhc)*
                           (class->clatoms_info.pi_beads-1) );
      }/*endif*/

/*==========================================================================*/
/*  Calculate c_nhc_div                                           */


   if(cp->cptherm_info.num_c_nhc != 0) {
     (*c_nhc_div) = (double) ((cp->cptherm_info.len_c_nhc)*
                              (cp->cptherm_info.num_c_nhc));
   }/* endif */

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void screen_write_cp_pimd(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,CP *cp,
                     double etot,double econv_now,double vtot, 
		     double avtot,double deth,double pnow, double apnow,
                     double nhc_div,double nhc_div_bead,
                     double a,double b,double c,
		     double tab,double tbc,double tac,double c_div,
                     double c_nhc_div)

/*==========================================================================*/
{/* begin routine */

 int npairs;
 int cp_lsda = cp->cpopts.cp_lsda;
 int anneal_opt = general_data->simopts.anneal_opt;
 static int count_norb=0;
 double atime,vol_div,atm_div;
 double updates_t,updates_true,updates_now; 

 double eelec;
 double aelec;

 double eu_conv=1.0;
 double dpi_beads = (double)(class->clatoms_info.pi_beads);
/*==========================================================================*/
/* Write to screen                                                          */

      atime = (double)(general_data->timeinfo.itime);
      if(general_data->filenames.iwrite_units==0){eu_conv=1.0;}
      if(general_data->filenames.iwrite_units==1){eu_conv=KCAL;}
      if(general_data->filenames.iwrite_units==2){eu_conv=BOLTZ;}

   if(general_data->timeinfo.itime==1) aelec=0.0;
   eelec =                   (general_data->stat_avg.cp_ehart
                            + general_data->stat_avg.cp_exc
                            + general_data->stat_avg.cp_eext
                            + general_data->stat_avg.cp_eke
                            + general_data->stat_avg.cp_enl);
   aelec =                  (general_data->stat_avg.acp_ehart
                            + general_data->stat_avg.acp_exc
                            + general_data->stat_avg.acp_eext
                            + general_data->stat_avg.acp_eke
                            + general_data->stat_avg.acp_enl);


/*==========================================================================*/
/*     A) Standard                                                   */
  printf("\n");
  printf("********************************************************\n");
  printf("QUANTITY          =  INSTANTANEOUS    AVERAGE           \n");
  printf("--------------------------------------------------------\n");
  if(general_data->ensopts.nve==1)  printf("Ensemble          = NVE     \n");
  if(general_data->ensopts.nvt==1)  printf("Ensemble          = NVT     \n");
  if(general_data->ensopts.npt_i==1)printf("Ensemble          = NPT-ISO \n");
  if(general_data->ensopts.npt_f==1)printf("Ensemble          = NPT-FLEX\n");
  printf("Time step         = %.10g\n",(double)(general_data->timeinfo.itime));
  printf("----------------- \n");
  printf("Econv             = %.10g %.10g\n",( econv_now ),
	 ( general_data->stat_avg.econv/atime));
  if(cp->cpopts.cp_isok_opt == 1)
     printf("CP KE conv        = %.10g %.10g \n",
	    fabs((general_data->stat_avg.kinet_cp - general_data->stat_avg.cp_kconv0)/
		 general_data->stat_avg.cp_kconv0),general_data->stat_avg.cp_kconv/atime);
  printf("----------------- \n");
  printf("Tot Energy        = %.10g %.10g\n",(general_data->stat_avg.kinet+vtot+eelec)*eu_conv,
         (general_data->stat_avg.akinet+avtot+aelec)/atime)*eu_conv;

  printf("----------------- \n");
  printf("Energy            = %.10g %.10g\n",(general_data->stat_avg.kinet+vtot)*eu_conv,
	 (general_data->stat_avg.akinet+avtot)/atime*eu_conv);
  printf("Total PE          = %.10g %.10g\n",(vtot)*eu_conv,(avtot/atime)*eu_conv);
  printf("Intermol PE       = %.10g %.10g\n",(general_data->stat_avg.vintert)*eu_conv,
	 (general_data->stat_avg.avintert/atime)*eu_conv);
  printf("Intramol PE       = %.10g %.10g\n",(general_data->stat_avg.vintrat)*eu_conv,
	 (general_data->stat_avg.avintrat/atime)*eu_conv);
  printf("Fict. Atm KE      = %.10g %.10g\n",(general_data->stat_avg.kinet)*eu_conv/dpi_beads,
	 (general_data->stat_avg.akinet/(atime*dpi_beads))*eu_conv);
  printf("Prim KE           = %.10g %.10g\n",(general_data->stat_avg.pi_ke_prim)*eu_conv,
	 (general_data->stat_avg.api_ke_prim/atime)*eu_conv);
  printf("Vir KE            = %.10g %.10g\n",(general_data->stat_avg.pi_ke_vir)*eu_conv,
	 (general_data->stat_avg.api_ke_vir/atime)*eu_conv);
  printf("----------------- \n");
  atm_div = (double)(class->clatoms_info.nfree);
  printf("Atm Deg. Free     = %.10g\n",atm_div); 
  atm_div = (double)(class->clatoms_info.nfree_pimd);
  printf("Fict. Atm Temperature   = %.10g %.10g \n",
	 (general_data->stat_avg.kinet*2.0*BOLTZ/(atm_div)),
	 (general_data->stat_avg.akinet*2.0*BOLTZ/(atm_div*atime)));

      printf("----------------- \n");
      printf("e-Energy          = %.10g %.10g\n",
                                (general_data->stat_avg.cp_ehart
			       + general_data->stat_avg.cp_exc
			       + general_data->stat_avg.cp_eext
			       + general_data->stat_avg.cp_eke
			       + general_data->stat_avg.cp_enl),
                                (general_data->stat_avg.acp_ehart
			       + general_data->stat_avg.acp_exc
			       + general_data->stat_avg.acp_eext
			       + general_data->stat_avg.acp_eke
			       + general_data->stat_avg.acp_enl)
                                 /atime);
      printf("e-Hartree + XC    = %.10g %.10g\n",(general_data->stat_avg.cp_ehart
				 	  + general_data->stat_avg.cp_exc),
                                           (general_data->stat_avg.acp_ehart
					  + general_data->stat_avg.acp_exc)
		                              /atime);
      printf("e-External PE     = %.10g %.10g\n",(general_data->stat_avg.cp_eext),
		                           (general_data->stat_avg.acp_eext)
		                              /atime);
      printf("e-Nonlocal PE     = %.10g %.10g\n",(general_data->stat_avg.cp_enl),
		                            (general_data->stat_avg.acp_enl)
		                              /atime);
      printf("e-Kinetic         = %.10g %.10g\n",(general_data->stat_avg.cp_eke),
		                            (general_data->stat_avg.acp_eke)
		                              /atime);
      printf("----------------- \n");
      printf("CP Fict KE        = %g %g\n",(general_data->stat_avg.kinet_cp),
                                           (general_data->stat_avg.akinet_cp/atime));

      printf("CP Temperature    = %.10g %.10g\n",(general_data->stat_avg.kinet_cp
					       *2.0*BOLTZ/c_div),
                                           (general_data->stat_avg.akinet_cp
 					       *2.0*BOLTZ/(c_div*atime)));

      if(anneal_opt == 0){
       if(general_data->stat_avg.kinet_cp*2.0*BOLTZ/c_div>10.0*cp->cpopts.te_ext){
         printf("\n$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
         printf("Your current CP temperature is greater than\n");
         printf("10 x the value set in the input file.\n");
         printf("Please reconsider your adiabaticity parameters\n");
         printf("and/or employ coefficient thermostats to control \n");
         printf("this quantity.\n");
         printf("$$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
         fflush(stdout);
       }/*endif*/
      }/* endif */

     printf("CP Up Force:   avg: %g max: %g\n",general_data->stat_avg.fc_mag_up,
                                             general_data->stat_avg.fc_max_up);
     if(cp_lsda==1){
      printf("CP Dn Force:   avg: %g max: %g\n",general_data->stat_avg.fc_mag_dn,
                                             general_data->stat_avg.fc_max_dn);
     }/*endif*/

  printf("----------------- \n");
/*==========================================================================*/
/*     B) Extended Class                                          */
     if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){
      if((general_data->ensopts.nvt + general_data->ensopts.npt_i
	  + general_data->ensopts.npt_f == 1) ){
	printf("NHC Temperature   = %.10g %.10g\n",
	       (general_data->stat_avg.kinet_nhc*2.0*BOLTZ/nhc_div),
	       (general_data->stat_avg.akinet_nhc*2.0*BOLTZ
		/(nhc_div*atime)));
        if(class->clatoms_info.pi_beads>1){
	 printf("Bead NHC Temperature   = %.10g %.10g\n",
	       (general_data->stat_avg.kinet_nhc_bead*2.0*BOLTZ/nhc_div_bead),
	       (general_data->stat_avg.akinet_nhc_bead*2.0*BOLTZ
		/(nhc_div_bead*atime)));
        }/*endif*/
      }/*endif*/
      if((general_data->ensopts.npt_i +general_data->ensopts.npt_f == 1)){
        vol_div = 1.0;
        if(general_data->ensopts.npt_f == 1) vol_div = 6.0;
	printf("Vol Temperature   = %.10g %.10g\n",
	       (general_data->stat_avg.kinet_v*2.0*BOLTZ/vol_div),
	       (general_data->stat_avg.akinet_v*2.0*BOLTZ/(atime*vol_div)));}
        printf("----------------- \n");
     if(cp->cptherm_info.num_c_nhc != 0) {
    printf("%.10g\n",c_nhc_div);
        printf("CP NHC Temp       = %.10g %.10g\n",
                (general_data->stat_avg.kinet_nhc_cp*2.0*BOLTZ)/c_nhc_div,
                (general_data->stat_avg.akinet_nhc_cp*2.0*BOLTZ)/(c_nhc_div*atime));
     } /* endif */
     } /* endif */
/*==========================================================================*/
/*     C)Pressure/Vol                                                      */
       if(((general_data->simopts.cp+general_data->simopts.cp_pimd) == 1)
          && cp->cpopts.cp_ptens_calc == 1
          &&(general_data->cell.iperd>=2)) {
       printf("Pressure          = %.10g %.10g\n",(pnow),( apnow));
       printf("Avg  P11,P22,P33  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten[1]/(PCONV*atime)),
              (general_data->stat_avg.apten[5]/(PCONV*atime)),
              (general_data->stat_avg.apten[9]/(PCONV*atime)));
       printf("Inst P11,P22,P33  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten_out[1]),
              (general_data->stat_avg.apten_out[5]),
              (general_data->stat_avg.apten_out[9]));
       printf("Avg  P12,P13,P23  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten[4])/(PCONV*atime),
              (general_data->stat_avg.apten[7])/(PCONV*atime),
              (general_data->stat_avg.apten[8])/(PCONV*atime));
       printf("Inst P12,P13,P23  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten_out[4]),
              (general_data->stat_avg.apten_out[7]),
              (general_data->stat_avg.apten_out[8]));
       printf("Avg  P21,P31,P32  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten[2])/(PCONV*atime),
              (general_data->stat_avg.apten[3])/(PCONV*atime),
              (general_data->stat_avg.apten[6])/(PCONV*atime));
       printf("Inst P21,P31,P32  = %.10g %.10g %.10g \n",
              (general_data->stat_avg.apten_out[2]),
              (general_data->stat_avg.apten_out[3]),
              (general_data->stat_avg.apten_out[6]));
       printf("----------------- \n");
       printf("Volume            = %.10g %.10g\n",(deth*BOHR*BOHR*BOHR),
              (general_data->stat_avg.avol/atime)*BOHR*BOHR*BOHR);
       printf("Inst cell lths    = %.10g %.10g %.10g\n",(a*BOHR),(b*BOHR),(c*BOHR));
       printf("Avg  cell lths    = %.10g %.10g %.10g\n",    
              (general_data->stat_avg.acella/atime)*BOHR,
              (general_data->stat_avg.acellb/atime)*BOHR,
              (general_data->stat_avg.acellc/atime)*BOHR);
       printf("Inst cell angs    = %.10g %.10g %.10g\n",(tab),(tac),(tbc));
       printf("Avg  cell angs    = %.10g %.10g %.10g\n",
              (general_data->stat_avg.acellab/atime),
              (general_data->stat_avg.acellac/atime),
              (general_data->stat_avg.acellbc/atime));
       printf("----------------- \n");
      }/* endif */

/*==========================================================================*/
/*     D)Constraint                                                  */

      if(bonded->constrnt.iconstrnt == 1) {
       if(general_data->simopts.cp==1||general_data->simopts.cp_pimd==1){
        if(bonded->bond.ncon > 0) {
         printf("Shake iter        = %.10g %.10g\n",
               (double)(general_data->stat_avg.iter_shake),
               (general_data->stat_avg.aiter_shake/atime));
         printf("Rattle iter       = %.10g %.10g\n",
               (double)(general_data->stat_avg.iter_ratl),
               (general_data->stat_avg.aiter_ratl/atime));
        }
        if(bonded->grp_bond_con.num_23 > 0) {
	 printf("Grp_23 shake iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_23,
		 general_data->stat_avg.aiter_23/atime);
         printf("Grp_23 ratl iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_23r,
		 general_data->stat_avg.aiter_23r/atime);
        }
        if(bonded->grp_bond_con.num_33 > 0) {
         printf("Grp_33 shake iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_33,
		 general_data->stat_avg.aiter_33/atime);
         printf("Grp_33 ratl iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_33r,
		 general_data->stat_avg.aiter_33r/atime);
        }
        if(bonded->grp_bond_con.num_46 > 0) {
         printf("Grp_46 shake iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_46,
		 general_data->stat_avg.aiter_46/atime);
         printf("Grp_46 ratl iter = %.10g %.10g\n",
		 general_data->stat_avg.iter_46r,
		 general_data->stat_avg.aiter_46r/atime);
        }
         printf("----------------- \n");
       }/*endif*/
      }/*endif*/
      if(!cp->cpopts.cp_norb) {
          printf("Orb shake iter    = %.10g %.10g\n",
                                    (double)(general_data->stat_avg.iter_shake_cp),
                                      (general_data->stat_avg.aiter_shake_cp/atime));
          printf("Orb rattle iter   = %.10g %.10g\n",
                                    (double) (general_data->stat_avg.iter_ratl_cp),
                                      (general_data->stat_avg.aiter_ratl_cp/atime));
      }
     if(cp->cpopts.cp_norb > 0) {
       if(cp->cpopts.cp_norb == 2 || cp->cpopts.cp_norb == 3)
          printf("Max off diag      = %.10g\n",
		                    cp->cpcoeffs_info.max_off_diag);
       if(cp->cpopts.cp_norb == 3)
          printf("Max diag          = %.10g\n",
		                    cp->cpcoeffs_info.max_diag);
      if(cp->cpcoeffs_info.max_off_diag > cp->cpconstrnt.c_tolnorb){
       ++count_norb;
      }
          printf("Rotations done    = %d\n",
		                    count_norb);
     }
/*==========================================================================*/
/*     D)Misc                                                         */

      if(class->nbr_list.iver == 1 
         && (general_data->simopts.cp+general_data->simopts.cp_pimd) == 1){
	  updates_t = general_data->stat_avg.updates;
	  if(updates_t == 0) updates_t = 1;
	  updates_now = (double ) (general_data->timeinfo.itime-
	 			 general_data->stat_avg.itime_update);
          updates_true = updates_t;
          npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
	       + class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
          printf("Inst steps/update = %.10g\n",updates_now);
          printf("Avg. steps/update = %.10g\n",(atime/updates_t));
          printf("Total list updates= %.10g\n",updates_true);
          printf("Number of pairs   = %d \n",npairs);
	  printf("----------------- \n");
      }/*endif*/
      printf("Cpu time          = %.10g %.10g\n",(general_data->stat_avg.cpu_now),
	     (general_data->stat_avg.acpu/atime));
      printf("--------------------------------------------------------\n");
      printf("\n"); 
      printf("********************************************************\n");
      fflush(stdout);

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_file_cp_pimd_class(CLASS *class,BONDED *bonded,
                                   GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/* begin routine */

 int i,j,ip,ncoef_up,ncoef_dn,izero,nwrite;
 int ichain,inhc,nstate;
 int pi_beads = class->clatoms_info.pi_beads;
 double deth;
 FILE *fp_dname;

/*==========================================================================*/
/* Open dump file                                                           */

  fp_dname = cfopen(general_data->filenames.dname,"o");

/*==========================================================================*/
/*     A)Atm positions                                               */

      fprintf(fp_dname,"natm_tot restart_typ itime\n");
      fprintf(fp_dname,"%d restart_all %d %d\n",class->clatoms_info.natm_tot,
	      general_data->timeinfo.itime,pi_beads);
      fprintf(fp_dname,"atm pos, atm_typ, mol_typ mol_num\n");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
	fprintf(fp_dname,"%g %g %g %s %s %s %d\n",
		class->clatoms_pos[ip].x[i],class->clatoms_pos[ip].y[i],
                                         class->clatoms_pos[ip].z[i],
		class->atommaps.atm_typ[class->atommaps.iatm_atm_typ[i]],
		class->atommaps.mol_typ[class->atommaps.iatm_mol_typ[i]],
		class->atommaps.res_typ[class->atommaps.iatm_res_typ[i]],
		class->atommaps.iatm_mol_num[i]);
       }/*endfor*/
      }/*endfor*/
/*==========================================================================*/
/*     B)Cell shape                                                 */

      if(general_data->ensopts.npt_i==0){
	deth = getdeth(general_data->cell.hmat);
	general_data->baro.x_lnv = log(deth)/3.0;
	general_data->baro.v_lnv = (general_data->par_rahman.vgmat[1]
			      + general_data->par_rahman.vgmat[5]
			      + general_data->par_rahman.vgmat[9])/3.0;
      }/*endif*/

      fprintf(fp_dname,"h matrix\n");
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat[1],
	      general_data->cell.hmat[4],general_data->cell.hmat[7]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat[2],
	      general_data->cell.hmat[5],general_data->cell.hmat[8]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat[3],
	      general_data->cell.hmat[6],general_data->cell.hmat[9]);
      fprintf(fp_dname,"h matrix for Ewald setup\n");
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_ewd[1],
	      general_data->cell.hmat_ewd[4],general_data->cell.hmat_ewd[7]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_ewd[2],
	      general_data->cell.hmat_ewd[5],general_data->cell.hmat_ewd[8]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_ewd[3],
	      general_data->cell.hmat_ewd[6],general_data->cell.hmat_ewd[9]);
      fprintf(fp_dname,"1/3 log(Vol)\n");
      fprintf(fp_dname,"%.12g\n",general_data->baro.x_lnv);

/*==========================================================================*/
/*    C)Atm and Atm NHC Velocities                                   */

      fprintf(fp_dname,"atm vel\n");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
	fprintf(fp_dname,"%g %g %g\n",class->clatoms_pos[ip].vx[i],
		class->clatoms_pos[ip].vy[i],class->clatoms_pos[ip].vz[i]);
       }/*endfor*/
      }/*endfor*/
      
      fprintf(fp_dname,"number of atm nhc, length of nhc\n");
      fprintf(fp_dname,"%d %d\n",class->therm_info_class.num_nhc,
	      class->therm_info_class.len_nhc);
      fprintf(fp_dname,"atm nhc velocities\n");
      for(j=1;j<=(class->therm_info_class.len_nhc);j++){
	for(i=1;i<=(class->therm_info_class.num_nhc);i++){
	  fprintf(fp_dname,"%g\n",class->therm_class.v_nhc[j][i]);
	} /*endfor*/
      }/*endfor*/

      fprintf(fp_dname,"number of bead nhc, length of nhc\n");
      fprintf(fp_dname,"%d %d\n",class->therm_info_bead.num_nhc,
	      class->therm_info_bead.len_nhc);
      fprintf(fp_dname,"bead nhc velocities\n");
      for(ip=2;ip<=pi_beads;ip++){
       for(j=1;j<=(class->therm_info_bead.len_nhc);j++){
	for(i=1;i<=(class->therm_info_bead.num_nhc);i++){
	  fprintf(fp_dname,"%g\n",class->therm_bead[ip].v_nhc[j][i]);
	} /*endfor*/
       }/*endfor*/
      }/*endfor*/

/*==========================================================================*/
/*    D)Vol and Vol NHC Velocities                             */

      fprintf(fp_dname,"vol velocities\n");
      fprintf(fp_dname,"%g %g %g\n",general_data->par_rahman.vgmat[1],
	      general_data->par_rahman.vgmat[4],general_data->par_rahman.vgmat[7]);
      fprintf(fp_dname,"%g %g %g\n",general_data->par_rahman.vgmat[2],
	      general_data->par_rahman.vgmat[5],general_data->par_rahman.vgmat[8]);
      fprintf(fp_dname,"%g %g %g\n",general_data->par_rahman.vgmat[3],
	      general_data->par_rahman.vgmat[6],general_data->par_rahman.vgmat[9]);
      fprintf(fp_dname,"log(vol) velocity\n");
      fprintf(fp_dname,"%g\n",general_data->baro.v_lnv);
      fprintf(fp_dname,"vol nhc velocities\n");
      for(i=1;i<=(class->therm_info_class.len_nhc);i++){
	fprintf(fp_dname,"%g\n",general_data->baro.v_vol_nhc[i]);
      }/*endfor*/

/*==========================================================================*/
/*    E)Misc                                                    */

      ncoef_up = cp->cpcoeffs_info.ncoef;
      ncoef_dn = cp->cpcoeffs_info.ncoef;
      if(cp->cpopts.cp_lda==1){ncoef_dn=0;}
      fprintf(fp_dname,"dt=%g\n",general_data->timeinfo.dt);
      fprintf(fp_dname,"nfree=%d\n",class->clatoms_info.nfree);
      fprintf(fp_dname,"nve=%d nvt=%d npt_i=%d npt_f=%d nst=%d\n",
	      general_data->ensopts.nve,  general_data->ensopts.nvt,  
	      general_data->ensopts.npt_i,  general_data->ensopts.npt_f,  
	      general_data->ensopts.nst);
      fprintf(fp_dname,"nbond_free=%d nbend_free=%d ntors_free=%d\n",
	      bonded->bond_free.num,bonded->bend_free.num,  
	      bonded->tors_free.num);
      fprintf(fp_dname,"t_ext=%g,pext=%g,stens_ext=%g\n",
              general_data->statepoint.t_ext,
	      general_data->statepoint.pext,
              general_data->statepoint.stens_ext);
      fprintf(fp_dname,"cp=%d cp_wave=%d cp_lda=%d cp_lsda=%d\n", 
                           general_data->simopts.cp,     general_data->simopts.cp_wave,  
		           cp->cpopts.cp_lda, cp->cpopts.cp_lsda);
      fprintf(fp_dname,"cp_sic=%d cp_norb=%d cp_nonint=%d cp_gga=%d\n",
                           cp->cpopts.cp_sic,    cp->cpopts.cp_norb,  
		           cp->cpopts.cp_nonint, cp->cpopts.cp_gga);
      fprintf(fp_dname,"nstate_up=%d nstate_dn=%d ncoef_up=%d ncoef_dn=%d\n",
                      cp->cpcoeffs_info.nstate_up,  cp->cpcoeffs_info.nstate_dn,
                      ncoef_up,ncoef_dn);
      fprintf(fp_dname,"vxc_typ=%s ggax_typ=%s\n",cp->pseudo.vxc_typ,
		                                     cp->pseudo.ggax_typ);
      fprintf(fp_dname,"ggac_typ=%s\n",cp->pseudo.ggac_typ);
      fflush(fp_dname); 
      fclose(fp_dname); 

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_config_files_cp_pimd_class(CLASS *class,BONDED *bonded,
                                      GENERAL_DATA *general_data,CP *cp)
 
/*==========================================================================*/

{/* begin routine */

  int i,ip,n,ncoef_up,ncoef_dn;
  int pi_beads = class->clatoms_info.pi_beads;
  FILE *fp_cpname, *fp_cvname,*fp_ccname, *fp_centname;

/*=====================================================================*/
/* 0) */
    ncoef_up = cp->cpcoeffs_info.ncoef;
    ncoef_dn = cp->cpcoeffs_info.ncoef;
    if(cp->cpopts.cp_lda==1){ncoef_dn=0;}

/*=====================================================================*/
/* I) Write to the atm config file                                     */ 

    if((general_data->timeinfo.itime % general_data->filenames.iwrite_confp) == 0 ){
     if(general_data->filenames.iwrite_conf_binary==0){
      fp_cpname  = cfopen(general_data->filenames.cpname,"a");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
	fprintf(fp_cpname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[ip].x[i],
		class->clatoms_pos[ip].y[i],class->clatoms_pos[ip].z[i]);
       }/*endfor*/
      }/*endfor*/
      for(i=1;i<=9;i+=3) 
	fprintf(fp_cpname,"%.13g %.13g %.13g\n",general_data->cell.hmat[i],
		general_data->cell.hmat[(i+1)],general_data->cell.hmat[(i+2)]);
      fflush(fp_cpname);
      fclose(fp_cpname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cpname = cfopen(general_data->filenames.cpname,"a");
    n=1;
   for(ip=1;ip<=pi_beads;ip++){
    for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
     fwrite(&(class->clatoms_pos[ip].x)[i],sizeof(double),n,fp_cpname);
     fwrite(&(class->clatoms_pos[ip].y)[i],sizeof(double),n,fp_cpname);
     fwrite(&(class->clatoms_pos[ip].z)[i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
   }/*endfor*/ 
    for(i=1;i<=9;i+=3){ 
      fwrite(&(general_data->cell.hmat)[i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[i+1],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[i+2],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
      fclose(fp_cpname);
   }/*endif*/
   }/*endif*/

  /*===================================================================*/
  /* I) Write to the partial atm config files                    */
  if((general_data->timeinfo.itime % general_data->filenames.iwrite_par_confp) == 0 ){

   if(general_data->filenames.iwrite_conf_binary==0){
    fp_cpname  = cfopen(general_data->filenames.cpparname,"a");
    for(ip=1;ip<=pi_beads;ip++){
     for(i=(general_data->filenames.low_lim_par);
       i<=(general_data->filenames.high_lim_par);i++){
       fprintf(fp_cpname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[ip].x[i],
	      class->clatoms_pos[ip].y[i],class->clatoms_pos[ip].z[i]);
     }/*endfor*/
    }/*endfor*/
    for(i=1;i<=9;i+=3) 
      fprintf(fp_cpname,"%.13g %.13g %.13g\n",general_data->cell.hmat[i],
	      general_data->cell.hmat[(i+1)],general_data->cell.hmat[(i+2)]);
    fflush(fp_cpname);
    fclose(fp_cpname);
   }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cpname = cfopen(general_data->filenames.cpparname,"a");
    n=1;
    for(ip=1;ip<=pi_beads;ip++){
     for(i=(general_data->filenames.low_lim_par);
      i<=(general_data->filenames.high_lim_par);i++){ 
      fwrite(&(class->clatoms_pos[ip].x)[i],sizeof(double),n,fp_cpname);
      fwrite(&(class->clatoms_pos[ip].y)[i],sizeof(double),n,fp_cpname);
      fwrite(&(class->clatoms_pos[ip].z)[i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
   }/*endfor*/ 
    for(i=1;i<=9;i+=3){ 
      fwrite(&(general_data->cell.hmat)[i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[i+1],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[i+2],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
      fclose(fp_cpname);
   }/*endif*/

  }/*endif*/
  /*=====================================================================*/

/*=====================================================================*/
 /* II) Write to the atm velocity config file                          */

  if((general_data->timeinfo.itime % general_data->filenames.iwrite_confv) == 0 ){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cvname  = cfopen(general_data->filenames.cvname,"a");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
	fprintf(fp_cvname,"%.12g  %.12g  %.12g\n",
                class->clatoms_pos[ip].vx[i],
		class->clatoms_pos[ip].vy[i],
                class->clatoms_pos[ip].vz[i]);
       }/*endfor*/
      }/*endfor*/
      for(i=1;i<=9;i+=3) 
	fprintf(fp_cvname,"%.13g %.13g %.13g\n",general_data->cell.hmat[i],
		general_data->cell.hmat[(i+1)],general_data->cell.hmat[(i+2)]);
      fflush(fp_cvname);
      fclose(fp_cvname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cvname = cfopen(general_data->filenames.cvname,"a");
    n=1;
  for(ip=1;ip<=pi_beads;ip++){
   for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
     fwrite(&(class->clatoms_pos[ip].vx)[i],sizeof(double),n,fp_cvname);
     fwrite(&(class->clatoms_pos[ip].vy)[i],sizeof(double),n,fp_cvname);
     fwrite(&(class->clatoms_pos[ip].vz)[i],sizeof(double),n,fp_cvname);
   }/*endfor*/ 
  }/*endfor*/ 
    for(i=1;i<=9;i+=3){ 
      fwrite(&(general_data->cell.hmat)[i],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[i+1],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[i+2],sizeof(double),n,fp_cvname);
    }/*endfor*/ 
      fclose(fp_cvname);
   }/*endif*/

  }/*endif*/
  
/*=====================================================================*/
 /* IV) Write to the centroid config file                             */

  if((general_data->timeinfo.itime % general_data->filenames.iwrite_path_cent) == 0 ){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_centname  = cfopen(general_data->filenames.centname,"a");
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
	fprintf(fp_centname,"%.12g  %.12g  %.12g\n",
                class->clatoms_pos[1].x[i],
		class->clatoms_pos[1].y[i],
                class->clatoms_pos[1].z[i]);
       }/*endfor*/
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
	fprintf(fp_centname,"%.12g  %.12g  %.12g\n",
                class->clatoms_pos[1].vx[i],
		class->clatoms_pos[1].vy[i],
                class->clatoms_pos[1].vz[i]);
       }/*endfor*/
      for(i=1;i<=9;i+=3) 
	fprintf(fp_centname,"%.13g %.13g %.13g\n",general_data->cell.hmat[i],
		general_data->cell.hmat[(i+1)],general_data->cell.hmat[(i+2)]);
      fflush(fp_centname);
      fclose(fp_centname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_centname = cfopen(general_data->filenames.centname,"a");
    n=1;
   for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
     fwrite(&(class->clatoms_pos[1].x)[i],sizeof(double),n,fp_centname);
     fwrite(&(class->clatoms_pos[1].y)[i],sizeof(double),n,fp_centname);
     fwrite(&(class->clatoms_pos[1].z)[i],sizeof(double),n,fp_centname);
   }/*endfor*/ 
   for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
     fwrite(&(class->clatoms_pos[1].vx)[i],sizeof(double),n,fp_centname);
     fwrite(&(class->clatoms_pos[1].vy)[i],sizeof(double),n,fp_centname);
     fwrite(&(class->clatoms_pos[1].vz)[i],sizeof(double),n,fp_centname);
   }/*endfor*/ 
    for(i=1;i<=9;i+=3){ 
      fwrite(&(general_data->cell.hmat)[i],sizeof(double),n,fp_centname);
      fwrite(&(general_data->cell.hmat)[i+1],sizeof(double),n,fp_centname);
      fwrite(&(general_data->cell.hmat)[i+2],sizeof(double),n,fp_centname);
    }/*endfor*/ 
      fclose(fp_centname);
   }/*endif*/

  }/*endif*/
  
/*==========================================================================*/
}/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_inst_file_cp_pimd(CLASS *class,GENERAL_DATA *general_data, CP *cp,
                        double etot,double a,double b,double c,
                        double tac,double tab,double tbc)

/*==========================================================================*/

{/* begin routine */

  int i;
  double inst_div;
  FILE *fp_iname;

/*==========================================================================*/
/*    A) Finish Averages                                            */

      inst_div    = (double)(general_data->filenames.iwrite_inst);    
      general_data->stat_avg.aikinet_cp     = general_data->stat_avg.aikinet_cp/inst_div;
      general_data->stat_avg.aikinet_nhc_cp = general_data->stat_avg.aikinet_nhc_cp
                                      /inst_div;
      general_data->stat_avg.aicp_ehart     = general_data->stat_avg.aicp_ehart/inst_div;
      general_data->stat_avg.aicp_eext      = general_data->stat_avg.aicp_eext/inst_div;
      general_data->stat_avg.aicp_exc       = general_data->stat_avg.aicp_exc/inst_div;
      general_data->stat_avg.aicp_eke       = general_data->stat_avg.aicp_eke/inst_div;
      general_data->stat_avg.aicp_enl       = general_data->stat_avg.aicp_enl/inst_div;
      general_data->stat_avg.aikinet     /= inst_div;
      general_data->stat_avg.aipi_ke_prim      /= inst_div;
      general_data->stat_avg.aipi_ke_vir       /= inst_div;
      general_data->stat_avg.aikin_harm        /= inst_div;
      general_data->stat_avg.aikinet_v   /= inst_div;
      general_data->stat_avg.aikinet_nhc /= inst_div;
      general_data->stat_avg.aikinet_nhc_bead /= inst_div;
      general_data->stat_avg.aivintert   /= inst_div;
      general_data->stat_avg.aivintrat   /= inst_div;
      general_data->stat_avg.aivol       /= inst_div;
      general_data->stat_avg.aicella     /= inst_div;
      general_data->stat_avg.aicellb     /= inst_div;
      general_data->stat_avg.aicellc     /= inst_div;
      general_data->stat_avg.aicellab    /= inst_div;
      general_data->stat_avg.aicellbc    /= inst_div;
      general_data->stat_avg.aicellac    /= inst_div;
      for(i=1;i<=9;i++){general_data->stat_avg.aipten[i] /= (inst_div*PCONV);}

/*==========================================================================*/
/*   B) Write Averages                                                */

      fp_iname = cfopen(general_data->filenames.iname,"a");
      fprintf(fp_iname,"\n");
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
	      general_data->stat_avg.aikinet,general_data->stat_avg.aipi_ke_prim,
              general_data->stat_avg.aipi_ke_vir,general_data->stat_avg.aikin_harm,
              general_data->stat_avg.aikinet_v,
              general_data->stat_avg.aikinet_nhc,general_data->stat_avg.aikinet_nhc_bead,
              general_data->stat_avg.aivintert,
	      general_data->stat_avg.aivintrat,general_data->stat_avg.aivol,etot);
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",
	      general_data->stat_avg.aicella,general_data->stat_avg.aicellb,
	      general_data->stat_avg.aicellc,general_data->stat_avg.aicellac,
	      general_data->stat_avg.aicellab,general_data->stat_avg.aicellbc);
      for(i=1;i<=9;i+=3) 
	fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.aipten[i],
		general_data->stat_avg.aipten[i+1],
		general_data->stat_avg.aipten[i+2]);
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
	      general_data->stat_avg.kinet,general_data->stat_avg.pi_ke_prim,
              general_data->stat_avg.pi_ke_vir,general_data->stat_avg.kin_harm,
              general_data->stat_avg.kinet_v,
	      general_data->stat_avg.kinet_nhc,general_data->stat_avg.kinet_nhc_bead,
              general_data->stat_avg.vintert,
	      general_data->stat_avg.vintrat,general_data->stat_avg.vol);
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",a,b,c,tac,tab,tbc);
      for(i=1;i<=9;i+=3) 
	fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.apten_out[i],
		general_data->stat_avg.apten_out[i+1],
		general_data->stat_avg.apten_out[i+2]);
       fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
               general_data->stat_avg.aikinet_cp,
               general_data->stat_avg.aikinet_nhc_cp,
               general_data->stat_avg.aicp_ehart,
               general_data->stat_avg.aicp_eext,
               general_data->stat_avg.aicp_exc,
               general_data->stat_avg.aicp_eke,
               general_data->stat_avg.aicp_enl);
      fflush(fp_iname);
      fclose(fp_iname);

/*==========================================================================*/
/*   C) Zero Averages                                               */

      general_data->stat_avg.aikinet_cp     = 0.0;
      general_data->stat_avg.aikinet_nhc_cp = 0.0;
      general_data->stat_avg.aicp_ehart     = 0.0;
      general_data->stat_avg.aicp_eext      = 0.0;
      general_data->stat_avg.aicp_exc       = 0.0;
      general_data->stat_avg.aicp_eke       = 0.0;
      general_data->stat_avg.aicp_enl       = 0.0;
      general_data->stat_avg.aikinet        = 0.0;
      general_data->stat_avg.aipi_ke_prim   = 0.0;
      general_data->stat_avg.aipi_ke_vir    = 0.0;
      general_data->stat_avg.aikin_harm     = 0.0;
      general_data->stat_avg.aikinet_v      = 0.0;
      general_data->stat_avg.aikinet_nhc    = 0.0;
      general_data->stat_avg.aikinet_nhc_bead    = 0.0;
      general_data->stat_avg.aivintert      = 0.0;
      general_data->stat_avg.aivintrat      = 0.0;
      general_data->stat_avg.aivol          = 0.0;
      general_data->stat_avg.aicella        = 0.0;
      general_data->stat_avg.aicellb        = 0.0;
      general_data->stat_avg.aicellc        = 0.0;
      general_data->stat_avg.aicellab       = 0.0;
      general_data->stat_avg.aicellbc       = 0.0;
      general_data->stat_avg.aicellac       = 0.0;
      for(i=1;i<=9;i++){general_data->stat_avg.aipten[i]  = 0.0;}

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*aaaaaa*/
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

  void write_dump_file_cp_pimd(CLASS *class,BONDED *bonded,
                               GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/* begin routine */

#include "../typ_defs/typ_mask.h"


 FILE *fp_dnamec;

 int i,j,ncoef_up,ncoef_dn,izero;
 int pi_beads = class->clatoms_info.pi_beads;
 int iproc,ioff_re,ioff_im;

 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo,iii;

 int ibinary = cp->cpopts.iwrite_coef_binary;
 
 MPI_Status stat;

/* Local pointers */
 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int nproc      = class->communicate.np;
 MPI_Comm world = class->communicate.world;


/*==========================================================================*/
/* Open cp dump file                                                        */


 if(myid == 0){
  fp_dnamec = cfopen(general_data->filenames.dnamec,"o");
 }

/*==========================================================================*/
/* Write the header                                                         */

 write_dump_header_cp_pimd(fp_dnamec,cp,general_data,myid,ibinary);

 if(nproc>1){Barrier(world);}

/*==========================================================================*/
/* write the occupation number                                              */

  write_dump_occ_cp_pimd(fp_dnamec,cp,myid,ibinary);

  if(nproc>1){Barrier(world);}

/*==========================================================================*/
/* Write the coefficients                                                   */

  write_dump_coef_cp_pimd(fp_dnamec,cp,class,ibinary);

/*==========================================================================*/
/* Write the coefficient velocities                                         */

  write_dump_vcoef_cp_pimd(fp_dnamec,cp,class,general_data,ibinary);

/*==========================================================================*/
/* Write the extended class stuff                                           */

  write_dump_extended_cp_pimd(fp_dnamec,cp,class,general_data,ibinary); 

/*==========================================================================*/
/* Done                                                                     */

 if(myid == 0){
   fflush(fp_dnamec); 
   fclose(fp_dnamec); 
 }


/*==========================================================================*/

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_conf_cp_pimd(CLASS *class,GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/* begin routine */

#include "../typ_defs/typ_mask.h"

 FILE *fp_ccname;


 int igo;
 int isoff,is,ip,iproc;
 int i,n,ncoef_up,ncoef_dn,ipoff;
 int cp_lsda = cp->cpopts.cp_lsda;
 int ncoef  = cp->cpcoeffs_info.ncoef;

 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 
 int upper,pflag;

 int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
 int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;

 MPI_Status stat;
 int pi_beads      = class->clatoms_info.pi_beads;
 int pi_beads_proc = class->clatoms_info.pi_beads_proc;


/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;

 int myid = class->communicate.myid;
 int nproc = class->communicate.np;
 MPI_Comm world = class->communicate.world;

/*=====================================================================*/
/* 0) */
    
    ncoef_up = cp->cpcoeffs_info.ncoef;
    ncoef_dn = cp->cpcoeffs_info.ncoef;
    if(cp->cpopts.cp_lda==1){ncoef_dn=0;}


/*==========================================================================*/
/* III) Write to the coefficient file                                       */
/*       Unformatted output is absolutely necessary here!!!!                */

     if((general_data->timeinfo.itime%
         general_data->filenames.iwrite_confc) == 0 ){
       if(myid==0) {fp_ccname = cfopen(general_data->filenames.ccname,"o");}
       if(nproc>1){Barrier(world);}
       n = ncoef;

  upper = pi_beads;
  pflag = 1;

  for(ip=1;ip<=upper;ip++){

       for(is=1;is<=nstate_up;is++){
         igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_up_st)*ncoef;
            for(i=1;i<=ncoef;i++){
              cre_up_tmp[i] = cp->cpcoeffs_pos[1].cre_up[(i+isoff)];
              cim_up_tmp[i] = cp->cpcoeffs_pos[1].cim_up[(i+isoff)];
            }/* endfor */
            igo = 1;
            if(myid!=0){
              Ssend(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
              Ssend(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
            }/*endif*/
         }/* endif */
         if(myid==0&&igo==0){
            Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,
                 MPI_ANY_TAG,world);
            Recv(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,
                 MPI_ANY_TAG,world);
          }/* endif */
          if(myid==0){
            fwrite(cre_up_tmp,sizeof(double),n,fp_ccname);
            fwrite(cim_up_tmp,sizeof(double),n,fp_ccname);
          }/* endif myid */
          if(nproc>1){Barrier(world);}
       }/*endfor:states*/
    }/*endfor:beads*/

 if(cp_lsda==1 && nstate_dn != 0){
    upper = pi_beads;
    pflag = 1;

  for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_dn;is++){
       igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_up_st)*ncoef;
             for(i=1;i<=ncoef;i++){
               cre_dn_tmp[i] = cp->cpcoeffs_pos[1].cre_dn[(i+isoff)];
               cim_dn_tmp[i] = cp->cpcoeffs_pos[1].cim_dn[(i+isoff)];
             }/* endfor */
             igo = 1;
             if(myid!=0){
               Ssend(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
               Ssend(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
             }/*endif*/
          }/* endif */
          if(myid==0&&igo==0){
             Recv(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,
                  MPI_ANY_TAG,world);
             Recv(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,
                  MPI_ANY_TAG,world);
           }/* endif */
           if(myid==0){
             fwrite(cre_dn_tmp,sizeof(double),n,fp_ccname);
             fwrite(cim_dn_tmp,sizeof(double),n,fp_ccname);
           }/* endif myid */
           if(nproc>1){Barrier(world);}
        }/*endfor:states*/
     }/*endfor:beads*/
       }/* endif lsda */
       if(myid==0){
        n=9;
        fwrite(general_data->cell.hmat,sizeof(double),n,fp_ccname);
        fclose(fp_ccname);
       }/* endif myid */
       if(nproc>1){Barrier(world);}
    }/*endif*/

/*==========================================================================*/

/*==========================================================================*/
/* Done                                                                     */

     if(nproc>1){Barrier(world);}

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_config_file_cp_pimd_quant(CP *cp,
                                   double *cre_temp, double *cim_temp,
                                   int ncoef,FILE *fp_ccname,int lsda,
                                   int ip,int pi_beads_proc,double *hmat)
 
/*==========================================================================*/

{/* begin routine */

  int i,n,ncoef_up,ncoef_dn;

/*=====================================================================*/
 /* III) Write to the coefficient file                                 */
/*       Unformatted output is absolutely necessary here!!!!           */


       if(lsda==0){
         n = (cp->cpcoeffs_info.nstate_up)*(ncoef);
       }else{
         n = (cp->cpcoeffs_info.nstate_dn)*(ncoef);
       }
       fwrite(cre_temp,sizeof(double),n,fp_ccname);
       fwrite(cim_temp,sizeof(double),n,fp_ccname);
       if((lsda==1)&&(ip==pi_beads_proc)){
         fwrite(hmat,sizeof(double),9,fp_ccname);
         fclose(fp_ccname);
       }

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_header_cp_pimd(FILE *fp_dnamec,CP *cp,GENERAL_DATA *general_data,int myid,
                               int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int izero,n;
 int csize = MAXWORD;
 int str_size,i;
 
 char *c_array1,*c_array2,*c_array3,*c_array4;
   
 if( (ibinary == 1) && (myid == 0) ){
   c_array1  = (char *) cmalloc(csize*sizeof(char));
   c_array2  = (char *) cmalloc(csize*sizeof(char));
   c_array3  = (char *) cmalloc(csize*sizeof(char));
   c_array4  = (char *) cmalloc(csize*sizeof(char));

/* initialize string  */
   for(i=0; i< csize-1; i++){
     c_array1[i] = ' ';
     c_array2[i] = ' ';
     c_array3[i] = ' ';
     c_array4[i] = ' ';
   }
     c_array1[csize-1] = '\0';
     c_array2[csize-1] = '\0';
     c_array3[csize-1] = '\0';
     c_array4[csize-1] = '\0';
 }/*endif*/


/*==========================================================================*/
/* Write the header                                                         */

if(myid==0){
 if(ibinary == 0){
    fprintf(fp_dnamec,"ncoef_up, ncoef_dn, nstate_up, nstate_dn, ");
    fprintf(fp_dnamec,"dft_typ, restart_typ, time of dump\n");
 }else{

#ifdef HP_VECLIB
   strcpy(c_array1,"hp_veclib:ncoef_up ");
#endif
#ifdef SGI_COMPLIB
   strcpy(c_array1,"sgi_complib:ncoef_up ");
#endif
#ifdef IBM_ESSL
   strcpy(c_array1,"ibm_essl:ncoef_up ");
#endif
#ifdef IBM_NOESSL
   strcpy(c_array1,"ibm_noessl:ncoef_up ");
#endif
#ifdef DEC_ALPHA
   strcpy(c_array1,"dec_alpha:ncoef_up ");
#endif
#ifdef _CRAY
   strcpy(c_array1,"cray:ncoef_up:ncoef_up ");
#endif
#ifdef SUN_COMPLIB
   strcpy(c_array1,"sun_complib:ncoef_up ");
#endif
#ifdef SUN_COMPLIB_OFF
   strcpy(c_array1,"sun_complib:ncoef_up ");
#endif
   strcpy(c_array2,"ncoef_dn, ");
   strcpy(c_array3,"nstate_up ");
   strcpy(c_array4,"nstate_dn ");

   fwrite(c_array1,sizeof(char),csize,fp_dnamec);
   fwrite(c_array2,sizeof(char),csize,fp_dnamec);
   fwrite(c_array3,sizeof(char),csize,fp_dnamec);
   fwrite(c_array4,sizeof(char),csize,fp_dnamec);

    strcpy(c_array1 ,"dft_typ,     ");
    strcpy(c_array2 ,"restart_typ  ");
    strcpy(c_array3 ,"time of dump ");

    fwrite(c_array1,sizeof(char),csize,fp_dnamec);
    fwrite(c_array2,sizeof(char),csize,fp_dnamec);
    fwrite(c_array3,sizeof(char),csize,fp_dnamec);  
 }/*endif ibinary*/


  if(cp->cpopts.cp_lda==1){
    izero = 0;
    if((general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd )==1){
      if(cp->cpopts.cp_norb>0){
       if(ibinary == 0){
       fprintf(fp_dnamec,"%d %d %d %d lda restart_all %d norb_on\n",
              cp->cpcoeffs_info.ncoef,izero,
              cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_up,
              general_data->timeinfo.itime);
       }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&izero,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);

      strcpy(c_array1,"lda");
      strcpy(c_array2,"restart_all");
      strcpy(c_array3,"norb_on");

      fwrite(c_array1,sizeof(char),csize,fp_dnamec);
      fwrite(c_array2,sizeof(char),csize,fp_dnamec);
      fwrite(&(general_data->timeinfo.itime),sizeof(int),n,fp_dnamec);
      fwrite(c_array3,sizeof(char),csize,fp_dnamec); 
      }/*endif ibinary*/

     }else{
      if(ibinary == 0){
       fprintf(fp_dnamec,"%d %d %d %d lda restart_all %d norb_off\n",
              cp->cpcoeffs_info.ncoef,izero,
              cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_up,
              general_data->timeinfo.itime);
      }else{
         n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&izero,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lda");
         strcpy(c_array2,"restart_all");
         strcpy(c_array3,"norb_off");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary*/
     }/*endif*/
    }else{

     if(cp->cpopts.cp_norb>0){

      if(ibinary == 0){
       fprintf(fp_dnamec,"%d %d %d %d lda restart_pos %d norb_on\n",
              cp->cpcoeffs_info.ncoef,izero,
              cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_up,
              general_data->timeinfo.itime);
      }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&izero,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);

	 strcpy(c_array1,"lda");
         strcpy(c_array2,"restart_pos");
         strcpy(c_array3,"norb_on");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec); 
      }/*endif ibinary*/
      }else{
     if(ibinary == 0){
       fprintf(fp_dnamec,"%d %d %d %d lda restart_pos %d norb_off\n",
              cp->cpcoeffs_info.ncoef,izero,
              cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_up,
              general_data->timeinfo.itime);
     }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&izero,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lda");
         strcpy(c_array2,"restart_pos");
         strcpy(c_array3,"norb_off");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
     }/*endif ibinary*/

      }/*endif*/
    }/*endif*/
  }/*endif*/
  if(cp->cpopts.cp_lsda==1){
    if((general_data->simopts.cp_pimd + general_data->simopts.cp_wave_pimd)==1){
     if(cp->cpopts.cp_norb>0){
      if(ibinary == 0){
      fprintf(fp_dnamec,"%d %d %d %d lsda restart_all %d norb_on\n",
         cp->cpcoeffs_info.ncoef,cp->cpcoeffs_info.ncoef,
         cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_dn,
         general_data->timeinfo.itime);
      }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_dn,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lsda");
         strcpy(c_array2,"restart_all");
         strcpy(c_array3,"norb_on");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary*/

     }else{
      if( ibinary == 0){
      fprintf(fp_dnamec,"%d %d %d %d lsda restart_all %d norb_off\n",
         cp->cpcoeffs_info.ncoef,cp->cpcoeffs_info.ncoef,
         cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_dn,
         general_data->timeinfo.itime);
      }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_dn,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lsda");
         strcpy(c_array2,"restart_all");
         strcpy(c_array3,"norb_off");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary */
     }/*endif*/
    }else{
     if(cp->cpopts.cp_norb>0){
      if(ibinary == 0){
      fprintf(fp_dnamec,"%d %d %d %d lsda restart_pos %d norb_on\n",
         cp->cpcoeffs_info.ncoef,cp->cpcoeffs_info.ncoef,
         cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_dn,
         general_data->timeinfo.itime);
      }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_dn,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lsda");
         strcpy(c_array2,"restart_pos");
         strcpy(c_array3,"norb_on");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary*/
     }else{
      if(ibinary == 0){
      fprintf(fp_dnamec,"%d %d %d %d lsda restart_pos %d norb_off\n",
         cp->cpcoeffs_info.ncoef,cp->cpcoeffs_info.ncoef,
         cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_dn,
         general_data->timeinfo.itime);
      }else{
        n = 1;
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.ncoef,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_up,sizeof(int),n,fp_dnamec);
        fwrite(&cp->cpcoeffs_info.nstate_dn,sizeof(int),n,fp_dnamec);

         strcpy(c_array1,"lsda");
         strcpy(c_array2,"restart_pos");
         strcpy(c_array3,"norb_off");

       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       fwrite(&general_data->timeinfo.itime,sizeof(int),n,fp_dnamec);
       fwrite(c_array3,sizeof(char),csize,fp_dnamec);
      }/*endif ibinary*/
     }/*endif*/
    }/*endif*/
  }/*endif*/
}/*endif:myid==0*/

 if( (ibinary == 1) && (myid == 0)){
   cfree(c_array1);
   cfree(c_array2);
   cfree(c_array3); 
   cfree(c_array4);  
 }/*endif*/


/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_occ_cp_pimd(FILE *fp_dnamec,CP *cp,int myid,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int i,n;
  int nstate;
  double occ_temp; /* for binary write option */

  char *c_array1,*c_array2;
  int csize = MAXWORD;


  if( (ibinary == 1) && (myid == 0) ){
    c_array1 = (char *) cmalloc(csize*sizeof(char ));
    c_array2 = (char *) cmalloc(csize*sizeof(char ));
  /* initialize arrays  */
    for(i=0 ; i< csize-1; i++){
      c_array1[i] = ' ';
      c_array2[i] = ' ';
    }
      c_array1[csize -1] = '\0';
      c_array2[csize -1] = '\0';
  }/*endif*/

/*==========================================================================*/
/* write the occupation number                                              */

 if(myid == 0){

  nstate  = MAX(cp->cpcoeffs_info.nstate_up,cp->cpcoeffs_info.nstate_dn);

  if(ibinary == 0){
   fprintf(fp_dnamec,"up occupation    dn occupation\n");
  }else{
    strcpy(c_array1,"up occupation");
    strcpy(c_array2,"dn occupation");

    fwrite(c_array1,sizeof(char),csize,fp_dnamec);
    fwrite(c_array2,sizeof(char),csize,fp_dnamec);

  }/*endif ibinary*/

  for(i=1;i<=nstate;i++){
   if(cp->cpopts.cp_lsda==1){
    if(ibinary == 0){
    fprintf(fp_dnamec,"%g %g\n",cp->cpopts.occ_up[i],cp->cpopts.occ_dn[i]);
    }else{
     n = 1;
     fwrite(&(cp->cpopts.occ_up[i]),sizeof(double),n,fp_dnamec);
     fwrite(&(cp->cpopts.occ_dn[i]),sizeof(double),n,fp_dnamec);
    }/*endif ibinary */

   }else{
    if(ibinary == 0){
    fprintf(fp_dnamec,"%g %g\n",(cp->cpopts.occ_up[i]-cp->cpopts.occ_dn[i]),
                                cp->cpopts.occ_dn[i]);
    }else{
     n = 1;
     occ_temp = (cp->cpopts.occ_up[i]-cp->cpopts.occ_dn[i]);
     fwrite(&occ_temp,sizeof(double),n,fp_dnamec);
     fwrite(&(cp->cpopts.occ_dn[i]),sizeof(double),n,fp_dnamec);

    }/*endif ibinary*/
   }/*endif*/
  }/*endfor*/

}/*endif:myid==0*/

/* free locally assigned memory  */
 if( (myid == 0) && (ibinary == 1)){
        cfree(c_array1);
        cfree(c_array2);  
 }

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_coef_cp_pimd(FILE *fp_dnamec,CP *cp,CLASS *class,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int i,j,izero,nwrite,pflag,lsda_flag,n,iii;
 int ichain,inhc,nstate,source,source_old;
 double deth;
 double *cre_up_temp,*cim_up_temp,*cre_dn_temp,*cim_dn_temp,*vc_nhc_temp;
 int pi_beads      = class->clatoms_info.pi_beads;
 int pi_beads_proc = class->clatoms_info.pi_beads_proc;
 int lsda          = cp->cpopts.cp_lsda;
 int myid          = class->communicate.myid;

 int upper;
 int igo;
 int isoff,ipoff,is,ip;

 int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
 int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
 int  ncoef  = cp->cpcoeffs_info.ncoef;
 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv;



 float cre_dum,cim_dum;

 char *c_array1,*c_array2;  /*for binary write */
 int csize = MAXWORD; 

 MPI_Comm world    = class->communicate.world;
 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;


/*==========================================================================*/
/* Initialize character arrays                                              */

  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(csize*sizeof(char ));
    c_array2 = (char *) cmalloc(csize*sizeof(char ));
  /* initialize array */
    for(i=0; i< csize-1; i++){
     c_array1[i] = ' ';
     c_array2[i] = ' ';
    }
     c_array1[csize -1] = '\0';
     c_array2[csize -1] = '\0';

  }/*endif */

/*==========================================================================*/
/*  I) Write Up states                                                      */

   if(myid==0){
     if(ibinary == 0){
       fprintf(fp_dnamec,"Up_State_Real         Up_State_Imag\n");
     }else{
       strcpy(c_array1,"Up_State_Real ");
       strcpy(c_array2,"Up_State_Imag ");
       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
     }/*endif ibinary */
   }/*endif myid*/
   Barrier(world);
/*-------------------------------------------------------------------------*/
/* Up states                                                               */

  upper = pi_beads;
  pflag = 1;

  for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_up;is++){
       igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_up_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_up_tmp[i] = cp->cpcoeffs_pos[ipoff].cre_up[(i+isoff)];
            cim_up_tmp[i] = cp->cpcoeffs_pos[ipoff].cim_up[(i+isoff)];
          }/* endfor */
          igo = 1;
          if(myid!=0){
            Ssend(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
            Ssend(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
          }/*endif*/
       }/* endif */
       if(myid==0&&igo==0){
         Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
         Recv(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
         n = 1;
         for(i=1;i<=ncoef;i++){
         if(ibinary == 0){
          fprintf(fp_dnamec," %.10e  %.10e \n",cre_up_tmp[i],cim_up_tmp[i]);
         }else{
           cre_dum = (float)cre_up_tmp[i];
           cim_dum = (float)cim_up_tmp[i];
           fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
           fwrite(&(cim_dum),sizeof(float),n,fp_dnamec);
         }/*endif ibinary */
         }/*endfor: write*/
        }/* endif myid */
        Barrier(world);
      }/*endfor:states*/
    }/*endfor:beads*/

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/* Dn states                                                               */

  if(cp->cpopts.cp_lsda==1){
   if(myid==0){
     if(ibinary == 0){
      fprintf(fp_dnamec,"Dn_State_Real         Dn_State_Imag\n");
     }else{
       strcpy(c_array1,"Dn_State_Real ");
       strcpy(c_array2,"Dn_State_Imag ");
       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
     }/*endif ibinary */
   }/*endif myid*/
   Barrier(world);

  upper = pi_beads;
  pflag = 1;

  for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_dn;is++){
       igo=0;
        if((is>=istate_dn_st) && (is<=istate_dn_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_dn_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_dn_tmp[i] = cp->cpcoeffs_pos[ipoff].cre_dn[(i+isoff)];
            cim_dn_tmp[i] = cp->cpcoeffs_pos[ipoff].cim_dn[(i+isoff)];
          }/* endfor */
          igo = 1;
          if(myid!=0){
            Ssend(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
            Ssend(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
          }/*endif*/
       }/* endif */
       if(myid==0&&igo==0){
         Recv(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
         Recv(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
          n = 1;
         for(i=1;i<=ncoef;i++){
        if(ibinary == 0){
          fprintf(fp_dnamec,"%.8g %.8g\n",(float)cre_dn_tmp[i],(float)cim_dn_tmp[i]);
         }else{
           cre_dum = (float) cre_dn_tmp[i];
           cim_dum = (float) cim_dn_tmp[i];
           fwrite(&(cre_dum),sizeof(double),n,fp_dnamec);
           fwrite(&(cim_dum),sizeof(double),n,fp_dnamec);
         }/*endif ibinary */
         }/*endfor: write*/
        }/* endif myid */
        Barrier(world);
     }/*endfor:states*/
   }/*endfor:beads*/
  }/*endif:lsda*/

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
        cfree(c_array1);
        cfree(c_array2); 
  }/*endif */

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_vcoef_cp_pimd(FILE *fp_dnamec,CP *cp,CLASS *class,
                         GENERAL_DATA *general_data,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int  ncoef  = cp->cpcoeffs_info.ncoef;
 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 int  cp_lsda       = cp->cpopts.cp_lsda;

 int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
 int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
 int pi_beads      = class->clatoms_info.pi_beads;
 int pi_beads_proc = class->clatoms_info.pi_beads_proc;

 int nrecv,isoff,ipoff;
 int is,igo,ip,pflag,upper;
 int i,n; 
 float cre_dum,cim_dum; 



 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;

 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 MPI_Comm world = class->communicate.world;

 char *c_array1,*c_array2;
 int csize = MAXWORD;


  if( (myid == 0) && (ibinary == 1)){
    c_array1 = (char *) cmalloc(csize*sizeof(char ));
    c_array2 = (char *) cmalloc(csize*sizeof(char ));
  /* initialize array */
   for(i=0 ; i< csize-1; i++){
    c_array1[i] = ' ';
    c_array2[i] = ' ';
   }
    c_array1[csize-1] = '\0';
    c_array2[csize-1] = '\0';
  }/*endif*/

/*==========================================================================*/
/* Write the coefficient velocities                                         */

/*-------------------------------------------------------------------------*/
/* Up states                                                               */

 if((general_data->simopts.cp_pimd+general_data->simopts.cp_wave_pimd)==1){

   if(myid==0){
     if(ibinary == 0){
     fprintf(fp_dnamec,"V_Up_State_Real         V_Up_State_Imag\n");
     }else{
       strcpy(c_array1,"V_UP_State_Real");
       strcpy(c_array2,"V_UP_State_Imag");
       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
     }/*endif ibinary */
   }/*endif myid*/
   Barrier(world);

  upper = pi_beads;
  pflag = 1;

  for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_up;is++){
       igo=0;
        if((is>=istate_up_st) && (is<=istate_up_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_up_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_up_tmp[i] = cp->cpcoeffs_pos[ipoff].vcre_up[(i+isoff)];
            cim_up_tmp[i] = cp->cpcoeffs_pos[ipoff].vcim_up[(i+isoff)];
          }/* endfor */
          igo = 1;
          if(myid!=0){
            Ssend(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
            Ssend(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
          }/*endif*/
       }/* endif */
       if(myid==0&&igo==0){
         Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
         Recv(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
         for(i=1;i<=ncoef;i++){
        if(ibinary == 0){
          fprintf(fp_dnamec,"%.8g %.8g\n",cre_up_tmp[i],cim_up_tmp[i]);
        }else{
          n = 1;
         cre_dum = (float) cre_up_tmp[i];
         cim_dum = (float) cim_up_tmp[i];
          fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
          fwrite(&(cim_dum),sizeof(float),n,fp_dnamec);
        }/*endif ibinary*/
         }/*endfor: write*/
        }/* endif myid */
        Barrier(world);
     }/*endfor:states*/
  }/*endfor beads*/

/*-------------------------------------------------------------------------*/
/* Dn states                                                               */

  if(cp->cpopts.cp_lsda==1){
    if(myid==0){
     if(ibinary == 0){
     fprintf(fp_dnamec,"V_Dn_State_Real         V_Dn_State_Imag\n");
     }else{
       strcpy(c_array1,"V_Dn_State_Real");
       strcpy(c_array2,"V_Dn_State_Imag");
       fwrite(c_array1,sizeof(char),csize,fp_dnamec);
       fwrite(c_array2,sizeof(char),csize,fp_dnamec);
     }/*endif ibinary*/
    }/*endif myid*/
    Barrier(world);


  upper = pi_beads;
  pflag = 1;

  for(ip=1;ip<=upper;ip++){
     for(is=1;is<=nstate_dn;is++){
       igo=0;
        if((is>=istate_dn_st) && (is<=istate_dn_end) && 
           (ip>=pi_beads_st)  && (ip<=pi_beads_end)){
          ipoff=(ip-pi_beads_st+1);
          isoff=(is-istate_dn_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_dn_tmp[i] = cp->cpcoeffs_pos[ipoff].vcre_dn[(i+isoff)];
            cim_dn_tmp[i] = cp->cpcoeffs_pos[ipoff].vcim_dn[(i+isoff)];
          }/* endfor */
          igo = 1;
          if(myid!=0){
            Ssend(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
            Ssend(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,0,0,world);
          }/*endif*/
       }/* endif */
       if(myid==0&&igo==0){
         Recv(&(cre_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
         Recv(&(cim_dn_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
         for(i=1;i<=ncoef;i++){
         if(ibinary == 0){
          fprintf(fp_dnamec,"%.8g %.8g\n",cre_dn_tmp[i],cim_dn_tmp[i]);
         }else{
           n = 1;
           cre_dum = (float) cre_dn_tmp[i];
           cim_dum = (float) cim_dn_tmp[i];
           fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
           fwrite(&(cim_dum),sizeof(float),n,fp_dnamec);
         }/*endif ibinary*/
         }/*endfor: write*/
        }/* endif myid */
        Barrier(world);
     }/*endfor:states*/
   }/*endfor:beads*/
  }/*endif:lsda*/

 }/*endif:write the velocities*/

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
     c_array1 = (char *) cmalloc(csize*sizeof(char ));
     c_array2 = (char *) cmalloc(csize*sizeof(char ));  
  }/*endif*/

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_extended_cp_pimd(FILE *fp_dnamec,CP *cp,CLASS *class,
                            GENERAL_DATA *general_data,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int iii;
 int i,j,ncoef_up,ncoef_dn;
 int ichain,inhc ;
 int iproc,ioff_re,ioff_im;
 int nstate_proc_now;
 int ip,ipoff;

 int  ncoef  = cp->cpcoeffs_info.ncoef;

 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;

 int  num_c_nhc_proc=cp->cptherm_info.num_c_nhc_proc;
 int  num_c_nhc     =cp->cptherm_info.num_c_nhc;
 int  len_c_nhc     =cp->cptherm_info.len_c_nhc;
 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo;

 int pi_beads_st        = cp->cpcoeffs_info.pi_beads_proc_st;
 int pi_beads_end       = cp->cpcoeffs_info.pi_beads_proc_end;
 int pi_beads      = class->clatoms_info.pi_beads;
 int pi_beads_proc = class->clatoms_info.pi_beads_proc;
 int myid_st          = cp->communicate.myid_state;
 int myid_bead          = class->communicate.myid_bead;

 float cre_dum,cim_dum;
 
 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;
 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int np_beads   = class->communicate.np_beads;
 MPI_Comm world = class->communicate.world;


 char *c_array1,*c_array2;
 int csize = MAXWORD;
 int n;
 int upper;

 int nsend;
         
  if( (myid == 0) && (ibinary == 1)){
     c_array1 = (char *) cmalloc(csize*sizeof(char ));
     c_array2 = (char *) cmalloc(csize*sizeof(char ));
  for(i=0; i< csize -1; i++){
     c_array1[i] = ' ';
     c_array2[i] = ' ';
  }  
     c_array1[csize -1] = '\0';
     c_array2[csize -1] = '\0';
  }/*endif*/


/*==========================================================================*/
/* Write the extended class stuff                                         */

  if((general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd)==1){
   if(myid==0){
    if(ibinary == 0){
    fprintf(fp_dnamec,"Number of chains    Length of chains \n");
    fprintf(fp_dnamec,"%d %d\n",cp->cptherm_info.num_c_nhc,
                                cp->cptherm_info.len_c_nhc);
    fprintf(fp_dnamec,"Chain velocities\n");
    }else{
     strcpy(c_array1,"Number of chains ");
     strcpy(c_array2,"Length of chains ");
     fwrite(c_array1,sizeof(char),csize,fp_dnamec);
     fwrite(c_array2,sizeof(char),csize,fp_dnamec);
       n = 1;
     fwrite(&(cp->cptherm_info.num_c_nhc),sizeof(int),n,fp_dnamec);
     fwrite(&(cp->cptherm_info.len_c_nhc),sizeof(int),n,fp_dnamec);
     strcpy(c_array1,"Chain velocities ");
     fwrite(c_array1,sizeof(char),csize,fp_dnamec);
   }/*endif ibinary */
   }/* endif : myid==0 */

/*==========================================================================*/


/*-----------------------------------------------------------------------*/
/* ii) Normal */

   if(cp->cptherm_info.istate_nhc_opt<=3){
          
    for(ichain=1;ichain<=len_c_nhc;ichain++){
       iii = 0;
     for(ip=1; ip<= pi_beads; ip++) {
     if(ip>=pi_beads_st && ip <= pi_beads_end ){
        ipoff=(ip-pi_beads_st+1);
      for(inhc=1;inhc<= num_c_nhc;inhc++){
         iii++;
        cim_up_tmp[iii] = cp->cptherm_pos[ipoff].vc_nhc[ichain][inhc];
      }
      } 
     } 

   nsend = pi_beads_proc*num_c_nhc;

   if(np_beads > 1 && myid_st==0){

   Gather(&(cim_up_tmp[1]),nsend,MPI_DOUBLE,
          &(cre_up_tmp[1]),nsend,MPI_DOUBLE,
          0,class->communicate.comm_beads);
   }else{
     for(i=1; i<= nsend ; i++){
        cre_up_tmp[i] = cim_up_tmp[i];
     }/*endfor*/
   }/*endif*/



  if(myid_bead == 0){ 
   for(i=1; i<= (num_c_nhc*pi_beads); i++){
     if(ibinary == 0){
       fprintf(fp_dnamec,"%.8g\n",cre_up_tmp[i]);
     }else{
       n = 1;
       cre_dum = (float) cre_up_tmp[i];
       fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
     }/*endif ibinary*/

   }
  }/*endif myid_bead*/

      
    
    }

   Barrier(world);

  }/*endif : istate_nhc_opt==3*/
/*-----------------------------------------------------------------------*/


   for(ichain=1;ichain<=len_c_nhc;ichain++){
     for(ip=1;ip<=pi_beads;ip++){

/*--------------------------------------------------------------------------*/
/* ii) Massive */

     if(cp->cptherm_info.istate_nhc_opt==4){


       for(iproc=0;iproc<np_states;iproc++){
         nstate_proc_now=0;
        if(ip>=pi_beads_st && ip <= pi_beads_end && myid_st == iproc){
           if(myid!=0){Ssend(&(nstate_up_proc),1,MPI_INT,0,0,world);}
           nstate_proc_now=nstate_up_proc;
        }/*endif*/
         if(myid==0&&iproc!=0){         
           Recv(&(nstate_proc_now),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }/*endif*/

        for(is=1;is<=nstate_proc_now;is++){
           if(myid==iproc){
             ipoff=(ip-pi_beads_st+1);
             ioff_re = (is-1)*ncoef;
             ioff_im = (nstate_proc_now+is-1)*ncoef;
             for(inhc=1;inhc<=ncoef;inhc++){
               cre_up_tmp[inhc] = 
                          cp->cptherm_pos[ipoff].vc_nhc[ichain][inhc+ioff_re];
             }/* endfor */
             for(inhc=1;inhc<=ncoef;inhc++){
               cre_up_tmp[inhc+ncoef] = 
                           cp->cptherm_pos[ipoff].vc_nhc[ichain][inhc+ioff_im];
             }/* endfor */
             if(myid!=0){
               Ssend(&(cre_up_tmp[0]),(2*ncoef+1),MPI_DOUBLE,0,0,world);
             } /* endif */
           }/* endif myid==iproc */
           if(myid==0 && iproc != 0){
             Recv(&(cre_up_tmp[0]),(2*ncoef+1),MPI_DOUBLE,
                                        MPI_ANY_SOURCE,MPI_ANY_TAG,world);
           }/* endif */
           if(myid==0){
             for(inhc=1;inhc<=2*ncoef;inhc++){
	      if(ibinary == 0){
               fprintf(fp_dnamec,"%.8g\n",cre_up_tmp[inhc]);
              }else{
               n = 1;
               cre_dum = (float) cre_up_tmp[inhc];
               fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);  
              }/*endif ibinary */
             }/*endfor inhc*/
           }/* endif myid*/
        }/*endfor : is*/

         if(cp_lsda==1){
           nstate_proc_now=0;
           if(myid==iproc){         
             if(myid!=0){Ssend(&(nstate_dn_proc),1,MPI_INT,0,0,world);}
             nstate_proc_now=nstate_dn_proc;
           }/*endif*/
           if(myid==0&&iproc!=0){         
             Recv(&(nstate_proc_now),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
                  world);
           }/*endif*/
           for(is=1;is<=nstate_proc_now;is++){
             if(myid==iproc){
               ipoff=(ip-pi_beads_st+1);
               ioff_re = (is-1+nstate_up_proc*2)*ncoef;
               ioff_im = (nstate_proc_now+is-1+nstate_up_proc*2)*ncoef;
               for(inhc=1;inhc<=ncoef;inhc++){
                 cre_up_tmp[inhc] = 
                     cp->cptherm_pos[ipoff].vc_nhc[ichain][inhc+ioff_re];
               }/* endfor */
               for(inhc=1;inhc<=ncoef;inhc++){
                 cre_up_tmp[inhc+ncoef] = 
                           cp->cptherm_pos[ipoff].vc_nhc[ichain][inhc+ioff_im];
               }/* endfor */
               if(myid!=0){
                 Ssend(&(cre_up_tmp[0]),(2*ncoef+1),MPI_DOUBLE,0,0,world);
               } /* endif */
             }/* endif myid==iproc */
             if(myid==0 && iproc != 0){
               Recv(&(cre_up_tmp[0]),(2*ncoef+1),MPI_DOUBLE,
                                          MPI_ANY_SOURCE,MPI_ANY_TAG,world);
             }/* endif */
             if(myid==0){
               for(inhc=1;inhc<=2*ncoef;inhc++){
                if(ibinary == 0){
                 fprintf(fp_dnamec,"%.8g\n",cre_up_tmp[inhc]);
                }else{
                 n = 1;
                 cre_dum = (float) cre_up_tmp[inhc];
                 fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
                }/*endif ibinary*/
               }/*endfor*/
             }/* endif */
           }/*endfor : is*/
         }/*endif : lsda*/
         Barrier(world);
       }/* endfor : iproc*/
     }/*endif:istate_nhc_opt*/


     }/*endfor:ip*/
   }/*endfor:ichain*/


  }/*endif:write extended*/

/* free locally assigned memory */
  if( (myid == 0) && (ibinary == 1)){
     cfree(c_array1);
     cfree(c_array2);  
  }/*endif*/

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/














