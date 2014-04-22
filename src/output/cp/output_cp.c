/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: output_cp                                    */
/*                                                                          */
/* This subprogram provides output for a CP on the                          */ 
/* ground state Born-Oppenheimer GGA-LDA surface                            */
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


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void output_cp(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,CP *cp)

/*==========================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  double etot,econv_now; 
  double vtot,avtot,pnow; 
  double apnow;
  double nhc_div;
  double a,b,c,tab,tac,tbc; 
  double deth;
  double c_nhc_div;
  double c_div = cp->cpcoeffs_info.cp_nfree;

/*         Local Pointers                             */
  int exit_flag         = general_data->timeinfo.exit_flag;
  int myid              = class->communicate.myid;
  int nproc             = class->communicate.np;
  MPI_Comm world        = class->communicate.world;
  int np_states         = cp->communicate.np_states;
  int igo;


/*=====================================================================*/
/* 0) Adjust form of CP vectors */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
   igo=1;
   if((general_data->timeinfo.itime % 
       general_data->filenames.iwrite_dump)==0 || exit_flag == 1){
    control_coef_transpose_bck(cp,3);
    igo=0;
   }/*endif*/

   if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_confc) == 0 && igo==1) {
    control_coef_transpose_bck(cp,1);
   }/*endif*/

 }/*endif*/

/*=====================================================================*/
/*  I)Open File option :                                               */

  if(myid == 0){   
    if(((general_data->timeinfo.itime)==0)&&((
         general_data->filenames.ifile_open))==1){
          initial_fopen_cp(class,general_data,bonded,cp);
    }/*endif*/
  }/* endif */
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  II) Write Initial Energies to screen                                 */

  if(myid == 0){   
    if(general_data->timeinfo.itime==0 && general_data->simopts.debug_cp == 0){
      initial_output_cp(class,general_data,bonded,cp);
    }/*endif*/
  }/* endif */
  if(nproc>1){Barrier(world);}

/*=======================================================================*/
/* II) Calculate some dinky quantities                                   */

  if(myid == 0){   
   if((general_data->timeinfo.itime)!=0){
     if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_screen) == 0 || 
        (general_data->timeinfo.itime % 
         general_data->filenames.iwrite_inst)==0 ||
         (exit_flag == 1)){
       get_cell(general_data->cell.hmat,&a,&b,&c,&tab,&tbc,&tac);
       dink_quant_calc_cp(class, general_data,cp,&etot,&econv_now,&vtot, 
                          &avtot,
                          &deth,&pnow, &apnow,&nhc_div,&c_nhc_div);
     } /*endif*/
  }/*endif*/
 }/* endif */
 if(nproc>1){Barrier(world);}

/*=======================================================================*/
/*  III) Write to the output to screen                                   */

  if(myid == 0){   
   if((general_data->timeinfo.itime) != 0){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_screen) == 0 ||
        (exit_flag == 1)){ 
     screen_write_cp(class,general_data,bonded,cp,
                     etot,econv_now,vtot,avtot,deth,pnow, apnow,
                     nhc_div,a,b,c,tab,tbc,tac,c_div,c_nhc_div);
    }/*endif*/
   }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* IV) Write to the output to dump file                                 */ 

  if((general_data->timeinfo.itime!=0)){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_dump)==0 ||
        (exit_flag == 1)){
      write_dump_file_cp(class,bonded,general_data,cp);
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* V) Write to the output to free energy                                */ 

  if(myid == 0){   
   if((general_data->timeinfo.itime)!=0&&(general_data->simopts.cp==1)){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_dump) == 0){
      write_free_energy_file(class,bonded,general_data);
    }/*endif*/
   }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* V) Write to the kseigs                                               */ 

  if((general_data->timeinfo.itime)!=0  && 
     (general_data->simopts.cp==1 || general_data->simopts.cp_wave==1 )){
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
   } /*endif*/
  }/* endif */
  if(nproc>1){Barrier(world);}

/*====================================================================*/
/* VII) Write to the config files                          */

  if((general_data->timeinfo.itime)!=0){
    write_config_files_cp(class,bonded,general_data,cp);
  }/*endif*/

/*======================================================================*/
/* VIII) Write to the inst avgs to inst file                            */

  if(myid == 0){   
   if((general_data->timeinfo.itime)!=0){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_inst) == 0 ){
      write_inst_file_cp(class,general_data,cp,etot,a,b,c,tac,tab,tbc);
    }/*endif*/
   }/*endif*/
  }/*endif*/
  if(nproc>1){Barrier(world);}

/*======================================================================*/
/* IX) Readjust form of CP vectors                                    */

 if(general_data->timeinfo.itime!=0 && np_states> 1){
   igo=1;
   if((general_data->timeinfo.itime % 
       general_data->filenames.iwrite_dump)==0 || exit_flag == 1){
    control_coef_transpose_fwd(cp,3);
    igo=0;
   }/*endif*/

   if((general_data->timeinfo.itime % 
         general_data->filenames.iwrite_confc) == 0 && igo==1) {
    control_coef_transpose_fwd(cp,1);
   }/*endif*/
 }/*endif*/

/*==========================================================================*/
     }/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_output_cp(CLASS *class,GENERAL_DATA *general_data,
                                              BONDED *bonded,CP *cp)

/*==========================================================================*/

{/* begin routine */

#include "../typ_defs/typ_mask.h"

  int pi_beads = class->clatoms_info.pi_beads;
  int ncons,np_tot,npairs,nsh_tot;
  double nhc_div,c_nhc_div;
  double c_div = cp->cpcoeffs_info.cp_nfree;
  double eu_conv = 1.0;
  if(general_data->filenames.iwrite_units==0){eu_conv = 1.0;}
  if(general_data->filenames.iwrite_units==1){eu_conv = KCAL;}
  if(general_data->filenames.iwrite_units==2){eu_conv = BOLTZ;}

/*==========================================================================*/
/* Output stuff */

  if((general_data->simopts.cp==1)||(general_data->simopts.debug_cp==1)){
    printf("Initial total      energy  %.10g\n",
                         (general_data->stat_avg.vintert +
                            + general_data->stat_avg.vintrat
                            + general_data->stat_avg.cp_ehart
                            + general_data->stat_avg.cp_exc
                            + general_data->stat_avg.cp_eext
                            + general_data->stat_avg.cp_eke
                            + general_data->stat_avg.cp_enl)*eu_conv);
    printf("Initial intra      energy  %g\n",
                general_data->stat_avg.vintrat*eu_conv);
    printf("Initial inter      energy  %.8g\n",
                general_data->stat_avg.vintert*eu_conv);
    printf("Initial kinetic    energy  %g\n",
                general_data->stat_avg.kinet*eu_conv);
    if(class->energy_ctrl.isep_vvdw == 1) {
      printf("Initial VDW        energy  %g\n",
             general_data->stat_avg.vvdw*eu_conv);
      printf("Initial Coul       energy  %g\n",
             general_data->stat_avg.vcoul*eu_conv);
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
    printf("Initial Torsion    energy  %g\n",
             general_data->stat_avg.vtorst*eu_conv);
    printf("Initial onefour    energy  %g\n",
             general_data->stat_avg.vonfot*eu_conv);
    printf("Initial surface    energy  %.10g\n",
             general_data->stat_avg.vsurft*eu_conv);
  }/*endif*/


    printf("Initial e-Energy   energy  %.10g\n",
                                (general_data->stat_avg.cp_ehart
                            + general_data->stat_avg.cp_exc
                            + general_data->stat_avg.cp_eext
                            + general_data->stat_avg.cp_eke
                            + general_data->stat_avg.cp_enl)*eu_conv);
    printf("Initial e-Hart+XC  energy  %.10g\n",
                            (general_data->stat_avg.cp_ehart
                            + general_data->stat_avg.cp_exc)*eu_conv);
    printf("Initial e-Hart+muXC-XC  energy  %.10g\n",
                                      (  general_data->stat_avg.cp_ehart
                                       + general_data->stat_avg.cp_muxc
                                 - general_data->stat_avg.cp_exc)*eu_conv);
    printf("Initial e-External energy  %.10g\n",
                                general_data->stat_avg.cp_eext*eu_conv);
    printf("Initial e-Nonlocal energy  %.10g\n",
                                general_data->stat_avg.cp_enl*eu_conv);
    printf("Initial e-Kinetic  energy  %.10g\n",
                                general_data->stat_avg.cp_eke*eu_conv);
    printf("Initial volume             %g\n",
                                   general_data->stat_avg.vol*BOHR*BOHR*BOHR);
   if(general_data->simopts.cp==1){
     printf("Initial Pressure      %g\n",
        (((general_data->ptens.tvten)[1]+(general_data->ptens.pvten_tot)[1]
      +(general_data->ptens.tvten)[5]+(general_data->ptens.pvten_tot)[5]
      +(general_data->ptens.tvten)[9]+(general_data->ptens.pvten_tot)[9])
                        /(3.0*(general_data->stat_avg.vol)*PCONV)));
     printf("Initial Atm temperature    %g\n",
            ((2.0*general_data->stat_avg.kinet*BOLTZ)/
           ((double)(class->clatoms_info.nfree))));
    }/*endif*/
    printf("Initial e-Temperature      %g\n",general_data->stat_avg.kinet_cp
                                             *2.0*BOLTZ/c_div);
/*==========================================================================*/
  /* NHC quantities                                                        */

  if(general_data->simopts.cp==1){
    if(general_data->ensopts.nvt==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc) );
      printf("Initial particle NHC temperature %g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      if(pi_beads>1){
        nhc_div = (double) ( class->therm_info_bead.len_nhc*
                            (class->therm_info_bead.num_nhc) );
        printf("Initial bead NHC temperature %g\n", 
              ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      }/*endif*/
    }/*endif*/

    if(general_data->ensopts.npt_i==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature %g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      nhc_div = (double) ( class->therm_info_bead.len_nhc*
                          (class->therm_info_bead.num_nhc +1));
      if(pi_beads>1){
       printf("Initial bead NHC temperature %g\n", 
             ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
       printf("Initial Volume temperature %g\n",
             (2.0*general_data->stat_avg.kinet_v*BOLTZ)); 
     }/*endif*/
    }/*endif*/

    if(general_data->ensopts.npt_f==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature %g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      nhc_div = (double) ( class->therm_info_bead.len_nhc*
                          (class->therm_info_bead.num_nhc +1));
      printf("Initial bead NHC temperature %g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      printf("Initial Volume temperature %g\n",      
            ((2.0*general_data->stat_avg.kinet_v*BOLTZ)/9.0));
    }/*endif*/
  }/*endif*/

   if(cp->cptherm_info.num_c_nhc != 0) {
    c_nhc_div = (double) ((cp->cptherm_info.len_c_nhc)*
                          (cp->cptherm_info.num_c_nhc));
    printf("CP NHC Temperature  = %g\n",general_data->stat_avg.kinet_nhc_cp
                                 *BOLTZ*2.0/c_nhc_div);
    } /* endif */

/*==========================================================================*/
  /* Neighbor list quantities                                              */

  if(general_data->simopts.cp==1){
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
      nsh_tot = ((class->nbr_list.lnklist.ncell_a)*
                 (class->nbr_list.lnklist.ncell_b)*
               (class->nbr_list.lnklist.ncell_c)-1)/2 +
                 (((class->nbr_list.lnklist.ncell_a)*
                   (class->nbr_list.lnklist.ncell_b)*
                   (class->nbr_list.lnklist.ncell_c)-1) % 2) + 1;
      printf("ncell_a = %d, ncell_b = %d, ncell_c = %d natm_cell_max %d \n",
            class->nbr_list.lnklist.ncell_a,class->nbr_list.lnklist.ncell_b,
            class->nbr_list.lnklist.ncell_c,
            class->nbr_list.lnklist.natm_cell_max);
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

void initial_fopen_cp(CLASS *class,GENERAL_DATA *general_data,
                                    BONDED *bonded,CP *cp)

/*==========================================================================*/

{/* begin routine */
/*==========================================================================*/

#include "../typ_defs/typ_mask.h"

  int n=1,iii,ibinary,iwrite_now;
  NAME file_typ;
  FILE *fp_bond_free,*fp_bend_free,*fp_tors_free,*fp_rbar_free;
  FILE *fp_iname, *fp_cpname, *fp_cvname,*fp_dname,*fp_ccname,*fp_kseigs;
  FILE *fp_cfname;

/*==========================================================================*/
/*     A) Open dump file                                          */

  fp_dname = cfopen(general_data->filenames.dname,"w");
  fclose(fp_dname);

/*==========================================================================*/
/*     A) Open inst avg file                                          */

 if(general_data->simopts.cp==1){ 
  ibinary    = 0; 
  iwrite_now = general_data->filenames.iwrite_inst;
  strcpy(file_typ,"ins_file");
  fp_iname   = cfopen(general_data->filenames.iname,"w");
  write_gen_header_cp(class,general_data,cp,fp_iname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_iname);
 }/*endif*/

/*==========================================================================*/
/*     B) Open atm vel conf file                                       */

  if(general_data->simopts.cp==1){ 
    ibinary    = general_data->filenames.iwrite_conf_binary;
    iwrite_now = general_data->filenames.iwrite_confv;
    strcpy(file_typ,"vel_file");
    fp_cvname  = cfopen(general_data->filenames.cvname,"w");
    write_gen_header_cp(class,general_data,cp,fp_cvname,ibinary,
                     iwrite_now,file_typ);
    fclose(fp_cvname); 
  }/*endif*/

/*==========================================================================*/
/*     C) Open force conf file                                           */

    ibinary    = general_data->filenames.iwrite_conf_binary;
    iwrite_now = general_data->filenames.iwrite_atm_for;
    strcpy(file_typ,"force_file");
    fp_cfname  = cfopen(general_data->filenames.forcename,"w");
    write_gen_header(class,general_data,fp_cfname,ibinary,
                     iwrite_now,file_typ);
    fclose(fp_cfname);

/*==========================================================================*/
/*     D) Open pos conf file                                           */

    ibinary    = general_data->filenames.iwrite_conf_binary;
    iwrite_now = general_data->filenames.iwrite_confp;
    strcpy(file_typ,"pos_file");
    fp_cpname  = cfopen(general_data->filenames.cpname,"w");
    write_gen_header_cp(class,general_data,cp,fp_cpname,ibinary,
                     iwrite_now,file_typ);
    fclose(fp_cpname); 

/*======================================================================*/
/*     E) Open partial pos conf file                                    */

 if(general_data->simopts.cp==1){ 
  if((general_data->filenames.low_lim_par<=
      general_data->filenames.high_lim_par)){
   ibinary    = general_data->filenames.iwrite_conf_binary;
   iwrite_now = general_data->filenames.iwrite_par_confp;
   strcpy(file_typ,"par_file");
   fp_cpname = cfopen(general_data->filenames.cpparname,"w"); 
   write_gen_header_cp(class,general_data,cp,fp_cpname,ibinary,
                   iwrite_now,file_typ);
   fclose(fp_cpname); 
  }/*endif*/
 }/*endif*/

/*==========================================================================*/
/*     F) Open bond free energy file                                   */

 if(general_data->simopts.cp==1){ 
    if(bonded->bond_free.num>0){
      fp_bond_free  = cfopen(bonded->bond_free.file,"w");
      fclose(fp_bond_free);
    }/*endif*/
 }/*endif*/

/*==========================================================================*/
/*     G) Open bend free energy file                                    */

 if(general_data->simopts.cp==1){ 
    if(bonded->bend_free.num>0){
      fp_bend_free  = cfopen(bonded->bend_free.file,"w");
      fclose(fp_bend_free);
    }/*endif*/
 }/*endif*/

/*==========================================================================*/
/*     H) Open tors free energy file                                    */

 if(general_data->simopts.cp==1){ 
    if(bonded->tors_free.num>0){
      fp_tors_free  = cfopen(bonded->tors_free.file,"w");
      fclose(fp_tors_free);
    }/*endif*/
 }/*endif*/    

/*==========================================================================*/
/*     I) Open rbar_sig free energy file                                    */

 if(general_data->simopts.cp==1){ 
    if(bonded->rbar_sig_free.nfree>0){
      fp_rbar_free  = cfopen(bonded->rbar_sig_free.file,"w");
      fclose(fp_rbar_free);
    }/*endif*/
 }/*endif*/        

/*==========================================================================*/
/*     J) Open CP conf file (Must be unformatted due to hugeness of         */
/*          wave functions)                                                 */

 if(general_data->simopts.cp == 1){
   ibinary    = 1;
   iwrite_now = general_data->filenames.iwrite_confc;
   strcpy(file_typ,"cof_file");
   fp_ccname = cfopen(general_data->filenames.ccname,"w");
   write_gen_header_cp(class,general_data,cp,fp_ccname,ibinary,
                   iwrite_now,file_typ);
   fclose(fp_ccname); 
 }
    
/*==========================================================================*/
/*  K) Open CP KS eigenvalues file                                          */

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

void dink_quant_calc_cp(CLASS *class, GENERAL_DATA *general_data,
                        CP *cp,double *etot,
                        double *econv_now,double *vtot, double *avtot,
                        double *deth,double *pnow, double *apnow,
                        double *nhc_div,double *c_nhc_div)

/*==========================================================================*/

{/* begin routine */

  int i,ncons;

#include "../typ_defs/typ_mask.h"

/*==========================================================================*/
 /*  Energy */

  (*etot) = general_data->stat_avg.kinet_cp + general_data->stat_avg.cp_ehart
          + general_data->stat_avg.cp_eext + general_data->stat_avg.cp_exc  
          + general_data->stat_avg.cp_eke + general_data->stat_avg.cp_enl;
  if(cp->cptherm_info.num_c_nhc > 0){
       (*etot) += ( general_data->stat_avg.kinet_nhc_cp 
                   +general_data->stat_avg.vpotnhc_cp );
  }/*endif:nhc_cp*/
  if(general_data->simopts.cp == 1) {
        (*etot) += ( (general_data->stat_avg.kinet) 
                    +(general_data->stat_avg.vintert) 
                    +(general_data->stat_avg.vintrat) );
        if((general_data->ensopts.nvt)==1)  
         {(*etot)     += ( (general_data->stat_avg.kinet_nhc) 
                          +(general_data->stat_avg.vpotnhc) );}
        if((general_data->ensopts.npt_i)==1)
         {(*etot)     += ( (general_data->stat_avg.kinet_nhc) 
                          +(general_data->stat_avg.vpotnhc)
                          +(general_data->stat_avg.kinet_v)   
                          +(general_data->stat_avg.vpot_v) );}
       if((general_data->ensopts.npt_f)==1)
         {(*etot)+= ( (general_data->stat_avg.kinet_nhc) 
                     +(general_data->stat_avg.vpotnhc)
                     +(general_data->stat_avg.kinet_v)   
                     +(general_data->stat_avg.vpot_v) );}

       if(((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==1)&&
           (general_data->stat_avg.iswit_vdw<=0))
         {(*etot)+= (general_data->stat_avg.vlong);}
  }/* endif  cp on*/

  (*econv_now) = fabs((*etot)-general_data->stat_avg.econv0)/
                 fabs(general_data->stat_avg.econv0);

 /*==========================================================================*/
 /*  Volume and pressure                                                    */

  if(general_data->simopts.cp == 1) {

    (*vtot) = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;
    (*avtot) = general_data->stat_avg.avintert+general_data->stat_avg.avintrat;
    (*deth) = getdeth(general_data->cell.hmat);
    for(i=1;i<=9;i++){
     general_data->stat_avg.apten_out[i]=(general_data->ptens.pvten_tot[i]
                            +general_data->ptens.tvten[i])/((*deth)*PCONV);
    }
    (*pnow)  = (general_data->stat_avg.apten_out[1]
            +  general_data->stat_avg.apten_out[5]
            +  general_data->stat_avg.apten_out[9])/(3.0);
    (*apnow) =  general_data->stat_avg.apress
               /(PCONV*((double)(general_data->timeinfo.itime)));

  }/*endif*/

 /*==========================================================================*/
 /*  NHC degrees of freedom                                                  */

  if(general_data->simopts.cp == 1) {
      (*nhc_div) = (double)((class->therm_info_class.len_nhc)*
                            (class->therm_info_class.num_nhc));
      if(general_data->ensopts.npt_f == 1 || general_data->ensopts.npt_i == 1){
       (*nhc_div) = (double) (class->therm_info_class.len_nhc*
                           (class->therm_info_class.num_nhc +1) );
      }/*endif*/
  }/*endif*/

/*==========================================================================*/
/*  Calculate c_nhc_div                                                     */


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

void screen_write_cp(CLASS *class,GENERAL_DATA *general_data,
                     BONDED *bonded,CP *cp,
                     double etot,double econv_now,double vtot, 
                     double avtot,double deth,double pnow, double apnow,
                     double nhc_div,double a,double b,double c,
                     double tab,double tbc,double tac,double c_div,
                     double c_nhc_div)

/*==========================================================================*/
{/* begin routine */
#include "../typ_defs/typ_mask.h"

 int npairs;
 int cp_lsda = cp->cpopts.cp_lsda;
 int iii;
 int anneal_opt = general_data->simopts.anneal_opt;
 static int count_norb=0;
 double atime,vol_div,atm_div;
 double updates_t,updates_true,updates_now; 
 double eelec;
 double aelec;
 double time_fact = general_data->timeinfo.dt * TIME_CONV / 1000.0;

/*==========================================================================*/
/* Write to screen                                                          */

   atime = (double)(general_data->timeinfo.itime);
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

  #define FMT "%18.10f" 

  printf("\n");
  printf("****************************************************************************\n");
  printf("\n");
  if(general_data->ensopts.nve==1)  printf("Ensemble          = NVE\n");
  if(general_data->ensopts.nvt==1)  printf("Ensemble          = NVT\n");
  if(general_data->ensopts.npt_i==1)printf("Ensemble          = NPT-ISO\n");
  if(general_data->ensopts.npt_f==1)printf("Ensemble          = NPT-FLEX\n");
  printf("Time step         = %d\n", general_data->timeinfo.itime);
  printf("Time              = %.4f ps\n", atime * time_fact);
  printf("\n");
  printf("QUANTITY                 INSTANTANEOUS            AVERAGE\n");
  printf("---------------------------------------------------------\n");
  printf("\n");
     printf("Econv             = "FMT" "FMT"\n",( econv_now ),
	 ( general_data->stat_avg.econv/atime));
  if(cp->cpopts.cp_isok_opt == 1)
     printf("CP KE conv        = "FMT" "FMT"\n",
	    fabs((general_data->stat_avg.kinet_cp - general_data->stat_avg.cp_kconv0)/
		 general_data->stat_avg.cp_kconv0),general_data->stat_avg.cp_kconv/atime);
   if(general_data->simopts.cp==1){
      printf("Total Energy      = "FMT" "FMT"\n",(general_data->stat_avg.kinet+vtot+eelec),
            (general_data->stat_avg.akinet+avtot+aelec)/atime);
      printf("Atm Energy        = "FMT" "FMT"\n",(general_data->stat_avg.kinet+vtot),
            (general_data->stat_avg.akinet+avtot)/atime);
      printf("Total Atm PE      = "FMT" "FMT"\n",(vtot),(avtot/atime));
      printf("Intermol Atm PE   = "FMT" "FMT"\n",(general_data->stat_avg.vintert),
            (general_data->stat_avg.avintert/atime));
      printf("Intramol Atm PE   = "FMT" "FMT"\n",(general_data->stat_avg.vintrat),
            (general_data->stat_avg.avintrat/atime));
      printf("Atm KE            = "FMT" "FMT"\n",(general_data->stat_avg.kinet),
            (general_data->stat_avg.akinet/atime));
      printf("\n");
      atm_div = (double)(class->clatoms_info.nfree);
      printf("Atm Deg. Free     = "FMT"\n",atm_div); 
      printf("Atm Temperature   = "FMT" "FMT"\n",
            (general_data->stat_avg.kinet*2.0*BOLTZ/atm_div),
            (general_data->stat_avg.akinet*2.0*BOLTZ/(atm_div*atime)));
     
   }/* endif */
   printf("\n");
   printf("e-Energy          = "FMT" "FMT"\n",
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
   printf("e-Hartree + XC    = "FMT" "FMT"\n",(general_data->stat_avg.cp_ehart
                                      + general_data->stat_avg.cp_exc),
                                           (general_data->stat_avg.acp_ehart
                                     + general_data->stat_avg.acp_exc)
                                            /atime);
   printf("e-External PE     = "FMT" "FMT"\n",(general_data->stat_avg.cp_eext),
                                         (general_data->stat_avg.acp_eext)
                                            /atime);
   printf("e-Nonlocal PE     = "FMT" "FMT"\n",(general_data->stat_avg.cp_enl),
                                          (general_data->stat_avg.acp_enl)
                                            /atime);
   printf("e-Kinetic         = "FMT" "FMT"\n",(general_data->stat_avg.cp_eke),
                                          (general_data->stat_avg.acp_eke)
                                            /atime);
   printf("\n");
   printf("CP Fict KE        = "FMT" "FMT"\n",(general_data->stat_avg.kinet_cp),
                                           (general_data->stat_avg.akinet_cp/atime));

   printf("CP Temperature    = "FMT" "FMT"\n",(general_data->stat_avg.kinet_cp
                                          *2.0*BOLTZ/c_div),
                                           (general_data->stat_avg.akinet_cp
                                          *2.0*BOLTZ/(c_div*atime)));

   if(anneal_opt == 0){
    if(general_data->stat_avg.kinet_cp*2.0*BOLTZ/c_div>100.0*cp->cpopts.te_ext){
     printf("\n$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
     printf("Your current CP temperature is greater than\n");
     printf("100 x the value set in the input file.\n");
     printf("Please reconsider your adiabaticity parameters\n");
     printf("and/or employ coefficient thermostats to control \n");
     printf("this quantity.\n");
     printf("$$$$$$$$$$$$$$$$-WARNING-$$$$$$$$$$$$$$$$$$$\n");
     fflush(stdout);
    }/*endif*/
   }/* endif */

   printf("CP Up Force:   avg: "FMT" max: "FMT"\n",general_data->stat_avg.fc_mag_up,
                                             general_data->stat_avg.fc_max_up);
   if(cp_lsda==1){
    printf("CP Dn Force:   avg: "FMT" max: "FMT"\n",general_data->stat_avg.fc_mag_dn,
                                             general_data->stat_avg.fc_max_dn);
   }/*endif*/

   printf("\n");

/*==========================================================================*/
/*     B.1) Extended Class                                          */

  if(general_data->simopts.cp==1){
   if((general_data->ensopts.nvt + general_data->ensopts.npt_i
      + general_data->ensopts.npt_f == 1)){
     printf("NHC Temperature   = "FMT" "FMT"\n",
           (general_data->stat_avg.kinet_nhc*2.0*BOLTZ/nhc_div),
            (general_data->stat_avg.akinet_nhc*2.0*BOLTZ
            /(nhc_div*atime)));
   }/*endif*/
   if((general_data->ensopts.npt_i +general_data->ensopts.npt_f == 1)) {
     vol_div = 1.0;
     if(general_data->ensopts.npt_f == 1) vol_div = 6.0;
     printf("Vol Temperature   = "FMT" "FMT"\n",
              (general_data->stat_avg.kinet_v*2.0*BOLTZ/vol_div),
              (general_data->stat_avg.akinet_v*2.0*BOLTZ/(atime*vol_div)));
   }/*endif*/
  }/* endif cp==1 */

/*==========================================================================*/
/*     B.2) Extended Coefs                                                  */

  if(cp->cptherm_info.num_c_nhc != 0) {
      printf("CP NHC Temp       = "FMT" "FMT"\n",
          (general_data->stat_avg.kinet_nhc_cp*2.0*BOLTZ)/c_nhc_div,
          (general_data->stat_avg.akinet_nhc_cp*2.0*BOLTZ)/(c_nhc_div*atime));
  }/* endif */

/*==========================================================================*/
/*     C)Pressure/Vol                                                      */

   if(general_data->simopts.cp == 1&&general_data->cell.iperd>=2
      && cp->cpopts.cp_ptens_calc == 1) {
       printf("Pressure          = "FMT" "FMT"\n",(pnow),( apnow));
       printf("Volume            = "FMT" "FMT"\n",(deth*BOHR*BOHR*BOHR),
              (general_data->stat_avg.avol/atime)*BOHR*BOHR*BOHR);
   }/* endif */
   printf("\n");

/*==========================================================================*/
/*     D.1)Constraint CLassical                                            */

   if(bonded->constrnt.iconstrnt == 1) {
       if(general_data->simopts.cp==1){
        if(bonded->bond.ncon > 0) {
         printf("Shake iter        = "FMT" "FMT"\n",
               (double)(general_data->stat_avg.iter_shake),
               (general_data->stat_avg.aiter_shake/atime));
         printf("Rattle iter       = "FMT" "FMT"\n",
               (double)(general_data->stat_avg.iter_ratl),
               (general_data->stat_avg.aiter_ratl/atime));
        }
        if(bonded->grp_bond_con.num_21 > 0) {
        printf("Grp_21 shake iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_21,
               general_data->stat_avg.aiter_21/atime);
         printf("Grp_21 ratl iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_21r,
               general_data->stat_avg.aiter_21r/atime);
        }
        if(bonded->grp_bond_con.num_23 > 0) {
        printf("Grp_23 shake iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_23,
               general_data->stat_avg.aiter_23/atime);
         printf("Grp_23 ratl iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_23r,
               general_data->stat_avg.aiter_23r/atime);
        }
        if(bonded->grp_bond_con.num_33 > 0) {
         printf("Grp_33 shake iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_33,
               general_data->stat_avg.aiter_33/atime);
         printf("Grp_33 ratl iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_33r,
               general_data->stat_avg.aiter_33r/atime);
        }
        if(bonded->grp_bond_con.num_46 > 0) {
         printf("Grp_46 shake iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_46,
               general_data->stat_avg.aiter_46/atime);
         printf("Grp_46 ratl iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_46r,
               general_data->stat_avg.aiter_46r/atime);
        }
        if(bonded->grp_bond_con.num_43 > 0) {
         printf("Grp_43 shake iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_43,
               general_data->stat_avg.aiter_43/atime);
         printf("Grp_43 ratl iter = "FMT" "FMT"\n",
               general_data->stat_avg.iter_43r,
               general_data->stat_avg.aiter_43r/atime);
        }
         printf("\n");
       }/*endif*/
   }/*endif*/

/*==========================================================================*/
/*     D.2)Constraint Coef                                                 */

   if(!cp->cpopts.cp_norb) {
       printf("Orb shake iter    = "FMT" "FMT"\n",
                        (double)(general_data->stat_avg.iter_shake_cp),
                                (general_data->stat_avg.aiter_shake_cp/atime));
      printf("Orb rattle iter   = "FMT" "FMT"\n",
                        (double) (general_data->stat_avg.iter_ratl_cp),
                                 (general_data->stat_avg.aiter_ratl_cp/atime));
   }/*endif*/
   if(cp->cpopts.cp_norb > 0) {
      if(cp->cpopts.cp_norb == 2 || cp->cpopts.cp_norb == 3){
          printf("Max off diag      = "FMT"\n",
                                  cp->cpcoeffs_info.max_off_diag);
      }/*endif*/
      if(cp->cpopts.cp_norb == 3){
          printf("Max diag          = "FMT"\n",
                                  cp->cpcoeffs_info.max_diag);
      }/*endif*/
      printf("Rotations done    = "FMT"\n",
                             general_data->stat_avg.count_diag_srot);
   }/*endif*/

/*==========================================================================*/
   /* timing */

   printf("Cpu time          = "FMT" "FMT"\n",
         general_data->stat_avg.cpu_now,
         general_data->stat_avg.acpu/atime);
   printf("\n");

/*==========================================================================*/
/*     D)Misc                                                         */

   if(class->nbr_list.iver == 1 && general_data->simopts.cp == 1){
      updates_t = general_data->stat_avg.updates;
      if(updates_t == 0) updates_t = 1;
      updates_now = (double ) (general_data->timeinfo.itime-
                               general_data->stat_avg.itime_update);
       updates_true = updates_t;
       npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
            + class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
       printf("Inst steps/update = "FMT"\n",updates_now);
       printf("Avg. steps/update = "FMT"\n",(atime/updates_t));
       printf("Total list updates= "FMT"\n",updates_true);
       printf("Number of pairs   = %d \n",npairs);
       printf("\n");
   }/*endif*/

  /*=====================*/
  /* pressure/vol matrix */

  if (general_data->cell.iperd >= 2) {
    printf("QUANTITY\n");
    printf("----------------------------------------------------------------------------\n");
    printf("\n");
       printf("Inst P11,P22,P33  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten_out[1]),
              (general_data->stat_avg.apten_out[5]),
              (general_data->stat_avg.apten_out[9]));
       printf("Avg  P11,P22,P33  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten[1]/(PCONV*atime)),
              (general_data->stat_avg.apten[5]/(PCONV*atime)),
              (general_data->stat_avg.apten[9]/(PCONV*atime)));
       printf("Inst P12,P13,P23  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten_out[4]),
              (general_data->stat_avg.apten_out[7]),
              (general_data->stat_avg.apten_out[8]));
       printf("Avg  P12,P13,P23  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten[4])/(PCONV*atime),
              (general_data->stat_avg.apten[7])/(PCONV*atime),
              (general_data->stat_avg.apten[8])/(PCONV*atime));
       printf("Inst P21,P31,P32  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten_out[2]),
              (general_data->stat_avg.apten_out[3]),
              (general_data->stat_avg.apten_out[6]));
       printf("Avg  P21,P31,P32  = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.apten[2])/(PCONV*atime),
              (general_data->stat_avg.apten[3])/(PCONV*atime),
              (general_data->stat_avg.apten[6])/(PCONV*atime));
       printf("\n");
       printf("Inst cell lths    = "FMT" "FMT" "FMT"\n",(a*BOHR),(b*BOHR),(c*BOHR));
       printf("Avg  cell lths    = "FMT" "FMT" "FMT"\n",    
              (general_data->stat_avg.acella/atime)*BOHR,
              (general_data->stat_avg.acellb/atime)*BOHR,
              (general_data->stat_avg.acellc/atime)*BOHR);
       printf("Inst cell angs    = "FMT" "FMT" "FMT"\n",(tab),(tac),(tbc));
       printf("Avg  cell angs    = "FMT" "FMT" "FMT"\n",
              (general_data->stat_avg.acellab/atime),
              (general_data->stat_avg.acellac/atime),
              (general_data->stat_avg.acellbc/atime));
       printf("\n");
}

   printf("----------------------------------------------------------------------------\n");
   printf("\n"); 
   printf("****************************************************************************\n");
   printf("\n"); 
   fflush(stdout);

/*==========================================================================*/
    }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_file_cp(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                        CP *cp)

/*==========================================================================*/
{/* begin routine */

#include "../typ_defs/typ_mask.h"

 int i,j,ncoef_up,ncoef_dn,izero;
 int pi_beads = class->clatoms_info.pi_beads;
 int iproc,ioff_re,ioff_im;
 double deth;
 FILE *fp_dname,*fp_dnamec;

 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo,iii;

 int ibinary = cp->cpopts.iwrite_coef_binary;
 int cp_dual_grid_opt = cp->cpopts.cp_dual_grid_opt;
 double sx,sy,sz;
 double cp_box_center_tmp[4];
 double *cp_box_center = general_data->cell.cp_box_center;
 double *hmati         = general_data->cell.hmati;

 MPI_Status stat;

/* Local pointers */
 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 MPI_Comm world = class->communicate.world;
 int nproc      = class->communicate.np;

/*==========================================================================*/

/*==========================================================================*/
/*==========================================================================*/
/* Open dump file                                                           */

 if(myid==0){

  fp_dname = cfopen(general_data->filenames.dname,"o");

/*==========================================================================*/
/*     A)Atm positions                                               */

      fprintf(fp_dname,"natm_tot restart_typ itime\n");
      fprintf(fp_dname,"%d restart_all %d %d\n",class->clatoms_info.natm_tot,
             general_data->timeinfo.itime,pi_beads);
      fprintf(fp_dname,"atm pos, atm_typ, mol_typ mol_num\n");
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
       fprintf(fp_dname,"%.12g %.12g %.12g %s %s %s %d\n",
              class->clatoms_pos[1].x[i],class->clatoms_pos[1].y[i],
                                         class->clatoms_pos[1].z[i],
              class->atommaps.atm_typ[class->atommaps.iatm_atm_typ[i]],
              class->atommaps.res_typ[class->atommaps.iatm_res_typ[i]],
              class->atommaps.mol_typ[class->atommaps.iatm_mol_typ[i]],
              class->atommaps.iatm_mol_num[i]);
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
    if( cp_dual_grid_opt >= 1){

      fprintf(fp_dname,"h matrix cp\n");
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_cp[1],
             general_data->cell.hmat_cp[4],general_data->cell.hmat_cp[7]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_cp[2],
             general_data->cell.hmat_cp[5],general_data->cell.hmat_cp[8]);
      fprintf(fp_dname,"%.12g %.12g %.12g\n",general_data->cell.hmat_cp[3],
             general_data->cell.hmat_cp[6],general_data->cell.hmat_cp[9]);

/* convert cp_box_center back into xtal coordinates*/
      fprintf(fp_dname,"cp box center\n");

    sx = cp_box_center[1];
    sy = cp_box_center[2];
    sz = cp_box_center[3];

    cp_box_center_tmp[1] =sx*hmati[1]+sy*hmati[4]+sz*hmati[7];
    cp_box_center_tmp[2] =sx*hmati[2]+sy*hmati[5]+sz*hmati[8];
    cp_box_center_tmp[3] =sx*hmati[3]+sy*hmati[6]+sz*hmati[9];

      fprintf(fp_dname,"%.12g %.12g %.12g\n",cp_box_center_tmp[1],cp_box_center_tmp[2],
                                    cp_box_center_tmp[3]);

    }/*endif cp_dual_grid_opt */

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
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
       fprintf(fp_dname,"%.12g %.12g %.12g\n",class->clatoms_pos[1].vx[i],
              class->clatoms_pos[1].vy[i],class->clatoms_pos[1].vz[i]);
       }/*endfor*/
      
      fprintf(fp_dname,"number of atm nhc, length of nhc\n");
      fprintf(fp_dname,"%d %d\n",class->therm_info_class.num_nhc,
             class->therm_info_class.len_nhc);
      fprintf(fp_dname,"atm nhc velocities\n");
      for(j=1;j<=(class->therm_info_class.len_nhc);j++){
       for(i=1;i<=(class->therm_info_class.num_nhc);i++){
         fprintf(fp_dname,"%.12g\n",class->therm_class.v_nhc[j][i]);
       } /*endfor*/
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
                      general_data->simopts.cp,general_data->simopts.cp_wave,  
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

 }/*endif : myid==0*/
 if(nproc>1){Barrier(world);}

/*==========================================================================*/
/*==========================================================================*/
/*==========================================================================*/
/* Open cp dump file                                                        */


if(myid ==0){
 if(ibinary == 0){
  fp_dnamec = fopen(general_data->filenames.dnamec,"w");
 }else{
  fp_dnamec = fopen(general_data->filenames.dnamec,"wb");
 }/*endif*/
}

/*==========================================================================*/
/*==========================================================================*/
/* Write the header                                                         */

/*==========================================================================*/
/* Write the header                                                         */

  write_dump_header_cp(fp_dnamec,cp,general_data,myid,ibinary);

  if(nproc>1){Barrier(world);}

/*==========================================================================*/
/* write the occupation number                                              */

  write_dump_occ_cp(fp_dnamec,cp,myid,ibinary);

  if(nproc>1){Barrier(world);}

/*==========================================================================*/
/*==========================================================================*/
/* Write the coefficients                                                  */

  write_dump_coef_cp(fp_dnamec,cp,class,ibinary);
 
/*==========================================================================*/
/* Write the coefficient velocities                                         */

  write_dump_vcoef_cp(fp_dnamec,cp,class,general_data,ibinary);

/*==========================================================================*/
/* Write the extended class stuff                                         */

 write_dump_extended_cp(fp_dnamec,cp,class,general_data,ibinary);

/*==========================================================================*/
/* Done                                                                     */

   if(myid==0) {
    fflush(fp_dnamec); 
    fclose(fp_dnamec); 
   };

/*==========================================================================*/
    }/* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_config_files_cp(CLASS *class,BONDED *bonded,
                           GENERAL_DATA *general_data, CP *cp)

/*==========================================================================*/

{/* begin routine */

#include "../typ_defs/typ_mask.h"

  int i,n,ncoef_up,ncoef_dn;
  int is,igo,iproc;

  int cp_on            = general_data->simopts.cp;

  int cp_lsda = cp->cpopts.cp_lsda;
  FILE *fp_cpname, *fp_cvname,*fp_ccname, *fp_cfname;
  int  ncoef  = cp->cpcoeffs_info.ncoef;

 int  istate_up_st  = cp->cpcoeffs_info.istate_up_st;
 int  istate_up_end = cp->cpcoeffs_info.istate_up_end;
 int  nstate_up     = cp->cpcoeffs_info.nstate_up;
 int  nstate_up_proc= cp->cpcoeffs_info.nstate_up_proc;

 int  istate_dn_st  = cp->cpcoeffs_info.istate_dn_st;
 int  istate_dn_end = cp->cpcoeffs_info.istate_dn_end;
 int  nstate_dn     = cp->cpcoeffs_info.nstate_dn;
 int  nstate_dn_proc= cp->cpcoeffs_info.nstate_dn_proc;
 int isoff;
 
 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;
 int myid = class->communicate.myid;
 MPI_Comm world = class->communicate.world;
 int nproc      = class->communicate.np;

/*=====================================================================*/
/* 0) */
    ncoef_up = cp->cpcoeffs_info.ncoef;
    ncoef_dn = cp->cpcoeffs_info.ncoef;
    if(cp->cpopts.cp_lda==1){ncoef_dn=0;}

/*=====================================================================*/
/* I) Write to the atm velocity config file                          */

  if(myid==0 ){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_confp) == 0 ){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cpname = cfopen(general_data->filenames.cpname,"a");

       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
       fprintf(fp_cpname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[1].x[i],
              class->clatoms_pos[1].y[i],class->clatoms_pos[1].z[i]);
       }/*endfor*/
      for(i=0;i<3;i++) 
       fprintf(fp_cpname,"%.12g %.12g %.12g\n",general_data->cell.hmat[1+i],
              general_data->cell.hmat[4+i],general_data->cell.hmat[7+i]);
      fflush(fp_cpname);
      fclose(fp_cpname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cpname = cfopen(general_data->filenames.cpname,"a");
    n=1;
    for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
     fwrite(&(class->clatoms_pos[1].x)[i],sizeof(double),n,fp_cpname);
     fwrite(&(class->clatoms_pos[1].y)[i],sizeof(double),n,fp_cpname);
     fwrite(&(class->clatoms_pos[1].z)[i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
    for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
      fflush(fp_cpname);
      fclose(fp_cpname);
   }/*endif*/
   }
  }/* endif myid==0 */
  if(nproc>1){Barrier(world); }

/*=====================================================================*/
 /* II) Write to the atm velocity config file                          */

  if(myid==0 && cp_on == 1){
    if((general_data->timeinfo.itime % 
        general_data->filenames.iwrite_confv) == 0 
        &&(general_data->simopts.cp==1)){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cvname = cfopen(general_data->filenames.cvname,"a");
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
       fprintf(fp_cvname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[1].vx[i],
              class->clatoms_pos[1].vy[i],class->clatoms_pos[1].vz[i]);
       }/*endfor*/
      for(i=0;i<3;i++) 
       fprintf(fp_cvname,"%.12g %.12g %.12g\n",general_data->cell.hmat[1+i],
              general_data->cell.hmat[4+i],general_data->cell.hmat[7+i]);
      fflush(fp_cvname);
      fclose(fp_cvname);
    }/*endif*/
   }/*endif binary*/

   if(general_data->filenames.iwrite_conf_binary==1){
     fp_cvname = cfopen(general_data->filenames.cvname,"a");
     n=1;
      for(i=1;i<=(class->clatoms_info.natm_tot);i++){ 
       fwrite(&(class->clatoms_pos[1].vx)[i],sizeof(double),n,fp_cvname);
       fwrite(&(class->clatoms_pos[1].vy)[i],sizeof(double),n,fp_cvname);
       fwrite(&(class->clatoms_pos[1].vz)[i],sizeof(double),n,fp_cvname);
      }/*endfor*/ 
     for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cvname);
     }/*endfor*/ 
     fflush(fp_cvname);
     fclose(fp_cvname);
   }/*endif:binary*/


  }/* endif myid */
  if(nproc>1){Barrier(world);}

/*=====================================================================*/
 /* III) Write to the coefficient file                                 */
/*       Unformatted output is absolutely necessary here!!!!           */

     if((general_data->timeinfo.itime%
         general_data->filenames.iwrite_confc) == 0 && cp_on == 1 ){
      if(myid==0 ) {fp_ccname = cfopen(general_data->filenames.ccname,"a"); }
       if(nproc>1){Barrier(world);}
       n = ncoef;
       for(is=1;is<=nstate_up;is++){
         igo=0;
         if((is>=istate_up_st) && (is<=istate_up_end)){
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
       if(cp_lsda==1 && nstate_dn != 0){
        for(is=1;is<=nstate_dn;is++){
          igo=0;
          if((is>=istate_dn_st) && (is<=istate_dn_end)){
             isoff=(is-istate_dn_st)*ncoef;
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
       }/* endif lsda */
       if(myid==0){
        n=9;
        fwrite(general_data->cell.hmat,sizeof(double),n,fp_ccname);
        fflush(fp_ccname);
        fclose(fp_ccname);
       }/* endif myid */
       if(nproc>1){Barrier(world);}
    }/*endif*/


/*=====================================================================*/
/* IV) Write to the atm force file                          */

  if(myid==0 ){
    if((general_data->timeinfo.itime %
        general_data->filenames.iwrite_atm_for) == 0 ){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cfname = cfopen(general_data->filenames.forcename,"a");

       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
       fprintf(fp_cfname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[1].fx[i],
              class->clatoms_pos[1].fy[i],class->clatoms_pos[1].fz[i]);
       }/*endfor*/
      for(i=0;i<3;i++) 
       fprintf(fp_cfname,"%.12g %.12g %.12g\n",general_data->cell.hmat[1+i],
              general_data->cell.hmat[4+i],general_data->cell.hmat[7+i]);
      fflush(fp_cfname);
      fclose(fp_cfname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cfname = cfopen(general_data->filenames.forcename,"a");
    n=1;
    for(i=1;i<=(class->clatoms_info.natm_tot);i++){
     fwrite(&(class->clatoms_pos[1].fx)[i],sizeof(double),n,fp_cfname);
     fwrite(&(class->clatoms_pos[1].fy)[i],sizeof(double),n,fp_cfname);
     fwrite(&(class->clatoms_pos[1].fz)[i],sizeof(double),n,fp_cfname);
    }/*endfor*/
    for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cfname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cfname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cfname);
    }/*endfor*/ 
      fflush(fp_cfname);
      fclose(fp_cfname);
   }/*endif*/
   }
  }/* endif myid==0 */
  if(nproc>1){Barrier(world); }



/*==========================================================================*/
   }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_inst_file_cp(CLASS *class,GENERAL_DATA *general_data, CP *cp,
                        double etot,double a,double b,double c,
                        double tac,double tab,double tbc)

/*==========================================================================*/

{/* begin routine */

#include "../typ_defs/typ_mask.h"

  int i;
  double inst_div;
  FILE *fp_iname;

/*==========================================================================*/
/*    A) Finish Averages                                            */


  inst_div    = (double)(general_data->filenames.iwrite_inst);    
  general_data->stat_avg.aikinet_cp     /= inst_div;
  general_data->stat_avg.aikinet_nhc_cp /= inst_div;
  general_data->stat_avg.aicp_ehart     /= inst_div;
  general_data->stat_avg.aicp_eext      /= inst_div;
  general_data->stat_avg.aicp_exc       /= inst_div;
  general_data->stat_avg.aicp_eke       /= inst_div;
  general_data->stat_avg.aicp_enl       /= inst_div;
  general_data->stat_avg.aikinet     /= inst_div;
  general_data->stat_avg.aikinet_v   /= inst_div;
  general_data->stat_avg.aikinet_nhc /= inst_div;
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
  fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
         general_data->stat_avg.aikinet,general_data->stat_avg.aikinet_v,
         general_data->stat_avg.aikinet_nhc,general_data->stat_avg.aivintert,
         general_data->stat_avg.aivintrat,general_data->stat_avg.aivol,etot);
  fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",
             general_data->stat_avg.aicella,general_data->stat_avg.aicellb,
             general_data->stat_avg.aicellc,general_data->stat_avg.aicellac,
             general_data->stat_avg.aicellab,general_data->stat_avg.aicellbc);
  for(i=1;i<=9;i+=3){
   fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.aipten[i],
              general_data->stat_avg.aipten[i+1],
              general_data->stat_avg.aipten[i+2]);
  }
  fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",
             general_data->stat_avg.kinet,general_data->stat_avg.kinet_v,
             general_data->stat_avg.kinet_nhc,general_data->stat_avg.vintert,
             general_data->stat_avg.vintrat,general_data->stat_avg.vol);
  fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",a,b,c,tac,tab,tbc);
  for(i=1;i<=9;i+=3){
    fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.apten_out[i],
              general_data->stat_avg.apten_out[i+1],
              general_data->stat_avg.apten_out[i+2]);
  }
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
  general_data->stat_avg.aikinet_v      = 0.0;
  general_data->stat_avg.aikinet_nhc    = 0.0;
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



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_header_cp(FILE *fp_dnamec,CP *cp,GENERAL_DATA *general_data,int myid,
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
/*==========================================================================*/
/* Write the header                                                         */

if(myid==0){
 if(ibinary == 0){
    fprintf(fp_dnamec,"ncoef_up, ncoef_dn, nstate_up, nstate_dn, ");
    fprintf(fp_dnamec,"dft_typ, restart_typ, time of dump\n");
 }else{

#ifdef IBM_ESSL
   strcpy(c_array1,"ibm_essl:ncoef_up ");
#endif
#ifdef IBM_NOESSL
   strcpy(c_array1,"ibm_noessl:ncoef_up ");
#endif
#ifdef LINUX
   strcpy(c_array1,"linux:ncoef_up ");
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
    if((general_data->simopts.cp+general_data->simopts.cp_wave)==1){
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
    if((general_data->simopts.cp+general_data->simopts.cp_wave)==1){
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

void write_dump_occ_cp(FILE *fp_dnamec,CP *cp,int myid,int ibinary)

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

void write_kseigs_file_cp(GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 int pi_beads        = cp->cpcoeffs_info.pi_beads;
 int nstate_up       = cp->cpcoeffs_info.nstate_up;
 int nstate_dn       = cp->cpcoeffs_info.nstate_dn;
 int myid            = cp->communicate.myid;
 int myid_state      = cp->communicate.myid_state;
 int ip_start        = cp->cpcoeffs_info.pi_beads_proc_st;
 int ip_end          = cp->cpcoeffs_info.pi_beads_proc_end;
 int cp_lsda         = cp->cpopts.cp_lsda;
 int ip,is,irecv;
 MPI_Comm world      = cp->communicate.world;
 FILE *fp_kseigs;
 double *kseig_tmp   = cp->cpscr.cpscr_ovmat.state_vec1;
 double *kseig;
 double ks_offset_tmp;


/*==========================================================================*/
/* Open dump file                                                           */

  if(myid==0){
   fp_kseigs = cfopen(general_data->filenames.ksname,"a");
  }/* endif myid */

/*==========================================================================*/
/* Write KS eigenvalues for up states                                       */

  for(ip=1;ip<=pi_beads;ip++){
    irecv = 1;
    if(ip >= ip_start && ip <= ip_end && myid_state == 0){
       irecv=0;
       kseig          = cp->cpcoeffs_pos[(ip-ip_start+1)].ksmat_eig_up;
       ks_offset_tmp  = cp->cpcoeffs_pos[(ip-ip_start+1)].ks_offset;
       for(is=1;is<=nstate_up;is++){
         kseig_tmp[is] = kseig[is];
       }/* endfor */
       if(myid != 0) {
         Ssend(&(kseig_tmp[1]),nstate_up,MPI_DOUBLE,0,0,world);
         Ssend(&(ks_offset_tmp),1,MPI_DOUBLE,0,0,world);
       }/* endif myid_bead */
    }/* endif ip */
    if(myid == 0){
     if(irecv != 0){
       Recv(&(kseig_tmp[1]),nstate_up,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
            world);
       Recv(&(ks_offset_tmp),1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
     }/* endif */
     for(is=1;is<=nstate_up;is++){
       fprintf(fp_kseigs,"%d  %.6g \n",is,kseig_tmp[is]);
     }/* endfor is */
     fprintf(fp_kseigs,"ks_offset %.6g\n",ks_offset_tmp);
   }/* endif myid */
  }/* endfor ip */


/*==========================================================================*/
/* Write KS eigenvalues for dn states                                       */

 if(cp_lsda == 1){

   for(ip=1;ip<=pi_beads;ip++){
     irecv = 1;
     if(ip >= ip_start && ip <= ip_end && myid_state == 0){
        irecv=0;
        kseig          = cp->cpcoeffs_pos[(ip-ip_start+1)].ksmat_eig_dn;
        ks_offset_tmp  = cp->cpcoeffs_pos[(ip-ip_start+1)].ks_offset;
        for(is=1;is<=nstate_dn;is++){
          kseig_tmp[is] = kseig[is];
        }/* endfor */
        if(myid != 0) {
          Ssend(&(kseig_tmp[1]),nstate_dn,MPI_DOUBLE,0,0,world);
          Ssend(&(ks_offset_tmp),1,MPI_DOUBLE,0,0,world);
        }/* endif myid_bead */
     }/* endif ip */
     if(myid == 0){
      if(irecv != 0){
        Recv(&(kseig_tmp[1]),nstate_dn,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
             world);
        Recv(&(ks_offset_tmp),1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
      }/* endif */
      for(is=1;is<=nstate_dn;is++){
        fprintf(fp_kseigs,"%d  %.6g \n",is,kseig_tmp[is]);
      }/* endfor is */
      fprintf(fp_kseigs,"%.6g\n",ks_offset_tmp);
    }/* endif myid */
   }/* endfor ip */

 }/* endif cp_lsda */

/*==========================================================================*/
/* Close file                                                              */

  if(myid==0){
    fflush(fp_kseigs);
    fclose(fp_kseigs);
  }/* endif */

/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_elf_file_cp(GENERAL_DATA *general_data,CP *cp)

/*==========================================================================*/
{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

 int index,ka,kb,kc;
 int count;
 int nkf1,nkf2,nkf3;
 int kb_str,kb_end;
 int skc_fft_ka_proc;
 int ekc_fft_ka_proc;
 int skb_fft_ka_proc;
 int ekb_fft_ka_proc;
 int skc_fft_ka_proc0;
 int ekc_fft_ka_proc0;
 int skb_fft_ka_proc0;
 int ekb_fft_ka_proc0;
 int skc_use;

 int nfft            =    cp->cp_para_fft_pkg3d_lg.nfft;
 int nfft_proc       =    cp->cp_para_fft_pkg3d_lg.nfft_proc;
 int nfft2           =    nfft/2;
 int nfft2_proc      =    nfft_proc/2;
 int nfft_proc_dens_cp_box;
 int nfft2_proc_dens_cp_box;
 int i;
 int iproc;
 static int iopen_flag=1;


 int myid             = cp->communicate.myid;
 int myid_state       = cp->communicate.myid_state;
 int np_states        = cp->communicate.np_states;
 int cp_lsda          = cp->cpopts.cp_lsda;
 int cp_ngrid_skip    = cp->cpopts.cp_ngrid_skip;
 int *kmax_cp;
 double *cp_elf_up    = cp->electronic_properties.cp_elf_up;
 double *cp_elf_dn    = cp->electronic_properties.cp_elf_dn;
 double *cp_elf_tmp   = cp->cpscr.cpscr_rho.rho_up;
 FILE *fp_elf;
 MPI_Comm world       = cp->communicate.world;
 MPI_Comm comm_states = cp->communicate.comm_states;
 

/*==========================================================================*/
/* Useful constants                                                         */

 if(cp->cpopts.cp_dual_grid_opt >= 1){
   kmax_cp         = cp->cpewald.kmax_cp_dens_cp_box;
   skc_fft_ka_proc = cp->cp_para_fft_pkg3d_dens_cp_box.skc_fft_ka_proc;
   ekc_fft_ka_proc = cp->cp_para_fft_pkg3d_dens_cp_box.ekc_fft_ka_proc;
   skb_fft_ka_proc = cp->cp_para_fft_pkg3d_dens_cp_box.skb_fft_ka_proc;
   ekb_fft_ka_proc = cp->cp_para_fft_pkg3d_dens_cp_box.ekb_fft_ka_proc;
 } else {
   kmax_cp         = cp->cpewald.kmax_cp;
   skc_fft_ka_proc = cp->cp_para_fft_pkg3d_lg.skc_fft_ka_proc;
   ekc_fft_ka_proc = cp->cp_para_fft_pkg3d_lg.ekc_fft_ka_proc;
   skb_fft_ka_proc = cp->cp_para_fft_pkg3d_lg.skb_fft_ka_proc;
   ekb_fft_ka_proc = cp->cp_para_fft_pkg3d_lg.ekb_fft_ka_proc;
 } /*endif dual grid opt on */

 /* These values are overwritten on the master processor */
 if((cp_lsda == 1)&&(myid_state==0)){
   skc_fft_ka_proc0 = skc_fft_ka_proc;
   ekc_fft_ka_proc0 = ekc_fft_ka_proc;
   skb_fft_ka_proc0 = skb_fft_ka_proc;
   ekb_fft_ka_proc0 = ekb_fft_ka_proc;
 }

 if(cp->cpopts.cp_dual_grid_opt >= 1){
   nfft_proc_dens_cp_box  = cp->cp_para_fft_pkg3d_dens_cp_box.nfft_proc;
   nfft2_proc_dens_cp_box = nfft_proc_dens_cp_box/2;
   nfft2_proc             = nfft2_proc_dens_cp_box;
 }/*endif cp_dual_grid_opt*/

  nkf1 = 4*(kmax_cp[1]+1);
  nkf2 = 4*(kmax_cp[2]+1);
  nkf3 = 4*(kmax_cp[3]+1);

/*==========================================================================*/
/* Check commensurability                                                   */

   if(((nkf1 % cp_ngrid_skip) != 0) || ((nkf2 % cp_ngrid_skip) != 0) ||
      ((nkf3 % cp_ngrid_skip) != 0)){
     if(myid_state==0){
       printf("@@@@@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("The grid skip value, cp_ngrid_skip, must be commensurate\n");
       printf("with the number of points along each direction\n");
       printf("nx = %d, ny = %d, nz = %d, grid skip = %d\n",nkf1,nkf2,nkf3,cp_ngrid_skip);
       printf("@@@@@@@@@@@@@@@@@@@@-ERROR-@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
     }/* endif myid */
       if(np_states>1){Barrier(comm_states);}
       Finalize();
       exit(1);
   }/* endif */
 

/*==========================================================================*/
/* Open elf file                                                           */

 if(iopen_flag == 1){
   if(myid_state==0){
      fp_elf = cfopen(general_data->filenames.elfname,"w");
   }/* endif myid */
   iopen_flag=0;
 } else {
   if(myid_state==0){
      fp_elf = cfopen(general_data->filenames.elfname,"a");
   }/* endif myid */
 }/* endif iopen_flag */

/*==========================================================================*/
/* Write the ELF for the up electrons                                       */

  if(np_states == 1){
    for(kc=1;kc<=nkf3;kc+=cp_ngrid_skip){
      for(kb=1;kb<=nkf2;kb+=cp_ngrid_skip){
        for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
          i = (ka-1) + (kb-1)*nkf1 + (kc-1)*nkf1*nkf2 + 1;
          fprintf(fp_elf,"%.5e\n",cp_elf_up[i]);
        }/* endfor */
     }/* endfor */
    }/* endfor */

   } else {

    for(iproc=0;iproc<np_states;iproc++){
      if(iproc==0){
        if(myid_state == 0){ 
           /* Assumes that skc_fft_ka_proc=1 for the 0th processor */
           i=1;
           for(kc=skc_fft_ka_proc;kc<=ekc_fft_ka_proc;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_fft_ka_proc ? skb_fft_ka_proc : 1);
            kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
             for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
              fprintf(fp_elf,"%.5e\n",cp_elf_up[i]);
              i += cp_ngrid_skip;
             }/* endfor */
             i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
           }/* endfor */
        }/* endif myid_state == 0 */
      }/* endif iproc == 0 */
      if(np_states>1){Barrier(comm_states);}
      if(iproc != 0){
	if(myid_state == iproc){
          Ssend(&cp_elf_up[1],nfft2_proc,MPI_DOUBLE,0,0,comm_states);
          Ssend(&skc_fft_ka_proc,1,MPI_INT,0,1,comm_states);
          Ssend(&ekc_fft_ka_proc,1,MPI_INT,0,2,comm_states);
          Ssend(&skb_fft_ka_proc,1,MPI_INT,0,3,comm_states);
          Ssend(&ekb_fft_ka_proc,1,MPI_INT,0,4,comm_states);
        } /* endif myid_state == iproc */
        if(myid_state == 0){
          Recv(&cp_elf_tmp[1],nfft2_proc,MPI_DOUBLE,MPI_ANY_SOURCE,0,comm_states);
          Recv(&skc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,1,comm_states);
          Recv(&ekc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,2,comm_states);
          Recv(&skb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,3,comm_states);
          Recv(&ekb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,4,comm_states);
          skc_use = (kb>nkf2 ? kc : kc-cp_ngrid_skip);
          i=(i%cp_ngrid_skip)+1;
          for(kc=skc_use;kc<=ekc_fft_ka_proc;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_use ? kb : 1);
            kb_str = (kb>nkf2 ? kb-nkf2 : kb);
            kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
             for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
              fprintf(fp_elf,"%.5e\n",cp_elf_tmp[i]);
              i += cp_ngrid_skip;
             }/* endfor */
             i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
          }/* endfor */
        }/* endif myid_state == 0 */
      }/* endif iproc != 0 */
      if(np_states>1){Barrier(comm_states);}
   }/* endfor iproc */
  }/* endif np_states */


/*==========================================================================*/
/* Write the ELF for the down electrons if necessary                        */

 if(cp_lsda == 1){
  if(np_states == 1){
    for(kc=1;kc<=nkf3;kc+=cp_ngrid_skip){
      for(kb=1;kb<=nkf2;kb+=cp_ngrid_skip){
        for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
          i = (ka-1) + (kb-1)*nkf1 + (kc-1)*nkf1*nkf2 + 1;
          fprintf(fp_elf,"%.5e\n",cp_elf_dn[i]);
        }/* endfor */
     }/* endfor */
    }/* endfor */

   } else {

    for(iproc=0;iproc<np_states;iproc++){
      if(iproc==0){
        if(myid_state == 0){ 
           i = 1;
           for(kc=skc_fft_ka_proc0;kc<=ekc_fft_ka_proc0;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_fft_ka_proc0 ? skb_fft_ka_proc0 : 1);
            kb_end = (kc==ekc_fft_ka_proc0 ? ekb_fft_ka_proc0 : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
             for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
              fprintf(fp_elf,"%.5e\n",cp_elf_dn[i]);
              i += cp_ngrid_skip;
             }/* endfor */
             i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
           }/* endfor */
        }/* endif myid_state == 0 */
      }/* endif iproc == 0 */
      if(np_states>1){Barrier(comm_states);}
      if(iproc != 0){
	if(myid_state == iproc){
          Ssend(&cp_elf_dn[1],nfft2_proc,MPI_DOUBLE,0,0,comm_states);
          Ssend(&skc_fft_ka_proc,1,MPI_INT,0,1,comm_states);
          Ssend(&ekc_fft_ka_proc,1,MPI_INT,0,2,comm_states);
          Ssend(&skb_fft_ka_proc,1,MPI_INT,0,3,comm_states);
          Ssend(&ekb_fft_ka_proc,1,MPI_INT,0,4,comm_states);
        } /* endif myid_state == iproc */
        if(myid_state == 0){
          Recv(&cp_elf_tmp[1],nfft2_proc,MPI_DOUBLE,MPI_ANY_SOURCE,0,comm_states);
          Recv(&skc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,1,comm_states);
          Recv(&ekc_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,2,comm_states);
          Recv(&skb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,3,comm_states);
          Recv(&ekb_fft_ka_proc,1,MPI_INT,MPI_ANY_SOURCE,4,comm_states);
          skc_use = (kb>nkf2 ? kc : kc-cp_ngrid_skip);
          i=(i%cp_ngrid_skip)+1;
          for(kc=skc_use;kc<=ekc_fft_ka_proc;kc+=cp_ngrid_skip){
            kb_str = (kc==skc_use ? kb : 1);
            kb_str = (kb>nkf2 ? kb-nkf2 : kb);
            kb_end = (kc==ekc_fft_ka_proc ? ekb_fft_ka_proc : nkf2);
            for(kb=kb_str;kb<=kb_end;kb+=cp_ngrid_skip){
             for(ka=1;ka<=nkf1;ka+=cp_ngrid_skip){
              fprintf(fp_elf,"%.5e\n",cp_elf_tmp[i]);
              i += cp_ngrid_skip; 
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1;
            }/* endfor */
            i += (cp_ngrid_skip-1)*nkf1*nkf2;
           }/* endfor */
        }/* endif myid_state == 0 */
      }/* endif iproc != 0 */
      if(np_states>1){Barrier(comm_states);}
   }/* endfor iproc */
  }/* endif np_states */
 }/* endif cp_lsda */

/*==========================================================================*/
/* Close the file                                                           */

  if(myid_state==0){fclose(fp_elf);}



/*==========================================================================*/
 }/* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_coef_cp(FILE *fp_dnamec,CP *cp,CLASS *class,int ibinary)

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
 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo;
 int i,n; 
 float cre_dum,cim_dum;

 char *c_array1,*c_array2;  /*for binary write */
 int csize = MAXWORD; 

 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;

 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int nproc      = class->communicate.np;
 MPI_Comm world = class->communicate.world;

/*=======================================================================*/
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
/* Write the coefficients                                                  */

/*-------------------------------------------------------------------------*/
/* Up states                                                               */

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
   if(nproc>1){Barrier(world);}
     for(is=1;is<=nstate_up;is++){
       igo=0;
       if((is>=istate_up_st) && (is<=istate_up_end)){
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
         Recv(&(cre_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
         Recv(&(cim_up_tmp[0]),(ncoef+1),MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,
               world);
        }/* endif */
        if(myid==0){
         n = 1;
         for(i=1;i<=ncoef;i++){
         if(ibinary == 0){
          fprintf(fp_dnamec,"%.8g %.8g\n",cre_up_tmp[i],cim_up_tmp[i]);
         }else{
           cre_dum = (float)cre_up_tmp[i];
           cim_dum = (float)cim_up_tmp[i];
           fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
           fwrite(&(cim_dum),sizeof(float),n,fp_dnamec);
         }/*endif ibinary */
         }/*endfor: write*/
        }/* endif myid */
        if(nproc>1){Barrier(world);}
     }/*endfor:states*/

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
   if(nproc>1){Barrier(world);}
     for(is=1;is<=nstate_dn;is++){
       igo=0;
       if((is>=istate_dn_st) && (is<=istate_dn_end)){
          isoff=(is-istate_dn_st)*ncoef;
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
           fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
           fwrite(&(cim_dum),sizeof(float),n,fp_dnamec);
         }/*endif ibinary */
         }/*endfor: write*/
        }/* endif myid */
        if(nproc>1){Barrier(world);}
     }/*endfor:states*/
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

void write_dump_vcoef_cp(FILE *fp_dnamec,CP *cp,CLASS *class,
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
 int cp_lsda        =cp->cpopts.cp_lsda;
 int nrecv,isoff;
 int is,igo;
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
 int nproc      = class->communicate.np;
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

 if((general_data->simopts.cp+general_data->simopts.cp_wave)==1){

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
   if(nproc>1){Barrier(world);}
     for(is=1;is<=nstate_up;is++){
       igo=0;
       if((is>=istate_up_st) && (is<=istate_up_end)){
          isoff=(is-istate_up_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_up_tmp[i] = cp->cpcoeffs_pos[1].vcre_up[(i+isoff)];
            cim_up_tmp[i] = cp->cpcoeffs_pos[1].vcim_up[(i+isoff)];
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
        if(nproc>1){Barrier(world);}
     }/*endfor:states*/

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
    if(nproc>1){Barrier(world);}
     for(is=1;is<=nstate_dn;is++){
       igo=0;
       if((is>=istate_dn_st) && (is<=istate_dn_end)){
          isoff=(is-istate_dn_st)*ncoef;
          for(i=1;i<=ncoef;i++){
            cre_dn_tmp[i] = cp->cpcoeffs_pos[1].vcre_dn[(i+isoff)];
            cim_dn_tmp[i] = cp->cpcoeffs_pos[1].vcim_dn[(i+isoff)];
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
        if(nproc>1){Barrier(world);}
     }/*endfor:states*/
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

void write_dump_extended_cp(FILE *fp_dnamec,CP *cp,CLASS *class,
                            GENERAL_DATA *general_data,int ibinary)

/*==========================================================================*/

{/* begin routine */
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

 int i,j,ncoef_up,ncoef_dn;
 int ichain,inhc ;
 int pi_beads = class->clatoms_info.pi_beads;
 int iproc,ioff_re,ioff_im;
 int nstate_proc_now;

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

 float cre_dum,cim_dum;
 
 MPI_Status stat;

/* Local pointers */
 double *cre_up_tmp = cp->cpscr.cpscr_wave.cre_up;
 double *cim_up_tmp = cp->cpscr.cpscr_wave.cim_up;
 double *cre_dn_tmp = cp->cpscr.cpscr_wave.cre_dn;
 double *cim_dn_tmp = cp->cpscr.cpscr_wave.cim_dn;
 int myid       = class->communicate.myid;
 int np_states  = class->communicate.np_states;
 int nproc      = class->communicate.np;
 MPI_Comm world = class->communicate.world;

 char *c_array1,*c_array2;
 int csize = MAXWORD;
 int n;


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

  if((general_data->simopts.cp+general_data->simopts.cp_wave)==1){
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

   for(ichain=1;ichain<=len_c_nhc;ichain++){

     if(cp->cptherm_info.istate_nhc_opt==4){
       for(iproc=0;iproc<np_states;iproc++){
         nstate_proc_now=0;
         if(myid==iproc){         
           if(myid!=0){Ssend(&(nstate_up_proc),1,MPI_INT,0,0,world);}
           nstate_proc_now=nstate_up_proc;
        }/*endif*/
         if(myid==0&&iproc!=0){         
           Recv(&(nstate_proc_now),1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,world);
        }/*endif*/
        for(is=1;is<=nstate_proc_now;is++){
           if(myid==iproc){
             ioff_re = (is-1)*ncoef;
             ioff_im = (nstate_proc_now+is-1)*ncoef;
             for(inhc=1;inhc<=ncoef;inhc++){
               cre_up_tmp[inhc] = 
                          cp->cptherm_pos[1].vc_nhc[ichain][inhc+ioff_re];
             }/* endfor */
             for(inhc=1;inhc<=ncoef;inhc++){
               cre_up_tmp[inhc+ncoef] = 
                           cp->cptherm_pos[1].vc_nhc[ichain][inhc+ioff_im];
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
             }/*endfor*/
           }/* endif */
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
               ioff_re = (is-1+nstate_up_proc*2)*ncoef;
               ioff_im = (nstate_proc_now+is-1+nstate_up_proc*2)*ncoef;
               for(inhc=1;inhc<=ncoef;inhc++){
                 cre_up_tmp[inhc] = 
                     cp->cptherm_pos[1].vc_nhc[ichain][inhc+ioff_re];
               }/* endfor */
               for(inhc=1;inhc<=ncoef;inhc++){
                 cre_up_tmp[inhc+ncoef] = 
                           cp->cptherm_pos[1].vc_nhc[ichain][inhc+ioff_im];
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
         if(nproc>1){Barrier(world);}
       }/* endfor : iproc*/
     }/*endif:istate_nhc_opt*/

     if(cp->cptherm_info.istate_nhc_opt<=3&&myid==0){
       for(inhc=1;inhc<=num_c_nhc;inhc++){
	 if(ibinary == 0){
         fprintf(fp_dnamec,"%.8g\n",cp->cptherm_pos[1].vc_nhc[ichain][inhc]);
         }else{
          n = 1;
         cre_dum = (float) cp->cptherm_pos[1].vc_nhc[ichain][inhc];
         fwrite(&(cre_dum),sizeof(float),n,fp_dnamec);
         }/*endif ibinary*/
       }/*endfor*/
     }/*endif:istate_nhc_opt*/

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






