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
/*               Header:                                                    */

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void output_min(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

/*==========================================================================*/
  {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int iii;

  int myid    = class->communicate.myid;
  int np_forc = class->communicate.np_forc;

/*=====================================================================*/
/*        Output_md routine                                            */
/*=====================================================================*/
/*  0)Open File option :                                               */

  if((general_data->timeinfo.itime==0)&&(myid==0)&&
     (general_data->filenames.ifile_open==1)){
    initial_fopen_md(class,general_data,bonded);
  }/*endif*/

/*=======================================================================*/
/*  II) Write Initial Energies to screen                                 */

  if( (general_data->timeinfo.itime==0) && (myid==0) &&
      (general_data->simopts.minimize == 1)){
     initial_output_min(class,general_data,bonded);
  }/*endif*/

/*=======================================================================*/
/*  IV) Write the output to screen                                       */

  if((general_data->timeinfo.itime != 0)&& (myid==0) ){
   if((general_data->timeinfo.itime%general_data->filenames.iwrite_screen)==0){
      screen_write_min(class,general_data,bonded);
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* V) Write to the output to dump file                                  */ 

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
   if((general_data->timeinfo.itime % general_data->filenames.iwrite_dump)==0){
    write_dump_file_md(class,bonded,general_data);
   }/*endif*/
  }/*endif*/

/*======================================================================*/
/* VII) Write to the config files                                      */

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
    write_config_files_md(class,bonded,general_data);
  }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_output_min(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

/*==========================================================================*/

{/* begin routine */

  int np_tot,npairs,nsh_tot;
  double nhc_div;
  double eu_conv=1.0;

  if(general_data->filenames.iwrite_units==0){eu_conv=1.0;}
  if(general_data->filenames.iwrite_units==1){eu_conv=KCAL;}
  if(general_data->filenames.iwrite_units==2){eu_conv=BOLTZ;}
     

/*==========================================================================*/
  printf("Initial intra        energy  %.10g\n",
                           general_data->stat_avg.vintrat*eu_conv);
  printf("Initial inter        energy  %.10g\n",
                           general_data->stat_avg.vintert*eu_conv);
  printf("Initial kinetic      energy  %.10g\n",
                           general_data->stat_avg.kinet*eu_conv);
  if(class->energy_ctrl.isep_vvdw == 1) {
    printf("Initial VDW          energy  %.10g\n",
                           general_data->stat_avg.vvdw*eu_conv);
    printf("Initial Coul         energy  %.10g\n",
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
  printf("Initial Torsion      energy  %.10g\n",
                                     general_data->stat_avg.vtorst*eu_conv);
  printf("Initial onefour      energy  %.10g\n",
                                     general_data->stat_avg.vonfot*eu_conv);
  printf("Initial surface      energy  %.10g\n",
                                     general_data->stat_avg.vsurft*eu_conv);
  printf("Initial volume               %.10g\n",general_data->stat_avg.vol
                                                 *BOHR*BOHR*BOHR);
  printf("Initial Atm temperature      %.10g\n",
        ((2.0*general_data->stat_avg.kinet*BOLTZ)/
         ((double)(class->clatoms_info.nfree))));

  if((general_data->cell.iperd==3)&&(general_data->simopts.md==1)){
    printf("Initial Pint (Cnstr-cntrb=0) %.10g\n",
       (((general_data->ptens.tvten)[1]+(general_data->ptens.pvten_tot)[1]
        +(general_data->ptens.tvten)[5]+(general_data->ptens.pvten_tot)[5]
        +(general_data->ptens.tvten)[9]+(general_data->ptens.pvten_tot)[9])
                   /(3.0*(general_data->stat_avg.vol)*PCONV)));
  }/*endif*/

/*==========================================================================*/
  /* NHC quantities                                                        */

    if(general_data->ensopts.nvt==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc) );
      printf("Initial NHC temperature      %.10g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      if(class->clatoms_info.pi_beads>1){
        nhc_div = (double) ( class->therm_info_bead.len_nhc*
                           (class->therm_info_bead.num_nhc) );
        printf("Initial bead NHC temperature %.10g\n", 
              ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      }/*endif*/
    }/*endif*/

    if(general_data->ensopts.npt_i==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature      %.10g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      nhc_div = (double) ( class->therm_info_bead.len_nhc*
                          (class->therm_info_bead.num_nhc +1));
      if(class->clatoms_info.pi_beads>1){
       printf("Initial bead NHC temperature %.10g\n", 
             ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      }/*endif*/
      printf("Initial Volume temperature   %.10g\n",
            (2.0*general_data->stat_avg.kinet_v*BOLTZ)); 
    }/*endif*/

    if(general_data->ensopts.npt_f==1){
      nhc_div = (double) ( class->therm_info_class.len_nhc*
                          (class->therm_info_class.num_nhc +1));
      printf("Initial NHC temperature %.10g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      nhc_div = (double) ( class->therm_info_bead.len_nhc*
                          (class->therm_info_bead.num_nhc +1));
      if(class->clatoms_info.pi_beads>1){
       printf("Initial bead NHC temperature %.10g\n", 
            ((2.0*general_data->stat_avg.kinet_nhc*BOLTZ)/nhc_div));
      }/*endif*/
      printf("Initial Volume temperature %.10g\n",      
            ((2.0*general_data->stat_avg.kinet_v*BOLTZ)/9.0));
    }/*endif*/

/*==========================================================================*/
  /* Neighbor list quantities                                              */

    if((class->nbr_list.iver)==1){
      npairs = class->nbr_list.verlist.nver_lst_now;
      np_tot = (class->clatoms_info.natm_tot)*
               (class->clatoms_info.natm_tot-1)/2 
       - (bonded->excl.nlst);
      printf("Number of pairs   = %d out of %d\n",npairs,np_tot);
      if(general_data->timeinfo.int_res_ter==1){
        npairs = class->nbr_list.verlist.nver_lst_now_res;
       printf("Number RESPA pairs = %d out of %d\n",npairs,np_tot);
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


/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_fopen_min(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

/*==========================================================================*/

{/* begin routine */
/*==========================================================================*/
/*   Local Variables                                                        */

  int n=1,iii,ibinary,iwrite_now;
  NAME file_typ;
  FILE *fp_cpname,*fp_dname;

/*==========================================================================*/
/*     A) Open dump file                                          */

  fp_dname = cfopen(general_data->filenames.dname,"w");
  fclose(fp_dname);

/*==========================================================================*/
/*     C) Open pos conf file                                           */

  ibinary    = general_data->filenames.iwrite_conf_binary;
  iwrite_now = general_data->filenames.iwrite_confp;
  strcpy(file_typ,"pos_file");
  fp_cpname  = cfopen(general_data->filenames.cpname,"w");
  write_gen_header(class,general_data,fp_cpname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_cpname); 

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void screen_write_min(CLASS *class,GENERAL_DATA *general_data,
                      BONDED *bonded)

/*==========================================================================*/
{/* begin routine */

 int npairs,ip;
 double atime,atm_div;
 double updates_t,updates_true,updates_now; 
 int pi_beads = class->clatoms_info.pi_beads;
 static int write_now = 1;
 double fc_mag_up = general_data->stat_avg.fc_mag_up;
 double fc_max_up = general_data->stat_avg.fc_max_up;
 double fc_mag_dn = general_data->stat_avg.fc_mag_dn;
 double fc_max_dn = general_data->stat_avg.fc_max_dn;
 double fatm_mag  = general_data->stat_avg.fatm_mag;
 double fatm_max  = general_data->stat_avg.fatm_max;
 double eu_conv;
 double vtot;

/*=======================================================================*/
/* Write to screen                                                       */
  
  atime = (double)(general_data->timeinfo.itime);
  if(general_data->filenames.iwrite_units==0){eu_conv=1.0;}
  if(general_data->filenames.iwrite_units==1){eu_conv=KCAL;}
  if(general_data->filenames.iwrite_units==2){eu_conv=BOLTZ;}

/*=======================================================================*/

   vtot = (general_data->stat_avg.vintert) + 
          (general_data->stat_avg.vintrat);

/*=======================================================================*/
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
     printf("Time step         = %g\n",(double)(general_data->timeinfo.itime));
     printf("----------------- \n");
     if(write_now==1){
      printf("Atm force        = %g           %g\n",fatm_mag,fatm_max);
     }/*endif*/
     if(write_now==1){
      printf("----------------- \n");
      printf("Atm Energy        = %g \n",(vtot));
      printf("Intermol Atm PE   = %g \n",(general_data->stat_avg.vintert));
      printf("Intramol Atm PE   = %g \n",(general_data->stat_avg.vintrat));

      if(class->energy_ctrl.isep_vvdw == 1) {
       printf("Initial VDW        energy  %g\n",general_data->stat_avg.vvdw);
       printf("Initial Coul       energy  %g\n",general_data->stat_avg.vcoul);
      }/*endif*/

      printf("Bond     energy  %g\n",general_data->stat_avg.vbondt);
      printf("Bend     energy  %g\n",general_data->stat_avg.vbendt);
      printf("Bendbnd  energy  %g\n",general_data->stat_avg.vbend_bndt);
      printf("Torsion  energy  %g\n",general_data->stat_avg.vtorst);
      printf("Onefour  energy  %g\n",general_data->stat_avg.vonfot);
    }/*endif*/

      printf("----------------- \n");
      atm_div = (double)(class->clatoms_info.nfree);
      printf("Atm Deg. Free     = %g\n",atm_div); 


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
      write_now=0;
      if(general_data->simopts.minimize==1) write_now=1;

  fflush(stdout);

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/



