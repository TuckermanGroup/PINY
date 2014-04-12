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

void output_md(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

/*==========================================================================*/
  {/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int iii;
  int exit_flag = general_data->timeinfo.exit_flag;
  double etot,econv_now; 
  double vtot,avtot,pnow; 
  double pnow_inter,apnow_inter,pnow_intra,apnow_intra;
  double pnow_kin,apnow_kin;
  double apnow;
  double nhc_div;
  double a,b,c,tab,tac,tbc; 
  double deth;

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
/* I) All gather velocities if necessary */ 

 if(general_data->timeinfo.itime!=0 && (np_forc>1)){
  if(((general_data->timeinfo.itime%general_data->filenames.iwrite_confv)==0)||
     ((general_data->timeinfo.itime%general_data->filenames.iwrite_dump)==0)){
    forc_level_vel_gather(class);
  }/*endif*/
 }/*endif*/

/*=======================================================================*/
/*  II) Write Initial Energies to screen                                 */

  if( (general_data->timeinfo.itime==0) && (myid==0) &&
      (general_data->simopts.debug == 0)){
     initial_output_md(class,general_data,bonded);
  }/*endif*/

/*=======================================================================*/
/* III) Calculate some dinky quantities                                   */

  if((general_data->timeinfo.itime!=0) && (myid==0) ){
   if(((general_data->timeinfo.itime%general_data->filenames.iwrite_screen)==0)
    ||((general_data->timeinfo.itime%general_data->filenames.iwrite_inst)==0)
    || (exit_flag == 1)){
      get_cell(general_data->cell.hmat,&a,&b,&c,&tab,&tbc,&tac);
      etot = general_data->stat_avg.econv_now;
      dink_quant_calc_md(class, general_data,&etot,&econv_now,&vtot, &avtot,
                         &deth,&pnow, &apnow, &pnow_inter, &apnow_inter,
                         &pnow_intra,&apnow_intra,&pnow_kin,&apnow_kin,
                         &nhc_div);
    }/*endif*/
  }/*endif*/

/*=======================================================================*/
/*  IV) Write to the output to screen                                   */

  if((general_data->timeinfo.itime != 0)&& (myid==0) ){
   if((general_data->timeinfo.itime%general_data->filenames.iwrite_screen)==0 ||
      (exit_flag == 1)){
      screen_write_md(class,general_data,bonded,
                      etot,econv_now,vtot,avtot,deth,pnow, apnow,
                      pnow_inter, apnow_inter, pnow_intra, apnow_intra,
                      pnow_kin, apnow_kin, 
                      nhc_div,a,b,c,tab,tbc,tac);
    }/*endif*/
  }/*endif*/

/*======================================================================*/
/* V) Write to the output to dump file                                  */ 

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
   if((general_data->timeinfo.itime % general_data->filenames.iwrite_dump)==0 ||
      (exit_flag == 1)){
    write_dump_file_md(class,bonded,general_data);
   }/*endif*/
  }/*endif*/

/*======================================================================*/
/* VI) Write to the output to free energy                                */ 

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
   if((general_data->timeinfo.itime % general_data->filenames.iwrite_dump)==0){
     write_free_energy_file(class,bonded,general_data);
   }/*endif*/
  }/*endif*/

/*====================================================================*/
/* VII) Write to the config files                                      */

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
    write_config_files_md(class,bonded,general_data);
  }/*endif*/

/*======================================================================*/
/* VIII) Write to the inst avgs to inst file                            */

  if((general_data->timeinfo.itime!=0)&& (myid==0) ){
   if((general_data->timeinfo.itime % general_data->filenames.iwrite_inst)==0){
      write_inst_file_md(class,general_data,etot,a,b,c,tac,tab,tbc);
    }/*endif*/
  }/*endif*/

/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void initial_output_md(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

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
            ((2.0*general_data->stat_avg.kinet_v*BOLTZ)/6.0));
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

void initial_fopen_md(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded)

/*==========================================================================*/

{/* begin routine */
/*==========================================================================*/
/*   Local Variables                                                        */

  int n=1,iii,ibinary,iwrite_now;
  NAME file_typ;
  FILE *fp_bond_free,*fp_bend_free,*fp_tors_free,*fp_rbar_free;
  FILE *fp_iname, *fp_cpname, *fp_cvname,*fp_dname,*fp_cfname;

/*==========================================================================*/
/*     A) Open dump file                                          */

  fp_dname = cfopen(general_data->filenames.dname,"w");
  fclose(fp_dname);

/*==========================================================================*/
/*     A) Open inst avg file                                          */

  ibinary    = 0; 
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
  write_gen_header(class,general_data,fp_cpname,ibinary,
                   iwrite_now,file_typ);
  fclose(fp_cpname); 

/*======================================================================*/
/*     E) Open partial pos conf file                                    */

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
/*     F) Open bond free energy file                                   */

    if(bonded->bond_free.num>0){
      fp_bond_free  = cfopen(bonded->bond_free.file,"w");
      fclose(fp_bond_free);
    }/*endif*/

/*==========================================================================*/
/*    G) Open bend free energy file                                    */

    if(bonded->bend_free.num>0){
      fp_bend_free  = cfopen(bonded->bend_free.file,"w");
      fclose(fp_bend_free);
    }/*endif*/

/*==========================================================================*/
/*     H) Open tors free energy file                                    */

    if(bonded->tors_free.num>0){
      fp_tors_free  = cfopen(bonded->tors_free.file,"w");
      fclose(fp_tors_free);
    }/*endif*/
    
/*==========================================================================*/
/*     I) Open rbar_sig free energy file                                    */

    if(bonded->rbar_sig_free.nfree>0){
      fp_rbar_free  = cfopen(bonded->rbar_sig_free.file,"w");
      fclose(fp_rbar_free);
    }/*endif*/
    
/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void dink_quant_calc_md(CLASS *class, GENERAL_DATA *general_data,double *etot,
                        double *econv_now,double *vtot, double *avtot,
                        double *deth,double *pnow, double *apnow, 
                        double *pnow_inter, double *apnow_inter,
                        double *pnow_intra, double *apnow_intra,
                        double *pnow_kin, double *apnow_kin,
                        double *nhc_div)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

  int i,iii;

/*=======================================================================*/
/*  Energy */

  (*econv_now) = fabs((*etot)-general_data->stat_avg.econv0)/
                 fabs(general_data->stat_avg.econv0);
  if(general_data->ensopts.nvt_isok==1) *econv_now=*etot;
/*=======================================================================*/
/*  Volume and pressure                                                    */
  
  (*vtot) = general_data->stat_avg.vintert + general_data->stat_avg.vintrat;
  (*avtot) = general_data->stat_avg.avintert + general_data->stat_avg.avintrat;
  (*deth) = getdeth(general_data->cell.hmat);
  for(i=1;i<=9;i++){
    general_data->stat_avg.apten_out[i]=(general_data->ptens.pvten_tot[i]
    +general_data->ptens.tvten[i])/((*deth)*PCONV);}
  (*pnow)  = (general_data->stat_avg.apten_out[1]
             +  general_data->stat_avg.apten_out[5]
             +  general_data->stat_avg.apten_out[9])/(3.0);
  (*apnow) =general_data->stat_avg.apress/(PCONV*((double)(general_data->timeinfo.itime)));
  (*pnow_inter)  = general_data->stat_avg.press_inter/PCONV;
  (*apnow_inter) = general_data->stat_avg.apress_inter
                       /(PCONV*((double)(general_data->timeinfo.itime)));
  (*pnow_intra)  = general_data->stat_avg.press_intra/PCONV;
  (*apnow_intra) = general_data->stat_avg.apress_intra
                       /(PCONV*((double)(general_data->timeinfo.itime)));
  (*pnow_kin)  = general_data->stat_avg.press_kin/PCONV;
  (*apnow_kin) = general_data->stat_avg.apress_kin
                       /(PCONV*((double)(general_data->timeinfo.itime)));

  /*======================================================================*/
  /*  NHC degrees of freedom                                              */
  
  (*nhc_div) = (double)((class->therm_info_class.len_nhc)*
                        (class->therm_info_class.num_nhc));
  if(general_data->ensopts.npt_f == 1 || general_data->ensopts.npt_i == 1){
    (*nhc_div) = (double) (class->therm_info_class.len_nhc*
                        (class->therm_info_class.num_nhc +1) );
  }/*endif*/
  
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void screen_write_md(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                     double etot,double econv_now,double vtot, 
                     double avtot,double deth,double pnow, double apnow,
                     double pnow_inter,double apnow_inter,
                     double pnow_intra,double apnow_intra,
                     double pnow_kin,double apnow_kin,
                     double nhc_div,double a,double b,double c,
                     double tab,double tbc,double tac)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

  int npairs,iprint;
  double atime,vol_div,atm_div,eu_conv=1.0;
  double updates_t,updates_true,updates_now; 

  char hat[3], hat_p[3];
  iprint = 0;strcpy(hat,"  ");strcpy(hat_p,"  ");
  if(general_data->timeinfo.iget_pe_real_inter_freq>1){
   iprint=1;strcpy(hat,"^*");
   if(general_data->timeinfo.int_res_ter==1){strcpy(hat_p,"^*");}
  }

/*=======================================================================*/
/* Write to screen                                                       */
  
  atime = (double)(general_data->timeinfo.itime);
  if(general_data->filenames.iwrite_units==0){eu_conv=1.0;}
  if(general_data->filenames.iwrite_units==1){eu_conv=KCAL;}
  if(general_data->filenames.iwrite_units==2){eu_conv=BOLTZ;}

/*=======================================================================*/
/*     A) Standard                                                   */

  printf("\n");
  printf("********************************************************\n");
  printf("QUANTITY          =  INSTANTANEOUS    AVERAGE           \n");
  printf("--------------------------------------------------------\n");
  if(general_data->ensopts.nve==1)  		printf("Ensemble          = NVE     \n");
  if(general_data->ensopts.nvt==1)  		printf("Ensemble          = NVT     \n");
  if(general_data->ensopts.nvt_isok==1)		printf("Ensemble          = NVT-ISOK\n");
  if(general_data->ensopts.npt_i==1)		printf("Ensemble          = NPT-ISO \n");
  if(general_data->ensopts.npt_f==1)		printf("Ensemble          = NPT-FLEX\n");
  printf("Time step         = %g\n",(double)(general_data->timeinfo.itime));
  printf("----------------- \n");
  atm_div = (double)(class->clatoms_info.nfree);
  if(general_data->ensopts.nvt_isok==1)	{
	/* printf("Isok Constraint%s = %g \n",hat,(  general_data->stat_avg.vpotnhc )); */

	 printf("Isok Conv%s       = %g   %g\n",hat,( general_data->stat_avg.isokconv_now ),
	        (( general_data->stat_avg.isokconv)/atime));
     printf("Ex. Sys. Temp. %s = %g       %g\n",hat,
		         (general_data->stat_avg.kinet_nhc)*eu_conv,
		         (general_data->stat_avg.akinet_nhc)/atime*eu_conv);
	 printf("----------------- \n");

  }
  else{
  printf("Econv%s           = %g %g\n",hat,( econv_now ),
        ( general_data->stat_avg.econv/atime));
  }
  printf("Energy%s          = %g %g\n",hat,
        (general_data->stat_avg.kinet+vtot)*eu_conv,
        (general_data->stat_avg.akinet+avtot)/atime*eu_conv);
  printf("Total PE%s        = %g %g\n",hat,
        (vtot)*eu_conv,(avtot/atime)*eu_conv);
  printf("Intermol PE%s     = %g %g\n",hat,
        (general_data->stat_avg.vintert)*eu_conv,
        (general_data->stat_avg.avintert/atime)*eu_conv);
  printf("Intramol PE       = %g %g\n",
        (general_data->stat_avg.vintrat)*eu_conv,
        (general_data->stat_avg.avintrat/atime)*eu_conv);
  printf("Atm KE            = %g %g\n",(general_data->stat_avg.kinet)*eu_conv,
        (general_data->stat_avg.akinet/atime)*eu_conv);
  printf("----------------- \n");
  printf("Atm Deg. Free     = %g\n",atm_div); 
  printf("Atm Temperature   = %g %g \n",
        (general_data->stat_avg.kinet*2.0*BOLTZ/atm_div),
        (general_data->stat_avg.akinet*2.0*BOLTZ/(atm_div*atime)));
  
/*======================================================================*/
/*     B) Extended Class                                          */

  if((general_data->ensopts.nvt + general_data->ensopts.npt_i
      + general_data->ensopts.npt_f == 1) && (general_data->simopts.md==1)){
    printf("NHC Temperature   = %g %g\n",
          (general_data->stat_avg.kinet_nhc*2.0*BOLTZ/nhc_div),
          (general_data->stat_avg.akinet_nhc*2.0*BOLTZ/(nhc_div*atime)));}
  if((general_data->ensopts.npt_i +general_data->ensopts.npt_f == 1) && 
     (general_data->simopts.md==1)) {
    vol_div = 1.0;
    if(general_data->ensopts.npt_f == 1) vol_div = 6.0;
    printf("Vol Temperature   = %g %g\n",
          (general_data->stat_avg.kinet_v*2.0*BOLTZ/vol_div),
          (general_data->stat_avg.akinet_v*2.0*BOLTZ/(atime*vol_div)));}
  
  printf("----------------- \n");

/*======================================================================*/
/*     C)Pressure/Vol                                                */

  if(general_data->cell.iperd>=2){
    printf("Total Pressure%s = %g %g\n",hat_p,(pnow),( apnow));
    printf("Inter Pressure%s = %g %g\n",hat_p,(pnow_inter),( apnow_inter));
    printf("Intra Pressure%s = %g %g\n",hat_p,(pnow_intra),( apnow_intra));
    printf("Kinetic Pressure = %g %g\n",(pnow_kin),( apnow_kin));
    printf("Avg  P11,P22,P33%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten[1]/(PCONV*atime),
          general_data->stat_avg.apten[5]/(PCONV*atime),
          general_data->stat_avg.apten[9]/(PCONV*atime));
    printf("Inst P11,P22,P33%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten_out[1],
          general_data->stat_avg.apten_out[5],
          general_data->stat_avg.apten_out[9]);
    printf("Avg  P12,P13,P23%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten[4]/(PCONV*atime),
          general_data->stat_avg.apten[7]/(PCONV*atime),
          general_data->stat_avg.apten[8]/(PCONV*atime));
    printf("Inst P12,P13,P23%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten_out[4],
          general_data->stat_avg.apten_out[7],
          general_data->stat_avg.apten_out[8]);
    printf("Avg  P21,P31,P32%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten[2]/(PCONV*atime),
          general_data->stat_avg.apten[3]/(PCONV*atime),
          general_data->stat_avg.apten[6]/(PCONV*atime));
    printf("Inst P21,P31,P32%s = %g %g %g \n",hat_p,
          general_data->stat_avg.apten_out[2],
          general_data->stat_avg.apten_out[3],
          general_data->stat_avg.apten_out[6]);
    printf("----------------- \n");
    printf("Volume            = %g %g\n",(deth*BOHR*BOHR*BOHR),
          ((general_data->stat_avg.avol*BOHR*BOHR*BOHR)/atime));
    printf("Inst cell lths    = %g %g %g\n",(a*BOHR),(b*BOHR),(c*BOHR));
    printf("Avg  cell lths    = %g %g %g\n",    
          ((general_data->stat_avg.acella*BOHR)/atime),
          ((general_data->stat_avg.acellb*BOHR)/atime),
          ((general_data->stat_avg.acellc*BOHR)/atime));
    printf("Inst cell angs    = %g %g %g\n",(tab),(tac),(tbc));
    printf("Avg  cell angs    = %g %g %g\n",
          (general_data->stat_avg.acellab/atime),
          (general_data->stat_avg.acellac/atime),
          (general_data->stat_avg.acellbc/atime));
    printf("----------------- \n");
  }/*endif*/
  
/*=======================================================================*/
/*     D)Constraint                                                  */

  if(bonded->constrnt.iconstrnt == 1) {
    if(general_data->simopts.md==1){
      if(bonded->bond.ncon > 0) {
       printf("Shake iter        = %g %g\n",
              (double)(general_data->stat_avg.iter_shake),
              (general_data->stat_avg.aiter_shake/atime));
       printf("Rattle iter       = %g %g\n",
              (double)(general_data->stat_avg.iter_ratl),
              (general_data->stat_avg.aiter_ratl/atime));
      }
      if(bonded->grp_bond_con.num_21 > 0) {
       printf("Grp_21 shake iter = %g %g\n",
              general_data->stat_avg.iter_21,
              general_data->stat_avg.aiter_21/atime);
       printf("Grp_21 ratl iter = %g %g\n",
              general_data->stat_avg.iter_21r,
              general_data->stat_avg.aiter_21r/atime);
      }
      if(bonded->grp_bond_con.num_23 > 0) {
       printf("Grp_23 shake iter = %g %g\n",
              general_data->stat_avg.iter_23,
              general_data->stat_avg.aiter_23/atime);
       printf("Grp_23 ratl iter = %g %g\n",
              general_data->stat_avg.iter_23r,
              general_data->stat_avg.aiter_23r/atime);
      }
      if(bonded->grp_bond_con.num_33 > 0) {
       printf("Grp_33 shake iter = %g %g\n",
              general_data->stat_avg.iter_33,
              general_data->stat_avg.aiter_33/atime);
       printf("Grp_33 ratl iter = %g %g\n",
              general_data->stat_avg.iter_33r,
              general_data->stat_avg.aiter_33r/atime);
      }
      if(bonded->grp_bond_con.num_43 > 0) {
       printf("Grp_43 shake iter = %g %g\n",
              general_data->stat_avg.iter_43,
              general_data->stat_avg.aiter_43/atime);
       printf("Grp_43 ratl iter = %g %g\n",
              general_data->stat_avg.iter_43r,
              general_data->stat_avg.aiter_43r/atime);
      }
      if(bonded->grp_bond_con.num_46 > 0) {
       printf("Grp_46 shake iter = %g %g\n",
              general_data->stat_avg.iter_46,
              general_data->stat_avg.aiter_46/atime);
       printf("Grp_46 ratl iter = %g %g\n",
              general_data->stat_avg.iter_46r,
              general_data->stat_avg.aiter_46r/atime);
      }
      printf("----------------- \n");
    }/*endif*/
  }/*endif*/
  
/*=======================================================================*/
/*     D)Misc                                                         */

  if(class->nbr_list.iver == 1){
    updates_t = general_data->stat_avg.updates;
    if(updates_t == 0) updates_t = 1;
    updates_now = (double ) (general_data->timeinfo.itime-
                          general_data->stat_avg.itime_update);
    updates_true = updates_t;
    npairs = class->nbr_list.verlist.nver_lst_now;
    printf("Inst steps/update = %g\n",updates_now);
    printf("Avg. steps/update = %g\n",(atime/updates_t));
    printf("Total list updates= %g\n",updates_true);
    printf("Number of pairs   = %d \n",npairs);
    if(general_data->timeinfo.int_res_ter==1){
      npairs = class->nbr_list.verlist.nver_lst_now_res;
      printf("Number RESPA pairs= %d \n",npairs);
    }/*endif*/
    printf("----------------- \n");
  }/*endif*/
  printf("Cpu time          = %g %g\n",(general_data->stat_avg.cpu_now),
        (general_data->stat_avg.acpu/atime));
  printf("--------------------------------------------------------\n");
  if(iprint==1){
   printf("^* Calculation frequency of position dependence = %d\n",
             general_data->timeinfo.iget_pe_real_inter_freq); 
  }/*endif*/
  printf("********************************************************\n");
  fflush(stdout);

/*--------------------------------------------------------------------------*/
  }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_dump_file_md(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

 int i,j,ip;
 int pi_beads = class->clatoms_info.pi_beads;
 double deth;
 FILE *fp_dname;

/*==========================================================================*/
/* Open dump file                                                           */

  fp_dname = cfopen(general_data->filenames.dname,"o");

/*==========================================================================*/
/*     A)Atm positions                                               */
      fprintf(fp_dname,"natm_tot restart_typ itime pi_beads\n");
      fprintf(fp_dname,"%d restart_all %d %d\n",class->clatoms_info.natm_tot,
             general_data->timeinfo.itime,pi_beads);
      fprintf(fp_dname,"atm pos, atm_typ, mol_typ mol_num\n");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
       fprintf(fp_dname,"%.13g %.13g %.13g %s %s %s %d %d\n",
              class->clatoms_pos[ip].x[i],class->clatoms_pos[ip].y[i],
                class->clatoms_pos[ip].z[i],
              class->atommaps.atm_typ[class->atommaps.iatm_atm_typ[i]],
              class->atommaps.res_typ[class->atommaps.iatm_res_typ[i]],
              class->atommaps.mol_typ[class->atommaps.iatm_mol_typ[i]],
              class->atommaps.iatm_mol_num[i],
              class->atommaps.iatm_res_num[i]);
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
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat[1],
             general_data->cell.hmat[4],general_data->cell.hmat[7]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat[2],
             general_data->cell.hmat[5],general_data->cell.hmat[8]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat[3],
             general_data->cell.hmat[6],general_data->cell.hmat[9]);
      fprintf(fp_dname,"h matrix for Ewald setup\n");
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat_ewd[1],
             general_data->cell.hmat_ewd[4],general_data->cell.hmat_ewd[7]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat_ewd[2],
             general_data->cell.hmat_ewd[5],general_data->cell.hmat_ewd[8]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->cell.hmat_ewd[3],
             general_data->cell.hmat_ewd[6],general_data->cell.hmat_ewd[9]);
      fprintf(fp_dname,"1/3 log(Vol)\n");
      fprintf(fp_dname,"%.13g\n",general_data->baro.x_lnv);
/*==========================================================================*/
/*    C)Atm and Atm NHC Velocities                                   */

      fprintf(fp_dname,"atm vel\n");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=class->clatoms_info.natm_tot;i++){
       fprintf(fp_dname,"%.13g %.13g %.13g\n",class->clatoms_pos[ip].vx[i],
              class->clatoms_pos[ip].vy[i],class->clatoms_pos[ip].vz[i]);
       }/*endfor*/
      }/*endfor*/
      
      fprintf(fp_dname,"number of atm nhc, length of nhc\n");
      fprintf(fp_dname,"%d %d\n",class->therm_info_class.num_nhc,
             class->therm_info_class.len_nhc);
      fprintf(fp_dname,"atm nhc velocities\n");
      for(j=1;j<=(class->therm_info_class.len_nhc);j++){
       for(i=1;i<=(class->therm_info_class.num_nhc);i++){
         fprintf(fp_dname,"%.13g\n",class->therm_class.v_nhc[j][i]);
       } /*endfor*/
      }/*endfor*/

/*==========================================================================*/
/*    D)Vol and Vol NHC Velocities                             */

      fprintf(fp_dname,"vol velocities\n");
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->par_rahman.vgmat[1],
             general_data->par_rahman.vgmat[4],general_data->par_rahman.vgmat[7]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->par_rahman.vgmat[2],
             general_data->par_rahman.vgmat[5],general_data->par_rahman.vgmat[8]);
      fprintf(fp_dname,"%.13g %.13g %.13g\n",general_data->par_rahman.vgmat[3],
             general_data->par_rahman.vgmat[6],general_data->par_rahman.vgmat[9]);
      fprintf(fp_dname,"log(vol) velocity\n");
      fprintf(fp_dname,"%.13g\n",general_data->baro.v_lnv);
      fprintf(fp_dname,"vol nhc velocities\n");
      for(i=1;i<=(class->therm_info_class.len_nhc);i++){
       fprintf(fp_dname,"%.13g\n",general_data->baro.v_vol_nhc[i]);
      }/*endfor*/
/*==========================================================================*/
/*    E)Misc                                                    */

      fprintf(fp_dname,"dt=%.13g\n",general_data->timeinfo.dt);
      fprintf(fp_dname,"nfree=%d\n",class->clatoms_info.nfree);
      fprintf(fp_dname,"nve=%d nvt=%d npt_i=%d npt_f=%d nst=%d\n",
             general_data->ensopts.nve,  general_data->ensopts.nvt,  
             general_data->ensopts.npt_i,  general_data->ensopts.npt_f,  
             general_data->ensopts.nst);
      fprintf(fp_dname,"nbond_free=%d nbend_free=%d ntors_free=%d\n",
             bonded->bond_free.num,bonded->bend_free.num,  
             bonded->tors_free.num);
      fprintf(fp_dname,"t_ext=%.13g,pext=%.13g,stens_ext=%.13g\n",
             general_data->statepoint.t_ext,
             general_data->statepoint.pext,
             general_data->statepoint.stens_ext);
/*==========================================================================*/
/*    F)CLose                                                               */

      fflush(fp_dname);
      fclose(fp_dname); 

/*==========================================================================*/
}/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_free_energy_file(CLASS *class,BONDED *bonded,
                            GENERAL_DATA *general_data)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

  int i,j,k,m,iii;
  FILE *fp_bond_free,*fp_bend_free,*fp_tors_free,*fp_rbar_free;

/*==========================================================================*/
/*    A) Bonds                                                       */

   if(bonded->bond_free.num >0 ){
       fp_bond_free  = cfopen(bonded->bond_free.file,"o");
       fprintf(fp_bond_free,"Itime\n");
       fprintf(fp_bond_free,"%d\n",general_data->timeinfo.itime);
       fprintf(fp_bond_free,"fk, eq, npow\n");
       fprintf(fp_bond_free,"%.13g %.13g %d\n",bonded->bond_free.fk,
              bonded->bond_free.eq,bonded->bond_free.npow);
       fprintf(fp_bond_free,"nhist,rmin,rmax,dr\n");
       fprintf(fp_bond_free,"%d %.13g %.13g %.13g\n",
              bonded->bond_free.nhist,bonded->bond_free.rmin,
              bonded->bond_free.rmax,bonded->bond_free.del);
       fprintf(fp_bond_free,"histogram\n");
       for(i=1;i<=(bonded->bond_free.nhist);i++){
         fprintf(fp_bond_free,"%.13g\n",bonded->bond_free.hist[i]);
       }/*endfor*/
        fflush(fp_bond_free); 
        fclose(fp_bond_free); 
   }/*endif*/

/*==========================================================================*/
/*    B) Bends                                                  */

   if(bonded->bend_free.num >0 ){
       fp_bend_free  = cfopen(bonded->bend_free.file,"o");
       fprintf(fp_bend_free,"Itime\n");
       fprintf(fp_bend_free,"%d\n",general_data->timeinfo.itime);
       fprintf(fp_bend_free,"fk, eq, npow\n");
       fprintf(fp_bend_free,"%.13g %.13g %d\n",bonded->bend_free.fk,
              bonded->bend_free.eq,bonded->bend_free.npow);
       fprintf(fp_bend_free,"nhist, dtheta\n");
       fprintf(fp_bend_free,"%d %.13g \n",
              bonded->bend_free.nhist,bonded->bend_free.del);
       fprintf(fp_bend_free,"Histogram\n");
       for(i=1;i<=(bonded->bend_free.nhist);i++){
         fprintf(fp_bend_free,"%.13g\n",bonded->bend_free.hist[i]);
       }/*endfor*/
        fflush(fp_bend_free); 
        fclose(fp_bend_free); 
   }/*endif*/

/*==========================================================================*/
/*    C) Tors                                                    */

   if(bonded->tors_free.num == 1){
       fp_tors_free  = cfopen(bonded->tors_free.file,"o");
       fprintf(fp_tors_free,"Itime\n");
       fprintf(fp_tors_free,"%d\n",general_data->timeinfo.itime);
       fprintf(fp_tors_free,"fk, eq, npow\n");
       fprintf(fp_tors_free,"%.13g %.13g %d\n",bonded->tors_free.fk,
              bonded->tors_free.eq[1],bonded->tors_free.npow);
       fprintf(fp_tors_free,"nhist, dtheta\n");
       fprintf(fp_tors_free,"%d -180.0 180.0 %.13g \n",
              bonded->tors_free.nhist,bonded->tors_free.del);
       fprintf(fp_tors_free,"Histogram\n");
       for(i=1;i<=(bonded->tors_free.nhist);i++){
         fprintf(fp_tors_free,"%.13g\n",bonded->tors_free.hist[i]);
       }/*endfor*/
        fflush(fp_tors_free); 
        fclose(fp_tors_free); 
   }/*endif*/

   if(bonded->tors_free.num == 2){
       fp_tors_free  = cfopen(bonded->tors_free.file,"o");
       fprintf(fp_tors_free,"Itime\n");
       fprintf(fp_tors_free,"%d\n",general_data->timeinfo.itime);
       fprintf(fp_tors_free,"fk, eq, fk eq npow\n");
       fprintf(fp_tors_free,"%.13g %.13g %.13g %.13g %d\n",
           bonded->tors_free.fk,bonded->tors_free.eq[1],
           bonded->tors_free.fk,bonded->tors_free.eq[2],
           bonded->tors_free.npow);
       fprintf(fp_tors_free,"nhist, dtheta\n");
       fprintf(fp_tors_free,"%d -180.0 180.0 %.13g \n",
              bonded->tors_free.nhist,bonded->tors_free.del);
       fprintf(fp_tors_free,"nhist, dtheta\n");
       fprintf(fp_tors_free,"%d -180.0 180.0 %.13g \n",
              bonded->tors_free.nhist,bonded->tors_free.del);
       fprintf(fp_tors_free,"Histogram\n");
       for(i=1;i<=(bonded->tors_free.nhist);i++){
       for(j=1;j<=(bonded->tors_free.nhist);j++){
        fprintf(fp_tors_free,"%d %d %.13g\n",j,i,
               bonded->tors_free.hist_2d[j][i]);
       }/*endfor*/
       }/*endfor*/
       fflush(fp_tors_free); 
       fclose(fp_tors_free); 
   }/*endif*/

/*==========================================================================*/
/*    D) rbar sigma                                                         */

   if(bonded->rbar_sig_free.nfree > 0){

        fp_rbar_free  = cfopen(bonded->rbar_sig_free.file,"o");

        fprintf(fp_rbar_free,"Itime\n");
        fprintf(fp_rbar_free,"%d\n",general_data->timeinfo.itime);
        fprintf(fp_rbar_free,"fk_bar, eq_bar, fk_sig, eq_sig nfree\n");
        fprintf(fp_rbar_free,"%.13g %.13g %.13g %.13g %d\n",
                             bonded->rbar_sig_free.fk_bar,
                             bonded->rbar_sig_free.eq_bar,
                             bonded->rbar_sig_free.fk_sigma,
                             bonded->rbar_sig_free.eq_sigma,
                             bonded->rbar_sig_free.nfree);

        fprintf(fp_rbar_free,"nhist_bar, rmin, rmax \n");
        fprintf(fp_rbar_free,"%d %.13g %.13g \n",
                             bonded->rbar_sig_free.nhist_bar,
                             bonded->rbar_sig_free.rmin,
                             bonded->rbar_sig_free.rmax);

        fprintf(fp_rbar_free,"nhist_sig, smin, smax \n");
        fprintf(fp_rbar_free,"%d %.13g %.13g \n",
                             bonded->rbar_sig_free.nhist_sig,
                             bonded->rbar_sig_free.smin,
                             bonded->rbar_sig_free.smax);

        fprintf(fp_rbar_free,"Histogram:ibar,isig,hist\n");
        for(i=1;i<=(bonded->rbar_sig_free.nhist_sig);i++){
         for(j=1;j<=(bonded->rbar_sig_free.nhist_bar);j++){
          fprintf(fp_rbar_free,"%d %d %.13g\n",
                       j,i,bonded->rbar_sig_free.hist[j][i]);
         }/*endfor*/
        }/*endfor*/
      
        for(m=1;m<=bonded->rbar_sig_free.nfree;m++){
         fprintf(fp_rbar_free,"Histogram %d:ibin,ibar,isig,hist\n",m);
         for(k=1;k<=bonded->rbar_sig_free.nhist_bar;k++){
               fprintf(fp_rbar_free,"%d %.13g\n",
                       k,bonded->rbar_sig_free.hist_rn[m][k]);
         }/*endfor*/
        }/*endfor*/

        fflush(fp_rbar_free); 
        fclose(fp_rbar_free); 

   }/*endif*/

/*==========================================================================*/
   }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_config_files_md(CLASS *class,BONDED *bonded,
                           GENERAL_DATA *general_data)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

  int i,ip,n,iii;
  int pi_beads = class->clatoms_info.pi_beads;
  FILE *fp_cpname, *fp_cvname, *fp_cfname;

/*=====================================================================*/
/* I) Write to the atm config file                          */

  if((general_data->timeinfo.itime % general_data->filenames.iwrite_confp)==0){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cpname  = cfopen(general_data->filenames.cpname,"a");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
       fprintf(fp_cpname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[ip].x[i],
              class->clatoms_pos[ip].y[i],class->clatoms_pos[ip].z[i]);
       }/*endfor*/
      }/*endfor*/
      for(i=0;i<3;i++) 
       fprintf(fp_cpname,"%.13g %.13g %.13g\n",general_data->cell.hmat[1+i],
              general_data->cell.hmat[(4+i)],general_data->cell.hmat[(7+i)]);
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
    for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
      fclose(fp_cpname);
   }/*endif*/
   }/*endif*/

/*===================================================================*/
/* II) Write to the partial atm config files                    */

if((general_data->timeinfo.itime%general_data->filenames.iwrite_par_confp)==0){

 if((general_data->filenames.iwrite_conf_binary==0)&&
  (general_data->filenames.low_lim_par<=general_data->filenames.high_lim_par)){
    fp_cpname  = cfopen(general_data->filenames.cpparname,"a");
    for(ip=1;ip<=pi_beads;ip++){
     for(i=(general_data->filenames.low_lim_par);
       i<=(general_data->filenames.high_lim_par);i++){
       fprintf(fp_cpname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[ip].x[i],
             class->clatoms_pos[ip].y[i],class->clatoms_pos[ip].z[i]);
     }/*endfor*/
    }/*endfor*/
    for(i=0;i<3;i++) 
      fprintf(fp_cpname,"%.13g %.13g %.13g\n",general_data->cell.hmat[1+i],
             general_data->cell.hmat[(4+i)],general_data->cell.hmat[(7+i)]);
    fflush(fp_cpname);
    fclose(fp_cpname);
   }/*endif*/

 if((general_data->filenames.iwrite_conf_binary==1)&&
  (general_data->filenames.low_lim_par<=general_data->filenames.high_lim_par)){
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
    for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cpname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cpname);
    }/*endfor*/ 
      fclose(fp_cpname);
   }/*endif*/

}/*endif*/

/*=====================================================================*/
/*=====================================================================*/
 /* III) Write to the atm velocity config file                          */

  if((general_data->timeinfo.itime % general_data->filenames.iwrite_confv)==0){

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
     for(i=0;i<3;i++){
       fprintf(fp_cvname,"%.13g %.13g %.13g\n",general_data->cell.hmat[1+i],
              general_data->cell.hmat[(4+i)],general_data->cell.hmat[(7+i)]);
     }/*endfor*/
     fflush(fp_cvname);
     fclose(fp_cvname);
   }/*endif:ascii*/

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
     for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cvname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cvname);
     }/*endfor*/ 
     fclose(fp_cvname);
   }/*endif:binary*/

  }/*endif:vel write*/

/*=====================================================================*/
/* IV) Write to the atm force file                          */

  if((general_data->timeinfo.itime % general_data->filenames.iwrite_atm_for)==0){
   if(general_data->filenames.iwrite_conf_binary==0){
      fp_cfname  = cfopen(general_data->filenames.forcename,"a");
      for(ip=1;ip<=pi_beads;ip++){
       for(i=1;i<=(class->clatoms_info.natm_tot);i++){
       fprintf(fp_cfname,"%.12g  %.12g  %.12g\n",class->clatoms_pos[ip].fx[i],
              class->clatoms_pos[ip].fy[i],class->clatoms_pos[ip].fz[i]);
       }/*endfor*/
      }/*endfor*/
      for(i=0;i<3;i++) 
        fprintf(fp_cfname,"%.13g %.13g %.13g\n",general_data->cell.hmat[1+i],
               general_data->cell.hmat[(4+i)],general_data->cell.hmat[(7+i)]);
      fflush(fp_cfname);
      fclose(fp_cfname);
    }/*endif*/

   if(general_data->filenames.iwrite_conf_binary==1){
    fp_cfname = cfopen(general_data->filenames.forcename,"a");
    n=1;
   for(ip=1;ip<=pi_beads;ip++){
    for(i=1;i<=(class->clatoms_info.natm_tot);i++){
     fwrite(&(class->clatoms_pos[ip].fx)[i],sizeof(double),n,fp_cfname);
     fwrite(&(class->clatoms_pos[ip].fy)[i],sizeof(double),n,fp_cfname);
     fwrite(&(class->clatoms_pos[ip].fz)[i],sizeof(double),n,fp_cfname);
    }/*endfor*/
   }/*endfor*/
    for(i=0;i<3;i++){ 
      fwrite(&(general_data->cell.hmat)[1+i],sizeof(double),n,fp_cfname);
      fwrite(&(general_data->cell.hmat)[4+i],sizeof(double),n,fp_cfname);
      fwrite(&(general_data->cell.hmat)[7+i],sizeof(double),n,fp_cfname);
    }/*endfor*/ 
      fclose(fp_cfname);
   }/*endif*/
   }/*endif*/

  
/*==========================================================================*/
  }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void write_inst_file_md(CLASS *class,GENERAL_DATA *general_data,
                        double etot,double a,double b,double c,
                        double tac,double tab,double tbc)

/*==========================================================================*/
{/* begin routine */
/*==========================================================================*/

  int i,iii;
  double inst_div;
  FILE *fp_iname;

/*==========================================================================*/
/*    A) Finish Averages                                            */

      inst_div    = (double)(general_data->filenames.iwrite_inst);    
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
      general_data->stat_avg.aipress_inter  /= (inst_div*PCONV);
      general_data->stat_avg.aipress_intra  /= (inst_div*PCONV);
      for(i=1;i<=9;i++){general_data->stat_avg.aipten[i] /= (inst_div*PCONV);}

/*==========================================================================*/
/*   B) Write Averages                                                */

      fp_iname = cfopen(general_data->filenames.iname,"a");
      fprintf(fp_iname,"\n");
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g\n",
           general_data->stat_avg.aikinet,general_data->stat_avg.aikinet_v,
           general_data->stat_avg.aikinet_nhc,general_data->stat_avg.aivintert,
           general_data->stat_avg.aivintrat,general_data->stat_avg.aivol,etot,
           general_data->stat_avg.aipress_inter,
           general_data->stat_avg.aipress_intra);
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",
             general_data->stat_avg.aicella,general_data->stat_avg.aicellb,
             general_data->stat_avg.aicellc,general_data->stat_avg.aicellac,
             general_data->stat_avg.aicellab,general_data->stat_avg.aicellbc);
      for(i=1;i<=9;i+=3){
       fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.aipten[i],
              general_data->stat_avg.aipten[(i+1)],
              general_data->stat_avg.aipten[(i+2)]);}
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",
             general_data->stat_avg.kinet,general_data->stat_avg.kinet_v,
             general_data->stat_avg.kinet_nhc,general_data->stat_avg.vintert,
             general_data->stat_avg.vintrat,general_data->stat_avg.vol);
      fprintf(fp_iname,"%.9g %.9g %.9g %.9g %.9g %.9g\n",a,b,c,tac,tab,tbc);
      for(i=1;i<=9;i+=3) {
       fprintf(fp_iname,"%.9g %.9g %.9g\n",general_data->stat_avg.apten_out[i],
              general_data->stat_avg.apten_out[(i+1)],
              general_data->stat_avg.apten_out[(i+2)]);}
      fflush(fp_iname);
      fclose(fp_iname);

/*==========================================================================*/
/*   C) Zero Averages                                               */

      general_data->stat_avg.aikinet     = 0.0;
      general_data->stat_avg.aikinet_v   = 0.0;
      general_data->stat_avg.aikinet_nhc = 0.0;
      general_data->stat_avg.aivintert   = 0.0;
      general_data->stat_avg.aivintrat   = 0.0;
      general_data->stat_avg.aivol       = 0.0;
      general_data->stat_avg.aicella     = 0.0;
      general_data->stat_avg.aicellb     = 0.0;
      general_data->stat_avg.aicellc     = 0.0;
      general_data->stat_avg.aicellab    = 0.0;
      general_data->stat_avg.aicellbc    = 0.0;
      general_data->stat_avg.aicellac    = 0.0;
      for(i=1;i<=9;i++){general_data->stat_avg.aipten[i]  = 0.0;}

/*==========================================================================*/
  }/* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void forc_level_vel_gather(CLASS *class)

/*==========================================================================*/
  {/* begin routine */
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"
 
  int i;

  double *vx          = class->clatoms_pos[1].vx;
  double *vy          = class->clatoms_pos[1].vy;
  double *vz          = class->clatoms_pos[1].vz;
  double *vxtemp      = class->ewd_scr.x;
  double *vytemp      = class->ewd_scr.y;
  double *vztemp      = class->ewd_scr.z;
  int *recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
  int *displs_atm     = class->class_comm_forc_pkg.displs_atm;
  int myid_forc       = class->class_comm_forc_pkg.myid;
  int myatm_start     = class->clatoms_info.myatm_start;
  int natm_tot        = class->clatoms_info.natm_tot;
  MPI_Comm comm_forc  = class->class_comm_forc_pkg.comm;

/*==========================================================================*/
/* Allgather the velocities */

  Allgatherv(&(vx[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vxtemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
  Allgatherv(&(vy[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vytemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);
  Allgatherv(&(vz[myatm_start]),recv_count_atm[(myid_forc+1)],
             MPI_DOUBLE,&(vztemp[1]),&recv_count_atm[1],&displs_atm[1],
             MPI_DOUBLE,0,comm_forc);

  for(i=1;i<=natm_tot;i++){
       vx[i] = vxtemp[i];
       vy[i] = vytemp[i];
       vz[i] = vztemp[i];
  }/*endfor*/

/*==========================================================================*/
  }/* end routine */
/*==========================================================================*/
