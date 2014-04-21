/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_sim_parms.c                              */
/*                                                                          */
/* This subprogram sets the simulation inputs                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"

#define DEBUG

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_list(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded,
                         CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index,ilimit;
  double real_key_arg;

/*========================================================================*/
  /* 1)\verlist_pad{#} */
       sscanf(dict[1].keyarg,"%lg",&real_key_arg);
       class->nbr_list.verlist.jver_pad = (int)(real_key_arg);
       ilimit = 50;
       if(class->communicate.np_forc>1){
       ilimit = 150;
       if(class->nbr_list.verlist.jver_pad==30){
          class->nbr_list.verlist.jver_pad = 100;
          }
        }
       index=1;
       if(class->nbr_list.verlist.jver_pad <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if((class->nbr_list.verlist.jver_pad < 10)||
                    (class->nbr_list.verlist.jver_pad > ilimit) ){ 
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a verlist_padding < 10 or > %d \n",ilimit);
       printf("Are you certain this is what you would like to do?   \n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\verlist_mem_safe{#} */
       sscanf(dict[2].keyarg,"%lg",&real_key_arg);
       class->nbr_list.verlist.mem_safe = (real_key_arg); 
      if(class->communicate.np_forc>1&&class->nbr_list.verlist.mem_safe==1.25){
       class->nbr_list.verlist.mem_safe = 1.5;
       }
       index=2;
       if(class->nbr_list.verlist.mem_safe <= 0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->nbr_list.verlist.mem_safe > 1.5){ 
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a verlist_mem_safe > 1.5      \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\verlist_mem_min{#} */
       sscanf(dict[3].keyarg,"%lg",&real_key_arg);
       class->nbr_list.verlist.nmem_min_lst = (int)(real_key_arg);
       index=3;
       if(class->nbr_list.verlist.nmem_min_lst <= 0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4)\verlist_skin{#} */
       sscanf(dict[4].keyarg,"%lg",&(class->interact.skin));
       index=4;
       if(class->interact.skin < 0.){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       }
       if(class->interact.skin>2.5){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a verlist skin > 2.5A         \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       class->interact.skin /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 8)\lnkcell_cell_divs{#} */
       sscanf(dict[8].keyarg,"%lg",&real_key_arg);
       class->nbr_list.lnklist.ncell_div_avg = (int)(real_key_arg);
       index=8;
       if(class->nbr_list.lnklist.ncell_div_avg <= 0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->nbr_list.lnklist.ncell_div_avg>40){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a lnk cell with > 40s per side \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\lnk_cell_excl_safe{#} */
       sscanf(dict[9].keyarg,"%lg",&real_key_arg);
       class->nbr_list.lnklist.lnk_excl_safe = (int)(real_key_arg);
       index=9;
       if(class->nbr_list.lnklist.lnk_excl_safe <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 10)\lnk_cell_vol_safe{#} */
       sscanf(dict[10].keyarg,"%lg",&real_key_arg);
       class->nbr_list.lnklist.lnk_vol_safe = (int)(real_key_arg);
       index=10;
       if(class->nbr_list.lnklist.lnk_vol_safe <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->nbr_list.lnklist.lnk_vol_safe > 1.1){ 
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a lnk_cell_vol_safe > 1.1  \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\lnkcell_force_odd{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[11].keyarg,"on")==0) {
       class->nbr_list.lnklist.lnk_for_odd = 1; ifound++;}
       if(strcasecmp(dict[11].keyarg,"off")==0){
          class->nbr_list.lnklist.lnk_for_odd = 0; ifound++;}
       index=11;
    if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 12)\shave_skin_opt{#} */
       ifound      = 0; 
       if(strcasecmp(dict[12].keyarg,"on")==0) {
       class->interact.ishave_opt = 1; ifound++;}
       if(strcasecmp(dict[12].keyarg,"off")==0){
          class->interact.ishave_opt = 0; ifound++;}
       index=12;
    if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /*  13)\neighbor_list{ver_list,lnk_list,no_list} */
       ifound = 0;
       class->nbr_list.iver = 0;
       class->nbr_list.ilnk = 0;
       class->nbr_list.nolst = 0;
       if(strcasecmp(dict[13].keyarg,"no_list")==0)  {
          class->nbr_list.nolst = 1;ifound++;}
       if(strcasecmp(dict[13].keyarg,"lnk_list")==0) {
          class->nbr_list.ilnk = 1;ifound++;}
       if(strcasecmp(dict[13].keyarg,"ver_list")==0) {
          class->nbr_list.iver = 1;ifound++;}
       index=13;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  14)\update_type{no_list,lnk_list} */
       ifound  = 0;
       if(strcasecmp(dict[14].keyarg,"no_list")==0){
          class->nbr_list.verlist.nolst_ver_update = 1;
          class->nbr_list.verlist.lnk_ver_update = 0;
          ifound++;
          }
       if(strcasecmp(dict[14].keyarg,"lnk_list")==0) {
          class->nbr_list.verlist.nolst_ver_update = 0;
          class->nbr_list.verlist.lnk_ver_update = 1;
          ifound++;
          }
       index=14;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/ 
  /* 5)\brnch_root_list_opt{#} */
       ifound      = 0;
       if(strcasecmp(dict[5].keyarg,"on")==0) {
          class->nbr_list.brnch_root_list_opt = 1; ifound++;}
       if(strcasecmp(dict[5].keyarg,"off")==0){
          class->nbr_list.brnch_root_list_opt = 0; ifound++;}
       if(strcasecmp(dict[5].keyarg,"pare_down")==0){
          class->nbr_list.brnch_root_list_opt = 2; ifound++;}
       index=5;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  if(class->nbr_list.iver ==0 && class->nbr_list.brnch_root_list_opt > 0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Brnch-root list option only implemented for verlet lists.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\brnch_root_list_skin{#} */
       sscanf(dict[6].keyarg,"%lg",&(class->interact.brnch_root_skin));
       index=6;
       if(class->interact.brnch_root_skin < 0.){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       }
       class->interact.brnch_root_skin /= BOHR;
       if(class->nbr_list.brnch_root_list_opt==0){
       class->interact.brnch_root_skin = 0.0;
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\brnch_root_cutoff{#} */
       sscanf(dict[7].keyarg,"%lg",&(class->interact.brnch_root_cut));
       index=7;
       if(class->interact.brnch_root_cut < 0.){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       }
       if(class->interact.brnch_root_cut > 1.112){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a brnch root cutoff > 1.112A \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       class->interact.brnch_root_cut /= BOHR;

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_cp(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                       CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                       FILENAME_PARSE *filename_parse,
                       DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */
int iii;

  int ifound,index,itemp,cp_orth_sum;
  double real_key_arg;

/*========================================================================*/

    cp->cpcoeffs_info.ncoef_l   = 0;
    cp->cpcoeffs_info.ncoef     = 0;
    cp->cpcoeffs_info.nstate_up = 0;
    cp->cpcoeffs_info.nstate_dn = 0;
    cp->cptherm_info.num_c_nhc  = 0;
    cp->cptherm_info.num_c_nhc_proc  = 0;
    cp->pseudo.n_ang_max   = 0;

/*========================================================================*/
  /*  8)\cp_dft_typ{lsda,lda,gga_lda,gga_lsda} */
        ifound  = 0;
        cp->cpopts.cp_lda  = 0;
        cp->cpopts.cp_lsda = 0;
        cp->cpopts.cp_gga  = 0;
        if(strcasecmp(dict[8].keyarg,"lsda")==0) {
        cp->cpopts.cp_lsda = 1;ifound++;}
        if(strcasecmp(dict[8].keyarg,"lda")==0)  {
        cp->cpopts.cp_lda  = 1;ifound++;}
        if(strcasecmp(dict[8].keyarg,"gga_lsda")==0){
        cp->cpopts.cp_lsda = 1;cp->cpopts.cp_gga=1; ifound++;}
        if(strcasecmp(dict[8].keyarg,"gga_lda")==0) {
        cp->cpopts.cp_lda  = 1;cp->cpopts.cp_gga=1; ifound++;}
        index=8;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 1)\cp_vxc_typ{pz_lda,pz_lsda,pw_lda,pw_lsda,pade_lda,pade_lsda}       */

       ifound = 0;
       index=1;
       if(strcasecmp(dict[1].keyarg,"pz_lda")==0) {ifound++;}
       if(strcasecmp(dict[1].keyarg,"pz_lsda")==0){ifound++;}
       if(strcasecmp(dict[1].keyarg,"pw_lda")==0) {ifound++;}
       if(strcasecmp(dict[1].keyarg,"pw_lsda")==0) {ifound++;}
       if(strcasecmp(dict[1].keyarg,"pade_lda")==0){ifound++;}
       if(strcasecmp(dict[1].keyarg,"pade_lsda")==0){ifound++;}
 
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       sscanf(dict[1].keyarg,"%s",cp->pseudo.vxc_typ);

  /*-----------------------------------------------------------------------*/ 
  /* 2)\cp_ggax_typ{becke,pw91x,debug97x,fila_1x,fila_2x,pbe_x,revpbe_x,rpbe_x,
       xpbe_x,brx89,brx2k,off}  */
       ifound = 0;
       itemp = cp->cpopts.cp_gga;
       cp->cpopts.cp_gga = 0;
       cp->cpcoeffs_info.cp_laplacian_on = 0;
       if(strcasecmp(dict[2].keyarg,"off")==0)    {
           ifound++;cp->cpopts.cp_gga=0;}
       if(strcasecmp(dict[2].keyarg,"becke")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"pw91x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"fila_1x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"fila_2x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"pbe_x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"revpbe_x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"rpbe_x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"xpbe_x")==0)  {
           ifound++;cp->cpopts.cp_gga=1;}
       if(strcasecmp(dict[2].keyarg,"brx89")==0 ||
          strcasecmp(dict[2].keyarg,"brx2k")==0) {
            ifound++;
            cp->cpopts.cp_gga=1;
            cp->cpcoeffs_info.cp_ke_dens_on = 1;
            cp->cpcoeffs_info.cp_tau_functional = 1;
            cp->cpcoeffs_info.cp_laplacian_on   = 1;
       }/* endif brx89 (a tau functional) */
       if(strcasecmp(dict[2].keyarg,"debug97x")==0)  {
        if(general_data->simopts.debug_cp != 1 || 
                           general_data->simopts.debug_cp_pimd != 1){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
         printf("You may not choose the debug97x functional unless you \n");
         printf("actually doing a debug_cp or debug_cp_pimd\n");
         printf("simulation. It isn't a real functional\n");
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
        }/* endif */
        ifound++;cp->cpopts.cp_gga=1;
       }/* endif */
       index=2;
       if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       sscanf(dict[2].keyarg,"%s",cp->pseudo.ggax_typ);
#ifdef JUNK
       if(itemp!=cp->cpopts.cp_gga)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
#endif

  /*-----------------------------------------------------------------------*/ 
  /* 3)\cp_ggac_typ{lyp,lypm1,pw91c,pbe_c,xpbe_c,tau1_c,off}               */
      ifound = 0;
      itemp = cp->cpopts.cp_gga;
      cp->cpopts.cp_gga = 0;
      index=3;
      if(strcasecmp(dict[3].keyarg,"off")==0)    {
         ifound++;cp->cpopts.cp_gga=1;}
      if(strcasecmp(dict[3].keyarg,"pw91c")==0)  {
         ifound++;cp->cpopts.cp_gga=1;}
      if(strcasecmp(dict[3].keyarg,"lyp")==0)    {
         ifound++;cp->cpopts.cp_gga=1;
#ifdef JUNK
         if(strcasecmp(cp->pseudo.vxc_typ,"pade_lda")==0){ 
           printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"); 
           printf("You may not choose the LYP functional with the Pade \n");
           printf("approximation to Perdew-Wang 1992 LDA functional\n");
           printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
         }/* endif */
         if(strcasecmp(cp->pseudo.vxc_typ,"pade_lsda")==0){ 
           printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"); 
           printf("You may not choose the LYP functional with the Pade \n");
           printf("approximation to Perdew-Wang 1992 LSDA functional\n");
           printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
         }/* endif */
#endif
      }/* endif */
      if(strcasecmp(dict[3].keyarg,"pbe_c")==0)    {
	ifound++;cp->cpopts.cp_gga=1;}
      if(strcasecmp(dict[3].keyarg,"xpbe_c")==0)   {
	ifound++;cp->cpopts.cp_gga=1;}
      if(strcasecmp(dict[3].keyarg,"lypm1")==0)    {
	ifound++;cp->cpopts.cp_gga=1;cp->cpcoeffs_info.cp_laplacian_on=1;}
      if(strcasecmp(dict[3].keyarg,"tau1_c")==0 ){
            ifound++;
            cp->cpopts.cp_gga=1;
            cp->cpcoeffs_info.cp_ke_dens_on     = 1;
            cp->cpcoeffs_info.cp_tau_functional = 1;
            cp->cpcoeffs_info.cp_laplacian_on   = 1;
       }/* endif tau1 (a tau- and Laplacian-dependent  functional) */
      if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      sscanf(dict[3].keyarg,"%s",cp->pseudo.ggac_typ);
      if(strcasecmp(dict[3].keyarg,"off")==0 && 
	 strcasecmp(dict[2].keyarg,"off")==0)    {
	cp->cpopts.cp_gga=0;}
      if(itemp!=cp->cpopts.cp_gga)
	keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/  
  /*  4)\cp_sic{on,off} */
      ifound  = 0;
      if(strcasecmp(dict[4].keyarg,"on")==0)  {cp->cpopts.cp_sic = 1;ifound++;}
      if(strcasecmp(dict[4].keyarg,"off")==0) {cp->cpopts.cp_sic = 0;ifound++;}
      index=4;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  5)\cp_e_e_interact{on,off} */
      ifound      = 0;
    if(strcasecmp(dict[5].keyarg,"on")==0) {cp->cpopts.cp_nonint=0;ifound++;}
    if(strcasecmp(dict[5].keyarg,"off")==0){cp->cpopts.cp_nonint=1;ifound++;}
      index=5;
    if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /*  6)\cp_norb{full_ortho,norm_only,no_constrnt,off} */
      ifound      = 0;
      if(strcasecmp(dict[6].keyarg,"full_ortho")==0) 
                             {cp->cpopts.cp_norb = 1;ifound++;}
      if(strcasecmp(dict[6].keyarg,"norm_only")==0) 
                             {cp->cpopts.cp_norb = 2;ifound++;}
    if(strcasecmp(dict[6].keyarg,"no_constrnt")==0) {
       cp->cpopts.cp_norb = 3;
       ifound++;
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Without constraints on the non-orthgonal orbitals\n");
       printf("they will probably develop large dependencies.\n");
       printf("Use at your own risk!");
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
    }
    if(strcasecmp(dict[6].keyarg,"off")==0){cp->cpopts.cp_norb = 0;ifound++;}
       index=6;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /*  7)\cp_gauss{on,off} */
        ifound      = 0;
     if(strcasecmp(dict[7].keyarg,"on")==0){cp->cpopts.cp_gauss=1;ifound++;}
     if(strcasecmp(dict[7].keyarg,"off")==0){cp->cpopts.cp_gauss=0;ifound++;}
        index=7;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 9)\cp_nl_list{on,off} */
        ifound      = 0;
        if(strcasecmp(dict[9].keyarg,"on")==0) {
        cp->pseudo.nl_cut_on = 1;ifound++;}
        if(strcasecmp(dict[9].keyarg,"off")==0){
        cp->pseudo.nl_cut_on = 0;ifound++;}
        index=9;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 10)\cp_mass_tau_def{#} */
        sscanf(dict[10].keyarg,"%lg",&(cp_parse->cp_mass_tau_def));
        index=10;
        if(cp_parse->cp_mass_tau_def<=0.)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_mass_cut_def{#} */
        sscanf(dict[11].keyarg,"%lg",&(cp_parse->cp_mass_cut_def));
        index=11;
        if(cp_parse->cp_mass_cut_def<=0.)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_energy_cut_def{#} */
        sscanf(dict[12].keyarg,"%lg",&(cp_parse->cp_ecut_def));
        index=12;
        if(cp_parse->cp_ecut_def <= 0. ) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\cp_fict_KE{#} */
        sscanf(dict[13].keyarg,"%lg",&(cp->cpopts.te_ext));
        index=13;
        if(cp->cpopts.te_ext <= 0. ) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        if((cp->cpopts.te_ext > 3000.0) || (cp->cpopts.te_ext <100.0)) {
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("The recommended fictious cp_fict_KE range is 100K - 3000K\n"); 
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /*  14)\cp_ptens{on,off} */
        ifound      = 0;
        if(strcasecmp(dict[14].keyarg,"on")==0)
                          {cp->cpopts.cp_ptens_calc=1;ifound++;}
        if(strcasecmp(dict[14].keyarg,"off")==0)
                          {cp->cpopts.cp_ptens_calc=0;ifound++;}
        index=14;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 15)\cp_init_orthog{on,off} */
        ifound  = 0;
        if(strcasecmp(dict[15].keyarg,"on")==0)  
           {cp->cpopts.cp_init_orthog = 1;ifound++;}
        if(strcasecmp(dict[15].keyarg,"off")==0) 
           {cp->cpopts.cp_init_orthog = 0;ifound++;}
        index=15;
  if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 16)\cp_cg_line_min_len{#} */
        sscanf(dict[16].keyarg,"%lg",&real_key_arg);
        general_data->minopts.cp_cg_line_min_len = (int)(real_key_arg);
        index=16;
        if(general_data->minopts.cp_cg_line_min_len<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_minimize_typ{min_std,min_cg,min_diis} */
        ifound = 0;
        cp->cpopts.cp_hess_calc=0;
        general_data->minopts.cp_min_std = 0;
        general_data->minopts.cp_min_cg = 0;
        general_data->minopts.cp_min_diis = 0;
        if(strcasecmp(dict[17].keyarg,"min_std")==0) {
           general_data->minopts.cp_min_std  = 1; ifound++;}
        if(strcasecmp(dict[17].keyarg,"min_cg")==0) {
           general_data->minopts.cp_min_cg   = 1; 
           cp->cpopts.cp_hess_calc=1; ifound++;}
        if(strcasecmp(dict[17].keyarg,"min_diis")==0){
           general_data->minopts.cp_min_diis = 1;
           cp->cpopts.cp_hess_calc=1; ifound++;}
        index=17;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
    if(general_data->simopts.cp_min == 0 || general_data->simopts.cp_wave_min == 0)
      cp->cpopts.cp_hess_calc = 0;
  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_diis_hist_len{#} */
        sscanf(dict[18].keyarg,"%lg",&real_key_arg);
        general_data->minopts.cp_diis_hist_len = (int)(real_key_arg);
        index=18;
        if(general_data->minopts.cp_diis_hist_len<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 19)\cp_orth_meth{gram_schmidt,lowdin,normalize,none} */
        ifound  = 0;
        cp->cpopts.cp_gs = 0;
        cp->cpopts.cp_low = 0;
        cp->cpopts.cp_normalize = 0;
        if(strcasecmp(dict[19].keyarg,"lowdin")==0) {
            cp->cpopts.cp_low = 1;ifound++;}
        if(strcasecmp(dict[19].keyarg,"gram_schmidt")==0) {
            cp->cpopts.cp_gs=1;ifound++;}
        if(strcasecmp(dict[19].keyarg,"normalize")==0) {
            cp->cpopts.cp_normalize=1;ifound++;}
        if(strcasecmp(dict[19].keyarg,"none")==0) {
            ifound++;}
        index=19;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        if(general_data->simopts.cp_min == 1 || 
              general_data->simopts.cp_wave_min == 1 ||
        general_data->simopts.cp_wave_min_pimd == 1 || 
              cp->cpopts.cp_init_orthog == 1){
        if(cp->cpopts.cp_normalize == 1 && cp->cpopts.cp_norb !=2){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
       printf("You may only use the normalization option under\n");
       printf("norb minimization or md with the norm_only option\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
       }/* endif */
       cp_orth_sum = cp->cpopts.cp_gs + cp->cpopts.cp_low 
               + cp->cpopts.cp_normalize;
       if(cp_orth_sum == 0 && cp->cpopts.cp_norb !=3){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
       printf("You may only use the none option under\n");
       printf("norb minimization or md with the none option\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
       }/* endif */
       }/* endif general_data */
  /*-----------------------------------------------------------------------*/ 
  /* 20)\cp_restart_type{initial,restart_pos,restart_posvel,restart_all}*/
       ifound = 0;
       if(strcasecmp(dict[20].keyarg,"gen_wave")==0)    {
          cp_parse->istart_cp = 0; ifound++;}
       if(strcasecmp(dict[20].keyarg,"initial")==0)    {
          cp_parse->istart_cp = 1; ifound++;}
       if(strcasecmp(dict[20].keyarg,"restart_pos")==0)  {
          cp_parse->istart_cp = 2; ifound++;}
       if(strcasecmp(dict[20].keyarg,"restart_posvel")==0){
          cp_parse->istart_cp= 3; ifound++;}
       if(strcasecmp(dict[20].keyarg,"restart_all")==0)  {
          cp_parse->istart_cp = 4; ifound++;}
       index=20;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 21)\diis_hist_len{#} */
        sscanf(dict[21].keyarg,"%lg",&real_key_arg);
        general_data->minopts.diis_hist_len = (int)(real_key_arg);
        index=21;
        if(general_data->minopts.diis_hist_len<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 22)\nlvps_list_skin{#}   */
        sscanf(dict[22].keyarg,"%lg",&real_key_arg);
        cp->pseudo.nlvps_skin = (real_key_arg);
        index=22;
        if(cp->pseudo.nlvps_skin<0){
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
        cp->pseudo.nlvps_skin /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 23)\gradient_cutoff{#}   */
       sscanf(dict[23].keyarg,"%lg",&real_key_arg);
       /* In order for LDA and LSDA calculations to be consistent with each
          other, gga_cut (lsda) = gga_cut/2.0 (lda) Search for rho_cut */
       if(cp->cpopts.cp_lsda==1){
        real_key_arg /= 2.0;
        sprintf(dict[23].keyarg,"%e",real_key_arg);
       }
       cp->pseudo.gga_cut = (real_key_arg);
       index=23;
       if(cp->pseudo.gga_cut<0){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /* 24)\zero_cp_vel{initial,periodic,no}*/
       ifound = 0;
       if(strcasecmp(dict[24].keyarg,"initial")==0)    {
          cp->cpopts.zero_cp_vel = 1; ifound++;}
       if(strcasecmp(dict[24].keyarg,"periodic")==0)    {
          cp->cpopts.zero_cp_vel = 2; ifound++;}
       if(strcasecmp(dict[24].keyarg,"no")==0)  {
          cp->cpopts.zero_cp_vel = 0; ifound++;}
       index=24;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 25)\cp_check_perd_size{#}   */
       ifound = 0;
       if(strcasecmp(dict[25].keyarg,"on")==0)    {
          cp->cpopts.icheck_perd_size = 1; ifound++;}
       if(strcasecmp(dict[25].keyarg,"off")==0)    {
          cp->cpopts.icheck_perd_size = 0; ifound++;}
       index=25;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 26)\cp_tol_edge_dist{#}   */
       sscanf(dict[26].keyarg,"%lg",&real_key_arg);
       cp->cpopts.tol_edge_dist = (real_key_arg);
       index=26;
       if(cp->cpopts.tol_edge_dist <= 0){
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

      if( (cp->cpopts.tol_edge_dist < 3.0) 
        || (cp->cpopts.tol_edge_dist > 6.0) ) {
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");    
        printf("The recommended value of the CP cluster size tolerance\n");
        printf("is 3-6 A\n"); 
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      }/*endif*/
      cp->cpopts.tol_edge_dist /= BOHR;

  /*-----------------------------------------------------------------------*/ 
  /* 27)\cp_para_typ{#}   */
       ifound = 0;
       if(strcasecmp(dict[27].keyarg,"hybrid")==0)    {
          cp->cpopts.cp_para_opt = 0; ifound++;}
       if(strcasecmp(dict[27].keyarg,"full_g")==0)    {
          cp->cpopts.cp_para_opt = 1; ifound++;}
       index=27;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 28)\cp_dual_grid_opt{#} */
       ifound = 0;
       if(strcasecmp(dict[28].keyarg,"off")==0)    {
          cp->cpopts.cp_dual_grid_opt = 0; ifound++;}
       if(strcasecmp(dict[28].keyarg,"prop")==0)    {
          cp->cpopts.cp_dual_grid_opt = 1; ifound++;}
       if(strcasecmp(dict[28].keyarg,"not_prop")==0)    {
          cp->cpopts.cp_dual_grid_opt = 2; ifound++;}
       index=28;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 29)\cp_check_perd_size{#}   */
       ifound = 0;
       if(strcasecmp(dict[29].keyarg,"on")==0)    {
          cp->cpopts.icheck_dual_size = 1; ifound++;}
       if(strcasecmp(dict[29].keyarg,"off")==0)    {
          cp->cpopts.icheck_dual_size = 0; ifound++;}
       index=29;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 30)\cp_move_dual_box_opt{#}   */
       ifound = 0;
       if(strcasecmp(dict[30].keyarg,"on")==0)    {
          general_data->cell.imov_cp_box = 1; ifound++;}
       if(strcasecmp(dict[30].keyarg,"off")==0)    {
          general_data->cell.imov_cp_box = 0; ifound++;}
       index=30;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 31)\cp_energy_cut_dual_grid_def{#} */
        sscanf(dict[31].keyarg,"%lg",&(cp_parse->cp_ecut_dual_grid_def));
        index=31;
        if(cp_parse->cp_ecut_dual_grid_def <= 0. ) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 32)\cp_alpha_conv_dual{#} */
        sscanf(dict[32].keyarg,"%lg",&(cp->pseudo.alpha_conv_dual));
        index=32;
        if(cp->pseudo.alpha_conv_dual <= 0. ) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if( cp->pseudo.alpha_conv_dual < 7.0 ||
           cp->pseudo.alpha_conv_dual > 10.0){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$$$\n");    
        printf("The recommended range of cp_alpha_conv_dual is 7-10 \n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       }
  /*-----------------------------------------------------------------------*/ 
  /* 33)\interp_pme_dual{#} */
       sscanf(dict[33].keyarg,"%lg",&real_key_arg);
       cp->pseudo.n_interp_pme_dual= (int)(real_key_arg);
       index=33;
       if((cp->pseudo.n_interp_pme_dual  <4 )||
       ((cp->pseudo.n_interp_pme_dual  %2)!=0))
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 34)\cp_elf_calc_frq{#} */
       sscanf(dict[34].keyarg,"%lg",&real_key_arg);
       cp->cpcoeffs_info.cp_elf_calc_frq =  (int)(real_key_arg);
       index=34;
       if((cp->cpcoeffs_info.cp_elf_calc_frq  <0 ))
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(cp->cpcoeffs_info.cp_elf_calc_frq > 0) 
         cp->cpcoeffs_info.cp_ke_dens_on = 1;

  /*-----------------------------------------------------------------------*/ 
  /* 35)\cp_ngrid_skip{#} */
       sscanf(dict[35].keyarg,"%lg",&real_key_arg);
       cp->cpopts.cp_ngrid_skip =  (int)(real_key_arg);
       index=35;
       if((cp->cpopts.cp_ngrid_skip  < 0 ))
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 36)\cp_isok_opt{#}   */
       ifound = 0;
       if(strcasecmp(dict[36].keyarg,"on")==0)    {
          cp->cpopts.cp_isok_opt = 1; ifound++;}
       if(strcasecmp(dict[36].keyarg,"off")==0)    {
          cp->cpopts.cp_isok_opt = 0; ifound++;}
       index=36;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 37)\cp_hess_cut{#} */
       sscanf(dict[37].keyarg,"%lg",&real_key_arg);
       cp->cpcoeffs_info.cp_hess_cut =  real_key_arg;
       index=37;
       if((cp->cpcoeffs_info.cp_hess_cut < 0.0 ))
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_gen(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                        CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse,
                        DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  double real_key_arg;
  int nproc_tot;

/*========================================================================*/
/*-----------------------------------------------------------------------*/ 
/*  1)\simulation_typ{md,minimize,pimd,cp,cp_wave,cp_min,cp_wave_min,
                       cp_wave_min_pimd,cp_pimd,cp_wave_pimd,
                       debug,debug_pimd,debug_cp} */

  ifound =0;

  general_data->simopts.minimize = 0;
  general_data->simopts.md = 0;
  general_data->simopts.pimd = 0;
  general_data->simopts.cp = 0;
  general_data->simopts.cp_wave = 0;
  general_data->simopts.cp_wave_pimd = 0;
  general_data->simopts.cp_pimd = 0; 
  general_data->simopts.cp_min = 0; 
  general_data->simopts.cp_wave_min = 0.;
  general_data->simopts.cp_wave_min_pimd = 0.;
  general_data->simopts.debug = 0;
  general_data->simopts.debug_pimd = 0;
  general_data->simopts.debug_cp = 0;
  general_data->simopts.debug_cp_pimd = 0; 

  if(strcasecmp(dict[1].keyarg,"minimize")==0)  {
    general_data->simopts.minimize = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"md")==0){
    general_data->simopts.md = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"pimd")==0){
    general_data->simopts.pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp")==0){
    general_data->simopts.cp = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_pimd")==0){
    general_data->simopts.cp_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave")==0){
    general_data->simopts.cp_wave = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_pimd")==0){
    general_data->simopts.cp_wave_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_min")==0){
    general_data->simopts.cp_min = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_min")==0){
    general_data->simopts.cp_wave_min=1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"cp_wave_min_pimd")==0){
    general_data->simopts.cp_wave_min_pimd=1;ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug")==0) {
    general_data->simopts.debug = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_cp")==0) {
    general_data->simopts.debug_cp = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_cp_pimd")==0) {
    general_data->simopts.debug_cp_pimd = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"debug_pimd")==0) {
    general_data->simopts.debug_pimd = 1; ifound++;}
  index=1;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /*  2)\ensemble_typ{nve,nvt,nvt_isok,npt_i,npt_f,nst} */
        ifound = 0;
        general_data->ensopts.nve      = 0;
        general_data->ensopts.nvt      = 0;
        general_data->ensopts.nvt_isok = 0;
        general_data->ensopts.npt_i    = 0;
        general_data->ensopts.npt_f    = 0;
        general_data->ensopts.nst      = 0;

        if(strcasecmp(dict[2].keyarg,"nve")==0)  {
           general_data->ensopts.nve = 1;ifound++;}
        if(strcasecmp(dict[2].keyarg,"nvt")==0)  {
           general_data->ensopts.nvt = 1;ifound++;}
        if(strcasecmp(dict[2].keyarg,"nvt_isok")==0)  {
           general_data->ensopts.nvt_isok = 1;ifound++;}
        if(strcasecmp(dict[2].keyarg,"npt_i")==0){
           general_data->ensopts.npt_i = 1;ifound++;}
        if(strcasecmp(dict[2].keyarg,"npt_f")==0){
           general_data->ensopts.npt_f = 1;ifound++;}
  
         index=2;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\num_time_step{#} */
        sscanf(dict[3].keyarg,"%lg",&real_key_arg);
        general_data->timeinfo.ntime = (int)(real_key_arg);
        index=3;
        if(general_data->timeinfo.ntime<0)
           keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        if(general_data->timeinfo.ntime>10000000){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("You have requested over 10 million time steps    \n");
        printf("Are you certain this is what you would like to do?\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 4)\time_step{#} */
       sscanf(dict[4].keyarg,"%lg",&real_key_arg);
       general_data->timeinfo.dt = real_key_arg;
       index=4;
       if(general_data->timeinfo.dt==0.0)
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->timeinfo.dt>50.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a time step greater than 50fs\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       general_data->timeinfo.dt /= TIME_CONV;

  /*-----------------------------------------------------------------------*/ 
  /* 5)\temperature{#} */
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       general_data->statepoint.t_ext = real_key_arg;
       index=5;
       if(general_data->statepoint.t_ext<0.0)
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->statepoint.t_ext>1000.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a temperature greater than 1000K\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       class->clatoms_info.pi_temperature = general_data->statepoint.t_ext;

  /*-----------------------------------------------------------------------*/ 
  /* 6)\pressure{#} */
       sscanf(dict[6].keyarg,"%lg",&real_key_arg);
       general_data->statepoint.pext = real_key_arg;
       index=6;
       if(general_data->statepoint.pext<0.0)
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->statepoint.pext>100000.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a pressure greater than 0.1Mbar\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       general_data->statepoint.pext *= PCONV;
  /*-----------------------------------------------------------------------*/  
  /* 7)\restart_type{initial,restart_pos,restart_posvel,restart_all}*/
      ifound = 0;
      if(strcasecmp(dict[7].keyarg,"initial")==0)    {
      class_parse->istart = 1; ifound++;}
      if(strcasecmp(dict[7].keyarg,"restart_pos")==0)  {
         class_parse->istart = 2; ifound++;}
      if(strcasecmp(dict[7].keyarg,"restart_posvel")==0){
         class_parse->istart= 3; ifound++;}
      if(strcasecmp(dict[7].keyarg,"restart_all")==0)  {
         class_parse->istart = 4; ifound++;}
      index=7;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 8)\minimize_typ{min_std,min_cg,min_diis} */
      ifound = 0;
      general_data->minopts.min_std = 0;
      general_data->minopts.min_cg = 0;
      general_data->minopts.min_diis = 0;
      if(strcasecmp(dict[8].keyarg,"min_std")==0) {
         general_data->minopts.min_std  = 1; ifound++;}
      if(strcasecmp(dict[8].keyarg,"min_cg")==0) {
         general_data->minopts.min_cg   = 1; ifound++;}
      if(strcasecmp(dict[8].keyarg,"min_diis")==0){
         general_data->minopts.min_diis = 1; ifound++;}
      index=8;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 9)\annealing_rate{#} */
      sscanf(dict[9].keyarg,"%lg",&real_key_arg);
      general_data->simopts.ann_rate = (double)(real_key_arg);
      index=9;
      if(general_data->simopts.ann_rate<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 10)\num_proc_beads{#} */
      sscanf(dict[10].keyarg,"%lg",&real_key_arg);
      class->communicate.np_beads = (int)(real_key_arg);
      index=10;
      if(class->communicate.np_beads<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 11)\num_proc_states{#} */
      sscanf(dict[11].keyarg,"%lg",&real_key_arg);
      class->communicate.np_states = (int)(real_key_arg);
      index=11;
      if(class->communicate.np_states<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 12)\num_proc_class_forc{#} */
      sscanf(dict[12].keyarg,"%lg",&real_key_arg);
      class->communicate.np_forc= (int)(real_key_arg);
      index=12;
      if(class->communicate.np_forc<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\num_proc_tot{#} */
      sscanf(dict[13].keyarg,"%lg",&real_key_arg);
      nproc_tot = (int)(real_key_arg);
      index=13;
      if(nproc_tot< 1)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class->communicate.np != nproc_tot){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The number of processors expected %d\n",class->communicate.np);
        printf("not equal to the number of processors given %d\n",nproc_tot);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");      
        fflush(stdout);
        exit(1);
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 14)\rndm_seed{#} */
      sscanf(dict[14].keyarg,"%d",&(class->vel_samp_class.iseed));
      sscanf(dict[14].keyarg,"%d",&(cp->vel_samp_cp.iseed));
      class->vel_samp_class.qseed = (double) class->vel_samp_class.iseed;
      cp->vel_samp_cp.qseed = (double) cp->vel_samp_cp.iseed;
  /*-----------------------------------------------------------------------*/ 
  /* 15)\rndm_seed2{#} */
      sscanf(dict[15].keyarg,"%d",&(class->vel_samp_class.iseed2));
      sscanf(dict[15].keyarg,"%d",&(cp->vel_samp_cp.iseed2));
  /*-----------------------------------------------------------------------*/ 
  /* 16)\generic_fft_opt{on,off} */
        ifound = 0;
        if(strcasecmp(dict[16].keyarg,"on")==0)  {
           cp->cp_sclr_fft_pkg3d_sm.igeneric_opt  = 1;
           cp->cp_para_fft_pkg3d_sm.igeneric_opt  = 1;
           cp->cp_sclr_fft_pkg3d_lg.igeneric_opt  = 1;
           cp->cp_para_fft_pkg3d_lg.igeneric_opt  = 1;
           cp->cp_sclr_fft_pkg3d_dens_cp_box.igeneric_opt  = 1;
           cp->cp_para_fft_pkg3d_dens_cp_box.igeneric_opt  = 1;
           general_data->pme_res_fft_pkg.igeneric_opt  = 1;
           general_data->pme_fft_pkg.igeneric_opt = 1;ifound++;
        }
        if(strcasecmp(dict[16].keyarg,"off")==0)  {
           cp->cp_sclr_fft_pkg3d_sm.igeneric_opt  = 0;
           cp->cp_para_fft_pkg3d_sm.igeneric_opt  = 0;
           cp->cp_sclr_fft_pkg3d_lg.igeneric_opt  = 0;
           cp->cp_para_fft_pkg3d_lg.igeneric_opt  = 0;
           cp->cp_sclr_fft_pkg3d_dens_cp_box.igeneric_opt  = 0;
           cp->cp_para_fft_pkg3d_dens_cp_box.igeneric_opt  = 0;
           general_data->pme_res_fft_pkg.igeneric_opt  = 0;
           general_data->pme_fft_pkg.igeneric_opt = 0;ifound++;
        }
        index=16;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 17)\alpha_clus{#}   */
       sscanf(dict[17].keyarg,"%lg",&real_key_arg);
       general_data->ewald.alp_clus = (real_key_arg);
       index=17;
       if(general_data->ewald.alp_clus<0){
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
      if( (general_data->ewald.alp_clus < 7.0) || 
          (general_data->ewald.alp_clus > 22.0) ) {
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("The recommended value of the dimensionless parameter\n");
        printf("alp_clus is 8-22\n"); 
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 18)\ecut_clus{#}   */
       sscanf(dict[18].keyarg,"%lg",&real_key_arg);
       general_data->ewald.ecut_clus = (real_key_arg);
       index=18;
       if(general_data->ewald.ecut_clus < 0){
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
  /* 19)\surf_tens{#} */
       sscanf(dict[19].keyarg,"%lg",&real_key_arg);
       general_data->statepoint.stens_ext= real_key_arg;
       index=19;
       if(general_data->statepoint.stens_ext>1000.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       printf("You have requested a surface tension greater than 1000 N/m\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       general_data->statepoint.stens_ext*= STENS_CONV;

  /*-----------------------------------------------------------------------*/
  /* 20)\min_num_atoms_per_proc{#}   */
       sscanf(dict[20].keyarg,"%lg",&real_key_arg);
       class->clatoms_info.natm_proc = (int)(real_key_arg);
       index=20;
       if(class->clatoms_info.natm_proc < 0){
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

  /*-----------------------------------------------------------------------*/
  /* 21)\num_proc_class_forc_src{#} */
      sscanf(dict[21].keyarg,"%lg",&real_key_arg);
      class->communicate.np_forc_src= (int)(real_key_arg);
      class->class_comm_forc_pkg.plimpton_ind.num_proc_source=
                             class->communicate.np_forc_src;
      index=21;
      if(class->communicate.np_forc_src<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 22)\num_proc_class_forc_trg{#} */
      sscanf(dict[22].keyarg,"%lg",&real_key_arg);
      class->communicate.np_forc_trg= (int)(real_key_arg);
      class->class_comm_forc_pkg.plimpton_ind.num_proc_target=
                             class->communicate.np_forc_trg;
      index=22;
      if(class->communicate.np_forc_trg<=0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);

#ifdef IVO_PLIMP
  /*-----------------------------------------------------------------------*/
      if((class->communicate.np_forc_trg)*(class->communicate.np_forc_src)
         !=class->communicate.np_forc){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("Num_proc_src*Num_proc_trg != Num_proc_forc\n");
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
  /*-----------------------------------------------------------------------*/
#else
      class->communicate.np_forc_src= class->communicate.np_forc;
      class->class_comm_forc_pkg.plimpton_ind.num_proc_source=
                                      class->communicate.np_forc;
#endif

  /*-----------------------------------------------------------------------*/
  /* 23)\annealing_opt{on,off} */

        ifound = 0;
        if(strcasecmp(dict[23].keyarg,"on")==0)  {
           general_data->simopts.anneal_opt = 1;           
           ifound++;
        }
        if(strcasecmp(dict[23].keyarg,"off")==0)  {
           general_data->simopts.anneal_opt = 0;           
           ifound++;
        }
        index=23;
      if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);

      if((general_data->simopts.anneal_opt == 1) && (general_data->simopts.ann_rate == 1.0)){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("Dude, you selected the annealing option but have an \n");
         printf("annealing rate of 1.0.  Maybe you want to reconsider \n");
         printf("something more sensible for the latter\n");
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/

  /*-----------------------------------------------------------------------*/
  /* 24)\ann_start_temperature{#} */

       sscanf(dict[24].keyarg,"%lg",&real_key_arg);
       general_data->simopts.ann_start_temp = real_key_arg;
       index=24;
       if(general_data->simopts.ann_start_temp < 0.0){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}


  /*-----------------------------------------------------------------------*/
  /* 25)\ann_target_temperature{#} */

       sscanf(dict[25].keyarg,"%lg",&real_key_arg);
       general_data->simopts.ann_target_temp = real_key_arg;
       index=25;
       if(general_data->simopts.ann_start_temp < 0.0){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_vpot(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded,
                         CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  double real_key_arg;

/*========================================================================*/
  /* 1)\shift_inter_pe{on,off} */
      ifound      = 0;
      class->energy_ctrl.iswit_vdw=0;
      if(strcasecmp(dict[1].keyarg,"swit")==0) {
         class_parse->ishift_pot = 2;class->energy_ctrl.iswit_vdw=1;ifound++;}
      if(strcasecmp(dict[1].keyarg,"on")==0) {
         class_parse->ishift_pot = 1;ifound++;}
      if(strcasecmp(dict[1].keyarg,"off")==0){
         class_parse->ishift_pot = 0;ifound++;}
      class->interact.iswit_vdw=class->energy_ctrl.iswit_vdw;
      general_data->stat_avg.iswit_vdw=class->energy_ctrl.iswit_vdw;
      index=1;
      if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class_parse->ishift_pot==0){
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
         printf("Potential shift is off.  Energy conservation will degrade.\n");
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
      if(class_parse->ishift_pot==1){
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
         printf("Potential shift is on.  Potential energy will change.\n");
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\inter_spline_pts{#} */
      sscanf(dict[2].keyarg,"%lg",&real_key_arg);
      class->interact.nsplin  = (int)(real_key_arg);
      bonded->ecor.nsplin     = (int)(real_key_arg);
      bonded->ecor.nsplin_m2  = bonded->ecor.nsplin-2;
      index=2;
      if(class->interact.nsplin <= 0) 
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class->interact.nsplin > 2000 ){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested more than 1000 spline points   \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\intra_block_min{} */
      sscanf(dict[3].keyarg,"%lg",&real_key_arg);
      class->energy_ctrl.nblock_min = (int)(real_key_arg);
      index=3;
      if(class->energy_ctrl.nblock_min<0)
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class->energy_ctrl.nblock_min<25){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested a memory blocking length less than 25.\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\pten_inter_respa{#} */
      index=4;
      sscanf(dict[4].keyarg,"%lg",&real_key_arg);
      class->interact.pten_inter_guess = real_key_arg*PCONV;
  /*-----------------------------------------------------------------------*/ 
  /* 5)\pten_kin_respa{#} */
      index=5;
      sscanf(dict[5].keyarg,"%lg",&real_key_arg);
      class->interact.pten_kin_guess = real_key_arg*PCONV;
  /*-----------------------------------------------------------------------*/ 
  /* 6)\pseud_spline_pts{#} */
       sscanf(dict[6].keyarg,"%lg",&real_key_arg);
       cp->pseudo.nsplin_g = (int)(real_key_arg);
       index=6;
       if(cp->pseudo.nsplin_g <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\ewald_interp_pme{#} */
       sscanf(dict[13].keyarg,"%lg",&real_key_arg);
       class->part_mesh.n_interp= (int)(real_key_arg);
       cp->cpscr.cpscr_atom_pme.n_interp = class->part_mesh.n_interp;
       index=13;
       if((class->part_mesh.n_interp  <4 )||
       ((class->part_mesh.n_interp  %2)!=0))
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\scratch_length{#} */
       sscanf(dict[7].keyarg,"%lg",&real_key_arg);
       bonded->intra_scr.nlen = (int)(real_key_arg);
       class->for_scr.nlen = bonded->intra_scr.nlen;
       class->part_mesh.nlen_pme = 
                      class->for_scr.nlen/class->part_mesh.n_interp;
       if(class->part_mesh.nlen_pme < 2*class->part_mesh.n_interp){
       class->part_mesh.nlen_pme = 2*class->part_mesh.n_interp;
       }/*endif*/
       cp->cpscr.cpscr_atom_pme.nlen_pme = class->part_mesh.nlen_pme;
       index=7;
       if(bonded->intra_scr.nlen<=0)
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(bonded->intra_scr.nlen < 50 ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a scratch length < 50   \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       if(bonded->intra_scr.nlen > 1500 ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a scratch length < 1500   \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 8)\ewald_alpha{#} */
       sscanf(dict[8].keyarg,"%lg",&(general_data->ewald.alp_ewd));
       index=8;
       if(general_data->ewald.alp_ewd <= 0.)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->ewald.alp_ewd>40){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested an alpha ewald parameter > 40\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\ewald_kmax{#} */
       sscanf(dict[9].keyarg,"%lg",&real_key_arg);
       class_parse->kmax_ewd = (int)(real_key_arg);
       index=9;
       if(class_parse->kmax_ewd <0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class_parse->kmax_ewd>40){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a kmax ewald parameter > 40\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       if((double)(class_parse->kmax_ewd)<general_data->ewald.alp_ewd){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a kmax ewald parameter < alpha ewald\n");
       printf("Are you certain this is what you would like to do?\n");
       printf("For well converged forces kmax ewald >= alpha ewald \n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\ewald_respa_kmax{#} */
       sscanf(dict[10].keyarg,"%lg",&real_key_arg);
       class_parse->kmax_res = (int)(real_key_arg);
       index=10;
       if((class_parse->kmax_res < 0)||
       (class_parse->kmax_res>class_parse->kmax_ewd)){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/ 
  /* 11)\ewald_pme_opt{#} */
       ifound      = 0;
       if(strcasecmp(dict[11].keyarg,"on")==0) {
       class->part_mesh.pme_on = 1; ifound++;}
       if(strcasecmp(dict[11].keyarg,"off")==0){
       class->part_mesh.pme_on = 0; ifound++;}
       cp->cpscr.cpscr_atom_pme.pme_on = class->part_mesh.pme_on;
       index=11;
    if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\ewald_kmax_pme{#} */
       sscanf(dict[12].keyarg,"%lg",&real_key_arg);
       class->part_mesh.kmax_pme = (int)(real_key_arg);
       index=12;
       if(class->part_mesh.kmax_pme  <0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->part_mesh.pme_on == 1){
        if(class->part_mesh.kmax_pme >40){
          printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
          printf("You have requested a PME kmax ewald parameter > 40\n");
          printf("Are you certain this is what you would like to do?\n");
          printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }/*endif*/
        if(class_parse->kmax_ewd>class->part_mesh.kmax_pme){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
          printf("The ewald PME kmax parameter < ewald kmax\n");
          printf("This is not allowed\n");
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }/*endif*/
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 14)\ewald_respa_pme_opt{#} */
       ifound      = 0;
       if(strcasecmp(dict[14].keyarg,"on")==0) {
       class->part_mesh.pme_res_on = 1; ifound++;}
       if(strcasecmp(dict[14].keyarg,"off")==0){
       class->part_mesh.pme_res_on = 0; ifound++;}
       index=14;
    if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 15)\ewald_respa_kmax_pme{#} */
       sscanf(dict[15].keyarg,"%lg",&real_key_arg);
       class->part_mesh.kmax_pme_res = (int)(real_key_arg);
       index=15;
       if(class->part_mesh.kmax_pme_res  <0) 
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->part_mesh.pme_res_on == 1){
        if(class->part_mesh.kmax_pme_res >40){
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
         printf("You have requested a RESPA PME kmax ewald parameter > 40\n");
         printf("Are you certain this is what you would like to do?\n");
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }/*endif*/
        if(class_parse->kmax_res>class->part_mesh.kmax_pme_res){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("The RESPA PME kmax parameter assigned < ewald respa kmax\n");
         printf("This is not allowed\n");      
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        }/*endif*/
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 16)\ewald_respa_interp_pme{#} */
       sscanf(dict[16].keyarg,"%lg",&real_key_arg);
       class->part_mesh.n_interp_res= (int)(real_key_arg);
       index=16;
       if((class->part_mesh.n_interp_res  <4 )||
       ((class->part_mesh.n_interp_res  %2)!=0))
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  /*-----------------------------------------------------------------------*/ 
  /* 17)\sep_VanderWaals{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[17].keyarg,"on")==0) {
       class->energy_ctrl.isep_vvdw = 1;ifound++;}
       if(strcasecmp(dict[17].keyarg,"off")==0){
       class->energy_ctrl.isep_vvdw = 0;ifound++;}
       index=17;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 18)\dielectric_opt{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[18].keyarg,"on")==0) {
       class->interact.dielectric_opt = 1;ifound++;}
       if(strcasecmp(dict[18].keyarg,"off")==0){
       class->interact.dielectric_opt = 0;ifound++;}
       index=18;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 19)\dielectric_rheal{#} */
       sscanf(dict[19].keyarg,"%lg",&real_key_arg);
       class->interact.dielectric_rheal = (double)(real_key_arg)/0.529177;
       index=19;
       if(class->interact.dielectric_rheal<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 20)\dielectric_cut{#} */
       sscanf(dict[20].keyarg,"%lg",&real_key_arg);
       class->interact.dielectric_cut = (double)(real_key_arg)/0.529177;
       index=20;
       if(class->interact.dielectric_cut<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 21)\dielectric_eps{#} */
       sscanf(dict[21].keyarg,"%lg",&real_key_arg);
       class->interact.dielectric_eps = (double)(real_key_arg);
       index=21;
       if(class->interact.dielectric_eps<1)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 22)\std_intra_block{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[22].keyarg,"on")==0) {
       class->energy_ctrl.block_std_on = 1; ifound++;}
       if(strcasecmp(dict[22].keyarg,"off")==0){
       class->energy_ctrl.block_std_on = 0; ifound++;}
       index=22;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
#ifndef NO_PRAGMA
       if(class->energy_ctrl.block_std_on == 1){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested the blocking option on a scalar machine \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
#endif
  /*-----------------------------------------------------------------------*/
  /* 23)\con_intra_block{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[23].keyarg,"on")==0) {
       class->energy_ctrl.block_con_on = 1; ifound++;}
       if(strcasecmp(dict[23].keyarg,"off")==0){
       class->energy_ctrl.block_con_on = 0; ifound++;}
       index=23;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 24)\inter_PE_calc_freq */
       sscanf(dict[24].keyarg,"%lg",&real_key_arg); 
       general_data->timeinfo.iget_pe_real_inter_freq = (int)(real_key_arg);
       class->energy_ctrl.iget_pe_real_inter_freq     = (int)(real_key_arg);
       class->energy_ctrl.itime              = 0;
       class->energy_ctrl.iget_pe_real_inter = 1;
       class->energy_ctrl.iget_pv_real_inter = 1;
       index=24;
       if(general_data->timeinfo.iget_pe_real_inter_freq < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->timeinfo.iget_pe_real_inter_freq > 10){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested an inter PE calc freq > 10 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  25)\pme_paralell_opt{#} */
       ifound = 0;
       if(strcasecmp(dict[25].keyarg,"hybrid")==0)    {
          class->part_mesh.pme_para_opt = 1; ifound++;}
       if(strcasecmp(dict[25].keyarg,"full_g")==0)    {
          class->part_mesh.pme_para_opt = 2; ifound++;}
       if(strcasecmp(dict[25].keyarg,"none")==0)    {
          class->part_mesh.pme_para_opt = 0; ifound++;}
       index=25;
       if(ifound != 1){
          keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
  /*-----------------------------------------------------------------------*/
/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_run(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                        CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse,
                        DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index,pimd_on;
  int iii;
  double real_key_arg;

  pimd_on = general_data->simopts.pimd + general_data->simopts.cp_pimd 
              + general_data->simopts.cp_wave_pimd 
              + general_data->simopts.debug_pimd 
              + general_data->simopts.debug_cp_pimd
              + general_data->simopts.cp_wave_min_pimd;


/*========================================================================*/
  /* 1)\init_resmpl_atm_vel{on,off} */
        ifound      = 0;
        if(strcasecmp(dict[1].keyarg,"on")==0) {
        class_parse->ivx_smpl = 1;ifound++;}
        if(strcasecmp(dict[1].keyarg,"off")==0){
        class_parse->ivx_smpl = 0;ifound++;}
        index=1;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 2)\init_resmpl_cp_vel{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[2].keyarg,"on")==0)   {
       cp_parse->ivc_smpl = 1;ifound++;}
       if(strcasecmp(dict[2].keyarg,"off")==0)  {
       cp_parse->ivc_smpl = 0;ifound++;}
       index=2;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\init_resmpl_cp_nhc{on,off} */
        ifound      = 0;
        if(strcasecmp(dict[3].keyarg,"on")==0)  {
        cp_parse->ivcnhc_smpl = 1;ifound++;}
        if(strcasecmp(dict[3].keyarg,"off")==0) {
        cp_parse->ivcnhc_smpl = 0;ifound++;}
        index=3;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4a)\resmpl_frq_atm_vel{#} */
        sscanf(dict[4].keyarg,"%lg",&real_key_arg);
        class->vel_samp_class.nvx_smpl = (int)(real_key_arg);
        index=4;
        if(class->vel_samp_class.nvx_smpl<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4b)\rescale_frq_atm_vel{#} */
        sscanf(dict[25].keyarg,"%lg",&real_key_arg);
        class->vel_samp_class.nvx_scale = (int)(real_key_arg);
        index=25;
        if(class->vel_samp_class.nvx_scale<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);


  /*-----------------------------------------------------------------------*/ 
  /* 5)\respa_steps_lrf{#} */
        sscanf(dict[5].keyarg,"%lg",&real_key_arg);
        general_data->timeinfo.nres_ter = (int)(real_key_arg);
        index=5;
        if(general_data->timeinfo.nres_ter<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        general_data->timeinfo.int_res_ter=1;
        class->energy_ctrl.int_res_ter=1;
 
        if(general_data->timeinfo.nres_ter==1){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("You have requested one long-range forces RESPA step.  \n");
        printf(
      "This will result in additional overhead with no gain in efficiency \n");
        printf("Are you certain this is what you would like to do?\n");
        printf("It would be better not to specify the key word or      \n");
        printf("give it the value zero.                                \n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }
        if(general_data->timeinfo.nres_ter==0){
        general_data->timeinfo.int_res_ter=0;
        class->energy_ctrl.int_res_ter=0;
        general_data->timeinfo.nres_ter=1;
        }
  /*-----------------------------------------------------------------------*/ 
  /* 6)\respa_steps_torsion{#} */
        sscanf(dict[6].keyarg,"%lg",&real_key_arg);
        general_data->timeinfo.nres_tor = (int)(real_key_arg);
        index=6;
        if(general_data->timeinfo.nres_tor<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        general_data->timeinfo.int_res_tor=1;
        if(general_data->timeinfo.nres_tor==0){
        general_data->timeinfo.int_res_tor=0;
        general_data->timeinfo.nres_tor=1;
        }
  /*-----------------------------------------------------------------------*/ 
  /* 7)\respa_steps_intra{#} */
        sscanf(dict[7].keyarg,"%lg",&real_key_arg);
        general_data->timeinfo.nres_tra = (int)(real_key_arg);
        index=7;
        if(general_data->timeinfo.nres_tra<0)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        general_data->timeinfo.int_res_tra=1;
        class->energy_ctrl.int_res_tra=1;
        if(general_data->timeinfo.nres_tra==0){
        general_data->timeinfo.int_res_tra=0;
        class->energy_ctrl.int_res_tra=0;
        general_data->timeinfo.nres_tra=1;
        }
  /*-----------------------------------------------------------------------*/ 
  /* 8)\respa_rheal{#} */
       sscanf(dict[8].keyarg,"%lg",&(class->interact.rheal_res));
       index=8;
       if(class->interact.rheal_res<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->interact.rheal_res>2.5){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested RESPA healing length > 2.5A   \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       class->interact.rheal_res  /= BOHR;
  /*-----------------------------------------------------------------------*/ 
  /* 9)\shake_tol{#} */
      sscanf(dict[9].keyarg,"%lg",&(bonded->constrnt.tolshake));
      index=9;
      if(bonded->constrnt.tolshake<=0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(bonded->constrnt.tolshake<1.0e-08) {
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      printf("You have requested a SHAKE tolerence < 1.0e-08    \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }
  /*-----------------------------------------------------------------------*/ 
  /* 10)\rattle_tol{#} */
      sscanf(dict[10].keyarg,"%lg",&(bonded->constrnt.tolratl));
      index=10;
      if(bonded->constrnt.tolratl <= 0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(bonded->constrnt.tolratl<1.0e-08) {
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      printf("You have requested a RATTLE tolerence < 1.0e-08    \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }
  /*-----------------------------------------------------------------------*/ 
  /* 11)\max_constrnt_iter{#} */
      sscanf(dict[11].keyarg,"%lg",&real_key_arg);
      bonded->constrnt.max_iter = (int)(real_key_arg);
      index=11;
      if(bonded->constrnt.max_iter <= 0) 
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      bonded->grp_bond_con.max_iter = bonded->constrnt.max_iter;
      if(bonded->constrnt.max_iter<100){ 
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested a max_constrnt_iter < 100  \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\init_rescale_atm_vel{on,off} */
      ifound      = 0;
      if(strcasecmp(dict[12].keyarg,"on")==0) {
      class_parse->ivx_scale = 1;ifound++;}
      if(strcasecmp(dict[12].keyarg,"off")==0){
      class_parse->ivx_scale = 0;ifound++;}
      index=12;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\init_rescale_atm_nhc{on,off} */
      ifound      = 0;
      if(strcasecmp(dict[13].keyarg,"on")==0) {
      class_parse->ivnhc_scale = 1;ifound++;}
      if(strcasecmp(dict[13].keyarg,"off")==0){
      class_parse->ivnhc_scale = 0;ifound++;}
      index=13;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 14)\init_rescale_cp_vel{on,off} */
      ifound      = 0;
      if(strcasecmp(dict[14].keyarg,"on")==0)   {
      cp_parse->ivc_scale = 1;ifound++;}
      if(strcasecmp(dict[14].keyarg,"off")==0)  {
      cp_parse->ivc_scale = 0;ifound++;}
      index=14;
      if((cp->cpopts.cp_isok_opt == 1) && (cp_parse->ivc_scale != 1)){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("When using the CP isokinetic method, it is advisable\n");
      printf("to do an initial rescale of your CP velocities\n");
      printf("Therefore, I am setting this flag for you.\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      }/* endif */
      if(cp->cpopts.cp_isok_opt == 1) cp_parse->ivc_scale = 1;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 15)\resmpl_frq_cp_vel{ # } */
      sscanf(dict[15].keyarg,"%lg",&real_key_arg);
      cp->vel_samp_cp.nvc_smpl = (int)(real_key_arg);
      index=15;
      if(cp->vel_samp_cp.nvc_smpl<0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 16)\group_con_tol{#} */
      sscanf(dict[16].keyarg,"%lg",&(bonded->grp_bond_con.tol));
      index=16;
      if(bonded->grp_bond_con.tol <= 0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(bonded->grp_bond_con.tol<1.0e-08) {
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      printf("You have requested a GROUP CONSTRAINT tolerence < 1.0e-08\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_norb_tol{#} */
      sscanf(dict[17].keyarg,"%lg",&(cp->cpconstrnt.c_tolnorb));
      index=17;
      if(cp->cpconstrnt.c_tolnorb<=0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);

      if(cp->cpconstrnt.c_tolnorb>0.01&&cp->cpopts.cp_norb==2){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
         printf("The tolerence must be quite low with the    \n");
         printf("norm_only norb option. The rotation does not\n");
         printf("preserve the individual norms only the sum.  \n");
         printf("Work is in progress to overcome this.        \n");
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_ks_rot_frq{#} */
      cp->cpcoeffs_info.ks_rot_on = 1;      
      sscanf(dict[18].keyarg,"%lg",&real_key_arg);
      cp->cpcoeffs_info.n_ks_rot = (int)(real_key_arg);
      index=18;
      if(cp->cpcoeffs_info.n_ks_rot<0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(cp->cpcoeffs_info.n_ks_rot==0){
      cp->cpcoeffs_info.ks_rot_on=0;
      cp->cpcoeffs_info.n_ks_rot=1;
      }
      if(cp->cpcoeffs_info.ks_rot_on == 1 && cp->cpopts.cp_norb > 0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");    
      printf("You are studying the Kohn-Sham eigenvalues under cp_norb \n");
      printf("Eigenvalues are not guaranteed to be correct\n");
      printf("Turn norb off to insure the correct eigenvalues\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
      }
   /*-----------------------------------------------------------------------*/ 
   /* 19)\cp_shake_tol{#} */
       sscanf(dict[19].keyarg,"%lg",&(cp->cpconstrnt.c_tolshake));
       index=19;
       if(cp->cpconstrnt.c_tolshake<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 20)\cp_rattle_tol{#} */
       sscanf(dict[20].keyarg,"%lg",&(cp->cpconstrnt.c_tolratl));
       index=20;
       if(cp->cpconstrnt.c_tolratl <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 24)\cp_min_tol{#} */
       sscanf(dict[24].keyarg,"%lg",&(general_data->minopts.tol_coef));
       index=24;
       if(general_data->minopts.tol_coef<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/

  /*-----------------------------------------------------------------------*/ 
  /* 21)\cp_run_tol{#} */
       sscanf(dict[21].keyarg,"%lg",&(cp->cpopts.tol_coef));
       index=21;
       if(cp->cpopts.tol_coef<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(cp->cpopts.tol_coef < general_data->minopts.tol_coef ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a cp run tolerence less than your \n");
       printf("minimization tolerence %g %g \n",cp->cpopts.tol_coef,
                                             general_data->minopts.tol_coef);
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       if(cp->cpopts.tol_coef > 10000.0*general_data->minopts.tol_coef ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a cp run tolerence greater than  \n");
       printf("10000 times the minimization tolerence %g %g \n",
                                             cp->cpopts.tol_coef,
                                             general_data->minopts.tol_coef);
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 22)\zero_com_vel{yes,no}*/
       ifound = 0;
       if(strcasecmp(dict[22].keyarg,"yes")==0){class_parse->zero_com_vel=1;ifound++;}
       if(strcasecmp(dict[22].keyarg,"no" )==0){class_parse->zero_com_vel=0;ifound++;}
       index = 22;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);  
  /*-----------------------------------------------------------------------*/
  /* 23)\min_tol{#} */
       sscanf(dict[23].keyarg,"%lg",&(general_data->minopts.tol_atom));
       index=23;
       if(general_data->minopts.tol_atom<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 26)hess_opt{full_an,full_num,unit} */
      ifound      = 0;
      if(strcasecmp(dict[26].keyarg,"full_num")==0)   {
      class->clatoms_info.hess_calc = 2;ifound++;}
      if(strcasecmp(dict[26].keyarg,"full_an")==0)   {
      class->clatoms_info.hess_calc = 3;ifound++;}
      if(strcasecmp(dict[26].keyarg,"diag")==0)   {
      class->clatoms_info.hess_calc = 1;ifound++;}
      if(strcasecmp(dict[26].keyarg,"unit")==0)  {
      class->clatoms_info.hess_calc = 0;ifound++;}
      index=26;
      if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(general_data->simopts.cp_min != 1 || general_data->simopts.debug_cp != 1) 
/*        class->clatoms_info.hess_calc=0; */
          iii=0;
      if((general_data->minopts.min_cg == 1) &&  class->clatoms_info.hess_calc > 0 ){
         printf("$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$\n");    
         printf("You have chosen to do a conjugate gradient minimization \n");
         printf("with a nonunit hessian preconditioner.                  \n");
         printf("There is no rigorous basis for this approach, and it    \n");
         printf("is likely not to work.  I will let you try as an        \n");
         printf("experiment, however, if you find that the minimizer is  \n");
         printf("is not behaving well, consider choosing the \"unit\" option\n");
         printf("for the \"hess_opt\" keyword\n");
         printf("$$$$$$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 27)\class_mass_scale_fact{#} */
       sscanf(dict[27].keyarg,"%lg",&(class->clatoms_info.mass_sc_fact));
       index=27;
       if(class->clatoms_info.mass_sc_fact<=0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if((general_data->simopts.cp_min + general_data->simopts.minimize == 0) &&
         class->clatoms_info.mass_sc_fact != 1.0 ){
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
         printf("Dude, if you are not doing a minimization run     \n");
         printf("then, like, why are you scaling your atom masses??\n");
         printf("It's cool!  I'll reset the scaling factor for you!\n");
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
         class->clatoms_info.mass_sc_fact = 1.0;
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 28)\hmat_int_typ{normal,upper_triangle} */
      ifound      = 0;
      if(strcasecmp(dict[28].keyarg,"normal")==0)   {
      general_data->cell.hmat_int_typ = 0;ifound++;}
      if(strcasecmp(dict[28].keyarg,"upper_triangle")==0)  {
      general_data->cell.hmat_int_typ = 1;ifound++;}
      index=28;

      if(general_data->ensopts.npt_f==1){
        if((pimd_on>0)&&(general_data->cell.hmat_int_typ==0)){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("You are running a path integral simulation under\n");
          printf("NPTF.  You must use the upper diagonal form of the\n");
          printf("cell matrix.\n");
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          Finalize();
          exit(1);
        }/*endif*/
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 29)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
      ifound      = 0;
      if(strcasecmp(dict[29].keyarg,"none")==0)   {
              general_data->cell.hmat_cons_typ = 0;ifound++;}
      if(strcasecmp(dict[29].keyarg,"ortho_rhom")==0)  {
              general_data->cell.hmat_cons_typ = 1;ifound++;}
      if(strcasecmp(dict[29].keyarg,"mono_clin")==0)  {
              general_data->cell.hmat_cons_typ = 2;ifound++;}

      index=29;
      if(ifound == 0){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 30)\hmat_cons_typ{none,ortho_rhom,mono_clin} */
      ifound      = 0;
      if(strcasecmp(dict[30].keyarg,"yes")==0)   {
              general_data->minopts.min_atm_com_fix_opt = 1;ifound++;}
      if(strcasecmp(dict[30].keyarg,"no")==0)  {
              general_data->minopts.min_atm_com_fix_opt = 0;ifound++;}
      index=30;
      if(ifound == 0){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     }/*endif*/

  /*-----------------------------------------------------------------------*/ 
  /* 31)\rescale_frq_cp_vel{ # } */
      sscanf(dict[31].keyarg,"%lg",&real_key_arg);
      cp->vel_samp_cp.nvc_scal = (int)(real_key_arg);
      index=31;
      if(cp->vel_samp_cp.nvc_scal<0)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 32)\auto_rescale_cp_vel{ on/off } */
      ifound = 0;
      if(strcasecmp(dict[32].keyarg,"on")==0)   {
              cp->vel_samp_cp.iauto_vc_scal_opt = 1;ifound++;}
      if(strcasecmp(dict[32].keyarg,"off")==0)  {
              cp->vel_samp_cp.iauto_vc_scal_opt = 0;ifound++;}
      index=32;
      if(ifound == 0){
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 33)\auto_rescale_cp_vel_tol{ # } */
      sscanf(dict[33].keyarg,"%lg",&real_key_arg);
      cp->vel_samp_cp.vc_scal_tol = real_key_arg;
      index=33;
      if(cp->vel_samp_cp.vc_scal_tol<=1.0){
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
      if(cp->vel_samp_cp.vc_scal_tol<3.0){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("Coef velocity Rescale tolerence set less than 3kT\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      }/*endif*/
      if(cp->vel_samp_cp.vc_scal_tol>10.0){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("Coef velocity Rescale tolerence set greater than 10kT\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      }/*endif*/

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_nhc(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                        CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse,
                        DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  double real_key_arg;

/*========================================================================*/
  /* 1)\respa_steps_nhc{#} */
      sscanf(dict[1].keyarg,"%lg",&real_key_arg);
      class->therm_info_bead.nres_nhc = (int)(real_key_arg);
      class->therm_info_class.nres_nhc = (int)(real_key_arg);
      index=1;
      if(class->therm_info_class.nres_nhc < 1) 
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class->therm_info_class.nres_nhc>5){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested more than 5 NHC RESPA steps    \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\yosh_steps_nhc{#} */
      sscanf(dict[2].keyarg,"%lg",&real_key_arg);
      class->therm_info_bead.nyosh_nhc = (int)(real_key_arg);
      class->therm_info_class.nyosh_nhc = (int)(real_key_arg);
      index=2;
      if(class->therm_info_class.nyosh_nhc != 1 && 
      class->therm_info_class.nyosh_nhc != 3 && 
      class->therm_info_class.nyosh_nhc != 5 && 
      class->therm_info_class.nyosh_nhc != 7 &&
      class->therm_info_class.nyosh_nhc != 9 )
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\atm_nhc_tau_def{#} */
      sscanf(dict[3].keyarg,"%lg",&(class_parse->tau_nhc_def));
      index=3;
      if(class_parse->tau_nhc_def <= 0) 
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class_parse->tau_nhc_def < 10.0 ){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested a NHC particle default tau < 10fs\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\cp_nhc_len{#} */
      sscanf(dict[4].keyarg,"%lg",&real_key_arg);
      cp->cptherm_info.len_c_nhc = (int)(real_key_arg);
      index=4;
      if(cp->cptherm_info.len_c_nhc <= 1) 
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 5)\resmpl_frq_cp_nhc{ #} */
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       cp->vel_samp_cp.nvcnhc_smpl = (int)(real_key_arg);
       index=5;
       if(cp->vel_samp_cp.nvcnhc_smpl<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 6)\init_resmpl_atm_nhc{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[6].keyarg,"on")==0) {
       class_parse->ivnhc_smpl = 1;ifound++;}
       if(strcasecmp(dict[6].keyarg,"off")==0){
       class_parse->ivnhc_smpl = 0;ifound++;}
       index=6;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\init_resmpl_cp_nhc{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[7].keyarg,"on")==0)  {
       cp_parse->ivcnhc_scale = 1;ifound++;}
       if(strcasecmp(dict[7].keyarg,"off")==0) {
       cp_parse->ivcnhc_scale = 0;ifound++;}
       index=7;
  if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 8)\resmpl_atm_nhc{#} */
       sscanf(dict[8].keyarg,"%lg",&real_key_arg);
       class->vel_samp_class.nvnhc_smpl = (int)(real_key_arg);
       index=8;
       if(class->vel_samp_class.nvnhc_smpl<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 9)\atm_nhc_len{#} */                       
       sscanf(dict[9].keyarg,"%lg",&real_key_arg);
       class->therm_info_bead.len_nhc  = (int)(real_key_arg);
       class->therm_info_class.len_nhc = (int)(real_key_arg);
       general_data->baro.len_nhc      = (int)(real_key_arg);
       index=9;
       if(class->therm_info_bead.len_nhc < 0){ 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);}
       if(class->therm_info_bead.therm_typ == 1
              && class->therm_info_bead.len_nhc > 3 ){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        printf("You have requested a NHC length > 3      \n");
        printf("Are you certain this is what you would like to do?\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 10)\cp_nhc_tau_def{#} */
       sscanf(dict[10].keyarg,"%lg",&(cp_parse->cp_tau_nhc_def));
       index=10;
       if(cp_parse->cp_tau_nhc_def<=0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_respa_steps_nhc{#} */
       sscanf(dict[11].keyarg,"%lg",&real_key_arg);
       cp->cptherm_info.nres_c_nhc = (int)(real_key_arg);
       index=11;
       if(cp->cptherm_info.nres_c_nhc < 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_yosh_steps_nhc{#} */
       sscanf(dict[12].keyarg,"%lg",&real_key_arg);
       cp->cptherm_info.nyosh_c_nhc = (int)(real_key_arg);
       index=12;
  if(cp->cptherm_info.nyosh_c_nhc != 1 && cp->cptherm_info.nyosh_c_nhc != 3 && 
       cp->cptherm_info.nyosh_c_nhc != 5 && cp->cptherm_info.nyosh_c_nhc != 7
       && cp->cptherm_info.nyosh_c_nhc != 9)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\respa_xi_opt{#} */
        sscanf(dict[13].keyarg,"%lg",&real_key_arg);
        general_data->timeinfo.ix_respa = (int)(real_key_arg);
        index=13;
        if(general_data->timeinfo.ix_respa <= 0 
                  || general_data->timeinfo.ix_respa > 5) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        if(general_data->timeinfo.ix_respa>1) {
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        printf("WARNING: xi respa option chosen!\n");
        printf("energy drifts may occur if the number \n");
        printf("of respa steps is taken too large\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }
  /*-----------------------------------------------------------------------*/
  /* 14)\thermostat_type{#} */

       if(strcasecmp(dict[14].keyarg,"NHC")==0)  {
         class->therm_info_bead.therm_typ = 1;
         class->therm_info_class.therm_typ = 1;
       }/*endif*/

       if(strcasecmp(dict[14].keyarg,"GGMT")==0)  {
         class->therm_info_bead.therm_typ = 2;
         class->therm_info_class.therm_typ = 2;
       }/*endif*/

       index=14;
       if(class->therm_info_class.therm_typ <=0
           || class->therm_info_class.therm_typ >= 3) {
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       }

       if(class->therm_info_bead.therm_typ == 2
         && class->therm_info_bead.len_nhc > 3 ){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have requested a GGMT length > 3      \n");
        printf("Available GGMT lengths are 2 or 3    \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        exit(1);
       }/*endif*/
   /*-----------------------------------------------------------------------*/ 
  /* 15)\cp_therm_heat_fact{#} */
     sscanf(dict[15].keyarg,"%lg",&real_key_arg);
     cp->cptherm_info.cp_therm_heat_fact = real_key_arg;
     index=15;
     if(real_key_arg<1.0){
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     }/*endif*/

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_vol(CLASS *class,GENERAL_DATA *general_data,BONDED *bonded,
                        CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse,
                        DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  double real_key_arg;
  int cp_dual_grid_opt_on = cp->cpopts.cp_dual_grid_opt;

/*========================================================================*/
  /* 1)\volume_tau{#} */
       sscanf(dict[1].keyarg,"%lg",&(class_parse->tau_vol));
       index=1;
       if(class_parse->tau_vol <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class_parse->tau_vol > 10000.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a volume tau > 10ps           \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       if(class_parse->tau_vol < 100.0 ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a volume tau < 100fs           \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\volume_nhc_tau{#} */
       sscanf(dict[2].keyarg,"%lg",&(class_parse->tau_vol_nhc));
       index=2;
       if(class_parse->tau_vol_nhc <= 0) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class_parse->tau_vol_nhc > 10000.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a NHC volume tau  > 10ps      \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
       if(class_parse->tau_vol_nhc < 100.0 ){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a NHC volume tau < 100fs      \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\periodicity{0,1,2,3} */
       if(strcasecmp(dict[3].keyarg,"0_ewald")==0)  {
         general_data->cell.iperd = 4;
       }else{
         sscanf(dict[3].keyarg,"%lg",&real_key_arg);
         general_data->cell.iperd = (int)(real_key_arg);
       }/*endif*/
       index=3;
       if(general_data->cell.iperd == 4 && cp_dual_grid_opt_on >= 1){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("cluster ewald  not supported for dualed systems yet. \n");
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          Finalize();
          exit(1);
       }/*endif*/
       if((general_data->cell.iperd<0)||(general_data->cell.iperd>4))
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 4)\intra_perds{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[4].keyarg,"on")==0) {
       general_data->cell.intra_perds = 1;ifound++;}
       if(strcasecmp(dict[4].keyarg,"off")==0){
       general_data->cell.intra_perds = 0;ifound++;}
       index=4;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_write(CLASS *class,GENERAL_DATA *general_data,
                          BONDED *bonded,
                          CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                          FILENAME_PARSE *filename_parse,
                          DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */
  int ifound,index;
  double real_key_arg;
  char *strip,*strip2;
/*========================================================================*/

  strip = (char *)cmalloc(MAXWORD*sizeof(char));
  strip2 = (char *)cmalloc(MAXWORD*sizeof(char));

/*========================================================================*/
  /* 1)\write_binary_cp_coef{on,off} */
       ifound  = 0;
       if(strcasecmp(dict[1].keyarg,"on")==0)  
          {cp->cpopts.iwrite_coef_binary = 1;ifound++;}
       if(strcasecmp(dict[1].keyarg,"off")==0) 
          {cp->cpopts.iwrite_coef_binary = 0;ifound++;}
       index=1;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 2)\write_force_freq{#} */
       sscanf(dict[2].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_atm_for = (int)(real_key_arg);
       index=2;
       if(general_data->filenames.iwrite_atm_for < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 3)\write_screen_freq */
       sscanf(dict[3].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_screen = (int)(real_key_arg);
       index=3;
       if(general_data->filenames.iwrite_screen < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_screen > 1000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a write screen freq > 1000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\write_dump_freq */
       sscanf(dict[4].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_dump = (int)(real_key_arg);
       index=4;
       if(general_data->filenames.iwrite_dump < 1) 
         keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_dump > 2000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a write dump freq > 2000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 5)\write_inst_freq */
        sscanf(dict[5].keyarg,"%lg",&real_key_arg);
        general_data->filenames.iwrite_inst = (int)(real_key_arg);
        index=5;
        if(general_data->filenames.iwrite_inst < 1) 
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
        if(general_data->filenames.iwrite_inst > 2000){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        printf("You have requested a write inst freq > 2000 steps \n");
        printf("Are you certain this is what you would like to do?\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\write_pos_freq */
       sscanf(dict[6].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_confp = (int)(real_key_arg);
       index=6;
       if(general_data->filenames.iwrite_confp < 1) 
           keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_confp > 2000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a write position freq > 2000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\write_vel_freq */
       sscanf(dict[7].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_confv = (int)(real_key_arg);
       index=7;
       if(general_data->filenames.iwrite_confv < 1) 
           keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_confv > 2000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a write velocity freq > 2000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 8)\write_cp_c_freq */
       sscanf(dict[8].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_confc = (int)(real_key_arg);
       index=8;
       if(general_data->filenames.iwrite_confc < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 9)\screen_output_units */
       ifound = 0;
       if(strcasecmp(dict[9].keyarg,"au")==0)    {
       general_data->filenames.iwrite_units=0; ifound++;}
       if(strcasecmp(dict[9].keyarg,"kcal_mol")==0)    {
       general_data->filenames.iwrite_units=1; ifound++;}
       if(strcasecmp(dict[9].keyarg,"kelvin")==0)    {
       general_data->filenames.iwrite_units=2; ifound++;}
       index=9;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 10)\conf_file_format */
       ifound = 0;
       if(strcasecmp(dict[10].keyarg,"binary")==0)    {
       general_data->filenames.iwrite_conf_binary=1; ifound++;}
       if(strcasecmp(dict[10].keyarg,"formatted")==0)  {
       general_data->filenames.iwrite_conf_binary=0; ifound++;}
       index=10;
    if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 11)\path_cent_file */
       sscanf(dict[11].keyarg,"%s",general_data->filenames.centname);
  /*-----------------------------------------------------------------------*/ 
  /* 12)\atm_force_file */
       sscanf(dict[12].keyarg,"%s",general_data->filenames.forcename);
  /*-----------------------------------------------------------------------*/ 
  /* 13)\cp_coef_file */
       sscanf(dict[13].keyarg,"%s",general_data->filenames.ccname);
  /*-----------------------------------------------------------------------*/ 
  /* 14)\conf_partial_freq */
       sscanf(dict[14].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_par_confp = (int)(real_key_arg);
       index=14;
       if(general_data->filenames.iwrite_par_confp < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_par_confp > 2000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a partial write pos freq > 2000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 15)\path_cent_freq */
       sscanf(dict[15].keyarg,"%lg",&real_key_arg);
       general_data->filenames.iwrite_path_cent = (int)(real_key_arg);
       index=15;
       if(general_data->filenames.iwrite_path_cent < 1) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->filenames.iwrite_path_cent > 2000){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a write centroid freq > 2000 steps \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 16)\conf_partial_limits */
       strcpy(strip,dict[16].keyarg);
       parse_part_lim(strip,strip2,&(general_data->filenames.low_lim_par),
                   &(general_data->filenames.high_lim_par));     
       index=16;
       if((general_data->filenames.low_lim_par <= 0) || 
        (general_data->filenames.high_lim_par < 0)) 
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if((general_data->filenames.low_lim_par > general_data->filenames.high_lim_par)&&
       (dict[16].iuset==1)){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested a lower limit higher than your \n");
       printf("upper limit in the conf_partial_limits keyarg     \n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/ 
   /* 17)\read_binary_cp_coef{on,off} */
       ifound  = 0;
       if(strcasecmp(dict[17].keyarg,"on")==0)  
          {cp->cpopts.iread_coef_binary = 1;ifound++;}
       if(strcasecmp(dict[17].keyarg,"off")==0) 
          {cp->cpopts.iread_coef_binary = 0;ifound++;}
       index=17;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/  
  /* 18)\sim_name */
       sscanf(dict[18].keyarg,"%s",filename_parse->simname);
  /*-----------------------------------------------------------------------*/ 
  /* 19)\out_restart_file */
       sscanf(dict[19].keyarg,"%s",general_data->filenames.dname);
  /*-----------------------------------------------------------------------*/ 
  /* 20)\in_restart_file */
       sscanf(dict[20].keyarg,"%s",filename_parse->dnamei);
  /*-----------------------------------------------------------------------*/ 
  /* 21)\instant_file */
       sscanf(dict[21].keyarg,"%s",general_data->filenames.iname);
  /*-----------------------------------------------------------------------*/ 
  /* 22)\atm_pos_file */
       sscanf(dict[22].keyarg,"%s",general_data->filenames.cpname);
  /*-----------------------------------------------------------------------*/ 
  /* 23)\atm_vel_file */
       sscanf(dict[23].keyarg,"%s",general_data->filenames.cvname);
  /*-----------------------------------------------------------------------*/ 
  /* 24)\conf_partial_file */
        sscanf(dict[24].keyarg,"%s",general_data->filenames.cpparname);
  /*-----------------------------------------------------------------------*/ 
  /* 25)\mol_set_file */
       sscanf(dict[25].keyarg,"%s",filename_parse->molsetname);
  /*-----------------------------------------------------------------------*/ 
  /* 26)\cp_restart_out_file */
       sscanf(dict[26].keyarg,"%s",general_data->filenames.dnamec);
  /*-----------------------------------------------------------------------*/ 
  /* 27)\cp_restart_in_file */
       sscanf(dict[27].keyarg,"%s",filename_parse->dnameci);
  /*-----------------------------------------------------------------------*/ 
  /* 28)\cp_kseigs_file */
        sscanf(dict[28].keyarg,"%s",general_data->filenames.ksname);
  /*-----------------------------------------------------------------------*/ 
  /* 29)\cp_elf_file */
        sscanf(dict[29].keyarg,"%s",general_data->filenames.elfname);

/*========================================================================*/

  cfree(strip);
  cfree(strip2);

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_pimd(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded,
                         CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  double real_key_arg;

/*========================================================================*/
  /* 1)\path_int_beads{#} */
      sscanf(dict[1].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.pi_beads = (int)(real_key_arg);
      cp->cpcoeffs_info.pi_beads = (int)(real_key_arg);
      general_data->simopts.pi_beads = (int)(real_key_arg);
      index=1;
      if((class->clatoms_info.pi_beads<0))
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      class->clatoms_info.pi_beads_proc = 
              (class->clatoms_info.pi_beads)/(class->communicate.np_beads);
      cp->cpcoeffs_info.pi_beads_proc =   class->clatoms_info.pi_beads_proc;

      if((class->clatoms_info.pi_beads%class->communicate.np_beads)!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The number of processors does not evenly divide \n");
      printf("the number of path integral beads.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      Finalize();
      exit(1);
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\path_int_gamma_adb{#} */
      sscanf(dict[2].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.gamma_adb = (real_key_arg);
      index=2;
      if((class->clatoms_info.gamma_adb<0.0))
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(class->clatoms_info.gamma_adb>1.0){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested a path integral                \n");
      printf("adiabaticity parameter more than one.             \n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\path_int_md_typ{staging,centroid} */
      ifound      = 0;
      if(strcasecmp(dict[3].keyarg,"staging")==0) {
      general_data->simopts.pi_md_typ = 1;ifound++;}
      if(strcasecmp(dict[3].keyarg,"centroid")==0){
      general_data->simopts.pi_md_typ = 2;ifound++;}
      index=3;
   if(ifound != 1) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 4)\pi_beads_level_full{#} */
      sscanf(dict[4].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.pi_beads_full_ter = (int)(real_key_arg);
      index=4;
      if(class->clatoms_info.pi_beads_full_ter<1||
      class->clatoms_info.pi_beads_full_ter>class->clatoms_info.pi_beads)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
      if(dict[4].iuset==0){
      class->clatoms_info.pi_beads_full_ter=class->clatoms_info.pi_beads;
      }/*endif*/
  /*-----------------------------------------------------------------------*/ 
  /* 5)\pi_beads_level_inter_short{#} */
      sscanf(dict[5].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.pi_beads_res_ter = (int)(real_key_arg);
      index=5;
      if(class->clatoms_info.pi_beads_res_ter<0||
      class->clatoms_info.pi_beads_res_ter>class->clatoms_info.pi_beads)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 6)\pi_beads_level_intra_res{#} */
      sscanf(dict[6].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.pi_beads_res_tra = (int)(real_key_arg);
      index=6;
      if(class->clatoms_info.pi_beads_res_tra<0||
      class->clatoms_info.pi_beads_res_tra>class->clatoms_info.pi_beads)
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 7)\pi_beads_level_intra{#} */
      sscanf(dict[7].keyarg,"%lg",&real_key_arg);
      class->clatoms_info.pi_beads_full_tra = (int)(real_key_arg);
      index=7;
      if(class->clatoms_info.pi_beads_full_tra<0||
        class->clatoms_info.pi_beads_full_tra>class->clatoms_info.pi_beads)
        keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/ 
  /* 8)\respa_steps_pimd{#} */
       sscanf(dict[8].keyarg,"%lg",&real_key_arg);
       general_data->timeinfo.nres_pimd = (int)(real_key_arg);
       index=8;
       if(general_data->timeinfo.nres_pimd<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->timeinfo.nres_pimd==0){
       general_data->timeinfo.nres_pimd=1;
       }
  /*-----------------------------------------------------------------------*/ 
  /* 9)\initial_spread_size{} */
       sscanf(dict[9].keyarg,"%lg",&real_key_arg);
       class->clatoms_info.rcut_spread = (real_key_arg)/BOHR;
       index=9;
       if(class->clatoms_info.rcut_spread<0)
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(class->clatoms_info.rcut_spread<1.0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested an initial spread size < 1 Bohr \n");
       printf("Are you certain this is what you would like to do?\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /* 10)\initial_spread_opt{on,off} */
       ifound      = 0;
       if(strcasecmp(dict[10].keyarg,"on")==0) {
           general_data->simopts.initial_spread_opt = 1; ifound++;}
       if(strcasecmp(dict[10].keyarg,"off")==0){
           general_data->simopts.initial_spread_opt = 0; ifound++;}
       index=10;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
       if(general_data->simopts.initial_spread_opt == 1 &&   
                           general_data->simopts.pimd == 0){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You have requested the initial spread option \n");
       printf("\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
       }/*endif*/
  /*-----------------------------------------------------------------------*/
  /*  11)\pimd_freeze_type{centroid,all_mode} */
       ifound      = 0;
       if(strcasecmp(dict[11].keyarg,"centroid")==0) {
           class->atommaps.pimd_freez_typ = 1; ifound++;}
       if(strcasecmp(dict[11].keyarg,"all_mode")==0){
           class->atommaps.pimd_freez_typ = 2; ifound++;}
       index=11;
   if(ifound == 0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_velo(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded, CP *cp,ANALYSIS *analysis,
                         CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  int iii;
  double real_key_arg;

/*========================================================================*/
  /*-----------------------------------------------------------------------*/
  /* 1) \corel_vovt{#} */
  index  = 1;
  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->velocorel.calcul_atm_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->velocorel.calcul_atm_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 2) \corel_vovt_ncor{#} */
  index = 2;
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.ncor = (int)(real_key_arg);
  if (analysis->velocorel.calcul_atm_on==1){
  if((analysis->velocorel.ncor<0)||
     (analysis->velocorel.ncor>general_data->timeinfo.ntime))
     {
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     };
   };
  /*-----------------------------------------------------------------------*/
  /* 3) \corel_vovt_njump{#} */
  index = 3;
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.njump = (int)(real_key_arg);
  if(analysis->velocorel.njump<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 4) \write_vovt */
  sscanf(dict[4].keyarg,"%s",analysis->velocorel.vovtname);
  /*-----------------------------------------------------------------------*/
  /* 5) \corel_vovt_normalize{#} */
  index  = 5;
  ifound = 0;
  if(strcasecmp(dict[5].keyarg,"on")==0) {
     analysis->velocorel.normalize = 1; ifound++;}
  if(strcasecmp(dict[5].keyarg,"off")==0){
     analysis->velocorel.normalize = 0; ifound++;}
  if (ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/
 /* 6) \corel_vovt_output */
  sscanf(dict[6].keyarg,"%s",analysis->velocorel.output_kind);
 /*-----------------------------------------------------------------------*/
  /* 7) \corel_vovt_periodic_output{#} */
  index = 7;
  sscanf(dict[7].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.nwrite = (int)(real_key_arg);
  if(analysis->velocorel.nwrite<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 8) \corel_vovt_com{#} */
  index  = 8;
  ifound = 0;
  if(strcasecmp(dict[8].keyarg,"on")==0) {
    analysis->velocorel.calcul_com_on = 1; ifound++;}
  if(strcasecmp(dict[8].keyarg,"off")==0){
    analysis->velocorel.calcul_com_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 9) \corel_vovt_com_ncor{#} */
  index = 9;
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.ncor_com = (int)(real_key_arg);
  if (analysis->velocorel.calcul_com_on==1){
  if((analysis->velocorel.ncor_com<0)||
     (analysis->velocorel.ncor_com>general_data->timeinfo.ntime))
     {
      keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     };
   };
  /*-----------------------------------------------------------------------*/
  /* 10) \corel_vovt_com_njump{#} */
  index = 10;
  sscanf(dict[10].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.njump_com = (int)(real_key_arg);
  if(analysis->velocorel.njump_com<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 11) \write_vovt_com */
  sscanf(dict[11].keyarg,"%s",analysis->velocorel.vovtname_com);
  /*-----------------------------------------------------------------------*/
  /* 12) \corel_vovt_com_normalize{#} */
  index  = 12;
  ifound = 0;
  if(strcasecmp(dict[12].keyarg,"on")==0) {
     analysis->velocorel.normalize_com = 1; ifound++;}
  if(strcasecmp(dict[12].keyarg,"off")==0){
     analysis->velocorel.normalize_com = 0; ifound++;}
  if (ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/
 /* 13) \corel_vovt_com_output */
  sscanf(dict[13].keyarg,"%s",analysis->velocorel.output_kind_com);
 /*-----------------------------------------------------------------------*/
  /* 14) \corel_vovt_com_periodic_output{#} */
  index = 14;
  sscanf(dict[14].keyarg,"%lg",&real_key_arg);
  analysis->velocorel.nwrite_com = (int)(real_key_arg);
  if(analysis->velocorel.nwrite_com<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/

/*========================================================================*/
    }/* end routine set_sim_params_velo */ 
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_msqd(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded, CP *cp,ANALYSIS *analysis,
                         CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  int iii;
  double real_key_arg;

/*========================================================================*/
/*========================================================================*/
  /*-----------------------------------------------------------------------*/
  /* 1) \corel_msqd{#} */
  index  = 1;
  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->msqdcorel.calcul_atm_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->msqdcorel.calcul_atm_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 2) \corel_msqd_ncor{#} */
  index = 2;
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  analysis->msqdcorel.ncor = (int)(real_key_arg);
  if (analysis->msqdcorel.calcul_atm_on==1){
  if((analysis->msqdcorel.ncor<0)||
     (analysis->msqdcorel.ncor>general_data->timeinfo.ntime))
     {
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     };
   };
  /*-----------------------------------------------------------------------*/
  /* 3) \corel_msqd_njump{#} */
  index = 3;
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  analysis->msqdcorel.njump = (int)(real_key_arg);
  if(analysis->msqdcorel.njump<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/
  /* 4) \write_msqd */
  sscanf(dict[4].keyarg,"%s",analysis->msqdcorel.msqdname);
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_msqd_output */
  sscanf(dict[5].keyarg,"%s",analysis->msqdcorel.output_kind);
 /*-----------------------------------------------------------------------*/

/*========================================================================*/
    }/* end routine set_sim_params_msqd */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_iikt_iso(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded, CP *cp,ANALYSIS *analysis,
                         CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index,count,i;
  int iii,iflag,natm_typ_now;
  char *tmp;
  char *tmp2;
  double real_key_arg;

  tmp   = (char *)cmalloc(1000*sizeof(char));
  tmp2  = (char *)cmalloc(1000*sizeof(char));

/*========================================================================*/
/*========================================================================*/
  /*-----------------------------------------------------------------------*/
  /* 1) \corel_iikt_iso{#} */
  index  = 1;
  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->iikt_iso_corel.calcul_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->iikt_iso_corel.calcul_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 2) \corel_iikt_iso_ncor{#} */
  index = 2;
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  analysis->iikt_iso_corel.ncor = (int)(real_key_arg);
  if (analysis->iikt_iso_corel.calcul_on==1){
  if((analysis->iikt_iso_corel.ncor<0)||
     (analysis->iikt_iso_corel.ncor>general_data->timeinfo.ntime))
     {
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     };
  };
  /*-----------------------------------------------------------------------*/
  /* 3) \corel_iikt_iso_njump{#} */
  index = 3;
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  analysis->iikt_iso_corel.njump = (int)(real_key_arg);
  if(analysis->iikt_iso_corel.njump<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/
  /* 4) \write_iikt_iso */
  sscanf(dict[4].keyarg,"%s",analysis->iikt_iso_corel.iikt_iso_name);
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_iikt_iso_output */
  sscanf(dict[5].keyarg,"%s",analysis->iikt_iso_corel.output_kind);
 /*-----------------------------------------------------------------------*/
 /* 6)\corel_iikt_iso_kmin{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  analysis->iikt_iso_corel.kmin_ask = (double)(real_key_arg);
  index=6;
  if(analysis->iikt_iso_corel.kmin_ask<0.0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/
 /* 7)\corel_iikt_iso_kmax{#} */
  sscanf(dict[7].keyarg,"%lg",&real_key_arg);
  analysis->iikt_iso_corel.kmax_ask = (double)(real_key_arg);
  index=7;
  if(analysis->iikt_iso_corel.kmax_ask<0.0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/
 /* 8)\corel_iikt_iso_nb_kvec{#} */
 index=8;
 sscanf(dict[8].keyarg,"%lg",&real_key_arg);
 analysis->iikt_iso_corel.nbk_ask = (int)(real_key_arg);
 if(analysis->iikt_iso_corel.nbk_ask<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 9) \corel_iikt_iso{#} */
  index  = 9;
  ifound = 0;
  if(strcasecmp(dict[9].keyarg,"on")==0) {
    analysis->iikt_iso_corel.eisf_on = 1; ifound++;}
  if(strcasecmp(dict[9].keyarg,"off")==0){
    analysis->iikt_iso_corel.eisf_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/
 /* 10) \corel_iikt_atoms_list{#} */
  index  = 10;
  ifound = 0;
  if(strcasecmp(dict[10].keyarg,"all")==0)
  {
   analysis->iikt_iso_corel.calc_all_atm_typ = 1; ifound++;
  }else{
   analysis->iikt_iso_corel.calc_all_atm_typ = 0; ifound++;
   strcpy(tmp,dict[10].keyarg);
   natm_typ_now = 1;

   analysis->iikt_iso_corel.list_atm_typ = 
      (NAME *)cmalloc(natm_typ_now*sizeof(NAME))-1;

   iflag = 0;
   count = 0;
   while(iflag != 1)
   {
    parse_atm_typ_site(tmp,tmp2,&natm_typ_now,&iflag,&count);
    sscanf(tmp2,"%s",analysis->iikt_iso_corel.list_atm_typ[natm_typ_now]);
    if(iflag==0)
    {
     natm_typ_now++;
     analysis->iikt_iso_corel.list_atm_typ = (NAME *) 
         crealloc(&(analysis->iikt_iso_corel.list_atm_typ[1]),
                  (natm_typ_now)*sizeof(NAME))-1;
    }; /*endif*/
    if(iflag==2) {break;};
   }; /*endwhile*/
  analysis->iikt_iso_corel.natm_typ_calc = natm_typ_now;
 }; /*end else*/

  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 /*-----------------------------------------------------------------------*/

  cfree(&tmp[0]);
  cfree(&tmp2[0]);
/*========================================================================*/
    }/* end routine set_sim_params_iikt_iso */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_ickt_iso(CLASS *class,GENERAL_DATA *general_data,
                         BONDED *bonded, CP *cp,ANALYSIS *analysis,
                         CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                         FILENAME_PARSE *filename_parse,
                         DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  int iii;
  double real_key_arg;

/*========================================================================*/
/*========================================================================*/
  /*-----------------------------------------------------------------------*/
  /* 1) \corel_ickt_iso{#} */
  index  = 1;
  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->ickt_iso_corel.calcul_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->ickt_iso_corel.calcul_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 2) \corel_ickt_iso_ncor{#} */
  index = 2;
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  analysis->ickt_iso_corel.ncor = (int)(real_key_arg);
  if(analysis->ickt_iso_corel.calcul_on==1){
  if((analysis->ickt_iso_corel.ncor<0)||
     (analysis->ickt_iso_corel.ncor>general_data->timeinfo.ntime))
     {
       keyarg_barf(dict,filename_parse->input_name,fun_key,index);
     };
   };
  /*-----------------------------------------------------------------------*/
  /* 3) \corel_ickt_iso_njump{#} */
  index = 3;
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  analysis->ickt_iso_corel.njump = (int)(real_key_arg);
  if(analysis->ickt_iso_corel.calcul_on==1){
  if(analysis->ickt_iso_corel.njump<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
   };
  /*-----------------------------------------------------------------------*/
  /*-----------------------------------------------------------------------*/
  /* 4) \write_ickt_iso */
  sscanf(dict[4].keyarg,"%s",analysis->ickt_iso_corel.ickt_iso_name);
 /*-----------------------------------------------------------------------*/
 /* 5)\corel_ickt_iso_kmin{#} */
  sscanf(dict[5].keyarg,"%lg",&real_key_arg);
  analysis->ickt_iso_corel.kmin_ask = (double)(real_key_arg);
  index=5;
  if(analysis->ickt_iso_corel.calcul_on==1){
  if(analysis->ickt_iso_corel.kmin_ask<0.0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  };
 /*-----------------------------------------------------------------------*/
 /* 6)\corel_ickt_iso_kmax{#} */
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  analysis->ickt_iso_corel.kmax_ask = (double)(real_key_arg);
  index=6;
  if(analysis->ickt_iso_corel.calcul_on==1){
  if(analysis->ickt_iso_corel.kmax_ask<0.0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  };
 /*-----------------------------------------------------------------------*/
 /* 7)\corel_ickt_iso_nb_kvec{#} */
 index=7;
 sscanf(dict[7].keyarg,"%lg",&real_key_arg);
 analysis->ickt_iso_corel.nbk_ask = (int)(real_key_arg);
 if(analysis->ickt_iso_corel.calcul_on==1){
 if(analysis->ickt_iso_corel.nbk_ask<0)
    keyarg_barf(dict,filename_parse->input_name,fun_key,index);
 };
  /*-----------------------------------------------------------------------*/

/*========================================================================*/
    }/* end routine set_sim_params_ickt_iso */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_rdf(CLASS *class,GENERAL_DATA *general_data,
                        BONDED *bonded, CP *cp,ANALYSIS *analysis,
                        CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FILENAME_PARSE *filename_parse,
                        DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  int iii;
  double real_key_arg;

/*========================================================================*/
/*========================================================================*/
  /*-----------------------------------------------------------------------*/
  /* 1) \corel_rdf{#} */
  index  = 1;
  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->rdf.calcul_gr_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->rdf.calcul_gr_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 2) \corel_rdf_intra{#} */
  index = 2;
  if(strcasecmp(dict[2].keyarg,"on")==0) {
    analysis->rdf.calcul_intra_on = 1; ifound++;}
  if(strcasecmp(dict[2].keyarg,"off")==0){
    analysis->rdf.calcul_intra_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 3) \corel_rdf_njump{#} */
  index = 3;
  sscanf(dict[3].keyarg,"%lg",&real_key_arg);
  analysis->rdf.njump = (int)(real_key_arg);
  if(analysis->rdf.njump<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 4) \corel_rdf_npts{#} */
  index = 4;
  sscanf(dict[4].keyarg,"%lg",&real_key_arg);
  analysis->rdf.number_of_points = (int)(real_key_arg);
  if(analysis->rdf.number_of_points<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 5) \corel_rdf_output_file{*} */
  sscanf(dict[5].keyarg,"%s",analysis->rdf.rdfname);
 /*-----------------------------------------------------------------------*/
  /* 6) \corel_rdf_periodic_output{#} */
  index = 6;
  sscanf(dict[6].keyarg,"%lg",&real_key_arg);
  analysis->rdf.periodic_write = (int)(real_key_arg);
  if(analysis->rdf.periodic_write<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 7) \corel_rdf_sq_from_gr{#} */
  index = 7;
  if(strcasecmp(dict[7].keyarg,"on")==0) {
    analysis->rdf.calcul_sq_from_gr_on = 1; ifound++;}
  if(strcasecmp(dict[7].keyarg,"off")==0){
    analysis->rdf.calcul_sq_from_gr_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 8) \corel_rdf_sq_nbk{#} */
  index = 8;
  sscanf(dict[8].keyarg,"%lg",&real_key_arg);
  analysis->rdf.nbk = (int)(real_key_arg);
  if(analysis->rdf.nbk<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/
  /* 9) \corel_rdf_sq_dk{#} */
  index = 9;
  sscanf(dict[9].keyarg,"%lg",&real_key_arg);
  analysis->rdf.dk = (double)(real_key_arg);
  if(analysis->rdf.dk<0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);
  /*-----------------------------------------------------------------------*/

/*========================================================================*/
} /* end routine set_sim_params_rdf */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_params_harmonic(CLASS *class,GENERAL_DATA *general_data,
                             BONDED *bonded, CP *cp,ANALYSIS *analysis,
                             CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                             FILENAME_PARSE *filename_parse,
                             DICT_WORD *dict,char *fun_key)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index;
  int iii;
  double real_key_arg;

  /*-----------------------------------------------------------------------*/
  /* 1) \harmonic_frequencies{#} */

  index = 1;  ifound = 0;
  if(strcasecmp(dict[1].keyarg,"on")==0) {
    analysis->harmonic_analysis.calcul_freq_on = 1; ifound++;}
  if(strcasecmp(dict[1].keyarg,"off")==0){
    analysis->harmonic_analysis.calcul_freq_on = 0; ifound++;}
  if(ifound==0) keyarg_barf(dict,filename_parse->input_name,fun_key,index);

  if(analysis->harmonic_analysis.calcul_freq_on == 1 ) {
     general_data->simopts.minimize = 0;
     general_data->simopts.md = 0;
     general_data->simopts.pimd = 0;
     general_data->simopts.cp = 0;
     general_data->simopts.cp_wave = 0;
     general_data->simopts.cp_wave_pimd = 0;
     general_data->simopts.cp_pimd = 0; 
     general_data->simopts.cp_min = 1; 
     general_data->simopts.cp_wave_min = 0.;
     general_data->simopts.cp_wave_min_pimd = 0.;
     general_data->simopts.debug = 0;
     general_data->simopts.debug_pimd = 0;
     general_data->simopts.debug_cp = 0;
     general_data->simopts.debug_cp_pimd = 0; 
     general_data->timeinfo.ntime = 100000;
     class->clatoms_info.hess_calc = 2;
   } /* endif */

  /*-----------------------------------------------------------------------*/
  /* 2) \finite_difference_displacement{#} */

  index = 2;
  sscanf(dict[2].keyarg,"%lg",&real_key_arg);
  analysis->harmonic_analysis.delta = (double)(real_key_arg);
  analysis->harmonic_analysis.delta /= BOHR;
  if(analysis->harmonic_analysis.delta < 0)
     keyarg_barf(dict,filename_parse->input_name,fun_key,index);

/*========================================================================*/
} /* end routine set_sim_params_harmonic */
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_sim_params_finale: Consistency checks                               */
/*==========================================================================*/

void set_sim_params_finale(CLASS *class,GENERAL_DATA *general_data,
                           BONDED *bonded,
                           CP *cp,CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                           FILENAME_PARSE *filename_parse)

/*=======================================================================*/
/*             Begin routine                                              */
{/*begin routine */
/*=======================================================================*/
/*             Local variable declarations                                */

  int ifound,index,itemp,pi_beads,cp_on,p1,p2,p3,p4,i1,i2,i3,i4,psum,ip,iii;
  int pimd_on,cppimd_on,cp_orth_sum,ilimit;
  double real_key_arg;
  int cp_min_on;
  int nproc_tot = class->communicate.np;
  int iperd= general_data->cell.iperd;
/*========================================================================*/
/* I) General Consistency Checks                                          */

  if(general_data->ensopts.nve==1){
    class->therm_info_bead.len_nhc = 0;
    class->therm_info_class.len_nhc = 0;
    class->therm_info_bead.num_nhc = 0;
    class->therm_info_class.num_nhc = 0;
  }

 pimd_on = general_data->simopts.pimd + general_data->simopts.cp_pimd 
            + general_data->simopts.cp_wave_pimd 
            + general_data->simopts.debug_pimd 
            + general_data->simopts.debug_cp_pimd
            + general_data->simopts.cp_wave_min_pimd;

 if((class->clatoms_info.pi_beads > 1)&&(pimd_on == 0)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The number of path integral beads is > 1, but \n");
       printf("you have not turned on a path integral simulation\n");
       printf("type.  Therefore, this run will terminate, now.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }/*endif*/

 if((general_data->ensopts.nve==1)&&(pimd_on == 1)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("You are doing a path integral calculation .  The \n");
       printf("microcanonical ensemble (NVE) is meaningless\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }/*endif*/

 if(class->nbr_list.verlist.nmem_min_lst < bonded->intra_scr.nlen){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Your scr_length is greater than your verlist_mem_min\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }

 if(class->energy_ctrl.nblock_min > bonded->intra_scr.nlen){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Your scr_length is less than your minimum block size \n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }

 if(pimd_on==0&&general_data->timeinfo.ix_respa==5){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("The xi_respa option must be 1,2,3, or 4 under \n");
       printf("classical MD \n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }

 if(class_parse->istart>2&&general_data->simopts.initial_spread_opt==1){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Spread option is supported only for restart_pos and\n");
       printf("initial restart types\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }

 if(class->interact.dielectric_opt==1 && general_data->cell.iperd != 0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Position dependent dielectric constant only    \n");
       printf("supported under cluster boundary conditions    \n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
 }

/*========================================================================*/
/* II) Respa consistency checks                                           */

  if(general_data->timeinfo.int_res_ter==1 && 
                     general_data->timeinfo.int_res_tra == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    printf("Having inter molecular respa on and intra molecular respa off\n");
    printf("doesn't makes sense.\n");
    printf("The program therefore insists that if int_res_ter=1\n");
    printf("then int_res_tra=1\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  if(general_data->timeinfo.int_res_ter==1 &&
                     general_data->timeinfo.int_res_tor == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    printf("Having inter molecular respa on and torsion respa off\n");
    printf("doesn'tmakes sense. \n");
    printf("The program therefore insists that if int_res_ter=1\n");
    printf("then int_res_tor=1\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  if(general_data->timeinfo.int_res_tor==1 && 
                     general_data->timeinfo.int_res_tra == 0){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    printf("Having torsional respa on and intra molecular respa off\n");
    printf("doesn't makes sense.\n");
    printf("The program therefore insists that if int_res_tor=1\n");
    printf("then int_res_tra=1\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  
/*========================================================================*/
/* II.V) Thermostat Consistency Checks                                          */

  if((class->therm_info_class.therm_typ == 2) && general_data->simopts.anneal_opt == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    printf("Annealing with GGMT thermostats not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  
  if((class->therm_info_bead.therm_typ == 2) && general_data->simopts.anneal_opt == 1){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@\n");
    printf("Annealing with GGMT thermostats not available\n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }
  
/*========================================================================*/
/* III) CP consistency checks                                             */

  if(general_data->simopts.cp==1 || general_data->simopts.cp_wave==1 || 
     general_data->simopts.cp_min == 1 || 
                  general_data->simopts.cp_wave_min==1){
    if(cp->cpopts.cp_lda == 1){
      if((strcasecmp(cp->pseudo.vxc_typ,"pz_lsda")==0)||
	 (strcasecmp(cp->pseudo.vxc_typ,"pw_lsda")==0)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("If you are doing a cp calculation using lda you\n");
       printf("cannot use an lsda vxc.\n"); 
       printf("you have selected vxc type %s\n",cp->pseudo.vxc_typ);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
       exit(1);
      } /*endif*/
    }/*endif*/

    if(general_data->simopts.cp_min==1 || 
                   general_data->simopts.cp_wave_min==1 || 
       general_data->simopts.cp_wave_min_pimd==1){
     if(cp->cpopts.cp_ptens_calc == 1){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("There really is not much use in computing the\n");
       printf("pressure tensor under CP minimization.  I will\n"); 
       printf("allow you to continue, but you might want to turn off\n");
       printf("the cp_ptens option.\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
        fflush(stdout);
      } /*endif*/
    }/*endif*/

    if(cp->cpopts.cp_lsda==1){
      if(strcasecmp(cp->pseudo.vxc_typ,"pz_lda")==0 ||
        strcasecmp(cp->pseudo.vxc_typ,"pw_lda")==0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("If you are doing a cp calculation using lsda you\n");
       printf("cannot use an lda vxc, you have selected vxc type %s\n",
              cp->pseudo.vxc_typ);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
      }/*endif*/
    }/*endif*/

    if(cp->cpopts.cp_sic==1 && cp->cpopts.cp_norb > 0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("If you are doing a CP calculation using sic\n");
      printf("you cannot use norb integration\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endif*/

  cp_min_on = general_data->simopts.cp_min 
            + general_data->simopts.cp_wave_min 
            + general_data->simopts.cp_wave_min_pimd;

  if(cp_min_on == 1 && cp->cpopts.cp_norb > 0) {
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Norb minimization not supported yet.  Sorry, but\n");
       printf("technical support is overworked and underpaid.  \n");
       printf("See Dawn, she'll give you the references!!\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
  }/*endif*/

  cp_on = general_data->simopts.cp
         +general_data->simopts.cp_wave
         +general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd
         +general_data->simopts.debug_cp+general_data->simopts.debug_cp_pimd;

  if( (cp_parse->istart_cp == 0) && (cp_on==1)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dude, you are doing a full CP calculation with a \n");
      printf("GEN_WAVE restart?????? NO WAY!!!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*enddif*/

 if(cp_on ==1 && cp->cpopts.cp_ptens_calc != 0){
   if(general_data->cell.iperd == 0 || general_data->cell.iperd == 4){ 
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dude, you are doing a CP calculation with cluster \n");
      printf("boundary conditions. Why are you calculating the\n");
      printf("the pressure tensor as well?\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
  }/*endif*/

  if(cp->cpcoeffs_info.ks_rot_on == 1){
     general_data->filenames.iwrite_kseigs = cp->cpcoeffs_info.n_ks_rot;
  } else {
    general_data->filenames.iwrite_kseigs = general_data->timeinfo.ntime+1;
  }/* endif */

  if( cp->cpcoeffs_info.cp_elf_calc_frq > 0){
     general_data->filenames.iwrite_elf = cp->cpcoeffs_info.cp_elf_calc_frq;
  } else {
    general_data->filenames.iwrite_elf = general_data->timeinfo.ntime+1;
  }/* endif */

  if( (cp->cpcoeffs_info.cp_elf_calc_frq > 0)&&(cp_min_on>0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("NO ELFing during cp-minimization \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if( (cp->cpcoeffs_info.cp_elf_calc_frq > 0)&& 
      (cp->cpopts.cp_dual_grid_opt >= 1)){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("NO ELFing with dual gridding without debugging! \n");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

  if( (cp->cpopts.icheck_perd_size != 1)  && (cp_on + cp_min_on > 0)) {
   if(general_data->cell.iperd <= 2 || general_data->cell.iperd == 4){ 
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("You are doing a CP calculation with reduced boundary\n");
       printf("conditions (< 3) and have elected NOT to check the\n"); 
       printf("inter-particle distances for instances of \n"); 
       printf("particles escaping the box. \n"); 
       printf("ARE YOU SURE THIS IS WHAT YOU WANT TO DO?????? \n"); 
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
   }/*endif*/
  }/*endif*/

  if( (general_data->cell.imov_cp_box == 1)  && (cp_on + cp_min_on > 0)) {
   if(cp->cpopts.cp_dual_grid_opt < 2){ 
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("If you are doing a CP calculation with the small box \n");
      printf("moving, you must use the not_prop option for the  \n");
      printf("cp_dual_grid_opt. Good bye, and have a nice day \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
  }/*endif*/

  if( (cp_on + cp_min_on > 0)) {
   if(cp->cpopts.cp_dual_grid_opt > 2 && cp->cpopts.cp_ptens_calc == 1){ 
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("If you would like to do a dual grid NPT CP calculation \n");
      printf("you will have to code up the pressure tensor \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
  }/*endif*/

 if((cp_on || cp_min_on) && (class->energy_ctrl.iget_pe_real_inter_freq != 1)){
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("For CP calculations, interatomic PE should be \n");
       printf("Calculated every step.  Setting the frequency to 1\n");
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n"); 
       class->energy_ctrl.iget_pe_real_inter_freq = 1;   
 }

  if(cp->cpopts.cp_dual_grid_opt >= 1 && cp->cpopts.cp_hess_calc != 0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dude, if you would like to do a dual grid with the hessian \n");
      printf("you will have to code it up \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if(class->clatoms_info.hess_calc == 1  && cp->cpopts.cp_ptens_calc != 0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Simulataneous calculation of pressure tensor and \n");
      printf("atomic Hessian under CP not currently available\n");
      printf("If you like REALLY need this option, see support staff\n");
      printf("and they will allocate the separate memory space\n");
      printf("required for this.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if( (cp_on + cp_min_on > 1) && (class->interact.dielectric_opt == 1)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dude, the dielectric option must be off for all  \n");
      printf("CP runs \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if( (cp_on + cp_min_on > 1)){
     if(cp->cpscr.cpscr_atom_pme.pme_on  == 1 && 
        cp->cpopts.cp_dual_grid_opt != 2){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Dude, you may only pme when you dual grid the \n");
      printf("not proportional way \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
     }/*endif*/
     if(cp->cpscr.cpscr_atom_pme.pme_on  == 0 && 
        cp->cpopts.cp_dual_grid_opt == 2){
      printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("Dude, why aren't you using pme when you are dual gridding\n");
      printf("the not proportional way. It really is a time saver. \n");
      printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
      fflush(stdout);
     }/*endif*/

  }/*endif*/

  if( (cp->cpopts.cp_dual_grid_opt >= 1) && 
      (general_data->simopts.cp_min== 1) &&
      (general_data->minopts.min_atm_com_fix_opt==0)){
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");
       printf("Dual gridding is great BUT you might want to fix\n");    
       printf("the com of the system using the min_atm_com_fix option\n");
       printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");   
  }/*endif*/

  if((cp->cpopts.cp_isok_opt == 1) && (cp->cpopts.cp_norb == 0)){
      printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Sorry, but if you want to use the CP isokinetic option,\n");
      printf("you need to turn on one of the nonorthogonal (cp_norb)\n");
      printf("options.  Try `full_ortho' if you are not sure what to do\n");
      printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/* endif */


/*========================================================================*/
/* Reduced Periodicity Warnings */

  if( (iperd!=3) && (cp_on + cp_min_on > 0) ){
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("Danger! Periodicity is not 3! The convergence of the\n");
    printf("electronic structure with box size in the non-periodic\n");
    printf("directions must be tested!\n");
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

  if( (iperd>0) && (iperd!=3) ){
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("Danger! Periodicity is not 3 or 0. The box side, L,\n");
    printf("in the non-periodic direction(s) must be in the range\n");
    printf("L>D where D is the relevent distance parameter, \n");
    printf("clusters, D=4R, surfaces D=2T, wires, D=2l. Here,\n");
    printf("R=cluster radius, T=surface thickness, l=wire length\n");
    printf("formed by the atomic positions \n");
    printf("$$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

/*========================================================================*/
/* IV) CP Ensemble checks                                                 */

  if(general_data->simopts.cp_wave==1||
                      general_data->simopts.cp_wave_pimd == 1){
      if((general_data->ensopts.npt_i + general_data->ensopts.npt_f == 1)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("If you are doing a CP calculation moving only the \n");
      printf("wave functions, constant pressure cannnot be employed \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
   }

  if(cp->cpopts.cp_ptens_calc==0 && cp_on ==1){
      if((general_data->ensopts.npt_i + general_data->ensopts.npt_f == 1)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("If you are doing a CP calculation under constant pressure \n");
      printf("You must turn on the CP pressure tensor calculation \n");
      printf("(cp_ptens{on} \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }
   }


/*========================================================================*/
/* V) Velocity sampling checks   */

  class->vel_samp_class.ivel_smpl_on=0;
  if(class->vel_samp_class.nvx_smpl>0){class->vel_samp_class.ivel_smpl_on=1;}
  if(class->vel_samp_class.nvnhc_smpl>0){class->vel_samp_class.ivel_smpl_on=1;}

  class->vel_samp_class.ivel_scale_on=0;
  if(class->vel_samp_class.nvx_scale>0){class->vel_samp_class.ivel_scale_on=1;}


  cp->vel_samp_cp.ivelc_smpl_on=0;
  if(cp->vel_samp_cp.nvc_smpl>0){cp->vel_samp_cp.ivelc_smpl_on=1;}
  if(cp->vel_samp_cp.nvcnhc_smpl>0){cp->vel_samp_cp.ivelc_smpl_on=1;}

  cp->vel_samp_cp.ivelc_scal_on=0;
  if(cp->vel_samp_cp.nvc_scal>0){cp->vel_samp_cp.ivelc_scal_on=1;}

/*========================================================================*/
/* V) Ensemble checks */

  if((general_data->cell.iperd<2) || (general_data->cell.iperd==4) ){
   if((general_data->ensopts.npt_f+general_data->ensopts.npt_i)==1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Constant Pressure is not implemented for\n");
      printf("periodicity in less than 2 dimensions   \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/
  }/*endif*/

  if((general_data->simopts.ann_rate!=1.0)&&
     ((general_data->ensopts.npt_f+general_data->ensopts.npt_i)==1)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Constant Pressure is not implemented using\n");
      printf("the temperature annealing option\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if(general_data->ensopts.npt_f == 1){
    if(general_data->cell.hmat_cons_typ==0){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested flexible constant pressure. \n");
      printf("This only is appropriate for 3D solid systems.\n");
      printf("For example, membranes won't work but protein crystals will.\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
    if(general_data->cell.hmat_cons_typ==2){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested flexible constant pressure\n");
      printf("with the mono_clinic constraint. For example,\n");
      printf("this ensemble will be stable for some membranes systems\n");
      printf("but not others. It is not a good idea for liquids.\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
    if(general_data->cell.hmat_cons_typ==1){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have requested flexible constant pressure\n");
      printf("with the orthorhombic constraint. This ensemble will be\n");
      printf("stable for membranes but not liquids.\n");
      printf("Are you certain this is what you would like to do?\n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    }/*endif*/
  }/*endif*/

/*========================================================================*/
/* VI) Path integral consistency checks */

  if((pimd_on==1)&&(general_data->simopts.pi_beads<=1)){
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
      printf("You have chosen the Path Integral option with  \n");
      printf("only one bead.                                 \n");
      printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      fflush(stdout);
  }/*endif*/

#ifdef NO_PATH_INT_ANNEALING
  if((pimd_on==1)&&(general_data->simopts.ann_rate!=1.0)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Temperature annealing with path integrals not allowed.\n");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
   }/*endif*/
#endif

/*========================================================================*/
/* VII) Bead RESPA checks*/

    p4 = class->clatoms_info.pi_beads_full_ter;
    p3 = class->clatoms_info.pi_beads_res_ter+1;
    p2 = class->clatoms_info.pi_beads_full_tra+1;
    p1 = class->clatoms_info.pi_beads_res_tra+1;
    psum = p4*p3*p2*p1;
    cppimd_on = general_data->simopts.cp_pimd
               +general_data->simopts.cp_wave_pimd
               +general_data->simopts.cp_wave_min_pimd;
    if(cppimd_on == 0){
      if(psum != class->clatoms_info.pi_beads){ 
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Your have improperly set up your RESPA beads   \n");
        printf("for a Path Integral MD simulation.  The product of \n");
        printf("bead_level parameters %d  must equal pi_beads %d. \n",psum,
                class->clatoms_info.pi_beads);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif*/
   class->clatoms_info.ip_lab = (int *)
                cmalloc(class->clatoms_info.pi_beads*sizeof(int))-1;
    ip=0;
    class->clatoms_info.pi_beads_full_ter_wght=0.0;
    class->clatoms_info.pi_beads_res_ter_wght=0.0;
    class->clatoms_info.pi_beads_full_tra_wght=0.0;
    for(i4=1;i4<=p4;i4++){
      for(i3=1;i3<=p3;i3++){
        for(i2=1;i2<=p2;i2++){
          for(i1=1;i1<=p1;i1++){
           ip++;
           if((i3==1)&&(i2==1)&&(i1==1)){class->clatoms_info.ip_lab[ip]=4;}
           if((i3!=1)&&(i2==1)&&(i1==1)){class->clatoms_info.ip_lab[ip]=3;}
           if((i3!=1)&&(i2!=1)&&(i1==1)){class->clatoms_info.ip_lab[ip]=2;}
           if((i3!=1)&&(i2!=1)&&(i1!=1)){class->clatoms_info.ip_lab[ip]=1;}
           if(class->clatoms_info.ip_lab[ip]==4){
             class->clatoms_info.pi_beads_full_ter_wght+=1.0;
             class->clatoms_info.pi_beads_res_ter_wght+=1.0;
             class->clatoms_info.pi_beads_full_tra_wght+=1.0;
           }/*endif*/
           if(class->clatoms_info.ip_lab[ip]==3){
             class->clatoms_info.pi_beads_res_ter_wght+=1.0;
             class->clatoms_info.pi_beads_full_tra_wght+=1.0;
           }/*endif*/
           if(class->clatoms_info.ip_lab[ip]==2){
             class->clatoms_info.pi_beads_full_tra_wght+=1.0;
           }/*endif*/
          }/*endfor*/
        }/*endfor*/
      }/*endfor*/
    }/*endfor*/
    class->clatoms_info.pi_beads_full_ter_wght/=
                     ((double)class->clatoms_info.pi_beads);
    class->clatoms_info.pi_beads_res_ter_wght/=
                     ((double)class->clatoms_info.pi_beads);
    class->clatoms_info.pi_beads_full_tra_wght/=
                     ((double)class->clatoms_info.pi_beads);
/*========================================================================*/
/* VIII) Parallel checks */

   if((cp_on+cp_min_on)==0){
    if(class->communicate.np !=
      class->communicate.np_forc*class->communicate.np_beads){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have not properly divided up processors among \n");
        printf("beads and classical forces: %d neq %d x %d\n",nproc_tot,
                class->communicate.np_forc,class->communicate.np_beads);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }/*endif*/
   }/*endif*/

   if((cp_on+cp_min_on)==1){
    if(class->communicate.np != 
      class->communicate.np_states*class->communicate.np_beads){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("You have not properly divided up processors among \n");
        printf("beads and states: %d neq %d x %d\n",nproc_tot,
               class->communicate.np_states,class->communicate.np_beads);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }
    if((class->communicate.np_forc != 1)&&
       (class->communicate.np_forc!=class->communicate.np_states)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Under CP you can only have np_forc=1 or np_forc=np_state.\n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }
   }

   if(class->communicate.np != 1){
     if(class_parse->istart!=4){
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       printf("Restart types other than restart_all employed\n");
       printf("in parallel may effect comparisons with simulations \n");
       printf("performed in scalar or with different numbers of procs\n");
       printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
     }/*endif*/
     if((cp_on+cp_min_on)==1){
       if(cp_parse->istart_cp!=4){
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
         printf("CP-Restart types other than restart_all employed\n");
         printf("in parallel may effect comparisons with simulations \n");
         printf("performed in scalar or with different numbers of procs\n");
         printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
       }/*endif*/
     }/*endif*/
   }/*endif*/

/*========================================================================*/
/*========================================================================*/
/* IX) PME checks                                                         */

if(general_data->cell.iperd > 0 && cp_on == 0){

   if((general_data->timeinfo.int_res_ter==0) &&
      (class->part_mesh.pme_res_on == 1)){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("PME respa option requires long range forces respa \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
   }/*endif*/

   if(class->part_mesh.pme_res_on == 1 && class->part_mesh.pme_on == 0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("No PME respa option without the PME option \n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
   }/*endif*/
   if((general_data->timeinfo.int_res_ter==1) &&
      (class->part_mesh.pme_res_on == 1)){ 
      if((class->part_mesh.n_interp_res  > class->part_mesh.n_interp)){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       printf("The RESPA PME n_interp parameter > PME n_interp \n");
       printf("This is not allowed\n");      
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);      
      }/*endif*/
      if((class->part_mesh.kmax_pme  < class->part_mesh.kmax_pme_res)){
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         printf("The RESPA PME kmax parameter > PME kmax parameter\n");
         printf("This is not allowed\n");      
         printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
   }/*endif*/

   if((general_data->timeinfo.int_res_ter==1) &&
      (class->part_mesh.pme_res_on == 0) && (class->part_mesh.pme_on == 1) ){
      if(class_parse->kmax_res==class_parse->kmax_ewd){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        printf("Inter respa with PME option on but PME respa opt is off.\n");
        printf("Placing full reciprocal space calculation \n");
        printf("inside the short range respa loop.\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
      }else{
       if(class_parse->kmax_res!=0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("PME option on but PME respa opt is off under LRF respa.\n");
        printf("This is only OK if kmax_res=kmax or kmax_res=0\n");
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
       }/*endif*/
      }/*endif*/
   }/*endif*/

   if(class->part_mesh.pme_res_on == 1 && class->part_mesh.pme_on == 1){
     if(class_parse->kmax_res==class_parse->kmax_ewd){
       if((class->part_mesh.kmax_pme_res == class->part_mesh.kmax_pme) &&
          (class->part_mesh.n_interp_res == class->part_mesh.n_interp)){
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        printf("PME Respa identical to full PME\n");
        printf("Placing full reciprocal space calculation \n");
        printf("inside the short range respa loop.\n");
        printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        class->part_mesh.pme_res_on = 0; 
       }/*endif*/
     }/*endif*/
   }/*endif*/

}/*endif:iperd > 0*/

/*========================================================================*/
    }/*end routine*/ 
/*========================================================================*/
