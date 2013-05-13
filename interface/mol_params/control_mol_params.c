/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_mol_parms.c                          */
/*                                                                          */
/* This subprogram reads in molecular setup input                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#define DEBUG

#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_mol_params(CLASS *class,GENERAL_DATA *general_data,
                        BONDED *bonded,CP *cp,
                        CLASS_PARSE *class_parse,CP_PARSE *cp_parse,
                        FREE_PARSE *free_parse,
                        FILENAME_PARSE *filename_parse)

/*==========================================================================*/
    { /*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int i,iii,jmol_typ,nmol_tot; 
  int nmol_typ,bond_free_num,bend_free_num,tors_free_num,rbar_sig_free_iopt;
  int ifirst,nline,nfun_key;
  int num_user,num_def,cp_on;
  double now_memory;
  char *fun_key;
  int nhist,nhist_bar,nhist_sig,nfree,nsurf;

  FILE *fp; 
  DICT_MOL dict_mol;                        /* Dictionaries and sizes */
  DICT_WORD *word;

  char *molsetname   = filename_parse->molsetname;
  double *text_mol;
  double *text_nhc_mol;
  int    *mol_freeze_opt;
  int    *imol_nhc_opt;
  int    *nmol_jmol_typ;
  double *tot_memory = &(class->tot_memory);
  int iperd          = general_data->cell.iperd;
  int pi_beads       = class->clatoms_info.pi_beads;
  int np_forc        = class->communicate.np_forc;
  int iextend,ipress;
  iextend            = (general_data->ensopts.nvt 
                       +general_data->ensopts.npt_i  
                       +general_data->ensopts.npt_f
                       +general_data->ensopts.nst);
  ipress             = (general_data->ensopts.npt_i
                       +general_data->ensopts.npt_f);

  cp_on  =    general_data->simopts.cp_min 
            + general_data->simopts.cp_wave_min 
            + general_data->simopts.cp_wave_min_pimd
         +general_data->simopts.cp
         +general_data->simopts.cp_wave
         +general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd
         +general_data->simopts.debug_cp+general_data->simopts.debug_cp_pimd;


/*========================================================================*/
/* 0) Output to screen */

  printf("\n");  PRINT_LINE_STAR;
  printf("Reading molecular set up file %s\n",molsetname);
  PRINT_LINE_DASH;printf("\n");

/*========================================================================*/
/* I) Set up dictionaries and malloc temporaries                          */

  ifirst = 1;

  set_molset_fun_dict(&dict_mol.fun_dict,&dict_mol.num_fun_dict);
  set_mol_dict(&dict_mol.mol_dict,&dict_mol.num_mol_dict,
               iextend,class_parse->tau_nhc_def,
               general_data->statepoint.t_ext,ifirst);

  set_wave_dict(&dict_mol.wave_dict,&dict_mol.num_wave_dict,cp_parse);
  set_bond_free_dict(&dict_mol.bond_free_dict,&dict_mol.num_bond_free_dict);
  set_bend_free_dict(&dict_mol.bend_free_dict,&dict_mol.num_bend_free_dict);
  set_tors_free_dict(&dict_mol.tors_free_dict,&dict_mol.num_tors_free_dict); 
  set_rbar_free_dict(&dict_mol.rbar_free_dict,&dict_mol.num_rbar_free_dict); 
  set_user_base_dict(&dict_mol.user_base_dict,&dict_mol.num_user_base_dict);   
  set_def_base_dict(&dict_mol.def_base_dict,&dict_mol.num_def_base_dict);
  set_surf_dict(&dict_mol.surface_dict,&dict_mol.num_surface_dict);
  num_user = dict_mol.num_user_base_dict;
  num_def  = dict_mol.num_def_base_dict;

  word = (DICT_WORD *)cmalloc(sizeof(DICT_WORD));
  fun_key = (char *)cmalloc(MAXWORD*sizeof(char));

/*========================================================================*/
/* II) Read the moleset file and count the molecule types                 */
/*       and the free energy stuff                                        */
  
  nmol_typ            = 0;
  nline               = 0;
  nfun_key            = 0;
  bond_free_num       = 0;
  bend_free_num       = 0;
  tors_free_num       = 0;
  rbar_sig_free_iopt  = 0;
  nsurf               = 0;

  fp = cfopen(molsetname,"r");
  while(get_fun_key_cnt(fp,fun_key,&nline,&nfun_key,molsetname)){
    if(!strcasecmp(fun_key,"molecule_def")) {nmol_typ+=1;}
    if(!strcasecmp(fun_key,"bond_free_def")){bond_free_num+=1;}
    if(!strcasecmp(fun_key,"bend_free_def")){bend_free_num+=1;}
    if(!strcasecmp(fun_key,"tors_free_def")){tors_free_num+=1;}
    if(!strcasecmp(fun_key,"rbar_sig_free_def")){rbar_sig_free_iopt+=1;}
    if(!strcasecmp(fun_key,"surface_def")){nsurf+=1;}
  }/*endwhile*/
  fclose(fp);

  class->atommaps.nmol_typ    = nmol_typ;
  bonded->bond_free.num       = bond_free_num;
  bonded->bend_free.num       = bend_free_num;
  bonded->tors_free.num       = tors_free_num;
  bonded->rbar_sig_free.iopt  = rbar_sig_free_iopt ;
  bonded->rbar_sig_free.nfree = 0;

  if((bond_free_num+bend_free_num+tors_free_num+rbar_sig_free_iopt)>1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Only one free energy definition permitted\n");
      printf("in set up file %s \n",molsetname);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

#ifdef JUNK
  if( (nsurf>0) && (iperd != 2)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Surface potentials permitted ONLY in systems with\n");
      printf("2D periodicity in set up file  %s \n",molsetname);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/
#endif

  class->surface.isurf_on = nsurf;

/*========================================================================*/
/* III) Malloc molecular data storage and initialize                      */
  
  now_memory = (double)((nmol_typ)*(sizeof(int)*2))*1.e-06;
  (*tot_memory) += now_memory;
  
  printf("Allocating molecular memory: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,(*tot_memory));
  
  class->atommaps.nmol_jmol_typ  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class->atommaps.mol_typ        = (NAME *)cmalloc(nmol_typ*sizeof(NAME))-1;
  class->atommaps.nres_1mol_jmol_typ  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class->atommaps.jres_jmol_typ_strt  = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class->atommaps.icons_jres_jmol_typ = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class_parse->ionfo_opt         = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class_parse->ires_bond_conv    = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class_parse->tau_nhc_mol       = (double*)cmalloc(nmol_typ*sizeof(double))-1;
  class_parse->text_nhc_mol      = (double*)cmalloc(nmol_typ*sizeof(double))-1;
  class_parse->imol_nhc_opt      = (int *)cmalloc(nmol_typ*sizeof(int))-1;  
  class_parse->mol_freeze_opt    = (int *)cmalloc(nmol_typ*sizeof(int))-1;  
  class_parse->mol_hydrog_mass_opt = (int *)cmalloc(nmol_typ*sizeof(int))-1;
  class_parse->mol_hydrog_mass_val = 
                               (double *)cmalloc(nmol_typ*sizeof(double))-1;
  class_parse->mol_hydrog_con_opt= (int *)cmalloc(nmol_typ*sizeof(int))-1;
  filename_parse->mol_param_name = (NAME *)cmalloc(nmol_typ*sizeof(NAME))-1;
  filename_parse->user_intra_name= (NAME *)cmalloc(num_user*sizeof(NAME))-1;
  filename_parse->def_intra_name = (NAME *)cmalloc(num_def*sizeof(NAME))-1;
  if(bond_free_num>0){
    bonded->bond_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    free_parse->imoltyp_bond_free= (int *) cmalloc(2*sizeof(int))-1;
    free_parse->imol_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
    free_parse->ires_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
    free_parse->iatm_bond_free   = (int *) cmalloc(2*sizeof(int))-1;
  }
  if(bend_free_num>0){
    bonded->bend_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    free_parse->imoltyp_bend_free= (int *) cmalloc(3*sizeof(int))-1;
    free_parse->imol_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
    free_parse->ires_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
    free_parse->iatm_bend_free   = (int *) cmalloc(3*sizeof(int))-1;
  }
  if(tors_free_num>0){
    bonded->tors_free.file       = (char *)cmalloc(MAXWORD*sizeof(char));
    free_parse->imoltyp_tors_free= (int *) cmalloc(8*sizeof(int))-1;
    free_parse->imol_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
    free_parse->ires_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
    free_parse->iatm_tors_free   = (int *) cmalloc(8*sizeof(int))-1;
  }

  if(rbar_sig_free_iopt>0){
    bonded->rbar_sig_free.file   = (char *)cmalloc(MAXWORD*sizeof(char));
  }
  class->clatoms_info.text_mol   = (double*)cmalloc(nmol_typ*sizeof(double))-1;

  text_mol       = class->clatoms_info.text_mol;
  text_nhc_mol   = class_parse->text_nhc_mol;
  mol_freeze_opt = class_parse->mol_freeze_opt;
  imol_nhc_opt   = class_parse->imol_nhc_opt;

/*========================================================================*/
/* IV) Get molecular/CP setup information */

  control_set_mol_params(&(class->atommaps),&(cp->cpopts),
                         &(cp->cpcoeffs_info),cp_parse,class_parse,
                         bonded,&(class->surface),
                         filename_parse,free_parse,
                         &dict_mol,word,
                         fun_key,&nfun_key,iextend,
                         general_data->statepoint.t_ext,
                         ifirst,pi_beads);
  if(tors_free_num==1){
    tors_free_num = bonded->tors_free.num; /* changed in above routine */
  }/*endif*/

  for(i=1; i<= nmol_typ; i++){
    text_mol[i] = text_nhc_mol[i];   
  }/*endfor*/

/*========================================================================*/
/* V) Free some memory and malloc some other                              */

  if(bond_free_num>0){
    nhist = (bonded->bond_free.nhist);
    bonded->bond_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/

  if(bend_free_num>0){
    nhist = (bonded->bend_free.nhist);
    bonded->bend_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/

  if(tors_free_num==1){
    nhist = (bonded->tors_free.nhist);
    bonded->tors_free.hist = (double *)cmalloc(nhist*sizeof(double))-1;
  }/*endif*/
  if(tors_free_num==2){
    nhist = (bonded->tors_free.nhist);
    bonded->tors_free.hist_2d = cmall_mat(1,nhist,1,nhist);
  }/*endif*/

  if(rbar_sig_free_iopt>0){
    nhist_bar = (bonded->rbar_sig_free.nhist_bar);
    nhist_sig = (bonded->rbar_sig_free.nhist_sig);
    nfree     = (bonded->rbar_sig_free.nfree);
    bonded->rbar_sig_free.hist    = cmall_mat(1,nhist_bar,1,nhist_sig);
    bonded->rbar_sig_free.hist_rn = cmall_mat(1,nfree,1,nhist_bar);
  }/*endif*/

  cfree(fun_key);
  cfree(word);
  cfree(&dict_mol.fun_dict[1]);
  cfree(&dict_mol.mol_dict[1]);
  cfree(&dict_mol.wave_dict[1]);
  cfree(&dict_mol.bond_free_dict[1]);
  cfree(&dict_mol.bend_free_dict[1]);
  cfree(&dict_mol.tors_free_dict[1]);
  cfree(&dict_mol.rbar_free_dict[1]);
  cfree(&dict_mol.user_base_dict[1]);
  cfree(&dict_mol.def_base_dict[1]);

  for(jmol_typ=1;jmol_typ<=nmol_typ;jmol_typ++){

    if((mol_freeze_opt[jmol_typ] > 0) &&  (ipress > 0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("No freezing allowed with NPT_I and NPT_F       \n");
      printf("Ensembles                                      \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/

    if((mol_freeze_opt[jmol_typ]==1)&&(imol_nhc_opt[jmol_typ] > 0)){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Do not use thermostating when freezing the     \n");
      printf("whole molecule.(Use none option)               \n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/

  }/*endfor*/

/*========================================================================*/
/* Error check */
 
  nmol_tot = 0;
  nmol_jmol_typ  = class->atommaps.nmol_jmol_typ;
  for(i=1;i<=nmol_typ;i++){nmol_tot+= nmol_jmol_typ[i];}

  if(nmol_tot < np_forc){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");   
      printf("Number of force level processors %d is greater than\n",np_forc);
      printf("the total number of molecules %d\n",nmol_tot);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

  if(cp_on==1 && cp->cpopts.cp_dual_grid_opt >= 1){
   if(cp_parse->cp_ecut_dual_grid > cp_parse->cp_ecut){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");   
      printf("The small dense grid cutoff is less than the large sparse");
      printf("grid cutoff. This might work, but I doubut it\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }
  }

/*========================================================================*/
/* VI) Done                                                               */

  printf("The class contains %d molecule types\n",nmol_typ);
  printf("\n"); 

  PRINT_LINE_DASH;
  printf("Completed reading molecular set up file %s\n",molsetname);
  PRINT_LINE_STAR;
  printf("\n");

/*------------------------------------------------------------------------*/
   }/*end routine*/
/*===========================================================================*/









