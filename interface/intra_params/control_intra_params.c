/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_intra_parms.c                        */
/*                                                                          */
/* This subprogram reads in molecular parameter files and sets              */
/* molecular and intramolecular data sets                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_intra_params(double *tot_memory,CLATOMS_INFO *clatoms_info,
                    CLATOMS_POS *clatoms_pos, 
                    GHOST_ATOMS *ghost_atoms,
                    ATOMMAPS *atommaps,BONDED *bonded,
                    FILENAME_PARSE *filename_parse,
                    FREE_PARSE *free_parse,CLASS_PARSE *class_parse,
                    NULL_INTER_PARSE *null_inter_parse,
                    SIMOPTS *simopts, COMMUNICATE *communicate,
                    int isurf_on)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/

/*========================================================================*/
/*     Local Variables                                                    */

  DICT_INTRA dict_intra;
  BUILD_INTRA build_intra;
  RESBOND_PARSE resbond_parse;
  START_INDEX start_index;

  int i,iii;
  int jmol_typ;
  int ifirst,nresidue,mol_or_res;
  int nmol_tot,nres_bond;
  int nres_tot,ncon_tot,nmass_unphys;
  int myid = communicate->myid;
  double cpu1,cpu2;

  char *filename,*fun_key;
  FILE *fp;

  int pi_beads = clatoms_info->pi_beads;

/*=======================================================================*/
/*   I) Start the routine                                                */

  printf("\n");
  PRINT_LINE_STAR;
  printf("Reading molecular parameter files\n");
  PRINT_LINE_DASH;printf("\n");
  cputime(&cpu1);


/*=======================================================================*/
/* II) Initialize/malloc intra stuff                                     */
/*      (init_intra_params.c)                                            */

  filename               = (char *)cmalloc(MAXWORD*sizeof(char));  
  fun_key                 = (char *)cmalloc(MAXWORD*sizeof(char));  

  init_intra_params(clatoms_info,ghost_atoms,atommaps,&build_intra,
                    bonded,null_inter_parse,&resbond_parse,
		    filename_parse);

/*=======================================================================*/
/* III) Set up the dictionaries                                          */
/*      (set_intra_dict.c)                                               */

  ifirst = 1;
  dict_intra.word         = (DICT_WORD *)cmalloc(sizeof(DICT_WORD));
  set_intra_fun_dict(&dict_intra.fun_dict,&dict_intra.num_fun_dict,ifirst);
  set_atm_dict(&dict_intra.atm_dict,&dict_intra.num_atm_dict,ifirst);
  set_intra_dict(&dict_intra.intra_dict,&dict_intra.num_intra_dict,ifirst);
  set_mol_name_dict(&dict_intra.mol_name_dict,&dict_intra.num_mol_name_dict,
                ifirst);
  set_res_name_dict(&dict_intra.res_name_dict,&dict_intra.num_res_name_dict,
                ifirst);
  set_res_def_dict(&dict_intra.res_def_dict,&dict_intra.num_res_def_dict,
                ifirst);
  set_res_bond_dict(&dict_intra.res_bond_dict,&dict_intra.num_res_bond_dict,
                ifirst);

  strcpy(atommaps->atm_typ[1],"HELP");

/*=======================================================================*/
/* V) Zero the list counters                                             */
 
   bonded->bond.npow             = 0;
   bonded->bond.ncon             = 0;
   null_inter_parse->nbond_nul   = 0;
   bonded->bend.npow             = 0;
   bonded->bend.ncon             = 0;
   null_inter_parse->nbend_nul   = 0;
   bonded->tors.npow             = 0;
   bonded->tors.ncon             = 0;
   bonded->tors.nimpr            = 0;
   null_inter_parse->ntors_nul   = 0;
   bonded->onfo.num              = 0;
   null_inter_parse->nonfo_nul   = 0;
   bonded->bend_bnd.num          = 0;
   clatoms_info->natm_tot        = 0;
   atommaps->nres_typ            = 0;
   atommaps->nfreeze             = 0;
   atommaps->natm_typ            = 0;
   ghost_atoms->nghost_tot       = 0;
   ghost_atoms->natm_comp_max    = 0;
   bonded->grp_bond_con.num_21   = 0;
   bonded->grp_bond_con.num_33   = 0;
   bonded->grp_bond_con.num_43   = 0;
   bonded->grp_bond_con.num_23   = 0;
   bonded->grp_bond_con.num_46   = 0;
   bonded->grp_bond_con.ntyp_21  = 0;
   bonded->grp_bond_con.ntyp_33  = 0;
   bonded->grp_bond_con.ntyp_43  = 0;
   bonded->grp_bond_con.ntyp_23  = 0;
   bonded->grp_bond_con.ntyp_46  = 0;
   bonded->grp_bond_watts.num_33 = 0;
   bonded->grp_bond_watts.ntyp_33= 0;

/*=======================================================================*/
/* IV) Loop over molecular parameter files                               */

  for(jmol_typ=1;jmol_typ<=atommaps->nmol_typ;jmol_typ++){

/*-----------------------------------------------------------------------*/
/*  0) Write to the screen                                               */

    printf("**************************************************************\n");
    printf("Reading from molecular parameter file %s\n",
          filename_parse->mol_param_name[jmol_typ]);
    printf("--------------------------------------------------------------\n");
    printf("\n");

/*-----------------------------------------------------------------------*/
/* 1) Store the present list counter values                              */

    start_index.nbond_pow    = bonded->bond.npow;
    start_index.nbond_con    = bonded->bond.ncon;
    start_index.nbond_nul    = null_inter_parse->nbond_nul;
    start_index.nbend_pow    = bonded->bend.npow;
    start_index.nbend_con    = bonded->bend.ncon;
    start_index.nbend_nul    = null_inter_parse->nbend_nul;
    start_index.ntors_pow    = bonded->tors.npow;
    start_index.ntors_con    = bonded->tors.ncon;
    start_index.ntors_nul    = null_inter_parse->ntors_nul;
    start_index.nonfo        = bonded->onfo.num;
    start_index.nonfo_nul    = null_inter_parse->nonfo_nul;
    start_index.nbend_bnd    = bonded->bend_bnd.num;
    start_index.natm         = clatoms_info->natm_tot;
    start_index.nfreeze      = atommaps->nfreeze;
    start_index.nghost_tot   = ghost_atoms->nghost_tot;
    start_index.ngrp_21      = bonded->grp_bond_con.num_21;
    start_index.ngrp_33      = bonded->grp_bond_con.num_33;
    start_index.ngrp_43      = bonded->grp_bond_con.num_43;
    start_index.ngrp_23      = bonded->grp_bond_con.num_23;
    start_index.ngrp_46      = bonded->grp_bond_con.num_46;
    start_index.ngrp_watt_33 = bonded->grp_bond_watts.num_33;
   
/*------------------------------------------------------------------------*/
/*  2) Count and error check the molecular parm file:                     */
/*     (fetch_residue.c)                                                  */

    strcpy(filename,filename_parse->mol_param_name[jmol_typ]);
    nresidue = atommaps->nres_1mol_jmol_typ[jmol_typ];/* spec in mol_set_file*/

    mol_or_res = 1;
    check_parmfile(filename,&(dict_intra.num_fun_dict),&(dict_intra.fun_dict),
                    fun_key,nresidue,&nres_bond,mol_or_res);
                                          /* in fetch_residue.c */
    resbond_parse.nres_bond = nres_bond; 
    resbond_parse.nresidue  = nresidue;
    resbond_parse_realloc(&resbond_parse);
                                          /* in manipulate_res_bonds.c */

/*-----------------------------------------------------------------------*/
/* 3) Read in the molecule name functional keyword: molecule type        */
/*     number of atoms and/or number of residues                         */
/*     (fetch_residue.c)                                                 */


    fetch_molname(filename,&dict_intra,atommaps,
                  fun_key,jmol_typ,nresidue);
/*-----------------------------------------------------------------------*/
/* 4) Read in the molecule residue definitions:                          */
/*       index,natom,parm_file,name                                      */
/*     (fetch_residue.c)                                                 */

    if(nresidue>0){
       fetch_residue_defs(filename,&dict_intra,
                          atommaps,filename_parse,
                          fun_key,nresidue,
                          jmol_typ,&build_intra);

    }else{
      fetch_residue_def0(atommaps,&build_intra,jmol_typ);  
      strcpy(filename_parse->res_param_name[1],
	     filename_parse->mol_param_name[jmol_typ]);
        /* assigns the molecule type to the residue type 
	   when there are no explicit residues           */
    }/*endif*/

/*-----------------------------------------------------------------------*/
/* 5) Read in the molecule resbonds then map them to the residues:       */
/*     Residue types(pairs),residue indices (pairs), bond sites(pairs),  */
/*     bond files(pairs), bond modifier and bond labels read in          */
/*     (fetch_residue.c and manipulate_res_bonds.c)                      */

    if(nres_bond>0){
       fetch_residue_bonds(filename,&dict_intra,fun_key,
                           resbond_parse.resbond_prm,
                           atommaps,nresidue,
                           nres_bond,jmol_typ,pi_beads); 
                    /* in fetch_residue.c */
       map_residue_bonds(&resbond_parse); /* in manipulate_res_bonds.c */
    }else{
       for(i=1;i<=MAX(nresidue,1);i++){
         resbond_parse.nres_bond1_jres[i]=0;resbond_parse.res_bond1_off[i]=0;
         resbond_parse.nres_bond2_jres[i]=0;resbond_parse.res_bond2_off[i]=0;
       }/*endfor*/  
    }/*endif*/
/*-----------------------------------------------------------------------*/
/* 6) Read in the residues: Set the atoms,bonds,bends,torsions,etc.      */
/*                          of each residue. Modificiations of residues  */
/*                          involved in bonds performed.                 */
/*                          Molecules with no residues treated here also */
/*     (control_res_params.c)                                           */


    control_res_params(tot_memory,clatoms_info,ghost_atoms,
                        atommaps,bonded,&resbond_parse,
			&build_intra,
                        filename_parse,free_parse,class_parse,
                        null_inter_parse,filename,
			&dict_intra,fun_key,jmol_typ);

/*------------------------------------------------------------------------*/
/* 7) Bond the residues together and create intramolecular interactions   */
/*    based on the molecular connectivity                                 */
/*    (residue_bond.c)                                                    */

    if(nres_bond>0){
       resbond_parse.ionfo  = class_parse->ionfo_opt[jmol_typ];
       resbond_parse.iconv  = class_parse->ires_bond_conv[jmol_typ];
       residue_bond(clatoms_info,clatoms_pos,
                    atommaps,bonded,&resbond_parse,
		    &build_intra,jmol_typ,
                    class_parse->mol_hydrog_con_opt[jmol_typ]);
    }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.5) Print out progress in intramolecular connectivity lists           */

#ifdef DEBUG
    printf("End Fetch: nbonds=%d,nbends=%d,ntors=%d,nimpr=%d,nonfo=%d\n",
	 bonded->bond.npow,bonded->bend_bnd.num,
	 bonded->tors.npow,bonded->tors.nimpr,bonded->onfo.num);
#endif

/*------------------------------------------------------------------------*/
/* 7.7) Implement the freeze option                                       */

 fetch_freeze(class_parse,atommaps,&build_intra,&start_index,clatoms_info,
              jmol_typ);

/*------------------------------------------------------------------------*/
/* 7.7.5) Check for consistency with zero_com_vel                          */

  if(atommaps->nfreeze > 0 && class_parse->zero_com_vel==1){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Frozen atoms not compatible with zeroing the center of mass\n");
    printf("If you want to freeze atoms, set the zero_com_vel to `no'\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

/*------------------------------------------------------------------------*/
/* 7.8) Implement the hydrog_mass option                                  */

 fetch_hydrog_mass(class_parse,atommaps,&build_intra,&start_index,
                   clatoms_info,jmol_typ);

/*-----------------------------------------------------------------------*/
/* 8) Replicate the molecule now that it has been constructed            */
/*    (replicate_mol.c)                                                    */

   if(atommaps->nmol_jmol_typ[jmol_typ] > 1){
     replicate_mol(clatoms_info,ghost_atoms,atommaps,&build_intra,bonded,
                   null_inter_parse,&start_index,jmol_typ);
   }/*endif*/

/*-----------------------------------------------------------------------*/

 }/*endfor:jmol_typ*/

/*=======================================================================*/
/*  Check for frozen atoms involved in constraints */

  freeze_con_check(atommaps,&(bonded->bond),&(bonded->bend),&(bonded->tors),
                   &(bonded->grp_bond_con));

/*=======================================================================*/
/*  DEBUG */

#ifdef DEBUG
  printf("Vomitting \n");mal_verify(1);
  vomit_intra_list(clatoms_info,ghost_atoms,atommaps,bonded,null_inter_parse);
#endif

/*=======================================================================*/
/*  VI)Tidy up                                                           */

 close_intra_params(clatoms_info,clatoms_pos,ghost_atoms,atommaps,
                    &build_intra,bonded,null_inter_parse,tot_memory,
                    simopts,communicate->np_forc);

/*=======================================================================*/
/*  VIII) Output to screen                                                */

  cputime(&cpu2);
  printf("\n"); 
  PRINT_LINE_DASH;
  printf("Completed reading molecular parameter files %g\n",cpu2-cpu1);
  PRINT_LINE_STAR;printf("\n");

/*========================================================================*/
/* Check for physical masses */

  nmass_unphys = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    if(clatoms_info->mass[i]<900.0){nmass_unphys++;}
  }/*endfor*/

#define NONEXPERT
#ifdef NONEXPERT
  if(nmass_unphys>0){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("There were %d atoms with masses less than 1/2 AMU\n",nmass_unphys);
    printf("This might be OK for experts, but not in general.\n");
    printf("If you are performing cp minimization, use the option\n");
    printf("class_mass_scale_fact in sim_run_def, instead.\n");
    printf("Expert users can disable this error by undefining the\n");
    printf("NONEXPERT ifdef in control_intra_params.c on line 358\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
#endif

/*========================================================================*/
/*  IX) Set the intramolecular potential params                          */

  set_intra_potent(bonded,&build_intra,
               filename_parse->def_intra_name,
               filename_parse->user_intra_name);
  
#ifdef DEBUG
   printf("Vomitting \n");mal_verify(1);
   vomit_intra_potent(bonded,&build_intra);
#endif

/*========================================================================*/
/*  X) Free memory                */                                     
  
  cfree(fun_key);
  cfree(filename);

  cfree(&dict_intra.atm_dict[1]);
  cfree(&dict_intra.intra_dict[1]);
  cfree(&dict_intra.mol_name_dict[1]);
  cfree(&dict_intra.word[0]);
  cfree(&dict_intra.fun_dict[1]);

  cfree(&build_intra.mask_atm[1]);
  cfree(&build_intra.bond_site[1]);
  cfree(&build_intra.index_atm[1]);
  cfree(&build_intra.iatm_ind_chk[1]);
  cfree(&build_intra.cbond_typ_pow[1]);
  cfree(&build_intra.cbond_typ_con[1]);
  cfree(build_intra.cbond_typ_now);
  cfree(&build_intra.cbend_typ_pow[1]);
  cfree(&build_intra.cbend_typ_con[1]);
  cfree(build_intra.cbend_typ_now);
  cfree(&build_intra.ctors_typ_pow[1]);
  cfree(&build_intra.ctors_typ_con[1]);
  cfree(build_intra.ctors_typ_now);
  cfree(&build_intra.confo_typ[1]);
  cfree(build_intra.confo_typ_now); 
  cfree(&build_intra.cbend_bnd_typ[1]);
  cfree(build_intra.cbend_bnd_typ_now);

/*=======================================================================*/
/*  X) Write out synopsis                                               */
  
  PRINT_LINE_STAR;
  printf("System Summary  \n");
  PRINT_LINE_DASH;printf("\n");

  nmol_tot=0;
  nres_tot=0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    nmol_tot+=(atommaps->nmol_jmol_typ)[i];
    nres_tot+=(atommaps->nres_1mol_jmol_typ)[i]*
              (atommaps->nmol_jmol_typ)[i];
  }/*endfor*/
  atommaps->nres_tot = nres_tot;   

  printf("There are %d molecular units      \n",nmol_tot);
  printf("There are %d molecular unit types \n",atommaps->nmol_typ);
  printf("There are %d residue units        \n",nres_tot);
  printf("There are %d residue unit types   \n",atommaps->nres_typ);
  printf("There are %d atoms                \n",clatoms_info->natm_tot);
  printf("There are %d degrees of freedom   \n",clatoms_info->nfree);
  printf("There are %d ghost atoms          \n",ghost_atoms->nghost_tot);
  printf("There are %d atom types           \n",atommaps->natm_typ);
  printf("There are %d charged atoms        \n",clatoms_info->nchrg);
  printf("There are %d power series bonds   \n",bonded->bond.npow);
  printf("There are %d constrained bonds    \n",bonded->bond.ncon);
  printf("There are %d null bonds           \n",null_inter_parse->nbond_nul);
  printf("There are %d power series bends   \n",bonded->bend.npow);
  printf("There are %d constrained bends    \n",bonded->bend.ncon);
  printf("There are %d null bends           \n",null_inter_parse->nbend_nul);
  printf("There are %d Urey-Bradley bends   \n",bonded->bend_bnd.num);
  printf("There are %d power series torsions\n",bonded->tors.npow
                                               -bonded->tors.nimpr);
  printf("There are %d improper torsions    \n",bonded->tors.nimpr);
  printf("There are %d constrained torsions \n",bonded->tors.ncon);
  printf("There are %d null torsions        \n",null_inter_parse->ntors_nul);
  printf("There are %d lj onefours          \n",bonded->onfo.num);
  printf("There are %d null onefours        \n",null_inter_parse->nonfo_nul);
  printf("There are %d free energy bonds    \n",bonded->bond_free.num);
  printf("There are %d free energy bends    \n",bonded->bend_free.num);
  printf("There are %d free energy torsions \n",bonded->tors_free.num);
  printf("There are %d free energy rbar-sig \n",bonded->rbar_sig_free.nfree);
  printf("There are %d 21 group constraints \n",bonded->grp_bond_con.num_21);
  printf("There are %d 23 group constraints \n",bonded->grp_bond_con.num_23);
  printf("There are %d 33 group constraints \n",bonded->grp_bond_con.num_33);
  printf("There are %d 43 group constraints \n",bonded->grp_bond_con.num_43);
  printf("There are %d 46 group constraints \n",bonded->grp_bond_con.num_46);
  printf("There are %d 33 group Watts       \n",bonded->grp_bond_watts.num_33);
  printf("There is  %d surface              \n",isurf_on);
  printf("\n");
  
  PRINT_LINE_DASH;
  printf("End System Summary\n");
  PRINT_LINE_STAR;printf("\n");
  if((simopts->debug+simopts->debug_cp+simopts->debug_pimd)==1){
    printf("Enter an integer ");scanf("%d",&iii);
  }/*endif*/

/*=======================================================================*/
/*   V) If np_forc > 1 and there are non-group constraints, die.         */

#ifdef DEVELOP
  if(((communicate->np_forc)>1)&&
      (bonded->bond.ncon+bonded->bend.ncon+bonded->tors.ncon>0)){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Classical force parallel routine implemented for group\n");
    printf("constraints only.\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
#endif
  
/*=======================================================================*/
/*   V) Assign mall variables  */

/*=======================================================================*/
/*  XI) No atm minimization with constraints                             */

  ncon_tot = bonded->grp_bond_con.num_21
           + bonded->grp_bond_con.num_23
           + bonded->grp_bond_con.num_33
           + bonded->grp_bond_con.num_43
           + bonded->grp_bond_con.num_46
           + bonded->bond.ncon
           + bonded->bend.ncon
           + bonded->tors.ncon;
  if(((simopts->cp_min)==1)&&(ncon_tot>0)){
    printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    printf("Atomic position minimization with constraints under CP \n");
    printf("not implemented\n");
    printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/ 

/*--------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/


#ifdef DEBUG

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void  vomit_intra_potent(BONDED *bonded, BUILD_INTRA *build_intra)

/*==========================================================================*/
{ /*begin routine */
/*==========================================================================*/

int i,iii,j;

/*==========================================================================*/
/* ii) Bonds */

printf("------\n");
printf("pow bonds typ: %d\n ",bonded->bond.ntyp_pow);
printf("------\n");
   for(i=1;i<=bonded->bond.ntyp_pow;i++){
     printf("%d atom1 %s atom2 %s eq %g \n",i,
             build_intra->cbond_typ_pow[i].atm1,
             build_intra->cbond_typ_pow[i].atm2,
             bonded->bond.eq_pow[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con bonds typ: %d\n ",bonded->bond.ntyp_con);
printf("------\n");
  for(i=1;i<=bonded->bond.ntyp_con;i++){
     printf("%d atom1 %s atom2 %s eq %g \n",i,
             build_intra->cbond_typ_con[i].atm1,
             build_intra->cbond_typ_con[i].atm2,
             bonded->bond.eq_con[i]);
  }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Bends */

printf("------\n");
printf("pow bends typ: %d\n ",bonded->bend.ntyp_pow);
printf("------\n");
   for(i=1;i<=bonded->bend.ntyp_pow;i++){
     printf("%d atom1 %s atom2 %s atom3 %s eq %g \n",i,
             build_intra->cbend_typ_pow[i].atm1,
             build_intra->cbend_typ_pow[i].atm2,
             build_intra->cbend_typ_pow[i].atm3,
             bonded->bend.eq_pow[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con bends typ: %d\n ",bonded->bend.ntyp_con);
printf("------\n");
   for(i=1;i<=bonded->bend.ntyp_con;i++){
     printf("%d atom1 %s atom2 %s atom3 %s eq %g \n",i,
             build_intra->cbend_typ_con[i].atm1,
             build_intra->cbend_typ_con[i].atm2,
             build_intra->cbend_typ_con[i].atm3,
             bonded->bend.eq_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Bend_bnds */

printf("------\n");
printf("bend_bnd typ: %d\n",bonded->bend_bnd.ntyp);
printf("------\n");
   for(i=1;i<=bonded->bend_bnd.ntyp;i++){
     printf("%d atom1 %s atom2 %s atom3 %s eq_bo %g eq_be %g \n",i,
             build_intra->cbend_bnd_typ[i].atm1,
             build_intra->cbend_bnd_typ[i].atm2,
             build_intra->cbend_bnd_typ[i].atm3,
             bonded->bend_bnd.eq_bond[i],
             bonded->bend_bnd.eq_bend[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Tors */

printf("------\n");
printf("pow torsions typ: %d\n",bonded->tors.ntyp_pow);
printf("------\n");
   for(i=1;i<=bonded->tors.ntyp_pow;i++){
     printf("%d atom1 %s atom2 %s atom3 %s atom4 %s eq %g \n",i,
             build_intra->ctors_typ_pow[i].atm1,
             build_intra->ctors_typ_pow[i].atm2,
             build_intra->ctors_typ_pow[i].atm3,
             build_intra->ctors_typ_pow[i].atm4,
             bonded->tors.eq_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con torsions: %d\n",bonded->tors.ncon);
printf("------\n");
   for(i=1;i<=bonded->tors.ntyp_con;i++){
     printf("%d atom1 %s atom2 %s atom3 %s atom4 %s eq %g \n",i,
             build_intra->ctors_typ_con[i].atm1,
             build_intra->ctors_typ_con[i].atm2,
             build_intra->ctors_typ_con[i].atm3,
             build_intra->ctors_typ_con[i].atm4,
             bonded->tors.eq_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Onfos */

printf("------\n");
printf("onfos typs: %d\n",bonded->onfo.ntyp);
printf("------\n");
   for(i=1;i<=bonded->onfo.ntyp;i++){
     printf("%d atom1 %s atom2 %s feps %g s6 %g \n",i,
             build_intra->confo_typ[i].atm1,
             build_intra->confo_typ[i].atm2,
             bonded->onfo.feps[i],bonded->onfo.s6[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 21 */

printf("------\n");mal_verify(1);
printf("Grp con 21 typ: %d \n",bonded->grp_bond_con.ntyp_21);
printf("------\n");
   for(i=1;i<=bonded->grp_bond_con.ntyp_21;i++){
     printf("%d atom1 %s atom2 %s eqs %g \n",i,
             build_intra->cgrp_typ_21[i].atm1,
             build_intra->cgrp_typ_21[i].atm2,
             bonded->grp_bond_con.eq_21[1][i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 23 */

printf("------\n");mal_verify(1);
printf("Grp con 23 typ: %d \n",bonded->grp_bond_con.ntyp_23);
printf("------\n");
   for(i=1;i<=bonded->grp_bond_con.ntyp_23;i++){
     printf("%d atom1 %s atom2 %s atom3 %s eqs %g %g\n",i,
             build_intra->cgrp_typ_23[i].atm1,
             build_intra->cgrp_typ_23[i].atm2,
             build_intra->cgrp_typ_23[i].atm3,
             bonded->grp_bond_con.eq_23[1][i],
             bonded->grp_bond_con.eq_23[2][i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 33 */

printf("------\n");mal_verify(1);
printf("Grp con 33 typ: %d \n",bonded->grp_bond_con.ntyp_33);
printf("------\n");
   for(i=1;i<=bonded->grp_bond_con.ntyp_33;i++){
     printf("%d atom1 %s atom2 %s atom3 %s eqs %g %g %g\n",i,
             build_intra->cgrp_typ_33[i].atm1,
             build_intra->cgrp_typ_33[i].atm2,
             build_intra->cgrp_typ_33[i].atm3,
             bonded->grp_bond_con.eq_33[1][i],
             bonded->grp_bond_con.eq_33[2][i],
             bonded->grp_bond_con.eq_33[3][i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 43 */

printf("------\n");mal_verify(1);
printf("Grp con 43 typ: %d \n",bonded->grp_bond_con.ntyp_43);
printf("------\n");
   for(i=1;i<=bonded->grp_bond_con.ntyp_43;i++){
     printf("%d atom1 %s atom2 %s atom3 %s atom4 %s\n",i,
             build_intra->cgrp_typ_43[i].atm1,
             build_intra->cgrp_typ_43[i].atm2,
             build_intra->cgrp_typ_43[i].atm3,
             build_intra->cgrp_typ_43[i].atm4);
     printf("%d eqs %g %g %g \n",i,
             bonded->grp_bond_con.eq_43[1][i],
             bonded->grp_bond_con.eq_43[2][i],
             bonded->grp_bond_con.eq_43[3][i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 46 */

printf("------\n");mal_verify(1);
printf("Grp con 46 typ: %d \n",bonded->grp_bond_con.ntyp_46);
printf("------\n");
   for(i=1;i<=bonded->grp_bond_con.ntyp_46;i++){
     printf("%d atom1 %s atom2 %s atom3 %s atom4 %s\n",i,
             build_intra->cgrp_typ_46[i].atm1,
             build_intra->cgrp_typ_46[i].atm2,
             build_intra->cgrp_typ_46[i].atm3,
             build_intra->cgrp_typ_46[i].atm4);
     printf("%d eqs %g %g %g %g %g %g \n",i,
             bonded->grp_bond_con.eq_46[1][i],
             bonded->grp_bond_con.eq_46[2][i],
             bonded->grp_bond_con.eq_46[3][i],
             bonded->grp_bond_con.eq_46[4][i],
             bonded->grp_bond_con.eq_46[5][i],
             bonded->grp_bond_con.eq_46[6][i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*-----------------------------------------------------------------------*/
}/*end routine */
/*=======================================================================*/

#endif

