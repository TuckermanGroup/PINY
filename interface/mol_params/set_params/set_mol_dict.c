/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_mol_dict.c                               */
/*                                                                          */
/* These subprograms sets the molecular setup dictionaries                  */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_mol_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*  set_molset_fun_dict:set up the molecular setup functional               */
/*                    key word dictionary                                   */ 
/*==========================================================================*/

void set_molset_fun_dict(DICT_WORD *fun_dict[],int *num_fun_dict)

/*==========================================================================*/
{
  int i;
/*=======================================================================*/
/* 0) Malloc the dictionary */

  *num_fun_dict   = 23;
  *fun_dict       = (DICT_WORD *)cmalloc(*num_fun_dict*sizeof(DICT_WORD))-1;  

/*=======================================================================*/
/* I) Assign the users set flags 0 */

  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].iuset = 0;}
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].iflag = 0;}
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].key_type = 1;}

/*=======================================================================*/
/* II) Fill the dictionary with words */
  /*--------------------------------------------------------------------*/ 
  /*  1) ~molecule_def[] */
  strcpy((*fun_dict)[1].keyword,"molecule_def");
  strcpy((*fun_dict)[1].keyarg,"");  
  strcpy((*fun_dict)[1].error_mes,"");
  (*fun_dict)[1].key_type=2; /*specify more than once*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) ~wavefunc_def[] */
  strcpy((*fun_dict)[2].keyword,"wavefunc_def");
  strcpy((*fun_dict)[2].keyarg,"");  
  strcpy((*fun_dict)[2].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  3) ~bond_free_def[] */
  strcpy((*fun_dict)[3].keyword,"bond_free_def");
  strcpy((*fun_dict)[3].keyarg,"");  
  strcpy((*fun_dict)[3].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  4) ~bend_free_def[] */
  strcpy((*fun_dict)[4].keyword,"bend_free_def");
  strcpy((*fun_dict)[4].keyarg,"");  
  strcpy((*fun_dict)[4].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  5) ~tors_free_def[] */
  strcpy((*fun_dict)[5].keyword,"tors_free_def");
  strcpy((*fun_dict)[5].keyarg,"");  
  strcpy((*fun_dict)[5].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  6) ~data_base_def[] */
  strcpy((*fun_dict)[6].keyword,"data_base_def");
  strcpy((*fun_dict)[6].keyarg,"");  
  strcpy((*fun_dict)[6].error_mes,
         "inter_file,vps_file,bond_file,bend_file,tors_file,onfo_file");
  /*-----------------------------------------------------------------------*/ 
  /*  7) ~user_data_base_def[] */
  strcpy((*fun_dict)[7].keyword,"user_data_base_def");
  strcpy((*fun_dict)[7].keyarg,"");  
  strcpy((*fun_dict)[7].error_mes,
    "inter_file,vps_file,bond_file,bend_file,tors_file,onfo_file,surf_file");
  /*-----------------------------------------------------------------------*/ 
  /*  8) ~rbar_sig_free_def[] */
  strcpy((*fun_dict)[8].keyword,"rbar_sig_free_def");
  strcpy((*fun_dict)[8].keyarg,"");  
  strcpy((*fun_dict)[8].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  9) ~surface_def[] */
  strcpy((*fun_dict)[9].keyword,"surface_def");
  strcpy((*fun_dict)[9].keyarg,"");  
  strcpy((*fun_dict)[9].error_mes,"");

  /*-----------------------------------------------------------------------*/ 
  /*  10) ~sim_list_def[ ] */
       strcpy((*fun_dict)[10].error_mes," ");
       strcpy((*fun_dict)[10].keyword,"sim_list_def");
       strcpy((*fun_dict)[10].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  11) ~sim_cp_def[] */
       strcpy((*fun_dict)[11].error_mes," ");
       strcpy((*fun_dict)[11].keyword,"sim_cp_def");
       strcpy((*fun_dict)[11].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  12) ~sim_gen_def[ ] */
       strcpy((*fun_dict)[12].error_mes," ");
       strcpy((*fun_dict)[12].keyword,"sim_gen_def");
       strcpy((*fun_dict)[12].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  13) ~sim_class_PE_def[ ] */
       strcpy((*fun_dict)[13].error_mes," ");
       strcpy((*fun_dict)[13].keyword,"sim_class_PE_def");
       strcpy((*fun_dict)[13].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  14) ~sim_run_def[] */
       strcpy((*fun_dict)[14].error_mes," ");
       strcpy((*fun_dict)[14].keyword,"sim_run_def");
       strcpy((*fun_dict)[14].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  15) ~sim_nhc_def[ ] */
       strcpy((*fun_dict)[15].error_mes," ");
       strcpy((*fun_dict)[15].keyword,"sim_nhc_def");
       strcpy((*fun_dict)[15].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  16) ~sim_vol_def[ ] */
       strcpy((*fun_dict)[16].error_mes," ");
       strcpy((*fun_dict)[16].keyword,"sim_vol_def");
       strcpy((*fun_dict)[16].keyarg," ");
  /*-----------------------------------------------------------------------*/
  /*  17) ~sim_write_def[] */
       strcpy((*fun_dict)[17].error_mes," ");
       strcpy((*fun_dict)[17].keyword,"sim_write_def");
       strcpy((*fun_dict)[17].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  18) ~sim_pimd_def[ ] */
       strcpy((*fun_dict)[18].error_mes," ");
       strcpy((*fun_dict)[18].keyword,"sim_pimd_def");
       strcpy((*fun_dict)[18].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  19) ~sim_velo_corel[ ] */
       strcpy((*fun_dict)[19].error_mes," ");
       strcpy((*fun_dict)[19].keyword,"sim_velo_corel");
       strcpy((*fun_dict)[19].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  20) ~sim_msqd_corel[ ] */
       strcpy((*fun_dict)[20].error_mes," ");
       strcpy((*fun_dict)[20].keyword,"sim_msqd_corel");
       strcpy((*fun_dict)[20].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  21) ~sim_iikt_iso_corel[ ] */
       strcpy((*fun_dict)[21].error_mes," ");
       strcpy((*fun_dict)[21].keyword,"sim_iikt_iso_corel");
       strcpy((*fun_dict)[21].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  22) ~sim_ickt_iso_corel[ ] */
       strcpy((*fun_dict)[22].error_mes," ");
       strcpy((*fun_dict)[22].keyword,"sim_ickt_iso_corel");
       strcpy((*fun_dict)[22].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  23) ~sim_rdf_corel[ ] */
       strcpy((*fun_dict)[23].error_mes," ");
       strcpy((*fun_dict)[23].keyword,"sim_rdf_corel");
       strcpy((*fun_dict)[23].keyarg," ");
  /*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/ 
   }  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_user_base_dict:set up the ~user_data_base_def[] dictionary          */ 
/*==========================================================================*/

void set_user_base_dict(DICT_WORD *user_base_dict[],int *num_user_base_dict)

/*==========================================================================*/
{
  int i;
  /*=======================================================================*/
  /* 0) Malloc the dictionary */

  *num_user_base_dict = 7;
  *user_base_dict=(DICT_WORD*)cmalloc(*num_user_base_dict*sizeof(DICT_WORD))-1;
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_user_base_dict;i++){(*user_base_dict)[i].iuset = 0;}
  for(i=1;i<=*num_user_base_dict;i++){(*user_base_dict)[i].key_type = 0;}
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*---------------------------------------------------------------------*/ 
  /* 1) /inter_file{} */
  strcpy((*user_base_dict)[1].keyword,"inter_file");
  strcpy((*user_base_dict)[1].keyarg,"");
  strcpy((*user_base_dict)[1].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 2) /vps_file{} */
  strcpy((*user_base_dict)[2].keyword,"vps_file");
  strcpy((*user_base_dict)[2].keyarg,"");
  strcpy((*user_base_dict)[2].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 3) /bond_file{} */
  strcpy((*user_base_dict)[3].keyword,"bond_file");
  strcpy((*user_base_dict)[3].keyarg,"");
  strcpy((*user_base_dict)[3].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 4) /bend_file{} */
  strcpy((*user_base_dict)[4].keyword,"bend_file");
  strcpy((*user_base_dict)[4].keyarg,"");
  strcpy((*user_base_dict)[4].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 5) /tors_file{} */
  strcpy((*user_base_dict)[5].keyword,"torsion_file");
  strcpy((*user_base_dict)[5].keyarg,"");
  strcpy((*user_base_dict)[5].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 6) /onefour_file{} */
  strcpy((*user_base_dict)[6].keyword,"onefour_file");
  strcpy((*user_base_dict)[6].keyarg,"");
  strcpy((*user_base_dict)[6].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 7) /surface_file{} */
  strcpy((*user_base_dict)[7].keyword,"surface_file");
  strcpy((*user_base_dict)[7].keyarg,"");
  strcpy((*user_base_dict)[7].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_def_base_dict:set up the ~data_base_def[] dictionary                */ 
/*==========================================================================*/

void set_def_base_dict(DICT_WORD *def_base_dict[],int *num_def_base_dict)

/*==========================================================================*/
{
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  *num_def_base_dict =  7;
  *def_base_dict=(DICT_WORD *)cmalloc(*num_def_base_dict*sizeof(DICT_WORD))-1;
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_def_base_dict;i++){(*def_base_dict)[i].iuset = 0;}
  for(i=1;i<=*num_def_base_dict;i++){(*def_base_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /* 1) /inter_file{} */
  strcpy((*def_base_dict)[1].keyword,"inter_file");
  strcpy((*def_base_dict)[1].keyarg,"pi_md.inter");
  strcpy((*def_base_dict)[1].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 2) /vps_file{} */
  strcpy((*def_base_dict)[2].keyword,"vps_file");
  strcpy((*def_base_dict)[2].keyarg,"pi_md.vps");
  strcpy((*def_base_dict)[2].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 3) /bond_file{} */
  strcpy((*def_base_dict)[3].keyword,"bond_file");
  strcpy((*def_base_dict)[3].keyarg,"pi_md.bond");
  strcpy((*def_base_dict)[3].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 4) /bend_file{} */
  strcpy((*def_base_dict)[4].keyword,"bend_file");
  strcpy((*def_base_dict)[4].keyarg,"pi_md.bend");
  strcpy((*def_base_dict)[4].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 5) /tors_file{} */
  strcpy((*def_base_dict)[5].keyword,"torsion_file");
  strcpy((*def_base_dict)[5].keyarg,"pi_md.tors");
  strcpy((*def_base_dict)[5].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 6) /onefour_file{} */
  strcpy((*def_base_dict)[6].keyword,"onefour_file");
  strcpy((*def_base_dict)[6].keyarg,"pi_md.onfo");
  strcpy((*def_base_dict)[6].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 7) /surface_file{} */
  strcpy((*def_base_dict)[7].keyword,"surface_file");
  strcpy((*def_base_dict)[7].keyarg,"pi_md.surf");
  strcpy((*def_base_dict)[7].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_wave_dict:set up the ~wavefunc_def[] keyword dictionary             */ 
/*==========================================================================*/

void set_wave_dict(DICT_WORD *wave_dict[],int *num_wave_dict,
		   CP_PARSE *cp_parse)

/*==========================================================================*/
{
  int i;
  /*=======================================================================*/
  /* 0) Malloc the dictionary */
  *num_wave_dict = 8;
  *wave_dict      = (DICT_WORD *)cmalloc(*num_wave_dict*sizeof(DICT_WORD))-1;  
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_wave_dict;i++){(*wave_dict)[i].iuset = 0;}
  for(i=1;i<=*num_wave_dict;i++){(*wave_dict)[i].key_type = 0;}
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /* 1) /nstate_up{} */
  strcpy((*wave_dict)[1].keyword,"nstate_up");
  strcpy((*wave_dict)[1].keyarg,"0");
  (*wave_dict)[1].key_type=1; /* must spec */
  strcpy((*wave_dict)[1].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 2) /nstate_dn{} */
  strcpy((*wave_dict)[2].keyword,"nstate_dn");
  strcpy((*wave_dict)[2].keyarg,"0");
  (*wave_dict)[2].key_type=1; /* must spec */
  strcpy((*wave_dict)[2].error_mes,
	 "a number > 0; under LDA nstate_dn=nstate_up");
  /*-----------------------------------------------------------------------*/ 
  /* 3) /cp_nhc_opt{} */
  strcpy((*wave_dict)[3].keyword,"cp_nhc_opt");
  strcpy((*wave_dict)[3].keyarg,"none");
  strcpy((*wave_dict)[3].error_mes,"none,global,glob_st,ind_st,mass_coef");
  /*-----------------------------------------------------------------------*/ 
  /* 4) /cp_tau_nhc{} */
  strcpy((*wave_dict)[4].keyword,"cp_tau_nhc");
  sprintf((*wave_dict)[4].keyarg,"%g",cp_parse->cp_tau_nhc_def);
  strcpy((*wave_dict)[4].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 5) /cp_energy_cut{} */
  strcpy((*wave_dict)[5].keyword,"cp_energy_cut");
  sprintf((*wave_dict)[5].keyarg,"%g",cp_parse->cp_ecut_def);
  strcpy((*wave_dict)[5].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 6) /cp_mass_cut{} */
  strcpy((*wave_dict)[6].keyword,"cp_mass_cut");
  sprintf((*wave_dict)[6].keyarg,"%g",cp_parse->cp_mass_cut_def);
  strcpy((*wave_dict)[6].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 7) /cp_mass_tau{} */
  strcpy((*wave_dict)[7].keyword,"cp_mass_tau");
  sprintf((*wave_dict)[7].keyarg,"%g",cp_parse->cp_mass_tau_def);
  strcpy((*wave_dict)[7].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/ 
  /* 8) /cp_energy_cut_dual_grid{} */
  strcpy((*wave_dict)[8].keyword,"cp_energy_cut_dual_grid");
  sprintf((*wave_dict)[8].keyarg,"%g",cp_parse->cp_ecut_dual_grid_def);
  strcpy((*wave_dict)[8].error_mes,"a number > 0");

  /*-----------------------------------------------------------------------*/ 

}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_mol_dict:set up the ~moleclule_def[] dictionary                     */ 
/*==========================================================================*/

void set_mol_dict(DICT_WORD *mol_dict[],int *num_mol_dict, int iextend,
                  double tau_nhc_def,double t_ext,int ifirst)

/*==========================================================================*/
{
  int i;
/*==========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_mol_dict = 13;
    *mol_dict      = (DICT_WORD *)cmalloc(*num_mol_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_mol_dict;i++){(*mol_dict)[i].iuset = 0;}
  for(i=1;i<=*num_mol_dict;i++){(*mol_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\mol_index{} */
  strcpy((*mol_dict)[1].keyword,"mol_index");
  strcpy((*mol_dict)[1].keyarg,"");
  strcpy((*mol_dict)[1].error_mes,"a number > 0");
  (*mol_dict)[1].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  2)\num_mol{} */
  strcpy((*mol_dict)[2].keyword,"num_mol");
  strcpy((*mol_dict)[2].keyarg,"");
  strcpy((*mol_dict)[2].error_mes,"a number > 0");
  (*mol_dict)[2].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  3)\mol_parm_file{} */
  strcpy((*mol_dict)[3].keyword,"mol_parm_file");
  strcpy((*mol_dict)[3].keyarg,"");
  strcpy((*mol_dict)[3].error_mes,"none");
  (*mol_dict)[3].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  4)\mol_nhc_opt{} */
  strcpy((*mol_dict)[4].keyword,"mol_opt_nhc");
  if(iextend==0)strcpy((*mol_dict)[4].keyarg,"none");
  if(iextend==1)strcpy((*mol_dict)[4].keyarg,"global");
  strcpy((*mol_dict)[4].error_mes,
	 "none,global,glob_mol,ind_mol,res_mol,atm_mol,mass_mol");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\mol_tau_nhc{} */
  strcpy((*mol_dict)[5].keyword,"mol_tau_nhc");
  sprintf((*mol_dict)[5].keyarg,"%g",tau_nhc_def);
  strcpy((*mol_dict)[5].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\mol_text_nhc{} */
  strcpy((*mol_dict)[6].keyword,"mol_text_nhc");
  sprintf((*mol_dict)[6].keyarg,"%g",t_ext);
  strcpy((*mol_dict)[6].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\mol_name{} */
  strcpy((*mol_dict)[7].keyword,"mol_name");
  strcpy((*mol_dict)[7].keyarg,"MOLECULE");
  strcpy((*mol_dict)[7].error_mes,"Molecule name (same as in parm file)");
  (*mol_dict)[7].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  8)\num_residue{} */
  strcpy((*mol_dict)[8].keyword,"num_residue");
  strcpy((*mol_dict)[8].keyarg,"0");
  strcpy((*mol_dict)[8].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\onefour_opt{} */
  strcpy((*mol_dict)[9].keyword,"onefour_opt");
  strcpy((*mol_dict)[9].keyarg,"off");
  strcpy((*mol_dict)[9].error_mes,"on,off");
/*-----------------------------------------------------------------------*/ 
  /*  10)\res_bond_convention{} */
  strcpy((*mol_dict)[10].keyword,"res_bond_convention");
  strcpy((*mol_dict)[10].keyarg,"pi_md");
  strcpy((*mol_dict)[10].error_mes,"pi_md,charmm,sybyl");
/*-----------------------------------------------------------------------*/ 
  /*  11)\mol_freeze_opt{} */
  strcpy((*mol_dict)[11].keyword,"mol_freeze_opt");
  strcpy((*mol_dict)[11].keyarg,"none");
  strcpy((*mol_dict)[11].error_mes,"none,all,heavy_atoms,backbone");
/*-----------------------------------------------------------------------*/ 
  /*  12)\hydrog_mass_opt{} */
  strcpy((*mol_dict)[12].keyword,"hydrog_mass_opt");
  strcpy((*mol_dict)[12].keyarg,"off,1.008");
  strcpy((*mol_dict)[12].error_mes,
                   "[off,all,backbone,sidechain],[number>0]");
/*-----------------------------------------------------------------------*/
  /*  13)\hydrog_con_opt{} */
  strcpy((*mol_dict)[13].keyword,"hydrog_con_opt");
  strcpy((*mol_dict)[13].keyarg,"off");
  strcpy((*mol_dict)[13].error_mes,"off,all,polar (NH,OH,SH)");
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bond_free_dict:set up the ~bond_free_def[] key word dictionary      */ 
/*==========================================================================*/
void set_bond_free_dict(DICT_WORD *bond_free_dict[],int *num_bond_free_dict)
{
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  *num_bond_free_dict = 15;
  *bond_free_dict=(DICT_WORD*)cmalloc(*num_bond_free_dict*sizeof(DICT_WORD))-1;
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_bond_free_dict;i++){(*bond_free_dict)[i].iuset = 0;}
  for(i=1;i<=*num_bond_free_dict;i++){(*bond_free_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom1_moltyp_ind{} */
  strcpy((*bond_free_dict)[1].keyword, "atom1_moltyp_ind");
  strcpy((*bond_free_dict)[1].keyarg,"");
  strcpy((*bond_free_dict)[1].error_mes,"a number > 0");
  (*bond_free_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2_moltyp_ind{} */
  strcpy((*bond_free_dict)[2].keyword, "atom2_moltyp_ind");
  strcpy((*bond_free_dict)[2].keyarg,"");
  strcpy((*bond_free_dict)[2].error_mes,"a number > 0");
  (*bond_free_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom1_mol_ind{}    */
  strcpy((*bond_free_dict)[3].keyword, "atom1_mol_ind");
  strcpy((*bond_free_dict)[3].keyarg,"");
  strcpy((*bond_free_dict)[3].error_mes,"a number > 0");
  (*bond_free_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom2_mol_ind{}    */
  strcpy((*bond_free_dict)[4].keyword, "atom2_mol_ind");
  strcpy((*bond_free_dict)[4].keyarg,"");
  strcpy((*bond_free_dict)[4].error_mes,"a number > 0");
  (*bond_free_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom1_atm_ind{}    */
  strcpy((*bond_free_dict)[5].keyword, "atom1_atm_ind");
  strcpy((*bond_free_dict)[5].keyarg,"");
  strcpy((*bond_free_dict)[5].error_mes,"a number > 0");
  (*bond_free_dict)[5].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\atom2_atm_ind{}     */
  strcpy((*bond_free_dict)[6].keyword,"atom2_atm_ind");
  strcpy((*bond_free_dict)[6].keyarg,"");
  strcpy((*bond_free_dict)[6].error_mes,"a number > 0");
  (*bond_free_dict)[6].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\eq{}               */
  strcpy((*bond_free_dict)[7].keyword,"eq");
  strcpy((*bond_free_dict)[7].keyarg,"");
  strcpy((*bond_free_dict)[7].error_mes,"a number > 0");
  (*bond_free_dict)[7].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 8)\fk{}               */
  strcpy((*bond_free_dict)[8].keyword,"fk");
  strcpy((*bond_free_dict)[8].keyarg,"");
  strcpy((*bond_free_dict)[8].error_mes,"a number > 0");
  (*bond_free_dict)[8].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\rmin_hist{}         */
  strcpy((*bond_free_dict)[9].keyword,"rmin_hist");
  strcpy((*bond_free_dict)[9].keyarg,"");
  strcpy((*bond_free_dict)[9].error_mes,"a number > 0");
  (*bond_free_dict)[9].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\rmax_hist{}         */
  strcpy((*bond_free_dict)[10].keyword,"rmax_hist");
  strcpy((*bond_free_dict)[10].keyarg,"");
  strcpy((*bond_free_dict)[10].error_mes,"a number > 0");
  (*bond_free_dict)[10].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\num_hist{}         */
  strcpy((*bond_free_dict)[11].keyword,"num_hist");
  strcpy((*bond_free_dict)[11].keyarg,"");
  strcpy((*bond_free_dict)[11].error_mes,"a number > 0");
  (*bond_free_dict)[11].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\power{}            */
  strcpy((*bond_free_dict)[12].keyword,"power");
  strcpy((*bond_free_dict)[12].keyarg,"2");
  strcpy((*bond_free_dict)[12].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 13)\hist_file{}        */
  strcpy((*bond_free_dict)[13].keyword,"hist_file");
  strcpy((*bond_free_dict)[13].keyarg,"bond_free_hist.out");
  strcpy((*bond_free_dict)[13].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 14)\atom1_residue_ind{} */
  strcpy((*bond_free_dict)[14].keyword, "atom1_residue_ind");
  strcpy((*bond_free_dict)[14].keyarg,"1");
  strcpy((*bond_free_dict)[14].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 15)\atom2_residue_ind{} */
  strcpy((*bond_free_dict)[15].keyword, "atom2_residue_ind");
  strcpy((*bond_free_dict)[15].keyarg,"1");
  strcpy((*bond_free_dict)[15].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_bend_free_dict:set up the ~bend_free_def[] key word dictionary      */ 
/*==========================================================================*/

void set_bend_free_dict(DICT_WORD *bend_free_dict[],int *num_bend_free_dict)

/*==========================================================================*/
{
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  *num_bend_free_dict = 17;
  *bend_free_dict = (DICT_WORD *)
    cmalloc(*num_bend_free_dict*sizeof(DICT_WORD))-1;
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_bend_free_dict;i++){(*bend_free_dict)[i].iuset = 0;}
  for(i=1;i<=*num_bend_free_dict;i++){(*bend_free_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom1_moltyp_ind{} */
  strcpy((*bend_free_dict)[1].keyword, "atom1_moltyp_ind");
  strcpy((*bend_free_dict)[1].keyarg,"");
  strcpy((*bend_free_dict)[1].error_mes,"a number > 0");
  (*bend_free_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2_moltyp_ind{} */
  strcpy((*bend_free_dict)[2].keyword, "atom2_moltyp_ind");
  strcpy((*bend_free_dict)[2].keyarg,"");
  strcpy((*bend_free_dict)[2].error_mes,"a number > 0");
  (*bend_free_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom3_moltyp_ind{} */
  strcpy((*bend_free_dict)[3].keyword, "atom3_moltyp_ind");
  strcpy((*bend_free_dict)[3].keyarg,"");
  strcpy((*bend_free_dict)[3].error_mes,"a number > 0");
  (*bend_free_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom1_mol_ind{}    */
  strcpy((*bend_free_dict)[4].keyword, "atom1_mol_ind");
  strcpy((*bend_free_dict)[4].keyarg,"");
  strcpy((*bend_free_dict)[4].error_mes,"a number > 0");
  (*bend_free_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom2_mol_ind{}    */
  strcpy((*bend_free_dict)[5].keyword, "atom2_mol_ind");
  strcpy((*bend_free_dict)[5].keyarg,"");
  strcpy((*bend_free_dict)[5].error_mes,"a number > 0");
  (*bend_free_dict)[5].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  6)\atom3_mol_ind{}    */
  strcpy((*bend_free_dict)[6].keyword, "atom3_mol_ind");
  strcpy((*bend_free_dict)[6].keyarg,"");
  strcpy((*bend_free_dict)[6].error_mes,"a number > 0");
  (*bend_free_dict)[6].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  7)\atom1_atm_ind{}    */
  strcpy((*bend_free_dict)[7].keyword, "atom1_atm_ind");
  strcpy((*bend_free_dict)[7].keyarg,"");
  strcpy((*bend_free_dict)[7].error_mes,"a number > 0");
  (*bend_free_dict)[7].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  8)\atom2_atm_ind{}    */
  strcpy((*bend_free_dict)[8].keyword,"atom2_atm_ind");
  strcpy((*bend_free_dict)[8].keyarg,"");
  strcpy((*bend_free_dict)[8].error_mes,"a number > 0");
  (*bend_free_dict)[8].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\atom3_atm_ind{}    */
  strcpy((*bend_free_dict)[9].keyword,"atom3_atm_ind");
  strcpy((*bend_free_dict)[9].keyarg,"");
  strcpy((*bend_free_dict)[9].error_mes,"a number > 0");
  (*bend_free_dict)[9].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\eq{}               */
  strcpy((*bend_free_dict)[10].keyword,"eq");
  strcpy((*bend_free_dict)[10].keyarg,"");
  strcpy((*bend_free_dict)[10].error_mes,"a number > 0");
  (*bend_free_dict)[10].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\fk{}               */
  strcpy((*bend_free_dict)[11].keyword,"fk");
  strcpy((*bend_free_dict)[11].keyarg,"");
  strcpy((*bend_free_dict)[11].error_mes,"a number > 0");
  (*bend_free_dict)[11].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\num_hist{}         */
  strcpy((*bend_free_dict)[12].keyword,"num_hist");
  strcpy((*bend_free_dict)[12].keyarg,"");
  strcpy((*bend_free_dict)[12].error_mes,"a number > 0");
  (*bend_free_dict)[12].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 13)\power{}            */
  strcpy((*bend_free_dict)[13].keyword,"power");
  strcpy((*bend_free_dict)[13].keyarg,"2");
  strcpy((*bend_free_dict)[13].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 14)\hist_file{}        */
  strcpy((*bend_free_dict)[14].keyword,"hist_file");
  strcpy((*bend_free_dict)[14].keyarg,"bend_free_hist.out");
  strcpy((*bend_free_dict)[14].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 15)\atom1_residue_ind{} */
  strcpy((*bend_free_dict)[15].keyword, "atom1_residue_ind");
  strcpy((*bend_free_dict)[15].keyarg,"1");
  strcpy((*bend_free_dict)[15].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 16)\atom2_residue_ind{} */
  strcpy((*bend_free_dict)[16].keyword, "atom2_residue_ind");
  strcpy((*bend_free_dict)[16].keyarg,"1");
  strcpy((*bend_free_dict)[16].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\atom3_residue_ind{} */
  strcpy((*bend_free_dict)[17].keyword, "atom3_residue_ind");
  strcpy((*bend_free_dict)[17].keyarg,"1");
  strcpy((*bend_free_dict)[17].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_tors_free_dict:set up the ~tors_free_def[] key word dictionary      */ 
/*==========================================================================*/

void set_tors_free_dict(DICT_WORD *tors_free_dict[],int *num_tors_free_dict)

/*==========================================================================*/
{/*begin routine*/

  int i,ind;
  /*=========================================================================*/
  /* 0) Malloc the dictionary */
  *num_tors_free_dict = 39;
  *tors_free_dict=(DICT_WORD*)cmalloc(*num_tors_free_dict*sizeof(DICT_WORD))-1;
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_tors_free_dict;i++){(*tors_free_dict)[i].iuset = 0;}
  for(i=1;i<=*num_tors_free_dict;i++){(*tors_free_dict)[i].key_type = 0;}

  /*=========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\atom1.1_moltyp_ind{} */
  strcpy((*tors_free_dict)[1].keyword, "atom1.1_moltyp_ind");
  strcpy((*tors_free_dict)[1].keyarg,"");
  strcpy((*tors_free_dict)[1].error_mes,"a number > 0");
  (*tors_free_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2)\atom2.1_moltyp_ind{} */
  strcpy((*tors_free_dict)[2].keyword, "atom2.1_moltyp_ind");
  strcpy((*tors_free_dict)[2].keyarg,"");
  strcpy((*tors_free_dict)[2].error_mes,"a number > 0");
  (*tors_free_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3)\atom3.1_moltyp_ind{} */
  strcpy((*tors_free_dict)[3].keyword, "atom3.1_moltyp_ind");
  strcpy((*tors_free_dict)[3].keyarg,"");
  strcpy((*tors_free_dict)[3].error_mes,"a number > 0");
  (*tors_free_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4)\atom4.1_moltyp_ind{} */
  strcpy((*tors_free_dict)[4].keyword, "atom4.1_moltyp_ind");
  strcpy((*tors_free_dict)[4].keyarg,"");
  strcpy((*tors_free_dict)[4].error_mes,"a number > 0");
  (*tors_free_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5)\atom1.1_mol_ind{}    */
  strcpy((*tors_free_dict)[5].keyword, "atom1.1_mol_ind");
  strcpy((*tors_free_dict)[5].keyarg,"");
  strcpy((*tors_free_dict)[5].error_mes,"a number > 0");
  (*tors_free_dict)[5].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  6)\atom2.1_mol_ind{}    */
  strcpy((*tors_free_dict)[6].keyword, "atom2.1_mol_ind");
  strcpy((*tors_free_dict)[6].keyarg,"");
  strcpy((*tors_free_dict)[6].error_mes,"a number > 0");
  (*tors_free_dict)[6].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  7)\atom3.1_mol_ind{}    */
  strcpy((*tors_free_dict)[7].keyword, "atom3.1_mol_ind");
  strcpy((*tors_free_dict)[7].keyarg,"");
  strcpy((*tors_free_dict)[7].error_mes,"a number > 0");
  (*tors_free_dict)[7].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  8)\atom4.1_mol_ind{}    */
  strcpy((*tors_free_dict)[8].keyword, "atom4.1_mol_ind");
  strcpy((*tors_free_dict)[8].keyarg,"");
  strcpy((*tors_free_dict)[8].error_mes,"a number > 0");
  (*tors_free_dict)[8].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  9)\atom1.1_atm_ind{}    */
  strcpy((*tors_free_dict)[9].keyword, "atom1.1_atm_ind");
  strcpy((*tors_free_dict)[9].keyarg,"");
  strcpy((*tors_free_dict)[9].error_mes,"a number > 0");
  (*tors_free_dict)[9].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\atom2.1_atm_ind{}    */
  strcpy((*tors_free_dict)[10].keyword,"atom2.1_atm_ind");
  strcpy((*tors_free_dict)[10].keyarg,"");
  strcpy((*tors_free_dict)[10].error_mes,"a number > 0");
  (*tors_free_dict)[10].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\atom3.1_atm_ind{}    */
  strcpy((*tors_free_dict)[11].keyword,"atom3.1_atm_ind");
  strcpy((*tors_free_dict)[11].keyarg,"");
  strcpy((*tors_free_dict)[11].error_mes,"a number > 0");
  (*tors_free_dict)[11].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\atom4.1_atm_ind{}    */
  strcpy((*tors_free_dict)[12].keyword,"atom4.1_atm_ind");
  strcpy((*tors_free_dict)[12].keyarg,"");
  strcpy((*tors_free_dict)[12].error_mes,"a number > 0");
  (*tors_free_dict)[12].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 13)\eq_1{}               */
  strcpy((*tors_free_dict)[13].keyword,"eq_1");
  strcpy((*tors_free_dict)[13].keyarg,"");
  strcpy((*tors_free_dict)[13].error_mes,"a number between -180 and 180");
  (*tors_free_dict)[13].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 14)\fk{}               */
  strcpy((*tors_free_dict)[14].keyword,"fk");
  strcpy((*tors_free_dict)[14].keyarg,"");
  strcpy((*tors_free_dict)[14].error_mes,"a number > 0");
  (*tors_free_dict)[14].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 15)\num_hist{}         */
  strcpy((*tors_free_dict)[15].keyword,"num_hist");
  strcpy((*tors_free_dict)[15].keyarg,"");
  strcpy((*tors_free_dict)[15].error_mes,"a number > 0");
  (*tors_free_dict)[15].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 16)\power{}            */
  strcpy((*tors_free_dict)[16].keyword,"power");
  strcpy((*tors_free_dict)[16].keyarg,"2");
  strcpy((*tors_free_dict)[16].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\hist_file{}        */
  strcpy((*tors_free_dict)[17].keyword,"hist_file");
  strcpy((*tors_free_dict)[17].keyarg,"tors_free_hist.out");
  strcpy((*tors_free_dict)[17].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 18)\atom1.1_residue_ind{} */
  strcpy((*tors_free_dict)[18].keyword, "atom1.1_residue_ind");
  strcpy((*tors_free_dict)[18].keyarg,"1");
  strcpy((*tors_free_dict)[18].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 19)\atom2.1_residue_ind{} */
  strcpy((*tors_free_dict)[19].keyword, "atom2.1_residue_ind");
  strcpy((*tors_free_dict)[19].keyarg,"1");
  strcpy((*tors_free_dict)[19].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 20)\atom3.1_residue_ind{} */
  strcpy((*tors_free_dict)[20].keyword, "atom3.1_residue_ind");
  strcpy((*tors_free_dict)[20].keyarg,"1");
  strcpy((*tors_free_dict)[20].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 21)\atom4.1_residue_ind{} */
  strcpy((*tors_free_dict)[21].keyword, "atom4.1_residue_ind");
  strcpy((*tors_free_dict)[21].keyarg,"1");
  strcpy((*tors_free_dict)[21].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /*  22)\atom1.2_moltyp_ind{} */
  ind = 22;
  strcpy((*tors_free_dict)[ind].keyword, "atom1.2_moltyp_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  23)\atom2.2_moltyp_ind{} */
  ind = 23;
  strcpy((*tors_free_dict)[ind].keyword, "atom2.2_moltyp_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  24)\atom3.2_moltyp_ind{} */
  ind = 24;
  strcpy((*tors_free_dict)[ind].keyword, "atom3.2_moltyp_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  25)\atom4.2_moltyp_ind{} */
  ind = 25;
  strcpy((*tors_free_dict)[ind].keyword, "atom4.2_moltyp_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  26)\atom1.2_mol_ind{}    */
  ind = 26;
  strcpy((*tors_free_dict)[ind].keyword, "atom1.2_mol_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  27)\atom2.2_mol_ind{}    */
  ind = 27;
  strcpy((*tors_free_dict)[ind].keyword, "atom2.2_mol_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  28)\atom3.2_mol_ind{}    */
  ind = 28;
  strcpy((*tors_free_dict)[ind].keyword, "atom3.2_mol_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  29)\atom4.2_mol_ind{}    */
  ind = 29;
  strcpy((*tors_free_dict)[ind].keyword, "atom4.2_mol_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  30)\atom1.2_atm_ind{}    */
  ind = 30;
  strcpy((*tors_free_dict)[ind].keyword, "atom1.2_atm_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 31)\atom2.2_atm_ind{}    */
  ind = 31;
  strcpy((*tors_free_dict)[ind].keyword,"atom2.2_atm_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 32)\atom3.2_atm_ind{}    */
  ind = 32;
  strcpy((*tors_free_dict)[ind].keyword,"atom3.2_atm_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 33)\atom4.2_atm_ind{}    */
  ind = 33;
  strcpy((*tors_free_dict)[ind].keyword,"atom4.2_atm_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"");
  strcpy((*tors_free_dict)[ind].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 34)\atom1.2_residue_ind{} */
  ind = 34;
  strcpy((*tors_free_dict)[ind].keyword, "atom1.2_residue_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"1");
  strcpy((*tors_free_dict)[ind].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 35)\atom2.2_residue_ind{} */
  ind = 35;
  strcpy((*tors_free_dict)[ind].keyword, "atom2.2_residue_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"1");
  strcpy((*tors_free_dict)[ind].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 36)\atom3.2_residue_ind{} */
  ind = 36;
  strcpy((*tors_free_dict)[ind].keyword, "atom3.2_residue_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"1");
  strcpy((*tors_free_dict)[ind].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/ 
  /* 37)\atom4.2_residue_ind{} */
  ind = 37;
  strcpy((*tors_free_dict)[ind].keyword, "atom4.2_residue_ind");
  strcpy((*tors_free_dict)[ind].keyarg,"1");
  strcpy((*tors_free_dict)[ind].error_mes,"a number >= 0");
  /*-----------------------------------------------------------------------*/
  /* 38)\eq_2{}               */
  ind = 38;
  strcpy((*tors_free_dict)[ind].keyword,"eq_2");
  strcpy((*tors_free_dict)[ind].keyarg,"0.0");
  strcpy((*tors_free_dict)[ind].error_mes,"a number between -180 and 180");
  /*-----------------------------------------------------------------------*/
  /* 39)\ntors{}               */
  ind = 39;
  strcpy((*tors_free_dict)[ind].keyword,"ntors");
  strcpy((*tors_free_dict)[ind].keyarg,"1");
  strcpy((*tors_free_dict)[ind].error_mes,"a number between 1 and 2");
  (*tors_free_dict)[ind].key_type = 1;  /* must spec*/

}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_rbar_free_dict:set up the ~rbar_sig_free_def[] key word dictionary  */ 
/*==========================================================================*/

void set_rbar_free_dict(DICT_WORD *rbar_free_dict[],int *num_rbar_free_dict)

/*==========================================================================*/
   {  /*begin routine */
/*==========================================================================*/
  int i;
/*========================================================================*/
/* 0) Malloc the dictionary */

  *num_rbar_free_dict = 13;
  *rbar_free_dict=(DICT_WORD*)cmalloc(*num_rbar_free_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/* I) Assign the users set flags 0 */

  for(i=1;i<=*num_rbar_free_dict;i++){(*rbar_free_dict)[i].iuset = 0;}
  for(i=1;i<=*num_rbar_free_dict;i++){(*rbar_free_dict)[i].key_type = 0;}

/*========================================================================*/
/* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /* 1)\num_bond{}     */
   strcpy((*rbar_free_dict)[1].keyword,"num_bond");
   strcpy((*rbar_free_dict)[1].keyarg,"0");
   strcpy((*rbar_free_dict)[1].error_mes,"a number > 0");
   (*rbar_free_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 2)\index_file{}     */
   strcpy((*rbar_free_dict)[2].keyword,"index_file");
   strcpy((*rbar_free_dict)[2].keyarg,"");
   strcpy((*rbar_free_dict)[2].error_mes,"none");
   (*rbar_free_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 3)\eq_rbar{}               */
   strcpy((*rbar_free_dict)[3].keyword,"eq_rbar");
   strcpy((*rbar_free_dict)[3].keyarg,"0");
   strcpy((*rbar_free_dict)[3].error_mes,"a number >= 0");
   (*rbar_free_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 4)\fk_rbar{}               */
   strcpy((*rbar_free_dict)[4].keyword,"fk_rbar");
   strcpy((*rbar_free_dict)[4].keyarg,"");
   strcpy((*rbar_free_dict)[4].error_mes,"a number >= 0");
   (*rbar_free_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 5)\eq_sig{}               */
   strcpy((*rbar_free_dict)[5].keyword,"eq_sig");
   strcpy((*rbar_free_dict)[5].keyarg,"0");
   strcpy((*rbar_free_dict)[5].error_mes,"a number >= 0");
   (*rbar_free_dict)[5].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 6)\fk_sig{}               */
   strcpy((*rbar_free_dict)[6].keyword,"fk_sig");
   strcpy((*rbar_free_dict)[6].keyarg,"");
   strcpy((*rbar_free_dict)[6].error_mes,"a number >= 0");
   (*rbar_free_dict)[6].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 7)\rmin_hist{}         */
   strcpy((*rbar_free_dict)[7].keyword,"rmin_hist");
   strcpy((*rbar_free_dict)[7].keyarg,"");
   strcpy((*rbar_free_dict)[7].error_mes,"a number > 0");
   (*rbar_free_dict)[7].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 8)\rmax_hist{}         */
   strcpy((*rbar_free_dict)[8].keyword,"rmax_hist");
   strcpy((*rbar_free_dict)[8].keyarg,"");
   strcpy((*rbar_free_dict)[8].error_mes,"a number > 0");
   (*rbar_free_dict)[8].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 9)\num_rbar_hist{}         */
   strcpy((*rbar_free_dict)[9].keyword,"num_rbar_hist");
   strcpy((*rbar_free_dict)[9].keyarg,"");
   strcpy((*rbar_free_dict)[9].error_mes,"a number > 0");
   (*rbar_free_dict)[9].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 10)\smin_hist{}         */
   strcpy((*rbar_free_dict)[10].keyword,"smin_hist");
   strcpy((*rbar_free_dict)[10].keyarg,"");
   strcpy((*rbar_free_dict)[10].error_mes,"a number > 0");
   (*rbar_free_dict)[10].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 11)\smax_hist{}         */
   strcpy((*rbar_free_dict)[11].keyword,"smax_hist");
   strcpy((*rbar_free_dict)[11].keyarg,"");
   strcpy((*rbar_free_dict)[11].error_mes,"a number > 0");
   (*rbar_free_dict)[11].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 12)\num_sig_hist{}         */
   strcpy((*rbar_free_dict)[12].keyword,"num_sig_hist");
   strcpy((*rbar_free_dict)[12].keyarg,"");
   strcpy((*rbar_free_dict)[12].error_mes,"a number > 0");
   (*rbar_free_dict)[12].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /* 13)\hist_file{}        */
   strcpy((*rbar_free_dict)[13].keyword,"hist_file");
   strcpy((*rbar_free_dict)[13].keyarg,"rbar_sig_free_hist.out");
   strcpy((*rbar_free_dict)[13].error_mes,"none");

/*-----------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_surface_dict : set up the ~surface_def[] key word dictionary        */ 
/*==========================================================================*/

void set_surf_dict(DICT_WORD *surface_dict[],int *num_surface_dict)

/*==========================================================================*/
   {  /*begin routine */
/*==========================================================================*/
  int i;
/*==========================================================================*/
/* 0) Malloc the dictionary                                                 */

  *num_surface_dict = 4;
  *surface_dict=(DICT_WORD*)cmalloc(*num_surface_dict*sizeof(DICT_WORD))-1;

/*==========================================================================*/
/* I) Assign the users set flags 0                                          */

  for(i=1;i<=*num_surface_dict;i++){(*surface_dict)[i].iuset = 0;}
  for(i=1;i<=*num_surface_dict;i++){(*surface_dict)[i].key_type = 0;}

/*==========================================================================*/
/* II) Fill the dictionary with words                                       */

  /* 1)\surface_type{}                                                      */
   strcpy((*surface_dict)[1].keyword,"surface_type");
   strcpy((*surface_dict)[1].keyarg,"");
   strcpy((*surface_dict)[1].error_mes,"none");
   (*surface_dict)[1].key_type = 1;                            /* must spec */
/*--------------------------------------------------------------------------*/
  /* 2)\surface_height{}                                                    */
   strcpy((*surface_dict)[2].keyword,"surface_height");
   strcpy((*surface_dict)[2].keyarg,"0.0");
   strcpy((*surface_dict)[2].error_mes,"none");
/*--------------------------------------------------------------------------*/
  /* 3)\num_spline_pts{}                                                    */
   strcpy((*surface_dict)[3].keyword,"num_spline_pts");
   strcpy((*surface_dict)[3].keyarg,"5000");
   strcpy((*surface_dict)[3].error_mes,"a number > 0");
/*--------------------------------------------------------------------------*/
  /* 4)\num_spline_pts{}                                                    */
   strcpy((*surface_dict)[4].keyword,"healing_length");
   strcpy((*surface_dict)[4].keyarg,"1.0");
   strcpy((*surface_dict)[4].error_mes,"a number > 0");
/*--------------------------------------------------------------------------*/
    }/*end routine*/
/*==========================================================================*/

