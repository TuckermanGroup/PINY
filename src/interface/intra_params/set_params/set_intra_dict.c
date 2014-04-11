/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_intra_dict.c                             */
/*                                                                          */
/* These subprograms set the intra parameter dictionaries                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_fun_dict:set up the intra parameter functional               */
/*                    key word dictionary                                   */ 
/*==========================================================================*/

void set_intra_fun_dict(DICT_WORD *fun_dict[],int *num_fun_dict,int ifirst)

/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*==========================================================================*/
/*       Local Variables */
  int i;
/*=======================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_fun_dict=15;
    *fun_dict   = (DICT_WORD *)cmalloc(*num_fun_dict*sizeof(DICT_WORD))-1;  
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].iuset = 0;}
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].key_type = 1;}
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1) ~molecule_name_def[] */
  strcpy((*fun_dict)[1].keyword,"molecule_name_def");
  strcpy((*fun_dict)[1].keyarg,"");  
  strcpy((*fun_dict)[1].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  2) ~atom_def[] */
  strcpy((*fun_dict)[2].keyword,"atom_def");
  strcpy((*fun_dict)[2].keyarg,"");  
  strcpy((*fun_dict)[2].error_mes,"");
  (*fun_dict)[2].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  3) ~bond_def[] */
  strcpy((*fun_dict)[3].keyword,"bond_def");
  strcpy((*fun_dict)[3].keyarg,"");  
  strcpy((*fun_dict)[3].error_mes,"");
  (*fun_dict)[3].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  4) ~bend_def[] */
  strcpy((*fun_dict)[4].keyword,"bend_def");
  strcpy((*fun_dict)[4].keyarg,"");  
  strcpy((*fun_dict)[4].error_mes,"");
  (*fun_dict)[4].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  5) ~tors_def[] */
  strcpy((*fun_dict)[5].keyword,"torsion_def");
  strcpy((*fun_dict)[5].keyarg,"");  
  strcpy((*fun_dict)[5].error_mes,"");
  (*fun_dict)[5].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  6) ~onfo_def[]     */
  strcpy((*fun_dict)[6].keyword,"onefour_def");
  strcpy((*fun_dict)[6].keyarg,"");  
  strcpy((*fun_dict)[6].error_mes,"");
  (*fun_dict)[6].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  7) ~residue_morph[]    */
  i = 7;
  strcpy((*fun_dict)[i].keyword,"residue_morph");
  strcpy((*fun_dict)[i].keyarg,"");  
  strcpy((*fun_dict)[i].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  8) ~atom_destroy[] */
  strcpy((*fun_dict)[8].keyword,"atom_destroy");
  strcpy((*fun_dict)[8].keyarg,"");  
  strcpy((*fun_dict)[8].error_mes,"");
  (*fun_dict)[8].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  9) ~atom_create[]  */
  strcpy((*fun_dict)[9].keyword,"atom_create");
  strcpy((*fun_dict)[9].keyarg,"");  
  strcpy((*fun_dict)[9].error_mes,"");
  (*fun_dict)[9].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 10) ~atom_morph[]   */
  strcpy((*fun_dict)[10].keyword,"atom_morph");
  strcpy((*fun_dict)[10].keyarg,"");  
  strcpy((*fun_dict)[10].error_mes,"");
  (*fun_dict)[10].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  11) ~bend_bnd_def[] */
  strcpy((*fun_dict)[11].keyword,"bend_bnd_def");
  strcpy((*fun_dict)[11].keyarg,"");  
  strcpy((*fun_dict)[11].error_mes,"");
  (*fun_dict)[11].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  12) ~residue_def[] */
  strcpy((*fun_dict)[12].keyword,"residue_def");
  strcpy((*fun_dict)[12].keyarg,""); 
  strcpy((*fun_dict)[12].error_mes,"");
  (*fun_dict)[12].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  13) ~residue_bond_def[] */
  strcpy((*fun_dict)[13].keyword,"residue_bond_def");
  strcpy((*fun_dict)[13].keyarg,"");  
  strcpy((*fun_dict)[13].error_mes,"");
  (*fun_dict)[13].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /*  14) ~residue_name_def[] */
  i = 14;
  strcpy((*fun_dict)[i].keyword,"residue_name_def");
  strcpy((*fun_dict)[i].keyarg,"");  
  strcpy((*fun_dict)[i].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  15) ~grp_bond_def[] */
  i = 15;
  strcpy((*fun_dict)[i].keyword,"grp_bond_def");
  strcpy((*fun_dict)[i].keyarg,"");  
  strcpy((*fun_dict)[i].error_mes,"");
  (*fun_dict)[i].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_atm_dict:set up the atom  key word dictionary                       */ 
/*==========================================================================*/

void set_atm_dict(DICT_WORD *atm_dict[],int *num_atm_dict, int ifirst)

/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*==========================================================================*/
/*       Local Variables */
  int i,ioff,joff;
/*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_atm_dict=17+NCOEF_GHOST_MAX;
    *atm_dict    = (DICT_WORD *)cmalloc(*num_atm_dict*sizeof(DICT_WORD))-1;  
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_atm_dict;i++){(*atm_dict)[i].iuset = 0;}
  for(i=1;i<=*num_atm_dict;i++){(*atm_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1) /atom_typ{} */
  strcpy((*atm_dict)[1].keyword,"atom_typ");
  strcpy((*atm_dict)[1].keyarg,"");
  strcpy((*atm_dict)[1].error_mes,"");
  (*atm_dict)[1].key_type = 1;  /* must spec*/
  /*------------------------------------------------------------------------*/ 
  /*  2) /atom_ind{} */
  strcpy((*atm_dict)[2].keyword,"atom_ind");
  strcpy((*atm_dict)[2].keyarg,"");
  strcpy((*atm_dict)[2].error_mes,"a number >0 <= natom");
  (*atm_dict)[2].key_type = 1;  /* must spec*/
  /*------------------------------------------------------------------------*/ 
  /*  3) /mass{} */
  strcpy((*atm_dict)[3].keyword,"mass");
  strcpy((*atm_dict)[3].keyarg,"");
  strcpy((*atm_dict)[3].error_mes,"a number =>0");
  (*atm_dict)[3].key_type = 1;  /* must spec*/
  /*------------------------------------------------------------------------*/ 
  /*  4) /charge{} */
  strcpy((*atm_dict)[4].keyword,"charge");
  strcpy((*atm_dict)[4].keyarg,"0.0");
  strcpy((*atm_dict)[4].error_mes,"");
  /*------------------------------------------------------------------------*/ 
  /*  5) /alpha_pol{} */
  strcpy((*atm_dict)[5].keyword,"alpha_pol");
  strcpy((*atm_dict)[5].keyarg,"0.0");
  strcpy((*atm_dict)[5].error_mes,"a number =>0");
  /*------------------------------------------------------------------------*/ 
  /*  6) /b_neut{} */
  strcpy((*atm_dict)[6].keyword,"b_neut");
  strcpy((*atm_dict)[6].keyarg,"0.0");
  strcpy((*atm_dict)[6].error_mes,"a number =>0");
  /*------------------------------------------------------------------------*/ 
  /*  7) /valence{} */
  strcpy((*atm_dict)[7].keyword,"valence");
  strcpy((*atm_dict)[7].keyarg,"1");
  strcpy((*atm_dict)[7].error_mes,"a number > 0  ");
  /*------------------------------------------------------------------------*/ 
  /*  8) /improper_def{} */
  strcpy((*atm_dict)[8].keyword,"improper_def");
  strcpy((*atm_dict)[8].keyarg,"0,0,0,0");
  strcpy((*atm_dict)[8].error_mes,"four integers expected like i1,i2,i3,i4");
  /*------------------------------------------------------------------------*/ 
  /*  9) /bond_site_1{} */
  strcpy((*atm_dict)[9].keyword,"bond_site_1");
  strcpy((*atm_dict)[9].keyarg,"null,-1,-1");
  strcpy((*atm_dict)[9].error_mes,
	 "a bondsite (char), branch (int), atoms from branch");
  /*------------------------------------------------------------------------*/ 
  /*  10) /bond_site_2{} */
  strcpy((*atm_dict)[10].keyword,"bond_site_2");
  strcpy((*atm_dict)[10].keyarg,"null,-1,-1");
  strcpy((*atm_dict)[10].error_mes,
	 "a bondsite (char), branch (int), atoms from branch");
  /*------------------------------------------------------------------------*/ 
  /*  11) /bond_site_3{} */
  strcpy((*atm_dict)[11].keyword,"bond_site_3");
  strcpy((*atm_dict)[11].keyarg,"null,-1,-1");
  strcpy((*atm_dict)[11].error_mes,
	 "a bondsite (char), branch (int), atoms from branch");
  /*------------------------------------------------------------------------*/ 
  /*  12) /bond_site_4{} */
  strcpy((*atm_dict)[12].keyword,"bond_site_4");
  strcpy((*atm_dict)[12].keyarg,"null,-1,-1");
  strcpy((*atm_dict)[12].error_mes,
	 "a bondsite (char), branch (int), atoms from branch");
  /*------------------------------------------------------------------------*/ 
  /*  13) /cp_valence_up{} */
  strcpy((*atm_dict)[13].keyword,"cp_valence_up");
  strcpy((*atm_dict)[13].keyarg,"0");
  strcpy((*atm_dict)[13].error_mes,"a number >= 0");
  /*------------------------------------------------------------------------*/ 
  /*  14) /cp_valence_dn{} */
  strcpy((*atm_dict)[14].keyword,"cp_valence_dn");
  strcpy((*atm_dict)[14].keyarg,"0");
  strcpy((*atm_dict)[14].error_mes,"a number >= 0");
  /*------------------------------------------------------------------------*/ 
  /*  15) /cp_atom{} */
  strcpy((*atm_dict)[15].keyword,"cp_atom");
  strcpy((*atm_dict)[15].keyarg,"no");
  strcpy((*atm_dict)[15].error_mes,"yes,no");
  /*------------------------------------------------------------------------*/ 
  /*  16 ... ) /def_ghost1{} */
  ioff = 15;
  joff = (int) '0';
  for(i=1;i<=NCOEF_GHOST_MAX;i++){
   strcpy((*atm_dict)[(i+ioff)].keyword,"def_ghost");
   (*atm_dict)[(i+ioff)].keyword[9] = (char) (i+joff);
   strcpy((*atm_dict)[(i+ioff)].keyarg,"0,0");
   strcpy((*atm_dict)[(i+ioff)].error_mes,"an index, a coefficient");
  }
  /*------------------------------------------------------------------------*/ 
  /*  16+NCOEF_GHOST_MAX)  /label{} */
  strcpy((*atm_dict)[(16+NCOEF_GHOST_MAX)].keyword,"label");
  strcpy((*atm_dict)[(16+NCOEF_GHOST_MAX)].keyarg,"standard");
  strcpy((*atm_dict)[(16+NCOEF_GHOST_MAX)].error_mes,
	 "backbone, sidechain, standard");
  /*------------------------------------------------------------------------*/ 
  /*  17+NCOEF_GHOST_MAX)  /pdb_typ{} */
  strcpy((*atm_dict)[(17+NCOEF_GHOST_MAX)].keyword,"pdb_typ");
  strcpy((*atm_dict)[(17+NCOEF_GHOST_MAX)].keyarg,"");
  strcpy((*atm_dict)[(17+NCOEF_GHOST_MAX)].error_mes,"");
/*==========================================================================*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_intra_dict:set up intramolecular interaction                        */
/*                    key word dictionary                                   */ 
/*==========================================================================*/

void set_intra_dict(DICT_WORD *intra_dict[],int *num_intra_dict, int ifirst)

/*==========================================================================*/
/*       Begin Routine */
   {/*begin routine*/
/*==========================================================================*/
/*       Local Variables */
  int i;
/*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_intra_dict=12;
    *intra_dict = (DICT_WORD *)cmalloc(*num_intra_dict*sizeof(DICT_WORD))-1;
    }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_intra_dict;i++){(*intra_dict)[i].iuset = 0;}
  for(i=1;i<=*num_intra_dict;i++){(*intra_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*------------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*intra_dict)[1].keyword,"atom1");
  strcpy((*intra_dict)[1].keyarg,"0");
  strcpy((*intra_dict)[1].error_mes,"a number >0 <=natom");
  (*intra_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */
  strcpy((*intra_dict)[2].keyword,"atom2");
  strcpy((*intra_dict)[2].keyarg,"0");
  strcpy((*intra_dict)[2].error_mes,"a number >0 <=natom");
  (*intra_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) \atom3{} */
  strcpy((*intra_dict)[3].keyword,"atom3");
  strcpy((*intra_dict)[3].keyarg,"0");
  strcpy((*intra_dict)[3].error_mes,"a number >0 <=natom");
  (*intra_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) \atom4{} */
  strcpy((*intra_dict)[4].keyword,"atom4");
  strcpy((*intra_dict)[4].keyarg,"0");
  strcpy((*intra_dict)[4].error_mes,"a number >0 <=natom");
  (*intra_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5) \modifier{} */
  strcpy((*intra_dict)[5].keyword,"modifier");
  strcpy((*intra_dict)[5].keyarg,"on");
  strcpy((*intra_dict)[5].error_mes,"con,on,off");
  /*---------------------------------------------------------------------*/ 
  /*  6) \label1{} */
  strcpy((*intra_dict)[6].keyword,"label1");
  strcpy((*intra_dict)[6].keyarg,"");
  strcpy((*intra_dict)[6].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /*  7) \label2{} */
  strcpy((*intra_dict)[7].keyword,"label2");
  strcpy((*intra_dict)[7].keyarg,"");
  strcpy((*intra_dict)[7].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /*  8) \label3{} */
  strcpy((*intra_dict)[8].keyword,"label3");
  strcpy((*intra_dict)[8].keyarg,"");
  strcpy((*intra_dict)[8].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /*  9) \label4{} */
  strcpy((*intra_dict)[9].keyword,"label4");
  strcpy((*intra_dict)[9].keyarg,"");
  strcpy((*intra_dict)[9].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /* 10) \label5{} */
  strcpy((*intra_dict)[10].keyword,"label5");
  strcpy((*intra_dict)[10].keyarg,"");
  strcpy((*intra_dict)[10].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /* 11) \label6{} */
  strcpy((*intra_dict)[11].keyword,"label6");
  strcpy((*intra_dict)[11].keyarg,"");
  strcpy((*intra_dict)[11].error_mes,"");
  /*---------------------------------------------------------------------*/ 
  /*  12) \group_typ{} */
  strcpy((*intra_dict)[12].keyword,"group_typ");
  strcpy((*intra_dict)[12].keyarg,"");
  strcpy((*intra_dict)[12].error_mes,"3x3,2x3,4x6,watts_3x3");
  /*---------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_mol_name_dict:set up the atom  key word dictionary                  */ 
/*==========================================================================*/

void set_mol_name_dict(DICT_WORD *mol_name_dict[],int *num_mol_name_dict,
                        int ifirst)
/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*========================================================================*/
/*       Local Variables */
  int i;
/*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_mol_name_dict=3;
    *mol_name_dict=(DICT_WORD*)cmalloc(*num_mol_name_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_mol_name_dict;i++){(*mol_name_dict)[i].iuset   = 0;}
  for(i=1;i<=*num_mol_name_dict;i++){(*mol_name_dict)[i].key_type = 0;}
  /*=========================================================================*/
  /* II) Fill the dictionary with words */
  /*-------------------------------------------------------------------------*/
  /* 1) \molecule_name{} */
  strcpy((*mol_name_dict)[1].keyword,"molecule_name");
  strcpy((*mol_name_dict)[1].keyarg,"MOLECULE");
  strcpy((*mol_name_dict)[1].error_mes,"the same as specified in molset file");
  (*mol_name_dict)[1].key_type = 1;  /* must spec */
  /*-------------------------------------------------------------------------*/
  /* 2) \nresidue{} */
  strcpy((*mol_name_dict)[2].keyword,"nresidue");
  strcpy((*mol_name_dict)[2].keyarg,"0");
  strcpy((*mol_name_dict)[2].error_mes,
	 "a number > 0 : = the same as specified in the molset file");
  /*-----------------------------------------------------------------------*/ 
  /* 3) \natom{} */
  strcpy((*mol_name_dict)[3].keyword,"natom");
  strcpy((*mol_name_dict)[3].keyarg,"0");
  strcpy((*mol_name_dict)[3].error_mes,
	 "a number > 0 : = the same as specified in the molset file");
  (*mol_name_dict)[3].key_type = 1;  /* must spec */
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_res_name_dict:set up the residue name key word dictionary          */ 
/*==========================================================================*/

void set_res_name_dict(DICT_WORD *res_name_dict[],int *num_res_name_dict,
                        int ifirst)

/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*========================================================================*/
/*       Local Variables */
  int i;
/*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_res_name_dict=3;
    *res_name_dict=(DICT_WORD*)cmalloc(*num_res_name_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_res_name_dict;i++){(*res_name_dict)[i].iuset   = 0;}
  for(i=1;i<=*num_res_name_dict;i++){(*res_name_dict)[i].key_type = 0;}
  /*=========================================================================*/
  /* II) Fill the dictionary with words */
  /*-------------------------------------------------------------------------*/
  /* 1) \res_name{} */
  strcpy((*res_name_dict)[1].keyword,"residue_name");
  strcpy((*res_name_dict)[1].keyarg,"RESIDUE");
  strcpy((*res_name_dict)[1].error_mes,
             "the same as specified in molecule parameter file");
  (*res_name_dict)[1].key_type = 1;  /* must spec */
  /*-------------------------------------------------------------------------*/
 /* 2) \nresidue{} */
  strcpy((*res_name_dict)[2].keyword,"nresidue");
  strcpy((*res_name_dict)[2].keyarg,"0");
  strcpy((*res_name_dict)[2].error_mes,
             "a number > 0 or same as specified in molecular parm file");
  (*res_name_dict)[2].key_type = 0;  /* must spec */
  /*-------------------------------------------------------------------------*/
 /* 3) \natom{} */
  strcpy((*res_name_dict)[3].keyword,"natom");
  strcpy((*res_name_dict)[3].keyarg,"0");
  strcpy((*res_name_dict)[3].error_mes,
             "a number > 0 or same as specified in molecular parm file");
  (*res_name_dict)[3].key_type = 1;  /* must spec */
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_res_def_dict:set up the ~res_def[] dictionary                       */ 
/*==========================================================================*/

void set_res_def_dict(DICT_WORD *res_def_dict[],int *num_res_def_dict,
		      int ifirst)

/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*========================================================================*/
/*       Local Variables */
  int i;
/*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
   *num_res_def_dict = 5;
   *res_def_dict = (DICT_WORD *)cmalloc(*num_res_def_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_res_def_dict;i++){(*res_def_dict)[i].iuset = 0;}
  for(i=1;i<=*num_res_def_dict;i++){(*res_def_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\residue_index{} */
  strcpy((*res_def_dict)[1].keyword,"residue_index");
  strcpy((*res_def_dict)[1].keyarg,"");
  strcpy((*res_def_dict)[1].error_mes,"a number > 0");
  (*res_def_dict)[1].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  2)\natom{} */
  strcpy((*res_def_dict)[2].keyword,"natom");
  strcpy((*res_def_dict)[2].keyarg,"");
  strcpy((*res_def_dict)[2].error_mes,"a number > 0");
  (*res_def_dict)[2].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  3)\residue_parm_file{} */
  strcpy((*res_def_dict)[3].keyword,"residue_parm_file");
  strcpy((*res_def_dict)[3].keyarg,"");
  strcpy((*res_def_dict)[3].error_mes,"none");
  (*res_def_dict)[3].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  4)\residue_name{} */
  strcpy((*res_def_dict)[4].keyword,"residue_name");
  strcpy((*res_def_dict)[4].keyarg,"RESIDUE");
  strcpy((*res_def_dict)[4].error_mes,"");
  (*res_def_dict)[4].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  5)\residue_fix_file{} */
  strcpy((*res_def_dict)[5].keyword,"residue_fix_file");
  strcpy((*res_def_dict)[5].keyarg,"\0");
  strcpy((*res_def_dict)[5].error_mes,"none");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_res_bond_dict:set up the ~res_bond[] key word dictionary            */ 
/*==========================================================================*/

void set_res_bond_dict(DICT_WORD *res_bond_dict[],int *num_res_bond_dict,
                       int ifirst)

/*==========================================================================*/
/*       Begin Routine */
  {/*begin routine*/
/*========================================================================*/
/*       Local Variables */
  int i;
/*========================================================================*/
/* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_res_bond_dict =  10;
    *res_bond_dict = (DICT_WORD *)
    cmalloc(*num_res_bond_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_res_bond_dict;i++){(*res_bond_dict)[i].iuset = 0;}
  for(i=1;i<=*num_res_bond_dict;i++){(*res_bond_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1)\res1_typ{} */
  strcpy((*res_bond_dict)[1].keyword,"res1_typ");
  strcpy((*res_bond_dict)[1].keyarg,"");
  strcpy((*res_bond_dict)[1].error_mes,"none");
  (*res_bond_dict)[1].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  2)\res2_typ{} */
  strcpy((*res_bond_dict)[2].keyword,"res2_typ");
  strcpy((*res_bond_dict)[2].keyarg,"");
  strcpy((*res_bond_dict)[2].error_mes,"none");
  (*res_bond_dict)[2].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  3)\res1_index{} */
  strcpy((*res_bond_dict)[3].keyword,"res1_index");
  strcpy((*res_bond_dict)[3].keyarg,"");
  strcpy((*res_bond_dict)[3].error_mes,"a number > 0");
  (*res_bond_dict)[3].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  4)\res2_index{} */
  strcpy((*res_bond_dict)[4].keyword,"res2_index");
  strcpy((*res_bond_dict)[4].keyarg,"");
  strcpy((*res_bond_dict)[4].error_mes,"a number > 0");
  (*res_bond_dict)[4].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  5)\res1_bond_site{} */
  strcpy((*res_bond_dict)[5].keyword,"res1_bond_site");
  strcpy((*res_bond_dict)[5].keyarg,"");
  strcpy((*res_bond_dict)[5].error_mes,"a number > 0");
  (*res_bond_dict)[5].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  6)\res2_bond_site{} */
  strcpy((*res_bond_dict)[6].keyword,"res2_bond_site");
  strcpy((*res_bond_dict)[6].keyarg,"");
  strcpy((*res_bond_dict)[6].error_mes,"a number > 0");
  (*res_bond_dict)[6].key_type=1; /* must spec */
  /*-----------------------------------------------------------------------*/ 
  /*  7)\res1_bondfile{} */
  strcpy((*res_bond_dict)[7].keyword,"res1_bondfile");
  strcpy((*res_bond_dict)[7].keyarg,"");
  strcpy((*res_bond_dict)[7].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\res2_bondfile{} */
  strcpy((*res_bond_dict)[8].keyword,"res2_bondfile");
  strcpy((*res_bond_dict)[8].keyarg,"");
  strcpy((*res_bond_dict)[8].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\modifier{} */
  strcpy((*res_bond_dict)[9].keyword,"modifier");
  strcpy((*res_bond_dict)[9].keyarg,"on");
  strcpy((*res_bond_dict)[9].error_mes,"on,con,null");
  /*-----------------------------------------------------------------------*/ 
  /*  10)\label{} */
  strcpy((*res_bond_dict)[10].keyword,"label");
  strcpy((*res_bond_dict)[10].keyarg,"");
  strcpy((*res_bond_dict)[10].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
}  /*end routine*/
/*==========================================================================*/









