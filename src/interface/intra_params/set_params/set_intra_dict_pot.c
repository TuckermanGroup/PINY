/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_intra_dict.c                             */
/*                                                                          */
/* These subprograms set the molecular parameter dictionaries               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*               Includes:                                                  */

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potbond_dict:bond potential key word dictionary                     */
/*==========================================================================*/
void set_potbond_dict(DICT_WORD *bond_dict[],int *num_bond_dict,int ifirst)
{/*begin routine*/
  int i;
  /*=======================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_bond_dict=16;
    *bond_dict = (DICT_WORD *)cmalloc(*num_bond_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*=======================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_bond_dict;i++){(*bond_dict)[i].iuset = 0;}
  for(i=1;i<=*num_bond_dict;i++){(*bond_dict)[i].key_type = 0;}
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*----------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*bond_dict)[1].keyword,"atom1");
  strcpy((*bond_dict)[1].keyarg,"");
  strcpy((*bond_dict)[1].error_mes,"");
  (*bond_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */ 
  strcpy((*bond_dict)[2].keyword,"atom2");
  strcpy((*bond_dict)[2].keyarg,"");
  strcpy((*bond_dict)[2].error_mes,"");
  (*bond_dict)[2].key_type = 1;  /* must spec*/
  /*----------------------------------------------------------------------*/ 
  /*  3) \pot_type{} */
  strcpy((*bond_dict)[3].keyword,"pot_type");
  strcpy((*bond_dict)[3].keyarg,"");
  strcpy((*bond_dict)[3].error_mes,"harmonic,power-series,morse,null");
  (*bond_dict)[3].key_type = 1;  /* must spec*/
  /*----------------------------------------------------------------------*/ 
  /*  4) \fk{} */
  strcpy((*bond_dict)[4].keyword,"fk");
  strcpy((*bond_dict)[4].keyarg,"0.0");
  strcpy((*bond_dict)[4].error_mes,"a number >= 0");
  /*----------------------------------------------------------------------*/ 
  /*  5) \eq{} */
  strcpy((*bond_dict)[5].keyword,"eq");
  strcpy((*bond_dict)[5].keyarg,"0.0");
  strcpy((*bond_dict)[5].error_mes,"a number");
  /*----------------------------------------------------------------------*/ 
  /*  16) \eq_res{} */
  strcpy((*bond_dict)[16].keyword,"eq_res");
  strcpy((*bond_dict)[16].keyarg,"0.0");
  strcpy((*bond_dict)[16].error_mes,"a number");
  /*----------------------------------------------------------------------*/ 
  /*  6) \p0{} */
  strcpy((*bond_dict)[6].keyword,"p0");
  strcpy((*bond_dict)[6].keyarg,"0.0");
  strcpy((*bond_dict)[6].error_mes,"a number");
  /*----------------------------------------------------------------------*/ 
  /*  7) \p1{} */
  strcpy((*bond_dict)[7].keyword,"p1");
  strcpy((*bond_dict)[7].keyarg,"0.0");
  strcpy((*bond_dict)[7].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  8) \p2{} */
  strcpy((*bond_dict)[8].keyword,"p2");
  strcpy((*bond_dict)[8].keyarg,"0.0");
  strcpy((*bond_dict)[8].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  9) \p3{} */
  strcpy((*bond_dict)[9].keyword,"p3");
  strcpy((*bond_dict)[9].keyarg,"0.0");
  strcpy((*bond_dict)[9].error_mes,"a number");
  /*----------------------------------------------------------------------*/ 
  /* 10) \p4{} */
  strcpy((*bond_dict)[10].keyword,"p4");
  strcpy((*bond_dict)[10].keyarg,"0.0");
  strcpy((*bond_dict)[10].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 11) \p5{} */
  strcpy((*bond_dict)[11].keyword,"p5");
  strcpy((*bond_dict)[11].keyarg,"0.0");
  strcpy((*bond_dict)[11].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 12) \p6{} */
  strcpy((*bond_dict)[12].keyword,"p6");
  strcpy((*bond_dict)[12].keyarg,"0.0");
  strcpy((*bond_dict)[12].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 13) \alpha} */
  strcpy((*bond_dict)[13].keyword,"alpha");
  strcpy((*bond_dict)[13].keyarg,"0.0");
  strcpy((*bond_dict)[13].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 14) \d0{} */
  strcpy((*bond_dict)[14].keyword,"d0");
  strcpy((*bond_dict)[14].keyarg,"0.0");
  strcpy((*bond_dict)[14].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 15) \label{} */
  strcpy((*bond_dict)[15].keyword,"label");
  strcpy((*bond_dict)[15].keyarg,"");
  strcpy((*bond_dict)[15].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*  set_potbend_dict:bend potential key word dictionary                     */
/*==========================================================================*/
void set_potbend_dict(DICT_WORD *bend_dict[],int *num_bend_dict,int ifirst)
{/*begin routine*/
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
     *num_bend_dict=14;
     *bend_dict = (DICT_WORD *)cmalloc(*num_bend_dict*sizeof(DICT_WORD))-1;
   }/*endif*/
  /*=======================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_bend_dict;i++){(*bend_dict)[i].iuset = 0;}
  for(i=1;i<=*num_bend_dict;i++){(*bend_dict)[i].key_type = 0;}
  /*=========================================================================*/
  /* II) Fill the dictionary with words */
  /*----------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*bend_dict)[1].keyword,"atom1");
  strcpy((*bend_dict)[1].keyarg,"");
  strcpy((*bend_dict)[1].error_mes,"");
  (*bend_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */ 
  strcpy((*bend_dict)[2].keyword,"atom2");
  strcpy((*bend_dict)[2].keyarg,"");
  strcpy((*bend_dict)[2].error_mes,"");
  (*bend_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) \atom3{} */ 
  strcpy((*bend_dict)[3].keyword,"atom3");
  strcpy((*bend_dict)[3].keyarg,"");
  strcpy((*bend_dict)[3].error_mes,"");
  (*bend_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) \pot_type{} */
  strcpy((*bend_dict)[4].keyword,"pot_type");
  strcpy((*bend_dict)[4].keyarg,"");
  strcpy((*bend_dict)[4].error_mes,"harmonic,power-series,freq_series,null");
  (*bend_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5) \fk{} */
  strcpy((*bend_dict)[5].keyword,"fk");
  strcpy((*bend_dict)[5].keyarg,"0.0");
  strcpy((*bend_dict)[5].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  6) \eq{} */
  strcpy((*bend_dict)[6].keyword,"eq");
  strcpy((*bend_dict)[6].keyarg,"0.0");
  strcpy((*bend_dict)[6].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  7) \p0{} */
  strcpy((*bend_dict)[7].keyword,"p0");
  strcpy((*bend_dict)[7].keyarg,"0.0");
  strcpy((*bend_dict)[7].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  8) \p1{} */
  strcpy((*bend_dict)[8].keyword,"p1");
  strcpy((*bend_dict)[8].keyarg,"0.0");
  strcpy((*bend_dict)[8].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  9) \p2{} */
  strcpy((*bend_dict)[9].keyword,"p2");
  strcpy((*bend_dict)[9].keyarg,"0.0");
  strcpy((*bend_dict)[9].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 10) \p3{} */
  strcpy((*bend_dict)[10].keyword,"p3");
  strcpy((*bend_dict)[10].keyarg,"0.0");
  strcpy((*bend_dict)[10].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 11) \p4{} */
  strcpy((*bend_dict)[11].keyword,"p4");
  strcpy((*bend_dict)[11].keyarg,"0.0");
  strcpy((*bend_dict)[11].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 12) \p5{} */
  strcpy((*bend_dict)[12].keyword,"p5");
  strcpy((*bend_dict)[12].keyarg,"0.0");
  strcpy((*bend_dict)[12].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 13) \p6{} */
  strcpy((*bend_dict)[13].keyword,"p6");
  strcpy((*bend_dict)[13].keyarg,"0.0");
  strcpy((*bend_dict)[13].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 14) \label{} */
  strcpy((*bend_dict)[14].keyword,"label");
  strcpy((*bend_dict)[14].keyarg,"");
  strcpy((*bend_dict)[14].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_pottors_dict:tors potential key word dictionary                     */
/*==========================================================================*/
void set_pottors_dict(DICT_WORD *tors_dict[],int *num_tors_dict,int ifirst)

{/*begin routine*/
  int i;
  /*=======================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
    *num_tors_dict=28;
    *tors_dict = (DICT_WORD *)cmalloc(*num_tors_dict*sizeof(DICT_WORD))-1;
  }/*endif*/
  /*=======================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_tors_dict;i++){(*tors_dict)[i].iuset = 0;}
  for(i=1;i<=*num_tors_dict;i++){(*tors_dict)[i].key_type = 0;}

  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*----------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*tors_dict)[1].keyword,"atom1");
  strcpy((*tors_dict)[1].keyarg,"");
  strcpy((*tors_dict)[1].error_mes,"");
  (*tors_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */ 
  strcpy((*tors_dict)[2].keyword,"atom2");
  strcpy((*tors_dict)[2].keyarg,"");
  strcpy((*tors_dict)[2].error_mes,"");
  (*tors_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) \atom3{} */ 
  strcpy((*tors_dict)[3].keyword,"atom3");
  strcpy((*tors_dict)[3].keyarg,"");
  strcpy((*tors_dict)[3].error_mes,"");
  (*tors_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) \atom4{} */ 
  strcpy((*tors_dict)[4].keyword,"atom4");
  strcpy((*tors_dict)[4].keyarg,"");
  strcpy((*tors_dict)[4].error_mes,"");
  (*tors_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5) \pot_type{} */
  strcpy((*tors_dict)[5].keyword,"pot_type");
  strcpy((*tors_dict)[5].keyarg,"");
  strcpy((*tors_dict)[5].error_mes,"harmonic,power-series,freq-series,null");
  (*tors_dict)[5].key_type = 1;  /* must spec*/
  /*----------------------------------------------------------------------*/ 
  /*  6) \fk{} */
  strcpy((*tors_dict)[6].keyword,"fk");
  strcpy((*tors_dict)[6].keyarg,"0.0");
  strcpy((*tors_dict)[6].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  7) \eq{} */
  strcpy((*tors_dict)[7].keyword,"eq");
  strcpy((*tors_dict)[7].keyarg,"0.0");
  strcpy((*tors_dict)[7].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  8) \p0{} */
  strcpy((*tors_dict)[8].keyword,"p0");
  strcpy((*tors_dict)[8].keyarg,"0.0");
  strcpy((*tors_dict)[8].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  9) \p1{} */
  strcpy((*tors_dict)[9].keyword,"p1");
  strcpy((*tors_dict)[9].keyarg,"0.0");
  strcpy((*tors_dict)[9].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  10) \p2{} */
  strcpy((*tors_dict)[10].keyword,"p2");
  strcpy((*tors_dict)[10].keyarg,"0.0");
  strcpy((*tors_dict)[10].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 11) \p3{} */
  strcpy((*tors_dict)[11].keyword,"p3");
  strcpy((*tors_dict)[11].keyarg,"0.0");
  strcpy((*tors_dict)[11].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 12) \p4{} */
  strcpy((*tors_dict)[12].keyword,"p4");
  strcpy((*tors_dict)[12].keyarg,"0.0");
  strcpy((*tors_dict)[12].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 13) \p5{} */
  strcpy((*tors_dict)[13].keyword,"p5");
  strcpy((*tors_dict)[13].keyarg,"0.0");
  strcpy((*tors_dict)[13].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 14) \p6{} */
  strcpy((*tors_dict)[14].keyword,"p6");
  strcpy((*tors_dict)[14].keyarg,"0.0");
  strcpy((*tors_dict)[14].error_mes,"a number");
  /*----------------------------------------------------------------------*/ 
  /* 15) \nfreq{} */
  strcpy((*tors_dict)[15].keyword,"nfreq");
  strcpy((*tors_dict)[15].keyarg,"0");
  strcpy((*tors_dict)[15].error_mes,"a number>0");
  /*-----------------------------------------------------------------------*/ 
  /* 16)  \a1{} */
  strcpy((*tors_dict)[16].keyword,"a1");
  strcpy((*tors_dict)[16].keyarg,"0.0");
  strcpy((*tors_dict)[16].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 17)  \c1{} */
  strcpy((*tors_dict)[17].keyword,"c1");
  strcpy((*tors_dict)[17].keyarg,"0.0");
  strcpy((*tors_dict)[17].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 18)  \d1{} */
  strcpy((*tors_dict)[18].keyword,"d1");
  strcpy((*tors_dict)[18].keyarg,"0.0");
  strcpy((*tors_dict)[18].error_mes,"0 or 180");
  /*-----------------------------------------------------------------------*/ 
  /* 19)  \a2{} */
  strcpy((*tors_dict)[19].keyword,"a2");
  strcpy((*tors_dict)[19].keyarg,"0.0");
  strcpy((*tors_dict)[19].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 20)  \c2{} */
  strcpy((*tors_dict)[20].keyword,"c2");
  strcpy((*tors_dict)[20].keyarg,"0.0");
  strcpy((*tors_dict)[20].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 21)  \d2{} */
  strcpy((*tors_dict)[21].keyword,"d2" );
  strcpy((*tors_dict)[21].keyarg,"0.0");
  strcpy((*tors_dict)[21].error_mes,"0 or 180");
  /*-----------------------------------------------------------------------*/ 
  /* 22)  \a3{} */
  strcpy((*tors_dict)[22].keyword,"a3");
  strcpy((*tors_dict)[22].keyarg,"0.0");
  strcpy((*tors_dict)[22].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 23)  \c3{} */
  strcpy((*tors_dict)[23].keyword,"c3");
  strcpy((*tors_dict)[23].keyarg,"0.0");
  strcpy((*tors_dict)[23].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 24)  \d3{} */
  strcpy((*tors_dict)[24].keyword,"d3");
  strcpy((*tors_dict)[24].keyarg,"0.0");
  strcpy((*tors_dict)[24].error_mes,"0 or 180");
  /*-----------------------------------------------------------------------*/ 
  /* 25) \label{} */
  strcpy((*tors_dict)[25].keyword,"label");
  strcpy((*tors_dict)[25].keyarg,"");
  strcpy((*tors_dict)[25].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /* 26)  \a4{} */
  strcpy((*tors_dict)[26].keyword,"a4");
  strcpy((*tors_dict)[26].keyarg,"0.0");
  strcpy((*tors_dict)[26].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 27)  \c4{} */
  strcpy((*tors_dict)[27].keyword,"c4");
  strcpy((*tors_dict)[27].keyarg,"0.0");
  strcpy((*tors_dict)[27].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 28)  \d4{} */
  strcpy((*tors_dict)[28].keyword,"d4");
  strcpy((*tors_dict)[28].keyarg,"0.0");
  strcpy((*tors_dict)[28].error_mes,"0 or 180");
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potonfo_dict:onfo potential key word dictionary                     */

void set_potonfo_dict(DICT_WORD *onfo_dict[],int *num_onfo_dict,int ifirst)

{/*begin routine*/
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
     *num_onfo_dict=7;
     *onfo_dict = (DICT_WORD *)cmalloc(*num_onfo_dict*sizeof(DICT_WORD))-1;
   }/*endif*/
  /*=======================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_onfo_dict;i++){(*onfo_dict)[i].iuset = 0;}
  for(i=1;i<=*num_onfo_dict;i++){(*onfo_dict)[i].key_type = 0;}
  /*=======================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*onfo_dict)[1].keyword,"atom1");
  strcpy((*onfo_dict)[1].keyarg,"");
  strcpy((*onfo_dict)[1].error_mes,"");
  (*onfo_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */ 
  strcpy((*onfo_dict)[2].keyword,"atom2");
  strcpy((*onfo_dict)[2].keyarg,"");
  strcpy((*onfo_dict)[2].error_mes,"");
  (*onfo_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) \pot_type{} */
  strcpy((*onfo_dict)[3].keyword,"pot_type");
  strcpy((*onfo_dict)[3].keyarg,"");
  strcpy((*onfo_dict)[3].error_mes,"lennard-jones,null");
  (*onfo_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) \sig{} */
  strcpy((*onfo_dict)[4].keyword,"sig");
  strcpy((*onfo_dict)[4].keyarg,"0.0");
  strcpy((*onfo_dict)[4].error_mes,"a number >= 0");
  /*------------------------------------------------------------------------*/ 
  /*  5) \eps{} */
  strcpy((*onfo_dict)[5].keyword,"eps");
  strcpy((*onfo_dict)[5].keyarg,"0.0");
  strcpy((*onfo_dict)[5].error_mes,"a number >= 0");
  /*------------------------------------------------------------------------*/ 
  /*  6) \scale{} */
  strcpy((*onfo_dict)[6].keyword,"scale");
  strcpy((*onfo_dict)[6].keyarg,"1.0");
  strcpy((*onfo_dict)[6].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 7) \label{} */
  strcpy((*onfo_dict)[7].keyword,"label");
  strcpy((*onfo_dict)[7].keyarg,"");
  strcpy((*onfo_dict)[7].error_mes,"");
  /*----------------------------------------------------------------------*/ 
}/*end routine*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_potfun_dict:intra potential fun key word dictionary                 */
/*==========================================================================*/
void set_potfun_dict(DICT_WORD *fun_dict[],int *num_fun_dict,int ifirst)
{/*begin routine*/
  int i;
  /*========================================================================*/
  /* 0) Malloc the dictionary */
  if(ifirst==1){
     *num_fun_dict=8;
     *fun_dict = (DICT_WORD *)cmalloc(*num_fun_dict*sizeof(DICT_WORD))-1;
   }/*endif*/
  /*========================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].iuset = 0;}
  for(i=1;i<=*num_fun_dict;i++){(*fun_dict)[i].key_type = 0;}
  /*========================================================================*/
  /* II) Fill the dictionary with words */
  /*-----------------------------------------------------------------------*/ 
  /* 1) ~inter_parm[] */
  strcpy((*fun_dict)[1].keyword,"inter_parm");
  strcpy((*fun_dict)[1].keyarg,"");
  strcpy((*fun_dict)[1].error_mes,"");
  (*fun_dict)[1].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 2) ~bond_parm[] */
  strcpy((*fun_dict)[2].keyword,"bond_parm");
  strcpy((*fun_dict)[2].keyarg,"");
  strcpy((*fun_dict)[2].error_mes,"");
  (*fun_dict)[2].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 3) ~bend_parm[] */
  strcpy((*fun_dict)[3].keyword,"bend_parm");
  strcpy((*fun_dict)[3].keyarg,"");
  strcpy((*fun_dict)[3].error_mes,"");
  (*fun_dict)[3].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 4) ~tors_parm[] */
  strcpy((*fun_dict)[4].keyword,"torsion_parm");
  strcpy((*fun_dict)[4].keyarg,"");
  strcpy((*fun_dict)[4].error_mes,"");
  (*fun_dict)[4].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 5) ~onfo_parm[] */
  strcpy((*fun_dict)[5].keyword,"onefour_parm");
  strcpy((*fun_dict)[5].keyarg,"");
  strcpy((*fun_dict)[5].error_mes,"");
  (*fun_dict)[5].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 6) ~pseudo_parm[] */
  strcpy((*fun_dict)[6].keyword,"pseudo_parm");
  strcpy((*fun_dict)[6].keyarg,"");
  strcpy((*fun_dict)[6].error_mes,"");
  (*fun_dict)[6].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 7) ~bend_bnd_parm[] */
  strcpy((*fun_dict)[7].keyword,"bend_bnd_parm");
  strcpy((*fun_dict)[7].keyarg,"");
  strcpy((*fun_dict)[7].error_mes,"");
  (*fun_dict)[7].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
  /* 8) ~surface_parm[] */
  strcpy((*fun_dict)[8].keyword,"surface_parm");
  strcpy((*fun_dict)[8].keyarg,"");
  strcpy((*fun_dict)[8].error_mes,"");
  (*fun_dict)[8].key_type = 2; /* specify more than once */
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*  set_potbend_bnd_dict:bend_bnd potential key word dictionary             */
/*==========================================================================*/

void set_potbend_bnd_dict(DICT_WORD *bend_bnd_dict[],int *num_bend_bnd_dict,
                          int ifirst)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
  int i;

/*========================================================================*/
/* 0) Malloc the dictionary */
  if(ifirst==1){
     *num_bend_bnd_dict=26;
     *bend_bnd_dict = (DICT_WORD *)cmalloc(*num_bend_bnd_dict*
                                            sizeof(DICT_WORD))-1;
   }/*endif*/
  /*=======================================================================*/
  /* I) Assign the users set flags 0 */
  for(i=1;i<=*num_bend_bnd_dict;i++){(*bend_bnd_dict)[i].iuset = 0;}
  for(i=1;i<=*num_bend_bnd_dict;i++){(*bend_bnd_dict)[i].key_type = 0;}
  /*=========================================================================*/
  /* II) Fill the dictionary with words */
  /*----------------------------------------------------------------------*/ 
  /*  1) \atom1{} */
  strcpy((*bend_bnd_dict)[1].keyword,"atom1");
  strcpy((*bend_bnd_dict)[1].keyarg,"");
  strcpy((*bend_bnd_dict)[1].error_mes,"");
  (*bend_bnd_dict)[1].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  2) \atom2{} */ 
  strcpy((*bend_bnd_dict)[2].keyword,"atom2");
  strcpy((*bend_bnd_dict)[2].keyarg,"");
  strcpy((*bend_bnd_dict)[2].error_mes,"");
  (*bend_bnd_dict)[2].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  3) \atom3{} */ 
  strcpy((*bend_bnd_dict)[3].keyword,"atom3");
  strcpy((*bend_bnd_dict)[3].keyarg,"");
  strcpy((*bend_bnd_dict)[3].error_mes,"");
  (*bend_bnd_dict)[3].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  4) \pot_type_bend{} */
  strcpy((*bend_bnd_dict)[4].keyword,"pot_type_bend");
  strcpy((*bend_bnd_dict)[4].keyarg,"");
  strcpy((*bend_bnd_dict)[4].error_mes,"harmonic,power-series,freq_series");
  (*bend_bnd_dict)[4].key_type = 1;  /* must spec*/
  /*-----------------------------------------------------------------------*/ 
  /*  5) \fk_bend{} */
  strcpy((*bend_bnd_dict)[5].keyword,"fk_bend");
  strcpy((*bend_bnd_dict)[5].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[5].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  6) \eq_bend{} */
  strcpy((*bend_bnd_dict)[6].keyword,"eq_bend");
  strcpy((*bend_bnd_dict)[6].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[6].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  7) \p0_bend{} */
  strcpy((*bend_bnd_dict)[7].keyword,"p0_bend");
  strcpy((*bend_bnd_dict)[7].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[7].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  8) \p1_bend{} */
  strcpy((*bend_bnd_dict)[8].keyword,"p1_bend");
  strcpy((*bend_bnd_dict)[8].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[8].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  9) \p2_bend{} */
  strcpy((*bend_bnd_dict)[9].keyword,"p2_bend");
  strcpy((*bend_bnd_dict)[9].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[9].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 10) \p3_bend{} */
  strcpy((*bend_bnd_dict)[10].keyword,"p3_bend");
  strcpy((*bend_bnd_dict)[10].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[10].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 11) \p4_bend{} */
  strcpy((*bend_bnd_dict)[11].keyword,"p4_bend");
  strcpy((*bend_bnd_dict)[11].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[11].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 12) \p5_bend{} */
  strcpy((*bend_bnd_dict)[12].keyword,"p5_bend");
  strcpy((*bend_bnd_dict)[12].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[12].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 13) \p6_bend{} */
  strcpy((*bend_bnd_dict)[13].keyword,"p6_bend");
  strcpy((*bend_bnd_dict)[13].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[13].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  14) \pot_type_bond{} */
  strcpy((*bend_bnd_dict)[14].keyword,"pot_type_bond");
  strcpy((*bend_bnd_dict)[14].keyarg,"Null");
  strcpy((*bend_bnd_dict)[14].error_mes,
                       "Null,harmonic,power-series,freq_series");
  /*-----------------------------------------------------------------------*/ 
  /*  15) \fk_bond{} */
  strcpy((*bend_bnd_dict)[15].keyword,"fk_bond");
  strcpy((*bend_bnd_dict)[15].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[15].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /*  16) \eq_bond{} */
  strcpy((*bend_bnd_dict)[16].keyword,"eq_bond");
  strcpy((*bend_bnd_dict)[16].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[16].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  17) \p0_bond{} */
  strcpy((*bend_bnd_dict)[17].keyword,"p0_bond");
  strcpy((*bend_bnd_dict)[17].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[17].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  18) \p1_bond{} */
  strcpy((*bend_bnd_dict)[18].keyword,"p1_bond");
  strcpy((*bend_bnd_dict)[18].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[18].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  19) \p2_bond{} */
  strcpy((*bend_bnd_dict)[19].keyword,"p2_bond");
  strcpy((*bend_bnd_dict)[19].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[19].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  20) \p3_bond{} */
  strcpy((*bend_bnd_dict)[20].keyword,"p3_bond");
  strcpy((*bend_bnd_dict)[20].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[20].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  21) \p4_bond{} */
  strcpy((*bend_bnd_dict)[21].keyword,"p4_bond");
  strcpy((*bend_bnd_dict)[21].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[21].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  22) \p5_bond{} */
  strcpy((*bend_bnd_dict)[22].keyword,"p5_bond");
  strcpy((*bend_bnd_dict)[22].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[22].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /*  23) \p6_bond{} */
  strcpy((*bend_bnd_dict)[23].keyword,"p6_bond");
  strcpy((*bend_bnd_dict)[23].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[23].error_mes,"a number");
  /*-----------------------------------------------------------------------*/ 
  /* 24) \alpha_bond} */
  strcpy((*bend_bnd_dict)[24].keyword,"alpha_bond");
  strcpy((*bend_bnd_dict)[24].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[24].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 25) \d0_bond{} */
  strcpy((*bend_bnd_dict)[25].keyword,"d0_bond");
  strcpy((*bend_bnd_dict)[25].keyarg,"0.0");
  strcpy((*bend_bnd_dict)[25].error_mes,"a number > 0");
  /*-----------------------------------------------------------------------*/ 
  /* 26) \label{} */
  strcpy((*bend_bnd_dict)[26].keyword,"label");
  strcpy((*bend_bnd_dict)[26].keyarg,"");
  strcpy((*bend_bnd_dict)[26].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
}/*end routine*/


