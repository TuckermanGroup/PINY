/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: bond_coef                                    */
/*                                                                          */
/* Subprogram to set up the power series coefficients for the bond          */
/* potential (6th order power series in (r-r0) is used)                     */
/*                                                                          */
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
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bond_coef(DICT_WORD *dict,char file_name[],char fun_key[],
                  DATA_BASE_ENTRIES *bond_base,CATM_LAB *cbond_base,int ibase)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */
       char *pot_typ;                             /* Intra potential type */
       int index;                                 /* Dictionary index     */
       double real_key_arg;                       /* Real key argument    */
       double fk_bond,eq_bond,deq_bond;           /* Harmonic bond        */
       double al_bond,d_bond;                     /* Morse bond           */
       double c0,c1,c2,c3,c4,c5,c6;               /* Power series         */
       int ifound=-1;                             /* Pot type match flag  */
       double a2,a3,a4,a5,a6;                     /* Powers of alpha      */

/*==========================================================================*/
/* 0) Fill atom types and label part of the data base     */
  strcpy(cbond_base[ibase].atm1,dict[1].keyarg);
  strcpy(cbond_base[ibase].atm2,dict[2].keyarg);
  strcpy(cbond_base[ibase].label,dict[15].keyarg);

/*==========================================================================*/
/* 0) Initial stuff                                                         */
   pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
   sscanf(dict[3].keyarg,"%s",pot_typ);

/*==========================================================================*/
/* 0) Constrained Bond                                                      */

      if(strcasecmp(pot_typ,"con") == 0) {
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       bond_base[ibase].eq = real_key_arg/BOHR;
       ifound = 0;
       cfree(pot_typ);
       return;
      }

/*==========================================================================*/
/* I) Harmonic Bond                                                         */

      if(strcasecmp(pot_typ,"harm") == 0) {
       sscanf(dict[4].keyarg,"%lg",&real_key_arg);
       fk_bond = real_key_arg;
       fk_bond *= (BOHR*BOHR/BOLTZ);
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       eq_bond = real_key_arg/BOHR;
       sscanf(dict[16].keyarg,"%lg",&real_key_arg);
       deq_bond = real_key_arg/BOHR;
       ifound = 0;
       bond_base[ibase].c_0 = 0.0;
       bond_base[ibase].c_1 = 0.0;
       bond_base[ibase].c_2 = 0.5*fk_bond;
       bond_base[ibase].c_3 = 0.0;
       bond_base[ibase].c_4 = 0.0;
       bond_base[ibase].c_5 = 0.0;
       bond_base[ibase].c_6 = 0.0;
       bond_base[ibase].eq  = eq_bond;
       bond_base[ibase].eq_res  = eq_bond + deq_bond;
     } /* endif harmonic */

/*==========================================================================*/
/* II) Morse                                                                */

      if(strcasecmp(pot_typ,"morse") == 0) {
       sscanf(dict[13].keyarg,"%lg",&real_key_arg);
       al_bond = real_key_arg;
       al_bond *= BOHR;
       sscanf(dict[14].keyarg,"%lg",&real_key_arg);
       d_bond = real_key_arg;
       d_bond /= BOLTZ;
       ifound = 0;
       a2 = al_bond*al_bond;
       a3 = al_bond*al_bond*al_bond;
       a4 = al_bond*al_bond*al_bond*al_bond;
       a5 = al_bond*al_bond*al_bond*al_bond*al_bond;
       a6 = al_bond*al_bond*al_bond*al_bond*al_bond*al_bond;
       bond_base[ibase].c_0 = 0.0;
       bond_base[ibase].c_1 = 0.0;
       bond_base[ibase].c_2 = a2*d_bond;
       bond_base[ibase].c_3 = -a3*d_bond;
       bond_base[ibase].c_4 = 7.0*a4*d_bond/12.0;
       bond_base[ibase].c_5 = -a5*d_bond/4.0;
       bond_base[ibase].c_6 = 62.0*a6*d_bond/720.0;
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       eq_bond = real_key_arg;
       eq_bond /= BOHR;
       bond_base[ibase].eq      = eq_bond;
       sscanf(dict[16].keyarg,"%lg",&real_key_arg);
       deq_bond = real_key_arg/BOHR;
       bond_base[ibase].eq_res  = eq_bond + deq_bond;
     } /* endif morse */

/*==========================================================================*/
/* III) Power Series                                                        */

      if(strcasecmp(pot_typ,"power") == 0) {
       sscanf(dict[6].keyarg,"%lg",&real_key_arg);
       c0 = real_key_arg;
       c0 /= BOLTZ;
       sscanf(dict[7].keyarg,"%lg",&real_key_arg);
       c1 = real_key_arg;
       c1 *= BOHR/BOLTZ;
       sscanf(dict[8].keyarg,"%lg",&real_key_arg);
       c2 = real_key_arg;
       c2 *= pow(BOHR,2.0)/BOLTZ;
       sscanf(dict[9].keyarg,"%lg",&real_key_arg);
       c3 = real_key_arg;
       c3 *= pow(BOHR,3.0)/BOLTZ;
       sscanf(dict[10].keyarg,"%lg",&real_key_arg);
       c4 = real_key_arg;
       c4 *= pow(BOHR,4.0)/BOLTZ;
       sscanf(dict[11].keyarg,"%lg",&real_key_arg);
       c5 = real_key_arg;
       c5 *= pow(BOHR,5.0)/BOLTZ;
       sscanf(dict[12].keyarg,"%lg",&real_key_arg);
       c6 = real_key_arg;
       c6 *= pow(BOHR,6.0)/BOLTZ;
       ifound = 0;
       bond_base[ibase].c_0 = c0;
       bond_base[ibase].c_1 = c1;
       bond_base[ibase].c_2 = c2;
       bond_base[ibase].c_3 = c3;
       bond_base[ibase].c_4 = c4;
       bond_base[ibase].c_5 = c5;
       bond_base[ibase].c_6 = c6;
       sscanf(dict[5].keyarg,"%lg",&real_key_arg);
       eq_bond = real_key_arg/BOHR;
       bond_base[ibase].eq  = eq_bond;
       sscanf(dict[16].keyarg,"%lg",&real_key_arg);
       deq_bond = real_key_arg/BOHR;
       bond_base[ibase].eq_res  = eq_bond + deq_bond;
     } /* endif power */

/*==========================================================================*/
/* IV) Null                                                                 */
      if(strcasecmp(pot_typ,"null") == 0) {
       ifound = 0;
       bond_base[ibase].c_0 = 0.0;
       bond_base[ibase].c_1 = 0.0;
       bond_base[ibase].c_2 = 0.0;
       bond_base[ibase].c_3 = 0.0;
       bond_base[ibase].c_4 = 0.0;
       bond_base[ibase].c_5 = 0.0;
       bond_base[ibase].c_6 = 0.0;
     } /* endif null */

/*==========================================================================*/
/* V) If pot type not found exit                                            */

       index=3;
       if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 

/*==========================================================================*/
/* VI) Calculate derivative coefficients                                     */

       bond_base[ibase].dc_0 = 0.0;
       bond_base[ibase].dc_1 =     bond_base[ibase].c_1;
       bond_base[ibase].dc_2 = 2.0*bond_base[ibase].c_2;
       bond_base[ibase].dc_3 = 3.0*bond_base[ibase].c_3;
       bond_base[ibase].dc_4 = 4.0*bond_base[ibase].c_4;
       bond_base[ibase].dc_5 = 5.0*bond_base[ibase].c_5;
       bond_base[ibase].dc_6 = 6.0*bond_base[ibase].c_6;

       cfree(pot_typ);
/*==========================================================================*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bend_coef(DICT_WORD *dict,char file_name[],char fun_key[],
                  DATA_BASE_ENTRIES *bend_base,CATM_LAB *cbend_base,int ibase)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  int ifound=-1;                            /* Pot type match flag     */
  int index;                                /* Dictionary index        */
  char *pot_typ;                            /* Potential type          */
  double real_key_arg;                      /* Value of key arg        */
  double fk_bend,theta_0;                   /* Harmonic bend           */
  double c0,c1,c2,c3,c4,c5,c6;              /* Power series            */
  double c2th0,c4th0,c6th0;                 /* Cos(n*theta_0)          */
  double s2th0,s4th0,s6th0;                 /* Sin(n*theta_0)          */
  static double cn1=49.0/72.0,cn2=3.0/4.0,cn3=3.0/40.0,cn4=1.0/180.0;
  static double cn5=3.0/2.0,cn6=3.0/5.0,cn7=1.0/10.0,cn8=4.0/15.0;
  static double cn9=8.0/45.0,cn10=3.0/10.0,cn11=1.0/30.0;
                                                 /* Useful constants        */

/*==========================================================================*/
/* 0) Fill atom types and label part of the data base     */
  strcpy(cbend_base[ibase].atm1,dict[1].keyarg);
  strcpy(cbend_base[ibase].atm2,dict[2].keyarg);
  strcpy(cbend_base[ibase].atm3,dict[3].keyarg);
  strcpy(cbend_base[ibase].label,dict[14].keyarg);

/*========================================================================*/
/* 0) Initial stuff                                                       */

  pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
  sscanf(dict[4].keyarg,"%s",pot_typ);

/*========================================================================*/
/* I) Harmonic Bend                                                       */
  
  if(strcasecmp(pot_typ,"harm") == 0) {
    sscanf(dict[5].keyarg,"%lg",&real_key_arg);
    fk_bend = real_key_arg;
    fk_bend /= BOLTZ;
    sscanf(dict[6].keyarg,"%lg",&real_key_arg);
    theta_0 = real_key_arg;
    theta_0 *= M_PI/180.0;
    ifound = 0;
    
    c2th0 = cos(2.0*theta_0);
    c4th0 = cos(4.0*theta_0);
    c6th0 = cos(6.0*theta_0);
    
    bend_base[ibase].c_0 =  0.5*fk_bend*(cn1 + cn2*c2th0 
                                       + cn3*c4th0 + cn4*c6th0);
    bend_base[ibase].c_1 = 0.0;
    bend_base[ibase].c_2 = 0.5*fk_bend*(-cn5*c2th0 - cn6*c4th0 - cn7*c6th0);
    bend_base[ibase].c_3 = 0.0;
    bend_base[ibase].c_4 = 0.5*fk_bend*(cn6*c4th0 + cn8*c6th0);
    bend_base[ibase].c_5 = 0.0;
    bend_base[ibase].c_6 = 0.5*fk_bend*(-cn9*c6th0);
    bend_base[ibase].eq   = theta_0;
    
    s2th0 = sin(2.0*theta_0);
    s4th0 = sin(4.0*theta_0);
    s6th0 = sin(6.0*theta_0);
    
    bend_base[ibase].s_0 = 0.0;
    bend_base[ibase].s_1 = 0.0;
    bend_base[ibase].s_2 = 0.5*fk_bend*(-cn5*s2th0 + cn10*s4th0 - cn11*s6th0);
    bend_base[ibase].s_3 = 0.0;
    bend_base[ibase].s_4 = 0.5*fk_bend*(-cn6*s4th0 + cn9*s6th0);
    bend_base[ibase].s_5 = 0.0;
    bend_base[ibase].s_6 = 0.5*fk_bend*(-cn9*s6th0);

  }  /* endif */
  
/*========================================================================*/
/* II) Power Series                                                      */
  
  if(strcasecmp(pot_typ,"power") == 0) {
    sscanf(dict[7].keyarg,"%lg",&real_key_arg);
    c0 = real_key_arg;
    c0 /= BOLTZ;
    sscanf(dict[8].keyarg,"%lg",&real_key_arg);
    c1 = real_key_arg;
    c1 /= BOLTZ;
    sscanf(dict[9].keyarg,"%lg",&real_key_arg);
    c2 = real_key_arg;
    c2 /= BOLTZ;
    sscanf(dict[10].keyarg,"%lg",&real_key_arg);
    c3 = real_key_arg;
    c3 /= BOLTZ;
    sscanf(dict[11].keyarg,"%lg",&real_key_arg);
    c4 = real_key_arg;
    c4 /= BOLTZ;
    sscanf(dict[12].keyarg,"%lg",&real_key_arg);
    c5 = real_key_arg;
    c5 /= BOLTZ;
    sscanf(dict[13].keyarg,"%lg",&real_key_arg);
    c6 = real_key_arg;
    c6 /= BOLTZ;

    ifound = 0;
    bend_base[ibase].c_0 = c0;
    bend_base[ibase].c_1 = c1;
    bend_base[ibase].c_2 = c2;
    bend_base[ibase].c_3 = c3;
    bend_base[ibase].c_4 = c4;
    bend_base[ibase].c_5 = c5;
    bend_base[ibase].c_6 = c6;

    bend_base[ibase].s_0 = 0.0;
    bend_base[ibase].s_1 = 0.0;
    bend_base[ibase].s_2 = 0.0;
    bend_base[ibase].s_3 = 0.0;
    bend_base[ibase].s_4 = 0.0;
    bend_base[ibase].s_5 = 0.0;
    bend_base[ibase].s_6 = 0.0;
    
  } /* endif power */

/*========================================================================*/
/* III) Null                                                             */

  if(strcasecmp(pot_typ,"null") == 0) {
    ifound = 0;
    bend_base[ibase].c_0 = 0.0;
    bend_base[ibase].c_1 = 0.0;
    bend_base[ibase].c_2 = 0.0;
    bend_base[ibase].c_3 = 0.0;
    bend_base[ibase].c_4 = 0.0;
    bend_base[ibase].c_5 = 0.0;
    bend_base[ibase].c_6 = 0.0;

    bend_base[ibase].s_0 = 0.0;
    bend_base[ibase].s_1 = 0.0;
    bend_base[ibase].s_2 = 0.0;
    bend_base[ibase].s_3 = 0.0;
    bend_base[ibase].s_4 = 0.0;
    bend_base[ibase].s_5 = 0.0;
    bend_base[ibase].s_6 = 0.0;
    
  } /* endif null */

/*========================================================================*/
/* IV) If pot type not found exit                                          */
  
  index = 4;
  if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 

/*========================================================================*/
/* V) Calculate derivative coefficients                                    */
  
  bend_base[ibase].dc_0 = 0.0;
  bend_base[ibase].dc_1 =     bend_base[ibase].c_1;
  bend_base[ibase].dc_2 = 2.0*bend_base[ibase].c_2;
  bend_base[ibase].dc_3 = 3.0*bend_base[ibase].c_3;
  bend_base[ibase].dc_4 = 4.0*bend_base[ibase].c_4;
  bend_base[ibase].dc_5 = 5.0*bend_base[ibase].c_5;
  bend_base[ibase].dc_6 = 6.0*bend_base[ibase].c_6;
  
  bend_base[ibase].ds_0 = 0.0;
  bend_base[ibase].ds_1 =     bend_base[ibase].s_1;
  bend_base[ibase].ds_2 =     bend_base[ibase].s_2;
  bend_base[ibase].ds_3 = 2.0*bend_base[ibase].s_3;
  bend_base[ibase].ds_4 = 3.0*bend_base[ibase].s_4;
  bend_base[ibase].ds_5 = 4.0*bend_base[ibase].s_5;
  bend_base[ibase].ds_6 = 5.0*bend_base[ibase].s_6;

  cfree(pot_typ);
/*=======================================================================*/
} /* end routine */
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void tors_coef(DICT_WORD *dict,char file_name[],char fun_key[],
                  DATA_BASE_ENTRIES *tors_base,CATM_LAB *ctors_base,int ibase)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  int ifound=-1;                            /* Pot type match flag     */
  int index;                                /* Dictionary index        */
  char *pot_typ;                            /* Potential type          */
  double real_key_arg;                      /* Value of key arg        */
  double fk_tors,theta_0;                   /* Harmonic torsion        */
  double c0,c1,c2,c3,c4,c5,c6;              /* Power series            */
  double d1,d2,d3,d4,a1,a2,a3,a4;int nfreq; /* Frequency Series     */
  double c[5],cos_d[5]; 
  double cp0[5],cp1[5],cp2[5],cp3[5],cp4[5],cp5[5],cp6[5]; 
  int n_tors[5];
  int ifreq;
  double fval,i_ptr;
  double c2th0,c4th0,c6th0;                 /* Cos(n*theta_0)          */
  double s2th0,s4th0,s6th0;                 /* Sin(n*theta_0)          */
  /* Useful constants        */
  static double cn1=49.0/72.0,cn2=3.0/4.0,cn3=3.0/40.0,cn4=1.0/180.0;
  static double cn5=3.0/2.0,cn6=3.0/5.0,cn7=1.0/10.0,cn8=4.0/15.0;
  static double cn9=8.0/45.0,cn10=3.0/10.0,cn11=1.0/30.0;
  static double small=1.0e-8;

/*==========================================================================*/
/* 0) Fill atom types and label part of the data base     */
  strcpy(ctors_base[ibase].atm1,dict[1].keyarg);
  strcpy(ctors_base[ibase].atm2,dict[2].keyarg);
  strcpy(ctors_base[ibase].atm3,dict[3].keyarg);
  strcpy(ctors_base[ibase].atm4,dict[4].keyarg);
  strcpy(ctors_base[ibase].label,dict[25].keyarg);

/*=======================================================================*/
/* 0) Initial stuff                                                      */

  pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
  sscanf(dict[5].keyarg,"%s",pot_typ);

/*=======================================================================*/
/* I) Harmonic tors                                                      */
  
  if(strcasecmp(pot_typ,"harm") == 0) {
    sscanf(dict[6].keyarg,"%lg",&real_key_arg);
    fk_tors = real_key_arg;
    fk_tors /= BOLTZ;
    sscanf(dict[7].keyarg,"%lg",&real_key_arg);
    theta_0 = real_key_arg;
    theta_0 *= M_PI/180.0;
    ifound = 0;

    c2th0 = cos(2.0*theta_0);
    c4th0 = cos(4.0*theta_0);
    c6th0 = cos(6.0*theta_0);

    tors_base[ibase].c_0 = 0.5*fk_tors*(cn1 + cn2*c2th0 + cn3*c4th0 
                                            + cn4*c6th0);
    tors_base[ibase].c_1 = 0.0;
    tors_base[ibase].c_2 = 0.5*fk_tors*(-cn5*c2th0 - cn6*c4th0 - cn7*c6th0);
    tors_base[ibase].c_3 = 0.0;
    tors_base[ibase].c_4 = 0.5*fk_tors*(cn6*c4th0 + cn8*c6th0);
    tors_base[ibase].c_5 = 0.0;
    tors_base[ibase].c_6 = 0.5*fk_tors*(-cn9*c6th0);
    tors_base[ibase].eq  = theta_0;
    
    s2th0 = sin(2.0*theta_0);
    s4th0 = sin(4.0*theta_0);
    s6th0 = sin(6.0*theta_0);
    
    tors_base[ibase].s_0 = 0.0;
    tors_base[ibase].s_1 = 0.0;
    tors_base[ibase].s_2 = 0.5*fk_tors*(-cn5*s2th0 + cn10*s4th0 - cn11*s6th0);
    tors_base[ibase].s_3 = 0.0;
    tors_base[ibase].s_4 = 0.5*fk_tors*(-cn6*s4th0 + cn9*s6th0);
    tors_base[ibase].s_5 = 0.0;
    tors_base[ibase].s_6 = 0.5*fk_tors*(-cn9*s6th0);

  }  /* endif */

/*========================================================================*/
/* II) Power Series                                                     */
  
  if(strcasecmp(pot_typ,"power") == 0) {
    sscanf(dict[8].keyarg,"%lg",&real_key_arg);
    c0 = real_key_arg;
    c0 /= BOLTZ;
    sscanf(dict[9].keyarg,"%lg",&real_key_arg);
    c1 = real_key_arg;
    c1 /= BOLTZ;
    sscanf(dict[10].keyarg,"%lg",&real_key_arg);
    c2 = real_key_arg;
    c2 /= BOLTZ;
    sscanf(dict[11].keyarg,"%lg",&real_key_arg);
    c3 = real_key_arg;
    c3 /= BOLTZ;
    sscanf(dict[12].keyarg,"%lg",&real_key_arg);
    c4 = real_key_arg;
    c4 /= BOLTZ;
    sscanf(dict[13].keyarg,"%lg",&real_key_arg);
    c5 = real_key_arg;
    c5 /= BOLTZ;
    sscanf(dict[14].keyarg,"%lg",&real_key_arg);
    c6 = real_key_arg;
    c6 /= BOLTZ;
    ifound = 0;
    tors_base[ibase].c_0 = c0;
    tors_base[ibase].c_1 = c1;
    tors_base[ibase].c_2 = c2;
    tors_base[ibase].c_3 = c3;
    tors_base[ibase].c_4 = c4;
    tors_base[ibase].c_5 = c5;
    tors_base[ibase].c_6 = c6;
    
    tors_base[ibase].s_0 = 0.0;
    tors_base[ibase].s_1 = 0.0;
    tors_base[ibase].s_2 = 0.0;
    tors_base[ibase].s_3 = 0.0;
    tors_base[ibase].s_4 = 0.0;
    tors_base[ibase].s_5 = 0.0;
    tors_base[ibase].s_6 = 0.0;
    
  } /* endif power */
  
/*========================================================================*/
/* III) Frequency Series                                                */

  if(strcasecmp(pot_typ,"freq-series") == 0) {
    sscanf(dict[15].keyarg,"%d",&nfreq);
    index =15;
    if((nfreq <= 0 || nfreq > 4) ){
      keyarg_barf(dict,file_name,fun_key,index);
    }
/* 1) 1st component */
    sscanf(dict[16].keyarg,"%lg",&real_key_arg);
    a1 = real_key_arg;
    index = 16;
    fval = modf(a1,&i_ptr);
    if(fval > small) keyarg_barf(dict,file_name,fun_key,index);

    sscanf(dict[17].keyarg,"%lg",&real_key_arg);
    c[1] = real_key_arg/BOLTZ;

    sscanf(dict[18].keyarg,"%lg",&real_key_arg);
    d1 = real_key_arg;
    index = 18;
    if(!(d1 == 0.0 || d1 == 180.0) ){
      keyarg_barf(dict,file_name,fun_key,index);
    }
/* 2) 2nd component */
    sscanf(dict[19].keyarg,"%lg",&real_key_arg);
    a2 = real_key_arg;
    index = 19;
    fval = modf(a1,&i_ptr);
    if(fval > small) keyarg_barf(dict,file_name,fun_key,index);

    sscanf(dict[20].keyarg,"%lg",&real_key_arg);
    c[2] = real_key_arg/BOLTZ;

    sscanf(dict[21].keyarg,"%lg",&real_key_arg);
    d2 = real_key_arg;
    index = 21;
    if(!(d2 == 0.0 || d2 == 180.0) ){
      keyarg_barf(dict,file_name,fun_key,index);
    }
/* 3) 3nd component */
    sscanf(dict[22].keyarg,"%lg",&real_key_arg);
    a3 = real_key_arg;
    index = 22;
    fval = modf(a1,&i_ptr);
    if(fval > small) keyarg_barf(dict,file_name,fun_key,index);

    sscanf(dict[23].keyarg,"%lg",&real_key_arg);
    c[3] = real_key_arg/BOLTZ;

    sscanf(dict[24].keyarg,"%lg",&real_key_arg);
    d3 = real_key_arg;
    index = 24;
    if(!(d3 == 0.0 || d3 == 180.0) ){
      keyarg_barf(dict,file_name,fun_key,index);
    }
/* 4) 4th component */
    sscanf(dict[26].keyarg,"%lg",&real_key_arg);
    a4 = real_key_arg;
    index = 26;
    fval = modf(a1,&i_ptr);
    if(fval > small) keyarg_barf(dict,file_name,fun_key,index);

    sscanf(dict[27].keyarg,"%lg",&real_key_arg);
    c[4] = real_key_arg/BOLTZ;

    sscanf(dict[28].keyarg,"%lg",&real_key_arg);
    d4 = real_key_arg;
    index = 24;
    if(!(d4 == 0.0 || d4 == 180.0) ){
      keyarg_barf(dict,file_name,fun_key,index);
    }


    n_tors[1] = (int) a1;
    n_tors[2] = (int) a2;
    n_tors[3] = (int) a3;
    n_tors[4] = (int) a4;
    d1 *= M_PI/180.0;
    d2 *= M_PI/180.0;
    d3 *= M_PI/180.0;
    d4 *= M_PI/180.0;
    cos_d[1] = cos(d1);
    cos_d[2] = cos(d2);
    cos_d[3] = cos(d3);
    cos_d[4] = cos(d4);

    ifound = 0;
    for(ifreq=1;ifreq<=nfreq;ifreq++){
      switch(n_tors[ifreq]) {
      case 0:
        cp0[ifreq] =  1.0;
        cp1[ifreq] =  0.0;
        cp2[ifreq] =  0.0;
        cp3[ifreq] =  0.0;
        cp4[ifreq] =  0.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  0.0;
        break;
      case 1:
        cp0[ifreq] =  0.0;
        cp1[ifreq] =  1.0;
        cp2[ifreq] =  0.0;
        cp3[ifreq] =  0.0;
        cp4[ifreq] =  0.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  0.0;
        break;
      case 2:
        cp0[ifreq] = -1.0;
        cp1[ifreq] =  0.0;
        cp2[ifreq] =  2.0;
        cp3[ifreq] =  0.0;
        cp4[ifreq] =  0.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  0.0;
        break;
      case 3:
        cp0[ifreq] =  0.0;
        cp1[ifreq] = -3.0;
        cp2[ifreq] =  0.0;
        cp3[ifreq] =  4.0;
        cp4[ifreq] =  0.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  0.0;
        break;
      case 4:
        cp0[ifreq] =  1.0;
        cp1[ifreq] =  0.0;
        cp2[ifreq] = -8.0;
        cp3[ifreq] =  0.0;
        cp4[ifreq] =  8.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  0.0;
        break;
      case 5:
        cp0[ifreq] =  0.0;
        cp1[ifreq] =  5.0;
        cp2[ifreq] =  0.0;
        cp3[ifreq] = -20.0;
        cp4[ifreq] =  0.0;
        cp5[ifreq] =  16.0;
        cp6[ifreq] =  0.0;
        break;
      case 6:
        cp0[ifreq] = -1.0;
        cp1[ifreq] =  0.0;
        cp2[ifreq] =  18.0;
        cp3[ifreq] =  0.0;
        cp4[ifreq] = -48.0;
        cp5[ifreq] =  0.0;
        cp6[ifreq] =  32.0;
        break;
      default:
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Frequency series with n > 6, n < 0\n");
        printf("or n not an integer\n");
        printf("not allowed\n");
        printf("n = %d\n",n_tors[ifreq]);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1); 
        break;
      } /* end switch(n_tors[ifreq]) */  
 
      cp0[ifreq] *= cos_d[ifreq];
      cp1[ifreq] *= cos_d[ifreq];
      cp2[ifreq] *= cos_d[ifreq];
      cp3[ifreq] *= cos_d[ifreq];
      cp4[ifreq] *= cos_d[ifreq];
      cp5[ifreq] *= cos_d[ifreq];
      cp6[ifreq] *= cos_d[ifreq];
    }/*endfor*/
    
    tors_base[ibase].c_0=0.0;
    tors_base[ibase].c_1=0.0;
    tors_base[ibase].c_2=0.0;
    tors_base[ibase].c_3=0.0;
    tors_base[ibase].c_4=0.0;
    tors_base[ibase].c_5=0.0;
    tors_base[ibase].c_6=0.0;
    for(ifreq=1;ifreq<=nfreq;ifreq++){
      tors_base[ibase].c_0 +=  c[ifreq];
    }
    for(ifreq=1;ifreq<=nfreq;ifreq++){
      tors_base[ibase].c_0 -=  c[ifreq]*cp0[ifreq];
      tors_base[ibase].c_1 -=  c[ifreq]*cp1[ifreq];
      tors_base[ibase].c_2 -=  c[ifreq]*cp2[ifreq];
      tors_base[ibase].c_3 -=  c[ifreq]*cp3[ifreq];
      tors_base[ibase].c_4 -=  c[ifreq]*cp4[ifreq];
      tors_base[ibase].c_5 -=  c[ifreq]*cp5[ifreq];
      tors_base[ibase].c_6 -=  c[ifreq]*cp6[ifreq];
    }


    
    tors_base[ibase].s_0 = 0.0;
    tors_base[ibase].s_1 = 0.0;
    tors_base[ibase].s_2 = 0.0;
    tors_base[ibase].s_3 = 0.0;
    tors_base[ibase].s_4 = 0.0;
    tors_base[ibase].s_5 = 0.0;
    tors_base[ibase].s_6 = 0.0;
    
  }/* endif frequency series */

/*========================================================================*/
/* III) Null                                                             */
  if(strcasecmp(pot_typ,"null") == 0) {
    ifound = 0;
    tors_base[ibase].c_0 = 0.0;
    tors_base[ibase].c_1 = 0.0;
    tors_base[ibase].c_2 = 0.0;
    tors_base[ibase].c_3 = 0.0;
    tors_base[ibase].c_4 = 0.0;
    tors_base[ibase].c_5 = 0.0;
    tors_base[ibase].c_6 = 0.0;

    tors_base[ibase].s_0 = 0.0;
    tors_base[ibase].s_1 = 0.0;
    tors_base[ibase].s_2 = 0.0;
    tors_base[ibase].s_3 = 0.0;
    tors_base[ibase].s_4 = 0.0;
    tors_base[ibase].s_5 = 0.0;
    tors_base[ibase].s_6 = 0.0;
    
  } /* endif null */

/*========================================================================*/
/* IV ) If pot type not found exit                                       */
  
  index = 5;
  if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 
  
/*========================================================================*/
/* V ) Calculate derivative coefficients                                 */
  
  tors_base[ibase].dc_0 = 0.0;
  tors_base[ibase].dc_1 =     tors_base[ibase].c_1;
  tors_base[ibase].dc_2 = 2.0*tors_base[ibase].c_2;
  tors_base[ibase].dc_3 = 3.0*tors_base[ibase].c_3;
  tors_base[ibase].dc_4 = 4.0*tors_base[ibase].c_4;
  tors_base[ibase].dc_5 = 5.0*tors_base[ibase].c_5;
  tors_base[ibase].dc_6 = 6.0*tors_base[ibase].c_6;
  
  tors_base[ibase].ds_0 = 0.0;
  tors_base[ibase].ds_1 =     tors_base[ibase].s_1;
  tors_base[ibase].ds_2 =     tors_base[ibase].s_2;
  tors_base[ibase].ds_3 = 2.0*tors_base[ibase].s_3;
  tors_base[ibase].ds_4 = 3.0*tors_base[ibase].s_4;
  tors_base[ibase].ds_5 = 4.0*tors_base[ibase].s_5;
  tors_base[ibase].ds_6 = 5.0*tors_base[ibase].s_6;
  
  cfree(pot_typ);

/*=======================================================================*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void onfo_coef(DICT_WORD *dict,char file_name[],char fun_key[],
                 DATA_BASE_ENTRIES *onfo_base,CATM_LAB *confo_base,int ibase)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*            Local variable declarations:                               */

  char *pot_typ;                             /* Intra potential type */
  int index;                                 /* Dictionary index     */
  double real_key_arg;                       /* Real key argument    */
  double sig,eps,sc;                         /* Lennard-Jones        */
  double feps,s6;
  int ifound=-1;                             /* Pot type match flag  */

/*==========================================================================*/
/* 0) Fill atom types and label part of the data base     */
  strcpy(confo_base[ibase].atm1,dict[1].keyarg);
  strcpy(confo_base[ibase].atm2,dict[2].keyarg);
  strcpy(confo_base[ibase].label,dict[7].keyarg);

/*=======================================================================*/
/* 0) Initial stuff                                                      */

  pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
  sscanf(dict[3].keyarg,"%s",pot_typ);

/*=======================================================================*/
/* I) Lennard-Jones                                                      */
  
  if(strcasecmp(pot_typ,"lennard-jones") == 0) {
    sscanf(dict[4].keyarg,"%lg",&real_key_arg);
    sig = real_key_arg;
    sig /= BOHR;
    sscanf(dict[5].keyarg,"%lg",&real_key_arg);
    eps = real_key_arg;
    eps /= BOLTZ;
    sscanf(dict[6].keyarg,"%lg",&real_key_arg);
    sc = real_key_arg;
    ifound = 0;
    feps = 4.0*eps;
    s6 = pow(sig,6.0);
    onfo_base[ibase].feps = feps;
    onfo_base[ibase].s6   = s6;
    onfo_base[ibase].sc   = sc;
  } /* endif lennard-jones */
  
/*=======================================================================*/
/* II) Null                                                             */

  if(strcasecmp(pot_typ,"null") == 0) {
    ifound = 0;
    onfo_base[ibase].feps = 0.0;
    onfo_base[ibase].s6   = 0.0;
    onfo_base[ibase].sc   = 0.0;
  } /* endif null */

/*=======================================================================*/
/* V) If pot type not found exit                                         */
  
  index=3;
  if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 
  
  cfree(pot_typ);

/*========================================================================*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void bend_bnd_coef(DICT_WORD *dict,char file_name[],char fun_key[],
                    DATA_BASE_ENTRIES *bend_bnd_base,CATM_LAB *cbend_bnd_base,
                    int ibase)

/*==========================================================================*/
/*               Begin subprogram:                                          */
      {/*begin routine*/
/*==========================================================================*/
/*               Local variable declarations:                               */
       char *pot_typ;                             /* Intra potential type */
       double real_key_arg;                       /* Real key argument    */
       double fk_bond,eq_bond;                    /* Harmonic bond        */
       double al_bond,d_bond;                     /* Morse bond           */
       double c0,c1,c2,c3,c4,c5,c6;               /* Power series         */
       int ifound=-1;                             /* Pot type match flag  */
       double a2,a3,a4,a5,a6;                     /* Powers of alpha      */
       double fk_bend,theta_0;
  double c2th0,c4th0,c6th0;                 /* Cos(n*theta_0)          */
  double s2th0,s4th0,s6th0;                 /* Sin(n*theta_0)          */
  static double cn1=49.0/72.0,cn2=3.0/4.0,cn3=3.0/40.0,cn4=1.0/180.0;
  static double cn5=3.0/2.0,cn6=3.0/5.0,cn7=1.0/10.0,cn8=4.0/15.0;
  static double cn9=8.0/45.0,cn10=3.0/10.0,cn11=1.0/30.0;
       int index;
/*==========================================================================*/
/*==========================================================================*/
/*  GET THE BOND PORTION OF THE BEND-BOND                                   */ 

/*==========================================================================*/
/* 0) Fill atom types and label part of the data base     */
  strcpy(cbend_bnd_base[ibase].atm1,dict[1].keyarg);
  strcpy(cbend_bnd_base[ibase].atm2,dict[2].keyarg);
  strcpy(cbend_bnd_base[ibase].atm3,dict[3].keyarg);
  strcpy(cbend_bnd_base[ibase].label,dict[26].keyarg);

/*==========================================================================*/
/* 0) Fetch Bond type                                                       */

  bend_bnd_base[ibase].ipure = 0;
   pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
   sscanf(dict[14].keyarg,"%s",pot_typ);

/*==========================================================================*/
/* 0) Eq Bond length                                                        */

   sscanf(dict[16].keyarg,"%lg",&real_key_arg);
   eq_bond = real_key_arg;
   eq_bond /= BOHR;
   bend_bnd_base[ibase].eqb  = eq_bond;
   index = 16;
   if(eq_bond < 0) keyarg_barf(dict,file_name,fun_key,index); 
      
/*==========================================================================*/
/* I) Harmonic Bond                                                         */

      if(strcasecmp(pot_typ,"harm") == 0) {
       sscanf(dict[15].keyarg,"%lg",&real_key_arg);
       fk_bond = real_key_arg;
       fk_bond *= (BOHR*BOHR/BOLTZ);
       ifound = 0;
       bend_bnd_base[ibase].cb_0 = 0.0;
       bend_bnd_base[ibase].cb_1 = 0.0;
       bend_bnd_base[ibase].cb_2 = 0.5*fk_bond;
       bend_bnd_base[ibase].cb_3 = 0.0;
       bend_bnd_base[ibase].cb_4 = 0.0;
       bend_bnd_base[ibase].cb_5 = 0.0;
       bend_bnd_base[ibase].cb_6 = 0.0;
     } /* endif harmonic */

/*==========================================================================*/
/* II) Morse                                                                */

      if(strcasecmp(pot_typ,"morse") == 0) {
       sscanf(dict[24].keyarg,"%lg",&real_key_arg);
       al_bond = real_key_arg;
       al_bond *= BOHR;
       sscanf(dict[25].keyarg,"%lg",&real_key_arg);
       d_bond = real_key_arg;
       d_bond /= BOLTZ;
       ifound = 0;
       a2 = al_bond*al_bond;
       a3 = al_bond*al_bond*al_bond;
       a4 = al_bond*al_bond*al_bond*al_bond;
       a5 = al_bond*al_bond*al_bond*al_bond*al_bond;
       a6 = al_bond*al_bond*al_bond*al_bond*al_bond*al_bond;
       bend_bnd_base[ibase].cb_0 = 0.0;
       bend_bnd_base[ibase].cb_1 = 0.0;
       bend_bnd_base[ibase].cb_2 = a2*d_bond;
       bend_bnd_base[ibase].cb_3 = -a3*d_bond;
       bend_bnd_base[ibase].cb_4 = 7.0*a4*d_bond/12.0;
       bend_bnd_base[ibase].cb_5 = -a5*d_bond/4.0;
       bend_bnd_base[ibase].cb_6 = 62.0*a6*d_bond/720.0;
     } /* endif morse */

/*==========================================================================*/
/* III) Power Series                                                        */

      if(strcasecmp(pot_typ,"power") == 0) {
       sscanf(dict[17].keyarg,"%lg",&real_key_arg);
       c0 = real_key_arg;
       c0 /= BOLTZ;
       sscanf(dict[18].keyarg,"%lg",&real_key_arg);
       c1 = real_key_arg;
       c1 *= BOHR/BOLTZ;
       sscanf(dict[19].keyarg,"%lg",&real_key_arg);
       c2 = real_key_arg;
       c2 *= pow(BOHR,2.0)/BOLTZ;
       sscanf(dict[20].keyarg,"%lg",&real_key_arg);
       c3 = real_key_arg;
       c3 *= pow(BOHR,3.0)/BOLTZ;
       sscanf(dict[21].keyarg,"%lg",&real_key_arg);
       c4 = real_key_arg;
       c4 *= pow(BOHR,4.0)/BOLTZ;
       sscanf(dict[22].keyarg,"%lg",&real_key_arg);
       c5 = real_key_arg;
       c5 *= pow(BOHR,5.0)/BOLTZ;
       sscanf(dict[23].keyarg,"%lg",&real_key_arg);
       c6 = real_key_arg;
       c6 *= pow(BOHR,6.0)/BOLTZ;
       ifound = 0;
       bend_bnd_base[ibase].cb_0 = c0;
       bend_bnd_base[ibase].cb_1 = c1;
       bend_bnd_base[ibase].cb_2 = c2;
       bend_bnd_base[ibase].cb_3 = c3;
       bend_bnd_base[ibase].cb_4 = c4;
       bend_bnd_base[ibase].cb_5 = c5;
       bend_bnd_base[ibase].cb_6 = c6;
     } /* endif power */

/*==========================================================================*/
/* IV) Null                                                                 */
      if(strcasecmp(pot_typ,"null") == 0) {
       ifound = 0;
       bend_bnd_base[ibase].cb_0 = 0.0;
       bend_bnd_base[ibase].cb_1 = 0.0;
       bend_bnd_base[ibase].cb_2 = 0.0;
       bend_bnd_base[ibase].cb_3 = 0.0;
       bend_bnd_base[ibase].cb_4 = 0.0;
       bend_bnd_base[ibase].cb_5 = 0.0;
       bend_bnd_base[ibase].cb_6 = 0.0;
       bend_bnd_base[ibase].ipure = 1;
     } /* endif null */

/*==========================================================================*/
/* V) If bond pot type not found exit                                       */

       index = 14;
       if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 

/*==========================================================================*/
/* VI) Calculate bond derivative coefficients                               */

       bend_bnd_base[ibase].dcb_0 = 0.0;
       bend_bnd_base[ibase].dcb_1 =     bend_bnd_base[ibase].cb_1;
       bend_bnd_base[ibase].dcb_2 = 2.0*bend_bnd_base[ibase].cb_2;
       bend_bnd_base[ibase].dcb_3 = 3.0*bend_bnd_base[ibase].cb_3;
       bend_bnd_base[ibase].dcb_4 = 4.0*bend_bnd_base[ibase].cb_4;
       bend_bnd_base[ibase].dcb_5 = 5.0*bend_bnd_base[ibase].cb_5;
       bend_bnd_base[ibase].dcb_6 = 6.0*bend_bnd_base[ibase].cb_6;

/*==========================================================================*/
/*==========================================================================*/
/*  GET THE BEND PORTION OF THE BEND-BOND                                   */ 


/*========================================================================*/
/* 0) Get potential type                                                  */

  pot_typ = (char *) cmalloc(MAXWORD*sizeof(char));
  sscanf(dict[4].keyarg,"%s",pot_typ);

/*========================================================================*/
/* I) Harmonic Bend                                                       */
  
  if(strcasecmp(pot_typ,"harm") == 0) {
    sscanf(dict[5].keyarg,"%lg",&real_key_arg);
    fk_bend = real_key_arg;
    fk_bend /= BOLTZ;
    sscanf(dict[6].keyarg,"%lg",&real_key_arg);
    theta_0 = real_key_arg;
    theta_0 *= M_PI/180.0;
    ifound = 0;
    
    c2th0 = cos(2.0*theta_0);
    c4th0 = cos(4.0*theta_0);
    c6th0 = cos(6.0*theta_0);
    
    bend_bnd_base[ibase].c_0 =  0.5*fk_bend*(cn1 + cn2*c2th0 + cn3*c4th0 
                                           + cn4*c6th0);
    bend_bnd_base[ibase].c_1 = 0.0;
    bend_bnd_base[ibase].c_2 = 0.5*fk_bend*(-cn5*c2th0 - cn6*c4th0 
                                           - cn7*c6th0);
    bend_bnd_base[ibase].c_3 = 0.0;
    bend_bnd_base[ibase].c_4 = 0.5*fk_bend*(cn6*c4th0 + cn8*c6th0);
    bend_bnd_base[ibase].c_5 = 0.0;
    bend_bnd_base[ibase].c_6 = 0.5*fk_bend*(-cn9*c6th0);
    bend_bnd_base[ibase].eq  = theta_0;
    
    s2th0 = sin(2.0*theta_0);
    s4th0 = sin(4.0*theta_0);
    s6th0 = sin(6.0*theta_0);
    
    bend_bnd_base[ibase].s_0 = 0.0;
    bend_bnd_base[ibase].s_1 = 0.0;
    bend_bnd_base[ibase].s_2 = 0.5*fk_bend*(-cn5*s2th0 + cn10*s4th0 
                                           - cn11*s6th0);
    bend_bnd_base[ibase].s_3 = 0.0;
    bend_bnd_base[ibase].s_4 = 0.5*fk_bend*(-cn6*s4th0 + cn9*s6th0);
    bend_bnd_base[ibase].s_5 = 0.0;
    bend_bnd_base[ibase].s_6 = 0.5*fk_bend*(-cn9*s6th0);

  }  /* endif */
  
/*========================================================================*/
/* II) Power Series                                                      */
  
  if(strcasecmp(pot_typ,"power") == 0) {
    sscanf(dict[7].keyarg,"%lg",&real_key_arg);
    c0 = real_key_arg;
    c0 /= BOLTZ;
    sscanf(dict[8].keyarg,"%lg",&real_key_arg);
    c1 = real_key_arg;
    c1 /= BOLTZ;
    sscanf(dict[9].keyarg,"%lg",&real_key_arg);
    c2 = real_key_arg;
    c2 /= BOLTZ;
    sscanf(dict[10].keyarg,"%lg",&real_key_arg);
    c3 = real_key_arg;
    c3 /= BOLTZ;
    sscanf(dict[11].keyarg,"%lg",&real_key_arg);
    c4 = real_key_arg;
    c4 /= BOLTZ;
    sscanf(dict[12].keyarg,"%lg",&real_key_arg);
    c5 = real_key_arg;
    c5 /= BOLTZ;
    sscanf(dict[13].keyarg,"%lg",&real_key_arg);
    c6 = real_key_arg;
    c6 /= BOLTZ;

    ifound = 0;
    bend_bnd_base[ibase].c_0 = c0;
    bend_bnd_base[ibase].c_1 = c1;
    bend_bnd_base[ibase].c_2 = c2;
    bend_bnd_base[ibase].c_3 = c3;
    bend_bnd_base[ibase].c_4 = c4;
    bend_bnd_base[ibase].c_5 = c5;
    bend_bnd_base[ibase].c_6 = c6;

    bend_bnd_base[ibase].s_0 = 0.0;
    bend_bnd_base[ibase].s_1 = 0.0;
    bend_bnd_base[ibase].s_2 = 0.0;
    bend_bnd_base[ibase].s_3 = 0.0;
    bend_bnd_base[ibase].s_4 = 0.0;
    bend_bnd_base[ibase].s_5 = 0.0;
    bend_bnd_base[ibase].s_6 = 0.0;
    
  } /* endif power */

/*========================================================================*/
/* III) Null                                                             */

  if(strcasecmp(pot_typ,"null") == 0) {
    ifound = 0;
    bend_bnd_base[ibase].c_0 = 0.0;
    bend_bnd_base[ibase].c_1 = 0.0;
    bend_bnd_base[ibase].c_2 = 0.0;
    bend_bnd_base[ibase].c_3 = 0.0;
    bend_bnd_base[ibase].c_4 = 0.0;
    bend_bnd_base[ibase].c_5 = 0.0;
    bend_bnd_base[ibase].c_6 = 0.0;

    bend_bnd_base[ibase].s_0 = 0.0;
    bend_bnd_base[ibase].s_1 = 0.0;
    bend_bnd_base[ibase].s_2 = 0.0;
    bend_bnd_base[ibase].s_3 = 0.0;
    bend_bnd_base[ibase].s_4 = 0.0;
    bend_bnd_base[ibase].s_5 = 0.0;
    bend_bnd_base[ibase].s_6 = 0.0;
    
  } /* endif null */

/*========================================================================*/
/* IV) If pot type not found exit                                          */
  
  index = 4;
  if(ifound < 0) keyarg_barf(dict,file_name,fun_key,index); 

/*========================================================================*/
/* V) Calculate derivative coefficients                                    */
  
  bend_bnd_base[ibase].dc_0 = 0.0;
  bend_bnd_base[ibase].dc_1 =     bend_bnd_base[ibase].c_1;
  bend_bnd_base[ibase].dc_2 = 2.0*bend_bnd_base[ibase].c_2;
  bend_bnd_base[ibase].dc_3 = 3.0*bend_bnd_base[ibase].c_3;
  bend_bnd_base[ibase].dc_4 = 4.0*bend_bnd_base[ibase].c_4;
  bend_bnd_base[ibase].dc_5 = 5.0*bend_bnd_base[ibase].c_5;
  bend_bnd_base[ibase].dc_6 = 6.0*bend_bnd_base[ibase].c_6;
  
  bend_bnd_base[ibase].ds_0 = 0.0;
  bend_bnd_base[ibase].ds_1 =     bend_bnd_base[ibase].s_1;
  bend_bnd_base[ibase].ds_2 =     bend_bnd_base[ibase].s_2;
  bend_bnd_base[ibase].ds_3 = 2.0*bend_bnd_base[ibase].s_3;
  bend_bnd_base[ibase].ds_4 = 3.0*bend_bnd_base[ibase].s_4;
  bend_bnd_base[ibase].ds_5 = 4.0*bend_bnd_base[ibase].s_5;
  bend_bnd_base[ibase].ds_6 = 5.0*bend_bnd_base[ibase].s_6;

/*==========================================================================*/
/* DONE */

       cfree(pot_typ);

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void assign_base_bond(DATA_BASE_ENTRIES *bond_base,int num_base,int *ifound,
                      BOND *bond,int *isearch,int nsearch,int icon_flag)
/*=======================================================================*/
{ /*begin routine */
  /*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,n;

  if(icon_flag==0){
    n = bond->ntyp_pow;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        bond->c_0[i] = bond_base[ibase].c_0;
        bond->c_1[i] = bond_base[ibase].c_1;
        bond->c_2[i] = bond_base[ibase].c_2;
        bond->c_3[i] = bond_base[ibase].c_3;
        bond->c_4[i] = bond_base[ibase].c_4;
        bond->c_5[i] = bond_base[ibase].c_5;
        bond->c_6[i] = bond_base[ibase].c_6;
        bond->dc_0[i] = bond_base[ibase].dc_0;
        bond->dc_1[i] = bond_base[ibase].dc_1;
        bond->dc_2[i] = bond_base[ibase].dc_2;
        bond->dc_3[i] = bond_base[ibase].dc_3;
        bond->dc_4[i] = bond_base[ibase].dc_4;
        bond->dc_5[i] = bond_base[ibase].dc_5;
        bond->dc_6[i] = bond_base[ibase].dc_6;
        bond->eq_pow[i] = bond_base[ibase].eq;
        bond->eq_pow_res[i] = bond_base[ibase].eq_res;
      }/*endif*/
    }/*endfor*/
  }else{
    n = bond->ntyp_con;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        bond->eq_con[i] = bond_base[ibase].eq;
      }/*endif*/
    }/*endfor*/
  }/*endelse*/
/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_base_bend(DATA_BASE_ENTRIES *bend_base,int num_base,int *ifound,
                      BEND *bend,int *isearch,int nsearch,int icon_flag)

/*=======================================================================*/
{ /*begin routine */
/*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,n;

  if(icon_flag==0){
    n = bend->ntyp_pow;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        bend->c_0[i] = bend_base[ibase].c_0;
        bend->c_1[i] = bend_base[ibase].c_1;
        bend->c_2[i] = bend_base[ibase].c_2;
        bend->c_3[i] = bend_base[ibase].c_3;
        bend->c_4[i] = bend_base[ibase].c_4;
        bend->c_5[i] = bend_base[ibase].c_5;
        bend->c_6[i] = bend_base[ibase].c_6;
        bend->dc_0[i] = bend_base[ibase].dc_0;
        bend->dc_1[i] = bend_base[ibase].dc_1;
        bend->dc_2[i] = bend_base[ibase].dc_2;
        bend->dc_3[i] = bend_base[ibase].dc_3;
        bend->dc_4[i] = bend_base[ibase].dc_4;
        bend->dc_5[i] = bend_base[ibase].dc_5;
        bend->dc_6[i] = bend_base[ibase].dc_6;
        bend->s_0[i] = bend_base[ibase].s_0;
        bend->s_1[i] = bend_base[ibase].s_1;
        bend->s_2[i] = bend_base[ibase].s_2;
        bend->s_3[i] = bend_base[ibase].s_3;
        bend->s_4[i] = bend_base[ibase].s_4;
        bend->s_5[i] = bend_base[ibase].s_5;
        bend->s_6[i] = bend_base[ibase].s_6;
        bend->ds_0[i] = bend_base[ibase].ds_0;
        bend->ds_1[i] = bend_base[ibase].ds_1;
        bend->ds_2[i] = bend_base[ibase].ds_2;
        bend->ds_3[i] = bend_base[ibase].ds_3;
        bend->ds_4[i] = bend_base[ibase].ds_4;
        bend->ds_5[i] = bend_base[ibase].ds_5;
        bend->ds_6[i] = bend_base[ibase].ds_6;
        bend->eq_pow[i] = bend_base[ibase].eq;
      }/*endif*/
    }/*endfor*/
  }else{
    n = bend->ntyp_con;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        bend->eq_con[i] = bend_base[ibase].eq;
      }/*endif*/
    }/*endfor*/
  }/*endelse*/
/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_base_tors(DATA_BASE_ENTRIES *tors_base,int num_base,int *ifound,
                      TORS *tors,int *isearch,int nsearch,int icon_flag)

/*=======================================================================*/
{ /*begin routine */
/*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,n;

  if(icon_flag==0){
    n = tors->ntyp_pow;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        tors->c_0[i] = tors_base[ibase].c_0;
        tors->c_1[i] = tors_base[ibase].c_1;
        tors->c_2[i] = tors_base[ibase].c_2;
        tors->c_3[i] = tors_base[ibase].c_3;
        tors->c_4[i] = tors_base[ibase].c_4;
        tors->c_5[i] = tors_base[ibase].c_5;
        tors->c_6[i] = tors_base[ibase].c_6;
        tors->dc_0[i] = tors_base[ibase].dc_0;
        tors->dc_1[i] = tors_base[ibase].dc_1;
        tors->dc_2[i] = tors_base[ibase].dc_2;
        tors->dc_3[i] = tors_base[ibase].dc_3;
        tors->dc_4[i] = tors_base[ibase].dc_4;
        tors->dc_5[i] = tors_base[ibase].dc_5;
        tors->dc_6[i] = tors_base[ibase].dc_6;
        tors->s_0[i] = tors_base[ibase].s_0;
        tors->s_1[i] = tors_base[ibase].s_1;
        tors->s_2[i] = tors_base[ibase].s_2;
        tors->s_3[i] = tors_base[ibase].s_3;
        tors->s_4[i] = tors_base[ibase].s_4;
        tors->s_5[i] = tors_base[ibase].s_5;
        tors->s_6[i] = tors_base[ibase].s_6;
        tors->ds_0[i] = tors_base[ibase].ds_0;
        tors->ds_1[i] = tors_base[ibase].ds_1;
        tors->ds_2[i] = tors_base[ibase].ds_2;
        tors->ds_3[i] = tors_base[ibase].ds_3;
        tors->ds_4[i] = tors_base[ibase].ds_4;
        tors->ds_5[i] = tors_base[ibase].ds_5;
        tors->ds_6[i] = tors_base[ibase].ds_6;
        tors->eq_pow[i] = tors_base[ibase].eq;
      }/*endif*/
    }/*endfor*/
  }else{
    n = tors->ntyp_con;
    for(i=1;i<=n;i++){
      if(ifound[i] > 0 && isearch[i]==nsearch){
        ibase = ifound[i];
        tors->eq_con[i] = tors_base[ibase].eq;
      }/*endif*/
    }/*endfor*/
  }/*endelse*/
/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void assign_base_bend_bnd(DATA_BASE_ENTRIES *bend_bnd_base,int num_base,
                      int *ifound,BEND_BND *bend_bnd,int *isearch,int nsearch,
                      int icon_flag,BUILD_INTRA *build_intra)
{ /*begin routine */
  /*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,ipure,n,ntyp,ntyp_pure,iii;

  n = bend_bnd->ntyp;
  for(i=1;i<=n;i++){
    if(ifound[i] > 0 && isearch[i]==nsearch){
      ibase = ifound[i];
      ipure = bend_bnd_base[ibase].ipure;
      build_intra->ibend_bnd_typ_pure[i] = ipure;
      bend_bnd->cbend_0[i] = bend_bnd_base[ibase].c_0;
      bend_bnd->cbend_1[i] = bend_bnd_base[ibase].c_1;
      bend_bnd->cbend_2[i] = bend_bnd_base[ibase].c_2;
      bend_bnd->cbend_3[i] = bend_bnd_base[ibase].c_3;
      bend_bnd->cbend_4[i] = bend_bnd_base[ibase].c_4;
      bend_bnd->cbend_5[i] = bend_bnd_base[ibase].c_5;
      bend_bnd->cbend_6[i] = bend_bnd_base[ibase].c_6;
      bend_bnd->dcbend_0[i] = bend_bnd_base[ibase].dc_0;
      bend_bnd->dcbend_1[i] = bend_bnd_base[ibase].dc_1;
      bend_bnd->dcbend_2[i] = bend_bnd_base[ibase].dc_2;
      bend_bnd->dcbend_3[i] = bend_bnd_base[ibase].dc_3;
      bend_bnd->dcbend_4[i] = bend_bnd_base[ibase].dc_4;
      bend_bnd->dcbend_5[i] = bend_bnd_base[ibase].dc_5;
      bend_bnd->dcbend_6[i] = bend_bnd_base[ibase].dc_6;
      bend_bnd->sbend_0[i] = bend_bnd_base[ibase].s_0;
      bend_bnd->sbend_1[i] = bend_bnd_base[ibase].s_1;
      bend_bnd->sbend_2[i] = bend_bnd_base[ibase].s_2;
      bend_bnd->sbend_3[i] = bend_bnd_base[ibase].s_3;
      bend_bnd->sbend_4[i] = bend_bnd_base[ibase].s_4;
      bend_bnd->sbend_5[i] = bend_bnd_base[ibase].s_5;
      bend_bnd->sbend_6[i] = bend_bnd_base[ibase].s_6;
      bend_bnd->dsbend_0[i] = bend_bnd_base[ibase].ds_0;
      bend_bnd->dsbend_1[i] = bend_bnd_base[ibase].ds_1;
      bend_bnd->dsbend_2[i] = bend_bnd_base[ibase].ds_2;
      bend_bnd->dsbend_3[i] = bend_bnd_base[ibase].ds_3;
      bend_bnd->dsbend_4[i] = bend_bnd_base[ibase].ds_4;
      bend_bnd->dsbend_5[i] = bend_bnd_base[ibase].ds_5;
      bend_bnd->dsbend_6[i] = bend_bnd_base[ibase].ds_6;
      bend_bnd->eq_bend[i] = bend_bnd_base[ibase].eq;
      if(ipure==0){
        bend_bnd->cbond_0[i] = bend_bnd_base[ibase].cb_0;
        bend_bnd->cbond_1[i] = bend_bnd_base[ibase].cb_1;
        bend_bnd->cbond_2[i] = bend_bnd_base[ibase].cb_2;
        bend_bnd->cbond_3[i] = bend_bnd_base[ibase].cb_3;
        bend_bnd->cbond_4[i] = bend_bnd_base[ibase].cb_4;
        bend_bnd->cbond_5[i] = bend_bnd_base[ibase].cb_5;
        bend_bnd->cbond_6[i] = bend_bnd_base[ibase].cb_6;
        bend_bnd->dcbond_0[i] = bend_bnd_base[ibase].dcb_0;
        bend_bnd->dcbond_1[i] = bend_bnd_base[ibase].dcb_1;
        bend_bnd->dcbond_2[i] = bend_bnd_base[ibase].dcb_2;
        bend_bnd->dcbond_3[i] = bend_bnd_base[ibase].dcb_3;
        bend_bnd->dcbond_4[i] = bend_bnd_base[ibase].dcb_4;
        bend_bnd->dcbond_5[i] = bend_bnd_base[ibase].dcb_5;
        bend_bnd->dcbond_6[i] = bend_bnd_base[ibase].dcb_6;
        bend_bnd->eq_bond[i] = bend_bnd_base[ibase].eqb;
      }/*endif*/
    }/*endif*/
  }/*endfor*/
/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void assign_base_onfo(DATA_BASE_ENTRIES *onfo_base,int num_base,int *ifound,
                      ONFO *onfo,int *isearch,int nsearch,int icon_flag)
{ /*begin routine */
  /*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,n;

  n = onfo->ntyp;
  for(i=1;i<=n;i++){
    if(ifound[i] > 0 && isearch[i]==nsearch){
      ibase = ifound[i];
      onfo->feps[i] = onfo_base[ibase].feps;
      onfo->s6[i] = onfo_base[ibase].s6;
      onfo->sc[i] = onfo_base[ibase].sc;
    }/*endif*/
  }/*endfor*/

/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void assign_base_grpbnd(DATA_BASE_ENTRIES *bond_base,int num_base,int *ifound,
                        GRP_BOND_CON *grp_bond_con,
                        GRP_BOND_WATTS *grp_bond_watts,
                        int *isearch,int nsearch,
                        int icon_flag)
/*=======================================================================*/
{ /*begin routine */
  /*=======================================================================*/
  /*           Local Variables                                             */
  int i,ibase,n;
  int j,k;
  double **eq_21 = grp_bond_con->eq_21;
  double **eq_23 = grp_bond_con->eq_23;      
  double **eq_33 = grp_bond_con->eq_33;      
  double **eq_43 = grp_bond_con->eq_43;
  double **eq_46 = grp_bond_con->eq_46;
  double **c_0_watt_33 = grp_bond_watts->c_0_33;
  double **c_1_watt_33 = grp_bond_watts->c_1_33;
  double **c_2_watt_33 = grp_bond_watts->c_2_33;
  double **c_3_watt_33 = grp_bond_watts->c_3_33;
  double **c_4_watt_33 = grp_bond_watts->c_4_33;
  double **c_5_watt_33 = grp_bond_watts->c_5_33;
  double **c_6_watt_33 = grp_bond_watts->c_6_33;
  double **dc_0_watt_33 = grp_bond_watts->dc_0_33;
  double **dc_1_watt_33 = grp_bond_watts->dc_1_33;
  double **dc_2_watt_33 = grp_bond_watts->dc_2_33;
  double **dc_3_watt_33 = grp_bond_watts->dc_3_33;
  double **dc_4_watt_33 = grp_bond_watts->dc_4_33;
  double **dc_5_watt_33 = grp_bond_watts->dc_5_33;
  double **dc_6_watt_33 = grp_bond_watts->dc_6_33;
  double **eq_watt_33   = grp_bond_watts->eq_33;
  double *cos_thet0_2 = grp_bond_watts->cos_thet0_2;
  double *sin_thet0_2 = grp_bond_watts->sin_thet0_2;

  switch(icon_flag){
    case 1: /* Group 21 [1-2,1-3]                  */
         n = grp_bond_con->ntyp_21;
         for(i=1;i<=n;i++){
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             eq_21[1][i] = bond_base[ibase].eq;
/*      printf("%d %d %g\n",ifound[i],i,bond_base[ibase].eq);*/
           }/*endif*/
         }/*endfor*/
         break;
    case 2: /* Group 23 [1-2,1-3]                  */
         n = 2*grp_bond_con->ntyp_23;
         j=0;k=0;
         for(i=1;i<=n;i++){
           if(((i-1)%2)==0){j++;k=0;}k++;
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             eq_23[k][j] = bond_base[ibase].eq;
           }/*endif*/
         }/*endfor*/
         break;
    case 3: /* Group 33 [1-2,1-3,2-3]              */
         n = 3*grp_bond_con->ntyp_33;
         j=0;k=0;
         for(i=1;i<=n;i++){
           if(((i-1)%3)==0){j++;k=0;}k++;
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             eq_33[k][j] = bond_base[ibase].eq;
/*      printf("%d %d %d %d %g\n",ifound[i],i,j,k,bond_base[ibase].eq);*/
           }/*endif*/
         }/*endfor*/
         break;
    case 4: /* Group 43 [1-2,1-3,2-3]              */
         n = 3*grp_bond_con->ntyp_43;
         j=0;k=0;
         for(i=1;i<=n;i++){
           if(((i-1)%3)==0){j++;k=0;}k++;
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             eq_43[k][j] = bond_base[ibase].eq;
           }/*endif*/
         }/*endfor*/
         break;
    case 5: /* Group 46 [1-2,1-3,1-4,2-3,2-4,3-4]  */
         n = 6*grp_bond_con->ntyp_46;
         j=0;k=0;
         for(i=1;i<=n;i++){
           if(((i-1)%6)==0){j++;k=0;}k++;
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             eq_46[k][j] = bond_base[ibase].eq;
           }/*endif*/
         }/*endfor*/
         break;
    case 6: /* Group 33 [1-2,1-3,2-3]              */
         n = 3*grp_bond_watts->ntyp_33;
         j=0;k=0;
         for(i=1;i<=n;i++){
           if(((i-1)%3)==0){j++;k=0;}k++;
           if(ifound[i] > 0 && isearch[i]==nsearch){
             ibase = ifound[i];
             c_0_watt_33[k][j] = bond_base[ibase].c_0;
             c_1_watt_33[k][j] = bond_base[ibase].c_1;
             c_2_watt_33[k][j] = bond_base[ibase].c_2;
             c_3_watt_33[k][j] = bond_base[ibase].c_3;
             c_4_watt_33[k][j] = bond_base[ibase].c_4;
             c_5_watt_33[k][j] = bond_base[ibase].c_5;
             c_6_watt_33[k][j] = bond_base[ibase].c_6;
             dc_0_watt_33[k][j] = bond_base[ibase].dc_0;
             dc_1_watt_33[k][j] = bond_base[ibase].dc_1;
             dc_2_watt_33[k][j] = bond_base[ibase].dc_2;
             dc_3_watt_33[k][j] = bond_base[ibase].dc_3;
             dc_4_watt_33[k][j] = bond_base[ibase].dc_4;
             dc_5_watt_33[k][j] = bond_base[ibase].dc_5;
             dc_6_watt_33[k][j] = bond_base[ibase].dc_6;
             eq_watt_33[k][j] = bond_base[ibase].eq;
             if(k==3){
               eq_watt_33[k][j] *= (BOHR*M_PI/180.0);
               cos_thet0_2[j] = cos(0.5*eq_watt_33[k][j]);
               sin_thet0_2[j] = sin(0.5*eq_watt_33[k][j]);
               c_1_watt_33[k][j] /= BOHR;
               c_2_watt_33[k][j] /= BOHR*BOHR;
               c_3_watt_33[k][j] /= BOHR*BOHR*BOHR;
               c_4_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR;
               c_5_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR*BOHR;
               c_6_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR*BOHR*BOHR;
               dc_1_watt_33[k][j] /= BOHR;
               dc_2_watt_33[k][j] /= BOHR*BOHR;
               dc_3_watt_33[k][j] /= BOHR*BOHR*BOHR;
               dc_4_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR;
               dc_5_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR*BOHR;
               dc_6_watt_33[k][j] /= BOHR*BOHR*BOHR*BOHR*BOHR*BOHR;
 	     }/*endif*/
           }/*endif*/
         }/*endfor*/
         break;
  }/*end switch*/
/*--------------------------------------------------------------------------*/
} /* end routine */
/*==========================================================================*/


