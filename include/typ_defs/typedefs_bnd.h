/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Structures: typedefs_bnd.h                             */
/*                                                                          */
/* The include file with the typedefs of all the bonded structures          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

/*==========================================================================*/
/*                Bond energy list info                                     */
/*             {Variables needed for mem allocation:                        */
/*                   nbond_ho,nbond_ho_typ,nbond_con,nbond_con_typ}         */
/*                                                                          */
typedef struct bond {
  int npow;                    /* Num: # pow series bonds             */
  int ntyp_pow;                /* Num: # pow series bond typs         */
  int ncon;                    /* Num: # constrained bonds            */
  int ntyp_con;                /* Num: # constrained bond typs        */
  int nbond_pow_mall;    
  int nbond_typ_pow_mall;
  int nbond_con_mall;    
  int nbond_typ_con_mall;
  int block_pow_on, block_con_on;

  int ncon_tot;
  int ncon_max;
  int nblock_pow, nblock_size_pow;
  int nblock_con, nblock_size_con; 

  int *j1_pow,*j2_pow;         /* Lst: Indices of atms in pow bond;
                                  Lth: npow                           */
  int *jtyp_pow;               /* Map: index of bond -> bond type;    
                                  Lth: npow                           */
  int *j1_con,*j2_con;         /* Lst: Indices of atms in cons bond;   
                                  Lth: ncon                           */
  int *jtyp_con;               /* Map: index of bond -> bond type;
                                  Lth: ncon                           */
  int *iblock_pow_size;
  int *iblock_pow_conflict_1,*iblock_pow_conflict_2;
  int *iblock_con_size;
  int *iblock_con_conflict_1,*iblock_con_conflict_2;

  double *eq_pow;              /* Lst: List of eq. bond lgths;        
                                  Lth: ntyp_pow                       */
  double *eq_pow_res;          /* Lst: List of RESPA eq. bond lgths;        
                                  Lth: ntyp_pow                       */
  double *c_0,*c_1,*c_2,*c_3,*c_4,*c_5,*c_6;            
  double *dc_0,*dc_1,*dc_2,*dc_3,*dc_4,*dc_5,*dc_6;        
                               /* Lst: Bond power series coefficients
                                  Lth: ntyp_pow                       */
  double *al_con;              /* Lst: List of cons bond multipliers;    
                                  Lth: ncon                           */
  double *eq_con;              /* Lst: List of eq. bond lgths;  
                                  Lth: ntyp_con                       */
} BOND;

/*==========================================================================*/
/*                Group Bond energy list info                               */
/*             {Variables needed for mem allocation:                        */
/*                num_33, num_23,num_46 }                                   */
/*                                                                          */
typedef struct grp_bond_con {
  int num_33;                  /* Num: # 3 atom - 3bond units     */
  int num_21;                  /* Num: # 3 atom - 3bond units     */
  int num_43;                  /* Num: # 3 atom - 3bond units     */
  int ntyp_33;                 /* Num: 33 types                   */
  int num_23;                  /* Num: # 3 atom - 2bond units     */
  int ntyp_23,ntyp_21,ntyp_43; /* Num: 33 types                   */
  int num_46;                  /* Lst: 4 atom - 6bond units lists */
  int ntyp_46;                 /* Num: 33 types                   */
  int max_iter;                /* Num: maximum iterations             */

  int ngrp_21_mall;
  int ngrp_33_mall;
  int ngrp_43_mall;
  int ngrp_23_mall;
  int ngrp_46_mall;    
  int ngrp_typ_21_mall;
  int ngrp_typ_33_mall;
  int ngrp_typ_43_mall;
  int ngrp_typ_23_mall;
  int ngrp_typ_46_mall;
  int num_21_tot,num_23_tot,num_33_tot,num_43_tot,num_46_tot;
  double tol;                  /* Num: bond_con tol               */

  int *j1_33,*j2_33,*j3_33;    /* Lst: 3 atom - 3bond units lists */
  int *j1_21,*j1_43,*j2_21;    /* Lst: 3 atom - 3bond units lists */
  int *j2_43,*j3_43,*j4_43;    /* Lst: 3 atom - 3bond units lists */
  int *jtyp_33,*jtyp_21,*jtyp_43; /* Map: Unit to unit type          */
  int *j1_23,*j2_23,*j3_23;    /* Lst: 3 atom - 3bond units lists */
  int *jtyp_23;                /* Map: Unit to unit type          */
  int *j1_46,*j2_46,*j3_46,*j4_46;
                               /* Lst: 4 atom - 6bond units lists */
  int *jtyp_46;                /* Map: Unit to unit type          */

  double **al_33,**eq_33;      /* Lst: Multipliers and eqbonds    */
  double **al_23,**eq_23;      /* Lst: Multipliers and eqbonds    */
  double **al_21,**al_43;      /* Lst: Multipliers and eqbonds    */
  double **eq_21,**eq_43;      /* Lst: Multipliers and eqbonds    */
  double **al_46,**eq_46;      /* Lst: Multipliers and eqbonds    */

  int ncons_gather_tot;        /* Num: length of gather list        */
  int ncons_gather;            /* Num: length of proc gather list        */
  int ncons_gather_max;        /* Num: max length of proc gather list    */
  int *ind_gather_fwd;         /* Map: particle list to proc gather list */
  int *ind_gather_bck;         /* Map: gather list to particle list */
  int *ind_gather_scr;         /* Map: gather list to gather list   */

} GRP_BOND_CON;

/*==========================================================================*/
/*                Group Bond Watts                                          */
/*                                                                          */
typedef struct grp_bond_watts {
  int num_33;                  /* Num: # 3 atom - 3bond units     */
  int num_33_tot;                  /* Num: # 3 atom - 3bond units     */
  int ntyp_33;                 /* Num: 33 types                   */

  int ngrp_33_mall;
  int ngrp_typ_33_mall;

  int *j1_33,*j2_33,*j3_33;    /* Lst: 3 atom - 3bond units lists */
  int *jtyp_33;                 /* Map: Unit to unit type          */

  double *cos_thet0_2,*sin_thet0_2;
  double **eq_33;      /* Lst: Multipliers and eqbonds    */
  double **c_0_33,**c_1_33,**c_2_33,**c_3_33;
  double **c_4_33,**c_5_33,**c_6_33;
  double **dc_0_33,**dc_1_33,**dc_2_33,**dc_3_33;
  double **dc_4_33,**dc_5_33,**dc_6_33;

} GRP_BOND_WATTS;

/*==========================================================================*/
/*                 Contraint control parameters                             */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */
typedef struct constrnt{
  int iconstrnt;               /* Opt: Atm constraints operative      */
  int iroll;                   /* Opt: Roll option                    */
  int max_iter;                /* Num: maximum iterations             */
  double tolshake;             /* Num: Atm shake tolerence            */
  double tolratl;              /* Num: Atm rattle tolerence           */
} CONSTRNT;


 /*==========================================================================*/
 /*                Bond free energy info                                     */
 /*             {Variables needed for mem allocation:                        */
 /*                   nhist_bond_free}                                       */
 /*                                                                          */
typedef struct bond_free {
  int num;                     /* Num: # of free energy bonds         */
  int j1,j2;                   /* Num: Indices of atms in free E bonds*/
  int npow;                    /* Num: Bond power                     */
  int nhist;                   /* Num: # of pts in histogram          */

  double fk;                   /* Num: Bond force constant            */ 
  double eq;                   /* Num: Equil bond length              */
  double del;                  /* Num: Bin width of histogram         */
  double rmin;                 /* Num: Min distance in histogram      */
  double rmax;                 /* Num: Max distance in histogram      */

  char *file;                  /* Chr: Bond free energy outputfile    */ 
  double *hist;                /* Lst: Bond free energy histogram  
                                  Lth:  nhist_bond_free               */
} BOND_FREE;

/*==========================================================================*/
/*                Bend energy info                                          */
/*             {Variables needed for mem allocation:                        */
/*                   nbend_ho,nbend_ho_typ,nbend_con,nbend_con_typ}         */
/*                                                                          */
typedef struct bend {
  int npow;                    /* Num: # pow series bends             */
  int ntyp_pow;                /* Num: # pow series bend typs         */ 
  int ncon;                    /* Num: # constrained bends            */
  int ntyp_con;                /* Num: # constrained bend typs        */ 
  int nbend_pow_mall;    
  int nbend_typ_pow_mall;
  int nbend_con_mall;    
  int nbend_typ_con_mall;
  int block_pow_on,block_con_on;

  int ncon_tot;
  int nblock_pow, nblock_size_pow;
  int nblock_con, nblock_size_con;

  int *j1_pow,*j2_pow,*j3_pow; /* Lst: Indices of atms in pow bends;   
                                  Lth: npow                           */
  int *jtyp_pow;               /* Map: index of bend -> bend type;
                                  Lth: npow                           */
  int *j1_con,*j2_con,*j3_con; /* Lst: Indices of atms in cons bend;
                                  Lth: ncon                           */
  int *jtyp_con;               /* Map: index of bend -> bend type;
                                  Lth: ncon                           */
  int *iblock_pow_size;
  int *iblock_pow_conflict_1,*iblock_pow_conflict_2;
  int *iblock_pow_conflict_3;
  int *iblock_con_size;
  int *iblock_con_conflict_1,*iblock_con_conflict_2;
  int *iblock_con_conflict_3;

  double *eq_pow;              /* Lst: List of eq. bend angs;    
                                  Lth: ntyp_pow                       */
  double *c_0,*c_1,*c_2,*c_3,*c_4,*c_5,*c_6;
  double *s_0,*s_1,*s_2,*s_3,*s_4,*s_5,*s_6;            
  double *dc_0,*dc_1,*dc_2,*dc_3,*dc_4,*dc_5,*dc_6;            
  double *ds_0,*ds_1,*ds_2,*ds_3,*ds_4,*ds_5,*ds_6;
                               /* Lst: Bend power series coefficients
                                  Lth: ntyp_pow                       */
  double *al_con;              /* Lst: List of cons bend multipliers;
                                  Lth: ncon                           */
  double *eq_con;              /* Lst: List of bend typ eq bend angs;
                                  Lth: ntyp_con                       */
} BEND;

/*==========================================================================*/
/*                Bend free energy                                          */
/*             {Variables needed for mem allocation:                        */
/*                   nhist_bend_free}                                       */
/*                                                                          */
typedef struct bend_free {
  int num;                     /* Num: # of free energy bends         */
  int j1,j2,j3;                /* Num: Indices of atms in Free E bend */
  int npow;                    /* Num: Bend power                     */
  int nhist;                   /* Num: # of pts in histogram          */

  double fk;                   /* Num: Bend force constant            */ 
  double eq;                   /* Num: Equil bend angle               */
  double del;                  /* Num: Bin width of histogram         */

  char *file;                  /* Chr: Bend free energy out filename  */ 
  double *hist;                /* Lst: Bond free energy histogram
                                  Lth:  nhist                         */
} BEND_FREE;

/*==========================================================================*/
/*                Bend_bnd energy                                           */
/*             {Variables needed for mem allocation:                        */
/*                   nhist_bend_free}                                       */
/*                                                                          */
typedef struct bend_bnd {
  int num;                     /* Num: # pow series bonds             */
  int ntyp;                    /* Num: # pow series bond typs         */
  int nbend_bnd_mall;    
  int nbend_bnd_typ_mall;
  int block_on;
  int nblock,nblock_size;

  int *j1,*j2,*j3;             /* Lst: Indices of atms in pow bond;
                                  Lth: num                            */
  int *jtyp;                   /* Map: index of bond -> bond type;    
                                  Lth: num                            */
  int *iblock_size;
  int *iblock_conflict_1,*iblock_conflict_2;
  int *iblock_conflict_3;

  double *eq_bond;             /* Lst: List of eq. bond lgths;        
                                  Lth: ntyp                           */
  double *eq_bend;             /* Lst: List of eq. bond lgths;        
                                  Lth: ntyp                           */
  double *cbond_0,*cbond_1,*cbond_2,*cbond_3,*cbond_4,*cbond_5,*cbond_6;   
  double *dcbond_0,*dcbond_1,*dcbond_2,*dcbond_3,*dcbond_4,*dcbond_5,*dcbond_6;
                               /* Lst: Bond power series coefficients
                                  Lth: ntyp                           */
  double *cbend_0,*cbend_1,*cbend_2,*cbend_3,*cbend_4,*cbend_5,*cbend_6;
  double *sbend_0,*sbend_1,*sbend_2,*sbend_3,*sbend_4,*sbend_5,*sbend_6;
  double *dcbend_0,*dcbend_1,*dcbend_2,*dcbend_3,*dcbend_4,*dcbend_5,*dcbend_6;
  double *dsbend_0,*dsbend_1,*dsbend_2,*dsbend_3,*dsbend_4,*dsbend_5,*dsbend_6;
                               /* Lst: Bend power series coefficients
                                  Lth: ntyp                           */
} BEND_BND;

/*==========================================================================*/
/*                Torsional energy info                                     */
/*             {Variables needed for mem allocation:                        */
/*                   ntors_ho,ntors_ho_typ,ntors_con,ntors_con_typ,         */
/*                   ntors_pow,ntors_pow_typ,ntors_frq,ntors_frq_typ}       */
/*                                                                          */
typedef struct tors {
  int npow;                    /* Num: # pow series torsions          */
  int ntyp_pow;                /* Num: # pow series torsion typs      */ 
  int nimpr;                   /* Num: # number of improper torsions  */
  int ncon;                    /* Num: # constrained torsions         */
  int ntyp_con;                /* Num: # constrained torsion typ      */

  int ntors_pow_mall;    
  int ntors_typ_pow_mall;
  int ntors_con_mall;    
  int ntors_typ_con_mall;
  int block_pow_on;
  int block_con_on;
  int ncon_tot;

  int nblock_pow,nblock_size_pow;
  int nblock_con,nblock_size_con;
  int *j1_pow,*j2_pow,*j3_pow,*j4_pow; 
                               /* Lst: Indices of atms in pow tors
                                  Lth: npow                           */
  int *jtyp_pow;               /* Map: index of tors -> tors type;
                                  Lth: npow                           */
  int *j1_con,*j2_con,*j3_con,*j4_con; 
                               /* Lst: Indices of atms in con tors;
                                  Lth: ntors_con                      */
  int *jtyp_con;               /* Map: index of tors -> tors type;
                                  Lth: ntors_con                      */
  int *iblock_pow_size;
  int *iblock_pow_conflict_1,*iblock_pow_conflict_2;
  int *iblock_pow_conflict_3,*iblock_pow_conflict_4;
  int *iblock_con_size;
  int *iblock_con_conflict_1,*iblock_con_conflict_2;
  int *iblock_con_conflict_3,*iblock_con_conflict_4;

  double *c_0,*c_1,*c_2,*c_3,*c_4,*c_5,*c_6;
  double *s_0,*s_1,*s_2,*s_3,*s_4,*s_5,*s_6;            
  double *dc_0,*dc_1,*dc_2,*dc_3,*dc_4,*dc_5,*dc_6;            
  double *ds_0,*ds_1,*ds_2,*ds_3,*ds_4,*ds_5,*ds_6;
                               /* Lst:  Tors pow series coeffs
                                  Lth: npow_typ                       */
  double *eq_pow;              /* Lst: List tors type eq. tors angs;  
                                  Lth: npow_typ                       */
  double  *al_con;             /* Lst: List of con tors multipliers;    
                                  Lth: ntors_con                      */
  double  *eq_con;             /* Lst: List of tor typ tors eq angles;
                                  Lth: ntors_con_typ                  */
} TORS;

/*==========================================================================*/
/*                 Torsional free energy                                    */
/*             {Variables needed for mem allocation:                        */
/*                   nhist_tors_free                                        */
/*                                                                          */
typedef struct tors_free {
  int num;                     /* Num: # of free energy tors          */
  int npow;                    /* Num: Tor power                      */
  int nhist;                   /* Num: # of pts in histogram          */

  int j1[3],j2[3],j3[3],j4[3]; /* Num: Indices of atms in free E tors;*/


  double fk;                   /* Num: Tors  force constant           */ 
  double eq[3];                /* Num: Equil tors angle               */
  double del;                  /* Num: Bin width of histogram         */

  char *file;                  /* Chr: Tors free energy output file   */
  double *hist;                /* Lst: Tors free energy histogram    
                                  Lth:  nhist_tors_free               */
  double **hist_2d;            /* Lst: Tors free energy histogram    
                                  Lth: nhist_tors_free*nhist_tors_free*/
} TORS_FREE;

/*==========================================================================*/
/*                 Rbar-sigma free energy                                   */
/*             {Variables needed for mem allocation:                        */
/*                   nfree,nhist_bar,nhist_sig                              */
/*                                                                          */
typedef struct rbar_sig_free {
  int nfree;                   /* Num: # of rbar_sigma free energy bonds */
  int nhist_bar,nhist_sig;     /* Num: # of pts in histogram          */
  int iopt;                    /* Num: 1/0 if rbar-sig is on/off      */
  int *j1,*j2;                 /* Num: atm index in ith bond          */
  double fk_bar;               /* Num: rbar  force constant           */ 
  double fk_sigma;             /* Num: sigma  force constant          */ 
  double eq_bar;               /* Num: Equil mean                     */
  double eq_sigma;             /* Num: Equil std                      */
  double del_bar,del_sig;      /* Num: Bin width of histogram         */
  double rmin,rmax;            /* Num: Min/max distance in histogram  */
  double smin,smax;            /* Num: Min/max distance in histogram  */
  double rnfree;               /* Num: # of bonds                     */
  double **hist;               /* Lst:  Free energy histogram    
                                  Lth:  nhist_sig*nhist_bar           */
  double **hist_rn;          /* Lst: Tors free energy histogram    
                                  Lth: nhist_sig*nhist_bar*nfree*nhist_bar */
  char *file;                  /* Chr: Free energy output file        */
} RBAR_SIG_FREE;

 /*==========================================================================*/
 /*                 One-four info                                            */
 /*             {Variables needed for mem allocation:                        */
 /*                   nonfo,nonfo_typ}                                       */
 /*                                                                          */
typedef struct onfo{
  int num;                     /* Num: # one-four interaction         */
  int ntyp;                    /* Num: # one-four typs                */ 

  int nonfo_mall;
  int nonfo_typ_mall;
  int block_on;
  int nblock, nblock_size;

  int *j1,*j2;                 /* Lst: Indices of atms in 1-4s;     
                                  Lth: num                            */
  int *jtyp;                   /* Map: index of 1-4 -> 1-4 type;
                                  Lth: num                            */
  int *iblock_size;
  int *iblock_conflict_1,*iblock_conflict_2;

  double *feps;                /* Lst: List of  one-four epsilons;
                                  Lth: ntyp                           */
  double *s6;                  /* Lst: List of 1-4 sigma^6's;     
                                  Lth: num                            */
  double *sc;                  /* Lst: List of 1-4 scaling factors;
                                  Lth: ntyp                           */
} ONFO;

/*==========================================================================*/
/*                Ecor energy list info                                     */
/*             {Variables needed for mem allocation:                        */
/*                   necor}                                                 */
/*                                                                          */
typedef struct ecor{
  int nsplin,nsplin_m2;        /* Num: spline coeffs                  */
  int num;                     /* Num: # ewald corrections            */
  int block_on;
  int nblock, nblock_size;
  int nktot_res;
  int necor_mall;

  double alp_ewd;              /* Num: Ewald convergence parameter    */
  double ecut,ecut_res;
  double rmin_spl;
  double rmax_spl;
  double dr_spl,dri_spl;

  int *j1,*j2;                 /* Lst: Indices of atms in ecors;
                                  Lth: num                            */
  int *iblock_size;
  int *iblock_conflict_1,*iblock_conflict_2;

  double *cv0,*cdv0,*cdv0_res;           /* Lst: spline coefs */


} ECOR;

/*==========================================================================*/
/*                Exclusion list stuff                                      */
/*             {Variables needed for mem allocation:                        */
/*                   num_excl_lst}                                          */
/*                                                                          */
typedef struct excl {
  int nlst;                    /* Num: Total # atm pairs in excl list */
  int nlst_root;               /* Num: Total # atm pairs in excl list */
  int num_cp;

  int *num,*num_root;          /* Lst: # of excls of each atm: 
                                  Exp: Atm 27 has 3 exclusions all of
                                       which have  indices < 27;
                                  Lth: natm_tot                        
                                  Can be eliminated using joff        */
  int *j,*j_root;              /* Lst: Indices of excluded atms.   
                                  Lth: num                            */
  int *j_off,*j_off_root;      /* Lst: Starting pts in *j of excls of 
                                      atms 2-natm_tot;
                                  Exp: Excls of atm  27 start at index 
                                      102 in *j
                                  Lth: natm_tot                       */
  int *j1_cp,*j2_cp;  /*Lst: Indices of atms needed for ewald summation */
                      /*  between ab initio and classical atoms         */


} EXCL;

 /*==========================================================================*/
 /*               Intra-molecular potential scratch memory                   */
 /*             {Variables needed for mem allocation:                        */
 /*                   nlen}                                                  */
 /*                                                                          */
typedef struct intra_scr{
  int nlen;                    /* Num: Lnth of temp vectors;          */

  double wght_ter;             /* Num: Intra force weight             */
  double wght_ter_res;         /* Num: Intra Respa force weight       */
  double wght_tra;             /* Num: Intra force weight             */
  double wght_tra_res;         /* Num: Intra Respa force weight       */
  double int_res_inter;

  int *iatm_typ;               /* Lst: Indx temps           Lth: nlen */
  int *jatm_typ;               /* Lst: Indx temps           Lth: nlen */
/*---------------------------------------------------------------------------*/
/* Used in all intra routines                                                */
  double *vpot;                /* Lst: PE temp;             Lth: nlen */
  double *p11,*p12,*p13;       /* Lst: Pres. tensor temp;   Lth: nlen */
  double *p21,*p22,*p23;      
  double *p31,*p32,*p33;  
  double *x1,*y1,*z1;          /* Lst: Atm position  temp;  Lth: nlen */
  double *x2,*y2,*z2;
  double *x3,*y3,*z3;
  double *x4,*y4,*z4; 
  double *x5,*y5,*z5;          /* Lst: Atm position  temp;  Lth: nlen */
  double *x6,*y6,*z6;
  double *x7,*y7,*z7;
  double *x8,*y8,*z8; 
  double *fx1,*fy1,*fz1;       /* Lst: Atm force temp;      Lth: nlen */
  double *fx2,*fy2,*fz2;
  double *fx3,*fy3,*fz3;
  double *fx4,*fy4,*fz4;      
  double *dx12,*dy12,*dz12;    /* Lst: Atm disp temp;       Lth: nlen */
  double *dx13,*dy13,*dz13;
  double *dx23,*dy23,*dz23;
  double *dx42,*dy42,*dz42;
  double *dx43,*dy43,*dz43; 
  double *dx56,*dy56,*dz56;    /* Lst: Atm disp temp;       Lth: nlen */
  double *dx57,*dy57,*dz57;
  double *dx67,*dy67,*dz67;
  double *dx86,*dy86,*dz86;
  double *dx87,*dy87,*dz87; 
/*---------------------------------------------------------------------------*/
/* Used only in bonds bends tors force_npol routines:                        */
  double *c_0,*c_1,*c_2,*c_3,*c_4,*c_5,*c_6;
  double *dc_0,*dc_1,*dc_2,*dc_3,*dc_4,*dc_5,*dc_6;
                               /* Lst: Pow series coef temp;Lth: nlen */
/*---------------------------------------------------------------------------*/
/* Used only in bonds bends tors routines                                    */
  double *ds_0,*ds_1,*ds_2,*ds_3,*ds_4,*ds_5,*ds_6;
  double *s_0,*s_1,*s_2,*s_3,*s_4,*s_5,*s_6;
                               /* Lst: Pow series coef temp;Lth: nlen */
/*---------------------------------------------------------------------------*/
/* Used in constraint routines                                           */
  double *eq;                  /* Lst: Eq. bond length  Pts: s_0      */
  double *m1,*m2;              /* Lst: Atm mass temps   Pts: s_1,s_2  */
  double *vx1,*vy1,*vz1;       /* Lst: Atm vel temps    Pts: s_3-s_5  */
  double *vx2,*vy2,*vz2;       /* Lst: Atm vel temps    Pts: ds_0-ds_2*/
  double *sc_1, *sc_2;         /*Lst scaling factors    Pts: ds_3,ds_4*/
/*---------------------------------------------------------------------------*/
/* Used in onefour routines                                                  */
  double *sc,*feps,*s6;        /* Lst: Onfo temps       Pts:s_0-s_2   */
/*---------------------------------------------------------------------------*/
/* Used in onefour,ecorr,force_npol routines                                 */
  double *q1,*q2;              /* Lst: Chrg temps       Pts:s_3,s_4   */
/*---------------------------------------------------------------------------*/
/* Used in force_npol routines                                               */
  double *rmin, *dr, *cutoff;  /* Lst:                  Pts:ds_0-ds_2 */
  double *dri; 
  double *vcut_coul,*r, *r3;   /* Lst: Dist temps       Pts:ds_3,s_5  */
  double *del_r,*del_rc;       /* Lst: Disp temps       Pts:s_5,ds_6  */
  double *swit,*cut,*spl_tmp;  /* Lst: temps            Pts:dc_4-dc_6 */
  double *cutoff_res;          /* Lst:                  Pts:ds_4      */ 
} INTRA_SCR;

/*==========================================================================*/
/*==========================================================================*/
/*                    Bonded                                                */

typedef struct bonded {
  BOND bond;
  GRP_BOND_CON grp_bond_con;
  GRP_BOND_WATTS grp_bond_watts;
  CONSTRNT constrnt;
  BOND_FREE bond_free;
  BEND bend;
  BEND_FREE bend_free;
  BEND_BND bend_bnd;
  TORS tors;
  TORS_FREE tors_free;
  RBAR_SIG_FREE rbar_sig_free;
  ONFO onfo;
  ECOR ecor;
  EXCL excl;
  INTRA_SCR intra_scr;
} BONDED;

