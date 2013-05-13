 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Structures: typedefs_par.h                             */
 /*                                                                          */
 /* The include file with the typedefs of all the interface structures       */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/


/*==========================================================================*/
/*            System-parameter setup                                        */

typedef struct class_parse{

  int ivx_smpl,ivnhc_smpl;     /* Num: Atm and Atm NHC initial
                                  velocity resampling opt             */
  int ivx_scale,ivnhc_scale;   /* Num: Atm and Atm NHC initial
                                  velocity resampling opt             */
  int kmax_ewd,kmax_res;       /* Num: Max k-vectors                  */
  int istart;                  /* Num: Type of start                  */
  int ishift_pot;              /* Flag to shift potential             */
  int initial_spread_opt;      /* Flag to spread a classical system out to
                                  a path integral system              */
  int zero_com_vel;            /* Flag to zero com velocities         */

  double tau_vol,tau_vol_nhc;  /* Num: Vol and Vol nhc time scales    */
  double tau_nhc_def;          /* Num: Atm nhc default time scale     */

  int *imol_nhc_opt;           /* Lst: Molecule NHC option;
                                  Lth: nmol_typ                       */
  int *ionfo_opt;              /* Lst: Molecule onefour option;
                                  Lth: nmol_typ                       */
  int *ires_bond_conv;         /* Lst: Molecule bonding convention 
                                  Lth: nmol_typ                       */
  int *mol_freeze_opt;         /* Lst: Molecule freezing option 
                                  Lth: nmol_typ                       */
  int *mol_hydrog_mass_opt;    /* Lst: Change hydrogen mass option
                                  Lth: nmol_typ                       */
  int *mol_hydrog_con_opt;     /* Lst: Constrain hydrogen option
                                  Lth: nmol_typ                       */

  double *tau_nhc_mol;         /* Lst: NHC time scale for each mol
                                  Lth:nmol_typ                        */
  double *text_nhc_mol;        /* Lst: NHC temperatures of each mol   
                                  Lth:nmol_typ                        */
  double *mol_hydrog_mass_val; /* Lst: Change hydrogen mass value
                                  Lth: nmol_typ                       */
} CLASS_PARSE;

/*==========================================================================*/
/*                    CP NHC time scales/options setup                      */

typedef struct cp_parse {
  int istart_cp;               /* Num: Type of start for CP           */
  int istate_nhc_opt;          /* Opt: CP NHC option                  */
  int ivc_smpl,ivcnhc_smpl;    /* Num: Coeff and coeff NHC initial
                                       velocity resampling opt        */

  int ivc_scale,ivcnhc_scale;    /* Num: Coeff and coeff NHC initial
                                       velocity resampling opt        */

  double cp_mass_tau_def;      /* Num: CP coeff default time scale    */
  double cp_mass_tau;          /* Num: CP coeff time scale            */
  double cp_tau_nhc_def;       /* Num: CP NHC default time scale      */
  double cp_tau_nhc;           /* Num: CP coeff time scale            */
  double cp_ecut_def, cp_ecut; /* Num: Coeff cutoffs                  */
  double cp_mass_cut_def,cp_mass_cut;/* Num: Coeff mass cutoff for 
                                             g dependent mass opt     */ 
  double cp_ecut_dual_grid_def; /*Num:Coef cutoffs for large sparse grid*/
  double cp_ecut_dual_grid;
} CP_PARSE;

/*==========================================================================*/
/*              Simulation/Molecule setup files                             */

typedef struct filename_parse{
  char  *simname;              /* Chr: Simulation param output file   */
  char  *molsetname;           /* Chr: Molecular set up file          */
  char  *dnamei;               /* Chr: Input dump file                */
  char  *dnameci;              /* Chr: CP Input dump file             */
  char  *input_name;           /* Chr: Input file name (PIMD.INPUTS)  */

  NAME def_vps_name;           /* Chr: Default pseudo pot lib file    */
  NAME user_vps_name;          /* Chr: user defined one               */
  NAME def_inter_name;         /* Chr: Default intermol lib file      */      
  NAME user_inter_name;        /* Chr: user defined one               */
  NAME def_surf_name;          /* Chr: Default surface lib file       */      
  NAME user_surf_name;         /* Chr: user defined one               */

  NAME *def_intra_name;        /* Chr: Default intramol lib file    
                                  Lst: Bondfile,Bendfile,
                                       Torsfile,Onfofile              */
  NAME  *mol_param_name;       /* Chr: Molecular param file names  
                                  Lth: nmol_typ                       */
  NAME  *res_param_name;       /* Chr: Residue param file names  
                                  Lth: nres_max                       */
  NAME  *res_fix_name;         /* Chr: Residue fix file names  
                                  Lth: nres_max                       */
  NAME *user_intra_name;       /* Chr: user defined one             
                                  Lst: Bondfile,Bendfile,
                                      Torsfile,Onfofile               */
  NAME *vps_name;              /* Chr: list of pseudo potential file  */
                               /*  names to be used in gen_wave       */
                               /* Lth: natm_typ                       */
} FILENAME_PARSE;

/*==========================================================================*/
/*                      Free energy setup                                   */

typedef struct free_parse {
/*---------------------------------------------------------------------------*/
/*            Bond  free energy indicies                                     */
  int *imoltyp_bond_free;      /* Num: moltyp index 1st atm in bond   */
  int *imol_bond_free;         /* Num: mol index of 1st atm in bond   */
  int *ires_bond_free;         /* Num: mol index of 1st atm in bond   */
  int *iatm_bond_free;         /* Num: Atm index of 1st atm in bond   */
/*---------------------------------------------------------------------------*/
/*          Bend  free energy indicies                                       */
  int *imoltyp_bend_free;      /* Num: moltyp index 1st atm in bend   */
  int *imol_bend_free;         /* Num: mol index of 1st atm in bend   */
  int *ires_bend_free;         /* Num: mol index of 1st atm in bend   */
  int *iatm_bend_free;         /* Num: Atm index of 1st atm in bend   */
/*---------------------------------------------------------------------------*/
/*         Tors  free energy indicies                                        */
  int *imoltyp_tors_free;      /* Num: moltyp index 1st atm in tors   */
  int *imol_tors_free;         /* Num: mol index of 1st atm in tors   */
  int *ires_tors_free;         /* Num: mol index of 1st atm in tors   */
  int *iatm_tors_free;         /* Num: Atm index of 1st atm in bend   */
/*---------------------------------------------------------------------------*/
/*     Rbar-sig free energy inicies                                          */
  int nbar_bond;
  int *imoltyp_rbar1_free;     /* Num: moltyp index 1st atm           */
  int *imoltyp_rbar2_free;     /* Num: moltyp index 2nd atm           */
  int *imol_rbar1_free;        /* Num: mol index 1st atm              */
  int *imol_rbar2_free;        /* Num: mol index 2nd atm              */
  int *ires_rbar1_free;        /* Num: res index 1st atm              */
  int *ires_rbar2_free;        /* Num: res index 2nd atm              */
  int *iatm_rbar1_free;        /* Num: atm index 1st atm              */    
  int *iatm_rbar2_free;        /* Num: atm index 2nd atm              */
} FREE_PARSE;

/*=========================================================================*/
/*               Spline scratch memory                                     */

typedef struct spline_parse {
  double *rspl_mem,*dspl_mem,*dgspl_mem;
                               /* Lst: Spline temps;
                                  Lth: max(nsplin,nsplin_g)            */
} SPLINE_PARSE;

/*==========================================================================*/
/*          Molecular bonding parameters                                  */
/*                                           C1-O2
                                             |
                                          N1-C0- res bond
                                             |
                                             O1-H1                           */

typedef struct resbond_prm{
  int opt;                                    /* OPT: residue bond is 
                                           con,off, or on              */
  int res1_index,res2_index;                  /* Num: indices of residues
                                                 to be bonded togethr   */

  int iatm_res1_site,iatm_res2_site;          /* Num: Index of atom at
                                                 residue bond site
                                             ex: iatm_res1_site = 27
                                                 C0 is atom 27         */
  int natm_res1_1back,natm_res2_1back;        /* Num: # atoms 1 atom back
                                                from bond site 
                                            ex: There are 3 atoms
                                                1 back from C0         */

  NAME label;                                 /* Chr: Bond label           */

  NAME res1_typ,res2_typ;                     /* Chr: residue types 
                                                to be bonded togethr   */
  NAME res1_site,res2_site;                   /* Chr: Names of sites at
                                                which residues are
                                                to be bonded
                                            ex: res1_site = "C-term"
                                                res2_site = "N-term"   */
  NAME res1_file,res2_file;                   /* Chr: Necessary modifications
                                                to the residues in 
                                                forming this bond are
                                                contained in the
                                                correction files       */

  int improp_res1[5],improp_res2[5];      /* Num: Indices of primary branch
                                                  atom to be used to form
                                                  an improper dihedral 
                                                  (123),(321),(213)... */
  int iatm_res1_1back[MAX_VALENCE1],
      iatm_res2_1back[MAX_VALENCE1];  /* Lst: Indices of atoms 1 atom
                                        back from bond site
                                    ex: iatm_res1_1back[1] = 36
                                        C1 is atom number 36           */
  int natm_res1_2back[MAX_VALENCE1],
      natm_res2_2back[MAX_VALENCE1];  /* Lst: List of number of atoms 2
                                        atoms back from bond site
                                        parsed into number of atoms
                                        1 atom back from atoms 1
                                        atom back from bond site
                                    ex: natm_res1_2back[1]=1
                                        There is 1 atom 2 back from
                                        C0 and 1 back from C1.          */
  int iatm_res1_2back[MAX_VALENCE1][MAX_VALENCE1],
      iatm_res2_2back[MAX_VALENCE1][MAX_VALENCE1];
                                              /* Lst: Indices of atoms 2
                                                atoms back from bond
                                                site, parsed into
                                                list of atoms 1 atom
                                                back from atoms 1
                                                atom back from bond site
                                           ex:  iatm_res1_2back[1][1]=42
                                                O2 is atom 42.         */
} RESBOND_PRM;

typedef struct resbond_parse {

  int nresidue,nresidue_max;    /* Num: Current/max number of residues  */
  int nres_bond,nres_bond_max;  /* Num: Current/max # of res_bonds      */
  int ionfo;                   /* one four option */
  int iconv;                   /* bonding convention */

  RESBOND_PRM *resbond_prm;  /* Lst: Residue bond parameters 
                                  Lth: nres_bond_max                   */
  int *nres_bond1_jres,*nres_bond2_jres;
                               /*  Lst: # times the residue is the
                                    first/second residue specified 
                                   in a residue bond
                                ex:nres_bond1_jres[24]=3
                                   Residue 24 is specified 3 times
                                   as the first residue in a res bond 
                                   Lth: nres_max                      */
  int *res_bond1_index,*res_bond2_index;  
                               /* Map: Map from residues to the 
                                       associated res_bonds
                               ex: res_bond1_index[3+res_bond1_off[24]]=15
                                   The 3rd time residue 24 is the
                                   1st residue specified in a 
                                   residue bond maps to the 15th residue
                                   bond given in the molecular
                                   param file.
                                 Lth: nres_bond_max                       */
  int *res_bond1_off,*res_bond2_off;  
                               /* Lst: Offset of res_bond_ind lsts  
                               ex: res_bond1_off[24]=200
                                   The list of residue bonds in which
                                   residue 24 is specified first
                                   starts at index 200+1 in the
                                   sorted list, res_bond1_index.
                                  Lth: nres_max                       */
} RESBOND_PARSE;

/*==========================================================================*/
/*            NULL interation temporaries                               */

typedef struct null_inter_parse {

  int nbond_nul,nbend_nul;                 /* Num: #s of null interacts */
  int ntors_nul,nonfo_nul;                 /* Num: #s of null interacts */

  int *jbond1_nul,*jbond2_nul;             /* Lst: null bonds Lth:nbond */
  int *jbend1_nul,*jbend2_nul,*jbend3_nul; /* Lst: null bends Lth:ntors */
  int *jtors1_nul,*jtors2_nul,*jtors3_nul,*jtors4_nul;
                                           /* Lst: null torss Lth:nonfo */
  int *jonfo1_nul,*jonfo2_nul;             /* Lst: null onfo  Lth:nonfo */
} NULL_INTER_PARSE;


/*==========================================================================*/
/*      Words in dictionaries                                               */

typedef struct dict_word{
  int iflag,iuset,key_type;    /* Opt: modifiers                      */
  NAME keyword,keyarg;         /* Chr: keyword and its arg            */
  LINE error_mes;              /* Chr: keyword error mess             */
} DICT_WORD;

/*==========================================================================*/
/*      Intra dictionaries                                                  */

typedef struct dict_intra {
      int num_fun_dict;          /* Num: # func key words(FKW)        */
      int num_mol_name_dict;
      int num_atm_dict;
      int num_intra_dict;
      int num_res_name_dict;
      int num_res_def_dict;
      int num_res_bond_dict;

      DICT_WORD *fun_dict;       
      DICT_WORD *atm_dict;
      DICT_WORD *intra_dict;
      DICT_WORD *mol_name_dict;
      DICT_WORD *res_name_dict;
      DICT_WORD *res_def_dict;
      DICT_WORD *res_bond_dict;
      DICT_WORD *word;
} DICT_INTRA;

/*==========================================================================*/
/*                 MOL_PARAMS DICTIONARIES                                  */

typedef struct dict_mol {
  int num_fun_dict;                         /* Num: # FKWs            */
  int num_mol_dict,num_wave_dict;           /* Num: # KWs             */
  int num_mol_bond_dict, num_bond_free_dict;/* Num: # KWs             */
  int num_bend_free_dict,num_tors_free_dict;/* Num: # KWs             */
  int num_user_base_dict,num_def_base_dict; /* Num: # KWs             */
  int num_rbar_free_dict;
  int num_surface_dict;                     /* Num: # KWs             */

  DICT_WORD *fun_dict;                      /* Dct: FKW dict Lth: num */
  DICT_WORD *mol_dict,*wave_dict;           /* Dct: KW dicts Lth: num */
  DICT_WORD *mol_bond_dict, *bond_free_dict;/* Dct: KW dicts Lth: num */
  DICT_WORD *bend_free_dict,*tors_free_dict;/* Dct: KW dicts Lth: num */
  DICT_WORD *user_base_dict,*def_base_dict; /* Dct: KW dicts Lth: num */
  DICT_WORD *rbar_free_dict;
  DICT_WORD *surface_dict;                  /* Dct: surface pot       */ 
} DICT_MOL;


/*==========================================================================*/
/*                 Build/Bond molecules  stuff                              */

typedef struct { NAME atm1;} CVPS;                  /* Chr: atms in VPS type */
typedef struct { NAME atm1,atm2,label;} CBOND;      /* Chr: atms in bond type*/
typedef struct { NAME atm1,atm2,atm3,label;} CBEND; /* Chr: atms in bend type*/
typedef struct { NAME atm1,atm2,atm3,atm4,label;} CTORS;  
typedef struct { NAME atm1,atm2,atm3,atm4,label;} CATM_LAB;  
typedef struct { NAME atm1,atm2,atm3,atm4,label[7];} CGRP_CONS;  
                                                    /* Chr: atms in tors type*/


/*                                          H2 C1-O2
                                            |  |
                          "N-term" --->    -N1-C0- <----- "C-term"
                                               |
                                               O1-H1                        
                                 
                                  Branch indexing convention
                                  ---------------------------

                                  C0:  (C-term 0 0) or (N-term 2 0)
                                  C1:  (C-term 1 0) or (N-term 2 1)
                                  O2:  (C-term 1 1) or too far from N-term
                                  O1:  (C-term 3 0) or (N-term 2 2)
				  H1:  (C-term 3 1) or too from N-term
				  N1   (C-term 2 0) or (Nterm 0 0)
				  H2:  (C-term 2 1) or (Nterm 1 0)

                                  too far is more than 2 atoms  back
                                  Initialization (Null -1 -1)              */

/* Length of this structure is natm_res (number of atoms in a residue)     */

typedef struct bond_site{
  NAME bond_site_name[MAX_BOND_SITE1];/* Chr: Names of bond sites this atom
                                     is involved in                 
                                 ex: bond_site[28].bond_site_name[1]="C-term"
                                     The first bond site in which
                                    C1, atom index 28,
                                    is involved in is called 
                                    "C-term"                       */
  int branch_1[MAX_BOND_SITE1];     /* Num: Branch of atm off bond site 
                                 ex: bond_site[28].branch_1[1]=1
                                     C1 is on first branch of "C-term"
                                 ex: bond_site[37].branch_1[2]=0
                                    N1, atom index 37 is the root
                                    of 2nd bond site "N-term"     */
  int branch_2[MAX_BOND_SITE1]; /* Num: branch of atm off branch_1
                                 ex: bond_site[28].branch_2[1]=0
                                     C1 is root of the first branch
                                    (i.e., it is atom 0 in the branch)
                                 ex: bond_site[40].branch_2[1]=1
                                     O2, atom index 40, is the first
                                    branch off C1, which
                                    is root of the first branch 
                                    off "C-term"                  
                                 ex: bond_site[42].branch_2[1]=0
                                     C0 is root of branch_1 and
                                    branch_2                          */
  int valence;                   /* Num: valence of the atm
                                 ex: bond_site[42].valence=4
                                     C0 has four connections
                                    Note: Chemist standard is
                                          to count double bonds
                                          as 2, we don't do this --
                                          double bonds count as 1. */
 int improper_ind[5];              /* Num: index of order of primary branch
                                        atoms to be used to form an 
                                        improper dihedral (123)or(321)... 
                                        Unspecifed means no improper        */
} BOND_SITE;

/*========================================================================*/
/* Intra molecular builder structure                                      */
typedef struct build_intra {
      BOND_SITE *bond_site;
      int *mask_atm;           /* Stf: Stuff for atoms                */
      int *index_atm;
      int *iatm_ind_chk;
      int *ires_ind_chk;
      int natm_1res_pure_now;
      int natm_tot_max;
      int nfreeze_max;
      int nfreeze_now;
      int natm_1res_now;
      int natm_1res_max;
      int natmind_1res_now;
      int natmind_1res_max;
      int nghost_now;
      int nghost_tot_max;      /* Stf: Stuff for ghosts               */
      int nbond_pow_max;       /* Stf: Stuff for bonds                */
      int nbond_con_max; 
      int nbond_nul_max;
      int nbond_typ_pow_max; 
      int nbond_typ_con_max;
      CBOND *cbond_typ_pow;
      CBOND *cbond_typ_con;
      CBOND *cbond_typ_now;
      int nbend_pow_max;       /* Stf: Stuff for bends                */
      int nbend_con_max; 
      int nbend_nul_max;
      int nbend_typ_pow_max; 
      int nbend_typ_con_max;
      CBEND *cbend_typ_pow;
      CBEND *cbend_typ_con;
      CBEND *cbend_typ_now;
      int ntors_pow_max;       /* Stf: Stuff for tors                 */
      int ntors_con_max; 
      int ntors_nul_max;
      int ntors_typ_pow_max; 
      int ntors_typ_con_max;
      CTORS *ctors_typ_pow;
      CTORS *ctors_typ_con;
      CTORS *ctors_typ_now;
      int nonfo_max;           /* Stf: Stuff for onfos                */
      int nonfo_nul_max;
      int nonfo_typ_max;
      CBOND *confo_typ;
      CBOND *confo_typ_now;
      int nbend_bnd_max;       /* Stf: Stuff for bend_bnd              */
      int nbend_bnd_typ_max; 
      CBEND *cbend_bnd_typ;
      CBEND *cbend_bnd_typ_now;
      int ngrp_33_max;       /* Stf: Stuff for grpcons              */
      int ngrp_watt_33_max;
      int ngrp_21_max; 
      int ngrp_43_max; 
      int ngrp_23_max; 
      int ngrp_46_max;
      int ngrp_typ_21_max; 
      int ngrp_typ_43_max; 
      int ngrp_typ_33_max; 
      int ngrp_typ_watt_33_max; 
      int ngrp_typ_23_max;
      int ngrp_typ_46_max;
      CGRP_CONS *cgrp_typ_now; 
      CGRP_CONS *cgrp_typ_21; 
      CGRP_CONS *cgrp_typ_43; 
      CGRP_CONS *cgrp_typ_33; 
      CGRP_CONS *cgrp_typ_watt_33; 
      CGRP_CONS *cgrp_typ_23; 
      CGRP_CONS *cgrp_typ_46; 
      int natm_typ_max;
      int nres_typ_max;
      char *strip1,*strip2;
      int *ibend_bnd_typ_pure;
      int *ibend_bnd_typ_map;
} BUILD_INTRA;

/*==========================================================================*/
/*  */

typedef struct start_index {
     int nbond_pow;
     int nbond_con; 
     int nbond_nul;
     int nbend_pow; 
     int nbend_con;
     int nbend_nul; 
     int ntors_pow;
     int ntors_con; 
     int nfreeze; 
     int ntors_nul; 
     int nonfo;     
     int nonfo_nul; 
     int nbend_bnd;        
     int natm;      
     int nghost_tot;
     int ngrp_21; 
     int ngrp_33; 
     int ngrp_watt_33; 
     int ngrp_43; 
     int ngrp_23;
     int ngrp_46;

} START_INDEX;

/*==========================================================================*/
/* type defs for search bases */

typedef struct {double eq,eq_res,eqb;
                double c_0,c_1,c_2,c_3,c_4,c_5,c_6;
                double dc_0,dc_1,dc_2,dc_3,dc_4,dc_5,dc_6;
                double s_0,s_1,s_2,s_3,s_4,s_5,s_6;
                double ds_0,ds_1,ds_2,ds_3,ds_4,ds_5,ds_6;
                double cb_0,cb_1,cb_2,cb_3,cb_4,cb_5,cb_6;
                double dcb_0,dcb_1,dcb_2,dcb_3,dcb_4,dcb_5,dcb_6;
                double feps,s6,sc;
                int ityp,ipure;
                double cutti,cutoff,cutoff_res;
                double sig,eps;
                double awill,bwill,c6m,c8m,c10m; 
                double cwill ,rm_swit, c9m; 
                int inter_label;
                int surf_label;
               } DATA_BASE_ENTRIES;

typedef struct { int atm1,atm2,atm3,atm4,label,good;} WILD;

/*  End type defs for search bases */
/*==========================================================================*/



