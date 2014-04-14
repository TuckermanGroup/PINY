 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Structures: typedefs_class.h                           */
 /*                                                                          */
 /* The include file with the typedefs of all the system structures          */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/

/*==========================================================================*/
/*                  Classical Communication variables                       */
/*                                                                          */

/*==========================================================================*/



typedef struct class_plimpton_ind {
  /*---------------------------*/
  /* Processor Numbers and IDs */
  
                
   int num_proc_source;        /* Num: Number of classical force level
                                       processors                         */
   
   int num_proc_target;        /* Num: Number of classical force level
                                       processors                         */
   int myid_source;            /* Num: Rank in the force communicator     */
   int myid_target;            /* Num: Rank in the force communicator     */
  int natoms_tot;             /* Num: Total # of particles                */
 
   /* Arrays with size O(N) */

   int *index_src;           /* Lst : Source atm index array      */
   int *index_src_even;      /* Lst : Even Source atm index array */
   int *index_src_odd;       /* Lst : Odd  Source atm index array */
   int *index_trg;           /* Lst : Target atm index array      */
   int *index_trg_even;      /* Lst : Even Target atm index array */
   int *index_trg_odd;       /* Lst : Odd Target atm index array  */

   /* Arrays with size O(numprocs) */

   int *displs_atm_trg,*counts_atm_trg; /* Lst: Displacements and 
                                         *  recieve counts of atoms 
                                         * into the recieve list  */
   int *displs_atm_src,*counts_atm_src; /* Lst: Displacements and recieve 
                                         * counts of atoms into 
                                         * the recieve list       */

   int *src_odd_end;         /* Lst : End indices in source_odd array */
   int *src_even_end;        /* Lst : End indices in source_even array*/
   int *trg_odd_end;         /* Lst : End indices in target_odd array */
   int *trg_even_end;        /* Lst : End indices in target_even array*/
   int *src_odd_begin;       /* Lst : Begin indices in source_odd array */
   int *src_even_begin;      /* Lst : Begin indices in source_even array*/
   int *trg_odd_begin;       /* Lst : Begin indices in target_odd array */
   int *trg_even_begin;      /* Lst : Begin indices in target_even array*/

   int my_src_odd_counts;         /* Lst : # indices in my_source_odd array */
   int my_src_even_counts;        /* Lst : # indices in my_source_even array*/
   int my_trg_odd_counts;         /* Lst : # indices in my_target_odd array */
   int my_trg_even_counts;        /* Lst : # indices in my_target_even array*/
   int *my_src_odd;               /* Lst : copy from src_ind_odd only myproc*/
                                  /* chunk into my_source_odd array         */
   int *my_src_even;              /* Lst : same for my_source_even array*/
   int *my_trg_odd;               /* Lst : same for my_target_odd array */
   int *my_trg_even;              /* Lst : same for my_target_even array*/



   int **jmax;               /* Lst : Adjusted loop bounds from function calls*/
   int *jmin;                /* Lst : Adjusted loop bounds from function calls*/
   int **imax;               /* Lst : Adjusted loop bounds from function calls*/

  /*-----------------------------*/
  /* Thermostats parallelization */

  /*---------------*/
  /* Communicators */

   MPI_Comm source;            /* Num: Force level communicator-source    */
   MPI_Comm target;            /* Num: Force level communicator-target    */

  

} CLASS_PLIMPTON_IND;

typedef struct class_comm_forc_pkg {
   int num_proc;               /* Num: Number of classical force level
                                       processors                         */
   int myid;                   /* Num: Rank in the force communicator     */
   int myatm_start,myatm_end;  /* Num: Start and ending index for atoms on 
                                       processors   */
   int mytherm_start,mytherm_end; /* Num: Start and end index of thermostat 
                                        on my processor    */
   int *displs_atm,*recv_count_atm; /* Lst: Displacements and recieve 
                                       counts of atoms into the recieve list*/
   int *displs_therm,*recv_count_therm; /* Map: Displacements of nhc's of 
                                          shared nhc's.  Number of 
                                          shared nhc's to recieve   
                                     Lth: np_forc    */
   MPI_Comm comm;              /* Num: Force level communicator           */
   MPI_Comm world;             /* Num: Global communicator                */
   double dbl_num_proc;        /* Num: num_proc as a double               */

  CLASS_PLIMPTON_IND plimpton_ind; /* holds everything related to creating */
                                     /* the plimpton indices */
} CLASS_COMM_FORC_PKG;

/*==========================================================================*/
/*                  Classical atoms/nuclei                                  */
/*             {Variables needed for mem allocation:                        */
/*                   nmol_typ,natm_typ,nmol_tot,                            */
/*                   natm_tot,atyplen,moltyplen,                            */
/*                   natm_nmol_nmol_typ,natm_typ_nmol_nmol_typ }            */
/*                                                                          */

typedef struct clatoms_info{
  int natm_tot;                /* Num: # of atm                          */
  int pi_beads;                /* Num: # of path integral beads          */
  int pi_beads_proc;
  int pi_beads_proc_st;        /* Num: Starting bead on this processor */
  int pi_beads_proc_end;       /* num: Ending bead on this processor   */
  int nfree;                   /* Num: # degrees of freedom              */
  int nfree_pimd;              /* Num: # degrees of freedom              */
  int nchrg;                   /* Num: # chrged atms       Lth:<natm_tot */
  int pi_beads_full_ter;       /* Num: full beads per break-up           */
  int pi_beads_res_ter;        /* Num: res_ter beads per break-up        */ 
  int pi_beads_full_tra;       /* Num: full_tra beads per res_ter break-up*/
  int pi_beads_res_tra;        /* Num: res_tra beads per full_tra break-up*/
  int cg_reset_flag;           /* Opt: Reset CG flag                     */
  int natm_mall;    
  int nchrg_mall;
  int hess_calc;               /* Num: Calculate the atomic hessian     */
  int natm_proc;               /* Num: # of atms per processor */
  int N_H2;                    /* Num: number of H2 molecules            */
  int i_H2_start;              /* Num: index of first atom of H2 block   */
  int nab_initio;              /* Num: # ab inito atms     Lth:<=natm_tot*/
  int myatm_start,myatm_end;   /* Num: Start and ending index for atoms on 
                                       processors   */

  double gamma_adb;            /* Num: Adiabaticity parameter            */
  double wght_pimd;            /* Num: bead harmonic RESPA wght          */
  double pi_temperature;       /* Num: Path integral temperature 
                                       (used in harmonic bead interaction) */
  double pi_beads_full_ter_wght;/* Num: full beads wght                   */
  double pi_beads_res_ter_wght; /* Num: res_ter beads wght                */ 
  double pi_beads_full_tra_wght; /* Num: full_tra beads wght              */
  double pi_beads_res_tra_wght;  /* Num: res_tra beads wght               */
  double rcut_spread;

  double mass_sc_fact;          /* Num: classical mass scaling factor
                                        (for minimization control)       */

  int *ip_lab;                 /* Lst: bead types          Lth: pi_beads */
  int *ichrg;                  /* Lst: Indx of chrged atms Lth:<natm_tot */
  int *cp_vlnc_up,*cp_vlnc_dn; /* Lst: cp atoms cp_valence > 0           */
  int *cp_atm_flag;            /* Lst: list of cp atoms                  */
  double *mass;                /* Lst: atm masses          Lth: natm_tot */
  double *q;                   /* Lst: atm chrges          Lth: natm_tot */
  double *alp_pol;             /* Lst: atm polarizability  Lth: natm_tot */
  double *b_neut;              /* Lst: atm scattering fact Lth: natm_tot */ 
  double *roll_sc;             /* Lst: roll scaling factor Lth: natm_tot */
  double *text_atm;            /* Lst: Atomic temperaures  Lth: natm_tot */
  double *text_mol;            /* Lst: mol  temperaures    Lth: nmol_typ */
  double *xold,*yold,*zold;    /* Lst: Atm positions temp  Lth: natm_tot */  
  double *prekf;               /* Lst: Atm prefks for path Lth:natm_tot  */ 
  double *xmod,*ymod,*zmod;    /* Lst: Atm mode bead;       Lth: natm_tot */
  double alp_ewd;              /* Num: Copy of Ewald alpha */

} CLATOMS_INFO;


typedef struct clatoms_tran{
/* Variables needed for path integral transformations                    */
  int pi_beads_proc;
  int pi_atm_proc;            /* Num: Number of atm per proc             */
  int pi_atm_proc_use;        /* Num: Number of atm per proc to use      */
  int pi_proc_rem;            /* Num: Remainder of natm_tot/numproc      */
  int nwork;                  /* Num: Size of FFT work space             */
  int *sendcounts,*recvcounts;
  int *senddspls,*recvdspls;
  int *ifax;                  /* Lst: Cray FFT work array                
                                  Lth: 13 (don't bother asking why)       */

  double *rat1_stag,*rat2_stag;/* Lst: bead convert factor Lth: pi_beads */
  double *path_eig;            /* Lst: List of bead eigenvalues Lth: pi_beads*/
  double *work,*work2;         /* Lst: FFT work space                     
                                  Lth: nwork;                             */
  double *work3,*work4;         /* Lst: FFT work space                     
                                  Lth: nwork;                             */
  double *x_trans,*y_trans,*z_trans;
                              /* Lst:  FFT transforms scratch arrays      
                                 Lth:  natm_tot                           */
  double *x_temp,*y_temp,*z_temp;
                              /* Lst:  FFT transforms scratch arrays      
                                 Lth:  natm_tot                           */
  double *xt_temp,*yt_temp,*zt_temp;
                              /* Lst:  FFT transforms scratch arrays      
                                 Lth:  natm_tot                           */
} CLATOMS_TRAN;

typedef struct clatoms_pos{
  double *mass;                /* Lst: atm masses          Lth: natm_tot */
  double *x,*y,*z;             /* Lst: Atm positions;      Lth: natm_tot */
  double *vx,*vy,*vz;          /* Lst: Atm velocity;       Lth: natm_tot */
  double *fx,*fy,*fz;          /* Lst: Atm force;          Lth: natm_tot */
/* Path integrals only */
  double *fxt,*fyt,*fzt;       /* Lst: Atm virial force;   Lth: natm_tot */
  double *fxm,*fym,*fzm;       /* Lst: Mode force;         Lth: natm_tot */
/*Minimization only */         /* Lst: Particle hessian    
                                               Lth:(natm_tot)*(natm_tot+1)/2*/
  double *hess_xx,*hess_xy,*hess_xz,*hess_yy,*hess_yz,*hess_zz;
} CLATOMS_POS;

/*==========================================================================*/
/*                          Ghosts                                          */

typedef struct ghost_atoms{
  int nghost_tot;              /* Num: # of ghost atoms                  */
  int natm_comp_max;           /* Num: Max # of atoms comprising a ghost */

  int nghost_mall;
  int nghost_old;
  int ncomp_mall;  
  int ncomp_old;  

  int *ighost_map;             /* Map: Lst of indicies of ghost atoms
                                       in atoms list. Ex: The 3rd ghost
                                       is atom number 27. Lth:nghost_tot */
  int *natm_comp;              /* Lst: # of atoms comprising each ghost  */

  int **iatm_comp;             /* Lst: Atoms comprising each ghost       */

  double **coef;               /* Lst: Coefficients of atoms comprising
                                       the ghost                         */
} GHOST_ATOMS;

/*==========================================================================*/
/*                           Maps                                           */

typedef struct atommaps{
  int nmol_typ;                /* Num: # of mol typs                     */
  int nres_typ;                /* Num: # of res typs                     */
  int natm_typ;                /* Num: # of atm typs                     */  
  int nfreeze;                 /* Num: number of frozen atoms  */

  int nfreeze_mall;
  int natm_typ_mall;
  int nres_typ_max;
  int nres_max;    
  int nres_tot;
  int nres_sum;
  int pimd_freez_typ;

  NAME *atm_typ;               /* Lst: Atm types as strgs; Lth: natm_typ */
  NAME *res_typ;               /* Lst: Res types as strgs; Lth: nres_typ */
  NAME *mol_typ;               /* Lst: Mol types as strgs; Lth: nmol_typ */

  int *nmol_jmol_typ;          /* Lst: # of mol of each Mol type  
                                  Lth: nmol_typ                          */
  int *nres_1mol_jmol_typ;     /* Lst: number of residues in 1 molecule
                                  of the jmol_typth  molecule type;         
                                  Lth: nmol_typ  */ 
  int *jatm_jmol_typ_strt;     /* Lst: start index in atm list where jmol_typth
                                  atoms are listed; Lth: nmol_typ */
  int *jres_jmol_typ_strt;     /* Lst: start index in residue lists 
				  where jmol_typth's residues are listed
                                  Lth: nmol_typ */
  int *ires_typ_jres_jmol_typ; /* Map: jth residue of the jth molecule type
				       to index of the residue type
				  Lth: sum(nres_1mol_jmol_typ) see below */
  int *jatm_jres_1mol_jmol_typ_strt;
                               /* Lst: start index in atm list where 
				  the atoms of the jth res
                                  of the 1st molecule of the 
				  jth molecule type is 
				  listed; Lth: sum(nres_1mol_jmol_typ) */
  int *natm_1mol_jmol_typ;     /* Lst: number of atms in 1 molecule
                                  of the jmol_typth  molecule type;         
                                  Lth: nmol_typ           */ 
  int *natm_jres_jmol_typ;     /* Lst: number of atms in the jth residue
                                  of the jmol_typth  molecule type;         
                                  Lth: sum(nres_1mol_jmol_typ) see below   */ 
  int *nfree_1mol_jmol_typ;     /* Lst: number of degrees of freedom
                                   in 1 molecule of the jmol_typth 
                                   molecule type; Lth: nmol_typ           */
  int *nfree_jres_jmol_typ;     /* Lst: number of degrees of freedom
                                   in the jth residue of the jmol_typth 
                                   molecule type; 
				   Lth: sum(nres_1mol_jmol_typ)           */
  int *icons_jmol_typ;     /* Lst: Flag to indicate if there are constraints
                                   in  the jmol_typth molecule type;         
                              Lth: nmol_typ           */
  int *icons_jres_jmol_typ;/* Lst: Flag to indicate if there are constraints
                                   in  the jth residue of the jmol_typth 
			           molecule type;         
                              Lth: sum(nres_1mol_jmol_typ)           */
  int *natm_jmol_typ;          /* Lst: # of atms in each molecule type
                                  Lth: nmol_typ           */
  int *iatm_mol_typ;           /* Map: Atm ind -> mol_typ ind;     
                                  Ex: 25th atm is of 27th mol_typ; 
                                  Lth: natm_tot           */
  int *iatm_res_typ;           /* Map: Atm ind -> res_typ ind;     
                                  Ex: 25th atm is in the 27th res_typ; 
                                  Lth: natm_tot           */
  int *iatm_atm_typ;           /* Map: Atm ind -> atm_typ ind;
                                  Ex: 25th atm is of 27th atm_typ; 
                                  Lth: natm_tot           */
  int *iatm_atm_typ_nl;        /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: 5th cp atm is of 3rd nonlocal atm_typ; 
                                  Lth: number of atoms with nonlocal pseuds */
  int *iatm_atm_typ_nl_rev;    /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: 5th cp atm is of 7th global atm_typ; */

  int *imap_atm_typ_nl;        /* Map: Atm ind -> cp nonlocal atm_typ ind;
                                  Ex: nonlocal type 1 is global type 3 */
  int *imap_atm_typ_nl_gh;     /* Map: Atm ind -> cp gauss-hermite (gh)
                                      nonlocal atm_typ ind;
                                  Ex: gh nonlocal type 3 is gh type 2 */

  int *iatm_mol_num;           /* Map: Atm ind -> mol ind 
                                  Ex: 25th atm is the 27th molecule of 
                                  mol typ specified by iatm_mol_typ;  
                                  Lth: natm_tot           */
  int *iatm_proc_num;           /* Map: Atm ind -> mol ind 
                                  Ex: 25th atm is on the ith processor
                                  Lth: natm_tot           */
  int *iatm_res_num;           /* Map: Atm ind -> res num
                                  Ex: 25th atm is in the 27th residue of 
                                  mol typ specified by iatm_mol_typ;  
                                  Lth: natm_tot           */
  int *ighost_flag;            /* Lst: 1/0 if atoms is/is not a ghost
                                       Lth: natm_tot                     */
  int *freeze_flag;           /* Lst: 1/0 if atoms is/is not frozen
                                       Lth: natm_tot                     */
  int *freeze_map;             /* Map: indices of frozen atoms 
                                       Lth: nfreeze  */
  int *atom_label;           /* Num: label of atom for freezing of length
                                 natm_tot.
                            Options: Standard(0), Backbone(1), Sidechain(2) */ 

  int *cp_atm_lst;          /*Map: indices of cp atoms used in gen_wave*/  
                            /* also used with dual option              */ 
} ATOMMAPS;

/*==========================================================================*/
/*                  Atom nose-hoover chains                                 */
/*     {Variables needed for mem allocation:  num_nhc,len_nhc}              */

typedef struct therm_info{
  int num_nhc;                 /* Num: # of NHC's                     */
  int len_nhc;                 /* Num: Length of NHC's                */
  int nres_nhc;                /* Num: # of RESPA NHC steps           */  
  int nyosh_nhc;               /* Num: # of Yosh NHC steps            */  
  int therm_typ;               /* Num: 1=NHC, 2=GGMT                  */
  int num_nhc_proc;            /* Num: # of NHC's on this proc        */
  int num_nhc_share;           /* Num: # of NHC's shared over processors */
  int mytherm_start,mytherm_end; /* Num: Start and end index of thermostat 
                                        on my processor    */
  double dt_nhc,dti_nhc,wght;  /* Num: NHC respa time steps           */

  int *inhc_x,*inhc_y,*inhc_z; /* Map: Atm index -> to atm NHC;
                                  Lth: natm_tot                       */
  int *itherm_nshare;           /* Map: Number of processors sharing 
                                        this thermostat
                                   Lth: num_nhc                   */
  double *ditherm_nshare_i;         /* Map: Number of processors sharing 
                                        this thermostat inverse
                                   Lth: num_nhc                   */
  int *map_share;              /* Map: Map of shared thermostat's on 
                                       different proc's
                                  Lth: Num_nhc_share        */
  double *text_nhc;            /* Lst: T_ext of NHC
                                  Lth: num_nhc                        */ 
  double *wdti,*wdti2,*wdti4,*wdti8,*wdti16;/* Lst: Yosh steps;  Lth:9        */

  double **mass_nhc,**gkt;     /* Lst: Mass,degs free of NHC's;
                                  Lth: num_nhc x len_nhc              */

} THERM_INFO;

typedef struct therm_pos{
  double **x_nhc,**v_nhc,**f_nhc;/* Lst: Atm NHC pos,vel,for;         
                                    Lth: num_nhc x len_nhc            */
  double x_nhc_tot;
} THERM_POS;

/*==========================================================================*/
/*                 Velocity Resampling info                                 */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */
typedef struct vel_samp_class{
  int ivel_smpl_on;            /* Opt: Periodic atm vel resampl opt   */
  int ivel_scale_on;            /* Opt: Periodic atm vel resampl opt   */
  int nvx_smpl;                /* Num: Freq of atm vel resampling     */
  int nvx_scale;               /* Num: Freq of atm vel rescaling      */
  int nvnhc_smpl;              /* Num: Freq of atm NHC vel resamp     */
  int iseed,iseed2;            /* Num: Random seeds                   */
  double qseed;                /* Num: Real seed for essl ran()       */
} VEL_SAMP_CLASS;

/*==========================================================================*/
/*                  Energy control options                                  */
/*             {Variables needed for mem allocation:}                       */
typedef struct energy_ctrl{
   int pme_on;                        /* Opt: Particle mesh ewald on    */
   int int_res_ter,int_res_tra;       /* Opt: Inter/Intra res on/off    */
   int isep_vvdw;                     /* Opt: Separate out Vander Waals */
                                      /*      and Coulomb energies      */
   int iswit_vdw;                     /* Opt: VdW  smoothly switched    */
   int block_std_on,block_con_on,nblock_min; 
                                      /* Opt: Intra blocking parameters */
   int iget_pe_real_inter_freq;       /* Opt: Freq calc of inter PE     */
   int iget_pe_real_inter,iget_pv_real_inter; 
                                      /* Opt: Calculate inter PE and PV */

   int itime;                         /* Num: Present time step;        */
   int iget_full_inter,iget_res_inter;/* Opt: Full or respa inter       */
   int iget_full_intra,iget_res_intra;/* Opt: Full or respa intra       */
   int iget_full_pimd,iget_res_pimd;  /* Opt: Full or respa intra       */
} ENERGY_CTRL;


/*==========================================================================*/
/*                Neighbor list                                             */
/*             {Variables needed for mem allocation:}                       */
/*                                                                          */

/*---------------------------------------------------------------------------*/
/*             Verlet list stuff                                             */

  typedef struct verlist{
    int iver_init;               /* Opt: 1st verlist fill               */
    int iver_fill,iver_count;    /* Opt: Fill and count options         */
    int nolst_ver_update,lnk_ver_update; 
                                 /* Opt: Update using nolst or lnk list */
    int jver_pad;                /* Num: Padding used in lnk_ver_update */
    int nmem_min_lst;            /* Num: */
    int nver_lst;                /* Num: # of atms in Verlet list       */
    int nver_lst_res;            /* Num: # of atms in RESPA Verlet list */ 
    int nver_lst_now;         /* Num: actual # of atms in Verlet list       */
    int nver_lst_now_res;     /* Num: actual # of atms in RESPA Verlet list */ 
    int *nter;                   /* Lst: # of neighbors of each atm: 
                                  Lth: natm_tot                         
                                  Can be eliminated using joff        */
    list_int *jter;           /* Lst: Indices of neighboring atms.  
                                  Lth: nver_lst                       */
    int *jver_off;               /* Lst: Starting pts of nghbors in jter 
                                  Lth: natm_tot                       */
    int *nter_res;               /* Lst: RESPA nter                        
                                  Lth: natm_tot                          
                                  Can be eliminated using joff_res    */
    list_int *jter_res;       /* Lst: RESPA jter
                                  Lth: nver_lst_res                   */
    int *jver_off_res;           /* Lst: RESPA jver      
                                  Lth: natm_tot                       */
    double mem_safe;                /* Num: */
  }VERLIST;

/*---------------------------------------------------------------------------*/
/*               Lnk list info                                               */

typedef struct lnklist {
    int ilnk_init;               /* Opt: 1st lnklst fill                */
    int lnk_for_odd;             /* Opt: Force the # of divisions along 
                                       each cell edge to be odd       */
    int ncell_div_avg;           /* Num: Input about how many divisions
                                      of each cell edge to use        */
    int lnk_excl_safe;           /* Num: */
    int lnk_vol_safe;            /* Num: */

    int ncell_a,ncell_b,ncell_c; /* Num: Cell divisions along a,b,c     */
    int natm_cell_max;           /* Num: Max # of atms in any cell      */
    int lnk_list_dim;            /* Num: Size of lnkcell                */
    int nshft_lnk;               /* Num: Total # of shifts needed to
                                        find the inter-mol interacts   */
    int nshft_lnk_res;           /* Num:  RESPA nshft_lnk               */


    double  vol_lnk_map;         /* Lst: Lnk volume                     */
    double rcut_max,rcut_max_res,rexl_max; /*Num: Lnk cell cutoffs      */

    list_int *lnk_list;               /* Lst: List of atms in each cell
                                  Lth: lnk_list_dim                   */
    int *ishft_a,*ishft_b,*ishft_c;/* Lst: Vector of each shift (a,b,c)
                                    Lth: nshft_lnk                    */
    int *iexl_chk;               /* Lst: Indicator of if a shift may 
                                      contains exlcs;
                                  Lth: nshft_lnk                      */
    int *ishft_a_res,*ishft_b_res,*ishft_c_res; 
                               /* Lst: RESPA shifts
                                  Num: nshft_lnk_res                  */
    int *iexl_chk_res;           /* Lst: Excl indicator   
                                  Lth: nshft_lnk_res                  */
    int *natm_cell;              /* Lst: # atms in each cell            */

    double *shft_wght;           /* Lst: Weight of each shift;    
                                  Lth: nshft_lnk_res                  */
    double *shft_wght_res;       /* Lst: Weight of each RESPA shift;
                                  Lth: nshft_lnk_res                  */
    double *hmat_lnk,*hmati_lnk; /* Lst: Lnk cell box shape             */
  }LNKLIST;


/*---------------------------------------------------------------------------*/
/*             Brnch_Root stuff                                              */

typedef struct brnch_root {

 int nbrnch_tot;               /* Num: Number of brnchs*/
 int nroot_tot;                /* Num: Number of roots*/
 int natm_add_tot;             /* Num: Number of additions*/
 int nbrnch_root_max;          /* Num: Max Number of brnch's off a root*/
 double brnch_root_dist_max;   /* Num: Max brnch root distance*/
 double cutdif_min_brnch_root; /* Num: MIN(root_root_cut - root_brnch_cut)*/
 double cutdif_min_brnch_brnch;/* Num: MIN(root_root_cut - brnch_brnch_cut)*/
 double cutdif_min_brnch_root_res;/* Num: MIN(root_root_cut - root_brnch_cut)*/
 double cutdif_min_brnch_brnch_res;/* Num: MIN(root_root_cut-brnch_brnch_cut)*/
 double r2_nocheck_dist;       /* Num: pare down no check distance*/
 double r2_nocheck_dist_res;   /* Num: pare down no check distance*/
 double cutskin_bb_min;        /* Num: minimum branch-branch cutskin*/
 double cutskin_bb_min_res;    /* Num: minimum branch-branch cutskin*/
 int *brnch_atm_list; /* Lst: List of brnch atoms
                         Exp: brnch_atm_list[3] = 25; 
                              The third brnch is atm 25 
                         Lth: nbrnch_tot */
 int *brnch_atm_root; /* Lst: List of brnch atom root
                         Exp: brnch_atm_roots[3] = 25; 
                              The third brnch's root is atm 25 
                         Lth: nbrnch_tot */
 int *root_atm_list;  /* Lst: List of root atoms
                         Exp: root_atm_list[3] = 25; 
                              The third root is atm 25 
                          Lth: nroot_tot */
 int *root_atm_map;   /* Map: Map of atoms into the root list 
                          Exp: root_atm_map[25] = 3; 
                              The 25th atom is the third root
                          Lth: natm_tot */
 int *brnch_atm_map;  /* Map: Map of atoms of big list into the branch list 
                          Exp: branch_atm_map[25] = 3; 
                              The 25th atom is the third branch
                          Lth: natm_tot */
 int *nbrnch_of_root_big;/* Lst: The number of brnchs of each atom in big
                                 list including self. Brnch atoms have 0. 
                          Exp: nbrnch_of_root[3] = 2; 
                              The 3rd atom has two brnchs.
                          Lth: nroot_tot */
 int *nbrnch_of_root;/* Lst: The number of brnchs of each root
                             doesn't not included itself.
                          Exp: nbrnch_of_root[3] = 2; 
                              The 3rd root has two brnchs
                          Lth: nroot_tot */
 int **ibrnch_of_root;/* Lst: The list of brnchs of each root
                          Exp: ibrnch_of_root[3][1] = 30; 
                              The first brnch of the 3rd root is atm 30.
                          Lth: nroot_tot x nbrnch_max */
 int *natm_add;       /* Lst: The number of atm's to be added to each atm's
                              neighbor list to correct for the root-brnch 
                              method
                          Exp: natm_add[3] = 2; 
                              The 3rd atm needs 2 atm's added to its list
                          Lth: natm_tot */
 int *iatm_add;       /* Lst: The list of atm's to be added to each atm's 
                              list
                          Exp: iatm_add[2+iatm_add_off[3]] = 27; 
                              The 2nd atm added to the list of the third 
                              atom is atm number 27
                          Lth: natm_add_tot */
 int *iatm_add_off;     /* Lst: The offset list
                          Exp: iatm_add[2+iatm_add_off[3]] = 27; 
                              The 2nd atm added to the list of the third 
                              atom is atm number 27
                          Lth: natm_add_tot */


} BRNCH_ROOT;

/*---------------------------------------------------------------------------*/
/*             Nbr list stuff                                                */

typedef struct nbr_list {

  int iver;                    /* Opt: Verlet list option             */ 
  int ilnk;                    /* Opt: Lnk list option                */
  int nolst;                   /* Opt: No list option                 */     
  int brnch_root_list_opt;     /* Num: Option to shave brnch interactions */
  double *x0,*y0,*z0;          /* Lst: Pos of the atms when the list
                                       is calculated;        
                                  Lth: natm_tot                       */

  VERLIST verlist;
  LNKLIST lnklist;
  BRNCH_ROOT brnch_root;

} NBR_LIST;

/*==========================================================================*/
/*                Atomic interaction info                                   */
/*             {Variables needed for mem allocation:                        */
/*                 nter_typ,nsplin,nsplin_tot       }                       */
/*                                                                          */
typedef struct interact {
  int nsplin;                  /* Num: # pts in spline of interact typ*/
  int iswit_vdw;                     /* Opt: Switch inter shift on/off */
  int dielectric_opt;          /* Flag to turn on distant dependent 
                                  dielectric constant                 */
  int ishave_opt;              /* Num: Option to shave skin in force_npol.c */

  int nter_typ;                /* Num: # inter-atm interaction 
                                       types=natm_typ*(natm_typ+1)/2  */
  int nsplin_tot;              /* Num: nsplin*nter_typ                */
  int ninter_mall;
  int nsplin_mall_tot;

  int *inter_map_index;        /* Num: length ninter maps into spline
                                       arrays for cv0 cdv0 */
  double dielectric_rheal;     /* Num: Dielectric constant healing length */ 
  double dielectric_cut;       /* Num: Distant Dependent Dielectric cutoff*/ 
  double dielectric_eps;       /* Num: Distant Dependent Dielectric constant*/

  double brnch_root_skin;      /* Num: Extra Skin needed to correctly obtain
                                       brnch-brnch pairs > 
                                       cut_root-2*brnch_root_cut           */
  double skin;                 /* Num: Skin depth used to generate ver_list 
                                       skin=skin_comp+brnch_root_skin      */
  double brnch_root_cut;       /* Num: An atom bonded to only one other atom 
                                       with eq_bond < brnch_root_cut is a 
                                       brnch */
  double rheal_res;            /* Num: RESPA healing length           */ 
  double pten_inter_guess;     /* Num: Estimate of intermolecular pressure */
  double pten_kin_guess;       /* Num: Estimate of kinetic pressure */

  double clong;                /* Num: Long range correction to 
                                       C_6/r^6 potentials             */
  double clong_res;            /* Num: RESPA long range correction to 
                                       C_6/r^6 pots                   */
  double spread_now, spread;   /* Num: Spread of the path integral beads */
  double spread_lnk;
  double brnch_root_skin_lnk;
  double cutoff_max;


  double *cutoff;              /* Lst: Inter-atm interact cutoff dist;
                                 Lth: nter_typ                        */ 
  double *cutoff_res;          /* Lst: Inter-atm interact RESPA 
                                       cutoff distance;
                                  Lth: nter_typ                       */ 
  double *cutti;               /* Lst: inner cutoff         
                                  Lth: nter_typ                       */ 
  double *cutskin;             /* Lst: Inter-atm interaction 
                                       cutoff+skin  distance;  
                                  Lth: nter_typ                       */ 
  double *cutskin_res;         /* Lst: Inter-atm interaction 
                                       RESPA cutoff+skin distance;  
                                  Lth: nter_typ                       */ 
  double *cutskin_root;        /* Lst: Inter-atm interaction
                                       cutoff+skin+branch root skin distance;
                                  Lth: nter_typ                       */
  double *cutskin_root_res;    /* Lst: Inter-atm interaction
                                       RESPA cutoff+skin+branch root skin
                                       distance;
                                  Lth: nter_typ                       */
  double *vcut_coul;           /* Lst: Value of coulomb potential at 
                                  interaction cutoff distance;        
                                  Lth: nter_typ                       */ 
  double *rmin_spl;            /* Lst: Min interact distance splined 
                                       of each atm-interact typ;
                                  Lth: nter_typ                       */ 
  double *dr_spl,*dri_spl;     /* Lst: Spacing of data pts used in 
                                       spline of each atm-interact typ
                                  Lth: nter_typ                       */ 
  double *cv0,*cv1,*cv2,*cv3;  /* Lst: Spline coef of atm-interact pot
                                           Lth: nsplin_tot                 */ 
  double *cdv0,*cdv1,*cdv2,*cdv3;/* Lst: Spline of  atm-interact force
                                             Lth: nsplin_tot               */ 
  double *cv0_c, *cv1_c,*cv2_c, *cv3_c;/* Lst: Spline coef of couloumb
                                                       -interact potential
                                                   Lth: nsplin_c      */ 
  double *cdv0_c,*cdv1_c,*cdv2_c,*cdv3_c;/* Lst: Spline coef of 
                                                     couloumb-interact force;
                                                     Lth: nsplin_c      */ 
} INTERACT;

/*==========================================================================*/
/*               Integrator scratch memory                                  */
/*             {Variables needed for mem allocation: }                      */
/*                                                                          */

typedef struct therm_scr {
  double *sc;                  /* Lst: Scaling factor of atm NHCs; 
                                    Lth: num_nhc+1                      */
  double *atm_kin;             /* Lst: Atm KE of each NHC; 
                                  Lth: num_nhc+1                      */
}THERM_SCR;

typedef struct int_scr {

  double v_lnv;                /* Num: Baro vel temp                  */

/*---------------------------------------------------------------------------*/
/*    Only needed when rolling under constant pressure                       */
  double *x,*y,*z;             /* Lst: Atm pos temps;   Lth: natm_tot */
                               /*      needed for min_CG              */
  double *vx,*vy,*vz;          /* Lst: Atm vel. temps;  Lth: natm_tot */
  double *gx,*gy,*gz;          /* Lst: Grad temp        Pts: vx,vy,vz */
  double **x_nhc,**v_nhc,**f_nhc;/* Lst: Atm NHC temps;  
                                    Lth: num_nhcxlen_nhc              */
  double *vgmat;               /* Lst: Baro temps                     */
  double *x_vol_nhc,*v_vol_nhc,*f_vol_nhc;
                               /* Lst: Vol nhc temp;  Lth: len_nhc    */
/*---------------------------------------------------------------------------*/
  double *sc;                  /* Lst: Scaling factor of atm NHCs; 
                                  Lth: num_nhc+1                      */
  double *sc_temp;             /* Lst: Scaling factor of atm NHCs; 
                                  Lth: num_nhc+1                      */
  double *atm_kin;             /* Lst: Atm KE of each NHC; 
                                  Lth: num_nhc+1                      */
  double *step_up,*step_dn;    /* Lst: Step lengths for CG min.       */

  THERM_SCR *therm_scr;

} INT_SCR;

/*==========================================================================*/
/*               Ewald sum scratch memory                                   */
/*             {Variables needed for mem allocation: }                      */
/*                                                                          */
typedef struct ewd_scr {
  double *q;                   /* Lst: Atm charge temp;    Lth: nchrg */
  double *x,*y,*z;             /* Lst: Atm position temps; Lth: nchrg */
  double *fx,*fy,*fz;          /* Lst: Atm force tmps;     Lth: nchrg */
/*---------------------------------------------------------------------------*/
/* Ewald or CP only */
  double *heli,*helr;          /* Lst: Struct factor temp; Lth: nchrg */
  double *cossc,*sinsc;        /* Lst: Struct factor temp; Lth: nchrg */
/*---------------------------------------------------------------------------*/
/* EWD RESPA OR CP */
  double *fx2,*fy2,*fz2;       /* Lst: Atm force tmps;     Lth: nchrg */
/*---------------------------------------------------------------------------*/
/* CP ONLY */
  double *temp;                /* Lst: Temporary variable  Lth: nchrg */
  double *helr_now,*heli_now,*vtemp_now;
                               /* Lst: Temporary variables Lth: nchrg */
} EWD_SCR;

/*==========================================================================*/
/*               Force scratch memory                                       */
/*             {Variables needed for mem allocation: }                      */
/*                                                                          */
typedef struct for_scr {
  int nlen;
  int nsplin_t;                /* Num:  max(nsplin,nsplin_g)          */
  int num_brk;                 /* Num: Current # of break pts         */
  int intact;                  /* Num: Current # of interactions      */

  double wght_lnk;             /* Num: link lst force weight          */
  double wght_ter;             /* Num: Inter force weight             */
  double wght_ter_res;         /* Num: Inter res force weight         */
  double wght_tra;
  double wght_tra_res;
  double wght_pimd;

  int *j_indext,*i_indext;     /* Lst: Atm index temps; 
                                  Lth: nlen                           */
  int *i_indext2;              /* Lst: Atm index temps; 
                                  Lth: nlen                           */
  list_int *j_lnk,*i_lnk,*k_lnk;     /* Lst: Atm index temp for link cell   */
                               /* Lth: lnk_list_dim */
  int *j_index,*i_index;       /* Lst: Atm index temps; Lth: nlen     */

  int *num_brk_i;              /* Lst: Atm force break point temps 
                                       for each break point; Lth:nlen */
  int *iexcl;                  /* Lst: Exclusion temp  Lth: natm_tot  */
  int *index_atm;              /* Lst: temp Lth : natm_typ */

} FOR_SCR;

/*==========================================================================*/
/*                  Particle mesh data                                      */
/*             {Variables needed for mem allocation:}                       */

typedef struct part_mesh{
/* NORMAL */
  int pme_on;                       /*Opt: PME on                         */
  int kmax_pme;                     /*Num: User PME Mesh estimate         */
  int n_interp;                     /*Num: Order of interpolation         */
  int pme_res_on;
  int kmax_pme_res;                 
  int n_interp_res;
  int pme_para_opt;
  int nlen_pme;                     /*Num: Scr lngth                      */

  int nktot_pme;                    /*Num: equal to ewald->nktot          */
  int ngrid_a,ngrid_b,ngrid_c;      /*Num: PME mesh                       */
  int nktot_pme_res;                /*Num: Equal to ewald->nktot_res      */
  int ngrid_a_res,ngrid_b_res,ngrid_c_res;      

  int ninterp_mall;    
  int ninterp_tot_mall;   

/* SCRATCH */

  int *iatemp,*ibtemp,*ictemp;      /*Lst: Lth: nlen_pme                  */
  int *nc,*ioff_c;                  /*Lst: Lth ngrid_c;                   */

  int **igrid_a,**igrid_b,**igrid_c;/*Lst: Lth: ninterp*nlen_pme          */
  int **igrid_now;                  /*Lst: Lth: ninterp*nlen_pme         */

  double ecut;                      /*Num: Cutoff in Hartree              */
  double ecut_res;                  /*Num: Cutoff in Hartree              */
  double *bweight_tot_res;          
  double *bweight_tot;              /*Lst: Lth: nktot                     */
  double *qgrid;                    /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_scr;                /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_tmp_real;           /*Lst: Lth: nktot                     */
  double *qgrid_tmp_imag;           /*Lst: Lth: nktot                     */
  double *frac_a,*frac_b,*frac_c;   /*Lst: Lth:nlen_pme                   */
  double *aj,*rn,*rn1;              /*Lst: Lth: ninterp*nlen_pme          */

  double **ua,**ub,**uc;            /*Lst: Lth:ninterp*nlen_pme           */
  double **mn_a,**mn_b,**mn_c;      /*Lst: Lth: ninterp*nlen_pme          */
  double **dmn_a,**dmn_b,**dmn_c;   /*Lst: Lth: ninterp*nlen_pme          */
  double **qgrid_now;               /*Lst: Lth: ninterp*nlen_pme         */
} PART_MESH;

/*==========================================================================*/
/*               Surface parameters                                         */
/*          {Variables needed for mem allocation: }                         */
/*                                                                          */
typedef struct surface {

  int    isurf_on;                  /* Opt : Surface potential on/off       */
  int    natm_typ;                  /* Num : Number of atom types           */
  int    nsplin_surf;               /* Num : spline points                  */
  int    nsplin_tot;                /* Num : natm_typ*nsplin                */
  
  double zheal;                     /* Num : Healing length                 */
  double surface_height;            /* Num : Height of the surface          */
  double *surface_pot;              /* Lst : surface pot                    */
  double *surface_forc;             /* Lst : surface force                  */
                                    /*     Lth: nsplin_tot                  */
  double *zcut_off;                 /* Lst : surface pot cutoff             */
                                    /*     Lth: natm_typ                    */
  double *dz_spl;                   /* Lst : surface pot dz                 */
                                    /*     Lth: natm_typ                    */
  double *dzi_spl;                  /* Lst : surface pot dz                 */
                                    /*     Lth: natm_typ                    */
  double *zmin_spl;                 /* Lst : surface pot zmin               */
                                    /*     Lth: natm_typ                    */
  NAME   surface_type;              /* Type of the surface                  */

} SURFACE; 

/*==========================================================================*/
/*==========================================================================*/
/*                       The system                                         */

typedef struct class {
  CLATOMS_INFO clatoms_info; 
  CLATOMS_TRAN clatoms_tran; 
  CLATOMS_POS *clatoms_pos; 
  GHOST_ATOMS ghost_atoms;
  ATOMMAPS atommaps;
  THERM_INFO therm_info_class; 
  THERM_INFO therm_info_bead; 
  THERM_POS therm_class; 
  THERM_POS *therm_bead; 
  VEL_SAMP_CLASS vel_samp_class;
  ENERGY_CTRL energy_ctrl;
  NBR_LIST nbr_list;
  INTERACT interact;
  INT_SCR int_scr;
  EWD_SCR ewd_scr;
  FOR_SCR for_scr;
  double tot_memory;
  PART_MESH part_mesh;
  SURFACE surface;
  COMMUNICATE communicate;
  CLASS_COMM_FORC_PKG class_comm_forc_pkg;
} CLASS;





