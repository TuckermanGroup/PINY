 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Structures: typedefs_cp.h                              */
 /*                                                                          */
 /* The include file with the typedefs of all the cp structures              */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/



/*==========================================================================*/
/*                  CP Communication variables                              */
/*                                                                          */


typedef struct cp_comm_state_pkg {
   int num_proc;               /* Num: Number of state level processors   */
   int myid;                   /* Num: Rank in the state communicator     */
   int nstate;                 /* Num: Number of states                   */
   int nstate_max;             /* Num: State padding (avoids Alltoallv)   */
   int ncoef;                  /* Num: Number of coefficients (in a state)*/
   int nstate_proc;            /* Num: Number of states on this processor */
   int nstate_proc_max;        /* Num: Max number of states on a processor*/
   int nstate_ncoef_proc_max;  /* Num: Max number of coefs/state on a 
                                       processor in transposed form       */
   int nstate_ncoef_proc;      /* Num: Number of coef/state on this 
                                       processor in transposed form       */
   int icoef_start;            /* Num: Starting pt of coeficients         */
   MPI_Comm comm;              /* Num: State level communicator           */
   MPI_Comm world;             /* Num: Global communicator           */
   int ioff_therm_norm[3],ioff_therm_tran[3];
                               /* Lst: Thermostat offsets for transposing */
} CP_COMM_STATE_PKG;

/*==========================================================================*/
/*                  CP simulation options                                   */
/*             {Variables needed for mem allocation:}                       */
/*                                                                          */
typedef struct cpopts{

  int cp_lda,cp_lsda;          /* Opt: Lda-lsda opt                   */
  int cp_sic;                  /* Opt: Sic opt                        */   
  int cp_gga;                  /* Opt: Gga opt                        */   
  int cp_nonint;               /* Opt: Non-interacting elec opt       */
  int cp_norb;                 /* Opt: Norb integration               */
  int cp_gauss;                /* Opt: Gaussian integration           */
  int cp_ptens_calc;           /* Opt: Calculate CP pressure tensor   */
  int cp_init_orthog;          /* Opt: Initially orthogonalize WFs    */
  int cp_gs,cp_low,cp_normalize;
                               /* Opt: Orthogonalization option for 
                                       minimization                   */
  int cp_ngrid_skip;           /* Num: How many grid points to skip
                                       when writing out grid-based quantities */
  int cp_hess_calc;            /* Opt: Calculate the diagonal Hessian
                                       for minimization               */
  int zero_cp_vel;             /* Opt: Zero all coefficient velocities */

  int iwrite_coef_binary;       /* binary write option for coef file */
  int iread_coef_binary;        /* binary read option for coef file */
  int icheck_perd_size;        /* Opt: Check atm dists under CBCs */
  int icheck_dual_size;        /* Opt: Check atm dists under CBCs */
  int cp_para_opt;             /* Opt: Type of CP parallelization 
                                        0=hybrid,1=full g space    */

  int cp_becke;                 /* cp GGA flags */
  int cp_pw91x;
  int cp_fila_1x; 
  int cp_fila_2x;
  int cp_pbe_x;
  int cp_revpbe_x;
  int cp_rpbe_x;
  int cp_xpbe_x;
  int cp_brx89;
  int cp_brx2k;
  int cp_lyp;  
  int cp_lypm1;  
  int cp_pw91c;
  int cp_pbe_c;
  int cp_xpbe_c;
  int cp_tau1_c;
  int cp_debug_xc;
  int cp_dual_grid_opt;              /*Opt: on = 1 off = 0 */
  int cp_isok_opt;                   /*Opt: Use Gaussian isokinetic for CP coeffs */

  double te_ext;               /* Num: Temperature of PW coef         */ 
  double tol_edge_dist;       /*  Num: Tolerance on max atom dist under CBC*/
  double tol_coef;             /* Num: cp and cp_wave fcoef tolerance*/

  double *occ_up,*occ_dn;      /* Lst: orbital occupation numbers     
                                   Lth: nstate_up                     */
  double *rocc_sum_up;            /* Lst: reciprocal sums of occ numbers 
                                   Lth: nstate_up*nstate_up           */
  double *rocc_sum_dn;            /* Lst: reciprocal sums of occ numbers 
                                   Lth: nstate_up*nstate_up           */
} CPOPTS;

/*==========================================================================*/
/*                  CP coeffs and states                                    */
/*             {Variables needed for mem allocation:                        */
/*                   ncoef_l,ncoef,(0,ncoef),nstate_up,nstate_dn}         */
/*                                                                          */
typedef struct cpcoeffs_info {
  int pi_beads;                /* Num: # of path integral beads        */
  int nstate_up,nstate_dn;     /* Num: # spin up/dn states             */
  int icmass_unif;             /* Opt: Uniform coeff PW masses         */
  int ks_rot_on;               /* Opt: Perioidic rotation to 
                                       KS states opt                   */
  int n_ks_rot;                /* Num: Freq of rotation to
                                       KS states                       */ 
  int cp_ke_dens_on;           /* Opt: Calculate the electron kinetic
                                       energy density                  */
  int cp_tau_functional;       /* Opt: A tau-dependent functional is 
                                       being used in this calculation  */
  int cp_elf_calc_frq;         /* Num: Compute the electron localization
                                       function evey ... steps         */
  int cp_laplacian_on;         /* Num: Compute the Laplacian of the density */

  int cp_nfree;                 /* Num: Number of fictitious electronic
				  degrees of freedom                   */
  int cp_nfree_up;              /* Num: Number of fictitious electronic
				  degrees of freedom                   */
  int cp_nfree_dn;              /* Num: Number of fictitious electronic
				  degrees of freedom                   */
  int pi_beads_proc;           /* Num: # of path integral beads on 
                                         this processor                */
  int pi_beads_proc_st;        /* Num: Starting bead on this processor */
  int pi_beads_proc_end;       /* num: Ending bead on this processor   */
  int nstate_up_proc,nstate_dn_proc;
                               /* Num: # spin up/dn states on this
                                       processor                       */
  int istate_up_st,istate_dn_st;  /*Num: Starting state on this processor*/
  int istate_up_end,istate_dn_end;/*Num: Ending state on this processor*/
  int cg_reset_flag;           /* Opt: Reset CG flag                   */
  int diis_reset_flag;         /* Opt: Reset DIIS flag                 */

  int ncoef;                   /* Num: # of PW coef on small 
                                       spherical cutoff grid dens_cp_box */  
  int ncoef_l;                 /* Num: # of coef on large 
                                       spherical cutoff grid sparse box  */  
  int ncoef_l_dens_cp_box;       /* Num: # of coef on small dense rho
                                       spherical cutoff cp_grid (hmat_cp)*/  
  /* All on dens_cp_box  */
  int icoef_start_up;          /* Num: Starting pt of coeficients         */
  int nstate_ncoef_proc_up;    /* Num: # coef/state on this 
                                       processor in transposed form       */
  int nstate_ncoef_proc_max_up;/* Num: Max # of coef/state on this 
                                       processor in transposed form       */
  int icoef_start_dn;          /* Num: Starting pt of coeficients         */
  int nstate_ncoef_proc_dn;    /* Num: # coef/state on this 
                                       processor in transposed form       */
  int nstate_ncoef_proc_max_dn;/* Num: Max # of coef/state on this 
                                       processor in transposed form       */

  int *ioff_up,*ioff_dn;       /* Lst: State Off sets in PW vector    */
  int *ioff_upt,*ioff_dnt;     /* Lst: State Off sets transposed PW vector */

  double ecut;                 /* Num: Energy cutoff large sparse grid    */
  double ecut_dens_cp_box;     /* Num: Energy cutoff cp_box               */
  double cp_hess_cut;          /* Num: Hessian cutoff and Hessian tau */
  double pseud_hess_loc;       /* Num: Local PP contribution to Hessian */
  double max_off_diag;         /* Num: maximum off-diagonal overlap
				       matrix element                 */
  double max_diag;             /* Num: maximum diagonal overlap
				       matrix element                 */

  double *cmass;                /* Lst: PW coef forces;
                                  Lth: ncoef                          */
  double k0_up,k0_dn;           /* Num: Fixed CP fictitious kinetic 
                                        energies for use in CP 
                                        isokinetic integration scheme  */
} CPCOEFFS_INFO;

typedef struct cpcoeffs_pos {
 /* All of this is dens_cp_box */
  double *cre_up,*cim_up;          /* Lst: Re/Im part of PW spin up 
                                      Lth: ncoef x nstate_up              */
  double *cre_dn,*cim_dn;          /* Lst: Re/Im part of PW spin dn 
                                      Lth: (0,ncoef)*nstate_dn            */
  double *vcre_up,*vcim_up;        /* Lst: Re/Im part PW spin up 
                                         coef velocities;
                                      Lth: ncoef x nstate_up              */
  double *vcre_dn,*vcim_dn;        /* Lst: Re/Im part of PW
                                           wave spin up coef velocities;
                                      Lth: (0,ncoef)x nstate_dn           */ 
  double *fcre_up,*fcim_up;        /* Lst: Re/Im part PW coeff forces;
                                      Lth: ncoef x nstate_up              */
  double *fcre_dn,*fcim_dn;        /* Lst: Re/Im part PW coeff forces;
                                      Lth: (0,ncoef) x nstate_dn          */
  double *cp_hess_re_up;           /* Lst: CP up diagonal Hessian;
                                      Lth: ncoef                          */
  double *cp_hess_im_up;           /* Lst: CP up diagonal Hessian;
                                      Lth: ncoef                          */
  double *cp_hess_re_dn;            /* Lst: CP dn diagonal Hessian;
                                      Lth: ncoef                          */
  double *cp_hess_im_dn;            /* Lst: CP dn diagonal Hessian;
                                   Lth: ncoef                          */
  double *ksmat_up,*ksmat_dn;      /* Lst: Up/dn Kohn-Sham matricies      */
                                   /* Lth: nstate x nstate                */

  double *ksmat_eig_up;            /* Lst: Up/dn Kohn-Sham eigenvalues    */
  double *ksmat_eig_dn;            /* Lth: nstate                         */
  double ks_offset;                /* Num: Energy correction to KS eigenvalues
                                            so that \sum_i f_i e_i = E_elec */
                                   
  double *norbmat_up,*norbmat_dn;  /* Lst: Up/dn St^{-1/2} matrix 
                                       St^{-1/2}*S*St^{-1/2} = diag_occ
                                      Lth: nstate x nstate                */
  double *norbmati_up,*norbmati_dn;/* Lst: Up/dn St^{1/2} matrix          
                                      Lth: nstate x nstate                */

  double *ovmat_eigv_up;           /* Lst: Up/dn S matrix eig vectors     */
  double *ovmat_eigv_dn;           /* Lth: nstate x nstate                */

  int icoef_form_up,icoef_form_dn;  /* Opt: PW coefs, vels, forces form    */
  int ivcoef_form_up,ivcoef_form_dn;/*      transposed(1) or normal(0)    */ 
  int ifcoef_form_up,ifcoef_form_dn;                 

  int icoef_orth_up,icoef_orth_dn;  /* Opt: PW coefs ortho(1)/nonortho(0)  */
  int ivcoef_orth_up,ivcoef_orth_dn;/* Opt: PW fcoefs ortho(1)/nonortho(0) */
  int ifcoef_orth_up,ifcoef_orth_dn;/* Opt: PW fcoefs ortho(1)/nonortho(0) */
  double max_diag,max_off_diag;     /* Num: max elements of ovlap matrix   */

} CPCOEFFS_POS;

/*==========================================================================*/
/*                  Electronic properties                                   */
/*             {Variables needed for mem allocation:                        */
/*                   nfft2_up,nfft2_up_proc}                                */
/*                                                                          */

typedef struct electronic_properties{
  double *cp_elf_up,*cp_elf_dn; /* Lst: Becke-Edgecomb electronic localization
                                        function
                                        Lth: nfft2                          */
} ELECTRONIC_PROPERTIES;

/*==========================================================================*/
/*                 CP coeff Nose'-Hoover chains                             */
/*             {Variables needed for mem allocation:                        */
/*                   num_c_nhc,len_c_nhc }                                  */
/*                                                                          */
typedef struct cptherm_info{
  int num_c_nhc;               /* Num: # of PW coeff NHC's            */
  int num_c_nhc_proc;          /* Num: # of PW coeff NHC's on proc    */
  int num_c_nhc_norm;          /* Num: # of PW coeff NHC's in normal form */
  int len_c_nhc;               /* Num: Lnth of PW coeff NHC's         */
  int nres_c_nhc;              /* Num: # of PW coeff RESPA NHC steps  */ 
  int nyosh_c_nhc;             /* Num: # of PW coeff Yosh NHC's steps */ 
  int massiv_flag;             /* opt: Massive thermo flag         */
  int istate_nhc_opt;          /* opt: NHC option                     */
  double cp_therm_heat_fact;   /* opt: Sample or Rescale Hot NHCs     */

  double dt_nhc,dti_nhc,wght;  /* Num: NHC respa time steps           */
  double c_gkt_massiv;         /* Num: gkt variable under massive  */
  double cmass_nhc_massiv;     /* Num: nhc mass under massive      */

  int *icmapup_nhc;            /* Map: Coeff index -> coeff NHC index;
                                  Lth: nstate_up              */ 
  int *icmapdn_nhc;            /* Map: Coeff index -> coeff NHC index; 
                                  Lth: nstate_dn                       */ 

  double *wdti,*wdti2,*wdti4,*wdti8,*wdti16;/* Lst: Yosh steps;  Lth:9        */

  double **c_gkt;               /* Lst: PW coeff NHC deg free;
                                  Lth: num_c_nhc x len_c_nhc          */
  double **cmass_nhc;           /* Lst: PW coeff NHC mass; 
                                  Lth: num_c_nhc x len_c_nhc          */  
                               

} CPTHERM_INFO;

typedef struct cptherm_pos{
  double c_nhc_massiv;           /* Num: c_nhc under massive         */

  double **c_nhc;               /* Lst: PW coef NHC pos; 
                                  Lth: num_c_nhc x len_c_nhc          */
  double **vc_nhc;              /* Lst: PW coef NHC vel;     
                                  Lth: num_c_nhc x len_c_nhc          */
  double **fc_nhc;              /* Lst: PW coef NHC for;     
                                  Lth: num_c_nhc x len_c_nhc          */
  int itherm_form_up;           /* Num: Thermostat form flag (massiv) */
  int itherm_form_dn;           /* Num: Thermostat form flag (massiv) */
} CPTHERM_POS;

/*==========================================================================*/
/*  CP constrnt stuff */

typedef struct cpconstrnt{
  double c_tolshake;           /* Num: PW coef shake tolerence        */
  double c_tolratl;            /* Num: PW coef rattle tolerence       */
  double c_tolnorb;            /* Num: PW coef norb tolerence         */

} CPCONSTRNT;

/*==========================================================================*/
/*  More CP/Ewald stuff */

typedef struct cpewald {
  int nktot_sm;                /* Num: # of PW coeffon small sphere
                                     cutoff g-space grid=ncoef  dens_cp_box*/
  int nktot_dens_cp_box;       /* When have 2 boxes  dens_cp_box
                                       Num: # of PW coeffon large sphere
                                       cutoff g-space cp_grid=ncoef_cp_l  */

  int nktot_cp_mall;   
  int nktot_cp_sm_mall;
  int nktot_dens_cp_box_mall;       /* if 2 boxes set by cp box */
  int box_rat;                      /* Box ratio for dual gridding */
                                    

  double dbox_rat;             /* Double value for box ratio */
  double gw_gmin,gw_gmax;
  double gw_gmin_dens_cp_box,gw_gmax_dens_cp_box;

  int *kmax_cp;                     /* Lst: Int cutoff in a,b,c directions */
  int *kmax_cp_dens_cp_box;         /* Lst: cutoff for the small box       */
                                    /*  dual=0 kmax_cp_dens_cp_box=kmax_cp */
                                    /*  dual=1 kmax_cp_dens_cp_box=kmax_cp */
                                    /*         but kmax_cp requires box_rat*/
                                    /*  dual=2 kmax_cp has true values     */

  int *kastr_sm,*kbstr_sm,*kcstr_sm; /* Lst: Small spherically cutoff 
                                             g-vectors; dens_cp_box
                                        Lth: nktot_sm                 */ 
  int *kastr_dens_cp_box,
      *kbstr_dens_cp_box,
      *kcstr_dens_cp_box;              /* Lst: large spherically cutoff 
                                          g-vectors in cp_box;
                                        Lth: nktot_cp_l                 */ 

  int *ibrk1_sm,*ibrk2_sm;     /* Lst: RESPA spherically cutoff 
                                        g-space grid break-pts;
                                  Lth: nktot_sm                       */ 
  int *ibrk1_dens_cp_box,
      *ibrk2_dens_cp_box;             /* Lst: RESPA spherically cutoff 
                                        g-space grid break-pts;
                                        Lth: nktot_sm                       */ 

  double *ak2,*ak2_sm;         /* Square of k vectors on big and small grid */
  double *ak2_dens_cp_box;      /* Square of k vectors on big and small grid */

} CPEWALD;


/*==========================================================================*/
/*                Pseudopotential interaction info                          */
/*             {Variables needed for mem allocation:                        */
/*                 vxctyplen,nsplin_g,nsplin_g_tot,n_ang_max,               */
/*                 num_nl_lst       }                                       */
/*                                                                          */
typedef struct pseudo{

  int n_ang_max;               /* Num: Max # of angular momentum 
                                       channels in an e-atm pseudopot */
  int n_ang_max_gh;            /* Num: Max # of angular momentum 
                                       channels in an e-atm pseudopot 
                                       only Gauss-Hermite type*/
  int n_ang_max_kb;            /* Num: Max # of angular momentum 
                                       channels in an e-atm pseudopot 
                                       KB and Goedecker type*/
  int n_rad_max;               /* Num: Max # radial channels in any
                                                      e-atm pseudopot */
  int nsplin_g;                /* Num: # of g-space spline pts in any 
                                       e-atm interaction channel      */
  int nsplin_g_tot;            /* Num: Total # of g-space spline pts  */
  int num_nl_lst;              /* Num: Total # of atms involved in the
                                       pseudopot angular momentum 
                                       channels =sum np_nl(j)         */
  int nl_cut_on;               /* Opt: Non-local cutoff scheme opt    */

  int n_ang_mall;      
  int nsplin_g_mall;   
  int nvpsnorm_mall;       
  int nlist_mall;       
  int n_interp_pme_dual;       /*Num: Order of interpolation         */
  int natm_typ_nl;             /*Num: number of nonlocal atom types  */
                               /*    this is a subset of iatm_atm_typ*/
  int natm_typ_gh;             /*Num: number of nonlocal atom types  */
  int ngh;                     /* Num: Number of gauss-hermite points for
                                       each nonlocal atom type */

  double gmin_true;            /* Num: Mag of smallest g-vector       */
  double gmin_spl;             /* Num: Min mag of g-vec in spline     */
  double gmax_spl;             /* Num: Max mag of g-vec in spline     */
  double dg_spl;               /* Num: Spacing between g-vectors pts
                                       in spline of e-atm pseudos     */
  double gga_cut;              /* Num: Gradient cutoff value          */
  double nlvps_skin;           /* Num: Nonlocal pseudopotential list
                                       skin length                     */
  double alpha_conv_dual;      /* Num: convergence factor for long     */
                               /*  and short range break up for dual grid*/
                                

  int *n_ang;                  /* Lst: # of angular momentum channels 
                                       in each e-atm pseudopot;
                                  Lth: natm_typ                       */ 
  int *nrad_0,*nrad_1,*nrad_2,*nrad_3; /* Lst: # rad channels in each
                                       l channel in each e-atm pseudopot;
                                  Lth: natm_typ                       */
  int *nrad_max_l;            /* Lst: max # rad channels in each
                                       l channel;
                                  Lth: natm_typ                       */
  int *loc_opt;                /* Lst: Angular momentum channel 
                                       chosen as local; Lth: natm_typ */ 
  int *ivps_label;             /* Lst: Type label of e-atm pseudopots;
                                  Lth: natm_typ                       */ 
  int *np_nl;                  /* Lst: # of atms in each  angular 
                                     momentum  channel except GAUSS-HERMITE
                                  Lth:  (n_ang_max+1)                  */
  int *np_nl_gh;               /* Lst: # of gauss-hermite atms in each 
                                       angular momentum  channel;               
                                  Lth:  (n_ang_max+1)                  */
  int **np_nl_rad_str;         /* Lst : where each rad channel strs    */
  int **np_nl_rad_end;         /* Lst : where each rad channel ends    */
  int *ip_nl;                  /* Lst: index of atms involved in each 
                                       angular momentum channel;
                                  Lth: num_nl_lst<natm_tot*(n_ang_max+1)*/
  int *ip_nl_rev;              /* Lst: index of atms involved in each 
                                       angular momentum channel;
                                  Lth: num_nl_lst<natm_tot*(n_ang_max+1)*/
  int *ip_nl_gh;               /* Lst: index of atms involved in each 
                                       angular momentum channel
                                       GAUSS-HERMITE
                                  Lth: num_nl_lst<natm_tot*(n_ang_max+1)*/
  int *ip_nl_rev_gh;           /* Lst: index of atms involved in each 
                                       angular momentum channel
                                       GAUSS-HERMITE
                                  Lth: num_nl_lst<natm_tot*(n_ang_max+1)*/
  int *map_nl;                 /* Lst: order atoms in decreasing order
                                       of # of radial channels          */

  int np_loc_cp_box;           /* Num: # of atms in small cp box       */ 
  int np_nonloc_cp_box;        /* Num: # of atms in small cp box       */ 
  int np_nonloc_cp_box_kb;     /* Num: # of kleinman-bylander nonlocal atoms*/
  int np_nonloc_cp_box_gh;     /* Num: # of gauss-hermite nonlocal atoms */
  int *ip_loc_cp_box;          /* Lst: index of atms on small cp box   */
                               /* Lth: natm_tot                        */
  int *ip_nonloc_cp_box;       /* Lst: index of atms on small cp box   */
                               /* Lth: natm_cp                         */
  double *vps0,*vps1,*vps2,*vps3;/* Lst: Spline coef of pseudopot 
                                    Lth: nsplin_g_tot                 */ 
  double *dvps0,*dvps1,*dvps2,*dvps3;/* Lst: Spline coef of pseudopot
                                    Lth: nsplin_g_tot                 */
  double *gzvps;               /* Lst: g=0 term of local part of 
                                       e-atm pseudopot; Lth: natm_typ */ 
  double *gzvps0;              /* Lst: g=0 term of l=0 channel of 
                                  e-atm pseudopot; Lth: natm_typ      */ 
  double *vpsnorm;             /* Lst: Norm of KB type non-local 
                                  pseudopot; Lth: natm_typ*(n_ang_max+1)  */
  double *rcut_nl;             /* Lst: Non-local cutoff distance;
                                  Lth: (natm_typ)                      */
  double *rgh,*wgh;            /* Lst: Gauss-Hermite nodes and weights 
                                  Lth: ngh  */
  double *x0w,*y0w,*z0w;       /* Lst: Position of atms when non-local
                                        list, ip_nl, is calculated;
                                  Lth: natm_tot                           
                                  Can be eliminated by writing to 
                                  disk and using vxg,vyg,vzg as a temp 
                                  to read it in                       */
  double *q_pseud;             /* charge associated with pseudo  */
                               /* Lth: natm_typ */

  char *vxc_typ;               /* Chr: Exchange-correlation type 
                                  Lth: MAXWORD                        */
  char *ggax_typ;               /* Chr: GGA-Exchange-correlation type 
                                  Lth: MAXWORD                        */
  char *ggac_typ;               /* Chr: GGA-Exchange-correlation type 
                                  Lth: MAXWORD                        */
} PSEUDO;

/*==========================================================================*/
/*               CP scratch memory                                          */
/*             {Variables needed for mem allocation: }                      */
/*                                                                          */

typedef struct cpscr_loc{
  double *vextr,*vexti;        /* Lst: External potential on large 
                                       sphere cutoff g-space grid;  
                                  Lth: ncoef_l large sparse grid       */
  double *vextr_dens_cp_box,*vexti_dens_cp_box;
                               /* Lst: External potential on large 
                                       sphere cutoff g-space grid;  
                                  Lth: ncoef_l_dens_cp_box         */
  double *dvextr,*dvexti;      /* Lst: G-deriv of external potential on large 
                                       sphere cutoff g-space grid;  
                                  Lth: ncoef_l                        */
}CPSCR_LOC;

typedef struct cpscr_nonloc{
 /* All done on the small dense grid  which is dens_cp_box */
  int natm_nls_max;            /* Num: Max # of atms in nl-pseudopot  */
  int nlscr_up;                /* Num: # non-local pseudopot interacts 
                                       for up states =  
                                       nstate_up_proc*(n_ang_max+1)^2*
                                       natm_nls_max                   */
  int nlscr_dn;                /* Num: # non-local e-atm pseudopot
                                       interactions for dn states = 
                                       nstate_dn_proc*(n_ang_max+1)^2*
                                       natm_nls_max                   */

  double *vnlre_up,*vnlim_up;  /* Lst: List of non-local pseudopot
                                      interacts of up states;      
                                  Lth: nlscr_up_proc                       */
  double *vnlre_dn,*vnlim_dn;  /* Lst: List of nl pseudopot 
                                       interacts for dn states;    
                                  Lth: nlscr_dn_proc                       */
  double *dvnlre_x_up,*dvnlre_y_up,*dvnlre_z_up;
  double *dvnlim_x_up,*dvnlim_y_up,*dvnlim_z_up;
                               /* Lst: List of div nl-pseudopot
                                       interacts for up states; 
                                  Lth: nlscr_up_proc                       */
  double *dvnlre_x_dn,*dvnlre_y_dn,*dvnlre_z_dn;
  double *dvnlim_x_dn,*dvnlim_y_dn,*dvnlim_z_dn;
                               /* Lst: List of div nl-pseudopot
                                       interacts for dn states; 
                                  Lth: nlscr_up_proc                       */
  double *dvnlre_gxgx_up,*dvnlim_gxgx_up;
  double *dvnlre_gygy_up,*dvnlim_gygy_up;
  double *dvnlre_gzgz_up,*dvnlim_gzgz_up;
  double *dvnlre_gxgy_up,*dvnlim_gxgy_up;
  double *dvnlre_gygz_up,*dvnlim_gygz_up;
  double *dvnlre_gxgz_up,*dvnlim_gxgz_up;
                               /* Lst: G-derivatives of nonlocal matrix
                                  Lth: nlscr_up_proc                       */
  double *dvnlre_gxgx_dn,*dvnlim_gxgx_dn;
  double *dvnlre_gygy_dn,*dvnlim_gygy_dn;
  double *dvnlre_gzgz_dn,*dvnlim_gzgz_dn;
  double *dvnlre_gxgy_dn,*dvnlim_gxgy_dn;
  double *dvnlre_gygz_dn,*dvnlim_gygz_dn;
  double *dvnlre_gxgz_dn,*dvnlim_gxgz_dn;
                               /* Lst: G-derivatives of nonlocal matrix
                                  Lth: nlscr_dn_proc                       */
}CPSCR_NONLOC;

typedef struct cpscr_rho{
  int iset_map_flag;           /* Opt: make map for dual gridding each step*/
                               /*      0 set once 1 set each step */
  int iset_map_upack_flag;     /* Opt: make map for dual gridding each step*/
                               /*      0 set once 1 set each step */
  int iset_map_flag_pme;       /* Opt: make map for pme dual gridding each step*/
                               /*      0 set once 1 set each step */
  int iset_map_upack_flag_pme; /* Opt: make map for pme dual gridding each step*/
                               /*      0 set once 1 set each step */

  int *map_dual;             /* Lst: array of indices for mapping real space */
                             /*  density from small grid to large grid in 
                                 parallel length nfft_proc
                                 eventually (nfft_proc/nkf1_cp_box)*/
  int *map_dual_upack;

  double *v_ks_up;             /* Lst: KS potential for up density 
                                       on dense cp box real-space grid;
                                  Lth: nfft                           */
  double *v_ks_dn;             /* Lst: KS potential for dn density 
                                       on dense cp box real-space grid;  
                                  Lth: nfft                           */

  double *v_ks_tau_up;         /* Lst: KS potential for up elec KE density 
                                       on dense cp box real-space grid;
                                  Lth: nfft                           */
  double *v_ks_tau_dn;         /* Lst: KS potential for dn elec KE density 
                                       on dense cp box real-space grid;  
                                  Lth: nfft                           */

  double *rho_up;              /* Lst: Spin up density on
                                       dense cp box real-space grid;
                                  Lth: nfft2                          */
  double *rho_dn;              /* Lst: Spin down density on
                                       dense cp box  real-space grid;
                                  Lth: nfft2                          */
  double *rhocr_up,*rhoci_up;  /* Lst: Spin up density on sphere 
                                       cutoff g-space grid;  
                                   large sparse grid
                                  Lth: ncoef_l_up                     */

  double *rhocr_scr,*rhoci_scr;/* Lst: Spin up scratch density on sphere 
                                       cutoff g-space grid;  
                                   large sparse grid
                                  Lth: ncoef_l_up                     */

  double *rhocr_up_dens_cp_box,*rhoci_up_dens_cp_box;
                                /* Lst: Spin up density on sphere 
                                       cutoff g-space grid;  
                                  dense cp box grid
                                  Lth: ncoef_dens_cp_box               */

  double *rhocr_dn,*rhoci_dn;  /* Lst: Spin dn density on sphere 
                                       cutoff g-space grid;    
                                   large sparse grid
                                  Lth: ncoef_l_dn        */
  double *rhocr_dn_dens_cp_box,*rhoci_dn_dens_cp_box;  /*PROBABLY DO NOT NEED*/
                                /* Lst: Spin up density on sphere 
                                       cutoff g-space grid;  
                                     dense cp box grid
                                  Lth: ncoef_dens_cp_box                 */

}CPSCR_RHO;

typedef struct cpscr_grho{
 /* All done on the small dense grid */
  double *d2_rho_up;           /* Lst: Del spin up density on
                                       large square real-space grid;
                                  Lth: nfft2                          */
  double *d2_rho_dn;           /* Lst: Del spin up density on
                                       large square real-space grid;
                                  Lth: nfft2                          */
  double *d2_rho_up_store;     /* Lst: Extra storage space for laplacian
                                       of density when laplacian AND ptens are on
                                  Lth: nfft2 */
  double *d2_rho_dn_store;     /* Lst: Extra storage space for laplacian
                                       of density when laplacian AND ptens are on
                                  Lth: nfft2 */
  double *d_rhox_up,*d_rhoy_up,*d_rhoz_up; 
                               /* Lst: Div spin up density on
                                       large square real-space grid;
                                  Lth: nfft2                          */
  double *d_rhox_dn,*d_rhoy_dn,*d_rhoz_dn; 
                               /* Lst: Div spin dn density on
                                       large square real-space grid;
                                  Lth: nfft2                          */
  double *dm_rhox_up,*dm_rhoy_up,*dm_rhoz_up; 
                               /* Lst: Div|Div| spin up density
                                       on large square real-space grid;
                                  Lth: nfft2                          */
  double *dm_rhox_dn,*dm_rhoy_dn,*dm_rhoz_dn; 
                               /* Lst: Div|Div| spin up density 
                                      on large square real-space grid;
                                  Lth: nfft2                          */
  double *elec_ke_dens_up,*elec_ke_dens_dn;
                                /* Lst: Spin up and Spin down electron
                                        kinetic energy densities
                                   Lth: nfft2                          */
  double *g_rhor_x,*g_rhor_y,*g_rhor_z;
                               /* Lst: Real part of FT of gradient 
                                       of spin density (up or down)
                                       cutoff g-space grid;    
                                  Lth: ncoef_l       */
  double *g_rhoi_x,*g_rhoi_y,*g_rhoi_z;
                               /* Lst: Real part FT of gradient 
                                       of spin density (up or down)
                                       cutoff g-space grid;    
                                  Lth: ncoef_l        */

  double *g2_rhor,*g2_rhoi;    /* Lst: Real/imag part of FT of Laplacian 
                                       of spin density (up or down)
                                       cutoff g-space grid;    
                                  Lth: ncoef_l        */
  double *rhocr_up_st,*rhoci_up_st;/* Lst: Real/imag part density stored only
                                       if laplacian and ptens_calc are on
                                  Lth: ncoef_l        */
  double *rhocr_dn_st,*rhoci_dn_st;/* Lst: Real/imag part density stored only
                                       if laplacian and ptens_calc are on
                                  Lth: ncoef_l        */
} CPSCR_GRHO;

typedef struct cpscr_ovmat{

  double *ovlap1,*ovlap2;
  double *ovlap3,*ovlap4;
  double *ovlap5,*ovlap6;
  double *ovlap7,*ovlap8;      /* Lst: Overlap matrices; Lth: nstate2 */
  double *state_vec1;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec2;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec3;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec4;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec5;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec6;          /* Lst: Eigenvalue vector;Lth: nstate  */
  double *state_vec7;          /* Lst: Eigenvalue vector;Lth: nstate  */

} CPSCR_OVMAT;

typedef struct cpscr_wave{
  /* CP dual grid opt = 0 : standard cp 1 grid */
  /* CP dual grid opt = 1 : dual proportional grids malloced to large
                            dense grid */
  /* CP dual grid opt = 0 : dual not proportional malloced to MAX 
                           of dense cp grid and sparse large grid */
  double *zfft,*zfft_tmp;      /* Lst: PW coeff: square g-space grid;
                                  Lth: nfft                           */
  int zfft_mall_size;     /* Length of zfft  and zfft_tmp*/
 /* Invariant to cp grid opt */
  double *cre_up,*cim_up;      /* Lst: Re/Im part PW coef;
                                  Lth: ncoef*nstate_up_proc            */
  double *cre_dn,*cim_dn;      /* Lst: Re/Im part PW coef;
                                  Lth: ncoef*nstate_up_proc            */
} CPSCR_WAVE;

typedef struct cpscr_therm{
  double *coef_kin;            /* Lst: PW coef KE of each NHC
                                  Lth: num_c_nhc+1                    */
  double *sc_cp;               /* Lst: Scaling of PW NHCs;            
                                  Lth: num_c_nhc+1                    */
} CPSCR_THERM;  

typedef struct cpscr_dual_pme{

 int *iatemp,*ibtemp,*ictemp;       /*Lth: nkf1(2,3)_dens_cp_box   */
 int **igrid_a,**igrid_b,**igrid_c; /*Lth: n_interp_dual_pme 
                                         X nkf1(2,3)_dens_cp_box*/

 int **igrid_now;                   /*Lth: n_interp_dual_pme X nkf1_cp_box */

 double *a_pme,*b_pme,*c_pme;      /*arrays to hold scaled grid points*/
                                   /*Lth: nkf1(2,3)_dens_cp_box       */
 double *frac_a,*frac_b,*frac_c;   /*Lth: nkf1_cp_box*/
 double *aj,*rn,*rn1;              /*Lth: n_interp_pme_dual*/

 double **mn_a,**mn_b,**mn_c;      /*Lth: n_interp_pme_dual X nkf1(2,3)_cp_box*/
 double **ua,**ub,**uc;            /*Lth: n_interp_pme_dual X nkf1(2,3)_cp_box*/
 double *bw_r,*bw_i;               /*Lth: nfft  weighting factors             */
                                   /*  for r->g space pme                     */
/* variables for parallel pme maps */
 int nfft_send_big_small;
 int nfft_recv_big_small;
 int num_rows_tot;
 int *ioff_lg_sm;
 int **joff_sm_lg;

} CPSCR_DUAL_PME;

typedef struct cpscr_atom_pme{

/* NORMAL */
  int pme_on;                       /*Opt: PME on                         */
  int n_interp;                     /*Num: Order of interpolation         */
  int pme_para_opt;
  int nlen_pme;                     /*Num: Scr lngth                      */

  int nktot_pme;                    /*Num: equal to ewald->nktot          */
  int nkf1,nkf2,nkf3;               /*Num: PME mesh same as large sparse grid */
  int ngrid_a,ngrid_b,ngrid_c;      /*Num: PME mesh same as large sparse grid */  

/* SCRATCH */

  int *iatemp,*ibtemp,*ictemp;      /*Lst: Lth: nlen_pme                  */
  int *nc,*ioff_c;                  /*Lst: Lth nkf3                       */

  int **igrid_a,**igrid_b,**igrid_c;/*Lst: Lth: ninterp*nlen_pme          */
  int **igrid_now;                  /*Lst: Lth: ninterp*nlen_pme          */

  double *frac_a,*frac_b,*frac_c;   /*Lst: Lth:nlen_pme                   */
  double *aj,*rn,*rn1;              /*Lst: Lth: ninterp*nlen_pme          */

  double **ua,**ub,**uc;            /*Lst: Lth:ninterp*nlen_pme           */
  double **mn_a,**mn_b,**mn_c;      /*Lst: Lth: ninterp*nlen_pme          */
  double **dmn_a,**dmn_b,**dmn_c;   /*Lst: Lth: ninterp*nlen_pme          */
  double **qgrid_now;               /*Lst: Lth: ninterp*nlen_pme          */
  double *qgrid;                    /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_scr;                /*Lst: Interpolation grid 
                                            Lth:2*ngrid_a*_b*_c           */
  double *qgrid_tmp_real;           /*Lst: Lth: nktot                     */
  double *qgrid_tmp_imag;           /*Lst: Lth: nktot                     */
  double *bweight_tot;              /*Lst: Lth: nktot : May not need?     */
  double *bw_r,*bw_i;               /*Lth: nfft  weighting factors        */
                                    /*  for r->g space pme                */

} CPSCR_ATOM_PME;

typedef struct cpscr{
  CPSCR_LOC    cpscr_loc;
  CPSCR_NONLOC cpscr_nonloc;
  CPSCR_RHO    cpscr_rho;
  CPSCR_GRHO   cpscr_grho;
  CPSCR_OVMAT  cpscr_ovmat;
  CPSCR_WAVE   cpscr_wave;
  CPSCR_THERM  cpscr_therm;  
  CPSCR_DUAL_PME cpscr_dual_pme;
  CPSCR_ATOM_PME cpscr_atom_pme;
}CPSCR;

/*==========================================================================*/
/*                 Velocity Resampling info                                 */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */

typedef struct vel_samp_cp{
  int ivelc_smpl_on;           /* Opt: Periodic PW coef vel resampl   */ 
  int ivelc_scal_on;           /* Opt: Periodic PW coef vel rescale   */ 
  int nvc_smpl;                /* Num: Freq of PW coef vel resampl    */
  int nvcnhc_smpl;             /* Num: Freq of PW coef NHC vel resamp */
  int nvc_scal;                /* Num: Freq of PW coef vel recale     */
  int iseed,iseed2;            /* Num: Random seeds                   */
  int iauto_vc_scal_opt;       /* Opt: auto rescale options on/off    */
  double qseed;                /* Num: Real seed for essl ran()       */
  double vc_scal_tol;          /* Num: tol of auto rscale of coef vel */
  double div_scal;
} VEL_SAMP_CP;

/*==========================================================================*/
/*                 Spherical Harmonic constants                             */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */

typedef struct ylm_cons {
  double rt_fpi,rt_thrfpi,rt_threpi;
  double hrt_fivfpi,rt_fiftepi,hrt_sevfpi;
  double hrt_toepi,hrt_ohffpi,hrt_tfepi;
} YLM_CONS;


/*==========================================================================*/
/*==========================================================================*/
/*                       CP                                                 */

typedef struct cp {
  CPOPTS cpopts;
  CPCOEFFS_INFO cpcoeffs_info;
  CPCOEFFS_POS *cpcoeffs_pos;
  CPTHERM_INFO cptherm_info;
  CPTHERM_POS *cptherm_pos;
  ELECTRONIC_PROPERTIES electronic_properties;
  CPEWALD cpewald; 
  CPCONSTRNT cpconstrnt;
  PSEUDO pseudo;
  CPSCR cpscr;
  COMMUNICATE communicate;
  PARA_FFT_PKG3D cp_sclr_fft_pkg3d_sm;
  PARA_FFT_PKG3D cp_para_fft_pkg3d_sm;
  PARA_FFT_PKG3D cp_sclr_fft_pkg3d_dens_cp_box;
  PARA_FFT_PKG3D cp_para_fft_pkg3d_dens_cp_box;
  PARA_FFT_PKG3D cp_sclr_fft_pkg3d_lg;
  PARA_FFT_PKG3D cp_para_fft_pkg3d_lg;
  CP_COMM_STATE_PKG cp_comm_state_pkg_up;
  CP_COMM_STATE_PKG cp_comm_state_pkg_dn;
  VEL_SAMP_CP vel_samp_cp;
} CP;



