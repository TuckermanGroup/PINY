 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/
 /*                                                                          */
 /*                         PI_MD:                                           */
 /*             The future of simulation technology                          */
 /*             ------------------------------------                         */
 /*                   Structures: typedefs_opt.h                             */
 /*                                                                          */
 /* The include file with the typedefs of all the option structures          */
 /*                                                                          */
 /*==========================================================================*/
 /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
 /*==========================================================================*/


/*==========================================================================*/
/*                  State point data                                        */
/*             {Variables needed for mem allocation:}                       */

typedef struct statepoint{
  double pext,t_ext, stens_ext;   /* Num: External press-temp, surface tens */ 
} STATEPOINT;

/*==========================================================================*/
/*                  Simulation cell data                                    */

typedef struct cell{
  int iperd;                   /* Opt: Periodicity of system(0-3)     */
  int intra_perds;             /* Opt: Intramol periodic flag         */
  int cubic_box_flag;          /* Opt: Cubic box flag for fast imaging */
  int hmat_int_typ;            /* Opt: Flag to distinguish between normal and
                                        upper triangular integration */
  int hmat_cons_typ;           /* Opt: Constraint on the box/h_mat 
                                       0=none,1=orthorhombic,2=monoclinic*/
  int imov_cp_box;             /* Opt: Flag to move center of cp_box */

  double vol;                  /* Num: The volume                     */
  double vol_cp;               /* Num: The volume                     */
  double vol0;                 /* Num: The initial volume            */
  double area;                 /* Num: The area                      */
  double alpha_conv_dual;

  double *hmat,*hmati;             /* Lst: Cell shape; Lth:9              */
  double *hmat_ewd,*hmat_ewd_cp;   /* Lst: Cell shape for g-vectors Lth:9 */
  double *hmat_cp,*hmati_cp;       /* Lst: Cell shape; Lth:9              */
  double *cp_box_center;           
  double *cp_box_center_rel;
  double *cp_vbox_center;          /*box center velocities */
  double *cp_fbox_center;          /*box center forces     */
} CELL;

/*==========================================================================*/
/*                  Timestep info                                           */
/*             {Variables needed for mem allocation:}                       */

typedef struct timeinfo{
  int ntime,itime;             /* Num: Time-minimization steps           */
  int ix_respa;                /* Opt: Extended system respa opt         */
  int int_res_ter,int_res_tra; /* Opt: Intra-inter respa opts            */
  int nres_ter,nres_tra;       /* Num: # intra-inter respa steps         */ 
  int int_res_tor,nres_tor;    /* Opt: Torsion respa opts                */
  int int_res_pimd,nres_pimd;  /* Opt: Path integral respa options       */
  int iget_pe_real_inter_freq;/*Opt: Freq calc of realspace inter PE     */
  int exit_flag;               /* Num: Flat used to tell code to exit.   */
                               /*      Used with annealing and auto exit */
                               /*      options.                          */
  double dt;                   /* Num: Time step                         */
} TIMEINFO;


/*==========================================================================*/
/*                  Andersen-Hoover NPT (NPT_I) info                        */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */
typedef struct baro {
  int len_nhc;                 /* Num: length of NHC                  */
  double x_lnv,v_lnv;          /* Num: log(vol),dlog(vol)/dt          */
  double v_lnv_glob;           /* Num: dlog(vol)/dt                   */
  double v_lnv_g;              /* Num: dlog(vol)/dt                   */
  double f_lnv_p,f_lnv_v;      /* Num: d^2log(vol)/dt^2               */
  double vol;                  /* Num: volume                         */
  double mass_lnv;             /* Num: Mass of log(vol)               */
  double c2_lnv;               /* Num: Useful constant                */
  double roll_scv;             /* Num: Roll scaling factor            */
  double roll_scg;             /* Num: Roll scaling factor            */
  double roll_scg0;            /* Num: Roll scaling factor            */
  double x_lnv_o;              /* Num: old log(vol)                   */
  double v_lnv_g_wght;         /* Num: weighted v_lnv_g               */
  double area;                 /* Num: area                           */

  double *x_vol_nhc,*v_vol_nhc;/* Lst: Volume NHCs Lth:len_nhc        */
  double *f_vol_nhc,*mass_vol_nhc,*gkt_vol;
  double *hmato;               /* Num: old cell matrix  :Lth 9        */
} BARO;

/*==========================================================================*/
/*                  Nuevo-Parrinello-Rahman info                            */
/*             {Variables needed for mem allocation:                        */
/*                                                  }                       */
/*                                                                          */
typedef struct par_rahman{
  double vol;
  double mass_hm;              /* Num: Mass of hh^-1 matrix           */
  double c1_hm;                /* Num: Useful constant                */
  double roll_scg;             /* Num: Roll scaling factor            */
  double roll_scg0;            /* Num: Roll scaling factor            */
  double area;                 /* Num: area                           */

  double *vgmat,*vgmat_g;      /* Lst: Velocity of hh^-1 matrix Lth:9 */
  double *vgmat_glob;          /* Lst: Velocity of hh^-1 matrix Lth:9 */
  double *fgmat_p,*fgmat_v;    /* Lst: Force on hh^-1 matrix    Lth:9 */
  double *roll_mtv,*roll_mtx;  /* Lst: Roll matrices            Lth:9 */
  double *roll_mtvv;           /* Lst: Roll matrices            Lth:9 */
  double *vtemps,*veigv,*veig; /* Lst: Integration temp         Lth:9 */
  double *vexpdt,*vsindt;      /* Lst: Integration temp         Lth:9 */
  double *vtempx, *vtempv;     /* Lst: Integration temp         Lth:9 */
  double *hmat_t, *hmato;      /* Lst: cell shape temps         Lth:9 */
  double *fv1,*fv2;            /* Lst: Rs temps                 Lth:3 */ 
  double *vgmat_g_wght;        /* Lst: weighted vgmat_g         Lth:9 */
} PAR_RAHMAN;

/*==========================================================================*/
/*               Pressure and Kinetic tensors                               */

typedef struct ptens {
  double *tvten;               /* Lst: KE  tensor           ; Lth: 9  */
  double *pvten,*pvten_tot;    /* Lst: PV  tensors          ; Lth: 9  */
  double *pvten_inc,*pvten_tmp;/* Lst: PV  tensors for SHAKE; Lth: 9  */
                               /* Lst: PV  tensors for RESPA; Lth: 9  */
  double *pvten_tmp_res       ;/* Lst: PV  tensors for RESPA; Lth: 9  */
  double *pvten_inc_t1,*pvten_inc_t2,*pvten_inc_t3,*pvten_inc_t4;
  double *pvten_inc_a;       
  double *pvten_inc_old;
  double *count_inc;           /* Lst: pvten_inc_ts counters; Lth: 4  */
  double *pvten_inc_glob,*pvten_inc_whol,*pvten_inc_std;
  double *pvten_mode,*pvten_mode_tot;
} PTENS;

/*=======================================================================*/
/*               Static averages                                         */
/*             {Variables needed for mem allocation:                     */
/*                                                                       */

typedef struct stat_avg {
  int nexclt;                        /* Num: Tot # of exclusions       */
  int iter_shake, iter_ratl; 
  int iswit_vdw;                     /* Opt: Switch inter shift on/off */
  int itime_update,itime_update_w;   /* Num: TIme of updates           */
  int iter_shake_cp,iter_ratl_cp;
  int write_cp_atm_flag;

  double vintert,aivintert,avintert; /* Num:Inst, and avgs inter PE    */
  double vintrat,aivintrat,avintrat; /* Num:Inst, and avgs intra PE    */
  double vsurft;                     /* Num: surface energy pot        */
  double vbondt,vbendt,vtorst,vonfot;/* Num: Bond,bend,tors, and onfo  */
  double vbend_bndt;                 /* NUM: Uri Bradleys              */
  double vbend_bnd_bond,vbend_bnd_bend;/* NUM: Uri Bradleys decomp   */
  double vbondt_watts,vbendt_watts,vtot_watts;/* NUM: Watts dcomp  */
  double vreal,vrecip;               /* Num: Inter mol PE              */
  double vvdw,vcoul;                 /* Num: Van der Waals and Coulomb */
  double vlong;                      /* Num: Long range correction to LJ*/
  double vbond_free,vbend_free,vtors_free; /* Num: Free intra pots     */
  double vbar_free;                  /* Num: more free energy pots     */
  double kinet,aikinet,akinet;       /* Num:Inst, and avg atm KE       */
  double kinet_v,aikinet_v,akinet_v; /* Num:Inst, and avg volume KE    */
  double vol,aivol,avol;             /* Num:Inst, and avg volume       */
  double kinet_nhc,aikinet_nhc,akinet_nhc;/*Num:Inst and avg Atm NHC KE*/
  double kinet_nhc_bead;
  double aikinet_nhc_bead;
  double akinet_nhc_bead;/*Num:Inst and avg Atm NHC KE*/
  double vpot_v;                     /* Num: Volume PE                 */
  double vpotnhc;                    /* Num: NHC PE                    */
  double aiter_shake,aiter_ratl;     /* Num:Inst and avg # Shake/Ratl  */ 
  double iter_23, iter_33, iter_46, iter_43, iter_21; 
  double aiter_23,aiter_33, aiter_46; /* Num:Inst and avg # grp Shake/Ratl  */ 
  double aiter_21,aiter_43;           /* Num:Inst and avg # grp Shake/Ratl  */ 
  double iter_23r, iter_33r, iter_46r, iter_43r, iter_21r; 
  double aiter_23r,aiter_33r, aiter_46r;/* Num:Inst and avg # grp Shake/Ratl*/
  double aiter_21r,aiter_43r;/* Num:Inst and avg # grp Shake/Ratl*/
  double acella,acellb,acellc;
  double aicella,aicellb,aicellc;    /* Num:Inst and avg cell lngth    */
  double acellab,acellbc,acellac;
  double aicellab,aicellbc,aicellac; /* Num:Inst and avg cell angles   */
  double apress,aipress;             /* Num: Avg, inst avg, pressure   */
  double press_inter,press_intra;    /* Num: Inter, intra pressure     */
  double apress_inter,aipress_inter; /* Num: Avg, inst avg, inter pressure   */
  double apress_intra,aipress_intra; /* Num: Avg, inst avg, intra pressure   */
  double press_kin;
  double apress_kin,aipress_kin; /* Num: Avg, inst avg, intra pressure   */
  double econv0,econv;               /* Num: Energy conservation       */
  double cp_kconv0,cp_kconv;         /* Num: CP kinetic energy conservation (isokinetic) */
  double cpu1,cpu2,acpu,cpu_now;     /* Num: Cpu time                  */
  double updates,updates_w;          /* Num: Ver-list,NL-list updates  */
  double kinet_cp_up,aikinet_cp_up,akinet_cp_up;  /* Num: Up CP fict KE */
  double kinet_cp_dn,aikinet_cp_dn,akinet_cp_dn;  /* Num: Up CP fict KE */
  double kinet_cp,aikinet_cp,akinet_cp;/* Num:Inst and avg CP-Class KE */
  double kinet_nhc_cp,aikinet_nhc_cp,akinet_nhc_cp;/* Num: CP-NHC KE   */
  double vpotnhc_cp;                 /* Num:Inst CP-NHC PE             */
  double cp_ehart,aicp_ehart,acp_ehart; /* Num:Inst, CP Hartree E      */
  double cp_eext,aicp_eext,acp_eext; /* Num:Inst  and avg CP Vext      */
  double cp_exc,aicp_exc,acp_exc;    /* Num:Inst and avg CP exc  E     */
  double cp_muxc;                    /* Num:Integral of xc pot times rho */
  double cp_eke,aicp_eke,acp_eke;    /* Num:Inst and avg CP KE         */
  double cp_enl,aicp_enl,acp_enl;    /* Num:Inst and avg CP V_nonlocal */
  double aiter_shake_cp,aiter_ratl_cp;/* Num: Inst and avg # Shake/Ratl 
                                              iterations for CP        */ 
  double maxfc,maxf;                 /* Num: MAX component of the force
                                             maxfc=coefs,maxf=atms     */ 
  double pi_ke_prim,pi_ke_vir;       /* Num: Quantum KE estimators     */
  double api_ke_prim,api_ke_vir;     /* Num: Average quantum KE est.   */
  double aipi_ke_prim,aipi_ke_vir;   /* Num: Inst. Avg. quantum KE est.*/
  double kin_harm,akin_harm,aikin_harm; /* Num: Harmonic KE             */

  double *apten, *aipten, *apten_out;/* Lst: Avg pressure tensor; Lth:9*/
  double fc_mag_up,fc_mag_dn;        /* Num: avg mag of coef force     */
  double fc_max_up,fc_max_dn;        /* Num: Max mag of coef force     */
  double fatm_max,fatm_mag;          /* Num: Max mag of atm  force     */
  double count_diag_srot;            /* Num: Number of rotations to 
                                             diagonal ovlap basis      */
  double econv_now;                       /* Num: The conserved quantity    */
  
} STAT_AVG;

/*==========================================================================*/
/*                 Indexing info for CP coeffs and Ewald                    */
/*             {Variables needed for mem allocation:                        */
/*                   nktot_sm,nktot,nktot_res}                              */
/*                                                                          */
typedef struct ewald {
  int nsplin_g;
  int nktot;                   /* Num: # of PW coeff on large sphere
                                      cutoff g-space grid =ncoef_l    */
  int nktot_res;               /* Num: # of PW coeffon RESPA sphere
                                     cutoff g-space grid <=ncoef_l    */

  int nktot_mall;      
  int nktot_res_mall;  
  int nkc_max;                 /* Num: Max value of kc, calculated by all */
  int nkb_max;                 /* Num: Max value of kb, calculated by all */

  double alp_ewd;              /* Num: Convergence param of Ewald sum */ 
  double self_erf;            
  double ecut,ecut_res;
  double ecut_clus;            /* Num: Energy cutoff for cluster BCs  */
  double alp_clus;             /* Num: Ewald alpha for cluster BCs */
  double gs_max;               /* Num: Max (gx^2+gy^2)^1/2 calc by all */

  int *kastr,*kbstr,*kcstr;    /* Lst: Large spherically cutoff g vecs
                                 Lth: nktot                           */ 
  int *ibrk1,*ibrk2,*ibrk3;    /* Lst: Large spherically cutoff 
                                       g-space grid break-pts;
                                  Lth: nktot                          */ 
  int *kastr_res,*kbstr_res,*kcstr_res; /* Lst: RESPA sphere cutoff 
                                                g-vectors; 
                                           Lth: nktot_res             */ 
  int *ibrk1_res,*ibrk2_res;   /* Lst: RESPA spherically cutoff 
                                       g-space grid break-pts;  
                                  Lth: nktot_res                      */
  double *clus_corr_r,*dclus_corr_r;
  double *clus_corr_r_res,*dclus_corr_r_res;
  double *cvgs_2d0,*cvgs_2d1,*cvgs_2d2,*cvgs_2d3,*vgc_2d;
  double *vgc_1d_y,*vgc_1d_z;
} EWALD;

/*==========================================================================*/
/*                  Simulation options                                      */
/*             {Variables needed for mem allocation:}                       */

typedef struct simopts {
  int md;                      /* Opt: Classical MD                         */
  int pimd;                    /* Opt: Path integral MD                     */
  int minimize;                /* Opt: Classical minimization               */
  int cp;                      /* Opt: Full CP                              */
  int cp_pimd;                 /* Opt: Full CP                              */
  int cp_wave;                 /* Opt: CP, wave function only               */
  int cp_wave_pimd;            /* Opt: CP, wave function only               */
  int cp_wave_min_pimd;        /* Opt: CP min, wave function only           */
  int cp_min;                  /* Opt: Full CP minimization                 */
  int cp_wave_min;             /* Opt: CP min, wave function only           */
  int debug;                   /* Opt: Internal use-backdoor chks           */
  int debug_pimd;              /* Opt: Internal use-backdoor chks for pimd  */
  int debug_cp;                /* Opt: Internal use-backdoor chks for CP    */
  int debug_cp_pimd;           /* Opt: Internal use-backdoor chks for CP    */
  int pi_beads;                /* Num: # of path integral descritizations   */
  int pi_md_typ;               /* Opt: Staging or normal modes              */
  int initial_spread_opt;      /* Opt: Spread coordinates for pimd          */
  int anneal_opt;              /* Opt: Do simulated annealing (on/off)      */
  double ann_rate;             /* Num: Annealing rate                       */
  double ann_start_temp;       /* Num: Annealing start temperature          */
  double ann_target_temp;      /* Num: Annealing final temperature          */
} SIMOPTS;

/*==========================================================================*/
/*                  Minimization options                                    */
/*             {Variables needed for mem allocation:}                       */

typedef struct minopts {
  int min_std;                 /* Opt: Steepest descent min           */
  int min_cg;                  /* Opt: Conjugate gradiant min         */
  int min_diis;                /* Opt: DIIS min                       */
  int diis_hist_len;           /* Num: Length of DIIS history         */
  int cp_min_std;              /* Opt: Steepest descent min           */
  int cp_min_cg;               /* Opt: Conjugate gradiant min         */
  int cp_min_diis;             /* Opt: DIIS min                       */
  int cp_diis_hist_len;        /* Num: Length of DIIS history         */
  int cp_cg_line_min_len;       /* Num: COMMENT ME */
  int min_atm_com_fix_opt;     /* Opt: Keep the com fixed             */

  double tol_coef;             /* Num: Tol on PW(plane wave) coef forces   */
  double tol_atom;             /* Num: Tol on atm forces              */
} MINOPTS;

/*==========================================================================*/
/*                  Ensemble options                                        */
/*             {Variables needed for mem allocation:}                       */

typedef struct ensopts{
  int nve,nvt,npt_i,npt_f,nst; /* Opt: Stat mech ensembles            */
} ENSOPTS;

/*==========================================================================*/
/*                 Input-output                                             */
/*             {Variables needed for mem allocation:                        */
/*                       file_len                           }               */
/*                                                                          */
typedef struct filenames {
  char *iname;                 /* Chr: Instananeous data file name    */
  char *dname;                 /* Chr: Atm dump file name             */
  char *dnamec;                /* Chr: PW coef dump file name         */
  char *cpname;                /* Chr: Atm pos conf file name         */
  char *ksname;                /* Chr: KS eigs file name              */
  char *elfname;               /* Chr: ELF file name                  */
  char *cvname;                /* Chr: Atm vel conf file name         */
  char *ccname;                /* Chr: PW coef conf file name         */
  char *cpparname;             /* Chr: Atm partial pos conf file name */
  char *centname;              /* Chr: Centroid conf file name        */
  char *forcename;             /* Chr: Force conf file name           */

  int low_lim_par;             /* Num: lower limit of partial conf write */
  int high_lim_par;            /* Num: upper limit of partial conf write */
  int iwrite_screen;           /* Num: Freq of screen writes          */
  int iwrite_inst;             /* Num: Freq of atm-pos config writes  */
  int iwrite_dump;             /* Num: Freq of atm and PW 
                                       coef dump file writes          */
  int iwrite_kseigs;           /* Num: Freq of KS eigs writes         */
  int iwrite_confv;            /* Num: Freq of atm-vel conf writes    */
  int iwrite_confp;            /* Num: Freq of atm-pos conf writes    */
  int iwrite_par_confp;        /* Num: Freq of atm-pos partial conf writes */
  int iwrite_confc;            /* Num: Freq of PW coef conf writes    */
  int iwrite_path_cent;        /* Num: Freq of centroid conf writes   */
  int ifile_open;              /* file open flag                      */
  int iwrite_conf_binary;      /* Opt: Write conf files in binary     */
  int iwrite_units;            /* Opt: Write screen output units      */
  int iwrite_atm_for;          /* Num: Freq of atm force              */
  int iwrite_elf; 
} FILENAMES;


/*==========================================================================*/
/*                  FFT Package variables                                   */
/*                                                                          */

typedef struct para_fft_pkg3d {

  /* General info */
  /*--------------*/
   int nktot,ncoef,ncoef_proc;          /* Num: Sphere cutoff coefs        */
   int ncoef_use,icoef_off,icoef_strt;  /* Num: Sphere cutoff coefs offset */
   int nfft_size;                       /* Num: Size of FFT grid required  */
   int nfft,nfft_proc;                  /* Num: Size real space and per proc*/
   int igeneric_opt;                    /* Num: Generic FFT option         */
   int scale_opt;                       /* Num: Back transform scale option*/
   int len_ifax;                        /* Num: Size of int scratch arrays */
   int nwork1,nwork2;                   /* Num: Size of FFT work arrays    */
   double ecut;                         /* Num: Energy cutoffs             */

  /* Communicator */
  /*--------------*/
   MPI_Comm comm;
   int num_proc,myid,myidp1;

  /* Maps for initial packing */
  /*--------------------------*/
   int *map_proc,*map_c_proc,*map_proc_post; /* Num: Map coefs to FFT order */

  /* Comm info for r-space allgather*/
  /*--------------------------------*/
  int *recv_counts_rho,*displs_rho;

  /* Comm info for g-space reduce_scatter*/
  /*--------------------------------*/
  int *recv_counts_coef,*recv_dspls_coef;

  /* Comm info for pme              */
  /*--------------------------------*/
  int *map_pme_dn,nrecv_dn;
  int *recvcounts_pme_dn,*recvdspls_pme_dn;
  int *sendcounts_pme_dn,*senddspls_pme_dn;
  int *map_pme_up,nrecv_up,nrecv_up_max;
  int *recvcounts_pme_up,*recvdspls_pme_up;
  int *sendcounts_pme_up,*senddspls_pme_up;

  /* FFT along kc  */
  /*--------------*/
   int ndata_kc;
   int nkf3,nfft_kc,nfft_kc_proc;            /* Num: Size and numof kc FFTs */
   int str_fft_kc_proc,end_fft_kc_proc;      /* Num: Str/End of my FFT's    */
   int ska_fft_kc_proc,eka_fft_kc_proc;      /* Num: my ka - str/end index  */
   int skb_fft_kc_proc,ekb_fft_kc_proc;      /* Num: my kb - str/end index  */
   int *sendcounts_fft_kc,*senddspls_fft_kc; /* Lst: Kc communication       */
   int *recvcounts_fft_kc,*recvdspls_fft_kc; /* Lst: Kc communication       */
   int *sum_fft_kc_proc;                     /* Lst: FFt's per ka value     */
   int nka_fft_kc, *nkb_fft_kc;              /* Lst: Ka vals and their kb's */
   int *ka_fft_kc,*kb_fft_kc,*ka_fft_kc_red; /* Lst: ka,kb values           */
   int *ifax_c_f,*ifax_c_r;                  /* Lst: work space             */
   double *work_1c_f,*work_2c_f;
   double *work_1c_r,*work_2c_r;

  /* FFT along kb  */
  /*---------------*/
   int ndata_kb;
   int nkf2, nfft_kb,nfft_kb_proc;           /* Num: Size and num of kb FFTs */
   int str_fft_kb_proc,end_fft_kb_proc;      /* Num: Str/End of my FFT's     */
   int ska_fft_kb_proc,eka_fft_kb_proc;      /* Num: my ka - str/end index   */
   int skc_fft_kb_proc,ekc_fft_kb_proc;      /* Num: my kb - str/end index   */
   int *ska_fft_kb_proc_all,*eka_fft_kb_proc_all;
                                            /* Lst: str/end values each proc */
   int *skc_fft_kb_proc_all,*ekc_fft_kb_proc_all;
   int *sendcounts_fft_kb,*senddspls_fft_kb;  /* Lst: Kb communication       */
   int *recvcounts_fft_kb,*recvdspls_fft_kb;  /* Lst: Kc communication       */
   int *sum_fft_kb_proc;                      /* Lst: FFt's per kc value     */
   int *ifax_b_f,*ifax_b_r;                   /* Lst: work space             */
   double *work_1b_f,*work_2b_f;
   double *work_1b_r,*work_2b_r;

  /* FFT along ka  */
  /*---------------*/
   int ndata_ka;
   int nkf1,nfft_ka,nfft_ka_proc;            /* Num: Size and num of ka FFTs */
   int str_fft_ka_proc,end_fft_ka_proc;      /* Num: Str/End of my FFT's     */
   int skb_fft_ka_proc,ekb_fft_ka_proc;      /* Num: my kb - str/end index   */
   int skc_fft_ka_proc,ekc_fft_ka_proc;      /* Num: my kc - str/end index   */
   int *skb_fft_ka_proc_all,*ekb_fft_ka_proc_all;
                                           /* Lst: str/end values each proc  */
   int *skc_fft_ka_proc_all,*ekc_fft_ka_proc_all;
   int *sendcounts_fft_ka,*senddspls_fft_ka;  /* Lst: Ka communication       */
   int *recvcounts_fft_ka,*recvdspls_fft_ka;  /* Lst: Ka communication       */
   int *ifax_a_f,*ifax_a_r;                   /* Lst: work space             */
   double *work_1a_f,*work_2a_f;
   double *work_1a_r,*work_2a_r;

  /* map variables for the dual grid option */
  /* used in map of the real space density from small grid to the large grid*/
   int *send_counts_row_big_small;          /* Length nproc */
   int *recv_counts_row_big_small;
   int *sdispls_row_big_small;
   int *rdispls_row_big_small;  /*Length nproc used in alltoallv in routine */

   int *send_counts_dual_map;          /* Length nproc */
   int *recv_counts_dual_map;
   int *sdispls_dual_map;
   int *rdispls_dual_map;   /*Length nproc used in alltoallv in routine */

   int *send_counts_row_big_small_upack;          /* Length nproc */
   int *recv_counts_row_big_small_upack;
   int *sdispls_row_big_small_upack;
   int *rdispls_row_big_small_upack;

   int *send_counts_dual_map_upack;          /* Length nproc */
   int *recv_counts_dual_map_upack;
   int *sdispls_dual_map_upack;
   int *rdispls_dual_map_upack;  /*Length nproc used in alltoallv in routine */

   int *send_counts_ioff_big_small;
   int *recv_counts_ioff_big_small;
   int *sdispls_ioff_big_small; /*Length nproc used in pme map parallel routine*/
   int *rdispls_ioff_big_small;


 } PARA_FFT_PKG3D;

/*==========================================================================*/
/*==========================================================================*/
/*                       The options                                        */

typedef struct general_data {
  SIMOPTS simopts;
  MINOPTS minopts;
  ENSOPTS ensopts;
  FILENAMES filenames;
  STATEPOINT statepoint;
  CELL cell;
  TIMEINFO timeinfo;
  BARO baro;
  PAR_RAHMAN par_rahman;
  PTENS ptens;
  STAT_AVG stat_avg;
  EWALD ewald;
  PARA_FFT_PKG3D pme_fft_pkg;
  PARA_FFT_PKG3D pme_res_fft_pkg;
  int error_check_on;
} GENERAL_DATA;

