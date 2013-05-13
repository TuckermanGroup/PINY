/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_sim_dict_list.c                          */
/*                                                                          */
/* This subprogram reads in simulation inputs                               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_sim_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_fun(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 24;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

  for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
  for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
  for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up list flags                                                */
/*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/ 
  /*  1)~sim_list_def[ ] */
        strcpy((*dict)[1].error_mes," ");
        strcpy((*dict)[1].keyword,"sim_list_def");
        strcpy((*dict)[1].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  2)~sim_cp_def[] */
        strcpy((*dict)[2].error_mes," ");
        strcpy((*dict)[2].keyword,"sim_cp_def");
        strcpy((*dict)[2].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  3)~sim_gen_def[ ] */
        strcpy((*dict)[3].error_mes," ");
        strcpy((*dict)[3].keyword,"sim_gen_def");
        strcpy((*dict)[3].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  4)~sim_class_PE_def[ ] */
        strcpy((*dict)[4].error_mes," ");
        strcpy((*dict)[4].keyword,"sim_class_PE_def");
        strcpy((*dict)[4].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  5)~sim_run_def[] */
        strcpy((*dict)[5].error_mes," ");
        strcpy((*dict)[5].keyword,"sim_run_def");
        strcpy((*dict)[5].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  6)~sim_nhc_def[ ] */
        strcpy((*dict)[6].error_mes," ");
        strcpy((*dict)[6].keyword,"sim_nhc_def");
        strcpy((*dict)[6].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  7)~sim_vol_def[ ] */
        strcpy((*dict)[7].error_mes," ");
        strcpy((*dict)[7].keyword,"sim_vol_def");
        strcpy((*dict)[7].keyarg," ");
  /*-----------------------------------------------------------------------*/
  /*  8)~sim_write_def[] */
        strcpy((*dict)[8].error_mes," ");
        strcpy((*dict)[8].keyword,"sim_write_def");
        strcpy((*dict)[8].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  9)~sim_pimd_def[ ] */
        strcpy((*dict)[9].error_mes," ");
        strcpy((*dict)[9].keyword,"sim_pimd_def");
        strcpy((*dict)[9].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  10)~sim_velo_corel[ ] */
        strcpy((*dict)[10].error_mes," ");
        strcpy((*dict)[10].keyword,"sim_velo_corel");
        strcpy((*dict)[10].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  11)~sim_msqd_corel[ ] */
        strcpy((*dict)[11].error_mes," ");
        strcpy((*dict)[11].keyword,"sim_msqd_corel");
        strcpy((*dict)[11].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  12)~sim_iikt_iso_corel[ ] */
        strcpy((*dict)[12].error_mes," ");
        strcpy((*dict)[12].keyword,"sim_iikt_iso_corel");
        strcpy((*dict)[12].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  13)~sim_ickt_iso_corel[ ] */
        strcpy((*dict)[13].error_mes," ");
        strcpy((*dict)[13].keyword,"sim_ickt_iso_corel");
        strcpy((*dict)[13].keyarg," ");
  /*-----------------------------------------------------------------------*/ 
  /*  14)~sim_rdf_corel[ ] */
        strcpy((*dict)[14].error_mes," ");
        strcpy((*dict)[14].keyword,"sim_rdf_corel");
        strcpy((*dict)[14].keyarg," ");
  /*--------------------------------------------------------------------*/ 
  /*  15) ~harmonic_analysis[] */
        strcpy((*dict)[15].keyword,"harmonic_analysis");
        strcpy((*dict)[15].keyarg,"");  
        strcpy((*dict)[15].error_mes,"");
  /*-----------------------------------------------------------------------*/
  /*  16) ~molecule_def[] */
        strcpy((*dict)[16].keyword,"molecule_def");
        strcpy((*dict)[16].keyarg,"");  
        strcpy((*dict)[16].error_mes,"");
        (*dict)[16].key_type=2;                    /*specify more than once*/
  /*-----------------------------------------------------------------------*/ 
  /*  17) ~wavefunc_def[] */
        strcpy((*dict)[17].keyword,"wavefunc_def");
        strcpy((*dict)[17].keyarg,"");  
        strcpy((*dict)[17].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  18) ~bond_free_def[] */
        strcpy((*dict)[18].keyword,"bond_free_def");
        strcpy((*dict)[18].keyarg,"");  
        strcpy((*dict)[18].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  19) ~bend_free_def[] */
        strcpy((*dict)[19].keyword,"bend_free_def");
        strcpy((*dict)[19].keyarg,"");  
        strcpy((*dict)[19].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  20) ~tors_free_def[] */
        strcpy((*dict)[20].keyword,"tors_free_def");
        strcpy((*dict)[20].keyarg,"");  
        strcpy((*dict)[20].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  21) ~data_base_def[] */
        strcpy((*dict)[21].keyword,"data_base_def");
        strcpy((*dict)[21].keyarg,"");  
        strcpy((*dict)[21].error_mes,
         "inter_file,vps_file,bond_file,bend_file,tors_file,onfo_file");
  /*-----------------------------------------------------------------------*/ 
  /*  22) ~user_data_base_def[] */
        strcpy((*dict)[22].keyword,"user_data_base_def");
        strcpy((*dict)[22].keyarg,"");  
        strcpy((*dict)[22].error_mes,
      "inter_file,vps_file,bond_file,bend_file,tors_file,onfo_file,surf_file");
  /*-----------------------------------------------------------------------*/ 
  /*  23) ~rbar_sig_free_def[] */
        strcpy((*dict)[23].keyword,"rbar_sig_free_def");
        strcpy((*dict)[23].keyarg,"");  
        strcpy((*dict)[23].error_mes,"");
  /*-----------------------------------------------------------------------*/ 
  /*  24) ~surface_def[] */
        strcpy((*dict)[24].keyword,"surface_def");
        strcpy((*dict)[24].keyarg,"");  
        strcpy((*dict)[24].error_mes,"");
  /*--------------------------------------------------------------------*/ 

/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_list(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 14;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up list flags                                                */
/*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/ 
  /*  1)\verlist_pad{#} */
        strcpy((*dict)[1].error_mes,"a number > 0 ");
        strcpy((*dict)[1].keyword,"verlist_pad");
        strcpy((*dict)[1].keyarg,"30");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\verlist_mem_safe{#} */
        strcpy((*dict)[2].error_mes,"a number > 0 ");
        strcpy((*dict)[2].keyword,"verlist_mem_safe");
        strcpy((*dict)[2].keyarg,"1.25");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\verlist_mem_min{#} */
        strcpy((*dict)[3].error_mes,"a number > 0 ");
        strcpy((*dict)[3].keyword,"verlist_mem_min");
        strcpy((*dict)[3].keyarg,"2000");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\verlist_skin{#} */
        strcpy((*dict)[4].error_mes,"a number > 0 ");
        strcpy((*dict)[4].keyword,"verlist_skin");
        strcpy((*dict)[4].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\brnch_root_list_opt{on,off} */
        strcpy((*dict)[5].error_mes,"on,off,pare_down");
        strcpy((*dict)[5].keyword,"brnch_root_list_opt");
        strcpy((*dict)[5].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\brnch_root_list_skin{#} */
        strcpy((*dict)[6].error_mes,"a number > 0 ");
        strcpy((*dict)[6].keyword,"brnch_root_list_skin");
        strcpy((*dict)[6].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\brnch_root_cutoff{#} */
        strcpy((*dict)[7].error_mes,"a number > 0 ");
        strcpy((*dict)[7].keyword,"brnch_root_cutoff");
        strcpy((*dict)[7].keyarg,"1.112");
  /*-----------------------------------------------------------------------*/
  /*  8)\lnkcell_cell_divs{#} */
        strcpy((*dict)[8].error_mes,"a number > 0 ");
        strcpy((*dict)[8].keyword,"lnkcell_cell_divs");
        strcpy((*dict)[8].keyarg,"7");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\lnk_cell_excl_safe{#} */
        strcpy((*dict)[9].error_mes,"a number > 0 ");
        strcpy((*dict)[9].keyword,"lnk_cell_excl_safe");
        strcpy((*dict)[9].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /* 10)\lnk_cell_vol_safe{#} */
        strcpy((*dict)[10].error_mes,"a number > 0 ");
        strcpy((*dict)[10].keyword,"lnk_cell_vol_safe");
        strcpy((*dict)[10].keyarg,"1.05");
  /*-----------------------------------------------------------------------*/ 
  /* 11)\lnkcell_force_odd{on,off} */
        strcpy((*dict)[11].error_mes,"on,off");
        strcpy((*dict)[11].keyword,"lnkcell_force_odd");
        strcpy((*dict)[11].keyarg,"on");
  /*-----------------------------------------------------------------------*/ 
  /* 12)\shave_skin_opt{#} */
        strcpy((*dict)[12].error_mes,"on or off");
        strcpy((*dict)[12].keyword,"shave_skin_opt");
#ifndef NO_PRAGMA
        strcpy((*dict)[12].keyarg,"off");
#else
        strcpy((*dict)[12].keyarg,"on");
#endif
  /*-----------------------------------------------------------------------*/
  /* 13)\neighbor_list{ver_list,lnk_list,no_list} */
        strcpy((*dict)[13].error_mes,"ver_list,lnk_list,no_list");
        strcpy((*dict)[13].keyword,"neighbor_list");
        strcpy((*dict)[13].keyarg,"ver_list");
  /*-----------------------------------------------------------------------*/ 
  /* 14)\update_type{no_list,lnk_list} */
        strcpy((*dict)[14].error_mes,"no_list,lnk_list");
        strcpy((*dict)[14].keyword,"update_type");
        strcpy((*dict)[14].keyarg,"lnk_list");
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_cp(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 37;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up cp flags                                                   */
/*------------------------------------------------------------------------*/ 
  /*  1)\cp_vxc_typ{pz_lda,pz_lsda,pw_lda,pw_lsda,pade_lda,pade_lsda}     */
        strcpy((*dict)[1].error_mes,
	       "pz_lda,pz_lsda,pw_lda,pw_lsda,pade_lda,pade_lsda");
        strcpy((*dict)[1].keyword,"cp_vxc_typ");
        strcpy((*dict)[1].keyarg,"pz_lda");
 /*-----------------------------------------------------------------------*/ 
 /*  2)\cp_ggax_typ{becke,pw91x,fila_1x,fil_2x,pbe_x,revpbe_x,rpbe_x,xpbe_x,brx89,brx2k,debug97x,off}*/
        strcpy((*dict)[2].error_mes,
        "becke,pw91x,fila_1x,fila_2x,brx89,brx2k,pbe_x,revpbe_x,rpbe_x,xpbe_x");
        strcpy((*dict)[2].keyword,"cp_ggax_typ");
        strcpy((*dict)[2].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\cp_ggac_typ{lyp,lypm1,pw91c,pbe_c,xpbe_c,tau1_c,off}                     */
        strcpy((*dict)[3].error_mes,
        "lyp,lypm1,pw91c,pbe_c,xpbe_c,tau1_c:consistent with gga_lda/gga_lsda opt");
        strcpy((*dict)[3].keyword,"cp_ggac_typ");
        strcpy((*dict)[3].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\cp_sic{on,off} */
        strcpy((*dict)[4].error_mes,"on,off");
        strcpy((*dict)[4].keyword,"cp_sic");
        strcpy((*dict)[4].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\cp_e_e_interact{on,off} */
        strcpy((*dict)[5].error_mes,"on,off");
        strcpy((*dict)[5].keyword,"cp_e_e_interact");
        strcpy((*dict)[5].keyarg,"on");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\cp_norb{full_ortho,norm_only,no_constrnt,off} */
        strcpy((*dict)[6].error_mes,"full_ortho,norm_only,no_constrnt,off");
        strcpy((*dict)[6].keyword,"cp_norb");
        strcpy((*dict)[6].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\cp_gauss{on,off} */
        strcpy((*dict)[7].error_mes,"on,off");
        strcpy((*dict)[7].keyword,"cp_gauss");
        strcpy((*dict)[7].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 8)\cp_dft_typ{lda,lsda,gga_lda,gga_lsda} */
        strcpy((*dict)[8].error_mes,"lda,lsda,gga_lda,gga_lsda");
        strcpy((*dict)[8].keyword,"cp_dft_typ");
        strcpy((*dict)[8].keyarg,"lda");
  /*-----------------------------------------------------------------------*/ 
  /* 9)\cp_nl_list{on,off} */
        strcpy((*dict)[9].error_mes,"on,off");
        strcpy((*dict)[9].keyword,"cp_nl_list");
        strcpy((*dict)[9].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 10)\cp_mass_tau_def{#} */
        strcpy((*dict)[10].error_mes,"a number > 0 ");
        strcpy((*dict)[10].keyword,"cp_mass_tau_def");
        strcpy((*dict)[10].keyarg,"25");
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_mass_cut_def{#} */
        strcpy((*dict)[11].error_mes,"a number > 0 ");
        strcpy((*dict)[11].keyword,"cp_mass_cut_def");
        strcpy((*dict)[11].keyarg,"2");
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_energy_cut_def{#} */
        strcpy((*dict)[12].error_mes,"a number > 0 ");
        strcpy((*dict)[12].keyword,"cp_energy_cut_def");
        strcpy((*dict)[12].keyarg,"2");
  /*-----------------------------------------------------------------------*/ 
  /* 13)\cp_fict_KE{#} */
        strcpy((*dict)[13].error_mes,"a number > 0 ");
        strcpy((*dict)[13].keyword,"cp_fict_KE");
        strcpy((*dict)[13].keyarg,"1000.0");
  /*-----------------------------------------------------------------------*/ 
  /* 14)\cp_ptens{on,off} */
        strcpy((*dict)[14].error_mes,"on,off");
        strcpy((*dict)[14].keyword,"cp_ptens");
        strcpy((*dict)[14].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 15)\cp_init_orthog{on,off} */
        strcpy((*dict)[15].error_mes,"on,off");
        strcpy((*dict)[15].keyword,"cp_init_orthog");
        strcpy((*dict)[15].keyarg,"off");
   /*-----------------------------------------------------------------------*/ 
   /* 16)\cp_cg_line_min_len{#} */
        strcpy((*dict)[16].error_mes,"a number >= 3");
        strcpy((*dict)[16].keyword,"cp_cg_line_min_len");
        strcpy((*dict)[16].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_minimize_typ{min_std,min_cg,min_diis} */
        strcpy((*dict)[17].error_mes,"min_std,min_cg,min_diis");
        strcpy((*dict)[17].keyword,"cp_minimize_typ");
        strcpy((*dict)[17].keyarg,"min_std");
  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_diis_hist_len{#} */
        strcpy((*dict)[18].error_mes,"a number >= 0");
        strcpy((*dict)[18].keyword,"cp_diis_hist_len");
        strcpy((*dict)[18].keyarg,"10");
  /*-----------------------------------------------------------------------*/ 
  /* 19)\cp_orth_meth{gram_schmidt,lowdin,normalize,none} */
        strcpy((*dict)[19].error_mes,"gram_schmidt,lowdin,normalize,none");
        strcpy((*dict)[19].keyword,"cp_orth_meth");
        strcpy((*dict)[19].keyarg,"gram_schmidt");
  /*-----------------------------------------------------------------------*/ 
  /* 20)\cp_restart_type{gen_wave,initial,restart_pos,restart_posvel,restart_all}*/
        strcpy((*dict)[20].error_mes,"initial");
        strcpy((*dict)[20].keyword,"cp_restart_type");
        strcpy((*dict)[20].keyarg,"gen_wave");
  /*-----------------------------------------------------------------------*/ 
  /*  21)\diis_hist_len{#} */
        strcpy((*dict)[21].error_mes,"a number >= 0");
        strcpy((*dict)[21].keyword,"diis_hist_len");
        strcpy((*dict)[21].keyarg,"10");
  /*-----------------------------------------------------------------------*/ 
  /* 22)\nlvps_list_skin{#}   */
        strcpy((*dict)[22].error_mes,"a number > 0");
        strcpy((*dict)[22].keyword,"nlvps_list_skin");
        strcpy((*dict)[22].keyarg,"50.0");
  /*-----------------------------------------------------------------------*/ 
  /* 23)\gradient_cutoff{#}   */
        strcpy((*dict)[23].error_mes,"a number > 0");
        strcpy((*dict)[23].keyword,"gradient_cutoff");
        strcpy((*dict)[23].keyarg,"5.0e-05");
  /*-----------------------------------------------------------------------*/
  /* 24)\zero_cp_vel{initial,periodic,no}   */
        strcpy((*dict)[24].error_mes,"initial periodic,,no");
        strcpy((*dict)[24].keyword,"zero_cp_vel");
        strcpy((*dict)[24].keyarg,"no");

  /*-----------------------------------------------------------------------*/
  /* 25)\cp_check_clus{#}   */
        strcpy((*dict)[25].error_mes,"on,off");
        strcpy((*dict)[25].keyword,"cp_check_perd_size");
        strcpy((*dict)[25].keyarg,"on");

  /*-----------------------------------------------------------------------*/
  /* 26)\cp_tol_edge_dist{#}   */
        strcpy((*dict)[26].error_mes,"a number > 0");
        strcpy((*dict)[26].keyword,"cp_tol_edge_dist");
        strcpy((*dict)[26].keyarg,"3");
  /*-----------------------------------------------------------------------*/
  /* 27)\cp_para_typ{#}   */
        strcpy((*dict)[27].error_mes,"hybrid,full_g");
        strcpy((*dict)[27].keyword,"cp_para_typ");
        strcpy((*dict)[27].keyarg,"full_g");
  /*-----------------------------------------------------------------------*/ 
  /* 28)\cp_dual_grid_opt{#} */
        strcpy((*dict)[28].error_mes,"prop,not_prop or off");
        strcpy((*dict)[28].keyword,"cp_dual_grid_opt");
        strcpy((*dict)[28].keyarg,"off");

  /*-----------------------------------------------------------------------*/
  /* 29)\cp_check_dual_size{#}   */
        strcpy((*dict)[29].error_mes,"on,off");
        strcpy((*dict)[29].keyword,"cp_check_dual_size");
        strcpy((*dict)[29].keyarg,"on");
  /*-----------------------------------------------------------------------*/
  /* 30)\cp_move_dual_box_opt{#}   */
        strcpy((*dict)[30].error_mes,"on,off");
        strcpy((*dict)[30].keyword,"cp_move_dual_box_opt");
        strcpy((*dict)[30].keyarg,"off");

  /*-----------------------------------------------------------------------*/ 
  /* 31)\cp_energy_cut_dual_grid_def{#} */
        strcpy((*dict)[31].error_mes,"a number > 0 ");
        strcpy((*dict)[31].keyword,"cp_energy_cut_dual_grid_def");
        strcpy((*dict)[31].keyarg,"2");
  /*-----------------------------------------------------------------------*/ 
  /* 32)\cp_alpha_conv_dual{#} */
        strcpy((*dict)[32].error_mes,"a number > 0 ");
        strcpy((*dict)[32].keyword,"cp_alpha_conv_dual");
        strcpy((*dict)[32].keyarg,"7");

  /*-----------------------------------------------------------------------*/ 
  /* 33)\interp_pme_dual{#} */
        strcpy((*dict)[33].error_mes,"an even number >= 4 ");
        strcpy((*dict)[33].keyword,"interp_pme_dual");
        strcpy((*dict)[33].keyarg,"4");
  /*-----------------------------------------------------------------------*/ 
  /* 34)\cp_elf_calc_frq{#} */
        strcpy((*dict)[34].error_mes,"a number >= 0");
        strcpy((*dict)[34].keyword,"cp_elf_calc_frq");
        strcpy((*dict)[34].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 35)\cp_ngrid_skip{#} */
        strcpy((*dict)[35].error_mes,"a number > 0");
        strcpy((*dict)[35].keyword,"cp_ngrid_skip");
        strcpy((*dict)[35].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /* 36)\cp_isok_opt{#} */
        strcpy((*dict)[36].error_mes,"on,off");
        strcpy((*dict)[36].keyword,"cp_isok_opt");
        strcpy((*dict)[36].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 37)\cp_hess_cut{#} */
        strcpy((*dict)[37].error_mes,"a number >= 0");
        strcpy((*dict)[37].keyword,"cp_hess_cut");
        strcpy((*dict)[37].keyarg,"1.5");
  /*-----------------------------------------------------------------------*/ 

/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_gen(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 25;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up gen flags                                                  */
/*------------------------------------------------------------------------*/ 
  /*  1)\simulation_typ{md,minimize,cp,cp_wave,cp_min,cp_wave_min,
                      debug,debug_cp}  */
        strcpy((*dict)[1].error_mes,
"minimize,md,pimd,cp,cp_pimd,cp_wave,cp_wave_pimd,cp_min,cp_wave_min,cp_wave_min_pimd");
        strcpy((*dict)[1].keyword,"simulation_typ");
        strcpy((*dict)[1].keyarg,"md");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\ensemble_typ{nve,nvt,npt_i,npt_f,nst} */
        strcpy((*dict)[2].error_mes,"nve,nvt,npt_i,npt_f,nst");
        strcpy((*dict)[2].keyword,"ensemble_typ");
        strcpy((*dict)[2].keyarg,"nve");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\num_time_step{#} */
        strcpy((*dict)[3].error_mes,"a number > 0 ");
        strcpy((*dict)[3].keyword,"num_time_step");
        strcpy((*dict)[3].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\time_step{#} */
        strcpy((*dict)[4].error_mes,"a number > 0 ");
        strcpy((*dict)[4].keyword,"time_step");
        strcpy((*dict)[4].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\temperature{#} */
        strcpy((*dict)[5].error_mes,"a number > 0 ");
        strcpy((*dict)[5].keyword,"temperature");
        strcpy((*dict)[5].keyarg,"300");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\pressure{#} */
        strcpy((*dict)[6].error_mes,"a number > 0 ");
        strcpy((*dict)[6].keyword,"pressure");
        strcpy((*dict)[6].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\restart_type{initial,restart_pos,restart_posvel,restart_all}*/
        strcpy((*dict)[7].error_mes,"initial");
        strcpy((*dict)[7].keyword,"restart_type");
        strcpy((*dict)[7].keyarg,"initial");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\minimize_typ{min_std,min_cg,min_diis} */
        strcpy((*dict)[8].error_mes,"min_std,min_cg,min_diis");
        strcpy((*dict)[8].keyword,"minimize_typ");
        strcpy((*dict)[8].keyarg,"min_std");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\annealing_rate{#} */
        strcpy((*dict)[9].error_mes,"a number > 0");
        strcpy((*dict)[9].keyword,"annealing_rate");
        strcpy((*dict)[9].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  10)\num_proc_beads{#} */
        strcpy((*dict)[10].error_mes,"a number > 0");
        strcpy((*dict)[10].keyword,"num_proc_beads");
        strcpy((*dict)[10].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  11)\num_proc_states{#} */
        strcpy((*dict)[11].error_mes,"a number > 0");
        strcpy((*dict)[11].keyword,"num_proc_states");
        strcpy((*dict)[11].keyarg,"1");
  /*-----------------------------------------------------------------------*/
  /*  12)\num_proc_class_forc{#} */
        strcpy((*dict)[12].error_mes,"a number > 0");
        strcpy((*dict)[12].keyword,"num_proc_class_forc");
        strcpy((*dict)[12].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  13)\num_proc_tot{#} */
        strcpy((*dict)[13].error_mes,
        "a number > 0 and equal to total number of processors to be used ");
        strcpy((*dict)[13].keyword,"num_proc_tot");
        strcpy((*dict)[13].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  14)\rndm_seed{#} */
        strcpy((*dict)[14].error_mes,"a number > 0 ");
        strcpy((*dict)[14].keyword,"iseed");
        strcpy((*dict)[14].keyarg,"19541");
  /*-----------------------------------------------------------------------*/ 
  /*  15)\rndm_seed2{#} */
        strcpy((*dict)[15].error_mes,"a number > 0 ");
        strcpy((*dict)[15].keyword,"iseed2");
        strcpy((*dict)[15].keyarg,"92571");
  /*-----------------------------------------------------------------------*/ 
  /* 16)\generic_fft_opt{on,off} */
        strcpy((*dict)[16].error_mes,"on or off");
        strcpy((*dict)[16].keyword,"generic_fft_opt");
        strcpy((*dict)[16].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\gen_alpha_clus{#}   */
        strcpy((*dict)[17].error_mes,"a number > 0");
        strcpy((*dict)[17].keyword,"gen_alpha_clus");
        strcpy((*dict)[17].keyarg,"7");
  /*-----------------------------------------------------------------------*/
  /* 18)\gen_ecut_clus{#}   */
        strcpy((*dict)[18].error_mes,"a number > 0");
        strcpy((*dict)[18].keyword,"gen_ecut_clus");
        strcpy((*dict)[18].keyarg,"20");
  /*-----------------------------------------------------------------------*/ 
  /*-----------------------------------------------------------------------*/
  /* 19)\surf_tens{#} */
        strcpy((*dict)[19].error_mes,"none");
        strcpy((*dict)[19].keyword,"surf_tens");
        strcpy((*dict)[19].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 20)\min_num_atoms_per_proc{#}   */
        strcpy((*dict)[20].error_mes,"to use option: a number > 0 otherwise = 0");
        strcpy((*dict)[20].keyword,"min_num_atoms_per_proc");
        strcpy((*dict)[20].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  21)\num_proc_class_forc_src{#} */
        strcpy((*dict)[21].error_mes,"a number > 0");
        strcpy((*dict)[21].keyword,"num_proc_class_forc_src");
        strcpy((*dict)[21].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  22)\num_proc_class_forc_trg{#} */
        strcpy((*dict)[22].error_mes,"a number > 0");
        strcpy((*dict)[22].keyword,"num_proc_class_forc_trg");
        strcpy((*dict)[22].keyarg,"1");

  /*-----------------------------------------------------------------------*/ 
  /*  23)\annealing_opt{on,off} */
        strcpy((*dict)[23].error_mes,"on or off");
        strcpy((*dict)[23].keyword,"annealing_opt");
        strcpy((*dict)[23].keyarg,"off");

  /*-----------------------------------------------------------------------*/ 
  /*  24)\annealing_start_temperature{#} */
        strcpy((*dict)[24].error_mes,"a number > 0");
        strcpy((*dict)[24].keyword,"annealing_start_temperature");
        strcpy((*dict)[24].keyarg,"0.1");

  /*-----------------------------------------------------------------------*/ 
  /*  25)\annealing_target_temperature{#} */
        strcpy((*dict)[25].error_mes,"a number > 0");
        strcpy((*dict)[25].keyword,"annealing_target_temperature");
        strcpy((*dict)[25].keyarg,"300");

/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_vpot(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 25;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up vpot flags                                                 */
/*------------------------------------------------------------------------*/ 
  /*  1)\shift_inter_pe{on,off,swit} */
        strcpy((*dict)[1].error_mes,"on,off,swit");
        strcpy((*dict)[1].keyword,"shift_inter_pe");
        strcpy((*dict)[1].keyarg,"swit");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\inter_spline_pts{#} */
        strcpy((*dict)[2].error_mes,"a number > 0 ");
        strcpy((*dict)[2].keyword,"inter_spline_pts");
        strcpy((*dict)[2].keyarg,"2000");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\intra_block_min{} */
        strcpy((*dict)[3].error_mes,"a number >= 0");
        strcpy((*dict)[3].keyword,"intra_block_min");
        strcpy((*dict)[3].keyarg,"50");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\pten_inter_respa{#} */
        strcpy((*dict)[4].error_mes,"a number");
        strcpy((*dict)[4].keyword,"pten_inter_respa");
        strcpy((*dict)[4].keyarg,"0.0");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\pten_kin_respa{#} */
        strcpy((*dict)[5].error_mes,"a number");
        strcpy((*dict)[5].keyword,"pten_kin_respa");
        strcpy((*dict)[5].keyarg,"0.0");
  /*-----------------------------------------------------------------------*/
  /*  6)\pseud_spline_pts{#} */
        strcpy((*dict)[6].error_mes,"a number > 0 ");
        strcpy((*dict)[6].keyword,"pseud_spline_pts");
        strcpy((*dict)[6].keyarg,"4000");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\scratch_length{#} */
        strcpy((*dict)[7].error_mes,"a number > 0 ");
        strcpy((*dict)[7].keyword,"scratch_length");
        strcpy((*dict)[7].keyarg,"1000");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\ewald_alpha{#} */
        strcpy((*dict)[8].error_mes,"a number > 0 ");
        strcpy((*dict)[8].keyword,"ewald_alpha");
        strcpy((*dict)[8].keyarg,"7");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\ewald_kmax{#} */
        strcpy((*dict)[9].error_mes,"a number > 0 ");
        strcpy((*dict)[9].keyword,"ewald_kmax");
        strcpy((*dict)[9].keyarg,"7");
  /*-----------------------------------------------------------------------*/ 
  /* 10)\ewald_respa_kmax{#} */
        strcpy((*dict)[10].error_mes,"a number > 0");
        strcpy((*dict)[10].keyword,"ewald_respa_kmax");
        strcpy((*dict)[10].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 11)\ewald_pme_opt{#} */
        strcpy((*dict)[11].error_mes,"on or off");
        strcpy((*dict)[11].keyword,"ewald_pme_opt");
        strcpy((*dict)[11].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 12)\ewald_kmax_pme{#} */
        strcpy((*dict)[12].error_mes,"a number > 0 ");
        strcpy((*dict)[12].keyword,"ewald_kmax_pme");
        strcpy((*dict)[12].keyarg,"7");
  /*-----------------------------------------------------------------------*/ 
  /* 13)\ewald_interp_pme{#} */
        strcpy((*dict)[13].error_mes,"an even number >= 4 ");
        strcpy((*dict)[13].keyword,"ewald_interp_pme");
        strcpy((*dict)[13].keyarg,"4");
  /*-----------------------------------------------------------------------*/ 
  /* 14)\ewald_respa_pme_opt{#} */
        strcpy((*dict)[14].error_mes,"on or off");
        strcpy((*dict)[14].keyword,"ewald_respa_pme_opt");
        strcpy((*dict)[14].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 15)\ewald_respa_kmax_pme{#} */
        strcpy((*dict)[15].error_mes,"a number > 0 ");
        strcpy((*dict)[15].keyword,"ewald_respa_kmax_pme");
        strcpy((*dict)[15].keyarg,"7");
  /*-----------------------------------------------------------------------*/ 
  /* 16)\ewald_respa_interp_pme{#} */
        strcpy((*dict)[16].error_mes,"an even number >= 4 ");
        strcpy((*dict)[16].keyword,"ewald_respa_interp_pme");
        strcpy((*dict)[16].keyarg,"4");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\sep_VanderWaals{on,off} */
        strcpy((*dict)[17].error_mes,"on,off");
        strcpy((*dict)[17].keyword,"sep_VanderWaals");
        strcpy((*dict)[17].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 18)\dielectric_opt{on,off} */
        strcpy((*dict)[18].error_mes,"on,off");
        strcpy((*dict)[18].keyword,"dielectric_opt");
        strcpy((*dict)[18].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 19)\dielectric_rheal{#} */
        strcpy((*dict)[19].error_mes,"a number > 0 ");
        strcpy((*dict)[19].keyword,"dielectric_rheal");
        strcpy((*dict)[19].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /* 20)\dielectric_cut{#} */
        strcpy((*dict)[20].error_mes,"a number > 0 ");
        strcpy((*dict)[20].keyword,"dielectric_cut");
        strcpy((*dict)[20].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /*  21)\dielectric_eps{#} */
        strcpy((*dict)[21].error_mes,"a number > 1 ");
        strcpy((*dict)[21].keyword,"dielectric_eps");
        strcpy((*dict)[21].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /*  22)\std_intra_block{on,off} */
        strcpy((*dict)[22].error_mes,"on,off");
        strcpy((*dict)[22].keyword,"std_intra_block");
        strcpy((*dict)[22].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  23)\con_intra_block{on,off} */
        strcpy((*dict)[23].error_mes,"on,off");
        strcpy((*dict)[23].keyword,"con_intra_block");
        strcpy((*dict)[23].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  24)\inter_PE_calc_freq{#} */
        strcpy((*dict)[24].error_mes,"a number > 0 ");
        strcpy((*dict)[24].keyword,"inter_PE_calc_freq");
        strcpy((*dict)[24].keyarg,"5");
  /*-----------------------------------------------------------------------*/ 
  /*  25)\pme_parallel_opt{#} */
        strcpy((*dict)[25].error_mes,"none,hybrid,full_g");
        strcpy((*dict)[25].keyword,"pme_parallel_opt");
        strcpy((*dict)[25].keyarg,"none");
  /*-----------------------------------------------------------------------*/ 
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_run(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 33;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up run flags                                                  */
/*------------------------------------------------------------------------*/ 
  /*  1)\init_resmpl_atm_vel{on,off} */
        strcpy((*dict)[1].error_mes,"on,off");
        strcpy((*dict)[1].keyword,"init_resmp_atm_vel");
        strcpy((*dict)[1].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 2)\init_resmpl_cp_vel{on,off} */
        strcpy((*dict)[2].error_mes,"on,off");
        strcpy((*dict)[2].keyword,"init_resmpl_cp_vel");
        strcpy((*dict)[2].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 3)\init_resmpl_cp_nhc{on,off} */
        strcpy((*dict)[3].error_mes,"on,off");
        strcpy((*dict)[3].keyword,"init_resmpl_cp_nhc");
        strcpy((*dict)[3].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 4)\resmpl_frq_atm_vel{#} */
        strcpy((*dict)[4].error_mes,"a number >= 0");
        strcpy((*dict)[4].keyword,"resmpl_frq_atm_vel");
        strcpy((*dict)[4].keyarg,"0");
  /*-----------------------------------------------------------------------*/
  /* 25)\rescale_frq_atm_vel{#} */
        strcpy((*dict)[25].error_mes,"a number >= 0");
        strcpy((*dict)[25].keyword,"rescale_frq_atm_vel");
        strcpy((*dict)[25].keyarg,"0");
  /*-----------------------------------------------------------------------*/
  /* 5)\respa_steps_lrf{#} */
        strcpy((*dict)[5].error_mes,"a number >=0");
        strcpy((*dict)[5].keyword,"respa_steps_lrf");
        strcpy((*dict)[5].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 6)\respa_steps_torsion{#} */
        strcpy((*dict)[6].error_mes,"a number >=0");
        strcpy((*dict)[6].keyword,"respa_steps_torsion");
        strcpy((*dict)[6].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\respa_steps_intra{#} */
        strcpy((*dict)[7].error_mes,"a number >= 0");
        strcpy((*dict)[7].keyword,"respa_steps_intra");
        strcpy((*dict)[7].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\respa_rheal{#} */
        strcpy((*dict)[8].error_mes,"a number > 0 ");
        strcpy((*dict)[8].keyword,"respa_rheal");
        strcpy((*dict)[8].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\shake_tol{#} */
        strcpy((*dict)[9].error_mes,"a number > 0 ");
        strcpy((*dict)[9].keyword,"shake_tol");
        strcpy((*dict)[9].keyarg,"1e-6");
  /*-----------------------------------------------------------------------*/ 
  /* 10)\rattle_tol{#} */
        strcpy((*dict)[10].error_mes,"a number > 0 ");
        strcpy((*dict)[10].keyword,"rattle_tol");
        strcpy((*dict)[10].keyarg,"1e-6");
  /*-----------------------------------------------------------------------*/ 
  /* 11)\max_constrnt_iter{#} */
        strcpy((*dict)[11].error_mes,"a number > 0 ");
        strcpy((*dict)[11].keyword,"max_constrnt_iter");
        strcpy((*dict)[11].keyarg,"200");
  /*-----------------------------------------------------------------------*/ 
  /* 12)\init_rescale_atm_vel{on,off} */
        strcpy((*dict)[12].error_mes,"on,off");
        strcpy((*dict)[12].keyword,"init_rescale_atm_vel");
        strcpy((*dict)[12].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 13)\init_rescale_atm_nhc{on,off} */
        strcpy((*dict)[13].error_mes,"on,off");
        strcpy((*dict)[13].keyword,"init_rescale_atm_nhc");
        strcpy((*dict)[13].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 14)\init_rscale_cp_vel{on,off} */
        strcpy((*dict)[14].error_mes,"on,off");
        strcpy((*dict)[14].keyword,"init_rescale_cp_vel");
        strcpy((*dict)[14].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /* 15)\resmpl_frq_cp_vel{ # } */
        strcpy((*dict)[15].error_mes,"a number >= 0");
        strcpy((*dict)[15].keyword,"resmpl_cp_vel");
        strcpy((*dict)[15].keyarg,"0");
  /*-----------------------------------------------------------------------*/
  /* 16)\group_con_tol{#} */
        strcpy((*dict)[16].error_mes,"a number > 0 ");
        strcpy((*dict)[16].keyword,"group_con_tol");
        strcpy((*dict)[16].keyarg,"1e-6");
  /*-----------------------------------------------------------------------*/ 
  /* 17)\cp_norb_tol */
        strcpy((*dict)[17].error_mes,"a number > 0 ");
        strcpy((*dict)[17].keyword,"cp_norb_tol");
        strcpy((*dict)[17].keyarg,"1.e-3");
  /*-----------------------------------------------------------------------*/ 
  /* 18)\cp_ks_rot_frq{#} */
        strcpy((*dict)[18].error_mes,"a number >=0");
        strcpy((*dict)[18].keyword,"cp_ks_rot");
        strcpy((*dict)[18].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 19)\cp_shake_tol{#} */
        strcpy((*dict)[19].error_mes,"a number > 0 ");
        strcpy((*dict)[19].keyword,"cp_shake_tol");
        strcpy((*dict)[19].keyarg,"1e-6");
  /*-----------------------------------------------------------------------*/ 
  /* 20)\cp_rattle_tol{#} */
        strcpy((*dict)[20].error_mes,"a number > 0 ");
        strcpy((*dict)[20].keyword,"cp_rattle_tol");
        strcpy((*dict)[20].keyarg,"1e-6");
  /*-----------------------------------------------------------------------*/ 
  /* 21)\cp_run_tol{#} */
        strcpy((*dict)[21].error_mes,"a number > 0 ");
        strcpy((*dict)[21].keyword,"cp_run_tol");
        strcpy((*dict)[21].keyarg,"2.0");
  /*-----------------------------------------------------------------------*/ 
  /* 22)\zero_com_vel{yes,no}*/
        strcpy((*dict)[22].error_mes,"yes,no");
        strcpy((*dict)[22].keyword,"zero_com_vel");
        strcpy((*dict)[22].keyarg,"no");
  /*-----------------------------------------------------------------------*/ 
  /* 23)\min_tol{#} */
        strcpy((*dict)[23].error_mes,"a number > 0 ");
        strcpy((*dict)[23].keyword,"min_tol");
        strcpy((*dict)[23].keyarg,"0.0002");
  /*-----------------------------------------------------------------------*/ 
  /* 24)\cp_min_tol{#} */
        strcpy((*dict)[24].error_mes,"a number > 0 ");
        strcpy((*dict)[24].keyword,"cp_min_tol");
        strcpy((*dict)[24].keyarg,"0.0002");
  /*-----------------------------------------------------------------------*/ 
  /* 26)\hess_opt{full_an,full_num,diag,unit} */
        strcpy((*dict)[26].error_mes,"full_an,full_num,diag,unit ");
        strcpy((*dict)[26].keyword,"hess_opt");
        strcpy((*dict)[26].keyarg,"unit");
  /*-----------------------------------------------------------------------*/ 
  /* 27)\class_mass_scale_fact{#} */
        strcpy((*dict)[27].error_mes,"a number > 0 ");
        strcpy((*dict)[27].keyword,"class_mass_scale_fact");
        strcpy((*dict)[27].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /* 28)\hmat_int_typ{#} */
        strcpy((*dict)[28].error_mes,"normal,upper_triangle");
        strcpy((*dict)[28].keyword,"hmat_int_typ");
        strcpy((*dict)[28].keyarg,"normal");
  /*-----------------------------------------------------------------------*/ 
  /* 29)\hmat_cons_typ{#} */
        strcpy((*dict)[29].error_mes,"none,ortho_rhom,mono_clin");
        strcpy((*dict)[29].keyword,"hmat_cons_typ");
        strcpy((*dict)[29].keyarg,"none");
  /*-----------------------------------------------------------------------*/ 
  /* 30)\min_atm_com_fix{#} */
        strcpy((*dict)[30].error_mes,"yes or no");
        strcpy((*dict)[30].keyword,"min_atm_com_fix");
        strcpy((*dict)[30].keyarg,"no");
  /*-----------------------------------------------------------------------*/ 
  /* 31)\rescale_frq_cp_vel{ # } */
        strcpy((*dict)[31].error_mes,"a number >= 0");
        strcpy((*dict)[31].keyword,"rescale_frq_cp_vel");
        strcpy((*dict)[31].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /* 32)\auto_rescale_cp_vel{ on/off } */
        strcpy((*dict)[32].error_mes,"on/off");
        strcpy((*dict)[32].keyword,"auto_rescale_cp_vel");
        strcpy((*dict)[32].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /* 33)\auto_rescale_cp_vel_tol{ # } */
        strcpy((*dict)[33].error_mes,"a number > 1");
        strcpy((*dict)[33].keyword,"auto_rescale_cp_vel_tol");
        strcpy((*dict)[33].keyarg,"9.0");
  /*-----------------------------------------------------------------------*/ 
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_nhc(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict =15;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up nhc flags                                                  */
/*------------------------------------------------------------------------*/ 
  /* 1)\respa_steps_nhc{#} */
        strcpy((*dict)[1].error_mes,"a number > 0 ");
        strcpy((*dict)[1].keyword,"respa_steps_nhc");
        strcpy((*dict)[1].keyarg,"2");
  /*-----------------------------------------------------------------------*/ 
  /* 2)\yosh_steps_nhc{#} */
        strcpy((*dict)[2].error_mes,"1,3,5,7");
        strcpy((*dict)[2].keyword,"yosh_steps_nhc");
        strcpy((*dict)[2].keyarg,"3");
  /*-----------------------------------------------------------------------*/ 
  /* 3)\atm_nhc_tau_def{#} */
        strcpy((*dict)[3].error_mes,"a number > 0 ");
        strcpy((*dict)[3].keyword,"atm_nhc_tau_def");
        strcpy((*dict)[3].keyarg,"1000");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\cp_nhc_len{#} */
        strcpy((*dict)[4].error_mes,"a number > 0 ");
        strcpy((*dict)[4].keyword,"cp_nhc_len");
        strcpy((*dict)[4].keyarg,"3");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\resmpl_frq_cp_nhc{ #} */
        strcpy((*dict)[5].error_mes,"on,off");
        strcpy((*dict)[5].keyword,"resmpl_cp_nhc");
        strcpy((*dict)[5].keyarg,"0");
  /*-----------------------------------------------------------------------*/
  /*  6)\init_resmpl_atm_nhc{on,off} */
        strcpy((*dict)[6].error_mes,"on,off");
        strcpy((*dict)[6].keyword,"init_resmp_atm_nhc");
        strcpy((*dict)[6].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /*  7)\init_rescale_cp_nhc{on,off} */
        strcpy((*dict)[7].error_mes,"on,off");
        strcpy((*dict)[7].keyword,"init_rescale_cp_nhc");
        strcpy((*dict)[7].keyarg,"off");
  /*-----------------------------------------------------------------------*/
  /*  8)\resmpl_atm_nhc{#} */
        strcpy((*dict)[8].error_mes,"a number >= 0");
        strcpy((*dict)[8].keyword,"resmpl_atm_nhc");
        strcpy((*dict)[8].keyarg,"0");
  /*-----------------------------------------------------------------------*/
  /*  9)\atm_nhc_len{#} */
        strcpy((*dict)[9].error_mes,"a number > 0 ");
        strcpy((*dict)[9].keyword,"atm_nhc_len");
        strcpy((*dict)[9].keyarg,"2");
  /*-----------------------------------------------------------------------*/ 
  /* 10)\cp_nhc_tau_def{#} */
        strcpy((*dict)[10].error_mes,"a number > 0 ");
        strcpy((*dict)[10].keyword,"cp_nhc_tau_def");
        strcpy((*dict)[10].keyarg,"25");
  /*-----------------------------------------------------------------------*/ 
  /* 11)\cp_respa_steps_nhc{#} */
        strcpy((*dict)[11].error_mes,"a number > 0 ");
        strcpy((*dict)[11].keyword,"cp_respa_steps_nhc");
        strcpy((*dict)[11].keyarg,"4");
  /*-----------------------------------------------------------------------*/ 
  /* 12)\cp_yosh_steps_nhc{#} */
        strcpy((*dict)[12].error_mes,"1,3,5,7 or 9");
        strcpy((*dict)[12].keyword,"cp_yosh_steps_nhc");
        strcpy((*dict)[12].keyarg,"9");
  /*-----------------------------------------------------------------------*/ 
  /* 13)\respa_xi_opt{#} */
        strcpy((*dict)[13].error_mes,"1,2,3,4");
        strcpy((*dict)[13].keyword,"respa_xi_opt");
        strcpy((*dict)[13].keyarg,"1");
  /*---------------------------------------------------------------------*/
  /* 14)\thermostat_type{#} */
        strcpy((*dict)[14].error_mes,"choose NHC or GGMT");
        strcpy((*dict)[14].keyword,"thermostat_type");
        strcpy((*dict)[14].keyarg,"NHC");
  /*---------------------------------------------------------------------*/
  /* 15)\cp_therm_heat_fact{#} */
        strcpy((*dict)[15].error_mes,"a number >= 1");
        strcpy((*dict)[15].keyword,"cp_therm_heat_fact");
        strcpy((*dict)[15].keyarg,"1.0");
  /*---------------------------------------------------------------------*/
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_vol(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 4;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up vol flags                                                  */
/*------------------------------------------------------------------------*/ 
  /*  1)\volume_tau{#} */
        strcpy((*dict)[1].error_mes,"a number > 0 ");
        strcpy((*dict)[1].keyword,"volume_tau");
        strcpy((*dict)[1].keyarg,"1000");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\volume_nhc_tau{#} */
        strcpy((*dict)[2].error_mes,"a number > 0 ");
        strcpy((*dict)[2].keyword,"volume_nhc_tau");
        strcpy((*dict)[2].keyarg,"1000");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\periodicity{#} */
        strcpy((*dict)[3].error_mes,"0,1,2,3,0_ewald");
        strcpy((*dict)[3].keyword,"periodicity");
        strcpy((*dict)[3].keyarg,"3");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\intra_perds{on,off} */
        strcpy((*dict)[4].error_mes,"on,off");
        strcpy((*dict)[4].keyword,"intra_perds");
        strcpy((*dict)[4].keyarg,"off");
  /*-----------------------------------------------------------------------*/
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_write(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 29;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up write flags                                                */
/*------------------------------------------------------------------------*/ 
  /*  1)\write_binary_cp_coef{on,off} */
        strcpy((*dict)[1].error_mes,"on,off");
        strcpy((*dict)[1].keyword,"write_binary_cp_coef");
        strcpy((*dict)[1].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\write_force_freq{#} */
        strcpy((*dict)[2].error_mes,"a number > 0 ");
        strcpy((*dict)[2].keyword,"write_force_freq");
        strcpy((*dict)[2].keyarg,"1000000");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\write_screen_freq{#} */
        strcpy((*dict)[3].error_mes,"a number > 0 ");
        strcpy((*dict)[3].keyword,"write_screen_freq");
        strcpy((*dict)[3].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  4)\write_dump_freq{#} */
        strcpy((*dict)[4].error_mes,"a number > 0 ");
        strcpy((*dict)[4].keyword,"write_dump_freq");
        strcpy((*dict)[4].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\write_inst_freq{#} */
        strcpy((*dict)[5].error_mes,"a number > 0 ");
        strcpy((*dict)[5].keyword,"write_inst_freq");
        strcpy((*dict)[5].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\write_pos_freq{#} */
        strcpy((*dict)[6].error_mes,"a number > 0 ");
        strcpy((*dict)[6].keyword,"write_pos_freq");
        strcpy((*dict)[6].keyarg,"100");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\write_vel_freq{#} */
        strcpy((*dict)[7].error_mes,"a number > 0 ");
        strcpy((*dict)[7].keyword,"write_vel_freq");
        strcpy((*dict)[7].keyarg,"100");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\write_cp_c_freq{#} */
        strcpy((*dict)[8].error_mes,"a number > 0 ");
        strcpy((*dict)[8].keyword,"write_cp_c_freq");
        strcpy((*dict)[8].keyarg,"100");
  /*-----------------------------------------------------------------------*/ 
  /*  9)\screen_output_units */
        strcpy((*dict)[9].error_mes,"au,kcal_mol,kelvin");
        strcpy((*dict)[9].keyword,"screen_output_units");
        strcpy((*dict)[9].keyarg,"au");
  /*-----------------------------------------------------------------------*/ 
  /*  10)\conf_file_format */
        strcpy((*dict)[10].error_mes,"binary,formatted");
        strcpy((*dict)[10].keyword,"conf_file_format");
        strcpy((*dict)[10].keyarg,"formatted");
  /*-----------------------------------------------------------------------*/ 
  /*  11)\path_cent_file */
        strcpy((*dict)[11].error_mes,"none");
        strcpy((*dict)[11].keyword,"path_cent_file");
        strcpy((*dict)[11].keyarg,"sim_path_cent.out");
  /*-----------------------------------------------------------------------*/ 
  /*  12)\atm_force_file */
        strcpy((*dict)[12].error_mes,"none");
        strcpy((*dict)[12].keyword,"atm_force_file");
        strcpy((*dict)[12].keyarg,"sim_atm_force.out");
  /*------------------------------------------------------------------------*/
  /*  13)\cp_coef_file */
        strcpy((*dict)[13].error_mes,"none");
        strcpy((*dict)[13].keyword,"cp_coef_file");
        strcpy((*dict)[13].keyarg,"sim_cp_coef.out");
  /*-----------------------------------------------------------------------*/ 
  /*  14)\conf_partial_freq */
        strcpy((*dict)[14].error_mes,"A number >= 0");
        strcpy((*dict)[14].keyword,"conf_partial_freq");
        strcpy((*dict)[14].keyarg,"100");
  /*-----------------------------------------------------------------------*/ 
  /*  15)\path_cent_freq */
        strcpy((*dict)[15].error_mes,"A number >= 0");
        strcpy((*dict)[15].keyword,"path_cent_freq");
        strcpy((*dict)[15].keyarg,"100");
  /*-----------------------------------------------------------------------*/ 
  /*  16)\conf_partial_limits */
        strcpy((*dict)[16].error_mes,"Two numbers > 0 (X,Y)");
        strcpy((*dict)[16].keyword,"conf_partial_limits");
        strcpy((*dict)[16].keyarg,"1,0");
  /*-----------------------------------------------------------------------*/ 
  /*  17)\read_binary_cp_coef{on,off} */
        strcpy((*dict)[17].error_mes,"on,off");
        strcpy((*dict)[17].keyword,"read_binary_cp_coef");
        strcpy((*dict)[17].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  18)\sim_name */
        strcpy((*dict)[18].error_mes,"none");
        strcpy((*dict)[18].keyword,"sim_name");
        strcpy((*dict)[18].keyarg,"sim_input.out");
  /*-----------------------------------------------------------------------*/ 
  /*  19)\out_restart_file */
        strcpy((*dict)[19].error_mes,"none");
        strcpy((*dict)[19].keyword,"out_restart_file");
        strcpy((*dict)[19].keyarg,"sim_restart.out");
  /*-----------------------------------------------------------------------*/ 
  /*  20)\in_restart_file */
        strcpy((*dict)[20].error_mes,"none");
        strcpy((*dict)[20].keyword,"in_restart_file");
        strcpy((*dict)[20].keyarg,"sim_restart.in");
  /*-----------------------------------------------------------------------*/ 
  /*  21)\instant_file */
        strcpy((*dict)[21].error_mes,"none");
        strcpy((*dict)[21].keyword,"instant_file");
        strcpy((*dict)[21].keyarg,"sim_instant.out");
  /*-----------------------------------------------------------------------*/ 
  /*  22)\atm_pos_file */
        strcpy((*dict)[22].error_mes,"none");
        strcpy((*dict)[22].keyword,"atm_pos_file");
        strcpy((*dict)[22].keyarg,"sim_atm_pos.out");
  /*-----------------------------------------------------------------------*/ 
  /*  23)\atm_vel_file */
        strcpy((*dict)[23].error_mes,"none");
        strcpy((*dict)[23].keyword,"atm_vel_file");
        strcpy((*dict)[23].keyarg,"sim_atm_vel.out");
  /*-----------------------------------------------------------------------*/ 
  /*  24)\conf_partial_file */
        strcpy((*dict)[24].error_mes,"none");
        strcpy((*dict)[24].keyword,"conf_partial_file");
        strcpy((*dict)[24].keyarg,"sim_atm_pos_part.out");
  /*-----------------------------------------------------------------------*/ 
  /*  25)\mol_set_file */
        strcpy((*dict)[25].error_mes,"none");
        strcpy((*dict)[25].keyword,"mol_set_file");
        strcpy((*dict)[25].keyarg,"sim_mol_set.in");
  /*-----------------------------------------------------------------------*/ 
  /*  26)\cp_restart_out_file */
        strcpy((*dict)[26].error_mes,"none");
        strcpy((*dict)[26].keyword,"cp_restart_out_file");
        strcpy((*dict)[26].keyarg,"sim_cp_restart.out");
  /*-----------------------------------------------------------------------*/ 
  /*  27)\cp_restart_in_file */
        strcpy((*dict)[27].error_mes,"none");
        strcpy((*dict)[27].keyword,"cp_restart_in_file");
        strcpy((*dict)[27].keyarg,"sim_cp_restart.in");
  /*-----------------------------------------------------------------------*/ 
  /*  28)\cp_kseigs_file */
         strcpy((*dict)[28].error_mes,"none");
         strcpy((*dict)[28].keyword,"cp_kseigs_file");
         strcpy((*dict)[28].keyarg,"sim_cp_kseigs.out");
  /*-----------------------------------------------------------------------*/ 
  /*  29)\cp_elf_file */
         strcpy((*dict)[29].error_mes,"none");
         strcpy((*dict)[29].keyword,"cp_elf_file");
         strcpy((*dict)[29].keyarg,"sim_cp_elf.out");
  /*-----------------------------------------------------------------------*/ 
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_pimd(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
   { /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                                */
  int i;

/*========================================================================*/
/*  0) Malloc the dictionary                                              */ 

  *num_dict = 11;
  *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

/*========================================================================*/
/*  I) Initialize the user set option(did the user set the key word       */ 

        for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
        for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

/*========================================================================*/ 
/*II) Set up the dictionary                                               */
/*========================================================================*/ 
/*   A) Set up pimd flags                                                 */
/*------------------------------------------------------------------------*/ 
  /*  1)\path_int_beads{#} */
        strcpy((*dict)[1].error_mes,"a number > 0 ");
        strcpy((*dict)[1].keyword,"path_int_beads");
        strcpy((*dict)[1].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  2)\path_int_gamma_adb{#} */
        strcpy((*dict)[2].error_mes,"a number > 0 ");
        strcpy((*dict)[2].keyword,"path_int_gamma_adb");
        strcpy((*dict)[2].keyarg,"1.0");
  /*-----------------------------------------------------------------------*/ 
  /*  3)\path_int_md_typ{staging,centroid} */
        strcpy((*dict)[3].error_mes,"staging,centroid");
        strcpy((*dict)[3].keyword,"path_int_md_typ");
        strcpy((*dict)[3].keyarg,"centroid");
  /*-----------------------------------------------------------------------*/
  /*  4)\pi_beads_level_full{#} */
        strcpy((*dict)[4].error_mes,"a number >= 0");
        strcpy((*dict)[4].keyword,"pi_beads_level_full");
        strcpy((*dict)[4].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  5)\pi_beads_level_inter_short{#} */
        strcpy((*dict)[5].error_mes,"a number >= 0");
        strcpy((*dict)[5].keyword,"pi_beads_level_inter_short");
        strcpy((*dict)[5].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  6)\pi_beads_level_intra_res{#} */
        strcpy((*dict)[6].error_mes,"a number >= 0");
        strcpy((*dict)[6].keyword,"pi_beads_level_intra_res");
        strcpy((*dict)[6].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  7)\pi_beads_level_intra{#} */
        strcpy((*dict)[7].error_mes,"a number >= 0");
        strcpy((*dict)[7].keyword,"pi_beads_level_intra");
        strcpy((*dict)[7].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*  8)\respa_steps_pimd{#} */
        strcpy((*dict)[8].error_mes,"a number >= 0");
        strcpy((*dict)[8].keyword,"respa_steps_pimd");
        strcpy((*dict)[8].keyarg,"0");
  /*-----------------------------------------------------------------------*/ 
  /*   9)\initial_spread_size{} */
        strcpy((*dict)[9].error_mes,"a number >= 0");
        strcpy((*dict)[9].keyword,"initial_spread_size");
        strcpy((*dict)[9].keyarg,"1");
  /*-----------------------------------------------------------------------*/ 
  /*  10)\initial_spread_opt{on,off} */
        strcpy((*dict)[10].error_mes,"on,off");
        strcpy((*dict)[10].keyword,"initial_spread_opt");
        strcpy((*dict)[10].keyarg,"off");
  /*-----------------------------------------------------------------------*/ 
  /*  11)\pimd_freeze_type{centroid,all_mode} */
        strcpy((*dict)[11].error_mes,"centroid,all_mode");
        strcpy((*dict)[11].keyword,"pimd_freeze_type");
        strcpy((*dict)[11].keyarg,"centroid");
/*========================================================================*/
/*------------------------------------------------------------------------*/
/*========================================================================*/
/*                   End Subprogram:                                      */
}/*end routine*/ 
/*========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_sim_dict_velo(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
   unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 14;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;

 /*========================================================================*/
 /*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */

 /* 1) \corel_vovt{#} */
    strcpy((*dict)[1].error_mes,"on or off");
    strcpy((*dict)[1].keyword,"corel_vovt");
    strcpy((*dict)[1].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 2) \corel_vovt_ncor{#} */
    strcpy((*dict)[2].error_mes,"nrun>ncor_velo>0");
    strcpy((*dict)[2].keyword,"corel_vovt_ncor");
    strcpy((*dict)[2].keyarg,"2");
 /*-----------------------------------------------------------------------*/
 /* 3) \corel_vovt_njump{#} */
    strcpy((*dict)[3].error_mes,"njump_vovt>0");
    strcpy((*dict)[3].keyword,"corel_vovt_njump");
    strcpy((*dict)[3].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 4) \write_vovt{#} */
    strcpy((*dict)[4].error_mes,"none");
    strcpy((*dict)[4].keyword,"write_vovt");
    strcpy((*dict)[4].keyarg,"vovt.data");
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_vovt_normalize{#} */
    strcpy((*dict)[5].error_mes,"on or off");
    strcpy((*dict)[5].keyword,"corel_vovt_normalize");
    strcpy((*dict)[5].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 6) \corel_vovt_output{#} */
    strcpy((*dict)[6].error_mes,"none");
    strcpy((*dict)[6].keyword,"corel_vovt_output");
    strcpy((*dict)[6].keyarg,"atoms");
 /*-----------------------------------------------------------------------*/
 /* 7) \corel_vovt_periodic_output{#} */
    strcpy((*dict)[7].error_mes,"nwrite_vovt>0");
    strcpy((*dict)[7].keyword,"corel_vovt_periodic_output");
    strcpy((*dict)[7].keyarg,"100000");
 /*-----------------------------------------------------------------------*/
 /* 8) \corel_vovt_com{#} */
    strcpy((*dict)[8].error_mes,"on or off");
    strcpy((*dict)[8].keyword,"corel_vovt_com");
    strcpy((*dict)[8].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 9) \corel_vovt_com_ncor{#} */
    strcpy((*dict)[9].error_mes,"nrun>ncor_velo>0");
    strcpy((*dict)[9].keyword,"corel_vovt_com_ncor");
    strcpy((*dict)[9].keyarg,"2");
 /*-----------------------------------------------------------------------*/
 /* 10) \corel_vovt_com_njump{#} */
    strcpy((*dict)[10].error_mes,"njump_vovt>0");
    strcpy((*dict)[10].keyword,"corel_vovt_com_njump");
    strcpy((*dict)[10].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 11) \write_vovt_com{#} */
    strcpy((*dict)[11].error_mes,"none");
    strcpy((*dict)[11].keyword,"write_vovt_com");
    strcpy((*dict)[11].keyarg,"vovt_com.data");
 /*-----------------------------------------------------------------------*/
 /* 12) \corel_vovt_com_normalize{#} */
    strcpy((*dict)[12].error_mes,"on or off");
    strcpy((*dict)[12].keyword,"corel_vovt_com_normalize");
    strcpy((*dict)[12].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 13) \corel_vovt_com_output{#} */
    strcpy((*dict)[13].error_mes,"none");
    strcpy((*dict)[13].keyword,"corel_vovt_com_output");
    strcpy((*dict)[13].keyarg,"molecules");
 /*-----------------------------------------------------------------------*/
 /* 14) \corel_vovt_com_periodic_output{#} */
    strcpy((*dict)[14].error_mes,"nwrite_vovt>0");
    strcpy((*dict)[14].keyword,"corel_vovt_com_periodic_output");
    strcpy((*dict)[14].keyarg,"100000");
/*-----------------------------------------------------------------------*/

 /*========================================================================*/
 } /* end routine set_sim_dict_velo */
 /*========================================================================*/

void set_sim_dict_msqd(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:                                          */
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
   unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 5;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;
 /*========================================================================*/
 /*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */

 /* 1) \corel_msqd{#} */
    strcpy((*dict)[1].error_mes,"on or off");
    strcpy((*dict)[1].keyword,"corel_msqd");
    strcpy((*dict)[1].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 2) \corel_msqd_ncor{#} */
    strcpy((*dict)[2].error_mes,"nrun>ncor_msqd>0");
    strcpy((*dict)[2].keyword,"corel_msqd_ncor");
    strcpy((*dict)[2].keyarg,"2");
 /*-----------------------------------------------------------------------*/
 /* 3) \corel_msqd_njump{#} */
    strcpy((*dict)[3].error_mes,"njump_msqd>0");
    strcpy((*dict)[3].keyword,"corel_msqd_njump");
    strcpy((*dict)[3].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 4) \write_msqd{#} */
    strcpy((*dict)[4].error_mes,"none");
    strcpy((*dict)[4].keyword,"write_msqd");
    strcpy((*dict)[4].keyarg,"msqd.data");
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_mqsd_output{#} */
    strcpy((*dict)[5].error_mes,"none");
    strcpy((*dict)[5].keyword,"corel_msqd_output");
    strcpy((*dict)[5].keyarg,"atoms");
 /*-----------------------------------------------------------------------*/

 /*========================================================================*/
 } /* end routine set_sim_dict_msqd */
 /*========================================================================*/

void set_sim_dict_iikt_iso(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:
*/
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
  unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 10;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;
 /*========================================================================*/
/*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */

 /* 1) \corel_iikt_iso{#} */
    strcpy((*dict)[1].error_mes,"on or off");
    strcpy((*dict)[1].keyword,"corel_iikt_iso");
    strcpy((*dict)[1].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 2) \corel_iikt_iso_ncor{#} */
    strcpy((*dict)[2].error_mes,"nrun>ncor_ikt_iso>0");
    strcpy((*dict)[2].keyword,"corel_iikt_iso_ncor");
    strcpy((*dict)[2].keyarg,"2");
 /*-----------------------------------------------------------------------*/
 /* 3) \corel_iikt_iso_njump{#} */
    strcpy((*dict)[3].error_mes,"njump_ikt_iso>0");
    strcpy((*dict)[3].keyword,"corel_iikt_iso_njump");
    strcpy((*dict)[3].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 4) \corel_iikt_iso_output_file{#} */
    strcpy((*dict)[4].error_mes,"none");
    strcpy((*dict)[4].keyword,"corel_iikt_iso_output_file");
    strcpy((*dict)[4].keyarg,"iikt_iso.data");
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_iikt_iso_output_kind{#} */
    strcpy((*dict)[5].error_mes,"none");
    strcpy((*dict)[5].keyword,"corel_iikt_iso_output_kind");
    strcpy((*dict)[5].keyarg,"atoms");
 /*-----------------------------------------------------------------------*/
 /* 6) \corel_iikt_iso_kmin{#} */
    strcpy((*dict)[6].error_mes,"a number > 0");
    strcpy((*dict)[6].keyword,"corel_iikt_iso_kmin");
    strcpy((*dict)[6].keyarg,"1.0");
 /*-----------------------------------------------------------------------*/
 /* 7) \corel_iikt_iso_kmax{#} */
    strcpy((*dict)[7].error_mes,"a number > 0");
    strcpy((*dict)[7].keyword,"corel_iikt_iso_kmax");
    strcpy((*dict)[7].keyarg,"5.0");
 /*-----------------------------------------------------------------------*/
 /* 8) \corel_iikt_iso_nb_kvec{#} */
    strcpy((*dict)[8].error_mes,"corel_iikt_iso_nb_kvec>2");
    strcpy((*dict)[8].keyword,"corel_iikt_iso_nb_kvec");
    strcpy((*dict)[8].keyarg,"3");
 /*-----------------------------------------------------------------------*/
 /* 9) \corel_iikt_eisf{#} */
    strcpy((*dict)[9].error_mes,"on or off");
    strcpy((*dict)[9].keyword,"corel_iikt_eisf");
    strcpy((*dict)[9].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 10) \corel_iikt_atoms_list{#} */
    strcpy((*dict)[10].error_mes,"list of atoms");
    strcpy((*dict)[10].keyword,"corel_iikt_atoms_list");
    strcpy((*dict)[10].keyarg,"all");
 /*-----------------------------------------------------------------------*/

 /*========================================================================*/
 } /* end routine set_sim_dict_iikt_iso */
 /*========================================================================*/

void set_sim_dict_ickt_iso(int *num_dict,DICT_WORD *dict[])

/*==========================================================================*/
/*               Begin subprogram:
*/
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
  unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 7;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;
 /*========================================================================*/
/*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */

 /* 1) \corel_ickt_iso{#} */
    strcpy((*dict)[1].error_mes,"on or off");
    strcpy((*dict)[1].keyword,"corel_ickt_iso");
    strcpy((*dict)[1].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 2) \corel_ickt_iso_ncor{#} */
    strcpy((*dict)[2].error_mes,"nrun>ncor_ikt_iso>0");
    strcpy((*dict)[2].keyword,"corel_ickt_iso_ncor");
    strcpy((*dict)[2].keyarg,"2");
 /*-----------------------------------------------------------------------*/
 /* 3) \corel_ickt_iso_njump{#} */
    strcpy((*dict)[3].error_mes,"njump_ikt_iso>0");
    strcpy((*dict)[3].keyword,"corel_ickt_iso_njump");
    strcpy((*dict)[3].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 4) \corel_ickt_iso_output_file{#} */
    strcpy((*dict)[4].error_mes,"none");
    strcpy((*dict)[4].keyword,"corel_ickt_iso_output_file");
    strcpy((*dict)[4].keyarg,"ickt_iso.data");
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_ickt_iso_kmin{#} */
    strcpy((*dict)[5].error_mes,"a number > 0");
    strcpy((*dict)[5].keyword,"corel_ickt_iso_kmin");
    strcpy((*dict)[5].keyarg,"1.0");
 /*-----------------------------------------------------------------------*/
 /* 6) \corel_ickt_iso_kmax{#} */
    strcpy((*dict)[6].error_mes,"a number > 0");
    strcpy((*dict)[6].keyword,"corel_ickt_iso_kmax");
    strcpy((*dict)[6].keyarg,"5.0");
 /*-----------------------------------------------------------------------*/
 /* 7) \corel_ickt_iso_nb_kvec{#} */
    strcpy((*dict)[7].error_mes,"corel_ickt_iso_nb_kvec>2");
    strcpy((*dict)[7].keyword,"corel_ickt_iso_nb_kvec");
    strcpy((*dict)[7].keyarg,"3");
 /*-----------------------------------------------------------------------*/

 /*========================================================================*/
 } /* end routine set_sim_dict_ickt_iso */
 /*========================================================================*/

void set_sim_dict_rdf(int *num_dict,DICT_WORD *dict[])
/*==========================================================================*/
/*               Begin subprogram:
*/
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
  unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 9;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;
 /*========================================================================*/
/*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */
 /*-----------------------------------------------------------------------*/
 /* 1) \corel_rdf{#} */
    strcpy((*dict)[1].error_mes,"on or off");
    strcpy((*dict)[1].keyword,"corel_rdf");
    strcpy((*dict)[1].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 2) \corel_rdf_intra{#} */
    strcpy((*dict)[2].error_mes,"on or off");
    strcpy((*dict)[2].keyword,"corel_rdf_intra");
    strcpy((*dict)[2].keyarg,"on");
 /*-----------------------------------------------------------------------*/
 /* 3) \corel_rdf_njump{#} */
    strcpy((*dict)[3].error_mes,"njump_ikt_iso>0");
    strcpy((*dict)[3].keyword,"corel_rdf_njump");
    strcpy((*dict)[3].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 4) \corel_rdf_npts{#} */
    strcpy((*dict)[4].error_mes,"a integer > 0");
    strcpy((*dict)[4].keyword,"corel_rdf_npts");
    strcpy((*dict)[4].keyarg,"200");
 /*-----------------------------------------------------------------------*/
 /* 5) \corel_rdf_output_file{#} */
    strcpy((*dict)[5].error_mes,"none");
    strcpy((*dict)[5].keyword,"corel_rdf_output_file");
    strcpy((*dict)[5].keyarg,"rdf.data");
 /*-----------------------------------------------------------------------*/
 /* 6) \corel_rdf_periodic_output{#} */
    strcpy((*dict)[6].error_mes,"none");
    strcpy((*dict)[6].keyword,"corel_rdf_periodic_output");
    strcpy((*dict)[6].keyarg,"1");
 /*-----------------------------------------------------------------------*/
 /* 7) \corel_rdf_sq_from_gr{#} */
    strcpy((*dict)[7].error_mes,"on or off");
    strcpy((*dict)[7].keyword,"corel_rdf_sq_from_gr");
    strcpy((*dict)[7].keyarg,"off");
 /*-----------------------------------------------------------------------*/
 /* 8) \corel_rdf_sq_nbk{#} */
    strcpy((*dict)[8].error_mes,"a integer > 0");
    strcpy((*dict)[8].keyword,"corel_rdf_sq_nbk");
    strcpy((*dict)[8].keyarg,"100");
 /*-----------------------------------------------------------------------*/
 /* 9) \corel_rdf_sq_dk{#} */
    strcpy((*dict)[9].error_mes,"a real > 0");
    strcpy((*dict)[9].keyword,"corel_rdf_sq_dk");
    strcpy((*dict)[9].keyarg,"0.1");
 /*-----------------------------------------------------------------------*/


 /*========================================================================*/
 } /* end routine set_sim_dict_rdf */
 /*========================================================================*/


void set_sim_dict_harmonic(int *num_dict,DICT_WORD *dict[])
/*==========================================================================*/
/*               Begin subprogram:
*/
{ /*begin routine*/
 /*=======================================================================*/
 /*             Local variable declarations                                */
  unsigned int i;
 /*========================================================================*/
 /*  0- Malloc the dictionary                                              */

   *num_dict = 2;
   *dict = (DICT_WORD *)cmalloc(*num_dict*sizeof(DICT_WORD))-1;
 /*========================================================================*/
/*========================================================================*/
 /*  I- Initialize the user set option(did the user set the key word       */

   for(i=1;i<=*num_dict;i++){(*dict)[i].iuset = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].iflag = 0;}
   for(i=1;i<=*num_dict;i++){(*dict)[i].key_type = 1;}

 /*========================================================================*/
 /* II- Set up the dictionary                                              */
 /*-----------------------------------------------------------------------*/
 /* 1) \harmonic_frequencies{on,off} */
    strcpy((*dict)[1].error_mes,"on,off");
    strcpy((*dict)[1].keyword,"harmonic_frequencies");
    strcpy((*dict)[1].keyarg,"off");

 /*-----------------------------------------------------------------------*/
 /* 2) \finite_difference_displacement{#} */
    strcpy((*dict)[2].error_mes,"A number > 0 for finite difference displacement");
    strcpy((*dict)[2].keyword,"finite_difference_displacement");
    strcpy((*dict)[2].keyarg,"0.001");

 /*========================================================================*/
 } /* end routine set_sim_dict_harmonic */
 /*========================================================================*/

