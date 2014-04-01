/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: parse                                        */
/*                                                                          */
/* This subprogram reads in all user inputs                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_parse_entry.h"
#include "../proto_defs/proto_parse_local.h"
#include "../proto_defs/proto_sim_params_entry.h"
#include "../proto_defs/proto_mol_params_entry.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_cp_ewald_entry.h"
#include "../proto_defs/proto_surf_params_entry.h"
#include "../proto_defs/proto_inter_params_entry.h"
#include "../proto_defs/proto_inter_params_local.h"
#include "../proto_defs/proto_vps_params_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_scratch_entry.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_entry.h"
#include "../proto_defs/proto_vel_sampl_cp_local.h"
#include "../proto_defs/proto_coords_cp_entry.h"
#include "../proto_defs/proto_coords_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_cpmin_entry.h"
#include "../proto_defs/proto_energy_cp_entry.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  Parse: note there is a noncommuting order of calls                      */
/*==========================================================================*/

void parse(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
           ANALYSIS *analysis,char *input_name)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,cp_on,pimd_on,cp_md,i,pi_beads_proc_st;
  int ivx_flag,ivnhc_flag,ipos_flag;
  int nchrg,iperd,myid,np_states,np_forc,pi_beads,pme_on;
  int myid_state,myid_forc,myatm_start,myatm_end;
  int cp_dual_grid_opt_on;        /*dualed option flag for CP */
  int ncons,ncons_dn;

  CLASS_PARSE      class_parse;
  CP_PARSE         cp_parse;
  FILENAME_PARSE   filename_parse;
  FREE_PARSE       free_parse;
  SPLINE_PARSE     spline_parse;
  NULL_INTER_PARSE null_inter_parse;

  int icontrol_proc  = general_data->error_check_on;
  int num_proc       = class->communicate.np;
  MPI_Comm world     = class->communicate.world;
  double *tot_memory = &(class->tot_memory);


/*========================================================================*/
/*   I) Zero the malloc size variables                                  */

    control_zero_mem(class,bonded,general_data,cp,&class_parse,
                     &null_inter_parse);

/*========================================================================*/
  if(icontrol_proc==1){
  /*========================================================================*/
  /*   II) Set the sim parameters: Done first to get input file names       */
  /*               (interface/sim_params/control_sim_params.c)              */
 
    filename_parse.input_name  = (char *)cmalloc(MAXWORD*sizeof(char));
    strcpy(filename_parse.input_name,input_name);
    control_sim_params(class,general_data,bonded,cp,analysis,
                       &class_parse,&cp_parse,&filename_parse);
  /*========================================================================*/
  /*   III) Read in atom and CP parameters: Done second to get info         */
  /*                                        needed for set_intra;           */
  /*                                        some atom mallocing             */
  /*               (interface/mol_params/control_mol_params.c)              */
  
    control_mol_params(class,general_data,bonded,cp,&class_parse,&cp_parse,
                       &free_parse,&filename_parse);
    
  /*========================================================================*/
  /*  IV) Read in atom, molecule connectivity data: Done before setting     */
  /*                                                 therms;                */
  /*                                                 majority atom mallocing*/
  /*                                                 intramol mallocing     */
  /*                                                 pressure mallocing     */
  /*               (interface/intra_params/control_intra_params.c)          */

    control_intra_params(tot_memory,
                         &(class->clatoms_info),(class->clatoms_pos),
                         &(class->ghost_atoms),&(class->atommaps),
                         bonded,&filename_parse,&free_parse,
                         &class_parse,&null_inter_parse,
                         &(general_data->simopts),&(class->communicate),
                         (class->surface.isurf_on));

  /*========================================================================*/
  }/*endif : icontrol_proc */
/*========================================================================*/
/*    V) Communicate class interface: done before proceeding further      */

 if(num_proc>1){
  Barrier(world);
  communicate_interface(class,bonded,cp,general_data,&null_inter_parse,
                        &class_parse,&cp_parse);
 }/*endif*/

/*========================================================================*/
/*   VI) Assign Flags                                                     */

  cp_dual_grid_opt_on = cp->cpopts.cp_dual_grid_opt;

  pimd_on = general_data->simopts.pimd + general_data->simopts.cp_pimd 
          + general_data->simopts.cp_wave_pimd 
          + general_data->simopts.cp_wave_min_pimd 
          + general_data->simopts.debug_pimd 
          + general_data->simopts.debug_cp_pimd;

  cp_on = general_data->simopts.cp_min  +general_data->simopts.cp_wave_min
         +general_data->simopts.cp 
         +general_data->simopts.cp_wave
         +general_data->simopts.cp_pimd +general_data->simopts.cp_wave_pimd
         +general_data->simopts.debug_cp 
         +general_data->simopts.debug_cp_pimd
         +general_data->simopts.cp_wave_min_pimd;


  cp_md = general_data->simopts.cp   
         +general_data->simopts.cp_wave
         +general_data->simopts.debug_cp
         +general_data->simopts.cp_pimd
         +general_data->simopts.cp_wave_pimd;

  nchrg       = class->clatoms_info.nchrg;
  iperd       = general_data->cell.iperd;
  myid        = class->communicate.myid;
  np_states   = class->communicate.np_states;
  np_forc     = class->communicate.np_forc;
  pi_beads    = class->clatoms_info.pi_beads;
  pme_on      = class->part_mesh.pme_on;

/*========================================================================*/
/*    VII) Read in hmat. Do before set_cp_ewald                           */
/*                (interface/coords/read_coord.c)                         */

  mall_coord(class,general_data);
  mall_pressure(class,general_data);  
  if(myid==0){
    read_hmat(class,general_data,&filename_parse,class_parse.istart,
              cp_dual_grid_opt_on,&(cp->cpewald.dbox_rat),
              &(cp->cpewald.box_rat));
  }/*endif*/
  if(num_proc>1){
    comm_cell_data(&(general_data->cell),
                   &(cp->cpewald.dbox_rat),
                   &(cp->cpewald.box_rat),world);

  }/*endif*/

/*========================================================================*/
/* VIII) Set up the ewald/cp: Done before setting intermol PE             */
/*                            CP/Ewald mallocing                          */
/*                (interface/cp_ewald/control_set_cp_ewald                */

/*--------------------------------------------------------------------------*/
  if((nchrg > 0 && iperd > 0) || cp_on==1){
/*--------------------------------------------------------------------------*/
/* Set up CP and Ewald stuff                                                */

     control_set_cp_ewald(&(general_data->simopts),&(general_data->cell),
                       &(cp->cpcoeffs_info),&(general_data->ewald),
                       &(cp->cpewald),&cp_parse,
                       &(cp->pseudo.gmin_true),
                       &(cp->pseudo.gmin_spl),
                       &(cp->pseudo.gmax_spl),
                       &(class->ewd_scr),(class_parse.kmax_ewd),
                       (class_parse.kmax_res),
                       tot_memory,general_data->timeinfo.int_res_ter,
                       &(class->part_mesh),&(bonded->ecor),myid,
                       cp->cpopts.cp_lsda,
                       general_data->minopts.cp_min_diis,
                       cp_dual_grid_opt_on); 

     class->clatoms_info.alp_ewd = general_data->ewald.alp_ewd;

/*--------------------------------------------------------------------------*/
/*  Calculate Number of CP fictitious degrees of freedom                    */

   if(cp_on) {
       calculate_cp_nfree(cp);

      cp->cpopts.te_ext /= (double) (cp->cpcoeffs_info.cp_nfree);
      cp->vel_samp_cp.div_scal = (double) cp->cpcoeffs_info.cp_nfree;

      if(myid==0){
         printf("\nYour fictious CP temperature is: %gK\n\n",cp->cpopts.te_ext);
      }/*endif*/
   }/* endif cp_on */

/*--------------------------------------------------------------------------*/

 }/*endif nchrg > 0 or cp is on */

/*========================================================================*/
/*   IX) Build communication groups: must be done in serial and parallel  */
/*                                    after call to set_cp_ewald          */

  control_group_communicators(class,cp,cp_on);
  if(num_proc>1){ Barrier(world); }

  myid_state  = class->communicate.myid_state;
  myid_forc   = class->communicate.myid_forc;

/*========================================================================*/
/*   IX) Create the FFT packages  */

  if((nchrg > 0 && iperd > 0 && pme_on==1) || cp_on==1){
    if(myid_state < num_proc && myid_state >= 0){
      control_fft_pkg(&(cp->cp_sclr_fft_pkg3d_sm),
                      &(cp->cp_para_fft_pkg3d_sm),
                      &(cp->cp_sclr_fft_pkg3d_dens_cp_box),
                      &(cp->cp_para_fft_pkg3d_dens_cp_box),
                      &(cp->cp_sclr_fft_pkg3d_lg),
                      &(cp->cp_para_fft_pkg3d_lg),
                      &(general_data->pme_fft_pkg),
                      &(general_data->pme_res_fft_pkg),
                      &(general_data->ewald),&(cp->cpewald),
                      &(class->part_mesh),&(cp->cpcoeffs_info),
                      &(class->communicate),cp_on,cp->cpopts.cp_lsda,
                      tot_memory,general_data->timeinfo.int_res_ter,
                      cp->cpopts.cp_para_opt,cp_dual_grid_opt_on);
    }/* endif */
  }/*endif*/

  if( (cp_on && iperd  == 0) || (iperd == 4)){
    get_coul_clus_corr(&(general_data->ewald),
                       &(class->communicate),
                       &(general_data->cell),(class_parse.kmax_ewd),
                       cp_on,cp->cpcoeffs_info.ecut,
                       cp_parse.cp_ecut_dual_grid,
                       cp->cp_para_fft_pkg3d_lg.igeneric_opt,
                       class->part_mesh.pme_para_opt,tot_memory,
                       cp_dual_grid_opt_on,cp->pseudo.nsplin_g);
  }/*endif*/

  if( iperd  == 2 ){
    get_coul_2D_corr(&(general_data->ewald),&(class->communicate),
                     &(general_data->cell),(cp->pseudo.nsplin_g),cp_on,
                     pme_on,
                     &(cp->cp_para_fft_pkg3d_lg),&(general_data->pme_fft_pkg),
                     tot_memory,1);

  }/*endif*/

  if( iperd  == 1){
     get_coul_1D_corr(&(general_data->ewald),&(class->communicate),
                      &(general_data->cell),(cp->pseudo.nsplin_g),cp_on,
                      pme_on,
                      &(cp->cp_para_fft_pkg3d_lg),&(general_data->pme_fft_pkg),
                      tot_memory,1);
   }/*endif*/

/*========================================================================*/
/*   XI) Setup intermolecular potential stuff: interspline mallocing      */
/*                (interface/inter_params/control_inter_params.c)         */

  control_inter_params(&(class->interact),&spline_parse,
                       &filename_parse,general_data->ewald.alp_ewd,
                       nchrg,class->clatoms_info.natm_tot,
                       class->atommaps.natm_typ,class->atommaps.atm_typ,
                       class->atommaps.iatm_atm_typ,iperd,
                       class_parse.ishift_pot,tot_memory,
                       general_data->timeinfo.int_res_ter,
                       myid,world,num_proc);

/*========================================================================*/
/*    XII) Setup the surface potential if needed                          */

  if(class->surface.isurf_on == 1){
     control_surf_params(&(class->surface), &(filename_parse),
                         class->atommaps.natm_typ, class->atommaps.atm_typ, 
                         tot_memory, myid, world, num_proc); 
  }/* endif */ 

/*========================================================================*/
/*    XII) Setup pseudopotential stuff: pseudospline mallocing            */
/*                (interface/interparams/control_vps_params.c)            */

  if(cp_on==1){
  /*---------------------------------------------------------*/
  /* Create a list of ab initio atoms                        */
  /* List is used in gen_wave and also cp_dual_check routine */

     make_cp_atom_list(&(class->atommaps),
                       class->clatoms_info.cp_atm_flag,
                       &(class->clatoms_info.nab_initio),
                       class->clatoms_info.natm_tot);


   if(myid_state<np_states){
    control_vps_params(&(cp->pseudo),&(general_data->cell),&filename_parse,
                       &spline_parse,class->atommaps.natm_typ,
                       class->atommaps.atm_typ,
                       tot_memory,class->clatoms_info.natm_tot,
                       class->clatoms_info.nab_initio,
                       cp->cpopts.cp_ptens_calc,cp_dual_grid_opt_on,
                       &(class->communicate),cp_parse.cp_ecut);
   }/*endif*/
  }/*endif*/


/*========================================================================*/
/*    XIII) Set particle exclusions and ewald corrections                 */
/*                (interface/lists/set_exclude.c)                         */
  /*THIS HAS ALSO BEEN MODIFIED FOR CLASSICAL HCA EQUILIBRATION */
  set_exclude(&(class->clatoms_info),&(class->ghost_atoms),bonded,
              &(bonded->excl),&null_inter_parse,
              iperd,tot_memory,
              general_data->ewald.alp_ewd,icontrol_proc);

/*========================================================================*/
/*   XV) Set thermostats: Done before reading the coordinates;           */
/*                        atm NHC mallocing                              */
/*                (interface/coords/set_atm_NHC.c                        */

  if( (myid_state==0) || (np_forc > 1) ){
    if((general_data->ensopts.nvt+general_data->ensopts.npt_i
       +general_data->ensopts.npt_f+general_data->ensopts.nvt_isok)==1){
      set_atm_NHC(&(general_data->ensopts),&(general_data->statepoint), 
                  &(general_data->simopts),
                  &(class->clatoms_info),&(class->ghost_atoms),
                  &(class->atommaps),
                  &(class->therm_info_bead),
                  &(class->therm_info_class),
                  (class->therm_bead),&(class->therm_class),
                  &(general_data->baro),&(general_data->par_rahman),
                  &class_parse,
                  &(bonded->constrnt.iconstrnt),tot_memory,pimd_on,
                  icontrol_proc, &(class->communicate),
                  (general_data->cell.hmat_cons_typ));

    }/*endif*/
  }/*endif*/

/*========================================================================*/
/*  XVI) Read in atm positions/velocities/NHCs:                           */
/*                (interface/coords/read_coord.c)                         */

   read_coord(class,general_data,&filename_parse,
              class_parse.istart,cp_dual_grid_opt_on);

/*========================================================================*/
/*  XVII) Spline the ewald corrections (needs particle positions)         */

  if((nchrg > 0 && iperd > 0)){
    splin_ecor(&(bonded->ecor),&(general_data->ewald),(class->clatoms_pos),
               pi_beads,icontrol_proc,tot_memory);
  }else{
    general_data->ewald.self_erf = 1.0;
  }/*endif*/

/*========================================================================*/
/* XVIII) Set thermostats: Done before reading the coeffs                 */
/*                        CP NHC mallocing                                */
/*                (interface/coords/set_coef_NHC.c)                       */

  if(cp_on==1){
   if(myid_state<np_states){
    mall_coef(cp,&(general_data->simopts),class->clatoms_info.pi_beads_proc);
    if(cp_md==1){
     set_coef_NHC(&(cp->cpopts),&(cp->cp_comm_state_pkg_up),
                  &(cp->cp_comm_state_pkg_dn),
                  (cp->cpcoeffs_info.pi_beads_proc),
                  &(cp->cptherm_info),(cp->cptherm_pos),
                  &cp_parse,tot_memory,&(class->communicate));
    } else {
     cp->cptherm_info.num_c_nhc      = 0;
     cp->cptherm_info.num_c_nhc_proc = 0;
     cp->cptherm_info.num_c_nhc_norm = 0;
     cp->cptherm_info.massiv_flag    = 0;
    }/*endif*/
   }/*endif*/
  }/*endif*/

/*========================================================================*/
/*   XIX)malloc scratch space                                            */
/*                (interface/scratch/mall_scratch.c)                     */

  control_mall_scratch(class,bonded,cp,general_data);

/*========================================================================*/
/* XX) Read in coeffs/velocities:                                         */
/*                (interface/coords/read_coef.c)                          */
/*     And tidy up the dual option                                        */


  if(cp_on==1){
   if(myid_state<np_states){
      read_coef(cp,general_data,class,&filename_parse,&cp_parse,tot_memory);
    if(myid == 0){cfree(&(filename_parse.vps_name[1]));} 
    if(cp->cpcoeffs_info.cp_elf_calc_frq > 0) mall_properties(cp);
   }/*endif*/


/* must be done after control_mall_scratch is called  */

   if((cp_dual_grid_opt_on == 1 && np_states >= 1)){
       control_init_dual_maps(&(cp->cpscr),&(general_data->cell),
                              (cp->cpewald.dbox_rat),
                              &(cp->cp_para_fft_pkg3d_dens_cp_box),
                              &(cp->cp_para_fft_pkg3d_lg),
                              cp_dual_grid_opt_on);
   }/*endif cp_dual_grid_opt*/

   if( cp_dual_grid_opt_on == 2){

      init_cp_dual_pme_bw(&(cp->cpscr.cpscr_dual_pme),&(general_data->ewald),
                          &(cp->cpewald),cp->pseudo.n_interp_pme_dual,
                          &(cp->cp_para_fft_pkg3d_lg));


       control_init_dual_maps_pme((cp->cpewald.dbox_rat),
                                  (cp->pseudo.n_interp_pme_dual),
                                   cp_dual_grid_opt_on,
                                  &(cp->cp_para_fft_pkg3d_dens_cp_box),
                                  &(cp->cp_para_fft_pkg3d_lg),
                                  &(cp->cpscr),&(general_data->cell));

    if(pme_on == 1){
      init_cp_atom_pme_bw(&(cp->cpscr.cpscr_atom_pme),&(general_data->ewald),
                          &(cp->cpewald),&(cp->cp_para_fft_pkg3d_lg));
    }/*endif*/

   }/*endif cp_dual_grid_opt*/


   if(cp_dual_grid_opt_on == 1){
     check_box_center(&(general_data->cell),&(cp->cp_para_fft_pkg3d_lg),myid);
   }/*endif cp_dual_grid_opt_on*/

  }/*endif*/
/*========================================================================*/
/* XXI) Set up branch root neighbor list data                             */

  if(class->nbr_list.brnch_root_list_opt>0){
   control_brnch_root_list(class,bonded);
  }/*endif*/

/*========================================================================*/
/* XXI) Control Molecular Decomposition       */

    control_molec_decomp(class,bonded,general_data);

/*========================================================================*/
/*   X) Initialize path integral transformations: after group communicators*/

  if(pi_beads>1||pimd_on==1){
   path_integral_init(&(class->clatoms_info),class->clatoms_pos, 
                      &(class->clatoms_tran),&(class->ghost_atoms), 
                      &(general_data->simopts),&(class->atommaps),
                      &(class->communicate));
  }/*endif*/

/*========================================================================*/
/*  XVII) Communicate classical coordinates to non-bead processors :      */
  /*                                         after path_integral_init     */

   if(np_states>1||np_forc>1){
    ipos_flag = 0;ivx_flag   = 1; ivnhc_flag = 1; /* just send the positions */
     Barrier(world);
     comm_coord_class_state(class,general_data,ivx_flag,ivnhc_flag,ipos_flag,
                            pimd_on);
   }/*endif : np_states*/

/*========================================================================*/
/*   XXII) malloc neigbor list memory                                     */
/*                (interface/lists/mall_make_lists.c)                     */
  
  get_cut_skin(class->interact.cutskin,
               class->interact.cutskin_root,
               class->interact.cutoff,
               class->interact.cutskin_res,
               class->interact.cutskin_root_res,
               class->interact.cutoff_res,
               class->interact.skin,
               class->interact.spread,
               class->interact.brnch_root_skin,
               class->interact.nter_typ);

   pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
   myatm_start = class->clatoms_info.myatm_start;
   myatm_end = class->clatoms_info.myatm_end;

  if(( (myid_state==0) || (np_forc==np_states!= 1)==0  )&&pimd_on==1){ /*DY*/
    control_pimd_trans_mode(class,general_data);

    control_pimd_trans_pos(class,general_data);

  }/*endif*/

  mall_make_lists(class,general_data,bonded,icontrol_proc);

/*========================================================================*/
/* XXIII) Orthogonalize coefficients                                      */

  if(myid_state<np_states){
    control_init_cp_orthog(general_data,cp,&cp_parse,cp_on,cp_md);
  }/*endif*/

/*========================================================================*/
/*  XXIV) Assign/resample the initial class velocities                    */
/*                (interface/vel_sampl)                                   */

  control_init_class_vsmpl(class,&(class_parse),bonded,general_data);

  if(np_states>1||np_forc>1){
    ivx_flag=0;ivnhc_flag=0;ipos_flag=1;  /* send it all */

    if((class_parse.istart <= 2)||(class_parse.ivx_smpl == 1) 
       ||(class_parse.zero_com_vel==1)){
      ivx_flag=1;
    }/*endif*/
    Barrier(world);

    comm_coord_class_state(class,general_data,ivx_flag,ivnhc_flag,ipos_flag,
                            pimd_on);

  }/*endif : np_states*/


/*========================================================================*/
/*  XXV) Assign/resample the initial coef velocities                      */
/*                (interface/vel_sampl)                                   */

  if(cp_md==1){
    if(myid_state<np_states){
     control_init_coef_vsmpl(general_data,cp,&cp_parse);
    }/*endif*/
  }/*endif*/


  /*=============*/
  /* XXVI) timer */
  /*=============*/

  if (myid == 0) {

    PRINT_LINE_STAR;
    printf("Code section timer\n");
    PRINT_LINE_DASH;
    printf("\n");

    #if defined(TIMER)
      printf("Code section timer is enabled.\n");
      #if defined(PARALLEL)
        printf("Parallel version, reporting from MPI rank 0.\n");
      #else
        printf("Serial version.\n");
      #endif
      printf("Printing timer log to standard error output.\n");
      printf("Timer resolution is %2.0e s.\n", timer_resolution());
    #else
      printf("Code section timer was not enabled at compile time.\n");
    #endif
    printf("\n");

    PRINT_LINE_DASH;
    printf("Code section timer - done\n");
    PRINT_LINE_STAR;
    printf("\n");

  }


  /*===============*/
  /* XXVII) PLUMED */

  if (myid == 0) {

    PRINT_LINE_STAR;
    printf("PLUMED interface\n");
    PRINT_LINE_DASH;
    printf("\n");

    #if defined(PLUMED)
      printf("Interface to PLUMED is enabled.\n");
      #if defined(PARALLEL)
        printf("Parallel version.\n");
        printf("Not yet implemented.\n");
        exit(1);
      #else
        printf("Serial version.\n");
      #endif
    #else
      printf("Interface to PLUMED was not enabled at compile time.\n");
    #endif
    printf("\n");

    PRINT_LINE_DASH;
    printf("PLUMED interface - done\n");
    PRINT_LINE_STAR;
    printf("\n");

  }


  /*===========================*/
  /* XXVIII) Flush the buffers */

  fflush(stdout);
  fflush(stderr);


/*------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_zero_mem(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,
                      CP *cp,CLASS_PARSE *class_parse,
                      NULL_INTER_PARSE *null_inter_parse)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/

 class->tot_memory = 0.0;
 zero_sys(class,general_data);
 zero_bnd(bonded);
 zero_cp(cp);
 zero_par(class_parse,null_inter_parse);

/*------------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_init_cp_orthog(GENERAL_DATA *general_data,CP *cp,
                            CP_PARSE *cp_parse,int cp_on,int cp_md)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
#include "../typ_defs/typ_mask.h"

  int iii,ip,is;

/*========================================================================*/
/* I) Find out if you need to orthog                                      */

  iii = cp->cpopts.cp_init_orthog;
  if(general_data->simopts.debug_cp == 1
      && cp->communicate.myid==0) {
   printf("Do you want to orthogonalize the wave functions? (1/0)\n");
   scanf("%d",&iii);
  }/*endif*/
  if(cp->communicate.np>1){Bcast(&iii,1,MPI_INT,0,cp->communicate.world);}

/*========================================================================*/
/* II) Orthogonalize */

  if((cp_on==1) && (iii == 1)){

    for(ip=1;ip<=cp->cpcoeffs_info.pi_beads_proc;ip++){

      cp_shuffle_states(cp,ip);
      orthog_control_cp(cp,ip);

      if((cp_parse->istart_cp > 2)&&(cp_parse->ivc_smpl == 0)&&(cp_md==1)){
        proj_velc(cp,ip);       
      }/*endif:proj_me*/

    }/* endfor:pi_beads */

  }/*endif:ortho*/

/*------------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_init_class_vsmpl(CLASS *class,CLASS_PARSE *class_parse,
                              BONDED *bonded,GENERAL_DATA *general_data)
                   
/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */
  int i;
  int igo;
  int anneal_opt = general_data->simopts.anneal_opt;
  double ann_start_temp = general_data->simopts.ann_start_temp;
  double scale_temp;

/*========================================================================*/
/* I) Particles                                                           */

  if((class_parse->istart <= 2)||(class_parse->ivx_smpl == 1)){
    control_vx_smpl(class,bonded,&(general_data->ensopts),
                    &(general_data->simopts),general_data,
                    general_data->error_check_on);
  }/*endif*/

  if(class_parse->ivx_scale == 1){   
    control_vx_scale(class,&(general_data->simopts),class_parse->text_nhc_mol); 
  }/* endif*/

  igo = 0;
  if( class->clatoms_info.pi_beads > 1 && class->communicate.myid_bead == 0){igo = 1;}
  if( class->clatoms_info.pi_beads == 1 ){igo = 1;}    

  if(class_parse->zero_com_vel==1 && igo == 1){   
     zero_com_vx(class); 
  }/*endif*/
  
  if(class->communicate.myid_bead<class->communicate.np_beads){

/*========================================================================*/
/* II) Extended system                                                    */

   if((general_data->ensopts.nvt+general_data->ensopts.npt_i
     +general_data->ensopts.npt_f+general_data->ensopts.nst+general_data->ensopts.nvt_isok)==1){

     if((class_parse->istart <= 3)||(class_parse->ivnhc_smpl == 1)){
        control_vnhc_smpl(class,general_data);
     }/*endif*/

     if(class_parse->ivnhc_scale == 1){
       scale_temp = (anneal_opt == 1 ? ann_start_temp : general_data->statepoint.t_ext);
       control_vnhc_scale(class,scale_temp);
     }/*endif*/

   }else{

    general_data->baro.v_lnv = 0;
    for(i=1;i<=9;i++){general_data->par_rahman.vgmat[i]=0.0;}

   }/*endif*/

  }/*endif*/

/*------------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/


    

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_init_coef_vsmpl(GENERAL_DATA *general_data,CP *cp,
                            CP_PARSE *cp_parse)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */


/*========================================================================*/
/* I) Coeficients                                                         */

      if((cp_parse->istart_cp <= 2) || (cp_parse->ivc_smpl == 1)){
        control_vc_smpl(general_data,cp);
      }/*endif*/

      if(cp_parse->ivc_scale == 1){
       control_vc_scale(cp);
      }/*endif*/

/*========================================================================*/
/* II) Extended system                                                    */

      if(cp_parse->istate_nhc_opt != 0) {
        if((cp_parse->istart_cp <= 3)||(cp_parse->ivcnhc_smpl == 1)){
         control_vcnhc_smpl(cp);
        }/*endif*/
        if(cp_parse->ivcnhc_scale == 1){
          control_vcnhc_scale(cp);
        }/*endif*/
      }/*endif*/

/*------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_cp_atom_list(ATOMMAPS *atommaps,int *cp_atm_flag,
                       int *nab_initio,int natm_tot)

/*========================================================================*/
/*             Begin Routine                                              */
         {/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

 int i,iii,j;
 int nab_init_cnt = 0;

/*------------------------------------------------------------------------*/

 for(i=1; i <= natm_tot; i++){
   if( cp_atm_flag[i] == 1 ){ nab_init_cnt++; }
 }

   *nab_initio = nab_init_cnt;

/*------------------------------------------------------------------------*/
/* malloc the list array that holds indices of ab initio atoms */

   atommaps->cp_atm_lst = (int *) cmalloc(nab_init_cnt*sizeof(int)) -1;

/*------------------------------------------------------------------------*/
/*  Fill the map with atom inidices of ab initio atoms */

     j = 1;
      
   for(i=1; i <= natm_tot; i++){
     if( cp_atm_flag[i] == 1){
       atommaps->cp_atm_lst[j] = i; j++; 
     }/*endif*/
   }/*endfor*/

/*------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calculate_cp_nfree(CP *cp)

/*========================================================================*/
/*             Begin Routine                                              */
{/*Begin subprogram: */
/*========================================================================*/
/*             Local variable declarations                                */

  int ncons,ncons_dn;

/*------------------------------------------------------------------------*/
/* Number of up degrees of freedom */

    if(cp->cpopts.cp_norb <= 1) {
        ncons = (cp->cpcoeffs_info.nstate_up)*
                     ((cp->cpcoeffs_info.nstate_up)+1)/2;
    }/* endif norb off or norb full_ortho */
    if(cp->cpopts.cp_norb == 2) {ncons = cp->cpcoeffs_info.nstate_up; } /* norm only */
    if(cp->cpopts.cp_norb == 3) {ncons = 0; }                           /* No constraint */
    cp->cpcoeffs_info.cp_nfree = (double) (2*cp->cpcoeffs_info.ncoef*
                                            cp->cpcoeffs_info.nstate_up-ncons);
    cp->cpcoeffs_info.cp_nfree_up = cp->cpcoeffs_info.cp_nfree;

/*------------------------------------------------------------------------*/
/* Number of down degrees of freedom */

  if(cp->cpopts.cp_lsda == 1) {
    if(cp->cpopts.cp_norb <= 1) {
        ncons = (cp->cpcoeffs_info.nstate_up)*
               ((cp->cpcoeffs_info.nstate_up)+1)/2+
                (cp->cpcoeffs_info.nstate_dn)*
               ((cp->cpcoeffs_info.nstate_dn)+1)/2;
        ncons_dn = (cp->cpcoeffs_info.nstate_dn)*
                  ((cp->cpcoeffs_info.nstate_dn)+1)/2;
    }/* endif norb off or norb full_ortho */
    if(cp->cpopts.cp_norb == 2) { 
             ncons = cp->cpcoeffs_info.nstate_up + cp->cpcoeffs_info.nstate_dn;
             ncons_dn = cp->cpcoeffs_info.nstate_dn;
	}/* endif norm only */
    if(cp->cpopts.cp_norb == 3) { 
     	  ncons = 0;
		  ncons_dn = 0;
	}/* endif no constraint */
    cp->cpcoeffs_info.cp_nfree = (double)(2*cp->cpcoeffs_info.ncoef*
                                         (cp->cpcoeffs_info.nstate_up+
                                          cp->cpcoeffs_info.nstate_dn)-ncons);
    cp->cpcoeffs_info.cp_nfree_dn = (double)(2*cp->cpcoeffs_info.ncoef*
                                            cp->cpcoeffs_info.nstate_dn-ncons_dn);
  } /* endif */


/*------------------------------------------------------------------------*/
  }/*end routine*/ 
/*==========================================================================*/
