/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_md                                   */
/*                                                                          */
/* This subprogram performs MD on a classical potential energy surface (PES)*/
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_cp_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_integrate_cp_entry.h"
#include "../proto_defs/proto_integrate_cp_local.h"
#include "../proto_defs/proto_output_cp_entry.h"
#include "../proto_defs/proto_output_cp_local.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_energy_ctrl_cp_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_debug_cp(CLASS *class,BONDED *bonded,
                              GENERAL_DATA *general_data,CP *cp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
#include "../typ_defs/typ_mask.h"
  int iii,is;
  int i,j,k,m,iflag;        
  int iflag_mass = 1;
  int npairs,np_tot;         /* Num: # pairs, Max # pairs         */
  double cpu1,cpu2,cpu; 
  int error_check_on = general_data->error_check_on;
  double kinet_cp_temp = 0.0;
  double kinet_nhc_cp_temp = 0.0;
  double cp_ehart_temp = 0.0;
  double cp_exc_temp = 0.0;
  double cp_muxc_temp = 0.0;
  double cp_eext_temp = 0.0;
  double cp_enl_temp = 0.0;
  double cp_eke_temp = 0.0;
  int myid_state = class->communicate.myid_state;
  int nstate_up = cp->cpcoeffs_info.nstate_up;
  MPI_Comm comm_states = class->communicate.comm_states;
  int num_proc         = class->communicate.np;

/*=======================================================================*/
/*   I) Write to Screen                                                  */

if(error_check_on==1){
  printf("\n");PRINT_LINE_STAR;
  printf("Debugging CP portion of pi_md \n");
  PRINT_LINE_DASH;printf("\n");
}/*endif*/
if(num_proc>1){Barrier(comm_states);}

/*=======================================================================*/
/* II) Hurl CP Scalars and Rolf CP Vectors */
   
/*=======================================================================*/
/* II) Hurl CP Scalars and Rolf CP Vectors */
   
#ifdef DEBUG
   printf("Do you want to hurl and rolf (1/0)");scanf("%d",&iii);
   if(iii==1){
     hurl_cp_scalars(cp);
     rolf_cp_vectors(cp,class->atommaps.natm_typ,class->clatoms_info.natm_tot);
   }/*endif*/
#endif

/*=======================================================================*/
/* III) Initialize */

  general_data->filenames.ifile_open = 1;

  general_data->timeinfo.itime = 0;
  class->energy_ctrl.itime     = 0;
  
  output_cp(class,general_data,bonded,cp);

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.acpu      = 0.0; 

/*========================================================================*/
/* IV) Initialize free energy stuff                                       */

  if(bonded->bond_free.num != 0){
    for(i=1;i<= (bonded->bond_free.num);i++){
      bonded->bond_free.hist[i] = 0.0;
    }
  }/*endif*/
  if(bonded->bend_free.num != 0){
    for(i=1;i<= (bonded->bend_free.num);i++){
      bonded->bend_free.hist[i] = 0.0;
    }
  }/*endif*/
  if(bonded->tors_free.num == 1){
    for(i=1;i<= bonded->tors_free.num;i++){
      bonded->tors_free.hist[i] = 0.0;
    }
  }/*endif*/
  if(bonded->tors_free.num == 2){
    for(i=1;i<= bonded->tors_free.nhist;i++){
    for(j=1;j<= bonded->tors_free.nhist;j++){
      bonded->tors_free.hist_2d[i][j] = 0.0;
    }/*endfor*/
    }/*endfor*/
  }/*endif*/


  if(bonded->rbar_sig_free.nfree != 0){
    for(i=1;i<=bonded->rbar_sig_free.nhist_sig;i++){
     for(j=1;j<=bonded->rbar_sig_free.nhist_bar;j++){
       bonded->rbar_sig_free.hist[j][i] = 0.0;
     }/*endfor*/
    }/*endfor*/
    for(m=1;m<=bonded->rbar_sig_free.nfree;m++){
     for(k=1;k<=bonded->rbar_sig_free.nhist_bar;k++){
         bonded->rbar_sig_free.hist_rn[m][k]= 0.0;
     }/*endfor*/
    }/*endfor*/
  }/*endif*/

/*=======================================================================*/
/* V) Set constraints                                                  */

  if(bonded->constrnt.iconstrnt==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/

    general_data->stat_avg.iter_shake   = 0;
    general_data->stat_avg.iter_ratl    = 0;
    general_data->stat_avg.iter_23    = 0;
    general_data->stat_avg.iter_33    = 0;
    general_data->stat_avg.iter_46    = 0;

/*=======================================================================*/
/* VI) Set energy flags                                                  */

  class->energy_ctrl.iget_full_inter = 1;
  class->energy_ctrl.iget_res_inter  = 0;
  if(general_data->timeinfo.int_res_ter==1){
    class->energy_ctrl.iget_res_inter=1;
  }
  class->energy_ctrl.iget_full_intra = 1;
  class->energy_ctrl.iget_res_intra  = 0;
  if((general_data->timeinfo.int_res_tra)==1){
    class->energy_ctrl.iget_res_intra=1;
  }
/*=======================================================================*/
/* VII) Test neighbor list update                                         */

  if((class->nbr_list.iver)==1){
    printf("======================================\n");
    if((class->nbr_list.verlist.nolst_ver_update)==1){
       printf("Checking verlist update type: no_lst\n");
    }else{
       printf("Checking verlist update type: lnk_lst\n");
    }/*endif*/ 
    printf("--------------------------------------\n\n");
    cputime(&cpu1);
    verlist_control(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),
                 &(bonded->excl),&(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),&(class->for_scr),
                 &(general_data->timeinfo),&(class->interact),error_check_on,
                 &(class->class_comm_forc_pkg));
    cputime(&cpu2);
    cpu = cpu2-cpu1;
    npairs = class->nbr_list.verlist.nter[(class->clatoms_info.natm_tot)] 
            + class->nbr_list.verlist.jver_off[(class->clatoms_info.natm_tot)];
    np_tot = (class->clatoms_info.natm_tot)*(class->clatoms_info.natm_tot-1)/2 
            - (bonded->excl.nlst);

    printf("Number of pairs = %d out of %d; cpu time %g\n",npairs,np_tot,cpu);
    printf("\n--------------------------------------\n");
    printf("Finished checking list update         \n");
    printf("======================================\n\n");
    printf("Enter an integer to continue:  ");
    scanf("%d",&iii); printf("\n");
  }/*endif*/

/*=======================================================================*/
/* VIII) Test energy                                                        */

  test_energy_cp(class,bonded,general_data,cp);

  get_tvten(&(class->clatoms_info),&(class->clatoms_pos[1]),
            &(general_data->stat_avg),&(general_data->ptens),
            &(general_data->cell));
  get_cpke(&(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
           &(general_data->stat_avg),cp->cpopts.cp_lsda,
           cp->communicate.np_states);
/*=======================================================================*/
/* IX)Get Atm NHC Forces                                                 */

  if((general_data->ensopts.nvt)==1 && myid_state==0){
     init_NHC_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
              &(class->therm_info_class),&(class->therm_class),
              &(class->int_scr),iflag_mass,&(class->class_comm_forc_pkg));
    iflag = 0;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),&(general_data->par_rahman),
                   &(general_data->stat_avg), 
                   &(general_data->statepoint),iflag,
                     class->class_comm_forc_pkg.myid);
  }/*endif*/
  if((general_data->ensopts.npt_i)==1){
     init_NHCPI_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
                &(class->therm_info_class),&(class->therm_class),
                &(general_data->baro),
                &(class->int_scr),&(class->class_comm_forc_pkg));
    iflag = 1;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
                   &(general_data->baro),&(general_data->par_rahman),
                   &(general_data->stat_avg), 
                   &(general_data->statepoint),iflag,
                     class->class_comm_forc_pkg.myid);
  }/*endif*/
  if((general_data->ensopts.npt_f)==1){
     init_NHCPF_par(&(class->clatoms_info),&(class->clatoms_pos[1]),
        &(class->therm_info_class),&(class->therm_class),
        &(general_data->baro),
        &(general_data->par_rahman),&(general_data->cell),&(class->int_scr),
        &(class->class_comm_forc_pkg));
    iflag = 2;
    nhc_vol_potkin(&(class->therm_info_class),&(class->therm_class),
        &(general_data->baro),&(general_data->par_rahman),
        &(general_data->stat_avg), 
        &(general_data->statepoint),iflag,
          class->class_comm_forc_pkg.myid);
  }/*endif*/

  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f)==0){
    general_data->stat_avg.kinet_v       = 0.0;
  }/*endif*/
  if((general_data->ensopts.npt_i+general_data->ensopts.npt_f+
                                  general_data->ensopts.nvt)==0){
    general_data->stat_avg.kinet_nhc     = 0.0;
  }/*endif*/

/*=======================================================================*/
/* X)Get CP NHC Forces                                                   */

  if(cp->cptherm_info.num_c_nhc > 0){
    init_cp_NHC(&(cp->cptherm_info),&(cp->cptherm_pos[1]),
                &(cp->cpcoeffs_info),&(cp->cpcoeffs_pos[1]),
                &(cp->cpscr),cp->cpopts.cp_lsda,comm_states,
                cp->communicate.np_states,cp->communicate.myid_state); 
    nhc_cp_potkin(&(cp->cptherm_info),&(cp->cptherm_pos[1]),
                  &(general_data->stat_avg),myid_state,comm_states);
  }/* endif */


/*=======================================================================*/
/*  Communication for initial output                                    */

if(num_proc>1){

   Reduce(&(general_data->stat_avg.kinet_cp), &kinet_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);
   Reduce(&(general_data->stat_avg.kinet_nhc_cp), &kinet_nhc_cp_temp,1,
                   MPI_DOUBLE,MPI_SUM,0,comm_states);

   general_data->stat_avg.kinet_cp = kinet_cp_temp;
   general_data->stat_avg.kinet_nhc_cp = kinet_nhc_cp_temp;

 }/*endif*/

/*=======================================================================*/
/*  XI) Write Energy to screen                                           */

   if(error_check_on==1){
     initial_output_cp(class,general_data,bonded,cp);
   }/*endif*/
   if(num_proc>1){Barrier(comm_states);}

/*=======================================================================*/
/* XII) Write to Screen         */

if(error_check_on==1){
  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed debug\n");
  PRINT_LINE_STAR;printf("\n");
}/*endif*/
if(num_proc>1){Barrier(comm_states);}

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



#ifdef DEBUG


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void hurl_cp_scalars(CP *cp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */
  int iii;

/*=======================================================================*/
/* 0) OUPUT */

   printf("----------------------------------------------------\n");
   printf("Writing out all CP related scalar quanitities      \n");
   printf("----------------------------------------------------\n\n");

/*=======================================================================*/
/* I) CPOPTS */

   printf("CPOPTS: %d %d %d %d %d %d %d %g \n",
           cp->cpopts.cp_lda,
           cp->cpopts.cp_lsda,
           cp->cpopts.cp_sic,
           cp->cpopts.cp_gga,
           cp->cpopts.cp_nonint,
           cp->cpopts.cp_norb,
           cp->cpopts.cp_gauss,
           cp->cpopts.te_ext); mal_verify(1);

/*=======================================================================*/
/* II) CPCOEFFS */

   printf("CPCOEFFS: %d %d %d %d %d %d %d \n",
           cp->cpcoeffs_info.ncoef,
           cp->cpcoeffs_info.ncoef_l,
           cp->cpcoeffs_info.nstate_up,
           cp->cpcoeffs_info.nstate_dn,
           cp->cpcoeffs_info.icmass_unif,
           cp->cpcoeffs_info.ks_rot_on,
           cp->cpcoeffs_info.n_ks_rot); mal_verify(1);

/*=======================================================================*/
/* III) CPCOEFFS */

   printf("CPTHERM: %d %d %d %d \n",   
           cp->cptherm_info.num_c_nhc,
           cp->cptherm_info.len_c_nhc,
           cp->cptherm_info.nres_c_nhc,
           cp->cptherm_info.nyosh_c_nhc); mal_verify(1);

/*=======================================================================*/
/* IV) CPCONSTRNT */

   printf("CPCONSTRNT: %g %g %g \n",   
           cp->cpconstrnt.c_tolshake,
           cp->cpconstrnt.c_tolratl,
           cp->cpconstrnt.c_tolnorb); mal_verify(1);

/*=======================================================================*/
/* V) CPCEWALD */

   printf("CPEWALD: %d %d %d %d \n",
          cp->cpewald.kmax_cp[1],
          cp->cpewald.kmax_cp[2],
          cp->cpewald.kmax_cp[3],
          cp->cpewald.nktot_sm); mal_verify(1);

/*=======================================================================*/
/* VI) PSEUDO */


   printf("PSEUDO: %d %d %d %g %g %g %g %d %d %d %s %s %s\n",
           cp->pseudo.n_ang_max,
           cp->pseudo.nsplin_g, 
           cp->pseudo.nsplin_g_tot,
           cp->pseudo.gmin_true,
           cp->pseudo.gmin_spl,
           cp->pseudo.gmax_spl,
           cp->pseudo.dg_spl,
           cp->pseudo.nvps_lab,
           cp->pseudo.num_nl_lst,
           cp->pseudo.nl_cut_on,
           cp->pseudo.vxc_typ,
           cp->pseudo.ggax_typ,
           cp->pseudo.ggac_typ); mal_verify(1);

/*=======================================================================*/
/* VII) DONE */

   printf("\n----------------------------------------------------\n");
   printf("Completed CP scalar debug                           \n");
   printf("----------------------------------------------------\n\n");

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void rolf_cp_vectors(CP *cp,int natm_typ,int natm_tot)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int imax_up,imax_dn,imax,iii;
  int ioff;
  double g;

/*=======================================================================*/
/* 0) OUPUT */

   printf("----------------------------------------------------\n");
   printf("Writing out selected elements of the CP vector stuff\n");
   printf("----------------------------------------------------\n\n");

/*=======================================================================*/
/* I) CPCOEFFS */

   printf("CPCOEFFS FORCE\n");
   imax_up = cp->cpcoeffs_info.ncoef*cp->cpcoeffs_info.nstate_up;
   cp->cpcoeffs_pos[1].fcre_up[1] = 0.;
   cp->cpcoeffs_pos[1].fcim_up[1] = 0.;
   cp->cpcoeffs_pos[1].fcre_up[imax_up] = 0.;
   cp->cpcoeffs_pos[1].fcim_up[imax_up] = 0.;mal_verify(1);
   printf("CPCOEFFS MASS: %g %g %g\n",
           cp->cpcoeffs_info.cmass[1],
           cp->cpcoeffs_info.cmass[2],
           cp->cpcoeffs_info.cmass[(cp->cpcoeffs_info.ncoef)]);
           mal_verify(1);
   printf("CPCOEFFS UP: %g %g %g %g %g %g %g %g \n",
   cp->cpcoeffs_pos[1].cre_up[1],
   cp->cpcoeffs_pos[1].cim_up[1],
   cp->cpcoeffs_pos[1].cre_up[imax_up],
   cp->cpcoeffs_pos[1].cim_up[imax_up],
   cp->cpcoeffs_pos[1].vcre_up[1],
   cp->cpcoeffs_pos[1].vcim_up[1],
   cp->cpcoeffs_pos[1].vcre_up[imax_up],
   cp->cpcoeffs_pos[1].vcim_up[imax_up]); mal_verify(1);
   if(cp->cpopts.cp_lsda==1){
      printf("CPCOEFFS FORCE\n");
      imax_dn = cp->cpcoeffs_info.ncoef*cp->cpcoeffs_info.nstate_dn;
      cp->cpcoeffs_pos[1].fcre_dn[1] = 0.;
      cp->cpcoeffs_pos[1].fcim_dn[1] = 0.;
      cp->cpcoeffs_pos[1].fcre_dn[imax_dn] = 0.;
      cp->cpcoeffs_pos[1].fcim_dn[imax_dn] = 0.;
      printf("CPCOEFFS DN: %g %g %g %g %g %g %g %g \n",
      cp->cpcoeffs_pos[1].cre_dn[1],
      cp->cpcoeffs_pos[1].cim_dn[1],
      cp->cpcoeffs_pos[1].cre_dn[imax_dn],
      cp->cpcoeffs_pos[1].cim_dn[imax_dn],
      cp->cpcoeffs_pos[1].vcre_dn[1],
      cp->cpcoeffs_pos[1].vcim_dn[1],
      cp->cpcoeffs_pos[1].vcre_dn[imax_dn],
      cp->cpcoeffs_pos[1].vcim_dn[imax_dn]); mal_verify(1);
   }/*endif*/

/*=======================================================================*/
/* II) CPCTHERM */

   if(cp->cptherm_info.num_c_nhc>0){
     imax_up = cp->cpcoeffs_info.ncoef*cp->cpcoeffs_info.nstate_up;
     imax_dn = cp->cpcoeffs_info.ncoef*cp->cpcoeffs_info.nstate_dn;
     printf("CPTHERM FORCE\n");
     cp->cptherm_pos[1].fc_nhc[1][1] = 0.0;
     cp->cptherm_pos[1].
         fc_nhc[cp->cptherm_info.len_c_nhc][cp->cptherm_info.num_c_nhc] = 0.0;
     printf("CPTHERM: %g %g %g %g %g %g %g %g \n",
          cp->cptherm_pos[1].c_nhc[1][1],
          cp->cptherm_pos[1].
              c_nhc[cp->cptherm_info.len_c_nhc][cp->cptherm_info.num_c_nhc],
          cp->cptherm_pos[1].vc_nhc[1][1],
          cp->cptherm_pos[1].
              vc_nhc[cp->cptherm_info.len_c_nhc][cp->cptherm_info.num_c_nhc],
          cp->cptherm_info.c_gkt[1][1],
          cp->cptherm_info.
              c_gkt[cp->cptherm_info.len_c_nhc][cp->cptherm_info.num_c_nhc],
          cp->cptherm_info.cmass_nhc[1][1],
          cp->cptherm_info.
          cmass_nhc[cp->cptherm_info.len_c_nhc][cp->cptherm_info.num_c_nhc]);
     printf("CPTHERM MAP UP: %d %d \n",
             cp->cptherm_info.icmapup_nhc[1],
             cp->cptherm_info.icmapup_nhc[imax_up]);mal_verify(1);
     if(cp->cpopts.cp_lsda==1){
       printf("CPTHERM MAP DN: %d %d \n",
             cp->cptherm_info.icmapdn_nhc[1],
             cp->cptherm_info.icmapdn_nhc[imax_dn]);mal_verify(1);
        mal_verify(1);
     }/*endif*/
   }/*endif*/

/*=======================================================================*/
/* III) CPCEWALD */

   imax = cp->cpewald.nktot_sm;
   printf("CPEWALD SM MAP: %d %d %d %d %d %d %d %d \n",
          cp->cpewald.map_cpr[1],
          cp->cpewald.map_cpi[1],
          cp->cpewald.map_cprc[1],
          cp->cpewald.map_cpic[1],
          cp->cpewald.map_cpr[imax],
          cp->cpewald.map_cpi[imax],
          cp->cpewald.map_cprc[imax],
          cp->cpewald.map_cpic[imax]);mal_verify(1);
   printf("CPEWALD SM: %d %d %d %d %d %d %d %d %d %d\n",
          cp->cpewald.kastr_sm[1],
          cp->cpewald.kbstr_sm[1],
          cp->cpewald.kcstr_sm[1],
          cp->cpewald.kastr_sm[imax],
          cp->cpewald.kbstr_sm[imax],
          cp->cpewald.kcstr_sm[imax],
          cp->cpewald.ibrk1_sm[1],
          cp->cpewald.ibrk2_sm[1],
          cp->cpewald.ibrk1_sm[imax],
          cp->cpewald.ibrk2_sm[imax]);mal_verify(1);
   imax = cp->cpcoeffs_info.ncoef_l-1;
   printf("CPEWALD LG MAP: %d %d %d %d %d %d %d %d \n",
          cp->cpewald.map_cpr_l[1],
          cp->cpewald.map_cpi_l[1],
          cp->cpewald.map_cprc_l[1],
          cp->cpewald.map_cpic_l[1],
          cp->cpewald.map_cpr_l[imax],
          cp->cpewald.map_cpi_l[imax],
          cp->cpewald.map_cprc_l[imax],
          cp->cpewald.map_cpic_l[imax]);mal_verify(1);

/*=======================================================================*/
/* IV) PSEUDO */

  printf("PSEUDO INFO: %d %d %d %d %d %d %g %g\n",
          cp->pseudo.n_ang[1],
          cp->pseudo.n_ang[natm_typ],
          cp->pseudo.loc_opt[1],
          cp->pseudo.loc_opt[natm_typ],
          cp->pseudo.ivps_label[1],
          cp->pseudo.ivps_label[natm_typ],
          cp->pseudo.rcut_nl[1],
          cp->pseudo.rcut_nl[natm_typ]);mal_verify(1);
  imax = natm_typ*(cp->pseudo.n_ang_max+1);
  printf("PSEUDO DATA: %g %g %g %g %g %g \n",
          cp->pseudo.gzvps[1],
          cp->pseudo.gzvps[natm_typ],
          cp->pseudo.gzvps0[1],
          cp->pseudo.gzvps0[natm_typ],
          cp->pseudo.vpsnorm[1],
          cp->pseudo.vpsnorm[imax]);mal_verify(1);

  printf("PSEUDO LIST\n");
  cp->pseudo.np_nl[1] = 0; cp->pseudo.np_nl[(cp->pseudo.n_ang_max+1)] = 0;
  imax = natm_tot*(cp->pseudo.n_ang_max+1);
  cp->pseudo.ip_nl[1] = 0; cp->pseudo.ip_nl[imax]=0;
  cp->pseudo.x0w[1] = 0.0; cp->pseudo.x0w[natm_tot] = 0.0;
  cp->pseudo.y0w[1] = 0.0; cp->pseudo.y0w[natm_tot] = 0.0;
  cp->pseudo.z0w[1] = 0.0; cp->pseudo.z0w[natm_tot] = 0.0; 
  printf("PSEUDO SPL: %g %g %g %g %g %g %g %g\n",
          cp->pseudo.vps0[1],
          cp->pseudo.vps0[cp->pseudo.nsplin_g_tot],
          cp->pseudo.vps1[1],
          cp->pseudo.vps1[cp->pseudo.nsplin_g_tot],
          cp->pseudo.vps2[1],
          cp->pseudo.vps2[cp->pseudo.nsplin_g_tot],
          cp->pseudo.vps3[1],
          cp->pseudo.vps3[cp->pseudo.nsplin_g_tot]);
   mal_verify(1);

   printf("\n----------------------------------------------------\n");
   printf("Completed CP vector debug                             \n");
   printf("----------------------------------------------------\n\n");

/*=======================================================================*/
/* V) Test VPS stuff */

   printf("----------------------------------------------------\n");
   printf("Selected elements of the non-local pseudopotential \n");
   printf("----------------------------------------------------\n\n");

   printf("   l=0\n");
   printf("================\n");
   g = cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[1]);
   g = cp->pseudo.dg_spl*((double) (cp->pseudo.nsplin_g-1))
     + cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[cp->pseudo.nsplin_g]);

   printf("   l=1\n");
   printf("================\n");
   ioff = cp->pseudo.nsplin_g;
   g = cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[(1+ioff)]);
   g = cp->pseudo.dg_spl*((double) (cp->pseudo.nsplin_g-1))
     + cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[cp->pseudo.nsplin_g+ioff]);

   printf("   l=2\n");
   printf("================\n");
   ioff = cp->pseudo.nsplin_g*2;
   g = cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[(1+ioff)]);
   g = cp->pseudo.dg_spl*((double) (cp->pseudo.nsplin_g-1))
     + cp->pseudo.gmin_spl;
   printf("g=%g,vps0[g]=%g\n",g,cp->pseudo.vps0[cp->pseudo.nsplin_g+ioff]);

   printf("\n----------------------------------------------------\n");
   printf("Completed non-local pseudopotential debug             \n");
   printf("----------------------------------------------------\n");

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



#endif





