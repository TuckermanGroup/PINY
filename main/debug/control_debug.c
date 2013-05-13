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
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_main_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_integrate_md_entry.h"
#include "../proto_defs/proto_integrate_md_local.h"
#include "../proto_defs/proto_output_entry.h"
#include "../proto_defs/proto_output_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_real_space_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_debug(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{   /*begin routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

#include "../typ_defs/typ_mask.h"

  FILE *fp;
  int i,iflag,iii,j,k,m,my_root;
  int iflag_mass = 1;
  int npairs,npairs_tmp,np_tot;       /* Num: # pairs, Max # pairs         */
  int num_proc         = class->communicate.np;
  double cpu1,cpu2,cpu; 
  int error_check_on = general_data->error_check_on;
  int rank           = class->communicate.myid;
  MPI_Comm world     = class->communicate.world;
  MPI_Comm comm_forc = class->communicate.comm_forc;
  double kinet_temp        = 0.0;
  double kinet_nhc_temp    = 0.0;
  double *ptens_tvten     = general_data->ptens.tvten;
  double *ptens_tvten_tmp = general_data->ptens.pvten_tmp;

/*=======================================================================*/
/*   I) Write to Screen                                                  */

  if(rank==0){
   PRINT_LINE_STAR;
   printf("Debugging classical portion of pi_md \n");
   PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*=======================================================================*/
/*  II) Initialize In-output                                             */

  general_data->filenames.ifile_open = 1;

  general_data->timeinfo.itime = 0;
  class->energy_ctrl.itime     = 0;
  
  if(rank==0){
   output_md(class,general_data,bonded);
  }/*endif*/

  general_data->stat_avg.updates   = 0.0;
  general_data->stat_avg.acpu      = 0.0;

/*========================================================================*/
/* III) Initialize free energy stuff                                       */

  if(bonded->bond_free.num != 0){
    for(i=1;i<= (bonded->bond_free.nhist);i++){
      bonded->bond_free.hist[i] = 0.0;
    }
  }/*endif*/
  if(bonded->bend_free.num != 0){
    for(i=1;i<= (bonded->bend_free.nhist);i++){
      bonded->bend_free.hist[i] = 0.0;
    }
  }/*endif*/
  if(bonded->tors_free.num == 1){
    for(i=1;i<= bonded->tors_free.nhist;i++){
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
         bonded->rbar_sig_free.hist_rn[m][k] = 0.0;
     }/*endfor*/
    }/*endfor*/
  }/*endif*/

/*=======================================================================*/
/* IV) Set constraints                                                  */

  if(bonded->constrnt.iconstrnt==1){
     init_constraint(bonded,&(general_data->ptens));
  }/*endif:init constraints, must be done before getting the energy*/

    general_data->stat_avg.iter_shake   = 0;
    general_data->stat_avg.iter_ratl    = 0;
    general_data->stat_avg.iter_23    = 0;
    general_data->stat_avg.iter_33    = 0;
    general_data->stat_avg.iter_46    = 0;

/*=======================================================================*/
/* V) Set energy flags                                                  */

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
/* VI) Test neighbor list update                                         */

  if((class->nbr_list.iver)==1){
    if(rank==0){
      printf("======================================\n");
      if((class->nbr_list.verlist.nolst_ver_update)==1){
         printf("Checking verlist update type: no_lst\n");
      }else{
         printf("Checking verlist update type: lnk_lst\n");
      }/*endif*/ 
      printf("--------------------------------------\n\n");
    }/*endif : rank=0*/
    cputime(&cpu1);
      verlist_control(&(class->clatoms_info),&(class->clatoms_pos[1]),
                 &(class->nbr_list),
                 &(bonded->excl),&(class->atommaps),&(general_data->cell),
                 &(bonded->intra_scr),&(class->for_scr),
                 &(general_data->timeinfo),&(class->interact),error_check_on,
                 &(class->class_comm_forc_pkg));
    cputime(&cpu2);
    cpu = cpu2-cpu1;
    npairs = class->nbr_list.verlist.nver_lst_now;
    np_tot=(class->clatoms_info.natm_tot)*(class->clatoms_info.natm_tot-1)/2
            - (bonded->excl.nlst);
    if(num_proc > 1){
     Reduce(&npairs,&npairs_tmp,1,MPI_INT,MPI_SUM,0,comm_forc);
     npairs = npairs_tmp;
     class->nbr_list.verlist.nver_lst_now = npairs;
    }/*endif*/
    if(rank==0){
      printf("Number of pairs   =%d out of %d; cpu time %g\n",
                                              npairs,np_tot,cpu);
    }
    if(general_data->timeinfo.int_res_ter==1){
        npairs = class->nbr_list.verlist.nver_lst_now_res;
        if(num_proc > 1){
         npairs_tmp=0;
         Reduce(&npairs,&npairs_tmp,1,MPI_INT,MPI_SUM,0,comm_forc);
         npairs = npairs_tmp;
         class->nbr_list.verlist.nver_lst_now_res = npairs;
        }/*endif*/
        if(rank==0){
         printf("Number RESPA pairs=%d out of %d; cpu time %g\n",
                 npairs,np_tot,cpu);
        }/*endif*/
    }
#ifdef DEBUG
   fp = fopen("junk","w");
   for(i=1;i<=class->clatoms_info.natm_tot;i++){
    fprintf(fp,"part %d has %d pairs\n",i,class->nbr_list.verlist.nter[i]);
    if(class->nbr_list.verlist.nter[i]>1){
      small_excl_sort(class->nbr_list.verlist.nter[i],
         &class->nbr_list.verlist.jter[class->nbr_list.verlist.jver_off[i]]);
    }
     for(j=1;j<=class->nbr_list.verlist.nter[i];j++){
      k=class->nbr_list.verlist.jter[j+class->nbr_list.verlist.jver_off[i]];
      fprintf(fp,"%d %d %s %s\n",i,k,
        class->atommaps.atm_typ[class->atommaps.iatm_atm_typ[i]],
        class->atommaps.atm_typ[class->atommaps.iatm_atm_typ[k]]);
     }/*endfor*/
   }/*endfor*/
   fclose(fp);
#endif
#ifdef DEBUG
   for(i=1;i<=5;i++){
     printf("Enter an atom number: ");scanf("%d",&j);
     my_root = j;
     if(class->nbr_list.brnch_root.brnch_atm_map[j]>0){
      my_root = class->nbr_list.brnch_root.brnch_atm_root[
         class->nbr_list.brnch_root.brnch_atm_map[j]];
     }   
      printf("The root atom is %d. It is the %dth root\n",my_root,
              class->nbr_list.brnch_root.root_atm_map[my_root]);
   }   
#endif

    if(rank==0){
      printf("\n--------------------------------------\n");
      printf("Finished checking list update         \n");
      printf("======================================\n\n");
      printf("Enter an integer to continue:  ");
      scanf("%d",&iii); printf("\n");
    }/*endif : rank=0*/
  }/*endif*/

/*=======================================================================*/
/* VI) Test energy                                                        */

  test_energy(class,bonded,general_data);
  get_tvten(&(class->clatoms_info),&(class->clatoms_pos[1]),
            &(general_data->stat_avg),&(general_data->ptens),
            &(general_data->cell));


/*=======================================================================*/
/* VIII)Get NHC Forces                                                     */

  if((general_data->ensopts.nvt)==1){
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
/* VIII) Kinetic Reduces */
 
  if(class->communicate.np_forc>1){

   Allreduce(&(ptens_tvten[1]), &(ptens_tvten_tmp[1]),9,MPI_DOUBLE,
              MPI_SUM,0,comm_forc);
   for(i=1;i<=9;i++){ptens_tvten[i] = ptens_tvten_tmp[i];}
   Reduce(&(general_data->stat_avg.kinet),&kinet_temp,1,
            MPI_DOUBLE,MPI_SUM,0,comm_forc);
   if((general_data->ensopts.nvt + general_data->ensopts.npt_i
      +general_data->ensopts.npt_f) > 0){
    Reduce(&(general_data->stat_avg.kinet_nhc), &kinet_nhc_temp,1,MPI_DOUBLE,
           MPI_SUM,0,comm_forc);
   }/*endif*/
   general_data->stat_avg.kinet          = kinet_temp;
   general_data->stat_avg.kinet_nhc      = kinet_nhc_temp;

  }/*endif*/


/*=======================================================================*/
/*  IX) Write Energy to screen                                           */

  if(rank==0){
    initial_output_md(class,general_data,bonded);
  }/*endif*/


/*=======================================================================*/
/*   X) Write to Screen         */

  if(rank==0){
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed debug\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif : rank=0*/

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/







