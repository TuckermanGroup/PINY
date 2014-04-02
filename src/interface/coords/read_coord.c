/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: read_coord                                   */
/*                                                                          */
/* This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        */
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_handle_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void read_coord(CLASS *class,GENERAL_DATA *general_data,
                FILENAME_PARSE *filename_parse,
                int istart,int cp_dual_grid_opt_on)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"

  int iii,upper,lower,igloc,ighost,igo;
  int i,j,ip,nmall,ip_now;       
  int natm_tot_now,istart_now,pi_beads_now;
  int imol_num_now,num_nhc_now,len_nhc_now;
  int itime_dump; 

  double vol,h1,h2,h3,deth;
  double *x_tmp,*y_tmp,*z_tmp;
  double *v_nhc_tmp;

  NAME restart_type_now,atm_typ_now,res_typ_now,mol_typ_now;
  NAME restart_type_spec;
  char *line;

  FILE *fp_dnamei;  

 /* Local Pointers */

  int iseed           = class->vel_samp_class.iseed;
  int iseed2          = class->vel_samp_class.iseed2;
  double qseed        = class->vel_samp_class.qseed;

  char *dnamei        = filename_parse->dnamei;

  int myid            = class->communicate.myid;
  int myid_state      = class->communicate.myid_state;
  int myid_forc       = class->communicate.myid_forc;
  int num_proc        = class->communicate.np;
  int num_proc_forc   = class->communicate.np_forc;
  MPI_Comm world      = class->communicate.world;
  MPI_Comm comm_forc  = class->communicate.comm_forc;

  int ip_start        = class->clatoms_info.pi_beads_proc_st;
  int ip_end          = class->clatoms_info.pi_beads_proc_end;
  int natm_tot        = class->clatoms_info.natm_tot;
  int pi_beads        = class->clatoms_info.pi_beads; 
  int pi_beads_proc   = class->clatoms_info.pi_beads_proc; 

  double *x,*y,*z;
  double *vx,*vy,*vz;

  int class_num_nhc_cast   = class->therm_info_class.num_nhc;
  int class_num_nhc        = class->therm_info_class.num_nhc;
  int class_len_nhc        = class->therm_info_class.len_nhc;
  int bead_len_nhc         = class->therm_info_bead.len_nhc;
  int bead_num_nhc         = class->therm_info_bead.num_nhc;
  int bead_num_nhc_cast    = class->therm_info_bead.num_nhc;

  double **v_nhc,**x_nhc;

  NAME *atm_typ       = class->atommaps.atm_typ;
  int *iatm_atm_typ   = class->atommaps.iatm_atm_typ;
  NAME *res_typ       = class->atommaps.res_typ;
  int *iatm_res_typ   = class->atommaps.iatm_res_typ;
  NAME *mol_typ       = class->atommaps.mol_typ;
  int *iatm_mol_typ   = class->atommaps.iatm_mol_typ;
  int *iatm_mol_num   = class->atommaps.iatm_mol_num;
  int nfreeze         = class->atommaps.nfreeze;
  int pimd_freez_typ  = class->atommaps.pimd_freez_typ;
  int *freeze_map     = class->atommaps.freeze_map;

  int nghost_tot      = class->ghost_atoms.nghost_tot;
  int *ighost_map     = class->ghost_atoms.ighost_map;

  int initial_spread_opt = general_data->simopts.initial_spread_opt;

  int nvt             = general_data->ensopts.nvt;
  int nvt_isok	      = general_data->ensopts.nvt_isok;
  int npt_i           = general_data->ensopts.npt_i;
  int npt_f           = general_data->ensopts.npt_f;

  double *vgmat       = general_data->par_rahman.vgmat;

  double *v_vol_nhc   = general_data->baro.v_vol_nhc;
  double *x_vol_nhc   = general_data->baro.x_vol_nhc;

  int hmat_int_typ    = general_data->cell.hmat_int_typ;
  int hmat_cons_typ   = general_data->cell.hmat_cons_typ;
  int iperd           = general_data->cell.iperd;
  double *hmat        = general_data->cell.hmat;
  int pimd_on;

  pimd_on = (general_data->simopts.pimd + general_data->simopts.cp_pimd 
           + general_data->simopts.cp_wave_pimd 
           + general_data->simopts.cp_wave_min_pimd  
           + general_data->simopts.debug_pimd 
           + general_data->simopts.debug_cp_pimd);

/*========================================================================*/
/* I)Write to screen:                                                   */

  if(myid==0){
    PRINT_LINE_STAR;
    printf("Reading user specified atm coordinate file %s\n",dnamei);
    if(istart==1){printf("using the `initial' restart option\n");}
    if(istart==2){printf("using the `restart_pos' restart option\n");}
    if(istart==3){printf("using the `restart_posvel' restart option\n");}
    if(istart==4){printf("using the `restart_all' restart option\n");}
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*========================================================================*/
/* II) Open the file and malloc:                                          */

  if(myid==0){
    fp_dnamei = cfopen(dnamei,"r");
    line = (char *)cmalloc(MAXLINE*sizeof(char));
  }/*endif*/

  if(num_proc>1){
    Bcast(&class_num_nhc_cast,1,MPI_INT,0,world);
    Bcast(&bead_num_nhc_cast,1,MPI_INT,0,world);
    Bcast(&qseed,1,MPI_DOUBLE,0,world);   /* syncronize random numbers */
    Bcast(&iseed,1,MPI_INT,0,world);
    Bcast(&iseed2,1,MPI_INT,0,world);
  }/*endif*/
  nmall     = MAX(class_num_nhc_cast,bead_num_nhc_cast);

  if( (num_proc_forc>1) && (myid_forc<num_proc_forc) ){
    Bcast(&ip_start,1,MPI_INT,0,comm_forc);
    Bcast(&ip_end,1,MPI_INT,0,comm_forc);
  }/*endif*/

  x_tmp     = (double *)cmalloc(natm_tot*sizeof(double))-1;
  y_tmp     = (double *)cmalloc(natm_tot*sizeof(double))-1;
  z_tmp     = (double *)cmalloc(natm_tot*sizeof(double))-1;
  v_nhc_tmp = (double *)cmalloc(nmall*sizeof(double))-1;
  if(myid_state!=0){
    v_vol_nhc = (double *)cmalloc(class_len_nhc*sizeof(double))-1;
  }/*endif*/

/*========================================================================*/
/*  III)Read in header                                                    */

  if(myid==0){

/*-------------------*/
/* A) Type 1 start : */
    if(istart==1){
      if(fscanf(fp_dnamei,"%d %d %d",&natm_tot_now,&istart_now,
                                   &pi_beads_now)!=3){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error reading start type and number of atoms \n");
        printf("in file %s\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      readtoendofline(fp_dnamei);
    }/*endif*/

/*---------------------*/
/* B) Type 2,3,4 start */
    if(istart>1){
      if(fgetc(fp_dnamei)==0){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error reading header information \n");
        printf("in file %s\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      readtoendofline(fp_dnamei);
      if(fscanf(fp_dnamei,"%d %s %d %d",&natm_tot_now,restart_type_now,
                                        &itime_dump,&pi_beads_now)!=4){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error reading start type and number of atoms \n");
        printf("in file %s\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      }/*endif*/
      readtoendofline(fp_dnamei);

      istart_now = 0;
      if(strcasecmp(restart_type_now,"initial") == 0)        istart_now = 1;
      if(strcasecmp(restart_type_now,"restart_pos") == 0)    istart_now = 2;
      if(strcasecmp(restart_type_now,"restart_posvel") == 0) istart_now = 3;
      if(strcasecmp(restart_type_now,"restart_all") == 0)    istart_now = 4;

      if(istart==1){strcpy(restart_type_spec,"initial");}
      if(istart==2){strcpy(restart_type_spec,"restart_pos");}
      if(istart==3){strcpy(restart_type_spec,"restart_posvel");}
      if(istart==4){strcpy(restart_type_spec,"restart_all");}

      if(istart_now==0){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Start up %s option in ",restart_type_now);
        printf("user specified coordinate file %s\n",dnamei);
        printf("not supported. Supported general_data are: \n");
        printf("initial, restart_pos, restart_posvel, restart_all. \n");     
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      } /* endif */
    } /* endif: istart */
  
/*---------------------*/ 
/* C) General checks   */
    if(istart_now < istart ) {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Start up option, %s, in ",restart_type_now);
      printf("user specified coordinate file %s\n",dnamei);
      printf("Incompatible with class setup,%s\n",restart_type_spec);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    } /* endif */
 
    if(natm_tot_now != natm_tot) {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of particles in\n");
      printf("user specified coordinate file %s \n",dnamei);
      printf("incompatible with class setup\n");
      printf("%d vs %d\n",natm_tot_now,natm_tot);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    } /* endif */
    if(pi_beads_now != pi_beads) {
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of path integral beads in\n");
      printf("user specified coordinate file %s \n",dnamei);
      printf("incompatible with class setup\n");
      printf("%d vs %d\n",pi_beads_now,pi_beads);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    } /* endif */
    if(istart>2 && pi_beads > 1 && initial_spread_opt == 1){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("The spread option is to be used only with \n");
      printf("the options : restart_pos or initial  \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/

  }/*endif:myid==0*/

/*========================================================================*/
/* IV)istart = 1 (initial)                                             */ 

  if(istart == 1) {

/*-----------------------------------------------------------------------*/
/* A)Atm positions                                                  */
    if(myid==0) printf("Reading in coordinates\n");
    upper = pi_beads;
    if(initial_spread_opt == 1){upper = 1;}
    for(ip=1;ip<=upper;ip++){
      if(myid==0){
        for(i=1;i<=natm_tot;i++){
          if(fscanf(fp_dnamei,"%lf %lf %lf",
                        &(x_tmp[i]),
                        &(y_tmp[i]),
                        &(z_tmp[i])) != 3) {
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Error while reading in the %d atom coordinate\n",i);
            printf("in file \"%s\"\n",dnamei);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          }/*endif*/
          readtoendofline(fp_dnamei);
          x_tmp[i] /= BOHR;
          y_tmp[i] /= BOHR;
          z_tmp[i] /= BOHR;
        }/*endfor:atoms*/
      }/*endif*/
      if(num_proc>1){
        Bcast(&(x_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
        Bcast(&(y_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
        Bcast(&(z_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
      }/*endif*/
      if(ip>=ip_start && ip <=ip_end){
         ip_now = ip-ip_start + 1;
         x = class->clatoms_pos[ip_now].x;
         y = class->clatoms_pos[ip_now].y;
         z = class->clatoms_pos[ip_now].z;
         for(i=1;i<=natm_tot;i++){
           x[i] = x_tmp[i];          
           y[i] = y_tmp[i];          
           z[i] = z_tmp[i];          
         }/*endfor*/
      }/*endif*/
    }/*endfor : beads*/

    if(initial_spread_opt == 1 && pi_beads>1){

      spread_coord(&(class->clatoms_info),(class->clatoms_pos),
                   x_tmp,y_tmp,z_tmp,&iseed,&iseed2,&qseed,
                   &(class->atommaps),myid);

      for(ip=1;ip<=pi_beads_proc;ip++){
        get_ghost_pos(&(class->clatoms_info), &(class->clatoms_pos[ip]),
                      &(class->ghost_atoms));
      }/* endfor */
    }/*endif*/
 
/*---------------------------------------------------------------------*/
/* B)Cell shape                                                   */

    if(myid==0){
      for(i=0;i<3;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Error while reading in the %d cell vector \n",i+1);
          printf("in file \"%s\"\n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
        readtoendofline(fp_dnamei);
      }/*endfor*/
    } /* endif */

  }/*endif: start=1 */

/*========================================================================*/
/* V)istart = 2 (restart_pos)                                            */ 

  if(istart >= 2) {

/*----------------------------------------------------------------------*/
/*     A)Atm positions                                                  */
    if(myid==0){
      printf("Reading in coordinates\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("EOF before particle coordinates \n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif : myid*/
    upper = pi_beads;
    if(initial_spread_opt == 1){upper = 1;}
    for(ip=1;ip<=upper;ip++){
      if(myid==0){
        for(i=1;i<=natm_tot;i++){
          if(fscanf(fp_dnamei,"%lf %lf %lf %s %s %s %d",
                            &(x_tmp[i]),
                            &(y_tmp[i]),
                            &(z_tmp[i]),
                            atm_typ_now,res_typ_now,mol_typ_now,
                            &imol_num_now) != 7) {
             printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
             printf("Error while reading in the %d atom coordinate\n",i);
             printf("in file \"%s\"\n",dnamei);
             printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
             fflush(stdout);
             exit(1);
          }/*endif*/
          readtoendofline(fp_dnamei);
          if(strcasecmp(atm_typ_now,atm_typ[iatm_atm_typ[i]]) != 0){
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Atom type mismatch for particle %d\n",i);
            printf("in user specified coordinate file %s \n",dnamei);
            printf("File says %s program expects %s\n",atm_typ_now,
                    atm_typ[iatm_atm_typ[i]]);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          } /* endif */
          if(strcasecmp(res_typ_now,res_typ[iatm_res_typ[i]]) != 0){
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Residue type mismatch for particle %d\n",i);
            printf("in user specified coordinate file %s \n",dnamei);
            printf("File says %s program expects %s \n",res_typ_now,
                    res_typ[iatm_res_typ[i]]);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          } /* endif */
          if(strcasecmp(mol_typ_now,mol_typ[iatm_mol_typ[i]]) != 0){
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Molecule type mismatch for particle %d\n",i);
            printf("in user specified coordinate file %s \n",dnamei);
            printf("File says %s program expects %s \n",mol_typ_now,
                    mol_typ[iatm_mol_typ[i]]);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          } /* endif */
          if(imol_num_now != iatm_mol_num[i]) {
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Molecule number mismatch for particle %d\n",i);
            printf("in user specified coordinate file %s \n",dnamei);
            printf("File says %d program expects %d \n",imol_num_now,
                    iatm_mol_num[i]);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          } /* endif */
        }/*endfor:atoms*/
      }/*endif* myid */
      if(num_proc>1){
        Bcast(&(x_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
        Bcast(&(y_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
        Bcast(&(z_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
      }/*endif*/
      if(ip>=ip_start && ip <=ip_end){
        ip_now = ip - ip_start + 1;
        x = class->clatoms_pos[ip_now].x;
        y = class->clatoms_pos[ip_now].y;
        z = class->clatoms_pos[ip_now].z;
        for(i=1;i<=natm_tot;i++){
          x[i] = x_tmp[i];          
          y[i] = y_tmp[i];          
          z[i] = z_tmp[i];          
        }/*endfor*/
      }/*endif*/
    }/*endfor:pi_beads*/

    if(initial_spread_opt == 1 && pi_beads>1){

      spread_coord(&(class->clatoms_info),(class->clatoms_pos),
                   x_tmp,y_tmp,z_tmp,&iseed,&iseed2,&qseed,
                   &(class->atommaps),myid);

      for(ip=1;ip<=pi_beads_proc;ip++){
        get_ghost_pos(&(class->clatoms_info), &(class->clatoms_pos[ip]),
                      &(class->ghost_atoms));
      }/* endfor */
    }/*endif*/

/*---------------------------------------------------------------------*/
/* B)Cell shape                                                        */

    if(myid==0){
      printf("Reading in cell shape information\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("EOF before cell vectors \n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      for(i=0;i<3;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf",&h1,&h2,&h3) != 3){
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Error reading in cell vector %d\n",i);
          printf("in file \"%s\"\n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
        readtoendofline(fp_dnamei);
      }/*endfor*/
   if(cp_dual_grid_opt_on >= 1){
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
   }/*endif cp_dual_grid_opt*/
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
      readtoendofline(fp_dnamei);
    }/* endif : myid==0 */

  }/*endif: start>=2*/

/*========================================================================*/
/* VI)Istart = 3 (restart_posvel)                                         */

  if(istart >= 3) {

/*----------------------------------------------------------------------*/
/* A)Atm Velocities                                                 */
    if(myid==0){
      printf("Reading in atom velocities\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("EOF before particle velocities \n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endif : myid == 0*/
    for(ip=1;ip<=pi_beads;ip++){
      if(myid==0){
        for(i=1;i<=natm_tot;i++){
          if(fscanf(fp_dnamei,"%lf %lf %lf",
                               &(x_tmp[i]),
                               &(y_tmp[i]),
                               &(z_tmp[i]) )!=3){
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Error while reading in the %d atom velocities\n",i);
            printf("in file \"%s\"\n",dnamei);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          }/*endif*/
          readtoendofline(fp_dnamei);
        }/*endfor:atoms*/
      }/*endif : myid==0*/
      if(num_proc>1){
         Bcast(&(x_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
         Bcast(&(y_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
         Bcast(&(z_tmp[1]),natm_tot,MPI_DOUBLE,0,world);
      }/*endif*/
      if(ip>=ip_start && ip <=ip_end){
         ip_now = ip-ip_start + 1;
         vx = class->clatoms_pos[ip_now].vx;
         vy = class->clatoms_pos[ip_now].vy;
         vz = class->clatoms_pos[ip_now].vz;
         for(i=1;i<=natm_tot;i++){
           vx[i] = x_tmp[i];          
           vy[i] = y_tmp[i];          
           vz[i] = z_tmp[i];          
         }/*endfor*/
      }/*endif*/
    }/*endfor : pi_beads */

/*----------------------------------------------------------------------*/
/* B) Zero velocities of ghost atoms if any                            */

    if(nghost_tot > 0) {
      for(ip=1;ip<=pi_beads_proc;ip++){
        vx = class->clatoms_pos[ip].vx;
        vy = class->clatoms_pos[ip].vy;
        vz = class->clatoms_pos[ip].vz;
        for(ighost=1;ighost <= nghost_tot;ighost++){
          igloc = ighost_map[ighost];
          vx[igloc] = 0.0;
          vy[igloc] = 0.0;
          vz[igloc] = 0.0;
        }/*endfor*/
      }/*endfor : pi_beads_proc*/
    }/*endif*/

/*----------------------------------------------------------------------*/
/* C) Zero velocities of freeze atoms if any                         */

    if(nfreeze > 0) {
      igo = 0;
      upper = (ip_start == 1 ? 1 : 0);
      if(pimd_freez_typ == 2){upper=pi_beads_proc;igo=1;}
      if(igo == 1){
	for(ip=1;ip<=upper;ip++){
        vx = class->clatoms_pos[ip].vx;
        vy = class->clatoms_pos[ip].vy;
        vz = class->clatoms_pos[ip].vz;
        for(i=1;i <= nfreeze;i++){
          igloc = freeze_map[i];
          vx[igloc] = 0.0;
          vy[igloc] = 0.0;
          vz[igloc] = 0.0;
        }/*endfor*/
        }/*endfor*/
      }/*endif*/
    }/*endif*/
  }/* endif : istart >=3*/

/*========================================================================*/
/*   VII)Istart = 4 (restart_all)                                         */

  if(istart == 4) {

/*------------------------------------------------------------------*/
/*  A)Atm NHC Velocities                                         */

    if(myid==0){
      printf("Reading Nose-Hoover chain (NHC) velocities\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("EOF before NHC information \n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
      if(fscanf(fp_dnamei,"%d %d",&num_nhc_now,&len_nhc_now)!=2){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error while reading NHC information\n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/ 
      readtoendofline(fp_dnamei);

      if(num_nhc_now != class_num_nhc) {
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Mismatched number of Nose-Hoover chains\n");
        printf("%d vs %d\n",class_num_nhc, num_nhc_now);
        printf("in user specified coordinate file %s \n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/* endif */
      if(len_nhc_now != class_len_nhc) {
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Mismatched length of Nose-Hoover chains\n");
        printf("%d vs %d\n",class_len_nhc,len_nhc_now);
        printf("in user specified coordinate file %s \n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      } /* endif */
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("EOF before NHC velocities \n");
         printf("in file \"%s\"\n",dnamei);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
    }/*endif : myid==0*/

    v_nhc = class->therm_class.v_nhc;
    for(j=1;j<=class_len_nhc;j++){
      if(myid==0){
        for(i=1;i<=class_num_nhc;i++){
          if(fscanf(fp_dnamei,"%lf",&(v_nhc_tmp[i]))!=1){
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            printf("Error while reading in the %d %d NHC velocity \n",j,i);
            printf("in file \"%s\"\n",dnamei);
            printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
          }/*endif*/
          readtoendofline(fp_dnamei);
        }/*endfor*/
      }/*endif*/
      if(num_proc>1){
        Bcast(&(v_nhc_tmp[1]),class_num_nhc_cast,MPI_DOUBLE,0,world);
      }/*endif*/
      if(ip_start==1){
        for(i=1;i<=class_num_nhc;i++){
          v_nhc[j][i] = v_nhc_tmp[i];
	}/*endfor*/
      }/*endif*/
    }/*endfor : length of nhc */

/*------------------------------------------------------------------*/
/*  B)Bead NHC Velocities                                         */

    if(pi_beads>1){
      if(myid==0){
        printf("Reading Nose-Hoover bead (NHC) velocities\n");
        if(fgets(line,MAXLINE,fp_dnamei)==NULL){
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("EOF before NHC information \n");
          printf("in file \"%s\"\n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
        fgets(line,MAXLINE,fp_dnamei);
        if(sscanf(line,"%d %d \n",&num_nhc_now,&len_nhc_now)!=2){
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Error while reading NHC information\n");
          printf("in file \"%s\"\n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/ 
        if(num_nhc_now != bead_num_nhc) {
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Mismatched number of Nose-Hoover chains\n");
          printf("%d vs %d\n",bead_num_nhc, num_nhc_now);
          printf("in user specified coordinate file %s \n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/* endif */
        if(len_nhc_now != bead_len_nhc) {
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Mismatched length of Nose-Hoover chains\n");
          printf("%d vs %d\n",bead_len_nhc,len_nhc_now);
          printf("in user specified coordinate file %s \n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        } /* endif */
        if(fgets(line,MAXLINE,fp_dnamei)==NULL){
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          printf("EOF before NHC velocities \n");
          printf("in file \"%s\"\n",dnamei);
          printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
      }/*endif : myid==0*/

      for(ip=2;ip<=pi_beads;ip++){
        for(j=1;j<=bead_len_nhc;j++){
          if(myid==0){
            for(i=1;i<=bead_num_nhc;i++){
              if(fgets(line,MAXLINE,fp_dnamei)==NULL){
                printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
                printf("EOF while reading in the %d %d NHC velocity \n",j,i);
                printf("in file \"%s\"\n",dnamei);
                printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
                fflush(stdout);
                exit(1);
              }/*endif*/
              if(sscanf(line,"%lf \n",&v_nhc_tmp[i])!=1){
                printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
                printf("Error reading in the %d %d bead NHC velocity \n",j,i);
                printf("in file \"%s\"\n",dnamei);
                printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
                fflush(stdout);
                exit(1);
              }/*endif*/
            }/*endfor:num_nhc*/
	  }/*endif : myid==0*/
          if(num_proc>1){
            Bcast(&(v_nhc_tmp[1]),bead_num_nhc_cast,MPI_DOUBLE,0,world);
	  }/*endif*/
          if(ip>=ip_start && ip <=ip_end){
            ip_now = ip-ip_start + 1;
            v_nhc  = class->therm_bead[ip_now].v_nhc;
            for(i=1;i<=bead_num_nhc;i++){
              v_nhc[j][i] = v_nhc_tmp[i];          
            }/*endfor*/
          }/*endif : ip in range*/
        }/*endfor:len_nhc*/
      }/*endfor:pi_beads*/

    }/*endif:pi_beads > 1*/

/*------------------------------------------------------------------*/
/*  C)Vol and Vol NHC Velocities                                    */

    if(myid==0){
      printf("Reading Cell velocities\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("EOF before Cell velocities \n");
         printf("in file \"%s\"\n",dnamei);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
      for(i=0;i<3;i++){
        if(fscanf(fp_dnamei,"%lf %lf %lf",
           &(vgmat[(1+i)]),&(vgmat[(4+i)]),&(vgmat[(7+i)])) != 3){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("EOF reading in cell vector velocity %d\n",i);
           printf("in file \"%s\"\n",dnamei);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
        readtoendofline(fp_dnamei);
      }/*endfor*/
      if(hmat_int_typ==0){
        if((vgmat[2]!=vgmat[4])&&
           (vgmat[3]!=vgmat[7])&&
           (vgmat[6]!=vgmat[8])){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("The Cell velocity matrix must be symmetric\n");
           printf("in file \"%s\"\n",dnamei);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
      }/*endif*/
      if(hmat_int_typ==1){
        if((vgmat[2]!=0.0)&&
           (vgmat[3]!=0.0)&&
           (vgmat[6]!=0.0)){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("The Cell velocity matrix must be upper triangular\n");
           printf("in file \"%s\"\n",dnamei);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
      }/*endif*/
     if(hmat_cons_typ==1){
      if( (vgmat[2] != 0.0) || (vgmat[3] != 0.0) ||
          (vgmat[4] != 0.0) || (vgmat[6] != 0.0) || 
          (vgmat[7] != 0.0) || (vgmat[8] != 0.0) ){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The Cell velocity matrix should contain no off-diagonal \n"); 
        printf("coupling with the orthorhombic constraint \n");
        printf("in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);      
      }/* endif */ 
    }/* endif */ 
    if(hmat_cons_typ==2){
      if( (vgmat[3] != 0.0) || (vgmat[6] != 0.0) ||
          (vgmat[7] != 0.0) || (vgmat[8] != 0.0) ){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The Cell velocity matrix should contain no c-coupling with\n");
        printf("the monoclinic constraint in file \"%s\"\n",dnamei);
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);      
      }/* endif */ 
    }/* endif */ 
      if(iperd==2){
        if((vgmat[7]!=0.0)&&
           (vgmat[8]!=0.0)&&
           (vgmat[3]!=0.0)&&
           (vgmat[6]!=0.0)){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("The Cell velocity matrix must not contain c-coupling\n");
           printf("in file \"%s\"\n",dnamei);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
      }/*endif*/
      printf("Reading Volume velocity\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("EOF before Volume velocity\n");
         printf("in file \"%s\"\n",dnamei);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
      if(fscanf(fp_dnamei,"%lf",&(general_data->baro.v_lnv)) != 1){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("Error reading in volume velocity\n");
         printf("in file \"%s\"\n",dnamei);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
      readtoendofline(fp_dnamei);
      printf("Reading Volume NHC velocities\n");
      if(fgets(line,MAXLINE,fp_dnamei)==NULL){
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         printf("EOF before Vol. NHC velocities \n");
         printf("in file \"%s\"\n",dnamei);
         printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
         fflush(stdout);
         exit(1);
      }/*endif*/
      for(i=1;i<=class_len_nhc;i++){
        if(fscanf(fp_dnamei,"%lf",&(v_vol_nhc[i])) != 1){
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           printf("Error reading in volume NHC velocity %d\n",i);
           printf("in file \"%s\"\n",dnamei);
           printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
        readtoendofline(fp_dnamei);
      }/*endfor*/
    } /* endif : myid==0*/

    if(num_proc>1){
      Bcast(&(vgmat[1]),9,MPI_DOUBLE,0,world);
      Bcast(&(v_vol_nhc[1]),class_len_nhc,MPI_DOUBLE,0,world);
      Bcast(&(general_data->baro.v_lnv),1,MPI_DOUBLE,0,world);
    }/*endif*/

  }/*endif : start==4 */

/*========================================================================*/
/*  VIII) Assign NHC positions                      */ 

  if( ( (nvt+npt_i+npt_f)==1 ) && (ip_start ==1) ){
     x_nhc = class->therm_class.x_nhc;
     for(j=1;j<=class_len_nhc;j++){
       for(i=1;i<=class_num_nhc;i++){
          x_nhc[j][i] = 0.0;
       }/*endfor*/
     }/*endfor*/
  }/*endif : bead and ensemble */

  if( pi_beads>1 ){
    lower = (ip_start==1 ? 2 : 1);
    for(ip=lower;ip<=pi_beads_proc;ip++){
      x_nhc = class->therm_bead[ip].x_nhc;
      for(j=1;j<=bead_len_nhc;j++){
        for(i=1;i<=bead_num_nhc;i++){
          x_nhc[j][i] = 0.0;
        }/*endfor:num_nhc*/
      }/*endfor:len_nhc*/
    }/*endfor:pi_beads*/
  }/*endif:pi_beads*/

  if( ((npt_i+npt_f)==1) && (ip_start==1) && (myid_state ==0) ){
    for(j=1;j<=class_len_nhc;j++){
      x_vol_nhc[j] = 0.0;
    }/*endfor*/
  }/*endif*/

/*========================================================================*/
/*  IX) Isokinetic restart_posvel/restart_all */

  if((nvt_isok==1) && (istart>=3)){
	  sampl_isok_restart(&(class->clatoms_info),&(class->clatoms_pos),&(class->therm_info_class),&(class->therm_class));
  }

/*========================================================================*/
/*  X) Calculate the spread                                          */

  class->interact.spread     = 0.0;
  class->interact.spread_now = 0.0;
  if( (pi_beads>1) || (pimd_on==1) ){
    get_pimd_spread(&(class->clatoms_info),class->clatoms_pos,
                    &(class->interact.spread_now),&(class->communicate));
    class->interact.spread = class->interact.spread_now;
  }/*endif*/

/*========================================================================*/
/*  XI) Close file, free arrays                                            */
  
  cfree(&x_tmp[1]);
  cfree(&y_tmp[1]);
  cfree(&z_tmp[1]);
  cfree(&v_nhc_tmp[1]);
  if(myid_state!=0){cfree(&v_vol_nhc[1]);}
  if(myid==0){ 
    fclose(fp_dnamei); 
    cfree(line);
  }/*endif*/

/*========================================================================*/
/*  XII) Output to the screen                                              */

  if(myid==0){

    printf("\n");
    PRINT_LINE_DASH;
    printf("Done reading user specified atm coordinate file %s\n",dnamei);
    PRINT_LINE_STAR;printf("\n");

  }/*endif myid=0*/

/*--------------------------------------------------------------------------*/
    }/* end routine */
/*==========================================================================*/



