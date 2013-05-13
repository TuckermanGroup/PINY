/*==================================  ===============================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                Module: control_vx_smpl.c                          */
/*                                                                   */
/* Control routine for velocity resampling                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_communicate_wrappers.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vx_smpl(CLASS *class,BONDED *bonded,ENSOPTS *ensopts,
		     SIMOPTS *simopts,GENERAL_DATA *general_data,
                     int error_check_on)

/*===================================================================*/
{/*begin routine*/
/*===================================================================*/
/*    Local variable declarations                                    */

#include "../typ_defs/typ_mask.h"

  int iproj_vel,igloc,ighost,ip,i,ivx_flag,ivnhc_flag,ipos_flag;
  int iup,igo;

  /*-------------------*/
  /*  Local pointers   */
 
  int myid_bead      = class->communicate.myid_bead;
  int np_beads       = class->communicate.np_beads;
  int np_states      = class->communicate.np_states;
  int np_forc        = class->communicate.np_forc;
  int np_tot         = class->communicate.np;
  MPI_Comm world     = class->communicate.world;
  int pi_beads_proc  = class->clatoms_info.pi_beads_proc;
  int pi_beads_proc_st= class->clatoms_info.pi_beads_proc_st;  
  int iconstrnt      = bonded->constrnt.iconstrnt;

  double *vx,*vy,*vz;  
  double temp;
  int nghost_tot     = class->ghost_atoms.nghost_tot;
  int *ighost_map    = class->ghost_atoms.ighost_map;
  int nfreeze        = class->atommaps.nfreeze;
  int pimd_freez_typ = class->atommaps.pimd_freez_typ;
  int *freeze_map    = class->atommaps.freeze_map;

  int npt_f          = ensopts->npt_f;
  int npt_i          = ensopts->npt_i;
  int nve            = ensopts->nve;
  int nvt            = ensopts->nvt;

  int pimd_on,debug_on;
  pimd_on            = general_data->simopts.pimd 
                     + general_data->simopts.cp_pimd 
                     + general_data->simopts.cp_wave_pimd 
                     + general_data->simopts.cp_wave_min_pimd 
                     + general_data->simopts.debug_pimd 
                     + general_data->simopts.debug_cp_pimd;
  debug_on           = simopts->debug + simopts->debug_pimd
                     + simopts->debug_cp + simopts->debug_cp_pimd;

/*===================================================================*/
/* 0) Write to screen                                                */

  if(error_check_on==1){
    PRINT_LINE_STAR;
    printf("Sampling atomic velocities\n");
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*===================================================================*/
/* I) Sample atom velocities                                          */

  if(myid_bead<np_beads){

   sampl_vx(&(class->clatoms_info),(class->clatoms_pos),&(general_data->simopts),
	   &(class->vel_samp_class.iseed),&(class->vel_samp_class.iseed2),
           &(class->vel_samp_class.qseed));

  }/*endif*/

/*====================================================================*/
/* II) Communicate velocities    */

  if( np_states>1 || np_forc>1){
     ivx_flag=0;ivnhc_flag=1;ipos_flag=1;  /* send velocities */
     Barrier(world);
     comm_coord_class_state(class,general_data,ivx_flag,ivnhc_flag,ipos_flag,
                            pimd_on);
  }/*endif : np_states*/

/*====================================================================*/
/* II) Project onto surface of constraint                             */

  if(iconstrnt==1){

    iproj_vel = 1;
    if(debug_on == 1) {
      if(error_check_on==1){
        printf("Do you wish to project velocities? (1 or 0)\n");
        scanf("%d",&iproj_vel);
      }/*endif*/
      if(np_tot>1){
        Bcast(&iproj_vel,1,MPI_INT,0,world);
      }
     }
     if(iproj_vel == 1) {
     if(npt_f==1){ proj_vel_rollf(class,bonded,general_data); }
     if(npt_i==1){ proj_vel_rolli(class,bonded,general_data); }
     if(nve  ==1 || nvt==1){proj_vel(class,bonded,general_data);}
    }/* endif iproj_vel */

  }/*endif iconstrnt*/
 
/*====================================================================*/
/* III) Zero velocities of ghost atoms if any                         */

  if(nghost_tot > 0) {

    for(ip=1;ip <= pi_beads_proc;ip++){
      vx = class->clatoms_pos[ip].vx;
      vy = class->clatoms_pos[ip].vy;
      vz = class->clatoms_pos[ip].vz;
      for(ighost=1;ighost <= nghost_tot;ighost++){
        igloc = ighost_map[ighost];
        vx[igloc] = 0.0;
        vy[igloc] = 0.0;
        vz[igloc] = 0.0;
      }/*endfor*/
    }/*endfor*/

  }/*endif*/

/*====================================================================*/
/* III) Zero velocities of freeze atoms if any                         */


  if(nfreeze > 0) {
    iup = 1; 
    igo=0;

    if(pimd_freez_typ==2){iup=pi_beads_proc;igo=1;}
    if(pimd_freez_typ==1 && pi_beads_proc_st==1){igo=1;}

    if(igo == 1){
      for(ip=1;ip<=iup;ip++){
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

/*====================================================================*/
/* III) Write to screen                                               */

  if(error_check_on==1){
    PRINT_LINE_DASH;
    printf("Atomic velocity sampling complete\n");
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/

/*====================================================================*/
   }/*end routine*/
/*====================================================================*/



/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*====================================================================*/

void zero_com_vx(CLASS *class)

/*===================================================================*/
/*  Begin Routine    */
  {/*begin routine*/
/*===================================================================*/
/*  Local Variables */

#include "../typ_defs/typ_mask.h"

  int i,ip;
  double vxcm,vycm,vzcm,mass_tot;
  double mass_tot_temp,vxcm_temp,vycm_temp,vzcm_temp;

  int pi_beads       = class->clatoms_info.pi_beads;
  int natm_tot       = class->clatoms_info.natm_tot;
  int myatm_start    = class->clatoms_info.myatm_start;
  int myatm_end      = class->clatoms_info.myatm_end;

  int np_forc        = class->class_comm_forc_pkg.num_proc;
  MPI_Comm comm_forc = class->class_comm_forc_pkg.comm;

  double *mass       = class->clatoms_info.mass;
  double *vx         = class->clatoms_pos[1].vx;
  double *vy         = class->clatoms_pos[1].vy;
  double *vz         = class->clatoms_pos[1].vz;
  int *ighost_flag   = class->atommaps.ighost_flag;

  if(pi_beads>1){mass = class->clatoms_pos[1].mass;}

/*===================================================================*/
/* Zero the com velocities                                           */

  mass_tot = 0.0;
  vxcm     = 0.0;
  vycm     = 0.0;
  vzcm     = 0.0;
  for(i=myatm_start;i<=myatm_end;i++){
    if (ighost_flag[i]==0) {
      vxcm     += vx[i]*mass[i];
      vycm     += vy[i]*mass[i];
      vzcm     += vz[i]*mass[i];
      mass_tot += mass[i];
    }
  }/*endfor*/

  if(np_forc > 1){
    vxcm_temp     = vxcm;
    vycm_temp     = vycm;
    vzcm_temp     = vzcm;
    mass_tot_temp = mass_tot;
    Allreduce(&(vxcm_temp),&(vxcm),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
    Allreduce(&(vycm_temp),&(vycm),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
    Allreduce(&(vzcm_temp),&(vzcm),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
    Allreduce(&(mass_tot_temp),&(mass_tot),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
  }/*endif*/
  vxcm /= mass_tot;
  vycm /= mass_tot;
  vzcm /= mass_tot;

  for(i=myatm_start;i<=myatm_end;i++){
    if (ighost_flag[i]==0) {
      vx[i] -= vxcm;
      vy[i] -= vycm;
      vz[i] -= vzcm;
    }
  }/*endfor*/

/*====================================================================*/
    }/*end routine*/
/*====================================================================*/
