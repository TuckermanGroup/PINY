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
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_vel_sampl_class_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#define DEBUG_OFF

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vx_scale(CLASS *class,SIMOPTS *simopts,double *text_mol)

/*===================================================================*/
{/*begin routine*/

#include "../typ_defs/typ_mask.h"

 int ip,imol_typ,iatm,iatm_st,iatm_end;
 double ake,sc,temp_now;
 double *vx,*vy,*vz,*mass;
 double ake_temp;
 int *jatm_jmol_typ_strt  = class->atommaps.jatm_jmol_typ_strt;
 int *natm_1mol_jmol_typ  = class->atommaps.natm_1mol_jmol_typ;
 int *nmol_jmol_typ       = class->atommaps.nmol_jmol_typ;
 int  nmol_typ            = class->atommaps.nmol_typ;
 int *nfree_1mol_jmol_typ = class->atommaps.nfree_1mol_jmol_typ;
 int nmol_have;
 int myatm_start          = class->class_comm_forc_pkg.myatm_start;
 int myatm_end            = class->class_comm_forc_pkg.myatm_end;
 int  pi_beads            = class->clatoms_info.pi_beads;
 int  pi_beads_proc       = class->clatoms_info.pi_beads_proc;
 int  myid                = class->communicate.myid;
 MPI_Comm comm_forc       = class->class_comm_forc_pkg.comm;
 int np_forc              = class->class_comm_forc_pkg.num_proc;
 int anneal_opt           = simopts->anneal_opt;
 double ann_start_temp    = simopts->ann_start_temp;
 double scale_temp;
 int iii;

/*===================================================================*/
/* 0) Write to screen                                                */

if(myid==0){
  PRINT_LINE_STAR;
  printf("Scaling atomic velocities\n");
  PRINT_LINE_DASH;printf("\n");
}/*endif : myid=0*/

/*====================================================================*/
/* I) Scale atom velocities on a molecule type by molecule type basis */

  for(ip=1;ip<=pi_beads_proc;ip++){

    vx = class->clatoms_pos[ip].vx;
    vy = class->clatoms_pos[ip].vy;
    vz = class->clatoms_pos[ip].vz;
    if(pi_beads>1){mass = class->clatoms_pos[ip].mass;}
    if(pi_beads==1){mass = class->clatoms_info.mass;}
    for(imol_typ=1;imol_typ<=nmol_typ;imol_typ++){

      iatm_st = jatm_jmol_typ_strt[imol_typ];
      iatm_end = nmol_jmol_typ[imol_typ]*
                 natm_1mol_jmol_typ[imol_typ] + iatm_st - 1;
      iatm_st   = MAX(iatm_st,myatm_start);
      iatm_end  = MIN(iatm_end,myatm_end);
      nmol_have = (iatm_end-iatm_st+1)/natm_1mol_jmol_typ[imol_typ];
      ake = 0.0;
      for(iatm=iatm_st;iatm<=iatm_end;iatm++){
        ake += (vx[iatm]*vx[iatm] + vy[iatm]*vy[iatm] + vz[iatm]*vz[iatm])
               *mass[iatm];
      }/*endfor:atoms*/
      if(np_forc > 1){
       ake_temp = ake;
       Allreduce(&(ake_temp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
      }/*endif*/
      if(ake>0.0 && nfree_1mol_jmol_typ[imol_typ]>0){
        temp_now = (BOLTZ*ake)/((double)(nfree_1mol_jmol_typ[imol_typ]
                                       *nmol_jmol_typ[imol_typ]));
        scale_temp = (anneal_opt == 1 ? ann_start_temp : text_mol[imol_typ]);
        sc = sqrt(scale_temp/temp_now);
        ake = 0.0;
        for(iatm=iatm_st;iatm<=iatm_end;iatm++){
          vx[iatm] *= sc;
          vy[iatm] *= sc;
          vz[iatm] *= sc;
          ake += (vx[iatm]*vx[iatm] + vy[iatm]*vy[iatm] + vz[iatm]*vz[iatm])
                *mass[iatm];
        }/*endfor:atoms*/
#ifdef DEBUG
        if(np_forc > 1){
          ake_temp = ake;
          Allreduce(&(ake_temp), &(ake),1,MPI_DOUBLE,MPI_SUM,0,comm_forc);
        }/*endif*/
        temp_now = (BOLTZ*ake)/((double)(nfree_1mol_jmol_typ[imol_typ]
                                        *nmol_jmol_typ[imol_typ]));
        printf("temp_now %d %d %g %d\n",ip,imol_typ,temp_now,
                                        class->class_comm_forc_pkg.myid);
        scanf("%d",&iii);
#endif
      }/*endif: something to rescale*/

    }/*endfor:molecules*/

  }/*endfor:beads*/

/*====================================================================*/
/* III) Write to screen                                               */

if(myid==0){
  PRINT_LINE_DASH;
  printf("Atomic velocity scaling complete\n");
  PRINT_LINE_STAR;printf("\n");
}/*endif : myid=0*/

/*====================================================================*/
}/*end routine*/
/*====================================================================*/





/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

void control_vnhc_scale(CLASS *class,double text)

/*===================================================================*/
{/*begin routine*/

#include "../typ_defs/typ_mask.h"

 int ip,inhc,ichain,num_nhc,len_nhc;
 double ake,sc,temp_now;
 double **mass_nhc,**vc_nhc;
 double *text_nhc  = class->therm_info_class.text_nhc;
 int pi_beads      = class->clatoms_info.pi_beads;
 int pi_beads_proc = class->clatoms_info.pi_beads_proc;
 int myid          = class->communicate.myid;

 int start_proc;
 MPI_Comm comm_beads = class->communicate.comm_beads;
      start_proc = (myid == 0 ? 2:1);

/*===================================================================*/
/* 0) Write to screen                                                */

if(myid==0){
  PRINT_LINE_STAR;
  printf("Scaling atomic NHC velocities\n");
  PRINT_LINE_DASH;printf("\n");
}/*endif : myid=0*/

/*====================================================================*/
/* I) Scale atom velocities on a molecule type by molecule type basis */
  ip=1;
  vc_nhc = class->therm_class.v_nhc;
  mass_nhc = class->therm_info_class.mass_nhc;
  num_nhc = class->therm_info_class.num_nhc;
  len_nhc = class->therm_info_class.len_nhc;
  for(inhc=1;inhc<=num_nhc;inhc++){
    ake = 0.0;
    for(ichain=1;ichain<=len_nhc;ichain++){
      ake += mass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                   *vc_nhc[ichain][inhc];     
    }/*endfor : ichain*/
    temp_now = (BOLTZ*ake)/(double)(len_nhc);
    sc = sqrt(text_nhc[inhc]/temp_now);
    ake = 0.0;
    for(ichain=1;ichain<=len_nhc;ichain++){
      vc_nhc[ichain][inhc] *= sc;
      ake += mass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                   *vc_nhc[ichain][inhc];     
    }/*endfor : ichain*/
    temp_now = (BOLTZ*ake)/(double)(len_nhc);
#ifdef DEBUG
      printf("temp_now %d %d %g\n",ip,inhc,temp_now);
#endif
  }/*endfor : inhc*/

  if(pi_beads>1){
    num_nhc = class->therm_info_bead.num_nhc;
    len_nhc = class->therm_info_bead.len_nhc;
    mass_nhc = class->therm_info_bead.mass_nhc;
  }/*endif*/
  for(ip=start_proc;ip<=pi_beads_proc;ip++){
    vc_nhc = class->therm_bead[ip].v_nhc;
    for(inhc=1;inhc<=num_nhc;inhc++){
      ake = 0.0;
      for(ichain=1;ichain<=len_nhc;ichain++){
        ake += mass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                     *vc_nhc[ichain][inhc];     
      }/*endfor : ichain*/
      temp_now = (BOLTZ*ake)/(double)(len_nhc);
      sc = sqrt(text/temp_now);
      ake = 0.0;
      for(ichain=1;ichain<=len_nhc;ichain++){
        vc_nhc[ichain][inhc] *= sc;
        ake += mass_nhc[ichain][inhc]*vc_nhc[ichain][inhc]
                                     *vc_nhc[ichain][inhc];     
      }/*endfor : ichain*/
      temp_now = (BOLTZ*ake)/(double)(len_nhc);
#ifdef DEBUG
      printf("temp_now %d %d %g\n",ip,inhc,temp_now);
#endif
    }/*endfor : inhc*/
  }/*endfor:beads*/

/*====================================================================*/
/* III) Write to screen                                               */

if(myid==0){
  PRINT_LINE_DASH;
  printf("Atomic velocity scaling complete\n");
  PRINT_LINE_STAR;printf("\n");
}/*endif : myid=0*/

/*====================================================================*/
}/*end routine*/
/*====================================================================*/
