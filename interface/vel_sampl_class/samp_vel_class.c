/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/*                                                                   */
/*                         PI_MD:                                    */
/*             The future of simulation technology                   */
/*             ------------------------------------                  */
/*                     Module: samp_vel.c                            */
/*                                                                   */
/* These subprograms sample the velocities                           */
/*                                                                   */
/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_class.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_vel_sampl_class_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* Particle Velocities */
/*===================================================================*/
void sampl_vx(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
              SIMOPTS *simopts,int *iseed,int *iseed2,double *qseed)

{/*begin routine*/
/* ==================================================================*/
/*               Local variable declarations                         */
    int iatm,ip,natm_tot=clatoms_info->natm_tot;
    int pi_beads = clatoms_info->pi_beads; 
    int pi_beads_proc = clatoms_info->pi_beads_proc;
    int anneal_opt = simopts->anneal_opt;
    double ann_start_temp = simopts->ann_start_temp;
    double width,*vx,*vy,*vz,*text_atm,*mass,start_temp;
/*-------------------------------------------------------------------*/
   for(ip=1;ip<=pi_beads_proc;ip++){
     vx = clatoms_pos[ip].vx;
     vy = clatoms_pos[ip].vy;
     vz = clatoms_pos[ip].vz;
     gaussran(natm_tot,iseed,iseed2,qseed,vx);
     gaussran(natm_tot,iseed,iseed2,qseed,vy);
     gaussran(natm_tot,iseed,iseed2,qseed,vz);
   }/*endfor*/
   text_atm = clatoms_info->text_atm;
   if(pi_beads==1){
     mass = clatoms_info->mass;
     for(ip=1;ip<=pi_beads_proc;ip++){
      vx = clatoms_pos[ip].vx;
      vy = clatoms_pos[ip].vy;
      vz = clatoms_pos[ip].vz;
      for(iatm=1;iatm<=natm_tot;iatm++){
        start_temp = (anneal_opt == 1 ? ann_start_temp : text_atm[iatm]);
        width = sqrt(start_temp/(mass[iatm]*BOLTZ));
        vx[iatm] *= width;
        vy[iatm] *= width;
        vz[iatm] *= width;
      } /*endfor*/
     } /*endfor*/
   }else{
     for(ip=1;ip<=pi_beads_proc;ip++){
       vx = clatoms_pos[ip].vx;
       vy = clatoms_pos[ip].vy;
       vz = clatoms_pos[ip].vz;
       mass = clatoms_pos[ip].mass;
       for(iatm=1;iatm<=natm_tot;iatm++){
         start_temp = (anneal_opt == 1 ? ann_start_temp:text_atm[iatm]);
         width = sqrt(start_temp/(mass[iatm]*BOLTZ));
         vx[iatm] *= width;
         vy[iatm] *= width;
         vz[iatm] *= width;
       }/*endfor*/
     }/*endfor*/
   }/*endif*/
}/*end routine*/
/*===================================================================*/

/*===================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*===================================================================*/
/* NHC Velocities */
/*===================================================================*/

void sampl_vnhc(THERM_INFO *therm_info_bead,THERM_POS *therm_bead,   
                THERM_INFO *therm_info_class,THERM_POS *therm_class,
                BARO *baro,PAR_RAHMAN *par_rahman,ENSOPTS *ensopts, 
                STATEPOINT *statepoint, INT_SCR *int_scr, int pi_beads,
                int *iseed, int *iseed2,double *qseed,int iperd,
                int rank,int pi_beads_proc,int hmat_int_typ, 
                int hmat_cons_typ)

/*==================================================================*/
/*            Begin routine*/
      {/*begin routine*/
/*==================================================================*/
/*                Local variable declarations                        */

   double temp[10],vgmat0[10];
   double width,sum;
   int i,j,ip,nnn;

   int npt_i               = ensopts->npt_i;
   int npt_f               = ensopts->npt_f;
   int nvt                 = ensopts->nvt;
   double *vgmat           = par_rahman->vgmat;
   double t_ext            = statepoint->t_ext;
   int class_len_nhc       = therm_info_class->len_nhc;
   int class_num_nhc       = therm_info_class->num_nhc;
   double *class_text_nhc  = therm_info_class->text_nhc;
   double mass_lnv         = baro->mass_lnv;
   double mass_hm          = par_rahman->mass_hm;
   double **class_mass_nhc = therm_info_class->mass_nhc;
   double **class_v_nhc    = therm_class->v_nhc;
   double *atm_kin         = int_scr->atm_kin;
   int bead_len_nhc        = therm_info_bead->len_nhc;
   int bead_num_nhc        = therm_info_bead->num_nhc;
   double *bead_text_nhc   = therm_info_bead->text_nhc;
   double **bead_mass_nhc  = therm_info_bead->mass_nhc;
   double *v_vol_nhc       = baro->v_vol_nhc;
   double *mass_vol_nhc    = baro->mass_vol_nhc;

   double **bead_v_nhc;
   double v_lnv;
  
   int start_proc;
       start_proc = (rank == 0 ? 2:1);

/*==================================================================*/
/* I) Get volume related velocities                                  */
/*------------------------------------------------------------------*/
/* A) NVT */
   if(npt_i==0 && npt_f==0){
     v_lnv= 0.0;
     for(i=1;i<=9;i++) vgmat[i]=0.0;
   }/*endif*/
/*------------------------------------------------------------------*/
/* B) NPTI */
   if(npt_i==1){
     nnn = 1;
     gaussran(nnn,iseed,iseed2,qseed,temp);
     width = sqrt(t_ext/(mass_lnv*BOLTZ));
     v_lnv= temp[1]*width;
     for(i=1;i<=9;i++){ vgmat[i]=0.0;}
     vgmat[1] = v_lnv/3.0;
     vgmat[5] = v_lnv/3.0;
     vgmat[9] = v_lnv/3.0;
   }/*endif*/
/*------------------------------------------------------------------*/
/* C) NPTF */
   if(npt_f==1){
     nnn = 6;
     gaussran(nnn,iseed,iseed2,qseed,temp);
     width = sqrt(t_ext/(mass_hm*BOLTZ));
     vgmat0[1] = temp[1]*width;
     vgmat0[5] = temp[2]*width;
     vgmat0[9] = temp[4]*width;
     vgmat0[2] = temp[3]*width;
     vgmat0[4] = temp[3]*width;
     vgmat0[3] = temp[5]*width;
     vgmat0[7] = temp[5]*width;
     vgmat0[6] = temp[6]*width;
     vgmat0[8] = temp[6]*width;
     constr_cell_mat(iperd,hmat_cons_typ,hmat_int_typ,vgmat0);
     v_lnv = vgmat0[1]+vgmat0[5]+vgmat0[9];
     for(i=1;i<=9;i++) vgmat[i]=vgmat0[i];
   }/*endif*/

   baro->v_lnv  = v_lnv;

/*==================================================================*/
/* II) Get nose-hoover chain velocities                    */ 
/*------------------------------------------------------------------*/
/*    A) Particle NHCs                                              */
   if(nvt==1 || npt_i==1 || npt_f==1) {
     for(i=1;i<=class_len_nhc;i++){
       gaussran(class_num_nhc,iseed,iseed2,qseed,atm_kin);
       for(j=1;j<=class_num_nhc;j++){
         width = sqrt(class_text_nhc[j]/(class_mass_nhc[i][j]*BOLTZ));
         class_v_nhc[i][j] = width*atm_kin[j];
       }/*endfor*/
     }/*endfor*/
   }/*endif*/
/*------------------------------------------------------------------*/
/*    B) Bead NHCs                                                  */
   if(pi_beads>1){
     if(nvt==1 || npt_i==1 || npt_f==1) {
       for(ip=start_proc;ip<=pi_beads_proc;ip++){
         bead_v_nhc = therm_bead[ip].v_nhc;
         for(i=1;i<=bead_len_nhc;i++){
           gaussran(bead_num_nhc,iseed,iseed2,qseed,atm_kin);
           for(j=1;j<=bead_num_nhc;j++){
             width = sqrt(bead_text_nhc[j]/(bead_mass_nhc[i][j]*BOLTZ));
             bead_v_nhc[i][j] = width*atm_kin[j];
           }/*endfor:num_nhc*/
         }/*endfor:len_nhc*/ 
       }/*endfor:pi_beads*/
      }/*endif:ensembles*/
   }/*endif:pi_beads*/
/*------------------------------------------------------------------*/
/*    C) Vol NHCs                                                   */
      nnn = class_len_nhc;
      if(npt_i==1 || npt_f==1){
        gaussran(nnn,iseed,iseed2,qseed,v_vol_nhc);
        for(i=1;i<=class_len_nhc;i++){
          width         = sqrt(t_ext/(mass_vol_nhc[i]*BOLTZ));
          v_vol_nhc[i] *= width;
        }/*endfor */
      } else {
        for(i=1;i<=class_len_nhc;i++){v_vol_nhc[i] = 0.0;}
      }/*endif*/


/*-------------------------------------------------------------------*/
    }/*end routine*/
/*===================================================================*/






