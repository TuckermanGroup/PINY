/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: int_utilities                                */
/*                                                                          */
/* This subprogram provides some integrator utility routines                */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_integrate_pimd_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_tvten_pimd(CLASS *class,GENERAL_DATA *general_data)
/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */
    int i,ipart,ip,iii;
    double akinet;
    double *ptens_tvten = general_data->ptens.tvten;
    double *clatoms_vx;
    double *clatoms_vy;
    double *clatoms_vz;
    double *clatoms_mass;
    int natm_tot         = class->clatoms_info.natm_tot;
    int myid_bead        = class->communicate.myid_bead;
    int pi_beads         = class->clatoms_info.pi_beads;
    int pi_beads_proc    = class->clatoms_info.pi_beads_proc;
    int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
    int myatm_start      = class->clatoms_info.myatm_start;
    int myatm_end        = class->clatoms_info.myatm_end;

/*========================================================================*/
           /* Tvten Matrix */
           /*  1:  2:  3   */
           /*  4:  5:  6   */
           /*  7:  8:  9   */
/*=======================================================================*/
/* II) Accumulate Tvten */

  for(i=1;i<=9;i++){
    (ptens_tvten)[i] = 0.0;
  /*endfor*/}

  if(pi_beads_proc_st==1){
    ip=1;
    clatoms_vx   = class->clatoms_pos[ip].vx;
    clatoms_vy   = class->clatoms_pos[ip].vy;
    clatoms_vz   = class->clatoms_pos[ip].vz;
    clatoms_mass = class->clatoms_pos[ip].mass;
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       ptens_tvten[1] += (clatoms_mass[ipart]*
                            clatoms_vx[ipart]*clatoms_vx[ipart]);
       ptens_tvten[5] += (clatoms_mass[ipart]*
                            clatoms_vy[ipart]*clatoms_vy[ipart]);
       ptens_tvten[9] += (clatoms_mass[ipart]*
                            clatoms_vz[ipart]*clatoms_vz[ipart]);
       ptens_tvten[2] += (clatoms_mass[ipart]*
                            clatoms_vx[ipart]*clatoms_vy[ipart]);
    /*endfor*/}
    (ptens_tvten)[4] = (ptens_tvten)[2];
    if((general_data->cell.iperd)==3){
        ip=1;
        clatoms_vx   = class->clatoms_pos[ip].vx;
        clatoms_vy   = class->clatoms_pos[ip].vy;
        clatoms_vz   = class->clatoms_pos[ip].vz;
        clatoms_mass = class->clatoms_pos[ip].mass;
        for(ipart=myatm_start;ipart<=myatm_end;ipart++){
          ptens_tvten[3] += (clatoms_mass[ipart]*
                               clatoms_vx[ipart]*clatoms_vz[ipart]);
          ptens_tvten[6] += (clatoms_mass[ipart]*
                               clatoms_vy[ipart]*clatoms_vz[ipart]);
        /*endfor*/}
      ptens_tvten[7] = ptens_tvten[3];
      ptens_tvten[8] = ptens_tvten[6];
    }/*endif : iperd*/
  }/*endif:myid_bead==0*/

/*=======================================================================*/
/* III) Kinetic Energy                                                   */

  akinet = 0.0;
  for(ip=1;ip<=pi_beads_proc;ip++){
    clatoms_vx   = class->clatoms_pos[ip].vx;
    clatoms_vy   = class->clatoms_pos[ip].vy;
    clatoms_vz   = class->clatoms_pos[ip].vz;
    clatoms_mass = class->clatoms_pos[ip].mass;
    for(ipart=myatm_start;ipart<=myatm_end;ipart++){
       akinet += clatoms_mass[ipart]*(
                            clatoms_vx[ipart]*clatoms_vx[ipart]
              +             clatoms_vy[ipart]*clatoms_vy[ipart]
              +             clatoms_vz[ipart]*clatoms_vz[ipart]);
      /*endfor*/}
  /*endfor*/}
  (general_data->stat_avg.kinet) = akinet/2.0;

/*-------------------------------------------------------------------------*/
/*end routine*/}
/*==========================================================================*/





/*==========================================================================*/
/*==========================================================================*/
void nhc_vol_potkin_pimd(CLASS *class, GENERAL_DATA *general_data,int iflag)

/*========================================================================*/
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

int i,ip,inhc,ichain,iii;
int len_nhc,num_nhc;
int pi_beads = class->clatoms_info.pi_beads;
int pi_beads_proc = class->clatoms_info.pi_beads_proc;
int rank = class->communicate.myid_bead;
int mytherm_start     = class->therm_info_class.mytherm_start;
int mytherm_end       = class->therm_info_class.mytherm_end;
int pi_beads_proc_st  = class->clatoms_info.pi_beads_proc_st;
int np_forc           = class->communicate.np_forc;
int myid_forc         = class->communicate.myid_forc;
int myid              = class->communicate.myid;
int start,num_nhc_proc;
double kinet_nhc,kinet_nhc_old,kinet_nhc_loop,vpotnhc,temp1,temp2;
double **therm_x_nhc,x_nhc_tot;
double **therm_v_nhc;
double **therm_mass_nhc;
double **therm_gkt,gkt;

/*==========================================================================*/
/* I) nhc stuff */

    vpotnhc   = 0.0;
    kinet_nhc = 0.0;
    kinet_nhc_old = 0.0;
   

  start = (pi_beads_proc_st == 1 ? 2:1);
  if(pi_beads_proc_st == 1){   
    len_nhc   = class->therm_info_class.len_nhc;
    num_nhc   = class->therm_info_class.num_nhc;
    therm_x_nhc = class->therm_class.x_nhc;
    therm_v_nhc = class->therm_class.v_nhc;
    therm_mass_nhc = class->therm_info_class.mass_nhc;
    therm_gkt = class->therm_info_class.gkt;
    for(ichain=1;ichain<=len_nhc;ichain++){
      for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
        kinet_nhc += (therm_mass_nhc[ichain][inhc]
               *therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]);
        vpotnhc += (therm_gkt[ichain][inhc]*
                              therm_x_nhc[ichain][inhc]);   
    /*endfor*/}
      if(iflag>0&&myid_forc==0){
        kinet_nhc += (general_data->baro.mass_vol_nhc[ichain]
               *general_data->baro.v_vol_nhc[ichain]
               *general_data->baro.v_vol_nhc[ichain]);
        vpotnhc   += (general_data->baro.gkt_vol[ichain]*
                              general_data->baro.x_vol_nhc[ichain]);
      /*endif*/}
     /*endfor*/}
  }/* endif rank */

    general_data->stat_avg.kinet_nhc = kinet_nhc/2.0;

    mytherm_start     = class->therm_info_bead.mytherm_start;
    mytherm_end       = class->therm_info_bead.mytherm_end;
    kinet_nhc = 0.0;
    len_nhc   = class->therm_info_bead.len_nhc;
    num_nhc   = class->therm_info_bead.num_nhc;
    therm_mass_nhc = class->therm_info_bead.mass_nhc;
    therm_gkt = class->therm_info_bead.gkt;
    gkt = class->therm_info_bead.gkt[1][1];
    num_nhc_proc = mytherm_end - mytherm_start + 1;
    for(ip=start;ip<=pi_beads_proc;ip++){
     x_nhc_tot   = class->therm_bead[ip].x_nhc_tot;
     therm_v_nhc = class->therm_bead[ip].v_nhc;
#ifdef FIX_OR_NUKE_ME
     vpotnhc += (gkt*x_nhc_tot);   
#endif 
     for(ichain=1;ichain<=len_nhc;ichain++){
       for(inhc=mytherm_start;inhc<=mytherm_end;inhc++){
         kinet_nhc += (therm_mass_nhc[ichain][inhc]
                *therm_v_nhc[ichain][inhc]*therm_v_nhc[ichain][inhc]);
         vpotnhc += gkt*class->therm_bead[ip].x_nhc[ichain][inhc];
       /*endfor*/}
     /*endfor*/}
    /*endfor*/}
    kinet_nhc /= 2.0;
    general_data->stat_avg.kinet_nhc_bead = kinet_nhc;
    general_data->stat_avg.vpotnhc = vpotnhc;

/*==========================================================================*/
/* II) vol stuff */


   if(pi_beads_proc_st == 1){
      if(iflag==1){
        general_data->stat_avg.kinet_v = general_data->baro.mass_lnv
                                        *general_data->baro.v_lnv*
                                   general_data->baro.v_lnv/2.0;
        general_data->stat_avg.vpot_v  = general_data->statepoint.pext
                                        *general_data->baro.vol
                                       - general_data->statepoint.stens_ext
                                        *general_data->baro.area;
      /*endif*/}
      if(iflag==2){
        general_data->stat_avg.kinet_v = 0.0;
        for(i=1;i<=9;i++){
          general_data->stat_avg.kinet_v += general_data->par_rahman.mass_hm
                                     *general_data->par_rahman.vgmat[i]
                                     *general_data->par_rahman.vgmat[i];
        /*endfor*/}
        general_data->stat_avg.kinet_v /= 2.0;
        general_data->stat_avg.vpot_v   = general_data->statepoint.pext
                                         *general_data->par_rahman.vol
                                        - general_data->statepoint.stens_ext
                                         *general_data->par_rahman.area;
      /*endif*/}
   }/* endif rank */

/*------------------------------------------------------------------------*/
/*end routine*/}
/*========================================================================*/






