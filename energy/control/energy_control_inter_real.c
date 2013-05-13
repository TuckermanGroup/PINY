/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                   Module: energy_control.c                               */
/*                                                                          */
/* This routine calls the required force and PE routines                    */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_energy_ctrl_entry.h"
#include "../proto_defs/proto_intra_entry.h"
#include "../proto_defs/proto_real_space_entry.h"
#include "../proto_defs/proto_recip3d_entry.h"
#include "../proto_defs/proto_energy_ctrl_local.h"
#include "../proto_defs/proto_intra_con_entry.h"
#include "../proto_defs/proto_pimd_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_output_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void energy_control_inter_real(CLASS *class, BONDED *bonded, 
                               GENERAL_DATA *general_data)

/*==========================================================================*/
   {/*Begin Routine*/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"
  
  double vreal,vlong,vvdw,vcoul;
  double vcoul_temp,vreal_temp;
  int iii,ip;

  double pext            = general_data->statepoint.pext;
  double text            = general_data->statepoint.t_ext;

  double vol             = general_data->cell.vol;
  double vol0            = general_data->cell.vol0;
  double *hmat           = general_data->cell.hmat;
  int iperd              = general_data->cell.iperd;

  int natm_tot           = class->clatoms_info.natm_tot;
  int nfree              = class->clatoms_info.nfree;
  int pi_beads           = class->clatoms_info.pi_beads;  
  int pi_beads_proc      = class->clatoms_info.pi_beads_proc;  

  double *x              = class->clatoms_pos[1].x;  
  double *y              = class->clatoms_pos[1].y;  
  double *z              = class->clatoms_pos[1].z;  

  int error_check_on     = general_data->error_check_on;
  MPI_Comm comm_beads    = class->communicate.comm_beads;
  int myid               = class->communicate.myid;
  int np_beads           = class->communicate.np_beads;

  double cutoff_max      = class->interact.cutoff_max;

  int iget_full_inter    = class->energy_ctrl.iget_full_inter;
  int iget_res_inter     = class->energy_ctrl.iget_res_inter;
  int iget_full_intra    = class->energy_ctrl.iget_full_intra;
  int iget_pe_real_inter = class->energy_ctrl.iget_pe_real_inter;


  int pimd_on,cp_on,nfree3;

/*======================================================================*/
/* 0) Get some constants */

  pimd_on = general_data->simopts.pimd
           +general_data->simopts.debug_pimd
           +general_data->simopts.cp_pimd
           +general_data->simopts.cp_wave_pimd
           +general_data->simopts.cp_wave_min_pimd
           +general_data->simopts.debug_cp_pimd;
  cp_on   = general_data->simopts.cp_min  
           +general_data->simopts.cp_wave_min
           +general_data->simopts.cp      
           +general_data->simopts.cp_wave
           +general_data->simopts.cp_pimd
           +general_data->simopts.cp_wave_pimd
           +general_data->simopts.debug_cp
           +general_data->simopts.debug_cp_pimd
           +general_data->simopts.cp_wave_min_pimd;

  nfree3  = nfree/3;
  if(pimd_on==1){nfree3=natm_tot;}

/*======================================================================*/
/* I) Initialize */

  vlong            = 0.0;
  vvdw             = 0.0;
  vcoul            = 0.0;
  vreal            = 0.0;

/*======================================================================*/
/* III) Get intermolecular real space force and  PE   */

  if( (iget_full_inter==1) || (iget_res_inter ==1) ){

    if(iperd>0 ){check_cutoff(myid,iperd,hmat,cutoff_max);}
    if(iperd==0){check_cutoff_clus(myid,x,y,z,natm_tot,cutoff_max);}

    nbr_list_control(&(class->clatoms_info),(class->clatoms_pos), 
                     &(class->for_scr),
                     &(class->atommaps),&(general_data->cell),
                     &(class->interact),&(general_data->timeinfo),
                     &(class->nbr_list),&(bonded->excl),    
                     &(bonded->intra_scr),&(general_data->stat_avg),
                     &(class->communicate),error_check_on,
                     &(class->class_comm_forc_pkg));

    for(ip=1;ip<=pi_beads_proc;ip++){
       force_control(&(class->clatoms_info), &(class->clatoms_pos[ip]),
                     &(class->for_scr),
                     &(class->atommaps),&(general_data->cell),
                     &(general_data->ptens),   &(class->interact),
                     &(class->energy_ctrl),&(class->nbr_list),
                     &(bonded->excl),    &(bonded->intra_scr),&vreal,
                     &vvdw,&vcoul,error_check_on,
                     &(class->class_comm_forc_pkg));
    }/*endfor*/

    vreal /=pi_beads;
    vvdw  /=pi_beads;
    vcoul /=pi_beads;

  }/*endif*/

/*======================================================================*/
/* II) Get the long range correction   */

  if( (iperd==3) && (myid==0) ){

     (class->interact.pten_kin_guess) = (((double)nfree3)*text)/(BOLTZ*vol0);
     long_range_corr(pi_beads,&vlong,vol,&(general_data->ptens),
                     &(class->for_scr),
                     &(class->interact),&(class->energy_ctrl),pext);

  }/*endif*/

/*======================================================================*/
/* II.V) Get the coloumb correction due to the CP-CP atoms              */
/*       which were included in exclusions list whose coloumb energy    */
/*       should be included                                             */
/*  NOTE: THIS HAS NOT BEEN IMPLEMENTED FOR PATH INTEGRALS              */

  if( bonded->excl.num_cp > 0 && pi_beads == 1){
     mix_coul_corr(&(class->clatoms_info),&(class->clatoms_pos[1]),
                   &(bonded->intra_scr),&(bonded->excl),
                   &(class->interact),
                   &(general_data->cell),&(general_data->ptens),
                   &vreal,general_data->ewald.alp_ewd); 
  }/*endif*/

/*======================================================================*/
/* III) Communicate and store  */

  if((cp_on==1)&&(np_beads>1)){
     vcoul_temp = 0.0;   vreal_temp       = 0.0;
     Allreduce(&(vreal), &(vreal_temp),1,MPI_DOUBLE,MPI_SUM,0,comm_beads);
     Allreduce(&(vcoul), &(vcoul_temp),1,MPI_DOUBLE,MPI_SUM,0,comm_beads);
     vcoul = vcoul_temp; vreal = vreal_temp;    
  }/*endif*/

  if( (iget_full_intra==1)|| (iget_full_inter==1) ){

    general_data->stat_avg.vlong          = vlong;
    if(iget_pe_real_inter==1){
      (general_data->stat_avg.vintert)    += (vreal+vlong);
      (general_data->stat_avg.vvdw)       += (vvdw +vlong);
      (general_data->stat_avg.vcoul)      += (vcoul);
    }/*endif*/

  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void long_range_corr(int pi_beads, double *vlong_ret, double vol,PTENS *ptens, 
                     FOR_SCR *for_scr,
                     INTERACT *interact, ENERGY_CTRL *energy_ctrl, double 
                     pext)

/*======================================================================*/
   {/*Begin Routine*/
/*======================================================================*/
/*    Local Variables */

  double vlong_now;
  double pdiag_ter_res,pdiag_tra_res,pdiag,pdiag_full;
  double pdiag_bead_res;
  double amult,amult_bead,dpi_beads;
  int iii;
  
  double *pvten          = ptens->pvten;
  double *pvten_tot      = ptens->pvten_tot;

  double clong           = interact->clong;
  double pten_inter_guess= interact->pten_inter_guess;
  double pkin_guess      = interact->pten_kin_guess;
  double clong_res       = interact->clong_res;

  double wght_ter        = for_scr->wght_ter;
  double wght_ter_res    = for_scr->wght_ter_res;
  double wght_tra        = for_scr->wght_tra;
  double wght_tra_res    = for_scr->wght_tra_res;
  
  int iswit_vdw          = energy_ctrl->iswit_vdw;
  int iget_res_intra     = energy_ctrl->iget_res_intra;
  int iget_res_inter     = energy_ctrl->iget_res_inter;
  int iget_full_inter    = energy_ctrl->iget_full_inter;
  int iget_pv_real_inter = energy_ctrl->iget_pv_real_inter;

/*======================================================================*/
/* 0) Set the multiplicative factor       */

  amult      = 2.0;  if(iswit_vdw==1){amult=1.0;}
  dpi_beads  = (double) pi_beads;
  amult_bead =  amult*dpi_beads;

/*======================================================================*/
/* I) Respa Bead step */

  if((pi_beads>1)&&(wght_tra_res!=1.0)){

     pdiag_bead_res = (pext-pkin_guess)*vol*dpi_beads;
     pvten[1]     += pdiag_bead_res;
     pvten[5]     += pdiag_bead_res;
     pvten[9]     += pdiag_bead_res;

  }else{ 

    pdiag_bead_res = 0.0;

  }/*endif*/

/*======================================================================*/
/* I) Respa Intra step */

  if(((iget_res_intra)==1)&&(wght_ter!=1.0)){

     pdiag_tra_res = (pten_inter_guess*vol)*dpi_beads;
     pvten[1]     += (pdiag_tra_res-pdiag_bead_res)*(wght_tra_res);
     pvten[5]     += (pdiag_tra_res-pdiag_bead_res)*(wght_tra_res);
     pvten[9]     += (pdiag_tra_res-pdiag_bead_res)*(wght_tra_res);

  }else{ 

    pdiag_tra_res = pdiag_bead_res;

  }/*endif*/

/*======================================================================*/
/* II) Respa Inter step */

  if( (iget_res_inter)==1){

     pdiag_ter_res = (clong_res/vol)*amult_bead;
     pvten[1]     += (pdiag_ter_res-pdiag_tra_res)*(wght_ter_res);
     pvten[5]     += (pdiag_ter_res-pdiag_tra_res)*(wght_ter_res);
     pvten[9]     += (pdiag_ter_res-pdiag_tra_res)*(wght_ter_res);

  }else{   

     pdiag_ter_res = pdiag_tra_res;

  }/*endif*/

/*======================================================================*/
/* III) Full inter respa */

  if( (iget_full_inter)==1){

     pdiag_full = (clong/vol)*amult_bead;
     pvten[1]     += (pdiag_full-pdiag_ter_res)*wght_ter;
     pvten[5]     += (pdiag_full-pdiag_ter_res)*wght_ter;
     pvten[9]     += (pdiag_full-pdiag_ter_res)*wght_ter;
     
     (*vlong_ret)  = clong/vol;
     if(iget_pv_real_inter==1){
      pvten_tot[1] += pdiag_full;
      pvten_tot[5] += pdiag_full;
      pvten_tot[9] += pdiag_full;
     }

   }/*endif*/

/*-------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_cutoff(int myid,int iperd,double *hmat,double cutoff_max)

/*=======================================================================*/
/*            Begin subprogram:                                          */
    {/*begin routine*/
/*==========================================================================*/
/* Local variables */

  int ierr;
  double a,b,c,tab,tbc,tac;

/*==========================================================================*/
/* Check the cutoff */

  get_cell(hmat,&a,&b,&c,&tab,&tbc,&tac);

  ierr = 0;
  if( (iperd>=1) && (cutoff_max>0.5*a) ){ierr++;}
  if( (iperd>=2) && (cutoff_max>0.5*b) ){ierr++;}
  if( (iperd>=3) && (cutoff_max>0.5*c) ){ierr++;}

  if( (ierr>0) && (myid==0) ){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("The max cutoff, %g, exceeds half the box edge, %g %g %g\n",
            cutoff_max,a,b,c);
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

/*------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_cutoff_clus(int myid,double *x,double *y,double *z, int natm_tot,
                       double cutoff_max)

/*==========================================================================*/
 { /*begin routine */
/*==========================================================================*/

  double dx,dy,dz,r2;
  double xcm,ycm,zcm;
  double rmax2,rmax;
  int i;

/*==========================================================================*/
/* I) Find the com */

  xcm = 0.0;
  ycm = 0.0;
  zcm = 0.0;
  for(i=1;i<=natm_tot;i++){
    xcm += x[i];
    ycm += y[i];
    zcm += z[i];
  }/*endfor*/
  xcm /= ((double)natm_tot);
  ycm /= ((double)natm_tot);
  zcm /= ((double)natm_tot);

/*==========================================================================*/
/* II) Find the maximum radial distance  */

  rmax2 = 0.0;
  for(i=1;i<=natm_tot;i++){
    dx = x[i]-xcm;
    dy = y[i]-ycm;
    dz = z[i]-zcm;
    r2 = dx*dx+dy*dy+dz*dz;
    rmax2 = MAX(rmax2,r2);
  }/*endfor*/
  rmax = sqrt(rmax2);

/*==========================================================================*/
/* III) Print the error */

  if( (cutoff_max<2.0*rmax) && (myid==0) ){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("The max cutoff, %g, is less twice the max radial distance, %g\n",
            cutoff_max,rmax);
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
  }/*endif*/

/*==========================================================================*/
  }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mix_coul_corr(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   INTRA_SCR *intra_scr,EXCL *excl,INTERACT *interact,
                   CELL *cell,PTENS *ptens,double *vcoul,double alp)

/*==========================================================================*/
 { /*begin routine */
/*==========================================================================*/

  int i,iii;
  int j1_cp_now,j2_cp_now;
  double v,dv,vnow,dvnow;
  double palp, ralp, talp2, gerfc, r2, tt, eee,dgerfc;
  double p=0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;
  double de1,de2,de3;
  double de4,de5,de6;
  double de7,de8,de9;
  double rt,s_eps,eps_r,deps_r;
  double fx,fy,fz;
  double dx12,dy12,dz12,r;
  double qij;
  double diele_cut    = interact->dielectric_cut;
  double diele_opt    = interact->dielectric_opt;
  double diele_rheal  = interact->dielectric_rheal;
  double diele_eps    = interact->dielectric_eps;

  double *pvten           = ptens->pvten;
  double *ptens_pvten_tmp = ptens->pvten_tmp;

  double *clatoms_q    = clatoms_info->q;

  double *clatoms_x    = clatoms_pos->x;
  double *clatoms_y    = clatoms_pos->y;
  double *clatoms_z    = clatoms_pos->z;

  double *clatoms_fx   = clatoms_pos->fx;
  double *clatoms_fy   = clatoms_pos->fy;
  double *clatoms_fz   = clatoms_pos->fz;

  int num_cp    = excl->num_cp;
  int *j1_cp    = excl->j1_cp;
  int *j2_cp    = excl->j2_cp; 

  int iperd     = cell->iperd;

  de1 = 1.0*e1; de2 = 2.0*e2; de3 = 3.0*e3;
  de4 = 4.0*e4; de5 = 5.0*e5; de6 = 6.0*e6;
  de7 = 7.0*e7; de8 = 8.0*e8; de9 = 9.0*e9;

/*==========================================================================*/

  for(i=1;i<=9;i++){(ptens_pvten_tmp)[i]=0;}

/*==========================================================================*/
/* Regular old Mr. c */

  if (iperd == 0) {

   if(diele_opt==0){
    for(i=1;i<= num_cp;i++){
      /* Make distance rij */
      j1_cp_now = j1_cp[i];   j2_cp_now = j2_cp[i];
     
      dx12  = clatoms_x[(j1_cp_now)] - clatoms_x[(j2_cp_now)];
      dy12  = clatoms_y[(j1_cp_now)] - clatoms_y[(j2_cp_now)];
      dz12  = clatoms_z[(j1_cp_now)] - clatoms_z[(j2_cp_now)];
       
      r2    = dx12*dx12 + dy12*dy12 + dz12*dz12;
      r     = sqrt(r2);
      qij   = clatoms_q[(j1_cp_now)]*clatoms_q[(j2_cp_now)];

      vnow  = qij/r;      
      dvnow = qij/(r*r*r);
      v  = vnow;
      dv = dvnow;

      *vcoul += v;  

      fx = dx12*dv;
      fy = dy12*dv;
      fz = dz12*dv;

      /* Particle i */
      clatoms_fx[(j1_cp_now)] += fx;
      clatoms_fy[(j1_cp_now)] += fy;
      clatoms_fz[(j1_cp_now)] += fz;


      /* Particle j */
      clatoms_fx[(j2_cp_now)] -= fx;
      clatoms_fy[(j2_cp_now)] -= fy;
      clatoms_fz[(j2_cp_now)] -= fz;

    }/*endfor*/
  }/*endif*/
 }/*endif*/


/*==========================================================================*/
/* Real space ewald */

   if (iperd > 0) { 
    talp2 = 2.0*alp*alp;
    palp  = p*alp;

    for(i=1; i <= num_cp ;i++){
      /* Make distance rij */

      j1_cp_now = j1_cp[i];   j2_cp_now = j2_cp[i];
      dx12  = clatoms_x[(j1_cp_now)] - clatoms_x[(j2_cp_now)];
      dy12  = clatoms_y[(j1_cp_now)] - clatoms_y[(j2_cp_now)];
      dz12  = clatoms_z[(j1_cp_now)] - clatoms_z[(j2_cp_now)];

      period_one(1,&dx12,&dy12,&dz12,cell);       

      r2    = dx12*dx12 + dy12*dy12 + dz12*dz12;
      r     = sqrt(r2);
      qij   = clatoms_q[(j1_cp_now)]*clatoms_q[(j2_cp_now)];

      r2     = r*r;
      ralp   = r * alp;
      eee    = exp(-ralp*ralp);
      tt     = 1.0/(1.0+p*ralp);
      gerfc  = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                      +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
      dgerfc = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt
                           +de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*palp
                           +talp2*gerfc*r;
      vnow   = qij * gerfc/r;               
      dvnow  = (gerfc/r2 + dgerfc/r)*qij/r; 
      v      = vnow;
      dv     = dvnow;

      *vcoul += v;

      fx = dx12*dv;
      fy = dy12*dv;
      fz = dz12*dv;

      /* Particle i */
      clatoms_fx[(j1_cp_now)] += fx;
      clatoms_fy[(j1_cp_now)] += fy;
      clatoms_fz[(j1_cp_now)] += fz;

      /* Particle j */
      clatoms_fx[(j2_cp_now)] -= fx;
      clatoms_fy[(j2_cp_now)] -= fy;
      clatoms_fz[(j2_cp_now)] -= fz; 

     if(iperd == 2 || iperd == 3) {
      ptens_pvten_tmp[1] += dx12*fx;   /*p11*/
      ptens_pvten_tmp[5] += dy12*fy;   /*p22*/
      ptens_pvten_tmp[9] += dz12*fz;   /*p33*/
      ptens_pvten_tmp[2] += dx12*fy;   /*p12*/
     }/*endif*/

     if(iperd == 3) {
      ptens_pvten_tmp[3] += dx12*fz;    /*p13*/
      ptens_pvten_tmp[6] += dy12*fz;    /*p23*/
     }/*endif*/


    }/*endfor*/
   }/*endif*/

/*==========================================================================*/
  }/*end routine */
/*==========================================================================*/








