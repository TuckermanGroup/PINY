/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: transform                                    */
/*                                                                          */
/* This subprogram transforms between cartesian and normal mode coords      */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_pimd_local.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_friend_lib_entry.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_force_stag(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran )
/*==========================================================================*/
{/*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iii;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *rat1 = clatoms_tran->rat1_stag;
  double *fxu,*fyu,*fzu;
  double *fxu_p,*fyu_p,*fzu_p;

/*==========================================================================*/
/* II) Get mode forces */



  if(pi_beads!=1){
    fxu_p  = clatoms_pos[1].fx;
    fyu_p  = clatoms_pos[1].fy;
    fzu_p  = clatoms_pos[1].fz;
    for(ip=2;ip<=pi_beads;ip++){
      fxu  = clatoms_pos[ip].fx;
      fyu  = clatoms_pos[ip].fy;
      fzu  = clatoms_pos[ip].fz;
      for(ipart=1;ipart<=natm_tot;ipart++){
        fxu_p[ipart] += fxu[ipart];
        fyu_p[ipart] += fyu[ipart];
        fzu_p[ipart] += fzu[ipart];
      }/*endfor*/
    }/*endfor*/
    for(ip=2;ip<=pi_beads-1;ip++){
      ip1 = ip+1;
      fxu    = clatoms_pos[ip].fx;
      fyu    = clatoms_pos[ip].fy;
      fzu    = clatoms_pos[ip].fz;
      fxu_p  = clatoms_pos[ip1].fx;
      fyu_p  = clatoms_pos[ip1].fy;
      fzu_p  = clatoms_pos[ip1].fz;
      for(ipart=1;ipart<=natm_tot;ipart++){
        fxu_p[ipart] += (rat1[ip]*fxu[ipart]);
        fyu_p[ipart] += (rat1[ip]*fyu[ipart]);
        fzu_p[ipart] += (rat1[ip]*fzu[ipart]);
      }/*endfor*/
    }/*endfor:ipart*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_mode_stag(
              CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              CLATOMS_TRAN *clatoms_tran )

/*==========================================================================*/
{ /*begin routine*/
/*======================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iii;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *rat1 = clatoms_tran->rat1_stag;
  double *rat2 = clatoms_tran->rat2_stag;
  double *xu,*yu,*zu;
  double *xu_p,*yu_p,*zu_p;
  double *xu_1,*yu_1,*zu_1;
  double *x,*y,*z;
  double xstar,ystar,zstar;

/*==========================================================================*/
/* I) Get mode stage */

  if(pi_beads!=1){
    xu_1   = clatoms_pos[1].x;
    yu_1   = clatoms_pos[1].y;
    zu_1   = clatoms_pos[1].z;
    for(ip=2;ip<=pi_beads-1;ip++){
      ip1 = ip+1;
      xu   = clatoms_pos[ip].x;
      yu   = clatoms_pos[ip].y;
      zu   = clatoms_pos[ip].z;
      xu_p = clatoms_pos[ip1].x;
      yu_p = clatoms_pos[ip1].y;
      zu_p = clatoms_pos[ip1].z;
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        xstar = rat1[ip]*xu_p[ipart] + rat2[ip]*xu_1[ipart];
        ystar = rat1[ip]*yu_p[ipart] + rat2[ip]*yu_1[ipart];
        zstar = rat1[ip]*zu_p[ipart] + rat2[ip]*zu_1[ipart];
        xu[ipart] -=  xstar;
        yu[ipart] -=  ystar;
        zu[ipart] -=  zstar;
      }/*endfor:ipart*/
    }/*endfor:ip*/
      ip = pi_beads;
      xu   = clatoms_pos[ip].x;
      yu   = clatoms_pos[ip].y;
      zu   = clatoms_pos[ip].z;
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        xstar = rat1[ip]*xu_1[ipart] + rat2[ip]*xu_1[ipart];
        ystar = rat1[ip]*yu_1[ipart] + rat2[ip]*yu_1[ipart];
        zstar = rat1[ip]*zu_1[ipart] + rat2[ip]*zu_1[ipart];
        xu[ipart] -=  xstar;
        yu[ipart] -=  ystar;
        zu[ipart] -=  zstar;
      }/*endfor:ipart*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void convert_pimd_pos_stag(
                CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
                CLATOMS_TRAN *clatoms_tran )

/*==========================================================================*/
{ /*begin routine*/
/*==========================================================================*/
/*           Local variable declarations                                */

  int ip,ipart,ip1,iii;
  int pi_beads = clatoms_info->pi_beads;
  int natm_tot = clatoms_info->natm_tot;
  int myatm_start = clatoms_info->myatm_start;
  int myatm_end = clatoms_info->myatm_end;
  double *rat1 = clatoms_tran->rat1_stag;
  double *rat2 = clatoms_tran->rat2_stag;
  double *xu,*yu,*zu;
  double *x,*y,*z;
  double *x_p,*y_p,*z_p;
  double *x_1,*y_1,*z_1;
  double xadd,yadd,zadd;

/*==========================================================================*/
/* II) Get cartesian positions */

  if(pi_beads!=1){
    x_1 = clatoms_pos[1].x;
    y_1 = clatoms_pos[1].y;
    z_1 = clatoms_pos[1].z;
    ip = pi_beads;
      x = clatoms_pos[ip].x;
      y = clatoms_pos[ip].y;
      z = clatoms_pos[ip].z;
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        x[ipart] += x_1[ipart];
        y[ipart] += y_1[ipart];
        z[ipart] += z_1[ipart];
      }/*endfor*/
    for(ip=pi_beads-1;ip>=2;ip--){
      ip1 = ip+1;
      x   = clatoms_pos[ip].x;
      y   = clatoms_pos[ip].y;
      z   = clatoms_pos[ip].z;
      x_p = clatoms_pos[ip1].x;
      y_p = clatoms_pos[ip1].y;
      z_p = clatoms_pos[ip1].z;
      for(ipart=myatm_start;ipart<=myatm_end;ipart++){
        xadd = rat1[ip]*x_p[ipart]+x_1[ipart]*rat2[ip];
        yadd = rat1[ip]*y_p[ipart]+y_1[ipart]*rat2[ip];
        zadd = rat1[ip]*z_p[ipart]+z_1[ipart]*rat2[ip];
        x[ipart] += xadd;
        y[ipart] += yadd;
        z[ipart] += zadd;
      }/*endfor*/
    }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void stage_1_part(double *x,double *gaussx,double amass,double text,
            double rcut,int pi_beads,int *iseed,int *iseed2,double *qseed,
            int *ierr_ret)
/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

      double tau,dsum;
      int nmove,iii,i,not_done,ierr,done;
      double endlx,endhx,prer,prerp,sig2,sig,ameanx,dist,xcm;
      double dpi_beads = (double)pi_beads;
/*======================================================================*/
/*               Spread out the particle                            */

   ierr = 0;
   x[1] = 0.0;
   if(pi_beads>1){
     nmove = pi_beads-1;
     endlx = x[1];
     endhx = x[1];
     prer  = (double)(nmove)/(double)(nmove+1);
     prerp = 1./(double)(nmove+1);
     tau   = BOLTZ/(text*dpi_beads);
     sig2  = prer*tau/amass;
     sig   = sqrt(sig2);
     ameanx = prer*endlx + prerp*endhx;
     gaussran(nmove,iseed,iseed2,qseed,gaussx);
     iii = 1;
     for(i=1;i<=nmove;i++){
       not_done=0;done=0;
       while(not_done<=5&&done==0){
         x[(i+1)] =  gaussx[iii]*sig + ameanx;
         iii++;
         dist   = fabs(x[i+1]-x[1]);
         if(iii>nmove){gaussran(nmove,iseed,iseed2,qseed,gaussx);iii=1;}
         if(dist<=rcut){done++;}else{not_done++;}
       }/*endwhile*/
       if(not_done>5){ierr++;}
       prer  = (double)(nmove-i)/(double)(nmove+1-i);
       prerp = 1./(double)(nmove+1-i);
       sig2  = prer*tau/amass;
       sig   = sqrt(sig2);
       ameanx = prer*x[i+1] + prerp*endhx;
     }/*endfor*/
     xcm = dsum1(pi_beads,x,1)/(double)(pi_beads);
     for(i=1;i<=pi_beads;i++){
       x[i] = x[i] - xcm;
     }/*endfor*/
   }/*endif*/
   *ierr_ret = ierr;
/*========================================================================*/
} /* end routine */
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void spread_coord(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  double *x,double *y,double *z,
                  int *iseed,int *iseed2,double *qseed,
                  ATOMMAPS *atommaps,int myid)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

  double *pos,*gauss;
  int ipart,ip,ierr,ierr_now,iii;
  double rcut;  

  int pi_beads       = clatoms_info->pi_beads;
  int pi_beads_proc  = clatoms_info->pi_beads_proc;
  int ip_start       = clatoms_info->pi_beads_proc_st;
  int natm_tot       = clatoms_info->natm_tot;
  int myatm_start    = clatoms_info->myatm_start;
  int myatm_end      = clatoms_info->myatm_end;
  double pi_temperature;
  double *mass       = clatoms_info->mass;

  int pimd_freez_typ = atommaps->pimd_freez_typ;
  int *freeze_flag   = atommaps->freeze_flag;

  int ip_off;
  ip_off = ip_start-1;

/*========================================================================*/
/*  I)Malloc local variables:                                             */

  pi_temperature = clatoms_info->pi_temperature*5.0;
  pos   = (double *) cmalloc(pi_beads*sizeof(double))-1;
  gauss = (double *) cmalloc(pi_beads*sizeof(double))-1;

/*========================================================================*/
/*  II)Spread out classical coords to path integral coords:               */

  ierr = 0;
  rcut = clatoms_info->rcut_spread;
  for(ipart=1;ipart<=natm_tot;ipart++){

    if( (pimd_freez_typ==1) || (freeze_flag[ipart]==0)){
      stage_1_part(pos,gauss,mass[ipart],pi_temperature,rcut, 
                   pi_beads,iseed,iseed2,qseed,&ierr_now);
      ierr += ierr_now;
      for(ip=1;ip<=pi_beads_proc;ip++){
        clatoms_pos[ip].x[ipart] = pos[(ip+ip_off)] + x[ipart];
      }/*endfor*/

      stage_1_part(pos ,gauss,mass[ipart],pi_temperature,rcut, 
                   pi_beads,iseed,iseed2,qseed,&ierr_now);
      ierr += ierr_now;
      for(ip=1;ip<=pi_beads_proc;ip++){
        clatoms_pos[ip].y[ipart] = pos[(ip+ip_off)] + y[ipart];
      }/*endfor*/

      stage_1_part(pos ,gauss,mass[ipart],pi_temperature,rcut, 
                   pi_beads,iseed,iseed2,qseed,&ierr_now);
      ierr += ierr_now;
      for(ip=1;ip<=pi_beads_proc;ip++){
        clatoms_pos[ip].z[ipart] = pos[(ip+ip_off)] + z[ipart];
      }/*endfor*/

    }else{

      for(ip=1;ip<=pi_beads_proc;ip++){
        clatoms_pos[ip].x[ipart] = x[ipart];
        clatoms_pos[ip].y[ipart] = y[ipart];
        clatoms_pos[ip].z[ipart] = z[ipart];
      }/*endfor*/

    }/*endif*/

  }/*endfor*/

/*========================================================================*/
/* III) Tuck variables back and free memory                               */

  cfree(&pos[1]);
  cfree(&gauss[1]);

/*========================================================================*/
/* IV) Error mesage */

  if((ierr>0) && (myid==0)){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("%d coordinates outside spread cutoff\n",ierr);
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

/*========================================================================*/
    }/* end routine */
/*==========================================================================*/






