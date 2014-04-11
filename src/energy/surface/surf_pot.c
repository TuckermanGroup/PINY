/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: surf_pot.c                                   */
/*                                                                          */
/* This routine computes the energy and forces from the surface.            */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_surf_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void surf_pot(CLATOMS_INFO *clatoms_info, CLATOMS_POS *clatoms_pos,
              ATOMMAPS *atommaps,SURFACE *surface, CELL *cell,
              INTRA_SCR *intra_scr, PTENS *ptens, double *vsurf_ret, 
              int iver_get, CLASS_COMM_FORC_PKG *class_comm_forc_pkg,
              int iget_pv_real_inter,int iget_pe_real_inter)

/*==========================================================================*/ 
  { /* begin routine */
/*==========================================================================*/
/* Local variable declarations                                              */

  int isurf,ibig,ibig1,ioff,iend,nnow;            /* Counters,offsets, etc  */
  int isoff;
  int nlen_now;
  int ktemp,kktemp;

  double dz,zsp,pvten_tmp,vsurf;
  double zmin_spl_now;
  double dz_spl_now;
  double dzi_spl_now;

  /* Define local pointers                                                 */
  double *scr_dz1         = intra_scr->dz12;
  double *scr_fz1         = intra_scr->fz1;
  double *scr_vpot        = intra_scr->vpot;
  double *scr_p33         = intra_scr->p33;
  double  wfor            = intra_scr->wght_tra;
  int nlen_use            = intra_scr->nlen;
  double *del_z           = intra_scr->del_r;
  int *index              = intra_scr->iatm_typ;
  double *spl_tmp         = intra_scr->spl_tmp;

  double *pvten           = ptens->pvten;
  double *pvten_tot       = ptens->pvten_tot;

  double *clatoms_z       = clatoms_pos->z;
  double *clatoms_fz      = clatoms_pos->fz;
  double *clatoms_fzt     = clatoms_pos->fzt;

  int *iatm_atm_typ       = atommaps->iatm_atm_typ;

  int iperd               = cell->iperd;

  int np_forc             = class_comm_forc_pkg->num_proc;
  int myid_forc           = class_comm_forc_pkg->myid;
  int myatm_start         = class_comm_forc_pkg->myatm_start;
  int myatm_end           = class_comm_forc_pkg->myatm_end;
  int ntot                = class_comm_forc_pkg->myatm_end;

  int nsplin              = surface->nsplin_surf;
  double z_surf           = surface->surface_height;
  double *surface_pot     = surface->surface_pot;
  double *surface_forc    = surface->surface_forc;
  double *zmin_spl        = surface->zmin_spl;
  double *dz_spl          = surface->dz_spl;
  double *dzi_spl         = surface->dzi_spl;

  int nsplin_m2;

      nsplin_m2           = nsplin-2;

/*=======================================================================*/
/* I) loop over all the atoms in steps of nlen to save memory            */

  pvten_tmp = 0.0;
  vsurf     = 0.0;

  for(ibig=myatm_start;ibig <= myatm_end;ibig += nlen_use) {
    /*---------------------------------------------------------------------*/
    /*  A) Offsets to save some memory by doing only nlen bonds at a time  */
    
    ibig1 = ibig-1;
    ioff  = -ibig1;
    iend  = MIN(ntot,ibig1+nlen_use);
    nnow  = iend-ibig1;

    /*---------------------------------------------------------------------*/
    /*  B) get the spline indicies                                         */

    for(isoff=1,isurf=ibig;isurf <= iend; ++isurf,++isoff) {       

      dz             = clatoms_z[isurf]-z_surf;

      ktemp          = iatm_atm_typ[isurf];
      zmin_spl_now   = zmin_spl[ktemp];
      dz_spl_now     = dz_spl[ktemp];
      dzi_spl_now    = dzi_spl[ktemp];

      kktemp         = (int) (((dz - zmin_spl_now)*dzi_spl_now)+0.5) + 3;
      kktemp         = MIN(kktemp,nsplin_m2);
      kktemp         = MAX(kktemp,3);
      zsp            = ((double)(kktemp-3))*dz_spl_now+zmin_spl_now;

      del_z[isoff]   = (dz - zsp)*dzi_spl_now;
      index[isoff]   = kktemp;
      scr_dz1[isoff] = dz;

    }/*endfor*/

    /*---------------------------------------------------------------------*/
    /*  C) Sum the potential energy                                        */

    if(iget_pe_real_inter==1){
      surf_vspl_fetch(nnow,del_z,spl_tmp,index,surface_pot);
      for(isoff=1;isoff<=nnow;++isoff){ vsurf += spl_tmp[isoff]; }
    }/*endif*/

    /*--------------------------------------------------------------------*/
    /*  D) Get the forces                                                 */

    surf_vspl_fetch(nnow,del_z,scr_fz1,index,surface_forc);

    /*--------------------------------------------------------------------*/
    /*  E) Scatter the forces                                             */

    if(iver_get==1){

      for(isoff=1,isurf=ibig;isurf <= iend; ++isurf,++isoff) {       
        clatoms_fzt[isurf] += scr_fz1[isoff];
      }/*endfor*/

    }/*endif*/

    for(isoff=1;isoff <= nnow; ++isoff) {
      pvten_tmp += scr_fz1[isoff]*scr_dz1[isoff];
    }/*endfor*/

    for(isoff=1;isoff <= nnow; ++isoff) {
      scr_fz1[isoff] *= wfor;
    }/*endfor*/

    for(isoff=1,isurf=ibig;isurf <= iend; ++isurf,++isoff) {       
      clatoms_fz[isurf]  +=  scr_fz1[isoff];
    }/*endfor*/

    /*--------------------------------------------------------------------*/
 }/*endfor ibig*/

/*==========================================================================*/
/* Increment the Pressure tensor */

  (*vsurf_ret) += vsurf;

  pvten[9]     += pvten_tmp*wfor;
  if(iget_pv_real_inter==1){    
    pvten_tot[9] += pvten_tmp;
  }/*endif*/

/*--------------------------------------------------------------------------*/
   }/* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void surf_vspl_fetch(int n,double del[], double spl_out[],
                     int i_index[],double c0_data[])

/*=======================================================================*/
/*              Begin Routine                                            */
    {/*Begin Routine*/
/*=======================================================================*/
/*            Local variable declarations                                */

  int i,ktemp0,ktemp1,ktemp2,ktemp3,iii;
  double del_tmp,c0,c1,c2,c3;
  double f0,fp1,fm1,fm2;
  static double oneth = (1.0/3.0);
  static double onesi = (1.0/6.0);

  /*========================================================================*/
  /* I) Fetch the coefs */

  for(i=1;i<=n;i++){
    ktemp0 = i_index[i];
    ktemp1 = ktemp0 + 1;
    ktemp2 = ktemp0 - 1;
    ktemp3 = ktemp0 - 2;
    c0  = c0_data[ktemp0];
    fp1 = c0_data[ktemp1];
    fm1 = c0_data[ktemp2];
    fm2 = c0_data[ktemp3];
    c1 = oneth*fp1+0.5*c0-fm1+onesi*fm2;
    c2 = 0.5*fp1-c0+0.5*fm1;
    c3 = onesi*fp1-0.5*c0+0.5*fm1-onesi*fm2;
    del_tmp = del[i];
    spl_out[i]  =  (c0 + del_tmp*(c1 + del_tmp* (c2 + del_tmp*c3)));
  }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/* end routine vspl_fetch */
/*==========================================================================*/


