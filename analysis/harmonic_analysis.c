/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: harmonic_analysis.c                          */
/*                                                                          */
/* Performs harmonic analysis of one or more configurations                 */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include <stddef.h>
#include "standard_include.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_analysis_local_entry.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_math.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void calcul_freqs(CLASS *class,GENERAL_DATA *general_data,ANALYSIS *analysis)

/*=======================================================================*/
{/*begin function*/

/*=======================================================================*/
/*            Local variable declarations                                */
/*=======================================================================*/

  double *scr_mat;
  double *scr_mat_mass;
  double *hess_vecs;
  double *rs1,*rs2;
  double *harm_freqs;
  double junk;
  double mass_factor;
  int ipart,jpart,indh,indm;
  int count;
  int natm_tot = class->clatoms_info.natm_tot;
  int itime = general_data->timeinfo.itime;
  int nhess,job,ierr;

  double *hess_xx = class->clatoms_pos[1].hess_xx;
  double *hess_xy = class->clatoms_pos[1].hess_xy;
  double *hess_xz = class->clatoms_pos[1].hess_xz;
  double *hess_yy = class->clatoms_pos[1].hess_yy;
  double *hess_yz = class->clatoms_pos[1].hess_yz;
  double *hess_zz = class->clatoms_pos[1].hess_zz;
  double *mass    = class->clatoms_info.mass;

/*========================================================================*/
/* I) Malloc up some local scratch                                        */


 scr_mat      = (double *) cmalloc(9*natm_tot*natm_tot*sizeof(double))-1;
 scr_mat_mass = (double *) cmalloc(9*natm_tot*natm_tot*sizeof(double))-1;
 hess_vecs    = (double *) cmalloc(9*natm_tot*natm_tot*sizeof(double))-1;
 rs1     = (double *) cmalloc(3*natm_tot*sizeof(double))-1;
 rs2     = (double *) cmalloc(3*natm_tot*sizeof(double))-1;

/*========================================================================*/
/* II) If this is the first step, malloc space for frequencies             */


  analysis->harmonic_analysis.harm_freqs = (double *) cmalloc(3*natm_tot*sizeof(double))-1;

/*========================================================================*/
/* III) Assign the local pointer                                          */


  harm_freqs = analysis->harmonic_analysis.harm_freqs;

/*========================================================================*/
/* IV) Pack the scratch matrices                                          */

/*-------------------------------------------------------------------------*/
/* i) First the Hessian */

 count = 1;
 for(ipart=1;ipart<=natm_tot;ipart++){
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] = hess_xx[indh]; count++;
      scr_mat[count] = hess_xy[indh]; count++;
      scr_mat[count] = hess_xz[indh]; count++;
   }
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] = hess_xy[indh]; count++;
      scr_mat[count] = hess_yy[indh]; count++;
      scr_mat[count] = hess_yz[indh]; count++;
   }         
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] = hess_xz[indh]; count++;
      scr_mat[count] = hess_yz[indh]; count++;
      scr_mat[count] = hess_zz[indh]; count++;
   }
 }

#ifdef JUNK
 count = 1;
 for(ipart=1;ipart<=natm_tot;ipart++){
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
   }
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
   }         
   for(jpart=1;jpart<=natm_tot;jpart++){
      indh = (ipart-1)*natm_tot + jpart;
      mass_factor = sqrt(mass[ipart]*mass[jpart]);
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
      scr_mat[count] /= mass_factor; count++;
   }
 }
#endif

/*-------------------------------------------------------------------------*/
/* ii) Next the diagonal matrix of masses */


  nhess = 3*natm_tot;
  for(ipart=1;ipart <= nhess*nhess; ipart++){
    scr_mat_mass[ipart] = 0.0;
  }


  for(ipart=1;ipart <= nhess; ipart++){
    indh = (ipart-1)*nhess + ipart;
    indm = (int) (((double) (ipart-1))/3.0) + 1;
    scr_mat_mass[indh] = mass[indm];
  }  


/*========================================================================*/
/* V) Diagonalize the Hessian                                             */

  
  job = 0;
  RSG(&nhess,&nhess,&(scr_mat[1]),&(scr_mat_mass[1]),&(harm_freqs[1]),
      &job,&(hess_vecs[1]),&(rs1[1]),&(rs2[1]),&ierr);

#ifdef JUNK
  RS(&nhess,&nhess,&(scr_mat[1]),&(harm_freqs[1]),
      &job,&(hess_vecs[1]),&(rs1[1]),&(rs2[1]),&ierr);
#endif

  for(ipart=1;ipart <= nhess; ipart++){
    if(harm_freqs[ipart] >= 0) 
         printf("Harmonic frequencies  %.12g\n",sqrt(harm_freqs[ipart])*219476.76);
    if(harm_freqs[ipart] < 0) 
         printf("Harmonic frequencies  %.12g i\n",sqrt(-harm_freqs[ipart])*219476.76);
  }


/*========================================================================*/
/*  Free scratch memory                                                   */

  cfree(&(scr_mat[1]));
  cfree(&(scr_mat_mass[1]));
  cfree(&(hess_vecs[1]));
  cfree(&(rs1[1]));
  cfree(&(rs2[1]));

/*=======================================================================*/
}/*end function*/
/*==========================================================================*/


