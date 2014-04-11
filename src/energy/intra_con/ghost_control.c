/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: ghost_control                                */
/*                                                                          */
/* This subprogram routines to calculate positions of ghost atoms and       */
/* force contributions due to ghost atoms                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_con_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void get_ghost_pos(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                   GHOST_ATOMS *ghost_atoms)
/*========================================================================*/
{/*begin routine*/

  int igloc,igpart,ighost,icoef,iii;
/* Local pointers */

   double *clatoms_x           = clatoms_pos->x;
   double *clatoms_y           = clatoms_pos->y;
   double *clatoms_z           = clatoms_pos->z;
   int nghost                  = ghost_atoms->nghost_tot;
   double **ghost_atoms_coef   = ghost_atoms->coef;
   int *ghost_atoms_ighost_map = ghost_atoms->ighost_map;
   int *ghost_atoms_natm_comp  = ghost_atoms->natm_comp;
   int **ghost_atoms_iatm_comp = ghost_atoms->iatm_comp;

/*==========================================================================*/
/* I) Recompute positions of ghost atoms                                    */

 if(nghost > 0) {
  for(ighost=1;ighost <= nghost;ighost++){
   igloc = ghost_atoms_ighost_map[ighost];
   clatoms_x[igloc] = 0.0;
   clatoms_y[igloc] = 0.0;
   clatoms_z[igloc] = 0.0;
   for(icoef=1;icoef<=ghost_atoms_natm_comp[ighost];icoef++) {
     igpart = ghost_atoms_iatm_comp[icoef][ighost];
     clatoms_x[igloc] += ghost_atoms_coef[icoef][ighost]*clatoms_x[igpart];
     clatoms_y[igloc] += ghost_atoms_coef[icoef][ighost]*clatoms_y[igpart];
     clatoms_z[igloc] += ghost_atoms_coef[icoef][ighost]*clatoms_z[igpart];
   }/*endfor icoef */
  }/*endfor ighost */
 }/* endif */
/*-------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void distrib_ghost_force(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                         GHOST_ATOMS *ghost_atoms, int iver_get)
/*==========================================================================*/

{/* begin routine */

 int igloc,igpart,ighost,icoef,nghost=ghost_atoms->nghost_tot;

/* Local pointers */

   double *clatoms_fx = clatoms_pos->fx;
   double *clatoms_fy = clatoms_pos->fy;
   double *clatoms_fz = clatoms_pos->fz;
   double *clatoms_fxt = clatoms_pos->fxt;
   double *clatoms_fyt = clatoms_pos->fyt;
   double *clatoms_fzt = clatoms_pos->fzt;
   double **ghost_atoms_coef = ghost_atoms->coef;
   int *ghost_atoms_ighost_map = ghost_atoms->ighost_map;
   int *ghost_atoms_natm_comp = ghost_atoms->natm_comp;
   int **ghost_atoms_iatm_comp = ghost_atoms->iatm_comp;

/*==========================================================================*/
/* I) Update forces due to ghost atoms and zero ghost atom forces         */

 if(nghost > 0) {
  for(ighost=1;ighost <= nghost;ighost++){
    igloc = ghost_atoms_ighost_map[ighost];
    for(icoef=1;icoef<=ghost_atoms_natm_comp[ighost];icoef++) {
      igpart = ghost_atoms_iatm_comp[icoef][ighost];
      clatoms_fx[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fx[igloc];
      clatoms_fy[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fy[igloc];
      clatoms_fz[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fz[igloc];
    } /* endfor icoef */
    clatoms_fx[igloc] = 0.0;
    clatoms_fy[igloc] = 0.0;
    clatoms_fz[igloc] = 0.0;
  } /* endfor ighost */
 }/* endif */
 if((nghost > 0)&&(iver_get==1)) {
  for(ighost=1;ighost <= nghost;ighost++){
    igloc = ghost_atoms_ighost_map[ighost];
    for(icoef=1;icoef<=ghost_atoms_natm_comp[ighost];icoef++) {
      igpart = ghost_atoms_iatm_comp[icoef][ighost];
      clatoms_fxt[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fxt[igloc];
      clatoms_fyt[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fyt[igloc];
      clatoms_fzt[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fzt[igloc];
    } /* endfor icoef */
    clatoms_fxt[igloc] = 0.0;
    clatoms_fyt[igloc] = 0.0;
    clatoms_fzt[igloc] = 0.0;
  } /* endfor ighost */
 }/* endif */
/*-------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void distrib_ghost_force_mode(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                              GHOST_ATOMS *ghost_atoms)
/*==========================================================================*/

{/* begin routine */

 int igloc,igpart,ighost,icoef,nghost=ghost_atoms->nghost_tot;

/* Local pointers */

   double *clatoms_fx = clatoms_pos->fx;
   double *clatoms_fy = clatoms_pos->fy;
   double *clatoms_fz = clatoms_pos->fz;
   double **ghost_atoms_coef = ghost_atoms->coef;
   int *ghost_atoms_ighost_map = ghost_atoms->ighost_map;
   int *ghost_atoms_natm_comp = ghost_atoms->natm_comp;
   int **ghost_atoms_iatm_comp = ghost_atoms->iatm_comp;

/*==========================================================================*/
/* I) Update forces due to ghost atoms and zero ghost atom forces         */

 if(nghost > 0) {
  for(ighost=1;ighost <= nghost;ighost++){
    igloc = ghost_atoms_ighost_map[ighost];
    for(icoef=1;icoef<=ghost_atoms_natm_comp[ighost];icoef++) {
      igpart = ghost_atoms_iatm_comp[icoef][ighost];
      clatoms_fx[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fx[igloc];
      clatoms_fy[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fy[igloc];
      clatoms_fz[igpart] += ghost_atoms_coef[icoef][ighost]*clatoms_fz[igloc];
    } /* endfor icoef */
    clatoms_fx[igloc] = 0.0;
    clatoms_fy[igloc] = 0.0;
    clatoms_fz[igloc] = 0.0;
  } /* endfor ighost */
 }/* endif */
/*-------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/






