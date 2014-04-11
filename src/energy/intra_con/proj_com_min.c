/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: proj_com_min                                 */
/*                                                                          */
/* This subprogram projects force along the com                             */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_intra_con_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void proj_com_out(int natm_tot, double *fx,double *fy, double *fz)

/*========================================================================*/
   {/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

    int i;
    double xcm,ycm,zcm;

/*==========================================================================*/
/* I) Get the force on the com */

   xcm = 0.0;
   ycm = 0.0;
   zcm = 0.0;
   for(i=1;i<=natm_tot;i++){
     xcm += fx[i];
     ycm += fy[i];
     zcm += fz[i];
   }/*endfor*/
   xcm /= (double)natm_tot;
   ycm /= (double)natm_tot;
   zcm /= (double)natm_tot;

/*==========================================================================*/
/* II) Remove it */

   for(i=1;i<=natm_tot;i++){
     fx[i] -= xcm;
     fy[i] -= ycm;
     fz[i] -= zcm;
   }/*endfor*/

/*--------------------------------------------------------------------------*/
  }/*end routine*/
/*==========================================================================*/
