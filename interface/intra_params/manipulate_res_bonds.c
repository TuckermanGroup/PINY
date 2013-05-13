/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: man_mol_parms.c                              */
/*                                                                          */
/* Manipulate molecular parameters                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void map_residue_bonds(RESBOND_PARSE *resbond_parse)

/*=======================================================================*/
{ /*begin routine*/
/*=======================================================================*/
/*             Local variable declarations                               */

  int i;
  int indx1,indx2,itemp1,itemp2;
  
/*========================================================================*/
/* I) Find the number of resbonds each residue is involved in             */

  for(i=1;i<=resbond_parse->nresidue;i++){
    resbond_parse->nres_bond1_jres[i]=0;
    resbond_parse->nres_bond2_jres[i]=0;
  }/*endfor*/

  for(i=1;i<=resbond_parse->nres_bond;i++){
   indx1 = resbond_parse->resbond_prm[i].res1_index;
   indx2 = resbond_parse->resbond_prm[i].res2_index;
   resbond_parse->nres_bond1_jres[indx1] += 1;
   resbond_parse->nres_bond2_jres[indx2] += 1;  
  }/*endfor*/

/*========================================================================*/
/* II) Find where in the residue bonding list the residue bonds of the    */
/*     jth residue will start                                             */

  resbond_parse->res_bond1_off[1]=0; 
  resbond_parse->res_bond2_off[1]=0;
  for(i=1;i<=resbond_parse->nresidue-1;i++){
    resbond_parse->res_bond1_off[(i+1)]=
      (resbond_parse->res_bond1_off[i]+resbond_parse->nres_bond1_jres[i]);
    resbond_parse->res_bond2_off[(i+1)]=
      (resbond_parse->res_bond2_off[i]+resbond_parse->nres_bond2_jres[i]);
  }/*endfor*/

/*=========================================================================*/
/* II) Create the residue bonding list: resbond_parse->res_bondk_index     */
/*                                                                         */
/*    It stores the indicies of the residue bonds in which the jth         */ 
/*    residue is involved. The list for the jth residue starts at          */
/*    resbond_parse->res_bondk_off[j]+1 and is completed at                */
/*    resbond_parse->res_bondk_off[j]+resbond_parse->nres_bondk_jres[j]    */

  for(i=1;i<=resbond_parse->nresidue;i++){
    resbond_parse->nres_bond1_jres[i]=0;
    resbond_parse->nres_bond2_jres[i]=0;
  }/*endfor*/
  
  for(i=1;i<=resbond_parse->nres_bond;i++){
    indx1 = resbond_parse->resbond_prm[i].res1_index;
    indx2 = resbond_parse->resbond_prm[i].res2_index;
    resbond_parse->nres_bond1_jres[indx1] += 1;
    resbond_parse->nres_bond2_jres[indx2] += 1;
    itemp1=(resbond_parse->nres_bond1_jres[indx1]
          + resbond_parse->res_bond1_off[indx1]);
    itemp2=(resbond_parse->nres_bond2_jres[indx2]
          + resbond_parse->res_bond2_off[indx2]);
    resbond_parse->res_bond1_index[itemp1] = i;
    resbond_parse->res_bond2_index[itemp2] = i;
  }/*endfor*/

/*========================================================================*/
}  /*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void resbond_parse_realloc(RESBOND_PARSE *resbond_parse)

/*=======================================================================*/
{ /*begin routine*/
/*=======================================================================*/
/* I) Realloc condition 1: Increased residue bonds                      */

  if(resbond_parse->nres_bond > resbond_parse->nres_bond_max){
     resbond_parse->nres_bond_max   += 
    MAX(NMEM_MIN,resbond_parse->nres_bond-resbond_parse->nres_bond_max);
     resbond_parse->resbond_prm     = (RESBOND_PRM *)
                    crealloc(&(resbond_parse->resbond_prm[1]),
                              (resbond_parse->nres_bond_max)
                              *sizeof(RESBOND_PRM))-1;
     resbond_parse->res_bond1_index = (int *)
                      crealloc(&(resbond_parse->res_bond1_index[1]),
                              (resbond_parse->nres_bond_max)*sizeof(int))-1; 
     resbond_parse->res_bond2_index = (int *)
                      crealloc(&(resbond_parse->res_bond2_index[1]),
                              (resbond_parse->nres_bond_max)*sizeof(int))-1; 
  }/*endif*/

/*=======================================================================*/
/* I) Realloc condition 2: Increased residues                            */

  if(resbond_parse->nresidue > resbond_parse->nresidue_max){
     resbond_parse->nresidue_max   += 
         MAX(NMEM_MIN,resbond_parse->nresidue-resbond_parse->nresidue_max);
     resbond_parse->res_bond1_off   = (int *)
                      crealloc(&(resbond_parse->res_bond1_off[1]),
                                (resbond_parse->nresidue_max)*sizeof(int))-1; 
     resbond_parse->res_bond2_off   = (int *)
                      crealloc(&(resbond_parse->res_bond2_off[1]),
                                (resbond_parse->nresidue_max)*sizeof(int))-1; 
     resbond_parse->nres_bond1_jres = (int *)
                      crealloc(&(resbond_parse->nres_bond1_jres[1]),
                                (resbond_parse->nresidue_max)*sizeof(int))-1; 
     resbond_parse->nres_bond2_jres = (int *)
                      crealloc(&(resbond_parse->nres_bond2_jres[1]),
                                (resbond_parse->nresidue_max)*sizeof(int))-1; 
  }/*endif*/
}/*end routine*/ 
/*==========================================================================*/




