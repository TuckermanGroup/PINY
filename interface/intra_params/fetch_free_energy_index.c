/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_res_parms.c                          */
/*                                                                          */
/* This subprogram reads in molecular parameter files and sets              */
/* molecular and intramolecular data sets                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_free_energy_index(BUILD_INTRA *build_intra,FREE_PARSE *free_parse,
                       BONDED *bonded,int jmol_typ,int iresidue,int natm_tot,
                       int natm_mol)

/*==========================================================================*/
/*      Begin Routine */
{ /* begin routine */
/*==========================================================================*/
/*      Local Variable Declarations */

  int itemp,i,iii;

  int bond_free_num       = bonded->bond_free.num;
  int bend_free_num       = bonded->bend_free.num;
  int tors_free_num       = bonded->tors_free.num;
  int rbar_sig_free_iopt  = bonded->rbar_sig_free.iopt;

  int *imoltyp_bond_free  = free_parse->imoltyp_bond_free;
  int *imol_bond_free     = free_parse->imol_bond_free;
  int *ires_bond_free     = free_parse->ires_bond_free;
  int *iatm_bond_free     = free_parse->iatm_bond_free;
  int *bond_free_j1       = &(bonded->bond_free.j1);
  int *bond_free_j2       = &(bonded->bond_free.j2);

  int *imoltyp_bend_free  = free_parse->imoltyp_bend_free;
  int *imol_bend_free     = free_parse->imol_bend_free;
  int *ires_bend_free     = free_parse->ires_bend_free;
  int *iatm_bend_free     = free_parse->iatm_bend_free;
  int *bend_free_j1       = &(bonded->bend_free.j1);
  int *bend_free_j2       = &(bonded->bend_free.j2);
  int *bend_free_j3       = &(bonded->bend_free.j3);

  int *imoltyp_tors_free  = free_parse->imoltyp_tors_free;
  int *imol_tors_free     = free_parse->imol_tors_free;
  int *ires_tors_free     = free_parse->ires_tors_free; 
  int *iatm_tors_free     = free_parse->iatm_tors_free;

  int *tors_free_j1       = bonded->tors_free.j1;
  int *tors_free_j2       = bonded->tors_free.j2;
  int *tors_free_j3       = bonded->tors_free.j3;
  int *tors_free_j4       = bonded->tors_free.j4;

  int *imoltyp_rbar1_free = free_parse->imoltyp_rbar1_free;
  int *imoltyp_rbar2_free = free_parse->imoltyp_rbar2_free;
  int *imol_rbar1_free    = free_parse->imol_rbar1_free;
  int *imol_rbar2_free    = free_parse->imol_rbar2_free;
  int *ires_rbar1_free    = free_parse->ires_rbar1_free;
  int *ires_rbar2_free    = free_parse->ires_rbar2_free;
  int *iatm_rbar1_free    = free_parse->iatm_rbar1_free;
  int *iatm_rbar2_free    = free_parse->iatm_rbar2_free;
  int nbar_bond           = free_parse->nbar_bond;

  int *rbar_sig_free_j1   = bonded->rbar_sig_free.j1;
  int *rbar_sig_free_j2   = bonded->rbar_sig_free.j2;

/*==========================================================================*/
/* I) Fix bonds */

  if(bond_free_num>0){

    if(imoltyp_bond_free[1]==jmol_typ){
      if(ires_bond_free[1]==iresidue){
        free_chk_fix(&itemp,iatm_bond_free[1],build_intra,jmol_typ,iresidue);
        (*bond_free_j1) = (imol_bond_free[1]-1)*natm_mol + itemp + natm_tot;
      }/*endif*/
    }/*endif*/
    if(imoltyp_bond_free[2]==jmol_typ){
      if(ires_bond_free[2]==iresidue){
        free_chk_fix(&itemp,iatm_bond_free[2],build_intra,jmol_typ,iresidue);
        (*bond_free_j2) = (imol_bond_free[2]-1)*natm_mol + itemp + natm_tot;
      }/*endif*/
   }/*endif*/

  }/*endif*/

/*==========================================================================*/
/* II) Fix bends */

  if(bend_free_num>0){

    if(imoltyp_bend_free[1]==jmol_typ){
      if(ires_bend_free[1]==iresidue){
        free_chk_fix(&itemp,iatm_bend_free[1],build_intra,jmol_typ,iresidue);
        (*bend_free_j1) = (imol_bend_free[1]-1)*natm_mol + itemp + natm_tot;
      }/*endif*/
    }/*endif*/
    if(imoltyp_bend_free[2]==jmol_typ){
      if(ires_bend_free[2]==iresidue){
        free_chk_fix(&itemp,iatm_bend_free[2],build_intra,jmol_typ,iresidue);
        (*bend_free_j2) = (imol_bend_free[2]-1)*natm_mol + itemp + natm_tot;
      }/*endif*/
    }/*endif*/
    if(imoltyp_bend_free[3]==jmol_typ){
      if(ires_bend_free[3]==iresidue){
        free_chk_fix(&itemp,iatm_bend_free[3],build_intra,jmol_typ,iresidue);
        (*bend_free_j3) = (imol_bend_free[3]-1)*natm_mol + itemp + natm_tot;
      }/*endif*/
    }/*endif*/

   }/*endif*/

/*==========================================================================*/
/* III) Fix tors */

  if(tors_free_num>0){

    if(imoltyp_tors_free[1]==jmol_typ){
       if(ires_tors_free[1]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[1],build_intra,jmol_typ,iresidue);
         tors_free_j1[1] = (imol_tors_free[1]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
    }/*endif*/
    if(imoltyp_tors_free[2]==jmol_typ){
       if(ires_tors_free[2]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[2],build_intra,jmol_typ,iresidue);
         tors_free_j2[1] = (imol_tors_free[2]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
    }/*endif*/
    if(imoltyp_tors_free[3]==jmol_typ){
       if(ires_tors_free[3]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[3],build_intra,jmol_typ,iresidue);
         tors_free_j3[1] = (imol_tors_free[3]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
    }/*endif*/
    if(imoltyp_tors_free[4]==jmol_typ){
       if(ires_tors_free[4]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[4],build_intra,jmol_typ,iresidue);
         tors_free_j4[1] = (imol_tors_free[4]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
    }/*endif*/

    if(tors_free_num==2){
      if(imoltyp_tors_free[5]==jmol_typ){
       if(ires_tors_free[5]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[5],build_intra,jmol_typ,iresidue);
         tors_free_j1[2] = (imol_tors_free[5]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
      }/*endif*/
      if(imoltyp_tors_free[6]==jmol_typ){
       if(ires_tors_free[6]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[6],build_intra,jmol_typ,iresidue);
         tors_free_j2[2] = (imol_tors_free[6]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
      }/*endif*/
      if(imoltyp_tors_free[7]==jmol_typ){
       if(ires_tors_free[7]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[7],build_intra,jmol_typ,iresidue);
         tors_free_j3[2] = (imol_tors_free[7]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
      }/*endif*/
      if(imoltyp_tors_free[8]==jmol_typ){
       if(ires_tors_free[8]==iresidue){
         free_chk_fix(&itemp,iatm_tors_free[8],build_intra,jmol_typ,iresidue);
         tors_free_j4[2] = (imol_tors_free[8]-1)*natm_mol + itemp + natm_tot;
       }/*endif*/
      }/*endif*/
    }/*endif : 2d*/

  }/*endif: free energy torsions*/

/*==========================================================================*/
/* IV) Fix rbar_sigma */

  if(rbar_sig_free_iopt>0){

   for(i=1;i<=nbar_bond;i++){

     if(imoltyp_rbar1_free[i]==jmol_typ){
       if(ires_rbar1_free[i]==iresidue){
         free_chk_fix(&itemp,iatm_rbar1_free[i],build_intra,jmol_typ,iresidue);
         rbar_sig_free_j1[i] = (imol_rbar1_free[i]-1)*natm_mol+itemp+natm_tot;
       }/*endif*/
     }/*endif*/
     if(imoltyp_rbar2_free[i]==jmol_typ){
       if(ires_rbar2_free[i]==iresidue){
         free_chk_fix(&itemp,iatm_rbar2_free[i],build_intra,jmol_typ,iresidue);
         rbar_sig_free_j2[i] = (imol_rbar2_free[i]-1)*natm_mol+itemp+natm_tot;
       }/*endif*/
     }/*endif*/

   }/*endfor*/

  }/*endif*/

/*==========================================================================*/
    }/*end routine */
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void free_chk_fix(int *itemp,int index,BUILD_INTRA *build_intra,
                  int jmol_typ, int iresidue)

/*==========================================================================*/
/*      Begin Routine */
{ /* begin routine */
/*==========================================================================*/
/*      Local Variable Declarations */

  int mask; 

  int natmind_1res_now = build_intra->natmind_1res_now;
  int *mask_atm        = build_intra->mask_atm;
  int *index_atm       = build_intra->index_atm;

/*==========================================================================*/
/* I) Check range of index */

  if(index>natmind_1res_now){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Free energy atom index out of range  \n");
      printf("in the residue number %d of molecule number %d\n",
              iresidue,jmol_typ);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

/*==========================================================================*/
/* II) Get and check mask of index */

  mask = mask_atm[index];
  if(mask==0){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Free energy atom has been nuked \n");
      printf("in residue number %d of molecule number %d \n",
              iresidue,jmol_typ);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
   }/*endif*/

/*==========================================================================*/
/* III) Get rejiggered index       */
 
  *itemp = index_atm[index];

/*==========================================================================*/
    }/*end routine */
/*==========================================================================*/







