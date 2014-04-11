/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: fetch_hydrog_mass.c                          */
/*                                                                          */
/* This subprogram sets the masses of selected hydrogen atoms               */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fetch_hydrog_mass(CLASS_PARSE *class_parse, ATOMMAPS *atommaps, 
                  BUILD_INTRA *build_intra, START_INDEX *start_index, 
                  CLATOMS_INFO *clatoms_info, int jmol_typ)

/*==========================================================================*/
/*begin routine*/{

/*==========================================================================*/
/* Define Local variable */

    int i,iii,nfreeze_now,iresidue_off,nresidue,count;
    int count_atm,ires;
    int mass_now;

/*==========================================================================*/
/* I) All hydrogen masses changed */ 
  if(class_parse->mol_hydrog_mass_opt[jmol_typ]==1){
       for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        mass_now = (int)(clatoms_info->mass[i]);
        if(mass_now<=2){
          clatoms_info->mass[i] = (class_parse->mol_hydrog_mass_val[jmol_typ]);
        }/*endif*/
      }/*endfor*/
   }/*endif*/

/*==========================================================================*/
/* II) Backbone hydrogen mass changes */ 

  if(class_parse->mol_hydrog_mass_opt[jmol_typ]==2){
     for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        mass_now = (int)(clatoms_info->mass[i]);
       if((atommaps->atom_label[i] == 1) && (mass_now<=2)){
         clatoms_info->mass[i] = (class_parse->mol_hydrog_mass_val[jmol_typ]);
       }/*endif*/
      }/*endfor*/
  }/*endif*/

/*==========================================================================*/
/* II) Sidechain hydrogen mass changes */

  if(class_parse->mol_hydrog_mass_opt[jmol_typ]==3){
     for(i=start_index->natm+1;i<=clatoms_info->natm_tot;i++){
        mass_now = (int)(clatoms_info->mass[i]);
       if((atommaps->atom_label[i] == 2) && (mass_now<=2)){
         clatoms_info->mass[i] = (class_parse->mol_hydrog_mass_val[jmol_typ]);
       }/*endif*/
      }/*endfor*/
  }/*endif*/

/*--------------------------------------------------------------------------*/
}/*end routine*/ 
/*==========================================================================*/
