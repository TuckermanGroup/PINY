/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: Build_map.c                                  */
/*                                                                          */
/* This subprogram builds some atoms maps                                   */
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

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void replicate_mol(CLATOMS_INFO *clatoms_info,GHOST_ATOMS *ghost_atoms,
                   ATOMMAPS *atommaps,
                   BUILD_INTRA *build_intra,BONDED *bonded,
                   NULL_INTER_PARSE *null_inter_parse,
		   START_INDEX *start_index,int jmol_typ)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/

/*========================================================================*/
/*     Local Variables                                                    */
   
   int imol,i,ishift,nmol,nmol1,k,iii,ig_shift,if_shift;
   int nbond_pow_add,nbond_con_add,nbond_null_add;
   int nbend_pow_add,nbend_con_add,nbend_null_add;
   int ntors_pow_add,ntors_con_add,ntors_null_add;
   int nonfo_add,nonfo_null_add,nbend_bnd_add;
   int natm_add;
   int nghost_add;
   int nfreeze_add;
   int ngrp_21_add; 
   int ngrp_33_add; 
   int ngrp_watt_33_add; 
   int ngrp_43_add; 
   int ngrp_23_add;
   int ngrp_46_add;
   int ibond_pow_off,ibond_con_off,ibond_null_off;
   int ibend_pow_off,ibend_con_off,ibend_null_off;
   int itors_pow_off,itors_con_off,itors_null_off;
   int ionfo_off,ionfo_null_off,ibend_bnd_off;
   int iatm_off;
   int ighost_off;
   int ifreeze_off;
   int igrp_21_off; 
   int igrp_33_off; 
   int igrp_watt_33_off; 
   int igrp_43_off; 
   int igrp_23_off;
   int igrp_46_off;
   int ibond_pow_offn,ibond_con_offn,ibond_null_offn;
   int ibend_pow_offn,ibend_con_offn,ibend_null_offn;
   int itors_pow_offn,itors_con_offn,itors_null_offn;
   int ionfo_offn,ionfo_null_offn,ibend_bnd_offn;
   int iatm_offn;
   int ighost_offn;
   int ifreeze_offn;
   int igrp_21_offn; 
   int igrp_33_offn; 
   int igrp_watt_33_offn; 
   int igrp_43_offn; 
   int igrp_23_offn;
   int igrp_46_offn;

/*========================================================================*/
/* I) Calculate how many of what to replicate                             */

   nbond_pow_add = bonded->bond.npow            - start_index->nbond_pow;
   nbond_con_add = bonded->bond.ncon            - start_index->nbond_con;
   nbond_null_add = null_inter_parse->nbond_nul - start_index->nbond_nul;
   nbend_pow_add = bonded->bend.npow            - start_index->nbend_pow;
   nbend_con_add = bonded->bend.ncon            - start_index->nbend_con;
   nbend_null_add = null_inter_parse->nbend_nul - start_index->nbend_nul;
   ntors_pow_add = bonded->tors.npow            - start_index->ntors_pow;
   ntors_con_add = bonded->tors.ncon            - start_index->ntors_con;
   ntors_null_add = null_inter_parse->ntors_nul - start_index->ntors_nul;
   nonfo_add     = bonded->onfo.num             - start_index->nonfo;
   nonfo_null_add = null_inter_parse->nonfo_nul - start_index->nonfo_nul;
   nbend_bnd_add = bonded->bend_bnd.num         - start_index->nbend_bnd;
   natm_add      = clatoms_info->natm_tot       - start_index->natm;
   nghost_add    = ghost_atoms->nghost_tot      - start_index->nghost_tot;
   nfreeze_add    = atommaps->nfreeze      - start_index->nfreeze;
   ngrp_21_add   = bonded->grp_bond_con.num_21         - start_index->ngrp_21;
   ngrp_23_add   = bonded->grp_bond_con.num_23         - start_index->ngrp_23;
   ngrp_43_add   = bonded->grp_bond_con.num_43         - start_index->ngrp_43;
   ngrp_33_add   = bonded->grp_bond_con.num_33         - start_index->ngrp_33;
   ngrp_watt_33_add   = bonded->grp_bond_watts.num_33  
                                     - start_index->ngrp_watt_33;
   ngrp_46_add   = bonded->grp_bond_con.num_46         - start_index->ngrp_46;


/*========================================================================*/
/* I) Rellocate                                                 */

   reallocate_intra_list(clatoms_info,ghost_atoms,
                         atommaps,build_intra,bonded,null_inter_parse,
			 jmol_typ,nbond_pow_add, nbond_con_add, nbond_null_add,
			 nbend_pow_add, nbend_con_add, nbend_null_add,
			 ntors_pow_add, ntors_con_add, ntors_null_add,
			 nonfo_add, nonfo_null_add, nbend_bnd_add,
		         natm_add,nghost_add, nfreeze_add,ngrp_43_add,
                         ngrp_33_add,ngrp_watt_33_add,
                         ngrp_21_add,ngrp_23_add, ngrp_46_add);

/*========================================================================*/
/* II) Store the offsets                                                 */

   ibond_pow_off  = start_index->nbond_pow;
   ibond_con_off  = start_index->nbond_con;
   ibond_null_off = start_index->nbond_nul;
   ibend_pow_off  = start_index->nbend_pow;
   ibend_con_off  = start_index->nbend_con;
   ibend_null_off = start_index->nbend_nul;
   itors_pow_off  = start_index->ntors_pow;
   itors_con_off  = start_index->ntors_con;
   itors_null_off = start_index->ntors_nul;
   ionfo_off      = start_index->nonfo;
   ionfo_null_off = start_index->nonfo_nul;
   ibend_bnd_off  = start_index->nbend_bnd;
   iatm_off       = start_index->natm;
   ighost_off     = start_index->nghost_tot;
   ifreeze_off    = start_index->nfreeze;
   igrp_21_off    = start_index->ngrp_21;
   igrp_23_off    = start_index->ngrp_23;
   igrp_33_off    = start_index->ngrp_33;
   igrp_watt_33_off    = start_index->ngrp_watt_33;
   igrp_43_off    = start_index->ngrp_43;
   igrp_46_off    = start_index->ngrp_46;

   ibond_pow_offn  = start_index->nbond_pow;
   ibond_con_offn  = start_index->nbond_con;
   ibond_null_offn = start_index->nbond_nul;
   ibend_pow_offn  = start_index->nbend_pow;
   ibend_con_offn  = start_index->nbend_con;
   ibend_null_offn = start_index->nbend_nul;
   itors_pow_offn  = start_index->ntors_pow;
   itors_con_offn  = start_index->ntors_con;
   itors_null_offn = start_index->ntors_nul;
   ionfo_offn      = start_index->nonfo;
   ionfo_null_offn = start_index->nonfo_nul;
   ibend_bnd_offn  = start_index->nbend_bnd;
   iatm_offn       = start_index->natm;
   ighost_offn     = start_index->nghost_tot;
   ifreeze_offn     = start_index->nfreeze;
   igrp_21_offn    = start_index->ngrp_21;
   igrp_23_offn    = start_index->ngrp_23;
   igrp_33_offn    = start_index->ngrp_33;
   igrp_watt_33_offn    = start_index->ngrp_watt_33;
   igrp_43_offn    = start_index->ngrp_43;
   igrp_46_offn    = start_index->ngrp_46;

/*========================================================================*/
/* III) Replicate the bonds                                               */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     ibond_pow_offn += nbond_pow_add;
     for(i=1;i<=nbond_pow_add;i++){
       bonded->bond.j1_pow[(i+ibond_pow_offn)] = 
                       bonded->bond.j1_pow[(i+ibond_pow_off)]+ishift;
       bonded->bond.j2_pow[(i+ibond_pow_offn)] = 
                       bonded->bond.j2_pow[(i+ibond_pow_off)]+ishift;
       bonded->bond.jtyp_pow[(i+ibond_pow_offn)] = 
                       bonded->bond.jtyp_pow[(i+ibond_pow_off)];
     }/*endfor*/
     ibond_con_offn += nbond_con_add;
     for(i=1;i<=nbond_con_add;i++){
       bonded->bond.j1_con[(i+ibond_con_offn)] = 
                       bonded->bond.j1_con[(i+ibond_con_off)]+ishift;
       bonded->bond.j2_con[(i+ibond_con_offn)] = 
                       bonded->bond.j2_con[(i+ibond_con_off)]+ishift;
       bonded->bond.jtyp_con[(i+ibond_con_offn)] = 
                       bonded->bond.jtyp_con[(i+ibond_con_off)];
     }/*endfor*/
     ibond_null_offn += nbond_null_add;
     for(i=1;i<=nbond_null_add;i++){
       null_inter_parse->jbond1_nul[(i+ibond_null_offn)] = 
                       null_inter_parse->jbond1_nul[(i+ibond_null_off)]+ishift;
       null_inter_parse->jbond2_nul[(i+ibond_null_offn)] = 
                       null_inter_parse->jbond2_nul[(i+ibond_null_off)]+ishift;
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* IV) Replicate the bends                                                */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     ibend_pow_offn += nbend_pow_add;
     for(i=1;i<=nbend_pow_add;i++){
       bonded->bend.j1_pow[(i+ibend_pow_offn)] = 
               bonded->bend.j1_pow[(i+ibend_pow_off)]+ishift;
       bonded->bend.j2_pow[(i+ibend_pow_offn)] = 
               bonded->bend.j2_pow[(i+ibend_pow_off)]+ishift;
       bonded->bend.j3_pow[(i+ibend_pow_offn)] = 
               bonded->bend.j3_pow[(i+ibend_pow_off)]+ishift;
       bonded->bend.jtyp_pow[(i+ibend_pow_offn)] = 
                       bonded->bend.jtyp_pow[(i+ibend_pow_off)];
     }/*endfor*/
     ibend_con_offn += nbend_con_add;
     for(i=1;i<=nbend_con_add;i++){
       bonded->bend.j1_con[(i+ibend_con_offn)] = 
                       bonded->bend.j1_con[(i+ibend_con_off)]+ishift;
       bonded->bend.j2_con[(i+ibend_con_offn)] = 
                       bonded->bend.j2_con[(i+ibend_con_off)]+ishift;
       bonded->bend.j3_con[(i+ibend_con_offn)] = 
                       bonded->bend.j3_con[(i+ibend_con_off)]+ishift;
       bonded->bend.jtyp_con[(i+ibend_con_offn)] = 
                       bonded->bend.jtyp_con[(i+ibend_con_off)];
     }/*endfor*/
     ibend_null_offn += nbend_null_add;
     for(i=1;i<=nbend_null_add;i++){
       null_inter_parse->jbend1_nul[(i+ibend_null_offn)] = 
                       null_inter_parse->jbend1_nul[(i+ibend_null_off)]+ishift;
       null_inter_parse->jbend2_nul[(i+ibend_null_offn)] = 
                       null_inter_parse->jbend2_nul[(i+ibend_null_off)]+ishift;
       null_inter_parse->jbend3_nul[(i+ibend_null_offn)] = 
                       null_inter_parse->jbend3_nul[(i+ibend_null_off)]+ishift;
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* V) Replicate the torsions                                              */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     itors_pow_offn += ntors_pow_add;
     for(i=1;i<=ntors_pow_add;i++){
       bonded->tors.j1_pow[(i+itors_pow_offn)] = 
               bonded->tors.j1_pow[(i+itors_pow_off)]+ishift;
       bonded->tors.j2_pow[(i+itors_pow_offn)] = 
               bonded->tors.j2_pow[(i+itors_pow_off)]+ishift;
       bonded->tors.j3_pow[(i+itors_pow_offn)] = 
               bonded->tors.j3_pow[(i+itors_pow_off)]+ishift;
       bonded->tors.j4_pow[(i+itors_pow_offn)] = 
               bonded->tors.j4_pow[(i+itors_pow_off)]+ishift;
       bonded->tors.jtyp_pow[(i+itors_pow_offn)] = 
                       bonded->tors.jtyp_pow[(i+itors_pow_off)];
     }/*endfor*/
     itors_con_offn += ntors_con_add;
     for(i=1;i<=ntors_con_add;i++){
       bonded->tors.j1_con[(i+itors_con_offn)] = 
                       bonded->tors.j1_con[(i+itors_con_off)]+ishift;
       bonded->tors.j2_con[(i+itors_con_offn)] = 
                       bonded->tors.j2_con[(i+itors_con_off)]+ishift;
       bonded->tors.j3_con[(i+itors_con_offn)] = 
                       bonded->tors.j3_con[(i+itors_con_off)]+ishift;
       bonded->tors.j4_con[(i+itors_con_offn)] = 
                       bonded->tors.j4_con[(i+itors_con_off)]+ishift;
       bonded->tors.jtyp_con[(i+itors_con_offn)] = 
                       bonded->tors.jtyp_con[(i+itors_con_off)];
     }/*endfor*/
     itors_null_offn += ntors_null_add;
     for(i=1;i<=ntors_null_add;i++){
       null_inter_parse->jtors1_nul[(i+itors_null_offn)] = 
                       null_inter_parse->jtors1_nul[(i+itors_null_off)]+ishift;
       null_inter_parse->jtors2_nul[(i+itors_null_offn)] = 
                       null_inter_parse->jtors2_nul[(i+itors_null_off)]+ishift;
       null_inter_parse->jtors3_nul[(i+itors_null_offn)] = 
                       null_inter_parse->jtors3_nul[(i+itors_null_off)]+ishift;
       null_inter_parse->jtors4_nul[(i+itors_null_offn)] = 
                       null_inter_parse->jtors4_nul[(i+itors_null_off)]+ishift;
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* VI) Replicate the onfos                                                */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     ionfo_offn += nonfo_add;
     for(i=1;i<=nonfo_add    ;i++){
       bonded->onfo.j1[(i+ionfo_offn)] = 
               bonded->onfo.j1[(i+ionfo_off)]+ishift;
       bonded->onfo.j2[(i+ionfo_offn)] = 
               bonded->onfo.j2[(i+ionfo_off)]+ishift;
       bonded->onfo.jtyp[(i+ionfo_offn)] = 
                       bonded->onfo.jtyp[(i+ionfo_off)];
     }/*endfor*/
     ionfo_null_offn += nonfo_null_add;
     for(i=1;i<=nonfo_null_add;i++){
       null_inter_parse->jonfo1_nul[(i+ionfo_null_offn)] = 
                       null_inter_parse->jonfo1_nul[(i+ionfo_null_off)]+ishift;
       null_inter_parse->jonfo2_nul[(i+ionfo_null_offn)] = 
                       null_inter_parse->jonfo2_nul[(i+ionfo_null_off)]+ishift;
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* VII) Replicate the bend_bnds                                           */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     ibend_bnd_offn += nbend_bnd_add;
     for(i=1;i<=nbend_bnd_add;i++){
       bonded->bend_bnd.j1[(i+ibend_bnd_offn)] = 
               bonded->bend_bnd.j1[(i+ibend_bnd_off)]+ishift;
       bonded->bend_bnd.j2[(i+ibend_bnd_offn)] = 
               bonded->bend_bnd.j2[(i+ibend_bnd_off)]+ishift;
       bonded->bend_bnd.j3[(i+ibend_bnd_offn)] = 
               bonded->bend_bnd.j3[(i+ibend_bnd_off)]+ishift;
       bonded->bend_bnd.jtyp[(i+ibend_bnd_offn)] = 
                       bonded->bend_bnd.jtyp[(i+ibend_bnd_off)];
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* VIII) Replicate the atoms                                              */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     iatm_offn += natm_add;
     ishift = imol*natm_add;
     ig_shift = (nghost_add)*imol;
     if_shift = (nfreeze_add)*imol;
     for(i=1;i<=natm_add     ;i++){
      clatoms_info->mass[(i+iatm_offn)]   = clatoms_info->mass[(i+iatm_off)];
      clatoms_info->q[(i+iatm_offn)]      = clatoms_info->q[(i+iatm_off)];
      clatoms_info->cp_vlnc_up[(i+iatm_offn)]= clatoms_info->cp_vlnc_up[(i+iatm_off)];
      clatoms_info->cp_vlnc_dn[(i+iatm_offn)]= clatoms_info->cp_vlnc_dn[(i+iatm_off)];
      clatoms_info->cp_atm_flag[(i+iatm_offn)]= clatoms_info->cp_atm_flag[(i+iatm_off)];
      clatoms_info->alp_pol[(i+iatm_offn)]= clatoms_info->alp_pol[(i+iatm_off)];
      clatoms_info->b_neut[(i+iatm_offn)] = clatoms_info->b_neut[(i+iatm_off)];
      clatoms_info->text_atm[(i+iatm_offn)] = 
                                     clatoms_info->text_atm[(i+iatm_off)];
      atommaps->iatm_atm_typ[(i+iatm_offn)] = 
                 atommaps->iatm_atm_typ[(i+iatm_off)];
      atommaps->iatm_mol_typ[(i+iatm_offn)] = jmol_typ;
      atommaps->iatm_mol_num[(i+iatm_offn)] = imol+1;
      atommaps->iatm_res_num[(i+iatm_offn)] = 
                                   atommaps->iatm_res_num[(i+iatm_off)];
      atommaps->iatm_res_typ[(i+iatm_offn)] = 
                                   atommaps->iatm_res_typ[(i+iatm_off)];
      atommaps->atom_label[(i+iatm_offn)] = 
                                   atommaps->atom_label[(i+iatm_off)];
      atommaps->freeze_flag[(i+iatm_offn)]  =  0;
      if(atommaps->freeze_flag[(i+iatm_off)] != 0) {
      atommaps->freeze_flag[(i+iatm_offn)]  =  
                              atommaps->freeze_flag[(i+iatm_off)]+if_shift;}
      atommaps->ighost_flag[(i+iatm_offn)]  =  0;
      if(atommaps->ighost_flag[(i+iatm_off)] != 0) {
      atommaps->ighost_flag[(i+iatm_offn)]  = 
                                atommaps->ighost_flag[(i+iatm_off)]+ig_shift;}
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* IX) Replicate the ghosts                                              */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ighost_offn += nghost_add;
     ishift = imol*natm_add;
     for(i=1;i<=nghost_add     ;i++){
       ghost_atoms->ighost_map[(i+ighost_offn)] = 
                             ghost_atoms->ighost_map[(i+ighost_off)] + ishift;
       ghost_atoms->natm_comp[(i+ighost_offn)] = 
                             ghost_atoms->natm_comp[(i+ighost_off)];
       for(k=1;k<=ghost_atoms->natm_comp[(i+ighost_off)];k++){
         ghost_atoms->iatm_comp[k][(i+ighost_offn)] = 
                          ghost_atoms->iatm_comp[k][(i+ighost_off)]+ishift;
         ghost_atoms->coef[k][(i+ighost_offn)] = 
                          ghost_atoms->coef[k][(i+ighost_off)];
       }/*endfor*/
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* IX) Replicate the frozen atoms                                         */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ifreeze_offn += nfreeze_add;
     ishift = imol*natm_add;
     for(i=1;i<=nfreeze_add     ;i++){
       atommaps->freeze_map[(i+ifreeze_offn)] = 
                             atommaps->freeze_map[(i+ifreeze_off)] + ishift;
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* X) Replicate the grp_cons                                              */

   nmol = atommaps->nmol_jmol_typ[jmol_typ];
   for(imol=1;imol<=(nmol-1);imol++){
     ishift = imol*natm_add;
     igrp_21_offn += ngrp_21_add;
     for(i=1;i<=ngrp_21_add;i++){
       bonded->grp_bond_con.j1_21[(i+igrp_21_offn)] = 
               bonded->grp_bond_con.j1_21[(i+igrp_21_off)]+ishift;
       bonded->grp_bond_con.j2_21[(i+igrp_21_offn)] = 
               bonded->grp_bond_con.j2_21[(i+igrp_21_off)]+ishift;
       bonded->grp_bond_con.jtyp_21[(i+igrp_21_offn)] = 
                       bonded->grp_bond_con.jtyp_21[(i+igrp_21_off)];
     }/*endfor*/
     igrp_23_offn += ngrp_23_add;
     for(i=1;i<=ngrp_23_add;i++){
       bonded->grp_bond_con.j1_23[(i+igrp_23_offn)] = 
               bonded->grp_bond_con.j1_23[(i+igrp_23_off)]+ishift;
       bonded->grp_bond_con.j2_23[(i+igrp_23_offn)] = 
               bonded->grp_bond_con.j2_23[(i+igrp_23_off)]+ishift;
       bonded->grp_bond_con.j3_23[(i+igrp_23_offn)] = 
               bonded->grp_bond_con.j3_23[(i+igrp_23_off)]+ishift;
       bonded->grp_bond_con.jtyp_23[(i+igrp_23_offn)] = 
                       bonded->grp_bond_con.jtyp_23[(i+igrp_23_off)];
     }/*endfor*/
     igrp_33_offn += ngrp_33_add;
     for(i=1;i<=ngrp_33_add;i++){
       bonded->grp_bond_con.j1_33[(i+igrp_33_offn)] = 
               bonded->grp_bond_con.j1_33[(i+igrp_33_off)]+ishift;
       bonded->grp_bond_con.j2_33[(i+igrp_33_offn)] = 
               bonded->grp_bond_con.j2_33[(i+igrp_33_off)]+ishift;
       bonded->grp_bond_con.j3_33[(i+igrp_33_offn)] = 
               bonded->grp_bond_con.j3_33[(i+igrp_33_off)]+ishift;
       bonded->grp_bond_con.jtyp_33[(i+igrp_33_offn)] = 
                       bonded->grp_bond_con.jtyp_33[(i+igrp_33_off)];
     }/*endfor*/
     igrp_watt_33_offn += ngrp_watt_33_add;
     for(i=1;i<=ngrp_watt_33_add;i++){
       bonded->grp_bond_watts.j1_33[(i+igrp_watt_33_offn)] = 
               bonded->grp_bond_watts.j1_33[(i+igrp_watt_33_off)]+ishift;
       bonded->grp_bond_watts.j2_33[(i+igrp_watt_33_offn)] = 
               bonded->grp_bond_watts.j2_33[(i+igrp_watt_33_off)]+ishift;
       bonded->grp_bond_watts.j3_33[(i+igrp_watt_33_offn)] = 
               bonded->grp_bond_watts.j3_33[(i+igrp_watt_33_off)]+ishift;
       bonded->grp_bond_watts.jtyp_33[(i+igrp_watt_33_offn)] = 
                       bonded->grp_bond_watts.jtyp_33[(i+igrp_watt_33_off)];
     }/*endfor*/
     igrp_43_offn += ngrp_43_add;
     for(i=1;i<=ngrp_43_add;i++){
       bonded->grp_bond_con.j1_43[(i+igrp_43_offn)] = 
               bonded->grp_bond_con.j1_43[(i+igrp_43_off)]+ishift;
       bonded->grp_bond_con.j2_43[(i+igrp_43_offn)] = 
               bonded->grp_bond_con.j2_43[(i+igrp_43_off)]+ishift;
       bonded->grp_bond_con.j3_43[(i+igrp_43_offn)] = 
               bonded->grp_bond_con.j3_43[(i+igrp_43_off)]+ishift;
       bonded->grp_bond_con.j4_43[(i+igrp_43_offn)] = 
               bonded->grp_bond_con.j4_43[(i+igrp_43_off)]+ishift;
       bonded->grp_bond_con.jtyp_43[(i+igrp_43_offn)] = 
                       bonded->grp_bond_con.jtyp_43[(i+igrp_43_off)];
     }/*endfor*/
     igrp_46_offn += ngrp_46_add;
     for(i=1;i<=ngrp_46_add;i++){
       bonded->grp_bond_con.j1_46[(i+igrp_46_offn)] = 
               bonded->grp_bond_con.j1_46[(i+igrp_46_off)]+ishift;
       bonded->grp_bond_con.j2_46[(i+igrp_46_offn)] = 
               bonded->grp_bond_con.j2_46[(i+igrp_46_off)]+ishift;
       bonded->grp_bond_con.j3_46[(i+igrp_46_offn)] = 
               bonded->grp_bond_con.j3_46[(i+igrp_46_off)]+ishift;
       bonded->grp_bond_con.j4_46[(i+igrp_46_offn)] = 
               bonded->grp_bond_con.j4_46[(i+igrp_46_off)]+ishift;
       bonded->grp_bond_con.jtyp_46[(i+igrp_46_offn)] = 
                       bonded->grp_bond_con.jtyp_46[(i+igrp_46_off)];
     }/*endfor*/
   }/*endfor*/

/*========================================================================*/
/* XI) Increase the list counters                                         */

   nmol1 = nmol-1;


   bonded->bond.npow += nbond_pow_add*nmol1;
   bonded->bond.ncon += nbond_con_add*nmol1;
   null_inter_parse->nbond_nul += nbond_null_add*nmol1;

   bonded->bend.npow += nbend_pow_add*nmol1;
   bonded->bend.ncon += nbend_con_add*nmol1;
   null_inter_parse->nbend_nul += nbend_null_add*nmol1;

   bonded->bend_bnd.num += nbend_bnd_add*nmol1;

   bonded->tors.npow += ntors_pow_add*nmol1;
   bonded->tors.ncon += ntors_con_add*nmol1;
   null_inter_parse->ntors_nul += ntors_null_add*nmol1;

   bonded->onfo.num += nonfo_add*nmol1;
   null_inter_parse->nonfo_nul += nonfo_null_add*nmol1;

   clatoms_info->natm_tot += natm_add*nmol1;
   ghost_atoms->nghost_tot += nghost_add*nmol1;
   atommaps->nfreeze += nfreeze_add*nmol1;

   bonded->grp_bond_con.num_21 += ngrp_21_add*nmol1;
   bonded->grp_bond_con.num_23 += ngrp_23_add*nmol1;
   bonded->grp_bond_con.num_33 += ngrp_33_add*nmol1;
   bonded->grp_bond_watts.num_33 += ngrp_watt_33_add*nmol1;
   bonded->grp_bond_con.num_43 += ngrp_43_add*nmol1;   
   bonded->grp_bond_con.num_46 += ngrp_46_add*nmol1;   


/*==========================================================================*/
}/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void reallocate_intra_list(CLATOMS_INFO *clatoms_info,GHOST_ATOMS *ghost_atoms,
                  ATOMMAPS *atommaps,
                  BUILD_INTRA *build_intra,BONDED *bonded,
		  NULL_INTER_PARSE *null_inter_parse,
                  int jmol_typ,
		  int nbond_pow_add,int nbond_con_add,int nbond_null_add,
		  int nbend_pow_add,int nbend_con_add,int nbend_null_add,
		  int ntors_pow_add,int ntors_con_add,int ntors_null_add,
		  int nonfo_add,int nonfo_null_add,int nbend_bnd_add,
		  int natm_add,int nghost_add,int nfreeze_add,int ngrp_43_add,
		  int ngrp_33_add,int ngrp_watt_33_add,int ngrp_21_add,
                  int ngrp_23_add,int ngrp_46_add)

/*==========================================================================*/
/*        Begin routine  */
 { /*begin routine */
/*==========================================================================*/
  int nmol,nmol1,mem_add_now,iii;
  int nfreeze_old,nfreeze_new;
  int nghost_new,nghost_old,ncomp_new,ncomp_old;
/*==========================================================================*/
/* I) Reallocate the atoms                                                  */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*natm_add;
  if(clatoms_info->natm_tot+nmol1 > build_intra->natm_tot_max){
     mem_add_now = clatoms_info->natm_tot+nmol1 - build_intra->natm_tot_max;
     build_intra->natm_tot_max += MAX(mem_add_now,NMEM_MIN);

      clatoms_info->mass = (double *)crealloc(&((clatoms_info->mass[1])),
                                         build_intra->natm_tot_max*
                                         sizeof(double))-1;
      clatoms_info->q        = (double *)crealloc(&((clatoms_info->q[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(double))-1;
      clatoms_info->cp_vlnc_up  = (int *)crealloc(&((clatoms_info->cp_vlnc_up[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(int))-1;
      clatoms_info->cp_vlnc_dn  = (int *)crealloc(&((clatoms_info->cp_vlnc_dn[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(int))-1;
      clatoms_info->cp_atm_flag  = (int *)crealloc(&((clatoms_info->cp_atm_flag[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(int))-1;
      clatoms_info->alp_pol  = (double *)crealloc(&((clatoms_info->alp_pol[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(double))-1;
      clatoms_info->b_neut   = (double *)crealloc(&((clatoms_info->b_neut[1])),
                                          build_intra->natm_tot_max*
                                          sizeof(double))-1;
      clatoms_info->text_atm =(double *)crealloc(&((clatoms_info->text_atm[1])),
                                          build_intra->natm_tot_max*
                                          sizeof(double))-1;
      atommaps->iatm_mol_typ  = (int *) crealloc(&(atommaps->iatm_mol_typ[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->iatm_atm_typ  = (int *) crealloc(&(atommaps->iatm_atm_typ[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->iatm_res_typ  = (int *) crealloc(&(atommaps->iatm_res_typ[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->iatm_mol_num  = (int *) crealloc(&(atommaps->iatm_mol_num[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->iatm_res_num  = (int *) crealloc(&(atommaps->iatm_res_num[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->ighost_flag   = (int *)crealloc(&(atommaps->ighost_flag[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->freeze_flag   = (int *)crealloc(&(atommaps->freeze_flag[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->atom_label   = (int *)crealloc(&(atommaps->atom_label[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
  }/*endif*/

/*==========================================================================*/
/* I) Reallocate the ghosts                                                 */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nghost_add;
  if(ghost_atoms->nghost_tot+nmol1 > build_intra->nghost_tot_max){
           nghost_old  = build_intra->nghost_tot_max;
           ncomp_old   = NCOEF_GHOST_MAX;
           build_intra->nghost_tot_max+=
    MAX((ghost_atoms->nghost_tot+nmol1-build_intra->nghost_tot_max),NMEM_MIN);
           nghost_new  = build_intra->nghost_tot_max;
           ncomp_new   = NCOEF_GHOST_MAX;
           ghost_atoms->ighost_map = (int *) crealloc(
                                     &(ghost_atoms->ighost_map)[1],
                                       nghost_new*sizeof(int))-1;
           ghost_atoms->natm_comp  = (int *) crealloc(  
                                     &(ghost_atoms->natm_comp)[1],
                                       nghost_new*sizeof(int))-1;
           ghost_atoms->iatm_comp  = creall_int_mat(ghost_atoms->iatm_comp,
                                                    1,ncomp_old,1,nghost_old,
                                                    1,ncomp_new,1,nghost_new);
           ghost_atoms->coef       = creall_mat(ghost_atoms->coef,
                                                1,ncomp_old,1,nghost_old,
                                                1,ncomp_new,1,nghost_new);
  }/*endif*/

/*==========================================================================*/
/* I) Reallocate the freeze                                                 */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nfreeze_add;
  if(atommaps->nfreeze+nmol1 > build_intra->nfreeze_max){
           nfreeze_old  = build_intra->nfreeze_max;
           build_intra->nfreeze_max+=
    MAX((atommaps->nfreeze+nmol1-build_intra->nfreeze_max),NMEM_MIN);
           nfreeze_new  = build_intra->nfreeze_max;
           atommaps->freeze_map = (int *) crealloc(
                                     &(atommaps->freeze_map)[1],
                                       nfreeze_new*sizeof(int))-1;
	 }/*endif*/

/*==========================================================================*/
/* I) Reallocate the grpcons                                                */


  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_21_add;

  if(bonded->grp_bond_con.num_21+nmol1 > build_intra->ngrp_21_max){
       build_intra->ngrp_21_max += 
  MAX((bonded->grp_bond_con.num_21+nmol1-build_intra->ngrp_21_max),NMEM_MIN);
       bonded->grp_bond_con.j1_21     = 
              (int *) crealloc(&(bonded->grp_bond_con.j1_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
       bonded->grp_bond_con.j2_21     = 
              (int *) crealloc(&(bonded->grp_bond_con.j2_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
       bonded->grp_bond_con.jtyp_21   = 
              (int *) crealloc(&(bonded->grp_bond_con.jtyp_21)[1],
                          build_intra->ngrp_21_max*sizeof(int))-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_23_add;
  if(bonded->grp_bond_con.num_23+nmol1 > build_intra->ngrp_23_max){
       build_intra->ngrp_23_max += 
  MAX((bonded->grp_bond_con.num_23+nmol1-build_intra->ngrp_23_max),NMEM_MIN);
       bonded->grp_bond_con.j1_23     = 
              (int *) crealloc(&(bonded->grp_bond_con.j1_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       bonded->grp_bond_con.j2_23     = 
              (int *) crealloc(&(bonded->grp_bond_con.j2_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       bonded->grp_bond_con.j3_23     = 
              (int *) crealloc(&(bonded->grp_bond_con.j3_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
       bonded->grp_bond_con.jtyp_23   = 
              (int *) crealloc(&(bonded->grp_bond_con.jtyp_23)[1],
                          build_intra->ngrp_23_max*sizeof(int))-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_33_add;
  if(bonded->grp_bond_con.num_33+nmol1 > build_intra->ngrp_33_max){
       build_intra->ngrp_33_max += 
   MAX((bonded->grp_bond_con.num_33+nmol1-build_intra->ngrp_33_max),NMEM_MIN);
       bonded->grp_bond_con.j1_33     = 
              (int *) crealloc(&(bonded->grp_bond_con.j1_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       bonded->grp_bond_con.j2_33     = 
              (int *) crealloc(&(bonded->grp_bond_con.j2_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       bonded->grp_bond_con.j3_33     = 
              (int *) crealloc(&(bonded->grp_bond_con.j3_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
       bonded->grp_bond_con.jtyp_33   = 
              (int *) crealloc(&(bonded->grp_bond_con.jtyp_33)[1],
                          build_intra->ngrp_33_max*sizeof(int))-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_watt_33_add;
  if(bonded->grp_bond_watts.num_33+nmol1 > build_intra->ngrp_watt_33_max){
       build_intra->ngrp_watt_33_max += 
   MAX((bonded->grp_bond_watts.num_33+nmol1-build_intra->ngrp_watt_33_max),
              NMEM_MIN);
       bonded->grp_bond_watts.j1_33     = 
              (int *) crealloc(&(bonded->grp_bond_watts.j1_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       bonded->grp_bond_watts.j2_33     = 
              (int *) crealloc(&(bonded->grp_bond_watts.j2_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       bonded->grp_bond_watts.j3_33     = 
              (int *) crealloc(&(bonded->grp_bond_watts.j3_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
       bonded->grp_bond_watts.jtyp_33   = 
              (int *) crealloc(&(bonded->grp_bond_watts.jtyp_33)[1],
                          build_intra->ngrp_watt_33_max*sizeof(int))-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_43_add;
  if(bonded->grp_bond_con.num_43+nmol1 > build_intra->ngrp_43_max){
       build_intra->ngrp_43_max += 
  MAX((bonded->grp_bond_con.num_43+nmol1-build_intra->ngrp_43_max),NMEM_MIN);
       bonded->grp_bond_con.j1_43     = 
              (int *) crealloc(&(bonded->grp_bond_con.j1_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       bonded->grp_bond_con.j2_43     = 
              (int *) crealloc(&(bonded->grp_bond_con.j2_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       bonded->grp_bond_con.j3_43     = 
              (int *) crealloc(&(bonded->grp_bond_con.j3_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       bonded->grp_bond_con.j4_43     =  
              (int *) crealloc(&(bonded->grp_bond_con.j4_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
       bonded->grp_bond_con.jtyp_43   = 
              (int *) crealloc(&(bonded->grp_bond_con.jtyp_43)[1],
                          build_intra->ngrp_43_max*sizeof(int))-1;
  }/*endif*/
  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*ngrp_46_add;
  if(bonded->grp_bond_con.num_46+nmol1 > build_intra->ngrp_46_max){
       build_intra->ngrp_46_max += 
  MAX((bonded->grp_bond_con.num_46+nmol1-build_intra->ngrp_46_max),NMEM_MIN);
       bonded->grp_bond_con.j1_46     = 
              (int *) crealloc(&(bonded->grp_bond_con.j1_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       bonded->grp_bond_con.j2_46     = 
              (int *) crealloc(&(bonded->grp_bond_con.j2_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       bonded->grp_bond_con.j3_46     = 
              (int *) crealloc(&(bonded->grp_bond_con.j3_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       bonded->grp_bond_con.j4_46     =  
              (int *) crealloc(&(bonded->grp_bond_con.j4_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
       bonded->grp_bond_con.jtyp_46   = 
              (int *) crealloc(&(bonded->grp_bond_con.jtyp_46)[1],
                          build_intra->ngrp_46_max*sizeof(int))-1;
  }/*endif*/

/*===========================================================================*/
/* I) Reallocate the bonds */

  nmol = atommaps->nmol_jmol_typ[jmol_typ];
  nmol1 = (nmol-1)*nbond_pow_add;
  if(bonded->bond.npow+nmol1 > build_intra->nbond_pow_max){
    build_intra->nbond_pow_max += 
          MAX((bonded->bond.npow+nmol1-build_intra->nbond_pow_max),NMEM_MIN);
    bonded->bond.j1_pow = (int *)crealloc(&(bonded->bond.j1_pow)[1],
	                           build_intra->nbond_pow_max*sizeof(int))-1;
    bonded->bond.j2_pow = (int *)crealloc(&(bonded->bond.j2_pow)[1],
				   build_intra->nbond_pow_max*sizeof(int))-1;
    bonded->bond.jtyp_pow=(int *)crealloc(&(bonded->bond.jtyp_pow)[1],
				   build_intra->nbond_pow_max*sizeof(int))-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbond_con_add;
  if(bonded->bond.ncon+nmol1 > build_intra->nbond_con_max){
    build_intra->nbond_con_max += 
          MAX((bonded->bond.ncon+nmol1-build_intra->nbond_con_max),NMEM_MIN);
    bonded->bond.j1_con = (int *)crealloc(&(bonded->bond.j1_con)[1],
	                           build_intra->nbond_con_max*sizeof(int))-1;
    bonded->bond.j2_con = (int *)crealloc(&(bonded->bond.j2_con)[1],
				   build_intra->nbond_con_max*sizeof(int))-1;
    bonded->bond.jtyp_con=(int *)crealloc(&(bonded->bond.jtyp_con)[1],
				   build_intra->nbond_con_max*sizeof(int))-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbond_null_add;
  if((null_inter_parse->nbond_nul)+nmol1 > build_intra->nbond_nul_max) {
    build_intra->nbond_nul_max   += 
 MAX((null_inter_parse->nbond_nul+nmol1-build_intra->nbond_nul_max),NMEM_MIN);
    null_inter_parse->jbond1_nul  = 
           (int *) crealloc(&(null_inter_parse->jbond1_nul)[1],
		            build_intra->nbond_nul_max*sizeof(int))-1;
    null_inter_parse->jbond2_nul  = 
	(int *) crealloc(&(null_inter_parse->jbond2_nul)[1],
			 build_intra->nbond_nul_max*sizeof(int))-1;
  }  /*endif*/

/*====================================================================*/
/* II) Reallocate the bends                                           */

  nmol1 = (nmol-1)*nbend_pow_add;
  if(bonded->bend.npow+nmol1 > build_intra->nbend_pow_max){
    build_intra->nbend_pow_max += 
        MAX((bonded->bend.npow+nmol1-build_intra->nbend_pow_max),NMEM_MIN);
    bonded->bend.j1_pow =(int*)crealloc(&(bonded->bend.j1_pow)[1], 
				 build_intra->nbend_pow_max*sizeof(int))-1;
    bonded->bend.j2_pow =(int*)crealloc(&(bonded->bend.j2_pow)[1],
				 build_intra->nbend_pow_max*sizeof(int))-1;
    bonded->bend.j3_pow =(int*)crealloc(&(bonded->bend.j3_pow)[1],
				 build_intra->nbend_pow_max*sizeof(int))-1;
    bonded->bend.jtyp_pow=(int*)crealloc(&(bonded->bend.jtyp_pow)[1],
				  build_intra->nbend_pow_max*sizeof(int))-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbend_con_add;
  if(bonded->bend.ncon+nmol1 > build_intra->nbend_con_max){    
    build_intra->nbend_con_max += 
          MAX((bonded->bend.ncon+nmol1-build_intra->nbend_con_max),NMEM_MIN);
    bonded->bend.j1_con = (int *)crealloc(&(bonded->bend.j1_con)[1],
				   build_intra->nbend_con_max*sizeof(int))-1;
    bonded->bend.j2_con = (int *)crealloc(&(bonded->bend.j2_con)[1],
				   build_intra->nbend_con_max*sizeof(int))-1;
    bonded->bend.jtyp_con=(int *)crealloc(&(bonded->bend.jtyp_con)[1],
				   build_intra->nbend_con_max*sizeof(int))-1;
  }/*endif*/

  nmol1 = (nmol-1)*nbend_null_add;
  if(null_inter_parse->nbend_nul+nmol1 > build_intra->nbend_nul_max){
    build_intra->nbend_nul_max += 
 MAX((null_inter_parse->nbend_nul+nmol1-build_intra->nbend_nul_max),NMEM_MIN);
    null_inter_parse->jbend1_nul     = 
	(int *) crealloc(&(null_inter_parse->jbend1_nul)[1],
			 build_intra->nbend_nul_max*sizeof(int))-1;
    null_inter_parse->jbend2_nul     = 
	(int *) crealloc(&(null_inter_parse->jbend2_nul)[1],
			 build_intra->nbend_nul_max*sizeof(int))-1;
    null_inter_parse->jbend3_nul     = 
	(int *) crealloc(&(null_inter_parse->jbend3_nul)[1],
			 build_intra->nbend_nul_max*sizeof(int))-1;
  }     /*endif*/

/*====================================================================*/
/* IV) Reallocate the torsions                                        */

  nmol1 = (nmol-1)*ntors_pow_add;
  if((bonded->tors.npow)+nmol1 > build_intra->ntors_pow_max){
    build_intra->ntors_pow_max += 
        MAX((bonded->tors.npow+nmol1-build_intra->ntors_pow_max),NMEM_MIN);
    bonded->tors.j1_pow =(int*)crealloc(&(bonded->tors.j1_pow)[1],
				 build_intra->ntors_pow_max*sizeof(int))-1;
    bonded->tors.j2_pow =(int*)crealloc(&(bonded->tors.j2_pow)[1],
				 build_intra->ntors_pow_max*sizeof(int))-1;
    bonded->tors.j3_pow =(int*)crealloc(&(bonded->tors.j3_pow)[1],
				 build_intra->ntors_pow_max*sizeof(int))-1;
    bonded->tors.j4_pow =(int*)crealloc(&(bonded->tors.j4_pow)[1],
				 build_intra->ntors_pow_max*sizeof(int))-1;
    bonded->tors.jtyp_pow=(int*)crealloc(&(bonded->tors.jtyp_pow)[1],
				  build_intra->ntors_pow_max*sizeof(int))-1;
  }/*endif*/
    
  nmol1 = (nmol-1)*ntors_con_add;
  if(bonded->tors.ncon+nmol1 > build_intra->ntors_con_max){    
    build_intra->ntors_con_max +=
         MAX((bonded->tors.ncon+nmol1-build_intra->ntors_con_max),NMEM_MIN);
    bonded->tors.j1_con = (int *) crealloc(&(bonded->tors.j1_con)[1],
			            build_intra->ntors_con_max*sizeof(int))-1;
    bonded->tors.j2_con = (int *) crealloc(&(bonded->tors.j2_con)[1],
				    build_intra->ntors_con_max*sizeof(int))-1;
    bonded->tors.j3_con = (int *) crealloc(&(bonded->tors.j3_con)[1],
				    build_intra->ntors_con_max*sizeof(int))-1;
    bonded->tors.j4_con = (int *) crealloc(&(bonded->tors.j4_con)[1],
				    build_intra->ntors_con_max*sizeof(int))-1;
    bonded->tors.jtyp_con=(int *) crealloc(&(bonded->tors.jtyp_con)[1],
				    build_intra->ntors_con_max*sizeof(int))-1;
  } /* endif */

  nmol1 = (nmol-1)*ntors_null_add;
  if(null_inter_parse->ntors_nul+nmol1 > build_intra->ntors_nul_max){
    build_intra->ntors_nul_max += 
  MAX((null_inter_parse->ntors_nul+nmol1-build_intra->ntors_nul_max),NMEM_MIN);
    null_inter_parse->jtors1_nul = 
      (int *) crealloc(&(null_inter_parse->jtors1_nul[1]),
		       build_intra->ntors_nul_max*sizeof(int))-1;
    null_inter_parse->jtors2_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors2_nul[1]),
		       build_intra->ntors_nul_max*sizeof(int))-1;
    null_inter_parse->jtors3_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors3_nul[1]),
		       build_intra->ntors_nul_max*sizeof(int))-1;
    null_inter_parse->jtors4_nul     = 
      (int *) crealloc(&(null_inter_parse->jtors4_nul[1]),
		       build_intra->ntors_nul_max*sizeof(int))-1;
  }      /*endif*/

/*====================================================================*/
/* V) Reallocate the torsions                                         */

  nmol1 = (nmol-1)*nonfo_add;
  if(bonded->onfo.num+nmol1 > build_intra->nonfo_max){
    build_intra->nonfo_max += 
          MAX((bonded->onfo.num+nmol1-build_intra->nonfo_max),NMEM_MIN);
    bonded->onfo.j1 =(int *) crealloc(&(bonded->onfo.j1)[1],
			       build_intra->nonfo_max*sizeof(int))-1;
    bonded->onfo.j2 =(int *) crealloc(&(bonded->onfo.j2)[1],
			       build_intra->nonfo_max*sizeof(int))-1;
    bonded->onfo.jtyp =(int *) crealloc(&(bonded->onfo.jtyp)[1],
				 build_intra->nonfo_max*sizeof(int))-1;
  }/*endif*/

  nmol1 = (nmol-1)*nonfo_null_add;
  if(null_inter_parse->nonfo_nul+nmol1 > build_intra->nonfo_nul_max){
    build_intra->nonfo_nul_max += 
  MAX((null_inter_parse->nonfo_nul+nmol1-build_intra->nonfo_nul_max),NMEM_MIN);
    null_inter_parse->jonfo1_nul = 
      (int *) crealloc(&(null_inter_parse->jonfo1_nul)[1],
		       build_intra->nonfo_nul_max*sizeof(int))-1;
    null_inter_parse->jonfo2_nul     = 
      (int *) crealloc(&(null_inter_parse->jonfo2_nul)[1],
		       build_intra->nonfo_nul_max*sizeof(int))-1;
  }/*endif*/

/*====================================================================*/
/* VI) Reallocate the Uri-Bradleys                                    */

  nmol1 = (nmol-1)*nbend_bnd_add;
  if(bonded->bend_bnd.num+nmol1 > build_intra->nbend_bnd_max){
    mem_add_now = bonded->bend_bnd.num+nmol1-build_intra->nbend_bnd_max;
    build_intra->nbend_bnd_max += MAX(mem_add_now,NMEM_MIN);
    bonded->bend_bnd.j1 =(int*)crealloc(&(bonded->bend_bnd.j1)[1], 
				 build_intra->nbend_bnd_max*sizeof(int))-1;
    bonded->bend_bnd.j2 =(int*)crealloc(&(bonded->bend_bnd.j2)[1],
				 build_intra->nbend_bnd_max*sizeof(int))-1;
    bonded->bend_bnd.j3 =(int*)crealloc(&(bonded->bend_bnd.j3)[1],
				 build_intra->nbend_bnd_max*sizeof(int))-1;
    bonded->bend_bnd.jtyp=(int*)crealloc(&(bonded->bend_bnd.jtyp)[1],
				  build_intra->nbend_bnd_max*sizeof(int))-1;
  }/*endif*/

/*==========================================================================*/
} /*end routine*/
/*==========================================================================*/





