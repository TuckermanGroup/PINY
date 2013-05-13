/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: init_intra_parms.c                           */
/*                                                                          */
/* This subprogram starts things up                                         */
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

void init_intra_params(CLATOMS_INFO *clatoms_info,
                       GHOST_ATOMS *ghost_atoms,
                       ATOMMAPS *atommaps, 
                       BUILD_INTRA *build_intra,BONDED *bonded,
                       NULL_INTER_PARSE *null_inter_parse,
		       RESBOND_PARSE *resbond_parse,
		       FILENAME_PARSE *filename_parse)

/*========================================================================*/
/*          Begin Routine                                                 */

   {/*begin routine */

/*========================================================================*/
/*          Local Variables                                               */

   int i,nres_max,nres_tot,iii,num;

/*=======================================================================*/
/* 0) Make jres_jmol_typ_strt map. Find max number of residues.          */
/*    Malloc residue index checker and parameter file memory             */ 

  atommaps->jres_jmol_typ_strt[1] = 1;
  nres_max = MAX(atommaps->nres_1mol_jmol_typ[1],1);
  nres_tot = MAX(atommaps->nres_1mol_jmol_typ[1],1);
  for(i=2;i<=atommaps->nmol_typ;i++){    
    atommaps->jres_jmol_typ_strt[i] = atommaps->jres_jmol_typ_strt[i-1] + 
                      MAX(atommaps->nres_1mol_jmol_typ[i-1],1); 
    nres_tot += MAX(atommaps->nres_1mol_jmol_typ[i],1);
    nres_max = MAX(atommaps->nres_1mol_jmol_typ[i],nres_max);
  }

  build_intra->ires_ind_chk       = (int *) cmalloc(nres_max*sizeof(int))-1;   
  filename_parse->res_param_name = 
               (NAME *) cmalloc(nres_max*sizeof(NAME))-1;   
  filename_parse->res_fix_name = 
               (NAME *) cmalloc(nres_max*sizeof(NAME))-1;   
  atommaps->jatm_jres_1mol_jmol_typ_strt = (int *) 
                                           cmalloc(nres_tot*sizeof(int))-1;
  atommaps->ires_typ_jres_jmol_typ       = (int *) 
                                           cmalloc(nres_tot*sizeof(int))-1;
  atommaps->natm_jres_jmol_typ           = (int *) 
                                           cmalloc(nres_tot*sizeof(int))-1;
  atommaps->nfree_jres_jmol_typ          = (int *) 
                                           cmalloc(nres_tot*sizeof(int))-1;
  atommaps->icons_jres_jmol_typ          = (int *) 
                                           cmalloc(nres_tot*sizeof(int))-1;   
  build_intra->strip1                    = (char *)
                                           cmalloc(MAXWORD*sizeof(char));  
  build_intra->strip2                    = (char *)
                                           cmalloc(MAXWORD*sizeof(char));  

/*========================================================================*/
/* I) Atm stuff */

  clatoms_info->natm_tot        = 0; 
  atommaps->nfreeze             = 0; 
  build_intra->natm_tot_max     = NMEM_MIN;
  build_intra->nfreeze_max      = NMEM_MIN;
  atommaps->natm_typ            = 0; 
  build_intra->natm_typ_max     = NMEM_MIN;
  build_intra->natmind_1res_max = NMEM_MIN;
  build_intra->natm_1res_max    = NMEM_MIN;
  
  atommaps->atm_typ         = (NAME *)  cmalloc(build_intra->natm_typ_max
                                                *sizeof(NAME))-1;
  clatoms_info->mass             = (double *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(double))-1;
  clatoms_info->q                = (double *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(double))-1;
  clatoms_info->cp_vlnc_up      = (int *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  clatoms_info->cp_vlnc_dn      = (int *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  clatoms_info->cp_atm_flag     = (int *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  clatoms_info->alp_pol          = (double *)cmalloc(build_intra->natm_tot_max
                                             *sizeof(double))-1;
  clatoms_info->b_neut           = (double *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(double))-1;
  clatoms_info->text_atm         = (double *)cmalloc(build_intra->natm_tot_max
                                                *sizeof(double))-1;
  atommaps->iatm_mol_typ    = (int *)   cmalloc(build_intra->natm_tot_max
                                               *sizeof(int))-1;
  atommaps->iatm_res_typ    = (int *)   cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  atommaps->iatm_atm_typ    = (int *)   cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  atommaps->iatm_mol_num    = (int *)   cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  atommaps->iatm_res_num    = (int *)   cmalloc(build_intra->natm_tot_max
                                                *sizeof(int))-1;
  atommaps->natm_1mol_jmol_typ=(int *)  cmalloc(atommaps->nmol_typ
                                                *sizeof(int))-1;
  atommaps->jatm_jmol_typ_strt=(int *)  cmalloc(atommaps->nmol_typ
                                                *sizeof(int))-1;
  atommaps->nfree_1mol_jmol_typ=(int *) cmalloc(atommaps->nmol_typ
                                                *sizeof(int))-1;
  atommaps->icons_jmol_typ  = (int *)   cmalloc(atommaps->nmol_typ
                                                *sizeof(int))-1;
  atommaps->ighost_flag     = (int *)   cmalloc(build_intra->natm_tot_max
                                               *sizeof(int))-1;
  atommaps->freeze_flag     = (int *)   cmalloc(build_intra->natm_tot_max
                                               *sizeof(int))-1;
  atommaps->freeze_map      = (int *)   cmalloc(build_intra->nfreeze_max
                                               *sizeof(int))-1;
  atommaps->atom_label      = (int *)   cmalloc(build_intra->natm_tot_max
                                               *sizeof(int))-1;
  build_intra->mask_atm     = (int *)   cmalloc(build_intra->natmind_1res_max
                                                *sizeof(int))-1;
  build_intra->index_atm    = (int *)   cmalloc(build_intra->natmind_1res_max
                                                *sizeof(int))-1;
  build_intra->iatm_ind_chk = (int *)   cmalloc(build_intra->natm_1res_max
                                                *sizeof(int))-1;
  build_intra->bond_site = (BOND_SITE *)
                                        cmalloc(build_intra->natm_1res_max
                                                *sizeof(BOND_SITE))-1;
/*========================================================================*/
/* I) Ghost_atm stuff */

  ghost_atoms->nghost_tot = 0;
  build_intra->nghost_tot_max = NMEM_MIN;
  ghost_atoms->ighost_map     = (int *) cmalloc(build_intra->nghost_tot_max
                                                *sizeof(int))-1;
  ghost_atoms->natm_comp_max  = NCOEF_GHOST_MAX;
  ghost_atoms->natm_comp      = (int *) cmalloc(build_intra->nghost_tot_max
                                                *sizeof(int))-1;
  ghost_atoms->ighost_map     = (int *) cmalloc(build_intra->nghost_tot_max
                                                *sizeof(int))-1;
  ghost_atoms->iatm_comp      = cmall_int_mat(1,ghost_atoms->natm_comp_max,
                                              1,build_intra->nghost_tot_max);
  ghost_atoms->coef           = cmall_mat(1,ghost_atoms->natm_comp_max,
                                          1,build_intra->nghost_tot_max);
/*========================================================================*/
/* II) Grp constrnts */
  bonded->grp_bond_con.num_21  = 0;
  bonded->grp_bond_con.num_33  = 0;
  bonded->grp_bond_con.num_43  = 0;
  bonded->grp_bond_con.num_23  = 0;
  bonded->grp_bond_con.num_46  = 0;
  bonded->grp_bond_watts.num_33  = 0;
  build_intra->ngrp_21_max     = NMEM_MIN;
  build_intra->ngrp_33_max     = NMEM_MIN;
  build_intra->ngrp_43_max     = NMEM_MIN;
  build_intra->ngrp_23_max     = NMEM_MIN;
  build_intra->ngrp_46_max     = NMEM_MIN;
  build_intra->ngrp_watt_33_max = NMEM_MIN;
  bonded->grp_bond_con.ntyp_21  = 0;
  bonded->grp_bond_con.ntyp_33  = 0;
  bonded->grp_bond_con.ntyp_43  = 0;
  bonded->grp_bond_con.ntyp_23  = 0;
  bonded->grp_bond_con.ntyp_46  = 0;
  bonded->grp_bond_watts.ntyp_33  = 0;
  build_intra->ngrp_typ_21_max  = NMEM_MIN;
  build_intra->ngrp_typ_33_max  = NMEM_MIN;
  build_intra->ngrp_typ_43_max  = NMEM_MIN;
  build_intra->ngrp_typ_23_max  = NMEM_MIN;
  build_intra->ngrp_typ_46_max  = NMEM_MIN;
  build_intra->ngrp_typ_watt_33_max  = NMEM_MIN;

  bonded->grp_bond_con.j1_21 = (int *)   cmalloc(build_intra->ngrp_21_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j2_21 = (int *)   cmalloc(build_intra->ngrp_21_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_21 = (int *) cmalloc(build_intra->ngrp_33_max
                                                *sizeof(int))-1;
  bonded->grp_bond_con.j1_33 = (int *)   cmalloc(build_intra->ngrp_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j2_33 = (int *)   cmalloc(build_intra->ngrp_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j3_33 = (int *)   cmalloc(build_intra->ngrp_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_33 = (int *) cmalloc(build_intra->ngrp_33_max
                                                *sizeof(int))-1;
  bonded->grp_bond_watts.j1_33 = (int *) cmalloc(build_intra->ngrp_watt_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_watts.j2_33 = (int *) cmalloc(build_intra->ngrp_watt_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_watts.j3_33 = (int *) cmalloc(build_intra->ngrp_watt_33_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_watts.jtyp_33 = 
                          (int *)cmalloc(build_intra->ngrp_watt_33_max
                                                *sizeof(int))-1;
  bonded->grp_bond_con.j1_23 = (int *)   cmalloc(build_intra->ngrp_23_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j2_23 = (int *)   cmalloc(build_intra->ngrp_23_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j3_23 = (int *)   cmalloc(build_intra->ngrp_23_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_23 = (int *) cmalloc(build_intra->ngrp_23_max
                                                *sizeof(int))-1;
  bonded->grp_bond_con.j1_43 = (int *)   cmalloc(build_intra->ngrp_43_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j2_43 = (int *)   cmalloc(build_intra->ngrp_43_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j3_43 = (int *)   cmalloc(build_intra->ngrp_43_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j4_43 = (int *)   cmalloc(build_intra->ngrp_43_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_43 = (int *) cmalloc(build_intra->ngrp_43_max
                                                *sizeof(int))-1;
  bonded->grp_bond_con.j1_46 = (int *)   cmalloc(build_intra->ngrp_46_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j2_46 = (int *)   cmalloc(build_intra->ngrp_46_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j3_46 = (int *)   cmalloc(build_intra->ngrp_46_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.j4_46 = (int *)   cmalloc(build_intra->ngrp_46_max
                                                *sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_46 = (int *) cmalloc(build_intra->ngrp_46_max
                                                *sizeof(int))-1;
  build_intra->cgrp_typ_21=(CGRP_CONS*) cmalloc(build_intra->ngrp_typ_21_max
                                             *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_33=(CGRP_CONS*) cmalloc(build_intra->ngrp_typ_33_max
                                             *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_43=(CGRP_CONS*) cmalloc(build_intra->ngrp_typ_43_max
                                             *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_23=(CGRP_CONS*) cmalloc(build_intra->ngrp_typ_23_max
                                             *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_46=(CGRP_CONS*) cmalloc(build_intra->ngrp_typ_46_max
                                             *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_watt_33=
          (CGRP_CONS*) cmalloc(build_intra->ngrp_typ_watt_33_max
                                               *sizeof(CGRP_CONS))-1;
  build_intra->cgrp_typ_now=(CGRP_CONS*) cmalloc(sizeof(CGRP_CONS));
/*========================================================================*/
/*   II) Bonds  */
  bonded->bond.npow             =0;
  bonded->bond.ncon             =0;
  null_inter_parse->nbond_nul   =0;
  build_intra->nbond_pow_max    =NMEM_MIN;
  build_intra->nbond_con_max    =NMEM_MIN;
  build_intra->nbond_nul_max    =NMEM_MIN;
  bonded->bond.ntyp_pow         =0;
  bonded->bond.ntyp_con         =0;
  build_intra->nbond_typ_pow_max=NMEM_MIN;
  build_intra->nbond_typ_con_max=NMEM_MIN;
  
  bonded->bond.j1_pow       = (int *)   cmalloc(build_intra->nbond_pow_max 
                                                *sizeof(int))-1;  
  bonded->bond.j2_pow       = (int *)   cmalloc(build_intra->nbond_pow_max
                                                *sizeof(int))-1;  
  bonded->bond.jtyp_pow     = (int *)   cmalloc(build_intra->nbond_pow_max
                                                *sizeof(int))-1;  
  bonded->bond.j1_con       = (int *)   cmalloc(build_intra->nbond_con_max
                                                *sizeof(int))-1;  
  bonded->bond.j2_con       = (int *)   cmalloc(build_intra->nbond_con_max
                                                *sizeof(int))-1;  
  bonded->bond.jtyp_con     = (int *)   cmalloc(build_intra->nbond_con_max
                                                *sizeof(int))-1;  

  num = bonded->rbar_sig_free.nfree;
  bonded->rbar_sig_free.j1 = (int *)cmalloc(num*sizeof(int))-1;
  bonded->rbar_sig_free.j2 = (int *)cmalloc(num*sizeof(int))-1;

  null_inter_parse->jbond1_nul= (int *) cmalloc(build_intra->nbond_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jbond2_nul= (int *) cmalloc(build_intra->nbond_nul_max
                                                *sizeof(int))-1;  
  build_intra->cbond_typ_pow=(CBOND*)cmalloc(build_intra->nbond_typ_pow_max
                                             *sizeof(CBOND))-1;
  build_intra->cbond_typ_con=(CBOND*)cmalloc(build_intra->nbond_typ_con_max
                                             *sizeof(CBOND))-1;
  build_intra->cbond_typ_now=(CBOND*)cmalloc(sizeof(CBOND));

/*========================================================================*/
/*  III) Bends  */

  bonded->bend.npow             =0;
  bonded->bend.ncon             =0;
  null_inter_parse->nbend_nul   =0;
  build_intra->nbend_pow_max    =NMEM_MIN;
  build_intra->nbend_con_max    =NMEM_MIN;
  build_intra->nbend_nul_max    =NMEM_MIN;
  bonded->bend.ntyp_pow         =0;
  bonded->bend.ntyp_con         =0;
  build_intra->nbend_typ_pow_max=NMEM_MIN;
  build_intra->nbend_typ_con_max=NMEM_MIN;
  
  bonded->bend.j1_pow       = (int *)   cmalloc(build_intra->nbend_pow_max
                                                *sizeof(int))-1;  
  bonded->bend.j2_pow       = (int *)   cmalloc(build_intra->nbend_pow_max
                                                *sizeof(int))-1;  
  bonded->bend.j3_pow       = (int *)   cmalloc(build_intra->nbend_pow_max
                                                *sizeof(int))-1;  
  bonded->bend.jtyp_pow     = (int *)   cmalloc(build_intra->nbend_pow_max
                                                *sizeof(int))-1;  
  bonded->bend.j1_con       = (int *)   cmalloc(build_intra->nbend_con_max
                                                *sizeof(int))-1;  
  bonded->bend.j2_con       = (int *)   cmalloc(build_intra->nbend_con_max
                                                *sizeof(int))-1;  
  bonded->bend.j3_con       = (int *)   cmalloc(build_intra->nbend_con_max
                                                *sizeof(int))-1;  
  bonded->bend.jtyp_con     = (int *)   cmalloc(build_intra->nbend_con_max
                                                *sizeof(int))-1;  
  null_inter_parse->jbend1_nul= (int *) cmalloc(build_intra->nbend_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jbend2_nul= (int *) cmalloc(build_intra->nbend_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jbend3_nul= (int *) cmalloc(build_intra->nbend_nul_max
                                                *sizeof(int))-1;  
  build_intra->cbend_typ_pow=(CBEND*)cmalloc(build_intra->nbend_typ_pow_max
                                             *sizeof(CBEND))-1;
  build_intra->cbend_typ_con=(CBEND*)cmalloc(build_intra->nbend_typ_con_max
                                             *sizeof(CBEND))-1;
  build_intra->cbend_typ_now=(CBEND*) cmalloc(sizeof(CBEND));  

/*=======================================================================*/
/*  III.V) Bend_bnd  */

  build_intra->nbend_bnd_max    =NMEM_MIN;
  bonded->bend_bnd.ntyp         =0;
  build_intra->nbend_bnd_typ_max=NMEM_MIN;
  bonded->bend_bnd.num          =0;

  bonded->bend_bnd.j1       = (int *)   cmalloc(build_intra->nbend_bnd_max
                                                *sizeof(int))-1;  
  bonded->bend_bnd.j2       = (int *)   cmalloc(build_intra->nbend_bnd_max
                                                *sizeof(int))-1;  
  bonded->bend_bnd.j3       = (int *)   cmalloc(build_intra->nbend_bnd_max
                                                *sizeof(int))-1;  
  bonded->bend_bnd.jtyp     = (int *)   cmalloc(build_intra->nbend_bnd_max
                                                *sizeof(int))-1;  
  build_intra->cbend_bnd_typ=(CBEND*)cmalloc(build_intra->nbend_bnd_typ_max
                                                *sizeof(CBEND))-1;
  build_intra->cbend_bnd_typ_now=(CBEND*) cmalloc(sizeof(CBEND));  


/*=======================================================================*/
/*  IV) Tors  */

  bonded->tors.npow             =0;
  bonded->tors.ncon             =0;
  bonded->tors.nimpr            =0;
  null_inter_parse->ntors_nul   =0;
  build_intra->ntors_pow_max    =NMEM_MIN;
  build_intra->ntors_con_max    =NMEM_MIN;
  build_intra->ntors_nul_max    =NMEM_MIN;
  bonded->tors.ntyp_pow         =0;
  bonded->tors.ntyp_con         =0;
  build_intra->ntors_typ_pow_max=NMEM_MIN;
  build_intra->ntors_typ_con_max=NMEM_MIN;
  
  bonded->tors.j1_pow       = (int *)   cmalloc(build_intra->ntors_pow_max
                                                *sizeof(int))-1;  
  bonded->tors.j2_pow       = (int *)   cmalloc(build_intra->ntors_pow_max
                                                *sizeof(int))-1;  
  bonded->tors.j3_pow       = (int *)   cmalloc(build_intra->ntors_pow_max
                                                *sizeof(int))-1;  
  bonded->tors.j4_pow       = (int *)   cmalloc(build_intra->ntors_pow_max
                                                *sizeof(int))-1;  
  bonded->tors.jtyp_pow     = (int *)   cmalloc(build_intra->ntors_pow_max
                                                *sizeof(int))-1;  
  bonded->tors.j1_con       = (int *)   cmalloc(build_intra->ntors_con_max
                                                *sizeof(int))-1;  
  bonded->tors.j2_con       = (int *)   cmalloc(build_intra->ntors_con_max
                                                *sizeof(int))-1;  
  bonded->tors.j3_con       = (int *)   cmalloc(build_intra->ntors_con_max
                                                *sizeof(int))-1;  
  bonded->tors.j4_con       = (int *)   cmalloc(build_intra->ntors_con_max
                                                *sizeof(int))-1;  
  bonded->tors.jtyp_con     = (int *)   cmalloc(build_intra->ntors_con_max
                                                *sizeof(int))-1;  
  null_inter_parse->jtors1_nul= (int *) cmalloc(build_intra->ntors_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jtors2_nul= (int *) cmalloc(build_intra->ntors_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jtors3_nul= (int *) cmalloc(build_intra->ntors_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jtors4_nul= (int *) cmalloc(build_intra->ntors_nul_max
                                                *sizeof(int))-1;  
  build_intra->ctors_typ_pow= (CTORS *) cmalloc(build_intra->ntors_typ_pow_max
                                                *sizeof(CTORS))-1;
  build_intra->ctors_typ_con= (CTORS *) cmalloc(build_intra->ntors_typ_con_max
                                                *sizeof(CTORS))-1;
  build_intra->ctors_typ_now= (CTORS *) cmalloc(sizeof(CTORS));

/*========================================================================*/
/*   V) Onfo  */

  bonded->onfo.num           =0;
  null_inter_parse->nonfo_nul=0;
  build_intra->nonfo_max     =NMEM_MIN;
  build_intra->nonfo_nul_max =NMEM_MIN;
  bonded->onfo.ntyp          =0;
  build_intra->nonfo_typ_max =NMEM_MIN;
  
  bonded->onfo.j1           = (int *)   cmalloc(build_intra->nonfo_max
                                                *sizeof(int))-1;  
  bonded->onfo.j2           = (int *)   cmalloc(build_intra->nonfo_max
                                                *sizeof(int))-1;  
  bonded->onfo.jtyp         = (int *)   cmalloc(build_intra->nonfo_max
                                                *sizeof(int))-1;  
  null_inter_parse->jonfo1_nul= (int *) cmalloc(build_intra->nonfo_nul_max
                                                *sizeof(int))-1;  
  null_inter_parse->jonfo2_nul= (int *) cmalloc(build_intra->nonfo_nul_max
                                                *sizeof(int))-1;  
  build_intra->confo_typ    = (CBOND *) cmalloc(build_intra->nonfo_typ_max
                                                *sizeof(CBOND))-1;
  build_intra->confo_typ_now= (CBOND *) cmalloc(sizeof(CBOND));

/*========================================================================*/
/* VI) Residues and residue bonds                                         */

  resbond_parse->nres_bond_max    = NMEM_MIN;
  resbond_parse->nresidue_max     = NMEM_MIN;
  resbond_parse->resbond_prm      = (RESBOND_PRM *)cmalloc(
                           resbond_parse->nres_bond_max*sizeof(RESBOND_PRM))-1;
  resbond_parse->res_bond1_index  = (int *)cmalloc(
                               resbond_parse->nres_bond_max*sizeof(int))-1; 
  resbond_parse->res_bond2_index  = (int *)cmalloc(
                               resbond_parse->nres_bond_max*sizeof(int))-1; 
  resbond_parse->res_bond1_off    = (int *)cmalloc(
                               resbond_parse->nresidue_max*sizeof(int))-1; 
  resbond_parse->res_bond2_off = (int *)cmalloc(
                               resbond_parse->nresidue_max*sizeof(int))-1; 
  resbond_parse->nres_bond1_jres = (int *)cmalloc(
                               resbond_parse->nresidue_max*sizeof(int))-1; 
  resbond_parse->nres_bond2_jres = (int *)cmalloc(
                               resbond_parse->nresidue_max*sizeof(int))-1; 

  build_intra->nres_typ_max     = NMEM_MIN;

  atommaps->res_typ             = (NAME *)cmalloc(
                                build_intra->nres_typ_max*sizeof(NAME))-1;

/*=======================================================================*/
/*   V) Assign mall variables  */

    atommaps->nres_max = nres_max;   

/*-----------------------------------------------------------------------*/
   }/*end routine*/
/*==========================================================================*/






