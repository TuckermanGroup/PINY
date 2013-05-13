/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: close_intra_parms.c                          */
/*                                                                          */
/* This subprogram tidies things up                                         */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../proto_defs/proto_intra_params_entry.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void close_intra_params(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                        GHOST_ATOMS *ghost_atoms,
                        ATOMMAPS *atommaps,
                        BUILD_INTRA *build_intra,
                        BONDED *bonded,
                        NULL_INTER_PARSE *null_inter_parse,
                        double *tot_memory,SIMOPTS *simopts,int np_forc)

/*=======================================================================*/
/*        Begin routine                                                 */
{/*begin routine*/
/*=======================================================================*/
/*           Local Variables                                             */

  int i,ntyp_pow,iii,nsec,inow,know;
  double now_memory;
  int natm_typ_mall;
  int natm_mall;
  int nghost_mall;
  int nfreeze_mall;
  int ncomp_mall;
  int nchrg_mall;
  int ngrp_21_mall;
  int ngrp_typ_21_mall;
  int ngrp_33_mall;
  int ngrp_typ_33_mall;
  int ngrp_watt_33_mall;
  int ngrp_typ_watt_33_mall;
  int ngrp_43_mall;
  int ngrp_typ_43_mall;
  int ngrp_23_mall;
  int ngrp_typ_23_mall;
  int ngrp_46_mall;
  int ngrp_typ_46_mall;
  int nbond_pow_mall;
  int nbond_typ_pow_mall;
  int nbond_con_mall;
  int nbond_typ_con_mall;
  int nbend_pow_mall;
  int nbend_typ_pow_mall;
  int nbend_con_mall;
  int nbend_typ_con_mall;
  int nbend_bnd_mall;
  int nbend_bnd_typ_mall;
  int ntors_pow_mall;
  int ntors_typ_pow_mall;
  int ntors_con_mall;
  int ntors_typ_con_mall;
  int nonfo_mall;
  int nonfo_typ_mall;
  int pimd_on;
  int nghost_old,nfreeze_old,ncomp_old;
  int pi_beads = clatoms_info->pi_beads; 
  int pimd_freez_typ = atommaps->pimd_freez_typ;

/*======================================================================*/
/* 0) Make the length of each vector an odd number                      */

  nghost_old          = build_intra->nghost_tot_max;
  nfreeze_old         = build_intra->nfreeze_max;
  ncomp_old           = NCOEF_GHOST_MAX;
  natm_mall           = clatoms_info->natm_tot;
  natm_typ_mall       = atommaps->natm_typ;
  nghost_mall         = ghost_atoms->nghost_tot;
  nfreeze_mall        = atommaps->nfreeze;
  ncomp_mall          = ghost_atoms->natm_comp_max;
  ngrp_21_mall        = bonded->grp_bond_con.num_21;
  ngrp_33_mall        = bonded->grp_bond_con.num_33;
  ngrp_watt_33_mall   = bonded->grp_bond_watts.num_33;
  ngrp_43_mall        = bonded->grp_bond_con.num_43;
  ngrp_23_mall        = bonded->grp_bond_con.num_23;
  ngrp_46_mall        = bonded->grp_bond_con.num_46;
  ngrp_typ_21_mall    = bonded->grp_bond_con.ntyp_21;
  ngrp_typ_33_mall    = bonded->grp_bond_con.ntyp_33;
  ngrp_typ_watt_33_mall  = bonded->grp_bond_watts.ntyp_33;
  ngrp_typ_43_mall    = bonded->grp_bond_con.ntyp_43;
  ngrp_typ_23_mall    = bonded->grp_bond_con.ntyp_23;
  ngrp_typ_46_mall    = bonded->grp_bond_con.ntyp_46;
  nbond_pow_mall      = bonded->bond.npow;
  nbond_typ_pow_mall  = bonded->bond.ntyp_pow;
  nbond_con_mall      = bonded->bond.ncon;
  nbond_typ_con_mall  = bonded->bond.ntyp_con;
  nbend_pow_mall      = bonded->bend.npow;
  nbend_typ_pow_mall  = bonded->bend.ntyp_pow;
  nbend_con_mall      = bonded->bend.ncon;
  nbend_typ_con_mall  = bonded->bend.ntyp_con;
  nbend_bnd_mall      = bonded->bend_bnd.num;
  nbend_bnd_typ_mall  = bonded->bend_bnd.ntyp;
  ntors_pow_mall      = bonded->tors.npow;
  ntors_typ_pow_mall  = bonded->tors.ntyp_pow;
  ntors_con_mall      = bonded->tors.ncon;
  ntors_typ_con_mall  = bonded->tors.ntyp_con;
  nonfo_mall          = bonded->onfo.num;
  nonfo_typ_mall      = bonded->onfo.ntyp;
  if((natm_mall           % 2)==0){natm_mall      += 1;}
  nsec=natm_mall/np_forc;
  if((natm_mall%np_forc)!=0){nsec++;}
  natm_mall = nsec*np_forc;
  if((natm_typ_mall       % 2)==0){natm_typ_mall  += 1;}
  if(((nghost_mall        % 2)==0)
    &&(nghost_mall        !=0)){nghost_mall++;}
  if(((nfreeze_mall        % 2)==0)
    &&(nfreeze_mall        !=0)){nfreeze_mall++;}
  if(((ncomp_mall      % 2)==0)
    &&(ncomp_mall     !=0)){ncomp_mall++;}
  if(((ngrp_21_mall        % 2)==0)
    &&(ngrp_21_mall       !=0)){ngrp_21_mall++;}
  if(((ngrp_33_mall        % 2)==0)
    &&(ngrp_33_mall       !=0)){ngrp_33_mall++;}
  if(((ngrp_watt_33_mall        % 2)==0)
    &&(ngrp_watt_33_mall       !=0)){ngrp_watt_33_mall++;}
  if(((ngrp_43_mall        % 2)==0)
    &&(ngrp_43_mall       !=0)){ngrp_43_mall++;}
  if(((ngrp_23_mall        % 2)==0)
    &&(ngrp_23_mall       !=0)){ngrp_23_mall++;}
  if(((ngrp_46_mall        % 2)==0)
    &&(ngrp_46_mall       !=0)){ngrp_46_mall++;}
  if(((ngrp_typ_21_mall    % 2)==0)
    &&(ngrp_typ_21_mall   !=0)){ngrp_typ_21_mall++;}
  if(((ngrp_typ_33_mall    % 2)==0)
    &&(ngrp_typ_33_mall   !=0)){ngrp_typ_33_mall++;}
  if(((ngrp_typ_watt_33_mall    % 2)==0)
    &&(ngrp_typ_watt_33_mall   !=0)){ngrp_typ_watt_33_mall++;}
  if(((ngrp_typ_43_mall    % 2)==0)
    &&(ngrp_typ_43_mall   !=0)){ngrp_typ_43_mall++;}
  if(((ngrp_typ_23_mall    % 2)==0)
    &&(ngrp_typ_23_mall   !=0)){ngrp_typ_23_mall++;}
  if(((ngrp_typ_46_mall    % 2)==0)
    &&(ngrp_typ_46_mall   !=0)){ngrp_typ_46_mall++;}
  if(((nbond_pow_mall     % 2)==0)
    &&(nbond_pow_mall     !=0)){nbond_pow_mall    += 1;}
  if(((nbond_typ_pow_mall % 2)==0)
    &&(nbond_typ_pow_mall !=0)){nbond_typ_pow_mall+= 1;}
  if(((nbond_con_mall     % 2)==0)
    &&(nbond_con_mall     !=0)){nbond_pow_mall    += 1;}
  if(((nbond_typ_con_mall % 2)==0)
    &&(nbond_typ_con_mall !=0)){nbond_typ_pow_mall+= 1;}
  if(((nbend_pow_mall     % 2)==0)
    &&(nbend_pow_mall     !=0)){nbend_pow_mall    += 1;}
  if(((nbend_typ_pow_mall % 2)==0)
    &&(nbend_typ_pow_mall !=0)){nbend_typ_pow_mall+= 1;}
  if(((nbend_con_mall     % 2)==0)
    &&(nbend_con_mall     !=0)){nbend_pow_mall    += 1;}
  if(((nbend_typ_con_mall % 2)==0)
    &&(nbend_typ_con_mall !=0)){nbend_typ_pow_mall+= 1;}
  if(((nbend_bnd_mall     % 2)==0)
    &&(nbend_bnd_mall     !=0)){nbend_bnd_mall    += 1;}
  if(((nbend_bnd_typ_mall % 2)==0)
    &&(nbend_bnd_typ_mall !=0)){nbend_bnd_typ_mall+= 1;}
  if(((ntors_pow_mall     % 2)==0)
    &&(ntors_pow_mall     !=0)){ntors_pow_mall    += 1;}
  if(((ntors_typ_pow_mall % 2)==0)
    &&(ntors_typ_pow_mall !=0)){ntors_typ_pow_mall+= 1;}
  if(((ntors_con_mall     % 2)==0)
    &&(ntors_con_mall     !=0)){ntors_pow_mall    += 1;}
  if(((ntors_typ_con_mall % 2)==0)
    &&(ntors_typ_con_mall !=0)){ntors_typ_pow_mall+= 1;}
  if(((nonfo_mall         % 2)==0)
    &&(nonfo_mall         !=0)){nonfo_mall        += 1;}
  if(((nonfo_typ_mall     % 2)==0)
    &&(nonfo_typ_mall     !=0)){nonfo_typ_mall    += 1;}


/*======================================================================*/
/* I) Atm stuff: Note constrain flag set here                           */
/*               and charged atoms are counted                          */

  bonded->constrnt.iconstrnt=0;
  if(bonded->bond.ncon > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_21 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_23 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_33 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_watts.num_33 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_43 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->grp_bond_con.num_46 > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->bend.ncon > 0) bonded->constrnt.iconstrnt=1;
  if(bonded->tors.ncon > 0) bonded->constrnt.iconstrnt=1;


  clatoms_info->ichrg      = (int *)cmalloc(natm_mall*sizeof(int))-1;
  clatoms_info->nfree      = 0;
  clatoms_info->nfree_pimd = 0;
  for(i=1;i<=atommaps->nmol_typ;i++){
    inow                 = atommaps->nfree_1mol_jmol_typ[i]
                          *atommaps->nmol_jmol_typ[i];
    know                 = 3*atommaps->natm_1mol_jmol_typ[i]
                            *atommaps->nmol_jmol_typ[i];
    clatoms_info->nfree += inow;
    if(pimd_freez_typ==2){
      clatoms_info->nfree_pimd += inow*pi_beads;
    }else{
      clatoms_info->nfree_pimd += (inow + know*(pi_beads-1));
    }/*endif*/
  }/*endfor*/

  clatoms_info->nchrg = 0;
  for(i=1;i<=clatoms_info->natm_tot;i++){
    clatoms_info->mass[i] *= PROT_MASS;
    if(clatoms_info->q[i] != 0) {
     clatoms_info->nchrg++;
     clatoms_info->ichrg[(clatoms_info->nchrg)] = i;
   }
  }
  nchrg_mall = clatoms_info->nchrg;
  if((nchrg_mall % 2)==0){nchrg_mall += 1;}


  now_memory = ((double)(((natm_mall)*sizeof(double)*13+sizeof(int)*4)+
			 natm_typ_mall*(sizeof(double)*0+sizeof(int)*25))
		*1.e-06);
  *tot_memory += now_memory;

  printf("Atom allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,*tot_memory);

  atommaps->atm_typ = (NAME *)crealloc(&(atommaps->atm_typ)[1],
				       natm_typ_mall*sizeof(NAME))-1;
  clatoms_info->mass     = (double *)crealloc(&(clatoms_info->mass)[1],
					 natm_mall*sizeof(double))-1;
  clatoms_info->q        = (double *)crealloc(&(clatoms_info->q)[1],
					 natm_mall*sizeof(double))-1;
  clatoms_info->cp_vlnc_up  = (int *)crealloc(&(clatoms_info->cp_vlnc_up)[1],
					 natm_mall*sizeof(int))-1;
  clatoms_info->cp_vlnc_dn  = (int *)crealloc(&(clatoms_info->cp_vlnc_dn)[1],
					 natm_mall*sizeof(int))-1;
  clatoms_info->cp_atm_flag  = (int *)crealloc(&(clatoms_info->cp_atm_flag)[1],
					 natm_mall*sizeof(int))-1;
  clatoms_info->ichrg  = (int *)crealloc(&(clatoms_info->ichrg)[1],
				    nchrg_mall*sizeof(int))-1;
  clatoms_info->alp_pol  = (double *)crealloc(&(clatoms_info->alp_pol)[1],
					 natm_mall*sizeof(double))-1;
  clatoms_info->b_neut   = (double *)crealloc(&(clatoms_info->b_neut)[1],
					 natm_mall*sizeof(double))-1;
  clatoms_info->text_atm = (double *)crealloc(&(clatoms_info->text_atm)[1],
					 natm_mall*sizeof(double))-1;
  atommaps->iatm_mol_typ  = (int *)crealloc(&(atommaps->iatm_mol_typ)[1],
					    natm_mall*sizeof(int))-1;
  atommaps->iatm_atm_typ  = (int *)crealloc(&(atommaps->iatm_atm_typ)[1],
					    natm_mall*sizeof(int))-1;
  atommaps->iatm_res_typ  = (int *)crealloc(&(atommaps->iatm_res_typ)[1],
					    natm_mall*sizeof(int))-1;
  atommaps->iatm_mol_num  = (int *)crealloc(&(atommaps->iatm_mol_num)[1],
					    natm_mall*sizeof(int))-1;
  atommaps->iatm_res_num  = (int *)crealloc(&(atommaps->iatm_res_num)[1],
					    natm_mall*sizeof(int))-1;
  pimd_on = simopts->pimd + simopts->cp_pimd 
          + simopts->cp_wave_pimd + simopts->cp_wave_min_pimd 
          + simopts->debug_pimd + simopts->debug_cp_pimd;

  clatoms_info->xold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  clatoms_info->yold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  clatoms_info->zold= (double *)cmalloc(natm_mall*sizeof(double))-1;
  if(pi_beads>1||pimd_on==1){
   clatoms_info->xmod = (double *)cmalloc(natm_mall*sizeof(double))-1;
   clatoms_info->ymod = (double *)cmalloc(natm_mall*sizeof(double))-1;
   clatoms_info->zmod = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }/*endif*/
  if(clatoms_info->pi_beads>1||pimd_on==1){
   clatoms_info->prekf = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }
  clatoms_info->roll_sc   = (double *)cmalloc(natm_mall*sizeof(double))-1;
/*========================================================================*/
/* I) Ghost_atm stuff */

  ghost_atoms->ighost_map     = (int *) crealloc(&(ghost_atoms->ighost_map)[1],
                                                 nghost_mall*sizeof(int))-1;
  ghost_atoms->natm_comp      = (int *) crealloc(&(ghost_atoms->natm_comp)[1],
                                                 nghost_mall*sizeof(int))-1;
  atommaps->ighost_flag       = (int *)   crealloc(&(atommaps->ighost_flag)[1],
						 natm_mall*sizeof(int))-1;

  /* need to set up reallocate matrix */
  ghost_atoms->iatm_comp      = creall_int_mat(ghost_atoms->iatm_comp,
                                               1,ncomp_old ,1,nghost_old,
                                               1,ncomp_mall,1,nghost_mall);
  ghost_atoms->coef           = creall_mat(ghost_atoms->coef,
                                           1,ncomp_old ,1,nghost_old,
                                           1,ncomp_mall,1,nghost_mall);
/*========================================================================*/
/* I) Freeze_atm stuff */

  atommaps->freeze_map     = (int *) crealloc(&(atommaps->freeze_map)[1],
                                                 nfreeze_mall*sizeof(int))-1;
  atommaps->freeze_flag       = (int *)   crealloc(&(atommaps->freeze_flag)[1],
						 natm_mall*sizeof(int))-1;
  atommaps->atom_label       = (int *)   crealloc(&(atommaps->atom_label)[1],
						 natm_mall*sizeof(int))-1;
/*========================================================================*/
/* II) Grp constrnts */



  bonded->grp_bond_con.j1_21 = 
                    (int *)crealloc(&(bonded->grp_bond_con.j1_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_21 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_21 = 
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_21)[1],
                                             ngrp_21_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_33 = 
                    (int *)crealloc(&(bonded->grp_bond_con.j1_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_33 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_33 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_33 = 
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_33)[1],
                                             ngrp_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j1_33 = 
                    (int *)crealloc(&(bonded->grp_bond_watts.j1_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.j2_33 = 
                    (int *)   crealloc(&(bonded->grp_bond_watts.j2_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;  
  bonded->grp_bond_watts.j3_33 = 
                    (int *)   crealloc(&(bonded->grp_bond_watts.j3_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;  
  bonded->grp_bond_watts.jtyp_33 = 
                    (int *) crealloc(&(bonded->grp_bond_watts.jtyp_33)[1],
                                        ngrp_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_23 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j1_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j2_23 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_23 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_23 = 
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_23)[1],
                                             ngrp_23_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_43 = 
                    (int *)crealloc(&(bonded->grp_bond_con.j1_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j2_43 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_43 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j4_43 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j4_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_43 = 
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_43)[1],
                                             ngrp_43_mall*sizeof(int))-1;
  bonded->grp_bond_con.j1_46 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j1_46)[1],
				             ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j2_46 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j2_46)[1],
                                              ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j3_46 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j3_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.j4_46 = 
                    (int *)   crealloc(&(bonded->grp_bond_con.j4_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;  
  bonded->grp_bond_con.jtyp_46 = 
                    (int *) crealloc(&(bonded->grp_bond_con.jtyp_46)[1],
                                             ngrp_46_mall*sizeof(int))-1;
  bonded->grp_bond_con.al_21   = cmall_mat(1,1,1,ngrp_21_mall);
  bonded->grp_bond_con.al_33   = cmall_mat(1,3,1,ngrp_33_mall);
  bonded->grp_bond_con.al_43   = cmall_mat(1,3,1,ngrp_43_mall);
  bonded->grp_bond_con.al_23   = cmall_mat(1,2,1,ngrp_23_mall);
  bonded->grp_bond_con.al_46   = cmall_mat(1,6,1,ngrp_46_mall);
  bonded->grp_bond_con.eq_21   = cmall_mat(1,1,1,ngrp_typ_21_mall);
  bonded->grp_bond_con.eq_33   = cmall_mat(1,3,1,ngrp_typ_33_mall);
  bonded->grp_bond_watts.eq_33 = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_con.eq_43   = cmall_mat(1,3,1,ngrp_typ_43_mall);
  bonded->grp_bond_con.eq_23   = cmall_mat(1,2,1,ngrp_typ_23_mall);
  bonded->grp_bond_con.eq_46   = cmall_mat(1,6,1,ngrp_typ_46_mall);

  bonded->grp_bond_watts.cos_thet0_2   = 
                (double *)malloc(ngrp_typ_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.sin_thet0_2   = 
                (double *)malloc(ngrp_typ_watt_33_mall*sizeof(int))-1;
  bonded->grp_bond_watts.c_0_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_1_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_2_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_3_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_4_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_5_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.c_6_33    = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_0_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_1_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_2_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_3_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_4_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_5_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);
  bonded->grp_bond_watts.dc_6_33   = cmall_mat(1,3,1,ngrp_typ_watt_33_mall);

/*======================================================================*/
/*   II) Bonds  */

  now_memory = ((nbond_pow_mall)*(sizeof(double)*0+sizeof(int)*3)+
		(nbond_con_mall)*(sizeof(double)*1 +  sizeof(int)*3)+
		(nbond_typ_pow_mall)*(sizeof(double)*14+ sizeof(int)*0 )+
		(nbond_typ_con_mall)*(sizeof(double)*1 +  sizeof(int)*0)+
		(nbend_pow_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
		(nbend_typ_pow_mall)*(sizeof(double)*29+  sizeof(int)*0)+
		(nbend_con_mall)*(sizeof(double)*1 +  sizeof(int)*4)+
		(nbend_typ_con_mall)*(sizeof(double)*1 +  sizeof(int)*0)+
		(nbend_bnd_mall)*(sizeof(double)*0 +  sizeof(int)*4)+
		(nbend_bnd_typ_mall)*(sizeof(double)*44 +  sizeof(int)*0)+
		(ntors_pow_mall)*(sizeof(double)*0 +  sizeof(int)*5)+
		(ntors_con_mall)*(sizeof(double)*1 +  sizeof(int)*5)+
		(ntors_typ_pow_mall)*(sizeof(double)*30 +  sizeof(int)*0)+
		(nonfo_mall)*(sizeof(double)*0 +  sizeof(int)*3)+
		(nonfo_typ_mall)*(sizeof(double)*3 +  sizeof(int)*0)
		)*1.e-06;

  *tot_memory += now_memory;
  printf("Intramolecular allocation: %g Mbytes; Total memory: %g Mbytes\n",
          now_memory,*tot_memory);
  ntyp_pow = nbond_typ_pow_mall;


  bonded->bond.j1_pow = (int *)    crealloc(&(bonded->bond.j1_pow)[1],
					    nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j2_pow = (int *)    crealloc(&(bonded->bond.j2_pow)[1],
					    nbond_pow_mall*sizeof(int))-1;
  bonded->bond.jtyp_pow   = (int *)crealloc(&(bonded->bond.jtyp_pow)[1],
					    nbond_pow_mall*sizeof(int))-1;
  bonded->bond.j1_con = (int *)    crealloc(&(bonded->bond.j1_con)[1],
					    nbond_con_mall*sizeof(int))-1;
  bonded->bond.j2_con = (int *)    crealloc(&(bonded->bond.j2_con)[1],
					    nbond_con_mall*sizeof(int))-1;
  bonded->bond.jtyp_con   = (int *)crealloc(&(bonded->bond.jtyp_con)[1],
					    nbond_con_mall*sizeof(int))-1;
  bonded->bond.al_con = (double *)cmalloc(nbond_con_mall*sizeof(double))-1;
  bonded->bond.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.eq_pow_res = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bond.eq_con = (double *)cmalloc(nbond_typ_con_mall
					  *sizeof(double))-1;

  null_inter_parse->jbond1_nul= (int *)
    crealloc(&(null_inter_parse->jbond1_nul)[1],
	     null_inter_parse->nbond_nul*sizeof(int))-1;
  null_inter_parse->jbond2_nul= (int *)
    crealloc(&(null_inter_parse->jbond2_nul)[1],
	     null_inter_parse->nbond_nul*sizeof(int))-1;


/*=======================================================================*/
/*  III) Bends  */

  ntyp_pow = nbend_typ_pow_mall;
  bonded->bend.j1_pow = (int *)   crealloc(&(bonded->bend.j1_pow)[1],
					   nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j2_pow = (int *)   crealloc(&(bonded->bend.j2_pow)[1],
					   nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j3_pow = (int *)   crealloc(&(bonded->bend.j3_pow)[1],
					   nbend_pow_mall*sizeof(int))-1;
  bonded->bend.jtyp_pow = (int *)crealloc(&(bonded->bend.jtyp_pow)[1],
					  nbend_pow_mall*sizeof(int))-1;
  bonded->bend.j1_con = (int *)   crealloc(&(bonded->bend.j1_con)[1],
					   nbend_con_mall*sizeof(int))-1;
  bonded->bend.j2_con = (int *)   crealloc(&(bonded->bend.j2_con)[1],
					   nbend_con_mall*sizeof(int))-1;
  bonded->bend.j3_con = (int *)   crealloc(&(bonded->bend.j3_con)[1],
					   nbend_con_mall*sizeof(int))-1;
  bonded->bend.jtyp_con   = (int *)crealloc(&(bonded->bend.jtyp_con)[1],
					    nbend_con_mall*sizeof(int))-1;
  bonded->bend.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->bend.eq_con = (double *)cmalloc(nbend_typ_con_mall
					  *sizeof(double))-1;
  bonded->bend.al_con = (double *)cmalloc(nbend_con_mall
					  *sizeof(double))-1;
  null_inter_parse->jbend1_nul= (int *)  
    crealloc(&(null_inter_parse->jbend1_nul)[1],
	     null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend2_nul= (int *)  
    crealloc(&(null_inter_parse->jbend2_nul)[1],
	     null_inter_parse->nbend_nul*sizeof(int))-1;
  null_inter_parse->jbend3_nul= (int *)  
    crealloc(&(null_inter_parse->jbend3_nul)[1],
	     null_inter_parse->nbend_nul*sizeof(int))-1;
/*=======================================================================*/
/*  III) Bend_bnds  */

  bonded->bend_bnd.j1 = (int *)   crealloc(&(bonded->bend_bnd.j1)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j2 = (int *)   crealloc(&(bonded->bend_bnd.j2)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.j3 = (int *)   crealloc(&(bonded->bend_bnd.j3)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.jtyp = (int *)crealloc(&(bonded->bend_bnd.jtyp)[1],
                                   nbend_bnd_mall*sizeof(int))-1;
  bonded->bend_bnd.eq_bond = (double *)cmalloc(nbend_bnd_typ_mall*
					       sizeof(double))-1;
  bonded->bend_bnd.eq_bend = (double *)cmalloc(nbend_bnd_typ_mall*
					       sizeof(double))-1;
  bonded->bend_bnd.cbond_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbond_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbond_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;

  bonded->bend_bnd.cbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.cbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_0    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_1    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_2    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_3    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_4    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_5    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.sbend_6    = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dcbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_0   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_1   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_2   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_3   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_4   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_5   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;
  bonded->bend_bnd.dsbend_6   = (double *)cmalloc(nbend_bnd_typ_mall
                                                  *sizeof(double))-1;

/*======================================================================*/
/*  IV) Tors  */

  ntyp_pow = ntors_typ_pow_mall;
  bonded->tors.j1_pow = (int *)   crealloc(&(bonded->tors.j1_pow)[1],
					   ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j2_pow = (int *)   crealloc(&(bonded->tors.j2_pow)[1],
					   ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j3_pow = (int *)   crealloc(&(bonded->tors.j3_pow)[1],
					   ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j4_pow = (int *)   crealloc(&(bonded->tors.j4_pow)[1],
					   ntors_pow_mall*sizeof(int))-1;
  bonded->tors.jtyp_pow = (int *)crealloc(&(bonded->tors.jtyp_pow)[1],
					    ntors_pow_mall*sizeof(int))-1;
  bonded->tors.j1_con = (int *)   crealloc(&(bonded->tors.j1_con)[1],
					   ntors_con_mall*sizeof(int))-1;
  bonded->tors.j2_con = (int *)   crealloc(&(bonded->tors.j2_con)[1],
					   ntors_con_mall*sizeof(int))-1;
  bonded->tors.j3_con = (int *)   crealloc(&(bonded->tors.j3_con)[1],
					   ntors_con_mall*sizeof(int))-1;
  bonded->tors.j4_con = (int *)   crealloc(&(bonded->tors.j4_con)[1],
					   ntors_con_mall*sizeof(int))-1;
  bonded->tors.jtyp_con = (int *)crealloc(&(bonded->tors.jtyp_con)[1],
					  ntors_con_mall*sizeof(int))-1;
  bonded->tors.al_con = (double *)cmalloc(ntors_con_mall
					  *sizeof(double))-1;
  bonded->tors.eq_pow = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.c_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_0    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_1    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_2    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_3    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_4    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_5    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.s_6    = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.dc_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_0   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_1   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_2   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_3   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_4   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_5   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.ds_6   = (double *)cmalloc(ntyp_pow*sizeof(double))-1;
  bonded->tors.eq_con = (double *)cmalloc(ntors_typ_pow_mall
                                            *sizeof(double))-1;
  null_inter_parse->jtors1_nul= (int *) 
    crealloc(&(null_inter_parse->jtors1_nul)[1],
	     null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors2_nul= (int *) 
    crealloc(&(null_inter_parse->jtors2_nul)[1],
	     null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors3_nul= (int *) 
    crealloc(&(null_inter_parse->jtors3_nul)[1],
	     null_inter_parse->ntors_nul*sizeof(int))-1;
  null_inter_parse->jtors4_nul= (int *) 
    crealloc(&(null_inter_parse->jtors4_nul)[1],
	     null_inter_parse->ntors_nul*sizeof(int))-1;

/*=======================================================================*/
/*   V) Onfo  */

  bonded->onfo.j1     = (int *)   crealloc(&(bonded->onfo.j1)[1],
					   nonfo_mall*sizeof(int))-1;
  bonded->onfo.j2     = (int *)   crealloc(&(bonded->onfo.j2)[1],
					   nonfo_mall*sizeof(int))-1;
  bonded->onfo.jtyp   = (int *)   crealloc(&(bonded->onfo.jtyp)[1],
					   nonfo_mall*sizeof(int))-1;
  bonded->onfo.feps   = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.s6 = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  bonded->onfo.sc = (double *)cmalloc(nonfo_typ_mall*sizeof(double))-1;
  null_inter_parse->jonfo1_nul= (int *) 
    crealloc(&(null_inter_parse->jonfo1_nul)[1],
	     null_inter_parse->nonfo_nul*sizeof(int))-1;
  null_inter_parse->jonfo2_nul  = (int *) 
    crealloc(&(null_inter_parse->jonfo2_nul)[1],
	     null_inter_parse->nonfo_nul*sizeof(int))-1;

/*=======================================================================*/
/*   V) Assign mall variables  */

    clatoms_info->natm_mall = natm_mall;   
    clatoms_info->nchrg_mall = nchrg_mall;   

    ghost_atoms->nghost_mall = nghost_mall;   
    ghost_atoms->ncomp_mall = ncomp_mall;   

    atommaps->nfreeze_mall = nfreeze_mall;   
    atommaps->natm_typ_mall = natm_typ_mall;   

    bonded->bond.nbond_pow_mall = nbond_pow_mall;    
    bonded->bond.nbond_typ_pow_mall = nbond_typ_pow_mall;
    bonded->bond.nbond_con_mall = nbond_con_mall;    
    bonded->bond.nbond_typ_con_mall = nbond_typ_con_mall;

    bonded->grp_bond_con.ngrp_21_mall = ngrp_21_mall;
    bonded->grp_bond_watts.ngrp_33_mall = ngrp_watt_33_mall;
    bonded->grp_bond_con.ngrp_33_mall = ngrp_33_mall;
    bonded->grp_bond_con.ngrp_43_mall = ngrp_43_mall;
    bonded->grp_bond_con.ngrp_23_mall = ngrp_23_mall;
    bonded->grp_bond_con.ngrp_46_mall = ngrp_46_mall;    
    bonded->grp_bond_con.ngrp_typ_21_mall = ngrp_typ_21_mall;
    bonded->grp_bond_con.ngrp_typ_33_mall = ngrp_typ_33_mall;
    bonded->grp_bond_watts.ngrp_typ_33_mall = ngrp_typ_watt_33_mall;
    bonded->grp_bond_con.ngrp_typ_43_mall = ngrp_typ_43_mall;
    bonded->grp_bond_con.ngrp_typ_23_mall = ngrp_typ_23_mall;
    bonded->grp_bond_con.ngrp_typ_46_mall = ngrp_typ_46_mall;

    bonded->bend.nbend_con_mall     = nbend_con_mall;
    bonded->bend.nbend_typ_con_mall = nbend_typ_con_mall;
    bonded->bend.nbend_pow_mall     = nbend_pow_mall;
    bonded->bend.nbend_typ_pow_mall = nbend_typ_pow_mall;

    bonded->bend_bnd.nbend_bnd_mall = nbend_bnd_mall;
    bonded->bend_bnd.nbend_bnd_typ_mall = nbend_bnd_typ_mall;

    bonded->tors.ntors_pow_mall = ntors_pow_mall;
    bonded->tors.ntors_typ_pow_mall = ntors_typ_pow_mall;
    bonded->tors.ntors_con_mall = ntors_con_mall;
    bonded->tors.ntors_typ_con_mall = ntors_typ_con_mall;

    bonded->onfo.nonfo_mall = nonfo_mall;
    bonded->onfo.nonfo_typ_mall = nonfo_typ_mall;

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/












