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
#include "../proto_defs/proto_friend_lib_entry.h"

#define DEBUG_OFF



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_res_params(double *tot_memory,CLATOMS_INFO *clatoms_info,
                        GHOST_ATOMS *ghost_atoms,
                        ATOMMAPS *atommaps,BONDED *bonded,
                        RESBOND_PARSE *resbond_parse,
                        BUILD_INTRA *build_intra,
                        FILENAME_PARSE *filename_parse,
                        FREE_PARSE *free_parse,CLASS_PARSE *class_parse,
                        NULL_INTER_PARSE *null_inter_parse,
                        char filename[],DICT_INTRA *dict_intra,
                        char fun_key[],int jmol_typ)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int iresidue,nresidue,iparm,iresidue_off;
  int i,istart,iend,jatm,iii;
  int mol_only,natm_mol,mol_or_res,nres_bond;

/*=======================================================================*/
/* 0) Initialize offsets and flags                                       */

  iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  nresidue     = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
  mol_or_res   = 2;
  mol_only     = 0;
  if(atommaps->nres_1mol_jmol_typ[jmol_typ]==0){mol_only=1;}
  atommaps->jatm_jmol_typ_strt[jmol_typ] = clatoms_info->natm_tot+1;

/*=======================================================================*/
/* I) Setup  residues of molecule jmol_typ                          */
  
  for(iresidue=1;iresidue<=nresidue;iresidue++){
    atommaps->jatm_jres_1mol_jmol_typ_strt[(iresidue+iresidue_off)]=
                                                 clatoms_info->natm_tot;
    build_intra->nghost_now = 0;
/*--------------------------------------------------------------------------*/
/*  0) Output                                                               */

#ifdef DEBUG
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    printf("**************************************************************\n");
    printf("Reading from residue parameter file %s\n",filename);
    printf("--------------------------------------------------------------\n");
    printf("\n");
#endif

/*--------------------------------------------------------------------------*/
/*  A) Get residue name and number of atoms if not a single molecule        */
    if(mol_only != 1){
/*    i) Read from residue parm file                                 */
        iparm = 1;
        strcpy(filename, filename_parse->res_param_name[iresidue]);
        check_parmfile(filename,&(dict_intra->num_fun_dict),
                       &(dict_intra->fun_dict),fun_key,
                       nresidue,&nres_bond,mol_or_res);
        fetch_residue_name(filename,dict_intra,atommaps,fun_key,build_intra,
                           jmol_typ,iresidue,iresidue_off,iparm,mol_only);
/*   ii) Read from residue bond file                                 */
        iparm = 0;
        istart = resbond_parse->res_bond1_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
           strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
           if(strlen(filename)!=0){
             check_parmfile(filename,&(dict_intra->num_fun_dict),
                            &(dict_intra->fun_dict),fun_key,
                            nresidue,&nres_bond,mol_or_res);
             fetch_residue_name(filename,dict_intra,atommaps,fun_key,
                    build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
           }/*endif*/
        }/*endfor*/
        iparm = 0;
        istart = resbond_parse->res_bond2_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
           strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
           if(strlen(filename)!=0){
             check_parmfile(filename,&(dict_intra->num_fun_dict),
                            &(dict_intra->fun_dict),fun_key,
                            nresidue,&nres_bond,mol_or_res);
             fetch_residue_name(filename,dict_intra,atommaps,fun_key,
                  build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
           }/*endif*/
        }/*endfor*/
/*   iii) Read from residue fix file                                 */
        iparm = 0;
        strcpy(filename,filename_parse->res_fix_name[iresidue]);
       if(strlen(filename)!=0){
        check_parmfile(filename,&(dict_intra->num_fun_dict),
                      &(dict_intra->fun_dict),fun_key,
                      nresidue,&nres_bond,mol_or_res);
        fetch_residue_name(filename,dict_intra,atommaps,fun_key,
                 build_intra,jmol_typ,iresidue,iresidue_off,iparm,mol_only);
       }/*endif*/
    }else{
        atommaps->natm_jres_jmol_typ[(iresidue_off+1)] = 
                            atommaps->natm_1mol_jmol_typ[jmol_typ];
    }/*endif*/
/*-------------------------------------------------------------------------*/
/*  B) Get atm indicies: set_params/set_atm_mask.c and fetch_residue.c     */
/*    i) Initialize build intra                                            */
        init_build_intra(build_intra,atommaps,iresidue,iresidue_off);
/*   ii) Read from residue parm file :                                     */
        iparm = 1;
        strcpy(filename, filename_parse->res_param_name[iresidue]);
        fetch_residue_atm_masks(filename,dict_intra,atommaps,
                            fun_key,build_intra,jmol_typ,iresidue, 
                            iresidue_off,iparm);
/*   iii) Check atom masks                                                 */
        check_atm_mask(build_intra,dict_intra,atommaps,fun_key,filename,
                       iresidue,iresidue_off);
/*   iv) Read from residue bond files                                      */
        iparm = 0;
        istart = resbond_parse->res_bond1_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
           strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
           if(strlen(filename)!=0){
             fetch_residue_atm_masks(filename,dict_intra,atommaps,
                               fun_key,build_intra,jmol_typ,iresidue, 
                               iresidue_off,iparm);
           }/*endif*/
        }/*endfor*/
        iparm = 0;
        istart = resbond_parse->res_bond2_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
           strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
           if(strlen(filename)!=0){
            fetch_residue_atm_masks(filename,dict_intra,atommaps,
                              fun_key,build_intra,jmol_typ,iresidue, 
                              iresidue_off,iparm);
           }/*endif*/
       }/*endfor*/
/*  v) Read from the fix file                                           */
        iparm = 0;
        strcpy(filename, filename_parse->res_fix_name[iresidue]);
       if(strlen(filename)!=0){
        fetch_residue_atm_masks(filename,dict_intra,atommaps,
                            fun_key,build_intra,jmol_typ,iresidue, 
                            iresidue_off,iparm);
       }/*endif*/
/*  vi) Make the atom indices                                               */
        create_atm_ind(clatoms_info,atommaps,build_intra,
                       dict_intra,fun_key,filename);
/*------------------------------------------------------------------------*/
/*  C) Get atms parameters                                                */

/*    i) Initialize the atom index checker                               */
        for(i=1;i<=build_intra->natm_1res_now;i++){
          (build_intra->iatm_ind_chk)[i]=0;
        }/*endfor*/
/*   ii) Read from residue parm file                                      */
        iparm = 1;
        strcpy(filename, filename_parse->res_param_name[iresidue]);
        fetch_residue_atm_prms(filename,dict_intra,atommaps,clatoms_info,
                   ghost_atoms,
                 fun_key,build_intra,jmol_typ,iresidue,iresidue_off,iparm);
/*  iii) Read residue bond files                                          */
        iparm = 0;
        istart = resbond_parse->res_bond1_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
          strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
          if(strlen(filename)!=0){
             fetch_residue_atm_prms(filename,dict_intra,atommaps,clatoms_info,
                    ghost_atoms,
                    fun_key,build_intra,jmol_typ,iresidue,iresidue_off,iparm);
          }/*endif*/
        }/*endfor*/
        iparm = 0;
        istart = resbond_parse->res_bond2_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
          strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
          if(strlen(filename)!=0){
             fetch_residue_atm_prms(filename,dict_intra,atommaps,clatoms_info,
                    ghost_atoms,
                    fun_key,build_intra,jmol_typ,iresidue,iresidue_off,iparm);
           }/*endif*/
        }/*endfor*/
/*  iv) Read the residue fix file                                    */
        iparm = 0;
        strcpy(filename, filename_parse->res_fix_name[iresidue]);
       if(strlen(filename)!=0){
        fetch_residue_atm_prms(filename,dict_intra,atommaps,
                   clatoms_info,ghost_atoms,
                 fun_key,build_intra,jmol_typ,iresidue,iresidue_off,iparm);
       }/*endif*/
/*    v) Check atm indicies and set atom numbers */
        check_atm_ind(build_intra);
        atommaps->natm_jres_jmol_typ[(iresidue_off+iresidue)]=
                           build_intra->natm_1res_now;
/*------------------------------------------------------------------------*/
/*  D) Get bond/bend/tors/onfo/grp_con params                             */

/*  i) Initial degrees of freedom and constraint options                  */
        atommaps->nfree_jres_jmol_typ[(iresidue_off+iresidue)] = 
                   3.0*(atommaps->natm_jres_jmol_typ[(iresidue_off+iresidue)]
                      - build_intra->nghost_now);
        atommaps->icons_jres_jmol_typ[(iresidue_off+iresidue)] = 0;
/*  ii) Read from residue parm file                                      */
        iparm = 1;
        strcpy(filename, filename_parse->res_param_name[iresidue]);
        fetch_residue_connectivity(filename,dict_intra,
                   atommaps,fun_key,clatoms_info,build_intra,bonded,
                   null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
                   class_parse->mol_hydrog_con_opt[jmol_typ]);
/*  iii) Read residue bond files                                          */
        iparm = 0;
        istart = resbond_parse->res_bond1_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond1_jres[iresidue];i++){
          strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res1_file);
          if(strlen(filename)!=0){
            fetch_residue_connectivity(filename,dict_intra,
                     atommaps,fun_key,clatoms_info,build_intra,bonded,
                     null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
                     class_parse->mol_hydrog_con_opt[jmol_typ]);
         }/*endif*/
        }/*endfor*/
        iparm = 0;
        istart = resbond_parse->res_bond2_off[iresidue];
        for(i=1;i<=resbond_parse->nres_bond2_jres[iresidue];i++){
          strcpy(filename,resbond_parse->resbond_prm[(istart+i)].res2_file);
          if(strlen(filename)!=0){
            fetch_residue_connectivity(filename,dict_intra,
                     atommaps,fun_key,clatoms_info, build_intra,bonded,
                     null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
                     class_parse->mol_hydrog_con_opt[jmol_typ]);
          }/*endif*/
        }/*endfor*/
/*  iv) Read residue fix files                                          */
        iparm = 0;
          strcpy(filename, filename_parse->res_fix_name[iresidue]);
          if(strlen(filename)!=0){
          fetch_residue_connectivity(filename,dict_intra,
                    atommaps,fun_key,clatoms_info,
                    build_intra,bonded,
                    null_inter_parse,jmol_typ,iresidue,iresidue_off,iparm,
                    class_parse->mol_hydrog_con_opt[jmol_typ]);
          }/*endif*/
/*------------------------------------------------------------------------*/
/*  E) Use connectivity information stored in the bondsite structure      */
/*     to fill the res_bondrpm structure for each residue bond to which   */
/*     the current residue is a member                                    */
        if(resbond_parse->nres_bond1_jres[iresidue]>0 ||
           resbond_parse->nres_bond2_jres[iresidue]>0){
          fetch_resbond_prm(resbond_parse,build_intra->bond_site,
                            jmol_typ,iresidue,build_intra->natm_1res_now,
                            clatoms_info->natm_tot);
        }/*endif*/
/*------------------------------------------------------------------------*/
/*  F) Set free energy indicies                                           */
     natm_mol = atommaps->natm_1mol_jmol_typ[jmol_typ];
     fetch_free_energy_index(build_intra,free_parse,bonded,
                           jmol_typ,iresidue,clatoms_info->natm_tot,natm_mol);
/*------------------------------------------------------------------------*/
/*  G) Increment total number of atoms in system                        */
        clatoms_info->natm_tot += build_intra->natm_1res_now;
/*------------------------------------------------------------------------*/
/*  H) Output                                                             */
#ifdef DEBUG
    printf("\n");
    printf("--------------------------------------------------------------\n");
    strcpy(filename, filename_parse->res_param_name[iresidue]);
    printf("Completed read from residue parameter file %s\n",filename);
    printf("**************************************************************\n");
    printf("\n");     
#endif

  }/*endfor:iresidue*/

/*=======================================================================*/
/* II) Get the number of degrees of freedom of each molecule             */

  iresidue_off = atommaps->jres_jmol_typ_strt[jmol_typ]-1;
  atommaps->nfree_1mol_jmol_typ[jmol_typ] = 0; 
  atommaps->natm_1mol_jmol_typ[jmol_typ]  = 0;
  atommaps->icons_jmol_typ[jmol_typ] = 0;
  nresidue = MAX(atommaps->nres_1mol_jmol_typ[jmol_typ],1);
  for(i=1;i<=nresidue;i++){
    atommaps->nfree_1mol_jmol_typ[jmol_typ] += 
                           atommaps->nfree_jres_jmol_typ[(iresidue_off+i)];
    atommaps->natm_1mol_jmol_typ[jmol_typ] += 
                           atommaps->natm_jres_jmol_typ[(iresidue_off+i)];
    if(atommaps->icons_jres_jmol_typ[(iresidue_off+i)]==1){
       atommaps->icons_jmol_typ[jmol_typ] = 1;
    }
#ifdef DEBUG
printf("jmol,ires,mol_free,res_free,natm_1mol,natm_res,icons_mol,icons_res\n");
printf(" %d   %d     %d        %d       %d          %d        %d      %d\n",
               jmol_typ,i,
               atommaps->nfree_1mol_jmol_typ[jmol_typ],
               atommaps->nfree_jres_jmol_typ[(iresidue_off+i)],
               atommaps->natm_1mol_jmol_typ[jmol_typ],
               atommaps->natm_jres_jmol_typ[(iresidue_off+i)],
               atommaps->icons_jmol_typ[jmol_typ],
               atommaps->icons_jres_jmol_typ[(iresidue_off+i)]);
mal_verify(1);
scanf("%d",&iii);
#endif
  }/*endfor*/
  if(natm_mol!=atommaps->natm_1mol_jmol_typ[jmol_typ]){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Error in number of atoms in molecule number %d\n",jmol_typ);
     printf("%d expected %d found\n",natm_mol,
                                     atommaps->natm_1mol_jmol_typ[jmol_typ]);
     printf("in file %s \n",filename_parse->mol_param_name[jmol_typ]);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
    }/*endif*/

/*=======================================================================*/
/* III) Build temperature stuff                                          */

    istart = atommaps->jatm_jmol_typ_strt[jmol_typ];
    iend = atommaps->natm_1mol_jmol_typ[jmol_typ]+istart-1;
    for(jatm=istart;jatm<=iend;jatm++){
      clatoms_info->text_atm[jatm] =  class_parse->text_nhc_mol[jmol_typ];
    }/*endfor*/
  
/*=======================================================================*/
/* IV) Debug                                                             */
#ifdef DEBUG
   vomit_intra_list(clatoms_info,ghost_atoms,atommaps,bonded,null_inter_parse);
#endif
/*-----------------------------------------------------------------------*/
} /*end routine */
/*=======================================================================*/


#ifdef DEBUG


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vomit_intra_list (CLATOMS_INFO *clatoms_info, GHOST_ATOMS *ghost_atoms, 
                       ATOMMAPS *atommaps, BONDED *bonded,
                       NULL_INTER_PARSE *null_inter_parse)

/*==========================================================================*/
{ /*begin routine */
/*==========================================================================*/

int i,iii,j;

/*==========================================================================*/
/* i) Atoms */


printf("------\n");mal_verify(1);
printf("Atoms: %d\n",clatoms_info->natm_tot);
printf("------\n");
   for(i=1;i<=clatoms_info->natm_tot;i++){
      printf("%d mass    %g\n",i,clatoms_info->mass[i]);
      printf("%d q       %g\n",i,clatoms_info->q[i]);            
      printf("%d alp_pol %g\n",i,clatoms_info->alp_pol[i]);     
      printf("%d b_neut  %g\n",i,clatoms_info->b_neut[i]);       
      printf("%d text    %g\n",i,clatoms_info->text_atm[i]);
      printf("%d mol num     %d\n",i,atommaps->iatm_mol_num[i]);
      printf("%d res num     %d\n",i,atommaps->iatm_res_num[i]);
      printf("%d ind mol_typ %d\n",i,atommaps->iatm_mol_typ[i]);
      printf("%d ind res_typ %d\n",i,atommaps->iatm_res_typ[i]);
      printf("%d ind atm_typ %d\n",i,atommaps->iatm_atm_typ[i]);
      printf("%d ind ghost_flag %d\n",i,atommaps->ighost_flag[i]);
      printf("%d ind mol_typ %s\n",i,
            atommaps->mol_typ[atommaps->iatm_mol_typ[i]]);
      printf("%d ind res_typ %s\n",i,
            atommaps->res_typ[atommaps->iatm_res_typ[i]]);
      printf("%d ind atm_typ %s\n",i,
            atommaps->atm_typ[atommaps->iatm_atm_typ[i]]);
      if((i%10)==0){scanf("%d",&iii);}
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Bonds */

printf("------\n");
printf("pow bonds: %d\n ",bonded->bond.npow);
printf("------\n");
   for(i=1;i<=bonded->bond.npow;i++){
     printf("%d atom1 %d atom2 %d type %d \n",i,
             bonded->bond.j1_pow[i],bonded->bond.j2_pow[i],
             bonded->bond.jtyp_pow[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con bonds: %d\n ",bonded->bond.ncon);
printf("------\n");
   for(i=1;i<=bonded->bond.ncon;i++){
     printf("%d atom1 %d atom2 %d type %d \n",i,
             bonded->bond.j1_con[i],bonded->bond.j2_con[i],
             bonded->bond.jtyp_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("null bonds: %d\n ",null_inter_parse->nbond_nul);
printf("------\n");
   for(i=1;i<=null_inter_parse->nbond_nul;i++){
     printf("%d atom1 %d atom2 %d \n",i,
             null_inter_parse->jbond1_nul[i],
             null_inter_parse->jbond2_nul[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Bends */

printf("------\n");
printf("pow bends: %d\n",bonded->bend.npow);
printf("------\n");
   for(i=1;i<=bonded->bend.npow;i++){
     printf("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
             bonded->bend.j1_pow[i],bonded->bend.j2_pow[i],
             bonded->bend.j3_pow[i],bonded->bend.jtyp_pow[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con bends: %d\n",bonded->bend.ncon);
printf("------\n");
   for(i=1;i<=bonded->bend.ncon;i++){
     printf("%d atom1 %d atom2 %d atom 3 %d type %d \n",i,
             bonded->bend.j1_con[i],bonded->bend.j2_con[i],
             bonded->bend.j3_con[i],bonded->bend.jtyp_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("null bends: %d\n ",null_inter_parse->nbend_nul);
printf("------\n");
   for(i=1;i<=null_inter_parse->nbend_nul;i++){
     printf("%d atom1 %d atom2 %d atom3 %d \n",i,
             null_inter_parse->jbend1_nul[i],
             null_inter_parse->jbend2_nul[i],
             null_inter_parse->jbend3_nul[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Bend_bnds */

printf("------\n");
printf("bend_bnd: %d\n",bonded->bend_bnd.num);
printf("------\n");
   for(i=1;i<=bonded->bend_bnd.num;i++){
     printf("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
             bonded->bend_bnd.j1[i],bonded->bend_bnd.j2[i],
             bonded->bend_bnd.j3[i],bonded->bend_bnd.jtyp[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Tors */
printf("------\n");
printf("pow torsions: %d\n",bonded->tors.npow);
printf("------\n");
   for(i=1;i<=bonded->tors.npow;i++){
     printf("%d atom1 %d atom2 %d atom3 %d atom 4 %d type %d \n",i,
             bonded->tors.j1_pow[i],bonded->tors.j2_pow[i],
             bonded->tors.j3_pow[i],bonded->tors.j4_pow[i],
             bonded->tors.jtyp_pow[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("con torsions: %d\n",bonded->tors.ncon);
printf("------\n");
   for(i=1;i<=bonded->tors.ncon;i++){
     printf("%d atom1 %d atom2 %d atom3 %d atom 4 %d type %d \n",i,
             bonded->tors.j1_con[i],bonded->tors.j2_con[i],
             bonded->tors.j3_con[i],bonded->tors.j4_con[i],
             bonded->tors.jtyp_con[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("null torsions: %d\n ",null_inter_parse->ntors_nul);
printf("------\n");
   for(i=1;i<=null_inter_parse->ntors_nul;i++){
     printf("%d atom1 %d atom2 %d atom3 %d atom4 %d\n",i,
             null_inter_parse->jtors1_nul[i],
             null_inter_parse->jtors2_nul[i],
             null_inter_parse->jtors3_nul[i],
             null_inter_parse->jtors4_nul[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Onfos */

printf("------\n");
printf("onfos: %d\n",bonded->onfo.num);
printf("------\n");
   for(i=1;i<=bonded->onfo.num;i++){
     printf("%d atom1 %d atom2 %d type %d \n",i,
             bonded->onfo.j1[i],bonded->onfo.j2[i],
             bonded->onfo.jtyp[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

printf("------\n");
printf("null onfos: %d\n ",null_inter_parse->nbond_nul);
printf("------\n");
   for(i=1;i<=null_inter_parse->nbond_nul;i++){
     printf("%d atom1 %d atom2 %d  \n",i,
             null_inter_parse->jbond1_nul[i],
             null_inter_parse->jbond2_nul[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Ghosts */

printf("------\n");
printf("ghost atoms: %d\n",ghost_atoms->nghost_tot);
printf("------\n");
  
   for(i=1;i<=ghost_atoms->nghost_tot;i++){
     printf("atm %d is ghost number %d and is formed from %d atms\n",
             ghost_atoms->ighost_map[i],i,ghost_atoms->natm_comp[i]);
     for(j=1;j<=ghost_atoms->natm_comp[i];j++){
       printf("atm %d atm_ind %d coef %g\n",
               j,ghost_atoms->iatm_comp[j][i],ghost_atoms->coef[j][i]);
     }/*endfor*/
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 21 */

printf("------\n");mal_verify(1);
printf("Grp con 21: %d \n",bonded->grp_bond_con.num_21);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_con.num_21;i++){
     printf("%d atom1 %d atom2 %d type %d \n",i,
             bonded->grp_bond_con.j1_21[i],bonded->grp_bond_con.j2_21[i],
             bonded->grp_bond_con.jtyp_21[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 23 */

printf("------\n");mal_verify(1);
printf("Grp con 23: %d \n",bonded->grp_bond_con.num_23);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_con.num_23;i++){
     printf("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
             bonded->grp_bond_con.j1_23[i],bonded->grp_bond_con.j2_23[i],
             bonded->grp_bond_con.j3_23[i],
             bonded->grp_bond_con.jtyp_23[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 33 */

printf("------\n");mal_verify(1);
printf("Grp con 33: %d \n",bonded->grp_bond_con.num_33);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_con.num_33;i++){
     printf("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
             bonded->grp_bond_con.j1_33[i],bonded->grp_bond_con.j2_33[i],
             bonded->grp_bond_con.j3_33[i],
             bonded->grp_bond_con.jtyp_33[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp Watts 33 */

printf("------\n");mal_verify(1);
printf("Grp Watts 33: %d \n",bonded->grp_bond_watts.num_33);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_watts.num_33;i++){
     printf("%d atom1 %d atom2 %d atom3 %d type %d \n",i,
             bonded->grp_bond_watts.j1_33[i],bonded->grp_bond_watts.j2_33[i],
             bonded->grp_bond_watts.j3_33[i],
             bonded->grp_bond_watts.jtyp_33[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 43 */

printf("------\n");mal_verify(1);
printf("Grp con 43: %d \n",bonded->grp_bond_con.num_43);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_con.num_43;i++){
     printf("%d atom1 %d atom2 %d atom3 %d atom4 %d type %d \n",i,
             bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j2_43[i],
             bonded->grp_bond_con.j3_43[i],
             bonded->grp_bond_con.j4_43[i],
             bonded->grp_bond_con.jtyp_43[i]);
   }/*endfor*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* i) Grp con 46 */

printf("------\n");mal_verify(1);
printf("Grp con 46: %d \n",bonded->grp_bond_con.num_46);
printf("------\n");

   for(i=1;i<=bonded->grp_bond_con.num_46;i++){
     printf("%d atom1 %d atom2 %d atom3 %d atom4 %d type %d \n",i,
             bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j2_46[i],
             bonded->grp_bond_con.j3_46[i],bonded->grp_bond_con.j4_46[i],
             bonded->grp_bond_con.jtyp_46[i]);
   }/*endfor*/

mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*==========================================================================*/
/* ii) Free energy */

printf("------\n");mal_verify(1);
printf("Bond free energy: %d \n",bonded->bond_free.num);
printf("------\n");
   if(bonded->bond_free.num>0){
    printf("file name %s atom1 %d atom2 %d\n",
            bonded->bond_free.file,bonded->bond_free.j1,bonded->bond_free.j2);
    printf("fk %g eq %g npow %d \n",
            bonded->bond_free.fk,bonded->bond_free.eq,bonded->bond_free.npow);
    printf("nhist %d del %g rmin %g rmax %g\n", 
            bonded->bond_free.nhist,bonded->bond_free.del,
            bonded->bond_free.rmin,bonded->bond_free.rmax);
   }/*endif*/

printf("------\n");mal_verify(1);
printf("Bend free energy: %d \n",bonded->bend_free.num);
printf("------\n");
   if(bonded->bend_free.num>0){
    printf("file name %s atom1 %d atom2 %d atom3 %d\n",
            bonded->bend_free.file,bonded->bend_free.j1,bonded->bend_free.j2,
            bonded->bend_free.j3);
    printf("fk %g eq %g npow %d \n",
            bonded->bend_free.fk,bonded->bend_free.eq,bonded->bend_free.npow);
    printf("nhist %d del %g \n", 
            bonded->bend_free.nhist,bonded->bend_free.del);
   }/*endif*/

printf("------\n");mal_verify(1);
printf("Tors free energy: %d \n",bonded->tors_free.num);
printf("------\n");
   if(bonded->tors_free.num>0){
    printf("file name %s atom1 %d atom2 %d atom3 %d atom4 %d\n",
            bonded->tors_free.file,bonded->tors_free.j1,bonded->tors_free.j2,
            bonded->tors_free.j3,bonded->tors_free.j4);
    printf("fk %g eq %g npow %d \n",
            bonded->tors_free.fk,bonded->tors_free.eq,bonded->tors_free.npow);
    printf("nhist %d del %g \n", 
            bonded->tors_free.nhist,bonded->tors_free.del);
   }/*endif*/
mal_verify(1);printf("Enter an integer ");scanf("%d",&iii);

/*=======================================================================*/
}/*end routine */
/*=======================================================================*/

void mal_verify(int i){}

#endif

