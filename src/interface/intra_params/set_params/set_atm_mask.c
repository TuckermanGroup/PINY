/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                  Mask routines                                           */
/*                                                                          */
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_handle_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_atm_mask(BUILD_INTRA *build_intra,DICT_INTRA *dict_intra,
                    ATOMMAPS *atommaps, char fun_key[], char filename[],
                    int iresidue, int iresidue_off)

/*========================================================================*/
/*     Begin routine                                                      */
   {/*begin routine*/
 int index,i;
/*========================================================================*/
/* Make sure you got the number of atoms in pure residue correct          */

    build_intra->natmind_1res_now = build_intra->natm_1res_now;
    if(build_intra->natm_1res_now!=
              atommaps->natm_jres_jmol_typ[iresidue_off+iresidue]){
      strcpy(fun_key,"res_name_def");
      index = 3;
      keyarg_barf(dict_intra->mol_name_dict,filename,fun_key,index);
    }/*endif*/

/*========================================================================*/
/* Make sure you got the atom indices of pure residue correct             */

    for(i=1;i<=build_intra->natm_1res_now;i++){
      if(build_intra->mask_atm[i]!=1){
        strcpy(fun_key,"res_name_def");
	index = 3;
	keyarg_barf(dict_intra->atm_dict,filename,fun_key,index);
      }/*endif*/
    }/*endfor*/

/*-----------------------------------------------------------------------*/
    }/*end routine*/
/*========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void create_atm_ind(CLATOMS_INFO *clatoms_info,ATOMMAPS *atommaps,
                    BUILD_INTRA *build_intra,
                    DICT_INTRA *dict_intra,char fun_key[],char filename[])

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/
int natm_now,i,index;
/*========================================================================*/
/* Check masks                                                            */

    natm_now = 0;
    for(i=1;i<=build_intra->natmind_1res_now;i++){
      if((build_intra->mask_atm[i]>1)||(build_intra->mask_atm[i]<0)){
        strcpy(fun_key,"atom_create or atom_destroy");
	index = 2;
	keyarg_barf(dict_intra->atm_dict,filename,fun_key,index);
      }/*endif*/
      natm_now += (build_intra->mask_atm)[i] ;
    }/*endfor*/

    if(natm_now!=build_intra->natm_1res_now){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("INTERNAL ERROR: atom masks incorrectly set ");
      printf("in control_res_params.c\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);  
      exit(1);
    }/*endfor*/

/*========================================================================*/
/* Build indicies                                                         */

    build_intra->index_atm[1] = (build_intra->mask_atm)[1];
    for(i=2;i<=build_intra->natmind_1res_now;i++){
      build_intra->index_atm[i] = (build_intra->mask_atm[i] + 
                                   build_intra->index_atm[(i-1)]);
    }/*endfor*/

/*========================================================================*/
/* Reallocate atom memory                                                 */

    if((clatoms_info->natm_tot+build_intra->natm_1res_now) > 
                                          build_intra->natm_tot_max){
      build_intra->natm_tot_max += MAX(NMEM_MIN,natm_now);
      clatoms_info->mass = (double *)crealloc(&((clatoms_info->mass[1])),
                                         build_intra->natm_tot_max*
                                         sizeof(double))-1;
      clatoms_info->q         = (double *)crealloc(&((clatoms_info->q[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(double))-1;
      clatoms_info->cp_vlnc_up = (int *)crealloc(&((clatoms_info->cp_vlnc_up[1])),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      clatoms_info->cp_vlnc_dn = (int *)crealloc(&((clatoms_info->cp_vlnc_dn[1])),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      clatoms_info->cp_atm_flag = (int *)crealloc(&((clatoms_info->cp_atm_flag[1])),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      clatoms_info->alp_pol =(double *)crealloc(&((clatoms_info->alp_pol[1])),
                                              build_intra->natm_tot_max*
                                              sizeof(double))-1;
      clatoms_info->b_neut  = (double *)crealloc(&((clatoms_info->b_neut[1])),
                                          build_intra->natm_tot_max*
                                          sizeof(double))-1;
      clatoms_info->text_atm = 
                            (double *)crealloc(&((clatoms_info->text_atm[1])),
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
      atommaps->ighost_flag  = (int *) crealloc(&(atommaps->ighost_flag[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->freeze_flag  = (int *) crealloc(&(atommaps->freeze_flag[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
      atommaps->atom_label  = (int *) crealloc(&(atommaps->atom_label[1]),
                                                 build_intra->natm_tot_max*
                                                 sizeof(int))-1;
    }/*endif*/
   
/*========================================================================*/
/* Reallocate bond_site memory                                            */

   if(build_intra->natm_1res_now>build_intra->natm_1res_max){
      build_intra->natm_1res_max += 
        MAX(NMEM_MIN,build_intra->natm_1res_now-build_intra->natm_1res_max);
      build_intra->bond_site  = (BOND_SITE *)
                                    crealloc(&(build_intra->bond_site[1]),
                                             build_intra->natm_1res_max*
                                             sizeof(BOND_SITE))-1;
      build_intra->iatm_ind_chk = (int *)
                                    crealloc(&(build_intra->iatm_ind_chk[1]),
                                             build_intra->natm_1res_max*
						sizeof(int))-1;
   }

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_atm_ind(BUILD_INTRA *build_intra)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/
 int i,iii;
/*========================================================================*/

  for(i=1;i<=build_intra->natm_1res_now;i++){
    if(build_intra->iatm_ind_chk[i]!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("INTERNAL ERROR: atom indices off ");
      printf("in control_res_params.c\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);  
      exit(1);
    }/*endif*/
  }/*endfor*/

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void init_build_intra(BUILD_INTRA *build_intra, ATOMMAPS *atommaps,
                       int iresidue, int iresidue_off)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/
int i;
/*========================================================================*/

     build_intra->natm_1res_now      = 0;
     build_intra->natmind_1res_now   = 0;
     build_intra->natm_1res_pure_now = 
                atommaps->natm_jres_jmol_typ[iresidue_off+iresidue];

     if(build_intra->natm_1res_pure_now > 
        build_intra->natmind_1res_max){
      build_intra->natmind_1res_max += 
            MAX(NMEM_MIN,(atommaps->natm_jres_jmol_typ)[iresidue_off+iresidue]
                           - (build_intra->natmind_1res_max));
      build_intra->mask_atm = (int *)
        crealloc(&(build_intra->mask_atm[1]),build_intra->natmind_1res_max*
                 sizeof(int))-1;
      build_intra->index_atm = (int *)
        crealloc(&(build_intra->index_atm[1]),build_intra->natmind_1res_max*
                 sizeof(int))-1;
    }/*endif*/

    for(i=1;i<=build_intra->natmind_1res_max;i++){build_intra->mask_atm[i]=0;}

/*-----------------------------------------------------------------------*/
}/*end routine*/
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_atm_mask(DICT_WORD atm_dict[],int num_atm_dict,
                  char fun_key[],char filename[],
                  BUILD_INTRA *build_intra,int num, int jres,int jmol_typ)

/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/
/*          Local Variables                                             */
  int index;
  int iatm_ind;
/*=======================================================================*/
/* Error check */

  if(num==9 || num==10){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Atom creates and atom destroys not permitted \n");
      printf("in the residue parameter file %s\n",filename);
      printf("of the %d residue of the %d molecule\n",jres,jmol_typ);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
  }/*endif*/

/*=======================================================================*/
/* Set the mask increment atoms counters                                 */

  sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
  index = 2;
  if((iatm_ind<0)||(iatm_ind>build_intra->natm_1res_pure_now)){
    keyarg_barf(atm_dict,filename,fun_key,index);
  }
  build_intra->mask_atm[iatm_ind]++;
  build_intra->natm_1res_now    ++;

/*=======================================================================*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/

void set_atm_mask_rb(DICT_WORD atm_dict[],int num_atm_dict,
		   char *fun_key,char *filename,
                   BUILD_INTRA *build_intra, int num, int jres,int jmol_typ)

/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/
/*              Local Variables                                             */
  int i,index;
  int iatm_ind,itemp;
/*==========================================================================*/

   if(num==2){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Atom definitions not permitted \n");
      printf("in the residue morphing file %s\n",filename);
      printf("of the %d residue of the %d molecule\n",jres,jmol_typ);
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }

/*==========================================================================*/
/* Get atom index                                                           */

    sscanf(atm_dict[2].keyarg,"%d",&iatm_ind);
    index = 2;

/*==========================================================================*/
/* Create new atom or destroy existing atoms                                */

    if(num==9){
       if((iatm_ind<0)||(iatm_ind<=build_intra->natm_1res_pure_now)){
            keyarg_barf(atm_dict,filename,fun_key,index);}
       build_intra->natm_1res_now += 1;
    }/*endif*/

    if(num==8){
       if((iatm_ind<0)||(iatm_ind>build_intra->natm_1res_pure_now)){
            keyarg_barf(atm_dict,filename,fun_key,index);}
      build_intra->natm_1res_now     -= 1;
    /*endif*/}

/*==========================================================================*/
/* Check array sizes                                                        */

    itemp = MAX(build_intra->natmind_1res_now,iatm_ind);
    build_intra->natmind_1res_now=itemp;
    if(build_intra->natmind_1res_now>build_intra->natmind_1res_max){
           build_intra->natmind_1res_max += 
    MAX(NMEM_MIN,build_intra->natmind_1res_now-build_intra->natmind_1res_max);
           build_intra->mask_atm = (int *)
              crealloc(&(build_intra->mask_atm[1]),
		       build_intra->natmind_1res_max*sizeof(int))-1;
           build_intra->index_atm = (int *)
              crealloc(&(build_intra->index_atm[1]),
		       build_intra->natmind_1res_max*sizeof(int))-1;
           build_intra->bond_site  = (BOND_SITE *)
              crealloc(&(build_intra->bond_site[1]),
		       build_intra->natmind_1res_max*sizeof(BOND_SITE))-1;
           build_intra->iatm_ind_chk= (int *)
              crealloc(&(build_intra->bond_site[1]),
		       build_intra->natmind_1res_max*sizeof(int))-1;
           for(i=build_intra->natmind_1res_now+1;
	       i<=build_intra->natmind_1res_max;i++)
	         {build_intra->mask_atm[i]=0;}
     /*endif*/}

/*==========================================================================*/
/* Set atom masks                                                           */

    if(strcasecmp(fun_key,"atom_create")==0){
     build_intra->mask_atm[iatm_ind]++;
    }/*endif*/
    if(strcasecmp(fun_key,"atom_destroy")==0){
      build_intra->mask_atm[iatm_ind]--;
    /*endif*/}

/*==========================================================================*/
/*end routine*/}
/*==========================================================================*/





