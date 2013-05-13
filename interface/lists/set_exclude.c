/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: set_exclude                                  */
/*                                                                          */
/* This subprogram sets exclusion and ewald corrections                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void set_exclude(CLATOMS_INFO *clatoms_info,GHOST_ATOMS *ghost_atoms,
                 BONDED *bonded,EXCL *excl,NULL_INTER_PARSE *null_inter_parse,
		 int iperd,double *tot_memory, double alp_ewd,
                 int error_check_on)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int num_excl_now;
  int *jtmp1,*jtmp2;
  int i,ioff,nghost_exl,k,iii;
  double now_memory;
  int kstart,kend,j,itemp;

  int bond_npow  = bonded->bond.npow;
  int bend_npow  = bonded->bend.npow;
  int nbond_minus = 0;
  int nbend_minus = 0;


/*========================================================================*/
/* 0) Output */

if(error_check_on==1){
  PRINT_LINE_STAR;
  printf("Finding the excluded interactions\n");
  PRINT_LINE_DASH;printf("\n");
}/*endif*/

/*========================================================================*/
/*  I) Allocate temporary memory */
  
  nghost_exl = 0;
  for(i=1;i<=ghost_atoms->nghost_tot;i++){
    nghost_exl +=   ghost_atoms->natm_comp[i];
  }/*endfor*/
  num_excl_now  =  (bond_npow + bonded->bond.ncon
		    + null_inter_parse->nbond_nul
		    + bend_npow+bonded->bend.ncon
		    + null_inter_parse->nbend_nul
		    + bonded->tors.npow+bonded->tors.ncon
		    + null_inter_parse->ntors_nul
		    + bonded->onfo.num+null_inter_parse->nonfo_nul
		    + bonded->bend_bnd.num
                    + 1*bonded->grp_bond_con.num_21
                    + 2*bonded->grp_bond_con.num_23
                    + 3*bonded->grp_bond_con.num_33
                    + 3*bonded->grp_bond_watts.num_33
                    + 3*bonded->grp_bond_con.num_43
                    + 6*bonded->grp_bond_con.num_46
                    + nghost_exl);

  jtmp1      = (int *)cmalloc(num_excl_now*sizeof(int))-1;
  jtmp2      = (int *)cmalloc(num_excl_now*sizeof(int))-1;

/*========================================================================*/
/* II) Set up exclusion list                                              */
/*------------------------------------------------------------------------*/
/*    1)Exclude Bonds */

  ioff = 0;

  for(i=1;i<= bond_npow;i++){  
    jtmp1[(i+ioff)] =  MAX(bonded->bond.j1_pow[i],bonded->bond.j2_pow[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->bond.j1_pow[i],bonded->bond.j2_pow[i]);
  }/*endfor*/

  ioff += bond_npow;  

  for(i=1;i<=bonded->bond.ncon;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->bond.j1_con[i],bonded->bond.j2_con[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->bond.j1_con[i],bonded->bond.j2_con[i]);
  }/*endfor*/
  ioff += bonded->bond.ncon;
  for(i=1;i<=null_inter_parse->nbond_nul;i++){
    jtmp1[(i+ioff)] =  MAX(null_inter_parse->jbond1_nul[i],
			 null_inter_parse->jbond2_nul[i]);
    jtmp2[(i+ioff)] =  MIN(null_inter_parse->jbond1_nul[i],
			 null_inter_parse->jbond2_nul[i]);
  }/*endfor*/
  ioff +=null_inter_parse->nbond_nul;

/*----------------------------------------------------------------------*/
/*    2)Exclude Bends */

  for(i=1;i<= bend_npow;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->bend.j1_pow[i],bonded->bend.j3_pow[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->bend.j1_pow[i],bonded->bend.j3_pow[i]);
  }/*endfor*/
  ioff += bend_npow;
  for(i=1;i<=bonded->bend.ncon;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->bend.j1_con[i],bonded->bend.j3_con[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->bend.j1_con[i],bonded->bend.j3_con[i]);
  }/*endfor*/
  ioff += bonded->bend.ncon;
  for(i=1;i<=null_inter_parse->nbend_nul;i++){
    jtmp1[(i+ioff)] =  MAX(null_inter_parse->jbend1_nul[i],
			 null_inter_parse->jbend3_nul[i]);
    jtmp2[(i+ioff)] =  MIN(null_inter_parse->jbend1_nul[i],
			 null_inter_parse->jbend3_nul[i]);
  }/*endfor*/
  ioff += null_inter_parse->nbend_nul;

 /*----------------------------------------------------------------------*/
 /*    3)Exclude Torsions */

  for(i=1;i<=bonded->tors.npow;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->tors.j1_pow[i],bonded->tors.j4_pow[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->tors.j1_pow[i],bonded->tors.j4_pow[i]);
  }/*endfor*/
  ioff += bonded->tors.npow;
  for(i=1;i<=bonded->tors.ncon;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->tors.j1_con[i],bonded->tors.j4_con[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->tors.j1_con[i],bonded->tors.j4_con[i]);
  }/*endfor*/
  ioff += bonded->tors.ncon;
  for(i=1;i<=null_inter_parse->ntors_nul;i++){
    jtmp1[(i+ioff)] =  MAX(null_inter_parse->jtors1_nul[i],
			 null_inter_parse->jtors4_nul[i]);
    jtmp2[(i+ioff)] =  MIN(null_inter_parse->jtors1_nul[i],
			 null_inter_parse->jtors4_nul[i]);
  }/*endfor*/
  ioff += null_inter_parse->ntors_nul;

/*---------------------------------------------------------------------*/
/*    4)Exclude Onefours */

  for(i=1;i<=bonded->onfo.num;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->onfo.j1[i],bonded->onfo.j2[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->onfo.j1[i],bonded->onfo.j2[i]);
  }/*endfor*/
  ioff += bonded->onfo.num;
  for(i=1;i<=null_inter_parse->nonfo_nul;i++){
    jtmp1[(i+ioff)] =  MAX(null_inter_parse->jonfo1_nul[i],
			 null_inter_parse->jonfo2_nul[i]);
    jtmp2[(i+ioff)] =  MIN(null_inter_parse->jonfo1_nul[i],
			 null_inter_parse->jonfo2_nul[i]);
  }/*endfor*/
  ioff += null_inter_parse->nonfo_nul;

/*---------------------------------------------------------------------*/
/*    5)Exclude bend_bonds */
  for(i=1;i<=bonded->bend_bnd.num;i++){
    jtmp1[(i+ioff)] =  MAX(bonded->bend_bnd.j1[i],bonded->bend_bnd.j3[i]);
    jtmp2[(i+ioff)] =  MIN(bonded->bend_bnd.j1[i],bonded->bend_bnd.j3[i]);
  }/*endfor*/
  ioff += bonded->bend_bnd.num;

/*---------------------------------------------------------------------*/
/*    6)Exclude ghost-ghost composition atoms */
  for(i=1;i<=ghost_atoms->nghost_tot;i++){
    for(k=1;k<=ghost_atoms->natm_comp[i];k++){
      jtmp1[(k+ioff)] =  
          MAX(ghost_atoms->ighost_map[i],ghost_atoms->iatm_comp[k][i]);
      jtmp2[(k+ioff)] =  
          MIN(ghost_atoms->ighost_map[i],ghost_atoms->iatm_comp[k][i]);
    }/*endfor*/
    ioff += ghost_atoms->natm_comp[i];
  }/*endfor*/

/*---------------------------------------------------------------------*/
/*  7)Exclude group constraints */

/*   i) 21 [1-2]                */
#ifdef DEBUG
  printf("number of group constrains %8d \n",bonded->grp_bond_con.num_21);
#endif
  for(i=1;i<=bonded->grp_bond_con.num_21;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_21[i],bonded->grp_bond_con.j2_21[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_21[i],bonded->grp_bond_con.j2_21[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_21;

/*   i) 23 [1-2,1-3]                */
  for(i=1;i<=bonded->grp_bond_con.num_23;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_23[i],bonded->grp_bond_con.j2_23[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_23[i],bonded->grp_bond_con.j2_23[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_23;
  for(i=1;i<=bonded->grp_bond_con.num_23;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_23[i],bonded->grp_bond_con.j3_23[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_23[i],bonded->grp_bond_con.j3_23[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_23;

/*   i) 33 [1-2,1-3,2-3]            */
#ifdef DEBUG
  printf("number of group constrains %8d \n",bonded->grp_bond_con.num_33);
#endif
  for(i=1;i<=bonded->grp_bond_con.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_33[i],bonded->grp_bond_con.j2_33[i]);
    jtmp2[(i+ioff)] =  

         MIN(bonded->grp_bond_con.j1_33[i],bonded->grp_bond_con.j2_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_33;
  for(i=1;i<=bonded->grp_bond_con.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_33[i],bonded->grp_bond_con.j3_33[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_33[i],bonded->grp_bond_con.j3_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_33;
  for(i=1;i<=bonded->grp_bond_con.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j2_33[i],bonded->grp_bond_con.j3_33[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j2_33[i],bonded->grp_bond_con.j3_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_33;

/*   i) 33 Watts [1-2,1-3,2-3]            */
#ifdef DEBUG
  printf("number of group constrains %8d \n",bonded->grp_bond_con.num_33);
#endif
  for(i=1;i<=bonded->grp_bond_watts.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_watts.j1_33[i],bonded->grp_bond_watts.j2_33[i]);
    jtmp2[(i+ioff)] =  

         MIN(bonded->grp_bond_watts.j1_33[i],bonded->grp_bond_watts.j2_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_watts.num_33;
  for(i=1;i<=bonded->grp_bond_watts.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_watts.j1_33[i],bonded->grp_bond_watts.j3_33[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_watts.j1_33[i],bonded->grp_bond_watts.j3_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_watts.num_33;
  for(i=1;i<=bonded->grp_bond_watts.num_33;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_watts.j2_33[i],bonded->grp_bond_watts.j3_33[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_watts.j2_33[i],bonded->grp_bond_watts.j3_33[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_watts.num_33;

/*   i) 43 [1-2,1-3,1-4]*/
  for(i=1;i<=bonded->grp_bond_con.num_43;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j2_43[i]);

    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j2_43[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_43;
  for(i=1;i<=bonded->grp_bond_con.num_43;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j3_43[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j3_43[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_43;
  for(i=1;i<=bonded->grp_bond_con.num_43;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j4_43[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_43[i],bonded->grp_bond_con.j4_43[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_43;

/*   i) 46 [1-2,1-3,1-4,2-3,2-4,3-4]*/
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j2_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j2_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j3_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j3_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j4_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j1_46[i],bonded->grp_bond_con.j4_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j2_46[i],bonded->grp_bond_con.j3_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j2_46[i],bonded->grp_bond_con.j3_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j2_46[i],bonded->grp_bond_con.j4_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j2_46[i],bonded->grp_bond_con.j4_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;
  for(i=1;i<=bonded->grp_bond_con.num_46;i++){
    jtmp1[(i+ioff)] =  
         MAX(bonded->grp_bond_con.j3_46[i],bonded->grp_bond_con.j4_46[i]);
    jtmp2[(i+ioff)] =  
         MIN(bonded->grp_bond_con.j3_46[i],bonded->grp_bond_con.j4_46[i]);
  }/*endfor*/
  ioff += bonded->grp_bond_con.num_46;

/*======================================================================*/
/* III) Sort the list                                                   */

  if(num_excl_now>1){
    exl_sort(&num_excl_now,jtmp1,jtmp2,clatoms_info->natm_tot);
  }/*endif*/ 

/*======================================================================*/
/* IV) Form exclusion list                                              */

  now_memory = (
                (num_excl_now)*(sizeof(double)*0 + sizeof(int)*1 )
               +(clatoms_info->natm_tot)*(sizeof(double)*0 + sizeof(int)*2 )
               )*1.e-06;
  *tot_memory += now_memory;
if(error_check_on==1){
  printf("Exclusion allocation: %g Mbytes; Total memory %g Mbytes\n",
           now_memory,*tot_memory);
}/*endif*/

  excl->nlst = num_excl_now;
  excl->j    =  (int *) cmalloc(excl->nlst*sizeof(int))-1;
  excl->num   = (int *) cmalloc(clatoms_info->natm_tot*sizeof(int))-1;
  excl->j_off = (int *) cmalloc(clatoms_info->natm_tot*sizeof(int))-1;   

  for(i=1;i<=clatoms_info->natm_tot;i++){
    excl->num[i]=0;
  }/*endfor*/
  for(i=1;i<=num_excl_now;i++){
    excl->j[i]=jtmp2[i];
    excl->num[jtmp1[i]]++;
  }/*endfor*/
  excl->j_off[1]=0;
  for(i=2;i<=clatoms_info->natm_tot;i++){
    excl->j_off[i] = excl->j_off[(i-1)]+excl->num[(i-1)];
  }/*endfor*/


/*========================================================================*/
/* VI) Sort the exclusions of each particle into ascending order */
/*     using a stupid sort because technical support got lazy     */

  for(i=2;i<=clatoms_info->natm_tot;i++){
    kstart = excl->j_off[i]+1;
    kend   = excl->j_off[i]+excl->num[i];
    for(k=kstart;k<=kend;k++){
      for(j=k+1;j<=kend;j++){
        if(excl->j[j]<excl->j[k] ){
         itemp = excl->j[k];
         excl->j[k] =  excl->j[j];
         excl->j[j] =  itemp;
        }/*endif*/
      }/*endfor*/
    }/*endfor*/
  }/*endfor*/

/*========================================================================*/
/* V) Output */

if(error_check_on==1){
  printf("Total number of exclusions, %d\n\n",num_excl_now);
  PRINT_LINE_DASH;
  printf("Completed exclusion generation\n"); 
  PRINT_LINE_STAR;printf("\n");
}/*endif*/

/*======================================================================*/
/*  VI) Form ewald corrections                                           */

  bonded->ecor.num = 0;

  if(iperd>0){

   if(error_check_on==1){
      PRINT_LINE_STAR;
      printf("Determining ewald corrections \n");
      PRINT_LINE_DASH;printf("\n");
   }/*endif*/

    bonded->ecor.alp_ewd = alp_ewd;
    for(i=1;i<=num_excl_now;i++){
      if((clatoms_info->q[(jtmp1[i])]!= 0.0)
       &&(clatoms_info->q[(jtmp2[i])]!= 0.0)
	  &&(
	  ((clatoms_info->cp_atm_flag[(jtmp1[i])] == 0)
          &&(clatoms_info->cp_atm_flag[(jtmp2[i])] == 0)))
       )
      {	bonded->ecor.num++; }
    }/*endfor*/

    now_memory      = (
                          (bonded->ecor.num)*
			  (sizeof(double)*0 + sizeof(int)*2 )
                                                            )*1.e-06;
    *tot_memory += now_memory;

   if(error_check_on==1){
       printf("Ecorr allocation: %g Mbytes; Total memory %g Mbytes\n",
               now_memory,*tot_memory);
   }/*endif*/

    bonded->ecor.j1 =(int *)cmalloc(bonded->ecor.num*sizeof(int))-1;
    bonded->ecor.j2 =(int *)cmalloc(bonded->ecor.num*sizeof(int))-1;
    bonded->ecor.num = 0;
    for(i=1;i<=num_excl_now;i++){
      if( (clatoms_info->q[(jtmp1[i])]!=0.0)
        &&(clatoms_info->q[(jtmp2[i])]!=0.0)
        &&(
	  ( (clatoms_info->cp_atm_flag[(jtmp1[i])] == 0)
           &&(clatoms_info->cp_atm_flag[(jtmp2[i])] == 0))) 
       )
     {
	bonded->ecor.num++;
	bonded->ecor.j1[(bonded->ecor.num)] = jtmp1[i];
	bonded->ecor.j2[(bonded->ecor.num)] = jtmp2[i]; 
      }/*endif*/
    }/*endfor*/


   if(error_check_on==1){
    printf("Total number of ewald corrections, %d \n\n",bonded->ecor.num);
    PRINT_LINE_DASH;
    printf("Ewald corrections determined\n");
    PRINT_LINE_STAR;printf("\n");
   }/*endif*/

  }/*endif*/

/*======================================================================*/
/*  VI) Form ewald corrections for interaction of CP atoms with         */
/*       classical atoms which were excluded                            */

   if(error_check_on == 1){
      PRINT_LINE_STAR;
      printf("Determining CP-CLASSICAL ewald corrections \n");
      PRINT_LINE_DASH;printf("\n");
   }/*endif*/

 /* CP-CP atoms are never excluded */
 /* Count number of pairs ab initio - ab initio */

	bonded->excl.num_cp = 0;
    for(i=1;i<=num_excl_now;i++){
      if( (clatoms_info->q[(jtmp1[i])]!= 0.0)
        &&(clatoms_info->q[(jtmp2[i])]!= 0.0)
        &&((clatoms_info->cp_atm_flag[(jtmp1[i])] == 1)
         &&(clatoms_info->cp_atm_flag[(jtmp2[i])] == 1))
        )
      {
	bonded->excl.num_cp++;
      }/*endif*/
    }/*endfor*/

    now_memory      = ( (bonded->excl.num_cp)*
			(sizeof(double)*0 + sizeof(int)*2 ))*1.e-06;
    *tot_memory += now_memory;

   if(error_check_on == 1){
       printf("MIXED CP CORR allocation: %g Mbytes; Total memory %g Mbytes\n",
               now_memory,*tot_memory);
   }/*endif*/

    bonded->excl.j1_cp =(int *)cmalloc((bonded->excl.num_cp)*sizeof(int))-1;
    bonded->excl.j2_cp =(int *)cmalloc((bonded->excl.num_cp)*sizeof(int))-1;

    bonded->excl.num_cp = 0;

    for(i=1;i<=num_excl_now;i++){
      if( ((clatoms_info->q[(jtmp1[i])]!= 0.0)
        &&(clatoms_info->q[(jtmp2[i])]!= 0.0))
       &&
         ( ((clatoms_info->cp_atm_flag[(jtmp1[i])] == 1)
            &&(clatoms_info->cp_atm_flag[(jtmp2[i])] == 1))
         ))
      {
	bonded->excl.num_cp++;
	bonded->excl.j1_cp[(bonded->excl.num_cp)] = jtmp1[i];
	bonded->excl.j2_cp[(bonded->excl.num_cp)] = jtmp2[i];
      }/*endif*/
    }/*endfor*/

   if(error_check_on==1){
     printf("Total number of MIXED-CP ewald corrections, %d \n",
                                                      bonded->excl.num_cp);
     PRINT_LINE_DASH;
     printf("MIXED-CP Ewald corrections determined\n");
     PRINT_LINE_STAR;printf("\n");
   }/*endif*/

  cfree(&jtmp1[1]);
  cfree(&jtmp2[1]);

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void splin_ecor(ECOR *ecor,EWALD *ewald,CLATOMS_POS *clatoms_pos,int pi_beads,
                int error_check_on,double *tot_memory)

/*==========================================================================*/ 
{ /*begin routine*/
/*==========================================================================*/
/*           Local Variables                                                */

  int i,ip,iii,k;
  double rmin,rmax,dr,dri,r12,r12i;
  double dx,dy,dz,r2;
  double talp2,palp,alp_ewd2;
  double ralp,eee,tt,temp,gerf,dgerf,dvecor;
  double pk,dk,kcut,kcut_res,kmax;
  double bess_j0,bess_j1,chi,dchi,arg;
  double *know,*fknow,*anode,*weight;
  int num_mem;
  double now_memory,self_erfc;

  double p  = 0.3614;
  double e1 = 0.2041422096422003, e2 = 0.1997535956961481;
  double e3 = 0.2213176596405576, e4 = 0.03360430734640255;
  double e5 = 0.4732592578721755, e6 =-0.509078520069735;
  double e7 = 0.6772631491947646, e8 =-0.369912979092217;
  double e9 = 0.06965131976970335;

  double de1,de2,de3;
  double de4,de5,de6;
  double de7,de8,de9;
  double ten;

/*           Local Pointers                                   */
  int nsplin      = ecor->nsplin;
  int nsplin_mal  = ecor->nsplin;
  int nk_use      = ecor->nsplin;
  int *j1         = ecor->j1;
  int *j2         = ecor->j2;
  int num         = ecor->num;
  int nktot_res   = ecor->nktot_res;
  double alp_ewd  = ecor->alp_ewd;
  double ecut     = ecor->ecut;
  double ecut_res = ecor->ecut_res;
  double *x,*y,*z;
  double *cv0,*cdv0,*cdv0_res;

  ten = 10.0/BOHR;
  de1 = 1.0*e1; de2 = 2.0*e2; de3 = 3.0*e3;
  de4 = 4.0*e4; de5 = 5.0*e5; de6 = 6.0*e6;
  de7 = 7.0*e7; de8 = 8.0*e8; de9 = 9.0*e9;

/*==========================================================================*/ 
/* 0) Output */

  if(error_check_on==1){
    PRINT_LINE_STAR;
    printf("Setup the finite k-space ewald corrections \n");
    PRINT_LINE_DASH;printf("\n");
  }/*endif*/

/*==========================================================================*/ 
/* I) Set the real space range and malloc the ecorr memory                  */

  rmax  = 0.0;
  for(ip=1;ip<=pi_beads;ip++){
   x = clatoms_pos[ip].x;
   y = clatoms_pos[ip].y;
   z = clatoms_pos[ip].z;
   for(i=1;i<=num;i++){
     dx   = x[j1[i]]-x[j2[i]];
     dy   = y[j1[i]]-y[j2[i]];
     dz   = z[j1[i]]-z[j2[i]];
     r2   = dx*dx+dy*dy+dz*dz;
     rmax = MAX(r2,rmax);
   }/*endfor*/
  }/*endfor*/
  rmax += 2.0/BOHR;      /* safety */
  rmax  = MAX(ten,rmax); /* safety */
  rmin  = 0.1;
  dr    = (rmax-rmin) /(double)(nsplin - 5);
  dri   = 1.0/dr;

  ecor->rmin_spl = rmin;
  ecor->rmax_spl = rmax;
  ecor->dr_spl   = dr;
  ecor->dri_spl  = dri;

  if((nsplin_mal % 2)==0){nsplin_mal++;}
  ecor->cv0      = (double *) cmalloc(nsplin_mal*sizeof(double))-1;
  ecor->cdv0     = (double *) cmalloc(nsplin_mal*sizeof(double))-1;
  cv0            = ecor->cv0;
  cdv0           = ecor->cdv0;
  num_mem        = 2*nsplin_mal;
  if(nktot_res>0){
    ecor->cdv0_res = (double *) cmalloc(nsplin_mal*sizeof(double))-1;
    cdv0_res       = ecor->cdv0_res;
    num_mem        = 3*nsplin_mal;
  }/*endif*/

/*==========================================================================*/ 
/* II) Get the infinte accuracy k-space ecorr */

  alp_ewd2 = alp_ewd*alp_ewd;
  talp2    = 2.0*(alp_ewd)*(alp_ewd);
  palp     = p*(alp_ewd);

  for (i = 1; i <= nsplin; ++i) {
    r12     = dr *(double)(i-3) + rmin;
    r12i    = 1.0/r12;
    ralp    = r12*(alp_ewd);
    eee     = exp(-ralp*ralp);
    tt      = 1.0/(1.0+p*ralp);
    temp    = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                           +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
    gerf    = 1.0-temp;
    dgerf   = ((((((((de9*tt+de8)*tt+de7)*tt+de6)*tt+de5)*tt
                            +de4)*tt+de3)*tt+de2)*tt+de1)*tt*tt*eee*palp
                            +talp2*temp*r12;
    cv0[i]  = -gerf*r12i;
    dvecor  = -(dgerf*r12i-gerf/(r12*r12));
    cdv0[i] = -dvecor*r12i;
  }/*endfor*/

  if(nktot_res>0){
    for (i = 1; i <= nsplin; ++i) {
      cdv0_res[i] = cdv0[i];
    }/*endfor*/
  }/*endif*/

/*==========================================================================*/ 
/* III) Add in the k-space cutoff dependent corrections                     */

 /*-------------------------------------------------------------*/
 /* i) Malloc the memory and get the gaussian quadrature points */
  nk_use = 128;
  know   = (double *) cmalloc(nk_use*sizeof(double))-1;
  fknow  = (double *) cmalloc(nk_use*sizeof(double))-1;
  anode  = (double *) cmalloc(nk_use*sizeof(double))-1;
  weight = (double *) cmalloc(nk_use*sizeof(double))-1;
#include "../proto_defs/weights_nodes_128.h"

 /*-----------------------------------------------------------*/
 /* ii) Compute the self term  correction                     */
  kcut = sqrt(2.0*ecut);
  arg  = kcut/(2.0*alp_ewd);
  eee  = exp(-arg*arg);
  tt   = 1.0/(1.0+p*arg);
  self_erfc = ((((((((e9*tt+e8)*tt+e7)*tt+e6)*tt+e5)*tt
                             +e4)*tt+e3)*tt+e2)*tt+e1)*tt*eee;
  ewald->self_erf  = 1.0-self_erfc;

 /*-----------------------------------------------------------*/
 /* iii) Compute the real space correction for standard ewald */
  kmax = 15.0*alp_ewd;
  if(kmax>kcut){
    dk   = (kmax-kcut);
    pk   = (kmax+kcut);
    for(k=1;k<=nk_use;k++){
      know[k]  = 0.5*(anode[k]*dk + pk);
      fknow[k] = (dk/M_PI)*exp(-0.25*know[k]*know[k]/alp_ewd2)*weight[k];
    }/*endfor*/
    for (i = 1; i <= nsplin; ++i) {
      r12     = dr *(double)(i-3) + rmin;
      chi = 0.0;
      dchi = 0.0;
      for(k=1;k<=nk_use;k++){
        arg      = know[k]*r12;
        bess_j0  = (sin(arg)/arg);
        bess_j1  = (sin(arg)/arg - cos(arg))/arg;
        chi     += (fknow[k]*bess_j0);
        dchi    += (know[k]*fknow[k]*bess_j1/r12);
      }/*endfor*/
      cv0[i]  += chi;
      cdv0[i] += dchi;
    }/*endfor:real space points*/
  }/*endif:kcut<kmax*/

 /*-----------------------------------------------------------*/
 /* iv) Compute the Real space correction for respa ewald    */
  if(nktot_res>0){
    kcut_res = sqrt(2.0*ecut_res);
    if(kmax>kcut_res){
      dk   = (kmax-kcut_res);
      pk   = (kmax+kcut_res);
      for(k=1;k<=nk_use;k++){
        know[k]  = 0.5*(anode[k]*dk + pk);
        fknow[k] = (dk/M_PI)*exp(-0.25*know[k]*know[k]/alp_ewd2)
                   *weight[k];
      }/*endfor*/
      for (i = 1; i <= nsplin; ++i) {
        r12     = dr *(double)(i-3) + rmin;
        dchi = 0.0;
        for(k=1;k<=nk_use;k++){
          arg     = know[k]*r12;
          bess_j1 = (sin(arg)/arg - cos(arg))/arg;
          dchi   += (know[k]*fknow[k]*bess_j1/r12);
        }/*endfor*/
        cdv0_res[i] += dchi;
      }/*endfor:real space points*/
    }/*endif:kcut_res<kmax*/
  }/*endif:irespa_on*/

 /*-----------------------------------------------------------*/
 /* v) Free the memory                                       */
  cfree(&know[1]);
  cfree(&fknow[1]);
  cfree(&anode[1]);
  cfree(&weight[1]);

/*==========================================================================*/ 
/* IV) Memory summary and Output */

  now_memory   = ( sizeof(double)*num_mem )*1.e-06;
  *tot_memory += now_memory;

  if(error_check_on==1){

    printf("The Ewald sum convergence factor, erfc(k_cut/2alp_ewd), is : %g\n",
            self_erfc);
    printf("Ecorr-spline allocation: %g Mbytes; Total memory %g Mbytes\n",
            now_memory,*tot_memory);

    printf("\n");PRINT_LINE_DASH;
    printf("Completed finite k-space ewald correction setup\n");
    PRINT_LINE_STAR;

  }/*endif*/

/*--------------------------------------------------------------------------*/
  }/*end routine */
/*==========================================================================*/


