
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  set_mol_bend_ind:                                                       */
/*==========================================================================*/

void fetch_resbond_prm(RESBOND_PARSE *resbond_parse,BOND_SITE *bond_site,
                         int jmol_typ,int iresidue,int natm_1res_now,
                         int natm_tot)

/*========================================================================*/
/*     Begin routine                                                      */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  NAME site;
  int i,j,k,m,ires_bond,index;
  int nres_bond_2,nres_bond_1;
  int ioff,ifound_site,num_root,iii;
  int num_1back,num_2back[MAX_VALENCE1],num_2back_tot,num_2back_val;
  int ichk_1back[MAX_VALENCE1];
  int ichk_2back[MAX_VALENCE1][MAX_VALENCE1];

/*========================================================================*/
/*  I) Loop over cases where this residue occurs first                    */

  nres_bond_1   = resbond_parse->nres_bond1_jres[iresidue]; 
  ioff          = resbond_parse->res_bond1_off[iresidue];
  for(ires_bond=1;ires_bond<=nres_bond_1;ires_bond++){
    index     = resbond_parse->res_bond1_index[(ires_bond+ioff)]; 
    strcpy(site,resbond_parse->resbond_prm[index].res1_site);
/*------------------------------------------------------------------------*/
/* A) Initialize number of bonds back                                     */
    ifound_site = num_root=0;
    num_1back=0;num_2back_tot=0;
    for(k=1;k<=MAX_VALENCE;k++){ichk_1back[k]=0;num_2back[k]=0;}
    for(k=1;k<=MAX_VALENCE;k++){
     for(m=1;m<=MAX_VALENCE;m++){ichk_2back[m][k]=0;}}
/*------------------------------------------------------------------------*/
/*B) Loop over all bond sites and  all atoms in the residue               */
/*       and check for bond site matches.                                 */
    for(i=1;i<=natm_1res_now;i++){      
      for(j=1;j<=MAX_BOND_SITE;j++){
       if(strcasecmp(site,bond_site[i].bond_site_name[j])==0){
        ifound_site=1;
#ifdef DEBUG
                  printf("%s %d\tj=%d\n",site,i,j);
                  printf("bond_site[%d].bond_site_name[%d] = %s %s\n",i,j,
                          bond_site[i].bond_site_name[j],site);
                  printf("bond_site[%d].branch_1[%d] = %d\n",i,j,
                          bond_site[i].branch_1[j]);
                  printf("bond_site[%d].branch_2[%d] = %d\n",i,j,
                          bond_site[i].branch_2[j]);
                  printf("bond_site[%d].valence = %d\n",i,
                          bond_site[i].valence);
                  scanf("%d",&iii);
#endif
     /*--------------------------------------------------------------------*/
     /* 1) Find the index of the the bond site atom (root atom)            */
     /*    Assign # of atoms in primary branch equal to root valence-1     */
     /*    Initialize index checker                                        */ 
         if(bond_site[i].branch_1[j]==0 && bond_site[i].branch_2[j]==0){
           num_root++; /* root atom found */
           resbond_parse->resbond_prm[index].iatm_res1_site  = i+natm_tot;
           resbond_parse->resbond_prm[index].natm_res1_1back = 
                                                    bond_site[i].valence-1;
           for(k=1;k<=4;k++){
             resbond_parse->resbond_prm[index].improp_res1[k] = 
                          bond_site[i].improper_ind[k];
           }/*endfor*/
         }/*endif*/
     /*--------------------------------------------------------------------*/
     /* 2) Get the indices of atoms 1 bond back (kth primary branch atoms) */
     /*    Assign # of atoms in the secondary branch off the kth atom      */
     /*    in the primary branch equal to the kth atom's valence -1        */
         if(bond_site[i].branch_1[j]>0 && bond_site[i].branch_2[j]==0){
           k = bond_site[i].branch_1[j]; /* kth atom in primary branch */
           num_1back++;ichk_1back[k]+=1;            
           resbond_parse->resbond_prm[index].iatm_res1_1back[k] = i+natm_tot;
           resbond_parse->resbond_prm[index].natm_res1_2back[k] =
                                                    bond_site[i].valence-1;
         }/*endif*/
      
     /*--------------------------------------------------------------------*/
     /* 3) Get the indices of atoms 2 bonds back (mth secdry branch atoms) */
         if(bond_site[i].branch_1[j]>0 && bond_site[i].branch_2[j]>0){
           k = bond_site[i].branch_1[j]; /*kth primary branch  */
           m = bond_site[i].branch_2[j]; /*mth atm in secondary branch */
           num_2back[k]++;ichk_2back[k][m]++;num_2back_tot++;
           resbond_parse->resbond_prm[index].iatm_res1_2back[k][m]=i+natm_tot;
         } /*endif*/
     /*--------------------------------------------------------------------*/
       } /*endif:bond sites match*/
      } /*endfor: each bond site in the atom */
    } /*endfor:each atm in the residue*/

/*------------------------------------------------------------------------*/
/*C) Consistancy checks                                                   */
    if(ifound_site==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
      printf("Bond site %s not found in residue %d\n",site,iresidue);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(num_root!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
      printf("Root atom at bond site %s found %d times in residue %d\n",
                           site,num_root,iresidue);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(num_1back!=resbond_parse->resbond_prm[index].natm_res1_1back){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
       printf("Valence of root atom at bond site %s \n",site);
       printf("not satisfied %d vs %d\n",num_1back,
                     resbond_parse->resbond_prm[index].natm_res1_1back);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
    }
    for(k=1;k<=num_1back;k++){
      if(ichk_1back[k]!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("%dth primary branch atom at bond site %s\n",k,site);
        printf("found %d times\n",ichk_1back[k]);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endfor*/
    num_2back_val=0;
    for(k=1;k<=num_1back;k++){
      num_2back_val += resbond_parse->resbond_prm[index].natm_res1_2back[k];
      if(num_2back[k]!=resbond_parse->resbond_prm[index].natm_res1_2back[k]){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("Valence of primary branch atom %d at bond site %s \n",k,site);
        printf("not satisfied %d vs %d\n",num_2back[k],
                     resbond_parse->resbond_prm[index].natm_res1_2back[k]);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endfor*/
    if(num_2back_val!=num_2back_tot){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("Total Valence of primary branch atoms at bond site %s\n",site);
        printf("not satisfied %d vs %d\n",num_2back_val,num_2back_tot);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }/*endif*/
    for(k=1;k<=num_1back;k++){
     for(m=1;m<=num_2back[k];m++){
        if(ichk_2back[k][m]!=1){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
          printf("%dth secondary branch atom off %dth primary branch atom\n",
                    k,m);
          printf("on bond site %s found %d times \n",site,ichk_2back[k][m]);
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
      }/*endfor:seconary*/
    }/*endfor:primary*/
/*------------------------------------------------------------------------*/
 } /*endfor:bonds with iresidue specified first */

/*========================================================================*/
/*  I) Loop over cases where this residue occurs first                    */

  nres_bond_2   = resbond_parse->nres_bond2_jres[iresidue]; 
  ioff          = resbond_parse->res_bond2_off[iresidue];
  for(ires_bond=1;ires_bond<=nres_bond_2;ires_bond++){      
    index     = resbond_parse->res_bond2_index[(ires_bond+ioff)]; 
    strcpy(site,resbond_parse->resbond_prm[index].res2_site);
/*------------------------------------------------------------------------*/
/* A) Initialize number of bonds back                                     */
    ifound_site = num_root=0;
    num_1back=0;num_2back_tot=0;
    for(k=1;k<=MAX_VALENCE;k++){ichk_1back[k]=0;num_2back[k]=0;}
    for(k=1;k<=MAX_VALENCE;k++){
     for(m=1;m<=MAX_VALENCE;m++){ichk_2back[m][k]=0;}}
/*------------------------------------------------------------------------*/
/*B) Loop over all bond sites and  all atoms in the residue               */
/*       and check for bond site matches.                                 */
    for(i=1;i<=natm_1res_now;i++){      
      for(j=1;j<=MAX_BOND_SITE;j++){
       if(strcasecmp(site,bond_site[i].bond_site_name[j])==0){
        ifound_site=1;
#ifdef DEBUG
                  printf("%s %d\tj=%d\n",site,i,j);
                  printf("bond_site[%d].bond_site_name[%d] = %s %s\n",i,j,
                          bond_site[i].bond_site_name[j],site);
                  printf("bond_site[%d].branch_1[%d] = %d\n",i,j,
                          bond_site[i].branch_1[j]);
                  printf("bond_site[%d].branch_2[%d] = %d\n",i,j,
                          bond_site[i].branch_2[j]);
                  printf("bond_site[%d].valence = %d\n",i,
                          bond_site[i].valence);
                  scanf("%d",&iii);
#endif
     /*--------------------------------------------------------------------*/
     /* 1) Find the index of the the bond site atom (root atom)            */
     /*    Assign # of atoms in primary branch equal to root valence-1     */
     /*    Initialize index checker                                        */ 
         if(bond_site[i].branch_1[j]==0 && bond_site[i].branch_2[j]==0){
           num_root++; /* root atom found */
           resbond_parse->resbond_prm[index].iatm_res2_site  = i+natm_tot;
           resbond_parse->resbond_prm[index].natm_res2_1back = 
                                                    bond_site[i].valence-1;
           for(k=1;k<=4;k++){
             resbond_parse->resbond_prm[index].improp_res2[k] = 
                          bond_site[i].improper_ind[k];
           }/*endfor*/
         }/*endif*/
     /*--------------------------------------------------------------------*/
     /* 2) Get the indices of atoms 1 bond back (kth primary branch atoms) */
     /*    Assign # of atoms in the secondary branch off the kth atom      */
     /*    in the primary branch equal to the kth atom's valence -1        */
         if(bond_site[i].branch_1[j]>0 && bond_site[i].branch_2[j]==0){
           k = bond_site[i].branch_1[j]; /* kth atom in primary branch */
           num_1back++;ichk_1back[k]+=1;            
           resbond_parse->resbond_prm[index].iatm_res2_1back[k] = i+natm_tot;
           resbond_parse->resbond_prm[index].natm_res2_2back[k] =
                                                    bond_site[i].valence-1;
         }/*endif*/
      
     /*--------------------------------------------------------------------*/
     /* 3) Get the indices of atoms 2 bonds back (mth secdry branch atoms) */
         if(bond_site[i].branch_1[j]>0 && bond_site[i].branch_2[j]>0){
           k = bond_site[i].branch_1[j]; /*kth primary branch  */
           m = bond_site[i].branch_2[j]; /*mth atm in secondary branch */
           num_2back[k]++;ichk_2back[k][m]++;num_2back_tot++;
           resbond_parse->resbond_prm[index].iatm_res2_2back[k][m]=i+natm_tot;
         } /*endif*/
     /*--------------------------------------------------------------------*/
       } /*endif:bond sites match*/
      } /*endfor: each bond site in the atom */
    } /*endfor:each atm in the residue*/

/*------------------------------------------------------------------------*/
/*C) Consistancy checks                                                   */
    if(ifound_site==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
      printf("Bond site %s not found in residue %d\n",site,iresidue);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(num_root!=1){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
      printf("Root atom at bond site %s found %d times in residue %d\n",
                           site,num_root,iresidue);
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
    if(num_1back!=resbond_parse->resbond_prm[index].natm_res2_1back){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
       printf("Valence of root atom at bond site %s \n",site);
       printf("not satisfied %d vs %d\n",num_1back,
                     resbond_parse->resbond_prm[index].natm_res2_1back);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
    }
    for(k=1;k<=num_1back;k++){
      if(ichk_1back[k]!=1){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("%dth primary branch atom at bond site %s\n",k,site);
        printf("found %d times\n",ichk_1back[k]);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endfor*/
    num_2back_val=0;
    for(k=1;k<=num_1back;k++){
      num_2back_val += resbond_parse->resbond_prm[index].natm_res2_2back[k];
      if(num_2back[k]!=resbond_parse->resbond_prm[index].natm_res2_2back[k]){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("Valence of primary branch atom %d at bond site %s \n",k,site);
        printf("not satisfied %d vs %d\n",num_2back[k],
                  resbond_parse->resbond_prm[index].natm_res2_2back[k]);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
      }/*endif*/
    }/*endfor*/
    if(num_2back_val!=num_2back_tot){
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
        printf("Total Valence of primary branch atoms at bond site %s\n",site);
        printf("not satisfied %d vs %d\n",num_2back_val,num_2back_tot);
        printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        fflush(stdout);
        exit(1);
    }/*endif*/
    for(k=1;k<=num_1back;k++){
      for(m=1;m<=num_2back[k];m++){
        if(ichk_2back[k][m]!=1){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Error in res_bond %d in molecule %d \n",index,jmol_typ);
          printf("%dth secondary branch atom off %dth primary branch atom\n",
                    k,m);
          printf("on bond site %s found %d times \n",site,ichk_2back[k][m]);
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          fflush(stdout);
          exit(1);
        }/*endif*/
      }/*endfor:seconary*/
    }/*endfor:primary*/
/*------------------------------------------------------------------------*/
 } /*endfor:bonds with iresidue specified first */

/*------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/
