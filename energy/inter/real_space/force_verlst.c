/*==========================================================================*/
/*                      Brief description                                   */

/* The (nter)[i] atoms that interact with atom, i,                */
/* are stored in (jter)[i]                                        */
/* starting at element 1+(jver_off[i])]                           */
/* ending   at element (nter)[i]+(jver_off[i])]         */ 
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated                                                           */

/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void force_verlst(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  FOR_SCR *for_scr,ATOMMAPS *atommaps,
                  CELL *cell,PTENS *ptens,INTERACT *interact,
                  ENERGY_CTRL *energy_ctrl, 
                  int nter[],list_int jter[],int jver_off[],
                  INTRA_SCR *intra_scr,double *vreal,double *vvdw,
                double *vcoul)
     
{/*Begin Routine*/
  /*========================================================================*/
  /*          Local Variable declarations               */
  
  int ipart, jpart;           /* ith particle, jth particle          */
  int intact_tot;             /* total number interactions           */
  int intact_save;             /* total number interactions           */
  int num_call;               /*  number of force calls              */
  int num_call_save;          /*  number of force calls              */
  int lower,upper;            /* lower and upper limits on for loop  */
  int ifor_call;              /* force call flag                     */
  int joff, jind_off;     /* off sets                            */
  int iii,itemp;
  int natm_tot,nlen;
  int lst_typ=1;

/* Local Pointers */
  int *for_scr_intact    = &(for_scr->intact);
  int *for_scr_num_brk   = &(for_scr->num_brk);
  int *for_scr_num_brk_i = for_scr->num_brk_i;
  int *for_scr_i_index   = for_scr->i_index;  
  int *for_scr_j_index   = for_scr->j_index;  

  /*========================================================================*/
  /*  I) Count the interactions and determine the number of force calls */

  intact_save = 0;
  natm_tot    =  clatoms_info->natm_tot;
  nlen        =  (for_scr->nlen);
  for(ipart=1;ipart<=natm_tot;ipart++){
    intact_save += nter[ipart];
  }/*endfor*/
  num_call_save = intact_save/(nlen);
  if( (intact_save % (nlen))!=0 ){num_call_save += 1;}
  
  /*========================================================================*/
  /*  II) Initialize */
  
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call         = 0;
  intact_tot       = 0;
  
  /*========================================================================*/
  /* III) Loop over the possible interactions:                              */
  
  for(ipart=2;ipart<=(natm_tot);ipart++){
    
    /*=======================================================================*/
    /*  A) Add interactions to the force routine interaction list            */
    /*     in chunks no greater than nlen: Keep track of the break points    */
    
    lower  = 1;
    upper  = 0;
    while(upper!=nter[ipart]){
      upper = nter[ipart];
      if(upper-lower+1+(*for_scr_intact)>(nlen)){
       upper = nlen - (*for_scr_intact) + lower-1;
      }/*endif*/
      jind_off = (*for_scr_intact)-lower+1;
      joff = jver_off[ipart];
      for(jpart=lower;jpart<=upper;jpart++){
       for_scr_j_index[(jpart+jind_off)]=jter[jpart+joff];
       for_scr_i_index[(jpart+jind_off)]=ipart;
      }/*endfor*/
      (*for_scr_intact) += upper-lower+1;
      (*for_scr_num_brk)++;
      for_scr_num_brk_i[(*for_scr_num_brk)] = upper-lower+1;
      lower = upper+1;
      
      /*=====================================================================*/
      /*  B) If enough interactions have accumulated (or you are done)       */
      /*       make a force call and then update and reinitialize            */
      /*       appropriate counters                                          */
      
      itemp = MIN(((nlen)-(*for_scr_intact)),1);
      if(!(itemp)){
       force_npol(clatoms_info,clatoms_pos,
                  for_scr,atommaps,
                  cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                  vvdw,vcoul,num_call,lst_typ);

       num_call        += 1;
       intact_tot      += (*for_scr_intact);
       (*for_scr_intact)  = 0;   
       (*for_scr_num_brk) = 0;
     }/*endif*/


     if((*for_scr_num_brk)>=nlen){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Internal Error: If statement to dump force not performed\n ");
       printf("in routine for_verlst \n ");
       printf("Contact product support \n ");
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
       exit(1);
      }/*endif*/
    }/*endwhile*/
  }/*endfor*/


  if((*for_scr_intact)!=0) {
     force_npol(clatoms_info,clatoms_pos,
                for_scr,atommaps,
                cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                vvdw,vcoul,num_call,lst_typ);
      num_call        += 1;
      intact_tot      += (*for_scr_intact);
      (*for_scr_intact)  = 0;   
      (*for_scr_num_brk) = 0;
  }/*endif*/
  
  /*========================================================================*/
  /* III) Check                */
  
  if(intact_tot!=intact_save){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error:\n ");
    printf("Incorrect number of interactions calculated\n ");
    printf("   in routine for_verlst \n ");
    printf("%d vs %d \n",intact_tot,intact_save);
    printf("Contact product support \n ");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
    /*endif*/}
   if(num_call!=num_call_save){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Internal Error:\n ");
     printf("Incorrect number of force calls\n ");
     printf("   in routine for_verlst \n ");
     printf("%d vs %d \n",num_call,num_call_save);
     printf("Contact product support \n ");
     fflush(stdout);
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     exit(1);
     /*endif*/}
  
  /*------------------------------------------------------------------------*/
}/*end routine */
/*==========================================================================*/




