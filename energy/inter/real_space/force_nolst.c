/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                      Brief description                                   */

/* All possible interactions are sent to the force routine in bundles of nlen*/
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated                                                           */

/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"


/*==========================================================================*/
void force_nolst(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                 FOR_SCR *for_scr,ATOMMAPS *atommaps,
                 CELL *cell,PTENS *ptens,INTERACT *interact,
                 ENERGY_CTRL *energy_ctrl, NBR_LIST *nbr_list, EXCL *excl,
                 INTRA_SCR *intra_scr, double *vreal,double *vvdw,
                 double *vcoul,CLASS_COMM_FORC_PKG *class_comm_forc_pkg)
{
  /*======================================================================*/
  /*            Local Variable declarations               */

  int ipart, jpart;           /* ith particle, jth particle          */
  int intact_tot;             /* total number interactions           */
  int intact_save;            /* total number interactions           */
  int num_call;               /*  number of force calls              */
  int num_call_save;          /*  number of force calls              */
  int lower,upper;            /* lower and upper limits on for loop  */
  int ifor_call;              /* force call flag                     */
  int intact_now;             /* interactions now                    */
  int jind_off;
  int nlen,natm_tot;
  int excl_nlst = excl->nlst;
  int lst_typ=1;
  int nskip,ipart_st;
/* Local Pointers */
  int *for_scr_intact    = &(for_scr->intact);
  int *for_scr_num_brk   = &(for_scr->num_brk);
  int *for_scr_num_brk_i = for_scr->num_brk_i;
  int *for_scr_i_index   = for_scr->i_index;  
  int *for_scr_j_index   = for_scr->j_index;  
  int *for_scr_iexcl     = for_scr->iexcl;
  int *excl_j            = excl->j;
  int *excl_j_off        = excl->j_off;
  int *excl_num          = excl->num;
  int myid_forc          = class_comm_forc_pkg->myid;
  int np_forc            = class_comm_forc_pkg->num_proc;

  /*========================================================================*/
  /* I) Count the interactions                                              */

  natm_tot =  clatoms_info->natm_tot;
  nlen     =  for_scr->nlen;

  ipart_st = 2+myid_forc;
  nskip    = np_forc;
  intact_save  = 0;
  for(ipart=ipart_st;ipart<=natm_tot;ipart+=nskip){
    intact_save += (ipart-1-excl_num[ipart]);
  }
  num_call_save = intact_save/(nlen);
  if( (intact_save % (nlen))!=0 ){num_call_save += 1;}
  
  /*=======================================================================*/
  /* II) Initialize                                                        */

  (*for_scr_intact)   = 0;
  (*for_scr_num_brk)  = 0;
  num_call          = 0;
  intact_tot        = 0;

  /*=======================================================================*/
  /* III) Loop over the possible interactions:                             */

  ipart_st = 2+myid_forc;
  nskip    = np_forc;
  for(ipart=ipart_st;ipart<=natm_tot;ipart+=nskip){
    
    /*=====================================================================*/
    /* A) Check for exclusions and create an interaction list              */

    if(excl_nlst!=0){
       for(jpart=1;jpart<=ipart-1;jpart++){
         (for_scr_iexcl)[jpart] = 1;
       }/*endfor*/
       lower = (excl_j_off)[ipart]+1;
       upper = (excl_num)[ipart]+(excl_j_off)[ipart];
       for(jpart=lower;jpart<=upper;jpart++){
         (for_scr_iexcl)[(excl_j)[jpart]] = 0;
       }/*endfor*/
       intact_now = 0;
       for(jpart=1;jpart<=ipart-1;jpart++){
         if((for_scr_iexcl)[jpart]==1){
              intact_now                  += 1;
              (for_scr_iexcl)[intact_now] = jpart;
         }/*endif*/
       }/*endfor*/
    }else{
       for(jpart=1;jpart<=ipart-1;jpart++){
           for_scr_iexcl[jpart] = jpart;
       }/*endfor*/
       intact_now = (ipart-1);
    }/*endif*/ 
    /*====================================================================*/
    /*  B) Add interactions to the force routine interaction list           */
    /*     in chunks no greater than nlen. Keep track of the break points   */

    lower  = 1;
    upper  = 0;
    while(upper!=intact_now){
      upper = intact_now;
      if(upper-lower+1+(*for_scr_intact)>nlen){
        upper = nlen-(*for_scr_intact)+lower-1;
      }/*endif*/
      jind_off = (*for_scr_intact)-lower+1;
      for(jpart=lower;jpart<=upper;jpart++){
        (for_scr_j_index)[(jpart+jind_off)]=(for_scr_iexcl)[jpart];
        (for_scr_i_index)[(jpart+jind_off)]=ipart;
      }/*endfor*/
      (*for_scr_intact)                       += upper-lower+1;
      (*for_scr_num_brk)                      += 1;
      (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
      lower                                    = upper+1;

      /*====================================================================*/
      /*  C) If enough interactions have accumulated (or you are done)      */
      /*       make a force call and then update and reinitialize           */
      /*       appropriate counters                                         */

      ifor_call = 0;
      if(((*for_scr_intact)==(nlen))){ifor_call=1;}
      if( (ipart==(natm_tot))&&
         (upper==intact_now)&&
         ((*for_scr_intact)!=0)   ){ifor_call=1;}
      if(ifor_call==1){
        force_npol(clatoms_info,clatoms_pos,
                   for_scr,atommaps,
                   cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                   vvdw,vcoul,num_call,lst_typ);
        num_call         += 1;
        intact_tot       += (*for_scr_intact);
        (*for_scr_intact)   = 0;
        (*for_scr_num_brk)  = 0;
      }/*endif*/
    }/*endwhile*/
  }/*endfor*/

  /*======================================================================*/
  /* Final dump */

  if((*for_scr_intact)!=0){
        force_npol(clatoms_info,clatoms_pos,
                   for_scr,atommaps,
                   cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                   vvdw,vcoul,num_call,lst_typ);
        num_call         += 1;
        intact_tot       += (*for_scr_intact);
        (*for_scr_intact)   = 0;
        (*for_scr_num_brk)  = 0;
  }/*endif*/

  /*======================================================================*/
  /* III) Check                */

  if(intact_tot!=intact_save){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error:\n ");
    printf("Incorrect number of interactions calculated\n ");
    printf("   in routine for_nolst \n ");
    printf("%d vs %d \n",intact_tot,intact_save);
    printf("Contact product support \n ");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
  if(num_call!=num_call_save){
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    printf("Internal Error:\n ");
    printf("Incorrect number of force calls\n ");
    printf("   in routine for_nolst \n ");
    printf("%d vs %d \n",num_call,num_call_save);
    printf("Contact product support \n ");
    printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/

}/*end routine */
/*==========================================================================*/










