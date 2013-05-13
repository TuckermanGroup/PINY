
#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_real_space_local.h"


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*==========================================================================*/
/*                      Brief description                                   */

/* lnk_list[ia][ib][ic][iatm] stores the indices of iatmth atm in           */
/*                           a cell labled ia,ib,ic                         */
/*                           (a pipette cut out of the full system cell).   */
/*                           It is packed with zeros for                    */
/*                           iatm greater than the number of atoms          */
/*                           in a given cell for computational              */
/*                           efficiency. Zero index atoms give              */
/*                           rise to ``null'' interactions that are         */
/*                           removed with logic                             */
/* ishft_a[i],ishft_b[i],ishft_c[i]                                      */
/*                           are vectors containing the allowed sets of     */
/*                           shifts along a,b,c                             */
/* iexl_check[i]             indicates whether or not exclusions needed     */
/*                           to be removed in a particular shift set        */
/* When the lnk_list is shifted interactions are generated between cells    */
/* of the shifted list and old unshifted list.                              */
/* These interactions are sent to the force routine in bundles of nlen.     */
/* Break points where memory conflicts can occur when summing the forces    */
/* are tablulated.                                                          */
/*==========================================================================*/

void force_lnklst(CLATOMS_INFO *clatoms_info,CLATOMS_POS *clatoms_pos,
                  FOR_SCR *for_scr,ATOMMAPS *atommaps,
                  CELL *cell,PTENS *ptens,INTERACT *interact,
                  ENERGY_CTRL *energy_ctrl, 
                  int ncell_a,int ncell_b,int ncell_c,
                  int natm_cell_max,list_int lnk_list[],
                  int nshft_lnk,double shft_wght[],
                  int ishft_a[],int ishft_b[],int ishft_c[],int iexl_chk[],
                  EXCL *excl,INTRA_SCR *intra_scr,double *vreal,double *vvdw,
                  double *vcoul,CLASS_COMM_FORC_PKG *class_comm_forc_pkg)

{/*Begin Routine*/
  /*=======================================================================*/
  /*          Local Variable declarations                                   */
  
  int iii;
  int itemp,nlen_use;              /* value of nlen to be used           */
  int ncell_tot, nlist_tot;        /* tot number of cells, size of list  */
  int natm_shft;                   /* number of atom shifts              */
  int iadd_shft;                   /* atm shift flag                     */
  int ishft, iatm_shft;           /* cell and atom shift loops indicies  */
  int ja,jb,jc,jatm;               /* the present shift                  */
  int icheck;                      /* excl check on                      */
  int intact_now, intact_tmp;      /* interaction temporaries            */
  int intact_tot, num_call;        /* tot interactions,num of force calls*/
  int i;                           /* loop counter                       */
  int imin,jmax;                   /* min/max atm ind in a interact pair */
  int nact;                        /* interaction is good flag           */
  int iact;                        /* excl interaction loop              */
  int lower,upper;                 /* limits on excl interaction loop    */
  int ifor_call;                   /* force call flag                    */
  int jind_off;                    /* integer offset                     */
  int ncell_af,ncell_bf,ncell_cf,natm_cellf;
  int mtemp,ntemp,ioff,natm_tot; 
  int intact_init,intact_half,nlist_tot2;
  int lst_typ=3;
  int nskip,ishft_st;
/* Local Pointers */
  int *for_scr_intact    = &(for_scr->intact);
  int *for_scr_num_brk   = &(for_scr->num_brk);
  int *for_scr_num_brk_i = for_scr->num_brk_i;
  int *for_scr_i_index   = for_scr->i_index;  
  int *for_scr_j_index   = for_scr->j_index;  
  list_int *for_scr_i_lnk     = for_scr->i_lnk;
  list_int *for_scr_j_lnk     = for_scr->j_lnk;
  int *excl_j            = excl->j;
  int *excl_j_off        = excl->j_off;
  int *excl_num          = excl->num;
  double *for_scr_wght_lnk   = &(for_scr->wght_lnk);
  int *for_scr_iexcl         = for_scr->iexcl;
  int myid_forc          = class_comm_forc_pkg->myid;
  int np_forc            = class_comm_forc_pkg->num_proc;

  /*========================================================================*/
  /* I) Get a nice value of nlen: a value divisible by natm_tot or vice versa*/
  
  natm_tot =  clatoms_info->natm_tot;
  if((for_scr->nlen)>(natm_tot)){
    nlen_use =  (for_scr->nlen) - ((for_scr->nlen)% (natm_tot));
  }else{
    itemp    = (natm_tot)/(for_scr->nlen) + 1;
    nlen_use = (natm_tot)/itemp;
  }/*endif*/
  
  /*=======================================================================*/
  /* II) Get some useful constants*/
  
  ncell_tot = (ncell_a)*(ncell_b)*(ncell_c);
  nlist_tot = ncell_tot*(natm_cell_max);

  /*=======================================================================*/
  /* III) Inititialize           */
  
  (*for_scr_intact)  = 0;   
  (*for_scr_num_brk) = 0;
  num_call           = 0;
  intact_tot         = 0;
  ncell_af = ncell_a;
  ncell_bf = ncell_b;
  ncell_cf = ncell_c;
  natm_cellf = natm_cell_max;

  intact_init = natm_tot;
  intact_half = 0; 
  nlist_tot2  = 4*natm_cellf*ncell_af*ncell_bf*ncell_cf
              - 4*ncell_af*ncell_bf*ncell_cf
              - 2*ncell_af*ncell_bf
              - ncell_af;
  for(i=1;i<=natm_tot;i++){
    if(for_scr_iexcl[i]<=nlist_tot2){intact_half=i;}
  }/*endfor*/
  
/*=========================================================================*/
/* IV) Loop over the allowed cell shifts: ja,jb,jc                         */
  
  ishft_st = 1 + myid_forc;
  nskip    = np_forc;
  for(ishft=ishft_st;ishft<=nshft_lnk;ishft+=nskip){
    ja                  = (ishft_a)[ishft];
    jb                  = (ishft_b)[ishft];
    jc                  = (ishft_c)[ishft];
    icheck              = (iexl_chk)[ishft];
    (*for_scr_wght_lnk) = (shft_wght)[ishft];
    
/*=======================================================================*/
/* V) Loop over all the atm shifts: jatm                                 */
/*     Bookeeping for the first cell shift to take care                  */
/*           of self interactions                                        */
    
    natm_shft = (natm_cell_max);
    iadd_shft = 0;
    if(ishft==1){
      iadd_shft  = (((natm_cell_max)-1) % 2);
      natm_shft  = ((natm_cell_max)-1)/2 + iadd_shft;
      /*endif*/}
    for(iatm_shft=1;iatm_shft<=natm_shft;iatm_shft++){
      jatm = iatm_shft;

/*==========================================================================*/
/*  B) Book keeping for last atm shift of 1st cell shift:                  */
/*         On the last iatm_shft of the 1st cell shift get rid of           */
/*         repeated interactions that occur if natm_cell_max is even.       */
/*         This simply means taking 1/2 of the interactions found           */

         intact_now = intact_init;
         if((ishft==1)&&(iadd_shft==1)&&(iatm_shft==natm_shft))
                 {intact_now = intact_half;}

/*==========================================================================*/
/*   C) Eliminate null interactions caused by packing the lnk_list          */
/*       with zero index atoms                                              */

         intact_tmp = intact_now;
         intact_now = 0;
         ioff =  ja + jb*ncell_af*2 + jc*ncell_af*ncell_bf*4 
               + jatm*ncell_af*ncell_bf*ncell_cf*8;
         for(i=1;i<=intact_tmp;i++){
           ntemp = for_scr_iexcl[i]+ioff;
           if(lnk_list[ntemp]!=0){
               intact_now += 1;
               mtemp = for_scr_iexcl[i];
               for_scr_i_lnk[intact_now] = lnk_list[mtemp];
               for_scr_j_lnk[intact_now] = lnk_list[ntemp];
           /*endif*/}
         /*endfor*/}

      /*====================================================================*/
      /*   E) Eliminate standard exclusions in shifts that might have them. */
      /*       Cells that are far enough away in space won't have any excls.*/
      
      if(icheck==1){
        intact_tmp = intact_now;
        intact_now  = 0;
        for(i=1;i<=intact_tmp;i++){
          imin  = MIN( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
          jmax  = MAX( (for_scr_i_lnk)[i],(for_scr_j_lnk)[i] );
          lower = (excl_j_off)[jmax]+1;
          upper = (excl_num)[jmax]+(excl_j_off)[jmax];
          nact  = 1;
          for(iact=lower;iact<=upper;iact++){
            if(imin==(excl_j)[iact]){nact = 0;}
          /*endfor*/}
          if(nact==1){
            intact_now += 1;
            (for_scr_i_lnk)[intact_now] =  (for_scr_i_lnk)[i];
            (for_scr_j_lnk)[intact_now] =  (for_scr_j_lnk)[i];
            /*endif*/}
          /*endfor*/}
        /*endif*/}
      
      /*====================================================================*/
      /* G) Add interactions to the temporary interaction list             */
      /*  in chunks no greater than nlen_use: Keep track of the break points*/
      
      lower = 1;
      upper = 0;
      while(upper<intact_now){
        upper = intact_now;
        if(upper-lower+1+(*for_scr_intact)>nlen_use){
          upper = nlen_use-(*for_scr_intact)+lower-1;
         /*endif*/}
        jind_off = (*for_scr_intact)-lower+1;
        for(i=lower;i<=upper;i++){
          mtemp = (i+jind_off);
          (for_scr_j_index)[mtemp]=(for_scr_j_lnk)[i];
          (for_scr_i_index)[mtemp]=(for_scr_i_lnk)[i];
          /*endfor*/}
        (*for_scr_intact)                       += upper-lower+1;
        (*for_scr_num_brk)                      += 1;
        (for_scr_num_brk_i)[(*for_scr_num_brk)] = upper-lower+1;
        lower                                    = upper+1;
        
        /*===================================================================*/
        /*  H) If enough interactions have accumulated (or you are done)     */
        /*       make a force call and then update and reinitialize          */
        /*       appropriate counters                                        */
        
        ifor_call = 0;
        if( (*for_scr_intact)>=nlen_use){ifor_call=1;}
        if(  (ishft==(nshft_lnk))  &&
             (iatm_shft==natm_shft)&&
             (upper==intact_now)   &&
            ((*for_scr_intact)!=0) ) {ifor_call=1;}
        if(ifor_call==1){
          force_npol(clatoms_info,clatoms_pos,
                     for_scr,atommaps,
                     cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                     vvdw,vcoul,num_call,lst_typ); 
          num_call          += 1;
          intact_tot        += (*for_scr_intact);
          (*for_scr_intact)  = 0;   
          (*for_scr_num_brk) = 0;
          /*endif*/}         
        /*endwhile*/}
      
/*====================================================================*/
      
      /*endfor:iatm_shft*/}
    /*endfor:ishft*/}
  
/*====================================================================*/
/* VI) Final dump */

   if((*for_scr_intact)>0){
      force_npol(clatoms_info,clatoms_pos,
                 for_scr,atommaps,
                 cell,ptens,interact,energy_ctrl,intra_scr,vreal,
                 vvdw,vcoul,num_call,lst_typ); 
       num_call          += 1;
       intact_tot        += (*for_scr_intact);
    }/*endif*/

  /*------------------------------------------------------------------------*/
  /*end routine */}
/*==========================================================================*/





