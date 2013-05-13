/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: search_base.c                                */
/*                                                                          */
/* This subprogram matches atom lists with database atom lists              */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_intra_params_local.h"
#include "../proto_defs/proto_search_entry.h"
#include "../proto_defs/proto_search_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define DEBUG_OFF

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  search_base:                                                            */
/*==========================================================================*/

void search_base(int nbase,int nbase2,CATM_LAB *cdata_base,
                 int n, CATM_LAB *cdata,int *igood, int *ifound,
                 int *isearch, int nsearch,int natm_srch,char *filename)
/*==========================================================================*/
/*    Begin Routine */
    {/* begin routine */ 
/*==========================================================================*/
/* Local Variables */
  int i,j,iii;
  int ityp,ibas;
  int comp;
  int *ind_fwd,*ind_bas;
  int ntame;
  int nwild_typ;
  NAME *wild_typ;
  WILD *wild;
  WILD *wild_bas;

/*==========================================================================*/
/* I) Malloc and initialize the wild card control data structures. */
  wild = (WILD *)malloc(n*sizeof(WILD))-1;
  wild_bas = (WILD *)malloc(nbase2*sizeof(WILD))-1;

/* Allow nwild_typ wild cards */
  nwild_typ = 2;
  wild_typ = (NAME *)malloc(nwild_typ*sizeof(NAME))-1;

/* Set the allowed wild cards */
  strcpy(wild_typ[1],"X");
  strcpy(wild_typ[2],"*");

/* Initialize data wild cards (data has none). */
  for(i=1;i<=n;i++){
    wild[i].atm1 = wild[i].atm2 = wild[i].atm3
                 = wild[i].atm4 = wild[i].label = 0;
    wild[i].good = 0;
  }/*endfor*/

/* Initialize data base wild card data. */
  for(i=1;i<=nbase2;i++){
   wild_bas[i].atm1 = wild_bas[i].atm2 = wild_bas[i].atm3
                    = wild_bas[i].atm4 = wild_bas[i].label = 0;
   for(j=1;j<=nwild_typ;j++){
    if(strcasecmp(cdata_base[i].atm1,wild_typ[j])==0){wild_bas[i].atm1 = 1;}
    if(strcasecmp(cdata_base[i].atm2,wild_typ[j])==0){wild_bas[i].atm2 = 1;}
    if(natm_srch > 2){
     if(strcasecmp(cdata_base[i].atm3,wild_typ[j])==0){wild_bas[i].atm3 = 1;}
    }/*endif*/
    if(natm_srch > 3){
      if(strcasecmp(cdata_base[i].atm4,wild_typ[j])==0){wild_bas[i].atm4 = 1;}
    }/*endif*/
    if(strcasecmp(cdata_base[i].label,wild_typ[j])==0){wild_bas[i].label = 1;}
    }/*endfor*/
    wild_bas[i].good = wild_bas[i].atm1 + wild_bas[i].atm2 + wild_bas[i].atm3
                     + wild_bas[i].atm4 + wild_bas[i].label;
  }/*endfor*/

/*==========================================================================*/
/* II) Sort data and base */
/*    A) Malloc the index arrays for sorting */
  ind_fwd = (int *)malloc(n*sizeof(int))-1;
  ind_bas = (int *)malloc(nbase2*sizeof(int))-1;

/*    B) Sort. */
  ind_fwd[1] = 1;
  if(n > 1){sort_atm_list(cdata,n,ind_fwd,wild,natm_srch);}
  ind_bas[1] = 1;
  if(nbase2 > 1){sort_atm_list(cdata_base,nbase2,ind_bas,wild_bas,
                                    natm_srch);}

/*    C) See how many tame (not wild) entries are in the base. */
  ntame = 0;
  for(i=1;i<=nbase2;i++){
    if((wild_bas[ind_bas[i]].good > 0)){break;}
    ntame = i;
  }/*endfor*/
/*    D) Check for multiple entries. */
  check_mult_base(cdata_base,ind_bas,nbase2,ntame,natm_srch,filename,
                  wild_bas);

/*==========================================================================*/
/* III) match 'em... */
  match_data_base(cdata,cdata_base,n,nbase,ntame,ind_fwd,ind_bas,ifound,igood,
                  isearch,nsearch,wild_bas,natm_srch);

  if(ntame < nbase2){
   match_wild_base(cdata,cdata_base,n,nbase,nbase2,ntame,ind_bas,ifound,
                  igood,isearch,nsearch,wild_bas,natm_srch); 
  }/*endif*/

/*==========================================================================*/
/* IV Free the memory.                                                      */
  cfree(&wild[1]);
  cfree(&wild_bas[1]);
  cfree(&wild_typ[1]);
  cfree(&ind_fwd[1]);
  cfree(&ind_bas[1]);
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  match_data_base:                                                        */
/*    match sorted atmlst lists to sorted base                              */
/*==========================================================================*/
void match_data_base(CATM_LAB *cdata, CATM_LAB *cdata_base, int n, int nbase,
                     int nbase2, int *ind,int *ind_bas,int *ifound,int *igood,
                     int *isearch,int nsearch,WILD *wild_bas,
                     int natm_srch)
{ /*begin routine */
/*==========================================================================*/
/* Local Variables */
  int iii,icount=0;
  int ityp,ibas;
  int comp;

/*==========================================================================*/
/* look for matches*/
  ibas = 1;
  for(ityp=1;ityp<=n;ityp++){
    if(ifound[ind[ityp]]==0 || igood[ind[ityp]] > 0){
      comp = 1;
      while((comp > 0) && (ibas <= nbase2)){
          match_atmlst(&comp,cdata,ind[ityp],cdata_base,
                         ind_bas[ibas],wild_bas,natm_srch);
        if(comp==0 && wild_bas[ind_bas[ibas]].good < igood[ind[ityp]]){
          ifound[ind[ityp]] = (ind_bas[ibas] <= nbase ?
                               ind_bas[ibas]:ind_bas[ibas]-nbase);
          igood[ind[ityp]] = 0;
          isearch[ind[ityp]] = nsearch;
        }/*endif*/
#ifdef DEBUG
          printf(
                "atom1: %s, atom 2: %s, atom 3: %s, atom 4: %s,label: %s \n",
                 cdata_base[ind_bas[ibas]].atm1,
                 cdata_base[ind_bas[ibas]].atm2,
                 cdata_base[ind_bas[ibas]].atm3,
                 cdata_base[ind_bas[ibas]].atm4,
                 cdata_base[ind_bas[ibas]].label);
          printf(
   "atom1: %s, atom 2: %s, atom 3: %s, atom 4: %s,label: %s %d \n",
                 cdata[ind[ityp]].atm1,
                 cdata[ind[ityp]].atm2,
                 cdata[ind[ityp]].atm3,
                 cdata[ind[ityp]].atm4,
                 cdata[ind[ityp]].label,ifound[ind[ityp]]);
          icount++;
          if((icount % 15) == 0){
            scanf("%d",&iii);
            icount = 0;
          }/*endif*/
#endif
        ibas++;
      }/*endwhile*/
      ibas--;
      if(ibas > nbase2){break;}
    }/*endif*/
  }/*endfor*/
/*--------------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void sort_atm_list(CATM_LAB *cdata,int ntyp,int *index, WILD *wild,
                        int natm_srch)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,comp;
  int k,isub;

/*=======================================================================*/
/* I) Setup                        */

  comp = 0;
  m  = ntyp/2+1;
  ir = ntyp;
  for(i=1;i<=ntyp;i++){
    index[i]=i;
  }/*endfor*/

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      index[ir]=index[1];
      ir--;
      if(ir==1){
	index[1]=rindex;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if(j<ir){
          comp_atmlst(&comp,cdata,wild,index[j],index[j+1],natm_srch);
      }/*endif*/
      if((j<ir) && (comp < 0)) j++;
      /*    b)demote */
          comp_atmlst(&comp,cdata,wild,rindex,index[j],natm_srch);
      if(comp < 0){
	index[i]=index[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*=======================================================================*/

/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/
/* comp_atmlst: compares atom list and base entries (forward compare) */
/*=======================================================================*/
void comp_atmlst(int *comp_ret,CATM_LAB *cdata,WILD *wild,int index,
                    int index1,int natm_srch)
{/*begin routine*/
/*=======================================================================*/

  int comp;

  comp = 0;
  if((wild[index].good + wild[index1].good)==0){
    comp = strcasecmp(cdata[index].atm1,cdata[index1].atm1);
    if(comp==0){
      comp = strcasecmp(cdata[index].atm2,cdata[index1].atm2);
    }/*endif*/
    if((comp==0)  && (natm_srch > 2)){
      comp = strcasecmp(cdata[index].atm3,cdata[index1].atm3);
    }/*endif*/
    if((comp==0) && (natm_srch > 3)){
      comp = strcasecmp(cdata[index].atm4,cdata[index1].atm4);
    }/*endif*/
    if(comp==0){
      comp = strcasecmp(cdata[index].label,cdata[index1].label);
    }/*endif*/
  }else{
    comp = wild[index].good - wild[index1].good;
  }/*endif*/
  *comp_ret = comp;
/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/
/* match_atmlst: matches lists and base entries (forward compare)    */
/*=======================================================================*/
void match_atmlst(int *comp_ret,CATM_LAB *cdata,
                     int ind,CATM_LAB *cdata_base,
                     int ind_bas, WILD *wild_bas,int natm_srch)
{/*begin routine*/
/*=======================================================================*/

  int comp;

  comp = 0;
  if(wild_bas[ind_bas].atm1==0){
    comp = strcasecmp(cdata[ind].atm1,cdata_base[ind_bas].atm1);
  }/*endif*/
  if((comp==0) && (wild_bas[ind_bas].atm2==0)){
    comp = strcasecmp(cdata[ind].atm2,cdata_base[ind_bas].atm2);
  }/*endif*/
  if((comp==0) && (wild_bas[ind_bas].atm3==0) && (natm_srch > 2)){
    comp = strcasecmp(cdata[ind].atm3,cdata_base[ind_bas].atm3);
  }/*endif*/
  if((comp==0) && (wild_bas[ind_bas].atm4==0) && (natm_srch > 3)){
    comp = strcasecmp(cdata[ind].atm4,cdata_base[ind_bas].atm4);
  }/*endif*/
  if((comp==0) && (wild_bas[ind_bas].label==0)){
    comp = strcasecmp(cdata[ind].label,cdata_base[ind_bas].label);
  }/*endif*/
  *comp_ret = comp;
/*-----------------------------------------------------------------------*/
}/*end routine*/
/*==========================================================================*/

/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/
/*  atmlst_not_found: die if entry were not found                        */
/*==========================================================================*/
void atmlst_not_found(int n,CATM_LAB *cdata,int *ifound,int natm,char *type)
/*=======================================================================*/
{/*begin routine*/
/*=======================================================================*/

int not_found,ityp;

/*==========================================================================*/
/* Print out entries that were not in the base */
  not_found = 0;
  for(ityp=1;ityp<=n;ityp++){
    if(ifound[ityp]==0){
      not_found++;
      if(not_found==1){
        printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
        printf("The following %ss were not in the data base:\n",type);
      }/*endif*/
      switch(natm){
        case 2: printf("atom 1: %s, atom 2: %s, label: %s \n",
                       cdata[ityp].atm1,
                       cdata[ityp].atm2,
                       cdata[ityp].label); break;
        case 3: printf("atom 1: %s, atom 2: %s, atom 3: %s, label: %s \n",
                       cdata[ityp].atm1,
                       cdata[ityp].atm2,
                       cdata[ityp].atm3,
                       cdata[ityp].label); break;
        case 4: printf(
                "atom 1: %s, atom 2: %s, atom 3: %s, atom 4: %s, label: %s \n",
                       cdata[ityp].atm1,
                       cdata[ityp].atm2,
                       cdata[ityp].atm3,
                       cdata[ityp].atm4,
                       cdata[ityp].label); break;
      }/*endswitch*/
    }/*endif*/
  }/*endfor*/
  if(not_found>0){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(1);
  }/*endif*/
/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  check_mult_base: Check for multiple entries in base                     */
/*==========================================================================*/
void check_mult_base(CATM_LAB *cdata_base,int *ind_bas,int nbase,
                     int ntame, int natm_srch, char *filename, WILD *wild_bas)
/*==========================================================================*/
{/*begin routine*/
/*==========================================================================*/

  int i;
  int nmult;

  nmult = 0;
  switch(natm_srch){
    case 2:
      for(i=1;i<=nbase-1;i++){
        if((strcasecmp(cdata_base[ind_bas[i]].atm1,
                      cdata_base[ind_bas[i+1]].atm1)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm2,
                      cdata_base[ind_bas[i+1]].atm2)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].label,
                      cdata_base[ind_bas[i+1]].label)==0)){
          nmult++;
          if(nmult==1){
            printf("$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$\n");
            printf("Multiple entry in data base file: %s\n",filename);
          }/*endif*/
          printf("atom1: %s, atom 2: %s, label: %s \n",
                 cdata_base[ind_bas[i]].atm1,
                 cdata_base[ind_bas[i]].atm2,
                 cdata_base[ind_bas[i]].label);
        }/*endif*/
      }/*endfor*/
      break;
    case 3:
      for(i=1;i<=nbase-1;i++){
        if((strcasecmp(cdata_base[ind_bas[i]].atm1,
                      cdata_base[ind_bas[i+1]].atm1)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm2,
                      cdata_base[ind_bas[i+1]].atm2)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm3,
                      cdata_base[ind_bas[i+1]].atm3)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].label,
                      cdata_base[ind_bas[i+1]].label)==0)){
          nmult++;
          if(nmult==1){
            printf("$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$\n");
            printf("Multiple entry in data base file: %s\n",filename);
          }/*endif*/
          printf("atom1: %s, atom 2: %s, atom 3: %s, label: %s \n",
                 cdata_base[ind_bas[i]].atm1,
                 cdata_base[ind_bas[i]].atm2,
                 cdata_base[ind_bas[i]].atm3,
                 cdata_base[ind_bas[i]].label);
        }/*endif*/
      }/*endfor*/
      break;
    case 4:
      for(i=1;i<=nbase-1;i++){
        if((strcasecmp(cdata_base[ind_bas[i]].atm1,
                      cdata_base[ind_bas[i+1]].atm1)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm2,
                      cdata_base[ind_bas[i+1]].atm2)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm3,
                      cdata_base[ind_bas[i+1]].atm3)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].atm4,
                      cdata_base[ind_bas[i+1]].atm4)==0)
        && (strcasecmp(cdata_base[ind_bas[i]].label,
                      cdata_base[ind_bas[i+1]].label)==0)){
          nmult++;
          if(nmult==1){
            printf("$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$\n");
            printf("Multiple entry in data base file: %s\n",filename);
          }/*endif*/
          printf(
                "atom1: %s, atom 2: %s, atom 3: %s, atom 4: %s,label: %s \n",
                 cdata_base[ind_bas[i]].atm1,
                 cdata_base[ind_bas[i]].atm2,
                 cdata_base[ind_bas[i]].atm3,
                 cdata_base[ind_bas[i]].atm4,
                 cdata_base[ind_bas[i]].label);
        }/*endif*/
      }/*endfor*/
      break;
    }/*endswitch*/
  if(nmult>0){
    printf("$$$$$$$$$$$$$$$$$$$_WARNING_$$$$$$$$$$$$$$$$$$$\n");
    fflush(stdout);
  }/*endif*/

  for(i=ntame+1;i<=nbase-1;i++){
    if(wild_bas[ind_bas[i]].good > wild_bas[ind_bas[i+1]].good ){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Internal error in search base: goodness not sorted \n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1);
    }/*endif*/
  }/*endfor*/

/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*  match_wld_base:                                                         */
/*    match sorted type lists  to sorted base                               */
/*==========================================================================*/
void match_wild_base(CATM_LAB *cdata, CATM_LAB *cdata_base, int n, int nbase,
                     int nbase2, int ntame, int *ind_bas,int *ifound,
                     int *igood,int *isearch,int nsearch,WILD *wild_bas,
                     int natm_srch)
/*==========================================================================*/
/*   Begin Routine */
    { /*begin routine */
/*==========================================================================*/
/* Local Variables */
  int iii;
  int ityp,ibas;
  int comp;
/*==========================================================================*/
/* Bubble sort the wild things */

  for(ityp=1;ityp<=n;ityp++){
     if(ifound[ityp]==0 || igood[ityp] > 0){
        for(ibas=ntame+1;ibas<=nbase2;ibas++){
          match_atmlst(&comp,cdata,ityp,cdata_base,
                         ind_bas[ibas],wild_bas,natm_srch);
          if(comp==0 && wild_bas[ind_bas[ibas]].good < igood[ityp]){
            ifound[ityp] = (ind_bas[ibas] <= nbase ?
                               ind_bas[ibas]:ind_bas[ibas]-nbase);
            igood[ityp] = 0;
            isearch[ityp] = nsearch;
            break;
          }/*endif*/  
          if(wild_bas[ind_bas[ibas]].good >= igood[ityp]){break;}
        }/*endfor*/  
     }/*endif*/
  }/*endfor*/

/*-----------------------------------------------------------------------*/
}  /*end routine*/
/*==========================================================================*/
