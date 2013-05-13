/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: sort_intra_lists                             */
/*                                                                          */
/* This subprogram blocks the intra lists to promote vectorization          */
/* reduce memory conflicts and improve shake/rattle                         */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/



#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void block_intra_lists(ATOMMAPS *atommaps,BONDED *bonded,
                       int block_std_on,int block_con_on, int nblock_min,
                       int nlen)

/*========================================================================*/
/* Explanation:                                                            */
/*  This routines calls a set of routines that sort an                     */
/*  interaction list into classes and then groups the classes into blocks  */
/*  in memory by rearranging the list.  This is done such that             */
/*  interactions of the same class are not in the same memory              */
/*  block. Interactions of the neighboring classes (1-2, 2-3, ...)         */
/*  are also excluded from occupying the same memory block.                */
/*  The classes are chosen to reduce bank conflicts by separating          */
/*  references to atoms in memory.                                         */
/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int nblock_tot,nconflict_tot,i;
  int nnow,iii,*iconf_now;
  double perconfs,block_size_avg,navg;

/*========================================================================*/
/* 0) Output */

  PRINT_LINE_STAR;
  printf("Blocking the intra molecular lists\n");
  PRINT_LINE_DASH;printf("\n");
   
  nblock_tot        = 0;
  nconflict_tot     = 0;
  block_size_avg    = 0.0;
  navg              = 0.0;

/*========================================================================*/
/* I) Pow Bonds */

  printf("Blocking the pow bonds\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->bond.npow>0)){
    bonded->bond.block_pow_on = 1;
    two_vector_block(bonded->bond.npow,&(bonded->bond.j1_pow),
                    &(bonded->bond.j2_pow),&(bonded->bond.jtyp_pow),
                    &(bonded->bond.nblock_pow),&(bonded->bond.nblock_size_pow),
                    &(bonded->bond.iblock_pow_size),
                    &(bonded->bond.iblock_pow_conflict_1),
                    &(bonded->bond.iblock_pow_conflict_2),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->bond.nblock_size_pow;
     nnow = bonded->bond.nblock_pow;
     nblock_tot += 2*nnow;
     iconf_now = bonded->bond.iblock_pow_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bond.iblock_pow_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* II) Con Bonds */

  printf("Blocking the con bonds\n");fflush(stdout);
  if((block_con_on==1)&&(bonded->bond.ncon>0)){
    bonded->bond.block_con_on = 1;
    two_vector_block(bonded->bond.ncon,&(bonded->bond.j1_con),
                    &(bonded->bond.j2_con),&(bonded->bond.jtyp_con),
                    &(bonded->bond.nblock_con),&(bonded->bond.nblock_size_con),
                    &(bonded->bond.iblock_con_size),
                    &(bonded->bond.iblock_con_conflict_1),
                    &(bonded->bond.iblock_con_conflict_2),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->bond.nblock_size_con;
     nnow = bonded->bond.nblock_con;
     nblock_tot += 2*nnow;
     iconf_now = bonded->bond.iblock_con_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bond.iblock_con_conflict_2;
  }/*endif*/

/*========================================================================*/
/* III) Pow Bends */

  printf("Blocking the pow bends\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->bend.npow>0)){
    bonded->bend.block_pow_on = 1;
    three_vector_block(bonded->bend.npow,&(bonded->bend.j1_pow),
                    &(bonded->bend.j2_pow),&(bonded->bend.j3_pow),
                    &(bonded->bend.jtyp_pow),
                    &(bonded->bend.nblock_pow),&(bonded->bend.nblock_size_pow),
                    &(bonded->bend.iblock_pow_size),
                    &(bonded->bend.iblock_pow_conflict_1),
                    &(bonded->bend.iblock_pow_conflict_2),
                    &(bonded->bend.iblock_pow_conflict_3),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->bend.nblock_size_pow;
     nnow = bonded->bend.nblock_pow;
     nblock_tot += 3*nnow;
     iconf_now = bonded->bend.iblock_pow_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend.iblock_pow_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend.iblock_pow_conflict_3;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* IV) Con Bends */

  printf("Blocking the con bends\n");fflush(stdout);
  if((block_con_on==1)&&(bonded->bend.ncon>0)){
    bonded->bend.block_con_on = 1;
    three_vector_block(bonded->bend.ncon,&(bonded->bend.j1_con),
                    &(bonded->bend.j2_con),&(bonded->bend.j3_con),
                    &(bonded->bend.jtyp_con),
                    &(bonded->bend.nblock_con),&(bonded->bend.nblock_size_con),
                    &(bonded->bend.iblock_con_size),
                    &(bonded->bend.iblock_con_conflict_1),
                    &(bonded->bend.iblock_con_conflict_2),
                    &(bonded->bend.iblock_con_conflict_3),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->bend.nblock_size_con;
     nnow = bonded->bend.nblock_con;
     nblock_tot += 3*nnow;
     iconf_now = bonded->bend.iblock_con_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend.iblock_con_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend.iblock_con_conflict_3;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* V) Pow Tors */

  printf("Blocking the pow torsions\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->tors.npow>0)){
    bonded->tors.block_pow_on = 1;
    four_vector_block(bonded->tors.npow,&(bonded->tors.j1_pow),
                    &(bonded->tors.j2_pow),&(bonded->tors.j3_pow),
                    &(bonded->tors.j4_pow),&(bonded->tors.jtyp_pow),
                    &(bonded->tors.nblock_pow),&(bonded->tors.nblock_size_pow),
                    &(bonded->tors.iblock_pow_size),
                    &(bonded->tors.iblock_pow_conflict_1),
                    &(bonded->tors.iblock_pow_conflict_2),
                    &(bonded->tors.iblock_pow_conflict_3),
                    &(bonded->tors.iblock_pow_conflict_4),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->tors.nblock_size_pow;
     nnow = bonded->tors.nblock_pow;
     nblock_tot += 4*nnow;
     iconf_now = bonded->tors.iblock_pow_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_pow_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_pow_conflict_3;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_pow_conflict_4;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* VI) Con Tors */

  printf("Blocking the con torsions\n");fflush(stdout);
  if((block_con_on==1)&&(bonded->tors.ncon>0)){
    bonded->tors.block_con_on = 1;
    four_vector_block(bonded->tors.ncon,&(bonded->tors.j1_con),
                    &(bonded->tors.j2_con),&(bonded->tors.j3_con),
                    &(bonded->tors.j4_con),&(bonded->tors.jtyp_con),
                    &(bonded->tors.nblock_con),&(bonded->tors.nblock_size_con),
                    &(bonded->tors.iblock_con_size),
                    &(bonded->tors.iblock_con_conflict_1),
                    &(bonded->tors.iblock_con_conflict_2),
                    &(bonded->tors.iblock_con_conflict_3),
                    &(bonded->tors.iblock_con_conflict_4),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->tors.nblock_size_con;
     nnow = bonded->tors.nblock_con;
     nblock_tot += 4*nnow;
     iconf_now = bonded->tors.iblock_con_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_con_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_con_conflict_3;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->tors.iblock_con_conflict_4;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* VII) Uri Bradleys */

  printf("Blocking the Uri-Bradley bends\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->bend_bnd.num>0)){
    bonded->bend_bnd.block_on = 1;
    three_vector_block(bonded->bend_bnd.num,&(bonded->bend_bnd.j1),
                    &(bonded->bend_bnd.j2),&(bonded->bend_bnd.j3),
                    &(bonded->bend_bnd.jtyp),
                    &(bonded->bend_bnd.nblock),&(bonded->bend_bnd.nblock_size),
                    &(bonded->bend_bnd.iblock_size),
                    &(bonded->bend_bnd.iblock_conflict_1),
                    &(bonded->bend_bnd.iblock_conflict_2),
                    &(bonded->bend_bnd.iblock_conflict_3),atommaps,
                    nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->bend_bnd.nblock_size;
     nnow = bonded->bend_bnd.nblock;
     nblock_tot += 3*nnow;
     iconf_now = bonded->bend_bnd.iblock_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend_bnd.iblock_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->bend_bnd.iblock_conflict_3;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* VIII) Onefours */

  printf("Blocking the onefours\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->onfo.num>0)){
   bonded->onfo.block_on = 1;
   two_vector_block(bonded->onfo.num,&(bonded->onfo.j1),
                   &(bonded->onfo.j2),&(bonded->onfo.jtyp),
                   &(bonded->onfo.nblock),&(bonded->onfo.nblock_size),
                   &(bonded->onfo.iblock_size),
                   &(bonded->onfo.iblock_conflict_1),
                   &(bonded->onfo.iblock_conflict_2),atommaps,
                   nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->onfo.nblock_size;
     nnow = bonded->onfo.nblock;
     nblock_tot += 2*nnow;
     iconf_now = bonded->onfo.iblock_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->onfo.iblock_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* IX) Ecorr */

  printf("Blocking the ecorrs\n");fflush(stdout);
  if((block_std_on==1)&&(bonded->ecor.num>0)){
   bonded->ecor.block_on = 1;
   two_vector_block_notyp(bonded->ecor.num,&(bonded->ecor.j1),
                        &(bonded->ecor.j2),
                        &(bonded->ecor.nblock),&(bonded->ecor.nblock_size),
                        &(bonded->ecor.iblock_size),
                        &(bonded->ecor.iblock_conflict_1),
                        &(bonded->ecor.iblock_conflict_2),atommaps,
                        nblock_min,nlen);

     navg+=1.0;
     block_size_avg += bonded->ecor.nblock_size;
     nnow = bonded->ecor.nblock;
     nblock_tot += 2*nnow;
     iconf_now = bonded->ecor.iblock_conflict_1;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
     iconf_now = bonded->ecor.iblock_conflict_2;
     for(i=1;i<=nnow;i++){nconflict_tot+=iconf_now[i];}
  }/*endif*/

/*========================================================================*/
/* 0) Output */

  block_size_avg/=MAX(navg,1.0);
  nblock_tot = MAX(nblock_tot,1);
  printf("The average block size is %g\n",block_size_avg);
  perconfs = 100.0*( ((double)nconflict_tot)/((double)nblock_tot) );
  printf("The percent conflicts is %g\n",perconfs);

  printf("\n");
  PRINT_LINE_DASH;
  printf("Completed intra molecular blocking procedure \n");
  PRINT_LINE_STAR;printf("\n");

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/

 



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void two_vector_block_notyp(int num,int **j1_now,int **j2_now,
                      int *nblock_now,int *nblock_size_now,
                      int **iblock_size_now,
                      int **iblock_conflict_1_now,int **iblock_conflict_2_now,
                      ATOMMAPS *atommaps,int nblock_size_min,
                      int nblock_size_max)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int i,j,lower,upper,nmem;
 int iblock_now,iblock_start,index;
 int nnow,nblock,nblock_size;
 int *iblock_size,*iblock_conflict_1,*iblock_conflict_2;
 int *j1,*j2;
 int nclass,iii;
 int *class_lst,*class_lst_num,*class_lst_off,*class_lst_map;
 int *j1_tmp,*j2_tmp;

/*==========================================================================*/
/* I) Malloc some nice memory and assign local pointers */

    class_lst     = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_num = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_off = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_map = (int *)cmalloc(num*sizeof(int))-1;
    j1_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j2_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j1            = *j1_now;
    j2            = *j2_now;

/*==========================================================================*/
/* II ) Sort the interactions into classes. The number of classes equals    */
/*      the number of blocks                                                */

    for(i=1;i<=num;i++){
      class_lst[i]     = j1[i];
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
    }/*endfor*/

    make_class_lst(num,&nclass,&nblock,class_lst,class_lst_num,class_lst_off,
                   class_lst_map,atommaps);

/*==========================================================================*/
/* IV ) Determine the maximum block size                                   */

    *iblock_size_now     = (int *)cmalloc(nblock*sizeof(int))-1;
    iblock_size       = *iblock_size_now;

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

    nblock_size = 0;
    for(i=1;i<=nblock;i++){nblock_size = MAX(iblock_size[i],nblock_size);}
    nblock_size = MIN(nblock_size,nblock_size_max);

/*==========================================================================*/
/* III ) Realloc the list and reassign local pointers                       */

    (*nblock_now)                    = nblock;
    (*nblock_size_now)               = nblock_size;
    nmem                             = nblock_size*nblock;
    if((nmem % 2)==0){nmem++;}
    (*j1_now)   = (int *)crealloc(&(*j1_now)[1],nmem*sizeof(int))-1;
    (*j2_now)   = (int *)crealloc(&(*j2_now)[1],nmem*sizeof(int))-1;
    (j1_tmp)   = (int *)crealloc(&(j1_tmp)[1],nmem*sizeof(int))-1;
    (j2_tmp)   = (int *)crealloc(&(j2_tmp)[1],nmem*sizeof(int))-1;
    *iblock_conflict_1_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_2_now = (int *)cmalloc(nblock*sizeof(int))-1;
    j1    = *j1_now;
    j2    = *j2_now;

/*==========================================================================*/
/* IV ) Block the list according to the classes                             */

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
       index       = (iblock_now-1)*nblock_size + iblock_size[iblock_now];
       j1[index]   = j1_tmp[class_lst_map[j]];
       j2[index]   = j2_tmp[class_lst_map[j]];
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

/* iv) Check for conflicts within the blocks */
    for(i=1;i<=nmem;i++){
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
    }/*endfor*/

    iblock_conflict_1   = *iblock_conflict_1_now;
    iblock_conflict_2   = *iblock_conflict_2_now;
    check_intra_conflict(nblock,nblock_size,j1_tmp,
                         iblock_size,iblock_conflict_1);
    check_intra_conflict(nblock,nblock_size,j2_tmp,
                         iblock_size,iblock_conflict_2);

/*========================================================================*/
/* IX) Free */

    free(&class_lst[1]);
    free(&class_lst_num[1]);
    free(&class_lst_off[1]);
    free(&class_lst_map[1]);
    free(&j1_tmp[1]);       
    free(&j2_tmp[1]);       

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void two_vector_block(int num,int **j1_now,int **j2_now,int **jtyp_now,
                      int *nblock_now,int *nblock_size_now,
                      int **iblock_size_now,
                      int **iblock_conflict_1_now,int **iblock_conflict_2_now,
                      ATOMMAPS *atommaps,int nblock_size_min,
                      int nblock_size_max)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int i,j,lower,upper,nmem,iii;
 int iblock_now,iblock_start,index;
 int nnow,nblock,nblock_size;
 int *iblock_size,*iblock_conflict_1,*iblock_conflict_2;
 int *j1,*j2,*jtyp;
 int nclass;
 int *class_lst,*class_lst_num,*class_lst_off,*class_lst_map;
 int *j1_tmp,*j2_tmp,*jtyp_tmp;

/*==========================================================================*/
/* I) Malloc some nice memory and assign local pointers */

    class_lst     = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_num = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_off = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_map = (int *)cmalloc(num*sizeof(int))-1;
    j1_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j2_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    jtyp_tmp      = (int *)cmalloc(num*sizeof(int))-1;
    j1            = *j1_now;
    j2            = *j2_now;
    jtyp          = *jtyp_now;

/*==========================================================================*/
/* II ) Sort the interactions into classes. The number of classes equals    */
/*      the number of blocks                                                */

    for(i=1;i<=num;i++){
      class_lst[i]     = j1[i];
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
      jtyp_tmp[i]      = jtyp[i];
    }/*endfor*/

    printf("Making the class list\n");fflush(stdout);
    make_class_lst(num,&nclass,&nblock,class_lst,class_lst_num,class_lst_off,
                   class_lst_map,atommaps);

/*==========================================================================*/
/* IV ) Determine the maximum block size                                   */

    *iblock_size_now     = (int *)cmalloc(nblock*sizeof(int))-1;
    iblock_size       = *iblock_size_now;

    printf("Determing the block size\n");fflush(stdout);
    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

    nblock_size = 0;
    for(i=1;i<=nblock;i++){
     nblock_size = MAX(iblock_size[i],nblock_size);
    }/*endfor*/
    nblock_size = MIN(nblock_size,nblock_size_max);

/*==========================================================================*/
/* III ) Realloc the list and reassign local pointers                       */

    (*nblock_now)                    = nblock;
    (*nblock_size_now)               = nblock_size;
    nmem                             = nblock_size*nblock;
    if((nmem % 2)==0){nmem++;}
    (*j1_now)   = (int *)crealloc(&(*j1_now)[1],nmem*sizeof(int))-1;
    (*j2_now)   = (int *)crealloc(&(*j2_now)[1],nmem*sizeof(int))-1;
    (*jtyp_now) = (int *)crealloc(&(*jtyp_now)[1],nmem*sizeof(int))-1;
    (j1_tmp)   = (int *)crealloc(&(j1_tmp)[1],nmem*sizeof(int))-1;
    (j2_tmp)   = (int *)crealloc(&(j2_tmp)[1],nmem*sizeof(int))-1;
    (jtyp_tmp) = (int *)crealloc(&(jtyp_tmp)[1],nmem*sizeof(int))-1;
    *iblock_conflict_1_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_2_now = (int *)cmalloc(nblock*sizeof(int))-1;
    j1    = *j1_now;
    j2    = *j2_now;
    jtyp  = *jtyp_now;

/*==========================================================================*/
/* IV ) Block the list according to the classes                             */

    printf("Filling the list\n");fflush(stdout);
    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
       index       = (iblock_now-1)*nblock_size + iblock_size[iblock_now];
       j1[index]   = j1_tmp[class_lst_map[j]];
       j2[index]   = j2_tmp[class_lst_map[j]];
       jtyp[index] = jtyp_tmp[class_lst_map[j]];
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

/* iv) Check for conflicts within the blocks */
    for(i=1;i<=nmem;i++){
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
    }/*endfor*/

    printf("Doing up the conflcits\n");fflush(stdout);
    iblock_conflict_1   = *iblock_conflict_1_now;
    iblock_conflict_2   = *iblock_conflict_2_now;
    check_intra_conflict(nblock,nblock_size,j1_tmp,
                         iblock_size,iblock_conflict_1);
    check_intra_conflict(nblock,nblock_size,j2_tmp,
                         iblock_size,iblock_conflict_2);

/*========================================================================*/
/* IX) Free */

    free(&class_lst[1]);
    free(&class_lst_num[1]);
    free(&class_lst_off[1]);
    free(&class_lst_map[1]);
    free(&j1_tmp[1]);       
    free(&j2_tmp[1]);       
    free(&jtyp_tmp[1]);     

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void three_vector_block(int num,int **j1_now,int **j2_now,
                      int **j3_now,
                      int **jtyp_now,
                      int *nblock_now,int *nblock_size_now,
                      int **iblock_size_now,
                      int **iblock_conflict_1_now,
                      int **iblock_conflict_2_now,
                      int **iblock_conflict_3_now,
                      ATOMMAPS *atommaps,int nblock_size_min,
                      int nblock_size_max)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int i,j,lower,upper,nmem,ntot;
 int iblock_now,iblock_start,index;
 int nnow,nblock,nblock_size;
 int *iblock_size,*iblock_conflict_1,*iblock_conflict_2;
 int *iblock_conflict_3;
 int *j1,*j2,*j3,*jtyp;
 int nclass;
 int *class_lst,*class_lst_num,*class_lst_off,*class_lst_map;
 int *j1_tmp,*j2_tmp,*j3_tmp,*jtyp_tmp;

/*==========================================================================*/
/* I) Malloc some nice memory and assign local pointers */

    class_lst     = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_num = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_off = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_map = (int *)cmalloc(num*sizeof(int))-1;
    j1_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j2_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j3_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    jtyp_tmp      = (int *)cmalloc(num*sizeof(int))-1;
    j1            = *j1_now;
    j2            = *j2_now;
    j3            = *j3_now;
    jtyp          = *jtyp_now;

/*==========================================================================*/
/* II ) Sort the interactions into classes. The number of classes equals    */
/*      the number of blocks                                                */

    for(i=1;i<=num;i++){
      class_lst[i]     = j1[i];
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
      j3_tmp[i]        = j3[i];
      jtyp_tmp[i]      = jtyp[i];
    }/*endfor*/

    make_class_lst(num,&nclass,&nblock,class_lst,class_lst_num,class_lst_off,
                   class_lst_map,atommaps);

/*==========================================================================*/
/* IV ) Determine the maximum block size                                   */

    *iblock_size_now     = (int *)cmalloc(nblock*sizeof(int))-1;
    iblock_size       = *iblock_size_now;

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

    nblock_size = 0;
    for(i=1;i<=nblock;i++){
      nblock_size = MAX(iblock_size[i],nblock_size);
    }
    nblock_size = MIN(nblock_size,nblock_size_max);

/*==========================================================================*/
/* III ) Realloc the list and reassign local pointers                       */

    (*nblock_now)                    = nblock;
    (*nblock_size_now)               = nblock_size;
    nmem                             = nblock_size*nblock;
    if((nmem % 2)==0){nmem++;}
    (*j1_now)   = (int *)crealloc(&(*j1_now)[1],nmem*sizeof(int))-1;
    (*j2_now)   = (int *)crealloc(&(*j2_now)[1],nmem*sizeof(int))-1;
    (*j3_now)   = (int *)crealloc(&(*j3_now)[1],nmem*sizeof(int))-1;
    (*jtyp_now) = (int *)crealloc(&(*jtyp_now)[1],nmem*sizeof(int))-1;
    (j1_tmp)   = (int *)crealloc(&(j1_tmp)[1],nmem*sizeof(int))-1;
    (j2_tmp)   = (int *)crealloc(&(j2_tmp)[1],nmem*sizeof(int))-1;
    (j3_tmp)   = (int *)crealloc(&(j3_tmp)[1],nmem*sizeof(int))-1;
    (jtyp_tmp) = (int *)crealloc(&(jtyp_tmp)[1],nmem*sizeof(int))-1;
    *iblock_conflict_1_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_2_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_3_now = (int *)cmalloc(nblock*sizeof(int))-1;
    j1    = *j1_now;
    j2    = *j2_now;
    j3    = *j3_now;
    jtyp  = *jtyp_now;

/*==========================================================================*/
/* IV ) Block the list according to the classes                             */

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
       index       = (iblock_now-1)*nblock_size + iblock_size[iblock_now];
       j1[index]   = j1_tmp[class_lst_map[j]];
       j2[index]   = j2_tmp[class_lst_map[j]];
       j3[index]   = j3_tmp[class_lst_map[j]];
       jtyp[index] = jtyp_tmp[class_lst_map[j]];
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

/* iv) Check for conflicts within the blocks */
    for(i=1;i<=nmem;i++){
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
      j3_tmp[i]        = j3[i];
    }/*endfor*/
    iblock_conflict_1   = *iblock_conflict_1_now;
    iblock_conflict_2   = *iblock_conflict_2_now;
    iblock_conflict_3   = *iblock_conflict_3_now;
    check_intra_conflict(nblock,nblock_size,j1_tmp,
                         iblock_size,iblock_conflict_1);
    check_intra_conflict(nblock,nblock_size,j2_tmp,
                         iblock_size,iblock_conflict_2);
    check_intra_conflict(nblock,nblock_size,j3_tmp,
                         iblock_size,iblock_conflict_3);

/*========================================================================*/
/* IX) Free */

    free(&class_lst[1]);
    free(&class_lst_num[1]);
    free(&class_lst_off[1]);
    free(&class_lst_map[1]);
    free(&j1_tmp[1]);       
    free(&j2_tmp[1]);       
    free(&j3_tmp[1]);       
    free(&jtyp_tmp[1]);     

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void four_vector_block(int num,int **j1_now,int **j2_now,int **j3_now,
                      int **j4_now,int **jtyp_now,
                      int *nblock_now,int *nblock_size_now,
                      int **iblock_size_now,
                      int **iblock_conflict_1_now,
                      int **iblock_conflict_2_now,
                      int **iblock_conflict_3_now,
                      int **iblock_conflict_4_now,
                      ATOMMAPS *atommaps,int nblock_size_min,
                      int nblock_size_max)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

 int i,j,lower,upper,nmem;
 int iblock_now,iblock_start,index;
 int nnow,nblock,nblock_size;
 int *iblock_size,*iblock_conflict_1,*iblock_conflict_2;
 int *iblock_conflict_3,*iblock_conflict_4;
 int *j1,*j2,*j3,*j4,*jtyp;
 int nclass;
 int *class_lst,*class_lst_num,*class_lst_off,*class_lst_map;
 int *j1_tmp,*j2_tmp,*j3_tmp,*j4_tmp,*jtyp_tmp;

/*==========================================================================*/
/* I) Malloc some nice memory and assign local pointers */

    class_lst     = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_num = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_off = (int *)cmalloc(num*sizeof(int))-1;
    class_lst_map = (int *)cmalloc(num*sizeof(int))-1;
    j1_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j2_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j3_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    j4_tmp        = (int *)cmalloc(num*sizeof(int))-1;
    jtyp_tmp      = (int *)cmalloc(num*sizeof(int))-1;
    j1            = *j1_now;
    j2            = *j2_now;
    j3            = *j3_now;
    j4            = *j4_now;
    jtyp          = *jtyp_now;

/*==========================================================================*/
/* II ) Sort the interactions into classes. The number of classes equals    */
/*      the number of blocks                                                */

    for(i=1;i<=num;i++){
      class_lst[i]     = j1[i];
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
      j3_tmp[i]        = j3[i];
      j4_tmp[i]        = j4[i];
      jtyp_tmp[i]      = jtyp[i];
    }/*endfor*/

    make_class_lst(num,&nclass,&nblock,class_lst,class_lst_num,class_lst_off,
                   class_lst_map,atommaps);

/*==========================================================================*/
/* IV ) Determine the maximum block size                                   */

    *iblock_size_now     = (int *)cmalloc(nblock*sizeof(int))-1;
    iblock_size       = *iblock_size_now;

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

    nblock_size = 0;
    for(i=1;i<=nblock;i++){nblock_size = MAX(iblock_size[i],nblock_size);}
    nblock_size = MIN(nblock_size,nblock_size_max);

/*==========================================================================*/
/* III ) Realloc the list and reassign local pointers                       */

    (*nblock_now)                    = nblock;
    (*nblock_size_now)               = nblock_size;
    nmem                             = nblock_size*nblock;
    if((nmem % 2)==0){nmem++;}
    (*j1_now)   = (int *)crealloc(&(*j1_now)[1],nmem*sizeof(int))-1;
    (*j2_now)   = (int *)crealloc(&(*j2_now)[1],nmem*sizeof(int))-1;
    (*j3_now)   = (int *)crealloc(&(*j3_now)[1],nmem*sizeof(int))-1;
    (*j4_now)   = (int *)crealloc(&(*j4_now)[1],nmem*sizeof(int))-1;
    (*jtyp_now) = (int *)crealloc(&(*jtyp_now)[1],nmem*sizeof(int))-1;
    (j1_tmp)   = (int *)crealloc(&(j1_tmp)[1],nmem*sizeof(int))-1;
    (j2_tmp)   = (int *)crealloc(&(j2_tmp)[1],nmem*sizeof(int))-1;
    (j3_tmp)   = (int *)crealloc(&(j3_tmp)[1],nmem*sizeof(int))-1;
    (j4_tmp)   = (int *)crealloc(&(j4_tmp)[1],nmem*sizeof(int))-1;
    (jtyp_tmp) = (int *)crealloc(&(jtyp_tmp)[1],nmem*sizeof(int))-1;
    *iblock_conflict_1_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_2_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_3_now = (int *)cmalloc(nblock*sizeof(int))-1;
    *iblock_conflict_4_now = (int *)cmalloc(nblock*sizeof(int))-1;
    j1    = *j1_now;
    j2    = *j2_now;
    j3    = *j3_now;
    j4    = *j4_now;
    jtyp  = *jtyp_now;

/*==========================================================================*/
/* IV ) Block the list according to the classes                             */

    for(i=1;i<=nblock;i++){iblock_size[i]=0;}
    iblock_now = 1;
    for(i=1;i<=nclass;i++){
/*  i) put the members of this class into the odd/even blocks */
      iblock_start = iblock_now;
      lower = class_lst_off[i]+1;
      upper = class_lst_off[i]+class_lst_num[i];
      for(j=lower;j<=upper;j++){
       iblock_size[iblock_now]++;
       index       = (iblock_now-1)*nblock_size + iblock_size[iblock_now];
       j1[index]   = j1_tmp[class_lst_map[j]];
       j2[index]   = j2_tmp[class_lst_map[j]];
       j3[index]   = j3_tmp[class_lst_map[j]];
       j4[index]   = j4_tmp[class_lst_map[j]];
       jtyp[index] = jtyp_tmp[class_lst_map[j]];
/*  ii) block wrap conditions */
       iblock_now+=2;
       if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
       if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       if(iblock_now==iblock_start){
         iblock_now++;
         if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
         if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
       }/*endif*/
      }/*endfor*/
/*  iii) change from odd/even filling to even/odd block filling */
      iblock_now++;
      if((iblock_now>nblock)&&((iblock_now % 2)==1)){iblock_now=1;}
      if((iblock_now>nblock)&&((iblock_now % 2)==0)){iblock_now=2;}
    }/*endfor*/

/* iv) Check for conflicts within the blocks */
    for(i=1;i<=nmem;i++){
      j1_tmp[i]        = j1[i];
      j2_tmp[i]        = j2[i];
      j3_tmp[i]        = j3[i];
      j4_tmp[i]        = j4[i];
    }/*endfor*/
    iblock_conflict_1   = *iblock_conflict_1_now;
    iblock_conflict_2   = *iblock_conflict_2_now;
    iblock_conflict_3   = *iblock_conflict_3_now;
    iblock_conflict_4   = *iblock_conflict_4_now;
    check_intra_conflict(nblock,nblock_size,j1_tmp,
                         iblock_size,iblock_conflict_1);
    check_intra_conflict(nblock,nblock_size,j2_tmp,
                         iblock_size,iblock_conflict_2);
    check_intra_conflict(nblock,nblock_size,j3_tmp,
                         iblock_size,iblock_conflict_3);
    check_intra_conflict(nblock,nblock_size,j4_tmp,
                         iblock_size,iblock_conflict_4);

/*========================================================================*/
/* IX) Free */

    free(&class_lst[1]);
    free(&class_lst_num[1]);
    free(&class_lst_off[1]);
    free(&class_lst_map[1]);
    free(&j1_tmp[1]);       
    free(&j2_tmp[1]);       
    free(&j3_tmp[1]);       
    free(&j4_tmp[1]);       
    free(&jtyp_tmp[1]);     

/*========================================================================*/
} /*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void make_class_lst(int num,int *nclass,int *nblock,int *class_lst,
                    int *class_lst_num,
                    int *class_lst_off,int *class_lst_map,ATOMMAPS *atommaps)

/*========================================================================*/
/*             Begin subprogram:                                          */
{/*begin routine*/
/*========================================================================*/
/*             Local variable declarations                                */

  int nclass_loc,nblock_loc,i,i1;
  int  imol_typ, imol_num, ires_num,ires_eff ;
  int  imol_typ_now, imol_num_now, ires_num_now;
  int *iatm_mol_typ = atommaps->iatm_mol_typ;
  int *iatm_mol_num = atommaps->iatm_mol_num;
  int *iatm_res_num = atommaps->iatm_res_num;

/*=======================================================================*/
/* I) Sort the class list according to atom index keeping track of the map */

   if(num>1){
    vec_sort_map(num,class_lst,class_lst_map);
   }/*endif*/

/*======================================================================= */
/* II) The atoms are laid out in memory by molecule_typ,molecule num,     */
/*     residue num. We just need to find the bondaries between classes    */
/*     In order to minimize class size let 10 be the max size before a   */
/*    break is forced*/

   nclass_loc = 1;
   nblock_loc = 1;
   class_lst_num[1] = 1;
   imol_typ = iatm_mol_typ[class_lst[1]];
   imol_num = iatm_mol_num[class_lst[1]];
   ires_num = iatm_res_num[class_lst[1]];
   ires_eff = 1;
   for(i=2;i<=num;i++){
    imol_typ_now = iatm_mol_typ[class_lst[i]];
    imol_num_now = iatm_mol_num[class_lst[i]];
    ires_num_now = iatm_res_num[class_lst[i]];
    if((imol_typ_now==imol_typ)&&(imol_num_now==imol_num)&&
       (ires_num_now==ires_num)&&(class_lst_num[nclass_loc]<10)){
       class_lst_num[nclass_loc]++;
       if((imol_typ_now==imol_typ)&&(imol_num_now==imol_num)){
          ires_eff++;
       }else{
          ires_eff=1;
       }/*endif*/
    }else{
       if(ires_eff>=2){ /* avoids potential conflicts */
        nblock_loc = MAX(nblock_loc,2*class_lst_num[nclass_loc]);
       }else{
        nblock_loc = MAX(nblock_loc,class_lst_num[nclass_loc]);
       }/*endif*/
       nclass_loc++;
       class_lst_num[nclass_loc] = 1;
    }/*endif*/
    imol_typ = imol_typ_now;
    imol_num = imol_num_now;
    ires_num = ires_num_now;
   }/*endfor*/
 
   class_lst_off[1] = 0;
   for(i = 2;i<=nclass_loc;i++){
    i1 = i-1;
    class_lst_off[i] = class_lst_off[i1]
                     + class_lst_num[i1];
  }/*endfor*/

   *nblock = nblock_loc;
   *nclass = nclass_loc;

/*========================================================================*/
} /*end routine*/ 
/*========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vec_sort_map(int n, int index[],int map[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex,rmap;
  int k,*kndex,*mask,isub;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;
  for(i=1;i<=n;i++){
   map[i] = i;
  }/*endfor*/

/*=======================================================================*/
/* II) Sort array index keeping track of the map                         */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rmap   = map[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex   = index[ir];
      rmap     = map[ir];
      index[ir]=index[1];
      map[ir]  = map[1];
      ir--;
      if(ir==1){
	index[1]=rindex;
        map[1]  = rmap;
	break;
      }/*endif*/
    }/*endif*/
/*---------------------------------------------------------------------*/
/*  C)put rindex in appropriate slot */
    i=m;
    j=2*m;
    while(j<=ir){
      /*    a)compare to rindex to underling */
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
	index[i]=index[j];
        map[i]  = map[j];
	i=j;
	j=2*j;
      }else{
	/*    c)if no demotations exit while */
	j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    map[i]   = rmap;
  }/*endfor*/

/*========================================================================*/
} /*end routine*/ 
/*========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void vec_sort(int n, int index[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex;
  int k,*kndex,*mask,isub;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index keeping track of the map                         */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex   = index[ir];
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
      if((j<ir) && (index[j]< index[(j+1)])) j++;
      /*    b)demote */
      if(rindex<index[j]){
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

/*========================================================================*/
} /*end routine*/ 
/*========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_intra_conflict(int nblock,int nblock_size,int *jindex,
                         int *iblock_size,int *iblock_conflict)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int i,istart,j,n;

/*=======================================================================    */
/* Sort the vector of indices by blocks and check for repeat indicies        */
/* within each block. If a repeating index is found in a block the conflict  */
/* flag for that block is set to one                                         */

  istart = 0;
  for(i=1;i<=nblock;i++){
    n = iblock_size[i];
    printf("sort i %d n %d\n",i,n);fflush(stdout);
    if(n>2){
     vec_sort(n,&jindex[istart]);
    }/*endif*/
    iblock_conflict[i] = 0;
    printf("checking conflicts %d %d\n",i,n);fflush(stdout);
    for(j=istart+2;j<=istart+n;j++){
     if(jindex[j]==jindex[(j-1)]){iblock_conflict[i]=1;}
    }/*endfor*/
    printf("DONE\n");fflush(stdout);
    istart += nblock_size;
  }/*endfor*/

/*========================================================================*/
} /*end routine*/ 
/*========================================================================*/
