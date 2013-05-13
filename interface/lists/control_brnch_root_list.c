/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: control_sim_parms.c                          */
/*                                                                          */
/*                                                                          */
/* This subprogram reads in user simulation params and echoes them          */
/* and the default parameters to a file                                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_lists_entry.h"
#include "../proto_defs/proto_lists_local.h"
#include "../proto_defs/proto_real_space_local.h"

#define DEBUG_EXCL_EXP_OFF
#define DEBUG_ATM_TYPE_OFF
#define DEBUG_CUTOFF_CHCK_OFF
#define DEBUG_BRNCH_ROOT_OFF
#define DEBUG_ADD_LIST_OFF
#define DEBUG_SORT_ADD_LIST_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/* Control_brnch_root_list:                                                 */
/* ------------------------                                                 */
/*    I) Finds all valence=1 atoms and assigns them to be ``branches'' of   */
/*       the ``root'' atoms to whom they are bonded if eq_bond < cut_bond.  */
/*   II) Assuming the root atoms are used to generate the pair interaction  */
/*       list, some interactions, will be left out due to root-root         */
/*       exclusion rules.  An ``addition'' list is generated for each       */
/*       atom which contains these missing interactions.                    */
/*  III) The valence = 1 schemes guarantees that there is no need to        */
/*       eliminate excluded interactions from the root-root list when       */
/*       expanding it to form the full list.                                */
/*       That is the root atoms will have the same or more excluded         */
/*       interactions thans its branch atoms barring some incredibly        */
/*       unchemical bonding scheme. In this case, the routine exits         */
/*       with an error. Also,the exclusions of a given root atom must be    */
/*       themselves root atoms or branches of root atoms which are          */
/*       themselves exclusions of the root atom in question. If this        */
/*       second condition is not met, an error is generated.                */
/*   IV) Cutoffs are checked to make sure that the root-root list when      */
/*       extended by the branches of each root and the addition list of     */
/*       each atom will contain all possible pairs. The code exits if all   */
/*       possible pairs within the specified cutoffs cannot be generated.   */
/*    V) 2 x 1 Group Cons not implemented. Ghost atoms are root atoms.      */
/*==========================================================================*/

void control_brnch_root_list(CLASS *class,BONDED *bonded)

/*=======================================================================*/
/*  Begin Routine */
  {/*begin routine*/ 
/*=======================================================================*/
/*          Local variable declarations                                  */

   int nexl_temp;
   int i,j,iii,ic,ind,ktemp;

   int natm_add_tot;
   int nbrnch_root_max,nbrnch_tot,nroot_tot;
   double brnch_root_dist_max,now_memory;
   double cutdif_min_brnch_root;
   double cutdif_min_brnch_brnch;
   double cutdif_min_brnch_root_res;
   double cutdif_min_brnch_brnch_res;
   double cutskin_bb_min;
   double r2_nocheck_dist;
   double cutskin_bb_min_res;
   double r2_nocheck_dist_res;
   int ires_flag;


   int *nvalence,*nbrnch_temp,*my_root;
   int *iatm_typ_brnch_flag;
   int *root_atm_typ;
   int *iexl_temp;                      /*    Scratch space         */
   int *jexl_temp;                      /* to be malloced and freed */
   int *nexl_atm;
   int *jtemp;
   int *iexl_atm_off;
   int *jatm_add;

/*=======================================================================*/
/*          Local Pointer declarations                                  */

   int *brnch_atm_list; 
   int *brnch_atm_root;
   int *root_atm_list;
   int *root_atm_map;
   int *brnch_atm_map;  /* Pointers into brnch_root: Assigned after malloc */
   int *nbrnch_of_root;
   int *nbrnch_of_root_big;
   int **ibrnch_of_root;
   int *natm_add;
   int *iatm_add_off;
   int *iatm_add;
   int *jexl_off_root;
   int *jexl_root;
   int *num_exl_root;

   int int_res_ter   = class->energy_ctrl.int_res_ter;
   int myid          = class->communicate.myid;
   int natm_tot      = class->clatoms_info.natm_tot;
   int natm_typ      = class->atommaps.natm_typ;
   int *iatm_atm_typ = class->atommaps.iatm_atm_typ;
   NAME *atm_typ     = class->atommaps.atm_typ;

   double *cutoff         = class->interact.cutoff;
   double skin            = class->interact.skin;
   double *cutoff_res     = class->interact.cutoff_res;
   double brnch_root_skin = class->interact.brnch_root_skin;
   double brnch_root_cut  = class->interact.brnch_root_cut;

   int *jexl    = bonded->excl.j;
   int *jexl_off    = bonded->excl.j_off;
   int *num_exl = bonded->excl.num;
   int nexl     = bonded->excl.nlst;
   

/*========================================================================*/
/*    I) Write to the screen                                              */

  if(myid==0){
    printf("\n");
    PRINT_LINE_STAR;
    printf("Performing preliminary tasks for branch_root_list_opt\n");
    PRINT_LINE_DASH; printf("\n");
  }/*endif*/

/*=========================================================================*/
/*    II) Find the Valence of each atom using the bond and group cons lists*/
/*        to provide the bonding pattern                                   */

  nvalence    = (int *)cmalloc(sizeof(int)*natm_tot)-1;

  count_valence(bonded,natm_tot,nvalence);

/*========================================================================*/
/*    III) Find and Count the Root and Brnch atoms:                       */
/*          An atom is a branch if it is valence 1 and bonded to another  */
/*          atom with a bond length less than brnch_root_cut.             */
/*          The atom it is bonded to is called a root.                    */
/*          my_root[i] is the index of the root of atom ``i''.            */
/*          nbrnch_temp[i] is the number of branch atoms assigned to      */
/*          each root atom. If ``i'' is a branch atom nbranch_temp[i]     */
/*          is zero.                                                      */
/*          nbrnch_tot is # of branch atoms                               */
/*          nroot_tot  is # of root atoms                                 */
/*          nbrnch_root_max is max number of branches off a root          */

  nbrnch_temp = (int *)cmalloc(sizeof(int)*natm_tot)-1;
  my_root     = (int *)cmalloc(sizeof(int)*natm_tot)-1;

  count_root_branch_data(bonded,natm_tot,brnch_root_cut,nvalence,
                         &nbrnch_tot,&nroot_tot,my_root,nbrnch_temp,
                         &nbrnch_root_max);


  if(nroot_tot==natm_tot && myid==0){
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");    
    printf("You have requested brnch_root_list_opt.  However, \n");
    printf("all the atoms are root atoms. This will result in\n");
    printf("additional overhead with no gain in efficiency.\n");
    printf("Are you certain this is what you would like to do?\n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
  }/*endif*/

  if(myid==0){
    printf("You have %d branch atoms and %d root atoms\n",nbrnch_tot,
                                                          nroot_tot);
  }/*endif*/

/*========================================================================*/
/*    IV) Malloc and Store the branch root data set                       */
/*          brnch_atm_map   = map from clatoms_pos to brnch_atm_list      */
/*          brnch_atm_list  = list of brnch atom                          */
/*          brnch_atm_root  = root of each brnch atom                     */
/*          root_atm_list   = list of root atoms                          */
/*          nbranch_of_root = number of branches of each root             */
/*          ibranch_of_root = branches of each root                       */

  class->nbr_list.brnch_root.nbrnch_root_max = nbrnch_root_max;
  class->nbr_list.brnch_root.nbrnch_tot    = nbrnch_tot;
  class->nbr_list.brnch_root.nroot_tot     = nroot_tot;
  class->nbr_list.brnch_root.brnch_atm_map = 
                           (int *)cmalloc(sizeof(int)*natm_tot)-1;
  class->nbr_list.brnch_root.nbrnch_of_root_big = 
                           (int *)cmalloc(sizeof(int)*natm_tot)-1;
  class->nbr_list.brnch_root.brnch_atm_list = 
                           (int *)cmalloc(sizeof(int)*nbrnch_tot)-1;
  class->nbr_list.brnch_root.brnch_atm_root = 
                           (int *)cmalloc(sizeof(int)*nbrnch_tot)-1;
  class->nbr_list.brnch_root.root_atm_map = 
                           (int *)cmalloc(sizeof(int)*natm_tot)-1;
  class->nbr_list.brnch_root.root_atm_list = 
                           (int *)cmalloc(sizeof(int)*nroot_tot)-1;
  class->nbr_list.brnch_root.nbrnch_of_root = 
                           (int *)cmalloc(sizeof(int)*nroot_tot)-1;
  class->nbr_list.brnch_root.ibrnch_of_root = 
                           cmall_int_mat(1,nroot_tot,1,nbrnch_root_max);

  now_memory =  (2*natm_tot + 2*nroot_tot + 2*nbrnch_tot
                +nroot_tot*nbrnch_root_max)*sizeof(int)*1.0e-06;

  brnch_atm_list = class->nbr_list.brnch_root.brnch_atm_list;
  brnch_atm_root = class->nbr_list.brnch_root.brnch_atm_root;
  brnch_atm_map  = class->nbr_list.brnch_root.brnch_atm_map;
  root_atm_list  = class->nbr_list.brnch_root.root_atm_list;
  root_atm_map   = class->nbr_list.brnch_root.root_atm_map;
  nbrnch_of_root = class->nbr_list.brnch_root.nbrnch_of_root;
  ibrnch_of_root = class->nbr_list.brnch_root.ibrnch_of_root;
  nbrnch_of_root_big = class->nbr_list.brnch_root.nbrnch_of_root_big;

  jtemp          = (int *)cmalloc(sizeof(int)*natm_tot)-1;

  store_branch_root_data(bonded,brnch_atm_list,brnch_atm_root,brnch_atm_map,
                         root_atm_list,root_atm_map,nbrnch_of_root,
                         ibrnch_of_root,brnch_root_cut,
                         natm_tot,nvalence,nbrnch_tot,nroot_tot,
                         my_root,nbrnch_temp,jtemp,
                         nbrnch_of_root_big,myid);

#ifdef DEBUG_BRNCH_ROOT

  if(myid==0){
    printf("I have this many branches %d\n",nbrnch_tot);
    printf("I have this many roots %d\n",nroot_tot);
    scanf("%d",&iii);

    printf("Branch atm list\n");
    scanf("%d",&iii);
    for(i=1;i<=nbrnch_tot;i++){
      printf("%d %d\n",brnch_atm_list[i],i);
    }/*endfor*/

    printf("Branch atm root\n");
    scanf("%d",&iii);
    for(i=1;i<=nbrnch_tot;i++){
      printf("%d %d\n",brnch_atm_root[i],i);
    }/*endfor*/

    printf("Root atm list\n");
    scanf("%d",&iii);
    for(i=1;i<=nroot_tot;i++){
      printf("%d %d\n",root_atm_list[i],i);
    }/*endfor*/

    printf("Root atm map\n");
    scanf("%d",&iii);
    for(i=1;i<=natm_tot;i++){
      printf("%d %d\n",root_atm_map[i],i);
    }/*endfor*/

    printf("Branch atm map\n");
    scanf("%d",&iii);
    for(i=1;i<=natm_tot;i++){
      printf("%d %d\n",brnch_atm_map[i],i);
    }/*endfor*/

    printf("Nbrnch of root\n");
    scanf("%d",&iii);
    for(i=1;i<=nroot_tot;i++){
      printf("%d %d\n",nbrnch_of_root[i],i);
    }/*endfor*/

    printf("ibrnch_of_root\n");
    scanf("%d",&iii);
    for(i=1;i<=nroot_tot;i++){
      for(j=1;j<=nbrnch_of_root[i];j++){
        printf("%d %d %d\n",ibrnch_of_root[i][j],i,j);
      }/*endfor*/
    }/*endfor*/
  }/*endif*/
#endif

/*=========================================================================*/
/*  V) Check the cutoffs of the Branch-Branch and Root-Branch              */
/*     interactions to make sure they are consistent with                  */
/*     the list generation scheme                                          */
  
  iatm_typ_brnch_flag = (int *)cmalloc(sizeof(int)*natm_typ)-1;
  root_atm_typ        = (int *)cmalloc(sizeof(int)*natm_typ)-1;

/*----------------------------------------------------------------------*/
/* i) Compute max brnch-root separation                                */
  find_max_brnch_root_dist(class->clatoms_pos,nbrnch_tot,brnch_atm_list,
                           brnch_atm_root,&(brnch_root_dist_max),
                           class->clatoms_info.pi_beads);
  class->nbr_list.brnch_root.brnch_root_dist_max = brnch_root_dist_max ;

  if(myid==0){ 
   printf("Your maximum brnch root separation is %g A\n",
                                  brnch_root_dist_max*BOHR);
  }/*endif*/

/*----------------------------------------------------------------------*/
/* ii) Check the cutoffs                                                */
  ires_flag = 0;
  check_root_branch_cutoffs(natm_tot,natm_typ,iatm_atm_typ,atm_typ,cutoff,
                            skin,brnch_root_cut,brnch_root_skin,
                            nbrnch_tot,brnch_atm_list,brnch_atm_root,
                            iatm_typ_brnch_flag,root_atm_typ,myid,
                            brnch_root_dist_max,
                            &cutdif_min_brnch_root,
                            &cutdif_min_brnch_brnch,
                            &cutskin_bb_min,&r2_nocheck_dist,
                            ires_flag);
  class->nbr_list.brnch_root.cutdif_min_brnch_root  = cutdif_min_brnch_root;
  class->nbr_list.brnch_root.cutdif_min_brnch_brnch = cutdif_min_brnch_brnch;
  class->nbr_list.brnch_root.cutskin_bb_min = cutskin_bb_min;
  class->nbr_list.brnch_root.r2_nocheck_dist = r2_nocheck_dist;

  if(myid==0){ 
   printf("Your brnch-root  safety margin is %g A\n",
      (cutdif_min_brnch_root -  brnch_root_dist_max+brnch_root_skin)*BOHR);
   printf("Your brnch-brnch safety margin is %g A\n",
      (cutdif_min_brnch_brnch-2*brnch_root_dist_max+brnch_root_skin)*BOHR);
   printf("Your pare down no-check distance is %g A\n",
      sqrt(r2_nocheck_dist)*BOHR);
  }/*endif*/

/*----------------------------------------------------------------------*/
/* iii) Check the RESPA cutoffs                                         */
  if(int_res_ter==1){

   ires_flag = 1;
   check_root_branch_cutoffs(natm_tot,natm_typ,iatm_atm_typ,atm_typ,cutoff_res,
                             skin,brnch_root_cut,brnch_root_skin,
                             nbrnch_tot,brnch_atm_list,brnch_atm_root,
                             iatm_typ_brnch_flag,root_atm_typ,myid,
                             brnch_root_dist_max,
                             &cutdif_min_brnch_root_res,
                             &cutdif_min_brnch_brnch_res,
                             &cutskin_bb_min_res,&r2_nocheck_dist_res,
                             ires_flag);
   class->nbr_list.brnch_root.cutdif_min_brnch_root_res  = 
                              cutdif_min_brnch_root_res;
   class->nbr_list.brnch_root.cutdif_min_brnch_brnch_res = 
                              cutdif_min_brnch_brnch_res;
   class->nbr_list.brnch_root.cutskin_bb_min_res = cutskin_bb_min_res;
   class->nbr_list.brnch_root.r2_nocheck_dist_res = r2_nocheck_dist_res;

   if(myid==0){ 
    printf("Your RESPA brnch-root  safety margin is %g A\n",
      (cutdif_min_brnch_root_res -  brnch_root_dist_max+brnch_root_skin)*BOHR);
    printf("Your RESPA brnch-brnch safety margin is %g A\n",
      (cutdif_min_brnch_brnch_res-2*brnch_root_dist_max+brnch_root_skin)*BOHR);
    printf("Your RESPA pare down no-check distance is %g A\n",
      sqrt(r2_nocheck_dist_res)*BOHR);
   }/*endif*/

  }/*endif:RESPA*/

/*=========================================================================*/
/* VI) Create a full form exclusion list to make comparisons easier        */
/*     i.e. 1-2 and 2-1 exclusions both appear                             */

  nexl_temp    = 2*nexl+natm_tot;
  iexl_temp    = (int *)cmalloc(sizeof(int)*nexl_temp)-1;
  jexl_temp    = (int *)cmalloc(sizeof(int)*nexl_temp)-1;
  nexl_atm     = (int *)cmalloc(sizeof(int)*natm_tot)-1;
  iexl_atm_off = (int *)cmalloc(sizeof(int)*natm_tot)-1;

  create_full_excl_list(natm_tot,nexl,num_exl,jexl,
                       nexl_temp,iexl_temp,jexl_temp,nexl_atm,  
                       iexl_atm_off,jtemp);

#ifdef DEBUG_EXCL_EXP
if(myid==0){
      printf("The total number of exclusions are %d\n",nexl_temp);
      scanf("%d",&iii);
      for(i=1;i<=205;i++){
        printf("Atom %d has %d excl and %d offset\n",i,nexl_atm[i],
                                                     iexl_atm_off[i]);
       for(j=1;j<=nexl_atm[i];j++){
        printf("Atom %d has excl %d \n",i,jexl_temp[j+iexl_atm_off[i]]);
       }/*endfor*/
      if(i%20==0){
        scanf("%d",&iii);
      }
      }/*endfor*/
}/*endif*/
#endif


/*=========================================================================*/
/* VII) Construct the addition list using the branches of the roots and    */
/*      the exclusion list.                                                */

/*-------------------------------------------------------------------------*/
/*  i) Count the addition list and check the scheme                        */

      class->nbr_list.brnch_root.natm_add = 
                           (int *)cmalloc(sizeof(int)*natm_tot)-1;
      class->nbr_list.brnch_root.iatm_add_off = 
                           (int *)cmalloc(sizeof(int)*natm_tot)-1;

      now_memory +=  (2*natm_tot)*sizeof(int)*1.0e-06;


      count_addition_list(&(class->nbr_list.brnch_root),
                          nexl_atm,iexl_atm_off,jexl_temp,natm_tot,jtemp,
                          myid);

/*-------------------------------------------------------------------------*/
/*  ii) Malloc, Local pointerize and Fill the addition list                */

      natm_add_tot = class->nbr_list.brnch_root.natm_add_tot;
      class->nbr_list.brnch_root.iatm_add = 
                           (int *)cmalloc(sizeof(int)*natm_add_tot)-1;
      iatm_add     =  class->nbr_list.brnch_root.iatm_add;
      natm_add     =  class->nbr_list.brnch_root.natm_add;
      iatm_add_off =  class->nbr_list.brnch_root.iatm_add_off;
      jatm_add     =  (int *)cmalloc(sizeof(int)*natm_add_tot)-1;

      fill_addition_list(&(class->nbr_list.brnch_root),nexl_atm,
                           iexl_atm_off,jexl_temp,natm_tot);

          
#ifdef DEBUG_ADD_LIST
    if(myid==0){
      printf("The total number of additions are %d\n",natm_add_tot);
      scanf("%d",&iii);
      for(i=1;i<=natm_tot;i++){
        printf("Atom %d has %d add and %d offset\n",i,natm_add[i],
                                                     iatm_add_off[i]);
       for(j=1;j<=natm_add[i];j++){
        printf("Atom %d has add %d \n",i,iatm_add[j+iatm_add_off[i]]);
       }/*endfor*/
       scanf("%d",&iii);
      }/*endfor*/
    }/*endif*/
#endif
      
/*-------------------------------------------------------------------------*/
/*  iii) Eliminate repeating pairs from the addition list                  */

   tidy_addition_list(natm_tot,&natm_add_tot,natm_add,iatm_add_off,
                      iatm_add,jatm_add,jtemp);

   if(myid==0){
     printf("You have %d additions\n",natm_add_tot);
   }/*endif*/

   class->nbr_list.brnch_root.natm_add_tot = natm_add_tot;
   class->nbr_list.brnch_root.iatm_add     = 
            (int *)crealloc(&(class->nbr_list.brnch_root.iatm_add[1]),
                                natm_add_tot*sizeof(int))-1;
   now_memory +=  (natm_add_tot)*sizeof(int)*1.0e-06;

#ifdef DEBUG_SORT_ADD_LIST
  if(myid==0){
   for(i=1;i<=natm_add_tot;i++){
     printf("list %d %d %d\n",i,jatm_add[i],iatm_add[i]);
     if(i%20==0){scanf("%d",&iii);}
   }/*endfor*/
  }/*endif*/
#endif

/*=========================================================================*/
/*   VIII) Calculate the root atom - root atom exclusion list              */
/*         with root atoms indexed from 1-nroot_tot                        */

  bonded->excl.j_off_root = (int *)cmalloc(sizeof(int)*nroot_tot)-1;
  bonded->excl.num_root   = (int *)cmalloc(sizeof(int)*nroot_tot)-1;
  jexl_off_root           = bonded->excl.j_off_root;
  num_exl_root            = bonded->excl.num_root;

  ic = 0;
  for(i=1;i<=nroot_tot;i++){
    jexl_off_root[i] = ic;
    num_exl_root[i]  = 0;
    ind              = root_atm_list[i];
    for(j=1;j<=num_exl[ind];j++){
      if(root_atm_map[jexl[(j+jexl_off[ind])]]>0){num_exl_root[i]++;}
    }/*endfor*/
    ic+=num_exl_root[i];
  }/*endfor*/

  bonded->excl.nlst_root = ic;
  bonded->excl.j_root    = (int *)cmalloc(sizeof(int)*ic)-1;
  jexl_root              = bonded->excl.j_root;

  ic = 0;
  for(i=1;i<=nroot_tot;i++){
    ind = root_atm_list[i];
    for(j=1;j<=num_exl[ind];j++){
      ktemp = root_atm_map[jexl[(j+jexl_off[ind])]];
      if(ktemp>0){
       ic++;
       jexl_root[ic] = root_atm_map[jexl[(j+jexl_off[ind])]];
      }/*endif*/
    }/*endfor*/
  }/*endfor*/

/*=========================================================================*/
/*   IX) Free memory                                                    */

  cfree(&(nbrnch_temp[1]));
  cfree(&(nvalence[1]));
  cfree(&(my_root[1]));
  cfree(&(iatm_typ_brnch_flag[1]));
  cfree(&(root_atm_typ[1])); 
  cfree(&(jtemp[1]));        
  cfree(&(iexl_temp[1])); 
  cfree(&(jatm_add[1])); 
  cfree(&(jexl_temp[1])); 
  cfree(&(nexl_atm[1])); 
  cfree(&(iexl_atm_off[1])); 

/*=========================================================================*/
/*   VIII) Write out to the screen                                         */

  class->tot_memory+= now_memory;
  if(myid==0){
    printf("Root-branch allocation: %g Mbytes; Total memory: %g Mbytes\n",
                      now_memory,class->tot_memory);
    printf("\n");
    PRINT_LINE_DASH;
    printf("Completed preliminary brnch_root_list_opt tasks.\n");
    PRINT_LINE_STAR;
    printf("\n");
  }/*endif*/

/*========================================================================*/
    }/*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void count_valence(BONDED *bonded,int natm_tot,int *nvalence)

/*==========================================================================*/
  { /*begin routine */
/*==========================================================================*/
/*     Local variable declarations */

   int i,iii;

/*    Local pointers              */
   int *j1_con    = bonded->bond.j1_con;
   int *j2_con    = bonded->bond.j2_con;
   int *j1_pow    = bonded->bond.j1_pow;
   int *j2_pow    = bonded->bond.j2_pow;
   int ncon       = bonded->bond.ncon;
   int npow       = bonded->bond.npow;

   int num_21     = bonded->grp_bond_con.num_21;
   int num_23     = bonded->grp_bond_con.num_23;
   int num_33     = bonded->grp_bond_con.num_33;
   int num_43     = bonded->grp_bond_con.num_43;
   int num_46     = bonded->grp_bond_con.num_46;
   int *j1_21     = bonded->grp_bond_con.j1_21;
   int *j1_23     = bonded->grp_bond_con.j1_23;
   int *j1_33     = bonded->grp_bond_con.j1_33;
   int *j1_43     = bonded->grp_bond_con.j1_43;
   int *j1_46     = bonded->grp_bond_con.j1_46;
   int *j2_21     = bonded->grp_bond_con.j2_21;
   int *j2_23     = bonded->grp_bond_con.j2_23;
   int *j2_33     = bonded->grp_bond_con.j2_33;
   int *j2_43     = bonded->grp_bond_con.j2_43;
   int *j2_46     = bonded->grp_bond_con.j2_46;
   int *j3_23     = bonded->grp_bond_con.j3_23;
   int *j3_33     = bonded->grp_bond_con.j3_33;
   int *j3_43     = bonded->grp_bond_con.j3_43;
   int *j3_46     = bonded->grp_bond_con.j3_46;
   int *j4_43     = bonded->grp_bond_con.j4_43;
   int *j4_46     = bonded->grp_bond_con.j4_46;

   int num_33_watts = bonded->grp_bond_watts.num_33;
   int *j1_33_watts = bonded->grp_bond_watts.j1_33;
   int *j2_33_watts = bonded->grp_bond_watts.j2_33;
   int *j3_33_watts = bonded->grp_bond_watts.j3_33;


/*==========================================================================*/
/* I) Initialize the valence                                                */

  for(i=1;i<=natm_tot;i++){
   nvalence[i]=0;
  }/*endfor*/

/*==========================================================================*/
/* II) Use connectivity to calculate the valence of each atom               */

/*----------------------------------------------------------------------*/
/*         i) Con Bonds                                                */

  for(i=1;i<=ncon;i++){
   nvalence[j1_con[i]]++;
   nvalence[j2_con[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         ii) Pow Bonds                                               */

  for(i=1;i<=npow;i++){
   nvalence[j1_pow[i]]++;
   nvalence[j2_pow[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         iii) Grp 21 Bonds                                            */

  for(i=1;i<=num_21;i++){
   nvalence[j1_21[i]]++;
   nvalence[j2_21[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         iii) Grp 23 Bonds                                            */

  for(i=1;i<=num_23;i++){
   nvalence[j1_23[i]] += 2;
   nvalence[j2_23[i]]++;
   nvalence[j3_23[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         iv) Grp 33 Bonds                                            */

  for(i=1;i<=num_33;i++){
   nvalence[j1_33[i]] += 2;
   nvalence[j2_33[i]]++;
   nvalence[j3_33[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         v) Grp 43 Bonds                                            */

  for(i=1;i<=num_43;i++){
   nvalence[j1_43[i]] +=3;
   nvalence[j2_43[i]]++;
   nvalence[j3_43[i]]++;
   nvalence[j4_43[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         vi) Grp 46 Bonds                                            */

  for(i=1;i<=num_46;i++){
   nvalence[j1_46[i]] +=3;
   nvalence[j2_46[i]]++;
   nvalence[j3_46[i]]++;
   nvalence[j4_46[i]]++;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         vi) Watts 33 Bonds                                            */

  for(i=1;i<=num_33_watts;i++){
   nvalence[j1_33_watts[i]] += 2;
   nvalence[j2_33_watts[i]]++;
   nvalence[j3_33_watts[i]]++;
  }/*endfor*/


/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void count_root_branch_data(BONDED *bonded,int natm_tot,double brnch_root_cut,
                            int *nvalence,
                            int *nbrnch_tot,int *nroot_tot,
                            int *my_root,int *nbrnch_temp,
                            int *nbrnch_root_max)

/*==========================================================================*/
  { /*begin routine */
/*==========================================================================*/
/*     Local variable declarations */
  
  int i,igo,iii;

/*    Local pointers              */
   double *eq_con = bonded->bond.eq_con;
   double *eq_pow = bonded->bond.eq_pow;
   int *j1_con    = bonded->bond.j1_con;
   int *j2_con    = bonded->bond.j2_con;
   int *j1_pow    = bonded->bond.j1_pow;
   int *j2_pow    = bonded->bond.j2_pow;
   int *jtyp_con  = bonded->bond.jtyp_con;
   int *jtyp_pow  = bonded->bond.jtyp_pow;
   int ncon       = bonded->bond.ncon;
   int npow       = bonded->bond.npow;

   int num_21     = bonded->grp_bond_con.num_21;
   int num_23     = bonded->grp_bond_con.num_23;
   int num_33     = bonded->grp_bond_con.num_33;
   int num_43     = bonded->grp_bond_con.num_43;
   int num_46     = bonded->grp_bond_con.num_46;
   int *j1_21     = bonded->grp_bond_con.j1_21;
   int *j1_23     = bonded->grp_bond_con.j1_23;
   int *j1_33     = bonded->grp_bond_con.j1_33;
   int *j1_43     = bonded->grp_bond_con.j1_43;
   int *j1_46     = bonded->grp_bond_con.j1_46;
   int *j2_21     = bonded->grp_bond_con.j2_21;
   int *j2_23     = bonded->grp_bond_con.j2_23;
   int *j2_33     = bonded->grp_bond_con.j2_33;
   int *j2_43     = bonded->grp_bond_con.j2_43;
   int *j2_46     = bonded->grp_bond_con.j2_46;
   int *j3_23     = bonded->grp_bond_con.j3_23;
   int *j3_33     = bonded->grp_bond_con.j3_33;
   int *j3_43     = bonded->grp_bond_con.j3_43;
   int *j3_46     = bonded->grp_bond_con.j3_46;
   int *j4_43     = bonded->grp_bond_con.j4_43;
   int *j4_46     = bonded->grp_bond_con.j4_46;
   double **eq_21 = bonded->grp_bond_con.eq_21;
   double **eq_23 = bonded->grp_bond_con.eq_23;
   double **eq_33 = bonded->grp_bond_con.eq_33;
   double **eq_43 = bonded->grp_bond_con.eq_43;
   double **eq_46 = bonded->grp_bond_con.eq_46;
   int *jtyp_21   = bonded->grp_bond_con.jtyp_21;
   int *jtyp_23   = bonded->grp_bond_con.jtyp_23;
   int *jtyp_33   = bonded->grp_bond_con.jtyp_33;
   int *jtyp_43   = bonded->grp_bond_con.jtyp_43;
   int *jtyp_46   = bonded->grp_bond_con.jtyp_46;

   int num_33_watts     = bonded->grp_bond_watts.num_33;
   int *j1_33_watts     = bonded->grp_bond_watts.j1_33;
   int *j2_33_watts     = bonded->grp_bond_watts.j2_33;
   int *j3_33_watts     = bonded->grp_bond_watts.j3_33;
   double **eq_33_watts = bonded->grp_bond_watts.eq_33;
   int *jtyp_33_watts   = bonded->grp_bond_watts.jtyp_33;

/*==========================================================================*/
/* 0) Initialize */

  for(i=1;i<=natm_tot;i++){
   nbrnch_temp[i]=1;
   my_root[i] = i;
  }/*endif*/

/*==========================================================================*/
/*  I) Use the connectivity to find roots and branches                      */

/*----------------------------------------------------------------------*/
/*         i) Con Bonds                                                 */

  (*nbrnch_tot) = 0;
  for(i=1;i<=ncon;i++){
   igo = 0;
   if((nvalence[j1_con[i]]==1)&&(eq_con[jtyp_con[i]]<=brnch_root_cut)){
    nbrnch_temp[j1_con[i]] = 0;
    nbrnch_temp[j2_con[i]] += 1;
    my_root[j1_con[i]] = j2_con[i];
    (*nbrnch_tot)++;
    igo = 1;
   }/*endif*/
   if(igo==0){
    if((nvalence[j2_con[i]]==1)&&(eq_con[jtyp_con[i]]<=brnch_root_cut)){
     nbrnch_temp[j2_con[i]] = 0;
     nbrnch_temp[j1_con[i]] += 1;
     my_root[j2_con[i]] = j1_con[i];
     (*nbrnch_tot)++;
    }/*endif*/
   }/*endif*/
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*         ii) Pow Bonds                                                */

  for(i=1;i<=npow;i++){
   igo = 0;
   if((nvalence[j1_pow[i]]==1)&&(eq_pow[jtyp_pow[i]]<=brnch_root_cut)){
    nbrnch_temp[j1_pow[i]] = 0;
    nbrnch_temp[j2_pow[i]] += 1;
    my_root[j1_pow[i]] = j2_pow[i];
    (*nbrnch_tot)++;
    igo = 1;
   }/*endif*/
   if(igo==0){
    if((nvalence[j2_pow[i]]==1)&&(eq_pow[jtyp_pow[i]]<=brnch_root_cut)){
     nbrnch_temp[j2_pow[i]] = 0;
     nbrnch_temp[j1_pow[i]] += 1;
     my_root[j2_pow[i]] = j1_pow[i];
     (*nbrnch_tot)++;
    }/*endif*/
   }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        iii) 21 Grp Bonds                                             */

  for(i=1;i<=num_21;i++){

   if((nvalence[j2_21[i]]==1)&&(eq_21[1][jtyp_21[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_21[i]] = 0;
    nbrnch_temp[j1_21[i]] += 1;
    my_root[j2_21[i]] = j1_21[i];
    (*nbrnch_tot)++;
   }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        iii) 23 Grp Bonds                                             */

  for(i=1;i<=num_23;i++){
   if((nvalence[j2_23[i]]==1)&&(eq_23[1][jtyp_23[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_23[i]] = 0;
    nbrnch_temp[j1_23[i]] += 1;
    my_root[j2_23[i]] = j1_23[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j3_23[i]]==1)&&(eq_23[2][jtyp_23[i]]<=brnch_root_cut)){
    nbrnch_temp[j3_23[i]] = 0;
    nbrnch_temp[j1_23[i]] += 1;
    my_root[j3_23[i]] = j1_23[i];
    (*nbrnch_tot)++;
   }/*endif*/
  }/*endfor*/


/*----------------------------------------------------------------------*/
/*        iv) 33 Grp Bonds                                              */

  for(i=1;i<=num_33;i++){
   if((nvalence[j2_33[i]]==1)&&(eq_33[1][jtyp_33[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_33[i]] = 0;
    nbrnch_temp[j1_33[i]] += 1;
    my_root[j2_33[i]] = j1_33[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j3_33[i]]==1)&&(eq_33[2][jtyp_33[i]]<=brnch_root_cut)){
    nbrnch_temp[j3_33[i]] = 0;
    nbrnch_temp[j1_33[i]] += 1;
    my_root[j3_33[i]] = j1_33[i];
    (*nbrnch_tot)++;
   }/*endif*/
  }/*endfor*/


/*----------------------------------------------------------------------*/
/*        v) 43 Grp Bonds                                              */

  for(i=1;i<=num_43;i++){
   if((nvalence[j2_43[i]]==1)&&(eq_43[1][jtyp_43[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_43[i]] = 0;
    nbrnch_temp[j1_43[i]] += 1;
    my_root[j2_43[i]] = j1_43[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j3_43[i]]==1)&&(eq_43[2][jtyp_43[i]]<=brnch_root_cut)){
    nbrnch_temp[j3_43[i]] = 0;
    nbrnch_temp[j1_43[i]] += 1;
    my_root[j3_43[i]] = j1_43[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j4_43[i]]==1)&&(eq_43[3][jtyp_43[i]]<=brnch_root_cut)){
    nbrnch_temp[j4_43[i]] = 0;
    nbrnch_temp[j1_43[i]] += 1;
    my_root[j4_43[i]] = j1_43[i];
    (*nbrnch_tot)++;
   }/*endif*/
  }/*endfor*/


/*----------------------------------------------------------------------*/
/*        vi) 46 Grp Bonds                                              */

  for(i=1;i<=num_46;i++){
   if((nvalence[j2_46[i]]==1)&&(eq_46[1][jtyp_46[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_46[i]] = 0;
    nbrnch_temp[j1_46[i]] += 1;
    my_root[j2_46[i]] = j1_46[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j3_46[i]]==1)&&(eq_46[2][jtyp_46[i]]<=brnch_root_cut)){
    nbrnch_temp[j3_46[i]] = 0;
    nbrnch_temp[j1_46[i]] += 1;
    my_root[j3_46[i]] = j1_46[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j4_46[i]]==1)&&(eq_46[3][jtyp_46[i]]<=brnch_root_cut)){
    nbrnch_temp[j4_46[i]] = 0;
    nbrnch_temp[j1_46[i]] += 1;
    my_root[j4_46[i]] = j1_46[i];
    (*nbrnch_tot)++;
   }/*endif*/
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        vii) 33 Grp watts Bonds                                       */

  for(i=1;i<=num_33_watts;i++){
   if((nvalence[j2_33_watts[i]]==1)&&
      (eq_33_watts[1][jtyp_33_watts[i]]<=brnch_root_cut)){
    nbrnch_temp[j2_33_watts[i]] = 0;
    nbrnch_temp[j1_33_watts[i]] += 1;
    my_root[j2_33_watts[i]] = j1_33_watts[i];
    (*nbrnch_tot)++;
   }/*endif*/

   if((nvalence[j3_33_watts[i]]==1)&&
      (eq_33_watts[2][jtyp_33_watts[i]]<=brnch_root_cut)){
    nbrnch_temp[j3_33_watts[i]] = 0;
    nbrnch_temp[j1_33_watts[i]] += 1;
    my_root[j3_33_watts[i]] = j1_33_watts[i];
    (*nbrnch_tot)++;
   }/*endif*/
  }/*endfor*/

/*=======================================================================*/
/* III) Count up roots and find max number of branches off a root        */

  (*nroot_tot) = 0;      
  (*nbrnch_root_max) = 0;
  for(i=1;i<=natm_tot;i++){
   if(nbrnch_temp[i]>0){(*nroot_tot)++;}
   (*nbrnch_root_max) = MAX((*nbrnch_root_max),nbrnch_temp[i]);
  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void store_branch_root_data(BONDED *bonded,int *brnch_atm_list,
                         int *brnch_atm_root,int *brnch_atm_map,
                         int *root_atm_list,int *root_atm_map,
                         int *nbrnch_of_root,
                         int **ibrnch_of_root,double brnch_root_cut,
                         int natm_tot,int *nvalence,int nbrnch_tot,
                         int nroot_tot,int *my_root,int *nbrnch_temp,
                         int *jtemp,int *nbrnch_of_root_big,int myid)

/*==========================================================================*/
/*  Begin routine */
   {/*begin routine */
/*=======================================================================*/
/*          Local variable declarations                                  */

  int ic,i,igo,ktemp,j,iii;


/*    Local pointers              */
   double *eq_con = bonded->bond.eq_con;
   double *eq_pow = bonded->bond.eq_pow;
   int *j1_con    = bonded->bond.j1_con;
   int *j2_con    = bonded->bond.j2_con;
   int *j1_pow    = bonded->bond.j1_pow;
   int *j2_pow    = bonded->bond.j2_pow;
   int *jtyp_con  = bonded->bond.jtyp_con;
   int *jtyp_pow  = bonded->bond.jtyp_pow;
   int ncon       = bonded->bond.ncon;
   int npow       = bonded->bond.npow;

   int num_21     = bonded->grp_bond_con.num_21;
   int num_23     = bonded->grp_bond_con.num_23;
   int num_33     = bonded->grp_bond_con.num_33;
   int num_43     = bonded->grp_bond_con.num_43;
   int num_46     = bonded->grp_bond_con.num_46;
   int *j1_21     = bonded->grp_bond_con.j1_21;
   int *j1_23     = bonded->grp_bond_con.j1_23;
   int *j1_33     = bonded->grp_bond_con.j1_33;
   int *j1_43     = bonded->grp_bond_con.j1_43;
   int *j1_46     = bonded->grp_bond_con.j1_46;
   int *j2_21     = bonded->grp_bond_con.j2_21;
   int *j2_23     = bonded->grp_bond_con.j2_23;
   int *j2_33     = bonded->grp_bond_con.j2_33;
   int *j2_43     = bonded->grp_bond_con.j2_43;
   int *j2_46     = bonded->grp_bond_con.j2_46;
   int *j3_23     = bonded->grp_bond_con.j3_23;
   int *j3_33     = bonded->grp_bond_con.j3_33;
   int *j3_43     = bonded->grp_bond_con.j3_43;
   int *j3_46     = bonded->grp_bond_con.j3_46;
   int *j4_43     = bonded->grp_bond_con.j4_43;
   int *j4_46     = bonded->grp_bond_con.j4_46;
   double **eq_21 = bonded->grp_bond_con.eq_21;
   double **eq_23 = bonded->grp_bond_con.eq_23;
   double **eq_33 = bonded->grp_bond_con.eq_33;
   double **eq_43 = bonded->grp_bond_con.eq_43;
   double **eq_46 = bonded->grp_bond_con.eq_46;
   int *jtyp_21   = bonded->grp_bond_con.jtyp_21;
   int *jtyp_23   = bonded->grp_bond_con.jtyp_23;
   int *jtyp_33   = bonded->grp_bond_con.jtyp_33;
   int *jtyp_43   = bonded->grp_bond_con.jtyp_43;
   int *jtyp_46   = bonded->grp_bond_con.jtyp_46;

   int num_33_watts     = bonded->grp_bond_watts.num_33;
   int *j1_33_watts     = bonded->grp_bond_watts.j1_33;
   int *j2_33_watts     = bonded->grp_bond_watts.j2_33;
   int *j3_33_watts     = bonded->grp_bond_watts.j3_33;
   double **eq_33_watts = bonded->grp_bond_watts.eq_33;
   int *jtyp_33_watts   = bonded->grp_bond_watts.jtyp_33;

/*========================================================================*/
/*    I) Create the list of branches and roots as well as                 */
/*        Map of Roots and Branches to clatoms_pos                        */

  ic = 0;
  for(i=1;i<=natm_tot;i++){
    brnch_atm_map[i] = 0;
   if(nbrnch_temp[i]==0){
    ic++;
    brnch_atm_list[ic] = i;
    brnch_atm_map[i] = ic;
    brnch_atm_root[ic] = my_root[i];
   }/*endif*/
  }/*endfor*/
  if(ic!=nbrnch_tot&&myid==0){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Internal Error in control_brnch_root_list.c\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

  ic = 0;
  for(i=1;i<=natm_tot;i++){
    root_atm_map[i] = 0;
   if(nbrnch_temp[i]>0){
    ic++;
    root_atm_list[ic] = i;
    root_atm_map[i] = ic;
   }/*endif*/
  }/*endfor*/
  if(ic!=nroot_tot&&myid==0){
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     printf("Internal Error in control_brnch_root_list.c\n");
     printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     fflush(stdout);
     exit(1);
  }/*endif*/

/*========================================================================*/
/*  II) Create a list of the root atoms branches                          */

/*----------------------------------------------------------------------*/
/*        i) Initialize                                                */

  for(i=1;i<=nroot_tot;i++){
    nbrnch_of_root[i] = 0;
  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        ii) Con Bonds                                                 */

  for(i=1;i<=ncon;i++){

       igo = 0;
       if((nvalence[j1_con[i]]==1)&&(eq_con[jtyp_con[i]]<=brnch_root_cut)){
        ktemp = root_atm_map[j2_con[i]];
        nbrnch_of_root[ktemp] += 1;
        ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j1_con[i];
        igo = 1;
       }/*endif*/

      if(igo==0){
       if((nvalence[j2_con[i]]==1)&&(eq_con[jtyp_con[i]]<=brnch_root_cut)){
        ktemp = root_atm_map[j1_con[i]];
        nbrnch_of_root[ktemp] += 1;
        ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_con[i];
       }/*endif*/
      }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        ii) Pow Bonds                                                 */

  for(i=1;i<=npow;i++){

      igo = 0;
      if((nvalence[j1_pow[i]]==1)&&(eq_pow[jtyp_pow[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j2_pow[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j1_pow[i];
       igo = 1;
      }/*endif*/

      if(igo==0){
       if((nvalence[j2_pow[i]]==1)&&(eq_pow[jtyp_pow[i]]<=brnch_root_cut)){
        ktemp = root_atm_map[j1_pow[i]];
        nbrnch_of_root[ktemp] += 1;
        ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_pow[i];
       }/*endif*/
      }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        iii) Grp 21 Bonds                                             */

  for(i=1;i<=num_21;i++){

      if((nvalence[j2_21[i]]==1)&&(eq_21[1][jtyp_21[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_21[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_21[i];
      }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        iii) Grp 23 Bonds                                             */

  for(i=1;i<=num_23;i++){

      if((nvalence[j2_23[i]]==1)&&(eq_23[1][jtyp_23[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_23[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_23[i];
      }/*endif*/

      if((nvalence[j3_23[i]]==1)&&(eq_23[2][jtyp_23[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_23[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j3_23[i];
      }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        iv) Grp 33 Bonds                                             */

  for(i=1;i<=num_33;i++){

      if((nvalence[j2_33[i]]==1)&&(eq_33[1][jtyp_33[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_33[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_33[i];
      }/*endif*/

      if((nvalence[j3_33[i]]==1)&&(eq_33[2][jtyp_33[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_33[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j3_33[i];
      }/*endif*/

  }/*endfor*/

/*----------------------------------------------------------------------*/
/*        v) Grp 43 Bonds                                             */

  for(i=1;i<=num_43;i++){

      if((nvalence[j2_43[i]]==1)&&(eq_43[1][jtyp_43[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_43[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_43[i];
      }/*endif*/

      if((nvalence[j3_43[i]]==1)&&(eq_43[2][jtyp_43[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_43[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j3_43[i];
      }/*endif*/

      if((nvalence[j4_43[i]]==1)&&(eq_43[3][jtyp_43[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_43[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j4_43[i];
      }/*endif*/

  }/*endfor*/


/*----------------------------------------------------------------------*/
/*        vi) Grp 46 Bonds                                             */

  for(i=1;i<=num_46;i++){

      if((nvalence[j2_46[i]]==1)&&(eq_46[1][jtyp_46[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_46[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_46[i];
      }/*endif*/

      if((nvalence[j3_46[i]]==1)&&(eq_46[2][jtyp_46[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_46[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j3_46[i];
      }/*endif*/

      if((nvalence[j4_46[i]]==1)&&(eq_46[3][jtyp_46[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_46[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j4_46[i];
      }/*endif*/

  }/*endfor*/

    
/*----------------------------------------------------------------------*/
/*        iv) Grp 33 watts Bonds                                        */

  for(i=1;i<=num_33_watts;i++){
      if((nvalence[j2_33_watts[i]]==1)&&
         (eq_33_watts[1][jtyp_33_watts[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_33_watts[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j2_33_watts[i];
      }/*endif*/

      if((nvalence[j3_33_watts[i]]==1)&&
         (eq_33_watts[2][jtyp_33_watts[i]]<=brnch_root_cut)){
       ktemp = root_atm_map[j1_33_watts[i]];
       nbrnch_of_root[ktemp] += 1;
       ibrnch_of_root[ktemp][nbrnch_of_root[ktemp]] = j3_33_watts[i];
      }/*endif*/

  }/*endfor*/

/*========================================================================*/
/*  III) Sort the branches of each atom                                   */

  for(i=1;i<=nroot_tot;i++){

    if(nbrnch_of_root[i] > 1){
      for(j=1;j<=nbrnch_of_root[i];j++){
       jtemp[j] = ibrnch_of_root[i][j];
      }/*endfor*/
      small_excl_sort(nbrnch_of_root[i],jtemp);
      for(j=1;j<=nbrnch_of_root[i];j++){
         ibrnch_of_root[i][j] = jtemp[j];
      }/*endfor*/
    }/*endif*/

  }/*endfor*/

/*========================================================================*/
/*  III) Fill big list of number of branches of each root                 */
/*       Here the atom itself counts as a branch. Branch atoms have zero  */
/*       branches */

  for(i=1;i<=natm_tot;i++){nbrnch_of_root_big[i]=0;}
  for(i=1;i<=nroot_tot;i++){
   ktemp  = root_atm_list[i];
   nbrnch_of_root_big[ktemp] = nbrnch_of_root[i]+1;
  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void check_root_branch_cutoffs(int natm_tot,int natm_typ,int *iatm_atm_typ,
                               NAME *atm_typ,double *cutoff,double skin,
                               double brnch_root_cut,double brnch_root_skin,
                               int nbrnch_tot,int *brnch_atm_list,
                               int *brnch_atm_root,
                               int *iatm_typ_brnch_flag,int *root_atm_typ,
                               int myid,double brnch_root_dist_max,
                               double *cutdif_min_brnch_root_ret,
                               double *cutdif_min_brnch_brnch_ret,
                               double *cutskin_bb_min_ret,
                               double *r2_nocheck_dist,
                               int ires_flag)


/*==========================================================================*/
/*  Begin routine */
   {/*begin routine */
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,j,ktemp,natm1_typ2,ityp,iitemp,jjtemp,kktemp,ind_typ,iii;
  double cuttrr,cuttij,diff;
  int ind_root_i,ind_root_j,ierr;
  double brnch_root_cut_use;
  double cutdif_min_brnch_root;
  double cutdif_min_brnch_brnch;
  double cutskin_bb_min;

  brnch_root_cut_use = MAX(brnch_root_cut,brnch_root_dist_max);

/*=======================================================================*/
/*  I) Find the root and branch atom types                                 */

  for(i=1;i<=natm_typ;i++){
    iatm_typ_brnch_flag[i] = 0;
    root_atm_typ[i] = i;
  }/*endfor*/

  for(i=1;i<=nbrnch_tot;i++){
     ktemp = brnch_atm_list[i];
     jjtemp = iatm_atm_typ[ktemp];
     iatm_typ_brnch_flag[jjtemp]++;
     root_atm_typ[jjtemp] = iatm_atm_typ[brnch_atm_root[i]];
  }/*endfor*/

#ifdef DEBUG_ATM_TYPE
if(myid==0){
  printf("Atm type, branch flag and root atm type\n");

  for(i=1;i<=natm_typ;i++){
      printf("%d %d %d\n",i,iatm_typ_brnch_flag[i],root_atm_typ[i]);
  }/*endfor*/
  scanf("%d",&iii);
}/*endif*/
#endif

/*=======================================================================*/
/* II) Loop over the interaction types and check the cutoffs             */

#ifdef DEBUG_CUTOFF_CHCK
if(myid==0){
  printf("Cutoff information\n");
}/*endif*/
#endif


  cutdif_min_brnch_root = 10.0;
  cutdif_min_brnch_brnch = 10.0;
  cutskin_bb_min = cutoff[1]+skin;
  natm1_typ2 = (natm_typ+1) * 2;
  ityp = 0;
  ierr = 0;
  for(i=1;i<=natm_typ;i++){
    for(j=i;j<=natm_typ;j++){

    /*----------------------------------------------------------------*/
    /* A) Find the present interaction cutoff and the root-root cutoff*/

      ityp++;        
      cuttij = cutoff[ityp];
      ind_root_i = root_atm_typ[i];
      ind_root_j = root_atm_typ[j];
      jjtemp      = MIN(ind_root_i,ind_root_j); 
      iitemp      = MAX(ind_root_i,ind_root_j); 
      kktemp      = ((jjtemp - 1) * (natm1_typ2 - jjtemp))/2; 
      ind_typ     = (kktemp + iitemp - jjtemp + 1); 
      cuttrr = cutoff[ind_typ];
      if(ind_root_j != j && ind_root_i != i){
        cutskin_bb_min = MIN(cutskin_bb_min,cuttij+skin);
      }/*endif*/
      diff = 0.0;
#ifdef DEBUG_CUTOFF_CHCK
      if(myid==0){
       printf("%d atm typ %s-%s, root_i %d root_j %d ind %d\n", 
              ityp,atm_typ[i],atm_typ[j],ind_root_i,ind_root_j,ind_typ);
      }/*endif*/
#endif
    /*-----------------------------------------------------------------*/
    /* B) Root-Branch  check                                           */

      if((iatm_typ_brnch_flag[i]==0)&&(iatm_typ_brnch_flag[j]>0)){
       diff = (cuttij - (cuttrr-brnch_root_cut_use+brnch_root_skin));
       cutdif_min_brnch_root = MIN((cuttrr - cuttij),cutdif_min_brnch_root);
       if(cuttij > (cuttrr-brnch_root_cut_use+brnch_root_skin)){
         diff*= BOHR;
         ierr++;
         if(myid==0){
          if(ierr==1){
           printf("\n@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          }/*endif*/
          if(ires_flag==0){
           printf("The %s-%s interaction %s-%s cutoff too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);
          }else{
           printf(
            "The %s-%s interaction %s-%s RESPA cutoff too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);

          }/*endif*/
         }/*endif*/
       }/*endif*/       
      }/*endif*/       

    /*-----------------------------------------------------------------*/
    /* C) Branch-Root  check                                           */

      if((iatm_typ_brnch_flag[i]>0)&&(iatm_typ_brnch_flag[j]==0)){
       diff = (cuttij - (cuttrr-brnch_root_cut_use+brnch_root_skin));
       cutdif_min_brnch_root = MIN((cuttrr - cuttij),cutdif_min_brnch_root);
       if(cuttij > (cuttrr-brnch_root_cut_use+brnch_root_skin)){        
         diff*= BOHR;
         ierr++;
         if(myid==0){
          if(ierr==1){
           printf("\n@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          }/*endif*/
          if(ires_flag==0){
           printf("The %s-%s interaction cutoff (%s-%s) too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);
          }else{
           printf(
             "The %s-%s interaction RESPA cutoff (%s-%s) too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);
          }/*endif*/
         }/*endif*/
       }/*endif*/       
      }/*endif*/       

    /*-----------------------------------------------------------------*/
    /* D) Branch-Branch check                                         */

      if((iatm_typ_brnch_flag[i]>0)&&(iatm_typ_brnch_flag[j]>0)){
       diff = (cuttij - (cuttrr-2*brnch_root_cut_use+brnch_root_skin));
       cutdif_min_brnch_brnch = MIN((cuttrr - cuttij),cutdif_min_brnch_brnch);
       if(cuttij > (cuttrr-2*brnch_root_cut_use+brnch_root_skin)){        
         diff*= BOHR;
         ierr++;
         if(myid==0){
          if(ierr==1){
           printf("\n@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          }/*endif*/
          if(ires_flag==0){
           printf("The %s-%s interaction cutoff (%s-%s) too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);
          }else{
           printf(
             "The %s-%s interaction RESPA (%s-%s) cutoff too large by %g A\n",
                                          atm_typ[i],atm_typ[j],
                                          atm_typ[ind_root_i],
                                          atm_typ[ind_root_j],
                                          diff);
          }/*endif*/
        }/*endif*/
       }/*endif*/       
      }/*endif*/       
    /*-----------------------------------------------------------------*/
    /* E) Increment safety                                             */
    }/*endfor:j*/
  }/*endfor:i*/

/*==========================================================================*/
/* III) Finish the error message                                            */

  if(ierr>0){
   if(myid==0){
      printf("-----------------------------------------------\n");
      printf("You may either increase the brnch_root_skin by \n");
      printf("the maximum difference or decrease each cutoff \n");
      printf("by the precise amount specified.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
   }/*endif*/       
      exit(1);
  }/*endif*/       

/*==========================================================================*/
/* IV) Return the safety                                                    */

  *cutdif_min_brnch_root_ret  = cutdif_min_brnch_root;
  *cutdif_min_brnch_brnch_ret = cutdif_min_brnch_brnch;
  *cutskin_bb_min_ret = cutskin_bb_min;
  if(cutskin_bb_min - 2*brnch_root_cut_use > 0.0){
    *r2_nocheck_dist = (cutskin_bb_min - 2*brnch_root_cut_use)
                      *(cutskin_bb_min - 2*brnch_root_cut_use);
  }else{
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    printf("Your minimum branch-branch cutoff is small: %g\n",cutskin_bb_min*BOHR);
    printf("This will result in additional overhead with no gain in efficiency.\n");
    printf("It is probably better to let the branches interact with a normal cutoff.\n");
    printf("Are you certain this is what you would like to do?\n");
    printf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
    *r2_nocheck_dist = 0.0;
  }/*endif*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/







/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void create_full_excl_list(int natm_tot,int nexl,int *num_exl,int *jexl,
                           int nexl_temp,int *iexl_temp,int *jexl_temp,
                           int *nexl_atm,int *iexl_atm_off,int *jtemp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
   {/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

int ioff,i,j,koff,iii;

/*=======================================================================*/
/* I) Creat a pair exclustion list with all pairs (i-j and j-i)          */

  ioff = 0;     
  for(i=1;i<=natm_tot;i++){
    for(j=ioff+1;j<=ioff+num_exl[i];j++){
      iexl_temp[j] = i;
      jexl_temp[j] = jexl[j];
    }/*endfor*/
    ioff += num_exl[i];
  }/*endfor*/

  ioff = 0;     
  koff = nexl;
  for(i=1;i<=natm_tot;i++){
    for(j=ioff+1;j<=ioff+num_exl[i];j++){
      iexl_temp[j+koff] = jexl[j];
      jexl_temp[j+koff] = i;
    }/*endfor*/
    ioff += num_exl[i];
  }/*endfor*/

  koff = 2*nexl;
  for(i=1;i<=natm_tot;i++){
    iexl_temp[i+koff] = i; 
    jexl_temp[i+koff] = i;
  }/*endfor*/ 

/*=======================================================================*/
/* II) Sort doubled list and count exclusions of each atom               */

  if(nexl_temp > 1){
    big_excl_sort(nexl_temp, iexl_temp,jexl_temp);
  }/*endif*/

  for(i=1;i<=natm_tot;i++){
    nexl_atm[i] = 0;
  }/*endfor*/
  for(i=1;i<=nexl_temp;i++){
    nexl_atm[iexl_temp[i]] += 1;
  }/*endfor*/

  iexl_atm_off[1] = 0;
  for(i=2;i<=natm_tot;i++){
    iexl_atm_off[i] = iexl_atm_off[i-1]+nexl_atm[i-1];
  }/*endfor*/

/*=======================================================================*/
/*  III) Sort the exclusions of each individual atom                     */

  for(i=1;i<=natm_tot;i++){

    if(nexl_atm[i] > 1){
      for(j=1;j<=nexl_atm[i];j++){
        jtemp[j] = jexl_temp[j+iexl_atm_off[i]];           
      }/*endfor*/
      small_excl_sort(nexl_atm[i], jtemp);
      for(j=1;j<=nexl_atm[i];j++){
        jexl_temp[j+iexl_atm_off[i]] = jtemp[j];           
      }/*endfor*/
    }/*endif*/

  }/*endfor*/

/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void count_addition_list(BRNCH_ROOT *brnch_root,int *nexl_atm ,int *jexl_off,
                         int *jexl,int natm_tot, int *ifound,int myid)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int ioff_add,i,k,j,my_root,ioff_me,ioff_root;
  int j_root,j_me,k_root,k_me,ind_now,j_brnch;
  int itemp,ind_root,ind_brnch,iii;
  int nfound;

/*          Local pointer declarations                                  */

  int *iatm_add_off    = brnch_root->iatm_add_off;
  int *iatm_add        = brnch_root->iatm_add;
  int *natm_add        = brnch_root->natm_add;
  int *brnch_atm_map   = brnch_root->brnch_atm_map;
  int *brnch_atm_root  = brnch_root->brnch_atm_root;
  int **ibrnch_of_root = brnch_root->ibrnch_of_root;
  int *root_atm_map    = brnch_root->root_atm_map;
  int *nbrnch_of_root  = brnch_root->nbrnch_of_root;

/*=======================================================================*/
/* I) Initialize the offset and number of additions per atom             */

  iatm_add_off[1] = 0;
  for(i=1;i<=natm_tot;i++){
     natm_add[i] = 0;
  }/*endfor*/

/*=======================================================================*/
/* II) Loop over the atoms and find the additions                        */

  for(i=1;i<=natm_tot;i++){

     my_root = i;
     ioff_me = jexl_off[i];
/*-------------------------------------------------------------------*/
/*  i) Branch atom additions only                                    */
   
     if(brnch_atm_map[i]>0){
         nfound = 0;
/*    a) Find my root      */
         ind_brnch = brnch_atm_map[i];
         my_root = brnch_atm_root[ind_brnch];      
/*    b) Compare my exclusions to those of my root and add extras to   */
/*       my addition list                                              */
         ioff_root = jexl_off[my_root];
         j_root = 1;
         j_me   = 1;
         while((j_me <= nexl_atm[i]) && (j_root <= nexl_atm[my_root])){
            if(jexl[j_root+ioff_root] >  jexl[j_me+ioff_me]){itemp=0;}
            if(jexl[j_root+ioff_root] <  jexl[j_me+ioff_me]){itemp=1;}
            if(jexl[j_root+ioff_root] == jexl[j_me+ioff_me]){itemp=2;}
            if(itemp == 0){j_me++;}
            if(itemp == 2){j_me++;j_root++;nfound++;}
            if(itemp == 1){natm_add[i]++;j_root++;}
         }/*endwhile:comparing exclusions*/
         if(nfound!=nexl_atm[i]&&myid==0){
            printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            printf("A branch atom has more excluded interactions   \n");
            printf("than its root. The brnch_root_list scheme will \n");
            printf("therefore fail as presently implemented.       \n");
            printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
            fflush(stdout);
            exit(1);
         }/*endif*/
         for(j=j_root;j<=nexl_atm[my_root];j++){
           natm_add[i]++;
         }/*endfor*/
     }/*endif:I am a branch atom*/

/*-------------------------------------------------------------------*/
/*  ii) Find exclusions of root atoms which are root atoms           */
/*      These will be excluded by root-root exclusion rules          */
/*      Branches of root atoms are also exclusions of the root       */

     if(root_atm_map[i]>0){

        for(j=1;j<=nexl_atm[i];j++){
          ifound[j]  = 0;
          if(root_atm_map[jexl[j+ioff_me]]){ifound[j]=1;}
        }/*endif*/
        ind_root = root_atm_map[i];
        j_brnch = 1;
        j_me   = 1;
        while((j_me <= nexl_atm[i]) && (j_brnch <= nbrnch_of_root[ind_root])){
           if(ibrnch_of_root[ind_root][j_brnch] > 
                               jexl[j_me+ioff_me]){itemp=0;}
           if(ibrnch_of_root[ind_root][j_brnch] < 
                               jexl[j_me+ioff_me]){itemp=1;}
           if(ibrnch_of_root[ind_root][j_brnch] == 
                               jexl[j_me+ioff_me]){itemp=2;}
           if(itemp == 0){j_me++;}
           if(itemp == 2){ifound[j_me]++;j_me++;j_brnch++;}
           if(itemp == 1){j_brnch++;}
        }/*endwhile:comparing exclusion and branche list*/

     }/*endif: i am a root*/

/*-------------------------------------------------------------------*/
/*  iii) Branch and root atom additions: Loop over the branches (k)  */
/*      of each exclusion (j) of my root atom and compare            */
/*      to my exclusion list. Add to addition list if necessary      */
/*      Check exclusion generation if I am a root atoms              */

     ioff_root = jexl_off[my_root];
     for(j=1;j<=nexl_atm[my_root];j++){

         ind_now = jexl[j+ioff_root];
         if((root_atm_map[ind_now] > 0)&&(ind_now != my_root)){
            ind_root = root_atm_map[ind_now];
            k_root = 1;
            k_me   = 1;
            while((k_root <= nbrnch_of_root[ind_root]) && 
                  (k_me <= nexl_atm[i])){
              if(ibrnch_of_root[ind_root][k_root] > 
                               jexl[k_me+ioff_me]){itemp=0;}
              if(ibrnch_of_root[ind_root][k_root] < 
                               jexl[k_me+ioff_me]){itemp=1;}
              if(ibrnch_of_root[ind_root][k_root] == 
                               jexl[k_me+ioff_me]){itemp=2;}
              if(itemp == 0){k_me++;}
              if(itemp == 2){
                if(root_atm_map[i]>0){ifound[k_me]++;}
                k_me++;k_root++;
              }/*endif*/
              if(itemp == 1){natm_add[i]++;k_root++;}
            }/*endwhile:comparing exclusions and branches*/
            for(k=k_root;k<=nbrnch_of_root[ind_root];k++){
              natm_add[i]++;
            }/*endfor*/
         }/*endif: jth exclusion is itself a root*/             

     }/*endfor:jth exclusion of my root */

     if(root_atm_map[i]>0){
      for(j=1;j<=nexl_atm[i];j++){
        if(ifound[j]!=1&&myid==0){
          printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
          printf("Cannot find all root atom excluded interactions\n");
          printf("by considering root atom-root atom exclusions  \n");
          printf("the branches of root atom exclusions and branches\n");
          printf("of the root atom itself. Therefore, the brnch_root_list\n");
          printf("scheme will fail as presently implemented: %d %d %d %d\n",
                  i,j,ifound[j],jexl[j+ioff_me]);
           printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
           fflush(stdout);
           exit(1);
        }/*endif*/
      }/*endfor*/
     }/*endif*/

/*-------------------------------------------------------------------*/
/*  iv) Increment addition list offset                               */

     if(i != natm_tot){iatm_add_off[i+1] = iatm_add_off[i]+natm_add[i];}

   }/*endfor: ith atom*/

/*-------------------------------------------------------------------*/
/*  v) Define the total number of additions                          */

   brnch_root->natm_add_tot = iatm_add_off[natm_tot] +  natm_add[natm_tot];

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void fill_addition_list(BRNCH_ROOT *brnch_root,int *nexl_atm ,int *jexl_off,
                         int *jexl,int natm_tot)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,k,j,my_root,ioff_me,ioff_root;
  int j_root,j_me,k_root,k_me,ind_now;
  int itemp,ind_root,ind_brnch;

/*         Local pointers                                               */
  int *iatm_add_off    = brnch_root->iatm_add_off;
  int *iatm_add        = brnch_root->iatm_add;
  int *natm_add        = brnch_root->natm_add;
  int *brnch_atm_map   = brnch_root->brnch_atm_map;
  int *brnch_atm_root  = brnch_root->brnch_atm_root;
  int **ibrnch_of_root = brnch_root->ibrnch_of_root;
  int *root_atm_map    = brnch_root->root_atm_map;
  int *nbrnch_of_root  = brnch_root->nbrnch_of_root;

/*=======================================================================*/
/* I) Initialize the offset and number of additions per atom             */

  iatm_add_off[1] = 0;
  for(i=1;i<=natm_tot;i++){
     natm_add[i] = 0;
  }/*endfor*/

/*=======================================================================*/
/* II) Loop over the atoms and find the additions                        */

  for(i=1;i<=natm_tot;i++){
     my_root = i;

/*-------------------------------------------------------------------*/
/*  i) Branch atom additions only                                    */

     if(brnch_atm_map[i]>0){

/*    a) Find my root      */
         ind_brnch = brnch_atm_map[i];
         my_root = brnch_atm_root[ind_brnch];      

/*    b) Compare my exclusions to those of my root and add them to my */
/*       addition list                                                */
         ioff_me = jexl_off[i];
         ioff_root = jexl_off[my_root];
         j_root = 1;
         j_me   = 1;
         while((j_me <= nexl_atm[i]) && (j_root <= nexl_atm[my_root])){
           if(jexl[j_root+ioff_root] >  jexl[j_me+ioff_me]){itemp=0;}
           if(jexl[j_root+ioff_root] <  jexl[j_me+ioff_me]){itemp=1;}
           if(jexl[j_root+ioff_root] == jexl[j_me+ioff_me]){itemp=2;}
           if(itemp == 0){j_me++;}
           if(itemp == 2){j_me++;j_root++;}
           if(itemp == 1){
            natm_add[i]++;
            iatm_add[natm_add[i]+iatm_add_off[i]] = jexl[j_root+ioff_root];
            j_root++;
           }/*endif*/
         }/*endwhile:comparing exclusions*/
         for(j=j_root;j<=nexl_atm[my_root];j++){
           natm_add[i]++;
           iatm_add[natm_add[i]+iatm_add_off[i]] = jexl[j+ioff_root];
         }/*endfor: add to addition list */

     }/*endif: I am a branch atom*/
           
/*-------------------------------------------------------------------*/
/*  ii) Branch and root atom additions: Loop over the branches (k)   */
/*      of each exclusion (j) of my root atom and compare            */
/*      to my exclusion list. Add to addition list if necessary      */

     ioff_me = jexl_off[i];
     ioff_root = jexl_off[my_root];
     for(j=1;j<=nexl_atm[my_root];j++){

         ind_now = jexl[j+ioff_root];
         if((root_atm_map[ind_now] > 0)&&(ind_now != my_root)){
            ind_root = root_atm_map[ind_now];
            k_me = 1;
            k_root = 1;
            while((k_root <= nbrnch_of_root[ind_root]) && 
                  (k_me <= nexl_atm[i])){
              if(ibrnch_of_root[ind_root][k_root] > 
                               jexl[k_me+ioff_me]){itemp=0;}
              if(ibrnch_of_root[ind_root][k_root] < 
                               jexl[k_me+ioff_me]){itemp=1;}
              if(ibrnch_of_root[ind_root][k_root] == 
                               jexl[k_me+ioff_me]){itemp=2;}
              if(itemp == 0){k_me++;}
              if(itemp == 2){k_me++;k_root++;}
              if(itemp == 1){
                natm_add[i]++;
                iatm_add[natm_add[i]+iatm_add_off[i]] = 
                      ibrnch_of_root[ind_root][k_root];
                k_root++;
              }/*endif*/
            }/*endwhile:comparing exclusions and branches*/
            for(k=k_root;k<=nbrnch_of_root[ind_root];k++){
             natm_add[i]++;
             iatm_add[natm_add[i]+iatm_add_off[i]] = 
                          ibrnch_of_root[ind_root][k];       
            }/*endfor: add to list */
         }/*endif:jth exclusion is a root*/             

     }/*endfor: j exclusion of my root*/

/*-------------------------------------------------------------------*/
/*  iii) Increment list offset                                       */

     if(i != natm_tot){iatm_add_off[i+1] = iatm_add_off[i]+natm_add[i];}

  }/*endfor: i atom*/

/*-----------------------------------------------------------------------*/
    }/*end routine*/ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void tidy_addition_list(int natm_tot,int *natm_add_tot,
                        int *natm_add,int *iatm_add_off,
                        int *iatm_add,int *jatm_add,int *jtemp)

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int i,j,iitemp,jjtemp,iii;

/*=======================================================================*/
/* I) Create the pair addition list                                      */

  for(i=1;i<=natm_tot;i++){
    for(j=1;j<=natm_add[i];j++){
      jjtemp = iatm_add[j+iatm_add_off[i]];
      iitemp = i;
      if(jjtemp >= iitemp){
         jatm_add[j+iatm_add_off[i]] = jjtemp;
         iatm_add[j+iatm_add_off[i]] = iitemp;
      }/*endif*/
      if(jjtemp < iitemp){
         iatm_add[j+iatm_add_off[i]] = jjtemp;
         jatm_add[j+iatm_add_off[i]] = iitemp;
      }/*endif*/
    }/*endfor*/
  }/*endfor*/

/*=======================================================================*/
/* II) Sort jatm_add keeping iatm_add commensurate and eliminate         */
/*         repeated pairs of additions                                   */
       
  if(*natm_add_tot>1){
    exl_sort(natm_add_tot,jatm_add,iatm_add,natm_tot);
  }/*endif*/

/*=======================================================================*/
/*  III) Recount the list and list offsets                               */

  for(i=1;i<=natm_tot;i++){
   natm_add[i] = 0;
  }/*endfor*/
  for(i=1;i<=(*natm_add_tot);i++){
   natm_add[jatm_add[i]]++;
  }/*endfor*/

  iatm_add_off[1]=0;
  for(i=2;i<=natm_tot;i++){
    iatm_add_off[i] = iatm_add_off[(i-1)]+natm_add[(i-1)];
  }/*endfor*/

/*=======================================================================*/
/*  IV) Sort the additions of each individual atom                       */

  for(i=1;i<=natm_tot;i++){

    if(natm_add[i] > 1){
      for(j=1;j<=natm_add[i];j++){
       jtemp[j] = iatm_add[j+iatm_add_off[i]];           
      }/*endfor*/
      small_excl_sort(natm_add[i], jtemp);
      for(j=1;j<=natm_add[i];j++){
       iatm_add[j+iatm_add_off[i]] = jtemp[j];           
      }/*endfor*/
    }/*endif*/

  }/*endfor*/


/*-----------------------------------------------------------------------*/
   }/*end routine*/ 
/*==========================================================================*/




/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void big_excl_sort(int n, int index[],int jndex[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,rjndex,iii;
  int k,*kndex,*mask,isub,temp;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

/*=======================================================================*/
/* II) Sort array index keeping jndex commensurrate */

  for(;;){

/*---------------------------------------------------------------------*/
/*  A)hire rindex */
    if(m>1){ 
      m--;
      rindex = index[m];
      rjndex = jndex[m];
/*--------------------------------------------------------------------*/
/*  B)retire/promote index[1] */
    }else{
      rindex = index[ir];
      rjndex = jndex[ir];
      index[ir]=index[1];
      jndex[ir]=jndex[1];
      ir--;
      if(ir==1){
       index[1]=rindex;
       jndex[1]=rjndex;
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
       jndex[i]=jndex[j];
       i=j;
       j=2*j;
      }else{
       /*    c)if no demotations exit while */
       j=ir+1;
      }/*endif*/
    } /*endwhile*/
    /*    d)slot rindex */
    index[i] = rindex;
    jndex[i] = rjndex;
  }/*endfor*/

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void small_excl_sort(int n, int index[])

/*=======================================================================*/
/*            Begin subprogram:                                          */
{/*begin routine*/
/*=======================================================================*/
/*          Local variable declarations                                  */

  int m,ir,i,j,rindex,iii;
  int k,*kndex,*mask,isub,temp;

/*=======================================================================*/
/* I) Setup                        */

  m  = n/2+1;
  ir = n;

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

/*-----------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/



