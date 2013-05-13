/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: comm_groups.c                                */
/*                                                                          */
/* Subprogram creates group communicators                                   */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_cp.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"

#define MAXPROCS 1000
#define DEBUG_OFF



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*  When CP=on                                                              */
/*     np_forc == np_state OR np_forc == 1                                  */
/*     if(np_forc==np_state) CP and Forc share processors                   */
/*     Path Integral commm is outer and State Level is inner                */
/*     np = np_state * np_bead ALWAYS                                       */
/*     forc processors are proc 0 only or equal to State                    */
/*     forc processors are then split into np_sourc*np_targ                 */
/*                                                                          */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  When CP=off                                                             */
/*     Path Integral commm is outer and Forc Level is inner                 */
/*     np = np_forc * np_bead ALWAYS                                        */
/*     forc processors are then split into np_sourc*np_targ                 */
/*     state processors are proc 0 only                                     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

  void build_communicate_groups(COMMUNICATE *communicate,int cp_on)

/*==========================================================================*/
/*             Begin routine                                                */
    {/*begin routine */
/*==========================================================================*/
#include "../typ_defs/typ_mask.h"

   int myid           = communicate->myid;
   int np_states      = communicate->np_states;
   int np_beads       = communicate->np_beads;
   int np             = communicate->np;
   int np_forc        = communicate->np_forc;
   int np_forc_src = communicate->np_forc_src;
   int np_forc_trg = communicate->np_forc_trg;

   int i,ii,loop,iii,ntemp,idiv,irem,icase,numcomm_bead;
   int *ranks;

   MPI_Comm world;
   MPI_Comm comm_forc;
   MPI_Group excl_group,temp;

/*=======================================================================*/
/* 0) Mallocs and Dups                                                   */

   numcomm_bead = np/np_beads;           /* number of bead communicators */
   if(numcomm_bead==0){numcomm_bead=1;}

   ranks = (int *) cmalloc(MAXPROCS*sizeof(int));

   Comm_dup(communicate->world,&(world));
   Barrier(world);

/*=======================================================================*/
/* I) Get path_integral or Bead level communicators                      */
/*     Bead communciators are OUTER communicators                        */
/*     Split world into forc  and bead pieces if cp is off               */
/*     Split world into state and bead pieces if cp is on                */
/*     Comm_beads_forc is a copy of Comm_beads                           */

   if(cp_on==1){icase=1;}
   if(cp_on==0){icase=2;}

   switch(icase){

     case 1:
        communicate->comm_beads = 
             build_grp_comm_outer(np,np_beads,np_states,myid,
                                &(communicate->myid_bead),ranks,world);
     break;

     case 2:
        communicate->comm_beads = 
             build_grp_comm_outer(np,np_beads,np_forc,myid,
                                &(communicate->myid_bead),ranks,world);
     break;

   }/*switch*/

   Comm_dup(communicate->comm_beads,&(communicate->comm_beads_forc));
   communicate->myid_bead_forc = communicate->myid_bead;

/*=======================================================================*/
/* II) Get state and force level communicators                           */
/*     State and Force comms are INNER communicators                     */
/*     CP is off : np = np_forc*np_bead                                  */

  if((cp_on==1)&&(np_forc==1))        {icase=1;}
  if((cp_on==1)&&(np_forc==np_states)){icase=2;}
  if((cp_on==0))                      {icase=3;}

  switch(icase){

  /*======================================================================*/
  /* CP IS ON : Np_forc = 1 */

   case 1:
 
     /*------------------------------------------------------*/
     /* i) The force level communicator is only proc 0       */

      ranks[0] = myid;
      Comm_group(world,&excl_group);
      Group_incl(excl_group,np_forc,ranks,&temp);
      Comm_create(world,temp,&communicate->comm_forc);
      Comm_dup(communicate->comm_forc,&(communicate->comm_forc_source));
      Comm_dup(communicate->comm_forc,&(communicate->comm_forc_target));
      communicate->myid_forc        = 0;
      communicate->myid_forc_source = 0;
      communicate->myid_forc_target = 0;
      Group_free(&excl_group);

     /*------------------------------------------------------*/
     /* ii) The state level communicator is an INNER         */
      communicate->comm_states =
           build_grp_comm_inner(np,np_beads,np_states,myid,
                              &(communicate->myid_state),ranks,world);

     /*------------------------------------------------*/
     /* iii) The myid_bead_prime and myid_bead store   */
     /*      the id of this proc in the comm_bead, if  */
     /*      it is the FIRST bead level communicator.  */
     /*      In subsequent comm_beads,myid_bead is out */
     /*      of range and myid_bead_prime=myid_state.  */
     /*      myid_bead_forc ALWAYS has carries         */
     /*      the rank of this proc in the bead comm    */
     /*      to which it is associated                 */

      communicate->myid_bead_prime = communicate->myid_bead;
      if(communicate->myid_state!=0){
        communicate->myid_bead       = communicate->np_beads;
      }/*endif*/

     break;

  /*======================================================================*/
   case 2:

     /*------------------------------------------------------*/
     /* i) The state level communicator is an INNER          */

      communicate->comm_states =
           build_grp_comm_inner(np,np_beads,np_states,myid,
                              &(communicate->myid_state),ranks,world);

  
     /*------------------------------------------------------*/
     /* ii) The forc level communicator is an INNER = state  */

       communicate->myid_forc = communicate->myid_state;
       Comm_dup(communicate->comm_states,&(communicate->comm_forc));
       Comm_dup(communicate->comm_states,&(comm_forc));

     /*------------------------------------------------*/
     /* iii) Split the force level communicator        */
     /*     Targets are OUTER and Sources are INNER    */

      communicate->comm_forc_target =
          build_grp_comm_outer(np_forc,np_forc_trg,np_forc_src,
                   communicate->myid_forc,&(communicate->myid_forc_target),
                   ranks,comm_forc);
      communicate->comm_forc_source =
          build_grp_comm_inner(np_forc,np_forc_trg,np_forc_src,
                   communicate->myid_forc,&(communicate->myid_forc_source),
                   ranks,comm_forc);
  
     /*------------------------------------------------*/
     /* iii) The myid_bead_prime and myid_bead store   */
     /*      the id of this proc in the comm_bead, if  */
     /*      it is the FIRST bead level communicator.  */
     /*      In subsequent comm_beads,myid_bead is out */
     /*      of range and myid_bead_prime=myid_state.  */
     /*      myid_bead_forc ALWAYS has carries         */
     /*      the rank of this proc in the bead comm    */
     /*      to which it is associated                 */

      communicate->myid_bead_prime = communicate->myid_bead;
      if(communicate->myid_state!=0){
        communicate->myid_bead       = communicate->np_beads;
      }/*endif*/

     break;

  /*======================================================================*/
   case 3: /* CP is off : force level only */

     /*------------------------------------------------*/
     /* i) The state communicator is only proc 0       */

      ranks[0] = myid;
      Comm_group(world,&excl_group);
      Group_incl(excl_group,np_states,ranks,&temp);
      Comm_create(world,temp,&communicate->comm_states);
      communicate->myid_state = 0;
      Group_free(&excl_group);

     /*------------------------------------------------*/
     /* ii) The force level communicator is an INNER   */

      communicate->comm_forc =
           build_grp_comm_inner(np,np_beads,np_forc,myid,
                              &(communicate->myid_forc),ranks,world);
      Comm_dup(communicate->comm_forc,&(comm_forc));

     /*------------------------------------------------*/
     /* iii) Split the force level communicator        */
     /*     Targets are OUTER and Sources are INNER    */

      communicate->comm_forc_target =
          build_grp_comm_outer(np_forc,np_forc_trg,np_forc_src,
                   communicate->myid_forc,&(communicate->myid_forc_target),
                   ranks,comm_forc);
      communicate->comm_forc_source =
          build_grp_comm_inner(np_forc,np_forc_trg,np_forc_src,
                   communicate->myid_forc,&(communicate->myid_forc_source),
                   ranks,comm_forc);

     /*------------------------------------------------*/
     /* iv) The myid_bead_prime and myid_bead store    */
     /*      the id of this proc in the comm_bead, if  */
     /*      it is the FIRST bead level communicator.  */
     /*      In subsequent comm_beads,myid_bead is out */
     /*      of range and myid_bead_prime=myid_forc.   */
     /*      myid_bead_forc ALWAYS has carries         */
     /*      the rank of this proc in the bead comm    */
     /*      to which it is associated                 */

      communicate->myid_bead_prime = communicate->myid_bead;
      if(communicate->myid_forc!=0){
        communicate->myid_bead       = communicate->np_beads;
      }/*endif*/

   break;

 }/*end switch*/

/*------------------------------------------------------------------------*/
/*  Free the memory                                                       */

      cfree(&ranks[0]);

/*------------------------------------------------------------------------*/
   } /*end routine*/ 
/*==========================================================================*/





/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

MPI_Comm build_grp_comm_outer(int np, int np_outer, int np_inner, int myid, 
                              int *myid_outer,int *ranks, MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                             */
{/*begin routine */
/*=======================================================================*/
/*  Local Variables */

  int iii;
  int myid_inner;

  int i,j,k;
  int ngroup;

  MPI_Group world_group,incl_group;
  MPI_Comm gen_comm;
  MPI_Comm junk_comm;

/*=====================================================================*/
/* Set local inner id's                                                */

  myid_inner = (myid % np_inner);

  if(np_outer*np_inner != np){
    printf("Incorrect number of procs %d vs %d\n",np,np_outer*np_inner);
    Finalize();
    exit(0);
  }/*endfor*/

/*=======================================================================*/
/*             Get rank of processor in new communicator                 */

  Comm_group(world,&world_group);

  for(j=0;j < np_inner;j++){
 /*-----------------------------------------------------------------------*/
 /* i) set the ranks   */

     for(i=0;i<np_outer;i++){
       ranks[i] = np_inner*i+j;
     }/*endfor*/

 /*-----------------------------------------------------------------------*/
 /* ii) Create the new communicator                                      */

     Group_incl(world_group,np_outer,ranks,&incl_group);
     Barrier(world);
     if(myid_inner==j){
       Comm_create(world,incl_group,&gen_comm);
       Barrier(world);
       Comm_rank(gen_comm,myid_outer);
       if(myid!= np_inner*(*myid_outer) + myid_inner){
         printf("Problems in building outer communicator\n");
         printf("ID expected %d not equal to ID given %d\n",myid,
                         np_inner*(*myid_outer) + myid_inner);
         printf("with myid_out %d and myid_in %d\n",*myid_outer,myid_inner);
       }/*endif*/
     }else{
       Comm_create(world,incl_group,&junk_comm);
       Barrier(world);
     }/*endif*/
     Group_free(&incl_group);
     Barrier(world);
   }/*endfor*/

   Group_free(&world_group);

   return gen_comm;

/*------------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/





/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

MPI_Comm build_grp_comm_inner(int np, int np_outer, int np_inner, int myid, 
                              int *myid_inner,int *ranks, MPI_Comm world)

/*=======================================================================*/
/*             Begin routine                                             */
{/*begin routine */
/*=======================================================================*/
/*  Local Variables */

  int iii;
  int myid_outer;

  int i,j,k;

  MPI_Group world_group,incl_group;
  MPI_Comm gen_comm;
  MPI_Comm junk_comm;

/*=====================================================================*/
/* Set local inner id's                                                */

  myid_outer = (myid / np_inner);

  if(np_outer*np_inner != np){
    printf("Incorrect number of procs %d vs %d\n",np,np_outer*np_inner);
    Finalize();
    exit(0);
  }/*endfor*/

/*=======================================================================*/
/*             Get rank of processor in new communicator                 */

  Comm_group(world,&world_group);

  for(j=0;j < np_outer;j++){
 /*-----------------------------------------------------------------------*/
 /* i) set the ranks   */

     for(i=0;i<np_inner;i++){
       ranks[i] = np_inner*j+i;
     }/*endfor*/

 /*-----------------------------------------------------------------------*/
 /* ii) Create the new communicator                                      */

     Group_incl(world_group,np_inner,ranks,&incl_group);
     Barrier(world);
     if(myid_outer==j){
       Comm_create(world,incl_group,&gen_comm);
       Barrier(world);
       Comm_rank(gen_comm,myid_inner);
       if(myid!= myid_outer*np_inner + *myid_inner){
         printf("Problems in building inner communicator\n");
         printf("ID expected %d not equal to ID given %d\n",myid,
                    myid_outer*np_inner + *myid_inner);
         printf("with myid_out %d and myid_in %d\n",myid_outer,*myid_inner);
       }/*endif*/
     }else{
       Comm_create(world,incl_group,&junk_comm);
       Barrier(world);
     }/*endif*/
     Group_free(&incl_group);
     Barrier(world);
   }/*endfor*/

   Group_free(&world_group);

   return gen_comm;

/*------------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/









