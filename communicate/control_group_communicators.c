/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                Module: control_communicate_groups                        */
/*                                                                          */
/* This subprogram reads in all user inputs                                 */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/


#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_communicate_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_communicate_local.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

 void control_group_communicators(CLASS *class,CP *cp,int cp_on)

/*==========================================================================*/
{/* begin routine */ 
/*==========================================================================*/
/*          Local variable declarations                                     */

 int iii;
 int num_proc = class->communicate.np;
 int myid     = class->communicate.myid;
 int comm_compare_result;
 MPI_Comm world;

/*==========================================================================*/
/* 0) Output */

 if(myid==0){
   PRINT_LINE_STAR;
   printf("Setting up communicators and communication pakages\n");
   PRINT_LINE_DASH;printf("\n");
 }/*endif*/

/*==========================================================================*/
/* I) Build different communicators                                        */

 if(num_proc>1){

   build_communicate_groups(&(class->communicate),cp_on);

 }else{

   class->communicate.myid_bead        = 0;
   class->communicate.myid_bead_prime  = 0;
   class->communicate.myid_state       = 0;
   class->communicate.myid_forc        = 0;
   class->communicate.myid_forc_source = 0;
   class->communicate.myid_forc_target = 0;

   Comm_dup(class->communicate.world,&class->communicate.comm_forc);
   Comm_dup(class->communicate.world,&class->communicate.comm_forc_source);
   Comm_dup(class->communicate.world,&class->communicate.comm_forc_target);
   Comm_dup(class->communicate.world,&class->communicate.comm_states);
   Comm_dup(class->communicate.world,&class->communicate.comm_beads);
   Comm_dup(class->communicate.world,&class->communicate.comm_faux);

 }/*endif*/

/*==========================================================================*/
/* II) Duplicate communicate into force package */

 /*-----------------------------------------------------------------------*/
 /* i) Force level communicator */

  class->class_comm_forc_pkg.myid     = class->communicate.myid_forc;
  class->class_comm_forc_pkg.num_proc = class->communicate.np_forc;
  Comm_dup(class->communicate.comm_forc,&class->class_comm_forc_pkg.comm);

  class->class_comm_forc_pkg.dbl_num_proc = 
                            (double)(class->communicate.np_forc);

 /*-----------------------------------------------------------------------*/
 /* ii) Source-Force level communicator */

  class->class_comm_forc_pkg.plimpton_ind.myid_source
                          = class->communicate.myid_forc_source;
  class->class_comm_forc_pkg.plimpton_ind.num_proc_source
                          = class->communicate.np_forc_src;
  Comm_dup(class->communicate.comm_forc_source,
          &(class->class_comm_forc_pkg.plimpton_ind.source));

 /*-----------------------------------------------------------------------*/
 /* iii) Target-Force level communicator */

  class->class_comm_forc_pkg.plimpton_ind.myid_target
                          = class->communicate.myid_forc_target;
  class->class_comm_forc_pkg.plimpton_ind.num_proc_target
                          = class->communicate.np_forc_trg;
  Comm_dup(class->communicate.comm_forc_target,
          &(class->class_comm_forc_pkg.plimpton_ind.target));

/*==========================================================================*/
/* II) Set up some Bead level parallel stuff                                */

  class->clatoms_info.pi_beads_proc_st    =
                     (class->communicate.myid_bead_prime)*
                     (class->clatoms_info.pi_beads_proc)+1;

  class->clatoms_info.pi_beads_proc_end   =
                     (class->communicate.myid_bead_prime+1)*
                     (class->clatoms_info.pi_beads_proc);

  cp->cpcoeffs_info.pi_beads_proc_st  = class->clatoms_info.pi_beads_proc_st;
  cp->cpcoeffs_info.pi_beads_proc_end = class->clatoms_info.pi_beads_proc_end;

/*==========================================================================*/
/* III) Duplicate class->communicate into cp->communicate                   */
/*      then set up the cp_comm_package                                     */

 /*-------------------------------------------------------------------------*/
 /* i) Duplicate */

  cp->communicate.np             = class->communicate.np;
  cp->communicate.np_beads       = class->communicate.np_beads;
  cp->communicate.np_states      = class->communicate.np_states;
  cp->communicate.np_forc        = class->communicate.np_forc;
  cp->communicate.np_forc_src    =  class->communicate.np_forc_src;
  cp->communicate.np_forc_trg    =  class->communicate.np_forc_trg;

  cp->communicate.myid             = class->communicate.myid;
  cp->communicate.myid_bead        = class->communicate.myid_bead;
  cp->communicate.myid_bead_forc   = class->communicate.myid_bead_forc;
  cp->communicate.myid_state       = class->communicate.myid_state;
  cp->communicate.myid_forc        = class->communicate.myid_forc;
  cp->communicate.myid_forc_source = class->communicate.myid_forc_source;
  cp->communicate.myid_forc_target = class->communicate.myid_forc_target;


  if(cp->communicate.np>1){

#ifdef DUPLICATE_NOT_CORRECT
    Comm_dup(class->communicate.world,&(cp->communicate.world));
    Comm_dup(class->communicate.comm_beads,&(cp->communicate.comm_beads));
    Comm_dup(class->communicate.comm_beads_forc,
              &(cp->communicate.comm_beads_forc));
    Comm_dup(class->communicate.comm_states,&(cp->communicate.comm_states));
    Comm_dup(class->communicate.comm_forc,&(cp->communicate.comm_forc));
    Comm_dup(class->communicate.comm_forc_source,
              &(cp->communicate.comm_forc_source));
    Comm_dup(class->communicate.comm_forc_target,
              &(cp->communicate.comm_forc_target));
#endif

    cp->communicate.world            = class->communicate.world;
    cp->communicate.comm_beads       = class->communicate.comm_beads;
    cp->communicate.comm_beads_forc  = class->communicate.comm_beads_forc;
    cp->communicate.comm_states      = class->communicate.comm_states;
    cp->communicate.comm_forc        = class->communicate.comm_forc;
    cp->communicate.comm_forc_source = class->communicate.comm_forc_source;
    cp->communicate.comm_forc_target = class->communicate.comm_forc_target;
    

  }else{

    cp->communicate.world            = class->communicate.world;
    cp->communicate.comm_beads       = class->communicate.comm_beads;
    cp->communicate.comm_beads_forc  = class->communicate.comm_beads_forc;
    cp->communicate.comm_states      = class->communicate.comm_states;
    cp->communicate.comm_forc        = class->communicate.comm_forc;
    cp->communicate.comm_forc_source = class->communicate.comm_forc_source;
    cp->communicate.comm_forc_target = class->communicate.comm_forc_target;
    cp->communicate.comm_faux        = class->communicate.comm_faux;

  }/*endif : cp->communicate.np>1*/

 /*-------------------------------------------------------------------------*/
 /* ii) Build the cp package */

  Comm_dup(class->communicate.world,&world);
  build_cp_comm_pkg(cp,world); /* Contains no src/target stuff */

/*==========================================================================*/
/* Output */

  if(myid==0){
    PRINT_LINE_DASH;
    printf("Completed communicator and package setup\n"); 
    PRINT_LINE_STAR;printf("\n");
  }/*endif*/

/*==========================================================================*/
   }/* end routine */ 
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
 void build_cp_comm_pkg(CP *cp,MPI_Comm world)
/*==========================================================================*/
{/* begin routine */ 
/*==========================================================================*/
/*          Local variable declarations                                     */

 int irem,idiv,iii;

/*==========================================================================*/
/* I) Up states                                                              */

 /*------------------------------------*/ 
 /* i) states per processor stuff      */

  idiv = cp->cpcoeffs_info.nstate_up/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_up % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_up_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_up_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_up_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_up_st = irem*(idiv+1) 
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
  cp->cpcoeffs_info.istate_up_end = cp->cpcoeffs_info.istate_up_st +
                                    cp->cpcoeffs_info.nstate_up_proc-1;

 /*------------------------------------*/
 /* ii) coefs per processor stuff      */

  cp->cp_comm_state_pkg_up.num_proc   = cp->communicate.np_states;
  cp->cp_comm_state_pkg_up.myid       = cp->communicate.myid_state;
  cp->cp_comm_state_pkg_up.nstate     = cp->cpcoeffs_info.nstate_up;
  cp->cp_comm_state_pkg_up.ncoef      = cp->cpcoeffs_info.ncoef;
  cp->cp_comm_state_pkg_up.nstate_proc= cp->cpcoeffs_info.nstate_up_proc;
  cp->cp_comm_state_pkg_up.world      = world;  

  if(cp->communicate.np_states > 1){
    Comm_dup(cp->communicate.comm_states,&(cp->cp_comm_state_pkg_up.comm));
  } else {
    cp->cp_comm_state_pkg_up.comm = cp->communicate.comm_states;
  }/* endif */

  irem             = (cp->cp_comm_state_pkg_up.nstate %
                      cp->cp_comm_state_pkg_up.num_proc);
  cp->cp_comm_state_pkg_up.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_up.nstate_max = (irem > 0 ? 
                          ((idiv+1)*cp->communicate.np_states) : 
                          (idiv*cp->communicate.np_states)) ;

  idiv                  = cp->cpcoeffs_info.ncoef/cp->communicate.np_states;
  irem                  = cp->cpcoeffs_info.ncoef % cp->communicate.np_states;
  cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_up.nstate_ncoef_proc   
                        = (cp->communicate.myid_state < irem ? idiv+1 : idiv);
  if(cp->communicate.myid_state <= irem) {
    cp->cp_comm_state_pkg_up.icoef_start = 
                             cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cp_comm_state_pkg_up.icoef_start = irem*(idiv+1) 
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/

  cp->cpcoeffs_info.nstate_ncoef_proc_up = 
                        cp->cp_comm_state_pkg_up.nstate_ncoef_proc;
  cp->cpcoeffs_info.nstate_ncoef_proc_max_up = 
                        cp->cp_comm_state_pkg_up.nstate_ncoef_proc_max;
  cp->cpcoeffs_info.icoef_start_up = 
                        cp->cp_comm_state_pkg_up.icoef_start;

/*==========================================================================*/
/* II) Down states                                                          */

 /*------------------------------------*/ 
 /* i) states per processor stuff      */

  idiv = cp->cpcoeffs_info.nstate_dn/cp->communicate.np_states;
  irem = (cp->cpcoeffs_info.nstate_dn % cp->communicate.np_states);
  cp->cpcoeffs_info.nstate_dn_proc = idiv;
  if(cp->communicate.myid_state < irem) {
     cp->cpcoeffs_info.nstate_dn_proc = idiv+1;
  }/*endif*/
  if(cp->communicate.myid_state <= irem) {
    cp->cpcoeffs_info.istate_dn_st = cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cpcoeffs_info.istate_dn_st = irem*(idiv+1) 
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/
  cp->cpcoeffs_info.istate_dn_end = cp->cpcoeffs_info.istate_dn_st +
                                    cp->cpcoeffs_info.nstate_dn_proc-1;

 /*------------------------------------*/
 /* ii) coefs per processor stuff      */

  cp->cp_comm_state_pkg_dn.num_proc   = cp->communicate.np_states;
  cp->cp_comm_state_pkg_dn.myid       = cp->communicate.myid_state;
  cp->cp_comm_state_pkg_dn.nstate     = cp->cpcoeffs_info.nstate_dn;
  cp->cp_comm_state_pkg_dn.ncoef      = cp->cpcoeffs_info.ncoef;
  cp->cp_comm_state_pkg_dn.nstate_proc= cp->cpcoeffs_info.nstate_dn_proc;
  cp->cp_comm_state_pkg_dn.world      = world;  
  if(cp->communicate.np_states > 1){
    Comm_dup(cp->communicate.comm_states,&(cp->cp_comm_state_pkg_dn.comm));
  } else {
    cp->cp_comm_state_pkg_dn.comm = cp->communicate.comm_states;
  }/* endif */

  irem             = (cp->cp_comm_state_pkg_dn.nstate %
                      cp->cp_comm_state_pkg_dn.num_proc);
  cp->cp_comm_state_pkg_dn.nstate_proc_max  = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_dn.nstate_max = (irem > 0 ? 
                          ((idiv+1)*cp->communicate.np_states) : 
                          (idiv*cp->communicate.np_states)) ;

  idiv                  = cp->cpcoeffs_info.ncoef/cp->communicate.np_states;
  irem                  = cp->cpcoeffs_info.ncoef % cp->communicate.np_states;
  cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max = (irem > 0 ? idiv+1 : idiv);
  cp->cp_comm_state_pkg_dn.nstate_ncoef_proc   
                        = (cp->communicate.myid_state < irem ? idiv+1 : idiv);
  if(cp->communicate.myid_state <= irem) {
    cp->cp_comm_state_pkg_dn.icoef_start = 
                             cp->communicate.myid_state*(idiv+1)+1;
  } else {
    cp->cp_comm_state_pkg_dn.icoef_start = irem*(idiv+1) 
                                   + (cp->communicate.myid_state-irem)*idiv+1;
  }/*endif*/


  cp->cpcoeffs_info.nstate_ncoef_proc_dn = 
                        cp->cp_comm_state_pkg_dn.nstate_ncoef_proc;
  cp->cpcoeffs_info.nstate_ncoef_proc_max_dn = 
                        cp->cp_comm_state_pkg_dn.nstate_ncoef_proc_max;
  cp->cpcoeffs_info.icoef_start_dn = 
                        cp->cp_comm_state_pkg_dn.icoef_start;

/*==========================================================================*/
   }/* end routine */
/*==========================================================================*/












