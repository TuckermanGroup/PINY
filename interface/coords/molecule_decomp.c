/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_coord                                   */
/*                                                                          */
/* This subprogram reads atm-atm_NHC vol-vol_NHC input for a MD on a        */ 
/* LD-classical potential energy surface (LD-PES)                           */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../proto_defs/proto_coords_entry.h"
#include "../proto_defs/proto_coords_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#define DEBUG_OFF


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void control_molec_decomp(CLASS *class,BONDED *bonded,
                          GENERAL_DATA *general_data)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */
     int i,iii;
     int nvt              = general_data->ensopts.nvt;
     int npt_i            = general_data->ensopts.npt_i;
     int npt_f            = general_data->ensopts.npt_f;
     int pimd_on; 

     pimd_on = general_data->simopts.pimd + general_data->simopts.cp_pimd 
              + general_data->simopts.cp_wave_pimd 
              + general_data->simopts.debug_pimd 
              + general_data->simopts.debug_cp_pimd
              + general_data->simopts.cp_wave_min_pimd;


/*========================================================================*/
/* I) Find the intramolecular stuff for each processor                  */ 

  if(class->communicate.np_forc>1){  
     assign_coord_forc(class);
  }else{
     class->clatoms_info.myatm_start = 1;
     class->clatoms_info.myatm_end   = class->clatoms_info.natm_tot;
     class->class_comm_forc_pkg.myatm_start = 1;
     class->class_comm_forc_pkg.myatm_end   = class->clatoms_info.natm_tot; 
  }/*endelse*/
  bonded->grp_bond_con.num_21_tot = bonded->grp_bond_con.num_21;
  bonded->grp_bond_con.num_23_tot = bonded->grp_bond_con.num_23;
  bonded->grp_bond_con.num_33_tot = bonded->grp_bond_con.num_33;
  bonded->grp_bond_watts.num_33_tot = bonded->grp_bond_watts.num_33;
  bonded->grp_bond_con.num_43_tot = bonded->grp_bond_con.num_43;
  bonded->grp_bond_con.num_46_tot = bonded->grp_bond_con.num_46;
  bonded->bond.ncon_tot = bonded->bond.ncon;
  bonded->bend.ncon_tot = bonded->bend.ncon;
  bonded->tors.ncon_tot = bonded->tors.ncon;
  if(class->communicate.np_forc>1){  
     assign_bonded_forc(bonded,class);
  }/*endif*/ 
  else{
   bonded->bond.ncon_max = bonded->bond.ncon;
  }/*endelse*/

/*========================================================================*/
/* II) Find the thermostat stuff for each processor                  */ 


if(class->clatoms_info.pi_beads_proc_st==1){
 if(nvt + npt_i + npt_f > 0){
  if(class->communicate.np_forc>1){  
     assign_thermo_forc(class);
  }else{
     class->therm_info_class.mytherm_start = 1;
     class->therm_info_class.mytherm_end   = class->therm_info_class.num_nhc;
     class->class_comm_forc_pkg.mytherm_start = 1;
     class->class_comm_forc_pkg.mytherm_end   = 
                                             class->therm_info_class.num_nhc;
     class->therm_info_class.itherm_nshare  = 
           (int *)cmalloc((class->therm_info_class.num_nhc+1)*sizeof(int))-1;
     class->therm_info_class.ditherm_nshare_i  = 
           (double *)cmalloc((class->therm_info_class.num_nhc+1)
                                              *sizeof(double))-1;

     for(i=1;i<=class->therm_info_class.num_nhc;i++){
      class->therm_info_class.itherm_nshare[i] = 1;
      class->therm_info_class.ditherm_nshare_i[i] = 1.0;
     }/*endfor*/

  }/*endelse*/
 }/*endif*/
}/*endif*/

if(pimd_on==1){
 if(nvt + npt_i + npt_f > 0){
  if(class->communicate.np_forc>1){  
     assign_bead_thermo_forc(class);
  }else{
     class->therm_info_bead.mytherm_start = 1;
     class->therm_info_bead.mytherm_end   = class->therm_info_bead.num_nhc;
     class->class_comm_forc_pkg.mytherm_start = 1;
     class->class_comm_forc_pkg.mytherm_end   = 
                                             class->therm_info_bead.num_nhc;
     class->therm_info_bead.itherm_nshare  = 
           (int *)cmalloc((class->therm_info_bead.num_nhc+1)*sizeof(int))-1;
     class->therm_info_bead.ditherm_nshare_i  = 
           (double *)cmalloc((class->therm_info_bead.num_nhc+1)
                                              *sizeof(double))-1;

     for(i=1;i<=class->therm_info_bead.num_nhc;i++){
      class->therm_info_bead.itherm_nshare[i] = 1;
      class->therm_info_bead.ditherm_nshare_i[i] = 1.0;
     }/*endfor*/
  }/*endif*/ 
 }/*endif*/
}/*endif*/

/*========================================================================*/
} /* end routine */
/*==========================================================================*/





/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_coord_forc(CLASS *class)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

  int icount,proc_now,my_mol,i,j,iii;   
  int natm_proc;

/*======================================================================*/
/*               Local Pointer declarations                            */

  MPI_Comm world       = class->communicate.world;
  int myid             = class->communicate.myid;
  int myid_forc        = class->communicate.myid_forc;
  int myid_bead        = class->communicate.myid_bead;
  int pi_beads_proc_st = class->clatoms_info.pi_beads_proc_st;
  int np_forc          = class->communicate.np_forc;
  int natm_tot         = class->clatoms_info.natm_tot;
  int *iatm_mol_num    = class->atommaps.iatm_mol_num;

  int *displs_atm;     /* defined after malloc */
  int *recv_count_atm; /* defined after malloc */
  int *iatm_proc_num;  /* defined after malloc */
  
  class->class_comm_forc_pkg.displs_atm = 
                            (int *) cmalloc(np_forc*sizeof(int))-1;
  class->class_comm_forc_pkg.recv_count_atm = 
                            (int *) cmalloc(np_forc*sizeof(int))-1;
  class->atommaps.iatm_proc_num  = 
                            (int *) cmalloc(natm_tot*sizeof(int))-1;

  displs_atm     = class->class_comm_forc_pkg.displs_atm;
  recv_count_atm = class->class_comm_forc_pkg.recv_count_atm;
  iatm_proc_num  = class->atommaps.iatm_proc_num;


/*======================================================================*/
/* I) Decompose the molecules among the processors                      */

  icount    = 0;
  proc_now  = 0;
  my_mol    = iatm_mol_num[1];
  natm_proc = ((class->clatoms_info.natm_proc != 0)
               ? class->clatoms_info.natm_proc : natm_tot/np_forc);
  if(myid_forc==0){
    class->clatoms_info.myatm_start = 1;
    class->class_comm_forc_pkg.myatm_start = 1;
  }/*endif*/

  for(i=1;i<=natm_tot;i++){

    icount++;
    if((icount > natm_proc) && (my_mol != iatm_mol_num[i]) &&
       (proc_now != np_forc-1) ){
      recv_count_atm[(proc_now+1)] = icount-1;
      icount = 1;
      if(proc_now==myid_forc){
        class->clatoms_info.myatm_end = i-1;
        class->class_comm_forc_pkg.myatm_end = i-1;
      }/*endif*/
      proc_now++;
      if(proc_now==myid_forc){
        class->clatoms_info.myatm_start = i;
        class->class_comm_forc_pkg.myatm_start = i;
      }/*endif*/
    }/*endif*/
    my_mol           = iatm_mol_num[i];
    iatm_proc_num[i] = proc_now;

  }/*endfor*/

  recv_count_atm[(proc_now+1)] = icount;
  if(proc_now==myid_forc){
    class->clatoms_info.myatm_end = natm_tot;
    class->class_comm_forc_pkg.myatm_end = natm_tot;
  }/*endif*/

 if(proc_now+1<np_forc){
    if(myid==0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Zero atoms on the processors > %d\n",proc_now);
      printf("Try reducing the number of atoms/proc from %d\n",natm_proc);
      printf("using the min_num_atoms_per_proc keyword in gen\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }/*endif*/
    exit(1);
 }/*endif*/

/*======================================================================*/
/* I) Calculate the displacement of atoms into atm array between processors */

  displs_atm[1] = 0;
  for(i=2;i<=np_forc;i++){
   displs_atm[i] = displs_atm[(i-1)] + recv_count_atm[(i-1)];
  }/*endfor*/

#ifdef DEBUG
  for(i=0;i<=np_forc-1;i++){
   Dbx_Barrier(world);
   if(myid_forc==i){
     printf("Myatm_start and Myatmend %d %d\n",
                         class->clatoms_info.myatm_start,
                         class->clatoms_info.myatm_end);
    for(j=1;j<=np_forc;j++){
     printf("Displs and Recv %d %d %d\n",displs_atm[j],recv_count_atm[j],j);
    }/*endfor*/
    scanf("%d",&iii);
#ifdef JUNK
      for(j=1;j<=natm_tot;j++){
       printf("Proc num %d %d\n",iatm_proc_num[j],j); 
       if(j%4000==0){scanf("%d",&iii);}
      }/*endfor*/
#endif
   }/*endif*/
  }/*endfor*/
#endif

/*========================================================================*/
   } /* end routine */
/*==========================================================================*/






/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_bonded_forc(BONDED *bonded,CLASS *class)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

#include "../typ_defs/typ_mask.h"    

    int iproc1,iproc2,iproc3,iproc4,ierr,ic,i,iii; 

    int *iatm_proc_num = class->atommaps.iatm_proc_num;
    int *bond_j1_pow        = bonded->bond.j1_pow;
    int *bond_j2_pow        = bonded->bond.j2_pow;
    int *bond_jtyp_pow      = bonded->bond.jtyp_pow;
    int bond_npow           = bonded->bond.npow;

    int *onfo_j1            = bonded->onfo.j1;
    int *onfo_j2            = bonded->onfo.j2;
    int *onfo_jtyp          = bonded->onfo.jtyp;
    int onfo_num            = bonded->onfo.num;

    int bond_num_free       = bonded->bond_free.num;
    int j1_bond_free        = bonded->bond_free.j1;
    int j2_bond_free        = bonded->bond_free.j2;
    int *bond_j1_con        = bonded->bond.j1_con;
    int *bond_j2_con        = bonded->bond.j2_con;
    int *bond_jtyp_con      = bonded->bond.jtyp_con;
    int bond_ncon           = bonded->bond.ncon;
    int *bend_j1_pow        = bonded->bend.j1_pow;
    int *bend_j2_pow        = bonded->bend.j2_pow;
    int *bend_j3_pow        = bonded->bend.j3_pow;
    int *bend_jtyp_pow      = bonded->bend.jtyp_pow;
    int bend_npow           = bonded->bend.npow;
    int bend_num_free       = bonded->bend_free.num;
    int j1_bend_free        = bonded->bend_free.j1;
    int j2_bend_free        = bonded->bend_free.j2;
    int j3_bend_free        = bonded->bend_free.j3;
    int *bend_j1_con        = bonded->bend.j1_con;
    int *bend_j2_con        = bonded->bend.j2_con;
    int *bend_j3_con        = bonded->bend.j3_con;
    int *bend_jtyp_con      = bonded->bend.jtyp_con;
    int bend_ncon           = bonded->bend.ncon;
    int *tors_j1_pow        = bonded->tors.j1_pow;
    int *tors_j2_pow        = bonded->tors.j2_pow;
    int *tors_j3_pow        = bonded->tors.j3_pow;
    int *tors_j4_pow        = bonded->tors.j4_pow;
    int *tors_jtyp_pow      = bonded->tors.jtyp_pow;
    int tors_npow           = bonded->tors.npow;
    int tors_num_free       = bonded->tors_free.num;
    int *j1_tors_free       = bonded->tors_free.j1;
    int *j2_tors_free       = bonded->tors_free.j2;
    int *j3_tors_free       = bonded->tors_free.j3;
    int *j4_tors_free       = bonded->tors_free.j4;
    double *eq_tors_free    = bonded->tors_free.eq;
    int *tors_j1_con        = bonded->tors.j1_con;
    int *tors_j2_con        = bonded->tors.j2_con;
    int *tors_j3_con        = bonded->tors.j3_con;
    int *tors_j4_con        = bonded->tors.j4_con;
    int *tors_jtyp_con      = bonded->tors.jtyp_con;
    int tors_ncon           = bonded->tors.ncon;
    int *ecor_j1            = bonded->ecor.j1;
    int *ecor_j2            = bonded->ecor.j2;
    int ecor_num            = bonded->ecor.num;
    int *excl_j1_cp         = bonded->excl.j1_cp;
    int *excl_j2_cp         = bonded->excl.j2_cp;
    int excl_num_cp         = bonded->excl.num_cp;

    int bend_bnd_num         = bonded->bend_bnd.num;
    int *bend_bnd_j1          = bonded->bend_bnd.j1;
    int *bend_bnd_j2          = bonded->bend_bnd.j2;
    int *bend_bnd_j3          = bonded->bend_bnd.j3;
    int *bend_bnd_jtyp        = bonded->bend_bnd.jtyp;
    int num_21     = bonded->grp_bond_con.num_21;
    int num_23     = bonded->grp_bond_con.num_23;
    int num_33     = bonded->grp_bond_con.num_33;
    int num_33_watts     = bonded->grp_bond_watts.num_33;
    int num_43     = bonded->grp_bond_con.num_43;
    int num_46     = bonded->grp_bond_con.num_46;
    int *jtyp_21   = bonded->grp_bond_con.jtyp_21;
    int *jtyp_23   = bonded->grp_bond_con.jtyp_23;
    int *jtyp_33   = bonded->grp_bond_con.jtyp_33;
    int *jtyp_33_watts   = bonded->grp_bond_watts.jtyp_33;
    int *jtyp_43   = bonded->grp_bond_con.jtyp_43;
    int *jtyp_46   = bonded->grp_bond_con.jtyp_46;
    int *j1_21     = bonded->grp_bond_con.j1_21;
    int *j1_23     = bonded->grp_bond_con.j1_23;
    int *j1_33     = bonded->grp_bond_con.j1_33;
    int *j1_33_watts     = bonded->grp_bond_watts.j1_33;
    int *j1_43     = bonded->grp_bond_con.j1_43;
    int *j1_46     = bonded->grp_bond_con.j1_46;
    int *j2_21     = bonded->grp_bond_con.j2_21;
    int *j2_23     = bonded->grp_bond_con.j2_23;
    int *j2_33     = bonded->grp_bond_con.j2_33;
    int *j2_33_watts     = bonded->grp_bond_watts.j2_33;
    int *j2_43     = bonded->grp_bond_con.j2_43;
    int *j2_46     = bonded->grp_bond_con.j2_46;
    int *j3_23     = bonded->grp_bond_con.j3_23;
    int *j3_33     = bonded->grp_bond_con.j3_33;
    int *j3_33_watts     = bonded->grp_bond_watts.j3_33;
    int *j3_43     = bonded->grp_bond_con.j3_43;
    int *j3_46     = bonded->grp_bond_con.j3_46;
    int *j4_43     = bonded->grp_bond_con.j4_43;
    int *j4_46     = bonded->grp_bond_con.j4_46;
    int rbar_num_free       = bonded->rbar_sig_free.nfree;
    int *j1_rbar_free        = bonded->rbar_sig_free.j1;
    int *j2_rbar_free        = bonded->rbar_sig_free.j2;

    int np_forc    = class->communicate.np_forc;
    int myid_forc  = class->communicate.myid_forc;
    MPI_Comm comm_forc = class->class_comm_forc_pkg.comm;
    MPI_Comm world = class->communicate.world;
    char *type;


    type = (char *)cmalloc(MAXWORD*sizeof(char));
   
/*======================================================================*/
/* I) Check to make sure bonds do not cross processors and fill new 
      Lists      */

/*   1) Pow Bonds   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=bond_npow;i++){
     iproc1 = iatm_proc_num[bond_j1_pow[i]];
     iproc2 = iatm_proc_num[bond_j2_pow[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      bond_j1_pow[ic] = bond_j1_pow[i];
      bond_j2_pow[ic] = bond_j2_pow[i];
      bond_jtyp_pow[ic] = bond_jtyp_pow[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->bond.npow = ic;
    if(ierr > 0){
     strcpy(type,"Power Series Bond");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Onefours   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=onfo_num;i++){
     iproc1 = iatm_proc_num[onfo_j1[i]];
     iproc2 = iatm_proc_num[onfo_j2[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      onfo_j1[ic] = onfo_j1[i];
      onfo_j2[ic] = onfo_j2[i];
      onfo_jtyp[ic] = onfo_jtyp[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->onfo.num = ic;
    if(ierr > 0){
     strcpy(type,"Onefour");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Ecors   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=ecor_num;i++){
     iproc1 = iatm_proc_num[ecor_j1[i]];
     iproc2 = iatm_proc_num[ecor_j2[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      ecor_j1[ic] = ecor_j1[i];
      ecor_j2[ic] = ecor_j2[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->ecor.num = ic;
    if(ierr > 0){
     strcpy(type,"Ecors");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) MIX Coul cors   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=excl_num_cp;i++){
     iproc1 = iatm_proc_num[excl_j1_cp[i]];
     iproc2 = iatm_proc_num[excl_j2_cp[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      excl_j1_cp[ic] = excl_j1_cp[i];
      excl_j2_cp[ic] = excl_j2_cp[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->excl.num_cp = ic;
    if(ierr > 0){
     strcpy(type,"MIX Coul cors");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Free Bonds   */
    ic = 0;   
    ierr = 0;
    if(bond_num_free>0){
     iproc1 = iatm_proc_num[j1_bond_free];
     iproc2 = iatm_proc_num[j2_bond_free];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
     }/*endif*/ 
    }/*endfor*/
    bonded->bond_free.num = ic;
    if(ierr > 0){
     strcpy(type,"Free energy Bond");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Con Bonds   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=bond_ncon;i++){
     iproc1 = iatm_proc_num[bond_j1_con[i]];
     iproc2 = iatm_proc_num[bond_j2_con[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      bond_j1_con[ic] = bond_j1_con[i];
      bond_j2_con[ic] = bond_j2_con[i];
      bond_jtyp_con[ic] = bond_jtyp_con[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->bond.ncon = ic;
    Allreduce(&(ic), &(bonded->bond.ncon_max),1,MPI_INT,
                   MPI_MAX,0,comm_forc);

    if(ierr > 0){
     strcpy(type,"Constrained Bond");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Pow Bends   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=bend_npow;i++){
     iproc1 = iatm_proc_num[bend_j1_pow[i]];
     iproc2 = iatm_proc_num[bend_j2_pow[i]];
     iproc3 = iatm_proc_num[bend_j3_pow[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      bend_j1_pow[ic] = bend_j1_pow[i];
      bend_j2_pow[ic] = bend_j2_pow[i];
      bend_j3_pow[ic] = bend_j3_pow[i];
      bend_jtyp_pow[ic] = bend_jtyp_pow[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->bend.npow = ic;
    if(ierr > 0){
     strcpy(type,"Power Series Bend");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Free Bends   */
    ic = 0;   
    ierr = 0;
    if(bend_num_free>0){
     iproc1 = iatm_proc_num[j1_bend_free];
     iproc2 = iatm_proc_num[j2_bend_free];
     iproc3 = iatm_proc_num[j3_bend_free];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
     }/*endif*/ 
    }/*endfor*/
    bonded->bend_free.num = ic;
    if(ierr > 0){
     strcpy(type,"Free energy bend");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/


/*   1) Bend_bnds   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=bend_bnd_num;i++){
     iproc1 = iatm_proc_num[bend_bnd_j1[i]];
     iproc2 = iatm_proc_num[bend_bnd_j2[i]];
     iproc3 = iatm_proc_num[bend_bnd_j3[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      bend_bnd_j1[ic] = bend_bnd_j1[i];
      bend_bnd_j2[ic] = bend_bnd_j2[i];
      bend_bnd_j3[ic] = bend_bnd_j3[i];
      bend_bnd_jtyp[ic] = bend_bnd_jtyp[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->bend_bnd.num = ic;
    if(ierr > 0){
     strcpy(type,"Uri-Bradley");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Con Bends   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=bend_ncon;i++){
     iproc1 = iatm_proc_num[bend_j1_con[i]];
     iproc2 = iatm_proc_num[bend_j2_con[i]];
     iproc3 = iatm_proc_num[bend_j3_con[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      bend_j1_con[ic] = bend_j1_con[i];
      bend_j2_con[ic] = bend_j2_con[i];
      bend_j3_con[ic] = bend_j3_con[i];
      bend_jtyp_con[ic] = bend_jtyp_con[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->bend.ncon = ic;
    if(ierr > 0){
     strcpy(type,"Constrained Bend");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Constrained Torsions   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=tors_ncon;i++){
     iproc1 = iatm_proc_num[tors_j1_con[i]];
     iproc2 = iatm_proc_num[tors_j2_con[i]];
     iproc3 = iatm_proc_num[tors_j3_con[i]];
     iproc4 = iatm_proc_num[tors_j4_con[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc1 != iproc4){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc2 != iproc4){ierr++;}
     if(iproc3 != iproc4){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      tors_j1_con[ic] = tors_j1_con[i];
      tors_j2_con[ic] = tors_j2_con[i];
      tors_j3_con[ic] = tors_j3_con[i];
      tors_j4_con[ic] = tors_j4_con[i];
      tors_jtyp_con[ic] = tors_jtyp_con[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->tors.ncon = ic;
    if(ierr > 0){
     strcpy(type,"Constrained Torsion");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Power Series Torsions   */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=tors_npow;i++){
     iproc1 = iatm_proc_num[tors_j1_pow[i]];
     iproc2 = iatm_proc_num[tors_j2_pow[i]];
     iproc3 = iatm_proc_num[tors_j3_pow[i]];
     iproc4 = iatm_proc_num[tors_j4_pow[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc1 != iproc4){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc2 != iproc4){ierr++;}
     if(iproc3 != iproc4){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      tors_j1_pow[ic] = tors_j1_pow[i];
      tors_j2_pow[ic] = tors_j2_pow[i];
      tors_j3_pow[ic] = tors_j3_pow[i];
      tors_j4_pow[ic] = tors_j4_pow[i];
      tors_jtyp_pow[ic] = tors_jtyp_pow[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->tors.npow = ic;
    if(ierr > 0){
     strcpy(type,"Power Series Torsion");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Free torsions   */
    ic = 0;   
    ierr = 0;
    if(tors_num_free>0){
     for(i=1;i<=tors_num_free;i++){
      iproc1 = iatm_proc_num[j1_tors_free[i]];
      iproc2 = iatm_proc_num[j2_tors_free[i]];
      iproc3 = iatm_proc_num[j3_tors_free[i]];
      iproc4 = iatm_proc_num[j4_tors_free[i]];
      if(iproc1 != iproc2){ierr++;}
      if(iproc1 != iproc3){ierr++;}
      if(iproc1 != iproc4){ierr++;}
      if(iproc2 != iproc3){ierr++;}
      if(iproc2 != iproc4){ierr++;}
      if(iproc3 != iproc4){ierr++;}
      if(iproc1 == myid_forc){
        ic++;
        j1_tors_free[ic] = j1_tors_free[i];
        j2_tors_free[ic] = j2_tors_free[i];
        j3_tors_free[ic] = j3_tors_free[i];
        j4_tors_free[ic] = j4_tors_free[i];
        eq_tors_free[ic] = eq_tors_free[i];
      }/*endif*/ 
     }/*endfor*/
     if( (tors_num_free!=ic) && (ic>0) ){ierr++;} 
     bonded->tors_free.num = ic;
     if(ierr > 0){
       strcpy(type,"Free energy Torsion");
       molecule_decomp_err(ierr,type,myid_forc);
     }/*endif*/
    }/*endif*/


/*   1) Grp_bond 21's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_21;i++){
     iproc1 = iatm_proc_num[j1_21[i]];
     iproc2 = iatm_proc_num[j2_21[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_21[ic] = j1_21[i];
      j2_21[ic] = j2_21[i];
      jtyp_21[ic] = jtyp_21[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_con.num_21 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond 21");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Grp_bond 23's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_23;i++){
     iproc1 = iatm_proc_num[j1_23[i]];
     iproc2 = iatm_proc_num[j2_23[i]];
     iproc3 = iatm_proc_num[j3_23[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_23[ic] = j1_23[i];
      j2_23[ic] = j2_23[i];
      j3_23[ic] = j3_23[i];
      jtyp_23[ic] = jtyp_23[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_con.num_23 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond 23");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Grp_bond 33's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_33;i++){
     iproc1 = iatm_proc_num[j1_33[i]];
     iproc2 = iatm_proc_num[j2_33[i]];
     iproc3 = iatm_proc_num[j3_33[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_33[ic] = j1_33[i];
      j2_33[ic] = j2_33[i];
      j3_33[ic] = j3_33[i];
      jtyp_33[ic] = jtyp_33[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_con.num_33 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond 33");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Grp_bond watts_33's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_33_watts;i++){
     iproc1 = iatm_proc_num[j1_33_watts[i]];
     iproc2 = iatm_proc_num[j2_33_watts[i]];
     iproc3 = iatm_proc_num[j3_33_watts[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_33_watts[ic] = j1_33_watts[i];
      j2_33_watts[ic] = j2_33_watts[i];
      j3_33_watts[ic] = j3_33_watts[i];
      jtyp_33_watts[ic] = jtyp_33_watts[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_watts.num_33 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond watts 33");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Grp_bond 43's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_43;i++){
     iproc1 = iatm_proc_num[j1_43[i]];
     iproc2 = iatm_proc_num[j2_43[i]];
     iproc3 = iatm_proc_num[j3_43[i]];
     iproc4 = iatm_proc_num[j4_43[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc1 != iproc4){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc2 != iproc4){ierr++;}
     if(iproc3 != iproc4){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_43[ic] = j1_43[i];
      j2_43[ic] = j2_43[i];
      j3_43[ic] = j3_43[i];
      j4_43[ic] = j4_43[i];
      jtyp_43[ic] = jtyp_43[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_con.num_43 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond 43");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Grp_bond 46's */
    ic = 0;   
    ierr = 0;
    for(i=1;i<=num_46;i++){
     iproc1 = iatm_proc_num[j1_46[i]];
     iproc2 = iatm_proc_num[j2_46[i]];
     iproc3 = iatm_proc_num[j3_46[i]];
     iproc4 = iatm_proc_num[j4_46[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 != iproc3){ierr++;}
     if(iproc1 != iproc4){ierr++;}
     if(iproc2 != iproc3){ierr++;}
     if(iproc2 != iproc4){ierr++;}
     if(iproc3 != iproc4){ierr++;}
     if(iproc1 == myid_forc){
      ic++;
      j1_46[ic] = j1_46[i];
      j2_46[ic] = j2_46[i];
      j3_46[ic] = j3_46[i];
      j4_46[ic] = j4_46[i];
      jtyp_46[ic] = jtyp_46[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->grp_bond_con.num_46 = ic;
    if(ierr > 0){
     strcpy(type,"Group bond 46");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/

/*   1) Rbar-Sigma Bonds   */

    ic = 0;   
    ierr = 0;
    for(i=1;i<=rbar_num_free;i++){
     iproc1 = iatm_proc_num[j1_rbar_free[i]];
     iproc2 = iatm_proc_num[j2_rbar_free[i]];
     if(iproc1 != iproc2){ierr++;}
     if(iproc1 == myid_forc){
       ic++;
       j1_rbar_free[ic] = j1_rbar_free[i];
       j2_rbar_free[ic] = j2_rbar_free[i];
     }/*endif*/ 
    }/*endfor*/
    bonded->rbar_sig_free.nfree = ic;
    if(ierr > 0){
     strcpy(type,"Rbar-Sigma free energy");
     molecule_decomp_err(ierr,type,myid_forc);
    }/*endif*/


    cfree(type);

#ifdef DEBUG

  for(i=0;i<=np_forc-1;i++){
   Dbx_Barrier(world);
   if(myid_forc==i){
    printf("On processor %d there are:\n",i);
    printf("There are %d power series bonds   \n",bonded->bond.npow);
    printf("There are %d constrained bonds    \n",bonded->bond.ncon);
    printf("There are %d power series bends   \n",bonded->bend.npow);
    printf("There are %d constrained bends    \n",bonded->bend.ncon);
    printf("There are %d Urey-Bradley bends   \n",bonded->bend_bnd.num);
    printf("There are %d Total torsions\n",bonded->tors.npow);
    printf("There are %d constrained torsions \n",bonded->tors.ncon);
    printf("There are %d lj onefours          \n",bonded->onfo.num);
    printf("There are %d free energy bonds    \n",bonded->bond_free.num);
    printf("There are %d free energy bends    \n",bonded->bend_free.num);
    printf("There are %d free energy torsions \n",bonded->tors_free.num);
    printf("There are %d 21 group constraints \n",bonded->grp_bond_con.num_21);
    printf("There are %d 23 group constraints \n",bonded->grp_bond_con.num_23);
    printf("There are %d 33 group constraints \n",bonded->grp_bond_con.num_33);
    printf("There are %d 33 group Watts \n",bonded->grp_bond_watts.num_33);
    printf("There are %d 43 group constraints \n",bonded->grp_bond_con.num_43);
    printf("There are %d 46 group constraints \n",bonded->grp_bond_con.num_46);
    printf("There are %d rbar-sigma bonds     \n",bonded->rbar_sigma.nfree);
    printf("\n");

   scanf("%d",&iii);
   }/*endif*/
  }/*endfor*/
#endif


/*========================================================================*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_thermo_forc(CLASS *class)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

  int icount,proc_now,my_mol,i,j,iii,ilast;   
  int natm_proc,ip,ic,is,itherm_share_min_old;
  
/*======================================================================*/
/*               Local Pointer declarations                            */

  MPI_Comm world    = class->communicate.world;
  int myid_forc     = class->communicate.myid_forc;
  int np_forc       = class->communicate.np_forc;
  int natm_tot      = class->clatoms_info.natm_tot;
  int *iatm_proc_num = class->atommaps.iatm_proc_num;
  int *itherm_proc,*itherm_share_min;
  int *map_share,*recv_count_therm,*displs_therm,*nshare;
  int num_nhc = class->therm_info_class.num_nhc;
  int num_nhc_share = class->therm_info_class.num_nhc_share;
  int num_nhc_proc = class->therm_info_class.num_nhc_proc;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int *inhc_x = class->therm_info_class.inhc_x;
  int *inhc_y = class->therm_info_class.inhc_y;
  int *inhc_z = class->therm_info_class.inhc_z;
  int *itherm_nshare,*itherm_share_last;
  double *ditherm_nshare_i;
  int mytherm_start,mytherm_end;
  double **v_nhc;

/*======================================================================*/
/* I) Find the thermostats on this processor which are shared */             

 /* Smallest processor number sharing this thermostat */
  itherm_share_min = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* 1/0 depending on whether this thermostat is on this processor */
  itherm_proc  = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* Largest processor number sharing this thermostat */
  itherm_share_last  = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* Number of processors sharing this thermostat */
  class->therm_info_class.itherm_nshare  = 
                       (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* 1/((double) itherm_nshare) */
  class->therm_info_class.ditherm_nshare_i  = 
                       (double *)cmalloc((num_nhc+1)*sizeof(double))-1;


  itherm_nshare     =   class->therm_info_class.itherm_nshare;
  ditherm_nshare_i  = class->therm_info_class.ditherm_nshare_i;

  for(i=1;i<=num_nhc;i++){
    itherm_nshare[i] = 0;
    itherm_share_min[i] = np_forc;
    itherm_proc[i]  = 0;
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){

   if(inhc_x[i] <= num_nhc){
     if(myid_forc == iatm_proc_num[i]){
       itherm_proc[inhc_x[i]]  = 1;
     }/*endif*/
     if(itherm_nshare[inhc_x[i]]==0){
       itherm_share_min[inhc_x[i]] = iatm_proc_num[i];
       itherm_nshare[inhc_x[i]] = 1;
       itherm_share_last[inhc_x[i]] = iatm_proc_num[i];
     } else {
      if(itherm_share_last[inhc_x[i]] != iatm_proc_num[i]){
        itherm_share_min[inhc_x[i]] = MIN(itherm_share_min[inhc_x[i]],
                                        iatm_proc_num[i]);
        itherm_share_last[inhc_x[i]] = iatm_proc_num[i];
        itherm_nshare[inhc_x[i]]++;
      }/*endif*/
     }/*endif*/
   }/*endif*/

   if(inhc_y[i] <= num_nhc){
    if(myid_forc == iatm_proc_num[i]){
     itherm_proc[inhc_y[i]]  = 1;
    }/*endif*/
    if(itherm_nshare[inhc_y[i]]==0){
     itherm_share_min[inhc_y[i]] = iatm_proc_num[i];
     itherm_nshare[inhc_y[i]] = 1;
     itherm_share_last[inhc_y[i]] = iatm_proc_num[i];
    }/*endif*/
    else{
     if(itherm_share_last[inhc_y[i]] != iatm_proc_num[i]){
      itherm_share_min[inhc_y[i]] = MIN(itherm_share_min[inhc_y[i]],
                                       iatm_proc_num[i]);
      itherm_share_last[inhc_y[i]] = iatm_proc_num[i];
      itherm_nshare[inhc_y[i]]++;
     }/*endif*/
    }/*endelse*/
   }/*endif*/

   if(inhc_z[i] <= num_nhc){
    if(myid_forc == iatm_proc_num[i]){
      itherm_proc[inhc_z[i]]  = 1;
    }/*endif*/
    if(itherm_nshare[inhc_z[i]]==0){
     itherm_share_min[inhc_z[i]] = iatm_proc_num[i];
     itherm_nshare[inhc_z[i]] = 1;
     itherm_share_last[inhc_z[i]] = iatm_proc_num[i];
    }/*endif*/
    else{
     if(itherm_share_last[inhc_z[i]] != iatm_proc_num[i]){
       itherm_share_min[inhc_z[i]] = MIN(itherm_share_min[inhc_z[i]],
                                       iatm_proc_num[i]);
       itherm_share_last[inhc_z[i]] = iatm_proc_num[i];
       itherm_nshare[inhc_z[i]]++;
     }/*endif*/
    }/*endelse*/
   }/*endif*/

  }/*endfor*/

/*======================================================================*/
/* II) Count the thermostats on this processor and shared thermos  */

   num_nhc_share = 0;
   num_nhc_proc  = 0;
  for(i=1;i<=num_nhc;i++){
   if(itherm_proc[i]==1){num_nhc_proc++;}
   if(itherm_nshare[i] > 1){num_nhc_share++;}
  }/*endfor*/

  for(i=1;i<=num_nhc;i++){
   ditherm_nshare_i[i] = 1.0/((double)itherm_nshare[i]);
  }/*endfor*/

   class->therm_info_class.map_share  = 
         (int *)cmalloc((num_nhc_share+1)*sizeof(int))-1;
   class->class_comm_forc_pkg.displs_therm = 
         (int *)cmalloc(np_forc*sizeof(int))-1;
   class->class_comm_forc_pkg.recv_count_therm  = 
         (int *)cmalloc(np_forc*sizeof(int))-1;
   class->int_scr.sc_temp  = 
         (double *)cmalloc((num_nhc_share+1)*sizeof(double))-1;

   map_share = class->therm_info_class.map_share;
   displs_therm = class->class_comm_forc_pkg.displs_therm;
   recv_count_therm  = class->class_comm_forc_pkg.recv_count_therm;
   class->therm_info_class.num_nhc_share = num_nhc_share;
   class->therm_info_class.num_nhc_proc = num_nhc_proc;

/*======================================================================*/
/* III) Compress the thermostat information and make the share map   */

   ic = 0;
   is = 0;
   mytherm_start = 1;
   mytherm_end = 0;
   itherm_share_min_old = -1;
  for(i=1;i<=num_nhc;i++){
   if(itherm_proc[i] == 1){
    ic++;
    if(ic==1){
     mytherm_start = i;
     ilast = i-1;
    }/*endif*/
     if(ilast != i-1){
      printf("@@@@@@@@@@error@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1); 
     }/*endif*/
     ilast = i;
     mytherm_end = i;
   }/*endif*/
   if(itherm_nshare[i] > 1){
    is++;
    if(itherm_share_min[i] < itherm_share_min_old){
     printf("@@@@@@@@@@error@@@@@@@@@@@@@@@\n");
     printf("%d %d\n",itherm_share_min[i],itherm_share_min_old);
     fflush(stdout);
     exit(1);
    }/*endif*/
    itherm_share_min_old = itherm_share_min[i];
    recv_count_therm[(itherm_share_min[i]+1)]++;
   }/*endif*/

   if(itherm_nshare[i] > 1 && itherm_proc[i] == 1){
     map_share[is] = i;
   } else {
     map_share[is] = num_nhc + 1;
   }/*endelse*/

  }/*endfor*/
  class->therm_info_class.mytherm_start = mytherm_start;
  class->therm_info_class.mytherm_end = mytherm_end;

  displs_therm[1] = 0;
  for(i=2;i<=np_forc;i++){
    displs_therm[i] = recv_count_therm[(i-1)] + displs_therm[(i-1)];
  }/*endfor*/


#ifdef DEBUG
  for(i=0;i<=np_forc-1;i++){
   Dbx_Barrier(world);
   if(myid_forc==i){
     printf("Mytherm_start and Mythermend %d %d\n",
                         class->therm_info_class.mytherm_start,
                         class->therm_info_class.mytherm_end);
    for(j=1;j<=np_forc;j++){
     printf("Displs and Recv %d %d %d\n",displs_therm[j],
                                         recv_count_therm[j],j);
    }/*endfor*/
     printf("num_nhc_proc %d\n",num_nhc_proc);
     printf("num_nhc_share %d\n",num_nhc_share);
    for(j=1;j<=num_nhc_share;j++){
     printf("Share map %d %d\n",map_share[j],j);
    }/*endfor*/
    for(j=1;j<=num_nhc;j++){
     printf("itherm's %d %d %d %d\n",itherm_share_min[j],
                               itherm_proc[j],itherm_nshare[j],j);
    }/*endfor*/
    
    scanf("%d",&iii);
   }/*endif*/
  }/*endfor*/
#endif

  cfree(&itherm_share_min[1]);
  cfree(&itherm_proc[1]);
  cfree(&itherm_share_last[1]);

/*========================================================================*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void assign_bead_thermo_forc(CLASS *class)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */

/*======================================================================*/
/*               Local variable declarations                            */

  int icount,proc_now,my_mol,i,j,iii,ilast;   
  int natm_proc,ip,ic,is,itherm_share_min_old;
  
/*======================================================================*/
/*               Local Pointer declarations                            */

  MPI_Comm world    = class->communicate.world;
  int myid_forc     = class->communicate.myid_forc;
  int np_forc       = class->communicate.np_forc;
  int natm_tot      = class->clatoms_info.natm_tot;
  int *iatm_proc_num = class->atommaps.iatm_proc_num;
  int *itherm_proc,*itherm_share_min;
  int *map_share,*recv_count_therm,*displs_therm,*nshare;
  int num_nhc = class->therm_info_bead.num_nhc;
  int num_nhc_share = class->therm_info_bead.num_nhc_share;
  int num_nhc_proc = class->therm_info_bead.num_nhc_proc;
  int pi_beads_proc = class->clatoms_info.pi_beads_proc;
  int *inhc_x = class->therm_info_bead.inhc_x;
  int *inhc_y = class->therm_info_bead.inhc_y;
  int *inhc_z = class->therm_info_bead.inhc_z;
  int *itherm_nshare,*itherm_share_last;
  double *ditherm_nshare_i;
  int mytherm_start,mytherm_end;
  double **v_nhc;

/*======================================================================*/
/* I) Find the thermostats on this processor which are shared */             

 /* Smallest processor number sharing this thermostat */
  itherm_share_min = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* 1/0 depending on whether this thermostat is on this processor */
  itherm_proc  = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* Largest processor number sharing this thermostat */
  itherm_share_last  = (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* Number of processors sharing this thermostat */
  class->therm_info_bead.itherm_nshare  = 
                       (int *)cmalloc((num_nhc+1)*sizeof(int))-1;

 /* 1/((double) itherm_nshare) */
  class->therm_info_bead.ditherm_nshare_i  = 
                       (double *)cmalloc((num_nhc+1)*sizeof(double))-1;


  itherm_nshare     =   class->therm_info_bead.itherm_nshare;
  ditherm_nshare_i  = class->therm_info_bead.ditherm_nshare_i;

  for(i=1;i<=num_nhc;i++){
    itherm_nshare[i] = 0;
    itherm_share_min[i] = np_forc;
    itherm_proc[i]  = 0;
  }/*endfor*/

  for(i=1;i<=natm_tot;i++){

   if(inhc_x[i] <= num_nhc){
     if(myid_forc == iatm_proc_num[i]){
       itherm_proc[inhc_x[i]]  = 1;
     }/*endif*/
     if(itherm_nshare[inhc_x[i]]==0){
       itherm_share_min[inhc_x[i]] = iatm_proc_num[i];
       itherm_nshare[inhc_x[i]] = 1;
       itherm_share_last[inhc_x[i]] = iatm_proc_num[i];
     } else {
      if(itherm_share_last[inhc_x[i]] != iatm_proc_num[i]){
        itherm_share_min[inhc_x[i]] = MIN(itherm_share_min[inhc_x[i]],
                                        iatm_proc_num[i]);
        itherm_share_last[inhc_x[i]] = iatm_proc_num[i];
        itherm_nshare[inhc_x[i]]++;
      }/*endif*/
     }/*endif*/
   }/*endif*/

   if(inhc_y[i] <= num_nhc){
    if(myid_forc == iatm_proc_num[i]){
     itherm_proc[inhc_y[i]]  = 1;
    }/*endif*/
    if(itherm_nshare[inhc_y[i]]==0){
     itherm_share_min[inhc_y[i]] = iatm_proc_num[i];
     itherm_nshare[inhc_y[i]] = 1;
     itherm_share_last[inhc_y[i]] = iatm_proc_num[i];
    }/*endif*/
    else{
     if(itherm_share_last[inhc_y[i]] != iatm_proc_num[i]){
      itherm_share_min[inhc_y[i]] = MIN(itherm_share_min[inhc_y[i]],
                                       iatm_proc_num[i]);
      itherm_share_last[inhc_y[i]] = iatm_proc_num[i];
      itherm_nshare[inhc_y[i]]++;
     }/*endif*/
    }/*endelse*/
   }/*endif*/

   if(inhc_z[i] <= num_nhc){
    if(myid_forc == iatm_proc_num[i]){
      itherm_proc[inhc_z[i]]  = 1;
    }/*endif*/
    if(itherm_nshare[inhc_z[i]]==0){
     itherm_share_min[inhc_z[i]] = iatm_proc_num[i];
     itherm_nshare[inhc_z[i]] = 1;
     itherm_share_last[inhc_z[i]] = iatm_proc_num[i];
    }/*endif*/
    else{
     if(itherm_share_last[inhc_z[i]] != iatm_proc_num[i]){
       itherm_share_min[inhc_z[i]] = MIN(itherm_share_min[inhc_z[i]],
                                       iatm_proc_num[i]);
       itherm_share_last[inhc_z[i]] = iatm_proc_num[i];
       itherm_nshare[inhc_z[i]]++;
     }/*endif*/
    }/*endelse*/
   }/*endif*/

  }/*endfor*/

/*======================================================================*/
/* II) Count the thermostats on this processor and shared thermos  */

   num_nhc_share = 0;
   num_nhc_proc  = 0;
  for(i=1;i<=num_nhc;i++){
   if(itherm_proc[i]==1){num_nhc_proc++;}
   if(itherm_nshare[i] > 1){num_nhc_share++;}
  }/*endfor*/

  for(i=1;i<=num_nhc;i++){
   ditherm_nshare_i[i] = 1.0/((double)itherm_nshare[i]);
  }/*endfor*/

   class->therm_info_bead.map_share  = 
         (int *)cmalloc((num_nhc_share+1)*sizeof(int))-1;
   class->class_comm_forc_pkg.displs_therm = 
         (int *)cmalloc(np_forc*sizeof(int))-1;
   class->class_comm_forc_pkg.recv_count_therm  = 
         (int *)cmalloc(np_forc*sizeof(int))-1;
   class->int_scr.sc_temp  = 
         (double *)cmalloc((num_nhc_share+1)*sizeof(double))-1;

   map_share = class->therm_info_bead.map_share;
   displs_therm = class->class_comm_forc_pkg.displs_therm;
   recv_count_therm  = class->class_comm_forc_pkg.recv_count_therm;
   class->therm_info_bead.num_nhc_share = num_nhc_share;
   class->therm_info_bead.num_nhc_proc = num_nhc_proc;

/*======================================================================*/
/* III) Compress the thermostat information and make the share map   */

   ic = 0;
   is = 0;
   mytherm_start = 1;
   mytherm_end = 0;
   itherm_share_min_old = -1;
  for(i=1;i<=num_nhc;i++){
   if(itherm_proc[i] == 1){
    ic++;
    if(ic==1){
     mytherm_start = i;
     ilast = i-1;
    }/*endif*/
     if(ilast != i-1){
      printf("@@@@@@@@@@error@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(1); 
     }/*endif*/
     ilast = i;
     mytherm_end = i;
   }/*endif*/
   if(itherm_nshare[i] > 1){
    is++;
    if(itherm_share_min[i] < itherm_share_min_old){
     printf("@@@@@@@@@@error@@@@@@@@@@@@@@@\n");
     printf("%d %d\n",itherm_share_min[i],itherm_share_min_old);
     fflush(stdout);
     exit(1);
    }/*endif*/
    itherm_share_min_old = itherm_share_min[i];
    recv_count_therm[(itherm_share_min[i]+1)]++;
   }/*endif*/

   if(itherm_nshare[i] > 1 && itherm_proc[i] == 1){
     map_share[is] = i;
   } else {
     map_share[is] = num_nhc + 1;
   }/*endelse*/

  }/*endfor*/

  class->therm_info_bead.mytherm_start = mytherm_start;
  class->therm_info_bead.mytherm_end = mytherm_end;

  displs_therm[1] = 0;
  for(i=2;i<=np_forc;i++){
    displs_therm[i] = recv_count_therm[(i-1)] + displs_therm[(i-1)];
  }/*endfor*/


#ifdef DEBUG
  for(i=0;i<=np_forc-1;i++){
   Dbx_Barrier(world);
   if(myid_forc==i){
     printf("Mytherm_start and Mythermend %d %d\n",
                         class->therm_info_bead.mytherm_start,
                         class->therm_info_bead.mytherm_end);
    for(j=1;j<=np_forc;j++){
     printf("Displs and Recv %d %d %d\n",displs_therm[j],
                                         recv_count_therm[j],j);
    }/*endfor*/
     printf("num_nhc_proc %d\n",num_nhc_proc);
     printf("num_nhc_share %d\n",num_nhc_share);
    for(j=1;j<=num_nhc_share;j++){
     printf("Share map %d %d\n",map_share[j],j);
    }/*endfor*/
    for(j=1;j<=num_nhc;j++){
     printf("itherm's %d %d %d %d\n",itherm_share_min[j],
                               itherm_proc[j],itherm_nshare[j],j);
    }/*endfor*/
    
    scanf("%d",&iii);
   }/*endif*/
  }/*endfor*/
#endif

  cfree(&itherm_share_min[1]);
  cfree(&itherm_proc[1]);
  cfree(&itherm_share_last[1]);

/*========================================================================*/
} /* end routine */
/*==========================================================================*/



/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void molecule_decomp_err(int ierr,char *type,int myid_forc)

/*======================================================================*/
/*                Begin Routine */
{   /*begin routine */
/*======================================================================*/
/*               Local variable declarations                            */

     if(myid_forc==0){
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       printf("Error in molecule decomp. \n");
       printf("%d atoms found spread across processors in %s\n",ierr,type);
       printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       fflush(stdout);
     }/*endif*/
       exit(1);

/*========================================================================*/
} /* end routine */
/*==========================================================================*/










