/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
/*                                                                          */
/*                         PI_MD:                                           */
/*             The future of simulation technology                          */
/*             ------------------------------------                         */
/*                     Module: mall_class.c                                */
/*                                                                          */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_par.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../proto_defs/proto_communicate_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_parse_entry.h"

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

void mall_class(CLASS *class,GENERAL_DATA *general_data,
                CLASS_PARSE *class_parse,int pimd_on)

/*========================================================================*/
/*     Begin routine                                                      */
    {/*begin routine*/
/*========================================================================*/
/* Local variables */

   int i, iii;

/*========================================================================*/
/* Local Pointer Sizes                                                    */

  int pi_beads        = class->clatoms_info.pi_beads;
  int natm_mall       = class->clatoms_info.natm_mall;
  int nghost_mall     = class->ghost_atoms.nghost_mall;
  int ncomp_mall      = class->ghost_atoms.ncomp_mall; 
  int nchrg_mall      = class->clatoms_info.nchrg_mall;
  int nres_typ_max    = class->atommaps.nres_typ_max;
  int natm_typ_mall   = class->atommaps.natm_typ_mall;
  int nfreeze_mall    = class->atommaps.nfreeze_mall; 
  int nres_tot        = class->atommaps.nres_tot;
  int nmol_typ        = class->atommaps.nmol_typ;
  int nres_sum        = class->atommaps.nres_sum;
  int nmol_typ_mall   = class->atommaps.nmol_typ;
  int ninter_mall     = class->interact.ninter_mall;
  int nsplin_mall_tot = class->interact.nsplin_mall_tot;

  size_t nine,four,three;
  nine                = (size_t)9;
  four                = (size_t)4;
  three               = (size_t)3;

#ifdef DEBUG_MALLOC
  printf("mall's %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
           natm_mall,nghost_mall,ncomp_mall,nchrg_mall,nres_typ_max,
           natm_typ_mall,nfreeze_mall,nres_tot,nmol_typ,nres_sum,
           nmol_typ_mall,nktot_mall,ninterp_mall,ninter_mall,
           nsplin_mall_tot);
  scanf("%d",&iii); 
#endif

/*========================================================================*/
/*     class->cell.cell_data                                             */
/*========================================================================*/

 general_data->cell.hmat        = (double *)cmalloc(nine*sizeof(double))-1;
 general_data->cell.hmati       = (double *)cmalloc(nine*sizeof(double))-1;
 general_data->cell.hmat_ewd    = (double *)cmalloc(nine*sizeof(double))-1;
 general_data->cell.hmat_ewd_cp = (double *)cmalloc(nine*sizeof(double))-1;

 general_data->cell.hmat_cp     = (double *)cmalloc(nine*sizeof(double))-1;
 general_data->cell.hmati_cp    = (double *)cmalloc(nine*sizeof(double))-1;
 general_data->cell.cp_box_center = (double *)cmalloc(three*sizeof(double))-1;
 general_data->cell.cp_box_center_rel
                                =(double *)cmalloc(three*sizeof(double))-1;

 general_data->cell.cp_vbox_center = (double *)cmalloc(three*sizeof(double))-1;
 general_data->cell.cp_fbox_center = (double *)cmalloc(three*sizeof(double))-1;

/*========================================================================*/
/*     class->clatoms_data                                               */
/*========================================================================*/

  class->clatoms_info.text_mol  = 
                        (double *) cmalloc(nmol_typ_mall*sizeof(double ))-1;

  class->clatoms_info.xold    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.yold    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.zold    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.roll_sc = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.mass    = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.q       = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.alp_pol = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.b_neut  = (double *)cmalloc(natm_mall*sizeof(double))-1;
  class->clatoms_info.text_atm= (double *)cmalloc(natm_mall*sizeof(double))-1;

  if(pi_beads>1||pimd_on==1){
   class->clatoms_info.xmod  = (double *)cmalloc(natm_mall*sizeof(double))-1;
   class->clatoms_info.ymod  = (double *)cmalloc(natm_mall*sizeof(double))-1;
   class->clatoms_info.zmod  = (double *)cmalloc(natm_mall*sizeof(double))-1;
   class->clatoms_info.prekf = (double *)cmalloc(natm_mall*sizeof(double))-1;
  }/*endif*/


/*========================================================================*/
/*     class->clatoms_list                                               */
/*========================================================================*/

  class->clatoms_info.ichrg  = (int *)cmalloc(natm_mall*sizeof(int))-1;
  class->clatoms_info.ip_lab = (int *)cmalloc(pi_beads*sizeof(int))-1;
  class->clatoms_info.cp_vlnc_up  = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->clatoms_info.cp_vlnc_dn  = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->clatoms_info.cp_atm_flag = (int *) cmalloc(natm_mall*sizeof(int))-1;

/*========================================================================*/
/*     class->ghost_atoms.ghost_atoms_list                               */
/*========================================================================*/

  class->ghost_atoms.ighost_map = (int *) cmalloc(nghost_mall*sizeof(int))-1;
  class->ghost_atoms.natm_comp  = (int *) cmalloc(nghost_mall*sizeof(int))-1;
  class->ghost_atoms.iatm_comp  = cmall_int_mat(1,ncomp_mall,1,nghost_mall);

/*========================================================================*/
/*     class->ghost_atoms.ghost_atoms_data                               */
/*========================================================================*/

  class->ghost_atoms.coef   = cmall_mat(1,ncomp_mall,1,nghost_mall);

/*========================================================================*/
/*     class->atommaps.atommaps_list                                     */
/*========================================================================*/

  class->atommaps.atm_typ    = (NAME *)cmalloc(natm_typ_mall*sizeof(NAME))-1;
  class->atommaps.res_typ    = (NAME *)cmalloc(nres_typ_max*sizeof(NAME))-1;
  class->atommaps.freeze_map = (int *) cmalloc(nfreeze_mall*sizeof(int))-1;

  class->atommaps.jatm_jres_1mol_jmol_typ_strt = 
                  (int *) cmalloc(nres_tot*sizeof(int))-1;
  class->atommaps.ires_typ_jres_jmol_typ       = 
                  (int *) cmalloc(nres_tot*sizeof(int))-1;
  class->atommaps.natm_jres_jmol_typ           = 
                  (int *) cmalloc(nres_tot*sizeof(int))-1;
  class->atommaps.nfree_jres_jmol_typ          = 
                  (int *) cmalloc(nres_tot*sizeof(int))-1;
  class->atommaps.icons_jres_jmol_typ          = 
                  (int *) cmalloc(nres_tot*sizeof(int))-1;   

  class->atommaps.ighost_flag  = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.freeze_flag  = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.atom_label   = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.iatm_mol_typ = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.iatm_res_typ = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.iatm_atm_typ = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.iatm_mol_num = (int *) cmalloc(natm_mall*sizeof(int))-1;
  class->atommaps.iatm_res_num = (int *) cmalloc(natm_mall*sizeof(int))-1;

  class->atommaps.mol_typ  = (NAME *)cmalloc(nmol_typ_mall*sizeof(NAME))-1;
  class->atommaps.natm_1mol_jmol_typ =
                  (int *) cmalloc(nmol_typ_mall*sizeof(int))-1;
  class->atommaps.jatm_jmol_typ_strt =
                  (int *) cmalloc(nmol_typ_mall*sizeof(int))-1;
  class->atommaps.nfree_1mol_jmol_typ =
                  (int *) cmalloc(nmol_typ_mall*sizeof(int))-1;
  class->atommaps.icons_jmol_typ      = 
                  (int *) cmalloc(nmol_typ_mall*sizeof(int))-1;

  class->atommaps.nmol_jmol_typ      =
                  (int *) cmalloc(nmol_typ*sizeof(int))-1; 
  class->atommaps.nres_1mol_jmol_typ =
                  (int *) cmalloc(nmol_typ*sizeof(int))-1; 
  class->atommaps.jres_jmol_typ_strt =
                  (int *) cmalloc(nmol_typ*sizeof(int))-1; 

/*========================================================================*/
/*     class->baro                                                       */
/*========================================================================*/

  general_data->baro.hmato   =  (double *)cmalloc(nine*sizeof(double))-1;

/*========================================================================*/
/*     class->par_rahman                                                 */
/*========================================================================*/

  general_data->par_rahman.vgmat     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vgmat_g   = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vgmat_glob= (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.fgmat_p   = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.fgmat_v   = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.roll_mtv  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.roll_mtx  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.roll_mtvv = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vtemps    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.veigv     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.veig      = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vexpdt    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vsindt    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vtempx    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.vtempv    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.hmat_t    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.hmato     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->par_rahman.fv1       = (double *)cmalloc(three*sizeof(double))-1;
  general_data->par_rahman.fv2       = (double *)cmalloc(three*sizeof(double))-1;

/*========================================================================*/
/*     class->ptens                                                      */
/*========================================================================*/

  general_data->ptens.tvten         = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten         = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_tot     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_tmp     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_tmp_res = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_t1  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_t2  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_t3  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_t4  = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_a   = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_old = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.count_inc     = (double *)cmalloc(four*sizeof(double))-1;
  general_data->ptens.pvten_inc_glob= (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_whol= (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_std = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_std = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_inc_std = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_mode    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->ptens.pvten_mode_tot= (double *)cmalloc(nine*sizeof(double))-1;

/*========================================================================*/
/*     class->stat_avg                                                   */
/*========================================================================*/

  general_data->stat_avg.apten     = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->stat_avg.aipten    = (double *)cmalloc(nine*sizeof(double))-1;
  general_data->stat_avg.apten_out = (double *)cmalloc(nine*sizeof(double))-1;

/*========================================================================*/
/*     class_parse->class_parse_data                                          */
/*========================================================================*/

  class_parse->tau_nhc_mol=
              (double *)cmalloc(class->atommaps.nmol_typ*sizeof(double))-1;
  class_parse->text_nhc_mol=
              (double *)cmalloc(class->atommaps.nmol_typ*sizeof(double))-1;
  class_parse->mol_hydrog_mass_val=
              (double *)cmalloc(class->atommaps.nmol_typ*sizeof(double))-1;

/*========================================================================*/
/*     class_parse->class_parse_list                                          */
/*========================================================================*/

  class_parse->imol_nhc_opt=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;  
  class_parse->ionfo_opt=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;
  class_parse->ires_bond_conv=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;
  class_parse->mol_freeze_opt=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;  
  class_parse->mol_hydrog_mass_opt=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;
  class_parse->mol_hydrog_con_opt=
             (int *)cmalloc(class->atommaps.nmol_typ*sizeof(int))-1;

/*------------------------------------------------------------------------*/
} /*end routine*/ 
/*==========================================================================*/













