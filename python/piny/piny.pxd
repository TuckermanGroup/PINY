"""PINY declarations relevant to the interface.

At the moment, only things that are needed are included.
More can be added by hand or using one of the pxd generation tools:

http://wiki.cython.org/AutoPxd

"""


ctypedef char NAME[50]


cdef extern from 'pentium_nopar/standard_include.h':

    pass


cdef extern from 'typ_defs/typedefs_class.h':

    ctypedef struct CLATOMS_INFO:
        int natm_tot
        int pi_beads
        int pi_beads_proc
        int pi_beads_proc_st
        int pi_beads_proc_end
        int nfree
        int nfree_pimd
        int nchrg
        int pi_beads_full_ter
        int pi_beads_res_ter
        int pi_beads_full_tra
        int pi_beads_res_tra
        int cg_reset_flag
        int natm_mall
        int nchrg_mall
        int hess_calc
        int natm_proc
        int nab_initio
        int myatm_start, myatm_end
        double gamma_adb
        double wght_pimd
        double pi_temperature
        double pi_beads_full_ter_wght
        double pi_beads_res_ter_wght
        double pi_beads_full_tra_wght
        double pi_beads_res_tra_wght
        double rcut_spread
        double mass_sc_fact
        int *ip_lab
        int *ichrg
        int *cp_vlnc_up, *cp_vlnc_dn
        int *cp_atm_flag
        double *mass
        double *q
        double *alp_pol
        double *b_neut
        double *roll_sc
        double *text_atm
        double *text_mol
        double *xold, *yold, *zold
        double *prekf
        double *xmod, *ymod, *zmod
        double alp_ewd

    ctypedef struct CLATOMS_TRAN:
        pass

    ctypedef struct CLATOMS_POS:
        double *mass
        double *x, *y, *z
        double *vx, *vy, *vz
        double *fx, *fy, *fz
        double *fxt, *fyt, *fzt
        double *fxm, *fym, *fzm
        double *hess_xx, *hess_xy, *hess_xz, *hess_yy, *hess_yz, *hess_zz

    ctypedef struct GHOST_ATOMS:
        pass

    ctypedef struct ATOMMAPS:
        int nmol_typ
        int nres_typ
        int natm_typ
        int nfreeze
        int nfreeze_mall
        int natm_typ_mall
        int nres_typ_max
        int nres_max
        int nres_tot
        int nres_sum
        int pimd_freez_typ
        NAME *atm_typ
        NAME *res_typ
        NAME *mol_typ
        int *nmol_jmol_typ
        int *nres_1mol_jmol_typ
        int *jatm_jmol_typ_strt
        int *jres_jmol_typ_strt
        int *ires_typ_jres_jmol_typ
        int *jatm_jres_1mol_jmol_typ_strt
        int *natm_1mol_jmol_typ
        int *natm_jres_jmol_typ
        int *nfree_1mol_jmol_typ
        int *nfree_jres_jmol_typ
        int *icons_jmol_typ
        int *icons_jres_jmol_typ
        int *natm_jmol_typ
        int *iatm_mol_typ
        int *iatm_res_typ
        int *iatm_atm_typ
        int *iatm_atm_typ_nl
        int *iatm_atm_typ_nl_rev
        int *imap_atm_typ_nl
        int *imap_atm_typ_nl_gh
        int *iatm_mol_num
        int *iatm_proc_num
        int *iatm_res_num
        int *ighost_flag
        int *freeze_flag
        int *freeze_map
        int *atom_label
        int *cp_atm_lst

    ctypedef struct THERM_INFO:
        pass

    ctypedef struct THERM_POS:
        pass

    ctypedef struct VEL_SAMP_CLASS:
        pass

    ctypedef struct ENERGY_CTRL:
        int pme_on
        int int_res_ter,int_res_tra
        int isep_vvdw
        int iswit_vdw
        int block_std_on, block_con_on, nblock_min
        int iget_pe_real_inter_freq
        int iget_pe_real_inter, iget_pv_real_inter
        int itime
        int iget_full_inter, iget_res_inter
        int iget_full_intra, iget_res_intra
        int iget_full_pimd, iget_res_pimd

    ctypedef struct NBR_LIST:
        pass

    ctypedef struct INTERACT:
        pass

    ctypedef struct INT_SCR:
        pass

    ctypedef struct EWD_SCR:
        pass

    ctypedef struct FOR_SCR:
        pass

    ctypedef struct PART_MESH:
        pass

    ctypedef struct SURFACE:
        pass

    ctypedef struct COMMUNICATE:
        pass

    ctypedef struct CLASS_COMM_FORC_PKG:
        pass

    ctypedef struct CLASS:
        CLATOMS_INFO clatoms_info
        CLATOMS_TRAN clatoms_tran
        CLATOMS_POS *clatoms_pos
        GHOST_ATOMS ghost_atoms
        ATOMMAPS atommaps
        THERM_INFO therm_info_class
        THERM_INFO therm_info_bead
        THERM_POS therm_class
        THERM_POS *therm_bead
        VEL_SAMP_CLASS vel_samp_class
        ENERGY_CTRL energy_ctrl
        NBR_LIST nbr_list
        INTERACT interact
        INT_SCR int_scr
        EWD_SCR ewd_scr
        FOR_SCR for_scr
        double tot_memory
        PART_MESH part_mesh
        SURFACE surface
        COMMUNICATE communicate
        CLASS_COMM_FORC_PKG class_comm_forc_pkg


cdef extern from 'typ_defs/typedefs_bnd.h':

    ctypedef struct BONDED:
        pass


cdef extern from 'typ_defs/typedefs_gen.h':

    ctypedef struct SIMOPTS:
        pass

    ctypedef struct MINOPTS:
        pass

    ctypedef struct ENSOPTS:
        pass

    ctypedef struct FILENAMES:
        pass

    ctypedef struct STATEPOINT:
        pass

    ctypedef struct CELL:
        pass

    ctypedef struct TIMEINFO:
        pass

    ctypedef struct BARO:
        pass

    ctypedef struct PAR_RAHMAN:
        pass

    ctypedef struct PTENS:
        pass

    ctypedef struct STAT_AVG:
        int nexclt
        int iter_shake, iter_ratl
        int iswit_vdw
        int itime_update, itime_update_w
        int iter_shake_cp, iter_ratl_cp
        int write_cp_atm_flag
        double vintert, aivintert, avintert
        double vintrat, aivintrat, avintrat
        double vsurft
        double vbondt, vbendt, vtorst, vonfot
        double vbend_bndt
        double vbend_bnd_bond, vbend_bnd_bend
        double vbondt_watts, vbendt_watts, vtot_watts
        double vreal, vrecip
        double vvdw, vcoul
        double vlong
        double vbond_free, vbend_free, vtors_free
        double vbar_free
        double kinet, aikinet, akinet
        double kinet_v, aikinet_v, akinet_v
        double vol, aivol, avol
        double kinet_nhc, aikinet_nhc, akinet_nhc
        double kinet_nhc_bead
        double aikinet_nhc_bead
        double akinet_nhc_bead
        double vpot_v
        double vpotnhc
        double aiter_shake, aiter_ratl
        double iter_23, iter_33, iter_46, iter_43, iter_21
        double aiter_23, aiter_33, aiter_46
        double aiter_21, aiter_43
        double iter_23r, iter_33r, iter_46r, iter_43r, iter_21r
        double aiter_23r, aiter_33r, aiter_46r
        double aiter_21r, aiter_43r
        double acella, acellb, acellc
        double aicella, aicellb, aicellc
        double acellab, acellbc, acellac
        double aicellab,aicellbc,aicellac
        double apress, aipress
        double press_inter, press_intra
        double apress_inter, aipress_inter
        double apress_intra, aipress_intra
        double press_kin
        double apress_kin, aipress_kin
        double econv0, econv
        double cp_kconv0, cp_kconv
        double cpu1, cpu2, acpu, cpu_now
        double updates, updates_w
        double kinet_cp_up, aikinet_cp_up, akinet_cp_up
        double kinet_cp_dn, aikinet_cp_dn, akinet_cp_dn
        double kinet_cp, aikinet_cp, akinet_cp
        double kinet_nhc_cp, aikinet_nhc_cp, akinet_nhc_cp
        double vpotnhc_cp
        double cp_ehart, aicp_ehart, acp_ehart
        double cp_eext, aicp_eext, acp_eext
        double cp_exc, aicp_exc, acp_exc
        double cp_muxc
        double cp_eke, aicp_eke, acp_eke
        double cp_enl, aicp_enl, acp_enl
        double aiter_shake_cp, aiter_ratl_cp
        double maxfc, maxf
        double pi_ke_prim, pi_ke_vir
        double api_ke_prim, api_ke_vir
        double aipi_ke_prim, aipi_ke_vir
        double kin_harm,akin_harm, aikin_harm
        double *apten, *aipten, *apten_out
        double fc_mag_up, fc_mag_dn
        double fc_max_up, fc_max_dn
        double fatm_max, fatm_mag
        double count_diag_srot
        double econv_now

    ctypedef struct EWALD:
        pass

    ctypedef struct PARA_FFT_PKG3D:
        pass

    ctypedef struct GENERAL_DATA:
        SIMOPTS simopts
        MINOPTS minopts
        ENSOPTS ensopts
        FILENAMES filenames
        STATEPOINT statepoint
        CELL cell
        TIMEINFO timeinfo
        BARO baro
        PAR_RAHMAN par_rahman
        PTENS ptens
        STAT_AVG stat_avg
        EWALD ewald
        PARA_FFT_PKG3D pme_fft_pkg
        PARA_FFT_PKG3D pme_res_fft_pkg
        int error_check_on


cdef extern from 'typ_defs/typedefs_cp.h':

    ctypedef struct CP:
        pass


cdef extern from 'typ_defs/typedefs_stat.h':

    ctypedef struct ANALYSIS:
        pass


cdef extern from 'proto_defs/proto_parse_entry.h':
    void parse(CLASS *, BONDED *, GENERAL_DATA *, CP *, ANALYSIS *, char *)


cdef extern from 'proto_defs/proto_main_local.h':
    void prelim_md(CLASS *, BONDED *, GENERAL_DATA *)


cdef extern from 'proto_defs/proto_energy_ctrl_entry.h':
    void energy_control(CLASS *, BONDED *, GENERAL_DATA *) nogil


cdef extern from 'proto_defs/proto_intra_con_entry.h':
    void get_ghost_pos(CLATOMS_INFO *, CLATOMS_POS *clatoms_pos, GHOST_ATOMS *)


cdef extern from 'proto_defs/proto_main_entry.h':
    void Init_PINY(int, char **, CLASS*, GENERAL_DATA*)
