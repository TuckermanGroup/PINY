#==================================================================
#               INTERFACE_INTRA_FILES
#==================================================================


#==================================================================
SET_PARM  = $(CODE)/interface/intra_params/set_params
DSET_PARM = $(DCODE)/interface/intra_params/set_params
#==================================================================


#==================================================================
close_intra_params.o :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/close_intra_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/close_intra_params.c

#------------------------------------------------------------------
control_intra_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) $(PIMD_LOC) \
                         $(MATH) \
                         $(CODE)/interface/intra_params/control_intra_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/control_intra_params.c

#------------------------------------------------------------------
control_res_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/control_res_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/control_res_params.c

#------------------------------------------------------------------
fetch_residue.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_residue.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_residue.c

#------------------------------------------------------------------
fetch_resbond_prm.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) \
                         $(CODE)/interface/intra_params/fetch_resbond_prm.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_resbond_prm.c

#------------------------------------------------------------------
fetch_free_energy_index.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) \
                       $(CODE)/interface/intra_params/fetch_free_energy_index.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_free_energy_index.c

#------------------------------------------------------------------
fetch_freeze.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_freeze.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_freeze.c

#------------------------------------------------------------------
init_intra_params.o  :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/init_intra_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/init_intra_params.c

#------------------------------------------------------------------
manipulate_res_bonds.o : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/manipulate_res_bonds.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/manipulate_res_bonds.c

#------------------------------------------------------------------
replicate_mol.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/replicate_mol.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/replicate_mol.c
#------------------------------------------------------------------
residue_bond.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_BND) $(TYP_CLASS) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/residue_bond.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/residue_bond.c

#==================================================================
set_atm_mask.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_mask.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_mask.c

#------------------------------------------------------------------
set_atm_morph.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_morph.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_morph.c

#------------------------------------------------------------------
set_atm_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_params.c

#------------------------------------------------------------------
set_bend_bnd_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_bnd_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bend_bnd_params.c

#------------------------------------------------------------------
set_bend_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bend_params.c

#------------------------------------------------------------------
set_bond_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bond_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bond_params.c

#------------------------------------------------------------------
set_intra_dict.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict.c

#------------------------------------------------------------------
set_intra_dict_pot.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict_pot.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict_pot.c

#------------------------------------------------------------------
set_intra_potent.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/set_intra_potent.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/set_intra_potent.c
#------------------------------------------------------------------
intra_coefs.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(CODE)/interface/intra_params/intra_coefs.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/intra_coefs.c

#------------------------------------------------------------------
set_mol_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_mol_name_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_mol_name_params.c

#------------------------------------------------------------------
set_onfo_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_onfo_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_onfo_params.c

#------------------------------------------------------------------
set_res_bond_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_bond_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_bond_params.c

#------------------------------------------------------------------
set_res_def_params.o   : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_def_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_def_params.c

#------------------------------------------------------------------
set_res_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_name_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_name_params.c

#------------------------------------------------------------------
set_res_morph_params.o : $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_morph_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_morph_params.c

#------------------------------------------------------------------
set_grp_con_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_grp_con_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_grp_con_params.c

#------------------------------------------------------------------
set_tors_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(SET_PARM)/set_tors_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_tors_params.c

#------------------------------------------------------------------
fetch_hydrog_mass.o :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_PAR) $(TYP_BND) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_hydrog_mass.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_hydrog_mass.c

#==================================================================
