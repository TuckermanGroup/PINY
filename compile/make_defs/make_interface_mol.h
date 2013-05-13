#==================================================================
#            INTERFACE_MOL_FILES
#==================================================================


#==================================================================
MOL_PARMS = $(CODE)/interface/mol_params
DMOL_PARMS = $(DCODE)/interface/mol_params
#==================================================================


#==================================================================
control_mol_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MOL_ENT) $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/control_mol_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/control_mol_params.c

#----------------------------------------------------------------
control_set_mol_params.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_GEN) $(TYP_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/control_set_mol_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/control_set_mol_params.c

#----------------------------------------------------------------
set_base_file_params.o :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(MOL_LOC) \
                         $(MOL_PARMS)/set_params/set_base_file_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_base_file_params.c

#----------------------------------------------------------------
set_surf_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(TYP_PAR) $(TYP_BND) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_surf_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_surf_params.c

#----------------------------------------------------------------
set_free_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_CP) \
                         $(TYP_PAR) $(TYP_BND) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_free_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_free_params.c

#----------------------------------------------------------------
set_mol_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_CLASS) \
                         $(TYP_BND) $(TYP_CP) \
                         $(MOL_LOC) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_dict.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_dict.c

#----------------------------------------------------------------
set_mol_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_PAR) $(TYP_BND) \
                         $(TYP_CP) $(TYP_CLASS) \
                         $(MOL_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_params.c
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_params.c

#==================================================================

