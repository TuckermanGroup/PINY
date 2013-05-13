#======================================================================
#            ENERGY_INTRA_FILES
#======================================================================


#======================================================================
bond.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/bond.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/bond.c

#------------------------------------------------------------------
bond_both.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/bond_both.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/bond_both.c

#------------------------------------------------------------------
bond_watts.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/bond_watts.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/bond_watts.c

#------------------------------------------------------------------
bend.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/bend.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/bend.c

#------------------------------------------------------------------
bend_bnd.o     :         $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/bend_bnd.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/bend_bnd.c

#------------------------------------------------------------------
tors.o     :             $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/tors.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/tors.c

#------------------------------------------------------------------
onefour.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/onefour.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/onefour.c

#------------------------------------------------------------------
ecorr.o     :            $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/ecorr.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/ecorr.c

#------------------------------------------------------------------
rbar_sigma.o  :          $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(INTRA_ENT) $(ENR_CTRL_ENT) \
                         $(CODE)/energy/intra/rbar_sigma.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/intra/rbar_sigma.c

#======================================================================





