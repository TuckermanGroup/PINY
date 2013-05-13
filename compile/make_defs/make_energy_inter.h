#======================================================================
#               ENERGY_INTER_FILES 
#======================================================================


#======================================================================
energy_control.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(MATH) \
                         $(CODE)/energy/control/energy_control.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control.c

#------------------------------------------------------------------
energy_control_initial.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(COMM_WRAP)    $(MATH) \
                         $(CODE)/energy/control/energy_control_initial.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_initial.c

#------------------------------------------------------------------
energy_control_intra.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)     $(MATH) \
                         $(CODE)/energy/control/energy_control_intra.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_intra.c

#------------------------------------------------------------------
energy_control_surf.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_LOC) $(SURF_ENT) \
                         $(CODE)/energy/control/energy_control_surf.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_surf.c

#------------------------------------------------------------------
energy_control_inter_real.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(COMM_WRAP) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)  $(MATH) \
                         $(CODE)/energy/control/energy_control_inter_real.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_inter_real.c

#------------------------------------------------------------------
energy_control_inter_recip.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT) \
                         $(REC_ENT)  $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(PIMD_ENT)  $(MATH) \
                         $(CODE)/energy/control/energy_control_inter_recip.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_inter_recip.c

#------------------------------------------------------------------
energy_control_final.o : $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(INTRA_ENT) $(REAL_ENT)  \
                         $(REC_ENT) $(ENR_CTRL_LOC) $(INTRA_CON_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) $(MATH) \
                         $(CODE)/energy/control/energy_control_final.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/energy_control_final.c

#------------------------------------------------------------------
test_energy.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) $(INTRA_ENT) \
                         $(INTRA_CON_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/control/test_energy.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/test_energy.c

#------------------------------------------------------------------
force_control.o   :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_ENT)  $(REAL_LOC) \
                         $(CODE)/energy/inter/real_space/force_control.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/force_control.c

#------------------------------------------------------------------
force_nolst.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
                         $(CODE)/energy/inter/real_space/force_nolst.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/force_nolst.c

#------------------------------------------------------------------
force_verlst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
                         $(CODE)/energy/inter/real_space/force_verlst.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/force_verlst.c

#------------------------------------------------------------------
force_lnklst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC) \
                         $(CODE)/energy/inter/real_space/force_lnklst.c
	$(ECHO) $@
	$(COBJ_FUSS) $(CODE)/energy/inter/real_space/force_lnklst.c

#------------------------------------------------------------------
force_npol.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_LOC) $(ENR_CTRL_ENT) $(MATH) \
                         $(CODE)/energy/inter/real_space/force_npol.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/force_npol.c

#------------------------------------------------------------------
period.o     :           $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(ENR_CTRL_ENT) \
                         $(CODE)/energy/control/period.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/control/period.c

#------------------------------------------------------------------
nbr_list_control.o   :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_ENT)  $(REAL_LOC) $(PIMD_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/inter/real_space/nbr_list_control.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/nbr_list_control.c

#------------------------------------------------------------------
verlist_control.o    :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/inter/real_space/verlist_control.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/verlist_control.c

#------------------------------------------------------------------
verlist_create.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) $(TYP_GEN) \
                         $(REAL_LOC)  $(ENR_CTRL_ENT) \
                         $(CODE)/energy/inter/real_space/verlist_create.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/real_space/verlist_create.c

#------------------------------------------------------------------
make_lnk_lst.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) \
                         $(CODE)/energy/inter/real_space/make_lnk_lst.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/energy/inter/real_space/make_lnk_lst.c

#------------------------------------------------------------------
make_lnk_map.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(FRND_ENT) \
                         $(CODE)/energy/inter/real_space/make_lnk_map.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/energy/inter/real_space/make_lnk_map.c

#------------------------------------------------------------------
lnk_lst_dis.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) $(TYP_BND) \
                         $(REAL_LOC)  $(MATH) \
                         $(CODE)/energy/inter/real_space/lnk_lst_dis.c
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/energy/inter/real_space/lnk_lst_dis.c

#------------------------------------------------------------------
ewald3d.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) \
                         $(MATH)    $(REC_ENT) \
                         $(CODE)/energy/inter/recip3d/ewald3d.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/recip3d/ewald3d.c

#------------------------------------------------------------------
ewald3d_both.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(REC_ENT) $(MATH) \
                         $(CODE)/energy/inter/recip3d/ewald3d_both.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/recip3d/ewald3d_both.c

#------------------------------------------------------------------
ewald3d_self_bgr.o  :    $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_GEN) \
                         $(REC_ENT) $(MATH) \
                         $(CODE)/energy/inter/recip3d/ewald3d_self_bgr.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/recip3d/ewald3d_self_bgr.c

#------------------------------------------------------------------
ewald3d_pme.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) \
                         $(MATH) $(REC_ENT) $(REC_LOC) $(COMM_WRAP) \
                         $(CODE)/energy/inter/recip3d/ewald3d_pme.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/inter/recip3d/ewald3d_pme.c
#------------------------------------------------------------------
surf_pot.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_CLASS) $(TYP_BND) \
                         $(TYP_GEN) $(SURF_ENT) \
                         $(CODE)/energy/surface/surf_pot.c
	$(ECHO) $@
	$(COBJ) $(CODE)/energy/surface/surf_pot.c

#======================================================================

