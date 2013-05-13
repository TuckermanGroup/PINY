#===================================================================
#               MATH_FILES
#===================================================================


#=================================================================
mathlib.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) \
                         $(MATH) $(COMM_WRAP) \
                         $(CODE)/mathlib/mathlib.c
	$(ECHO) $@
	$(COBJ) $(CODE)/mathlib/mathlib.c

#------------------------------------------------------------------
blas_wrappers.o       :  $(DEFINES) $(STANDARD) $(MATH) \
                         $(CODE)/mathlib/blas_wrappers.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/blas_wrappers.f

#------------------------------------------------------------------
fft_generic.o       :    $(DEFINES) $(STANDARD) $(MATH) \
                         $(CODE)/mathlib/fft_generic.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/fft_generic.f

#------------------------------------------------------------------
math_generic.o       :   $(DEFINES) $(STANDARD) $(MATH) \
                         $(CODE)/mathlib/math_generic.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/math_generic.f

#------------------------------------------------------------------
math_rs.o       :   $(DEFINES) $(STANDARD) $(MATH) \
                         $(CODE)/mathlib/math_rs.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/math_rs.f

#------------------------------------------------------------------
math_rs_v2.o       :   $(DEFINES) $(STANDARD) $(MATH) \
                         $(CODE)/mathlib/math_rs_v2.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/math_rs_v2.f

#------------------------------------------------------------------
fft_package.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(MATH) $(FRND_ENT) $(ENR_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/mathlib/fft_package.c
	$(ECHO) $@
	$(COBJ) $(CODE)/mathlib/fft_package.c

#------------------------------------------------------------------
fft_create_package.o     :      $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CP) $(TYP_CLASS) \
                         $(MATH) $(FRND_ENT) $(ENR_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/mathlib/fft_create_package.c
	$(ECHO) $@
	$(COBJ_NOOPT) $(CODE)/mathlib/fft_create_package.c

#===================================================================
