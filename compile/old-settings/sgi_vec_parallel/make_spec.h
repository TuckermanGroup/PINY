
#SPEC_FILES = z3dfft_sg.o 
#spec.a : spec.a(z3dfft_sg.o) 
#	 $(AR) $(ARFLAGS) spec.a $(SPEC_FILES)

spec.a(z3dfft_sg.o) : z3dfft_sg.o
z3dfft_sg.o : $(PROTO) $(CODE)/mathlib/z3dfft_sg.f
	$(ECHO) $@
	$(FOBJ_CARE) $(CODE)/mathlib/z3dfft_sg.f
