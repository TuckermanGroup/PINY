#SPEC_FILES = z3dfft_dec.o 
#spec.a : spec.a(z3dfft_dec.o) 
#	 $(AR) $(ARFLAGS) spec.a $(SPEC_FILES)

spec.a(z3dfft_dec_noimsl.o) : z3dfft_dec_noimsl.o
z3dfft_dec_noimsl.o : $(PROTO) $(CODE)/mathlib/z3dfft_dec_noimsl.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/z3dfft_dec_noimsl.f

spec.a(math_rs.o) : math_rs.o
math_rs.o : $(PROTO) $(CODE)/mathlib/math_rs.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/math_rs.f


