
#SPEC_FILES = z3dfft_ibm.o 
#spec.a : spec.a(z3dfft_ibm.o) 
#	 $(AR) $(ARFLAGS) spec.a $(SPEC_FILES)

spec.a(z3dfft_ibm.o) : z3dfft_ibm.o
z3dfft_ibm.o : $(PROTO) $(CODE)/mathlib/z3dfft_ibm.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/z3dfft_ibm.f

spec.a(math_rs.o) : math_rs.o
math_rs.o : $(PROTO) $(CODE)/mathlib/math_rs.f
	$(ECHO) $@
	$(FOBJ) $(CODE)/mathlib/math_rs.f
