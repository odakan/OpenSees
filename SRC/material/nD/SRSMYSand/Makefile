include ../../../../Makefile.def

OBJS       = CTensor.o \
			MaterialStateVariables.o \
			NestedSurface.o \
			SRSMYSand.o \
			OPS_SRSMYSand.o

all:	$(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



