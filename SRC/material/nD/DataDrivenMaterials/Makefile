include ../../../../Makefile.def

OBJS       = MaterialFileOperations.o \
            MaterialTensorOperations.o \
            MaterialDataPoint.o \
			MaterialResponseDatabase.o \
			MaterialStateVariables.o \
			DataDrivenNestedSurfaces.o \
			YieldSurfacePackage.o \
			MultiYieldSurfaceHardeningSoftening.o \
			VonMisesDMM.o \
			DruckerPragerDMM.o \
			MatsuokaNakaiDMM.o \
			OPS_VonMisesDMM.o \
			OPS_DruckerPragerDMM.o \
			OPS_MatsuokaNakaiDMM.o


all:	$(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.


