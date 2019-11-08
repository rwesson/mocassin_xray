# standard
FC = mpif90
LD = mpif90

#IBM
#FC = mpxlf90_r
#LD = mpxlf90_r
#FFLAGS += -O3 -qstrict -qrealsize=4 -bmaxdata:0x40000000 -qmaxmem=-1 -cpp -DPREFIX=\"$(PREFIX)\"

#SUN
#FC = mpf90
#LD = mpf90
#FFLAGS += -fast -xarch=v9b -ftrap=%none -lmpi -cpp -DPREFIX=\"$(PREFIX)\"

#SGI
#FC = f90
#LD = f90
#FFLAGS += -64 -C -mpio -OPT:Olimit=3495 -O3 -lmpi -cpp -DPREFIX=\"$(PREFIX)\"

# set prefix depending on OS
OS := $(shell uname)
ifeq ($(OS),Darwin)
  PREFIX=/usr/local
else
  PREFIX=/usr
endif

# get version from changelog if debian package, or git log otherwise
VERSION := $(shell if [ -e debian/ ]; then dpkg-parsechangelog -S version; else git describe --always --tags --dirty || echo "3.04.00"; fi)

# set flags

ifeq ($(FC),ifort)
  FFLAGS += -cpp -DPREFIX=\"$(PREFIX)\" -DVERSION=\"$(VERSION)\" -module source/
else
  FFLAGS += -cpp -Jsource/ -ffree-line-length-0 -lm -DPREFIX=\"$(PREFIX)\" -DVERSION=\"$(VERSION)\" -DCOMPILER=\"$(FC)\" -DCO=\"$(CO)\" -I/usr/include
endif

MANDIR=$(DESTDIR)$(PREFIX)/share/man/man1
SOURCES = source/constants_mod.o source/vector_mod.o source/common_mod.o source/interpolation_mod.o \
	source/set_input_mod.o source/hydro_mod.o source/ph_mod.o source/composition_mod.o \
	source/continuum_mod.o source/ionization_mod.o source/pathIntegration_mod.o \
	source/grid_mod.o source/dust_mod.o source/emission_mod.o source/photon_mod.o  \
	source/update_mod.o source/output_mod.o source/iteration_mod.o source/readdata_mod.o

ifeq ($(CO),debug) #to show all compiler warnings
  FFLAGS += -fbounds-check -Wall -Wuninitialized -g -pg #-ffpe-trap=zero,overflow,invalid,underflow,denormal -fbacktrace -fcheck=all
else ifeq ($(CO),valgrind)
  FFLAGS += -g
else ifeq ($(CO),gprof)
  FFLAGS += -pg
else
  FFLAGS += -O2
endif

.PHONY: all clean new install uninstall

all: mocassinX mocassinXWarm mocassinXOutput mocassinXPlot mocassinXFluorescence

new: clean all

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

mocassinX: $(SOURCES) source/mocassinX.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinXWarm: $(SOURCES) source/mocassinXWarm.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinXOutput: $(SOURCES) source/mocassinXOutput.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinXPlot: $(SOURCES) source/mocassinXPlot.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

mocassinXFluorescence: $(SOURCES) source/mocassinXFluorescence.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^

clean:
	/bin/rm -f source/*.o *~ source/*.mod mocassinX mocassinXWarm mocassinXOutput mocassinXPlot mocassinXFluorescence

install: mocassinX mocassinXWarm mocassinXOutput mocassinXPlot
	test -e $(DESTDIR)$(PREFIX)/share/mocassinX || mkdir -p $(DESTDIR)$(PREFIX)/share/mocassinX
	test -e $(DESTDIR)$(PREFIX)/bin || mkdir -p $(DESTDIR)$(PREFIX)/bin
	test -e $(MANDIR) || mkdir -p $(MANDIR)
	cp -R data $(DESTDIR)$(PREFIX)/share/mocassinX
	cp -R dustData $(DESTDIR)$(PREFIX)/share/mocassinX
	cp -R benchmarks $(DESTDIR)$(PREFIX)/share/mocassinX
	cp -R examples $(DESTDIR)$(PREFIX)/share/mocassinX
#	install -m 644 man/mocassinX.1 $(MANDIR)
#	gzip -f $(MANDIR)/mocassinX.1
#	ln -s -f $(MANDIR)/mocassinX.1.gz $(MANDIR)/mocassinXWarm.1.gz
#	ln -s -f $(MANDIR)/mocassinX.1.gz $(MANDIR)/mocassinXOutput.1.gz
#	ln -s -f $(MANDIR)/mocassinX.1.gz $(MANDIR)/mocassinXPlot.1.gz
	install mocassinX $(DESTDIR)$(PREFIX)/bin
	install mocassinXWarm $(DESTDIR)$(PREFIX)/bin
	install mocassinXPlot $(DESTDIR)$(PREFIX)/bin
	install mocassinXOutput $(DESTDIR)$(PREFIX)/bin
	install mocassinXFluorescence $(DESTDIR)$(PREFIX)/bin

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/mocassinX*
	rm -f $(MANDIR)/mocassinX*.1.gz
	rm -rf $(DESTDIR)$(PREFIX)/share/mocassinX
