
UNAME_S := $(shell uname -s)


all: spectrum turbospectrum moog isochrones

spectrum: ispec/synthesizer.so synthesizer/spectrum/spectrum

turbospectrum: synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu

moog: synthesizer/moog/MOOGSILENT

ares: synthesizer/ARES/bin/ARES

isochrones: isochrones/YYmix2


ispec/synthesizer.so: synthesizer/synthesizer.pyx synthesizer/synthesizer_func.c synthesizer/spectrum/*.c
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	cd synthesizer/ ; python setup.py build_ext --inplace
	mv -f synthesizer/synthesizer*.so ispec/synthesizer.so
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/

synthesizer/spectrum/spectrum: synthesizer/spectrum/*.c
	rm -f synthesizer/spectrum/*.o
	$(MAKE) -C synthesizer/spectrum/

isochrones/YYmix2: isochrones/YYmix2.f
	rm -f isochrones/YYmix2
	gfortran -L/usr/lib -o isochrones/YYmix2 isochrones/YYmix2.f

synthesizer/ARES/bin/ARES:	synthesizer/ARES/src/ARES_v2.c synthesizer/ARES/src/*.h
	mkdir -p synthesizer/ARES/bin/
	gcc -o synthesizer/ARES/bin/ARES synthesizer/ARES/src/ARES_v2.c -lcfitsio -lgsl -lgslcblas -lm -lgomp -fopenmp ${CPPFLAGS} ${LDFLAGS}

synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu: synthesizer/turbospectrum/source/*.f
	rm -f synthesizer/turbospectrum/exec-gf/*.o
	rm -f synthesizer/turbospectrum/exec-gf/*.mod
	rm -f synthesizer/turbospectrum/exec-gf/babsma_lu
	rm -f synthesizer/turbospectrum/exec-gf/bsyn_lu
	$(MAKE) -C synthesizer/turbospectrum/exec-gf/
	mkdir -p synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf/babsma_lu synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf/bsyn_lu synthesizer/turbospectrum/bin/

synthesizer/moog/MOOGSILENT: synthesizer/moog/*.f
ifeq ($(UNAME_S),Linux)
	sed -i 's/machine = "mac"/machine = "pcl"/' synthesizer/moog/Moogsilent.f
	sed -i 's/machine = "uni"/machine = "pcl"/' synthesizer/moog/Moogsilent.f
endif
ifeq ($(UNAME_S),Darwin)
	sed -i.bak 's/machine = "pcl"/machine = "mac"/' synthesizer/moog/Moogsilent.f
	sed -i.bak 's/machine = "uni"/machine = "mac"/' synthesizer/moog/Moogsilent.f
	rm -f synthesizer/moog/Moogsilent.f.bak
endif
ifeq ($(UNAME_S),Solaris)
	sed -i 's/machine = "mac"/machine = "uni"/' synthesizer/moog/Moogsilent.f
	sed -i 's/machine = "pcl"/machine = "uni"/' synthesizer/moog/Moogsilent.f
endif
	rm -f synthesizer/moog/*.o
	rm -f synthesizer/moog/MOOG
	rm -f synthesizer/moog/MOOGSILENT
	$(MAKE) -C synthesizer/moog/fake_sm-2.4.35/ -f Makefile
	$(MAKE) -C synthesizer/moog/ -f Makefile.rhsilent
ifeq ($(UNAME_S),Darwin)
	sed -i.bak 's/machine = "mac"/machine = "pcl"/' synthesizer/moog/Moogsilent.f
	rm -f synthesizer/moog/Moogsilent.f.bak
endif
ifeq ($(UNAME_S),Solaris)
	sed -i 's/machine = "uni"/machine = "pcl"/' synthesizer/moog/Moogsilent.f
endif


.PHONY: clean


clean:
	rm -f ispec.log
	rm -f synthesizer/moog/fake_sm-2.4.35/lib/libfakesm.a
	rm -f synthesizer/moog/fake_sm-2.4.35/src/fakesm.o
	rm -f synthesizer/moog/*.o
	rm -f synthesizer/moog/MOOG
	rm -f synthesizer/moog/MOOGSILENT
	rm -f synthesizer/turbospectrum/exec-gf/*.o
	rm -f synthesizer/turbospectrum/exec-gf/*.mod
	rm -f synthesizer/turbospectrum/bin/babsma_lu 
	rm -f synthesizer/turbospectrum/bin/bsyn_lu 
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	rm -f synthesizer/synthesizer*.so
	rm -f synthesizer/spectrum/*.o
	rm -f synthesizer/spectrum/spectrum
	rm -f ispec/synthesizer.so
	rm -f isochrones/YYmix2
	rm -f synthesizer/ARES/bin/ARES

