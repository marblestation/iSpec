
all: spectrum turbospectrum moog isochrones

spectrum: ispec/synthesizer.so

turbospectrum: synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu synthesizer/turbospectrum/bin/eqwidt_lu

moog: synthesizer/moog/MOOGSILENT

ares: synthesizer/ARES/bin/ARES

isochrones: isochrones/YYmix2


ispec/synthesizer.so: synthesizer/synthesizer.pyx synthesizer/synthesizer_func.c synthesizer/spectrum276e/*.c
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	cd synthesizer/ ; python setup.py build_ext --inplace
	mv -f synthesizer/synthesizer.so ispec/
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/

isochrones/YYmix2: isochrones/YYmix2.f
	rm -f isochrones/YYmix2
	gfortran -o isochrones/YYmix2 isochrones/YYmix2.f

synthesizer/ARES/bin/ARES:	synthesizer/ARES/src/ARES_v2.c synthesizer/ARES/src/*.h
	mkdir -p synthesizer/ARES/bin/
	gcc -o synthesizer/ARES/bin/ARES synthesizer/ARES/src/ARES_v2.c -lcfitsio -lgsl -lgslcblas -lm -lgomp -fopenmp ${CPPFLAGS} ${LDFLAGS}

synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu synthesizer/turbospectrum/bin/eqwidt_lu: synthesizer/turbospectrum/src/*.f
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/*.o
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/babsma_lu
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/bsyn_lu
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/eqwidt_lu
	$(MAKE) -C synthesizer/turbospectrum/exec-gf-v15.1/
	mkdir -p synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf-v15.1/babsma_lu synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf-v15.1/bsyn_lu synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf-v15.1/eqwidt_lu synthesizer/turbospectrum/bin/

synthesizer/moog/MOOGSILENT: synthesizer/moog/*.f
	rm -f synthesizer/moog/*.o
	rm -f synthesizer/moog/MOOG
	rm -f synthesizer/moog/MOOGSILENT
	$(MAKE) -C synthesizer/moog/fake_sm-2.4.35/ -f Makefile
	$(MAKE) -C synthesizer/moog/ -f Makefile.rhsilent


.PHONY: clean


clean:
	rm -f ispec.log
	rm -f synthesizer/moog/fake_sm-2.4.35/lib/libfakesm.a
	rm -f synthesizer/moog/fake_sm-2.4.35/src/fakesm.o
	rm -f synthesizer/moog/*.o
	rm -f synthesizer/moog/MOOG
	rm -f synthesizer/moog/MOOGSILENT
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/*.o
	rm -f synthesizer/turbospectrum/bin/babsma_lu 
	rm -f synthesizer/turbospectrum/bin/bsyn_lu 
	rm -f synthesizer/turbospectrum/bin/eqwidt_lu
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	rm -f synthesizer/synthesizer.so
	rm -f ispec/synthesizer.so
	rm -f isochrones/YYmix2
	rm -f synthesizer/ARES/bin/ARES

