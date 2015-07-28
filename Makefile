
all: synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu synthesizer/turbospectrum/bin/eqwidt_lu ispec/synthesizer.so isochrones/YYmix2


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

synthesizer/turbospectrum/bin/babsma_lu synthesizer/turbospectrum/bin/bsyn_lu synthesizer/turbospectrum/bin/eqwidt_lu: synthesizer/turbospectrum/src/*.f
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/*.o
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/babsma_lu
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/bsyn_lu
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/eqwidt_lu
	$(MAKE) -C synthesizer/turbospectrum/exec-gf-v15.1/
	mv synthesizer/turbospectrum/exec-gf-v15.1/babsma_lu synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf-v15.1/bsyn_lu synthesizer/turbospectrum/bin/
	mv synthesizer/turbospectrum/exec-gf-v15.1/eqwidt_lu synthesizer/turbospectrum/bin/


.PHONY: clean


clean:
	rm -f ispec.log
	rm -f synthesizer/turbospectrum/exec-gf-v15.1/*.o
	rm -f synthesizer/turbospectrum/bin/babsma_lu 
	rm -f synthesizer/turbospectrum/bin/bsyn_lu 
	rm -f synthesizer/turbospectrum/bin/eqwidt_lu
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	rm -f synthesizer/synthesizer.so
	rm -f ispec/synthesizer.so
	rm -f isochrones/YYmix2

