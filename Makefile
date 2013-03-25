
all: sve/synthesizer.so


sve/synthesizer.so: synthesizer/synthesizer.pyx synthesizer/synthesizer_func.c
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	cd synthesizer/ ; python setup.py build_ext --inplace
	cp -f synthesizer/synthesizer.so sve/

.PHONY: clean pyinstaller docs


clean:
	rm -f synthesizer/synthesizer.c
	rm -rf synthesizer/build/
	rm -f synthesizer/synthesizer.so
	rm -f sve/synthesizer.so
	# pyinstaller
	rm -rf pyinstaller/sve/build/
	rm -rf pyinstaller/sve/dist/
	# docs
	cd docs-sphinx/; make clean

pyinstaller:
	rm -rf pyinstaller/sve/build/
	rm -rf pyinstaller/sve/dist/
	cd pyinstaller/; python utils/Build.py sve/sve.spec
	rm -rf pyinstaller/sve/dist/input/spectra/examples/
	rm -rf pyinstaller/sve/dist/*.command
	cp -rf *.command ./pyinstaller/sve/dist/
	mkdir -p pyinstaller/sve/dist/input/spectra/examples/
	cp -rf input/spectra/examples/*.s.gz pyinstaller/sve/dist/input/spectra/examples/
	rm -f pyinstaller/sve/dist/input/spectra/examples/*_norm.s.gz
	mkdir -p pyinstaller/sve/dist/input/regions/
	cp -rf input/regions/*.txt pyinstaller/sve/dist/input/regions/
	rm -rf pyinstaller/sve/dist/input/spectra/binaries/
	mkdir -p pyinstaller/sve/dist/input/spectra/binaries/
	cp -f input/spectra/binaries/*.s.gz pyinstaller/sve/dist/input/spectra/binaries/
	echo "-----------------------------------------------"
	echo " SVE distribution in './pyinstaller/sve/dist/'"
	echo "-----------------------------------------------"

docs:
	cd docs-sphinx/; make latexpdf
	echo "----------------------------------------------------------"
	echo " SVE distribution in './docs-sphinx/_build/latex/SVE.pdf'"
	echo "----------------------------------------------------------"

