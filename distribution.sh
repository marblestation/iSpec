#!/bin/bash
cd pyinstaller/
rm -rf sve/build/
rm -rf sve/dist/
python utils/Build.py sve/sve.spec
cd ..
rm -rf ./pyinstaller/sve/dist/input/spectra/examples/ ./pyinstaller/sve/dist/*.command
cp -rf *.command ./pyinstaller/sve/dist/
mkdir -p ./pyinstaller/sve/dist/input/spectra/examples/
cp -rf input/spectra/examples/*.s.gz ./pyinstaller/sve/dist/input/spectra/examples/
rm -f ./pyinstaller/sve/dist/input/spectra/examples/*_norm.s.gz
mkdir -p ./pyinstaller/sve/dist/input/regions/
cp -rf input/regions/*.txt ./pyinstaller/sve/dist/input/regions/
rm -rf ./pyinstaller/sve/dist/input/spectra/binaries/
mkdir -p ./pyinstaller/sve/dist/input/spectra/binaries/
cp -f input/spectra/binaries/*.s.gz ./pyinstaller/sve/dist/input/spectra/binaries/
echo "-----------------------------------------------"
echo " SVE distribution in './pyinstaller/sve/dist/'"
echo "-----------------------------------------------"
cd ..
cd docs-sphinx/
make latexpdf
echo "----------------------------------------------------------"
echo " SVE distribution in './docs-sphinx/_build/latex/SVE.pdf'"
echo "----------------------------------------------------------"

