#!/bin/bash
rm -f synthesizer.so
cd synthesizer/
rm -f synthesizer.so
rm -f synthesizer.c
rm -rf build/
python setup.py build_ext --inplace
#python setup.py build --build_lib .
cp -f synthesizer.so ../ispec/
cd ..

rm -f isochrones/YYmix2
gfortran -o isochrones/YYmix2 isochrones/YYmix2.f
