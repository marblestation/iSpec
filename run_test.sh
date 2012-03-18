#!/bin/bash
echo "------------------------------------------------------------------"
echo "     Normalized / R = 47000 / Wavelength range = 480 - 680 nm"
echo "------------------------------------------------------------------"
echo "1) Sun       (Metal Rich Dwarf) - Teff 5777 logg 4.44 [M/H]  0.00"
echo "2) Mu Cas A  (Metal Poor Dwarf) - Teff 5308 logg 4.41 [M/H] -0.89"
echo "3) Arcturus  (Metal Poor Giant) - Teff 4247 logg 1.59 [M/H] -0.54"
echo "4) Mu Leo    (Metal Rich Giant) - Teff 4433 logg 2.50 [M/H]  0.29"
echo "5) All of them simultaneously"
echo "6) Exit"
echo "------------------------------------------------------------------"
echo "Enter the spectra number you want to use : "
read i;

case $i in
1) python interactive.py --continuum=input/test/continuum_regions.txt --lines=input/test/line_masks.txt --segments=input/test/segments.txt input/test/observed_sun.s.gz ;;
2) python interactive.py --continuum=input/test/continuum_regions.txt --lines=input/test/line_masks.txt --segments=input/test/segments.txt input/test/observed_mu_cas_a.s.gz ;;
3) python interactive.py --continuum=input/test/continuum_regions.txt --lines=input/test/line_masks.txt --segments=input/test/segments.txt input/test/observed_arcturus.s.gz ;;
4) python interactive.py --continuum=input/test/continuum_regions.txt --lines=input/test/line_masks.txt --segments=input/test/segments.txt input/test/observed_mu_leo.s.gz ;;
5) python interactive.py --continuum=input/test/continuum_regions.txt --lines=input/test/line_masks.txt --segments=input/test/segments.txt input/test/observed_sun.s.gz input/test/observed_mu_leo.s.gz input/test/observed_mu_cas_a.s.gz input/test/observed_arcturus.s.gz ;;
6) exit
esac
