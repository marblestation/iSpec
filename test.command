#!/bin/bash
echo "----------------------------------------------------------------------------"
echo "                 Resolving power per instrument"
echo "----------------------------------------------------------------------------"
echo "  Narval   ~70,000 - 90,000"
echo "  ESPaDOnS ~68,000 - 81,000"
echo "  HARPS    ~115,000"
echo "----------------------------------------------------------------------------"
echo "                 Wavelength range = 480 - 680 nm"
echo "----------------------------------------------------------------------------"
echo "1) Procyon   (Metal Rich Dwarf) - Teff 6545 logg 3.99 [M/H] -0.02 - HARPS"
echo "2) Sun       (Metal Rich Dwarf) - Teff 5777 logg 4.44 [M/H]  0.00 - Narval"
echo "3) Mu Cas A  (Metal Poor Dwarf) - Teff 5308 logg 4.41 [M/H] -0.89 - Narval"
echo "4) Arcturus  (Metal Poor Giant) - Teff 4247 logg 1.59 [M/H] -0.54 - Narval"
echo "5) Mu Leo    (Metal Rich Giant) - Teff 4433 logg 2.50 [M/H]  0.29 - ESPaDOnS"
echo "............................................................................"
echo "6) All of them simultaneously (1 to 5)"
echo "7) Sun with example region masks"
echo "8) Exit"
echo "----------------------------------------------------------------------------"
echo "Enter the spectra number you want to use : "
read i;

# Directory of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Default exeuction
SVE='/usr/bin/env python '$DIR'/interactive.py'

if [ ! -f interactive.py ]; then
    platform=`uname`
    if [[ "$platform" == 'Linux' && -f $DIR/sve.linux2/sve.linux2 ]]; then
        SVE=$DIR'/sve.linux2/sve.linux2'
    elif [[ "$platform" == 'Darwin' && -f $DIR/sve.darwin/sve.darwin ]]; then
        SVE=$DIR'/sve.darwin/sve.darwin'
    fi
fi

case $i in
1) $SVE $DIR/input/spectra/examples/harps_procyon.s.gz ;;
2) $SVE $DIR/input/spectra/examples/narval_sun.s.gz ;;
3) $SVE $DIR/input/spectra/examples/narval_mu_cas.s.gz ;;
4) $SVE $DIR/input/spectra/examples/narval_arcturus.s.gz ;;
5) $SVE $DIR/input/spectra/examples/espadons_mu_leo.s.gz ;;
6) $SVE $DIR/input/spectra/examples/espadons_mu_leo.s.gz $DIR/input/spectra/examples/narval_arcturus.s.gz $DIR/input/spectra/examples/harps_procyon.s.gz $DIR/input/spectra/examples/narval_mu_cas.s.gz $DIR/input/spectra/examples/narval_sun.s.gz ;;
7) $SVE --continuum=$DIR/input/regions/continuum_regions.txt --lines=$DIR/input/regions/line_masks.txt --segments=$DIR/input/regions/segments.txt $DIR/input/spectra/examples/narval_sun.s.gz ;;
8) exit
esac
