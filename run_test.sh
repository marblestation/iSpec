#!/bin/bash
echo "1) MRD Sun"
echo "2) MRG Mu Leo"
echo "3) MPD Mu Cas A"
echo "4) MPG Arcturus"
echo "5) RV test with Mu Cas A"
echo "6) Exit"
echo "Enter the spectra number you want to use : "
read i;

case $i in
1) python interactive.py --continuum=input/LUMBA/UVES_MRD_sun_cmask.txt --lines=input/LUMBA/UVES_MRD_sun_Fe-linelist.txt --segments=input/LUMBA/UVES_MRD_sun_segments.txt input/LUMBA/UVES_MRD_sun_official.s.gz ;;
2) python interactive.py --continuum=input/LUMBA/UVES_MRG_mu_leo_cmask.txt --lines=input/LUMBA/UVES_MRG_mu_leo_Fe-linelist.txt --segments=input/LUMBA/UVES_MRG_mu_leo_segments.txt input/LUMBA/UVES_MRG_mu_leo_official.s.gz ;;
3) python interactive.py --continuum=input/LUMBA/UVES_MPD_mu_cas_a_cmask.txt --lines=input/LUMBA/UVES_MPD_mu_cas_a_Fe-linelist.txt --segments=input/LUMBA/UVES_MPD_mu_cas_a_segments.txt input/LUMBA/UVES_MPD_mu_cas_a_official.s.gz ;;
4) python interactive.py --continuum=input/LUMBA/UVES_MPG_arcturus_cmask.txt --lines=input/LUMBA/UVES_MPG_arcturus_Fe-linelist.txt --segments=input/LUMBA/UVES_MPG_arcturus_segments.txt input/LUMBA/UVES_MPG_arcturus_official.s.gz ;;
5) python interactive.py --continuum=input/LUMBA/UVES_MPD_mu_cas_a_cmask.txt --lines=input/RV_linelist.200.txt --segments=input/LUMBA/UVES_MPD_mu_cas_a_segments.txt input/LUMBA/UVES_MPD_mu_cas_a_official.s.gz ;;
6) exit
esac
