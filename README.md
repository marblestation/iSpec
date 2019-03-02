
# iSpec

[iSpec](https://www.blancocuaresma.com/s/) is a tool for the treatment and analysis of stellar spectra. Some of the main functionalities for spectra treatment are the following:

- Cosmic rays removal
- Continuum normalization
- Resolution degradation
- Radial velocity determination and correction
- Telluric lines identification
- Re-sampling

## Atmospheric parameters and chemical abundances.

iSpec is capable of determining atmospheric parameters (i.e effective temperature, surface gravity, metallicity, micro/macroturbulence, rotation) and individual chemical abundances for AFGKM stars by using two different approaches: synthetic spectra fitting technique or equivalent widths method. iSpec integrates [MARCS](http://marcs.astro.uu.se) and [ATLAS](http://kurucz.harvard.edu/) model atmospheres together with the following radial transfer codes:

- [SPECTRUM](http://www.appstate.edu/~grayro/spectrum/spectrum.html) R. O. Gray
- [Turbospectrum](http://www.appstate.edu/~grayro/spectrum/spectrum.html) Bertrand Plez
- [SME](http://www.stsci.edu/~valenti/sme.html) Valenti & Piskunov
- [MOOG](http://www.as.utexas.edu/~chris/moog.html) Chris Sneden
- [Synthe/WIDTH9](http://atmos.obspm.fr) Kurucz/Atmos


## Python powered automatized analyses.

The user-friendly interface is perfect for learning and testing. However, to take advantage of the full potential, iSpec can be used from Python. This is the recommended way to use iSpec for complex scientific studies, it ensures reproducibility and give access to a wider range of functionalities and options.


## Citations

If you use iSpec, we thank you to cite these two articles: [A&A (2014)](http://adsabs.harvard.edu/abs/2014A%26A...569A.111B) and [MNRAS (2019)](http://adsabs.harvard.edu/abs/2019arXiv190209558B).

[More information and documentation](https://www.blancocuaresma.com/s/)

