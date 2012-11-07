# python setup.py build_ext --inplace
import os
import sys
import numpy
from distutils.core import setup
from distutils.extension import Extension
from distutils.sysconfig import get_config_vars
from Cython.Distutils import build_ext
from Cython.Build import cythonize


#os.environ['CC'] = 'gcc'
#get_config_vars()['OPT'] = ''
#'OPT': '-DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes'
#get_config_vars()['CFLAGS'] = '-fno-strict-aliasing -fno-common -dynamic -pipe -fwrapv  -DNDEBUG -g -fwrapv -O3 '
#'CFLAGS': '-fno-strict-aliasing -fno-common -dynamic -pipe -O2 -fwrapv  -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes'

ext_modules=[
    ## Requires external C routine:
    Extension(
                name="synthesizer",
                version='1.0',
                description='Python integration of SPECTRUM a Stellar Spectral Synthesis Program (C) Richard O. Gray 1992 - 2010 Version 2.76e',
                author='Sergi Blanco Cuaresma',
                url='http://www.marblestation.com',
                sources=["synthesizer.pyx", "synthesizer_func.c", "spectrum276e/abund2.c", "spectrum276e/al1op3.c", "spectrum276e/autoion3.c", "spectrum276e/balmer8.c", "spectrum276e/brackett.c", "spectrum276e/broad12.c", "spectrum276e/c1op_av.c", "spectrum276e/ca1op_av.c", "spectrum276e/capnu.c", "spectrum276e/chop.c", "spectrum276e/coolop5.c", "spectrum276e/density9.c", "spectrum276e/depth.c", "spectrum276e/eqtaukap.c", "spectrum276e/fe1op2.c", "spectrum276e/flux.c", "spectrum276e/fluxflx2.c", "spectrum276e/getisotope.c", "spectrum276e/he12.c", "spectrum276e/he13.c", "spectrum276e/he14a.c", "spectrum276e/he15a.c", "spectrum276e/he16a.c", "spectrum276e/he17a.c", "spectrum276e/he1op_av.c", "spectrum276e/he313.c", "spectrum276e/he314a.c", "spectrum276e/he315a.c", "spectrum276e/he316a.c", "spectrum276e/he317a.c", "spectrum276e/he617a.c", "spectrum276e/helines.c", "spectrum276e/helium6.c", "spectrum276e/heprof4.c", "spectrum276e/hotdensity.c", "spectrum276e/hprofl5b.c", "spectrum276e/humphreys.c", "spectrum276e/inatom2.c", "spectrum276e/infix.c", "spectrum276e/inisotope.c", "spectrum276e/inline8.c", "spectrum276e/inmodel6.c", "spectrum276e/integ4.c", "spectrum276e/intensit.c", "spectrum276e/interva4.c", "spectrum276e/invelgrad.c", "spectrum276e/isorelabun.c", "spectrum276e/linelst12b.c", "spectrum276e/lukeop2.c", "spectrum276e/lyman3.c", "spectrum276e/maxcharge.c", "spectrum276e/mg1op_av.c", "spectrum276e/mghop.c", "spectrum276e/ohop.c", "spectrum276e/opacity6.c", "spectrum276e/optdepth.c", "spectrum276e/opttrap.c", "spectrum276e/partfn5.c", "spectrum276e/paschen3.c", "spectrum276e/pfinit5.c", "spectrum276e/pfunctio.c", "spectrum276e/pfund.c", "spectrum276e/planck.c", "spectrum276e/pop13.c", "spectrum276e/qround.c", "spectrum276e/setreset.c", "spectrum276e/si1op3.c", "spectrum276e/spaux.c", "spectrum276e/strong7.c", "spectrum276e/tauflx2.c", "spectrum276e/taukap7.c", "spectrum276e/tauref.c", "spectrum276e/tauwave.c", "spectrum276e/trapez.c", "spectrum276e/unified.c", "spectrum276e/veryhotdensity.c", "spectrum276e/voigt.c", "spectrum276e/xi7.c"],

                include_dirs = [numpy.get_include()],  # .../site-packages/numpy/core/include
                extra_compile_args = [],
                extra_link_args = [],
                language="c",
    )
]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
