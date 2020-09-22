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
                description='Python integration of SPECTRUM a Stellar Spectral Synthesis Program (C) Richard O. Gray 1992 - 2010 Version 2.77',
                author='Sergi Blanco Cuaresma',
                url='http://www.marblestation.com',
                sources=["synthesizer.pyx", "synthesizer_func.c", "spectrum/abund2.c", "spectrum/al1op3.c", "spectrum/autoion3.c", "spectrum/balmer8.c", "spectrum/brackett.c", "spectrum/broad12.c", "spectrum/c1op_av.c", "spectrum/ca1op_av.c", "spectrum/capnu.c", "spectrum/chop.c", "spectrum/coolop5.c", "spectrum/density9.c", "spectrum/depth.c", "spectrum/eqtaukap.c", "spectrum/fe1op2.c", "spectrum/flux.c", "spectrum/fluxflx2.c", "spectrum/getisotope.c", "spectrum/he12.c", "spectrum/he13.c", "spectrum/he14a.c", "spectrum/he15a.c", "spectrum/he16a.c", "spectrum/he17a.c", "spectrum/he1op_av.c", "spectrum/he313.c", "spectrum/he314a.c", "spectrum/he315a.c", "spectrum/he316a.c", "spectrum/he317a.c", "spectrum/he617a.c", "spectrum/helines.c", "spectrum/helium6.c", "spectrum/heprof4.c", "spectrum/hotdensity.c", "spectrum/hprofl5.c", "spectrum/humphreys.c", "spectrum/inatom2.c", "spectrum/infix.c", "spectrum/inisotope.c", "spectrum/inline8.c", "spectrum/inmodel6.c", "spectrum/integ4.c", "spectrum/intensit.c", "spectrum/interva4.c", "spectrum/invelgrad.c", "spectrum/isorelabun.c", "spectrum/linelst12b.c", "spectrum/lline6.c", "spectrum/lukeop2.c", "spectrum/lyman3.c", "spectrum/maxcharge.c", "spectrum/mg1op_av.c", "spectrum/mghop.c", "spectrum/ohop.c", "spectrum/opacity6.c", "spectrum/optdepth.c", "spectrum/opttrap.c", "spectrum/partfn5.c", "spectrum/paschen3.c", "spectrum/pfinit5.c", "spectrum/pfunctio.c", "spectrum/pfund.c", "spectrum/planck.c", "spectrum/pop13.c", "spectrum/qround.c", "spectrum/setreset.c", "spectrum/si1op3.c", "spectrum/spaux.c", "spectrum/strong8.c", "spectrum/tauflx2.c", "spectrum/taukap7.c", "spectrum/tauref.c", "spectrum/tauwave.c", "spectrum/trapez.c", "spectrum/unified.c", "spectrum/veryhotdensity.c", "spectrum/voigt.c", "spectrum/xi7.c"],

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
