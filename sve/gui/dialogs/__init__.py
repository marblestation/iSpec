"""
    This file is part of Spectra Visual Editor (SVE).
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com

    SVE is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SVE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with SVE. If not, see <http://www.gnu.org/licenses/>.
"""
from __future__ import division

import os
import sys

if sys.hexversion < 0x02050000:
    raise RuntimeError("SVE requires at least Python 2.5")

from AbundancesDialog import *
from AddNoiseDialog import *
from AdjustLinesDialog import *
from CleanSpectrumDialog import *
from CleanTelluricsDialog import *
from CombineSpectraDialog import *
from CorrectVelocityDialog import *
from CustomDialog import *
from CutSpectrumDialog import *
from DegradeResolutionDialog import *
from DetermineBarycentricCorrectionDialog import *
from EstimateErrorsDialog import *
from EstimateSNRDialog import *
from ExampleDialog import *
from FindContinuumDialog import *
from FindLinesDialog import *
from FindSegmentsDialog import *
from FitContinuumDialog import *
from FitLinesDialog import *
from OperateSpectrumDialog import *
from ResampleSpectrumDialog import *
from SendSpectrumDialog import *
from SyntheticSpectrumDialog import *
from VelocityProfileDialog import *
from SolverDialog import *

