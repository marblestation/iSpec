#
#    This file is part of Spectra Visual Editor (SVE).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    SVE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SVE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with SVE. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import log
import logging

def read_segment_regions(segment_regions_filename):
    """
    Read segment regions.
    The specified file should be plain text with **tab** character as column delimiter.
    Two columns should exists: 'wave_base' and 'wave_top' (the first line should contain those header names).
    They indicate the beginning and end of each region (one per line). For instance:
    ::

        wave_base       wave_top
        480.6000        480.6100
        481.1570        481.1670
        491.2240        491.2260
        492.5800        492.5990
    """
    segment_regions = np.array([tuple(seg.rstrip('\r\n').split("\t")) for seg in open(segment_regions_filename,)][1:], dtype=[('wave_base', float),('wave_top', float)])

    if np.any(segment_regions['wave_top'] - segment_regions['wave_base'] <= 0):
        logging.error("Segments regions where wave_base is equal or bigger than wave_top")
        raise Exception("Incompatible format")
    return segment_regions

def write_segment_regions(segment_regions, segment_regions_filename):
    """
    Write segment regions file with the following format:
    ::

        wave_base       wave_top
        480.6000        480.6100
        481.1570        481.1670
        491.2240        491.2260
        492.5800        492.5990
    """
    out = open(segments_regions_filename, "w")
    out.write("wave_base\twave_top\n")
    out.write("\n".join(["\t".join(map(str, (seg['wave_base'], seg['wave_top']))) for seg in segments_regions]))
    out.close()


