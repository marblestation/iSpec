#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from . import log
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
    out = open(segment_regions_filename, "w")
    out.write("wave_base\twave_top\n")
    out.write("\n".join(["\t".join(map(str, (seg['wave_base'], seg['wave_top']))) for seg in segment_regions]))
    out.close()

def create_segments_around_lines(linemasks, margin=0.5):
    linemasks.sort(order='wave_peak')
    dirty_segments = np.recarray((len(linemasks),),  dtype=[('wave_base', float), ('wave_top', float)])
    for i, line_mask in enumerate(linemasks):
        dirty_segments['wave_base'][i] = line_mask['wave_base'] - margin
        dirty_segments['wave_top'][i] = line_mask['wave_top'] + margin

    # Given a group of segments of a spectrum, merge those that are
    # consecutive.
    cleaned_segments = []
    i = 0
    # For all regions (except the last one), check the next one is consecutive in wavelengths
    while i < len(dirty_segments):
        j = 0
        # While wave_top of the current is equal to wave_base of the next...
        while ((j < len(dirty_segments) - 1 - i) and (dirty_segments[j+i]['wave_top'] >= dirty_segments[j+i+1]['wave_base'])):
            j += 1

        wave_base = dirty_segments[i]['wave_base']
        wave_top = dirty_segments[j+i]['wave_top']

        cleaned_segments.append((wave_base, wave_top))
        i += j + 1 # Skip the regions that have been merged
    return np.array(cleaned_segments,  dtype=[('wave_base', float), ('wave_top', float)])

def merge_overlapping_regions(regions):
    """
    Merges partially or totally overlapping regions in a NumPy recarray.

    The regions are defined by 'wave_base' and 'wave_top'. The function can handle
    two types of recarrays:
    1. With fields ('wave_base', 'wave_top')
    2. With fields ('wave_peak', 'wave_base', 'wave_top')

    If 'wave_peak' exists and regions are merged, the new 'wave_peak'
    is set to the midpoint of the new merged region's base and top.

    Args:
        regions (np.recarray): A recarray containing the regions to merge.
                               Must have 'wave_base' and 'wave_top' fields.

    Returns:
        np.recarray: A new recarray of the same type as the input, containing
                     the merged regions sorted by 'wave_base'.
    """
    # If there are 0 or 1 regions, there's nothing to merge.
    if len(regions) < 2:
        return regions.copy()

    # Sort regions by the starting point ('wave_base')
    sorted_regions = np.sort(regions, order='wave_base', kind='mergesort')

    # Check if 'wave_peak' exists and determine indices for clarity
    has_peak = 'wave_peak' in regions.dtype.names
    if has_peak:
        peak_idx, base_idx, top_idx = 0, 1, 2
    else:
        # Assign dummy index for peak to prevent errors if used
        peak_idx, base_idx, top_idx = -1, 0, 1

    # Use a list to store the merged regions. Start with the first region.
    # Convert to a list to make it mutable.
    merged_list = [list(sorted_regions[0].tolist())]

    for current_region in sorted_regions[1:]:
        last_merged = merged_list[-1]

        # Always check against the *current* top of the last_merged region,
        # as it may have been updated by a previous merge.
        if current_region['wave_base'] <= last_merged[top_idx]:
            # There is an overlap, so we merge by extending the top.
            last_merged[top_idx] = max(last_merged[top_idx], current_region['wave_top'])

            # If wave_peak exists, recalculate it based on the new boundaries.
            if has_peak:
                # The base is last_merged[base_idx], the new top is last_merged[top_idx]
                last_merged[peak_idx] = (last_merged[base_idx] + last_merged[top_idx]) / 2.0
        else:
            # No overlap, so we start a new merged region
            merged_list.append(list(current_region.tolist()))

    # Convert the list of lists back to a recarray with the original dtype
    merged_list = [tuple(m) for m in merged_list] # np.rec.fromrecords expects tuples not lists
    merged_regions = np.rec.fromrecords(merged_list, dtype=regions.dtype)

    return merged_regions
