"""
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
    
    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
cimport numpy as np

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b
#cdef inline int float_max(float a, float b): return a if a >= b else b
#cdef inline int float_min(float a, float b): return a if a <= b else b

from cpython cimport bool

cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
def find_local_max_values(np.ndarray[np.double_t,ndim=1] x):
    """
    For an array of values, find the position of local maximum values considering only
    the next and previous elements, except they have the same value.
    In that case, the next/previous different value is checked. Therefore,
    ::

        find_local_max([1,2,3,3,2,1,4,3])

    would return:
    ::

        [2, 3, 6]
    """
    ret = []
    cdef int i
    cdef int j
    cdef int n = len(x)
    cdef int l_min
    cdef int r_min
    cdef bool is_min
    for i from 0 <= i < n by 1:
        l_min = int_max(i-1, 0)
        #l_max = i-1
        #r_min = i+1
        #r_max = np.min([i+1, n-1])
        r_min = int_min(i+1, n-1)
        is_max = True

        # left side
        j = l_min
        # If value is equal, search for the last different value
        while j >= 0 and x[j] == x[i]:
            j -= 1

        if (j < 0 or x[j] > x[i]) and i > 0:
            is_max = False

        # right side
        if is_max:
            j = r_min
            # If value is equal, search for the next different value
            while j < n and x[j] == x[i]:
                j += 1
            if (j >= n or x[j] > x[i]) and i < n-1:
                is_max = False

        if is_max:
            ret.append(i)
    return np.asarray(ret)


@cython.boundscheck(False)
@cython.wraparound(False)
def find_local_min_values(np.ndarray[np.double_t,ndim=1] x):
    """
    For an array of values, find the position of local maximum values considering only
    the next and previous elements, except they have the same value.
    In that case, the next/previous different value is checked. Therefore,
    ::

        find_local_max([10,9,3,3,9,10,4,30])

    would return:
    ::

        [2, 3, 6]
    """
    ret = []
    cdef int i
    cdef int j
    cdef int n = len(x)
    cdef int l_min
    cdef int r_min
    cdef bool is_min
    for i from 0 <= i < n by 1:
        l_min = int_max(i-1, 0)
        #l_max = i-1
        #r_min = i+1
        #r_max = np.min([i+1, n-1])
        r_min = int_min(i+1, n-1)
        is_min = True
        # left side
        j = l_min
        # If value is equal, search for the last different value
        while j >= 0 and x[j] == x[i]:
            j -= 1

        if j < 0 or x[j] < x[i]:
            is_min = False

        # right side
        if is_min:
            j = r_min
            # If value is equal, search for the next different value
            while j < n and x[j] == x[i]:
                j += 1

            if j >= n or x[j] < x[i]:
                is_min = False

        if is_min:
            ret.append(i)
    return np.asarray(ret)



