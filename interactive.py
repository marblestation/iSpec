#!/usr/bin/env python
import sys
import os
import getopt
import numpy as np
import sve
from sve.gui import *

## Print usage
def usage():
    print "Usage:"
    print sys.argv[0], "[--continuum=file] [--lines=file] [--segments=file] [spectrum_file]"

## Interpret arguments
def get_arguments():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "cls", ["continuum=", "lines=", "segments="])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    continuum_file = None
    lines_file = None
    segments_file = None
    spectrum_file = None
    for o, a in opts:
        if o in ("-c", "--continuum"):
            continuum_file = a
            if not os.path.exists(continuum_file):
                print "Continuum file", args[0], "does not exists!"
                sys.exit(2)
        elif o in ("-l", "--lines"):
            lines_file = a
            if not os.path.exists(lines_file):
                print "Lines file", args[0], "does not exists!"
                sys.exit(2)
        elif o in ("-s", "--segments"):
            segments_file = a
            if not os.path.exists(segments_file):
                print "Segments file", args[0], "does not exists!"
                sys.exit(2)
        else:
            print "Argument", o, "not recognized!"
            usage()
            sys.exit(2)

    filenames = {}
    filenames['spectra'] = []
    filenames['continuum'] = continuum_file
    filenames['lines'] = lines_file
    filenames['segments'] = segments_file

    # Open spectra
    for arg in args:
        spectrum_file = arg
        if not os.path.exists(spectrum_file):
            print "Spectrum file", arg, "does not exists!"
            sys.exit(2)
        filenames['spectra'].append(spectrum_file)

    return filenames


#~ Example:
#~ python interactive.py --continuum=input/regions/continuum_regions.txt --lines=input/regions/line_masks.txt --segments=input/regions/segments.txt input/spectra/examples/narval_sun.s.gz

if __name__ == '__main__':
    filenames = get_arguments()

    #### Read files

    spectra = []
    for path in filenames['spectra']:
        try:
            spectrum = sve.read_spectrum(path)
            #wfilter = (spectrum['waveobs'] >= 516.0) & (spectrum['waveobs'] <= 519.0)
            #spectrum = spectrum[wfilter]
        except Exception as e:
            print "Spectrum file", path, "has an incompatible format!"
            sys.exit(2)
        spectra.append(spectrum)

    if filenames['continuum'] != None:
        try:
            continuum = sve.read_continuum_regions(filenames['continuum'])
        except Exception as e:
            print "Continuum file", filenames['continuum'], "has an incompatible format!"
            sys.exit(2)

        ## Validations
        if np.any((continuum['wave_top'] - continuum['wave_base']) < 0):
            raise Exception("Continuum: wave_top cannot be smaller than wave_base")
    else:
        continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])

    if filenames['lines'] != None:
        try:
            lines = sve.read_line_regions(filenames['lines'])
        except Exception as e:
            print "Lines file", filenames['lines'], "has an incompatible format!"
            sys.exit(2)


        ## Validations
        if np.any((lines['wave_top'] - lines['wave_base']) < 0):
            raise Exception("Lines: wave_top cannot be smaller than wave_base")

        if np.any((lines['wave_top'] - lines['wave_peak']) < 0):
            raise Exception("Lines: wave_top cannot be smaller than wave_peak")

        if np.any((lines['wave_peak'] - lines['wave_base']) < 0):
            raise Exception("Lines: wave_peak cannot be smaller than wave_base")
    else:
        lines = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])

    if filenames['segments'] != None:
        try:
            segments = sve.read_segment_regions(filenames['segments'])
        except Exception as e:
            print "Segments file", filenames['segments'], "has an incompatible format!"
            sys.exit(2)

        ## Validations
        if np.any((segments['wave_top'] - segments['wave_base']) < 0):
            raise Exception("Segments: wave_top cannot be smaller than wave_base")
    else:
        segments = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])



    regions = {}
    regions['continuum'] = continuum
    regions['lines'] = lines
    regions['segments'] = segments

    ## Force locale to English. If for instance Spanish is the default,
    ## it may make that functions like atof(), called from external
    ## C libraries through Cython, behave interpreting comma as decimal separator
    ## and brokes the code that expects to have dots as decimal separators
    ## (i.e. SPECTRUM and its external files like linelists)
    os.environ['LANG'] = u'en_US.UTF-8'

    #app = wx.App(redirect=False) # False to avoid additional windows for stdout/stderr
    #app.frame = SpectraFrame(spectra, regions, filenames)
    #app.frame.Show()
    #app.MainLoop()


    app = SVEBaseApp(spectra, regions, filenames)
    app.mainloop()
