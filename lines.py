import asciitable
import numpy as np


# Load a VALD linelist (short format) and filter it
# - depth_limit: filter out all the lines with a depth less than this percentage (0.05 = 5%)
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def load_and_filter_VALD(vald_file, depth_limit=0.0, data_end=None):
    # Original VALD linelist
    if data_end == None:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "Rad", "Stark", "Waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "Rad", "Stark", "Waals", "factor", "Reference"], guess=False)
    else:
        vald = asciitable.read(vald_file, delimiter=",", quotechar="'", data_start=3, data_end=data_end, names=["element", "wave (A)", "lower state (eV)", "Vmic (km/s)", "log(gf)", "Rad", "Stark", "Waals", "factor", "depth", "Reference"], exclude_names=["Vmic (km/s)", "Rad", "Stark", "Waals", "factor", "Reference"], guess=False)

    if depth_limit <= 0.0:
        return vald
    else:
        # Filter
        vald_limited = vald[vald['depth'] >= depth_limit]
        return vald_limited


# Convert a VALD linelist (short format) to a format that can be used with SPECTRUM
# - depth_limit: filter out all the lines with a depth less than this percentage (0.05 = 5%)
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def VALD_to_SPECTRUM_format(vald_file, output_file, depth_limit=0.0, data_end=None):
    # Planck constant
    h = 6.626068 * 10e-34 # m^2 kg / s
    # Light speed in vacuum
    c = 299792458.0 # m/s

    # Original VALD linelist
    vald_limited = load_and_filter_VALD(vald_file, depth_limit=depth_limit, data_end=data_end)
    
    # Periodic table
    table = asciitable.read("input/abundances/chemical_elements_symbols.dat", delimiter="\t")

    # Prepare resulting structure
    linelist = np.recarray((len(vald_limited), ), dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'),('transition type', '|S10'), ('note', '|S100')])

    i = 0
    for line in vald_limited:
        linelist[i]['wave (A)'] = line['wave (A)']
        
        element = line['element'].split(" ")
        symbol = element[0]
        try:
            element.remove('') # Make sure there are not additional spaces between the symbol and the ionization state
            element.remove('')
            element.remove('')
        except ValueError:
            pass
        ionization = str(int(element[1]) - 1)
        
        tfilter = (table['symbol'] == symbol)
        linelist[i]['species'] = str(table[tfilter]["atomic_num"][0]) + "." + ionization
        
        #print linelist[i]['species']
        linelist[i]['lower state (cm^-1)'] = int(line['lower state (eV)'] * 8065.73) #cm-1
        # Wavelength
        l = (line[ "wave (A)"] / 10) * 10e-9 # m
        # Frequency
        f = c/l # Hz
        # Energy
        E = h * f # Joules
        E = E * 6.24150974e18 # electron Volt (eV)
        E = E * 8065.73 #cm-1
        linelist[i]['upper state (cm^-1)'] = int(linelist[i]['lower state (cm^-1)'] + E)
        linelist[i]['log(gf)'] = line['log(gf)'] *1.5
        linelist[i]['fudge factor'] = 1.0
        linelist[i]['transition type'] = "99"
        linelist[i]['note'] = line["element"].replace(" ", "_")
        i += 1

    asciitable.write(linelist, output=output_file, Writer=asciitable.FixedWidthNoHeader, delimiter=None, bookend=False, formats={'wave (A)': '%4.3f', })



# Convert a VALD linelist (short format) to a format that can be used to measure radial velocities
# and select the top N deepest lines every wave_step armstrongs (1 nm).
# - data_end can be negative in order to ignore the last nth rows (VALD usually 
#            adds references to the end of the file that should be ignored)
def VALD_top_N_to_RV_format(vald_file, output_file, top=1, wave_step=10, data_end=None):

    vald_limited = load_and_filter_VALD(vald_file, depth_limit=0, data_end=data_end)
    
    wave_base = np.min(vald_limited["wave (A)"])
    wave_top = np.max(vald_limited["wave (A)"])
    
    # Prepare resulting structure
    linelist = np.recarray((top*np.ceil((wave_top - wave_base) / wave_step), ), dtype=[('wave_peak', '<f8'), ('depth', '<f8'), ('element', '|S100')])
    
    wave_current = wave_base
    i = 0
    # For each segment
    while wave_current < wave_top:
        wfilter = (vald_limited["wave (A)"] >= wave_current) & (vald_limited["wave (A)"] < wave_current + wave_step)
        vald_filtered = vald_limited[wfilter]
        vald_filtered.sort(order="depth")
        # Select the top 3 deepest lines
        for j in np.arange(top):
            pos = -1*(j+1)
            linelist[i+j]["wave_peak"] = vald_filtered[pos]["wave (A)"] / 10 #nm
            linelist[i+j]["depth"] = vald_filtered[pos]["depth"]
            linelist[i+j]["element"] = vald_filtered[pos]["element"]
        wave_current += wave_step
        i += top
    
    linelist.sort(order="wave_peak")
    asciitable.write(linelist, output=output_file, delimiter="\t")

if __name__ == '__main__':
    #VALD_to_SPECTRUM_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/linelists/VALD.300_1100nm.lst", depth_limit=0.0)

    VALD_top_3_to_RV_format("input/linelists/original/VALD.300_1100nm_teff_5770.0_logg_4.40.lst", "input/rv/VALD.300_1100nm.rv.lst", top=1, wave_step=10)

