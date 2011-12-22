#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def plot_spectra(spectra, filename=None, continuum=None, grid = True, title = None, ylabel = 'Flux', xlabel = 'Wavelength (nm)'):
    figure = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    
    if grid:
        ax1.grid(True, which="both")

    if title != None:
        ax1.set_title(title, fontsize="10")

    if xlabel != None:
        ax1.set_xlabel(xlabel, fontsize="10")

    if ylabel != None:
        ax1.set_ylabel(ylabel, fontsize="10")
    
    #~ if filename != None:
        #~ major_tick = np.round((np.max(spectra['waveobs']) - np.min(spectra['waveobs'])) / 100, decimals=2)
        #~ minor_tick = major_tick / 2
        #~ ax1.xaxis.set_major_locator(MultipleLocator(major_tick))
        #~ ax1.xaxis.set_minor_locator(MultipleLocator(minor_tick))
        #~ ax1.xaxis.set_major_formatter(FormatStrFormatter("%1.2f"))
    
    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
    linestyles = ['', 'steps', 'steps:', '-', '--', ':']
    
    if not type(spectra) in [list, np.ndarray]:
        spectra = [spectra]
    
    color = 0
    for spec in spectra:
        ax1.plot(spec['waveobs'], spec['flux'], lw=1, color=colors[color % 7], linestyle=linestyles[3], marker='', markersize=1, markeredgewidth=0, markerfacecolor=colors[0])
        color = (color + 1) % len(colors)
    
    if continuum != None:
        for c in continuum:
            ax1.axvspan(c['wave_base'], c['wave_top'], facecolor='g', alpha=0.5)
    
    if filename == None:
        plt.show()
    else:
        plt.savefig(filename+".png")

