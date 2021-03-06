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
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def plot_spectra(spectra, filename=None, grid = True, title = None, ylabel = 'Flux', xlabel = 'Wavelength (nm)'):
    """
        Plot a spectrum or array of spectra.
        If filename (i.e. plot.png) is specified, then the plot is saved into that file
        and not shown on the screen.
    """
    figure = plt.figure()
    ax1 = plt.subplot(1, 1, 1)

    if grid:
        ax1.grid(True, which="both")

    if title is not None:
        ax1.set_title(title, fontsize="10")

    if xlabel is not None:
        ax1.set_xlabel(xlabel, fontsize="10")

    if ylabel is not None:
        ax1.set_ylabel(ylabel, fontsize="10")

    #~ if filename is not None:
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

    if filename is None:
        plt.show()
    else:
        if filename[-4:] == ".png":
            plt.savefig(filename)
        else:
            plt.savefig(filename+".png")


def show_histogram(x, xlabel='Units', nbins=50):
    """
    Build a histogram from 'x' values and plot it.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # the histogram of the data
    n, bins, patches = ax.hist(x, nbins, normed=0, facecolor='green', alpha=0.75)
    # Mean
    l = plt.axvline(x = np.mean(x), linewidth=1, color='red')
    ax.annotate('Mean', xy=(np.mean(x), np.max(n)),  xycoords='data',
        xytext=(10, 20), textcoords='offset points',
        size=8,
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                        edgecolor='black'),
        horizontalalignment='right', verticalalignment='top',
        )
    # Median
    l = plt.axvline(x = np.median(x), linewidth=1, color='orange')
    ax.annotate('Median', xy=(np.median(x), np.max(n)),  xycoords='data',
        xytext=(10, 35), textcoords='offset points',
        size=8,
        bbox=dict(boxstyle="round", fc="0.8"),
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                        edgecolor='black'),
        horizontalalignment='right', verticalalignment='top',
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Counts')
    ax.grid(True)
    plt.show()


