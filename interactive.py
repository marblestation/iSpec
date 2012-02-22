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
#!/usr/bin/env python
#################
# Run with ipython -pdb -c "%run interactive.py"
#################
import ipdb
import os
import pprint
import random
import wx
import asciitable
import numpy as np

import threading
import sys
import getopt

if os.path.exists("synthesizer.so"):
    import synthesizer
    from atmospheres import *


# The recommended way to use wx with mpl is with the WXAgg backend.
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from common import *
from continuum import *
from fitting import *
from radial_velocity import *


class InfoRegionDialog(wx.Dialog):
    def add_parameter(self, name, value, editable=False):
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        param_text = wx.StaticText(self, -1, name, style=wx.ALIGN_LEFT)
        parameter = wx.TextCtrl(self, -1, str(value),  style=wx.TE_RIGHT)
        parameter.SetEditable(editable)
        
        hbox.AddSpacer(10)
        hbox.Add(param_text, 0, border=3, flag=wx.ALIGN_LEFT | wx.ALL | wx.EXPAND | wx.ALIGN_CENTER_VERTICAL)
        hbox.Add(parameter, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        return hbox, parameter
        
    def __init__(self, parent, id, title, spectra, region):
        wx.Dialog.__init__(self, parent, id, title)
        self.region = region
        
        self.update = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Info
        wave_base = region.get_wave_base()
        wave_top = region.get_wave_top()
        wave_filter = (spectra['waveobs'] >= wave_base) & (spectra['waveobs'] <= wave_top)
        spectra_window = spectra[wave_filter]
        std = np.std(spectra_window['flux'])
        num = len(spectra_window['flux'])
        
        hbox, self.element_type = self.add_parameter("Type:", region.element_type)
        self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        hbox, self.wave_base = self.add_parameter("Min. wavelength:", wave_base)
        self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        hbox, self.wave_top = self.add_parameter("Max. wavelength:", wave_top)
        self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        hbox, self.wave_top = self.add_parameter("Flux standard deviation:", std)
        self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        hbox, self.wave_top = self.add_parameter("Number of points:", num)
        self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        if region.element_type == "lines":
            note = region.get_note_text()
            wave_peak = region.get_wave_peak()
            hbox, self.wave_top = self.add_parameter("Wave peak:", wave_peak)
            self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
            hbox, self.wave_top = self.add_parameter("Note:", note)
            self.vbox.Add(hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.update = True
        self.Close()

class FitContinuumDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_nknots = wx.StaticText(self, -1, "Number of knots for a uniform spline fit: ", style=wx.ALIGN_LEFT)
        self.nknots = wx.TextCtrl(self, -1, '20',  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_nknots, 0, border=3, flag=flags)
        self.hbox.Add(self.nknots, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Fit using: ", style=wx.ALIGN_LEFT)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_where, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.radio_button_spectra = wx.RadioButton(self, -1, 'The whole spectra', style=wx.RB_GROUP)
        self.vbox2.Add(self.radio_button_spectra, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_continuum = wx.RadioButton(self, -1, 'Only continuum regions')
        self.vbox2.Add(self.radio_button_continuum, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_spectra.SetValue(True)
        self.hbox.Add(self.vbox2, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()


class FindContinuumDialog(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Max step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_fixed_wave_step = wx.StaticText(self, -1, "Check for regions of minimum size: ", style=wx.ALIGN_LEFT)
        self.fixed_wave_step = wx.TextCtrl(self, -1, '0.01',  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_fixed_wave_step, 0, border=3, flag=flags)
        self.hbox.Add(self.fixed_wave_step, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_sigma = wx.StaticText(self, -1, "Maximum standard deviation: ", style=wx.ALIGN_LEFT)
        self.sigma = wx.TextCtrl(self, -1, '0.001',  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_sigma, 0, border=3, flag=flags)
        self.hbox.Add(self.sigma, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Where to look
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Look for continuum regions in: ", style=wx.ALIGN_LEFT)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_where, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.radio_button_spectra = wx.RadioButton(self, -1, 'The whole spectra', style=wx.RB_GROUP)
        self.vbox2.Add(self.radio_button_spectra, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_segments = wx.RadioButton(self, -1, 'Only inside segments')
        self.vbox2.Add(self.radio_button_segments, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_spectra.SetValue(True)
        self.hbox.Add(self.vbox2, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()

class CorrectRVDialog(wx.Dialog):
    def __init__(self, parent, id, title, rv):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### RV
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_rv = wx.StaticText(self, -1, "Radial velocity (km/s): ", style=wx.ALIGN_LEFT)
        self.rv = wx.TextCtrl(self, -1, str(rv),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv, 0, border=3, flag=flags)
        self.hbox.Add(self.rv, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Where to look
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Apply correction on: ", style=wx.ALIGN_LEFT)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_where, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.radio_button_spectra = wx.RadioButton(self, -1, 'Spectra', style=wx.RB_GROUP)
        self.vbox2.Add(self.radio_button_spectra, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_regions = wx.RadioButton(self, -1, 'Regions')
        self.vbox2.Add(self.radio_button_regions, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_spectra.SetValue(True)
        self.hbox.Add(self.vbox2, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()

class DetermineRVDialog(wx.Dialog):
    def __init__(self, parent, id, title, rv_limit, rv_step):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### RV limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_rv_limit = wx.StaticText(self, -1, "Radial velocity limit abs(km/s): ", style=wx.ALIGN_LEFT)
        self.rv_limit = wx.TextCtrl(self, -1, str(rv_limit),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_limit, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_limit, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### RV step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_rv_step = wx.StaticText(self, -1, "Radial velocity steps (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_step = wx.TextCtrl(self, -1, str(rv_step),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_step, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_step, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()


class RVProfileDialog(wx.Dialog):
    def __init__(self, parent, id, title, xcoord, fluxes, model):
        wx.Dialog.__init__(self, parent, id, title, size=(600, 600))
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Plot
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((3.0, 3.0), dpi=self.dpi)
        self.canvas = FigCanvas(self, -1, self.fig)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        self.axes = self.fig.add_subplot(1, 1, 1)
        
        self.toolbar = NavigationToolbar(self.canvas)
        
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        
        ### Stats
        self.stats = wx.ListCtrl(self, -1, style=wx.LC_REPORT)
        self.stats.InsertColumn(0, 'Property')
        self.stats.InsertColumn(1, 'Value')
        self.stats.SetColumnWidth(0, 300)
        self.stats.SetColumnWidth(1, 300)
        
        #self.vbox.Add(self.stats, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP | wx.EXPAND)
        self.vbox.Add(self.stats, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        
        ### RV limit
        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        #self.text_rv_limit = wx.StaticText(self, -1, "Radial velocity limit abs(km/s): ", style=wx.ALIGN_LEFT)
        #self.rv_limit = wx.TextCtrl(self, -1, str(rv_limit),  style=wx.TE_RIGHT)
        
        #self.hbox.AddSpacer(10)
        #self.hbox.Add(self.text_rv_limit, 0, border=3, flag=flags)
        #self.hbox.Add(self.rv_limit, 0, border=3, flag=flags)
        
        #self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        sizer =  self.CreateButtonSizer(wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        
        ## Draw
        self.axes.plot(xcoord, fluxes, lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)
        self.axes.plot(xcoord, model(xcoord), lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)
        
        self.axes.grid(True, which="both")
        self.axes.set_title("Profile", fontsize="10")
        self.axes.set_xlabel("radial velocity (km/s)", fontsize="10")
        self.axes.set_ylabel("average flux", fontsize="10")
        
        ## Stats
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, "Mean (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(model.mu, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Sigma (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(model.sig, 2)))
        #num_items += 1
        #self.stats.InsertStringItem(num_items, "A")
        #self.stats.SetStringItem(num_items, 1, str(np.round(model.A, 2)))
        fwhm = model.sig * (2*np.sqrt(2*np.log(2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "FWHM (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(fwhm, 2)))

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()

class SyntheticSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top, resolution, teff, logg, MH, microturbulence_vel):
        wx.Dialog.__init__(self, parent, id, title, size=(450,450))
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        
        ### Teff
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_teff = wx.StaticText(self, -1, "Effective temperature (K): ", style=wx.ALIGN_LEFT)
        self.teff = wx.TextCtrl(self, -1, str(teff),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_teff, 0, border=3, flag=flags)
        self.hbox.Add(self.teff, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Log g
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_logg = wx.StaticText(self, -1, "Gravity (log g): ", style=wx.ALIGN_LEFT)
        self.logg = wx.TextCtrl(self, -1, str(logg),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_logg, 0, border=3, flag=flags)
        self.hbox.Add(self.logg, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### MH
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_MH = wx.StaticText(self, -1, "Metallicity [M/H]: ", style=wx.ALIGN_LEFT)
        self.MH = wx.TextCtrl(self, -1, str(MH),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_MH, 0, border=3, flag=flags)
        self.hbox.Add(self.MH, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### MH
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_microturbulence_vel = wx.StaticText(self, -1, "Microturbulence velocity (km/s): ", style=wx.ALIGN_LEFT)
        self.microturbulence_vel = wx.TextCtrl(self, -1, str(microturbulence_vel),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_microturbulence_vel, 0, border=3, flag=flags)
        self.hbox.Add(self.microturbulence_vel, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        
        ### Resolution
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_resolution = wx.StaticText(self, -1, "Resolution: ", style=wx.ALIGN_LEFT)
        self.resolution = wx.TextCtrl(self, -1, str(resolution),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_resolution, 0, border=3, flag=flags)
        self.hbox.Add(self.resolution, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Wave base
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_base = wx.StaticText(self, -1, "Wavelength min (nm): ", style=wx.ALIGN_LEFT)
        self.wave_base = wx.TextCtrl(self, -1, str(wave_base),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_base, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_base, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Wave top
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_top = wx.StaticText(self, -1, "Wavelength max (nm): ", style=wx.ALIGN_LEFT)
        self.wave_top = wx.TextCtrl(self, -1, str(wave_top),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_top, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_top, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        self.vbox.AddSpacer(20)
        
        ### Wavelength range
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Generate spectrum for: ", style=wx.ALIGN_LEFT)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_where, 0, border=3, flag=flags)
        
        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.radio_button_custom_range = wx.RadioButton(self, -1, 'Custom range (defined above)', style=wx.RB_GROUP)
        self.vbox2.Add(self.radio_button_custom_range, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_segments = wx.RadioButton(self, -1, 'Segments')
        self.vbox2.Add(self.radio_button_segments, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_custom_range.SetValue(True)
        self.hbox.Add(self.vbox2, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(30)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()


class CustomizableRegion:
    min_width = 0.002 # nm
    mark = None
    
    def __init__(self, frame, element_type, axvspan, mark=None, note=None):
        self.frame = frame
        self.axvspan = axvspan
        self.element_type = element_type
        self.press = None
        self.original_facecolor = self.axvspan.get_facecolor()
        self.original_edgecolor = self.axvspan.get_edgecolor()
        ## Line region specific properties:
        self.mark = mark
        self.note = note
        # Fit line properties, dictionary of spectrum (to vinculate to different spectra):
        self.line_plot_id = {}
        self.line_model = {}
        self.continuum_base_level = {}
    
    def connect(self):
        # Connect to all the events
        self.cid_press = self.axvspan.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.axvspan.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.axvspan.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
    
    
    def update_size(self, event):
        # Update the position of the edge left or right
        xy = self.axvspan.get_xy()
        
        # If left button clicked, modify left edge
        if event.button == 1:
            # Check condition to allow modification
            if self.mark != None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata < self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True
            
            # Do not modify if the region will become slimmer than...
            big_enough = xy[2,0] - event.xdata > self.min_width
            
            if big_enough and compatible_with_mark:
                xy[0,0] = event.xdata
                xy[1,0] = event.xdata
                xy[4,0] = event.xdata
                self.frame.status_message("Moving left edge to %.4f" % xy[0,0])
                self.frame.canvas.draw()
        elif event.button == 3:
            # Check condition to allow modification
            if self.mark != None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata > self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True
            
            # Do not modify if the region will become slimmer than...
            big_enough = event.xdata - xy[0,0] > self.min_width
            
            if big_enough and compatible_with_mark:
                xy[2,0] = event.xdata
                xy[3,0] = event.xdata
                self.frame.status_message("Moving right edge to %.4f" % xy[2,0])
                self.frame.canvas.draw()
        
    def update_mark(self, event):
        # Position of the edge left or right
        xy = self.axvspan.get_xy()
        
        # If left button clicked, modify mark
        if event.button == 1 and self.mark != None:
            x = self.mark.get_xdata()
            
            inside_region = (event.xdata > xy[0,0]) and (event.xdata < xy[2,0]) 
            if inside_region:
                x[0] = event.xdata
                x[1] = event.xdata
                self.mark.set_xdata(x)
                if self.note != None:
                    self.note.xy = (event.xdata, 1)
                self.frame.status_message("Moving mark to %.4f" % x[0])
                self.frame.canvas.draw()
    
    def update_mark_note(self, event):
        # Position of the mark
        x = self.mark.get_xdata()
        
        if self.note != None:
            note_text = self.note.get_text()
        else:
            note_text = ""
        
        note_text = self.frame.ask('Note for the new line region:', 'Note', note_text)
        
        if note_text != None and note_text != "":
            # Set
            if self.note == None:
                # New
                self.note = self.frame.axes.annotate(note_text, xy=(x[0], 1),  xycoords='data',
                    xytext=(-10, 20), textcoords='offset points',
                    size=8,
                    bbox=dict(boxstyle="round", fc="0.8"),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                    edgecolor='black'),
                    horizontalalignment='right', verticalalignment='top',
                    )
            else:
                # Update
                self.note.set_text(note_text)
                self.frame.canvas.draw()
        elif note_text == "" and self.note != None:
            # Remove
            self.note.set_visible(False)
            self.note = None
            self.frame.canvas.draw()
        
        
    def on_press(self, event):
        # Validate that the click is on this axis/object
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        contains, attrd = self.axvspan.contains(event)
        if not contains: return
        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return
        
        # When regions overlap two or more can receive the click event, so
        # let's use a lock to allow the modification of one of them
        if self.frame.lock.acquire(False):
            if self.frame.action == "Remove":
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Remove note
                    if self.note != None:
                        self.note.set_visible(False)
                        self.note = None
                        self.frame.canvas.draw()
                else:
                    ## Remove region
                    self.disconnect_and_hide()
                    self.frame.region_widgets[self.element_type].remove(self)
                    self.frame.regions_changed(self.element_type)
                    self.frame.flash_status_message("Removed region from " + "%.4f" % self.axvspan.get_xy()[0,0] + " to " + "%.4f" % self.axvspan.get_xy()[2,0])
                    self.frame.canvas.draw()
                    self.press = event.button, event.x, event.xdata
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Modify":
                ## Modify region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata
                
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    if event.button == 1:
                        # Left button, modify mark position
                        self.update_mark(event)
                    elif event.button == 2:
                        pass
                        # Modification of the note is managed on_release
                        # because we cannot show modal windows asking information
                        # while on_press or the mouse stops working
                else:
                    self.update_size(event)
                self.frame.regions_changed(self.element_type)
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                ## Create note but handel it on_release
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata
            elif self.frame.action == "Stats":
                ## Statistics about the region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata
                
            
    def disconnect_and_hide(self):
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
        self.axvspan.set_visible(False)
        if self.mark != None:
            self.mark.set_visible(False)
        if self.note != None:
            self.note.set_visible(False)
        if self.line_plot_id != None:
            for plot_id in self.line_plot_id.values():
                if v != None:
                    self.frame.axes.lines.remove(plot_id)
        
    
    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        
        # If it is in remove mode, the lock should be released
        # NOTE: The selected region has already been removed, so it will not  
        # receive the "release" event.
        if self.frame.action == "Remove" and self.frame.lock.locked():
            self.frame.lock.release()
        elif self.press != None:
            # In modification mode, if it is the current selected widget
            self.press = None
            if self.frame.action == "Modify":
                # If it was a right click when elements is lines and subelements marks
                # => Change note
                if event.button == 3 and self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Right button, modify note
                    if self.frame.action == "Modify" or self.frame.action == "Create":
                        self.update_mark_note(event)
                    else:
                        # Note removal is managed on_press
                        pass
            elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                # Create a note (left or right click)
                if (event.button == 1 or event.button == 3) and self.note == None:
                    self.update_mark_note(event)
            elif self.frame.action == "Stats" and self.frame.lock.locked():
                self.frame.update_stats(self)
                #~ dlg = InfoRegionDialog(self.frame, -1, "Info", spectra, self)
                #~ dlg.ShowModal()
                #~ if not dlg.action_accepted:
                #~ dlg.Destroy()
                #~ resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
            
            self.axvspan.set_facecolor(self.original_facecolor)
            self.axvspan.set_edgecolor(self.original_edgecolor)
            self.axvspan.set_alpha(0.30)
            self.frame.status_message("")
            self.frame.lock.release()
            self.frame.canvas.draw()
    
    def on_motion(self, event):
        # Validate that the object has been pressed and the click is on this axis
        if self.press is None: return
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes.get_navigate_mode() != None: return
        if self.frame.action == "Stats": return
        
#        button, x, xpress = self.press
        if self.frame.subelements == "marks":
            self.update_mark(event)
        else:
            self.update_size(event)
            
            
    def get_note_text(self):
        note_text = ""
        if self.element_type == "lines" and self.note != None:
            note_text = self.note.get_text()
        return note_text
    
    def get_wave_peak(self):
        if self.element_type == "lines":
            return self.mark.get_xdata()[0]
        else:
            return None
    
    def get_wave_base(self):
        return self.axvspan.get_xy()[0,0]
    
    def get_wave_top(self):
        return self.axvspan.get_xy()[2,0]
    

class Spectrum():
    def __init__(self, data, name, color='b', not_saved=False, plot_id=None, continuum_model=None, continuum_data=None, continuum_plot_id=None, rv=0, path=None):
        self.data = data
        self.name = name
        self.color = color
        self.not_saved = not_saved
        self.plot_id = plot_id
        self.continuum_model = continuum_model
        self.continuum_data = continuum_data
        self.continuum_plot_id = continuum_plot_id
        self.rv = rv # Radial velocity (km/s)
        self.path = path


class SpectraFrame(wx.Frame):
    title = 'Spectra Visual Editor '
    
    def ipython(self):
        import IPython
        self.embedshell = IPython.Shell.IPShellEmbed(argv=['--pdb'])
        self.embedshell()
    
    # regions should be a dictionary with 'continuum', 'lines' and 'segments' keys
    # filenames should be a dictionary with 'spectra', 'continuum', 'lines' and 'segments' keys
    def __init__(self, spectra=None, name=None, regions=None, filenames=None):
        
        self.embedshell = None
        self.ipython_thread = None
        #self.ipython_thread = threading.Thread(target=self.ipython)
        #self.ipython_thread.setDaemon(True)
        #self.ipython_thread.start()
        
        self.spectra_color = ('#3100FF', '#0074FF', '#E49400', '#8b6914', '#00DEFF', '#02AA4D', '#49CD2E', '#A4DC3A', '#FF69D2',)
        self.spectra_color_id = 0
        
        if spectra == None:
            self.spectra = []
            self.active_spectrum = None
        else:
            self.active_spectrum = Spectrum(spectra, name, path=filenames["spectra"], color=self.spectra_color[self.spectra_color_id])
            self.spectra_color_id = (self.spectra_color_id + 1) % len(self.spectra_color)
            self.spectra = [self.active_spectrum]
        
        
        if regions == None:
            continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
            lines = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
            segments = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])

            self.regions = {}
            self.regions['continuum'] = continuum
            self.regions['lines'] = lines
            self.regions['segments'] = segments
        else:
            self.regions = regions
        
        if filenames == None:
            self.filenames = {}
            self.filenames['continuum'] = None
            self.filenames['lines'] = None
            self.filenames['segments'] = None
        else:
            self.filenames = filenames
        
        self.region_widgets = {} # It is important not to lose the reference to the regions or the garbage 
        self.region_widgets['continuum'] = []
        self.region_widgets['lines'] = []
        self.region_widgets['segments'] = []
        
        self.action = "Stats"
        self.elements = "continuum"
        self.subelements = None
        self.lock = threading.Lock()
        self.not_saved = {}
        self.not_saved["continuum"] = False
        self.not_saved["lines"] = False
        self.not_saved["segments"] = False
        self.rv_limit = 200 # km/s
        self.rv_step = 0.5 # km/s
        self.modeled_layers_pack = None # Synthesize spectrum (atmospheric models)
        
        self.timeroff = None
        
        wx.Frame.__init__(self, None, -1, self.title)
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.draw_figure()
    
    def ask(self, text, title, default_value):
            dlg = wx.TextEntryDialog(self, text, title)
            
            dlg.SetValue(default_value)
            if dlg.ShowModal() == wx.ID_OK:
                response = dlg.GetValue()
            else:
                response = None
                
            dlg.Destroy()
            return response
    
    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Open spectra\tCtrl-O", "Open spectra file")
        self.Bind(wx.EVT_MENU, self.on_open_spectra, m_expt)
        m_expt = menu_file.Append(-1, "Open co&ntinuum regions\tCtrl-N", "Open continuum regions file")
        self.Bind(wx.EVT_MENU, self.on_open_continuum, m_expt)
        m_expt = menu_file.Append(-1, "Open l&ine regions\tCtrl-I", "Open line regions file")
        self.Bind(wx.EVT_MENU, self.on_open_lines, m_expt)
        m_expt = menu_file.Append(-1, "Open se&gments\tCtrl-G", "Open segments file")
        self.Bind(wx.EVT_MENU, self.on_open_segments, m_expt)
        menu_file.AppendSeparator()
        m_expt = menu_file.Append(-1, "&Save plot image as\tCtrl-S", "Save plot image to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        m_expt = menu_file.Append(-1, "Save s&pectra as\tCtrl-P", "Save active spectra to file")
        self.Bind(wx.EVT_MENU, self.on_save_spectra, m_expt)
        m_expt = menu_file.Append(-1, "Save &continuum regions as\tCtrl-C", "Save continuum regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_continuum_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save &line regions as\tCtrl-L", "Save line regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_line_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save s&egments as\tCtrl-E", "Save segments to file")
        self.Bind(wx.EVT_MENU, self.on_save_segments, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_edit = wx.Menu()
        self.menu_active_spectrum = wx.Menu()
        self.update_menu_active_spectrum()
        menu_edit.AppendMenu(-1, 'Select spectra', self.menu_active_spectrum)
        
        m_fit_continuum = menu_edit.Append(-1, "Fit &continuum", "Fit continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_continuum, m_fit_continuum)
        m_fit_lines = menu_edit.Append(-1, "Fit &lines", "Fit lines using the fitted continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_lines, m_fit_lines)
        m_find_continuum = menu_edit.Append(-1, "&Find continuum regions", "Find continuum regions")
        self.Bind(wx.EVT_MENU, self.on_find_continuum, m_find_continuum)
        m_determine_rv = menu_edit.Append(-1, "&Determine RV", "Determine radial velocity")
        self.Bind(wx.EVT_MENU, self.on_determine_rv, m_determine_rv)
        m_correct_rv = menu_edit.Append(-1, "&Correct RV", "Correct spectra by using its radial velocity")
        self.Bind(wx.EVT_MENU, self.on_correct_rv, m_correct_rv)
        
        if sys.modules.has_key('synthesizer'):
            m_synthesize = menu_edit.Append(-1, "&Synthesize spectrum", "Synthesize spectrum")
            self.Bind(wx.EVT_MENU, self.on_synthesize, m_synthesize)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the visual editor")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(1, 1, 1)
        
        # Make space for the legend
        #box = self.axes.get_position()
        #self.axes.set_position([box.x0, box.y0, box.width * 0.80, box.height])
        
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
                
        text_action = wx.StaticText(self.panel, -1, "Action: ", style=wx.ALIGN_LEFT)
        self.radio_button_stats = wx.RadioButton(self.panel, -1, 'Stats', style=wx.RB_GROUP)
        self.radio_button_create = wx.RadioButton(self.panel, -1, 'Create')
        self.radio_button_modify = wx.RadioButton(self.panel, -1, 'Modify')
        self.radio_button_remove = wx.RadioButton(self.panel, -1, 'Remove')
        self.radio_button_stats.SetValue(True)
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_create.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_modify.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_remove.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_stats.GetId())
        
        text_element = wx.StaticText(self.panel, -1, "Elements: ", style=wx.ALIGN_LEFT)
        self.radio_button_continuum = wx.RadioButton(self.panel, -1, 'Continuum', style=wx.RB_GROUP)
        self.radio_button_lines = wx.RadioButton(self.panel, -1, 'Lines')
        self.radio_button_linemarks = wx.RadioButton(self.panel, -1, 'Line marks')
        self.radio_button_segments = wx.RadioButton(self.panel, -1, 'Segments')
        self.radio_button_continuum.SetValue(True)
        self.radio_button_linemarks.Enable(False) # Because "Stats" is selected as an action by default
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_continuum.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_lines.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_linemarks.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_segments.GetId())

        # Progress bar
        self.gauge = wx.Gauge(self.panel, range=100)
        
        self.stats = wx.ListCtrl(self.panel, -1, style=wx.LC_REPORT)
        self.stats.InsertColumn(0, 'Property')
        self.stats.InsertColumn(1, 'Value')
        self.stats.SetColumnWidth(0, 300)
        self.stats.SetColumnWidth(1, 300)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)
        
        #
        # Layout with box sizers
        #
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)
        
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        self.hbox.Add(text_action, 0, border=3, flag=flags)
        
        self.hbox.Add(self.radio_button_stats, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_create, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_modify, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_remove, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)
        
        self.hbox.Add(text_element, 0, border=3, flag=flags)
        
        self.hbox.Add(self.radio_button_continuum, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_lines, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_linemarks, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_segments, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)
        
        self.parent_hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.parent_hbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP)
        
        # Progress bar
        self.parent_hbox.Add(self.gauge, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        
        self.vbox.Add(self.parent_hbox, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP)
        
        self.vbox.Add(self.stats, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP | wx.EXPAND)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
    
   
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_active_spectrum(self):
        if self.active_spectrum == None:
            return
        
        # Remove spectrum plot if exists
        if self.active_spectrum.plot_id != None:
            self.axes.lines.remove(self.active_spectrum.plot_id)
        
        # zorder = 1, always in the background
        self.active_spectrum.plot_id = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'], lw=1, color=self.active_spectrum.color, linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor=self.active_spectrum.color, zorder=1, label="[A] "+self.active_spectrum.name)[0]
        
        # Put a legend to the right of the current axis
        self.update_legend()
        
    def remove_drawn_continuum_spectra(self):
        if self.active_spectrum != None and self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
            self.active_spectrum.continuum_plot_id = None
            self.canvas.draw()
    
    def draw_continuum_spectra(self):
        if self.active_spectrum == None:
            return
        
        # Remove continuum plot if exists
        if self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
        
        # zorder = 1, always in the background
        if self.active_spectrum.continuum_data != None:
            self.active_spectrum.continuum_plot_id = self.axes.plot(self.active_spectrum.continuum_data['waveobs'], self.active_spectrum.continuum_data['flux'], lw=1, color='green', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]
    
    # Draw all elements in regions array and creates widgets in regions widget array
    # Elements cann be "continuum", "lines" or "segments"
    def draw_regions(self, elements):
        for r in self.region_widgets[elements]:
            r.disconnect_and_hide()
        self.region_widgets[elements] = []
        
        if elements == "continuum":
            color = 'green'
        elif elements == "lines":
            color = 'yellow'
        else:
            color = 'grey'
        
        for r in self.regions[elements]:
            axvs = self.axes.axvspan(r['wave_base'], r['wave_top'], facecolor=color, alpha=0.30)
            axvs.zorder = 3
            
            # Segments always in the background but above spectra
            if elements == "segments":
                axvs.zorder = 2
            
            if elements == "lines":
                # http://matplotlib.sourceforge.net/examples/pylab_examples/annotation_demo.html
                # http://matplotlib.sourceforge.net/examples/pylab_examples/annotation_demo2.html
                if r['note'] != None and r['note'] != "":
                    note = self.axes.annotate(r['note'], xy=(r['wave_peak'], 1),  xycoords='data',
                        xytext=(-10, 20), textcoords='offset points',
                        size=8,
                        bbox=dict(boxstyle="round", fc="0.8"),
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                        edgecolor='black'),
                        horizontalalignment='right', verticalalignment='top',
                        )
                else:
                    note = None
                axvline = self.axes.axvline(x = r['wave_peak'], linewidth=1, color='orange')
                region = CustomizableRegion(self, "lines", axvs, mark=axvline, note=note)
            else:
                region = CustomizableRegion(self, elements, axvs)
            region.connect()
            self.region_widgets[elements].append(region)
    
    def draw_figure(self):
        """ Redraws the figure
        """
        self.axes.grid(True, which="both")
        self.axes.set_title("Spectra", fontsize="10")
        self.axes.set_xlabel("wavelength", fontsize="10")
        self.axes.set_ylabel("flux", fontsize="10")
        
        self.draw_active_spectrum()
        self.draw_regions("continuum")
        self.draw_regions("lines")
        self.draw_regions("segments")
        self.canvas.draw()
    
    def update_title(self):
        title = self.title
        
        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"] or spectra_not_saved:
            title += "("
            if self.not_saved["continuum"]:
                title += "*continuum"
            if self.not_saved["lines"]:
                title += "*lines"
            if self.not_saved["segments"]:
                title += "*segments"
            for spec in self.spectra:
                if spec.not_saved:
                    title += "*" + spec.name
            title += "* not saved)"
        self.SetTitle(title)
    
    # Change saved status and title for the concrete element
    def regions_changed(self, element_type):
        self.not_saved[element_type] = True
        self.update_title()
    
    # Change saved status and title for the concrete element
    def regions_saved(self, element_type):
        self.not_saved[element_type] = False
        self.update_title()
        
    
    def on_close(self, event):
        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"] or spectra_not_saved:
            dlg = wx.MessageDialog(self, 'Are you sure you want to exit without saving the regions/spectra?', 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                self.Destroy()
        else:
            self.Destroy()
        
        # Embeded ipython terminal
        if self.ipython_thread != None:
            self.ipython_thread.join()
        

    def on_action_change(self, event):
        self.enable_elements()
        self.stats.DeleteAllItems()
        
        if self.radio_button_create.GetId() == event.GetId():
            self.action = "Create"
        elif self.radio_button_modify.GetId() == event.GetId():
            self.action = "Modify"
        elif self.radio_button_remove.GetId() == event.GetId():
            self.action = "Remove"
        else:
            self.action = "Stats"
            self.radio_button_linemarks.Enable(False)
            if self.radio_button_linemarks.GetValue():
                self.radio_button_continuum.SetValue(True)
                self.elements = "continuum"
            
    def disable_elements(self):
        self.radio_button_continuum.Enable(False)
        self.radio_button_lines.Enable(False)
        self.radio_button_linemarks.Enable(False)
        self.radio_button_segments.Enable(False)
    
    def enable_elements(self):
        self.radio_button_continuum.Enable(True)
        self.radio_button_lines.Enable(True)
        self.radio_button_linemarks.Enable(True)
        self.radio_button_segments.Enable(True)
    
    def on_element_change(self, event):
        if self.radio_button_lines.GetId() == event.GetId():
            self.elements = "lines"
            self.subelements = None
        elif self.radio_button_linemarks.GetId() == event.GetId():
            self.elements = "lines"
            self.subelements = "marks"
        elif self.radio_button_segments.GetId() == event.GetId():
            self.elements = "segments"
            self.subelements = None
        else:
            self.elements = "continuum"
            self.subelements = None
    
    def on_motion(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        
        self.status_message("Cursor on wavelength %.4f" % event.xdata + " and flux %.4f" % event.ydata)

    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        
        new_halfwidth = 0.05
        # If the left button is clicked when action "create" is active, 
        # create a new region of the selected element type
        # NOTE: line marks only can be created from a line region
        if self.action == "Create" and self.subelements != "marks" and event.button == 1 and event.key == None:
            if self.elements == "continuum":
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='green', alpha=0.30)
                axvs.zorder = 3
                region = CustomizableRegion(self, "continuum", axvs)
                region.connect()
                self.region_widgets['continuum'].append(region)
            elif self.elements == "lines":
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='yellow', alpha=0.30)
                axvs.zorder = 3
                axvline = self.axes.axvline(x = event.xdata, linewidth=1, color='orange')
                
                dlg = wx.TextEntryDialog(self, 'Note for the new line region:','Note')
                dlg.SetValue("")
                if dlg.ShowModal() == wx.ID_OK:
                    note_text = dlg.GetValue()
                    if note_text != "":
                        note = self.axes.annotate(note_text, xy=(event.xdata, 1),  xycoords='data',
                            xytext=(-10, 20), textcoords='offset points',
                            size=8,
                            bbox=dict(boxstyle="round", fc="0.8"),
                            arrowprops=dict(arrowstyle="->",
                                            connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                            edgecolor='black'),
                            horizontalalignment='right', verticalalignment='top',
                            )
                    else:
                        note = None
                else:
                    note = None
                dlg.Destroy()
                
                region = CustomizableRegion(self, "lines", axvs, mark=axvline, note=note)
                region.connect()
                self.region_widgets['lines'].append(region)
            elif self.elements == "segments":
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='grey', alpha=0.30)
                # Segments always in the background but above spectra
                axvs.zorder = 2
                region = CustomizableRegion(self, "segments", axvs)
                region.connect()
                self.region_widgets['segments'].append(region)
            self.canvas.draw()
            self.regions_changed(self.elements)
            self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - new_halfwidth) + " to " + "%.4f" % (event.xdata + new_halfwidth))
            ## Automatically change from create to modify, because this will be
            ## the general will of the user
            #~ self.radio_button_create.SetValue(False)
            #~ self.radio_button_modify.SetValue(True)
            #~ self.action = "Modify"
    
    def update_progress(self, value):
        wx.CallAfter(self.gauge.SetValue, pos=value)

    def update_menu_active_spectrum(self):
        # Remove everything
        for i in np.arange(self.menu_active_spectrum.MenuItemCount):
            self.menu_active_spectrum.Remove(i)
        
        # Add as many options as spectra
        for i in np.arange(len(self.spectra)):
            if self.spectra[i] == None:
                continue
            
            spec_element = self.menu_active_spectrum.Append(i, self.spectra[i].name, kind=wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.on_change_active_spectrum, spec_element)

            if self.active_spectrum == self.spectra[i]:
                spec_element.Check()
    
    def update_legend(self):
        leg = self.axes.legend(loc='center right', bbox_to_anchor=(1.1, 1), ncol=1, shadow=False)
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='8')
    
    def add_stats(self, k, v):
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, k)
        self.stats.SetStringItem(num_items, 1, str(v))
        
    def update_stats(self, region):
        self.stats.DeleteAllItems()
        
        wave_base = region.get_wave_base()
        wave_top = region.get_wave_top()
        
        wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
        spectra_window = self.active_spectrum.data[wave_filter]
        num_points = len(spectra_window['flux'])
        
        self.add_stats("Number of measures", num_points)
        
        self.add_stats("Wavelength min.", "%.4f" % wave_base)
        self.add_stats("Wavelength max.", "%.4f" % wave_top)
        wave_diff = wave_top - wave_base
        self.add_stats("Wavelength difference", "%.4f" % wave_diff)
        
        if region.element_type == "lines":
            self.add_stats("Wavelength peak", "%.4f" % region.get_wave_peak())
            note = region.get_note_text()
            if note != "":
                self.add_stats("Note", note)
        
        if num_points > 0:
            self.add_stats("Flux min.", "%.6f" % np.min(spectra_window['flux']))
            self.add_stats("Flux max.", "%.6f" % np.max(spectra_window['flux']))
            self.add_stats("Flux mean", "%.6f" % np.mean(spectra_window['flux']))
            self.add_stats("Flux median", "%.6f" % np.median(spectra_window['flux']))
            self.add_stats("Flux standard deviation", "%.6f" % np.std(spectra_window['flux']))
        
        if region.element_type == "lines" and region.line_model != None:
            self.add_stats("Gaussian mean (mu)", "%.4f" % region.line_model.mu)
            self.add_stats("Gaussian width (A)", "%.4f" % region.line_model.A)
            self.add_stats("Gaussian standard deviation (sigma)", "%.4f" % region.line_model.sig)
            self.add_stats("Gaussian base level (mean continuum)", "%.4f" % region.continuum_base_level)
            rms = np.mean(region.line_model.residuals()) + np.std(region.line_model.residuals())
            self.add_stats("Gaussian fit root mean squeare (RMS)", "%.4f" % rms)
        
        if self.active_spectrum.continuum_model != None:
            if num_points > 0:
                mean_continuum = np.mean(self.active_spectrum.continuum_model(spectra_window['waveobs']))
                self.add_stats("Continuum mean for the region", "%.4f" % mean_continuum)
            rms = np.mean(self.active_spectrum.continuum_model.residuals()) + np.std(self.active_spectrum.continuum_model.residuals())
            self.add_stats("Continuum fit root mean square (RMS)", "%.4f" % rms)
            
        
    
    def open_file(self, elements):
        file_choices = "All|*"
        
        if elements != "spectra" and self.not_saved[elements]:
            dlg = wx.MessageDialog(self, 'Are you sure you want to open a new %s file without saving the current regions?' % elements, 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() != wx.ID_YES:
                self.flash_status_message("Discarded.")
                return
        
        
        if elements == "spectra":
            if self.active_spectrum != None and self.active_spectrum.path != None:
                filename = self.active_spectrum.path.split('/')[-1]
                filename_length = len(filename)
                dirname = self.active_spectrum.path[:-filename_length]
            elif self.active_spectrum != None:
                filename = self.active_spectrum.name
                dirname = os.getcwd()
            else:
                filename = ""
                dirname = os.getcwd()
        else:
            if self.filenames[elements] != None:
                filename = self.filenames[elements].split('/')[-1]
                filename_length = len(filename)
                dirname = self.filenames[elements][:-filename_length]
            else:
                filename = ""
                dirname = os.getcwd()
        # create a new thread when a button is pressed
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Open %s..." % elements,
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.OPEN)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if not os.path.exists(path):
                    dlg_error = wx.MessageDialog(self, 'File %s does not exist.' % dlg.GetFilename(), 'File does not exist', wx.OK|wx.ICON_ERROR)
                    dlg_error.ShowModal()
                    dlg_error.Destroy()
                    continue # Give the oportunity to select another file name
                
                try:
                    if elements == "spectra":
                        new_spectra_data = read_spectra(path)
                        
                        # Remove current continuum from plot if exists
                        self.remove_drawn_continuum_spectra()
                        
                        # Remove current drawn fitted lines if they exist
                        self.remove_drawn_fitted_lines()
                        
                        # Remove "[A]  " from spectra name (legend) if it exists
                        if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
                        
                        self.active_spectrum = Spectrum(new_spectra_data, path.split('/')[-1], path = path, color=self.spectra_color[self.spectra_color_id])
                        self.spectra_color_id = (self.spectra_color_id + 1) % len(self.spectra_color)
                        
                        self.spectra.append(self.active_spectrum)
                        self.update_menu_active_spectrum()                        
                        self.draw_active_spectrum()
                    elif elements == "continuum":
                        self.regions[elements] = read_continuum_regions(path)
                        self.draw_regions(elements)
                        self.not_saved[elements] = False
                        self.update_title()
                    elif elements == "lines":
                        self.regions[elements] = read_line_regions(path)
                        self.draw_regions(elements)
                        self.not_saved[elements] = False
                        self.update_title()
                    else:
                        # 'segments'
                        self.regions[elements] = read_segment_regions(path)
                        self.draw_regions(elements)
                        self.not_saved[elements] = False
                        self.update_title()
                    self.filenames[elements] = path
                    self.canvas.draw()
                    self.flash_status_message("Opened file %s" % path)
                    action_ended = True
                except Exception:
                    dlg_error = wx.MessageDialog(self, 'File %s does not have a compatible format.' % dlg.GetFilename(), 'File format incompatible', wx.OK | wx.ICON_ERROR)
                    dlg_error.ShowModal()
                    dlg_error.Destroy()
                    continue
            else:
                self.flash_status_message("Discarded.")
                action_ended = True
    
    def on_open_spectra(self, event):
        self.open_file("spectra")
    
    def on_open_continuum(self, event):
        self.open_file("continuum")
    
    def on_open_lines(self, event):
        self.open_file("lines")
    
    def on_open_segments(self, event):
        self.open_file("segments")
    
    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        elements = "spectra"
        if self.filenames[elements] != None:
            filename = self.filenames[elements].split('/')[-1]
            filename_length = len(filename)
            dirname = self.filenames[elements][:-filename_length]
            filename = filename + ".png"
        else:
            filename = elements + ".png"
            dirname = os.getcwd()
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Save plot as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    dlg_confirm = wx.MessageDialog(self, 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename(), 'File already exists', wx.YES|wx.NO|wx.ICON_QUESTION)
                    if dlg_confirm.ShowModal() != wx.ID_YES:
                        continue # Give the oportunity to select a new file name
                self.canvas.print_figure(path, dpi=self.dpi)
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                action_ended = True
    
    
    def on_save_spectra(self, event):
        file_choices = "All|*"
        
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            dlg_error = wx.MessageDialog(self, "There is no spectrum to be saved", 'Spectrum not present', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        if self.active_spectrum.path != None:
            filename = self.active_spectrum.path.split('/')[-1]
            filename_length = len(filename)
            dirname = self.active_spectrum.path[:-filename_length]
        else:
            filename = self.active_spectrum.name + ".txt"
            dirname = os.getcwd()
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Save spectrum as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    dlg_confirm = wx.MessageDialog(self, 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename(), 'File already exists', wx.YES|wx.NO|wx.ICON_QUESTION)
                    if dlg_confirm.ShowModal() != wx.ID_YES:
                        continue # Give the oportunity to select a new file name
                
                self.status_message("Saving %s..." % path)
                # Save, compress if the filename ends with ".gz"
                write_spectra(self.active_spectrum.data, path, compress=(path.split('.')[-1] == "gz"))
                self.active_spectrum.not_saved = False
                self.update_title()
                
                # Change name and path
                self.active_spectrum.path = path
                self.active_spectrum.name = path.split('/')[-1]
                self.active_spectrum.plot_id.set_label("[A] " + self.active_spectrum.name)
                self.update_legend()
                self.canvas.draw()
                
                # Menu active spectra
                for item in self.menu_active_spectrum.GetMenuItems():
                    if item.IsChecked():
                        item.SetItemLabel(self.active_spectrum.name)
                        break
                
                saved = True
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                saved = False
                action_ended = True
        
        return saved
        
    
    
    ## Save to a file
    # elements can be "continuum", "lines" or "segments"
    def save_regions(self, elements):
        file_choices = "All|*"
        saved = False
        elements = elements.lower()
        
        if len(self.region_widgets[elements]) == 0:
            dlg_error = wx.MessageDialog(self, "There is no regions to be saved", 'Empty regions', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        if self.filenames[elements] != None:
            filename = self.filenames[elements].split('/')[-1]
            filename_length = len(filename)
            dirname = self.filenames[elements][:-filename_length]
        else:
            filename = elements + ".txt"
            dirname = os.getcwd()
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Save regions as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    dlg_confirm = wx.MessageDialog(self, 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename(), 'File already exists', wx.YES|wx.NO|wx.ICON_QUESTION)
                    if dlg_confirm.ShowModal() != wx.ID_YES:
                        continue # Give the oportunity to select a new file name
                self.filenames[elements] = path
                output = open(path, "w")
                # Write header
                if elements == "lines":
                    output.write("wave_peak\twave_base\twave_top\tnote\n")
                else:
                    output.write("wave_base\twave_top\n")
                
                ## Visible regions: maybe they have been modified, removed or created
                for region in self.region_widgets[elements]:
                    # If the widget is not visible, it has been deleted by the user
                    if region.axvspan.get_visible():
                        if elements == "lines":
                            if region.note != None:
                                note = region.note.get_text()
                            else:
                                note = ""
                            output.write("%.4f" % region.get_wave_peak() + "\t" + "%.4f" % region.get_wave_base() + "\t" + "%.4f" % region.get_wave_top() + "\t" + note + "\n")
                        else:
                            output.write("%.4f" % region.get_wave_base() + "\t" + "%.4f" % region.get_wave_top() + "\n")
                
                output.close()
                self.regions_saved(elements)
                saved = True
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                saved = False
                action_ended = True
        
        return saved
    
    
    def on_save_continuum_regions(self, event):
        self.save_regions("continuum")
    
    def on_save_line_regions(self, event):
        self.save_regions("lines")
    
    def on_save_segments(self, event):
        self.save_regions("segments")
    
    def text2float(self, value, msg):
        result = None
        try:
            result = float(value)
        except ValueError:
            dlg_error = wx.MessageDialog(self, msg, 'Bad value format', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
        
        return result
    
    def on_fit_lines(self, event):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            dlg_error = wx.MessageDialog(self, "Spectra not loaded, there is nothing to fit.", 'Spectra not loaded', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        if self.active_spectrum.continuum_model == None:
            dlg_error = wx.MessageDialog(self, "Please, execute a general continuum fit first.", 'Continuum model not fitted', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        # Remove drawd lines if they exist
        self.remove_drawn_fitted_lines()
        
        for region in self.region_widgets["lines"]:
            wave_base = region.get_wave_base()
            wave_top = region.get_wave_top()
            loc = region.get_wave_peak()
            wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
            spectra_window = self.active_spectrum.data[wave_filter]
        
            try:
                # Remove old possible results
                region.line_model[self.active_spectrum] = None
                region.continuum_base_level[self.active_spectrum] = None
                
                # Fit
                line_model, continuum_value, lineflux, ew = fit_line(spectra_window, loc, continuum_model=self.active_spectrum.continuum_model)
                
                # Save results
                region.line_model[self.active_spectrum] = line_model
                region.continuum_base_level[self.active_spectrum] = np.mean(continuum_value)
            except Exception, e:
                print "Error:", wave_base, wave_top, e
        # Draw new lines
        self.draw_fitted_lines()
        self.flash_status_message("Lines fitted.")
    
    def remove_drawn_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_plot_id.has_key(self.active_spectrum) and region.line_plot_id[self.active_spectrum] != None:
                self.axes.lines.remove(region.line_plot_id[self.active_spectrum])
                region.line_plot_id[self.active_spectrum] = None
        self.canvas.draw()
        
    
    def draw_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum) and region.line_model[self.active_spectrum] != None:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()
                wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
                spectra_window = self.active_spectrum.data[wave_filter]
                
                spectra_line = get_spectra_from_model(region.line_model[self.active_spectrum], spectra_window['waveobs'])
                # Add the continuum base because the line_model has substracted it
                spectra_line['flux'] += region.continuum_base_level[self.active_spectrum]
                
                # zorder = 4, above the line region
                line_plot_id = self.axes.plot(spectra_line['waveobs'], spectra_line['flux'], lw=1, color='red', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=4)[0]
                region.line_plot_id[self.active_spectrum] = line_plot_id
        self.canvas.draw()
    
    def on_change_active_spectrum(self, event):
        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectra()
        
        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()
        
        # Determine the new spectra and draw its continuum if exists
        i = 0
        for item in self.menu_active_spectrum.GetMenuItems():
            if item.IsChecked():
                # Remove "[A]  " from spectra name (legend)
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
                # Change active spectrum
                self.active_spectrum = self.spectra[i]
                self.draw_active_spectrum()
                
                self.draw_continuum_spectra()
                self.draw_fitted_lines()
                self.canvas.draw()
                self.flash_status_message("Active spectrum: %s" % self.active_spectrum.name)
                break
            i += 1
        
    
    def on_fit_continuum(self, event):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            dlg_error = wx.MessageDialog(self, "Spectra not loaded, there is nothing to fit.", 'Spectra not loaded', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        dlg = FitContinuumDialog(self, -1, "Properties for fitting continuum")
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        
        nknots = self.text2float(dlg.nknots.GetValue(), 'Number of knots value is not a valid one.')
        in_continuum = dlg.radio_button_continuum.GetValue()
        dlg.Destroy()
        
        if nknots == None:
            self.flash_status_message("Bad value.")
            return
        
        self.update_progress(25)
        thread = threading.Thread(target=self.on_fit_continuum_thread, args=(nknots,), kwargs={'in_continuum':in_continuum})
        thread.setDaemon(True)
        thread.start()
        
    def on_fit_continuum_thread(self, nknots, in_continuum=False):
        if in_continuum:
            # Select from the spectra the regions that should be used to fit the continuum
            spectra_regions = None
            for region in self.region_widgets["continuum"]:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()
                wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
                new_spectra_region = self.active_spectrum.data[wave_filter]
                if spectra_regions == None:
                    spectra_regions = new_spectra_region
                else:
                    spectra_regions = np.hstack((spectra_regions, new_spectra_region))
        else:
            spectra_regions = self.active_spectrum.data
        
        if spectra_regions != None:
            self.active_spectrum.continuum_model = fit_continuum(spectra_regions, nknots=nknots)
            self.active_spectrum.continuum_data = get_spectra_from_model(self.active_spectrum.continuum_model, self.active_spectrum.data['waveobs'])
            wx.CallAfter(self.on_fit_continuum_finish, nknots)
        else:
            wx.CallAfter(self.flash_status_message, "No continuum regions found.")
    
    def on_fit_continuum_finish(self, nknots):
        self.draw_continuum_spectra()
        self.canvas.draw()
        self.flash_status_message("Continuum fitted with %s knots uniform Spline model." % str(nknots))
    
    def on_find_continuum(self, event):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            dlg_error = wx.MessageDialog(self, "Spectra not loaded, there is nothing to find.", 'Spectra not loaded', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        if self.active_spectrum.continuum_model == None:
            dlg_error = wx.MessageDialog(self, "Please, execute a general continuum fit first.", 'Continuum model not fitted', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        if self.not_saved["continuum"]:
            dlg = wx.MessageDialog(self, 'Are you sure you want to find new continuum regions without saving the current ones?', 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() != wx.ID_YES:
                self.flash_status_message("Discarded.")
                return
        
        dlg = FindContinuumDialog(self, -1, "Properties for finding continuum regions")
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        #~ resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
        resolution = None
        fixed_wave_step = self.text2float(dlg.fixed_wave_step.GetValue(), 'Maximum standard deviation value is not a valid one.')
        sigma = self.text2float(dlg.sigma.GetValue(), 'Maximum standard deviation value is not a valid one.')
        in_segments = dlg.radio_button_segments.GetValue()
        dlg.Destroy()
        
        if fixed_wave_step == None or sigma == None:
            wx.CallAfter(self.flash_status_message, "Bad value.")
            return
        
        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            wx.CallAfter(self.flash_status_message, "No segments found.")
            return
    
        thread = threading.Thread(target=self.on_find_continuum_thread, args=(resolution, sigma,), kwargs={'fixed_wave_step':fixed_wave_step, 'in_segments':in_segments, })
        thread.setDaemon(True)
        thread.start()
    
    def update_numpy_arrays_from_widgets(self, elements):
        total_regions = len(self.region_widgets[elements])
        if elements == "lines":
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                    self.regions[elements]['wave_peak'][i] = region.get_wave_peak()
                    self.regions[elements]['note'][i] = region.get_note_text()
                i += 1
        else:
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_base', float),('wave_top', float)])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                i += 1
    
    def on_find_continuum_thread(self, resolution, sigma, fixed_wave_step=None, in_segments=False):
        wx.CallAfter(self.status_message, "Finding continuum regions...")
        if in_segments:
            # TODO: Update self.regions["segments"] with self.region_widgets["segments"]
            self.update_numpy_arrays_from_widgets("segments")
            continuum_regions = find_continuum_on_regions(self.active_spectrum.data, resolution, self.regions["segments"], log_filename=None, max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, fixed_wave_step=fixed_wave_step, frame=self)
        else:
            continuum_regions = find_continuum(self.active_spectrum.data, resolution, log_filename=None, max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, fixed_wave_step=fixed_wave_step, frame=self)
        continuum_regions = merge_regions(self.active_spectrum.data, continuum_regions)
        
        wx.CallAfter(self.on_find_continuum_finish, continuum_regions)
    
    def on_find_continuum_finish(self, continuum_regions):
        elements = "continuum"
        self.regions[elements] = continuum_regions
        self.draw_regions(elements)
        self.not_saved[elements] = True
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Automatic finding of continuum regions ended.")
                    
    def on_determine_rv(self, event):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            dlg_error = wx.MessageDialog(self, "Spectra not loaded, there is nothing to determine.", 'Spectra not loaded', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        if self.active_spectrum.continuum_model == None:
            dlg_error = wx.MessageDialog(self, "Please, execute a general continuum fit first.", 'Continuum model not fitted', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        self.update_numpy_arrays_from_widgets("lines")
        
        if len(self.regions["lines"]) == 0:
            dlg_error = wx.MessageDialog(self, "Please, create or load a line list first.", 'Lines not present', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
            return
        
        dlg = DetermineRVDialog(self, -1, "Radial velocity determination", rv_limit=self.rv_limit, rv_step=self.rv_step)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        rv_limit = self.text2float(dlg.rv_limit.GetValue(), 'Radial velocity limit is not a valid one.')
        rv_step = self.text2float(dlg.rv_step.GetValue(), 'Radial velocity step is not a valid one.')
        renormalize = False #dlg.renormalize.GetValue()
        dlg.Destroy()
        
        if rv_limit == None or rv_step == None:
            self.flash_status_message("Bad value.")
            return
        self.rv_limit = rv_limit
        self.rv_step = rv_step
        
        thread = threading.Thread(target=self.on_determine_rv_thread, args=(renormalize,))
        thread.setDaemon(True)
        thread.start()

    def on_determine_rv_thread(self, renormalize):
        wx.CallAfter(self.status_message, "Determining RV...")
        xcoord, fluxes = build_radial_velocity_profile(self.active_spectrum.data, self.active_spectrum.continuum_model, self.regions["lines"], rv_limit=self.rv_limit, rv_step=self.rv_step, frame=self)
        wx.CallAfter(self.on_determine_rv_finish, xcoord, fluxes, renormalize)
    
    def on_determine_rv_finish(self, xcoord, fluxes, renormalize):
        # Modelize
        model = model_radial_velocity_profile(xcoord, fluxes, renormalize=renormalize)
        self.active_spectrum.rv = np.round(model.mu, 2) # km/s
        self.flash_status_message("Radial velocity determined: " + str(self.active_spectrum.rv) + " km/s")
        
        dlg = RVProfileDialog(self, -1, "Radial velocity profile", xcoord, fluxes, model)
        dlg.ShowModal()

    def on_correct_rv(self, event):
        if self.active_spectrum != None:
            rv = self.active_spectrum.rv
        else:
            rv = 0
        dlg = CorrectRVDialog(self, -1, "Radial velocity correction", rv)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        radial_vel = self.text2float(dlg.rv.GetValue(), 'Radial velocity value is not a valid one.')
        in_regions = dlg.radio_button_regions.GetValue()
        dlg.Destroy()
        
        if radial_vel == None:
            self.flash_status_message("Bad value.")
            return
        
        if in_regions:
            for elements in ["lines", "continuum", "segments"]:
                self.update_numpy_arrays_from_widgets(elements)
                if len(self.regions[elements]) > 0:
                    if elements == "lines":
                        self.regions[elements] = correct_radial_velocity_regions(self.regions[elements], radial_vel, with_peak=True)
                    else:
                        self.regions[elements] = correct_radial_velocity_regions(self.regions[elements], radial_vel)
                    self.draw_regions(elements)
                    self.not_saved[elements] = False
        else:
            if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
                dlg_error = wx.MessageDialog(self, "Spectra not loaded, there is nothing to determine.", 'Spectra not loaded', wx.OK | wx.ICON_ERROR)
                dlg_error.ShowModal()
                dlg_error.Destroy()
                return
            self.active_spectrum.data = correct_radial_velocity(self.active_spectrum.data, radial_vel)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
            if self.active_spectrum.continuum_model != None:
                self.active_spectrum.continuum_data = correct_radial_velocity(self.active_spectrum.continuum_data, radial_vel)
                self.draw_continuum_spectra()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a radial velocity correction of %s." % radial_vel)
    
    def on_synthesize(self, event):
        
        if sys.modules.has_key('synthesizer'):
            if self.active_spectrum != None:
                wave_base = np.round(np.min(self.active_spectrum.data['waveobs']), 2)
                wave_top = np.round(np.max(self.active_spectrum.data['waveobs']), 2)
            else:
                wave_base = 515.0 # Magnesium triplet region
                wave_top = 525.0
            teff = 5770.0
            logg = 4.40
            MH = 0.02
            microturbulence_vel = 2.0
            resolution = 47000
            dlg = SyntheticSpectrumDialog(self, -1, "Synthetic spectrum generator", wave_base, wave_top, resolution, teff, logg, MH, microturbulence_vel)
            dlg.ShowModal()
            
            if not dlg.action_accepted:
                dlg.Destroy()
                return
            
            teff = self.text2float(dlg.teff.GetValue(), 'Effective temperature value is not a valid one.')
            logg = self.text2float(dlg.logg.GetValue(), 'Gravity (log g) value is not a valid one.')
            MH = self.text2float(dlg.MH.GetValue(), 'Metallicity [M/H] value is not a valid one.')
            microturbulence_vel = self.text2float(dlg.microturbulence_vel.GetValue(), 'Microturbulence velocity value is not a valid one.')
            resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
            wave_base = self.text2float(dlg.wave_base.GetValue(), 'Wavelength min. value is not a valid one.')
            wave_top = self.text2float(dlg.wave_top.GetValue(), 'Wavelength max. value is not a valid one.')
            
            in_segments = dlg.radio_button_segments.GetValue() # else in spectrum
            dlg.Destroy()
            
            if teff == None or logg == None or MH == None or microturbulence_vel == None or resolution == None or wave_base == None or wave_top == None:
                self.flash_status_message("Bad value.")
                return
                
            
            
            ## TODO: Control wavelength margins, they should not be bigger/lower thant the linelist
            
            if self.modeled_layers_pack == None:
                print "Loading modeled atmospheres..."
                self.status_message("Loading modeled atmospheres...")
                self.modeled_layers_pack = load_modeled_layers_pack()
            
            if not valid_objective(self.modeled_layers_pack, teff, logg, MH):
                dlg_error = wx.MessageDialog(self, "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models.", 'Out of the atmospheric models', wx.OK | wx.ICON_ERROR)
                dlg_error.ShowModal()
                dlg_error.Destroy()
                self.flash_status_message("Bad values.")
                return
            
            # Prepare atmosphere model
            print "Interpolating atmosphere model..."
            self.status_message("Interpolating atmosphere model...")
            layers = interpolate_atmosphere_layers(self.modeled_layers_pack, teff, logg, MH)
            atm_filename = write_atmosphere(teff, logg, MH, layers)
            
            # Generate
            if not in_segments:
                if wave_base >= wave_top:
                    dlg_error = wx.MessageDialog(self, "Bad wavelength range definition, maximum value cannot be lower than minimum value.", 'Wavelength range', wx.OK | wx.ICON_ERROR)
                    dlg_error.ShowModal()
                    dlg_error.Destroy()
                    self.flash_status_message("Bad values.")
                    return
                
                waveobs = generate_wavelength_grid(wave_base, wave_top, resolution, points_per_fwhm = 3)
            else:
                #in_segments
                if len(self.region_widgets["segments"]) == 0:
                    self.flash_status_message("No segments present for synthetic spectrum generation.")
                    return
                
                # Build wavelength points from regions
                waveobs = None
                for region in self.region_widgets["segments"]:
                    wave_base = region.get_wave_base()
                    wave_top = region.get_wave_top()
                    
                    new_waveobs = generate_wavelength_grid(wave_base, wave_top, resolution, points_per_fwhm = 3)
                    if waveobs == None:
                        waveobs = new_waveobs
                    else:
                        waveobs = np.hstack((waveobs, new_waveobs))
            
            # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectra
            if np.min(waveobs) <= 300.0 or np.max(waveobs) >= 1100.0:
                # luke.300_1000nm.lst
                dlg_error = wx.MessageDialog(self, "Wavelength range is outside line list for spectrum generation.", 'Wavelength range', wx.OK | wx.ICON_ERROR)
                dlg_error.ShowModal()
                dlg_error.Destroy()
                self.flash_status_message("Bad values.")
                return
            
            total_points = len(waveobs)
            
            if total_points < 2:
                dlg_error = wx.MessageDialog(self, "Wavelength range too narrow.", 'Wavelength range', wx.OK | wx.ICON_ERROR)
                dlg_error.ShowModal()
                dlg_error.Destroy()
                self.flash_status_message("Bad values.")
                return
            
            print "Synthesizing spectrum..."
            self.status_message("Synthesizing spectrum...")
            
            synth_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
            synth_spectra['waveobs'] = waveobs
            
            synth_spectra['flux'] = synthesizer.spectrum(synth_spectra['waveobs']*10.0, atm_filename, linelist_file = "input/luke.300_1100nm.lst", abundances_file = "input/stdatom.dat", microturbulence_vel = microturbulence_vel, verbose=1)
                
            
            synth_spectra.sort(order='waveobs') # Make sure it is ordered by wavelength
            
            # Remove atmosphere model temporary file
            os.remove(atm_filename)
            
            # Remove current continuum from plot if exists
            self.remove_drawn_continuum_spectra()
            
            # Remove current drawn fitted lines if they exist
            self.remove_drawn_fitted_lines()
            
            # Remove "[A]  " from spectra name (legend) if it exists
            if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
            
            # Count how many synthetic spectra there are with the same name and find the max identificator number
            name = "Synth_" + str(teff) + "_" + str(logg) + "_"  + str(MH) + "_" + str(microturbulence_vel)
            num_synthetic = 0
            max_num = 0
            for spec in self.spectra:
                if spec.name.startswith(name):
                    try:
                        num = int(spec.name.split("-")[-1])
                        if num > max_num:
                            max_num = num
                    except ValueError:
                        pass
                    num_synthetic += 1
            
            if num_synthetic > 0: # There are repeated names
                name = name + "-" + str(max_num+1) # Add identificator number + 1
            
            self.active_spectrum = Spectrum(synth_spectra, name, color=self.spectra_color[self.spectra_color_id])
            self.spectra_color_id = (self.spectra_color_id + 1) % len(self.spectra_color)
            
            self.spectra.append(self.active_spectrum)
            self.active_spectrum.not_saved = True
            self.update_title()
            self.update_menu_active_spectrum()                        
            self.draw_active_spectrum()
            
            self.canvas.draw()

            self.flash_status_message("Synthetic spectra generated!")

    
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ Spectra Visual Editor is a tool for the definition of continuum regions, line masks and segments needed for the execution of the SME tool (Spectroscopy Made Easy).
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=3000, progress=True):
        self.timeroff = wx.Timer(self)
        
        self.statusbar.SetStatusText(msg)
        if progress:
            self.gauge.SetValue(pos=100)
        
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
        
    def status_message(self, msg):
        if self.timeroff != None:
            self.timeroff.Stop()
            self.on_flash_status_off(None)
        
        self.statusbar.SetStatusText(msg)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
        # Progress bar to zero
        self.gauge.SetValue(pos=0)
        self.timeroff = None


## Print usage
def usage():
    print "Usage:"
    print sys.argv[0], "[--continuum=file] [--lines=file] [--segments=file] [spectra_file]"

## Interpret arguments
def get_arguments():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "cls", ["continuum=", "lines=", "segments="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
    continuum_file = None
    lines_file = None
    segments_file = None
    spectra_file = None
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
              
    if (len(args) == 1):
        spectra_file = args[0]
        if not os.path.exists(spectra_file):
            print "Spectra file", args[0], "does not exists!"
            sys.exit(2)
    
    filenames = {}
    filenames['spectra'] = spectra_file
    filenames['continuum'] = continuum_file
    filenames['lines'] = lines_file
    filenames['segments'] = segments_file

    return filenames

#~ Example:
#~ python interactive.py --continuum=input/LUMBA/UVES_MRD_sun_cmask.txt --lines=input/LUMBA/UVES_MRD_sun_Fe-linelist.txt --segments=input/LUMBA/UVES_MRD_sun_segments.txt input/LUMBA/UVES_MRD_sun_official.s.gz

## Input files should be '\t' separated, comments with # and the first line is the header
if __name__ == '__main__':
    filenames = get_arguments()    
    
    #### Read files

    if filenames['spectra'] != None:
        try:
            name = filenames['spectra'].split('/')[-1]
            spectra = read_spectra(filenames['spectra'])
        except Exception:
            print "Spectra file", filenames['spectra'], "has an incompatible format!"
            sys.exit(2)
    else:
        spectra = None #np.zeros((0,), dtype=[('waveobs', '<f8'), ('flux', '<f8'), ('err', '<f8')])
        name = None
    
    if filenames['continuum'] != None:
        try:
            continuum = read_continuum_regions(filenames['continuum'])
        except Exception:
            print "Continuum file", filenames['continuum'], "has an incompatible format!"
            sys.exit(2)
        
        ## Validations
        if np.any((continuum['wave_top'] - continuum['wave_base']) < 0):
            raise Exception("Continuum: wave_top cannot be smaller than wave_base")
    else:
        continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
    
    if filenames['lines'] != None:
        try:
            lines = read_line_regions(filenames['lines'])
        except Exception:
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
            segments = read_segment_regions(filenames['segments'])
        except Exception:
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
    ## wxPython makes that functions like atof(), called from external
    ## C libraries through Cython, behave interpreting comma as decimal separator
    ## and brokes the code that expects to have dots as decimal separators
    ## (i.e. SPECTRUM and its external files like linelists)
    locales = {
        u'en' : (wx.LANGUAGE_ENGLISH, u'en_US.UTF-8'),
        u'es' : (wx.LANGUAGE_SPANISH, u'es_ES.UTF-8'),
        u'fr' : (wx.LANGUAGE_FRENCH, u'fr_FR.UTF-8'),
        }
    os.environ['LANG'] = locales[u'en'][1]

    app = wx.PySimpleApp()
    app.frame = SpectraFrame(spectra, name, regions, filenames)
    app.frame.Show()
    app.MainLoop()


#### IGNORE: Code for original files
#~ continuum = asciitable.read(table=filenames['continuum'], delimiter=' ', comment=';', names=['wave_base', 'wave_top'], data_start=0)
#~ continuum['wave_base'] = continuum['wave_base'] / 10 # nm
#~ continuum['wave_top'] = continuum['wave_top'] / 10 # nm

#~ lines = asciitable.read(table=filenames['lines'], delimiter=' ', comment=';', names=['wave_peak', 'wave_base', 'wave_top'], data_start=0)
#~ lines['wave_base'] = lines['wave_base'] / 10 # nm
#~ lines['wave_peak'] = lines['wave_peak'] / 10 # nm
#~ lines['wave_top'] = lines['wave_top'] / 10 # nm

#~ segments = asciitable.read(table=filenames['segments'], delimiter=' ', comment=';', names=['wave_base', 'wave_top'], data_start=0)
#~ segments['wave_base'] = segments['wave_base'] / 10 # nm
#~ segments['wave_top'] = segments['wave_top'] / 10 # nm

