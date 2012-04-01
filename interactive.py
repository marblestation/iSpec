#!/usr/bin/env python
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

if os.path.exists("synthesizer.so") and os.path.exists("atmospheres.py"):
    try:
        import synthesizer
        from atmospheres import *
    except ImportError as e:
        print "Spectra synthesizer not available"
else:
    #print "Spectra synthesizer not available"
    pass


# The recommended way to use wx with mpl is with the WXAgg backend.
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from common import *
from continuum import *
from lines import *
from radial_velocity import *
from convolve import *


class FitContinuumDialog(wx.Dialog):
    def __init__(self, parent, id, title, nknots=1):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_nknots = wx.StaticText(self, -1, "Number of knots for a uniform spline fit: ", style=wx.ALIGN_LEFT)
        self.nknots = wx.TextCtrl(self, -1, str(nknots),  style=wx.TE_RIGHT)
        
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
    def __init__(self, parent, id, title, rv_upper_limit, rv_lower_limit, rv_step):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### RV lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_rv_lower_limit = wx.StaticText(self, -1, "Radial velocity lower limit (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_lower_limit = wx.TextCtrl(self, -1, str(int(rv_lower_limit)),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_lower_limit, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_lower_limit, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### RV upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_rv_upper_limit = wx.StaticText(self, -1, "Radial velocity upper limit (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_upper_limit = wx.TextCtrl(self, -1, str(int(rv_upper_limit)),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_upper_limit, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_upper_limit, 0, border=3, flag=flags)
        
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

class DetermineBarycentricCorrectionDialog(wx.Dialog):
    def __init__(self, parent, id, title, day, month, year, hours, minutes, seconds, ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        # Date
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_date = wx.StaticText(self, -1, "Date: ", style=wx.ALIGN_LEFT)
        self.day = wx.TextCtrl(self, -1, str(int(day)),  style=wx.TE_RIGHT)
        self.month = wx.TextCtrl(self, -1, str(int(month)),  style=wx.TE_RIGHT)
        self.year = wx.TextCtrl(self, -1, str(int(year)),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_date, 0, border=3, flag=flags)
        self.hbox.Add(self.day, 0, border=3, flag=flags)
        self.hbox.Add(self.month, 0, border=3, flag=flags)
        self.hbox.Add(self.year, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        # Time
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_time = wx.StaticText(self, -1, "Time: ", style=wx.ALIGN_LEFT)
        self.hours = wx.TextCtrl(self, -1, str(int(hours)),  style=wx.TE_RIGHT)
        self.minutes = wx.TextCtrl(self, -1, str(int(minutes)),  style=wx.TE_RIGHT)
        self.seconds = wx.TextCtrl(self, -1, str(seconds),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_time, 0, border=3, flag=flags)
        self.hbox.Add(self.hours, 0, border=3, flag=flags)
        self.hbox.Add(self.minutes, 0, border=3, flag=flags)
        self.hbox.Add(self.seconds, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        # Right Ascension
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_ra = wx.StaticText(self, -1, "Right Ascension: ", style=wx.ALIGN_LEFT)
        self.ra_hours = wx.TextCtrl(self, -1, str(int(ra_hours)),  style=wx.TE_RIGHT)
        self.ra_minutes = wx.TextCtrl(self, -1, str(int(ra_minutes)),  style=wx.TE_RIGHT)
        self.ra_seconds = wx.TextCtrl(self, -1, str(ra_seconds),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_ra, 0, border=3, flag=flags)
        self.hbox.Add(self.ra_hours, 0, border=3, flag=flags)
        self.hbox.Add(self.ra_minutes, 0, border=3, flag=flags)
        self.hbox.Add(self.ra_seconds, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        # Declination
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_dec = wx.StaticText(self, -1, "Declination: ", style=wx.ALIGN_LEFT)
        self.dec_degrees = wx.TextCtrl(self, -1, str(int(dec_degrees)),  style=wx.TE_RIGHT)
        self.dec_minutes = wx.TextCtrl(self, -1, str(int(dec_minutes)),  style=wx.TE_RIGHT)
        self.dec_seconds = wx.TextCtrl(self, -1, str(dec_seconds),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_dec, 0, border=3, flag=flags)
        self.hbox.Add(self.dec_degrees, 0, border=3, flag=flags)
        self.hbox.Add(self.dec_minutes, 0, border=3, flag=flags)
        self.hbox.Add(self.dec_seconds, 0, border=3, flag=flags)
        
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
    def __init__(self, parent, id, title, xcoord, fluxes, model, num_used_lines, rv_step):
        wx.Dialog.__init__(self, parent, id, title, size=(600, 600))
        
        self.recalculate = False
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
        
        self.text_question = wx.StaticText(self, -1, "Recalculate again? ", style=wx.ALIGN_CENTER)
        
        sizer =  self.CreateButtonSizer(wx.YES_NO | wx.NO_DEFAULT)
        self.vbox.Add(self.text_question, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_yes, id=wx.ID_YES)
        self.Bind(wx.EVT_BUTTON, self.on_no, id=wx.ID_NO)
        
        ## Draw
        self.axes.plot(xcoord, fluxes+1, lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)
        if rv_step >= 0.1:
            xcoord_mod = np.arange(np.min(xcoord), np.max(xcoord), 0.1)
            self.axes.plot(xcoord_mod, model(xcoord_mod)+1, lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)
        else:
            self.axes.plot(xcoord, model(xcoord)+1, lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)
        
        self.axes.grid(True, which="both")
        self.axes.set_title("Profile", fontsize="10")
        self.axes.set_xlabel("radial velocity (km/s)", fontsize="10")
        self.axes.set_ylabel("relative intensity", fontsize="10")
        
        ## Stats
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, "Mean (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(model.mu, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Min. error (+/- km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(rv_step, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Sigma (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(model.sig, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "A (rel. intensity)")
        self.stats.SetStringItem(num_items, 1, str(np.round(model.A+1, 2)))
        fwhm = model.sig * (2*np.sqrt(2*np.log(2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "FWHM (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(fwhm, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Estimated resolution/resolving power (R)")
        c = 299792458.0 # m/s
        self.stats.SetStringItem(num_items, 1, str(np.int(c/(1000.0*np.round(fwhm, 2)))))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Number of lines used")
        self.stats.SetStringItem(num_items, 1, str(num_used_lines))

    def on_no(self, event):
        self.recalculate = False
        self.Close()
    
    def on_yes(self, event):
        self.recalculate = True
        self.Close()

class CutSpectraDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Wave lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_base = wx.StaticText(self, -1, "Base wavelength: ", style=wx.ALIGN_LEFT)
        self.wave_base = wx.TextCtrl(self, -1, str(wave_base),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_base, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_base, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Wave upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_top = wx.StaticText(self, -1, "Top wavelength: ", style=wx.ALIGN_LEFT)
        self.wave_top = wx.TextCtrl(self, -1, str(wave_top),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_top, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_top, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        self.vbox.AddSpacer(10)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        

    def on_ok(self, event):
        self.action_accepted = True
        self.Close()

class HomogenizeCombineSpectraDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top, resolution):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        ### Wave lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_base = wx.StaticText(self, -1, "Base wavelength: ", style=wx.ALIGN_LEFT)
        self.wave_base = wx.TextCtrl(self, -1, str(wave_base),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_base, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_base, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Wave upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_wave_top = wx.StaticText(self, -1, "Top wavelength: ", style=wx.ALIGN_LEFT)
        self.wave_top = wx.TextCtrl(self, -1, str(wave_top),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_top, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_top, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### Resolution
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_resolution = wx.StaticText(self, -1, "Resolution of the current spectra: ", style=wx.ALIGN_LEFT)
        self.resolution = wx.TextCtrl(self, -1, str(resolution),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_resolution, 0, border=3, flag=flags)
        self.hbox.Add(self.resolution, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        self.vbox.AddSpacer(10)
        
        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)
        

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

class DegradeResolutionDialog(wx.Dialog):
    def __init__(self, parent, id, title, from_resolution, to_resolution):
        wx.Dialog.__init__(self, parent, id, title)
        
        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        
        
        ### From resolution
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_from_resolution = wx.StaticText(self, -1, "Initial resolution: ", style=wx.ALIGN_LEFT)
        self.from_resolution = wx.TextCtrl(self, -1, str(from_resolution),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_from_resolution, 0, border=3, flag=flags)
        self.hbox.Add(self.from_resolution, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        ### To resolution
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        self.text_to_resolution = wx.StaticText(self, -1, "Final resolution: ", style=wx.ALIGN_LEFT)
        self.to_resolution = wx.TextCtrl(self, -1, str(to_resolution),  style=wx.TE_RIGHT)
        
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_to_resolution, 0, border=3, flag=flags)
        self.hbox.Add(self.to_resolution, 0, border=3, flag=flags)
        
        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)
        
        self.vbox.AddSpacer(10)
        
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
        
        note_text = self.frame.ask_value('Note for the new line region:', 'Note', note_text)
        
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
        
        if self.frame.operation_in_progress:
            return
        
        # When regions overlap two or more can receive the click event, so
        # let's use a lock to allow the modification of one of them
        if self.frame.lock.acquire(False):
            if self.frame.action == "Remove":
                # The real removal is in the "on_release" event
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    self.press = event.button, event.x, event.xdata
                else:
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
                
            
    def disconnect_and_remove(self):
        #self.hide()
        #self.frame.canvas.draw()
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
        self.axvspan.remove()
        if self.mark != None:
            self.mark.remove()
        if self.note != None:
            self.note.remove()
    
    def hide(self):
        self.axvspan.set_visible(False)
        if self.mark != None:
            self.mark.set_visible(False)
        if self.note != None:
            self.note.set_visible(False)
        if self.line_plot_id != None:
            for plot_id in self.line_plot_id.values():
                if plot_id != None:
                    self.frame.axes.lines.remove(plot_id)
        
    
    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        
        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return
        
        if self.press != None:
            # In modification mode, if it is the current selected widget
            self.press = None
            if self.frame.action == "Remove":
                self.frame.lock.release() # Release now because the next commands will 
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Remove note
                    if self.note != None:
                        self.note.set_visible(False)
                        self.note.remove()
                        self.note = None
                        self.frame.canvas.draw()
                else:
                    # Remove from the figure now (except for line marks)
                    #  (if we do it on_press, this event would not be trigered and the lock not released)
                    self.frame.region_widgets[self.element_type].remove(self)
                    self.frame.regions_changed(self.element_type)
                    self.frame.flash_status_message("Removed region from " + "%.4f" % self.axvspan.get_xy()[0,0] + " to " + "%.4f" % self.axvspan.get_xy()[2,0])
                    self.disconnect_and_remove()
                    self.hide()
                    self.frame.canvas.draw()
            else:
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
        if self.frame.action == "Modify":
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
        self.rv_profile_xcoord = None
        self.rv_profile_fluxes = None
        self.rv_profile_model = None
        self.rv_profile_num_used_lines = None
        self.rv_profile_rv_step = None
        self.resolution = 47000 # default

    
class SpectraFrame(wx.Frame):
    title = 'Spectra Visual Editor '
    
    def ipython(self):
        import IPython
        self.embedshell = IPython.Shell.IPShellEmbed(argv=['--pdb'])
        self.embedshell()
    
    # regions should be a dictionary with 'continuum', 'lines' and 'segments' keys
    # filenames should be a dictionary with 'spectra', 'continuum', 'lines' and 'segments' keys
    def __init__(self, spectra=None, regions=None, filenames=None):
        self.embedshell = None
        self.ipython_thread = None
        #self.ipython_thread = threading.Thread(target=self.ipython)
        #self.ipython_thread.setDaemon(True)
        #self.ipython_thread.start()
        
        self.spectra_colors = ('#0000FF', '#A52A2A', '#A020F0', '#34764A', '#000000', '#90EE90', '#FFA500', '#1E90FF',   '#FFC0CB', '#7F7F7F', '#00FF00',)
        
        self.spectra = []
        self.active_spectrum = None
        for i in np.arange(len(spectra)):
            path = filenames["spectra"][i]
            name = path.split('/')[-1]
            name = self.get_name(name) # If it already exists, add a suffix
            color = self.get_color()
            self.active_spectrum = Spectrum(spectra[i], name, path=path, color=color)
            self.spectra.append(self.active_spectrum)
        
        
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
            self.filenames = {}
            self.filenames['continuum'] = filenames['continuum']
            self.filenames['lines'] = filenames['lines']
            self.filenames['segments'] = filenames['segments']
        
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
        self.rv_lower_limit = -200 # km/s
        self.rv_upper_limit = 200 # km/s
        self.rv_step = 5.0 # km/s
        self.linelist_rv = None
        self.modeled_layers_pack = None # Synthesize spectrum (atmospheric models)
        
        # Barycentric velocity determination (default params)
        self.day = 15
        self.month = 2
        self.year = 2012
        self.hours = 0
        self.minutes = 0
        self.seconds = 0
        self.ra_hours = 19
        self.ra_minutes = 50
        self.ra_seconds = 46.99
        self.dec_degrees = 8
        self.dec_minutes = 52
        self.dec_seconds = 5.96
        self.barycentric_vel = 0.0 # km/s
        
        self.operation_in_progress = False
        
        wx.Frame.__init__(self, None, -1, self.title)
        self.icon = wx.Icon("images/SVE.ico", wx.BITMAP_TYPE_ICO)
        self.tbicon = wx.TaskBarIcon()
        self.tbicon.SetIcon(self.icon, "")
        wx.Frame.SetIcon(self, self.icon)
                
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.timeroff = wx.Timer(self)
        self.draw_figure_for_first_time()
    
    ####################################################################
    def error(self, title, msg):
        dlg_error = wx.MessageDialog(self, msg, title, wx.OK|wx.ICON_ERROR)
        dlg_error.ShowModal()
        dlg_error.Destroy()
    
    def question(self, title, msg):
        dlg_question = wx.MessageDialog(self, msg, title, wx.YES|wx.NO|wx.ICON_QUESTION)
        answer_yes = (dlg_question.ShowModal() == wx.ID_YES)
        dlg_question.Destroy()
        if not answer_yes:
            self.flash_status_message("Discarded")
        return answer_yes
    
    def ask_value(self, text, title, default_value):
        dlg_question = wx.TextEntryDialog(self, text, title)
        
        dlg_question.SetValue(default_value)
        if dlg_question.ShowModal() == wx.ID_OK:
            response = dlg_question.GetValue()
        else:
            response = None
            
        dlg_question.Destroy()
        return response
    ####################################################################
    
    def get_color(self):
        # Look for a free color
        free_color_found = None
        if len(self.spectra) < len(self.spectra_colors):
            for color in self.spectra_colors:
                discard_color = False
                for spec in self.spectra:
                    if spec.color == color:
                        discard_color = True
                        break
                if discard_color:
                    continue
                else:
                    free_color_found = color
                    break
        
        if free_color_found == None:
            good_color = False
            while not good_color:
                # Random number converted to hexadecimal
                random_num = np.random.randint(0, 16777215)
                random_color = "#%x" % random_num
                # Test that the generated color is correct
                try:
                    rgba = matplotlib.colors.colorConverter.to_rgba(random_color)
                    good_color = True
                except ValueError as e:
                    pass
            
            free_color_found = random_color
            
        return free_color_found
    ####################################################################
    
    # Check if exists a spectrum with that name, in that case add a suffix
    def get_name(self, name):
        num_repeated = 0
        max_num = 0
        for spec in self.spectra:
            if spec.name.startswith(name):
                try:
                    # Does it has already a suffix?
                    num = int(spec.name.split("-")[-1])
                    if num > max_num:
                        max_num = num
                except ValueError as e:
                    pass
                num_repeated += 1
            
        if num_repeated > 0: # There are repeated names
            name = name + "-" + str(max_num+1) # Add identificator number + 1
        return name
    
    
    def check_operation_in_progress(self):
        if self.operation_in_progress:
            msg = "This action cannot be done while there is another operation in progress"
            title = "Operation in progress"
            self.error(title, msg)
            return True
        return False
    
    def check_active_spectrum_exists(self):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            msg = "There is no active spectrum"
            title = "Spectrum not present"
            self.error(title, msg)
            return False
        return True
    
    def check_continuum_model_exists(self):
        if self.active_spectrum.continuum_model == None:
            msg = "Please, execute a general continuum fit first"
            title = "Continuum model not fitted"
            self.error(title, msg)
            return False
        return True
    ####################################################################
    
    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        self.spectrum_function_items = []
        
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
        self.spectrum_function_items.append(m_expt)
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
        m_expt = menu_edit.AppendMenu(-1, 'Select spectra', self.menu_active_spectrum)
        self.spectrum_function_items.append(m_expt)
        
        m_expt = m_close_spectrum = menu_edit.Append(-1, "Close spectrum", "Close current active spectrum")
        self.Bind(wx.EVT_MENU, self.on_close_spectrum, m_close_spectrum)
        self.spectrum_function_items.append(m_expt)
        menu_edit.AppendSeparator()
        
        m_fit_continuum = menu_edit.Append(-1, "Fit &continuum", "Fit continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_continuum, m_fit_continuum)
        self.spectrum_function_items.append(m_fit_continuum)
        m_fit_lines = menu_edit.Append(-1, "Fit &lines", "Fit lines using the fitted continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_lines, m_fit_lines)
        self.spectrum_function_items.append(m_fit_lines)
        menu_edit.AppendSeparator()
        
        m_remove_fitted_continuum = menu_edit.Append(-1, "Clear fitted continuum", "Remove the fitted continuum")
        self.Bind(wx.EVT_MENU, self.on_remove_fitted_continuum, m_remove_fitted_continuum)
        self.spectrum_function_items.append(m_remove_fitted_continuum)
        m_remove_fitted_lines = menu_edit.Append(-1, "Clear fitted lines", "Remove fitted lines")
        self.Bind(wx.EVT_MENU, self.on_remove_fitted_lines, m_remove_fitted_lines)
        self.spectrum_function_items.append(m_remove_fitted_lines)
        
        m_remove_continuum_regions = menu_edit.Append(-1, "Clear continuum regions", "Clear continuum regions")
        self.Bind(wx.EVT_MENU, self.on_remove_continuum_regions, m_remove_continuum_regions)
        m_remove_line_masks = menu_edit.Append(-1, "Clear line masks", "Clear line masks")
        self.Bind(wx.EVT_MENU, self.on_remove_line_masks, m_remove_line_masks)
        m_remove_segments = menu_edit.Append(-1, "Clear segments", "Clear segments")
        self.Bind(wx.EVT_MENU, self.on_remove_segments, m_remove_segments)
        menu_edit.AppendSeparator()
        
        m_find_continuum = menu_edit.Append(-1, "&Find continuum regions", "Find continuum regions")
        self.Bind(wx.EVT_MENU, self.on_find_continuum, m_find_continuum)
        self.spectrum_function_items.append(m_find_continuum)
        menu_edit.AppendSeparator()
        
        m_determine_barycentric_vel = menu_edit.Append(-1, "&Determine barycentric velocity", "Determine baricentryc velocity")
        self.Bind(wx.EVT_MENU, self.on_determine_barycentric_vel, m_determine_barycentric_vel)
        m_determine_rv = menu_edit.Append(-1, "&Determine radial velocity", "Determine radial velocity")
        self.Bind(wx.EVT_MENU, self.on_determine_rv, m_determine_rv)
        self.spectrum_function_items.append(m_determine_rv)
        m_correct_barycentric_vel = menu_edit.Append(-1, "Correct barycentric velocity", "Correct spectra by using its radial velocity")
        self.Bind(wx.EVT_MENU, self.on_correct_barycentric_vel, m_correct_barycentric_vel)
        self.spectrum_function_items.append(m_correct_barycentric_vel)
        m_correct_rv = menu_edit.Append(-1, "&Correct radial velocity", "Correct spectra by using its radial velocity")
        self.Bind(wx.EVT_MENU, self.on_correct_rv, m_correct_rv)
        self.spectrum_function_items.append(m_correct_rv)
        menu_edit.AppendSeparator()
        
        m_degrade_resolution = menu_edit.Append(-1, "Degrade resolution", "Degread spectrum resolution")
        self.Bind(wx.EVT_MENU, self.on_degrade_resolution, m_degrade_resolution)
        self.spectrum_function_items.append(m_degrade_resolution)
        m_normalize_spectrum = menu_edit.Append(-1, "Continuum normalization", "Normalize spectrum")
        self.Bind(wx.EVT_MENU, self.on_continuum_normalization, m_normalize_spectrum)
        self.spectrum_function_items.append(m_normalize_spectrum)
        m_cut_spectra = menu_edit.Append(-1, "Wavelength range reduction", "Reduce the wavelength range")
        self.Bind(wx.EVT_MENU, self.on_cut_spectra, m_cut_spectra)
        self.spectrum_function_items.append(m_cut_spectra)
        m_convert_to_nm = menu_edit.Append(-1, "Convert to nanometers", "Divide wavelength by 10")
        self.Bind(wx.EVT_MENU, self.on_convert_to_nm, m_convert_to_nm)
        self.spectrum_function_items.append(m_convert_to_nm)
        m_combine_spectra = menu_edit.Append(-1, "Combine all spectra", "Combine all spectra into one")
        self.Bind(wx.EVT_MENU, self.on_combine_spectra, m_combine_spectra)
        self.spectrum_function_items.append(m_combine_spectra)
        
        if sys.modules.has_key('synthesizer'):
            menu_edit.AppendSeparator()
            m_synthesize = menu_edit.Append(-1, "&Synthesize spectrum", "Synthesize spectrum")
            self.Bind(wx.EVT_MENU, self.on_synthesize, m_synthesize)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the visual editor")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)
        self.update_menu_active_spectrum()

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
        box = self.axes.get_position()
        self.axes.set_position([box.x0, box.y0, box.width * 0.80, box.height])
        
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
        self.gauge = wx.Gauge(self.panel, range=100, size=(100, 20), style=wx.GA_HORIZONTAL)
        
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
        if self.active_spectrum != None:
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
    
    def remove_continuum_spectra(self):
        self.remove_drawn_continuum_spectra()
        self.active_spectrum.continuum_model = None
        self.active_spectrum.continuum_data = None
    
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
            r.disconnect_and_remove()
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
    
    def draw_figure_for_first_time(self):
        """ Redraws the figure
        """
        self.axes.grid(True, which="both")
        self.axes.set_title("Spectra", fontsize="10")
        self.axes.set_xlabel("wavelength", fontsize="10")
        self.axes.set_ylabel("flux", fontsize="10")
        
        for spec in self.spectra:
            # Remove "[A]  " from spectra name (legend) if it exists
            if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
            self.active_spectrum = spec
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
        if self.operation_in_progress:
            msg = "There is an operation in progress, are you sure you want to exit anyway?"
            title = "Operation in progress"
            if not self.question(title, msg):
                return
        
        # Check if there is any spectrum not saved
        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"] or spectra_not_saved:
            msg = "Are you sure you want to exit without saving the regions/spectra?"
            title = "Changes not saved"
            if self.question(title, msg):
                self.tbicon.Destroy()
                self.Destroy()
        else:
            self.tbicon.Destroy()
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
        if self.timeroff.IsRunning(): return
        if self.operation_in_progress: return
        self.status_message("Cursor on wavelength %.4f" % event.xdata + " and flux %.4f" % event.ydata)

    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        if self.check_operation_in_progress():
            return
        
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
            self.menu_active_spectrum.Remove(i+1)
        
        # Add as many options as spectra
        for i in np.arange(len(self.spectra)):
            if self.spectra[i] == None:
                continue
            
            spec_element = self.menu_active_spectrum.Append(i+1, self.spectra[i].name, kind=wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.on_change_active_spectrum, spec_element)

            if self.active_spectrum == self.spectra[i]:
                spec_element.Check()
        
        if len(self.spectra) > 0:
            for item in self.spectrum_function_items:
                item.Enable(True)
        else:
            for item in self.spectrum_function_items:
                item.Enable(False)
    
    def update_legend(self):
        if len(self.spectra) == 0:
            self.axes.legend_ = None
        else:
            leg = self.axes.legend(loc='upper right', bbox_to_anchor=(1.4, 1.1), ncol=1, shadow=False)
            ltext  = leg.get_texts()
            plt.setp(ltext, fontsize='8')
    
    def add_stats(self, k, v):
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, k)
        self.stats.SetStringItem(num_items, 1, str(v))
        
    def update_stats(self, region):
        if not self.check_active_spectrum_exists():
            return
        
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
        
        if region.element_type == "lines" and region.line_plot_id.has_key(self.active_spectrum) and region.line_model[self.active_spectrum] != None:
            self.add_stats("Gaussian mean (mu)", "%.4f" % region.line_model[self.active_spectrum].mu)
            self.add_stats("Gaussian width (A)", "%.4f" % region.line_model[self.active_spectrum].A)
            self.add_stats("Gaussian standard deviation (sigma)", "%.4f" % region.line_model[self.active_spectrum].sig)
            self.add_stats("Gaussian base level (mean continuum)", "%.4f" % region.continuum_base_level[self.active_spectrum])
            rms = np.mean(np.abs(region.line_model[self.active_spectrum].residuals())) + np.std(np.abs(region.line_model[self.active_spectrum].residuals()))
            self.add_stats("Gaussian fit root mean squeare (RMS)", "%.4f" % rms)
        
        if self.active_spectrum.continuum_model != None:
            if num_points > 0:
                mean_continuum = np.mean(self.active_spectrum.continuum_model(spectra_window['waveobs']))
                self.add_stats("Continuum mean for the region", "%.4f" % mean_continuum)
            rms = np.mean(np.abs(self.active_spectrum.continuum_model.residuals())) + np.std(np.abs(self.active_spectrum.continuum_model.residuals()))
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
                        
                        name = self.get_name(path.split('/')[-1]) # If it already exists, add a suffix
                        color = self.get_color()
                        self.active_spectrum = Spectrum(new_spectra_data, name, path = path, color=color)
                        
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
                except Exception as e:
                    dlg_error = wx.MessageDialog(self, 'File %s does not have a compatible format.' % dlg.GetFilename(), 'File format incompatible', wx.OK | wx.ICON_ERROR)
                    dlg_error.ShowModal()
                    dlg_error.Destroy()
                    continue
            else:
                self.flash_status_message("Discarded.")
                action_ended = True
    
    def on_open_spectra(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("spectra")
    
    def on_open_continuum(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("continuum")
    
    def on_open_lines(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("lines")
    
    def on_open_segments(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("segments")
    
    def on_close_spectrum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        if self.active_spectrum.not_saved:
            msg = "Are you sure you want to close the spectrum without saving it?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        
        self.spectra.remove(self.active_spectrum)
        self.axes.lines.remove(self.active_spectrum.plot_id)
        
        # Remove fitted continuum if it exists
        if self.active_spectrum != None and self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
            self.active_spectrum.continuum_plot_id = None
            self.active_spectrum.continuum_model = None
            self.active_spectrum.continuum_data = None
        # Remove fitted lines if they exist
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum):
                if region.line_plot_id[self.active_spectrum] != None:
                    self.axes.lines.remove(region.line_plot_id[self.active_spectrum])
                del region.line_plot_id[self.active_spectrum]
                del region.line_model[self.active_spectrum]
                del region.continuum_base_level[self.active_spectrum]
        if len(self.spectra) == 0:
            self.active_spectrum = None
        else:
            self.active_spectrum = self.spectra[0]
        
        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectra()
        self.draw_fitted_lines()
        self.canvas.draw()
        self.flash_status_message("Spectrum closed.")
    
    def on_save_plot(self, event):
        if self.check_operation_in_progress():
            return
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
        if self.check_operation_in_progress():
            return
        file_choices = "All|*"
        
        if not self.check_active_spectrum_exists():
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
        if self.check_operation_in_progress():
            return
        self.save_regions("continuum")
    
    def on_save_line_regions(self, event):
        if self.check_operation_in_progress():
            return
        self.save_regions("lines")
    
    def on_save_segments(self, event):
        if self.check_operation_in_progress():
            return
        self.save_regions("segments")
    
    def text2float(self, value, msg):
        result = None
        try:
            result = float(value)
        except ValueError as e:
            dlg_error = wx.MessageDialog(self, msg, 'Bad value format', wx.OK | wx.ICON_ERROR)
            dlg_error.ShowModal()
            dlg_error.Destroy()
        
        return result
    
    def on_remove_fitted_lines(self, event):
        if self.check_operation_in_progress():
            return
        self.remove_fitted_lines()
        
    def remove_regions(self, elements):
        if self.not_saved[elements]:
            msg = "Are you sure you want to remove this regions without saving them?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        if elements == "continuum":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        elif elements == "lines":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
        else:
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        self.draw_regions(elements)
        self.not_saved[elements] = False
        self.update_title()
        self.filenames[elements] = None
        self.canvas.draw()
    
    def on_remove_continuum_regions(self, event):
        if self.check_operation_in_progress():
            return
        elements = "continuum"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned continuum regions")
        
    
    def on_remove_line_masks(self, event):
        if self.check_operation_in_progress():
            return
        elements = "lines"
        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements)
        self.flash_status_message("Cleaned line masks")
    
    def on_remove_segments(self, event):
        if self.check_operation_in_progress():
            return
        elements = "segments"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned segments")
    
    def on_degrade_resolution(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = DegradeResolutionDialog(self, -1, "Degrade spectrum resolution", self.active_spectrum.resolution, 10000)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        from_resolution = self.text2float(dlg.from_resolution.GetValue(), 'Initial resolution value is not a valid one.')
        to_resolution = self.text2float(dlg.to_resolution.GetValue(), 'Final resolution value is not a valid one.')
        dlg.Destroy()
        
        if from_resolution == None or to_resolution == None or from_resolution <= to_resolution:
            self.flash_status_message("Bad value.")
            return
        
        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to degrade its resolution now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_degrade_resolution_thread, args=(from_resolution, to_resolution,))
        thread.setDaemon(True)
        thread.start()
    
    def on_degrade_resolution_thread(self, from_resolution, to_resolution):
        wx.CallAfter(self.status_message, "Degrading spectrum resolution...")
        convolved_spectra = convolve_spectra(self.active_spectrum.data, from_resolution, to_resolution, frame=self)
        wx.CallAfter(self.on_degrade_resolution_finnish, convolved_spectra, from_resolution, to_resolution)
    
    def on_degrade_resolution_finnish(self, convolved_spectra, from_resolution, to_resolution):
        self.active_spectrum.data = convolved_spectra
        self.active_spectrum.not_saved = True
        self.active_spectrum.resolution = to_resolution
        
        # Remove current continuum from plot if exists
        self.remove_continuum_spectra()
        
        # Remove current drawn fitted lines if they exist
        self.remove_fitted_lines()
        
        self.draw_active_spectrum()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Spectrum degreaded from resolution %.2f to %.2f" % (from_resolution, to_resolution))
        self.operation_in_progress = False
    
    def on_cut_spectra(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = CutSpectraDialog(self, -1, "Wavelength range reduction", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2))
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        wave_base = self.text2float(dlg.wave_base.GetValue(), 'Base wavelength value is not a valid one.')
        wave_top = self.text2float(dlg.wave_top.GetValue(), 'Top wavelength value is not a valid one.')
        dlg.Destroy()
        
        if wave_base == None or wave_top == None or wave_top <= wave_base:
            self.flash_status_message("Bad value.")
            return
        
        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to cut it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        
        self.status_message("Cutting spectra...")
        wfilter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
        self.active_spectrum.data = self.active_spectrum.data[wfilter]
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()
        
        self.remove_continuum_spectra()
        self.remove_fitted_lines()
        
        self.update_title()
        # Autoscale
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw()
        self.flash_status_message("Spectrum cutted.")
        
    def on_combine_spectra(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = HomogenizeCombineSpectraDialog(self, -1, "Homogenize & combine spectrum", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2), self.active_spectrum.resolution)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        wave_base = self.text2float(dlg.wave_base.GetValue(), 'Base wavelength value is not a valid one.')
        wave_top = self.text2float(dlg.wave_top.GetValue(), 'Top wavelength value is not a valid one.')
        resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
        dlg.Destroy()
        
        if wave_base == None or wave_top == None or wave_top <= wave_base:
            self.flash_status_message("Bad value.")
            return
        # It is not necessary to check if base and top are out of the current spectra,
        # the resample_spectra function can deal it
        
        # Check if there is any spectrum not saved
        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if spectra_not_saved:
            msg = "Some of the spectra have not been saved, are you sure you want to combine them all without saving previously?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_combine_spectra_thread, args=(wave_base, wave_top, resolution,))
        thread.setDaemon(True)
        thread.start()
        
    def on_combine_spectra_thread(self, wave_base, wave_top, resolution):
        # Homogenize
        i = 0
        total = len(self.spectra)
        for spec in self.spectra:
            wx.CallAfter(self.status_message, "Homogenizing spectrum %i of %i (%s)..." % (i+1, total, spec.name))
            xaxis = generate_wavelength_grid(wave_base, wave_top, resolution)
            resampled_spectra_data = resample_spectra(spec.data, xaxis, frame=self)
            spec.data = resampled_spectra_data
            
            # Remove plot
            self.axes.lines.remove(spec.plot_id)
            # Remove fitted lines if they exist
            for region in self.region_widgets["lines"]:
                if region.line_plot_id.has_key(spec):
                    if region.line_plot_id[spec] != None:
                        self.axes.lines.remove(region.line_plot_id[spec])
                    region.line_plot_id.remove(spec)
                    region.line_model.remove(spec)
                    region.continuum_base_level.remove(spec)
            i += 1
        
        # Combine
        wx.CallAfter(self.status_message, "Combining spectra...")
        wx.CallAfter(self.update_progress, 10)
        total_spectra = len(self.spectra)
        total_wavelengths = len(self.active_spectrum.data['waveobs'])
        matrix = np.empty((total_spectra, total_wavelengths))
        
        # Create a matrix for an easy combination
        i = 0
        for spec in self.spectra:
            matrix[i] = spec.data['flux']
            i += 1
        
        # Standard deviation
        std = np.std(matrix, axis=0)
        
        # Mean fluxes
        cumulative_mean_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        cumulative_mean_spectra['waveobs'] = self.active_spectrum.data['waveobs'] # Wavelengths from the reference spectra
        cumulative_mean_spectra['flux'] = np.mean(matrix, axis=0)
        cumulative_mean_spectra['err'] = std
        
        # Median fluxes
        cumulative_median_spectra = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        cumulative_median_spectra['waveobs'] = self.active_spectrum.data['waveobs']
        cumulative_median_spectra['flux'] = np.median(matrix, axis=0)
        cumulative_median_spectra['err'] = std
        
        wx.CallAfter(self.on_combine_spectra_finnish, cumulative_mean_spectra, cumulative_median_spectra)
    
    def on_combine_spectra_finnish(self, cumulative_mean_spectra, cumulative_median_spectra):
        # Add plots
        self.spectra = []
        
        name = self.get_name("Cumulative_mean_spectra") # If it already exists, add a suffix
        color = self.get_color()
        self.active_spectrum = Spectrum(cumulative_mean_spectra, name, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.draw_active_spectrum()
        
        # Remove "[A]  " from spectra name (legend)
        self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
        
        name = self.get_name("Cumulative_median_spectra") # If it already exists, add a suffix
        color = self.get_color()
        self.active_spectrum = Spectrum(cumulative_median_spectra, name, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.draw_active_spectrum()
        
        self.update_menu_active_spectrum()                        
        
        self.update_title()
        # Autoscale
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw()
        self.flash_status_message("Spectra combined.")
        self.operation_in_progress = False
        
        
    def on_continuum_normalization(self, event):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return
        
        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to divide all the fluxes by the fitted continuum anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            msg = "This operation is going to divide all the fluxes of the active spectrum by the fitted continuum. Are you sure?"
            title = "Confirmation"
            if not self.question(title, msg):
                return
        
        self.status_message("Continuum normalization...")
        self.active_spectrum.data['flux'] /= self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'])
        self.active_spectrum.not_saved = True
        
        # Remove current continuum from plot if exists
        self.remove_continuum_spectra()
        
        # Remove current drawn fitted lines if they exist
        self.remove_fitted_lines()
        
        self.draw_active_spectrum()
        self.update_title()
        # Autoscale
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw()
        self.flash_status_message("Spectrum normalized")
        pass
    
    def on_fit_lines(self, event):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return
        
        # Remove drawd lines if they exist
        self.remove_drawn_fitted_lines()
        
        self.operation_in_progress = True
        self.status_message("Fitting lines...")
        thread = threading.Thread(target=self.on_fit_lines_thread)
        thread.setDaemon(True)
        thread.start()
    
    def on_fit_lines_thread(self):
        total_regions = len(self.region_widgets["lines"])
        i = 0
        for region in self.region_widgets["lines"]:
            wave_base = region.get_wave_base()
            wave_top = region.get_wave_top()
            mu = region.get_wave_peak()
            wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
            spectra_window = self.active_spectrum.data[wave_filter]
        
            try:
                # Remove old possible results
                region.line_model[self.active_spectrum] = None
                region.continuum_base_level[self.active_spectrum] = None
                
                # Fit
                line_model = fit_line(spectra_window, self.active_spectrum.continuum_mode, mu)
                continuum_value = continuum_model(spectra_window['waveobs'])
                
                # Save results
                region.line_model[self.active_spectrum] = line_model
                region.continuum_base_level[self.active_spectrum] = np.mean(continuum_value)
            except Exception as e:
                print "Error:", wave_base, wave_top, e
            
            current_work_progress = (i*1.0 / total_regions) * 100
            wx.CallAfter(self.update_progress, current_work_progress)
            i += 1
        
        wx.CallAfter(self.on_fit_lines_finnish)
    
    def on_fit_lines_finnish(self):
        # Draw new lines
        self.draw_fitted_lines()
        self.operation_in_progress = False
        self.flash_status_message("Lines fitted.")

    
    def remove_drawn_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_plot_id.has_key(self.active_spectrum) and region.line_plot_id[self.active_spectrum] != None:
                self.axes.lines.remove(region.line_plot_id[self.active_spectrum])
                region.line_plot_id[self.active_spectrum] = None
        self.canvas.draw()
        
    def remove_fitted_lines(self):
        self.remove_drawn_fitted_lines()
        # Remove fitted lines if they exist
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum):
                del region.line_plot_id[self.active_spectrum]
                del region.line_model[self.active_spectrum]
                del region.continuum_base_level[self.active_spectrum]
    
    def draw_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum) and region.line_model[self.active_spectrum] != None:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()
                wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
                spectra_window = self.active_spectrum.data[wave_filter]
                
                # Get fluxes from model
                line_fluxes = region.line_model[self.active_spectrum](spectra_window['waveobs'])
                # Add the continuum base because the line_model has substracted it
                line_fluxes += region.continuum_base_level[self.active_spectrum]
                
                # zorder = 4, above the line region
                line_plot_id = self.axes.plot(spectra_window['waveobs'], line_fluxes, lw=1, color='red', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=4)[0]
                region.line_plot_id[self.active_spectrum] = line_plot_id
        self.canvas.draw()
    
    def on_change_active_spectrum(self, event):
        if self.check_operation_in_progress():
            return
        
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
    
    def on_remove_fitted_continuum(self, event):
        if self.check_operation_in_progress():
            return
        self.remove_continuum_spectra()
    
    def on_fit_continuum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        # Initial recommendation: 1 knot every 10 nm
        nknots = np.max([1, int((np.max(self.active_spectrum.data['waveobs']) - np.min(self.active_spectrum.data['waveobs'])) / 10)])
        
        dlg = FitContinuumDialog(self, -1, "Properties for fitting continuum", nknots=nknots)
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
        
        
        self.operation_in_progress = True
        self.status_message("Fitting continuum...")
        self.update_progress(10)
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
            try:
                self.active_spectrum.continuum_model = fit_continuum(spectra_regions, nknots=nknots)
                self.active_spectrum.continuum_data = get_spectra_from_model(self.active_spectrum.continuum_model, self.active_spectrum.data['waveobs'])
                wx.CallAfter(self.on_fit_continuum_finish, nknots)
            except Exception as e:
                self.operation_in_progress = False
                wx.CallAfter(self.flash_status_message, "Not enough data points to fit, reduce the number of nknots or increase the spectra regions.")
        else:
            self.operation_in_progress = False
            wx.CallAfter(self.flash_status_message, "No continuum regions found.")
    
    def on_fit_continuum_finish(self, nknots):
        self.draw_continuum_spectra()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Continuum fitted with %s knots uniform Spline model." % str(nknots))
    
    def on_find_continuum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return
        
        if self.not_saved["continuum"]:
            msg = "Are you sure you want to find new continuum regions without saving the current ones?"
            title = "Changes not saved"
            if not self.question(title, msg):
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
        
        self.operation_in_progress = True
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
        self.operation_in_progress = False
        self.flash_status_message("Automatic finding of continuum regions ended.")
    
    def on_determine_barycentric_vel(self, event):
        dlg = DetermineBarycentricCorrectionDialog(self, -1, "Barycentric velocity determination", self.day, self.month, self.year, self.hours, self.minutes, self.seconds, self.ra_hours, self.ra_minutes, self.ra_seconds, self.dec_degrees, self.dec_minutes, self.dec_seconds)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        day = self.text2float(dlg.day.GetValue(), 'Day value is not a valid one.')
        month = self.text2float(dlg.month.GetValue(), 'Month value is not a valid one.')
        year = self.text2float(dlg.year.GetValue(), 'Year value is not a valid one.')
        
        hours = self.text2float(dlg.hours.GetValue(), 'Hours value is not a valid one.')
        minutes = self.text2float(dlg.minutes.GetValue(), 'Minutes value is not a valid one.')
        seconds = self.text2float(dlg.seconds.GetValue(), 'Seconds value is not a valid one.')
        
        ra_hours = self.text2float(dlg.ra_hours.GetValue(), 'RA hours value is not a valid one.')
        ra_minutes = self.text2float(dlg.ra_minutes.GetValue(), 'RA minutes value is not a valid one.')
        ra_seconds = self.text2float(dlg.ra_seconds.GetValue(), 'RA seconds value is not a valid one.')
        
        dec_degrees = self.text2float(dlg.dec_degrees.GetValue(), 'DEC degrees value is not a valid one.')
        dec_minutes = self.text2float(dlg.dec_minutes.GetValue(), 'DEC minutes value is not a valid one.')
        dec_seconds = self.text2float(dlg.dec_seconds.GetValue(), 'DEC seconds value is not a valid one.')
        
        if None in [day, month, year, hours, minutes, seconds, ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds]:
            self.flash_status_message("Bad value.")
            return
        
        self.day = day
        self.month = month
        self.year = year
        self.hours = hours
        self.minutes = minutes
        self.seconds = seconds
        self.ra_hours = ra_hours
        self.ra_minutes = ra_minutes
        self.ra_seconds = ra_seconds
        self.dec_degrees = dec_degrees
        self.dec_minutes = dec_minutes
        self.dec_seconds = dec_seconds
        
        vh, vb = baryvel((year, month, day, hours, minutes, seconds))

        ra = (ra_hours + ra_minutes/60 + ra_seconds/(60*60)) # hours
        ra = ra * 360/24 # degrees
        ra = ra * ((2*np.pi) / 360) # radians
        dec = (dec_degrees + dec_minutes/60 + dec_seconds/(60*60)) # degrees
        dec = dec * ((2*np.pi) / 360) # radians
        
        # Project velocity toward star
        self.barycentric_vel = vb[0]*np.cos(dec)*np.cos(ra) + vb[1]*np.cos(dec)*np.sin(ra) + vb[2]*np.sin(dec) # km/s
        self.barycentric_vel = np.round(self.barycentric_vel, 2) # km/s
        
        self.flash_status_message("Barycentric velocity determined: " + str(self.barycentric_vel) + " km/s")

        
    def on_determine_rv(self, event, show_previous_results=True):
        if not self.check_active_spectrum_exists():
            return
        
        if self.active_spectrum.rv_profile_model != None and show_previous_results:
            dlg = RVProfileDialog(self, -1, "Radial velocity profile", self.active_spectrum.rv_profile_xcoord, self.active_spectrum.rv_profile_fluxes, self.active_spectrum.rv_profile_model, self.active_spectrum.rv_profile_num_used_lines, self.active_spectrum.rv_profile_rv_step)
            dlg.ShowModal()
            recalculate = dlg.recalculate
            dlg.Destroy()
            if not recalculate:
                return
        
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = DetermineRVDialog(self, -1, "Radial velocity determination", rv_lower_limit=self.rv_lower_limit, rv_upper_limit=self.rv_upper_limit, rv_step=self.rv_step)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        rv_lower_limit = self.text2float(dlg.rv_lower_limit.GetValue(), 'Radial velocity lower limit is not a valid one.')
        rv_upper_limit = self.text2float(dlg.rv_upper_limit.GetValue(), 'Radial velocity upper limit is not a valid one.')
        rv_step = self.text2float(dlg.rv_step.GetValue(), 'Radial velocity step is not a valid one.')
        renormalize = True #dlg.renormalize.GetValue()
        dlg.Destroy()
        
        if rv_lower_limit == None or rv_upper_limit == None or rv_step == None:
            self.flash_status_message("Bad value.")
            return
        
        if rv_lower_limit >= rv_upper_limit:
            msg = "Upper radial velocity limit should be greater than lower limit"
            title = "Radial velocity error"
            self.error(title, msg)
            return
        
        if (np.abs(rv_lower_limit) + np.abs(rv_upper_limit)) <= 4*rv_step:
            msg = "Radial velocity step too small for the established limits"
            title = "Radial velocity error"
            self.error(title, msg)
            return
        
        self.rv_lower_limit = rv_lower_limit
        self.rv_upper_limit = rv_upper_limit
        self.rv_step = rv_step
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_determine_rv_thread, args=(renormalize,))
        thread.setDaemon(True)
        thread.start()

    def on_determine_rv_thread(self, renormalize):
        wx.CallAfter(self.status_message, "Determining RV...")
        if self.linelist_rv == None:
            self.linelist_rv = asciitable.read("input/rv/default.300_1100nm.rv.lst")
        xcoord, fluxes, num_used_lines = build_radial_velocity_profile(self.active_spectrum.data, self.active_spectrum.continuum_model, self.linelist_rv, rv_lower_limit=self.rv_lower_limit, rv_upper_limit=self.rv_upper_limit, rv_step=self.rv_step, frame=self)
        wx.CallAfter(self.on_determine_rv_finish, xcoord, fluxes, renormalize, num_used_lines)
    
    def on_determine_rv_finish(self, xcoord, fluxes, renormalize, num_used_lines):
        # Modelize
        model, fluxes = model_radial_velocity_profile(xcoord, fluxes, renormalize=renormalize)
        self.active_spectrum.rv = np.round(model.mu, 2) # km/s
        # A positive radial velocity indicates the distance between the objects is or was increasing;
        # A negative radial velocity indicates the distance between the source and observer is or was decreasing.
        self.flash_status_message("Radial velocity determined: " + str(self.active_spectrum.rv) + " km/s")
        self.operation_in_progress = False
        
        self.active_spectrum.rv_profile_xcoord = xcoord
        self.active_spectrum.rv_profile_fluxes = fluxes
        self.active_spectrum.rv_profile_model = model
        self.active_spectrum.rv_profile_num_used_lines = num_used_lines
        self.active_spectrum.rv_profile_rv_step = self.rv_step
        
        dlg = RVProfileDialog(self, -1, "Radial velocity profile", xcoord, fluxes, model, num_used_lines, self.rv_step)
        dlg.ShowModal()
        recalculate = dlg.recalculate
        dlg.Destroy()
        if recalculate:
            self.on_determine_rv(None, show_previous_results=False)


    def on_correct_barycentric_vel(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = CorrectRVDialog(self, -1, "Barycentric velocity correction", self.barycentric_vel)
        dlg.ShowModal()
        
        if not dlg.action_accepted:
            dlg.Destroy()
            return
        
        barycentric_vel = self.text2float(dlg.rv.GetValue(), 'Radial velocity value is not a valid one.')
        in_regions = dlg.radio_button_regions.GetValue()
        dlg.Destroy()
        
        if barycentric_vel == None:
            self.flash_status_message("Bad value.")
            return
        
        self.status_message("Correcting barycentric velocity...")
        if in_regions:
            for elements in ["lines", "continuum", "segments"]:
                self.update_numpy_arrays_from_widgets(elements)
                if len(self.regions[elements]) > 0:
                    if elements == "lines":
                        self.regions[elements] = correct_radial_velocity_regions(self.regions[elements], barycentric_vel, with_peak=True)
                    else:
                        self.regions[elements] = correct_radial_velocity_regions(self.regions[elements], barycentric_vel)
                    self.draw_regions(elements)
                    self.not_saved[elements] = False
        else:
            if not self.check_active_spectrum_exists():
                return
            self.active_spectrum.data = correct_radial_velocity(self.active_spectrum.data, barycentric_vel)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
        self.remove_continuum_spectra()
        self.remove_fitted_lines()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a radial velocity correction of %s." % barycentric_vel)
    
    
    def on_correct_rv(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        dlg = CorrectRVDialog(self, -1, "Radial velocity correction", self.active_spectrum.rv)
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
        
        self.status_message("Correcting radial velocity...")
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
            if not self.check_active_spectrum_exists():
                return
            self.active_spectrum.data = correct_radial_velocity(self.active_spectrum.data, radial_vel)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
        self.remove_continuum_spectra()
        self.remove_fitted_lines()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a radial velocity correction of %s." % radial_vel)
    
    def on_convert_to_nm(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        
        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to divide all the wavelength by 10 anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            msg = "This operation is going to divide all the wavelengths of the active spectrum by 10. Are you sure?"
            title = "Confirmation"
            if not self.question(title, msg):
                return
        
        self.status_message("Converting to nanometers...")
        self.active_spectrum.data['waveobs'] = self.active_spectrum.data['waveobs'] / 10
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()
        
        self.remove_continuum_spectra()
        self.remove_fitted_lines()
        
        self.update_title()
        # Autoscale
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw()
        self.flash_status_message("Wavelength divided by 10.")
    
    def on_synthesize(self, event):
        if self.check_operation_in_progress():
            return
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
            
            
            if self.modeled_layers_pack == None:
                print "Loading modeled atmospheres..."
                self.status_message("Loading modeled atmospheres...")
                self.modeled_layers_pack = load_modeled_layers_pack(filename='input/atmospheres/default.modeled_layers_pack.dump')
            
            if not valid_objective(self.modeled_layers_pack, teff, logg, MH):
                dlg_error = wx.MessageDialog(self, "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models.", 'Out of the atmospheric models', wx.OK | wx.ICON_ERROR)
                dlg_error.ShowModal()
                dlg_error.Destroy()
                self.flash_status_message("Bad values.")
                return
            
            # Prepare atmosphere model
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
            
            ## TODO: Control wavelength margins, they should not be bigger/lower thant the linelist
            
            self.operation_in_progress = True
            self.status_message("Synthesizing spectrum...")
            self.update_progress(10)
            thread = threading.Thread(target=self.on_synthesize_thread, args=(waveobs, atm_filename, teff, logg, MH, microturbulence_vel,))
            thread.setDaemon(True)
            thread.start()
        
    def on_synthesize_thread(self, waveobs, atm_filename, teff, logg, MH, microturbulence_vel):
        total_points = len(waveobs)
        
        synth_spectra = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        synth_spectra['waveobs'] = waveobs
        
        synth_spectra['flux'] = synthesizer.spectrum(synth_spectra['waveobs']*10.0, atm_filename, linelist_file = "input/linelists/default.300_1100nm.lst", abundances_file = "input/abundances/default.stdatom.dat", microturbulence_vel = microturbulence_vel, verbose=1, update_progress_func=self.update_progress)
            
        
        synth_spectra.sort(order='waveobs') # Make sure it is ordered by wavelength
        
        # Remove atmosphere model temporary file
        os.remove(atm_filename)
        wx.CallAfter(self.on_synthesize_finnish, synth_spectra, teff, logg, MH, microturbulence_vel)
    
    def on_synthesize_finnish(self, synth_spectra, teff, logg, MH, microturbulence_vel):
        # Draw new lines
        self.draw_fitted_lines()
        self.operation_in_progress = False
        self.flash_status_message("Lines fitted.")
        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectra()
        
        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()
        
        # Remove "[A]  " from spectra name (legend) if it exists
        if self.active_spectrum != None and self.active_spectrum.plot_id != None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
        
        # Name: If it already exists, add a suffix
        name = self.get_name("Synth_" + str(teff) + "_" + str(logg) + "_"  + str(MH) + "_" + str(microturbulence_vel))
        color = self.get_color()
        self.active_spectrum = Spectrum(synth_spectra, name, color=color)
        
        self.spectra.append(self.active_spectrum)
        self.active_spectrum.not_saved = True
        self.update_title()
        self.update_menu_active_spectrum()                        
        self.draw_active_spectrum()
        
        self.canvas.draw()

        self.operation_in_progress = False
        self.flash_status_message("Synthetic spectra generated!")

    
    def on_exit(self, event):
        self.on_close(event)
        
    def on_about(self, event):
        description = """Spectra Visual Editor is a tool for the treatment of spectrum files in order to identify lines, continuum regions and determine radial velocities among other options.
"""
        if sys.modules.has_key('synthesizer'):
            description += """
The generation of synthetic spectrum is done thanks to:

SPECTRUM a Stellar Spectral Synthesis Program 
(C) Richard O. Gray 1992 - 2010 Version 2.76e
"""

        licence = """Spectra Visual Editor is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Spectra Visual Editor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Spectra Visual Editor.  If not, see <http://www.gnu.org/licenses/>."""


        info = wx.AboutDialogInfo()

        info.SetIcon(wx.Icon('images/SVE.png', wx.BITMAP_TYPE_PNG))
        info.SetName('Spectra Visual Editor')
        #info.SetVersion('2012.02.22')
        info.SetDescription(description)
        info.SetCopyright('(C) 2011 - 2012 Sergi Blanco Cuaresma')
        info.SetWebSite('http://www.marblestation.com')
        info.SetLicence(licence)
        info.AddDeveloper('Sergi Blanco Cuaresma')
        #info.AddDocWriter('Jan Bodnar')
        #info.AddArtist('The Tango crew')
        #info.AddTranslator('Jan Bodnar')

        wx.AboutBox(info)
    
    def flash_status_message(self, msg, flash_len_ms=3000, progress=True):
        if self.timeroff.IsRunning():
            self.timeroff.Stop()
        
        self.statusbar.SetStatusText(msg)
        if progress:
            self.gauge.SetValue(pos=100)
        
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
        
    def status_message(self, msg):
        if self.timeroff.IsRunning():
            self.timeroff.Stop()
            # Progress bar to zero
            self.gauge.SetValue(pos=0)
        
        self.statusbar.SetStatusText(msg)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
        # Progress bar to zero
        self.gauge.SetValue(pos=0)


## Print usage
def usage():
    print "Usage:"
    print sys.argv[0], "[--continuum=file] [--lines=file] [--segments=file] [spectra_file]"

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
    
    filenames = {}
    filenames['spectra'] = []
    filenames['continuum'] = continuum_file
    filenames['lines'] = lines_file
    filenames['segments'] = segments_file
    
    # Open spectra
    for arg in args:
        spectra_file = arg
        if not os.path.exists(spectra_file):
            print "Spectra file", arg, "does not exists!"
            sys.exit(2)
        filenames['spectra'].append(spectra_file)

    return filenames

#~ Example:
#~ python interactive.py --continuum=input/LUMBA/UVES_MRD_sun_cmask.txt --lines=input/LUMBA/UVES_MRD_sun_Fe-linelist.txt --segments=input/LUMBA/UVES_MRD_sun_segments.txt input/LUMBA/UVES_MRD_sun_official.s.gz

## Input files should be '\t' separated, comments with # and the first line is the header
if __name__ == '__main__':
    filenames = get_arguments()    
    
    #### Read files

    spectra = []
    for path in filenames['spectra']:
        try:
            spectrum = read_spectra(path)
            #wfilter = (spectrum['waveobs'] >= 516.0) & (spectrum['waveobs'] <= 519.0)
            #spectrum = spectrum[wfilter]
        except Exception as e:
            print "Spectra file", path, "has an incompatible format!"
            sys.exit(2)
        spectra.append(spectrum)
    
    if filenames['continuum'] != None:
        try:
            continuum = read_continuum_regions(filenames['continuum'])
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
            lines = read_line_regions(filenames['lines'])
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
            segments = read_segment_regions(filenames['segments'])
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
    app.frame = SpectraFrame(spectra, regions, filenames)
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

