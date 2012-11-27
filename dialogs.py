#!/usr/bin/env python
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
#################
# Run with ipython -pdb -c "%run interactive.py"
#################
#import ipdb

import os
import sys

import wx
import numpy as np
import scipy.stats as stats

# The recommended way to use wx with mpl is with the WXAgg backend.
import matplotlib
#matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.ticker import ScalarFormatter

class FitContinuumDialog(wx.Dialog):
    def __init__(self, parent, id, title, nknots=1, median_wave_range=0.1, max_wave_range=1):
        wx.Dialog.__init__(self, parent, id, title, size=(400,300))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Splines
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_nknots = wx.StaticText(self, -1, "Number of splines: ", style=wx.ALIGN_LEFT)
        self.nknots = wx.TextCtrl(self, -1, str(nknots),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_nknots, 0, border=3, flag=flags)
        self.hbox.Add(self.nknots, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)


        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_median_wave_range = wx.StaticText(self, -1, "Wavelength step for median selection: ", style=wx.ALIGN_LEFT)
        self.median_wave_range = wx.TextCtrl(self, -1, str(median_wave_range),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_median_wave_range, 0, border=3, flag=flags)
        self.hbox.Add(self.median_wave_range, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)


        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_max_wave_range = wx.StaticText(self, -1, "Wavelength step for max selection: ", style=wx.ALIGN_LEFT)
        self.max_wave_range = wx.TextCtrl(self, -1, str(max_wave_range),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_max_wave_range, 0, border=3, flag=flags)
        self.hbox.Add(self.max_wave_range, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.hbox.AddSpacer(10)

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
        self.EndModal(wx.ID_OK)

class EstimateSNRDialog(wx.Dialog):
    def __init__(self, parent, id, title, num_points=10):
        wx.Dialog.__init__(self, parent, id, title)

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Number of points
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_num_points = wx.StaticText(self, -1, "Number of points: ", style=wx.ALIGN_LEFT)
        self.num_points = wx.TextCtrl(self, -1, str(num_points),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_num_points, 0, border=3, flag=flags)
        self.hbox.Add(self.num_points, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Estimate SNR: ", style=wx.ALIGN_LEFT)
        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_where, 0, border=3, flag=flags)

        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.radio_button_from_errors = wx.RadioButton(self, -1, 'Directly from reported errors', style=wx.RB_GROUP)
        self.vbox2.Add(self.radio_button_from_errors, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_from_flux = wx.RadioButton(self, -1, 'From fluxes in blocks of N points')
        self.vbox2.Add(self.radio_button_from_flux, 0, border=3, flag=wx.LEFT | wx.TOP | wx.GROW)
        self.radio_button_from_flux.SetValue(True)
        self.hbox.Add(self.vbox2, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)


class FindContinuumDialog(wx.Dialog):
    def __init__(self, parent, id, title, fixed_wave_step=0.05, sigma=0.001, max_continuum_diff=1.0):
        wx.Dialog.__init__(self, parent, id, title, size=(450,350))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Max step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_fixed_wave_step = wx.StaticText(self, -1, "Check for regions of minimum size: ", style=wx.ALIGN_LEFT)
        self.fixed_wave_step = wx.TextCtrl(self, -1, str(fixed_wave_step),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_fixed_wave_step, 0, border=3, flag=flags)
        self.hbox.Add(self.fixed_wave_step, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_sigma = wx.StaticText(self, -1, "Maximum standard deviation: ", style=wx.ALIGN_LEFT)
        self.sigma = wx.TextCtrl(self, -1, str(sigma),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_sigma, 0, border=3, flag=flags)
        self.hbox.Add(self.sigma, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_max_continuum_diff = wx.StaticText(self, -1, "Maximum fitted continuum difference (%): ", style=wx.ALIGN_LEFT)
        self.max_continuum_diff = wx.TextCtrl(self, -1, str(max_continuum_diff),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_max_continuum_diff, 0, border=3, flag=flags)
        self.hbox.Add(self.max_continuum_diff, 0, border=3, flag=flags)

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
        self.EndModal(wx.ID_OK)

class FindLinesDialog(wx.Dialog):
    def __init__(self, parent, id, title, min_depth=0.05, max_depth=1.0, vel_atomic=0.0, vel_telluric=0.0, resolution=300000, elements="Fe 1, Fe 2"):
        wx.Dialog.__init__(self, parent, id, title, size=(450,450))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Min Depth
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_min_depth = wx.StaticText(self, -1, "Minimum depth (% of the continuum): ", style=wx.ALIGN_LEFT)
        self.min_depth = wx.TextCtrl(self, -1, str(min_depth),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_min_depth, 0, border=3, flag=flags)
        self.hbox.Add(self.min_depth, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Max Depth
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_max_depth = wx.StaticText(self, -1, "Maximum depth (% of the continuum): ", style=wx.ALIGN_LEFT)
        self.max_depth = wx.TextCtrl(self, -1, str(max_depth),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_max_depth, 0, border=3, flag=flags)
        self.hbox.Add(self.max_depth, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Elements
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_elements = wx.StaticText(self, -1, "Select elements (comma separated): ", style=wx.ALIGN_LEFT)
        self.elements = wx.TextCtrl(self, -1, str(elements),  style=wx.TE_LEFT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_elements, 0, border=3, flag=flags)
        self.hbox.Add(self.elements, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Resolution
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_resolution = wx.StaticText(self, -1, "Resolution: ", style=wx.ALIGN_LEFT)
        self.resolution = wx.TextCtrl(self, -1, str(resolution),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_resolution, 0, border=3, flag=flags)
        self.hbox.Add(self.resolution, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Velocity respect to atomic data
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_vel_atomic = wx.StaticText(self, -1, "Velocity respect to atomic lines (km/s): ", style=wx.ALIGN_LEFT)
        self.vel_atomic = wx.TextCtrl(self, -1, str(vel_atomic),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_vel_atomic, 0, border=3, flag=flags)
        self.hbox.Add(self.vel_atomic, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Velocity respect to telluric data
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_vel_telluric = wx.StaticText(self, -1, "Velocity respect to telluric lines (km/s): ", style=wx.ALIGN_LEFT)
        self.vel_telluric = wx.TextCtrl(self, -1, str(vel_telluric),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_vel_telluric, 0, border=3, flag=flags)
        self.hbox.Add(self.vel_telluric, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Standard deviation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.discard_tellurics = wx.CheckBox(self, -1, 'Discard affected by tellurics', style=wx.ALIGN_LEFT)
        self.discard_tellurics.SetValue(True)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.discard_tellurics, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)


        ### Where to look
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_where = wx.StaticText(self, -1, "Look for line masks in: ", style=wx.ALIGN_LEFT)
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
        self.EndModal(wx.ID_OK)

class CorrectVelocityDialog(wx.Dialog):
    def __init__(self, parent, id, vel_type, rv):
        title = vel_type.capitalize() + " velocity correction"
        wx.Dialog.__init__(self, parent, id, title)

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### RV
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_rv = wx.StaticText(self, -1, vel_type.capitalize() + " velocity (km/s): ", style=wx.ALIGN_LEFT)
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
        self.radio_button_spectra = wx.RadioButton(self, -1, 'Spectrum', style=wx.RB_GROUP)
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
        self.EndModal(wx.ID_OK)

class DetermineVelocityDialog(wx.Dialog):
    def __init__(self, parent, id, title, rv_upper_limit, rv_lower_limit, rv_step, templates):
        wx.Dialog.__init__(self, parent, id, title)

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### RV lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_rv_lower_limit = wx.StaticText(self, -1, "Velocity lower limit (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_lower_limit = wx.TextCtrl(self, -1, str(int(rv_lower_limit)),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_lower_limit, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_lower_limit, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### RV upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_rv_upper_limit = wx.StaticText(self, -1, "Velocity upper limit (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_upper_limit = wx.TextCtrl(self, -1, str(int(rv_upper_limit)),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_upper_limit, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_upper_limit, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### RV step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_rv_step = wx.StaticText(self, -1, "Velocity steps (km/s): ", style=wx.ALIGN_LEFT)
        self.rv_step = wx.TextCtrl(self, -1, str(rv_step),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv_step, 0, border=3, flag=flags)
        self.hbox.Add(self.rv_step, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Template
        if len(templates) > 0:
            self.hbox = wx.BoxSizer(wx.HORIZONTAL)

            self.text_templates = wx.StaticText(self, -1, "Cross-correlate with: ", style=wx.ALIGN_LEFT)
            self.templates = wx.ComboBox (self, wx.ID_ANY, templates[0], choices=templates, style=wx.CB_READONLY)

            self.hbox.AddSpacer(10)
            self.hbox.Add(self.text_templates, 0, border=3, flag=flags)
            self.hbox.Add(self.templates, 0, border=3, flag=flags)

            self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

            self.vbox.AddSpacer(30)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

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
        self.EndModal(wx.ID_OK)


class VelocityProfileDialog(wx.Dialog):
    def __init__(self, parent, id, title, xcoord, fluxes, errors, models, rv_step, telluric_fwhm=0.0, snr=0.0, template=None):
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
      	self.axes.set_ylim([0,1])
        # Avoid using special notation that are not easy to understand in axis for big zoom
        myyfmt = ScalarFormatter(useOffset=False)
    	self.axes.get_xaxis().set_major_formatter(myyfmt)
    	self.axes.get_yaxis().set_major_formatter(myyfmt)

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
        self.axes.plot(xcoord, fluxes, lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)
        self.axes.fill_between(xcoord, fluxes+errors, fluxes-errors, color='#CCCCCC')

        if rv_step >= 0.1:
            xcoord_mod = np.arange(np.min(xcoord), np.max(xcoord), 0.1)
            for model in models:
                self.axes.plot(xcoord_mod, model(xcoord_mod), lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)
        else:
            for model in models:
                self.axes.plot(xcoord, model(xcoord), lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)

        self.axes.grid(True, which="both")
        self.axes.set_title("Profile", fontsize="10")
        self.axes.set_xlabel("velocity (km/s)", fontsize="10")
        self.axes.set_ylabel("relative intensity", fontsize="10")

        ## Stats
        num_items = self.stats.GetItemCount()
        if snr != 0:
            self.stats.InsertStringItem(num_items, "Estimated local SNR")
            self.stats.SetStringItem(num_items, 1, str(np.round(snr, 2)))
            num_items += 1
        if template != None:
            self.stats.InsertStringItem(num_items, "Template")
            self.stats.SetStringItem(num_items, 1, template)
            num_items += 1
        for i, model in enumerate(models):
            self.stats.InsertStringItem(num_items, "Mean (km/s)")
            self.stats.SetStringItem(num_items, 1, str(np.round(model.mu(), 2)))
            num_items += 1

            self.stats.InsertStringItem(num_items, "Min. error (+/- km/s)")
            self.stats.SetStringItem(num_items, 1, str(np.round(rv_step/2, 4)))
            num_items += 1

            self.stats.InsertStringItem(num_items, "Baseline")
            self.stats.SetStringItem(num_items, 1, str(np.round(model.baseline(), 2)))
            num_items += 1
            self.stats.InsertStringItem(num_items, "A (rel. intensity)")
            self.stats.SetStringItem(num_items, 1, str(np.round(model.A(), 2)))
            num_items += 1
            self.stats.InsertStringItem(num_items, "Sigma (km/s)")
            self.stats.SetStringItem(num_items, 1, str(np.round(model.sig(), 2)))
            num_items += 1

            try:
                # If model is VoigtModel
                self.stats.InsertStringItem(num_items, "Gamma")
                self.stats.SetStringItem(num_items, 1, str(np.round(model.gamma(), 2)))
                num_items += 1
            except AttributeError:
                # model is GaussianModel
                pass

            fwhm = model.fwhm()[0] # km/s (because xcoord is already velocity)
            self.stats.InsertStringItem(num_items, "Measured FWHM (km/s)")
            self.stats.SetStringItem(num_items, 1, str(np.round(fwhm, 2)))
            num_items += 1
            if telluric_fwhm != 0:
                self.stats.InsertStringItem(num_items, "Theoretical telluric FWHM (km/s)")
                self.stats.SetStringItem(num_items, 1, str(np.round(telluric_fwhm, 2)))
                num_items += 1
                self.stats.InsertStringItem(num_items, "Corrected FWHM (km/s)")
                self.stats.SetStringItem(num_items, 1, str(np.round(fwhm - telluric_fwhm, 2)))
                num_items += 1
            self.stats.InsertStringItem(num_items, "Estimated resolving power (R)")
            c = 299792458.0 # m/s
            R = np.int(c/(1000.0*fwhm - telluric_fwhm))
            self.stats.SetStringItem(num_items, 1, str(np.round(R, 2)))
            num_items += 1

            rms = model.rms
            self.stats.InsertStringItem(num_items, "RMS")
            self.stats.SetStringItem(num_items, 1, str(np.round(rms, 5)))
            num_items += 1
            self.stats.InsertStringItem(num_items, "-------------------------------------------------------")
            self.stats.SetStringItem(num_items, 1, "---------------")
            num_items += 1

    def on_no(self, event):
        self.recalculate = False
        self.EndModal(wx.ID_NO)

    def on_yes(self, event):
        self.recalculate = True
        self.EndModal(wx.ID_YES)

class CleanSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, flux_base, flux_top, err_base, err_top):
        wx.Dialog.__init__(self, parent, id, title)

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Flux filter
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.filter_by_flux = wx.CheckBox(self, -1, 'Filter by flux', style=wx.ALIGN_LEFT)
        self.filter_by_flux.SetValue(True)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.filter_by_flux, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### flux lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_flux_base = wx.StaticText(self, -1, "Base flux: ", style=wx.ALIGN_LEFT)
        self.flux_base = wx.TextCtrl(self, -1, str(flux_base),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_flux_base, 0, border=3, flag=flags)
        self.hbox.Add(self.flux_base, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### flux upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_flux_top = wx.StaticText(self, -1, "Top flux: ", style=wx.ALIGN_LEFT)
        self.flux_top = wx.TextCtrl(self, -1, str(flux_top),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_flux_top, 0, border=3, flag=flags)
        self.hbox.Add(self.flux_top, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Error filter
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.filter_by_error = wx.CheckBox(self, -1, 'Filter by error', style=wx.ALIGN_LEFT)
        self.filter_by_error.SetValue(True)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.filter_by_error, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### err lower limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_err_base = wx.StaticText(self, -1, "Base error: ", style=wx.ALIGN_LEFT)
        self.err_base = wx.TextCtrl(self, -1, str(err_base),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_err_base, 0, border=3, flag=flags)
        self.hbox.Add(self.err_base, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### err upper limit
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_err_top = wx.StaticText(self, -1, "Top error: ", style=wx.ALIGN_LEFT)
        self.err_top = wx.TextCtrl(self, -1, str(err_top),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_err_top, 0, border=3, flag=flags)
        self.hbox.Add(self.err_top, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(10)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)


    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class CleanTelluricsDialog(wx.Dialog):
    def __init__(self, parent, id, title, rv, min_vel, max_vel, min_depth):
        wx.Dialog.__init__(self, parent, id, title)

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        ### RV
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_rv = wx.StaticText(self, -1, "Radial velocity: ", style=wx.ALIGN_LEFT)
        self.rv = wx.TextCtrl(self, -1, str(rv),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_rv, 0, border=3, flag=flags)
        self.hbox.Add(self.rv, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Min vel
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_min_vel = wx.StaticText(self, -1, "Minimum velocity: ", style=wx.ALIGN_LEFT)
        self.min_vel = wx.TextCtrl(self, -1, str(min_vel),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_min_vel, 0, border=3, flag=flags)
        self.hbox.Add(self.min_vel, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Max vel
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_max_vel = wx.StaticText(self, -1, "Maximum velocity: ", style=wx.ALIGN_LEFT)
        self.max_vel = wx.TextCtrl(self, -1, str(max_vel),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_max_vel, 0, border=3, flag=flags)
        self.hbox.Add(self.max_vel, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Min depth
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_min_depth = wx.StaticText(self, -1, "Minimum tellurics depth: ", style=wx.ALIGN_LEFT)
        self.min_depth = wx.TextCtrl(self, -1, str(min_depth),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_min_depth, 0, border=3, flag=flags)
        self.hbox.Add(self.min_depth, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(10)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)


    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)


class CutSpectrumDialog(wx.Dialog):
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
        self.EndModal(wx.ID_OK)

class ResampleSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top, wave_step):
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

        ### Wave step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_wave_step = wx.StaticText(self, -1, "Wavelength step: ", style=wx.ALIGN_LEFT)
        self.wave_step = wx.TextCtrl(self, -1, str(wave_step),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_step, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_step, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(10)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)


    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class CombineSpectraDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top, wave_step):
        wx.Dialog.__init__(self, parent, id, title, size=(400,450))

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

        ### Wave step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_wave_step = wx.StaticText(self, -1, "Wavelength step: ", style=wx.ALIGN_LEFT)
        self.wave_step = wx.TextCtrl(self, -1, str(wave_step),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_step, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_step, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        # Operation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_operation = wx.StaticText(self, -1, "Method: ", style=wx.ALIGN_LEFT)
        self.operation = wx.ComboBox (self, wx.ID_ANY, "Median", choices=["Median", "Mean", "Subtract", "Add", "Divide"], style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_operation, 0, border=3, flag=flags)
        self.hbox.Add(self.operation, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(30)


        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)


    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class OperateSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, operations_description, waveobs, flux, err):
        wx.Dialog.__init__(self, parent, id, title, size=(700,200))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        # Operation
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_operation = wx.StaticText(self, -1, "Available functions: ", style=wx.ALIGN_LEFT)
        self.operation = wx.ComboBox (self, wx.ID_ANY, operations_description[0], choices=operations_description, style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_operation, 0, border=3, flag=flags)
        self.hbox.Add(self.operation, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)


        #### Operate waveobs
        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.operate_waveobs = wx.CheckBox(self, -1, 'Modify wavelengths', style=wx.ALIGN_LEFT)
        #self.operate_waveobs.SetValue(True)

        #self.hbox.AddSpacer(10)
        #self.hbox.Add(self.operate_waveobs, 0, border=3, flag=flags)

        #self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Value
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_waveobs = wx.StaticText(self, -1, "Waves\t = ", style=wx.ALIGN_LEFT)
        self.waveobs = wx.TextCtrl(self, -1, str(waveobs),  style=wx.TE_LEFT, size = wx.Size(500, -1))

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_waveobs, 0, border=3, flag=flags)
        self.hbox.Add(self.waveobs, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        #### Flux operate
        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.operate_flux = wx.CheckBox(self, -1, 'Modify fluxes', style=wx.ALIGN_LEFT)
        #self.operate_flux.SetValue(True)

        #self.hbox.AddSpacer(10)
        #self.hbox.Add(self.operate_flux, 0, border=3, flag=flags)

        #self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Value
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_flux = wx.StaticText(self, -1, "Fluxes\t = ", style=wx.ALIGN_LEFT)
        self.flux = wx.TextCtrl(self, -1, str(flux),  style=wx.TE_LEFT, size = wx.Size(500, -1))

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_flux, 0, border=3, flag=flags)
        self.hbox.Add(self.flux, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        #### Err operate
        #self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        #self.operate_err = wx.CheckBox(self, -1, 'Modify errors', style=wx.ALIGN_LEFT)
        #self.operate_err.SetValue(True)

        #self.hbox.AddSpacer(10)
        #self.hbox.Add(self.operate_err, 0, border=3, flag=flags)

        #self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Value
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_err = wx.StaticText(self, -1, "Errors\t = ", style=wx.ALIGN_LEFT)
        self.err = wx.TextCtrl(self, -1, str(err),  style=wx.TE_LEFT, size = wx.Size(500, -1))

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_err, 0, border=3, flag=flags)
        self.hbox.Add(self.err, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(20)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)


    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)


class SendSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, applications):
        wx.Dialog.__init__(self, parent, id, title, size=(450,150))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        ### Application
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_application = wx.StaticText(self, -1, "Application: ", style=wx.ALIGN_LEFT)
        self.application = wx.ComboBox (self, wx.ID_ANY, applications[0], choices=applications, style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_application, 0, border=3, flag=flags)
        self.hbox.Add(self.application, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(30)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class SyntheticSpectrumDialog(wx.Dialog):
    def __init__(self, parent, id, title, wave_base, wave_top, wave_step, resolution, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff):
        wx.Dialog.__init__(self, parent, id, title, size=(450,650))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        ### Model atmosphere
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_atmospheres = wx.StaticText(self, -1, "Model atmosphere: ", style=wx.ALIGN_LEFT)
        self.atmospheres = wx.ComboBox (self, wx.ID_ANY, 'MARCS', choices=['Kurucz', 'Castelli', 'MARCS'], style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_atmospheres, 0, border=3, flag=flags)
        self.hbox.Add(self.atmospheres, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Linelist
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_linelist = wx.StaticText(self, -1, "Line list: ", style=wx.ALIGN_LEFT)
        self.linelist = wx.ComboBox (self, wx.ID_ANY, 'VALD', choices=['SPECTRUM', 'VALD'], style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_linelist, 0, border=3, flag=flags)
        self.hbox.Add(self.linelist, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

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

        ### Microturbulence
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_microturbulence_vel = wx.StaticText(self, -1, "Microturbulence velocity (km/s): ", style=wx.ALIGN_LEFT)
        self.microturbulence_vel = wx.TextCtrl(self, -1, str(microturbulence_vel),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_microturbulence_vel, 0, border=3, flag=flags)
        self.hbox.Add(self.microturbulence_vel, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Macro
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_macroturbulence = wx.StaticText(self, -1, "Macroturbulence velocity (km/s): ", style=wx.ALIGN_LEFT)
        self.macroturbulence = wx.TextCtrl(self, -1, str(macroturbulence),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_macroturbulence, 0, border=3, flag=flags)
        self.hbox.Add(self.macroturbulence, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### vsini
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_vsini = wx.StaticText(self, -1, "Rotation (v sin(i)) (km/s): ", style=wx.ALIGN_LEFT)
        self.vsini = wx.TextCtrl(self, -1, str(vsini),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_vsini, 0, border=3, flag=flags)
        self.hbox.Add(self.vsini, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Limb
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_limb_darkening_coeff = wx.StaticText(self, -1, "Limb darkening coeff: ", style=wx.ALIGN_LEFT)
        self.limb_darkening_coeff = wx.TextCtrl(self, -1, str(limb_darkening_coeff),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_limb_darkening_coeff, 0, border=3, flag=flags)
        self.hbox.Add(self.limb_darkening_coeff, 0, border=3, flag=flags)

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

        ### Wave step
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_wave_step = wx.StaticText(self, -1, "Wavelength step (nm): ", style=wx.ALIGN_LEFT)
        self.wave_step = wx.TextCtrl(self, -1, str(wave_step),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_wave_step, 0, border=3, flag=flags)
        self.hbox.Add(self.wave_step, 0, border=3, flag=flags)

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
        self.EndModal(wx.ID_OK)

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
        self.text_from_resolution_explanation = wx.StaticText(self, -1, "If initial resolution is set to zero, spectrum will be just smoothed by a gaussian determined by the final resolution (FWHM = wavelengths / final resolution)", style=wx.ALIGN_CENTER)
        self.vbox.Add(self.text_from_resolution_explanation, 1,  wx.CENTER | wx.TOP | wx.GROW)

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
        self.EndModal(wx.ID_OK)


class EstimateErrorsDialog(wx.Dialog):
    def __init__(self, parent, id, title, snr):
        wx.Dialog.__init__(self, parent, id, title, size=(250,200))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL


        ### SNR
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_snr = wx.StaticText(self, -1, "SNR: ", style=wx.ALIGN_LEFT)
        self.snr = wx.TextCtrl(self, -1, str(snr),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_snr, 0, border=3, flag=flags)
        self.hbox.Add(self.snr, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(10)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)


class FitLinesDialog(wx.Dialog):
    def __init__(self, parent, id, title, vel_atomic, vel_telluric):
        wx.Dialog.__init__(self, parent, id, title, size=(450,250))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Velocity respect to atomic data
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_vel_atomic = wx.StaticText(self, -1, "Velocity respect to atomic lines (km/s): ", style=wx.ALIGN_LEFT)
        self.vel_atomic = wx.TextCtrl(self, -1, str(vel_atomic),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_vel_atomic, 0, border=3, flag=flags)
        self.hbox.Add(self.vel_atomic, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        ### Velocity respect to telluric data
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.text_vel_telluric = wx.StaticText(self, -1, "Velocity respect to telluric lines (km/s): ", style=wx.ALIGN_LEFT)
        self.vel_telluric = wx.TextCtrl(self, -1, str(vel_telluric),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_vel_telluric, 0, border=3, flag=flags)
        self.hbox.Add(self.vel_telluric, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        self.vbox.AddSpacer(10)

        sizer =  self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class DetermineAbundancesDialog(wx.Dialog):
    def __init__(self, parent, id, title, teff, logg, MH, microturbulence_vel):
        wx.Dialog.__init__(self, parent, id, title, size=(450,450))

        self.action_accepted = False

        self.vbox = wx.BoxSizer(wx.VERTICAL)

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        ### Model atmosphere
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_atmospheres = wx.StaticText(self, -1, "Model atmosphere: ", style=wx.ALIGN_LEFT)
        self.atmospheres = wx.ComboBox (self, wx.ID_ANY, 'MARCS', choices=['Kurucz', 'Castelli', 'MARCS'], style=wx.CB_READONLY)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_atmospheres, 0, border=3, flag=flags)
        self.hbox.Add(self.atmospheres, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

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

        ### Microturbulence
        self.hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.text_microturbulence_vel = wx.StaticText(self, -1, "Microturbulence velocity (km/s): ", style=wx.ALIGN_LEFT)
        self.microturbulence_vel = wx.TextCtrl(self, -1, str(microturbulence_vel),  style=wx.TE_RIGHT)

        self.hbox.AddSpacer(10)
        self.hbox.Add(self.text_microturbulence_vel, 0, border=3, flag=flags)
        self.hbox.Add(self.microturbulence_vel, 0, border=3, flag=flags)

        self.vbox.Add(self.hbox, 1,  wx.LEFT | wx.TOP | wx.GROW)

        sizer = self.CreateButtonSizer(wx.CANCEL|wx.OK)
        self.vbox.Add(sizer, 0, wx.ALIGN_CENTER)
        self.vbox.AddSpacer(10)
        self.SetSizer(self.vbox)
        self.Bind(wx.EVT_BUTTON, self.on_ok, id=wx.ID_OK)

    def on_ok(self, event):
        self.action_accepted = True
        self.EndModal(wx.ID_OK)

class AbundancesDialog(wx.Dialog):
    def __init__(self, parent, id, title, linemasks, abundances, teff, logg, MH, microturbulence_vel):
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
        # Avoid using special notation that are not easy to understand in axis for big zoom
        myyfmt = ScalarFormatter(useOffset=False)
    	self.axes.get_xaxis().set_major_formatter(myyfmt)
    	self.axes.get_yaxis().set_major_formatter(myyfmt)

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

        ## Stats
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, "Temperature effective")
        self.stats.SetStringItem(num_items, 1, str(np.round(teff, 1)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Gravity")
        self.stats.SetStringItem(num_items, 1, str(np.round(logg, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Metallicity")
        self.stats.SetStringItem(num_items, 1, str(np.round(MH, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Microturbulence velocity (km/s)")
        self.stats.SetStringItem(num_items, 1, str(np.round(microturbulence_vel, 2)))
        num_items += 1
        self.stats.InsertStringItem(num_items, "Total number of lines")
        self.stats.SetStringItem(num_items, 1, str(len(abundances)))
        num_items += 1

        ## Draw
        elements = np.unique(linemasks['element'])
        for element in elements:
            flines = linemasks['element'] == element
            element_linemasks = linemasks[flines]
            element_abundances = abundances[flines]
            self.axes.plot(element_linemasks['VALD_wave_peak'], element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)
            self.stats.InsertStringItem(num_items, element + " abundance (%)")
            self.stats.SetStringItem(num_items, 1, str(np.round(100.*np.power(10, np.median(element_abundances)), 5)))
            num_items += 1
            self.stats.InsertStringItem(num_items, element + " median abundance in log(N/Ntot) (dex)")
            self.stats.SetStringItem(num_items, 1, str(np.round(np.median(element_abundances), 2)))
            num_items += 1
            self.stats.InsertStringItem(num_items, element + " mean abundance in log(N/Ntot) (dex)")
            self.stats.SetStringItem(num_items, 1, str(np.round(np.mean(element_abundances), 2)))
            num_items += 1
            self.stats.InsertStringItem(num_items, element + " standard deviation in abundance (dex)")
            self.stats.SetStringItem(num_items, 1, str(np.round(np.std(element_abundances), 2)))
            num_items += 1
            self.stats.InsertStringItem(num_items, element + " lines number")
            self.stats.SetStringItem(num_items, 1, str(np.round(len(element_abundances), 0)))
            num_items += 1

        leg = self.axes.legend(loc='upper right', shadow=False)
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='8')

        self.axes.grid(True, which="both")
        self.axes.set_title("Abundances", fontsize="10")
        self.axes.set_xlabel("wavelength (nm)", fontsize="10")
        self.axes.set_ylabel("abundance (dex)", fontsize="10")


    def on_no(self, event):
        self.recalculate = False
        self.EndModal(wx.ID_NO)

    def on_yes(self, event):
        self.recalculate = True
        self.EndModal(wx.ID_YES)
