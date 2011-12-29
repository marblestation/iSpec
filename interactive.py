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

# The recommended way to use wx with mpl is with the WXAgg backend.
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from common import *


class CustomizableRegion:
    min_width = 0.002 # nm
    mark = None
    
    def __init__(self, frame, element_type, axvspan, mark=None, note=None):
        self.frame = frame
        self.axvspan = axvspan
        self.element_type = element_type
        self.press = None
        self.mark = mark
        self.note = note
        self.original_facecolor = self.axvspan.get_facecolor()
        self.original_edgecolor = self.axvspan.get_edgecolor()
    
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
        if self.frame.action == "Create": return
        
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
            
    def disconnect_and_hide(self):
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
        self.axvspan.set_visible(False)
        if self.mark != None:
            self.mark.set_visible(False)
        if self.note != None:
            self.note.set_visible(False)
    
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
            # If it was a right click when elements is lines and subelements marks
            # => Change note
            if event.button == 3 and self.frame.elements == "lines" and self.frame.subelements == "marks":
                # Right button, modify note
                if self.frame.action == "Modify" or self.frame.action == "Create":
                    self.update_mark_note(event)
                else:
                    # Note removal is managed on_press
                    pass
            
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
        
#        button, x, xpress = self.press
        if self.frame.subelements == "marks":
            self.update_mark(event)
        else:
            self.update_size(event)
        


class SpectraFrame(wx.Frame):
    title = 'Spectra: Regions definition'
    
    # regions should be a dictionary with 'continuum', 'lines' and 'segments' keys
    # filenames should be a dictionary with 'spectra', 'continuum', 'lines' and 'segments' keys
    def __init__(self, spectra=None, regions=None, filenames=None):
        if spectra == None:
            self.spectra = np.zeros((0,), dtype=[('waveobs', '<f8'), ('flux', '<f8'), ('err', '<f8')])
        else:
            self.spectra = spectra
        self.spectra_id = None
        
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
            self.filenames['spectra'] = None
            self.filenames['continuum'] = None
            self.filenames['lines'] = None
            self.filenames['segments'] = None
        else:
            self.filenames = filenames
        
        self.region_widgets = {} # It is important not to lose the reference to the regions or the garbage 
        self.region_widgets['continuum'] = []
        self.region_widgets['lines'] = []
        self.region_widgets['segments'] = []
        
        self.action = "Modify"
        self.elements = "continuum"
        self.subelements = None
        self.lock = threading.Lock()
        self.not_saved = {}
        self.not_saved["continuum"] = False
        self.not_saved["lines"] = False
        self.not_saved["segments"] = False
        
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
        m_expt = menu_file.Append(-1, "Save &continuum regions as\tCtrl-C", "Save continuum regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_continuum_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save &line regions as\tCtrl-L", "Save line regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_line_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save s&egments as\tCtrl-E", "Save segments to file")
        self.Bind(wx.EVT_MENU, self.on_save_segments, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the visual editor")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
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
        self.fig = Figure((10.0, 6.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(1, 1, 1)
        
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
                
        text_action = wx.StaticText(self.panel, -1, "Action: ", style=wx.ALIGN_LEFT)
        self.radio_button_create = wx.RadioButton(self.panel, -1, 'Create', style=wx.RB_GROUP)
        self.radio_button_modify = wx.RadioButton(self.panel, -1, 'Modify')
        self.radio_button_remove = wx.RadioButton(self.panel, -1, 'Remove')
        self.radio_button_modify.SetValue(True)
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_create.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_modify.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_remove.GetId())
        
        text_element = wx.StaticText(self.panel, -1, "Elements: ", style=wx.ALIGN_LEFT)
        self.radio_button_continuum = wx.RadioButton(self.panel, -1, 'Continuum', style=wx.RB_GROUP)
        self.radio_button_lines = wx.RadioButton(self.panel, -1, 'Lines')
        self.radio_button_linemarks = wx.RadioButton(self.panel, -1, 'Line marks')
        self.radio_button_segments = wx.RadioButton(self.panel, -1, 'Segments')
        self.radio_button_continuum.SetValue(True)
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_continuum.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_lines.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_linemarks.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_segments.GetId())
        
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
        
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
    
   
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_spectra(self):
        # Remove current spectra plot if exists
        if self.spectra_id != None:
            self.axes.lines.remove(self.spectra_id[0])
        
        self.spectra_id = self.axes.plot(self.spectra['waveobs'], self.spectra['flux'], lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b')
    
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
            
            # Segments always in the background
            if elements == "segments":
                axvs.zorder = 0
            
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
        
        self.draw_spectra()
        self.draw_regions("continuum")
        self.draw_regions("lines")
        self.draw_regions("segments")
        self.canvas.draw()
    
    def update_title(self):
        title = "Spectra: Regions definition"
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"]:
            title += "("
            if self.not_saved["continuum"]:
                title += "*continuum"
            if self.not_saved["lines"]:
                title += "*lines"
            if self.not_saved["segments"]:
                title += "*segments"
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
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"]:
            dlg = wx.MessageDialog(self, 'Are you sure you want to exit without saving the regions?', 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                self.Destroy()
        else:
            self.Destroy()

    def on_action_change(self, event):
        if self.radio_button_create.GetId() == event.GetId():
            self.action = "Create"
        elif self.radio_button_modify.GetId() == event.GetId():
            self.action = "Modify"
        else:
            self.action = "Remove"
    
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
                region = CustomizableRegion(self, "continuum", axvs)
                region.connect()
                self.region_widgets['continuum'].append(region)
            elif self.elements == "lines":
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='yellow', alpha=0.30)
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
                # Segments always in the background
                axvs.zorder = 0
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
    
    def open_file(self, elements):
        file_choices = "All|*"
        
        if elements != "spectra" and self.not_saved[elements]:
            dlg = wx.MessageDialog(self, 'Are you sure you want to open a new %s file without saving the current regions?' % elements, 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() != wx.ID_YES:
                self.flash_status_message("Discarded.")
                return
        
        if self.filenames[elements] != None:
            filename = self.filenames[elements].split('/')[-1]
            filename_length = len(filename)
            dirname = self.filenames[elements][:-filename_length]
        else:
            filename = ""
            dirname = os.getcwd()
        
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
                        self.spectra = read_spectra(path)
                        self.draw_spectra()
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
                    self.flash_status_message("Opened file %s" % path, flash_len_ms=10000)
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
    
    ## Save to a file
    # elements can be "continuum", "lines" or "segments"
    def save_regions(self, elements):
        file_choices = "All|*"
        saved = False
        elements = elements.lower()
        
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
                self.filenames[elements] = dirname + filename
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
                            output.write("%.4f" % region.mark.get_xdata()[0] + "\t" + "%.4f" % region.axvspan.get_xy()[0,0] + "\t" + "%.4f" % region.axvspan.get_xy()[2,0] + "\t" + note + "\n")
                        else:
                            output.write("%.4f" % region.axvspan.get_xy()[0,0] + "\t" + "%.4f" % region.axvspan.get_xy()[2,0] + "\n")
                
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
    
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ Spectra visualizer for interactive regions definition:
        
         * Create action:
            - Left click to create a new region
            - It will create the kind of region selected in the elements options
         * Modify action:
            - Left/Right click on a region to modify its left/right limit
            - Middle click on a line region to modify its mark
            - It will modify only regions of the kind selected in the elements options
         * Remove action:
            - Click on the region to be removed
            - It will remove only regions of the kind selected in the elements options
         * Zoom mode:
            - Click and drag to define the zone to be zoomed
         * Pan mode:
            - Click and drag to move the plot
            
        Create, modify and remove action are only effective when zoom/pan modes are not be active.
        
        The user is responsible for creating overlapping regions.
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
    
    def flash_status_message(self, msg, flash_len_ms=3000):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
        
    def status_message(self, msg):
        self.statusbar.SetStatusText(msg)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')


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
#~ python interactive.py --continuum=input/LUMBA/UVES_MRD_cmask.txt --lines=input/LUMBA/UVES_MRD_Fe-linelist.txt --segments=input/LUMBA/UVES_MRD_segments.txt input/LUMBA/UVES_MRD_sun_official.s.gz

## Input files should be '\t' separated, comments with # and the first line is the header
if __name__ == '__main__':
    filenames = get_arguments()    
    
    #### Read files

    if filenames['spectra'] != None:
        try:
            spectra = read_spectra(filenames['spectra'])
        except Exception:
            print "Spectra file", args[0], "has an incompatible format!"
            sys.exit(2)
    else:
        spectra = np.zeros((0,), dtype=[('waveobs', '<f8'), ('flux', '<f8'), ('err', '<f8')])
    
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

