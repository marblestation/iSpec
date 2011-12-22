#!/usr/bin/env python
#################
# Run with ipython -pdb -c "%run interactive.py"
#################
#import ipdb
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
    min_width = 0.01 # nm
    
    def __init__(self, frame, axvspan):
        self.frame = frame
        self.axvspan = axvspan
        self.press = None
    
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
            # Do not modify if the region will become slimmer than...
            if xy[2,0] - event.xdata > self.min_width:
                xy[0,0] = event.xdata
                xy[1,0] = event.xdata
                xy[4,0] = event.xdata
                self.frame.status_message("Modifying region from %.4f" % xy[0,0] + " to " + "%.4f" % xy[2,0])
                self.axvspan.figure.canvas.draw()
        elif event.button == 3:
            # Do not modify if the region will become slimmer than...
            if  event.xdata - xy[0,0] > self.min_width:
                xy[2,0] = event.xdata
                xy[3,0] = event.xdata
                self.frame.status_message("Modifying region from %.4f" % xy[0,0] + " to " + "%.4f" % xy[2,0])
                self.axvspan.figure.canvas.draw()
    
    def on_press(self, event):
        # Validate that the click is on this axis/object
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        contains, attrd = self.axvspan.contains(event)
        if not contains: return
        
        # When regions overlap two or more can receive the click event, so
        # let's use a lock to allow the modification of one of them
        if self.frame.lock.acquire(False):
            if self.frame.remove_mode:
                ## Remove region
                self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
                self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
                self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
                self.axvspan.set_visible(False)
                self.frame.regions_changed()
                self.frame.flash_status_message("Removed region from " + "%.4f" % self.axvspan.get_xy()[0,0] + " to " + "%.4f" % self.axvspan.get_xy()[2,0])
                self.axvspan.figure.canvas.draw()
                self.press = event.button, event.x, event.xdata
                # Do not free the lock until the user releases the click
            else:
                ## Modify region
                self.axvspan.set_color('r')
                self.axvspan.set_alpha(0.5)
                self.press = event.button, event.x, event.xdata
                self.update_size(event)
                self.frame.regions_changed()
                # Do not free the lock until the user releases the click
            
    
    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        
        # If it is in remove mode, the lock should be released
        # NOTE: The selected region has already been removed, so it will not  
        # receive the "release" event.
        if self.frame.remove_mode and self.frame.lock.locked():
            self.frame.lock.release()
        elif self.press != None:
            # In modification mode, if it is the current selected widget
            self.press = None
            self.axvspan.set_color('g')
            self.axvspan.set_alpha(0.5)
            self.frame.status_message("")
            self.frame.lock.release()
            self.axvspan.figure.canvas.draw()
        
    
    def on_motion(self, event):
        # Validate that the object has been pressed and the click is on this axis
        if self.press is None: return
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes.get_navigate_mode() != None: return
        
#        button, x, xpress = self.press
        self.update_size(event)
        


class SpectraFrame(wx.Frame):
    title = 'Spectra: Regions definition'
    
    def __init__(self, spectra, regions, regions_file, wave_base, wave_top):
        self.spectra = spectra
        self.regions = regions
        self.regions_file = regions_file
        self.remove_mode = False
        self.lock = threading.Lock()
        self.changes_not_saved = False
        
        if wave_base == None:
            self.wave_base = spectra['waveobs'].min()
        else:
            self.wave_base = wave_base
        if wave_top == None:
            self.wave_top = spectra['waveobs'].max()
        else:
            self.wave_top = wave_top
        
        wave_filter = (self.spectra['waveobs'] >= self.wave_base) & (self.spectra['waveobs'] < self.wave_top)    
        self.spectra_visible = self.spectra[wave_filter]
        
        wave_filter = (self.regions['wave_base'] >= self.wave_base) & (self.regions['wave_base'] < self.wave_top) & (self.regions['wave_top'] >= self.wave_base) & (self.regions['wave_top'] < self.wave_top)
        self.regions_visible = self.regions[wave_filter]
        
        wx.Frame.__init__(self, None, -1, self.title)
        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.draw_figure()
    
    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        m_expt = menu_file.Append(-1, "&Save regions\tCtrl-S", "Save regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_regions, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
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
        
        self.canvas.mpl_connect('button_press_event', self.on_press)
                
        text_mode = wx.StaticText(self.panel, -1, "Mode: ", style=wx.ALIGN_LEFT)
        self.radio_button_modify = wx.RadioButton(self.panel, -1, 'Modify', style=wx.RB_GROUP)
        self.radio_button_remove = wx.RadioButton(self.panel, -1, 'Remove')
        self.radio_button_modify.SetValue(True)
        self.Bind(wx.EVT_RADIOBUTTON, self.on_mode_change, id=self.radio_button_modify.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_mode_change, id=self.radio_button_remove.GetId())
        
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
        
        self.hbox.Add(text_mode, 0, border=3, flag=flags)
        
        self.hbox.Add(self.radio_button_modify, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_remove, 0, border=3, flag=flags)
        #self.hbox.AddSpacer(30)
        
        self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
    
   
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_figure(self):
        """ Redraws the figure
        """
        grid = True
        title = "Spectra"
        xlabel = "wavelength"
        ylabel = "flux"
        xdata = self.spectra_visible['waveobs']
        ydata = self.spectra_visible['flux']

        if grid:
            self.axes.grid(True, which="both")

        if title != None:
            self.axes.set_title(title, fontsize="10")

        if xlabel != None:
            self.axes.set_xlabel(xlabel, fontsize="10")

        if ylabel != None:
            self.axes.set_ylabel(ylabel, fontsize="10")

        #major_tick = np.round((xdata.max() - xdata.min()) / 100, decimals=2)
        #minor_tick = major_tick / 2
        #ax1.xaxis.set_major_locator(MultipleLocator(major_tick))
        #ax1.xaxis.set_minor_locator(MultipleLocator(minor_tick))
        #ax1.xaxis.set_major_formatter(FormatStrFormatter("%1.2f"))

        colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
        linestyles = ['', 'steps', 'steps:', '-', '--', ':']

        self.axes.plot(xdata, ydata, lw=1, color=colors[0], linestyle=linestyles[3], marker='', markersize=1, markeredgewidth=0, markerfacecolor=colors[0])
        
        self.region_widgets = [] # It is important not to lose the reference to the regions or the garbage collector will destroy them
        for c in self.regions_visible:
            axvs = self.axes.axvspan(c['wave_base'], c['wave_top'], facecolor='g', alpha=0.5)
            region = CustomizableRegion(self, axvs)
            region.connect()
            self.region_widgets.append(region)
        
        self.canvas.draw()
    
    def regions_changed(self):
        self.SetTitle('Spectra: Regions definition (*not saved)')
        self.changes_not_saved == True
    
    def regions_saved(self):
        self.SetTitle('Spectra: Regions definition')
        self.changes_not_saved == True
    
    def on_close(self, event):
        if self.changes_not_saved == True:
            dlg = wx.MessageDialog(self, 'Are you sure you want to exit without saving the regions?', 'Changes not saved', wx.YES|wx.NO|wx.ICON_QUESTION)
            if dlg.ShowModal() == wx.ID_YES:
                self.Destroy()
        else:
            self.Destroy()

    def on_mode_change(self, event):
        if self.radio_button_modify.GetId() == event.GetId():
            self.remove_mode = False
        else:
            self.remove_mode = True
    
    def on_press(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        
        new_halfwidth = 0.05
        # If the middle button, create a new region
        if event.button == 2 and event.key == None:
            axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='g', alpha=0.5)
            region = CustomizableRegion(self, axvs)
            region.connect()
            self.region_widgets.append(region)
            self.canvas.draw()
            self.regions_changed()
            self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - new_halfwidth) + " to " + "%.4f" % (event.xdata + new_halfwidth))
    

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Save plot as...",
                defaultDir=os.getcwd(),
                defaultFile=self.regions_file + ".png",
                wildcard=file_choices,
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    dlg_confirm = wx.MessageDialog(self, 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename(), 'File already exists', wx.YES|wx.NO|wx.ICON_QUESTION)
                    if dlg_confirm.ShowModal() == wx.ID_NO:
                        continue # Give the oportunity to select a new file name
                self.canvas.print_figure(path, dpi=self.dpi)
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                action_ended = True
    
    def on_save_regions(self, event):
        file_choices = "All|*"
        
        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self, 
                message="Save regions as...",
                defaultDir=os.getcwd(),
                defaultFile=self.regions_file,
                wildcard=file_choices,
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    dlg_confirm = wx.MessageDialog(self, 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename(), 'File already exists', wx.YES|wx.NO|wx.ICON_QUESTION)
                    if dlg_confirm.ShowModal() == wx.ID_NO:
                        continue # Give the oportunity to select a new file name
                self.regions_file = dlg.GetFilename()
                output = open(path, "w")
                output.write("wave_base\twave_top\n")
                
                ## Regions not visible that are below the visible range
                wave_filter = (self.regions['wave_base'] < self.wave_base) & (self.regions['wave_top'] < self.wave_base)
                for region in self.regions[wave_filter]:
                    output.write("%.4f" % region['wave_base'] + "\t" + "%.4f" % region['wave_top'] + "\n")
                
                ## Visible regions: maybe they have been modified, removed or created
                for region in self.region_widgets:
                    # If the widget is not visible, it has been deleted by the user
                    if region.axvspan.get_visible():
                        output.write("%.4f" % region.axvspan.get_xy()[0,0] + "\t" + "%.4f" % region.axvspan.get_xy()[2,0] + "\n")
                
                ## Regions not visible that are above the visible range
                wave_filter = (self.regions['wave_base'] > self.wave_top) & (self.regions['wave_top'] > self.wave_top)
                for region in self.regions[wave_filter]:
                    output.write("%.4f" % region['wave_base'] + "\t" + "%.4f" % region['wave_top'] + "\n")

                output.close()
                self.regions_saved()
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                action_ended = True
        
    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ Spectra visualizer for interactive regions definition:
        
         * Modify mode:
            - Click on a space without regions to create a new one
            - Left/Right click on a region modifies its left/right limit
         * Remove mode:
            - Click on the region to be removed
         * Zoom mode:
            - Click and drag to define the zone to be zoomed
         * Pan mode:
            - Click and drag to move the plot
            
        Modify/remove modes only work when zoom/pan modes are not be active.
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
    print sys.argv[0], "[--wave_base=#] [--wave_top=#] spectra_file regions_file"

## Interpret arguments
def get_arguments():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "bt", ["wave_base=", "wave_top="])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    
    wave_base = None
    wave_top = None
    create_regions = False
    spectra_file = None
    regions_file = None
    try:
        for o, a in opts:
            if o in ("-b", "--wave_base"):
                wave_base = float(a)
            elif o in ("-t", "--wave_top"):
                wave_top = float(a)
            else:
                print "Argument", o, "not recognized!"
                usage()
                sys.exit(2)
    except ValueError, err:
        print "Wave values should be numbers!"
        usage()
        sys.exit(2)
              
    if (len(args) == 0):
        print "Missing spectra and regions files as arguments!"
        usage()
        sys.exit(2)
    elif not os.path.exists(args[0]):
        print "File", args[0], "does not exists!"
        sys.exit(2)
    elif len(args) >= 2 and not os.path.exists(args[1]):
        print "File", args[1], "does not exists!"
        sys.exit(2)
    
    spectra_file = args[0]
    if len(args) == 2:
        regions_file = args[1]
    else:
        regions_file = spectra_file + ".txt"
        create_regions = True
    
    return wave_base, wave_top, create_regions, spectra_file, regions_file


if __name__ == '__main__':
    wave_base, wave_top, create_regions, spectra_file, regions_file = get_arguments()    
    
    #~ spectra = asciitable.read(table=spectra_file, delimiter=' ', data_start=2, names=['waveobs', 'flux', 'err'])
    #~ spectra.sort(order='waveobs') # Make sure it is ordered by wavelength
    spectra = read_spectra(spectra_file)
    
    if not create_regions:
        regions = asciitable.read(table=regions_file, delimiter='\t')
    else:
        regions = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
    
    app = wx.PySimpleApp()
    app.frame = SpectraFrame(spectra, regions, regions_file, wave_base, wave_top)
    app.frame.Show()
    app.MainLoop()

