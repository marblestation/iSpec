#!/usr/bin/env python
#################
# Run with ipython -pdb -c "%run visualize.py"
#################
# notebook.py

import asciitable
import numpy as np
import wx
import wx.lib.sheet as sheet
import ipdb

class MySheet(sheet.CSheet):
    def __init__(self, parent, num_rows, num_cols):
        sheet.CSheet.__init__(self, parent)
        self.SetNumberRows(num_rows)
        self.SetNumberCols(num_cols)

class Notebook(wx.Frame):

    def __init__(self, parent, id, title, data):
        wx.Frame.__init__(self, parent, id, title, size=(600, 500))
        menubar = wx.MenuBar()
        file = wx.Menu()
        file.Append(101, 'Quit', '' )
        menubar.Append(file, '&File')
        self.SetMenuBar(menubar)

        wx.EVT_MENU(self, 101, self.OnQuit)

        nb = wx.Notebook(self, -1, style=wx.NB_BOTTOM)
        
        self.sheet1 = self.create_sheet(nb, data)
        
        nb.AddPage(self.sheet1, 'Sheet1')
        
        self.sheet1.SetFocus()
        self.StatusBar()
        self.Centre()
        self.Show()

    def create_sheet(self, nb, data):
        if type(data) == list:
            return self.create_sheet_from_array(nb, np.array(data))
        elif type(data) == np.ndarray and data.dtype.names == None:
            return self.create_sheet_from_array(nb, data)
        elif (type(data) == np.core.records.recarray) or (type(data) == np.ndarray and data.dtype.names != None):
            return self.create_sheet_from_record(nb, data)
        else:
            raise Exception("Type %s not suported" % type(data))
    
    # Numpy records
    def create_sheet_from_record(self, nb, data):
        num_rows = len(data)
        num_cols = len(data.dtype.names)
        
        sheet = MySheet(nb, num_rows, num_cols)
        
        # Fill the spreadsheet with data
        for row in np.arange(num_rows):
            for col in np.arange(num_cols):
                sheet.SetCellValue(row, col, str(data[row][col]))
        ipdb.set_trace() 
        # Column names
        for i in np.arange(num_cols):
            sheet.SetColLabelValue(i, data.dtype.names[i])
        
        return sheet
    
    # Numpy array
    def create_sheet_from_array(self, nb, data):
        num_rows = data.shape[0]
        # If there is not 2 columns, consider just one
        if len(data.shape) != 2:
            num_cols = 1
        else:
            num_cols = data.shape[1]
        
        sheet = MySheet(nb, num_rows, num_cols)
        
        # Fill the spreadsheet with data
        for row in np.arange(num_rows):
            for col in np.arange(num_cols):
                sheet.SetCellValue(row, col, str(data[row][col]))
        
        return sheet

    def StatusBar(self):
        self.statusbar = self.CreateStatusBar()

    def OnQuit(self, event):
        self.Close()

if __name__ == '__main__':
    regions_file = "input/L082N03_spec_norm/10oct08/sp2_Normal/hd4614_001.s.txt"
    regions = asciitable.read(table=regions_file, delimiter='\t')
#    regions = np.array([(1,2,3), (3,4,5)], dtype=np.dtype([('x', 'int32'), ('y', 'int32'), ('z', 'int32')]))
#    regions = np.array([[1,2,3], [3,4,5]])
#    regions = [[1,2,3], [3,4,5]]
   
#    app = wx.PySimpleApp()
    app = wx.App()
    n = Notebook(None, -1, 'visualize.py', regions)
    app.MainLoop()
    
