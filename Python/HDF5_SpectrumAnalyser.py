# -*- coding: utf-8 -*-
"""
Created on Mon May 09 15:26:26 2016

@author: Ana Andres
"""

#from PyQt4.uic import loadUiType

from nplab.utils.gui import *
from nplab.utils.notified_property import DumbNotifiedProperty
import nplab.ui.hdf5_browser as hdf5_browser
from PyQt4 import uic
from PyQt4 import QtGui, QtCore #TODO: I think these should be wrapped by nplab.utils.gui? rwb


import numpy as np
import h5py

import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Qt4Agg')
#from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.figure import Figure
#from nplab.ui.data_renderers import suitable_renderers
from nplab.ui.ui_tools import UiTools

import os


class HDF5SpectrumAnalyser(QtGui.QMainWindow, UiTools):
    """
    GUI which analyses the timelapse spectra recorded with Spectrometerself.timelapse.py    
    """
    def __init__(self, data_file, parent=None):
        super(HDF5SpectrumAnalyser, self).__init__(parent)
        self.data_file = data_file
        uic.loadUi(os.path.join(os.path.dirname(__file__), 'HDF5_SpectrumAnalyser_GUI.ui'), self)
        
#        self.tree = HDF5TreeWidget(self.data_file,
#                             parent=self,
#                             treeWidget=self.FileContentsTreeView, 
#                             refresh_button=self.RefreshTreePushButton,
#                             )
#        self.tree.treeWidget.itemClicked.connect(self.on_click)
                             
        self.treeWidget = hdf5_browser.HDF5TreeWidget(self.data_file,
                                                     parent=self,
                                                     )
        self.replace_widget(self.TreeWidgetContainer, self.FileContentsTreeView, self.treeWidget)        
                             
        self.treeWidget.selectionModel().selectionChanged.connect(self.selection_changed)
                             
#        self.viewer = HDF5ItemViewer(parent=self, 
#                                     show_controls=True,
#                                     )
               
        self.viewer = hdf5_browser.HDF5ItemViewer(parent=self, 
                                                 figure_widget=self.RendererFigureWidget,
                                                 show_controls=False, 
                                                 renderer_combobox = self.HDF5RendererComboBox,
                                                 default_button=self.DefaultRendererPushButton,
                                                 refresh_button=self.RefreshFigurePushButton,
                                                 copy_button=self.CopyFigurePushButton,
                                                 )             
        self.replace_widget(self.RendererWidgetContainer, self.RendererFigureWidget, self.viewer)
        
        self.NewFilePushButton.clicked.connect(self.new_file)
        self.RefreshTreePushButton.clicked.connect(self.refresh_tree)
        self.AnalysePlotPushButton.clicked.connect(self.analyse_and_plot_data)
        
        """Populate Default Combobox Options"""
        
#        self.BackgroundComboBox.addItem('None', 0)
#        self.ReferenceComboBox.addItem('None', 0)
        self.NormalisationComboBox.addItem('None', 0)
        
        self.TwoDPlotColourComboBox.addItem('Default Colours', 0)        
        self.SpectraColourComboBox.addItem('Default Colours', 0)        
        self.TraceColourComboBox.addItem('Default Colours', 0)        
        
        self.TwoDPlotFigureComboBox.addItem('New Figure', 0)
        self.SpectraFigureComboBox.addItem('New Figure', 0)
        self.TraceFigureComboBox.addItem('New Figure', 0)
        
        self.auto_connect_by_name()
    
    """"Connecting the checkboxes to the properties"""
    TwoDPlot = DumbNotifiedProperty()
    Spectra = DumbNotifiedProperty()
    Trace = DumbNotifiedProperty()
        
    
    def selection_changed(self, selected, deselected):
        """Callback function to update the displayed item when the tree selection changes."""
        self.viewer.data = self.treeWidget.selected_h5item()
        
#    def on_click(self, item, column):
#        """Handle clicks on items in the tree."""
##        item.setExpanded(True) # auto expand the item upon click
#        if len(self.treeWidget.selectedItems())>1: 
#            self.viewer.data = DummyHDF5Group({treeitem.data(column, QtCore.Qt.UserRole).name : treeitem.data(column, QtCore.Qt.UserRole) \
#                                                for treeitem in self.tree.treeWidget.selectedItems() })
#        else:
#            self.viewer.data = item.data(column, QtCore.Qt.UserRole)
        
#    def on_click(self, item, column):        
#        """Handle clicks on items in the tree."""
#        item.setExpanded(True)
#        if len(self.tree.treeWidget.selectedItems())>1: 
##            self.viewer.data = [treeitem.data(column, QtCore.Qt.UserRole) for treeitem in self.tree.treeWidget.selectedItems() ]
#            self.viewer.data = DummyHDF5Group({treeitem.data(column, QtCore.Qt.UserRole).name : treeitem.data(column, QtCore.Qt.UserRole) \
#                                                for treeitem in self.tree.treeWidget.selectedItems() })
#        else:
#            self.viewer.data = item.data(column, QtCore.Qt.UserRole)
            
        self.item_names = self.viewer.data.name
        print self.viewer.data.name

        self.timelapse = self.viewer.data
        self.Spectra = self.timelapse.items()

        raw = []
        self.spectra_name = []
        for Spectrum in self.Spectra:
            raw.append(np.array(self.timelapse.get(Spectrum[0])))
            self.spectra_name.append(Spectrum[0])
        self.intensity = raw            
        
        self.attributes = self.timelapse.get(self.Spectra[0][0]).attrs.items()
        
        """ Repopulate the Comboboxes"""
        # only checks the attributes of the last one: fix this!
        self.BackgroundComboBox.clear()
        self.BackgroundComboBox.addItem('None', 0)
        self.ReferenceComboBox.clear()
        self.ReferenceComboBox.addItem('None', 0)     
        self.BackgroundComboBox.addItem('From File', 3)
        self.ReferenceComboBox.addItem('From File', 3)
        for attribute in self.attributes:
            if attribute[0] == 'background':
                self.BackgroundComboBox.addItem('From Metadata', 3)
            if attribute[0] == 'reference':
                self.ReferenceComboBox.addItem('From Metadata', 3)                
                    
                        
        
        
    def new_file(self):
        """Select a new HDF5 file and populate the file contents tree widget."""
        new_file = nplab.datafile.open_file()
        assert new_file is not None, "Error opening the new file!  We'll stick with the old one :)"
        print "closing old file"
        self.data_file.close()
        self.data_file = new_file
        print "new file loaded:"
        print self.data_file.file
        print self.data_file.name
        print "changing tree widget's file"
        self.treeWidget.model.data_group = self.data_file
        print "refreshing tree"
        self.treeWidget.model.refresh_tree()
    
    def refresh_tree(self):
        """Refresh the file contents tree widget."""
        self.treeWidget.model.refresh_tree()
    
    def analyse_and_plot_data(self):
        
        """
        This is reading the attributes of the last spectra only.
        We assume that all spectra have the same attributes
        """
        
        for attribute in self.attributes:
        #    print attribute[0]
            if attribute[0] == 'wavelengths':
#                print 'Read Wavelengths'
                Wavelengths = attribute[1]
#            if attribute[0] == 'reference':
#                print 'Read Reference'
#                Reference = attribute[1]
#            if attribute[0] == 'background':
#                print 'Read Background'
#                Background = attribute[1]
            if attribute[0] == 'information':
#                print 'Read Information'
                Information = attribute[1]     
        
        plt.rcParams.update({'font.size': 32})
        
        plt.figure(figsize=[25,15])
        
        
        AxisSpectra = plt.subplot(111)
        PlotsSpectra = []
        for s in range(0,len(self.intensity),1):
            PlotsSpectra.append(
                               AxisSpectra.plot(
                                               Wavelengths,self.intensity[s], 
                                               '-', 
                                               label=self.spectra_name[s], 
                                               linewidth=4.0,
                                               )
                               )
                               
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity (a.u.)')
        plt.title(self.timelapse.name + ' // ' + Information)
        AxisSpectra.legend(numpoints = 1, loc='upper right')


if __name__ == '__main__':
    import sys, h5py, os, numpy as np
    import nplab
    
    print os.getcwd()
    app = get_qt_app()
#    data_file = h5py.File('C:/Users/Ana Andres/Documents/Python Scripts/2016-05-17.h5', 'r')
    data_file = nplab.datafile.open_file()
    ui = HDF5SpectrumAnalyser(data_file)
    ui.show()
