# -*- coding: utf-8 -*-
"""
@author: Ana Andres-Arroyo
GUI which controls a Thorlabs camera
"""

# pyuic4 CameraGUIdesign.ui -o CameraGUIdesign.py

import sys
import CameraGUIdesign
from PyQt4 import QtGui

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from MatplotlibSettings import *

#import matplotlib.image as mpimg
#import matplotlib.pyplot as plt

from instrumental import list_instruments, instrument


class ClassCameraGUI(QtGui.QMainWindow, CameraGUIdesign.Ui_CameraGUI):
    """
    GUI which controls a Thorlabs camera
    """
    def __init__(self, parent=None):
        super(ClassCameraGUI, self).__init__(parent)
        self.setupUi(self)
        
        # set starting parameters
        self.ExposureTimeNumberBox.setValue(10)
      
        # Connect GUI elements
        self.OpenCameraPushButton.clicked.connect(self.open_camera)
        self.CloseCameraPushButton.clicked.connect(self.close_camera)
        self.TakeImagePushButton.clicked.connect(self.take_image)
#        self.ExposureTimeNumberBox.valueChanged.connect(self.take_image)
        self.AutoExposureCheckBox.stateChanged.connect(self.set_auto_exposure)
        self.LiveViewCheckBox.stateChanged.connect(self.live_view)
        
        self.figure_liveview = Figure()
        self.canvas_liveview = FigureCanvas(self.figure_liveview)
        self.toolbar_liveview = NavigationToolbar(self.canvas_liveview, self)
        self.LiveViewWidgetContainer.addWidget(self.toolbar_liveview)
        self.LiveViewWidgetContainer.addWidget(self.canvas_liveview)    
        self.figure_liveview.clear()
        self.ax = self.figure_liveview.add_subplot(111)
                
        print list_instruments()
        self.open_camera()
        self.take_image()
        
    
    def open_camera(self):
        """Connect to a Thorlabs camera.""" 
        self.camera = instrument('uc480')
        print 'Camera connection started.'     
#        print self.camera.color_mode
    
    def close_camera(self):
        """Close the Thorlabs camera connection.""" 
        self.camera.close()
        print 'Camera connection closed.'
    
    def take_image(self):
        """Grab an image and display it."""
        self.grab_image()
        self.display_image()
        
    def display_image(self):
        """Display an image."""
        self.ax.imshow(self.image)
        self.canvas_liveview.draw()

    def set_auto_exposure(self):
        """Enable or disable the auto exposure shutter."""
        self.camera.set_auto_exposure(self.AutoExposureCheckBox.checkState())
        if self.AutoExposureCheckBox.checkState():
            print 'Auto exposure is ON'
        else:
            print 'Auto exposure is OFF'
        
    def grab_image(self):
        """Grab an image with the camera
        If live view is on it will take the latest frame."""
        if self.AutoExposureCheckBox.checkState():
            print 'Exposure time = AUTO'
            self.image = self.camera.grab_image()            
        else:
            exposure_time = "{} millisecond".format(str(self.ExposureTimeNumberBox.value()))
            print 'Exposure time = ' + str(exposure_time)
            self.image = self.camera.grab_image(exposure_time=exposure_time)
        print 'Image grabbed'

    def live_view(self):
        """Start/stop the live view."""
        print self.LiveViewCheckBox.checkState()
#        if self.LiveViewCheckBox.checkState():
#            self.take_image()
#            live_view = self.LiveViewCheckBox.checkState()
#            if live_view:
#                self.live_view()
            
#            self.camera.start_live_video()
#            i = 0
#            live_view = True
#            while live_view and i<10:
#                i += 1
#                print i
#                self.take_image()
#                live_view = self.LiveViewCheckBox.checkState()
#                timeout = "{} millisecond".format(str(1000))
#                image_exists = self.camera.wait_for_frame(timeout=timeout)
#                print image_exists
#                if image_exists:
##                    self.image = self.camera.get_captured_image()
#                    self.image = self.camera.latest_frame()
#                    self.display_image()
#            self.camera.stop_live_video()
        
        
        
        
if __name__ == '__main__':
    
    app = QtGui.QApplication(sys.argv)
    gui = ClassCameraGUI()
    gui.show()
    app.exec_()