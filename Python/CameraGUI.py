# -*- coding: utf-8 -*-
"""
@author: Ana Andres-Arroyo
GUI which controls a Thorlabs camera
"""
# create/update GUI file:
# pyuic4 CameraGUIdesign.ui -o CameraGUIdesign.py

# documentation:
# http://instrumental-lib.readthedocs.io/en/latest/uc480-cameras.html

import sys
#import CameraGUIdesign
from qtpy import QtGui, QtCore, QtWidgets, uic
import pyqtgraph as pg
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from MatplotlibSettings import *

from nplab.ui.ui_tools import UiTools

from instrumental import list_instruments, instrument

#CameraGUIdesign.Ui_CameraGUI,UiTools
class ClassCameraGUI(QtWidgets.QMainWindow,UiTools):
    """
    GUI which controls a Thorlabs camera.
    """
    
    def __init__(self):
        super(self.__class__, self).__init__()
        ui_file = 'C:\Users\Ana Andres\Documents\GitHub\optofluidics\Python\CameraGUIdesign.ui'
        uic.loadUi(ui_file, self)
     #   self.setupUi(self)
        
        # set starting parameters
        self.ExposureTimeNumberBox.setValue(5)
      
        # Connect GUI elements
        self.OpenCameraPushButton.clicked.connect(self.open_camera)
        self.CloseCameraPushButton.clicked.connect(self.close_camera)
        self.TakeImagePushButton.clicked.connect(self.take_image)
#        self.ExposureTimeNumberBox.valueChanged.connect(self.take_image)
        self.AutoExposureCheckBox.stateChanged.connect(self.set_auto_exposure)
#        self.LiveViewCheckBox.stateChanged.connect(self.live_view)
        self.StartLiveViewPushButton.clicked.connect(self.live_view)
        self.StopLiveViewPushButton.setEnabled(False)
        
#        self.figure_liveview = Figure()
#        self.canvas_liveview = FigureCanvas(self.figure_liveview)
#        self.toolbar_liveview = NavigationToolbar(self.canvas_liveview, self)
#        self.LiveViewWidgetContainer.addWidget(self.toolbar_liveview)
#        self.LiveViewWidgetContainer.addWidget(self.canvas_liveview)    
#        self.figure_liveview.clear()
#        self.ax = self.figure_liveview.add_subplot(111)
        
        # initialise empty dictionary
        self.capture_params = dict()
        self.video_params = dict()
                
        print list_instruments()
        self.open_camera()
        self.imv = pg.ImageView()
        self.replace_widget(self.verticalLayout,self.LiveViewWidget,self.imv)
     #   self.imv.layout()
        
#        self.imv.show()
        self.take_image()
        
        # initialise live view thread        
#        self.LiveView = LiveView(self.camera, self.LiveViewWidgetContainer)
#        self.LiveView = LiveView(self.camera)
    
    def open_camera(self):
        """Connect to a Thorlabs camera.""" 
        self.camera = instrument('uc480')
        print 'Camera connection started.'  
        self.camera.start_live_video()

        self.ROIWidthNumberBox.setValue(self.camera.max_width)
        self.ROIHeightNumberBox.setValue(self.camera.max_height)
        self.CameraWidthLabel.setText(str(self.camera.max_width))
        self.CameraHeightLabel.setText(str(self.camera.max_height))
        
        
        
    def close_camera(self):
        """Close the Thorlabs camera connection.""" 
        print 'close_camera not implemented yet.'
#        self.camera.close()
#        print 'Camera connection closed.'        
    
    def take_image(self):
        """Grab an image and display it."""
        self.qimage = self.grab_image()
        self.display_image(self.qimage)
        
    def display_image(self, image):
        """Display the latest captured image."""
     #   print image
        
        self.imv.setImage(image, autoLevels=True, autoHistogramRange=True)
     #   self.imv.show()
        
        
#        imv.setImage(pg.imageToArray(qimage))        
#        self.pixmap = QtGui.QPixmap.fromImage(qimage)
#        pixmap.detach()
#        pixmap_image = QtGui.QPixmap(pixmap)
#        self.label = QtGui.QLabel(self.LiveViewWidget)
#        label = QtGui.QLabel()
#        self.label.setPixmap(self.pixmap)
#        label.setAlignment(QtCore.Qt.AlignCenter)
#        label.setScaledContents(True)
#        label.setMinimumSize(1,1)
#        self.label.show()
#        label.app.processEvents()
#        self.LiveViewWidgetContainer.addWidget(label)        
#        if type(images) == tuple:
#            image = images[-1]
#        else:
#            image = images
#        self.ax.imshow(image)
#        self.canvas_liveview.draw()
            
    def set_capture_parameters(self):
        """Read capture parameters from the GUI."""
        self.capture_params['n_frames'] = int(self.FramesNumberBox.value())
#        self.capture_params['framerate'] = self.FrameRateNumberBox.value()
        self.capture_params = self.set_exposure_time(self.capture_params)
        self.capture_params = self.set_ROI(self.capture_params)
        print self.capture_params
    
    def set_ROI(self, params_dict):
        """Read ROI coordinates from the GUI."""
        # TODO: width and height must be multiples of 100: why?
        # TODO: it takes to captures to update the ROI: why?
        
        ROI_dict = {'width':[self.ROIWidthCheckBox, self.ROIWidthNumberBox],
                    'height':[self.ROIHeightCheckBox, self.ROIHeightNumberBox],
                    'left':[self.ROILeftEdgeCheckBox, self.ROILeftEdgeNumberBox],
                    'right':[self.ROIRightEdgeCheckBox, self.ROIRightEdgeNumberBox],
                    'top':[self.ROITopEdgeCheckBox, self.ROITopEdgeNumberBox],
                    'bottom':[self.ROIBottomEdgeCheckBox, self.ROIBottomEdgeNumberBox],
                    'cx':[self.ROICentreXCheckBox, self.ROICentreXNumberBox],
                    'cy':[self.ROICentreYCheckBox, self.ROICentreYNumberBox],
                    }
        
        for item in ROI_dict.keys():
            if item in params_dict.keys():
                del params_dict[item]

        if self.ROICheckBox.checkState():            
            for item in ROI_dict.keys():
                if ROI_dict[item][0].checkState():
                    params_dict[item] = int(ROI_dict[item][1].value())
        
        return params_dict
        
    def set_exposure_time(self, params_dict):
        """Read exposure time from the GUI."""
        if self.AutoExposureCheckBox.checkState():
            if 'exposure_time' in params_dict.keys():
                del params_dict['exposure_time']
#            print 'Exposure time = AUTO'          
        else:
            exposure_time = "{} millisecond".format(str(self.ExposureTimeNumberBox.value()))
            params_dict['exposure_time'] = exposure_time
#            print 'Exposure time = ' + str(exposure_time)
        return params_dict
            
    def set_auto_exposure(self):
        """Enable or disable the auto exposure shutter."""
        self.camera.set_auto_exposure(self.AutoExposureCheckBox.checkState())
#        if self.AutoExposureCheckBox.checkState():
#            print 'Auto exposure is ON'
#        else:
#            print 'Auto exposure is OFF'
        
    def grab_image(self):
        """Grab an image with the camera."""
        self.set_capture_parameters()
        image = self.camera.grab_image(**self.capture_params)
#        qimage = QtGui.QImage(image, image.shape[0], image.shape[1], QtGui.QImage.Format_RGB32)
        print 'Image(s) grabbed'

        self.CurrentWidthLabel.setText(str(self.camera.width))
        self.CurrentHeightLabel.setText(str(self.camera.height))
        self.MaxWidthLabel.setText(str(self.camera.max_width))
        self.MaxHeightLabel.setText(str(self.camera.max_height))
        self.MaxFrameRateLabel.setText(str(round(self.camera.framerate.magnitude,3)))
        self.CurrentExposureLabel.setText(str(round(self.camera._get_exposure().magnitude,3)))
        
        return image
        
    def set_video_parameters(self):
        """Read video parameters from the GUI."""
        self.video_params['framerate'] = "{} hertz".format(str(self.FrameRateNumberBox.value()))
        self.video_params['timeout'] = "{} millisecond".format(str(self.TimeoutNumberBox.value()))
        self.video_params = self.set_exposure_time(self.video_params)
        self.video_params = self.set_ROI(self.video_params)
        print self.video_params

    def live_view(self):
        """Start/stop the live view."""
#        print 'live_view not implemented yet.'
        self.set_video_parameters()       
        
#        if self.LiveViewCheckBox.checkState():
        self.LiveView = LiveViewThread(self.camera)
        self.LiveView.display_signal.connect(self.display_image)
#        self.connect(self.LiveView, QtCore.SIGNAL("display_image(QImage)"), self.display_image)
        self.LiveView.finished.connect(self.done)
    #    self.connect(self.LiveView, QtCore.SIGNAL("finished()"), self.done)
        self.StopLiveViewPushButton.setEnabled(True)
        self.StopLiveViewPushButton.clicked.connect(self.LiveView.terminate)
        self.StartLiveViewPushButton.setEnabled(False)
        
        self.LiveView.start()
        print "Started live view"
#        else:
#            self.LiveView.terminate()
#            print "Terminated live view"
    
    def done(self):
        print "terminated live view"
        self.StopLiveViewPushButton.setEnabled(False)
        self.StartLiveViewPushButton.setEnabled(True)
            
        
class LiveViewThread(QtCore.QThread):
    """
    Thread wich allows live view of the camera.
    """
    display_signal = QtCore.Signal(np.ndarray)
    
    def __init__(self, camera):
        QtCore.QThread.__init__(self)
       
        self.camera = camera        
        
        
    def __del__(self):
        self.wait()
        
    def run(self):
        """Live view."""
        # TODO: connect framerate from the GUI
        # TODO: set exposure time for the live video
        framerate = "{} hertz".format(str(20))
        self.camera.start_live_video(framerate=framerate)
        timeout = "{} millisecond".format(str(1000))

        i = 0
#        while  i < 100 :
        while not self.isFinished():
#            i += 1
#            print i
            if self.camera.wait_for_frame(timeout=timeout):
                image = self.camera.latest_frame()
#                qimage = QtGui.QImage(image, image.shape[0], image.shape[1], QtGui.QImage.Format_RGB32)
                self.display_signal.emit(image)
             #   QtCore.Signal()
        
#        
#def main():        #sys.argv
#    app = QtWidgets.QApplication([])
#    gui = ClassCameraGUI()
#    gui.show()
#    gui.activateWindow()
#    
#    app.exec_()
    
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    gui = ClassCameraGUI()
    gui.show()
    gui.activateWindow()
    
  #  app.exec_()
