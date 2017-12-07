# -*- coding: utf-8 -*-
"""
@author: Ana Andres-Arroyo
GUI which controls a Thorlabs camera
"""
# create/update GUI file:
# pyuic4 CameraGUIdesign.ui -o CameraGUIdesign.py

# documentation:
# http://instrumental-lib.readthedocs.io/en/latest/uc480-cameras.html

from qtpy import QtCore, QtWidgets, uic
from scipy.misc import imsave
import pyqtgraph as pg
import numpy as np
from MatplotlibSettings import *
from nplab.ui.ui_tools import UiTools
from instrumental import list_instruments, instrument

class ClassCameraGUI(QtWidgets.QMainWindow,UiTools):
    """
    GUI which controls a Thorlabs camera.
    """
    
    def __init__(self):
        super(self.__class__, self).__init__()
        ui_file = 'C:\Users\Ana Andres\Documents\GitHub\optofluidics\Python\CameraGUIdesign.ui'
        uic.loadUi(ui_file, self)
        
        # set starting parameters
        self.file_path = ''
        self.ExposureTimeNumberBox.setValue(5)
      
        # Connect GUI elements
        self.TakeImagePushButton.clicked.connect(self.take_image)
        self.SaveImagePushButton.clicked.connect(self.save_image)
#        self.ExposureTimeNumberBox.valueChanged.connect(self.take_image)
        self.AutoExposureCheckBox.stateChanged.connect(self.set_auto_exposure)
        self.LiveViewCheckBox.stateChanged.connect(self.live_view)
        
        # initialise empty dictionary for capture and video parameters
        self.capture_params = dict()
        self.video_params = dict()
                
        print list_instruments()
        self.open_camera()
        
        # create live view widget
#        self.imv = pg.ImageView()        
#        self.imv.ui.roiBtn.hide()
#        self.imv.ui.menuBtn.hide()
#        self.imv.ui.histogram.hide()
#        self.replace_widget(self.verticalLayout, self.LiveViewWidget, self.imv)
        image_widget = pg.GraphicsLayoutWidget()
        self.replace_widget(self.verticalLayout, self.LiveViewWidget, image_widget)
        view_box = image_widget.addViewBox(row=1,col=1)        
        self.imv = pg.ImageItem()
        self.imv.setOpts(axisOrder='row-major')
        view_box.addItem(self.imv)
        view_box.setAspectLocked(True)

        self.take_image()

    
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
        self.camera.close()
        print 'Camera connection closed.'        
    
    def take_image(self):
        """Grab an image and display it.
        First check if live view is running and terminate it"""
        image = self.grab_image()
        self.display_image(image)
        
    def display_image(self, image):
        """Display the latest captured image."""
        # make a copy of the data so it can be accessed when saving an image
        self.image = image
        # set levels to [0,255] because otherwise it autoscales when plotting
        self.imv.setImage(image, autoDownsample=True, levels=[0,255])        

            
    def set_capture_parameters(self):
        """Read capture parameters from the GUI."""
        # TODO: fix errors that ocurr when n_frames > 1
        self.capture_params['n_frames'] = int(self.FramesNumberBox.value())
        self.capture_params = self.set_exposure_time(self.capture_params)
        self.capture_params = self.set_ROI(self.capture_params)
        print self.capture_params
    
    def set_ROI(self, params_dict):
        """Read ROI coordinates from the GUI."""
        # TODO: not all widths and heights are allowed: why?        
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
        print 'Image(s) grabbed'

        self.CurrentWidthLabel.setText(str(self.camera.width))
        self.CurrentHeightLabel.setText(str(self.camera.height))
        self.MaxWidthLabel.setText(str(self.camera.max_width))
        self.MaxHeightLabel.setText(str(self.camera.max_height))
        self.display_framerate(self.camera.framerate.magnitude)
        self.display_exposure_time(self.camera._get_exposure().magnitude)
        return image
        
    def save_image(self):
        """Save the latest image as a .png file."""
        # make a copy of the image so the saved image is the one that was on the 
        # screen when the button was pressed, not when the file name was chosen
        image = self.image
        # user input to choose file name
        self.file_path = QtWidgets.QFileDialog.getSaveFileName(self, 'Save image', self.file_path, "PNG files (*.png)")
        # save image
        imsave(self.file_path, np.flip(image, axis=0))
        print "Image saved: " + self.file_path
        
    def set_video_parameters(self):
        """Read video parameters from the GUI."""
        self.timeout = "{} millisecond".format(str(self.TimeoutNumberBox.value()))
        self.video_params['framerate'] = "{} hertz".format(str(self.FrameRateNumberBox.value()))
        self.video_params = self.set_exposure_time(self.video_params)
        self.video_params = self.set_ROI(self.video_params)
        print self.video_params
    
    def display_framerate(self, framerate):
        """Display the current framerate on the GUI."""
        self.MaxFrameRateLabel.setText(str(round(framerate,3)))
    
    def display_exposure_time(self, exposure_time):
        """Display the current exposure time on the gui."""
        self.CurrentExposureLabel.setText(str(round(exposure_time,3)))
        
    def live_view(self):
        """Start/stop the live view."""
        if self.LiveViewCheckBox.isChecked():
            
            # start thread
            self.LiveView = LiveViewThread(self.camera)
            
            # connect signals
            self.LiveView.display_signal.connect(self.display_image)
            self.LiveView.framerate_signal.connect(self.display_framerate)
            self.LiveView.exposure_time_signal.connect(self.display_exposure_time)
            
            # connect buttons and functions            
            self.LiveViewCheckBox.stateChanged.connect(self.LiveView.terminate)
            self.TakeImagePushButton.setEnabled(False)
            self.LiveView.finished.connect(self.done)
            
            print "Starting live view."
            self.set_video_parameters()
            self.LiveView.set_video_parameters(self.video_params, self.timeout)
            self.LiveView.start()
            
    
    def done(self):
        print "Terminated live view."
        self.LiveViewCheckBox.setChecked(False)
        self.TakeImagePushButton.setEnabled(True)
            
        
class LiveViewThread(QtCore.QThread):
    """
    Thread wich allows live view of the camera.
    """
    display_signal = QtCore.Signal(np.ndarray)
    framerate_signal = QtCore.Signal(float)
    exposure_time_signal = QtCore.Signal(float)
    
    def __init__(self, camera):
        QtCore.QThread.__init__(self)       
        self.camera = camera                
        
    def __del__(self):
        self.wait()
    
    def set_video_parameters(self, video_params, timeout):
        """Receive video parameters from main GUI."""
        self.timeout = timeout        
        self.video_params = video_params
        
    def run(self):
        """Start live view and continuously acquire frames."""
        self.camera.start_live_video(**self.video_params)
        self.framerate_signal.emit(self.camera.framerate.magnitude)
        self.exposure_time_signal.emit(self.camera._get_exposure().magnitude)

        while not self.isFinished():
            if self.camera.wait_for_frame(timeout=self.timeout):
                image = self.camera.latest_frame()
                self.display_signal.emit(image)

    
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    gui = ClassCameraGUI()
    gui.show()
    gui.activateWindow()
