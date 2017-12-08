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
        self.capture_parameters = dict()
        self.video_parameters = dict()
                
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
        self.capture_parameters['n_frames'] = int(self.FramesNumberBox.value())
        self.capture_parameters = self.set_exposure_time(self.capture_parameters)
        self.capture_parameters = self.set_ROI(self.capture_parameters)
        print "Capture parameters:"
        print self.capture_parameters
    
    def set_ROI(self, parameters_dict):
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
            if item in parameters_dict.keys():
                del parameters_dict[item]

        if self.ROICheckBox.checkState():            
            for item in ROI_dict.keys():
                if ROI_dict[item][0].checkState():
                    parameters_dict[item] = int(ROI_dict[item][1].value())
        
        return parameters_dict
        
    def set_exposure_time(self, parameters_dict):
        """Read exposure time from the GUI."""
        if self.AutoExposureCheckBox.checkState():
            if 'exposure_time' in parameters_dict.keys():
                del parameters_dict['exposure_time']        
        else:
            exposure_time = "{} millisecond".format(str(self.ExposureTimeNumberBox.value()))
            parameters_dict['exposure_time'] = exposure_time
        return parameters_dict
            
    def set_auto_exposure(self):
        """Enable or disable the auto exposure shutter."""
        self.camera.set_auto_exposure(self.AutoExposureCheckBox.checkState())
        
    def grab_image(self):
        """Grab an image with the camera."""
        # set the desired capture parameters
        self.set_capture_parameters() # populate self.capture_parameters
        # grab the image
        image = self.camera.grab_image(**self.capture_parameters)
        print 'Image(s) grabbed.'
        # get camera parameters and display them on the GUI
        camera_parameters = self.get_camera_parameters()
        self.display_camera_parameters(camera_parameters)        
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
        self.video_parameters['framerate'] = "{} hertz".format(str(self.FrameRateNumberBox.value()))
        self.video_parameters = self.set_exposure_time(self.video_parameters)
        self.video_parameters = self.set_ROI(self.video_parameters)
        print "Video parameters:"
        print self.video_parameters

    def get_camera_parameters(self):
        """Get camera parameters."""
        camera_parameters = dict()
        camera_parameters['framerate'] = self.camera.framerate
        camera_parameters['exposure_time'] = self.camera._get_exposure()
        camera_parameters['width'] = self.camera.width
        camera_parameters['max_width'] = self.camera.max_width
        camera_parameters['height'] = self.camera.height
        camera_parameters['max_height'] = self.camera.max_height
        return camera_parameters
    
    def display_camera_parameters(self, camera_parameters):
        """Display the current camera parameters on the GUI."""
        self.CurrentFrameRateLabel.setText(str(round(camera_parameters['framerate'].magnitude,3)))
        self.CurrentExposureLabel.setText(str(round(camera_parameters['exposure_time'].magnitude,3)))
        self.CurrentWidthLabel.setText(str(camera_parameters['width']))
        self.CurrentHeightLabel.setText(str(camera_parameters['height']))
        self.MaxWidthLabel.setText(str(camera_parameters['max_width']))
        self.MaxHeightLabel.setText(str(camera_parameters['max_height']))        
        
    def live_view(self):
        """Start/stop the live view."""
        if self.LiveViewCheckBox.isChecked():
            
            # start thread
            self.LiveView = LiveViewThread(self.camera)
            
            # connect signals
            self.LiveView.display_signal.connect(self.display_image)
            self.LiveView.camera_parameters_signal.connect(self.display_camera_parameters)
            
            # connect buttons and functions            
            self.LiveViewCheckBox.stateChanged.connect(self.LiveView.terminate)
            self.TakeImagePushButton.setEnabled(False)
            self.LiveView.finished.connect(self.done)
                        
            self.set_video_parameters() # populate self.video_parameters and self.timeout
            print "Starting live view."
            self.LiveView.start_live_view(self.video_parameters, self.timeout)
            camera_parameters = self.get_camera_parameters()
            self.display_camera_parameters(camera_parameters)
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
    camera_parameters_signal = QtCore.Signal(dict)
    
    def __init__(self, camera):
        QtCore.QThread.__init__(self)       
        self.camera = camera                
        
    def __del__(self):
        self.wait()
    
    def set_video_parameters(self, video_parameters, timeout):
        """Set video parameters."""
        self.timeout = timeout        
        self.video_parameters = video_parameters
            
    def start_live_view(self, video_parameters, timeout):
        """Start live view with the video parameters received from the main GUI."""
        self.set_video_parameters(video_parameters, timeout)
        self.camera.start_live_video(**self.video_parameters)
        
    def run(self):
        """Continuously acquire frames."""
        while not self.isFinished():
            if self.camera.wait_for_frame(timeout=self.timeout):
                # capture the latest frame
                image = self.camera.latest_frame()
                # send image to the main GUI
                self.display_signal.emit(image)

    
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    gui = ClassCameraGUI()
    gui.show()
    gui.activateWindow()
