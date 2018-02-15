# -*- coding: utf-8 -*-
"""
@author: Ana Andres-Arroyo
GUI which controls a Thorlabs camera
"""
# documentation:
# http://instrumental-lib.readthedocs.io/en/latest/uc480-cameras.html

import os
import datetime
import time
from qtpy import QtCore, QtWidgets, uic
from scipy.misc import imsave
import pyqtgraph as pg
import numpy as np
import nplab
from nplab.ui.ui_tools import UiTools
from instrumental import list_instruments, instrument

def get_camera_parameters(camera):
    """Get camera parameters."""
    camera_parameters = dict()
    camera_parameters['framerate'] = camera.framerate.magnitude
    camera_parameters['exposure_time'] = camera._get_exposure().magnitude
    camera_parameters['width'] = camera.width
    camera_parameters['max_width'] = camera.max_width
    camera_parameters['height'] = camera.height
    camera_parameters['max_height'] = camera.max_height        
    camera_parameters['master_gain'] = camera.master_gain
    camera_parameters['gain_boost'] = camera.gain_boost
    camera_parameters['auto_blacklevel'] = camera.auto_blacklevel
    camera_parameters['blacklevel_offset'] = camera.blacklevel_offset
    try:
        camera_parameters['gamma'] = camera.gamma
    except:
#        print "WARNING: Can't read gamma value from the camera!!!"
        pass
    
    return camera_parameters

class uc480(QtWidgets.QMainWindow, UiTools):
    """
    GUI which controls a Thorlabs camera.
    """
    
    def __init__(self):
        super(self.__class__, self).__init__()
        # Get the path of this file in case we are calling this class from another location
        file_path = os.path.dirname(__file__)
        ui_file = file_path + '\uc480_gui_design.ui'
        uic.loadUi(ui_file, self)
        
        # Maximise GUI window
#        self.showMaximized()
        
        # Set initial tabs to display
        self.SettingsTabWidget.setCurrentIndex(0) 
        
        # Set initial splitter sizes
        self.splitter.setSizes([50,60000])
        
        # Disable stop video button
        self.StopVideoPushButton.setEnabled(False)     
      
        # Connect GUI elements
        self.AutoExposurePushButton.clicked.connect(self.auto_exposure)
        self.TakeImagePushButton.clicked.connect(self.take_image)
        self.SaveImagePushButton.clicked.connect(self.save_image)
        self.NewFilePushButton.clicked.connect(self.new_hdf5_file)
        self.LiveViewCheckBox.stateChanged.connect(self.start_live_view)
        self.StartVideoPushButton.clicked.connect(self.save_video)
        
        # Initialise empty dictionary for capture and video parameters
        self.capture_parameters = dict()
        self.video_parameters = dict()
        self.camera_parameters = dict()
        self.attributes = dict()
                        
        # Create live view widget
        image_widget = pg.GraphicsLayoutWidget()
        self.replace_widget(self.verticalLayout, self.LiveViewWidget, image_widget)
        view_box = image_widget.addViewBox(row=1,col=1)        
        self.imv = pg.ImageItem()
        self.imv.setOpts(axisOrder='row-major')
        view_box.addItem(self.imv)
        view_box.setAspectLocked(True)
        
        # Populate image format combobox
        self.ImageFormatComboBox.addItem('hdf5',0)
        self.ImageFormatComboBox.addItem('png',1)
        self.ImageFormatComboBox.addItem('tiff',2)
        self.ImageFormatComboBox.addItem('jpg',3)
        self.ImageFormatComboBox.setCurrentIndex(0)    
        
        # Populate video format combobox
        self.VideoFormatComboBox.addItem('hdf5',0)
        self.VideoFormatComboBox.setCurrentIndex(0)    

        # Open hdf5 file show hdf5 browser gui
        # TODO: only get hdf5 file the first time it's needed to save data
        self.df = False
        self.new_hdf5_file()

        # open camera
        print list_instruments()
        self.open_camera()
        
        # Set initial parameters
        self.file_path = ''
        self.ExposureTimeNumberBox.setValue(3)
        self.FramerateNumberBox.setValue(30)
        self.DisplayFramerateNumberBox.setValue(10)
        self.GainNumberBox.setValue(0)
        self.GammaNumberBox.setValue(1)
        self.BlacklevelNumberBox.setValue(255)       
#        self.ROICheckBox.setChecked(True)
#        self.ROIWidthCheckBox.setChecked(True)
#        self.ROIHeightCheckBox.setChecked(True)
#        self.ROIWidthNumberBox.setValue(700)
#        self.ROIHeightNumberBox.setValue(50)
        
        # take image and calculate the best exposure time
        self.auto_exposure()

    
    def open_camera(self):
        """Connect to a Thorlabs camera.""" 
        self.camera = instrument('uc480')
        print 'Camera connection started.'  

        self.ROIWidthNumberBox.setValue(self.camera.max_width)
        self.ROIHeightNumberBox.setValue(self.camera.max_height)
        self.CameraWidthLabel.setText(str(self.camera.max_width))
        self.CameraHeightLabel.setText(str(self.camera.max_height))                 

        
    def close_camera(self):
        """Close the Thorlabs camera connection.""" 
        self.camera.close()
        print 'Camera connection closed.'  
        
    def closeEvent(self,event):
        self.close_camera()
        self.df.close()
        self.df_gui.close()
    
    def take_image(self):
        """Grab an image and display it."""
        image = self.grab_image()
        self.display_image(image)
        return image

    def get_max_grayscale(self, image):
        max_grayscale = np.amax(image)
        self.CurrentMaxGrayLabel.setText(str(max_grayscale))
        return max_grayscale
        
    def auto_exposure(self):
        min_gray = self.MinGrayNumberBox.value()
        max_gray = self.MaxGrayNumberBox.value()
        precision = self.ExposureTimePrecisionNumberBox.value()
        self.set_auto_exposure(min_gray=min_gray, max_gray=max_gray, precision=precision)
    
    def set_auto_exposure(self, min_gray=200, max_gray=250, precision=1, max_attempts=10):
        image = self.take_image()
        max_image = self.get_max_grayscale(image)
        okay = True
        attempt = 0
        
        while (max_image > max_gray or max_image < min_gray) and okay:
            attempt += 1
            current_exposure = float(self.CurrentExposureLabel.text())
            
            if max_image > max_gray:
                print "REDUCE exposure time...\n"
                new_exposure = current_exposure/2
            elif max_image < min_gray:
                print "INCREASE exposure time...\n"
                new_exposure = current_exposure/max_image*max_gray*0.99
                
            self.ExposureTimeNumberBox.setValue(new_exposure)
            image = self.take_image()
            max_image = self.get_max_grayscale(image)            
                        
            previous_exposure = current_exposure
            current_exposure = float(self.CurrentExposureLabel.text())
            print "Exposure time = %f ms" %current_exposure
            print "Brightest pix = %d\n" %max_image
            if np.abs(previous_exposure - current_exposure) < precision:
                okay = False 
            if attempt > max_attempts:
                okay = False # make sure we're not trying forever

            # make sure the GUI is updated            
            QtWidgets.qApp.processEvents()

        
    def display_image(self, image):
        """Display the latest captured image."""
        # make a copy of the data so it can be accessed when saving an image
        self.image = image
        # set levels to [0,255] because otherwise it autoscales when plotting
        self.imv.setImage(image, autoDownsample=True, levels=[0,255])   
    
    def display_camera_parameters(self, camera_parameters):
        """Display the current camera parameters on the GUI."""
        self.CurrentFramerateLabel.setText(str(camera_parameters['framerate']))
        self.CurrentExposureLabel.setText(str(camera_parameters['exposure_time']))
        self.CurrentWidthLabel.setText(str(camera_parameters['width']))
        self.CurrentHeightLabel.setText(str(camera_parameters['height']))
        self.MaxWidthLabel.setText(str(camera_parameters['max_width']))
        self.MaxHeightLabel.setText(str(camera_parameters['max_height']))  
        self.CurrentMasterGainLabel.setText(str(camera_parameters['master_gain']))
        self.CurrentGainBoostLabel.setText(str(camera_parameters['gain_boost']))
        self.CurrentAutoBlacklevelLabel.setText(str(camera_parameters['auto_blacklevel']))
        self.CurrentBlacklevelLabel.setText(str(camera_parameters['blacklevel_offset']))
        try:
            self.CurrentGammaLabel.setText(str(camera_parameters['gamma']))
        except:
#            print "WARNING: Gamma not in camera parameters dictionary!!!"
            pass
        
            
    def set_capture_parameters(self):
        """Read capture parameters from the GUI."""
        self.capture_parameters = self.set_exposure_time(self.capture_parameters)
        self.capture_parameters = self.set_ROI(self.capture_parameters)
        self.capture_parameters['gain'] = float(self.GainNumberBox.value())
        self.capture_parameters['vbin'] = int(self.VBinNumberBox.value())
        self.capture_parameters['hbin'] = int(self.HBinNumberBox.value())
        self.capture_parameters['vsub'] = int(self.VSubNumberBox.value())
        self.capture_parameters['hsub'] = int(self.HSubNumberBox.value())
        self.set_camera_properties()        
    
    def set_camera_properties(self):
        """Read capture parameters from the GUI and set the corresponding camera properties."""
        self.camera.auto_blacklevel = self.AutoBlacklevelCheckBox.checkState()
        self.camera.blacklevel_offset = int(self.BlacklevelNumberBox.value())         
        self.camera.gain_boost = self.GainBoostCheckBox.checkState()
        try:
            self.camera.gamma = int(self.GammaNumberBox.value())
            # TODO: some camera models don't have a gamma. account for this in a neater way than try/except
        except:
#            print "WARNING: Can't set gamma!!!"
            pass
    
    def set_ROI(self, parameters_dict):
        """Read ROI coordinates from the GUI."""
        ROI_dict = {'width':[self.ROIWidthCheckBox, self.ROIWidthNumberBox],
                    'height':[self.ROIHeightCheckBox, self.ROIHeightNumberBox],
                    'left':[self.ROILeftEdgeCheckBox, self.ROILeftEdgeNumberBox],
                    'right':[self.ROIRightEdgeCheckBox, self.ROIRightEdgeNumberBox],
                    'top':[self.ROITopEdgeCheckBox, self.ROITopEdgeNumberBox],
                    'bottom':[self.ROIBottomEdgeCheckBox, self.ROIBottomEdgeNumberBox],
                    'cx':[self.ROICentreXCheckBox, self.ROICentreXNumberBox],
                    'cy':[self.ROICentreYCheckBox, self.ROICentreYNumberBox],
                    }
        
        # clear all of the old ROI parameters
        for item in ROI_dict.keys():
            if item in parameters_dict.keys():
                del parameters_dict[item]
        
        # use maximum width and height available
        parameters_dict['width'] = self.camera.max_width
        parameters_dict['height'] = self.camera.max_height
        
        # repopulate ROI parameters with the selected ones
        if self.ROICheckBox.checkState():            
            for item in ROI_dict.keys():
                if ROI_dict[item][0].checkState():
                    parameters_dict[item] = int(ROI_dict[item][1].value())
        
        return parameters_dict
        
    def set_exposure_time(self, parameters_dict):
        """Read exposure time from the GUI."""
        exposure_time = "{} millisecond".format(str(self.ExposureTimeNumberBox.value()))
        parameters_dict['exposure_time'] = exposure_time
        return parameters_dict            
        
    def grab_image(self):
        """Grab an image with the camera."""
        # set the desired capture parameters
        self.set_capture_parameters() # populate self.capture_parameters
        # get the capture_timestamp with millisecond precision
        # insert a T to match the creation_timestamp formatting
        self.attributes['capture_timestamp'] = str(datetime.datetime.now()).replace(' ', 'T')
        # grab the image
        image = self.camera.grab_image(**self.capture_parameters)
        print 'Image grabbed.\n'
        # get camera parameters and display them on the GUI
        self.camera_parameters = get_camera_parameters(self.camera)
        # update the attributes dictionary
        self.attributes.update(self.capture_parameters)
        self.attributes.update(self.camera_parameters)
        self.display_camera_parameters(self.camera_parameters)    
        self.get_max_grayscale(image)
        return image
        
    def get_info(self):
        """Get info from the GUI."""
        info = dict()
        info['description'] = self.DescriptionLineEdit.text()
        return info
        
    def save_image(self, group_name=False):
        """Save the latest image."""
        # make a copy of the image so the saved image is the one that was on the 
        # screen when the save button was pressed, not when the file name was chosen
        image = self.image
        image_format = self.ImageFormatComboBox.currentText()
        
        # when the save button is pressed group_name=False even if a default value is specified in the function definition
        # so we need this statement to make 'images' be the default
        if not group_name:
            group_name = 'images'
        
        if image_format == 'hdf5':
            # update the attributes dictionary
            self.attributes.update(self.get_info())
            # get the datafile
            df = nplab.current_datafile()
            # write data in the "images" group within the datafile
            dg = df.require_group(name=group_name)
            # write data to the file
            dg.create_dataset("image_%d", data=image, attrs=self.attributes)
            dg.file.flush()
            print "Image saved to the hdf5 file.\n"
            
        else:
            # user input to choose file name
            self.file_path = QtWidgets.QFileDialog.getSaveFileName(self, 'Save image', 
                                                                   self.file_path, 
                                                                   "(*."+self.ImageFormatComboBox.currentText()+")")
            if len(self.file_path):        
                # save image            
                imsave(self.file_path, np.flip(image, axis=0))
                print "Image saved: " + self.file_path + "\n"
            else:
                print "WARNING: Image wasn't saved.\n" 

    def new_hdf5_file(self):
        if self.df:
            # close the datafile
            self.df.close() 
            # open new datafile
            self.df = nplab.current_datafile() 
            # update the datafile for the gui
            self.df_gui.treeWidget.model.data_group = self.df 
            # refresh the tree
            self.df_gui.refresh_tree_button.click() 
        else:
            # open new datafile
            self.df = nplab.current_datafile()
            # open gui
            self.df_gui = self.df.show_gui(blocking=False)                    
        self.FilePathLineEdit.setText(self.df.filename)
        return self.df
        
    def set_video_parameters(self):
        """Read video parameters from the GUI."""
        self.timeout = "{} millisecond".format(str(self.TimeoutNumberBox.value()))
        self.video_parameters['framerate'] = "{} hertz".format(str(self.FramerateNumberBox.value()))
        self.video_parameters = self.set_exposure_time(self.video_parameters)
        self.video_parameters = self.set_ROI(self.video_parameters)
        self.video_parameters['gain'] = float(self.GainNumberBox.value())
        self.video_parameters['vbin'] = int(self.VBinNumberBox.value())
        self.video_parameters['hbin'] = int(self.HBinNumberBox.value())
        self.video_parameters['vsub'] = int(self.VSubNumberBox.value())
        self.video_parameters['hsub'] = int(self.HSubNumberBox.value())
        self.set_camera_properties()        
        
    def start_live_view(self):
        """Start/stop the live view."""
        if self.LiveViewCheckBox.isChecked():
                        
            # enable/disable gui buttons
            self.TakeImagePushButton.setEnabled(False)
            self.AutoExposurePushButton.setEnabled(False)
            self.StartVideoPushButton.setEnabled(False)
            self.StopVideoPushButton.setEnabled(False)
            
            # start thread
            self.LiveView = LiveViewThread(self.camera)
            # connect thread
            self.LiveViewCheckBox.stateChanged.connect(self.LiveView.terminate)
            self.LiveView.finished.connect(self.terminated_live_view)
            # live view
            print "Live view..."
            self.live_view(save=False)                       

            
    def save_video(self):
        """Save a video."""            
        # connect buttons and functions              
        self.LiveViewCheckBox.setEnabled(False)
        self.TakeImagePushButton.setEnabled(False)
        self.AutoExposurePushButton.setEnabled(False)
        self.StartVideoPushButton.setEnabled(False)
        self.StopVideoPushButton.setEnabled(True)
        
        # start thread
        self.LiveView = LiveViewThread(self.camera)
        # connect thread
        self.StopVideoPushButton.clicked.connect(self.LiveView.terminate)
        self.LiveView.finished.connect(self.terminated_live_view)
        # live view
        self.live_view(save=True)        
        
    def live_view(self, save=False):                
        # connect signals
        self.LiveView.display_signal.connect(self.display_image)
        self.LiveView.attributes_signal.connect(self.update_attributes)
        
        # populate self.video_parameters and self.timeout
        self.set_video_parameters() 
                
        # starting live view
        self.LiveView.start_live_view(self.video_parameters, 
                                      self.timeout, 
                                      save=save,
                                      total_frames=self.TotalFramesNumberBox.value(),
                                      display_framerate=self.DisplayFramerateNumberBox.value(),
                                      )
        
        self.LiveView.attributes.update(self.video_parameters)
        self.camera_parameters = get_camera_parameters(self.camera)    
        self.LiveView.attributes.update(self.camera_parameters)
        self.LiveView.attributes.update(self.get_info())
        self.display_camera_parameters(self.camera_parameters)
        self.LiveView.start()
    
    def update_attributes(self, attributes):
        self.attributes.update(attributes)   
        self.display_camera_parameters(self.attributes)                             
    
    def terminated_live_view(self):
        print "Terminated live view.\n"
        
        if self.LiveView.save:
            print "Saving video, please wait for a while..."
            
            # disable the GUI whilst the data is saved
            self.setEnabled(False)
            
            QtWidgets.qApp.processEvents()
            # get the datafile
            df = nplab.current_datafile()
            # get the "videos" group within the datafile
            datagroup = df.require_group("videos")
            # save video to the datafile
            # TODO: check wether this works with RGB camera too
            datagroup.create_dataset("video_%d", 
                                     data=self.LiveView.image_array,#.astype(int), 
                                     attrs=self.LiveView.attributes)
            
            del self.LiveView.image_array
            # flushing at the end
            datagroup.file.flush() 
            print "Done!\n"
            
            # enable the GUI again
            self.setEnabled(True)
            
            # remove capture_time_sec so it doesn't get recorded in still images
            del self.attributes['capture_time_sec']            
            del self.attributes['total_frames']    
        
        self.StopVideoPushButton.setEnabled(False)
        self.LiveViewCheckBox.setEnabled(True)
        self.LiveViewCheckBox.setChecked(False)
        self.TakeImagePushButton.setEnabled(True)
        self.AutoExposurePushButton.setEnabled(True)
        self.StartVideoPushButton.setEnabled(True)
        
            
        
class LiveViewThread(QtCore.QThread):
    """
    Thread wich allows live view of the camera.
    """
    display_signal = QtCore.Signal(np.ndarray)
    attributes_signal = QtCore.Signal(dict)
    
    def __init__(self, camera):
        QtCore.QThread.__init__(self)       
        self.camera = camera
        self.attributes = dict()               
        
    def __del__(self):
        self.wait()
    
    def set_video_parameters(self, video_parameters, timeout):
        """Set video parameters."""
        self.timeout = timeout        
        self.video_parameters = video_parameters
            
    def start_live_view(self, video_parameters, timeout, 
                        save=False, total_frames=2000,
                        display_framerate = 20):
        """Start live view with the video parameters received from the main GUI."""

        self.timeout = timeout        
        self.video_parameters = video_parameters
        self.save = save
        self.total_frames = total_frames
               
        if save:
            print "Recording video..."            
            self.capture_timestamp_array = np.empty(self.total_frames)
            self.image_array = np.empty([self.total_frames,
                                         self.video_parameters['height'],
                                         self.video_parameters['width']], 
                                         dtype='uint8')
        
        # get the capture_timestamp with millisecond precision
        # insert a T to match the creation_timestamp formatting
        self.attributes['capture_timestamp'] = str(datetime.datetime.now()).replace(' ', 'T')
        self.attributes['total_frames'] = total_frames
        
        # start timer with microsecond precision
        self.high_precision_time = HighPrecisionWallTime()     
        
        # start live video
        self.camera.start_live_video(**self.video_parameters)       
        
        # calculate when we need to emit the image to the gui        
        capture_framerate = self.camera.framerate.magnitude
        self.frame_multiple = int(capture_framerate / display_framerate)
        # if display_framerate > capture_framerate then frame_multiple < 1
        # since we cannot emit each image more than once, frame_multiple must be >= 1
        if self.frame_multiple < 1:
            self.frame_multiple = 1
        
        print "Framerate = %d" %capture_framerate
        print "Total frames = %d" %self.total_frames
        print "Frame multiple = %d\n" %self.frame_multiple
        
    def run(self):
        """Continuously acquire frames."""
        frame_number = 0
        while not self.isFinished() and frame_number < self.total_frames:
            if self.camera.wait_for_frame(timeout=self.timeout):
                # get the capture_time with microsecond precision
                capture_time_sec = self.high_precision_time.sample()                
                
                # capture the latest frame
                image = self.camera.latest_frame()
                if self.save:
                    self.save_frame(image, capture_time_sec, frame_number)
                
                if frame_number % self.frame_multiple == 0:
                    # emit signals to the main gui
                    self.attributes_signal.emit(self.attributes)
                    self.display_signal.emit(image) # crashes more often - maybe?               
                frame_number += 1
                
        
    def save_frame(self, image, capture_time_sec, frame_number):
        """Save the frame to memory."""
        self.capture_timestamp_array[frame_number] = capture_time_sec
        self.attributes['capture_time_sec'] = self.capture_timestamp_array        
        self.image_array[frame_number,:,:] = image
        

class HighPrecisionWallTime():
    def __init__(self,):
        self._wall_time_0 = time.time()
        self._clock_0 = time.clock()

    def sample(self,):
        dc = time.clock()-self._clock_0
        return self._wall_time_0 + dc
    
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    camera = uc480()
    camera.show()
    camera.activateWindow()

