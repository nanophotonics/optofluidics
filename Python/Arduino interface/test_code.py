#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 13:46:11 2017

@author: philippkoehler
"""

from qtpy import QtWidgets, uic
import pyqtgraph as pg
from nplab.ui.ui_tools import UiTools
import arduino_due_optofluidics as arduino



class arduinoGUI(QtWidgets.QMainWindow,UiTools):
    def __init__(self):    
        super(self.__class__, self).__init__()
        ui_file='arduino_GUI_design.ui'
        uic.loadUi(ui_file, self)
        self.StartMeasurementPushButton.clicked.connect(self.run_measurement)
        self.SinglePulseResponsePushButton.clicked.connect(self.run_SinglePulseResponse)
         #self.MeasuremenTypeComboBox.addItem('Delay Trace',0)
         #self.MeasuremenTypeComboBox.addItem('Intensity Trace',1)
         #self.MeasuremenTypeComboBox.setCurrentIndex(1)
        self.TimeUnitComboBox.addItem('min',0)
        self.TimeUnitComboBox.addItem('sec',1)
        self.TimeUnitComboBox.setCurrentIndex(0)
#        self.SelectNewFolderPushButton.clicked.connect(self.select_new_folder)
        self.LaserOnPushButton.clicked.connect(self.laser_on)
        self.LaserOffPushButton.clicked.connect(self.laser_off)
        self.ClearGUIPushButton.clicked.connect(self.clearGUI)
        #self.NumberOfPointsDoubleSpinBox.setValue(200)
        
        # open arduino        
        self.due = arduino.ArduinoDueOptofluidics("COM9")
        
        # show data browser
        self.due.data_file.show_gui(blocking=False)
        
#        plot_widget = pg.GraphicsLayoutWidget()
        self.replace_widget(self.PlotVerticalLayout, self.SignalPlotWidget, self.due.gui)
        
        #self.select_new_folder()
    
    def laser_on(self):
        #self.send("set_laser_on",True)
        self.due
        print 'Laser on'
        
    def laser_off(self):
        #self.send("set_laser_on",False)
        self.due.laser_on(False)
        print 'Laser off'
    
    def clearGUI(self):
        self.due.gui.clear()
        print 'GUI cleared'
        
    def run_SinglePulseResponse(self):
        print 'running single pulse response'

        parametersSS = dict()
#        
        NumberOfPoints=int(self.SSNumberOfPointsDoubleSpinBox.value())
        TimeIncrement=int(float(self.SSTimeIncrementDoubleSpinBox.value()))
        Averaging=int(self.SSAveragingDoubleSpinBox.value())
        
        parametersSS['num_points']=NumberOfPoints
        parametersSS['num_pulses']=Averaging
        parametersSS['time_increment_in_sec']=float(TimeIncrement)*1e-6
        print parametersSS
        
               
        concentration = self.ConcentrationLineEdit.text()
        sample = self.SampleNameLineEdit.text()

        # create attributes dictionary
        self.due.attributes = dict(parametersSS)
        self.due.attributes['sample'] = sample
        self.due.attributes['concentration_uM'] = concentration
        
        self.due.delay_trace(**parametersSS)

        
        
    def run_measurement(self):
        TypeOfMeasurement = 'Intensity Trace'
        TimeUnit = self.TimeUnitComboBox.currentText()
        
        print 'running ' + TypeOfMeasurement
        
        # buildin the parameter dictionary
        
        parameters = dict()
        
        frequency=int(self.FrequencyDoubleSpinBox.value())
                        
        if TimeUnit == 'min':
            time_sec=int(self.MeasurementTimeDoubleSpinBox.value())*60
        if TimeUnit == 'sec':
            time_sec=int(self.MeasurementTimeDoubleSpinBox.value())
        
        if int(frequency*float(self.TimeDelayDoubleSpinBox.value())) == 0:
            self.FrequencyDoubleSpinBox.setValue(int(1/float(self.TimeDelayDoubleSpinBox.value())))
            frequency=int(self.FrequencyDoubleSpinBox.value())
            
       
        parameters['num_points']=int(time_sec/float(self.TimeDelayDoubleSpinBox.value()))
            
        parameters['num_pulses']=int(frequency*float(self.TimeDelayDoubleSpinBox.value()))
        
        parameters['time_delay_in_sec']=self.TimeDelayDoubleSpinBox.value()
        parameters['points_per_pulse']=int(self.PointsPerPulseSpinBox.value())
        print parameters
        
               
        concentration = self.ConcentrationLineEdit.text()
        sample = self.SampleNameLineEdit.text()
#        print 'concentration = ' + concentration
#        print 'sample = ' + sample

        # create attributes dictionary
        self.due.attributes = dict(parameters)
        self.due.attributes['sample'] = sample
        self.due.attributes['concentration_uM'] = concentration
        
        
        # running the measurement
        self.due.pulse_intensity_trace(**parameters)

            
#    def select_new_folder(self):
#        print 'select new folder'
#        self.folder_path = self.FileNameLineEdit.text()
#        self.folder_path = QtWidgets.QFileDialog.getExistingDirectory(directory=self.folder_path)
#        self.FileNameLineEdit.clear()
#        self.FileNameLineEdit.insert(self.folder_path)
#        print self.folder_path
        
        
    def delay_trace(self, num_points, num_pulses, time_increment_in_sec):
        print 'running delay trace'
        print num_points
        print num_pulses
        print time_increment_in_sec

    
#    def pulse_intensity_trace(self, num_points, num_pulses, points_per_pulse, time_delay_in_sec):
#        print num_points
#        print num_pulses
#        print points_per_pulse
#        print time_delay_in_sec
#        self.due.pulse_intensity_trace(num_points, num_pulses, points_per_pulse, time_delay_in_sec)
        
        
if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    gui = arduinoGUI()
    gui.show()
    #due = arduino.ArduinoDueOptofluidics("COM5")
    gui.activateWindow()
    
