#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 13:46:11 2017

@author: philippkoehler
"""

from qtpy import QtWidgets, uic
import pyqtgraph as pg
from nplab.ui.ui_tools import UiTools



class arduinoGUI(QtWidgets.QMainWindow,UiTools):
    def __init__(self):    
        super(self.__class__, self).__init__()
        ui_file='arduino_GUI_design.ui'
        uic.loadUi(ui_file, self)
        self.StartMeasurementPushButton.clicked.connect(self.run_measurement)
        self.MeasuremenTypeComboBox.addItem('Delay Trace',0)
        self.MeasuremenTypeComboBox.addItem('Intensity Trace',1)
        self.MeasuremenTypeComboBox.setCurrentIndex(1)
        self.TimeUnitComboBox.addItem('min',0)
        self.TimeUnitComboBox.addItem('sec',1)
        self.TimeUnitComboBox.setCurrentIndex(0)
        self.SelectNewFolderPushButton.clicked.connect(self.select_new_folder)
        #self.NumberOfPointsDoubleSpinBox.setValue(200)
        
        
        
        
        plot_widget = pg.GraphicsLayoutWidget()
        #self.replace_widget(self.VerticalLayout, self.SignalPlotWidget, plot_widget)
        
        #self.select_new_folder()
        
    
    def run_measurement(self):
        TypeOfMeasurement = self.MeasuremenTypeComboBox.currentText()
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
        
        

        if TypeOfMeasurement == 'Delay Trace':
            parameters['time_increment_in_sec']=1
        if TypeOfMeasurement == 'Intensity Trace':
            parameters['time_delay_in_sec']=self.TimeDelayDoubleSpinBox.value()
            parameters['points_per_pulse']=int(self.PointsPerPulseSpinBox.value())
        print parameters
        
               
        concentration = self.ConcentrationLineEdit.text()
        sample = self.SampleNameLineEdit.text()
        print 'concentration = ' + concentration
        print 'sample = ' + sample

        # running the measurement
        if TypeOfMeasurement == 'Delay Trace':
            self.delay_trace(**parameters)
        elif TypeOfMeasurement == 'Intensity Trace':
            self.pulse_intensity_trace(**parameters)

            
    def select_new_folder(self):
        print 'select new folder'
        self.folder_path = self.FileNameLineEdit.text()
        self.folder_path = QtWidgets.QFileDialog.getExistingDirectory(directory=self.folder_path)
        self.FileNameLineEdit.clear()
        self.FileNameLineEdit.insert(self.folder_path)
        print self.folder_path
        
        
    def delay_trace(self, num_points, num_pulses, time_increment_in_sec):
        print 'running delay trace'
        print num_points
        print num_pulses
        print time_increment_in_sec

    
    def pulse_intensity_trace(self, num_points, num_pulses, points_per_pulse, time_delay_in_sec):
        print num_points
        print num_pulses
        print points_per_pulse
        print time_delay_in_sec
        
        
if __name__ == '__main__':
    app =QtWidgets.QApplication([])
    gui = arduinoGUI()
    gui.show()
    gui.activateWindow()
    
