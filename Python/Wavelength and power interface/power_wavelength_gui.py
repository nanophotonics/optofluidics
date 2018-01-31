# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 13:59:49 2018

@author: Stijn de Graaf (sjd87)
"""

from qtpy import QtCore, QtWidgets, uic
from nplab.ui.ui_tools import UiTools
import visa
from ThorlabsPM100 import ThorlabsPM100
import wlm
import time
from nplab.instrument.stage import PyAPT
from nplab.instrument.camera import uc480
import pandas as pd
import numpy as np
import pyqtgraph as pg
import SolsTiS

class wavelength_controller(QtWidgets.QMainWindow, UiTools):
    
    def __init__(self):
        super(self.__class__, self).__init__()
        ui_file = 'power_wavelength_gui_design.ui'
        uic.loadUi(ui_file, self)
        
        # connect gui buttons
        self.SetWavelengthPushButton.clicked.connect(self.button_set_laser_wavelength)
        self.ReadWavelengthPushButton.clicked.connect(self.read_laser_wavelength)        
        self.ReadPowerPushButton.clicked.connect(self.read_power)
        self.SetWaveplateAnglePushButton.clicked.connect(self.button_set_waveplate_angle)
        self.ReadWaveplateAnglePushButton.clicked.connect(self.read_waveplate_angle)
        self.WavelengthSweepPushButton.clicked.connect(self.button_wavelength_sweep)
        self.WaveplateAngleSweepPushButton.clicked.connect(self.button_waveplate_sweep)
        

        # set gui parameters
        self.min_angle_step = 0.01
        self.min_wavelength_step = 0.01

        
        # set up waveplate
        print "Connecting to waveplate..."
        self.waveplate = PyAPT.APTMotor(SerialNum=27500609)  
        
        # set up powermeter
        rm = visa.ResourceManager()
        inst = rm.open_resource('USB0::0x1313::0x8078::P0011774::INSTR',timeout=1)
        self.power_meter=ThorlabsPM100(inst=inst)
        self.power_meter.system.beeper.immediate()
        
        # set up wavemeter
        self.wavemeter = wlm.Wavemeter(verbosity=False)
        self.wavemeter.active = 1
        time.sleep(1)
        
        # start camera        
        self.camera_gui = uc480.uc480()
        self.camera_gui.show()
        self.camera_gui.activateWindow() 
        
        # start laser        
        self.laser = SolsTiS.SolsTiS(('172.24.60.15', 39933))
        
        # read wavelength
        wavelength = self.read_laser_wavelength()
        self.WavelengthDoubleSpinBox.setValue(wavelength)
        # read power
        self.read_power()
        # read waveplate angle
        waveplate_angle = self.read_waveplate_angle()
        self.WaveplateAngleDoubleSpinBox.setValue(waveplate_angle)      
        
        
    def button_set_laser_wavelength(self):
        wavelength = self.WavelengthDoubleSpinBox.value()
        self.set_laser_wavelength(wavelength)
        self.read_power()
        
    def wait_for_wavelength(self, std=0.001, length=5, sleep=0.1):
        wavelength_std = std*10
        wavelength_list = []        
        while wavelength_std > std:
            time.sleep(sleep)
            wavelength = self.read_laser_wavelength()
            wavelength_list.append(wavelength)
            if len(wavelength_list) > length:
                wavelength_list.pop(0)
                wavelength_std = np.std(wavelength_list)  
        return wavelength
        
    def set_laser_wavelength(self, wavelength=800, precision=1, max_attempts = 5):
        # NOTE: if the wavemeter is linked to the laser the change_wavelength will fail
        self.laser.change_wavelength(wavelength)
        current_wavelength = self.wait_for_wavelength()
        offset = wavelength - current_wavelength
        attempt = 1
        while np.abs(offset) > precision and attempt < max_attempts:  
            attempt += 1            
            self.laser.change_wavelength(wavelength + offset)
            current_wavelength = self.wait_for_wavelength()
            offset = wavelength - current_wavelength
        print attempt
        return current_wavelength
        
    def read_laser_wavelength(self):
        wavelength = self.wavemeter.wavelength # get wavelength from power meter
        self.LaserWavelengthLabel.setText(str(wavelength)) # update gui
        self.set_powermeter_wavelength(wavelength) # set powermeter_wavelength        
        return wavelength        
        
    def read_power(self):
        # TODO: sometimes the value is of the order of e-37. figure out why
        power = self.power_meter.read
        self.PowerLabel.setText(str(power))
        return power
        
    def set_powermeter_wavelength(self, wavelength):
        self.power_meter.sense.correction.wavelength = wavelength # set powermeter wavelength
        powermeter_wavelength = self.power_meter.sense.correction.wavelength # read powermeter wavelength
        self.PowermeterWavelengthLabel.setText(str(powermeter_wavelength)) # update gui
        return powermeter_wavelength

    def button_set_waveplate_angle(self):
        waveplate_angle = self.WaveplateAngleDoubleSpinBox.value()
        self.set_waveplate_angle(waveplate_angle)
        self.read_power()

    def set_waveplate_angle(self, waveplate_angle):
        self.waveplate.mbAbs(waveplate_angle)
        waveplate_angle = self.read_waveplate_angle()
        return waveplate_angle
        
    def read_waveplate_angle(self):    
        waveplate_angle = self.waveplate.getPos()
        self.WaveplateAngleLabel.setText(str(waveplate_angle))
        return waveplate_angle

    def button_wavelength_sweep(self):
        start = self.WavelengthStartDoubleSpinBox.value()
        end = self.WavelengthEndDoubleSpinBox.value()
        step = self.WavelengthStepDoubleSpinBox.value()
        sample_description = self.SampleDescriptionLineEdit.text()
        self.wavelength_sweep(start, end, step, sample_description)
        
    def wavelength_sweep(self, start, end, step, sample_description):
        wavelength_range = np.arange(start, end + self.min_wavelength_step, step)
        sweep_dict = {'wavelength_nm':[], 'waveplate_angle_deg':[], 'power_w':[]}
        group_name = 'wavelength_scan_%d
        dg = self.camera_gui.df.require_group(group_name)
        group_name = dg.name
        for target_wavelength in wavelength_range:
            
            wavelength = self.set_laser_wavelength(target_wavelength) # set laser wavelength
            waveplate_angle = self.read_waveplate_angle() # read waveplate angle
            
            power = self.read_power() # read power from power meter
            
            # set camera attributes
            self.camera_gui.attributes['wavelength_nm'] = wavelength
            self.camera_gui.attributes['waveplate_angle_deg'] = waveplate_angle
            self.camera_gui.attributes['power_w'] = power
            self.camera_gui.attributes['sample_description'] = sample_description
            
            self.camera_gui.take_image() # take image
            self.camera_gui.save_image(group_name=group_name) # save image
            
            # write data to dictionary
            sweep_dict['wavelength_nm'].append(wavelength)
            sweep_dict['waveplate_angle_deg'].append(waveplate_angle)
            sweep_dict['power_w'].append(power)
            
        sweep_df = pd.DataFrame(data=sweep_dict)
        print sweep_df
        
        
    def button_waveplate_sweep(self):
        start = self.WaveplateStartDoubleSpinBox.value()
        end = self.WaveplateEndDoubleSpinBox.value()
        step = self.WaveplateStepDoubleSpinBox.value()
        sample_description = self.SampleDescriptionLineEdit.text()
        self.waveplate_sweep(start, end, step, sample_description)
        
    def waveplate_sweep(self, start, end, step, sample_description):
        # TODO: the current waveplate angle label is not being updated in the gui during the sweep: fix it
        waveplate_range = np.arange(start, end + self.min_angle_step, step) # range includes end value
        sweep_dict = {'wavelength_nm':[], 'waveplate_angle_deg':[], 'power_w':[]}
        group_name = 'waveplate_scan_%d
        dg = self.camera_gui.df.require_group(group_name)
        group_name = dg.name
        for target_angle in waveplate_range:
            
            waveplate_angle = self.set_waveplate_angle(target_angle) # set waveplate angle
            wavelength = self.read_laser_wavelength() # read wavelength
            power = self.read_power() # read power from power meter
            
            # set camera attributes
            self.camera_gui.attributes['wavelength_nm'] = wavelength
            self.camera_gui.attributes['waveplate_angle_deg'] = waveplate_angle
            self.camera_gui.attributes['power_w'] = power
            self.camera_gui.attributes['sample_description'] = sample_description
            
            self.camera_gui.take_image() # take image
            self.camera_gui.save_image(group_name=group_name) # save image
            
            # write data to dictionary
            sweep_dict['wavelength_nm'].append(wavelength)
            sweep_dict['waveplate_angle_deg'].append(waveplate_angle)
            sweep_dict['power_w'].append(power)
            
        sweep_df = pd.DataFrame(data=sweep_dict)
        print sweep_df
        
            
        
    def closeEvent(self, event):
        print "Closing GUIs..."
        self.waveplate.cleanUpAPT()
        self.camera_gui.close()
        # NOTE: should we close the wavemeter here or not?
        print "Done!"
    
if __name__=='__main__':
    app = QtWidgets.QApplication([])
    c = wavelength_controller()
    c.show()
    c.activateWindow()
