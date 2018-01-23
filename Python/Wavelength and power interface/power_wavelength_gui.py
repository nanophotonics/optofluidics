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



class wavelength_controller(QtWidgets.QMainWindow, UiTools):
    
    def __init__(self):
        super(self.__class__, self).__init__()
        ui_file = 'power_wavelength_gui_design.ui'
        uic.loadUi(ui_file, self)
        
        # connect gui buttons
        self.SetWavelengthPushButton.clicked.connect(self.set_laser_wavelength)
        self.ReadPowerPushButton.clicked.connect(self.read_power)
        self.SetWaveplateAnglePushButton.clicked.connect(self.set_waveplate_angle)
        
        # set up waveplate
        print "Connecting to waveplate..."
        self.waveplate = PyAPT.APTMotor(SerialNum=27500609)  
        
        # set up powermeter
        rm = visa.ResourceManager()
        inst = rm.open_resource('USB0::0x1313::0x8078::P0011774::INSTR',timeout=1)
        self.power_meter=ThorlabsPM100(inst=inst)
        self.power_meter.system.beeper.immediate()
        
        # set up wavemeter
        self.wavemeter = wlm.Wavemeter()
        self.wavemeter.active = 1
        time.sleep(1)
        
        # start camer        
        self.camera_gui = uc480.uc480()
        self.camera_gui.show()
        self.camera_gui.activateWindow() 
        
        # read wavelength
        wavelength = self.read_laser_wavelength()
        self.WavelengthDoubleSpinBox.setValue(wavelength)
        # read power
        self.read_power()
        # read waveplate angle
        waveplate_angle = self.read_waveplate_angle()
        self.WaveplateAngleDoubleSpinBox.setValue(waveplate_angle)
        
        
    
    def set_laser_wavelength(self):
        self.wavelength = self.WavelengthDoubleSpinBox.value()
        self.read_laser_wavelength()
        
    def read_laser_wavelength(self):
        wavelength = self.wavemeter.wavelength # get wavelength from power meter
        self.LaserWavelengthLabel.setText(str(wavelength)) # update gui
        self.set_powermeter_wavelength(wavelength) # set powermeter_wavelength
        
        # TODO: update with latest camera gui
        self.camera_gui.CurrentExposureLabel.setText(str(wavelength)) 
        
        return wavelength        
        
    def read_power(self):
        power = self.power_meter.read
        self.PowerLabel.setText(str(power))
        return power
        
    def set_powermeter_wavelength(self, wavelength):
        self.power_meter.sense.correction.wavelength = wavelength # set powermeter wavelength
        powermeter_wavelength = self.power_meter.sense.correction.wavelength # read powermeter wavelength
        self.PowermeterWavelengthLabel.setText(str(powermeter_wavelength)) # update gui
        return powermeter_wavelength
        
    def set_waveplate_angle(self):
        waveplate_angle = self.WaveplateAngleDoubleSpinBox.value()
        self.waveplate.mAbs(waveplate_angle)
        self.read_waveplate_angle()
        
    def read_waveplate_angle(self):    
        # TODO: make button for this
        waveplate_angle = self.waveplate.getPos()
        self.WaveplateAngleLabel.setText(str(waveplate_angle))
        return waveplate_angle
        
    def closeEvent(self, event):
        print 'Closing GUIs'
        self.waveplate.cleanUpAPT()
        self.camera.close()
    
if __name__=='__main__':
    app = QtWidgets.QApplication([])
    c = wavelength_controller()
    c.show()
    c.activateWindow()
