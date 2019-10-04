# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:37:31 2016

@author: Ana Andres
"""
#modules
import nplab
import visa
import wlm
import matplotlib.pyplot as plt
import scipy
import os
import signal
import time
import threading
import numpy as np
import pandas as pd
import pyqtgraph as pg
from traits.api import HasTraits, Float, Int, String, Button, Bool
from traitsui.api import View, Group, HGroup, VGroup, Item
from nplab.utils.thread_utils import locked_action, background_action

#intruments
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
import nplab.instrument.spectrometer.seabreeze as seabreeze
#from ThorlabsPM100 import ThorlabsPM100
#from nplab.instrument.camera import uc480
from nplab.instrument.shutter import thorlabs_sc10

rm = visa.ResourceManager()

class SpectrometerTimelapse(nplab.instrument.Instrument, HasTraits):
    """An Instrument that saves spectra periodically."""
    number_of_spectra = Int(5)
    time_interval_seconds = Float(0.5)
    button_start = Button("Start")
    button_stop = Button("Stop")
    button_new_file = Button("New File")
    group_name = String("Experiment001")
    taking_spectra = Bool(False)
    group_spectra = Bool(False)
    info_string = String("Info...")

    # Set up shutter
    self.shutter = thorlabs_sc10.ThorLabsSC10()
    self.shutter.close_shutter()

    traits_view = View(VGroup(
                             HGroup(
                                 VGroup(
                                       Item("group_name", label = "Group Name"),
                                       Item("number_of_spectra", label = "Number of Spectra"),
                                       Item("time_interval_seconds", label = "Time Interval (s)"),
                                       Item("taking_spectra", label = "Taking Spectra"),
                                       ),
                                VGroup(
                                        Item("uv_start_delay", label = "UV Start Delay (# spectra)"),
                                        Item("uv_on_period", label = "UV On Period (# spectra)"),
                                        Item("uv_off_period", label = "UV Off Period (# spectra)"),
                                        ),
                                 VGroup(
                                       Item("button_start", show_label = False),
                                       Item("button_stop", show_label = False),
                                       Item("button_new_file", show_label = False),
                             ))),
                        Item("info_string", label = "Information", resizable=True, springy=False),
                        )
    traits_view.title = "Spectrometer Timelapse"
    traits_view.x = 200
    traits_view.y = 200
    traits_view.resizable = True

    def __init__(self,spectrometer):
        "Making a new timelapse GUI"
        super(SpectrometerTimelapse,self).__init__() # initialise the parent class
        self.spectrometer = spectrometer
        print "New timelapse GUI."

    def _button_start_fired(self):
        """This is called when "button_start" is clicked."""
        assert self.taking_spectra == False, "Error: already taking spectra"

        df = nplab.current_datafile()
        dg = df.require_group(self.group_name)

        datagroup = dg.create_group("timelapse_%d")
        self.acquisition_thread = threading.Thread(target=self.take_spectra, args=[datagroup],)
        self.acquisition_thread.start()



    def _button_stop_fired(self):
        """Stop the current acquisition"""
        self.taking_spectra = False

    def _button_new_file_fired(self):
        df = nplab.current_datafile()
        df.close()
        df = nplab.current_datafile()

    def take_spectra(self, datagroup):
        """Take the defined number of spectra, and save to the file."""
        N = 0
        M1 = 1
        M2 = 1
        self.taking_spectra = True
        try:
            while N < self.number_of_spectra and self.taking_spectra:
                if N!=0: # starting data collection immediately
                    time.sleep(self.time_interval_seconds)
                ds = datagroup.create_dataset("spectrum_%d",
                                 data=self.spectrometer.read_spectrum(bundle_metadata=True),
                                 attrs=self.spectrometer.get_metadata(),
                                 )
                ds.attrs.create("time_interval", self.time_interval_seconds)
                ds.attrs.create("information", self.info_string)
                datagroup.file.flush()
                print "Spectra %d of %d recorded" % (N,self.number_of_spectra)
                if N=uv_start_delay:
                    self.shutter.open_shutter()
                if N=uv_start_delay+(M1*uv_on_period):
                    self.shutter.close_shutter()
                    M1 +=1
                if N=uv_start_delay+((M1-1)*uv_on_period)+(M2*uv_off_period):
                    self.shutter.open_shutter()
                    M2 += 1
            print "Done!\n"
        finally:
            self.taking_spectra = False

class DummySpectrometer():
    def read_spectrum(self, bundle_metadata=False):
        return np.random.random(100)
    def get_metadata(self):
        return {'a':2,'b':3}

# example code:
if __name__ == "__main__":
    try:
#        spectrometer = DummySpectrometer()
        spectrometer = OceanOpticsSpectrometer(0)
        timelapse = SpectrometerTimelapse(spectrometer)
        nplab.utils.gui.show_guis([
                    spectrometer,
                    nplab.current_datafile(),
                    timelapse,
                    ])
    except seabreeze.OceanOpticsError as error:
        print "An error occurred with the spectrometer: %s" % error
    finally:
        try:
            seabreeze.shutdown_seabreeze()
        except:  # of course, if there's no spectrometer this will fail, hence the error handling
            print """The spectrometer did not reset nicely."""
