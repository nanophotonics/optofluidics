__author__ = 'chrisgrosse'


from PyCmdMessenger import CmdMessengerThreaded, ArduinoDueBoard

import numpy as np
import pyqtgraph as pg
from collections import deque

from nplab.instrument import Instrument
from nplab.utils.gui import QtWidgets, get_qt_app, QtCore
from nplab import datafile


class ArduinoDueOptofluidics(CmdMessengerThreaded, Instrument, QtCore.QObject):

    # definition of CmdMessenger commands
    commands = [["check_arduino",""],  
                ["ret_fast_intensity","ff"],
                ["sweep_delay","iif"], # param: number of data points, averaged data points, delay increment in seconds
                ["set_delay","f"],
                ["get_delay","f"],
                ["get_intensity_trace","iiif"],
                ["osci_mode","i"],
                ["set_laser_on","?"],
                ["get_laser_on","?"],
                ["ret_string","s"],
                ["ret_int","i"],
                ["ret_char","c"],
                ["ret_float","f"],
                ["ret_double","d"],
                ["ret_bool","?"],
                ["ret_byte","b"]]

    intensity_signal = QtCore.Signal(float,float)

    def __init__(self, port):
        # get Arduino board/serial instance
        self.DueBoard = ArduinoDueBoard(port, baud_rate=115200)
        # initialize PyCmdMessenger and start serial reading thread
        CmdMessengerThreaded.__init__(self, self.DueBoard, self.commands)
        QtCore.QObject.__init__(self)
        self.check_arduino_response = True
        self.incomming_data = []  
        # open a data file and set it as "current"
        self.data_file = datafile.current()  
        self.get_qt_ui()
#        self.data = np.zeros(0)    
        self.attributes = dict()      

    def __del__(self):
        if isinstance(self.data_file, datafile.DataFile):
            self.data_file.close()
        self.close()
        
    def init_scan(self, num_points):
        self.data_index=0
        self.data_delay = np.zeros(num_points)
        self.data_signal1 = np.zeros(num_points)
        self.data_delay.fill(np.nan)
        self.data_signal1.fill(np.nan)        
            
    def close_scan(self):
        curr_scan = self.datafile_group.create_group('scan_%d', attrs=self.attributes)
        curr_scan.create_dataset('delay', data=self.data_delay)
        curr_scan.create_dataset('signal1', data=self.data_signal1)
        
    def osci_trace(self, num_points, description=None):
        self.init_scan(num_points)
        self.send("osci_mode", num_points)
        self.datafile_group = self.data_file.require_group('oscilloscope_traces') # create new group in data file        
#        self.curr_dataset = datafile_group.create_dataset('scan_%d', attrs=dict(description=description)) 
#        self.curr_dataset = datafile_group.create_group('scan_%d', attrs=dict(description=description))

    def delay_trace(self, num_points, num_pulses, time_increment_in_sec):
        self.init_scan(num_points)
        self.send("sweep_delay", num_points, num_pulses, time_increment_in_sec)
        self.datafile_group = self.data_file.require_group('delay_traces')
    
    def pulse_intensity_trace(self, num_points, num_pulses, points_per_pulse, time_delay_in_sec):
        self.init_scan(num_points)
        self.send("get_intensity_trace", num_points, num_pulses, points_per_pulse, time_delay_in_sec)
        self.datafile_group = self.data_file.require_group('pulse_intensity_traces')
        
    def laser_on(self, status):
        self.send("set_laser_on",status)
              

    def response_to_command(self, cmd_name, msg, message_time):
        if cmd_name == 'ret_fast_intensity':
            self.data_delay[self.data_index]=msg[1]
            self.data_signal1[self.data_index]=msg[0]

            self.intensity_signal.emit(msg[0],msg[1])
            self.data_index+=1
            if self.data_delay.shape[0] == self.data_index:
                self.close_scan()
        else:
            print cmd_name, msg
            
    def get_qt_ui(self):
        self.gui=ArduinoUI(self)
        return self.gui

class ArduinoUI(QtWidgets.QWidget):
    def __init__(self, arduino_instance):
        assert isinstance(arduino_instance, ArduinoDueOptofluidics) ,\
            "experiment must be an instance of ArduinoDueOptofluidics"
        super(ArduinoUI, self).__init__()
        self.due =arduino_instance
        self.signal1_deque = deque(maxlen=500)#2e3)
#        self.signal1_deque.append(0.0)
        self.signal2_deque = deque(maxlen=500)
#        self.signal2_deque.append(0.0)
        self.xaxis_deque = deque(maxlen=500)
        self.init_ui()
        self.update_counter=0

    def init_ui(self):
        self.setWindowTitle('signal monitor')
        self.signal1_plot = pg.PlotWidget(labels = {'left':'intensity','bottom':'time'},symbol='o', pen=None)
        self.signal1_plot.enableAutoRange()
        self.signal2_plot = pg.PlotWidget(labels = {'left':'intensity','bottom':'time'})
        self.signal2_plot.enableAutoRange()

        # create grid layout to manage widgets size
        layout = QtWidgets.QGridLayout()
        self.setLayout(layout)
        layout.addWidget(self.signal1_plot,0,0)
        layout.addWidget(self.signal2_plot,1,0)
#        self.signal1_plot.plot(y=self.signal1_deque, x=self.xaxis_deque)
#        self.signal2_plot.plot(y=self.signal2_deque, x=self.xaxis_deque)

        # event handling
        self.due.intensity_signal.connect(self.update_gui)

    def update_gui(self,value1,value2):
        self.signal1_deque.append(value1)
#        self.signal2_deque.append(value2)
        self.xaxis_deque.append(value2)
#        if self.update_counter == 9:
        self.signal1_plot.clear()
        self.signal1_plot.plot(self.xaxis_deque, self.signal1_deque, symbol='o')
        self.signal2_plot.clear()
        self.signal2_plot.plot(self.xaxis_deque, symbol='o')
        self.update_counter=0
        
    def clear(self):
        self.signal1_deque=[]
        self.signal2_deque=[]
        self.xaxis_deque=[]

if __name__ == '__main__':

    import sys
    from nplab.utils.gui import get_qt_app 
    due = ArduinoDueOptofluidics("COM9")
    due.show_gui(blocking=False)
    app = get_qt_app()