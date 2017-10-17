__author__ = 'chrisgrosse'


from PyCmdMessenger import CmdMessengerThreaded, ArduinoDueBoard

import pyqtgraph as pg
from collections import deque

from nplab.instrument import Instrument
from nplab.utils.gui import QtWidgets, get_qt_app, QtCore


class ArduinoDueOptofluidics(CmdMessengerThreaded, Instrument, QtCore.QObject):

    # definition of CmdMessenger commands
    commands = [["check_arduino",""],                
                ["get_intensity",""],
                ["ret_intensity","fff"],
                ["ret_fast_intensity","ff"],
                ["set_cont_signal_readout","?"],
                ["get_cont_signal_readout",""],
                ["set_daq_freq","i"],
                ["get_daq_freq",""],
                ["sweep_delay","iif"], # param: number of data points, averaged data points, delay increment in seconds
                ["fast_ADC_read","i"],
                ["ret_string","s*"],
                ["ret_int","i*"],
                ["ret_char","c*"],
                ["ret_float","f*"],
                ["ret_double","d*"],
                ["ret_bool","?*"],
                ["ret_byte","b*"]]

    intensity_signal = QtCore.Signal(float,float)

    def __init__(self, port):
        # get Arduino board/serial instance
        self.DueBoard = ArduinoDueBoard(port, baud_rate=115200)
        # initialize PyCmdMessenger and start serial reading thread
        CmdMessengerThreaded.__init__(self, self.DueBoard, self.commands)
        QtCore.QObject.__init__(self)
        self.check_arduino_response = True

    def __del__(self):
        self.close()

    def response_to_command(self, cmd_name, msg, message_time):
#        print "cmd_name:", cmd_name
#        print "msg:", msg
#      if cmd_name == 'ret_intensity':  # received intensity values
#            self.intensity_signal.emit(msg[0],msg[1])
        if cmd_name == 'ret_fast_intensity':
#            print "msg:",msg
            self.intensity_signal.emit(msg[0],msg[1])
        else:
            print cmd_name, msg[0]

    def get_qt_ui(self):
        return ArduinoUI(self)

class ArduinoUI(QtWidgets.QWidget):
    def __init__(self, arduino_instance):
        assert isinstance(arduino_instance, ArduinoDueOptofluidics) ,\
            "experiment must be an instance of ArduinoDueOptofluidics"
        super(ArduinoUI, self).__init__()
        self.due =arduino_instance
        self.signal1_deque = deque(maxlen=100)#2e3)
#        self.signal1_deque.append(0.0)
        self.signal2_deque = deque(maxlen=2e3)
#        self.signal2_deque.append(0.0)
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
        self.signal1_plot.plot(self.signal1_deque)
        self.signal2_plot.plot(self.signal2_deque)

        # event handling
        self.due.intensity_signal.connect(self.update_gui)

    def update_gui(self,value1,value2):
#        print "value1:", value1, value2
        self.signal1_deque.append(value1)
        self.signal2_deque.append(value2)
#        if self.update_counter == 9:
        self.signal1_plot.clear()
        self.signal1_plot.plot(self.signal1_deque,symbol='o')
        self.signal2_plot.clear()
        self.signal2_plot.plot(self.signal2_deque,symbol='o')
        self.update_counter=0

if __name__ == '__main__':

    import sys
    from nplab.utils.gui import get_qt_app 
    due = ArduinoDueOptofluidics("COM11")
    due.show_gui(blocking=False)
    app = get_qt_app()