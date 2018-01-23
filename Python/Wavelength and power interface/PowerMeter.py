# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 12:17:13 2018

@author: Hera
"""

import visa
from ThorlabsPM100 import ThorlabsPM100
rm = visa.ResourceManager()
inst = rm.open_resource('USB0::0x1313::0x8078::P0011774::INSTR',timeout=1)
power_meter=ThorlabsPM100(inst=inst)
power_meter.system.beeper.immediate()
print(power_meter.read)