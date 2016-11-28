# -*- coding: utf-8 -*-
"""
@author: Ana Andres

A script to scan the wavelength range of the SolsTiS laser. \
The precise wavelength is read with the HighFinesse wavemeter. \
The intensity of the photodiode is monitored with the PicoScope 5444B.

Not implemented:
- Linking/unking the wavemeter and the SolsTiS

"""

#import sys
import time
from wlm import wavemeter
import SolsTiS
#from picoscope import PicoScope5000a

if __name__ == "__main__":
    print(__doc__)
    
    """Create the laser class"""
    laser = SolsTiS.SolsTiS(('192.168.1.222', 39933))
    
    """Make sure all the wavelength locking mechanisms are off"""
    laser.cavity_lock('off')
    laser.etalon_lock('off')
    laser.wave_lock('off')
    
    """Create the wavelength meter class"""
    wlm = wavemeter.Wavemeter()
    wlm.version
    
    """Start the wavelength meter"""
    wlm.active = 1
    
    """Specify the wavelength range in nm"""
    wavelength_start = 725.0
    wavelength_end =   975.0
    wavelength_step =   50.0
    
    # Bring the laser to the intial wavelength of the scan
    laser.change_wavelength(wavelength_start)
    time.sleep(3) # seconds
    
    wlm.verbosity = False
    wavelength = wavelength_start
    set_wavelengths = []
    measured_wavelengths = []    
    while wavelength <= wavelength_end:     
        set_wavelengths.append(wavelength)    
        
        """Set the laser wavelength"""
        laser.change_wavelength(wavelength)
        
        """Wait for the laser to change its wavelength"""
        time.sleep(0.5) # seconds
               
        """Read the wavelength from the wavemeter"""
        current_wavelength = wlm.wavelength        
        measured_wavelengths.append(current_wavelength)       
        
        laser.system_status()
        print 'set wavelength (nm): ' + str(wavelength)
        print 'TiS wavelength (nm): ' + str(laser.laser_status['wavelength'])
        print 'wlm wavelength (nm): ' + str(current_wavelength) + '\n'
#        print str(set_wavelength) + '\t\t' + str(current_wavelength) + '\n'
        
        wavelength = wavelength + wavelength_step 
    
    time.sleep(3) # seconds
    wlm.verbosity = True
    
    """Stop the wavelength meter"""
    wlm.active = 0
