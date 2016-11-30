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
import numpy
import datetime
#from picoscope import PicoScope5000a

error_log = [[],[],[]]
error_log[0] = ['Timestamp']
error_log[1] = ['Error']
error_log[2] = ['Info']

if __name__ == "__main__":
    print(__doc__)
    verbosity = True
    error_log = [[]]
    error_log.append('Timestamp')
    error_log.append('Error')
    error_log.append('Info')
    
    """Create the laser class"""
    laser = SolsTiS.SolsTiS(('192.168.1.222', 39933))
    
    """Make sure all the wavelength locking mechanisms are off"""
    laser.cavity_lock('off')
    laser.etalon_lock('off')
    laser.wave_lock('off')
    
    """Create the wavelength meter class"""
    wlm = wavemeter.Wavemeter()
    j = 0
    while j <= 5:
        try:
            wlm.version
            j = 99;
        except wavemeter.WavemeterException as e:
            error_log[0].append(datetime.datetime.now())
            error_log[1].append(e.message)
            error_log[2].append('wlm.version')
            if e.message == "ErrWlmMissing":
                print "ERROR! " + e.message + ': try again\n'
            else:
                raise e
            # wait for the wavemeter to finish initialising
            time.sleep(1) # seconds
            j = j + 1
    
    """Start the wavelength meter"""
    wlm.active = 1
    
    """Specify the wavelength parameters in nm"""
    wavelength_start =   725.0
    wavelength_end =     975.0
    wavelength_step =       1.0    
    wavelength_precision = 0.001 
    wavelength_accuracy  = 2.0
    
    wavelength_precision_changed = False    
    wavelength_accuracy_changed = False
   
#    """Specify the time to wait after setting the wavelength"""
#    speed = 40 # nm/s
#    wait_time = 4 * wavelength_step / speed # seconds
    
    if verbosity:
        print '\n'
        print 'Wavelength start (nm):     ' + str(wavelength_start)
        print 'Wavelength end (nm):       ' + str(wavelength_end)
        print 'Wavelength step (nm):      ' + str(wavelength_step)
        print 'Wavelength precision (nm): ' + str(wavelength_precision)
        print 'Wavelength accuracy (nm):  ' + str(wavelength_accuracy)
    #    print 'Speed  (nm/s): ' + str(speed)
    #    print 'Wait time (s): ' + str(wait_time)
        print '\n'
    
    wlm.verbosity = False
    laser.verbosity = False
    verbosity = False
    
#    """ Bring the laser to the intial wavelength of the scan"""
#    laser.change_wavelength(wavelength_start)
#    while abs(wavelength_start - wlm.wavelength) > wavelength_accuracy:
#        time.sleep(0.5) # seconds
    current_wavelength = laser.laser_status['wavelength']
    wavelength = wavelength_start
    set_wavelengths = []
    measured_wavelengths = []    
    while wavelength <= wavelength_end:   
        verbosity = False
        
        set_wavelengths.append(wavelength)    
#        current_wavelength = wlm.wavelength 
        maximum_iterations = max(int(abs(current_wavelength - wavelength)/2), 20)
        if verbosity:
            print 'maximum iterations allowed: ' + str(maximum_iterations)
            
        verbosity = False
        
        """Set the laser wavelength"""
        laser.change_wavelength(wavelength)
        
        """Wait for the laser to change its wavelength"""
#        time.sleep(wait_time) # seconds            
        previous_wavelength = []
        while len(previous_wavelength) <= 3:
            previous_wavelength.append(0.0)
#        current_wavelength = wlm.wavelength    
        i = 0
        while (abs(current_wavelength - wavelength) \
              >= wavelength_accuracy or \
              abs(current_wavelength - numpy.mean(previous_wavelength)) \
              >= wavelength_precision) and \
              i <= maximum_iterations:                        
                  
            del previous_wavelength[0]
            previous_wavelength.append(current_wavelength)
            time.sleep(0.1) # seconds
            
            """Read the wavelength from the wavemeter"""
            j = 0
            while j <= 5:
                try:
                    current_wavelength = wlm.wavelength
                    j = 99;
                except wavemeter.WavemeterException as e:
                    error_log[0].append(datetime.datetime.now())
                    error_log[1].append(e.message)
                    error_log[2].append('wavelength = ' + str(wavelength) + 'nm')
                    if e.message == "ErrBigSignal":
                        print "ERROR! " + e.message + ': try again\n'
                        print datetime.datetime.now()
                    else:
                        raise e
                    time.sleep(0.1) # seconds
                    j = j + 1
                    
            if verbosity:
                print str(i) + ': ' + str(current_wavelength)
                
            i = i + 1
            
        if abs(current_wavelength - wavelength) \
           <= wavelength_accuracy:
            if verbosity:
                print 'reached accuracy'
            
        if abs(current_wavelength - numpy.mean(previous_wavelength)) \
           <= wavelength_precision:
            if verbosity:
                print 'reached precision'
            
        verbosity = False
        if verbosity:
            print 'no. of iterations required: ' + str(i)
            
        verbosity = True
        if i >= maximum_iterations:
            if abs(current_wavelength - wavelength) \
               >= wavelength_accuracy:
                wavelength_accuracy = wavelength_accuracy + 0.5
                wavelength_accuracy_changed = False
                if verbosity:
                    print 'reached maximum iterations'
                    print 'change wavelength accuracy (nm): ' + str(wavelength_accuracy)
                
            if abs(current_wavelength - numpy.mean(previous_wavelength)) \
               >= wavelength_precision:
                wavelength_precision = wavelength_precision * 10
                wavelength_precision_changed = False    
                if verbosity:
                    print 'reached maximum iterations'
                    print 'change wavelenth precision (nm): ' + str(wavelength_precision)
            
        verbosity = True
        if verbosity:
            print 'set wavelength (nm): ' + str(wavelength)
            print 'wlm wavelength (nm): ' + str(current_wavelength) + '\n'
        
        # for some reason the laser.system_status() takes ver long to update
#        laser.system_status()
#        print 'TiS wavelength (nm): ' + str(laser.laser_status['wavelength'])
        
        measured_wavelengths.append(current_wavelength)   
        wavelength = wavelength + wavelength_step 
        
   
#    print wlm.wavelength
#    print laser.laser_status['wavelength']
#    laser.system_status()
#    print laser.laser_status['wavelength']

    
    """Stop the wavelength meter"""
    wlm.verbosity = True    
    wlm.active = 0
    wlm.verbosity = False
    
    verbosity = True
    time.sleep(2) # seconds
    if verbosity:
        if wavelength_precision_changed:
            print 'Wavelength precision (nm): ' + str(wavelength_precision)
        if wavelength_accuracy_changed:        
            print 'Wavelength accuracy (nm):  ' + str(wavelength_accuracy)
        print error_log
