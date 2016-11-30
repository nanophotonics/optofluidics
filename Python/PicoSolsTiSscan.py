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
#import os
import time
from wlm import wavemeter
import SolsTiS
import numpy
import datetime
from picoscope import PicoScope5000a



if __name__ == "__main__":
    print(__doc__)
    verbosity = True
    error_log = [[],[],[]]
    error_log[0] = ['Timestamp']
    error_log[1] = ['Error']
    error_log[2] = ['Info']
    
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

    """Create the picoscope class"""
    ps = PicoScope5000a.PicoScope5000a()    
    
    """Set up the picoscope channels"""
    ps_channels = ['A', 'B']
    ps_channels_range = []
    ps_data = []
    for item in ps_channels:
        channelRange = ps.setChannel(channel=item,                                  
                                     coupling="DC", 
                                     VRange=10.0, 
                                     VOffset=0.0,
                                     enabled=True,                                  
                                     BWLimited=False,
                                     probeAttenuation=1.0,
                                     )        
        ps_channels_range.append(channelRange)    
        ps_data.append([])                           
    
    """ Set picoscope sampling interval and timebase"""
    waveform_duration = 1 # seconds
    number_of_samples = 4096
    sampling_interval = waveform_duration / number_of_samples
    (actualSamplingInterval, nSamples, maxSamples) = \
        ps.setSamplingInterval(sampling_interval, waveform_duration)

    """ Set picoscope trigger"""
    ps.setSimpleTrigger(trigSrc=ps_channels[0], 
                        threshold_V=0.5, 
                        direction='Rising', 
                        delay=0,
                        timeout_ms=100, 
                        enabled=True
                        )
        
    if verbosity:
#        for i in range(len(ps_channels)
#            print ps_channels[i] + ' channel range = ' + str(ps_channels_range[i]) + ' V'
        print("Sampling interval = %f ns" % (actualSamplingInterval * 1E9))
        print("Taking  samples = %d" % nSamples)
        print("Maximum samples = %d" % maxSamples)
        print '\n'
                        
    
    """Start the wavelength meter"""
    wlm.active = 1
    
    """Specify the wavelength parameters in nm"""
    wavelength_start = 725.0
    wavelength_end = 975.0
    wavelength_step = 10.0    
    wavelength_precision = 0.001 
    wavelength_accuracy = 2.0
    
    wavelength_precision_changed = False    
    wavelength_accuracy_changed = False
   
   
    if verbosity:
        print '\n'
        print 'Wavelength start (nm):     ' + str(wavelength_start)
        print 'Wavelength end (nm):       ' + str(wavelength_end)
        print 'Wavelength step (nm):      ' + str(wavelength_step)
        print 'Wavelength precision (nm): ' + str(wavelength_precision)
        print 'Wavelength accuracy (nm):  ' + str(wavelength_accuracy)
        print '\n'
    
    wlm.verbosity = False
    laser.verbosity = False
    verbosity = False
    

    current_wavelength = laser.laser_status['wavelength']
    wavelength = wavelength_start
    set_wavelengths = []
    measured_wavelengths = []    
    while wavelength <= wavelength_end:   
        verbosity = False
        
        set_wavelengths.append(wavelength)    
        maximum_iterations = max(int(abs(current_wavelength - wavelength)/2), 20)
        if verbosity:
            print 'maximum iterations allowed: ' + str(maximum_iterations)
            
        verbosity = False
        
        """Set the laser wavelength"""
        laser.change_wavelength(wavelength)
        
        """Wait for the laser to change its wavelength"""           
        previous_wavelength = []
        while len(previous_wavelength) <= 3:
            previous_wavelength.append(0.0) 
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
        
        # for some reason the laser.system_status() takes very long to update
#        laser.system_status()
#        print 'TiS wavelength (nm): ' + str(laser.laser_status['wavelength'])
        
        measured_wavelengths.append(current_wavelength)   
        wavelength = wavelength + wavelength_step 
        
        """Collect data from picoscope"""
        ps.runBlock()
        ps.waitReady()
        
        for i in range(len(ps_channels)):
            (data, numSamplesReturnedA, overflowA) = \
                ps.getDataRaw(channel=ps_channels[i], 
                              numSamples=0, 
                              startIndex=0, 
                              downSampleRatio=1,
                              downSampleMode=0, 
                              segmentIndex=0, 
                              data=None,
                              )
            ps_data[i].append(data)
                                                            
    # end of wavelength for loop

    
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
        if len(error_log[0]) > 1:            
            for i in range(len(error_log[0])):
                for j in range(len(error_log)):
                    print error_log[j][i]
                print '\t'

    """Write data to files"""
    
#    directory = 'R:\\aa938\\'
    directory = 'C:\\Users\\Ana Andres\\Documents\\Python Scripts\\'

    for i in range(len(ps_data)): # each channel
        file_name = 'test' + ps_channels[i] + '.txt'
        file = open(directory + file_name,'w')
        
        file.write('SolsTiS IP = ' + laser.computerIP + '\n')    
        
        file.write('Wavelength meter units = nm\n')    
        
        file.write('PicoScope model  = ' + ps.get_model_name() + '\n')
        file.write('Channel = ' + ps_channels[i] + '\n')
        file.write('Units = V\n')
        file.write('Waveform duration = ' + str(waveform_duration) + ' s\n')
        file.write('Number of samples = ' + str(nSamples) + '\n')
        file.write('Sampling interval = ' + str(actualSamplingInterval * 1E9) + ' ns\n')
        
        for item in measured_wavelengths: # each wavelength
            file.write(str(item))
            file.write('\t')
        file.write('\n')    
            
        for j in range(len(ps_data[i][0])): # each time
            for k in range(len(ps_data[i])): # each wavelength
                file.write(str(ps_data[i][k][j]))
                file.write('\t')
            file.write('\n')    
        file.close()
