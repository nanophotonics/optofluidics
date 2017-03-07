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
import matplotlib.pyplot as plt
import datetime
from picoscope import PicoScope5000a



if __name__ == "__main__":
    print(__doc__)
    verbosity = True
    error_log = [[],[],[]]
    error_log[0] = ['Timestamp']
    error_log[1] = ['Error']
    error_log[2] = ['Info']
    
    """Instantiate the laser object"""
    laser = SolsTiS.SolsTiS(('192.168.1.222', 39933))
    
    """Make sure all the wavelength locking mechanisms are off"""
    laser.cavity_lock('off')
    laser.etalon_lock('off')
    laser.wave_lock('off')
    
    """Instantiate the wavelength meter object"""
    wlm = wavemeter.Wavemeter()
    j = 0
    while j <= 5:
        try:
            wlm.version
            break
        except wavemeter.WavemeterException as e:
            error_log[0].append(datetime.datetime.now())
            error_log[1].append(e.message)
            error_log[2].append('wlm.version')
            if e.message == "ErrWlmMissing":
                print "ERROR! " + e.message + ': try again\n'
            else:
                raise e
            if j == 5:
                raise e
            # wait for the wavemeter to finish initialising
            time.sleep(1) # seconds
            j = j + 1

    """Instantiate the picoscope object"""
    ps = PicoScope5000a.PicoScope5000a()   
    
    """Specify which picoscope channels to use."""
    ps_channels = ['A', 'B']
#    ps_channels = ['B'] 
    ps_channels_range = []
    ps_data = []    
    ps_average = []    
    ps_autorange = []
    for i in range(len(ps_channels)):
        ps_channels_range.append(0.2)
        ps_data.append([])
        ps_average.append([])
        ps_autorange.append(True)
                            
    
    """ Set picoscope sampling interval and timebase"""
    waveform_duration = 1.0 # seconds
    number_of_samples = 2**10
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
        print("Waveform duration = %f s" % (nSamples * actualSamplingInterval))
        print("Sampling interval = %f ns" % (actualSamplingInterval * 1E9))
        print("Taking  samples = %d" % nSamples)
        print("Maximum samples = %d" % maxSamples)
        print('\n')
                        
    
    """Start the wavelength meter"""
    wlm.active = 1
    
    """Specify the wavelength parameters in nm"""
    wavelength_start = 725.0
    wavelength_end = 975.0
    wavelength_step = 5.0    
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
    
    """Make sure that the wavelength limits are within range for the SolsTiS"""
    wavelength_start = max(wavelength_start, 725)
    wavelength_end   = min(wavelength_end, 975)

    
    current_wavelength = laser.laser_status['wavelength']
    wavelength = wavelength_start
    set_wavelengths = []
    measured_wavelengths = []    
    while wavelength <= wavelength_end:  
        verbosity = True
        if verbosity:
            print 'set wavelength (nm): ' + str(wavelength)
            
        verbosity = False        
        set_wavelengths.append(wavelength)    
        maximum_iterations = max(int(abs(current_wavelength - wavelength)/2), 20)
        if verbosity:
            print 'maximum iterations allowed: ' + str(maximum_iterations)
            
                
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
            j_max = 9
            while j <= j_max:
                try:
                    current_wavelength = wlm.wavelength
                    break
                except wavemeter.WavemeterException as e:
                    error_log[0].append(datetime.datetime.now())
                    error_log[1].append(e.message)
                    error_log[2].append('wavelength = ' + str(wavelength) + 'nm')
                    if e.message == "ErrBigSignal" or \
                       e.message == "ErrLowSignal" or \
                       e.message == "ErrNoValue":                           
                        print "ERROR! " + e.message + ': try again (' + str(j) + ')'
                    else:
                        raise e
                    if j == j_max:
                        raise e
                    time.sleep(0.1) # seconds
                    j = j + 1
                    
            verbosity = False
            if verbosity:
                print str(i) + ': ' + str(current_wavelength)
                
            i = i + 1

        verbosity = False
        if verbosity:
            print 'no. of iterations required: ' + str(i)            

        verbosity = True        
        if verbosity:
            print 'wlm wavelength (nm): ' + str(current_wavelength) + '\n'
        
        """Append measured wavelength to the list"""
        measured_wavelengths.append(current_wavelength)   
        wavelength = wavelength + wavelength_step             
            
        verbosity = False
        if abs(current_wavelength - wavelength) \
           <= wavelength_accuracy:
            if verbosity:
                print 'Acchieved correct accuracy'
        if abs(current_wavelength - numpy.mean(previous_wavelength)) \
           <= wavelength_precision:
            if verbosity:
               print 'Achieved correct precision'          
            
        if i >= maximum_iterations:
            if abs(current_wavelength - wavelength) \
               >= wavelength_accuracy:
                wavelength_accuracy = wavelength_accuracy + 0.5
                wavelength_accuracy_changed = False
                print 'Reached maximum iterations'
                print 'Change wavelength accuracy (nm): ' + str(wavelength_accuracy)
                
            if abs(current_wavelength - numpy.mean(previous_wavelength)) \
               >= wavelength_precision:
                wavelength_precision = wavelength_precision * 10
                wavelength_precision_changed = False    
                print 'Reached maximum iterations'
                print 'Change wavelenth precision (nm): ' + str(wavelength_precision)
                
 
        verbosity = True
        """Collect data from picoscope"""   
        for i in range(len(ps_channels)):                                                                     
            ps_data[i].append([])
            ps_average[i].append([])
            ps_autorange[i] = True
        while any(ps_autorange):
            for i in range(len(ps_channels)):    
                channelRange = ps_channels_range[i]
                """Set up channel range."""
                channelRange = ps.setChannel(channel=ps_channels[i],
                                     coupling="DC", 
                                     VRange=channelRange, 
                                     VOffset=0.0,
                                     enabled=True,                                  
                                     BWLimited=False,
                                     probeAttenuation=1.0,
                                     ) 
                if verbosity:
                    print 'Channel ' + ps_channels[i] + ' range = ' + str(channelRange) + ' V'
                ps_channels_range[i] = channelRange
            if verbosity:
                print ''
            
            """Collect data."""
            ps.runBlock()
            ps.waitReady()
                      
            for i in range(len(ps_channels)):
                """Read data from the buffer."""
                (data, numSamplesReturned, overflow) = \
                    ps.getDataRaw(channel=ps_channels[i])
                    
                """Convert data to volts."""
                data = data/32512.0*ps_channels_range[i]
                if verbosity:
                    print 'Channel ' + ps_channels[i] + ' max = ' + str(max(data)) + ' V'
                
                """Modify chanel range if data is too small/big."""
                if max(abs(data)) < ps_channels_range[i]*0.3 and ps_channels_range[i] > 0.01:
                    ps_channels_range[i] = ps_channels_range[i] / 2.5
                elif max(abs(data)) > ps_channels_range[i]*0.95 and ps_channels_range[1] < 20:
                    ps_channels_range[i] = ps_channels_range[i] * 2     
                else:
                    ps_autorange[i] = False
                        
                if max(abs(data)) >= 12:
                    print 'WARNING!!! Channel ' + ps_channels[i] + ' is saturated!!!'
            
                """Append data to list."""
                ps_data[i][len(ps_data[i])-1] = data # units = volts
                ps_average[i][len(ps_average[i])-1] = numpy.mean(data)
                
            if verbosity:
                print ''
                

                                                            
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
        
    """Plot data."""        
    plt.hold(False)
    for i in range(len(ps_channels)):     
        plt.plot(measured_wavelengths, ps_average[i], 
                 linewidth = 2.0, 
                 label = 'Channel ' + ps_channels[i],
                 )
        plt.hold(True)
        
    plt.grid()    
    plot_font_size = 32
    plt.legend(fontsize = plot_font_size, loc = 1)
    now = datetime.datetime.now()
    plt.title(now, fontsize = plot_font_size)
#    plt.title('Picoscope average vs. wavelength', fontsize = plot_font_size)
    plt.xlabel('Wavemeter wavelength (nm)', fontsize = plot_font_size)
    plt.ylabel('Picoscope average (V)', fontsize = plot_font_size)
    plt.tick_params(axis='both', which='major', 
                    length = 20, width = 2, labelsize = plot_font_size)    
    plt.show()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
        
    """Write data to files"""
    save_data = True
    
    if save_data:
        directory = 'R:\\3-Temporary\\aa938\\'
#        now = datetime.datetime.now()
        file_name_core = "%04i" % now.year + '.' + \
                         "%02i" % now.month + '.' + \
                         "%02i" % now.day + '_' + \
                         "%02i" % now.hour + '.' + \
                         "%02i" % now.minute + '.' + \
                         "%02i" % now.second
    
        for i in range(len(ps_data)): # each channel
            file_name = file_name_core + '_' + ps_channels[i] + '.txt'        
            file = open(directory + file_name,'w')
            
            print 'Writing channel ' + ps_channels[i] + ' data to file...'
            
            file.write('SolsTiS IP = ' + laser.computerIP + '\n')    
            
            file.write('Wavelength meter units = nm\n')    
            
            file.write('PicoScope model  = ' + ps.get_model_name() + '\n')
            file.write('Channel = ' + ps_channels[i] + '\n')
            file.write('Units = V\n')
            file.write('Waveform duration = ' + str(nSamples*actualSamplingInterval) + ' s\n')
            file.write('Number of samples = ' + str(nSamples) + '\n')
            file.write('Sampling interval = ' + str(actualSamplingInterval * 1E9) + ' ns\n')
            
            file.write('\n')    
            
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

    ps.close()