# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 13:10:04 2018

@author: Ana Andres
"""

import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from MatplotlibSettings import *

#from 'C:\Users\Ana Andres\Documents\GitHub\Python-GarminDataAnalyser\DatabaseRead.py' import ElapsedTime
#from Python-GarminDataAnalyser.DatabaseRead import ElapsedTime

def ElapsedTime(timestamp, units_t='sec', mode='start'):    
    """Calculate the elapsed time from a timestamp pandas Series.
        
    Arguments:
    timestamp : timestamp pandas Series
        Timestamp values
    units_d: string
        Units of the calculated time, e.g. 'sec' (default), 'min', or 'h'
    mode : string
        If 'start' calculate the elapsed time between all points of the array and the first point.
        If 'previous' calculate the elapsed time between all points of the array the previous point.
    
    Output:
    elapsed_time: float pandas Series
        Contains the calculated elapsed time in units of units_t
    """
    
    # The Garmin Forerunner 35 takes data every 1 second

    origin_time = np.empty(timestamp.shape, dtype=type(timestamp))
    if mode == 'start':
        origin_time[:] = timestamp[0]
    
    elif mode == 'previous':
        origin_time[0] = timestamp[0]
        for i, time in enumerate(timestamp[0:-1]):
            origin_time[i+1] = time
            # TODO: change to origin_time.append(time)
            
    else:
        raise ValueError('Unable to recognise the mode.')  

    timedelta = timestamp-origin_time
    elapsed_time = timedelta.astype('timedelta64[ms]')/1000 # in seconds
    
    if units_t == 's':
        pass    
    elif units_t == 'sec':
        pass
    elif units_t == 'min':
        # Convert seconds to minutes
        elapsed_time = elapsed_time/60 
    elif units_t == 'h':
        # Convert seconds to hours
        elapsed_time = elapsed_time/60/60 
    else:
        raise ValueError('Unable to recognise the units for the time.')    
    return elapsed_time

def get_item_names(data_group):
    items = []
    for index in range(len(data_group.items())):
        items.append(data_group.items()[index][0])
    return items
    
def elapsed_time(series):
    origin_time = [series[0]]
    for time in series[0:-1]:
        origin_time.append(time)
    timedelta = series - pd.Series(origin_time)
    elapsed_time = timedelta.astype('timedelta64[ms]')/1000 # in seconds
    return elapsed_time

def populate_plot_options(ax, x='capture_timestamp', y='power_w', 
                          kind='line', label=False):    
    plot_options=dict()
    
    plot_options['ax'] = ax
    plot_options['kind'] = kind
    
    if plot_options['kind'] == 'line':
        plot_options['style'] = '.-'
#        plot_options['markersize'] = 30
        plot_options['x'] = x
        plot_options['y'] = y
        
    elif plot_options['kind'] == 'scatter':
        plot_options['x'] = x
        plot_options['y'] = y
    
    if label:
        plot_options['label'] = label
        
    return plot_options
    
    
def get_attributes_dataframe(items):        
    # NOTE: there will be errors if not all items have the same attributes
    attributes_dict = dict()    
    for i, item in enumerate(items):
        attributes = item[1].attrs        
        if i == 0: attributes_dict['item'] = []
        attributes_dict['item'].append(item[0])
        for key in attributes.keys():
            if i == 0: attributes_dict[key] = []
            attributes_dict[key].append(attributes[key])
    df = pd.DataFrame.from_dict(attributes_dict)
    timestamp_columns = ['creation_timestamp', 'capture_timestamp']
    for column in timestamp_columns:
        if column in df.columns:
            df[column] = pd.to_datetime(df[column])
    return df
    
    
def plot_video(data_group):
    all_video_names = get_item_names(data_group)
    desired_video_names = all_video_names
#    desired_video_names = ['video_1']    
    video_names = np.sort(list(set(all_video_names) & set(desired_video_names)))
    
    timestamp_columns = [
#                         'creation_timestamp', 
#                         'capture_timestamp', 
#                         'pre_capture_timestamp', 
#                         'post_capture_timestamp',
                         ]
    
#    fig1 = plt.figure()
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    for i, video_name in enumerate(video_names):
        print video_name
        
        video_group = data_group[video_name]
        frames = video_group.items()
        data_frame = get_attributes_dataframe(frames)
#        data_frame.sort_values(by='capture_timestamp', inplace=True)
#        data_frame.sort_values(by='pre_capture_timestamp', inplace=True)
        data_frame.sort_values(by='pre_capture_time', inplace=True)
        data_frame['i'] = range(len(data_frame))
        
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
#        subplot = int(str(len(video_names)) + '1' + str(i+1))
#        ax = fig1.add_subplot(subplot)               
        
        for column in list(set(timestamp_columns) & set(data_frame.columns)):
            data_frame[column] = pd.to_datetime(data_frame[column])
            data_frame[column+'_elapsed_time'] = ElapsedTime(data_frame[column])  
            delta = ElapsedTime(data_frame[column], mode='previous')  
#            delta = elapsed_time(data_frame[column])
            data_frame[column+'_delta_sec'] = delta
            data_frame[column+'_framerate'] = 1/delta            
            
            plot_options = populate_plot_options(ax=ax,
#                                                 x=column, 
#                                                 x='i',
                                                 x=column+'_elapsed_time', 
                                                 y=column+'_framerate',
                                                 )
            data_frame.plot(**plot_options)                
          
#        time_reference = data_frame['pre_capture_time']
#        delta = data_frame['pre_capture_time'] - data_frame.loc[0,'pre_capture_time']
#        a = data_frame['pre_capture_time']
#        delta = [x - a[i - 1] for i, x in enumerate(a) if i > 0]
#        data_frame.loc[1:,'pre_capture_framerate'] = 1/np.array(delta)
        delta = np.ediff1d(data_frame['pre_capture_time'])
        data_frame.loc[1:,'pre_capture_framerate'] = 1/delta
        plot_options = populate_plot_options(ax=ax,
#                                             x='capture_timestamp_elapsed_time', 
#                                             x='pre_capture_timestamp_elapsed_time', 
                                             x='pre_capture_time', 
                                             y='pre_capture_framerate')
        data_frame.plot(**plot_options)
        
        plot_options = populate_plot_options(ax=ax,
#                                             x='capture_timestamp_elapsed_time', 
#                                             x='pre_capture_timestamp_elapsed_time', 
                                             x='pre_capture_time', 
                                             y='framerate')
        data_frame.plot(**plot_options)
#        ax.set_xlabel('timestamp')       
        ax.set_xlabel('time (s)')
        ax.set_title(video_name)
#        ax.set_xlim(0.1,0.2)
#        ax.set_ylim(120,180)
        
#        difference = data_frame['creation_timestamp'] - data_frame['capture_timestamp']
#        difference = data_frame['post_capture_timestamp'] - data_frame['pre_capture_timestamp']
#        difference = data_frame['creation_timestamp'] - data_frame['pre_capture_timestamp']
        difference = data_frame['post_capture_time'] - data_frame['pre_capture_time']
        data_frame['time_difference_sec'] = difference # in seconds
#        data_frame['time_difference_sec'] = difference.dt.microseconds/1e6 # in seconds
        
        plot_options = populate_plot_options(ax=ax2,
#                                             x='capture_timestamp_elapsed_time', 
#                                             x='pre_capture_timestamp_elapsed_time',
                                             x='pre_capture_time',
                                             y='time_difference_sec',
                                             label=video_name)                
        data_frame.plot(**plot_options)
#        ax2.set_ylabel('time difference (sec): creation - capture')
        ax2.set_ylabel('capture time diff (sec): post - pre')
        ax2.set_xlabel('time (s)')
#        ax2.set_xlim(0,1)        
    
    return data_frame
    
def plot_waveplate_scan(data_group, ax):    
    frames = data_group.items()
    data_frame = get_attributes_dataframe(frames)
    
    x = 'waveplate_angle_deg'
    y = 'power_w'

#    x= 'capture_timestamp'
#    y = ['waveplate_angle_deg','wavelength_nm', 'power_w']    
    
    data_frame.sort_values(by=x, inplace=True)
    plot_options = populate_plot_options(ax=ax, x=x, y=y)
    ax = data_frame.plot(subplots=True, **plot_options)
    return data_frame           
     
    

if __name__ == "__main__":
#    data_file = h5py.File('R:/3-Temporary/aa938/2018.01.26 - sweep tests/2018-01-26.h5', 'r')
#    data_file = h5py.File('R:/3-Temporary/aa938/2018.01.26 - camera framerate tests/2018-01-26.h5', 'r')
    data_file = h5py.File('R:/3-Temporary/aa938/2018-01-31-new.h5', 'r')

    group_names = get_item_names(data_file)
    
    # plot data from videos
    if 'videos' in group_names:
        data_group = data_file['videos']
        video_df = plot_video(data_group)
    
    # plot data from waveplate scans
    all_waveplate_scans = [name for name in group_names if 'waveplate_scan' in name]
    desired_waveplate_scans = ['waveplate_scan_1']
    available_waveplate_scans = list(set(all_waveplate_scans) & set(desired_waveplate_scans))

    for name in available_waveplate_scans:
        if name in group_names:
            data_group = data_file[name]
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
            wp_scan_df = plot_waveplate_scan(data_group=data_group, ax=ax)   
            
            

            

    