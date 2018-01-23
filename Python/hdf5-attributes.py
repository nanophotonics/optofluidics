# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 13:10:04 2018

@author: Ana Andres
"""

#import nplab
import h5py
import pandas as pd
import matplotlib.pyplot as plt
from MatplotlibSettings import *

def get_items(data_group):
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

if __name__ == "__main__":
    data_file = h5py.File('C:/Users/Ana Andres/Documents/GitHub/optofluidics/Python/camera tests/2018-01-22.h5', 'r')

    group_names = get_items(data_file)
    
    if 'videos' in group_names:
        videos_info = dict()
        data_group = data_file['videos']
        video_names = get_items(data_group)
        video_attributes_dict = dict()
        desired_keys = ['creation_timestamp', 'capture_timestamp', 'framerate']
        timestamp_columns = ['creation_timestamp', 'capture_timestamp']
        
        fig = plt.figure()
        
        for i, video_name in enumerate(video_names):
            video_group = data_group[video_name]
            frame_names = get_items(video_group)
            
            attributes_dict = dict()
            for key in desired_keys:
                attributes_dict[key] = []
                
            for frame_name in frame_names:
                data = video_group[frame_name]
                attributes = data.attrs      
                for key in desired_keys:
                    if key in attributes.keys():
                        attributes_dict[key].append(attributes[key])
                    else:
                        attributes_dict[key].append('nan')
                        
            subplot = int(str(len(video_names)) + '1' + str(i+1))
            ax = fig.add_subplot(subplot)
            
            plot_options=dict()
            plot_options['ax'] = ax

            plot_options['kind'] = 'line'
#            plot_options['kind'] = 'scatter'
            
            if plot_options['kind'] == 'line':
                plot_options['style'] = '.-'
                plot_options['markersize'] = 20
                plot_options['x'] = ['creation_timestamp']
                plot_options['y'] = ['creation_timestamp_framerate']
                
            elif plot_options['kind'] == 'scatter':
                plot_options['x'] = 'creation_timestamp'
                plot_options['y'] = 'framerate'
                
            data_frame = pd.DataFrame.from_dict(attributes_dict)
            for column in list(set(timestamp_columns) & set(data_frame.columns)):
                data_frame[column] = pd.to_datetime(data_frame[column], errors='ignore')
                delta = elapsed_time(data_frame[column])
                data_frame[column+'_delta_sec'] = delta
                data_frame[column+'_framerate'] = 1/delta
                
                plot_options['x'] = [column]
                plot_options['y'] = [column+'_framerate']
                
#                data_frame.plot(**plot_options)
            
            difference = data_frame['creation_timestamp'] - data_frame['capture_timestamp']
            data_frame['time_difference_sec'] = difference.dt.microseconds/1e6 # inseconds
            plot_options['y'] = 'time_difference_sec'
                
            data_frame.plot(**plot_options)
                
                        
            

            

    