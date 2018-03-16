# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 14:22:20 2018

@author: Ana Andres
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from MatplotlibSettings import *

def get_item_names(data_group):
    items = []
    for index in range(len(data_group.items())):
        items.append(data_group.items()[index][0])
    return items

folder_path = 'C:/Users/Ana Andres/Desktop/'
file_name = '2018-03-13.h5'
group_name = 'waveplate_scan_2'
y_pixel = 600

file_path = folder_path + file_name
data_file = h5py.File(file_path, 'r')

data_group = data_file[group_name]
#item_names = get_item_names(data_group)
item_names = []
for i in range(len(get_item_names(data_group))):
    item_names.append('image_' + str(i))

cmap = plt.get_cmap('jet')
colours = cmap(np.linspace(0,1,len(item_names)+1))
colour_dict = dict(zip(item_names,colours))

fig = plt.figure()
ax = fig.add_subplot(111)

for item_name in item_names:    
    item = data_group[item_name]
    image = item.value
    attributes = item.attrs
    exposure_time = attributes['exposure_time'] 
    
#    plt.imshow(image)
#    normalisation = exposure_time
    normalisation = float(max(image[y_pixel,:]))
    
    ax.plot(image[y_pixel,:]/normalisation, c=colour_dict[item_name], label=item_name)

ax.legend()
ax.set_title(group_name)
ax.set_xlabel('pixels')
ax.set_ylabel('intensity (a.u.)')
ax.set_xlim(300,750)
    
    