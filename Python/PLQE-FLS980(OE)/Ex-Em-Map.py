
# coding: utf-8

# In[1]:


import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import plotly.plotly as py
import plotly.graph_objs as go
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from IPython.display import HTML, display
import tabulate
import pandas as pd
import os
color = cm.viridis
jet = plt.get_cmap('viridis')
values = range(100)
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='nipy_spectral')
TUM=['#005c96','#a0b239','#ed762a']
OPT=['#347f3a','#36358b','#e47327']
matplotlib.pyplot.set_cmap('prism')
get_cmap()


# In[2]:


path='C://Users//tl457//OneDrive - University Of Cambridge 1//Measurements//180914-N-CND(gr)-Cuvette-PLQE//'
plotpath='C://Users//tl457//OneDrive - University Of Cambridge 1//Plots//'


# In[3]:


def SG(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
# precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
# pad the signal at the extremes with
# values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    return np.convolve( m[::-1], y, mode='valid')

def WAVE(key,number):
    dataset=data[key]
    wavelengths = dataset['spectrum_{}'.format(number)].attrs['wavelengths']
    return wavelengths

def SGW(array, window):
    return array[window/2:-window/2]


# In[4]:


read_row_start=24 #start reading the txt file at
read_column_end=8 #end reading the txt file at

with open(path+'n-cnd(gr)_EEMAP.txt') as f:
    reader = csv.reader(f, delimiter="\t")
    results = list(reader)

d1=[]
exwav=[]
emwav=[]

for i in range(read_row_start,len(results)):
    emwav.append(float(results[i][0]))
for j in range(1,len(results[read_column_end])-1):
    exwav.append(float(results[read_column_end][j]))
for i in range(24,len(results)):
    for j in range(1,len(results[read_row_start])-1):
        d1.append(results[i][j].replace("E", "e"))

d1=map(float,d1)
d1=map(abs,d1)
d2=np.array(d1).reshape(len(emwav),len(exwav))


# In[5]:


#plot parameters
SG_Window=7 #odd number needed
SG_Order=2

cpsMAX=1.05*(np.amax(d2))
cpsMIN=-0
emextentmin=340
emextentmax=800
exextentmin=340
exextentmax=600

#cross section parameters
exavg=400
emMIN=497
emMAX=503

em_avg=(emMIN+emMAX)/2


# In[6]:


with rc_context(fname='plotting'):
    fig = plt.figure(figsize=(6, 6))
    grid = plt.GridSpec(4, 4, hspace=0.4, wspace=0.4, right=0.75, bottom=0.25)
    main_ax = fig.add_subplot(grid[:-1, 1:])
    y_plot = fig.add_subplot(grid[:-1, 0], sharey=main_ax)
    x_plot = fig.add_subplot(grid[-1, 1:], sharex=main_ax)

# main axes
with rc_context(fname='plotting'):
    cb=main_ax.imshow(d2,cmap='jet',extent=(np.amin(exwav), np.amax(exwav), np.amin(emwav), np.amax(emwav)),vmin=cpsMIN,vmax=cpsMAX,origin='lower',interpolation='spline16',aspect='auto')
    main_ax.axvline(exavg,color='white',linewidth=1,linestyle='-')
    main_ax.axhline(em_avg,color='white',linewidth=1,linestyle='-')
    main_ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

# trace at specific emission wavelength
trace=[]
for i in range(0,len(exwav)):
    if int(exwav[i])==int(exavg):
        number=i
    else:
        pass
for j in range (0,len(emwav)):
    trace.append(d2[j][number])

trace_SG=SG(trace,SG_Window,SG_Order)
waves_SG=SGW(emwav,SG_Window-1)

with rc_context(fname='plotting'):
    y_plot.plot(trace_SG,waves_SG,color=cm.jet(60),label='{}'.format('( Ex: '+str(exavg)+') nm'))
    y_plot.set_ylim(emextentmin,emextentmax)
    #y_plot.legend(fontsize=6)
    y_plot.set_xlim(cpsMIN,cpsMAX)
    #y_plot.tick_params(rotation=45)
    y_plot.axhline(em_avg,color='grey',linewidth=1,linestyle='-')
    y_plot.set_ylabel('Emission Wavelength / nm')


# trace at specific excitation wavelength
for i in range(0,len(emwav)):
    if int(round(emwav[i]))==int(emMIN):
        numberMIN=i
    if int(round(emwav[i]))==int(emMAX):
        numberMAX=i
    else:
        pass

ValuesToBeAveraged=[]
trace2=[]
for i in range(numberMIN,numberMAX):
    ValuesToBeAveraged.append(d2[i])

trace2=mean(ValuesToBeAveraged,axis=0)
waves2=exwav

with rc_context(fname='plotting'):
    x_plot.plot(waves2,trace2,color=cm.jet(60),label='{}'.format('Em: ('+str(emMIN)+'-'+str(emMAX)+') nm'))
    x_plot.legend(fontsize=6)
    x_plot.set_ylim(cpsMIN,cpsMAX)
    x_plot.set_xlim(exextentmin,exextentmax)
    x_plot.axvline(exavg,color='grey',linewidth=1,linestyle='-')
    x_plot.set_xlabel('Excitation Wavelength / nm')
    #y_plot.tick_params(rotation=45)

with rc_context(fname='plotting'):
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    cbax = fig.add_axes([0.8, 0.428, 0.02, 0.443])
    cbar=fig.colorbar(cb,cax=cbax, orientation='vertical')
    cbar.set_label('Counts / #/s ',rotation=270,labelpad=15)

for i in ['pdf','jpg']:
    fig.savefig(plotpath+'N-CND(gr)-EmExMap.{}'.format(i))
plt.show()
