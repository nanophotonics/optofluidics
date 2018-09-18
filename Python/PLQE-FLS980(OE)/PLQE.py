
# coding: utf-8

# In[1]:


import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
from IPython.display import HTML, display
import tabulate
import pandas as pd
import matplotlib.cm as cm
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import matplotlib.gridspec as gridspec
import plotly.plotly as py
import plotly.graph_objs as go
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


c_h=float(6.62607004e-34) #Plank's constant
c_c=float(299792458) #speed of light
c_e=float(1.60217662e-19) #elementary charge

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

def wav2en(x):
    return (c_h*c_c)/(x*c_e*10**-9) #wavelength in nm


# In[4]:


read_row_start=22 #start reading the txt file at

with open(path+'n-cnd(gr)_PLQE_blank.txt') as a: #reference
    reader1 = csv.reader(a, delimiter="\t")
    blank = list(reader1)
with open(path+'n-cnd(gr)_PLQE_sample.txt') as b: #sample
    reader2 = csv.reader(b, delimiter="\t")
    sample = list(reader2)

#reference
d1=[]
emwav=[]
for i in range(read_row_start,len(blank)):
    emwav.append(float(blank[i][0]))
for i in range(read_row_start,len(blank)):
    for j in range(1,len(blank[read_row_start])-1):
        d1.append(blank[i][j].replace("E", "e"))
d1=map(float,d1)
d1=np.array(map(abs,d1))

emwav=np.array(emwav)
energy=[]
for i in range(0,len(emwav)):
    energy.append(wav2en(emwav[i]))

#sample
d2=[]
for i in range(22,len(sample)):
    for j in range(1,len(sample[22])-1):
        d2.append(sample[i][j].replace("E", "e"))
d2=map(float,d2)
d2=np.array(map(abs,d2))


# In[5]:


#excitation integrals
wav1=float(375)
wav2=float(385)

for i in range(0,len(emwav)):
    if int(round(emwav[i]))==int(wav1):
        numberMIN=i
    if int(round(emwav[i]))==int(wav2):
        numberMAX=i
    else:
        pass

splice_blank=[]
splice_sample=[]
splice_emwav=[]

for i in range(numberMIN,numberMAX):
    splice_blank.append(d1[i])
    splice_sample.append(d2[i])
    splice_emwav.append(energy[i])

mean_blank = sum(ar(splice_blank)*ar(splice_emwav))/sum(ar(splice_blank))
mean_sample = sum(ar(splice_sample)*ar(splice_emwav))/sum(ar(splice_sample))
sigma_blank = np.sqrt(sum(ar(splice_blank)*(ar(splice_emwav) - mean_blank)**2)/sum(ar(splice_blank)))
sigma_sample = np.sqrt(sum(ar(splice_sample)*(ar(splice_emwav) - mean_sample)**2)/sum(ar(splice_sample)))

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gauss_blank(x):
    return gauss_function(x,*popt_blank)

def gauss_sample(x):
    return gauss_function(x,*popt_sample)

#curve fitting reference
popt_blank, pcov_blank = curve_fit(gauss_function, splice_emwav, splice_blank, p0 = [max(splice_blank), mean_blank, sigma_blank])
#curve fitting sample
popt_sample, pcov_sample = curve_fit(gauss_function, splice_emwav, splice_sample, p0 = [max(splice_sample), mean_sample, sigma_sample])

A_blank,err_blank=integrate.quad(gauss_blank, wav2en(wav1), wav2en(wav2))
A_sample,err_sample=integrate.quad(gauss_sample, wav2en(wav1), wav2en(wav2))

#emission integral
wav3=float(400)
wav4=float(650)

for i in range(0,len(emwav)):
    if int(round(emwav[i]))==int(wav3):
        numberMIN1=i
    if int(round(emwav[i]))==int(wav4):
        numberMAX1=i
    else:
        pass

splice_sample2=[]
splice_emwav2=[]

for i in range(numberMIN1,numberMAX1):
    splice_sample2.append(d2[i])
    splice_emwav2.append(energy[i])

mean_sample2 = sum(ar(splice_sample2)*ar(splice_emwav2))/sum(ar(splice_sample2))
sigma_sample2 = np.sqrt(sum(ar(splice_sample2)*(ar(splice_emwav2) - mean_sample2)**2)/sum(ar(splice_sample2)))

#curve fitting emission
popt_sample2, pcov_sample2 = curve_fit(gauss_function, splice_emwav2, splice_sample2, p0 = [max(splice_sample2), mean_sample2, sigma_sample2])

def gauss_sample2(x):
    return gauss_function(x,*popt_sample2)

A_sample2,err_sample2=integrate.quad(gauss_sample2, wav2en(wav3), wav2en(wav4))

print 'Area (1) sample emission (2) sample excitation (3) reference excitation:'
print A_sample2
print A_sample
print A_blank

PLQE=A_sample2/(A_blank-A_sample)
print 'PLQE:'
print PLQE


# In[6]:


#plot parameters
SG_Window=11 #odd number needed
SG_Order=2

cpsMAX=1.05*(np.amax(d1))
cpsMIN=-0
emextentmin=wav2en(800)
emextentmax=wav2en(350)

d3=SG(d2,SG_Window,SG_Order)
emwav3=SGW(energy,SG_Window-1)


# In[7]:


with rc_context(fname='plotting'):
    fig, ax1 = plt.subplots()
    left, bottom, width, height = [0.25, 0.4, 0.35, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax1.semilogy(energy,d1,color=cm.jet(60),label='{}'.format('water'))
    #ax1.plot(energy,gauss_blank(energy),color='green',linewidth=2,label='{}'.format('blank_fit'))
    ax1.semilogy(energy,d2,color='red',label='{}'.format('sample'))
    ax2.plot(emwav3,d3,color='red',linewidth=1.5,label='{}'.format('sample_magnify'))
    ax2.plot(emwav3,gauss_sample2(emwav3),color='green',linestyle='--',linewidth=1.5,label='{}'.format('sample_fit'))
    ax2.set_xlim(1.72,3.18)
    ax2.set_ylim(60,500)
    ax1.legend()
    ax1.set_ylim(90,cpsMAX*1.1)
    ax1.set_xlim(emextentmin,emextentmax)
    ax1.set_xlabel('Emission Energy / eV')
    ax1.set_ylabel('Counts per Second / #/s')
    ax1.set_title('N-CND(gr) PL (Excitation: 380 nm)',y=1.03)
    ax1.text(2.3, 300000,'{}'.format('PLQE: '+str(round(PLQE*100,2))+' %'),fontsize=10,color='green')

for i in ['pdf','jpg']:
    fig.savefig(plotpath+'N-CND(gr)-PLQE.{}'.format(i))
plt.show()
