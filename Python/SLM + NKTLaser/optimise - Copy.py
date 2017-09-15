#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:34:45 2017

@author: Andrei
"""
import matplotlib.image as mpimg
import scipy.misc
from autofunct import optsetup, compare
from scipy.ndimage.interpolation import zoom, rotate
import numpy as np
from camera import initcam, autoexposure, closecam, takeimagebw, takeimagecolour
from functions import SLMinit, useSLM
from beamshapes import LG, square, gausscomp, nothing, pic, test
import time
import matplotlib.pyplot as plt

slmpix=[512,512]
spfrq=0.15
wavelength=675
distribution=test(70,slmpix)
slm, slmstate =SLMinit(1)
tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(0.1)

cam=initcam(100)
currentimg=takeimagebw(cam)

mask, cutmask, hexcor=optsetup(currentimg, 'mask.bmp', 'cutmask.bmp')
fig=plt.imshow(mask, cmap='gray')
plt.show()

valmin=999999999
angmin=10
for rotip in np.linspace(0.4, 2.5, 22):
    distribution=np.real((LG(3,2,180,slmpix))*np.exp(1.j*rotip))
    tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
    plt.pause(0.05)
    
    autoexposure(cam)
    currentimg=takeimagebw(cam)
    
    #scipy.misc.imsave("{}.png".format(str(rotip)), currentimg)
    val, currentimg=compare(currentimg, mask, cutmask, hexcor)
    if val>valmin:
        valmin=val
        angmin=rotip
    print(rotip, val)
    fig.set_data(currentimg)
    plt.draw()
print(angmin)
closecam(cam)
slm.close()
'''
angmin=1.25
valmin=999999999
sizmin=0
for sizip in np.linspace(60, 300, 20):
    distribution=np.real((LG(3,2,sizip,slmpix))*np.exp(1.j*angmin))
    tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
    plt.pause(0.05)
    autoexposure(cam)
    currentimg=takeimagebw(cam)
    currentimg=cropcurrent(currentimg, cent, dist, angle)
    fig.set_data(currentimg)
    plt.draw()
    val, maskout, cutmaskout=compare(currentimg, mask, cutmask)
    if val<valmin:
        valmin=val
        sizmin=sizip
    print(val)    
print(sizmin)
'''


#scipy.misc.imsave('tmp2.bmp',cutmask)
#scipy.misc.imsave('tmp1.bmp',mask)
#scipy.misc.imsave('tmp.bmp',currentimg)

