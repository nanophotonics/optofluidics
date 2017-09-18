#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:34:45 2017

@author: Andrei
"""
import matplotlib.image as mpimg
import scipy.misc
import scipy.signal
from autofunct import optsetup, compare
from scipy.ndimage.interpolation import zoom, rotate
from beamshapes import test, LG, displace, angled
from functions import SLMinit, useSLM
from camera import takeimagebw, initcam, closecam, autoexposure
import numpy as np
import time
import matplotlib.pyplot as plt
import cv2

'''
Put identifier mode onto SLM.
Take a piture.
'''
intaim=np.load('calib.npy')
cam=initcam(4)
slm, slmstate =SLMinit(1)
slmpix=[512,512]
spfrq=0.15
wavelength=675

distribution=test(70, slmpix)*intaim
useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
img=takeimagebw(cam)
mask, cutmask, hexcor=optsetup(img, 'mask0.bmp','cutmask.bmp')
fig=plt.imshow(mask, cmap='gray')
plt.show()
xbest=0
ybest=0
xangbest=0
yangbest=0
optpara=0
distribution=LG(0, 0, 80, slmpix)

for xpos in np.linspace(-50,50,7):
    for xangle in np.linspace(-6*np.pi,6*np.pi,10):
        trydist=displace(angled(distribution, xangle, yangbest, slmpix),xpos,ybest,slmpix)
        useSLM(trydist*intaim,spfrq,wavelength,slm,slmstate,slmpix)
        #print(autoexposure(cam))
        img=takeimagebw(cam)
        val, currentimg =compare(img, mask,cutmask, hexcor)
        fig.set_data(currentimg)
        plt.draw()
        plt.pause(0.05)
        print(xpos, xangle, val)
        if val>optpara:
            xbest=xpos
            xangbest=xangle
            optpara=val

for ypos in np.linspace(-50,50,7):
    for yangle in np.linspace(-6*np.pi,6*np.pi,10):
        trydist=displace(angled(distribution, xangbest, yangle, slmpix),xbest,ypos,slmpix)
        useSLM(trydist*intaim,spfrq,wavelength,slm,slmstate,slmpix)
        #autoexposure(cam)
        img=takeimagebw(cam)
        val, currentimg=compare(img, mask,cutmask, hexcor)
        fig.set_data(currentimg)
        plt.draw()
        plt.pause(0.05)
        print(ypos, yangle, val)
        if val>optpara:
            ybest=ypos
            yangbest=yangle
            optpara=val
            
slm.close()
closecam(cam)
