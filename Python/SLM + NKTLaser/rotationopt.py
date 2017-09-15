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
from cameraclass import camera
from functions import SLMinit, useSLM
from beamshapes import LG, square, gausscomp, nothing, pic, test, displace, angled, rotator
import time
import matplotlib.pyplot as plt

intaim=np.load('calib.npy')
savefolder='23rd March/optexrot/'
slmpix=[512,512]
spfrq=0.15
wavelength=675
distribution=test(70,slmpix)*intaim
slm, slmstate =SLMinit(1)
tosend=useSLM(distribution,spfrq,wavelength,slm,slmstate,slmpix)
time.sleep(0.1)

cam=camera(20)
currentimg=cam.takeimagebw()
mask, cutmask, hexcor=optsetup(currentimg, '32.png', 'cutmask.bmp',1)
fig=plt.imshow(mask, cmap='gray')
plt.show()

valmin=999999999
angmin=10
valh=[]
for rotang in np.linspace(0, 60, 31):
    distribution=np.real((LG(3,2,180,slmpix)))
    distribution=displace(angled(rotator(distribution,rotang), 0, 0, slmpix), 0, 0, slmpix)
    scipy.misc.imsave(savefolder+"in{}.png".format(str(rotang)), np.abs(distribution))
    tosend=useSLM(distribution*intaim,spfrq,wavelength,slm,slmstate,slmpix)
    cam.autoexposurecut(hexcor[0], hexcor[1], hexcor[2])
    currentimg=cam.takeimagebw()
    plt.pause(0.05)    
    val, currentimg=compare(currentimg, mask, cutmask, hexcor)
    scipy.misc.imsave(savefolder+"out{}.png".format(str(rotang)), currentimg)
    valh.append((rotang, val))
    print(rotang, val)
    currentimg=np.expand_dims(currentimg,2)
    currentimg=np.concatenate((currentimg,currentimg,currentimg),axis=2)    
    fig.set_data(currentimg.astype(float)/255)
    plt.draw()
cam.closecam()
slm.close()


#scipy.misc.imsave('tmp2.bmp',cutmask)
#scipy.misc.imsave('tmp1.bmp',mask)
#scipy.misc.imsave('tmp.bmp',currentimg)

