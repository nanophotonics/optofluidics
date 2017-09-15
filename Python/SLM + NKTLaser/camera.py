# -*- coding: utf-8 -*-
"""
Created on Thu Dec 01 13:55:58 2016

@author: Hera
"""

import sys
import time
#from PySide.QtGui import *
from instrumental import instrument, gui, errors
import scipy.misc
import numpy as np
from autofunct import cropcurrent

def initcam(exp):
    cam = instrument('uc480')
    cam.start_live_video(exposure_time="{} millisecond".format(str(exp)))
    return cam

def closecam(cam):
        cam.stop_live_video()
        cam.close()


def takeimagecolour(cam):
    frame_ready = cam.wait_for_frame()
    if frame_ready:
        a = cam.latest_frame()
        print(cam._get_exposure())
        a=a[:,:,0:3]
    else:
        a=0
        raise NameError('Failed to get frame.')
    return a
        
def takeimagebw(cam):
    frame_ready = cam.wait_for_frame()
    if frame_ready:
        a = cam.latest_frame()
        a=a[:,:,0:3]
        a=np.sum(a,2)
    else:
        a=0
        raise NameError('Failed to get frame.')
    return a/3

def autoexposure(cam):
    emax=20
    emin=0.01
    while (emax-emin>=0.05):
        cam.start_live_video(exposure_time="{} millisecond".format(str((emax+emin)/2)))
        a=takeimagecolour(cam)
        if np.sum(a==255)>10:
            emax=(emax+emin)/2
        else:
            emin=(emax+emin)/2
    cam.start_live_video(exposure_time="{} millisecond".format(str(emin)))
    a=takeimagecolour(cam)
    return emin
    
def autoexposurecut(cam, cent, dist, angle):
    emax=50
    emin=0.01
    while (emax-emin>=0.05):
        cam.start_live_video(exposure_time="{} millisecond".format(str((emax+emin)/2)))
        a=takeimagecolour(cam)
        a=cropcurrent(a, cent, dist, angle)
        if np.sum(a==255)>10:
            emax=(emax+emin)/2
        else:
            emin=(emax+emin)/2
    cam.start_live_video(exposure_time="{} millisecond".format(str(emin)))
    a=takeimagecolour(cam)
    return emin    
    
def exposure(cam, expo):    
    cam.start_live_video(exposure_time="{} millisecond".format(str(expo)))
#print(np.shape(a))
#print(np.max(np.max(a)))